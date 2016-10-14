!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module momentum_diagnostics

  use boundary_conditions
  use coriolis_module, only : two_omega => coriolis
  use diagnostic_source_fields
  use diagnostic_fields
  use field_derivatives
  use field_options
  use fields
  use fldebug
  use geostrophic_pressure
  use global_parameters, only : OPTION_PATH_LEN
  use multimaterial_module
  use solvers
  use sparsity_patterns_meshes
  use spud
  use state_fields_module
  use sediment, only : get_n_sediment_fields, get_sediment_item
  use state_module
  use initialise_fields_module
  use parallel_tools
  
  implicit none
  
  private
  
  public :: calculate_strain_rate, calculate_bulk_viscosity, calculate_strain_rate_second_invariant, &
	    calculate_sediment_concentration_dependent_viscosity, &
	    calculate_buoyancy, calculate_coriolis, calculate_tensor_second_invariant, &
	    calculate_imposed_material_velocity_source, &
            calculate_imposed_material_velocity_absorption, &
            calculate_scalar_potential, calculate_projection_scalar_potential, &
            calculate_geostrophic_velocity, calculate_atmosphere_forcing_vector, calculate_viscous_dissipation
           
contains

  subroutine calculate_strain_rate(state, t_field)
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: t_field
    
    type(vector_field), pointer :: source_field
    type(vector_field), pointer :: positions

    positions => extract_vector_field(state, "Coordinate")
    source_field => vector_source_field(state, t_field)

    call check_source_mesh_derivative(source_field, "strain_rate")

    call strain_rate(source_field, positions, t_field)

  end subroutine calculate_strain_rate

  subroutine calculate_strain_rate_second_invariant(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: positions
    type(vector_field), pointer :: velocity

    type(tensor_field) :: strain_rate_tensor

    ewrite(1,*) 'In calculate_strain_rate_second_invariant'
    positions => extract_vector_field(state, "Coordinate")
    velocity  => extract_vector_field(state, "IteratedVelocity")

    ! Allocate strain_rate tensor:
    call allocate(strain_rate_tensor, s_field%mesh, name="strain_rate_II")

    call check_source_mesh_derivative(velocity, "strain_rate_second_invariant")

    ! Calculate strain_rate and second invariant:
    call strain_rate(velocity, positions, strain_rate_tensor)
    call tensor_second_invariant(strain_rate_tensor, s_field)

    ! Clean-up:
    call deallocate(strain_rate_tensor)

    ! Prin min and max:
    ewrite_minmax(s_field) 

  end subroutine calculate_strain_rate_second_invariant

  subroutine calculate_sediment_concentration_dependent_viscosity(state, t_field)
    ! calculates viscosity based upon total sediment concentration
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: t_field
    
    type(scalar_field_pointer), dimension(:), allocatable :: sediment_concs
    type(tensor_field), pointer :: zero_conc_viscosity
    type(scalar_field) :: rhs
    integer :: sediment_classes, i
    character(len = FIELD_NAME_LEN) :: field_name
    
    ewrite(1,*) 'In calculate_sediment_concentration_dependent_viscosity'

    sediment_classes = get_n_sediment_fields()

    if (sediment_classes > 0) then
	allocate(sediment_concs(sediment_classes))
	
	call get_sediment_item(state, 1, sediment_concs(1)%ptr)
	
	call allocate(rhs, sediment_concs(1)%ptr%mesh, name="Rhs")
	call set(rhs, 1.0)
	
	! get sediment concentrations and remove c/0.65 from rhs
	do i=1, sediment_classes
	   call get_sediment_item(state, i, sediment_concs(i)%ptr)
	   call addto(rhs, sediment_concs(i)%ptr, scale=-(1.0/0.65))
	end do
	
	! raise rhs to power of -1.625
	do i = 1, node_count(rhs)
	   call set(rhs, i, node_val(rhs, i)**(-1.625))
	end do
	
	! check for presence of ZeroSedimentConcentrationViscosity field
	if (.not. has_tensor_field(state, "ZeroSedimentConcentrationViscosity")) then
	   FLExit("You must specify an zero sediment concentration viscosity to be able &
		&to calculate sediment concentration dependent viscosity field values")
	endif
	zero_conc_viscosity => extract_tensor_field(state, 'ZeroSedimentConcentrationViscosity')
	
	call set(t_field, zero_conc_viscosity)
	call scale(t_field, rhs)
	ewrite_minmax(t_field) 

	deallocate(sediment_concs)
	call deallocate(rhs)
    else
	ewrite(1,*) 'No sediment in problem definition'
    end if  
  end subroutine calculate_sediment_concentration_dependent_viscosity
  
  subroutine calculate_tensor_second_invariant(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(tensor_field), pointer :: source_field

    source_field => tensor_source_field(state, s_field)

    call tensor_second_invariant(source_field, s_field)

  end subroutine calculate_tensor_second_invariant

  subroutine calculate_viscous_dissipation(state, s_field)
    ! A routine to calculate the viscous dissipation. Currently
    ! assumes a constant viscosity tensor:
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: positions
    type(vector_field), pointer :: velocity
    type(tensor_field), pointer :: viscosity

    type(scalar_field) :: velocity_divergence
    type(scalar_field) :: viscosity_component, viscosity_component_remap
    type(tensor_field) :: strain_rate_tensor

    integer :: dim1, dim2, node
    real :: val

    ewrite(1,*) 'In calculate_viscous_dissipation'

    ! Extract velocity field from state - will be used to calculate strain-
    ! rate tensor:
    velocity => extract_vector_field(state, "NonlinearVelocity")
    ! Check velocity field is not on a discontinous mesh:
    call check_source_mesh_derivative(velocity, "Viscous_Dissipation")

    ! Extract positions field from state:
    positions => extract_vector_field(state, "Coordinate")

    ! Allocate and initialize strain rate tensor:
    call allocate(strain_rate_tensor, s_field%mesh, "Strain_Rate_VD")
    call zero(strain_rate_tensor)

    ! Calculate strain rate tensor:
    call strain_rate(velocity, positions, strain_rate_tensor)

    ! Calculate velocity divergence for correct definition of stress:
    call allocate(velocity_divergence, s_field%mesh, 'Velocity_divergence')
    call div(velocity, positions, velocity_divergence)
    ewrite_minmax(velocity_divergence)

    ! Extract viscosity from state and remap to s_field mesh:
    viscosity => extract_tensor_field(state, "Viscosity")
    ! Extract first component of viscosity tensor from full tensor:
    !*** This is not ideal - only valid for constant viscosity tensors
    !*** though they can still vary spatially and temporally.
    viscosity_component = extract_scalar_field(viscosity,1,1)  
    call allocate(viscosity_component_remap, s_field%mesh, "RemappedViscosityComponent")
    call remap_field(viscosity_component, viscosity_component_remap)

    ! Calculate viscous dissipation (scalar s_field):
    do node=1,node_count(s_field)
       val = 0.
       do dim1 = 1, velocity%dim
	  do dim2 = 1, velocity%dim
	     if(dim1==dim2) then
		! Add divergence of velocity term to diagonal only: 
		val = val + 2.*node_val(viscosity_component_remap, node) * & 
		     & (node_val(strain_rate_tensor,dim1,dim2,node)	 - &
		     & 1./3. * node_val(velocity_divergence, node))**2
	     else
		val = val + 2.*node_val(viscosity_component_remap, node) * & 
		     & node_val(strain_rate_tensor,dim1,dim2,node)**2	
	     end if
	  end do
       end do
       call set(s_field, node, val)
    end do

    ewrite_minmax(s_field)

    ! Deallocate:
    call deallocate(strain_rate_tensor)
    call deallocate(viscosity_component_remap)
    call deallocate(velocity_divergence)

  end subroutine calculate_viscous_dissipation

  subroutine calculate_bulk_viscosity(states, t_field)
    type(state_type), dimension(:), intent(inout) :: states
    type(tensor_field), intent(inout) :: t_field
    character(len = OPTION_PATH_LEN) :: mean_type
    
    call get_option(trim(complete_field_path(trim(t_field%option_path))) // &
		    "/algorithm[0]/mean/name", mean_type, default="arithmetic")

    call calculate_bulk_property(states, t_field, "MaterialViscosity", &
      & mean_type = mean_type, momentum_diagnostic = .true.)

  end subroutine calculate_bulk_viscosity
    
  subroutine calculate_atmosphere_forcing_vector(state, v_field, dt)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    real, intent(in) :: dt
    
    integer :: i, j, stat, ele, dim
    real  :: time_scale, z_base, x_base_r, x_base_l, y_base_r, y_base_l, z_max, x_max, x_min, y_max, y_min, xi, X_val, coeff
    real, allocatable, dimension(:) :: absorption_valx, absorption_valz, absorption_valy, time_scale_tab
    type(vector_field), pointer :: X
    type(vector_field) :: X_local
    character(len=OPTION_PATH_LEN) :: absorption_path
    
    ewrite(2,*) 'Inside calculate_atmosphere_forcing_vector '//trim(v_field%name)
  
    call zero(v_field)
    
    X=>extract_vector_field(state,"Coordinate")
    call allocate(X_local,X%dim,v_field%mesh,"LocalMesh")
    X_local = get_remapped_coordinates(X, v_field%mesh)

    absorption_path = trim(v_field%option_path)//"/diagnostic/algorithm::atmosphere_forcing_vector"
    
    if (have_option(trim(absorption_path)//"/sponge_layer_velocity_absorption")) then
    
      call get_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/coefficient',coeff)
    
      if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_right')) &
          call get_option (trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_right',x_base_r)
      if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_left')) &
          call get_option (trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_left',x_base_l)
      if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_right')) &
          call get_option (trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_right',y_base_r)
      if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_left')) &
          call get_option (trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_left',y_base_l)
      if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/z_sponge')) &
          call get_option (trim(absorption_path)//'/sponge_layer_velocity_absorption/z_sponge',z_base)

      allocate(absorption_valx(1:X_local%dim), absorption_valy(1:X_local%dim), absorption_valz(1:X_local%dim))

      dim=X%dim
      x_max=maxval(X_local%val(1,:))	  
      x_min=minval(X_local%val(1,:))	  
      y_max=maxval(X_local%val(2,:))	  
      y_min=minval(X_local%val(2,:))	  
      z_max=maxval(X_local%val(dim,:))

      call allmax(x_max)
      call allmax(y_max)
      call allmax(z_max)
      call allmin(x_min)
      call allmin(y_min)

      do ele=1,node_count(v_field)
        absorption_valx=0.
        absorption_valy=0.
        absorption_valz=0.
        
        X_val=node_val(X_local,1,ele)
        if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_right') .and. X_val > x_base_r) then
          xi=(X_val - x_base_r)/(x_max - x_base_r)
          if (X_val >= x_base_r .and. xi >= 0. .and. xi < 0.75) then
            absorption_valx=1./dt*sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
          else if (X_val >= x_base_r .and. xi > 0.75 .and. xi <= 1.) then
            absorption_valx=1./dt
          else
            absorption_valx=0.
          endif
        endif

        if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_left') .and. X_val < x_base_l) then
          xi=(x_base_l - X_val)/(x_base_l - x_min)
          if (X_val <= x_base_l .and. xi >= 0. .and. xi < 0.75) then
            absorption_valx=1./dt*sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
          else if (X_val <= x_base_l .and. xi > 0.75 .and. xi <= 1.) then
            absorption_valx=1./dt
          else
            absorption_valx=0.
          endif
        endif	
        
        if (dim > 2) then
          X_val=node_val(X_local,2,ele)
          if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_right') .and. X_val > y_base_r) then
            xi=(X_val - y_base_r)/(y_max - y_base_r)
            if (X_val >= y_base_r .and. xi >= 0. .and. xi < 0.75) then
              absorption_valy=1./dt*sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
            else if (X_val >= y_base_r .and. xi > 0.75 .and. xi <= 1.) then
              absorption_valy=1./dt
            else
              absorption_valy=0.
            endif
          endif

          if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_left') .and. X_val < y_base_l) then
            xi=(y_base_l - X_val)/(y_base_l - y_min)
            if (X_val <= y_base_l .and. xi >= 0. .and. xi < 0.75) then
              absorption_valy=1./dt*sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
            else if (X_val <= x_base_l .and. xi > 0.75 .and. xi <= 1.) then
              absorption_valy=1./dt
            else
              absorption_valy=0.
            endif
          endif   
        endif
        
        X_val=node_val(X_local,X_local%dim,ele)
        if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/z_sponge') .and. X_val > z_base) then
          xi=(X_val - z_base)/(z_max - z_base)
          if (X_val >= z_base .and. xi >= 0. .and. xi < 0.75) then
            absorption_valz=1./dt*sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
          else if (X_val >= z_base .and. xi > 0.75 .and. xi <= 1.) then
            absorption_valz=1./dt
          else
            absorption_valz=0.
          endif
        endif  

        absorption_valx=coeff*sqrt((absorption_valx**2. + absorption_valy**2. + absorption_valz**2.)/real(dim))
        call addto(v_field, ele, absorption_valx)
      enddo
    
      deallocate (absorption_valx, absorption_valy, absorption_valz)
    end if 
    
    if (have_option(trim(absorption_path)//"/velocity_nudging")) then

      allocate(time_scale_tab(1:X_local%dim))
      call get_option(trim(absorption_path)//"/velocity_nudging/time_scale", time_scale)
      do j=1,X_local%dim
        time_scale_tab(j)=1./time_scale
      enddo

      call addto(v_field, time_scale_tab)    
      deallocate(time_scale_tab)
    
    end if
    
    call deallocate (X_local)
    
  end subroutine calculate_atmosphere_forcing_vector

  subroutine calculate_imposed_material_velocity_source(states, state_index, v_field)
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(inout) :: v_field
    
    logical :: prescribed
    integer :: i, stat
    type(vector_field), pointer :: absorption, mat_vel
  
    call zero(v_field)
    
    do i = 1, size(states)
  
      mat_vel => extract_vector_field(states(i), "MaterialVelocity", stat)
      
      if(stat==0) then
      
        call add_scaled_material_property(states(i), v_field, mat_vel, &
                                          momentum_diagnostic=.true.)
                                          
      else
        ! alternatively use the Velocity field from the state
        
        mat_vel => extract_vector_field(states(i), "Velocity", stat)
        
        if(stat==0) then
          prescribed = have_option(trim(mat_vel%option_path)//"/prescribed")
          
          if(prescribed.and.(.not.aliased(mat_vel))) then
            ! but make sure it's prescribed and not aliased
          
            call add_scaled_material_property(states(i), v_field, mat_vel, &
                                              momentum_diagnostic=.true.)
                                              
          end if
        
        end if
        
      end if

    end do

    absorption => extract_vector_field(states(state_index), "VelocityAbsorption")
    call scale(v_field, absorption)
  
  end subroutine calculate_imposed_material_velocity_source

  subroutine calculate_imposed_material_velocity_absorption(states, v_field)
    type(state_type), dimension(:), intent(inout) :: states
    type(vector_field), intent(inout) :: v_field
    
    logical :: prescribed
    integer :: i, stat
    real :: dt
    real, dimension(v_field%dim) :: factor
    type(vector_field) :: temp_abs
    type(vector_field), pointer :: mat_vel
    
    call get_option("/timestepping/timestep", dt)
    call get_option(trim(complete_field_path(trim(v_field%option_path))) // &
                    "/algorithm[0]/relaxation_factor", factor, default=spread(1.0, 1, v_field%dim))
    
    call allocate(temp_abs, v_field%dim, v_field%mesh, "TemporaryAbsorption", &
                  field_type=FIELD_TYPE_CONSTANT)
    call set(temp_abs, factor/dt)
        
    call zero(v_field)
    
    do i = 1, size(states)
  
      mat_vel => extract_vector_field(states(i), "MaterialVelocity", stat)
      
      if(stat==0) then
      
        call add_scaled_material_property(states(i), v_field, temp_abs, &
                                          momentum_diagnostic=.true.)

      else
        ! alternatively use the Velocity field from the state
        
        mat_vel => extract_vector_field(states(i), "Velocity", stat)
        
        if(stat==0) then
          prescribed = have_option(trim(mat_vel%option_path)//"/prescribed")
          
          if(prescribed.and.(.not.aliased(mat_vel))) then
            ! but make sure it's prescribed and not aliased
          
            call add_scaled_material_property(states(i), v_field, temp_abs, &
                                              momentum_diagnostic=.true.)
            
          end if
        
        end if
        
      end if

    end do
    
    call deallocate(temp_abs)
    
  end subroutine calculate_imposed_material_velocity_absorption

  subroutine calculate_buoyancy(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    integer :: i, stat
    real :: gravity_magnitude
    type(scalar_field), pointer :: buoyancy_density
    type(vector_field), pointer :: gravity
  
    ewrite(1, *) "In calculate_buoyancy"
    
    buoyancy_density => extract_scalar_field(state, "VelocityBuoyancyDensity", stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "Warning: Cannot calculate Buoyancy without VelocityBuoyancyDensity field"
      call zero(v_field)
      ewrite(1, *) "Exiting calculate_buoyancy"
      return
    end if    
    ewrite_minmax(buoyancy_density)
    
    gravity => extract_vector_field(state, "GravityDirection", stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "Warning: Cannot calculate Buoyancy without GravityDirection field"
      call zero(v_field)
      ewrite(1, *) "Exiting calculate_buoyancy"
      return
    end if    
    ewrite_minmax(gravity)
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    ewrite(2, *) "Gravity magnitude = ", gravity_magnitude
    
    if(.not. v_field%mesh == buoyancy_density%mesh) then
      ewrite(-1, *) "VelocityBuoyancyDensity mesh: " // trim(buoyancy_density%mesh%name)
      FLExit("Buoyancy must be on the VelocityBuoyancyDensity mesh")
    end if
    
    do i = 1, node_count(v_field)
      call set(v_field, i, node_val(gravity, i) * node_val(buoyancy_density, i) * gravity_magnitude)
    end do
    
    ewrite(1, *) "Exiting calculate_buoyancy"
  
  end subroutine calculate_buoyancy
  
  subroutine calculate_coriolis(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
   
    character(len = OPTION_PATH_LEN) :: base_path
    
    base_path = trim(complete_field_path(v_field%option_path)) // "/algorithm"

    if(have_option(trim(base_path) // "/consistent_interpolation")) then
      call compute_coriolis_ci(state, v_field)
    else if(have_option(trim(base_path) // "/galerkin_projection")) then
      if(have_option(trim(base_path) // "/galerkin_projection/lump_mass")) then
        call compute_coriolis_gp_lumped(state, v_field)
      else
        call compute_coriolis_gp(state, v_field, option_path = trim(base_path) // "/galerkin_projection")
      end if
    else
      FLAbort("Failed to determine interpolation method")
    end if  
      
  end subroutine calculate_coriolis
  
  subroutine compute_coriolis_ci(state, coriolis)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: coriolis
  
    integer :: i
    type(vector_field) :: positions, velocity_remap
    type(vector_field), pointer :: velocity
    
    positions = get_nodal_coordinate_field(state, coriolis%mesh)
    velocity => extract_vector_field(state, "Velocity")
    
    if(velocity%mesh == coriolis%mesh) then
      velocity_remap = velocity
      call incref(velocity_remap)
    else
      call allocate(velocity_remap, velocity%dim, coriolis%mesh, "VelocityRemap")
      call remap_field(velocity, velocity_remap)
    end if
    
    do i = 1, node_count(coriolis)
      call set(coriolis, i, coriolis_val(node_val(positions, i), node_val(velocity_remap, i)))
    end do
    
    call deallocate(positions)
    call deallocate(velocity_remap)
    
  end subroutine compute_coriolis_ci
  
  subroutine compute_coriolis_gp(state, coriolis, option_path)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: coriolis
    character(len = *), optional, intent(in) :: option_path
  
    integer :: i
    type(csr_matrix), pointer :: mass
    type(vector_field) :: rhs
    type(vector_field), pointer :: positions, velocity
    
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    
    mass => get_mass_matrix(state, coriolis%mesh)
    call allocate(rhs, coriolis%dim, coriolis%mesh, "CoriolisRhs")

    call zero(rhs)
    do i = 1, ele_count(rhs)
      call assemble_coriolis_ele(i, positions, velocity, rhs)
    end do

    call petsc_solve(coriolis, mass, rhs, option_path = option_path)

    call deallocate(rhs)
  
  end subroutine compute_coriolis_gp
  
  subroutine compute_coriolis_gp_lumped(state, coriolis)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: coriolis
    
    integer :: i
    type(scalar_field), pointer :: masslump
    type(vector_field), pointer :: positions, velocity
    
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    
    masslump => get_lumped_mass(state, coriolis%mesh)

    call zero(coriolis)
    do i = 1, ele_count(coriolis)
      call assemble_coriolis_ele(i, positions, velocity, coriolis)
    end do
    
    do i = 1, coriolis%dim
      coriolis%val(i,:) = coriolis%val(i,:) / masslump%val
    end do
    
  end subroutine compute_coriolis_gp_lumped
    
  subroutine assemble_coriolis_ele(ele, positions, velocity, rhs)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    type(vector_field), intent(inout) :: rhs

    real, dimension(ele_ngi(rhs, ele)) :: detwei

    call transform_to_physical(positions, ele, detwei = detwei)

    call addto(rhs, ele_nodes(rhs, ele), &
      & shape_vector_rhs(ele_shape(rhs, ele), &
        & coriolis_val(ele_val_at_quad(positions, ele), ele_val_at_quad(velocity, ele)), &
      & detwei))

  end subroutine assemble_coriolis_ele
  
  subroutine calculate_scalar_potential(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: source_field

    source_field => vector_source_field(state, s_field)
    call geopressure_decomposition(state, source_field, s_field, &
      & option_path = trim(complete_field_path(s_field%option_path)) // "/algorithm")

  end subroutine calculate_scalar_potential
  
  subroutine calculate_projection_scalar_potential(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    character(len = OPTION_PATH_LEN) :: bcfield_name, path
    type(scalar_field), pointer :: gp
    type(vector_field), pointer :: bcfield, source_field

    source_field => vector_source_field(state, s_field, index = 1)
    
    path = trim(complete_field_path(s_field%option_path)) // "/algorithm"

    if(have_option(trim(path) // "/bc_field")) then
      call get_option(trim(path) // "/bc_field/name", bcfield_name)
      bcfield => extract_vector_field(state, bcfield_name)
    else
      bcfield => source_field
    end if
    
    if(have_option(trim(path) // "/source_field_2_name")) then
      gp => scalar_source_field(state, s_field, index = 2)
      call projection_decomposition(state, source_field, s_field, &
        & bcfield = bcfield, gp = gp, option_path = path)
    else
      call projection_decomposition(state, source_field, s_field, &
        & bcfield = bcfield, option_path = path)
    end if

  end subroutine calculate_projection_scalar_potential

  subroutine calculate_geostrophic_velocity(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path
    integer :: stat
    real :: scale_factor
    type(scalar_field), pointer :: source_field
    type(cmc_matrices) :: matrices
    type(vector_field), pointer :: velocity
    
    source_field => scalar_source_field(state, v_field)
    velocity => extract_vector_field(state, "Velocity")
    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    call allocate(matrices, state, velocity, source_field, option_path = path, add_cmc = .false.)
    
    call geostrophic_velocity(matrices, state, v_field, source_field) 
    
    call deallocate(matrices)
    
    call get_option(trim(path) // "/scale_factor", scale_factor, stat = stat)
    if(stat == SPUD_NO_ERROR) call scale(v_field, scale_factor)
  
  end subroutine calculate_geostrophic_velocity

end module momentum_diagnostics
