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

module diagnostic_fields_surface

  use fields
  use fldebug
  use spud
  use field_options
  use parallel_tools
  use global_parameters, only:FIELD_NAME_LEN, current_time, OPTION_PATH_LEN
  use diagnostic_fields, only: safe_set
  use boundary_conditions
  use state_module
  use equation_of_state
  use filter_diagnostics
  use column_module
  use node_ownership
  use populate_state_module
  use vtk_interfaces, only: vtk_write_fields
  use write_state_module, only: include_scalar_field_in_vtu

  implicit none
  
  private
  public :: calculate_diagnostic_variable_surface, write_surface
  
  interface calculate_diagnostic_variable_surface
    module procedure diagnostic_variable_surface_single, 		&
        diagnostic_variable_surface_multiple_non_indexed, 		&
	diagnostic_variable_surface_multiple_indexed
  end interface calculate_diagnostic_variable_surface
  
contains
  
  subroutine diagnostic_variable_surface_single (state, s_field_name)
  
    type(state_type), intent(inout) :: state
    character(len = *), intent(in) :: s_field_name
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call calculate_diagnostic_variable_surface(states, 1, s_field_name)
    state = states(1)
  
  end subroutine diagnostic_variable_surface_single
  
  subroutine diagnostic_variable_surface_multiple_non_indexed(states, s_field_name)
  
    type(state_type), dimension(:), target, intent(inout) :: states
    character(len = *), intent(inout) :: s_field_name
    
    call calculate_diagnostic_variable_surface(states, 1, s_field_name)
  
  end subroutine diagnostic_variable_surface_multiple_non_indexed

  subroutine diagnostic_variable_surface_multiple_indexed (states, state_index, sfield_name)
    implicit none
    type(state_type), dimension(:), target, intent(inout) :: states
    character(len = *), intent(in) :: sfield_name
    integer, intent(in) :: state_index

    type(state_type), pointer :: state => null()
    type(scalar_field), pointer :: surface_field, pfield
    type(vector_field), pointer :: position
    type(vector_field) :: surface_position
    type(mesh_type), pointer :: surface_mesh
    integer, dimension(:), pointer:: surface_element_list   
    integer, dimension(:), allocatable:: surface_ids
    integer :: stat, shape_option(2)
    character(len=OPTION_PATH_LEN) :: bc_path, bc_type, option_path
    character(len=FIELD_NAME_LEN) :: parent_name
    logical :: recalculate

    state => states(state_index)

    ! Surface scalars:        
    option_path='/material_phase[0]/scalar_field::'//trim(sfield_name)//'/diagnostic'
    if (have_option(trim(option_path)//'/boundary_conditions::diagnostic')) then
      
      bc_path = trim(option_path)//'/boundary_conditions::diagnostic'
      call get_option(trim(bc_path)//"/parent_field_name", parent_name)
      call get_option(trim(bc_path)//"/type[0]/name", bc_type)
    
      pfield => extract_scalar_field(state, trim(parent_name))
      position => extract_vector_field(state,"Coordinate")
      surface_field=>extract_surface_field(pfield, trim(bc_type), name=trim(bc_type))
    
      call get_boundary_condition(pfield, trim(bc_type), surface_mesh=surface_mesh, &
                    surface_element_list=surface_element_list)

      surface_position=get_coordinates_remapped_to_surface(position, surface_mesh, surface_element_list) 

      select case(trim(bc_type))
        
	case("surface_temperature")
	  call get_surface_temperature(state, surface_field, pfield, surface_element_list, stat)
        
	case("surface_precipitation")
	  call get_surface_precipitation(state, surface_field, pfield, surface_element_list, stat)
        
	case("cumulated_surface_precipitation")
	  call get_cumulated_surface_precipitation(state, surface_field, pfield, surface_mesh, surface_element_list, stat)

	case("cloud_cover")
	  call get_cloud_cover(state, surface_field, surface_position, pfield, surface_element_list, option_path, stat)

	case("liquid_water_path")
	  call get_liquid_water_path(state, surface_field, surface_position, pfield, surface_element_list, option_path, stat)

	case("sensible_flux")
	  call get_sensible_flux(state, surface_field, pfield, surface_element_list, stat)

	case("latent_flux")
	  call get_latent_flux(state, surface_field, pfield, surface_element_list, stat)

	case("rain_flux")
	  call get_rain_flux(state, surface_field, pfield, surface_element_list, stat)
      
        case default
        ! Nothing to be done here
      
      end select
      
      call deallocate(surface_position)
      
    end if

  end subroutine diagnostic_variable_surface_multiple_indexed

  subroutine get_sensible_flux(state, s_field, p_field, surface_element_list, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    integer :: nbcs, i
    character(len=OPTION_PATH_LEN) :: bc_path
    character(len=FIELD_NAME_LEN) :: bc_name, bc_type
    type(scalar_field), pointer :: flux
    
    call zero(s_field)
    
    bc_path="/material_phase[0]/scalar_field::"//trim(p_field%name)//"/prognostic/boundary_conditions"
    nbcs=option_count(trim(bc_path))

    do i=0, nbcs-1

      call get_option(trim(bc_path)//"["//int2str(i)//"]/type[0]/name", bc_type)
      
      if (trim(bc_type) == "surface_ocean_COARE3") then
        call get_option(trim(bc_path)//"["//int2str(i)//"]/name", bc_name)
        flux => extract_surface_field(p_field, bc_name, "sensible_heat_flux", stat=stat)
    
        call set(s_field,flux)
      endif
      
    enddo
    
  end subroutine get_sensible_flux

  subroutine get_latent_flux(state, s_field, p_field, surface_element_list, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    integer :: nbcs, i
    character(len=OPTION_PATH_LEN) :: bc_path
    character(len=FIELD_NAME_LEN) :: bc_name, bc_type
    type(scalar_field), pointer :: flux
    
    call zero(s_field)
    
    bc_path="/material_phase[0]/scalar_field::"//trim(p_field%name)//"/prognostic/boundary_conditions"
    nbcs=option_count(trim(bc_path))

    do i=0, nbcs-1

      call get_option(trim(bc_path)//"["//int2str(i)//"]/type[0]/name", bc_type)
      
      if (trim(bc_type) == "surface_ocean_COARE3") then
        call get_option(trim(bc_path)//"["//int2str(i)//"]/name", bc_name)
        flux => extract_surface_field(p_field, bc_name, "latent_heat_flux", stat=stat)
    
        call set(s_field,flux)
      endif
      
    enddo
    
  end subroutine get_latent_flux

  subroutine get_rain_flux(state, s_field, p_field, surface_element_list, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    integer :: nbcs, i
    character(len=OPTION_PATH_LEN) :: bc_path
    character(len=FIELD_NAME_LEN) :: bc_name, bc_type
    type(scalar_field), pointer :: flux
    
    call zero(s_field)
    
    bc_path="/material_phase[0]/scalar_field::"//trim(p_field%name)//"/prognostic/boundary_conditions"
    nbcs=option_count(trim(bc_path))

    do i=0, nbcs-1

      call get_option(trim(bc_path)//"["//int2str(i)//"]/type[0]/name", bc_type)
      
      if (trim(bc_type) == "surface_ocean_COARE3") then
        call get_option(trim(bc_path)//"["//int2str(i)//"]/name", bc_name)
        flux => extract_surface_field(p_field, bc_name, "rain_heat_flux", stat=stat)
    
        call set(s_field,flux)
      endif
      
    enddo
    
  end subroutine get_rain_flux

  subroutine get_surface_temperature(state, s_field, p_field, surface_element_list, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    type(scalar_field) :: temperature
    type(scalar_field), pointer :: local
    
    call zero(s_field)
    
    call allocate(temperature,p_field%mesh,"LocalTemperature")
    call compressible_eos(state, temperature=temperature)
    
    call remap_field_to_surface (temperature, s_field, surface_element_list)
    
    call deallocate(temperature)

  end subroutine get_surface_temperature

  subroutine get_surface_precipitation(state, s_field, p_field, surface_element_list, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    integer :: qstat, vstat
    real :: rho_w=1000.
    type(scalar_field) :: precip, density_remap
    type(scalar_field), pointer :: density, QR, QR_sink
    
    call zero(s_field)
    
    call allocate(precip,p_field%mesh,"LocalPrecipitation")
    call allocate(density_remap,p_field%mesh,"RemappedDensity")

    density=>extract_scalar_field(state, "Density")
    QR=>extract_scalar_field(state, "Qrain", qstat)
    QR_sink=>extract_scalar_field(state, "QrainSinkingVelocity", vstat)
    
    call safe_set(state,density_remap,density)
    
    stat=0
    if (qstat/=0 .or. vstat/=0) stat=-1
    if (stat == 0) then      
      call set(precip,QR)
      call scale(precip,QR_sink)
      call scale(precip,density_remap)
      call scale(precip,3.6e6/rho_w)  !Conversion to mm/h
    
      call remap_field_to_surface (precip, s_field, surface_element_list)
    endif
    
    call deallocate(precip)
    call deallocate(density_remap)
    
  end subroutine get_surface_precipitation

  subroutine get_cumulated_surface_precipitation(state, s_field, p_field, surface_mesh, surface_element_list, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    type(mesh_type), intent(inout) :: surface_mesh
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    integer :: qstat, vstat
    real :: dt, rho_w=1000.
    type(scalar_field) :: precip, s_precip, density_remap
    type(scalar_field), pointer :: density, QR, QR_sink
    
    call allocate(precip,p_field%mesh,"LocalPrecipitation")
    call allocate(density_remap,p_field%mesh,"RemappedDensity")
    call allocate(s_precip,s_field%mesh,"LocalSurfacePrecipitation")

    call zero(s_precip)
       
    call get_option("/timestepping/timestep", dt)

    density=>extract_scalar_field(state, "Density")
    QR=>extract_scalar_field(state, "Qrain", qstat)
    QR_sink=>extract_scalar_field(state, "QrainSinkingVelocity", vstat)
    
    call safe_set(state,density_remap,density)
    
    stat=0
    if (present(stat) .and. qstat/=0 .or. vstat/=0) stat=-1
    if (stat == 0) then      
      call set(precip,QR)
      call scale(precip,QR_sink)
      call scale(precip,density_remap)
      call scale(precip,1000.*dt/rho_w)  !Conversion to mm
    
      call remap_field_to_surface(precip, s_precip, surface_element_list)

      call addto(s_field, s_precip)
    endif
    
    call deallocate(precip)    
    call deallocate(density_remap)
    call deallocate(s_precip)    
    
  end subroutine get_cumulated_surface_precipitation

  subroutine get_cloud_cover(state, s_field, position, p_field, surface_element_list, bc_path, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    type(vector_field), intent(in) :: position
    character(len=OPTION_PATH_LEN) :: bc_path
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    integer :: i, rstat, cstat, nits
    real :: threshold, alpha, height, top
    type(scalar_field) :: local, localmean
    type(scalar_field), pointer :: QR, QC, tmp
    
    call allocate(local,p_field%mesh,"Local")
    call allocate(localmean,p_field%mesh,"LocalMean")
    call zero(localmean)

    QR=>extract_scalar_field(state, "Qrain", rstat)
    QC=>extract_scalar_field(state, "Qdrop", cstat)
    
    stat=0
    if (cstat/=0 .or. rstat/=0) stat=-1
    if (stat == 0) then      
      call set(local,QC)
      call addto(local,QR)
      
      call get_option(trim(bc_path)//'/number_of_iterations',nits)
      call get_option(trim(bc_path)//'/alpha',alpha)
      call get_option(trim(bc_path)//'/domain_top',top)
      call get_option(trim(bc_path)//'/limit',threshold)
      
      call calculate_horizontal_filter(state, local, localmean, nits=nits, alpha=alpha, &
         base_path=p_field%option_path, vertical=.true.)
    
      call remap_field_to_surface(localmean, s_field, surface_element_list)
      
      do i =1, node_count(s_field)
        height=top-position%val(position%dim,i)
        if (s_field%val(i) < threshold*height) then
	  s_field%val(i) = 0.0
	else
	  s_field%val(i) = 1.0
	endif
      enddo
    endif
    
    tmp=>extract_scalar_field(state, "Temporary", stat=stat)
    if (stat==0) call set(tmp,localmean)
    
    call deallocate(local)    
    call deallocate(localmean)
    
  end subroutine get_cloud_cover

  subroutine get_liquid_water_path(state, s_field, position, p_field, surface_element_list, bc_path, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: p_field
    type(vector_field), intent(inout) :: position
    character(len=OPTION_PATH_LEN) :: bc_path
    integer, dimension(:), pointer, intent(in) :: surface_element_list    
    integer, intent(out), optional :: stat 
    
    integer :: i, rstat, cstat, nits, dim
    real :: ql, ql_old, alpha, top
    type(scalar_field) :: local_qc, local_qr, localmean, density_remap
    type(scalar_field), pointer :: QR, QC, density, sfield
    type(vector_field), pointer :: X

    QR=>extract_scalar_field(state, "Qrain", rstat)
    QC=>extract_scalar_field(state, "Qdrop", cstat)
    density=>extract_scalar_field(state, "Density")
    X=>extract_vector_field(state, "Coordinate")
    
    call allocate(density_remap,X%mesh,"RemapDensity")
    call allocate(local_qc,X%mesh,"LocalQC")
    call allocate(local_qr,X%mesh,"LocalQR")
    call allocate(localmean,X%mesh,"LocalMean")    
    call safe_set(state,density_remap,density)

    stat=0
    if (cstat/=0 .or. rstat/=0) stat=-1
    if (stat == 0) then
    
      call safe_set(state,local_qc,QC)
      call safe_set(state,local_qr,QR)
      call addto(local_qc,local_qr)
      call scale(local_qc,density_remap)
      
      call get_option(trim(bc_path)//'/number_of_iterations',nits)
      call get_option(trim(bc_path)//'/alpha',alpha)
      call get_option(trim(bc_path)//'/domain_top',top)
      
      call calculate_horizontal_filter(state, local_qc, localmean, nits=nits, alpha=alpha, &
         base_path=p_field%option_path, vertical=.true.)   
      
      dim=X%dim
      do i =1, node_count(s_field)
        s_field%val(i)=(top-position%val(dim,i))*localmean%val(i)
      enddo

    endif
    
    call deallocate(localmean)
    call deallocate(local_qc)
    call deallocate(local_qr)    
    call deallocate(density_remap)   
    
  end subroutine get_liquid_water_path

  subroutine write_surface(dump_no, state)
  
    integer, intent(in) :: dump_no
    type(state_type), dimension(:), intent(inout) :: state
    
    type(vector_field) :: surface_coordinate
    type(vector_field), pointer :: model_coordinate
    type(mesh_type), pointer :: model_mesh, surface_mesh

    type(scalar_field), pointer :: pfield
    type(scalar_field), dimension(:), allocatable :: lsfields
    character(len = OPTION_PATH_LEN) :: filename, bc_path, dump_format
    character(len = FIELD_NAME_LEN) :: field_name, mesh_name, bc_name, parent_name
    integer, dimension(:), pointer:: surface_element_list    
    integer :: i, f, counter, max_dump_no, stat, index
    logical :: multi_state, write_region_ids=.false.
    
    ewrite(1, *) "In write_surface"

    call get_option("/simulation_name", filename)
    call get_option("/io/max_dump_file_count", max_dump_no, stat, default = huge(0))
    
    index = modulo(dump_no, max_dump_no)
        
    call get_option("/io/output_mesh[0]/name", mesh_name)
    model_mesh => extract_mesh(state(1), mesh_name)
    
    call get_option("/io/dump_format", dump_format)
    select case(trim(dump_format))
    case("vtk")
    
    ewrite(2, *) "Writing output " // int2str(index) // " to vtu"
    
    multi_state = size(state) > 1
        
    ! count number of scalar fields in output:
    counter = 0
    do i = 1, size(state)
      if (associated(state(i)%scalar_fields)) then
        do f = 1, option_count("/material_phase["//int2str(i-1)//"]/scalar_field")
          call get_option("/material_phase["//int2str(i-1)//"]/scalar_field["//int2str(f-1)//"]/name",field_name)
	  if (have_option("/material_phase["//int2str(i-1)//"]/scalar_field::" &
	      //trim(field_name)//"/diagnostic/boundary_conditions::diagnostic") .and. &
	      include_scalar_field_in_vtu(state, i, field_name)) then
            counter = counter + 1
          end if
        end do
      end if
    end do
    
    ! collect scalar fields:
    allocate(lsfields(1:counter))
    counter = 0
    do i = 1, size(state)
      if (associated(state(i)%scalar_fields)) then
        do f = 1, option_count("/material_phase["//int2str(i-1)//"]/scalar_field")
          call get_option("/material_phase["//int2str(i-1)//"]/scalar_field["//int2str(f-1)//"]/name",field_name)
	  if (have_option("/material_phase["//int2str(i-1)//"]/scalar_field::" &
	      //trim(field_name)//"/diagnostic/boundary_conditions::diagnostic") .and. &
	      include_scalar_field_in_vtu(state, i, field_name)) then
            counter = counter + 1
      
            bc_path='/material_phase['//int2str(i-1)//']/scalar_field::'//trim(field_name)//'/diagnostic/boundary_conditions'
            call get_option(trim(bc_path)//"::diagnostic/parent_field_name", parent_name)
            call get_option(trim(bc_path)//"::diagnostic/type[0]/name", bc_name)
    
            pfield=>extract_scalar_field(state(i), trim(parent_name))
            lsfields(counter)=extract_surface_field(pfield, trim(bc_name), name=trim(bc_name))
            if (multi_state) then
              lsfields(counter)%name = trim(state(i)%name)//'::'//trim(field_name)
            end if
          end if
        end do
      end if
    end do
    
    if (counter > 0) then
      call get_boundary_condition(pfield, trim(bc_name), surface_mesh=surface_mesh, &
    	   surface_element_list=surface_element_list)
    
      model_coordinate=>get_external_coordinate_field(state(1), model_mesh)
      surface_coordinate=get_coordinates_remapped_to_surface(model_coordinate,surface_mesh,surface_element_list) 
    
      if (counter < 1 .or. node_count(surface_coordinate) < 1) return
      ewrite(2, *) "Writing using mesh " // trim(mesh_name)
      ewrite(2, "(a,i0,a)") "Writing ", size(lsfields), " surface field(s)"
    
      call vtk_write_fields('surface_'//filename, &
           index=index, &
           position=surface_coordinate, &
           model=surface_mesh,  &
           sfields=lsfields, &
           write_region_ids=write_region_ids)
    endif
    
    case default
    FLAbort("Unrecognised dump file format.")      
    end select
         
    ewrite(1, *) "Exiting write_surface"
    
  end subroutine write_surface
          
end module diagnostic_fields_surface
