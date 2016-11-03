!    Copyright (C) 2008 Imperial College London and others.
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
!    License as published by the Free Software Foundation
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

module les_module
  !!< This module contains several subroutines and functions used to implement LES models
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use vector_tools
  use fetools
  use field_derivatives
  use equation_of_state
  use fields
  use state_module
  use field_options
  use solvers
  use smoothing_module
  use state_fields_module
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
  public les_init_diagnostic_fields, les_assemble_diagnostic_fields, les_solve_diagnostic_fields, &
         leonard_tensor, les_strain_rate, tke_sgs, smagorinsky_2nd_sgs, smagorinsky_4th_sgs, non_linear_sgs, &
	 wale_sgs, dynamic_smagorinsky_sgs, calculate_richardson, get_les_options
  public smagorinsky_coefficient, backscatter_coefficient, prandtl_number, &
         length_scale_type, buoyancy_correction, les_second_order, les_tke, &
	 les_non_linear, les_fourth_order, wale, dynamic_les, use_hydrostatic

  
  logical :: buoyancy_correction, use_hydrostatic, les_tke, les_second_order, les_non_linear, les_fourth_order, wale, dynamic_les
  real :: smagorinsky_coefficient, backscatter_coefficient, prandtl_number
  character(len=OPTION_PATH_LEN) :: length_scale_type, correction_model
  
  real, parameter :: kappa=0.4 !vonKarman constant

contains

  subroutine get_les_options (les_option_path, momentum)
  
    implicit none
    logical :: momentum
    character(len=OPTION_PATH_LEN), intent(in) :: les_option_path

       buoyancy_correction=.false.; les_tke=.false.; les_second_order=.false.; les_non_linear=.false.; les_fourth_order=.false.; wale=.false.; dynamic_les=.false.
       les_second_order=have_option(trim(les_option_path)//"/second_order")
       les_fourth_order=have_option(trim(les_option_path)//"/fourth_order")
       les_non_linear=have_option(trim(les_option_path)//"/non_linear")
       wale=have_option(trim(les_option_path)//"/wale")
       dynamic_les=have_option(trim(les_option_path)//"/dynamic_les")
       les_tke=have_option(trim(les_option_path)//"/tke")

       if (les_second_order) then
    	  call get_option(trim(les_option_path)//"/second_order/smagorinsky_coefficient", &
    	       smagorinsky_coefficient)
    	  call get_option(trim(les_option_path)//"/second_order/Prandtl_turbulent", &
    	       prandtl_number)
          call get_option(trim(les_option_path)//"/second_order/length_scale_type", &
	       length_scale_type)
  	  buoyancy_correction=have_option(trim(les_option_path)//"/second_order/buoyancy_correction")
	  if (buoyancy_correction) then
	     call get_option(trim(les_option_path)//"/second_order/buoyancy_correction/name",correction_model)
	     use_hydrostatic=have_option(trim(les_option_path)//"/second_order/buoyancy_correction::"//trim(correction_model)//"/use_hydrostatic_state")
	  endif
       
       else if (les_fourth_order) then
    	  call get_option(trim(les_option_path)//"/fourth_order/smagorinsky_coefficient", &
    	       smagorinsky_coefficient)

       else if (les_non_linear) then
    	  call get_option(trim(les_option_path)//"/non_linear/backscatter_coefficient", &
    	       backscatter_coefficient)
    	  call get_option(trim(les_option_path)//"/non_linear/Prandtl_turbulent", &
    	       prandtl_number)
          call get_option(trim(les_option_path)//"/non_linear/length_scale_type", &
	       length_scale_type)
          smagorinsky_coefficient=backscatter_coefficient
	  
       else if (wale) then
    	  call get_option(trim(les_option_path)//"/wale/smagorinsky_coefficient", &
    	       smagorinsky_coefficient)
       
       else if(dynamic_les) then
          call get_option(trim(les_option_path)//"/dynamic_les/length_scale_type", &
	       length_scale_type)
       
       else if(les_tke) then
    	  call get_option(trim(les_option_path)//"/non_linear/Prandtl_turbulent", &
    	       prandtl_number)
          call get_option(trim(les_option_path)//"/dynamic_les/length_scale_type", &
	       length_scale_type)
       end if        
  
  end subroutine get_les_options

  subroutine les_init_diagnostic_fields(state, have_eddy_viscosity, have_non_linear, &
                           have_filter_width, have_coeff, have_richardson)

    ! Arguments
    type(state_type), intent(inout)             :: state
    logical, intent(in)                         :: have_filter_width, have_coeff, have_richardson, have_eddy_viscosity, have_non_linear
    
    ! Local variables
    logical, dimension(3)                       :: have_diagnostic_tfield
    logical, dimension(2)                       :: have_diagnostic_sfield
    character(len=FIELD_NAME_LEN), dimension(3) :: diagnostic_tfield_names
    character(len=FIELD_NAME_LEN), dimension(2) :: diagnostic_sfield_names
    type(tensor_field), pointer                 :: tfield
    type(scalar_field), pointer                 :: sfield
    integer                                     :: i, stat

    ewrite(2,*) "Initialising optional LES diagnostic fields"
    
    have_diagnostic_tfield = (/have_eddy_viscosity, have_non_linear, have_filter_width/)
    diagnostic_tfield_names(1) = "EddyViscosity"
    diagnostic_tfield_names(2) = "SFSLeonard"
    diagnostic_tfield_names(3) = "FilterWidth"
    
    diagnostic_tfield_loop: do i = 1, size(diagnostic_tfield_names)
      if(have_diagnostic_tfield(i)) then
         tfield => extract_tensor_field(state, diagnostic_tfield_names(i), stat=stat)
         if (stat==0) call zero(tfield)
      end if
    end do diagnostic_tfield_loop

    have_diagnostic_sfield = (/have_coeff, have_richardson/)
    diagnostic_sfield_names(1) = "SmagorinskyCoefficient"
    diagnostic_sfield_names(2) = "RichardsonNumber"

    diagnostic_sfield_loop: do i = 1, size(diagnostic_sfield_names)
      if(have_diagnostic_sfield(i)) then
         sfield => extract_scalar_field(state, diagnostic_sfield_names(i))
         call zero(sfield)
      end if
    end do diagnostic_sfield_loop

  end subroutine les_init_diagnostic_fields

  subroutine les_assemble_diagnostic_fields(state, X, u, Rho, viscosity, leonard, ele, &
                 mesh_size_gi, les_tensor_gi, les_leonard_gi, les_coef_gi, &
                 have_sfs_leonard, have_filter_width, have_coeff)

    ! Arguments
    type(state_type), intent(inout)                             :: state
    type(vector_field), intent(in)                              :: u, X
    type(scalar_field), intent(in)				:: Rho
    type(tensor_field), intent(inout)				:: viscosity, leonard
    integer, intent(in)                                         :: ele
    real, dimension(ele_ngi(u,ele)), intent(in)                 :: les_coef_gi
    real, dimension(u%dim,u%dim,ele_ngi(u,ele)), intent(in)     :: mesh_size_gi, les_tensor_gi, les_leonard_gi
    logical, intent(in) :: have_filter_width, have_coeff, have_sfs_leonard
    
    ! Local variables
    type(element_type), pointer 				:: shape
    type(tensor_field), pointer                                 :: tfield
    type(scalar_field), pointer                                 :: sfield
    real, dimension(u%dim,u%dim,ele_loc(u,ele))                 :: tensor_loc
    real, dimension(ele_loc(u,ele))                             :: scalar_loc
    real, dimension(ele_ngi(Rho,ele))				:: Rho_q
    real, dimension(ele_ngi(u,ele))	                        :: detwei
    
    Rho_q=ele_val_at_quad(Rho, ele)

    ! Eddy viscosity
    shape=>ele_shape(viscosity,ele)
    call transform_to_physical(X, ele, detwei=detwei)
    call set(viscosity, ele_nodes(viscosity, ele), shape_tensor_rhs(shape, les_tensor_gi, detwei))

    ! SFS strain
    if (have_sfs_leonard) then
       shape=>ele_shape(leonard,ele)
       call transform_to_physical(X, ele, detwei=detwei)
       call set(leonard, ele_nodes(leonard, ele), shape_tensor_rhs(shape, les_leonard_gi, detwei))
    end if

    ! Filter width
    if(have_filter_width) then
       tfield=>extract_tensor_field(state, "FilterWidth")
       shape=>ele_shape(tfield,ele)
       call transform_to_physical(X, ele, detwei=detwei)
       call set(tfield, ele_nodes(tfield, ele), shape_tensor_rhs(shape, mesh_size_gi, detwei))
    end if

    ! Smagorinsky Coefficient
    if(have_coeff) then
       sfield=>extract_scalar_field(state, "SmagorinskyCoefficient")
       shape=>ele_shape(sfield,ele)
       call transform_to_physical(X, ele, detwei=detwei)
       call set(sfield, ele_nodes(sfield, ele), shape_rhs(shape, les_coef_gi*detwei))
    end if

  end subroutine les_assemble_diagnostic_fields
  
  subroutine les_solve_diagnostic_fields(state, have_filter_width, have_coeff)

    ! Arguments
    type(state_type), intent(inout) :: state
    logical, intent(in) :: have_filter_width, have_coeff
    
    ! Local variables
    logical, dimension(2)                       :: have_diagnostic_tfield
    logical, dimension(1)                       :: have_diagnostic_sfield
    character(len=FIELD_NAME_LEN), dimension(2) :: diagnostic_tfield_names
    character(len=FIELD_NAME_LEN), dimension(1) :: diagnostic_sfield_names
    type(tensor_field), pointer                 :: tfield
    type(scalar_field), pointer                 :: sfield
    integer                                     :: i
    type(vector_field), pointer                 :: u
    type(csr_matrix), pointer                   :: mass_matrix
    type(scalar_field), pointer                 :: lumped_mass
    type(scalar_field)                          :: inv_lumped_mass
    logical                                     :: lump_mass = .false.
    logical                                     :: use_submesh = .false.
    
    ewrite(2,*) "Solving for optional LES diagnostic fields"
        
    u => extract_vector_field(state, "Velocity")
    
    have_diagnostic_tfield = (/.true., have_filter_width/)
    diagnostic_tfield_names(1) = "EddyViscosity"
    diagnostic_tfield_names(2) = "FilterWidth"
    
    diagnostic_tfield_loop: do i = 1, size(diagnostic_tfield_names)
      if(have_diagnostic_tfield(i)) then
         tfield => extract_tensor_field(state, diagnostic_tfield_names(i))
         lump_mass = have_option(trim(tfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix")
         use_submesh = have_option(trim(tfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix/use_submesh") ! For P2 meshes.
            
         if(lump_mass) then
            if(use_submesh) then
               lumped_mass => get_lumped_mass_on_submesh(state, tfield%mesh)
            else
               lumped_mass => get_lumped_mass(state, tfield%mesh)
            end if
            call allocate(inv_lumped_mass, tfield%mesh)
            call invert(lumped_mass, inv_lumped_mass)
            call scale(tfield, inv_lumped_mass)
            call deallocate(inv_lumped_mass)
         else
            mass_matrix => get_mass_matrix(state, tfield%mesh)
            call petsc_solve(tfield, mass_matrix, tfield, option_path=u%option_path)
         end if
      end if
    end do diagnostic_tfield_loop
    
    have_diagnostic_sfield = (/have_coeff/)
    diagnostic_sfield_names(1) = "SmagorinskyCoefficient"
    
    diagnostic_sfield_loop: do i = 1, size(diagnostic_sfield_names)
      if(have_diagnostic_sfield(i)) then
         sfield => extract_scalar_field(state, diagnostic_sfield_names(i))
         lump_mass = have_option(trim(sfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix")
         use_submesh = have_option(trim(sfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix/use_submesh") ! For P2 meshes.
            
         if(lump_mass) then
            if(use_submesh) then
               lumped_mass => get_lumped_mass_on_submesh(state, sfield%mesh)
            else
               lumped_mass => get_lumped_mass(state, sfield%mesh)
            end if
            call allocate(inv_lumped_mass, sfield%mesh)
            call invert(lumped_mass, inv_lumped_mass)
            call scale(sfield, inv_lumped_mass)
            call deallocate(inv_lumped_mass)
         else
            mass_matrix => get_mass_matrix(state, sfield%mesh)
            call petsc_solve(sfield, mass_matrix, sfield, option_path=u%option_path)
         end if
      end if
    end do diagnostic_sfield_loop

  end subroutine les_solve_diagnostic_fields
  
  subroutine leonard_tensor(u, positions, fnu, tnu, leonard, strainprod, alpha, gamma, path)

    ! Unfiltered velocity
    type(vector_field), pointer                           :: u
    type(vector_field), intent(in)                        :: positions
    ! Filtered velocities
    type(vector_field), pointer                           :: fnu, tnu
    ! Leonard tensor and strain product
    type(tensor_field), pointer                           :: leonard, strainprod
    ! Scale factors
    real, intent(in)                                      :: alpha, gamma
    character(len=OPTION_PATH_LEN), intent(in)            :: path
    ! Local quantities
    type(vector_field), dimension(positions%dim)	  :: nu_grad
    type(tensor_field), pointer                           :: ui_uj, tui_tuj
    character(len=OPTION_PATH_LEN)                        :: lpath
    integer                                               :: i, d,ele, gi
    real, dimension(:), allocatable                       :: u_loc
    real, dimension(:,:), allocatable                     :: t_loc
    real, dimension(ele_loc(u,1), ele_ngi(u,1), u%dim)    :: du_t
    real, dimension(ele_ngi(u,1))                         :: detwei
    real, dimension(u%dim, u%dim, ele_ngi(u,1))           :: strain_gi, strain_prod_gi
    type(element_type)                                    :: shape_nu
    real, dimension(positions%dim, ele_ngi(u,1))		  :: nu_ele
    real, dimension(positions%dim, positions%dim, ele_ngi(u,1))	  :: nu_grad_ele 

    ! Path is to level above solver options
    lpath = (trim(path)//"/dynamic_les")
    ewrite(2,*) "filter factor alpha: ", alpha
    ewrite(2,*) "filter factor gamma: ", gamma

    ! First filter operator returns u^f:
    call anisotropic_smooth_vector(u, positions, fnu, alpha, lpath)
    ! Test filter operator needs the ratio of test filter to mesh size and returns u^ft:
    call anisotropic_smooth_vector(fnu, positions, tnu, alpha*gamma, lpath)
    ewrite_minmax(u)
    ewrite_minmax(fnu)
    ewrite_minmax(tnu)

    ! Velocity products (ui*uj)
    allocate(ui_uj); allocate(tui_tuj)
    call allocate(ui_uj, u%mesh, "NonlinearVelocityProduct")
    call allocate(tui_tuj, u%mesh, "TestNonlinearVelocityProduct")
    call zero(ui_uj); call zero(tui_tuj)

    ! Other local variables
    allocate(u_loc(u%dim)); allocate(t_loc(u%dim, u%dim))
    u_loc=0.0; t_loc=0.0

    ! Get cross products of velocities
    do i=1, node_count(u)
      u_loc = node_val(fnu,i)
      t_loc = outer_product(u_loc, u_loc)
      call set( ui_uj, i, t_loc )
      u_loc = node_val(tnu,i)
      ! Calculate (test-filtered velocity) products: (ui^ft*uj^ft)
      t_loc = outer_product(u_loc, u_loc)
      call set( tui_tuj, i, t_loc )
    end do

    ! Calculate test-filtered (velocity products): (ui^f*uj^f)^t
    call anisotropic_smooth_tensor(ui_uj, positions, leonard, alpha*gamma, lpath)

    ! Leonard tensor field
    call addto( leonard, tui_tuj, -1.0 )

    ! Zero tensor field for reuse in strain product assembly
    call zero(ui_uj)
    
    ! Calculate velocity gradient
    call grad(u, positions, nu_grad)

    do i=1, element_count(u)
      shape_nu = ele_shape(u, i)
      nu_ele = ele_val_at_quad(u, i)
      do d=1,u%dim
        nu_grad_ele(i,:,:) = ele_val_at_quad(nu_grad(d), i)
      enddo
      ! Assuming no FE stabilisation is used with LES so we can use velocity shape.
      call transform_to_physical(positions, i, shape_nu, dshape=du_t, detwei=detwei)
      ! Strain rate of first filtered velocity S1^f
      strain_gi = les_strain_rate(nu_ele, nu_grad_ele)
      do gi=1, ele_ngi(u, ele)
        ! Strain product = strain modulus*strain rate: |S1^f|S1^f
        strain_prod_gi(:,:,gi) = sqrt(2*sum(strain_gi(:,:,gi)*strain_gi(:,:,gi))) * strain_gi(:,:,gi)
      end do
      ! Assemble local tensor field
      call addto(ui_uj, ele_nodes(u,i), shape_tensor_rhs(shape_nu, strain_prod_gi, detwei))
    end do

    ! Filter strain product with test filter: (|S1^f|S1^f)^t
    call anisotropic_smooth_tensor(ui_uj, positions, strainprod, alpha*gamma, lpath)

    ! Deallocates
    deallocate(u_loc, t_loc)
    call deallocate(ui_uj)
    call deallocate(tui_tuj)
    deallocate(ui_uj); deallocate(tui_tuj)

  end subroutine leonard_tensor

  function les_strain_rate(nu_ele, nu_grad_ele)
    !! Computes the strain rate
    real, dimension(:,:) :: nu_ele
    real, dimension(:,:,:) :: nu_grad_ele
    
    real, dimension( size(nu_ele,1),size(nu_ele,1),size(nu_ele,2) ):: les_strain_rate
    real, dimension( size(nu_ele,1),size(nu_ele,1) ):: s
    integer dim, ngi, gi, di

    ngi=size(nu_ele,2)
    dim=size(nu_ele,1)
       
    do gi=1, ngi
    
      s=0.5*nu_grad_ele(:,:,gi)
      les_strain_rate(:,:,gi)=s+transpose(s)
      
    enddo
    
  end function les_strain_rate

  function les_viscosity_strength(ele, nu_ele, nu_grad_ele, correction, ri, pr)
    !! Computes the strain rate modulus for the LES model
    integer, intent(in) :: ele
    real, dimension(:,:) :: nu_ele
    real, dimension(:,:,:) :: nu_grad_ele
    !! Richardson number
    real, dimension(:), intent(in), optional :: ri
    !! Prandtl number
    real, intent(in), optional :: pr
    character(len=OPTION_PATH_LEN), intent(in), optional :: correction

    real, dimension(size(nu_ele,1),size(nu_ele,1)) :: s
    real, dimension(size(nu_ele,2)) :: les_viscosity_strength

    real vis, corr, sii
    integer dim, ngi, gi, di, i
    real, parameter :: b=0.5, c=4., r=4., a=16., d=1.15, e=1.5, f=1.2, g=0.0623

    ngi=size(nu_ele,2)
    dim=size(nu_ele,1)
    les_viscosity_strength=0.0
    
    do gi=1,ngi

       s=0.5*nu_grad_ele(:,:,gi)
       s=s+transpose(s)
       vis=sqrt( 2.*sum( s**2 ) )
       
       sii=0.0
       do i = 1, size(nu_ele,1)
         sii=sii+s(i,i)
       enddo
       sii=abs(sii)
       
       ! Buoyancy correction
       if (present(ri).and.present(correction)) then
         assert(size(nu_ele,2) == size(ri))
         corr=1.
	 select case (trim(correction))
	   case("Lilly") 
	     if (ri(gi) <= 0.25 .and. ri(gi) >= 0.) then
	       corr=sqrt(1. - ri(gi)/pr)
	     else if (ri(gi) <= 0.) then
	       corr=0.
	     endif
	   case("Brown")
	     if (ri(gi) <= 0.) then
	       corr=(1. - a*max(ri(gi)/pr,-0.0625))**b
	     else
	       corr=(1. - c*min(ri(gi)/pr,0.25))**r
	     endif
	   case("Cheng")
	     if (ri(gi) <= 0.25 .and. ri(gi) >= 0.) then
	       corr=(1. + f*ri(gi)/pr)**(-b)
	     else if (ri(gi) > 0.25) then
	       corr=(1. + g*min(ri(gi)/pr,10.))
	     else if (ri(gi) <= 0.) then
	       corr=0.
	     endif
	     corr=corr*(1. + d*(sqrt( 2.*sum( s**2 ) )/sii)**2.)**(-e)
	   case default
	     corr=1.
	 end select
	 vis=vis*corr
       endif

       les_viscosity_strength(gi)=vis

    end do

  end function les_viscosity_strength
  
  function wale_viscosity_strength(nu_ele, nu_grad_ele)
    !! Computes the traceless symmetric part of the square of
    !! the resolved velocity gradient tensor for the LES model
    !! See a WALE paper for more (G_{ij})
    real, dimension(:,:) :: nu_ele
    real, dimension(:,:,:) :: nu_grad_ele

    real, dimension( size(nu_ele,1) ):: wale_viscosity_strength

    real, dimension( size(nu_ele,1),size(nu_ele,1) ):: s
    real, dimension( size(nu_ele,1),size(nu_ele,1) ):: g
    real vis
    integer dim, ngi, gi, di, i

    ngi=size(nu_ele,2)
    dim=size(nu_ele,1)

    do gi=1, ngi

       s=nu_grad_ele(:,:,gi)
       g=0.5*matmul(s,s)
       g=g+transpose(g)
       forall(i=1:dim) g(i,i)=0.
       
       vis=sqrt( 2.*sum( g**2 ) )

       wale_viscosity_strength(gi)=vis

    end do

  end function wale_viscosity_strength
  
  subroutine calculate_richardson (state, ele, X, relu, grad_u, pottem, grad_pt, gravity, richardson)
  
    integer, intent(in) :: ele
    type(state_type), intent(inout) :: state
    type(vector_field), intent(in) :: X, relu
    type(scalar_field), intent(in) :: pottem
    type(vector_field), intent(in) :: grad_pt
    type(vector_field), dimension(X%dim) :: grad_u
    type(scalar_field), intent(inout) :: richardson
    real, intent(in) :: gravity
    
    integer :: dim, gi, ngi, i
    type(element_type), pointer :: ri_shape
    real, dimension(ele_ngi(pottem,ele)) :: pottem_ele
    real, dimension(X%dim, ele_ngi(pottem,ele)) :: grad_pt_ele
    real, dimension(ele_ngi(relu,ele)) :: s2, n2, detwei
    real, dimension(X%dim, X%dim, ele_ngi(relu,ele)) :: grad_u_ele
    real, dimension(ele_loc(relu,ele), ele_ngi(relu,ele), x%dim) :: dri
    type(scalar_field), pointer :: sfield
    logical, save :: test=.true.

    dim=relu%dim
    ngi=ele_ngi(relu, ele)
    assert(ngi == ele_ngi(richardson,ele))
    assert(ngi == ele_ngi(pottem,ele))
    
    ri_shape=>ele_shape(richardson,ele)
    call transform_to_physical(X, ele, ri_shape, dshape=dri, detwei=detwei)
    
    pottem_ele = ele_val_at_quad(pottem, ele)
    do i=1,dim
      grad_u_ele(i,:,:)=ele_val_at_quad(grad_u(i), ele)
    enddo
    grad_pt_ele=ele_val_at_quad(grad_pt, ele)
    
    ! Calculate denominator:
    do i=1,dim-1
      s2=s2+grad_u_ele(i,dim,:)**2.
    enddo
    s2=sign(1.,s2)*max(abs(s2),1.e-10)
    
    ! Calculate Buoyancy fequency
    n2 = gravity * grad_pt_ele(dim,:)/pottem_ele
    
    ! set Richardson
    call set(richardson, ele_nodes(richardson, ele), shape_rhs(ri_shape, n2/s2*detwei))
    
    ! Richardson number
    if (has_scalar_field(state,"RichardsonNumber")) then
      sfield=>extract_scalar_field(state, "RichardsonNumber")
      call set(sfield, ele_nodes(sfield, ele), shape_rhs(ri_shape, n2/s2*detwei))
    end if
            
  end subroutine calculate_richardson
  
  subroutine tke_sgs (ele, g, x, t, u, nu_grad, du_t, rho, hb, richardson, tke, les_tensor_gi, tke_source_gi)
  
    integer, intent(in) :: ele
    real, intent(in)    :: g
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: x
    type(scalar_field), intent(in) :: t, rho, hb, richardson, tke
    type(vector_field), dimension(x%dim), intent(in) :: nu_grad
    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), mesh_dim(u)), intent(in) :: du_t
    
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)), intent(inout) :: les_tensor_gi
    real, dimension(ele_ngi(u, ele)), optional, intent(inout) :: tke_source_gi
    
    integer					   :: dim, ngi, gi, i
    real					   :: length_scale, coefficient_e, coefficient_eps, scale, pow, s2
    real, dimension(ele_ngi(u, ele))	 	   :: rho_gi,hb_gi,ri_gi,tke_gi,tke_prod_gi,tke_diss_gi, tke_buoy_gi
    real, dimension(x%dim,ele_ngi(u, ele))         :: nu_ele
    real, dimension(x%dim,x%dim,ele_ngi(u, ele))   :: nu_grad_ele, les_strain_gi 
    real, dimension(x%dim,x%dim)                   :: s
    
    ngi=ele_ngi(u, ele)
    dim=u%dim
    coefficient_e=0.15
    coefficient_eps=0.7
  
    ri_gi=ele_val(richardson, ele)    
    tke_gi=ele_val(tke, ele) 
    rho_gi=ele_val(Rho, ele) 
    hb_gi=ele_val(hb, ele) 
       
    nu_ele=ele_val_at_quad(u, ele)
    do i = 1,dim
      nu_grad_ele(i,:,:)=ele_val_at_quad(nu_grad(i), ele)
    enddo
    
    select case(length_scale_type)
       case("scalar")
	  ! Length scale is the cube root of the element's volume in 3D.
	  ! In 2D, it is the square root of the element's area.
	  do gi = 1, ngi   
	     les_tensor_gi(:,:,gi)=2.*les_tensor_gi(:,:,gi)
             s=0.5*nu_grad_ele(:,:,gi)
             les_strain_gi(:,:,gi)=s+transpose(s)
	     s2=2.*sum( les_strain_gi(:,:,gi)*les_strain_gi(:,:,gi) )
	     
	     length_scale=tke_length_scale(ele, x, tke_gi(gi), ri_gi(gi), s2)
	     les_tensor_gi(:,:,gi)=coefficient_e*length_scale*sqrt(tke_gi(gi))
	     
	     if (trim(t%name) == "TurbulentKineticEnergy" .and. present(tke_source_gi)) then	       
	       tke_prod_gi(gi)=les_tensor_gi(1,1,gi)*sum( les_strain_gi(:,:,gi)*nu_grad_ele(:,:,gi) ) 
	       tke_diss_gi(gi)=coefficient_eps/length_scale*tke_gi(gi)**(3./2.)
	       tke_buoy_gi(gi)=g*nu_ele(dim,gi)*(rho_gi(gi) - hb_gi(gi))
	       
	       tke_source_gi(gi)=tke_prod_gi(gi)+tke_diss_gi(gi)+tke_buoy_gi(gi)
	     endif
	  end do
	  
       case default
	  FLExit("Unknown length scale type")
    end select  

contains

    real function tke_length_scale(ele, x, tke_gi, ri_gi, s2)
      integer :: ele
      type(vector_field) :: x
      real :: tke_gi, ri_gi, s2
      real :: delta, deardorff, hunt
      character(len=OPTION_PATH_LEN) :: correction="Kosovic"
      
      delta=1./length_scale_scalar(x, ele)
      
      select case (trim(correction))
        case("Kosovic") 
          deardorff=max(0.577*tke_gi/(max(ri_gi*s2,1.e-6)),0.0) !0.577=0.76**2
          hunt=max(7.6176*tke_gi/max(s2,1.e-6),0.0)		!7.6167=2.76**2
      
          tke_length_scale=delta/sqrt(1. + delta**2./deardorff + delta**2./hunt)
	  
        case("Moeng")  
	  tke_length_scale=delta
          if (ri_gi <= 0.25) &
	      tke_length_scale=sqrt(max(0.577*tke_gi/(max(ri_gi*s2,1.e-6)),0.0))
	  
	case default
	  tke_length_scale=delta
      end select
      
      return
    end function
        
  end subroutine tke_sgs
    
  subroutine smagorinsky_2nd_sgs (ele, x, u, nu_grad, du_t, richardson, les_coef_gi, les_tensor_gi)
  
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: x
    type(scalar_field), intent(in) :: richardson
    type(vector_field), dimension(x%dim), intent(in) :: nu_grad
    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), mesh_dim(u)), intent(in) :: du_t
    
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)), intent(inout) :: les_tensor_gi
    real, dimension(ele_ngi(u, ele)), intent(inout)		  :: les_coef_gi
    
    integer					   :: dim, gi, i
    real, dimension(ele_ngi(u, ele))	 	   :: ri_gi
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)) :: length_scale_full
    real, dimension(x%dim, x%dim)		   :: length_scale
    real, dimension(x%dim,ele_ngi(u, ele))         :: nu_ele
    real, dimension(x%dim,x%dim,ele_ngi(u, ele))   :: nu_grad_ele 
    
    dim=x%dim
  
    ri_gi=ele_val_at_quad(richardson, ele)    
    nu_ele=ele_val_at_quad(u, ele)
    do i = 1, u%dim
      nu_grad_ele(i,:,:)=ele_val_at_quad(nu_grad(i), ele)
    enddo
    
    les_coef_gi = les_viscosity_strength(ele, nu_ele, nu_grad_ele, correction=correction_model, ri=ri_gi, pr=prandtl_number)
    
    select case(length_scale_type)
       case("scalar")
	  ! Length scale is the cube root of the element's volume in 3D.
	  ! In 2D, it is the square root of the element's area.
	  do gi = 1, size(les_coef_gi)
	     length_scale = 1./length_scale_scalar(x, ele) !+1./position_gi(dim, gi)
	     les_tensor_gi(:,:,gi) = 1.0/(length_scale)*(smagorinsky_coefficient**2)*les_coef_gi(gi)
	  end do
	  
       case("tensor")
	  ! This uses a tensor length scale metric from the adaptivity process
	  ! to better handle anisotropic elements.
	  length_scale_full = length_scale_tensor(du_t, ele_shape(u, ele))
	  do gi = 1, size(les_coef_gi)
	     length_scale = 1./length_scale_full(:,:,gi)
	     les_tensor_gi(:,:,gi) = 1.0/(length_scale)*(smagorinsky_coefficient**2)*les_coef_gi(gi)
	  end do
	  
       case default
	  FLExit("Unknown length scale type")
    end select  
  
  end subroutine smagorinsky_2nd_sgs
    
  subroutine non_linear_sgs (ele, x, u, nu_grad, du_t, les_tensor_gi, les_leonard_gi)
  
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: x
    type(vector_field), dimension(x%dim), intent(in) :: nu_grad
    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), mesh_dim(u)), intent(in) :: du_t
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)), intent(inout) :: les_tensor_gi, les_leonard_gi
    
    real, dimension(x%dim, ele_ngi(u, ele))	   :: nu_ele
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)) :: nu_grad_ele
    real, dimension(ele_ngi(u, ele))		   :: les_coef_gi
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)) :: les_nl_coef_gi
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)) :: length_scale_full
    real, dimension(x%dim, x%dim)		   :: length_scale
    real, dimension(ele_ngi(u, ele)) 		   :: isotropic_part
    
    integer :: dim, gi, i
    real :: c1, c2, cs
    
    dim=x%dim
    cs=sqrt((1.+backscatter_coefficient)/40.)
    c1=31.*backscatter_coefficient/(3.5*(1.+backscatter_coefficient))
    c2=c1
       
    nu_ele = ele_val_at_quad(u, ele)
    do i =1, u%dim
      nu_grad_ele(i,:,:) = ele_val_at_quad(nu_grad(i), ele)
    enddo
    
    ! Calculate linear part: smagorinsky eddy diffusivity model
    les_coef_gi = les_viscosity_strength(ele, nu_ele, nu_grad_ele)
    
    ! Calculate non-linear part
    les_nl_coef_gi = c1*non_linear_part(nu_ele, nu_grad_ele) + c2*rotational_part(nu_ele, nu_grad_ele)
    
    ! Multiply by length scale squared
    select case(length_scale_type)
       case("scalar")
	  ! Length scale is the cube root of the element's volume in 3D.
	  ! In 2D, it is the square root of the element's area.
	  do gi = 1, size(les_coef_gi)
	     length_scale = 1./length_scale_scalar(x, ele) 
	     les_tensor_gi(:,:,gi) = 1.0/(length_scale)*(cs**2)*les_coef_gi(gi)
	     les_leonard_gi(:,:,gi) = 1.0/(length_scale)*(cs**2)*les_nl_coef_gi(:,:,gi)
	  end do
	  
       case("tensor")
	  ! This uses a tensor length scale metric from the adaptivity process
	  ! to better handle anisotropic elements.
	  length_scale_full = length_scale_tensor(du_t, ele_shape(u, ele))
	  do gi = 1, size(les_coef_gi)
	     length_scale = 1./(length_scale_full(:,:,gi))
	     les_tensor_gi(:,:,gi) = 1.0/(length_scale)*(cs**2)*les_coef_gi(gi)
	     les_leonard_gi(:,:,gi) = 1.0/(length_scale)*(cs**2)*les_nl_coef_gi(:,:,gi)
	  end do
	  
       case default
	  FLExit("Unknown length scale type")
    end select
        
    contains
    
    function non_linear_part (nu_ele, nu_grad_ele) 
    
      real, dimension(:,:,:), intent(in) :: nu_grad_ele
      real, dimension(:,:), intent(in) :: nu_ele
      
      real, dimension(size(nu_ele,1),size(nu_ele,1),size(nu_ele,2)):: non_linear_part
      real, dimension(size(nu_ele,1),size(nu_ele,1)):: s
      integer :: ngi, dim, i, j
      real :: s2

      ngi=size(nu_ele,2)
      dim=size(nu_ele,1)
      do gi=1, ngi

         s=0.5*nu_grad_ele(:,:,gi)
         s=s+transpose(s)
         s2=sum( 2.*s*s )

	 do i = 1, dim
	   do j = 1, dim
             non_linear_part(i,j,gi) = sum( s(i,:)*s(:,j) )
	   enddo
	   non_linear_part(i,i,gi) = non_linear_part(i,i,gi) - s2/3.
	 enddo
	 
      enddo
      
    end function non_linear_part
    
    function rotational_part (nu_ele, nu_grad_ele) 
    
      real, dimension(:,:,:), intent(in) :: nu_grad_ele
      real, dimension(:,:), intent(in) :: nu_ele
      
      real, dimension(size(nu_ele,1),size(nu_ele,1),size(nu_ele,2)):: rotational_part
      real, dimension(size(nu_ele,1),size(nu_ele,1)):: s, r
      integer :: ngi, dim, i, j

      ngi=size(nu_ele,2)
      dim=size(nu_ele,1)
      
      do gi=1, ngi

         s=0.5*nu_grad_ele(:,:,gi)
	 r=s-transpose(s)
         s=s+transpose(s)

	 do i = 1, dim
	   do j = 1, dim
             rotational_part(i,j,gi) = sum(s(i,:)*r(:,j)) - sum(r(i,:)*s(:,j))
	   enddo
	 enddo	   
	 
      enddo
      
    end function rotational_part
  
  end subroutine non_linear_sgs
	    
  subroutine wale_sgs (ele, x, u, nu_grad, du_t, les_tensor_gi)

    integer, intent(in) :: ele
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: x
    type(vector_field), dimension(x%dim), intent(in) :: nu_grad
    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)), intent(inout) :: les_tensor_gi

    integer					   :: gi, i
    real, dimension(x%dim, ele_ngi(u, ele))	   :: nu_ele
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)) :: nu_grad_ele
    real, dimension(ele_ngi(u, ele))		   :: les_coef_gi, wale_coef_gi

    les_tensor_gi=length_scale_tensor(du_t, ele_shape(u, ele))
    nu_ele = ele_val_at_quad(u, ele)
    do i=1, u%dim
      nu_grad_ele(i,:,:) = ele_val_at_quad(nu_grad(i), ele)
    enddo
    
    les_coef_gi=les_viscosity_strength(ele, nu_ele, nu_grad_ele)
    wale_coef_gi=wale_viscosity_strength(nu_ele, nu_grad_ele)
    
    do gi=1, size(les_coef_gi)
    
      les_tensor_gi(:,:,gi)=4.*les_tensor_gi(:,:,gi)* &
            wale_coef_gi(gi)**3 * smagorinsky_coefficient**2 / &
            max(les_coef_gi(gi)**5 + wale_coef_gi(gi)**2.5, 1.e-10)
	    
    end do
  
  end subroutine wale_sgs
  
  subroutine smagorinsky_4th_sgs (ele, x, u, nu_grad, du_t, detwei, les_tensor_gi, rhs_addto)
  
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: x
    type(vector_field), dimension(x%dim), intent(in) :: nu_grad
    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
    real, dimension(ele_ngi(u, ele)), intent(in) :: detwei

    real, dimension(x%dim, x%dim, ele_ngi(u, ele)), intent(inout) :: les_tensor_gi
    real, dimension(u%dim, ele_loc(u, ele)), intent(inout)  :: rhs_addto

    integer :: dim, iloc
    real, dimension(x%dim, ele_loc(u,ele), ele_loc(u,ele)) :: div_les_viscosity
    real, dimension(ele_ngi(u, ele))		  :: les_coef_gi
    real, dimension(u%dim, ele_loc(u, ele))	  :: nu_ele
    real, dimension(u%dim, ele_ngi(u, ele))	  :: nu_q
    real, dimension(x%dim, x%dim, ele_loc(u,ele)) :: grad_u_nodes
    real, dimension(x%dim, x%dim, ele_ngi(u,ele)) :: grad_u_q
       
    les_tensor_gi=length_scale_tensor(du_t, ele_shape(u, ele))
    nu_ele = ele_val(u, ele)
    nu_q = ele_val_at_quad(u, ele)
    do dim=1, u%dim
      grad_u_nodes(dim,:,:)=ele_val(nu_grad(dim), ele)
      grad_u_q(dim,:,:)=ele_val_at_quad(nu_grad(dim), ele)
    enddo
    
    les_coef_gi=les_viscosity_strength(ele, nu_q, grad_u_q)
    div_les_viscosity=dshape_dot_tensor_shape(du_t, les_tensor_gi, ele_shape(u, ele), detwei)

    do dim=1, u%dim
       do iloc=1, ele_loc(u, ele)
          rhs_addto(dim,iloc)=rhs_addto(dim,iloc)+ &
               sum(div_les_viscosity(:,:,iloc)*grad_u_nodes(:,dim,:))
       end do
    end do  
  
  end subroutine smagorinsky_4th_sgs
  
  subroutine dynamic_smagorinsky_sgs (ele, x, u, fnu, tnu, gfnu, gtnu, du_t, strainprod, leonard, alpha, gamma, & 
                                  mesh_size_gi, les_coef_gi, les_tensor_gi)
  
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: x
    type(vector_field), intent(in) :: fnu, tnu
    type(vector_field), dimension(x%dim), intent(in) :: gfnu, gtnu
    real, intent(in)  		   :: alpha, gamma
    type(tensor_field), intent(in)    :: leonard, strainprod
    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
    
    real, dimension(x%dim, x%dim, ele_ngi(u, ele)), intent(inout) :: les_tensor_gi
    real, dimension(ele_ngi(u, ele)), intent(inout)		  :: les_coef_gi
  
    integer :: gi, i
    type(element_type)  	   :: shape_nu
    integer, dimension(:), pointer :: nodes_nu
    real, dimension(x%dim, ele_ngi(u, ele))        :: fnu_ele, tnu_ele
    real, dimension(x%dim, x%dim, ele_ngi(u,ele))  :: gfnu_ele, gtnu_ele
    real, dimension(x%dim, x%dim, ele_ngi(u,ele))  :: mesh_size_gi, leonard_gi
    real, dimension(x%dim, x%dim, ele_ngi(u,ele))  :: strain_gi, t_strain_gi, strainprod_gi
    real, dimension(x%dim, x%dim, ele_ngi(u,ele))  :: f_tensor, t_tensor
    real, dimension(ele_ngi(u, ele))		   :: f_scalar, t_scalar
    real, dimension(ele_ngi(u, ele))		   :: strain_mod, t_strain_mod
    real, dimension(x%dim, x%dim)                  :: mij

    shape_nu = ele_shape(u, ele)
    nodes_nu => ele_nodes(u, ele)
    
    fnu_ele=ele_val_at_quad(fnu, ele)
    tnu_ele=ele_val_at_quad(tnu, ele)
    do i=1, x%dim
      gfnu_ele(i,:,:)=ele_val_at_quad(gfnu(i), ele)
      gtnu_ele(i,:,:)=ele_val_at_quad(gtnu(i), ele)
    enddo

    ! Get strain S1 for first-filtered velocity
    strain_gi = les_strain_rate(fnu_ele, gfnu_ele)
    
    ! Get strain S2 for test-filtered velocity
    t_strain_gi = les_strain_rate(tnu_ele, gtnu_ele)
    
    ! Mesh size (untis length^2)
    mesh_size_gi = length_scale_tensor(du_t, ele_shape(u, ele))
    
    ! Leonard tensor and strain product at Gauss points
    leonard_gi = ele_val_at_quad(leonard, ele)
    strainprod_gi = ele_val_at_quad(strainprod, ele)

    do gi=1, ele_ngi(u, ele)
       ! Get strain modulus |S1| for first-filtered velocity
       strain_mod(gi) = sqrt( 2*sum(strain_gi(:,:,gi)*strain_gi(:,:,gi) ) )
       
       ! Get strain modulus |S2| for test-filtered velocity
       t_strain_mod(gi) = sqrt( 2*sum(t_strain_gi(:,:,gi)*t_strain_gi(:,:,gi) ) )
    end do

    select case(length_scale_type)
    
       case("scalar")
	  ! Scalar first filter width G1 = alpha^2*meshsize (units length^2)
	  f_scalar = alpha**2*length_scale_scalar(x, ele)
	  ! Combined width G2 = (1+gamma^2)*G1
	  t_scalar = (1.0+gamma**2)*f_scalar
	  
	  do gi=1, ele_ngi(u, ele)
	    ! Tensor M_ij = (|S2|*S2)G2 - ((|S1|S1)^f2)G1
	    mij = t_strain_mod(gi)*t_strain_gi(:,:,gi)*t_scalar(gi) - strainprod_gi(:,:,gi)*f_scalar(gi)
	    ! Model coeff C_S = -(L_ij M_ij) / 2(M_ij M_ij)
	    les_coef_gi(gi) = -0.5*sum(leonard_gi(:,:,gi)*mij) / sum(mij*mij)
	    ! Constrain C_S to be between 0 and 0.04.
	    les_coef_gi(gi) = min(max(les_coef_gi(gi),0.0), 0.04)
	    ! Isotropic tensor dynamic eddy viscosity = -2C_S|S1|.alpha^2.G1
	    les_tensor_gi(:,:,gi) = 2*alpha**2*les_coef_gi(gi)*strain_mod(gi)*f_scalar(gi)
	  end do
	  
       case("tensor")
	  ! First filter width G1 = alpha^2*mesh size (units length^2)
	  f_tensor = alpha**2*mesh_size_gi
	  ! Combined width G2 = (1+gamma^2)*G1
	  t_tensor = (1.0+gamma**2)*f_tensor
	  
	  do gi=1, ele_ngi(u, ele)
	    ! Tensor M_ij = (|S2|*S2).G2 - ((|S1|S1)^f2).G1
	    mij = t_strain_mod(gi)*t_strain_gi(:,:,gi)*t_tensor(:,:,gi) - strainprod_gi(:,:,gi)*f_tensor(:,:,gi)
	    ! Model coeff C_S = -(L_ij M_ij) / 2(M_ij M_ij)
	    les_coef_gi(gi) = -0.5*sum(leonard_gi(:,:,gi)*mij) / sum(mij*mij)
	    ! Constrain C_S to be between 0 and 0.04.
	    les_coef_gi(gi) = min(max(les_coef_gi(gi),0.0), 0.04)
	    ! Anisotropic tensor dynamic eddy viscosity m_ij = -2C_S|S1|.alpha^2.G1
	    les_tensor_gi(:,:,gi) = 2*alpha**2*les_coef_gi(gi)*strain_mod(gi)*f_tensor(:,:,gi)
	  end do
	  
    end select
  
  end subroutine dynamic_smagorinsky_sgs
    
  subroutine velocity_gradient (ele, q_mesh, u, x, g_nl, velocity_bc)
  
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: u, x
    type(mesh_type), intent(in) :: q_mesh
    type(vector_field), intent(in), optional :: velocity_bc
    real, dimension(mesh_dim(u), mesh_dim(u), ele_loc(u,ele)), intent(inout) :: g_nl
    
    type(element_type), pointer :: u_shape, q_shape
    real, dimension(ele_loc(q_mesh,ele), ele_ngi(q_mesh,ele), mesh_dim(U)) :: dq_t
    real, dimension(ele_ngi(u,ele)) :: detwei
    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim) :: du_t
    real, dimension(ele_loc(u,ele), ele_loc(u,ele)) :: M_inv
    real, dimension(U%dim, ele_loc(q_mesh,ele), ele_and_faces_loc(U,ele)) :: Grad_u_mat_q
    integer, dimension(:), pointer :: neigh
    integer :: ni, dim1, dim2, ele_2, face, start, finish, loc
    
    loc = ele_loc(u, ele)
    
    neigh=>ele_neigh(U, ele)
    u_shape=>ele_shape(U,ele)
    call transform_to_physical(X, ele, u_shape, dshape=du_t, detwei=detwei)
    
    q_shape=>ele_shape(q_mesh, ele)
    call transform_to_physical(X, ele, q_shape , dshape=dq_t)
    Grad_U_mat_q(:, :, :loc) = -dshape_shape(dq_t, u_shape, detwei)

    ! get inverse mass
    M_inv = shape_shape(u_shape, u_shape, detwei)
    call invert(M_inv)
    
    ! Compute gradient of non-linear velocity
    do dim1=1,mesh_dim(u)
      do dim2=1,mesh_dim(u)
    	! interior contribution
    	g_nl(dim1,dim2,:)=matmul(grad_U_mat_q(dim2,:,:loc), ele_val(u,dim1,ele))

    	! boundary contributions (have to be done seperately as we need to apply bc's at boundaries)
    	! local node map counter.
    	start=loc+1
    	do ni=1,size(neigh)
    	  ! get neighbour ele, corresponding faces, and complete local node map
    	  ele_2=neigh(ni)

    	  if (ele_2>0) then
    	    ! obtain corresponding faces, and complete local node map
    	    face=ele_face(U, ele_2, ele)
    	    finish=start+face_loc(U, face)-1  
    	    ! for interior faces we use the face values  
    	    g_nl(dim1,dim2,:)=g_nl(dim1,dim2,:)+matmul(grad_U_mat_q(dim2,:,start:finish), face_val(u,dim1,face))
    	  else
    	    ! obtain corresponding faces, and complete local node map
    	    face=ele_face(U, ele, ele_2)
    	    finish=start+face_loc(U, face)-1 
    	    ! for boundary faces the value we use depends upon if a weak bc is applied
    	    if (present(velocity_bc)) then
    	      ! weak bc! use the bc value
    	      g_nl(dim1,dim2,:)=g_nl(dim1,dim2,:)+matmul(grad_U_mat_q(dim2,:,start:finish), ele_val(velocity_bc,dim1,face))
    	    else
    	      ! no weak bc, use node values on internal face
    	      g_nl(dim1,dim2,:)=g_nl(dim1,dim2,:)+matmul(grad_U_mat_q(dim2,:,start:finish), face_val(u,dim1,face))
    	    end if
    	  end if

    	  ! update node map counter
    	  start=start+face_loc(U, face)
    	end do

    	! apply inverse mass
    	g_nl(dim1,dim2,:)=matmul(M_inv, g_nl(dim1,dim2,:))
      end do
    end do
        
  end subroutine velocity_gradient

end module les_module
