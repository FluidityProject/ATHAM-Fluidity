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

module helmholtz_projection

  use elements
  use sparse_tools
  use fetools
  use dgtools
  use fields
  use fefields
  use state_module
  use shape_functions
  use transform_elements
  use fldebug
  use field_options
  use field_derivatives
  use fields_calculations
  use spud
  use petsc_solve_state_module
  use sparsity_patterns
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN
  use state_fields_module
  
  public :: project_irrotational_velocity, project_irrotational_velocity_simple, &
            project_velocity_divergence
  
  contains
      
    subroutine project_irrotational_velocity(state, T, divergence, U_new, X, lvelocity_name)
   	 	      
       type(state_type), intent(inout) :: state
       type(vector_field), intent(inout) :: X, U_new
       type(scalar_field), intent(inout) :: T
       type(scalar_field), intent(inout) :: divergence
  
       integer, dimension(:), pointer :: test_nodes, neigh
       real, dimension(:,:,:), allocatable :: dtest_t, sdtest_t, ddensity_t, dt_t
       real, dimension(:), allocatable :: detwei, skew_div_at_quad, density_at_quad, olddensity_at_quad
       real, dimension(:,:), allocatable :: stiff_mat, rhs_mat
       real, dimension(:,:), allocatable :: skew_grad_at_quad, density_grad_at_quad, um_at_quad, u_at_quad, oldu_at_quad

       type(mesh_type), pointer :: test_mesh
       type(scalar_field) :: rhs_d, rhs_s, potential, stream, component
       type(vector_field) :: U_div, U_mean
       type(csr_matrix) :: stiff_m, skew_m
       
       type(vector_field), pointer :: v_field, U, U_old, U_proj
       type(scalar_field), pointer :: s_field, density, old_density
       type(element_type), pointer :: test_shape, d_shape, t_shape
       type(csr_sparsity), pointer :: m_sparsity
       
       character(len=FIELD_NAME_LEN) :: lvelocity_name
       integer :: ele, face, ni, ng, dim, stat
       real :: dt, umean
       
       ewrite(2,*) "Enter project_irrotational_velocity"
       
       call zero(U_new)
       call zero(divergence)
       
       density=>extract_scalar_field(state, "Density")
       old_density=>extract_scalar_field(state, "OldDensity")
       U_proj=>extract_vector_field(state, trim(lvelocity_name))
       test_mesh=>extract_mesh(state, trim(U_new%mesh%name))
       m_sparsity => get_csr_sparsity_firstorder(state, test_mesh, test_mesh)
       
       call get_option("/timestepping/timestep", dt)
       
       ! aux. field to store increment between subcycles
       call allocate(rhs_s, test_mesh, "RHS_s")
       call zero(rhs_s)
       call allocate(rhs_d, test_mesh, "RHS_d")
       call zero(rhs_d)
       call allocate(stream, test_mesh, "Stream")
       call zero(stream)
       call allocate(potential, test_mesh, "Potential")
       call zero(potential)
       call allocate(stiff_m, m_sparsity, name="Stiffness_m")
       call zero(stiff_m)
       call allocate(skew_m, m_sparsity, name="Skew_m")
       call zero(skew_m)
       call allocate(U_div, U_new%dim, U_new%mesh, "DivergenceVelocity")
       call zero(U_div)
       call allocate(U_mean, U_new%dim, U_new%mesh, "MeanVelocity")
       call zero(U_mean)
       
       allocate(detwei(ele_ngi(test_mesh, 1)), u_at_quad(U_proj%dim, ele_ngi(U_proj, 1)), um_at_quad(U_proj%dim, ele_ngi(U_proj, 1)), &
                skew_div_at_quad(ele_ngi(U_proj, 1)), skew_grad_at_quad(X%dim, ele_ngi(U_proj, 1)), & 
                density_at_quad(ele_ngi(density, 1)), olddensity_at_quad(ele_ngi(density, 1)), &
		density_grad_at_quad(X%dim, ele_ngi(density,1)), &
   	 	dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), x%dim), sdtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), x%dim), &
   	 	ddensity_t(ele_loc(density, 1), ele_ngi(density, 1), x%dim), &
   	 	rhs_mat(ele_loc(test_mesh, 1), ele_loc(test_mesh, 1)), &
		stiff_mat(ele_loc(test_mesh, 1), ele_loc(test_mesh, 1)))

       ! Compute mean velocity
       do dim=1,U_proj%dim
         component=extract_scalar_field_from_vector_field(U_proj, dim)
	 umean = mean(component)
	 call allmean(umean)
	 call set(U_mean, dim, umean)
       enddo
       ewrite_minmax(U_mean)
       
       ! Compute 2D stream function
       do ele=1, ele_count(test_mesh)

   	 test_nodes => ele_nodes(test_mesh, ele)
   	 test_shape => ele_shape(test_mesh, ele)
   	 call transform_to_physical(X, ele, test_shape, dshape=dtest_t, detwei=detwei)	 
	 sdtest_t(:,:,1)=-dtest_t(:,:,2)
	 sdtest_t(:,:,2)=dtest_t(:,:,1)
	 
	 ! Assemble skew stiffness matrix
	 call addto(skew_m, test_nodes, test_nodes, dshape_dot_dshape(sdtest_t, sdtest_t, detwei))

   	 ! Assemble RHS_s
	 skew_div_at_quad=0.0
         do dim=1,U_proj%dim
           skew_div_at_quad=skew_div_at_quad+matmul(ele_val(U_proj, dim, ele),sdtest_t(:,:,dim))
         end do
	 call addto(rhs_s, test_nodes, -shape_rhs(test_shape, detwei*skew_div_at_quad))

       enddo
       
       ! Boundary conditions: zero skew gradient
       do ni = 1, surface_element_count(rhs_s)
         call addto_diag(skew_m, face_global_nodes(rhs_s, ni), spread(INFINITY, 1, face_loc(rhs_s, ni)))
       end do

       ! Solve for stream function
       call petsc_solve(stream, skew_m, rhs_s, option_path=T%option_path)
       call halo_update(stream)

       s_field=>extract_scalar_field(state, "ScalarStream", stat=stat)
       if (stat==0) call set(s_field,stream)
       
       call compute_grad(state, stream, X, T, U_new, .true.)
       
       ! Correct divergence free velocity
       call zero(skew_m)
       call zero(rhs_s)
       do ele=1, ele_count(test_mesh)
   	 test_nodes => ele_nodes(test_mesh, ele)
   	 test_shape => ele_shape(test_mesh, ele)
   	 call transform_to_physical(X, ele, test_shape, dshape=dtest_t, detwei=detwei)	 
	 
	 ! Assemble skew stiffness matrix
	 call addto(skew_m, test_nodes, test_nodes, dshape_dot_dshape(dtest_t, dtest_t, detwei))
	 call addto(rhs_s, test_nodes, -shape_rhs(test_shape, detwei*ele_div_at_quad(U_new, ele, dtest_t)))
       enddo
       
       ! Solve for error potential
       call petsc_solve(potential, skew_m, rhs_s, option_path=T%option_path)
       call halo_update(potential)
       
       call compute_grad(state, potential, X, T, U_div, .false.)
       call addto(U_new, U_div, scale=-1.0)
       call addto(U_new, U_mean)
  
       ! Calculate potential
       call zero(potential)
       do ele=1, ele_count(test_mesh)

   	 test_nodes => ele_nodes(test_mesh, ele)
   	 test_shape => ele_shape(test_mesh, ele)
   	 call transform_to_physical(X, ele, test_shape, dshape=dtest_t, detwei=detwei)
	 sdtest_t(:,:,1)=-dtest_t(:,:,2)
	 sdtest_t(:,:,2)=dtest_t(:,:,1)

   	 d_shape => ele_shape(density, ele)
   	 call transform_to_physical(X, ele, d_shape, dshape=ddensity_t)

   	 u_at_quad = ele_val_at_quad(U_proj, ele)
   	 um_at_quad = ele_val_at_quad(U_mean, ele)
   	 density_at_quad = ele_val_at_quad(density, ele)
   	 olddensity_at_quad = ele_val_at_quad(old_density, ele)
         do dim=1,U_proj%dim
	   skew_grad_at_quad(dim,:) = matmul(ele_val(stream, ele),sdtest_t(:,:,dim))
	 enddo
	 
   	 density_grad_at_quad = theta_rho*(ele_grad_at_quad(density, ele, ddensity_t))+ &
   	 		     (1-theta_rho)*(ele_grad_at_quad(old_density, ele, ddensity_t))
  
   	 ! Assemble RHS_d
	 call addto(rhs_d, test_nodes, shape_rhs(test_shape, (1./dt)*detwei*(density_at_quad - olddensity_at_quad)))
         call addto(rhs_d, test_nodes, shape_rhs(test_shape, detwei*sum(skew_grad_at_quad*density_grad_at_quad,1)))
         call addto(rhs_d, test_nodes, shape_rhs(test_shape, detwei*sum(um_at_quad*density_grad_at_quad,1)))
	 
	 ! Assemble stiffness matrix
	 stiff_mat = dshape_dot_dshape(dtest_t, dtest_t, detwei*(theta_rho*density_at_quad+(1-theta_rho)*olddensity_at_quad))
	 call addto(stiff_m, test_nodes, test_nodes, stiff_mat)
	 	 
       enddo

	! Solve for scalar potential
       call petsc_solve(potential, stiff_m, rhs_d, option_path=T%option_path)
       call halo_update(potential)

       s_field=>extract_scalar_field(state,"ScalarPotential",stat=stat)
       if (stat==0) call set(s_field,potential)
              
       ! Reassemble projected velocity
       call zero(U_div)
       call compute_grad(state, potential, X, T, U_div, .false.)
       call addto(U_new, U_div)
       
       call compute_lap(state, potential, X, T, divergence)
       
       call deallocate(rhs_d)
       call deallocate(rhs_s)
       call deallocate(potential)
       call deallocate(stream)
       call deallocate(stiff_m)
       call deallocate(skew_m)
       call deallocate(U_div)
       call deallocate(U_mean)
       
    end subroutine project_irrotational_velocity
    
    subroutine project_irrotational_velocity_simple(state, T, divergence, X, lvelocity_name)
   	 	      
       type(state_type), intent(inout) :: state
       type(vector_field), intent(inout) :: X
       type(scalar_field), intent(inout) :: T
       type(scalar_field), intent(inout) :: divergence
  
       integer, dimension(:), pointer :: test_nodes, neigh
       real, dimension(:,:,:), allocatable :: dtest_t, ddensity_t, dt_t
       real, dimension(:), allocatable :: detwei, density_at_quad, olddensity_at_quad
       real, dimension(:,:), allocatable :: density_grad_at_quad, u_at_quad

       type(mesh_type), pointer :: test_mesh
       type(scalar_field) :: rhs_d
       type(csr_matrix) :: mass_m
       
       type(vector_field), pointer :: U_proj
       type(scalar_field), pointer :: density, old_density
       type(element_type), pointer :: test_shape, d_shape, t_shape
       type(csr_sparsity), pointer :: m_sparsity
       
       character(len=FIELD_NAME_LEN) :: lvelocity_name
       integer :: ele, face, ni, ng, dim, stat
       real :: dt, umean
       
       ewrite(2,*) "Enter project_irrotational_velocity_simple"
       
       call zero(divergence)
       
       density=>extract_scalar_field(state, "Density")
       old_density=>extract_scalar_field(state, "OldDensity")
       U_proj=>extract_vector_field(state, trim(lvelocity_name))
       test_mesh=>extract_mesh(state, trim(divergence%mesh%name))
       m_sparsity => get_csr_sparsity_firstorder(state, test_mesh, test_mesh)
       
       call get_option("/timestepping/timestep", dt)
       
       ! aux. field to store increment between subcycles
       call allocate(rhs_d, test_mesh, "RHS_d")
       call zero(rhs_d)
       call allocate(mass_m, m_sparsity, name="Mass_m")
       call zero(mass_m)
       
       allocate(detwei(ele_ngi(test_mesh, 1)), u_at_quad(U_proj%dim, ele_ngi(U_proj, 1)), &
                density_at_quad(ele_ngi(density, 1)), olddensity_at_quad(ele_ngi(density, 1)), &
		density_grad_at_quad(X%dim, ele_ngi(density,1)), &
   	 	dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), x%dim), &
   	 	ddensity_t(ele_loc(density, 1), ele_ngi(density, 1), x%dim))
  
       ! Calculate potential
       do ele=1, ele_count(test_mesh)

   	 test_nodes => ele_nodes(test_mesh, ele)
   	 test_shape => ele_shape(test_mesh, ele)
   	 call transform_to_physical(X, ele, test_shape, dshape=dtest_t, detwei=detwei)

   	 d_shape => ele_shape(density, ele)
   	 call transform_to_physical(X, ele, d_shape, dshape=ddensity_t)

   	 u_at_quad = ele_val_at_quad(U_proj, ele)
   	 density_at_quad = ele_val_at_quad(density, ele)
   	 olddensity_at_quad = ele_val_at_quad(old_density, ele)
	 
   	 density_grad_at_quad = theta_rho*(ele_grad_at_quad(density, ele, ddensity_t))+ &
   	 		     (1-theta_rho)*(ele_grad_at_quad(old_density, ele, ddensity_t))
  
   	 ! Assemble RHS_d
	 call addto(rhs_d, test_nodes, shape_rhs(test_shape, (1./dt)*detwei*(density_at_quad - olddensity_at_quad)))
         call addto(rhs_d, test_nodes, shape_rhs(test_shape, detwei*sum(u_at_quad*density_grad_at_quad,1)))
	 
	 ! Assemble stiffness matrix
	 call addto(mass_m, test_nodes, test_nodes, -shape_shape(test_shape, test_shape, detwei*(theta_rho*density_at_quad+(1-theta_rho)*olddensity_at_quad)))
	 	 
       enddo

	! Solve for scalar potential
       call petsc_solve(divergence, mass_m, rhs_d, option_path=T%option_path)
       
       call deallocate(rhs_d)
       call deallocate(mass_m)
       
    end subroutine project_irrotational_velocity_simple
    
    subroutine project_velocity_divergence(state, T, divergence, X, lvelocity_name)
   	 	      
       type(state_type), intent(inout) :: state
       type(vector_field), intent(inout) :: X
       type(scalar_field), intent(in) :: T
       type(scalar_field), intent(inout) :: divergence
       character(len=FIELD_NAME_LEN), intent(in) :: lvelocity_name
  
       integer, dimension(:), pointer :: test_nodes
       real, dimension(:,:,:), allocatable :: rhs_mat, dtest_t, du_t
       real, dimension(:), allocatable :: detwei, detwei_u

       type(mesh_type), pointer :: test_mesh, u_mesh
       type(scalar_field) :: rhs_d
       type(csr_matrix) :: mass_m
       
       type(vector_field), pointer :: U_proj
       type(element_type), pointer :: test_shape, u_shape
       type(csr_sparsity), pointer :: m_sparsity
       integer :: ele, dim
       
       ewrite(2,*) "Enter project_velocity_divergence"
       
       call zero(divergence)
       
       U_proj=>extract_vector_field(state, trim(lvelocity_name))
       test_mesh=>extract_mesh(state, trim(divergence%mesh%name))
       u_mesh=>extract_mesh(state, trim(U_proj%mesh%name))
       m_sparsity => get_csr_sparsity_firstorder(state, test_mesh, test_mesh)

       ! aux. field to store increment between subcycles
       call allocate(rhs_d, test_mesh, "RHS_d")
       call zero(rhs_d)
       call allocate(mass_m, m_sparsity, name="Mass_m")
       call zero(mass_m)
	
       allocate(detwei(ele_ngi(test_mesh, 1)), detwei_u(ele_ngi(U_proj, 1)), &
       		rhs_mat(x%dim, ele_loc(test_mesh, 1), ele_loc(u_mesh, 1)), &
   	 	du_t(ele_loc(u_mesh, 1), ele_ngi(u_mesh, 1), x%dim), &
   	 	dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), x%dim))
  
       ! Calculate potential
       do ele=1, ele_count(test_mesh)

   	 test_nodes => ele_nodes(test_mesh, ele)
   	 test_shape => ele_shape(test_mesh, ele)
   	 call transform_to_physical(X, ele, test_shape, dshape=dtest_t, detwei=detwei)

   	 u_shape => ele_shape(test_mesh, ele)
   	 call transform_to_physical(X, ele, u_shape, dshape=du_t, detwei=detwei_u)
  
   	 ! Assemble RHS_d
	 rhs_mat=shape_dshape(test_shape, du_t, detwei)
	 do dim=1,U_proj%dim
	   call addto(rhs_d, test_nodes, matmul(rhs_mat(dim,:,:),ele_val(U_proj, dim, ele)))
	 enddo
	 
	 ! Assemble stiffness matrix
	 call addto(mass_m, test_nodes, test_nodes, shape_shape(test_shape, test_shape, detwei))
	 	 
       enddo

	! Solve for scalar potential
       call petsc_solve(divergence, mass_m, rhs_d, option_path=T%option_path)
       
       call deallocate(rhs_d)
       call deallocate(mass_m)
       
    end subroutine project_velocity_divergence

    subroutine compute_div(state, vector, X, T, div)
      
      type(state_type)   :: state
      type(scalar_field) :: T, div
      type(scalar_field) :: component_grad, component_vect
      type(vector_field) :: X, vector, grad_loc
      integer :: dim
      
      call allocate(grad_loc, X%dim, div%mesh, 'ScalarGradient')
            
      call zero(div)
      do dim = 1, X%dim
        call zero(grad_loc)
        component_vect=extract_scalar_field_from_vector_field(vector, dim)
        call compute_grad(state, component_vect, X, T, grad_loc, .false.)
	
        component_grad=extract_scalar_field_from_vector_field(grad_loc, dim)
        call addto(div, component_grad)
      enddo
      
      call deallocate(grad_loc)
      
    end subroutine compute_div

    subroutine compute_grad(state, scalar, X, T, grad_loc, skew)
      
      type(state_type)   :: state
      type(scalar_field) :: T, scalar
      type(vector_field) :: X, grad_loc
      logical :: skew
      
      integer :: ele
      type(element_type), pointer :: test_shape
      integer, dimension(:), pointer :: test_nodes, neigh
      real, dimension(:,:,:), allocatable :: dtest_t, sdtest_t
      real, dimension(:), allocatable :: detwei

      type(mesh_type), pointer :: test_mesh
      type(csr_sparsity), pointer :: m_sparsity
      type(vector_field) :: rhs
      type(csr_matrix) :: mass
       
      test_mesh=>extract_mesh(state, trim(scalar%mesh%name))
      m_sparsity => get_csr_sparsity_firstorder(state, test_mesh, test_mesh)
      
      allocate(detwei(ele_ngi(test_mesh, 1)), &
   	       dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), x%dim), &
	       sdtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), x%dim))
      
      call allocate(rhs, X%dim, test_mesh, 'RHS_matrix')
      call allocate(mass, m_sparsity, name='Mass_matrix')
      
      call zero(mass)
      call zero(rhs)
      do ele=1, ele_count(scalar)

        test_nodes => ele_nodes(scalar, ele)
        test_shape => ele_shape(scalar, ele)
        call transform_to_physical(X, ele, test_shape, dshape=dtest_t, detwei=detwei)
	
	if (skew) then
          sdtest_t(:,:,1)=-dtest_t(:,:,2)
          sdtest_t(:,:,2)=dtest_t(:,:,1)
        else
	  sdtest_t=dtest_t
	endif
  
        ! Assemble RHS
        call addto(rhs, test_nodes, shape_vector_rhs(test_shape, ele_grad_at_quad(scalar, ele, dtest_t), detwei))
        call addto(mass, test_nodes, test_nodes, shape_shape(test_shape, test_shape, detwei))
        	
      enddo

       ! Solve for scalar potential
      call petsc_solve(grad_loc, mass, rhs, option_path=T%option_path)
      
      call deallocate(rhs)
      call deallocate(mass)
      
    end subroutine compute_grad

    subroutine compute_lap(state, scalar, X, T, laplace)
      
      type(state_type)   :: state
      type(scalar_field) :: T, scalar, laplace
      type(vector_field) :: X
      logical :: skew
      
      integer :: ele, j
      type(element_type), pointer :: test_shape
      integer, dimension(:), pointer :: test_nodes, neigh
      real, dimension(:,:,:), allocatable :: dtest_t
      real, dimension(:,:), allocatable :: stiff_mat, grad_gi
      real, dimension(:), allocatable :: detwei

      type(mesh_type), pointer :: test_mesh
      type(csr_sparsity), pointer :: m_sparsity
      type(csr_matrix) :: mass
      type(scalar_field) :: rhs
       
      test_mesh=>extract_mesh(state, trim(scalar%mesh%name))
      m_sparsity => get_csr_sparsity_firstorder(state, test_mesh, test_mesh)
      
      allocate(detwei(ele_ngi(test_mesh, 1)), grad_gi(X%dim, ele_ngi(test_mesh, 1)), &
               stiff_mat(ele_loc(test_mesh, 1), ele_loc(test_mesh, 1)), &
   	       dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), x%dim))
      
      call allocate(rhs, test_mesh, name='RHS_matrix')
      call allocate(mass, m_sparsity, name='Mass_matrix')
      
      call zero(mass)
      call zero(rhs)
      do ele=1, ele_count(scalar)

        test_nodes => ele_nodes(scalar, ele)
        test_shape => ele_shape(scalar, ele)
        call transform_to_physical(X, ele, test_shape, dshape=dtest_t, detwei=detwei)
  
        ! Assemble RHS
        call addto(rhs, test_nodes, -dshape_dot_vector_rhs(dtest_t, ele_grad_at_quad(scalar, ele, dtest_t), detwei))
        call addto(mass, test_nodes, test_nodes, shape_shape(test_shape, test_shape, detwei))
        	
      enddo

       ! Solve for scalar potential
      call petsc_solve(laplace, mass, rhs, option_path=T%option_path)
      
      call deallocate(rhs)
      call deallocate(mass)
      
    end subroutine compute_lap

end module helmholtz_projection
