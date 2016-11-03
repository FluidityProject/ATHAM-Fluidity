!    Copyright (C) 2009 Imperial College London and others.
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

module slope_limiters_dg
use fefields
use fldebug_parameters
use ieee_arithmetic
use spud
use elements
use eventcounter
use transform_elements
use fields
use state_module
use vtk_interfaces
use state_fields_module
use bound_field_module
use node_boundary
use global_parameters, only: FIELD_NAME_LEN
use diagnostic_fields, only : safe_set

implicit none

private
public limit_slope_dg, limit_fpn, limit_vb

integer, parameter :: LIMITER_MINIMAL=1
integer, parameter :: LIMITER_COCKBURN=2
integer, parameter :: LIMITER_HERMITE_WENO=3
integer, parameter :: LIMITER_FPN=4
integer, parameter :: LIMITER_VB=5
integer, parameter :: LIMITER_VBJ=6

public :: LIMITER_MINIMAL, LIMITER_COCKBURN, LIMITER_HERMITE_WENO,&
     & LIMITER_FPN, LIMITER_VB, LIMITER_VBJ

!!CockburnShuLimiter stuff
real :: TVB_factor=0.1
logical :: TVB_factor_absolute=.false., TVB_factor_relative=.false.
real :: Limit_factor=1.1, tolerance=1.e-5
real, dimension(:,:,:), pointer :: alpha => null()
real, dimension(:,:), pointer :: dx2 => null(), dxn2 => null(), A => null()
integer :: CSL_adapt_counter = -666
logical :: CSL_initialised = .false.
logical :: tolerate_negative_weights

!!Hermite Weno limiter stuff
real :: gam0 !power coefficient in weights
real :: eps_o !relative/absolute tolerance threshold for oscillation indicator
real :: eps_w !relative/absolute tolerance threshold for WENO weights
real :: disc_tol !Value for discontinuity test
real :: limit_tol !Do not limit if infinity norm of tracer is less than
!this value on an element
logical :: debugging !Switch to bung out lots of debugging output
integer, parameter :: IGNORE_MISSING_POLYS=1
integer, parameter :: REPLACE_MISSING_POLYS=2
integer, parameter :: LOWER_ORDER=3
integer :: missing_polys
logical :: leave_out_hermite_polynomials
logical :: has_discontinuity_detector_field
type(scalar_field), pointer :: discontinuity_detector_field
integer :: limit_count

contains

  subroutine limit_slope_dg(T, U, X, state, limiter, not_initialised)
    !! Assume 1D linear elements
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(inout) :: X, U
    type(state_type), intent(inout) :: state
    integer, intent(inout) :: limiter
    logical, intent(in), optional :: not_initialised

    integer :: ele, stat
    type(scalar_field) :: T_limit
    character(len=FIELD_NAME_LEN) :: limiter_name

    !assert(mesh_dim(coordinate)==1)
    !assert(field%mesh%continuity<0)
    !assert(field%mesh%shape%degree==1)

    ewrite(2,*) 'subroutiune limit_slope_dg'
    
    if (present(not_initialised)) then
       ! Note unsafe for mixed element meshes
       if (element_degree(T,1)==0) then
          FLExit("Slope limiters make no sense for degree 0 fields")
       end if

       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/slope_limiter/name",limiter_name)

       select case(trim(limiter_name))
       case("Cockburn_Shu")
          limiter=LIMITER_COCKBURN
       case("Hermite_Weno")
          limiter=LIMITER_HERMITE_WENO
       case("minimal")
          limiter=LIMITER_MINIMAL
       case("FPN")
          limiter=LIMITER_FPN
       case("Vertex_Based")
          limiter=LIMITER_VB
       case("Vertex_Based_Julien")
          limiter=LIMITER_VBJ
       case default
          FLAbort('No such limiter')
       end select
    endif
    
    select case (limiter)
    case (LIMITER_MINIMAL)
       T_limit=extract_scalar_field(state, trim(T%name)//"Limiter", stat=stat)
       
       do ele=1,element_count(T)
          
          if (stat==0) then
             call limit_slope_ele_dg(ele, T, X, T_limit)
          else
             call limit_slope_ele_dg(ele, T, X)
          end if

       end do

    case (LIMITER_COCKBURN)

       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Cockburn_Shu/TVB_factor_absolute", &
            &TVB_factor,stat=stat)
       TVB_factor_absolute=.true.
       if (stat/=0) then
         call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Cockburn_Shu/TVB_factor_relative", &
            &TVB_factor,stat=stat)
         TVB_factor_relative=.true.; TVB_factor_absolute=.false.
       endif
       
       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Cockburn_Shu/limit_factor", &
            &limit_factor)

       tolerate_negative_weights = &
            &have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Cockburn_Shu/tolerate_negative_weights")

       call cockburn_shu_setup(T, X)
              
       do ele=1,element_count(T)
          
          call limit_slope_ele_cockburn_shu(ele, T, X)
          
       end do

    case (LIMITER_HERMITE_WENO)

       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/power_coeffi&
            &cient", &
       & gam0)
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/tolerance_th&
            &reshold_oscillations", &
            &eps_o) 
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/tolerance_th&
            &reshold_weights", &
            &eps_w) 
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/discontinuit&
            &y_tolerance",disc_tol) 
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/limit_tolera&
            &nce",limit_tol) 
       debugging = have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/debugging")
       missing_polys = IGNORE_MISSING_POLYS
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &boundary_treatment::ignore_missing_polys"))&
            & missing_polys = IGNORE_MISSING_POLYS
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &boundary_treatment::replace_missing_polys"))&
            & missing_polys = REPLACE_MISSING_POLYS
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &boundary_treatment::lower_order")) then
          missing_polys = LOWER_ORDER
       end if

       leave_out_hermite_polynomials = .false.
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &leave_out_hermite_polynomials")) &
            & leave_out_hermite_polynomials = .true.

       call allocate(T_limit, T%mesh, name="NewT") 
       T_limit%val = T%val
       
       limit_count = 0.0

       has_discontinuity_detector_field = has_scalar_field( &
            state, "DiscontinuityDetector")
       if(has_discontinuity_detector_field) then
          discontinuity_detector_field &
               => extract_scalar_field(state, "DiscontinuityDetector")
          discontinuity_detector_field%val = 0.0
       end if

       do ele = 1, element_count(T)
          
          call limit_slope_ele_hermite_weno(ele, T, T_limit, X, U)

       end do

       ewrite(3,*) 'Limit count = ',limit_count

       T%val = T_limit%val
       call deallocate(T_limit)

    case (LIMITER_VB)
       call limit_VB(state, T)

    case (LIMITER_VBJ)
       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Vertex_Based_Julien/Limit_factor", &
            &limit_factor,stat=stat)
       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Vertex_Based_Julien/tolerance", &
            &tolerance,stat=stat)

       call limit_VB_julien(state, X, T)

    case (LIMITER_FPN)
       call limit_fpn(state, T)      

    case default
       ewrite(-1,*) 'limiter = ', limiter
       FLAbort('no such limiter exists')
    end select

    ewrite(2,*) 'END subroutiune limit_slope_dg'

  end subroutine limit_slope_dg

  subroutine limit_slope_ele_dg(ele, T, X, T_limit)
    
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(inout), optional :: T_limit
    integer, dimension(:), pointer :: neigh, neigh2, x_neigh, T_ele, face_nodes
    real, dimension(X%dim) :: ele_centre, face2_centre
    real :: ele_mean, miss_val, lim_value_l, lim_value_r
    integer :: ele_2, ni, ni2, face, face2, d, deg, i, j, jj, ii, miss
    real, dimension(X%dim, ele_loc(X,ele)) :: X_val, X_val_2
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    real, dimension(X%dim, ele_face_count(T,ele)) :: neigh_centre, face_centre
    real, dimension(ele_face_count(T,ele)) :: neigh_mean, face_mean, new_val_f
    integer, dimension(T%mesh%shape%numbering%vertices) :: ele_vertices      
    real, dimension(ele_loc(T,ele)) :: b, new_val
    logical :: limit
    real :: alpha_l=1.
    
    if(associated(A)) then
       deallocate(A)
       A => null()
    end if

    X_val=ele_val(X, ele)
    T_val=ele_val(T, ele)
    
    ele_vertices = local_vertices(T%mesh%shape)

    ele_centre=sum(X_val,2)/size(X_val,2)
    
    ele_mean=sum(T_val(ele_vertices))/size(ele_vertices)
    
    neigh=>ele_neigh(T, ele)
    face_mean=0.0
    x_neigh=>ele_neigh(X, ele)

    limit=.false.

    searchloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       
       ! These finding routines are outside the inner loop so as to allow
       ! for local stack variables of the right size in
       ! construct_add_diff_interface_dg.

       ele_2=neigh(ni)
    
       ! Note that although face is calculated on field U, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(T, ele, ele_2)
       face_nodes=>face_local_nodes(T%mesh, face)

       face_centre(:,ni) = sum(face_val(X,face),2)/size(face_val(X,face),2)
       
       face_mean(ni)=0.
       do j = 1,size(face_nodes)
         if (any(ele_vertices==face_nodes(j))) then    
           face_mean(ni)=face_mean(ni)+T_val(face_nodes(j))/real(size(ele_vertices)-1)
         endif
       enddo
       
       if (ele_2<=0) then
          ! External face.
          cycle
       end if

       X_val_2=ele_val(X, ele_2)
       T_val_2=ele_val(T, ele_2)

       neigh_centre(:,ni)=sum(X_val_2,2)/size(X_val_2,2)
       
       neigh_mean(ni)=sum(T_val_2(ele_vertices))/size(ele_vertices)
    
       if ((face_mean(ni)-ele_mean)*(face_mean(ni)-neigh_mean(ni))>0.0) then
          ! Limit if face_mean does not lie between ele_mean and neigh_mean
          limit=.true.
          
	  lim_value_l=(1.-alpha_l)*ele_mean + alpha_l*min(ele_mean,neigh_mean(ni))
	  lim_value_r=(1.-alpha_l)*ele_mean + alpha_l*max(ele_mean,neigh_mean(ni))
	  face_mean(ni)=max(min(face_mean(ni),lim_value_r),lim_value_l)

       end if


    end do searchloop

    if (present(T_limit)) then
       T_ele=>ele_nodes(T_limit,ele)
       call set(T_limit, T_ele, ele_mean+0.0*T_ele)
    end if

    if (.not.limit) then
       return
    end if
    
    d=mesh_dim(T)
    new_val_f=ele_mean

    do miss=1,d+1
       
       ! If the missed side is a boundary, it is not possible to limit in
       ! this direction without violating the boundary condition.
       if (neigh(miss)<=0) cycle

       A=0.0
       b(1)=ele_mean

       do i=1, d+1
          ! Enforce preservation of the element mean value.
          A(1,i)=1.0/(d+1) 
          
          jj=1
          do j=1,d+1
             if (j==miss) cycle
             jj=jj+1
             
             if (i/=j) then
                A(jj,i)=1.0/d
             else
                b(jj)=face_mean(j)
             end if
             
          end do
          
       end do

       call invert(A)
       b=matmul(A,b)
       
       if (maxval(abs(b-ele_mean))>maxval(abs(new_val_f-ele_mean))) then
          !! The slope is larger than the current best guess.
          
          miss_val=0.0
          do ni=1, d+1
             if (ni==miss) cycle

             miss_val=miss_val+b(ni)/d
          end do
             
          if ((miss_val-ele_mean)*(miss_val-neigh_mean(miss))<=0.0) then
             ! The slope is legal.
             
             new_val_f=b

          end if

       end if

    end do
    
    new_val=T_val
    do ni=1,size(neigh)
      ele_2=neigh(ni)
      face=ele_face(T, ele, ele_2)
      face_nodes=>face_local_nodes(T%mesh, face)
      
      do j = 1,size(face_nodes)
        if (.not.any(ele_vertices==face_nodes(j))) &    
              new_val(face_nodes(j))=face_mean(ni)
      enddo
    enddo
    
    jj=0
    do j=1,ele_loc(T,ele)
      if (any(ele_vertices==j)) then
        jj=jj+1
        new_val(j)=new_val_f(jj)
      endif
    enddo
   
    ! Success or non-boundary failure.
    T_ele=>ele_nodes(T,ele)
    
    call set(T, T_ele, new_val)

    if (present(T_limit)) then
       T_ele=>ele_nodes(T_limit, ele)
       
       if (all(new_val==ele_mean)) then
          call set(T_limit, T_ele, 1.0+T_ele*0.0)
       else
          call set(T_limit, T_ele, -1.0+0.0*T_ele)
       end if

    end if

  end subroutine limit_slope_ele_dg

  subroutine cockburn_shu_setup(T,X)
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    !
    logical :: do_setup
    integer :: cnt, d, i, j, ele, deg, ele_2, face
    integer, dimension(:), pointer :: face_nodes, neigh
    integer, dimension(T%mesh%shape%numbering%vertices) :: ele_vertices  
    
    do_setup = .false.
    if(.not.CSL_initialised) then
       CALL GetEventCounter(EVENT_ADAPTIVITY, csl_adapt_counter)
       do_setup = .true.
       CSL_initialised = .true.
    else 
       CALL GetEventCounter(EVENT_ADAPTIVITY, CNT)
       if(cnt.ne.csl_adapt_counter) then
          do_setup= .true.
          csl_adapt_counter = cnt
       end if
    end if

    if(do_setup) then

       if(associated(alpha)) then
          deallocate(alpha)
          alpha => null()
       end if
       if(associated(dx2)) then
          deallocate(dx2)
          dx2 => null()
       end if
       if(associated(A)) then
          deallocate(A)
          A => null()
       end if

       !!ATTENTION: This assumes that all elements have the same number of faces
       allocate(alpha(element_count(T),ele_face_count(T,1)&
            &,ele_face_count(T,1)))
       allocate(dx2(element_count(T),ele_face_count(T,1)))

       d=mesh_dim(T)
       allocate(A(d+1,d+1))

       ! Initialise A with the change from face centre values to node values.
       do i=1, size(A,1)
          do j=1,size(A,2)
             if (i==j) then
                A(i,j)=0.0
             else
                A(i,j)=1.0/d
             end if
          end do
       end do
       
       call invert(A)
       
       do ele = 1, element_count(T)

          call cockburn_shu_setup_ele(ele,T,X)
	  
       end do

    end if    

  end subroutine cockburn_shu_setup

  subroutine cockburn_shu_setup_ele(ele, T, X)
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    
    integer, dimension(:), pointer :: neigh, x_neigh, face_nodes
    real, dimension(X%dim) :: ele_centre, face_2_centre
    real :: max_alpha, min_alpha, neg_alpha
    integer :: ele_2, ni, nj, face, face_2, i, nk, ni_skip, info, nl, j
    real, dimension(X%dim, ele_loc(X,ele)) :: X_val, X_val_2
    real, dimension(X%dim, ele_face_count(T,ele)) :: neigh_centre, face_centre
    real, dimension(X%dim) :: alpha1, alpha2
    real, dimension(X%dim,X%dim) :: alphamat
    real, dimension(X%dim,X%dim+1) :: dx_f, dx_c
    integer, dimension(T%mesh%shape%numbering%vertices) :: ele_vertices  
    
    ele_vertices = local_vertices(T%mesh%shape)
    
    X_val=ele_val(X, ele)
    
    ele_centre=sum(X_val,2)/size(X_val,2)
    
    neigh=>ele_neigh(T, ele)
    ! x_neigh/=t_neigh only on periodic boundaries.
    x_neigh=>ele_neigh(X, ele)

    searchloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       ele_2=neigh(ni)
    
       ! Note that although face is calculated on field U, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(T, ele, ele_2)
       face_nodes=>face_local_nodes(X%mesh, face)

       face_centre(:,ni) = sum(X_val(:,face_nodes),2)/size(face_nodes)

       if (ele_2<=0) then
          ! External face.
          neigh_centre(:,ni)=face_centre(:,ni)
          cycle
       end if

       X_val_2=ele_val(X, ele_2)
       
       neigh_centre(:,ni)=sum(X_val_2,2)/size(X_val_2,2)
       if (ele_2/=x_neigh(ni)) then
          ! Periodic boundary case. We have to cook up the coordinate by
          ! adding vectors to the face from each side.
          face_2=ele_face(T, ele_2, ele)
          face_nodes=>face_local_nodes(X%mesh, face_2)
	  
          face_2_centre = sum(X_val_2(:,face_nodes),2)/size(face_nodes)
          neigh_centre(:,ni)=face_centre(:,ni) + &
               (neigh_centre(:,ni) - face_2_centre)
       end if

    end do searchloop

    do ni = 1, size(neigh)
       dx_c(:,ni)=neigh_centre(:,ni)-ele_centre !Vectors from ni centres to
                                                !ele centre
       dx_f(:,ni)=face_centre(:,ni)-ele_centre !Vectors from ni face centres
                                               !to ele centre
    end do
    
    alpha_construction_loop: do ni = 1, size(neigh)
       !Loop for constructing Delta v(m_i,K_0) as described in C&S
       alphamat(:,1) = dx_c(:,ni)

       max_alpha = -1.0
       ni_skip = 0

       choosing_best_other_face_loop: do nj = 1, size(neigh)
          !Loop over the other faces to choose best one to use
          !for linear basis across face

          if(nj==ni) cycle
          
          !Construct a linear basis using all faces except for nj
          nl = 1
          do nk = 1, size(neigh)
             if(nk==nj.or.nk==ni) cycle
             nl = nl + 1
             alphamat(:,nl) = dx_c(:,nk)
          end do
          
          !Solve for basis coefficients alpha
          alpha2 = dx_f(:,ni)
          call solve(alphamat,alpha2,info)

          if((.not.any(alpha2<0.0)).and.alpha2(1)/norm2(alpha2)>max_alpha) &
               & then
             alpha1 = alpha2
             ni_skip = nj
             max_alpha = alpha2(1)/norm2(alpha2)
          end if

       end do choosing_best_other_face_loop

       if(max_alpha<0.0) then
          if(tolerate_negative_weights) then
             min_alpha = huge(0.0)
             ni_skip = 0
             choosing_best_other_face_neg_weights_loop: do nj = 1, size(neigh)
                !Loop over the other faces to choose best one to use
                !for linear basis across face
                
                if(nj==ni) cycle
                
                !Construct a linear basis using all faces except for nj
                nl = 1
                do nk = 1, size(neigh)
                   if(nk==nj.or.nk==ni) cycle
                   nl = nl + 1
                   alphamat(:,nl) = dx_c(:,nk)
                end do
                
                !Solve for basis coefficients alpha
                alpha2 = dx_f(:,ni)
                call solve(alphamat,alpha2,info)

                neg_alpha = 0.0
                do i = 1, size(alpha2)
                   if(alpha2(i)<0.0) then
                      neg_alpha = neg_alpha + alpha2(i)**2
                   end if
                end do
                neg_alpha = sqrt(neg_alpha)

                if(min_alpha>neg_alpha) then
                   alpha1 = alpha2
                   ni_skip = nj
                   min_alpha = neg_alpha
                end if
             end do choosing_best_other_face_neg_weights_loop
          else
             FLAbort('solving for alpha failed')
          end if
       end if
       
       alpha(ele,ni,:) = 0.0
       alpha(ele,ni,ni) = alpha1(1)
       nl = 1
       do nj = 1, size(neigh)
          if(nj==ni.or.nj==ni_skip) cycle
          nl = nl + 1
          alpha(ele,ni,nj) = alpha1(nl)
       end do
       
       dx2(ele,ni) = norm2(dx_c(:,ni))

    end do alpha_construction_loop

  end subroutine cockburn_shu_setup_ele

  subroutine limit_slope_ele_cockburn_shu(ele, T, X)
    !!< Slope limiter according to Cockburn and Shu (2001) 
    !!< http://dx.doi.org/10.1023/A:1012873910884
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    
    integer, dimension(:), pointer :: neigh, x_neigh, T_ele, face_nodes
    real :: ele_mean
    real :: pos, neg, SumUp
    logical :: limit
    integer :: ele_2, ni, ni2, face, d, deg, i, j, jj
    real, dimension(X%dim,ele_loc(X,ele)) :: X_val
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    real, dimension(ele_loc(T,ele)) :: new_val
    real, dimension(ele_face_count(T,ele)) :: neigh_mean, face_mean
    real, dimension(mesh_dim(T)+1) :: delta_v
    real, dimension(mesh_dim(T)+1) :: Delta, new_val_f
    integer, dimension(T%mesh%shape%numbering%vertices) :: ele_vertices  
      
    T_val=ele_val(T, ele)
    X_val=ele_val(X, ele)
    
    ele_vertices = local_vertices(T%mesh%shape)
    
    ele_mean=sum(T_val(ele_vertices))/size(ele_vertices)
    
    neigh=>ele_neigh(T, ele)
    face_mean=0.0
    ! x_neigh/=t_neigh only on periodic boundaries.
    x_neigh=>ele_neigh(X, ele)
    
    limit=.false.

    searchloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       ele_2=neigh(ni)
    
       ! Note that although face is calculated on field U, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(T, ele, ele_2)
       face_nodes=>face_local_nodes(T%mesh, face)
    
       face_mean(ni)=0.
       do j = 1,size(face_nodes)
         if (any(ele_vertices==face_nodes(j))) then    
           face_mean(ni)=face_mean(ni)+T_val(face_nodes(j))/mesh_dim(T)
         endif
       enddo
       
       if (ele_2<=0) then
          ! External face.
          neigh_mean(ni)=face_mean(ni)
          cycle
       end if

       T_val_2=ele_val(T, ele_2)

       neigh_mean=sum(T_val_2(ele_vertices))/size(ele_vertices)
    
    end do searchloop

    delta_v = matmul(alpha(ele,:,:),neigh_mean-ele_mean)

    delta_loop: do ni=1,size(neigh)

       Delta(ni)=TVB_minmod(face_mean(ni)-ele_mean, Limit_factor*delta_v(ni), dx2(ele,ni))

    end do delta_loop
        
    ! Apply limiting in the element only if actually needed
    SumUp=0.0
    do ni = 1, size(neigh)
      SumUp = SumUp+Delta(ni)+ele_mean
    enddo
    limit=(sum(face_mean)/=SumUp)
    
    if (.not.limit) then
       return
    end if

    if (abs(sum(Delta))>1000.0*epsilon(0.0)) then
       ! Coefficients do not sum to 0.0

       pos=sum(max(0.0, Delta))
       neg=sum(max(0.0, -Delta))
       
       Delta = min(1.0,neg/pos)*max(0.0,Delta) &
             - min(1.0,pos/neg)*max(0.0,-Delta)
       
    end if

    new_val_f=matmul(A,Delta+ele_mean)
    
    new_val=T_val
    do ni=1,size(ele_vertices)
      new_val(ele_vertices(ni))=new_val_f(ni)
    enddo
      
    ! Success or non-boundary failure.
    T_ele=>ele_nodes(T,ele)
    
    call set(T, T_ele, new_val)

  end subroutine limit_slope_ele_cockburn_shu

  function TVB_minmod(a1,a2, dx)
    real :: TVB_minmod
    real, intent(in) :: a1, a2, dx

    if (TVB_factor_absolute .and. abs(a1)<TVB_factor) then
       TVB_minmod=a1
    else if (TVB_factor_relative .and. abs(a1)<TVB_factor*dx**2) then
       TVB_minmod=a1
    else if (abs(a1)<abs(a2)) then
       TVB_minmod=a1
    else
       TVB_minmod=a2
    end if

  end function TVB_minmod

  !11:25 <Guest54276>     do ele_A=1,ele_count(old_position)
  !11:25 <Guest54276>       call local_coords_matrix(old_position, ele_A, 
  !                   inversion_matrices_A(:, :, ele_A))
  !11:25 <Guest54276>     end do
  
  !subroutine local_coords_matrix(positions, ele, mat)
  !inputs global coordinates
  !outputs local coordinates

  subroutine limit_slope_ele_hermite_weno(ele, T, T_limit, X, U)
    !!< Hermite Weno Slope limiter
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: T
    type(scalar_field), intent(inout) :: T_limit
    type(vector_field), intent(in) :: X, U

    integer, dimension(:), pointer :: neigh, x_neigh, T_ele, face_nodes
    real :: ele_mean, ele_mean_2
    real, dimension(ele_face_count(T,ele)) :: ele_means
    real :: residual
    integer :: ele_2, ni, nj, face, face_2,i, nk, info, nl
    integer :: l_face, l_face_2
    logical :: limit_slope
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    real, dimension(face_loc(T,1)) :: T_val_face
    real, dimension(face_ngi(T,1)) :: T_face_quad
    real, dimension(ele_face_count(T,ele),ele_loc(X,ele)) :: T_vals
    real, dimension(ele_face_count(T,ele),X%dim, ele_loc(X,ele)) :: X_vals
    real, dimension(X%dim, ele_loc(X,ele)) :: X_val
    real, dimension(ele_face_count(T,ele)) :: neigh_mean, face_mean
    real, dimension(ele_loc(T,ele)) :: new_val
    
    type(element_type), pointer :: shape_T
    real, dimension(ele_loc(T, ele), ele_ngi(T, ele), mesh_dim(T)) :: du_t
    real, dimension(ele_ngi(T,ele)) :: detwei
    real, dimension(ele_ngi(T,ele)) :: p_quad, T_quad
    real, dimension(1+2*ele_loc(T,ele),ele_loc(T,ele)) :: Polys
    real, dimension(ele_loc(T,ele)*2+1) :: Polys_o, Polys_w
    logical, dimension(ele_face_count(T,ele)) :: boundaries
    logical, dimension(ele_face_count(T,ele)) :: construct_Lagrange
    real, dimension(ele_loc(T,ele),ele_loc(T,ele)) :: Imat
    real, dimension(ele_loc(T,ele)) :: Irhs
    real, dimension(ele_loc(X,ele)) :: local_coords
    integer, dimension(face_loc(T,1)) :: l_face_list,l_face_list_2
    real, dimension(mesh_dim(T),ele_ngi(T,ele)) :: dp_quad
    integer, dimension(T%mesh%shape%numbering%vertices) :: ele_vertices  

    real :: Discontinuity_indicator, inflow_integral, h
    real, dimension(ele_loc(T,ele)) :: ones
    integer :: discontinuity_option 
    real :: face_max, face_min

    if(debugging) then
       ewrite(2,*) 'Limit_slope_Hermite_weno_ele'
    end if

    limit_slope = .false.

    boundaries = .false.
    construct_Lagrange = .true.

    T_val=ele_val(T, ele)
    X_val=ele_val(X, ele)

    ele_mean=sum(T_val)/size(T_val)

    neigh=>ele_neigh(T, ele)
    ! x_neigh/=t_neigh only on periodic boundaries.
    x_neigh=>ele_neigh(X, ele)
    
    ele_vertices = local_vertices(T%mesh%shape)

    discontinuity_option = 2

    select case(discontinuity_option)

    case (1)
       !=========================================================
       !Discontinuity detector using TVB condition
       !Checks solution on each face is between mean values of 
       !ele and ele_2
       !=========================================================
       do ni=1,size(neigh)

          !--------------------------------------------------------------------
          ! Find the relevant faces.
          !--------------------------------------------------------------------
          ele_2=neigh(ni)

          if(ele_2<0) cycle
          
          T_val_2 = ele_val(T,ele_2)
          face=ele_face(T, ele, ele_2)
          T_val_face = face_val(T, face)
          T_face_quad = face_val_at_quad(T,face)
          !face_max = maxval(T_val_face)
          !face_min = minval(T_val_face)
          face_max = maxval(T_face_quad)
          face_min = minval(T_face_quad)
          ele_mean_2 = sum(T_val_2)/size(T_val_2)

          if(face_max>max(ele_mean,ele_mean_2)+disc_tol) limit_slope = .true.
          if(face_min<min(ele_mean,ele_mean_2)-disc_tol) limit_slope = .true.

          if(has_discontinuity_detector_field) then
             if(limit_slope) then
                ewrite(3,*) 'cjc limit_slope', ele
                ones = 1.0
                T_ele=>ele_nodes(Discontinuity_Detector_field,ele)
                call set(Discontinuity_detector_field,T_ele,ones)
             end if
          end if
       end do

    case (2)

       !=================================================================
       !DISCONTINUITY INDICATOR,
       !from http://www.gce.ucl.ac.be/~remacle/pdf/detect.pdf
       !We compute the jump of the solution on upwind boundaries
       !=================================================================

       !Initial value of integral of jump of solution on inflow boundaries
       Discontinuity_indicator = 0.0
       !Initial value of inflow area/length
       Inflow_integral = 0.0
       !We are going to increment these

       do ni=1,size(neigh)

          !--------------------------------------------------------------------
          ! Find the relevant faces.
          !--------------------------------------------------------------------
          ele_2=neigh(ni)

          if(ele_2<0) cycle

          face=ele_face(T, ele, ele_2)
          face_2=ele_face(T, ele_2, ele)

          call Discontinuity_indicator_face(Discontinuity_indicator, &
               & Inflow_integral, &
               & U,T,X,ele,face,face_2)

       end do

       discontinuity_indicator = abs(discontinuity_indicator)
       inflow_integral = abs(inflow_integral)

       !Compute h
       h = get_H(ele_vertices, X_val)

       !Get max norm in element of T
       T_quad = ele_val_at_quad(T,ele)

       if(Discontinuity_Indicator>disc_tol*Inflow_integral&
            &*maxval(abs(T_quad))*h) limit_slope = .true.

       if(has_discontinuity_detector_field) then
          ones = 1.0
          T_ele=>ele_nodes(Discontinuity_Detector_field,ele)
          call set(Discontinuity_detector_field,T_ele&
               &,Discontinuity_Indicator*ones/inflow_integral/&
               maxval(abs(T_quad)+limit_tol)/h)
       end if

    case default
       FLExit('no such discontinuity option')
    end select

    if(limit_slope) then
       limit_count = limit_count + 1

       !Apply HWENO limiter

       setuploop: do ni=1,size(neigh)

          ele_2=neigh(ni)

          if (ele_2<=0) then
             ! External face.
             neigh_mean(ni)=face_mean(ni)
             boundaries(ni) = .true.

             do nj = 1, size(neigh)
                if(ni==nj) cycle
                construct_Lagrange(nj) = .false.
             end do
             cycle

          end if

          ! Note that although face is calculated on field U, it is in fact
          ! applicable to any field which shares the same mesh topology.
          face=ele_face(T, ele, ele_2)
          face_2=ele_face(T, ele_2, ele)
          face_nodes=>face_local_nodes(T%mesh, face)

          face_mean(ni) = sum(T_val(face_nodes))/size(face_nodes)

          T_val_2=ele_val(T, ele_2)

          T_vals(ni,:) = T_val_2
          X_vals(ni,:,:) = ele_val(X, ele_2)

          neigh_mean(ni)=sum(T_val_2)/size(T_val_2)
          ele_means(ni) = neigh_mean(ni)

       end do setuploop

       if(any(boundaries).and.(missing_polys==LOWER_ORDER)) then
          !On boundary, with this option, just project to p(n-1)
          !We have only coded P1 so this projects to P0

          new_val = sum(T_val)/size(T_val)

       else

          Polys = 0.0

          if(debugging) then
             ewrite(2,*) 'Limiting slope.'
          end if

          !We store transformations (du_t, detwei) in the following way:
          ! 1:size(neigh) : transformations for neighbouring element
          ! size(neigh)+1 : transformations for this element

          shape_T=>ele_shape(T,ele)

          !Construct transformations for this element
          call transform_to_physical(X, ele,&
               & shape_T , dshape=du_t, detwei=detwei)

          !Polynomials are stored in the following way:
          ! i = 1:size(neigh) : Lagrange polynomials obtained by
          !                     missing out the i-th neighbour
          ! i = size(neigh)+1 : The existing polynomial representation
          ! j = size(neigh)+1 + i, i = 1:size(neigh) : The function with same mean
          !                                            as existing polynomial
          !                                            with slope taken from 
          !                                            i-th neighbour
          ! The latter representations are Hermite polynomials using 
          ! gradient information

          !Construct Lagrange polys
          !Fails if non-flat elements
          LagrangeP_loop: do ni = 1, size(neigh)
             if(.not.construct_Lagrange(ni)) then
                Polys(ni,:) = T_val
             else
                nl = 0
                do nj = 1, size(neigh)
                   if(nj==ni) cycle
                   nl = nl + 1
                   !
                   !This row requires that the mean of the polynomial over
                   !neighbour element nj be equal to the mean of the unlimited
                   !solution in that element.
                   ! 
                   ! This is done by computing the local coordinates of the
                   ! centre of the neighbour element (only works for P1) which
                   ! gives you the coefficients of the local expansion which
                   ! give you the polynomial value at that point which is then
                   ! required to be equal to the current element mean.
                   !
                   do nk = 1, size(neigh)
                      Imat(nl,:) = local_coords_interpolation(X,ele&
                           &,sum(X_vals(nj,:,:),2)/size(X_vals,3))
                   end do
                   Irhs(nl) = ele_means(nj)
                end do
                !Last column sets the mean value
                Imat(size(neigh),:) = 1.0/size(Imat,2)
                Irhs(size(neigh)) = ele_mean

                !Solve for the Polynomial
                call solve(Imat,Irhs,info)
                Polys(ni,:) = Irhs

                if(debugging) then
                   !Do some checking
                   !Compute the polynomial at the quad points
                   !Check polynomial has mean value ele_mean in this element
                   if(abs(sum(Polys(ni,:)-ele_mean))>1.0e-5) then
                      FLAbort('failed to get the correct mean value in this element')
                   end if
                   !Check polynomial has mean value ele_means in other two elements
                   do nj = 1, size(neigh)
                      if(nj==ni) cycle
                      local_coords = local_coords_interpolation(X,ele&
                           &,sum(X_vals(nj,:,:),2)/size(X_vals,3))
                      residual = sum(Polys(ni,:)*local_coords) - ele_means(nj)
                      if(abs(residual)>1.0e-5) then
                         FLAbort('failed to get the correct mean value in neighbour')
                      end if
                   end do
                end if

             end if
          end do LagrangeP_loop

          !Construct Hermite polys
          !Fails if non-flat elements

          !Original polynomial
          Polys(size(neigh)+1,:) = T_val
          HermiteP_loop: do ni = 1, size(neigh)
             !Mean of original values with slope of neighbouring element
             nk = size(neigh)+1+ni

             ele_2=neigh(ni)

             if(ele_2<0) then
                Polys(nk,:) = T_val
                cycle
             end if

             !face number in ele
             face=ele_face(T, ele, ele_2)
             !face number in ele_2
             face_2=ele_face(T, ele_2, ele)

             !local face number in ele
             l_face = local_face_number(T, face)
             !local face number in ele_2
             l_face_2 = local_face_number(T, face_2)

             !Local face list in ele
             l_face_list = face_local_nodes(T, face)
             !Local face list in ele_2
             l_face_list_2 = face_local_nodes(T, face_2)

             !T values in ele_2
             T_val_2=ele_val(T, ele_2)

             !First we "continue" the polynomial in ele_2 into ele

             !Polynomial takes same values on shared nodes
             Polys(nk,l_face_list) = &
                  & T_val_2(l_face_list_2)

             !Compute local coordinates (relative to ele)
             !of vertex in ele_2 which is opposite the face
             local_coords = local_coords_interpolation(X,ele,X_vals(ni,:,l_face_2))

             !Solve 1D linear system to get polynomial value
             !from
             !T_val_2(l_face_2) = sum(Poly*local_coords)
             !we have already computed the values on the face
             !so we rearrange.
             Polys(nk,l_face) = (T_val_2(l_face_2) - &
                  & sum(Polys(nk,l_face_list)*local_coords(l_face_list)))/ &
                  local_coords(l_face)

             !ADD SOME DEBUGGING TESTS

             !Second we adjust the mean so that it is the same as the 
             !mean of the current solution in ele
             Polys(nk,:) = Polys(nk,:)- &
                  & sum(Polys(nk,:))/size(T_val) + &
                  & sum(T_val)/size(T_val)         

          end do HermiteP_loop

          if(debugging) then
             ewrite(2,*) 'Dumping polynomials'
             do ni = 1, size(neigh)*2 + 1
                ewrite(2,*) Polys(ni,:)
             end do
          end if

          !Compute oscillatory indicators

          do ni = 1, size(neigh)*2 + 1
             !construct the ni-th polynomial at the quadrature points
             P_quad=matmul(Polys(ni,:), shape_T%n)
             !construct the gradient of the ni-th polynomial at the quad points
             do i = 1, mesh_dim(T)
                dP_quad(i,:) = &
                     &matmul(Polys(ni,:), du_t(:,:,i))
             end do

             !construct the oscillator index of the ni-th polynomial
             Polys_o(ni) = 0.0
             do i = 1, mesh_dim(T)
                Polys_o(ni) = Polys_o(ni) + &
                     & sum(detwei*dP_quad(i,:)**2)
             end do
             Polys_o(ni) = Polys_o(ni)/&
                  &sum(detwei*(eps_o + P_quad)**2)
          end do

          if(debugging) then
             ewrite(2,*) 'Dumping oscillatory indicators'
             ewrite(2,*) Polys_o
          end if

          !Compute weights
          do ni = 1, size(neigh)*2 + 1
             Polys_w(ni) = (eps_w + Polys_o(ni))**(-gam0)
          end do

          if(missing_polys==IGNORE_MISSING_POLYS.and.any(boundaries)) then
             do ni = 1, size(neigh)
                if(boundaries(ni)) then
                   do nj = 1, size(neigh)
                      if(ni==nj) cycle
                      Polys_w(nj) = 0.0
                   end do
                   Polys_w(size(neigh)+1+ni) = 0.0
                end if
             end do
          end if

          if(leave_out_hermite_polynomials) then
             Polys_w(size(neigh)+1:2*size(neigh)+1) = 0.0
          end if

          Polys_w = Polys_w/sum(Polys_w)

          if(debugging) then
             ewrite(2,*) 'Dumping weights'
             ewrite(2,*) Polys_w
          end if

          new_val = 0.
          do ni = 1, size(neigh)*2 + 1
             new_val = new_val + Polys_w(ni)*Polys(ni,:)
          end do

          if(debugging) then
             ewrite(2,*) 'new val is'
             ewrite(2,*) new_val

             ewrite(2,*) 'old slope was'
             do i = 1, mesh_dim(T)
                ewrite(2,*) maxval(matmul(T_val, du_t(:,:,i)))
             end do

             ewrite(2,*) 'new slope is'
             do i = 1, mesh_dim(T)
                ewrite(2,*) maxval(matmul(new_val, du_t(:,:,i)))
             end do
          end if
       end if

       T_ele=>ele_nodes(T,ele)
       call set(T_limit, T_ele, new_val)

    end if

  end subroutine limit_slope_ele_hermite_weno

  subroutine Discontinuity_indicator_face(Discontinuity_indicator, &
       & Inflow_integral, &
       & U,T,X,ele,face,face_2)
    real, intent(inout) :: Discontinuity_indicator, Inflow_integral
    type(vector_field), intent(in) :: U,X
    type(scalar_field), intent(in) :: T
    integer, intent(in) :: ele, face, face_2
    !
    real, dimension(face_ngi(T,face)) :: detwei
    real, dimension(mesh_dim(T),face_ngi(T,face)) :: normal
    real, dimension(mesh_dim(U),face_ngi(U,face)) :: U_flux
    integer, dimension(face_ngi(U,face)) :: inflow
    real :: Area

    call transform_facet_to_physical( X, face,&
    	 & detwei_f=detwei,normal=normal)
    
    U_flux = 0.5*(face_val_at_quad(U,face)+ &
    	 & face_val_at_quad(U,face_2))
    
    !We only compute on inflow boundaries    
    inflow = merge(1.0,0.0,sum(U_flux*normal,1)<0.0)

    Discontinuity_indicator = &
    	 & Discontinuity_indicator + &
    	 abs(sum( (face_val_at_quad(T, face) &
    	 - face_val_at_quad(T,face_2))*detwei*inflow ))

    Area = abs(sum(detwei*inflow))
    Inflow_integral = Inflow_integral + Area

  end subroutine Discontinuity_indicator_face

  function get_H(vertices, X) result (h)
    real, dimension(:,:), intent(in) :: X
    integer, dimension(:), intent(in) :: vertices
    real :: h
    !
    integer, dimension(size(X,1)+1) :: ind
    integer :: i,j,dim
    real :: a,b,c

    dim = size(X,1)
    
    do i = 1,dim+1
      ind(i) = vertices(i)
    enddo

    select case(dim)
    case (1)
       !Just take the difference
       h = abs(X(1,ind(1))-X(1,ind(2)))
    case (2)
       !Circumradius
       a = sqrt(sum((X(:,ind(2))-X(:,ind(1)))**2))
       b = sqrt(sum((X(:,ind(3))-X(:,ind(2)))**2))
       c = sqrt(sum((X(:,ind(1))-X(:,ind(3)))**2))
       h = a*b*c/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))
    case (3)
       !This should be circumradius too but I didn't code it
       h = 0.0
       do i = 1, size(X,ind(2))
          do j = 2, size(X,ind(2))
             h = max(h,sqrt(sum( (X(:,ind(i))-X(:,ind(j)))**2 )))
          end do
       end do
    case default
       FLExit('dont know that dimension.')
    end select
    
  end function get_H

  subroutine limit_slope_ele_dg_1d(ele, T, X)
    
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    
    integer, dimension(:), pointer :: neigh, T_ele
    real, dimension(mesh_dim(X)) :: ele_centre, ele_2_centre
    real :: ele_mean, ele_2_mean, dx
    real :: ele_slope, old_ele_slope, ele_2_slope
    integer :: ele_2, ni
    real, dimension(mesh_dim(X), ele_loc(X,ele)) :: X_val, X_val_2
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    
    X_val=ele_val(X, ele)
    T_val=ele_val(T, ele)
    

    ele_centre=sum(X_val,2)/size(X_val,2)

    ele_mean=sum(T_val)/size(T_val)
    
    dx=X_val(1,2)-X_val(1,1)
    
    ele_slope=(T_val(2)-T_val(1))/dx
    
    old_ele_slope=ele_slope

    neigh=>ele_neigh(T, ele)

    neighbourloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       
       ! These finding routines are outside the inner loop so as to allow
       ! for local stack variables of the right size in
       ! construct_add_diff_interface_dg.

       ele_2=neigh(ni)
           
       if (ele_2<=0) then
          ! External face.
          cycle
       end if

       X_val_2=ele_val(X, ele_2)
       T_val_2=ele_val(T, ele_2)

       ele_2_centre=sum(X_val_2,2)/size(X_val_2,2)

       ele_2_mean=sum(T_val_2)/size(T_val_2)
    
       ele_2_slope=(ele_2_mean-ele_mean)/sum(ele_2_centre-ele_centre)
       
       if (ele_slope*ele_2_slope<0.0) then
          ! Slope sign changes
          ele_slope=0.0
          exit neighbourloop
       end if

       ele_slope=sign(min(abs(ele_slope),abs(ele_2_slope)), ele_slope)

    end do neighbourloop

    if (old_ele_slope/=ele_slope) then
       
       ! Remove high order stuff here.
       T_ele=>ele_nodes(T,ele)

       call set(T, T_ele(1), ele_mean-0.5*dx*ele_slope)
       call set(T, T_ele(2), ele_mean+0.5*dx*ele_slope)

    end if

  end subroutine limit_slope_ele_dg_1d

  subroutine limit_vb(state, t)
    !Vertex-based (not Victoria Bitter) limiter from
    !Kuzmin, J. Comp. Appl. Math., 2010
    ! doi:10.1016/j.cam.2009.05.028
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: t
    !
    ! This is the limited version of the field, we have to make a copy
    type(scalar_field) :: T_limit, T_max, T_min
    type(mesh_type), pointer :: vertex_mesh
    ! counters
    integer :: ele, node
    ! local numbers
    integer, dimension(:), pointer :: T_ele
    ! gradient scaling factor
    real :: alpha
    ! local field values
    real, dimension(ele_loc(T,1)) :: T_val, T_val_slope, T_val_min,T_val_max
    real :: Tbar

    if (.not. element_degree(T%mesh, 1)==1 .or. continuity(T%mesh)>=0) then
      FLExit("The vertex based slope limiter only works for P1DG fields.")
    end if
    
    ! Allocate copy of field
    call allocate(T_limit, T%mesh,trim(T%name)//"Limited")
    call set(T_limit, T)
    
    ! returns linear version of T%mesh (if T%mesh is periodic, so is vertex_mesh)
    call find_linear_parent_mesh(state, T%mesh, vertex_mesh)

    call allocate(T_max, vertex_mesh, trim(T%name)//"LimitMax")
    call allocate(T_min, vertex_mesh, trim(T%name)//"LimitMin")
 
    call set(T_max, -huge(0.0))
    call set(T_min, huge(0.0))

    ! for each vertex in the mesh store the min and max values of the P1DG nodes directly surrounding it
    do ele = 1, ele_count(T)
       T_ele => ele_nodes(T,ele)
       T_val = ele_val(T,ele)
       Tbar = sum(T_val)/size(T_val)
       ! we assume here T is P1DG and vertex_mesh is linear
       assert( size(T_ele)==ele_loc(vertex_mesh,ele) )
       
       ! do maxes
       T_val_max = ele_val(T_max,ele)
       do node = 1, size(T_val)
          T_val_max(node) = max(T_val_max(node), Tbar)
       end do
       call set(T_max, ele_nodes(T_max, ele), T_val_max)
       
       ! do mins
       T_val_min = ele_val(T_min,ele)
       do node = 1, size(T_val)
          T_val_min(node) = min(T_val_min(node), Tbar)
       end do
       call set(T_min, ele_nodes(T_min,ele), T_val_min)
    end do

    ! now for each P1DG node make sure the field value is between the recorded vertex min and max
    ! this is done without changing the element average (Tbar)
    do ele = 1, ele_count(T)
       !Set slope factor to 1
       alpha = 1.
       !Get local node lists
       T_ele=>ele_nodes(T,ele)
       
       T_val = ele_val(T,ele)
       Tbar = sum(T_val)/size(T_val)
       T_val_slope = T_val - Tbar
       T_val_max = ele_val(T_max,ele)
       T_val_min = ele_val(T_min,ele)

       !loop over nodes, adjust alpha
       do node = 1, size(T_val)
	 !check whether to use max or min, and avoid floating point algebra errors due to round-off and underflow
	 if(T_val(node)>Tbar*(1.0+sign(1.0e-12,Tbar)) .and. T_val(node)-Tbar > tiny(0.0)*1e10) then
	   alpha = min(alpha,(T_val_max(node)-Tbar)/(T_val(node)-Tbar))
	 else if(T_val(node)<Tbar*(1.0-sign(1.0e-12,Tbar)) .and. T_val(node)-Tbar < -tiny(0.0)*1e10) then
	   alpha = min(alpha,(T_val_min(node)-Tbar)/(T_val(node)-Tbar))
	 end if
       end do

       call set(T_limit, T_ele, Tbar + alpha*T_val_slope)
    end do


    !Deallocate copy of field
    call set(T, T_limit)
    call halo_update(T)
    call deallocate(T_limit)
    call deallocate(T_max)
    call deallocate(T_min)

  end subroutine limit_vb

  subroutine limit_vb_julien(state, X, T)
    implicit none
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: X
    type(scalar_field), intent(inout) :: T
    !
    ! This is the limited version of the field, we have to make a copy
    type(scalar_field) :: DT_limit, T_proj, T_rem
    type(vector_field) :: X_proj
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(mesh_type), pointer :: vmesh
    
    ! counters
    integer :: deg, stat, ele, ele_2, face, ni, nj, nl, nf, nv, l
    logical :: limit
    
    ! local
    integer, dimension(:), pointer :: neigh, nodes, nodes_x, face_nodes, face_nodes_x, ele_neighbours
    integer, dimension(T%mesh%shape%numbering%vertices) :: ele_vertices 
    integer, dimension(X%mesh%faces%shape%loc) :: face_global
    integer, dimension(ele_loc(T,1)) :: flag
             
    real, dimension(ele_loc(T,1)) :: T_neigh_max, T_neigh_min, T_val, T_val_n, new_val
    real, dimension(X%dim,ele_loc(T,1)) :: X_proj_val
    real, dimension(X%dim,ele_loc(X,1)) :: X_val, X_val_n
    real, dimension(ele_loc(X,1)) :: base, Delta, Delta_low, T_proj_val, T_rem_val, T_low
    real, dimension(X%dim) :: X_neigh, X_mean
    real :: pos, neg, sum, T_mean, T_neigh, alpha

    deg=T%mesh%shape%degree
    ele_vertices = local_vertices(T%mesh%shape) 
    
    vfield=>extract_vector_field(state, "Velocity")
    vmesh=>vfield%mesh

    call allocate(X_proj, X%dim, T%mesh, "ProjectedMesh")
    call allocate(DT_limit, T%mesh, "DeltaLimitedScalar")
    call allocate(T_proj, vmesh, "ProjectedScalar")
    call allocate(T_rem, vmesh, "RemappedScalar")
    
    call remap_field(X, X_proj)
    call safe_set(state, T_proj, T)
    call zero(DT_limit)

    ! Find min/max values on each node
    do ele = 1, ele_count(T)

      neigh => ele_neigh(T,ele)      
      nodes => ele_nodes(T,ele)
      nodes_x => ele_nodes(X,ele)
      
      ! Exclude boundary elements
      if ( any(neigh <= 0) ) cycle
      
      T_val=ele_val(T,ele)
      T_proj_val=ele_val(T_proj,ele)
      T_mean=sum(T_proj_val)/size(T_proj_val)

      X_val=ele_val(X, ele)	  
      X_proj_val=ele_val(X_proj, ele)
      X_mean=sum(X_val,2)/size(X_val)
                  
      ! Find extrema at nodes
      flag=0
      T_neigh_max=-huge(0.0)
      T_neigh_min=huge(0.0)
      do nf = 1, size(neigh)
        ele_2=neigh(nf)
        face=ele_face(T, ele, ele_2)
        face_nodes=>face_local_nodes(T%mesh, face)
		
	nv=0	
	do nj = 1, size(face_nodes)
	  ni=face_nodes(nj)	  
	  if ( flag(ni) == 1 .or. .not.any(ele_vertices == ni) ) cycle
	  flag(ni)=1

	  nv=nv+1
          face_global=face_global_nodes(X%mesh, face)
          ele_neighbours=>node_neigh(X%mesh, face_global(nv))
          
          ! Get min/max limiting values
          do l = 1, size(ele_neighbours)
            if ( ele_neighbours(l) == ele ) cycle
            T_val_n=ele_val(T,ele_neighbours(l))
            T_neigh_max(ni)=max(T_neigh_max(ni),maxval(T_val_n))
            T_neigh_min(ni)=min(T_neigh_min(ni),minval(T_val_n))
     	  enddo
        enddo
	
	do nj = 1, size(face_nodes)
	  ni=face_nodes(nj)	  
	  if ( flag(ni) == 0 .and. .not.any(ele_vertices == ni) ) then
	    flag(ni)=1
	    T_neigh_max(ni)=maxval(T_neigh_max(face_nodes))
	    T_neigh_min(ni)=minval(T_neigh_min(face_nodes))
	  endif
	enddo
      enddo
     
      ! Limiter indicator based on higher-order terms
      limit=.false.
      nj=0
      do ni = 1, size(T_val)
        if ( abs(T_val(ni)-T_mean)/abs(T_mean) <= tolerance ) cycle
      
     	if (T_neigh_max(ni) < T_neigh_min(ni)) then
	  T_neigh_max(ni)=T_val(ni); T_neigh_min(ni)=T_val(ni)
	endif
	
	if (deg > 1) then
	  if ( T_val(ni)-T_mean < (0.333*(limit_factor-1)+1)*(T_neigh_min(ni)-T_mean) &
	  .or. T_val(ni)-T_mean > (0.333*(limit_factor-1)+1)*(T_neigh_max(ni)-T_mean) ) &
	     limit=.true.
	else
	  limit=.true.
	endif
      enddo
      
      if ( .not.limit ) cycle
      
      nj=0
      do ni = 1, size(nodes)
        if ( any(ele_vertices == ni) ) then
          nj=nj+1
	  Delta_low(nj)=minmod(T_neigh_max(ni)-T_mean,T_neigh_min(ni)-T_mean,T_proj_val(nj)-T_mean,limit_factor)
        endif
      enddo   
      
      if (sum(T_mean+Delta_low) == sum(T_proj_val)) cycle
      
      ! Enforce mass conservation on low order limiter
      if (abs(sum(Delta_low))>1000.0*epsilon(0.0)) then
        pos=sum(max(0.0, Delta_low))
        neg=sum(max(0.0, -Delta_low))
      
        Delta_low = min(1.0,neg/pos)*max(0.0,Delta_low)   &
              - min(1.0,pos/neg)*max(0.0,-Delta_low)
      end if
            
      ! Apply limiter to vertices
      nj=0
      do ni = 1, size(nodes)
        if ( any(ele_vertices == ni) ) then
          nj=nj+1
	  new_val(ni) = T_mean + Delta_low(nj)
	  T_low(nj) = new_val(ni)
        endif
      enddo   
	
      ! Apply limiter to mid-face nodes
      do nf = 1, size(neigh)
        ele_2=neigh(nf)
        face=ele_face(T, ele, ele_2)
        face_nodes=>face_local_nodes(T%mesh, face)
		
	do nj = 1, size(face_nodes)
	  ni=face_nodes(nj)	  
          if ( .not.any(ele_vertices == ni) ) then
            base=create_basis(face,X_proj_val(:,ni),X_val,X_mean)
	    new_val(ni) = sum(base*T_low)
          endif
	enddo
      enddo 
      
      call set(DT_limit, ele_nodes(T,ele), new_val-T_val)
       
    enddo 
    ewrite_minmax(DT_limit)

    !Deallocate copy of field
    call addto(T, DT_limit)
    call halo_update(T)
 
    if (trim(T%name)=="PotentialTemperature") then  
      sfield=>extract_scalar_field(state, "DeltaLimitPotentialTemperature", stat)
      if (stat==0) call set(sfield, DT_limit)
    endif
    
    call deallocate(X_proj)
    call deallocate(T_proj)
    call deallocate(DT_limit)

contains

    function minmod(phi_max,phi_min,lim,limit_factor)
      real :: phi_max, phi_min, lim
      real :: minmod
      real :: limit_factor
      
      if ( sign(1.0,phi_max) == sign(1.0,phi_min) .and. sign(1.0,phi_max) == sign(1.0,lim) ) then
        minmod=sign(1.0,lim)*min(limit_factor*max(abs(phi_max),abs(phi_min)),abs(lim))
      else if ( sign(1.0,phi_max) /= sign(1.0,phi_min) .and. sign(1.0,phi_max) == sign(1.0,lim) ) then
        minmod=min(limit_factor*phi_max,lim)
      else if ( sign(1.0,phi_max) /= sign(1.0,phi_min) .and. sign(1.0,phi_min) == sign(1.0,lim) ) then
        minmod=max(limit_factor*phi_min,lim)
      else
	minmod=0.0
      endif
      
      return
    end function minmod  

    function sweby(phi_max,phi_min,lim,limit_factor)
      real :: phi_max, phi_min, lim, dx2
      real :: minmod, sweby
      real :: limit_factor
      
      if (sign(1.0,phi_max) == sign(1.0,phi_min)) then
        minmod=min(limit_factor*max(phi_max,phi_min),abs(lim))
      else 
        minmod=0.0
      endif
      
      if ( minmod == abs(lim) ) return
      
      if (sign(1.0,phi_max) == sign(1.0,phi_min) .and. sign(1.0,phi_max) == sign(1.0,lim)) then
        sweby=min(limit_factor*abs(phi_max),limit_factor*abs(phi_min),abs(lim))
        sweby=sign(1.0,lim)*max(0.,sweby,minmod)
      else 
        sweby=0.0
      endif
      
      return
    end function sweby  

    function create_basis(face,X_proj,X_val,X_mean)
      real, dimension(:), intent(in) :: X_mean,X_proj
      real, dimension(:,:), intent(in) :: X_val
      real, dimension(size(X_val,2)) :: create_basis
      integer :: face
      
      integer, dimension(:), pointer :: face_nodes
      integer :: info, ni
      real, dimension(size(X_val,1),size(X_val,1)) :: phimat
      real, dimension(size(X_val,1)) :: phi
      real :: dx_max

      face_nodes=>face_local_nodes(X%mesh, face)

      do ni=1,size(face_nodes)
        phimat(:,ni)=X_val(:,face_nodes(ni))-X_mean
      enddo
      phi=X_proj-X_mean

      !Solve for basis coefficients
      call solve(phimat,phi,info)

      create_basis=0.0
      do ni=1,size(face_nodes)
        create_basis(face_nodes(ni))=phi(ni)
      enddo   

    end function create_basis

  end subroutine limit_vb_julien
  
  subroutine limit_fpn(state, t)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: t
    
    type(scalar_field), pointer :: limiting_t, lumped_mass
    type(scalar_field) :: lowerbound, upperbound, inverse_lumped_mass
    type(csr_matrix), pointer :: mass

    type(csr_sparsity), pointer :: eelist
    integer :: ele, i, j, k, row, column
    integer :: rows, columns
    real :: node_max, node_min, extra_val, extra_val2

    integer, dimension(:), pointer :: nodelist, faces, neighbouring_ele_nodes
    integer, dimension(:), allocatable :: face_nodes, neighbouring_nodes
    integer :: neighbouring_face, neighbouring_ele
    logical, save :: first=.true.
    logical :: midpoint, extrapolate, pre_dist_mass

    real :: beta=1.0, mean_val
    type(vector_field), pointer :: position
    type(vector_field) :: dg_position
    real, dimension(:), allocatable :: e_vec_1
    real, dimension(:,:,:), allocatable :: dt_t
    real, dimension(:,:), allocatable :: grad_t
    real :: grad, e_dist

    real, dimension(ele_loc(t,1)) :: weight, tracer_val
    logical, dimension(ele_loc(t,1)) :: nweight, pweight
    real :: nodeval, nodemin, nodemax, adjust

    integer, dimension(:,:,:), allocatable, save :: nodes_array
!     real, dimension(2,4) :: local_values
!     real, dimension(2) :: line_max, line_min
    real, dimension(:,:), allocatable :: local_values
    real, dimension(:), allocatable :: line_max, line_min
    integer :: node, adjacent_node, local_face

    integer, dimension(face_loc(t,1)) :: fnodes, neighbouring_face_nodes

    type(vector_field), pointer :: u
    real :: mod_u, tol
    real, dimension(2) :: u_node, e_vec_2, e_vec_3
    logical :: upwind=.false.
    real :: cor1, cor2, cor3, cor4
    integer :: problem_dimension, values

    call get_option('/geometry/dimension', problem_dimension)

    rows=problem_dimension ! The number of 'lines' to look along
    columns=3 ! The number of nodes on each line. Should always be three
    if (upwind) then
      values=columns+1
    else
      values=columns
    end if
    tol=0.25

    allocate(local_values(rows,values),line_max(rows),line_min(rows))

    midpoint=have_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter::FPN/mid-point_scheme")
    if (midpoint) then
      call get_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter::FPN/mid-point_scheme/beta", beta, default=1.0)
      extrapolate=have_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
           &discontinuous_galerkin/slope_limiter::FPN/mid-point_scheme/extrapolate")
    end if

!     pre_dist_mass=have_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
!          &discontinuous_galerkin/slope_limiter::FPN/pre_distribute_mass")

    mass => get_mass_matrix(state, t%mesh)
    lumped_mass => get_lumped_mass(state, t%mesh)
    call allocate(inverse_lumped_mass, lumped_mass%mesh, "InverseLumpedMass")
    inverse_lumped_mass%val = 1.0/lumped_mass%val

    limiting_t => extract_scalar_field(state, trim(t%name))

!     eelist => extract_eelist(t%mesh)

    call allocate(lowerbound, t%mesh, "LowerBound")
    call allocate(upperbound, t%mesh, "UpperBound")
    call zero(lowerbound); call zero(upperbound)

    allocate (neighbouring_nodes(ele_loc(limiting_t,1)))

!     allocate (face_nodes(face_loc(limiting_t,1)), neighbouring_nodes(face_loc(limiting_t,1)))

    if (extrapolate) then
      position => extract_vector_field(state, "Coordinate")
      call allocate(dg_position, position%dim, t%mesh, name="DG_Coordinate")
      call remap_field(position, dg_position)
      allocate (e_vec_1(position%dim),dt_t(ele_loc(limiting_t, 1), ele_ngi(limiting_t, 1), mesh_dim(limiting_t)))
      allocate (grad_t(mesh_dim(limiting_t), limiting_t%mesh%shape%ngi))
    end if

    if (upwind) then
      u => extract_vector_field(state, "Velocity")
    end if

    ! Loop to construct an array containing the global node numbers required to compute the limiting values
    ! at a node i. Only evaluated on the first timestep (and after every adapt for adaptive runs).
    if (first) then
      allocate (nodes_array(node_count(t),rows,columns))
      first=.false.
      do node=1,node_count(limiting_t)
        ele=node_ele(limiting_t, node)
        nodelist => ele_nodes(limiting_t, ele)
        faces => ele_faces(limiting_t, ele)
        row=0
        do i=1,size(nodelist)
          if (nodelist(i)==node) cycle
          row=row+1
          fnodes=face_global_nodes(limiting_t,faces(i))
          do j=1,size(fnodes)
            if (fnodes(j)==node) adjacent_node=j
          end do
          neighbouring_face = face_neigh(limiting_t, faces(i))
          neighbouring_face_nodes = face_global_nodes(limiting_t,neighbouring_face)
!           secnd_val=neighbouring_face_nodes(adjacent_node) ! 2nd node we want
          local_face=local_face_number(limiting_t,neighbouring_face)
          neighbouring_ele = face_ele(limiting_t, neighbouring_face)
          neighbouring_nodes = ele_nodes(limiting_t, neighbouring_ele)
!           thrid_val=neighbouring_nodes(local_face)
          nodes_array(node,row,1)=nodelist(i)
          nodes_array(node,row,2)=neighbouring_face_nodes(adjacent_node)
          nodes_array(node,row,3)=neighbouring_nodes(local_face)
        end do
      end do
    end if

    ! Loop through the nodes and calculate the bounds for each node
    do node=1,node_count(limiting_t)
      ! Calculate the av. value of the tracer within the element
      ele=node_ele(limiting_t, node)
      nodelist => ele_nodes(limiting_t, ele)
      do i=1, size(nodelist)
        tracer_val(i)=node_val(limiting_t,nodelist(i))
      end do
      mean_val=sum(tracer_val)/float(size(nodelist))
      ! Get the values needed for calculating the bounds
      do row=1,rows
        do column=1,columns
          local_values(row,column)=node_val(limiting_t, nodes_array(node,row,column))
        end do
        ! Adjust values depending on options
        if (midpoint.and.(.not.extrapolate)) then
          local_values(row,3) = (1.0-beta)*local_values(row,2)+beta*local_values(row,3)
        else if (midpoint.and.extrapolate) then
          local_values(row,3) = (1.0-beta)*local_values(row,2)+beta*local_values(row,3)
          ! Extrapolate using the gradients of the neighbouring element to form our extra value
                
          ! 1st, work out the direction in which we want to extrapolate, e_vec_1
          e_vec_1=node_val(dg_position,node)-node_val(dg_position,nodes_array(node,row,1))
          ! Work out the distance to exprapolate
          e_dist=sqrt(sum(e_vec_1(:)**2))
          ! Turn this into a unit vector
          e_vec_1=e_vec_1/e_dist

          call transform_to_physical(dg_position, node_ele(limiting_t, nodes_array(node,row,2)), &
                                     ele_shape(limiting_t,node_ele(limiting_t, nodes_array(node,row,2))), dshape=dt_t)

          grad_t=ele_grad_at_quad(limiting_t, node_ele(limiting_t, node), dt_t)

          ! Calculate the gradient in the desired direction
          ! Note that grad_t will be the same at all gauss points in the linear element case
          grad=dot_product(grad_t(1,:),e_vec_1)

          local_values(row,4) = local_values(row,2)+beta*grad*e_dist 
        end if
        if (upwind) then
          u_node=node_val(u, node)
          mod_u=sqrt(dot_product(u_node,u_node))
          u_node=u_node/mod_u
          e_vec_2=e_vec_1
          e_vec_1=-e_vec_1
          e_vec_3=node_val(dg_position,nodes_array(node,row,3))-node_val(dg_position,node)
          e_vec_3=e_vec_3/sqrt(dot_product(e_vec_3,e_vec_3))
          cor1=dot_product(u_node,e_vec_1)
          cor2=dot_product(u_node,e_vec_2)
          cor3=dot_product(u_node,e_vec_3)
          cor4=dot_product(u_node,e_vec_3)
          if (cor1<tol) then
            local_values(row,1)=node_val(limiting_t,node)
          end if
          if (cor2<tol) then
            local_values(row,2)=node_val(limiting_t,node)
          end if
          if (cor3<tol) then
            local_values(row,3)=node_val(limiting_t,node)
          end if
          if (cor4<tol) then
            local_values(row,4)=node_val(limiting_t,node)
          end if
        end if
      end do
      ! Calculate and set the bounds
      line_max(1)=maxval(local_values(1,:))
      line_min(1)=minval(local_values(1,:))
      line_max(2)=maxval(local_values(2,:))
      line_min(2)=minval(local_values(2,:))
      node_max=minval(line_max)
      node_min=maxval(line_min)
      if (node_max<mean_val) node_max=mean_val
      if (node_min>mean_val) node_min=mean_val
      call set(lowerbound, node, node_min)
      call set(upperbound, node, node_max)
    end do

    call bound_field_diffuse(t, upperbound, lowerbound, mass, lumped_mass, inverse_lumped_mass)

    deallocate (neighbouring_nodes)
    call deallocate(inverse_lumped_mass)
    call deallocate(upperbound)
    call deallocate(lowerbound)
    if (extrapolate) then
      call deallocate(dg_position)
      deallocate (e_vec_1, dt_t, grad_t)
    end if

  end subroutine limit_fpn

end module slope_limiters_dg
