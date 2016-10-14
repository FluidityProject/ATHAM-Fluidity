!    Copyright (C) 2007 Imperial College London and others.
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
module filter_diagnostics

  use field_options
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use parallel_tools
  use spud
  use state_module
  use diagnostic_fields, only: safe_set
  use smoothing_module
  
  private
  public :: calculate_horizontal_filter, calculate_sponge_coefficient_vector, calculate_sponge_coefficient_scalar
  
  interface calculate_horizontal_filter
    module procedure calculate_horizontal_filter_scalar, calculate_horizontal_filter_vector
  end interface

contains

  subroutine calculate_horizontal_filter_scalar(state,s_field,f_field,nits,alpha,base_path,horizontal,vertical)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), intent(inout) :: f_field
    logical, intent(in), optional :: horizontal, vertical

    integer, optional :: nits
    real, optional ::  alpha
    character(len = OPTION_PATH_LEN), optional :: base_path
    
    type(vector_field), pointer :: positions

    integer :: its
    character(len = OPTION_PATH_LEN) :: direction
    real :: a=1.e-6, r=1.e-10, diff_max, rhs_max
    type(scalar_field) :: rhs, out
    
    assert(element_count(s_field) == element_count(f_field))
    
    positions => extract_vector_field(state, "Coordinate")

    call allocate(rhs, positions%mesh, "RHS")
    call allocate(out, positions%mesh, "OutField")
    call safe_set(state,rhs,s_field)

    if (.not.present(base_path)) base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"
    if (.not.present(nits)) call get_option(trim(base_path) // "/number_of_iterations",nits)
    if (.not.present(alpha)) call get_option(trim(base_path) // "/alpha",alpha)

    if (present(horizontal)) then
      direction='horizontal'
    else
      direction='vertical'
    endif

    lexit=.false.
    do its=1,nits
       call horizontal_vertical_smooth_scalar(rhs,positions,out,alpha,base_path,direction)
       
       diff_max=maxval(abs(rhs%val))-maxval(abs(out%val))
       rhs_max=maxval(abs(rhs%val))
       call allmax(diff_max)
       call allmax(rhs_max)

       call set(rhs,out)
       call halo_update(rhs)       
       
       ewrite(2,*) 'Vertical smoothing '//int2str(its)
       if ( abs(diff_max) < a*rhs_max+r ) exit
    end do

    call safe_set(state,f_field,rhs)
    call deallocate(rhs)
    call deallocate(out)

  end subroutine calculate_horizontal_filter_scalar

  subroutine calculate_horizontal_filter_vector(state,v_field,f_field,nits,alpha,base_path,horizontal,vertical)

    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    type(vector_field), intent(out) :: f_field
    logical, intent(in), optional :: horizontal, vertical

    integer, optional :: nits
    real, optional ::  alpha
    character(len = OPTION_PATH_LEN), optional :: base_path
    
    type(vector_field), pointer :: positions

    integer :: its, i
    character(len = OPTION_PATH_LEN) :: direction

    type(scalar_field) :: rhs, f_rhs
    
    positions => extract_vector_field(state, "Coordinate")

    call allocate(rhs, v_field%mesh, "RHS")
    call allocate(f_rhs, v_field%mesh, "FilteredRHS")

    if (.not.present(base_path)) base_path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    if (.not.present(nits)) call get_option(trim(base_path) // "/number_of_iterations",nits)
    if (.not.present(alpha)) call get_option(trim(base_path) // "/alpha",alpha)

    if (present(horizontal)) then
      direction='horizontal'
    else
      direction='vertical'
    endif

    do i = 1, v_field%dim

      call set(rhs,v_field,i)
      do its=1,nits
         call horizontal_vertical_smooth_scalar(rhs,positions,f_rhs,alpha,base_path,direction)
         call set(rhs,f_rhs)
      end do
      call set(f_field,i,f_rhs)

    enddo

    call deallocate(rhs)
    call deallocate(f_rhs)

  end subroutine calculate_horizontal_filter_vector
    
  subroutine calculate_sponge_coefficient_scalar(s_field, X, option_path)
    type(scalar_field), intent(inout) :: s_field
    type(vector_field), intent(inout) :: X
    character(len=OPTION_PATH_LEN), intent(in) :: option_path    
    
    integer :: i, stat, ele, dim
    real  :: z_base, x_base_r, x_base_l, y_base_r, y_base_l, z_max, x_max, x_min, y_max, y_min, xi, X_val
    real  :: absorption_valx, absorption_valy, absorption_valz
    type(vector_field) :: X_local
    character(len=OPTION_PATH_LEN) :: absorption_path
    
    call zero(s_field)
  
    absorption_path = trim(option_path)//&
    	"/prognostic/scalar_field::Absorption/diagnostic/algorithm::atmosphere_forcing_scalar"

    if (have_option(trim(absorption_path)//"/sponge_layer_scalar_absorption")) then
      
      call allocate (X_local,X%dim,s_field%mesh,"LocalMesh")
      call remap_field(X,X_local)
      
      if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/x_sponge_right')) &
    	  call get_option (trim(absorption_path)//'/sponge_layer_scalar_absorption/x_sponge_right',x_base_r)
      if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/x_sponge_left')) &
    	  call get_option (trim(absorption_path)//'/sponge_layer_scalar_absorption/x_sponge_left',x_base_l)
      if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/y_sponge_right')) &
    	  call get_option (trim(absorption_path)//'/sponge_layer_scalar_absorption/y_sponge_right',y_base_r)
      if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/y_sponge_left')) &
    	  call get_option (trim(absorption_path)//'/sponge_layer_scalar_absorption/y_sponge_left',y_base_l)
      if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/z_sponge')) &
    	  call get_option (trim(absorption_path)//'/sponge_layer_scalar_absorption/z_sponge',z_base)

      dim=X_local%dim
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

      do ele=1,node_count(s_field)
        absorption_valx=0.
        absorption_valy=0.
        absorption_valz=0.
        
        X_val=node_val(X_local,1,ele)
        if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/x_sponge_right') .and. X_val >= x_base_r) then
           xi=(X_val - x_base_r)/(x_max - x_base_r)
           if (xi >= 0. .and. xi < 0.75) then
             absorption_valx=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
           else if (xi > 0.75 .and. xi <= 1.) then
             absorption_valx=1.
           else
             absorption_valx=0.
           endif
    	endif
        
        if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/x_sponge_left') .and. X_val <= x_base_l) then
           xi=(x_base_l - X_val)/(x_base_l - x_min)  
           if (xi >= 0. .and. xi < 0.75) then
             absorption_valx=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
           else if (xi > 0.75 .and. xi <= 1.) then
             absorption_valx=1.
           else
             absorption_valx=0.
           endif
        endif
        
        if (dim > 2) then
          X_val=node_val(X_local,2,ele)
          if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/y_sponge_right') .and. X_val >= y_base_r) then
             xi=(X_val - y_base_r)/(y_max - y_base_r)
             if (xi >= 0. .and. xi < 0.75) then
               absorption_valy=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
             else if (xi > 0.75 .and. xi <= 1.) then
               absorption_valy=1.
             else
               absorption_valy=0.
             endif
    	  endif
          if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/y_sponge_left') .and. X_val <= y_base_l) then
             xi=(y_base_l - X_val)/(y_base_l - y_min)
             if (xi >= 0. .and. xi < 0.75) then
               absorption_valy=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
             else if (xi > 0.75 .and. xi <= 1.) then
               absorption_valy=1.
             else
               absorption_valy=0.
             endif
          endif
        endif
        
        X_val=node_val(X_local,X_local%dim,ele)
        if (have_option(trim(absorption_path)//'/sponge_layer_scalar_absorption/z_sponge') .and. X_val >= z_base) then
           xi=(X_val - z_base)/(z_max - z_base)
           if (xi >= 0. .and. xi < 0.75) then
             absorption_valz=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
           else if (xi > 0.75 .and. xi <= 1.) then
             absorption_valz=1.
           else
             absorption_valz=0.
           endif
    	endif
        
        absorption_valx = 1. - min(sqrt(absorption_valx**2. + absorption_valy**2. + absorption_valz**2.),1.)
        call addto(s_field, ele, absorption_valx)
      enddo
    
      call deallocate (X_local)
    end if 

  end subroutine calculate_sponge_coefficient_scalar
    
  subroutine calculate_sponge_coefficient_vector(s_field, X, option_path)
    type(scalar_field), intent(inout) :: s_field
    type(vector_field), intent(inout) :: X
    character(len=OPTION_PATH_LEN), intent(in) :: option_path
    
    integer :: j, stat, ele, dim
    real  :: z_base, x_base_r, x_base_l, y_base_r, y_base_l, z_max, x_max, x_min, y_max, y_min, xi, X_val
    real :: absorption_valx, absorption_valz, absorption_valy
    type(vector_field) :: X_local
    character(len=OPTION_PATH_LEN) :: absorption_path
  
    call zero(s_field)

    absorption_path = trim(option_path)//  &
    	 "/prognostic/vector_field::Absorption/diagnostic/algorithm::atmosphere_forcing_vector"
    
    if (have_option(trim(absorption_path)//"/sponge_layer_velocity_absorption")) then
    
      call allocate(X_local,X%dim,s_field%mesh,"LocalMesh")
      call remap_field(X,X_local)
    
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

      dim=X_local%dim
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

      do ele=1,node_count(s_field)
        absorption_valx=0.
        absorption_valy=0.
        absorption_valz=0.
        
        X_val=node_val(X_local,1,ele)
        if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_right') .and. X_val >= x_base_r) then
           xi=(X_val - x_base_r)/(x_max - x_base_r)
           if (xi >= 0. .and. xi < 0.75) then
             absorption_valx=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
           else if (xi > 0.75 .and. xi <= 1.) then
             absorption_valx=1.
           else
             absorption_valx=0.
           endif
    	endif
        
        if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/x_sponge_left') .and. X_val <= x_base_l) then
           xi=(x_base_l - X_val)/(x_base_l - x_min)
           if (xi >= 0. .and. xi < 0.75) then
             absorption_valx=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
           else if (xi > 0.75 .and. xi <= 1.) then
             absorption_valx=1.
           else
             absorption_valx=0.
           endif	    
        endif
        
        if (dim > 2) then
          X_val=node_val(X_local,2,ele)
          if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_right') .and. X_val >= y_base_r) then
             xi=(X_val - y_base_r)/(y_max - y_base_r)
             if (xi >= 0. .and. xi < 0.75) then
               absorption_valy=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
             else if (xi > 0.75 .and. xi <= 1.) then
               absorption_valy=1.
             else
               absorption_valy=0.
             endif
    	  endif
          
          if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/y_sponge_left') .and. X_val <= y_base_l) then
             xi=(y_base_l - X_val)/(y_base_l - y_min)
             if (xi >= 0. .and. xi < 0.75) then
               absorption_valy=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
             else if (xi > 0.75 .and. xi <= 1.) then
               absorption_valy=1.
             else
               absorption_valy=0.
             endif
    	  endif
        endif
        
        X_val=node_val(X_local,X_local%dim,ele)
        if (have_option(trim(absorption_path)//'/sponge_layer_velocity_absorption/z_sponge') .and. X_val >= z_base) then
           xi=(X_val - z_base)/(z_max - z_base)
           if (xi >= 0. .and. xi < 0.75) then
             absorption_valz=sin(3.1416/2.*xi/0.75)*sin(3.1416/2.*xi/0.75)
           else if (xi > 0.75 .and. xi <= 1.) then
             absorption_valz=1.
           else
             absorption_valz=0.
           endif
        endif  

        absorption_valx = 1. - min(sqrt(absorption_valx**2. + absorption_valy**2. + absorption_valz**2.),1.)
        call addto(s_field, ele, absorption_valx)
      enddo
    
      call deallocate (X_local)
    end if 

  end subroutine calculate_sponge_coefficient_vector


end module filter_diagnostics
