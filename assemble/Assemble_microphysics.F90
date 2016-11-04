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
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.f
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module assemble_microphysics
  !!< This module implements a warm cloud microphysics model in Fluidity

  use spud
  use state_module
  use fields
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN
  use slope_limiters_dg
  use microphysics
  use futils, only : int2str
  use fldebug

  implicit none

  public :: assemble_split_microphysics

  contains

  subroutine assemble_split_microphysics (state, dt)
  
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: dt
    character(len=OPTION_PATH_LEN) :: mom_path, micro_path, micro_name
    logical :: have_ndrop, have_nrain, have_cold, have_ccn
    integer :: j, nmom
        
    ewrite(1,*) 'Assemble split microphysics'
    
    do j = 1, size(state)
    
    mom_path = '/material_phase['//int2str(j-1)//']/cloud_microphysics'
    if (have_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics')) then
      nmom = 2
      call get_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics/name', &
    		      micro_name)
      mom_path=trim(mom_path)//'/fortran_microphysics/two_moment_microphysics::'//trim(micro_name)
    else if (have_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics')) then
      nmom = 1
      call get_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics/name', &
    		      micro_name)
      mom_path=trim(mom_path)//'/fortran_microphysics/one_moment_microphysics::'//trim(micro_name)
    endif
    
    have_cold=have_option(trim(mom_path)//'/cold_microphysics')
    have_ndrop=have_option(trim(mom_path)//'/scalar_field::Ndrop/prognostic')
    have_nrain=have_option(trim(mom_path)//'/scalar_field::Nrain/prognostic')
    have_ccn=have_option(trim(mom_path)//'/scalar_field::CCN/prognostic')

    micro_path='/material_phase['//int2str(j-1)//']/cloud_microphysics'
    if (have_option(trim(micro_path)//'/time_integration::Splitting') .or. &
        have_option(trim(micro_path)//'/time_integration::Strang') ) then
      call advance_microphysics(state(j),dt,"Qdrop")
      call advance_microphysics(state(j),dt,"Qrain")
      if (have_ndrop) call advance_microphysics(state(j),dt,"Ndrop")
      if (have_nrain) call advance_microphysics(state(j),dt,"Nrain")  
      if (have_ccn) call advance_microphysics(state(j),dt,"CCN")

      if (have_cold) then
        call advance_microphysics(state(j),dt,"Qice")
        call advance_microphysics(state(j),dt,"Qgrau")
        call advance_microphysics(state(j),dt,"Qsnow")
        call advance_microphysics(state(j),dt,"Nice")
        if (nmom>1) call advance_microphysics(state(j),dt,"Ngrau")    
        if (nmom>1) call advance_microphysics(state(j),dt,"Nsnow")    
      endif
      	      
      if (has_scalar_field(state(j),"Temperature")) then
         call advance_microphysics(state(j),dt,"Temperature")
      end if
      if (has_scalar_field(state(j),"PotentialTemperature")) then
         call advance_microphysics(state(j),dt,"PotentialTemperature")
      end if
      if (has_scalar_field(state(j),"ConservedPotentialTemperature")) then
         call advance_microphysics(state(j),dt,"ConservedPotentialTemperature")
      end if
      if (has_scalar_field(state(j),"VapourWaterQ")) then
         call advance_microphysics(state(j),dt,"VapourWaterQ")
      end if
    endif
    
    enddo

  end subroutine assemble_split_microphysics
  
  subroutine advance_microphysics(state,dt,name)
    type(state_type),intent(inout) :: state
    character(len=*) :: name
    real :: dt
    type(scalar_field), pointer :: sfield, mfield
    type(vector_field), pointer :: U, X
    integer :: stat, limiter
    logical :: limit_after_advance
    character(len=OPTION_PATH_LEN) :: micro_path
        
    ewrite(2,*) 'Advance split microphysics for scalar: '//trim(name)

    sfield=>extract_scalar_field(state,name,stat=stat)
    if (stat /= 0) then
       ewrite(3,*) ('No field '//name//' in state to store!')
       return
    end if

    mfield=>extract_scalar_field(state,trim(name)//'MicrophysicsSource',stat=stat)
    if (stat /= 0) then
       ewrite(3,*) ('No field '//name//' in state to store!')
       return
    end if

    !Advance sfield
    call addto(sfield,mfield,scale=dt)
    
    !Store old
    call store_microphysics_source(state,trim(mfield%name))
    
    !Limit if necessary
    micro_path='/material_phase[0]/cloud_microphysics'
    if (have_option(trim(micro_path)//'/time_integration::Splitting/limit_after_advance') .or. &
        have_option(trim(micro_path)//'/time_integration::Strang/limit_after_advance').and. &
         &have_option(trim(sfield%option_path)//"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/slope_limiter")) then
       U=>extract_vector_field(state, "NonlinearVelocity", stat)
       if (stat/=0) U=>extract_vector_field(state, "Velocity", stat)
       X=>extract_vector_field(state, "Coordinate", stat)       
       call limit_slope_dg(sfield, U, X, state, limiter, not_initialised=.true.)
    endif

  end subroutine advance_microphysics
  
end module assemble_microphysics
