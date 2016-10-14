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

module diagnostic_fields_eos

  use fields
  use fldebug
  use spud
  use field_options
  use global_parameters, only:FIELD_NAME_LEN, current_time, OPTION_PATH_LEN
  use diagnostic_fields, only: safe_set
  use state_module
  use equation_of_state
  use microphysics
  
  implicit none
  
  private
  public :: calculate_diagnostic_variable_eos  
  
  interface calculate_diagnostic_variable_eos
    module procedure diagnostic_variable_eos_single, 		&
        diagnostic_variable_eos_multiple_non_indexed, 		&
	diagnostic_variable_eos_multiple_indexed
  end interface calculate_diagnostic_variable_eos
  
contains
  
  subroutine diagnostic_variable_eos_single (state, s_field, s_field_name)
  
    type(state_type), intent(inout) :: state
    character(len = *), intent(in) :: s_field_name
    type(scalar_field), intent(inout) :: s_field
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call calculate_diagnostic_variable_eos(states, 1, s_field, s_field_name)
    state = states(1)
  
  end subroutine diagnostic_variable_eos_single
  
  subroutine diagnostic_variable_eos_multiple_non_indexed(states, s_field, s_field_name)
  
    type(state_type), dimension(:), target, intent(inout) :: states
    type(scalar_field), intent(inout) :: s_field
    character(len = *), intent(inout) :: s_field_name
    
    call calculate_diagnostic_variable_eos(states, 1, s_field, s_field_name)
  
  end subroutine diagnostic_variable_eos_multiple_non_indexed

  subroutine diagnostic_variable_eos_multiple_indexed (states, state_index, s_field, s_field_name)
    !!< Calculate the specified scalar diagnostic field s_field_name from state
    !!< and return the field in s_field.

    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: s_field_name
    type(scalar_field), intent(inout) ::s_field

    integer :: stat
    type(state_type), pointer :: state => null()
    type(mesh_type), pointer :: mesh
    type(scalar_field) :: tmp
    character(len=FIELD_NAME_LEN) :: mesh_name

    state => states(state_index)
    call get_option(trim(s_field%option_path)//'/diagnostic/mesh/name',mesh_name,stat)
    if (stat==0) mesh=>extract_mesh(state, trim(mesh_name))    

    select case(s_field_name)
    
      case("ExnerPressure")
        call get_exner_pressure(state, s_field, stat)
	
      case("SoundSpeed")
        call compressible_eos(state, sound_speed=s_field)
	
      case("VapourWaterQ")
        call get_water_vapour(state, s_field, stat)

      case("Saturation")
        tmp=extract_scalar_field(state, "VapourWaterQ", stat=stat)
	if (stat /= 0) then
          call allocate(tmp,s_field%mesh,"WaterVapour")
	  call zero(tmp)
          call get_water_vapour(state, tmp, stat)
        else
	  call incref(tmp)
	endif
	
	call compressible_eos(state, saturation=s_field)
	
	call invert(s_field)
	call scale(s_field,tmp)
	call addto(s_field,-1.)
	call scale(s_field,100.)
	
	call deallocate(tmp)

      case("Reflectivity")
        call get_reflectivity(state, s_field, stat)

      case default
        ! Nothing to be done: All other variables are already treated in Diagnostic_Fields.

    end select
    
    ewrite_minmax(s_field)

  end subroutine diagnostic_variable_eos_multiple_indexed

  subroutine get_exner_pressure(state, s_field, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    integer, intent(out), optional :: stat 
    
    real :: p_0
    type(scalar_field) :: pressure, c_p, c_v

    call allocate(pressure,s_field%mesh,"LocalPressure")
    call allocate(c_p,s_field%mesh,"LocalCP")
    call allocate(c_v,s_field%mesh,"LocalCV")
    
    call zero(s_field)

    call get_option('/material_phase[0]/equation_of_state/compressible/giraldo/reference_pressure',p_0, default=1.0e5)
    
    call compressible_eos(state, pressure=pressure, qc_p=c_p, qc_v=c_v)
    
    call addto(s_field,pressure,scale=1./p_0)
    s_field%val=(s_field%val)**((c_p%val-c_v%val)/c_p%val)
    
    call deallocate(pressure)
    call deallocate(c_p)
    call deallocate(c_v)

  end subroutine get_exner_pressure

  subroutine get_water_vapour(state, s_field, stat)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    integer, intent(out), optional :: stat 
    
    type(scalar_field), pointer :: QT,QC,QR,QI,QG,QS
    integer :: ele

    QT=>extract_scalar_field(state, "TotalWaterQ", stat)
    
    if (stat == 0) then
    
    call set (s_field, QT)
    
    if (has_scalar_field(state,'Qdrop')) then
      QC=>extract_scalar_field(state, "Qdrop", stat)
      call addto (s_field, QC, scale=-1.)
    endif
    if (has_scalar_field(state,'Qrain')) then
      QR=>extract_scalar_field(state, "Qrain", stat)
      call addto (s_field, QR, scale=-1.)
    endif
    if (has_scalar_field(state,'Qice')) then
      QI=>extract_scalar_field(state, "Qice", stat)
      call addto (s_field, QI, scale=-1.)
    endif
    if (has_scalar_field(state,'Qgrau')) then
      QG=>extract_scalar_field(state, "Qgrau", stat)
      call addto (s_field, QG, scale=-1.)
    endif
    if (has_scalar_field(state,'Qsnow')) then
      QS=>extract_scalar_field(state, "Qsnow", stat)
      call addto (s_field, QS, scale=-1.)
    endif
    
    endif

  end subroutine get_water_vapour
  
  subroutine get_reflectivity(state, s_field, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    integer, intent(out), optional :: stat 
    
    type(scalar_field), pointer :: nr,ni,ng,ns,lambda
    type(scalar_field), target :: dummyscalar
    type(scalar_field) :: temperature,n_remap
    real :: param, kappa, ze, eps=1.e-3
    character(len=OPTION_PATH_LEN) :: mom_path
    integer :: ele
    
    call zero(s_field)
      
    if (have_option('/material_phase[0]/cloud_microphysics')) then
    
      call allocate(dummyscalar,s_field%mesh,"DummyScalar")
      call allocate(temperature,s_field%mesh,"TemperatureLocal")
      call allocate(n_remap,s_field%mesh,"NRemap")
      
      call compressible_eos(state,temperature=temperature)
      call zero(dummyscalar)
      
      ! Rain
      if (has_scalar_field(state,"Nrain")) then
        param = (nu(2)+6.)*(nu(2)+5.)*(nu(2)+4.)*(nu(2)+3.)*(nu(2)+2.)*(nu(2)+1.)
        nr=>extract_scalar_field(state, "Nrain", stat)
	lambda=>dummyscalar
	
	call safe_set(state,n_remap,nr)
        call calculate_lambda_two(state,"rain",lambda)
	
	do ele = 1, node_count(lambda)
	  ze=node_val(lambda,ele)
	  if (node_val(n_remap,ele) > eps .and. ze > eps) then
  	    ze=param*node_val(n_remap,ele)*ze**(-6.)
	    call addto(s_field,ele,ze)
	  endif
	enddo
      end if
      
      ! Ice
      if (has_scalar_field(state,"Nice")) then
        param = 0.224*(nu(3)+6.)*(nu(3)+5.)*(nu(3)+4.)*(nu(3)+3.)*(nu(3)+2.)*(nu(3)+1.)
        ni=>extract_scalar_field(state, "Nice", stat)
	lambda=>dummyscalar

	call safe_set(state,n_remap,ni)
        call calculate_lambda_two(state,"ice",lambda)
	
	do ele = 1, node_count(lambda)
	  ze=node_val(lambda,ele)
	  if (node_val(n_remap,ele) > eps .and. ze > eps .and. node_val(temperature,ele) < 273.15) then
  	    ze=0.224*param*node_val(n_remap,ele)*ze**(-6.)
	    call addto(s_field,ele,ze)
	  endif
	enddo
      end if
      
      ! Graupel
      if (has_scalar_field(state,"Ngrau")) then
        param = (nu(4)+6.)*(nu(4)+5.)*(nu(4)+4.)*(nu(4)+3.)*(nu(4)+2.)*(nu(4)+1.)
        ng=>extract_scalar_field(state, "Ngrau", stat)
	lambda=>dummyscalar

	call safe_set(state,n_remap,ng)
        call calculate_lambda_two(state,"grau",lambda)
	
	do ele = 1, node_count(lambda)
	  ze=node_val(lambda,ele)
	  if (node_val(n_remap,ele) > eps .and. ze > eps) then
	    kappa=1.
	    if (node_val(temperature,ele) < 273.15) kappa=0.224
  	    ze=kappa*param*node_val(n_remap,ele)*ze**(-6.)
	    call addto(s_field,ele,ze)
	  endif
	enddo
      end if
      
      ! Snow
      if (has_scalar_field(state,"Nsnow")) then
        param =(nu(5)+6.)*(nu(5)+5.)*(nu(5)+4.)*(nu(5)+3.)*(nu(5)+2.)*(nu(5)+1.)
        ns=>extract_scalar_field(state, "Nsnow", stat)
	lambda=>dummyscalar

	call safe_set(state,n_remap,ns)
        call calculate_lambda_two(state,"snow",lambda)
	
	do ele = 1, node_count(lambda)
	  ze=node_val(lambda,ele)
	  if (node_val(n_remap,ele) > eps .and. ze > eps) then
	    kappa=1.
	    if (node_val(temperature,ele) < 273.15) kappa=0.224
  	    ze=kappa*param*node_val(n_remap,ele)*ze**(-6.)
	    call addto(s_field,ele,ze)
	  endif
	enddo
      end if
      
      call deallocate (dummyscalar)
      call deallocate (temperature)
      call deallocate (n_remap) 
    endif

  end subroutine get_reflectivity
          
end module diagnostic_fields_eos
