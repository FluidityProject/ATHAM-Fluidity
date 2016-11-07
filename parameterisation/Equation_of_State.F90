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

module equation_of_state
  !!< This module contains functions used to evaluate the equation of state.

  use fldebug
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use futils, only: int2str
  use spud
  use fields
  use state_module
  use diagnostic_fields, only: safe_set
  use sediment, only: get_n_sediment_fields, get_sediment_item
  
  implicit none

  real :: c_p,c_v,c_v_v,c_p_v,c_v_l,c_p_l,c_v_i,c_p_i
  
  interface compressible_eos
      module procedure compressible_eos_1mat, compressible_eos_mmat
  end interface compressible_eos
  
  private
  public :: calculate_perturbation_density, mcD_J_W_F2002, &
            compressible_eos, compressible_material_eos, &
            initialise_from_eos, set_EOS_pressure_and_temperature, &
            make_giraldo_quantities_1mat, extract_entropy_variable, scale_pressure, &
	    c_p,c_v,c_p_v,c_v_v,c_v_l,c_p_l, &
	    get_thermo_variable, get_cp_cv, cal_flv, cal_fls, cal_flm, &
	    cal_esat, cal_esati, cal_qsat, cal_qsati, cal_dqsatdt, cal_dqsatidt

contains

  subroutine calculate_perturbation_density(state, density, reference_density)
    !!< Calculates the perturbation density (i.e. the reference density is already subtracted)
    !!< of a state with equation_of_state fluids/linear or 
    !!< fluids/ocean_pade_approximation.
    type(state_type), intent(in):: state
    type(scalar_field), intent(inout) :: density
    real, intent(out), optional :: reference_density
    
    type(vector_field), pointer:: u
    type(scalar_field), pointer:: T, S, oldT, oldS, topdis
    type(scalar_field) :: DeltaT, DeltaS, remapT, remapS, fluidconcentration,&
         & sedimentdensity
    character(len=OPTION_PATH_LEN) option_path, dep_option_path, sediment_field_name, class_name, sfield_name
    logical, dimension(:), allocatable:: done
    logical include_depth_below
    real T0, S0, gamma, rho_0, salt, temp, dist, dens, theta
    integer, dimension(:), pointer :: density_nodes
    integer ele, i, node, n_sediment_fields, f
    
    ewrite(1,*) 'In calculate_perturbation_density'
    
    u => extract_vector_field(state, "Velocity")
    
    call zero(density)
    
    option_path='/material_phase::'//trim(state%name)//'/equation_of_state/fluids'
    
    call get_option(trim(u%option_path)//'/prognostic/temporal_discretisation/relaxation', &
                    theta, default = 1.0)
    
    rho_0 = 0.0
    
    if (have_option(trim(option_path)//'/linear')) then
    
       option_path=trim(option_path)//'/linear'
       
       if (have_option(trim(option_path)//'/temperature_dependency')) then
          dep_option_path=trim(option_path)//'/temperature_dependency'
          call get_option(trim(dep_option_path)//'/reference_temperature', T0)
          call get_option(trim(dep_option_path)//'/thermal_expansion_coefficient', gamma)
          T => extract_scalar_field(state, "Temperature")
          oldT => extract_scalar_field(state, "OldTemperature")
          call allocate(deltaT, density%mesh, "DeltaT")
          call allocate(remapT, density%mesh, "RemapT")
          
          ! deltaT=theta*T+(1-theta)*oldT-T0
          call remap_field(T, remapT)
          call set(deltaT, remapT)
          call scale(deltaT, theta)
          
          call remap_field(oldT, remapT)
          call addto(deltaT, remapT, 1.0-theta)
          call addto(deltaT, -T0)
          ! density=density-gamma*deltaT
          call addto(density, deltaT, scale=-gamma)
          call deallocate(deltaT)
          call deallocate(remapT)
       end if
       
       if (have_option(trim(option_path)//'/salinity_dependency')) then
          dep_option_path=trim(option_path)//'/salinity_dependency'
          call get_option(trim(dep_option_path)//'/reference_salinity', S0)
          call get_option(trim(dep_option_path)//'/saline_contraction_coefficient', gamma)
          S => extract_scalar_field(state, "Salinity")
          oldS => extract_scalar_field(state, "OldSalinity")
          call allocate(deltaS, density%mesh, "DeltaS")
          call allocate(remapS, density%mesh, "RemapS")
          
          ! deltaS=theta*S+(1-theta)*oldS-S0
          call remap_field(S, remapS)
          call set(deltaS, remapS)
          call scale(deltaS, theta)
          
          call remap_field(oldS, remapS)
          call addto(deltaS, remapS, 1.0-theta)
          call addto(deltaS, -S0)
          ! density=density+gamma*deltaS
          call addto(density, deltaS, scale=gamma)
          call deallocate(deltaS)
          call deallocate(remapS)
       end if

       if (have_option(trim(option_path)//'/generic_scalar_field_dependency')) then
	  do f = 1, option_count(trim(option_path)//'/generic_scalar_field_dependency')
	     dep_option_path=trim(option_path)//'/generic_scalar_field_dependency['//int2str(f-1)//']'
	     call get_option(trim(dep_option_path)//'/name', sfield_name)
	     call get_option(trim(dep_option_path)//'/reference_value', T0)
	     call get_option(trim(dep_option_path)//'/expansion_coefficient', gamma)
	     T => extract_scalar_field(state, trim(sfield_name))
	     oldT => extract_scalar_field(state, "Old"//trim(sfield_name))
	     call allocate(deltaT, density%mesh, "DeltaT")
	     call allocate(remapT, density%mesh, "RemapT")

	     ! deltaT=theta*T+(1-theta)*oldT-T0
	     call remap_field(T, remapT)
	     call set(deltaT, remapT)
	     call scale(deltaT, theta)

	     call remap_field(oldT, remapT)
	     call addto(deltaT, remapT, 1.0-theta)
	     call addto(deltaT, -T0)
	     ! density=density-gamma*deltaT
	     call addto(density, deltaT, scale=-gamma)
	     call deallocate(deltaT)
	     call deallocate(remapT)
	  end do
       end if
       
       call get_option(trim(option_path)//'/reference_density', rho_0)
       call scale(density, rho_0)
       
    elseif (have_option(trim(option_path)//'/ocean_pade_approximation')) then
      
       option_path=trim(option_path)//'/ocean_pade_approximation'
       
       include_depth_below=have_option(trim(option_path)//'/include_depth_below_surface')
       
       T => extract_scalar_field(state, "Temperature")
       oldT => extract_scalar_field(state, "OldTemperature")
       S => extract_scalar_field(state, "Salinity")
       oldS => extract_scalar_field(state, "OldSalinity")
       if (include_depth_below) then
          topdis => extract_scalar_field(state, "DistanceToTop")
       endif
      
       allocate( done(1:node_count(density)) )
       done=.false.
       
       do ele=1, element_count(density)

          density_nodes => ele_nodes(density, ele)

          do i=1,size(density_nodes)
             node=density_nodes(i)
             ! In the continuous case ensure we only do each calculation once.
             if (done(node)) cycle
             done(node)=.true.
            
             salt=theta*node_val(S, node)+(1-theta)*node_val(oldS, node)
             temp=node_val(T, node)+(1-theta)*node_val(oldT, node)
             if (include_depth_below) then
                dist=node_val(topdis, node)
             else
                dist=0.0
             end if            
               
             call mcD_J_W_F2002(dens,temp,salt,dist)
             call addto(density, node, dens)
          end do
           
       end do
         
       ! reference density is assumed 1 for the pade approximation
       rho_0=1.0

    end if
    
    if (have_option('/material_phase::'//trim(state%name)//'/sediment')) then 

       call allocate(deltaS, density%mesh, "DeltaS")
       call allocate(remapS, density%mesh, "RemapS")
       call allocate(sedimentdensity, density%mesh, "SedimentDensity")
       call zero(sedimentdensity)

       n_sediment_fields = get_n_sediment_fields()

       do i=1,n_sediment_fields
       
	  call get_sediment_item(state, i, S)
	  call get_sediment_item(state, i, 'submerged_specific_gravity', gamma)
	  gamma = gamma * rho_0

          oldS => extract_scalar_field(state, &
               "Old"//trim(S%name))
          
          ! deltaS=theta*S+(1-theta)*oldS-S0
          call remap_field(S, remapS)
          call set(deltaS, remapS)
          call scale(deltaS, theta)
          
          call remap_field(oldS, remapS)
          call addto(deltaS, remapS, 1.0-theta)
          ! density=density+gamma*deltaS
          call addto(sedimentdensity, deltaS, scale=gamma)
       end do
       
       call addto(density,sedimentdensity)

       call deallocate(deltaS)
       call deallocate(remapS)
       call deallocate(sedimentdensity)
       
    end if

    if(present(reference_density)) then
      reference_density = rho_0
    end if
    
  end subroutine calculate_perturbation_density

  subroutine mcD_J_W_F2002(density,T,Salinity,distance_to_top)
    !!<  function to evaluate density from the 2002 McDougall, Jackett,
    !!< Wright and Feistel equation of state using Pade approximation.  
    real, intent(out) :: density
    real, intent(in) :: T,Salinity,distance_to_top
   
    real :: p,p1,p2,S

    ! Salinity can be negitive because it's numerically diffused,
    ! some regions may be initialised with zero salinity, and
    ! therefore undershoot may occur.
    S = max(Salinity, 0.0)

    ! calculate pressure in decibars from hydrostatic pressure
    ! using reference density 1000 kg m^-2

    p = 9.81*1000.0*distance_to_top*1.0e-4

    !     evaluate top and bottom of Pade approximant
      
    p1 = 9.99843699e2 &
         + 7.35212840*T - 5.45928211e-2*(T**2) + 3.98476704e-4*(T**3) &
         + 2.96938239*S - 7.23268813e-3*S*T + 2.12382341e-3*(S**2) &
         + 1.04004591e-2*p + 1.03970529e-7*p*(T**2) &
         + 5.18761880e-6*p*S - 3.24041825e-8*(p**2) &
         - 1.23869360e-11*(p**2)*(t**2)

    p2 = 1.0 &
         + 7.28606739e-3*T - 4.60835542e-5*(T**2) + 3.68390573e-7*(T**3) &
         + 1.80809186e-10*(T**4) &
         + 2.14691708e-3*S - 9.27062484e-6*S*T - 1.78343643e-10*S*(T**3) &
         + 4.76534122e-6*(S**1.5) + 1.63410736e-9*(S**1.5)*(T**2) &
         + 5.30848875e-6*p -3.03175128e-16*(p**2)*(t**3) &
         - 1.27934137e-17*(p**3)*T
    
    ! calculate the resulting density
    
    density = p1/p2
    
    ! the perturbation density
    
    density = (density-1000.0)/1000.0
    
  end subroutine mcD_J_W_F2002  
  
  subroutine compressible_eos_1mat(state,full_pressure,full_density,energy,density,pressure,	&
  				   drhodp,saturation,supersaturation,dqsaturation,temperature,	&
				   potentialtem,density_pottem,qc_p,qc_v,sound_speed, getold)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), target, optional :: energy, full_pressure, full_density
    type(scalar_field), intent(inout), optional :: density, drhodp, saturation, supersaturation, dqsaturation
    type(scalar_field), intent(inout), optional :: pressure, temperature, potentialtem, density_pottem, qc_p, qc_v,sound_speed
    logical, intent(in), optional :: getold

    type(scalar_field), pointer :: output, local_density, local_pressure, thermal, lp
    type(scalar_field), pointer :: q_v,q_c,q_r,q_i,q_g,q_s
    type(scalar_field), target  :: dummyscalar
    
    character(len=OPTION_PATH_LEN) :: eos_path
    type(state_type), dimension(1) :: states
    type(scalar_field) :: drhodp_local
    logical :: getoldlocal=.false., have_vapour=.true., have_liquid=.true., have_ice=.true.
    integer :: thermal_variable, stat

    if (present(getold)) then
       getoldlocal=getold
    else
       getoldlocal=.false.
    end if
    
    ewrite(1,*) 'Entering compressible_eos'
    
    if (present(drhodp)) then
      drhodp_local=drhodp
      if (present(density)) then
         assert(drhodp%mesh==density%mesh)
      end if
      if (present(pressure)) then
         assert(drhodp%mesh==pressure%mesh)
      end if
      if (present(full_pressure)) then
         assert(drhodp%mesh==full_pressure%mesh)
      end if
    else if (present(density)) then
      call allocate(drhodp_local, density%mesh, 'Localdrhop')
    else if (present(pressure)) then
      call allocate(drhodp_local, pressure%mesh, 'Localdrhop')
    else if (present(temperature)) then
      call allocate(drhodp_local, temperature%mesh, 'Localdrhop')
    else if (present(potentialtem)) then
      call allocate(drhodp_local, potentialtem%mesh, 'Localdrhop')
    else if (present(density_pottem)) then
      call allocate(drhodp_local, density_pottem%mesh, 'Localdrhop')
    else if (present(saturation)) then
      call allocate(drhodp_local, saturation%mesh, 'Localdrhop')
    else if (present(sound_speed)) then
      call allocate(drhodp_local, sound_speed%mesh, 'Localdrhodp')
    else
      FLAbort("No point in being in here if you don't want anything out.")
    end if
    
   call get_cp_cv (c_p,c_v,c_p_v=c_p_v,c_v_v=c_v_v,c_p_l=c_p_l,c_v_l=c_v_l,c_p_i=c_p_i,c_v_i=c_v_i)
   
   eos_path = '/material_phase::'//trim(state%name)//'/equation_of_state'
   
   if(have_option(trim(eos_path)//'/compressible')) then
	 
      ! each of the following compressible_eos_XXX() routines should always calculate drhodp
      ! (zero if density does not depend on pressure) and calculate density and
      ! pressure if present
      
      if(have_option(trim(eos_path)//'/compressible/stiffened_gas')) then
         
         ! standard stiffened gas eos
         
         if (present(pressure))then
            call compressible_eos_stiffened_gas(state, eos_path,&
                 drhodp_local, density=density, pressure=pressure)
         else 
            call compressible_eos_stiffened_gas(state, eos_path,&
                 drhodp_local, density=density, pressure=full_pressure)
         end if

      elseif(have_option(trim(eos_path)//'/compressible/foam')) then
        
        ! eos used in foam modelling
        
         call compressible_eos_foam(state, eos_path, drhodp_local, &
              density=density, pressure=pressure)
         
      else if(have_option(trim(eos_path)//'/compressible/giraldo')) then
         
        ! Eq. of state commonly used in atmospheric applications. See
        ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
        ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
        
	 if (present(full_pressure)) then
	   local_pressure => full_pressure
	 else
           local_pressure => extract_scalar_field(state,"Pressure",stat=stat)
         endif
        
	 if (present(full_density)) then
	   local_density => full_density
	 else
           local_density => extract_scalar_field(state,"Density",stat=stat)
         endif
        
	 if (present(energy)) then
	   thermal => energy
	   if (trim(energy%name)=='PotentialTemperature' .or. &
	       trim(energy%name)=='LocalPotentialTemperature' .or. &
	       trim(energy%name)=='HydrostaticReferencePotentialTemperature') then
	       thermal_variable=3
	   else if (trim(energy%name)=='Temperature') then
	       thermal_variable=2
	   else if (trim(energy%name)=='ConservedPotentialTemperature') then
	       thermal_variable=5
	   endif
	 else
           call get_thermo_variable(state,thermal,index=thermal_variable)
         endif
	 
         call allocate (dummyscalar,thermal%mesh,"DummyScalar")
	 call zero(dummyscalar)
    
         if (has_scalar_field(state,"TotalWaterQ")) then
            q_v=>extract_scalar_field(state,"TotalWaterQ")
         else if (has_scalar_field(state,"VapourWaterQ")) then
            q_v=>extract_scalar_field(state,"VapourWaterQ")
         else
	    have_vapour=.false.
            q_v=>dummyscalar
         end if
         if (has_scalar_field(state,"Qdrop")) then
            q_c=>extract_scalar_field(state,"Qdrop")
         else
	    have_liquid=.false.
            q_c=>dummyscalar
         end if
         if (has_scalar_field(state,"Qrain")) then
            q_r=>extract_scalar_field(state,"Qrain")
         else
	    have_liquid=.false.
            q_r=>dummyscalar
         end if
         if (has_scalar_field(state,"Qice")) then
            q_i=>extract_scalar_field(state,"Qice")
         else
	    have_ice=.false.
            q_i=>dummyscalar
         end if         
	 if (has_scalar_field(state,"Qgrau")) then
            q_g=>extract_scalar_field(state,"Qgrau")
         else
	    have_ice=.false.
            q_g=>dummyscalar
         end if         
	 if (has_scalar_field(state,"Qsnow")) then
            q_s=>extract_scalar_field(state,"Qsnow")
         else
	    have_ice=.false.
            q_s=>dummyscalar
         end if  
	 
         if (present(pressure))then
           ewrite(1,*) 'compressible_eos_giraldo: compute pressure'
           call compressible_eos_giraldo_1mat(state,eos_path,q_v,q_c,q_r,q_i,q_g,q_s,   &
                local_pressure,local_density,thermal,thermal_variable,drhodp_local,	&
		have_vapour, have_liquid, have_ice, pressure=pressure)
         endif

         if (present(density))then
           ewrite(1,*) 'compressible_eos_giraldo: compute density'
           call compressible_eos_giraldo_1mat(state,eos_path,q_v,q_c,q_r,q_i,q_g,q_s,   &
                local_pressure,local_density,thermal,thermal_variable,drhodp_local,	&
		have_vapour, have_liquid, have_ice, density=density)
         endif
	 
	 if (present(temperature).or.present(potentialtem))then
           ewrite(1,*) 'compressible_eos_giraldo: compute temperature'
           call compressible_eos_giraldo_1mat(state,eos_path,q_v,q_c,q_r,q_i,q_g,q_s,   &
                local_pressure,local_density,thermal,thermal_variable,drhodp_local,	&
		have_vapour, have_liquid, have_ice, temperature=temperature, potentialtem=potentialtem)
         end if
	 
	 if (present(density_pottem))then
           ewrite(1,*) 'compressible_eos_giraldo: compute density potential temperature'
           call compressible_eos_giraldo_1mat(state,eos_path,q_v,q_c,q_r,q_i,q_g,q_s,   &
                local_pressure,local_density,thermal,thermal_variable,drhodp_local,	&
		have_vapour, have_liquid, have_ice, density_pottem=density_pottem)
         end if

	 if (present(saturation)) then
           ewrite(1,*) 'compressible_eos_giraldo: compute saturation'
	   call compressible_eos_giraldo_1mat (state,eos_path,q_v,q_c,q_r,q_i,q_g,q_s,  &
	        local_pressure,local_density,thermal,thermal_variable,drhodp_local,	&
		have_vapour, have_liquid, have_ice, saturation=saturation,dqsaturation=dqsaturation)
		
	   if (present(supersaturation)) then
	     call set(supersaturation,q_v)
	     if (has_scalar_field(state,"TotalWaterQ")) then
	       call addto(supersaturation,q_c,scale=-1.)
	       call addto(supersaturation,q_r,scale=-1.)
	       call addto(supersaturation,q_i,scale=-1.)
	       call addto(supersaturation,q_g,scale=-1.)
	       call addto(supersaturation,q_s,scale=-1.)
	     endif
	     call addto(supersaturation,saturation,-1.)	     
	   endif
	 endif

	 if (present(qc_p).or.present(qc_v)) then
           ewrite(1,*) 'compressible_eos_giraldo: compute cp cv'
	   call compressible_eos_giraldo_1mat (state,eos_path,q_v,q_c,q_r,q_i,q_g,q_s,  &
	        local_pressure,local_density,thermal,thermal_variable,drhodp_local,	&
		have_vapour, have_liquid, have_ice, qc_p=qc_p,qc_v=qc_v)
	 endif

	 if (present(sound_speed)) then
           ewrite(1,*) 'compressible_eos_giraldo: compute sound speed'
	   call compressible_eos_giraldo_1mat (state,eos_path,q_v,q_c,q_r,q_i,q_g,q_s,  &
	        local_pressure,local_density,thermal,thermal_variable,drhodp_local,	&
		have_vapour, have_liquid, have_ice, sound_speed=sound_speed)
	 endif
	 	 
	 call deallocate (dummyscalar)	 

      else if(have_option(trim(eos_path)//'/compressible/ATHAM')) then
 
         call get_thermo_variable(state,thermal,index=thermal_variable)

         if (present(pressure))then
           output=>extract_scalar_field(state,"Pressure")
	 else if (present(density)) then
           output=>extract_scalar_field(state,"Density")
	 endif
	 
	 states=(/state/)
         if (present(pressure))then
            call compressible_eos_atham(states,&
	        pressure,"Pressure",thermal,thermal_variable,output2=drhodp_local,getold=getold)
         else if (present(density))then
            call compressible_eos_atham(states,&
	        density,"Density",thermal,thermal_variable,output2=drhodp_local,getold=getold)
         endif
	 
	 if (present(temperature).or.present(potentialtem))then
            call compressible_eos_atham(states,&
	 	temperature,"Temperature",thermal,thermal_variable,getold=getold)
            call compressible_eos_atham(states,&
	 	potentialtem,"PotentialTemperature",thermal,thermal_variable,getold=getold)
         end if

	 if (present(saturation)) then
	   call compressible_eos_atham(states,&
	 	saturation,"Saturation",thermal,thermal_variable,output2=dqsaturation,qc_p=qc_p,qc_v=qc_v,getold=getold)
	 endif

	 if (present(qc_p).or.present(qc_v)) then
	   call compressible_eos_atham(states,&
	 	density,"Density",thermal,thermal_variable,qc_p=qc_p,qc_v=qc_v,getold=getold)
	 endif
	 
	 state=states(1)
      else
      FLAbort('Gone into compressible_eos without having equation_of_state/compressible')
   end if
   
   end if

    if(present(density)) then
      ewrite_minmax(density)
    end if

    if(present(pressure)) then
      ewrite_minmax(pressure)
    end if

    if(present(temperature)) then
      ewrite_minmax(temperature)
    end if

    if(present(potentialtem)) then
      ewrite_minmax(potentialtem)
    end if  

    if(present(density_pottem)) then
      ewrite_minmax(density_pottem)
    end if  

    if(present(saturation)) then
      ewrite_minmax(saturation)
    end if  

    if(present(supersaturation)) then
      ewrite_minmax(supersaturation)
    end if  
    
    if(present(full_pressure)) then
      ewrite_minmax(full_pressure%val)
    end if

    if(present(drhodp)) then      
      ewrite_minmax(drhodp)
    else
      call deallocate(drhodp_local)
    end if

    ewrite(1,*) 'Exiting compressible_eos'

  end subroutine compressible_eos_1mat

  subroutine compressible_eos_mmat(state,full_pressure,full_density,energy,density,pressure,	& 
       		drhodp,saturation,dqsaturation,temperature,potentialtem,density_pottem,qc_p,qc_v,sound_speed,getold)

    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), intent(inout), optional :: full_pressure,full_density,energy
    type(scalar_field), intent(inout), optional :: density, pressure, drhodp
    type(scalar_field), intent(inout), optional :: saturation, dqsaturation, qc_p, qc_v, sound_speed
    type(scalar_field), intent(inout), optional :: temperature, potentialtem, density_pottem
    logical, intent(in), optional :: getold

    integer :: stat, thermal_variable
    
    type(scalar_field), pointer :: output,thermal
    character(len=OPTION_PATH_LEN) :: eos_path
    logical :: getoldlocal = .false.
    
    if (present(getold)) then
       getoldlocal=getold
    else
       getoldlocal=.false.
    end if

    ewrite(1,*) 'Entering compressible_eos_mmat'

    if (size(state)==1) then !<--------
       call compressible_eos_1mat(state(1), full_pressure=full_pressure, full_density=full_density, energy=energy, &
            density=density, pressure=pressure,drhodp=drhodp,saturation=saturation,dqsaturation=dqsaturation,&
	    temperature=temperature,potentialtem=potentialtem,density_pottem=density_pottem,qc_p=qc_p,qc_v=qc_v,getold=getold)
       return
    else
    
    if (present(drhodp)) then
       if (present(density)) then
          assert(drhodp%mesh==density%mesh)
       end if
       if (present(pressure)) then
          assert(drhodp%mesh==pressure%mesh)
       end if
    end if
    
    if (.not. (present(drhodp) .or. present(density) .or.&
         present(pressure) .or. present(temperature))) then
       FLAbort("No point in being in here if you don't want anything out.")
    end if
   
    eos_path = trim(state(1)%option_path)//'/equation_of_state'
    
    if(have_option(trim(eos_path)//'/compressible')) then       
       if(have_option(trim(eos_path)//'/compressible/ATHAM')) then

	 call get_thermo_variable(state(1),thermal,index=thermal_variable)
         if (present(pressure))then
          output=>extract_scalar_field(state(1),"Pressure")
	 else if (present(density)) then
          output=>extract_scalar_field(state(1),"Density")
	 endif
 	         
         !ATHAM equation of state
         if (present(pressure)) then
            if (present(drhodp)) then
               call compressible_eos_atham(state,&
	           pressure,"Pressure",thermal,thermal_variable,output2=drhodp,getold=getold)
            else
               call compressible_eos_atham(state,&
	           pressure,"Pressure",thermal,thermal_variable,getold=getold)
            end if
         else if (present(density)) then
            if (present(drhodp)) then
               call compressible_eos_atham(state,&
	           density,"Density",thermal,thermal_variable,output2=drhodp,getold=getold)
            else
               call compressible_eos_atham(state,&
	           density,"Density",thermal,thermal_variable,getold=getold)
            end if
         else if (present(temperature)) then
            call compressible_eos_atham(state,&
	 	temperature,"Temperature",thermal,thermal_variable,getold=getold)
            call compressible_eos_atham(state,&
	 	potentialtem,"PotentialTemperature",thermal,thermal_variable,getold=getold)
         end if
	 
         if (present(saturation)) then
            call compressible_eos_atham(state,saturation,"Saturation", &
	        thermal,thermal_variable,output2=dqsaturation,qc_p=qc_p,qc_v=qc_v,getold=getold)
         end if
	 
         if (present(qc_p).or.present(qc_v)) then
            call compressible_eos_atham(state,density,"Density", &
	         thermal,thermal_variable,qc_p=qc_p,qc_v=qc_v,getold=getold)
         end if
       else

          FLAbort('Gone into multimaterial compressible_eos without having equation_of_state/compressible/ATHAM')

       end if
                 
    else
       FLAbort('Gone into compressible_eos without having equation_of_state/compressible')
    end if
    
    endif !<--------

    if(present(density)) then
       ewrite_minmax(density)
    end if

    if(present(pressure)) then
       ewrite_minmax(pressure)
    end if

    if(present(drhodp)) then      
      ewrite_minmax(drhodp)
   end if

   if(present(temperature)) then      
      ewrite_minmax(temperature)
   end if
   
   if(present(saturation)) then      
      ewrite_minmax(saturation)
   end if
   
 end subroutine compressible_eos_mmat
    
  subroutine compressible_eos_stiffened_gas(state, eos_path, drhodp, &
    density, pressure)
    ! Standard stiffened gas equation
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure
    
    !locals
    integer :: stat, gstat, cstat
    type(scalar_field), pointer :: pressure_local, energy_local, density_local, hp
    real :: reference_density, ratio_specific_heats
    real :: bulk_sound_speed_squared, atmospheric_pressure
    type(scalar_field) :: energy_remap, pressure_remap, density_remap
    logical :: incompressible

    character(len = *), parameter:: hp_name = "HydrostaticReferencePressure"
    character(len = *), parameter:: hpg_name= "HydrostaticPressureGradient"
    
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/reference_density', &
                        reference_density, default=0.0)
        
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/ratio_specific_heats', &
                    ratio_specific_heats, stat=gstat)
    if(gstat/=0) then
      ratio_specific_heats=1.4
    end if
    
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/bulk_sound_speed_squared', &
                    bulk_sound_speed_squared, stat=cstat)
    if(cstat/=0) then
      bulk_sound_speed_squared=0.0
    end if
    
    incompressible = ((gstat/=0).and.(cstat/=0))
    if(incompressible) then
      ewrite(0,*) "Selected compressible eos but not specified a bulk_sound_speed_squared or a ratio_specific_heats."
    end if
    
    call zero(drhodp)
    
    if(.not.incompressible) then
      energy_local=>extract_scalar_field(state,'InternalEnergy',stat=stat)
      ! drhodp = 1.0/( bulk_sound_speed_squared + (ratio_specific_heats - 1.0)*energy )
      if((stat==0).and.(gstat==0)) then   ! we have an internal energy field and we want to use it
        call allocate(energy_remap, drhodp%mesh, 'RemappedInternalEnergy')
        call safe_set(state,energy_remap,energy_local)
        
        call addto(drhodp, energy_remap, (ratio_specific_heats-1.0))
        
        call deallocate(energy_remap)
      end if
      call addto(drhodp, bulk_sound_speed_squared)
      call invert(drhodp)
    end if

    if(present(density)) then
      ! calculate the density
      ! density may equal density in state depending on how this
      ! subroutine is called
      if(incompressible) then
        ! density = reference_density
        call set(density, reference_density)
      else
        pressure_local=>extract_scalar_field(state,'Pressure',stat=stat)
        if (has_scalar_field(state,hp_name)) then
           hp => extract_scalar_field(state,hp_name)
           call addto(pressure_local,hp)
        end if
        if (stat==0) then
          assert(density%mesh==drhodp%mesh)
        
          ! density = drhodp*(pressure_local + atmospheric_pressure
          !                  + bulk_sound_speed_squared*reference_density)
          call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                          atmospheric_pressure, default=0.0)
          
          call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
          call safe_set(state,pressure_remap,pressure_local)
          
          call set(density, reference_density*bulk_sound_speed_squared + atmospheric_pressure)
          call addto(density, pressure_remap)
          call scale(density, drhodp)
          
          call deallocate(pressure_remap)
        else
          FLExit('No Pressure in material_phase::'//trim(state%name))
        end if
      end if
    end if

    if(present(pressure)) then
      if(incompressible) then
        ! pressure is unrelated to density in this case
        call zero(pressure)
      else
        ! calculate the pressure using the eos and the calculated (probably prognostic)
        ! density
        density_local=>extract_scalar_field(state,'Density',stat=stat)
        if (stat==0) then
          assert(pressure%mesh==drhodp%mesh)
          
          ! pressure = density_local/drhodp &
          !          - bulk_sound_speed_squared*reference_density
          
          call allocate(density_remap, drhodp%mesh, "RemappedDensity")
          call safe_set(state,density_remap,density_local)
          
          call set(pressure, drhodp)
          call invert(pressure)
          call scale(pressure, density_remap)
          call addto(pressure, -bulk_sound_speed_squared*reference_density)
          
          call deallocate(density_remap)
        else
          FLExit('No Density in material_phase::'//trim(state%name))
        end if
      end if
    end if

  end subroutine compressible_eos_stiffened_gas
  
  subroutine compressible_eos_giraldo_1mat(state, eos_path, q_v, q_c, q_r, q_i, q_g, q_s,	&
    pressure_local,density_local,thermal_local, thermal_variable, drhodp, &
    have_vapour, have_liquid, have_ice, density, pressure, temperature, potentialtem, &
    density_pottem, saturation, dqsaturation, qc_p, qc_v, sound_speed)
    ! Eq. of state commonly used in atmospheric applications. See
    ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
    ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
    
    integer, intent(in) :: thermal_variable
    logical, intent(in) :: have_vapour, have_ice, have_liquid
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp,q_v,q_c,q_r,q_i,q_g,q_s
    type(scalar_field), intent(inout) :: pressure_local,density_local,thermal_local
    type(scalar_field), intent(inout), optional :: density, pressure, temperature, potentialtem, density_pottem
    type(scalar_field), intent(inout), optional :: saturation, dqsaturation, qc_p, qc_v, sound_speed
    
    ! locals
    type(scalar_field) :: energy_remap, pressure_remap, density_remap, &
    	 drhodp_remap, qg, temperature_remap, potentialtem_remap, &
	 density_pottem_remap,saturation_remap,dqsaturation_remap, &
	 incompfix,incompfix2,qc_p_remap,qc_v_remap,sound_speed_remap, &
	 qv_remap,qc_remap,qr_remap,qi_remap,qg_remap,qs_remap
    real :: reference_density, p_0, r_cp, r, gamm
    real :: drhodp_node, rhog_node, temperature_node, pottem_node, dpt_node, density_node, c_node, &
            pressure_node, energy_node,esat_node,saturation_node,dqsaturation_node,i_node,i2_node, qt_node
    real :: exn, epsilon, rho_w=1000., rho_i=920.
    integer :: node
    logical :: constant_cp_cv
        
    call get_option(trim(eos_path)//'/compressible/giraldo/reference_pressure',p_0, default=1.0e5)
    constant_cp_cv = have_option(trim(eos_path)//'/compressible/giraldo/constant_cp_cv')
    
    call allocate(incompfix,drhodp%mesh,"IncompressibleScaleFactor")          
    call allocate(incompfix2,drhodp%mesh,"IncompressibleScaleFactor2")          
    call allocate(energy_remap, drhodp%mesh, "RemappedThermal")
    call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
    call allocate(density_remap, drhodp%mesh, "RemappedDensity")
    call allocate(qc_p_remap, drhodp%mesh, "RemappedCP")
    call allocate(qc_v_remap, drhodp%mesh, "RemappedCV")
    call allocate(sound_speed_remap, drhodp%mesh, "RemappedSoundSpeed")
    call allocate(potentialtem_remap, drhodp%mesh, "RemappedPotentialTemperature")
    call allocate(temperature_remap, drhodp%mesh, "RemappedTemperature")
    
    call safe_set(state,density_remap,density_local)   
    call safe_set(state,pressure_remap,pressure_local)   
    call safe_set(state,energy_remap,thermal_local)    
    
    call zero(incompfix)
    call zero(incompfix2)
    call zero(drhodp) 
    call zero(potentialtem_remap)
    call zero(temperature_remap)
    
    if (have_vapour) then     
      call allocate(qv_remap, drhodp%mesh, "RemappedQV")
      call safe_set(state,qv_remap,q_v)	
    endif
    if (have_liquid) then
      call allocate(qc_remap, drhodp%mesh, "RemappedQC")
      call safe_set(state,qc_remap,q_c)	
      call allocate(qr_remap, drhodp%mesh, "RemappedQR")
      call safe_set(state,qr_remap,q_r)	
    endif
    if (have_ice) then
      call allocate(qi_remap, drhodp%mesh, "RemappedQI")
      call safe_set(state,qi_remap,q_i)	
      call allocate(qg_remap, drhodp%mesh, "RemappedQG")
      call safe_set(state,qg_remap,q_g)	
      call allocate(qs_remap, drhodp%mesh, "RemappedQS")
      call safe_set(state,qs_remap,q_s)	
    endif
    
    if (has_scalar_field(state,"TotalWaterQ") .and. (have_liquid.or.have_ice)) then
      if (have_liquid) call addto (qv_remap,qc_remap,scale=-1.)
      if (have_liquid) call addto (qv_remap,qr_remap,scale=-1.)
      if (have_ice) call addto (qv_remap,qi_remap,scale=-1.)
      if (have_ice) call addto (qv_remap,qg_remap,scale=-1.)
      if (have_ice) call addto (qv_remap,qs_remap,scale=-1.)
    endif
 
    if (constant_cp_cv) then
      call make_giraldo_quantities_1mat_cst (state,qc_p_remap,qc_v_remap)
    else
      call make_giraldo_quantities_1mat (state,have_vapour,have_liquid,have_ice, &
                   qv_remap,qc_remap,qr_remap,qi_remap,qg_remap,qs_remap, &
		   qc_p_remap,qc_v_remap)
    endif
    
    do node=1,node_count(drhodp)
      r=node_val(qc_p_remap,node)-node_val(qc_v_remap,node)
      r_cp=r / node_val(qc_p_remap,node)
      gamm=node_val(qc_p_remap,node)/node_val(qc_v_remap,node)
      exn=(node_val(pressure_remap,node)/p_0)**(r_cp)
      
      select case (thermal_variable)
      case(1)
       temperature_node=node_val(energy_remap,node)/node_val(qc_v_remap,node)
       pottem_node=temperature_node/exn
      case(2)
       temperature_node=node_val(energy_remap,node)
       pottem_node=temperature_node/exn
      case(3)
       temperature_node=node_val(energy_remap,node)*exn
       pottem_node=node_val(energy_remap,node)
      case(5)
       pottem_node=node_val(energy_remap,node)/node_val(density_remap,node)
       temperature_node=pottem_node*exn
      end select 
       
      drhodp_node=r*temperature_node
      select case (thermal_variable)
      case(3, 5)
    	drhodp_node=drhodp_node*gamm
      end select  
      rhog_node=node_val(pressure_remap,node)/(temperature_node*r)
      c_node=sqrt(gamm*r*temperature_node)
      
      i_node=1.
      if (have_liquid) i_node=i_node+(node_val(qc_remap,node)+node_val(qr_remap,node))*rhog_node/rho_w
      if (have_ice) i_node=i_node+(node_val(qi_remap,node)+node_val(qg_remap,node)+node_val(qs_remap,node))*rhog_node/rho_i
      i2_node=1.
      if (have_liquid) i2_node=i2_node-(node_val(qc_remap,node)+node_val(qr_remap,node))*node_val(density_remap,node)/rho_w
      if (have_ice) i2_node=i2_node-(node_val(qi_remap,node)+node_val(qg_remap,node)+node_val(qs_remap,node))*node_val(density_remap,node)/rho_i
      
      call set(drhodp, node, drhodp_node)
      call set(incompfix, node, i_node)
      call set(incompfix2, node, i2_node)
      call set(sound_speed_remap, node, c_node)
      call set(temperature_remap, node, temperature_node)
      call set(potentialtem_remap, node, pottem_node)
    end do

    call scale(drhodp,incompfix)
    call scale(drhodp,incompfix)
    call invert(drhodp)
    
    if (present(temperature)) &
    	call safe_set(state,temperature,temperature_remap)
    if (present(potentialtem)) &
    	call safe_set(state,potentialtem,potentialtem_remap)
    
    if (present(density_pottem)) then
      call allocate(density_pottem_remap, drhodp%mesh, "RemappedDensityPotentialTemperature")
    
      do node=1,node_count(drhodp)
        pottem_node=node_val(potentialtem_remap,node)
	qt_node=node_val(qv_remap,node)
	if (have_liquid) qt_node=qt_node+node_val(qc_remap,node)+node_val(qr_remap,node)
	if (have_ice) qt_node=qt_node+node_val(qi_remap,node)+node_val(qg_remap,node)+node_val(qs_remap,node)

        dpt_node=pottem_node*(1. + node_val(qv_remap,node)*(c_p_v-c_v_v)/(c_p-c_v))/(1. + qt_node)

        call set(density_pottem_remap, node, dpt_node)	      
      end do

      call safe_set(state,density_pottem,density_pottem_remap)        
      call deallocate(density_pottem_remap)
    end if  
    
    if (present(saturation)) then
      call allocate(saturation_remap, drhodp%mesh, "RemappedSaturation")
      call allocate(dqsaturation_remap, drhodp%mesh, "RemappedDqSaturation")
    
      do node=1,node_count(drhodp)
        temperature_node=node_val(temperature_remap,node)
        pressure_node=node_val(pressure_remap,node)

        saturation_node=cal_qsat(temperature_node,pressure_node)
        dqsaturation_node=cal_dqsatdt(temperature_node,pressure_node)

        call set(saturation_remap, node, saturation_node)	      
        call set(dqsaturation_remap, node, dqsaturation_node)	
      end do

      call safe_set(state,saturation,saturation_remap)        
      if (present(dqsaturation)) call safe_set(state,dqsaturation,dqsaturation_remap)
      call deallocate(saturation_remap)
      call deallocate(dqsaturation_remap)		      
    end if  

    if (present(density)) then
      assert(density%mesh==drhodp%mesh)

      do node=1,node_count(drhodp)
        pressure_node=node_val(pressure_remap,node)
	temperature_node=node_val(temperature_remap,node)
	rhog_node=pressure_node / (temperature_node*(node_val(qc_p_remap,node)-node_val(qc_v_remap,node)))
	i_node=node_val(incompfix,node)
        
        density_node=rhog_node/i_node  
        call set(density, node, density_node)	    
      end do
    endif
    
    if (present(pressure)) then
      assert(pressure%mesh==drhodp%mesh)

      do node=1,node_count(drhodp)
        density_node=node_val(density_remap,node)
        temperature_node=node_val(temperature_remap,node)
	i2_node=node_val(incompfix2,node)
	
        rhog_node=density_node/i2_node
        pressure_node=rhog_node*temperature_node*(node_val(qc_p_remap,node)-node_val(qc_v_remap,node))
        call set(pressure, node, pressure_node)       
      end do
    endif
      
    if(present(qc_p)) call safe_set(state,qc_p,qc_p_remap)
    if(present(qc_v)) call safe_set(state,qc_v,qc_v_remap)
    if(present(sound_speed)) call safe_set(state,sound_speed,sound_speed_remap)
    
    call deallocate(temperature_remap)
    call deallocate(potentialtem_remap)
    call deallocate(energy_remap)
    call deallocate(pressure_remap)
    call deallocate(density_remap)
    
    call deallocate(incompfix)
    call deallocate(incompfix2)
    call deallocate(qc_p_remap)
    call deallocate(qc_v_remap)
    call deallocate(sound_speed_remap)
    
    if (have_vapour) then
      call deallocate(qv_remap)
    endif
    if (have_liquid) then
      call deallocate(qc_remap)
      call deallocate(qr_remap)
    endif
    if (have_ice) then
      call deallocate(qi_remap)
      call deallocate(qg_remap)
      call deallocate(qs_remap)
    endif
      
  end subroutine compressible_eos_giraldo_1mat
  
  subroutine make_giraldo_quantities_1mat_cst (state,qc_p,qc_v)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: qc_p,qc_v
        
    call set (qc_p,c_p)
    call set (qc_v,c_v)
          
  end subroutine make_giraldo_quantities_1mat_cst

  subroutine make_giraldo_quantities_1mat(state,have_qv,have_ql,have_qi,q_v,q_c,q_r,q_i,q_g,q_s, &
                               qc_p,qc_v)

    ! Here, qc_p and qc_v are mixture averaged heat capacities: qc_p=Sum_i(q_i*c_p_i) with i the consituents of the mixture.
   
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: qc_p,qc_v,q_v,q_c,q_r,q_i,q_g,q_s
    logical, intent(in) :: have_qv, have_ql,have_qi
    type(scalar_field) :: qc_p_local,qc_v_local, tmp
    type(scalar_field) :: ql_local, qi_local, qv_local, qd_local
    character(len=OPTION_PATH_LEN) :: eos_path
    integer :: i,stat
    
    call zero(qc_p)
    call zero(qc_v)

    ! Initialize local arrays
    call allocate(qd_local,qc_p%mesh,"DryAirFraction")
    call set(qd_local,1.)
    
    if (have_qv) then
      call allocate(qv_local,qc_p%mesh,"VapourFraction")    
      call safe_set(state,qv_local,q_v)			! Water vapour fraction
      call addto(qd_local,qv_local,scale=-1.0)
    endif
    
    if (have_ql) then
      call allocate(ql_local,qc_p%mesh,"CloudFraction")    
      call safe_set(state,ql_local,q_c)	       
      
      call allocate(tmp,qc_p%mesh,"Tmp")    
      call safe_set(state,tmp,q_r)	       
      call addto(ql_local,tmp)	
      call deallocate(tmp)  
      call addto(qd_local,ql_local,scale=-1.0)
    endif
    
    if (have_qi) then
      call allocate(qi_local,qc_p%mesh,"IceFraction")    
      call safe_set(state,qi_local,q_i)	
             
      call allocate(tmp,qc_p%mesh,"Tmp")    
      call safe_set(state,tmp,q_g)
      call addto(qi_local,tmp)
      call safe_set(state,tmp,q_s)
      call addto(qi_local,tmp)
      call deallocate(tmp)  
      call addto(qd_local,qi_local,scale=-1.0)
    endif

    ! Compute cp cv
    call addto(qc_p,qd_local,scale=c_p)
    call addto(qc_v,qd_local,scale=c_v)
    
    if (have_qv) then
      call addto(qc_p,qv_local,scale=c_p_v)
      call addto(qc_v,qv_local,scale=c_v_v)
    endif
    
    if (have_ql) then
      call addto(qc_p,ql_local,scale=c_p_l)
      call addto(qc_v,ql_local,scale=c_v_l)
    endif

    if (have_qi) then
      call addto(qc_p,qi_local,scale=c_p_i)
      call addto(qc_v,qi_local,scale=c_v_i)
    endif
    
    call deallocate(qd_local)
    if (have_qv) call deallocate(qv_local)
    if (have_ql) call deallocate(ql_local)
    if (have_qi) call deallocate(qi_local)
    
  end subroutine make_giraldo_quantities_1mat

  subroutine compressible_eos_atham(state, &
  		output,output_name,thermal,thermal_variable,output2,qc_p,qc_v,getold)
    ! Routines to invert the ATHAM equation of state in various
    ! permuations
    integer, intent(in) :: thermal_variable
    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), intent(inout) :: output, thermal
    character(len=*), intent(in) :: output_name
    type(scalar_field), intent(inout), optional :: output2, qc_p, qc_v
    logical, intent(in), optional :: getold
    
    ! locals
    integer :: stat
    type(scalar_field), pointer :: pressure,density
    type(scalar_field) :: thermal_remap, pressure_remap, density_remap, output_remap, &
         qc_solid,scaledq,qc_p_local,qc_v_local

    real :: p0,TV_node,saturation_node,dqsaturation_node,esat_node,rho_node,p_node,qc_p_node, &
         qc_v_node,qc_solid_node, scaledq_node, drhodp_node,temp_node
    integer :: node
    character (len=10) :: old="          "

    if (present(getold)) then
       if (getold) then
          old="Old"
       else
          old="   "
       end if
    else
       old="   "
    end if

    call get_option(trim(state(1)%option_path)//'/compressible/ATHAM/reference_pressure', &
                    p0, default=1.0e5)

    density=>extract_scalar_field(state(1),trim(old)//"Density",stat=stat)
    if (stat /=0 ) then
       FLAbort("No density in bulk state for ATHAM EoS")
    end if
    pressure=>extract_scalar_field(state(1),trim(old)//"Pressure",stat=stat)
    if (stat /=0 ) then
       FLAbort("No pressure in bulk state for ATHAM EoS")
    end if
 
    call allocate(thermal_remap,output%mesh,"Remeshed"//trim(thermal%name))
    call allocate(density_remap,output%mesh,"RemeshedDensity")
    call allocate(pressure_remap,output%mesh,"RemeshedPressure")

    call set(pressure_remap,p0)
    call set(thermal_remap,275.0)
    call set(density_remap,1.0)

    call safe_set(state(1),density_remap,density)
    call safe_set(state(1),pressure_remap,pressure)
    call safe_set(state(1),thermal_remap,thermal)

    call allocate(qc_p_local,output%mesh,"LocalCP")
    call allocate(qc_v_local,output%mesh,"LocalCV")
    call allocate(qc_solid,output%mesh,"LocalCSolid")
    call allocate(scaledq,output%mesh,"LocalScaledQ")
    
    call make_atham_quantities &
	     (state,thermal,qc_p_local,qc_v_local,qc_solid=qc_solid,scaledq=scaledq,getold=getold)

    if (output_name=="Pressure") &
         call pressure_atham(output,qc_p_local,qc_v_local,qc_solid,scaledq,density_remap,thermal_remap)
    if (output_name=="Density") &
         call density_atham(output,qc_p_local,qc_v_local,qc_solid,scaledq,pressure_remap,thermal_remap)
    if (output_name=="Temperature") &
         call temperature_atham(output,qc_p_local,qc_v_local,qc_solid,scaledq,thermal_remap,pressure_remap)
    if (output_name=="Saturation" .and. present(output2)) &
         call saturation_atham(output,output2,qc_p_local,qc_v_local,qc_solid,scaledq,thermal_remap,pressure_remap)
    if (present(output2)) then
       if (output_name=="Pressure") then
          call drhodp_atham(output2,qc_p_local,qc_v_local,qc_solid,scaledq,density_remap,pressure_remap,&
               thermal_remap)
       else if (output_name=="Density") then
          call drhodp_atham(output2,qc_p_local,qc_v_local,qc_solid,scaledq,density_remap,pressure_remap,&
               thermal_remap)
       else
          call drhodp_atham(output2,qc_p_local,qc_v_local,qc_solid,scaledq,density_remap,pressure_remap,&
               thermal_remap)
       end if
    end if
    
    if(present(qc_p)) call safe_set(state(1),qc_p,qc_p_local)
    if(present(qc_v)) call safe_set(state(1),qc_v,qc_v_local)   

    call deallocate(pressure_remap)
    call deallocate(density_remap)
    call deallocate(thermal_remap)
    call deallocate(qc_p_local)
    call deallocate(qc_v_local)
    call deallocate(qc_solid)
    call deallocate(scaledq)

  contains

      subroutine pressure_atham(p,qc_p,qc_v,qc_solid,scaledq,rho,TV)
        type(scalar_field), intent(inout) :: p
        type(scalar_field), intent(inout):: rho,TV,qc_p,qc_v,qc_solid,scaledq
	type(scalar_field) :: qc_p_remap, qc_v_remap, qc_s_remap, scaledq_remap

        type(scalar_field) :: rhs

        real :: fp_node

        integer :: nlin
	
        call allocate(qc_p_remap,output%mesh,"RemappedCP")
	call safe_set(state(1),qc_p_remap,qc_p)
        call allocate(qc_v_remap,output%mesh,"RemappedCV")
	call safe_set(state(1),qc_v_remap,qc_v)
        call allocate(qc_p_remap,output%mesh,"RemappedCS")
	call safe_set(state(1),qc_s_remap,qc_solid)
        call allocate(scaledq_remap,output%mesh,"RemappedDQ")
	call safe_set(state(1),scaledq_remap,scaledq)

        select case(thermal_variable)
        case(1)
           do node=1,node_count(p)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              p_node=rho_node*(qc_p_node-qc_v_node)*TV_node&
                   /((1.0-rho_node*scaledq_node)*(qc_v_node+qc_solid_node))
              
              call set(p,node,p_node)
	    enddo
	   
        case(2)
           call allocate(rhs,p%mesh,"RHS")

           do node=1,node_count(p)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              p_node=p0**(-(qc_p_node-qc_v_node)/qc_v_node)&
                      *(rho_node*(qc_p_node-qc_v_node)*TV_node)&
                      **(qc_p_node/qc_v_node)

              if (abs(scaledq_node)>1e-8) then

                 temp_node=rho_node*(qc_p_node-qc_v_node)&
                   *(qc_p_node+qc_solid_node)&
                   *TV_node/(1.0-rho_node*scaledq_node)
                 fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)

                 do nlin=1,200
                    p_node=p_node&
                         -(fp_node-temp_node)&
                         /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_v_node+qc_solid_node)
                    
                    fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_p_node+qc_solid_node)
                    
                    if (sqrt((fp_node-temp_node)**2)<1.0e-10) then
                       exit
                    end if
                 
                 end do
              end if
              call set(p,node,p_node)
              call set(rhs,node,temp_node)

           end do

           call deallocate(rhs)
	   
        case(3)
           
           call allocate(rhs,p%mesh,"RHS")

           do node=1,node_count(p)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              p_node=rho_node*(qc_p_node-qc_v_node)*TV_node	&
		 /(1.0-rho_node*scaledq_node)

              call set(p,node,p_node)
              call set(rhs,node,temp_node)

           end do

           call deallocate(rhs)

        case(4)

           call allocate(rhs,p%mesh,"RHS")

           do node=1,node_count(p)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              p_node=p0**(-(qc_p_node-qc_v_node)/qc_v_node)&
                      *(rho_node*(qc_p_node-qc_v_node)*TV_node/(qc_p_node+qc_solid_node))&
                      **(qc_p_node/qc_v_node)

              if (abs(scaledq_node)>1e-8) then

                 temp_node=rho_node*(qc_p_node-qc_v_node)&
                   *TV_node/(1.0-rho_node*scaledq_node)
                 fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)

                 do nlin=1,200
                    p_node=p_node&
                         -(fp_node-temp_node)&
                         /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_v_node+qc_solid_node)
                    
                    fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_p_node+qc_solid_node)
                    
                    if (sqrt((fp_node-temp_node)**2)<1.0e-10) then
                       exit
                    end if
                 
                 end do
              end if
              call set(p,node,p_node)
              call set(rhs,node,temp_node)


           end do

           call deallocate(rhs)
        end select
	
	call deallocate(qc_p_remap)
	call deallocate(qc_v_remap)
	call deallocate(qc_s_remap)
	call deallocate(scaledq_remap)

      end subroutine pressure_atham

      subroutine density_atham(rho,qc_p,qc_v,qc_solid,scaledq,p,TV)
        type(scalar_field), intent(inout) :: rho
        type(scalar_field), intent(inout):: p,TV,qc_p,qc_v,qc_solid,scaledq
	type(scalar_field) :: qc_p_remap, qc_v_remap, qc_s_remap, scaledq_remap
	
        call allocate(qc_p_remap,output%mesh,"RemappedCP")
	call safe_set(state(1),qc_p_remap,qc_p)
        call allocate(qc_v_remap,output%mesh,"RemappedCV")
	call safe_set(state(1),qc_v_remap,qc_v)
        call allocate(qc_p_remap,output%mesh,"RemappedCS")
	call safe_set(state(1),qc_s_remap,qc_solid)
        call allocate(scaledq_remap,output%mesh,"RemappedDQ")
	call safe_set(state(1),scaledq_remap,scaledq)

        select case(thermal_variable)
        case(1)
           do node=1,node_count(rho)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              rho_node=(qc_v_node+qc_solid_node)*p_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +(qc_v_node+qc_solid_node)*p_node*scaledq_node)
              
              call set(rho,node,rho_node)
           end do

        case(2)
           do node=1,node_count(rho)
              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              rho_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)&
                   /((qc_p_node-qc_v_node)*(qc_p_node+qc_solid_node)*TV_node&
                   +p_node*scaledq_node&
                   *((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node))
              
              call set(rho,node,rho_node)

           end do
	   
        case(3)
           do node=1,node_count(rho)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              rho_node=p_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +(qc_v_node+qc_solid_node)*p_node*scaledq_node)
              
              call set(rho,node,rho_node)
           end do

        case(4)
           do node=1,node_count(rho)
              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              rho_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +p_node*scaledq_node&
                   *((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node))
              
              call set(rho,node,rho_node)

           end do

        end select
	
	call deallocate(qc_p_remap)
	call deallocate(qc_v_remap)
	call deallocate(qc_s_remap)
	call deallocate(scaledq_remap)

      end subroutine density_atham

      subroutine temperature_atham(T,qc_p,qc_v,qc_solid,scaledq,TV,p)
        type(scalar_field), intent(inout) :: T
        type(scalar_field), intent(inout):: p,TV,qc_p,qc_v,qc_solid,scaledq
	type(scalar_field) :: qc_p_remap, qc_v_remap, qc_s_remap, scaledq_remap
	
        call allocate(qc_p_remap,output%mesh,"RemappedCP")
	call safe_set(state(1),qc_p_remap,qc_p)
        call allocate(qc_v_remap,output%mesh,"RemappedCV")
	call safe_set(state(1),qc_v_remap,qc_v)
        call allocate(qc_p_remap,output%mesh,"RemappedCS")
	call safe_set(state(1),qc_s_remap,qc_solid)
        call allocate(scaledq_remap,output%mesh,"RemappedDQ")
	call safe_set(state(1),scaledq_remap,scaledq)

        select case(thermal_variable)
        case(1)

           call set(T,qc_v_remap)
           call addto(T,qc_s_remap)
           call invert(T)
           call scale(T,TV)
	   
        case(2)

           call set(T,TV)
	   
        case(3)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)

              temp_node=TV_node*(qc_p_node+qc_solid_node)&
                   /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node&
                   +qc_solid_node)
                   
              call set(T,node,temp_node)
              
           end do
	   
        case(4)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)

              temp_node=TV_node&
                   /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node&
                   +qc_solid_node)
                   

              call set(T,node,temp_node)
              
           end do

        end select
	
        call deallocate(qc_p_remap)
        call deallocate(qc_v_remap)
        call deallocate(qc_s_remap)
        call deallocate(scaledq_remap)

      end subroutine temperature_atham

      subroutine potential_temperature_atham(T,qc_p,qc_v,qc_solid,scaledq,TV,p)
        type(scalar_field), intent(inout) :: T
        type(scalar_field), intent(inout):: p,TV,qc_p,qc_v,qc_solid,scaledq
	type(scalar_field) :: qc_p_remap, qc_v_remap, qc_s_remap, scaledq_remap
	
        call allocate(qc_p_remap,output%mesh,"RemappedCP")
	call safe_set(state(1),qc_p_remap,qc_p)
        call allocate(qc_v_remap,output%mesh,"RemappedCV")
	call safe_set(state(1),qc_v_remap,qc_v)
        call allocate(qc_p_remap,output%mesh,"RemappedCS")
	call safe_set(state(1),qc_s_remap,qc_solid)
        call allocate(scaledq_remap,output%mesh,"RemappedDQ")
	call safe_set(state(1),scaledq_remap,scaledq)

        select case(thermal_variable)
        case(1)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)

              temp_node=TV_node*&
                   ((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node &
                   +qc_solid_node)/((qc_v_node+qc_solid_node)*(qc_p_node+qc_solid_node))
                   
              call set(T,node,temp_node)
              
           end do

        case(2)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)

              temp_node=TV_node*&
                   ((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node&
                   +qc_solid_node)/(qc_p_node+qc_solid_node)
                   
              call set(T,node,temp_node)
              
           end do
	   
        case(3)

           call set(T,TV)
	   
        case(4)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)

              temp_node=TV_node*(qc_p_node+qc_solid_node)
                   
              call set(T,node,temp_node)
              
           end do

        end select
	
        call deallocate(qc_p_remap)
        call deallocate(qc_v_remap)
        call deallocate(qc_s_remap)
        call deallocate(scaledq_remap)
	
      end subroutine potential_temperature_atham

      subroutine saturation_atham(T,dT,qc_p,qc_v,qc_solid,scaledq,TV,p)
        type(scalar_field), intent(inout) :: T,dT
        type(scalar_field), intent(inout):: p,TV,qc_p,qc_v,qc_solid,scaledq
	type(scalar_field) :: qc_p_remap, qc_v_remap, qc_s_remap, scaledq_remap
	
        call allocate(qc_p_remap,output%mesh,"RemappedCP")
	call safe_set(state(1),qc_p_remap,qc_p)
        call allocate(qc_v_remap,output%mesh,"RemappedCV")
	call safe_set(state(1),qc_v_remap,qc_v)
        call allocate(qc_p_remap,output%mesh,"RemappedCS")
	call safe_set(state(1),qc_s_remap,qc_solid)
        call allocate(scaledq_remap,output%mesh,"RemappedDQ")
	call safe_set(state(1),scaledq_remap,scaledq)

        select case(thermal_variable)
        case(1)

           do node=1,node_count(T)
              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
 
              temp_node=TV_node/(qc_v_node+qc_solid_node)
		   
	      saturation_node=cal_qsat(temp_node,p_node)
	      dqsaturation_node= cal_dqsatdt(temp_node,p_node)
		                      
              call set(T,node,saturation_node)
              call set(dT,node,dqsaturation_node)	      
	   enddo
	   
        case(2)

           do node=1,node_count(T)
              temp_node=node_val(TV,node)
              p_node=node_val(p,node)
		   
	      saturation_node=cal_qsat(temp_node,p_node)
	      dqsaturation_node= cal_dqsatdt(temp_node,p_node)
	                        
              call set(T,node,saturation_node)
              call set(dT,node,dqsaturation_node)  
	   enddo
	   	
        case(3)

           do node=1,node_count(T)
              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)

              temp_node=TV_node*(qc_p_node+qc_solid_node)&
                   /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node&
                   +qc_solid_node)
		   
	      saturation_node=cal_qsat(temp_node,p_node)
	      dqsaturation_node= cal_dqsatdt(temp_node,p_node)
                   
              call set(T,node,saturation_node)
              call set(dT,node,dqsaturation_node)
           end do
	   
        case(4)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)

              temp_node=TV_node&
                   /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node&
                   +qc_solid_node)
		   
	      saturation_node=cal_qsat(temp_node,p_node)
	      dqsaturation_node= cal_dqsatdt(temp_node,p_node)
                   
              call set(T,node,saturation_node)
              call set(dT,node,dqsaturation_node)		   
           end do

        end select
	
        call deallocate(qc_p_remap)
        call deallocate(qc_v_remap)
        call deallocate(qc_s_remap)
        call deallocate(scaledq_remap)

      end subroutine saturation_atham

      subroutine drhodp_atham(drhodp,qc_p,qc_v,qc_solid,scaledq,rho,p,TV)
        type(scalar_field), intent(inout) :: drhodp
        type(scalar_field), intent(inout):: rho,p,TV,qc_p,qc_v,qc_solid,scaledq
	type(scalar_field) :: qc_p_remap, qc_v_remap, qc_s_remap, scaledq_remap
	
        call allocate(qc_p_remap,output%mesh,"RemappedCP")
	call safe_set(state(1),qc_p_remap,qc_p)
        call allocate(qc_v_remap,output%mesh,"RemappedCV")
	call safe_set(state(1),qc_v_remap,qc_v)
        call allocate(qc_p_remap,output%mesh,"RemappedCS")
	call safe_set(state(1),qc_s_remap,qc_solid)
        call allocate(scaledq_remap,output%mesh,"RemappedDQ")
	call safe_set(state(1),scaledq_remap,scaledq)
              
        select case(thermal_variable)
	
        case(1)
           do node=1,node_count(drhodp)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)


              drhodp_node=(qc_v_node+qc_solid_node)&
                   *(qc_p_node-qc_v_node)*TV_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +(qc_v_node+qc_solid_node)*p_node*scaledq_node)**2

              call set(drhodp, node, drhodp_node)
           end do
        case(2)
           do node=1,node_count(drhodp)


              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              temp_node=(p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)
              drhodp_node=(temp_node*qc_v_node+qc_solid_node)*(qc_p_node-qc_v_node)&
                   *(qc_p_node+qc_solid_node)*TV_node&
                   /((qc_p_node-qc_v_node)*(qc_p_node+qc_solid_node)*TV_node&
                   +p_node*(temp_node*qc_p_node+qc_solid_node)*scaledq_node)**2

              call set(drhodp, node, drhodp_node)
           end do
        case(3)
           do node=1,node_count(drhodp)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)


              drhodp_node=(qc_p_node-qc_v_node)*TV_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +(qc_v_node+qc_solid_node)*p_node*scaledq_node)**2

              call set(drhodp, node, drhodp_node)
           end do
        case(4)
           do node=1,node_count(drhodp)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p_remap,node)
              qc_v_node=node_val(qc_v_remap,node)
              qc_solid_node=node_val(qc_s_remap,node)
              scaledq_node=node_val(scaledq_remap,node)

              temp_node=(p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)
              drhodp_node=(temp_node*qc_v_node+qc_solid_node)*(qc_p_node-qc_v_node)&
                   *TV_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +p_node*(temp_node*qc_p_node+qc_solid_node)*scaledq_node)**2

              call set(drhodp, node, drhodp_node)
           end do

        end select
	
        call deallocate(qc_p_remap)
        call deallocate(qc_v_remap)
        call deallocate(qc_s_remap)
        call deallocate(scaledq_remap)

      end subroutine drhodp_atham

  end subroutine compressible_eos_atham

  subroutine make_atham_quantities(state,thermal,qc_p,qc_v,qc_solid,scaledq,getold)

    logical, intent(in), optional :: getold
    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), intent(inout) :: qc_p,qc_v,qc_solid,scaledq,thermal

    type(scalar_field), pointer :: fraction, fractionc, fractionr, fraction_dry
    type(scalar_field) :: qc_p_local,qc_v_local,&
         qc_solid_local,scaledq_local
    character(len=OPTION_PATH_LEN) :: eos_path, micro_path
    character(len=FIELD_NAME_LEN)  :: mom_name
    character (len=10) :: old="          "
    integer :: i,stat
    real :: c_vi,c_pi,rho

    if (present(getold)) then
       if (getold) then
          old="Old"
       else
          old="   "
       end if
    else
       old="   "
    end if

    call allocate(qc_p_local,thermal%mesh,"LocalqC_p")
    call allocate(qc_v_local,thermal%mesh,"LocalqC_v")
    call allocate(qc_solid_local,thermal%mesh,"LocalqC_solid")
    call allocate(scaledq_local,thermal%mesh,"LocalScaledQ")
    call allocate(fraction_dry,thermal%mesh,"DryAirQ")

    call zero(qc_p_local)
    call zero(qc_v_local)
    call zero(qc_solid_local)
    call zero(scaledq_local)
    call set (fraction_dry,1.0)

    do i=2,size(state)
       if (has_scalar_field(state(i),"MassFraction")) then
          fraction=>extract_scalar_field(state(i),trim(old)//"MassFraction")
	  call addto(fraction_dry,fraction,scale=-1.)
	  
          eos_path = trim(state(i)%option_path)//'/equation_of_state'
          call get_option(trim(eos_path)//'/compressible/ATHAM/C_P', &
               c_pi, stat=stat)
          if (stat /= 0) then
             ewrite(0,*) "Selected compressible eos but not specified C_P in "//state(i)%name//"."
             cycle
          end if
          call get_option(trim(eos_path)//'/compressible/ATHAM/C_V', &
               c_vi, stat=stat)
          if (stat /= 0) then 
             ewrite(0,*) "Selected compressible eos but not specified C_V in "//state(i)%name//"."
             cycle
          end if

          if (have_option(trim(eos_path)//'/compressible/ATHAM/density')) then
             call get_option(trim(eos_path)//'/compressible/ATHAM/density',rho,stat)
             call addto(qc_solid_local,fraction,scale=c_vi)
             call addto(scaledq_local,fraction,scale=1.0/rho)
          else
             call addto(qc_v_local,fraction,scale=c_vi)
             call addto(qc_p_local,fraction,scale=c_pi)
          end if
       end if
    end do
    call addto(qc_v_local,fraction_dry,scale=c_v)
    call addto(qc_p_local,fraction_dry,scale=c_p)

    call safe_set(state(1),qc_p,qc_p_local)
    call safe_set(state(1),qc_v,qc_v_local)
    call safe_set(state(1),qc_solid,qc_solid_local)
    call safe_set(state(1),scaledq,scaledq_local)

    call deallocate(qc_p_local)
    call deallocate(qc_v_local)
    call deallocate(qc_solid_local)
    call deallocate(scaledq_local)
    call deallocate(fraction_dry)
    
  end subroutine make_atham_quantities
  
  subroutine compressible_eos_foam(state, eos_path, drhodp, &
    density, pressure)
    ! Foam EoS Used with compressible simulations of liquid drainage in foams.
    ! It describes the liquid content in the foam as the product of the  Plateau 
    ! border cross sectional area and the local Plateau  border length per unit volume (lambda).
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure

    ! locals
    integer :: pstat, dstat
    type(scalar_field), pointer :: pressure_local, density_local, drainagelambda_local
    real :: atmospheric_pressure
    type(scalar_field) :: pressure_remap, density_remap, drainagelambda_remap

    call zero(drhodp)

    pressure_local => extract_scalar_field(state,'Pressure', stat=pstat)

    drainagelambda_local => extract_scalar_field(state,'DrainageLambda')

    call allocate(drainagelambda_remap, drhodp%mesh, 'RemappedDrainageLambda')
    call remap_field(drainagelambda_local, drainagelambda_remap)

    call addto(drhodp, drainagelambda_remap)

    call deallocate(drainagelambda_remap)

    if(present(density)) then
      if (pstat==0) then
        assert(density%mesh==drhodp%mesh)

        call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
        call remap_field(pressure_local, pressure_remap)

        call set(density, atmospheric_pressure)
        call addto(density, pressure_remap)
        call scale(density, drhodp)

        call deallocate(pressure_remap)
      else
        FLExit('No Pressure in material_phase::'//trim(state%name))
      end if
    end if

    if(present(pressure)) then
      density_local=>extract_scalar_field(state,'Density',stat=dstat)
      if (dstat==0) then
        assert(pressure%mesh==drhodp%mesh)

        call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call allocate(density_remap, drhodp%mesh, "RemappedDensity")
        call remap_field(density_local, density_remap)

        call set(pressure, drhodp)
        call invert(pressure)
        call scale(pressure, density_remap)

        call deallocate(density_remap)
      else
        FLExit('No Density in material_phase::'//trim(state%name))
      end if
    end if
        

  end subroutine compressible_eos_foam
        
  subroutine compressible_material_eos(state,materialdensity,materialpressure,materialdrhodp)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: materialdensity, &
                                             materialpressure, materialdrhodp

    !locals
    integer :: stat, gstat, cstat
    type(scalar_field), pointer :: pressure, materialenergy, materialdensity_local
    character(len=4000) :: thismaterial_phase, eos_path
    real :: reference_density, ratio_specific_heats
    real :: bulk_sound_speed_squared, atmospheric_pressure
    type(scalar_field) :: drhodp

    ewrite(1,*) 'Entering compressible_material_eos'

    if (present(materialdensity)) then
      call allocate(drhodp, materialdensity%mesh, 'Gradient of density wrt pressure')
    else if (present(materialpressure)) then
      call allocate(drhodp, materialpressure%mesh, 'Gradient of density wrt pressure')
    else if (present(materialdrhodp)) then
      call allocate(drhodp, materialdrhodp%mesh, 'Gradient of density wrt pressure')
    else
      FLAbort("No point in being in here if you don't want anything out.")
    end if

    thismaterial_phase = '/material_phase::'//trim(state%name)
    eos_path = trim(thismaterial_phase)//'/equation_of_state'

    if(have_option(trim(eos_path)//'/compressible')) then

      if(have_option(trim(eos_path)//'/compressible/stiffened_gas')) then
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/reference_density', &
                        reference_density, default=0.0)
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/ratio_specific_heats', &
                        ratio_specific_heats, stat=gstat)
        if(gstat/=0) then
          ratio_specific_heats=1.0
        end if
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/bulk_sound_speed_squared', &
                        bulk_sound_speed_squared, stat = cstat)
        if(cstat/=0) then
          bulk_sound_speed_squared=0.0
        end if
        if((gstat/=0).and.(cstat/=0)) then
          FLExit("Must set either a bulk_sound_speed_squared or a ratio_specific_heats.")
        end if
        materialenergy=>extract_scalar_field(state,'MaterialInternalEnergy',stat=stat)
        if(stat==0) then   ! we have an internal energy field
          drhodp%val = 1.0/( bulk_sound_speed_squared + (ratio_specific_heats - 1.0)*materialenergy%val )
        else               ! we don't have an internal energy field
          call set(drhodp, 1.0/bulk_sound_speed_squared)
        end if

        if(present(materialdensity)) then
          ! calculate the materialdensity
          ! materialdensity can equal materialdensity in state depending on how this
          ! subroutine is called
          pressure=>extract_scalar_field(state,'Pressure',stat=stat)
          if (stat==0) then
            call get_option(trim(pressure%option_path)//'/prognostic/atmospheric_pressure', &
                            atmospheric_pressure, default=0.0)
            materialdensity%val = drhodp%val*(pressure%val + atmospheric_pressure &
                                   + bulk_sound_speed_squared*reference_density)
          else
            FLExit('No Pressure in material_phase::'//trim(state%name))
          end if
        end if

        if(present(materialpressure)) then
          ! calculate the materialpressure using the eos and the calculated (probably prognostic)
          ! materialdensity
          ! materialpressure /= bulk pressure
          materialdensity_local=>extract_scalar_field(state,'MaterialDensity',stat=stat)
          if (stat==0) then
            materialpressure%val = materialdensity_local%val/drhodp%val &
                                  - bulk_sound_speed_squared*reference_density
          else
            FLExit('No MaterialDensity in material_phase::'//trim(state%name))
          end if
        end if

!       else
!       ! place other compressible material eos here

      end if

!     else
!     ! an incompressible option?

    end if      

    if(present(materialdensity)) then
      ewrite_minmax(materialdensity)
    end if

    if(present(materialpressure)) then
      ewrite_minmax(materialpressure)
    end if

    if(present(materialdrhodp)) then
      materialdrhodp%val=drhodp%val
      ewrite_minmax(materialdrhodp)
    end if

    call deallocate(drhodp)

  end subroutine compressible_material_eos

  subroutine set_EOS_pressure_and_temperature(state,mesh,have_temperature,have_pottemperature,&
      have_EOSPressure,have_EOSDensity,extras)
    
    type(state_type), intent(inout) :: state
    type(scalar_field), pointer :: temperature,pottemperature,thermal, pressure, density
    type(mesh_type), pointer, intent(inout) :: mesh
    integer :: i, stat
    logical :: have_temperature,have_pottemperature,have_EOSPressure,have_EOSDensity
    type(scalar_field), dimension(:), intent(inout), target, optional :: extras

    ewrite(3,*) 'Entering set_EOS_pressure_and_temperature'

    ! Here there will be magic done to put the Equation of state pressure
    ! and insitu bulk Temperature on the density mesh in the input state

    i=1
    if (have_temperature) then
       temperature=>extract_scalar_field(state,"Temperature",stat=stat)
       if (stat == 0) call compressible_eos(state,temperature=temperature)
    else
       call allocate(extras(i),mesh,"IteratedTemperature")
       temperature=>extras(i)
       call zero(temperature)
       call compressible_eos(state,temperature=temperature)
       call allocate(extras(i+1),mesh,"OldTemperature")
       temperature=>extras(i+1)
       call zero(temperature)
       call compressible_eos(state,temperature=temperature,getold=.true.)
       i=i+2
       stat=0
    end if
    
    if (have_pottemperature) then
       pottemperature=>extract_scalar_field(state,"PotentialTemperature",stat=stat)
       if (stat == 0) call compressible_eos(state,potentialtem=pottemperature)
    else
       call allocate(extras(i),mesh,"IteratedPotentialTemperature")
       pottemperature=>extras(i)
       call zero(pottemperature)
       call compressible_eos(state,potentialtem=pottemperature)
       call allocate(extras(i+1),mesh,"OldPotentialTemperature")
       pottemperature=>extras(i+1)
       call zero(pottemperature)
       call compressible_eos(state,potentialtem=pottemperature,getold=.true.)
       i=i+2
       stat=0
    end if

    if (have_EOSPressure) then
       pressure=>extract_scalar_field(state,"EOSPressure",stat)
       if (stat == 0) call compressible_eos(state,pressure=pressure)
    else
       call allocate(extras(i),mesh,"IteratedEOSPressure")
       pressure=>extras(i)
       call zero(pressure)
       call compressible_eos(state,pressure=pressure)
       call allocate(extras(i+1),mesh,"OldEOSPressure")
       pressure=>extras(i+1)
       call zero(pressure)
       call compressible_eos(state,pressure=pressure,getold=.true.)
       i=i+2
    end if
    
    if (have_EOSDensity) then
       density=>extract_scalar_field(state,"IteratedEOSDensity",stat)
       if (stat == 0) call compressible_eos(state,density=density)
    else
       call allocate(extras(i),mesh,"IteratedEOSDensity")
       density=>extras(i)
       call zero(density)
       call compressible_eos(state,density=density)
       call allocate(extras(i+1),mesh,"OldEOSDensity")
       density=>extras(i+1)
       call zero(density)
       call compressible_eos(state,density=density,getold=.true.)
       i=i+2
    end if

    ewrite(3,*) 'End set_EOS_pressure_and_temperature'

  end subroutine set_EOS_pressure_and_temperature

  subroutine initialise_from_eos(state,field,path)
    implicit none
    
    character(len=*), intent(in) :: path
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: field
    type(scalar_field), pointer :: pressure, density, thermal
    integer :: value, nvalues
    
    nvalues=option_count(trim(path))

    do value = 0, nvalues-1
      if (have_option(trim(path)//'['//int2str(value)//']'//'/from_equation_of_state')) then

        call zero(field)
        if (trim(field%name) == "Density") then
           pressure=>extract_scalar_field(state, "Pressure")
           call get_thermo_variable(state, thermal)

           call compressible_eos (state, full_pressure=pressure, energy=thermal, &
        	      density=field)
		      
        else if (trim(field%name) == "Pressure") then
           density=>extract_scalar_field(state, "Density")
           call get_thermo_variable(state, thermal)

           call compressible_eos (state, full_density=density, energy=thermal, &
        	      pressure=field)
	
	else if(trim(field%name) == "HydrostaticReferenceDensity") then
           pressure=>extract_scalar_field(state, "HydrostaticReferencePressure")
           thermal=>extract_scalar_field(state, "HydrostaticReferencePotentialTemperature")

           call compressible_eos (state, full_pressure=pressure, energy=thermal, &
        	      density=field)
		      
        else if (trim(field%name) == "HydrostaticReferencePressure") then
           density=>extract_scalar_field(state, "HydrostaticReferenceDensity")
           call get_thermo_variable(state, thermal)

           call compressible_eos (state, full_density=density, energy=thermal, &
        	      pressure=field)
        endif
      
      endif
    end do

  end subroutine initialise_from_eos
  
  function extract_entropy_variable(state,prefix,suffix,stat) result(sfield)
    type(state_type) :: state
    character(len=*), optional :: prefix, suffix
    type(scalar_field), pointer :: sfield
    integer, optional :: stat


    character(len=OPTION_PATH_LEN) ::lprefix,lsuffix
    integer :: i,j,lstat

    character(len=*), dimension(5), parameter ::&
         entropy_names= (/ "PotentialTemperature         ",&
                           "ConservedPotentialTemperature",&
                           "InternalEnergy               ",&
                           "Temperature                  ",&
                           "ATHAMInternalEnergy          "/)

    ! Routine searches (in order) for one of the entropy names above 
    ! in the states given to it, then returns the field
                           
    if (present(prefix)) then
       lprefix=prefix
    else
       lprefix=""
    end if

    if (present(suffix)) then
       lsuffix=suffix
    else
       lsuffix=""
    end if

    do j=1,size(entropy_names)
       sfield=>extract_scalar_field(state,trim(lprefix)//trim(entropy_names(j))//trim(lsuffix),stat=lstat)
       if (present(stat)) stat=lstat
       if (lstat==0) return
    end do
    
    if (.not. present(stat)) then 
       FLAbort('Failed to find entropy field when using Atham Equation of state')
    end if

  end function extract_entropy_variable

  subroutine scale_pressure(t, state, pressure, pressure_scale)

    type(scalar_field), intent(inout) :: t, pressure, pressure_scale
    type(state_type), intent(inout) :: state
        
    integer  :: cstat, i
    real  :: c_v, c_p
    type(scalar_field) :: qc_p,qc_v
    character(len=OPTION_PATH_LEN) :: eos_path,eos_option

    call allocate(qc_p,pressure%mesh,"qc_p")
    call allocate(qc_v,pressure%mesh,"qc_v")
    
    do i = 1, size(qc_v%val)
      call set (qc_v,i,c_v)
      call set (qc_p,i,c_p)
    enddo  
        
    select case (trim(t%name))
      case('InternalEnergy')
        do i = 1, size(qc_v%val)
	  pressure_scale%val(i) = 1.
	enddo  

      case('PotentialTemperature') 
        do i = 1, size(qc_v%val)
	  pressure_scale%val(i) = 0.
	enddo  

      case('ConservedPotentialTemperature') 
        do i = 1, size(qc_v%val)
	  pressure_scale%val(i) = 0.
	enddo  
	  
      case('Temperature') 
        do i = 1, size(qc_v%val)
	  pressure_scale%val(i) = 1./qc_v%val(i)
	enddo  

      case('ATHAMInternalEnergy') 
        do i = 1, size(qc_v%val)
	  pressure_scale%val(i) = 1.
	enddo  
       
      case default
        do i = 1, size(qc_v%val)
	  pressure_scale%val(i) = 1.
	enddo 
    end select

    call deallocate(qc_p)
    call deallocate(qc_v)
    
  end subroutine scale_pressure
   
  subroutine get_thermo_variable (state,field,field_name,index)
   
   type(state_type), intent(in) :: state
   type(scalar_field), pointer, intent(out) :: field
   integer, optional, intent(out) :: index
   character(len=*), optional, intent(out) :: field_name
   
   integer :: stat
       
   field=>extract_scalar_field(state,"PotentialTemperature",stat=stat)
   if(present(index)) index=3
   if (stat /= 0) then
       field=>extract_scalar_field(state,"InternalEnergy",stat=stat)
       if(present(index)) index=1
   endif
   if (stat /= 0) then
       field=>extract_scalar_field(state,"Temperature",stat=stat)
       if(present(index)) index=2
   endif
   if (stat /= 0) then
       field=>extract_scalar_field(state,"ATHAMInternalEnergy",stat=stat)
       if(present(index)) index=4
   endif
   if (stat /= 0) then
       field=>extract_scalar_field(state,"ConservedPotentialTemperature",stat=stat)
       if(present(index)) index=5
   endif
   if (stat /= 0) then
       FLAbort("Can t find thermodynamic variable in state!")
   endif
   
   if (present(field_name)) field_name=field%name

  end subroutine get_thermo_variable
      
  subroutine get_cp_cv (c_p,c_v,c_p_v,c_v_v,c_p_l,c_v_l,c_p_i,c_v_i)
    
    real, intent(out) :: c_p,c_v
    real, intent(out), optional :: c_p_v, c_v_v, c_p_l, c_v_l, c_p_i, c_v_i
    
    character(len=OPTION_PATH_LEN) :: eos_path
    integer :: nphases, stat    
    
    nphases = option_count('/material_phase')
    
    if (nphases == 1) then
      eos_path = '/material_phase[0]/equation_of_state'
    
      if (have_option(trim(eos_path)//'/compressible')) then
     
      if (have_option(trim(eos_path)//'/compressible/ATHAM')) then
        call get_option(trim(eos_path)//'/compressible/ATHAM/C_P', &
     	     c_p, stat=stat)
        if (stat /= 0) then
           ewrite(0,*) "Selected compressible eos but not specified C_P in state #0"
        end if
        call get_option(trim(eos_path)//'/compressible/ATHAM/C_V', &
     	     c_v, stat=stat)
        if (stat /= 0) then 
           ewrite(0,*) "Selected compressible eos but not specified C_V in state #0" 
        end if
      else if (have_option(trim(eos_path)//'/compressible/giraldo')) then
        call get_option(trim(eos_path)//'/compressible/giraldo/C_P', &
         	   c_p, stat=stat)
        if (stat /= 0) then
           ewrite(0,*) "Selected compressible eos but not specified C_P in state #0" 
        end if
        call get_option(trim(eos_path)//'/compressible/giraldo/C_V', &
       	     c_v, stat=stat)
        if (stat /= 0) then 
           ewrite(0,*) "Selected compressible eos but not specified C_V in state #0" 
        end if
      else 
        ewrite(0,*) "C_P and C_V are not explicitly defined. Default values are selected"
        c_p = 1004.64
        c_v = 717.6
      endif
      
      endif

      if (present(c_v_v)) c_v_v=1423.5
      if (present(c_p_v)) c_p_v=1885.
      if (present(c_v_l)) c_v_l=4186.
      if (present(c_p_l)) c_p_l=4186.
      if (present(c_v_i)) c_v_i=2050.
      if (present(c_p_i)) c_p_i=2050.
      
    else
       FLAbort("get_cp_cv work only with a single material phase")      
    endif
    
  end subroutine get_cp_cv
  
  function cal_esat (tem)
    real :: tem, es, cal_esat
    real, parameter :: k1=54.842763, k2=-6763.22, k3=-4.21, k4=0.000367
    real, parameter :: k5=0.0415, k6=53.878, k7=-1331.22, k8=-9.44523, k9=0.014025
    
    es=k1 + k2/tem + k3*log(tem) + k4*tem + tanh(k5*(tem-218.8))*(k6 + k7/tem + k8*log(tem) + k9*tem)
    cal_esat=exp(es)
    
    return
  end function
  
  function cal_esati (tem)
    real :: tem, es, cal_esati
    real, parameter :: k1=9.550426, k2=-5723.265, k3=3.53068, k4=-0.00728332
    
    es=k1 + k2/tem + k3*log(tem) + k4*tem
    cal_esati=exp(es)
    
    return
  end function
  
  function cal_qsat (tem,pp)
    real :: pp, tem, cal_qsat
    
    cal_qsat=0.622*cal_esat(tem) / (pp - cal_esat(tem))
    
    return
  end function
  
  function cal_qsati (tem,pp)
    real :: pp, tem, cal_qsati
    
    cal_qsati=0.622*cal_esati(tem) / (pp - cal_esati(tem))
    
    return
  end function
  
  function cal_dqsatdt (tem,pp)
    real :: pp, tem, desdt, cal_dqsatdt
    real, parameter :: k1=54.842763, k2=-6763.22, k3=-4.21, k4=0.000367
    real, parameter :: k5=0.0415, k6=53.878, k7=-1331.22, k8=-9.44523, k9=0.014025
    
    desdt=-k2/(tem*tem) + k3/tem + k4 + k5*(1.-tanh(k5*(tem-218.8))*tanh(k5*(tem-218.8)))* &
        (k6 + k7/tem + k8*log(tem) + k9*tem) + tanh(k5*(tem-218.8))*(-k7/(tem*tem) + k8/tem + k9)
    desdt=desdt*cal_esat(tem)
    
    cal_dqsatdt=0.622*pp/(pp-cal_esat(tem))**2. * desdt

    return
  end function
  
  function cal_dqsatidt (tem,pp)
    real :: pp, tem, desdt, cal_dqsatidt
    real, parameter :: k1=9.550426, k2=-5723.265, k3=3.53068, k4=-0.00728332
    
    desdt=-k2/(tem*tem) + k3/tem + k4
    desdt=desdt*cal_esati(tem)
    
    cal_dqsatidt=0.622*pp/(pp-cal_esati(tem))**2. * desdt
    
    return
  end function
  
  function cal_flv (tem)
    real :: tem, cal_flv
    real, parameter :: T0=273.15
    
    cal_flv=(2500.8-2.36*(tem-T0)+0.0016*(tem-T0)**2.-0.00006*(tem-T0)**3.)*1000.
!     cal_flv=2501000.
    
    return
  end function
  
  function cal_fls (tem)
    real :: tem, cal_fls
    real, parameter :: T0=273.15
    
    cal_fls=(2834.1-0.29*(tem-T0)-0.004*(tem-T0)**2.)*1000.
    
    return
  end function
  
  function cal_flm (tem)
    real :: tem, cal_flm
    real :: flv,fls
    
    flv=cal_flv(tem)
    fls=cal_fls(tem)
    cal_flm=fls-flv
    
    return
  end function
	   
end module equation_of_state
