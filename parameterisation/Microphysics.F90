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

module microphysics
  !!< This module implements a warm cloud microphysics model in Fluidity

  use spud
  use state_module
  use fields
  use field_derivatives
  use sparse_tools
  use fefields
  use fetools
  use boundary_conditions
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN,&
       PYTHON_FUNC_LEN, CURRENT_DEBUG_LEVEL
  use solvers
  use python_state
  use sparsity_patterns_meshes
  use diagnostic_fields, only: safe_set
  use equation_of_state
  use filter_diagnostics
  use slope_limiters_dg

  implicit none
  
  real, parameter :: xmin=1.e-15, xnmin=0.1, xactiv=4.1888e-15, rho0=1.225, rhow=1000., kt=2.4e-2, pi=3.1415927
  real, parameter :: N0r = 1.e7, N0g = 4.e6, N0s = 2.e7
  real :: rauto, mauto, k0, rd0, sd0, mindt, rtol, relaxation
  character(len=OPTION_PATH_LEN) :: condensation_evaporation
  logical :: constant_cp_cv, sat_adj=.false., drop_adj=.false.

  real, dimension(5), parameter  :: nu=(/5., 0., 0., 0., 0./),			 &
  				    nu_sb=(/1., -2./3., 1., 1., 1./),		 &
				    mu_sb=(/1., 1./3., 1./3., 1./3., 1./3./),	 &
  				    am=(/3., 3., 3., 3., 3./),			 &
  				    cm=(/523.6, 523.6, 523.6, 523.6, 523.6/),	 &
				    btv=(/2., 0.9, 0.9, 0.9, 0.9/),		 &
				    ctv=(/3.025e+7, 1860., 1860., 1860., 1860./),&	
				    al_sb=(/0., 159., 0., 0., 0./),		 &
				    be_sb=(/0., 0.266, 0.667, 0.667, 0.667/)	
  type microphysics_field
    type(scalar_field), dimension(3) :: data
    type(scalar_field), pointer :: forcing, sinking_velocity
    logical :: has_forcing,has_sinking_velocity
  end type microphysics_field

  type logic_array
     logical :: have_temp=.false.,&
         have_pottemp=.false.,&
         have_qv=.false.,&
         have_EOSDensity=.false.,&
         have_EOSPressure=.false.
  end type logic_array

  private
  public :: calculate_microphysics_forcings, initialise_microphysics, store_microphysics_source, &
       calculate_diagnostic_microphysics, calculate_lambda_two, calculate_lambda_one, limit_microphysics
  public :: nu, am, cm, btv, ctv, cal_gamma 

contains

  subroutine store_microphysics_source(state,name)
    type(state_type),intent(inout) :: state
    character(len=*) :: name
    character(len=FIELD_NAME_LEN) :: parent
    type(vector_field), pointer :: X
    type(scalar_field), pointer :: sfield, sfield1, sfield2
    type(scalar_field) :: sponge
    integer :: stat, ind

    sfield1=>extract_scalar_field(state,name,stat=stat)
    if (stat /= 0) then
       ewrite(3,*) ('No field '//name//' in state to store!')
       return
    end if
    
    ind = INDEX(name,'MicrophysicsSource')
    parent = name(1:ind-1)
    sfield => extract_scalar_field(state, trim(parent), stat=stat)
    if (stat /= 0) &
        call get_thermo_variable(state,sfield)
    
    if ( have_option(trim(sfield%option_path)//&
    	 "/prognostic/scalar_field::Absorption/diagnostic/algorithm::&
    	  atmosphere_forcing_scalar/sponge_layer_scalar_absorption") ) then
       X => extract_vector_field(state, "Coordinate")
       call allocate(sponge, sfield1%mesh, "SpongeLayer")
       call calculate_sponge_coefficient_scalar (sponge, X, sfield%option_path)

       call scale(sfield1, sponge)
       call deallocate(sponge)
    endif
    
    sfield2=>extract_scalar_field(state,'Old'//name,stat=stat)
    if (stat /= 0) then
       ewrite(3,*) ('No old field '//name//' in state to store in!')
       return
    end if

    call set(sfield2,sfield1)

  end subroutine store_microphysics_source

  subroutine reset_microphysics_source(state,name)
    type(state_type),intent(inout) :: state
    character(len=*) :: name
    type(scalar_field), pointer :: sfield
    integer :: stat

    sfield=>extract_scalar_field(state,name,stat=stat)
    if (stat /= 0) then
       ewrite(3,*) ('No field '//name//' in state to store!')
       return
    end if

    call zero(sfield)
    
  end subroutine reset_microphysics_source

  subroutine calculate_microphysics_forcings(state,current_time,dt,adapt,reset)
    type(state_type), intent(inout) :: state(:)
    type(scalar_field), pointer :: sfield1,sfield2
    type(mesh_type), pointer :: mesh
    real, intent(in) :: current_time
    real, intent(in) :: dt
    logical, intent(in), optional :: adapt,reset

    integer :: i, j, nmom, stat
    logical :: have_cold, have_ccn, have_ndrop, have_nrain
    character(len=OPTION_PATH_LEN) :: mom_path, micro_path, micro_name
    type(scalar_field), pointer :: Qdrop
    type(logic_array) :: logic

    ewrite(1,*) 'Entering calculate_microphysics_forcings'
    
    do j=1,size(state)	!<------------------------------
    
    micro_path="/material_phase["//int2str(j-1)//"]/cloud_microphysics/fortran_microphysics"
    if (have_option(trim(micro_path))) then
       if (.not. present(adapt)) then
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
	   
	  Qdrop=>extract_scalar_field(state(j), 'Qdrop')
          mesh=>extract_mesh(state(j), trim(Qdrop%mesh%name), stat=stat) 
	  
          if (stat /= 0) then
             ewrite(-1,*) "Mesh:"
             ewrite(-1,*) trim(mesh%name)
             FLExit("does not exit. Stopping!")
          endif

          if (size(state)>1) then
             if (have_option(trim(state(j)%option_path)//'/equation_of_state/compressible/ATHAM')) &
             call store_microphysics_source(state(j),"MassFractionMicrophysicsSource")
             call reset_microphysics_source(state(j),"MassFractionMicrophysicsSource")
          else
             call store_microphysics_source(state(j),"QdropMicrophysicsSource")
             call store_microphysics_source(state(j),"QrainMicrophysicsSource")
             if (have_ndrop) call store_microphysics_source(state(j),"NdropMicrophysicsSource")
             if (have_nrain) call store_microphysics_source(state(j),"NrainMicrophysicsSource")	 
             if (have_ccn) call store_microphysics_source(state(j),"CCNMicrophysicsSource")	 

             call reset_microphysics_source(state(j),"QdropMicrophysicsSource")
             call reset_microphysics_source(state(j),"QrainMicrophysicsSource")
             if (have_ndrop) call reset_microphysics_source(state(j),"NdropMicrophysicsSource")
             if (have_nrain) call reset_microphysics_source(state(j),"NrainMicrophysicsSource")	     
	     if (have_ccn) call reset_microphysics_source(state(j),"CCNMicrophysicsSource")	     

	     if (have_cold) then
               call store_microphysics_source(state(j),"QiceMicrophysicsSource")
               call store_microphysics_source(state(j),"QgrauMicrophysicsSource")
               call store_microphysics_source(state(j),"QsnowMicrophysicsSource")
               call store_microphysics_source(state(j),"NiceMicrophysicsSource")
               if (nmom>1) call store_microphysics_source(state(j),"NgrauMicrophysicsSource")	 
               if (nmom>1) call store_microphysics_source(state(j),"NsnowMicrophysicsSource")	 
	       
               call reset_microphysics_source(state(j),"QiceMicrophysicsSource")
               call reset_microphysics_source(state(j),"QgrauMicrophysicsSource")
               call reset_microphysics_source(state(j),"QsnowMicrophysicsSource")
               call reset_microphysics_source(state(j),"NiceMicrophysicsSource")
               if (nmom>1) call reset_microphysics_source(state(j),"NgrauMicrophysicsSource")	 
               if (nmom>1) call reset_microphysics_source(state(j),"NsnowMicrophysicsSource")	 
	     endif    	     
          end if
	  
          i=4
          if (has_scalar_field(state(j),"Temperature")) then
             call store_microphysics_source(state(j),"TemperatureMicrophysicsSource")
             call reset_microphysics_source(state(j),"TemperatureMicrophysicsSource")
             logic%have_temp=.true.
             i=i-1
          end if
          if (has_scalar_field(state(j),"EOSPressure")) then
             call store_microphysics_source(state(j),"EOSPressure")
             call reset_microphysics_source(state(j),"EOSPressure")
             logic%have_EOSPressure=.true.
             i=i-1
          end if
          if (has_scalar_field(state(j),"EOSDensity")) then
             call store_microphysics_source(state(j),"EOSDensity")
             call reset_microphysics_source(state(j),"EOSDensity")
             logic%have_EOSDensity=.true.
             i=i-1
          end if
          if (has_scalar_field(state(j),"PotentialTemperature")) then
             call store_microphysics_source(state(j),"PotentialTemperatureMicrophysicsSource")
             call reset_microphysics_source(state(j),"PotentialTemperatureMicrophysicsSource")
             logic%have_pottemp=.true.
             i=i-1
          end if
          if (has_scalar_field(state(j),"ConservedPotentialTemperature")) then
             call store_microphysics_source(state(j),"ConservedPotentialTemperatureMicrophysicsSource")
             call reset_microphysics_source(state(j),"ConservedPotentialTemperatureMicrophysicsSource")
             logic%have_pottemp=.true.
             i=i-1
          end if
          if (has_scalar_field(state(j),"VapourWaterQ")) then
             call store_microphysics_source(state(j),"VapourWaterQMicrophysicsSource")
             call reset_microphysics_source(state(j),"VapourWaterQMicrophysicsSource")
             logic%have_qv=.true.
             i=i-1
          end if
       end if
       
       if (present(reset)) exit
 
       ! Call main microphysics
       if (.not.have_option(trim(micro_path)//'/no_sources')) 	&
           call calculate_microphysics_from_fortran(state(j),mesh,micro_name,mom_path,current_time,dt)

       if (present(adapt)) then
          if (size(state)>1) then
             call store_microphysics_source(state(j),"MassFractionMicrophysicsSource")
          else
             call store_microphysics_source(state(j),"QdropMicrophysicsSource")
             call store_microphysics_source(state(j),"QrainMicrophysicsSource")
             if (have_ndrop) call store_microphysics_source(state(j),"NdropMicrophysicsSource")
             if (have_nrain) call store_microphysics_source(state(j),"NrainMicrophysicsSource")
             if (have_ccn) call store_microphysics_source(state(j),"CCNMicrophysicsSource")
	     
	     if (have_cold) then
               call store_microphysics_source(state(j),"QiceMicrophysicsSource")
               call store_microphysics_source(state(j),"QgrauMicrophysicsSource")
               call store_microphysics_source(state(j),"QsnowMicrophysicsSource")
               call store_microphysics_source(state(j),"NiceMicrophysicsSource")
               if (nmom>1) call store_microphysics_source(state(j),"NgrauMicrophysicsSource")
               if (nmom>1) call store_microphysics_source(state(j),"NsnowMicrophysicsSource")
             endif
	  end if
	  
          if (has_scalar_field(state(j),"Temperature")) &
              call store_microphysics_source(state(j),"TemperatureMicrophysicsSource")
          if (has_scalar_field(state(j),"PotentialTemperature")) &
              call store_microphysics_source(state(j),"PotentialTemperatureMicrophysicsSource")
          if (has_scalar_field(state(j),"ConservedPotentialTemperature")) &
              call store_microphysics_source(state(j),"ConservedPotentialTemperatureMicrophysicsSource")
          if (has_scalar_field(state(j),"VapourWaterQ")) &
              call store_microphysics_source(state(j),"VapourWaterQMicrophysicsSource")
       end if
    end if
       
    enddo		!<-----------------------

    ewrite(1,*) 'End calculate_microphysics_forcings'

  end subroutine calculate_microphysics_forcings

  subroutine calculate_microphysics_from_fortran(state,mesh,mname,mom_path,current_time,dt)
    ! Set microphysical  source terms from python.
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout), pointer :: mesh
    real, intent(in) :: current_time
    real, intent(in) :: dt
    character(len=FIELD_NAME_LEN):: option_path, mom_path, mname

    logical :: have_qt=.false.
    integer :: stat,i,j,ele,thermal_variable
    real :: qs, dqsdt, cp, cv, ccn, g, w, gsat
    real :: c_p,c_v,c_p_v,c_v_v,c_p_l,c_v_l,c_p_i,c_v_i

    type(microphysics_field) :: ptem,pth,pq_t,pq_v,pq_r,pq_c,pq_i,pq_g,pq_s,pn_r,pn_c,pn_i,pn_g,pn_s,pn_ccn,pp,prho
    type(scalar_field) :: qsat, dqsatdt, qc_p, qc_v
    type(scalar_field) :: temperature, eosdensity, supersat, pressure_remap
    type(vector_field) :: grad_saturation, pvelocity, pposition
    type(scalar_field), pointer :: ccn_concentration, pressure, thermal
    type(scalar_field), pointer :: tmp
    type(vector_field), pointer :: velocity, position

    real, dimension(5) :: tem,th,pt,q_t,q_v,q_r,q_c,q_i,q_g,q_s,n_r,n_c,n_i,n_g,n_s,n_ccn,p,rho

    ! data arrays for the pass into the fortran routine
    ! These have the following order
    ! (1) Old Timelevel Value (input, projected)
    ! (2) Previous nonlinear iteration (input, projected)
    ! (3) Current nonlinear iteration (input, projected)
    ! (4) Sinking Velocity (output, on Microphysics mesh)
    ! (5) Microphysics forcing (output, on Microphysics mesh)
    
    ewrite(1,*) 'In calculate_microphysics_from_fortran ', trim(mname)

    ! Extract microphysics and thermo constants
    call extract_constants (state)

    sat_adj=have_option('/material_phase[0]/cloud_microphysics/saturation_adjustment')
    if (.not.sat_adj) then
      call get_option('/material_phase[0]/cloud_microphysics/condensation_evaporation/name',condensation_evaporation)
      if (trim(condensation_evaporation)=='Adaptive') then
        call get_option('/material_phase[0]/cloud_microphysics/condensation_evaporation::Adaptive/min_step_size',mindt)
        call get_option('/material_phase[0]/cloud_microphysics/condensation_evaporation::Adaptive/tolerance',rtol)
      else
        drop_adj = have_option('/material_phase[0]/cloud_microphysics/fortran_microphysics/two_moment_microphysics::'//trim(mname)//'/droplet_adjustment')
      endif      
    endif
    constant_cp_cv=have_option('/material_phase[0]/equation_of_state/compressible/giraldo/constant_cp_cv')
    
    call get_option('/material_phase[0]/cloud_microphysics/relaxation',relaxation,default=1.)

    call get_option("/physical_parameters/gravity/magnitude", g, stat)
    have_qt = has_scalar_field(state,"TotalWaterQ")

    call allocate (temperature,mesh,"Temperature")
    call allocate (pressure_remap,mesh,"PressureRemap")
    call allocate (eosdensity,mesh,"EOSDensity")
    call allocate (qsat,mesh,"QSaturation")
    call allocate (dqsatdt,mesh,"dQSaturationdT")
    call allocate (supersat,mesh,"Supersaturation")
    call allocate (qc_p,mesh,"CpLocal")
    call allocate (qc_v,mesh,"CvLocal")
    
    call get_thermo_variable(state, thermal, index=thermal_variable)
    assert(mesh==thermal%mesh)
    
    call get_cp_cv(c_p,c_v,c_p_v=c_p_v,c_v_v=c_v_v,c_p_l=c_p_l,c_v_l=c_v_l,c_p_i=c_p_i,c_v_i=c_v_i)
    
    call compressible_eos(state,saturation=qsat,supersaturation=supersat,dqsaturation=dqsatdt,density=eosdensity,&
                          temperature=temperature,qc_p=qc_p,qc_v=qc_v)
    
    pressure=>extract_scalar_field(state,"Pressure")
    call safe_set(state,pressure_remap,pressure)
    
    call extract_and_project(state,mesh,pp,   	"PressureRemap", local=pressure_remap)
    call extract_and_project(state,mesh,ptem,   "Temperature", local=temperature)  
    call extract_and_project(state,mesh,prho, 	"EOSDensity", local=eosdensity)
    call extract_and_project(state,mesh,pth,     trim(thermal%name), local=thermal)  
    
    if (have_qt) then 
      call extract_and_project(state,mesh,pq_t,"TotalWaterQ")    
    else
      call extract_and_project(state,mesh,pq_v,"VapourWaterQ")    
    endif
    call extract_and_project(state,mesh,pq_c,       "Qdrop")
    call extract_and_project(state,mesh,pq_r,       "Qrain")
    call extract_and_project(state,mesh,pq_i,       "Qice")
    call extract_and_project(state,mesh,pq_g,       "Qgrau")
    call extract_and_project(state,mesh,pq_s,       "Qsnow")
    call extract_and_project(state,mesh,pn_c,       "Ndrop")
    call extract_and_project(state,mesh,pn_r,       "Nrain")
    call extract_and_project(state,mesh,pn_i,       "Nice")
    call extract_and_project(state,mesh,pn_g,       "Ngrau")
    call extract_and_project(state,mesh,pn_s,       "Nsnow")
    call extract_and_project(state,mesh,pn_ccn,     "CCN")
    
    position => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    call allocate(pvelocity,position%dim,mesh,"ProjVelocity")
    call allocate(grad_saturation,position%dim,mesh,"GradientSupersaturation")
    call allocate(pposition,position%dim,mesh,"ProjPosition")

    call remap_field(position, pposition)    
    if (trim(velocity%mesh%name) /= trim(mesh%name)) then
      call remap_field(velocity, pvelocity)
    else
      call set(pvelocity, velocity)
    endif
   
    call grad(supersat,pposition,grad_saturation)
    
    do ele=1,node_count(pp%data(1))

       ! Get element arrays
       call set_local_array(q_c,pq_c,ele)
       call set_local_array(q_r,pq_r,ele)
       call set_local_array(q_i,pq_i,ele)
       call set_local_array(q_g,pq_g,ele)
       call set_local_array(q_s,pq_s,ele)
       call set_local_array(n_c,pn_c,ele)
       call set_local_array(n_r,pn_r,ele)
       call set_local_array(n_i,pn_i,ele)
       call set_local_array(n_g,pn_g,ele)
       call set_local_array(n_s,pn_s,ele)
       call set_local_array(p,pp,ele)
       call set_local_array(rho,prho,ele)
       call set_local_array(tem,ptem,ele)
       call set_local_array(th,pth,ele)
       call set_local_array(n_ccn,pn_ccn,ele)
       if (have_qt) then
         call set_local_array(q_t,pq_t,ele)
         q_v=q_t-q_c-q_r-q_i-q_s-q_g
       else
         call set_local_array(q_v,pq_v,ele)
         q_t=q_v+q_c+q_r+q_i+q_s+q_g
       endif
       
       w=node_val(pvelocity,position%dim,ele)
       cp=node_val(qc_p,ele)
       cv=node_val(qc_v,ele)
       qs=node_val(qsat,ele)
       gsat=node_val(grad_saturation,position%dim,ele)

       !   Call external routines
       select case (trim(mname))
       
         case('morrison', 'seifert_beheng')
         call microphysics_2mom(state,sat_adj,mom_path,mname,current_time,dt,	&
         	   w,qs,gsat,g,cp,cv,q_t,q_v,q_r,q_c,q_i,q_g,q_s,n_r,n_c,n_i,n_g,n_s,n_ccn,tem,p,rho)

         case('thompson')
	 call microphysics_1mom(state,sat_adj,current_time,dt,		&
         	   qs,cp,q_t,q_v,q_r,q_c,q_i,q_g,q_s,n_r,n_c,n_i,n_g,n_s,tem,p,rho)
       
       end select
       
       ! Call store results
       call store_result(pq_c,q_c,ele)
       call store_result(pq_r,q_r,ele)
       call store_result(pq_i,q_i,ele)
       call store_result(pq_g,q_g,ele)
       call store_result(pq_s,q_s,ele)
       call store_result(pn_c,n_c,ele)
       call store_result(pn_r,n_r,ele)
       call store_result(pn_i,n_i,ele)
       call store_result(pn_g,n_g,ele)
       call store_result(pn_s,n_s,ele)
       call store_result(pn_ccn,n_ccn,ele)
       if (have_qt) then
         call store_result(pq_t,q_t,ele)
       else
         call store_result(pq_v,q_v,ele)
       endif
       call store_thermal(thermal_variable,pth,th,tem,p,q_v,cp,cv,ele)
    end do
    
    if (have_qt) then
      call clean_up(pq_t)
    else
      call clean_up(pq_v)
    endif
    call clean_up(pq_c)
    call clean_up(pq_r)
    call clean_up(pq_i)
    call clean_up(pq_g)
    call clean_up(pq_s)
    call clean_up(pn_c)
    call clean_up(pn_r)
    call clean_up(pn_i)
    call clean_up(pn_g)
    call clean_up(pn_s)
    call clean_up(pp)
    call clean_up(prho)
    call clean_up(ptem)
    call clean_up(pth)
    call clean_up(pn_ccn)
    
    call deallocate (temperature)
    call deallocate (pressure_remap)
    call deallocate (eosdensity)
    call deallocate (qsat)
    call deallocate (dqsatdt)
    call deallocate (supersat)
    call deallocate (qc_p)
    call deallocate (qc_v)
    call deallocate (grad_saturation)
    call deallocate (pvelocity)
    call deallocate (pposition)
    
    ewrite(1,*) 'End calculate_microphysics_from_fortran'

contains 

    subroutine clean_up(field)
      type(microphysics_field), intent(inout) :: field
      integer :: i

      do i=1,size(field%data)
         call deallocate(field%data(i))
      end do

      nullify(field%forcing)
      nullify(field%sinking_velocity)
         
    end subroutine clean_up

    subroutine set_local_array(local_array,fields,n)
      real, dimension(:), intent(inout) :: local_array
      type(microphysics_field), intent(in) :: fields
      integer, intent(in) :: n
      integer :: i
      
      local_array=0.0
      do i=1,size(fields%data)
         local_array(i)=node_val(fields%data(i),n)
      end do
      
    end subroutine set_local_array
    
    subroutine store_result(field,vloc,n)
      type(microphysics_field), intent(inout) :: field
      real, intent(in), dimension(5) :: vloc
      integer, intent(in) :: n

      if (field%has_sinking_velocity)&
           call set(field%sinking_velocity,n,vloc(4))
      if (field%has_forcing)&
           call set(field%forcing,n,vloc(5))

    end subroutine store_result
    
    subroutine store_thermal(index,field,vloc,tem,p,q_v,cp,cv,n)
      type(microphysics_field), intent(inout) :: field
      real, intent(in), dimension(5) :: vloc
      real, intent(in), dimension(5) :: p, tem, q_v
      integer, intent(in) :: n, index
      real :: source, precip, cp, cv, pp
      
      if (index==1) then
        source = cv*tem(5)
      else if (index==2) then
        source = tem(5)
      else if (index==3) then
        pp = (p(1)/1.E+05)**((cp-cv)/cp)
        source = tem(5) * vloc(3)/tem(3)
        if (.not.constant_cp_cv) then
          source = source - vloc(3)*log(pp)*((c_p_v-c_v_v)/(cp-cv) - (c_p_v-c_p_l)/cp)*q_v(5)
        endif
      endif

      if (field%has_forcing)&
           call set(field%forcing,n,source)
      
    end subroutine store_thermal
        
    subroutine extract_constants (state)
      implicit none
      type(state_type), intent(inout) :: state
      character(len=OPTION_PATH_LEN) :: option_path,eos_path,mic_path
      character(len=FIELD_NAME_LEN) :: mic_name
      integer :: stats(4)

      !
      !  Set thermodynamic constants
      !
      option_path='/material_phase::'//trim(state%name)
      eos_path = trim(option_path)//'/equation_of_state'

      !
      !  Set microphysical constants
      !
      mic_path =trim(option_path)//'/cloud_microphysics'
      if(have_option(trim(mic_path)//'/fortran_microphysics/two_moment_microphysics')) then
        call get_option(trim(mic_path)//'/fortran_microphysics/two_moment_microphysics/name',mic_name)
        mic_path = trim(mic_path)//'/fortran_microphysics/two_moment_microphysics::'//trim(mic_name)
	  
	call get_option(trim(mic_path)//'/autoconversion_radius', &
        		rauto, stat=stats(1))
        
	if (have_option(trim(mic_path)//'/detailed_activation')) then
          call get_option(trim(mic_path)//'/detailed_activation/CCN_kappa',         &
                          k0, stat=stats(2))
          call get_option(trim(mic_path)//'/detailed_activation/CCN_mean_radius',   &
                          rd0, stat=stats(3))
          call get_option(trim(mic_path)//'/detailed_activation/CCN_dev_radius',    &
                          sd0, stat=stats(4))

          if(stats(1)/=0) then
            rauto = 40.e-6
          end if
          if(stats(2)/=0) then
            k0 = 0.7
          endif    
          if(stats(3)/=0) then
            rd0 = 0.3e-6
          endif    
          if(stats(4)/=0) then
            sd0 = 1.4
          endif    
        endif
	
      else if (have_option(trim(mic_path)//'/fortran_microphysics/one_moment_microphysics')) then
        call get_option(trim(mic_path)//'/fortran_microphysics/one_moment_microphysics/name',mic_name)
        mic_path = trim(mic_path)//'/fortran_microphysics/one_moment_microphysics::'//trim(mic_name)
        call get_option(trim(mic_path)//'/mass_threshold', &
                        mauto, stat=stats(1))

        if(stats(1)/=0) then
          mauto = 1.e-3
        end if
      endif
  
    end subroutine
 
  end subroutine calculate_microphysics_from_fortran

subroutine microphysics_2mom(state,sat_adj,mom_path,mname,time,dt, &
        w,qs,gsat,g,cp,cv,lq_t,lq_v,lq_r,lq_c,lq_i,lq_g,lq_s, &
	ln_r,ln_c,ln_i,ln_g,ln_s,ln_ccn,ltem,lp,lrho)

  implicit none
  type(state_type),intent(inout) :: state
  logical :: sat_adj, have_ndrop, have_nrain, have_ccn
  real :: time, dt, qs, cp, cv, g, w, gsat
  real :: lq_c0, lq_r0, ln_c0, ln_r0, lq_v0, ltem0, ln_ccn0, lp0, lrho0
  real :: lq_c1, lq_r1, ln_c1, ln_r1, lq_v1, ltem1, ln_ccn1
  real, dimension(5) :: lq_t,lq_v,lq_r,lq_c,lq_i,lq_g,lq_s,ln_r,ln_c,ln_i,ln_g,ln_s,ln_ccn,ltem,lp,lrho
  character(len=FIELD_NAME_LEN):: mom_path, mname
  
  have_ndrop = have_option(trim(mom_path)//'/scalar_field::Ndrop/prognostic')
  have_nrain = have_option(trim(mom_path)//'/scalar_field::Nrain/prognostic')
  have_ccn = have_option(trim(mom_path)//'/scalar_field::CCN/prognostic')

  ! data arrays for the pass into the fortran routine
  ! These have the following order
  
  ! (1) Old Timelevel Value (input, projected)
  ! (2) Previous nonlinear iteration (input, projected)
  ! (3) Current nonlinear iteration (input, projected)
  ! (4) Sinking Velocity (output, on Microphysics mesh)
  ! (5) Microphysics forcing (output, on Microphysics mesh)
    
  lq_v(5) = 0.
  lq_c(5) = 0.
  lq_r(5) = 0.
  ln_c(5) = 0.
  ln_r(5) = 0.
  ltem(5) = 0.
  ln_ccn(5) = 0.
  
  lp0   = lp(3)
  lrho0 = lrho(3)
  lq_v0 = lq_v(3); lq_v1 = lq_v(3)
  lq_c0 = lq_c(3); lq_c1 = lq_c(3)
  lq_r0 = lq_r(3); lq_r1 = lq_r(3)
  ln_c0 = ln_c(3); ln_c1 = ln_c(3)
  ln_r0 = ln_r(3); ln_r1 = ln_r(3)
  ltem0 = ltem(3); ltem1 = ltem(3)
  ln_ccn0 = ln_ccn(3); ln_ccn1 = ln_ccn(3)

  !
  !  1.: Evaluate microphysical processes: diffusion/evaporation
  !
  select case(trim(condensation_evaporation))
  
    case('Analytic')
      call diff_analytic (dt, w, qs, gsat, g, cp, lq_v, lq_c, ln_c, lq_r, ln_r, ln_ccn, lrho, lp, ltem)
      
    case('Adaptive')
      call diff_adaptive_rkf (dt, w, qs, g, cp, lq_v, lq_c, ln_c, lq_r, ln_r, ln_ccn, lrho, lp, ltem)
  
    case default
      call diff_simple (dt, qs, g, cp, lq_v, lq_c, ln_c, lq_r, ln_r, ln_ccn, lrho, lp, ltem)
    
  end select
  
  lq_v1 = lq_v1 + dt*lq_v(5)
  lq_c1 = lq_c1 + dt*lq_c(5)
  lq_r1 = lq_r1 + dt*lq_r(5)
  ln_c1 = ln_c1 + dt*ln_c(5)
  ln_r1 = ln_r1 + dt*ln_r(5)
  ltem1 = ltem1 + dt*ltem(5)
  ln_ccn1 = ln_ccn1 + dt*ln_ccn(5)
  
  lq_v(5) = 0.
  lq_c(5) = 0.
  lq_r(5) = 0.
  ln_c(5) = 0.
  ln_r(5) = 0.
  ltem(5) = 0.
  ln_ccn(5) = 0.

  !
  !  2.: Evaluate microphysical processes: collision-coalescence
  !
  if (trim(mname) == 'morrison') then 
    ! Calculate autoconversion and self-collection rates: c -> r & c+c -> c or r
    call aucr_k (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  
    ! Calculate self-collection rates: r+r -> r
    if (have_nrain) call scr_k (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  
    ! Calculate accretion rates: c+r -> r
    call ccr_k (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  
  else if (trim(mname) == 'seifert_beheng') then 
    ! Calculate autoconversion and self-collection rates: c -> r & c+c -> c or r
    call aucr_sb (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  
    ! Calculate self-collection rates: r+r -> r
    if (have_nrain) call scr_sb (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  
    ! Calculate accretion rates: c+r -> r
    call ccr_sb (dt, lrho, lq_c, lq_r, ln_c, ln_r)
    
  endif
  
  lq_c1 = lq_c1 + dt*lq_c(5)
  lq_r1 = lq_r1 + dt*lq_r(5)
  ln_c1 = ln_c1 + dt*ln_c(5)
  ln_r1 = ln_r1 + dt*ln_r(5)
  
  lq_c(5) = 0.
  lq_r(5) = 0.
  ln_c(5) = 0.
  ln_r(5) = 0.

  !
  !  3.: Activate droplets 
  !  
  call activ (dt, w, cp, cv, g, gsat, ln_ccn, lq_v, lq_c, ln_c, lp, ltem, mom_path)  
  
  lq_c1 = lq_c1 + dt*lq_c(5)
  ln_c1 = ln_c1 + dt*ln_c(5)
  lq_v1 = lq_v1 + dt*lq_v(5)
  ltem1 = ltem1 + dt*ltem(5)
  ln_ccn1 = ln_ccn1 + dt*ln_ccn(5)
  
  !
  !  4.: Sedimentation velocities
  !
  call sedim (dt, 2, lq_r0, ln_r0, lq_r1, ln_r1, lq_r(4), ln_r(4))
  
  !
  !  5.: Assemble microphysics sources 
  !  
  call assemble_micro (dt, cp, lq_v, lq_c, ln_c, lq_r, ln_r, ln_ccn, ltem)

contains

subroutine diff_analytic (dt, w, qs, gsat, g, cp, lqv, lqc, lnc, lqr, lnr, lnccn, lrho, lp, ltem)

  implicit none
  integer :: i
  real  :: dt, fl, fv, vis, dif, gamm, tau, fac, aa
  real  :: dql, ql, qc, nc, qr, nr, qv, ss
  real  :: w, qs, es, dqsdt, g, cp, tem, pp, den, sat0, sat, gsat
  real, dimension(5) :: lqv,lqr,lnr,lqc,lnc,lp,ltem,lrho,lnccn,cd

  pp   = lp0
  den  = lrho0
  tem  = relaxation*ltem1+(1.-relaxation)*ltem0
  qv   = relaxation*lq_v1+(1.-relaxation)*lq_v0
  qc   = relaxation*lq_c1+(1.-relaxation)*lq_c0
  nc   = relaxation*ln_c1+(1.-relaxation)*ln_c0
  qr   = relaxation*lq_r1+(1.-relaxation)*lq_r0
  nr   = relaxation*ln_r1+(1.-relaxation)*ln_r0

  call diff_cd (dt, qs, cp, qv, qc, nc, qr, nr, den, pp, tem, cd)

  ql  = qc + qr
  dql = cd(1) + cd(2)
  
  if ( dql > 1.e-6 ) then

    dqsdt= cal_dqsatdt(ltem1,pp)
    es   = cal_esat(ltem1)
    sat0 = lq_v1 - cal_qsat(ltem1,pp)
    ss   = lq_v1 / cal_qsat(ltem1,pp)
    
    aa   = g/cp*dqsdt - den*g*qs/(pp - es)
    gamm = 1. + ss*cal_flv(tem)/cp*dqsdt
    tau  = max(1./(gamm*dql),1.e-15)
    
    if (have_option('/material_phase[0]/cloud_microphysics/time_integration::Splitting') .or. &
        have_option('/material_phase[0]/cloud_microphysics/time_integration::Strang')) then
      sat = tau/dt*sat0*(1. - exp(-dt/tau))
    else
      aa  = g/cp*dqsdt - den*g*qs/(pp - es)
      sat = w*aa*tau + tau/dt*(sat0 - w*aa*tau - w*gsat*dt)*(1. - exp(-dt/tau))
    endif
    
    if ((lq_c1 > xmin .and. ln_c1 > xnmin) .and. (.not.drop_adj)) then
      lqc(5) = lqc(5) + max(min(cd(1)*sat, qv/dt), -lq_c1/dt)
      if (have_ndrop) lnc(5) = lnc(5) + min(max(cd(1)*sat*ln_c0/lq_c0,-ln_c1/dt),0.)
    endif
    
    if (lq_r1 > xmin .and. ln_r1 > xnmin) then    
      lqr(5) = lqr(5) + max(min(cd(2)*sat, qv/dt-lqc(5)), -lq_r1/dt)
      if (have_nrain) lnr(5) = lnr(5) + max(min(cd(2)*sat*ln_r0/lq_r0, 0.), -ln_r1/dt)
    endif
        
    lqv(5) = lqv(5) - lqc(5) - lqr(5)
    ltem(5) = ltem(5) + (lqc(5) + lqr(5))*cal_flv(tem)/cp
    
  endif
  
end subroutine
    
subroutine diff_simple (dt, qs, g, cp, lqv, lqc, lnc, lqr, lnr, lnccn, lrho, lp, ltem)
  implicit none
  integer :: i, nit
  real  :: dt, gamm, dqsdt, beta, dql
  real  :: qc0, qr0, qv0, tem0, nc0, nr0
  real  :: qs, ss, sat0, dqsatt, g, cp, pp, den, dqc, dqr
  real, dimension(5) :: lqv,lqr,lnr,lqc,lnc,lp,ltem,lrho,lnccn,cd

  beta  = 1.0

  pp    = lp0
  den   = lrho0
  tem0  = relaxation*ltem1+(1.-relaxation)*ltem0
  qv0   = relaxation*lq_v1+(1.-relaxation)*lq_v0
  qc0   = relaxation*lq_c1+(1.-relaxation)*lq_c0
  qr0   = relaxation*lq_r1+(1.-relaxation)*lq_r0
  nc0   = relaxation*ln_c1+(1.-relaxation)*ln_c0
  nr0   = relaxation*ln_r1+(1.-relaxation)*ln_r0

  call diff_cd (dt, qs, cp, qv0, qc0, nc0, qr0, nr0, den, pp, tem0, cd)
  
  dql = cd(1) + cd(2)
  
  if (dql > 1.e-6) then

    sat0   = lq_v1 - cal_qsat(ltem1,pp)
    ss     = lq_v1 / cal_qsat(ltem1,pp)
    dqsatt = cal_dqsatdt(ltem1,pp)
    gamm   = dql*(1. + ss*dqsatt*cal_flv(ltem0)/cp)
    
    if ((lq_c1 > xmin .and. ln_c1 > xnmin) .and. (.not.drop_adj)) then    
      dqc  = max(min(cd(1)*sat0 / (1. + beta*dt*gamm),qv0/dt),-lq_c1/dt)
      lqc(5) = lqc(5) + dqc
      if (have_ndrop) lnc(5)   = lnc(5) + max(min(ln_c0/lq_c0*dqc,0.),-ln_c1/dt)
!      if (have_ccn) lnccn(5) = lnccn(5) - min(ln_c0/lq_c0*(qc - lq_c0)/dt,0.)
    endif
    
    if (lq_r1 > xmin .and. ln_r1 > xnmin) then    
      dqr  = max(min(cd(2)*sat0 / (1. + beta*dt*gamm),max(qv0/dt-dqc,0.)),-lq_r1/dt)
      lqr(5) = lqr(5) + dqr
      if (have_nrain) lnr(5) = lnr(5) + max(min(ln_r0/lq_r0*dqr,0.),-ln_r1/dt)
    endif
   
    lqv(5) = lqv(5) - (lqc(5) + lqr(5))
    ltem(5) = ltem(5) - lqv(5)*cal_flv(tem0)/cp
  endif

end subroutine
    
subroutine diff_adaptive_rkf (dt, w, qs, g, cp, lqv, lqc, lnc, lqr, lnr, lnccn, lrho, lp, ltem)
  implicit none
  integer :: i, nit
  real  :: dt, gamm, dqsdt, beta, dql
  real  :: qc0, qr0, qv0, tem0, nc0, nr0
  real  :: w, qs, es, g, cp, tem, pp, den, qv, qc, qr
  real, dimension(5) :: lqv,lqr,lnr,lqc,lnc,lp,ltem,lrho,lnccn,cd

  pp    = lp0
  den   = lrho0
  tem0  = relaxation*ltem1+(1.-relaxation)*ltem0; tem = ltem0
  qv0   = relaxation*lq_v1+(1.-relaxation)*lq_v0;  qv = lq_v0
  qc0   = relaxation*lq_c1+(1.-relaxation)*lq_c0;  qc = lq_c0
  qr0   = relaxation*lq_r1+(1.-relaxation)*lq_r0;  qr = lq_r0
  nc0   = relaxation*ln_c1+(1.-relaxation)*ln_c0
  nr0   = relaxation*ln_r1+(1.-relaxation)*ln_r0

  call diff_cd (dt, qs, cp, qv0, qc0, nc0, qr0, nr0, den, pp, tem0, cd)
  
  dql = cd(1)+cd(2)
  
  if (dql > 1.e-6) then

    call direct_solve_rkf (dt, cd, cp, pp, tem0, qv0, qc0, qr0, tem, qv, qc, qr)

    if ((lq_c0 > xmin .and. ln_c0 > xnmin) .and. (.not.drop_adj)) then    
      lqc(5) = lqc(5) + (qc - lq_c0)/dt
      if (have_ndrop) lnc(5)   = lnc(5) + min(ln_c0/lq_c0*(qc - lq_c0)/dt,0.)
!      if (have_ccn) lnccn(5) = lnccn(5) - min(ln_c0/lq_c0*(qc - lq_c0)/dt,0.)
    endif
    
    if (lq_r0 > xmin .and. ln_r0 > xnmin) then    
      lqr(5) = lqr(5) + (qr - lq_r0)/dt
      if (have_nrain) lnr(5) = lnr(5) + min(ln_r0/lq_r0*(qr - lq_r0)/dt,0.)
    endif
    
    lqv(5) = lqv(5) - lqc(5) - lqr(5)
    ltem(5) = ltem(5) + (lqc(5) + lqr(5))*cal_flv(tem0)/cp
  endif

end subroutine

subroutine droplet_adjustment(qc, qr, qv, pp, tem)
  implicit none
  
  integer :: i
  real :: sat, pp, tem, qv, qc, qr
  real :: epsilon, cp, cv, dq_c, d_tem, qc_new
  real, parameter :: tol=1.e-8
  integer, parameter :: max_iter=15
    
  if ( qc > xmin ) then
	
    i=0
    epsilon=1.
    do while (abs(epsilon) > tol .and. i < max_iter)
      i=i+1

      cp=c_p
      cv=c_v
      if (.not.constant_cp_cv) then
  	cp = cp + qv*c_p_v + (qc+qr)*c_p_l
  	cv = cv + qv*c_v_v + (qc+qr)*c_v_l
      endif
      sat=cal_qsat(tem,pp)

      dq_c = (qv - sat) / (1. + cal_flv(tem)**2.*sat/(cp*(cp-cv) * tem**2.))
      qc_new=min(max(qc+dq_c,0.),qv)
      dq_c=qc_new-qc
    
      tem = tem + dq_c*cal_flv(tem)/cp
      qc  = qc + dq_c
      qv  = min(max(qv-dq_c,0.),qv+qc)
      
      epsilon = (qv-sat)/sat
  	
    enddo
    
  endif

end subroutine

subroutine aucr_sb (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  ! autoconversion of cloud droplets into rain drops
  ! plus droplet self-collection
  ! following Seifert and Beheng (2001)
  implicit none
  real :: xav, xr, xnc, xnr
  real :: aun, auq, scnc, aunr
  real :: xc, nuc, xprim, tau, fautau, dens
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,ln_r,ln_c,lrho  
  real, parameter :: kc=9.44e+09
  
  xprim = 1000.*4./3.*pi*rauto**3.
  
  dens = lrho0
  xc  = dens*(relaxation*lq_c1+(1.-relaxation)*lq_c0)
  xr  = dens*(relaxation*lq_r1+(1.-relaxation)*lq_r0)
  xnc = dens*(relaxation*ln_c1+(1.-relaxation)*ln_c0)
  xnr = dens*(relaxation*ln_r1+(1.-relaxation)*ln_r0)

  if (xnc > xnmin .and. xc > xmin) then
    xav = min(2.e-9,max(2.e-14,xc/xnc))
    tau = 1. - xc/(xc + xr)
    tau = min(max(tau,1.e-9),0.98)
    fautau = 600.*tau**0.68 * (1. - tau**0.68)**3.
  
    ! autoconversion rate (droplet conversion)
    auq = kc/(20.*xprim)*(nu_sb(1)+2.)*(nu_sb(1)+4.)/(nu_sb(1)+1.)**2.	&
        * xc**2.*xav**2.*(1. + fautau/(1.-tau)**2.) * rho0/dens
    
    ! autoconversion (rain formation)
    aunr = auq/xprim
    
    ! autoconversion + self-collection (droplet removal)
    scnc = kc*xc**2.*(nu_sb(1)+2.)/(nu_sb(1)+1.) * rho0/dens
    
    lq_c(5) = lq_c(5) - min(auq/dens,lq_c1/dt)
    lq_r(5) = lq_r(5) + min(auq/dens,lq_c1/dt)
    ln_c(5) = ln_c(5) - min(scnc/dens,ln_c1/dt)
    ln_r(5) = ln_r(5) + min(aunr/dens,ln_c1/dt)
  endif
  
end subroutine

subroutine scr_sb (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  ! Rain self-collection and break-up following Seifert and Beheng (2006)
  implicit none
  real :: lambda, avd
  real :: dens, xc, xr, xnc, xnr, xav
  real :: crrn, phibr
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,ln_r,ln_c,lrho
  real, parameter :: deq=1.3e-3, krr=7.12, kar=60.7
  
  dens = lrho0
  xc  = dens*(relaxation*lq_c1+(1.-relaxation)*lq_c0)
  xr  = dens*(relaxation*lq_r1+(1.-relaxation)*lq_r0)
  xnc = dens*(relaxation*ln_c1+(1.-relaxation)*ln_c0)
  xnr = dens*(relaxation*ln_r1+(1.-relaxation)*ln_r0)

  if (xnr > xnmin .and. xr > xmin) then
    lambda = cal_lambda_sb (2, xr, xnr)
    
    avd = 0.12407/lambda
    if (avd < 0.45e-3) then
      phibr = -1.
    else if (avd >= 0.45e-3) then
      phibr = 1000.*(avd - deq)
    endif
    
    crrn = -phibr*krr*xr*xnr * (1. + kar/lambda)**(-9.) * sqrt(rho0/dens)
    
    ln_r(5) = ln_r(5) - min(crrn/dens,ln_r1/dt)
  endif  
  
end subroutine

subroutine ccr_sb (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  ! Accretion of cloud droplets by rain drops (Seifert and Beheng 2006)
  implicit none
  real  :: dens, xc, xr, xnc, xnr
  real  :: tau, factau
  real  :: clr, clrn
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,ln_r,ln_c,lrho
  real, parameter :: kr=5.25
  
  dens = lrho0
  xc  = dens*(relaxation*lq_c1+(1.-relaxation)*lq_c0)
  xr  = dens*(relaxation*lq_r1+(1.-relaxation)*lq_r0)
  xnc = dens*(relaxation*ln_c1+(1.-relaxation)*ln_c0)
  xnr = dens*(relaxation*ln_r1+(1.-relaxation)*ln_r0)

  if ((xnc > xnmin .and. xc > xmin) .and. (xnr > xnmin .and. xr > xmin)) then
    tau = 1. - xc/(xc+xr)
    tau = min(max(tau,1.e-9),0.98)
    factau = (tau / (tau + 5.e-5))**4.
    
    ! accretion of cloud droplets by rain (rain formation)
    clr  = kr*xc*xr*factau * sqrt(rho0/dens)
    
    ! accretion (droplet removal)
    clrn = kr*xnc*xr*factau * sqrt(rho0/dens)
    
    lq_c(5) = lq_c(5) - min(clr/dens,lq_c1/dt)
    lq_r(5) = lq_r(5) + min(clr/dens,lq_c1/dt)
    ln_c(5) = ln_c(5) - min(clrn/dens,ln_c1/dt)
  endif
  
end subroutine

subroutine aucr_k (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  implicit none
  real :: xr, xnc, xnr
  real :: aunc, auq
  real :: xc, nuc, xprim, xref, k1, dens
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,ln_r,ln_c,lrho  
  real, parameter :: kau=7.98e10
  
  xprim = 1000.*4./3.*pi*(rauto)**3.
  xref  = 1000.*4./3.*pi*(40.e-6)**3.
  k1    = xref/xprim
  
  ! Convert to: kg/kg and cm-3
  dens = lrho0
  xc  = (relaxation*lq_c1+(1.-relaxation)*lq_c0)/dens
  xnc = (relaxation*ln_c1+(1.-relaxation)*ln_c0)/1000000.

  if (1000000.*xnc > xnmin .and. dens*xc > xmin) then
    auq  = log(k1*kau) + 4.22*log(xc) - 3.01*log(xnc)
    auq  = exp(auq)
    aunc = dens*auq/xprim
    
    lq_c(5) = lq_c(5) - min(dens*auq,lq_c1/dt)
    lq_r(5) = lq_r(5) + min(dens*auq,lq_c1/dt)
    if (have_ndrop) ln_c(5) = ln_c(5) - min(aunc,ln_c1/dt)
    if (have_nrain) ln_r(5) = ln_r(5) + min(aunc,ln_c1/dt)
  endif
  
end subroutine

subroutine scr_k (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  implicit none
  real :: nuc
  real :: dens, xc, xr, xnc, xnr
  real :: crrn, crcn
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,ln_r,ln_c,lrho
  real, parameter :: kscr=205.

  dens = lrho0
  xc  = (relaxation*lq_c1+(1.-relaxation)*lq_c0)/dens
  xr  = (relaxation*lq_r1+(1.-relaxation)*lq_r0)/dens
  xnc = (relaxation*ln_c1+(1.-relaxation)*ln_c0)/1000000.
  xnr = (relaxation*ln_r1+(1.-relaxation)*ln_r0)/1000000.

  if (1000000.*xnr > xnmin .and. dens*xr > xmin) then    
    crrn = log(kscr) + 0.6*log(xnr) + 1.55*log(xr)
    crrn = 1000000.*exp(crrn)
    
    ln_r(5) = ln_r(5) - min(crrn,ln_c1/dt)
  endif  
  
end subroutine

subroutine ccr_k (dt, lrho, lq_c, lq_r, ln_c, ln_r)
  implicit none
  real  :: dens, xc, xr, xnc, xnr
  real  :: tau, factau
  real  :: clr, clrn
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,ln_r,ln_c,lrho
  real, parameter :: kr=8.53
  
  dens = lrho0
  xc  = (relaxation*lq_c1+(1.-relaxation)*lq_c0)/dens
  xr  = (relaxation*lq_r1+(1.-relaxation)*lq_r0)/dens
  xnc = (relaxation*ln_c1+(1.-relaxation)*ln_c0)/1000000.
  xnr = (relaxation*ln_r1+(1.-relaxation)*ln_r0)/1000000.

  if ((1000000.*xnc > xnmin .and. dens*xc > xmin) .and. (1000000.*xnr > xnmin .and. dens*xr > xmin)) then
    clr  = log(kr) + 1.05*log(xc) + 0.98*log(xr)
    clr  = dens*exp(clr)
    clrn = clr * 1000000.*xnc/(dens*xc)
    
    lq_c(5) = lq_c(5) - min(clr,lq_c1/dt)
    lq_r(5) = lq_r(5) + min(clr,lq_c1/dt)
    if (have_ndrop) ln_c(5) = ln_c(5) - min(clrn,ln_c1/dt)
  endif
  
end subroutine

subroutine diff_cd (dt, qs, cp, qv, qc, nc, qr, nr, den, pp, tem, cd)

  implicit none
  integer :: i
  real  :: dt, fl, fv, vis, dif, lambda, gamm, fac
  real  :: cd0, dq, qc, nc, qr, nr, qv, qs, cp, tem, pp, den, esw, esi, qsi, sat, dav
  real, dimension(5) :: cd

  cd  = 0.0
  esw = cal_esat(tem)
  esi = cal_esati(tem)
  dif = cal_xfkd(tem,pp)
  vis = cal_xfmu(tem)
  sat = qv - qs

  cd0 = qs * ((c_p_v-c_v_v)*tem/(dif*esw) + (cal_flv(tem)/(kt*tem))*(cal_flv(tem)/((c_p_v-c_v_v)*tem) - 1.))

  ! Cloud drops
  if ((qc > xmin .or. nc > xnmin) .and. cd0 > 1.e-9) then
    fv=1.
    dav = (1./cm(1))**(1./am(1)) * cal_moment_sb(1, qc, nc, 1./am(1))
    call ventilation (1, vis, qc, nc, den, fv)
    
    cd(1) = 2.*pi*nc*dav*fv / cd0
  endif

  ! Rain drops
  if ((qr > xmin .or. nr > xnmin) .and. cd0 > 1.e-9) then
    fv=1.
    dav = (1./cm(2))**(1./am(2)) * cal_moment_sb(2, qr, nr, 1./am(2))
    if (sat < 0.) call ventilation (2, vis, qr, nr, den, fv)
    
    cd(2) = 2.*pi*nr*dav*fv / cd0
  endif

end subroutine

subroutine direct_solve_trapeze (dtt, cd, cp, pp, tt_1, qv_1, qc_1, qr_1, beta)

  implicit none
  real :: dtt, pp, cp, beta, cd(2)
  real :: tt_1, qv_1, qc_1, qr_1
  real :: gamm, sat0, dqsatt, dql, ss, dqc, dqr
!    
    dqc    = 0.0
    dqr    = 0.0
    dql    = cd(1) + cd(2)

    sat0   = qv_1 - cal_qsat(tt_1,pp)
    ss     = qv_1/cal_qsat(tt_1,pp)
    dqsatt = cal_dqsatdt(tt_1,pp)
    gamm   = dql*(1. + ss*dqsatt*cal_flv(tt_1)/cp)

    if (qc_1 > xmin) then 
      dqc  = max(min(dtt*cd(1)*sat0 / (1. + beta*dtt*gamm),qv_1),-qc_1)
      qc_1 = qc_1 + dqc
    endif
    
    if (qr_1 > xmin) then  
      dqr  = max(min(dtt*cd(2)*sat0 / (1. + beta*dtt*gamm),qv_1-dqc),-qr_1)
      qr_1 = qr_1 + dqr
    endif
    
    qv_1 = qv_1 - dqc - dqr
    tt_1 = tt_1 + (dqc + dqr)*cal_flv(tt_1)/cp
  
end subroutine

subroutine direct_solve_rkf (dtt, cd, cp, pp, tt_0, qv_0, qc_0, qr_0, tt_1, qv_1, qc_1, qr_1)

  implicit none
  real :: dt0, dt1, dtt, pp, cp, cd(2)
  real :: tt_0, qv_0, qc_0, qr_0, tt_1, qv_1, qc_1, qr_1
  real :: tt_y, qv_y, qc_y, qr_y, tt_z, qv_z, qc_z, qr_z
  real :: qc_rk(6), qr_rk(6), qv_rk(6), tt_rk(6)
  real :: tt_m, qv_m, sat_m
  real :: gamm, dqsatt, dql, ss
  real :: c(6), d(6), a(6,6)
  integer :: it0, it1, m, nit0, nit1

  c(1)=25./216.; c(2)=0.; c(3)=1408./2565.; c(4)=2197./4101.; c(5)=-1./5.; c(6)=0.
  d(1)=16./135.; d(2)=0.; d(3)=6656./12825.; d(4)=28561./56430.; d(5)=-9./50.; d(6)=2./55.
  a(1,1)=0.; a(1,2)=0.; a(1,3)=0.; a(1,4)=0.; a(1,5)=0.; a(1,6)=0.
  a(2,1)=0.25; a(2,2)=0.; a(2,3)=0.; a(2,4)=0.; a(2,5)=0.; a(2,6)=0.
  a(3,1)=3./32.; a(3,2)=9./32.; a(3,3)=0.; a(3,4)=0.; a(3,5)=0.; a(3,6)=0.
  a(4,1)=1932./2197.; a(4,2)=-7200./2197.; a(4,3)=7296./2197.; a(4,4)=0.; a(4,5)=0.; a(4,6)=0.
  a(5,1)=439./216.; a(5,2)=-8.; a(5,3)=3680./513.; a(5,4)=-845./4104.; a(5,5)=0.; a(5,6)=0.
  a(6,1)=-8./27.; a(6,2)=2.; a(6,3)=-3544./2565.; a(6,4)=1859./4104.; a(6,5)=-11./40.; a(6,6)=0.
    
  qv_1 = qv_0
  qc_1 = qc_0
  qr_1 = qr_0
  tt_1 = tt_0
  
  dql	 = cd(1) + cd(2)
  dqsatt = cal_dqsatdt(tt_0,pp)
  ss     = qv_0/cal_qsat(tt_0,pp)
  dt0    = 0.5/(dql*(1. + ss*dqsatt*cal_flv(tt_0)/cp))
  nit0   = floor(dtt/min(max(dt0,mindt),dtt))
  dt0    = dtt/real(nit0)
  
  do it0 = 1, nit0

    qv_rk=0.0
    qc_rk=0.0
    qr_rk=0.0
    tt_rk=0.0
    do m = 1, 6
      qv_m  = qv_1 + sum(a(m,:)*qv_rk)
      tt_m  = tt_1 + sum(a(m,:)*tt_rk)
      sat_m = qv_m - cal_qsat(tt_m,pp)
      
      if (qc_1 > xmin) then    
        qc_rk(m) = dt0*max(min(cd(1)*sat_m,qv_1/dt0),-qc_1/dt0)
      endif
      if (qr_1 > xmin) then    
        qr_rk(m) = dt0*max(min(cd(2)*sat_m,qv_1/dt0),-qr_1/dt0)
      endif
      
      qv_rk(m) = - qc_rk(m) - qr_rk(m)
      tt_rk(m) = (qc_rk(m) + qr_rk(m))*cal_flv(tt_0)/cp
    enddo
    qc_y = qv_1 + sum(c*qc_rk)
    qr_y = qc_1 + sum(c*qr_rk) 
    qv_y = qr_1 + sum(c*qv_rk) 
    tt_y = tt_1 + sum(c*tt_rk) 
    qc_z = qv_1 + sum(d*qc_rk)
    qr_z = qc_1 + sum(d*qr_rk) 
    qv_z = qr_1 + sum(d*qv_rk) 
    tt_z = tt_1 + sum(d*tt_rk)
    
    dt1  = 0.84*dt0*(rtol*dt0 / max(abs(qv_z-qv_y),abs(qc_z-qc_y),abs(qr_z-qr_y),abs(tt_z-tt_y)) )**(1./4.)
    nit1 = floor(dt0/min(max(dt1,mindt),dtt))
    dt1  = dt0/real(nit1)
    
    do it1 = 1, nit1
      qv_rk=0.0
      qc_rk=0.0
      qr_rk=0.0
      tt_rk=0.0
      do m = 1, 6
    	qv_m = qv_1 + sum(a(m,:)*qv_rk)
    	tt_m = tt_1 + sum(a(m,:)*tt_rk)
    	sat_m = qv_m-cal_qsat(tt_m,pp)
    	
    	if (qc_1 > xmin) then
    	  qc_rk(m) = dt1*max(min(cd(1)*sat_m,qv_1/dt1),-qc_1/dt1)
    	endif
    	if (qr_1 > xmin) then	 
    	  qr_rk(m) = dt1*max(min(cd(2)*sat_m,qv_1/dt1),-qr_1/dt1)
    	endif
	
	qv_rk(m) = - qc_rk(m) - qr_rk(m)
    	tt_rk(m) = (qc_rk(m) + qr_rk(m))*cal_flv(tt_0)/cp
      enddo
      qc_1 = qv_1 + sum(c*qc_rk)
      qr_1 = qc_1 + sum(c*qr_rk) 
      qv_1 = qr_1 + sum(c*qv_rk) 
      tt_1 = tt_1 + sum(c*tt_rk) 
    enddo
    
  enddo
!  
end subroutine

subroutine activ (dt, w, cp, cv, g, gsat, lnccn, lq_v, lq_c, ln_c, lp, ltem, mom_path)

  implicit none
  real  :: dt
  real  :: qs, ss, w, g, cp, cv, gsat, ssmax, ssmin, nc0, dnc, dqc
  real  :: qc, nc, qv, dens, ccn, tem, pp, gamm, taucd
  real  :: ak, rdcr, mu, erf, xxx, beta=0.5
  real, dimension(5) :: lq_v,lq_c,ln_c,lp,ltem,lnccn
  character(len=FIELD_NAME_LEN):: mom_path
  
  real :: eps, dtem, tem_new, tem_old, qv_new, qc_new, qs_new
  integer :: iter
  
  pp  = lp0
  tem = relaxation*ltem1+(1.-relaxation)*ltem0
  qv  = relaxation*lq_v1+(1.-relaxation)*lq_v0
  qc  = relaxation*lq_c1+(1.-relaxation)*lq_c0
  nc  = relaxation*ln_c1+(1.-relaxation)*ln_c0
  ccn = relaxation*ln_ccn1+(1.-relaxation)*ln_ccn0

  !Supersaturation arbitrarily limited to 5%
  qs = cal_qsat(tem,pp)
  ss = (qv - qs)/qs
  ss = min(ss,0.05)
    
  if (ss > 0.) then
  
    if (have_ndrop) then
    
      if (have_option(trim(mom_path)//'/detailed_activation')) then
        ak = 2.*cal_sigmawv(tem) / ((c_p_v-c_v_v)*tem*rhow)
        rdcr = sqrt(4*ak**3./(27.*k0))
        rdcr = (rdcr/ss)**(1./(beta + 1.))
      
        mu  = log(rd0) + log(sd0)*log(sd0)
        xxx = (log(rdcr) - mu) / (sqrt(2.)*log(sd0))
        erf = max(min(calerf(xxx,0),1.),-1.)
      
        nc0 = ccn/2. * (1. - erf)
      else if (have_option(trim(mom_path)//'/simple_activation')) then
        nc0 = ccn
      endif
      
      dnc = 1./dt*min(max((nc0 - nc),0.),ln_ccn1,ss*qs/xactiv)
     
      lq_c(5) = lq_c(5)  + xactiv*dnc
      lq_v(5) = lq_v(5)  - xactiv*dnc
      ln_c(5) = ln_c(5)  + dnc
      ltem(5) = ltem(5)  + xactiv*dnc*cal_flv(tem)/cp
!      lnccn(5)= lnccn(5) - dnc
    
    else
    
      qc_new = lq_c1
      qv_new = lq_v1
      tem_new = ltem1
      if (qc < xmin .and. qv > xmin) then
	iter=0
	eps=1.
        do while (eps > 1.e-6 .and. iter < 15)
          tem_old = tem_new
	  qs_new = cal_qsat(tem_new,pp)
          dqc = (qv_new-qs_new) / (1. + cal_flv(tem)**2.*qs_new/(cp * (cp-cv)*tem_new**2.))
	  tem_new = tem_new+max(dqc,-qc_new)*cal_flv(tem)/cp
          qc_new = max(qc_new+dqc,0.)
          qv_new = max(qv_new-dqc,0.)

          eps=abs(tem_old-tem_new)
	  iter=iter+1
        enddo
	dqc = max(min((qc_new - lq_c1)/dt,lq_v1/dt),0.0)
      
        lq_c(5) = lq_c(5) + dqc
        lq_v(5) = lq_v(5) - dqc
        ltem(5) = ltem(5) + dqc*cal_flv(tem)/cp
      endif
      
    endif
    
  endif
  
end subroutine

subroutine ventilation (i, vis, q, n, dens, fv)
  implicit none
  integer :: i, k
  real  :: dens, q, n, vis, fv
  real  :: beta, re, const, dav, velav
  real, parameter :: avent=0.78, bvent=0.308
  
  fv = 1.
  beta=0.5*(be_sb(i)+3./am(i))
  if (i /= 3 .and. i /= 1) then
    dav = (1./cm(i))**(1./am(i))*cal_moment_sb(i, q, n, 1./am(i))
    velav = al_sb(i)*cal_moment_sb(i, q, n, be_sb(i))*sqrt(rho0/dens)
    re=dav*velav/vis

    const=(cal_gamma((beta+nu_sb(i)+1.)/mu_sb(i))*cal_gamma((nu_sb(i)+1.)/mu_sb(i))) / &
          (cal_gamma((1./am(i)+nu_sb(i)+1.)/mu_sb(i))**(3./2.)*cal_gamma((be_sb(i)+nu_sb(i)+1.)/mu_sb(i))**(1./2.))
    
    fv = avent + const*bvent*0.892*sqrt(re)
  endif
  
end subroutine

subroutine sedim (dt, i, lq0, ln0, lq1, ln1, lqv, lnv)
  implicit none
  integer :: i
  real :: q, n, den, vel1, vel2, lambda
  
  real :: dt, lq1, ln1, lq0, ln0, lqv, lnv
  
  den = lrho0
  q = relaxation*lq1+(1.-relaxation)*lq0
  n = relaxation*ln1+(1.-relaxation)*ln0

  if (q > xmin .or. n > xnmin) then
    lambda = cal_lambda_sb (i, q, n)
    
    vel1 = 9.65 - 10.3*(1.+600./lambda)**(-4.)
    vel2 = 9.65 - 10.3*(1.+600./lambda)**(-1.)
 
    lqv = max(vel1,0.) * sqrt(rho0/den)
    lnv = max(vel2,0.) * sqrt(rho0/den)
  endif  
  
end subroutine

function cal_moment_sb (i, q, n, k)
  integer :: i
  real :: q, n, k
  real :: lambda, cal_moment_sb
  
  if (q > xmin .or. n > xnmin) then
    lambda=cal_lambda_sb(i, q, n)
    cal_moment_sb=cal_gamma((k+nu_sb(i)+1.)/mu_sb(i))/cal_gamma((nu_sb(i)+1.)/mu_sb(i)) * lambda**(-k/mu_sb(i))
  endif

  return

end function

function cal_lambda_sb (i, q, n)
  integer :: i
  real :: q, n
  real :: lambda, cal_lambda_sb
  
  cal_lambda_sb=0.
  if (q > xmin .or. n > xnmin) then
    lambda=(cal_gamma((nu_sb(i)+1.)/mu_sb(i))/cal_gamma((nu_sb(i)+2.)/mu_sb(i))*(q+xmin)/(n+xnmin))**(-mu_sb(i))
    
    if (i == 2) cal_lambda_sb=min(max(lambda,100.),1.e+6)
    if (i == 1) cal_lambda_sb=min(max(lambda,1.e+10),1.e+16)
  endif
  
  return 

end function

subroutine assemble_micro (dt, cp, lqv, lqc, lnc, lqr, lnr, lnccn, ltem)

  integer :: i
  real  :: dt, cp
  real  :: fqc, fqr, tem
  real, dimension(5) :: lqv,lqr,lnr,lqc,lnc,lnccn,ltem

  if (lq_v1 < xmin) then
    fqc = lq_c1/(lq_c1+lq_r1)
    fqr = lq_r1/(lq_c1+lq_r1)
    lq_c1 = lq_c1 + max(fqc*lq_v1,0.0)
    lq_r1 = lq_r1 + max(fqr*lq_v1,0.0)
    lq_v1 = 0.0
  endif
  
  lq_v(5) = (max(lq_v1,0.) - lq_v0)/dt
  
  lq_c(5) = (max(lq_c1,0.) - lq_c0)/dt  

  lq_r(5) = (max(lq_r1,0.) - lq_r0)/dt

  ln_c(5) = (max(ln_c1,0.) - ln_c0)/dt

  ln_r(5) = (max(ln_r1,0.) - ln_r0)/dt

!  ln_ccn(5) = (max(ln_ccn1,0.) - ln_ccn0)/dt

  tem = relaxation*ltem1 +(1.-relaxation)*ltem0
  ltem(5) = -lq_v(5)*cal_flv(tem)/cp
  
end subroutine

end subroutine microphysics_2mom
  
subroutine microphysics_1mom(state,sat_adj,time,dt,qs,cp,lq_t,lq_v,lq_r,lq_c,lq_i,lq_g,lq_s,ln_r,ln_c,ln_i,ln_g,ln_s,ltem,lp,lrho)
  implicit none
  type(state_type),intent(inout) :: state
  logical :: sat_adj
  real :: time, dt, qs, cp
  real :: lq_v0, lq_c0, lq_r0, lq_i0, lq_g0, lq_s0, ltem0, ln_c0, ln_i0
  real :: lq_v1, lq_c1, lq_r1, lq_i1, lq_g1, lq_s1, ltem1, ln_c1, ln_i1
  real, dimension(5) :: lq_t,lq_v,lq_r,lq_c,lq_i,lq_g,lq_s,ln_r,ln_c,ln_i,ln_g,ln_s,ltem,lp,lrho
      
  ! data arrays for the pass into the fortran routine
  ! These have the following order
  
  ! (1) Old Timelevel Value (input, projected)
  ! (2) Previous nonlinear iteration (input, projected)
  ! (3) Current nonlinear iteration (input, projected)
  ! (4) Sinking Velocity (output, on Microphysics mesh)
  ! (5) Microphysics forcing (output, on Microphysics mesh)
    
  lq_v(5) = 0.
  lq_c(5) = 0.
  lq_r(5) = 0.
  lq_i(5) = 0.
  lq_g(5) = 0.
  lq_s(5) = 0.
  ln_c(5) = 0.
  ln_i(5) = 0.
  ltem(5) = 0.
  
  lq_v0 = lq_v(1); lq_v1 = lq_v(3)
  lq_c0 = lq_c(1); lq_c1 = lq_c(3)
  lq_r0 = lq_r(1); lq_r1 = lq_r(3)
  lq_i0 = lq_i(1); lq_i1 = lq_i(3)
  lq_g0 = lq_g(1); lq_g1 = lq_g(3)
  lq_s0 = lq_s(1); lq_s1 = lq_s(3)
  ln_c0 = ln_c(1); ln_c1 = ln_c(3)
  ln_i0 = ln_i(1); ln_i1 = ln_i(3)
  ltem0 = ltem(1); ltem1 = ltem(3)
  
  !
  !  1.: Sedimentation velocities
  !
    
  call sedim (dt, 2, lq_r1, lq_r, ln_r, lrho)
 
  if (has_scalar_field(state,'Qice')) then
    call sedim (dt, 3, lq_i1, lq_i, ln_i, lrho)
    call sedim (dt, 4, lq_g1, lq_g, ln_g, lrho)
    call sedim (dt, 5, lq_s1, lq_s, ln_s, lrho)
  endif
  
  !
  !  2.: Evaluate microphysical processes: diffusion/evaporation/melting...
  !

  ! condensation/evaporation rain: v <-> c
  if (.not.sat_adj) call diff (dt, 1, cp, lq_v, lq_c, lq_c1, ln_c, lrho, lp, ltem)
  
  ! condensation/evaporation rain: v <-> r
  call diff (dt, 2, cp, lq_v, lq_r, lq_r1, ln_r, lrho, lp, ltem)
 
  if (has_scalar_field(state,'Qice')) then
  
    ! condensation/evaporation ice: v <-> i (including melted i)
    call diff (dt, 3, cp, lq_v, lq_i, lq_i1, ln_i, lrho, lp, ltem)
  
    ! condensation/evaporation grau: v <-> g (including melted g)
    call diff (dt, 4, cp, lq_v, lq_g, lq_g1, ln_g, lrho, lp, ltem)
  
    ! condensation/evaporation snow: v <-> s (including melted s)
    call diff (dt, 5, cp, lq_v, lq_s, lq_s1, ln_s, lrho, lp, ltem)
    
    ! melting ice: i -> c
    call melt_ice (dt, cp, lq_i, ln_i, lq_c, ltem)
  
    ! melting graupel: g -> r
    call melt (dt, 4, cp, lq_r, lq_g, lq_g1, ln_g, lrho, ltem)
  
    ! melting snow: s -> r
    call melt (dt, 5, cp, lq_r, lq_s, lq_s1, ln_s, lrho, ltem)
  
  endif

  lq_v1 = lq_v1 + dt*lq_v(5)
  lq_c1 = lq_c1 + dt*lq_c(5)
  lq_r1 = lq_r1 + dt*lq_r(5)
  lq_i1 = lq_i1 + dt*lq_i(5)
  lq_g1 = lq_g1 + dt*lq_g(5)
  lq_s1 = lq_s1 + dt*lq_s(5)
  ltem1 = ltem0 + dt*ltem(5)
  
  lq_v(5) = 0.
  lq_c(5) = 0.
  lq_r(5) = 0.
  lq_i(5) = 0.
  lq_g(5) = 0.
  lq_s(5) = 0.
  ltem(5) = 0.
  
  !
  !  3.: Evaluate microphysical processes: collision-coalescence
  !
   
  ! Calculate autoconversion rates: c -> r
  call auto_drop (dt, lrho, lq_c, lq_r)
  
  ! Calculate accretion rates: c+r -> r
  call acc_cr (dt, lrho, lq_c, lq_r)
  
  if (has_scalar_field(state,'Qice')) then
   
    ! Calculate autoconversion rates: i -> s
    call auto_ice (dt, lrho, lq_i, lq_s, ltem)
  
    ! Calculate accretion rates: c+i -> g
!    call acc_ci (dt, lrho, lq_c, lq_i, lq_g, lq_r, ltem)
  
    ! Calculate accretion rates: c+g -> g (includes riming enhanced melting) 
    call acc_cg (dt, lrho, lq_c, ln_c, lq_g, lq_r, ltem)
  
    ! Calculate accretion rates: c+s -> s (includes riming enhanced melting) 
    call acc_cs (dt, lrho, lq_c, ln_c, lq_s, lq_r, ltem)
  
    ! Calculate accretion rates: r+i -> g
    call acc_ri (dt, lrho, lq_r, lq_i, ln_i, lq_g, ltem)
  
    ! Calculate accretion rates: r+s -> g (includes riming enhanced melting) 
    call acc_rs (dt, lrho, lq_r, lq_s, lq_g, ltem)
  
    ! Calculate accretion rates: r+g -> g (includes riming enhanced melting) 
    call acc_rg (dt, lrho, lq_r, lq_g, ltem)
  
    ! Calculate accretion rates: i+g -> g
    call acc_ig (dt, lrho, lq_g, lq_i, ltem)
  
    ! Calculate accretion rates: s+g -> g
    call acc_sg (dt, lrho, lq_s, lq_g, ltem)
  
    ! Calculate accretion rates: i+s -> s
    call acc_is (dt, lrho, lq_i, ln_i, lq_s, ltem)
  
  endif

  lq_c1 = lq_c1 + dt*lq_c(5)
  lq_r1 = lq_r1 + dt*lq_r(5)
  lq_i1 = lq_i1 + dt*lq_i(5)
  lq_g1 = lq_g1 + dt*lq_g(5)
  lq_s1 = lq_s1 + dt*lq_s(5)
  ln_c1 = ln_c1 + dt*ln_c(5)
  ln_i1 = ln_i1 + dt*ln_i(5)
  
  lq_c(5) = 0.
  lq_r(5) = 0.
  lq_i(5) = 0.
  lq_g(5) = 0.
  lq_s(5) = 0.
  ln_c(5) = 0.
  ln_i(5) = 0.
  
  !
  !  4.: Production: freezing/activation
  !
 
  if (has_scalar_field(state,'Qice')) then

    ! Freezing heterogeneousd
    call freezing_D10 (dt, ltem, lrho, lq_i, ln_i)
  
    ! Freezing homogeneous: c -> i and r -> g
    call freezing_hom (dt, ltem, lrho, lq_i, ln_i, lq_g, lq_c, ln_c, lq_r)
  
  endif
  
  lq_c1 = lq_c1 + dt*lq_c(5)
  lq_r1 = lq_r1 + dt*lq_r(5)
  lq_i1 = lq_i1 + dt*lq_i(5)
  lq_g1 = lq_g1 + dt*lq_g(5)
  ln_c1 = ln_c1 + dt*ln_c(5)
  ln_i1 = ln_i1 + dt*ln_i(5)
  
  lq_c(5) = 0.
  lq_r(5) = 0.
  lq_i(5) = 0.
  lq_g(5) = 0.
  ln_c(5) = 0.
  ln_i(5) = 0.
    
  !
  !  5.: Assemble microphysics sources 
  !  
  
  call assemble_micro (dt, cp, lq_v, lq_c, lq_r, lq_i, lq_g, lq_s, ln_c, ln_i, ltem)
      
  
  ! Treated in calculate_diagnostic_microphysics
  
contains

subroutine freezing_D10(dt, ltem, lrho, lq_i, ln_i)
  implicit none
  real :: ni, qi, tem, den
  
  real :: dt
  real :: dq, dn, ni0, qi0
  real, dimension(5) :: lq_i,ln_i,lrho,ltem
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qi  = relaxation*lq_i1+(1.-relaxation)*lq_i(3)
  ni  = relaxation*ln_i1+(1.-relaxation)*ln_i(3)
   
  if (ltem(1) < 268.15) then
    ni0 = exp(32. - 0.125*tem)
    qi0 = 1.e-12*ni0/den
    dq	= qi0 - qi
    dn  = ni0 - ni

    ln_i(5) = ln_i(5) + max(dn,0.)
    lq_i(5) = lq_i(5) + max(dq/dt,0.)
  endif
  
end subroutine freezing_D10

subroutine freezing_hom(dt, ltem, lrho, lq_i, ln_i, lq_g, lq_c, ln_c, lq_r)
  implicit none
  real :: dn, qc, nc, qr, tem, den
  
  real :: dt
  real, dimension(5) :: lq_i,ln_i,lq_g,lq_c,ln_c,lq_r,lrho,ltem
 
  if (lq_c1 > xmin .and. ltem1 < 235.) then
    lq_i(5) = lq_i(5) + max(lq_c1/dt,0.)
    lq_c(5) = lq_c(5) - max(lq_c1/dt,0.)
    ln_i(5) = ln_i(5) + max(ln_c1,0.)
  endif
 
  if (lq_r1 > xmin .and. ltem1 < 235.) then
    lq_g(5) = lq_g(5) + max(lq_r1/dt,0.)
    lq_r(5) = lq_r(5) - max(lq_r1/dt,0.)
  endif
  
end subroutine freezing_hom

subroutine auto_drop (dt, lrho, lq_c, lq_r)
  implicit none
  real :: qc, qr, den, Pauto
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,lrho 
  real, parameter :: k1=5.5e-4	!1.e-3
  
  den = lrho(3)
  qc  = relaxation*lq_c1+(1.-relaxation)*lq_c(3)
  qr  = relaxation*lq_r1+(1.-relaxation)*lq_r(3)

  if (lq_c1 > mauto/1000.) then
    Pauto=k1*(qc-mauto/1000.)
  else
    Pauto=0.
  endif
    
  lq_c(5) = lq_c(5) - min(max(Pauto,0.0),lq_c1/dt)
  lq_r(5) = lq_r(5) + min(max(Pauto,0.0),lq_c1/dt)
  
end subroutine

subroutine auto_ice (dt, lrho, lq_i, lq_s, ltem)
  implicit none
  real :: qi, qs, den, tem, Pauto
  
  real :: dt
  real, dimension(5) :: lq_i,lq_s,lrho,ltem
  real, parameter :: k1=1.e-3, miauto=8.e-5
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qi  = relaxation*lq_i1+(1.-relaxation)*lq_i(3)
  qs  = relaxation*lq_s1+(1.-relaxation)*lq_s(3)

  if (lq_i1 > miauto .and. ltem1 < 273.15) then
    Pauto=k1*(qi-miauto)
  else
    Pauto=0.
  endif
    
  lq_i(5) = lq_i(5) - min(max(Pauto,0.0),lq_i1/dt)
  lq_s(5) = lq_s(5) + min(max(Pauto,0.0),lq_i1/dt)
  
end subroutine

subroutine acc_cr (dt, lrho, lq_c, lq_r)
  implicit none
  real :: qc, qr, den, lambda, Paccc
  
  real :: dt
  real, dimension(5) :: lq_r,lq_c,lrho
  real, parameter :: Eacc=1.
  
  !
  !  Accretion: c + r --> r
  !
  
  den = lrho(3)
  qr  = relaxation*lq_r1+(1.-relaxation)*lq_r(3)
  qc  = relaxation*lq_c1+(1.-relaxation)*lq_c(3)

  if (lq_c1 > xmin .and. lq_r1 > xmin) then
    lambda = (cm(2)*N0r*cal_gamma(nu(2)+am(2)+1.)/(den*qr))**(1./(nu(2)+am(2)+1.))
    Paccc = 0.25*pi*Eacc*N0r*ctv(2)*cal_gamma(3.+btv(2)+nu(2)) * sqrt(rho0/den)
    Paccc = Paccc*qc*lambda**(-3.-btv(2)-nu(2))
    
    lq_c(5) = lq_c(5) - min(max(Paccc,0.0),lq_c1/dt)
    lq_r(5) = lq_r(5) + min(max(Paccc,0.0),lq_c1/dt)
  endif  
  
end subroutine

subroutine acc_cg (dt, lrho, lq_c, ln_c, lq_g, lq_r, ltem)
  implicit none
  real :: qc, nc, qg, den, tem, lambda, Paccc, Paccg
  
  real :: dt
  real, dimension(5) :: lq_g,lq_c,ln_c,lq_r,lrho,ltem
  real, parameter :: Eacc=1.
  
  !
  !  Accretion: c + g --> g
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  nc  = relaxation*ln_c1+(1.-relaxation)*ln_c(3)
  qc  = relaxation*lq_c1+(1.-relaxation)*lq_c(3)
  qg  = relaxation*lq_g1+(1.-relaxation)*lq_g(3)

  if (lq_c1 > xmin .and. lq_g1 > xmin) then
    lambda = (cm(4)*N0g*cal_gamma(nu(4)+am(4)+1.)/(den*qg))**(1./(nu(4)+am(4)+1.))
    
    Paccc = 0.25*pi*Eacc*N0g*ctv(4)*cal_gamma(3.+btv(4)+nu(4)) * sqrt(rho0/den)
    Paccc = Paccc*qc*lambda**(-3.-btv(4)-nu(4))
    
    Paccg = 0.25*pi*Eacc*N0g*ctv(4)*cm(4)*cal_gamma(3.+am(4)+btv(4)+nu(4)) * sqrt(rho0/den)
    Paccg = Paccg*nc*lambda**(-3.-am(4)-btv(4)-nu(4))
    
    lq_c(5) = lq_c(5) - min(max(Paccc,0.0),lq_c1/dt)
    if (ltem1 >= 273.15) then
      Paccg   = -c_p_l/cal_flm(tem)*(tem - 273.15) * Paccg
      lq_g(5) = lq_g(5) - min(max(Paccg,0.0),lq_g1/dt)
      lq_r(5) = lq_r(5) + min(max(Paccg,0.0),lq_g1/dt) + min(max(Paccc,0.0),lq_c1/dt)
    else
      lq_g(5) = lq_g(5) + min(max(Paccc,0.0),lq_c1/dt)
    endif
  endif  
  
end subroutine

subroutine acc_cs (dt, lrho, lq_c, ln_c, lq_s, lq_r, ltem)
  implicit none
  real :: qc, nc, qs, den, tem, lambda, Paccc, Paccs
  
  real :: dt
  real, dimension(5) :: lq_s,lq_c,ln_c,lq_r,lrho,ltem
  real, parameter :: Eacc=1.
  
  !
  !  Accretion: c + s --> s
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  nc  = relaxation*ln_c1+(1.-relaxation)*ln_c(3)
  qc  = relaxation*lq_c1+(1.-relaxation)*lq_c(3)
  qs  = relaxation*lq_s1+(1.-relaxation)*lq_g(3)

  if (lq_c1 > xmin .and. lq_s1 > xmin) then
    lambda = (cm(5)*N0s*cal_gamma(nu(5)*am(5)+1.)/(den*qs))**(1./(nu(5)*am(5)+1.))
    
    Paccc = 0.25*pi*Eacc*N0s*ctv(5)*cal_gamma(3.+btv(5)+nu(5)) * sqrt(rho0/den)
    Paccc = Paccc*qc*lambda**(-3.-btv(5)-nu(5))

    Paccs = 0.25*pi*Eacc*N0s*ctv(5)*cm(5)*cal_gamma(3.+am(5)+btv(5)+nu(5)) * sqrt(rho0/den)
    Paccs = Paccs*nc*lambda**(-3.-am(5)-btv(5)-nu(5))

    lq_c(5) = lq_c(5) - min(max(Paccc,0.0),lq_c1/dt)
    if (ltem1 >= 273.15) then
      Paccs   = -c_p_l/cal_flm(tem)*(tem - 273.15) * Paccs     
      lq_s(5) = lq_s(5) - min(max(Paccs,0.0),lq_s1/dt)
      lq_r(5) = lq_r(5) + min(max(Paccs,0.0),lq_s1/dt) + min(max(Paccc,0.0),lq_c1/dt)
    else
      lq_s(5) = lq_s(5) + min(max(Paccc,0.0),lq_c1/dt)
    endif
  endif  
  
end subroutine

subroutine acc_ri (dt, lrho, lq_r, lq_i, ln_i, lq_g, ltem)
  implicit none
  real :: qr, qi, ni, den, tem, lambda, Pacci, Paccr
  
  real :: dt
  real, dimension(5) :: lq_r,lq_i,ln_i,lq_g,lrho,ltem
  real, parameter :: Eacc=1.
  
  !
  !  Accretion: r + i --> g
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  ni  = relaxation*ln_i1+(1.-relaxation)*ln_i(3)
  qi  = relaxation*lq_i1+(1.-relaxation)*lq_i(3)
  qr  = relaxation*lq_r1+(1.-relaxation)*lq_r(3)

  if (lq_r1 > xmin .and. lq_i1 > xmin .and. ltem1 < 273.15) then
    lambda = (cm(2)*N0r*cal_gamma(nu(2)+am(2)+1.)/(den*qr))**(1./(nu(2)+am(2)+1.))
    
    Pacci = 0.25*pi*Eacc*N0r*ctv(2)*cal_gamma(3.+btv(2)+nu(2)) * sqrt(rho0/den)
    Pacci = Pacci*qi*lambda**(-3.-btv(2)-nu(2))
    
    Paccr = 0.25*pi*Eacc*N0r*cm(2)*ctv(2)*cal_gamma(3.+am(2)+btv(2)+nu(2)) * sqrt(rho0/den)
    Paccr = Paccr*ni*lambda**(-3.-am(2)-btv(2)-nu(2))
    
    lq_r(5) = lq_r(5) - min(max(Paccr,0.0),lq_r1/dt)
    lq_i(5) = lq_i(5) - min(max(Pacci,0.0),lq_i1/dt)
    lq_g(5) = lq_g(5) + min(max(Paccr,0.0),lq_r1/dt) + min(max(Pacci,0.0),lq_i1/dt)
  endif  
  
end subroutine

subroutine acc_rg (dt, lrho, lq_r, lq_g, ltem)
  implicit none
  real :: qr, qg, vr, vg, den, tem, lambdar, lambdag, Paccr, Paccg
  
  real :: dt
  real, dimension(5) :: lq_r,lq_g,lrho,ltem
  real, parameter :: Eacc=1.
  
  !
  !  Accretion: r + g --> g
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qr  = relaxation*lq_r1+(1.-relaxation)*lq_r(3)
  qg  = relaxation*lq_g1+(1.-relaxation)*lq_g(3)
  vr = lq_r(4)
  vg = lq_g(4)

  if (qr > xmin .and. qg > xmin) then
    lambdar = (cm(2)*N0r*cal_gamma(nu(2)+am(2)+1.)/(den*qr))**(1./(nu(2)+am(2)+1.))
    lambdag = (cm(4)*N0g*cal_gamma(nu(4)+am(4)+1.)/(den*qg))**(1./(nu(4)+am(4)+1.))
    
    Paccr = 0.25*pi*Eacc*N0g*abs(vr - vg)
    Paccr = Paccr*qr*(cal_gamma(nu(4)+3.)*lambdag**(-nu(4)-3.) +		&
    	  	      (am(2)+nu(2)+1.)*cal_gamma(nu(4)+2.)*lambdag**(-nu(4)-2.)*lambdar**(-1.) +	&
		      (am(2)+nu(2)+2.)*(am(2)+nu(2)+1.)*cal_gamma(nu(4)+1.)*lambdag**(-nu(4)-1.)*lambdar**(-2.))
    
    Paccg = 0.25*pi*Eacc*N0r*abs(vr - vg)
    Paccg = Paccg*qg*(cal_gamma(nu(2)+3.)*lambdar**(-nu(2)-3.) +		&
    	  	      (am(4)+nu(4)+1.)*cal_gamma(nu(2)+2.)*lambdar**(-nu(2)-2.)*lambdag**(-1.) +	&
		      (am(4)+nu(4)+2.)*(am(4)+nu(4)+1.)*cal_gamma(nu(2)+1.)*lambdar**(-nu(2)-1.)*lambdag**(-2.))
    
    if (ltem1 >= 273.15) then
      Paccg   = -c_p_l/cal_flm(tem)*(tem - 273.15) * Paccg   
      lq_r(5) = lq_r(5) + min(max(Paccg,0.0),lq_g1/dt)
      lq_g(5) = lq_g(5) - min(max(Paccg,0.0),lq_g1/dt)    
    else
      lq_r(5) = lq_r(5) - min(max(Paccr,0.0),lq_r1/dt)
      lq_g(5) = lq_g(5) + min(max(Paccr,0.0),lq_r1/dt)
    endif
  endif  
  
end subroutine

subroutine acc_rs (dt, lrho, lq_r, lq_s, lq_g, ltem)
  implicit none
  real :: qr, qs, vr, vs, den, tem, lambdar, lambdas, Paccs, Paccr
  
  real :: dt
  real, dimension(5) :: lq_r,lq_s,lq_g,lrho,ltem
  real, parameter :: Eacc=1.
  
  !
  !  Accretion: r + s --> g
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qr  = relaxation*lq_r1+(1.-relaxation)*lq_r(3)
  qs  = relaxation*lq_s1+(1.-relaxation)*lq_s(3)
  vr = lq_r(4)
  vs = lq_s(4)

  if (lq_r1 > xmin .and. lq_s1 > xmin) then
    lambdar = (cm(2)*N0r*cal_gamma(nu(2)+am(2)+1.)/(den*qr))**(1./(nu(2)+am(2)+1.))
    lambdas = (cm(5)*N0s*cal_gamma(nu(5)+am(5)+1.)/(den*qs))**(1./(nu(5)+am(5)+1.))
    
    Paccr = 0.25*pi*Eacc*N0s*abs(vr - vs)
    Paccr = Paccr*qr*(cal_gamma(nu(5)+3.)*lambdas**(-nu(5)-3.) +		&
    	  	      (am(2)+nu(2)+1.)*cal_gamma(nu(5)+2.)*lambdas**(-nu(5)-2.)*lambdar**(-1.) +	&
		      (am(2)+nu(2)+2.)*(am(2)+nu(2)+1.)*cal_gamma(nu(5)+1.)*lambdas**(-nu(5)-1.)*lambdar**(-2.))

    Paccs = 0.25*pi*Eacc*N0r*abs(vr - vs) 
    Paccs = Paccs*qs*(cal_gamma(nu(2)+3.)*lambdar**(-nu(2)-3.) +		&
    		      (am(5)+nu(5)+1.)*cal_gamma(nu(2)+2.)*lambdar**(-nu(2)-2.)*lambdas**(-1.) +	&
		      (am(5)+nu(5)+2.)*(am(5)+nu(5)+1.)*cal_gamma(nu(2)+1.)*lambdar**(-nu(2)-1.)*lambdas**(-2.))

    if (ltem1 >= 273.15) then    
      Paccs   = -c_p_l/cal_flm(tem)*(tem - 273.15) * Paccs          
      lq_s(5) = lq_s(5) - min(max(Paccs,0.0),lq_s1/dt)
      lq_r(5) = lq_r(5) + min(max(Paccs,0.0),lq_s1/dt)
    else
      lq_r(5) = lq_r(5) - min(max(Paccr,0.0),lq_r1/dt)
      lq_s(5) = lq_s(5) - min(max(Paccs,0.0),lq_s1/dt)
      lq_g(5) = lq_g(5) + min(max(Paccr,0.0),lq_r1/dt) + min(max(Paccs,0.0),lq_s1/dt)
    endif
  endif  
  
end subroutine

subroutine acc_ig (dt, lrho, lq_g, lq_i, ltem)
  implicit none
  real :: qg, qi, den, tem, lambda, Pacci
  
  real :: dt
  real, dimension(5) :: lq_i,lq_g,lrho,ltem
  real :: Eacc=0.8
  
  !
  !  Accretion: i + g --> g
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qg  = relaxation*lq_g1+(1.-relaxation)*lq_g(3)
  qi  = relaxation*lq_i1+(1.-relaxation)*lq_i(3)
  
  if (lq_g1 > xmin .and. lq_i1 > xmin .and. ltem1 < 273.15) then
    Eacc=min(exp(0.07*(tem-273.15)),1.)
    lambda = (cm(4)*N0g*cal_gamma(nu(4)+am(4)+1.)/(den*qg))**(1./(nu(4)+am(4)+1.))
    
    Pacci = 0.25*pi*Eacc*N0g*ctv(4)*cal_gamma(3.+btv(4)+nu(4)) * sqrt(rho0/den)
    Pacci = Pacci*qi*lambda**(-3.-btv(4)-nu(4))
    
    lq_i(5) = lq_i(5) - min(max(Pacci,0.0),lq_i1/dt)
    lq_g(5) = lq_g(5) + min(max(Pacci,0.0),lq_i1/dt)
  endif  
  
end subroutine

subroutine acc_sg (dt, lrho, lq_s, lq_g, ltem)
  implicit none
  real :: qs, qg, vs, vg, den, tem, lambdas, lambdag, Paccs
  
  real :: dt
  real, dimension(5) :: lq_s,lq_g,lrho,ltem
  real :: Eacc=0.8
  
  !
  !  Accretion: s + g --> g
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qg  = relaxation*lq_g1+(1.-relaxation)*lq_g(3)
  qs  = relaxation*lq_s1+(1.-relaxation)*lq_s(3)
  
  if (lq_s1 > xmin .and. lq_g1 > xmin) then
    Eacc=min(exp(0.09*(tem-273.15)),1.)
    lambdas = (cm(5)*N0s*cal_gamma(nu(5)+am(5)+1.)/(den*qs))**(1./(nu(5)+am(5)+1.))
    lambdag = (cm(4)*N0g*cal_gamma(nu(4)+am(4)+1.)/(den*qg))**(1./(nu(4)+am(4)+1.))
    
    Paccs = 0.25*pi*Eacc*N0g*abs(vs - vg)
    Paccs = Paccs*qs*(cal_gamma(nu(4)+3.)*lambdag**(-nu(4)-3.) +		&
    		      (am(5)+nu(5)+1.)*cal_gamma(nu(4)+2.)*lambdag**(-nu(4)-2.)*lambdas**(-1.) +	&
		      (am(5)+nu(5)+2.)*(am(5)+nu(5)+1.)*cal_gamma(nu(4)+1.)*lambdag**(-nu(4)-1.)*lambdas**(-2.))
    
    lq_s(5) = lq_s(5) - min(max(Paccs,0.0),lq_s1/dt)
    lq_g(5) = lq_g(5) + min(max(Paccs,0.0),lq_s1/dt)
  endif  
  
end subroutine

subroutine acc_is (dt, lrho, lq_i, ln_i, lq_s, ltem)
  implicit none
  real :: vi, vs, qi, ni, qs, den, tem, lambdai, lambdas, Pacci
  
  real :: dt
  real, dimension(5) :: lq_i,ln_i,lq_s,lrho,ltem
  real :: Eacc=0.8
  
  !
  !  Accretion: i + s --> s
  !
  
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qi  = relaxation*lq_i1+(1.-relaxation)*lq_i(3)
  qs  = relaxation*lq_s1+(1.-relaxation)*lq_s(3)
  ni  = relaxation*ln_i1+(1.-relaxation)*ln_i(3)
  vi = lq_i(4)
  vs = lq_s(4)

  if (lq_i1 > xmin .and. lq_s1 > xmin .and. ltem1 < 273.15) then
    Eacc = min(exp(0.07*(tem-273.15)),1.)
    lambdai = (cm(3)*ni*cal_gamma(nu(3)+am(3)+1.)/(den*cal_gamma(nu(3)+1.)*qi))**(1./am(3))
    lambdas = (cm(5)*N0s*cal_gamma(nu(5)+am(5)+1.)/(den*qs))**(1./(nu(5)+am(5)+1.))
    
    Pacci = 0.25*pi*Eacc*N0s*abs(vi - vs)
    Pacci = Pacci*qi*(cal_gamma(nu(5)+3.)*lambdas**(-nu(5)-3.) +		&
    	              (am(3)+nu(3)+1.)*cal_gamma(nu(3)+2.)*lambdas**(-nu(5)-2.)*lambdai**(-1.) +	&
		      (am(3)+nu(3)+2.)*(am(3)+nu(3)+1.)*cal_gamma(nu(5)+1.)*lambdas**(-nu(5)-1.)*lambdai**(-2.))
    
    lq_i(5) = lq_i(5) - min(max(Pacci,0.0),lq_i1/dt)
    lq_s(5) = lq_s(5) + min(max(Pacci,0.0),lq_i1/dt)
  endif  
  
end subroutine

subroutine diff (dt, i, cp, lq_v, lq, lq_1, ln, lrho, lp, ltem)
  implicit none
  integer :: i
  real  :: dt, lq_1, fl, fv, cd, vis, dif, lambda, gamm, beta
  real  :: dav, dq, q, n, qv, qs, cp, tem, pp, den, es, dqs, sat
  real, dimension(5) :: lq_v,lq,ln,lp,ltem,lrho

  pp  = lp(3)
  den = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  qv  = relaxation*lq_v1+(1.-relaxation)*lq_v(3)
  q   = relaxation*lq_1+(1.-relaxation)*lq(3)
  n   = ln(3)

  if (lq_1 > xmin .and. n > xnmin .and. .not.(i==3 .and. ltem1>=273.15)) then

    if (i >= 3 .and. ltem1 < 273.15) then
      qs=cal_qsati(tem,pp)
      dqs=cal_dqsatidt(tem,pp)
      es=cal_esati(tem)
      fl=cal_fls(tem)
    else
      qs=cal_qsat(tem,pp)
      dqs=cal_dqsatdt(tem,pp)
      es=cal_esat(tem)
      fl=cal_flv(tem)
    endif
    sat=qv/qs - 1.
  
    lambda=(cm(i)*n*cal_gamma(nu(i)+am(i)+1.)/(den*cal_gamma(nu(i)+1.)*q))**(1./(am(i)))
    lambda=min(max(lambda,10.),100000000.)
    dif = cal_xfkd(tem,pp)
    vis = cal_xfmu(tem)
    
    call ventilation (i, vis, lambda, den, fv)
    
    cd  = qs * ((c_p_v-c_v_v)*tem/(dif*es) + (fl/(kt*tem))*(fl/((c_p_v-c_v_v)*tem) - 1.))
    dav = cal_gamma(nu(i)+2.)*lambda**(-nu(i)-2.)
    cd  = 2.*pi*n*fv*dav/cd
  
    if ( cd > 1.e-9 .and. lq_1 > xmin ) then
      gamm = cd*(1. + (sat+1.)*dqs*fl/cp)
      beta = 1.

      ! Implicit formulation
      dq = max(min(dt*cd*(qv - qs) / (1. + beta*dt*gamm),lq_v1),-lq_1)

      if (lq_1 > xmin .and. n > xnmin) then    
        lq(5) = lq(5) + dq/dt
        ln(3) = ln(3) + min(n*dq/lq_1,0.)
      endif
   
      lq_v(5) = lq_v(5) - lq(5)
      ltem(5) = ltem(5) + lq(5)*fl/cp
    endif
    
  endif
  
end subroutine

subroutine ventilation (i, vis, lambda, dens, fv)
  implicit none
  integer :: i
  real  :: fv
  real  :: dens, lambda, vis
  real, parameter :: avent=0.78, bvent=0.308
  
  if (i /= 3 .and. i /= 1) then
    fv = bvent*0.84343*sqrt(ctv(i)*dens/vis)*(rho0/dens)**0.25    ! 0.84343=Sc**(1/3)
    fv = avent + fv*cal_gamma(nu(i)+2.+(btv(i)+1.)/2.)/cal_gamma(nu(i)+2.)*lambda**(-(btv(i)+1.)/2.)
  else
    fv = 1.
  endif
  
end subroutine

subroutine melt_ice (dt, cp, lq_i, ln_i, lq_c, ltem)
  implicit none
  real  :: dt, cp, tem
  real, dimension(5) :: ltem,lq_c,lq_i,ln_i

  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)

  if (lq_i1 > xmin .and. ltem1 >= 273.15) then
    ln_i(3) = 0.
    lq_i(5) = lq_i(5) - max(lq_i1/dt,0.)
    lq_c(5) = lq_c(5) + max(lq_i1/dt,0.)
    ltem(5) = ltem(5) - max(lq_i1/dt,0.) * cal_flm(tem)/cp
  endif
  
end subroutine

subroutine melt (dt, i, cp, lq_l, lq, lq_1, ln, lrho, ltem)
  implicit none
  integer :: i
  real  :: dt, cme, fv, vis, lambda
  real  :: dq, q, n, cp, tem, den, lq_1
  real, dimension(5) :: lq_l,lq,ln,ltem,lrho

  den  = lrho(3)
  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  q   = relaxation*lq_1+(1.-relaxation)*lq(3)
  n   = ln(3)

  if (lq_1 > xmin .and. n > xnmin .and. ltem1 >= 273.15) then
  
    lambda=(cm(i)*n*cal_gamma(nu(i)+am(i)+1.)/(den*cal_gamma(nu(i)+1.)*q))**(1./(am(i)))
    vis=cal_xfmu(tem)
    
    call ventilation (i,vis,lambda,den,fv)
    
    cme = -2.*pi*kt*(tem - 273.15)/cal_flm(tem)
    dq  = cme*n*fv*cal_gamma(nu(i)+2.)*lambda**(-nu(i)-2.)
    
    lq(5)   = lq(5)   - min(dq,lq_1/dt)
    lq_l(5) = lq_l(5) + min(dq,lq_1/dt)
    ltem(5) = ltem(5) - min(dq,lq_1/dt) * cal_flm(tem)/cp
    
  endif
  
end subroutine

subroutine sedim (dt, i, lq1, lq, ln, lrho)
  implicit none
  integer :: i
  real :: q, n, den, vel1, vel2, lambda
  
  real :: dt, lq1
  real, dimension(5) :: lq,ln,lrho
  
  den = lrho(3)
  q   = relaxation*lq1+(1.-relaxation)*lq(3)
  n   = ln(3)

  if (q > xmin .and. n > xnmin) then
    lambda=(cm(i)*n*cal_gamma(nu(i)+am(i)+1.)/(cal_gamma(nu(i)+1.)*q))**(1./(am(i)))
    lambda=max(10.,min(100000000.,lambda))
    
    vel1 = cal_gamma(am(i)+btv(i)+nu(i)+1.)/cal_gamma(am(i)+nu(i)+1.)
    vel2 = cal_gamma(btv(i)+nu(i)+1.)/cal_gamma(nu(i)+1.)
    
    lq(4) = vel1*ctv(i)*lambda**(-btv(i)) * sqrt(rho0/den)
    ln(4) = vel2*ctv(i)*lambda**(-btv(i)) * sqrt(rho0/den)
  endif  
  
end subroutine

subroutine assemble_micro (dt, cp, lqv, lqc, lqr, lqi, lqg, lqs, lnc, lni, ltem)

  integer :: i
  real  :: dt, cp
  real  :: fqc, fqr, tem
  real, dimension(5) :: lqv,lqr,lqc,lqi,lqg,lqs,lnc,lni,ltem

  tem = relaxation*ltem1+(1.-relaxation)*ltem(3)
  if (lq_v1 < xmin) then
    fqc = lq_c1/(lq_c1+lq_r1)
    fqr = lq_r1/(lq_c1+lq_r1)
    lq_c1 = lq_c1 + max(fqc*lq_v1,0.0)
    lq_r1 = lq_r1 + max(fqr*lq_v1,0.0)
    lq_v1 = 0.0
  endif
  
  lq_v(5) = (max(lq_v1,0.) - lq_v0)/dt
  
  lq_c(5) = (max(lq_c1,0.) - lq_c0)/dt  

  lq_r(5) = (max(lq_r1,0.) - lq_r0)/dt

  lq_i(5) = (max(lq_i1,0.) - lq_i0)/dt

  lq_g(5) = (max(lq_g1,0.) - lq_g0)/dt

  lq_s(5) = (max(lq_s1,0.) - lq_s0)/dt

  ln_c(5) = (max(ln_c1,0.) - ln_c0)/dt

  ln_i(5) = (max(ln_i1,0.) - ln_i0)/dt

  ltem(5) = (lq_c(5)+lq_r(5))*cal_flv(tem)/cp + (lq_i(5)+lq_g(5)+lq_s(5))*cal_fls(tem)/cp
  
end subroutine

end subroutine microphysics_1mom

subroutine calculate_diagnostic_microphysics(state,current_time,dt,init)

    type(state_type), intent(inout) :: state(:)
    type(mesh_type), pointer :: mesh
    real, intent(in) :: current_time
    real, intent(in) :: dt
    logical, intent(in), optional :: init

    type(scalar_field), pointer :: n
    character(len=OPTION_PATH_LEN) :: mom_path, sat_path
    character(len=FIELD_NAME_LEN) :: mesh_name, mom_name, mic_name
    logical :: one_mom, linit
    integer :: i, j
    
    if (present(init)) then
      linit=init
    else
      linit=.false.
    endif
    
    do j =1, size(state)	!<-------------------------------

    mom_path = '/material_phase['//int2str(j-1)//']/cloud_microphysics'
    if (have_option(trim(mom_path)//'/fortran_microphysics')) then
      one_mom = have_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics')
      
      if (have_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics')) then
        call get_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics/name',mic_name)
	mom_path=trim(mom_path)//'/fortran_microphysics/two_moment_microphysics::'//trim(mic_name)
      else if (have_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics')) then
        call get_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics/name',mic_name)
	mom_path=trim(mom_path)//'/fortran_microphysics/one_moment_microphysics::'//trim(mic_name)
      endif
    
    else if (have_option(trim(mom_path)//'/diagnostic')) then
      one_mom=.false. 
      mom_path=trim(mom_path)//'/diagnostic'
    endif
    
    sat_path = '/material_phase['//int2str(j-1)//']/cloud_microphysics/saturation_adjustment'
    sat_adj = have_option(trim(sat_path))
    drop_adj = have_option('/material_phase[0]/cloud_microphysics/fortran_microphysics/two_moment_microphysics::'//trim(mic_name)//'/droplet_adjustment')
 
    ! Reference mesh is the one used for microphysics variables
    if (have_option(trim(mom_path)//'/scalar_field::Qdrop/prognostic'))then
    	call get_option(trim(mom_path)//'/scalar_field::Qdrop/prognostic/mesh/name',mesh_name)
    else if (have_option(trim(mom_path)//'/scalar_field::Qdrop/diagnostic'))then
    	call get_option(trim(mom_path)//'/scalar_field::Qdrop/diagnostic/mesh/name',mesh_name)
    else if (have_option(trim(mom_path)//'/scalar_field::Qdrop/prescribed'))then
    	call get_option(trim(mom_path)//'/scalar_field::Qdrop/prescribed/mesh/name',mesh_name)
    endif  
    
    mesh=>extract_mesh(state(j), trim(mesh_name)) 
       
    ! Calculate saturation adjustment
    if (sat_adj .or. drop_adj .or. linit) then 
      call saturation_adjustment (state(j), mesh, sat_path, dt, linit)
    end if
    
    ! One moment scheme: diagnose number concentrations
    if (one_mom) then
      if (has_scalar_field(state(j),"Nrain")) then
        n=>extract_scalar_field(state(j),"Nrain")
        call calculate_n_one(state(j),"rain",n)
      end if
      if (has_scalar_field(state(j),"Ngrau")) then
        n=>extract_scalar_field(state(j),"Ngrau")
        call calculate_n_one(state(j),"grau",n)
      end if
      if (has_scalar_field(state(j),"Qsnow")) then
        n=>extract_scalar_field(state(j),"Nsnow")
        call calculate_n_one(state(j),"snow",n)
      end if
    endif
    
    enddo  	!<----------------------------------

end subroutine calculate_diagnostic_microphysics
  
subroutine saturation_adjustment(state, mesh, sat_path, dt, linit)

    type(state_type), intent(inout) :: state
    type(mesh_type), pointer, intent(inout) :: mesh
    real, intent(in) :: dt
    character(len=*), intent(in) :: sat_path
    logical, intent(in) :: linit

    type(scalar_field), pointer :: scal,thermal,density,pressure,sat,q_c,q_r,q_i,q_g,q_s,q_v,fraction
    type(scalar_field) :: pressure_local,temperature,eosdensity,qc_p,qc_v
    type(scalar_field) :: thermal_old,q_c_old,epsilon
    type(scalar_field), target :: dummyscalar
    real :: thermal_node, saturation, exn, qv_node, qc_node, tem, therm, pp, cp_node, cv_node
    real :: dq_c, d_therm, dtem, qc_new, qv_new, thermal_new, tol
    real :: c_p, c_v, c_p_v, c_v_v, c_p_l, c_v_l, c_p_i, c_v_i
    character(len=OPTION_PATH_LEN) :: eos_path, phase_name
    integer :: i,node,stat,thermal_variable, max_iter
    logical :: have_qv, have_ql, have_qi, constant_cp_cv

    ewrite(1,*) 'Entering saturation adjustment'
    
    have_qv=.false.; have_ql=.false.; have_qi=.false.
    call get_option(trim(sat_path)//"/tolerance",tol,default=1.e-6)
    call get_option(trim(sat_path)//"/max_iterations",max_iter,default=20)
    constant_cp_cv=have_option('/material_phase[0]/equation_of_state/compressible/giraldo/constant_cp_cv')
    
    call allocate (dummyscalar,mesh,"DummyScalar")
    call zero(dummyscalar)
        
    call get_thermo_variable(state, thermal, index=thermal_variable)
    assert(thermal%mesh==mesh)
      
    eos_path = trim(state%option_path)//'/equation_of_state'
    if(have_option(trim(eos_path)//'/compressible')) then
    
      call get_cp_cv (c_p,c_v,c_p_v=c_p_v,c_v_v=c_v_v, &
                   c_p_l=c_p_l,c_v_l=c_v_l,c_p_i=c_p_i,c_v_i=c_v_i)
    
      call allocate(eosdensity,mesh,"EOSDensity")
      call allocate(temperature,mesh,"Temperature")
      call allocate(pressure_local,mesh,"PressureLocal")
      call allocate(qc_p,mesh,"LocalCp")
      call allocate(qc_v,mesh,"LocalCv")
      call allocate(epsilon,mesh,"Epsilon")
	          
      if(have_option(trim(eos_path)//'/compressible/giraldo')) then
      
        if (has_scalar_field(state,"TotalWaterQ")) then
          q_v=>extract_scalar_field(state,"TotalWaterQ")
	  have_qv=.true.
        else if (has_scalar_field(state,"VapourWaterQ")) then
          q_v=>extract_scalar_field(state,"VapourWaterQ")
	  have_qv=.true.
        else
          q_v=>dummyscalar
        end if
        if (has_scalar_field(state,"Qdrop")) then
          q_c=>extract_scalar_field(state,"Qdrop")
	  have_ql=.true.
        else
          q_c=>dummyscalar
        end if
        if (has_scalar_field(state,"Qrain")) then
          q_r=>extract_scalar_field(state,"Qrain")
        else
          q_r=>dummyscalar
        end if
        if (has_scalar_field(state,"Qice")) then
          q_i=>extract_scalar_field(state,"Qice")
	  have_qi=.true.
        else
          q_i=>dummyscalar
        end if
        if (has_scalar_field(state,"Qgrau")) then
          q_g=>extract_scalar_field(state,"Qgrau")
        else
          q_g=>dummyscalar
        end if
        if (has_scalar_field(state,"Qsnow")) then
          q_s=>extract_scalar_field(state,"Qsnow")
        else
          q_s=>dummyscalar
        end if
	  
        if (has_scalar_field(state,"TotalWaterQ") .or. has_scalar_field(state,"VapourWaterQ")) then
           assert(q_v%mesh==mesh)
        end if
        if (has_scalar_field(state,"Qdrop")) then
           assert(q_c%mesh==mesh)
        end if
        if (has_scalar_field(state,"Qrain")) then
           assert(q_r%mesh==mesh)
        end if
        if (has_scalar_field(state,"Qice")) then
           assert(q_i%mesh==mesh)
        end if
        if (has_scalar_field(state,"Qgrau")) then
           assert(q_g%mesh==mesh)
        end if
        if (has_scalar_field(state,"Qsnow")) then
           assert(q_s%mesh==mesh)
        end if
	
	pressure=>extract_scalar_field(state,"Pressure")
	call safe_set(state,pressure_local,pressure)
	
      end if
       
      call compressible_eos(state,density=eosdensity,temperature=temperature)
	
      i=0
      call set(epsilon,1.)
      do while (max(abs(maxval(epsilon%val)),abs(minval(epsilon%val))) > tol .and. i < max_iter)
        i=i+1
	
        call make_giraldo_quantities_1mat(state,have_qv,have_ql,have_qi,q_v,q_c,q_r,q_i,q_g,q_s,qc_p,qc_v)

	do node=1,node_count(q_c)
          tem=node_val(temperature,node)
	  therm=node_val(thermal,node)
	  pp=node_val(pressure_local,node)
	  saturation=cal_qsat(tem,pp)
	  qc_node=node_val(q_c,node)
	  qv_node=node_val(q_v,node)
	  cp_node=node_val(qc_p,node)
	  cv_node=node_val(qc_v,node)
	  if (has_scalar_field(state,"TotalWaterQ")) qv_node=qv_node-node_val(q_c,node)-node_val(q_r,node) &
	  					            -node_val(q_i,node)-node_val(q_g,node)-node_val(q_s,node)

          dq_c=0.
	  d_therm=0.
	  if ( qv_node > saturation .or. qc_node > xmin ) then
  	    dq_c=(qv_node-saturation) / (1. + cal_flv(tem)**2.*saturation/(node_val(qc_p,node)*   &
	         (node_val(qc_p,node)-node_val(qc_v,node)) * tem**2.))
		 
	    qc_new=min(max(qc_node+dq_c,0.),qv_node)
	    dq_c=qc_new-qc_node
	    qv_new=min(max(qv_node-dq_c,0.),qv_node+qc_node)
   
            dtem=cal_flv(tem)/node_val(qc_p,node)
	    if (thermal_variable == 1) then
	      if (.not.constant_cp_cv) then
	        d_therm = node_val(qc_v,node)*dtem - tem*(c_p_v-c_p_l)
	      else
	        d_therm = node_val(qc_v,node)*dtem
	      endif
	    else if (thermal_variable == 2) then
	      d_therm = dtem
            else if (thermal_variable == 3) then
	      exn = (pp/1.E+05)**((cp_node-cv_node)/cp_node)
              d_therm = dtem * therm/tem
	      if (.not.constant_cp_cv) then
                d_therm = d_therm + therm*log(exn)*((c_p_v-c_v_v)/(node_val(qc_p,node)-node_val(qc_v,node)) &
		        - (c_p_v-c_p_l)/node_val(qc_p,node))
	      endif
	    endif
	    thermal_new=therm+d_therm*dq_c
	    tem=tem+dtem*dq_c
	    
            call set (thermal,node,thermal_new)
	    call set (temperature,node,tem)
            call set (q_c,node,qc_new)
	    if (.not.has_scalar_field(state,"TotalWaterQ")) call set (q_v,node,qv_new)
	  endif
          call set (epsilon,node,d_therm*dq_c)
	  
	enddo
      enddo
      
      ! Update density
      density=>extract_scalar_field(state,'Density')
      call compressible_eos(state, density=density)
 
      call deallocate (eosdensity)
      call deallocate (temperature)
      call deallocate (pressure_local)
      call deallocate (qc_p)
      call deallocate (qc_v)
      call deallocate (epsilon)
      call deallocate (dummyscalar)
     
    end if
      
contains
    
    subroutine store_result(field,vloc1,vloc2,dt)
      type(microphysics_field), intent(inout) :: field
      type(scalar_field) :: vloc1,vloc2
      real, intent(in) :: dt

      if (field%has_forcing) then
         call addto (vloc2,vloc1,scale=-1.)
	 call scale (vloc2,-1./dt)
         call set(field%forcing,vloc2)	 
      endif
    end subroutine store_result

    subroutine clean_up(field)
      type(microphysics_field), intent(inout) :: field
      integer :: i

      do i=1,size(field%data)
         call deallocate(field%data(i))
      end do

      nullify(field%forcing)
      nullify(field%sinking_velocity)
    end subroutine clean_up

end subroutine saturation_adjustment

subroutine limit_microphysics(state)

    type(state_type), intent(inout) :: state(:)
    integer :: j
    
    ewrite(2,*) 'In limit_microphysics'
    
    do j = 1,size(state)
    
       call limit(state(j), "drop")
       call limit(state(j), "rain")
       call limit(state(j), "TotalWaterQ")
       call limit(state(j), "VapourWaterQ")
    
    enddo

contains

  subroutine limit(state,field_name)
    type(state_type), intent(in) :: state
    character(len=*), intent(in) :: field_name
    type(scalar_field), pointer :: sfield, nfield, qfield
    real :: delta
    integer :: n,stat,stat0
    logical :: have_n=.false., have_q=.false., have_s=.false.
    
    sfield=>extract_scalar_field(state,trim(field_name),stat0)
    if (stat0 /= 0) then
      sfield=>extract_scalar_field(state,'VapourWaterQ',stat)
      if (stat==0 .and. have_option(trim(sfield%option_path)//'/prognostic')) have_s=.true.
      nfield=>extract_scalar_field(state,'N'//trim(field_name),stat)
      if (stat==0 .and. have_option(trim(nfield%option_path)//'/prognostic')) have_n=.true.
      qfield=>extract_scalar_field(state,'Q'//trim(field_name),stat)
      if (stat==0 .and. have_option(trim(qfield%option_path)//'/prognostic')) have_q=.true.
    else
      if (have_option(trim(sfield%option_path)//'/prognostic')) have_s=.true.
    endif
    
    if (stat0==0 .and. have_s) then

      do n = 1, size(sfield%val)
        sfield%val(n)=max(sfield%val(n),0.0)
      enddo
    
    else if (stat==0 .and. (have_n.or.have_q)) then
    
      do n = 1, size(qfield%val)
        
	if (qfield%val(n)<xmin .or. nfield%val(n)<xnmin) then
!	  if (have_s .and. have_q) sfield%val(n)=sfield%val(n)+qfield%val(n)
	  if (have_n) nfield%val(n)=0.
	  if (have_q) qfield%val(n)=0.
	endif
	
      enddo      
    
    endif
    
  end subroutine limit

end subroutine limit_microphysics

subroutine extract_and_project(lstate,lmesh,mfield,fname,local,entropy)

    type(state_type), intent(inout) :: lstate
    type(mesh_type),intent(inout), pointer :: lmesh
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: X
    character(len=*), intent(in) :: fname
    type(microphysics_field) :: mfield
    type(scalar_field), target, optional, intent(in) :: local
    type(scalar_field), target :: dummyscalar
     
    logical, intent(in), optional :: entropy
    logical :: lentropy
    
    integer :: i,stat
    character(len=*), dimension(3), parameter::&
    	 old=(/ "Old     ","Iterated","        " /)
	 
    ewrite(3,*) 'In extract_and_project'
	
    if (present(entropy)) then
       lentropy=entropy
    else
       lentropy=.false.
    end if

    X=>extract_vector_field(lstate,"Coordinate")

    call allocate(dummyscalar,lmesh,name="Dummy")
    call zero(dummyscalar)
    
    do i=1,size(mfield%data)
       call allocate(mfield%data(i),lmesh,name=trim(old(i))//trim(fname))

       if (lentropy) then
    	  sfield=>extract_entropy_variable(lstate,prefix=trim(old(i)))
       else
          if (has_scalar_field(lstate,trim(fname))) then
    	    sfield=>extract_scalar_field(lstate,trim(old(i))//trim(fname))	    
          else if (present(local)) then      
            sfield=>local
          else
            sfield=>dummyscalar
          endif
       endif
       
       call project_field(sfield,mfield%data(i),X)
    end do

    mfield%has_forcing=.false.
    mfield%has_sinking_velocity=.false.

    if (lentropy) then
       mfield%forcing=>extract_entropy_variable(lstate,&
    	    suffix="MicrophysicsSource",stat=stat)
       if (stat==0) mfield%has_forcing=.true.       
       mfield%sinking_velocity=>extract_entropy_variable(lstate,&
    	    suffix="SinkingVelocity",stat=stat)
       if (stat==0) mfield%has_sinking_velocity=.true.
    else
       if (present(local) .and. trim(fname)=="Temperature") then 
         mfield%forcing=>extract_entropy_variable(lstate,&
    	      suffix="MicrophysicsSource",stat=stat)
         if (stat==0) mfield%has_forcing=.true.       
!         mfield%has_sinking_velocity=.true.
       else
         mfield%forcing=>extract_scalar_field(lstate,&
    	      trim(fname)//"MicrophysicsSource",stat=stat)
         if (stat==0) mfield%has_forcing=.true.
         mfield%sinking_velocity=>extract_scalar_field(lstate,&
    	      trim(fname)//"SinkingVelocity",stat=stat)
         if (stat==0) mfield%has_sinking_velocity=.true.
      endif
    end if
    call deallocate(dummyscalar)
        
end subroutine extract_and_project

subroutine initialise_microphysics(state,current_time,dt)
  type(state_type), intent(inout) :: state(:)
  type(mesh_type), pointer :: mesh
  type(scalar_field), pointer :: gas_density
  real, intent(in) :: current_time
  real, intent(in) :: dt
  character(len=OPTION_PATH_LEN) :: mom_path
  character(len=FIELD_NAME_LEN) :: mesh_name, mic_name
  
  logical :: microphysics_on
  integer :: i, j

  ! In case it's needed later, so we can ensure everything we need exists
  
  do j= 1, size(state)
  
  microphysics_on=have_option('/material_phase::'//trim(state(j)%name)//'/cloud_microphysics/fortran_microphysics')

  if (microphysics_on) then
      mom_path = '/material_phase['//int2str(j-1)//']/cloud_microphysics'
      if (have_option(trim(mom_path)//'/fortran_microphysics')) then
        if (have_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics')) then
          call get_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics/name',mic_name)
	  mom_path=trim(mom_path)//'/fortran_microphysics/two_moment_microphysics::'//trim(mic_name)
	else if (have_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics')) then
          call get_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics/name',mic_name)
	  mom_path=trim(mom_path)//'/fortran_microphysics/one_moment_microphysics::'//trim(mic_name)
	endif
      endif
 
    if (have_option(trim(mom_path)//'/scalar_field::Qdrop/prognostic'))then
   	call get_option(trim(mom_path)//'/scalar_field::Qdrop/prognostic/mesh/name',mesh_name)    
          else if (have_option(trim(mom_path)//'/scalar_field::Qdrop/diagnostic'))then
        call get_option(trim(mom_path)//'/scalar_field::Qdrop/diagnostic/mesh/name',mesh_name)	
          else if (have_option(trim(mom_path)//'/scalar_field::Qdrop/prescribed'))then
        call get_option(trim(mom_path)//'/scalar_field::Qdrop/prescribed/mesh/name',mesh_name)	
    endif  
    
    mesh=>extract_mesh(state(j), trim(mesh_name))
  end if
  
  enddo

end subroutine initialise_microphysics

!
!  Useful microphysics functions 
!

subroutine calculate_lambda_two(state, hydro_name, lambda)

  type(state_type), intent(inout) :: state
  character(len=*), intent(in) :: hydro_name
  type(scalar_field), pointer, intent(inout) :: lambda
  
  type(mesh_type), pointer :: mesh
  type(scalar_field), pointer :: q, n, density
  type(scalar_field), target :: dummyscalar
  type(scalar_field) :: density_remap
  integer :: h, node
  real :: lambda_node
  real :: param, eps=1.e-12
 
  mesh=>lambda%mesh
  call allocate (dummyscalar,mesh,"DummyScalar")
  call zero(dummyscalar)
  
  density=>extract_scalar_field(state,"Density")
  call allocate(density_remap,mesh,"DensityRemap")
  call safe_set(state,density_remap,density)
  
  if (has_scalar_field(state,'Q'//trim(hydro_name))) then
    q=>extract_scalar_field(state,'Q'//trim(hydro_name))
  else
    q=>dummyscalar
  end if
  if (has_scalar_field(state,'N'//trim(hydro_name))) then
    n=>extract_scalar_field(state,'N'//trim(hydro_name))
  else
    n=>dummyscalar
  end if
  
  assert(mesh==q%mesh)
  assert(mesh==n%mesh)
  
  if (trim(hydro_name) == "rain") then
    h = 2
  else if (trim(hydro_name) == "ice") then
    h = 3
  else if (trim(hydro_name) == "grau") then
    h = 4
  else if (trim(hydro_name) == "snow") then
    h = 5
  endif
  
  param = cm(h)*cal_gamma(am(h)+nu(h)+1.)/cal_gamma(nu(h)+1.)
  
  call zero(lambda)
  
  do node = 1, node_count(q)
    if (node_val(q,node) > xmin .and. node_val(n,node) > xnmin) then
      lambda_node=(param*node_val(n,node)/(node_val(density_remap,node)*node_val(q,node)))**(1./am(h))
      call set(lambda,node,lambda_node)
    endif
  enddo
  
  call deallocate (dummyscalar)
  call deallocate (density_remap)
  
end subroutine calculate_lambda_two

subroutine calculate_lambda_one(state, hydro_name, lambda)

  type(state_type), intent(inout) :: state
  character(len=*), intent(in) :: hydro_name
  type(scalar_field), pointer, intent(inout) :: lambda
  
  type(mesh_type), pointer :: mesh
  type(scalar_field), pointer :: q, density
  type(scalar_field), target :: dummyscalar 
  type(scalar_field) :: density_remap
  integer :: h, node
  real :: lambda_node
  real :: N0, param
 
  mesh=>lambda%mesh
  call allocate (dummyscalar,mesh,"DummyScalar")
  call set(dummyscalar,0.0)
  
  density=>extract_scalar_field(state,"Density")
  call allocate(density_remap,mesh,"DensityRemap")
  call safe_set(state,density_remap,density)
  
  if (has_scalar_field(state,'Q'//trim(hydro_name))) then
    q=>extract_scalar_field(state,'Q'//trim(hydro_name))
  else
    q=>dummyscalar
  end if
  
  if (trim(hydro_name) == "rain") then
    h = 2
    N0=N0r
  else if (trim(hydro_name) == "grau") then
    h = 4
    N0=N0g
  else if (trim(hydro_name) == "snow") then
    h = 5
    N0=N0s
  endif
  
  assert(mesh==q%mesh)
  
  param = cm(h)*N0*cal_gamma(am(h)+nu(h)+1.)
  
  call zero(lambda)
  
  do node = 1, node_count(q)
    if (node_val(q,node) > xmin) then
      lambda_node=(param/(node_val(density_remap,node)*node_val(q,node)))**(1./(am(h)+nu(h)+1.))
      call set(lambda,node,lambda_node)
    endif
  enddo
  
  call deallocate (dummyscalar)
  call deallocate (density_remap)
  
end subroutine calculate_lambda_one

subroutine calculate_n_one(state, hydro_name, n)

  type(state_type), intent(inout) :: state
  character(len=*), intent(in) :: hydro_name
  type(scalar_field), pointer, intent(inout) :: n
  
  type(mesh_type), pointer :: mesh
  type(scalar_field), pointer :: lambda
  type(scalar_field), target :: dummyscalar
  integer :: h
  real :: N0, param
 
  mesh=>n%mesh
  call allocate (dummyscalar,mesh,"DummyScalar")
  call set(dummyscalar,0.0)
  
  lambda=>dummyscalar
  call calculate_lambda_one(state,trim(hydro_name),lambda)
  
  if (trim(hydro_name) == "rain") then
    h = 2
    N0=N0r
  else if (trim(hydro_name) == "grau") then
    h = 4
    N0=N0g
  else if (trim(hydro_name) == "snow") then
    h = 5
    N0=N0s
  endif
  
  param = N0*cal_gamma(nu(h)+1.)
  
  call set(n,lambda)
  where(n%val > xnmin) n%val=n%val**(-(nu(h)+1.))
  call scale(n,param)
  
  call deallocate(dummyscalar)
  
end subroutine calculate_n_one

function cal_xfkd (tem,p)
  implicit none
  real :: tem,p
  real :: cal_xfkd
  
  cal_xfkd = 10./p * (1.53e-3*tem - 0.192)
  
  return
end function

function cal_xfmu (tem)
  implicit none
  real :: tem
  real :: cal_xfmu
  
  cal_xfmu = 1.72e-5*(392./(tem + 120.)) * (tem/273.)**(1.5)
  
  return
end function

real function cal_sigmawv( T )
  implicit none
  real	:: T, Tc

  Tc = T - 273.15
  cal_sigmawv = (75.93 + Tc*(0.115 + Tc*(6.818e-2 + Tc*(6.511e-3 + Tc*(2.933e-4 + Tc*(6.283e-6 + Tc*5.285e-8)))))) * 1.e-3

  return                                                                   
end function

function cal_gamma (x)
  implicit none
  real :: x, xx, cal_gamma
  real, parameter :: p0=1.000000000190015,	&
  		     p1=76.18009172947146,	&
		     p2=-86.50532032941677,	&
		     p3=24.01409824083091,	&
		     p4=-1.231739572450155,	&
		     p5=1.208650973866179e-3,	&
		     p6=-5.395239384953e-6
		     
  xx = p0 + p1/(x+1.) + p2/(x+2.) + p3/(x+3.) + p4/(x+4.) + p5/(x+5.) + p6/(x+6.)
  cal_gamma = sqrt(2.*pi) * xx/x * (x+5.5)**(x+0.5) * exp(-(x+5.5))
  return
end function

real function erreur(tt, qv, qc, qr, tt_1, qv_1, qc_1, qr_1, p)
  
  implicit none
  real :: tt, qv, qc, qr, tt_1, qv_1, qc_1, qr_1
  real :: err(4)
  integer :: p
  
  err(1) = real(fac(p))/(real(p)+1) * abs((tt-tt_1)/tt)
  err(2) = real(fac(p))/(real(p)+1) * abs((qv-qv_1)/qv)
  err(3) = real(fac(p))/(real(p)+1) * abs((qc-qc_1)/qc)
  err(4) = real(fac(p))/(real(p)+1) * abs((qr-qr_1)/qr)

  erreur=(rtol/maxval(err))**(1./real(p+1))

  return
end function

integer function fac( p )
  integer :: i, p
  
  fac=1
  do i = 1, p  
    fac=fac*p
  enddo
  
  return
end function

function calerf ( arg, jint )

!*****************************************************************************80
!
!! CALERF computes various forms of the error function.
!
!  Discussion:
!
!    This routine evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
!    for a real argument x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev Approximations for the Error Function,
!    Mathematics of Computation,
!    Volume 23, Number 107, July 1969, pages 631-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT is 1, the
!    argument must be less than XBIG.  If JINT is 2, the argument
!    must lie between XNEG and XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = erf(x);
!    1, RESULT = erfc(x) = 1 - erf(x);
!    2, RESULT = exp(x*x)*erfc(x) = exp(x*x) - erf(x*x)*erf(x).
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, erf(x);
!    1, erfc(x);
!    2, exp(x*x)*erfc(x).
!
  implicit none

  real a(5)
  real arg
  real b(4)
  real c(9)
  real d(8)
  real del
  integer i
  integer jint
  real p(6)
  real q(5)
  real result
  real sixten
  real sqrpi
  real thresh
  real x
  real xbig
  real xden
  real xhuge
  real xinf
  real xmax
  real xneg
  real xnum
  real xsmall
  real y
  real ysq
  real calerf
!
!  Mathematical constants
!
  data sqrpi / 5.6418958354775628695d-1 /
  data thresh / 0.46875d0 /
  data sixten / 16.0d0 /
!
!  Machine-dependent constants
!
  data xinf /1.79d308 /
  data xneg / -26.628d0 /
  data xsmall /1.11d-16/
  data xbig /26.543d0 /
  data xhuge /6.71d7/
  data xmax /2.53d307/
!
!  Coefficients for approximation to  erf  in first interval
!
  data a/3.16112374387056560d00,1.13864154151050156d02, &
         3.77485237685302021d02,3.20937758913846947d03, &
         1.85777706184603153d-1/
  data b/2.36012909523441209d01,2.44024637934444173d02, &
         1.28261652607737228d03,2.84423683343917062d03/
!
!  Coefficients for approximation to  erfc  in second interval
!
  data c/5.64188496988670089d-1,8.88314979438837594d0, &
         6.61191906371416295d01,2.98635138197400131d02, &
         8.81952221241769090d02,1.71204761263407058d03, &
         2.05107837782607147d03,1.23033935479799725d03, &
         2.15311535474403846d-8/
  data d/1.57449261107098347d01,1.17693950891312499d02, &
         5.37181101862009858d02,1.62138957456669019d03, &
         3.29079923573345963d03,4.36261909014324716d03, &
         3.43936767414372164d03,1.23033935480374942d03/
!
!  Coefficients for approximation to  erfc  in third interval
!
  data p/3.05326634961232344d-1,3.60344899949804439d-1, &
         1.25781726111229246d-1,1.60837851487422766d-2, &
         6.58749161529837803d-4,1.63153871373020978d-2/
  data q/2.56852019228982242d00,1.87295284992346047d00, &
         5.27905102951428412d-1,6.05183413124413191d-2, &
         2.33520497626869185d-3/

  x = arg
  y = abs ( x )
!
!  Evaluate erf for |X| <= 0.46875.
!
  if ( y <= thresh ) then

    result = -1.0D+00

    ysq = 0.0D+00
    if ( xsmall < y ) then
      ysq = y * y
    end if

    xnum = a(5) * ysq
    xden = ysq

    do i = 1, 3
      xnum = ( xnum + a(i) ) * ysq
      xden = ( xden + b(i) ) * ysq
    end do

    result = x * ( xnum + a(4) ) / ( xden + b(4) )

    if ( jint /= 0 ) then
      result = 1.0D+00 - result
    end if

    if ( jint == 2 ) then
      result = exp ( ysq ) * result
    end if

    calerf = result

    return
!
!  Evaluate erfc for 0.46875 <= |X| <= 4.0.
!
   else if ( y <= 4.0D+00 ) then

     xnum = c(9) * y
     xden = y

     do i = 1, 7
       xnum = ( xnum + c(i) ) * y
       xden = ( xden + d(i) ) * y
     end do

     result = ( xnum + c(8) ) / ( xden + d(8) )

     if ( jint /= 2 ) then
       ysq = aint ( y * sixten ) / sixten
       del = ( y - ysq ) * ( y + ysq )
       result = exp ( -ysq * ysq ) * exp ( -del ) * result
     end if
!
!  Evaluate erfc for 4.0 < |X|.
!
   else

     result = 0.0D+00

     if ( xbig <= y ) then

       if ( jint /= 2 .or. xmax <= y ) then
         go to 300
       end if

       if ( xhuge <= y ) then
         result = sqrpi / y
         go to 300
       end if

     end if

     ysq = 1.0D+00 / ( y * y )
     xnum = p(6) * ysq
     xden = ysq
     do i = 1, 4
       xnum = ( xnum + p(i) ) * ysq
       xden = ( xden + q(i) ) * ysq
      end do

      result = ysq * ( xnum + p(5) ) / ( xden + q(5) )
      result = ( sqrpi -  result ) / y

      if ( jint /= 2 ) then
        ysq = aint ( y * sixten ) / sixten
        del = ( y - ysq ) * ( y + ysq )
        result = exp ( -ysq * ysq ) * exp ( -del ) * result
      end if

  end if
!
!  Fix up for negative argument, erf, etc.
!
  300 continue

  if ( jint == 0 ) then

    result = ( 0.5D+00 - result ) + 0.5D+00
    if ( x < 0.0D+00 ) then
      result = -result
    end if

  else if ( jint == 1 ) then

    if ( x < 0.0D+00 ) then
      result = 2.0D+00 - result
    end if

  else

    if ( x < 0.0D+00 ) then

      if ( x < xneg ) then
        result = xinf
      else
        ysq = aint ( x * sixten ) / sixten
        del = ( x - ysq ) * ( x + ysq )
        y = exp ( ysq * ysq ) * exp ( del )
        result = ( y + y ) - result
      end if

    end if

  end if
  
  calerf = result

  return
end function

end module microphysics
