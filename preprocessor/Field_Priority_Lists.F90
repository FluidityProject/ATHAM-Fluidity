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

module field_priority_lists

  use fldebug
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use futils, only: int2str
  use spud
  use fields
  use state_module
  use sediment, only: get_n_sediment_fields, get_sediment_item

  implicit none
  
  !! Field name list for tracers (from 1 to NTSOL)
  character(len=FIELD_NAME_LEN), save, &
       dimension(:), allocatable :: field_name_list
  !! Field list for tracters (from 1 to NTSOL)
  type(scalar_field), dimension(:), allocatable, save :: field_list
  !! Options path list for tracers (from 1 to NTSOL)
  character(len=OPTION_PATH_LEN), save, &
       dimension(:), allocatable :: field_optionpath_list
  !! State list for tracers (from 1 to NTSOL)
  integer, save, dimension(:), allocatable :: field_state_list

  private
  public :: field_name_list, field_list, field_optionpath_list,&
       & field_state_list, initialise_field_lists_from_options,&
       & get_ntsol 
       
contains

  subroutine initialise_field_lists_from_options(state, ntsol)
    type(state_type), dimension(:), intent(in) :: state
    integer, intent(in) :: ntsol

    logical, save:: initialised=.false.
    integer :: nsol, nphases,nfields,ncars,p,f,i, tmpint
    character(len=FIELD_NAME_LEN) :: tmpstring
    logical :: aliased, pressure, surface

    integer, dimension(:), allocatable :: priority
    !! Field list for tracers (from 1 to NTSOL)
    character(len=FIELD_NAME_LEN), save, &
        dimension(:), allocatable :: temp_field_name_list
    !! Options path list for tracers (from 1 to NTSOL)
    character(len=OPTION_PATH_LEN), save, &
        dimension(:), allocatable :: temp_field_optionpath_list
    !! State list for tracers (from 1 to NTSOL)
    integer, save, dimension(:), allocatable :: temp_field_state_list
    character(len=FIELD_NAME_LEN):: micro_optionpath, micro_name
    
    ewrite(1,*) 'In initialise_field_lists_from_options'
    
    
    ! if called for the second time return immediately
    if (.not.initialised) then
    
       allocate( field_name_list(ntsol), &
            field_state_list(ntsol), &
            field_optionpath_list(ntsol),&
            priority(ntsol), &
            field_list(ntsol), &
            temp_field_name_list(ntsol), &
            temp_field_state_list(ntsol), &
            temp_field_optionpath_list(ntsol) )

       nsol = 0

       nphases = option_count('/material_phase')  
       do p = 0, nphases-1
          nfields = option_count('/material_phase[' &
               //int2str(p)//']/scalar_field')
          do f = 0,nfields-1
             aliased = have_option('/material_phase['// &
                  int2str(p)//']/scalar_field['//int2str(f)//']/aliased')
             surface = have_option('/material_phase['// &
                  int2str(p)//']/scalar_field['//int2str(f)//']/diagnostic/boundary_conditions')
             call get_option('/material_phase['// &
                  int2str(p)// &
                  ']/scalar_field['//int2str(f)//']/name', &
                  tmpstring)
             call get_option('/material_phase['// &
                  int2str(p)// &
                  ']/scalar_field['//int2str(f)//']/&
                  &prognostic/priority', &
                  tmpint, default=0)
             pressure = (trim(tmpstring)=='Pressure')

             if (.not. aliased .and. .not. pressure .and. .not. surface) then
                nsol = nsol + 1
                temp_field_name_list(nsol) = tmpstring
                temp_field_optionpath_list(nsol) = '/material_phase['// &
                     int2str(p)// &
                     ']/scalar_field::'//trim(tmpstring)
                temp_field_state_list(nsol) = p+1
                priority(nsol) = tmpint
             end if
          end do

          ! prognostic sediment fields
          if (have_option('/material_phase['//int2str(p)//']/sediment')) then
             nfields = get_n_sediment_fields()
             do f = 1, nfields
                nsol=nsol+1
                
                call get_sediment_item(state(p+1), f, temp_field_name_list(nsol))
                temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                     ']/sediment/scalar_field['//int2str(f-1)//']'
                temp_field_state_list(nsol) = p+1
                call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                     tmpint, default=0)
                priority(nsol) = tmpint
             end do
          end if

          ! prognostic Mellor Yamada fields:
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::KineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "KineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::KineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol) = tmpint
          end if
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::TurbulentLengthScalexKineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "TurbulentLengthScalexKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::TurbulentLengthScalexKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol) = tmpint
          end if

          ! Check for GLS - we need to make sure these fields are solved *after*
          ! everything else, so set to a big negative value. In addition, the
          ! Psi solve *must* come after the TKE solve, so make sure the priority
          ! is set such that this happens
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "GLSTurbulentKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "GLSGenericSecondQuantity"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if

          ! Check for k-epsilon - we need to make sure these fields are solved *after*
          ! everything else, so set to a big negative value. In addition, the
          ! TurbulentDissipation (Epsilon) solve *must* come after the TKE solve,
          ! so make sure the priority is set such that this happens.
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "TurbulentKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "TurbulentDissipation"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if
          ! Check for subgrid-scale kinetic energy equation
          ! - we need to make sure this is solved *after*
          ! everything else, so set to a big negative value.
          if(have_option('/material_phase['//int2str(p)// &
             ']/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/scalar_field::SubgridKineticEnergy')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "SubgridKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
             ']/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/scalar_field::SubgridKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
!!! Melt rate should be the last thing to calculate, Sb
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Sb/diagnostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "Sb"
             temp_field_optionpath_list(nsol)='/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Sb'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/diagnostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
!!! Melt rate should be the last thing to calculate, Tb
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Tb/diagnostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "Tb"
             temp_field_optionpath_list(nsol)='/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Tb'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/diagnostic/priority', &
             tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
!!!/ocean_forcing/iceshelf_meltrate/Holland08
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::MeltRate/diagnostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "MeltRate"
             temp_field_optionpath_list(nsol)='/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::MeltRate'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/diagnostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
	  
!!!/microphysics
       
         if (have_option('/material_phase['//int2str(p)//']/cloud_microphysics')) then
	 
	 micro_optionpath='/material_phase['//int2str(p)//']/cloud_microphysics'
         if (have_option(trim(micro_optionpath)//'/fortran_microphysics')) then
	   if (have_option(trim(micro_optionpath)//'/fortran_microphysics/two_moment_microphysics')) then
	     call get_option(trim(micro_optionpath)//'/fortran_microphysics/two_moment_microphysics/name', &
	                     micro_name)
	     micro_optionpath=trim(micro_optionpath)//'/fortran_microphysics/two_moment_microphysics::'//trim(micro_name)
	   else if (have_option(trim(micro_optionpath)//'/fortran_microphysics/one_moment_microphysics')) then
	     call get_option(trim(micro_optionpath)//'/fortran_microphysics/one_moment_microphysics/name', &
	                     micro_name)
	     micro_optionpath=trim(micro_optionpath)//'/fortran_microphysics/one_moment_microphysics::'//trim(micro_name)
	   endif
	   
         else if (have_option(trim(micro_optionpath)//'/diagnostic')) then
           micro_optionpath=trim(micro_optionpath)//'/diagnostic'
	 endif
	 
	 if (have_option(trim(micro_optionpath)//'/scalar_field::Ndrop')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Ndrop"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/scalar_field::Ndrop'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint
         endif
	 if (have_option(trim(micro_optionpath)//'/scalar_field::Qdrop')) then
	   nsol=nsol+1
           temp_field_name_list(nsol) = "Qdrop"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/scalar_field::Qdrop'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint	      
	 endif
	 if (have_option(trim(micro_optionpath)//'/scalar_field::Nrain')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Nrain"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/scalar_field::Nrain'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint	      
	 endif
	 if (have_option(trim(micro_optionpath)//'/scalar_field::Qrain')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Qrain"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/scalar_field::Qrain'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint   
	 endif
	 if (have_option(trim(micro_optionpath)//'/cold_microphysics/scalar_field::Nice')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Nice"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/cold_microphysics/scalar_field::Nice'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint	      
	 endif
	 if (have_option(trim(micro_optionpath)//'/cold_microphysics/scalar_field::Qice')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Qice"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/cold_microphysics/scalar_field::Qice'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint   
	 endif	  
	 if (have_option(trim(micro_optionpath)//'/cold_microphysics/scalar_field::Ngrau')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Ngrau"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/cold_microphysics/scalar_field::Ngrau'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint	      
	 endif
	 if (have_option(trim(micro_optionpath)//'/cold_microphysics/scalar_field::Qgrau')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Qgrau"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/cold_microphysics/scalar_field::Qgrau'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint   
	 endif	 
	 if (have_option(trim(micro_optionpath)//'/cold_microphysics/scalar_field::Nsnow')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Nsnow"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/cold_microphysics/scalar_field::Nsnow'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint	      
	 endif
	 if (have_option(trim(micro_optionpath)//'/cold_microphysics/scalar_field::Qsnow')) then
           nsol=nsol+1
           temp_field_name_list(nsol) = "Qsnow"
           temp_field_optionpath_list(nsol)=trim(micro_optionpath)//'/cold_microphysics/scalar_field::Qsnow'
           temp_field_state_list(nsol) = p+1
           priority(nsol) = tmpint   
	 endif	 
	 	  	 
         end if
	 
       end do

       ! make sure we have found all ntsol scalar fields:
       assert(nsol==ntsol)

       nsol=0
       do p=maxval(priority),minval(priority),-1
          do f=1,ntsol
             if (priority(f)==p) then
                nsol = nsol + 1
                field_name_list(nsol) = temp_field_name_list(f)
                field_optionpath_list(nsol) = temp_field_optionpath_list(f)
                field_state_list(nsol) = temp_field_state_list(f)
             end if
          end do
       end do

       deallocate( priority, &
            temp_field_name_list, &
            temp_field_state_list, &
            temp_field_optionpath_list )

       initialised = .true.
       
    end if ! End of if(initialised)

    ! Point the list of fields. This has to be done every adapt as the
    ! field structures will be reallocated.
    
    ! Note that we use borrowed references for this so as not to interfere
    ! with adaptivity.
    do f=1,ntsol
       field_list(f) = extract_scalar_field(state(field_state_list(f)),&
            &                               field_name_list(f))
    end do
          
  end subroutine initialise_field_lists_from_options

  subroutine get_ntsol(ntsol)
    integer, intent(out) :: ntsol
    integer :: nphases,nfields,ncars,p,f
    character(len=FIELD_NAME_LEN) :: tmpstring
    logical :: aliased, pressure, surface
    
    ewrite(1,*) 'In get_ntsol'

    ntsol = 0

    nphases = option_count('/material_phase')  
    do p = 0, nphases-1
       nfields = option_count('/material_phase[' &
            //int2str(p)//']/scalar_field')
       do f = 0, nfields-1
          aliased = have_option('/material_phase['// &
               int2str(p)//']/scalar_field['//int2str(f)//']/aliased')
          surface = have_option('/material_phase['// &
               int2str(p)//']/scalar_field['//int2str(f)//']/diagnostic/boundary_conditions')
          call get_option('/material_phase['// &
               int2str(p)//']/scalar_field['//int2str(f)//']/name', tmpstring)
          pressure = (trim(tmpstring)=='Pressure')

          if (.not. aliased .and. .not. pressure .and. .not. surface) then
             ntsol = ntsol + 1
          end if
       end do
       ! prognostic scalar fields for Mellor Yamada:
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::KineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::TurbulentLengthScalexKineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity/prognostic')) then
          ntsol=ntsol + 1
       end if
       ! prognostic scalar fields for k-epsilon turbulence model:
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation/prognostic')) then
          ntsol=ntsol + 1
       end if
       ! prognostic scalar fields for subgrid-scale kinetic energy model:
       if(have_option('/material_phase['//int2str(p)// &
            ']/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/scalar_field::SubgridKineticEnergy')) then
          ntsol=ntsol + 1
       end if
       !Melting
       if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Tb/diagnostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Sb/diagnostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::MeltRate/diagnostic')) then
          ntsol=ntsol + 1
       end if
       !Sediments
       if (have_option('/material_phase['//int2str(p)//']/sediment')) then
          ntsol=ntsol + get_n_sediment_fields()
       end if

       ! scalar fields for cloud microphysics:
       if(have_option('/material_phase['//int2str(p)//']/cloud_microphysics')) then
	  call get_microphysics_ntsol(p,ntsol)
       end if
       
    end do

  end subroutine get_ntsol

  subroutine get_microphysics_ntsol(p,ntsol)
      integer, intent(in) :: p
      integer, intent(inout) :: ntsol
      
      character(len=OPTION_PATH_LEN) :: mom_path  
      character(len=FIELD_NAME_LEN) :: micro_name  
  
      ! In case it's needed later, so we can ensure everything we need exists
      
      mom_path = '/material_phase['//int2str(p)//']/cloud_microphysics'
      if (have_option(trim(mom_path)//'/fortran_microphysics')) then
    	if (have_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics')) then
	  call get_option(trim(mom_path)//'/fortran_microphysics/two_moment_microphysics/name',micro_name)
    	  mom_path=trim(mom_path)//'/fortran_microphysics/two_moment_microphysics::'//trim(micro_name)
    	else if (have_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics')) then
	  call get_option(trim(mom_path)//'/fortran_microphysics/one_moment_microphysics/name',micro_name)
    	  mom_path=trim(mom_path)//'/fortran_microphysics/one_moment_microphysics::'//trim(micro_name)
    	endif
	
      else if (have_option(trim(mom_path)//'/diagnostic')) then	
        mom_path=trim(mom_path)//'/diagnostic'
      endif
     
      if (have_option(trim(mom_path)//'/scalar_field::Qdrop'))then
    	ntsol = ntsol + 1  
      endif  
      if (have_option(trim(mom_path)//'/scalar_field::Ndrop'))then
    	ntsol = ntsol + 1  
      endif  
      if (have_option(trim(mom_path)//'/scalar_field::Qrain'))then
    	ntsol = ntsol + 1  
      endif  
      if (have_option(trim(mom_path)//'/scalar_field::Nrain'))then
    	ntsol = ntsol + 1  
      endif  
      
      if (have_option(trim(mom_path)//'/cold_microphysics')) then
    	mom_path = trim(mom_path)//'/cold_microphysics'
    	if (have_option(trim(mom_path)//'/scalar_field::Qice'))then
    	  ntsol = ntsol + 1  
    	endif  
    	if (have_option(trim(mom_path)//'/scalar_field::Nice'))then
    	  ntsol = ntsol + 1  
    	endif  
    	if (have_option(trim(mom_path)//'/scalar_field::Qsnow'))then
    	  ntsol = ntsol + 1  
    	endif  
    	if (have_option(trim(mom_path)//'/scalar_field::Nsnow'))then
    	  ntsol = ntsol + 1  
    	endif  
    	if (have_option(trim(mom_path)//'/scalar_field::Qgrau'))then
    	  ntsol = ntsol + 1  
    	endif  
    	if (have_option(trim(mom_path)//'/scalar_field::Ngrau'))then
    	  ntsol = ntsol + 1  
    	endif	  
      endif
  
  end subroutine get_microphysics_ntsol

end module field_priority_lists
