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
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
!
#include "fdebug.h"

module Radiation_interface

use spud
use state_module
use fields
use sparse_tools
use fefields
use fldebug
use fetools
use node_ownership
use interpolation_module
use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN,&
      PYTHON_FUNC_LEN, CURRENT_DEBUG_LEVEL
use diagnostic_fields, only: safe_set
use boundary_conditions
use equation_of_state
use futils, only : int2str, present_and_true

implicit none

  !--------------------------------------------------------------------!
  ! radiation input data                                               !
  !   nsw    [1]     number of shortwave intervals                     !
  !   nlw    [1]     number of longwave intervals                      !
  !   mm_co2 [kg/kg] specific concentration of co2                     !
  !   mm_o2  [kg/kg] specific concentation of o2                       !
  !   mm_o3  [kg/kg] specific concentation of o2                       !
  !   slofrc [w/m2]  fraction of solar irradiance per wavelength       !
  !                     interval                                       !
  !   fsol   [w/m2]  solar constant                                    !
  !   wsw            weights of shortwave spectral intervals           !
  !   abswvap[m2/kg] absorption coefficient of water vapour (shortwave)!
  !   abswo3 [m2/kg] absorption coefficient of o3 (shortwave)          !
  !   abswco2[m2/kg] absorption coeffcient of co2 (shortwave)          !
  !   abswo2 [m2/kg] absorption coefficient of o2 (shortwave)          !
  !   tauray         optical depth for rayleigh scattering (shortwave) !
  !   asl, bsl, csl  optical properties of cloudwater (shortwave)      !
  !   dsl, esl, fsl                                                    !
  !   wavenum[1/m]   wavenumber of longwave in 1/micrometre            !
  !   delwav [m]     wave length interval in micrometre                !
  !   coszen         cosine of zenith angle                            !
  !   sf     [/]     solar angle factor (influence on albedo)          !
  !   surface_type   surface type for calculating albedo               !
  !   flxdn_*[W/m2]  short wave (sw) and long wave (lw) ground (gr)    !
  !                  irradiance total (tot) and direct (dir)           !
  !--------------------------------------------------------------------!
  
  ! SW radiation variables
  integer :: nsw
  real    :: mm_co2, mm_o2, n2o, ch4, cfc11, cfc12, ccl4, cfc22, rad_ice, rad_wat, thresh
  real    :: fsol, pss, hscale, dtaumax, dtaumin, ompdiff, expmax, al0, sf
  real, allocatable, dimension(:) :: mm_o3, fdiff, fdir
  real, allocatable, dimension(:) :: wsw, solfrc, tauray, abswo2, abswco2,  &
       abswo3, abswvap, asl, bsl, csl, dsl, esl, fsl, wavnum, delwav, wlw,  &
       bl0, bl1, bl2, cl0, cl1, cl2

  ! LW radiation variables
  integer :: npdiff, nplevel
  integer, allocatable, dimension(:) :: mw
  real, allocatable, dimension(:) :: xkw, xke, aw, bw, pm
  real, allocatable, dimension(:,:) :: fkw, gkw, aib, awb, aiw, aww, aig, awg, cb
  
  integer :: ip, iw
  integer :: nv, surface_type  
  integer, parameter :: nxp=26, no=21, nc=30, nh=31, pi=3.14159
  
  real, dimension(nxp, nc) :: c1, c2, c3
  real, dimension(nxp, no) :: o1, o2, o3
  real, dimension(nxp, nh) :: h11, h12, h13, h21, h22, h23, h81, h82, h83

  ! Ozone profile
  integer :: nmar
  real, dimension(750) :: zmar, omar
  logical, save :: o3flg

  logical :: sw_radiation, lw_radiation 
  logical :: sw_radiation_input, swrad_old, lwrad_old  

! "radiation.in" contains look-up tables for co2, o3 and h2o:
! Parameters c1, c2, c3, o1, o2, o3, h11, h12, h13, h21, h22, h23, h81, h82, h83
#include "radiation.in"  
  
  private

  public :: radiation_init, calculate_radiation

contains

  subroutine radiation_init(state)
    
    type(state_type), intent(inout), dimension(:) :: state
    type(mesh_type), pointer :: radiation_mesh
        
    character(len=OPTION_PATH_LEN) :: rad_path, date_path
    real :: co2, oxy
    integer :: i, nlevel, stat

    namelist /radiation_setup/                         &
         rad_ice, rad_wat, co2, oxy, n2o, ch4, cfc11,  &
         cfc12, cfc22, ccl4, thresh 
	 
    ewrite(1,*) 'In radiation_init'
    
    ! Check if radiation mesh is present 
    radiation_mesh => extract_mesh(state(1),"RadiationMesh",stat=stat)
    if (stat/=0) then
      FLAbort('Trying to initialise Radiation but no mesh specified!')
    endif
    
    rad_path="/radiation_model"
    if (have_option("/radiation_model/simple") .or. &
        have_option("/radiation_model/from_file")) then
      ! Nothing to be done here
    
    else if (have_option("/radiation_model/rrtm")) then
    
    if(have_option('/physical_parameters/location_date'))then
      date_path='/physical_parameters/location_date/julian_day'
    else
      FLAbort("No location and date selected to compute radiation!")
    endif
    sw_radiation=.not.have_option(trim(rad_path)//'/no_short_wave')
    lw_radiation=.not.have_option(trim(rad_path)//'/no_long_wave')
    
    !------------------------------------------------------------------------!
    ! atmospheric trace gas concentrations (part/part)                       !
    ! NOTE: rrtm requires by default volume mixing ratios (Documentation)    !
    ! rrtm code was changed to mm by activating predefined coeffs in rrtm    !
    !------------------------------------------------------------------------!
    rad_ice = 50.
    rad_wat = 10.
    thresh  = 1.0e-9
    co2     = 350.e-6
    oxy     = 0.2095
    n2o     = 0.28e-6
    ch4     = 1.75e-6
    cfc11   = 0.3e-9   ! original rrtm_lw is capable of processing full  
    cfc12   = 0.5e-9   ! matrices of trace gas concentrations: code was  
    cfc22   = 0.2e-9   ! changed assuming that cfcs are well mixed -->
    ccl4    = 0.1e-9   ! only one value 
    mm_o2   = 0.
    mm_co2  = 0.
    o3flg   = .true.
   
    ! read namelist file INPUT_radiation     
    open(56,file='src/radiation/INPUT_radiation',form='formatted',status='old',err=2889)
    read(56,radiation_setup,end=2889)
      
    goto 2888
    
2889 ewrite(3,*) 'NO namelist from INPUT_radiation read, NO DEFAULTS!!!'

2888 close(56)     

    call get_option('geometry/mesh::RadiationMesh/from_mesh/extrude/regions[0]/vertical_levels',nlevel)
    
    if (sw_radiation) then
      mm_o2  = oxy * 1.104801	  ! conversion factor from vm to mm
      mm_co2 = co2* 1.519149	  ! taken from rrtm_lw
      
      call sw_radiation_atham_init(nlevel)
    endif

    if (lw_radiation) then
      call lw_radiation_nasa_init(nlevel)
    end if
    
    endif
   
  end subroutine radiation_init

  subroutine sw_radiation_atham_init(nlevel)                   

    integer, intent(in) :: nlevel                             

    real :: tal0  
    real, dimension(99) :: twsw,tsolfrc, ttauray, tabswo2,tabswco2, tabswo3, tabswvap
    real, dimension(99) :: tasl, tbsl, tcsl, tdsl, tesl, tfsl
    real, dimension(99) :: tbl0, tbl1, tbl2, tcl0, tcl1, tcl2
    real, dimension(240,3,18) :: rad_down_sw 
    integer :: i, a, b, c, level_in

    namelist /radiation_setup_sw/ nsw, dtaumax,dtaumin,ompdiff,expmax,        &
                             pss,hscale,tsolfrc,fsol,                         &
                             twsw, tabswvap, tabswo3, tabswco2,  tabswo2,     & 
                             ttauray, tasl, tbsl, tcsl, tdsl, tesl, tfsl,     &
                             tal0, tbl0, tbl1, tbl2, tcl0, tcl1, tcl2,        &
                             sw_radiation_input, swrad_old

    sw_radiation_input =.false.
    swrad_old =.false.
    
    ! read namelist file INPUT_radiation    
    open(56,file='src/radiation/INPUT_radiation_sw_atham',form='formatted',status='old',err=2999)
    read(56,radiation_setup_sw,end=2999)

    allocate(wsw(nsw), solfrc(nsw), tauray(nsw), abswo2(nsw),           &
             abswco2(nsw), abswo3(nsw), abswvap(nsw))
    allocate(asl(nsw), bsl(nsw), csl(nsw), dsl(nsw), esl(nsw), fsl(nsw))
    allocate(bl0(nsw), bl1(nsw), bl2(nsw), cl0(nsw), cl1(nsw), cl2(nsw))
   
    solfrc  = tsolfrc(1:nsw)
    wsw     = twsw(1:nsw)
    abswvap = tabswvap(1:nsw)
    abswo3  = tabswo3(1:nsw)
    abswco2 = tabswco2(1:nsw)
    abswo2  = tabswo2(1:nsw)
    tauray  = ttauray(1:nsw)
    asl     = tasl(1:nsw)
    bsl     = tbsl(1:nsw)
    csl     = tcsl(1:nsw)
    dsl     = tdsl(1:nsw)
    esl     = tesl(1:nsw)
    fsl     = tfsl(1:nsw)
    bl0     = tbl0(1:nsw)
    bl1     = tbl1(1:nsw)
    bl2     = tbl2(1:nsw)
    cl0     = tcl0(1:nsw)
    cl1     = tcl1(1:nsw)
    cl2     = tcl2(1:nsw)
    al0     = tal0

    goto 2998
    
2999 ewrite(3,*) 'ERROR: NO namelist from INPUT_rad_sw read, NO DEFAULTS!!!'

2998 close(56) 

    !------------------------------------------------------------------!
    ! dertermine radiation input per sw band at model top              !
    !------------------------------------------------------------------!
    allocate (fdiff(nsw), fdir(nsw))
    fdiff = 0.
    fdir  = 0.

    if(sw_radiation_input) then

       level_in = NINT((nlevel-125.)/250.)
       if (level_in<0) then
    	  FLAbort("Can't find thermodynamic variable in state!")
       endif

       open(56,file='src/radiation/SWRIN',form='unformatted')
       READ(56,END=2899) (((rad_down_sw(a,b,c),a=1,240),b=1,3),c=1,18)

       fdiff = rad_down_sw(level_in,2,:)
       fdir  = rad_down_sw(level_in,3,:)

       goto 2898

2899   ewrite(3,*) 'ERROR: NO Radiation Input'

2898   close(56)
    end if

  end subroutine sw_radiation_atham_init
  
  subroutine lw_radiation_nasa_init(nlevel)

    integer, intent(in) :: nlevel

    real :: txkw(9), txke(9), taw(9), tbw(9), tpm(9), tfkw(6,9), tgkw(6,3)
    real :: taib(3,10)
    real, dimension(4,10) :: tawb, taiw, taww, taig, tawg
    real, dimension(6,10) :: tcb
    integer :: tmw(9)

    namelist /radiation_setup_lw/                                             &
                             txkw, txke, tmw, taw, tbw, tpm, tfkw, tgkw,      &
                             taib, tawb, taiw, taww, taig, tawg, tcb,         &
                             npdiff, lwrad_old

    lwrad_old = .false.

    ! read namelist file INPUT_radiation   
    open(56,file='src/radiation/INPUT_radiation_lw_nasa',form='formatted',status='old',err=2999)
    read(56,radiation_setup_lw,end=2999)

    allocate(xkw(9), xke(9), aw(9), bw(9), pm(9), fkw(6,9), gkw(6,3), mw(9))
    allocate(aib(3,10), awb(4,10), aiw(4,10), aww(4,10), aig(4,10), awg(4,10))
    allocate(cb(6,10))

    xkw = txkw
    xke = txke
    mw  = tmw
    aw  = taw
    bw  = tbw
    pm  = tpm
    fkw = tfkw
    gkw = tgkw
    mw  = tmw
    aib = taib
    awb = tawb
    aiw = taiw
    aww = taww
    aig = taig
    awg = tawg
    cb  = tcb
   
    goto 2998
    
2999 ewrite(3,*) 'ERROR: NO namelist from INPUT_rad_nasa read, NO DEFAULTS!!!'

2998 close(56) 

    nplevel = (nlevel + npdiff)
   
  end subroutine lw_radiation_nasa_init

  subroutine define_ozone(nz,z)

    ! local variables 
    integer, intent(in)		    :: nz
    real, dimension(nz), intent(in) :: z
    
    integer              :: k, km
    real		 :: zh1,zh2,dh,dhr,dom
    real, dimension(750) :: zmar, omar

    ! read ozone data 
    if (o3flg) then 
      allocate(mm_o3(1:nz))
      mm_o3=0.
      
      open(9,file='input/INPUT_profile_ozone', form='formatted', status='old',err=999)
      read(9,100,end=999) nmar 
      do k=1,nmar
         read(9,111,end=999) zmar(k), omar(k)
      enddo
      close(9)
      
      o3flg=.false.
    endif
    
    ! interpolate to model layer     
    do km=1,nmar-1
       do k=1,nz
          if (z(k)>=zmar(km) .and. z(k)<zmar(km+1))then
             zh1=zmar(km)
             zh2=zmar(km+1)
             dh =zh2-zh1
             dhr=z(k)-zh1
             dom=(-omar(km)+omar(km+1))/dh
             mm_o3(k)=omar(km)+dhr*dom
          endif
       enddo
    enddo

    return
    
100 format(20x,i4,/,/,/,/)
111 format(f8.1,42x,f10.4)
999 write(*,*) 'No OZONE profile!!'
    o3flg=.false.
     
  end subroutine define_ozone
  
  subroutine calculate_radiation (state,timestep,current_time,dt)
    type(state_type), intent(inout) :: state(:)
    real, intent(in) :: current_time, dt
    integer, intent(in) :: timestep
    integer :: radiation_period
    logical :: do_radiation
    
    call get_option("/radiation_model/radiation_period",radiation_period)
    radiation_period=max(radiation_period,1)
    if (mod(timestep,radiation_period)==0) then
      do_radiation=.true.
    else
      do_radiation=.false.
    endif
    
    call reset_radiation(state(1))
    
    if (do_radiation) then	!<----------------------------------------

      if (have_option("/radiation_model/from_file")) then
        call radiation_from_file(state(1),current_time,dt)
      else if (have_option("/radiation_model/rrtm")) then
        call calculate_radiation_rrtm(state(1),current_time,dt)
      else if (have_option("/radiation_model/simple")) then
        call calculate_radiation_simple(state(1))
      endif
    
    endif
    
    contains
    
    subroutine reset_radiation (state)
      type(state_type), intent(inout) :: state
      type(scalar_field), pointer :: thermal, flux, source
      integer :: stat
      
      call get_thermo_variable (state, thermal)
      source=>extract_scalar_field(state, trim(thermal%name)//'RadiationSource', stat=stat)
      if (stat /= 0) &
          FLAbort("Radiation model selected but no Radiation Source field present")
      call zero(source)
	
      if (have_option('/radiation_model/simple') .or. &
          have_option('/radiation_model/from_file')) then
        flux=>extract_scalar_field(state, 'RadiationFlux', stat=stat)
        if (stat /= 0) &
            FLAbort("Radiation model selected but no Radiation Flux field present")
        call zero(flux)
      else
        flux=>extract_scalar_field(state, 'RadiationLWFlux', stat=stat)
        if (stat /= 0) &
            FLAbort("Radiation model selected but no LW Radiation Flux field present")
        call zero(flux)
	
        flux=>extract_scalar_field(state, 'RadiationSWFlux', stat=stat)
        if (stat /= 0) &
            FLAbort("Radiation model selected but no SW Radiation Flux field present")
        call zero(flux)
      endif
    
    end subroutine reset_radiation
	     
  end subroutine calculate_radiation
  
  subroutine radiation_from_file (state,current_time,dt)

    type(state_type), intent(inout) :: state
    real, intent(in) :: current_time, dt
        
    type(scalar_field), pointer :: pressure, thermal, flux, radiation_source
    type(scalar_field) :: pressure_remap
    
    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN), allocatable :: filename(:)
    integer, allocatable :: time(:)
    
    real :: frac, time_frac, p_node, p0, p1, hr, hr0, hr1, flx, flx0, flx1
    real :: p01, p11, flx01, flx11, hr01, hr11, up_flx0, dn_flx0, up_flx1, dn_flx1
    integer :: i, k, l, nfile, lev, stat, lunit1, lunit2
    logical :: lopened
    
    ewrite(1,*) 'Start read_radiation_file'
    
    call get_thermo_variable (state, thermal)
    
    pressure=>extract_scalar_field(state,"Pressure")
    flux=>extract_scalar_field(state, 'RadiationFlux', stat=stat)
    radiation_source=>extract_scalar_field(state, trim(thermal%name)//'RadiationSource', stat=stat)
    
    call allocate(pressure_remap, flux%mesh, "PressureRemap")
    call safe_set(state,pressure_remap,pressure)

    call get_option("/radiation_model/from_file/file_name", name)
    call time_dependent_filenames(name, nfile, time, filename)
    
    if (nfile == 1) then
      l=0
      stat=1
      do while (stat/=0)
	lunit1=20+l
	inquire(unit=lunit1, opened=lopened, iostat=stat)
	l=l+1
      enddo  
      open (UNIT=lunit1, FILE=trim(name), IOSTAT=stat)

    else
      i=1
      do i = 1, size(time)-1
	if (current_time>=time(i).and.current_time<time(i+1)) then
	
	  l=0
	  stat=1
	  do while (stat/=0)
	    lunit1=20+l
	    lunit2=40+l
	    inquire(unit=lunit1, opened=lopened, iostat=stat)
	    if(stat==0) inquire(unit=lunit2, opened=lopened, iostat=stat)
	    l=l+1
	  enddo  
	  
	  open (UNIT=lunit1, FILE=trim(trim(filename(i))), IOSTAT=stat)	
	  open (UNIT=lunit2, FILE=trim(trim(filename(i+1))), IOSTAT=stat)
	  
	  time_frac=(current_time-real(time(i)))/(real(time(i+1))-real(time(i)))
	  exit
	  
        endif
      enddo
    endif
    
    ! Main loop
    do i = 1, node_count(flux)
      stat = 0
      p_node=node_val(pressure_remap,i)

      ! Read first lines
      do k = 1, 3
        read(lunit1,*)
      enddo
      read(lunit1,*) lev, p0, up_flx0, dn_flx0, flx0, hr0
      p0  = p0*100.     !From mb to Pa
      hr0 = hr0/86400.  !From K/day to K/s
      if (nfile > 1) then
        do k = 1, 3
          read(lunit1,*)
        enddo
        read(lunit2,*) lev, p1, up_flx1, dn_flx1, flx1, hr1
        p1  = p1*100.   	!From mb to Pa
        hr1 = hr1/86400.  !From K/day to K/s
      endif
      
      do while (stat == 0)
    	read(lunit1,*,iostat=stat) lev, p01, up_flx0, dn_flx0, flx01, hr01
    	p01  = p01*100.	    !From mb to Pa
    	hr01 = hr01/86400.    !From K/day to K/s
        
    	frac = (p_node - p0)/(p01 - p0)
	flx  = flx0 + frac*(flx01 - flx0)
	hr   = hr0 + frac*(hr01 - hr0)
	
        if (nfile > 1) then
          read(lunit2,*,iostat=stat) lev, p11, up_flx1, dn_flx1, flx11, hr11
          p11  = p11*100.     !From mb to Pa
          hr11 = hr11/86400.  !From K/day to K/s

          frac = (p_node - p1)/(p11 - p1)
          flx  = flx + time_frac*((flx1+frac*(flx11-flx1)) - flx)
          hr   = hr + time_frac*((hr1+frac*(hr11-hr1)) - hr)
        endif

        ! Set scalar values
	if (frac >= 0. .and. frac <= 1.) then
          call set(radiation_source,i,hr)
          call set(flux,i,flx)
	  exit
	else
	  p0 = p01
	  hr0 = hr01
	  flx0 = flx01
	  if (nfile > 1) then
	    p1 = p11
	    hr1 = hr11
	    flx1 = flx11
	  endif
	  cycle
	endif
      enddo
      
      rewind(lunit1)
      if (nfile > 1) rewind(lunit2)

    enddo
    
    call deallocate(pressure_remap)
     
    ewrite(1,*) 'End read_radiation_file'
  
  contains
  
    subroutine time_dependent_filenames(name, nfiles, time, filename)
    character(len=*), intent(in) :: name
    integer, allocatable, intent(out) :: time(:)
    character(len=OPTION_PATH_LEN), allocatable, intent(out) :: filename(:)
    integer, intent(out) :: nfiles
    
    logical :: ex
    integer :: i, j
    real, dimension(:), allocatable :: time_old
    character(len=OPTION_PATH_LEN), dimension(:), allocatable :: name_old
    character(len=OPTION_PATH_LEN) :: filename2
    
    call get_option("/radiation_model/from_file/time_dependent/number_of_files", nfiles)
    
    if (nfiles > 1) then

    allocate(time(1), filename(1))
    time=0
    j=0
    i=0
    do while (j < nfiles)
      if(i<10)then
        filename2=trim(name)//".0000"//int2str(i)
      else if(i>=10.and.i<100)then
        filename2=trim(name)//".000"//int2str(i)
      else if(i>=100.and.i<1000)then
        filename2=trim(name)//".00"//int2str(i)
      else if(i>=1000.and.i<10000)then
        filename2=trim(name)//".0"//int2str(i)
      else if(i>=10000.and.i<100000)then
        filename2=trim(name)//"."//int2str(i)        
      else
        FLExit("Error in initialise_scalar_field: Number of files is too large.")	  
      endif
      inquire(FILE=trim(filename2),EXIST=ex)
      
      if (ex) then
        if(j/=0) then
	  allocate(time_old(1:j), name_old(1:j))
	  time_old=time
	  name_old=filename
	  deallocate(time, filename)
	  
	  allocate(time(1:j+1), filename(1:j+1))
	  time(1:j)=time_old
	  filename(1:j)=name_old
	  deallocate(time_old, name_old)
	endif
	
	j=j+1
        time(j)=i
	filename(j)=filename2
        if (j==nfiles) exit
      endif
      i=i+1
    enddo
    
    endif

    end subroutine time_dependent_filenames
  
  end subroutine radiation_from_file
    
  subroutine calculate_radiation_simple (state)

    type(state_type), intent(inout) :: state
        
    type(mesh_type), pointer :: radiation_mesh, thermal_mesh
    type(scalar_field), target :: dummyscalar
    type(scalar_field), pointer :: thermal, eospressure, eosdensity, qc_p, qc_v, n_c, q_c, q_r, q_
    
    type(scalar_field) :: pressure_remap, density_remap, temperature_remap, cp_remap, &
    			  cv_remap, qc_remap, qr_remap, nc_remap, radiation_source_remap, flux_remap
    type(scalar_field), pointer :: radiation_source, flux
    type(vector_field), pointer :: radiation_coord, coord
    
    type(scalar_field), dimension(:), allocatable :: new_scalar, old_scalar
    real, allocatable, dimension(:) :: dens, plev, tlev, cplev, cvlev, nclev, wclev, wrlev, z, heatrate, flx
    real :: F0, F1, kappa, lwp, lwp_z, rc, pcom
    integer :: stat, nlevel, ncolumns, radiation_period
    integer :: i, ilev, inode, inode0, thermal_variable
    
    ewrite(1,*) 'Start calculate_radiation_simple'

    call get_thermo_variable (state, thermal, index=thermal_variable)

    thermal_mesh => extract_mesh(state, trim(thermal%mesh%name))
    radiation_mesh=>extract_mesh(state, "RadiationMesh", stat=stat)
    if (stat/=0) then
      FLAbort("Trying to compute radiation, but no RadiationMesh exists")
    endif

    call compressible_eos(state, pressure=eospressure, density=eosdensity, &
    		qc_p=qc_p, qc_v=qc_v)

    call allocate(dummyscalar,thermal_mesh,name="DummyScalar")
    call zero(dummyscalar)

    if (has_scalar_field(state,"Qdrop")) then
      q_c=>extract_scalar_field(state,"Qdrop")
    else
      q_c=>dummyscalar
    end if
    if (has_scalar_field(state,"Ndrop")) then
      n_c=>extract_scalar_field(state,"Ndrop")
    else
      n_c=>dummyscalar
    end if
    if (has_scalar_field(state,"Qrain")) then
      q_r=>extract_scalar_field(state,"Qrain")
    else
      q_r=>dummyscalar
    end if
    
    ! Get options
    call get_option('/radiation_model/simple/F0',F0)
    call get_option('/radiation_model/simple/F1',F1)
    call get_option('/radiation_model/simple/kappa',kappa)

    call get_option('/geometry/mesh::RadiationMesh/from_mesh/extrude/regions[0]/vertical_levels',nlevel)
    nlevel=nlevel+1
    ncolumns=size(radiation_mesh%columns)/real(nlevel)
        
    flux=>extract_scalar_field(state, 'RadiationFlux', stat=stat)
    radiation_source=>extract_scalar_field(state, trim(thermal%name)//'RadiationSource', stat=stat)
    
    ! Allocate scalars on radiation mesh
    call allocate(pressure_remap,radiation_mesh,name="RemappedPressure")
    call allocate(density_remap,radiation_mesh,name="RemappedDensity")
    call allocate(cp_remap,radiation_mesh,name="RemappedCP")
    call allocate(cv_remap,radiation_mesh,name="RemappedCV")
    call allocate(qc_remap,radiation_mesh,name="RemappedQC")
    call allocate(qr_remap,radiation_mesh,name="RemappedQR")
    call allocate(nc_remap,radiation_mesh,name="RemappedNC")
    call allocate(flux_remap,radiation_mesh,name="RemappedFlux")
    call allocate(radiation_source_remap,radiation_mesh,name="RemappedHeatRate")

    call zero(flux_remap)
    call zero(radiation_source_remap)
    
    ! Get coordinates
    coord=>extract_vector_field(state,'Coordinate',stat=stat)
    radiation_coord=>extract_vector_field(state,trim(radiation_mesh%name)//'Coordinate',stat=stat)
    
    ! Time to interpolate scalar fields on the extruded RadiationMesh
    allocate(new_scalar(7))
    allocate(old_scalar(7))
    old_scalar(1) = eospressure
    old_scalar(2) = eosdensity
    old_scalar(3) = qc_p
    old_scalar(4) = qc_v
    old_scalar(5) = q_c
    old_scalar(6) = n_c
    old_scalar(7) = q_r

    call interpolate_on_radiation_mesh (coord,radiation_coord,old_scalar,new_scalar,input=.true.)

    pressure_remap = new_scalar(1)
    density_remap = new_scalar(2)
    cp_remap = new_scalar(3)
    cv_remap = new_scalar(4)
    qc_remap = new_scalar(5)
    nc_remap = new_scalar(6)
    qr_remap = new_scalar(7)
    deallocate(new_scalar)
    deallocate(old_scalar)
    
    ! Allocate columns
    allocate(flx(1:nlevel), z(1:nlevel), heatrate(1:nlevel))    
    allocate(plev(1:nlevel), tlev(1:nlevel), dens(1:nlevel), cplev(1:nlevel), nclev(1:nlevel), wclev(1:nlevel), wrlev(1:nlevel))

    ! Loop over columns
    inode=1
 column_loop: do i=1,ncolumns
    
    z=0.
    inode0=inode
    lwp=0.
    do ilev=1,nlevel
      z(ilev)=node_val(radiation_coord,radiation_coord%dim,inode)
      plev(ilev)=node_val(pressure_remap,inode)
      dens(ilev)=node_val(density_remap,inode)
      cplev(ilev)=node_val(cp_remap,inode)
      cvlev(ilev)=node_val(cv_remap,inode)
      wclev(ilev)=node_val(qc_remap,inode)
      nclev(ilev)=node_val(nc_remap,inode)
      wrlev(ilev)=node_val(qr_remap,inode)
      
      if (ilev > 1) lwp=lwp + (wclev(ilev)+wrlev(ilev)+wclev(ilev-1)+wrlev(ilev-1)) &
                       * 0.5*(z(ilev) - z(ilev-1))
      inode=inode+1
    enddo
    
    lwp_z=0.
    do ilev=1,nlevel
      if (ilev > 1) lwp_z=lwp_z + (wclev(ilev)+wrlev(ilev)+wclev(ilev-1)+wrlev(ilev-1)) &
                       * 0.5*(z(ilev) - z(ilev-1))
      rc=0.
      if (nclev(ilev)>100.) rc=(3./(4.*pi*1000.)*wclev(ilev)/nclev(ilev))**(1./3.)
      
      flx(ilev)=F0*exp(kappa*(lwp-lwp_z)) + F1*exp(-kappa*(lwp-lwp_z))    
      heatrate(ilev)=-kappa/cplev(ilev)*rc*(F0*exp(-kappa*(lwp-lwp_z)) - F1*exp(-kappa*lwp_z))
      if (thermal_variable==1) then        
        heatrate(ilev)=cvlev(ilev)*heatrate(ilev)
      else if (thermal_variable==3) then
        pcom=(plev(ilev)/1.e+05)**((cplev(ilev)-cvlev(ilev))/cplev(ilev))    
        heatrate(ilev)=heatrate(ilev)/pcom
      endif
    enddo
    
    ! 
    inode=inode0
    do ilev=1,nlevel
      call set(flux_remap,inode,flx(ilev))
      call set(radiation_source_remap,inode,heatrate(ilev))
      inode=inode+1
    enddo
    
  enddo column_loop
  
  ! Interpolate fluxes back onto main mesh
  allocate(new_scalar(2))
  new_scalar(1) = flux
  new_scalar(2) = radiation_source

  allocate(old_scalar(2))
  old_scalar(1) = flux_remap
  old_scalar(3) = radiation_source_remap

  call interpolate_on_radiation_mesh (radiation_coord,coord,old_scalar,new_scalar,output=.true.)

  flux = new_scalar(1)
  radiation_source = new_scalar(2)
  deallocate(new_scalar)
  deallocate(old_scalar)

  ! deallocate
  call deallocate(pressure_remap)
  call deallocate(density_remap)
  call deallocate(cp_remap)
  call deallocate(cv_remap)
  call deallocate(qc_remap)
  call deallocate(qr_remap)
  call deallocate(nc_remap)
  call deallocate(flux_remap)
  call deallocate(radiation_source_remap)
  call deallocate(dummyscalar)
  
  deallocate (z, plev, wclev, wrlev, dens, cplev, cvlev, flx, heatrate)

  ewrite(1,*) 'End calculate_radiation_simple'
  
  end subroutine calculate_radiation_simple
    
  subroutine calculate_radiation_rrtm (state,current_time,dt)

    type(state_type), intent(inout) :: state
    real, intent(in) :: current_time,dt
    
    type(mesh_type), pointer :: radiation_mesh, mesh
    real, allocatable, dimension(:) :: z, frcld, heatrate_lw, heatrate_sw, tetflx, sw, lw
    real, allocatable, dimension(:) :: plev, dens, tlev, cplev, cvlev, wclev, wrlev, wvlev, icld
    
    character(len=FIELD_NAME_LEN) :: mesh_name 
    integer :: stat, allocated, nmeshes, nlevel, ncolumns, radiation_period, landseamask
    integer :: i, k, icol, inode, inode0, ilev, iradmesh, klow
    real :: pcom, coszen, albedo, tsfc, emis, dayinv, flxdn_sw_gr_dir, flxdn_sw_gr_tot, flxdn_lw_gr_tot
    
    type(scalar_field), target :: dummyscalar
    type(scalar_field) :: pressure, density, eospressure, eosdensity, temperature, qc_p, qc_v
    type(scalar_field) :: pressure_remap, density_remap, temperature_remap, cp_remap, &
    			  cv_remap, qc_remap, qr_remap, qv_remap
    type(scalar_field) :: radiation_source_remap, sw_flux_remap, lw_flux_remap
    type(scalar_field), dimension(:), allocatable :: new_scalar, old_scalar
    type(scalar_field), pointer :: reference_pressure, q_c, q_r, q_v
    type(scalar_field), pointer :: radiation_source, sw_flux, lw_flux
    type(vector_field), pointer :: radiation_coord, coord
    
    logical :: surface_model, use_dynamic_osa, not_have_pottemp
    
    ewrite(1,*) 'Start calculate_radiation'

    ! Get meshes : Temperature mesh and radiation mesh
    if (have_option('/material_phase[0]/scalar_field::PotentialTemperature/prognostic')) then
    	call get_option('/material_phase[0]/scalar_field::PotentialTemperature/prognostic/mesh/name',mesh_name)
    else if (have_option('/material_phase[0]/scalar_field::Temperature/diagnostic'))then
    	call get_option('/material_phase[0]/scalar_field::Temperature/diagnostic/mesh/name',mesh_name)
    else if (have_option('/material_phase[0]/scalar_field::InternalEnergy/prescribed'))then
    	call get_option('/material_phase[0]/scalar_field::InternalEnergy/prescribed/mesh/name',mesh_name)
    endif  
    
    if (trim(mesh_name)=='PressureMesh') then
      mesh=>extract_mesh(state, "PressureMesh")
    else if (trim(mesh_name)=='CoordinateMesh') then 
      mesh=>extract_mesh(state, "CoordinateMesh")
    else if (trim(mesh_name)=='VelocityMesh') then 
      mesh=>extract_mesh(state, "VelocityMesh")
    endif  
    
    reference_pressure=>extract_scalar_field(state,"Pressure")

    ! Compute and remap thermodynamic variables
    call allocate(eospressure,reference_pressure%mesh,"EOSPressure")
    call allocate(eosdensity,reference_pressure%mesh,"EOSDensity")
    
    call allocate(qc_p,mesh,name="CP")
    call allocate(qc_v,mesh,name="CV")
    call allocate(pressure,mesh,name="Pressure")
    call allocate(density,mesh,name="Density")
    call allocate(temperature,mesh,name="Temperature")
    call allocate(dummyscalar,mesh,name="DummyScalar")

    call compressible_eos(state, pressure=eospressure, density=eosdensity, temperature=temperature, qc_p=qc_p, qc_v=qc_v)
    
    call safe_set(state,pressure,eospressure)
    call safe_set(state,density,eosdensity)
    
    ! Get radiation mesh options
    radiation_mesh=>extract_mesh(state, "RadiationMesh", stat=stat)

    nmeshes=option_count("/geometry/mesh")
    iradmesh=-1
    do i=0,nmeshes-1
      call get_option('/geometry/mesh['//int2str(i)//']/name',mesh_name)
      if (trim(mesh_name)=='RadiationMesh') iradmesh = i
    enddo
    if (iradmesh==-1) then
      FLAbort("Trying to compute radiation, but no RadiationMesh exists!")
    endif
    
    call get_option('/geometry/mesh['//int2str(iradmesh)//']/from_mesh/extrude/regions[0]/vertical_levels',nlevel)
    call get_option('/radiation_model/fixed_surface_albedo',albedo,default=0.4)
    call get_option('/radiation_model/surface_type',landseamask,default=0)
    call get_option('/radiation_model/surface_temperature',tsfc,default=0.)
    nlevel=nlevel+1
    ncolumns=size(radiation_mesh%columns)/real(nlevel)

    ! Extract fluxes	
    call zero(dummyscalar)
    
    if (.not.have_option('/radiation_model/no_short_wave')) then
      sw_flux=>extract_scalar_field(state, 'RadiationSWFlux')
    else
      sw_flux=>dummyscalar
    endif
    if (.not.have_option('/radiation_model/no_long_wave')) then
      lw_flux=>extract_scalar_field(state, 'RadiationLWFlux')
    else
      lw_flux=>dummyscalar
    endif
    
    radiation_source=>extract_scalar_field(state, 'PotentialTemperatureRadiationSource', stat=stat)
    if(stat/=0) then
      not_have_pottemp=.true.
      radiation_source=>extract_scalar_field(state, 'TemperatureRadiationSource', stat=stat)
    endif
    
    ! Get coordinates
    coord=>extract_vector_field(state,'Coordinate',stat=stat)
    radiation_coord=>extract_vector_field(state,trim(radiation_mesh%name)//'Coordinate',stat=stat)
    
    ! Allocate columns
    allocate(tetflx(1:nlevel), sw(1:nlevel), lw(1:nlevel))
    allocate(z(1:nlevel), frcld(1:nlevel), heatrate_sw(1:nlevel), heatrate_lw(1:nlevel))    
    allocate(plev(1:nlevel), tlev(1:nlevel), dens(1:nlevel), cplev(1:nlevel), cvlev(1:nlevel), wclev(1:nlevel), wrlev(1:nlevel), wvlev(1:nlevel), icld(1:nlevel))

    if (has_scalar_field(state,"TotalWaterQ")) then
      q_v=>extract_scalar_field(state,"TotalWaterQ")
    else
      q_v=>dummyscalar
    end if
    if (has_scalar_field(state,"Qdrop")) then
      q_c=>extract_scalar_field(state,"Qdrop")
    else
      q_c=>dummyscalar
    end if
    if (has_scalar_field(state,"Qrain")) then
      q_r=>extract_scalar_field(state,"Qrain")
    else
      q_r=>dummyscalar
    end if
    call addto(q_v,q_c,scale=-1.)
    call addto(q_v,q_r,scale=-1.)
    
    ! Allocate scalars on radiation mesh
    call allocate(pressure_remap,radiation_mesh,name="RemappedPressure")
    call allocate(density_remap,radiation_mesh,name="RemappedDensity")
    call allocate(temperature_remap,radiation_mesh,name="RemappedTemperature")
    call allocate(cp_remap,radiation_mesh,name="RemappedCP")
    call allocate(cv_remap,radiation_mesh,name="RemappedCV")
    call allocate(qc_remap,radiation_mesh,name="RemappedQC")
    call allocate(qr_remap,radiation_mesh,name="RemappedQR")
    call allocate(qv_remap,radiation_mesh,name="RemappedQV")

    call allocate(sw_flux_remap,radiation_mesh,name="RemappedSWFlux")
    call allocate(lw_flux_remap,radiation_mesh,name="RemappedLWFlux")
    call allocate(radiation_source_remap,radiation_mesh,name="RemappedHeatRate")

    call zero(sw_flux_remap)
    call zero(lw_flux_remap)
    call zero(radiation_source_remap)
    
    ! Time to interpolate scalar fields on the extruded RadiationMesh
    allocate(new_scalar(8))
    new_scalar(1) = pressure_remap
    new_scalar(2) = density_remap
    new_scalar(3) = temperature_remap
    new_scalar(4) = cp_remap
    new_scalar(5) = cv_remap
    new_scalar(6) = qc_remap
    new_scalar(7) = qr_remap
    new_scalar(8) = qv_remap

    allocate(old_scalar(8))
    old_scalar(1) = pressure
    old_scalar(2) = density
    old_scalar(3) = temperature
    old_scalar(4) = qc_p
    old_scalar(5) = qc_v
    old_scalar(6) = q_c
    old_scalar(7) = q_r
    old_scalar(8) = q_v

    call interpolate_on_radiation_mesh (coord,radiation_coord,old_scalar,new_scalar,input=.true.)

    pressure_remap = new_scalar(1)
    density_remap = new_scalar(2)
    temperature_remap = new_scalar(3)
    cp_remap = new_scalar(4)
    cv_remap = new_scalar(5)
    qc_remap = new_scalar(6)
    qr_remap = new_scalar(7)
    qv_remap = new_scalar(8)
    deallocate(new_scalar)
    deallocate(old_scalar)
      
    call solar_angle(current_time,coszen)

    icld = 0.
    inode=1
    dayinv  = 1./86400.
    surface_model=.false.   ! Presence of a surface model (temporary)
    use_dynamic_osa=.false. ! If surface model, albedo computed internally (temporary)

    ! Loop over columns
 column_loop: do i=1,ncolumns
    
    inode0      =inode
    frcld       =1.
    heatrate_sw =0.
    heatrate_lw =0.
    tetflx      =0.
    sw          =0.
    lw          =0.
    
    z=0.
    do ilev=1,nlevel
      z(ilev)=node_val(radiation_coord,radiation_coord%dim,inode)
      plev(ilev)=node_val(pressure_remap,inode)
      tlev(ilev)=node_val(temperature_remap,inode)
      dens(ilev)=node_val(density_remap,inode)
      cplev(ilev)=node_val(cp_remap,inode)
      cvlev(ilev)=node_val(cv_remap,inode)
      wclev(ilev)=node_val(qc_remap,inode)
      wrlev(ilev)=node_val(qr_remap,inode)
      wvlev(ilev)=node_val(qv_remap,inode)
      inode=inode+1
    enddo
    
    ! Interpolate ozone profile
    call define_ozone (nlevel,z)
    
    if (surface_model) then
!       emis = es_atham     ! ground emissivity from surfacemodel
!       tsfc = ts_atham	   ! surface temp from surfacemodel  
    else
       emis = 1. 	                     ! black body emissivity
       if (tsfc == 0) tsfc=tlev(1)        ! use first model layer as tsfc
    endif

    where (icld .le. thresh) 
       icld = 0.
       frcld= 0.
    end where
    where (wclev .le. thresh) 
       wclev= 0. 
       frcld= 0. 
    end where
  
    ! Calculate short-wave radiation
    if (.not.have_option ('/radiation_model/no_short_wave')) then
!       print*, 'before SW radiation'
!       call shortwave_radiation_atham (dt, nlevel, z, dens, 		 & 
!       				  cplev, cvlev, plev,    	 &
!                                  heatrate_sw, coszen,        	  	 &
!                                  wclev, wvlev, icld,   	    	 &
!                                  mm_o2, mm_co2, mm_o3,                          &
!                                  landseamask, surface_model, use_dynamic_osa,   &
!                                  not_have_pottemp, albedo, sw,      	          &
!				  flxdn_sw_gr_dir, flxdn_sw_gr_tot )
!
      do k=1,nlevel
        if (not_have_pottemp) then        
          tetflx(k)= tetflx(k) + heatrate_sw(k)
        else
          pcom=(plev(k)/1.e+05)**((cplev(k)-cvlev(k))/cplev(k))	   
          tetflx(k)= tetflx(k) + heatrate_sw(k)/pcom
        endif
      enddo
    endif
    
    ! Calculate long-wave radiation
    if (.not.have_option ('/radiation_model/no_long_wave')) then
!       print*, 'before LW radiation'
!       call longwave_radiation_nasa (nlevel, tlev, plev, wvlev,   		       &
!                                mm_co2, n2o, ch4,mm_o2, cfc11, cfc12,cfc22,ccl4,      &
!                                mm_o3, frcld, icld, wclev, wrlev,   		       &
!				rad_ice, rad_wat, emis, tsfc, 		   	      &
!				surface_model, heatrate_lw,   		   	      &
!                                lw, flxdn_lw_gr_tot )
!
      do  k=1,nlevel
        if (not_have_pottemp) then        
          tetflx(k)= tetflx(k) - heatrate_lw(k)*dayinv
        else
          pcom=(plev(k)/1.e+05)**((cplev(k)-cvlev(k))/cplev(k))	   
          tetflx(k)= tetflx(k) - heatrate_lw(k)*dayinv/pcom
        endif
      enddo
    end if
    
    ! 
    inode=inode0
    do ilev=1,nlevel
      call set(sw_flux_remap,inode,sw(ilev))
      call set(lw_flux_remap,inode,lw(ilev))
      call set(radiation_source_remap,inode,tetflx(ilev))
      inode=inode+1
    enddo
    
  enddo column_loop
  
  ! Interpolate fluxes back onto main mesh
  allocate(new_scalar(3))
  new_scalar(1) = sw_flux
  new_scalar(2) = lw_flux
  new_scalar(3) = radiation_source

  allocate(old_scalar(3))
  old_scalar(1) = sw_flux_remap
  old_scalar(2) = lw_flux_remap
  old_scalar(3) = radiation_source_remap

  call interpolate_on_radiation_mesh (radiation_coord,coord,old_scalar,new_scalar,output=.true.)

  sw_flux = new_scalar(1)
  lw_flux = new_scalar(2)
  radiation_source = new_scalar(3)
  deallocate(new_scalar)
  deallocate(old_scalar)

  ! deallocate
  call deallocate(pressure_remap)
  call deallocate(density_remap)
  call deallocate(temperature_remap)
  call deallocate(cp_remap)
  call deallocate(cv_remap)
  call deallocate(qv_remap)
  call deallocate(qr_remap)
  call deallocate(qc_remap)
  
  call deallocate(sw_flux_remap)
  call deallocate(lw_flux_remap)
  call deallocate(radiation_source_remap)

  call deallocate(eospressure)
  call deallocate(eosdensity)
  call deallocate(pressure)
  call deallocate(density)
  call deallocate(temperature)
  call deallocate(qc_p)
  call deallocate(qc_v)    
  call deallocate(dummyscalar)
  
  deallocate (frcld, wvlev, wclev, wrlev, icld, dens, cplev, cvlev, tetflx, sw, lw)
  deallocate (z, plev, tlev, heatrate_sw, heatrate_lw)

  ewrite(1,*) 'End calculate_radiation'
  
  end subroutine calculate_radiation_rrtm

  subroutine solar_angle(current_time,coszen)

    real, intent(in) :: current_time
    real, intent(inout) :: coszen
    
    integer :: current_jday, hour_start
    real :: deglat
    real :: dechr, decday, sin_lat, cos_lat, decmonth
    real :: sun_azim, sun_lat
    real :: pi=3.14159
    
    call get_option('/physical_parameters/location_date/latitude',deglat)
    call get_option('/physical_parameters/location_date/julian_day',current_jday)
    call get_option('/physical_parameters/location_date/hour_start',hour_start)
    
    dechr    = (current_time+real(hour_start))/24./3600.
    decmonth = 1.+((real(current_jday)+dechr)/365.*12.)

    sin_lat=sin(pi/180.*deglat)
    cos_lat=cos(pi/180.*deglat)  

    sun_azim=2.*pi*dechr  !sun_azim=2.*pi*(decday + deglon/360.-time_zone/24.) 
    sun_lat=(pi*23./180.)*COS((decmonth-6.75)*pi/6.) ! radians
    
    coszen = max(-cos(sun_lat)*cos(sun_azim)*cos_lat + sin(sun_lat)*sin_lat, 1.e-9)

  end subroutine solar_angle

  subroutine shortwave_radiation_atham                                      &
                                 (dt, nlevel, zv, density, cptot, cvtot,    &
                                  pnew, heatrate,coszen, 	   	    &
                                  watcnew, wetnew, icenew,      	         &
                                  oxy2, carbo2, ozone,   	               	         &
                                  landseamask, surface_model, use_dynamic_osa,   &
                                  not_have_theta, albedo, flux, 		 &
				  flxdn_sw_gr_dir, flxdn_sw_gr_tot)
 
    logical, intent(in) :: surface_model, use_dynamic_osa, not_have_theta
    integer, intent(in) :: nlevel, landseamask

    real, intent(in)   :: dt, oxy2, carbo2, coszen
    real, DIMENSION(:), INTENT(in) :: density, cptot, cvtot, pnew, wetnew, watcnew, icenew, ozone, zv
    
    real, dimension(:), intent(inout) :: heatrate, flux 
    real, intent(inout)   :: flxdn_sw_gr_tot, flxdn_sw_gr_dir, albedo

    ! Local
    real, dimension(nlevel) :: flxup, flxdn, dtaup, omp, gp
    REAL :: dzv, tempflx, cp, pcom, albedo_new

    integer :: k, l, ku, nv

    !-------------------------------------------------------------------------!
    ! Set surface_type                                                        !
    ! Currently set to landseamask                                            !
    !   0 = ocean                                                             !
    !   1 = Tall/medium grassland, evergreen shrubland as defined             !
    !       in Briegleb 1992                                                  !
    ! Possible to add further definitions!                                    !
    !-------------------------------------------------------------------------!
    surface_type = landseamask

    nv = 2*nlevel-1
    do k=1,nlevel
      flxup(k)=0.
      flxdn(k)=0.
    enddo
    flxdn_sw_gr_dir=0.
    
    ! perform short wave calculations for nsw wave length intervals
    do l=1,nsw
       call opticsw(l, dtaup, omp, gp, density, pnew, 		 	&
    		    albedo_new, albedo, surface_model, oxy2, carbo2,   	&
    		    ozone, watcnew, wetnew, icenew,			&
    		    nlevel, zv, coszen, use_dynamic_osa )
		    
       call fluxsw(l, dtaup, omp, gp, flxup, flxdn, albedo_new, density, 	&
    		   watcnew, wetnew, nlevel, flxdn_sw_gr_dir, coszen )
    enddo
    
    do k=2,nlevel
       ku=k-1
       dzv=zv(k)-zv(ku)
       flux(k)=flxdn(k)-flxup(k)    
       if (k==2) flxdn_sw_gr_tot = flxdn(k)

       heatrate(k)=((flxdn(k)-flxup(k)) - (flxdn(ku)-flxup(ku))) / (density(k)*cptot(k)*dzv)  
    enddo
      
  end subroutine shortwave_radiation_atham

  subroutine opticsw(l,dtaup,omp,gp, density, pnew, albedo,              &
                     al_ATHAM, surface_model,                            &
                     oxy2, carbo2, ozone, watcnew, wetnew, icenew,          &
                     nlevel, zv, coszen, use_dynamic_osa)
		     
    !=======================================================================!
    !    calculate optical properties:                                      !
    !    for wave length l in the short wave region                         !
    !                                                                       !
    !    dtaup -- (scaled) optical depth                                    !
    !    omp   -- (scaled) single scattering albedo                         !
    !    gp    -- (scaled) asymmetry factor                                 !
    !=======================================================================!
    
    real, dimension(nlevel), intent(inout) :: dtaup, omp, gp
    real, intent(inout)     	       :: albedo
    real, intent(in)                   :: oxy2, carbo2, coszen, al_ATHAM
    real, dimension(:), intent(in)     :: ozone, zv
    real, dimension(:), intent(in)     :: density, pnew, watcnew, wetnew, icenew

    LOGICAL, INTENT(in)                :: surface_model, use_dynamic_osa
    
    integer, intent(in) :: l, nlevel
    
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !                                                                  !
    !    aerosols:                                                     !
    !    scaaer - scattering coefficient                               !
    !             ( =  specific scattering coefficient   in m^2/kg     !
    !                * specific concentration of aerosol in kg/kg)     !
    !             dependent on wave lentgth interval l                 !
    !    extaer - extinction coefficient in 1/m                        !
    !             (see comment to scaaer)                              !
    !    gaer   - asymmetry factor, dimensionless                      !
    !                                                                  !
    !    optical properties for cloud droplets:                        !
    !    droplet radius in meter                                       !
    !    cloud cover    [0.,1.]                                        !
    !------------------------------------------------------------------!
    
    real :: scaaer=0., extaer=0., gaer  =0.
    real :: radwatc, clcover, radiusw, extwat, omwat, gwat, scawat
    real :: radice, radiusi, extice, omice, gice, scaice
    real :: pssinv, rayleigh, pfac, scamie, g, f, absgas, extinc
    real :: dtau, om, omdiff
    REAL :: dz, co, oxy, fac

    integer ::k
    
    radwatc=5.32e-6
    radice = 50e-6
    clcover=1.
    radiusw=radwatc*1.e6
    radiusi=radice*1.e6
    
    extwat=1000.*(asl(l)+bsl(l)/radiusw)*clcover*sqrt(clcover)
    omwat =1.-csl(l)-dsl(l)*radiusw
    gwat  =   esl(l)+fsl(l)*radiusw
    scawat=extwat*omwat

    extice= 1000.*(al0/radiusi)*clcover*sqrt(clcover)
    omice = 1. - bl0(l) - bl1(l)*radiusi - bl2(l)*radiusi*radiusi
    gice  = cl0(l) + cl1(l) * radiusi + cl2(l) *radiusi*radiusi
    scaice= extice*omice

    !-----------------------------------------------------------------!
    ! start loop over xz-slice                                        !
    !-----------------------------------------------------------------!
    
    !marker 
   
    IF (swrad_old) THEN
       co = 4.558e-4
       oxy=.2314
       fac = 50.
    else
       co = carbo2
       oxy= oxy2
       fac= 1.
    END IF

    pssinv=1./pss
    
    do k=2,nlevel
       dz=zv(k)-zv(k-1)

       pfac=pnew(k)*pssinv
       scamie=scawat*watcnew(k) + scaaer + scaice*icenew(k)
       rayleigh=tauray(l)*pfac/(hscale*density(k))
       g=(gwat*scawat*watcnew(k) + gice*scaice*icenew(k) + gaer*scaaer) &
            /max(1.e-20,scamie+rayleigh) 
       f=g*g
       
       absgas= abswvap(l) *					&
              (wetnew(k)*pfac + 0.0008*sqrt(coszen*wetnew(k)))+ &
              abswco2(l)*sqrt(coszen*co)		      + &
              abswo2 (l)*sqrt(coszen*oxy)		      + &
              abswo3 (l)* fac *ozone(k)     

       extinc=extwat*watcnew(k) + extaer + rayleigh + absgas + extice*icenew(k)
       dtau=extinc*density(k)*dz

       dtau=max(dtaumin,dtau)
       om=(scamie+rayleigh)/max(1.e-20,extinc)
       dtaup(k)=(1.-om*f)*dtau
       omp(k)  =(1.-f)*om/(1.-om*f)

       omp(k)=min(omp(k),1.-ompdiff)
       gp(k) =(g-f)/(1.-f)
    enddo

    !----------------------------------------------------------------!
    ! Calculate albedo                                               !
    ! sf is solar angle factor (now a public variable)               !
    !----------------------------------------------------------------!
    
    sf = 1.4 / (1.0 + 0.4 * coszen)
    if (surface_model) then
       if (surface_type .eq. 0) then
    	  if (use_dynamic_osa) then
    	     albedo = min(al_ATHAM,1.)     ! sf included in dynamic_osa
    	  else
    	     albedo = min(al_ATHAM * sf,1.)
    	  endif

       elseif (surface_type .eq. 1) then
    	  albedo = min(al_ATHAM * sf,1.)
       endif
    else
       if (surface_type .eq. 0) then
    	  albedo = 0.06 * sf
       else if (surface_type .eq. 1) then
    	  albedo = 0.19 * sf
       else 
    	  albedo = 0.0
    	  ewrite(3,*) 'surface type undefined'	    
       endif
    endif

  end subroutine opticsw

!=============================================================================¬

  subroutine fluxsw(l,dtaup,omp,gp,flxup,flxdn, albedo, density,	 &
                    watcnew, wetnew, nlevel, flxdn_sw_gr_dir, coszen)

    !-----------------------------------------------------------------------!
    ! prepare matrix equation for delta-eddington approximation,            !
    ! compute solution and fluxes                                           !
    !                                                                       !
    ! input:   coszen -- cosine of zenith angle                             ! 
    !          albedo -- ground albedo                                      !
    !          dtaup  -- opt depth for each layer                           !
    !          omp    -- sing scat albedo for each layer                    !
    !          gp     -- asym factor for each layer                         !
    !                                                                       !
    !          fsun   -- solar beam flux                                    !
    !          fdinc  -- incident diffuse down-flux at level nlevel (top)       !
    !          fuinc  -- incident diffuse up-flux at level 1 (bottom)       !
    !                                                                       !
    ! output:  flxup, flxdn -- up/down fluxes in w/m^2                      ! 
    !-----------------------------------------------------------------------!
    
    real, intent(inout) :: flxdn_sw_gr_dir
    real, dimension(nlevel), intent(in) :: omp, gp, density, watcnew, wetnew 
    real, intent(in) :: albedo
    real, dimension(nlevel), intent(inout) :: dtaup, flxup, flxdn
    real, intent(in)               :: coszen
    integer, intent(in) :: l, nlevel

    ! local
    integer ::  k, km1, k2, k2p1, k2m1, klow
    real :: fdinc, fuinc, fsun, dir
    real :: h1, h2, h3, xk, xk2, expon, smu 
    real :: pp, pm
    real, dimension(nlevel) :: p, ex, extau
    real, dimension(nv,6) :: pen
    real, dimension(nv) :: solut
    real, dimension(nlevel) :: alph, beta
    character(len=50) :: message
 
    if (sw_radiation_input) then 
       fdinc = fdiff(l)*coszen
       fsun  = fdir(l)
       fuinc = 0.
    else
       fdinc=0.
       fuinc=0.
       fsun =wsw(l)*solfrc(l)*fsol
    endif
    
    if (nv /= 2*nlevel-1) then
       ewrite(3,*) 'error: nv=', nv, '> 2*nlevel-1=', 2*nlevel-1
       return
    endif
    
    !-----------------------------------------------------------------------!
    ! precalculate coefficients                                             !
    !-----------------------------------------------------------------------!

    do k=2,nlevel
       h1 =1.-omp(k)
       h2 =1.-omp(k)*gp(k)
       xk2=3.*h1*h2
       h3 =1./(coszen*coszen)-xk2

    !-----------------------------------------------------------------------!
    ! cut in h3                                                             !
    !-----------------------------------------------------------------------!
    	 h3 =sign(1.,h3)*max(1.e-20,abs(h3))
    	 xk =sqrt(xk2)
    	 p(k)=2.*xk/(3.*h2)
    	 alph(k)=0.75*fsun*omp(k)*(1.+gp(k)*h1)/h3
    	 beta(k)=0.5*fsun*omp(k)*(1./coszen+3.*coszen*gp(k)*h1)/h3

    !----------------------------------------------------------------------!
    ! cut in dtaup							   !
    !----------------------------------------------------------------------!
       dtaup(k)=min(xk*dtaup(k),dtaumax)/xk
       expon=xk*dtaup(k)
       ex(k)=exp(expon)
    enddo

    extau(nlevel)=dtaup(nlevel)
    do k=nlevel-1,2,-1
       extau(k)=extau(k+1)+dtaup(k)
    enddo
    
    smu=coszen
    do k=2,nlevel
       expon=min(expmax,extau(k)/smu)
       extau(k)=1./exp(expon)
    enddo

    !--------------------------------------------------------------------!
    ! calculate matrix: upper boundary  				 !
    !--------------------------------------------------------------------!

    pen(nv,1) = 0.
    pen(nv,2) = 0.
    pen(nv,3) = 1.
    pen(nv,4) = 0.
    pen(nv,5) = 0.
    pen(nv,6) = 0.

    pen(nv-1,5)=0.
    pen(nv-1,4)=0.
    pen(nv-1,3)=(1.+p(nlevel))*ex(nlevel)
    pen(nv-1,2)=(1.-p(nlevel))/ex(nlevel)
    pen(nv-1,1)=0.
    pen(nv-1,6)=alph(nlevel)+beta(nlevel)+fdinc

    !--------------------------------------------------------------------!
    ! atmosphere							 !
    !--------------------------------------------------------------------!

    do k=2,nlevel
       km1=k-1
       k2=2*k-4
       k2p1=k2+1
       pen(k2  ,5)= p(k)
       pen(k2  ,4)=-p(k)
       pen(k2  ,3)=-p(km1)*ex(km1)
       pen(k2  ,2)= p(km1)/ex(km1)
       pen(k2  ,1)= 0.
       pen(k2  ,6)= (beta(k)-beta(km1))*extau(k)
    		   
       pen(k2p1,5)= 0.
       pen(k2p1,4)=pen(k2,5)-1.
       pen(k2p1,3)=pen(k2,4)-1.
       pen(k2p1,2)=pen(k2,3)+ex(km1)
       pen(k2p1,1)=pen(k2,2)+1./ex(km1)
       pen(k2p1,6)=pen(k2,6)-(alph(k)-alph(km1))*extau(k)
    enddo

    !---------------------------------------------------------------------!
    ! lower boundary							  !
    !---------------------------------------------------------------------!

    k=klow
    k2=2*k-4
    k2p1=k2+1
    pp=1.+p(k)
    pm=1.-p(k)
    pen(k2p1,5)=0.
    pen(k2p1,4)=pm-pp*albedo
    pen(k2p1,3)=pp-pm*albedo
    pen(k2p1,2)=0.
    pen(k2p1,1)=0.
    pen(k2p1,6)=fuinc+(alph(k)-beta(k)+albedo*(fsun*coszen-alph(k)-beta(k)))*extau(k)

    !----------------------------------------------------------------------!
    ! ground								   !
    !----------------------------------------------------------------------!

    do k=3,klow
       k2=2*k-4
       k2m1=k2-1
       pen(k2  ,5)=0.
       pen(k2  ,4)=0.
       pen(k2  ,3)=1.
       pen(k2  ,2)=0.
       pen(k2  ,1)=0.
       pen(k2  ,6)=0.
    		   
       pen(k2m1,5)=0.
       pen(k2m1,4)=0.
       pen(k2m1,3)=1.
       pen(k2m1,2)=0.
       pen(k2m1,1)=0.
       pen(k2m1,6)=0.
    enddo
    
    !-----------------------------------------------------------------------!
    ! solve matrix system						    !
    !-----------------------------------------------------------------------!

    call pendia(pen,solut)

    !-----------------------------------------------------------------------!
    ! calculate fluxes  						    !
    !-----------------------------------------------------------------------!

    dir=fsun*coszen
    flxdn(nlevel)=flxdn(nlevel)+fdinc+dir
    flxup(nlevel)=flxup(nlevel)			     &
    		 +(1.+p(nlevel))/ex(nlevel)*solut(nv-2)      &
    		 +(1.-p(nlevel))*ex(nlevel)*solut(nv-1)      &
    		 -alph(nlevel)+beta(nlevel)

    do k=2,nlevel
       k2=2*k-3
       k2p1=k2+1
       dir=fsun*coszen*extau(k)
       flxdn(k-1)= flxdn(k-1)			     &
    		     +(1.-p(k))*solut(k2)	     &
    		     +(1.+p(k))*solut(k2p1)	     &
    		     -(alph(k)+beta(k))*extau(k)     &
    		     +dir
       flxup(k-1)= flxup(k-1)			     &
    		     +(1.+p(k))*solut(k2)	     &
    		     +(1.-p(k))*solut(k2p1)	     &
    		     -(alph(k)-beta(k))*extau(k)     
    enddo
    	     
    flxdn_sw_gr_dir = flxdn_sw_gr_dir + fsun*coszen*extau(1)

  end subroutine fluxsw

  subroutine pendia(pen,solut)

    !=======================================================================!
    !                                                                       !
    !    purpose : parallel solution of nx pentadiagonal linear equations   !
    !                                                                       !
    !    method  : gaussian elimination                                     !
    !    author  : michael herzog                                           !
    !    date    : 21 june 1994                                             !
    !                                                                       ! 
    !=======================================================================!
    
    integer :: l
    real, dimension(nv) :: solut(nv)
    real, dimension(nv,6) :: pen
    real :: fac1, fac2
    
    ! normalization of the main diagonal      
    do l=1,nv
       fac1=1./pen(l,3)
       pen(l,1)=pen(l,1)*fac1
       pen(l,2)=pen(l,2)*fac1
       pen(l,3)=pen(l,3)*fac1
       pen(l,4)=pen(l,4)*fac1
       pen(l,5)=pen(l,5)*fac1
       pen(l,6)=pen(l,6)*fac1
    enddo
    
    ! triangularization 
    do l=2,nv-1
       fac1=pen(l  ,2)/pen(l-1,3)
       fac2=pen(l+1,1)/pen(l-1,3)
       
       pen(l  ,6) = pen(l  ,6) - pen(l-1,6)*fac1
       pen(l  ,3) = pen(l  ,3) - pen(l-1,4)*fac1
       pen(l  ,4) = pen(l  ,4) - pen(l-1,5)*fac1
       
       pen(l+1,6) = pen(l+1,6) - pen(l-1,6)*fac2
       pen(l+1,2) = pen(l+1,2) - pen(l-1,4)*fac2
       pen(l+1,3) = pen(l+1,3) - pen(l-1,5)*fac2
    enddo
    
    ! calculation of the unknowns
    solut(nv  ) =  pen(nv  ,6)/pen(nv  ,3)
    solut(nv-1) = (pen(nv-1,6)-pen(nv-1,4)*solut(nv)) /pen(nv-1,3)
    
    do l=nv-2,1,-1
       solut(l)=(pen(l,6) -pen(l,4)*solut(l+1) -pen(l,5)*solut(l+2))/pen(l,3)
    enddo

  end subroutine pendia

  subroutine longwave_radiation_nasa                                       &
                               (nz, tempnew, pnew, wetnew,	           &
                                co2, n2o, ch4, o2, cfc11, cfc12,           &
                                cfc22, ccl4, ozone,                        &
                                frcld, icld, wcld, wrain, radice, radwat,  &
                                emis, tsfc, surface_model, heatlwrate,     &
                                flxlwtmp, flxdn_lw_gr_tot)
				
    !=======================================================================!
    !                                                                       !
    ! solution of the radiation transfer equation                           !
    ! in delta-eddington approxomation for longwave                         !
    !                                                                       !
    ! calculation of heatingrates                                           !
    !                                                                       !
    !   THE EQUATION NUMBERS noted in this code follows the latest          !
    !    version (July 2002) of the NASA Tech. Memo. (2001), A Thermal      !
    !    Infrared Radiation Parameterization for Atmospheric Studies, where !
    !    a detailed description of the radiation routine can be found.      !
    !                                                                       !
    !                                                                       !
    !   The IR spectrum is divided into nine bands:                         !
    !                                                                       !
    !   band     wavenumber (/cm)   absorber                                !
    !                                                                       !
    !    1           0 - 340           h2o                                  !
    !    2         340 - 540           h2o                                  !
    !    3         540 - 800       h2o,cont,co2                             !
    !    4         800 - 980       h2o,cont                                 !
    !                              co2,f11,f12,f22                          !
    !    5         980 - 1100      h2o,cont,o3                              !
    !                              co2,f11                                  !
    !    6        1100 - 1215      h2o,cont                                 !
    !                              n2o,ch4,f12,f22                          !
    !    7        1215 - 1380      h2o,cont                                 !
    !                              n2o,ch4                                  !
    !    8        1380 - 1900          h2o                                  !
    !    9        1900 - 3000          h2o                                  !
    !                                                                       !
    !   In addition, a narrow band in the 17 micrometer region (Band 10) is !
    !   added to compute flux reduction due to n2o                          !
    !                                                                       !
    !    10        540 - 620       h2o,cont,co2,n2o                         !
    !                                                                       !
    !   Band 3 (540-800/cm) is further divided into 3 sub-bands :           !
    !                                                                       !
    !   subband   wavenumber (/cm)                                          !
    !                                                                       !
    !    3a        540 - 620                                                !
    !    3b        620 - 720                                                !
    !    3c        720 - 800                                                !
    !                                                                       !
    !   August 2009: several changes made to code in order to account for   !
    !   topography by Tobias Gerken (see (*))                               !
    !                                                                       !
    !=======================================================================!

    logical, intent(in) :: surface_model
    integer, intent(in) :: nz
    real, intent(in) :: co2, n2o, ch4, o2, cfc11, cfc12, cfc22, ccl4, &
                     	radice, radwat, emis, tsfc
    real, dimension(nz), intent(in) :: ozone, frcld, icld, wcld, wrain,  &
                                       tempnew, pnew, wetnew

    real, intent(inout) :: flxdn_lw_gr_tot
    real, dimension(nz), intent(out):: heatlwrate, flxlwtmp

    ! temporary arrays   
    real, dimension(10) :: eg
    real, dimension(nplevel) :: pl, flx
    real, dimension(nplevel-1) :: ta, wa, oa, fcld
    real, dimension(nplevel-1, 3) :: cwc, reff
    real :: xpi, tg, tb
    real :: avg_pl, deltap, radwat2, radice2
    real :: co2lw, ch4lw, n2olw, cfc11lw, cfc12lw, cfc22lw

    integer :: ibin, k, klow
    
    klow = 1
    xpi =4./3.*pi

    ! set ground emissivity, eg, read from process_data if surface  
    ! module is  coupled and initialized   
    do ibin = 1 , 10
      eg(ibin) = emis
    enddo

    ! assign ground temperature, tg 
    ! perturbation : dtempg is a randomly selected number with max
    ! amplitude   
    tg = tsfc

    ! assign air temperature at ground 
    if (surface_model) then
       tb = tempnew(1)
    else
       tb = tg
    endif
  
    ! assign variables in pressure levels
    IF (.NOT. lwrad_old) THEN
      radice2 = radice
      radwat2 = radwat
      co2lw   = co2
      n2olw   = n2o
      ch4lw   = ch4
      cfc11lw = cfc11
      cfc12lw = cfc12
      cfc22lw = cfc22
    ELSE
      radice2 = 50.
      radwat2 = 10.
      co2lw = 350.e-6
      n2olw=0.28e-6
      ch4lw=1.75e-6
      cfc11lw=0.3e-9
      cfc12lw=0.5e-9
      cfc22lw=0.2e-9
    END IF
    
    DO k=2,nz
      pl(nplevel+2-k) = pnew(k)       
      ta(nplevel+1-k) = tempnew(k)
      wa(nplevel+1-k) = wetnew(k)
      oa(nplevel+1-k) = ozone(k)       

      cwc(nplevel+1-k,1) = icld(k)
      cwc(nplevel+1-k,2) = wcld(k)	      !cloud water
      cwc(nplevel+1-k,3) = wrain(k)	      !rain water
      fcld(nplevel+1-k)  = frcld(k)	      !cloud fraction

      flx(nplevel+2-k)   = 0.
    ENDDO

    k=nz+1
    pl(nplevel+2-k)=pnew(k-1)*(pnew(k-1)-pl(nplevel+3-k))
    flx(nplevel+2-k)=0.

    ! convert pressure level unit from Pa to hPa   
    do k=nplevel,npdiff+1,-1
       pl(k)=pl(k)/100.0
    enddo
    
    avg_pl = pl(nplevel+1-nz)
    deltap = (avg_pl - 1)/npdiff
    deltap = min(deltap, 4.0)

    ! top boundary layer, k from 6 to 1  
    do k=npdiff,1,-1
!       print*, npdiff, k, pl(k), deltap
       pl(k)=pl(k+1)-deltap
       ta(k)=ta(k+1)
       wa(k)=wa(k+1)
       oa(k)=oa(k+1)
       cwc(k,1)=0.
       cwc(k,2)=0.
       cwc(k,3)=0.
       fcld(k)=0.
       flx(k) =0.
    enddo

    ! set effective cloud particle sizes, 1 = ice, 2 = cloud water,
    !  3 = rain water (micrometres) 
    do k= 1, nplevel-1
       !reff(k,1)  = 80.
       reff(k,1) = radice2
       reff(k,2) = radwat2
       reff(k,3) = 1000.
    enddo

    ! set the cloud droplet (cloud water) effective radius  
!!$    do k=2,nz
!!$ 	  if(watcnew(k).gt.watcri .and. tnumnew(k,1).gt.tnumcri)  then
!!$ 	     reff(nplevel+1-k,2)=min(2.5d-5,max(1.0d-6,		      &
!!$ 	     	  (watcnew(k)/tnumnew(k,1)/rhowat/xpi)**0.33))
!!$ 	     reff(nplevel+1-k,2) =1.0d6*reff(nplevel+1-k,2)!m to um
!!$ 	     reff(nplevel+1-k,2) =1.08*reff(nplevel+1-k,2) !rv-->re
!!$ 	  endif
!!$    enddo

    ! calculate infrared radiation flux  
    call irrad(pl,ta,wa,oa,tb,tg,eg,cwc,reff,fcld,flx,      &
    	       co2lw, n2olw, ch4lw, cfc11lw, cfc12lw,       &
    	       cfc22lw, flxdn_lw_gr_tot)
    print*, 'after irrad'
 
    ! calculate infrared cooling rate	
    do k=npdiff+1,nplevel+1-klow
       heatlwrate(nplevel+1-k)=(flx(k+1)-flx(k))*8.441874/(pl(k+1)-pl(k))
       flxlwtmp(nplevel-k)=-flx(k+1)
    enddo
    k=npdiff
    flxlwtmp(nplevel-k)=-flx(k+1)
         
  end subroutine longwave_radiation_nasa

 subroutine irrad (pl,ta,wa,oa,tb,tg,eg,cwc,reff,fcld,flx,       &
                   co2lw, n2o, ch4, cfc11, cfc12, cfc22,         &
                   flxdn_lw_gr_tot)
    
    real, intent(in)  :: co2lw, n2o, ch4, cfc11, cfc12, cfc22
    real, intent(in)  :: tg, tb
    
    real, dimension(nplevel), intent(in)      :: pl
    real, dimension(nplevel-1), intent(in)    :: ta, wa, oa, fcld
    real, dimension(nplevel-1, 3), intent(in) :: cwc, reff
    real, dimension(10), intent(in):: eg
    
    real, intent(inout) :: flxdn_lw_gr_tot
    real, dimension(nplevel), intent(inout) :: flx

    real, parameter :: T250 = 250.0, C789 = 789.0

    real, dimension(nplevel-1) :: pa, dt250, dp
    real, dimension(nplevel-1) :: dh2o,dco2,do3,dn2o,dch4,dcont,df11,df12,df22
    real, dimension(nplevel-1, 3) :: cwp, taucl
    real, dimension(nplevel) :: transfc
    real, dimension(0:nplevel) ::blayer
    real, dimension(nplevel) :: blevel
    real, dimension(nplevel-1) :: tcldlyr
    real, dimension(nplevel-1, 6)   :: h2oexp, comexp
    real, dimension(nplevel-1, 3)   :: conexp
    real, dimension(nplevel-1, 6,2) :: co2exp
    real, dimension(nplevel-1, 4)   :: n2oexp, ch4exp
    real, dimension(nplevel-1)      :: f11exp, f12exp, f22exp
    real, dimension(nplevel) :: flxu, flxd
    real, dimension(0:nplevel) :: bu, bd
    real, dimension(6)    :: th2o, tcom
    real, dimension(6, 2) :: tco2
    real, dimension(4)    :: tn2o, tch4
    real, dimension(nplevel) :: fclr, trant
    real, dimension(3)    :: tcon
    
    real :: tx, xlayer, bs, rflxs, flag
    real :: tf11, tf12, tf22
    real :: x1, x2, x3

    real :: xx, tauc,reff1, reff2, w1, w2, w3, ww, g1, g2, g3, ff, gg,        &
            a1, b1, fk1, a2, b2, fk2, p1, dwe, dpe, yy

    logical :: oznbnd,co2bnd,h2otbl,conbnd,n2obnd
    logical :: ch4bnd,combnd,f11bnd,f12bnd,f22bnd,b10bnd, h2Oon, co2on, n2Oon

    integer :: k,ne,k1,k2,ik,isb, ibn, klow

    klow=1
     
    ! compute layer pressure (pa) and layer temperature minus 250K (dt250) 
    do k=1,nplevel-1
       pa(k)=0.5 * (pl(k) + pl(k+1))
       dt250(k)=ta(k) - T250
    enddo

    !-------------------------------------------------------------------------!
    !  compute layer absorber amount                                          !
    !                                                                         !
    !     dh2o : water vapor amount (g/cm**2)                                 !
    !     dcont: scaled water vapor amount for continuum absorption (g/cm**2) !
    !     dco2 : co2 amount (cm-atm)stp                                       !
    !     do3  : o3 amount (cm-atm)stp                                        !
    !     dn2o : n2o amount (cm-atm)stp                                       !
    !     dch4 : ch4 amount (cm-atm)stp                                       !
    !     df11 : cfc11 amount (cm-atm)stp                                     !
    !     df12 : cfc12 amount (cm-atm)stp                                     ! 
    !     df22 : cfc22 amount (cm-atm)stp                                     !
    !     the factor 1.02 is equal to 1000/980                                !
    !     factors 789 and 476 are for unit conversion                         !
    !     the factor 0.001618 is equal to 1.02/(.622*1013.25)                 !
    !     the factor 6.081 is equal to 1800/296                               !
    !-------------------------------------------------------------------------!

    do k=1, nplevel-1
       dp   (k) = pl(k+1)-pl(k)
       
       dh2o(k) = 1.02 * wa(k)*dp(k)
       dh2o(k) = max(dh2o(k),1.e-8)
       do3 (k) = 476.0* oa(k)*dp(k)
       do3 (k) = max(do3 (k),1.e-6)
       dco2(k) = C789 * co2lw*dp(k)
       dco2(k) = max(dco2(k),1.e-4)
       
       dch4(k) = C789 *ch4*dp(k)
       dn2o(k) = C789 *n2o*dp(k)
       df11(k) = C789 *cfc11*dp(k)
       df12(k) = C789 *cfc12*dp(k)
       df22(k) = C789 *cfc22*dp(k)

       ! compute scaled water vapor amount for h2o continuum absorption   
       !  following eq. (4.21).  
       xx=pa(k)*0.001618*wa(k)*wa(k)*dp(k)
       dcont(k) = xx*exp(1800.0/ta(k)-6.081)
    enddo

    ! compute layer cloud water amount (gm/m**2)
    !  index is 1 for ice, 2 for waterdrops and 3 for raindrops.
    do  k=1,nplevel + 1 -klow
      xx=1.02*10000.0 *(pl(k+1)-pl(k))
      cwp(k,1)=xx*cwc(k,1)
      cwp(k,2)=xx*cwc(k,2)
      cwp(k,3)=xx*cwc(k,3)
    enddo

    ! The surface (nplevel) is treated as a layer filled with black clouds.  !
    ! transfc is the transmittance between the surface and a pressure level. !
    ! trantcr is the clear-sky transmittance between the surface and a       !
    ! pressure level.                                                        !
    transfc(nplevel-klow+2)=1.

    ! initialize fluxes  
    do k=1,nplevel
       flx(k) = 0.
    enddo
    
    ! reset total ground irradiance for x-line 
    flxdn_lw_gr_tot = 0.
  
    !-----------------------------------------------------------------------!
    ! integration over spectral bands                                       !
    !                                                                       !
    !     if h2otbl, compute h2o (line) transmittance using table look-up.  !
    !     if conbnd, compute h2o (continuum) transmittance in bands 2-7.    !
    !     if co2bnd, compute co2 transmittance in band 3.                   !
    !     if oznbnd, compute  o3 transmittance in band 5.                   !
    !     if n2obnd, compute n2o transmittance in bands 6 and 7.            !
    !     if ch4bnd, compute ch4 transmittance in bands 6 and 7.            !
    !     if combnd, compute co2-minor transmittance in bands 4 and 5.      ! 
    !     if f11bnd, compute cfc11 transmittance in bands 4 and 5.          !
    !     if f12bnd, compute cfc12 transmittance in bands 4 and 6.          !
    !     if f22bnd, compute cfc22 transmittance in bands 4 and 6.          !
    !     if b10bnd, compute flux reduction due to n2o in band 10.          !
    !-----------------------------------------------------------------------!

    print*, 'start loop 1000'
    do 1000 ibn=1,10
       h2Oon=.true.
       co2on=.true.
       n2Oon=.true.
       h2otbl=.false.
       conbnd=ibn.ge.2.and.ibn.le.7
       co2bnd=ibn.eq.3
       oznbnd=ibn.eq.5
       n2obnd=ibn.eq.6.or.ibn.eq.7
       ch4bnd=ibn.eq.6.or.ibn.eq.7
       combnd=ibn.eq.4.or.ibn.eq.5
       f11bnd=ibn.eq.4.or.ibn.eq.5
       f12bnd=ibn.eq.4.or.ibn.eq.6
       f22bnd=ibn.eq.4.or.ibn.eq.6
       b10bnd=ibn.eq.10

       ! blayer is the spectrally integrated planck flux of the mean layer 
       ! temperature derived from eq. (3.11) 
       ! The fitting for the planck flux is valid for the range 160-345 K. 
       do k=1,nplevel-1
          tx=ta(k)
          call planck(ibn,tx,xlayer)
          blayer(k)=xlayer
       enddo

       ! Index "0" is the layer above the top of the atmosphere.
       blayer(0)=0.
         
       ! Surface emission and reflectivity. See Section 9.  
       call sfcflux (ibn,tg,eg,bs,rflxs) 
       blayer(nplevel-klow+2)=bs

       ! interpolate Planck function at model levels (linear in p) 
       do  k=2,nplevel-klow+1
          blevel(k)=(blayer(k-1)*dp(k)+blayer(k)*dp(k-1))/(dp(k-1)+dp(k))
       enddo

       ! Extrapolate blevel(i,1) from blayer(i,2) and blayer(i,1)    
       blevel(1)=blayer(1)+(blayer(1)-blayer(2))*dp(1)/(dp(1)+dp(2))
         
       ! If the surface air temperature tb is known, compute blevel(i,nplevel)
       call planck(ibn,tb,xlayer)
       blevel(nplevel-klow+2)=xlayer

       !---------------------------------------------------------------------!
       ! Otherwise, extrapolate blevel(nplevel) from blayer(nplevel-2)       !
       ! and blayer                                                          !
       ! (nplevel-1)                                                         !
       !       do i=nxl,nxr                                                  !
       !           blevel(i,nplevel)=blayer(i,nplevel-1)               &     !
       !                        +(blayer(i,nplevel-1)-blayer(i,nplevel-2))   &    !
       !                        *dp(i,nplevel-1)/(dp(i,nplevel-1)+dp(i,nplevel-2))!
       !       enddo                                                         !
       !---------------------------------------------------------------------!

       ! Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
       ! Rain optical thickness is set to 0.00307 /(gm/m**2).  
       ! It is for a specific drop size distribution provided by Q. Fu. 
       do  k=1,nplevel-klow+1
          taucl(k,1)=cwp(k,1)*(aib(1,ibn)+aib(2,ibn)  	       &
        	        /reff(k,1)**aib(3,ibn))
          taucl(k,2)=cwp(k,2)*(awb(1,ibn)+(awb(2,ibn) 	       &
        	        +(awb(3,ibn)+awb(4,ibn)*reff(k,2))*reff(k,2))&
        	        *reff(k,2))
          taucl(k,3)=0.00307*cwp(k,3)
       enddo

       !----------------------------------------------------------------------
       ! Compute cloud single-scattering albedo and asymmetry factor for      !
       !     a mixture of ice particles and liquid drops following            !
       !     Eqs. (6.5), (6.6), (6.9) and (6.10).                             !
       !     Single-scattering albedo and asymmetry factor of rain are set    !
       !     to 0.54 and 0.95, respectively, based on the information provided!
       !     by Prof. Qiang Fu.                                               !
       !----------------------------------------------------------------------!

       do  k=1,nplevel-klow+1 

          tcldlyr(k) = 1.
          tauc=taucl(k,1)+taucl(k,2)+taucl(k,3)
          if (tauc.gt.0.02 .and. fcld(k).gt.0.01) then
             reff1=min(reff(k,1),130.0*1.)
             reff2=min(reff(k,2),20.0*1.)
             
             w1=taucl(k,1)*(aiw(1,ibn)+(aiw(2,ibn)+(aiw(3,ibn)      &
        	    +aiw(4,ibn)*reff1)*reff1)*reff1)
             w2=taucl(k,2)*(aww(1,ibn)+(aww(2,ibn)+(aww(3,ibn)      &
        	    +aww(4,ibn)*reff2)*reff2)*reff2) 
             w3=taucl(k,3)*0.54
             ww=(w1+w2+w3)/tauc

             g1=w1*(aig(1,ibn)+(aig(2,ibn)+(aig(3,ibn)  	      &
        	    +aig(4,ibn)*reff1)*reff1)*reff1)
             g2=w2*(awg(1,ibn)+(awg(2,ibn)+(awg(3,ibn)  	      &
        	    +awg(4,ibn)*reff2)*reff2)*reff2)
             g3=w3*0.95

             gg=(g1+g2+g3)/(w1+w2+w3)

             ! Parameterization of LW scattering following Eqs. (6.11) and (6.12)
             ff=0.5+(0.3739+(0.0076+0.1185*gg)*gg)*gg
             tauc=(1.-ww*ff)*tauc

             ! compute cloud diffuse transmittance. It is approximated  
             ! by using a diffusivity factor of 1.66.
             tcldlyr(k)=exp(-1.66*tauc)
          endif
            
       enddo

       ! Compute the exponential terms (Eq. 8.21) at each layer due to 
       ! water vapor line absorption when k-distribution is used  
       if (.not.h2otbl .and. .not.b10bnd .and. h2Oon) then
          call h2oexps(ibn,dh2o,pa,dt250,h2oexp)
       endif

       ! compute the exponential terms (Eq. 4.24) at each layer due to   
       ! water vapor continuum absorption. 
       ! ne is the number of terms used in each band to compute water  
       ! vapor continuum transmittance (Table 9). 
       print*, 'before calc exponentials'
       ne=0
       if (conbnd) then
          ne=1
          if (ibn.eq.3) ne=3
          call conexps(ibn,dcont,conexp)
       endif

       ! compute the exponential terms (Eq. 8.21) at each layer due to   
       ! co2 absorption  
       if (co2bnd) then
          call co2exps(dco2,pa,dt250,co2exp)
       endif

       ! compute the exponential terms at each layer due to n2o absorption
       if (n2obnd) then
          call n2oexps(ibn,dn2o,pa,dt250,n2oexp)
       endif

       ! compute the exponential terms at each layer due to ch4 absorption
       if (ch4bnd) then
          call ch4exps(ibn,dch4,pa,dt250,ch4exp)
       endif

       ! Compute the exponential terms due to co2 minor absorption     
       if (combnd) then
          call comexps(ibn,dco2,dt250,comexp)
       endif

       ! Compute the exponential terms due to cfc11 absorption.
       ! The values of the parameters are given in Table 7.  
       if (f11bnd) then
          a1  = 1.26610e-3
          b1  = 3.55940e-6
          fk1 = 1.89736e+1
          a2  = 8.19370e-4
          b2  = 4.67810e-6
          fk2 = 1.01487e+1
          call cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,df11,dt250,f11exp)
       endif

       ! Compute the exponential terms due to cfc12 absorption.
       if (f12bnd) then
          a1  = 8.77370e-4
          b1  =-5.88440e-6
          fk1 = 1.58104e+1
          a2  = 8.62000e-4
          b2  =-4.22500e-6
          fk2 = 3.70107e+1
          call cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,df12,dt250,f12exp)
       endif

       ! Compute the exponential terms due to cfc22 absorption.  
       if (f22bnd) then
          a1  = 9.65130e-4
          b1  = 1.31280e-5
          fk1 = 6.18536e+0
          a2  =-3.00010e-5
          b2  = 5.25010e-7
          fk2 = 3.27912e+1
          call cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,df22,dt250,f22exp)
       endif

       ! Compute the exponential terms at each layer in band 10 due to 
       !     h2o line and continuum, co2, and n2o absorption   
       if (b10bnd) then
          call b10exps(dh2o,dcont,dco2,dn2o,pa,dt250,       &
                        h2oexp,conexp,co2exp,n2oexp, h2Oon, co2on, n2Oon)
       endif

       ! initialize fluxes                                                !
       do k=1,nplevel
         flxu(k) = 0.
         flxd(k) = 0.
       enddo

       ! For a given level, k1, compute the transmittance between this level  !
       ! all the levels below, trant(i,k2).                                   !
       ! Also, compute the upward and doward blackbody emissions of a layer,  !
       ! bu and bd.                                                           !

       print*, 'before loop 2000'
       bd(0) = 0.
       do 2000 k1=1,nplevel-1
          if (h2Oon .and. .not. h2otbl) then
             do ik=1,6
                th2o(ik)=1.
             enddo
          endif
          tcon(1:3) = 0.
          if (conbnd) then
             do ik=1,3
                tcon(ik)=1.
             enddo
          endif
          if (co2bnd) then
             do isb=1,2
                do ik=1,6
                   tco2(ik,isb)=1.
                enddo
             enddo
          endif
          if (n2obnd) then
             do ik=1,4
                tn2o(ik)=1.
             enddo
          endif
          if (ch4bnd) then
             do ik=1,4
                tch4(ik)=1.
             enddo
          endif
          if (combnd) then
             do ik=1,6
                tcom(ik)=1.
             enddo
          endif
          if (f11bnd) then
                tf11=1.
          endif
          if (f12bnd) then
             tf12=1.
          endif
          if (f22bnd) then
             tf22=1.
          endif
          if (b10bnd) then
             do ik=1,5
                th2o(ik)=1.
             enddo
             do ik=1,6
                tco2(ik,1)=1.
             enddo
             tcon(1)=1.
                          
             do ik=1,2
                tn2o(ik)=1.
             enddo
          endif

          x1=0.
          x2=0.
          x3=0.
          do k=1,nplevel
             fclr(k)=1.
          enddo

          ! loop over the bottom level of the region (k2)    
          do 3000 k2=k1+1,nplevel
             trant(k2)=1.
             
             if (h2otbl) then
                w1=-8.0
                p1=-2.0
                dwe=0.3
                dpe=0.2

                if (h2Oon) then
                   if (ibn.eq.1) then
                      call tablup(k2,nxp,nh,dh2o,pa,dt250,x1,x2,x3,   &
                           w1,p1,dwe,dpe,nh,h11,h12,h13,trant)
                   endif
                   
                   if (ibn.eq.2) then
                      call tablup(k2,nxp,nh,dh2o,pa,dt250,x1,x2,x3,       &
                           w1,p1,dwe,dpe,nh,h21,h22,h23,trant)
                   endif
                   if (ibn.eq.8) then
                      call tablup(k2,nxp,nh,dh2o,pa,dt250,x1,x2,x3,       &
                           w1,p1,dwe,dpe,nh,h81,h82,h83,trant)
                   endif
                endif

                if (conbnd) then
                   tcon(1)=tcon(1)*conexp(k2-1,1)
                   trant(k2)=trant(k2)*tcon(1)
                endif
             else
                if (h2Oon .and. .not.b10bnd) then
                   call h2okdis(ibn,k2-1,ne,h2oexp,conexp,th2o,tcon,trant)
                endif
             endif

             if (co2bnd) then
                call co2kdis(k2-1,co2exp,tco2,trant)
             endif

             if (oznbnd) then
                w1=-6.0
                p1=-2.0
                dwe=0.3
                dpe=0.2
                call tablup(k2,nxp,nh,do3,pa,dt250,x1,x2,x3,            &
                     w1,p1,dwe,dpe,no,o1,o2,o3,trant)
             endif

             if (n2obnd) then
                call n2okdis(ibn,k2-1,n2oexp,tn2o,trant)
             endif

             if (ch4bnd) then
                call ch4kdis(ibn,k2-1,ch4exp,tch4,trant)
             endif

             if (combnd) then
                call comkdis(ibn,k2-1,comexp,tcom,trant)
             endif

             if (f11bnd) then
                call cfckdis(k2-1,f11exp,tf11,trant)
             endif

             if (f12bnd) then
                call cfckdis(k2-1,f12exp,tf12,trant)
             endif

             if (f22bnd) then
                call cfckdis(k2-1,f22exp,tf22,trant)
             endif

             ! Compute transmittance in band 10 using k-distribution method.  
             ! For band 10, trant is the change in transmittance due to n2o absorption   !
             if (b10bnd) then
                call b10kdis(k2-1,h2oexp,conexp,co2exp,n2oexp,        &
                             th2o,tcon,tco2,tn2o,trant, h2Oon,co2on, n2Oon)

             endif
             fclr(k2)=fclr(k2-1)*tcldlyr(k2-1)

             ! Compute upward and downward blackbody emission of a layer
             if (k1.eq.1) then
                xx=(blayer(k2-1)-blevel(k2-1)) * (blayer(k2-1)-blevel(k2))

                if (xx.gt.0.0) then
                   ! If xx>0, there is a local temperature minimum or maximum.  
                   ! Computations of bd and bu follow Eq. (8.20).  
                   bd(k2-1)=0.5*blayer(k2-1)+0.25*(blevel(k2-1)+blevel(k2))
                   bu(k2-1)=bd(k2-1)
                else
                   ! Computations of bd and bu following Eqs.(8.17) and (8.18).  
                   ! The effect of clouds on the transmission of a layer is taken 
                   ! into account, following Eq. (8.19).	   
                   xx=(fcld(k2-1)*tcldlyr(k2-1)+(1.-fcld(k2-1))) * trant(k2)
                   yy=min(0.9999,xx)
                   yy=max(0.00001,yy)
                   xx=(blevel(k2-1)-blevel(k2))/log(yy)
                   bd(k2-1)=(blevel(k2)-blevel(k2-1)*yy)/(1.-yy)-xx
                   bu(k2-1)=(blevel(k2-1)+blevel(k2))-bd(k2-1)
                endif
                bu(nplevel-klow+2)=blayer(nplevel-klow+2)
 
             endif
3000      end do

          ! upward and downward flux calculations. 
          do 4000 k2=k1+1,nplevel
             if (k2.eq.k1+1 .and. ibn .ne. 10) then
                ! The first terms on the rhs of Eqs. (8.15) and (8.16) 
                flxu(k1)=flxu(k1)-bu(k1)
                flxd(k2)=flxd(k2)+bd(k1)
             endif

             ! The summation terms on the rhs of Eqs. (8.15) and (8.16).
             ! Also see Eqs. (5.4) and (5.5) for Band 10. 
             xx=trant(k2)*(bu(k2-1)-bu(k2))
             flxu(k1)=flxu(k1)+xx*fclr(k2)
             xx=trant(k2)*(bd(k1-1)-bd(k1))
             flxd(k2)=flxd(k2)+xx*fclr(k2)
4000      end do

          ! Here, fclr and trant are, respectively, the clear line-of-sight 
          ! and the transmittance between k1 and the surface. 
          transfc(k1) =trant(nplevel-klow+2)*fclr(nplevel-klow+2)

2000   end do

       if (.not. b10bnd) then
          ! For surface emission. 
          ! Note: blayer(i,nplevel) and dbs include the surface emissivity effect.  
          flxu(nplevel-klow+2)=-blayer(nplevel-klow+2)

          ! Add the flux reflected by the surface.(Second term on the rhs of Eq. 8.16)
          do k=1,nplevel-klow+1    
             flxu(k)=flxu(k)- flxd(nplevel-klow+2)*transfc(k)*rflxs  
          enddo
       endif

       ! Summation of fluxes over spectral bands 
       do k=1,nplevel-klow+2
          flx(k)=flx(k)+flxd(k)+flxu(k)

          ! extract and add total ground irradiance (for each band)  
          if (k==nplevel-klow+2) flxdn_lw_gr_tot = flxdn_lw_gr_tot + flxd(k)
       enddo

1000 end do
      
  end subroutine irrad

  subroutine planck(ibn,t,xlayer)

    ! Compute spectrally integrated Planck flux   
    integer, intent(in) :: ibn          ! spectral band index
    real, intent(in) :: t               ! temperature (K)
    real, intent(inout) ::  xlayer      ! planck flux (w/m2)

    xlayer=(t*(t*(t*(t*(t*cb(6,ibn)+cb(5,ibn))  		&
          +cb(4,ibn))+cb(3,ibn))+cb(2,ibn)))+cb(1,ibn)

  end subroutine planck

  subroutine sfcflux (ibn,tg,eg,bs,rflxs)

    !------------------------------------------------------------------------!
    ! Compute emission and reflection by an inhomogeneous surface with       !
    ! vegetation cover.                                                      !
    !                                                                        !
    ! Input parameters                                                       !
    ! ibn   = index for the spectral band                                    !
    ! nx    = number of grid boxs                                            !
    ! ns    = number of sub-grid box                                         !
    ! fs    = fractional cover of sub-grid box                               !
    ! tg    = sub-grid ground temperature                                    !
    ! eg    = sub-grid ground emissivity                                     !
    ! tv    = sub-grid vegetation temperature                                !
    ! ev    = sub-grid vegetation emissivity                                 !
    ! rv    = sub-grid vegetation reflectivity (rv)                          !
    ! if there is vegetation cover, vege=.true.                              !
    !                                                                        !
    ! Output parameters                                                      !
    ! bs    = Emission by the surface (ground+vegetation)                    !
    ! dbs   = Derivative of bs rwt temperature                               !
    ! rflxs = Reflection by the surface                                      !
    !------------------------------------------------------------------------!

    integer, intent(in) :: ibn
    real, intent(in) :: tg,eg(10)
    real, intent(inout) ::  bs,rflxs
    real :: tx

    ! For homogeneous surface without vegetation following Eqs. (9.4), (9.5), and (3.13) !
    tx=tg
    call planck(ibn,tx,bs)
    rflxs=1.-eg(ibn)
       
  end subroutine sfcflux

  subroutine h2oexps(ibn,dh2o,pa,dt250,h2oexp)

    !-------------------------------------------------------------------------!
    ! Compute exponentials for water vapor line absorption                    !
    ! in individual layers using Eqs. (8.21) and (8.22).                      !
    !                                                                         !
    ! Input parameters                                                        !
    ! dh2o    = layer water vapor amount for line absorption                  !
    ! pa      = layer pressure                                                !
    ! dt250   = layer temperature minus 250K                                  !
    ! xkw     = absorption coefficients for the first k-distribution function !
    !             due to h2o line absorption                                  !
    ! aw, bw, pm =  coefficients for the temperature and pressure scaling     ! 
    ! mw      =  ratios between neighboring absorption coefficients for       !
    !            h2o line absorption                                          !
    !                                                                         !
    ! Output parameters                                                       !
    ! h2oexp  = exponentials for each layer                                   !
    !-------------------------------------------------------------------------!

    integer, intent(in) :: ibn
    real, intent(in) :: dh2o(nplevel-1), pa(nplevel-1), dt250(nplevel-1)
    real, intent(inout), dimension(nplevel-1, 6) :: h2oexp

    real, parameter :: pref = 500.
    integer :: k, ik, klow
    real :: xh

    ! Note that the 3 sub-bands in band 3 use the same set of xkw, aw, 
    ! and bw,  therefore, h2oexp for these sub-bands are identical.  

    klow = 1
    do k = 1, nplevel-klow+1
       ! xh is the scaled water vapor amount for line absorption computed from Eq. (4.4). 
       xh = dh2o(k)*(pa(k)/pref)**pm(ibn)	      &
    	    *(1.+(aw(ibn)+bw(ibn)* dt250(k))*dt250(k))

       ! h2oexp is the water vapor transmittance of the layer k due to line absorption  
       h2oexp(k,1) = exp(-xh*xkw(ibn))
    enddo

    ! compute transmittances from Eq. (8.22)                                    !
    do ik=2,6

       if (mw(ibn).eq.6) then
          do k = 1, nplevel-klow+1
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            h2oexp(k,ik) = xh*xh*xh
          enddo
       elseif (mw(ibn).eq.8) then
          do k = 1, nplevel-klow+1
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            xh = xh*xh
            h2oexp(k,ik) = xh*xh
          enddo
       elseif (mw(ibn).eq.9) then
          do k = 1, nplevel-klow+1
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)*h2oexp(k,ik-1)
            h2oexp(k,ik) = xh*xh*xh
          enddo
       else
          do k = 1, nplevel-klow+1
            xh = h2oexp(k,ik-1)*h2oexp(k,ik-1)
            xh = xh*xh
            xh = xh*xh
            h2oexp(k,ik) = xh*xh
          enddo
       endif
       
    enddo

  end subroutine h2oexps

  subroutine conexps(ibn,dcont,conexp)

    !-------------------------------------------------------------------------------!
    ! Compute exponentials for continuum absorption in individual layers.           !
    !                                                                               !
    ! Input parameters                                                              !
    ! dcont = layer scaled water vapor amount for continuum absorption              !
    ! xke   = absorption coefficients for the first k-distribution function due to  !
    !         water vapor continuum absorption (xke)                                !
    !                                                                               !
    ! Output parameters                                                             !
    ! conexp = 1 or 3 exponentials for each layer                                   !
    !-------------------------------------------------------------------------------!

    real, intent(in) :: dcont(nplevel-1)
    integer, intent(in) :: ibn
    real, intent(inout) :: conexp(nplevel-1,3)

    integer :: i, k, klow
      
    klow = 1
    do k = 1, nplevel-klow+1
       conexp(k,1) = exp(-dcont(k)*xke(ibn))
    enddo

    if (ibn .eq. 3) then
       ! The absorption coefficients for sub-bands 3b and 3a are, respectively,   
       ! two and four times the absorption coefficient for sub-band 3c (Table 9).   
       ! Note that conexp(i,k,3) is for sub-band 3a. 
       klow = 1
       do k = 1, nplevel-klow+1
          conexp(k,2) = conexp(k,1) *conexp(k,1)
          conexp(k,3) = conexp(k,2) *conexp(k,2)
       enddo
    endif

  end subroutine conexps

  subroutine co2exps(dco2,pa,dt250,co2exp)

    !----------------------------------------------------------------------------!
    ! Compute co2 exponentials for individual layers.                            !
    !                                                                            !
    ! Input parameters                                                           !
    ! dco2  = layer co2 amount                                                   !
    ! pa    = layer pressure                                                     !
    !                                                                            !
    ! Output parameters                                                          !
    ! co2exp = 6 exponentials for each layer                                     !
    !----------------------------------------------------------------------------!
      
    real, intent(in) ::  dco2(nplevel-1), pa(nplevel-1), dt250(nplevel-1)
    real, intent(inout) :: co2exp(nplevel-1,6,2)
    
    integer :: k, klow
    real :: xc

    klow = 1
    do k = 1, nplevel-klow+1
       ! The scaling parameters are given in Table 3, and values of the absorption !
       ! coefficient are given in Table 10.					   !

       ! Scaled co2 amount for band-wings (sub-bands 3a and 3c) 
       xc = dco2(k)*((pa(k)/300.0)**0.5)*(1.+(0.0182+1.07e-4*dt250(k))*dt250(k))

       ! six exponentials by powers of 8 (See Eqs. 8.21, 8.22 and Table 10).
       co2exp(k,1,1)=exp(-xc*2.656e-5)
       
       xc=co2exp(k,1,1)*co2exp(k,1,1)
       xc=xc*xc
       co2exp(k,2,1)=xc*xc

       xc=co2exp(k,2,1)*co2exp(k,2,1)
       xc=xc*xc
       co2exp(k,3,1)=xc*xc

       xc=co2exp(k,3,1)*co2exp(k,3,1)
       xc=xc*xc
       co2exp(k,4,1)=xc*xc

       xc=co2exp(k,4,1)*co2exp(k,4,1)
       xc=xc*xc
       co2exp(k,5,1)=xc*xc

       xc=co2exp(k,5,1)*co2exp(k,5,1)
       xc=xc*xc
       co2exp(k,6,1)=xc*xc

       ! For band-center region (sub-band 3b)	
       xc = dco2(k)*(pa(k)/30.0)**0.85*(1.+(0.0042+2.00e-5*dt250(k))*dt250(k))

       co2exp(k,1,2)=exp(-xc*2.656e-3)
    	 
       xc=co2exp(k,1,2)*co2exp(k,1,2)
       xc=xc*xc
       co2exp(k,2,2)=xc*xc
    	 
       xc=co2exp(k,2,2)*co2exp(k,2,2)
       xc=xc*xc
       co2exp(k,3,2)=xc*xc

       xc=co2exp(k,3,2)*co2exp(k,3,2)
       xc=xc*xc
       co2exp(k,4,2)=xc*xc

       xc=co2exp(k,4,2)*co2exp(k,4,2)
       xc=xc*xc
       co2exp(k,5,2)=xc*xc

       xc=co2exp(k,5,2)*co2exp(k,5,2)
       xc=xc*xc
       co2exp(k,6,2)=xc*xc
    enddo

  end subroutine co2exps

  subroutine n2oexps(ibn,dn2o,pa,dt250,n2oexp)

    !------------------------------------------------------------------------!
    ! Compute n2o exponentials for individual layers                         !
    !                                                                        !
    ! Output parameters: n2oexp  =  2 or 4 exponentials for each layer       !
    !------------------------------------------------------------------------!
      
    integer, intent(in) :: ibn
    real, intent(in) ::  dn2o(nplevel-1), pa(nplevel-1), dt250(nplevel-1)
    real, intent(inout) :: n2oexp(nplevel-1,4)
    
    integer :: k, klow
    real :: xc, xc1, xc2

    ! Scaling and absorption data are given in Table 5.  
    ! Transmittances are computed using Eqs. (8.21) and (8.22).  
    klow = 1
    do k = 1, nplevel-klow+1
       ! Four exponential by powers of 21 for band 6.	
       if (ibn .eq. 6) then
    	  xc=dn2o(k)*(1.+(1.9297e-3+4.3750e-6*dt250(k))*dt250(k))
    	  n2oexp(k,1)=exp(-xc*6.31582e-2)
    	  
    	  xc=n2oexp(k,1)*n2oexp(k,1)*n2oexp(k,1)
    	  xc1=xc*xc
    	  xc2=xc1*xc1
    	  n2oexp(k,2)=xc*xc1*xc2	  
       else
    	  ! four exponential by powers of 8 for band 7 
    	  xc=dn2o(k)*(pa(k)/500.0)**0.48 *    &
    		 (1.+(1.3804e-3+7.4838e-6*dt250(k))*dt250(k))
    	  n2oexp(k,1)=exp(-xc*5.35779e-2)

    	  xc=n2oexp(k,1)*n2oexp(k,1)
    	  xc=xc*xc
    	  n2oexp(k,2)=xc*xc
    	  xc=n2oexp(k,2)*n2oexp(k,2)
    	  xc=xc*xc
    	  n2oexp(k,3)=xc*xc
    	  xc=n2oexp(k,3)*n2oexp(k,3)
    	  xc=xc*xc
    	  n2oexp(k,4)=xc*xc
       endif
    enddo

  end subroutine n2oexps

  subroutine ch4exps(ibn,dch4,pa,dt250,ch4exp)

    !-------------------------------------------------------------------------!
    ! Compute ch4 exponentials for individual layers                          !
    !                                                                         !
    ! Output parameters : ch2exp = 1 or 4 exponentials for each layer         !
    !-------------------------------------------------------------------------!
      
    integer, intent(in) :: ibn
    real, intent(in) ::  dch4(nplevel-1), pa(nplevel-1), dt250(nplevel-1)
    real, intent(inout) :: ch4exp(nplevel-1,4)
    
    integer :: k, klow
    real :: xc

    ! Scaling and absorpton data are given in Table 5   
    klow = 1
    do k = 1, nplevel-klow+1
       ! four exponentials for band 6 
       if (ibn .eq. 6) then
    	  xc=dch4(k)*(1.+(1.7007e-2+1.5826e-4*dt250(k))*dt250(k))
    	  ch4exp(k,1)=exp(-xc*5.80708e-3)
       else
    	  ! four exponentials by powers of 12 for band 7   
    	  xc=dch4(k)*(pa(k)/500.0)**0.65	    &
    		   *(1.+(5.9590e-4-2.2931e-6*dt250(k))*dt250(k))
    	  ch4exp(k,1)=exp(-xc*6.29247e-2)

    	  xc=ch4exp(k,1)*ch4exp(k,1)*ch4exp(k,1)
    	  xc=xc*xc
    	  ch4exp(k,2)=xc*xc

    	  xc=ch4exp(k,2)*ch4exp(k,2)*ch4exp(k,2)
    	  xc=xc*xc
    	  ch4exp(k,3)=xc*xc

    	  xc=ch4exp(k,3)*ch4exp(k,3)*ch4exp(k,3)
    	  xc=xc*xc
    	  ch4exp(k,4)=xc*xc
       endif
    enddo

  end subroutine ch4exps

  subroutine comexps(ibn,dcom,dt250,comexp)
    
    !------------------------------------------------------------------------!
    ! Compute co2-minor exponentials for individual layers using Eqs. (8.21) !
    ! and (8.22).                                                            !
    !                                                                        !
    ! Output parameters : comexp =  6 exponentials for each layer            !
    !------------------------------------------------------------------------!

    integer, intent(in) :: ibn
    real, intent(in) ::  dcom(nplevel-1), dt250(nplevel-1)
    real, intent(inout) :: comexp(nplevel-1,6)
    
    integer :: k, ik, klow 
    real :: xc

    ! Scaling and absorption data are given in Table 6 
    klow = 1
    do k = 1, nplevel-klow+1
       if (ibn .eq. 4) then
    	  xc=dcom(k)*(1.+(3.5775e-2 + 4.0447e-4 * dt250(k))*dt250(k))
       endif

       if (ibn.eq.5) then
    	  xc=dcom(k)*(1.+(3.4268e-2 + 3.7401e-4 * dt250(k))*dt250(k))
       endif

       comexp(k,1)=exp(-xc*1.922e-7)

       do ik=2,6
    	  xc=comexp(k,ik-1)*comexp(k,ik-1)
    	  xc=xc*xc
    	  comexp(k,ik)=xc*comexp(k,ik-1)
       enddo
    enddo

  end subroutine comexps

  subroutine cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,dcfc,dt250,cfcexp)

    !-------------------------------------------------------------------------!
    ! compute cfc(-11, -12, -22) exponentials for individual layers.          !
    !                                                                         !
    ! Input parameters                                                        !
    ! a1,b1,a2,b2 = parameters for computing the scaled cfc amounts           !
    !                  for temperature scaling                                !
    ! fk1, fk2    = the absorption coefficients for the first k-distribution  !
    !                  function due to cfcs                                   !
    !                                                                         !
    ! Output parameters: cfcexp = 1 exponential for each layer                !
    !-------------------------------------------------------------------------!

    integer, intent(in) :: ibn
    real, intent(in) ::  dcfc(nplevel-1), dt250(nplevel-1)
    real, intent(in) :: a1, b1, fk1, a2, b2, fk2
    real, intent(inout) :: cfcexp(nplevel-1)
    
    integer :: k, klow
    real :: xf

    klow = 1
    do k = 1, nplevel-klow+1
       ! compute the scaled cfc amount (xf) and exponential (cfcexp)
       if (ibn .eq. 4) then
    	  xf=dcfc(k)*(1.+(a1+b1*dt250(k))*dt250(k))
    	  cfcexp(k)=exp(-xf*fk1)
       else
    	  xf=dcfc(k)*(1.+(a2+b2*dt250(k))*dt250(k))
    	  cfcexp(k)=exp(-xf*fk2)
       endif

    enddo

  end subroutine cfcexps

  subroutine b10exps(dh2o,dcont,dco2,dn2o,pa,dt250,h2oexp,conexp,co2exp,n2oexp, h2Oon,  &
                     co2on, n2Oon)

    !------------------------------------------------------------------------!
    ! Compute band3a exponentials for individual layers                      !
    !                                                                        !
    ! Output parameters: h2oexp,conexp,co2exp,n2oexp = exponentials for each !
    !layer                                                                   !
    !------------------------------------------------------------------------!

    real, intent(in), dimension(nplevel-1) :: dh2o, dcont, dco2, dn2o
    real, intent(in), dimension(nplevel-1) :: pa, dt250
    real, intent(inout) :: h2oexp(nplevel-1,6), conexp(nplevel-1,3),  &
                                  co2exp(nplevel-1,6,2), n2oexp(nplevel-1,4)
    logical, intent(in) :: h2Oon, co2on, n2Oon

    integer :: k, klow
    real :: xx, xx1, xx2, xx3
     
    klow = 1
    do k = 1, nplevel-klow+1
      
      ! Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3). 
      if (h2Oon) then
         xx=dh2o(k)*(pa(k)/500.0) * (1.+(0.0149+6.20e-5*dt250(k))*dt250(k))

         ! six exponentials by powers of 8
         h2oexp(k,1)=exp(-xx*0.10624)

         xx=h2oexp(k,1)*h2oexp(k,1)
         xx=xx*xx
         h2oexp(k,2)=xx*xx

         xx=h2oexp(k,2)*h2oexp(k,2)
         xx=xx*xx
         h2oexp(k,3)=xx*xx

         xx=h2oexp(k,3)*h2oexp(k,3)
         xx=xx*xx
         h2oexp(k,4)=xx*xx

         xx=h2oexp(k,4)*h2oexp(k,4)
         xx=xx*xx
         h2oexp(k,5)=xx*xx

         ! One exponential of h2o continuum for sub-band 3a (Table 9).  
         conexp(k,1)=exp(-dcont(k)*109.0)
      endif

      ! Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).	
      if (co2on) then
         xx=dco2(k)*(pa(k)/300.0)**0.5 * (1.+(0.0179+1.02e-4*dt250(k))*dt250(k))

         ! six exponentials by powers of 8	
         co2exp(k,1,1)=exp(-xx*2.656e-5)

         xx=co2exp(k,1,1)*co2exp(k,1,1)
         xx=xx*xx
         co2exp(k,2,1)=xx*xx

         xx=co2exp(k,2,1)*co2exp(k,2,1)
         xx=xx*xx
         co2exp(k,3,1)=xx*xx

         xx=co2exp(k,3,1)*co2exp(k,3,1)
         xx=xx*xx
         co2exp(k,4,1)=xx*xx

         xx=co2exp(k,4,1)*co2exp(k,4,1)
         xx=xx*xx
         co2exp(k,5,1)=xx*xx

         xx=co2exp(k,5,1)*co2exp(k,5,1)
         xx=xx*xx
         co2exp(k,6,1)=xx*xx
      endif

      ! Compute the scaled n2o amount for Band 10 (Table 5).
      if (n2Oon) then
         xx=dn2o(k)*(1.+(1.4476e-3 + 3.6656e-6 * dt250(k))*dt250(k))

         ! Two exponentials by powers of 58   
         n2oexp(k,1)=exp(-xx*0.25238)

         xx=n2oexp(k,1)*n2oexp(k,1)
         xx1=xx*xx
         xx1=xx1*xx1
         xx2=xx1*xx1
         xx3=xx2*xx2
         n2oexp(k,2)=xx*xx1*xx2*xx3
      endif
      
    enddo

  end subroutine b10exps

  subroutine tablup(k2,nxp,nh,dw,pa,dt250,x1,x2,x3,w1,p1,          &
                    dwe,dpe,ncoef,coef1,coef2,coef3,trant)

    !-----------------------------------------------------------------------------!
    ! Compute water vapor, co2 and o3 transmittances between level k1 and and     !
    ! level k2 for m soundings, using table look-up.                              !
    ! Calculations follow Eq. (4.16).                                             !
    !                                                                             !
    ! Input variables:                                                            !
    ! k2  = index for level                                                       !
    ! nxp = number of pressure intervals in the table                             !
    ! nh  = number of absorber amount intervals in the table                      !
    ! dw  = layer absorber amount                                                 !
    ! w1  = first value of absorber amount (log10) in the table                   !
    ! p1  = first value of pressure (log10) in the table                          !
    ! dwe = size of the interval of absorber amount (log10) in the table          !
    ! dpe = size of the interval of pressure (log10) in the table                 !
    ! coef1, coef2, coef3 = pre-computed coefficients                             !
    ! Output variables:                                                           !
    ! x1 = column integrated absorber amount                                      !
    ! x2 = absorber-weighted column pressure                                      !
    ! x3 = absorber-weighted column temperature                                   !
    ! trant = transmittance                                                       !
    !   Note: Units of x1: g/cm**2 for water vapor and (cm-atm)stp for co2 and o3.!
    !-----------------------------------------------------------------------------!
 
    integer, intent(in) :: k2, nxp, nh, ncoef
    real, intent(in) :: w1, p1, dwe, dpe
    real, intent(in), dimension(nplevel-1) :: dw, pa, dt250
    real, intent(in), dimension(nxp,ncoef) ::  coef1, coef2, coef3
    real, intent(inout) :: x1, x2, x3, trant(nplevel)
    
    integer :: iw, ip
    real :: effp, efft, xx, we, pe, fw, fp, pra, prb, prc, ax,  &
            ba, bb, t1, ca, cb, t2

    ! Compute effective pressure (x2)effp and temperature (x3)efft following  !
    ! Eqs. (8.28) and (8.29)                                                  !
    x1=x1+dw(k2-1)
    x2=x2+pa(k2-1)*dw(k2-1)
    x3=x3+dt250(k2-1)*dw(k2-1)

    xx=x1
    effp=x2/x1
    efft=x3/x1

    ! normalize we and pe
    we=(log10(xx)-w1)/dwe
    pe=(log10(effp)-p1)/dpe

    ! Restrict the magnitudes of the normalized we and pe.
    we=min(we,real(nh-1))
    pe=min(pe,real(nxp-1))

    ! assign iw and ip and compute the distance of we and pe from iw and ip.
    iw=int(we+1.)
    iw=min(iw,nh-1)
    iw=max(iw, 2)
    fw=we-real(iw-1)

    ip=int(pe+1.)
    ip=min(ip,nxp-1)
    ip=max(ip, 1)
    fp=pe-real(ip-1)

    ! linear interpolation in pressure  
    pra = coef1(ip,iw-1)*(1.-fp)+coef1(ip+1,iw-1)*fp
    prb = coef1(ip,  iw)*(1.-fp)+coef1(ip+1,  iw)*fp
    prc = coef1(ip,iw+1)*(1.-fp)+coef1(ip+1,iw+1)*fp

    ! quadratic interpolation in absorber amount for coef1	
    ax = (-pra*(1.-fw)+prc*(1.+fw)) *fw*0.5 + prb*(1.-fw*fw)

    ! linear interpolation in absorber amount for coef2 and coef3
    ba = coef2(ip,  iw)*(1.-fp)+coef2(ip+1,  iw)*fp
    bb = coef2(ip,iw+1)*(1.-fp)+coef2(ip+1,iw+1)*fp
    t1 = ba*(1.-fw) + bb*fw

    ca = coef3(ip,  iw)*(1.-fp)+coef3(ip+1,  iw)*fp
    cb = coef3(ip,iw+1)*(1.-fp)+coef3(ip+1,iw+1)*fp
    t2 = ca*(1.-fw) + cb*fw

    ! update the total transmittance between levels k1 and k2	
    trant(k2)= (ax + (t1 + (t2*effp)) * efft) *trant(k2)
     
    ! Max/Min transmission values? commented out in original	
    !	    trant(k2)=min(trant(k2),0.9999999)
    !	    trant(k2)=max(trant(k2),0.0000001)			    

  end subroutine tablup

  subroutine h2okdis(ibn,k,ne,h2oexp,conexp,th2o,tcon,trant)

    !-------------------------------------------------------------------------------!
    ! Compute water vapor transmittance between levels k1 and k2 for m soundings,   !
    ! using the k-distribution method.                                              !
    !                                                                               !
    ! Input parameters:                                                             !
    ! k  = current pressure level                                                   !
    ! ne = no. of terms in each band to compute watervapour continuum transmittance !
    ! h2oexp = exponentials for line absorption                                     !
    ! conexp = exponentials for continuum absorption                                !
    ! Output parameters:                                                            !
    ! th2o = transmittance between levels k1 & k2 due to wv line absorption         !
    ! tcon = transmittance between levels k1 & k2 due to wv continuum absorption    !
    ! trant = total transmittance                                                   !
    !-------------------------------------------------------------------------------!
    
    integer, intent(in) :: ibn, k, ne
    real, intent(in) :: conexp(nplevel-1,3), h2oexp(nplevel-1,6) 
    real, intent(inout) :: th2o(6), tcon(3), trant(nplevel)
    real :: trnth2o

    ! tco2 are the six exp factors between levels k1 and k2                         !
    ! tran is the updated total transmittance between levels k1 and k2              !
    ! h2o are 6 exp factors between levels k1 & k2 from h2o line absorption         !
    ! tcon are 3 exp factors between levels k1 & k2 from h2o continuum absorption   !
    ! trnth2o is the total transmittance between levels k1 and k2 due to both line  !
    !        and continuum absorption.                                              !

    ! Compute th2o following Eq. (8.23).  
    th2o(1) = th2o(1)*h2oexp(k,1)
    th2o(2) = th2o(2)*h2oexp(k,2)
    th2o(3) = th2o(3)*h2oexp(k,3)
    th2o(4) = th2o(4)*h2oexp(k,4)
    th2o(5) = th2o(5)*h2oexp(k,5)
    th2o(6) = th2o(6)*h2oexp(k,6)

    if (ne .eq. 0) then
      ! Compute trnh2o following Eq. (8.25). fkw is given in Table 4.
      trnth2o	   =(fkw(1,ibn)*th2o(1)	           &
        		 + fkw(2,ibn)*th2o(2)	   &
        		 + fkw(3,ibn)*th2o(3)	   &
        		 + fkw(4,ibn)*th2o(4)	   &
        		 + fkw(5,ibn)*th2o(5)	   &
        		 + fkw(6,ibn)*th2o(6))

      trant(k+1)=trant(k+1)*trnth2o

    elseif (ne .eq. 1) then
      ! Compute trnh2o following Eqs. (8.25) and (4.27).
      tcon(1)= tcon(1)*conexp(k,1)
      trnth2o	   =(fkw(1,ibn)*th2o(1)	        &
        	      + fkw(2,ibn)*th2o(2)	&
        	      + fkw(3,ibn)*th2o(3)	&
        	      + fkw(4,ibn)*th2o(4)	&
        	      + fkw(5,ibn)*th2o(5)	&
        	      + fkw(6,ibn)*th2o(6))*tcon(1)

      trant(k+1)=trant(k+1)*trnth2o

    else
      ! For band 3. This band is divided into 3 subbands.
      tcon(1)= tcon(1)*conexp(k,1)
      tcon(2)= tcon(2)*conexp(k,2)
      tcon(3)= tcon(3)*conexp(k,3)
      
      ! Compute trnh2o following Eqs. (4.29) and (8.25).		    
      trnth2o	   = (  gkw(1,1)*th2o(1)		  &
        		 + gkw(2,1)*th2o(2)		  &
        		 + gkw(3,1)*th2o(3)		  &
        		 + gkw(4,1)*th2o(4)		  &
        		 + gkw(5,1)*th2o(5)		  &
        		 + gkw(6,1)*th2o(6))*tcon(1)      &
        		 + (  gkw(1,2)*th2o(1)	          &
        		 + gkw(2,2)*th2o(2)		  &
        		 + gkw(3,2)*th2o(3)		  &
        		 + gkw(4,2)*th2o(4)		  &
        		 + gkw(5,2)*th2o(5)		  &
        		 + gkw(6,2)*th2o(6))*tcon(2)      &
        		 + (  gkw(1,3)*th2o(1)	          &
        		 + gkw(2,3)*th2o(2)		  &
        		 + gkw(3,3)*th2o(3)		  &
        		 + gkw(4,3)*th2o(4)		  &
        		 + gkw(5,3)*th2o(5)		  &
        		 + gkw(6,3)*th2o(6))*tcon(3)

      trant(k+1)=trant(k+1)*trnth2o
    endif

  end subroutine h2okdis

  subroutine co2kdis(k,co2exp,tco2,trant)

    !------------------------------------------------------------------------------!
    ! Compute co2 transmittances between levels k1 and k2 for m soundings, using   !
    ! the k-distribution method with linear pressure scaling.                      !
    !                                                                              !
    ! Output variables:                                                            !
    ! tco2  = transmittance between levels k1 and k2 due to co2 absorption for     !
    !         the various values of the absorption coefficient                     !
    ! trant = total transmittance                                                  !
    !------------------------------------------------------------------------------!

    integer, intent(in) :: k
    real, intent(in) :: co2exp(nplevel-1,6,2)
    real, intent(inout) :: tco2(6,2), trant(nplevel)
    real :: xc

    ! tco2 is the 6 exp factors between levels k1 & k2 computed from Eqs.(8.23) &  !
    ! (8.25), (See Eq.(4.30)). k-distribution functions are given in Table 10.     !

    tco2(1,1)=tco2(1,1)*co2exp(k,1,1)
    xc=   0.1395 *tco2(1,1)

    tco2(2,1)=tco2(2,1)*co2exp(k,2,1)
    xc=xc+0.1407 *tco2(2,1)

    tco2(3,1)=tco2(3,1)*co2exp(k,3,1)
    xc=xc+0.1549 *tco2(3,1)

    tco2(4,1)=tco2(4,1)*co2exp(k,4,1)
    xc=xc+0.1357 *tco2(4,1)

    tco2(5,1)=tco2(5,1)*co2exp(k,5,1)
    xc=xc+0.0182 *tco2(5,1)

    tco2(6,1)=tco2(6,1)*co2exp(k,6,1)
    xc=xc+0.0220 *tco2(6,1)

    ! Band-center region 
    tco2(1,2)=tco2(1,2)*co2exp(k,1,2)
    xc=xc+0.0766 *tco2(1,2)

    tco2(2,2)=tco2(2,2)*co2exp(k,2,2)
    xc=xc+0.1372 *tco2(2,2)

    tco2(3,2)=tco2(3,2)*co2exp(k,3,2)
    xc=xc+0.1189 *tco2(3,2)

    tco2(4,2)=tco2(4,2)*co2exp(k,4,2)
    xc=xc+0.0335 *tco2(4,2)

    tco2(5,2)=tco2(5,2)*co2exp(k,5,2)
    xc=xc+0.0169 *tco2(5,2)

    tco2(6,2)=tco2(6,2)*co2exp(k,6,2)
    xc=xc+0.0059 *tco2(6,2)

    trant(k+1)=trant(k+1)*xc

  end subroutine co2kdis

  subroutine n2okdis(ibn,k,n2oexp,tn2o,trant)

    !-----------------------------------------------------------------------------!
    ! Compute n2o transmittances between levels k1 and k2 for m soundings, using  !
    ! the k-distribution method with linear pressure scaling.                     !
    !                                                                             !
    ! Output parameters                                                           !
    ! tn2o  = transmittance between levels k1 and k2 due to n2o absorption for    !
    !           the various values of the absorption coefficient                  !
    ! trant = total transmittance                                                 !
    !-----------------------------------------------------------------------------!

    integer, intent(in) :: ibn, k
    real, intent(in) :: n2oexp(nplevel-1,4)
    real, intent(inout) :: tn2o(4), trant(nplevel)
    real :: xc

    ! tn2o is computed from Eq. (8.23).                                           !
    ! xc is the total n2o transmittance computed from (8.25)                      !
    ! The k-distribution functions are given in Table 5.                          !

    if (ibn .eq. 6) then
      ! band 6
      tn2o(1)=tn2o(1)*n2oexp(k,1)
      xc= 0.940414*tn2o(1)

      tn2o(2)=tn2o(2)*n2oexp(k,2)
      xc=xc + 0.059586*tn2o(2)
    else
      ! band 7
      tn2o(1)=tn2o(1)*n2oexp(k,1)
      xc= 0.561961*tn2o(1)
      
      tn2o(2)=tn2o(2)*n2oexp(k,2)
      xc=xc + 0.138707*tn2o(2)

      tn2o(3)=tn2o(3)*n2oexp(k,3)
      xc=xc + 0.240670*tn2o(3)

      tn2o(4)=tn2o(4)*n2oexp(k,4)
      xc=xc + 0.058662*tn2o(4)
    endif

    trant(k+1)=trant(k+1)*xc

  end subroutine n2okdis
   
  subroutine ch4kdis(ibn,k,ch4exp,tch4,trant)

    !-----------------------------------------------------------------------------!
    ! compute ch4 transmittances between levels k1 and k2 for m soundings, using  !
    ! the k-distribution method with linear pressure scaling.                     !
    !                                                                             !
    ! Output parameters                                                           !
    ! tch4  = transmittance between levels k1 and k2 due to ch4 absorption for    !
    !           the various values of the absorption coefficient                  !
    ! trant = total transmittance                                                 !
    !-----------------------------------------------------------------------------!

    integer, intent(in) :: ibn, k
    real, intent(in) :: ch4exp(nplevel-1,4)
    real, intent(inout) :: tch4(4), trant(nplevel)
    real :: xc

    ! tch4 is computed from Eq. (8.23). xc is total ch4 transmittance computed    !
    ! from (8.25). The k-distribution functions are given in Table 5.             !

    if (ibn .eq. 6) then
      ! band 6   
       tch4(1)=tch4(1)*ch4exp(k,1)
       xc= tch4(1)
    else
       ! band 7 
       tch4(1)=tch4(1)*ch4exp(k,1)
       xc= 0.610650*tch4(1)
       
       tch4(2)=tch4(2)*ch4exp(k,2)
       xc=xc + 0.280212*tch4(2)

       tch4(3)=tch4(3)*ch4exp(k,3)
       xc=xc + 0.107349*tch4(3)

       tch4(4)=tch4(4)*ch4exp(k,4)
       xc=xc + 0.001789*tch4(4)
    endif

    trant(k+1)=trant(k+1)*xc

  end subroutine ch4kdis

  subroutine comkdis(ibn,k,comexp,tcom,trant)

    !-----------------------------------------------------------------------------!
    ! compute co2-minor transmittances between levels k1 and k2 for m soundings,  !
    ! using the k-distribution method with linear pressure scaling.               !
    !                                                                             !
    ! Output parameters                                                           !
    ! tcom  = transmittance between levels k1 and k2 due to co2-minor absorption  !
    !           for the various values of the absorption coefficient              !
    ! trant = total transmittance                                                 !
    !-----------------------------------------------------------------------------!

    integer, intent(in) :: ibn, k
    real, intent(in) :: comexp(nplevel-1,6)
    real, intent(inout) :: tcom(6), trant(nplevel)
    real :: xc

    ! tcom is computed from Eq. (8.23). xc is total co2 transmittance computed  
    ! from (8.25). The k-distribution functions are given in Table 6. 

    if (ibn .eq. 4) then
       ! band 4  
       tcom(1)=tcom(1)*comexp(k,1)
       xc= 0.12159*tcom(1)
       tcom(2)=tcom(2)*comexp(k,2)
       xc=xc + 0.24359*tcom(2)
       tcom(3)=tcom(3)*comexp(k,3)
       xc=xc + 0.24981*tcom(3)
       tcom(4)=tcom(4)*comexp(k,4)
       xc=xc + 0.26427*tcom(4)
       tcom(5)=tcom(5)*comexp(k,5)
       xc=xc + 0.07807*tcom(5)
       tcom(6)=tcom(6)*comexp(k,6)
       xc=xc + 0.04267*tcom(6)
    else
       ! band 5 
       tcom(1)=tcom(1)*comexp(k,1)
       xc= 0.06869*tcom(1)
       tcom(2)=tcom(2)*comexp(k,2)
       xc=xc + 0.14795*tcom(2)
       tcom(3)=tcom(3)*comexp(k,3)
       xc=xc + 0.19512*tcom(3)
       tcom(4)=tcom(4)*comexp(k,4)
       xc=xc + 0.33446*tcom(4)
       tcom(5)=tcom(5)*comexp(k,5)
       xc=xc + 0.17199*tcom(5)
       tcom(6)=tcom(6)*comexp(k,6)
       xc=xc + 0.08179*tcom(6)
    endif

    trant(k+1)=trant(k+1)*xc

  end subroutine comkdis

  subroutine cfckdis(k,cfcexp,tcfc,trant)

    !-----------------------------------------------------------------------------!
    ! compute cfc-(11,12,22) transmittances between levels k1 and k2 for m        !
    ! soundings, using the k-distribution method with linear pressure scaling.    !
    !                                                                             !
    ! Output parameters                                                           !
    ! tcfc  = transmittance between levels k1 and k2 due to cfc absorption        !
    !           for the various values of the absorption coefficient              !
    ! trant = total transmittance                                                 !
    !-----------------------------------------------------------------------------!
      
    integer, intent(in) :: k
    real, intent(in) :: cfcexp(nplevel-1)
    real, intent(inout) :: tcfc, trant(nplevel)
    
    tcfc=tcfc*cfcexp(k)
    trant(k+1)=trant(k+1)*tcfc

  end subroutine cfckdis

  subroutine b10kdis(k,h2oexp,conexp,co2exp,n2oexp,th2o,tcon,tco2,tn2o,trant, h2Oon,  &
                     co2on, n2Oon)
 
    !-----------------------------------------------------------------------------!
    ! compute h20(line & continuum),co2,n2o transmittances between levels k1 & k2 !
    ! for m soundings, using k-distribution method with linear pressure scaling.  !
    !                                                                             !
    ! Output parameters                                                           !
    ! t(h2o) etc = transmittance between levels k1 and k2 due to (h2o) absorption !
    !           for the various values of the absorption coefficient              !
    ! trant = total transmittance                                                 !
    !-----------------------------------------------------------------------------!

    integer, intent(in) :: k
    real, intent(in) :: h2oexp(nplevel-1,6), conexp(nplevel-1,3), co2exp(nplevel-1,6,2), n2oexp(nplevel-1,4)
    real, intent(inout) :: th2o(6), tcon(3), tco2(6,2), tn2o(4), trant(nplevel)    
    logical, intent(in) :: h2Oon, co2on, n2Oon
    real :: xx

    ! For h2o line. The k-distribution functions are given in Table 4.    
    trant(k+1) = 1.

    if (h2Oon) then
      th2o(1)=th2o(1)*h2oexp(k,1)
      xx= 0.3153*th2o(1)

      th2o(2)=th2o(2)*h2oexp(k,2)
      xx=xx + 0.4604*th2o(2)

      th2o(3)=th2o(3)*h2oexp(k,3)
      xx=xx + 0.1326*th2o(3)

      th2o(4)=th2o(4)*h2oexp(k,4)
      xx=xx + 0.0798*th2o(4)
 
      th2o(5)=th2o(5)*h2oexp(k,5)
      xx=xx + 0.0119*th2o(5)

      trant(k+1)=trant(k+1) * xx

    ! For h2o continuum. Note that conexp(k,3) is for subband 3a.   
      tcon(1)=tcon(1)*conexp(k,3)
      trant(k+1)=trant(k+1)*tcon(1)
    endif
 
    ! For co2 (Table 6)      
    if (co2on) then
      tco2(1,1)=tco2(1,1)*co2exp(k,1,1)
      xx= 0.2673*tco2(1,1)
      
      tco2(2,1)=tco2(2,1)*co2exp(k,2,1)
      xx=xx + 0.2201*tco2(2,1)
      
      tco2(3,1)=tco2(3,1)*co2exp(k,3,1)
      xx=xx + 0.2106*tco2(3,1)

      tco2(4,1)=tco2(4,1)*co2exp(k,4,1)
      xx=xx + 0.2409*tco2(4,1)

      tco2(5,1)=tco2(5,1)*co2exp(k,5,1)
      xx=xx + 0.0196*tco2(5,1)
      
      tco2(6,1)=tco2(6,1)*co2exp(k,6,1)
      xx=xx + 0.0415*tco2(6,1)

      trant(k+1)=trant(k+1)*xx
    endif

    ! For n2o (Table 5)   
    if (n2Oon) then
      tn2o(1)=tn2o(1)*n2oexp(k,1)
      xx= 0.970831*tn2o(1)

      tn2o(2)=tn2o(2)*n2oexp(k,2)
      xx=xx + 0.029169*tn2o(2)
      trant(k+1)=trant(k+1)*(xx-1.)
    endif
       
  end subroutine b10kdis
  
  subroutine interpolate_on_radiation_mesh (old_coordinate, new_coordinate, old_fields, new_fields, input, output)
    
    type(vector_field), intent(in) :: old_coordinate, new_coordinate
    type(scalar_field), dimension(:), intent(in)    :: old_fields
    type(scalar_field), dimension(:), intent(inout) :: new_fields
    logical, optional, intent(in) :: input, output
    
    type(state_type), dimension(1) :: old_state, new_state
    type(vector_field) :: old_position, new_position
    integer, dimension(node_count(new_fields(1))) :: map
    integer :: i, stat, node
    
    call allocate(old_position,old_coordinate%dim,old_fields(1)%mesh,'Coordinate')
    call allocate(new_position,old_coordinate%dim,new_fields(1)%mesh,'Coordinate')

    if (present_and_true(input)) then
      call remap_field(old_coordinate,old_position)
      call set(new_position,new_coordinate)
    else if (present_and_true(output)) then
      call set(old_position,old_coordinate)
      call remap_field(new_coordinate,new_position)
    endif
    new_position%name='Coordinate'
    old_position%name='Coordinate'
    
    call insert(old_state(1), old_position%mesh, name=old_position%mesh%name)
    call insert(new_state(1), new_position%mesh, name=new_position%mesh%name)
    
    call insert(old_state(1), old_position, name="Coordinate")
    call insert(new_state(1), new_position, name="Coordinate")

    ! Insert fields and meshes in new and old states
    do i = 1, size(new_fields)
      call zero(new_fields(i))
      call insert(old_state(1), old_fields(i), name=old_fields(i)%name)
      call insert(new_state(1), new_fields(i), name=old_fields(i)%name)
    enddo
    
    ! Get node ownership
    map = get_element_mapping(old_position, new_position, different_domains=.true.)

    ! Perform interpolations
    call linear_interpolate_states(old_state, new_state, map = map)
!    call interpolate(old_state,new_state,map=map)    

    call deallocate(old_position)
    call deallocate(new_position)
    call deallocate(old_state)
    call deallocate(new_state)

  end subroutine interpolate_on_radiation_mesh
  
  subroutine mesh_to_column (mesh,scalar,column)
    
    type(mesh_type), intent(in) :: mesh
    type(scalar_field), intent(in) :: scalar
    real, dimension(:,:), intent(inout) :: column
    
    integer :: icolumn, ilevel, node
    
    icolumn=0
    ilevel=0
    do node=1,node_count(scalar)
      if (mesh%columns(node) /= icolumn) ilevel=0
      icolumn=mesh%columns(node)
      ilevel=ilevel+1      
      column(icolumn,ilevel)=node_val(scalar,node)
    enddo    
    
    end subroutine mesh_to_column
  
    subroutine column_to_mesh (mesh,column,scalar)
    
    type(mesh_type), intent(in) :: mesh
    real, dimension(:,:), intent(in) :: column
    type(scalar_field), intent(inout) :: scalar
    
    real :: scalar_node
    integer :: icolumn, ilevel, node
    
    icolumn=0
    ilevel=0
    do node=1,node_count(scalar)
      if (mesh%columns(node) /= icolumn) ilevel=0
      icolumn=mesh%columns(node)
      ilevel=ilevel+1      
      scalar_node=column(icolumn,ilevel)
      call set(scalar, node, scalar_node)
    enddo    
    
  end subroutine column_to_mesh
  
end module Radiation_interface  
