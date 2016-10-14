!-*- F90 -*- so emacs thinks this is an f90 file
module sw_radiation_atham

  !--------------------------------------------------------------------!
  ! author:   Michael Herzog  (short wave)                             !
  !                                                                    !
  ! modified by : Rachel White                                         !
  ! date:         Jan 2008                                             !
  !--------------------------------------------------------------------!
  use precision,    only: kint, kreal
  implicit none
  private

  public :: shortwave_radiation_atham, sw_radiation_atham_init
  public :: flxdn_sw_gr_tot, flxdn_sw_gr_dir

 !--------------------------------------------------------------------!
  ! radiation input data                                               !
  !   nsw    [1]     number of shortwave intervals                     !
  !   nlw    [1]     number of longwave intervals                      !
  !   co2    [kg/kg] specific concentration of co2                     !
  !   oxy2   [kg/kg] specific concentation of o2                       !
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
  INTEGER(kint), SAVE :: nsw, nv
  real(kreal), save :: fsol, co2, oxy2, pss, hscale, dtaumax, dtaumin, &
                       ompdiff, expmax, sf, al0
  real(kreal), save, allocatable, dimension(:) :: wsw, solfrc, tauray, &
       abswo2, abswco2, abswo3, abswvap, asl, bsl, csl, dsl, esl, fsl, &
       wavnum, delwav, wlw, bl0, bl1, bl2, cl0, cl1, cl2
  real(kreal), save, allocatable, dimension(:,:) :: flxdn_sw_gr_tot,   &
       flxdn_sw_gr_dir, flxdn_lw_gr_tot
   real(kreal), dimension(:), allocatable, save  :: fdiff, fdir
  integer(kint), save, dimension(6)              :: start_time
  integer(kint), save, allocatable, dimension(:,:) :: surface_type
  LOGICAL, SAVE                                  :: sw_radiation_input,&
                                                    swrad_old   

contains

  !--------------------------------------------------------------------------

  subroutine sw_radiation_atham_init(nx,ny,nz)                   

    !------------------------------------------------------------------!
    ! Allocates arrays                                                 !
    ! Reads in data from INPUT_radiation and INPUT_profile_ozone       !
    !------------------------------------------------------------------!
    use atham_module,   only: cylindric_geometry
    use atham_module,   only: myrank
    use phys_constants, only: r0
    USE atham_module,   ONLY: atham_stop, ztotal, zstart

    integer(kint), intent(in) :: nx, ny,nz                             

    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(kreal) :: tal0  
    real(kreal), dimension(99) :: twsw,tsolfrc, ttauray,          &
         tabswo2,tabswco2, tabswo3, tabswvap
    real(kreal), dimension(99) :: tasl, tbsl, tcsl, tdsl, tesl, tfsl
    real(kreal), dimension(99) :: tbl0, tbl1, tbl2, tcl0, tcl1, tcl2
    REAL(kreal), DIMENSION(240,3,18) :: rad_down_sw 
    INTEGER(kint)                    :: i, a, b, c, level_in

    namelist /radiation_setup_sw/ nsw, dtaumax,dtaumin,ompdiff,expmax,        &
                             pss,hscale,tsolfrc,fsol,                         &
                             twsw, tabswvap, tabswo3, tabswco2,  tabswo2,     & 
                             ttauray, tasl, tbsl, tcsl, tdsl, tesl, tfsl,     &
                             tal0, tbl0, tbl1, tbl2, tcl0, tcl1, tcl2,        &
                             sw_radiation_input, swrad_old

    ! set defaults
    sw_radiation_input =.false.
    swrad_old =.false.

    !------------------------------------------------------------------!
    ! read namelist file INPUT_radiation                               !
    !------------------------------------------------------------------!
    open(56,file='input/INPUT_radiation_sw_atham',form='formatted',status='old',err=2999)
    read(56,radiation_setup_sw,end=2999)

    !------------------------------------------------------------------!
    ! namelist variables cannot have allocatable type                  !
    !------------------------------------------------------------------!
    allocate(wsw(nsw), solfrc(nsw), tauray(nsw), abswo2(nsw),           &
             abswco2(nsw), abswo3(nsw), abswvap(nsw))
    allocate(asl(nsw), bsl(nsw), csl(nsw), dsl(nsw), esl(nsw), fsl(nsw))
    allocate(surface_type(nx,ny))
    allocate(flxdn_sw_gr_tot(nx,ny), flxdn_sw_gr_dir(nx,ny),            &
             flxdn_lw_gr_tot(nx,ny))
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
2999 if (myrank==0) print 2200, 'NO namelist from INPUT_rad_sw read, NO DEFAULTS!!!'
2200 format ('***********************************************',/,a,/,    &
            '***********************************************')
2998 close(56) 


    !------------------------------------------------------------------!
    ! dertermine radiation input per sw band at model top              !
    !------------------------------------------------------------------!
    allocate (fdiff(nsw), fdir(nsw))
    fdiff(:) = r0
    fdir(:)  = r0

    if(sw_radiation_input) then

       level_in = NINT((ztotal+zstart-125._kreal)/250._kreal)

       open(56,file='input/SWRIN',form='unformatted')

       READ(56,END=2899) (((rad_down_sw(a,b,c),a=1,240),b=1,3),c=1,18)

          fdiff(:) = rad_down_sw(level_in,2,:)
          fdir(:) = rad_down_sw(level_in,3,:)

       goto 2898

2899   if (myrank ==0) print 2300, 'NO Radiation Input at model top, END ATHAM!'
2300   format  ('***********************************************',/,a,/,    &
            '***********************************************')
       call atham_stop(' ') 

2898   close(56)
    end if

  end subroutine sw_radiation_atham_init

  !--------------------------------------------------------------------------

  subroutine shortwave_radiation_atham                                     &
                                 (nx, ny, nz, nyv, nyh, nxl, nxr, zv, xv,   &
                                  icenter, density, cptot, &
                                  tetflx, cptgas, rgasnew, tracnew, ntrac,  &
                                  tgasnew, ntgas, p0, pnew, iflgs, ifeld,   &
                                  landseamask,coszen,                       &
                                  watcnew, wetnew, icenew,                  &
                                  oxy2, co2, ozone, dt,                     &
                                  surface_model, use_dynamic_osa, al_ATHAM, &
                                  flxdn_sw_gr_dir, flxdn_sw_gr_tot, totalK)
 
    !=======================================================================!
    !                                                                       !
    ! solution of the radiation transfer equation                           !
    ! in delta-eddington approxomation for shortwave                        !
    !                                                                       !
    ! calculation of heatingrates                                           !
    !                                                                       !
    !=======================================================================!

    use phys_constants, only: r0, r1, gasmin, epsmach, ps0, cpair
    
    integer(kint), intent(in) :: nx, ny, nz, nyv, nyh, nxl, nxr,          &
                                 icenter, ntrac, ntgas
    integer(kint), dimension(:,:), intent(in)   ::  ifeld, landseamask
    integer(kint), dimension(:,:,:), intent(in) :: iflgs 

    real(kreal), intent(in)   :: dt, oxy2, co2, coszen
    REAL(kreal), DIMENSION(:,:,:), INTENT(in) :: density, cptot, &
                                 rgasnew, pnew, wetnew,        &
                                 watcnew, icenew
    real(kreal), dimension(:,:,:,:), intent(in)  :: tracnew, tgasnew
    real(kreal), dimension(:), intent(in) :: ozone, xv, zv, cptgas, p0
    real(kreal), dimension(:,:,:), intent(inout) :: tetflx 
    real(kreal), dimension(:), intent(inout)     :: totalK
    real(kreal), dimension(:,:), intent(inout)   :: flxdn_sw_gr_tot,        &   
                                                    flxdn_sw_gr_dir, al_ATHAM
    logical, intent(in)                          :: surface_model, use_dynamic_osa
    !-----------------------------------------------------------------------!
    ! Local variables                                                       !
    !-----------------------------------------------------------------------!

    real(kreal), dimension(nx, nz) :: flxup, flxdn, dtaup, omp, gp
    REAL(kreal) :: dzv, dxv, tempflx, gasnew, drynew, totgas, totinv, cp, pcom
    real(kreal), dimension(nx) :: albedo

    integer(kint) :: i, j, k, l, ku

    !-------------------------------------------------------------------------!
    ! Set surface_type                                                        !
    ! Currently set to landseamask                                            !
    !   0 = ocean                                                             !
    !   1 = Tall/medium grassland, evergreen shrubland as defined             !
    !       Briegleb 1992                                                     !
    ! Possible to add further definitions!                                    !
    !-------------------------------------------------------------------------!
    surface_type = landseamask

    nv = 2*nz -1

    do j=nyv,nyh
       !--------------------------------------------------------------------!
       ! reset up/down fluxes for xz-slice                                  !
       ! reset direct ground irradiance for x-line                          !
       !--------------------------------------------------------------------!
       do k=1,nz
          do i=nxl,nxr
             flxup(i,k)=r0
             flxdn(i,k)=r0
             flxdn_sw_gr_dir(i,j)=r0
          enddo
       enddo
       !--------------------------------------------------------------------!
       ! perform short wave calculations for nsw wave length intervals      !
       !--------------------------------------------------------------------!
       do l=1,nsw
          call opticsw(j,l,dtaup,omp,gp, density(:,:,:), pnew(:,:,:), p0(:)&
                       , albedo, al_ATHAM(:,:), surface_model, oxy2, co2,   &
                       ozone, watcnew(:,:,:), wetnew(:,:,:), icenew(:,:,:), &
                       nx, ny, nz, nxl, nxr, zv, coszen, use_dynamic_osa)
          call fluxsw(j,l,dtaup,omp,gp,flxup,flxdn, albedo, density(:,:,:), &
                      pnew(:,:,:), p0(:), watcnew(:,:,:), wetnew(:,:,:), &
                      iflgs(:,:,:), ifeld(:,:), nx, ny, nz, nxl, nxr,    &
                      flxdn_sw_gr_dir, coszen)
 

     enddo
       !--------------------------------------------------------------------!
       ! compute heating rates                                              !
       !--------------------------------------------------------------------!
       do k=2,nz
          ku=k-1
          dzv=zv(k)-zv(ku)
          do i=nxl,nxr
             dxv = xv(i) - xv(i-1)

             gasnew=r1-sum(tracnew(i,j,k,1:ntrac))
             gasnew=max(gasmin,min(gasnew,r1))
             drynew=gasnew-sum(tgasnew(i,j,k,1:ntgas))
             drynew=max(epsmach,min(drynew,r1))
             
             totgas=drynew+sum(tgasnew(i,j,k,1:ntgas))
             totinv=r1/totgas

             cp=(cpair*drynew+sum(cptgas(1:ntgas)*tgasnew(i,j,k,1:ntgas)))*totinv

             pcom=((p0(k)+pnew(i,j,k))/ps0)**(rgasnew(i,j,k)/cp)	  

             tempflx=(flxdn(i,k)-flxup(i,k) -flxdn(i,ku)+flxup(i,ku))       &
                         /(density(i,j,k)*cptot(i,j,k)*dzv) *iflgs(i,j,k)
       

             tetflx(i,j,k)=tetflx(i,j,k) + dt*tempflx/pcom*iflgs(i,j,k)

             if (i .EQ. icenter) totalK(k) = tempflx*86169.0_kreal
             
             !--------------------------------------------------------------!
             ! extract and assign total ground irradiance                   !
             ! direct ground irradiance calculated in fluxsw (not exported) !
             !--------------------------------------------------------------!             
             if (k==ifeld(i,j)) flxdn_sw_gr_tot(i,j) = flxdn(i,k)

          enddo
       enddo
    enddo
      
  end subroutine shortwave_radiation_atham
!============================================================================!
  subroutine opticsw(j,l,dtaup,omp,gp, density, pnew, p0, albedo,        &
                     al_ATHAM, surface_model,                            &
                     oxy2, co2, ozone, watcnew, wetnew, icenew,          &
                     nx, ny, nz, nxl, nxr, zv, coszen, use_dynamic_osa)
    !=======================================================================!
    !    calculate optical properties:                                      !
    !    for wave length l in the short wave region                         !
    !                                                                       !
    !    dtaup -- (scaled) optical depth                                    !
    !    omp   -- (scaled) single scattering albedo                         !
    !    gp    -- (scaled) asymmetry factor                                 !
    !=======================================================================!

    use phys_constants,only : r0, r1, epsilon
    
    real(kreal), dimension(nx, nz), intent(inout) :: dtaup, omp, gp 
    real(kreal), dimension(nx), intent(inout)     :: albedo
    real(kreal), intent(in)                       :: oxy2, co2, coszen
    real(kreal), dimension(:), intent(in)         :: ozone, p0
    real(kreal), dimension(:,:,:), intent(in)     :: density, pnew, watcnew,  &
                                                     wetnew, icenew
    real(kreal), dimension(:,:), intent(in)       :: al_ATHAM
    real(kreal), dimension(nz)                    :: zv, xv
    LOGICAL, INTENT(in)                           :: surface_model, use_dynamic_osa
    
    integer(kint), intent(in) :: j, l, nx, ny, nz, nxl, nxr
    
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
    real(kreal) :: scaaer=r0, extaer=r0, gaer  =r0
    real(kreal) :: radwatc, clcover, radiusw, extwat, omwat, gwat, scawat
    real(kreal) :: radice, radiusi, extice, omice, gice, scaice
    real(kreal) :: pssinv, rayleigh, pfac, scamie, g, f, absgas, extinc
    real(kreal) :: dtau, om, omdiff
    REAL(kreal) :: dz, co, oxy, fac

    integer(kint) :: i, k
    
    radwatc=5.32e-6_kreal
    radice = 50e-6_kreal
    clcover=r1
    radiusw=radwatc*1.e6_kreal
    radiusi=radice*1.e6_kreal
    
    extwat=1000._kreal*(asl(l)+bsl(l)/radiusw)*clcover*sqrt(clcover)
    omwat =r1-csl(l)-dsl(l)*radiusw
    gwat  =   esl(l)+fsl(l)*radiusw
    scawat=extwat*omwat

    extice= 1000._kreal*(al0/radiusi)*clcover*sqrt(clcover)
    omice = r1 - bl0(l) - bl1(l)*radiusi - bl2(l)*radiusi*radiusi
    gice  = cl0(l) + cl1(l) * radiusi + cl2(l) *radiusi*radiusi
    scaice= extice*omice

    !-----------------------------------------------------------------!
    ! start loop over xz-slice                                        !
    !-----------------------------------------------------------------!
    
    !marker 
   
    IF (swrad_old) THEN
       co = 4.558e-4_kreal
       oxy=.2314_kreal
       fac = 50._kreal
    else
       co = co2
       oxy= oxy2
       fac= r1
    END IF

    pssinv=r1/pss
    
    do k=2,nz
       dz=zv(k)-zv(k-1)
       do i=nxl,nxr
          pfac=(p0(k)+pnew(i,j,k))*pssinv
          scamie=scawat*watcnew(i,j,k) + scaaer + scaice*icenew(i,j,k)
          rayleigh=tauray(l)*pfac/(hscale*density(i,j,k))
          g=(gwat*scawat*watcnew(i,j,k) + gice*scaice*icenew(i,j,k) + gaer*scaaer) &
               /max(epsilon,scamie+rayleigh) 
          f=g*g
          
          absgas= abswvap(l) *                                             &
                 (wetnew(i,j,k)*pfac + 0.0008*sqrt(coszen*wetnew(i,j,k)))+ &
                 abswco2(l)*sqrt(coszen*co)                             + &
                 abswo2 (l)*sqrt(coszen*oxy)                            + &
                 abswo3 (l)* fac *ozone(k)     

          extinc=extwat*watcnew(i,j,k) + extaer + rayleigh + absgas + extice*icenew(i,j,k)
          dtau=extinc*density(i,j,k)*dz

          !------------------------------------------------------------------!
          ! cut in dtau                                                      !
          !------------------------------------------------------------------!
          
          dtau=max(dtaumin,dtau)
          om=(scamie+rayleigh)/max(epsilon,extinc)
          dtaup(i,k)=(r1-om*f)*dtau
          omp(i,k)  =(r1-f)*om/(r1-om*f)

          !------------------------------------------------------------------!
          ! cut in omp                                                       !
          !------------------------------------------------------------------!

          omp(i,k)=min(omp(i,k),r1-ompdiff)
          gp(i,k) =(g-f)/(r1-f)
          
       enddo
    enddo

    !----------------------------------------------------------------!
    ! Calculate albedo                                               !
    ! sf is solar angle factor (now a public variable)               !
    !----------------------------------------------------------------!
    
    sf = 1.4_kreal / (1.0_kreal + 0.4_kreal * coszen)

    do i = nxl, nxr
       if (surface_model) then
          if (surface_type(i,j) .eq. 0) then
             if (use_dynamic_osa) then
                albedo(i) = min(al_ATHAM(i,j),r1)     ! sf included in dynamic_osa
             else
                albedo(i) = min(al_ATHAM(i,j) * sf,r1)
             endif

          elseif (surface_type(i,j) .eq. 1) then
             albedo(i) = min(al_ATHAM(i,j) * sf,r1)
          endif
       else
          if (surface_type(i,j) .eq. 0) then
             albedo(i) = 0.06_kreal * sf
          else if (surface_type(i,j) .eq. 1) then
             albedo(i) = 0.19_kreal * sf
          else 
             albedo(i) = 0.0_kreal
             write(*,*) 'surface type undefined'       
          endif
       endif
    enddo

  end subroutine opticsw

!=============================================================================¬

  subroutine fluxsw(j,l,dtaup,omp,gp,flxup,flxdn, albedo, density, pnew, p0, &
                    watcnew, wetnew, iflgs, ifeld, nx, ny, nz, nxl, nxr,     &
                    flxdn_sw_gr_dir, coszen)

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
    !          fdinc  -- incident diffuse down-flux at level nz (top)       !
    !          fuinc  -- incident diffuse up-flux at level 1 (bottom)       !
    !                                                                       !
    ! output:  flxup, flxdn -- up/down fluxes in w/m^2                      ! 
    !-----------------------------------------------------------------------!
    use phys_constants, only: r0,r1, r2, r1h, r3q, epsilon
    use atham_module,  only : atham_stop 
    
    real(kreal), dimension(:,:), intent(inout) :: flxdn_sw_gr_dir
    real(kreal), dimension(nx, nz), intent(in) :: omp, gp, density, pnew, &
                                                  p0, watcnew, wetnew 
    real(kreal), dimension(nx), intent(in) :: albedo
    real(kreal), dimension(nx, nz), intent(inout) :: dtaup, flxup, flxdn
    real(kreal), intent(in)               :: coszen
    integer(kint), intent(in) :: j, l, nx, ny, nz, nxl, nxr
    integer(kint), intent(in) :: ifeld(:,:), iflgs(:,:,:) 

    !-----------------------------------------------------------------------!
    ! local variables                                                       !
    !-----------------------------------------------------------------------!
    integer(kint) ::  i, k, km1, k2, k2p1, k2m1, klow
    real(kreal) :: fdinc, fuinc, fsun, dir
    real(kreal) :: dblflag
    real(kreal) :: h1, h2, h3, xk, xk2, expon, smu 
    real(kreal) :: pp, pm
    real(kreal), dimension(nz) :: p, ex, extau
    real(kreal), dimension(nv,6) :: pen
    real(kreal), dimension(nv) :: solut
    real(kreal), dimension(nz) :: alph, beta
    character(len=50) :: message

    real(kreal), parameter :: r3=3._kreal
 
    if (sw_radiation_input) then 
       fdinc = fdiff(l)*coszen
       fsun  = fdir(l)
       fuinc = r0
    else
       fdinc=r0
       fuinc=r0
       fsun =wsw(l)*solfrc(l)*fsol
    endif
    
    message = 'nv is not equal to 2*nz -1 further down module'

    if (nv .NE. 2*nz - 1) call atham_stop(message)
    
    !-----------------------------------------------------------------------!
    ! precalculate coefficients                                             !
    !-----------------------------------------------------------------------!

    do i=nxl,nxr
       do k=2,nz
          dblflag=real(iflgs(i,j,k),kreal)
          
          h1 =r1-omp(i,k)
          h2 =r1-omp(i,k)*gp(i,k)
          xk2=r3*h1*h2
          h3 =r1/(coszen*coszen)-xk2

    !-----------------------------------------------------------------------!
    ! cut in h3                                                             !
    !-----------------------------------------------------------------------!
          h3 =sign(r1,h3)*max(epsilon,abs(h3))
          xk =sqrt(xk2)
          p(k)=r2*xk/(r3*h2)*dblflag
          alph(k)=r3q*fsun*omp(i,k)*(r1+gp(i,k)*h1)/h3*dblflag
          beta(k)=r1h*fsun*omp(i,k)*dblflag*                  &
               (r1/coszen+r3*coszen*gp(i,k)*h1)/h3

     !----------------------------------------------------------------------!
     ! cut in dtaup                                                         !
     !----------------------------------------------------------------------!
          dtaup(i,k)=min(xk*dtaup(i,k),dtaumax)/xk
          expon=xk*dtaup(i,k)
          ex(k)=exp(expon)
       enddo

       extau(nz)=dtaup(i,nz)
       do k=nz-1,2,-1
          extau(k)=extau(k+1)+dtaup(i,k)
       enddo
       smu=coszen
       do k=2,nz
          dblflag=real(iflgs(i,j,k),kreal)
          expon=min(expmax,extau(k)/smu)
          extau(k)=dblflag/exp(expon)
       enddo

       !--------------------------------------------------------------------!
       ! calculate matrix: upper boundary                                   !
       !--------------------------------------------------------------------!

       pen(nv,1) = r0
       pen(nv,2) = r0
       pen(nv,3) = r1
       pen(nv,4) = r0
       pen(nv,5) = r0
       pen(nv,6) = r0

       pen(nv-1,5)=r0
       pen(nv-1,4)=r0
       pen(nv-1,3)=(r1+p(nz))*ex(nz)
       pen(nv-1,2)=(r1-p(nz))/ex(nz)
       pen(nv-1,1)=r0
       pen(nv-1,6)=alph(nz)+beta(nz) +fdinc

       !--------------------------------------------------------------------!
       ! atmosphere                                                         !
       !--------------------------------------------------------------------!

       klow=ifeld(i,j)
       do k=klow+1,nz
          km1=k-1
          k2=2*k-4
          k2p1=k2+1
          pen(k2  ,5)= p(k)
          pen(k2  ,4)=-p(k)
          pen(k2  ,3)=-p(km1)*ex(km1)
          pen(k2  ,2)= p(km1)/ex(km1)
          pen(k2  ,1)= r0
          pen(k2  ,6)= (beta(k)-beta(km1))*extau(k)
                      
          pen(k2p1,5)= r0
          pen(k2p1,4)=pen(k2,5)-r1
          pen(k2p1,3)=pen(k2,4)-r1
          pen(k2p1,2)=pen(k2,3)+   ex(km1)
          pen(k2p1,1)=pen(k2,2)+r1/ex(km1)
          pen(k2p1,6)=pen(k2,6)-(alph(k)-alph(km1))*extau(k)
       enddo

       !---------------------------------------------------------------------!
       ! lower boundary                                                      !
       !---------------------------------------------------------------------!

       k=klow
       k2=2*k-4
       k2p1=k2+1
       pp=r1+p(k)
       pm=r1-p(k)
       pen(k2p1,5)=r0
       pen(k2p1,4)=pm-pp*albedo(i)
       pen(k2p1,3)=pp-pm*albedo(i)
       pen(k2p1,2)=r0
       pen(k2p1,1)=r0
       pen(k2p1,6)=fuinc                                                             &
                  +(alph(k)-beta(k)+albedo(i)*(fsun*coszen-alph(k)-beta(k)))  &
                  *extau(k)

       !----------------------------------------------------------------------!
       ! ground                                                               !
       !----------------------------------------------------------------------!

       do k=3,klow
          k2=2*k-4
          k2m1=k2-1
          pen(k2  ,5)=r0
          pen(k2  ,4)=r0
          pen(k2  ,3)=r1
          pen(k2  ,2)=r0
          pen(k2  ,1)=r0
          pen(k2  ,6)=r0
                      
          pen(k2m1,5)=r0
          pen(k2m1,4)=r0
          pen(k2m1,3)=r1
          pen(k2m1,2)=r0
          pen(k2m1,1)=r0
          pen(k2m1,6)=r0
       enddo
     
       !-----------------------------------------------------------------------!
       ! solve matrix system                                                   !
       !-----------------------------------------------------------------------!

       call pendia(pen,solut)

       !-----------------------------------------------------------------------!
       ! calculate fluxes                                                      !
       !-----------------------------------------------------------------------!

       dir=fsun*coszen
       flxdn(i,nz)=flxdn(i,nz)+fdinc+dir
       flxup(i,nz)=flxup(i,nz)                          &
                    +(r1+p(nz))/ex(nz)*solut(nv-2)      &
                    +(r1-p(nz))*ex(nz)*solut(nv-1)      &
                    -alph(nz)+beta(nz)

       do k=2,nz
          k2=2*k-3
          k2p1=k2+1
          dir=fsun*coszen*extau(k)
          flxdn(i,k-1)= flxdn(i,k-1)                    &
                        +(r1-p(k))*solut(k2)            &
                        +(r1+p(k))*solut(k2p1)          &
                        -(alph(k)+beta(k))*extau(k)     &
                        +dir
          flxup(i,k-1)= flxup(i,k-1)                    &
                        +(r1+p(k))*solut(k2)            &
                        +(r1-p(k))*solut(k2p1)          &
                        -(alph(k)-beta(k))*extau(k)     
       enddo

       !-----------------------------------------------------------------------!
       ! calculate ground irradiance (direct)                                  !
       !-----------------------------------------------------------------------!
                
       flxdn_sw_gr_dir(i,j) = flxdn_sw_gr_dir(i,j) + fsun*coszen*extau(ifeld(i,j))

    enddo

  end subroutine fluxsw
  
!==============================================================================!

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
      
    use atham_module,   only : nz
    use phys_constants, only : r1
    
    integer(kint) :: l
    real(kreal), dimension(nv) :: solut(nv)
    real(kreal), dimension(nv,6) :: pen
    real(kreal) :: fac1, fac2
    
    !-----------------------------------------------------------------------!
    ! normalization of the main diagonal                                    !
    !-----------------------------------------------------------------------!
    
    do l=1,nv
       fac1=r1/pen(l,3)
       pen(l,1)=pen(l,1)*fac1
       pen(l,2)=pen(l,2)*fac1
       pen(l,3)=pen(l,3)*fac1
       pen(l,4)=pen(l,4)*fac1
       pen(l,5)=pen(l,5)*fac1
       pen(l,6)=pen(l,6)*fac1
    enddo
    
    !------------------------------------------------------------------------!
    ! triangularization                                                      !
    !------------------------------------------------------------------------!
    
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
    
    !------------------------------------------------------------------------!
    ! calculation of the unknowns                                            !
    !------------------------------------------------------------------------!

    solut(nv  ) =  pen(nv  ,6)/pen(nv  ,3)
    solut(nv-1) = (pen(nv-1,6)-pen(nv-1,4)*solut(nv)) /pen(nv-1,3)
    
    do l=nv-2,1,-1
       solut(l)=(pen(l,6)                       &
                   -pen(l,4)*solut(l+1)         &
                   -pen(l,5)*solut(l+2))        &
                   /pen(l,3)
    enddo

  end subroutine pendia

end module sw_radiation_atham
