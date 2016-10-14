!-*- F90 -*- so emacs thinks this is an f90 file

module lw_radiation_nasa

  !--------------------------------------------------------------------!
  !           NASA ...        (long wave)                              !
  !                                                                    !
  ! modified by : Rachel White                                         !
  ! date:         Jan 2008                                             !
  !--------------------------------------------------------------------!
  use precision,    only: kint, kreal

  implicit none
  private

  public :: longwave_radiation_nasa, lw_radiation_nasa_init

  !--------------------------------------------------------------------!
  ! radiation input LW                                                 !
  !   wavenum[1/m]   wavenumber of longwave in 1/micrometre            !
  !   delwav [m]     wave length interval in micrometre                !
  !   wlw            weights of longwave spectral intervals            !
  !   ablwvap[m2/kg] absorption coefficient of water vapour (longwave) !
  !   ablwo3 [m2/kg] absorption coefficient of o3 (longwave)           !
  !   ablwco2[m2/kg] absorption coefficient of co2 (longwave)          !
  !   ablwch4[m2/kg] absorption coefficient of ch4 (longwave)          !
  !                                                                    !
  !   coszen         cosine of zenith angle                            !
  !   sf     [/]     solar angle factor (influence on albedo)          !
  !   surface_type   surface type for calculating albedo               !
  !   flxdn_*[W/m2]  short wave (sw) and long wave (lw) ground (gr)    !
  !                  irradiance total (tot) and direct (dir)           !
  !--------------------------------------------------------------------!

  real(kreal), save, allocatable, dimension(:) :: xkw, xke, aw, bw, pm
  real(kreal), save, allocatable, dimension(:,:) :: fkw, gkw, aib, awb,    &
                     aiw, aww, aig, awg, cb
  integer(kint), save, allocatable, dimension(:) :: mw
  integer(kint), save :: nplevel, npdiff

  !----------------------------------------------------------------------!
  ! local definitions                                                    !
  !                                                                      !
  ! parameters defining the size of the pre-computed tables for          !
  ! transmittance using table look-up.                                   !
  !                                                                      !
  ! nxp = number of pressure intervals                                   !
  ! no  = number of o3 intervals                                         !
  ! nc  = number of co2 intervals                                        !
  ! nh  = number of h20 intervals                                        !
  !----------------------------------------------------------------------!
  integer(kint), parameter :: nxp=26, no=21, nc=30, nh=31
  !------------------------------------------------------------------------!
  !     include tables used in the table look-up for co2 (band 3),         !
  !     o3 (band 5), and h2o (bands 1, 2, and 7) transmission functions    !
  !     "co2.tran4" is the co2 transmission table applicable to a large    !
  !     range of co2 amount (up to 100 times of the present-time value).   !
  !------------------------------------------------------------------------!
  real(kreal), save, dimension(nxp, nc) :: c1, c2, c3
  real(kreal), save, dimension(nxp, no) :: o1, o2, o3
  real(kreal), save, dimension(nxp, nh) :: h11, h12, h13,    &
                                           h21, h22, h23,    &
                                           h81, h82, h83
  integer(kint) ip,iw
  LOGICAL, save                         :: lwrad_old 

#include "include/h2o.tran3"
#include "include/co2.tran4"
#include "include/o3.tran3"


!===========================================================================

contains 
  subroutine lw_radiation_nasa_init(nx,ny,nz)

    !------------------------------------------------------------------!
    ! Allocates arrays                                                 !
    ! Reads in data from INPUT_radiation and INPUT_profile_ozone       !
    !------------------------------------------------------------------!
    use atham_module,   only: myrank
    use phys_constants, only: r0

    integer(kint), intent(in) :: nx, ny, nz

    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(kreal) :: txkw(9), txke(9), taw(9), tbw(9),        &
         tpm(9), tfkw(6,9), tgkw(6,3)
    real(kreal) :: taib(3,10)
    real(kreal), dimension(4,10) :: tawb, taiw, taww, taig, tawg
    real(kreal), dimension(6,10) :: tcb
    integer(kint) :: tmw(9)

    namelist /radiation_setup_lw/                                             &
                             txkw, txke, tmw, taw, tbw, tpm, tfkw, tgkw,      &
                             taib, tawb, taiw, taww, taig, tawg, tcb,         &
                             npdiff, lwrad_old

    lwrad_old = .false.

    !------------------------------------------------------------------!
    ! read namelist file INPUT_radiation                               !
    !------------------------------------------------------------------!
    open(56,file='input/INPUT_radiation_lw_nasa',form='formatted',status='old',err=2999)
    read(56,radiation_setup_lw,end=2999)

    !------------------------------------------------------------------!
    ! namelist variables cannot have allocatable type                  !
    !------------------------------------------------------------------!

    allocate(xkw(9), xke(9), aw(9), bw(9), pm(9), fkw(6,9), gkw(6,3),   &
              mw(9))
    allocate(aib(3,10), awb(4,10), aiw(4,10), aww(4,10), aig(4,10), awg(4,10))
    allocate(cb(6,10))

   
    xkw     = txkw
    xke     = txke
    mw      = tmw
    aw      = taw
    bw      = tbw
    pm      = tpm
    fkw     = tfkw
    gkw     = tgkw
    mw      = tmw
    aib = taib
    awb = tawb
    aiw = taiw
    aww = taww
    aig = taig
    awg = tawg
    cb  = tcb
   
    goto 2998
2999 if (myrank==0) print 2200, 'NO namelist from INPUT_rad_nasa read, NO DEFAULTS!!!'
2200 format ('***********************************************',/,a,/,    &
            '***********************************************')
2998 close(56) 
    !----------------------------------------------------------------------!
    ! set nplevel
    !----------------------------------------------------------------------!
    nplevel = (nz + npdiff)
   
  end subroutine lw_radiation_nasa_init

!----------------------------------------------------------------------------

  subroutine longwave_radiation_nasa                                       &
                               (nx, ny, nz, nyv, nyh, nxl, nxr, z, zv, xv, &
                                tempnew, p0, pnew, plevel, iflgs, ifeld,   &
                                wetnew, watcnew, icenew,                    &
                                co2, n2o, ch4, o2, cfc11, cfc12,           &
                                cfc22, ccl4, ozone,                        &
                                frcld, icld, wcld, radice, radwat,         &
                                emis, tsfc, surface_model, heatlwrate,     &
                                flxdn_lw_gr_tot)
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

    use phys_constants, only: r0, r1, pi, ps0, grav
    use phys_constants, only: gasmin, epsmach
    use process_data,   only: watpnew, rhowat, granew

    integer(kint), intent(in) ::  nx, ny, nz, nyv, nyh, nxl, nxr
    integer(kint), intent(in) :: ifeld(:,:), iflgs(:,:,:)
    real(kreal), intent(in)   :: co2, n2o, ch4, o2, cfc11, cfc12, cfc22, ccl4, &
                                 radice, radwat
    real(kreal), dimension(nz), intent(in)       :: ozone, p0, z, zv, xv
    real(kreal), dimension(nx,ny), intent(in)    :: emis, tsfc
    real(kreal), dimension(nx,ny,nz), intent(in) :: frcld, icld, wcld,       &
                                                    tempnew, pnew,           &
                                                    wetnew, icenew, watcnew,  &
                                                    plevel
    logical, intent(in)                          :: surface_model

    real(kreal), dimension(nx,ny), intent(inout) :: flxdn_lw_gr_tot
    real(kreal), dimension(nx,ny,nz), intent(out):: heatlwrate

    real(kreal), parameter :: gasmax = r1, r3=3._kreal, r4=4._kreal
    real(kreal), parameter :: watcri = 1.0e-5
    
    !---------------------------------------------------------------------!
    ! temporary arrays                                                    !
    !------------------------------------------------------  -------------!

    real(kreal), dimension(nx, 10) :: eg
    real(kreal), dimension(nx) :: tg, tb
    real(kreal), dimension(nx, nplevel) :: pl, flx
    real(kreal), dimension(nx, nplevel-1) :: ta, wa, oa, fcld
    real(kreal), dimension(nx, nplevel-1, 3) :: cwc, reff
    REAL(kreal), DIMENSION(nx, ny, nz) ::  flxlwtmp, heatlwrate2
    real(kreal) :: xpi
    REAL(kreal) :: avg_pl, deltap, radwat2, radice2
    REAL(kreal) :: co2lw, ch4lw, n2olw, cfc11lw, cfc12lw, cfc22lw

    integer(kint) :: i, ibin, j, k, klow
    
    xpi =r4/r3*pi

   do j=nyv,nyh

       !!!$$$ PUT INTO SETUP??!
       !------------------------- ---------------------------------------!
       ! set ground emissivity, eg, read from process_data if surface    !
       ! module is  coupled and initialized                              !
       !-----------------------------------------------------------------!

       do ibin = 1 , 10
          DO i=nxl, nxr
             eg(i, ibin) = emis(i,j)
          END DO
       enddo
       !-------------------------------------------------------------------!
       ! assign ground temperature, tg                                     !
       ! perturbation : dtempg is a randomly selected number with max      !
       ! amplitude                                                         ! 
       !-------------------------------------------------------------------!

       !#ifdef PERTURBATION
       !       do i=nxl,nxr
       !          tg(i)=tempground1+dtempg(i,j)
       !       enddo
       !#else
       DO i = nxl, nxr
          tg(i) = tsfc(i,j)
       ENDDO

       !#endif

       !------------------------------------------------------------------!
       ! assign air temperature at ground                                 !
       !------------------------------------------------------------------!
       do i=nxl,nxr
          if (surface_model) then
             tb(i) = tempnew(i,j,ifeld(i,j))
          else
             tb(i) = tg(i)
          endif
       enddo
  
       !-------------------------------------------------------------------!
       ! assign variables in pressure levels!                              !
       !-------------------------------------------------------------------!
      
       IF (.NOT. lwrad_old) THEN
         radice2 = radice
         radwat2 = radwat
         co2lw   = co2
         n2olw   = n2o
         ch4lw   = ch4
         cfc11lw = cfc11
         cfc12lw = cfc12
         cfc22lw = cfc22
          DO k=2,nz
             DO i=nxl,nxr
                pl(i,nplevel+2-k) = plevel(i,j,k)           
                ta(i,nplevel+1-k) = tempnew(i,j,k)
                wa(i,nplevel+1-k) = wetnew(i,j,k)
                oa(i,nplevel+1-k)   = ozone(k)

                cwc(i,nplevel+1-k,1) = icld(i,j,k)
                cwc(i,nplevel+1-k,2) = wcld(i,j,k)            !cloud water
                cwc(i,nplevel+1-k,3) = watpnew(i,j,k)         !rain water
                fcld(i,nplevel+1-k)  = frcld(i,j,k)           !cloud fraction
             ENDDO
          ENDDO
       ELSE
         radice2 = 50._kreal
         radwat2 = 10._kreal
         co2lw = 350.e-6
         n2olw=0.28e-6
         ch4lw=1.75e-6
         cfc11lw=0.3e-9
         cfc12lw=0.5e-9
         cfc22lw=0.2e-9
          DO k=2,nz
             DO i=nxl,nxr
                pl(i,nplevel+2-k) = plevel(i,j,k)           
                ta(i,nplevel+1-k) = tempnew(i,j,k)
                wa(i,nplevel+1-k) = wetnew(i,j,k)
                oa(i,nplevel+1-k)   = ozone(k)       
                cwc(i,nplevel+1-k,1) = icenew(i,j,k) +granew(i,j,k) 

                IF (watcnew(i,j,k) .GE. watcri) THEN
                   cwc(i,nplevel+1-k,2) = watcnew(i,j,k)            
                   fcld(i,nplevel+1-k)= r1                          
                ELSE
                   cwc(i,nplevel+1-k,2) = r0 
                   fcld(i,nplevel+1-k)  = r0
                ENDIF
                cwc(i,nplevel+1-k,3) = watpnew(i,j,k)               !rain water

                flx(i,nplevel+2-k)   = r0
             ENDDO
          enddo
       END IF

       k=nz+1
       do i=nxl,nxr
          pl(i,nplevel+2-k) = p0(k-1)+pnew(i,j,k-1)                          &
                                 + (zv(k-1)-z(k-1))/(z(k-1)-zv(k-2))         &      
                                 *(p0(k-1)+pnew(i,j,k-1)-pl(i,nplevel+3-k))
          flx(i,nplevel+2-k) = r0
       enddo

       !------------------------------------------------------------------!
       ! convert pressure level unit from Pa to hPa                       !
       !------------------------------------------------------------------!
       do k=nplevel,npdiff+1,-1
          do i=nxl,nxr
             pl(i,k)=pl(i,k)/100.0_kreal
          enddo
       enddo
       
       avg_pl = (SUM(pl(nxl:nxr,nplevel+1-nz))/(nxr-nxl+1))
       deltap = (avg_pl - 1)/npdiff
       deltap = min(deltap, 4.0_kreal)

       !-------------------------------------------------------------------!
       ! top boundary layer, k from 6 to 1                                 !
       !-------------------------------------------------------------------!
       do k=npdiff,1,-1
          do i=nxl,nxr
             pl(i,k)=pl(i,k+1)-deltap
             ta(i,k)=ta(i,k+1)
             wa(i,k)=wa(i,k+1)
             oa(i,k)=oa(i,k+1)
             cwc(i,k,1)=r0
             cwc(i,k,2)=r0
             cwc(i,k,3)=r0
             fcld(i,k)=r0
             flx(i,k) =r0
          enddo
       enddo
       !-------------------------------------------------------------------!
       ! set effective cloud particle sizes, 1 = ice, 2 = cloud water,     ! 
       !  3 = rain water (micrometres)                                     !
       !-------------------------------------------------------------------!
       do k= 1, nplevel-1
          do i= nxl,nxr
             !reff(i,k,1)  = 80._kreal
             reff(i,k,1) =  radice2
             reff(i,k,2) =  radwat2
             reff(i,k,3) = 1000._kreal
          enddo
       enddo

       !-------------------------------------------------------------------!
       ! set the cloud droplet (cloud water) effective radius              !
       !-------------------------------------------------------------------!
!!$       do k=2,nz
!!$          do i=nxl,nxr
!!$             if(watcnew(i,j,k).gt.watcri .and. tnumnew(i,j,k,1).gt.tnumcri)  then
!!$                reff(i,nplevel+1-k,2)=min(2.5d-5,max(1.0d-6,                     &
!!$                     (watcnew(i,j,k)/tnumnew(i,j,k,1)/rhowat/xpi)**0.33))
!!$                reff(i,nplevel+1-k,2) =1.0d6*reff(i,nplevel+1-k,2)!m to um
!!$                reff(i,nplevel+1-k,2) =1.08*reff(i,nplevel+1-k,2) !rv-->re
!!$             endif
!!$          enddo
!!$       enddo

       !-------------------------------------------------------------------!
       ! calculate infrared radiation flux                                 !
       !-------------------------------------------------------------------!
       call irrad(pl,ta,wa,oa,tb, tg,eg,cwc,reff,fcld,flx,j, nxl, nxr, nx, &
        ny, p0, ifeld, iflgs, co2lw, n2olw, ch4lw, cfc11lw, cfc12lw,       &
        cfc22lw, flxdn_lw_gr_tot)

       !-------------------------------------------------------------------!
       ! calculate infrared cooling rate                                   !
       !-------------------------------------------------------------------!
       do  i=nxl,nxr                     !(*)
          klow = ifeld(i,j)
          do k=npdiff+1,nplevel+1-klow
             heatlwrate(i,j,nplevel+1-k)=(flx(i,k+1)-flx(i,k))         &
                                           *8.441874_kreal/(pl(i,k+1)-pl(i,k))

             flxlwtmp(i,j,nplevel-k)=-flx(i,k+1)

             heatlwrate2(i,j,nplevel+1-k)=(flx(i,k+1)-flx(i,k))/(pl(i,k+1)-pl(i,k))
      
          enddo
       enddo

       k=npdiff
       do i=nxl,nxr
          flxlwtmp(i,j,nplevel-k)=-flx(i,k+1)
       enddo

    enddo
    
         
  end subroutine longwave_radiation_nasa

!=======================================================================

 subroutine irrad (pl,ta,wa,oa,tb, tg,eg,cwc,reff,fcld, flx, j, nxl, nxr, nx, &
        ny, p0, ifeld, iflgs, co2lw, n2o, ch4, cfc11, cfc12, cfc22,           &
        flxdn_lw_gr_tot)

    use phys_constants,only: r1, r0, r1h, r1q

    !-----------------------------------------------------------------------!
    ! input parameters                                                      !
    !-----------------------------------------------------------------------!
    
    real(kreal), intent(in)  :: co2lw, n2o, ch4, cfc11, cfc12, cfc22
    real(kreal), dimension(nx, 10), intent(in)           :: eg
    real(kreal), dimension(nx), intent(in)               :: tg, tb
    real(kreal), dimension(nx, nplevel), intent(in)      :: pl
    real(kreal), dimension(nx, nplevel-1), intent(in)    :: ta, wa, oa, fcld
    real(kreal), dimension(nx, nplevel-1, 3), intent(in) :: cwc, reff
    real(kreal), dimension(:), intent(in)                :: p0
    integer(kint), intent(in)                            :: j, nxl, nxr, nx, &
                                                            ny
    integer, intent(in)                             :: ifeld(:,:), iflgs(:,:,:)
    real(kreal), dimension(nx,ny), intent(inout)         :: flxdn_lw_gr_tot
    
    !----------------------------------------------------------------------!
    ! output parameters                                                    !
    !----------------------------------------------------------------------!
    
    real(kreal), dimension(nx, nplevel), intent(inout) :: flx

    !----------------------------------------------------------------------!
    ! local tracer definitions                                             !
    !----------------------------------------------------------------------!

    real(kreal), parameter :: T250 = 250.0, C789 = 789.0
   
    !----------------------------------------------------------------------!
    ! temporary arrays                                                     !
    !----------------------------------------------------------------------!

    real(kreal), dimension(nx, nplevel-1) :: pa, dt250, dp
    real(kreal), dimension(nx, nplevel-1) :: dh2o,dco2,do3,dn2o,dch4,dcont, &
                                             df11,df12,df22
    real(kreal), dimension(nx, nplevel-1, 3) :: cwp, taucl
    real(kreal), dimension(nx, nplevel) :: transfc
    real(kreal), dimension(nx) :: tx, xlayer, bs, rflxs, flag
    real(kreal), dimension(nx, 0:nplevel) ::blayer
    real(kreal), dimension(nx, nplevel) :: blevel
    real(kreal), dimension(nx, nplevel-1) :: tcldlyr
    real(kreal), dimension(nx, nplevel-1, 6)   :: h2oexp, comexp
    real(kreal), dimension(nx, nplevel-1, 3)   :: conexp
    real(kreal), dimension(nx, nplevel-1, 6,2) :: co2exp
    real(kreal), dimension(nx, nplevel-1, 4)   :: n2oexp, ch4exp
    real(kreal), dimension(nx, nplevel-1)      :: f11exp, f12exp, f22exp
    real(kreal), dimension(nx, nplevel) :: flxu, flxd
    real(kreal), dimension(nx, 0:nplevel) :: bu, bd
    real(kreal), dimension(nx, 6)    :: th2o, tcom
    real(kreal), dimension(nx, 6, 2) :: tco2
    real(kreal), dimension(nx, 4)    :: tn2o, tch4
    real(kreal), dimension(nx)       :: tf11, tf12, tf22
    real(kreal), dimension(nx, 3)    :: tcon
    real(kreal), dimension(nx) :: x1, x2, x3
    real(kreal), dimension(nx, nplevel) :: fclr, trant



    real(kreal) :: xx, tauc,                                                &
                   reff1, reff2, w1, w2, w3, ww, g1, g2, g3, ff, gg,        &
                   a1, b1, fk1, a2, b2, fk2, p1, dwe, dpe,                  &
                   yy

    logical :: oznbnd,co2bnd,h2otbl,conbnd,n2obnd
    logical :: ch4bnd,combnd,f11bnd,f12bnd,f22bnd,b10bnd, h2Oon, co2on, n2Oon

    integer(kint) :: i,k,ne,k1,k2,ik,isb, ibn, klow

    !---------------------------------------------------------------------------!
    ! compute layer pressure (pa) and layer temperature minus 250K (dt250)      !
    !-------------------------------------------------------------------------!

     do k=1,nplevel-1
       do i=nxl,nxr
          pa(i,k)=r1h * (pl(i,k) + pl(i,k+1))
          dt250(i,k)=ta(i,k) - T250
       enddo
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

    do  i=nxl,nxr
       do k=1, nplevel-1

          dp   (i,k) = pl(i,k+1)-pl(i,k)
          
          dh2o (i,k) = 1.02_kreal * wa(i,k)*dp(i,k)
          dh2o (i,k) = max(dh2o(i,k),1.e-8_kreal)
          do3  (i,k) = 476.0_kreal* oa(i,k)*dp(i,k)
          do3  (i,k) = max(do3 (i,k),1.e-6_kreal)
          dco2 (i,k) = C789 * co2lw*dp(i,k)
          dco2 (i,k) = max(dco2(i,k),1.e-4_kreal)
          
          dch4 (i,k) = C789 *ch4*dp(i,k)
          dn2o (i,k) = C789 *n2o*dp(i,k)
          df11 (i,k) = C789 *cfc11*dp(i,k)
          df12 (i,k) = C789 *cfc12*dp(i,k)
          df22 (i,k) = C789 *cfc22*dp(i,k)

          !-------------------------------------------------------------------!
          ! compute scaled water vapor amount for h2o continuum absorption    !
          !  following eq. (4.21).                                            !
          !-------------------------------------------------------------------!

          xx=pa(i,k)*0.001618_kreal*wa(i,k)*wa(i,k)*dp(i,k)
          dcont(i,k) = xx*exp(1800.0_kreal/ta(i,k)-6.081_kreal)

       enddo
    enddo

    !------------------------------------------------------------------------!
    ! compute layer cloud water amount (gm/m**2)                             !  
    !  index is 1 for ice, 2 for waterdrops and 3 for raindrops.             !
    !------------------------------------------------------------------------!

    do i=nxl,nxr       !(*)
       klow = ifeld(i,j)
       do  k=1,nplevel + 1 -klow
          xx=1.02_kreal*10000.0_kreal *(pl(i,k+1)-pl(i,k))
          cwp(i,k,1)=xx*cwc(i,k,1)
          cwp(i,k,2)=xx*cwc(i,k,2)
          cwp(i,k,3)=xx*cwc(i,k,3)
       enddo
    enddo

    !------------------------------------------------------------------------!
    ! The surface (nplevel) is treated as a layer filled with black clouds.  !
    ! transfc is the transmittance between the surface and a pressure level. !
    ! trantcr is the clear-sky transmittance between the surface and a       !
    ! pressure level.                                                        !
    !------------------------------------------------------------------------!
 
    do i=nxl,nxr       !(*)
        klow = ifeld(i,j)
       transfc(i,nplevel-klow+2)=r1
    enddo

    !------------------------------------------------------------------------!
    ! initialize fluxes                                                      !
    !------------------------------------------------------------------------!

    do k=1,nplevel
       flx(nxl:nxr,k) = r0
    enddo
    
    !------------------------------------------------------------------------!
    ! reset total ground irradiance for x-line                               !
    !------------------------------------------------------------------------!

    flxdn_lw_gr_tot(:,j) = r0
  
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

       !---------------------------------------------------------------------!
       ! blayer is the spectrally integrated planck flux of the mean layer   !
       ! temperature derived from eq. (3.11)                                 !
       ! The fitting for the planck flux is valid for the range 160-345 K.   !
       !---------------------------------------------------------------------!

       do k=1,nplevel-1

          tx(nxl:nxr)=ta(nxl:nxr,k)
                       
          call planck(ibn,tx,xlayer, nxl, nxr, nx)

          blayer(nxl:nxr,k)=xlayer(nxl:nxr)
            
       enddo

       !---------------------------------------------------------------------!
       ! Index "0" is the layer above the top of the atmosphere.             !
       !---------------------------------------------------------------------!
         
       blayer(nxl:nxr,0)=r0
         
       !---------------------------------------------------------------------!
       ! Surface emission and reflectivity. See Section 9.                   !
       !---------------------------------------------------------------------!

       call sfcflux (ibn,tg,eg,bs,rflxs, nxl, nxr, nx) 
       
       do i= nxl, nxr       !(*)
          klow = ifeld(i,j)
          blayer(i,nplevel-klow+2)=bs(i)
       end do
       !--------------------------------------------------------------------!
       ! interpolate Planck function at model levels (linear in p)          ! 
       !--------------------------------------------------------------------!

       do i=nxl,nxr       !(*)
          klow = ifeld(i,j)
          do  k=2,nplevel-klow+1
             blevel(i,k)=(blayer(i,k-1)*dp(i,k)+blayer(i,k)*dp(i,k-1))     &
                          /(dp(i,k-1)+dp(i,k))
          enddo
       enddo

       !--------------------------------------------------------------------!
       ! Extrapolate blevel(i,1) from blayer(i,2) and blayer(i,1)           !
       !--------------------------------------------------------------------!

       do i=nxl,nxr
          blevel(i,1)=blayer(i,1)+(blayer(i,1)-blayer(i,2))*dp(i,1)        &
                       /(dp(i,1)+dp(i,2))
       enddo
         
       !----------------------------------------------------------------------!
       ! If the surface air temperature tb is known, compute blevel(i,nplevel)!
       !----------------------------------------------------------------------!

       call planck(ibn,tb,xlayer, nxl, nxr, nx)
         
       do i= nxl,nxr       !(*)
          klow = ifeld(i,j)
          blevel(i,nplevel-klow+2)=xlayer(i)
       enddo
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


       !----------------------------------------------------------------------!
       ! Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)    !
       ! Rain optical thickness is set to 0.00307 /(gm/m**2).                 !
       ! It is for a specific drop size distribution provided by Q. Fu.       !
       !---------------------------------------------------------------- -----!

       do i=nxl,nxr       !(*)
          klow = ifeld(i,j)
          do  k=1,nplevel-klow+1
             taucl(i,k,1)=cwp(i,k,1)*(aib(1,ibn)+aib(2,ibn)                 &
                             /reff(i,k,1)**aib(3,ibn))
             taucl(i,k,2)=cwp(i,k,2)*(awb(1,ibn)+(awb(2,ibn)                &
                             +(awb(3,ibn)+awb(4,ibn)*reff(i,k,2))*reff(i,k,2))&
                             *reff(i,k,2))
             taucl(i,k,3)=0.00307_kreal*cwp(i,k,3)
          enddo
       enddo

       !----------------------------------------------------------------------
       ! Compute cloud single-scattering albedo and asymmetry factor for      !
       !     a mixture of ice particles and liquid drops following            !
       !     Eqs. (6.5), (6.6), (6.9) and (6.10).                             !
       !     Single-scattering albedo and asymmetry factor of rain are set    !
       !     to 0.54 and 0.95, respectively, based on the information provided!
       !     by Prof. Qiang Fu.                                               !
       !----------------------------------------------------------------------!

       do i=nxl,nxr              !(*)
          klow = ifeld(i,j)
          do  k=1,nplevel-klow+1 

             tcldlyr(i,k) = r1
             tauc=taucl(i,k,1)+taucl(i,k,2)+taucl(i,k,3)

             if (tauc.gt.0.02 .and. fcld(i,k).gt.0.01) then
                reff1=min(reff(i,k,1),130.0_kreal*r1)
                reff2=min(reff(i,k,2),20.0_kreal*r1)
                
                w1=taucl(i,k,1)*(aiw(1,ibn)+(aiw(2,ibn)+(aiw(3,ibn)      &
                       +aiw(4,ibn)*reff1)*reff1)*reff1)
                w2=taucl(i,k,2)*(aww(1,ibn)+(aww(2,ibn)+(aww(3,ibn)      &
                       +aww(4,ibn)*reff2)*reff2)*reff2) 
                w3=taucl(i,k,3)*0.54_kreal
                ww=(w1+w2+w3)/tauc

                g1=w1*(aig(1,ibn)+(aig(2,ibn)+(aig(3,ibn)                &
                       +aig(4,ibn)*reff1)*reff1)*reff1)
                g2=w2*(awg(1,ibn)+(awg(2,ibn)+(awg(3,ibn)                &
                       +awg(4,ibn)*reff2)*reff2)*reff2)
                g3=w3*0.95_kreal

                gg=(g1+g2+g3)/(w1+w2+w3)

                !-----------------------------------------------------------!
                ! Parameterization of LW scattering following Eqs. (6.11)   !
                ! and (6.12).                                               !
                !-----------------------------------------------------------!

                ff=r1h+(0.3739_kreal+(0.0076_kreal+0.1185_kreal*gg)*gg)*gg
                tauc=(r1-ww*ff)*tauc

                !-----------------------------------------------------------!
                ! compute cloud diffuse transmittance. It is approximated   !
                ! by using                                                  !
                ! a diffusivity factor of 1.66.                             !
                !-----------------------------------------------------------!

                tcldlyr(i,k)=exp(-1.66_kreal*tauc)
             endif
               
          enddo
       enddo

       !---------------------------------------------------------------------!
       ! Compute the exponential terms (Eq. 8.21) at each layer due to       !
       ! water vapor line absorption when k-distribution is used             !
       !---------------------------------------------------------------------!

       if (.not.h2otbl .and. .not.b10bnd .and. h2Oon) then
          call h2oexps(ibn,dh2o,pa,dt250,h2oexp,j)
       endif

       !--------------------------------------------------------------------!
       ! compute the exponential terms (Eq. 4.24) at each layer due to      !
       ! water vapor continuum absorption.                                  !
       ! ne is the number of terms used in each band to compute water       !
       ! vapor continuum transmittance (Table 9).                           !
       !--------------------------------------------------------------------!

       ne=0
       if (conbnd) then

          ne=1
          if (ibn.eq.3) ne=3
          call conexps(ibn,dcont,conexp,j)

       endif

       !------------------------------------------------------------------!
       ! compute the exponential terms (Eq. 8.21) at each layer due to    !
       ! co2 absorption                                                   !
       !------------------------------------------------------------------!

       if (co2bnd) then
          call co2exps(dco2,pa,dt250,co2exp,j)
       endif

       !------------------------------------------------------------------!
       ! ***** for trace gases *****                                      !
       !------------------------------------------------------------------!


       !------------------------------------------------------------------!
       ! compute the exponential terms at each layer due to n2o absorption!
       !------------------------------------------------------------------!

       if (n2obnd) then
          call n2oexps(ibn,dn2o,pa,dt250,n2oexp,j)
       endif

       !------------------------------------------------------------------!
       ! compute the exponential terms at each layer due to ch4 absorption!
       !------------------------------------------------------------------!

       if (ch4bnd) then
          call ch4exps(ibn,dch4,pa,dt250,ch4exp,j)
       endif
       !------------------------------------------------------------------!
       ! Compute the exponential terms due to co2 minor absorption        !     
       !------------------------------------------------------------------!
       
       if (combnd) then
          call comexps(ibn,dco2,dt250,comexp,j)
       endif

       !-------------------------------------------------------- ---------!
       ! Compute the exponential terms due to cfc11 absorption.           !
       ! The values of the parameters are given in Table 7.               !
       !------------------------------------------------------------------!

       if (f11bnd) then
          a1  = 1.26610e-3_kreal
          b1  = 3.55940e-6_kreal
          fk1 = 1.89736e+1_kreal
          a2  = 8.19370e-4_kreal
          b2  = 4.67810e-6_kreal
          fk2 = 1.01487e+1_kreal
          call cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,df11,dt250,f11exp,j)
       endif

       !------------------------------------------------------------------!
       ! Compute the exponential terms due to cfc12 absorption.           !
       !------------------------------------------------------------------!

       if (f12bnd) then
          a1  = 8.77370e-4_kreal
          b1  =-5.88440e-6_kreal
          fk1 = 1.58104e+1_kreal
          a2  = 8.62000e-4_kreal
          b2  =-4.22500e-6_kreal
          fk2 = 3.70107e+1_kreal
          call cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,df12,dt250,f12exp,j)
       endif

       !-----------------------------------------------------------------!
       ! Compute the exponential terms due to cfc22 absorption.          !
       !-----------------------------------------------------------------!

       if (f22bnd) then
          a1  = 9.65130e-4_kreal
          b1  = 1.31280e-5_kreal
          fk1 = 6.18536e+0_kreal
          a2  =-3.00010e-5_kreal
          b2  = 5.25010e-7_kreal
          fk2 = 3.27912e+1_kreal
          call cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,df22,dt250,f22exp,j)
       endif

       !------------------------------------------------------------------!
       ! Compute the exponential terms at each layer in band 10 due to    !
       !     h2o line and continuum, co2, and n2o absorption              !
       !------------------------------------------------------------------!

       if (b10bnd) then
          call b10exps(dh2o,dcont,dco2,dn2o,pa,dt250,       &
                        h2oexp,conexp,co2exp,n2oexp, h2Oon, co2on, n2Oon,j)
       endif


       !------------------------------------------------------------------!
       ! initialize fluxes                                                !
       !------------------------------------------------------------------!

       do k=1,nplevel
             flxu(nxl:nxr,k) = r0
             flxd(nxl:nxr,k) = r0
       enddo

       !----------------------------------------------------------------------!
       ! For a given level, k1, compute the transmittance between this level  !
       ! all the levels below, trant(i,k2).                                   !
       ! Also, compute the upward and doward blackbody emissions of a layer,  !
       ! bu and bd.                                                           !
       !----------------------------------------------------------------------!

       bd(:,0) = r0    ! (*)

       
       do 2000 k1=1,nplevel-1

          !-----------------------------------------------------------------!
          ! initialization...                                               !
          !-----------------------------------------------------------------!
          
          !-----------------------------------------------------------------!
          ! for h2o line transmission                                       ! 
          !-----------------------------------------------------------------!
          
          if (h2Oon .and. .not. h2otbl) then
             do ik=1,6
                th2o(nxl:nxr,ik)=r1
             enddo
          endif

          !-----------------------------------------------------------------!
          ! for h2o continuum transmission                                  !
          !-----------------------------------------------------------------!
          tcon(nxl:nxr,1:3) = r0
          if (conbnd) then
             do ik=1,3
                tcon(nxl:nxr,ik)=r1
             enddo
          endif

          !-----------------------------------------------------------------!
          ! for co2 transmission using k-distribution method.               !
          ! band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c     !
          ! are combined in computing the co2 transmittance.                !
          !-----------------------------------------------------------------!

          if (co2bnd) then
             do isb=1,2
                do ik=1,6
                   tco2(nxl:nxr,ik,isb)=r1
                enddo
             enddo
          endif

          !-----------------------------------------------------------------!
          ! ***** for trace gases *****                                     !
          !-----------------------------------------------------------------!

          !-----------------------------------------------------------------!  
          ! for n2o transmission using k-distribution method.               !
          !-----------------------------------------------------------------!

          if (n2obnd) then
             do ik=1,4
                tn2o(nxl:nxr,ik)=r1
             enddo
          endif

          !-------------------------------------------------------------!
          ! for ch4 transmission using k-distribution method.           !
          !-------------------------------------------------------------!

          if (ch4bnd) then
             do ik=1,4
                   tch4(nxl:nxr,ik)=r1
             enddo
          endif

          !----------------------------------------------------------------!
          ! for co2-minor transmission using k-distribution method.        !
          !----------------------------------------------------------------!

          if (combnd) then
             do ik=1,6
                tcom(nxl:nxr,ik)=r1
             enddo
          endif

          !----------------------------------------------------------------!
          ! or cfc-11 transmission using k-distribution method.            !
          !----------------------------------------------------------------!

          if (f11bnd) then
                tf11(nxl:nxr)=r1
          endif

          !----------------------------------------------------------------!
          ! for cfc-12 transmission using k-distribution method.           !
          !----------------------------------------------------------------!

          if (f12bnd) then
             tf12(nxl:nxr)=r1
          endif

          !----------------------------------------------------------------!
          ! for cfc-22 transmission when using k-distribution method.      !
          !----------------------------------------------------------------!

          if (f22bnd) then
             tf22(nxl:nxr)=r1
          endif

          !----------------------------------------------------------------!
          ! or the transmission in band 10 using k-distribution method.    !
          !----------------------------------------------------------------!

          if (b10bnd) then
             do ik=1,5
                th2o(nxl:nxr,ik)=r1
             enddo
             
             do ik=1,6
                tco2(nxl:nxr,ik,1)=r1
             enddo
             
             tcon(nxl:nxr,1)=r1
                          
             do ik=1,2
                tn2o(nxl:nxr,ik)=r1
             enddo
          endif

          !----------------------------------------------------------------!
          ! ***** end trace gases *****                                    !
          !----------------------------------------------------------------!

          x1(nxl:nxr)=r0
          x2(nxl:nxr)=r0
          x3(nxl:nxr)=r0
     

          do k=1,nplevel
             fclr(nxl:nxr,k)=r1
          enddo

          !--------------------------------------------------------------!
          ! loop over the bottom level of the region (k2)                !   
          !--------------------------------------------------------------!

          do 3000 k2=k1+1,nplevel

             !-------------------------------------------------------------!
             ! trant is the total transmittance between levels k1 and k2.  !
             !-------------------------------------------------------------!

             trant(nxl:nxr,k2)=r1
             
             if (h2otbl) then

                !----------------------------------------------------------!
                ! Compute water vapor transmittance using table look-up.   !
                ! The following values are taken from Table 8.             !
                !----------------------------------------------------------!

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
                   do i=nxl,nxr
                      tcon(i,1)=tcon(i,1)*conexp(i,k2-1,1)
                      trant(i,k2)=trant(i,k2)*tcon(i,1)
                   enddo
                endif

             else

                !-------------------------------------------------------------!
                ! compute water vapor transmittance using k-distribution      !
                !-------------------------------------------------------------!

                if (h2Oon .and. .not.b10bnd) then
                   call h2okdis(ibn,k2-1,ne,h2oexp,conexp,     &
                                th2o,tcon,trant)
                endif
               
             endif


             if (co2bnd) then

                !--------------------------------------------------------------!
                ! Compute co2 transmittance using table k-distribution method. !
                !--------------------------------------------------------------!

                call co2kdis(k2-1,co2exp,tco2,trant)

             endif

             !---------------------------------------------------------------------------!
             ! Always use table look-up to compute o3 transmittance.                     !
             ! The following values are taken from Table 8.                              !
             !---------------------------------------------------------------------------!

             if (oznbnd) then
                w1=-6.0
                p1=-2.0
                dwe=0.3
                dpe=0.2
                call tablup(k2,nxp,nh,do3,pa,dt250,x1,x2,x3,            &
                             w1,p1,dwe,dpe,no,o1,o2,o3,trant)

             endif



             !---------------------------------------------------------------------------!
             ! ***** for trace gases *****                                               !
             !---------------------------------------------------------------------------!

             !---------------------------------------------------------------------------!
             ! compute n2o transmittance using k-distribution method                     !
             !---------------------------------------------------------------------------!
             if (n2obnd) then
                call n2okdis(ibn,k2-1,n2oexp,tn2o,trant)
             endif

             !---------------------------------------------------------------------------!
             ! compute ch4 transmittance using k-distribution method                     !
             !---------------------------------------------------------------------------!
             if (ch4bnd) then
                call ch4kdis(ibn,k2-1,ch4exp,tch4,trant)
             endif

             !---------------------------------------------------------------------------!
             ! compute co2-minor transmittance using k-distribution method               !
             !---------------------------------------------------------------------------!
             if (combnd) then
                call comkdis(ibn,k2-1,comexp,tcom,trant)
             endif

             !---------------------------------------------------------------------------!
             ! compute cfc11 transmittance using k-distribution method                   !
             !---------------------------------------------------------------------------!

             if (f11bnd) then
                call cfckdis(k2-1,f11exp,tf11,trant)
             endif

             !---------------------------------------------------------------------------!
             ! compute cfc12 transmittance using k-distribution method                   !
             !---------------------------------------------------------------------------!

             if (f12bnd) then
                call cfckdis(k2-1,f12exp,tf12,trant)
             endif

             !---------------------------------------------------------------------------!
             ! compute cfc22 transmittance using k-distribution method                   !
             !---------------------------------------------------------------------------!

             if (f22bnd) then
                call cfckdis(k2-1,f22exp,tf22,trant)
             endif

             !---------------------------------------------------------------------------!
             ! Compute transmittance in band 10 using k-distribution method.             !
             ! For band 10, trant is the change in transmittance due to n2o absorption   !
             !---------------------------------------------------------------------------!

             if (b10bnd) then
                call b10kdis(k2-1,h2oexp,conexp,co2exp,n2oexp,        &
                             th2o,tcon,tco2,tn2o,trant, h2Oon,co2on, n2Oon)

             endif

             !---------------------------------------------------------------------------!
             ! *****   end trace gases  *****                                            !
             !---------------------------------------------------------------------------!

             do i=nxl,nxr
                !------------------------------------------------------------------------!
                ! > original code:                                                       !	
                !        fclr(i,k2)=fclr(i,k2)*tcldlyr(i,k2-1)                           !
                !> after Yang Chen's discovery:                                          !
                !------------------------------------------------------------------------!

                fclr(i,k2)=fclr(i,k2-1)*tcldlyr(i,k2-1)
             enddo

             !-----------------------------------------------------------------------!
             ! Compute upward and downward blackbody emission of a layer             !
             !-----------------------------------------------------------------------!

             if (k1.eq.1) then
                do i=nxl,nxr

                   xx=(blayer(i,k2-1)-blevel(i,k2-1)) * (blayer(i,k2-1)-blevel(i,k2))

                   if (xx.gt.0.0) then

                      !--------------------------------------------------------------------------!
                      ! If xx>0, there is a local temperature minimum or maximum.                !
                      ! Computations of bd and bu follow Eq. (8.20).                             !
                      !--------------------------------------------------------------------------!

                      bd(i,k2-1)=r1h*blayer(i,k2-1)+r1q*(blevel(i,k2-1)+blevel(i,k2))
                      bu(i,k2-1)=bd(i,k2-1)

                   else

                      !--------------------------------------------------------------------------!
                      ! Computations of bd and bu following Eqs.(8.17) and (8.18).               !
                      ! The effect of clouds on the transmission of a layer is taken             !
                      ! into account, following Eq. (8.19).                                      !
                      !--------------------------------------------------------------------------!

                      xx=(fcld(i,k2-1)*tcldlyr(i,k2-1)+(r1-fcld(i,k2-1))) * trant(i,k2)
                      yy=min(0.9999_kreal,xx)
                      yy=max(0.00001_kreal,yy)
                      xx=(blevel(i,k2-1)-blevel(i,k2))/log(yy)
                      bd(i,k2-1)=(blevel(i,k2)-blevel(i,k2-1)*yy)/(r1-yy)-xx
                      bu(i,k2-1)=(blevel(i,k2-1)+blevel(i,k2))-bd(i,k2-1)

                   endif

                enddo

                do i=nxl,nxr         !(*)
                   klow = ifeld(i,j)
                   bu(i,nplevel-klow+2)=blayer(i,nplevel-klow+2)
                enddo
 
             endif
            
3000      end do

          !-------------------------------------------------------------------------!
          ! upward and downward flux calculations.                                  !
          !-------------------------------------------------------------------------!

          do 4000 k2=k1+1,nplevel

             if (k2 <= npdiff+1) then      ! (*)
                flag(nxl:nxr) = r1
             else
                flag(nxl:nxr) = iflgs(nxl:nxr,j,nplevel-k2+2)
             end if

             if (k2.eq.k1+1 .and. ibn .ne. 10) then

                !----------------------------------------------------------------------------!
                ! The first terms on the rhs of Eqs. (8.15) and (8.16)                       !
                !----------------------------------------------------------------------------!
              
                do i=nxl,nxr
                   flxu(i,k1)=flxu(i,k1)-bu(i,k1)*flag(i) ! (*)
                   flxd(i,k2)=flxd(i,k2)+bd(i,k1)*flag(i)
                   
                enddo

             endif

             !---------------------------------------------------------------------------!
             ! The summation terms on the rhs of Eqs. (8.15) and (8.16).                 !
             ! Also see Eqs. (5.4) and (5.5) for Band 10.                                !
             !---------------------------------------------------------------------------!

             do i=nxl,nxr
                xx=trant(i,k2)*(bu(i,k2-1)-bu(i,k2))*flag(i) ! (*)
                flxu(i,k1) =flxu(i,k1)+xx*fclr(i,k2)
                xx=trant(i,k2)*(bd(i,k1-1)-bd(i,k1))*flag(i)
                flxd(i,k2) =flxd(i,k2)+xx*fclr(i,k2)

             enddo

4000      end do

          !---------------------------------------------------------------------------!
          ! Here, fclr and trant are, respectively, the clear line-of-sight           !
          ! and the transmittance between k1 and the surface.                         !
          !---------------------------------------------------------------------------!

          do i=nxl,nxr    ! (*)
             klow = ifeld(i,j)
             transfc(i,k1) =trant(i,nplevel-klow+2)*fclr(i,nplevel-klow+2)
          enddo


2000   end do

       if (.not. b10bnd) then

          !---------------------------------------------------------------------------!
          ! For surface emission.                                                     !
          ! Note: blayer(i,nplevel) and dbs include the surface emissivity effect.    !
          !---------------------------------------------------------------------------!

          do i=nxl,nxr
             klow = ifeld(i,j)    ! (*)
             flxu(i,nplevel-klow+2)=-blayer(i,nplevel-klow+2)
          enddo

          !---------------------------------------------------------------------------!
          ! Add the flux reflected by the surface.(Second term on the rhs of Eq. 8.16)!
          !---------------------------------------------------------------------------!

          do i=nxl,nxr
             klow = ifeld(i,j)    
             do k=1,nplevel-klow+1    
                flxu(i,k)=flxu(i,k)- flxd(i,nplevel-klow+2)*transfc(i,k)*rflxs(i)  
             enddo
          enddo

       endif

       !------------------------------------------------------------------!
       ! Summation of fluxes over spectral bands                          !
       !-------------------------------------------------------- ---------!

       do i=nxl,nxr
          klow = ifeld(i,j)   ! (*)  
          do k=1,nplevel-klow+2
             flx(i,k)=flx(i,k)+flxd(i,k)+flxu(i,k)
             !--------------------------------------------------------------!
             ! extract and add total ground irradiance (for each band)      !
             !--------------------------------------------------------------!  
             if (k==nplevel-klow+2) flxdn_lw_gr_tot(i,j) = &
                                     flxdn_lw_gr_tot(i,j) + flxd(i,k)

          enddo
       enddo

1000 end do
      
  end subroutine irrad

!============================================================================================!

  subroutine planck(ibn,t,xlayer, nxl, nxr, nx)

    !------------------------------------------------------------------------!
    ! Compute spectrally integrated Planck flux                              !
    !------------------------------------------------------------------------!


    integer(kint), intent(in) :: nxl, nxr, nx
    integer(kint), intent(in) :: ibn                          ! spectral band index
    real(kreal), intent(in), dimension(nx) :: t               ! temperature (K)
    real(kreal), intent(inout), dimension(nx) ::  xlayer      ! planck flux (w/m2)

    integer(kint) :: i

    do i=nxl,nxr
       xlayer(i)=(t(i)*(t(i)*(t(i)*(t(i)*(t(i)*cb(6,ibn)+cb(5,ibn))  &
                 +cb(4,ibn))+cb(3,ibn))+cb(2,ibn)))+cb(1,ibn)
    enddo

  end subroutine planck

!=============================================================================================!

  subroutine sfcflux (ibn,tg,eg,bs,rflxs, nxl, nxr, nx)

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

    use phys_constants, only: r1

    integer(kint), intent(in) :: nxl, nxr, nx
    integer(kint), intent(in) :: ibn
    real(kreal), intent(in) :: tg(nx),eg(nx,10)

    real(kreal), intent(inout) ::  bs(nx),rflxs(nx)

    integer(kint) :: i
    real(kreal) :: tx(nx)

    !------------------------------------------------------------------------------------!
    ! For homogeneous surface without vegetation following Eqs. (9.4), (9.5), and (3.13) !
    !------------------------------------------------------------------------------------!

    do i=nxl,nxr
       tx(i)=tg(i)
    enddo

    call planck(ibn,tx,bs,nxl, nxr, nx)

    do i=nxl,nxr
       rflxs(i)=r1-eg(i,ibn)
    enddo
       
  end subroutine sfcflux

!==============================================================================================!

  subroutine h2oexps(ibn,dh2o,pa,dt250,h2oexp,j)

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

    use phys_constants, only: r1
    use atham_module ,  only: nxl,nxr
    use atham_module,   only: nx, ifeld

    integer(kint), intent(in) :: ibn, j
    real(kreal), intent(in) :: dh2o(nx,nplevel-1), pa(nx,nplevel-1), dt250(nx,nplevel-1)
      
    real(kreal), intent(inout), dimension(nx, nplevel-1, 6) :: h2oexp

    real(kreal), parameter :: pref = 500._kreal
    integer(kint) :: i, k, ik, klow
    real(kreal) :: xh

    !--------------------------------------------------------------------------!
    ! Note that the 3 sub-bands in band 3 use the same set of xkw, aw,         !
    ! and bw,  therefore, h2oexp for these sub-bands are identical.            !
    !--------------------------------------------------------------------------!

       do i=nxl,nxr
          klow = ifeld(i,j)      ! (*)
          do k = 1, nplevel-klow+1
          !----------------------------------------------------------------------------------!
          ! xh is the scaled water vapor amount for line absorption computed from Eq. (4.4). !
          !----------------------------------------------------------------------------------!

          xh = dh2o(i,k)*(pa(i,k)/pref)**pm(ibn)         &
               *(r1+(aw(ibn)+bw(ibn)* dt250(i,k))*dt250(i,k))

          !----------------------------------------------------------------------------------!
          ! h2oexp is the water vapor transmittance of the layer k due to line absorption    !
          !----------------------------------------------------------------------------------!

          h2oexp(i,k,1) = exp(-xh*xkw(ibn))

       enddo
    enddo

    !---------------------------------------------------------------------------!
    ! compute transmittances from Eq. (8.22)                                    !
    !---------------------------------------------------------------------------!

    do ik=2,6

       if (mw(ibn).eq.6) then

       do i=nxl,nxr
          klow = ifeld(i,j)   ! (*)
          do k = 1, nplevel-klow+1
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                h2oexp(i,k,ik) = xh*xh*xh
             enddo
          enddo

       elseif (mw(ibn).eq.8) then

       do i=nxl,nxr
          klow = ifeld(i,j)    ! (*)
          do k = 1, nplevel-klow+1
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                xh = xh*xh
                h2oexp(i,k,ik) = xh*xh
             enddo
          enddo

       elseif (mw(ibn).eq.9) then

       do i=nxl,nxr
          klow = ifeld(i,j)   ! (*)
          do k = 1, nplevel-klow+1
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                h2oexp(i,k,ik) = xh*xh*xh
             enddo
          enddo

       else

       do i=nxl,nxr
          klow = ifeld(i,j)    ! (*)
          do k = 1, nplevel-klow+1
                xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
                xh = xh*xh
                xh = xh*xh
                h2oexp(i,k,ik) = xh*xh
             enddo
          enddo

       endif
    enddo

  end subroutine h2oexps

!==========================================================================================!

  subroutine conexps(ibn,dcont,conexp,j)

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

    use atham_module,  only: nx, nxl,nxr, ifeld

    real(kreal), intent(in) :: dcont(nx,nplevel-1)
    integer(kint), intent(in) :: ibn,j
    
    real(kreal), intent(inout) :: conexp(nx,nplevel-1,3)

    integer(kint) :: i, k, klow
      
       do i=nxl,nxr
          klow = ifeld(i,j)   ! (*)
          do k = 1, nplevel-klow+1
          conexp(i,k,1) = exp(-dcont(i,k)*xke(ibn))
       enddo
    enddo

    if (ibn .eq. 3) then

       !------------------------------------------------------------------------------!
       ! The absorption coefficients for sub-bands 3b and 3a are, respectively,       !
       ! two and four times the absorption coefficient for sub-band 3c (Table 9).     !
       ! Note that conexp(i,k,3) is for sub-band 3a.                                  !
       !------------------------------------------------------------------------------!
       
       do i=nxl,nxr
          klow = ifeld(i,j)    ! (*)
          do k = 1, nplevel-klow+1
             conexp(i,k,2) = conexp(i,k,1) *conexp(i,k,1)
             conexp(i,k,3) = conexp(i,k,2) *conexp(i,k,2)
          enddo
       enddo
       
    endif

  end subroutine conexps

!===========================================================================================!

  subroutine co2exps(dco2,pa,dt250,co2exp,j)

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

    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx, ifeld
      
    real(kreal), intent(in) ::  dco2(nx,nplevel-1), pa(nx,nplevel-1), dt250(nx,nplevel-1)
    real(kreal), intent(inout) :: co2exp(nx,nplevel-1,6,2)
    integer(kint), intent(in) :: j
    
    integer(kint) :: i, k, klow
    real(kreal) :: xc

       do i=nxl,nxr
          klow = ifeld(i,j)    ! (*)
          do k = 1, nplevel-klow+1

          !---------------------------------------------------------------------------!
          ! The scaling parameters are given in Table 3, and values of the absorption !
          ! coefficient are given in Table 10.                                        !
          !                                                                           !
          ! Scaled co2 amount for band-wings (sub-bands 3a and 3c)                    !
          !---------------------------------------------------------------------------!

          xc = dco2(i,k)*((pa(i,k)/300.0_kreal)**r1h)            &
               *(r1+(0.0182_kreal+1.07e-4_kreal*dt250(i,k))*dt250(i,k))

          !---------------------------------------------------------------------------!
          ! six exponentials by powers of 8 (See Eqs. 8.21, 8.22 and Table 10).       !
          !---------------------------------------------------------------------------!

          co2exp(i,k,1,1)=exp(-xc*2.656e-5_kreal)
          
          xc=co2exp(i,k,1,1)*co2exp(i,k,1,1)
          xc=xc*xc
          co2exp(i,k,2,1)=xc*xc

          xc=co2exp(i,k,2,1)*co2exp(i,k,2,1)
          xc=xc*xc
          co2exp(i,k,3,1)=xc*xc

          xc=co2exp(i,k,3,1)*co2exp(i,k,3,1)
          xc=xc*xc
          co2exp(i,k,4,1)=xc*xc

          xc=co2exp(i,k,4,1)*co2exp(i,k,4,1)
          xc=xc*xc
          co2exp(i,k,5,1)=xc*xc

          xc=co2exp(i,k,5,1)*co2exp(i,k,5,1)
          xc=xc*xc
          co2exp(i,k,6,1)=xc*xc

          !-------------------------------------------------------------------------!
          ! For band-center region (sub-band 3b)                                    !
          !-------------------------------------------------------------------------!

          xc = dco2(i,k)*(pa(i,k)/30.0_kreal)**0.85_kreal   &
               *(r1+(0.0042_kreal+2.00e-5_kreal*dt250(i,k))*dt250(i,k))

          co2exp(i,k,1,2)=exp(-xc*2.656e-3_kreal)
            
          xc=co2exp(i,k,1,2)*co2exp(i,k,1,2)
          xc=xc*xc
          co2exp(i,k,2,2)=xc*xc
            
          xc=co2exp(i,k,2,2)*co2exp(i,k,2,2)
          xc=xc*xc
          co2exp(i,k,3,2)=xc*xc

          xc=co2exp(i,k,3,2)*co2exp(i,k,3,2)
          xc=xc*xc
          co2exp(i,k,4,2)=xc*xc

          xc=co2exp(i,k,4,2)*co2exp(i,k,4,2)
          xc=xc*xc
          co2exp(i,k,5,2)=xc*xc

          xc=co2exp(i,k,5,2)*co2exp(i,k,5,2)
          xc=xc*xc
          co2exp(i,k,6,2)=xc*xc

       enddo
    enddo

  end subroutine co2exps

!=====================================================================================!

  subroutine n2oexps(ibn,dn2o,pa,dt250,n2oexp,j)

    !------------------------------------------------------------------------!
    ! Compute n2o exponentials for individual layers                         !
    !                                                                        !
    ! Output parameters: n2oexp  =  2 or 4 exponentials for each layer       !
    !------------------------------------------------------------------------!
    use phys_constants,only: r1
    use atham_module,  only: nxl, nxr, nx, ifeld
      
    integer(kint), intent(in) :: ibn, j
    real(kreal), intent(in) ::  dn2o(nx,nplevel-1), pa(nx,nplevel-1), dt250(nx,nplevel-1)

    real(kreal), intent(inout) :: n2oexp(nx,nplevel-1,4)
    
    integer(kint) :: i, k, klow
    real(kreal) :: xc, xc1, xc2

    !-------------------------------------------------------------------------!
    ! Scaling and absorption data are given in Table 5.                       !
    ! Transmittances are computed using Eqs. (8.21) and (8.22).               !
    !-------------------------------------------------------------------------!

       do i=nxl,nxr
          klow = ifeld(i,j)    ! (*)
          do k = 1, nplevel-klow+1

          !---------------------------------------------------------------!
          ! Four exponential by powers of 21 for band 6.                  !
          !---------------------------------------------------------------!

          if (ibn .eq. 6) then

             xc=dn2o(i,k)*(r1+(1.9297e-3_kreal+4.3750e-6_kreal*dt250(i,k))*dt250(i,k))
             n2oexp(i,k,1)=exp(-xc*6.31582e-2_kreal)
             
             xc=n2oexp(i,k,1)*n2oexp(i,k,1)*n2oexp(i,k,1)
             xc1=xc*xc
             xc2=xc1*xc1
             n2oexp(i,k,2)=xc*xc1*xc2        
          else

             !------------------------------------------------------------!
             ! four exponential by powers of 8 for band 7                 !
             !------------------------------------------------------------!

             xc=dn2o(i,k)*(pa(i,k)/500.0_kreal)**0.48_kreal *    &
                    (r1+(1.3804e-3_kreal+7.4838e-6_kreal*dt250(i,k))*dt250(i,k))
             n2oexp(i,k,1)=exp(-xc*5.35779e-2_kreal)

             xc=n2oexp(i,k,1)*n2oexp(i,k,1)
             xc=xc*xc
             n2oexp(i,k,2)=xc*xc
             xc=n2oexp(i,k,2)*n2oexp(i,k,2)
             xc=xc*xc
             n2oexp(i,k,3)=xc*xc
             xc=n2oexp(i,k,3)*n2oexp(i,k,3)
             xc=xc*xc
             n2oexp(i,k,4)=xc*xc

          endif

       enddo
    enddo

  end subroutine n2oexps

!=============================================================================!

  subroutine ch4exps(ibn,dch4,pa,dt250,ch4exp,j)

    !-------------------------------------------------------------------------!
    ! Compute ch4 exponentials for individual layers                          !
    !                                                                         !
    ! Output parameters : ch2exp = 1 or 4 exponentials for each layer         !
    !-------------------------------------------------------------------------!
    use phys_constants,only: r1
    use atham_module,  only: nxl, nxr, nx, ifeld
      
    integer(kint), intent(in) :: ibn,j
    real(kreal), intent(in) ::  dch4(nx,nplevel-1), pa(nx,nplevel-1), dt250(nx,nplevel-1)

    real(kreal), intent(inout) :: ch4exp(nx,nplevel-1,4)
    
    integer(kint) :: i, k, klow
    real(kreal) :: xc

    !-------------------------------------------------------------------------!
    ! Scaling and absorpton data are given in Table 5                         !
    !-------------------------------------------------------------------------!

       do i=nxl,nxr
          klow = ifeld(i,j)    ! (*)
          do k = 1, nplevel-klow+1


          !------------------------------------------------------------------!
          ! four exponentials for band 6                                     !
          !------------------------------------------------------------------!

          if (ibn .eq. 6) then

             xc=dch4(i,k)*(r1+(1.7007e-2_kreal+1.5826e-4_kreal*dt250(i,k))*dt250(i,k))
             ch4exp(i,k,1)=exp(-xc*5.80708e-3_kreal)
         
          else
             !---------------------------------------------------------------!
             ! four exponentials by powers of 12 for band 7                  !
             !---------------------------------------------------------------!

             xc=dch4(i,k)*(pa(i,k)/500.0_kreal)**0.65_kreal            &
                      *(r1+(5.9590e-4_kreal-2.2931e-6_kreal*dt250(i,k))*dt250(i,k))
             ch4exp(i,k,1)=exp(-xc*6.29247e-2_kreal)

             xc=ch4exp(i,k,1)*ch4exp(i,k,1)*ch4exp(i,k,1)
             xc=xc*xc
             ch4exp(i,k,2)=xc*xc

             xc=ch4exp(i,k,2)*ch4exp(i,k,2)*ch4exp(i,k,2)
             xc=xc*xc
             ch4exp(i,k,3)=xc*xc

             xc=ch4exp(i,k,3)*ch4exp(i,k,3)*ch4exp(i,k,3)
             xc=xc*xc
             ch4exp(i,k,4)=xc*xc

          endif

       enddo
    enddo

  end subroutine ch4exps

!=============================================================================!

  subroutine comexps(ibn,dcom,dt250,comexp,j)
    
    !------------------------------------------------------------------------!
    ! Compute co2-minor exponentials for individual layers using Eqs. (8.21) !
    ! and (8.22).                                                            !
    !                                                                        !
    ! Output parameters : comexp =  6 exponentials for each layer            !
    !------------------------------------------------------------------------!

    use phys_constants,only: r1
    use atham_module,  only: nxl, nxr, nx, ifeld
      
    integer(kint), intent(in) :: ibn,j
    real(kreal), intent(in) ::  dcom(nx,nplevel-1), dt250(nx,nplevel-1)

    real(kreal), intent(inout) :: comexp(nx,nplevel-1,6)
    
    integer(kint) :: i, k, ik, klow 
    real(kreal) :: xc

    !--------------------------------------------------------------------!
    ! Scaling and absorption data are given in Table 6                   !
    !--------------------------------------------------------------------!

       do i=nxl,nxr
          klow = ifeld(i,j)    ! (*)
          do k = 1, nplevel-klow+1

          if (ibn .eq. 4) then
             xc=dcom(i,k)*(r1+(3.5775e-2_kreal + 4.0447e-4_kreal * dt250(i,k))*dt250(i,k))
          endif

          if (ibn.eq.5) then
             xc=dcom(i,k)*(r1+(3.4268e-2_kreal + 3.7401e-4_kreal * dt250(i,k))*dt250(i,k))
          endif

          comexp(i,k,1)=exp(-xc*1.922e-7_kreal)

          do ik=2,6
             xc=comexp(i,k,ik-1)*comexp(i,k,ik-1)
             xc=xc*xc
             comexp(i,k,ik)=xc*comexp(i,k,ik-1)
          enddo

       enddo
    enddo

  end subroutine comexps

!==========================================================================!

  subroutine cfcexps(ibn,a1,b1,fk1,a2,b2,fk2,dcfc,dt250,cfcexp,j)

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

    use phys_constants,only: r1
    use atham_module,  only: nxl, nxr, nx, ifeld
      
    integer(kint), intent(in) :: ibn, j
    real(kreal), intent(in) ::  dcfc(nx,nplevel-1), dt250(nx,nplevel-1)
    real(kreal), intent(in) :: a1, b1, fk1, a2, b2, fk2

    real(kreal), intent(inout) :: cfcexp(nx,nplevel-1)
    
    integer(kint) :: i, k, klow
    real(kreal) :: xf

       do i=nxl,nxr  
          klow = ifeld(i,j)  ! (*)
          do k = 1, nplevel-klow+1

          !--------------------------------------------------------------------------!
          ! compute the scaled cfc amount (xf) and exponential (cfcexp)              !
          !--------------------------------------------------------------------------!

          if (ibn .eq. 4) then
             xf=dcfc(i,k)*(r1+(a1+b1*dt250(i,k))*dt250(i,k))
             cfcexp(i,k)=exp(-xf*fk1)
          else
             xf=dcfc(i,k)*(r1+(a2+b2*dt250(i,k))*dt250(i,k))
             cfcexp(i,k)=exp(-xf*fk2)
          endif

       enddo
    enddo

  end subroutine cfcexps

!=======================================================================================!

  subroutine b10exps(dh2o,dcont,dco2,dn2o,pa,dt250,h2oexp,conexp,co2exp,n2oexp, h2Oon,  &
                     co2on, n2Oon,j)

    !------------------------------------------------------------------------!
    ! Compute band3a exponentials for individual layers                      !
    !                                                                        !
    ! Output parameters: h2oexp,conexp,co2exp,n2oexp = exponentials for each !
    !layer                                                                   !
    !------------------------------------------------------------------------!
    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx, ifeld

    integer(kint), intent(in) :: j   
    real(kreal), intent(in), dimension(nx, nplevel-1) ::  dh2o, dcont, dco2, dn2o
    real(kreal), intent(in), dimension(nx, nplevel-1) :: pa, dt250
    real(kreal), intent(inout) :: h2oexp(nx,nplevel-1,6), conexp(nx,nplevel-1,3),  &
                                  co2exp(nx, nplevel-1,6,2), n2oexp(nx, nplevel-1,4)
    logical, intent(in) :: h2Oon, co2on, n2Oon

    integer(kint) :: i, k, klow
    real(kreal) :: xx, xx1, xx2, xx3
     
       do i=nxl,nxr
          klow = ifeld(i,j)  ! (*)
          do k = 1, nplevel-klow+1
          
          !-------------------------------------------------------------------!
          ! Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3). !
          !-------------------------------------------------------------------!
          if (h2Oon) then

             xx=dh2o(i,k)*(pa(i,k)/500.0_kreal)         &
                  * (r1+(0.0149_kreal+6.20e-5_kreal*dt250(i,k))*dt250(i,k))

             !--------------------------------------------------------------!
             ! six exponentials by powers of 8                              !
             !--------------------------------------------------------------!
             
             h2oexp(i,k,1)=exp(-xx*0.10624_kreal)

             xx=h2oexp(i,k,1)*h2oexp(i,k,1)
             xx=xx*xx
             h2oexp(i,k,2)=xx*xx

             xx=h2oexp(i,k,2)*h2oexp(i,k,2)
             xx=xx*xx
             h2oexp(i,k,3)=xx*xx

             xx=h2oexp(i,k,3)*h2oexp(i,k,3)
             xx=xx*xx
             h2oexp(i,k,4)=xx*xx

             xx=h2oexp(i,k,4)*h2oexp(i,k,4)
             xx=xx*xx
             h2oexp(i,k,5)=xx*xx

             !----------------------------------------------------------------!
             ! One exponential of h2o continuum for sub-band 3a (Table 9).    !
             !----------------------------------------------------------------!
       
             conexp(i,k,1)=exp(-dcont(i,k)*109.0_kreal)

          endif


          !------------------------------------------------------------------!
          ! Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).     !
          !------------------------------------------------------------------!
          if (co2on) then

             xx=dco2(i,k)*(pa(i,k)/300.0_kreal)**r1h       &
                  *(r1+(0.0179_kreal+1.02e-4_kreal*dt250(i,k))*dt250(i,k))

             !---------------------------------------------------------------!
             ! six exponentials by powers of 8                               !
             !---------------------------------------------------------------!

             co2exp(i,k,1,1)=exp(-xx*2.656e-5_kreal)

             xx=co2exp(i,k,1,1)*co2exp(i,k,1,1)
             xx=xx*xx
             co2exp(i,k,2,1)=xx*xx

             xx=co2exp(i,k,2,1)*co2exp(i,k,2,1)
             xx=xx*xx
             co2exp(i,k,3,1)=xx*xx

             xx=co2exp(i,k,3,1)*co2exp(i,k,3,1)
             xx=xx*xx
             co2exp(i,k,4,1)=xx*xx

             xx=co2exp(i,k,4,1)*co2exp(i,k,4,1)
             xx=xx*xx
             co2exp(i,k,5,1)=xx*xx

             xx=co2exp(i,k,5,1)*co2exp(i,k,5,1)
             xx=xx*xx
             co2exp(i,k,6,1)=xx*xx
          endif

             !----------------------------------------------------------------!
             ! Compute the scaled n2o amount for Band 10 (Table 5).           !
             !----------------------------------------------------------------!

          if (n2Oon) then
             xx=dn2o(i,k)*(r1+(1.4476e-3_kreal + 3.6656e-6_kreal * dt250(i,k))*dt250(i,k))

             !---------------------------------------------------------------!
             ! Two exponentials by powers of 58                              !
             !---------------------------------------------------------------!

             n2oexp(i,k,1)=exp(-xx*0.25238_kreal)

             xx=n2oexp(i,k,1)*n2oexp(i,k,1)
             xx1=xx*xx
             xx1=xx1*xx1
             xx2=xx1*xx1
             xx3=xx2*xx2
             n2oexp(i,k,2)=xx*xx1*xx2*xx3
             
          endif
       enddo
    enddo

  end subroutine b10exps

!=============================================================================!

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
 
    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx
      
    integer(kint), intent(in) :: k2, nxp, nh, ncoef
    real(kreal), intent(in) :: w1, p1, dwe, dpe
    real(kreal), intent(in), dimension(nx, nplevel-1) :: dw, pa, dt250
    real(kreal), intent(in), dimension(nxp,ncoef) ::  coef1, coef2, coef3

    real(kreal), intent(inout) :: x1(nx), x2(nx), x3(nx), trant(nx, nplevel)
    
    integer(kint) :: i, iw, ip
    real(kreal) :: effp, efft, xx, we, pe, fw, fp, pra, prb, prc, ax,  &
                   ba, bb, t1, ca, cb, t2

    !-------------------------------------------------------------------------!
    ! Compute effective pressure (x2)effp and temperature (x3)efft following  !
    ! Eqs. (8.28) and (8.29)                                                  !
    !-------------------------------------------------------------------------!

    do i=nxl,nxr

       x1(i)=x1(i)+dw(i,k2-1)
       x2(i)=x2(i)+pa(i,k2-1)*dw(i,k2-1)
       x3(i)=x3(i)+dt250(i,k2-1)*dw(i,k2-1)

       xx=x1(i)
       effp=x2(i)/x1(i)
       efft=x3(i)/x1(i)

       !---------------------------------------------------------------------!
       ! normalize we and pe                                                 !
       !---------------------------------------------------------------------!
       
       we=(log10(xx)-w1)/dwe
       pe=(log10(effp)-p1)/dpe

       !---------------------------------------------------------------------!
       ! Restrict the magnitudes of the normalized we and pe.                !
       !---------------------------------------------------------------------!

       we=min(we,real(nh-1_kint,kreal))
       pe=min(pe,real(nxp-1_kint,kreal))

       !-----------------------------------------------------------------------!
       ! assign iw and ip and compute the distance of we and pe from iw and ip.!
       !-----------------------------------------------------------------------!

       iw=int(we+r1)
       iw=min(iw,nh-1_kint)
       iw=max(iw, 2_kint)
       fw=we-real(iw-1_kint,kreal)

       ip=int(pe+r1)
       ip=min(ip,nxp-1_kint)
       ip=max(ip, 1_kint)
       fp=pe-real(ip-1_kint,kreal)

       !---------------------------------------------------------------------!
       ! linear interpolation in pressure                                    !
       !---------------------------------------------------------------------!

       pra = coef1(ip,iw-1)*(r1-fp)+coef1(ip+1,iw-1)*fp
       prb = coef1(ip,  iw)*(r1-fp)+coef1(ip+1,  iw)*fp
       prc = coef1(ip,iw+1)*(r1-fp)+coef1(ip+1,iw+1)*fp

       !---------------------------------------------------------------------!
       ! quadratic interpolation in absorber amount for coef1                !
       !---------------------------------------------------------------------!

       ax = (-pra*(r1-fw)+prc*(r1+fw)) *fw*r1h + prb*(r1-fw*fw)

       !---------------------------------------------------------------------!
       ! linear interpolation in absorber amount for coef2 and coef3         !
       !---------------------------------------------------------------------!

       ba = coef2(ip,  iw)*(r1-fp)+coef2(ip+1,  iw)*fp
       bb = coef2(ip,iw+1)*(r1-fp)+coef2(ip+1,iw+1)*fp
       t1 = ba*(r1-fw) + bb*fw

       ca = coef3(ip,  iw)*(r1-fp)+coef3(ip+1,  iw)*fp
       cb = coef3(ip,iw+1)*(r1-fp)+coef3(ip+1,iw+1)*fp
       t2 = ca*(r1-fw) + cb*fw

       !---------------------------------------------------------------------!
       ! update the total transmittance between levels k1 and k2             !
       !---------------------------------------------------------------------!

       trant(i,k2)= (ax + (t1 + (t2*effp)) * efft) *trant(i,k2)
        
       !---------------------------------------------------------------------!
       ! Max/Min transmission values? commented out in original              !
       !       trant(i,k2)=min(trant(i,k2),0.9999999)                        !
       !       trant(i,k2)=max(trant(i,k2),0.0000001)                          
       !---------------------------------------------------------------------!

    enddo

  end subroutine tablup

!=============================================================================================!

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
    
    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx
      
    integer(kint), intent(in) :: ibn, k, ne
    real(kreal), intent(in) :: conexp(nx, nplevel-1,3), h2oexp(nx, nplevel-1,6) 
     
    real(kreal), intent(inout) :: th2o(nx,6), tcon(nx,3), trant(nx, nplevel)
    
    integer(kint) :: i
    real(kreal) :: trnth2o

    !-------------------------------------------------------------------------------!
    ! tco2 are the six exp factors between levels k1 and k2                         !
    ! tran is the updated total transmittance between levels k1 and k2              !
    ! h2o are 6 exp factors between levels k1 & k2 from h2o line absorption         !
    ! tcon are 3 exp factors between levels k1 & k2 from h2o continuum absorption   !
    ! trnth2o is the total transmittance between levels k1 and k2 due to both line  !
    !        and continuum absorption.                                              !
    !-------------------------------------------------------------------------------!

    !------------------------------------------------------------------------!
    ! Compute th2o following Eq. (8.23).                                     !
    !------------------------------------------------------------------------!

    do i=nxl,nxr
       th2o(i,1) = th2o(i,1)*h2oexp(i,k,1)
       th2o(i,2) = th2o(i,2)*h2oexp(i,k,2)
       th2o(i,3) = th2o(i,3)*h2oexp(i,k,3)
       th2o(i,4) = th2o(i,4)*h2oexp(i,k,4)
       th2o(i,5) = th2o(i,5)*h2oexp(i,k,5)
       th2o(i,6) = th2o(i,6)*h2oexp(i,k,6)
    enddo


    if (ne .eq. 0) then

       !----------------------------------------------------------------------!
       ! Compute trnh2o following Eq. (8.25). fkw is given in Table 4.        !
       !----------------------------------------------------------------------!

       do i=nxl,nxr

          trnth2o      =(fkw(1,ibn)*th2o(i,1)         &
                          + fkw(2,ibn)*th2o(i,2)      &
                          + fkw(3,ibn)*th2o(i,3)      &
                          + fkw(4,ibn)*th2o(i,4)      &
                          + fkw(5,ibn)*th2o(i,5)      &
                          + fkw(6,ibn)*th2o(i,6))

          trant(i,k+1)=trant(i,k+1)*trnth2o

       enddo

    elseif (ne .eq. 1) then

       !-------------------------------------------------------------------!
       ! Compute trnh2o following Eqs. (8.25) and (4.27).                  !
       !-------------------------------------------------------------------!

       do i=nxl,nxr

          tcon(i,1)= tcon(i,1)*conexp(i,k,1)

          trnth2o      =(fkw(1,ibn)*th2o(i,1)         &
                          + fkw(2,ibn)*th2o(i,2)      &
                          + fkw(3,ibn)*th2o(i,3)      &
                          + fkw(4,ibn)*th2o(i,4)      &
                          + fkw(5,ibn)*th2o(i,5)      &
                          + fkw(6,ibn)*th2o(i,6))*tcon(i,1)

          trant(i,k+1)=trant(i,k+1)*trnth2o

       enddo

    else

       !------------------------------------------------------------------!
       ! For band 3. This band is divided into 3 subbands.                !
       !------------------------------------------------------------------!

       do i=nxl,nxr

          tcon(i,1)= tcon(i,1)*conexp(i,k,1)
          tcon(i,2)= tcon(i,2)*conexp(i,k,2)
          tcon(i,3)= tcon(i,3)*conexp(i,k,3)
          
          !----------------------------------------------------------------!
          ! Compute trnh2o following Eqs. (4.29) and (8.25).                    
          !----------------------------------------------------------------!

          trnth2o      = (  gkw(1,1)*th2o(i,1)                &
                             + gkw(2,1)*th2o(i,2)             &
                             + gkw(3,1)*th2o(i,3)             &
                             + gkw(4,1)*th2o(i,4)             &
                             + gkw(5,1)*th2o(i,5)             &
                             + gkw(6,1)*th2o(i,6))*tcon(i,1)  &
                             + (  gkw(1,2)*th2o(i,1)          &
                             + gkw(2,2)*th2o(i,2)             &
                             + gkw(3,2)*th2o(i,3)             &
                             + gkw(4,2)*th2o(i,4)             &
                             + gkw(5,2)*th2o(i,5)             &
                             + gkw(6,2)*th2o(i,6))*tcon(i,2)  &
                             + (  gkw(1,3)*th2o(i,1)          &
                             + gkw(2,3)*th2o(i,2)             &
                             + gkw(3,3)*th2o(i,3)             &
                             + gkw(4,3)*th2o(i,4)             &
                             + gkw(5,3)*th2o(i,5)             &
                             + gkw(6,3)*th2o(i,6))*tcon(i,3)

          trant(i,k+1)=trant(i,k+1)*trnth2o

       enddo

    endif


  end subroutine h2okdis

!===========================================================================!

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

    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx
        
    integer(kint), intent(in) :: k
    real(kreal), intent(in) :: co2exp(nx, nplevel-1,6,2)
     
    real(kreal), intent(inout) :: tco2(nx,6,2), trant(nx, nplevel)
    
    integer(kint) :: i
    real(kreal) :: xc


    !------------------------------------------------------------------------------!
    ! tco2 is the 6 exp factors between levels k1 & k2 computed from Eqs.(8.23) &  !
    ! (8.25), (See Eq.(4.30)). k-distribution functions are given in Table 10.     !
    !------------------------------------------------------------------------------!

    do i=nxl,nxr

       tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
       xc=   0.1395_kreal *tco2(i,1,1)

       tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
       xc=xc+0.1407_kreal *tco2(i,2,1)

       tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
       xc=xc+0.1549_kreal *tco2(i,3,1)

       tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
       xc=xc+0.1357_kreal *tco2(i,4,1)

       tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
       xc=xc+0.0182_kreal *tco2(i,5,1)

       tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
       xc=xc+0.0220_kreal *tco2(i,6,1)

       !-----------------------------------------------------------------------------!
       ! Band-center region                                                          !
       !-----------------------------------------------------------------------------!

       tco2(i,1,2)=tco2(i,1,2)*co2exp(i,k,1,2)
       xc=xc+0.0766_kreal *tco2(i,1,2)

       tco2(i,2,2)=tco2(i,2,2)*co2exp(i,k,2,2)
       xc=xc+0.1372_kreal *tco2(i,2,2)

       tco2(i,3,2)=tco2(i,3,2)*co2exp(i,k,3,2)
       xc=xc+0.1189_kreal *tco2(i,3,2)

       tco2(i,4,2)=tco2(i,4,2)*co2exp(i,k,4,2)
       xc=xc+0.0335_kreal *tco2(i,4,2)

       tco2(i,5,2)=tco2(i,5,2)*co2exp(i,k,5,2)
       xc=xc+0.0169_kreal *tco2(i,5,2)

       tco2(i,6,2)=tco2(i,6,2)*co2exp(i,k,6,2)
       xc=xc+0.0059_kreal *tco2(i,6,2)

       trant(i,k+1)=trant(i,k+1)*xc

    enddo

  end subroutine co2kdis

!===========================================================================================!

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

    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx
      
    integer(kint), intent(in) :: ibn, k
    real(kreal), intent(in) :: n2oexp(nx, nplevel-1,4)
     
    real(kreal), intent(inout) :: tn2o(nx,4), trant(nx, nplevel)
    
    integer(kint) :: i
    real(kreal) :: xc


    !-----------------------------------------------------------------------------!
    ! tn2o is computed from Eq. (8.23).                                           !
    ! xc is the total n2o transmittance computed from (8.25)                      !
    ! The k-distribution functions are given in Table 5.                          !
    !-----------------------------------------------------------------------------!

    do i=nxl,nxr

       !------------------------------------------------------------------------------!
       ! band 6                                                                       !
       !------------------------------------------------------------------------------!

       if (ibn .eq. 6) then

          tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
          xc= 0.940414_kreal*tn2o(i,1)

          tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
          xc=xc + 0.059586_kreal*tn2o(i,2)

       else
          
          !--------------------------------------------------------------------------------!
          ! band 7                                                                         !
          !--------------------------------------------------------------------------------!

          tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
          xc= 0.561961_kreal*tn2o(i,1)
          
          tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
          xc=xc + 0.138707_kreal*tn2o(i,2)

          tn2o(i,3)=tn2o(i,3)*n2oexp(i,k,3)
          xc=xc + 0.240670_kreal*tn2o(i,3)

          tn2o(i,4)=tn2o(i,4)*n2oexp(i,k,4)
          xc=xc + 0.058662_kreal*tn2o(i,4)
           
       endif

       trant(i,k+1)=trant(i,k+1)*xc
       
    enddo

  end subroutine n2okdis
   
!============================================================================================!

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

    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx
      
    integer(kint), intent(in) :: ibn, k
    real(kreal), intent(in) :: ch4exp(nx, nplevel-1,4)
     
    real(kreal), intent(inout) :: tch4(nx,4), trant(nx, nplevel)
    
    integer(kint) :: i
    real(kreal) :: xc

    !-----------------------------------------------------------------------------!
    ! tch4 is computed from Eq. (8.23). xc is total ch4 transmittance computed    !
    ! from (8.25). The k-distribution functions are given in Table 5.             !
    !-----------------------------------------------------------------------------!

    do i=nxl,nxr

       !-----------------------------------------------------------------------------!
       ! band 6                                                                      !
       !-----------------------------------------------------------------------------!

       if (ibn .eq. 6) then

          tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
          xc= tch4(i,1)

       else
          
          !--------------------------------------------------------------------------------!
          ! band 7                                                                         !
          !--------------------------------------------------------------------------------!
          
          tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
          xc= 0.610650_kreal*tch4(i,1)
          
          tch4(i,2)=tch4(i,2)*ch4exp(i,k,2)
          xc=xc + 0.280212_kreal*tch4(i,2)

          tch4(i,3)=tch4(i,3)*ch4exp(i,k,3)
          xc=xc + 0.107349_kreal*tch4(i,3)

          tch4(i,4)=tch4(i,4)*ch4exp(i,k,4)
          xc=xc + 0.001789_kreal*tch4(i,4)

       endif

       trant(i,k+1)=trant(i,k+1)*xc

    enddo

  end subroutine ch4kdis

!============================================================================================!


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

    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx
      
    integer(kint), intent(in) :: ibn, k
    real(kreal), intent(in) :: comexp(nx, nplevel-1,6)
     
    real(kreal), intent(inout) :: tcom(nx,6), trant(nx, nplevel)
    
    integer(kint) :: i
    real(kreal) :: xc

  
    !-----------------------------------------------------------------------------!
    ! tcom is computed from Eq. (8.23). xc is total co2 transmittance computed    !
    ! from (8.25). The k-distribution functions are given in Table 6.             !
    !-----------------------------------------------------------------------------!

    do i=nxl,nxr

       !-----------------------------------------------------------------------------!
       ! band 4                                                                      !
       !-----------------------------------------------------------------------------!

       if (ibn .eq. 4) then

          tcom(i,1)=tcom(i,1)*comexp(i,k,1)
          xc= 0.12159_kreal*tcom(i,1)
          tcom(i,2)=tcom(i,2)*comexp(i,k,2)
          xc=xc + 0.24359_kreal*tcom(i,2)
          tcom(i,3)=tcom(i,3)*comexp(i,k,3)
          xc=xc + 0.24981_kreal*tcom(i,3)
          tcom(i,4)=tcom(i,4)*comexp(i,k,4)
          xc=xc + 0.26427_kreal*tcom(i,4)
          tcom(i,5)=tcom(i,5)*comexp(i,k,5)
          xc=xc + 0.07807_kreal*tcom(i,5)
          tcom(i,6)=tcom(i,6)*comexp(i,k,6)
          xc=xc + 0.04267_kreal*tcom(i,6)
        
       else
          
          !--------------------------------------------------------------------------------!
          ! band 5                                                                         !
          !--------------------------------------------------------------------------------!

          tcom(i,1)=tcom(i,1)*comexp(i,k,1)
          xc= 0.06869_kreal*tcom(i,1)
          tcom(i,2)=tcom(i,2)*comexp(i,k,2)
          xc=xc + 0.14795_kreal*tcom(i,2)
          tcom(i,3)=tcom(i,3)*comexp(i,k,3)
          xc=xc + 0.19512_kreal*tcom(i,3)
          tcom(i,4)=tcom(i,4)*comexp(i,k,4)
          xc=xc + 0.33446_kreal*tcom(i,4)
          tcom(i,5)=tcom(i,5)*comexp(i,k,5)
          xc=xc + 0.17199_kreal*tcom(i,5)
          tcom(i,6)=tcom(i,6)*comexp(i,k,6)
          xc=xc + 0.08179_kreal*tcom(i,6)
       endif

       trant(i,k+1)=trant(i,k+1)*xc
    enddo

  end subroutine comkdis

!===============================================================================================!

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

    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx
      
    integer(kint), intent(in) :: k
    real(kreal), intent(in) :: cfcexp(nx, nplevel-1)
     
    real(kreal), intent(inout) :: tcfc(nx), trant(nx, nplevel)
    
    integer(kint) :: i
    
    do i=nxl,nxr

       tcfc(i)=tcfc(i)*cfcexp(i,k)
       trant(i,k+1)=trant(i,k+1)*tcfc(i)
  
    enddo
  

  end subroutine cfckdis

!===============================================================================================!

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

    use phys_constants,only: r1, r1h
    use atham_module,  only: nxl, nxr, nx

    integer(kint), intent(in) :: k
    real(kreal), intent(in) :: h2oexp(nx,nplevel-1,6), conexp(nx,nplevel-1,3),      &
                                  co2exp(nx, nplevel-1,6,2), n2oexp(nx, nplevel-1,4)
    real(kreal), intent(inout) :: th2o(nx,6), tcon(nx,3), tco2(nx,6,2), tn2o(nx,4), &
                                    trant(nx, nplevel)    

    logical, intent(in) :: h2Oon, co2on, n2Oon

    integer(kint) :: i
    real(kreal) :: xx

    !-----------------------------------------------------------------------------!
    ! For h2o line. The k-distribution functions are given in Table 4.            !
    !-----------------------------------------------------------------------------!
    trant(nxl:nxr,k+1) = r1

    if (h2Oon) then
       do i=nxl,nxr

          th2o(i,1)=th2o(i,1)*h2oexp(i,k,1)
          xx= 0.3153_kreal*th2o(i,1)

          th2o(i,2)=th2o(i,2)*h2oexp(i,k,2)
          xx=xx + 0.4604_kreal*th2o(i,2)

          th2o(i,3)=th2o(i,3)*h2oexp(i,k,3)
          xx=xx + 0.1326_kreal*th2o(i,3)

          th2o(i,4)=th2o(i,4)*h2oexp(i,k,4)
          xx=xx + 0.0798_kreal*th2o(i,4)
 
          th2o(i,5)=th2o(i,5)*h2oexp(i,k,5)
          xx=xx + 0.0119_kreal*th2o(i,5)

          trant(i,k+1)=trant(i,k+1) * xx
       enddo
    !-----------------------------------------------------------------------------!
    ! For h2o continuum. Note that conexp(i,k,3) is for subband 3a.               !
    !-----------------------------------------------------------------------------!

       do i=nxl,nxr

          tcon(i,1)=tcon(i,1)*conexp(i,k,3)
          trant(i,k+1)=trant(i,k+1)*tcon(i,1)
       enddo

    endif
 
    !-----------------------------------------------------------------------------!
    ! For co2 (Table 6)                                                           !
    !-----------------------------------------------------------------------------!
    
    if (co2on) then
       do i=nxl,nxr

          tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
          xx= 0.2673_kreal*tco2(i,1,1)
          
          tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
          xx=xx + 0.2201_kreal*tco2(i,2,1)
       
          tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
          xx=xx + 0.2106_kreal*tco2(i,3,1)

          tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
          xx=xx + 0.2409_kreal*tco2(i,4,1)

          tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
          xx=xx + 0.0196_kreal*tco2(i,5,1)
       
          tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
          xx=xx + 0.0415_kreal*tco2(i,6,1)

          trant(i,k+1)=trant(i,k+1)*xx
       enddo
    endif

    !-----------------------------------------------------------------------------!
    ! For n2o (Table 5)                                                           !
    !-----------------------------------------------------------------------------!

    if (n2Oon) then
       do i=nxl,nxr

          tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
          xx= 0.970831_kreal*tn2o(i,1)

          tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
          xx=xx + 0.029169_kreal*tn2o(i,2)
          trant(i,k+1)=trant(i,k+1)*(xx-r1)
       enddo
    endif
       
  end subroutine b10kdis

!-----------------------------------------------------------------------------------
 
end module lw_radiation_nasa
