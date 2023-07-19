!>  \file sfc_diff.f
!! This file contains the surface roughness length formulation based on
!! the surface sublayer scheme in
!! Zeng and Dickinson (1998) \cite zeng_and_dickinson_1998.

!> This module contains the CCPP-compliant GFS surface layer scheme.
      module sfc_diff

      use machine , only : kind_phys

      implicit none

      public :: sfc_diff_run
      public :: stability

      private

      real (kind=kind_phys), parameter :: ca=0.4_kind_phys  ! ca - von karman constant

      contains

!> \defgroup GFS_diff_main GFS Surface Layer Module
!> This module calculates surface roughness length.
!> @{
!! This subroutine includes the surface roughness length formulation
!! based on the surface sublayer scheme in
!! Zeng and Dickinson (1998) \cite zeng_and_dickinson_1998.
!> \section arg_table_sfc_diff_run Argument Table
!! \htmlinclude sfc_diff_run.html
!!
!>  \section general_diff GFS Surface Layer Scheme General Algorithm
!! - Calculate the thermal roughness length formulation over the ocean (see eq. (25) and (26)
!!  in Zeng et al. (1998) \cite zeng_et_al_1998).
!! - Calculate Zeng's momentum roughness length formulation over land and sea ice.
!! - Calculate the new vegetation-dependent formulation of thermal roughness length
!! (Zheng et al.(2012) \cite zheng_et_al_2012).
!! Zheng et al. (2012) \cite zheng_et_al_2012 proposed a new formulation on
!! \f$ln(Z_{0m}^,/Z_{0t})\f$ as follows:
!! \f[
!!  ln(Z_{0m}^,/Z_{0t})=(1-GVF)^2C_{zil}k(u*Z_{0g}/\nu)^{0.5}
!! \f]
!! where \f$Z_{0m}^,\f$ is the effective momentum roughness length
!! computed in the following equation for each grid, \f$Z_{0t}\f$
!! is the roughness lenghth for heat, \f$C_{zil}\f$ is a coefficient
!! (taken as 0.8), k is the Von Karman constant (0.4),
!! \f$\nu=1.5\times10^{-5}m^{2}s^{-1}\f$ is the molecular viscosity,
!! \f$u*\f$ is the friction velocity, and \f$Z_{0g}\f$ is the bare
!! soil roughness length for momentum (taken as 0.01).
!! \n In order to consider the convergence of \f$Z_{0m}\f$ between
!! fully vegetated and bare soil, the effective \f$Z_{0m}^,\f$ is
!! computed:
!! \f[
!!  \ln(Z_{0m}^,)=(1-GVF)^{2}\ln(Z_{0g})+\left[1-(1-GVF)^{2}\right]\ln(Z_{0m})
!!\f]
!! - Calculate the exchange coefficients:\f$cm\f$, \f$ch\f$, and \f$stress\f$ as inputs of other \a sfc schemes.
!!
      subroutine sfc_diff_run (im,rvrdm1,eps,epsm1,grav,                &  !intent(in)
     &                    ps,t1,q1,z1,garea,wind,                       &  !intent(in)
     &                    prsl1,prslki,prsik1,prslk1,                   &  !intent(in)
     &                    sigmaf,vegtype,shdmax,ivegsrc,                &  !intent(in)
     &                    z0pert,ztpert,                                &  ! mg, sfc-perts !intent(in)
     &                    flag_iter,redrag,                             &  !intent(in)
     &                    u10m,v10m,sfc_z0_type,                        &  !hafs,z0 type !intent(in)
     &                    wet,dry,icy,                                  &  !intent(in)
     &                    thsfc_loc,                                    &  !intent(in)
     &                    tskin_wat, tskin_lnd, tskin_ice,              &  !intent(in)
     &                    tsurf_wat, tsurf_lnd, tsurf_ice,              &  !intent(in)
     &                     z0rl_wat,  z0rl_lnd,  z0rl_ice,              &  !intent(inout)
     &                     z0rl_wav,                                    &  !intent(inout)
     &                    ustar_wat, ustar_lnd, ustar_ice,              &  !intent(inout)
     &                       cm_wat,    cm_lnd,    cm_ice,              &  !intent(inout)
     &                       ch_wat,    ch_lnd,    ch_ice,              &  !intent(inout)
     &                       rb_wat,    rb_lnd,    rb_ice,              &  !intent(inout)
     &                   stress_wat,stress_lnd,stress_ice,              &  !intent(inout)
     &                       fm_wat,    fm_lnd,    fm_ice,              &  !intent(inout)
     &                       fh_wat,    fh_lnd,    fh_ice,              &  !intent(inout)
     &                     fm10_wat,  fm10_lnd,  fm10_ice,              &  !intent(inout)
     &                      fh2_wat,   fh2_lnd,   fh2_ice,              &  !intent(inout)
     &                    ztmax_wat, ztmax_lnd, ztmax_ice,              &  !intent(inout)
     &                    zvfun,                                        &  !intent(out)
     &                    errmsg, errflg)                                  !intent(out)
!
      implicit none
!
      integer, parameter  :: kp = kind_phys
      integer, intent(in) :: im, ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean

      integer, dimension(:), intent(in) :: vegtype

      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      logical, dimension(:), intent(in) :: flag_iter, dry, icy
      logical, dimension(:), intent(inout) :: wet

      logical, intent(in) :: thsfc_loc ! Flag for reference pressure in theta calculation

      real(kind=kind_phys), dimension(:), intent(in)    :: u10m,v10m
      real(kind=kind_phys), intent(in) :: rvrdm1, eps, epsm1, grav
      real(kind=kind_phys), dimension(:), intent(in)    ::              &
     &                    ps,t1,q1,z1,garea,prsl1,prslki,prsik1,prslk1, &
     &                    wind,sigmaf,shdmax,                           &
     &                    z0pert,ztpert ! mg, sfc-perts
      real(kind=kind_phys), dimension(:), intent(in)    ::              &
     &                    tskin_wat, tskin_lnd, tskin_ice,              &
     &                    tsurf_wat, tsurf_lnd, tsurf_ice

      real(kind=kind_phys), dimension(:), intent(in)    :: z0rl_wav
      real(kind=kind_phys), dimension(:), intent(inout) ::              &
     &                     z0rl_wat,  z0rl_lnd,  z0rl_ice,              &
     &                    ustar_wat, ustar_lnd, ustar_ice,              &
     &                       cm_wat,    cm_lnd,    cm_ice,              &
     &                       ch_wat,    ch_lnd,    ch_ice,              &
     &                       rb_wat,    rb_lnd,    rb_ice,              &
     &                   stress_wat,stress_lnd,stress_ice,              &
     &                       fm_wat,    fm_lnd,    fm_ice,              &
     &                       fh_wat,    fh_lnd,    fh_ice,              &
     &                     fm10_wat,  fm10_lnd,  fm10_ice,              &
     &                      fh2_wat,   fh2_lnd,   fh2_ice,              &
     &                    ztmax_wat, ztmax_lnd, ztmax_ice
      real(kind=kind_phys), dimension(:), intent(out)    :: zvfun
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!     locals
!
      integer   i
!
      real(kind=kind_phys) :: rat, tv1, thv1, restar, wind10m,
     &                        czilc, tem1, tem2, virtfac
!

      real(kind=kind_phys) :: tvs, z0, z0max, ztmax, gdx
!
      real(kind=kind_phys), parameter :: z0lo=0.1, z0up=1.0
!
      real(kind=kind_phys), parameter ::
     &        one=1.0_kp, zero=0.0_kp, half=0.5_kp, qmin=1.0e-8_kp
     &,       charnock=.018_kp, z0s_max=.317e-2_kp                      &! a limiting value at high winds over sea
     &,       zmin=1.0e-6_kp                                            &
     &,       vis=1.4e-5_kp, rnu=1.51e-5_kp, visi=one/vis               &
     &,       log01=log(0.01_kp), log05=log(0.05_kp), log07=log(0.07_kp)

!     parameter (charnock=.014,ca=.4)!c ca is the von karman constant
!     parameter (alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041)
!     parameter (a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5)

!     real(kind=kind_phys) aa1,bb1,bb2,cc,cc1,cc2,arnu
!     parameter (aa1=-1.076,bb1=.7045,cc1=-.05808)
!     parameter (bb2=-.1954,cc2=.009999)
!     parameter (arnu=.135*rnu)
!
!    z0s_max=.196e-2 for u10_crit=25 m/s
!    z0s_max=.317e-2 for u10_crit=30 m/s
!    z0s_max=.479e-2 for u10_crit=35 m/s
!
! mbek -- toga-coare flux algorithm
!     parameter (rnu=1.51e-5,arnu=0.11*rnu)

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals, wind is wind speed,
!  surface roughness length is converted to m from cm
!
!       write(0,*)'in sfc_diff, sfc_z0_type=',sfc_z0_type

      do i=1,im
        if(flag_iter(i)) then

          ! Need to initialize ztmax arrays
          ztmax_lnd(i) = 1. ! log(1) = 0
          ztmax_ice(i) = 1. ! log(1) = 0
          ztmax_wat(i) = 1. ! log(1) = 0

          virtfac = one + rvrdm1 * max(q1(i),qmin)

          tv1 = t1(i) * virtfac ! Virtual temperature in middle of lowest layer
          if(thsfc_loc) then ! Use local potential temperature
            thv1  = t1(i) * prslki(i) * virtfac
          else ! Use potential temperature reference to 1000 hPa
            thv1    = t1(i) / prslk1(i) * virtfac
          endif

          zvfun(i) = zero
          gdx = sqrt(garea(i))

!  compute stability dependent exchange coefficients
!  this portion of the code is presently suppressed
!
          if (dry(i)) then ! Some land

            if(thsfc_loc) then ! Use local potential temperature
              tvs   = half * (tsurf_lnd(i)+tskin_lnd(i)) * virtfac
            else ! Use potential temperature referenced to 1000 hPa
              tvs   = half * (tsurf_lnd(i)+tskin_lnd(i))/prsik1(i)
     &                     * virtfac
            endif

            z0max = max(zmin, min(0.01_kp * z0rl_lnd(i), z1(i)))
!** xubin's new z0  over land
            tem1  = one - shdmax(i)
            tem2  = tem1 * tem1
            tem1  = one  - tem2

            if( ivegsrc == 1 ) then

              if (vegtype(i) == 10) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 6) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 7) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01_kp
              elseif (vegtype(i) == 16) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01_kp
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            elseif (ivegsrc == 2 ) then

              if (vegtype(i) == 7) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 8) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 9) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01_kp
              elseif (vegtype(i) == 11) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01_kp
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            endif
! mg, sfc-perts: add surface perturbations to z0max over land
            if (z0pert(i) /= zero ) then
              z0max = z0max * (10.0_kp**z0pert(i))
            endif

            z0max = max(z0max, zmin)

!!          czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
!           czilc = 0.8_kp
!
!           tem1  = 1.0_kp - sigmaf(i)
!           ztmax_lnd(i) = z0max*exp( - tem1*tem1
!    &              * czilc*ca*sqrt(ustar_lnd(i)*(0.01/1.5e-05)))
!
            czilc = 10.0_kp ** (- 4.0_kp * z0max) ! Trier et al. (2011,WAF)
            czilc = max(min(czilc, 0.8_kp), 0.08_kp)
            tem1 = 1.0_kp - sigmaf(i)
            czilc = czilc * tem1 * tem1
            ztmax_lnd(i) = z0max * exp( - czilc * ca
     &            * 258.2_kp * sqrt(ustar_lnd(i)*z0max) )
!
! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
            if (ztpert(i) /= zero) then
              ztmax_lnd(i) = ztmax_lnd(i) * (10.0_kp**ztpert(i))
            endif
            ztmax_lnd(i) = max(ztmax_lnd(i), zmin)
!
! compute a function of surface roughness & green vegetation fraction (zvfun)       
!
            tem1 = (z0max - z0lo) / (z0up - z0lo)
            tem1 = min(max(tem1, zero), 1.0_kp)
            tem2 = max(sigmaf(i), 0.1_kp)
            zvfun(i) = sqrt(tem1 * tem2)
!
            call stability
!  ---  inputs:
     &       (z1(i), zvfun(i), gdx, tv1, thv1, wind(i),
     &        z0max, ztmax_lnd(i), tvs, grav, thsfc_loc,
!  ---  outputs:
     &        rb_lnd(i), fm_lnd(i), fh_lnd(i), fm10_lnd(i), fh2_lnd(i),
     &        cm_lnd(i), ch_lnd(i), stress_lnd(i), ustar_lnd(i))
          endif ! Dry points

          if (icy(i)) then ! Some ice

            zvfun(i) = zero

            if(thsfc_loc) then ! Use local potential temperature
              tvs   = half * (tsurf_ice(i)+tskin_ice(i)) * virtfac
            else ! Use potential temperature referenced to 1000 hPa
              tvs   = half * (tsurf_ice(i)+tskin_ice(i))/prsik1(i)
     &                     * virtfac 
            endif

            z0max = max(zmin, min(0.01_kp * z0rl_ice(i), z1(i)))
!** xubin's new z0  over land and sea ice
            tem1  = one - shdmax(i)
            tem2  = tem1 * tem1
            tem1  = one  - tem2

! Removed the following lines by W. Zheng, for effective z0m (z0max) is applied only
! for land.
!wz         if( ivegsrc == 1 ) then
!wz
!wz           z0max = exp( tem2*log01 + tem1*log(z0max) )
!wz         elseif (ivegsrc == 2 ) then
!wz           z0max = exp( tem2*log01 + tem1*log(z0max) )
!wz         endif

            z0max = max(z0max, zmin)

!!          czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height
!           dependance of czil
!           czilc = 0.8_kp
!
!           tem1  = 1.0_kp - sigmaf(i)
!           ztmax_ice(i) = z0max*exp( - tem1*tem1
!    &              * czilc*ca*sqrt(ustar_ice(i)*(0.01/1.5e-05)))
!
            czilc = 10.0_kp ** (- 4.0_kp * z0max)
            czilc = max(min(czilc, 0.8_kp), 0.08_kp)
            tem1 = 1.0_kp - sigmaf(i)
            czilc = czilc * tem1 * tem1
            ztmax_ice(i) = z0max * exp( - czilc * ca
     &            * 258.2_kp * sqrt(ustar_ice(i)*z0max) )
!
            ztmax_ice(i) = max(ztmax_ice(i), 1.0e-6)
!
            call stability
!  ---  inputs:
     &     (z1(i), zvfun(i), gdx, tv1, thv1, wind(i),
     &      z0max, ztmax_ice(i), tvs, grav, thsfc_loc,
!  ---  outputs:
     &      rb_ice(i), fm_ice(i), fh_ice(i), fm10_ice(i), fh2_ice(i),
     &      cm_ice(i), ch_ice(i), stress_ice(i), ustar_ice(i))
      endif ! Icy points

! BWG: Everything from here to end of subroutine was after
!      the stuff now put into "stability"

          if (wet(i)) then ! Some open ocean
  
            zvfun(i) = zero

            if(thsfc_loc) then ! Use local potential temperature
              tvs        = half * (tsurf_wat(i)+tskin_wat(i)) * virtfac
            else
              tvs        = half * (tsurf_wat(i)+tskin_wat(i))/prsik1(i)
     &                          * virtfac
            endif

            z0           = 0.01_kp * z0rl_wat(i)
            z0max        = max(zmin, min(z0,z1(i)))
!           ustar_wat(i) = sqrt(grav * z0 / charnock)
            wind10m      = sqrt(u10m(i)*u10m(i)+v10m(i)*v10m(i))

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar_wat(i)*z0max*visi, 0.000001_kp)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat   = min(7.0_kp, 2.67_kp * sqrt(sqrt(restar)) - 2.57_kp)
            ztmax_wat(i) = max(z0max * exp(-rat), zmin)
!
            if (sfc_z0_type == 6) then
              call znot_t_v6(wind10m, ztmax_wat(i))   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type == 7) then
              call znot_t_v7(wind10m, ztmax_wat(i))   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type > 0) then
              write(0,*)'no option for sfc_z0_type=',sfc_z0_type
              errflg = 1
              errmsg = 'ERROR(sfc_diff_run): no option for sfc_z0_type'
              return
            endif
!
            call stability
!  ---  inputs:
     &       (z1(i), zvfun(i), gdx, tv1, thv1, wind(i),
     &        z0max, ztmax_wat(i), tvs, grav, thsfc_loc,
!  ---  outputs:
     &        rb_wat(i), fm_wat(i), fh_wat(i), fm10_wat(i), fh2_wat(i),
     &        cm_wat(i), ch_wat(i), stress_wat(i), ustar_wat(i))
!
!  update z0 over ocean
!
            if (sfc_z0_type >= 0) then
              if (sfc_z0_type == 0) then
!               z0 = (charnock / grav) * ustar_wat(i) * ustar_wat(i)
                tem1 = 0.11 * vis / ustar_wat(i)
                z0 = tem1 + (charnock/grav)*ustar_wat(i)*ustar_wat(i)


! mbek -- toga-coare flux algorithm
!               z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!               cc = ustar(i) * z0 / rnu
!               pp = cc / (1. + cc)
!               ff = grav * arnu / (charnock * ustar(i) ** 3)
!               z0 = arnu / (ustar(i) * ff ** pp)

                if (redrag) then
                  z0rl_wat(i) = 100.0_kp * max(min(z0, z0s_max),        &
     &                                                 1.0e-7_kp)
                else
                  z0rl_wat(i) = 100.0_kp * max(min(z0,0.1_kp), 1.e-7_kp)
                endif

              elseif (sfc_z0_type == 6) then   ! wang
                 call znot_m_v6(wind10m, z0)   ! wind, m/s, z0, m
                 z0rl_wat(i) = 100.0_kp * z0   ! cm
              elseif (sfc_z0_type == 7) then   ! wang
                 call znot_m_v7(wind10m, z0)   ! wind, m/s, z0, m
                 z0rl_wat(i) = 100.0_kp * z0   ! cm
              else
                 z0rl_wat(i) = 1.0e-4_kp
              endif

            elseif (z0rl_wav(i) <= 1.0e-7_kp .or.                       &
     &              z0rl_wav(i) > 1.0_kp) then
!             z0 = (charnock / grav) * ustar_wat(i) * ustar_wat(i)
              tem1 = 0.11 * vis / ustar_wat(i)
              z0 = tem1 + (charnock/grav)*ustar_wat(i)*ustar_wat(i)

              if (redrag) then
                z0rl_wat(i) = 100.0_kp * max(min(z0, z0s_max),1.0e-7_kp)
              else
                z0rl_wat(i) = 100.0_kp * max(min(z0,0.1_kp), 1.0e-7_kp)
              endif
            endif

          endif              ! end of if(open ocean)
!
        endif                ! end of if(flagiter) loop
      enddo

      return
      end subroutine sfc_diff_run

!----------------------------------------
!>\ingroup GFS_diff_main
      subroutine stability                                              &
!  ---  inputs:
     &     ( z1, zvfun, gdx, tv1, thv1, wind, z0max, ztmax, tvs, grav,  &
     &       thsfc_loc,                                                 &
!  ---  outputs:
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar)
!-----

      integer, parameter :: kp = kind_phys
!  ---  inputs:
      real(kind=kind_phys), intent(in) ::                               &
     &       z1, zvfun, gdx, tv1, thv1, wind, z0max, ztmax, tvs, grav
      logical,              intent(in) :: thsfc_loc

!  ---  outputs:
      real(kind=kind_phys), intent(out) ::                              &
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar

!  ---  locals:
      real(kind=kind_phys), parameter :: alpha=5.0_kp, a0=-3.975_kp     &
     &,         a1=12.32_kp,   alpha4=4.0_kp*alpha                      &
     &,         b1=-7.755_kp,  b2=6.041_kp                              &
     &,         xkrefsqr=0.3_kp,  xkmin=0.05_kp                         & 
     &,         xkgdx=3000.0_kp                                         &
     &,         a0p=-7.941_kp, a1p=24.75_kp, b1p=-8.705_kp, b2p=7.899_kp&
     &,         zolmin=-10.0_kp, zero=0.0_kp, one=1.0_kp

      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2,
     &                     z1i,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
     &                     tem1,   tem2,   zolmax

      real(kind=kind_phys) xkzo

          z1i = one / z1

!
!  set background diffusivities with one for gdx >= xkgdx and 
!   as a function of horizontal grid size for gdx < xkgdx 
!   (i.e., gdx/xkgdx for gdx < xkgdx)
!
          if(gdx >= xkgdx) then
            xkzo = one
          else
            xkzo = gdx / xkgdx
          endif

          tem1 = tv1 - tvs
          if(tem1 > zero) then
            tem2 = xkzo * zvfun
            xkzo = min(max(tem2, xkmin), xkzo)
          endif

          zolmax = xkrefsqr / sqrt(xkzo)

!  compute stability indices (rb and hlinf)

          dtv     = thv1 - tvs
          adtv    = max(abs(dtv),0.001_kp)
          dtv     = sign(1.0_kp,dtv) * adtv

          if(thsfc_loc) then ! Use local potential temperature
            rb      = max(-5000.0_kp, (grav+grav) * dtv * z1
     &              / ((thv1 + tvs) * wind * wind))
          else ! Use potential temperature referenced to 1000 hPa
            rb      = max(-5000.0_kp, grav * dtv * z1
     &              / (tv1 * wind * wind))
          endif

          tem1    = one / z0max
          tem2    = one / ztmax
          fm      = log((z0max+z1)  * tem1)
          fh      = log((ztmax+z1)  * tem2)
          fm10    = log((z0max+10.0_kp) * tem1)
          fh2     = log((ztmax+2.0_kp)  * tem2)
          hlinf   = rb * fm * fm / fh
          hlinf   = min(max(hlinf,zolmin),zolmax)
!
!  stable case
!
          if (dtv >= zero) then
            hl1 = hlinf
            if(hlinf > 0.25_kp) then
              tem1   = hlinf * z1i
              hl0inf = z0max * tem1
              hltinf = ztmax * tem1
              aa     = sqrt(one + alpha4 * hlinf)
              aa0    = sqrt(one + alpha4 * hl0inf)
              bb     = aa
              bb0    = sqrt(one + alpha4 * hltinf)
              pm     = aa0 - aa + log( (aa + one)/(aa0 + one) )
              ph     = bb0 - bb + log( (bb + one)/(bb0 + one) )
              fms    = fm - pm
              fhs    = fh - ph
              hl1    = fms * fms * rb / fhs
              hl1    = min(hl1, zolmax)
            endif
!
!  second iteration
!
            tem1  = hl1 * z1i
            hl0   = z0max * tem1
            hlt   = ztmax * tem1
            aa    = sqrt(one + alpha4 * hl1)
            aa0   = sqrt(one + alpha4 * hl0)
            bb    = aa
            bb0   = sqrt(one + alpha4 * hlt)
            pm    = aa0 - aa + log( (one+aa)/(one+aa0) )
            ph    = bb0 - bb + log( (one+bb)/(one+bb0) )
            hl110 = hl1 * 10.0_kp * z1i
            aa    = sqrt(one + alpha4 * hl110)
            pm10  = aa0 - aa + log( (one+aa)/(one+aa0) )
            hl12  = (hl1+hl1) * z1i
!           aa    = sqrt(one + alpha4 * hl12)
            bb    = sqrt(one + alpha4 * hl12)
            ph2   = bb0 - bb + log( (one+bb)/(one+bb0) )
!
!  unstable case - check for unphysical obukhov length
!
          else                          ! dtv < 0 case
            olinf = z1 / hlinf
            tem1  = 50.0_kp * z0max
            if(abs(olinf) <= tem1) then
              hlinf = -z1 / tem1
              hlinf = max(hlinf, zolmin)
            endif
!
!  get pm and ph
!
            if (hlinf >= -0.5_kp) then
              hl1   = hlinf
              pm    = (a0  + a1*hl1)  * hl1   / (one+ (b1+b2*hl1)  *hl1)
              ph    = (a0p + a1p*hl1) * hl1   / (one+ (b1p+b2p*hl1)*hl1)
              hl110 = hl1 * 10.0_kp * z1i
              pm10  = (a0 + a1*hl110) * hl110/(one+(b1+b2*hl110)*hl110)
              hl12  = (hl1+hl1) * z1i
              ph2   = (a0p + a1p*hl12) * hl12/(one+(b1p+b2p*hl12)*hl12)
            else                       ! hlinf < 0.05
              hl1   = -hlinf
              tem1  = one / sqrt(hl1)
              pm    = log(hl1) + 2.0_kp * sqrt(tem1) - .8776_kp
              ph    = log(hl1) + 0.5_kp * tem1 + 1.386_kp
!             pm    = log(hl1) + 2.0 * hl1 ** (-.25) - .8776
!             ph    = log(hl1) + 0.5 * hl1 ** (-.5) + 1.386
              hl110 = hl1 * 10.0_kp * z1i
              pm10  = log(hl110) + 2.0_kp/sqrt(sqrt(hl110)) - 0.8776_kp
!             pm10  = log(hl110) + 2. * hl110 ** (-.25) - .8776
              hl12  = (hl1+hl1) * z1i
              ph2   = log(hl12) + 0.5_kp / sqrt(hl12) + 1.386_kp
!             ph2   = log(hl12) + .5 * hl12 ** (-.5) + 1.386
            endif

          endif          ! end of if (dtv >= 0 ) then loop
!
!  finish the exchange coefficient computation to provide fm and fh
!
          fm        = fm - pm
          fh        = fh - ph
          fm10      = fm10 - pm10
          fh2       = fh2 - ph2
          cm        = ca * ca / (fm * fm)
          ch        = ca * ca / (fm * fh)
          tem1      = 0.00001_kp/z1
          cm        = max(cm, tem1)
          ch        = max(ch, tem1)
          stress    = cm * wind * wind
          ustar     = sqrt(stress)

      return
!.................................
      end subroutine stability
!---------------------------------


!! add fitted z0,zt curves for hurricane application (used in HWRF/HMON)
!! Weiguo Wang, 2019-0425

      SUBROUTINE znot_m_v6(uref, znotm)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate areodynamical roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
! For high winds, try to fit available observational data
!
! Bin Liu, NOAA/NCEP/EMC 2017
! 
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znotm
      real(kind=kind_phys), parameter  :: p13 = -1.296521881682694e-02,
     &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,
     &      p10 = -8.396975715683501e+00,

     &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,
     &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,
     &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,

     &      p35 =  1.840430200185075e-07, p34 = -2.793849676757154e-05,
     &      p33 =  1.735308193700643e-03, p32 = -6.139315534216305e-02,
     &      p31 =  1.255457892775006e+00, p30 = -1.663993561652530e+01,

     &      p40 =  4.579369142033410e-04
  

       if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp(p10 + uref * (p11 + uref * (p12 + uref*p13))) 
       elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23
     &              + uref * (p24 + uref * p25))))
       elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33
     &                   + uref * (p34 + uref * p35)))))
       elseif ( uref > 53.0) then
         znotm = p40
       else
        print*, 'Wrong input uref value:',uref
       endif

      END SUBROUTINE znot_m_v6

      SUBROUTINE znot_t_v6(uref, znott)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate scalar roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
!
! Bin Liu, NOAA/NCEP/EMC 2017
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znott
      real(kind=kind_phys), parameter  :: p00 =  1.100000000000000e-04,
     &      p15 = -9.144581627678278e-10, p14 =  7.020346616456421e-08,
     &      p13 = -2.155602086883837e-06, p12 =  3.333848806567684e-05,
     &      p11 = -2.628501274963990e-04, p10 =  8.634221567969181e-04,

     &      p25 = -8.654513012535990e-12, p24 =  1.232380050058077e-09,
     &      p23 = -6.837922749505057e-08, p22 =  1.871407733439947e-06,
     &      p21 = -2.552246987137160e-05, p20 =  1.428968311457630e-04,

     &      p35 =  3.207515102100162e-12, p34 = -2.945761895342535e-10,
     &      p33 =  8.788972147364181e-09, p32 = -3.814457439412957e-08,
     &      p31 = -2.448983648874671e-06, p30 =  3.436721779020359e-05,

     &      p45 = -3.530687797132211e-11, p44 =  3.939867958963747e-09,
     &      p43 = -1.227668406985956e-08, p42 = -1.367469811838390e-05,
     &      p41 =  5.988240863928883e-04, p40 = -7.746288511324971e-03,

     &      p56 = -1.187982453329086e-13, p55 =  4.801984186231693e-11,
     &      p54 = -8.049200462388188e-09, p53 =  7.169872601310186e-07,
     &      p52 = -3.581694433758150e-05, p51 =  9.503919224192534e-04,
     &      p50 = -1.036679430885215e-02,

     &      p60 =  4.751256171799112e-05

      if (uref >= 0.0 .and. uref < 5.9 ) then
         znott = p00
      elseif (uref >= 5.9 .and. uref <= 15.4) then
         znott = p10 + uref * (p11 + uref * (p12 + uref * (p13
     &               + uref * (p14 + uref * p15))))
      elseif (uref > 15.4 .and. uref <= 21.6) then
         znott = p20 + uref * (p21 + uref * (p22 + uref * (p23 
     &               + uref * (p24 + uref * p25))))
      elseif (uref > 21.6 .and. uref <= 42.2) then
         znott = p30 + uref * (p31 + uref * (p32 + uref * (p33 
     &               + uref * (p34 + uref * p35))))
      elseif ( uref > 42.2 .and. uref <= 53.3) then
         znott = p40 + uref * (p41 + uref * (p42 + uref * (p43 
     &               + uref * (p44 + uref * p45))))
      elseif ( uref > 53.3 .and. uref <= 80.0) then
         znott = p50 + uref * (p51 + uref * (p52 + uref * (p53 
     &               + uref * (p54 + uref * (p55 + uref * p56)))))
      elseif ( uref > 80.0) then
         znott = p60
      else
         print*, 'Wrong input uref value:',uref
      endif

      END SUBROUTINE znot_t_v6


      SUBROUTINE znot_m_v7(uref, znotm)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate areodynamical roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
! For high winds, try to fit available observational data
! Comparing to znot_t_v6, slightly decrease Cd for higher wind speed
!
! Bin Liu, NOAA/NCEP/EMC 2018
!
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znotm

      real(kind=kind_phys), parameter  :: p13 = -1.296521881682694e-02,
     &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,
     &      p10 = -8.396975715683501e+00,

     &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,
     &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,
     &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,

     &      p35 =  1.897534489606422e-07, p34 = -3.019495980684978e-05,
     &      p33 =  1.931392924987349e-03, p32 = -6.797293095862357e-02,
     &      p31 =  1.346757797103756e+00, p30 = -1.707846930193362e+01,

     &      p40 =  3.371427455376717e-04

      if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp( p10 + uref * (p11 + uref * (p12 + uref * p13)))
      elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23
     &              + uref * (p24 + uref * p25))))
      elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33
     &                   + uref * (p34 + uref * p35)))))
      elseif ( uref > 53.0) then
        znotm = p40
      else
        print*, 'Wrong input uref value:',uref
      endif

      END SUBROUTINE znot_m_v7
      SUBROUTINE znot_t_v7(uref, znott)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate scalar roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
! To be compatible with the slightly decreased Cd for higher wind speed
!
! Bin Liu, NOAA/NCEP/EMC 2018
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znott

      real(kind=kind_phys), parameter  :: p00 =  1.100000000000000e-04,

     &      p15 = -9.193764479895316e-10, p14 =  7.052217518653943e-08,
     &      p13 = -2.163419217747114e-06, p12 =  3.342963077911962e-05,
     &      p11 = -2.633566691328004e-04, p10 =  8.644979973037803e-04,

     &      p25 = -9.402722450219142e-12, p24 =  1.325396583616614e-09,
     &      p23 = -7.299148051141852e-08, p22 =  1.982901461144764e-06,
     &      p21 = -2.680293455916390e-05, p20 =  1.484341646128200e-04,

     &      p35 =  7.921446674311864e-12, p34 = -1.019028029546602e-09,
     &      p33 =  5.251986927351103e-08, p32 = -1.337841892062716e-06,
     &      p31 =  1.659454106237737e-05, p30 = -7.558911792344770e-05,

     &      p45 = -2.694370426850801e-10, p44 =  5.817362913967911e-08,
     &      p43 = -5.000813324746342e-06, p42 =  2.143803523428029e-04,
     &      p41 = -4.588070983722060e-03, p40 =  3.924356617245624e-02,

     &      p56 = -1.663918773476178e-13, p55 =  6.724854483077447e-11,
     &      p54 = -1.127030176632823e-08, p53 =  1.003683177025925e-06,
     &      p52 = -5.012618091180904e-05, p51 =  1.329762020689302e-03,
     &      p50 = -1.450062148367566e-02, p60 =  6.840803042788488e-05

        if (uref >= 0.0 .and. uref < 5.9 ) then
           znott = p00
        elseif (uref >= 5.9 .and. uref <= 15.4) then
           znott = p10 + uref * (p11 + uref * (p12 + uref * (p13
     &                 + uref * (p14 + uref * p15))))
        elseif (uref > 15.4 .and. uref <= 21.6) then
           znott = p20 + uref * (p21 + uref * (p22 + uref * (p23
     &                 + uref * (p24 + uref * p25))))
        elseif (uref > 21.6 .and. uref <= 42.6) then
           znott = p30 + uref * (p31 + uref * (p32 + uref * (p33
     &                 + uref * (p34 + uref * p35))))
        elseif ( uref > 42.6 .and. uref <= 53.0) then
           znott = p40 + uref * (p41 + uref * (p42 + uref * (p43
     &                 + uref * (p44 + uref * p45))))
        elseif ( uref > 53.0 .and. uref <= 80.0) then
           znott = p50 + uref * (p51 + uref * (p52 + uref * (p53
     &                 + uref * (p54 + uref * (p55 + uref * p56)))))
        elseif ( uref > 80.0) then
           znott = p60
        else
           print*, 'Wrong input uref value:',uref
        endif

        END SUBROUTINE znot_t_v7
!> @}
      end module sfc_diff
