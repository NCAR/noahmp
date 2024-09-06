#define CCPP
!>  \file noahmpdrv.F90
!!  This file contains the NoahMP land surface scheme driver.

!>\defgroup NoahMP_LSM NoahMP LSM Model
!! \brief This is the NoahMP LSM driver module, with the functionality of 
!! preparing variables to run the NoahMP LSM subroutine noahmp_sflx(), calling NoahMP LSM and post-processing
!! variables for return to the parent model suite including unit conversion, as well
!! as diagnotics calculation.

!> This module contains the CCPP-compliant NoahMP land surface model driver.
      module noahmpdrv

      use module_sf_noahmplsm

      implicit none

      integer, parameter :: psi_opt = 0 ! 0: MYNN or 1:GFS

      private

      public :: noahmpdrv_init, noahmpdrv_run

      contains

!> \ingroup NoahMP_LSM
!! \brief This subroutine is called during the CCPP initialization phase and calls set_soilveg() to 
!! initialize soil and vegetation parameters for the chosen soil and vegetation data sources.
!! \section arg_table_noahmpdrv_init Argument Table
!! \htmlinclude noahmpdrv_init.html
!!
      subroutine noahmpdrv_init(lsm, lsm_noahmp, me, isot, ivegsrc, &
                                nlunit, pores, resid,               &
                                do_mynnsfclay,do_mynnedmf,          &
                                errmsg, errflg)

        use machine,          only: kind_phys
        use set_soilveg_mod,  only: set_soilveg
        use namelist_soilveg
        use noahmp_tables

        implicit none
        integer,              intent(in) :: lsm
        integer,              intent(in) :: lsm_noahmp    
        integer,              intent(in)  :: me, isot, ivegsrc, nlunit

        real (kind=kind_phys), dimension(:), intent(out) :: pores, resid

        logical,              intent(in) :: do_mynnsfclay
        logical,              intent(in) :: do_mynnedmf


        character(len=*),     intent(out) :: errmsg
        integer,              intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Consistency checks
        if (lsm/=lsm_noahmp) then
          write(errmsg,'(*(a))') 'Logic error: namelist choice of ',   &
       &       'LSM is different from Noah'
          errflg = 1
          return
        end if

        if (ivegsrc /= 1) then
          errmsg = 'The NOAHMP LSM expects that the ivegsrc physics '// &
                   'namelist parameter is 1. Exiting...'
          errflg = 1
          return
        end if
        if (isot /= 1) then
          errmsg = 'The NOAHMP LSM expects that the isot physics '// &
                   'namelist parameter is 1. Exiting...'
          errflg = 1
          return
        end if

        if ( do_mynnsfclay .and. .not. do_mynnedmf) then
          errmsg = 'Problem : do_mynnsfclay = .true.' // &
                   'but mynnpbl is .false.. Exiting ...'
          errflg = 1
          return
        end if


        !--- initialize soil vegetation
        call set_soilveg(me, isot, ivegsrc, nlunit, errmsg, errflg)

        !--- read in noahmp table
        call read_mp_table_parameters(errmsg, errflg)

        ! initialize psih and psim 

        if ( do_mynnsfclay ) then
        call psi_init(psi_opt,errmsg,errflg)
        endif

        pores (:) = maxsmc (:)
        resid (:) = drysmc (:)

      end subroutine noahmpdrv_init

!> \ingroup NoahMP_LSM
!! \brief This subroutine is the main CCPP entry point for the NoahMP LSM.
!! \section arg_table_noahmpdrv_run Argument Table
!! \htmlinclude noahmpdrv_run.html
!!
!! \section general_noahmpdrv NoahMP Driver General Algorithm
!!  @{
!!    - Initialize CCPP error handling variables.
!!    - Set a flag to only continue with each grid cell if the fraction of land is non-zero.
!!    - This driver may be called as part of an iterative loop. If called as the first "guess" run, 
!!        save land-related prognostic fields to restore.
!!    - Initialize output variables to zero and prepare variables for input into the NoahMP LSM.
!!    - Call transfer_mp_parameters() to fill a derived datatype for input into the NoahMP LSM.
!!    - Call noahmp_options() to set module-level scheme options for the NoahMP LSM.
!!    - If the vegetation type is ice for the grid cell, call noahmp_options_glacier() to set 
!!        module-level scheme options for NoahMP Glacier and call noahmp_glacier().
!!    - For other vegetation types, call noahmp_sflx(), the entry point of the NoahMP LSM.
!!    - Set output variables from the output of noahmp_glacier() and/or noahmp_sflx().
!!    - Call penman() to calculate potential evaporation.
!!    - Calculate the surface specific humidity and convert surface sensible and latent heat fluxes in W m-2 from their kinematic values.
!!    - If a "guess" run, restore the land-related prognostic fields.
  subroutine noahmpdrv_run                                       &
!...................................
!  ---  inputs:
    ( im, km, lsnowl, itime, ps, u1, v1, t1, q1, soiltyp,soilcol,&
      vegtype, sigmaf, dlwflx, dswsfc, snet, delt, tg3, cm, ch,  &
      prsl1, prslk1, prslki, prsik1, zf,pblh, dry, wind, slopetyp,&
      shdmin, shdmax, snoalb, sfalb, flag_iter,con_g,            &
      idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc, iopt_frz,   &
      iopt_inf, iopt_rad, iopt_alb, iopt_snf, iopt_tbot,iopt_stc,&
      iopt_trs,iopt_diag,xlatin, xcoszin, iyrlen, julian, garea, &
      rainn_mp, rainc_mp, snow_mp, graupel_mp, ice_mp, rhonewsn1,&
      con_hvap, con_cp, con_jcal, rhoh2o, con_eps, con_epsm1,    &
      con_fvirt, con_rd, con_hfus, thsfc_loc, cpllnd, cpllnd2atm,&

!  ---  in/outs:
      weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
      canopy, trans, tsurf, zorl,                                &
      rb1, fm1, fh1, ustar1, stress1, fm101, fh21,               &
      rmol1,flhc1,flqc1,do_mynnsfclay,                           &

! --- Noah MP specific

      snowxy, tvxy, tgxy, canicexy, canliqxy, eahxy, tahxy, cmxy,&
      chxy, fwetxy, sneqvoxy, alboldxy, qsnowxy, wslakexy, zwtxy,&
      waxy, wtxy, tsnoxy, zsnsoxy, snicexy, snliqxy, lfmassxy,   &
      rtmassxy, stmassxy, woodxy, stblcpxy, fastcpxy, xlaixy,    &
      xsaixy, taussxy, smoiseq, smcwtdxy, deeprechxy, rechxy,    &
      albdvis, albdnir,  albivis,  albinir,emiss,                &

!  ---  outputs:
      sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
      cmm, chh, evbs, evcw, sbsno, pah, ecan, etran, edir, snowc,&
      stm, snohf,smcwlt2, smcref2, wet1, t2mmp, q2mp,zvfun,      &
      ztmax, rca, errmsg, errflg,                                &
      canopy_heat_storage_ccpp,                                  &
      rainfall_ccpp,                                             &
      sw_absorbed_total_ccpp,                                    &
      sw_reflected_total_ccpp,                                   &
      lw_absorbed_total_ccpp,                                    &
      temperature_bare_grd_ccpp,                                 &
      temperature_veg_grd_ccpp,                                  &
      temperature_veg_2m_ccpp,                                   &
      temperature_bare_2m_ccpp,                                  &
      spec_humidity_veg_2m_ccpp,                                 &
      spec_humidity_bare_2m_ccpp,                                &
      sw_absorbed_veg_ccpp,                                      &
      sw_absorbed_ground_ccpp,                                   &
      snowmelt_out_ccpp,                                         &
      snowmelt_shallow_ccpp,                                     &
      albedo_direct_snow_ccpp,                                   &
      albedo_diffuse_snow_ccpp,                                  &
      ch_vegetated_ccpp,                                         &
      ch_bare_ground_ccpp,                                       &
      sensible_heat_grd_veg_ccpp,                                &
      sensible_heat_leaf_ccpp,                                   &
      sensible_heat_grd_bar_ccpp,                                &
      latent_heat_grd_veg_ccpp,                                  &
      latent_heat_grd_bare_ccpp,                                 &
      ground_heat_veg_ccpp,                                      &
      ground_heat_bare_ccpp,                                     &
      lw_absorbed_grd_veg_ccpp,                                  &
      lw_absorbed_leaf_ccpp,                                     &
      lw_absorbed_grd_bare_ccpp,                                 &
      latent_heat_trans_ccpp,                                    &
      latent_heat_leaf_ccpp,                                     &
      ch_leaf_ccpp,                                              &
      ch_below_canopy_ccpp,                                      &
      ch_vegetated_2m_ccpp,                                      &
      ch_bare_ground_2m_ccpp,                                    &
      precip_adv_heat_veg_ccpp,                                  &
      precip_adv_heat_grd_v_ccpp,                                &
      precip_adv_heat_grd_b_ccpp,                                &
      spec_humid_sfc_veg_ccpp,                                   &
      spec_humid_sfc_bare_ccpp                                   &
      )

  use machine ,   only : kind_phys
  use funcphys,   only : fpvs

  use module_sf_noahmplsm,   only : gfs_stability
  use module_sf_noahmp_glacier
  use noahmp_tables

  implicit none
      
  real(kind=kind_phys), parameter  :: a2      = 17.2693882
  real(kind=kind_phys), parameter  :: a3      = 273.16
  real(kind=kind_phys), parameter  :: a4      = 35.86
  real(kind=kind_phys), parameter  :: a23m4   = a2*(a3-a4)
  real(kind=kind_phys), intent(in) :: con_g 
      
  real, parameter                  :: undefined  =  9.99e20_kind_phys

  integer, parameter               :: nsoil   = 4   ! hardwired to Noah
  integer, parameter               :: nsnow   = 3   ! max. snow layers

  integer, parameter               :: iz0tlnd = 0   ! z0t treatment option

  real(kind=kind_phys), save  :: zsoil(nsoil)
  data zsoil  / -0.1, -0.4, -1.0, -2.0 /

!
!  ---  CCPP interface fields (in call order)
!

  integer                                , intent(in)    :: im         ! horiz dimension and num of used pts
  integer                                , intent(in)    :: km         ! vertical soil layer dimension
  integer                                , intent(in)    :: lsnowl     ! lower bound for snow level arrays
  integer                                , intent(in)    :: itime      ! NOT USED
  real(kind=kind_phys), dimension(:)     , intent(in)    :: ps         ! surface pressure [Pa]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: u1         ! u-component of wind [m/s]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: v1         ! u-component of wind [m/s]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: t1         ! layer 1 temperature [K]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: q1         ! layer 1 specific humidity [kg/kg]
  integer             , dimension(:)     , intent(in)    :: soiltyp    ! soil type (integer index)
  integer             , dimension(:)     , intent(in)    :: soilcol    ! soil color (integer index)
  integer             , dimension(:)     , intent(in)    :: vegtype    ! vegetation type (integer index)
  real(kind=kind_phys), dimension(:)     , intent(in)    :: sigmaf     ! areal fractional cover of green vegetation
  real(kind=kind_phys), dimension(:)     , intent(in)    :: dlwflx     ! downward longwave radiation [W/m2]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: dswsfc     ! downward shortwave radiation [W/m2]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: snet       ! total sky sfc netsw flx into ground[W/m2]
  real(kind=kind_phys)                   , intent(in)    :: delt       ! time interval [s]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: tg3        ! deep soil temperature [K]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: cm         ! surface exchange coeff for momentum [-]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: ch         ! surface exchange coeff heat & moisture[-]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: prsl1      ! sfc layer 1 mean pressure [Pa]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: prslk1     !  exner_function_at lowest model layer

  real(kind=kind_phys), dimension(:)     , intent(in)    :: prslki     ! Exner function bt midlayer and interface at 1st layer
  real(kind=kind_phys), dimension(:)     , intent(in)    :: prsik1     ! Exner function at the ground surfac

  real(kind=kind_phys), dimension(:)     , intent(in)    :: zf         ! height of bottom layer [m]

  logical                                , intent(in)    :: do_mynnsfclay !flag for MYNN sfc layer scheme

  real(kind=kind_phys), dimension(:)     , intent(in)    :: pblh       ! height of pbl
  real(kind=kind_phys), dimension(:)     , intent(inout) :: rmol1      !
  real(kind=kind_phys), dimension(:)     , intent(inout) :: flhc1      !
  real(kind=kind_phys), dimension(:)     , intent(inout) :: flqc1      !


  logical             , dimension(:)     , intent(in)    :: dry        ! = T if a point with any land
  real(kind=kind_phys), dimension(:)     , intent(in)    :: wind       ! wind speed [m/s]
  integer             , dimension(:)     , intent(in)    :: slopetyp   ! surface slope classification
  real(kind=kind_phys), dimension(:)     , intent(in)    :: shdmin     ! min green vegetation coverage [fraction]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: shdmax     ! max green vegetation coverage [fraction]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: snoalb     ! upper bound on max albedo over deep snow
  real(kind=kind_phys), dimension(:)     , intent(inout) :: sfalb      ! mean surface albedo [fraction]
  logical             , dimension(:)     , intent(in)    :: flag_iter  !
  integer                                , intent(in)    :: idveg      ! option for dynamic vegetation
  integer                                , intent(in)    :: iopt_crs   ! option for canopy stomatal resistance
  integer                                , intent(in)    :: iopt_btr   ! option for soil moisture factor for stomatal resistance
  integer                                , intent(in)    :: iopt_run   ! option for runoff and groundwater
  integer                                , intent(in)    :: iopt_sfc   ! option for surface layer drag coeff (ch & cm)
  integer                                , intent(in)    :: iopt_frz   ! option for supercooled liquid water (or ice fraction)
  integer                                , intent(in)    :: iopt_inf   ! option for frozen soil permeability
  integer                                , intent(in)    :: iopt_rad   ! option for radiation transfer
  integer                                , intent(in)    :: iopt_alb   ! option for ground snow surface albedo
  integer                                , intent(in)    :: iopt_snf   ! option for partitioning  precipitation into rainfall & snowfall
  integer                                , intent(in)    :: iopt_tbot  ! option for lower boundary condition of soil temperature
  integer                                , intent(in)    :: iopt_stc   ! option for snow/soil temperature time scheme (only layer 1)
  integer                                , intent(in)    :: iopt_trs   ! option for thermal roughness scheme
  integer                                , intent(in)    :: iopt_diag  ! option for surface diagnose approach
  real(kind=kind_phys), dimension(:)     , intent(in)    :: xlatin     ! latitude
  real(kind=kind_phys), dimension(:)     , intent(in)    :: xcoszin    ! cosine of zenith angle
  integer                                , intent(in)    :: iyrlen     ! year length [days]
  real(kind=kind_phys)                   , intent(in)    :: julian     ! julian day of year
  real(kind=kind_phys), dimension(:)     , intent(in)    :: garea      ! area of the grid cell
  real(kind=kind_phys), dimension(:)     , intent(in), optional :: rainn_mp   ! microphysics non-convective precipitation [mm]
  real(kind=kind_phys), dimension(:)     , intent(in), optional :: rainc_mp   ! microphysics convective precipitation [mm]
  real(kind=kind_phys), dimension(:)     , intent(in), optional :: snow_mp    ! microphysics snow [mm]
  real(kind=kind_phys), dimension(:)     , intent(in), optional :: graupel_mp ! microphysics graupel [mm]
  real(kind=kind_phys), dimension(:)     , intent(in), optional :: ice_mp     ! microphysics ice/hail [mm]
  real(kind=kind_phys), dimension(:)     , intent(in)    :: rhonewsn1  ! precipitation ice density (kg/m^3)
  real(kind=kind_phys)                   , intent(in)    :: con_hvap   ! latent heat condensation [J/kg]
  real(kind=kind_phys)                   , intent(in)    :: con_cp     ! specific heat air [J/kg/K] 
  real(kind=kind_phys)                   , intent(in)    :: con_jcal   ! joules per calorie (not used)
  real(kind=kind_phys)                   , intent(in)    :: rhoh2o     ! density of water [kg/m^3]
  real(kind=kind_phys)                   , intent(in)    :: con_eps    ! Rd/Rv 
  real(kind=kind_phys)                   , intent(in)    :: con_epsm1  ! Rd/Rv - 1
  real(kind=kind_phys)                   , intent(in)    :: con_fvirt  ! Rv/Rd - 1
  real(kind=kind_phys)                   , intent(in)    :: con_rd     ! gas constant air [J/kg/K]
  real(kind=kind_phys)                   , intent(in)    :: con_hfus   ! lat heat H2O fusion  [J/kg]

  logical                                , intent(in)    :: thsfc_loc  ! Flag for reference pressure in theta calculation

  logical                                , intent(in)    :: cpllnd     ! Flag for land coupling (atm->lnd)
  logical                                , intent(in)    :: cpllnd2atm ! Flag for land coupling (lnd->atm)

  real(kind=kind_phys), dimension(:)     , intent(inout) :: weasd      ! water equivalent accumulated snow depth [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: snwdph     ! snow depth [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: tskin      ! ground surface skin temperature [K]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: tprcp      ! total precipitation [m]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: srflag     ! snow/rain flag for precipitation
  real(kind=kind_phys), dimension(:,:)   , intent(inout) :: smc        ! total soil moisture content [m3/m3]
  real(kind=kind_phys), dimension(:,:)   , intent(inout) :: stc        ! soil temp [K]
  real(kind=kind_phys), dimension(:,:)   , intent(inout) :: slc        ! liquid soil moisture [m3/m3]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: canopy     ! canopy moisture content [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: trans      ! total plant transpiration [m/s]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: tsurf      ! surface skin temperature [K]
  real(kind=kind_phys), dimension(:)     , intent(inout) :: zorl       ! surface roughness [cm]

  real(kind=kind_phys), dimension(:)     , intent(inout) :: rb1        ! bulk richardson #
  real(kind=kind_phys), dimension(:)     , intent(inout) :: fm1        !  Monin_Obukhov_silarity_function for momentum
  real(kind=kind_phys), dimension(:)     , intent(inout) :: fh1        !  Monin_Obukhov_silarity_function for heat
  real(kind=kind_phys), dimension(:)     , intent(inout) :: ustar1     !  friction velocity m s-1
  real(kind=kind_phys), dimension(:)     , intent(inout) :: stress1    ! Wind stress m2 S-2
  real(kind=kind_phys), dimension(:)     , intent(inout) :: fm101      ! MOS function for momentum evaulated @ 10 m
  real(kind=kind_phys), dimension(:)     , intent(inout) :: fh21       ! MOS function for heat evaulated @ 2m

  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: snowxy     ! actual no. of snow layers
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: tvxy       ! vegetation leaf temperature [K]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: tgxy       ! bulk ground surface temperature [K]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: canicexy   ! canopy-intercepted ice [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: canliqxy   ! canopy-intercepted liquid water [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: eahxy      ! canopy air vapor pressure [Pa]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: tahxy      ! canopy air temperature [K]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: cmxy       ! bulk momentum drag coefficient [m/s]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: chxy       ! bulk sensible heat exchange coefficient [m/s]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: fwetxy     ! wetted or snowed fraction of the canopy [-]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: sneqvoxy   ! snow mass at last time step[mm h2o]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: alboldxy   ! snow albedo at last time step [-]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: qsnowxy    ! snowfall on the ground [mm/s]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: wslakexy   ! lake water storage [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: zwtxy      ! water table depth [m]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: waxy       ! water in the "aquifer" [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: wtxy       ! groundwater storage [mm]
  real(kind=kind_phys), dimension(:,lsnowl:), intent(inout), optional :: tsnoxy  ! snow temperature [K]
  real(kind=kind_phys), dimension(:,lsnowl:), intent(inout), optional :: zsnsoxy ! snow/soil layer depth [m]
  real(kind=kind_phys), dimension(:,lsnowl:), intent(inout), optional :: snicexy ! snow layer ice [mm]
  real(kind=kind_phys), dimension(:,lsnowl:), intent(inout), optional :: snliqxy ! snow layer liquid water [mm]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: lfmassxy   ! leaf mass [g/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: rtmassxy   ! mass of fine roots [g/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: stmassxy   ! stem mass [g/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: woodxy     ! mass of wood (incl. woody roots) [g/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: stblcpxy   ! stable carbon in deep soil [g/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: fastcpxy   ! short-lived carbon, shallow soil [g/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: xlaixy     ! leaf area index [m2/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: xsaixy     ! stem area index [m2/m2]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: taussxy    ! snow age factor [-]
  real(kind=kind_phys), dimension(:,:)   , intent(inout), optional :: smoiseq    ! eq volumetric soil moisture [m3/m3]
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: smcwtdxy   ! soil moisture content in the layer to the water table when deep
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: deeprechxy ! recharge to the water table when deep
  real(kind=kind_phys), dimension(:)     , intent(inout), optional :: rechxy     ! recharge to the water table
  real(kind=kind_phys), dimension(:)     , intent(out)   :: albdvis    ! albedo - direct  visible [fraction]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: albdnir    ! albedo - direct  NIR     [fraction]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: albivis    ! albedo - diffuse visible [fraction]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: albinir    ! albedo - diffuse NIR     [fraction]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: emiss      ! sfc lw emissivity [fraction]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: sncovr1    ! snow cover over land [fraction]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: qsurf      ! specific humidity at sfc [kg/kg]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: gflux      ! soil heat flux [W/m2]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: drain      ! subsurface runoff [mm/s]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: evap       ! total latent heat flux [W/m2]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: hflx       ! sensible heat flux [W/m2]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: ep         ! potential evaporation [mm/s?]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: runoff     ! surface runoff [mm/s]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: cmm        ! cm*U [m/s]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: chh        ! ch*U*rho [kg/m2/s]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: evbs       ! direct soil evaporation [m/s]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: evcw       ! canopy water evaporation [m/s]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: sbsno      ! sublimation/deposit from snopack [W/m2]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: pah        ! precipitation advected heat - total (w/m2) 
  real(kind=kind_phys), dimension(:)     , intent(out)   :: ecan       ! evaporation of intercepted water (mm/s)
  real(kind=kind_phys), dimension(:)     , intent(out)   :: etran      ! transpiration rate (mm/s) 
  real(kind=kind_phys), dimension(:)     , intent(out)   :: edir       ! soil surface evaporation rate (mm/s)
  real(kind=kind_phys), dimension(:)     , intent(out)   :: snowc      ! fractional snow cover [-]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: stm        ! total soil column moisture content [mm]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: snohf      ! snow/freezing-rain latent heat flux [W/m2]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: smcwlt2    ! dry soil moisture threshold [m3/m3]
  real(kind=kind_phys), dimension(:)     , intent(out)   :: smcref2    ! soil moisture threshold [m3/m3]
  real(kind=kind_phys), dimension(:)     , intent(out), optional :: wet1       ! normalized surface soil saturated fraction
  real(kind=kind_phys), dimension(:)     , intent(out), optional :: t2mmp      ! combined T2m from tiles
  real(kind=kind_phys), dimension(:)     , intent(out), optional :: q2mp       ! combined q2m from tiles
  real(kind=kind_phys), dimension(:)     , intent(out)   :: zvfun      ! 
  real(kind=kind_phys), dimension(:)     , intent(out)   :: ztmax      ! thermal roughness length
  real(kind=kind_phys), dimension(:)     , intent(out), optional :: rca        ! total canopy/stomatal resistance (s/m)

  character(len=*)    ,                    intent(out)   :: errmsg
  integer             ,                    intent(out)   :: errflg

  real(kind=kind_phys), dimension(:)  , intent(out), optional :: canopy_heat_storage_ccpp ! within-canopy heat [W/m2]
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: rainfall_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: sw_absorbed_total_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: sw_reflected_total_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: lw_absorbed_total_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: temperature_bare_grd_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: temperature_veg_grd_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: temperature_veg_2m_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: temperature_bare_2m_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: spec_humidity_veg_2m_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: spec_humidity_bare_2m_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: sw_absorbed_veg_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: sw_absorbed_ground_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: snowmelt_out_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: snowmelt_shallow_ccpp
  real(kind=kind_phys), dimension(:,:), intent(out), optional :: albedo_direct_snow_ccpp
  real(kind=kind_phys), dimension(:,:), intent(out), optional :: albedo_diffuse_snow_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ch_vegetated_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ch_bare_ground_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: sensible_heat_grd_veg_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: sensible_heat_leaf_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: sensible_heat_grd_bar_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: latent_heat_grd_veg_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: latent_heat_grd_bare_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ground_heat_veg_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ground_heat_bare_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: lw_absorbed_grd_veg_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: lw_absorbed_leaf_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: lw_absorbed_grd_bare_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: latent_heat_trans_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: latent_heat_leaf_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ch_leaf_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ch_below_canopy_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ch_vegetated_2m_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: ch_bare_ground_2m_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: precip_adv_heat_veg_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: precip_adv_heat_grd_v_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: precip_adv_heat_grd_b_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: spec_humid_sfc_veg_ccpp
  real(kind=kind_phys), dimension(:)  , intent(out), optional :: spec_humid_sfc_bare_ccpp

!
!  ---  some new options, hard code for now
!

  integer    :: iopt_rsf  = 4 ! option for surface resistance
  integer    :: iopt_soil = 1 ! option for soil parameter treatment
  integer    :: iopt_pedo = 1 ! option for pedotransfer function
  integer    :: iopt_crop = 0 ! option for crop model
  integer    :: iopt_gla  = 2 ! option for glacier treatment
  integer    :: iopt_z0m  = 1 ! option for z0m treatment

!
!  ---  local inputs to noah-mp and glacier subroutines; listed in order in noah-mp call
!
                                                                            ! intent
  integer                                          :: i_location            ! in    | grid index
  integer                                          :: j_location            ! in    | grid index (not used in ccpp)
  real (kind=kind_phys)                            :: latitude              ! in    | latitude [radians]
  integer                                          :: year_length           ! in    | number of days in the current year
  real (kind=kind_phys)                            :: julian_day            ! in    | julian day of year [floating point]
  real (kind=kind_phys)                            :: cosine_zenith         ! in    | cosine solar zenith angle [-1,1]
  real (kind=kind_phys)                            :: timestep              ! in    | time step [sec]
  real (kind=kind_phys)                            :: spatial_scale         ! in    | spatial scale [m] (not used in noah-mp)
  real (kind=kind_phys)                            :: atmosphere_thickness  ! in    | thickness of lowest atmo layer [m] (not used in noah-mp)
  integer                                          :: soil_levels           ! in    | soil levels
  real (kind=kind_phys), dimension(       1:nsoil) :: soil_interface_depth  ! in    | soil layer-bottom depth from surface [m]
  integer                                          :: max_snow_levels       ! in    | maximum number of snow levels
  real (kind=kind_phys)                            :: vegetation_frac       ! in    | vegetation fraction [0.0-1.0]
  real (kind=kind_phys)                            :: area_grid             ! in    | 
  real (kind=kind_phys)                            :: max_vegetation_frac   ! in    | annual maximum vegetation fraction [0.0-1.0]
  integer                                          :: vegetation_category   ! in    | vegetation category
  integer                                          :: ice_flag              ! in    | ice flag (1->ice)
  integer                                          :: surface_type          ! in    | surface type flag 1->soil; 2->lake
  integer                                          :: crop_type             ! in    | crop type category
  real (kind=kind_phys), dimension(       1:nsoil) :: eq_soil_water_vol     ! in    | (opt_run=5) equilibrium soil water content [m3/m3]
  real (kind=kind_phys)                            :: temperature_forcing   ! in    | forcing air temperature [K]
  real (kind=kind_phys)                            :: air_pressure_surface  ! in    | surface air pressure [Pa]
  real (kind=kind_phys)                            :: air_pressure_forcing  ! in    | forcing air pressure [Pa]
  real (kind=kind_phys)                            :: uwind_forcing         ! in    | forcing u-wind [m/s]
  real (kind=kind_phys)                            :: vwind_forcing         ! in    | forcing v-wind [m/s]
  real (kind=kind_phys)                            :: spec_humidity_forcing ! in    | forcing mixing ratio [kg/kg]
  real (kind=kind_phys)                            :: cloud_water_forcing   ! in    | cloud water mixing ratio [kg/kg] (not used in noah-mp)
  real (kind=kind_phys)                            :: sw_radiation_forcing  ! in    | forcing downward shortwave radiation [W/m2]
  real (kind=kind_phys)                            :: radiation_lw_forcing  ! in    | forcing downward longwave radiation [W/m2]
  real (kind=kind_phys)                            :: precipitation_forcing ! in    | total precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_convective     ! in    | convective precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_non_convective ! in    | non-convective precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_sh_convective  ! in    | shallow convective precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_snow           ! in    | snow precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_graupel        ! in    | graupel precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_hail           ! in    | hail precipitation [mm/s]
  real (kind=kind_phys)                            :: temperature_soil_bot  ! in    | soil bottom boundary condition temperature [K]
  real (kind=kind_phys)                            :: co2_air               ! in    | atmospheric co2 concentration [Pa]
  real (kind=kind_phys)                            :: o2_air                ! in    | atmospheric o2 concentration [Pa]
  real (kind=kind_phys)                            :: foliage_nitrogen      ! in    | foliage nitrogen [%] [1-saturated]
  real (kind=kind_phys), dimension(-nsnow+1:    0) :: snow_ice_frac_old     ! in    | snow ice fraction at last timestep [-]
  real (kind=kind_phys)                            :: forcing_height        ! inout | forcing height [m]
  real (kind=kind_phys)                            :: snow_albedo_old       ! inout | snow albedo at last time step (class option) [-]
  real (kind=kind_phys)                            :: snow_water_equiv_old  ! inout | snow water equivalent at last time step [mm]
  real (kind=kind_phys), dimension(-nsnow+1:nsoil) :: temperature_snow_soil ! inout | snow/soil temperature [K]
  real (kind=kind_phys), dimension(       1:nsoil) :: soil_liquid_vol       ! inout | volumetric liquid soil moisture [m3/m3]
  real (kind=kind_phys), dimension(       1:nsoil) :: soil_moisture_vol     ! inout | volumetric soil moisture (ice + liq.) [m3/m3]

  real (kind=kind_phys)                            :: surface_temperature   !  out  | surface aerodynamic temp

  real (kind=kind_phys)                            :: temperature_canopy_air! inout | canopy air tmeperature [K]
  real (kind=kind_phys)                            :: vapor_pres_canopy_air ! inout | canopy air vapor pressure [Pa]
  real (kind=kind_phys)                            :: canopy_wet_fraction   ! inout | wetted or snowed fraction of canopy (-)
  real (kind=kind_phys)                            :: canopy_liquid         ! inout | canopy intercepted liquid [mm]
  real (kind=kind_phys)                            :: canopy_ice            ! inout | canopy intercepted ice [mm]
  real (kind=kind_phys)                            :: temperature_leaf      ! inout | leaf temperature [K]
  real (kind=kind_phys)                            :: temperature_ground    ! inout | grid ground surface temperature [K]
  real (kind=kind_phys)                            :: spec_humidity_surface ! inout | surface specific humidty [kg/kg]
  real (kind=kind_phys)                            :: snowfall              ! inout | land model partitioned snowfall [mm/s]
  real (kind=kind_phys)                            :: rainfall              ! inout | land model partitioned rainfall [mm/s]
  integer                                          :: snow_levels           ! inout | active snow levels [-]
  real (kind=kind_phys), dimension(-nsnow+1:nsoil) :: interface_depth       ! inout | layer-bottom depth from snow surf [m]
  real (kind=kind_phys)                            :: snow_depth            ! inout | snow depth [m]
  real (kind=kind_phys)                            :: snow_water_equiv      ! inout | snow water equivalent [mm]
  real (kind=kind_phys), dimension(-nsnow+1:    0) :: snow_level_ice        ! inout | snow level ice [mm]
  real (kind=kind_phys), dimension(-nsnow+1:    0) :: snow_level_liquid     ! inout | snow level liquid [mm]
  real (kind=kind_phys)                            :: depth_water_table     ! inout | depth to water table [m]
  real (kind=kind_phys)                            :: aquifer_water         ! inout | water storage in aquifer [mm]
  real (kind=kind_phys)                            :: saturated_water       ! inout | water in aquifer+saturated soil [mm]
  real (kind=kind_phys)                            :: lake_water            ! inout | lake water storage (can be neg.) [mm]
  real (kind=kind_phys)                            :: leaf_carbon           ! inout | leaf mass [g/m2]
  real (kind=kind_phys)                            :: root_carbon           ! inout | mass of fine roots [g/m2]
  real (kind=kind_phys)                            :: stem_carbon           ! inout | stem mass [g/m2]
  real (kind=kind_phys)                            :: wood_carbon           ! inout | mass of wood (incl. woody roots) [g/m2]
  real (kind=kind_phys)                            :: soil_carbon_stable    ! inout | stable soil carbon [g/m2]
  real (kind=kind_phys)                            :: soil_carbon_fast      ! inout | short-lived soil carbon [g/m2]
  real (kind=kind_phys)                            :: leaf_area_index       ! inout | leaf area index [-]
  real (kind=kind_phys)                            :: stem_area_index       ! inout | stem area index [-]
  real (kind=kind_phys)                            :: cm_noahmp             ! inout | grid momentum drag coefficient [m/s]
  real (kind=kind_phys)                            :: ch_noahmp             ! inout | grid heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: snow_age              ! inout | non-dimensional snow age [-]
  real (kind=kind_phys)                            :: grain_carbon          ! inout | grain mass [g/m2]
  real (kind=kind_phys)                            :: growing_deg_days      ! inout | growing degree days [-]
  integer                                          :: plant_growth_stage    ! inout | plant growing stage [-]
  real (kind=kind_phys)                            :: soil_moisture_wtd     ! inout | (opt_run=5) soil water content between bottom of the soil and water table [m3/m3]
  real (kind=kind_phys)                            :: deep_recharge         ! inout | (opt_run=5) recharge to or from the water table when deep [m]
  real (kind=kind_phys)                            :: recharge              ! inout | (opt_run=5) recharge to or from the water table when shallow [m] (diagnostic)
  real (kind=kind_phys)                            :: z0_total              !   out | weighted z0 sent to coupled model [m]
  real (kind=kind_phys)                            :: z0h_total             !   out | weighted z0h sent to coupled model [m]

  real (kind=kind_phys)                            :: sw_absorbed_total     !   out | total absorbed solar radiation [W/m2]
  real (kind=kind_phys)                            :: sw_reflected_total    !   out | total reflected solar radiation [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_total     !   out | total net lw rad [W/m2]  [+ to atm]
  real (kind=kind_phys)                            :: sensible_heat_total   !   out | total sensible heat [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: ground_heat_total     !   out | ground heat flux [W/m2]   [+ to soil]
  real (kind=kind_phys)                            :: latent_heat_canopy    !   out | canopy evaporation heat flux [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: latent_heat_ground    !   out | ground evaporation heat flux [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: transpiration_heat    !   out | transpiration heat flux [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: evaporation_canopy    !   out | canopy evaporation [mm/s]
  real (kind=kind_phys)                            :: transpiration         !   out | transpiration [mm/s]
  real (kind=kind_phys)                            :: evaporation_soil      !   out | soil surface evaporation [mm/s]
  real (kind=kind_phys)                            :: temperature_radiative !   out | surface radiative temperature [K]
  real (kind=kind_phys)                            :: temperature_bare_grd  !   out | bare ground surface temperature [K]
  real (kind=kind_phys)                            :: temperature_veg_grd   !   out | below_canopy ground surface temperature [K]
  real (kind=kind_phys)                            :: temperature_veg_2m    !   out | vegetated 2-m air temperature [K]
  real (kind=kind_phys)                            :: temperature_bare_2m   !   out | bare ground 2-m air temperature [K]
  real (kind=kind_phys)                            :: spec_humidity_veg_2m  !   out | vegetated 2-m air specific humidity [K]
  real (kind=kind_phys)                            :: spec_humidity_bare_2m !   out | bare ground 2-m air specfic humidity [K]
  real (kind=kind_phys)                            :: runoff_surface        !   out | surface runoff [mm/s] 
  real (kind=kind_phys)                            :: runoff_baseflow       !   out | baseflow runoff [mm/s]
  real (kind=kind_phys)                            :: par_absorbed          !   out | absorbed photosynthesis active radiation [W/m2]
  real (kind=kind_phys)                            :: photosynthesis        !   out | total photosynthesis [umol CO2/m2/s] [+ out]
  real (kind=kind_phys)                            :: sw_absorbed_veg       !   out | solar radiation absorbed by vegetation [W/m2]
  real (kind=kind_phys)                            :: sw_absorbed_ground    !   out | solar radiation absorbed by ground [W/m2]
  real (kind=kind_phys)                            :: snow_cover_fraction   !   out | snow cover fraction on the ground [-]
  real (kind=kind_phys)                            :: net_eco_exchange      !   out | net ecosystem exchange [g/m2/s CO2]
  real (kind=kind_phys)                            :: global_prim_prod      !   out | global primary production [g/m2/s C]
  real (kind=kind_phys)                            :: net_prim_prod         !   out | net primary productivity [g/m2/s C]
  real (kind=kind_phys)                            :: vegetation_fraction   !   out | vegetation fraction [0.0-1.0]
  real (kind=kind_phys)                            :: albedo_total          !   out | total surface albedo [-]
  real (kind=kind_phys)                            :: snowmelt_out          !   out | snowmelt out bottom of pack [mm/s]
  real (kind=kind_phys)                            :: snowmelt_shallow      !   out | shallow snow melt [mm]
  real (kind=kind_phys)                            :: snowmelt_shallow_1    !   out | additional shallow snow melt [mm]
  real (kind=kind_phys)                            :: snowmelt_shallow_2    !   out | additional shallow snow melt [mm]
  real (kind=kind_phys)                            :: rs_sunlit             !   out | sunlit leaf stomatal resistance [s/m]
  real (kind=kind_phys)                            :: rs_shaded             !   out | shaded leaf stomatal resistance [s/m]
  real (kind=kind_phys), dimension(1:2)            :: albedo_direct         !   out | direct vis/nir albedo [-]
  real (kind=kind_phys), dimension(1:2)            :: albedo_diffuse        !   out | diffuse vis/nir albedo [-]
  real (kind=kind_phys), dimension(1:2)            :: albedo_direct_snow    !   out | direct vis/nir snow albedo [-]
  real (kind=kind_phys), dimension(1:2)            :: albedo_diffuse_snow   !   out | diffuse vis/nir snow albedo [-]
  real (kind=kind_phys)                            :: canopy_gap_fraction   !   out | between canopy gap fraction [-]
  real (kind=kind_phys)                            :: incanopy_gap_fraction !   out | within canopy gap fraction for beam [-]
  real (kind=kind_phys)                            :: ch_vegetated          !   out | vegetated heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_bare_ground        !   out | bare-ground heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: emissivity_total      !   out | grid emissivity [-]
  real (kind=kind_phys)                            :: sensible_heat_grd_veg !   out | below-canopy ground sensible heat flux [W/m2]
  real (kind=kind_phys)                            :: sensible_heat_leaf    !   out | leaf-to-canopy sensible heat flux [W/m2]
  real (kind=kind_phys)                            :: sensible_heat_grd_bar !   out | bare ground sensible heat flux [W/m2]
  real (kind=kind_phys)                            :: latent_heat_grd_veg   !   out | below-canopy ground evaporation heat flux [W/m2]
  real (kind=kind_phys)                            :: latent_heat_grd_bare  !   out | bare ground evaporation heat flux [W/m2]
  real (kind=kind_phys)                            :: ground_heat_veg       !   out | below-canopy ground heat flux [W/m2]
  real (kind=kind_phys)                            :: ground_heat_bare      !   out | bare ground heat flux [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_grd_veg   !   out | below-canopy ground absorbed longwave radiation [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_leaf      !   out | leaf absorbed longwave radiation [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_grd_bare  !   out | bare ground net longwave radiation [W/m2]
  real (kind=kind_phys)                            :: latent_heat_trans     !   out | transpiration [W/m2]
  real (kind=kind_phys)                            :: latent_heat_leaf      !   out | leaf evaporation [W/m2]
  real (kind=kind_phys)                            :: ch_leaf               !   out | leaf exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_below_canopy       !   out | below-canopy exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_vegetated_2m       !   out | 2-m vegetated  heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_bare_ground_2m     !   out | 2-m bare-ground heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: precip_frozen_frac    !   out | precipitation snow fraction [-]
  real (kind=kind_phys)                            :: precip_adv_heat_veg   !   out | precipitation advected heat - vegetation net [W/m2]
  real (kind=kind_phys)                            :: precip_adv_heat_grd_v !   out | precipitation advected heat - below-canopy net [W/m2]
  real (kind=kind_phys)                            :: precip_adv_heat_grd_b !   out | precipitation advected heat - bare ground net [W/m2]
  real (kind=kind_phys)                            :: precip_adv_heat_total !   out | precipitation advected heat - total [W/m2)
  real (kind=kind_phys)                            :: snow_sublimation      !   out | snow sublimation [W/m2]
  real (kind=kind_phys)                            :: lai_sunlit            !   out | sunlit leaf area index [m2/m2]
  real (kind=kind_phys)                            :: lai_shaded            !   out | shaded leaf area index [m2/m2]
  real (kind=kind_phys)                            :: leaf_air_resistance   !   out | leaf boundary layer resistance [s/m]

  real (kind=kind_phys)                            :: canopy_heat_storage   !   out | within-canopy heat [W/m2]
  real (kind=kind_phys)                            :: spec_humid_sfc_veg    !   out | surface specific humidty over vegetation [kg/kg]
  real (kind=kind_phys)                            :: spec_humid_sfc_bare   !   out | surface specific humidty over bare soil [kg/kg]
  
  real (kind=kind_phys)                            :: ustarx                !  inout |surface friction velocity
  real (kind=kind_phys)                            :: prslkix               !  in exner function
  real (kind=kind_phys)                            :: prsik1x               !  in exner function
  real (kind=kind_phys)                            :: prslk1x               !  in exner function

  real (kind=kind_phys)                            :: ch2
  real (kind=kind_phys)                            :: cq2
  real (kind=kind_phys)                            :: qfx
  real (kind=kind_phys)                            :: wspd1                 !  wind speed with all components
  real (kind=kind_phys)                            :: pblhx                 !  height of pbl
   integer                                         :: mnice

  real (kind=kind_phys)                            :: rah_total             !
  real (kind=kind_phys)                            :: cah_total             !


!
!  ---  local variable
!

  integer :: soil_category(nsoil)
  integer :: slope_category
  integer :: soil_color_category
  character(len=256)                     :: dataset_identifier

  real (kind=kind_phys) :: spec_humidity_sat      ! saturation specific humidity
  real (kind=kind_phys) :: vapor_pressure_sat     ! saturation vapor pressure
  real (kind=kind_phys) :: latent_heat_total        ! total latent heat flux [W/m2]
  real (kind=kind_phys) :: density                ! air density
  real (kind=kind_phys) :: virtual_temperature    ! used for penman calculation and density
  real (kind=kind_phys) :: potential_evaporation  ! used for penman calculation
  real (kind=kind_phys) :: potential_temperature  ! used for penman calculation
  real (kind=kind_phys) :: penman_radiation       ! used for penman calculation
  real (kind=kind_phys) :: dqsdt                  ! used for penman calculation
  real (kind=kind_phys) :: precip_freeze_frac_in  ! used for penman calculation
 
  real (kind=kind_phys) :: virtfac1               ! virtual factor
  real (kind=kind_phys) :: tflux                  ! surface flux temp
  real (kind=kind_phys) :: tvs1                   ! surface virtual temp
  real (kind=kind_phys) :: vptemp                 ! virtual potential temp

  real(kind=kind_phys) ::  tem1,tem2,gdx
  real(kind=kind_phys), parameter :: z0lo=0.1, z0up=1.0

  logical               :: is_snowing             ! used for penman calculation
  logical               :: is_freeze_rain         ! used for penman calculation
  integer :: i, k
      
!
!  --- local derived constants:
!
      
  type(noahmp_parameters) :: parameters

!
!  --- end declaration
!     

!
!  --- Initialize CCPP error handling variables
!
  errmsg = ''
  errflg = 0

!
!  --- Just return if external land component is activated for two-way interaction
!
  if (cpllnd .and. cpllnd2atm) return

  do i = 1, im

    if (flag_iter(i) .and. dry(i)) then

!
!  --- variable checks and derived fields
!

      if(vegtype(i) == isice_table ) then
        if(weasd(i) < 0.1) then
          weasd(i)  = 0.1
        end if
      end if                                       

!
!  --- noah-mp input variables (except snow_ice_frac_old done later)
!

      dataset_identifier    = "modified_igbp_modis_noah"

      i_location            = i
      j_location            = -9999
      latitude              = xlatin(i)
      year_length           = iyrlen
      julian_day            = julian
      cosine_zenith         = xcoszin(i)
      timestep              = delt
      spatial_scale         = -9999.0
      atmosphere_thickness  = -9999.0
      soil_levels           = km
      soil_interface_depth  = zsoil
      max_snow_levels       = nsnow
      vegetation_frac       = sigmaf(i)
      max_vegetation_frac   = shdmax(i)
      vegetation_category   = vegtype(i)
      surface_type          = 1
      crop_type             = 0
      eq_soil_water_vol     = smoiseq(i,:) ! only need for run=5
      temperature_forcing   = t1(i) 
      air_pressure_surface  = ps(i)
      air_pressure_forcing  = prsl1(i)
      uwind_forcing         = u1(i)
      vwind_forcing         = v1(i)
      area_grid             = garea(i)

      pblhx                 = pblh(i)

      prslkix               = prslki(i)
      prsik1x               = prsik1(i)
      prslk1x               = prslk1(i)

      spec_humidity_forcing = max(q1(i), 1.e-8)                            ! specific humidity at level 1 (kg/kg)
      virtual_temperature   = temperature_forcing * &
                               (1.0 + con_fvirt * spec_humidity_forcing)   ! virtual temperature
      vapor_pressure_sat    = fpvs( temperature_forcing )                  ! sat. vapor pressure at level 1 (Pa)
      spec_humidity_sat     = con_eps*vapor_pressure_sat / &
                              (prsl1(i) + con_epsm1*vapor_pressure_sat)    ! sat. specific humidity at level 1 (kg/kg)
      spec_humidity_sat     = max(spec_humidity_sat, 1.e-8)                ! lower limit sat. specific humidity (kg/kg)
      spec_humidity_forcing = min(spec_humidity_sat,spec_humidity_forcing) ! limit specific humidity at level 1 (kg/kg)

      cloud_water_forcing   = -9999.0
      sw_radiation_forcing  = dswsfc(i)
      radiation_lw_forcing  = dlwflx(i)
      precipitation_forcing = 1000.0 * tprcp(i) / delt 
      precip_convective     = rainc_mp(i)
      precip_non_convective = rainn_mp(i)
      precip_sh_convective  = 0.
      precip_snow           = snow_mp(i)
      precip_graupel        = graupel_mp(i)
      precip_hail           = ice_mp(i)
      temperature_soil_bot  = tg3(i)
      co2_air               = co2_table * air_pressure_forcing
      o2_air                = o2_table  * air_pressure_forcing
      foliage_nitrogen      = 1.0
      
!
!  --- noah-mp inout variables
!

      forcing_height               = zf(i)
      snow_albedo_old              = alboldxy(i)
      snow_water_equiv_old         = sneqvoxy(i)
      temperature_snow_soil(-2: 0) = tsnoxy(i,:)
      temperature_snow_soil( 1:km) = stc(i,:)
      soil_liquid_vol              = slc(i,:)
      soil_moisture_vol            = smc(i,:)
      temperature_canopy_air       = tahxy(i)
      vapor_pres_canopy_air        = eahxy(i)
      canopy_wet_fraction          = fwetxy(i)
      canopy_liquid                = canliqxy(i)
      canopy_ice                   = canicexy(i)
      temperature_leaf             = tvxy(i)
      temperature_ground           = tgxy(i)
      spec_humidity_surface        = undefined  ! doesn't need inout; should be out
      snowfall                     = qsnowxy(i) ! doesn't need inout; should be out
      rainfall                     = -9999.0    ! doesn't need inout; should be out
      snow_levels                  = nint(snowxy(i))
      interface_depth              = zsnsoxy(i,:)
      snow_depth                   = snwdph(i) * 0.001         ! convert from mm to m
      snow_water_equiv             = weasd(i)
          if (snow_water_equiv /= 0.0 .and. snow_depth == 0.0) then
            snow_depth = 10.0 * snow_water_equiv /1000.0
          endif
      snow_level_ice               = snicexy(i,:)
      snow_level_liquid            = snliqxy(i,:)
      depth_water_table            = zwtxy(i)
      aquifer_water                = waxy(i)
      saturated_water              = wtxy(i)
      lake_water                   = wslakexy(i)
      leaf_carbon                  = lfmassxy(i)
      root_carbon                  = rtmassxy(i)
      stem_carbon                  = stmassxy(i)
      wood_carbon                  = woodxy(i)
      soil_carbon_stable           = stblcpxy(i)
      soil_carbon_fast             = fastcpxy(i)
      leaf_area_index              = xlaixy(i)
      stem_area_index              = xsaixy(i)
      cm_noahmp                    = cmxy(i)
      ch_noahmp                    = chxy(i)
      snow_age                     = taussxy(i)
!      grain_carbon          ! new variable
!      growing_deg_days      ! new variable
!      plant_growth_stage    ! new variable
      soil_moisture_wtd            = smcwtdxy(i)
      deep_recharge                = deeprechxy(i)
      recharge                     = rechxy(i)
      
      ustarx                       = ustar1(i)

      snow_ice_frac_old = 0.0
      do k = snow_levels+1, 0
        if(snow_level_ice(k) > 0.0 ) &
          snow_ice_frac_old(k) = snow_level_ice(k) /(snow_level_ice(k)+snow_level_liquid(k)) 
      end do


       if (snow_depth .gt. 0.1 .or. vegetation_category  == isice_table ) then
         mnice = 1
       else    
         mnice = 0
       endif   

!
!  --- some outputs for atm model?
!
      density = air_pressure_forcing / (con_rd * virtual_temperature)
      chh(i) = ch(i)  * wind(i) * density
      cmm(i) = cm(i)  * wind(i)
!
!  --- noah-mp additional variables
!

      soil_category       = soiltyp(i)
      slope_category      = slopetyp(i)
      soil_color_category = soilcol(i)
!     soil_color_category = 4


      call transfer_mp_parameters(vegetation_category, soil_category, &
                        slope_category, soil_color_category, crop_type,parameters)
      parameters%prcpiceden = rhonewsn1(i)
      call noahmp_options(idveg ,iopt_crs, iopt_btr , iopt_run, iopt_sfc,  &
                                 iopt_frz, iopt_inf , iopt_rad, iopt_alb,  &
                                 iopt_snf, iopt_tbot, iopt_stc, iopt_rsf,  &
			         iopt_soil,iopt_pedo, iopt_crop,iopt_trs,  &
                                 iopt_diag,iopt_z0m)

      if ( vegetation_category == isice_table )  then

        if (precipitation_forcing > 0.0) then
          if (srflag(i) > 0.0) then
            snowfall = srflag(i) * precipitation_forcing ! need snowfall for glacier snow age
          endif
        endif

        ice_flag = -1
        temperature_soil_bot = min(temperature_soil_bot,263.15)

        call noahmp_options_glacier(iopt_alb, iopt_snf, iopt_tbot, iopt_stc, iopt_gla, &
                                    iopt_sfc ,iopt_trs)
        vegetation_frac = 0.0
        call noahmp_glacier (                                                                      &
          i_location           ,1                    ,cosine_zenith        ,nsnow                , &
          nsoil                ,timestep             ,                                             &
          temperature_forcing  ,air_pressure_forcing ,uwind_forcing        ,vwind_forcing        , &
          spec_humidity_forcing,sw_radiation_forcing ,precipitation_forcing,radiation_lw_forcing , &
          temperature_soil_bot ,forcing_height       ,snow_ice_frac_old    ,zsoil                , &
          thsfc_loc            ,prslkix              ,prsik1x              ,prslk1x              , &
          air_pressure_surface ,pblhx                ,iz0tlnd              ,itime                , &
	  vegetation_frac      ,area_grid            ,psi_opt                                    , &
          con_fvirt            ,con_eps              ,con_epsm1            ,con_cp               , &
          snowfall             ,snow_water_equiv_old ,snow_albedo_old      ,                       &
          cm_noahmp            ,ch_noahmp            ,snow_levels          ,snow_water_equiv     , &
          soil_moisture_vol    ,interface_depth      ,snow_depth           ,snow_level_ice       , &
          snow_level_liquid    ,temperature_ground   ,temperature_snow_soil,soil_liquid_vol      , &
          snow_age             ,spec_humidity_surface,sw_absorbed_total    ,sw_reflected_total   , &
          lw_absorbed_total    ,sensible_heat_total  ,latent_heat_ground   ,ground_heat_total    , &
          temperature_radiative,evaporation_soil     ,runoff_surface       ,runoff_baseflow      , &
          sw_absorbed_ground   ,albedo_total         ,snowmelt_out         ,snowmelt_shallow     , &
          snowmelt_shallow_1   ,snowmelt_shallow_2   ,temperature_bare_2m  ,spec_humidity_bare_2m, &
          z0h_total                                                                              , &
          emissivity_total     ,precip_frozen_frac   ,ch_bare_ground_2m    ,snow_sublimation     , &
#ifdef CCPP
          albedo_direct        ,albedo_diffuse,       errmsg               ,errflg )
#else
          albedo_direct        ,albedo_diffuse)
#endif

#ifdef CCPP
        if (errflg /= 0) return
#endif

!
! set some non-glacier fields over the glacier
!

        snow_cover_fraction    = 1.0
        temperature_leaf       = undefined  
        canopy_ice             = undefined
        canopy_liquid          = undefined
        vapor_pres_canopy_air  = undefined
        temperature_canopy_air = undefined
        canopy_wet_fraction    = undefined
        lake_water             = undefined
        depth_water_table      = undefined
        aquifer_water          = undefined
        saturated_water        = undefined
        leaf_carbon            = undefined
        root_carbon            = undefined
        stem_carbon            = undefined
        wood_carbon            = undefined
        soil_carbon_stable     = undefined
        soil_carbon_fast       = undefined
        leaf_area_index        = undefined
        stem_area_index        = undefined
        evaporation_canopy     = undefined
        transpiration          = undefined
        aquifer_water          = undefined
        precip_adv_heat_total  = undefined
        soil_moisture_wtd      = 0.0
        recharge               = 0.0
        deep_recharge          = 0.0
        eq_soil_water_vol      = soil_moisture_vol
        transpiration_heat     = undefined
        latent_heat_canopy     = undefined
        z0_total               = 0.002
        latent_heat_total      = latent_heat_ground
        t2mmp(i)               = temperature_bare_2m
        q2mp(i)                = spec_humidity_bare_2m

        tskin(i)               = temperature_radiative
        tflux                  = temperature_ground
        surface_temperature    = temperature_ground
        vegetation_fraction    = vegetation_frac
        ch_vegetated           = 0.0
        ch_bare_ground         = ch_noahmp
        canopy_heat_storage    = 0.0
        lai_sunlit             = 0.0
        lai_shaded             = 0.0
        rs_sunlit              = 0.0
        rs_shaded              = 0.0

      else  ! not glacier

        ice_flag = 0 

        call noahmp_sflx (parameters                                          , &
          i_location            ,j_location            ,latitude              , &
          year_length           ,julian_day            ,cosine_zenith         , &
          timestep              ,spatial_scale         ,atmosphere_thickness  , &
          soil_levels           ,soil_interface_depth  ,max_snow_levels       , &
          vegetation_frac       ,max_vegetation_frac   ,vegetation_category   , &
          ice_flag              ,surface_type          ,crop_type             , &
          eq_soil_water_vol     ,temperature_forcing   ,air_pressure_forcing  , &
          air_pressure_surface  ,uwind_forcing         ,vwind_forcing         , &
          spec_humidity_forcing ,area_grid             ,cloud_water_forcing   , &
          sw_radiation_forcing  ,radiation_lw_forcing  ,thsfc_loc             , &
          prslkix               ,prsik1x               ,prslk1x               , &
          pblhx                 ,iz0tlnd               ,itime                 , &
          psi_opt                                                             , &
          precip_convective                                                   , &
          precip_non_convective ,precip_sh_convective  ,precip_snow           , &
          precip_graupel        ,precip_hail           ,temperature_soil_bot  , &
          co2_air               ,o2_air                ,foliage_nitrogen      , &
          snow_ice_frac_old     ,forcing_height                               , &
          con_fvirt             ,con_eps, con_epsm1    ,con_cp                , &
          snow_albedo_old       ,snow_water_equiv_old                         , &
          temperature_snow_soil ,soil_liquid_vol       ,soil_moisture_vol     , &
          temperature_canopy_air,vapor_pres_canopy_air ,canopy_wet_fraction   , &
          canopy_liquid         ,canopy_ice            ,temperature_leaf      , &
          temperature_ground    ,spec_humidity_surface ,snowfall              , &
          rainfall              ,snow_levels           ,interface_depth       , &
          snow_depth            ,snow_water_equiv      ,snow_level_ice        , &
          snow_level_liquid     ,depth_water_table     ,aquifer_water         , &
          saturated_water       ,                                               &
          lake_water            ,leaf_carbon           ,root_carbon           , &
          stem_carbon           ,wood_carbon           ,soil_carbon_stable    , &
          soil_carbon_fast      ,leaf_area_index       ,stem_area_index       , &
          cm_noahmp             ,ch_noahmp             ,snow_age              , &
          grain_carbon          ,growing_deg_days      ,plant_growth_stage    , &
          soil_moisture_wtd     ,deep_recharge         ,recharge,ustarx       , &
          z0_total              ,z0h_total             ,surface_temperature   , &
          sw_absorbed_total     ,sw_reflected_total                           , &
          lw_absorbed_total     ,sensible_heat_total   ,ground_heat_total     , &
          latent_heat_canopy    ,latent_heat_ground    ,transpiration_heat    , &
          evaporation_canopy    ,transpiration         ,evaporation_soil      , &
          temperature_radiative ,temperature_bare_grd  ,temperature_veg_grd   , &
          temperature_veg_2m    ,temperature_bare_2m   ,spec_humidity_veg_2m  , &
          spec_humidity_bare_2m ,runoff_surface        ,runoff_baseflow       , &
          par_absorbed          ,photosynthesis        ,sw_absorbed_veg       , &
          sw_absorbed_ground    ,snow_cover_fraction   ,net_eco_exchange      , &
          global_prim_prod      ,net_prim_prod         ,vegetation_fraction   , &
          albedo_total          ,snowmelt_out          ,snowmelt_shallow      , &
          snowmelt_shallow_1    ,snowmelt_shallow_2    ,rs_sunlit             , &
          rs_shaded             ,albedo_direct         ,albedo_diffuse        , &
          albedo_direct_snow    ,albedo_diffuse_snow                          , &
          canopy_gap_fraction                                                 , &
          incanopy_gap_fraction ,ch_vegetated          ,ch_bare_ground        , &
          emissivity_total      ,sensible_heat_grd_veg ,sensible_heat_leaf    , &
          sensible_heat_grd_bar ,latent_heat_grd_veg   ,latent_heat_grd_bare  , &
          ground_heat_veg       ,ground_heat_bare      ,lw_absorbed_grd_veg   , &
          lw_absorbed_leaf      ,lw_absorbed_grd_bare  ,latent_heat_trans     , &
          latent_heat_leaf      ,ch_leaf               ,ch_below_canopy       , &
          ch_vegetated_2m       ,ch_bare_ground_2m     ,precip_frozen_frac    , &
          precip_adv_heat_veg   ,precip_adv_heat_grd_v ,precip_adv_heat_grd_b , &
          precip_adv_heat_total ,snow_sublimation      ,canopy_heat_storage   , &
          lai_sunlit            ,lai_shaded            ,leaf_air_resistance   , &
#ifdef CCPP
          spec_humid_sfc_veg    ,spec_humid_sfc_bare   ,                        &
          errmsg                ,errflg                )
#else
          spec_humid_sfc_veg    ,spec_humid_sfc_bare   )
#endif
        
#ifdef CCPP
        if (errflg /= 0) return
#endif

        latent_heat_total  = latent_heat_canopy + latent_heat_ground + transpiration_heat

        t2mmp(i)  = temperature_veg_2m   * vegetation_fraction + &
                   temperature_bare_2m   * (1-vegetation_fraction) 
         q2mp(i)  = spec_humidity_veg_2m * vegetation_fraction + &
                   spec_humidity_bare_2m * (1-vegetation_fraction)

         tskin(i)               = temperature_radiative
         tflux                  = surface_temperature

      endif          ! glacial split ends

!
!  --- noah-mp inout and out variables
!

      tsnoxy    (i,:) = temperature_snow_soil(-2: 0)
      stc       (i,:) = temperature_snow_soil( 1:km)
      hflx      (i)   = sensible_heat_total !note unit change below
      evap      (i)   = latent_heat_total   !note unit change below
      evbs      (i)   = latent_heat_ground
      evcw      (i)   = latent_heat_canopy
      trans     (i)   = transpiration_heat
      gflux     (i)   = -1.0*ground_heat_total                          ! opposite sign to be consistent with noah
      snohf     (i)   = snowmelt_out * con_hfus         ! only snow that exits pack
      sbsno     (i)   = snow_sublimation
      pah       (i)   = precip_adv_heat_total

      cmxy      (i)   = cm_noahmp
      chxy      (i)   = ch_noahmp
      zorl      (i)   = z0_total * 100.0  ! convert to cm
      ztmax     (i)   = z0h_total 
      
      !LAI-scale canopy resistance based on weighted sunlit shaded fraction
      if(rs_sunlit .le. 0.0 .or. rs_shaded .le. 0.0 .or. &
          lai_sunlit .eq. 0.0 .or. lai_shaded .eq. 0.0) then
        rca(i) = parameters%rsmax
      else !calculate LAI-scale canopy conductance (1/Rs)
        rca(i) = ((1.0/(rs_sunlit+leaf_air_resistance)*lai_sunlit) + &
                 ((1.0/(rs_shaded+leaf_air_resistance))*lai_shaded))
        rca(i) = max((1.0/rca(i)),parameters%rsmin) !resistance
      end if
      
      smc       (i,:) = soil_moisture_vol
      slc       (i,:) = soil_liquid_vol
      snowxy    (i)   = float(snow_levels)
      weasd     (i)   = snow_water_equiv
      snicexy   (i,:) = snow_level_ice
      snliqxy   (i,:) = snow_level_liquid
      snwdph    (i)   = snow_depth * 1000.0       ! convert from mm to m
      canopy    (i)   = canopy_ice + canopy_liquid
      canliqxy  (i)   = canopy_liquid
      canicexy  (i)   = canopy_ice
      zwtxy     (i)   = depth_water_table
      waxy      (i)   = aquifer_water
      wtxy      (i)   = saturated_water
      qsnowxy   (i)   = snowfall
      ecan      (i)   = evaporation_canopy
      etran     (i)   = transpiration
      edir      (i)   = evaporation_soil
      drain     (i)   = runoff_baseflow
      runoff    (i)   = runoff_surface

      lfmassxy  (i)   = leaf_carbon
      rtmassxy  (i)   = root_carbon
      stmassxy  (i)   = stem_carbon
      woodxy    (i)   = wood_carbon
      stblcpxy  (i)   = soil_carbon_stable
      fastcpxy  (i)   = soil_carbon_fast
      xlaixy    (i)   = leaf_area_index
      xsaixy    (i)   = stem_area_index

      snowc     (i)   = snow_cover_fraction
      sncovr1   (i)   = snow_cover_fraction

      qsurf     (i)   = spec_humidity_surface
      tsurf     (i)   = tskin(i)

      tvxy      (i)   = temperature_leaf
      tgxy      (i)   = temperature_ground
      tahxy     (i)   = temperature_canopy_air
      eahxy     (i)   = vapor_pres_canopy_air
      emiss     (i)   = emissivity_total

      if(albedo_total > 0.0) then
        sfalb(i)      = albedo_total
        albdvis(i)    = albedo_direct(1)
        albdnir(i)    = albedo_direct(2)
        albivis(i)    = albedo_diffuse(1)
        albinir(i)    = albedo_diffuse(2)
      end if

      zsnsoxy   (i,:) = interface_depth

      if(present(canopy_heat_storage_ccpp  )) canopy_heat_storage_ccpp  (i)   = canopy_heat_storage
      if(present(rainfall_ccpp             )) rainfall_ccpp             (i)   = rainfall
      if(present(sw_absorbed_total_ccpp    )) sw_absorbed_total_ccpp    (i)   = sw_absorbed_total
      if(present(sw_reflected_total_ccpp   )) sw_reflected_total_ccpp   (i)   = sw_reflected_total
      if(present(lw_absorbed_total_ccpp    )) lw_absorbed_total_ccpp    (i)   = lw_absorbed_total
      if(present(temperature_bare_grd_ccpp )) temperature_bare_grd_ccpp (i)   = temperature_bare_grd
      if(present(temperature_veg_grd_ccpp  )) temperature_veg_grd_ccpp  (i)   = temperature_veg_grd
      if(present(temperature_veg_2m_ccpp   )) temperature_veg_2m_ccpp   (i)   = temperature_veg_2m
      if(present(temperature_bare_2m_ccpp  )) temperature_bare_2m_ccpp  (i)   = temperature_bare_2m
      if(present(spec_humidity_veg_2m_ccpp )) spec_humidity_veg_2m_ccpp (i)   = spec_humidity_veg_2m
      if(present(spec_humidity_bare_2m_ccpp)) spec_humidity_bare_2m_ccpp(i)   = spec_humidity_bare_2m
      if(present(sw_absorbed_veg_ccpp      )) sw_absorbed_veg_ccpp      (i)   = sw_absorbed_veg
      if(present(sw_absorbed_ground_ccpp   )) sw_absorbed_ground_ccpp   (i)   = sw_absorbed_ground
      if(present(snowmelt_out_ccpp         )) snowmelt_out_ccpp         (i)   = snowmelt_out
      if(present(snowmelt_shallow_ccpp     )) snowmelt_shallow_ccpp     (i)   = snowmelt_shallow
      if(present(albedo_direct_snow_ccpp   )) albedo_direct_snow_ccpp   (i,:) = albedo_direct_snow
      if(present(albedo_diffuse_snow_ccpp  )) albedo_diffuse_snow_ccpp  (i,:) = albedo_diffuse_snow
      if(present(ch_vegetated_ccpp         )) ch_vegetated_ccpp         (i)   = ch_vegetated
      if(present(ch_bare_ground_ccpp       )) ch_bare_ground_ccpp       (i)   = ch_bare_ground
      if(present(sensible_heat_grd_veg_ccpp)) sensible_heat_grd_veg_ccpp(i)   = sensible_heat_grd_veg
      if(present(sensible_heat_leaf_ccpp   )) sensible_heat_leaf_ccpp   (i)   = sensible_heat_leaf
      if(present(sensible_heat_grd_bar_ccpp)) sensible_heat_grd_bar_ccpp(i)   = sensible_heat_grd_bar
      if(present(latent_heat_grd_veg_ccpp  )) latent_heat_grd_veg_ccpp  (i)   = latent_heat_grd_veg
      if(present(latent_heat_grd_bare_ccpp )) latent_heat_grd_bare_ccpp (i)   = latent_heat_grd_bare
      if(present(ground_heat_veg_ccpp      )) ground_heat_veg_ccpp      (i)   = ground_heat_veg
      if(present(ground_heat_bare_ccpp     )) ground_heat_bare_ccpp     (i)   = ground_heat_bare
      if(present(lw_absorbed_grd_veg_ccpp  )) lw_absorbed_grd_veg_ccpp  (i)   = lw_absorbed_grd_veg
      if(present(lw_absorbed_leaf_ccpp     )) lw_absorbed_leaf_ccpp     (i)   = lw_absorbed_leaf
      if(present(lw_absorbed_grd_bare_ccpp )) lw_absorbed_grd_bare_ccpp (i)   = lw_absorbed_grd_bare
      if(present(latent_heat_trans_ccpp    )) latent_heat_trans_ccpp    (i)   = latent_heat_trans
      if(present(latent_heat_leaf_ccpp     )) latent_heat_leaf_ccpp     (i)   = latent_heat_leaf
      if(present(ch_leaf_ccpp              )) ch_leaf_ccpp              (i)   = ch_leaf
      if(present(ch_below_canopy_ccpp      )) ch_below_canopy_ccpp      (i)   = ch_below_canopy
      if(present(ch_vegetated_2m_ccpp      )) ch_vegetated_2m_ccpp      (i)   = ch_vegetated_2m
      if(present(ch_bare_ground_2m_ccpp    )) ch_bare_ground_2m_ccpp    (i)   = ch_bare_ground_2m
      if(present(precip_adv_heat_veg_ccpp  )) precip_adv_heat_veg_ccpp  (i)   = precip_adv_heat_veg
      if(present(precip_adv_heat_grd_v_ccpp)) precip_adv_heat_grd_v_ccpp(i)   = precip_adv_heat_grd_v
      if(present(precip_adv_heat_grd_b_ccpp)) precip_adv_heat_grd_b_ccpp(i)   = precip_adv_heat_grd_b
      if(present(spec_humid_sfc_veg_ccpp   )) spec_humid_sfc_veg_ccpp   (i)   = spec_humid_sfc_veg
      if(present(spec_humid_sfc_bare_ccpp  )) spec_humid_sfc_bare_ccpp  (i)   = spec_humid_sfc_bare

      wslakexy  (i)   = lake_water         ! not active
      fwetxy    (i)   = canopy_wet_fraction
      taussxy   (i)   = snow_age
      alboldxy  (i)   = snow_albedo_old
      sneqvoxy  (i)   = snow_water_equiv_old

      smcwtdxy  (i)   = soil_moisture_wtd  ! only need for run=5
      deeprechxy(i)   = deep_recharge      ! only need for run=5
      rechxy    (i)   = recharge           ! only need for run=5
      smoiseq   (i,:) = eq_soil_water_vol  ! only need for run=5; listed as in

      stm       (i)   = (0.1*soil_moisture_vol(1) + &
                         0.3*soil_moisture_vol(2) + &
                         0.6*soil_moisture_vol(3) + &        ! clean up and use depths above
                         1.0*soil_moisture_vol(4))*1000.0    ! unit conversion from m to kg m-2

      wet1   (i) = soil_moisture_vol(1) /  smcmax_table(soil_category(1))
      smcwlt2(i) = smcdry_table(soil_category(1))   !!!change to wilt?
      smcref2(i) = smcref_table(soil_category(1))

      virtfac1  = 1.0 +  con_fvirt * max(q1(i), 1.e-8)          !from forcing

      if(thsfc_loc) then ! Use local potential temperature
            vptemp    =temperature_forcing * prslki(i)*virtfac1       !virtual potential temperature @zlvl 1
        else ! Use potential temperature reference to 1000 hPa
            vptemp    =temperature_forcing /prslk1(i) * virtfac1
       endif

       if(thsfc_loc) then ! Use local potential temperature
              tvs1   = tflux * virtfac1
         else ! Use potential temperature referenced to 1000 hPa
              tvs1   = tflux/prsik1(i) * virtfac1
       endif

      z0_total  = max(min(z0_total,forcing_height),1.0e-6)
      z0h_total = max(z0h_total,1.0e-6)


            tem1 = (z0_total - z0lo) / (z0up - z0lo)
            tem1 = min(max(tem1, 0.0_kind_phys), 1.0_kind_phys)
            tem2 = max(vegetation_fraction, 0.1_kind_phys)
            zvfun(i) = sqrt(tem1 * tem2)
            gdx=sqrt(garea(i))

!      if ( .not. do_mynnsfclay) then   !GFS sfcdiff
       if ( iopt_sfc .ne. 4 ) then   !GFS sfcdiff

      call       gfs_stability                                                               &
        (zf(i), zvfun(i), gdx, virtual_temperature, vptemp,wind(i), z0_total, z0h_total, & 
         tvs1, con_g, thsfc_loc,                                                         &
         rb1(i), fm1(i), fh1(i), fm101(i), fh21(i), cm(i), ch(i), stress1(i), ustar1(i))

       rmol1(i) = undefined  !not used in GFS sfcdif -> to satsify output
       flhc1(i) = undefined
       flqc1(i) = undefined

        rah_total = max(1.0,1.0/( ch(i)*wind(i)) )
        cah_total = density * con_cp /rah_total
!       tskin(i) = sensible_heat_total/cah_total + temperature_forcing ! test to use combined ch and SH to backout Ts

!        ch(i) = ch_vegetated * vegetation_frac + ch_bare_ground*(1.0-vegetation_frac)

      else    ! MYNN - note the GFS option is the same as sfcdif3; so removed.

             qfx = evap(i) / con_hvap         ! use flux from output

           call sfcdif4(i_location  ,j_location  ,uwind_forcing ,vwind_forcing ,           &
                        temperature_forcing, air_pressure_forcing ,air_pressure_surface  , &
                        pblhx,gdx,z0_total,con_fvirt,con_eps,con_cp,itime,snwdph(i),mnice, &
                        psi_opt,surface_temperature,                                       &
                        spec_humidity_forcing,forcing_height,iz0tlnd,spec_humidity_surface,&
                        sensible_heat_total,qfx,cm(i),ch(i),ch2,cq2,rmol1(i),ustar1(i),    &
                        rb1(i),fm1(i),fh1(i),stress1(i),fm101(i),fh21(i),wspd1,flhc1(i),   &
                        flqc1(i) )

              ch(i)=ch(i)/wspd1
              cm(i)=cm(i)/wspd1

              ch(i) = ch_vegetated * vegetation_fraction + ch_bare_ground*(1.0-vegetation_fraction)

          rah_total = max(1.0,1.0/( ch(i)*wind(i)) )
          cah_total = density * con_cp /rah_total

!          tskin(i) = sensible_heat_total/cah_total + temperature_forcing !

       endif



      cmxy(i) = cm(i)
      chxy(i) = ch(i)

      chh       (i)   = chxy(i)  * wind(i) * density
      cmm       (i)   = cmxy(i)  * wind(i)

      snwdph    (i)   = snow_depth * 1000.0       ! convert from m to mm; wait after the stability call
!     qsurf     (i)   = q1(i) + evap(i)/(con_hvap*density*ch(i)*wind(i))

!      
!  --- change units for output
!
      hflx(i) = hflx(i) / density / con_cp
      evap(i) = evap(i) / density / con_hvap

!
!  --- calculate potential evaporation using noah code
!
      potential_temperature = temperature_forcing * prslki(i)
      virtual_temperature   = temperature_forcing * (1.0 + 0.61*spec_humidity_forcing)
      penman_radiation      = sw_absorbed_total + radiation_lw_forcing
      dqsdt                 = spec_humidity_sat * a23m4/(temperature_forcing-a4)**2

      precip_freeze_frac_in = srflag(i)
      is_snowing            = .false.          
      is_freeze_rain        = .false.
      if (precipitation_forcing > 0.0) then
        if (precip_freeze_frac_in > 0.0) then                   ! rain/snow flag, one condition is enough?
          is_snowing = .true.
        else
          if (temperature_forcing <= 275.15) is_freeze_rain = .true.  
        end if
      end if

!
! using new combined ch output to compute ep
!
      ch_noahmp = chxy(i) * wind(i)
      
      call penman (temperature_forcing, air_pressure_forcing , ch_noahmp            , &
                   virtual_temperature, potential_temperature, precipitation_forcing, &
                   penman_radiation   , ground_heat_total    , spec_humidity_forcing, &
                   spec_humidity_sat  , potential_evaporation, is_snowing           , &
                   is_freeze_rain     , precip_freeze_frac_in, dqsdt                , &
                   emissivity_total   , snow_cover_fraction  )

      ep(i) = potential_evaporation

    end if ! flag_iter(i) .and. dry(i)

  end do ! im loop

      return

      end subroutine noahmpdrv_run
!> @}
!-----------------------------------

!> \ingroup NoahMP_LSM
!! \brief This subroutine fills in a derived data type of type noahmp_parameters with data
!! from the module \ref noahmp_tables.
      subroutine transfer_mp_parameters (vegtype,soiltype,slopetype,    &
                                       soilcolor,croptype,parameters)
     
        use noahmp_tables
        use module_sf_noahmplsm
      
        implicit none
      
        integer, intent(in)    :: vegtype
        integer, intent(in)    :: soiltype(4)
        integer, intent(in)    :: slopetype
        integer, intent(in)    :: soilcolor
        integer, intent(in)    :: croptype
          
        type (noahmp_parameters), intent(out) :: parameters
          
        real    :: refdk
        real    :: refkdt
        real    :: frzk
        real    :: frzfact
        integer :: isoil
      
        parameters%iswater   =  iswater_table
        parameters%isbarren  =  isbarren_table
        parameters%isice     =  isice_table
        parameters%iscrop    =  iscrop_table
        parameters%eblforest =  eblforest_table
      
!-----------------------------------------------------------------------&
        parameters%urban_flag = .false.
        if( vegtype == isurban_table .or. vegtype == 31                 &
     &         .or.vegtype  == 32 .or. vegtype == 33) then
           parameters%urban_flag = .true.
        endif
      
!------------------------------------------------------------------------------------------!
! transfer veg parameters
!------------------------------------------------------------------------------------------!
      
        parameters%ch2op  =  ch2op_table(vegtype)       !maximum intercepted h2o per unit lai+sai (mm)
        parameters%dleaf  =  dleaf_table(vegtype)       !characteristic leaf dimension (m)
        parameters%z0mvt  =  z0mvt_table(vegtype)       !momentum roughness length (m)
        parameters%hvt    =    hvt_table(vegtype)       !top of canopy (m)
        parameters%hvb    =    hvb_table(vegtype)       !bottom of canopy (m)
        parameters%z0mhvt = z0mhvt_table(vegtype)       !momentum roughness length (m)
        parameters%den    =    den_table(vegtype)       !tree density (no. of trunks per m2)
        parameters%rc     =     rc_table(vegtype)       !tree crown radius (m)
        parameters%mfsno  =  mfsno_table(vegtype)       !snowmelt m parameter ()
	parameters%scffac = scffac_table(vegtype)       !snow cover factor
        parameters%cbiom  =  cbiom_table(vegtype)       !canopy biomass heat capacity parameter (m)
        parameters%saim   =   saim_table(vegtype,:)     !monthly stem area index, one-sided
        parameters%laim   =   laim_table(vegtype,:)     !monthly leaf area index, one-sided
        parameters%sla    =    sla_table(vegtype)       !single-side leaf area per kg [m2/kg]
        parameters%dilefc = dilefc_table(vegtype)       !coeficient for leaf stress death [1/s]
        parameters%dilefw = dilefw_table(vegtype)       !coeficient for leaf stress death [1/s]
        parameters%fragr  =  fragr_table(vegtype)       !fraction of growth respiration  !original was 0.3 
        parameters%ltovrc = ltovrc_table(vegtype)       !leaf turnover [1/s]
      
        parameters%c3psn  =  c3psn_table(vegtype)       !photosynthetic pathway: 0. = c4, 1. = c3
        parameters%kc25   =   kc25_table(vegtype)       !co2 michaelis-menten constant at 25c (pa)
        parameters%akc    =    akc_table(vegtype)       !q10 for kc25
        parameters%ko25   =   ko25_table(vegtype)       !o2 michaelis-menten constant at 25c (pa)
        parameters%ako    =    ako_table(vegtype)       !q10 for ko25
        parameters%vcmx25 = vcmx25_table(vegtype)       !maximum rate of carboxylation at 25c (umol co2/m**2/s)
        parameters%avcmx  =  avcmx_table(vegtype)       !q10 for vcmx25
        parameters%bp     =     bp_table(vegtype)       !minimum leaf conductance (umol/m**2/s)
        parameters%mp     =     mp_table(vegtype)       !slope of conductance-to-photosynthesis relationship
        parameters%qe25   =   qe25_table(vegtype)       !quantum efficiency at 25c (umol co2 / umol photon)
        parameters%aqe    =    aqe_table(vegtype)       !q10 for qe25
        parameters%rmf25  =  rmf25_table(vegtype)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
        parameters%rms25  =  rms25_table(vegtype)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
        parameters%rmr25  =  rmr25_table(vegtype)       !root maintenance respiration at 25c (umol co2/kg bio/s)
        parameters%arm    =    arm_table(vegtype)       !q10 for maintenance respiration
        parameters%folnmx = folnmx_table(vegtype)       !foliage nitrogen concentration when f(n)=1 (%)
        parameters%tmin   =   tmin_table(vegtype)       !minimum temperature for photosynthesis (k)
      
        parameters%xl     =     xl_table(vegtype)       !leaf/stem orientation index
        parameters%rhol   =   rhol_table(vegtype,:)     !leaf reflectance: 1=vis, 2=nir
        parameters%rhos   =   rhos_table(vegtype,:)     !stem reflectance: 1=vis, 2=nir
        parameters%taul   =   taul_table(vegtype,:)     !leaf transmittance: 1=vis, 2=nir
        parameters%taus   =   taus_table(vegtype,:)     !stem transmittance: 1=vis, 2=nir
      
        parameters%mrp    =    mrp_table(vegtype)       !microbial respiration parameter (umol co2 /kg c/ s)
        parameters%cwpvt  =  cwpvt_table(vegtype)       !empirical canopy wind parameter
      
        parameters%wrrat  =  wrrat_table(vegtype)       !wood to non-wood ratio
        parameters%wdpool = wdpool_table(vegtype)       !wood pool (switch 1 or 0) depending on woody or not [-]
        parameters%tdlef  =  tdlef_table(vegtype)       !characteristic t for leaf freezing [k]
      
        parameters%nroot  =  nroot_table(vegtype)       !number of soil layers with root present
        parameters%rgl    =    rgl_table(vegtype)       !parameter used in radiation stress function
        parameters%rsmin  =     rs_table(vegtype)       !minimum stomatal resistance [s m-1]
        parameters%hs     =     hs_table(vegtype)       !parameter used in vapor pressure deficit function
        parameters%topt   =   topt_table(vegtype)       !optimum transpiration air temperature [k]
        parameters%rsmax  =  rsmax_table(vegtype)       !maximal stomatal resistance [s m-1]
      
!------------------------------------------------------------------------------------------!
! transfer rad parameters
!------------------------------------------------------------------------------------------!
      
         parameters%albsat    = albsat_table(soilcolor,:)
         parameters%albdry    = albdry_table(soilcolor,:)
         parameters%albice    = albice_table
         parameters%alblak    = alblak_table               
         parameters%omegas    = omegas_table
         parameters%betads    = betads_table
         parameters%betais    = betais_table
         parameters%eg        = eg_table
      
!------------------------------------------------------------------------------------------!
! Transfer crop parameters
!------------------------------------------------------------------------------------------!

      if(croptype > 0) then
        parameters%pltday    =    pltday_table(croptype)    ! planting date
        parameters%hsday     =     hsday_table(croptype)    ! harvest date
        parameters%plantpop  =  plantpop_table(croptype)    ! plant density [per ha] - used?
        parameters%irri      =      irri_table(croptype)    ! irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
        parameters%gddtbase  =  gddtbase_table(croptype)    ! base temperature for gdd accumulation [c]
        parameters%gddtcut   =   gddtcut_table(croptype)    ! upper temperature for gdd accumulation [c]
        parameters%gdds1     =     gdds1_table(croptype)    ! gdd from seeding to emergence
        parameters%gdds2     =     gdds2_table(croptype)    ! gdd from seeding to initial vegetative 
        parameters%gdds3     =     gdds3_table(croptype)    ! gdd from seeding to post vegetative 
        parameters%gdds4     =     gdds4_table(croptype)    ! gdd from seeding to intial reproductive
        parameters%gdds5     =     gdds5_table(croptype)    ! gdd from seeding to pysical maturity 
        parameters%c3c4      =      c3c4_table(croptype)    ! photosynthetic pathway:  1. = c3 2. = c4
        parameters%aref      =      aref_table(croptype)    ! reference maximum co2 assimulation rate 
        parameters%psnrf     =     psnrf_table(croptype)    ! co2 assimulation reduction factor(0-1) (e.g.pests, weeds)
        parameters%i2par     =     i2par_table(croptype)    ! fraction of incoming solar radiation to photosynthetically active radiation
        parameters%tassim0   =   tassim0_table(croptype)    ! minimum temperature for co2 assimulation [c]
        parameters%tassim1   =   tassim1_table(croptype)    ! co2 assimulation linearly increasing until temperature reaches t1 [c]
        parameters%tassim2   =   tassim2_table(croptype)    ! co2 assmilation rate remain at aref until temperature reaches t2 [c]
        parameters%k         =         k_table(croptype)    ! light extinction coefficient
        parameters%epsi      =      epsi_table(croptype)    ! initial light use efficiency
        parameters%q10mr     =     q10mr_table(croptype)    ! q10 for maintainance respiration
        parameters%foln_mx   =   foln_mx_table(croptype)    ! foliage nitrogen concentration when f(n)=1 (%)
        parameters%lefreez   =   lefreez_table(croptype)    ! characteristic t for leaf freezing [k]
        parameters%dile_fc   =   dile_fc_table(croptype,:)  ! coeficient for temperature leaf stress death [1/s]
        parameters%dile_fw   =   dile_fw_table(croptype,:)  ! coeficient for water leaf stress death [1/s]
        parameters%fra_gr    =    fra_gr_table(croptype)    ! fraction of growth respiration
        parameters%lf_ovrc   =   lf_ovrc_table(croptype,:)  ! fraction of leaf turnover  [1/s]
        parameters%st_ovrc   =   st_ovrc_table(croptype,:)  ! fraction of stem turnover  [1/s]
        parameters%rt_ovrc   =   rt_ovrc_table(croptype,:)  ! fraction of root tunrover  [1/s]
        parameters%lfmr25    =    lfmr25_table(croptype)    ! leaf maintenance respiration at 25c [umol co2/m**2  /s]
        parameters%stmr25    =    stmr25_table(croptype)    ! stem maintenance respiration at 25c [umol co2/kg bio/s]
        parameters%rtmr25    =    rtmr25_table(croptype)    ! root maintenance respiration at 25c [umol co2/kg bio/s]
        parameters%grainmr25 = grainmr25_table(croptype)    ! grain maintenance respiration at 25c [umol co2/kg bio/s]
        parameters%lfpt      =      lfpt_table(croptype,:)  ! fraction of carbohydrate flux to leaf
        parameters%stpt      =      stpt_table(croptype,:)  ! fraction of carbohydrate flux to stem
        parameters%rtpt      =      rtpt_table(croptype,:)  ! fraction of carbohydrate flux to root
        parameters%grainpt   =   grainpt_table(croptype,:)  ! fraction of carbohydrate flux to grain
        parameters%bio2lai   =   bio2lai_table(croptype)    ! leaf are per living leaf biomass [m^2/kg]
      end if

!------------------------------------------------------------------------------------------!
! transfer global parameters
!------------------------------------------------------------------------------------------!
      
         parameters%co2          =          co2_table
         parameters%o2           =           o2_table
         parameters%timean       =       timean_table
         parameters%fsatmx       =       fsatmx_table
         parameters%z0sno        =        z0sno_table
         parameters%ssi          =          ssi_table
         parameters%snow_ret_fac = snow_ret_fac_table
         parameters%swemx        =        swemx_table
         parameters%tau0         =         tau0_table
         parameters%grain_growth = grain_growth_table
         parameters%extra_growth = extra_growth_table
         parameters%dirt_soot    =    dirt_soot_table
         parameters%bats_cosz    =    bats_cosz_table
         parameters%bats_vis_new = bats_vis_new_table
         parameters%bats_nir_new = bats_nir_new_table
         parameters%bats_vis_age = bats_vis_age_table
         parameters%bats_nir_age = bats_nir_age_table
         parameters%bats_vis_dir = bats_vis_dir_table
         parameters%bats_nir_dir = bats_nir_dir_table
         parameters%rsurf_snow   =   rsurf_snow_table
         parameters%rsurf_exp    =    rsurf_exp_table
         parameters%snow_emis    =    snow_emis_table
      
! ----------------------------------------------------------------------
!  transfer soil parameters
! ----------------------------------------------------------------------
      
      do isoil = 1, size(soiltype)
        parameters%bexp(isoil)   = bexp_table   (soiltype(isoil))
        parameters%dksat(isoil)  = dksat_table  (soiltype(isoil))
        parameters%dwsat(isoil)  = dwsat_table  (soiltype(isoil))
        parameters%psisat(isoil) = psisat_table (soiltype(isoil))
        parameters%quartz(isoil) = quartz_table (soiltype(isoil))
        parameters%smcdry(isoil) = smcdry_table (soiltype(isoil))
        parameters%smcmax(isoil) = smcmax_table (soiltype(isoil))
        parameters%smcref(isoil) = smcref_table (soiltype(isoil))
        parameters%smcwlt(isoil) = smcwlt_table (soiltype(isoil))
      end do
          
      parameters%f1     = f1_table(soiltype(1))
      parameters%refdk  = refdk_table
      parameters%refkdt = refkdt_table

! ----------------------------------------------------------------------
! transfer genparm parameters
! ----------------------------------------------------------------------
          parameters%csoil  = csoil_table
          parameters%zbot   = zbot_table
          parameters%czil   = czil_table
      
          frzk   = frzk_table
          parameters%kdt    = parameters%refkdt * parameters%dksat(1) / parameters%refdk
          parameters%slope  = slope_table(slopetype)
      
          if(parameters%urban_flag)then  ! hardcoding some urban parameters for soil
             parameters%smcmax = 0.45 
             parameters%smcref = 0.42 
             parameters%smcwlt = 0.40 
             parameters%smcdry = 0.40 
             parameters%csoil  = 3.e6
          endif
      
      ! adjust frzk parameter to actual soil type: frzk * frzfact
      
!-----------------------------------------------------------------------&
          if(soiltype(1) /= 14) then
            frzfact = (parameters%smcmax(1) / parameters%smcref(1)) * (0.412 / 0.468)
            parameters%frzx = frzk * frzfact
          end if
      
       end subroutine transfer_mp_parameters

!> \ingroup NoahMP_LSM
!! \brief This subroutine uses a pedotransfer method to calculate soil properties.
SUBROUTINE PEDOTRANSFER_SR2006(nsoil,sand,clay,orgm,parameters)

  use module_sf_noahmplsm
  use noahmp_tables
        
  implicit none
        
  integer,                    intent(in   ) :: nsoil     ! number of soil layers
  real, dimension( 1:nsoil ), intent(inout) :: sand
  real, dimension( 1:nsoil ), intent(inout) :: clay
  real, dimension( 1:nsoil ), intent(inout) :: orgm
    
  real, dimension( 1:nsoil ) :: theta_1500t
  real, dimension( 1:nsoil ) :: theta_1500
  real, dimension( 1:nsoil ) :: theta_33t
  real, dimension( 1:nsoil ) :: theta_33
  real, dimension( 1:nsoil ) :: theta_s33t
  real, dimension( 1:nsoil ) :: theta_s33
  real, dimension( 1:nsoil ) :: psi_et
  real, dimension( 1:nsoil ) :: psi_e
    
  type(noahmp_parameters), intent(inout) :: parameters
  integer :: k

  do k = 1,4
    if(sand(k) <= 0 .or. clay(k) <= 0) then
      sand(k) = 0.41
      clay(k) = 0.18
    end if
    if(orgm(k) <= 0 ) orgm(k) = 0.0
  end do
        
  theta_1500t =   sr2006_theta_1500t_a*sand       &
                + sr2006_theta_1500t_b*clay       &
                + sr2006_theta_1500t_c*orgm       &
                + sr2006_theta_1500t_d*sand*orgm  &
                + sr2006_theta_1500t_e*clay*orgm  &
                + sr2006_theta_1500t_f*sand*clay  &
                + sr2006_theta_1500t_g

  theta_1500  =   theta_1500t                      &
                + sr2006_theta_1500_a*theta_1500t  &
                + sr2006_theta_1500_b

  theta_33t   =   sr2006_theta_33t_a*sand       &
                + sr2006_theta_33t_b*clay       &
                + sr2006_theta_33t_c*orgm       &
                + sr2006_theta_33t_d*sand*orgm  &
                + sr2006_theta_33t_e*clay*orgm  &
                + sr2006_theta_33t_f*sand*clay  &
                + sr2006_theta_33t_g

  theta_33    =   theta_33t                              &
                + sr2006_theta_33_a*theta_33t*theta_33t  &
                + sr2006_theta_33_b*theta_33t            &
                + sr2006_theta_33_c

  theta_s33t  =   sr2006_theta_s33t_a*sand      &
                + sr2006_theta_s33t_b*clay      &
                + sr2006_theta_s33t_c*orgm      &
                + sr2006_theta_s33t_d*sand*orgm &
                + sr2006_theta_s33t_e*clay*orgm &
                + sr2006_theta_s33t_f*sand*clay &
                + sr2006_theta_s33t_g

  theta_s33   = theta_s33t                       &
                + sr2006_theta_s33_a*theta_s33t  &
                + sr2006_theta_s33_b

  psi_et      =   sr2006_psi_et_a*sand           &
                + sr2006_psi_et_b*clay           &
                + sr2006_psi_et_c*theta_s33      &
                + sr2006_psi_et_d*sand*theta_s33 &
                + sr2006_psi_et_e*clay*theta_s33 &
                + sr2006_psi_et_f*sand*clay      &
                + sr2006_psi_et_g
 
  psi_e       =   psi_et                        &
                + sr2006_psi_e_a*psi_et*psi_et  &
                + sr2006_psi_e_b*psi_et         &
                + sr2006_psi_e_c
    
  parameters%smcwlt = theta_1500
  parameters%smcref = theta_33
  parameters%smcmax =   theta_33    &
                      + theta_s33            &
                      + sr2006_smcmax_a*sand &
                      + sr2006_smcmax_b

  parameters%bexp   = 3.816712826 / (log(theta_33) - log(theta_1500) )
  parameters%psisat = psi_e
  parameters%dksat  = 1930.0 * (parameters%smcmax - theta_33) ** (3.0 - 1.0/parameters%bexp)
  parameters%quartz = sand
    
! Units conversion
    
  parameters%psisat = max(0.1,parameters%psisat)     ! arbitrarily impose a limit of 0.1kpa
  parameters%psisat = 0.101997 * parameters%psisat   ! convert kpa to m
  parameters%dksat  = parameters%dksat / 3600000.0   ! convert mm/h to m/s
  parameters%dwsat  = parameters%dksat * parameters%psisat *parameters%bexp / parameters%smcmax  ! units should be m*m/s
  parameters%smcdry = parameters%smcwlt
  
! Introducing somewhat arbitrary limits (based on SOILPARM) to prevent bad things
  
  parameters%smcmax = max(0.32 ,min(parameters%smcmax,             0.50 ))
  parameters%smcref = max(0.17 ,min(parameters%smcref,parameters%smcmax ))
  parameters%smcwlt = max(0.01 ,min(parameters%smcwlt,parameters%smcref ))
  parameters%smcdry = max(0.01 ,min(parameters%smcdry,parameters%smcref ))
  parameters%bexp   = max(2.50 ,min(parameters%bexp,               12.0 ))
  parameters%psisat = max(0.03 ,min(parameters%psisat,             1.00 ))
  parameters%dksat  = max(5.e-7,min(parameters%dksat,              1.e-5))
  parameters%dwsat  = max(1.e-6,min(parameters%dwsat,              3.e-5))
  parameters%quartz = max(0.05 ,min(parameters%quartz,             0.95 ))
    
 END SUBROUTINE PEDOTRANSFER_SR2006

!-----------------------------------------------------------------------&

!> \ingroup NoahMP_LSM
!! brief Calculate potential evaporation for the current point. Various
!! partial sums/products are also calculated and passed back to the
!! calling routine for later use.
      subroutine penman (sfctmp,sfcprs,ch,t2v,th2,prcp,fdown,ssoil,     &
     &                   q2,q2sat,etp,snowng,frzgra,ffrozp,             &
     &                   dqsdt2,emissi_in,sncovr)
 
! etp is calcuated right after ssoil

! ----------------------------------------------------------------------
! subroutine penman
! ----------------------------------------------------------------------
      use machine, only: kind_phys
      implicit none
      logical, intent(in)     :: snowng, frzgra
      real(kind=kind_phys), intent(in)        :: ch, dqsdt2,fdown,prcp,ffrozp,          &
     &                           q2, q2sat,ssoil, sfcprs, sfctmp,       &
     &                           t2v, th2,emissi_in,sncovr
      real(kind=kind_phys), intent(out)       :: etp
      real(kind=kind_phys)                    :: epsca,flx2,rch,rr,t24
      real(kind=kind_phys)                    :: a, delta, fnet,rad,rho,emissi,elcp1,lvs

      real(kind=kind_phys), parameter :: elcp = 2.4888e+3, lsubc = 2.501000e+6,cp = 1004.6
      real(kind=kind_phys), parameter :: lsubs = 2.83e+6, rd = 287.05, cph2o = 4.1855e+3
      real(kind=kind_phys), parameter :: cpice = 2.106e+3, lsubf   = 3.335e5  
      real(kind=kind_phys), parameter :: sigma = 5.6704e-8

! ----------------------------------------------------------------------
! executable code begins here:
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! prepare partial quantities for penman equation.
! ----------------------------------------------------------------------
        emissi=emissi_in
!       elcp1  = (1.0-sncovr)*elcp  + sncovr*elcp*lsubs/lsubc
        lvs    = (1.0-sncovr)*lsubc + sncovr*lsubs

      flx2 = 0.0
      delta = elcp * dqsdt2
!     delta = elcp1 * dqsdt2
      t24 = sfctmp * sfctmp * sfctmp * sfctmp
       rr = t24 * 6.48e-8 / (sfcprs * ch) + 1.0
!     rr = emissi*t24 * 6.48e-8 / (sfcprs * ch) + 1.0
      rho = sfcprs / (rd * t2v)

! ----------------------------------------------------------------------
! adjust the partial sums / products with the latent heat
! effects caused by falling precipitation.
! ----------------------------------------------------------------------
      rch = rho * cp * ch
      if (.not. snowng) then
         if (prcp >  0.0) rr = rr + cph2o * prcp / rch
      else
! ---- ...  fractional snowfall/rainfall
        rr = rr + (cpice*ffrozp+cph2o*(1.-ffrozp))                      &
     &       *prcp/rch
      end if

! ----------------------------------------------------------------------
! include the latent heat effects of frzng rain converting to ice on
! impact in the calculation of flx2 and fnet.
! ----------------------------------------------------------------------
!      fnet = fdown - sigma * t24- ssoil
      fnet = fdown -  emissi*sigma * t24- ssoil
      if (frzgra) then
         flx2 = - lsubf * prcp
         fnet = fnet - flx2
! ----------------------------------------------------------------------
! finish penman equation calculations.
! ----------------------------------------------------------------------
      end if
      rad = fnet / rch + th2- sfctmp
       a = elcp * (q2sat - q2)
!     a = elcp1 * (q2sat - q2)
      epsca = (a * rr + rad * delta) / (delta + rr)
       etp = epsca * rch / lsubc
!     etp = epsca * rch / lvs

! ----------------------------------------------------------------------
      end subroutine penman

      end module noahmpdrv
