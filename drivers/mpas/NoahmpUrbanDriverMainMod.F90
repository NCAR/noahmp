module NoahmpUrbanDriverMainMod

!--------------------------------------------------------------------------
!  SUBROUTINE noahmp_urban seperated from module_sf_noahmpdrv.F
!  P. Valayamkunnath C. He & refactor team (April 08 2022)
!  The entire urban code is not refactored.
!  Sep 24, 2026: adapted to MPAS compatibility, Cenlin He (NCAR)
!--------------------------------------------------------------------------

  use module_sf_urban,        only: urban
! use module_sf_bep,          only: bep
! use module_sf_bep_bem,      only: bep_bem
  use mpas_atmphys_date_time, only: cal_mon_day
  use mpas_kind_types
  use mpas_constants,         only: cp
  use mpas_atmphys_constants, only: karman,xlv,rcp
  use NoahmpIOVarType,        only: NoahmpIO_type
 
  implicit none
  
contains

  subroutine noahmp_urban(NoahmpIO, sf_urban_physics,  NSOIL, IVGTYP,  ITIMESTEP,           & ! IN : Model configuration 
                                 DT,     COSZ_URB2D,     XLAT_URB2D,                        & ! IN : Tite/Space-related
                                T3D,           QV3D,          U_PHY,      V_PHY,   SWDOWN,  & ! IN : Forcing
                             SWDDIR,         SWDDIF,                                        &
		                GLW,          P8W3D,         RAINBL,       DZ8W,      ZNT,  & ! IN : Forcing
                                TSK,            HFX,            QFX,         LH,   GRDFLX,  & ! IN/OUT : LSM 
		             ALBEDO,          EMISS,           QSFC,                        & ! IN/OUT : LSM 
                            its,ite,        kts,kte,                                        &
                         cmr_sfcdif,     chr_sfcdif,     cmc_sfcdif,                        &
	                 chc_sfcdif,    cmgr_sfcdif,    chgr_sfcdif,                        &
                           tr_urb2d,       tb_urb2d,       tg_urb2d,                        & !H urban
	                   tc_urb2d,       qc_urb2d,       uc_urb2d,                        & !H urban
                         xxxr_urb2d,     xxxb_urb2d,     xxxg_urb2d, xxxc_urb2d,            & !H urban
                          trl_urb3d,      tbl_urb3d,      tgl_urb3d,                        & !H urban
                           sh_urb2d,       lh_urb2d,        g_urb2d,   rn_urb2d,  ts_urb2d, & !H urban
                         psim_urb2d,     psih_urb2d,      u10_urb2d,  v10_urb2d,            & !O urban
                       GZ1OZ0_urb2d,     AKMS_URB2D,                                        & !O urban
                          th2_urb2d,       q2_urb2d,      ust_urb2d,                        & !O urban
                         declin_urb,      omg_urb2d,                                        & !I urban
                    num_roof_layers,num_wall_layers,num_road_layers,                        & !I urban
                                dzr,            dzb,            dzg,                        & !I urban
                         cmcr_urb2d,      tgr_urb2d,     tgrl_urb3d,  smr_urb3d,            & !H urban
                        drelr_urb2d,    drelb_urb2d,    drelg_urb2d,                        & !H urban
                      flxhumr_urb2d,  flxhumb_urb2d,  flxhumg_urb2d,                        & !H urban
                       tvg_urb2d,      xxxvg_urb2d,    tvgl_urb3d,     smg_urb3d, tt_urb2d, & !H urban
                      cmcg_urb2d, flxhumvg_urb2d, flxhumt_urb2d,                            & !H urban
                             julian,          julyr,                                        & !H urban
                          frc_urb2d,    utype_urb2d,                                        & !I urban
                                chs,           chs2,           cqs2,                        & !H
                      num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd,          & !I multi-layer urban
                       urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,          & !I multi-layer urban
                      urban_map_gbd,  urban_map_fbd, urban_map_zgrd,                        & !I multi-layer urban
                       num_urban_hi,                                                        & !I multi-layer urban
                          trb_urb4d,      tw1_urb4d,      tw2_urb4d,  tgb_urb4d,            & !H multi-layer urban
                         tlev_urb3d,     qlev_urb3d,                                        & !H multi-layer urban
                       tw1lev_urb3d,   tw2lev_urb3d,                                        & !H multi-layer urban
                        tglev_urb3d,    tflev_urb3d,                                        & !H multi-layer urban
                        sf_ac_urb3d,    lf_ac_urb3d,    cm_ac_urb3d,                        & !H multi-layer urban
                       sfvent_urb3d,   lfvent_urb3d,                                        & !H multi-layer urban
                       sfwin1_urb3d,   sfwin2_urb3d,                                        & !H multi-layer urban
                         sfw1_urb3d,     sfw2_urb3d,      sfr_urb3d,  sfg_urb3d,            & !H multi-layer urban
                        ep_pv_urb3d,     t_pv_urb3d,                                        & !RMS
                          trv_urb4d,       qr_urb4d,      qgr_urb3d,  tgr_urb3d,            & !RMS
                        drain_urb4d,  draingr_urb3d,     sfrv_urb3d, lfrv_urb3d,            & !RMS
                          dgr_urb3d,       dg_urb3d,      lfr_urb3d,  lfg_urb3d,            & !RMS
                           lp_urb2d,       hi_urb2d,       lb_urb2d,  hgt_urb2d,            & !H multi-layer urban
                           mh_urb2d,     stdh_urb2d,       lf_urb2d,                        & !SLUCM
                         lf_urb2d_s,       z0_urb2d,         vegfra,                        & !SLUCM
                             th_phy,            rho,          p_phy,        ust,            & !I multi-layer urban
                                gmt,         julday,          xlong,       xlat,            & !I multi-layer urban
                            a_u_bep,        a_v_bep,        a_t_bep,    a_q_bep,            & !O multi-layer urban
                            a_e_bep,        b_u_bep,        b_v_bep,                        & !O multi-layer urban
                            b_t_bep,        b_q_bep,        b_e_bep,    dlg_bep,            & !O multi-layer urban
                           dl_u_bep,         sf_bep,         vl_bep,                        & !O multi-layer urban
                                 T2,             Q2,           PSFC,        TH2,            &
                                U10,            V10,           PSIM,       PSIH,            &
                             GZ1OZ0, distributed_ahe_opt,     xtime,  ahe_urb2d             )

!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

    type(NoahmpIO_type),           intent(in   ) :: NoahmpIO
    INTEGER,                       INTENT(IN   ) :: sf_urban_physics   ! urban physics option
    INTEGER,                       INTENT(IN   ) :: NSOIL              ! number of soil layers
    INTEGER,                       INTENT(IN   ) :: ITIMESTEP          ! titestep number
    INTEGER,                       INTENT(IN   ) :: its,ite, kts,kte
    INTEGER,                       INTENT(IN   ) :: distributed_ahe_opt
    INTEGER, DIMENSION( its:ite ), INTENT(IN   ) :: IVGTYP             ! vegetation type
    REAL(kind=RKIND),              INTENT(IN   ) :: DT                 ! titestep [s]
    REAL(kind=RKIND), OPTIONAL,    INTENT(IN   ) :: xtime
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  COSZ_URB2D
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  XLAT_URB2D
    REAL(kind=RKIND),    DIMENSION( kts:kte, its:ite ), INTENT(IN   ) ::  T3D       ! 3D atmospheric temperature valid at mid-levels [K]
    REAL(kind=RKIND),    DIMENSION( kts:kte, its:ite ), INTENT(IN   ) ::  QV3D      ! 3D water vapor mixing ratio [kg/kg_dry]
    REAL(kind=RKIND),    DIMENSION( kts:kte, its:ite ), INTENT(IN   ) ::  U_PHY     ! 3D U wind component [m/s]
    REAL(kind=RKIND),    DIMENSION( kts:kte, its:ite ), INTENT(IN   ) ::  V_PHY     ! 3D V wind component [m/s]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  SWDOWN    ! solar down at surface [W m-2]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  SWDDIF    ! solar down at surface [W m-2]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  SWDDIR    ! solar down at surface [W m-2]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  GLW       ! longwave down at surface [W m-2]
    REAL(kind=RKIND),    DIMENSION( kts:kte, its:ite ), INTENT(IN   ) ::  P8W3D     ! 3D pressure, valid at interface [Pa]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  RAINBL    ! total input precipitation [mm]
    REAL(kind=RKIND),    DIMENSION( kts:kte, its:ite ), INTENT(IN   ) ::  DZ8W      ! thickness of atmo layers [m]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  ZNT       ! combined z0 sent to coupled model
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  TSK       ! surface radiative temperature [K]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  HFX       ! sensible heat flux [W m-2]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  QFX       ! latent heat flux [kg s-1 m-2]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  LH        ! latent heat flux [W m-2]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  GRDFLX    ! ground/snow heat flux [W m-2]
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  ALBEDO    ! total grid albedo []
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  EMISS     ! surface bulk emissivity
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  QSFC      ! bulk surface mixing ratio
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  Q2        ! 2-m mixing ratio (kg/kg)
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  T2        ! 2-m air temperature (K)
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(IN   ) ::  PSFC      ! surface pressure (pa)
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  TH2       ! 2-m potential temperature (K)
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  U10       ! 10-m wind (m/s)
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  V10       ! 10-m wind (m/s)
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  PSIM      ! similarity stability function for momentum
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  PSIH      ! similarity stability function for heat
    REAL(kind=RKIND),    DIMENSION( its:ite          ), INTENT(INOUT) ::  GZ1OZ0    ! log(z/z0) where z0 is roughness length

    ! input variables surface_driver --> lsm

     INTEGER,                       INTENT(IN   ) :: julian, julyr 
     INTEGER,                       INTENT(IN   ) :: num_roof_layers
     INTEGER,                       INTENT(IN   ) :: num_wall_layers
     INTEGER,                       INTENT(IN   ) :: num_road_layers
     INTEGER, DIMENSION( its:ite ), INTENT(IN   ) :: UTYPE_URB2D

     REAL(kind=RKIND),           DIMENSION( its:ite ),          INTENT(IN   ) :: FRC_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION(1:num_roof_layers),  INTENT(IN   ) :: DZR
     REAL(kind=RKIND), OPTIONAL, DIMENSION(1:num_wall_layers),  INTENT(IN   ) :: DZB
     REAL(kind=RKIND), OPTIONAL, DIMENSION(1:num_road_layers),  INTENT(IN   ) :: DZG
     REAL(kind=RKIND), OPTIONAL,                                INTENT(IN   ) :: DECLIN_URB
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),          INTENT(IN   ) :: OMG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ), INTENT(IN   ) :: TH_PHY
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ), INTENT(IN   ) :: P_PHY
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ), INTENT(IN   ) :: RHO
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),          INTENT(INOUT) :: UST
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),          INTENT(INOUT) :: CHS, CHS2, CQS2
     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:288, its:ite ),   INTENT(IN   ) :: ahe_urb2d

     ! local variables lsm --> urban

     INTEGER             :: UTYPE_URB    ! urban type
     LOGICAL             :: LSOLAR_URB
     REAL(kind=RKIND)    :: TA_URB       ! potential temp at 1st atmospheric level [K]
     REAL(kind=RKIND)    :: QA_URB       ! mixing ratio at 1st atmospheric level  [kg/kg]
     REAL(kind=RKIND)    :: UA_URB       ! wind speed at 1st atmospheric level    [m/s]
     REAL(kind=RKIND)    :: U1_URB       ! u at 1st atmospheric level             [m/s]
     REAL(kind=RKIND)    :: V1_URB       ! v at 1st atmospheric level             [m/s]
     REAL(kind=RKIND)    :: SSG_URB      ! downward total short wave radiation    [W/m/m]
     REAL(kind=RKIND)    :: LLG_URB      ! downward long wave radiation           [W/m/m]
     REAL(kind=RKIND)    :: RAIN_URB     ! precipitation                          [mm/h]
     REAL(kind=RKIND)    :: RHOO_URB     ! air density                            [kg/m^3]
     REAL(kind=RKIND)    :: ZA_URB       ! first atmospheric level                [m]
     REAL(kind=RKIND)    :: DELT_URB     ! tite step                              [s]
     REAL(kind=RKIND)    :: SSGD_URB     ! downward direct short wave radiation   [W/m/m]
     REAL(kind=RKIND)    :: SSGQ_URB     ! downward diffuse short wave radiation  [W/m/m]
     REAL(kind=RKIND)    :: XLAT_URB     ! latitude                               [deg]
     REAL(kind=RKIND)    :: COSZ_URB     ! cosz
     REAL(kind=RKIND)    :: OMG_URB      ! hour angle
     REAL(kind=RKIND)    :: ZNT_URB      ! roughness length                       [m]
     REAL(kind=RKIND)    :: TR_URB
     REAL(kind=RKIND)    :: TB_URB
     REAL(kind=RKIND)    :: TG_URB
     REAL(kind=RKIND)    :: TC_URB
     REAL(kind=RKIND)    :: QC_URB
     REAL(kind=RKIND)    :: UC_URB
     REAL(kind=RKIND)    :: XXXR_URB
     REAL(kind=RKIND)    :: XXXB_URB
     REAL(kind=RKIND)    :: XXXG_URB
     REAL(kind=RKIND)    :: XXXC_URB
     REAL(kind=RKIND), DIMENSION(1:num_roof_layers) :: TRL_URB  ! roof layer temp [K]
     REAL(kind=RKIND), DIMENSION(1:num_wall_layers) :: TBL_URB  ! wall layer temp [K]
     REAL(kind=RKIND), DIMENSION(1:num_road_layers) :: TGL_URB  ! road layer temp [K]

     !===hydrological variable for single layer UCM===

     INTEGER             :: jmonth, jday, ihour
     REAL(kind=RKIND)    :: DRELR_URB
     REAL(kind=RKIND)    :: DRELB_URB
     REAL(kind=RKIND)    :: DRELG_URB
     REAL(kind=RKIND)    :: FLXHUMR_URB
     REAL(kind=RKIND)    :: FLXHUMB_URB
     REAL(kind=RKIND)    :: FLXHUMG_URB
     REAL(kind=RKIND)    :: TVG_URB
     REAL(kind=RKIND)    :: TT_URB
     REAL(kind=RKIND)    :: XXXVG_URB
     REAL(kind=RKIND)    :: CMCG_URB
     REAL(kind=RKIND)    :: FLXHUMVG_URB
     REAL(kind=RKIND)    :: FLXHUMT_URB
     REAL(kind=RKIND)    :: CMCR_URB
     REAL(kind=RKIND)    :: TGR_URB

     REAL(kind=RKIND), DIMENSION(1:num_roof_layers) :: SMR_URB  ! green roof layer moisture
     REAL(kind=RKIND), DIMENSION(1:num_roof_layers) :: TGRL_URB ! green roof layer temp [K]
     REAL(kind=RKIND), DIMENSION(1:num_road_layers) :: TVGL_URB ! vegetated ground layer temp [K]
     REAL(kind=RKIND), DIMENSION(1:num_road_layers) :: SMG_URB  ! vegetated ground layer moisture

     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: DRELR_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: DRELB_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: DRELG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: FLXHUMR_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: FLXHUMB_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: FLXHUMG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TVG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TT_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: XXXVG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMCG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: FLXHUMVG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: FLXHUMT_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMCR_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TGR_URB2D

     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_roof_layers, its:ite ), INTENT(INOUT) :: TGRL_URB3D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_roof_layers, its:ite ), INTENT(INOUT) :: SMR_URB3D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_road_layers, its:ite ), INTENT(INOUT) :: TVGL_URB3D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_road_layers, its:ite ), INTENT(INOUT) :: SMG_URB3D

     ! state variable surface_driver <--> lsm <--> urban

     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TR_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TB_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TC_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: QC_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: UC_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: XXXR_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: XXXB_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: XXXG_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: XXXC_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: SH_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: LH_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: G_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: RN_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: TS_URB2D

     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_roof_layers, its:ite ), INTENT(INOUT) :: TRL_URB3D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_wall_layers, its:ite ), INTENT(INOUT) :: TBL_URB3D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_road_layers, its:ite ), INTENT(INOUT) :: TGL_URB3D

     ! output variable lsm --> surface_driver

     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: PSIM_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: PSIH_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: GZ1OZ0_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: U10_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: V10_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: TH2_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: Q2_URB2D
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(OUT) :: AKMS_URB2D
     REAL(kind=RKIND),           DIMENSION( its:ite ), INTENT(OUT) :: UST_URB2D


     ! output variables urban --> lsm

     REAL(kind=RKIND) :: TS_URB           ! surface radiative temperature    [K]
     REAL(kind=RKIND) :: QS_URB           ! surface humidity                 [-]
     REAL(kind=RKIND) :: SH_URB           ! sensible heat flux               [W/m/m]
     REAL(kind=RKIND) :: LH_URB           ! latent heat flux                 [W/m/m]
     REAL(kind=RKIND) :: LH_KINEMATIC_URB ! latent heat flux, kinetic  [kg/m/m/s]
     REAL(kind=RKIND) :: SW_URB           ! upward short wave radiation flux [W/m/m]
     REAL(kind=RKIND) :: ALB_URB          ! tite-varying albedo            [fraction]
     REAL(kind=RKIND) :: LW_URB           ! upward long wave radiation flux  [W/m/m]
     REAL(kind=RKIND) :: G_URB            ! heat flux into the ground        [W/m/m]
     REAL(kind=RKIND) :: RN_URB           ! net radiation                    [W/m/m]
     REAL(kind=RKIND) :: PSIM_URB         ! shear f for momentum             [-]
     REAL(kind=RKIND) :: PSIH_URB         ! shear f for heat                 [-]
     REAL(kind=RKIND) :: GZ1OZ0_URB       ! shear f for heat                 [-]
     REAL(kind=RKIND) :: U10_URB          ! wind u component at 10 m         [m/s]
     REAL(kind=RKIND) :: V10_URB          ! wind v component at 10 m         [m/s]
     REAL(kind=RKIND) :: TH2_URB          ! potential temperature at 2 m     [K]
     REAL(kind=RKIND) :: Q2_URB           ! humidity at 2 m                  [-]
     REAL(kind=RKIND) :: CHS_URB
     REAL(kind=RKIND) :: CHS2_URB
     REAL(kind=RKIND) :: UST_URB
     REAL(kind=RKIND) :: EMISS_URB

     ! NUDAPT Parameters urban --> lam

     REAL(kind=RKIND) :: mh_urb
     REAL(kind=RKIND) :: stdh_urb
     REAL(kind=RKIND) :: lp_urb
     REAL(kind=RKIND) :: hgt_urb
     REAL(kind=RKIND), DIMENSION(4) :: lf_urb

     ! Distributed aerodynamics parameters

     REAL(kind=RKIND) :: lf_urb_s
     REAL(kind=RKIND) :: z0_urb
     REAL(kind=RKIND) :: vegfrac

     ! Local variables

     INTEGER :: I,J,K
     REAL(kind=RKIND) :: Q1

     ! Noah UA changes

     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMR_SFCDIF
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CHR_SFCDIF
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMGR_SFCDIF
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CHGR_SFCDIF
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CMC_SFCDIF
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(INOUT) :: CHC_SFCDIF

     ! Variables for multi-layer UCM

     REAL(kind=RKIND), OPTIONAL,                       INTENT(IN   ) :: gmt
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: JULDAY
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ), INTENT(IN   ) :: XLAT, XLONG
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: num_urban_ndm
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_zrd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_zwd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_gd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_zd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_zdf
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_bd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_wd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_gbd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_fbd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: urban_map_zgrd
     INTEGER, OPTIONAL,                                INTENT(IN   ) :: NUM_URBAN_HI

     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_urban_hi, its:ite ), INTENT(IN   ) :: hi_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: lp_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: lb_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: hgt_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: mh_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: stdh_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( 4,its:ite ),               INTENT(IN   ) :: lf_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: lf_urb2d_s
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: z0_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                 INTENT(IN   ) :: vegfra

     ! temporary solution: making BEP/BEM vars optional (cenlinhe,10/18/2026)
     REAL(kind=RKIND), OPTIONAL, DIMENSION(1:,its:), INTENT(IN   ) :: hi_urb2d
     REAL(kind=RKIND), OPTIONAL, DIMENSION(1:,its:), INTENT(INOUT) :: &
                                            trb_urb4d, tw1_urb4d, tw2_urb4d, tgb_urb4d, &
                                            tlev_urb3d, qlev_urb3d, tw1lev_urb3d, tw2lev_urb3d, &
                                            tglev_urb3d, tflev_urb3d, &
                                            sfwin1_urb3d, sfwin2_urb3d, sfw1_urb3d, sfw2_urb3d, &
                                            sfr_urb3d, sfg_urb3d, t_pv_urb3d, &
                                            trv_urb4d, qr_urb4d, drain_urb4d, &
                                            sfrv_urb3d, lfrv_urb3d, dgr_urb3d, dg_urb3d, &
                                            lfr_urb3d, lfg_urb3d

     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zrd, its:ite ), INTENT(INOUT) :: trb_urb4d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zwd, its:ite ), INTENT(INOUT) :: tw1_urb4d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zwd, its:ite ), INTENT(INOUT) :: tw2_urb4d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_gd , its:ite ), INTENT(INOUT) :: tgb_urb4d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_bd , its:ite ), INTENT(INOUT) :: tlev_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_bd , its:ite ), INTENT(INOUT) :: qlev_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_wd , its:ite ), INTENT(INOUT) :: tw1lev_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_wd , its:ite ), INTENT(INOUT) :: tw2lev_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_gbd, its:ite ), INTENT(INOUT) :: tglev_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_fbd, its:ite ), INTENT(INOUT) :: tflev_urb3d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: lf_ac_urb3d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: sf_ac_urb3d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: cm_ac_urb3d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: sfvent_urb3d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: lfvent_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_wd , its:ite ), INTENT(INOUT) :: sfwin1_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_wd , its:ite ), INTENT(INOUT) :: sfwin2_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zd , its:ite ), INTENT(INOUT) :: sfw1_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zd , its:ite ), INTENT(INOUT) :: sfw2_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zdf, its:ite ), INTENT(INOUT) :: sfr_urb3d
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_urban_ndm, its:ite ), INTENT(INOUT) :: sfg_urb3d
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: ep_pv_urb3d   !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zdf , its:ite ),INTENT(INOUT) :: t_pv_urb3d    !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zgrd, its:ite ),INTENT(INOUT) :: trv_urb4d     !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zgrd, its:ite ),INTENT(INOUT) :: qr_urb4d      !GRZ
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: qgr_urb3d     !GRZ
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: tgr_urb3d     !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zdf, its:ite ), INTENT(INOUT) :: drain_urb4d   !GRZ
     REAL(kind=RKIND), OPTIONAL, DIMENSION( its:ite ),                  INTENT(INOUT) :: draingr_urb3d !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zdf, its:ite ), INTENT(INOUT) :: sfrv_urb3d    !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zdf, its:ite ), INTENT(INOUT) :: lfrv_urb3d    !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zdf, its:ite ), INTENT(INOUT) :: dgr_urb3d     !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_urban_ndm, its:ite ), INTENT(INOUT) :: dg_urb3d      !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:urban_map_zdf, its:ite ), INTENT(INOUT) :: lfr_urb3d     !GRZ
     !REAL(kind=RKIND), OPTIONAL, DIMENSION( 1:num_urban_ndm, its:ite ), INTENT(INOUT) :: lfg_urb3d     !GRZ
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: a_u_bep       !Implicit momemtum component X-direction
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: a_v_bep       !Implicit momemtum component Y-direction
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: a_t_bep       !Implicit component pot. temperature
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: a_q_bep       !Implicit momemtum component X-direction
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: a_e_bep       !Implicit component TKE
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: b_u_bep       !Explicit momentum component X-direction
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: b_v_bep       !Explicit momentum component Y-direction
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: b_t_bep       !Explicit component pot. temperature
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: b_q_bep       !Implicit momemtum component Y-direction
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: b_e_bep       !Explicit component TKE
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: vl_bep        !Fraction air volume in grid cell
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: dlg_bep       !Height above ground
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: sf_bep        !Fraction air at the face of grid cell
     REAL(kind=RKIND), OPTIONAL, DIMENSION( kts:kte, its:ite ),         INTENT(INOUT) :: dl_u_bep      !Length scale

     ! Local variables for multi-layer UCM

     REAL(kind=RKIND), DIMENSION( its:ite ) :: HFX_RURAL,GRDFLX_RURAL,UST_RURAL ! RN_RURAL
     REAL(kind=RKIND), DIMENSION( its:ite ) :: QFX_RURAL,LH_RURAL,QSFC_RURAL    ! UMOM_RURAL,VMOM_RURAL
     REAL(kind=RKIND), DIMENSION( its:ite ) :: ALB_RURAL,EMISS_RURAL,TSK_RURAL,Q2_RURAL,T2_RURAL
     REAL(kind=RKIND), DIMENSION( its:ite ) :: HFX_URB,UMOM_URB,VMOM_URB
     REAL(kind=RKIND), DIMENSION( its:ite ) :: QFX_URB
     REAL(kind=RKIND), DIMENSION( its:ite ) :: RL_UP_URB
     REAL(kind=RKIND), DIMENSION( its:ite ) :: RS_ABS_URB
     REAL(kind=RKIND), DIMENSION( its:ite ) :: GRDFLX_URB

     REAL(kind=RKIND) :: SIGMA_SB,RL_UP_RURAL,RL_UP_TOT,RS_ABS_TOT,UMOM,VMOM
     REAL(kind=RKIND) :: r1,r2,r3
     REAL(kind=RKIND) :: CMR_URB, CHR_URB, CMC_URB, CHC_URB, CMGR_URB, CHGR_URB
     REAL(kind=RKIND) :: frc_urb,lb_urb
     REAL(kind=RKIND) :: check

    character(len=80) :: message

    DO I=ITS,ITE
      HFX_RURAL(I)    = HFX(I)
      QFX_RURAL(I)    = QFX(I)
      GRDFLX_RURAL(I) = GRDFLX(I)
      EMISS_RURAL(I)  = EMISS(I)
      TSK_RURAL(I)    = TSK(I)
      ALB_RURAL(I)    = ALBEDO(I)
      LH_RURAL(I)     = LH(I)
      QSFC_RURAL(I)   = QSFC(I)
      UST_RURAL(I)    = UST(I)
      Q2_RURAL(I)     = Q2(I)
      T2_RURAL(I)     = T2(I)
    END DO

IF (SF_URBAN_PHYSICS == 1 ) THEN         ! Beginning of UCM CALL if block

!--------------------------------------
! URBAN CANOPY MODEL START
!--------------------------------------

ILOOP : DO I = its, ite


  IF ( IVGTYP(I) == NoahmpIO%ISURBAN_TABLE .or. IVGTYP(I) > NoahmpIO%URBTYPE_beg ) THEN

  IF ((FRC_URB2D(I) > 0.) .and. (FRC_URB2D(I) <= 1.)) THEN

    UTYPE_URB = UTYPE_URB2D(I) !urban type (low, high or industrial)

    TA_URB    = T3D(1,I)                                ! [K]            
    QA_URB    = QV3D(1,I)/(1.0+QV3D(1,I))               ! [kg/kg]       
    UA_URB    = SQRT(U_PHY(1,I)**2.+V_PHY(1,I)**2.)
    U1_URB    = U_PHY(1,I)
    V1_URB    = V_PHY(1,I)
    IF(UA_URB < 1.0) UA_URB=1.                           ! [m/s]
    SSG_URB   = SWDOWN(I)                               ! [W/m/m]      
    SSGD_URB  = 0.8*SWDOWN(I)                           ! [W/m/m]     
    SSGQ_URB  = SSG_URB-SSGD_URB                        ! [W/m/m]
    LLG_URB   = GLW(I)                                  ! [W/m/m]
    RAIN_URB  = RAINBL(I) / DT * 3600.0                 ! [mm/hr]       
    RHOO_URB  = (P8W3D(KTS+1,I)+P8W3D(KTS,I))*0.5 / (287.04 * TA_URB * (1.0+ 0.61 * QA_URB)) ![kg/m/m/m] 
    ZA_URB    = 0.5*DZ8W(1,I)                           ! [m]         
    DELT_URB  = DT                                      ! [sec]
    XLAT_URB  = XLAT_URB2D(I)                           ! [deg]
    COSZ_URB  = COSZ_URB2D(I) 
    OMG_URB   = OMG_URB2D(I)
    ZNT_URB   = ZNT(I)

    LSOLAR_URB = .FALSE.

    TR_URB = TR_URB2D(I)
    TB_URB = TB_URB2D(I)
    TG_URB = TG_URB2D(I)
    TC_URB = TC_URB2D(I)
    QC_URB = QC_URB2D(I)
    UC_URB = UC_URB2D(I)

    TGR_URB     = TGR_URB2D(I)
    CMCR_URB    = CMCR_URB2D(I)
    FLXHUMR_URB = FLXHUMR_URB2D(I)
    FLXHUMB_URB = FLXHUMB_URB2D(I)
    FLXHUMG_URB = FLXHUMG_URB2D(I)
    TVG_URB     = TVG_URB2D(I)
    TT_URB      = TT_URB2D(I)
    XXXVG_URB   = XXXVG_URB2D(I)
    CMCG_URB    = CMCG_URB2D(I)
    FLXHUMVG_URB= FLXHUMVG_URB2D(I)
    FLXHUMT_URB = FLXHUMT_URB2D(I)
    DRELR_URB   = DRELR_URB2D(I)
    DRELB_URB   = DRELB_URB2D(I)
    DRELG_URB   = DRELG_URB2D(I)

    DO K = 1,num_roof_layers
      TRL_URB(K) = TRL_URB3D(K,I)
      SMR_URB(K) = SMR_URB3D(K,I)
      TGRL_URB(K)= TGRL_URB3D(K,I)
    END DO

    DO K = 1,num_wall_layers
      TBL_URB(K) = TBL_URB3D(K,I)
    END DO

    DO K = 1,num_road_layers
      TGL_URB(K) = TGL_URB3D(K,I)
      TVGL_URB(K)= TVGL_URB3D(K,I)
      SMG_URB(K) = SMG_URB3D(K,I)
    END DO

    XXXR_URB = XXXR_URB2D(I)
    XXXB_URB = XXXB_URB2D(I)
    XXXG_URB = XXXG_URB2D(I)
    XXXC_URB = XXXC_URB2D(I)

    ! Limits to avoid dividing by small number
    IF (CHS(I) < 1.0E-02) THEN
      CHS(I) = 1.0E-02
    ENDIF
    IF (CHS2(I) < 1.0E-02) THEN
      CHS2(I) = 1.0E-02
    ENDIF
    IF (CQS2(I) < 1.0E-02) THEN
      CQS2(I) = 1.0E-02
    ENDIF

    CHS_URB  = CHS(I)
    CHS2(I)  = CQS2(I)      
    CHS2_URB = CHS2(I)
    IF (PRESENT(CMR_SFCDIF)) THEN
      CMR_URB  = CMR_SFCDIF(I)
      CHR_URB  = CHR_SFCDIF(I)
      CMGR_URB = CMGR_SFCDIF(I)
      CHGR_URB = CHGR_SFCDIF(I)
      CMC_URB  = CMC_SFCDIF(I)
      CHC_URB  = CHC_SFCDIF(I)
    ENDIF

    ! NUDAPT for SLUCM

    MH_URB   = MH_URB2D(I)
    STDH_URB = STDH_URB2D(I)
    LP_URB   = LP_URB2D(I)
    HGT_URB  = HGT_URB2D(I)
    LF_URB   = 0.0
    DO K = 1,4
      LF_URB(K) = LF_URB2D(K,I)
    ENDDO
    FRC_URB  = FRC_URB2D(I)
    LB_URB   = LB_URB2D(I)
    CHECK    = 0
    !IF (I.EQ.73)THEN
      CHECK = 0
    !END IF
    ! Distributed aerodynamics
    lf_urb_s = lf_urb2d_s(I)
    z0_urb   = z0_urb2d(I)
    vegfrac  = vegfra(I) / 100

    ! Call urban

    CALL cal_mon_day(julian,julyr,jmonth,jday) ! julday in MPAS used for julian
    CALL urban(LSOLAR_URB,                                                             & ! I
          num_roof_layers, num_wall_layers, num_road_layers,                           & ! C
                DZR,        DZB,        DZG,                                           & ! C
          UTYPE_URB,     TA_URB,     QA_URB,     UA_URB,   U1_URB,  V1_URB, SSG_URB,   & ! I
           SSGD_URB,   SSGQ_URB,    LLG_URB,   RAIN_URB, RHOO_URB,                     & ! I
             ZA_URB, DECLIN_URB,   COSZ_URB,    OMG_URB,                               & ! I
           XLAT_URB,   DELT_URB,    ZNT_URB,                                           & ! I
            CHS_URB,   CHS2_URB,                                                       & ! I
             TR_URB,     TB_URB,     TG_URB,     TC_URB,   QC_URB,   UC_URB,           & ! H
            TRL_URB,    TBL_URB,    TGL_URB,                                           & ! H
           XXXR_URB,   XXXB_URB,   XXXG_URB,   XXXC_URB,                               & ! H
             TS_URB,     QS_URB,     SH_URB,     LH_URB, LH_KINEMATIC_URB,             & ! O
             SW_URB,    ALB_URB,     LW_URB,      G_URB,   RN_URB, PSIM_URB, PSIH_URB, & ! O
         GZ1OZ0_URB,                                                                   & ! O
            CMR_URB,    CHR_URB,    CMC_URB,    CHC_URB,                               &
            U10_URB,    V10_URB,    TH2_URB,     Q2_URB,                               & ! O
            UST_URB,     mh_urb,   stdh_urb,     lf_urb,   lp_urb,                     & ! O
            hgt_urb,    frc_urb,     lb_urb,      check, CMCR_URB,TGR_URB,             & ! H
           TGRL_URB,    SMR_URB,   CMGR_URB,   CHGR_URB,   jmonth,                     & ! H
          DRELR_URB,  DRELB_URB,                                                       & ! H
          DRELG_URB,FLXHUMR_URB,FLXHUMB_URB,FLXHUMG_URB,                               &
          TVG_URB,TT_URB,XXXVG_URB,TVGL_URB,SMG_URB,CMCG_URB,FLXHUMVG_URB,FLXHUMT_URB, &
          lf_urb_s, z0_urb, vegfrac, EMISS_URB)

    TS_URB2D(I) = TS_URB

    ALBEDO(I)   = FRC_URB2D(I) * ALB_URB + (1-FRC_URB2D(I)) * ALB_RURAL(I)           ![-]      
    HFX(I)      = FRC_URB2D(I) * SH_URB  + (1-FRC_URB2D(I)) * HFX_RURAL(I)           ![W/m/m] 
    QFX(I)      = FRC_URB2D(I) * LH_KINEMATIC_URB + (1-FRC_URB2D(I))* QFX_RURAL(I)   ![kg/m/m/s] 
    LH(I)       = FRC_URB2D(I) * LH_URB  + (1-FRC_URB2D(I)) * LH_RURAL(I)            ![W/m/m]   
    GRDFLX(I)   = FRC_URB2D(I) * (G_URB*(-1.0)) + (1-FRC_URB2D(I)) * GRDFLX_RURAL(I) ![W/m/m] positive: downward direction
    TSK(I)      = FRC_URB2D(I) * TS_URB  + (1-FRC_URB2D(I)) * TSK_RURAL(I)           ![K]
    EMISS(I)    = FRC_URB2D(I) * EMISS_URB + (1-FRC_URB2D(I)) * EMISS_RURAL(I)     
    !Q1         = QSFC(I)/(1.0+QSFC(I))                                         
    !Q1         = FRC_URB2D(I) * QS_URB  + (1-FRC_URB2D(I)) * Q1               ![-]

   ! Convert QSFC back to mixing ratio

   !QSFC(I) = Q1/(1.0-Q1)
    QSFC(I)     = FRC_URB2D(I) * QS_URB  + (1-FRC_URB2D(I)) * QSFC_RURAL(I)    !!   QSFC(I)=QSFC1D
    UST(I)      = FRC_URB2D(I) * UST_URB + (1-FRC_URB2D(I)) * UST_RURAL(I)     ![m/s]

   ! Renew Urban State Variables

    TR_URB2D(I) = TR_URB
    TB_URB2D(I) = TB_URB
    TG_URB2D(I) = TG_URB
    TC_URB2D(I) = TC_URB
    QC_URB2D(I) = QC_URB
    UC_URB2D(I) = UC_URB

    TGR_URB2D(I)     = TGR_URB
    CMCR_URB2D(I)    = CMCR_URB
    FLXHUMR_URB2D(I) = FLXHUMR_URB
    FLXHUMB_URB2D(I) = FLXHUMB_URB
    FLXHUMG_URB2D(I) = FLXHUMG_URB
    TVG_URB2D(I)     = TVG_URB
    TT_URB2D(I)      = TT_URB
    XXXVG_URB2D(I)   = XXXVG_URB
    CMCG_URB2D(I)    = CMCG_URB
    FLXHUMVG_URB2D(I)= FLXHUMVG_URB
    FLXHUMT_URB2D(I) = FLXHUMT_URB
    DRELR_URB2D(I)   = DRELR_URB
    DRELB_URB2D(I)   = DRELB_URB
    DRELG_URB2D(I)   = DRELG_URB

    DO K = 1,num_roof_layers
      TRL_URB3D(K,I) = TRL_URB(K)
      SMR_URB3D(K,I) = SMR_URB(K)
      TGRL_URB3D(K,I)= TGRL_URB(K)
    END DO
    DO K = 1,num_wall_layers
      TBL_URB3D(K,I) = TBL_URB(K)
    END DO
    DO K = 1,num_road_layers
      TGL_URB3D(K,I) = TGL_URB(K)
      TVGL_URB3D(K,I)= TVGL_URB(K)
      SMG_URB3D(K,I) = SMG_URB(K)
    END DO

    XXXR_URB2D(I)    = XXXR_URB
    XXXB_URB2D(I)    = XXXB_URB
    XXXG_URB2D(I)    = XXXG_URB
    XXXC_URB2D(I)    = XXXC_URB
    SH_URB2D(I)      = SH_URB
    LH_URB2D(I)      = LH_URB
    G_URB2D(I)       = G_URB         
    RN_URB2D(I)      = RN_URB
    PSIM_URB2D(I)    = PSIM_URB
    PSIH_URB2D(I)    = PSIH_URB
    GZ1OZ0_URB2D(I)  = GZ1OZ0_URB
    U10_URB2D(I)     = U10_URB
    V10_URB2D(I)     = V10_URB
    TH2_URB2D(I)     = TH2_URB
    Q2_URB2D(I)      = Q2_URB
    UST_URB2D(I)     = UST_URB
    AKMS_URB2D(I)    = KARMAN * UST_URB2D(I)/(GZ1OZ0_URB2D(I)-PSIM_URB2D(I))
    IF (PRESENT(CMR_SFCDIF)) THEN
      CMR_SFCDIF(I)  = CMR_URB
      CHR_SFCDIF(I)  = CHR_URB
      CMGR_SFCDIF(I) = CMGR_URB
      CHGR_SFCDIF(I) = CHGR_URB
      CMC_SFCDIF(I)  = CMC_URB
      CHC_SFCDIF(I)  = CHC_URB
    ENDIF

    ! cenlin(10/17/2026): update diagnostic host model variables (moved from WRF surface driver to here)
    Q2(I)     = FRC_URB2D(I) * Q2_URB + (1.-FRC_URB2D(I)) * Q2_RURAL(I)
    T2(I)     = FRC_URB2D(I) * (TH2_URB/((1.E5/PSFC(I))**RCP)) + (1.-FRC_URB2D(I)) * T2_RURAL(I)
    TH2(I)    = T2(I) * (1.E5/PSFC(I))**RCP
    U10(I)    = U10_URB
    V10(I)    = V10_URB  
    PSIM(I)   = PSIM_URB
    PSIH(I)   = PSIH_URB
    GZ1OZ0(I) = GZ1OZ0_URB 
    !AKHS(I)   = CHS(I)        ! only needed for MYJ scheme, not available in MPAS
    !AKMS(I)   = AKMS_URB2D(I) ! only needed for MYJ scheme, not available in MPAS

    ! cenlin(10/17/2026): move anthropogenic heat treatment from  WRF surface driver to here
    IF (distributed_ahe_opt == 2) THEN
       ihour  = 1 + (jmonth-1) * 24 + MOD(INT(gmt+xtime/60.0),24) ! convert ihour from 0-based to 1-based
       HFX(I) = HFX(I) + ahe_urb2d(ihour,I)
    ENDIF

  ELSE ! FRC_URB2D <=0 or >1

    SH_URB2D(I)      = 0.
    LH_URB2D(I)      = 0.
    G_URB2D(I)       = 0.
    RN_URB2D(I)      = 0.
    TS_URB2D(I)      = 0.

  ENDIF    ! FRC_URB2D

  ENDIF      ! urban land used type block

ENDDO ILOOP  ! of I loop

ENDIF        ! sf_urban_physics = 1 block

!--------------------------------------
! URBAN CANOPY MODEL END
!--------------------------------------

!--------------------------------------
! URBAN BEP and BEM MODEL BEGIN
!--------------------------------------

!IF (SF_URBAN_PHYSICS == 2) THEN
!
!DO I=ITS,ITE
!
!  EMISS_URB(I)       = 0.
!  RL_UP_URB(I)       = 0.
!  RS_ABS_URB(I)      = 0.
!  GRDFLX_URB(I)      = 0.
!  B_Q_BEP(KTS:KTE,I) = 0.
!
!END DO
!
!  CALL BEP(frc_urb2d,  utype_urb2d, ititestep,       dz8w,         &
!                  dt,        u_phy,     v_phy,                     &
!              th_phy,          rho,     p_phy,     swdown,    glw, &
!                 gmt,       julday,     xlong,       xlat,         &
!          declin_urb,   cosz_urb2d, omg_urb2d,                     &
!       num_urban_ndm, urban_map_zrd, urban_map_zwd, urban_map_gd,  &
!        urban_map_zd, urban_map_zdf,  urban_map_bd, urban_map_wd,  &
!       urban_map_gbd, urban_map_fbd,  num_urban_hi,                &
!           trb_urb4d,    tw1_urb4d, tw2_urb4d,  tgb_urb4d,         &
!          sfw1_urb3d,   sfw2_urb3d, sfr_urb3d,  sfg_urb3d,         &
!            lp_urb2d,     hi_urb2d,  lb_urb2d,  hgt_urb2d,         &
!             a_u_bep,      a_v_bep,   a_t_bep,                     &
!             a_e_bep,      b_u_bep,   b_v_bep,                     &
!             b_t_bep,      b_e_bep,   b_q_bep,    dlg_bep,         &
!            dl_u_bep,       sf_bep,    vl_bep,                     &
!           rl_up_urb,   rs_abs_urb, emiss_urb, grdflx_urb,         &
!         ids,ide, jds,jde, kds,kde,                                &
!         ims,ime, jms,jme, kms,kme,                                &
!         its,ite, jts,jte, kts,kte )
!
!ENDIF ! SF_URBAN_PHYSICS == 2
!
!IF (SF_URBAN_PHYSICS == 3) THEN
!
!DO I=ITS,ITE
!
!  EMISS_URB(I)       = 0.
!  RL_UP_URB(I)       = 0.
!  RS_ABS_URB(I)      = 0.
!  GRDFLX_URB(I)      = 0.
!  B_Q_BEP(KTS:KTE,I) = 0.
!
!END DO
!
!  CALL BEP_BEM( frc_urb2d,  utype_urb2d,    ititestep,         dz8w,       &
!                       dt,        u_phy,        v_phy,                     &
!                   th_phy,          rho,        p_phy,       swdown,  glw, &
!                      gmt,       julday,        xlong,         xlat,       &
!               declin_urb,   cosz_urb2d,    omg_urb2d,                     &
!            num_urban_ndm, urban_map_zrd, urban_map_zwd, urban_map_gd,     &
!             urban_map_zd, urban_map_zdf,  urban_map_bd, urban_map_wd,     &
!            urban_map_gbd, urban_map_fbd,  urban_map_zgrd,num_urban_hi,    &
!                trb_urb4d,    tw1_urb4d,    tw2_urb4d,    tgb_urb4d,       &
!               tlev_urb3d,   qlev_urb3d, tw1lev_urb3d, tw2lev_urb3d,       &
!              tglev_urb3d,  tflev_urb3d,  sf_ac_urb3d,  lf_ac_urb3d,       &
!              cm_ac_urb3d, sfvent_urb3d, lfvent_urb3d,                     &
!             sfwin1_urb3d, sfwin2_urb3d,                                   &
!               sfw1_urb3d,   sfw2_urb3d,    sfr_urb3d,    sfg_urb3d,       &
!              ep_pv_urb3d,   t_pv_urb3d,                                   & !RMS
!                trv_urb4d,     qr_urb4d,    qgr_urb3d,   tgr_urb3d,        & !RMS
!              drain_urb4d,draingr_urb3d,   sfrv_urb3d,  lfrv_urb3d,        & !RMS
!                dgr_urb3d,     dg_urb3d,    lfr_urb3d,   lfg_urb3d,        & !RMS
!                   rainbl,       swddir,       swddif,                     &
!                 lp_urb2d,     hi_urb2d,     lb_urb2d,    hgt_urb2d,       &
!                  a_u_bep,      a_v_bep,      a_t_bep,                     &
!                  a_e_bep,      b_u_bep,      b_v_bep,                     &
!                  b_t_bep,      b_e_bep,      b_q_bep,      dlg_bep,       &
!                 dl_u_bep,       sf_bep,       vl_bep,                     &
!                rl_up_urb,   rs_abs_urb,    emiss_urb,   grdflx_urb, qv3d, &
!             ids,ide, jds,jde, kds,kde,                                    &
!             ims,ime, jms,jme, kms,kme,                                    &
!             its,ite, jts,jte, kts,kte )
!
!ENDIF ! SF_URBAN_PHYSICS == 3
!
!IF((SF_URBAN_PHYSICS == 2).OR.(SF_URBAN_PHYSICS == 3))THEN 
!
!  sigma_sb=5.67e-08
!  do i = its, ite
!    UMOM_URB(I)     = 0.
!    VMOM_URB(I)     = 0.
!    HFX_URB(I)      = 0.
!    QFX_URB(I)      = 0.
!
!    do k=kts,kte
!      a_u_bep(k,i) = a_u_bep(k,i)*frc_urb2d(i)
!      a_v_bep(k,i) = a_v_bep(k,i)*frc_urb2d(i)
!      a_t_bep(k,i) = a_t_bep(k,i)*frc_urb2d(i)
!      a_q_bep(k,i) = 0.
!      a_e_bep(k,i) = 0.
!      b_u_bep(k,i) = b_u_bep(k,i)*frc_urb2d(i)
!      b_v_bep(k,i) = b_v_bep(k,i)*frc_urb2d(i)
!      b_t_bep(k,i) = b_t_bep(k,i)*frc_urb2d(i)
!      b_q_bep(k,i) = b_q_bep(k,i)*frc_urb2d(i)
!      b_e_bep(k,i) = b_e_bep(k,i)*frc_urb2d(i)
!      HFX_URB(I)   = HFX_URB(I) + B_T_BEP(K,I)*RHO(K,I)*CP*DZ8W(K,I)*VL_BEP(K,I)
!      QFX_URB(I)   = QFX_URB(I) + B_Q_BEP(K,I)*DZ8W(K,I)*VL_BEP(K,I)
!      UMOM_URB(I)  = UMOM_URB(I)+ (A_U_BEP(K,I)*U_PHY(K,I)+B_U_BEP(K,I))*DZ8W(K,I)*VL_BEP(K,I)
!      VMOM_URB(I)  = VMOM_URB(I)+ (A_V_BEP(K,I)*V_PHY(K,I)+B_V_BEP(K,I))*DZ8W(K,I)*VL_BEP(K,I)
!      vl_bep(k,i)  = (1.-frc_urb2d(i)) + vl_bep(k,i)*frc_urb2d(i)
!      sf_bep(k,i)  = (1.-frc_urb2d(i)) + sf_bep(k,i)*frc_urb2d(i)
!    end do
!
!    a_u_bep(1,i)   = (1.-frc_urb2d(i))*(-ust(I)*ust(I))/dz8w(1,i)/   &
!                          ((u_phy(1,i)**2+v_phy(1,i)**2.)**.5)+a_u_bep(1,i)
!
!    a_v_bep(1,i)   = (1.-frc_urb2d(i))*(-ust(I)*ust(I))/dz8w(1,i)/   &
!                          ((u_phy(1,i)**2+v_phy(1,i)**2.)**.5)+a_v_bep(1,i)
!
!    b_t_bep(1,i)   = (1.-frc_urb2d(i))*hfx_rural(i)/dz8w(1,i)/rho(1,i)/CP+ & 
!                           b_t_bep(1,i)
!
!    b_q_bep(1,i)   = (1.-frc_urb2d(i))*qfx_rural(i)/dz8w(1,i)/rho(1,i)+b_q_bep(1,i)
!
!    umom             = (1.-frc_urb2d(i))*ust(i)*ust(i)*u_phy(1,i)/               &
!                         ((u_phy(1,i)**2+v_phy(1,i)**2.)**.5)+umom_urb(i)
!
!    vmom             = (1.-frc_urb2d(i))*ust(i)*ust(i)*v_phy(1,i)/               &
!                         ((u_phy(1,i)**2+v_phy(1,i)**2.)**.5)+vmom_urb(i)
!    sf_bep(1,i)    = 1.
!
!  ! using the emissivity and the total longwave upward radiation estimate the averaged skin temperature
!
!  IF (FRC_URB2D(I).GT.0.) THEN
!    rl_up_rural   = -emiss_rural(i)*sigma_sb*(tsk_rural(i)**4.)-(1.-emiss_rural(i))*glw(i)
!    rl_up_tot     = (1.-frc_urb2d(i))*rl_up_rural     + frc_urb2d(i)*rl_up_urb(i)
!    emiss(i)    = (1.-frc_urb2d(i))*emiss_rural(i)+ frc_urb2d(i)*emiss_urb(i)
!    ts_urb2d(i) = (max(0.,(-rl_up_urb(i)-(1.-emiss_urb(i))*glw(i))/emiss_urb(i)/sigma_sb))**0.25
!    tsk(i)      = (max(0., (-1.*rl_up_tot-(1.-emiss(i))*glw(i) )/emiss(i)/sigma_sb))**.25
!    rs_abs_tot    = (1.-frc_urb2d(i))*swdown(i)*(1.-albedo(i))+frc_urb2d(i)*rs_abs_urb(i)
!
!    if(swdown(i) > 0.)then
!      albedo(i) = 1.-rs_abs_tot/swdown(i)
!    else
!      albedo(i) = alb_rural(i)
!    endif
!
!    ! rename *_urb to sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d
!
!    grdflx(i)   = (1.-frc_urb2d(i))*grdflx_rural(i)+ frc_urb2d(i)*grdflx_urb(i)*(-1.0) ! positive: downward direction
!    qfx(i)      = (1.-frc_urb2d(i))*qfx_rural(i)   + qfx_urb(i)
!    lh(i)       = qfx(i)*xlv
!    hfx(i)      = hfx_urb(i) + (1-frc_urb2d(i))*hfx_rural(i)      ![W/m/m]
!    sh_urb2d(i) = hfx_urb(i)/frc_urb2d(i)
!    lh_urb2d(i) = qfx_urb(i)*xlv/frc_urb2d(i)
!    g_urb2d(i)  = grdflx_urb(i)
!    rn_urb2d(i) = rs_abs_urb(i)+emiss_urb(i)*glw(i)-rl_up_urb(i)
!    ust(i)      = (umom**2.+vmom**2.)**.25
!
!  ELSE
!
!    sh_urb2d(i)    = 0.
!    lh_urb2d(i)    = 0.
!    g_urb2d(i)     = 0.
!    rn_urb2d(i)    = 0.
!
!  ENDIF
!
!  enddo ! iloop
!
!ENDIF ! SF_URBAN_PHYSICS == 2 or 3

!--------------------------------------
! URBAN BEP and BEM MODEL END
!--------------------------------------

  end subroutine noahmp_urban

end module NoahmpUrbanDriverMainMod
