module NoahmpWRFmainMod 

! -------------------------------------------------------------
! this is the interface for NoahMP and WRF variable remapping
! and calling the main NoahMP driver: NoahmpDriverMain(NoahmpIO)
! adapted from original module_sf_noahmpdrv.F file
!
! Coder: Cenlin He (NCAR), December 2025
! -------------------------------------------------------------

contains

  subroutine NoahmpWRFmain(NoahmpIO, ITIMESTEP, YR, JULIAN, COSZIN, XLAT, XLONG, & ! IN : Time/Space-related
                   DZ8W,          DT,        DZS,    NSOIL,       DX,            & ! IN : Model configuration 
                   IVGTYP,    ISLTYP,     VEGFRA,   VEGMAX,      TMN,            & ! IN : Vegetation/Soil characteristics
                   XLAND,       XICE, XICE_THRES,  CROPCAT,                      & ! IN : Vegetation/Soil characteristics
                   PLANTING, HARVEST, SEASON_GDD,                                & ! IN : agriculture input
                   IDVEG,   IOPT_CRS, IOPT_BTR, IOPT_RUNSUB,IOPT_SFC, IOPT_FRZ,  & ! IN : User options
                   IOPT_INF,IOPT_RAD,   IOPT_ALB, IOPT_SNF,IOPT_TBOT, IOPT_STC,  & ! IN : User options
                   IOPT_GLA,IOPT_RSF,  IOPT_SOIL,IOPT_PEDO,IOPT_CROP, IOPT_IRR,  & ! IN : User options
                   IOPT_IRRM,IOPT_INFDV,IOPT_TDRN, soiltstep,                    & ! IN : User options
                   IOPT_RUNSRF, IOPT_TKSNO, IOPT_COMPACT, IOPT_SCF, IOPT_WETLAND,& ! IN : User options                 
                   IZ0TLND, SF_URBAN_PHYSICS,                                    & ! IN : User options
                   SOILCOMP, SOILCL1,    SOILCL2,  SOILCL3,  SOILCL4,            & ! IN : User options
                   T3D, QV3D,  U_PHY,      V_PHY,   SWDOWN,   SWDDIR,            & ! IN : forcing
                   SWDDIF,       GLW,                                            & ! IN : Forcing
                   P8W3D,PRECIP_IN,           SR,                                & ! IN : Forcing
                   IRFRACT,  SIFRACT,    MIFRACT,  FIFRACT,                      & ! IN : irrigation
                   TSK,          HFX,        QFX,       LH,   GRDFLX,   SMSTAV,  & ! IN/OUT LSM eqv
                   SMSTOT, SFCRUNOFF,   UDRUNOFF,   ALBEDO,    SNOWC,    SMOIS,  & ! IN/OUT LSM eqv
                   SH2O,        TSLB,       SNOW,    SNOWH,   CANWAT,   ACSNOM,  & ! IN/OUT LSM eqv
                   ACSNOW,     EMISS,       QSFC,                                & ! IN/OUT LSM eqv
                   Z0,      ZNT,                                                 & ! IN/OUT LSM eqv
                   IRNUMSI,  IRNUMMI,    IRNUMFI,  IRWATSI,  IRWATMI,  IRWATFI,  & ! IN/OUT irrigation
                   IRELOSS,  IRSIVOL,    IRMIVOL,  IRFIVOL,  IRRSPLH, LLANDUSE,  & ! IN/OUT irrigation
                   ISNOWXY,     TVXY,       TGXY, CANICEXY, CANLIQXY,    EAHXY,  & ! IN/OUT Noah MP only
                   TAHXY,       CMXY,       CHXY,   FWETXY, SNEQVOXY, ALBOLDXY,  & ! IN/OUT Noah MP only
                   QSNOWXY,  QRAINXY,   WSLAKEXY,    ZWTXY,WAXY,WTXY,   TSNOXY,  & ! IN/OUT Noah MP only
                   ZSNSOXY,  SNICEXY,    SNLIQXY, LFMASSXY, RTMASSXY, STMASSXY,  & ! IN/OUT Noah MP only
                   WOODXY,  STBLCPXY,   FASTCPXY,   XLAIXY,   XSAIXY,  TAUSSXY,  & ! IN/OUT Noah MP only
                   SMOISEQ, SMCWTDXY, DEEPRECHXY,   RECHXY,GRAINXY,GDDXY,PGSXY,  & ! IN/OUT Noah MP only
                   QTDRAIN,   TD_FRACTION,                                       & ! IN/OUT tile drainage
                   T2MVXY,    T2MBXY,     Q2MVXY,   Q2MBXY,                      & ! OUT Noah MP only
                   TRADXY,     NEEXY,      GPPXY,    NPPXY,   FVEGXY,  RUNSFXY,  & ! OUT Noah MP only
                   RUNSBXY,   ECANXY,     EDIRXY,  ETRANXY,    FSAXY,   FIRAXY,  & ! OUT Noah MP only
                   APARXY,     PSNXY,      SAVXY,    SAGXY,  RSSUNXY,  RSSHAXY,  & ! OUT Noah MP only
                   BGAPXY,    WGAPXY,      TGVXY,    TGBXY,    CHVXY,    CHBXY,  & ! OUT Noah MP only
                   SHGXY,      SHCXY,      SHBXY,    EVGXY,    EVBXY,    GHVXY,  & ! OUT Noah MP only
                   GHBXY,      IRGXY,      IRCXY,    IRBXY,     TRXY,    EVCXY,  & ! OUT Noah MP only
                   CHLEAFXY,  CHUCXY,     CHV2XY,   CHB2XY,       RS,            & ! OUT Noah MP only
                   QINTSXY,  QINTRXY,   QDRIPSXY,                                & ! OUT Noah MP only
                   QDRIPRXY,QTHROSXY,   QTHRORXY,                                & ! OUT Noah MP only
                   QSNSUBXY,QSNFROXY,    QSUBCXY,                                & ! OUT Noah MP only
                   QFROCXY,  QEVACXY,    QDEWCXY,  QFRZCXY, QMELTCXY,            & ! OUT Noah MP only
                   QSNBOTXY, QMELTXY,  PONDINGXY,  PAHXY,PAHGXY,PAHVXY, PAHBXY,  & ! OUT Noah MP only
                   FPICEXY,RAINLSM,SNOWLSM,FORCTLSM,FORCQLSM,FORCPLSM,FORCZLSM,  & ! OUT Noah MP only
                   FORCWLSM,ACC_SSOILXY,ACC_QINSURXY,ACC_QSEVAXY, ACC_ETRANIXY,  & ! IN/OUT Noah MP
                   EFLXBXY, SOILENERGY, SNOWENERGY, CANHSXY,                     & ! OUT Noah MP only
                   ACC_DWATERXY, ACC_PRCPXY, ACC_ECANXY,ACC_ETRANXY,ACC_EDIRXY,  & ! IN/OUT Noah MP
                   FSATXY, WSURFXY,                                              & ! IN/OUT Noah MP
                   SNICAR_BANDNUMBER_OPT, SNICAR_SOLARSPEC_OPT,                  & ! SNICAR variable
                   SNICAR_SNOWOPTICS_OPT, SNICAR_DUSTOPTICS_OPT,                 & ! SNICAR variable
                   SNICAR_RTSOLVER_OPT, SNICAR_SNOWSHAPE_OPT,                    & ! SNICAR variable
                   SNICAR_USE_AEROSOL, SNICAR_SNOWBC_INTMIX,                     & ! SNICAR variable
                   SNICAR_SNOWDUST_INTMIX, SNICAR_USE_OC,                        & ! SNICAR variable
                   SNICAR_AEROSOL_READTABLE, SNRDSXY, SNFRXY, BCPHIXY, BCPHOXY,  & ! SNICAR variable
                   OCPHIXY, OCPHOXY, DUST1XY, DUST2XY, DUST3XY, DUST4XY, DUST5XY,& ! SNICAR variable
                   MassConcBCPHIXY, MassConcBCPHOXY, MassConcOCPHIXY,            & ! SNICAR variable
                   MassConcOCPHOXY, MassConcDUST1XY, MassConcDUST2XY,            & ! SNICAR variable
                   MassConcDUST3XY, MassConcDUST4XY, MassConcDUST5XY,            & ! SNICAR variable
                   ALBSOILDIRXY, ALBSOILDIFXY,                                   & ! SNICAR variable
                 ! BEXP_3D,SMCDRY_3D,SMCWLT_3D,SMCREF_3D,SMCMAX_3D,              & ! placeholders to activate 3D soil
                 ! DKSAT_3D,DWSAT_3D,PSISAT_3D,QUARTZ_3D,                        & ! placeholders to activate 3D soil
                 ! REFDK_2D,REFKDT_2D,                                           & ! placeholders to activate 3D soil
                 ! IRR_FRAC_2D,IRR_HAR_2D,IRR_LAI_2D,IRR_MAD_2D,FILOSS_2D,       & ! placeholders to activate 3D soil
                 ! SPRIR_RATE_2D,MICIR_RATE_2D,FIRTFAC_2D,IR_RAIN_2D,            & ! placeholders to activate 3D soil
                 ! BVIC_2D,AXAJ_2D,BXAJ_2D,XXAJ_2D,BDVIC_2D,GDVIC_2D,BBVIC_2D,   & ! placeholders to activate 3D soil
                 ! KLAT_FAC,TDSMC_FAC,TD_DC,TD_DCOEF,TD_DDRAIN,TD_RADI,TD_SPAC,  & ! placeholders to activate 3D soil
#ifdef WRF_HYDRO
                   sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D,              & ! OUT WRF-Hydro only
#endif
                   ids,ide,  jds,jde,  kds,kde,                                  & ! IN: WRF dimension
                   ims,ime,  jms,jme,  kms,kme,                                  & ! IN: WRF dimension
                   its,ite,  jts,jte,  kts,kte,                                  & ! IN: WRF dimension
                   MP_RAINC,MP_RAINNC,MP_SHCV,MP_SNOW,MP_GRAUP,MP_HAIL           ) ! IN: WRF forcing

!----------------------------------------------------------------

    use NoahmpIOVarType
    use NoahmpIOVarInitMod
    use NoahmpReadTableMod
    use SnowInputSnicarMod
    use NoahmpDriverMainMod
    use module_sf_urban,    only: IRI_SCHEME

    implicit none

    type(NoahmpIO_type),                             intent(inout) ::  NoahmpIO

    ! input
    INTEGER,                                         INTENT(IN   ) ::  ids,ide, jds,jde, kds,kde,  &  ! d -> domain
                                                                       ims,ime, jms,jme, kms,kme,  &  ! m -> memory
                                                                       its,ite, jts,jte, kts,kte      ! t -> tile
    INTEGER,                                         INTENT(IN   ) ::  ITIMESTEP    ! timestep number
    INTEGER,                                         INTENT(IN   ) ::  YR           ! 4-digit year
    REAL,                                            INTENT(IN   ) ::  JULIAN       ! Julian day
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  COSZIN       ! cosine zenith angle
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAT         ! latitude [rad]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLONG        ! latitude [rad]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  DZ8W         ! thickness of atmo layers [m]
    REAL,                                            INTENT(IN   ) ::  DT           ! timestep [s]
    REAL,    DIMENSION(1:nsoil),                     INTENT(IN   ) ::  DZS          ! thickness of soil layers [m]
    INTEGER,                                         INTENT(IN   ) ::  NSOIL        ! number of soil layers
    REAL,                                            INTENT(IN   ) ::  DX           ! horizontal grid spacing [m]
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  IVGTYP       ! vegetation type
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  ISLTYP       ! soil type
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  VEGFRA       ! vegetation fraction []
    REAL,    DIMENSION( ims:ime ,         jms:jme ), INTENT(IN   ) ::  VEGMAX       ! annual max vegetation fraction []
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  TMN          ! deep soil temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XLAND        ! =2 ocean; =1 land/seaice
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  XICE         ! fraction of grid that is seaice
    REAL,                                            INTENT(IN   ) ::  XICE_THRES   ! fraction of grid determining seaice
    INTEGER,                                         INTENT(IN   ) ::  IDVEG        ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1      
    INTEGER,                                         INTENT(IN   ) ::  IOPT_CRS     ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_BTR     ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RUNSUB  ! subsurface runoff and groundwater (currently keep the same as surface runoff option)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RUNSRF  ! surface runoff (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XianAnJiang; 8->DynVIC)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_COMPACT ! snowpack compaction (1->Anderson1976; 2->Abolafia-Rosenzweig2024)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_TKSNO   ! snow thermal conductivity: 1 -> Stieglitz(yen,1965) scheme (default), 2 -> Anderson, 1976 scheme, 3 -> constant, 4 -> Verseghy (1991) scheme, 5 -> Douvill(Yen, 1981) scheme
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SCF     ! snow cover fraction (1->NiuYang07; 2->Abolafia-Rosenzweig2025)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_WETLAND ! wetland model option (0->off; 1->Zhang2022 fixed parameter; 2->Zhang2022 read in 2D parameter)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SFC     ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_FRZ     ! supercooled liquid water (1-> NY06; 2->Koren99)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_INF     ! frozen soil permeability (1-> NY06; 2->Koren99)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RAD     ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_ALB     ! snow surface albedo (1->BATS; 2->CLASS; 3->SNICAR)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SNF     ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_TBOT    ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_STC     ! snow/soil temperature time scheme
    INTEGER,                                         INTENT(IN   ) ::  IOPT_GLA     ! glacier option (1->phase change; 2->simple)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RSF     ! surface resistance (1->Sakaguchi/Zeng; 2->Seller; 3->mod Sellers; 4->1+snow)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SOIL    ! soil configuration option
    INTEGER,                                         INTENT(IN   ) ::  IOPT_PEDO    ! soil pedotransfer function option
    INTEGER,                                         INTENT(IN   ) ::  IOPT_CROP    ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_IRR     ! irrigation scheme (0->none; >1 irrigation scheme ON)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_IRRM    ! irrigation method
    INTEGER,                                         INTENT(IN   ) ::  IOPT_INFDV   ! infiltration options for dynamic VIC infiltration (1->Philip; 2-> Green-Ampt;3->Smith-Parlange)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_TDRN    ! tile drainage (0-> no tile drainage; 1-> simple tile drainage;2->Hooghoudt's)
    REAL,                                            INTENT(IN   ) ::  soiltstep    ! soil timestep (s), default:0->same as main model timestep
    INTEGER,                                         INTENT(IN   ) ::  IZ0TLND      ! option of Chen adjustment of Czil (not used)
    INTEGER,                                         INTENT(IN   ) ::  sf_urban_physics ! urban physics option
    REAL,    DIMENSION( ims:ime, nsoil*2, jms:jme ), INTENT(IN   ) ::  SOILCOMP     ! soil sand and clay percentage
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL1      ! soil texture in layer 1
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL2      ! soil texture in layer 2
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL3      ! soil texture in layer 3
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SOILCL4      ! soil texture in layer 4
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  T3D          ! 3D atmospheric temperature valid at mid-levels [K]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  QV3D         ! 3D water vapor mixing ratio [kg/kg_dry]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  U_PHY        ! 3D U wind component [m/s]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  V_PHY        ! 3D V wind component [m/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDOWN       ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIF       ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SWDDIR       ! solar down at surface [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  GLW          ! longwave down at surface [W m-2]
    REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::  P8W3D        ! 3D pressure, valid at interface [Pa]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PRECIP_IN    ! total input precipitation [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SR           ! frozen precipitation ratio [-]
    ! Optional Detailed Precipitation Partitioning Inputs
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINC  ! convective precipitation entering land model [mm] ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_RAINNC ! large-scale precipitation entering land model [mm]! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SHCV   ! shallow conv precip entering land model [mm]      ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_SNOW   ! snow precipitation entering land model [mm]       ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_GRAUP  ! graupel precipitation entering land model [mm]    ! MB/AN : v3.7
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ), OPTIONAL ::  MP_HAIL   ! hail precipitation entering land model [mm]       ! MB/AN : v3.7
    ! Crop Model
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  CROPCAT     ! crop catagory
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  PLANTING    ! planting date
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  HARVEST     ! harvest date
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SEASON_GDD  ! growing season GDD
    ! Tile drain variables    
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  TD_FRACTION
    !2D inout irrigation variables 
    CHARACTER(LEN=256),                              INTENT(IN   ) ::  LLANDUSE    ! landuse data name (USGS or MODIS_IGBP)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  IRFRACT     ! irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  SIFRACT     ! sprinkler irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  MIFRACT     ! micro irrigation fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) ::  FIFRACT     ! flood irrigation fraction   
    ! SNICAR snow albedo variables
    INTEGER,                                         INTENT(IN   ) ::  SNICAR_BANDNUMBER_OPT    ! number of wavelength bands used in SNICAR
    INTEGER,                                         INTENT(IN   ) ::  SNICAR_SOLARSPEC_OPT     ! type of downward solar radiation spectrum for SNICAR
    INTEGER,                                         INTENT(IN   ) ::  SNICAR_SNOWOPTICS_OPT    ! snow optics type using different refractive index databases in SNICAR
    INTEGER,                                         INTENT(IN   ) ::  SNICAR_DUSTOPTICS_OPT    ! dust optics type for SNICAR snow albedo calculation
    INTEGER,                                         INTENT(IN   ) ::  SNICAR_RTSOLVER_OPT      ! option for two different SNICAR radiative transfer solver
    INTEGER,                                         INTENT(IN   ) ::  SNICAR_SNOWSHAPE_OPT     ! option for snow grain shape in SNICAR
    LOGICAL,                                         INTENT(IN   ) ::  SNICAR_USE_AEROSOL       ! option to turn on/off aerosol deposition flux effect in snow in SNICAR
    LOGICAL,                                         INTENT(IN   ) ::  SNICAR_SNOWBC_INTMIX     ! option to activate BC-snow internal mixing in SNICAR
    LOGICAL,                                         INTENT(IN   ) ::  SNICAR_SNOWDUST_INTMIX   ! option to activate dust-snow internal mixing in SNICAR 
    LOGICAL,                                         INTENT(IN   ) ::  SNICAR_USE_OC            ! option to activate OC in snow in SNICAR
    LOGICAL,                                         INTENT(IN   ) ::  SNICAR_AEROSOL_READTABLE ! option to read aerosol deposition fluxes from table (on) or NetCDF forcing file (off)

    ! placeholders for 2D/3D soil
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  BEXP_3D       ! C-H B exponent
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCDRY_3D     ! Soil Moisture Limit: Dry
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCWLT_3D     ! Soil Moisture Limit: Wilt
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCREF_3D     ! Soil Moisture Limit: Reference
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  SMCMAX_3D     ! Soil Moisture Limit: Max
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  DKSAT_3D      ! Saturated Soil Conductivity
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  DWSAT_3D      ! Saturated Soil Diffusivity
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  PSISAT_3D     ! Saturated Matric Potential
    ! REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(IN) ::  QUARTZ_3D     ! Soil quartz content
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFDK_2D      ! Reference Soil Conductivity
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFKDT_2D     ! Soil Infiltration Parameter
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BVIC_2D       ! VIC model infiltration parameter [-] for opt_run=6
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  AXAJ_2D       ! Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BXAJ_2D       ! Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  XXAJ_2D       ! Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BDVIC_2D      ! VIC model infiltration parameter [-]  
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  GDVIC_2D      ! Mean Capillary Drive (m) for infiltration models
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  BBVIC_2D      ! DVIC heterogeniety paramater [-]

    ! placeholders for 2D irrigation parameters
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_FRAC_2D   ! irrigation Fraction
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_HAR_2D    ! number of days before harvest date to stop irrigation 
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_LAI_2D    ! Minimum lai to trigger irrigation
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IRR_MAD_2D    ! management allowable deficit (0-1)
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  FILOSS_2D     ! fraction of flood irrigation loss (0-1) 
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  SPRIR_RATE_2D ! mm/h, sprinkler irrigation rate
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  MICIR_RATE_2D ! mm/h, micro irrigation rate
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  FIRTFAC_2D    ! flood application rate factor
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  IR_RAIN_2D    ! maximum precipitation to stop irrigation trigger

    ! placeholders for 2D tile drainage parameters
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  KLAT_FAC      ! factor multiplier to hydraulic conductivity
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TDSMC_FAC     ! factor multiplier to field capacity
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DC         ! drainage coefficient for simple
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DCOEF      ! drainge coefficient for Hooghoudt 
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_DDRAIN     ! depth of drain
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_RADI       ! tile radius
    ! REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  TD_SPAC       ! tile spacing

    ! INOUT (with generic LSM equivalent)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TSK          ! surface radiative temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  HFX          ! sensible heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QFX          ! latent heat flux [kg s-1 m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LH           ! latent heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRDFLX       ! ground/snow heat flux [W m-2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTAV       ! soil moisture avail. [not used]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMSTOT       ! total soil water [mm][not used]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SFCRUNOFF    ! accumulated surface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  UDRUNOFF     ! accumulated sub-surface runoff [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBEDO       ! total grid albedo []
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWC        ! snow cover fraction []
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOIS        ! volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SH2O         ! volumetric liquid soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  TSLB         ! soil temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOW         ! snow water equivalent [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNOWH        ! physical snow depth [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANWAT       ! total canopy water + ice [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOM       ! accumulated snow melt (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACSNOW       ! accumulated snow on grid
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EMISS        ! surface bulk emissivity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSFC         ! bulk surface specific humidity
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  Z0           ! combined z0 sent to coupled model
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZNT          ! combined z0 sent to coupled model
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RS           ! Total stomatal resistance (s/m)
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ISNOWXY      ! actual no. of snow layers
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TVXY         ! vegetation leaf temperature
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TGXY         ! bulk ground surface temperature
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANICEXY     ! canopy-intercepted ice (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CANLIQXY     ! canopy-intercepted liquid water (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  EAHXY        ! canopy air vapor pressure (pa)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAHXY        ! canopy air temperature (k)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CMXY         ! bulk momentum drag coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  CHXY         ! bulk sensible heat exchange coefficient
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FWETXY       ! wetted or snowed fraction of the canopy (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SNEQVOXY     ! snow mass at last time step(mm h2o)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ALBOLDXY     ! snow albedo at last time step (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QSNOWXY      ! snowfall on the ground [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QRAINXY      ! rainfall on the ground [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WSLAKEXY     ! lake water storage [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ZWTXY        ! water table depth [m]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WAXY         ! water in the "aquifer" [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WTXY         ! groundwater storage [mm]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  TSNOXY       ! snow temperature [K]
    REAL,    DIMENSION( ims:ime,-2:NSOIL, jms:jme ), INTENT(INOUT) ::  ZSNSOXY      ! snow layer depth [m]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNICEXY      ! snow layer ice [mm]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNLIQXY      ! snow layer liquid water [mm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  LFMASSXY     ! leaf mass [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RTMASSXY     ! mass of fine roots [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STMASSXY     ! stem mass [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WOODXY       ! mass of wood (incl. woody roots) [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  STBLCPXY     ! stable carbon in deep soil [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FASTCPXY     ! short-lived carbon, shallow soil [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XLAIXY       ! leaf area index
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  XSAIXY       ! stem area index
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  TAUSSXY      ! snow age factor
    REAL,    DIMENSION( ims:ime, 1:nsoil, jms:jme ), INTENT(INOUT) ::  SMOISEQ      ! eq volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  SMCWTDXY     ! soil moisture content in the layer to the water table when deep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  DEEPRECHXY   ! recharge to the water table when deep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  RECHXY       ! recharge to the water table (diagnostic) 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GRAINXY      ! mass of grain XING [g/m2]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  GDDXY        ! growing degree days XING (based on 10C) 
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  PGSXY        ! growing stage
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_SSOILXY  ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QINSURXY ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_QSEVAXY  ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime, 1:NSOIL, jms:jme ), INTENT(INOUT) ::  ACC_ETRANIXY ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_DWATERXY ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_PRCPXY   ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ECANXY   ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_ETRANXY  ! m/s * soil_dt/main_dt
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  ACC_EDIRXY   ! m/s * soil_dt/main_dt
    !2D inout irrigation variables 
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRNUMSI      ! irrigation event number, Sprinkler
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRNUMMI      ! irrigation event number, Micro
    INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRNUMFI      ! irrigation event number, Flood 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRWATSI      ! irrigation water amount [m] to be applied, Sprinkler
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRWATMI      ! irrigation water amount [m] to be applied, Micro
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRWATFI      ! irrigation water amount [m] to be applied, Flood
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRELOSS      ! loss of irrigation water to evaporation,sprinkler [m/timestep]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRSIVOL      ! amount of irrigation by sprinkler (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRMIVOL      ! amount of irrigation by micro (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRFIVOL      ! amount of irrigation by micro (mm)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  IRRSPLH      ! latent heating from sprinkler evaporation (w/m2)    
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  QTDRAIN      ! Tile drain
    ! wetland state varible
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  FSATXY       ! saturated fraction of the grid (-)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  WSURFXY      ! wetland water storage [mm]
    ! SNICAR state variable
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNRDSXY         ! snow layer effective grain radius [microns, m-6]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  SNFRXY          ! snow layer rate of snow freezing [mm/s]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  BCPHIXY         ! mass of hydrophillic Black Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  BCPHOXY         ! mass of hydrophobic Black Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  OCPHIXY         ! mass of hydrophillic Organic Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  OCPHOXY         ! mass of hydrophobic Organic Carbon in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST1XY         ! mass of dust species 1 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST2XY         ! mass of dust species 2 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST3XY         ! mass of dust species 3 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST4XY         ! mass of dust species 4 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  DUST5XY         ! mass of dust species 5 in snow [kg/m2]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcBCPHIXY ! mass concentration of hydrophillic Black Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcBCPHOXY ! mass concentration of hydrophobic Black Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcOCPHIXY ! mass concentration of hydrophillic Organic Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcOCPHOXY ! mass concentration of hydrophobic Organic Carbon in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST1XY ! mass concentration of dust species 1 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST2XY ! mass concentration of dust species 2 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST3XY ! mass concentration of dust species 3 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST4XY ! mass concentration of dust species 4 in snow [kg/kg]
    REAL,    DIMENSION( ims:ime,-2:0,     jms:jme ), INTENT(INOUT) ::  MassConcDUST5XY ! mass concentration of dust species 5 in snow [kg/kg]
    REAL,    DIMENSION(ims:ime,1:2,       jms:jme),  INTENT(INOUT) ::  ALBSOILDIRXY    ! soil albedo direct
    REAL,    DIMENSION(ims:ime,1:2,       jms:jme),  INTENT(INOUT) ::  ALBSOILDIFXY    ! soil albedo diffuse
#ifdef WRF_HYDRO
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) ::  sfcheadrt,INFXSRT,soldrain,qtiledrain,ZWATBLE2D   ! for WRF-Hydro
#endif

    ! OUT (with no Noah LSM equivalent)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MVXY       ! 2m temperature of vegetation part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  T2MBXY       ! 2m temperature of bare ground part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MVXY       ! 2m mixing ratio of vegetation part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  Q2MBXY       ! 2m mixing ratio of bare ground part
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRADXY       ! surface radiative temperature (k)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NEEXY        ! net ecosys exchange (g/m2/s CO2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GPPXY        ! gross primary assimilation [g/m2/s C]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  NPPXY        ! net primary productivity [g/m2/s C]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FVEGXY       ! Noah-MP vegetation fraction [-]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSFXY      ! surface runoff [mm] per soil timestep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RUNSBXY      ! subsurface runoff [mm] per soil timestep
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ECANXY       ! evaporation of intercepted water (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EDIRXY       ! soil surface evaporation rate (mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  ETRANXY      ! transpiration rate (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FSAXY        ! total absorbed solar radiation (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FIRAXY       ! total net longwave rad (w/m2) [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  APARXY       ! photosyn active energy by canopy (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PSNXY        ! total photosynthesis (umol co2/m2/s) [+]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAVXY        ! solar rad absorbed by veg. (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SAGXY        ! solar rad absorbed by ground (w/m2)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSUNXY      ! sunlit leaf stomatal resistance (s/m)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RSSHAXY      ! shaded leaf stomatal resistance (s/m)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  BGAPXY       ! between gap fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  WGAPXY       ! within gap fraction
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGVXY        ! under canopy ground temperature[K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TGBXY        ! bare ground temperature [K]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHVXY        ! sensible heat exchange coefficient vegetated
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHBXY        ! sensible heat exchange coefficient bare-ground
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHGXY        ! veg ground sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHCXY        ! canopy sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SHBXY        ! bare sensible heat [w/m2]     [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVGXY        ! veg ground evap. heat [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVBXY        ! bare soil evaporation [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHVXY        ! veg ground heat flux [w/m2]  [+ to soil]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  GHBXY        ! bare ground heat flux [w/m2] [+ to soil]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRGXY        ! veg ground net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRCXY        ! canopy net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  IRBXY        ! bare net longwave rad. [w/m2] [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  TRXY         ! transpiration [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EVCXY        ! canopy evaporation heat [w/m2]  [+ to atm]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHLEAFXY     ! leaf exchange coefficient 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHUCXY       ! under canopy exchange coefficient 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHV2XY       ! veg 2m exchange coefficient 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CHB2XY       ! bare 2m exchange coefficient
    ! additional output variables
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHXY        ! precipitation advected heat 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHGXY       ! precipitation advected heat 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHBXY       ! precipitation advected heat 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PAHVXY       ! precipitation advected heat 
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTSXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QINTRXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPSXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDRIPRXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHROSXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QTHRORXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNSUBXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNFROXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSUBCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFROCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QEVACXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QDEWCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QFRZCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTCXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QSNBOTXY     ! total liquid water (snowmelt + rain through pack)out of snowpack bottom [mm/s]
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  QMELTXY      ! snowmelt due to phase change (mm/s)
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  PONDINGXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FPICEXY      ! fraction of ice in precip
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  RAINLSM      ! rain rate                   (mm/s)  AJN
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWLSM      ! liquid equivalent snow rate (mm/s)  AJN
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCTLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCQLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCPLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCZLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  FORCWLSM
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  EFLXBXY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SOILENERGY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  SNOWENERGY
    REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(OUT  ) ::  CANHSXY

! ----------------------------------------------------------------------

    !--------- Input variable mapping start ---------

    ! input WRF variables mapped to NoahmpIO variables
    NoahmpIO%ids                = ids
    NoahmpIO%ide                = ide
    NoahmpIO%jds                = jds
    NoahmpIO%jde                = jde
    NoahmpIO%kds                = kds
    NoahmpIO%kde                = kde
    NoahmpIO%ims                = ims
    NoahmpIO%ime                = ime
    NoahmpIO%jms                = jms
    NoahmpIO%jme                = jme
    NoahmpIO%kms                = kms
    NoahmpIO%kme                = kme
    NoahmpIO%its                = its
    NoahmpIO%ite                = ite
    NoahmpIO%jts                = jts
    NoahmpIO%jte                = jte
    NoahmpIO%kts                = kts
    NoahmpIO%kte                = kte
    NoahmpIO%xstart             = ims
    NoahmpIO%xend               = ime
    NoahmpIO%ystart             = jms
    NoahmpIO%yend               = jme    
    NoahmpIO%YR                 = YR
    NoahmpIO%JULIAN             = JULIAN
    NoahmpIO%DTBL               = DT
    NoahmpIO%NSOIL              = NSOIL
    NoahmpIO%DX                 = DX
    NoahmpIO%DY                 = DX
    NoahmpIO%IOPT_DVEG          = IDVEG
    NoahmpIO%IOPT_CRS           = IOPT_CRS
    NoahmpIO%IOPT_BTR           = IOPT_BTR
    NoahmpIO%IOPT_SFC           = IOPT_SFC
    NoahmpIO%IOPT_FRZ           = IOPT_FRZ
    NoahmpIO%IOPT_INF           = IOPT_INF
    NoahmpIO%IOPT_RAD           = IOPT_RAD
    NoahmpIO%IOPT_ALB           = IOPT_ALB
    NoahmpIO%IOPT_SNF           = IOPT_SNF
    NoahmpIO%IOPT_TBOT          = IOPT_TBOT
    NoahmpIO%IOPT_STC           = IOPT_STC
    NoahmpIO%IOPT_GLA           = IOPT_GLA
    NoahmpIO%IOPT_RSF           = IOPT_RSF
    NoahmpIO%IOPT_SOIL          = IOPT_SOIL
    NoahmpIO%IOPT_PEDO          = IOPT_PEDO
    NoahmpIO%IOPT_CROP          = IOPT_CROP
    NoahmpIO%IOPT_IRR           = IOPT_IRR
    NoahmpIO%IOPT_IRRM          = IOPT_IRRM
    NoahmpIO%IOPT_INFDV         = IOPT_INFDV
    NoahmpIO%IOPT_TDRN          = IOPT_TDRN
    NoahmpIO%IOPT_RUNSRF        = IOPT_RUNSRF
    NoahmpIO%IOPT_RUNSUB        = IOPT_RUNSUB
    NoahmpIO%IOPT_TKSNO         = IOPT_TKSNO
    NoahmpIO%IOPT_COMPACT       = IOPT_COMPACT
    NoahmpIO%IOPT_SCF           = IOPT_SCF
    NoahmpIO%IOPT_WETLAND       = IOPT_WETLAND
    NoahmpIO%SF_URBAN_PHYSICS   = SF_URBAN_PHYSICS
    NoahmpIO%IZ0TLND            = IZ0TLND
    NoahmpIO%LLANDUSE           = LLANDUSE
    NoahmpIO%SOILTSTEP          = SOILTSTEP
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%SNICAR_BANDNUMBER_OPT    = SNICAR_BANDNUMBER_OPT
       NoahmpIO%SNICAR_SOLARSPEC_OPT     = SNICAR_SOLARSPEC_OPT
       NoahmpIO%SNICAR_SNOWOPTICS_OPT    = SNICAR_SNOWOPTICS_OPT
       NoahmpIO%SNICAR_DUSTOPTICS_OPT    = SNICAR_DUSTOPTICS_OPT
       NoahmpIO%SNICAR_RTSOLVER_OPT      = SNICAR_RTSOLVER_OPT
       NoahmpIO%SNICAR_SNOWSHAPE_OPT     = SNICAR_SNOWSHAPE_OPT
       NoahmpIO%SNICAR_USE_AEROSOL       = SNICAR_USE_AEROSOL
       NoahmpIO%SNICAR_SNOWBC_INTMIX     = SNICAR_SNOWBC_INTMIX
       NoahmpIO%SNICAR_SNOWDUST_INTMIX   = SNICAR_SNOWDUST_INTMIX
       NoahmpIO%SNICAR_USE_OC            = SNICAR_USE_OC
       NoahmpIO%SNICAR_AEROSOL_READTABLE = SNICAR_AEROSOL_READTABLE
    endif

    ! initialze all NoahmpIO variables with default values
    call NoahmpIOVarInitDefault(NoahmpIO)

    ! read in Noahmp table parameters
    call NoahmpReadTable(NoahmpIO)

    ! read in SNICAR parameter netcdif file
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%snicar_optic_flnm = "snicar_optics_5bnd_c013122.nc"
       NoahmpIO%snicar_age_flnm   = "snicar_drdt_bst_fit_60_c070416.nc"
       call SnowInputSnicar(NoahmpIO)
    endif

    if ( NoahmpIO%IOPT_SOIL > 1 ) then
       NoahmpIO%SOILCOMP        = SOILCOMP
       NoahmpIO%SOILCL1         = SOILCL1
       NoahmpIO%SOILCL2         = SOILCL2
       NoahmpIO%SOILCL3         = SOILCL3
       NoahmpIO%SOILCL4         = SOILCL4
    endif
    NoahmpIO%IVGTYP             = IVGTYP
    NoahmpIO%ISLTYP             = ISLTYP
    NoahmpIO%VEGFRA             = VEGFRA
    NoahmpIO%GVFMAX             = VEGMAX
    NoahmpIO%TMN                = TMN
    NoahmpIO%XLAND              = XLAND
    NoahmpIO%XICE               = XICE
    NoahmpIO%XICE_THRESHOLD     = XICE_THRES
    NoahmpIO%CROPCAT            = CROPCAT
    NoahmpIO%PLANTING           = PLANTING
    NoahmpIO%HARVEST            = HARVEST
    NoahmpIO%SEASON_GDD         = SEASON_GDD
    NoahmpIO%DZS                = DZS
    NoahmpIO%XLAT               = XLAT
    NoahmpIO%XLONG              = XLONG
    NoahmpIO%COSZEN             = COSZIN
    NoahmpIO%IRI_URBAN          = IRI_SCHEME
    NoahmpIO%ITIMESTEP          = ITIMESTEP
    NoahmpIO%DZ8W               = DZ8W
    NoahmpIO%T_PHY              = T3D
    NoahmpIO%QV_CURR            = QV3D
    NoahmpIO%U_PHY              = U_PHY
    NoahmpIO%V_PHY              = V_PHY
    NoahmpIO%SWDOWN             = SWDOWN
    NoahmpIO%SWDDIR             = SWDDIR
    NoahmpIO%SWDDIF             = SWDDIF
    NoahmpIO%GLW                = GLW
    NoahmpIO%P8W                = P8W3D
    NoahmpIO%RAINBL             = PRECIP_IN
    NoahmpIO%SR                 = SR
    NoahmpIO%IRFRACT            = IRFRACT
    NoahmpIO%SIFRACT            = SIFRACT
    NoahmpIO%MIFRACT            = MIFRACT
    NoahmpIO%FIFRACT            = FIFRACT
    NoahmpIO%TD_FRACTION        = TD_FRACTION
    if (present(MP_RAINC) .and. present(MP_RAINNC) .and. &
        present(MP_SHCV)  .and. present(MP_SNOW)   .and. &
        present(MP_GRAUP) .and. present(MP_HAIL) ) then
       NoahmpIO%MP_RAINC        = MP_RAINC
       NoahmpIO%MP_RAINNC       = MP_RAINNC
       NoahmpIO%MP_SHCV         = MP_SHCV
       NoahmpIO%MP_SNOW         = MP_SNOW
       NoahmpIO%MP_GRAUP        = MP_GRAUP
       NoahmpIO%MP_HAIL         = MP_HAIL
    endif
    ! NoahmpIO%BEXP_3D          = BEXP_3D
    ! NoahmpIO%SMCDRY_3D        = SMCDRY_3D
    ! NoahmpIO%SMCWLT_3D        = SMCWLT_3D
    ! NoahmpIO%SMCREF_3D        = SMCREF_3D
    ! NoahmpIO%SMCMAX_3D        = SMCMAX_3D
    ! NoahmpIO%DKSAT_3D         = DKSAT_3D
    ! NoahmpIO%DWSAT_3D         = DWSAT_3D
    ! NoahmpIO%PSISAT_3D        = PSISAT_3D
    ! NoahmpIO%QUARTZ_3D        = QUARTZ_3D
    ! NoahmpIO%REFDK_2D         = REFDK_2D
    ! NoahmpIO%REFKDT_2D        = REFKDT_2D
    ! NoahmpIO%IRR_FRAC_2D      = IRR_FRAC_2D
    ! NoahmpIO%IRR_HAR_2D       = IRR_HAR_2D
    ! NoahmpIO%IRR_LAI_2D       = IRR_LAI_2D
    ! NoahmpIO%IRR_MAD_2D       = IRR_MAD_2D
    ! NoahmpIO%FILOSS_2D        = FILOSS_2D
    ! NoahmpIO%SPRIR_RATE_2D    = SPRIR_RATE_2D
    ! NoahmpIO%MICIR_RATE_2D    = MICIR_RATE_2D
    ! NoahmpIO%FIRTFAC_2D       = FIRTFAC_2D
    ! NoahmpIO%IR_RAIN_2D       = IR_RAIN_2D
    ! NoahmpIO%BVIC_2D          = BVIC_2D
    ! NoahmpIO%AXAJ_2D          = AXAJ_2D
    ! NoahmpIO%BXAJ_2D          = BXAJ_2D
    ! NoahmpIO%XXAJ_2D          = XXAJ_2D
    ! NoahmpIO%BDVIC_2D         = BDVIC_2D
    ! NoahmpIO%GDVIC_2D         = GDVIC_2D
    ! NoahmpIO%BBVIC_2D         = BBVIC_2D
    ! NoahmpIO%KLAT_FAC         = KLAT_FAC
    ! NoahmpIO%TDSMC_FAC        = TDSMC_FAC
    ! NoahmpIO%TD_DC            = TD_DC
    ! NoahmpIO%TD_DCOEF         = TD_DCOEF
    ! NoahmpIO%TD_DDRAIN        = TD_DDRAIN
    ! NoahmpIO%TD_RADI          = TD_RADI
    ! NoahmpIO%TD_SPAC          = TD_SPAC
    
    ! in/out WRF variables mapped to NoahmpIO variables
    NoahmpIO%TSK                = TSK
    NoahmpIO%HFX                = HFX
    NoahmpIO%QFX                = QFX
    NoahmpIO%LH                 = LH
    NoahmpIO%GRDFLX             = GRDFLX
    NoahmpIO%SMSTAV             = SMSTAV
    NoahmpIO%SMSTOT             = SMSTOT
    NoahmpIO%SFCRUNOFF          = SFCRUNOFF
    NoahmpIO%UDRUNOFF           = UDRUNOFF
    NoahmpIO%ALBEDO             = ALBEDO
    NoahmpIO%SNOWC              = SNOWC
    NoahmpIO%SMOIS              = SMOIS
    NoahmpIO%SH2O               = SH2O
    NoahmpIO%TSLB               = TSLB
    NoahmpIO%SNOW               = SNOW
    NoahmpIO%SNOWH              = SNOWH
    NoahmpIO%CANWAT             = CANWAT
    NoahmpIO%CANICEXY           = CANICEXY
    NoahmpIO%CANLIQXY           = CANLIQXY
    NoahmpIO%ACSNOM             = ACSNOM
    NoahmpIO%ACSNOW             = ACSNOW
    NoahmpIO%EMISS              = EMISS
    NoahmpIO%QSFC               = QSFC
    NoahmpIO%Z0                 = Z0
    NoahmpIO%ZNT                = ZNT
    NoahmpIO%IRNUMSI            = IRNUMSI
    NoahmpIO%IRNUMMI            = IRNUMMI
    NoahmpIO%IRNUMFI            = IRNUMFI
    NoahmpIO%IRWATSI            = IRWATSI
    NoahmpIO%IRWATMI            = IRWATMI
    NoahmpIO%IRWATFI            = IRWATFI
    NoahmpIO%IRELOSS            = IRELOSS
    NoahmpIO%IRSIVOL            = IRSIVOL
    NoahmpIO%IRMIVOL            = IRMIVOL
    NoahmpIO%IRFIVOL            = IRFIVOL
    NoahmpIO%IRRSPLH            = IRRSPLH
    NoahmpIO%ISNOWXY            = ISNOWXY
    NoahmpIO%TVXY               = TVXY
    NoahmpIO%TGXY               = TGXY
    NoahmpIO%EAHXY              = EAHXY
    NoahmpIO%TAHXY              = TAHXY
    NoahmpIO%CMXY               = CMXY
    NoahmpIO%CHXY               = CHXY
    NoahmpIO%FWETXY             = FWETXY
    NoahmpIO%SNEQVOXY           = SNEQVOXY
    NoahmpIO%ALBOLDXY           = ALBOLDXY
    NoahmpIO%QSNOWXY            = QSNOWXY
    NoahmpIO%QRAINXY            = QRAINXY
    NoahmpIO%WSLAKEXY           = WSLAKEXY
    NoahmpIO%ZWTXY              = ZWTXY
    NoahmpIO%WAXY               = WAXY
    NoahmpIO%WTXY               = WTXY
    NoahmpIO%TSNOXY             = TSNOXY
    NoahmpIO%ZSNSOXY            = ZSNSOXY
    NoahmpIO%SNICEXY            = SNICEXY
    NoahmpIO%SNLIQXY            = SNLIQXY
    NoahmpIO%LFMASSXY           = LFMASSXY
    NoahmpIO%RTMASSXY           = RTMASSXY
    NoahmpIO%STMASSXY           = STMASSXY
    NoahmpIO%WOODXY             = WOODXY
    NoahmpIO%STBLCPXY           = STBLCPXY
    NoahmpIO%FASTCPXY           = FASTCPXY
    NoahmpIO%LAI                = XLAIXY
    NoahmpIO%XSAIXY             = XSAIXY
    NoahmpIO%TAUSSXY            = TAUSSXY
    NoahmpIO%SMOISEQ            = SMOISEQ
    NoahmpIO%SMCWTDXY           = SMCWTDXY
    NoahmpIO%DEEPRECHXY         = DEEPRECHXY
    NoahmpIO%RECHXY             = RECHXY
    NoahmpIO%GRAINXY            = GRAINXY
    NoahmpIO%GDDXY              = GDDXY
    NoahmpIO%PGSXY              = PGSXY
    NoahmpIO%QTDRAIN            = QTDRAIN
    NoahmpIO%RS                 = RS
    NoahmpIO%ACC_SSOILXY        = ACC_SSOILXY
    NoahmpIO%ACC_QINSURXY       = ACC_QINSURXY
    NoahmpIO%ACC_QSEVAXY        = ACC_QSEVAXY
    NoahmpIO%ACC_ETRANIXY       = ACC_ETRANIXY
    NoahmpIO%ACC_DWATERXY       = ACC_DWATERXY
    NoahmpIO%ACC_PRCPXY         = ACC_PRCPXY
    NoahmpIO%ACC_ECANXY         = ACC_ECANXY
    NoahmpIO%ACC_ETRANXY        = ACC_ETRANXY
    NoahmpIO%ACC_EDIRXY         = ACC_EDIRXY
    NoahmpIO%ALBSOILDIRXY       = ALBSOILDIRXY
    NoahmpIO%ALBSOILDIFXY       = ALBSOILDIFXY
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       NoahmpIO%FSATXY          = FSATXY
       NoahmpIO%WSURFXY         = WSURFXY
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%SNRDSXY         = SNRDSXY
       NoahmpIO%SNFRXY          = SNFRXY
       NoahmpIO%BCPHIXY         = BCPHIXY
       NoahmpIO%BCPHOXY         = BCPHOXY
       NoahmpIO%OCPHIXY         = OCPHIXY
       NoahmpIO%OCPHOXY         = OCPHOXY
       NoahmpIO%DUST1XY         = DUST1XY
       NoahmpIO%DUST2XY         = DUST2XY
       NoahmpIO%DUST3XY         = DUST3XY
       NoahmpIO%DUST4XY         = DUST4XY
       NoahmpIO%DUST5XY         = DUST5XY
       NoahmpIO%MassConcBCPHIXY = MassConcBCPHIXY
       NoahmpIO%MassConcBCPHOXY = MassConcBCPHOXY
       NoahmpIO%MassConcOCPHIXY = MassConcOCPHIXY
       NoahmpIO%MassConcOCPHOXY = MassConcOCPHOXY
       NoahmpIO%MassConcDUST1XY = MassConcDUST1XY
       NoahmpIO%MassConcDUST2XY = MassConcDUST2XY
       NoahmpIO%MassConcDUST3XY = MassConcDUST3XY
       NoahmpIO%MassConcDUST4XY = MassConcDUST4XY
       NoahmpIO%MassConcDUST5XY = MassConcDUST5XY
    endif
#ifdef WRF_HYDRO
    NoahmpIO%sfcheadrt          = sfcheadrt
    NoahmpIO%INFXSRT            = INFXSRT
    NoahmpIO%soldrain           = soldrain
    NoahmpIO%qtiledrain         = qtiledrain
    NoahmpIO%ZWATBLE2D          = ZWATBLE2D              
#endif

    !--------- Input variable mapping end ---------


    !-------- call main Noah-MP driver ------------
    call NoahmpDriverMain(NoahmpIO)


    !--------- Output variable mapping start ---------

    ! in/out NoahmpIO variables mapped to WRF variables
    TSK          = NoahmpIO%TSK
    HFX          = NoahmpIO%HFX
    QFX          = NoahmpIO%QFX
    LH           = NoahmpIO%LH
    GRDFLX       = NoahmpIO%GRDFLX
    SMSTAV       = NoahmpIO%SMSTAV
    SMSTOT       = NoahmpIO%SMSTOT
    SFCRUNOFF    = NoahmpIO%SFCRUNOFF
    UDRUNOFF     = NoahmpIO%UDRUNOFF
    ALBEDO       = NoahmpIO%ALBEDO
    SNOWC        = NoahmpIO%SNOWC
    SMOIS        = NoahmpIO%SMOIS
    SH2O         = NoahmpIO%SH2O
    TSLB         = NoahmpIO%TSLB
    SNOW         = NoahmpIO%SNOW
    SNOWH        = NoahmpIO%SNOWH
    CANWAT       = NoahmpIO%CANWAT
    CANICEXY     = NoahmpIO%CANICEXY
    CANLIQXY     = NoahmpIO%CANLIQXY
    ACSNOM       = NoahmpIO%ACSNOM
    ACSNOW       = NoahmpIO%ACSNOW
    EMISS        = NoahmpIO%EMISS
    QSFC         = NoahmpIO%QSFC
    Z0           = NoahmpIO%Z0
    ZNT          = NoahmpIO%ZNT
    IRNUMSI      = NoahmpIO%IRNUMSI
    IRNUMMI      = NoahmpIO%IRNUMMI
    IRNUMFI      = NoahmpIO%IRNUMFI
    IRWATSI      = NoahmpIO%IRWATSI
    IRWATMI      = NoahmpIO%IRWATMI
    IRWATFI      = NoahmpIO%IRWATFI
    IRELOSS      = NoahmpIO%IRELOSS
    IRSIVOL      = NoahmpIO%IRSIVOL
    IRMIVOL      = NoahmpIO%IRMIVOL
    IRFIVOL      = NoahmpIO%IRFIVOL
    IRRSPLH      = NoahmpIO%IRRSPLH
    ISNOWXY      = NoahmpIO%ISNOWXY
    TVXY         = NoahmpIO%TVXY
    TGXY         = NoahmpIO%TGXY
    EAHXY        = NoahmpIO%EAHXY
    TAHXY        = NoahmpIO%TAHXY
    CMXY         = NoahmpIO%CMXY
    CHXY         = NoahmpIO%CHXY
    FWETXY       = NoahmpIO%FWETXY
    SNEQVOXY     = NoahmpIO%SNEQVOXY
    ALBOLDXY     = NoahmpIO%ALBOLDXY
    QSNOWXY      = NoahmpIO%QSNOWXY
    QRAINXY      = NoahmpIO%QRAINXY
    WSLAKEXY     = NoahmpIO%WSLAKEXY
    ZWTXY        = NoahmpIO%ZWTXY
    WAXY         = NoahmpIO%WAXY
    WTXY         = NoahmpIO%WTXY
    TSNOXY       = NoahmpIO%TSNOXY
    ZSNSOXY      = NoahmpIO%ZSNSOXY
    SNICEXY      = NoahmpIO%SNICEXY
    SNLIQXY      = NoahmpIO%SNLIQXY
    LFMASSXY     = NoahmpIO%LFMASSXY
    RTMASSXY     = NoahmpIO%RTMASSXY
    STMASSXY     = NoahmpIO%STMASSXY
    WOODXY       = NoahmpIO%WOODXY
    STBLCPXY     = NoahmpIO%STBLCPXY
    FASTCPXY     = NoahmpIO%FASTCPXY
    XLAIXY       = NoahmpIO%LAI
    XSAIXY       = NoahmpIO%XSAIXY
    TAUSSXY      = NoahmpIO%TAUSSXY
    SMOISEQ      = NoahmpIO%SMOISEQ
    SMCWTDXY     = NoahmpIO%SMCWTDXY
    DEEPRECHXY   = NoahmpIO%DEEPRECHXY
    RECHXY       = NoahmpIO%RECHXY
    GRAINXY      = NoahmpIO%GRAINXY
    GDDXY        = NoahmpIO%GDDXY
    PGSXY        = NoahmpIO%PGSXY
    QTDRAIN      = NoahmpIO%QTDRAIN
    RS           = NoahmpIO%RS
    ACC_SSOILXY  = NoahmpIO%ACC_SSOILXY
    ACC_QINSURXY = NoahmpIO%ACC_QINSURXY
    ACC_QSEVAXY  = NoahmpIO%ACC_QSEVAXY
    ACC_ETRANIXY = NoahmpIO%ACC_ETRANIXY
    ACC_DWATERXY = NoahmpIO%ACC_DWATERXY
    ACC_PRCPXY   = NoahmpIO%ACC_PRCPXY
    ACC_ECANXY   = NoahmpIO%ACC_ECANXY
    ACC_ETRANXY  = NoahmpIO%ACC_ETRANXY
    ACC_EDIRXY   = NoahmpIO%ACC_EDIRXY
    ALBSOILDIRXY = NoahmpIO%ALBSOILDIRXY
    ALBSOILDIFXY = NoahmpIO%ALBSOILDIFXY
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       FSATXY    = NoahmpIO%FSATXY
       WSURFXY   = NoahmpIO%WSURFXY
    endif
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       SNRDSXY         = NoahmpIO%SNRDSXY
       SNFRXY          = NoahmpIO%SNFRXY
       BCPHIXY         = NoahmpIO%BCPHIXY
       BCPHOXY         = NoahmpIO%BCPHOXY
       OCPHIXY         = NoahmpIO%OCPHIXY
       OCPHOXY         = NoahmpIO%OCPHOXY
       DUST1XY         = NoahmpIO%DUST1XY
       DUST2XY         = NoahmpIO%DUST2XY
       DUST3XY         = NoahmpIO%DUST3XY
       DUST4XY         = NoahmpIO%DUST4XY
       DUST5XY         = NoahmpIO%DUST5XY
       MassConcBCPHIXY = NoahmpIO%MassConcBCPHIXY
       MassConcBCPHOXY = NoahmpIO%MassConcBCPHOXY
       MassConcOCPHIXY = NoahmpIO%MassConcOCPHIXY
       MassConcOCPHOXY = NoahmpIO%MassConcOCPHOXY
       MassConcDUST1XY = NoahmpIO%MassConcDUST1XY
       MassConcDUST2XY = NoahmpIO%MassConcDUST2XY
       MassConcDUST3XY = NoahmpIO%MassConcDUST3XY
       MassConcDUST4XY = NoahmpIO%MassConcDUST4XY
       MassConcDUST5XY = NoahmpIO%MassConcDUST5XY
    endif
#ifdef WRF_HYDRO
    sfcheadrt    = NoahmpIO%sfcheadrt
    INFXSRT      = NoahmpIO%INFXSRT
    soldrain     = NoahmpIO%soldrain
    qtiledrain   = NoahmpIO%qtiledrain
    ZWATBLE2D    = NoahmpIO%ZWATBLE2D
#endif

    ! output NoahmpIO variables mapped to WRF variables
    T2MVXY       = NoahmpIO%T2MVXY
    T2MBXY       = NoahmpIO%T2MBXY
    Q2MVXY       = NoahmpIO%Q2MVXY
    Q2MBXY       = NoahmpIO%Q2MBXY
    TRADXY       = NoahmpIO%TRADXY
    NEEXY        = NoahmpIO%NEEXY
    GPPXY        = NoahmpIO%GPPXY
    NPPXY        = NoahmpIO%NPPXY
    FVEGXY       = NoahmpIO%FVEGXY
    RUNSFXY      = NoahmpIO%RUNSFXY
    RUNSBXY      = NoahmpIO%RUNSBXY
    ECANXY       = NoahmpIO%ECANXY
    EDIRXY       = NoahmpIO%EDIRXY
    ETRANXY      = NoahmpIO%ETRANXY
    FSAXY        = NoahmpIO%FSAXY
    FIRAXY       = NoahmpIO%FIRAXY
    APARXY       = NoahmpIO%APARXY
    PSNXY        = NoahmpIO%PSNXY
    SAVXY        = NoahmpIO%SAVXY
    SAGXY        = NoahmpIO%SAGXY
    RSSUNXY      = NoahmpIO%RSSUNXY
    RSSHAXY      = NoahmpIO%RSSHAXY
    BGAPXY       = NoahmpIO%BGAPXY
    WGAPXY       = NoahmpIO%WGAPXY
    TGVXY        = NoahmpIO%TGVXY
    TGBXY        = NoahmpIO%TGBXY
    CHVXY        = NoahmpIO%CHVXY
    CHBXY        = NoahmpIO%CHBXY
    SHGXY        = NoahmpIO%SHGXY
    SHCXY        = NoahmpIO%SHCXY
    SHBXY        = NoahmpIO%SHBXY
    EVGXY        = NoahmpIO%EVGXY
    EVBXY        = NoahmpIO%EVBXY
    GHVXY        = NoahmpIO%GHVXY
    GHBXY        = NoahmpIO%GHBXY
    IRGXY        = NoahmpIO%IRGXY 
    IRCXY        = NoahmpIO%IRCXY
    IRBXY        = NoahmpIO%IRBXY
    TRXY         = NoahmpIO%TRXY
    EVCXY        = NoahmpIO%EVCXY
    CHLEAFXY     = NoahmpIO%CHLEAFXY
    CHUCXY       = NoahmpIO%CHUCXY
    CHV2XY       = NoahmpIO%CHV2XY
    CHB2XY       = NoahmpIO%CHB2XY
    QINTSXY      = NoahmpIO%QINTSXY
    QINTRXY      = NoahmpIO%QINTRXY
    QDRIPSXY     = NoahmpIO%QDRIPSXY
    QDRIPRXY     = NoahmpIO%QDRIPRXY
    QTHROSXY     = NoahmpIO%QTHROSXY
    QTHRORXY     = NoahmpIO%QTHRORXY
    QSNSUBXY     = NoahmpIO%QSNSUBXY
    QSNFROXY     = NoahmpIO%QSNFROXY
    QSUBCXY      = NoahmpIO%QSUBCXY
    QFROCXY      = NoahmpIO%QFROCXY
    QEVACXY      = NoahmpIO%QEVACXY
    QDEWCXY      = NoahmpIO%QDEWCXY
    QFRZCXY      = NoahmpIO%QFRZCXY
    QMELTCXY     = NoahmpIO%QMELTCXY
    QSNBOTXY     = NoahmpIO%QSNBOTXY
    QMELTXY      = NoahmpIO%QMELTXY
    PONDINGXY    = NoahmpIO%PONDINGXY
    PAHXY        = NoahmpIO%PAHXY
    PAHGXY       = NoahmpIO%PAHGXY
    PAHVXY       = NoahmpIO%PAHVXY
    PAHBXY       = NoahmpIO%PAHBXY
    FPICEXY      = NoahmpIO%FPICEXY
    RAINLSM      = NoahmpIO%RAINLSM
    SNOWLSM      = NoahmpIO%SNOWLSM
    FORCTLSM     = NoahmpIO%FORCTLSM
    FORCQLSM     = NoahmpIO%FORCQLSM
    FORCPLSM     = NoahmpIO%FORCPLSM
    FORCZLSM     = NoahmpIO%FORCZLSM
    FORCWLSM     = NoahmpIO%FORCWLSM
    EFLXBXY      = NoahmpIO%EFLXBXY
    SOILENERGY   = NoahmpIO%SOILENERGY
    SNOWENERGY   = NoahmpIO%SNOWENERGY
    CANHSXY      = NoahmpIO%CANHSXY

    !--------- Output variable mapping end ---------

  end subroutine NoahmpWRFmain

end module NoahmpWRFmainMod
