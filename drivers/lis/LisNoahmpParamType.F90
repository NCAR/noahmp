module LisNoahmpParamType

!!! Define LIS-Noah-MP table parameter variables

! ------------------------ Code history -----------------------------------
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine

  implicit none
  save
  private

  integer, private, parameter :: MBAND  = 2
  integer, private, parameter :: NSOIL  = 4
  integer, private, parameter :: NSTAGE = 8

  type, public :: LisNoahmpParam_type

!----------------------------------------------------------------
! Noahmp Parameters Table 
!----------------------------------------------------------------

    ! vegetation parameters
    logical                     :: URBAN_FLAG          ! urban flag
    integer                     :: ISWATER             ! water flag
    integer                     :: ISBARREN            ! barren ground flag
    integer                     :: ISICE               ! ice flag
    integer                     :: ISCROP              ! cropland flag
    integer                     :: EBLFOREST           ! evergreen broadleaf forest flag
    real(kind=kind_noahmp)      :: CH2OP               ! maximum intercepted h2o per unit lai+sai (mm)
    real(kind=kind_noahmp)      :: DLEAF               ! characteristic leaf dimension (m)
    real(kind=kind_noahmp)      :: Z0MVT               ! momentum roughness length (m)
    real(kind=kind_noahmp)      :: HVT                 ! top of canopy (m)
    real(kind=kind_noahmp)      :: HVB                 ! bottom of canopy (m)
    real(kind=kind_noahmp)      :: DEN                 ! tree density (no. of trunks per m2)
    real(kind=kind_noahmp)      :: RC                  ! tree crown radius (m)
    real(kind=kind_noahmp)      :: MFSNO               ! snowmelt curve parameter
    real(kind=kind_noahmp)      :: SCFFAC              ! snow cover factor (m) (replace original hard-coded 2.5*z0 in SCF formulation)
    real(kind=kind_noahmp)      :: CBIOM               ! canopy biomass heat capacity parameter (m) 
    real(kind=kind_noahmp)      :: SAIM(12)            ! monthly stem area index, one-sided
    real(kind=kind_noahmp)      :: LAIM(12)            ! monthly leaf area index, one-sided
    real(kind=kind_noahmp)      :: SLA                 ! single-side leaf area per Kg [m2/kg]
    real(kind=kind_noahmp)      :: DILEFC              ! coeficient for leaf stress death [1/s]
    real(kind=kind_noahmp)      :: DILEFW              ! coeficient for leaf stress death [1/s]
    real(kind=kind_noahmp)      :: FRAGR               ! fraction of growth respiration  !original was 0.3 
    real(kind=kind_noahmp)      :: LTOVRC              ! leaf turnover [1/s]
    real(kind=kind_noahmp)      :: C3PSN               ! photosynthetic pathway: 0. = c4, 1. = c3
    real(kind=kind_noahmp)      :: KC25                ! co2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp)      :: AKC                 ! q10 for kc25
    real(kind=kind_noahmp)      :: KO25                ! o2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp)      :: AKO                 ! q10 for ko25
    real(kind=kind_noahmp)      :: VCMX25              ! maximum rate of carboxylation at 25c (umol co2/m2/s)
    real(kind=kind_noahmp)      :: AVCMX               ! q10 for vcmx25
    real(kind=kind_noahmp)      :: BP                  ! minimum leaf conductance (umol/m2/s)
    real(kind=kind_noahmp)      :: MP                  ! slope of conductance-to-photosynthesis relationship
    real(kind=kind_noahmp)      :: QE25                ! quantum efficiency at 25c (umol co2 / umol photon)
    real(kind=kind_noahmp)      :: AQE                 ! q10 for qe25
    real(kind=kind_noahmp)      :: RMF25               ! leaf maintenance respiration at 25c (umol co2/m2/s)
    real(kind=kind_noahmp)      :: RMS25               ! stem maintenance respiration at 25c (umol co2/kg bio/s)
    real(kind=kind_noahmp)      :: RMR25               ! root maintenance respiration at 25c (umol co2/kg bio/s)
    real(kind=kind_noahmp)      :: ARM                 ! q10 for maintenance respiration
    real(kind=kind_noahmp)      :: FOLNMX              ! foliage nitrogen concentration when f(n)=1 (%)
    real(kind=kind_noahmp)      :: TMIN                ! minimum temperature for photosynthesis (k)
    real(kind=kind_noahmp)      :: XL                  ! leaf/stem orientation index
    real(kind=kind_noahmp)      :: RHOL(MBAND)         ! leaf reflectance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: RHOS(MBAND)         ! stem reflectance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: TAUL(MBAND)         ! leaf transmittance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: TAUS(MBAND)         ! stem transmittance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: MRP                 ! microbial respiration parameter (umol co2 /kg c/ s)
    real(kind=kind_noahmp)      :: CWPVT               ! empirical canopy wind parameter
    real(kind=kind_noahmp), allocatable, dimension(:)      :: WRRAT_TABLE               ! wood to non-wood ratio
    real(kind=kind_noahmp), allocatable, dimension(:)      :: WDPOOL_TABLE              ! wood pool (switch 1 or 0) depending on woody or not [-]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TDLEF_TABLE               ! characteristic T for leaf freezing [K]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: NROOT_TABLE               ! number of soil layers with root present
    real(kind=kind_noahmp), allocatable, dimension(:)      :: RGL_TABLE                 ! Parameter used in radiation stress function
    real(kind=kind_noahmp), allocatable, dimension(:)      :: RS_TABLE                  ! Minimum stomatal resistance [s m-1]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: HS_TABLE                  ! Parameter used in vapor pressure deficit function
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TOPT_TABLE                ! Optimum transpiration air temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: RSMAX_TABLE               ! Maximal stomatal resistance [s m-1]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: RTOVRC_TABLE              ! root turnover coefficient [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: RSWOODC_TABLE             ! wood respiration coeficient [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BF_TABLE                  ! parameter for present wood allocation [-]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: WSTRC_TABLE               ! water stress coeficient [-]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: LAIMIN_TABLE              ! minimum leaf area index [m2/m2]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: XSAMIN_TABLE              ! minimum stem area index [m2/m2]

    ! radiation parameters
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: ALBSAT_TABLE              ! saturated soil albedos: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: ALBDRY_TABLE              ! dry soil albedos: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:)      :: ALBICE_TABLE              ! albedo land ice: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:)      :: ALBLAK_TABLE              ! albedo frozen lakes: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:)      :: OMEGAS_TABLE              ! two-stream parameter omega for snow
    real(kind=kind_noahmp)                                 :: BETADS_TABLE              ! two-stream parameter betad for snow
    real(kind=kind_noahmp)                                 :: BETAIS_TABLE              ! two-stream parameter betad for snow
    real(kind=kind_noahmp), allocatable, dimension(:)      :: EG_TABLE                  ! emissivity soil surface
    real(kind=kind_noahmp)                                 :: EICE_TABLE                ! ice surface emissivity

    ! global parameters
    real(kind=kind_noahmp)                                 :: CO2_TABLE                 ! co2 partial pressure
    real(kind=kind_noahmp)                                 :: O2_TABLE                  ! o2 partial pressure
    real(kind=kind_noahmp)                                 :: TIMEAN_TABLE              ! gridcell mean topgraphic index (global mean)
    real(kind=kind_noahmp)                                 :: FSATMX_TABLE              ! maximum surface saturated fraction (global mean)
    real(kind=kind_noahmp)                                 :: Z0SNO_TABLE               ! snow surface roughness length (m) (0.002)
    real(kind=kind_noahmp)                                 :: SSI_TABLE                 ! liquid water holding capacity for snowpack (m3/m3) (0.03)
    real(kind=kind_noahmp)                                 :: SNOW_RET_FAC_TABLE        ! snowpack water release timescale factor (1/s)
    real(kind=kind_noahmp)                                 :: SNOW_EMIS_TABLE           ! snow emissivity
    real(kind=kind_noahmp)                                 :: SWEMX_TABLE               ! new snow mass to fully cover old snow (mm)
    real(kind=kind_noahmp)                                 :: TAU0_TABLE                ! tau0 from Yang97 eqn. 10a
    real(kind=kind_noahmp)                                 :: GRAIN_GROWTH_TABLE        ! growth from vapor diffusion Yang97 eqn. 10b
    real(kind=kind_noahmp)                                 :: EXTRA_GROWTH_TABLE        ! extra growth near freezing Yang97 eqn. 10c
    real(kind=kind_noahmp)                                 :: DIRT_SOOT_TABLE           ! dirt and soot term Yang97 eqn. 10d
    real(kind=kind_noahmp)                                 :: BATS_COSZ_TABLE           ! zenith angle snow albedo adjustment; b in Yang97 eqn. 15
    real(kind=kind_noahmp)                                 :: BATS_VIS_NEW_TABLE        ! new snow visible albedo
    real(kind=kind_noahmp)                                 :: BATS_NIR_NEW_TABLE        ! new snow NIR albedo
    real(kind=kind_noahmp)                                 :: BATS_VIS_AGE_TABLE        ! age factor for diffuse visible snow albedo Yang97 eqn. 17
    real(kind=kind_noahmp)                                 :: BATS_NIR_AGE_TABLE        ! age factor for diffuse NIR snow albedo Yang97 eqn. 18
    real(kind=kind_noahmp)                                 :: BATS_VIS_DIR_TABLE        ! cosz factor for direct visible snow albedo Yang97 eqn. 15
    real(kind=kind_noahmp)                                 :: BATS_NIR_DIR_TABLE        ! cosz factor for direct NIR snow albedo Yang97 eqn. 16
    real(kind=kind_noahmp)                                 :: RSURF_SNOW_TABLE          ! surface resistance for snow(s/m)
    real(kind=kind_noahmp)                                 :: RSURF_EXP_TABLE           ! exponent in the shape parameter for soil resistance option 1
    real(kind=kind_noahmp)                                 :: C2_SNOWCOMPACT_TABLE      ! overburden snow compaction parameter (m3/kg)
    real(kind=kind_noahmp)                                 :: C3_SNOWCOMPACT_TABLE      ! snow desctructive metamorphism compaction parameter1 [1/s]
    real(kind=kind_noahmp)                                 :: C4_SNOWCOMPACT_TABLE      ! snow desctructive metamorphism compaction parameter2 [1/k]
    real(kind=kind_noahmp)                                 :: C5_SNOWCOMPACT_TABLE      ! snow desctructive metamorphism compaction parameter3
    real(kind=kind_noahmp)                                 :: DM_SNOWCOMPACT_TABLE      ! upper Limit on destructive metamorphism compaction [kg/m3]
    real(kind=kind_noahmp)                                 :: ETA0_SNOWCOMPACT_TABLE    ! snow viscosity coefficient [kg-s/m2]
    real(kind=kind_noahmp)                                 :: SNLIQMAXFRAC_TABLE        ! maximum liquid water fraction in snow
    real(kind=kind_noahmp)                                 :: SWEMAXGLA_TABLE           ! Maximum SWE allowed at glaciers (mm)
    real(kind=kind_noahmp)                                 :: WSLMAX_TABLE              ! maximum lake water storage (mm)
    real(kind=kind_noahmp)                                 :: ROUS_TABLE                ! specific yield [-] for Niu et al. 2007 groundwater scheme
    real(kind=kind_noahmp)                                 :: CMIC_TABLE                ! microprore content (0.0-1.0), 0.0: close to free drainage
    real(kind=kind_noahmp)                                 :: SNOWDEN_MAX_TABLE         ! maximum fresh snowfall density (kg/m3)
    real(kind=kind_noahmp)                                 :: CLASS_ALB_REF_TABLE       ! reference snow albedo in CLASS scheme
    real(kind=kind_noahmp)                                 :: CLASS_SNO_AGE_TABLE       ! snow aging e-folding time (s) in CLASS albedo scheme
    real(kind=kind_noahmp)                                 :: CLASS_ALB_NEW_TABLE       ! fresh snow albedo in CLASS scheme
    real(kind=kind_noahmp)                                 :: PSIWLT_TABLE              ! soil metric potential for wilting point (m)
    real(kind=kind_noahmp)                                 :: Z0SOIL_TABLE              ! Bare-soil roughness length (m) (i.e., under the canopy)
    real(kind=kind_noahmp)                                 :: Z0LAKE_TABLE              ! Lake surface roughness length (m)

    ! irrigation parameters
    integer                                                :: IRR_HAR_TABLE             ! number of days before harvest date to stop irrigation 
    real(kind=kind_noahmp)                                 :: IRR_FRAC_TABLE            ! irrigation Fraction
    real(kind=kind_noahmp)                                 :: IRR_LAI_TABLE             ! Minimum lai to trigger irrigation
    real(kind=kind_noahmp)                                 :: IRR_MAD_TABLE             ! management allowable deficit (0-1)
    real(kind=kind_noahmp)                                 :: FILOSS_TABLE              ! factor of flood irrigation loss
    real(kind=kind_noahmp)                                 :: SPRIR_RATE_TABLE          ! mm/h, sprinkler irrigation rate
    real(kind=kind_noahmp)                                 :: MICIR_RATE_TABLE          ! mm/h, micro irrigation rate
    real(kind=kind_noahmp)                                 :: FIRTFAC_TABLE             ! flood application rate factor
    real(kind=kind_noahmp)                                 :: IR_RAIN_TABLE             ! maximum precipitation to stop irrigation trigger

    ! tile drainage parameters
    integer                                                :: DRAIN_LAYER_OPT_TABLE     ! tile drainage layer
    integer               , allocatable, dimension(:)      :: TD_DEPTH_TABLE            ! tile drainage depth (layer number) from soil surface
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TDSMC_FAC_TABLE           ! tile drainage soil moisture factor
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TD_DC_TABLE               ! tile drainage coefficient [mm/d]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TD_DCOEF_TABLE            ! tile drainage coefficient [mm/d]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TD_D_TABLE                ! depth to impervious layer from drain water level [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TD_ADEPTH_TABLE           ! actual depth of impervious layer from land surface [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TD_RADI_TABLE             ! effective radius of drain tubes [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TD_SPAC_TABLE             ! distance between two drain tubes or tiles [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TD_DDRAIN_TABLE           ! tile drainage depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: KLAT_FAC_TABLE            ! hydraulic conductivity mutiplification factor

    ! crop parameters
    integer                                                :: DEFAULT_CROP_TABLE        ! Default crop index
    integer               , allocatable, dimension(:)      :: PLTDAY_TABLE              ! Planting date
    integer               , allocatable, dimension(:)      :: HSDAY_TABLE               ! Harvest date
    real(kind=kind_noahmp), allocatable, dimension(:)      :: PLANTPOP_TABLE            ! Plant density [per ha] - used?
    real(kind=kind_noahmp), allocatable, dimension(:)      :: IRRI_TABLE                ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDTBASE_TABLE            ! Base temperature for GDD accumulation [C]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDTCUT_TABLE             ! Upper temperature for GDD accumulation [C]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS1_TABLE               ! GDD from seeding to emergence
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS2_TABLE               ! GDD from seeding to initial vegetative 
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS3_TABLE               ! GDD from seeding to post vegetative 
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS4_TABLE               ! GDD from seeding to intial reproductive
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS5_TABLE               ! GDD from seeding to pysical maturity 
    real(kind=kind_noahmp), allocatable, dimension(:)      :: C3PSNI_TABLE              ! photosynthetic pathway: 0. = c4, 1. = c3 ! Zhe Zhang 2020-07-03
    real(kind=kind_noahmp), allocatable, dimension(:)      :: KC25I_TABLE               ! co2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: AKCI_TABLE                ! q10 for kc25
    real(kind=kind_noahmp), allocatable, dimension(:)      :: KO25I_TABLE               ! o2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: AKOI_TABLE                ! q10 for ko25
    real(kind=kind_noahmp), allocatable, dimension(:)      :: VCMX25I_TABLE             ! maximum rate of carboxylation at 25c (umol co2/m2/s)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: AVCMXI_TABLE              ! q10 for vcmx25
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BPI_TABLE                 ! minimum leaf conductance (umol/m2/s)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: MPI_TABLE                 ! slope of conductance-to-photosynthesis relationship
    real(kind=kind_noahmp), allocatable, dimension(:)      :: QE25I_TABLE               ! quantum efficiency at 25c (umol co2 / umol photon)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: FOLNMXI_TABLE             ! foliage nitrogen concentration when f(n)=1 (%)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: AREF_TABLE                ! reference maximum CO2 assimulation rate 
    real(kind=kind_noahmp), allocatable, dimension(:)      :: PSNRF_TABLE               ! CO2 assimulation reduction factor(0-1) (caused by non-modeled part, pest,weeds)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: I2PAR_TABLE               ! Fraction of incoming solar radiation to photosynthetically active radiation
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TASSIM0_TABLE             ! Minimum temperature for CO2 assimulation [C]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TASSIM1_TABLE             ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: TASSIM2_TABLE             ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: K_TABLE                   ! light extinction coefficient
    real(kind=kind_noahmp), allocatable, dimension(:)      :: EPSI_TABLE                ! initial light use efficiency
    real(kind=kind_noahmp), allocatable, dimension(:)      :: Q10MR_TABLE               ! q10 for maintainance respiration
    real(kind=kind_noahmp), allocatable, dimension(:)      :: LEFREEZ_TABLE             ! characteristic T for leaf freezing [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: DILE_FC_TABLE             ! coeficient for temperature leaf stress death [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: DILE_FW_TABLE             ! coeficient for water leaf stress death [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: FRA_GR_TABLE              ! fraction of growth respiration
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: LF_OVRC_TABLE             ! fraction of leaf turnover  [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: ST_OVRC_TABLE             ! fraction of stem turnover  [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: RT_OVRC_TABLE             ! fraction of root tunrover  [1/s]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: LFMR25_TABLE              ! leaf maintenance respiration at 25C [umol CO2/m2/s]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: STMR25_TABLE              ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: RTMR25_TABLE              ! root maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GRAINMR25_TABLE           ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: LFPT_TABLE                ! fraction of carbohydrate flux to leaf
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: STPT_TABLE                ! fraction of carbohydrate flux to stem
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: RTPT_TABLE                ! fraction of carbohydrate flux to root
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: GRAINPT_TABLE             ! fraction of carbohydrate flux to grain
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: LFCT_TABLE                ! fraction of carbohydrate translocation from leaf to grain 
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: STCT_TABLE                ! fraction of carbohydrate translocation from stem to grain
    real(kind=kind_noahmp), allocatable, dimension(:,:)    :: RTCT_TABLE                ! fraction of carbohydrate translocation from root to grain
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BIO2LAI_TABLE             ! leaf area per living leaf biomass [m2/kg]

    ! soil parameters
    integer                                                :: SLCATS_TABLE              ! number of soil categories
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BEXP_TABLE                ! soil B parameter
    real(kind=kind_noahmp), allocatable, dimension(:)      :: SMCDRY_TABLE              ! dry soil moisture threshold
    real(kind=kind_noahmp), allocatable, dimension(:)      :: SMCMAX_TABLE              ! porosity, saturated value of soil moisture (volumetric)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: SMCREF_TABLE              ! reference soil moisture (field capacity) (volumetric)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: PSISAT_TABLE              ! saturated soil matric potential
    real(kind=kind_noahmp), allocatable, dimension(:)      :: DKSAT_TABLE               ! saturated soil hydraulic conductivity
    real(kind=kind_noahmp), allocatable, dimension(:)      :: DWSAT_TABLE               ! saturated soil hydraulic diffusivity
    real(kind=kind_noahmp), allocatable, dimension(:)      :: SMCWLT_TABLE              ! wilting point soil moisture (volumetric)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: QUARTZ_TABLE              ! soil quartz content
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BVIC_TABLE                ! VIC model infiltration parameter (-) for opt_run=6
    real(kind=kind_noahmp), allocatable, dimension(:)      :: AXAJ_TABLE                ! Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BXAJ_TABLE                ! Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)      :: XXAJ_TABLE                ! Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BDVIC_TABLE               ! VIC model infiltration parameter (-)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: GDVIC_TABLE               ! mean capilary drive (m)
    real(kind=kind_noahmp), allocatable, dimension(:)      :: BBVIC_TABLE               ! heterogeniety parameter for DVIC infiltration [-]

    ! general parameters
    real(kind=kind_noahmp), allocatable, dimension(:)      :: SLOPE_TABLE               ! slope factor for soil drainage
    real(kind=kind_noahmp)                                 :: CSOIL_TABLE               ! Soil heat capacity [J m-3 K-1]
    real(kind=kind_noahmp)                                 :: REFDK_TABLE               ! Parameter in the surface runoff parameterization
    real(kind=kind_noahmp)                                 :: REFKDT_TABLE              ! Parameter in the surface runoff parameterization
    real(kind=kind_noahmp)                                 :: FRZK_TABLE                ! Frozen ground parameter
    real(kind=kind_noahmp)                                 :: ZBOT_TABLE                ! Depth [m] of lower boundary soil temperature
    real(kind=kind_noahmp)                                 :: CZIL_TABLE                ! Parameter used in the calculation of the roughness length for heat

    ! optional parameters
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_a_TABLE      ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_b_TABLE      ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_c_TABLE      ! orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_d_TABLE      ! sand*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_e_TABLE      ! clay*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_f_TABLE      ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_g_TABLE      ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500_a_TABLE       ! theta_1500t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500_b_TABLE       ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_a_TABLE        ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_b_TABLE        ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_c_TABLE        ! orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_d_TABLE        ! sand*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_e_TABLE        ! clay*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_f_TABLE        ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_g_TABLE        ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_33_a_TABLE         ! theta_33t*theta_33t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33_b_TABLE         ! theta_33t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33_c_TABLE         ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_a_TABLE       ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_b_TABLE       ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_c_TABLE       ! orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_d_TABLE       ! sand*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_e_TABLE       ! clay*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_f_TABLE       ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_g_TABLE       ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33_a_TABLE        ! theta_s33t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33_b_TABLE        ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_a_TABLE           ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_b_TABLE           ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_c_TABLE           ! theta_s33 coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_d_TABLE           ! sand*theta_s33 coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_e_TABLE           ! clay*theta_s33 coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_f_TABLE           ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_g_TABLE           ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_psi_e_a_TABLE            ! psi_et*psi_et coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_e_b_TABLE            ! psi_et coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_e_c_TABLE            ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_smcmax_a_TABLE           ! sand adjustment
    real(kind=kind_noahmp)                                 :: sr2006_smcmax_b_TABLE           ! constant adjustment

  end type LisNoahmpParam_type

end module LisNoahmpParamType
