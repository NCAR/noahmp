module lnd_comp_types 

  use ESMF           , only : ESMF_Grid, ESMF_Mesh, ESMF_RouteHandle
  use machine        , only : kp => kind_phys
  use mpp_domains_mod, only : domain2d

  !----------------------------------------------------------------------------
  ! Land component specific data types
  !----------------------------------------------------------------------------
  public

  ! data type for initial conditions
  type initial_type
     real(kind=kp), allocatable :: snow_water_equivalent(:)
     real(kind=kp), allocatable :: snow_depth(:)
     real(kind=kp), allocatable :: canopy_water(:)
     real(kind=kp), allocatable :: skin_temperature(:)
     real(kind=kp), allocatable :: soil_temperature(:,:)
     real(kind=kp), allocatable :: soil_moisture(:,:)
     real(kind=kp), allocatable :: soil_liquid(:,:)
     real(kind=kp), allocatable :: surface_roughness(:)
     real(kind=kp), allocatable :: friction_velocity(:)
  end type initial_type

  ! data type for forcing
  type forcing_type
    real(kind=kp),  allocatable :: t1         (:) ! air temperature (K)
    real(kind=kp),  allocatable :: q1         (:) ! mixing ratio/specific humidty at lowest model layer (kg/kg)
    real(kind=kp),  allocatable :: u1         (:) ! u-component of wind (m/s)
    real(kind=kp),  allocatable :: v1         (:) ! v-component of wind (m/s)
    real(kind=kp),  allocatable :: ps         (:) ! surface pressure (Pa)
    real(kind=kp),  allocatable :: pbot       (:) ! bottom layer pressure (Pa)
    real(kind=kp),  allocatable :: tskin      (:) ! skin temperature (K)
    real(kind=kp),  allocatable :: dlwflx     (:) ! downward longwave radiation (W/m2)
    real(kind=kp),  allocatable :: dswsfc     (:) ! downward shortwave radiation (W/m2)
    real(kind=kp),  allocatable :: snet       (:) ! net shortwave radiation (W/m2)
    real(kind=kp),  allocatable :: wind       (:) ! wind speed (m/s)
    real(kind=kp),  allocatable :: tprcp      (:) ! total precipitation (mm/s)
    real(kind=kp),  allocatable :: tprcpc     (:) ! convective component of precipitation (mm/s)
    real(kind=kp),  allocatable :: tprcpl     (:) ! large-scale component of precipitation (mm/s)
    real(kind=kp),  allocatable :: snow       (:) ! snow fall (mm/s)
    real(kind=kp),  allocatable :: snowc      (:) ! convective component of snow fall (mm/s)
    real(kind=kp),  allocatable :: snowl      (:) ! large-scale component of snow fall (mm/s)
    real(kind=kp),  allocatable :: vegfrac    (:) ! vegetation fraction (unitless, 0-1)
    real(kind=kp),  allocatable :: hgt        (:) ! forcing height (m)
    real(kind=kp),  allocatable :: prslk1     (:) ! dimensionless Exner function at the lowest model layer
    real(kind=kp),  allocatable :: ustar1     (:) ! friction velocity (m/s)
    real(kind=kp),  allocatable :: zorl       (:) ! surface roughness (m)
  end type forcing_type

  ! data type for static information provided by nems.configure
  type static_type
    integer                     :: im             ! horiz dimension and number of used points
    integer                     :: km             ! vertical soil layer dimension
    integer                     :: lsnowl = -2    ! lower bound of vertical dimension of surface snow
    integer                     :: itime = huge(1)! not used, set to something huge by default
    integer                     :: isot           ! sfc soil type data source zobler or statsgo
    integer                     :: ivegsrc        ! sfc veg type data source umd or igbp
    integer                     :: idveg          ! option for dynamic vegetation
    real(kind=kp)               :: delt           ! time interval (s) 
    integer                     :: iopt_crs       ! option for canopy stomatal resistance
    integer                     :: iopt_btr       ! option for soil moisture factor for stomatal resistance
    integer                     :: iopt_run       ! option for runoff and groundwater
    integer                     :: iopt_sfc       ! option for surface layer drag coeff (ch & cm)
    integer                     :: iopt_frz       ! option for supercooled liquid water (or ice fraction)
    integer                     :: iopt_inf       ! option for frozen soil permeability
    integer                     :: iopt_rad       ! option for radiation transfer
    integer                     :: iopt_alb       ! option for ground snow surface albedo
    integer                     :: iopt_snf       ! option for partitioning  precipitation into rainfall & snowfall
    integer                     :: iopt_tbot      ! option for lower boundary condition of soil temperature
    integer                     :: iopt_stc       ! option for snow/soil temperature time scheme (only layer 1)
    integer                     :: iopt_rsf       ! option for surface resistent to evaporation/sublimation
    integer                     :: iopt_gla       ! option for glacier treatment
    integer                     :: iopt_trs       ! option for surface thermal roughness option
    logical                     :: do_mynnedmf    ! option for MYNN-EDMF
    logical                     :: do_mynnsfclay  ! option for MYNN surface layer scheme
    character(len=128)          :: errmsg         ! error message
    integer                     :: errflg         ! error flag 
  end type static_type

  ! data type for model internal data structures
  type model_type
     ! scalar variables
     integer                    :: iyrlen           ! year length
     real(kind=kp)              :: julian           ! julian day of year
     character(len=19)          :: reference_date   ! reference date
     logical                    :: thsfc_loc        ! flag for reference pressure theta
     ! variables in dimension im
     real(kind=kp), allocatable :: u1         (:)   ! u-component of wind (m/s)
     real(kind=kp), allocatable :: v1         (:)   ! v-component of wind (m/s)
     integer      , allocatable :: soiltyp    (:)   ! soil type (integer index)
     integer      , allocatable :: vegtype    (:)   ! vegetation type (integer index)
     real(kind=kp), allocatable :: sigmaf     (:)   ! areal fractional cover of green vegetation
     real(kind=kp), allocatable :: emiss      (:)   ! sfc lw emissivity (fraction)
     real(kind=kp), allocatable :: albdvis    (:)   ! albedo - direct  visible (fraction)
     real(kind=kp), allocatable :: albdnir    (:)   ! albedo - direct  NIR (fraction)
     real(kind=kp), allocatable :: albivis    (:)   ! albedo - diffuse visible (fraction)
     real(kind=kp), allocatable :: albinir    (:)   ! albedo - diffuse NIR (fraction)
     real(kind=kp), allocatable :: snet       (:)   ! total sky sfc netsw flx into ground (W/m^2)    NOT USED
     real(kind=kp), allocatable :: tg3        (:)   ! deep soil temperature (K)
     real(kind=kp), allocatable :: cm         (:)   ! surface exchange coeff for momentum (m/s)
     real(kind=kp), allocatable :: ch         (:)   ! surface exchange coeff heat & moisture (m/s)
     real(kind=kp), allocatable :: prsl1      (:)   ! sfc layer 1 mean pressure (Pa)
     real(kind=kp), allocatable :: prslki     (:)   ! ? 
     real(kind=kp), allocatable :: prslk1     (:)   ! dimensionless exner function at surface adjacent layer
     real(kind=kp), allocatable :: prsik1     (:)   ! surface dimensionless exner function
     real(kind=kp), allocatable :: zf         (:)   ! height of bottom layer (m)
     real(kind=kp), allocatable :: pblh       (:)   ! PBL thickness (m)
     logical      , allocatable :: dry        (:)   ! = T if a point with any land
     integer      , allocatable :: slopetyp   (:)   ! class of sfc slope (integer index)
     real(kind=kp), allocatable :: alb_monthly(:,:) ! surface albedo
     real(kind=kp), allocatable :: gvf_monthly(:,:) ! fractional coverage of green veg
     real(kind=kp), allocatable :: shdmin     (:)   ! min fractional coverage of green veg           NOT USED
     real(kind=kp), allocatable :: shdmax     (:)   ! max fractnl cover of green veg                 NOT USED
     real(kind=kp), allocatable :: snoalb     (:)   ! upper bound on max albedo over deep snow       NOT USED
     real(kind=kp), allocatable :: sfalb      (:)   ! mean sfc diffused sw albedo (fractional)       NOT USED
     logical      , allocatable :: flag_iter  (:)   ! ? 
     real(kind=kp), allocatable :: xlatin     (:)   ! latitude
     real(kind=kp), allocatable :: xcoszin    (:)   ! cosine of zenith angle
     real(kind=kp), allocatable :: rainn_mp   (:)   ! microphysics non-convective precipitation (mm)
     real(kind=kp), allocatable :: rainc_mp   (:)   ! microphysics convective precipitation (mm)
     real(kind=kp), allocatable :: snow_mp    (:)   ! microphysics snow (mm)
     real(kind=kp), allocatable :: graupel_mp (:)   ! microphysics graupel (mm)
     real(kind=kp), allocatable :: ice_mp     (:)   ! microphysics ice/hail (mm)
     real(kind=kp), allocatable :: tprcp      (:)   ! total precipitation (mm/s)
     real(kind=kp), allocatable :: weasd      (:)   ! water equivalent accumulated snow depth (mm)
     real(kind=kp), allocatable :: snwdph     (:)   ! snow depth (water equiv) over land
     real(kind=kp), allocatable :: tskin      (:)   ! ground surface skin temperature (K)
     real(kind=kp), allocatable :: srflag     (:)   ! snow/rain flag for precipitation
     real(kind=kp), allocatable :: canopy     (:)   ! canopy moisture content (m)
     real(kind=kp), allocatable :: trans      (:)   ! total plant transpiration (m/s)
     real(kind=kp), allocatable :: tsurf      (:)   ! surface skin temperature (after iteration)
     real(kind=kp), allocatable :: zorl       (:)   ! surface roughness
     real(kind=kp), allocatable :: rb1        (:)   ! composite bulk richardson number
     real(kind=kp), allocatable :: fm1        (:)   ! composite momemtum stability
     real(kind=kp), allocatable :: fh1        (:)   ! composite heat/moisture stability
     real(kind=kp), allocatable :: ustar1     (:)   ! composite friction velocity
     real(kind=kp), allocatable :: stress1    (:)   ! composite surface stress
     real(kind=kp), allocatable :: fm101      (:)   ! composite 2-meter momemtum stability
     real(kind=kp), allocatable :: fh21       (:)   ! composite 10-meter heat/moisture stability
     real(kind=kp), allocatable :: rmol1      (:)   ! one over obukhov length
     real(kind=kp), allocatable :: flhc1      (:)   ! surface exchange coefficient for heat
     real(kind=kp), allocatable :: flqc1      (:)   ! surface exchange coefficient for moisture
     logical                    :: do_mynnsfclay    ! flag to activate MYNN surface layer
     real(kind=kp), allocatable :: snowxy     (:)   ! actual no. of snow layers
     real(kind=kp), allocatable :: tvxy       (:)   ! vegetation leaf temperature (K)
     real(kind=kp), allocatable :: tgxy       (:)   ! bulk ground surface temperature (K)
     real(kind=kp), allocatable :: canicexy   (:)   ! canopy-intercepted ice (mm)
     real(kind=kp), allocatable :: canliqxy   (:)   ! canopy-intercepted liquid water (mm)
     real(kind=kp), allocatable :: eahxy      (:)   ! canopy air vapor pressure (Pa)
     real(kind=kp), allocatable :: tahxy      (:)   ! canopy air temperature (K)
     real(kind=kp), allocatable :: cmxy       (:)   ! bulk momentum drag coefficient (m/s)
     real(kind=kp), allocatable :: chxy       (:)   ! bulk sensible heat exchange coefficient (m/s(
     real(kind=kp), allocatable :: fwetxy     (:)   ! wetted or snowed fraction of the canopy (fractional)
     real(kind=kp), allocatable :: sneqvoxy   (:)   ! snow mass at last time step (mm h2o)
     real(kind=kp), allocatable :: alboldxy   (:)   ! snow albedo at last time step (fractional)
     real(kind=kp), allocatable :: qsnowxy    (:)   ! snowfall on the ground (mm/s)
     real(kind=kp), allocatable :: wslakexy   (:)   ! lake water storage (mm)
     real(kind=kp), allocatable :: zwtxy      (:)   ! water table depth (m)
     real(kind=kp), allocatable :: waxy       (:)   ! water in the "aquifer" (mm)
     real(kind=kp), allocatable :: wtxy       (:)   ! groundwater storage (mm)
     real(kind=kp), allocatable :: lfmassxy   (:)   ! leaf mass (g/m^2)
     real(kind=kp), allocatable :: rtmassxy   (:)   ! mass of fine roots (g/m^2)
     real(kind=kp), allocatable :: stmassxy   (:)   ! stem mass (g/m^2)
     real(kind=kp), allocatable :: woodxy     (:)   ! mass of wood (incl. woody roots) (g/m^2)
     real(kind=kp), allocatable :: stblcpxy   (:)   ! stable carbon in deep soil (g/m^2)
     real(kind=kp), allocatable :: fastcpxy   (:)   ! short-lived carbon, shallow soil (g/m^2)
     real(kind=kp), allocatable :: xlaixy     (:)   ! leaf area index
     real(kind=kp), allocatable :: xsaixy     (:)   ! stem area index
     real(kind=kp), allocatable :: taussxy    (:)   ! snow age factor
     real(kind=kp), allocatable :: smcwtdxy   (:)   ! soil moisture content in the layer to the water table when deep
     real(kind=kp), allocatable :: deeprechxy (:)   ! recharge to the water table when deep
     real(kind=kp), allocatable :: rechxy     (:)   ! recharge to the water table (diagnostic) 
     real(kind=kp), allocatable :: sncovr1    (:)   ! snow cover over land (fractional)
     real(kind=kp), allocatable :: qsurf      (:)   ! specific humidity at sfc
     real(kind=kp), allocatable :: gflux      (:)   ! soil heat flux (W/m^2)
     real(kind=kp), allocatable :: drain      (:)   ! subsurface runoff (mm/s)
     real(kind=kp), allocatable :: evap       (:)   ! evaperation from latent heat flux
     real(kind=kp), allocatable :: hflx       (:)   ! sensible heat flux
     real(kind=kp), allocatable :: ep         (:)   ! potential evaporation
     real(kind=kp), allocatable :: runoff     (:)   ! surface runoff (m/s)
     real(kind=kp), allocatable :: cmm        (:)   ! ?
     real(kind=kp), allocatable :: chh        (:)   ! ?
     real(kind=kp), allocatable :: evbs       (:)   ! direct soil evaporation (m/s)
     real(kind=kp), allocatable :: evcw       (:)   ! canopy water evaporation (m/s)
     real(kind=kp), allocatable :: sbsno      (:)   ! sublimation/deposit from snopack (m/s)
     real(kind=kp), allocatable :: pah        (:)   ! precipitation advected heat - total (W/m2)
     real(kind=kp), allocatable :: ecan       (:)   ! canopy evaporation (mm/s)
     real(kind=kp), allocatable :: etran      (:)   ! transpiration (mm/s)
     real(kind=kp), allocatable :: edir       (:)   ! soil surface evaporation (mm/s)
     real(kind=kp), allocatable :: snowc      (:)   ! fractional snow cover 
     real(kind=kp), allocatable :: stm        (:)   ! total soil column moisture content (m)
     real(kind=kp), allocatable :: rhonewsn1  (:)   ! precipitation ice density (kg/m^3) 
     real(kind=kp), allocatable :: snohf      (:)   ! snow/freezing-rain latent heat flux (W/m^2)
     real(kind=kp), allocatable :: smcwlt2    (:)   ! dry soil moisture threshold
     real(kind=kp), allocatable :: smcref2    (:)   ! soil moisture threshold
     real(kind=kp), allocatable :: wet1       (:)   ! normalized soil wetness
     real(kind=kp), allocatable :: t2mmp      (:)   ! combined T2m from tiles
     real(kind=kp), allocatable :: q2mp       (:)   ! combined q2m from tiles
     real(kind=kp), allocatable :: zvfun      (:)   ! some function of vegetation used for gfs stability
     real(kind=kp), allocatable :: ztmax      (:)   ! bounded surface roughness length for heat over land
     real(kind=kp), allocatable :: rho        (:)   ! air density
     real(kind=kp), allocatable :: pores      (:)   ! max soil moisture for a given soil type for land surface model
     real(kind=kp), allocatable :: resid      (:)   ! min soil moisture for a given soil type for land surface model
     ! variables in dimensions im and km
     real(kind=kp), allocatable :: smc      (:,:)   ! total soil moisture content (fractional)
     real(kind=kp), allocatable :: stc      (:,:)   ! soil temp (K)
     real(kind=kp), allocatable :: slc      (:,:)   ! liquid soil moisture (?)
     real(kind=kp), allocatable :: tsnoxy   (:,:)   ! snow temperature (K)
     real(kind=kp), allocatable :: zsnsoxy  (:,:)   ! snow layer depth (m)
     real(kind=kp), allocatable :: snicexy  (:,:)   ! snow layer ice (mm)
     real(kind=kp), allocatable :: snliqxy  (:,:)   ! snow layer liquid water (mm)
     real(kind=kp), allocatable :: smoiseq  (:,:)   ! eq volumetric soil moisture (m^3/m^3)
  end type model_type

  ! data type for domain related variables
  type domain_type
     type(ESMF_Grid)            :: grid           ! ESMF grid object, only for mosaic case
     type(ESMF_Mesh)            :: mesh           ! ESMF mesh object
     type(ESMF_RouteHandle)     :: rh_grid2mesh   ! ESMF RouteHandle for redist from grid to mesh
     type(ESMF_RouteHandle)     :: rh_mesh2grid_r8! ESMF RouteHandle for redist from mesh to grid
     type(domain2d)             :: mosaic_domain  ! domain object created by FMS
     logical                    :: global         ! flag for global vs. regional domain
     integer                    :: ntiles         ! number of tiles in case of having CS grid
     integer                    :: ncontacts      ! number of contacts in case of having CS grid
     integer      , allocatable :: tile1      (:) ! list of tile numbers in tile 1 of each contact
     integer      , allocatable :: tile2      (:) ! list of tile numbers in tile 2 of each contact
     integer      , allocatable :: istart1    (:) ! list of starting i-index in tile 1 of each contact
     integer      , allocatable :: iend1      (:) ! list of ending i-index in tile 1 of each contact
     integer      , allocatable :: jstart1    (:) ! list of starting j-index in tile 1 of each contact
     integer      , allocatable :: jend1      (:) ! list of ending j-index in tile 1 of each contact
     integer      , allocatable :: istart2    (:) ! list of starting i-index in tile 2 of each contact
     integer      , allocatable :: iend2      (:) ! list of ending i-index in tile 2 of each contact
     integer      , allocatable :: jstart2    (:) ! list of starting j-index in tile 2 of each contact
     integer      , allocatable :: jend2      (:) ! list of ending j-index in tile 2 of each contact
     integer                    :: ni             ! global size in i direction, only for mesh as input
     integer                    :: nj             ! global size in j direction, only for mesh as input
     integer      , allocatable :: nit        (:) ! size of tile in i direction
     integer      , allocatable :: njt        (:) ! size of tile in j direction
     real(kind=kp), allocatable :: latt     (:,:) ! mosaic latitude
     real(kind=kp), allocatable :: lont     (:,:) ! mosaic longitude
     integer                    :: begl           ! starting index of size
     integer                    :: endl           ! ending index of size 
     integer                    :: layout     (2) ! layout for domain decomposition
     real(kind=kp), allocatable :: hgt        (:) ! topography
     integer      , allocatable :: mask       (:) ! mesh land mask: 1 = land, 0 = ocean
     real(kind=kp), allocatable :: frac       (:) ! mesh fractional land
     real(kind=kp), allocatable :: lats       (:) ! mesh latitude
     real(kind=kp), allocatable :: lons       (:) ! mesh longitude
     real(kind=kp), allocatable :: garea      (:) ! mesh cell area
  end type domain_type

  ! data type for namelist options
  ! TODO: this could be merged with static_type
  type namelist_type
     character*100              :: case_name                         ! name of case
     character*255              :: mosaic_file                       ! name of mosaic file
     character*255              :: input_dir                         ! input directory for tiled files
     character*255              :: restart_dir                       ! restart directory
     character*255              :: restart_file                      ! restart file name
     character*255              :: ic_type                           ! source of initial conditions, custom vs. sfc
     logical                    :: restart_run                       ! flag for restart run
     integer                    :: num_soil_levels                   ! number of soil levels
     real(kind=kp), allocatable :: soil_level_thickness(:)           ! soil level thicknesses (m)
     real(kind=kp), allocatable :: soil_level_nodes(:)               ! soil level centroids from surface (m)
     integer                    :: dynamic_vegetation_option         ! choice for dynamic vegetation option
     integer                    :: canopy_stomatal_resistance_option ! choice for canopy stomatal resistance option
     integer                    :: soil_wetness_option               ! 
     integer                    :: runoff_option
     integer                    :: surface_exchange_option
     integer                    :: supercooled_soilwater_option
     integer                    :: frozen_soil_adjust_option
     integer                    :: radiative_transfer_option
     integer                    :: snow_albedo_option
     integer                    :: precip_partition_option
     integer                    :: soil_temp_lower_bdy_option
     integer                    :: soil_temp_time_scheme_option
     integer                    :: surface_evap_resistance_option
     integer                    :: glacier_option
     integer                    :: surface_thermal_roughness_option
     integer                    :: output_freq                       ! model output interval
     logical                    :: has_export                        ! enable/disable export fields
     logical                    :: calc_snet                         ! enable/disable calculating net shortwave rad. internally
     integer                    :: soil_type_category                ! soil type (category)
     integer                    :: veg_type_category                 ! vegetation type (category)
     real(kind=kp)              :: initial_emiss                     ! initial value for the emissivity (constant in everywhere)
     real(kind=kp)              :: initial_albedo                    ! initial value for the monthly albedo (constant in everywhere)
  end type namelist_type

  type noahmp_type
     type(initial_type)  :: init
     type(forcing_type)  :: forc
     type(static_type)   :: static
     type(model_type)    :: model
     type(domain_type)   :: domain
     type(namelist_type) :: nmlist

  contains

     procedure, public  :: Initialize
     procedure, private :: InitializeAllocate
     procedure, private :: InitializeDefault
     procedure, public  :: InitializeStates

  end type noahmp_type

  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
     logical :: connected = .false.
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type(fld_list_type)    :: fldsToLnd(fldsMax)
  type(fld_list_type)    :: fldsFrLnd(fldsMax)

contains

  subroutine Initialize(this, begl, endl, km, lsnowl)

    class(noahmp_type) :: this
    integer :: begl, endl, km

    call this%InitializeAllocate(begl, endl, km, lsnowl)
    call this%InitializeDefault()

  end subroutine Initialize 

  subroutine InitializeAllocate(this, begl, endl, km, lsnowl)

    class(noahmp_type) :: this
    integer :: begl, endl

    allocate(this%init%snow_water_equivalent(begl:endl))
    allocate(this%init%snow_depth(begl:endl))
    allocate(this%init%canopy_water(begl:endl))
    allocate(this%init%skin_temperature(begl:endl))
    allocate(this%init%soil_temperature(begl:endl,km))
    allocate(this%init%soil_moisture(begl:endl,km))
    allocate(this%init%soil_liquid(begl:endl,km))
    allocate(this%init%surface_roughness(begl:endl))
    allocate(this%init%friction_velocity(begl:endl))

    allocate(this%forc%t1     (begl:endl))
    allocate(this%forc%q1     (begl:endl))
    allocate(this%forc%u1     (begl:endl))
    allocate(this%forc%v1     (begl:endl))
    allocate(this%forc%ps     (begl:endl))
    allocate(this%forc%pbot   (begl:endl))
    allocate(this%forc%tskin  (begl:endl))
    allocate(this%forc%dlwflx (begl:endl))
    allocate(this%forc%dswsfc (begl:endl))
    allocate(this%forc%snet   (begl:endl))
    allocate(this%forc%wind   (begl:endl))
    allocate(this%forc%tprcp  (begl:endl))
    allocate(this%forc%tprcpc (begl:endl))
    allocate(this%forc%tprcpl (begl:endl))
    allocate(this%forc%snow   (begl:endl))
    allocate(this%forc%snowc  (begl:endl))
    allocate(this%forc%snowl  (begl:endl))
    allocate(this%forc%vegfrac(begl:endl))
    allocate(this%forc%hgt    (begl:endl))
    allocate(this%forc%prslk1 (begl:endl))
    allocate(this%forc%ustar1 (begl:endl))
    allocate(this%forc%zorl   (begl:endl))

    allocate(this%model%u1         (begl:endl))
    allocate(this%model%v1         (begl:endl))
    allocate(this%model%soiltyp    (begl:endl))
    allocate(this%model%vegtype    (begl:endl))
    allocate(this%model%sigmaf     (begl:endl))
    allocate(this%model%emiss      (begl:endl))
    allocate(this%model%albdvis    (begl:endl))
    allocate(this%model%albdnir    (begl:endl))
    allocate(this%model%albivis    (begl:endl))
    allocate(this%model%albinir    (begl:endl))
    allocate(this%model%snet       (begl:endl))
    allocate(this%model%tg3        (begl:endl))
    allocate(this%model%cm         (begl:endl))
    allocate(this%model%ch         (begl:endl))
    allocate(this%model%prsl1      (begl:endl))
    allocate(this%model%prslki     (begl:endl))
    allocate(this%model%prslk1     (begl:endl))
    allocate(this%model%prsik1     (begl:endl))
    allocate(this%model%zf         (begl:endl))
    allocate(this%model%pblh       (begl:endl))
    allocate(this%model%dry        (begl:endl))
    allocate(this%model%slopetyp   (begl:endl))
    allocate(this%model%alb_monthly(begl:endl,12))
    allocate(this%model%gvf_monthly(begl:endl,12))
    allocate(this%model%shdmin     (begl:endl))
    allocate(this%model%shdmax     (begl:endl))
    allocate(this%model%snoalb     (begl:endl))
    allocate(this%model%sfalb      (begl:endl))
    allocate(this%model%flag_iter  (begl:endl))
    allocate(this%model%xlatin     (begl:endl))
    allocate(this%model%xcoszin    (begl:endl))
    allocate(this%model%rainn_mp   (begl:endl))
    allocate(this%model%rainc_mp   (begl:endl))
    allocate(this%model%snow_mp    (begl:endl))
    allocate(this%model%graupel_mp (begl:endl))
    allocate(this%model%ice_mp     (begl:endl))
    allocate(this%model%tprcp      (begl:endl))
    allocate(this%model%weasd      (begl:endl))
    allocate(this%model%snwdph     (begl:endl))
    allocate(this%model%tskin      (begl:endl))
    allocate(this%model%srflag     (begl:endl))
    allocate(this%model%canopy     (begl:endl))
    allocate(this%model%trans      (begl:endl))
    allocate(this%model%tsurf      (begl:endl))
    allocate(this%model%zorl       (begl:endl))
    allocate(this%model%rb1        (begl:endl))
    allocate(this%model%fm1        (begl:endl))
    allocate(this%model%fh1        (begl:endl))
    allocate(this%model%ustar1     (begl:endl))
    allocate(this%model%stress1    (begl:endl))
    allocate(this%model%fm101      (begl:endl))
    allocate(this%model%fh21       (begl:endl))
    allocate(this%model%rmol1      (begl:endl))
    allocate(this%model%flhc1      (begl:endl))
    allocate(this%model%flqc1      (begl:endl))
    allocate(this%model%snowxy     (begl:endl))
    allocate(this%model%tvxy       (begl:endl))
    allocate(this%model%tgxy       (begl:endl))
    allocate(this%model%canicexy   (begl:endl))
    allocate(this%model%canliqxy   (begl:endl))
    allocate(this%model%eahxy      (begl:endl))
    allocate(this%model%tahxy      (begl:endl))
    allocate(this%model%cmxy       (begl:endl))
    allocate(this%model%chxy       (begl:endl))
    allocate(this%model%fwetxy     (begl:endl))
    allocate(this%model%sneqvoxy   (begl:endl))
    allocate(this%model%alboldxy   (begl:endl))
    allocate(this%model%qsnowxy    (begl:endl))
    allocate(this%model%wslakexy   (begl:endl))
    allocate(this%model%zwtxy      (begl:endl))
    allocate(this%model%waxy       (begl:endl))
    allocate(this%model%wtxy       (begl:endl))
    allocate(this%model%lfmassxy   (begl:endl))
    allocate(this%model%rtmassxy   (begl:endl))
    allocate(this%model%stmassxy   (begl:endl))
    allocate(this%model%woodxy     (begl:endl))
    allocate(this%model%stblcpxy   (begl:endl))
    allocate(this%model%fastcpxy   (begl:endl))
    allocate(this%model%xlaixy     (begl:endl))
    allocate(this%model%xsaixy     (begl:endl))
    allocate(this%model%taussxy    (begl:endl))
    allocate(this%model%smcwtdxy   (begl:endl))
    allocate(this%model%deeprechxy (begl:endl))
    allocate(this%model%rechxy     (begl:endl))
    allocate(this%model%sncovr1    (begl:endl))
    allocate(this%model%qsurf      (begl:endl))
    allocate(this%model%gflux      (begl:endl))
    allocate(this%model%drain      (begl:endl))
    allocate(this%model%evap       (begl:endl))
    allocate(this%model%hflx       (begl:endl))
    allocate(this%model%ep         (begl:endl))
    allocate(this%model%runoff     (begl:endl))
    allocate(this%model%cmm        (begl:endl))
    allocate(this%model%chh        (begl:endl))
    allocate(this%model%evbs       (begl:endl))
    allocate(this%model%evcw       (begl:endl))
    allocate(this%model%sbsno      (begl:endl))
    allocate(this%model%pah        (begl:endl))
    allocate(this%model%ecan       (begl:endl))
    allocate(this%model%etran      (begl:endl))
    allocate(this%model%edir       (begl:endl))
    allocate(this%model%snowc      (begl:endl))
    allocate(this%model%stm        (begl:endl))
    allocate(this%model%rhonewsn1  (begl:endl))
    allocate(this%model%snohf      (begl:endl))
    allocate(this%model%smcwlt2    (begl:endl))
    allocate(this%model%smcref2    (begl:endl))
    allocate(this%model%wet1       (begl:endl))
    allocate(this%model%t2mmp      (begl:endl))
    allocate(this%model%q2mp       (begl:endl))
    allocate(this%model%zvfun      (begl:endl))
    allocate(this%model%ztmax      (begl:endl))
    allocate(this%model%rho        (begl:endl))
    allocate(this%model%pores      (30))
    allocate(this%model%resid      (30))

    allocate(this%model%smc        (begl:endl,km))
    allocate(this%model%stc        (begl:endl,km))
    allocate(this%model%slc        (begl:endl,km))
    allocate(this%model%smoiseq    (begl:endl,km))
    allocate(this%model%tsnoxy     (begl:endl,lsnowl:0))
    allocate(this%model%zsnsoxy    (begl:endl,lsnowl:km))
    allocate(this%model%snicexy    (begl:endl,lsnowl:0))
    allocate(this%model%snliqxy    (begl:endl,lsnowl:0))

  end subroutine InitializeAllocate

  subroutine InitializeDefault(this)

    class(noahmp_type) :: this

    this%init%snow_water_equivalent = 0.0_kp
    this%init%snow_depth            = 0.0_kp
    this%init%canopy_water          = 0.0_kp
    this%init%skin_temperature      = 0.0_kp
    this%init%soil_temperature      = 0.0_kp
    this%init%soil_moisture         = 0.0_kp
    this%init%soil_liquid           = 0.0_kp
    this%init%surface_roughness     = 0.0_kp
    this%init%friction_velocity     = 0.0_kp

    this%forc%t1      = 0.0_kp
    this%forc%q1      = 0.0_kp
    this%forc%u1      = 0.0_kp
    this%forc%v1      = 0.0_kp
    this%forc%ps      = 0.0_kp
    this%forc%pbot    = 0.0_kp
    this%forc%tskin   = 0.0_kp
    this%forc%dlwflx  = 0.0_kp
    this%forc%dswsfc  = 0.0_kp
    this%forc%snet    = 0.0_kp
    this%forc%wind    = 0.0_kp
    this%forc%tprcp   = 0.0_kp
    this%forc%tprcpc  = 0.0_kp
    this%forc%tprcpl  = 0.0_kp
    this%forc%snow    = 0.0_kp
    this%forc%snowc   = 0.0_kp
    this%forc%snowl   = 0.0_kp
    this%forc%vegfrac = 0.0_kp
    this%forc%hgt     = 0.0_kp
    this%forc%prslk1  = 0.0_kp
    this%forc%ustar1  = 0.0_kp
    this%forc%zorl    = 0.0_kp

    this%model%iyrlen      = 0
    this%model%julian      = 0.0_kp
    this%model%u1          = 0.0_kp
    this%model%v1          = 0.0_kp
    this%model%soiltyp     = 0
    this%model%vegtype     = 0
    this%model%sigmaf      = 0.0_kp
    this%model%emiss       = 0.0_kp
    this%model%albdvis     = 0.0_kp
    this%model%albdnir     = 0.0_kp
    this%model%albivis     = 0.0_kp
    this%model%albinir     = 0.0_kp
    this%model%snet        = 0.0_kp
    this%model%tg3         = 0.0_kp
    this%model%cm          = 0.0_kp
    this%model%ch          = 0.0_kp
    this%model%prsl1       = 0.0_kp
    this%model%prslki      = 0.0_kp
    this%model%prslk1      = 0.0_kp
    this%model%prsik1      = 0.0_kp
    this%model%zf          = 0.0_kp
    this%model%pblh        = 0.0_kp
    this%model%dry         = .false.
    this%model%slopetyp    = 0
    this%model%alb_monthly = 0.0_kp
    this%model%gvf_monthly = 0.0_kp
    this%model%shdmin      = 0.0_kp
    this%model%shdmax      = 0.0_kp
    this%model%snoalb      = 0.0_kp
    this%model%sfalb       = 0.0_kp
    this%model%flag_iter   = .false.
    this%model%xlatin      = 0.0_kp
    this%model%xcoszin     = 0.0_kp
    this%model%rainn_mp    = 0.0_kp
    this%model%rainc_mp    = 0.0_kp
    this%model%snow_mp     = 0.0_kp
    this%model%graupel_mp  = 0.0_kp
    this%model%ice_mp      = 0.0_kp
    this%model%tprcp       = 0.0_kp
    this%model%weasd       = 0.0_kp
    this%model%snwdph      = 0.0_kp
    this%model%tskin       = 0.0_kp
    this%model%srflag      = 0.0_kp
    this%model%canopy      = 0.0_kp
    this%model%trans       = 0.0_kp
    this%model%tsurf       = 0.0_kp
    this%model%zorl        = 0.0_kp
    this%model%rb1         = 0.0_kp  
    this%model%fm1         = 0.0_kp   
    this%model%fh1         = 0.0_kp
    this%model%ustar1      = 0.0_kp
    this%model%stress1     = 0.0_kp
    this%model%fm101       = 0.0_kp
    this%model%fh21        = 0.0_kp
    this%model%rmol1       = 0.0_kp
    this%model%flhc1       = 0.0_kp
    this%model%flqc1       = 0.0_kp
    this%model%snowxy      = 0.0_kp
    this%model%tvxy        = 0.0_kp
    this%model%tgxy        = 0.0_kp
    this%model%canicexy    = 0.0_kp
    this%model%canliqxy    = 0.0_kp
    this%model%eahxy       = 0.0_kp
    this%model%tahxy       = 0.0_kp
    this%model%cmxy        = 0.0_kp
    this%model%chxy        = 0.0_kp
    this%model%fwetxy      = 0.0_kp
    this%model%sneqvoxy    = 0.0_kp
    this%model%alboldxy    = 0.0_kp
    this%model%qsnowxy     = 0.0_kp
    this%model%wslakexy    = 0.0_kp
    this%model%zwtxy       = 0.0_kp
    this%model%waxy        = 0.0_kp
    this%model%wtxy        = 0.0_kp
    this%model%lfmassxy    = 0.0_kp
    this%model%rtmassxy    = 0.0_kp
    this%model%stmassxy    = 0.0_kp
    this%model%woodxy      = 0.0_kp
    this%model%stblcpxy    = 0.0_kp
    this%model%fastcpxy    = 0.0_kp
    this%model%xlaixy      = 0.0_kp
    this%model%xsaixy      = 0.0_kp
    this%model%taussxy     = 0.0_kp
    this%model%smoiseq     = 0.0_kp
    this%model%smcwtdxy    = 0.0_kp
    this%model%deeprechxy  = 0.0_kp
    this%model%rechxy      = 0.0_kp
    this%model%sncovr1     = 0.0_kp
    this%model%qsurf       = 0.0_kp
    this%model%gflux       = 0.0_kp
    this%model%drain       = 0.0_kp
    this%model%evap        = 0.0_kp
    this%model%hflx        = 0.0_kp
    this%model%ep          = 0.0_kp
    this%model%runoff      = 0.0_kp
    this%model%cmm         = 0.0_kp
    this%model%chh         = 0.0_kp
    this%model%evbs        = 0.0_kp
    this%model%evcw        = 0.0_kp
    this%model%sbsno       = 0.0_kp
    this%model%pah         = 0.0_kp
    this%model%ecan        = 0.0_kp
    this%model%etran       = 0.0_kp
    this%model%edir        = 0.0_kp
    this%model%snowc       = 0.0_kp
    this%model%stm         = 0.0_kp
    this%model%rhonewsn1   = 0.0_kp
    this%model%snohf       = 0.0_kp
    this%model%smcwlt2     = 0.0_kp
    this%model%smcref2     = 0.0_kp
    this%model%wet1        = 0.0_kp
    this%model%t2mmp       = 0.0_kp
    this%model%q2mp        = 0.0_kp
    this%model%zvfun       = 0.0_kp
    this%model%ztmax       = 0.0_kp
    this%model%rho         = 0.0_kp
    this%model%pores       = 0.0_kp
    this%model%resid       = 0.0_kp
    this%model%smc         = 0.0_kp
    this%model%stc         = 0.0_kp
    this%model%slc         = 0.0_kp
    this%model%smoiseq     = 0.0_kp
    this%model%tsnoxy      = 0.0_kp
    this%model%zsnsoxy     = 0.0_kp
    this%model%snicexy     = 0.0_kp
    this%model%snliqxy     = 0.0_kp

  end subroutine InitializeDefault

  subroutine InitializeStates(this, namelist, static, month)
    
    ! out of necessity, this is extracted from FV3GFS_io.F90
    
    use noahmp_tables, only: isbarren_table, isice_table, isurban_table
    use noahmp_tables, only: iswater_table, laim_table, saim_table, sla_table
      
    class(noahmp_type)  :: this
    type(namelist_type) :: namelist
    type(static_type)   :: static

    integer             :: iloc, ilevel, isnow
    real                :: masslai, masssai, depthm
    real                :: dzsno(static%lsnowl:static%km)
    
    do iloc = 1, static%im
      this%model%tvxy (iloc) = this%model%tskin(iloc)
      this%model%tgxy (iloc) = this%model%tskin(iloc)
      this%model%tahxy(iloc) = this%model%tskin(iloc)
  
      if (this%model%snwdph(iloc) > 0.01 .and. this%model%tskin(iloc) > 273.15 ) then
         this%model%tvxy (iloc) = 273.15
         this%model%tgxy (iloc) = 273.15
         this%model%tahxy(iloc) = 273.15
      end if
    
      this%model%canicexy(iloc) = 0.0
      this%model%canliqxy(iloc) = this%model%canopy(iloc)
      this%model%eahxy(iloc)    = 2000.0
      this%model%cmxy    (iloc) = 0.0
      this%model%chxy    (iloc) = 0.0
      this%model%fwetxy  (iloc) = 0.0
      this%model%sneqvoxy(iloc) = this%model%weasd(iloc)
      this%model%alboldxy(iloc) = 0.65
      this%model%qsnowxy (iloc) = 0.0
      this%model%wslakexy(iloc) = 0.0
      this%model%taussxy (iloc) = 0.0
      this%model%waxy    (iloc) = 4900.0 ! this assumes water table is at 2.5m
      this%model%wtxy    (iloc) = this%model%waxy(iloc)
      this%model%zwtxy   (iloc) = (25.0 + 2.0)-this%model%waxy(iloc)/1000.0/0.2
  
      if ((this%model%vegtype(iloc) == isbarren_table) .or. (this%model%vegtype(iloc) == isice_table) .or. &
          (this%model%vegtype(iloc) == isurban_table)  .or. (this%model%vegtype(iloc) == iswater_table) .or. &
          (this%model%vegtype(iloc) < 0)) then
         this%model%xlaixy  (iloc) = 0.0
         this%model%xsaixy  (iloc) = 0.0
         this%model%lfmassxy(iloc) = 0.0
         this%model%stmassxy(iloc) = 0.0
         this%model%rtmassxy(iloc) = 0.0
         this%model%woodxy  (iloc) = 0.0       
         this%model%stblcpxy(iloc) = 0.0      
         this%model%fastcpxy(iloc) = 0.0     
      else
         this%model%xlaixy  (iloc) = max(laim_table(this%model%vegtype(iloc), month),0.05)
         this%model%xsaixy  (iloc) = max(this%model%xlaixy(iloc)*0.1,0.05)
         masslai             = 1000.0/max(sla_table(this%model%vegtype(iloc)),1.0)
         this%model%lfmassxy(iloc) = this%model%xlaixy(iloc)*masslai
         masssai             = 1000.0/3.0
         this%model%stmassxy(iloc) = this%model%xsaixy(iloc)*masssai
         this%model%rtmassxy(iloc) = 500.0      
         this%model%woodxy  (iloc) = 500.0       
         this%model%stblcpxy(iloc) = 1000.0      
         this%model%fastcpxy(iloc) = 1000.0     
      end if
  
      if (this%model%vegtype(iloc)  == isice_table )  then
         do ilevel = 1, namelist%num_soil_levels
            this%model%stc(iloc,ilevel) = min(this%model%stc(iloc,ilevel),min(this%model%tg3(iloc),263.15))
         end do
         this%model%smc(iloc,:) = 1.0
         this%model%slc(iloc,:) = 0.0
      end if
  
      depthm = this%model%snwdph(iloc)/1000.0  ! snow depth in meters
  
      if (this%model%weasd(iloc) /= 0.0 .and. depthm == 0.0 ) then
        depthm = this%model%weasd(iloc)/1000.0*10.0 ! assume 10/1 ratio
      endif
  
      if (this%model%vegtype(iloc) == isice_table) then  ! land ice in MODIS/IGBP
         if (this%model%weasd(iloc) < 0.1) then
            this%model%weasd(iloc) = 0.1
            depthm = this%model%weasd(iloc)/1000.0*10.0 
         end if
      end if
  
      dzsno = 0.0
      if (depthm < 0.025 ) then
         this%model%snowxy(iloc) = 0.0
         dzsno(-2:0) = 0.0
      elseif (depthm >= 0.025 .and. depthm <= 0.05 ) then
         this%model%snowxy(iloc) = -1.0
         dzsno(0) = depthm
      elseif (depthm > 0.05 .and. depthm <= 0.10 ) then
         this%model%snowxy(iloc) = -2.0
         dzsno(-1) = 0.5*depthm
         dzsno(0) = 0.5*depthm
      elseif (depthm > 0.10 .and. depthm <= 0.25 ) then
         this%model%snowxy(iloc) = -2.0
         dzsno(-1) = 0.05
         dzsno(0) = depthm - 0.05
      elseif (depthm > 0.25 .and. depthm <= 0.45 ) then
         this%model%snowxy(iloc) = -3.0
         dzsno(-2) = 0.05
         dzsno(-1) = 0.5*(depthm-0.05)
         dzsno(0) = 0.5*(depthm-0.05)
      elseif (depthm > 0.45) then 
         this%model%snowxy(iloc) = -3.0
         dzsno(-2) = 0.05
         dzsno(-1) = 0.20
         dzsno(0) = depthm - 0.05 - 0.20
      endif
  
      ! Now we have the snowxy field
      ! snice + snliq + tsno allocation and compute them from what we have
      this%model%tsnoxy (iloc,-2:0) = 0.0
      this%model%snicexy(iloc,-2:0) = 0.0
      this%model%snliqxy(iloc,-2:0) = 0.0
      this%model%zsnsoxy(iloc,-2:namelist%num_soil_levels) = 0.0
  
      isnow = nint(this%model%snowxy(iloc))+1 ! snowxy <=0.0, dzsno >= 0.0
  
      do ilevel = isnow, 0
         this%model%tsnoxy (iloc,ilevel) = this%model%tgxy(iloc)
         this%model%snliqxy(iloc,ilevel) = 0.0
         this%model%snicexy(iloc,ilevel) = dzsno(ilevel)/depthm*this%model%weasd(iloc)
      end do
  
      this%model%zsnsoxy(iloc,isnow) = -1.0*dzsno(isnow)
      do ilevel = isnow+1, 0
         this%model%zsnsoxy(iloc,ilevel) = this%model%zsnsoxy(iloc,ilevel-1)-dzsno(ilevel)
      end do
      do ilevel = 1, namelist%num_soil_levels
         this%model%zsnsoxy(iloc,ilevel) = this%model%zsnsoxy(iloc,ilevel-1)-namelist%soil_level_thickness(ilevel)
      end do
  
      ! TODO: Should not need to initialize these
      this%model%smoiseq(iloc,:)  = 0.0
      this%model%smcwtdxy(iloc)   = 0.0
      this%model%deeprechxy(iloc) = 0.0
      this%model%rechxy(iloc)     = 0.0
    end do ! iloc
    
  end subroutine InitializeStates
  
end module lnd_comp_types
