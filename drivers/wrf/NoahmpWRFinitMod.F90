module NoahmpWRFinitMod

! --------------------------------------------------------------------------
! this is for NoahmpIO variable mapping & initialization in WRF physics_init
! and calling the main NoahMP initialization module: NoahmpInitMain(NoahmpIO) 
! adapted from NOAHMP_INIT in original module_sf_noahmpdrv.F file
!
! Coder: Cenlin He (NCAR), December 2025
! ---------------------------------------------------------------------------

contains

  subroutine NoahmpWRFinit(NoahmpIO, MMINLU, SNOW, SNOWH, CANWAT, ISLTYP, IVGTYP, XLAT, &
                   TSLB,  SMOIS, SH2O,   DZS, FNDSOILW, FNDSNOWH,                       &
                   TSK, isnowxy, tvxy,  tgxy, canicexy,      TMN,   XICE,               &
                   canliqxy,    eahxy, tahxy,     cmxy,     chxy,                       &
                   fwetxy, sneqvoxy, alboldxy, qsnowxy, qrainxy, wslakexy, zwtxy, waxy, &
                   wtxy, tsnoxy, zsnsoxy, snicexy, snliqxy, lfmassxy, rtmassxy,         &
                   stmassxy, woodxy, stblcpxy, fastcpxy, xsaixy, lai,                   &
                   grainxy,   gddxy,                                                    &
                   croptype, cropcat,                                                   &
                   irnumsi, irnummi, irnumfi, irwatsi,                                  &
                   irwatmi, irwatfi, ireloss, irsivol,                                  &
                   irmivol, irfivol, irrsplh,                                           &
                   t2mvxy,   t2mbxy, chstarxy, fsatxy, wsurfxy,                         &
                   snrdsxy, snfrxy, bcphixy, bcphoxy, ocphixy, ocphoxy, dust1xy,        &
                   dust2xy, dust3xy, dust4xy, dust5xy, massconcbcphixy, massconcbcphoxy,&
                   massconcocphixy, massconcocphoxy, massconcdust1xy, massconcdust2xy,  &
                   massconcdust3xy, massconcdust4xy, massconcdust5xy,                   &
                   ALBSOILDIRXY, ALBSOILDIFXY,                                          &
                   NSOIL,   restart,                                                    &
                   allowed_to_read , IOPT_RUNSUB, IOPT_CROP, IOPT_IRR, IOPT_IRRM,       &
                   SF_URBAN_PHYSICS, IOPT_SOIL, IOPT_ALB, IOPT_WETLAND,                 &
                   SNICAR_SNOWOPTICS_OPT, SNICAR_DUSTOPTICS_OPT, SNICAR_SOLARSPEC_OPT,  & ! optional SNICAR option
                   SNICAR_BANDNUMBER_OPT,                                               & ! optional SNICAR option
                   ids,ide, jds,jde, kds,kde,                                           &
                   ims,ime, jms,jme, kms,kme,                                           &
                   its,ite, jts,jte, kts,kte,                                           &
                   smoiseq,smcwtdxy,rechxy,deeprechxy,qtdrain,areaxy,dx,dy,msftx,msfty, & ! Optional groundwater
                   wtddt,   stepwtd, dt, qrfsxy, qspringsxy, qslatxy,                   & ! Optional groundwater
                   fdepthxy, ht, riverbedxy, eqzwt, rivercondxy, pexpxy, rechclim       ) ! Optional groundwater

! ---------------------------------------------------------------------------------------

    use NoahmpIOVarType
    use NoahmpIOVarInitMod
    use NoahmpReadTableMod
    use NoahmpInitMainMod
    use SnowInputSnicarMod

    implicit none

    type(NoahmpIO_type), intent(inout)                         :: NoahmpIO

    ! input only
    INTEGER, INTENT(IN)                                        :: ids,ide, jds,jde, kds,kde,  &
                                                                  ims,ime, jms,jme, kms,kme,  &
                                                                  its,ite, jts,jte, kts,kte
    INTEGER, INTENT(IN)                                        :: NSOIL, IOPT_RUNSUB, IOPT_CROP, &
                                                                  IOPT_IRR, IOPT_IRRM, IOPT_SOIL,&
                                                                  IOPT_ALB, IOPT_WETLAND
    LOGICAL, INTENT(IN)                                        :: restart, allowed_to_read
    LOGICAL, INTENT(IN)                                        :: FNDSOILW, FNDSNOWH
    INTEGER, INTENT(IN)                                        :: SF_URBAN_PHYSICS                      
    CHARACTER(LEN=*),                    INTENT(IN)            :: MMINLU
    REAL,    DIMENSION(NSOIL), INTENT(IN)                      :: DZS                 ! Thickness of the soil layers [m]
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: ISLTYP, IVGTYP
    REAL   , DIMENSION(ims:ime,5,jms:jme),INTENT(IN)           :: croptype            ! crop type fraction
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: XLAT                ! latitude
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: TSK                 ! skin temperature (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN)            :: XICE                ! sea ice fraction
    REAL,                                INTENT(IN), OPTIONAL  :: DT, WTDDT
    REAL,                                INTENT(IN), OPTIONAL  :: DX, DY
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: FDEPTHXY            ! efolding depth for transmissivity (m)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: HT                  ! terrain height (m)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: MSFTX, MSFTY
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(IN), OPTIONAL  :: rechclim
    INTEGER, INTENT(IN),                             OPTIONAL  :: SNICAR_SNOWOPTICS_OPT, &
                                                                  SNICAR_DUSTOPTICS_OPT, &
                                                                  SNICAR_SOLARSPEC_OPT,  &
                                                                  SNICAR_BANDNUMBER_OPT
    ! in/out
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: TMN                 ! deep soil temperature (k)
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: isnowxy             ! actual no. of snow layers
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(INOUT) :: SMOIS, SH2O, TSLB
    REAL,    DIMENSION(ims:ime, jms:jme), INTENT(INOUT)        :: SNOW, SNOWH, CANWAT
    REAL,    DIMENSION(ims:ime,-2:NSOIL,jms:jme),INTENT(INOUT) :: zsnsoxy             ! snow layer depth [m]
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: tsnoxy              ! snow temperature [K]
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snicexy             ! snow layer ice [mm]
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snliqxy             ! snow layer liquid water [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tvxy                ! vegetation canopy temperature
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tgxy                ! ground surface temperature
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: canicexy            ! canopy-intercepted ice (mm)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: canliqxy            ! canopy-intercepted liquid water (mm)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: eahxy               ! canopy air vapor pressure (pa)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: tahxy               ! canopy air temperature (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cmxy                ! momentum drag coefficient
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chxy                ! sensible heat exchange coefficient
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fwetxy              ! wetted or snowed fraction of the canopy (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: sneqvoxy            ! snow mass at last time step(mm h2o)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: alboldxy            ! snow albedo at last time step (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qsnowxy             ! snowfall on the ground [mm/s]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qrainxy             ! rainfall on the ground [mm/s]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wslakexy            ! lake water storage [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: zwtxy               ! water table depth [m]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: waxy                ! water in the "aquifer" [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wtxy                ! groundwater storage [mm]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: lfmassxy            ! leaf mass [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: rtmassxy            ! mass of fine roots [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: stmassxy            ! stem mass [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: woodxy              ! mass of wood (incl. woody roots) [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: grainxy             ! mass of grain [g/m2] !XING
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gddxy               ! growing degree days !XING
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: stblcpxy            ! stable carbon in deep soil [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fastcpxy            ! short-lived carbon, shallow soil [g/m2]
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: xsaixy              ! stem area index
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: lai                 ! leaf area index
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: qtdrain             ! tile drainage (mm)
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnumsi             ! irrigation number
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnummi             ! irrigation number
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irnumfi             ! irrigation number
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatsi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatmi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irwatfi             ! irrigation water amount
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: ireloss             ! irrigation loss
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irsivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irmivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irfivol             ! irrigation water volume
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: irrsplh             ! irrigation evaporation heat
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: t2mvxy              ! 2m temperature vegetation part (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: t2mbxy              ! 2m temperature bare ground part (k)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chstarxy            ! dummy
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: fsatxy              ! saturation fraction of grid (-)
    REAL,    DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wsurfxy             ! wetland water storage (mm)
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snrdsxy             ! SNICAR snow radius
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: snfrxy              ! SNICAR snow freezing rate
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: bcphixy             ! SNICAR BCPHI mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: bcphoxy             ! SNICAR BCPHO mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: ocphixy             ! SNICAR OCPHI mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: ocphoxy             ! SNICAR OCPHO mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust1xy             ! SNICAR DUST1 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust2xy             ! SNICAR DUST2 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust3xy             ! SNICAR DUST3 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust4xy             ! SNICAR DUST4 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: dust5xy             ! SNICAR DUST5 mass in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcbcphixy     ! SNICAR BCPHI mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcbcphoxy     ! SNICAR BCPHO mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcocphixy     ! SNICAR OCPHI mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcocphoxy     ! SNICAR OCPHO mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust1xy     ! SNICAR DUST1 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust2xy     ! SNICAR DUST2 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust3xy     ! SNICAR DUST3 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust4xy     ! SNICAR DUST4 mass conc in snow
    REAL,    DIMENSION(ims:ime,-2:0,jms:jme), INTENT(INOUT)    :: massconcdust5xy     ! SNICAR DUST5 mass conc in snow
    REAL,    DIMENSION(ims:ime,1:2,jms:jme),  INTENT(INOUT)    :: ALBSOILDIRXY        ! soil albedo direct
    REAL,    DIMENSION(ims:ime,1:2,jms:jme),  INTENT(INOUT)    :: ALBSOILDIFXY        ! soil albedo diffuse
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: smcwtdxy     ! deep soil moisture content [m3m-3]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: deeprechxy   ! deep recharge [m]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: rechxy       ! accumulated recharge [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qrfsxy       ! accumulated flux from groundwater to rivers [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qspringsxy   ! accumulated seeping water [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: qslatxy      ! accumulated lateral flow [mm]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: areaxy       ! grid cell area [m2]
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: RIVERBEDXY   ! riverbed depth (m)
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: EQZWT        ! equilibrium water table depth (m)
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: RIVERCONDXY  ! river conductance
    REAL,    DIMENSION(ims:ime,jms:jme),         INTENT(INOUT), OPTIONAL :: PEXPXY       ! factor for river conductance
    REAL,    DIMENSION(ims:ime,1:NSOIL,jms:jme), INTENT(INOUT), OPTIONAL :: smoiseq      ! equilibrium soil moisture content [m3m-3]
    ! output only
    INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(OUT)           :: cropcat             ! crop type
    INTEGER,                             INTENT(OUT), OPTIONAL :: STEPWTD
 
! ----------------------------------------------------------------------------------

    ! initialize NoahmpIO dimension and key config variables
    NoahmpIO%xstart           = ims
    NoahmpIO%xend             = ime
    NoahmpIO%ystart           = jms
    NoahmpIO%yend             = jme
    NoahmpIO%ids              = ids
    NoahmpIO%ide              = ide
    NoahmpIO%jds              = jds
    NoahmpIO%jde              = jde
    NoahmpIO%kds              = kds
    NoahmpIO%kde              = kde
    NoahmpIO%ims              = ims
    NoahmpIO%ime              = ime
    NoahmpIO%jms              = jms
    NoahmpIO%jme              = jme
    NoahmpIO%kms              = kms
    NoahmpIO%kme              = kme
    NoahmpIO%its              = its
    NoahmpIO%ite              = ite
    NoahmpIO%jts              = jts
    NoahmpIO%jte              = jte
    NoahmpIO%kts              = kts
    NoahmpIO%kte              = kte
    NoahmpIO%NSOIL            = NSOIL
    NoahmpIO%LLANDUSE         = MMINLU
    NoahmpIO%IOPT_CROP        = IOPT_CROP
    NoahmpIO%IOPT_IRR         = IOPT_IRR
    NoahmpIO%IOPT_IRRM        = IOPT_IRRM
    NoahmpIO%SF_URBAN_PHYSICS = SF_URBAN_PHYSICS
    NoahmpIO%IOPT_SOIL        = IOPT_SOIL
    NoahmpIO%IOPT_ALB         = IOPT_ALB
    NoahmpIO%IOPT_WETLAND     = IOPT_WETLAND
    NoahmpIO%IOPT_RUNSUB      = IOPT_RUNSUB
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%SNICAR_SNOWOPTICS_OPT = SNICAR_SNOWOPTICS_OPT
       NoahmpIO%SNICAR_DUSTOPTICS_OPT = SNICAR_DUSTOPTICS_OPT
       NoahmpIO%SNICAR_SOLARSPEC_OPT  = SNICAR_SOLARSPEC_OPT
       NoahmpIO%SNICAR_BANDNUMBER_OPT = SNICAR_BANDNUMBER_OPT
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

    !--------- WRF variables mapped to NoahmpIO variables
    ! input only
    NoahmpIO%restart_flag       = restart
    NoahmpIO%IVGTYP             = IVGTYP
    NoahmpIO%ISLTYP             = ISLTYP
    NoahmpIO%FNDSNOWH           = FNDSNOWH
    NoahmpIO%XLAT               = XLAT    
    NoahmpIO%TSK                = TSK
    NoahmpIO%DZS                = DZS
    NoahmpIO%XICE               = XICE    
    NoahmpIO%CROPTYPE           = CROPTYPE
    NoahmpIO%DTBL               = DT
    NoahmpIO%WTDDT              = WTDDT
    NoahmpIO%DX                 = DX
    NoahmpIO%DY                 = DY
    NoahmpIO%FDEPTHXY           = FDEPTHXY
    NoahmpIO%MSFTX              = MSFTX
    NoahmpIO%MSFTY              = MSFTY
    NoahmpIO%TERRAIN            = HT
    NoahmpIO%RECHCLIM           = RECHCLIM

    ! in/out variables
    NoahmpIO%SMOIS              = SMOIS
    NoahmpIO%SH2O               = SH2O
    NoahmpIO%TSLB               = TSLB
    NoahmpIO%SNOW               = SNOW    
    NoahmpIO%SNOWH              = SNOWH   
    NoahmpIO%CANWAT             = CANWAT  
    NoahmpIO%CANICEXY           = CANICEXY
    NoahmpIO%CANLIQXY           = CANLIQXY
    NoahmpIO%TMN                = TMN
    NoahmpIO%ISNOWXY            = ISNOWXY
    NoahmpIO%ZSNSOXY            = ZSNSOXY 
    NoahmpIO%TSNOXY             = TSNOXY  
    NoahmpIO%SNICEXY            = SNICEXY 
    NoahmpIO%SNLIQXY            = SNLIQXY 
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
    NoahmpIO%LFMASSXY           = LFMASSXY
    NoahmpIO%RTMASSXY           = RTMASSXY
    NoahmpIO%STMASSXY           = STMASSXY
    NoahmpIO%WOODXY             = WOODXY
    NoahmpIO%GRAINXY            = GRAINXY
    NoahmpIO%GDDXY              = GDDXY
    NoahmpIO%STBLCPXY           = STBLCPXY
    NoahmpIO%FASTCPXY           = FASTCPXY
    NoahmpIO%LAI                = LAI
    NoahmpIO%XSAIXY             = XSAIXY
    NoahmpIO%QTDRAIN            = QTDRAIN
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
    NoahmpIO%T2MVXY             = T2MVXY
    NoahmpIO%T2MBXY             = T2MBXY
    NoahmpIO%SMCWTDXY           = SMCWTDXY
    NoahmpIO%DEEPRECHXY         = DEEPRECHXY
    NoahmpIO%RECHXY             = RECHXY
    NoahmpIO%QRFSXY             = QRFSXY
    NoahmpIO%QSPRINGSXY         = QSPRINGSXY
    NoahmpIO%QSLATXY            = QSLATXY
    NoahmpIO%AREAXY             = AREAXY
    NoahmpIO%RIVERBEDXY         = RIVERBEDXY
    NoahmpIO%EQZWT              = EQZWT
    NoahmpIO%RIVERCONDXY        = RIVERCONDXY
    NoahmpIO%PEXPXY             = PEXPXY
    NoahmpIO%SMOISEQ            = SMOISEQ
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

    !--------- WRF -> NoahmpIO variables mapping ends

    !--------- main Noahmp initialization module
    call NoahmpInitMain(NoahmpIO)
    !---------

    !--------- initialized NoahmpIO variable mapped to WRF variables
    ! in/out variables
    SMOIS        = NoahmpIO%SMOIS 
    SH2O         = NoahmpIO%SH2O
    TSLB         = NoahmpIO%TSLB
    SNOW         = NoahmpIO%SNOW
    SNOWH        = NoahmpIO%SNOWH 
    CANWAT       = NoahmpIO%CANWAT
    CANICEXY     = NoahmpIO%CANICEXY
    CANLIQXY     = NoahmpIO%CANLIQXY
    TMN          = NoahmpIO%TMN
    ISNOWXY      = NoahmpIO%ISNOWXY
    ZSNSOXY      = NoahmpIO%ZSNSOXY
    TSNOXY       = NoahmpIO%TSNOXY
    SNICEXY      = NoahmpIO%SNICEXY
    SNLIQXY      = NoahmpIO%SNLIQXY
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
    LFMASSXY     = NoahmpIO%LFMASSXY
    RTMASSXY     = NoahmpIO%RTMASSXY
    STMASSXY     = NoahmpIO%STMASSXY
    WOODXY       = NoahmpIO%WOODXY
    GRAINXY      = NoahmpIO%GRAINXY
    GDDXY        = NoahmpIO%GDDXY 
    STBLCPXY     = NoahmpIO%STBLCPXY
    FASTCPXY     = NoahmpIO%FASTCPXY
    LAI          = NoahmpIO%LAI 
    XSAIXY       = NoahmpIO%XSAIXY
    QTDRAIN      = NoahmpIO%QTDRAIN
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
    T2MVXY       = NoahmpIO%T2MVXY
    T2MBXY       = NoahmpIO%T2MBXY
    SMCWTDXY     = NoahmpIO%SMCWTDXY
    DEEPRECHXY   = NoahmpIO%DEEPRECHXY
    RECHXY       = NoahmpIO%RECHXY
    QRFSXY       = NoahmpIO%QRFSXY
    QSPRINGSXY   = NoahmpIO%QSPRINGSXY
    QSLATXY      = NoahmpIO%QSLATXY
    AREAXY       = NoahmpIO%AREAXY 
    RIVERBEDXY   = NoahmpIO%RIVERBEDXY
    EQZWT        = NoahmpIO%EQZWT
    RIVERCONDXY  = NoahmpIO%RIVERCONDXY
    PEXPXY       = NoahmpIO%PEXPXY
    SMOISEQ      = NoahmpIO%SMOISEQ
    CHSTARXY     = 0.1 ! dummy
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

    ! out variables only
    CROPCAT      = NoahmpIO%CROPCAT
    STEPWTD      = NoahmpIO%STEPWTD

    !--------- NoahmpIO -> WRF variables mapping ends

  end subroutine NoahmpWRFinit

end module NoahmpWRFinitMod
