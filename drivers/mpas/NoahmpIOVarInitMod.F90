module NoahmpIOVarInitMod

!!! Initialize Noah-MP input/output variables
!!! Input/Output variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Sep 13, 2026: NoahmpIO%xx change to 1-D vector for MPAS, Cenlin He (NCAR)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type

  implicit none

contains

!=== initialize with default values

  subroutine NoahmpIOVarInitDefault(NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
   
! ------------------------------------------------- 
    associate(                               &
              ITS     =>  NoahmpIO%ITS      ,&
              ITE     =>  NoahmpIO%ITE      ,&
              KTS     =>  NoahmpIO%KTS      ,&
              KTE     =>  NoahmpIO%KTE      ,&
              NSOIL   =>  NoahmpIO%NSOIL    ,&
              NSNOW   =>  NoahmpIO%NSNOW    ,&
              NUMRAD  =>  NoahmpIO%NUMRAD    &
             )
! -------------------------------------------------

    ! Input variables
    if ( .not. allocated (NoahmpIO%COSZEN)    ) allocate ( NoahmpIO%COSZEN     (ITS:ITE         ) ) ! cosine zenith angle
    if ( .not. allocated (NoahmpIO%XLAT)      ) allocate ( NoahmpIO%XLAT       (ITS:ITE         ) ) ! latitude [radians] 
    if ( .not. allocated (NoahmpIO%DX)        ) allocate ( NoahmpIO%DX         (ITS:ITE         ) ) ! grid spacing [m] 
    if ( .not. allocated (NoahmpIO%DY)        ) allocate ( NoahmpIO%DY         (ITS:ITE         ) ) ! grid spacing [m]  
    if ( .not. allocated (NoahmpIO%DZS)       ) allocate ( NoahmpIO%DZS        (1:NSOIL         ) ) ! thickness of soil layers [m]
    if ( .not. allocated (NoahmpIO%ZSOIL)     ) allocate ( NoahmpIO%ZSOIL      (1:NSOIL         ) ) ! depth to soil interfaces [m] 
    if ( .not. allocated (NoahmpIO%IVGTYP)    ) allocate ( NoahmpIO%IVGTYP     (ITS:ITE         ) ) ! vegetation type
    if ( .not. allocated (NoahmpIO%ISLTYP)    ) allocate ( NoahmpIO%ISLTYP     (ITS:ITE         ) ) ! soil type
    if ( .not. allocated (NoahmpIO%VEGFRA)    ) allocate ( NoahmpIO%VEGFRA     (ITS:ITE         ) ) ! vegetation fraction []
    if ( .not. allocated (NoahmpIO%TMN)       ) allocate ( NoahmpIO%TMN        (ITS:ITE         ) ) ! deep soil temperature [K]
    if ( .not. allocated (NoahmpIO%XLAND)     ) allocate ( NoahmpIO%XLAND      (ITS:ITE         ) ) ! =2 ocean; =1 land/seaice
    if ( .not. allocated (NoahmpIO%XICE)      ) allocate ( NoahmpIO%XICE       (ITS:ITE         ) ) ! fraction of grid that is seaice
    if ( .not. allocated (NoahmpIO%SWDOWN)    ) allocate ( NoahmpIO%SWDOWN     (ITS:ITE         ) ) ! solar down at surface [W m-2]
    if ( .not. allocated (NoahmpIO%SWDDIR)    ) allocate ( NoahmpIO%SWDDIR     (ITS:ITE         ) ) ! solar down at surface [W m-2] for new urban solar panel
    if ( .not. allocated (NoahmpIO%SWDDIF)    ) allocate ( NoahmpIO%SWDDIF     (ITS:ITE         ) ) ! solar down at surface [W m-2] for new urban solar panel
    if ( .not. allocated (NoahmpIO%GLW)       ) allocate ( NoahmpIO%GLW        (ITS:ITE         ) ) ! longwave down at surface [W m-2]
    if ( .not. allocated (NoahmpIO%RAINBL)    ) allocate ( NoahmpIO%RAINBL     (ITS:ITE         ) ) ! total precipitation entering land model [mm] per time step
    if ( .not. allocated (NoahmpIO%SNOWBL)    ) allocate ( NoahmpIO%SNOWBL     (ITS:ITE         ) ) ! snow entering land model [mm] per time step
    if ( .not. allocated (NoahmpIO%SR)        ) allocate ( NoahmpIO%SR         (ITS:ITE         ) ) ! frozen precip ratio entering land model [-]
    if ( .not. allocated (NoahmpIO%RAINCV)    ) allocate ( NoahmpIO%RAINCV     (ITS:ITE         ) ) ! convective precip forcing [mm]
    if ( .not. allocated (NoahmpIO%RAINNCV)   ) allocate ( NoahmpIO%RAINNCV    (ITS:ITE         ) ) ! non-convective precip forcing [mm]
    if ( .not. allocated (NoahmpIO%RAINSHV)   ) allocate ( NoahmpIO%RAINSHV    (ITS:ITE         ) ) ! shallow conv. precip forcing [mm]
    if ( .not. allocated (NoahmpIO%SNOWNCV)   ) allocate ( NoahmpIO%SNOWNCV    (ITS:ITE         ) ) ! non-covective snow forcing (subset of rainncv) [mm]
    if ( .not. allocated (NoahmpIO%GRAUPELNCV)) allocate ( NoahmpIO%GRAUPELNCV (ITS:ITE         ) ) ! non-convective graupel forcing (subset of rainncv) [mm]
    if ( .not. allocated (NoahmpIO%HAILNCV)   ) allocate ( NoahmpIO%HAILNCV    (ITS:ITE         ) ) ! non-convective hail forcing (subset of rainncv) [mm]
    if ( .not. allocated (NoahmpIO%MP_RAINC)  ) allocate ( NoahmpIO%MP_RAINC   (ITS:ITE         ) ) ! convective precip forcing [mm]
    if ( .not. allocated (NoahmpIO%MP_RAINNC) ) allocate ( NoahmpIO%MP_RAINNC  (ITS:ITE         ) ) ! non-convective precip forcing [mm]
    if ( .not. allocated (NoahmpIO%MP_SHCV)   ) allocate ( NoahmpIO%MP_SHCV    (ITS:ITE         ) ) ! shallow conv. precip forcing [mm]
    if ( .not. allocated (NoahmpIO%MP_SNOW)   ) allocate ( NoahmpIO%MP_SNOW    (ITS:ITE         ) ) ! non-covective snow (subset of rainnc) [mm]
    if ( .not. allocated (NoahmpIO%MP_GRAUP)  ) allocate ( NoahmpIO%MP_GRAUP   (ITS:ITE         ) ) ! non-convective graupel (subset of rainnc) [mm]
    if ( .not. allocated (NoahmpIO%MP_HAIL)   ) allocate ( NoahmpIO%MP_HAIL    (ITS:ITE         ) ) ! non-convective hail (subset of rainnc) [mm]
    if ( .not. allocated (NoahmpIO%SEAICE)    ) allocate ( NoahmpIO%SEAICE     (ITS:ITE         ) ) ! seaice fraction
    if ( .not. allocated (NoahmpIO%DZ8W)      ) allocate ( NoahmpIO%DZ8W       (ITS:ITE,KTS:KTE ) ) ! thickness of atmo layers [m]
    if ( .not. allocated (NoahmpIO%T_PHY)     ) allocate ( NoahmpIO%T_PHY      (ITS:ITE,KTS:KTE ) ) ! 3D atmospheric temperature valid at mid-levels [K]
    if ( .not. allocated (NoahmpIO%QV_CURR)   ) allocate ( NoahmpIO%QV_CURR    (ITS:ITE,KTS:KTE ) ) ! 3D water vapor mixing ratio [kg/kg_dry]
    if ( .not. allocated (NoahmpIO%U_PHY)     ) allocate ( NoahmpIO%U_PHY      (ITS:ITE,KTS:KTE ) ) ! 3D U wind component [m/s]
    if ( .not. allocated (NoahmpIO%V_PHY)     ) allocate ( NoahmpIO%V_PHY      (ITS:ITE,KTS:KTE ) ) ! 3D V wind component [m/s]
    if ( .not. allocated (NoahmpIO%P8W)       ) allocate ( NoahmpIO%P8W        (ITS:ITE,KTS:KTE ) ) ! 3D pressure, valid at interface [Pa]
 
    ! spatial varying parameter map
    if ( NoahmpIO%IOPT_SOIL > 1 ) then
       if ( .not. allocated (NoahmpIO%soilcomp)) allocate ( NoahmpIO%soilcomp (ITS:ITE,1:2*NSOIL) ) ! Soil sand and clay content [fraction]
       if ( .not. allocated (NoahmpIO%soilcl1) ) allocate ( NoahmpIO%soilcl1  (ITS:ITE          ) ) ! Soil texture class with depth
       if ( .not. allocated (NoahmpIO%soilcl2) ) allocate ( NoahmpIO%soilcl2  (ITS:ITE          ) ) ! Soil texture class with depth
       if ( .not. allocated (NoahmpIO%soilcl3) ) allocate ( NoahmpIO%soilcl3  (ITS:ITE          ) ) ! Soil texture class with depth
       if ( .not. allocated (NoahmpIO%soilcl4) ) allocate ( NoahmpIO%soilcl4  (ITS:ITE          ) ) ! Soil texture class with depth
    endif
    if ( NoahmpIO%IOPT_SOIL == 4 ) then
       if ( .not. allocated (NoahmpIO%bexp_3d)      ) allocate ( NoahmpIO%bexp_3d       (ITS:ITE,1:NSOIL) ) ! C-H B exponent
       if ( .not. allocated (NoahmpIO%smcdry_3D)    ) allocate ( NoahmpIO%smcdry_3D     (ITS:ITE,1:NSOIL) ) ! Soil Moisture Limit: Dry
       if ( .not. allocated (NoahmpIO%smcwlt_3D)    ) allocate ( NoahmpIO%smcwlt_3D     (ITS:ITE,1:NSOIL) ) ! Soil Moisture Limit: Wilt
       if ( .not. allocated (NoahmpIO%smcref_3D)    ) allocate ( NoahmpIO%smcref_3D     (ITS:ITE,1:NSOIL) ) ! Soil Moisture Limit: Reference
       if ( .not. allocated (NoahmpIO%smcmax_3D)    ) allocate ( NoahmpIO%smcmax_3D     (ITS:ITE,1:NSOIL) ) ! Soil Moisture Limit: Max
       if ( .not. allocated (NoahmpIO%dksat_3D)     ) allocate ( NoahmpIO%dksat_3D      (ITS:ITE,1:NSOIL) ) ! Saturated Soil Conductivity
       if ( .not. allocated (NoahmpIO%dwsat_3D)     ) allocate ( NoahmpIO%dwsat_3D      (ITS:ITE,1:NSOIL) ) ! Saturated Soil Diffusivity
       if ( .not. allocated (NoahmpIO%psisat_3D)    ) allocate ( NoahmpIO%psisat_3D     (ITS:ITE,1:NSOIL) ) ! Saturated Matric Potential
       if ( .not. allocated (NoahmpIO%quartz_3D)    ) allocate ( NoahmpIO%quartz_3D     (ITS:ITE,1:NSOIL) ) ! Soil quartz content
       if ( .not. allocated (NoahmpIO%refdk_2D)     ) allocate ( NoahmpIO%refdk_2D      (ITS:ITE        ) ) ! Reference Soil Conductivity
       if ( .not. allocated (NoahmpIO%refkdt_2D)    ) allocate ( NoahmpIO%refkdt_2D     (ITS:ITE        ) ) ! Soil Infiltration Parameter
       if ( .not. allocated (NoahmpIO%irr_frac_2D)  ) allocate ( NoahmpIO%irr_frac_2D   (ITS:ITE        ) ) ! irrigation Fraction
       if ( .not. allocated (NoahmpIO%irr_har_2D)   ) allocate ( NoahmpIO%irr_har_2D    (ITS:ITE        ) ) ! number of days before harvest date to stop irrigation 
       if ( .not. allocated (NoahmpIO%irr_lai_2D)   ) allocate ( NoahmpIO%irr_lai_2D    (ITS:ITE        ) ) ! Minimum lai to trigger irrigation
       if ( .not. allocated (NoahmpIO%irr_mad_2D)   ) allocate ( NoahmpIO%irr_mad_2D    (ITS:ITE        ) ) ! management allowable deficit (0-1)
       if ( .not. allocated (NoahmpIO%filoss_2D)    ) allocate ( NoahmpIO%filoss_2D     (ITS:ITE        ) ) ! fraction of flood irrigation loss (0-1) 
       if ( .not. allocated (NoahmpIO%sprir_rate_2D)) allocate ( NoahmpIO%sprir_rate_2D (ITS:ITE        ) ) ! mm/h, sprinkler irrigation rate
       if ( .not. allocated (NoahmpIO%micir_rate_2D)) allocate ( NoahmpIO%micir_rate_2D (ITS:ITE        ) ) ! mm/h, micro irrigation rate
       if ( .not. allocated (NoahmpIO%firtfac_2D)   ) allocate ( NoahmpIO%firtfac_2D    (ITS:ITE        ) ) ! flood application rate factor
       if ( .not. allocated (NoahmpIO%ir_rain_2D)   ) allocate ( NoahmpIO%ir_rain_2D    (ITS:ITE        ) ) ! maximum precipitation to stop irrigation trigger
       if ( .not. allocated (NoahmpIO%bvic_2D)      ) allocate ( NoahmpIO%bvic_2D       (ITS:ITE        ) ) ! VIC model infiltration parameter [-]
       if ( .not. allocated (NoahmpIO%axaj_2D)      ) allocate ( NoahmpIO%axaj_2D       (ITS:ITE        ) ) ! Tension water distribution inflection parameter [-]
       if ( .not. allocated (NoahmpIO%bxaj_2D)      ) allocate ( NoahmpIO%bxaj_2D       (ITS:ITE        ) ) ! Tension water distribution shape parameter [-]
       if ( .not. allocated (NoahmpIO%xxaj_2D)      ) allocate ( NoahmpIO%xxaj_2D       (ITS:ITE        ) ) ! Free water distribution shape parameter [-]
       if ( .not. allocated (NoahmpIO%bdvic_2D)     ) allocate ( NoahmpIO%bdvic_2D      (ITS:ITE        ) ) ! DVIC model infiltration parameter [-]
       if ( .not. allocated (NoahmpIO%gdvic_2D)     ) allocate ( NoahmpIO%gdvic_2D      (ITS:ITE        ) ) ! Mean Capillary Drive (m) for infiltration models
       if ( .not. allocated (NoahmpIO%bbvic_2D)     ) allocate ( NoahmpIO%bbvic_2D      (ITS:ITE        ) ) ! DVIC heterogeniety parameter for infiltration [-]
       if ( .not. allocated (NoahmpIO%KLAT_FAC)     ) allocate ( NoahmpIO%KLAT_FAC      (ITS:ITE        ) ) ! factor multiplier to hydraulic conductivity
       if ( .not. allocated (NoahmpIO%TDSMC_FAC)    ) allocate ( NoahmpIO%TDSMC_FAC     (ITS:ITE        ) ) ! factor multiplier to field capacity
       if ( .not. allocated (NoahmpIO%TD_DC)        ) allocate ( NoahmpIO%TD_DC         (ITS:ITE        ) ) ! drainage coefficient for simple
       if ( .not. allocated (NoahmpIO%TD_DCOEF)     ) allocate ( NoahmpIO%TD_DCOEF      (ITS:ITE        ) ) ! drainge coefficient for Hooghoudt 
       if ( .not. allocated (NoahmpIO%TD_DDRAIN)    ) allocate ( NoahmpIO%TD_DDRAIN     (ITS:ITE        ) ) ! depth of drain
       if ( .not. allocated (NoahmpIO%TD_RADI)      ) allocate ( NoahmpIO%TD_RADI       (ITS:ITE        ) ) ! tile radius
       if ( .not. allocated (NoahmpIO%TD_SPAC)      ) allocate ( NoahmpIO%TD_SPAC       (ITS:ITE        ) ) ! tile spacing
    endif

    ! INOUT (with generic LSM equivalent) (as defined in WRF)
    if ( .not. allocated (NoahmpIO%TSK)      ) allocate ( NoahmpIO%TSK       (ITS:ITE        ) ) ! surface radiative temperature [K]
    if ( .not. allocated (NoahmpIO%HFX)      ) allocate ( NoahmpIO%HFX       (ITS:ITE        ) ) ! sensible heat flux [W m-2]
    if ( .not. allocated (NoahmpIO%QFX)      ) allocate ( NoahmpIO%QFX       (ITS:ITE        ) ) ! latent heat flux [kg s-1 m-2]
    if ( .not. allocated (NoahmpIO%LH)       ) allocate ( NoahmpIO%LH        (ITS:ITE        ) ) ! latent heat flux [W m-2]
    if ( .not. allocated (NoahmpIO%GRDFLX)   ) allocate ( NoahmpIO%GRDFLX    (ITS:ITE        ) ) ! ground/snow heat flux [W m-2]
    if ( .not. allocated (NoahmpIO%SMSTAV)   ) allocate ( NoahmpIO%SMSTAV    (ITS:ITE        ) ) ! soil moisture avail. [not used]
    if ( .not. allocated (NoahmpIO%SMSTOT)   ) allocate ( NoahmpIO%SMSTOT    (ITS:ITE        ) ) ! total soil water [mm][not used]
    if ( .not. allocated (NoahmpIO%SFCRUNOFF)) allocate ( NoahmpIO%SFCRUNOFF (ITS:ITE        ) ) ! accumulated surface runoff [m]
    if ( .not. allocated (NoahmpIO%UDRUNOFF) ) allocate ( NoahmpIO%UDRUNOFF  (ITS:ITE        ) ) ! accumulated sub-surface runoff [m]
    if ( .not. allocated (NoahmpIO%ALBEDO)   ) allocate ( NoahmpIO%ALBEDO    (ITS:ITE        ) ) ! total grid albedo []
    if ( .not. allocated (NoahmpIO%SNOWC)    ) allocate ( NoahmpIO%SNOWC     (ITS:ITE        ) ) ! snow cover fraction []
    if ( .not. allocated (NoahmpIO%SNOW)     ) allocate ( NoahmpIO%SNOW      (ITS:ITE        ) ) ! snow water equivalent [mm]
    if ( .not. allocated (NoahmpIO%SNOWH)    ) allocate ( NoahmpIO%SNOWH     (ITS:ITE        ) ) ! physical snow depth [m]
    if ( .not. allocated (NoahmpIO%CANWAT)   ) allocate ( NoahmpIO%CANWAT    (ITS:ITE        ) ) ! total canopy water + ice [mm]
    if ( .not. allocated (NoahmpIO%ACSNOM)   ) allocate ( NoahmpIO%ACSNOM    (ITS:ITE        ) ) ! accumulated snow melt leaving pack
    if ( .not. allocated (NoahmpIO%ACSNOW)   ) allocate ( NoahmpIO%ACSNOW    (ITS:ITE        ) ) ! accumulated snow on grid
    if ( .not. allocated (NoahmpIO%EMISS)    ) allocate ( NoahmpIO%EMISS     (ITS:ITE        ) ) ! surface bulk emissivity
    if ( .not. allocated (NoahmpIO%QSFC)     ) allocate ( NoahmpIO%QSFC      (ITS:ITE        ) ) ! bulk surface specific humidity
    if ( .not. allocated (NoahmpIO%SMOISEQ)  ) allocate ( NoahmpIO%SMOISEQ   (ITS:ITE,1:NSOIL) ) ! equilibrium volumetric soil moisture [m3/m3]
    if ( .not. allocated (NoahmpIO%SMOIS)    ) allocate ( NoahmpIO%SMOIS     (ITS:ITE,1:NSOIL) ) ! volumetric soil moisture [m3/m3]
    if ( .not. allocated (NoahmpIO%SH2O)     ) allocate ( NoahmpIO%SH2O      (ITS:ITE,1:NSOIL) ) ! volumetric liquid soil moisture [m3/m3]
    if ( .not. allocated (NoahmpIO%TSLB)     ) allocate ( NoahmpIO%TSLB      (ITS:ITE,1:NSOIL) ) ! soil temperature [K]

    ! INOUT (with no Noah LSM equivalent) (as defined in WRF)
    if ( .not. allocated (NoahmpIO%ISNOWXY)   ) allocate ( NoahmpIO%ISNOWXY    (ITS:ITE               ) ) ! actual no. of snow layers
    if ( .not. allocated (NoahmpIO%TVXY)      ) allocate ( NoahmpIO%TVXY       (ITS:ITE               ) ) ! vegetation leaf temperature
    if ( .not. allocated (NoahmpIO%TGXY)      ) allocate ( NoahmpIO%TGXY       (ITS:ITE               ) ) ! bulk ground surface temperature
    if ( .not. allocated (NoahmpIO%CANICEXY)  ) allocate ( NoahmpIO%CANICEXY   (ITS:ITE               ) ) ! canopy-intercepted ice (mm)
    if ( .not. allocated (NoahmpIO%CANLIQXY)  ) allocate ( NoahmpIO%CANLIQXY   (ITS:ITE               ) ) ! canopy-intercepted liquid water (mm)
    if ( .not. allocated (NoahmpIO%EAHXY)     ) allocate ( NoahmpIO%EAHXY      (ITS:ITE               ) ) ! canopy air vapor pressure (pa)
    if ( .not. allocated (NoahmpIO%TAHXY)     ) allocate ( NoahmpIO%TAHXY      (ITS:ITE               ) ) ! canopy air temperature (k)
    if ( .not. allocated (NoahmpIO%CMXY)      ) allocate ( NoahmpIO%CMXY       (ITS:ITE               ) ) ! bulk momentum drag coefficient
    if ( .not. allocated (NoahmpIO%CHXY)      ) allocate ( NoahmpIO%CHXY       (ITS:ITE               ) ) ! bulk sensible heat exchange coefficient
    if ( .not. allocated (NoahmpIO%FWETXY)    ) allocate ( NoahmpIO%FWETXY     (ITS:ITE               ) ) ! wetted or snowed fraction of the canopy (-)
    if ( .not. allocated (NoahmpIO%SNEQVOXY)  ) allocate ( NoahmpIO%SNEQVOXY   (ITS:ITE               ) ) ! snow mass at last time step(mm h2o)
    if ( .not. allocated (NoahmpIO%ALBOLDXY)  ) allocate ( NoahmpIO%ALBOLDXY   (ITS:ITE               ) ) ! snow albedo at last time step (-)
    if ( .not. allocated (NoahmpIO%QSNOWXY)   ) allocate ( NoahmpIO%QSNOWXY    (ITS:ITE               ) ) ! snowfall on the ground [mm/s]
    if ( .not. allocated (NoahmpIO%QRAINXY)   ) allocate ( NoahmpIO%QRAINXY    (ITS:ITE               ) ) ! rainfall on the ground [mm/s]
    if ( .not. allocated (NoahmpIO%WSLAKEXY)  ) allocate ( NoahmpIO%WSLAKEXY   (ITS:ITE               ) ) ! lake water storage [mm]
    if ( .not. allocated (NoahmpIO%ZWTXY)     ) allocate ( NoahmpIO%ZWTXY      (ITS:ITE               ) ) ! water table depth [m]
    if ( .not. allocated (NoahmpIO%WAXY)      ) allocate ( NoahmpIO%WAXY       (ITS:ITE               ) ) ! water in the "aquifer" [mm]
    if ( .not. allocated (NoahmpIO%WTXY)      ) allocate ( NoahmpIO%WTXY       (ITS:ITE               ) ) ! groundwater storage [mm]
    if ( .not. allocated (NoahmpIO%SMCWTDXY)  ) allocate ( NoahmpIO%SMCWTDXY   (ITS:ITE               ) ) ! soil moisture below the bottom of the column (m3m-3)
    if ( .not. allocated (NoahmpIO%DEEPRECHXY)) allocate ( NoahmpIO%DEEPRECHXY (ITS:ITE               ) ) ! recharge to the water table when deep (m)
    if ( .not. allocated (NoahmpIO%RECHXY)    ) allocate ( NoahmpIO%RECHXY     (ITS:ITE               ) ) ! recharge to the water table (diagnostic) (mm)
    if ( .not. allocated (NoahmpIO%LFMASSXY)  ) allocate ( NoahmpIO%LFMASSXY   (ITS:ITE               ) ) ! leaf mass [g/m2]
    if ( .not. allocated (NoahmpIO%RTMASSXY)  ) allocate ( NoahmpIO%RTMASSXY   (ITS:ITE               ) ) ! mass of fine roots [g/m2]
    if ( .not. allocated (NoahmpIO%STMASSXY)  ) allocate ( NoahmpIO%STMASSXY   (ITS:ITE               ) ) ! stem mass [g/m2]
    if ( .not. allocated (NoahmpIO%WOODXY)    ) allocate ( NoahmpIO%WOODXY     (ITS:ITE               ) ) ! mass of wood (incl. woody roots) [g/m2]
    if ( .not. allocated (NoahmpIO%GRAINXY)   ) allocate ( NoahmpIO%GRAINXY    (ITS:ITE               ) ) ! mass of grain XING [g/m2]
    if ( .not. allocated (NoahmpIO%GDDXY)     ) allocate ( NoahmpIO%GDDXY      (ITS:ITE               ) ) ! growing degree days XING FOUR
    if ( .not. allocated (NoahmpIO%STBLCPXY)  ) allocate ( NoahmpIO%STBLCPXY   (ITS:ITE               ) ) ! stable carbon in deep soil [g/m2]
    if ( .not. allocated (NoahmpIO%FASTCPXY)  ) allocate ( NoahmpIO%FASTCPXY   (ITS:ITE               ) ) ! short-lived carbon, shallow soil [g/m2]
    if ( .not. allocated (NoahmpIO%LAI)       ) allocate ( NoahmpIO%LAI        (ITS:ITE               ) ) ! leaf area index
    if ( .not. allocated (NoahmpIO%XSAIXY)    ) allocate ( NoahmpIO%XSAIXY     (ITS:ITE               ) ) ! stem area index
    if ( .not. allocated (NoahmpIO%TAUSSXY)   ) allocate ( NoahmpIO%TAUSSXY    (ITS:ITE               ) ) ! snow age factor
    if ( .not. allocated (NoahmpIO%TSNOXY)    ) allocate ( NoahmpIO%TSNOXY     (ITS:ITE,-NSNOW+1:0    ) ) ! snow temperature [K]
    if ( .not. allocated (NoahmpIO%ZSNSOXY)   ) allocate ( NoahmpIO%ZSNSOXY    (ITS:ITE,-NSNOW+1:NSOIL) ) ! snow layer depth [m]
    if ( .not. allocated (NoahmpIO%SNICEXY)   ) allocate ( NoahmpIO%SNICEXY    (ITS:ITE,-NSNOW+1:0    ) ) ! snow layer ice [mm]
    if ( .not. allocated (NoahmpIO%SNLIQXY)   ) allocate ( NoahmpIO%SNLIQXY    (ITS:ITE,-NSNOW+1:0    ) ) ! snow layer liquid water [mm]

    ! irrigation
    if ( .not. allocated (NoahmpIO%IRFRACT) ) allocate ( NoahmpIO%IRFRACT (ITS:ITE) ) ! irrigation fraction
    if ( .not. allocated (NoahmpIO%SIFRACT) ) allocate ( NoahmpIO%SIFRACT (ITS:ITE) ) ! sprinkler irrigation fraction
    if ( .not. allocated (NoahmpIO%MIFRACT) ) allocate ( NoahmpIO%MIFRACT (ITS:ITE) ) ! micro irrigation fraction
    if ( .not. allocated (NoahmpIO%FIFRACT) ) allocate ( NoahmpIO%FIFRACT (ITS:ITE) ) ! flood irrigation fraction   
    if ( .not. allocated (NoahmpIO%IRNUMSI) ) allocate ( NoahmpIO%IRNUMSI (ITS:ITE) ) ! irrigation event number, Sprinkler
    if ( .not. allocated (NoahmpIO%IRNUMMI) ) allocate ( NoahmpIO%IRNUMMI (ITS:ITE) ) ! irrigation event number, Micro
    if ( .not. allocated (NoahmpIO%IRNUMFI) ) allocate ( NoahmpIO%IRNUMFI (ITS:ITE) ) ! irrigation event number, Flood 
    if ( .not. allocated (NoahmpIO%IRWATSI) ) allocate ( NoahmpIO%IRWATSI (ITS:ITE) ) ! irrigation water amount [m] to be applied, Sprinkler
    if ( .not. allocated (NoahmpIO%IRWATMI) ) allocate ( NoahmpIO%IRWATMI (ITS:ITE) ) ! irrigation water amount [m] to be applied, Micro
    if ( .not. allocated (NoahmpIO%IRWATFI) ) allocate ( NoahmpIO%IRWATFI (ITS:ITE) ) ! irrigation water amount [m] to be applied, Flood
    if ( .not. allocated (NoahmpIO%IRELOSS) ) allocate ( NoahmpIO%IRELOSS (ITS:ITE) ) ! loss of irrigation water to evaporation,sprinkler [mm]
    if ( .not. allocated (NoahmpIO%IRSIVOL) ) allocate ( NoahmpIO%IRSIVOL (ITS:ITE) ) ! amount of irrigation by sprinkler (mm)
    if ( .not. allocated (NoahmpIO%IRMIVOL) ) allocate ( NoahmpIO%IRMIVOL (ITS:ITE) ) ! amount of irrigation by micro (mm)
    if ( .not. allocated (NoahmpIO%IRFIVOL) ) allocate ( NoahmpIO%IRFIVOL (ITS:ITE) ) ! amount of irrigation by micro (mm)
    if ( .not. allocated (NoahmpIO%IRRSPLH) ) allocate ( NoahmpIO%IRRSPLH (ITS:ITE) ) ! latent heating from sprinkler evaporation (w/m2)
    if ( .not. allocated (NoahmpIO%LOCTIM)  ) allocate ( NoahmpIO%LOCTIM  (ITS:ITE) ) ! local time
  
    ! OUT (with no Noah LSM equivalent) (as defined in WRF)   
    if ( .not. allocated (NoahmpIO%T2MVXY)     ) allocate ( NoahmpIO%T2MVXY      (ITS:ITE) ) ! 2m temperature of vegetation part [K]
    if ( .not. allocated (NoahmpIO%T2MBXY)     ) allocate ( NoahmpIO%T2MBXY      (ITS:ITE) ) ! 2m temperature of bare ground part [K]
    if ( .not. allocated (NoahmpIO%T2MXY )     ) allocate ( NoahmpIO%T2MXY       (ITS:ITE) ) ! 2m temperature grid mean [K]
    if ( .not. allocated (NoahmpIO%Q2MVXY)     ) allocate ( NoahmpIO%Q2MVXY      (ITS:ITE) ) ! 2m mixing ratio of vegetation part [kg/kg]
    if ( .not. allocated (NoahmpIO%Q2MBXY)     ) allocate ( NoahmpIO%Q2MBXY      (ITS:ITE) ) ! 2m mixing ratio of bare ground part [kg/kg]
    if ( .not. allocated (NoahmpIO%Q2MXY )     ) allocate ( NoahmpIO%Q2MXY       (ITS:ITE) ) ! 2m mixing ratio grid mean [kg/kg]
    if ( .not. allocated (NoahmpIO%TRADXY)     ) allocate ( NoahmpIO%TRADXY      (ITS:ITE) ) ! surface radiative temperature (k)
    if ( .not. allocated (NoahmpIO%NEEXY)      ) allocate ( NoahmpIO%NEEXY       (ITS:ITE) ) ! net ecosys exchange (g/m2/s CO2)
    if ( .not. allocated (NoahmpIO%GPPXY)      ) allocate ( NoahmpIO%GPPXY       (ITS:ITE) ) ! gross primary assimilation [g/m2/s C]
    if ( .not. allocated (NoahmpIO%NPPXY)      ) allocate ( NoahmpIO%NPPXY       (ITS:ITE) ) ! net primary productivity [g/m2/s C]
    if ( .not. allocated (NoahmpIO%FVEGXY)     ) allocate ( NoahmpIO%FVEGXY      (ITS:ITE) ) ! Noah-MP vegetation fraction [-]
    if ( .not. allocated (NoahmpIO%RUNSFXY)    ) allocate ( NoahmpIO%RUNSFXY     (ITS:ITE) ) ! surface runoff [mm per soil timestep]
    if ( .not. allocated (NoahmpIO%RUNSBXY)    ) allocate ( NoahmpIO%RUNSBXY     (ITS:ITE) ) ! subsurface runoff [mm per soil timestep]
    if ( .not. allocated (NoahmpIO%ECANXY)     ) allocate ( NoahmpIO%ECANXY      (ITS:ITE) ) ! evaporation of intercepted water (mm/s)
    if ( .not. allocated (NoahmpIO%EDIRXY)     ) allocate ( NoahmpIO%EDIRXY      (ITS:ITE) ) ! soil surface evaporation rate (mm/s]
    if ( .not. allocated (NoahmpIO%ETRANXY)    ) allocate ( NoahmpIO%ETRANXY     (ITS:ITE) ) ! transpiration rate (mm/s)
    if ( .not. allocated (NoahmpIO%FSAXY)      ) allocate ( NoahmpIO%FSAXY       (ITS:ITE) ) ! total absorbed solar radiation (w/m2)
    if ( .not. allocated (NoahmpIO%FIRAXY)     ) allocate ( NoahmpIO%FIRAXY      (ITS:ITE) ) ! total net longwave rad (w/m2) [+ to atm]
    if ( .not. allocated (NoahmpIO%APARXY)     ) allocate ( NoahmpIO%APARXY      (ITS:ITE) ) ! photosyn active energy by canopy (w/m2)
    if ( .not. allocated (NoahmpIO%PSNXY)      ) allocate ( NoahmpIO%PSNXY       (ITS:ITE) ) ! total photosynthesis (umol co2/m2/s) [+]
    if ( .not. allocated (NoahmpIO%SAVXY)      ) allocate ( NoahmpIO%SAVXY       (ITS:ITE) ) ! solar rad absorbed by veg. (w/m2)
    if ( .not. allocated (NoahmpIO%SAGXY)      ) allocate ( NoahmpIO%SAGXY       (ITS:ITE) ) ! solar rad absorbed by ground (w/m2)
    if ( .not. allocated (NoahmpIO%RSSUNXY)    ) allocate ( NoahmpIO%RSSUNXY     (ITS:ITE) ) ! sunlit leaf stomatal resistance (s/m)
    if ( .not. allocated (NoahmpIO%RSSHAXY)    ) allocate ( NoahmpIO%RSSHAXY     (ITS:ITE) ) ! shaded leaf stomatal resistance (s/m)
    if ( .not. allocated (NoahmpIO%BGAPXY)     ) allocate ( NoahmpIO%BGAPXY      (ITS:ITE) ) ! between gap fraction
    if ( .not. allocated (NoahmpIO%WGAPXY)     ) allocate ( NoahmpIO%WGAPXY      (ITS:ITE) ) ! within gap fraction
    if ( .not. allocated (NoahmpIO%TGVXY)      ) allocate ( NoahmpIO%TGVXY       (ITS:ITE) ) ! under canopy ground temperature[K]
    if ( .not. allocated (NoahmpIO%TGBXY)      ) allocate ( NoahmpIO%TGBXY       (ITS:ITE) ) ! bare ground temperature [K]
    if ( .not. allocated (NoahmpIO%CHVXY)      ) allocate ( NoahmpIO%CHVXY       (ITS:ITE) ) ! sensible heat exchange coefficient vegetated
    if ( .not. allocated (NoahmpIO%CHBXY)      ) allocate ( NoahmpIO%CHBXY       (ITS:ITE) ) ! sensible heat exchange coefficient bare-ground
    if ( .not. allocated (NoahmpIO%SHGXY)      ) allocate ( NoahmpIO%SHGXY       (ITS:ITE) ) ! veg ground sen. heat [w/m2]   [+ to atm]
    if ( .not. allocated (NoahmpIO%SHCXY)      ) allocate ( NoahmpIO%SHCXY       (ITS:ITE) ) ! canopy sen. heat [w/m2]   [+ to atm]
    if ( .not. allocated (NoahmpIO%SHBXY)      ) allocate ( NoahmpIO%SHBXY       (ITS:ITE) ) ! bare sensible heat [w/m2]  [+ to atm]
    if ( .not. allocated (NoahmpIO%EVGXY)      ) allocate ( NoahmpIO%EVGXY       (ITS:ITE) ) ! veg ground evap. heat [w/m2]  [+ to atm]
    if ( .not. allocated (NoahmpIO%EVBXY)      ) allocate ( NoahmpIO%EVBXY       (ITS:ITE) ) ! bare soil evaporation [w/m2]  [+ to atm]
    if ( .not. allocated (NoahmpIO%GHVXY)      ) allocate ( NoahmpIO%GHVXY       (ITS:ITE) ) ! veg ground heat flux [w/m2]  [+ to soil]
    if ( .not. allocated (NoahmpIO%GHBXY)      ) allocate ( NoahmpIO%GHBXY       (ITS:ITE) ) ! bare ground heat flux [w/m2] [+ to soil]
    if ( .not. allocated (NoahmpIO%IRGXY)      ) allocate ( NoahmpIO%IRGXY       (ITS:ITE) ) ! veg ground net LW rad. [w/m2] [+ to atm]
    if ( .not. allocated (NoahmpIO%IRCXY)      ) allocate ( NoahmpIO%IRCXY       (ITS:ITE) ) ! canopy net LW rad. [w/m2] [+ to atm]
    if ( .not. allocated (NoahmpIO%IRBXY)      ) allocate ( NoahmpIO%IRBXY       (ITS:ITE) ) ! bare net longwave rad. [w/m2] [+ to atm]
    if ( .not. allocated (NoahmpIO%TRXY)       ) allocate ( NoahmpIO%TRXY        (ITS:ITE) ) ! transpiration [w/m2]  [+ to atm]
    if ( .not. allocated (NoahmpIO%EVCXY)      ) allocate ( NoahmpIO%EVCXY       (ITS:ITE) ) ! canopy evaporation heat [w/m2]  [+ to atm]
    if ( .not. allocated (NoahmpIO%CHLEAFXY)   ) allocate ( NoahmpIO%CHLEAFXY    (ITS:ITE) ) ! leaf exchange coefficient 
    if ( .not. allocated (NoahmpIO%CHUCXY)     ) allocate ( NoahmpIO%CHUCXY      (ITS:ITE) ) ! under canopy exchange coefficient 
    if ( .not. allocated (NoahmpIO%CHV2XY)     ) allocate ( NoahmpIO%CHV2XY      (ITS:ITE) ) ! veg 2m exchange coefficient 
    if ( .not. allocated (NoahmpIO%CHB2XY)     ) allocate ( NoahmpIO%CHB2XY      (ITS:ITE) ) ! bare 2m exchange coefficient 
    if ( .not. allocated (NoahmpIO%RS)         ) allocate ( NoahmpIO%RS          (ITS:ITE) ) ! Total stomatal resistance (s/m)
    if ( .not. allocated (NoahmpIO%Z0)         ) allocate ( NoahmpIO%Z0          (ITS:ITE) ) ! roughness length output to WRF 
    if ( .not. allocated (NoahmpIO%ZNT)        ) allocate ( NoahmpIO%ZNT         (ITS:ITE) ) ! roughness length output to WRF 
    if ( .not. allocated (NoahmpIO%QTDRAIN)    ) allocate ( NoahmpIO%QTDRAIN     (ITS:ITE) ) ! tile drainage (mm)
    if ( .not. allocated (NoahmpIO%TD_FRACTION)) allocate ( NoahmpIO%TD_FRACTION (ITS:ITE) ) ! tile drainage fraction
    if ( .not. allocated (NoahmpIO%XLONG)      ) allocate ( NoahmpIO%XLONG       (ITS:ITE) ) ! longitude
    if ( .not. allocated (NoahmpIO%TERRAIN)    ) allocate ( NoahmpIO%TERRAIN     (ITS:ITE) ) ! terrain height
    if ( .not. allocated (NoahmpIO%GVFMIN)     ) allocate ( NoahmpIO%GVFMIN      (ITS:ITE) ) ! annual minimum in vegetation fraction
    if ( .not. allocated (NoahmpIO%GVFMAX)     ) allocate ( NoahmpIO%GVFMAX      (ITS:ITE) ) ! annual maximum in vegetation fraction

    ! additional output variables
    if ( .not. allocated (NoahmpIO%PAHXY)       ) allocate ( NoahmpIO%PAHXY        (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%PAHGXY)      ) allocate ( NoahmpIO%PAHGXY       (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%PAHBXY)      ) allocate ( NoahmpIO%PAHBXY       (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%PAHVXY)      ) allocate ( NoahmpIO%PAHVXY       (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QINTSXY)     ) allocate ( NoahmpIO%QINTSXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QINTRXY)     ) allocate ( NoahmpIO%QINTRXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QDRIPSXY)    ) allocate ( NoahmpIO%QDRIPSXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QDRIPRXY)    ) allocate ( NoahmpIO%QDRIPRXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QTHROSXY)    ) allocate ( NoahmpIO%QTHROSXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QTHRORXY)    ) allocate ( NoahmpIO%QTHRORXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QSNSUBXY)    ) allocate ( NoahmpIO%QSNSUBXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QSNFROXY)    ) allocate ( NoahmpIO%QSNFROXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QSUBCXY)     ) allocate ( NoahmpIO%QSUBCXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QFROCXY)     ) allocate ( NoahmpIO%QFROCXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QEVACXY)     ) allocate ( NoahmpIO%QEVACXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QDEWCXY)     ) allocate ( NoahmpIO%QDEWCXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QFRZCXY)     ) allocate ( NoahmpIO%QFRZCXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QMELTCXY)    ) allocate ( NoahmpIO%QMELTCXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QSNBOTXY)    ) allocate ( NoahmpIO%QSNBOTXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QMELTXY)     ) allocate ( NoahmpIO%QMELTXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%PONDINGXY)   ) allocate ( NoahmpIO%PONDINGXY    (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%FPICEXY)     ) allocate ( NoahmpIO%FPICEXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%RAINLSM)     ) allocate ( NoahmpIO%RAINLSM      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%SNOWLSM)     ) allocate ( NoahmpIO%SNOWLSM      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%FORCTLSM)    ) allocate ( NoahmpIO%FORCTLSM     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%FORCQLSM)    ) allocate ( NoahmpIO%FORCQLSM     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%FORCPLSM)    ) allocate ( NoahmpIO%FORCPLSM     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%FORCZLSM)    ) allocate ( NoahmpIO%FORCZLSM     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%FORCWLSM)    ) allocate ( NoahmpIO%FORCWLSM     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%EFLXBXY)     ) allocate ( NoahmpIO%EFLXBXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%SOILENERGY)  ) allocate ( NoahmpIO%SOILENERGY   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%SNOWENERGY)  ) allocate ( NoahmpIO%SNOWENERGY   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%CANHSXY)     ) allocate ( NoahmpIO%CANHSXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_DWATERXY)) allocate ( NoahmpIO%ACC_DWATERXY (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_PRCPXY)  ) allocate ( NoahmpIO%ACC_PRCPXY   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_ECANXY)  ) allocate ( NoahmpIO%ACC_ECANXY   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_ETRANXY) ) allocate ( NoahmpIO%ACC_ETRANXY  (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_EDIRXY)  ) allocate ( NoahmpIO%ACC_EDIRXY   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_SSOILXY) ) allocate ( NoahmpIO%ACC_SSOILXY  (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_QINSURXY)) allocate ( NoahmpIO%ACC_QINSURXY (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_QSEVAXY) ) allocate ( NoahmpIO%ACC_QSEVAXY  (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ACC_ETRANIXY)) allocate ( NoahmpIO%ACC_ETRANIXY (ITS:ITE,1:NSOIL) )
    if ( .not. allocated (NoahmpIO%ACC_GLAFLWXY)) allocate ( NoahmpIO%ACC_GLAFLWXY (ITS:ITE) )

    ! Needed for MMF_RUNOFF (IOPT_RUNSUB = 5); not part of MP driver in WRF
    if ( .not. allocated (NoahmpIO%MSFTX)      ) allocate ( NoahmpIO%MSFTX       (ITS:ITE) ) 
    if ( .not. allocated (NoahmpIO%MSFTY)      ) allocate ( NoahmpIO%MSFTY       (ITS:ITE) ) 
    if ( .not. allocated (NoahmpIO%EQZWT)      ) allocate ( NoahmpIO%EQZWT       (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%RIVERBEDXY) ) allocate ( NoahmpIO%RIVERBEDXY  (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%RIVERCONDXY)) allocate ( NoahmpIO%RIVERCONDXY (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%PEXPXY)     ) allocate ( NoahmpIO%PEXPXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%FDEPTHXY)   ) allocate ( NoahmpIO%FDEPTHXY    (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%AREAXY)     ) allocate ( NoahmpIO%AREAXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QRFSXY)     ) allocate ( NoahmpIO%QRFSXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QSPRINGSXY) ) allocate ( NoahmpIO%QSPRINGSXY  (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QRFXY)      ) allocate ( NoahmpIO%QRFXY       (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QSPRINGXY)  ) allocate ( NoahmpIO%QSPRINGXY   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QSLATXY)    ) allocate ( NoahmpIO%QSLATXY     (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%QLATXY)     ) allocate ( NoahmpIO%QLATXY      (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%RECHCLIM)   ) allocate ( NoahmpIO%RECHCLIM    (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%RIVERMASK)  ) allocate ( NoahmpIO%RIVERMASK   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%NONRIVERXY) ) allocate ( NoahmpIO%NONRIVERXY  (ITS:ITE) )

    ! Needed for SNICAR SNOW ALBEDO (IOPT_ALB = 3)
    if ( NoahmpIO%IOPT_ALB == 3 ) then

       if ( NoahmpIO%SNICAR_BANDNUMBER_OPT == 1 ) then
          NoahmpIO%snicar_numrad_snw = 5
       elseif ( NoahmpIO%SNICAR_BANDNUMBER_OPT == 2 ) then
          NoahmpIO%snicar_numrad_snw = 480
       endif

       if ( .not. allocated (NoahmpIO%ss_alb_snw_drc)      ) allocate ( NoahmpIO%ss_alb_snw_drc      (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_snw_drc)     ) allocate ( NoahmpIO%asm_prm_snw_drc     (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_snw_drc) ) allocate ( NoahmpIO%ext_cff_mss_snw_drc (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_snw_dfs)      ) allocate ( NoahmpIO%ss_alb_snw_dfs      (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_snw_dfs)     ) allocate ( NoahmpIO%asm_prm_snw_dfs     (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )         
       if ( .not. allocated (NoahmpIO%ext_cff_mss_snw_dfs) ) allocate ( NoahmpIO%ext_cff_mss_snw_dfs (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_bc1)      )     allocate ( NoahmpIO%ss_alb_bc1          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_bc1)     )     allocate ( NoahmpIO%asm_prm_bc1         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_bc1) )     allocate ( NoahmpIO%ext_cff_mss_bc1     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_bc2)      )     allocate ( NoahmpIO%ss_alb_bc2          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_bc2)     )     allocate ( NoahmpIO%asm_prm_bc2         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_bc2) )     allocate ( NoahmpIO%ext_cff_mss_bc2     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_oc1)      )     allocate ( NoahmpIO%ss_alb_oc1          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_oc1)     )     allocate ( NoahmpIO%asm_prm_oc1         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_oc1) )     allocate ( NoahmpIO%ext_cff_mss_oc1     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_oc2)      )     allocate ( NoahmpIO%ss_alb_oc2          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_oc2)     )     allocate ( NoahmpIO%asm_prm_oc2         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_oc2) )     allocate ( NoahmpIO%ext_cff_mss_oc2     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst1)     )     allocate ( NoahmpIO%ss_alb_dst1         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst1)    )     allocate ( NoahmpIO%asm_prm_dst1        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst1))     allocate ( NoahmpIO%ext_cff_mss_dst1    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst2)     )     allocate ( NoahmpIO%ss_alb_dst2         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst2)    )     allocate ( NoahmpIO%asm_prm_dst2        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst2))     allocate ( NoahmpIO%ext_cff_mss_dst2    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst3)     )     allocate ( NoahmpIO%ss_alb_dst3         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst3)    )     allocate ( NoahmpIO%asm_prm_dst3        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst3))     allocate ( NoahmpIO%ext_cff_mss_dst3    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst4)     )     allocate ( NoahmpIO%ss_alb_dst4         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst4)    )     allocate ( NoahmpIO%asm_prm_dst4        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst4))     allocate ( NoahmpIO%ext_cff_mss_dst4    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst5)     )     allocate ( NoahmpIO%ss_alb_dst5         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst5)    )     allocate ( NoahmpIO%asm_prm_dst5        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst5))     allocate ( NoahmpIO%ext_cff_mss_dst5    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%flx_wgt_dir)     )     allocate ( NoahmpIO%flx_wgt_dir         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%flx_wgt_dif)     )     allocate ( NoahmpIO%flx_wgt_dif         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%snowage_tau)     )     allocate ( NoahmpIO%snowage_tau         (NoahmpIO%idx_rhos_max,NoahmpIO%idx_Tgrd_max,NoahmpIO%idx_T_max) )
       if ( .not. allocated (NoahmpIO%snowage_kappa)   )     allocate ( NoahmpIO%snowage_kappa       (NoahmpIO%idx_rhos_max,NoahmpIO%idx_Tgrd_max,NoahmpIO%idx_T_max) )
       if ( .not. allocated (NoahmpIO%snowage_drdt0)   )     allocate ( NoahmpIO%snowage_drdt0       (NoahmpIO%idx_rhos_max,NoahmpIO%idx_Tgrd_max,NoahmpIO%idx_T_max) )
       if ( .not. allocated (NoahmpIO%SNRDSXY) )             allocate ( NoahmpIO%SNRDSXY             (ITS:ITE,-NSNOW+1:0) ) ! snow layer effective grain radius [microns, m-6]
       if ( .not. allocated (NoahmpIO%SNFRXY)  )             allocate ( NoahmpIO%SNFRXY              (ITS:ITE,-NSNOW+1:0) ) ! snow layer rate of snow freezing [mm/s]
       if ( .not. allocated (NoahmpIO%BCPHIXY) )             allocate ( NoahmpIO%BCPHIXY             (ITS:ITE,-NSNOW+1:0) ) 
       if ( .not. allocated (NoahmpIO%BCPHOXY) )             allocate ( NoahmpIO%BCPHOXY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%OCPHIXY) )             allocate ( NoahmpIO%OCPHIXY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%OCPHOXY) )             allocate ( NoahmpIO%OCPHOXY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%DUST1XY) )             allocate ( NoahmpIO%DUST1XY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%DUST2XY) )             allocate ( NoahmpIO%DUST2XY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%DUST3XY) )             allocate ( NoahmpIO%DUST3XY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%DUST4XY) )             allocate ( NoahmpIO%DUST4XY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%DUST5XY) )             allocate ( NoahmpIO%DUST5XY             (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcBCPHIXY) )     allocate ( NoahmpIO%MassConcBCPHIXY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcBCPHOXY) )     allocate ( NoahmpIO%MassConcBCPHOXY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcOCPHIXY) )     allocate ( NoahmpIO%MassConcOCPHIXY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcOCPHOXY) )     allocate ( NoahmpIO%MassConcOCPHOXY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcDUST1XY) )     allocate ( NoahmpIO%MassConcDUST1XY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcDUST2XY) )     allocate ( NoahmpIO%MassConcDUST2XY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcDUST3XY) )     allocate ( NoahmpIO%MassConcDUST3XY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcDUST4XY) )     allocate ( NoahmpIO%MassConcDUST4XY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%MassConcDUST5XY) )     allocate ( NoahmpIO%MassConcDUST5XY     (ITS:ITE,-NSNOW+1:0) )
       if ( .not. allocated (NoahmpIO%DepBChydrophoXY) )     allocate ( NoahmpIO%DepBChydrophoXY     (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepBChydrophiXY) )     allocate ( NoahmpIO%DepBChydrophiXY     (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepOChydrophoXY) )     allocate ( NoahmpIO%DepOChydrophoXY     (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepOChydrophiXY) )     allocate ( NoahmpIO%DepOChydrophiXY     (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepDust1XY)      )     allocate ( NoahmpIO%DepDust1XY          (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepDust2XY)      )     allocate ( NoahmpIO%DepDust2XY          (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepDust3XY)      )     allocate ( NoahmpIO%DepDust3XY          (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepDust4XY)      )     allocate ( NoahmpIO%DepDust4XY          (ITS:ITE           ) )
       if ( .not. allocated (NoahmpIO%DepDust5XY)      )     allocate ( NoahmpIO%DepDust5XY          (ITS:ITE           ) )
    endif

    if ( .not. allocated (NoahmpIO%ALBSNOWDIRXY) ) allocate ( NoahmpIO%ALBSNOWDIRXY (ITS:ITE,1:NUMRAD) ) ! snow albedo (direct)
    if ( .not. allocated (NoahmpIO%ALBSNOWDIFXY) ) allocate ( NoahmpIO%ALBSNOWDIFXY (ITS:ITE,1:NUMRAD) ) ! snow albedo (diffuse)
    if ( .not. allocated (NoahmpIO%ALBSFCDIRXY)  ) allocate ( NoahmpIO%ALBSFCDIRXY  (ITS:ITE,1:NUMRAD) ) ! surface albedo (direct)
    if ( .not. allocated (NoahmpIO%ALBSFCDIFXY)  ) allocate ( NoahmpIO%ALBSFCDIFXY  (ITS:ITE,1:NUMRAD) ) ! surface albedo (diffuse)
    if ( .not. allocated (NoahmpIO%ALBSOILDIRXY) ) allocate ( NoahmpIO%ALBSOILDIRXY (ITS:ITE,1:NUMRAD) ) ! soil albedo (direct)
    if ( .not. allocated (NoahmpIO%ALBSOILDIFXY) ) allocate ( NoahmpIO%ALBSOILDIFXY (ITS:ITE,1:NUMRAD) ) ! soil albedo (diffuse)
    if ( .not. allocated (NoahmpIO%RadSwVisFrac) ) allocate ( NoahmpIO%RadSwVisFrac (ITS:ITE         ) ) ! downward solar radation visible fraction
    if ( .not. allocated (NoahmpIO%RadSwDirFrac) ) allocate ( NoahmpIO%RadSwDirFrac (ITS:ITE         ) ) ! downward solar radation direct fraction

    ! Needed for crop model (OPT_CROP=1)
    if ( .not. allocated (NoahmpIO%PGSXY)     ) allocate ( NoahmpIO%PGSXY      (ITS:ITE  ) )
    if ( .not. allocated (NoahmpIO%CROPCAT)   ) allocate ( NoahmpIO%CROPCAT    (ITS:ITE  ) )
    if ( .not. allocated (NoahmpIO%PLANTING)  ) allocate ( NoahmpIO%PLANTING   (ITS:ITE  ) )
    if ( .not. allocated (NoahmpIO%HARVEST)   ) allocate ( NoahmpIO%HARVEST    (ITS:ITE  ) )
    if ( .not. allocated (NoahmpIO%SEASON_GDD)) allocate ( NoahmpIO%SEASON_GDD (ITS:ITE  ) )
    if ( .not. allocated (NoahmpIO%CROPTYPE)  ) allocate ( NoahmpIO%CROPTYPE   (ITS:ITE,5) )

    ! Needed for Zhang et al. 2022 wetland model (OPT_WETLAND=1 or 2)
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       if ( .not. allocated (NoahmpIO%FSATXY) ) allocate ( NoahmpIO%FSATXY     (ITS:ITE) ) ! saturated fraction of the grid (-)
       if ( .not. allocated (NoahmpIO%WSURFXY)) allocate ( NoahmpIO%WSURFXY    (ITS:ITE) ) ! wetland water storage [mm]
    endif
    if ( NoahmpIO%IOPT_WETLAND == 2 ) then
       if ( .not. allocated (NoahmpIO%FSATMX) ) allocate ( NoahmpIO%FSATMX     (ITS:ITE) ) ! maximum saturated fraction
       if ( .not. allocated (NoahmpIO%WCAP)   ) allocate ( NoahmpIO%WCAP       (ITS:ITE) ) ! maximum wetland capacity [m]
    endif

    !---- For MPAS, urban scheme will be coupled later, so deactivated here
    if (0 == 1) then 
    ! Single- and Multi-layer Urban Models
    if ( NoahmpIO%SF_URBAN_PHYSICS > 0 ) then

       if ( .not. allocated (NoahmpIO%sh_urb2d)   ) allocate ( NoahmpIO%sh_urb2d    (ITS:ITE) ) 
       if ( .not. allocated (NoahmpIO%lh_urb2d)   ) allocate ( NoahmpIO%lh_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%g_urb2d)    ) allocate ( NoahmpIO%g_urb2d     (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%rn_urb2d)   ) allocate ( NoahmpIO%rn_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%ts_urb2d)   ) allocate ( NoahmpIO%ts_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%HRANG)      ) allocate ( NoahmpIO%HRANG       (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%frc_urb2d)  ) allocate ( NoahmpIO%frc_urb2d   (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%utype_urb2d)) allocate ( NoahmpIO%utype_urb2d (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%lp_urb2d)   ) allocate ( NoahmpIO%lp_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%lb_urb2d)   ) allocate ( NoahmpIO%lb_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%hgt_urb2d)  ) allocate ( NoahmpIO%hgt_urb2d   (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%ust)        ) allocate ( NoahmpIO%ust         (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%qc_urb2d)   ) allocate ( NoahmpIO%qc_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%dzr)        ) allocate ( NoahmpIO%dzr         (1:NSOIL) )
       if ( .not. allocated (NoahmpIO%dzb)        ) allocate ( NoahmpIO%dzb         (1:NSOIL) )
       if ( .not. allocated (NoahmpIO%dzg)        ) allocate ( NoahmpIO%dzg         (1:NSOIL) )
       if ( .not. allocated (NoahmpIO%xxxr_urb2d) ) allocate ( NoahmpIO%xxxr_urb2d  (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%xxxb_urb2d) ) allocate ( NoahmpIO%xxxb_urb2d  (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%xxxg_urb2d) ) allocate ( NoahmpIO%xxxg_urb2d  (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%xxxc_urb2d) ) allocate ( NoahmpIO%xxxc_urb2d  (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%tr_urb2d)   ) allocate ( NoahmpIO%tr_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%tb_urb2d)   ) allocate ( NoahmpIO%tb_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%tg_urb2d)   ) allocate ( NoahmpIO%tg_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%tc_urb2d)   ) allocate ( NoahmpIO%tc_urb2d    (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%trl_urb3d)  ) allocate ( NoahmpIO%trl_urb3d   (ITS:ITE,1:NSOIL) )
       if ( .not. allocated (NoahmpIO%tbl_urb3d)  ) allocate ( NoahmpIO%tbl_urb3d   (ITS:ITE,1:NSOIL) )
       if ( .not. allocated (NoahmpIO%tgl_urb3d)  ) allocate ( NoahmpIO%tgl_urb3d   (ITS:ITE,1:NSOIL) )
       if ( .not. allocated (NoahmpIO%sf_ac_urb3d)) allocate ( NoahmpIO%sf_ac_urb3d (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%lf_ac_urb3d)) allocate ( NoahmpIO%lf_ac_urb3d (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%cm_ac_urb3d)) allocate ( NoahmpIO%cm_ac_urb3d (ITS:ITE) )
       if ( .not. allocated (NoahmpIO%sfvent_urb3d))allocate ( NoahmpIO%sfvent_urb3d(ITS:ITE) )
       if ( .not. allocated (NoahmpIO%lfvent_urb3d))allocate ( NoahmpIO%lfvent_urb3d(ITS:ITE) )
         
       if ( NoahmpIO%SF_URBAN_PHYSICS == 1 ) then  ! single layer urban model  
          if ( .not. allocated (NoahmpIO%cmr_sfcdif)   ) allocate ( NoahmpIO%cmr_sfcdif    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%chr_sfcdif)   ) allocate ( NoahmpIO%chr_sfcdif    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%cmc_sfcdif)   ) allocate ( NoahmpIO%cmc_sfcdif    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%chc_sfcdif)   ) allocate ( NoahmpIO%chc_sfcdif    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%cmgr_sfcdif)  ) allocate ( NoahmpIO%cmgr_sfcdif   (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%chgr_sfcdif)  ) allocate ( NoahmpIO%chgr_sfcdif   (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%uc_urb2d)     ) allocate ( NoahmpIO%uc_urb2d      (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%psim_urb2d)   ) allocate ( NoahmpIO%psim_urb2d    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%psih_urb2d)   ) allocate ( NoahmpIO%psih_urb2d    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%u10_urb2d)    ) allocate ( NoahmpIO%u10_urb2d     (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%v10_urb2d)    ) allocate ( NoahmpIO%v10_urb2d     (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%GZ1OZ0_urb2d) ) allocate ( NoahmpIO%GZ1OZ0_urb2d  (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%AKMS_URB2D)   ) allocate ( NoahmpIO%AKMS_URB2D    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%th2_urb2d)    ) allocate ( NoahmpIO%th2_urb2d     (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%q2_urb2d)     ) allocate ( NoahmpIO%q2_urb2d      (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%ust_urb2d)    ) allocate ( NoahmpIO%ust_urb2d     (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%cmcr_urb2d)   ) allocate ( NoahmpIO%cmcr_urb2d    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%tgr_urb2d)    ) allocate ( NoahmpIO%tgr_urb2d     (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%drelr_urb2d)  ) allocate ( NoahmpIO%drelr_urb2d   (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%drelb_urb2d)  ) allocate ( NoahmpIO%drelb_urb2d   (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%drelg_urb2d)  ) allocate ( NoahmpIO%drelg_urb2d   (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%flxhumr_urb2d)) allocate ( NoahmpIO%flxhumr_urb2d (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%flxhumb_urb2d)) allocate ( NoahmpIO%flxhumb_urb2d (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%flxhumg_urb2d)) allocate ( NoahmpIO%flxhumg_urb2d (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%chs)          ) allocate ( NoahmpIO%chs           (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%chs2)         ) allocate ( NoahmpIO%chs2          (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%cqs2)         ) allocate ( NoahmpIO%cqs2          (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%mh_urb2d)     ) allocate ( NoahmpIO%mh_urb2d      (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%stdh_urb2d)   ) allocate ( NoahmpIO%stdh_urb2d    (ITS:ITE) )
          if ( .not. allocated (NoahmpIO%lf_urb2d)     ) allocate ( NoahmpIO%lf_urb2d      (ITS:ITE,4) )
          if ( .not. allocated (NoahmpIO%tgrl_urb3d)   ) allocate ( NoahmpIO%tgrl_urb3d    (ITS:ITE,1:NSOIL) )
          if ( .not. allocated (NoahmpIO%smr_urb3d)    ) allocate ( NoahmpIO%smr_urb3d     (ITS:ITE,1:NSOIL) )
       endif ! SLUCM

       if ( NoahmpIO%SF_URBAN_PHYSICS == 2 .or. NoahmpIO%SF_URBAN_PHYSICS == 3 ) then  ! BEP or BEM urban models
          if ( .not. allocated (NoahmpIO%trb_urb4d)  ) allocate ( NoahmpIO%trb_urb4d   (ITS:ITE,NoahmpIO%urban_map_zrd) )
          if ( .not. allocated (NoahmpIO%tw1_urb4d)  ) allocate ( NoahmpIO%tw1_urb4d   (ITS:ITE,NoahmpIO%urban_map_zwd) )
          if ( .not. allocated (NoahmpIO%tw2_urb4d)  ) allocate ( NoahmpIO%tw2_urb4d   (ITS:ITE,NoahmpIO%urban_map_zwd) )
          if ( .not. allocated (NoahmpIO%tgb_urb4d)  ) allocate ( NoahmpIO%tgb_urb4d   (ITS:ITE,NoahmpIO%urban_map_gd ) )
          if ( .not. allocated (NoahmpIO%sfw1_urb3d) ) allocate ( NoahmpIO%sfw1_urb3d  (ITS:ITE,NoahmpIO%urban_map_zd ) )
          if ( .not. allocated (NoahmpIO%sfw2_urb3d) ) allocate ( NoahmpIO%sfw2_urb3d  (ITS:ITE,NoahmpIO%urban_map_zd ) )
          if ( .not. allocated (NoahmpIO%sfr_urb3d)  ) allocate ( NoahmpIO%sfr_urb3d   (ITS:ITE,NoahmpIO%urban_map_zdf) )
          if ( .not. allocated (NoahmpIO%sfg_urb3d)  ) allocate ( NoahmpIO%sfg_urb3d   (ITS:ITE,NoahmpIO%num_urban_ndm) )
          if ( .not. allocated (NoahmpIO%hi_urb2d)   ) allocate ( NoahmpIO%hi_urb2d    (ITS:ITE,NoahmpIO%num_urban_hi ) )
          if ( .not. allocated (NoahmpIO%theta_urban)) allocate ( NoahmpIO%theta_urban (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%u_urban)    ) allocate ( NoahmpIO%u_urban     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%v_urban)    ) allocate ( NoahmpIO%v_urban     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%dz_urban)   ) allocate ( NoahmpIO%dz_urban    (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%rho_urban)  ) allocate ( NoahmpIO%rho_urban   (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%p_urban)    ) allocate ( NoahmpIO%p_urban     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%a_u_bep)    ) allocate ( NoahmpIO%a_u_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%a_v_bep)    ) allocate ( NoahmpIO%a_v_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%a_t_bep)    ) allocate ( NoahmpIO%a_t_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%a_q_bep)    ) allocate ( NoahmpIO%a_q_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%a_e_bep)    ) allocate ( NoahmpIO%a_e_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%b_u_bep)    ) allocate ( NoahmpIO%b_u_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%b_v_bep)    ) allocate ( NoahmpIO%b_v_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%b_t_bep)    ) allocate ( NoahmpIO%b_t_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%b_q_bep)    ) allocate ( NoahmpIO%b_q_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%b_e_bep)    ) allocate ( NoahmpIO%b_e_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%dlg_bep)    ) allocate ( NoahmpIO%dlg_bep     (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%dl_u_bep)   ) allocate ( NoahmpIO%dl_u_bep    (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%sf_bep)     ) allocate ( NoahmpIO%sf_bep      (ITS:ITE,KTS:KTE) )
          if ( .not. allocated (NoahmpIO%vl_bep)     ) allocate ( NoahmpIO%vl_bep      (ITS:ITE,KTS:KTE) )
       endif ! BEP/BEP-BEM

       if ( NoahmpIO%SF_URBAN_PHYSICS == 3 ) then  ! BEM urban model
          if ( .not. allocated (NoahmpIO%tlev_urb3d)   ) allocate ( NoahmpIO%tlev_urb3d    (ITS:ITE,NoahmpIO%urban_map_bd   ) )
          if ( .not. allocated (NoahmpIO%qlev_urb3d)   ) allocate ( NoahmpIO%qlev_urb3d    (ITS:ITE,NoahmpIO%urban_map_bd   ) )
          if ( .not. allocated (NoahmpIO%tw1lev_urb3d) ) allocate ( NoahmpIO%tw1lev_urb3d  (ITS:ITE,NoahmpIO%urban_map_wd   ) )
          if ( .not. allocated (NoahmpIO%tw2lev_urb3d) ) allocate ( NoahmpIO%tw2lev_urb3d  (ITS:ITE,NoahmpIO%urban_map_wd   ) )
          if ( .not. allocated (NoahmpIO%tglev_urb3d)  ) allocate ( NoahmpIO%tglev_urb3d   (ITS:ITE,NoahmpIO%urban_map_gbd  ) )
          if ( .not. allocated (NoahmpIO%tflev_urb3d)  ) allocate ( NoahmpIO%tflev_urb3d   (ITS:ITE,NoahmpIO%urban_map_fbd  ) )
          if ( .not. allocated (NoahmpIO%sfwin1_urb3d) ) allocate ( NoahmpIO%sfwin1_urb3d  (ITS:ITE,NoahmpIO%urban_map_wd   ) )
          if ( .not. allocated (NoahmpIO%sfwin2_urb3d) ) allocate ( NoahmpIO%sfwin2_urb3d  (ITS:ITE,NoahmpIO%urban_map_wd   ) )
          if ( .not. allocated (NoahmpIO%ep_pv_urb3d)  ) allocate ( NoahmpIO%ep_pv_urb3d   (ITS:ITE                         ) )
          if ( .not. allocated (NoahmpIO%t_pv_urb3d)   ) allocate ( NoahmpIO%t_pv_urb3d    (ITS:ITE,NoahmpIO%urban_map_zdf  ) )
          if ( .not. allocated (NoahmpIO%trv_urb4d)    ) allocate ( NoahmpIO%trv_urb4d     (ITS:ITE,NoahmpIO%urban_map_zgrd ) )
          if ( .not. allocated (NoahmpIO%qr_urb4d)     ) allocate ( NoahmpIO%qr_urb4d      (ITS:ITE,NoahmpIO%urban_map_zgrd ) )
          if ( .not. allocated (NoahmpIO%qgr_urb3d)    ) allocate ( NoahmpIO%qgr_urb3d     (ITS:ITE                         ) )
          if ( .not. allocated (NoahmpIO%tgr_urb3d)    ) allocate ( NoahmpIO%tgr_urb3d     (ITS:ITE                         ) )
          if ( .not. allocated (NoahmpIO%drain_urb4d)  ) allocate ( NoahmpIO%drain_urb4d   (ITS:ITE,NoahmpIO%urban_map_zdf  ) )
          if ( .not. allocated (NoahmpIO%draingr_urb3d)) allocate ( NoahmpIO%draingr_urb3d (ITS:ITE                         ) )
          if ( .not. allocated (NoahmpIO%sfrv_urb3d)   ) allocate ( NoahmpIO%sfrv_urb3d    (ITS:ITE,NoahmpIO%urban_map_zdf  ) )
          if ( .not. allocated (NoahmpIO%lfrv_urb3d)   ) allocate ( NoahmpIO%lfrv_urb3d    (ITS:ITE,NoahmpIO%urban_map_zdf  ) )
          if ( .not. allocated (NoahmpIO%dgr_urb3d)    ) allocate ( NoahmpIO%dgr_urb3d     (ITS:ITE,NoahmpIO%urban_map_zdf  ) )
          if ( .not. allocated (NoahmpIO%dg_urb3d)     ) allocate ( NoahmpIO%dg_urb3d      (ITS:ITE,NoahmpIO%num_urban_ndm  ) )
          if ( .not. allocated (NoahmpIO%lfr_urb3d)    ) allocate ( NoahmpIO%lfr_urb3d     (ITS:ITE,NoahmpIO%urban_map_zdf  ) )
          if ( .not. allocated (NoahmpIO%lfg_urb3d)    ) allocate ( NoahmpIO%lfg_urb3d     (ITS:ITE,NoahmpIO%num_urban_ndm  ) )
       endif ! BEM 

    endif ! urban physics
    endif ! deactivate urban for MPAS

    ! needs to be updated for MPAS-Hydro later
#ifdef WRF_HYDRO
    if ( .not. allocated (NoahmpIO%infxsrt)   ) allocate ( NoahmpIO%infxsrt    (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%sfcheadrt) ) allocate ( NoahmpIO%sfcheadrt  (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%soldrain)  ) allocate ( NoahmpIO%soldrain   (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%qtiledrain)) allocate ( NoahmpIO%qtiledrain (ITS:ITE) )
    if ( .not. allocated (NoahmpIO%ZWATBLE2D) ) allocate ( NoahmpIO%ZWATBLE2D  (ITS:ITE) )
#endif    

    !-------------------------------------------------------------------
    ! Initialize variables with default values 
    !-------------------------------------------------------------------
    
    NoahmpIO%ICE             = undefined_int
    NoahmpIO%IVGTYP          = undefined_int
    NoahmpIO%ISLTYP          = undefined_int
    NoahmpIO%ISNOWXY         = undefined_int
    NoahmpIO%COSZEN          = undefined_real
    NoahmpIO%XLAT            = undefined_real
    NoahmpIO%DX              = undefined_real 
    NoahmpIO%DY              = undefined_real 
    NoahmpIO%DZ8W            = undefined_real
    NoahmpIO%DZS             = undefined_real
    NoahmpIO%ZSOIL           = undefined_real
    NoahmpIO%VEGFRA          = undefined_real
    NoahmpIO%TMN             = undefined_real
    NoahmpIO%XLAND           = undefined_real
    NoahmpIO%XICE            = undefined_real
    NoahmpIO%T_PHY           = undefined_real
    NoahmpIO%QV_CURR         = undefined_real
    NoahmpIO%U_PHY           = undefined_real
    NoahmpIO%V_PHY           = undefined_real
    NoahmpIO%SWDOWN          = undefined_real
    NoahmpIO%SWDDIR          = undefined_real
    NoahmpIO%SWDDIF          = undefined_real
    NoahmpIO%GLW             = undefined_real
    NoahmpIO%P8W             = undefined_real
    NoahmpIO%RAINBL          = undefined_real
    NoahmpIO%SNOWBL          = undefined_real
    NoahmpIO%SR              = undefined_real
    NoahmpIO%RAINCV          = undefined_real
    NoahmpIO%RAINNCV         = undefined_real
    NoahmpIO%RAINSHV         = undefined_real
    NoahmpIO%SNOWNCV         = undefined_real
    NoahmpIO%GRAUPELNCV      = undefined_real
    NoahmpIO%HAILNCV         = undefined_real
    NoahmpIO%QSFC            = undefined_real
    NoahmpIO%TSK             = undefined_real
    NoahmpIO%QFX             = undefined_real
    NoahmpIO%SMSTAV          = undefined_real
    NoahmpIO%SMSTOT          = undefined_real
    NoahmpIO%SMOIS           = undefined_real
    NoahmpIO%SH2O            = undefined_real
    NoahmpIO%TSLB            = undefined_real
    NoahmpIO%SNOW            = undefined_real
    NoahmpIO%SNOWH           = undefined_real
    NoahmpIO%CANWAT          = undefined_real
    NoahmpIO%SMOISEQ         = undefined_real
    NoahmpIO%ALBEDO          = undefined_real
    NoahmpIO%TVXY            = undefined_real
    NoahmpIO%TGXY            = undefined_real
    NoahmpIO%CANICEXY        = undefined_real
    NoahmpIO%CANLIQXY        = undefined_real
    NoahmpIO%EAHXY           = undefined_real
    NoahmpIO%TAHXY           = undefined_real
    NoahmpIO%CMXY            = undefined_real
    NoahmpIO%CHXY            = undefined_real
    NoahmpIO%FWETXY          = undefined_real
    NoahmpIO%SNEQVOXY        = undefined_real
    NoahmpIO%ALBOLDXY        = undefined_real
    NoahmpIO%QSNOWXY         = undefined_real
    NoahmpIO%QRAINXY         = undefined_real
    NoahmpIO%WSLAKEXY        = undefined_real
    NoahmpIO%ZWTXY           = undefined_real
    NoahmpIO%WAXY            = undefined_real
    NoahmpIO%WTXY            = undefined_real
    NoahmpIO%TSNOXY          = undefined_real
    NoahmpIO%SNICEXY         = undefined_real
    NoahmpIO%SNLIQXY         = undefined_real
    NoahmpIO%LFMASSXY        = undefined_real
    NoahmpIO%RTMASSXY        = undefined_real
    NoahmpIO%STMASSXY        = undefined_real
    NoahmpIO%WOODXY          = undefined_real
    NoahmpIO%STBLCPXY        = undefined_real
    NoahmpIO%FASTCPXY        = undefined_real
    NoahmpIO%LAI             = undefined_real
    NoahmpIO%XSAIXY          = undefined_real
    NoahmpIO%XLONG           = undefined_real
    NoahmpIO%SEAICE          = undefined_real
    NoahmpIO%SMCWTDXY        = undefined_real
    NoahmpIO%ZSNSOXY         = undefined_real
    NoahmpIO%GRDFLX          = undefined_real
    NoahmpIO%HFX             = undefined_real
    NoahmpIO%LH              = undefined_real
    NoahmpIO%EMISS           = undefined_real
    NoahmpIO%SNOWC           = undefined_real
    NoahmpIO%T2MVXY          = undefined_real
    NoahmpIO%T2MBXY          = undefined_real
    NoahmpIO%Q2MVXY          = undefined_real
    NoahmpIO%Q2MBXY          = undefined_real
    NoahmpIO%TRADXY          = undefined_real
    NoahmpIO%NEEXY           = undefined_real
    NoahmpIO%GPPXY           = undefined_real
    NoahmpIO%NPPXY           = undefined_real
    NoahmpIO%FVEGXY          = undefined_real
    NoahmpIO%RUNSFXY         = undefined_real
    NoahmpIO%RUNSBXY         = undefined_real
    NoahmpIO%ECANXY          = undefined_real
    NoahmpIO%EDIRXY          = undefined_real
    NoahmpIO%ETRANXY         = undefined_real
    NoahmpIO%FSAXY           = undefined_real
    NoahmpIO%FIRAXY          = undefined_real
    NoahmpIO%APARXY          = undefined_real
    NoahmpIO%PSNXY           = undefined_real
    NoahmpIO%SAVXY           = undefined_real
    NoahmpIO%SAGXY           = undefined_real
    NoahmpIO%RSSUNXY         = undefined_real
    NoahmpIO%RSSHAXY         = undefined_real
    NoahmpIO%BGAPXY          = undefined_real
    NoahmpIO%WGAPXY          = undefined_real
    NoahmpIO%TGVXY           = undefined_real
    NoahmpIO%TGBXY           = undefined_real
    NoahmpIO%CHVXY           = undefined_real
    NoahmpIO%CHBXY           = undefined_real
    NoahmpIO%SHGXY           = undefined_real
    NoahmpIO%SHCXY           = undefined_real
    NoahmpIO%SHBXY           = undefined_real
    NoahmpIO%EVGXY           = undefined_real
    NoahmpIO%EVBXY           = undefined_real
    NoahmpIO%GHVXY           = undefined_real
    NoahmpIO%GHBXY           = undefined_real
    NoahmpIO%IRGXY           = undefined_real
    NoahmpIO%IRCXY           = undefined_real
    NoahmpIO%IRBXY           = undefined_real
    NoahmpIO%TRXY            = undefined_real
    NoahmpIO%EVCXY           = undefined_real
    NoahmpIO%CHLEAFXY        = undefined_real
    NoahmpIO%CHUCXY          = undefined_real
    NoahmpIO%CHV2XY          = undefined_real
    NoahmpIO%CHB2XY          = undefined_real
    NoahmpIO%RS              = undefined_real
    NoahmpIO%CANHSXY         = undefined_real
    NoahmpIO%Z0              = undefined_real
    NoahmpIO%ZNT             = undefined_real
    NoahmpIO%ALBSNOWDIRXY    = undefined_real
    NoahmpIO%ALBSNOWDIFXY    = undefined_real
    NoahmpIO%ALBSFCDIRXY     = undefined_real
    NoahmpIO%ALBSFCDIFXY     = undefined_real
    NoahmpIO%ALBSOILDIRXY    = 0.0
    NoahmpIO%ALBSOILDIFXY    = 0.0
    NoahmpIO%RadSwVisFrac    = 0.5
    NoahmpIO%RadSwDirFrac    = 0.7
    NoahmpIO%TAUSSXY         = 0.0
    NoahmpIO%DEEPRECHXY      = 0.0
    NoahmpIO%RECHXY          = 0.0
    NoahmpIO%ACSNOM          = 0.0
    NoahmpIO%ACSNOW          = 0.0
    NoahmpIO%MP_RAINC        = 0.0
    NoahmpIO%MP_RAINNC       = 0.0
    NoahmpIO%MP_SHCV         = 0.0
    NoahmpIO%MP_SNOW         = 0.0
    NoahmpIO%MP_GRAUP        = 0.0
    NoahmpIO%MP_HAIL         = 0.0
    NoahmpIO%SFCRUNOFF       = 0.0
    NoahmpIO%UDRUNOFF        = 0.0

    ! additional output
    NoahmpIO%PAHXY           = undefined_real
    NoahmpIO%PAHGXY          = undefined_real
    NoahmpIO%PAHBXY          = undefined_real
    NoahmpIO%PAHVXY          = undefined_real
    NoahmpIO%QINTSXY         = undefined_real
    NoahmpIO%QINTRXY         = undefined_real
    NoahmpIO%QDRIPSXY        = undefined_real
    NoahmpIO%QDRIPRXY        = undefined_real
    NoahmpIO%QTHROSXY        = undefined_real
    NoahmpIO%QTHRORXY        = undefined_real
    NoahmpIO%QSNSUBXY        = undefined_real
    NoahmpIO%QSNFROXY        = undefined_real
    NoahmpIO%QSUBCXY         = undefined_real
    NoahmpIO%QFROCXY         = undefined_real
    NoahmpIO%QEVACXY         = undefined_real
    NoahmpIO%QDEWCXY         = undefined_real
    NoahmpIO%QFRZCXY         = undefined_real
    NoahmpIO%QMELTCXY        = undefined_real
    NoahmpIO%QSNBOTXY        = undefined_real
    NoahmpIO%QMELTXY         = undefined_real
    NoahmpIO%FPICEXY         = undefined_real
    NoahmpIO%RAINLSM         = undefined_real
    NoahmpIO%SNOWLSM         = undefined_real
    NoahmpIO%FORCTLSM        = undefined_real
    NoahmpIO%FORCQLSM        = undefined_real
    NoahmpIO%FORCPLSM        = undefined_real
    NoahmpIO%FORCZLSM        = undefined_real
    NoahmpIO%FORCWLSM        = undefined_real
    NoahmpIO%EFLXBXY         = undefined_real
    NoahmpIO%SOILENERGY      = undefined_real
    NoahmpIO%SNOWENERGY      = undefined_real
    NoahmpIO%PONDINGXY       = 0.0
    NoahmpIO%ACC_SSOILXY     = 0.0
    NoahmpIO%ACC_QINSURXY    = 0.0
    NoahmpIO%ACC_QSEVAXY     = 0.0
    NoahmpIO%ACC_ETRANIXY    = 0.0
    NoahmpIO%ACC_DWATERXY    = 0.0
    NoahmpIO%ACC_PRCPXY      = 0.0
    NoahmpIO%ACC_ECANXY      = 0.0
    NoahmpIO%ACC_ETRANXY     = 0.0
    NoahmpIO%ACC_EDIRXY      = 0.0
    NoahmpIO%ACC_GLAFLWXY    = 0.0

    ! MMF Groundwater
    NoahmpIO%TERRAIN         = undefined_real
    NoahmpIO%GVFMIN          = undefined_real
    NoahmpIO%GVFMAX          = undefined_real
    NoahmpIO%MSFTX           = undefined_real
    NoahmpIO%MSFTY           = undefined_real
    NoahmpIO%EQZWT           = undefined_real
    NoahmpIO%RIVERBEDXY      = undefined_real
    NoahmpIO%RIVERCONDXY     = undefined_real
    NoahmpIO%PEXPXY          = undefined_real
    NoahmpIO%FDEPTHXY        = undefined_real
    NoahmpIO%AREAXY          = undefined_real
    NoahmpIO%QRFSXY          = undefined_real
    NoahmpIO%QSPRINGSXY      = undefined_real
    NoahmpIO%QRFXY           = undefined_real
    NoahmpIO%QSPRINGXY       = undefined_real
    NoahmpIO%QSLATXY         = undefined_real
    NoahmpIO%QLATXY          = undefined_real

    ! SNICAR snow albedo
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%ss_alb_snw_drc           = undefined_real
       NoahmpIO%asm_prm_snw_drc          = undefined_real
       NoahmpIO%ext_cff_mss_snw_drc      = undefined_real
       NoahmpIO%ss_alb_snw_dfs           = undefined_real
       NoahmpIO%asm_prm_snw_dfs          = undefined_real
       NoahmpIO%ext_cff_mss_snw_dfs      = undefined_real
       NoahmpIO%ss_alb_bc1               = undefined_real
       NoahmpIO%asm_prm_bc1              = undefined_real
       NoahmpIO%ext_cff_mss_bc1          = undefined_real
       NoahmpIO%ss_alb_bc2               = undefined_real
       NoahmpIO%asm_prm_bc2              = undefined_real
       NoahmpIO%ext_cff_mss_bc2          = undefined_real
       NoahmpIO%ss_alb_oc1               = undefined_real
       NoahmpIO%asm_prm_oc1              = undefined_real
       NoahmpIO%ext_cff_mss_oc1          = undefined_real
       NoahmpIO%ss_alb_oc2               = undefined_real
       NoahmpIO%asm_prm_oc2              = undefined_real
       NoahmpIO%ext_cff_mss_oc2          = undefined_real
       NoahmpIO%ss_alb_dst1              = undefined_real
       NoahmpIO%asm_prm_dst1             = undefined_real
       NoahmpIO%ext_cff_mss_dst1         = undefined_real
       NoahmpIO%ss_alb_dst2              = undefined_real
       NoahmpIO%asm_prm_dst2             = undefined_real
       NoahmpIO%ext_cff_mss_dst2         = undefined_real
       NoahmpIO%ss_alb_dst3              = undefined_real
       NoahmpIO%asm_prm_dst3             = undefined_real
       NoahmpIO%ext_cff_mss_dst3         = undefined_real
       NoahmpIO%ss_alb_dst4              = undefined_real
       NoahmpIO%asm_prm_dst4             = undefined_real
       NoahmpIO%ext_cff_mss_dst4         = undefined_real
       NoahmpIO%ss_alb_dst5              = undefined_real
       NoahmpIO%asm_prm_dst5             = undefined_real
       NoahmpIO%ext_cff_mss_dst5         = undefined_real
       NoahmpIO%flx_wgt_dir              = undefined_real
       NoahmpIO%flx_wgt_dif              = undefined_real
       NoahmpIO%snowage_tau              = undefined_real
       NoahmpIO%snowage_kappa            = undefined_real
       NoahmpIO%snowage_drdt0            = undefined_real
       NoahmpIO%SNRDSXY                  = undefined_real
       NoahmpIO%SNFRXY                   = undefined_real
       NoahmpIO%BCPHOXY                  = undefined_real
       NoahmpIO%BCPHIXY                  = undefined_real
       NoahmpIO%OCPHOXY                  = undefined_real
       NoahmpIO%OCPHIXY                  = undefined_real
       NoahmpIO%DUST1XY                  = undefined_real
       NoahmpIO%DUST2XY                  = undefined_real
       NoahmpIO%DUST3XY                  = undefined_real
       NoahmpIO%DUST4XY                  = undefined_real
       NoahmpIO%DUST5XY                  = undefined_real
       NoahmpIO%MassConcBCPHOXY          = undefined_real
       NoahmpIO%MassConcBCPHIXY          = undefined_real
       NoahmpIO%MassConcOCPHOXY          = undefined_real
       NoahmpIO%MassConcOCPHIXY          = undefined_real
       NoahmpIO%MassConcDUST1XY          = undefined_real
       NoahmpIO%MassConcDUST2XY          = undefined_real
       NoahmpIO%MassConcDUST3XY          = undefined_real
       NoahmpIO%MassConcDUST4XY          = undefined_real
       NoahmpIO%MassConcDUST5XY          = undefined_real
       NoahmpIO%DepBChydrophoXY          = undefined_real
       NoahmpIO%DepBChydrophiXY          = undefined_real
       NoahmpIO%DepOChydrophoXY          = undefined_real
       NoahmpIO%DepOChydrophiXY          = undefined_real
       NoahmpIO%DepDust1XY               = undefined_real
       NoahmpIO%DepDust2XY               = undefined_real
       NoahmpIO%DepDust3XY               = undefined_real
       NoahmpIO%DepDust4XY               = undefined_real
       NoahmpIO%DepDust5XY               = undefined_real
    endif

    ! crop model
    NoahmpIO%PGSXY           = undefined_int
    NoahmpIO%CROPCAT         = undefined_int
    NoahmpIO%PLANTING        = undefined_real
    NoahmpIO%HARVEST         = undefined_real
    NoahmpIO%SEASON_GDD      = undefined_real
    NoahmpIO%CROPTYPE        = undefined_real

    ! tile drainage
    NoahmpIO%QTDRAIN         = 0.0
    NoahmpIO%TD_FRACTION     = undefined_real

    ! irrigation
    NoahmpIO%IRFRACT         = 0.0
    NoahmpIO%SIFRACT         = 0.0
    NoahmpIO%MIFRACT         = 0.0
    NoahmpIO%FIFRACT         = 0.0
    NoahmpIO%IRNUMSI         = 0
    NoahmpIO%IRNUMMI         = 0
    NoahmpIO%IRNUMFI         = 0
    NoahmpIO%IRWATSI         = 0.0
    NoahmpIO%IRWATMI         = 0.0
    NoahmpIO%IRWATFI         = 0.0
    NoahmpIO%IRELOSS         = 0.0
    NoahmpIO%IRSIVOL         = 0.0
    NoahmpIO%IRMIVOL         = 0.0
    NoahmpIO%IRFIVOL         = 0.0
    NoahmpIO%IRRSPLH         = 0.0
    NoahmpIO%LOCTIM          = undefined_real

    ! wetland model (Zhang et al. 2022)
    if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       NoahmpIO%FSATXY       = undefined_real
       NoahmpIO%WSURFXY      = undefined_real
    endif
    if ( NoahmpIO%IOPT_WETLAND == 2 ) then
       NoahmpIO%FSATMX       = undefined_real
       NoahmpIO%WCAP         = undefined_real
    endif

    ! spatial varying soil texture
    if ( NoahmpIO%IOPT_SOIL > 1 ) then
       NoahmpIO%SOILCL1      = undefined_real
       NoahmpIO%SOILCL2      = undefined_real
       NoahmpIO%SOILCL3      = undefined_real
       NoahmpIO%SOILCL4      = undefined_real
       NoahmpIO%SOILCOMP     = undefined_real
    endif

    !--- For WRF, this is not needed, so deactivated here
    if (0 == 1) then
    ! urban model 
    if ( NoahmpIO%SF_URBAN_PHYSICS > 0 ) then
       NoahmpIO%JULDAY        = undefined_int_neg
       NoahmpIO%IRI_URBAN     = undefined_int_neg
       NoahmpIO%utype_urb2d   = undefined_int_neg
       NoahmpIO%HRANG         = undefined_real_neg
       NoahmpIO%DECLIN        = undefined_real_neg
       NoahmpIO%sh_urb2d      = undefined_real_neg
       NoahmpIO%lh_urb2d      = undefined_real_neg
       NoahmpIO%g_urb2d       = undefined_real_neg
       NoahmpIO%rn_urb2d      = undefined_real_neg
       NoahmpIO%ts_urb2d      = undefined_real_neg
       NoahmpIO%GMT           = undefined_real_neg
       NoahmpIO%frc_urb2d     = undefined_real_neg
       NoahmpIO%lp_urb2d      = undefined_real_neg
       NoahmpIO%lb_urb2d      = undefined_real_neg
       NoahmpIO%hgt_urb2d     = undefined_real_neg
       NoahmpIO%ust           = undefined_real_neg
       NoahmpIO%qc_urb2d      = undefined_real_neg
       NoahmpIO%dzr           = undefined_real_neg
       NoahmpIO%dzb           = undefined_real_neg
       NoahmpIO%dzg           = undefined_real_neg
       NoahmpIO%xxxr_urb2d    = undefined_real_neg
       NoahmpIO%xxxb_urb2d    = undefined_real_neg
       NoahmpIO%xxxg_urb2d    = undefined_real_neg
       NoahmpIO%xxxc_urb2d    = undefined_real_neg
       NoahmpIO%tr_urb2d      = undefined_real_neg
       NoahmpIO%tb_urb2d      = undefined_real_neg
       NoahmpIO%tg_urb2d      = undefined_real_neg
       NoahmpIO%tc_urb2d      = undefined_real_neg
       NoahmpIO%trl_urb3d     = undefined_real_neg
       NoahmpIO%tbl_urb3d     = undefined_real_neg
       NoahmpIO%tgl_urb3d     = undefined_real_neg
       NoahmpIO%sf_ac_urb3d   = undefined_real_neg
       NoahmpIO%lf_ac_urb3d   = undefined_real_neg
       NoahmpIO%cm_ac_urb3d   = undefined_real_neg
       NoahmpIO%sfvent_urb3d  = undefined_real_neg
       NoahmpIO%lfvent_urb3d  = undefined_real_neg

       if ( NoahmpIO%SF_URBAN_PHYSICS == 1 ) then  ! single layer urban model
          NoahmpIO%cmr_sfcdif    = 1.0e-4
          NoahmpIO%chr_sfcdif    = 1.0e-4
          NoahmpIO%cmc_sfcdif    = 1.0e-4
          NoahmpIO%chc_sfcdif    = 1.0e-4
          NoahmpIO%cmgr_sfcdif   = 1.0e-4
          NoahmpIO%chgr_sfcdif   = 1.0e-4
          NoahmpIO%uc_urb2d      = undefined_real_neg
          NoahmpIO%psim_urb2d    = undefined_real_neg
          NoahmpIO%psih_urb2d    = undefined_real_neg
          NoahmpIO%u10_urb2d     = undefined_real_neg
          NoahmpIO%v10_urb2d     = undefined_real_neg
          NoahmpIO%GZ1OZ0_urb2d  = undefined_real_neg
          NoahmpIO%AKMS_URB2D    = undefined_real_neg
          NoahmpIO%th2_urb2d     = undefined_real_neg
          NoahmpIO%q2_urb2d      = undefined_real_neg
          NoahmpIO%ust_urb2d     = undefined_real_neg
          NoahmpIO%cmcr_urb2d    = undefined_real_neg
          NoahmpIO%tgr_urb2d     = undefined_real_neg
          NoahmpIO%tgrl_urb3d    = undefined_real_neg
          NoahmpIO%smr_urb3d     = undefined_real_neg
          NoahmpIO%drelr_urb2d   = undefined_real_neg
          NoahmpIO%drelb_urb2d   = undefined_real_neg
          NoahmpIO%drelg_urb2d   = undefined_real_neg
          NoahmpIO%flxhumr_urb2d = undefined_real_neg
          NoahmpIO%flxhumb_urb2d = undefined_real_neg
          NoahmpIO%flxhumg_urb2d = undefined_real_neg
          NoahmpIO%chs           = 1.0e-4
          NoahmpIO%chs2          = 1.0e-4
          NoahmpIO%cqs2          = 1.0e-4
          NoahmpIO%mh_urb2d      = undefined_real_neg
          NoahmpIO%stdh_urb2d    = undefined_real_neg
          NoahmpIO%lf_urb2d      = undefined_real_neg
       endif ! SLUCM
       if ( NoahmpIO%SF_URBAN_PHYSICS == 2 .or. NoahmpIO%SF_URBAN_PHYSICS == 3 ) then  ! BEP or BEM urban models
          NoahmpIO%trb_urb4d     = undefined_real_neg
          NoahmpIO%tw1_urb4d     = undefined_real_neg
          NoahmpIO%tw2_urb4d     = undefined_real_neg
          NoahmpIO%tgb_urb4d     = undefined_real_neg
          NoahmpIO%sfw1_urb3d    = undefined_real_neg
          NoahmpIO%sfw2_urb3d    = undefined_real_neg
          NoahmpIO%sfr_urb3d     = undefined_real_neg
          NoahmpIO%sfg_urb3d     = undefined_real_neg
          NoahmpIO%hi_urb2d      = undefined_real_neg
          NoahmpIO%theta_urban   = undefined_real_neg
          NoahmpIO%u_urban       = undefined_real_neg
          NoahmpIO%v_urban       = undefined_real_neg
          NoahmpIO%dz_urban      = undefined_real_neg
          NoahmpIO%rho_urban     = undefined_real_neg
          NoahmpIO%p_urban       = undefined_real_neg
          NoahmpIO%a_u_bep       = undefined_real_neg
          NoahmpIO%a_v_bep       = undefined_real_neg
          NoahmpIO%a_t_bep       = undefined_real_neg
          NoahmpIO%a_q_bep       = undefined_real_neg
          NoahmpIO%a_e_bep       = undefined_real_neg
          NoahmpIO%b_u_bep       = undefined_real_neg
          NoahmpIO%b_v_bep       = undefined_real_neg
          NoahmpIO%b_t_bep       = undefined_real_neg
          NoahmpIO%b_q_bep       = undefined_real_neg
          NoahmpIO%b_e_bep       = undefined_real_neg
          NoahmpIO%dlg_bep       = undefined_real_neg
          NoahmpIO%dl_u_bep      = undefined_real_neg
          NoahmpIO%sf_bep        = undefined_real_neg
          NoahmpIO%vl_bep        = undefined_real_neg
       endif ! BEP/BEP-BEM
       if ( NoahmpIO%SF_URBAN_PHYSICS == 3 ) then  ! BEM urban model
          NoahmpIO%tlev_urb3d    = undefined_real_neg
          NoahmpIO%qlev_urb3d    = undefined_real_neg
          NoahmpIO%tw1lev_urb3d  = undefined_real_neg
          NoahmpIO%tw2lev_urb3d  = undefined_real_neg
          NoahmpIO%tglev_urb3d   = undefined_real_neg
          NoahmpIO%tflev_urb3d   = undefined_real_neg
          NoahmpIO%sfwin1_urb3d  = undefined_real_neg
          NoahmpIO%sfwin2_urb3d  = undefined_real_neg
          NoahmpIO%ep_pv_urb3d   = undefined_real_neg
          NoahmpIO%t_pv_urb3d    = undefined_real_neg
          NoahmpIO%trv_urb4d     = undefined_real_neg
          NoahmpIO%qr_urb4d      = undefined_real_neg
          NoahmpIO%qgr_urb3d     = undefined_real_neg
          NoahmpIO%tgr_urb3d     = undefined_real_neg
          NoahmpIO%drain_urb4d   = undefined_real_neg
          NoahmpIO%draingr_urb3d = undefined_real_neg
          NoahmpIO%sfrv_urb3d    = undefined_real_neg
          NoahmpIO%lfrv_urb3d    = undefined_real_neg
          NoahmpIO%dgr_urb3d     = undefined_real_neg
          NoahmpIO%dg_urb3d      = undefined_real_neg
          NoahmpIO%lfr_urb3d     = undefined_real_neg
          NoahmpIO%lfg_urb3d     = undefined_real_neg
       endif ! BEM 
    endif ! urban physics
    endif ! deactivate urban for WRF

    NoahmpIO%XLAND             = 1.0      ! water = 2.0, land = 1.0
    NoahmpIO%XICE              = 0.0      ! fraction of grid that is seaice
    NoahmpIO%XICE_THRESHOLD    = 0.5      ! fraction of grid determining seaice (from WRF)
    NoahmpIO%SLOPETYP          = 1        ! soil parameter slope type
    NoahmpIO%soil_update_steps = 1        ! number of model time step to update soil proces
    NoahmpIO%calculate_soil    = .false.  ! index for if do soil process
    NoahmpIO%ITIMESTEP         = 0        ! model time step count

#ifdef WRF_HYDRO
    NoahmpIO%infxsrt         = 0.0
    NoahmpIO%sfcheadrt       = 0.0 
    NoahmpIO%soldrain        = 0.0
    NoahmpIO%qtiledrain      = 0.0
    NoahmpIO%ZWATBLE2D       = 0.0
#endif 
   
    end associate
 
  end subroutine NoahmpIOVarInitDefault

end module NoahmpIOVarInitMod
