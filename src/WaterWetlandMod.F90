module WaterWetlandMod

!!! Calculate surface runoff based on first-layer soil moisture with wetland scheme (Zhang et al., 2020)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine WaterWetland(noahmp,TimeStep)

! ------------------------ Code history --------------------------------------------------
! Originally embeded in SOILWATER subroutine instead of as a separate subroutine
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout)      :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep                      ! timestep (may not be the same as model timestep)

! local variables
    real(kind=kind_noahmp)           :: WCAP    ! surface wetland storage capacity
    real(kind=kind_noahmp)           :: LATHEA  ! latent heat vap./sublimation
    real(kind=kind_noahmp)           :: SC      ! slope of saturation vapor pressure curve
    real(kind=kind_noahmp)           :: QE      ! qevaporation from P-T method, energy flux [W/m2]
    real(kind=kind_noahmp)           :: GAM     ! gamma psychrometic constant
    real(kind=kind_noahmp)           :: FWAT    ! evaporation heat from surface [W/m2]
    real(kind=kind_noahmp)           :: EWAT    ! evaporation water from surface [mm/s]
! --------------------------------------------------------------------
    associate(                                                           &
              TemperatureSfc    => noahmp%energy%state%TemperatureSfc      ,& ! in, surface temperature [K]
              RadSwAbsSfc       => noahmp%energy%flux%RadSwAbsSfc          ,& ! in,  total absorbed solar radiation [W/m2]
              RadSwReflSfc      => noahmp%energy%flux%RadSwReflSfc         ,& ! in,  total reflected solar radiation [W/m2]
              RadLwNetSfc       => noahmp%energy%flux%RadLwNetSfc          ,& ! in,  total net longwave rad [W/m2] (+ to atm)
              HeatGroundTot     => noahmp%energy%flux%HeatGroundTot        ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
              SoilSaturateFrac  => noahmp%water%state%SoilSaturateFrac      ,& ! in,  fractional saturated area for soil moisture
              HeatSensibleSfc   => noahmp%energy%flux%HeatSensibleSfc       ,& ! in,  total sensible heat [W/m2] (+ to atm)
              WetlandCapMax     => noahmp%water%param%WetlandCapMax         ,& ! in,  maximum wetland capacity [mm]
              RunoffSurface     => noahmp%water%flux%RunoffSurface          ,& ! inout, surface runoff [m/s]
              WaterStorageWetland => noahmp%water%state%WaterStorageWetland ,& ! inout,  
              EvapGroundNet       => noahmp%water%flux%EvapGroundNet        ,& ! inout, accumulated net ground evaporation per soil timestep [mm]
              HeatLatentGrd       => noahmp%energy%flux%HeatLatentGrd       ,& ! inout, ground evaporation heat flux [W/m2] (+ to atm)
              HeatLatentCanEvap   => noahmp%energy%flux%HeatLatentCanEvap   ,& ! inout, canopy evaporation heat flux [W/m2] (+ to atm)
              HeatLatentCanTransp => noahmp%energy%flux%HeatLatentCanTransp  & ! inout, canopy transpiration heat flux [W/m2] (+ to atm)
             )
! ----------------------------------------------------------------------

! set initial value
 EWAT = 0.0
 FWAT = 0.0
 WCAP = WetlandCapMax * 1000.0

! set psychrometric constant and SC
 IF (TemperatureSfc .GT. ConstFreezePoint) THEN 
    LATHEA = ConstLatHeatEvap
 ELSE
    LATHEA = ConstLatHeatSublim
 ENDIF

! determine gamma
 IF (TemperatureSfc .LT. 273.15+26.85) THEN
    GAM = 0.00040
 ELSE 
    GAM = 0.00041
 ENDIF 

! calculate slope SC based on surface temperature
 if (TemperatureSfc .LT. 273.15+6.85) then
    SC = 0.00022
 elseif (TemperatureSfc .LT. 273.15+16.85) then
    SC = 0.00042
 elseif (TemperatureSfc .LT. 273.15+26.85) then
    SC = 0.00078 
 else 
    SC = 0.00132
 endif 
 

 QE = 1.26*SC*(RadSwAbsSfc-RadLwNetSfc-HeatGroundTot)/(SC+GAM)           ! QE Potential latent heat W/m2 P-T method
 WaterStorageWetland = WaterStorageWetland + RunoffSurface     
 IF (WaterStorageWetland .GT. QE*TimeStep/LATHEA*SoilSaturateFrac) THEN  ! if current wetland storage is larger than PET rate 
    EWAT = max(QE/LATHEA*SoilSaturateFrac,0.0)                           ! EWAT mm/s
 ELSE                                                                    ! if current wetland storage is less than PET rate
    EWAT = max(WaterStorageWetland/TimeStep, 0.0)                        ! use all of it
 ENDIF 
 WaterStorageWetland = max(WaterStorageWetland-EWAT*TimeStep,0.0)        ! adjust surface wetland storage

 ! adjuest energy and water balance                            
 FWAT = EWAT * LATHEA                                                    ! convert evaporation to latent heat flux 
 HeatSensibleSfc  = HeatSensibleSfc  - FWAT                              ! reduce sensible heat flux
 HeatLatentGrd    = HeatLatentGrd + FWAT                                 ! increase direct evaporation
 RunoffSurface = max(WaterStorageWetland-WetlandCapMax*1000.0, 0.0)      ! excessive storage becomes runoff
 EvapGroundNet = EvapGroundNet + EWAT                                    ! increase direct evaporation
 WaterStorageWetland = min(WetlandCapMax*1000.0,WaterStorageWetland)

    end associate

  end subroutine WaterWetland

end module WaterWetlandMod
