module DynaRootGrowDieMod

!!! compute soil water transpiration factor that will be used for
!!! stomata resistance and evapotranspiration calculations

  use Machine
  use NoahmpVarType

  implicit none

contains

  subroutine DynaRootGrowDie(noahmp)

! ------------------------ Code history -----------------------------------
! Original code: Carolina Bieri
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type), intent(inout)  :: noahmp

! local variables
    integer                            :: IndSoil         ! loop index
    integer                            :: WaterTableLayer ! layer that contains water table
    integer                            :: AboveWaterTableLayer ! layer above water table
    integer                            :: DynaRootLayer   ! bottom layer of root zone (root depth)
    integer                            :: SoilIceFactor   ! multiplication factor corresponding to IndicatorIceSfc
    real(kind=kind_noahmp)             :: LayerMidpoint   ! depth to midpoint of soil layer [m]
    real(kind=kind_noahmp)             :: MaximumInactiveDays ! max number of days without active roots until layer is inactive [s]
    real(kind=kind_noahmp)             :: MaximumEaseFunction ! maximum column-wise ease function value
    real(kind=kind_noahmp)             :: SumEaseFunction     ! column-wise sum of ease function
    real(kind=kind_noahmp)             :: EffectiveCanopyHeight ! canopy height used in ease function calculation [m]
    real(kind=kind_noahmp), parameter  :: LeafWiltMatPotential = -150.0 ! leaf wilting point matric potential [m]

    real(kind=kind_noahmp), allocatable, dimension(:) :: DynaRootMask    ! indicator for presence of active roots (1=roots active; 0=roots not active)
    real(kind=kind_noahmp), allocatable, dimension(:) :: ThicknessSoilLayerTmp   ! ThicknessSoilLayer accounting for location of water table, only for root scheme [mm]    
    real(kind=kind_noahmp), allocatable, dimension(:) :: InactiveDays ! number of days without active roots [s]
! --------------------------------------------------------------------
    associate(                                                                             &
              ThicknessSoilLayer        => noahmp%config%domain%ThicknessSoilLayer        ,& ! in,  soil layer thickness [m]
              DepthSoilLayer            => noahmp%config%domain%DepthSoilLayer            ,& ! in,  layer-bottom depth from soil surface [m]
              NumSoilLayer              => noahmp%config%domain%NumSoilLayer              ,& ! in,  number of soil layers
              MainTimeStep              => noahmp%config%domain%MainTimeStep              ,& ! in,  main Noah-MP model time step [s]
              IndicatorIceSfc           => noahmp%config%domain%IndicatorIceSfc           ,& ! in,  indicator for sea ice and land ice (1=sea ice; -1=land ice; 0=non-ice)
              WaterTableDepth           => noahmp%water%state%WaterTableDepth             ,& ! in,  depth to water table [m]
              SoilMatPotential          => noahmp%water%state%SoilMatPotential            ,& ! in,  soil matrix potential [m]
              HeightCanopyTop           => noahmp%energy%param%HeightCanopyTop            ,& ! in,  height of canopy top [m]
              InactiveTimeSteps         => noahmp%water%state%InactiveTimeSteps           ,& ! inout, number of time steps without active roots               
              NumSoilLayerRoot          => noahmp%water%param%NumSoilLayerRoot            ,& ! inout,  number of soil layers with root present
              EaseFunction              => noahmp%water%state%EaseFunction                ,& ! out, ease function
              RootActivity              => noahmp%water%state%RootActivity                 & ! out, root activity function
              )
! ----------------------------------------------------------------------

    if (.not. allocated(DynaRootMask)  ) allocate(DynaRootMask  (1:NumSoilLayer))
    if (.not. allocated(ThicknessSoilLayerTmp)  ) allocate(ThicknessSoilLayerTmp  (1:NumSoilLayer))
    if (.not. allocated(InactiveDays)  ) allocate(InactiveDays  (1:NumSoilLayer))

    ! initialize
    DynaRootMask          = 0.0
    ThicknessSoilLayerTmp = undefined_real
    MaximumInactiveDays   = 31536000.0 ! seconds in a year 
    EaseFunction          = 0.0
    RootActivity          = 0.0
    InactiveDays          = 0.0

    ThicknessSoilLayerTmp = ThicknessSoilLayer(1:NumSoilLayer)   ! soil thickness variable for use in root scheme

    ! determine WaterTableLayer
    IF(WaterTableDepth + 1.0E-6 .lt. DepthSoilLayer(NumSoilLayer)) THEN  ! set WaterTableLayer to NumSoilLayer if water table below resolved layers
        WaterTableLayer = NumSoilLayer
    ELSE
        DO IndSoil=NumSoilLayer,1,-1                                 ! find WaterTableLayer if within resolved layers
                IF(WaterTableDepth + 1.0E-6 < DepthSoilLayer(IndSoil)) EXIT
        ENDDO
        AboveWaterTableLayer = IndSoil                               ! layer above water table
        WaterTableLayer = MIN(AboveWaterTableLayer+1,NumSoilLayer)  

        ! Adjust layer thickness accordingly
        IF (AboveWaterTableLayer .ge. 1) ThicknessSoilLayerTmp(WaterTableLayer) = -1*(WaterTableDepth-DepthSoilLayer(AboveWaterTableLayer))
    END IF

    ! determine DynaRootLayer - bottom layer of root zone
    InactiveDays(1:NumSoilLayer) = InactiveTimeSteps(1:NumSoilLayer)*MainTimeStep ! Convert InactiveTimeSteps to InactiveDays (in seconds)
    DO IndSoil = NumSoilLayer,1,-1
        IF(InactiveDays(IndSoil) .le. MaximumInactiveDays)EXIT
    END DO
    DynaRootLayer = MIN(MAX(IndSoil,1),NumSoilLayer)
    print *, "DynaRootLayer:", DynaRootLayer

    EffectiveCanopyHeight = HeightCanopyTop*(2./3.)
    ! calculation of ease function and root activity
    DO IndSoil = 1, MIN(MIN(DynaRootLayer+1,NumSoilLayer), WaterTableLayer)
        IF (IndSoil==1) THEN            ! get layer midpoints
                LayerMidpoint = -0.05
        ELSE
                LayerMidpoint =  0.5 * (DepthSoilLayer(IndSoil-1) + DepthSoilLayer(IndSoil))
        ENDIF

        IF(InactiveDays(IndSoil) .le. MaximumInactiveDays) DynaRootMask(IndSoil) = 1.    ! set root mask accordingly

        ! account for frozen soil
        IF(IndicatorIceSfc.eq.0) THEN
                SoilIceFactor = 1.
        ELSE
                SoilIceFactor = 0.
        ENDIF

        ! calculate ease function
        ! EaseFunction(IndSoil) = MAX(-( LeafWiltMatPotential - SoilMatPotential(IndSoil) )*SoilIceFactor / ( EffectiveCanopyHeight-LayerMidpoint ), 0.)
        EaseFunction(IndSoil) = -( LeafWiltMatPotential - SoilMatPotential(IndSoil) )*SoilIceFactor / ( EffectiveCanopyHeight-LayerMidpoint )
        print *, "Ease function:", EaseFunction(IndSoil)
        print *, "Matric potential:", SoilMatPotential(IndSoil)
    
    END DO

    ! to grow roots anew, the layer has to be easiest to get water from than the
    ! current active layers with roots
    MaximumEaseFunction = MAXVAL(EaseFunction, DynaRootMask == 1.)

    ! eliminate small root activity
    WHERE(EaseFunction .lt. 0.001*MaximumEaseFunction) EaseFunction = 0.

    ! allow root water uptake to extend deeper only if needed
    DO IndSoil = 1, MIN(DynaRootLayer+1,NumSoilLayer)
        IF(InactiveDays(IndSoil) .gt. MaximumInactiveDays .and. EaseFunction(IndSoil) .lt. MaximumEaseFunction) EaseFunction(IndSoil) = 0.
    ENDDO
    ! update root zone index
    IF(DynaRootLayer .lt. NumSoilLayer)THEN
        IF(EaseFunction(DynaRootLayer+1) .gt. 0.) DynaRootLayer = DynaRootLayer+1
    ENDIF

    ! calculate root activity
    SumEaseFunction = SUM(EaseFunction*ThicknessSoilLayerTmp(1:NumSoilLayer))
    IF(SumEaseFunction .eq. 0.)THEN
        RootActivity = 0.
    ELSE
        RootActivity = MIN(MAX((EaseFunction*ThicknessSoilLayerTmp(1:NumSoilLayer)) / SumEaseFunction , 0. ), 1. )
    ENDIF

    ! Increment or restart inactive timestep count based on root activity
    DO IndSoil = 1, NumSoilLayer
        IF(EaseFunction(IndSoil) .eq. 0.)THEN
                InactiveTimeSteps(IndSoil) = InactiveTimeSteps(IndSoil) + 1
        ELSE
                InactiveTimeSteps(IndSoil) = 0
        ENDIF
    ENDDO

         !INACTIVEDAYS = MIN(INACTIVEDAYS,MaximumInactiveDays+1)

    NumSoilLayerRoot = DynaRootLayer

   end associate
  
  end subroutine DynaRootGrowDie 

end module DynaRootGrowDieMod
