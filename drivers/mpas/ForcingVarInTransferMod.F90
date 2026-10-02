module ForcingVarInTransferMod

!!! Transfer input 2-D NoahmpIO Forcing variables to 1-D column variable
!!! 1-D variables should be first defined in /src/ForcingVarType.F90
!!! 2-D variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Sep 13, 2026: NoahmpIO%xx change to 1-D vector for MPAS, Cenlin He (NCAR)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpVarType

  implicit none

contains

!=== initialize with input data or table values

  subroutine ForcingVarInTransfer(noahmp, NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    type(noahmp_type),   intent(inout) :: noahmp
    
    ! local variables
    real(kind=kind_noahmp)              :: PrecipOtherRefHeight  ! other precipitation, e.g. fog [mm/s] at reference height
    real(kind=kind_noahmp)              :: PrecipTotalRefHeight  ! total precipitation [mm/s] at reference height

! ---------------------------------------------------------------
    associate(                                           &
              I      => noahmp%config%domain%GridIndexI  &
             )
! ---------------------------------------------------------------

    noahmp%forcing%TemperatureAirRefHeight = NoahmpIO%T_PHY(I,1)
    noahmp%forcing%WindEastwardRefHeight   = NoahmpIO%U_PHY(I,1)
    noahmp%forcing%WindNorthwardRefHeight  = NoahmpIO%V_PHY(I,1)
    noahmp%forcing%SpecHumidityRefHeight   = NoahmpIO%QV_CURR(I,1)/(1.0+NoahmpIO%QV_CURR(I,1))  ! convert from mixing ratio to specific humidity
    noahmp%forcing%PressureAirRefHeight    = (NoahmpIO%P8W(I,1) + NoahmpIO%P8W(I,2)) * 0.5 ! air pressure at middle point of lowest atmos model layer
    noahmp%forcing%PressureAirSurface      = NoahmpIO%P8W      (I,1)
    noahmp%forcing%RadLwDownRefHeight      = NoahmpIO%GLW      (I)
    noahmp%forcing%RadSwDownRefHeight      = NoahmpIO%SWDOWN   (I)
    noahmp%forcing%TemperatureSoilBottom   = NoahmpIO%TMN      (I)

    ! treat different precipitation types
    PrecipTotalRefHeight                   = NoahmpIO%RAINBL   (I) / NoahmpIO%DTBL                ! convert precip unit from mm/timestep to mm/s
    noahmp%forcing%PrecipConvRefHeight     = NoahmpIO%MP_RAINC (I) / NoahmpIO%DTBL
    noahmp%forcing%PrecipNonConvRefHeight  = NoahmpIO%MP_RAINNC(I) / NoahmpIO%DTBL
    noahmp%forcing%PrecipShConvRefHeight   = NoahmpIO%MP_SHCV  (I) / NoahmpIO%DTBL
    noahmp%forcing%PrecipSnowRefHeight     = NoahmpIO%MP_SNOW  (I) / NoahmpIO%DTBL
    noahmp%forcing%PrecipGraupelRefHeight  = NoahmpIO%MP_GRAUP (I) / NoahmpIO%DTBL
    noahmp%forcing%PrecipHailRefHeight     = NoahmpIO%MP_HAIL  (I) / NoahmpIO%DTBL
    ! treat other precipitation (e.g. fog) contained in total precipitation
    PrecipOtherRefHeight                   = PrecipTotalRefHeight - noahmp%forcing%PrecipConvRefHeight - &
                                             noahmp%forcing%PrecipNonConvRefHeight - noahmp%forcing%PrecipShConvRefHeight
    PrecipOtherRefHeight                   = max(0.0, PrecipOtherRefHeight)
    noahmp%forcing%PrecipNonConvRefHeight  = noahmp%forcing%PrecipNonConvRefHeight + PrecipOtherRefHeight
    noahmp%forcing%PrecipSnowRefHeight     = noahmp%forcing%PrecipSnowRefHeight + PrecipOtherRefHeight * NoahmpIO%SR(I)

    ! downward solar radiation direct/diffuse and visible/NIR partition
    noahmp%forcing%RadSwDirFrac            = NoahmpIO%RadSwDirFrac(I)
    noahmp%forcing%RadSwVisFrac            = NoahmpIO%RadSwVisFrac(I)

    ! SNICAR aerosol deposition flux forcing
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( noahmp%config%nmlist%FlagSnicarAerosolReadTable .eqv. .true. ) then 
          noahmp%forcing%DepBChydropho     = NoahmpIO%DepBChydropho_TABLE
          noahmp%forcing%DepBChydrophi     = NoahmpIO%DepBChydrophi_TABLE
          noahmp%forcing%DepOChydropho     = NoahmpIO%DepOChydropho_TABLE
          noahmp%forcing%DepOChydrophi     = NoahmpIO%DepOChydrophi_TABLE
          noahmp%forcing%DepDust1          = NoahmpIO%DepDust1_TABLE
          noahmp%forcing%DepDust2          = NoahmpIO%DepDust2_TABLE
          noahmp%forcing%DepDust3          = NoahmpIO%DepDust3_TABLE
          noahmp%forcing%DepDust4          = NoahmpIO%DepDust4_TABLE
          noahmp%forcing%DepDust5          = NoahmpIO%DepDust5_TABLE
       else
          noahmp%forcing%DepBChydropho     = NoahmpIO%DepBChydrophoXY(I)
          noahmp%forcing%DepBChydrophi     = NoahmpIO%DepBChydrophiXY(I)
          noahmp%forcing%DepOChydropho     = NoahmpIO%DepOChydrophoXY(I)
          noahmp%forcing%DepOChydrophi     = NoahmpIO%DepOChydrophiXY(I)
          noahmp%forcing%DepDust1          = NoahmpIO%DepDust1XY(I)
          noahmp%forcing%DepDust2          = NoahmpIO%DepDust2XY(I)
          noahmp%forcing%DepDust3          = NoahmpIO%DepDust3XY(I)
          noahmp%forcing%DepDust4          = NoahmpIO%DepDust4XY(I)
          noahmp%forcing%DepDust5          = NoahmpIO%DepDust5XY(I)
       endif
    endif

    end associate
 
  end subroutine ForcingVarInTransfer

end module ForcingVarInTransferMod
