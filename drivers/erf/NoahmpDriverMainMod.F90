module NoahmpDriverMainMod

  use Machine
  use NoahmpVarType
  use NoahmpIOVarType
  use ConfigVarInitMod
  use EnergyVarInitMod
  use ForcingVarInitMod
  use WaterVarInitMod
  use BiochemVarInitMod
  use ConfigVarInTransferMod
  use EnergyVarInTransferMod
  use ForcingVarInTransferMod
  use WaterVarInTransferMod
  use BiochemVarInTransferMod
  use ConfigVarOutTransferMod
  use ForcingVarOutTransferMod
  use EnergyVarOutTransferMod
  use WaterVarOutTransferMod
  use BiochemVarOutTransferMod
  use NoahmpMainMod
  use NoahmpMainGlacierMod
  use NoahmpFatalMod, only: NoahmpIO_abort

  implicit none

  ! Days per month, non-leap. Shared so CAL_MON_DAY and NoahmpCalendarAdvance cannot
  ! drift; February is adjusted per-call against NoahmpYearLength, never in place.
  integer, parameter, private :: MONTH_DAYS(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  
contains  

  subroutine NoahmpDriverMain(NoahmpIO)
  
! Code history: original Noah-MP subroutine noahmplsm (Niu et al. 2011);
! refactored by C. He, P. Valayamkunnath & team (He et al. 2023)
 
    implicit none 
    
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO
    
    ! local variables
    type(noahmp_type)                   :: noahmp
    integer                             :: I
    integer                             :: J
    integer                             :: K
    integer                             :: JMONTH, JDAY
    integer                             :: CAL_YR
    real(kind=kind_noahmp)              :: SOLAR_TIME
    real(kind=kind_noahmp)              :: CAL_JULIAN

      ! ERF provides one atmospheric level; WRF physics expects two -- duplicate layer 1 into layer 2.
      NoahmpIO%P8W(:, 2, :) = NoahmpIO%P8W(:, 1, :)
      NoahmpIO%T_PHY(:, 2, :) = NoahmpIO%T_PHY(:, 1, :)
      NoahmpIO%U_PHY(:, 2, :) = NoahmpIO%U_PHY(:, 1, :)
      NoahmpIO%V_PHY(:, 2, :) = NoahmpIO%V_PHY(:, 1, :)
      NoahmpIO%QV_CURR(:, 2, :) = NoahmpIO%QV_CURR(:, 1, :)
      NoahmpIO%DZ8W(:, 2, :) = NoahmpIO%DZ8W(:, 1, :)
      ! Zero the channels ERF lacks (convective, shallow) and unused SNOWBL.
      NoahmpIO%SNOWBL  = 0.0
      NoahmpIO%RAINCV  = 0.0
      NoahmpIO%RAINSHV = 0.0

      ! The host (ERF) stages DZ8W per column, per step, from its own MOST reference
      ! height -- see stage_forcing in ERF_NOAHMP_Advance.cpp. Only fall back to the
      ! ZLVL broadcast (2x for the WRF model-level convention) if nothing was staged,
      ! which should not happen when driven from ERF.
      if (any(NoahmpIO%DZ8W == undefined_real)) then
         NoahmpIO%DZ8W = 2*NoahmpIO%ZLVL
         if (NoahmpIO%rank == 0) write(*,'(" ***** Noah-MP: host did not stage DZ8W; using ZLVL ",F0.3," m.")') NoahmpIO%ZLVL
      endif

      NoahmpIO%SWDDIR = NoahmpIO%SWDOWN*0.7                    ! following noahmplsm ATM 70% direct radiation
      NoahmpIO%SWDDIF = NoahmpIO%SWDOWN*0.3                    ! following noahmplsm ATM 30% diffuse radiation

      IF (NoahmpIO%ITIMESTEP == 1) THEN
         if (NoahmpIO%rank == 0) write(*,'("Noah-MP ITIMESTEP == 1 setting initial guess for variables")')
         NoahmpIO%EAHXY = (NoahmpIO%P8W(:, 1, :)*NoahmpIO%QV_CURR(:, 1, :))/(0.622+NoahmpIO%QV_CURR(:, 1, :)) ! Initial guess only.
         NoahmpIO%TAHXY = NoahmpIO%T_PHY(:, 1, :)                                                         ! Initial guess only.
         NoahmpIO%CHXY = 0.1
         NoahmpIO%CMXY = 0.1
      END IF

      IF (NoahmpIO%ITIMESTEP > 0) THEN
         if (NoahmpIO%rank == 0) write(*,'("Noah-MP running physical processes")')
         ! Convective/shallow absent in ERF; MP_HAIL is ERF-supplied (not set here).
         NoahmpIO%MP_RAINC = NoahmpIO%RAINCV
         NoahmpIO%MP_SHCV = NoahmpIO%RAINSHV

    !  Treatment of Noah-MP soil timestep
    NoahmpIO%CALCULATE_SOIL    = .false.
    NoahmpIO%SOIL_UPDATE_STEPS = nint(NoahmpIO%SOILTSTEP / NoahmpIO%DTBL)
    NoahmpIO%SOIL_UPDATE_STEPS = max(NoahmpIO%SOIL_UPDATE_STEPS,1)

    if ( NoahmpIO%SOIL_UPDATE_STEPS == 1 ) then
       NoahmpIO%ACC_SSOILXY  = 0.0
       NoahmpIO%ACC_QINSURXY = 0.0
       NoahmpIO%ACC_QSEVAXY  = 0.0
       NoahmpIO%ACC_ETRANIXY = 0.0
       NoahmpIO%ACC_DWATERXY = 0.0
       NoahmpIO%ACC_PRCPXY   = 0.0
       NoahmpIO%ACC_ECANXY   = 0.0
       NoahmpIO%ACC_ETRANXY  = 0.0
       NoahmpIO%ACC_EDIRXY   = 0.0
       NoahmpIO%ACC_GLAFLWXY = 0.0
    endif

    if ( NoahmpIO%SOIL_UPDATE_STEPS > 1 ) then
       if ( mod(NoahmpIO%ITIMESTEP, NoahmpIO%SOIL_UPDATE_STEPS) == 1 ) then
          NoahmpIO%ACC_SSOILXY  = 0.0
          NoahmpIO%ACC_QINSURXY = 0.0
          NoahmpIO%ACC_QSEVAXY  = 0.0
          NoahmpIO%ACC_ETRANIXY = 0.0
          NoahmpIO%ACC_DWATERXY = 0.0
          NoahmpIO%ACC_PRCPXY   = 0.0
          NoahmpIO%ACC_ECANXY   = 0.0
          NoahmpIO%ACC_ETRANXY  = 0.0
          NoahmpIO%ACC_EDIRXY   = 0.0
          NoahmpIO%ACC_GLAFLWXY = 0.0
       end if
    endif

    ! Set directly (not in the if above) to avoid stale calculate_soil across cpu threads
    NoahmpIO%CALCULATE_SOIL = mod(NoahmpIO%ITIMESTEP, NoahmpIO%SOIL_UPDATE_STEPS) == 0

    !  Prepare Noah-MP driver

    ! Wall-clock date of this call: the namelist start date advanced by the elapsed
    ! model time. ITIMESTEP is 1-based, so firing i sits at (i-1)*DTBL seconds.
    call NoahmpCalendarAdvance(NoahmpIO, real(NoahmpIO%ITIMESTEP-1, kind=kind_noahmp) * NoahmpIO%DTBL, &
                               CAL_YR, CAL_JULIAN)
    NoahmpIO%YR     = CAL_YR
    NoahmpIO%JULIAN = CAL_JULIAN

    ! Announce the clock once. In the ERF driver the calendar start now comes
    ! directly from the wrfinput/WPS file header.
    if ( (NoahmpIO%ITIMESTEP == 1) .and. (NoahmpIO%rank == 0) ) then
       write(*,'(" ***** Noah-MP calendar start (from wrfinput/WPS file): ",       &
                 &I0,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2," UTC")')                   &
             NoahmpIO%start_year, NoahmpIO%start_month, NoahmpIO%start_day,          &
             max(NoahmpIO%start_hour, 0), max(NoahmpIO%start_min, 0)
    endif

    ! find length of year for phenology (also S Hemisphere)
    NoahmpIO%YEARLEN = NoahmpYearLength(NoahmpIO%YR)

    ! depth to soil interfaces (<0) [m]
    NoahmpIO%ZSOIL(1) = -NoahmpIO%DZS(1)
    do K = 2, NoahmpIO%NSOIL
       NoahmpIO%ZSOIL(K) = -NoahmpIO%DZS(K) + NoahmpIO%ZSOIL(K-1)
    enddo
    
    JLOOP : do J = NoahmpIO%JTS, NoahmpIO%JTE

       NoahmpIO%J = J
       if ( NoahmpIO%ITIMESTEP == 1 ) then
          do I = NoahmpIO%ITS, NoahmpIO%ITE
             if ( (NoahmpIO%XLAND(I,J)-1.5) >= 0.0 ) then  ! Open water point
                if ( (NoahmpIO%XICE(I,J) == 1.0) .and. (NoahmpIO%rank == 0) )                &
                   write(*,'(" ***** Noah-MP: sea-ice at water point, I=",I0," J=",I0)') I, J
                NoahmpIO%SMSTAV(I,J) = 1.0
                NoahmpIO%SMSTOT(I,J) = 1.0
                do K = 1, NoahmpIO%NSOIL
                   NoahmpIO%SMOIS(I,K,J) = 1.0
                   NoahmpIO%TSLB(I,K,J)  = 273.16
                enddo
             else
                if ( NoahmpIO%XICE(I,J) == 1.0 ) then      ! Sea-ice case
                   NoahmpIO%SMSTAV(I,J) = 1.0
                   NoahmpIO%SMSTOT(I,J) = 1.0
                   do K = 1, NoahmpIO%NSOIL
                      NoahmpIO%SMOIS(I,K,J) = 1.0
                   enddo
                endif
             endif
          enddo
       endif  ! end of initialization over ocean

       ILOOP : do I = NoahmpIO%ITS, NoahmpIO%ITE

          NoahmpIO%I = I
          if ( NoahmpIO%XICE(I,J) >= NoahmpIO%XICE_THRESHOLD ) then  ! Sea-ice point
             NoahmpIO%ICE                        = 1
             NoahmpIO%SH2O(I,1:NoahmpIO%NSOIL,J) = 1.0
             NoahmpIO%LAI (I,J)                  = 0.01
             cycle ILOOP                                             ! Skip any sea-ice points
          else
             if ( (NoahmpIO%XLAND(I,J)-1.5) >= 0.0 ) cycle ILOOP     ! Skip any open water points

             ! ICE is a scalar carried across iterations, so without this a land point
             ! after a sea-ice point inherits ICE=1 (and the first sees undefined_int).
             NoahmpIO%ICE = 0                                        ! Non-ice land point

             !  Initialize data types and transfer inputs from 2-D to 1-D column variables
             call ConfigVarInitDefault  (noahmp)
             call ConfigVarInTransfer   (noahmp, NoahmpIO)
             call ForcingVarInitDefault (noahmp)
             call ForcingVarInTransfer  (noahmp, NoahmpIO)
             call EnergyVarInitDefault  (noahmp)
             call EnergyVarInTransfer   (noahmp, NoahmpIO)
             call WaterVarInitDefault   (noahmp)
             call WaterVarInTransfer    (noahmp, NoahmpIO)
             call BiochemVarInitDefault (noahmp)
             call BiochemVarInTransfer  (noahmp, NoahmpIO)

             !  Urban vegetation hydrology: irrigate only in urban area, MAY-SEP, 9-11pm.
             !  TODO: separate urban-specific logic out of the Noah-MP driver.
             if ( (NoahmpIO%IVGTYP(I,J) == NoahmpIO%ISURBAN_TABLE) .or. &
                  (NoahmpIO%IVGTYP(I,J) > NoahmpIO%URBTYPE_beg) ) then
                if ( (NoahmpIO%SF_URBAN_PHYSICS > 0) .and. (NoahmpIO%IRI_URBAN == 1) ) then
                   SOLAR_TIME = (NoahmpIO%JULIAN - int(NoahmpIO%JULIAN))*24 + NoahmpIO%XLONG(I,J)/15.0
                   if ( SOLAR_TIME < 0.0 ) SOLAR_TIME = SOLAR_TIME + 24.0
                   ! JULIAN is 0-based; CAL_MON_DAY takes 1..YEARLEN.
                   call CAL_MON_DAY(int(NoahmpIO%JULIAN)+1, NoahmpIO%YR, JMONTH, JDAY)
                   if ( (SOLAR_TIME >= 21.0) .and. (SOLAR_TIME <= 23.0) .and. &
                        (JMONTH >= 5) .and. (JMONTH <= 9) ) then
                       noahmp%water%state%SoilMoisture(1) = &
                              max(noahmp%water%state%SoilMoisture(1),noahmp%water%param%SoilMoistureFieldCap(1))
                       noahmp%water%state%SoilMoisture(2) = &
                              max(noahmp%water%state%SoilMoisture(2),noahmp%water%param%SoilMoistureFieldCap(2))
                   endif
                endif
             endif

             !  Call 1D Noah-MP LSM

             ! glacier ice
             if (noahmp%config%domain%VegType == noahmp%config%domain%IndexIcePoint ) then
                 noahmp%config%domain%IndicatorIceSfc = -1  ! Land-ice point      
                 noahmp%forcing%TemperatureSoilBottom = min(noahmp%forcing%TemperatureSoilBottom,263.15) ! set deep glaicer temp to >= -10C
                 call NoahmpMainGlacier(noahmp)
             ! non-glacier land
             else
                 noahmp%config%domain%IndicatorIceSfc = 0   ! land soil point.
                 call NoahmpMain(noahmp)
             endif ! glacial split ends

             !  Transfer 1-D Noah-MP column variables to 2-D output variables
             call ConfigVarOutTransfer (noahmp, NoahmpIO)
             call ForcingVarOutTransfer(noahmp, NoahmpIO)
             call EnergyVarOutTransfer (noahmp, NoahmpIO)
             call WaterVarOutTransfer  (noahmp, NoahmpIO)
             call BiochemVarOutTransfer(noahmp, NoahmpIO) 

          endif     ! land-sea split ends

       enddo ILOOP  ! I loop
    enddo  JLOOP    ! J loop
    end if
 
  end subroutine NoahmpDriverMain

  subroutine CAL_MON_DAY(JULDAY, julyr, Jmonth, Jday)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JULDAY, julyr
      INTEGER, INTENT(OUT) :: Jmonth, Jday
      LOGICAL :: NOT_FIND_DATE
      INTEGER :: MONTH(12), itmpday, i
      NOT_FIND_DATE = .true.

      ! Re-seed every call: a DATA-initialized MONTH is implicitly SAVEd, so one
      ! leap-year call used to leave MONTH(2)=29 set for the rest of the run.
      MONTH   = MONTH_DAYS
      itmpday = JULDAY
      IF (NoahmpYearLength(julyr) .EQ. 366) MONTH(2) = 29

      do i = 1, 12
         IF (itmpday .GT. MONTH(i)) THEN
            itmpday = itmpday-MONTH(i)
         ELSE
            Jday = itmpday
            Jmonth = i
            NOT_FIND_DATE = .false.
            exit
         END IF
      end do
      if (NOT_FIND_DATE) then
         write(*,'(" ***** CAL_MON_DAY: day-of-year ",I0," out of range 1..",I0,          &
                   &" for year ",I0)') JULDAY, NoahmpYearLength(julyr), julyr
         call NoahmpIO_abort()
      endif

   end subroutine CAL_MON_DAY

   ! ---------------------------------------------------------------------------
   !  Calendar. The host owns ITIMESTEP; the wall-clock date follows from the
   !  namelist start date plus the elapsed model time. Before this existed, YR and
   !  JULIAN sat frozen at their NoahmpIOVarInitMod defaults (2000, Jan 1) for the
   !  whole run, silently pinning phenology and the urban irrigation window.
   ! ---------------------------------------------------------------------------

   ! Days in a year, proleptic Gregorian.
   pure integer function NoahmpYearLength(year) result(ylen)

      implicit none
      integer, intent(in) :: year

      ylen = 365
      if (mod(year,4) == 0) then
         ylen = 366
         if (mod(year,100) == 0) then
            ylen = 365
            if (mod(year,400) == 0) ylen = 366
         endif
      endif

   end function NoahmpYearLength

   ! Advance the namelist start date by elapsed_sec. Returns the calendar year and
   ! the day-of-year carrying the fraction of the day, 0-based (Jan 1 00Z -> 0.0)
   ! so 0 <= JULIAN < YEARLEN. That is the convention DayJulianInYear is consumed
   ! in: PhenologyMainMod interpolates the monthly table LAI/SAI with
   ! 12*DayCurrent/NumDayInYear, and the crop and irrigation date checks compare
   ! against it. Consumers needing a 1-based day convert at the call site.
   !
   ! NoahmpIO is read-only here: the caller assigns the results, so YR/JULIAN are
   ! never aliased against the intent(in) dummy.
   subroutine NoahmpCalendarAdvance(NoahmpIO, elapsed_sec, YR, JULIAN)

      implicit none
      type(NoahmpIO_type),    intent(in)  :: NoahmpIO
      real(kind=kind_noahmp), intent(in)  :: elapsed_sec   ! since the start date [s]
      integer,                intent(out) :: YR            ! 4-digit calendar year
      real(kind=kind_noahmp), intent(out) :: JULIAN        ! 0-based day-of-year + day fraction

      integer, parameter     :: i8 = selected_int_kind(18)
      integer                :: hh, mm, imon, doy, ylen, dmon
      integer(kind=i8)       :: isec, sec_of_year, sec_in_year
      real(kind=kind_noahmp) :: fsec

      ! Only start_year/month/day are mandatory in namelist.erf, so an unset
      ! hour/minute means midnight rather than the -9999 sentinel.
      hh = NoahmpIO%start_hour
      mm = NoahmpIO%start_min
      if (hh == undefined_int) hh = 0
      if (mm == undefined_int) mm = 0

      YR   = NoahmpIO%start_year
      ylen = NoahmpYearLength(YR)

      dmon = 0
      if ( (NoahmpIO%start_month >= 1) .and. (NoahmpIO%start_month <= 12) ) then
         dmon = MONTH_DAYS(NoahmpIO%start_month)
         if ( (NoahmpIO%start_month == 2) .and. (ylen == 366) ) dmon = 29
      endif

      if ( (YR <= 0) .or. (NoahmpIO%start_month < 1) .or. (NoahmpIO%start_month > 12) .or. &
           (NoahmpIO%start_day < 1) .or. (NoahmpIO%start_day > dmon) .or.                  &
           (hh < 0) .or. (hh > 23) .or. (mm < 0) .or. (mm > 59) ) then
         if (NoahmpIO%rank == 0) then
            write(*,'(" ***** Noah-MP: namelist.erf start date is not a valid date: ",     &
                      &"year=",I0," month=",I0," day=",I0," hour=",I0," min=",I0)')        &
                  NoahmpIO%start_year, NoahmpIO%start_month, NoahmpIO%start_day, hh, mm
         endif
         call NoahmpIO_abort()
      endif

      ! Day-of-year of the start date (1-based).
      doy = NoahmpIO%start_day
      do imon = 1, NoahmpIO%start_month - 1
         if ( (imon == 2) .and. (ylen == 366) ) then
            doy = doy + 29
         else
            doy = doy + MONTH_DAYS(imon)
         endif
      end do

      ! Carry the rollover in whole seconds so it stays exact in a single-precision
      ! build; only the sub-second remainder goes through the real accumulator.
      isec        = floor(elapsed_sec, kind=i8)
      fsec        = elapsed_sec - real(isec, kind=kind_noahmp)
      sec_of_year = int(doy-1,i8)*86400_i8 + int(hh,i8)*3600_i8 + int(mm,i8)*60_i8 + isec

      do
         sec_in_year = int(NoahmpYearLength(YR),i8) * 86400_i8
         if (sec_of_year < sec_in_year) exit
         sec_of_year = sec_of_year - sec_in_year
         YR = YR + 1
      end do
      do while (sec_of_year < 0_i8)
         YR = YR - 1
         sec_of_year = sec_of_year + int(NoahmpYearLength(YR),i8)*86400_i8
      end do

      JULIAN = (real(sec_of_year, kind=kind_noahmp) + fsec)          &
             / 86400.0_kind_noahmp

   end subroutine NoahmpCalendarAdvance

end module NoahmpDriverMainMod  
