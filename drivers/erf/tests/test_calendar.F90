! ===========================================================================
! test_calendar -- NoahmpCalendarAdvance's day-of-year convention.
!
! DayJulianInYear is 0-based: PhenologyMainMod interpolates the monthly table
! LAI/SAI with 12*DayCurrent/NumDayInYear over 0<=DayCurrent<NumDayInYear, and
! the crop and irrigation date checks compare against it. A 1-based day reads
! the table one day late and pushes MonthCurrent past 12 on 31 December.
!
!   * Phase A -- the convention: Jan 1 00Z is 0.0, the day fraction rides along,
!     and 0 <= JULIAN < YEARLEN including on 31 December.
!   * Phase B -- CAL_MON_DAY takes a 1-based day, so int(JULIAN)+1 round-trips
!     to the right month and day (leap day included).
!   * Phase C -- the phenology month interpolation reproduces the
!     hand-interpolated table value for a known date.
!
! Pure Fortran, no NetCDF, no domain fixture -- fits the Tier-1 test style.
! Exit code 0 = pass.
! ===========================================================================
program test_calendar

  use Machine,             only : kind_noahmp
  use NoahmpIOVarType,     only : NoahmpIO_type
  use NoahmpDriverMainMod, only : NoahmpCalendarAdvance, NoahmpYearLength, CAL_MON_DAY

  implicit none

  type(NoahmpIO_type) :: blk
  integer             :: nfail
  real(kind=kind_noahmp) :: jul
  integer                :: yr

  real(kind=kind_noahmp), parameter :: DAY = 86400.0_kind_noahmp
  real(kind=kind_noahmp), parameter :: HOUR = 3600.0_kind_noahmp

  nfail = 0

  ! Only the error path reads rank, but it is a pointer component: give it a
  ! target so an unexpected abort prints rather than dereferencing null.
  allocate(blk%rank); blk%rank = 0

  ! =====================================================================
  ! Phase A -- the 0-based convention
  ! =====================================================================

  ! Jan 1 00Z is day 0.0, not 1.0.
  call set_start(2023, 1, 1, 0, 0)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_int (yr,  2023,              "Jan 1 00Z year")
  call expect_close(jul, 0.0_kind_noahmp,  "Jan 1 00Z is day 0")

  ! The fraction of the day rides on top of the integer day.
  call advance(12.0_kind_noahmp*HOUR, yr, jul)
  call expect_close(jul, 0.5_kind_noahmp,  "Jan 1 12Z is day 0.5")

  ! An unset start hour/minute means midnight, so a start time also lands on .0.
  call set_start(2023, 1, 2, 0, 0)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_close(jul, 1.0_kind_noahmp,  "Jan 2 00Z is day 1")

  ! A non-leap year: Aug 5 is the 217th day counting from 1, so 216 from 0.
  call set_start(2023, 8, 5, 0, 0)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_close(jul, 216.0_kind_noahmp, "2023-08-05 00Z (non-leap)")

  ! A leap year shifts everything after February by one day.
  call set_start(2024, 8, 5, 0, 0)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_int  (NoahmpYearLength(2024), 366, "2024 is a leap year")
  call expect_close(jul, 217.0_kind_noahmp, "2024-08-05 00Z (leap)")

  ! The last day of a 365-day year stays below YEARLEN; a 1-based day reached
  ! 365.96 and pushed the phenology month past 12.
  call set_start(2023, 12, 31, 23, 0)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_close(jul, 364.0_kind_noahmp + 23.0_kind_noahmp/24.0_kind_noahmp, &
                    "2023-12-31 23Z")
  call expect_lt(jul, real(NoahmpYearLength(2023), kind_noahmp), &
                 "31 December stays below YEARLEN")
  call expect_lt(12.0_kind_noahmp*jul/real(NoahmpYearLength(2023), kind_noahmp), &
                 12.0_kind_noahmp, "phenology MonthCurrent stays below 12")

  ! Rolling forward past 31 December resets the day to 0 in the new year.
  call set_start(2023, 12, 31, 0, 0)
  call advance(DAY, yr, jul)
  call expect_int  (yr, 2024,              "rollover year")
  call expect_close(jul, 0.0_kind_noahmp,  "rollover lands on day 0")

  ! Rolling backward lands on the last day of the previous year, not past it.
  call set_start(2023, 1, 1, 0, 0)
  call advance(-DAY, yr, jul)
  call expect_int  (yr, 2022,              "backward rollover year")
  call expect_close(jul, 364.0_kind_noahmp, "backward rollover is day 364")

  ! =====================================================================
  ! Phase B -- CAL_MON_DAY takes a 1-based day: int(JULIAN)+1
  ! =====================================================================

  call set_start(2023, 1, 1, 0, 0)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_date(jul, yr, 1, 1, "Jan 1 round-trips")

  call set_start(2023, 8, 5, 6, 30)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_date(jul, yr, 8, 5, "Aug 5 round-trips (non-leap)")

  ! The leap day, and the day after it.
  call set_start(2024, 2, 28, 0, 0)
  call advance(DAY, yr, jul)
  call expect_date(jul, yr, 2, 29, "Feb 29 round-trips")
  call advance(2.0_kind_noahmp*DAY, yr, jul)
  call expect_date(jul, yr, 3, 1, "Mar 1 after the leap day")

  ! int(JULIAN)+1 must equal YEARLEN, the largest day CAL_MON_DAY accepts.
  call set_start(2023, 12, 31, 23, 59)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_int  (int(jul)+1, NoahmpYearLength(2023), "last day is YEARLEN")
  call expect_date(jul, yr, 12, 31, "Dec 31 round-trips")

  ! =====================================================================
  ! Phase C -- the phenology month interpolation on the table LAI
  ! =====================================================================
  !
  ! MODIS grassland (class 10) has a monthly table LAI of 3.5 in July and 1.5 in
  ! August, and on 2024-08-05 00Z the two are weighted by 12*DayCurrent/366:
  !
  !   0-based day 217 -> LAI 2.2705
  !   1-based day 218 -> LAI 2.2049    (the table of one day later)
  call set_start(2024, 8, 5, 0, 0)
  call advance(0.0_kind_noahmp, yr, jul)
  call expect_close(table_lai(jul, NoahmpYearLength(yr)), 2.2705_kind_noahmp, &
                    "grassland table LAI on 2024-08-05")

  if (nfail == 0) then
     write(*,'(A)') "PASS: test_calendar"
     call exit(0)
  else
     write(*,'(A,I0,A)') "FAILED: test_calendar (", nfail, " check(s))"
     call exit(1)
  end if

contains

  subroutine set_start(y, mo, d, h, mi)
    integer, intent(in) :: y, mo, d, h, mi
    blk%start_year  = y
    blk%start_month = mo
    blk%start_day   = d
    blk%start_hour  = h
    blk%start_min   = mi
  end subroutine set_start

  subroutine advance(elapsed_sec, y, j)
    real(kind=kind_noahmp), intent(in)  :: elapsed_sec
    integer,                intent(out) :: y
    real(kind=kind_noahmp), intent(out) :: j
    call NoahmpCalendarAdvance(blk, elapsed_sec, y, j)
  end subroutine advance

  ! PhenologyMainMod's monthly interpolation, Northern Hemisphere, apart from
  ! the hard-wired grassland row. IntpMonth* are integers there too; the int()
  ! is that same truncation, written out.
  function table_lai(j, yearlen) result(lai)
    real(kind=kind_noahmp), intent(in) :: j
    integer,                intent(in) :: yearlen
    real(kind=kind_noahmp) :: lai, MonthCurrent, IntpWgt1, IntpWgt2
    integer                :: IntpMonth1, IntpMonth2
    ! MODIS class 10, grassland: monthly table LeafAreaIndex.
    real(kind=kind_noahmp), parameter :: LAIM(12) = (/ &
         0.4_kind_noahmp, 0.5_kind_noahmp, 0.6_kind_noahmp, 0.7_kind_noahmp, &
         1.2_kind_noahmp, 3.0_kind_noahmp, 3.5_kind_noahmp, 1.5_kind_noahmp, &
         0.7_kind_noahmp, 0.6_kind_noahmp, 0.5_kind_noahmp, 0.4_kind_noahmp /)

    MonthCurrent = 12.0_kind_noahmp * j / real(yearlen, kind_noahmp)
    IntpMonth1   = int(MonthCurrent + 0.5_kind_noahmp)
    IntpMonth2   = IntpMonth1 + 1
    IntpWgt1     = (IntpMonth1 + 0.5_kind_noahmp) - MonthCurrent
    IntpWgt2     = 1.0_kind_noahmp - IntpWgt1
    if ( IntpMonth1 <  1 ) IntpMonth1 = 12
    if ( IntpMonth2 > 12 ) IntpMonth2 = 1
    lai = IntpWgt1 * LAIM(IntpMonth1) + IntpWgt2 * LAIM(IntpMonth2)
  end function table_lai

  ! CAL_MON_DAY takes a 1-based day, as int(JULIAN)+1.
  subroutine expect_date(j, y, want_mon, want_day, name)
    real(kind=kind_noahmp), intent(in) :: j
    integer,                intent(in) :: y, want_mon, want_day
    character(*),           intent(in) :: name
    integer :: jmonth, jday
    call CAL_MON_DAY(int(j)+1, y, jmonth, jday)
    if ( (jmonth /= want_mon) .or. (jday /= want_day) ) then
       write(0,'(A,A,A,I0,A,I0,A,I0,A,I0)') "  FAIL: ", name, " got=", jmonth, "/", jday, &
             " want=", want_mon, "/", want_day
       nfail = nfail + 1
    end if
  end subroutine expect_date

  subroutine expect_int(got, want, name)
    integer,      intent(in) :: got, want
    character(*), intent(in) :: name
    if (got /= want) then
       write(0,'(A,A,A,I0,A,I0)') "  FAIL: ", name, " got=", got, " want=", want
       nfail = nfail + 1
    end if
  end subroutine expect_int

  subroutine expect_lt(got, bound, name)
    real(kind=kind_noahmp), intent(in) :: got, bound
    character(*),           intent(in) :: name
    if (.not. (got < bound)) then
       write(0,'(A,A,A,ES14.6,A,ES14.6)') "  FAIL: ", name, " got=", got, " not < ", bound
       nfail = nfail + 1
    end if
  end subroutine expect_lt

  subroutine expect_close(got, want, name)
    real(kind=kind_noahmp), intent(in) :: got, want
    character(*),           intent(in) :: name
    real(kind=kind_noahmp) :: tol
    tol = 1.0e-4_kind_noahmp * max(1.0_kind_noahmp, abs(want))
    if (abs(got - want) > tol) then
       write(0,'(A,A,A,ES14.6,A,ES14.6)') "  FAIL: ", name, " got=", got, " want=", want
       nfail = nfail + 1
    end if
  end subroutine expect_close

end program test_calendar
