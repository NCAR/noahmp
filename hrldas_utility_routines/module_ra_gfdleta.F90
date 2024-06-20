!
! An attempt to mimic the GFDL code without requiring all the unnecessary extras.


MODULE MODULE_RA_GFDLETA

    contains

    !---------------------------------------------------------------------
        SUBROUTINE CAL_MON_DAY(JULDAY,julyr,Jmonth,Jday)
    !---------------------------------------------------------------------
        IMPLICIT NONE
    !-----------------------------------------------------------------------
        INTEGER, INTENT(IN) :: JULDAY,julyr
        INTEGER, INTENT(OUT) :: Jmonth,Jday
        LOGICAL :: LEAP,NOT_FIND_DATE
        INTEGER :: MONTH (12),itmpday,itmpmon,i
    !-----------------------------------------------------------------------
        DATA MONTH/31,28,31,30,31,30,31,31,30,31,30,31/
    !***********************************************************************
        NOT_FIND_DATE = .true.

        itmpday = JULDAY
        itmpmon = 1
        LEAP=.FALSE.
        IF(MOD(julyr,4).EQ.0)THEN
          MONTH(2)=29
          LEAP=.TRUE.
        ENDIF

        i = 1
        DO WHILE (NOT_FIND_DATE)
           IF(itmpday.GT.MONTH(i))THEN
             itmpday=itmpday-MONTH(i)
           ELSE
             Jday=itmpday
             Jmonth=i
             NOT_FIND_DATE = .false.
           ENDIF
           i = i+1
        END DO

        END SUBROUTINE CAL_MON_DAY
    !!================================================================================

          END MODULE module_RA_GFDLETA
