module lnd_comp_shr

  ! Coupling related shared routines
  use ESMF,  only : operator(-), operator(+), operator(<=), operator(*)
  use ESMF,  only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_FieldDestroy
  use ESMF,  only : ESMF_Mesh, ESMF_MeshGet, ESMF_MeshCreate, ESMF_VM
  use ESMF,  only : ESMF_FILEFORMAT_ESMFMESH, ESMF_SUCCESS, ESMF_REDUCE_SUM
  use ESMF,  only : ESMF_DistGrid, ESMF_DistGridGet, ESMF_FieldGet
  use ESMF,  only : ESMF_FieldRegrid, ESMF_REGION_TOTAL, ESMF_TERMORDER_SRCSEQ
  use ESMF,  only : ESMF_VMReduce, ESMF_VMAllReduce, ESMF_VMBroadcast
  use ESMF,  only : ESMF_RouteHandle, ESMF_Field, ESMF_Array, ESMF_ArrayCreate
  use ESMF,  only : ESMF_FieldCreate, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
  use ESMF,  only : ESMF_FieldRegridStore, ESMF_REGRIDMETHOD_CONSERVE
  use ESMF,  only : ESMF_NORMTYPE_DSTAREA, ESMF_UNMAPPEDACTION_IGNORE
  use ESMF,  only : ESMF_MeshSet, ESMF_FAILURE, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF,  only : ESMF_VMGet, ESMF_VMGetCurrent, ESMF_State, ESMF_StateGet
  use ESMF,  only : ESMF_DistGridGet, ESMF_Clock, ESMF_Alarm, ESMF_Time
  use ESMF,  only : ESMF_TimeInterval, ESMF_ClockGet, ESMF_TimeGet
  use ESMF,  only : ESMF_TimeIntervalSet, ESMF_TimeSet, ESMF_Calendar
  use ESMF,  only : ESMF_AlarmCreate, ESMF_LOGMSG_ERROR, ESMF_GridComp
  use ESMF,  only : ESMF_GridCompGet, ESMF_Grid, ESMF_GridCreateMosaic
  use ESMF,  only : ESMF_GridDestroy, ESMF_GridAddItem, ESMF_GridGetItem
  use ESMF,  only : ESMF_KIND_I4, ESMF_GRIDITEM_MASK, ESMF_STAGGERLOC_CENTER
  use ESMF,  only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
  use ESMF,  only : ESMF_FieldWriteVTK, ESMF_FieldWrite, ESMF_Decomp_Flag, ESMF_DECOMP_SYMMEDGEMAX
  use ESMF,  only : ESMF_INDEX_GLOBAL, ESMF_KIND_R4, ESMF_Field, ESMF_FieldGet
  use ESMF,  only : ESMF_LogFoundNetCDFError, ESMF_ArraySpec, ESMF_ArraySpecSet, ESMF_TYPEKIND_R4
  use ESMF,  only : ESMF_RouteHandle, ESMF_FieldRegridStore, ESMF_FieldRedist
  use NUOPC, only : NUOPC_CompAttributeGet

  use lnd_comp_kind  , only : r8 => shr_kind_r8 
  use lnd_comp_kind  , only : cl => shr_kind_cl
  use lnd_comp_types , only : noahmp_type

  implicit none
  private

  public :: alarm_init
  public :: chkerr 
  public :: chkerrnc
  public :: read_namelist
  public :: shr_string_listGetName
  public :: shr_string_listGetNum
  public :: shr_string_countChar
  public :: time_init

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  ! clock and alarm options
  character(len=*), private, parameter :: &
       optNONE           = "none"      , &
       optNever          = "never"     , &
       optNSteps         = "nsteps"    , &
       optNStep          = "nstep"     , &
       optNSeconds       = "nseconds"  , &
       optNSecond        = "nsecond"   , &
       optNMinutes       = "nminutes"  , &
       optNMinute        = "nminute"   , &
       optNHours         = "nhours"    , &
       optNHour          = "nhour"     , &
       optNDays          = "ndays"     , &
       optNDay           = "nday"      , &
       optNMonths        = "nmonths"   , &
       optNMonth         = "nmonth"    , &
       optNYears         = "nyears"    , &
       optNYear          = "nyear"     , &
       optMonthly        = "monthly"   , &
       optYearly         = "yearly"    , &
       optEnd            = "end"       , &
       optDate           = "date"

  integer :: master_task = 0
  integer :: dbug = 1
  integer :: iswater = 17
  character(len=1), save  :: listDel  = ":"
  character(*), parameter :: modName =  "(comp_shr)"

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  logical function ChkErr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    ChkErr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       ChkErr = .true.
    endif
  end function ChkErr

  !===============================================================================
  logical function ChkErrNc(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    ChkErrNc = .false.
    lrc = rc
    if (ESMF_LogFoundNetCDFError(lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       ChkErrNc = .true.
    endif
  end function ChkErrNc

  !===============================================================================
  subroutine alarm_init( clock, alarm, option, &
       opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)

    ! Setup an alarm in a clock
    ! Notes: The ringtime sent to AlarmCreate MUST be the next alarm
    ! time.  If you send an arbitrary but proper ringtime from the
    ! past and the ring interval, the alarm will always go off on the
    ! next clock advance and this will cause serious problems.  Even
    ! if it makes sense to initialize an alarm with some reference
    ! time and the alarm interval, that reference time has to be
    ! advance forward to be >= the current time.  In the logic below
    ! we set an appropriate "NextAlarm" and then we make sure to
    ! advance it properly based on the ring interval.

    ! input/output variables
    type(ESMF_Clock)            , intent(inout) :: clock     ! clock
    type(ESMF_Alarm)            , intent(inout) :: alarm     ! alarm
    character(len=*)            , intent(in)    :: option    ! alarm option
    integer          , optional , intent(in)    :: opt_n     ! alarm freq
    integer          , optional , intent(in)    :: opt_ymd   ! alarm ymd
    integer          , optional , intent(in)    :: opt_tod   ! alarm tod (sec)
    type(ESMF_Time)  , optional , intent(in)    :: RefTime   ! ref time
    character(len=*) , optional , intent(in)    :: alarmname ! alarm name
    integer                     , intent(inout) :: rc        ! Return code

    ! local variables
    type(ESMF_Calendar)     :: cal              ! calendar
    integer                 :: lymd             ! local ymd
    integer                 :: ltod             ! local tod
    integer                 :: cyy,cmm,cdd,csec ! time info
    character(len=64)       :: lalarmname       ! local alarm name
    logical                 :: update_nextalarm ! update next alarm
    type(ESMF_Time)         :: CurrTime         ! Current Time
    type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    integer                 :: sec
    character(len=*), parameter :: subname = '(set_alarmInit): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lalarmname = 'alarm_unknown'
    if (present(alarmname)) lalarmname = trim(alarmname)
    ltod = 0
    if (present(opt_tod)) ltod = opt_tod
    lymd = -1
    if (present(opt_ymd)) lymd = opt_ymd

    call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initial guess of next alarm, this will be updated below
    if (present(RefTime)) then
       NextAlarm = RefTime
    else
       NextAlarm = CurrTime
    endif

    ! Determine calendar
    call ESMF_ClockGet(clock, calendar=cal)

    ! Determine inputs for call to create alarm
    selectcase (trim(option))

    case (optNONE)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optNever)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optEnd)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optDate)
       if (.not. present(opt_ymd)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_ymd", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (lymd < 0 .or. ltod < 0) then
          call ESMF_LogWrite(trim(subname)//": "//trim(option)//'opt_ymd, opt_tod invalid', ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call time_init(NextAlarm, lymd, cal, ltod, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optNSteps)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNStep)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSeconds)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSecond)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinutes)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinute)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHours)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHour)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDays)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDay)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonths)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonth)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optMonthly)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case (optNYears)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNYear)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(trim(subname)//": requires opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(trim(subname)//": invalid opt_n", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optYearly)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.
       
    case default
       call ESMF_LogWrite(trim(subname)//": unknown option "//trim(option), ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return

    end select

    ! --------------------------------------------------------------------------------
    ! --- AlarmInterval and NextAlarm should be set ---
    ! --------------------------------------------------------------------------------

    ! --- advance Next Alarm so it won't ring on first timestep for
    ! --- most options above. go back one alarminterval just to be careful

    if (update_nextalarm) then
       NextAlarm = NextAlarm - AlarmInterval
       do while (NextAlarm <= CurrTime)
          NextAlarm = NextAlarm + AlarmInterval
       enddo
    endif
    alarm = ESMF_AlarmCreate( name=lalarmname, clock=clock, ringTime=NextAlarm, &
         ringInterval=AlarmInterval, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine alarm_init

  !===============================================================================
  subroutine time_init( Time, ymd, cal, tod, rc)

    ! Create the ESMF_Time object corresponding to the given input time,
    ! given in YMD (Year Month Day) and TOD (Time-of-day) format.
    ! Set the time by an integer as YYYYMMDD and integer seconds in the day

    ! input/output parameters:
    type(ESMF_Time)     , intent(inout) :: Time ! ESMF time
    integer             , intent(in)    :: ymd  ! year, month, day YYYYMMDD
    type(ESMF_Calendar) , intent(in)    :: cal  ! ESMF calendar
    integer             , intent(in)    :: tod  ! time of day in seconds
    integer             , intent(out)   :: rc

    ! local variables
    integer :: year, mon, day ! year, month, day as integers
    integer :: tdate          ! temporary date
    character(len=*), parameter :: subname='(time_init)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if ( (ymd < 0) .or. (tod < 0) .or. (tod > 86400) )then
       call ESMF_LogWrite(trim(subname)//": yymmdd is a negative number or time-of-day out of bounds", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    tdate = abs(ymd)
    year = int(tdate/10000)
    if (ymd < 0) year = -year
    mon = int( mod(tdate,10000)/  100)
    day = mod(tdate,  100)

    call ESMF_TimeSet( Time, yy=year, mm=mon, dd=day, s=tod, calendar=cal, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine time_init

  !===============================================================================
  subroutine read_namelist(gcomp, noahmp, rc)

    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_GridComp), intent(in)    :: gcomp
    type(noahmp_type)  , intent(inout) :: noahmp 
    integer, intent(out) :: rc

    ! local variables
    integer :: n
    character(len=cl) :: cname
    character(len=cl) :: cvalue
    character(len=cl) :: msg 
    character(len=cl), allocatable :: valueList(:)
    logical :: isPresent, isSet
    character(len=*),parameter :: subname=trim(modName)//':(read_namelist) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! get case name
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%case_name
    else
       noahmp%nmlist%case_name = 'ufs'
    end if

    ! forcing height
    call NUOPC_CompAttributeGet(gcomp, name='forcing_height', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%forcing_height
    else
       noahmp%nmlist%forcing_height = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : forcing_height = '//trim(cvalue), ESMF_LOGMSG_INFO)

    ! get num_soil_levels
    call NUOPC_CompAttributeGet(gcomp, name='num_soil_levels', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%num_soil_levels
       ! allocate variables
       if (.not. allocated(noahmp%nmlist%soil_level_thickness)) then
          allocate(noahmp%nmlist%soil_level_thickness(noahmp%nmlist%num_soil_levels))
       end if
       if (.not. allocated(noahmp%nmlist%soil_level_nodes)) then
          allocate(noahmp%nmlist%soil_level_nodes(noahmp%nmlist%num_soil_levels))
       end if
    else
       noahmp%nmlist%num_soil_levels = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : num_soil_levels = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%km = noahmp%nmlist%num_soil_levels

    ! get soil_level_thickness
    call NUOPC_CompAttributeGet(gcomp, name='soil_level_thickness', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       do n = 1, noahmp%nmlist%num_soil_levels
          call shr_string_listGetName(cvalue, n, cname, rc) 
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          read(cname,*) noahmp%nmlist%soil_level_thickness(n)
          write(msg, fmt='(A,I2.2,A)') trim(subname)//' : soil_level_thickness(',n,') = '//trim(cname)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
        end do
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in soil_level_thickness", ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return      
    end if

    ! get soil_level_nodes
    call NUOPC_CompAttributeGet(gcomp, name='soil_level_nodes', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       do n = 1, noahmp%nmlist%num_soil_levels
          call shr_string_listGetName(cvalue, n, cname, rc) 
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          read(cname,*) noahmp%nmlist%soil_level_nodes(n)
          write(msg, fmt='(A,I2.2,A)') trim(subname)//' : soil_level_nodes(',n,') = '//trim(cname)
          call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
       end do
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in soil_level_nodes", ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return 
    end if

    ! get dynamic_vegetation_option
    call NUOPC_CompAttributeGet(gcomp, name='dynamic_vegetation_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%dynamic_vegetation_option 
    else
       noahmp%nmlist%dynamic_vegetation_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : dynamic_vegetation_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%idveg = noahmp%nmlist%dynamic_vegetation_option

    ! get canopy_stomatal_resistance_option
    call NUOPC_CompAttributeGet(gcomp, name='canopy_stomatal_resistance_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%canopy_stomatal_resistance_option
    else
       noahmp%nmlist%canopy_stomatal_resistance_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : canopy_stomatal_resistance_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_crs = noahmp%nmlist%canopy_stomatal_resistance_option

    ! get soil_wetness_option
    call NUOPC_CompAttributeGet(gcomp, name='soil_wetness_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%soil_wetness_option
    else
       noahmp%nmlist%soil_wetness_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : soil_wetness_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_btr = noahmp%nmlist%soil_wetness_option

    ! get runoff_option
    call NUOPC_CompAttributeGet(gcomp, name='runoff_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%runoff_option
    else
       noahmp%nmlist%runoff_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : runoff_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_run = noahmp%nmlist%runoff_option

    ! get surface_exchange_option
    call NUOPC_CompAttributeGet(gcomp, name='surface_exchange_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%surface_exchange_option
    else
       noahmp%nmlist%surface_exchange_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : surface_exchange_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_sfc = noahmp%nmlist%surface_exchange_option

    ! get supercooled_soilwater_option
    call NUOPC_CompAttributeGet(gcomp, name='supercooled_soilwater_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%supercooled_soilwater_option
    else
       noahmp%nmlist%supercooled_soilwater_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : supercooled_soilwater_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_frz = noahmp%nmlist%supercooled_soilwater_option

    ! get frozen_soil_adjust_option
    call NUOPC_CompAttributeGet(gcomp, name='frozen_soil_adjust_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%frozen_soil_adjust_option
    else
       noahmp%nmlist%frozen_soil_adjust_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : frozen_soil_adjust_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_inf = noahmp%nmlist%frozen_soil_adjust_option

    ! get radiative_transfer_option
    call NUOPC_CompAttributeGet(gcomp, name='radiative_transfer_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%radiative_transfer_option
    else
       noahmp%nmlist%radiative_transfer_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : radiative_transfer_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_rad = noahmp%nmlist%radiative_transfer_option

    ! get snow_albedo_option
    call NUOPC_CompAttributeGet(gcomp, name='snow_albedo_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%snow_albedo_option
    else
       noahmp%nmlist%snow_albedo_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : snow_albedo_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_alb = noahmp%nmlist%snow_albedo_option

    ! get precip_partition_option
    call NUOPC_CompAttributeGet(gcomp, name='precip_partition_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%precip_partition_option
    else
       noahmp%nmlist%precip_partition_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : precip_partition_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_snf = noahmp%nmlist%precip_partition_option

    ! get soil_temp_lower_bdy_option
    call NUOPC_CompAttributeGet(gcomp, name='soil_temp_lower_bdy_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%soil_temp_lower_bdy_option
    else
       noahmp%nmlist%soil_temp_lower_bdy_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : soil_temp_lower_bdy_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_tbot = noahmp%nmlist%soil_temp_lower_bdy_option

    ! get soil_temp_time_scheme_option
    call NUOPC_CompAttributeGet(gcomp, name='soil_temp_time_scheme_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%soil_temp_time_scheme_option
    else
       noahmp%nmlist%soil_temp_time_scheme_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : soil_temp_time_scheme_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_stc = noahmp%nmlist%soil_temp_time_scheme_option

    ! get surface_evap_resistance_option
    call NUOPC_CompAttributeGet(gcomp, name='surface_evap_resistance_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%surface_evap_resistance_option
    else
       noahmp%nmlist%surface_evap_resistance_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : surface_evap_resistance_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_rsf = noahmp%nmlist%surface_evap_resistance_option

    ! get glacier_option
    call NUOPC_CompAttributeGet(gcomp, name='glacier_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%glacier_option
    else
       noahmp%nmlist%glacier_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : glacier_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_gla = noahmp%nmlist%glacier_option

    ! get surface_thermal_roughness_option
    call NUOPC_CompAttributeGet(gcomp, name='surface_thermal_roughness_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%surface_thermal_roughness_option
    else
       noahmp%nmlist%surface_thermal_roughness_option = -999
    end if
    call ESMF_LogWrite(trim(subname)//' : surface_thermal_roughness_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_trs = noahmp%nmlist%surface_thermal_roughness_option

    ! get surface_diagnose_approach_option
    call NUOPC_CompAttributeGet(gcomp, name='surface_diagnose_approach_option', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%surface_diagnose_approach_option
    else
       noahmp%nmlist%surface_diagnose_approach_option = 2
    end if
    call ESMF_LogWrite(trim(subname)//' : surface_diagnose_approach_option = '//trim(cvalue), ESMF_LOGMSG_INFO)
    noahmp%static%iopt_diag = noahmp%nmlist%surface_diagnose_approach_option

    ! output frequency
    call NUOPC_CompAttributeGet(gcomp, name='output_freq', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%output_freq
    else
       noahmp%nmlist%output_freq = 6
    end if
    call ESMF_LogWrite(trim(subname)//' : output_freq = '//trim(cvalue), ESMF_LOGMSG_INFO)

    ! output mode (high, low, debug)
    call NUOPC_CompAttributeGet(gcomp, name='output_mode', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%output_mode
       if (trim(noahmp%nmlist%output_mode) == 'all' .or. &
           trim(noahmp%nmlist%output_mode) == 'mid' .or. &
           trim(noahmp%nmlist%output_mode) == 'low') then
       else
         call ESMF_LogWrite(trim(subname)//": ERROR in output_mode. Only 'all', 'mid' and 'low' are allowed!", ESMF_LOGMSG_INFO)
         rc = ESMF_FAILURE
         return 
       end if
    else
       noahmp%nmlist%output_mode = 'all'
    end if

    ! restart frequency
    call NUOPC_CompAttributeGet(gcomp, name='restart_freq', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%restart_freq
    else
       noahmp%nmlist%restart_freq = noahmp%nmlist%output_freq
    end if
    write(msg, fmt='(A,I6)') trim(subname)//': restart_freq = ', noahmp%nmlist%restart_freq
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! MYNN-EDMF
    call NUOPC_CompAttributeGet(gcomp, name='do_mynnedmf', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%static%do_mynnedmf = .false.
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.' .or. trim(cvalue) .eq. 'true') noahmp%static%do_mynnedmf = .true.
    end if
    write(msg, fmt='(A,L)') trim(subname)//': do_mynnedmf = ', noahmp%static%do_mynnedmf
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! MYNN Surface Layer Scheme
    call NUOPC_CompAttributeGet(gcomp, name='do_mynnsfclay', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    noahmp%static%do_mynnsfclay = .false.
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.' .or. trim(cvalue) .eq. 'true') noahmp%model%do_mynnsfclay = .true.
    end if
    noahmp%model%do_mynnsfclay = noahmp%static%do_mynnsfclay
    if (noahmp%static%iopt_sfc == 4) noahmp%model%do_mynnsfclay = .true.
    write(msg, fmt='(A,L)') trim(subname)//': do_mynnsfclay = ', noahmp%static%do_mynnsfclay
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! Option for soil type category
    call NUOPC_CompAttributeGet(gcomp, name='soil_type_category', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%soil_type_category
    else
       noahmp%nmlist%soil_type_category = 1
    end if
    write(msg, fmt='(A,I1)') trim(subname)//' : soil_type_category (isot) = ', noahmp%nmlist%soil_type_category
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    noahmp%static%isot = noahmp%nmlist%soil_type_category

    ! Option for vegetation type category
    call NUOPC_CompAttributeGet(gcomp, name='veg_type_category', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%veg_type_category
    else
       noahmp%nmlist%veg_type_category = 1
    end if
    write(msg, fmt='(A,I1)') trim(subname)//' : veg_type_category (ivegsrc) = ', noahmp%nmlist%veg_type_category
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    noahmp%static%ivegsrc = noahmp%nmlist%veg_type_category

    ! Initial value of emissivity (constant in everywhere)
    call NUOPC_CompAttributeGet(gcomp, name='initial_emiss', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%initial_emiss
    else
       noahmp%nmlist%initial_emiss = 0.95
    end if
    write(msg, fmt='(A,F6.2)') trim(subname)//' : initial_emiss = ', noahmp%nmlist%initial_emiss
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! Initial value of monthly albedo (constant in everywhere)
    call NUOPC_CompAttributeGet(gcomp, name='initial_albedo', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) noahmp%nmlist%initial_albedo
    else
       noahmp%nmlist%initial_albedo = 0.2
    end if
    write(msg, fmt='(A,F6.2)') trim(subname)//' : initial_albedo = ', noahmp%nmlist%initial_albedo
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine read_namelist

  !===============================================================================
  subroutine shr_string_listGetName(list, k, name, rc)

    ! ----------------------------------------------
    ! Get name of k-th field in list
    ! It is adapted from CDEPS, shr_string_listGetName
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*)     , intent(in)  :: list    ! list/string
    integer          , intent(in)  :: k       ! index of field
    character(*)     , intent(out) :: name    ! k-th name in list
    integer          , intent(out) :: rc

    ! local variables
    integer :: i,n     ! generic indecies
    integer :: kFlds   ! number of fields in list
    integer :: i0,i1   ! name = list(i0:i1)
    character(*), parameter :: subName = '(shr_string_listGetName)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    !--- check that this is a valid index ---
    kFlds = shr_string_listGetNum(list)
    if (k < 1 .or. kFlds < k) then
      call ESMF_LogWrite(trim(subname)//": ERROR invalid index ", ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
    end if

    !--- start with whole list, then remove fields before and after desired field ---
    i0 = 1
    i1 = len_trim(list)

    !--- remove field names before desired field ---
    do n=2,k
       i = index(list(i0:i1),listDel)
       i0 = i0 + i
    end do

    !--- remove field names after desired field ---
    if ( k < kFlds ) then
       i = index(list(i0:i1),listDel)
       i1 = i0 + i - 2
    end if

    !--- copy result into output variable ---
    name = list(i0:i1)//"   "

  end subroutine shr_string_listGetName

  !===============================================================================
  integer function shr_string_listGetNum(str)

    ! ----------------------------------------------
    ! Get number of fields in a string list
    ! It is adapted from CDEPS, shr_string_listGetNum
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*), intent(in) :: str   ! string to search

    ! local variables
    integer :: count ! counts occurances of char
    character(*), parameter :: subName = '(shr_string_listGetNum'
    ! ----------------------------------------------

    shr_string_listGetNum = 0

    if (len_trim(str) > 0) then
       count = shr_string_countChar(str,listDel)
       shr_string_listGetNum = count + 1
    endif

  end function shr_string_listGetNum

  !===============================================================================
  integer function shr_string_countChar(str,char,rc)

    ! ----------------------------------------------
    ! Count number of occurances of a character
    ! It is adapted from CDEPS, shr_string_countChar
    ! ----------------------------------------------

    implicit none

    ! input/output variables
    character(*), intent(in)       :: str   ! string to search
    character(1), intent(in)       :: char  ! char to search for
    integer, intent(out), optional :: rc    ! return code

    ! local variables
    integer :: count    ! counts occurances of char
    integer :: n        ! generic index
    character(*), parameter :: subName = '(shr_string_countChar)'
    ! ----------------------------------------------

    count = 0
    do n = 1, len_trim(str)
      if (str(n:n) == char) count = count + 1
    end do
    shr_string_countChar = count

  end function shr_string_countChar

end module lnd_comp_shr
