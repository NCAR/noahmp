!> @file
!> @brief Manage cpl_scalars
!> @author mvertens@ucar.edu
!> @author modified for NOAHMP by Denise.Worthen@noaa.gov @date 03-03-2024
!>
!> Manage scalars in import and export states. Called at realization to set the
!> required scalar data into a state. The scalar_value will be set into a field
!> with name flds_scalar_name. The scalar_id identifies which dimension in the
!> scalar field is given by the scalar_value. The number of scalars is used to
!> ensure that the scalar_id is within the bounds of the scalar field

module lnd_comp_cplscalars

  use NUOPC
  use ESMF, only : ESMF_Field, ESMF_Distgrid, ESMF_Grid, ESMF_State
  use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet, ESMF_VMBroadCast
  use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_FAILURE, ESMF_SUCCESS
  use ESMF, only : ESMF_LOGMSG_INFO, ESMF_LOGWRITE, ESMF_TYPEKIND_R8, ESMF_KIND_R8
  use ESMF, only : ESMF_GridCreate, ESMF_FieldCreate, ESMF_StateGet, ESMF_DistGridCreate
  use ESMF, only : ESMF_FieldGet

  implicit none

  private
  public SetScalarField
  public State_SetScalar
  public State_GetScalar

  integer, public           :: flds_scalar_num, flds_scalar_index_nx
  integer, public           :: flds_scalar_index_ny, flds_scalar_index_ntile
  character(len=80), public :: flds_scalar_name

contains

  !================================================================================
  !> Create a scalar field
  !>
  !> @param[inout]   field            an ESMF_Field
  !> @param[in]      flds_scalar_name the name of the scalar
  !> @param[in]      flds_scalar_num  the number of scalars
  !> @param[inout]   rc               a return code
  subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)

    type(ESMF_Field) , intent(inout) :: field
    character(len=*) , intent(in)    :: flds_scalar_name
    integer          , intent(in)    :: flds_scalar_num
    integer          , intent(inout) :: rc

    ! local variables
    type(ESMF_Distgrid) :: distgrid
    type(ESMF_Grid)     :: grid

    character(len=*), parameter :: subname='(lnd_comp_cplscalars:SetScalarField)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
    distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    grid = ESMF_GridCreate(distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
         ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc) ! num of scalar values
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetScalarField

  !================================================================================
  !> Set scalar data into a state
  !>
  !> @param[inout]   State            an ESMF_State
  !> @param[in]      scalar_value     the value of the scalar
  !> @param[in]      scalar_id        the identity of the scalar
  !> @param[in]      flds_scalar_name the name of the scalar
  !> @param[in]      flds_scalar_num  the number of scalars
  !> @param[inout]   rc               a return code
  subroutine State_SetScalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

    ! input/output arguments
    real(ESMF_KIND_R8), intent(in)   :: scalar_value
    integer,          intent(in)     :: scalar_id
    type(ESMF_State), intent(inout)  :: State
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: lfield
    type(ESMF_VM)     :: vm
    real(ESMF_KIND_R8), pointer :: farrayptr(:,:)

    character(len=*), parameter :: subname = ' (lnd_comp_cplscalars:state_setscalar) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=lfield, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (mytask == 0) then
       call ESMF_FieldGet(lfield, farrayPtr = farrayptr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       endif
       farrayptr(scalar_id,1) = scalar_value
    endif

  end subroutine State_SetScalar

  !===============================================================================
  !> Get scalar data from a state
  !>
  !> @details Obtain the field flds_scalar_name from a State and broadcast and
  !> it to all PEs
  !>
  !> @param[in]    State            an ESMF_State
  !> @param[in]    scalar_value     the value of the scalar
  !> @param[in]    scalar_id        the identity of the scalar
  !> @param[in]    flds_scalar_name the name of the scalar
  !> @param[in]    flds_scalar_num  the number of scalars
  !> @param[out]   rc               a return code
  subroutine State_GetScalar(state, scalar_id, scalar_value, flds_scalar_name, flds_scalar_num, rc)

    ! ----------------------------------------------
    ! Get scalar data from State for a particular name and broadcast it to all other pets
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)            :: state
    integer,          intent(in)            :: scalar_id
    real(ESMF_KIND_R8),         intent(out) :: scalar_value
    character(len=*), intent(in)            :: flds_scalar_name
    integer,          intent(in)            :: flds_scalar_num
    integer,          intent(inout)         :: rc

    ! local variables
    integer                     :: mytask, ierr, icount
    type(ESMF_VM)               :: vm
    type(ESMF_Field)            :: field
    real(ESMF_KIND_R8), pointer :: farrayptr(:,:)
    real(ESMF_KIND_R8)          :: tmp(1)

    character(len=*), parameter :: subname = ' (lnd_comp_cplscalars:state_getscalar) '
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! check item exist or not?
    call ESMF_StateGet(State, itemSearch=trim(flds_scalar_name), itemCount=icount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (icount > 0) then
      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__)
          rc = ESMF_FAILURE
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif
        tmp(:) = farrayptr(scalar_id,:)
      endif
      call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      scalar_value = tmp(1)
    else
      scalar_value = 0.0_ESMF_KIND_R8
      call ESMF_LogWrite(trim(subname)//": no ESMF_Field found named: "//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
    end if

  end subroutine State_GetScalar
end module lnd_comp_cplscalars
