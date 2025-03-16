module NoahmpReadLandMod

  use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE, C_PTR, C_CHAR
  use netcdf
  use Machine
  use NoahmpIOVarType

  implicit none

  public :: NoahmpReadLandHeader, NoahmpReadLandMain
  private :: FATAL, NOT_FATAL, get_2d_netcdf, error_handler, &
             get_2d_netcdf_cfloat, get_2d_netcdf_ffloat, get_2d_netcdf_finteger

  logical, parameter :: FATAL = .TRUE.
  logical, parameter :: NOT_FATAL = .FALSE.

  interface get_2d_netcdf
    module procedure get_2d_netcdf_cfloat
    module procedure get_2d_netcdf_ffloat
    module procedure get_2d_netcdf_finteger
  end interface get_2d_netcdf

contains

subroutine NoahmpReadLandHeader(NoahmpIO)

    implicit none
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

    integer :: ncid, dimid, varid, ierr
    real, allocatable, dimension(:,:) :: dum2d
    character(len=256) :: units
    integer :: i
    integer :: rank

    write(*,'("erfinput_flnm: ''", A, "''")') trim(NoahmpIO%erf_setup_file)

    ierr = nf90_open(NoahmpIO%erf_setup_file, NF90_NOWRITE, ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file))

    !!ierr = nf90_inq_dimid(ncid, "NX", dimid)
    !!call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding dimension 'NX'")

    !!ierr = nf90_inquire_dimension(ncid, dimid, len=NoahmpIO%xend)
    !!call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding dimension length for 'NX'")

    !!ierr = nf90_inq_dimid(ncid, "NY", dimid)
    !!call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding dimension 'NY'")

    !!ierr = nf90_inquire_dimension(ncid, dimid, len=NoahmpIO%yend)
    !!call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding dimension length for 'NY'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "DX", NoahmpIO%dx)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'DX'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "DY", NoahmpIO%dy)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'DY'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISWATER", NoahmpIO%iswater)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISWATER'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISLAKE", NoahmpIO%islake)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISLAKE'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISURBAN", NoahmpIO%isurban)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISURBAN'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISICE", NoahmpIO%isice)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISICE'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "LLANDUSE", NoahmpIO%llanduse)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'LLANDUSE'")   
 
    ! IBM XLF seems to need something like this:
    do i = 1, 256
       if (ichar(NoahmpIO%llanduse(i:i)) == 0) NoahmpIO%llanduse(i:i) = " "
    enddo

    ierr = nf90_close(ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems closing NetCDF file.")

end subroutine NoahmpReadLandHeader

subroutine NoahmpReadLandMain(NoahmpIO)
    implicit none
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

    character(len=256) :: units
    integer :: ierr
    integer :: ncid
    integer :: i

    real, dimension(NoahmpIO%xstart:NoahmpIO%xend, NoahmpIO%ystart:NoahmpIO%yend) :: soildummy

    write(*,'("erfinput_flnm: ''", A, "''")') trim(NoahmpIO%erf_setup_file)

    ierr = nf90_open(NoahmpIO%erf_setup_file, NF90_NOWRITE, ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file))

    call get_2d_netcdf("TERRAIN", ncid, NoahmpIO%TERRAIN, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)
    call get_2d_netcdf("ISLTYP", ncid, NoahmpIO%ISLTYP, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)

    call get_2d_netcdf("LAI", ncid, NoahmpIO%LAI, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)
    call get_2d_netcdf("VEGFRA", ncid, NoahmpIO%VEGFRA, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)

    call get_2d_netcdf("SMOIS", ncid, soildummy, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)
    do i=1,NoahmpIO%nsoil
       NoahmpIO%SMOIS(:,i,:) = soildummy
    end do

    call get_2d_netcdf("TSLB", ncid, soildummy, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)
    do i=1,NoahmpIO%nsoil
       NoahmpIO%TSLB(:,i,:) = soildummy
    end do

    call get_2d_netcdf("SNOW", ncid, NoahmpIO%SNOW, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)
    call get_2d_netcdf("SNOWH", ncid, NoahmpIO%SNOWH, units, NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, FATAL, ierr)
    NoahmpIO%FNDSNOWH = .TRUE.
  
    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= 0) stop "MODULE_NOAHLSM_ERF_INPUT:  READLAND_ERF:  NF90_CLOSE"

end subroutine NoahmpReadLandMain

subroutine get_2d_netcdf_cfloat(name, ncid, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real(c_double), dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=*), intent(out) :: units
    integer :: iret, varid
    ! FATAL_IF_ERROR:  an input code value:
    !      .TRUE. if an error in reading the data should stop the program.
    !      Otherwise the, IERR error flag is set, but the program continues.
    logical, intent(in) :: fatal_if_error 
    integer, intent(out) :: ierr
 
    units = " "
    
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid = ', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file.")
       else
          ierr = iret
         return
       endif
    endif

    iret = nf90_get_att(ncid, varid, "units", units)
    if (iret /= 0) units = "units unknown"

    iret = nf90_get_var(ncid, varid, values=array, start=(/xstart+1,ystart+1/), count=(/xend-xstart+1,yend-ystart+1/))

    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid =', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem retrieving variable '"//trim(name)//"' from NetCDF file.")
       else
          ierr = iret
          return
       endif
    endif
    ierr = 0;
end subroutine get_2d_netcdf_cfloat

subroutine get_2d_netcdf_ffloat(name, ncid, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real, dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=*), intent(out) :: units
    integer :: iret, varid
    ! FATAL_IF_ERROR:  an input code value:
    !      .TRUE. if an error in reading the data should stop the program.
    !      Otherwise the, IERR error flag is set, but the program continues.
    logical, intent(in) :: fatal_if_error 
    integer, intent(out) :: ierr
 
    units = " "
    
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid = ', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file.")
       else
          ierr = iret
         return
       endif
    endif

    iret = nf90_get_att(ncid, varid, "units", units)
    if (iret /= 0) units = "units unknown"

    iret = nf90_get_var(ncid, varid, values=array, start=(/xstart+1,ystart+1/), count=(/xend-xstart+1,yend-ystart+1/))

    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid =', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem retrieving variable '"//trim(name)//"' from NetCDF file.")
       else
          ierr = iret
          return
       endif
    endif
    ierr = 0;
end subroutine get_2d_netcdf_ffloat

subroutine get_2d_netcdf_finteger(name, ncid, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    integer, dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=*), intent(out) :: units
    integer :: iret, varid
    ! FATAL_IF_ERROR:  an input code value:
    !      .TRUE. if an error in reading the data should stop the program.
    !      Otherwise the, IERR error flag is set, but the program continues.
    logical, intent(in) :: fatal_if_error 
    integer, intent(out) :: ierr
 
    units = " "
    
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid = ', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file.")
       else
          ierr = iret
         return
       endif
    endif

    iret = nf90_get_att(ncid, varid, "units", units)
    if (iret /= 0) units = "units unknown"

    iret = nf90_get_var(ncid, varid, values=array, start=(/xstart+1,ystart+1/), count=(/xend-xstart+1,yend-ystart+1/))

    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid =', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem retrieving variable '"//trim(name)//"' from NetCDF file.")
       else
          ierr = iret
          return
       endif
    endif
    ierr = 0;
end subroutine get_2d_netcdf_finteger

subroutine error_handler(status, failure, success)
    !
    ! Check the error flag from a NetCDF function call, and print appropriate
    ! error message.
    !
    implicit none
    integer,                    intent(in) :: status
    character(len=*), optional, intent(in) :: failure
    character(len=*), optional, intent(in) :: success

    if (status .ne. NF90_NOERR) then
       write(*,'(/,A)') nf90_strerror(status)
       if (present(failure)) then
          write(*,'(/," ***** ", A,/)') failure
       endif
       stop 'Stopped'
    endif

    if (present(success)) then
       write(*,'(A)') success
    endif
end subroutine error_handler

end module NoahmpReadLandMod
