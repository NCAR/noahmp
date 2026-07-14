module NoahmpReadLandMod

  use, intrinsic :: iso_c_binding, only: C_INT, C_PTR, C_CHAR
  use netcdf
  use Machine
  use NoahmpIOVarType
  use NoahmpFatalMod, only : NoahmpIO_abort

  implicit none

  public :: NoahmpReadLandHeader, NoahmpReadLandMain

  private :: FATAL, NOT_FATAL, get_2d_netcdf, get_2d_netcdf_c, error_handler, &
             get_landuse_netcdf, get_soilcat_netcdf, get_netcdf_soillevel, init_interp

  logical, parameter :: FATAL = .TRUE.
  logical, parameter :: NOT_FATAL = .FALSE.

  ! get_2d_netcdf reads 2-D fields into kind_noahmp arrays; get_2d_netcdf_c reads
  ! the ERF C++-owned fields (XLAT, TSK) with the C-interop kind c_kind_noahmp
  ! (== kind_noahmp in every build, so the readers are otherwise identical).

contains

subroutine NoahmpReadLandHeader(NoahmpIO)

    implicit none
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

    integer :: ncid, dimid, varid, ierr
    real(kind_noahmp), allocatable, dimension(:,:) :: dum2d  ! scratch to extract lat1/lon1 scalars
    character(len=256) :: units
    integer :: i
    integer :: rank
    integer :: ilev, is, js, ratio, xoffset, yoffset

    if (NoahmpIO%rank == 0) write(*,'("Noah-MP reading ''", A, "'' headers")') trim(NoahmpIO%erf_setup_file_lev)

    ierr = nf90_open(NoahmpIO%erf_setup_file_lev, NF90_NOWRITE, ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_lev))

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION", NoahmpIO%xsglobal)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'WEST-EAST_GRID_DIMENSION'")
    NoahmpIO%xsglobal = NoahmpIO%xsglobal-1

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", NoahmpIO%ysglobal)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'SOUTH-NORTH_GRID_DIMENSION'")
    NoahmpIO%ysglobal = NoahmpIO%ysglobal-1

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "DX", NoahmpIO%dx)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'DX'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "DY", NoahmpIO%dy)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'DY'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1", NoahmpIO%truelat1)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'TRUELAT1'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2", NoahmpIO%truelat2)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'TRUELAT2'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", NoahmpIO%cen_lon)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'STAND_LON'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "MAP_PROJ", NoahmpIO%mapproj)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'MAP_PROJ'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "GRID_ID", NoahmpIO%igrid)
    if (ierr /= 0) then
       ierr = nf90_get_att(ncid, NF90_GLOBAL, "grid_id", NoahmpIO%igrid)
       call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'GRID_ID' or 'grid_id'")
    endif

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISWATER", NoahmpIO%iswater)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISWATER'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISLAKE", NoahmpIO%islake)
    if(ierr /= 0) then
      if (NoahmpIO%rank == 0) write(*,*) "Problems finding global attribute: ISLAKE; setting to -1"
      NoahmpIO%islake = -1
    end if

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISURBAN", NoahmpIO%isurban)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISURBAN'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISICE", NoahmpIO%isice)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISICE'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "MMINLU", NoahmpIO%llanduse)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'MMINLU'")
    
    ! IBM XLF seems to need something like this:
    do i = 1, 256
       if (ichar(NoahmpIO%llanduse(i:i)) == 0) NoahmpIO%llanduse(i:i) = " "
    enddo

    allocate(dum2d(NoahmpIO%xstart:NoahmpIO%xstart,NoahmpIO%ystart:NoahmpIO%ystart))
    call get_2d_netcdf("XLAT", ncid, dum2d,  units, NoahmpIO%xstart, NoahmpIO%xstart, NoahmpIO%ystart, NoahmpIO%ystart, FATAL, ierr)
    NoahmpIO%lat1 = dum2d(NoahmpIO%xstart,NoahmpIO%ystart)

    call get_2d_netcdf("XLONG", ncid, dum2d,  units, NoahmpIO%xstart, NoahmpIO%xstart, NoahmpIO%ystart, NoahmpIO%ystart, FATAL, ierr)
    NoahmpIO%lon1 = dum2d(NoahmpIO%xstart,NoahmpIO%ystart)
    deallocate (dum2d)

    ierr = nf90_close(ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems closing NetCDF file.")

    NoahmpIO%xoffset = 0
    NoahmpIO%yoffset = 0

    do ilev=0,NoahmpIO%level

      select case(ilev)
      case(0)
        ierr = nf90_open(NoahmpIO%erf_setup_file_01, NF90_NOWRITE, ncid)
        call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_01))
      case(1)
        ierr = nf90_open(NoahmpIO%erf_setup_file_02, NF90_NOWRITE, ncid)
        call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_02))
      case(2)
        ierr = nf90_open(NoahmpIO%erf_setup_file_03, NF90_NOWRITE, ncid)
        call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_03)) 
      case default
        if (NoahmpIO%rank == 0) print *, "Error: unsupported level: ", ilev
        call NoahmpIO_abort()
      end select

      ierr = nf90_get_att(ncid, NF90_GLOBAL, "I_PARENT_START", is)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'I_PARENT_START'")

      ierr = nf90_get_att(ncid, NF90_GLOBAL, "J_PARENT_START", js)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'J_PARENT_START'")

      ierr = nf90_get_att(ncid, NF90_GLOBAL, "PARENT_GRID_RATIO", ratio)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'PARENT_GRID_RATIO'")

      ierr = nf90_close(ncid)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems closing NetCDF file.")        

      NoahmpIO%xoffset = ratio*(NoahmpIO%xoffset+is-1)
      NoahmpIO%yoffset = ratio*(NoahmpIO%yoffset+js-1)

    end do

end subroutine NoahmpReadLandHeader

subroutine NoahmpReadLandMain(NoahmpIO)
    implicit none
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

    character(len=256) :: units
    integer :: ierr
    integer :: ncid
    real, dimension(NoahmpIO%xstart-NoahmpIO%xoffset:NoahmpIO%xend-NoahmpIO%xoffset, &
                    NoahmpIO%ystart-NoahmpIO%yoffset:NoahmpIO%yend-NoahmpIO%yoffset) :: xdum
    character(len=256) :: llanduse

    integer :: ierr_snodep, varid
    integer :: isoil
    real(kind_noahmp), dimension(100) :: layer_bottom
    real(kind_noahmp), dimension(100) :: layer_top
    real(kind_noahmp), dimension(NoahmpIO%nsoil)   :: dzs

    real(kind_noahmp), dimension(NoahmpIO%xstart-NoahmpIO%xoffset:NoahmpIO%xend-NoahmpIO%xoffset, &
                    NoahmpIO%nsoil, &
                    NoahmpIO%ystart-NoahmpIO%yoffset:NoahmpIO%yend-NoahmpIO%yoffset) :: soildummy

    integer :: ierr_vegfra
    integer :: ierr_lai

    integer :: i, j
    integer :: xstart, ystart, xend, yend

    if (NoahmpIO%rank == 0) write(*,'("Noah-MP reading ''", A, "'' variables")') trim(NoahmpIO%erf_setup_file_lev)

    ierr = nf90_open(NoahmpIO%erf_setup_file_lev, NF90_NOWRITE, ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_lev))

    xstart = NoahmpIO%xstart-NoahmpIO%xoffset
    ystart = NoahmpIO%ystart-NoahmpIO%yoffset

    xend = NoahmpIO%xend-NoahmpIO%xoffset
    yend = NoahmpIO%yend-NoahmpIO%yoffset

    call get_2d_netcdf_c("XLAT", ncid, NoahmpIO%xlat,  units, xstart, xend, ystart, yend, FATAL, ierr)  ! c_kind_noahmp

    call get_2d_netcdf("XLONG", ncid, NoahmpIO%xlong, units, xstart, xend, ystart, yend, FATAL, ierr)

    call get_2d_netcdf("XLAND", ncid, NoahmpIO%xland, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)

    call get_2d_netcdf("SEAICE", ncid, NoahmpIO%seaice, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)

    call get_2d_netcdf("HGT", ncid, NoahmpIO%terrain, units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Deep-layer soil temperature
    call get_2d_netcdf("TMN", ncid, NoahmpIO%TMN, units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Map factors, only needed for iopt_run=5
    call get_2d_netcdf("MAPFAC_MX", ncid, NoahmpIO%msftx, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)
    if (ierr /= 0) print*, 'Did not find MAPFAC_MX, only needed for iopt_run=5'

    call get_2d_netcdf("MAPFAC_MY", ncid, NoahmpIO%msfty, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)
    if (ierr /= 0) print*, 'Did not find MAPFAC_MY, only needed for iopt_run=5'

    ! Dominant land-use category
    call get_landuse_netcdf(ncid, xdum , units, xstart, xend, ystart, yend)
    NoahmpIO%ivgtyp = nint(xdum)

    ! Dominant top-layer soil-type category
    call get_soilcat_netcdf(ncid, xdum , units, xstart, xend, ystart, yend)
    NoahmpIO%isltyp = nint(xdum)

    where (NoahmpIO%SEAICE > 0.0) NoahmpIO%XICE = 1.0
 
    NoahmpIO%CROPTYPE   = 0       ! no crops by default

    NoahmpIO%TD_FRACTION = 0.0

    NoahmpIO%SLOPETYP  =  1 ! matches noahmpdrv
    NoahmpIO%DZS       =  NoahmpIO%SOIL_THICK_INPUT(1:NoahmpIO%NSOIL)
    NoahmpIO%ITIMESTEP = 1
    NoahmpIO%restart_flag = .false.

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "MMINLU", llanduse)
    if (ierr /= 0) then
       if (NoahmpIO%rank == 0) write(*,'("WARNING:  Noah-MP input file does not have MMINLU attribute.")')
       if (NoahmpIO%rank == 0) write(*,'("          This probably means that the file  is from an older release.")')
       if (NoahmpIO%rank == 0) write(*,'("          I assume you know what you are doing.")')
    else
       if (NoahmpIO%rank == 0) write(*,'("MMNINLU attribute: ", A)') llanduse
    endif

    call get_2d_netcdf("CANWAT", ncid, NoahmpIO%canwat, units, xstart, xend, ystart, yend, FATAL, ierr)
    call get_2d_netcdf_c("TSK",  ncid, NoahmpIO%tsk, units, xstart, xend, ystart, yend, FATAL, ierr)  ! c_kind_noahmp
    call get_2d_netcdf("SNOW",   ncid, NoahmpIO%snow, units, xstart, xend, ystart, yend, FATAL, ierr)
    call get_2d_netcdf("SNOWC",  ncid, NoahmpIO%snowc, units, xstart, xend, ystart, yend, FATAL, ierr)

    NoahmpIO%snowh = 0.0
    call get_2d_netcdf("SNOWH", ncid, NoahmpIO%snowh, units, xstart, xend, ystart, yend, NOT_FATAL, ierr_snodep)
    NoahmpIO%fndsnowh = .true.
    if (ierr_snodep /= 0) NoahmpIO%fndsnowh = .false.
   
    ierr = nf90_inq_varid(ncid,  "DZS",  varid)
    call error_handler(ierr, "READLAND_ERF:  Problem finding variable 'DZS' in the wrfinput file.")
    ierr = nf90_get_var(ncid, varid, values=dzs, start=(/1/), count=(/NoahmpIO%nsoil/))
    call error_handler(ierr, "READLAND_ERF:  Problem retrieving variable 'DZS' from the wrfinput file.")

    layer_top(1) = 0.0
    layer_bottom(1) = dzs(1)
    do isoil = 2, NoahmpIO%nsoil
      layer_top(isoil) = layer_bottom(isoil-1)
      layer_bottom(isoil) = layer_top(isoil) + dzs(isoil)
    end do

    call get_netcdf_soillevel("TSLB", ncid, NoahmpIO%nsoil, soildummy, units,  xstart, xend, ystart, yend, FATAL, ierr)

    call init_interp(NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, NoahmpIO%nsoil, &
                     NoahmpIO%dzs, NoahmpIO%tslb, NoahmpIO%nsoil, soildummy, layer_bottom(1:NoahmpIO%nsoil), layer_top(1:NoahmpIO%nsoil), NoahmpIO%rank)

    call get_netcdf_soillevel("SMOIS", ncid, NoahmpIO%nsoil, soildummy, units,  xstart, xend, ystart, yend, FATAL, ierr)

    call init_interp(NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, NoahmpIO%nsoil, &
    NoahmpIO%dzs, NoahmpIO%smois, NoahmpIO%nsoil, soildummy, layer_bottom(1:NoahmpIO%nsoil), layer_top(1:NoahmpIO%nsoil), NoahmpIO%rank)

    NoahmpIO%VEGFRA =  0.0

    call get_2d_netcdf("VEGFRA", ncid, NoahmpIO%vegfra, units, xstart, xend, ystart, yend, NOT_FATAL, ierr_vegfra)
    call get_2d_netcdf("LAI", ncid, NoahmpIO%lai, units, xstart, xend, ystart, yend, NOT_FATAL, ierr_lai)

    ! Min/max green vegetation fraction
    call get_2d_netcdf("SHDMIN", ncid, NoahmpIO%gvfmin, units, xstart, xend, ystart, yend, FATAL, ierr)

    call get_2d_netcdf("SHDMAX", ncid, NoahmpIO%gvfmax, units, xstart, xend, ystart, yend, FATAL, ierr)

    ierr = nf90_close(ncid)
    call error_handler(ierr, "MODULE_NOAHLSM_ERF_INPUT:  READLAND_ERF:  NF90_CLOSE")

end subroutine NoahmpReadLandMain

subroutine get_2d_netcdf(name, ncid, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real(kind_noahmp), dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=*), intent(out) :: units
    integer :: iret, varid
    logical, intent(in) :: fatal_if_error  ! .TRUE. aborts on error; else set ierr and return
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

end subroutine get_2d_netcdf

! Variant for the ERF C++-owned fields (XLAT, TSK) with C-interop kind c_kind_noahmp;
! identical to get_2d_netcdf, kept distinct to mark the C boundary (cf. get2dd).
subroutine get_2d_netcdf_c(name, ncid, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real(c_kind_noahmp), dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=*), intent(out) :: units
    integer :: iret, varid
    logical, intent(in) :: fatal_if_error  ! .TRUE. aborts on error; else set ierr and return
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

end subroutine get_2d_netcdf_c

subroutine error_handler(status, failure, success)
    ! Abort with a message if a NetCDF status flag indicates failure.
    implicit none
    integer,                    intent(in) :: status
    character(len=*), optional, intent(in) :: failure
    character(len=*), optional, intent(in) :: success

    if (status .ne. NF90_NOERR) then
       write(*,'(/,A)') nf90_strerror(status)
       if (present(failure)) then
          write(*,'(/," ***** ", A,/)') failure
       endif
       call NoahmpIO_abort()
    endif

    if (present(success)) then
       write(*,'(A)') success
    endif
end subroutine error_handler

subroutine get_landuse_netcdf(ncid, array, units, xstart, xend, ystart, yend)
    implicit none
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real, dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=256), intent(out) :: units
    integer :: iret, varid
    character(len=24), parameter :: name = "IVGTYP"

    units = " "

    iret = nf90_inq_varid(ncid,  trim(name),  varid)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       call error_handler(iret, "MODULE_NOAHLSM_ERF_INPUT:  get_landuse_netcdf:  nf90_inq_varid")
    endif

    iret = nf90_get_var(ncid, varid, array, (/xstart+1, ystart+1/), (/xend-xstart+1, yend-ystart+1/))
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       call error_handler(iret, "MODULE_NOAHLSM_ERF_INPUT:  get_landuse_netcdf:  nf90_get_var")
    endif
end subroutine get_landuse_netcdf

subroutine get_soilcat_netcdf(ncid, array, units, xstart, xend, ystart, yend)
    implicit none
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real, dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=256), intent(out) :: units
    integer :: iret, varid
    character(len=24), parameter :: name = "ISLTYP"

    units = " "

    iret = nf90_inq_varid(ncid,  trim(name),  varid)
    call error_handler(iret, "Problem finding variable '"//trim(name)//"' in the wrfinput file.")

    iret = nf90_get_var(ncid, varid, array, (/xstart+1, ystart+1/), (/xend-xstart+1, yend-ystart+1/))
    call error_handler(iret, "Problem retrieving variable "//trim(name)//" from the wrfinput file.")

end subroutine get_soilcat_netcdf

subroutine get_netcdf_soillevel(name, ncid, nsoil, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: nsoil
    integer, intent(in) :: xstart, xend, ystart, yend
    real(kind_noahmp), dimension(xstart:xend,nsoil,ystart:yend), intent(out) :: array
    character(len=256), intent(out) :: units
    logical, intent(in) :: fatal_if_error
    integer, intent(out) :: ierr

    integer :: iret, varid, isoil
    real(kind_noahmp):: insoil(xstart:xend,ystart:yend,nsoil)

    units = " "

    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'name = "', trim(name)//'"'
          call error_handler(iret, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_2d_netcdf:  nf90_inq_varid")
       else
          ierr = iret
          return
       endif
    endif

    iret = nf90_get_att(ncid, varid, "units", units)
    if (iret /= 0) units = "units unknown"

    iret = nf90_get_var(ncid, varid, values=insoil, start=(/xstart+1,ystart+1,1,1/), count=(/xend-xstart+1,yend-ystart+1,nsoil,1/))
    ! Check the read before reshaping so a non-fatal failure does not copy garbage.
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'name = "', trim(name)//'"'
          print*, 'varid =', varid
          call error_handler(iret, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_2d_netcdf:  nf90_get_var")
       else
          ierr = iret
          return
       endif
    endif

    do isoil = 1,nsoil
      array(:,isoil,:) = insoil(:,:,isoil)
    end do

    ierr = 0;
end subroutine get_netcdf_soillevel

subroutine init_interp(xstart, xend, ystart, yend, nsoil, sldpth, var, nvar, src, layer_bottom, layer_top, rank)
    implicit none
    integer, intent(in)    :: xstart, xend, ystart, yend, nsoil, nvar
    real(kind_noahmp), dimension(nsoil) :: sldpth ! the thickness of each layer
    real(kind_noahmp), dimension(xstart:xend, nsoil, ystart:yend), intent(out) :: var
    real(kind_noahmp), dimension(xstart:xend, nvar, ystart:yend ), intent(in)  :: src
    real(kind_noahmp), dimension(nvar),               intent(in)  :: layer_bottom ! The depth from the surface of each layer bottom.
    real(kind_noahmp), dimension(nvar),               intent(in)  :: layer_top    ! The depth from the surface of each layer top.
    integer :: i, j, k, kk, ktop, kbottom
    real(kind_noahmp), dimension(nsoil) :: dst_centerpoint
    real(kind_noahmp), dimension(nvar)  :: src_centerpoint
    real(kind_noahmp) :: fraction
    integer :: ierr
    integer, intent(in) :: rank

    do k = 1, nsoil
       if (k==1) then
          dst_centerpoint(k) = sldpth(k)/2.
       else
          dst_centerpoint(k) = sldpth(k)/2. + sum(sldpth(1:k-1))
       endif
    enddo

    do k = 1, nvar
       src_centerpoint(k) = 0.5*(layer_bottom(k)+layer_top(k))
    enddo

    KLOOP : do k = 1, nsoil

       if (dst_centerpoint(k) < src_centerpoint(1)) then
          ! Destination center shallower than the topmost source: use topmost source
          var(:,k,:) = src(:,1,:)
          cycle KLOOP
       endif

       if (dst_centerpoint(k) > src_centerpoint(nvar)) then
          ! Destination center deeper than the deepest source: use deepest source
          var(:,k,:) = src(:,nvar,:)
          cycle KLOOP
       endif

       ! If the destination center is "close" to a source center, use that layer
       do kk = 1, nvar
          if (abs(dst_centerpoint(k)-src_centerpoint(kk)) < 0.01) then
             var(:,k,:) = src(:,kk,:)
             cycle KLOOP
          endif
       enddo

       ! Otherwise, do a linear interpolation

       ! ktop: top bracketing source layer (first, from the bottom up, shallower
       ! than the destination level)
       ktop = -99999
       TOPLOOP : do kk = nvar,1,-1
          if (src_centerpoint(kk) < dst_centerpoint(k)) then
             ktop = kk
             exit TOPLOOP
          endif
       enddo TOPLOOP
       if (ktop < -99998) then
          if (rank == 0) write(*,'("***** ERROR: ktop problem in soil layer interpolation")')
          call NoahmpIO_abort()
       endif

       ! kbottom: bottom bracketing source layer (first, from the top down, deeper
       ! than the destination level)
       kbottom = -99999
       BOTTOMLOOP : do kk = 1, nvar
          if ( src_centerpoint(kk) > dst_centerpoint(k) ) then
             kbottom = kk
             exit BOTTOMLOOP
          endif
       enddo BOTTOMLOOP
       if (kbottom < -99998) then
          if (rank == 0) write(*,'("***** ERROR: kbottom problem in soil layer interpolation")')
          call NoahmpIO_abort()
       endif

       fraction = (src_centerpoint(kbottom)-dst_centerpoint(k)) / (src_centerpoint(kbottom)-src_centerpoint(ktop))

       var(:,k,:) = (src(:,ktop,:)*fraction) + (src(:,kbottom,:)*(1.0-fraction))

    enddo KLOOP     
end subroutine init_interp


end module NoahmpReadLandMod
