module NoahmpReadRestartMod

! Read the NoahMP prognostic state back into the NoahmpIO arrays on restart,
! overwriting cold-init state with the checkpointed values for a bitwise restart.
! Optional fields are read only if present in the file and allocated (a missing
! varid is skipped). Layer counts are asserted against the file.

   use mpi
   use netcdf
   use Machine, only : kind_noahmp, c_kind_noahmp
   use NoahmpIOVarType
   use NoahmpFatalMod, only : check_nc, NoahmpIO_abort

   implicit none

   integer, save, private :: ncid

contains

   subroutine NoahmpReadRestart(NoahmpIO, dir, maxblocks)

      implicit none

      type(NoahmpIO_type), intent(inout) :: NoahmpIO
      character(len=*),    intent(in)    :: dir
      integer,             intent(in)    :: maxblocks

      integer :: ierr, start(2), count(2)
      integer :: file_nsoil, file_nsnow
      character(len=1)   :: lev_str
      character(len=512) :: filename

      if (NoahmpIO%blkid == 0) then
         write (lev_str, '(I1.1)') NoahmpIO%LEVEL
         filename = trim(dir)//"/Level_"//trim(lev_str)//".nc"

         ierr = nf90_open(trim(filename), IOR(NF90_NOWRITE, NF90_MPIIO), &
                          ncid, comm=NoahmpIO%comm, info=MPI_INFO_NULL)
         if (ierr /= nf90_noerr) then
            print *, "NoahmpReadRestart: NetCDF open failed: ", trim(nf90_strerror(ierr))
            call NoahmpIO_abort()
         end if

         ! Assert the layer geometry matches; mismatch => corrupt restore.
         call check_nc(nf90_get_att(ncid, NF90_GLOBAL, "NSOIL", file_nsoil), "get_att NSOIL")
         call check_nc(nf90_get_att(ncid, NF90_GLOBAL, "NSNOW", file_nsnow), "get_att NSNOW")
         if (file_nsoil /= NoahmpIO%NSOIL .or. file_nsnow /= NoahmpIO%NSNOW) then
            print *, "NoahmpReadRestart: layer mismatch. file NSOIL/NSNOW=", &
                     file_nsoil, file_nsnow, " run=", NoahmpIO%NSOIL, NoahmpIO%NSNOW
            call NoahmpIO_abort()
         end if
      end if

      start = (/NoahmpIO%xstart-NoahmpIO%xoffset+1, NoahmpIO%ystart-NoahmpIO%yoffset+1/)
      count = (/NoahmpIO%xend-NoahmpIO%xstart+1,    NoahmpIO%yend-NoahmpIO%ystart+1/)

      ! --- soil
      call get3d(ncid, "TSLB",  NoahmpIO%TSLB,  start, count, NoahmpIO%NSOIL, .true.)
      call get3d(ncid, "SMOIS", NoahmpIO%SMOIS, start, count, NoahmpIO%NSOIL, .true.)
      call get3d(ncid, "SH2O",  NoahmpIO%SH2O,  start, count, NoahmpIO%NSOIL, .true.)
      if (allocated(NoahmpIO%SMOISEQ)) &
         call get3d(ncid, "SMOISEQ", NoahmpIO%SMOISEQ, start, count, NoahmpIO%NSOIL, .false.)

      ! --- soil albedo, banded (carried state; .false. so pre-fix checkpoints still restart -> cold-init 0.0) (#3446)
      call get3d(ncid, "ALBSOILDIRXY", NoahmpIO%ALBSOILDIRXY, start, count, NoahmpIO%NUMRAD, .false.)
      call get3d(ncid, "ALBSOILDIFXY", NoahmpIO%ALBSOILDIFXY, start, count, NoahmpIO%NUMRAD, .false.)

      ! --- snow layers
      call get3d(ncid, "TSNOXY",  NoahmpIO%TSNOXY,  start, count, NoahmpIO%NSNOW, .true.)
      call get3d(ncid, "SNICEXY", NoahmpIO%SNICEXY, start, count, NoahmpIO%NSNOW, .true.)
      call get3d(ncid, "SNLIQXY", NoahmpIO%SNLIQXY, start, count, NoahmpIO%NSNOW, .true.)
      call get3d(ncid, "ZSNSOXY", NoahmpIO%ZSNSOXY, start, count, NoahmpIO%NSNOW+NoahmpIO%NSOIL, .true.)

      ! --- snowpack scalars
      call get2d (ncid, "SNOW",   NoahmpIO%SNOW,    start, count, .true.)
      call get2d (ncid, "SNOWH",  NoahmpIO%SNOWH,   start, count, .true.)
      call get2d (ncid, "SNOWC",  NoahmpIO%SNOWC,   start, count, .true.)
      call get2di(ncid, "ISNOWXY",NoahmpIO%ISNOWXY, start, count, .true.)
      call get2d (ncid, "CANWAT", NoahmpIO%CANWAT,  start, count, .true.)
      call get2d (ncid, "ACSNOM", NoahmpIO%ACSNOM,  start, count, .true.)
      call get2d (ncid, "ACSNOW", NoahmpIO%ACSNOW,  start, count, .true.)

      ! --- canopy / surface
      call get2d(ncid, "TVXY",    NoahmpIO%TVXY,     start, count, .true.)
      call get2d(ncid, "TGXY",    NoahmpIO%TGXY,     start, count, .true.)
      call get2d(ncid, "CANICEXY",NoahmpIO%CANICEXY, start, count, .true.)
      call get2d(ncid, "CANLIQXY",NoahmpIO%CANLIQXY, start, count, .true.)
      call get2d(ncid, "EAHXY",   NoahmpIO%EAHXY,    start, count, .true.)
      call get2d(ncid, "TAHXY",   NoahmpIO%TAHXY,    start, count, .true.)
      call get2d(ncid, "CMXY",    NoahmpIO%CMXY,     start, count, .true.)
      call get2d(ncid, "CHXY",    NoahmpIO%CHXY,     start, count, .true.)
      call get2d(ncid, "FWETXY",  NoahmpIO%FWETXY,   start, count, .true.)
      call get2d(ncid, "QSFC",    NoahmpIO%QSFC,     start, count, .true.)
      call get2dd(ncid, "TSK",    NoahmpIO%TSK,      start, count, .true.)  ! c_kind_noahmp
      call get2d(ncid, "QSNOWXY", NoahmpIO%QSNOWXY,  start, count, .true.)
      call get2d(ncid, "QRAINXY", NoahmpIO%QRAINXY,  start, count, .true.)

      ! --- albedo history
      call get2d(ncid, "SNEQVOXY",NoahmpIO%SNEQVOXY, start, count, .true.)
      call get2d(ncid, "ALBOLDXY",NoahmpIO%ALBOLDXY, start, count, .true.)
      call get2d(ncid, "TAUSSXY", NoahmpIO%TAUSSXY,  start, count, .true.)
      call get2d(ncid, "ALBEDO",  NoahmpIO%ALBEDO,   start, count, .true.)

      ! --- aquifer / groundwater
      call get2d(ncid, "ZWTXY",     NoahmpIO%ZWTXY,      start, count, .true.)
      call get2d(ncid, "WAXY",      NoahmpIO%WAXY,       start, count, .true.)
      call get2d(ncid, "WTXY",      NoahmpIO%WTXY,       start, count, .true.)
      call get2d(ncid, "SMCWTDXY",  NoahmpIO%SMCWTDXY,   start, count, .true.)
      call get2d(ncid, "DEEPRECHXY",NoahmpIO%DEEPRECHXY, start, count, .true.)
      call get2d(ncid, "RECHXY",    NoahmpIO%RECHXY,     start, count, .true.)

      ! --- phenology
      call get2d(ncid, "LAI",    NoahmpIO%LAI,    start, count, .true.)
      call get2d(ncid, "XSAIXY", NoahmpIO%XSAIXY, start, count, .true.)

      ! --- accumulators / carried state
      call get2d(ncid, "SFCRUNOFF",NoahmpIO%SFCRUNOFF, start, count, .true.)
      call get2d(ncid, "UDRUNOFF", NoahmpIO%UDRUNOFF,  start, count, .true.)
      call get2d(ncid, "SMSTAV",   NoahmpIO%SMSTAV,    start, count, .true.)
      call get2d(ncid, "SMSTOT",   NoahmpIO%SMSTOT,    start, count, .true.)
      call get2dd(ncid,"EMISS",    NoahmpIO%EMISS,     start, count, .true.)  ! c_kind_noahmp
      call get2d(ncid, "GRDFLX",   NoahmpIO%GRDFLX,    start, count, .true.)

      ! --- soil-cycle accumulators (mid-cycle carry when SOIL_UPDATE_STEPS>1); optional
      ! so a pre-fix checkpoint still restarts, leaving them at cold-init 0.
      call get2d(ncid, "ACC_SSOILXY",  NoahmpIO%ACC_SSOILXY,  start, count, .false.)
      call get2d(ncid, "ACC_QINSURXY", NoahmpIO%ACC_QINSURXY, start, count, .false.)
      call get2d(ncid, "ACC_QSEVAXY",  NoahmpIO%ACC_QSEVAXY,  start, count, .false.)
      call get2d(ncid, "ACC_DWATERXY", NoahmpIO%ACC_DWATERXY, start, count, .false.)
      call get2d(ncid, "ACC_PRCPXY",   NoahmpIO%ACC_PRCPXY,   start, count, .false.)
      call get2d(ncid, "ACC_ECANXY",   NoahmpIO%ACC_ECANXY,   start, count, .false.)
      call get2d(ncid, "ACC_ETRANXY",  NoahmpIO%ACC_ETRANXY,  start, count, .false.)
      call get2d(ncid, "ACC_EDIRXY",   NoahmpIO%ACC_EDIRXY,   start, count, .false.)
      call get3d(ncid, "ACC_ETRANIXY", NoahmpIO%ACC_ETRANIXY, start, count, NoahmpIO%NSOIL, .false.)
      call get2d(ncid, "ACC_GLAFLWXY", NoahmpIO%ACC_GLAFLWXY, start, count, .false.)

      ! --- optional carbon / lake (only if allocated; missing var is skipped)
      if (allocated(NoahmpIO%LFMASSXY)) call get2d(ncid, "LFMASSXY", NoahmpIO%LFMASSXY, start, count, .false.)
      if (allocated(NoahmpIO%RTMASSXY)) call get2d(ncid, "RTMASSXY", NoahmpIO%RTMASSXY, start, count, .false.)
      if (allocated(NoahmpIO%STMASSXY)) call get2d(ncid, "STMASSXY", NoahmpIO%STMASSXY, start, count, .false.)
      if (allocated(NoahmpIO%WOODXY))   call get2d(ncid, "WOODXY",   NoahmpIO%WOODXY,   start, count, .false.)
      if (allocated(NoahmpIO%GRAINXY))  call get2d(ncid, "GRAINXY",  NoahmpIO%GRAINXY,  start, count, .false.)
      if (allocated(NoahmpIO%GDDXY))    call get2d(ncid, "GDDXY",    NoahmpIO%GDDXY,    start, count, .false.)
      if (allocated(NoahmpIO%WSLAKEXY)) call get2dd(ncid, "WSLAKEXY", NoahmpIO%WSLAKEXY, start, count, .false.)  ! c_kind_noahmp

      ! --- optional-mode carried state (restored only if this run enabled the mode
      ! AND the checkpoint has it; pre-fix / mode-off checkpoints leave cold-init values).
      if (allocated(NoahmpIO%FASTCPXY)) call get2d (ncid, "FASTCPXY", NoahmpIO%FASTCPXY, start, count, .false.)
      if (allocated(NoahmpIO%STBLCPXY)) call get2d (ncid, "STBLCPXY", NoahmpIO%STBLCPXY, start, count, .false.)
      if (allocated(NoahmpIO%PGSXY))    call get2di(ncid, "PGSXY",    NoahmpIO%PGSXY,    start, count, .false.)
      if (allocated(NoahmpIO%IRNUMSI))  call get2di(ncid, "IRNUMSI",  NoahmpIO%IRNUMSI,  start, count, .false.)
      if (allocated(NoahmpIO%IRNUMMI))  call get2di(ncid, "IRNUMMI",  NoahmpIO%IRNUMMI,  start, count, .false.)
      if (allocated(NoahmpIO%IRNUMFI))  call get2di(ncid, "IRNUMFI",  NoahmpIO%IRNUMFI,  start, count, .false.)
      if (allocated(NoahmpIO%IRWATSI))  call get2d (ncid, "IRWATSI",  NoahmpIO%IRWATSI,  start, count, .false.)
      if (allocated(NoahmpIO%IRWATMI))  call get2d (ncid, "IRWATMI",  NoahmpIO%IRWATMI,  start, count, .false.)
      if (allocated(NoahmpIO%IRWATFI))  call get2d (ncid, "IRWATFI",  NoahmpIO%IRWATFI,  start, count, .false.)
      if (allocated(NoahmpIO%WSURFXY))  call get2d (ncid, "WSURFXY",  NoahmpIO%WSURFXY,  start, count, .false.)
      if (allocated(NoahmpIO%FSATXY))   call get2d (ncid, "FSATXY",   NoahmpIO%FSATXY,   start, count, .false.)
      if (allocated(NoahmpIO%SNRDSXY))  call get3d(ncid, "SNRDSXY", NoahmpIO%SNRDSXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%BCPHIXY))  call get3d(ncid, "BCPHIXY", NoahmpIO%BCPHIXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%BCPHOXY))  call get3d(ncid, "BCPHOXY", NoahmpIO%BCPHOXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%OCPHIXY))  call get3d(ncid, "OCPHIXY", NoahmpIO%OCPHIXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%OCPHOXY))  call get3d(ncid, "OCPHOXY", NoahmpIO%OCPHOXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%DUST1XY))  call get3d(ncid, "DUST1XY", NoahmpIO%DUST1XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%DUST2XY))  call get3d(ncid, "DUST2XY", NoahmpIO%DUST2XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%DUST3XY))  call get3d(ncid, "DUST3XY", NoahmpIO%DUST3XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%DUST4XY))  call get3d(ncid, "DUST4XY", NoahmpIO%DUST4XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%DUST5XY))  call get3d(ncid, "DUST5XY", NoahmpIO%DUST5XY, start, count, NoahmpIO%NSNOW, .false.)
      ! SNICAR freeze rate + mass concentrations: recomputed every step from the
      ! fields above, but only after this restart's first Advance -- restore them
      ! so that first step isn't wrong too.
      if (allocated(NoahmpIO%SNFRXY))          call get3d(ncid, "SNFRXY",          NoahmpIO%SNFRXY,          start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcBCPHIXY)) call get3d(ncid, "MassConcBCPHIXY", NoahmpIO%MassConcBCPHIXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcBCPHOXY)) call get3d(ncid, "MassConcBCPHOXY", NoahmpIO%MassConcBCPHOXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcOCPHIXY)) call get3d(ncid, "MassConcOCPHIXY", NoahmpIO%MassConcOCPHIXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcOCPHOXY)) call get3d(ncid, "MassConcOCPHOXY", NoahmpIO%MassConcOCPHOXY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcDUST1XY)) call get3d(ncid, "MassConcDUST1XY", NoahmpIO%MassConcDUST1XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcDUST2XY)) call get3d(ncid, "MassConcDUST2XY", NoahmpIO%MassConcDUST2XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcDUST3XY)) call get3d(ncid, "MassConcDUST3XY", NoahmpIO%MassConcDUST3XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcDUST4XY)) call get3d(ncid, "MassConcDUST4XY", NoahmpIO%MassConcDUST4XY, start, count, NoahmpIO%NSNOW, .false.)
      if (allocated(NoahmpIO%MassConcDUST5XY)) call get3d(ncid, "MassConcDUST5XY", NoahmpIO%MassConcDUST5XY, start, count, NoahmpIO%NSNOW, .false.)

      if (NoahmpIO%blkid == (maxblocks-1)) then
         call check_nc(nf90_close(ncid), "close")
      end if

   end subroutine NoahmpReadRestart

   ! helpers

   subroutine get2d(nc, name, arr, start, count, required)
      integer,          intent(in)    :: nc
      character(len=*), intent(in)    :: name
      real(kind=kind_noahmp), intent(inout) :: arr(:,:)
      integer,          intent(in)    :: start(2), count(2)
      logical,          intent(in)    :: required
      integer :: vid, ierr
      ierr = nf90_inq_varid(nc, name, vid)
      if (ierr /= nf90_noerr) then
         if (required) then
            print *, "NoahmpReadRestart: missing required var ", trim(name)
            call NoahmpIO_abort()
         end if
         return
      end if
      ierr = nf90_get_var(nc, vid, arr, start=start, count=count)
      if (ierr /= nf90_noerr) then
         print *, "NoahmpReadRestart: read failed for ", trim(name), ": ", trim(nf90_strerror(ierr))
         call NoahmpIO_abort()
      end if
   end subroutine get2d

   ! Variant for C-boundary 2D fields (TSK, EMISS, WSLAKEXY) declared with the
   ! C-interop kind c_kind_noahmp (== kind_noahmp in every build).
   subroutine get2dd(nc, name, arr, start, count, required)
      integer,          intent(in)    :: nc
      character(len=*), intent(in)    :: name
      real(kind=c_kind_noahmp), intent(inout) :: arr(:,:)
      integer,          intent(in)    :: start(2), count(2)
      logical,          intent(in)    :: required
      integer :: vid, ierr
      ierr = nf90_inq_varid(nc, name, vid)
      if (ierr /= nf90_noerr) then
         if (required) then
            print *, "NoahmpReadRestart: missing required var ", trim(name)
            call NoahmpIO_abort()
         end if
         return
      end if
      ierr = nf90_get_var(nc, vid, arr, start=start, count=count)
      if (ierr /= nf90_noerr) then
         print *, "NoahmpReadRestart: read failed for ", trim(name), ": ", trim(nf90_strerror(ierr))
         call NoahmpIO_abort()
      end if
   end subroutine get2dd

   subroutine get2di(nc, name, arr, start, count, required)
      integer,          intent(in)    :: nc
      character(len=*), intent(in)    :: name
      integer,          intent(inout) :: arr(:,:)
      integer,          intent(in)    :: start(2), count(2)
      logical,          intent(in)    :: required
      integer :: vid, ierr
      ierr = nf90_inq_varid(nc, name, vid)
      if (ierr /= nf90_noerr) then
         if (required) then
            print *, "NoahmpReadRestart: missing required var ", trim(name)
            call NoahmpIO_abort()
         end if
         return
      end if
      ierr = nf90_get_var(nc, vid, arr, start=start, count=count)
      if (ierr /= nf90_noerr) then
         print *, "NoahmpReadRestart: read failed for ", trim(name), ": ", trim(nf90_strerror(ierr))
         call NoahmpIO_abort()
      end if
   end subroutine get2di

   subroutine get3d(nc, name, arr, start, count, nk, required)
      integer,          intent(in)    :: nc
      character(len=*), intent(in)    :: name
      real(kind=kind_noahmp), intent(inout) :: arr(:,:,:)
      integer,          intent(in)    :: start(2), count(2), nk
      logical,          intent(in)    :: required
      integer :: vid, ierr
      ierr = nf90_inq_varid(nc, name, vid)
      if (ierr /= nf90_noerr) then
         if (required) then
            print *, "NoahmpReadRestart: missing required var ", trim(name)
            call NoahmpIO_abort()
         end if
         return
      end if
      ierr = nf90_get_var(nc, vid, arr, start=(/start(1),1,start(2)/), &
                                        count=(/count(1),nk,count(2)/))
      if (ierr /= nf90_noerr) then
         print *, "NoahmpReadRestart: read failed for ", trim(name), ": ", trim(nf90_strerror(ierr))
         call NoahmpIO_abort()
      end if
   end subroutine get3d

end module NoahmpReadRestartMod
