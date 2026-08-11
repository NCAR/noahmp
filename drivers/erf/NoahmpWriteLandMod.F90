module NoahmpWriteLandMod

   use mpi
   use netcdf
   use NoahmpIOVarType
   use NoahmpFatalMod, only : check_nc, NoahmpIO_abort

   implicit none

   integer, save, private :: ncid, terrain, snowh, shbxy, evbxy, &
                             vegfra, gvfmin, gvfmax, tsk, emiss, &
                             albsfcdirxy, albsfcdifxy, &
                             savxy, sagxy, pahxy, firaxy, hfx, &
                             lh, grdflx, ghbxy, canhsxy, tslb, smois, &
                             tau_ew, tau_ns, &
                             ! core surface diagnostics
                             t2mvxy, t2mbxy, q2mvxy, q2mbxy, tradxy, &
                             fvegxy, tgvxy, tgbxy, shgxy, shcxy, &
                             evgxy, evcxy, trxy, runsfxy, runsbxy, &
                             ecanxy, edirxy, etranxy, fsaxy, rs, z0, znt

contains

   subroutine NoahmpWriteLand(NoahmpIO, filenum, maxblocks)

      implicit none

      type(NoahmpIO_type), intent(inout) :: NoahmpIO
      integer, intent(in) :: filenum, maxblocks

      ! local variables
      integer :: ierr, start(2), count(2), nx, ny, comp2d, nsoil
      character(len=32) :: ts_str   ! wide enough to avoid overflow for large step counts
      character(len=1) :: lev_str
      character(len=100) :: dir, filename
      logical :: ex

      if (NoahmpIO%blkid == 0) then
         write (ts_str, '(I12.5)') filenum   ! zero-pad to >=5 digits, grows for larger counts
         ts_str = adjustl(ts_str)            ! drop leading blanks, keep leading zeros
         write (lev_str, '(I1.1)') NoahmpIO%LEVEL

         dir = "lnd"//trim(ts_str)
         inquire (file=trim(dir), exist=ex)
         if (.not. ex) then
            call execute_command_line("mkdir -p "//trim(dir), exitstat=ierr)
            if (ierr /= 0) then
               print *, "Failed to create directory: ", trim(dir)
               call NoahmpIO_abort()
            end if
         end if

         filename = trim(dir)//"/Level_"//trim(lev_str)//".nc"
         call check_nc(nf90_create(trim(filename), IOR(NF90_CLOBBER, IOR(NF90_NETCDF4, NF90_MPIIO)), &
                       ncid, comm=NoahmpIO%comm, info=MPI_INFO_NULL), "create "//trim(filename))

         ! Define dimensions and variable
         call check_nc(nf90_def_dim(ncid, "NX", NoahmpIO%xsglobal, nx), "def_dim NX")
         call check_nc(nf90_def_dim(ncid, "NY", NoahmpIO%ysglobal, ny), "def_dim NY")
         call check_nc(nf90_def_dim(ncid, "COMP2D", 2, comp2d), "def_dim COMP2D")
         call check_nc(nf90_def_dim(ncid, "NSOIL", NoahmpIO%NSOIL, nsoil), "def_dim NSOIL")
         call check_nc(nf90_def_var(ncid, "TERRAIN", NF90_FLOAT, (/nx, ny/), terrain), "def_var TERRAIN")
         call check_nc(nf90_def_var(ncid, "SNOWH", NF90_FLOAT, (/nx, ny/), snowh), "def_var SNOWH")
         call check_nc(nf90_def_var(ncid, "SHBXY", NF90_FLOAT, (/nx, ny/), shbxy), "def_var SHBXY")
         call check_nc(nf90_def_var(ncid, "EVBXY", NF90_FLOAT, (/nx, ny/), evbxy), "def_var EVBXY")
         call check_nc(nf90_def_var(ncid, "VEGFRA", NF90_FLOAT, (/nx, ny/), vegfra), "def_var VEGFRA")
         call check_nc(nf90_def_var(ncid, "GVFMIN", NF90_FLOAT, (/nx, ny/), gvfmin), "def_var GVFMIN")
         call check_nc(nf90_def_var(ncid, "GVFMAX", NF90_FLOAT, (/nx, ny/), gvfmax), "def_var GVFMAX")
         call check_nc(nf90_def_var(ncid, "TSK", NF90_FLOAT, (/nx, ny/), tsk), "def_var TSK")
         call check_nc(nf90_def_var(ncid, "EMISS", NF90_FLOAT, (/nx, ny/), emiss), "def_var EMISS")
         call check_nc(nf90_def_var(ncid, "ALBSFCDIRXY", NF90_FLOAT, (/nx, comp2d, ny/), albsfcdirxy), "def_var ALBSFCDIRXY")
         call check_nc(nf90_def_var(ncid, "ALBSFCDIFXY", NF90_FLOAT, (/nx, comp2d, ny/), albsfcdifxy), "def_var ALBSFCDIFXY")
         call check_nc(nf90_def_var(ncid, "SAVXY", NF90_FLOAT, (/nx, ny/), savxy), "def_var SAVXY")
         call check_nc(nf90_def_var(ncid, "SAGXY", NF90_FLOAT, (/nx, ny/), sagxy), "def_var SAGXY")
         call check_nc(nf90_def_var(ncid, "PAHXY", NF90_FLOAT, (/nx, ny/), pahxy), "def_var PAHXY")
         call check_nc(nf90_def_var(ncid, "FIRAXY", NF90_FLOAT, (/nx, ny/), firaxy), "def_var FIRAXY")
         call check_nc(nf90_def_var(ncid, "HFX", NF90_FLOAT, (/nx, ny/), hfx), "def_var HFX")
         call check_nc(nf90_def_var(ncid, "LH", NF90_FLOAT, (/nx, ny/), lh), "def_var LH")
         call check_nc(nf90_def_var(ncid, "GRDFLX", NF90_FLOAT, (/nx, ny/), grdflx), "def_var GRDFLX")
         call check_nc(nf90_def_var(ncid, "GHBXY", NF90_FLOAT, (/nx, ny/), ghbxy), "def_var GHBXY")
         call check_nc(nf90_def_var(ncid, "CANHSXY", NF90_FLOAT, (/nx, ny/), canhsxy), "def_var CANHSXY")
         call check_nc(nf90_def_var(ncid, "TSLB", NF90_FLOAT, (/nx, nsoil, ny/), tslb), "def_var TSLB")
         call check_nc(nf90_def_var(ncid, "SMOIS", NF90_FLOAT, (/nx, nsoil, ny/), smois), "def_var SMOIS")
         call check_nc(nf90_def_var(ncid, "TAU_EW", NF90_FLOAT, (/nx, ny/), tau_ew), "def_var TAU_EW")
         call check_nc(nf90_def_var(ncid, "TAU_NS", NF90_FLOAT, (/nx, ny/), tau_ns), "def_var TAU_NS")
         ! core surface diagnostics
         call check_nc(nf90_def_var(ncid, "T2MVXY", NF90_FLOAT, (/nx, ny/), t2mvxy), "def_var T2MVXY")
         call check_nc(nf90_def_var(ncid, "T2MBXY", NF90_FLOAT, (/nx, ny/), t2mbxy), "def_var T2MBXY")
         call check_nc(nf90_def_var(ncid, "Q2MVXY", NF90_FLOAT, (/nx, ny/), q2mvxy), "def_var Q2MVXY")
         call check_nc(nf90_def_var(ncid, "Q2MBXY", NF90_FLOAT, (/nx, ny/), q2mbxy), "def_var Q2MBXY")
         call check_nc(nf90_def_var(ncid, "TRADXY", NF90_FLOAT, (/nx, ny/), tradxy), "def_var TRADXY")
         call check_nc(nf90_def_var(ncid, "FVEGXY", NF90_FLOAT, (/nx, ny/), fvegxy), "def_var FVEGXY")
         call check_nc(nf90_def_var(ncid, "TGVXY", NF90_FLOAT, (/nx, ny/), tgvxy), "def_var TGVXY")
         call check_nc(nf90_def_var(ncid, "TGBXY", NF90_FLOAT, (/nx, ny/), tgbxy), "def_var TGBXY")
         call check_nc(nf90_def_var(ncid, "SHGXY", NF90_FLOAT, (/nx, ny/), shgxy), "def_var SHGXY")
         call check_nc(nf90_def_var(ncid, "SHCXY", NF90_FLOAT, (/nx, ny/), shcxy), "def_var SHCXY")
         call check_nc(nf90_def_var(ncid, "EVGXY", NF90_FLOAT, (/nx, ny/), evgxy), "def_var EVGXY")
         call check_nc(nf90_def_var(ncid, "EVCXY", NF90_FLOAT, (/nx, ny/), evcxy), "def_var EVCXY")
         call check_nc(nf90_def_var(ncid, "TRXY", NF90_FLOAT, (/nx, ny/), trxy), "def_var TRXY")
         call check_nc(nf90_def_var(ncid, "RUNSFXY", NF90_FLOAT, (/nx, ny/), runsfxy), "def_var RUNSFXY")
         call check_nc(nf90_def_var(ncid, "RUNSBXY", NF90_FLOAT, (/nx, ny/), runsbxy), "def_var RUNSBXY")
         call check_nc(nf90_def_var(ncid, "ECANXY", NF90_FLOAT, (/nx, ny/), ecanxy), "def_var ECANXY")
         call check_nc(nf90_def_var(ncid, "EDIRXY", NF90_FLOAT, (/nx, ny/), edirxy), "def_var EDIRXY")
         call check_nc(nf90_def_var(ncid, "ETRANXY", NF90_FLOAT, (/nx, ny/), etranxy), "def_var ETRANXY")
         call check_nc(nf90_def_var(ncid, "FSAXY", NF90_FLOAT, (/nx, ny/), fsaxy), "def_var FSAXY")
         call check_nc(nf90_def_var(ncid, "RS", NF90_FLOAT, (/nx, ny/), rs), "def_var RS")
         call check_nc(nf90_def_var(ncid, "Z0", NF90_FLOAT, (/nx, ny/), z0), "def_var Z0")
         call check_nc(nf90_def_var(ncid, "ZNT", NF90_FLOAT, (/nx, ny/), znt), "def_var ZNT")
         call check_nc(nf90_enddef(ncid), "enddef")
      end if

      ! Select portion to write
      start = (/NoahmpIO%xstart-NoahmpIO%xoffset+1, NoahmpIO%ystart-NoahmpIO%yoffset+1/)
      count = (/NoahmpIO%xend-NoahmpIO%xstart+1, NoahmpIO%yend-NoahmpIO%ystart+1/)

      ! Write data
      call check_nc(nf90_put_var(ncid, terrain, NoahmpIO%TERRAIN, start=start, count=count), "put_var TERRAIN")
      call check_nc(nf90_put_var(ncid, snowh, NoahmpIO%SNOWH, start=start, count=count), "put_var SNOWH")
      call check_nc(nf90_put_var(ncid, shbxy, NoahmpIO%SHBXY, start=start, count=count), "put_var SHBXY")
      call check_nc(nf90_put_var(ncid, evbxy, NoahmpIO%EVBXY, start=start, count=count), "put_var EVBXY")
      call check_nc(nf90_put_var(ncid, vegfra, NoahmpIO%VEGFRA, start=start, count=count), "put_var VEGFRA")
      call check_nc(nf90_put_var(ncid, gvfmin, NoahmpIO%GVFMIN, start=start, count=count), "put_var GVFMIN")
      call check_nc(nf90_put_var(ncid, gvfmax, NoahmpIO%GVFMAX, start=start, count=count), "put_var GVFMAX")
      call check_nc(nf90_put_var(ncid, tsk, NoahmpIO%TSK, start=start, count=count), "put_var TSK")
      call check_nc(nf90_put_var(ncid, emiss, NoahmpIO%EMISS, start=start, count=count), "put_var EMISS")
      call check_nc(nf90_put_var(ncid, albsfcdirxy, NoahmpIO%ALBSFCDIRXY, start=(/start(1), 1, start(2)/), &
                                                                   count=(/count(1), 2, count(2)/)), "put_var ALBSFCDIRXY")
      call check_nc(nf90_put_var(ncid, albsfcdifxy, NoahmpIO%ALBSFCDIFXY, start=(/start(1), 1, start(2)/), &
                                                                   count=(/count(1), 2, count(2)/)), "put_var ALBSFCDIFXY")
      call check_nc(nf90_put_var(ncid, savxy, NoahmpIO%SAVXY, start=start, count=count), "put_var SAVXY")
      call check_nc(nf90_put_var(ncid, sagxy, NoahmpIO%SAGXY, start=start, count=count), "put_var SAGXY")
      call check_nc(nf90_put_var(ncid, pahxy, NoahmpIO%PAHXY, start=start, count=count), "put_var PAHXY")
      call check_nc(nf90_put_var(ncid, firaxy, NoahmpIO%FIRAXY, start=start, count=count), "put_var FIRAXY")
      call check_nc(nf90_put_var(ncid, hfx, NoahmpIO%HFX, start=start, count=count), "put_var HFX")
      call check_nc(nf90_put_var(ncid, lh, NoahmpIO%LH, start=start, count=count), "put_var LH")
      call check_nc(nf90_put_var(ncid, grdflx, NoahmpIO%GRDFLX, start=start, count=count), "put_var GRDFLX")
      call check_nc(nf90_put_var(ncid, ghbxy, NoahmpIO%GHBXY, start=start, count=count), "put_var GHBXY")
      call check_nc(nf90_put_var(ncid, canhsxy, NoahmpIO%CANHSXY, start=start, count=count), "put_var CANHSXY")
      call check_nc(nf90_put_var(ncid, tslb, NoahmpIO%TSLB, start=(/start(1), 1, start(2)/), &
                                                     count=(/count(1), NoahmpIO%NSOIL, count(2)/)), "put_var TSLB")
      call check_nc(nf90_put_var(ncid, smois, NoahmpIO%SMOIS, start=(/start(1), 1, start(2)/), &
                                                       count=(/count(1), NoahmpIO%NSOIL, count(2)/)), "put_var SMOIS")
      call check_nc(nf90_put_var(ncid, tau_ew, NoahmpIO%TAU_EW, start=start, count=count), "put_var TAU_EW")
      call check_nc(nf90_put_var(ncid, tau_ns, NoahmpIO%TAU_NS, start=start, count=count), "put_var TAU_NS")
      ! core surface diagnostics
      call check_nc(nf90_put_var(ncid, t2mvxy, NoahmpIO%T2MVXY, start=start, count=count), "put_var T2MVXY")
      call check_nc(nf90_put_var(ncid, t2mbxy, NoahmpIO%T2MBXY, start=start, count=count), "put_var T2MBXY")
      call check_nc(nf90_put_var(ncid, q2mvxy, NoahmpIO%Q2MVXY, start=start, count=count), "put_var Q2MVXY")
      call check_nc(nf90_put_var(ncid, q2mbxy, NoahmpIO%Q2MBXY, start=start, count=count), "put_var Q2MBXY")
      call check_nc(nf90_put_var(ncid, tradxy, NoahmpIO%TRADXY, start=start, count=count), "put_var TRADXY")
      call check_nc(nf90_put_var(ncid, fvegxy, NoahmpIO%FVEGXY, start=start, count=count), "put_var FVEGXY")
      call check_nc(nf90_put_var(ncid, tgvxy, NoahmpIO%TGVXY, start=start, count=count), "put_var TGVXY")
      call check_nc(nf90_put_var(ncid, tgbxy, NoahmpIO%TGBXY, start=start, count=count), "put_var TGBXY")
      call check_nc(nf90_put_var(ncid, shgxy, NoahmpIO%SHGXY, start=start, count=count), "put_var SHGXY")
      call check_nc(nf90_put_var(ncid, shcxy, NoahmpIO%SHCXY, start=start, count=count), "put_var SHCXY")
      call check_nc(nf90_put_var(ncid, evgxy, NoahmpIO%EVGXY, start=start, count=count), "put_var EVGXY")
      call check_nc(nf90_put_var(ncid, evcxy, NoahmpIO%EVCXY, start=start, count=count), "put_var EVCXY")
      call check_nc(nf90_put_var(ncid, trxy, NoahmpIO%TRXY, start=start, count=count), "put_var TRXY")
      call check_nc(nf90_put_var(ncid, runsfxy, NoahmpIO%RUNSFXY, start=start, count=count), "put_var RUNSFXY")
      call check_nc(nf90_put_var(ncid, runsbxy, NoahmpIO%RUNSBXY, start=start, count=count), "put_var RUNSBXY")
      call check_nc(nf90_put_var(ncid, ecanxy, NoahmpIO%ECANXY, start=start, count=count), "put_var ECANXY")
      call check_nc(nf90_put_var(ncid, edirxy, NoahmpIO%EDIRXY, start=start, count=count), "put_var EDIRXY")
      call check_nc(nf90_put_var(ncid, etranxy, NoahmpIO%ETRANXY, start=start, count=count), "put_var ETRANXY")
      call check_nc(nf90_put_var(ncid, fsaxy, NoahmpIO%FSAXY, start=start, count=count), "put_var FSAXY")
      call check_nc(nf90_put_var(ncid, rs, NoahmpIO%RS, start=start, count=count), "put_var RS")
      call check_nc(nf90_put_var(ncid, z0, NoahmpIO%Z0, start=start, count=count), "put_var Z0")
      call check_nc(nf90_put_var(ncid, znt, NoahmpIO%ZNT, start=start, count=count), "put_var ZNT")

      if (NoahmpIO%blkid == (maxblocks-1)) then
         call check_nc(nf90_close(ncid), "close")
      end if

   end subroutine NoahmpWriteLand

end module NoahmpWriteLandMod
