module NoahmpWriteRestartMod

!------------------------------------------------------------------------------
! Write the full NoahMP prognostic state to a NetCDF restart file so that an
! ERF restart reproduces a cold-start trajectory bitwise. This mirrors the
! collective per-block MPI-IO pattern of NoahmpWriteLandMod, but:
!   * serializes the PROGNOSTIC state (soil, snow incl. layers, canopy/veg,
!     albedo history, aquifer, phenology, accumulators) rather than the small
!     diagnostic/coupling subset the lnd plotfile carries;
!   * stores reals at the model's working precision (kind_noahmp) so there is
!     no lossy conversion -- NF90_DOUBLE when kind_noahmp==8, else NF90_REAL;
!   * stores ISNOWXY (the active snow-layer count) as NF90_INT, which is
!     required to interpret the negative-indexed snow-layer arrays on restart;
!   * guards every optional (conditionally allocated) field with allocated().
! Companion reader: NoahmpReadRestartMod.
!------------------------------------------------------------------------------

   use mpi
   use netcdf
   use iso_c_binding, only : C_DOUBLE
   use Machine, only : kind_noahmp
   use NoahmpIOVarType

   implicit none

   integer, save, private :: ncid

   ! varids -- soil
   integer, save, private :: id_tslb, id_smois, id_sh2o, id_smoiseq
   ! varids -- snowpack scalars
   integer, save, private :: id_snow, id_snowh, id_snowc, id_isnow, &
                             id_canwat, id_acsnom, id_acsnow
   ! varids -- snow layers (negative-indexed)
   integer, save, private :: id_tsno, id_zsnso, id_snice, id_snliq
   ! varids -- canopy / surface
   integer, save, private :: id_tv, id_tg, id_canice, id_canliq, id_eah, &
                             id_tah, id_cm, id_ch, id_fwet, id_qsfc, id_tsk, &
                             id_qsnow, id_qrain
   ! varids -- albedo history
   integer, save, private :: id_sneqvo, id_albold, id_tauss, id_albedo
   ! varids -- aquifer / groundwater
   integer, save, private :: id_zwt, id_wa, id_wt, id_smcwtd, id_deeprech, id_rech
   ! varids -- phenology
   integer, save, private :: id_lai, id_xsai
   ! varids -- accumulators / misc carried state
   integer, save, private :: id_sfcrunoff, id_udrunoff, id_smstav, id_smstot, &
                             id_emiss, id_grdflx
   ! varids -- optional carbon / lake
   integer, save, private :: id_lfmass, id_rtmass, id_stmass, id_wood, &
                             id_grain, id_gdd, id_wslake

contains

   subroutine NoahmpWriteRestart(NoahmpIO, dir, maxblocks)

      implicit none

      type(NoahmpIO_type), intent(inout) :: NoahmpIO
      character(len=*),    intent(in)    :: dir
      integer,             intent(in)    :: maxblocks

      ! local variables
      integer :: ierr, start(2), count(2)
      integer :: nx, ny, nsoil_d, nsnow_d, nsnso_d
      integer :: rtype
      character(len=1)   :: lev_str
      character(len=512) :: filename
      logical :: ex

      ! NetCDF real type chosen to match the in-memory kind so the round-trip
      ! is bit-exact (no float<->double conversion).
      rtype = NF90_REAL
      if (kind_noahmp == 8) rtype = NF90_DOUBLE

      if (NoahmpIO%blkid == 0) then
         write (lev_str, '(I1.1)') NoahmpIO%LEVEL

         inquire (file=trim(dir), exist=ex)
         if (.not. ex) then
            call execute_command_line("mkdir -p "//trim(dir), exitstat=ierr)
            if (ierr /= 0) then
               print *, "NoahmpWriteRestart: failed to create directory: ", trim(dir)
               stop
            end if
         end if

         filename = trim(dir)//"/Level_"//trim(lev_str)//".nc"
         ierr = nf90_create(trim(filename), IOR(NF90_CLOBBER, IOR(NF90_NETCDF4, NF90_MPIIO)), &
                            ncid, comm=NoahmpIO%comm, info=MPI_INFO_NULL)
         if (ierr /= nf90_noerr) then
            print *, "NoahmpWriteRestart: NetCDF create failed: ", trim(nf90_strerror(ierr))
            stop
         end if

         ! Dimensions
         ierr = nf90_def_dim(ncid, "NX",    NoahmpIO%xsglobal, nx)
         ierr = nf90_def_dim(ncid, "NY",    NoahmpIO%ysglobal, ny)
         ierr = nf90_def_dim(ncid, "NSOIL", NoahmpIO%NSOIL,    nsoil_d)
         ierr = nf90_def_dim(ncid, "NSNOW", NoahmpIO%NSNOW,    nsnow_d)
         ierr = nf90_def_dim(ncid, "NSNSO", NoahmpIO%NSNOW+NoahmpIO%NSOIL, nsnso_d)

         ! Record the layer counts as global attributes for the read-side assert.
         ierr = nf90_put_att(ncid, NF90_GLOBAL, "NSOIL", NoahmpIO%NSOIL)
         ierr = nf90_put_att(ncid, NF90_GLOBAL, "NSNOW", NoahmpIO%NSNOW)

         ! --- soil (NX, NSOIL, NY)
         ierr = nf90_def_var(ncid, "TSLB",     rtype, (/nx, nsoil_d, ny/), id_tslb)
         ierr = nf90_def_var(ncid, "SMOIS",    rtype, (/nx, nsoil_d, ny/), id_smois)
         ierr = nf90_def_var(ncid, "SH2O",     rtype, (/nx, nsoil_d, ny/), id_sh2o)
         if (allocated(NoahmpIO%SMOISEQ)) &
            ierr = nf90_def_var(ncid, "SMOISEQ", rtype, (/nx, nsoil_d, ny/), id_smoiseq)

         ! --- snow layers (NX, NSNOW, NY) and snow+soil (NX, NSNSO, NY)
         ierr = nf90_def_var(ncid, "TSNOXY",  rtype, (/nx, nsnow_d, ny/), id_tsno)
         ierr = nf90_def_var(ncid, "SNICEXY", rtype, (/nx, nsnow_d, ny/), id_snice)
         ierr = nf90_def_var(ncid, "SNLIQXY", rtype, (/nx, nsnow_d, ny/), id_snliq)
         ierr = nf90_def_var(ncid, "ZSNSOXY", rtype, (/nx, nsnso_d, ny/), id_zsnso)

         ! --- snowpack scalars (NX, NY)
         ierr = nf90_def_var(ncid, "SNOW",    rtype,     (/nx, ny/), id_snow)
         ierr = nf90_def_var(ncid, "SNOWH",   rtype,     (/nx, ny/), id_snowh)
         ierr = nf90_def_var(ncid, "SNOWC",   rtype,     (/nx, ny/), id_snowc)
         ierr = nf90_def_var(ncid, "ISNOWXY", NF90_INT,  (/nx, ny/), id_isnow)
         ierr = nf90_def_var(ncid, "CANWAT",  rtype,     (/nx, ny/), id_canwat)
         ierr = nf90_def_var(ncid, "ACSNOM",  rtype,     (/nx, ny/), id_acsnom)
         ierr = nf90_def_var(ncid, "ACSNOW",  rtype,     (/nx, ny/), id_acsnow)

         ! --- canopy / surface
         ierr = nf90_def_var(ncid, "TVXY",    rtype, (/nx, ny/), id_tv)
         ierr = nf90_def_var(ncid, "TGXY",    rtype, (/nx, ny/), id_tg)
         ierr = nf90_def_var(ncid, "CANICEXY",rtype, (/nx, ny/), id_canice)
         ierr = nf90_def_var(ncid, "CANLIQXY",rtype, (/nx, ny/), id_canliq)
         ierr = nf90_def_var(ncid, "EAHXY",   rtype, (/nx, ny/), id_eah)
         ierr = nf90_def_var(ncid, "TAHXY",   rtype, (/nx, ny/), id_tah)
         ierr = nf90_def_var(ncid, "CMXY",    rtype, (/nx, ny/), id_cm)
         ierr = nf90_def_var(ncid, "CHXY",    rtype, (/nx, ny/), id_ch)
         ierr = nf90_def_var(ncid, "FWETXY",  rtype, (/nx, ny/), id_fwet)
         ierr = nf90_def_var(ncid, "QSFC",    rtype, (/nx, ny/), id_qsfc)
         ! TSK, EMISS are declared C_DOUBLE in NoahmpIO -> always NF90_DOUBLE.
         ierr = nf90_def_var(ncid, "TSK",     NF90_DOUBLE, (/nx, ny/), id_tsk)
         ierr = nf90_def_var(ncid, "QSNOWXY", rtype, (/nx, ny/), id_qsnow)
         ierr = nf90_def_var(ncid, "QRAINXY", rtype, (/nx, ny/), id_qrain)

         ! --- albedo history
         ierr = nf90_def_var(ncid, "SNEQVOXY",rtype, (/nx, ny/), id_sneqvo)
         ierr = nf90_def_var(ncid, "ALBOLDXY",rtype, (/nx, ny/), id_albold)
         ierr = nf90_def_var(ncid, "TAUSSXY", rtype, (/nx, ny/), id_tauss)
         ierr = nf90_def_var(ncid, "ALBEDO",  rtype, (/nx, ny/), id_albedo)

         ! --- aquifer / groundwater
         ierr = nf90_def_var(ncid, "ZWTXY",     rtype, (/nx, ny/), id_zwt)
         ierr = nf90_def_var(ncid, "WAXY",      rtype, (/nx, ny/), id_wa)
         ierr = nf90_def_var(ncid, "WTXY",      rtype, (/nx, ny/), id_wt)
         ierr = nf90_def_var(ncid, "SMCWTDXY",  rtype, (/nx, ny/), id_smcwtd)
         ierr = nf90_def_var(ncid, "DEEPRECHXY",rtype, (/nx, ny/), id_deeprech)
         ierr = nf90_def_var(ncid, "RECHXY",    rtype, (/nx, ny/), id_rech)

         ! --- phenology
         ierr = nf90_def_var(ncid, "LAI",     rtype, (/nx, ny/), id_lai)
         ierr = nf90_def_var(ncid, "XSAIXY",  rtype, (/nx, ny/), id_xsai)

         ! --- accumulators / carried state
         ierr = nf90_def_var(ncid, "SFCRUNOFF",rtype, (/nx, ny/), id_sfcrunoff)
         ierr = nf90_def_var(ncid, "UDRUNOFF", rtype, (/nx, ny/), id_udrunoff)
         ierr = nf90_def_var(ncid, "SMSTAV",   rtype, (/nx, ny/), id_smstav)
         ierr = nf90_def_var(ncid, "SMSTOT",   rtype, (/nx, ny/), id_smstot)
         ierr = nf90_def_var(ncid, "EMISS",    NF90_DOUBLE, (/nx, ny/), id_emiss)
         ierr = nf90_def_var(ncid, "GRDFLX",   rtype, (/nx, ny/), id_grdflx)

         ! --- optional carbon / dveg (only if allocated)
         if (allocated(NoahmpIO%LFMASSXY)) &
            ierr = nf90_def_var(ncid, "LFMASSXY", rtype, (/nx, ny/), id_lfmass)
         if (allocated(NoahmpIO%RTMASSXY)) &
            ierr = nf90_def_var(ncid, "RTMASSXY", rtype, (/nx, ny/), id_rtmass)
         if (allocated(NoahmpIO%STMASSXY)) &
            ierr = nf90_def_var(ncid, "STMASSXY", rtype, (/nx, ny/), id_stmass)
         if (allocated(NoahmpIO%WOODXY)) &
            ierr = nf90_def_var(ncid, "WOODXY",   rtype, (/nx, ny/), id_wood)
         if (allocated(NoahmpIO%GRAINXY)) &
            ierr = nf90_def_var(ncid, "GRAINXY",  rtype, (/nx, ny/), id_grain)
         if (allocated(NoahmpIO%GDDXY)) &
            ierr = nf90_def_var(ncid, "GDDXY",    rtype, (/nx, ny/), id_gdd)
         if (allocated(NoahmpIO%WSLAKEXY)) &
            ierr = nf90_def_var(ncid, "WSLAKEXY", NF90_DOUBLE, (/nx, ny/), id_wslake)

         ierr = nf90_enddef(ncid)
      end if

      ! Hyperslab for this block within the global domain.
      start = (/NoahmpIO%xstart-NoahmpIO%xoffset+1, NoahmpIO%ystart-NoahmpIO%yoffset+1/)
      count = (/NoahmpIO%xend-NoahmpIO%xstart+1,    NoahmpIO%yend-NoahmpIO%ystart+1/)

      ! --- soil
      ierr = nf90_put_var(ncid, id_tslb,  NoahmpIO%TSLB,  start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/))
      ierr = nf90_put_var(ncid, id_smois, NoahmpIO%SMOIS, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/))
      ierr = nf90_put_var(ncid, id_sh2o,  NoahmpIO%SH2O,  start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/))
      if (allocated(NoahmpIO%SMOISEQ)) &
         ierr = nf90_put_var(ncid, id_smoiseq, NoahmpIO%SMOISEQ, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/))

      ! --- snow layers
      ierr = nf90_put_var(ncid, id_tsno,  NoahmpIO%TSNOXY,  start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW,count(2)/))
      ierr = nf90_put_var(ncid, id_snice, NoahmpIO%SNICEXY, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW,count(2)/))
      ierr = nf90_put_var(ncid, id_snliq, NoahmpIO%SNLIQXY, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW,count(2)/))
      ierr = nf90_put_var(ncid, id_zsnso, NoahmpIO%ZSNSOXY, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW+NoahmpIO%NSOIL,count(2)/))

      ! --- snowpack scalars
      ierr = nf90_put_var(ncid, id_snow,   NoahmpIO%SNOW,    start=start, count=count)
      ierr = nf90_put_var(ncid, id_snowh,  NoahmpIO%SNOWH,   start=start, count=count)
      ierr = nf90_put_var(ncid, id_snowc,  NoahmpIO%SNOWC,   start=start, count=count)
      ierr = nf90_put_var(ncid, id_isnow,  NoahmpIO%ISNOWXY, start=start, count=count)
      ierr = nf90_put_var(ncid, id_canwat, NoahmpIO%CANWAT,  start=start, count=count)
      ierr = nf90_put_var(ncid, id_acsnom, NoahmpIO%ACSNOM,  start=start, count=count)
      ierr = nf90_put_var(ncid, id_acsnow, NoahmpIO%ACSNOW,  start=start, count=count)

      ! --- canopy / surface
      ierr = nf90_put_var(ncid, id_tv,     NoahmpIO%TVXY,     start=start, count=count)
      ierr = nf90_put_var(ncid, id_tg,     NoahmpIO%TGXY,     start=start, count=count)
      ierr = nf90_put_var(ncid, id_canice, NoahmpIO%CANICEXY, start=start, count=count)
      ierr = nf90_put_var(ncid, id_canliq, NoahmpIO%CANLIQXY, start=start, count=count)
      ierr = nf90_put_var(ncid, id_eah,    NoahmpIO%EAHXY,    start=start, count=count)
      ierr = nf90_put_var(ncid, id_tah,    NoahmpIO%TAHXY,    start=start, count=count)
      ierr = nf90_put_var(ncid, id_cm,     NoahmpIO%CMXY,     start=start, count=count)
      ierr = nf90_put_var(ncid, id_ch,     NoahmpIO%CHXY,     start=start, count=count)
      ierr = nf90_put_var(ncid, id_fwet,   NoahmpIO%FWETXY,   start=start, count=count)
      ierr = nf90_put_var(ncid, id_qsfc,   NoahmpIO%QSFC,     start=start, count=count)
      ierr = nf90_put_var(ncid, id_tsk,    NoahmpIO%TSK,      start=start, count=count)
      ierr = nf90_put_var(ncid, id_qsnow,  NoahmpIO%QSNOWXY,  start=start, count=count)
      ierr = nf90_put_var(ncid, id_qrain,  NoahmpIO%QRAINXY,  start=start, count=count)

      ! --- albedo history
      ierr = nf90_put_var(ncid, id_sneqvo, NoahmpIO%SNEQVOXY, start=start, count=count)
      ierr = nf90_put_var(ncid, id_albold, NoahmpIO%ALBOLDXY, start=start, count=count)
      ierr = nf90_put_var(ncid, id_tauss,  NoahmpIO%TAUSSXY,  start=start, count=count)
      ierr = nf90_put_var(ncid, id_albedo, NoahmpIO%ALBEDO,   start=start, count=count)

      ! --- aquifer / groundwater
      ierr = nf90_put_var(ncid, id_zwt,      NoahmpIO%ZWTXY,      start=start, count=count)
      ierr = nf90_put_var(ncid, id_wa,       NoahmpIO%WAXY,       start=start, count=count)
      ierr = nf90_put_var(ncid, id_wt,       NoahmpIO%WTXY,       start=start, count=count)
      ierr = nf90_put_var(ncid, id_smcwtd,   NoahmpIO%SMCWTDXY,   start=start, count=count)
      ierr = nf90_put_var(ncid, id_deeprech, NoahmpIO%DEEPRECHXY, start=start, count=count)
      ierr = nf90_put_var(ncid, id_rech,     NoahmpIO%RECHXY,     start=start, count=count)

      ! --- phenology
      ierr = nf90_put_var(ncid, id_lai,    NoahmpIO%LAI,    start=start, count=count)
      ierr = nf90_put_var(ncid, id_xsai,   NoahmpIO%XSAIXY, start=start, count=count)

      ! --- accumulators / carried state
      ierr = nf90_put_var(ncid, id_sfcrunoff, NoahmpIO%SFCRUNOFF, start=start, count=count)
      ierr = nf90_put_var(ncid, id_udrunoff,  NoahmpIO%UDRUNOFF,  start=start, count=count)
      ierr = nf90_put_var(ncid, id_smstav,    NoahmpIO%SMSTAV,    start=start, count=count)
      ierr = nf90_put_var(ncid, id_smstot,    NoahmpIO%SMSTOT,    start=start, count=count)
      ierr = nf90_put_var(ncid, id_emiss,     NoahmpIO%EMISS,     start=start, count=count)
      ierr = nf90_put_var(ncid, id_grdflx,    NoahmpIO%GRDFLX,    start=start, count=count)

      ! --- optional carbon / lake
      if (allocated(NoahmpIO%LFMASSXY)) ierr = nf90_put_var(ncid, id_lfmass, NoahmpIO%LFMASSXY, start=start, count=count)
      if (allocated(NoahmpIO%RTMASSXY)) ierr = nf90_put_var(ncid, id_rtmass, NoahmpIO%RTMASSXY, start=start, count=count)
      if (allocated(NoahmpIO%STMASSXY)) ierr = nf90_put_var(ncid, id_stmass, NoahmpIO%STMASSXY, start=start, count=count)
      if (allocated(NoahmpIO%WOODXY))   ierr = nf90_put_var(ncid, id_wood,   NoahmpIO%WOODXY,   start=start, count=count)
      if (allocated(NoahmpIO%GRAINXY))  ierr = nf90_put_var(ncid, id_grain,  NoahmpIO%GRAINXY,  start=start, count=count)
      if (allocated(NoahmpIO%GDDXY))    ierr = nf90_put_var(ncid, id_gdd,    NoahmpIO%GDDXY,    start=start, count=count)
      if (allocated(NoahmpIO%WSLAKEXY)) ierr = nf90_put_var(ncid, id_wslake, NoahmpIO%WSLAKEXY, start=start, count=count)

      if (NoahmpIO%blkid == (maxblocks-1)) then
         ierr = nf90_close(ncid)
      end if

   end subroutine NoahmpWriteRestart

end module NoahmpWriteRestartMod
