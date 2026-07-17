module NoahmpWriteRestartMod

! Write the full NoahMP prognostic state to a NetCDF restart file for bit-exact
! ERF restart, using the collective per-block MPI-IO pattern of NoahmpWriteLandMod.
! State is serialized at working precision (NF90_DOUBLE when kind_noahmp==8, else
! NF90_REAL); ISNOWXY is NF90_INT (needed to interpret the negative-indexed
! snow-layer arrays). Companion reader: NoahmpReadRestartMod.

   use mpi
   use netcdf
   use Machine, only : kind_noahmp
   use NoahmpIOVarType
   use NoahmpFatalMod, only : check_nc, NoahmpIO_abort

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
   ! varids -- soil-cycle accumulators (carried across steps when SOIL_UPDATE_STEPS>1)
   integer, save, private :: id_acc_ssoil, id_acc_qinsur, id_acc_qseva, id_acc_dwater, &
                             id_acc_prcp, id_acc_ecan, id_acc_etran, id_acc_edir, &
                             id_acc_etrani, id_acc_glaflw

contains

   subroutine NoahmpWriteRestart(NoahmpIO, dir, maxblocks)

      implicit none

      type(NoahmpIO_type), intent(inout) :: NoahmpIO
      character(len=*),    intent(in)    :: dir
      integer,             intent(in)    :: maxblocks

      integer :: ierr, start(2), count(2)
      integer :: nx, ny, nsoil_d, nsnow_d, nsnso_d
      integer :: rtype
      character(len=1)   :: lev_str
      character(len=512) :: filename
      logical :: ex

      ! Match the NetCDF real type to the in-memory kind for a bit-exact round-trip.
      rtype = NF90_REAL
      if (kind_noahmp == 8) rtype = NF90_DOUBLE

      if (NoahmpIO%blkid == 0) then
         write (lev_str, '(I1.1)') NoahmpIO%LEVEL

         inquire (file=trim(dir), exist=ex)
         if (.not. ex) then
            call execute_command_line("mkdir -p "//trim(dir), exitstat=ierr)
            if (ierr /= 0) then
               print *, "NoahmpWriteRestart: failed to create directory: ", trim(dir)
               call NoahmpIO_abort()
            end if
         end if

         filename = trim(dir)//"/Level_"//trim(lev_str)//".nc"
         call check_nc(nf90_create(trim(filename), IOR(NF90_CLOBBER, IOR(NF90_NETCDF4, NF90_MPIIO)), &
                       ncid, comm=NoahmpIO%comm, info=MPI_INFO_NULL), "create "//trim(filename))

         ! Dimensions
         call check_nc(nf90_def_dim(ncid, "NX",    NoahmpIO%xsglobal, nx),      "def_dim NX")
         call check_nc(nf90_def_dim(ncid, "NY",    NoahmpIO%ysglobal, ny),      "def_dim NY")
         call check_nc(nf90_def_dim(ncid, "NSOIL", NoahmpIO%NSOIL,    nsoil_d), "def_dim NSOIL")
         call check_nc(nf90_def_dim(ncid, "NSNOW", NoahmpIO%NSNOW,    nsnow_d), "def_dim NSNOW")
         call check_nc(nf90_def_dim(ncid, "NSNSO", NoahmpIO%NSNOW+NoahmpIO%NSOIL, nsnso_d), "def_dim NSNSO")

         ! Layer counts as global attributes for the read-side assert.
         call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "NSOIL", NoahmpIO%NSOIL), "put_att NSOIL")
         call check_nc(nf90_put_att(ncid, NF90_GLOBAL, "NSNOW", NoahmpIO%NSNOW), "put_att NSNOW")

         ! --- soil (NX, NSOIL, NY)
         call check_nc(nf90_def_var(ncid, "TSLB",     rtype, (/nx, nsoil_d, ny/), id_tslb),  "def_var TSLB")
         call check_nc(nf90_def_var(ncid, "SMOIS",    rtype, (/nx, nsoil_d, ny/), id_smois), "def_var SMOIS")
         call check_nc(nf90_def_var(ncid, "SH2O",     rtype, (/nx, nsoil_d, ny/), id_sh2o),  "def_var SH2O")
         if (allocated(NoahmpIO%SMOISEQ)) &
            call check_nc(nf90_def_var(ncid, "SMOISEQ", rtype, (/nx, nsoil_d, ny/), id_smoiseq), "def_var SMOISEQ")

         ! --- snow layers (NX, NSNOW, NY) and snow+soil (NX, NSNSO, NY)
         call check_nc(nf90_def_var(ncid, "TSNOXY",  rtype, (/nx, nsnow_d, ny/), id_tsno),  "def_var TSNOXY")
         call check_nc(nf90_def_var(ncid, "SNICEXY", rtype, (/nx, nsnow_d, ny/), id_snice), "def_var SNICEXY")
         call check_nc(nf90_def_var(ncid, "SNLIQXY", rtype, (/nx, nsnow_d, ny/), id_snliq), "def_var SNLIQXY")
         call check_nc(nf90_def_var(ncid, "ZSNSOXY", rtype, (/nx, nsnso_d, ny/), id_zsnso), "def_var ZSNSOXY")

         ! --- snowpack scalars (NX, NY)
         call check_nc(nf90_def_var(ncid, "SNOW",    rtype,     (/nx, ny/), id_snow),   "def_var SNOW")
         call check_nc(nf90_def_var(ncid, "SNOWH",   rtype,     (/nx, ny/), id_snowh),  "def_var SNOWH")
         call check_nc(nf90_def_var(ncid, "SNOWC",   rtype,     (/nx, ny/), id_snowc),  "def_var SNOWC")
         call check_nc(nf90_def_var(ncid, "ISNOWXY", NF90_INT,  (/nx, ny/), id_isnow),  "def_var ISNOWXY")
         call check_nc(nf90_def_var(ncid, "CANWAT",  rtype,     (/nx, ny/), id_canwat), "def_var CANWAT")
         call check_nc(nf90_def_var(ncid, "ACSNOM",  rtype,     (/nx, ny/), id_acsnom), "def_var ACSNOM")
         call check_nc(nf90_def_var(ncid, "ACSNOW",  rtype,     (/nx, ny/), id_acsnow), "def_var ACSNOW")

         ! --- canopy / surface
         call check_nc(nf90_def_var(ncid, "TVXY",    rtype, (/nx, ny/), id_tv),     "def_var TVXY")
         call check_nc(nf90_def_var(ncid, "TGXY",    rtype, (/nx, ny/), id_tg),     "def_var TGXY")
         call check_nc(nf90_def_var(ncid, "CANICEXY",rtype, (/nx, ny/), id_canice), "def_var CANICEXY")
         call check_nc(nf90_def_var(ncid, "CANLIQXY",rtype, (/nx, ny/), id_canliq), "def_var CANLIQXY")
         call check_nc(nf90_def_var(ncid, "EAHXY",   rtype, (/nx, ny/), id_eah),    "def_var EAHXY")
         call check_nc(nf90_def_var(ncid, "TAHXY",   rtype, (/nx, ny/), id_tah),    "def_var TAHXY")
         call check_nc(nf90_def_var(ncid, "CMXY",    rtype, (/nx, ny/), id_cm),     "def_var CMXY")
         call check_nc(nf90_def_var(ncid, "CHXY",    rtype, (/nx, ny/), id_ch),     "def_var CHXY")
         call check_nc(nf90_def_var(ncid, "FWETXY",  rtype, (/nx, ny/), id_fwet),   "def_var FWETXY")
         call check_nc(nf90_def_var(ncid, "QSFC",    rtype, (/nx, ny/), id_qsfc),   "def_var QSFC")
         ! TSK, EMISS, WSLAKEXY use c_kind_noahmp == kind_noahmp, so rtype applies.
         call check_nc(nf90_def_var(ncid, "TSK",     rtype, (/nx, ny/), id_tsk),    "def_var TSK")
         call check_nc(nf90_def_var(ncid, "QSNOWXY", rtype, (/nx, ny/), id_qsnow),  "def_var QSNOWXY")
         call check_nc(nf90_def_var(ncid, "QRAINXY", rtype, (/nx, ny/), id_qrain),  "def_var QRAINXY")

         ! --- albedo history
         call check_nc(nf90_def_var(ncid, "SNEQVOXY",rtype, (/nx, ny/), id_sneqvo), "def_var SNEQVOXY")
         call check_nc(nf90_def_var(ncid, "ALBOLDXY",rtype, (/nx, ny/), id_albold), "def_var ALBOLDXY")
         call check_nc(nf90_def_var(ncid, "TAUSSXY", rtype, (/nx, ny/), id_tauss),  "def_var TAUSSXY")
         call check_nc(nf90_def_var(ncid, "ALBEDO",  rtype, (/nx, ny/), id_albedo), "def_var ALBEDO")

         ! --- aquifer / groundwater
         call check_nc(nf90_def_var(ncid, "ZWTXY",     rtype, (/nx, ny/), id_zwt),      "def_var ZWTXY")
         call check_nc(nf90_def_var(ncid, "WAXY",      rtype, (/nx, ny/), id_wa),       "def_var WAXY")
         call check_nc(nf90_def_var(ncid, "WTXY",      rtype, (/nx, ny/), id_wt),       "def_var WTXY")
         call check_nc(nf90_def_var(ncid, "SMCWTDXY",  rtype, (/nx, ny/), id_smcwtd),   "def_var SMCWTDXY")
         call check_nc(nf90_def_var(ncid, "DEEPRECHXY",rtype, (/nx, ny/), id_deeprech), "def_var DEEPRECHXY")
         call check_nc(nf90_def_var(ncid, "RECHXY",    rtype, (/nx, ny/), id_rech),     "def_var RECHXY")

         ! --- phenology
         call check_nc(nf90_def_var(ncid, "LAI",     rtype, (/nx, ny/), id_lai),  "def_var LAI")
         call check_nc(nf90_def_var(ncid, "XSAIXY",  rtype, (/nx, ny/), id_xsai), "def_var XSAIXY")

         ! --- accumulators / carried state
         call check_nc(nf90_def_var(ncid, "SFCRUNOFF",rtype, (/nx, ny/), id_sfcrunoff), "def_var SFCRUNOFF")
         call check_nc(nf90_def_var(ncid, "UDRUNOFF", rtype, (/nx, ny/), id_udrunoff),  "def_var UDRUNOFF")
         call check_nc(nf90_def_var(ncid, "SMSTAV",   rtype, (/nx, ny/), id_smstav),    "def_var SMSTAV")
         call check_nc(nf90_def_var(ncid, "SMSTOT",   rtype, (/nx, ny/), id_smstot),    "def_var SMSTOT")
         call check_nc(nf90_def_var(ncid, "EMISS",    rtype, (/nx, ny/), id_emiss),     "def_var EMISS")
         call check_nc(nf90_def_var(ncid, "GRDFLX",   rtype, (/nx, ny/), id_grdflx),    "def_var GRDFLX")

         ! --- soil-cycle accumulators (mid-cycle carry when SOIL_UPDATE_STEPS>1)
         call check_nc(nf90_def_var(ncid, "ACC_SSOILXY", rtype, (/nx, ny/), id_acc_ssoil),  "def_var ACC_SSOILXY")
         call check_nc(nf90_def_var(ncid, "ACC_QINSURXY",rtype, (/nx, ny/), id_acc_qinsur), "def_var ACC_QINSURXY")
         call check_nc(nf90_def_var(ncid, "ACC_QSEVAXY", rtype, (/nx, ny/), id_acc_qseva),  "def_var ACC_QSEVAXY")
         call check_nc(nf90_def_var(ncid, "ACC_DWATERXY",rtype, (/nx, ny/), id_acc_dwater), "def_var ACC_DWATERXY")
         call check_nc(nf90_def_var(ncid, "ACC_PRCPXY",  rtype, (/nx, ny/), id_acc_prcp),   "def_var ACC_PRCPXY")
         call check_nc(nf90_def_var(ncid, "ACC_ECANXY",  rtype, (/nx, ny/), id_acc_ecan),   "def_var ACC_ECANXY")
         call check_nc(nf90_def_var(ncid, "ACC_ETRANXY", rtype, (/nx, ny/), id_acc_etran),  "def_var ACC_ETRANXY")
         call check_nc(nf90_def_var(ncid, "ACC_EDIRXY",  rtype, (/nx, ny/), id_acc_edir),   "def_var ACC_EDIRXY")
         call check_nc(nf90_def_var(ncid, "ACC_ETRANIXY",rtype, (/nx, nsoil_d, ny/), id_acc_etrani), "def_var ACC_ETRANIXY")
         call check_nc(nf90_def_var(ncid, "ACC_GLAFLWXY",rtype, (/nx, ny/), id_acc_glaflw), "def_var ACC_GLAFLWXY")

         ! --- optional carbon / dveg (only if allocated)
         if (allocated(NoahmpIO%LFMASSXY)) &
            call check_nc(nf90_def_var(ncid, "LFMASSXY", rtype, (/nx, ny/), id_lfmass), "def_var LFMASSXY")
         if (allocated(NoahmpIO%RTMASSXY)) &
            call check_nc(nf90_def_var(ncid, "RTMASSXY", rtype, (/nx, ny/), id_rtmass), "def_var RTMASSXY")
         if (allocated(NoahmpIO%STMASSXY)) &
            call check_nc(nf90_def_var(ncid, "STMASSXY", rtype, (/nx, ny/), id_stmass), "def_var STMASSXY")
         if (allocated(NoahmpIO%WOODXY)) &
            call check_nc(nf90_def_var(ncid, "WOODXY",   rtype, (/nx, ny/), id_wood),   "def_var WOODXY")
         if (allocated(NoahmpIO%GRAINXY)) &
            call check_nc(nf90_def_var(ncid, "GRAINXY",  rtype, (/nx, ny/), id_grain),  "def_var GRAINXY")
         if (allocated(NoahmpIO%GDDXY)) &
            call check_nc(nf90_def_var(ncid, "GDDXY",    rtype, (/nx, ny/), id_gdd),    "def_var GDDXY")
         if (allocated(NoahmpIO%WSLAKEXY)) &
            call check_nc(nf90_def_var(ncid, "WSLAKEXY", rtype, (/nx, ny/), id_wslake), "def_var WSLAKEXY")

         call check_nc(nf90_enddef(ncid), "enddef")
      end if

      ! Hyperslab for this block within the global domain.
      start = (/NoahmpIO%xstart-NoahmpIO%xoffset+1, NoahmpIO%ystart-NoahmpIO%yoffset+1/)
      count = (/NoahmpIO%xend-NoahmpIO%xstart+1,    NoahmpIO%yend-NoahmpIO%ystart+1/)

      ! --- soil
      call check_nc(nf90_put_var(ncid, id_tslb,  NoahmpIO%TSLB,  start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/)), "put_var TSLB")
      call check_nc(nf90_put_var(ncid, id_smois, NoahmpIO%SMOIS, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/)), "put_var SMOIS")
      call check_nc(nf90_put_var(ncid, id_sh2o,  NoahmpIO%SH2O,  start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/)), "put_var SH2O")
      if (allocated(NoahmpIO%SMOISEQ)) &
         call check_nc(nf90_put_var(ncid, id_smoiseq, NoahmpIO%SMOISEQ, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/)), "put_var SMOISEQ")

      ! --- snow layers
      call check_nc(nf90_put_var(ncid, id_tsno,  NoahmpIO%TSNOXY,  start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW,count(2)/)), "put_var TSNOXY")
      call check_nc(nf90_put_var(ncid, id_snice, NoahmpIO%SNICEXY, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW,count(2)/)), "put_var SNICEXY")
      call check_nc(nf90_put_var(ncid, id_snliq, NoahmpIO%SNLIQXY, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW,count(2)/)), "put_var SNLIQXY")
      call check_nc(nf90_put_var(ncid, id_zsnso, NoahmpIO%ZSNSOXY, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSNOW+NoahmpIO%NSOIL,count(2)/)), "put_var ZSNSOXY")

      ! --- snowpack scalars
      call check_nc(nf90_put_var(ncid, id_snow,   NoahmpIO%SNOW,    start=start, count=count), "put_var SNOW")
      call check_nc(nf90_put_var(ncid, id_snowh,  NoahmpIO%SNOWH,   start=start, count=count), "put_var SNOWH")
      call check_nc(nf90_put_var(ncid, id_snowc,  NoahmpIO%SNOWC,   start=start, count=count), "put_var SNOWC")
      call check_nc(nf90_put_var(ncid, id_isnow,  NoahmpIO%ISNOWXY, start=start, count=count), "put_var ISNOWXY")
      call check_nc(nf90_put_var(ncid, id_canwat, NoahmpIO%CANWAT,  start=start, count=count), "put_var CANWAT")
      call check_nc(nf90_put_var(ncid, id_acsnom, NoahmpIO%ACSNOM,  start=start, count=count), "put_var ACSNOM")
      call check_nc(nf90_put_var(ncid, id_acsnow, NoahmpIO%ACSNOW,  start=start, count=count), "put_var ACSNOW")

      ! --- canopy / surface
      call check_nc(nf90_put_var(ncid, id_tv,     NoahmpIO%TVXY,     start=start, count=count), "put_var TVXY")
      call check_nc(nf90_put_var(ncid, id_tg,     NoahmpIO%TGXY,     start=start, count=count), "put_var TGXY")
      call check_nc(nf90_put_var(ncid, id_canice, NoahmpIO%CANICEXY, start=start, count=count), "put_var CANICEXY")
      call check_nc(nf90_put_var(ncid, id_canliq, NoahmpIO%CANLIQXY, start=start, count=count), "put_var CANLIQXY")
      call check_nc(nf90_put_var(ncid, id_eah,    NoahmpIO%EAHXY,    start=start, count=count), "put_var EAHXY")
      call check_nc(nf90_put_var(ncid, id_tah,    NoahmpIO%TAHXY,    start=start, count=count), "put_var TAHXY")
      call check_nc(nf90_put_var(ncid, id_cm,     NoahmpIO%CMXY,     start=start, count=count), "put_var CMXY")
      call check_nc(nf90_put_var(ncid, id_ch,     NoahmpIO%CHXY,     start=start, count=count), "put_var CHXY")
      call check_nc(nf90_put_var(ncid, id_fwet,   NoahmpIO%FWETXY,   start=start, count=count), "put_var FWETXY")
      call check_nc(nf90_put_var(ncid, id_qsfc,   NoahmpIO%QSFC,     start=start, count=count), "put_var QSFC")
      call check_nc(nf90_put_var(ncid, id_tsk,    NoahmpIO%TSK,      start=start, count=count), "put_var TSK")
      call check_nc(nf90_put_var(ncid, id_qsnow,  NoahmpIO%QSNOWXY,  start=start, count=count), "put_var QSNOWXY")
      call check_nc(nf90_put_var(ncid, id_qrain,  NoahmpIO%QRAINXY,  start=start, count=count), "put_var QRAINXY")

      ! --- albedo history
      call check_nc(nf90_put_var(ncid, id_sneqvo, NoahmpIO%SNEQVOXY, start=start, count=count), "put_var SNEQVOXY")
      call check_nc(nf90_put_var(ncid, id_albold, NoahmpIO%ALBOLDXY, start=start, count=count), "put_var ALBOLDXY")
      call check_nc(nf90_put_var(ncid, id_tauss,  NoahmpIO%TAUSSXY,  start=start, count=count), "put_var TAUSSXY")
      call check_nc(nf90_put_var(ncid, id_albedo, NoahmpIO%ALBEDO,   start=start, count=count), "put_var ALBEDO")

      ! --- aquifer / groundwater
      call check_nc(nf90_put_var(ncid, id_zwt,      NoahmpIO%ZWTXY,      start=start, count=count), "put_var ZWTXY")
      call check_nc(nf90_put_var(ncid, id_wa,       NoahmpIO%WAXY,       start=start, count=count), "put_var WAXY")
      call check_nc(nf90_put_var(ncid, id_wt,       NoahmpIO%WTXY,       start=start, count=count), "put_var WTXY")
      call check_nc(nf90_put_var(ncid, id_smcwtd,   NoahmpIO%SMCWTDXY,   start=start, count=count), "put_var SMCWTDXY")
      call check_nc(nf90_put_var(ncid, id_deeprech, NoahmpIO%DEEPRECHXY, start=start, count=count), "put_var DEEPRECHXY")
      call check_nc(nf90_put_var(ncid, id_rech,     NoahmpIO%RECHXY,     start=start, count=count), "put_var RECHXY")

      ! --- phenology
      call check_nc(nf90_put_var(ncid, id_lai,    NoahmpIO%LAI,    start=start, count=count), "put_var LAI")
      call check_nc(nf90_put_var(ncid, id_xsai,   NoahmpIO%XSAIXY, start=start, count=count), "put_var XSAIXY")

      ! --- accumulators / carried state
      call check_nc(nf90_put_var(ncid, id_sfcrunoff, NoahmpIO%SFCRUNOFF, start=start, count=count), "put_var SFCRUNOFF")
      call check_nc(nf90_put_var(ncid, id_udrunoff,  NoahmpIO%UDRUNOFF,  start=start, count=count), "put_var UDRUNOFF")
      call check_nc(nf90_put_var(ncid, id_smstav,    NoahmpIO%SMSTAV,    start=start, count=count), "put_var SMSTAV")
      call check_nc(nf90_put_var(ncid, id_smstot,    NoahmpIO%SMSTOT,    start=start, count=count), "put_var SMSTOT")
      call check_nc(nf90_put_var(ncid, id_emiss,     NoahmpIO%EMISS,     start=start, count=count), "put_var EMISS")
      call check_nc(nf90_put_var(ncid, id_grdflx,    NoahmpIO%GRDFLX,    start=start, count=count), "put_var GRDFLX")

      ! --- soil-cycle accumulators
      call check_nc(nf90_put_var(ncid, id_acc_ssoil,  NoahmpIO%ACC_SSOILXY,  start=start, count=count), "put_var ACC_SSOILXY")
      call check_nc(nf90_put_var(ncid, id_acc_qinsur, NoahmpIO%ACC_QINSURXY, start=start, count=count), "put_var ACC_QINSURXY")
      call check_nc(nf90_put_var(ncid, id_acc_qseva,  NoahmpIO%ACC_QSEVAXY,  start=start, count=count), "put_var ACC_QSEVAXY")
      call check_nc(nf90_put_var(ncid, id_acc_dwater, NoahmpIO%ACC_DWATERXY, start=start, count=count), "put_var ACC_DWATERXY")
      call check_nc(nf90_put_var(ncid, id_acc_prcp,   NoahmpIO%ACC_PRCPXY,   start=start, count=count), "put_var ACC_PRCPXY")
      call check_nc(nf90_put_var(ncid, id_acc_ecan,   NoahmpIO%ACC_ECANXY,   start=start, count=count), "put_var ACC_ECANXY")
      call check_nc(nf90_put_var(ncid, id_acc_etran,  NoahmpIO%ACC_ETRANXY,  start=start, count=count), "put_var ACC_ETRANXY")
      call check_nc(nf90_put_var(ncid, id_acc_edir,   NoahmpIO%ACC_EDIRXY,   start=start, count=count), "put_var ACC_EDIRXY")
      call check_nc(nf90_put_var(ncid, id_acc_etrani, NoahmpIO%ACC_ETRANIXY, start=(/start(1),1,start(2)/), count=(/count(1),NoahmpIO%NSOIL,count(2)/)), "put_var ACC_ETRANIXY")
      call check_nc(nf90_put_var(ncid, id_acc_glaflw, NoahmpIO%ACC_GLAFLWXY, start=start, count=count), "put_var ACC_GLAFLWXY")

      ! --- optional carbon / lake
      if (allocated(NoahmpIO%LFMASSXY)) call check_nc(nf90_put_var(ncid, id_lfmass, NoahmpIO%LFMASSXY, start=start, count=count), "put_var LFMASSXY")
      if (allocated(NoahmpIO%RTMASSXY)) call check_nc(nf90_put_var(ncid, id_rtmass, NoahmpIO%RTMASSXY, start=start, count=count), "put_var RTMASSXY")
      if (allocated(NoahmpIO%STMASSXY)) call check_nc(nf90_put_var(ncid, id_stmass, NoahmpIO%STMASSXY, start=start, count=count), "put_var STMASSXY")
      if (allocated(NoahmpIO%WOODXY))   call check_nc(nf90_put_var(ncid, id_wood,   NoahmpIO%WOODXY,   start=start, count=count), "put_var WOODXY")
      if (allocated(NoahmpIO%GRAINXY))  call check_nc(nf90_put_var(ncid, id_grain,  NoahmpIO%GRAINXY,  start=start, count=count), "put_var GRAINXY")
      if (allocated(NoahmpIO%GDDXY))    call check_nc(nf90_put_var(ncid, id_gdd,    NoahmpIO%GDDXY,    start=start, count=count), "put_var GDDXY")
      if (allocated(NoahmpIO%WSLAKEXY)) call check_nc(nf90_put_var(ncid, id_wslake, NoahmpIO%WSLAKEXY, start=start, count=count), "put_var WSLAKEXY")

      if (NoahmpIO%blkid == (maxblocks-1)) then
         call check_nc(nf90_close(ncid), "close")
      end if

   end subroutine NoahmpWriteRestart

end module NoahmpWriteRestartMod
