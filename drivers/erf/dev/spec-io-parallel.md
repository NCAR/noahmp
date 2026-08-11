# Spec: parallel land-output I/O

> Status: living document · Owns: `NoahmpWriteLandMod.F90` — the per-timestep,
> collective, parallel-NetCDF land output. · Companion:
> [`spec-io-restart.md`](spec-io-restart.md) (same I/O pattern, prognostic
> state), [`spec-overview.md`](spec-overview.md) §3 (level/block geometry).

## 1. Goal

Write the small **diagnostic / coupling** subset of Noah-MP state to NetCDF every
output step, **collectively over MPI** so every block (MPI rank's tile) writes
its own hyperslab of one shared file per AMR level — no gather-to-root, no
per-rank files to post-process.

Entry point: `NoahmpIO_type::WriteLand(int filenum)` → `NoahmpWriteLand_fi` →
`NoahmpWriteLand(NoahmpIO, filenum, maxblocks)`. `maxblocks =
SIZE(NoahmpIO_vect(level)%NoahmpIO)`.

## 2. On-disk layout

One directory per output step, one NetCDF-4 file per level:

```
lnd<NNNNN>/Level_<L>.nc      # NNNNN = zero-padded filenum (step); L = 1-digit level
```

Dimensions:

| Dim | Size | Meaning |
|-----|------|---------|
| `NX` | `xsglobal` | global grid points in x |
| `NY` | `ysglobal` | global grid points in y |
| `NSOIL` | `NSOIL` | soil layers |
| `COMP2D` | `2` | 2 spectral/albedo components |

Variable shapes — written from Fortran, so the **fastest-varying axis is first**:

- **2-D surface** `(NX, NY)`: `TERRAIN, SNOWH, SHBXY, EVBXY, VEGFRA, GVFMIN,
  GVFMAX, TSK, EMISS, SAVXY, SAGXY, PAHXY, FIRAXY, HFX, LH, GRDFLX, GHBXY,
  CANHSXY, TAU_EW, TAU_NS` plus the core surface diagnostics `T2MVXY, T2MBXY,
  Q2MVXY, Q2MBXY, TRADXY, FVEGXY, TGVXY, TGBXY, SHGXY, SHCXY, EVGXY, EVCXY,
  TRXY, RUNSFXY, RUNSBXY, ECANXY, EDIRXY, ETRANXY, FSAXY, RS, Z0, ZNT`
- **3-D soil** `(NX, NSOIL, NY)`: `TSLB, SMOIS`
- **3-D 2-component** `(NX, COMP2D, NY)`: `ALBSFCDIRXY, ALBSFCDIFXY`

All output is `NF90_FLOAT` (single precision — this is a diagnostic plotfile, not
a bit-exact checkpoint; contrast [`spec-io-restart.md`](spec-io-restart.md)).
There are currently no coordinate variables for `NX`/`NY` (index space only).

> This list is exactly what `NoahmpWriteLandMod.F90` defines today. Consumers
> should **discover** variables/dims at run time, not hardcode this table, so a
> tool keeps working when a field is added (see §6).

## 3. The collective write pattern

```
if (blkid == 0) then            ! the FIRST block opens + defines the file
   nf90_create(filename, NF90_CLOBBER | NF90_NETCDF4 | NF90_MPIIO,
               ncid, comm=NoahmpIO%comm, info=MPI_INFO_NULL)
   nf90_def_dim(...) ; nf90_def_var(...) ; nf90_enddef(ncid)
end if

start = (/ xstart-xoffset+1, ystart-yoffset+1 /)     ! this block's hyperslab origin
count = (/ xend-xstart+1,    yend-ystart+1    /)     ! this block's extent
nf90_put_var(ncid, varid, NoahmpIO%FIELD, start=start, count=count)   ! every block writes its slab

if (blkid == maxblocks-1) then  ! the LAST block closes
   nf90_close(ncid)
end if
```

Key points:

- **`NF90_MPIIO`** + `comm=NoahmpIO%comm` → parallel NetCDF-4. `comm` is the
  host's MPI communicator, passed in across the boundary (Noah-MP holds no MPI
  dependency of its own — see [`spec-memory-safety.md`](spec-memory-safety.md)).
- The file `ncid` is a **module `save` variable**: opened by `blkid==0`, reused by
  every block on the rank, closed by `blkid==maxblocks-1`. (This assumes blocks
  are visited in `blkid` order on each rank — see §5.)
- **Hyperslab math**: `start` converts this block's global `(xstart, ystart)` into
  a 1-based NetCDF offset within the global `(NX, NY)` dataset using
  `xoffset/yoffset`. `count` is the owned patch size. 3-D fields expand `start`/
  `count` to 3 entries with the layer axis in the middle
  (`(/start(1),1,start(2)/)`, `(/count(1),NLAYER,count(2)/)`).
- Directory creation is done once by `blkid==0` via `execute_command_line("mkdir
  -p …")`.

### 3.1 Device→host refresh (once the GPU offload lands)

Under the GPU offload ([`plan-cpp-interface.md`](plan-cpp-interface.md)) the
written fields are **device-resident** — the host copy `nf90_put_var` reads is
stale between steps. Before the collective write (i.e. before the `blkid==0`
define/first `put_var`), refresh host from device for the output set:

```fortran
!$acc update host(NoahmpIO%TSK, NoahmpIO%HFX, ...)   ! only the fields being written
```

The write itself stays **host-side** (parallel NetCDF over `NoahmpIO%comm`); we do
not do GPU-direct I/O. This is off the hot path (only output steps), so the copy
cost is acceptable. Omitting it does not crash — it silently writes stale data
(see [`spec-memory-safety.md`](spec-memory-safety.md) §7.3). The output subset is
lossy `NF90_FLOAT` anyway, so GPU/host round-off differences here are immaterial
(contrast the restart path, which must stay bit-exact —
[`spec-io-restart.md`](spec-io-restart.md) §4).

## 4. Adding a field to the output

`NoahmpWriteLandMod.F90` is **not** generated — edit it directly, following the
`TSK` / `SMOIS` pattern:

1. Add a NetCDF var id to the module-level
   `integer, save, private :: …` list.
2. Inside the `blkid==0` define block (before `nf90_enddef`):
   - 2-D: `nf90_def_var(ncid, "NAME", NF90_FLOAT, (/nx, ny/), nameid)`
   - 3-D soil: `(/nx, nsoil, ny/)`; 2-component: `(/nx, comp2d, ny/)`
3. After `nf90_enddef`, write it:
   - 2-D: `nf90_put_var(ncid, nameid, NoahmpIO%NAME, start=start, count=count)`
   - 3-D: expand `start`/`count` to 3 entries as `TSLB`/`SMOIS` do.

If the field is new to the boundary too, do
[`spec-add-coupled-variable.md`](spec-add-coupled-variable.md) first.

## 5. Known constraints / hardening backlog

These are **current limitations** worth tracking, not necessarily bugs:

- **NetCDF error handling is `print` + `stop`.** A serial `stop` terminates only
  the calling image and can deadlock peers in the next collective. For
  consistency with the rest of the driver these should route through
  `NoahmpIO_abort()` (see [`spec-memory-safety.md`](spec-memory-safety.md) §5).
  The restart modules share this and should be fixed together.
- **`ncid` is module-global and `save`.** Correct only if (a) blocks are visited
  in `blkid` order and (b) `WriteLand` for different levels does not interleave on
  one rank. Both hold today; a refactor that breaks either must carry `ncid` in
  state instead.
- **No coordinate variables / metadata.** Index space only; no `XLAT`/`XLONG`,
  units, or `_FillValue`. Consider adding CF-style coordinates and attributes.
- **Hardcoded single precision.** Fine for a plotfile; do not copy this choice
  into the restart path (which is precision-exact by design).

## 6. Consumer side (plotting)

A reader just opens `lnd<NNNNN>/Level_<L>.nc` with `xarray`/`netCDF4`. Because
Fortran writes fastest-axis-first, a reader sees axes as listed in §2 — transpose
2-D fields so the image reads x-horizontal / y-vertical. A plotting tool should:

- discover variables/dims at run time (do not hardcode §2);
- locate `NSOIL`/`COMP2D` **by name** and slice 3-D fields by `--layer`/
  `--component`;
- mask sentinel/undefined values (e.g. `-9999`) so they don't flatten the color
  scale;
- fail clearly on a missing variable, listing what is present.

No reference plotting tool ships in `tools/` yet (it holds only
`NoahmpMacro.py`) — build one to this contract.

## 7. Acceptance checklist

- [ ] A multi-rank run produces one `lnd<NNNNN>/Level_<L>.nc` per level per step.
- [ ] Each block's hyperslab lands in the right place (no overlap/gap) — a stitched
      field is continuous across block boundaries.
- [ ] A new field appears with correct dims and values after the §4 steps.
- [ ] File opens cleanly in `ncdump -h` / `xarray`.
