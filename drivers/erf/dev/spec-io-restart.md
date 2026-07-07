# Spec: checkpoint / restart I/O

> Status: living document · Owns: `NoahmpWriteRestartMod.F90` /
> `NoahmpReadRestartMod.F90` — full prognostic-state checkpoint and restore.
> Companion: [`spec-io-parallel.md`](spec-io-parallel.md) (shares the collective
> write pattern). Origin: ERF #3255 / ERF restart issue #3255.

## 1. Goal

Serialize the **full Noah-MP prognostic state** so an ERF restart reproduces a
cold-start trajectory **bitwise**. This is distinct from the land output
([`spec-io-parallel.md`](spec-io-parallel.md)) in three ways:

| | Land output | Restart |
|---|---|---|
| Scope | small diagnostic/coupling subset | full prognostic state |
| Precision | `NF90_FLOAT` (lossy) | `NF90_DOUBLE` when `kind_noahmp==8`, else `NF90_REAL` (**bit-exact**) |
| Cadence | every output step | at checkpoints |

Entry points: `NoahmpIO_type::WriteRestart(dir)` / `ReadRestart(dir)` →
`NoahmpWriteRestart_fi` / `NoahmpReadRestart_fi` →
`NoahmpWriteRestart` / `NoahmpReadRestart(NoahmpIO, dir, maxblocks)`.

## 2. On-disk layout

One file per AMR level under the caller-supplied directory:

```
<dir>/Level_<L>.nc           # L = 1-digit level
```

Same collective parallel-NetCDF-4 pattern as the land output (open/define on
`blkid==0`, each block writes its hyperslab, close on `blkid==maxblocks-1`;
`comm=NoahmpIO%comm`, `NF90_MPIIO`). The string `dir` crosses the boundary as
`(const char*, int len)` with an overflow-guarded length (see
[`spec-fc-api.md`](spec-fc-api.md) §3).

### NetCDF real type — precision exactness

```fortran
rtype = NF90_REAL
if (kind_noahmp == 8) rtype = NF90_DOUBLE
```

so the round trip has **no float↔double conversion**. The few fields declared
with the C-boundary kind `c_kind_noahmp` (`TSK`, `EMISS`, `WSLAKEXY`) use the
same `rtype` because `c_kind_noahmp == kind_noahmp` in every build; the reader
marks them with a distinct helper (`get2dd`) only to flag them as boundary
fields.

### Dimensions and layer-geometry guard

| Dim | Size |
|-----|------|
| `NX` | `xsglobal` |
| `NY` | `ysglobal` |
| `NSOIL` | `NSOIL` |
| `NSNOW` | `NSNOW` |
| `NSNSO` | `NSNOW + NSOIL` (combined snow+soil column) |

`NSOIL` and `NSNOW` are also written as **global attributes**. On read,
`NoahmpReadRestart` asserts the file's `NSOIL`/`NSNOW` match the run and aborts on
mismatch — a different layer geometry means a corrupt restore.

## 3. State captured

Grouped as in the source (all per-block hyperslabs):

- **Soil** `(NX,NSOIL,NY)`: `TSLB, SMOIS, SH2O`, and `SMOISEQ` (optional).
- **Snow layers** `(NX,NSNOW,NY)`: `TSNOXY, SNICEXY, SNLIQXY`; `ZSNSOXY`
  `(NX,NSNSO,NY)`.
- **Snowpack scalars** `(NX,NY)`: `SNOW, SNOWH, SNOWC, CANWAT, ACSNOM, ACSNOW`,
  and **`ISNOWXY` as `NF90_INT`** — the active snow-layer count, required to
  interpret the negative-indexed snow-layer arrays on restore.
- **Canopy / surface**: `TVXY, TGXY, CANICEXY, CANLIQXY, EAHXY, TAHXY, CMXY,
  CHXY, FWETXY, QSFC, TSK, QSNOWXY, QRAINXY`.
- **Albedo history**: `SNEQVOXY, ALBOLDXY, TAUSSXY, ALBEDO`.
- **Aquifer / groundwater**: `ZWTXY, WAXY, WTXY, SMCWTDXY, DEEPRECHXY, RECHXY`.
- **Phenology**: `LAI, XSAIXY`.
- **Accumulators / carried state**: `SFCRUNOFF, UDRUNOFF, SMSTAV, SMSTOT, EMISS,
  GRDFLX`.
- **Optional carbon / dveg / lake** (only if allocated): `LFMASSXY, RTMASSXY,
  STMASSXY, WOODXY, GRAINXY, GDDXY, WSLAKEXY`.

### Optional-field discipline

- **Write**: every conditionally-allocated field is guarded with `allocated(...)`
  on both `nf90_def_var` and `nf90_put_var`, so a build that doesn't allocate it
  simply omits it.
- **Read**: a field is read only if it is **both** present in the file *and*
  allocated in `NoahmpIO`. Each `get*` helper takes a `required` flag:
  - `required=.true.` → a missing varid aborts (mandatory state).
  - `required=.false.` → a missing varid is **silently skipped**, preserving the
    cold-init value (optional state).

Helpers: `get2d` (real `kind_noahmp`), `get2dd` (real `c_kind_noahmp`, the
boundary fields), `get2di` (integer, for `ISNOWXY`), `get3d` (real, with a layer
count `nk`).

## 4. Restart ordering (how a bitwise restart happens)

`NoahmpReadRestart` runs **after** the normal cold init
(`ReadLandMain` / `InitMain`), so it *overwrites* the table/`wrfinput`-initialized
state with the checkpointed values. On the first `DriverMain`, the per-step
In-transfer (`*VarInTransfer`) pulls those values into the 1-D physics state,
giving a bitwise-identical continuation. Host sequence per block:

```
ScalarInitDefault → ReadNamelist → ReadTable → ReadLandHeader → ReadLandMain
  → VarInitDefault → InitMain
  → ReadRestart(dir)          # only on a restart; overwrites cold-init state
  → (time loop) DriverMain ...
```

### 4.1 Host↔device reconciliation (once the GPU offload lands)

Restart state is **device-resident** under the offload
([`plan-cpp-interface.md`](plan-cpp-interface.md)), so the two `update`s must
bracket the host-side NetCDF I/O:

- **Write** — refresh host from device *before* `WriteRestart` serializes it:
  `!$acc update host(<the full prognostic set>)`.
- **Read** — push the restored values to the device *after* `ReadRestart`
  overwrites the host arrays, *before* the first `DriverMain` reads them on device:
  `!$acc update device(<the full prognostic set>)`.

Because bit-exactness is the whole point of restart, the **host Fortran path stays
the bitwise oracle**: a restart on the GPU path is validated within the documented
GPU tolerance against it, not bit-for-bit (see
[`plan-cpp-interface.md`](plan-cpp-interface.md) §4, §5.9). Getting the `update`
placement wrong reads/writes stale state and silently corrupts the restart — see
[`spec-memory-safety.md`](spec-memory-safety.md) §7.3.

## 5. Known constraints / hardening backlog

- **Error handling is `print` + `stop`** in both modules (shared with the land
  writer). A serial `stop` can deadlock peers at the next collective; route
  through `NoahmpIO_abort()` instead. Fix the land and restart paths together
  (see [`spec-io-parallel.md`](spec-io-parallel.md) §5).
- **`ncid` is module-global `save`** — same ordering assumption as the land
  writer (blocks visited in `blkid` order, no cross-level interleave on a rank).
- **Layer geometry is the only structural assert.** Grid geometry
  (`xsglobal`/`ysglobal`) is trusted to match; a domain reshape between write and
  read is not detected. Consider asserting global extents too.
- **`mkdir -p` via `execute_command_line`** is done by `blkid==0`; a failure
  prints and `stop`s.

## 6. Adding a prognostic field to the checkpoint

The reader and writer are **hand-kept in sync** — edit both:

1. **Writer** (`NoahmpWriteRestartMod.F90`): add a `varid` to the right module
   list; `nf90_def_var` it with `rtype` (or `NF90_INT`) and the correct dims
   inside the `blkid==0` define block; `nf90_put_var` it with the hyperslab. Guard
   with `allocated()` if optional.
2. **Reader** (`NoahmpReadRestartMod.F90`): add a matching `get2d/get2dd/get2di/
   get3d` call with the correct `required` flag.
3. If the field changes the layer geometry, update the dim/attribute asserts.

Keep the two files' field lists identical (modulo the `required`/`allocated`
guards) — a field written but not read silently drops state on restart.

## 7. Acceptance checklist

- [ ] Cold run to step N, checkpoint; restart from the checkpoint; continue to
      step M. Compare against an uninterrupted cold run to step M — **bitwise
      identical**.
- [ ] `NF90_DOUBLE` in the file for a double-precision build (`ncdump -h`).
- [ ] Layer-mismatch file is rejected with a clear message.
- [ ] An optional field absent from the build round-trips without error
      (skipped on read).
- [ ] Negative-indexed snow arrays restore correctly (cross-check `ISNOWXY`).
