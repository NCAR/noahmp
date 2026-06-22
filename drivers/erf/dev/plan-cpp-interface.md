# Plan: port the Noah-MP ERF driver to C++ for GPU execution

> Status: **proposal / roadmap** (not yet started) · Audience: contributors
> planning the GPU effort. · Prereqs read: [`spec-overview.md`](spec-overview.md),
> [`spec-fc-api.md`](spec-fc-api.md), [`spec-memory-safety.md`](spec-memory-safety.md).

## 0. TL;DR

The **boundary** is already C++ (`NoahmpIO_type`, `NoahmpArray`,
`noahmp::fatal`). What is still Fortran and **CPU-only** is the **physics**: the
per-column solver (`NoahmpMain` / `NoahmpMainGlacier`) and the `i,j` driver loop
that packs/unpacks it (`NoahmpDriverMainMod.F90`). This plan ports that column
solver to **device-callable C++** so ERF can run it as an AMReX `ParallelFor`
over each block's box — on GPU where available, threaded CPU otherwise — while
keeping the existing ABI as the fallback and the migration's safety net.

The end state: `NoahmpIO_type::DriverMain()` dispatches to a C++ `ParallelFor`
that calls a `__host__ __device__` `noahmp_column(...)` per land cell, instead of
crossing into Fortran. The Fortran path stays buildable as the reference oracle.

## 1. Where we are vs. where we're going

```
  TODAY                                            TARGET
  ─────                                            ──────
  ERF (C++)                                        ERF (C++)
    └ NoahmpIO_type::DriverMain()                    └ NoahmpIO_type::DriverMain()
        └ NoahmpDriverMain_fi  (bind C)                  └ amrex::ParallelFor(box, [=] (i,j) {
            └ NoahmpDriverMain  (Fortran)                      noahmp_column(state, params, i, j);  // C++ device fn
                └ i,j loop                                 });
                    └ *VarInTransfer (2D→1D)            (no language boundary on the hot path)
                    └ NoahmpMain    (Fortran column)
                    └ *VarOutTransfer (1D→2D)
```

- **Reused as-is:** the AMR level/block model, `NoahmpArray{1,2,3}D` (already
  Fortran-layout, bounds-checked, `__host__ __device__`-friendly once annotated),
  the precision discipline (`noahmp_real`), the fatal-handler abstraction, the
  parallel NetCDF I/O (stays host-side, runs after a device→host copy).
- **Ported:** the column physics and the 2-D↔1-D transfer.
- **Kept as oracle:** the Fortran driver, behind a build/runtime switch.

## 2. GPU execution model & constraints

Target AMReX's portability layer (CUDA / HIP / SYCL / OpenMP) via `ParallelFor`
and device lambdas, so one source runs on every backend. A function on the hot
path must be:

1. **`AMREX_GPU_HOST_DEVICE`** (≈ `__host__ __device__`) and self-contained — it
   may call only other device functions.
2. **Allocation-free.** No `new`/`malloc`/`std::vector`/Fortran `allocate` inside
   the kernel. The 1-D column working set (today the nested `noahmp_type`) must
   become **fixed-size stack/local storage** sized by compile-time-or-launch
   constants (`NSOIL`, `NSNOW`, `numrad`, …). This is the single biggest design
   change.
3. **I/O-free, intrinsic-safe.** No `print`, no file ops, no host-only library
   calls in the kernel. Errors cannot `abort` per-thread — use a returned status
   / flagged cell, reduced and reported on the host (`noahmp::fatal`).
4. **Capture-by-value friendly.** Kernels capture `NoahmpArray` views (pointer +
   bounds, trivially copyable) and a parameter struct by value; large read-only
   tables live in device memory.

## 3. Principal challenges (size the work honestly)

| Challenge | Why it's hard | Approach |
|-----------|---------------|----------|
| **`noahmp_type` is a deep nested derived type** (`config/energy/forcing/water/biochem`, each with `state`/`param`/`flux`…) | It's the entire column working set; allocatable members | Mirror as a flat, fixed-size C++ POD per column (stack-local in the kernel) or a struct-of-arrays in device scratch. No allocations. |
| **Dynamic allocation inside physics** | Many routines allocate local arrays sized by `NSOIL`/`NSNOW` | Convert to fixed-size locals (`amrex::GpuArray<noahmp_real, NSOIL_MAX>`); cap layer counts with compile-time maxima. |
| **2-D↔1-D transfer** (`*VarInTransfer`/`*VarOutTransfer`) | 10 Fortran modules of hand copies | Fold the gather/scatter into the kernel: read `NoahmpArray(i,*,j)` into the local column, run, write back. |
| **Table parameters** (`NoahmpReadTableMod`, `NoahmpTable.TBL`) | Read on host in Fortran; needed on device | Read on host (keep Fortran/host reader), pack into a device-resident parameter struct/table, capture by value. |
| **Translating ~Fortran column code to C++** | Volume + index-base (1-based) + column-major semantics | `NoahmpArray` already gives 1-based, column-major indexing; port routine-by-routine, oracle-checked. |
| **Glacier vs. land split + sea-ice/water masking** | Branchy per-cell control flow | Preserve the masks (`XLAND`, `XICE`, `XICE_THRESHOLD`, veg `IndexIcePoint`) as in-kernel branches; divergence is acceptable initially. |
| **Determinism / bitwise parity** | Reductions, transcendental functions, FMA differ host vs. device | Accept tolerance-based parity for GPU; keep the Fortran path for bitwise restart validation (see [`spec-io-restart.md`](spec-io-restart.md)). |

## 4. Phased task list

### Phase 0 — groundwork (no behavior change)
- [ ] Adopt AMReX (or a thin shim) as the portability dependency for the driver;
      decide the kernel macro set (`AMREX_GPU_HOST_DEVICE`, `GpuArray`, launch).
- [ ] Annotate `NoahmpArray{1,2,3}D` access operators `AMREX_GPU_HOST_DEVICE`;
      confirm they remain trivially copyable and that `noahmparray_check` is
      device-safe or compiled out on device. (See
      [`spec-memory-safety.md`](spec-memory-safety.md) §4.)
- [ ] Define compile-time/launch maxima (`NSOIL_MAX`, `NSNOW_MAX`, `NRAD_MAX`).
- [ ] Stand up a **column oracle harness**: capture inputs/outputs of
      `NoahmpMain` for a set of representative cells (land, glacier, sea-ice,
      water) to regression-test the C++ port against.

### Phase 1 — data model
- [ ] Design the flat C++ column state POD mirroring `noahmp_type` (no
      allocatables; fixed-size members). Document the field mapping.
- [ ] Design the device-resident **parameter table** struct; write a host packer
      that fills it from the existing Fortran/host table reader output.
- [ ] Design the per-block kernel signature: which `NoahmpArray` views + scalars
      it captures.

### Phase 2 — port the column solver
- [ ] Port leaf physics routines bottom-up (energy → water → biochem → glacier),
      each as `AMREX_GPU_HOST_DEVICE`, each validated against the oracle within
      tolerance.
- [ ] Assemble `noahmp_column(...)` (the C++ analogue of `NoahmpMain`/
      `NoahmpMainGlacier`) calling the ported leaves on stack-local state.
- [ ] Replace per-thread aborts with a status flag + host-side reduction/report.

### Phase 3 — driver integration
- [ ] Implement the C++ `DriverMain` body: the pre-loop setup currently at the top
      of `NoahmpDriverMain` (WRF layer fill, `RAINBL*DTBL`, `SWDDIR/SWDDIF`,
      `YEARLEN`, `ZSOIL`, soil-timestep accumulators) in C++ on the box.
- [ ] Fold `*VarInTransfer` / `*VarOutTransfer` into the kernel's gather/scatter.
- [ ] Launch `ParallelFor` over `(its:ite, jts:jte)` with the masks; call
      `noahmp_column` per land cell.
- [ ] Add a **dispatch switch**: `NoahmpIO_type::DriverMain()` chooses C++/GPU vs.
      the Fortran `*_fi` path (build flag + runtime env override).

### Phase 4 — validation & performance
- [ ] Whole-step parity: C++/CPU vs. Fortran on a real case (field-by-field,
      tolerance documented).
- [ ] GPU correctness across backends (CUDA/HIP at minimum).
- [ ] Profile; remove host↔device copies from the hot path (keep state resident on
      device across steps; copy to host only for I/O / restart).
- [ ] Restart parity: the Fortran path stays the bitwise oracle; document GPU
      tolerance vs. it.

### Phase 5 — consolidation
- [ ] Decide the long-term status of the Fortran physics path (keep as oracle /
      retire). Update [`spec-overview.md`](spec-overview.md) and
      [`spec-fc-api.md`](spec-fc-api.md) to reflect the dual path.
- [ ] If the column physics no longer crosses into Fortran, prune the now-unused
      `bind(C)` driver shim and transfer modules from the ABI surface (keep I/O
      and init shims, which stay host/Fortran).

## 5. Design principles to hold throughout

1. **The ABI stays valid the whole way.** Each phase is independently
   buildable/testable with the Fortran path as fallback; never a big-bang switch.
2. **Precision discipline unchanged.** `noahmp_real` everywhere on device too; no
   literal `double`/`float`. (See [`spec-fc-api.md`](spec-fc-api.md) §2.)
3. **No allocation, no I/O, no MPI on the device.** Errors flag-and-reduce, then
   `noahmp::fatal` on the host (see [`spec-memory-safety.md`](spec-memory-safety.md) §5).
4. **State resident on device.** Forcing in / fluxes out and I/O are the only
   host↔device transfers; the column working set never leaves the kernel.
5. **Oracle-checked, routine by routine.** Port small, diff against Fortran, keep
   the diff harness in CI.

## 6. Open questions

- Generate the flat column POD and the 2-D↔1-D gather/scatter from a `@NoahmpMacro`
  region (extending [`spec-fc-api.md`](spec-fc-api.md)'s generator), to avoid
  hand-maintaining a second large field list?
- Per-cell stack pressure on GPU: does the full column state fit in registers/
  local memory, or is per-block device scratch (struct-of-arrays) needed?
- Variable layer counts across cells/levels vs. compile-time maxima — pad to max,
  or specialize?
- How much physics divergence (glacier/sea-ice/water branches) is acceptable
  before it's worth sorting cells by type?
