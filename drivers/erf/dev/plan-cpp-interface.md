# Plan: GPU-enable Noah-MP under the ERF driver via Fortran offload

> Status: **proposal / roadmap** (not yet started) · Audience: contributors
> planning the GPU effort. · Prereqs read: [`spec-overview.md`](spec-overview.md),
> [`spec-fc-api.md`](spec-fc-api.md), [`spec-memory-safety.md`](spec-memory-safety.md).
>
> Filename kept (`plan-cpp-interface.md`) for link stability; the approach is
> **Fortran GPU offload**, not a C++ rewrite of the physics. See §0.

## 0. TL;DR

The **boundary** is already C++ (`NoahmpIO_type`, `NoahmpArray`,
`noahmp::fatal`), and it **stays** that way: `drivers/erf` remains the C++ host
ERF calls. The **physics in `src/` stays Fortran** — it is *scalar, per-column*
math, and we do **not** rewrite it in C++.

What changes is the *compilation target* of that physics. To run Noah-MP on the
GPU, the per-column solver (`NoahmpMain` / `NoahmpMainGlacier`) and the `i,j`
driver loop that drives it (`NoahmpDriverMainMod.F90`) are compiled as **GPU
device code using a Fortran offload model** — OpenACC as the primary, OpenMP
target as the portability hedge. The `i,j` loop becomes an offload region; each
iteration calls `NoahmpMain` marked as a device routine. The host C++
`DriverMain()` still just calls (via `bind(C)`) the Fortran driver entry point;
the loop runs on the device underneath it.

> **Why not a C++ `ParallelFor`?** An AMReX device lambda is GPU code and can only
> call functions *compiled for the device*. Host-compiled Fortran cannot be called
> from a kernel — not a complexity issue, a compilation-target one. So to keep the
> physics in Fortran *and* run on GPU, the **Fortran itself** must be compiled for
> the device and the loop offloaded **in Fortran**. That is option B.

The end state: `NoahmpIO_type::DriverMain()` → `bind(C)` → a Fortran
`NoahmpDriverMain` whose `i,j` loop is an `!$acc parallel loop` over the block's
box, calling a device-compiled `NoahmpMain` per land cell, with state resident on
the device. The unmodified-offload (host) Fortran build stays the reference oracle.

> **Deliberate trade-off (read this):** Fortran offload is **NVHPC/NVIDIA-first**.
> AMD (OpenMP target via the AMD/Cray Fortran compilers) is possible but less
> mature; Intel GPUs (SYCL) have **no** Fortran offload path. This is the price of
> keeping the physics in Fortran instead of porting to C++ — we lose AMReX's "one
> source, every backend" portability. The host CPU Fortran path always builds with
> any compiler.

## 1. Where we are vs. where we're going

```
  TODAY                                            TARGET
  ─────                                            ──────
  ERF (C++)                                        ERF (C++)
    └ NoahmpIO_type::DriverMain()                    └ NoahmpIO_type::DriverMain()
        └ NoahmpDriverMain_fi  (bind C)                  └ NoahmpDriverMain_fi  (bind C)
            └ NoahmpDriverMain  (Fortran)                    └ NoahmpDriverMain  (Fortran)
                └ i,j loop (host, serial/threaded)               └ !$acc parallel loop (i,j)  ← on GPU
                    └ *VarInTransfer (2D→1D)                          └ private column state (fixed-size)
                    └ NoahmpMain    (Fortran column, host)            └ NoahmpMain  (!$acc routine seq, on device)
                    └ *VarOutTransfer (1D→2D)                     (state resident on device; host copy for I/O only)
```

- **Reused as-is:** the AMR level/block model, the precision discipline
  (`noahmp_real` / `kind_noahmp`), the fatal-handler abstraction (host-side), the
  parallel NetCDF I/O (stays host-side, runs after a device→host copy), and **the
  entire `src/` physics source** — *annotated*, not rewritten.
- **Changed:** `src/` routines gain device annotations and lose dynamic
  allocation; the `noahmp_type` column working set becomes device-friendly
  (private per thread, fixed-size, no allocatable components); the `i,j` loop
  becomes an offload region; data-residency clauses are added.
- **Kept as oracle:** the host Fortran path (compile without offload flags).

> **Note on `NoahmpArray`:** unlike a C++ port, option B does **not** put the C++
> `NoahmpArray` views on the device hot path. The offloaded Fortran reads the
> Fortran `allocatable` arrays (the `NoahmpIO_type` storage) directly. So the
> `AMREX_GPU_HOST_DEVICE` annotation of `NoahmpArray` is **not** needed here — that
> task belonged to the C++-port approach and is dropped.

## 2. GPU execution model & constraints

Target a **Fortran offload model**. Pick one primary; keep a second compiling:

- **OpenACC** (`!$acc`) — *recommended primary*. Most mature in NVHPC; the
  de-facto weather/climate Fortran-on-GPU path. Best control over routine marking
  and data residency.
- **OpenMP target** (`!$omp target` / `declare target`) — *portability hedge*.
  Standard, and the realistic route to AMD GPUs (AMD/Cray Fortran). Historically
  less mature in `nvfortran` but improving.
- **`do concurrent` + `-stdpar=gpu`** — cleanest source (no pragmas) but the least
  control (mandates managed memory, limited routine-marking story). Noted as an
  option, not the plan of record.

A routine on the hot path must be:

1. **Device-marked.** *Every* routine transitively reachable from the loop needs
   `!$acc routine seq` (or `!$omp declare target`). Noah-MP's call tree is deep
   (energy → water → biochem → glacier), so this is a **mechanical sweep over all
   of `src/`** plus the `drivers/erf` transfer modules.
2. **Allocation-free on device.** No `allocate`/`deallocate` inside the kernel. The
   in-physics locals sized by `NSOIL`/`NSNOW` become **fixed-size automatic
   arrays** capped by compile-time maxima.
3. **I/O-free, intrinsic-safe.** No `print`, no file ops, no MPI in the kernel.
   Errors cannot `abort` per-thread — use a returned status / flagged cell,
   reduced and reported on the host (`noahmp::fatal`).
4. **Operating on per-thread-private column state.** The `noahmp_type` working set
   must be `private`/`firstprivate` per loop iteration and **free of allocatable
   components on device** (deep-copy of allocatable-component derived types to the
   device is the classic offload trap).

Plus a **data-residency** contract: the large `NoahmpIO` state arrays must live on
the device *across* timesteps (managed memory, or explicit `enter data` /
`target enter data`), copied to host only for I/O / restart.

## 3. Principal challenges (size the work honestly)

| Challenge | Why it's hard | Approach |
|-----------|---------------|----------|
| **`noahmp_type` has allocatable components** (`config/energy/forcing/water/biochem`, each with `state`/`param`/`flux`…) | OpenACC/OMP deep-copy of allocatable-component derived types to device is painful; a per-thread `private` copy must not device-allocate | Make the per-column working set a **fixed-size** derived type (no allocatables); declare it `private`/`firstprivate` in the loop. Same fixed-size redesign a C++ POD would need — but it stays Fortran. |
| **Dynamic allocation inside physics** | `allocate` in a `seq` device routine is slow/unsupported | Convert local allocatables to fixed-size automatic arrays sized by `NSOIL_MAX`/`NSNOW_MAX`/`NRAD_MAX`. |
| **Device-marking the whole call tree** | Every transitively-called routine needs `!$acc routine`/`declare target` | Mechanical sweep across `src/` (and the transfer modules); a CI check that no unmarked routine is reached. Consider generating the markers (see §6). |
| **AMReX ↔ Fortran-offload memory interop** | AMReX owns/manages the arrays; the Fortran offload runtime must see the *same* device memory | Use **one vendor runtime** across both sides (NVHPC: `nvc++`/`nvcc` for AMReX + `nvfortran` for the offload) and a shared memory model (CUDA managed memory, or hand specific device pointers across via `c_loc`). Document the contract. |
| **Table parameters** (`NoahmpReadTableMod`, `NoahmpTable.TBL`) | Read on host; needed on device | Keep the host Fortran reader; `enter data` the tables to the device once; mark read-only. |
| **Glacier vs. land split + sea-ice/water masking** | Branchy per-cell control flow | Keep the masks (`XLAND`, `XICE`, `XICE_THRESHOLD`, `IndexIcePoint`) as in-kernel branches; divergence acceptable initially. |
| **Determinism / bitwise parity** | FMA / transcendentals differ host vs. device | Tolerance-based parity on GPU; keep the host Fortran path as the bitwise oracle (see [`spec-io-restart.md`](spec-io-restart.md)). |
| **Toolchain lock-in** | Fortran offload is NVHPC-first | Accept NVIDIA-first; keep the OpenMP-target variant compiling for AMD; the CPU path builds with any compiler. |

## 4. Phased task list

### Phase 0 — groundwork (no behavior change)
- [ ] Choose the offload model (**recommend OpenACC** primary; OpenMP-target hedge)
      and pin the toolchain (NVHPC `nvfortran` + `nvc++`/`nvcc` for AMReX). Decide
      the memory model (managed vs. explicit data regions).
- [ ] Get `src/` **and** `drivers/erf` compiling under the chosen compiler with
      **offload OFF** (host build) — flush out plain compiler-compatibility issues
      before any GPU work.
- [ ] Define compile-time maxima (`NSOIL_MAX`, `NSNOW_MAX`, `NRAD_MAX`).
- [ ] Stand up a **column oracle harness**: capture inputs/outputs of `NoahmpMain`
      for representative cells (land, glacier, sea-ice, water) to regression-test
      every device-readiness refactor against.

### Phase 1 — make the column device-ready (still host, still CPU)
- [ ] Redesign the `noahmp_type` per-column working set for the device: remove
      allocatable components; fixed-size members. **CPU semantics identical**,
      oracle-checked.
- [ ] Convert in-physics local allocatables to fixed-size automatic arrays.
- [ ] Remove/guard `print`/I/O in the physics; route errors to a status flag
      instead of aborting (see [`spec-memory-safety.md`](spec-memory-safety.md) §5).
- [ ] Validate the whole refactor on CPU against the oracle — *before* any offload.

### Phase 2 — device-enable the physics
- [ ] Sweep `src/`: add `!$acc routine seq` (and/or `!$omp declare target`) to
      every routine in the `NoahmpMain` / `NoahmpMainGlacier` call tree.
- [ ] Device-mark the parameter tables and module constants.
- [ ] Compile the call tree *for the device*; resolve unsupported constructs
      (residual allocations, I/O, non-device intrinsics).

### Phase 3 — offload the driver loop
- [ ] Turn the `i,j` loop in `NoahmpDriverMain` into an offload region
      (`!$acc parallel loop` / `collapse`), with the column state `private` per
      iteration.
- [ ] Fold the masks (`XICE`/`XLAND`/`IndexIcePoint`) into in-kernel branches.
- [ ] Device-mark the `*VarInTransfer` / `*VarOutTransfer` modules (they're scalar
      gather/scatter copies) so they run inside the kernel on device data.
- [ ] Add data residency: `enter data` the `NoahmpIO` state + tables once; keep
      resident across steps; refresh the host only for I/O / restart.
- [ ] Add a **dispatch switch**: `NoahmpIO_type::DriverMain()` chooses the offload
      (GPU) path vs. the host Fortran path (build flag + runtime env override).

### Phase 4 — validation & performance
- [ ] Whole-step parity: GPU vs. host Fortran on a real case (field-by-field,
      tolerance documented).
- [ ] Verify data residency: no per-step host↔device copies except I/O; profile.
- [ ] Validate AMReX ↔ Fortran-runtime memory interop (no double-copies, no
      host/device aliasing bugs).
- [ ] Restart parity: the host Fortran path stays the bitwise oracle; document the
      GPU tolerance vs. it (see [`spec-io-restart.md`](spec-io-restart.md)).

### Phase 5 — consolidation & portability
- [ ] Bring up the OpenMP-target variant for AMD (if required); document the
      supported compilers/GPUs and the memory-model contract.
- [ ] Decide the long-term status of the host-only path (keep as oracle / retire).
- [ ] Update [`spec-overview.md`](spec-overview.md) and
      [`spec-fc-api.md`](spec-fc-api.md) to document the dual path (host Fortran /
      device Fortran) and the toolchain + memory-model contract.

## 5. Design principles to hold throughout

1. **The ABI stays valid the whole way.** Each phase is independently
   buildable/testable with the host Fortran path as fallback; never a big-bang
   switch.
2. **The physics stays Fortran — annotated, not rewritten.** No C++ translation of
   `src/`. The diff against upstream Noah-MP stays a set of pragmas + fixed-size
   refactors, not a reimplementation.
3. **Precision discipline unchanged.** `noahmp_real` / `kind_noahmp` everywhere on
   device too. (See [`spec-fc-api.md`](spec-fc-api.md) §2.)
4. **No allocation, no I/O, no MPI on the device.** Errors flag-and-reduce, then
   `noahmp::fatal` on the host (see [`spec-memory-safety.md`](spec-memory-safety.md) §5).
5. **State resident on device.** Forcing in / fluxes out and I/O are the only
   host↔device transfers; the column working set never leaves the kernel.
6. **One vendor runtime across the boundary.** AMReX (C++) and Noah-MP (Fortran)
   share the device memory model — same toolchain family, same managed-memory /
   device-pointer contract.
7. **Oracle-checked.** Every device-readiness refactor (Phase 1) is validated on
   CPU first; the GPU path is then checked within tolerance.

## 6. Open questions

- **OpenACC vs. OpenMP target** as primary? OpenACC = maturity on NVIDIA; OpenMP
  target = the path to AMD. Is `do concurrent` + `-stdpar=gpu` (cleanest source,
  but managed-memory-mandatory and weak routine-marking control) viable instead?
- **Managed (unified) memory vs. explicit data clauses** for the `NoahmpIO`
  arrays — which gives acceptable performance without hand-maintaining `map`/`data`
  lists as the coupled-variable set grows (see
  [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md))?
- Can AMReX's arena allocator and the Fortran offload runtime **share device
  allocations** cleanly, or do we hand specific arrays across via `c_loc` / device
  pointers?
- **Stack/register pressure** of the per-thread fixed-size column state: does it
  fit in device local memory, or do we need block-shared scratch?
- Generate the **device annotations** and the **fixed-size column type** from a
  `@NoahmpMacro` region (extending [`spec-fc-api.md`](spec-fc-api.md)'s generator),
  to avoid hand-maintaining the sweep and a second large field list?
- How much physics divergence (glacier/sea-ice/water branches) is acceptable
  before it's worth sorting cells by type?
