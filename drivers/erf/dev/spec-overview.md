# Spec: ERF ↔ Noah-MP driver — architecture overview

> Status: living document · Audience: contributors to `drivers/erf` · Read this
> before the other specs.

## 1. Purpose

ERF is a C++/AMReX atmospheric model. Noah-MP is a Fortran land-surface model.
This driver is the boundary that lets ERF own the parallel runtime (MPI domain
decomposition, AMR, time stepping) while Noah-MP runs as a *column physics
library* underneath it. The driver's job is to:

1. Present Noah-MP to ERF as a **C++ class** (`NoahmpIO_type`) with ordinary
   methods (`InitMain()`, `DriverMain()`, `WriteLand()`, …) — no Fortran in
   sight from the caller's side.
2. Share the large state arrays **without copying** them across the language
   boundary every step: ERF/Noah-MP work on the same memory through pointers.
3. Carry **zero parallel-runtime dependency inside Noah-MP** — Noah-MP never
   calls MPI directly; the host owns MPI and AMReX.

## 2. The two sides of the boundary

```
  ┌──────────────── C++ (ERF host) ────────────────┐   ┌──────── Fortran (Noah-MP) ────────┐
  │  NoahmpIO_type            (public API class)    │   │  NoahmpIO_type   (storage type)    │
  │    └ fptr : NoahmpIO_type_fi  (flat pointers)   │◄─►│  NoahmpIO_type_fi  bind(C) mirror  │
  │  NoahmpIO_vector  (per-AMR-level block array)   │   │  NoahmpIO_vect(level)%NoahmpIO(:)  │
  │  NoahmpArray{1,2,3}D<noahmp_real>  (array view) │   │  allocatable real arrays (owner)   │
  └─────────────────────────────────────────────────┘   └────────────────────────────────────┘
                          ▲                                              ▲
                 NoahmpArray.H / NoahmpIO.H               NoahmpIOVarType.F90 / *Mod.F90
```

- **C++ `NoahmpIO_type`** is the API the host calls. It holds a private member
  `fptr` of type `NoahmpIO_type_fi` — a flat struct of opaque pointers that *is*
  the ABI.
- **Fortran `NoahmpIO_type`** is the storage type: it owns the `allocatable`
  arrays where the physics state actually lives.
- **`NoahmpIO_type_fi`** exists on *both* sides (a C++ struct and a Fortran
  `bind(c)` derived type) with identical member count/order/precision. It is the
  wire format: every member is a pointer handle.
- Scalars are owned by C++ (Fortran points back at them); arrays are owned by
  Fortran (C++ `NoahmpArray` views point into them via `C_LOC`).

See [`spec-fc-api.md`](spec-fc-api.md) for the exact wiring and
[`spec-memory-safety.md`](spec-memory-safety.md) for why this ownership split is
safe.

## 3. AMR levels and blocks

ERF is block-structured AMR. The driver mirrors that shape:

- **Level** — an AMR refinement level, `NLEVEL_MIN..NLEVEL_MAX` (`0..2` today,
  in `NoahmpIO_fi.F90-mc`). Each level has its own block array.
- **Block (`blkid`)** — one AMReX box / MPI-decomposed tile within a level,
  `0..NBlocks-1`. Each block is one `NoahmpIO_type` instance and runs Noah-MP
  over its own `(xstart:xend, ystart:yend)` patch.

On the C++ side this is `NoahmpIO_vector` (one per level); on the Fortran side
`NoahmpIO_vect(level)%NoahmpIO(0:NBlocks-1)`. The array is **sized exactly once**
(`NoahmpIO_vector::resize(size, level)` → `NoahmpIOTypeVectInit_fi`) and never
resized again — element relocation would dangle the self-referential pointers.

Global-vs-local grid geometry per block:

- `xstart..xend`, `ystart..yend` — this block's owned patch (1-based, global
  indices).
- `xsglobal`, `ysglobal` — global domain extents (the NetCDF dataset size).
- `xoffset`, `yoffset` — origin offset used to turn global indices into a NetCDF
  hyperslab `start` (see the I/O specs).

## 4. Lifecycle of a run

Per level, the host:

1. `NoahmpIO_vector::resize(NBlocks, level)` — size the block array once.
2. Per block, in order, call the methods that map 1:1 to `bind(C)` entry points:
   - `ScalarInitDefault()` — wire each block's scalar pointers to C++ memory
     (must run first; later calls assert it ran).
   - `ReadNamelist()`, `ReadTable()` — configuration.
   - `ReadLandHeader()`, `ReadLandMain()` — static land/soil inputs.
   - `VarInitDefault()` + `InitMain()` — allocate arrays, set initial state.
   - (restart only) `ReadRestart(dir)` — overwrite cold-init state with a
     checkpoint (see [`spec-io-restart.md`](spec-io-restart.md)).
3. Per time step, per block: fill forcing arrays, call `DriverMain()`, read back
   fluxes. `DriverMain` loops `j=jts..jte`, `i=its..ite`, packs each land column
   into the 1-D `noahmp_type`, runs `NoahmpMain` (or `NoahmpMainGlacier`), and
   unpacks the result (`NoahmpDriverMainMod.F90`).
4. Output: `WriteLand(filenum)` each output step, `WriteRestart(dir)` at
   checkpoints (see [`spec-io-parallel.md`](spec-io-parallel.md)).

## 5. Design tenets (do not regress)

1. **One source of truth for the ABI.** The member list in the
   `@NoahmpMacro:Source` block of `NoahmpIO.H-mc` generates *both* sides. Count,
   order, and precision match by construction. See [`spec-fc-api.md`](spec-fc-api.md).
2. **No parallel-runtime dependency in Noah-MP.** Noah-MP never calls MPI/AMReX.
   Fatal errors route through `NoahmpIO_fatal()` → a host-installed handler
   (`amrex::Abort` → `MPI_Abort`). See [`spec-memory-safety.md`](spec-memory-safety.md).
3. **Precision is host-owned.** `DOUBLE_PREC` selects `noahmp_real` /
   `c_kind_noahmp` so the boundary always matches `amrex::Real`. A run-time guard
   aborts on a mismatch between the two separately-invoked compilers.
4. **No silent relocation.** Both `NoahmpIO_type` and `NoahmpIO_vector` forbid
   the std-container operations that would move elements and dangle the
   self-referential pointers — compile errors, not run-time surprises.
5. **Collective, parallel I/O.** All NetCDF output is written collectively over
   `NoahmpIO%comm` (the host's MPI communicator) as parallel NetCDF-4.

## 6. File map

Generated/templated coupling glue (see [`spec-fc-api.md`](spec-fc-api.md)):

| Template (`-mc`, tracked) | Generated target (gitignored) | Role |
|---------------------------|-------------------------------|------|
| `NoahmpIO.H-mc` | `NoahmpIO.H` | C++ API class + `fi` struct + ABI guards; **holds the Source block** |
| `NoahmpIO.cpp-mc` | `NoahmpIO.cpp` | C++ method bodies, fatal handler, ABI precision guard |
| `NoahmpIO_fi.F90-mc` | `NoahmpIO_fi.F90` | Fortran `bind(C)` entry points, block resolution |
| `NoahmpIOVarType.F90-mc` | `NoahmpIOVarType.F90` | Fortran storage type (`NoahmpIO_type`) |
| `NoahmpIOVarInitMod.F90-mc` | `NoahmpIOVarInitMod.F90` | Fortran array `allocate()` |

Hand-written (not generated):

The public C++ client API is exposed at global scope (`NoahmpIO_type`,
`NoahmpIO_vector`, `NoahmpArray{2,3}D`, `noahmp_real`,
`NoahmpIO_fatal` / `NoahmpIO_set_fatal_handler`); see [`spec-fc-api.md`](spec-fc-api.md) §3.

| File | Role |
|------|------|
| `NoahmpArray.H` | `NoahmpArray{1,2,3}D` Fortran-layout array views + `noahmp_real` |
| `NoahmpFatal.H` | `NoahmpIO_fatal` / `NoahmpIO_set_fatal_handler` abstraction |
| `NoahmpDriverMainMod.F90` | Per-step driver: 2-D ↔ 1-D transfer + `NoahmpMain` |
| `NoahmpWriteLandMod.F90` | Parallel land-output writer ([`spec-io-parallel.md`](spec-io-parallel.md)) |
| `NoahmpWrite/ReadRestartMod.F90` | Checkpoint/restart ([`spec-io-restart.md`](spec-io-restart.md)) |
| `Noahmp{Init,ReadNamelist,ReadTable,ReadLand}*Mod.F90` | Init / config / static input |
| `tools/NoahmpMacro.py` | The code generator |

## 7. Glossary

This is the shared glossary for all the specs; the others link here on first use
rather than re-defining terms.

- **ABI** — the binary contract between the C++ `NoahmpIO_type_fi` and the
  Fortran `bind(c)` mirror: same member count, order, all-pointer types. ("API"
  is what the *programmer* sees — method names, arguments; "ABI" is what the
  *compiled binaries* must agree on — byte layout. The boundary needs both.)
- **`fi`** — "Fortran interop"; the boundary mirror struct and `*_fi` shims.
- **handle / opaque pointer** — a bare memory address with no type information
  attached; a member of `NoahmpIO_type_fi` is just such an address.
- **self-referential struct** — an object that stores the addresses of its *own*
  members, so a second language can read/write them. Powerful, but it means the
  object must never be moved or copied (see
  [`spec-memory-safety.md`](spec-memory-safety.md)).
- **standard-layout / "no padding"** — the struct is a plain flat array of
  pointers with no hidden bytes inserted by the compiler, so the C++ and Fortran
  pictures of it line up exactly.
- **hyperslab** — a NetCDF term: a rectangular sub-region of a larger array, given
  as an origin (`start`) and a size (`count`). Each block writes its own
  hyperslab of the shared file (see [`spec-io-parallel.md`](spec-io-parallel.md)).
- **collective I/O** — every MPI rank participates in one shared write to one file,
  instead of writing private per-rank files that need stitching afterward.
- **oracle** — a trusted reference implementation you check a new one against; here,
  the CPU Fortran path is the oracle for the GPU path
  (see [`plan-cpp-interface.md`](plan-cpp-interface.md)).
- **column / 1-D physics** — Noah-MP's per-grid-cell solver (`noahmp_type`,
  `NoahmpMain`), the part destined for the GPU (see
  [`plan-cpp-interface.md`](plan-cpp-interface.md)).
