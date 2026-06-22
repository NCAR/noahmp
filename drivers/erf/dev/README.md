# Noah-MP ERF driver — developer specs (`dev/`)

This directory holds the design specifications and forward-looking plans for the
**ERF driver of Noah-MP** (`drivers/erf`): the C++ ↔ Fortran coupling layer that
lets [ERF](https://github.com/erf-model/ERF) (an AMReX-based atmospheric model)
drive the Noah-MP land-surface model.

These are *living design documents*, not user manuals. They describe how the
boundary is built, why it is built that way, and where it is going. Public,
user-facing documentation lives at
<https://erf.readthedocs.io/en/latest/CouplingToNoahMP.html>.

## Index

| File | What it covers |
|------|----------------|
| [`spec-overview.md`](spec-overview.md) | Architecture, lifecycle, design tenets, file map. **Start here.** |
| [`spec-fc-api.md`](spec-fc-api.md) | The Fortran ↔ C ABI/API and the `@NoahmpMacro` code generator. |
| [`spec-memory-safety.md`](spec-memory-safety.md) | Self-referential pointer handling, copy/move bans, bounds checking, fatal-handler abstraction. |
| [`spec-io-parallel.md`](spec-io-parallel.md) | Parallel (collective MPI) NetCDF land-output implementation. |
| [`spec-io-restart.md`](spec-io-restart.md) | Full prognostic-state checkpoint/restart I/O. |
| [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md) | Workflow: expose a new variable across the boundary. |
| [`plan-cpp-interface.md`](plan-cpp-interface.md) | Plan/task list to port the driver to C++ and run Noah-MP physics on GPU. |

## Conventions used across these docs

- **Templates vs targets.** Files ending in `-mc` are hand-edited *templates*.
  The build runs `tools/NoahmpMacro.py` to expand them into *generated targets*
  (`NoahmpIO.H`, `NoahmpIO.cpp`, `NoahmpIO_fi.F90`, `NoahmpIOVarType.F90`,
  `NoahmpIOVarInitMod.F90`), which are gitignored. **Never hand-edit a target.**
- **`noahmp_real` / `c_kind_noahmp`.** The coupling precision, selected by the
  `DOUBLE_PREC` macro so it tracks `amrex::Real`. Never hardcode `double` /
  `C_DOUBLE`.
- **`fi` suffix.** Marks the Fortran-interop boundary mirror (`NoahmpIO_type_fi`
  and the `*_fi` `bind(C)` entry points).
