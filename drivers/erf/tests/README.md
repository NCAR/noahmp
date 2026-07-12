# Noah-MP ERF driver — regression tests

Fast, dependency-free tests that link the built `noahmp` library and exercise
the C++ ↔ Fortran coupling boundary through **both** interfaces. Plain CTest with
a tiny assert helper (`test_util.H`) — no GoogleTest, matching the repo's
minimal-dependency style. See `../dev/spec-overview.md` and `../dev/spec-fc-api.md`
for the machinery under test.

## Running

Tests are built when `NOAHMP_ENABLE_TESTS` is ON — the default for a **standalone**
`cmake .../drivers/erf` build, and OFF when ERF consumes this directory as a
subproject (so ERF's build is never perturbed).

```sh
# under the site environment (compilers + spack netcdf-fortran):
source <project>/sites/<site>/environment.sh
export NETCDF_HOME=<...>/netcdf-install-<site>

cmake -DNETCDF_DIR=$NETCDF_HOME -DNOAHMP_ENABLE_TESTS=ON <...>/drivers/erf
make -j
ctest --output-on-failure
```

The `tests/NOAH-MP` harness (`build.sh` / `job.submit`) does exactly this.

## Coverage

**Tier 1 — coupling machinery (no external data files)**

| Test | Interface | What it guards |
|------|-----------|----------------|
| `array_view` | C++ | `NoahmpArray{1,2,3}D` Fortran-column-major index math with non-1 lower bounds |
| `abi_guard` | C++/Fortran | precision guard (`NoahmpRealSize_fi` == `sizeof(noahmp_real)`), flat-pointer layout, `NoahmpIO_AssertAbi` |
| `fatal_handler` | C++ | `NoahmpIO_fatal` host-handler dispatch + `std::abort` fallback (forked) |
| `lifecycle` | C++ | `NoahmpIO_vector`/`NoahmpIO_type` cold-init order; array-view bounds derived from C++-owned scalars (proves scalar wiring) |
| `sharing` | C++ + Fortran | **zero-copy proof**: C++ view writes seen by the Fortran module-global block and vice versa (2-D + 3-D column-major) |
| `fortran_alloc` | Fortran | `NoahmpIOVarInitDefault` array shapes/bounds; `NoahmpRealSize_fi`; idempotent `NoahmpIOTypeVectInit_fi` |
| `water_flux_units` | Fortran | `WaterVarOutTransfer` units/accumulation contract: `QTDRAIN`/`SFCRUNOFF`/`UDRUNOFF` accumulate a depth `[mm per soil timestep]`, `RUNSFXY`/`RUNSBXY` snapshot, glacier sentinel zeroes `TileDrain` and scales runoff by `MainTimeStep` |
| `cpp_guard_*` | C++ | `NoahmpIO_vector` size-once / non-zero guards (WILL_FAIL) |
| `fortran_guard_*` | Fortran | `NoahmpIOTypeVectInit_fi` level / NBlocks / re-init validation (WILL_FAIL) |

**Tier 2 — config smoke (shipped `NoahmpTable.TBL` + tiny `namelist.erf`, no NetCDF)**

| Test | Interface | What it guards |
|------|-----------|----------------|
| `config_smoke` | C++ + Fortran | `ReadNamelist` (asserts parsed coupling scalars) and `ReadTable` (asserts known general-parameter table values) via the public C++ API |

**Tier 3 — NetCDF I/O + driver (parallel NetCDF-4; each test creates its own MPI world)**

These drive the real NetCDF modules and the physics driver end to end. To stay
committable they generate their own fixtures instead of depending on a large
land file: `make_wrfinput` (in `test_io_support.F90`) synthesizes a tiny
wrfinput/WPS file carrying every attribute/variable the reader needs, and the
write/restart tests round-trip their own output. The read and driver tests also
accept **any real local WPS/wrfinput file** via `NOAHMP_TEST_WRFINPUT=<path>`
(standard `west_east`/`south_north`/`soil_layers_stag` dims), skipping the
synthetic fixture.

| Test | Interface | What it guards |
|------|-----------|----------------|
| `io_readland` | Fortran | `ReadLandHeader`/`ReadLandMain`: header globals (`xsglobal`/offsets), and exact recovery of per-cell fields (`XLAT`, `TERRAIN`, `TMN`, `TSK`, `IVGTYP`, `ISLTYP`, soil `TSLB`/`SMOIS`) written to the fixture |
| `io_writeland` | Fortran | `NoahmpWriteLand`: reopens the produced `lnd*/Level_0.nc` and checks dims + written `TERRAIN`/`TSK`/`HFX`/`TSLB` values |
| `io_restart_roundtrip` | Fortran | `WriteRestart`→`ReadRestart` **bit-exact** round trip over the full prognostic state, plus the checkpoint precision (`NF90_DOUBLE`/`REAL`) and `NSOIL`/`NSNOW`/`ISNOWXY` metadata contract |
| `io_restart_mismatch` | Fortran | restart layer-geometry guard aborts on an `NSOIL` mismatch (WILL_FAIL) |
| `io_driver` | Fortran | full cold-init chain (`ReadNamelist`→`ReadLandHeader`→`VarInitDefault`→`ReadTable`→`ReadLandMain`→`InitMain`) then `NoahmpDriverMain` over two steps; asserts finite, physical surface state |
| `io_driver_cpp` | C++ + Fortran | same full step through the **public C++ API in ERF's exact call order** (`ERF_NOAHMP_Init.cpp` / `ERF_NOAHMP_Advance.cpp`), incl. `WriteLand(0)` |
| `io_abort_check_ok` | Fortran | `NoahmpFatalMod::check_nc` returns on `NF90_NOERR` |
| `io_abort_check_bad` | Fortran | `check_nc` aborts on a NetCDF error status (WILL_FAIL) |
| `io_abort_direct` | Fortran | `NoahmpIO_abort` terminates (WILL_FAIL) |

Tier 3 needs a **parallel-I/O-capable** NetCDF (the write/restart paths use
`NF90_MPIIO`); point `-DNETCDF_DIR` at the project's parallel `netcdf-fortran`
install (`nc-config --has-parallel4` = yes), as `tests/NOAH-MP/build.sh` does.

## Files

- `test_util.H` — `CHECK` / `CHECK_EQ` / `CHECK_CLOSE` / `TEST_SUMMARY`.
- `test_support.F90` — `bind(C)` helpers compiled **into the tests** (not the
  library): set safe `IOPT_*` (stand-in for `ReadNamelist`), read/write
  Fortran-owned arrays, probe table scalars.
- `test_abort_handler.cpp` — routes the fatal `std::abort` to a non-zero
  `_Exit` in the WILL_FAIL executables so CTest inverts the result reliably.
- `namelist.erf` — minimal valid `NOAHLSM_OFFLINE` fixture for `config_smoke`
  and the Tier-3 driver tests.
- `test_io_support.F90` — Tier-3 shared helpers: synthetic wrfinput generator
  (`make_wrfinput`), single-block setup, MPI-world init, deterministic reference
  fields, assert harness, and `bind(C)` shims used by `io_driver_cpp`.

## Notes

- **Precision**: tests use `noahmp_real` / `c_kind_noahmp`, never literal
  `double`/`C_DOUBLE`, and inherit the library's precision. For a double-precision
  build pass the same flag to **both** compilers (as ERF/the Makefile do), e.g.
  `-DCMAKE_Fortran_FLAGS=-DDOUBLE_PREC -DCMAKE_CXX_FLAGS=-DDOUBLE_PREC`.
- Coupled dimension scalars (`xstart`, `nsoil`, `itimestep`, …) are C++-owned
  **pointer** components on the Fortran side, normally wired by
  `ScalarInitDefault`. The pure-Fortran `fortran_alloc` test therefore
  `allocate`s the ones `NoahmpIOVarInitDefault` touches before use.
