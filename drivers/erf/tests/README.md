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
| `cpp_guard_*` | C++ | `NoahmpIO_vector` size-once / non-zero guards (WILL_FAIL) |
| `fortran_guard_*` | Fortran | `NoahmpIOTypeVectInit_fi` level / NBlocks / re-init validation (WILL_FAIL) |

**Tier 2 — config smoke (shipped `NoahmpTable.TBL` + tiny `namelist.erf`, no NetCDF)**

| Test | Interface | What it guards |
|------|-----------|----------------|
| `config_smoke` | C++ + Fortran | `ReadNamelist` (asserts parsed coupling scalars) and `ReadTable` (asserts known general-parameter table values) via the public C++ API |

Deferred to a future integration tier (need domain-sized NetCDF fixtures + MPI
launch): `ReadLand*`, `InitMain`, `DriverMain`, `WriteLand`, restart round-trip.

## Files

- `test_util.H` — `CHECK` / `CHECK_EQ` / `CHECK_CLOSE` / `TEST_SUMMARY`.
- `test_support.F90` — `bind(C)` helpers compiled **into the tests** (not the
  library): set safe `IOPT_*` (stand-in for `ReadNamelist`), read/write
  Fortran-owned arrays, probe table scalars.
- `test_abort_handler.cpp` — routes the fatal `std::abort` to a non-zero
  `_Exit` in the WILL_FAIL executables so CTest inverts the result reliably.
- `namelist.erf` — minimal valid `NOAHLSM_OFFLINE` fixture for `config_smoke`.

## Notes

- **Precision**: tests use `noahmp_real` / `c_kind_noahmp`, never literal
  `double`/`C_DOUBLE`, and inherit the library's precision. For a double-precision
  build pass the same flag to **both** compilers (as ERF/the Makefile do), e.g.
  `-DCMAKE_Fortran_FLAGS=-DDOUBLE_PREC -DCMAKE_CXX_FLAGS=-DDOUBLE_PREC`.
- Coupled dimension scalars (`xstart`, `nsoil`, `itimestep`, …) are C++-owned
  **pointer** components on the Fortran side, normally wired by
  `ScalarInitDefault`. The pure-Fortran `fortran_alloc` test therefore
  `allocate`s the ones `NoahmpIOVarInitDefault` touches before use.
