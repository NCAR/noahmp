# ===========================================================================
# Regression tests for the Noah-MP ERF driver.
#
# Added via `add_subdirectory(tests)` from the parent when NOAHMP_ENABLE_TESTS
# is ON, so the `noahmp` target and its generated Fortran .mod files (in
# ${CMAKE_BINARY_DIR}) already exist. Tests link the library and exercise the
# C++ <-> Fortran coupling boundary through both interfaces. Run with `ctest`.
# ===========================================================================

# Where the library's generated Fortran modules (noahmpio_fi.mod, ...) live.
set(NOAHMP_MOD_DIR ${CMAKE_BINARY_DIR})

# Fortran runtime + implicit libs needed when a C++-driven executable pulls in
# Fortran objects from libnoahmp (libgfortran, libquadmath, ...).
set(NOAHMP_FORTRAN_RT ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

# ---------------------------------------------------------------------------
# Fortran test-support helpers (bind(C)), compiled once into their own static
# lib and linked into the C++ tests that need to peek at Fortran-owned state.
# NOT part of libnoahmp. A private module dir avoids .mod clashes.
# ---------------------------------------------------------------------------
add_library(noahmp_testsupport STATIC test_support.F90)
add_dependencies(noahmp_testsupport noahmp)
target_link_libraries(noahmp_testsupport PUBLIC noahmp)
target_include_directories(noahmp_testsupport PRIVATE ${NOAHMP_MOD_DIR})
set_target_properties(noahmp_testsupport PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tsmod)
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(noahmp_testsupport PRIVATE "-ffree-line-length-none")
endif()

# ---------------------------------------------------------------------------
# Helper: a C++ regression test executable linking the library (+ Fortran rt).
#   noahmp_add_cxx_test(<name> [WITH_SUPPORT])
# ---------------------------------------------------------------------------
function(noahmp_add_cxx_test name)
    cmake_parse_arguments(ARG "WITH_SUPPORT" "" "" ${ARGN})
    add_executable(${name} ${name}.cpp)
    target_compile_features(${name} PRIVATE cxx_std_17)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    if(ARG_WITH_SUPPORT)
        target_link_libraries(${name} PRIVATE noahmp_testsupport)
    endif()
    target_link_libraries(${name} PRIVATE noahmp ${NOAHMP_FORTRAN_RT})
endfunction()

# ---------------------------------------------------------------------------
# Tier 1 -- coupling machinery (no external data)
# ---------------------------------------------------------------------------
noahmp_add_cxx_test(test_array_view)
noahmp_add_cxx_test(test_abi_guard)
noahmp_add_cxx_test(test_fatal_handler)
noahmp_add_cxx_test(test_lifecycle   WITH_SUPPORT)
noahmp_add_cxx_test(test_sharing     WITH_SUPPORT)

add_test(NAME array_view    COMMAND test_array_view)
add_test(NAME abi_guard     COMMAND test_abi_guard)
add_test(NAME fatal_handler COMMAND test_fatal_handler)
add_test(NAME lifecycle     COMMAND test_lifecycle)
add_test(NAME sharing       COMMAND test_sharing)

# C++ sizing-guard aborts. test_abort_handler.cpp turns the fatal std::abort()
# into a non-zero std::_Exit so CTest's WILL_FAIL inverts it reliably (a SIGABRT
# is reported as "Subprocess aborted" and is NOT inverted by WILL_FAIL).
add_executable(test_cpp_guards test_cpp_guards.cpp test_abort_handler.cpp)
target_compile_features(test_cpp_guards PRIVATE cxx_std_17)
target_include_directories(test_cpp_guards PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_cpp_guards PRIVATE noahmp ${NOAHMP_FORTRAN_RT})
add_test(NAME cpp_guard_double_resize COMMAND test_cpp_guards double)
add_test(NAME cpp_guard_zero_resize   COMMAND test_cpp_guards zero)
set_tests_properties(cpp_guard_double_resize cpp_guard_zero_resize
    PROPERTIES WILL_FAIL TRUE)

# ---- Fortran-side executables ----
add_executable(test_fortran_alloc test_fortran_alloc.F90)
add_dependencies(test_fortran_alloc noahmp)
target_include_directories(test_fortran_alloc PRIVATE ${NOAHMP_MOD_DIR})
target_link_libraries(test_fortran_alloc PRIVATE noahmp)
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(test_fortran_alloc PRIVATE "-ffree-line-length-none")
endif()
add_test(NAME fortran_alloc COMMAND test_fortran_alloc)

# Water output-transfer units/accumulation contract: drives the real
# WaterVarOutTransfer and pins the [mm per soil timestep] depth semantics
# (QTDRAIN/SFCRUNOFF accumulate a depth, RUNSFXY/RUNSBXY snapshot, glacier
# sentinel zeroes TileDrain + scales runoff by MainTimeStep).
add_executable(test_water_flux_units test_water_flux_units.F90)
add_dependencies(test_water_flux_units noahmp)
target_include_directories(test_water_flux_units PRIVATE ${NOAHMP_MOD_DIR})
target_link_libraries(test_water_flux_units PRIVATE noahmp)
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(test_water_flux_units PRIVATE "-ffree-line-length-none")
endif()
add_test(NAME water_flux_units COMMAND test_water_flux_units)

# test_abort_handler.cpp routes NoahmpIO_abort's std::abort() to a non-zero
# std::_Exit so the WILL_FAIL cases below invert cleanly (see note above).
add_executable(test_fortran_guards test_fortran_guards.F90 test_abort_handler.cpp)
add_dependencies(test_fortran_guards noahmp)
target_include_directories(test_fortran_guards PRIVATE ${NOAHMP_MOD_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_fortran_guards PRIVATE noahmp)
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(test_fortran_guards PRIVATE "-ffree-line-length-none")
endif()
add_test(NAME fortran_guard_level   COMMAND test_fortran_guards level)
add_test(NAME fortran_guard_nblocks COMMAND test_fortran_guards nblocks)
add_test(NAME fortran_guard_reinit  COMMAND test_fortran_guards reinit)
set_tests_properties(fortran_guard_level fortran_guard_nblocks fortran_guard_reinit
    PROPERTIES WILL_FAIL TRUE)

# ---------------------------------------------------------------------------
# Tier 2 -- config smoke (shipped TBL + tiny namelist, no NetCDF)
# Stage the data files the parsers read from CWD into the test binary dir.
# NoahmpTable.TBL is a symlink in the driver dir; file(COPY) copies its target.
# ---------------------------------------------------------------------------
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/namelist.erf
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# NoahmpTable.TBL is a symlink in the driver dir; FOLLOW_SYMLINK_CHAIN copies the
# real file so the staged copy is not a dangling link in the build tree.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NoahmpTable.TBL
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FOLLOW_SYMLINK_CHAIN)

noahmp_add_cxx_test(test_config_smoke WITH_SUPPORT)
add_test(NAME config_smoke COMMAND test_config_smoke)
set_tests_properties(config_smoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ---------------------------------------------------------------------------
# Tier 3 -- NetCDF I/O (read / write / restart / checkpoint) + fatal path.
#
# These drive the real NetCDF modules end to end. The read path is serial
# NetCDF against a wrfinput/WPS file; the write and restart paths use collective
# parallel NetCDF-4 (NF90_MPIIO), so each test PROGRAM creates its own MPI world
# (MPI_COMM_WORLD) and runs single-rank. They therefore need the Fortran NetCDF
# headers/libs and MPI directly (libnoahmp links NetCDF PRIVATE and MPI PUBLIC).
#
# Genericity: the read test synthesizes its own tiny wrfinput fixture at run time
# (no committed multi-hundred-MB data file), so it runs anywhere; point it at any
# real local WPS/wrfinput file instead with NOAHMP_TEST_WRFINPUT=<path>.
# ---------------------------------------------------------------------------

# Shared I/O test-support (fixture generation + block setup + MPI world), built
# once into its own static lib with a private module dir (as noahmp_testsupport).
add_library(noahmp_io_testsupport STATIC test_io_support.F90)
add_dependencies(noahmp_io_testsupport noahmp)
target_link_libraries(noahmp_io_testsupport PUBLIC noahmp ${NETCDF_LIBRARIES_F90})
target_include_directories(noahmp_io_testsupport PRIVATE
    ${NOAHMP_MOD_DIR} ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_F90})
set_target_properties(noahmp_io_testsupport PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/iomod)
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(noahmp_io_testsupport PRIVATE "-ffree-line-length-none")
endif()

# Helper: a Fortran I/O test linking the library, NetCDF, and the I/O support lib.
# Extra sources (e.g. test_abort_handler.cpp for WILL_FAIL abort cases) via ARGN.
function(noahmp_add_io_test name)
    add_executable(${name} ${name}.F90 ${ARGN})
    add_dependencies(${name} noahmp noahmp_io_testsupport)
    target_include_directories(${name} PRIVATE
        ${NOAHMP_MOD_DIR} ${CMAKE_CURRENT_BINARY_DIR}/iomod
        ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_F90} ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(${name} PRIVATE
        noahmp_io_testsupport noahmp ${NETCDF_LIBRARIES_F90})
    if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
        target_compile_options(${name} PRIVATE "-ffree-line-length-none")
    endif()
endfunction()

# Read path (serial NetCDF): synthesize a wrfinput fixture and read it back.
noahmp_add_io_test(test_io_readland)
add_test(NAME io_readland COMMAND test_io_readland)
set_tests_properties(io_readland PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Land-output write (parallel NetCDF): write then reopen and compare.
noahmp_add_io_test(test_io_writeland)
add_test(NAME io_writeland COMMAND test_io_writeland)
set_tests_properties(io_writeland PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Checkpoint/restart (parallel NetCDF): bit-exact round trip + a WILL_FAIL
# layer-mismatch abort. test_abort_handler.cpp routes the mismatch abort to a
# non-zero _Exit so WILL_FAIL inverts it cleanly.
noahmp_add_io_test(test_io_restart test_abort_handler.cpp)
add_test(NAME io_restart_roundtrip COMMAND test_io_restart roundtrip)
add_test(NAME io_restart_mismatch  COMMAND test_io_restart mismatch)
set_tests_properties(io_restart_roundtrip io_restart_mismatch PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(io_restart_mismatch PROPERTIES WILL_FAIL TRUE)

# Full cold-init chain + NoahmpDriverMain over two steps, from a synthetic (or
# external) wrfinput. Needs the staged namelist.erf + NoahmpTable.TBL in CWD.
noahmp_add_io_test(test_io_driver)
add_test(NAME io_driver COMMAND test_io_driver)
set_tests_properties(io_driver PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Same full step, but driven through the PUBLIC C++ API in ERF's exact call order
# (ERF_NOAHMP_Init.cpp / ERF_NOAHMP_Advance.cpp). C++ main linking the library +
# the Fortran I/O support lib (MPI world, fixture, Fortran-only field shims).
add_executable(test_io_driver_cpp test_io_driver_cpp.cpp)
add_dependencies(test_io_driver_cpp noahmp noahmp_io_testsupport)
target_compile_features(test_io_driver_cpp PRIVATE cxx_std_17)
target_include_directories(test_io_driver_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_io_driver_cpp PRIVATE
    noahmp_io_testsupport noahmp ${NETCDF_LIBRARIES_F90} ${NOAHMP_FORTRAN_RT})
add_test(NAME io_driver_cpp COMMAND test_io_driver_cpp)
set_tests_properties(io_driver_cpp PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# NoahmpFatalMod: check_nc success returns; error status and NoahmpIO_abort
# terminate (WILL_FAIL via the same _Exit-on-fatal handler).
noahmp_add_io_test(test_io_abort test_abort_handler.cpp)
add_test(NAME io_abort_check_ok  COMMAND test_io_abort ok)
add_test(NAME io_abort_check_bad COMMAND test_io_abort checkbad)
add_test(NAME io_abort_direct    COMMAND test_io_abort abort)
set_tests_properties(io_abort_check_bad io_abort_direct PROPERTIES WILL_FAIL TRUE)
