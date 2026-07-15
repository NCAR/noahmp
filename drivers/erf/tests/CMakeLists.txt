# Regression tests for the Noah-MP ERF driver: link the library and exercise the
# C++ <-> Fortran coupling boundary through both interfaces. Run with `ctest`.

# Tee each test's stdout+stderr to console.log in its WORKING_DIRECTORY; pipefail
# preserves the real exit status and stdout for CTest PASS/FAIL regex scanning.
function(add_test)
    cmake_parse_arguments(PARSE_ARGV 0 _AT "" "NAME" "COMMAND")
    list(POP_FRONT _AT_COMMAND _exe)
    if(TARGET ${_exe})
        set(_exe $<TARGET_FILE:${_exe}>)
    endif()
    _add_test(NAME ${_AT_NAME} COMMAND
        bash -c "set -o pipefail; \"$@\" 2>&1 | tee console.log" wrap ${_exe} ${_AT_COMMAND})
endfunction()

# Where the library's generated Fortran modules (noahmpio_fi.mod, ...) live.
set(NOAHMP_MOD_DIR ${CMAKE_BINARY_DIR})

# Fortran runtime + implicit libs needed when a C++-driven executable pulls in
# Fortran objects from libnoahmp (libgfortran, libquadmath, ...).
set(NOAHMP_FORTRAN_RT ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

# GNU gfortran needs unlimited free-form line length for these long-lined sources.
function(noahmp_fortran_flags target)
    if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
        target_compile_options(${target} PRIVATE "-ffree-line-length-none")
    endif()
endfunction()

# Fortran test-support helpers (bind(C)) in their own static lib for C++ tests
# that peek at Fortran-owned state. Not part of libnoahmp; private module dir.
add_library(noahmp_testsupport STATIC test_support.F90)
add_dependencies(noahmp_testsupport noahmp)
target_link_libraries(noahmp_testsupport PUBLIC noahmp)
target_include_directories(noahmp_testsupport PRIVATE ${NOAHMP_MOD_DIR})
set_target_properties(noahmp_testsupport PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tsmod)
noahmp_fortran_flags(noahmp_testsupport)

# Helpers: noahmp_add_cxx_test(<name> [WITH_SUPPORT] [SOURCES ...]) -- C++ + rt,
# noahmp_add_fortran_test(<name> [SOURCES ...]) -- Fortran + noahmp.
function(noahmp_add_cxx_test name)
    cmake_parse_arguments(ARG "WITH_SUPPORT" "" "SOURCES" ${ARGN})
    add_executable(${name} ${name}.cpp ${ARG_SOURCES})
    target_compile_features(${name} PRIVATE cxx_std_17)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    if(ARG_WITH_SUPPORT)
        target_link_libraries(${name} PRIVATE noahmp_testsupport)
    endif()
    target_link_libraries(${name} PRIVATE noahmp ${NOAHMP_FORTRAN_RT})
endfunction()

function(noahmp_add_fortran_test name)
    cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})
    add_executable(${name} ${name}.F90 ${ARG_SOURCES})
    add_dependencies(${name} noahmp)
    target_include_directories(${name} PRIVATE ${NOAHMP_MOD_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(${name} PRIVATE noahmp)
    noahmp_fortran_flags(${name})
endfunction()

# Per-test artifact isolation: each test runs in its own artifacts/<name>/ dir.
# CONFIG_INPUTS stages namelist.erf + NoahmpTable.TBL; call AFTER add_test().
set(NOAHMP_TEST_ARTIFACTS_DIR ${CMAKE_CURRENT_BINARY_DIR}/artifacts)
function(noahmp_test_workdir test_name)
    cmake_parse_arguments(ARG "CONFIG_INPUTS" "" "" ${ARGN})
    set(_wd ${NOAHMP_TEST_ARTIFACTS_DIR}/${test_name})
    file(MAKE_DIRECTORY ${_wd})
    if(ARG_CONFIG_INPUTS)
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/namelist.erf DESTINATION ${_wd})
        # NoahmpTable.TBL is a symlink in the driver dir; FOLLOW_SYMLINK_CHAIN
        # copies the real file so the staged copy is not a dangling link.
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../NoahmpTable.TBL
             DESTINATION ${_wd} FOLLOW_SYMLINK_CHAIN)
    endif()
    set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY ${_wd})
endfunction()

# Tier 1 -- coupling machinery (no external data).
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

# C++ sizing-guard aborts, detected by OUTPUT not exit code: each must print its
# guard diagnostic (PASS) and not hit the fall-through (FAIL). Stricter than WILL_FAIL.
noahmp_add_cxx_test(test_cpp_guards SOURCES test_abort_handler.cpp)
add_test(NAME cpp_guard_double_resize COMMAND test_cpp_guards double)
add_test(NAME cpp_guard_zero_resize   COMMAND test_cpp_guards zero)
set_tests_properties(cpp_guard_double_resize PROPERTIES
    PASS_REGULAR_EXPRESSION "already been initialized and cannot be resized"
    FAIL_REGULAR_EXPRESSION "expected abort did NOT occur")
set_tests_properties(cpp_guard_zero_resize PROPERTIES
    PASS_REGULAR_EXPRESSION "size must be > 0"
    FAIL_REGULAR_EXPRESSION "expected abort did NOT occur")

# ---- Fortran-side executables ----
noahmp_add_fortran_test(test_fortran_alloc)
add_test(NAME fortran_alloc COMMAND test_fortran_alloc)

# Water output-transfer units/accumulation contract: drives WaterVarOutTransfer
# and pins the [mm per soil timestep] depth semantics.
noahmp_add_fortran_test(test_water_flux_units)
add_test(NAME water_flux_units COMMAND test_water_flux_units)

# test_abort_handler.cpp routes NoahmpIO_abort's std::abort() to a clean _Exit so
# the diagnostic flushes; cases below are detected by output, not exit code.
noahmp_add_fortran_test(test_fortran_guards SOURCES test_abort_handler.cpp)
add_test(NAME fortran_guard_level   COMMAND test_fortran_guards level)
add_test(NAME fortran_guard_nblocks COMMAND test_fortran_guards nblocks)
add_test(NAME fortran_guard_reinit  COMMAND test_fortran_guards reinit)
# Output-based detection (see the cpp-guard note above): each case must print its
# validation diagnostic and not fall through. Stricter than WILL_FAIL.
set_tests_properties(fortran_guard_level PROPERTIES
    PASS_REGULAR_EXPRESSION "NoahmpIOTypeVectInit_fi: level out of range"
    FAIL_REGULAR_EXPRESSION "expected abort did NOT occur")
set_tests_properties(fortran_guard_nblocks PROPERTIES
    PASS_REGULAR_EXPRESSION "NBlocks must be"
    FAIL_REGULAR_EXPRESSION "expected abort did NOT occur")
set_tests_properties(fortran_guard_reinit PROPERTIES
    PASS_REGULAR_EXPRESSION "refusing to re-init"
    FAIL_REGULAR_EXPRESSION "expected abort did NOT occur")

# Tier 2 -- config smoke (shipped TBL + tiny namelist, no NetCDF).
# CONFIG_INPUTS stages namelist.erf + NoahmpTable.TBL, read from CWD.
noahmp_add_cxx_test(test_config_smoke WITH_SUPPORT)
add_test(NAME config_smoke COMMAND test_config_smoke)

# Tier 3 -- NetCDF I/O (read / write / restart / checkpoint) + fatal path. Drives
# the real NetCDF modules end to end; each test builds its own single-rank MPI world.

# Shared I/O test-support (fixture generation + block setup + MPI world), built
# once into its own static lib with a private module dir (as noahmp_testsupport).
add_library(noahmp_io_testsupport STATIC test_io_support.F90)
add_dependencies(noahmp_io_testsupport noahmp)
target_link_libraries(noahmp_io_testsupport PUBLIC noahmp ${NETCDF_LIBRARIES_F90})
target_include_directories(noahmp_io_testsupport PRIVATE
    ${NOAHMP_MOD_DIR} ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_F90})
set_target_properties(noahmp_io_testsupport PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/iomod)
noahmp_fortran_flags(noahmp_io_testsupport)

# Helper: a Fortran I/O test linking the library, NetCDF, and the I/O support lib.
# Extra sources (e.g. test_abort_handler.cpp for expected-abort cases) via ARGN.
function(noahmp_add_io_test name)
    add_executable(${name} ${name}.F90 ${ARGN})
    add_dependencies(${name} noahmp noahmp_io_testsupport)
    target_include_directories(${name} PRIVATE
        ${NOAHMP_MOD_DIR} ${CMAKE_CURRENT_BINARY_DIR}/iomod
        ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_F90} ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(${name} PRIVATE
        noahmp_io_testsupport noahmp ${NETCDF_LIBRARIES_F90})
    noahmp_fortran_flags(${name})
endfunction()

# Read path (serial NetCDF): synthesize a wrfinput fixture and read it back.
noahmp_add_io_test(test_io_readland)
add_test(NAME io_readland COMMAND test_io_readland)

# Land-output write (parallel NetCDF): write then reopen and compare.
noahmp_add_io_test(test_io_writeland)
add_test(NAME io_writeland COMMAND test_io_writeland)

# Checkpoint/restart (parallel NetCDF): bit-exact round trip + a layer-mismatch
# abort detected by output (see cpp-guard note); test_abort_handler.cpp linked.
noahmp_add_io_test(test_io_restart test_abort_handler.cpp)
add_test(NAME io_restart_roundtrip COMMAND test_io_restart roundtrip)
add_test(NAME io_restart_mismatch  COMMAND test_io_restart mismatch)
# The "layer mismatch" line is stdout; _Exit(7) skips libgfortran's flush, so
# GFORTRAN_UNBUFFERED_ALL keeps the PASS token from being dropped before CTest scans.
set_tests_properties(io_restart_mismatch PROPERTIES
    ENVIRONMENT "GFORTRAN_UNBUFFERED_ALL=y"
    PASS_REGULAR_EXPRESSION "layer mismatch"
    FAIL_REGULAR_EXPRESSION "expected layer-mismatch abort did NOT occur")

# Full cold-init chain + NoahmpDriverMain over two steps, from a synthetic
# wrfinput. Needs the staged namelist.erf + NoahmpTable.TBL in CWD (CONFIG_INPUTS).
noahmp_add_io_test(test_io_driver)
add_test(NAME io_driver COMMAND test_io_driver)

# Checkpoint/restart REPRODUCIBILITY across time steps: drive NoahmpDriverMain for
# a continuous run (steps 1..N) vs a checkpoint-at-(N-1) + restart + step-N run,
# and require the final state to match bit-for-bit. Catches prognostic state the
# driver evolves but the restart set fails to persist (test_io_restart alone
# cannot). One process sweeps every Noah-MP subsystem (land/rain/snow/carbon/
# glacier/groundwater), each driving its state to life and asserting (activity
# guard) it actually evolved. Needs staged namelist.erf + NoahmpTable.TBL
# (CONFIG_INPUTS). A single scenario can be run for debugging: `... <name>`.
noahmp_add_io_test(test_io_restart_step)
add_test(NAME io_restart_step COMMAND test_io_restart_step all)

# Same full step driven through the PUBLIC C++ API in ERF's exact call order.
# C++ main linking the library + Fortran I/O support lib.
add_executable(test_io_driver_cpp test_io_driver_cpp.cpp)
add_dependencies(test_io_driver_cpp noahmp noahmp_io_testsupport)
target_compile_features(test_io_driver_cpp PRIVATE cxx_std_17)
target_include_directories(test_io_driver_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_io_driver_cpp PRIVATE
    noahmp_io_testsupport noahmp ${NETCDF_LIBRARIES_F90} ${NOAHMP_FORTRAN_RT})
add_test(NAME io_driver_cpp COMMAND test_io_driver_cpp)

# NoahmpFatalMod: check_nc success returns; error status and NoahmpIO_abort
# terminate. Abort cases detected by output, not exit code (see cpp-guard note).
noahmp_add_io_test(test_io_abort test_abort_handler.cpp)
add_test(NAME io_abort_check_ok  COMMAND test_io_abort ok)
add_test(NAME io_abort_check_bad COMMAND test_io_abort checkbad)
add_test(NAME io_abort_direct    COMMAND test_io_abort abort)
# check_nc must abort on error status, printing the forced-error context (stdout).
# GFORTRAN_UNBUFFERED_ALL keeps it from being dropped by _Exit(7) before CTest scans.
set_tests_properties(io_abort_check_bad PROPERTIES
    ENVIRONMENT "GFORTRAN_UNBUFFERED_ALL=y"
    PASS_REGULAR_EXPRESSION "forced NetCDF error"
    FAIL_REGULAR_EXPRESSION "did NOT abort")
# NoahmpIO_abort must not return: the pre-call marker must appear and the
# post-call "returned unexpectedly" line must not.
set_tests_properties(io_abort_direct PROPERTIES
    PASS_REGULAR_EXPRESSION "reached direct-abort call"
    FAIL_REGULAR_EXPRESSION "returned unexpectedly")

# Give every test its own artifacts/<name>/ working dir. Runs after add_test();
# the test list comes from the TESTS property so it never drifts.
set(NOAHMP_CONFIG_INPUT_TESTS config_smoke io_driver io_driver_cpp io_restart_step)
get_property(_noahmp_tests DIRECTORY PROPERTY TESTS)
foreach(t ${_noahmp_tests})
    if(t IN_LIST NOAHMP_CONFIG_INPUT_TESTS)
        noahmp_test_workdir(${t} CONFIG_INPUTS)
    else()
        noahmp_test_workdir(${t})
    endif()
endforeach()
