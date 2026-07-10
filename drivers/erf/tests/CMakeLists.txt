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
