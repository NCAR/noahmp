cmake_minimum_required(VERSION 3.17)
project(NoahMP LANGUAGES Fortran CXX C)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

set(NOAHMP_MASTER_PROJECT OFF)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(NOAHMP_MASTER_PROJECT ON)
endif()

set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

find_package(MPI COMPONENTS Fortran REQUIRED)

# ---------------------------------------------------------------------------
# Coupling code generation (configure time). The @NoahmpIOCoupling:<region>
# blocks in NoahmpIO.{H,cpp}, NoahmpIO_fi.F90 and NoahmpIOVarType.F90 are
# expanded from the @NoahmpIOCoupling:source block in the tracked *-mc
# templates. The generated files are gitignored (not committed), so they MUST be
# produced before the file(GLOB) below discovers them as sources. The generator
# is stdlib-only Python and idempotent. Editing a template or the generator
# re-triggers configure (CMAKE_CONFIGURE_DEPENDS), which regenerates them.
# ---------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(NOAHMP_CODEGEN ${CMAKE_CURRENT_SOURCE_DIR}/tools/NoahmpIOCoupling.py)
execute_process(
    COMMAND ${Python3_EXECUTABLE} ${NOAHMP_CODEGEN}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE NOAHMP_CODEGEN_RC)
if(NOT NOAHMP_CODEGEN_RC EQUAL 0)
    message(FATAL_ERROR "Noah-MP coupling codegen failed (tools/NoahmpIOCoupling.py). "
                        "Generated NoahmpIO.{H,cpp}, NoahmpIO_fi.F90 and "
                        "NoahmpIOVarType.F90 are required to build this library.")
endif()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${NOAHMP_CODEGEN}
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIO.H-mc
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIO.cpp-mc
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIO_fi.F90-mc
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIOVarType.F90-mc)

file(GLOB NoahMP_Sources
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../../utility/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)
add_library(noahmp ${NoahMP_Sources})
add_library(NoahMP::noahmp ALIAS noahmp)

# Floating-point precision (the DOUBLE_PREC macro -> kind_noahmp) is selected by
# the PARENT project, not here. ERF compiles this library with -DDOUBLE_PREC
# because the ERF<->Noah-MP C++ interface (NoahArray<double>, NoahmpIO_type) is
# double-precision only; see the ERF root CMakeLists.txt. When this directory is
# built standalone, precision falls back to Noah-MP's upstream default (single);
# add the definition explicitly for a double-precision standalone build:
#   target_compile_definitions(noahmp PRIVATE DOUBLE_PREC)

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(noahmp PRIVATE "-ffree-line-length-none")
endif()

target_include_directories(noahmp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(noahmp PRIVATE $<BUILD_INTERFACE:${NETCDF_INCLUDES}>)
target_link_libraries(noahmp PRIVATE ${NETCDF_LIBRARIES_F90})
target_link_libraries(noahmp PUBLIC MPI::MPI_Fortran)
target_include_directories(noahmp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# ---------------------------------------------------------------------------
# Convenience codegen targets. The glue is already (re)generated at configure
# time above; these let you regenerate or drift-check on demand:
#   cmake --build . --target noahmp_codegen        # regenerate in place
#   cmake --build . --target noahmp_codegen_check   # fail on drift (for CI)
# ---------------------------------------------------------------------------
add_custom_target(noahmp_codegen
    COMMAND ${Python3_EXECUTABLE} ${NOAHMP_CODEGEN}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Regenerating Noah-MP coupling glue from @NoahmpIOCoupling:source")
add_custom_target(noahmp_codegen_check
    COMMAND ${Python3_EXECUTABLE} ${NOAHMP_CODEGEN} --check
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Checking Noah-MP coupling glue is in sync with @NoahmpIOCoupling:source")

# Install the library and export the target
install(TARGETS noahmp
    EXPORT NoahMPTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)

install(EXPORT NoahMPTargets
    FILE NoahMPTargets.cmake
    NAMESPACE NoahMP::
    DESTINATION lib/cmake/noahmp
)
