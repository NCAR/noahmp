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

# Coupling codegen (configure time): expand *-mc templates into the gitignored
# generated sources BEFORE the file(GLOB) below discovers them. Idempotent.
find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(NOAHMP_CODEGEN ${CMAKE_CURRENT_SOURCE_DIR}/tools/NoahmpMacro.py)
execute_process(
    COMMAND ${Python3_EXECUTABLE} ${NOAHMP_CODEGEN}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE NOAHMP_CODEGEN_RC)
if(NOT NOAHMP_CODEGEN_RC EQUAL 0)
    message(FATAL_ERROR "Noah-MP coupling codegen failed (tools/NoahmpMacro.py). "
                        "Generated NoahmpIO.{H,cpp}, NoahmpIO_fi.F90, "
                        "NoahmpIOVarType.F90 and NoahmpIOVarInitMod.F90 are "
                        "required to build this library.")
endif()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${NOAHMP_CODEGEN}
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIO.H-mc
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIO.cpp-mc
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIO_fi.F90-mc
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIOVarType.F90-mc
    ${CMAKE_CURRENT_SOURCE_DIR}/NoahmpIOVarInitMod.F90-mc)

file(GLOB NoahMP_Sources
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../../utility/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)
add_library(noahmp ${NoahMP_Sources})
add_library(NoahMP::noahmp ALIAS noahmp)

# Precision (DOUBLE_PREC -> kind_noahmp) is set by the parent project; ERF uses
# -DDOUBLE_PREC. Standalone falls back to upstream single unless you add it.

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(noahmp PRIVATE "-ffree-line-length-none")
endif()

target_include_directories(noahmp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(noahmp PRIVATE $<BUILD_INTERFACE:${NETCDF_INCLUDES}>)
if(NETCDF_INCLUDES_F90)
    target_include_directories(noahmp PRIVATE $<BUILD_INTERFACE:${NETCDF_INCLUDES_F90}>)
endif()
target_link_libraries(noahmp PRIVATE ${NETCDF_LIBRARIES_F90})
target_link_libraries(noahmp PUBLIC MPI::MPI_Fortran)
target_include_directories(noahmp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Convenience codegen targets: regenerate in place, or drift-check for CI.
add_custom_target(noahmp_codegen
    COMMAND ${Python3_EXECUTABLE} ${NOAHMP_CODEGEN}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Regenerating Noah-MP coupling glue from @NoahmpMacro:source")
add_custom_target(noahmp_codegen_check
    COMMAND ${Python3_EXECUTABLE} ${NOAHMP_CODEGEN} --check
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Checking Noah-MP coupling glue is in sync with @NoahmpMacro:source")

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

# Regression tests (drivers/erf/tests). Default ON standalone, OFF as an ERF
# subproject; override with -DNOAHMP_ENABLE_TESTS=ON.
option(NOAHMP_ENABLE_TESTS "Build Noah-MP ERF driver regression tests" ${NOAHMP_MASTER_PROJECT})
if(NOAHMP_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
