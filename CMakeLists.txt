# Shared library for the Python package (ssfortran) and other C callers.
# The Fortran-only build and the test suite use fpm (fpm.toml); this build
# compiles the same sources with -fPIC into libstatespace, with L-BFGS-B
# fetched at the commit fpm.toml pins.
#
#   cmake -S . -B build/cmake -G Ninja && cmake --build build/cmake
#
# Python wheels are built through scikit-build-core (pyproject.toml).
cmake_minimum_required(VERSION 3.24)
project(statespace VERSION 0.1.0 LANGUAGES Fortran)

option(STATESPACE_OPENMP "Build fit_many with OpenMP" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)
FetchContent_Declare(
  lbfgsb
  GIT_REPOSITORY https://github.com/jacobwilliams/lbfgsb.git
  GIT_TAG ce1a9322016d6f31b652393e000a914b63af6fef
  SOURCE_SUBDIR no-cmake  # only the sources are used
)
FetchContent_MakeAvailable(lbfgsb)

file(GLOB STATESPACE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)
set(LBFGSB_SOURCES
  ${lbfgsb_SOURCE_DIR}/src/lbfgsb_kinds_module.F90
  ${lbfgsb_SOURCE_DIR}/src/lbfgsb_blas_module.F90
  ${lbfgsb_SOURCE_DIR}/src/lbfgsb_linpack_module.f90
  ${lbfgsb_SOURCE_DIR}/src/lbfgsb.f90
)

add_library(statespace SHARED ${STATESPACE_SOURCES} ${LBFGSB_SOURCES})
set_target_properties(statespace PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(statespace PRIVATE
    $<$<CONFIG:Release>:-O3 -funroll-loops>
    -ffree-line-length-none
  )
endif()

find_package(LAPACK REQUIRED)
target_link_libraries(statespace PRIVATE ${LAPACK_LIBRARIES})

if(STATESPACE_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(statespace PRIVATE OpenMP::OpenMP_Fortran)
  endif()
endif()

# The wheel puts the library inside the Python package.
if(SKBUILD)
  install(TARGETS statespace LIBRARY DESTINATION ssfortran RUNTIME DESTINATION ssfortran)
else()
  install(TARGETS statespace)
endif()
