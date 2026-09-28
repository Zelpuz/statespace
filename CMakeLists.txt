# Shared library for the Python package (ssfortran) and other C callers.
# The Fortran-only build and the test suite use fpm (fpm.toml); this build
# compiles the same sources with -fPIC into libstatespace, with the bundled
# L-BFGS-B that fpm also uses.
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

# L-BFGS-B, bundled in third_party/lbfgsb (see its README), so the build needs
# no network access.
set(LBFGSB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lbfgsb/src)
file(GLOB STATESPACE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)
set(LBFGSB_SOURCES
  ${LBFGSB_DIR}/lbfgsb_kinds_module.F90
  ${LBFGSB_DIR}/lbfgsb_blas_module.F90
  ${LBFGSB_DIR}/lbfgsb_linpack_module.f90
  ${LBFGSB_DIR}/lbfgsb.f90
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
