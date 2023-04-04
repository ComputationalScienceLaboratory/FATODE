# FATODE
[![License](https://img.shields.io/github/license/ComputationalScienceLaboratory/FATODE.svg)](https://github.com/ComputationalScienceLaboratory/FATODE/blob/master/LICENSE)
[![CI Status](https://github.com/ComputationalScienceLaboratory/FATODE/actions/workflows/test.yml/badge.svg)](https://github.com/ComputationalScienceLaboratory/FATODE/actions/workflows/test.yml)

FATODE version 1.2
by Hong Zhang and Adrian Sandu
Released April 30, 2013

FATODE version 1.1
by Hong Zhang and Adrian Sandu
Released November 11, 2012

FATODE contains a set of integrators to solve ordinary differential system
y'=f(t,y) with capabilities of direct and adjoint sensitivity analysis.
FATODE is implemented in FORTRAN90, and has been tested with the following
compilers: Portland group's pgf90, Intel's ifort, Lahey's lf95, Sun's sunf90,
gfortran, Absoft.

FATODE contains the following directory structure:

    FATODE/README      general instructions
    FATODE/FWD         forward model integrators   
    FATODE/ADJ         adjoint model integrators
    FATODE/TLM         tangent linear model integrators
    FATODE/LSS_LIBS    lib files of linear solvers
    FATODE/test        example programs that also function as tests
    FATODE/DOC         documentation (user's guide)

FATODE implementation requires BLAS and LAPACK libraries. For solution of
linear systems arised in the integration, FATODE provides options to use BLAS/
LAPACK routines (with macro definition -DFULL_ALGEBRA) and to use sparse linear
solvers such as UMFPACK (-DSPARSE_UMF) and SuperLU (-DSPARSE_LU). If sparse
linear solvers are needed, the FATODE code must be linked with corresponding
(static) libraries for the linear solvers. For installation instuctions of
UMFPACK and SuperLU, please refer to their websites.
UMFPACK can be accessed at:
http://www.cise.ufl.edu/research/sparse/umfpack/
SuperLU can be accessed at:
http://crd.lbl.gov/~xiaoye/SuperLU/

Since UMFPACK and SuperLU are written in C, wrappers for calling C routines from
FORTRAN are needed. They can be found in the original packages of UMFPACK and
SuperLU. For convenience, they are also included in FATODE/LSS_LIBS directory as well as makefiles to generate object files which can be linked to FATODE code. Make sure to modify the installation PATH first in the makefiles if you need to regenerate the object files.

For usage of FATODE, please refer to the user's guide, as well as the example programs containing drivers and makefiles which can be used as templates.


## Changes

### Changes from Version 1.2.1
- The library is now built with CMake
- Dockerfiles that build the library with different algebra schemes are included

### Changes from Version 1.1
Bug fixes:
- Fixed a possible initilization problem of method coefficients for windows users.
- Fixed improper use of drdy in the adjoint Rosenbrock solver.
Complier options:
- The option -mp for ifort compiler is deprecated. Use -fltconsistency instead.


### Changes from Version 1.0
New example:
- CBM-IV: a stiff chemical system used to demonstrate the usage of FATODE for computing sensivities with respect to system parameters.
Bug fixes:
- Fixed several bugs in Rosenbrock solvers
Complier options:
- ifort sometimes reulsts in inconsistent output with other compilers due to floating-point calculation consistency. The option -mp can solve this problem.
Documentation:
- Add a user's guide for FATODE


# Installation

1. Clone the repository
2. Install the depencencies for the linear algebra library you would like to use. The options are SuiteSparse, SuperLU, or a full algebra solve with lapack. All versions require lapack.
3. Make a build directory `cd FATODE && mkdir build && cd build`
4. Configure the project with the options that you want. The most important options are the type of linear algebra library, SuiteSparse, SuperLU, or full algebra (driven by lapack).
    - SuiteSparse: `cmake -D USE_SUITESPARSE_UMF=ON ..`
    - SuperLU: `cmake -D USE_SUPER_LU=ON ..`
    - Full Algebra `cmake -D USE_FULL_ALGEBRA=ON ..`
5. Build `make`
6. Install `make install`

## Usage
After installing, you can find the libray with cmake like this

```
cmake_minimum_required(VERSION 3.21)

project(
  example
  VERSION 0.0
  LANGUAGES Fortran
)

find_package(LAPACK REQUIRED)

find_package(FATODE 1.2.1 REQUIRED)

add_executable(test
  cbm4_rk_dr.F90 cbm4_parameters.F90
)

target_link_libraries(test
  PUBLIC
    csl::FATODE
)
```
