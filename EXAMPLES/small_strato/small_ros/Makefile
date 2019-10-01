############################################################################
#
##  Module:          Makefile
#
##  Purpose:         Top-level Makefile
#
##  Creation date:   Nov 27, 2010
#
##  Modified:
#
##  Send bug reports, comments or suggestions to zhang@vt.edu
#
#############################################################################
export ARCH=ar -cr
export RANLIB=ranlib

# Configuration flags for linear solvers. Must choose one from the following three options.
#  -DFULL_ALGEBRA            use BLAS and LAPACK library (full algebra)
#  -DSPARSE_UMF              use UMFpack (sparse format) 
#  -DSPARSE_LU               use SuperLU (sparse format)
# 
LS_CONFIG  =  -DFULL_ALGEBRA

#~~~> NAG fortran
#export FC=nagfor
#export FFLAGS = -O2 $(LS_CONFIG)
#BLAS = ~/nagfor_libs/libblas.a
#LAPACK = ~/nagfor_libs/liblapack.a

#~~~> absoft fortran
#export FC=f90
#export FFLAGS = -O2 $(LS_CONFIG)
#BLAS = /Applications/Absoft10.1/lib/libblas.a
#LAPACK = /Applications/Absoft10.1/lib/liblapack.a

#~~~> intel fortran
#export FC=ifort
#export FFLAGS= -cpp -O2 -nogen-interface -fltconsitency $(LS_CONFIG)
#BLAS = /opt/ifort_libs/libblas.a
#LAPACK = /opt/ifort_libs/liblapack.a


#~~~> sunf90
#export FC = /opt/oracle/solstudio12.2/bin/sunf90
#export FFLAGS = -fpp -O2 -free $(LS_CONFIG)
#BLAS =  -xlic_lib=sunperf
#LAPACK =

#~~~> gfortran (GNU FORTRAN Compiler)
export FC = gfortran
export FFLAGS = -cpp -O2 -ffree-line-length-none $(LS_CONFIG)
BLAS=/opt/lapack/lib/libblas.a
LAPACK=/opt/lapack/lib/liblapack.a


#~~~> PGF90 (Portland Group Compiler)
#export FC=pgf90
#export FFLAGS=-O2 -fastsse -Mcache_align -tp=penryn-64 -Mflushz -Minform=warn -Mpreprocess $(LS_CONFIG)
#BLAS=/opt/pgi/linux86-64/7.2-5/lib/libblas.a
#LAPACK=/opt/pgi/linux86-64/7.2-5/lib/liblapack.a
#other libraries
#XERBLA = -lpgftnrtl -pgf90libs


#~~~> LAHEY
#export FC  = lf95
#export FFLAGS = -O2 -Cpp $(LS_CONFIG)
#BLAS=-lssl2mt
#LAPACK=


FATDIR = ../../..
MODEL   = FWD
FAMILY  = ROS

APP = small_strato 
PAR = User_Parameters.o

LSSOLVER=LS_Solver.o

INTEGRATOR=$(FAMILY)_f90_Integrator.o

default: driver

User_Parameters.o: User_Parameters.F90
	$(FC) $(FFLAGS) -c $<

LS_Solver.o: $(FATDIR)/$(MODEL)/$(FAMILY)/LS_Solver.F90
	$(FC) $(FFLAGS) -c $<

$(FAMILY)_f90_Integrator.o: $(FATDIR)/$(MODEL)/$(FAMILY)/$(FAMILY)_f90_Integrator.F90
	$(FC) $(FFLAGS) -c $<

small_strato_dr.o: small_strato_dr.F90 $(PAR) $(LSSOLVER) $(INTEGRATOR) 
	$(FC) $(FFLAGS) -c $<


LIB = -lm

#UMFPACK = $(FATDIR)/LSS_LIBS/libumfpack.a $(FATDIR)/LSS_LIBS/libamd.a
#UMFWRAP = $(FATDIR)/LSS_LIBS/umf4_f77wrapper.o

#LUPACK = $(FATDIR)/LSS_LIBS/libsuperlu_4.0.a
#LUWRAP = $(FATDIR)/LSS_LIBS/c_fortran_dgssv.o

default: driver clean

driver: small_strato_dr.o $(PAR) $(LSSOLVER) $(INTEGRATOR)
	$(FC) $(FFLAGS) -o $(APP) $< $(LUWRAP) $(LUPACK) $(UMFWRAP) $(UMFPACK) $(INTEGRATOR) $(LSSOLVER) $(PAR) $(LAPACK) $(BLAS) $(XERBLA) $(LIB)

purge:  clean

clean:
	rm -f *~ *.o *.mod 
	rm -f driver
help:
	@$(ECHO) "usage: make ?"

