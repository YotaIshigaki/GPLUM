#======================================================================
#   Numerical Libraries and Compilers
#======================================================================
QD_INC = -I/nfshome/namekata/gnu/include
QD_LIB = -L/nfshome/namekata/gnu/lib -lqd -lm -lrt
# /nfshome/namekata/gnu/bin/qd-config --cxxflags
# /nfshome/namekata/gnu/bin/qd-config --libs
GSL_INC = -I/nfshome/namekata/gnu/gsl-2.1/include
GSL_LIB = -L/nfshome/namekata/gnu/gsl-2.1/lib -lgsl -lgslcblas
EIGEN_INC = -I/nfshome/namekata/gnu/include/eigen3
EIGEN_LIB = 
FFTW_INC = -I/nfshome/namekata/gnu/include
FFTW_LIB = -L/nfshome/namekata/gnu/lib -lfftw3_mpi -lfftw3_omp -lfftw3 -lm
FDPS_INC = -I../../../../src -I../../../../src/particle_mesh
FDPS_LIB = -L../../../../src/particle_mesh -lpm
PMMM_INC = -I/nfshome/namekata/codes/fdps/sandbox/namekata_box/codes/PMMM
PMMM_LIB = 
LDFLAGS = $(QD_LIB) $(GSL_LIB) $(EIGEN_LIB) $(FFTW_LIB) $(FDPS_LIB) $(PMMM_LIB)

CXXFLAGS_COMMON = -std=c++11 -O3 -mavx -mfma4 -ffast-math -funroll-loops $(QD_INC) $(GSL_INC) $(EIGEN_INC) $(FFTW_INC) $(FDPS_INC) $(PMMM_INC)
#CXXFLAGS_COMMON = -std=c++11 -O0 -Wall -Wextra -ftrapv -fexceptions -g3 $(QD_INC) $(GSL_INC) $(EIGEN_INC) $(FFTW_INC) $(FDPS_INC) $(PMMM_INC)
#CXXFLAGS_COMMON = -std=c++11 -O2 -mavx -mfma4 -Wall -Wextra -ftrapv -fexceptions -g3 $(QD_INC) $(GSL_INC) $(EIGEN_INC) $(FFTW_INC) $(FDPS_INC) $(PMMM_INC)
# [1] Serial
#CXX = g++
#CXXFLAGS = $(CFLAGS_COMMON) 
#CXXFLAGS = $(CFLAGS_COMMON) -DPARTICLE_SIMULATOR_DEBUG_PRINT
# [2] OpenMP
#CXX = g++
#CXXFLAGS = $(CFLAGS_COMMON) -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# [3] MPI
#CXX = mpicxx
#CXXFLAGS = $(CFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL 
# [4] MPI + OpenMP
CXX = mpicxx
CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

#----------------------------------------------------------------------
#   WORKDIR
#----------------------------------------------------------------------
WORKDIR = /work2/namekata
