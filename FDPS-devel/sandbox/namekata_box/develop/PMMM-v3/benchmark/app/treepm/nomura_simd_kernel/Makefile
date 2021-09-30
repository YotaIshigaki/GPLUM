#======================================================================
#   Numerical Libraries and Compilers
#======================================================================

FFTW_LOC = /usr/local
FFTW_INC = -I$(FFTW_LOC)/include
FFTW_LIB = -L$(FFTW_LOC)/lib  -lfftw3f_mpi -lfftw3f_omp -lfftw3f -lm

FDPS_PATH ?= $(HOME)/FDPS
PM_PATH = $(FDPS_PATH)/src/particle_mesh
FDPS_INC = -I$(FDPS_PATH)/src -I$(PM_PATH)
FDPS_LIB = -L$(PM_PATH) -lpm

CXXFLAGS_COMMON = -std=c++11 -O2 -funroll-loops $(FDPS_INC) $(FFTW_INC) -I../inc
LDFLAGS += $(FDPS_LIB) $(FFTW_LIB)
#CXXFLAGS_COMMON += -g

CXXFLAGS_COMMON += -mavx512f -mavx512dq -DENABLE_AVX512
CXXFLAGS_COMMON += -DENABLE_SINGLE
CXXFLAGS_COMMON += -DTUNING # for time profile

# [1] Serial
CXX = g++
CXXFLAGS = $(CXXFLAGS_COMMON)
#CXXFLAGS+= -DPARTICLE_SIMULATOR_DEBUG_PRINT

# [2] OpenMP
#CXX = g++
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

# [3] MPI
MPI_PATH = /usr/local/opt/open-mpi
CXX = g++
CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -L$(MPI_PATH)/lib -I$(MPI_PATH)/include -lmpi

# [4] MPI + OpenMP
#CXX = mpicxx
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#LDFLAGS += -L$(MPI_PATH)/lib64
# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

#----------------------------------------------------------------------
#   Source and Object files
#----------------------------------------------------------------------
SRCS	= main.cpp
HEDS    = user_defined_class.h  kernel.hpp

TARGET = main.out

$(TARGET): Makefile  $(SRCS) $(HEDS)
	$(CXX) $(SRCS) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS) 
	
clean:
	rm -f $(TARGET) *.o

distclean: clean
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
