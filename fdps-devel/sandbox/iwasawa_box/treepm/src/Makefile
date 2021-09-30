#CC := mpicc
CC := time mpifccpx
#CXX := mpicxx
CXX = time mpiFCCpx

CFLAGS = -Kfast
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -x32
CFLAGS += -Xg
CFLAGS += -DFAST_ALL_TO_ALL_FOR_K
CFLAGS += -DFAST_WALK_K
CFLAGS += -DFFTW_DOUBLE

#CFLAGS = -O3 -ffast-math -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL


CXXFLAGS := $(CFLAGS)

all: treepm makeglass

TREEPM_OBJ = treepm.o make_directory.o timing.o
TREEPM_DEP := $(TREEPM_OBJ) treepm.hpp cosmology.hpp run_param.hpp

MAKEGLASS_OBJ = makeglass.o make_directory.o timing.o
MAKEGLASS_DEP := $(MAKEGLASS_OBJ) treepm.hpp cosmology.hpp run_param.hpp

#FDPS_ROOT=/home/masaki/project/fdps
#FDPS_ROOT=../..
#FDPS_ROOT=../../../../fdps/
FDPS_ROOT=/home/ra000008/a03080/project/fdps
#PM_ROOT=$(FDPS_ROOT)/src/particle_mesh
PM_ROOT=$(FDPS_ROOT)/src/particle_meshD

CXXFLAGS += -I$(FDPS_ROOT)/src
CXXFLAGS += -I$(PM_ROOT)

#FFTW_ROOT = /usr/local/fftw-3
#FFTW_ROOT = ../../../../../../
#FFTW_ROOT = /home/masaki/opt
#FFTW_LIBS = -L$(FFTW_ROOT)/lib -lfftw3f_mpi -lfftw3f
FFTW_ROOT = /home/apps/fftw/3.3.3
FFTW_LIBS = -L$(FFTW_ROOT)/lib64 -lfftw3_mpi -lfftw3_omp -lfftw3

CXXFLAGS += -I$(FFTW_ROOT)/include
CFLAGS += -I$(FFTW_ROOT)/include
CLIBS +=  -L$(PM_ROOT) -lpm $(FFTW_LIBS)

$(PM_ROOT)/libpm.a: $(PM_ROOT)/param_fdps.h
	cd $(PM_ROOT) && $(MAKE) allclean libpm.a

treepm.o: treepm.hpp cosmology.hpp run_param.hpp

treepm:  $(TREEPM_DEP) $(PM_ROOT)/libpm.a
	$(CXX) $(CXXFLAGS) -o $@ $(TREEPM_OBJ) $(CLIBS) -lm 

makeglass: $(MAKEGLASS_DEP) $(PM_ROOT)/libpm.a
	$(CXX) $(CXXFLAGS) -o $@ $(MAKEGLASS_OBJ) $(CLIBS) -lm

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp -o $*.o

clean: 
	-rm -rf *.o

distclean: clean
	-rm -rf *~



