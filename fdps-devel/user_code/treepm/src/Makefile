#CC := mpicc
#CXX := mpicxx
CC := cc
CXX := CC
CFLAGS = -O3 -ffast-math -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS := $(CFLAGS)

all: treepm makeglass santabarbara

TREEPM_OBJ = treepm.o make_directory.o timing.o
TREEPM_DEP := $(TREEPM_OBJ) treepm.hpp cosmology.hpp run_param.hpp

MAKEGLASS_OBJ = makeglass.o make_directory.o timing.o
MAKEGLASS_DEP := $(MAKEGLASS_OBJ) treepm.hpp cosmology.hpp run_param.hpp

SANTABARBARA_OBJ = santabarbara.o make_directory.o timing.o
SANTABARBARA_DEP := $(SANTABARBARA_OBJ) treepm.hpp cosmology.hpp run_param.hpp

FDPS_ROOT=/home/iwaswams/project/fdps
PM_ROOT=$(FDPS_ROOT)/src/particle_mesh

CXXFLAGS += -I$(FDPS_ROOT)/src
CXXFLAGS += -I$(PM_ROOT)

#FFTW_ROOT = /usr/local/fftw-3
FFTW_LIBS = -L$(FFTW_ROOT)/lib -lfftw3f_mpi -lfftw3f

CXXFLAGS += -I$(FFTW_ROOT)/include
CFLAGS += -I$(FFTW_ROOT)/include
CFLAGS += -I$(FDPS_ROOT)/src
CLIBS +=  -L$(PM_ROOT) -lpm $(FFTW_LIBS)

$(PM_ROOT)/libpm.a: $(PM_ROOT)/param_fdps.h
	cd $(PM_ROOT) && $(MAKE) allclean libpm.a

treepm.o: treepm.hpp cosmology.hpp run_param.hpp

treepm:  $(TREEPM_DEP) $(PM_ROOT)/libpm.a
	$(CXX) $(CXXFLAGS) -o $@ $(TREEPM_OBJ) $(CLIBS) -lm 

makeglass: $(MAKEGLASS_DEP) $(PM_ROOT)/libpm.a
	$(CXX) $(CXXFLAGS) -o $@ $(MAKEGLASS_OBJ) $(CLIBS) -lm

santabarbara: $(SANTABARBARA_DEP) $(PM_ROOT)/libpm.a
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $@ $(SANTABARBARA_OBJ) $(CLIBS) -lm

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp -o $*.o

clean: 
	-rm -rf *.o

distclean: clean
	-rm -rf *~



