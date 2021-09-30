#PS_PATH = -I../../../src/
#PS_PATH = -I../src_sc/
PS_PATH = ../../../src/
INC  = -I$(PS_PATH)

#CXX = time mpiFCCpx
#CXX = time FCCpx
#CXXFLAGS = -Kfast,openmp -NRtrap
#CXXFLAGS = -Kfast
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -x32
#CXXFLAGS += -Xg
#CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
#CXXFLAGS += -Krestp=all
#CXXFLAGS += -DFAST_WALK_K
#CXXFLAGS += -Nfjcex
#CXXFLAGS = -Kfast,openmp,restp=all

CXX = mpicxx
#CXX = mpiicpc
CXX = CC
#CCC = g++
CXXFLAGS = -O2
#CXXFLAGS += -march=core-avx2
CXXFLAGS += -Wall 
#CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CXXFLAGS += -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

#CXXFLAGS += -DUSEPHANTOMGRAPE
#MULEXP = -DQUADRUPOLE
MULEXP = MONOPOLE

use_gpu_cuda = yes

ifeq ($(use_gpu_cuda),yes)
#NVCC = time nvcc
#CUDA_HOME = /gwfefs/opt/x86_64/cuda/7.5
CUDA_HOME = /usr/local/cuda
INC  += -I$(CUDA_HOME)/samples/common/inc/
DEF  = -DUSE_GPU_CUDA -D$(MULEXP)
#COMP = --compiler-bindir=g++ --compiler-options=-fopenmp
#COMP = --compiler-bindir=mpicxx --compiler-options=-fopenmp --compiler-options=-O3
#COMP = -ccbin=mpicxx -Xcompiler='-O3 -ffast-math -funroll-loops -fopenmp'
#COMP = -ccbin=mpicxx -Xcompiler='-O2 -ffast-math -funroll-loops'
NVCC = time $(CUDA_HOME)/bin/nvcc -Xcompiler="-Wall -O2 -march=core2 -mavx"
#COMP = --compiler-bindir mpicxx
endif

SRC = nbody_disk.cc

all:nbody_disk_gpu.out

nbody_disk_gpu.out:nbody_disk_gpu.o force.o
	$(CXX) $(DEF) $(CXXFLAGS) $(INC) -L $(CUDA_HOME)/lib64 -o $@ $^ -lcudart

nbody_disk_gpu.o:nbody_disk.cc
	$(CXX) $(DEF) $(CXXFLAGS) $(INC) -L $(CUDA_HOME)/lib64 -c -o $@ $^ -lcudart

#force.o:force_float.cu force.hpp
force.o:force_kn.cu force.hpp
	$(NVCC) $(INC) $(DEF) $(COMP) -c -o $@ $<

nbody_disk.out:$(SRC)
	$(CXX) $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_kfull.out:$(SRC)
	$(CXX) -DUSE_KFULL $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_k2.out:$(SRC)
	$(CXX) -DUSE_K2 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_snapshot.out:$(SRC)
	$(CXX) -DGET_SNAPSHOT $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm nbody_disk_gpu.out *.o
