FDPS_HOME = ../../../src/

#CXX = time mpiFCCpx
#CXXFLAGS = -Kfast,openmp -NRtrap
#CXXFLAGS = -Kfast
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -x32
#CXXFLAGS += -Xg
#CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
#CXXFLAGS += -Nfjcex
#CXXFLAGS += -Krestp=all
#CXXFLAGS = -Kfast,openmp,restp=all

#CXX = time mpicxx -Wall
#CXX = CC
#CCC = g++
#CXXFLAGS = -O3 -march=core-avx2
#CXXFLAGS += -ffast-math -funroll-loops
#CXXFLAGS += -std=c++11
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CXXFLAGS += -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

#CXXFLAGS += -DUSEPHANTOMGRAPE
#CXXFLAGS += -DFORDEBUG

#CXX = time mpiicpc
CXX = time mpicxx
#NVCC = time nvcc
#CUDA_HOME = /gwfefs/opt/x86_64/cuda/7.5
CUDA_HOME = /usr/local/cuda
#NVCC = time $(CUDA_HOME)/bin/nvcc -arch=sm_52 -Xcompiler="-Wall -O2 -march=core2 -mavx"
#NVCC = time $(CUDA_HOME)/bin/nvcc -arch=sm_35 -Xcompiler="-Wall -O2 -march=core2 -mavx"
NVCC = time $(CUDA_HOME)/bin/nvcc -arch=sm_30 -Xcompiler="-Wall -O2 -march=core2"
#NVCC = time $(CUDA_HOME)/bin/nvcc -arch=sm_35 -Xcompiler="-Wall -O2 -march=core-avx2"
INC = -I$(CUDA_HOME)/samples/common/inc/ -I$(FDPS_HOME)
CXXFLAGS = -O2
CXXFLAGS+= -Wall
CXXFLAGS+= -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL 
LIB = -lgomp
#COMP = --compiler-bindir=g++ --compiler-options=-fopenmp
COMP = --compiler-bindir=mpicxx --compiler-options=-fopenmp
#COMP = --compiler-bindir mpiicpc --compiler-options=-openmp
#CXXFLAGS += --compiler-bindir mpiicpc


nbody_plummer.out:nbody_plummer.o force.o
	$(NVCC) $(COMP) -o $@ $^ -lgomp
#	$(CXX) $(CXXFLAGS) -L $(CUDA_HOME)/lib64 -o $@ $^ -lcudart

test.out:test.o force.o
	$(NVCC) $(INC) $(CXXFLAGS) $(COMP) -o test.out test.o force.o -lgomp

test.o:test.cc
	$(NVCC) -DGPU $(INC) $(CXXFLAGS) $(COMP) -c -o test.o test.cc

#force.o:force_float.cu
#force.o:force.cu
force.o:force_kn.cu
	$(NVCC) -DGPU $(INC) -c -o $@ $<

nbody_plummer.o:nbody_plummer.cc
	$(CXX) -DGPU -DMONOPOLE $(CXXFLAGS) $(INC) -c -o $@ $<

force.cubin:force_float.cu
	$(NVCC) -DGPU --cubin $(INC) -o $@ $<

force.sass:force.cubin
	$(CUDA_HOME)/bin/cuobjdump -sass $< > $@

clean:
	rm *.o *.out

