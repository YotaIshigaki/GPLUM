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
NVCC = time $(CUDA_HOME)/bin/nvcc 
INC = -I$(CUDA_HOME)/samples/common/inc/ -I$(FDPS_HOME)
CXXFLAGS = -O3
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL 
LIB = -lgomp
#COMP = --compiler-bindir=g++ --compiler-options=-fopenmp
COMP = --compiler-bindir mpicxx --compiler-options=-fopenmp
#COMP = --compiler-bindir mpiicpc --compiler-options=-openmp
#CXXFLAGS += --compiler-bindir mpiicpc


#test.out:test.o force.o
#	$(CXX) -DGPU $(INC) $(CXXFLAGS) -fopenmp -o $@ test.o force.o $(LIB)

test.out:test.o force.o
	$(NVCC) $(INC) $(CXXFLAGS) $(COMP) -o test.out test.o force.o -lgomp

test.o:test.cc
	$(NVCC) -DGPU $(INC) $(CXXFLAGS) $(COMP) -c -o test.o test.cc

#force.o:force.cu
#	$(NVCC) -DGPU $(INC) -c -o force.o force.cu
force.o:force_float.cu
	$(NVCC) -DGPU $(INC) -c -o $@ $<

#test.o:test.cc
#	$(CXX) $(INC) $(CXXFLAGS) $^ -c -o $@

#force.o:force.cu
#	$(NVCC) -DGPU $(INC) $^ -c -o $@

#force.o:force.cc Makefile
#	$(CXX) $(INC) -c -o $@ $<


clean:
	rm *.o *.out

