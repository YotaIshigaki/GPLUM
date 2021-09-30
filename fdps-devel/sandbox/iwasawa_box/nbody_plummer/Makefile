FDPS_HOME = ../../../src/

CXX = time mpiFCCpx
CXXFLAGS = -Kfast
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
CXXFLAGS += -Krestp=all
CXXFLAGS += -DFAST_WALK_K
CXXFLAGS += -Nfjcex


#CXX = time mpicxx
#CUDA_HOME = /usr/local/cuda
#NVCC = time $(CUDA_HOME)/bin/nvcc -Xcompiler="-Wall -O2 -march=core2 -mavx"
#INC = -I$(CUDA_HOME)/samples/common/inc/ -I$(FDPS_HOME)
INC = -I$(FDPS_HOME)
#CXXFLAGS = -O2
#CXXFLAGS+= -Wall
#CXXFLAGS+= -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL 
#LIB = -lgomp
#COMP = --compiler-bindir=mpicxx --compiler-options=-fopenmp


nbody_plummer.out:nbody_plummer.cc
	$(CXX) -DMONOPOLE $(CXXFLAGS) $(INC) -o $@ $^ 

#nbody_plummer.out:nbody_plummer.o force.o
#	$(CXX) $(CXXFLAGS) $(INC) -L $(CUDA_HOME)/lib64 -o $@ $^ -lcudart

test.out:test.o force.o
	$(NVCC) $(INC) $(CXXFLAGS) $(COMP) -o test.out test.o force.o -lgomp

test.o:test.cc
	$(NVCC) -DUSE_GPU $(INC) $(CXXFLAGS) $(COMP) -c -o test.o test.cc


#force.o:force.cu
#force.o:force_kn.cu
#force.o:force_float.cu
#	$(NVCC) -DUSE_GPU_CUDA $(INC) -c -o $@ $<

nbody_plummer.o:nbody_plummer.cc
	$(CXX) -DUSE_GPU_CUDA -DMONOPOLE $(CXXFLAGS) $(INC) -c -o $@ $<

force.cubin:force_float.cu
	$(NVCC) -DUSE_GPU_CUDA --cubin $(INC) -o $@ $<

force.sass:force.cubin
	$(CUDA_HOME)/bin/cuobjdump -sass $< > $@

clean:
	rm *.o *.out

