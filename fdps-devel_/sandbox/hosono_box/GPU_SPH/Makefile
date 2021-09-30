FDPS_HOME = ../../../src/

#CXX = time mpiicpc
CXX = time mpic++
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

PROGRAM = sph.out
CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CUOBJS  = $(patsubst %.cu, %.o, $(wildcard *.cu))
CPPHDRS = $(wildcard *.h)

all:	$(CPPOBJS) $(CUOBJS) $(CPPHDRS)
	@echo "Linking object files..."
#	@$(NVCC) $(CXXFLAGS) $(INC) $(CPPOBJS) $(CUOBJS) -o $(PROGRAM) -lgomp
	@$(CXX) $(CXXFLAGS) $(INC) $(CPPOBJS) $(CUOBJS) -L $(CUDA_HOME)/lib64 -o $(PROGRAM) -lcudart -lgomp
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
#	@$(NVCC) -DGPU $(INC) -c $< $(CXXFLAGS) $(COMP)
	@$(CXX) -DGPU $(CXXFLAGS) $(INC) -c -o $@ $<
	@echo "[CPP $< OK]"
	@echo ""

%.o:	%.cu $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(NVCC) -DGPU $(INC) -c $< -o $@
	@echo "[CUDA $< OK]"
	@echo ""

clean:
	rm *.o *.out

