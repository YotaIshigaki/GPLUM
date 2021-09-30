PS_PATH = -I../../../../src/

CXX = mpicxx
#CXX = CC
#CCC = g++
CXXFLAGS = -O3
CXXFLAGS += -std=c++11
#CXXFLAGS += -ffast-math -funroll-loops -Wall 
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

SRC = test_exptcl.cpp
PROGRAM = test_exptcl.out

$(PROGRAM):$(SRC) 
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

test_exptcl2.out:test_exptcl2.cpp
	$(CXX) $(PS_PATH) $(CXXFLAGS) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp -o $@ $<

clean:
	rm $(PROGRAM) 
