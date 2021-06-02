PS_PATH = -I../../my_trunk/src/

CXX = mpicxx
CXXFLAGS = -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -Wall
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

SRC = test.cc
PROGRAM = test.out

$(PROGRAM):$(SRC) 
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm $(PROGRAM) 
