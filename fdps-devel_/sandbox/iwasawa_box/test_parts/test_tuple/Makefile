PS_PATH = -I../../../../src/

CXX = mpicxx
CXXFLAGS = -O3 -Wall -march=core-avx2
CXXFLAGS = -std=c++11
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DTEST_VARIADIC_TEMPLATE

SRC = test_set_particle_local_tree.cc
PROGRAM = test_set_particle_local_tree.out

$(PROGRAM):$(SRC)
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm $(PROGRAM) 
