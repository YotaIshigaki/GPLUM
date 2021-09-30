PS_PATH = -I../../../../src/

GTEST_PATH = -I/home/masaki/project/git/googletest/googletest/include/
GTEST_LIB = -L/home/masaki/project/git/googletest/googletest/lib -lgtest

CXX = mpicxx
CXXFLAGS += -O3 -Wall
CXXFLAGS += -std=c++11
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

SRC = test_key.cpp
PROGRAM = test_key.out

$(PROGRAM):$(SRC)
	$(CXX) $(PS_PATH) $(GTEST_PATH) $(CXXFLAGS) -DSANITY_CHECK_REALLOCATABLE_ARRAY=1 -o $@ $< $(GTEST_LIB)

clean:
	rm $(PROGRAM) *~

