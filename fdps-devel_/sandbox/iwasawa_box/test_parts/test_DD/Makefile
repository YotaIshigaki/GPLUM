#PS_PATH = -I../../develop/src_remove_a2a_in_exLET/
PS_PATH = -I../../../../src/

#GTEST_PATH = -I/home/masaki/local/include/
#GTEST_LIB = -L/home/masaki/local/lib -lgtest

GTEST_PATH = -I/home/masaki/project/git/googletest/googletest/include/
GTEST_LIB = -L/home/masaki/project/git/googletest/googletest/lib -lgtest

CXX = mpicxx
CXXFLAGS += -O3 -Wall
CXXFLAGS += -std=c++11
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

SRC = test_DD.cpp
PROGRAM = test_DD.out

$(PROGRAM):$(SRC)
	$(CXX) $(PS_PATH) $(GTEST_PATH) $(CXXFLAGS) -DSANITY_CHECK_REALLOCATABLE_ARRAY=1 -o $@ $< $(GTEST_LIB)

clean:
	rm $(PROGRAM) *~

