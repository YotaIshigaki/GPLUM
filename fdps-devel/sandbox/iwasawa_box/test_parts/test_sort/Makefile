#PS_PATH = -I../../../../src/
#PS_PATH = -I../../my_trunk/src/
PS_PATH = -I../../../../../galaxy-challenge/src
GTEST_PATH = -I/home/masaki/project/git/googletest/googletest/include/
GTEST_LIB = -L/home/masaki/project/git/googletest/googletest/lib -lgtest


CXX = mpicxx
CXXFLAGS = -O3
CXXFLAGS += -std=c++14
CXXFLAGS += -Wall
CXXFLAGS += --pedantic-errors
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DUSE_128BIT_KEY
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY=1
#CXXFLAGS += -DGOOGLE_TEST

#SRC = test.cc
#PROGRAM = test.out
SRC = test_sort.cpp
PROGRAM = test_sort.out

$(PROGRAM):$(SRC)
	$(CXX) $(PS_PATH) $(GTEST_PATH)  $(CXXFLAGS) -o $@ $< $(GTEST_LIB)

test_sort.s:$(SRC) 
	$(CXX)   -S $(PS_PATH) $(GTEST_PATH)  $(CXXFLAGS) -o $@ $< $(GTEST_LIB)

clean:
	rm $(PROGRAM) 
