PS_PATH = -I../../../../src/
#PS_PATH = -I../../my_trunk/src/
#PS_PATH = -I../../../../../../project/git/FDPS/src/
#PS_PATH = -I../../develop/src_remove_a2a_in_exLET/
#PS_PATH = -I../../develop/src2/

CXX = mpicxx
#CXX = g++
#CXXFLAGS = -O1
CXXFLAGS += -std=c++11
CXXFLAGS += -Wall
#CXXFLAGS += -march=core-avx2
#CXXFLAGS += -ffast-math -funroll-loops -Wall 
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY=1
CXXFLAGS += -DSHRINK_SIZE_OF_SPJ_SORTED
#CXXFLAGS += -DUSE_ALLGATHER_EXLET

SRC = test.cc
PROGRAM = test.out

$(PROGRAM):$(SRC) 
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm $(PROGRAM) 

