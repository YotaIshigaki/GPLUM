FDPS_HOME = ../../../../src/

CXX = time mpicxx
INC = -I$(FDPS_HOME)
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DSET_CALLER_INFO_REALLOCATABLE_ARRAY
CXXFLAGS += -O3 -Wall

test_reallocatable_array.out:test_reallocatable_array.cpp
	$(CXX) $(INC) $(CXXFLAGS) -o $@ $<
