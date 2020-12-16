CCC = mpicxx
#CCC = mpiFCCpx
#PS = ../../develop/src_remove_a2a_in_exLET/
#PS = ../../develop/src2/
PS = ../../../../src/

#CFLAGS = -O3 -Wall -DNDEBUG
CFLAGS = -O3 -ffast-math -funroll-loops -Wall
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

GTEST_PATH = -I/home/masaki/project/git/googletest/googletest/include/
GTEST_LIB = -L/home/masaki/project/git/googletest/googletest/lib -lgtest


SRC = test_comm.cc

PROGRAM = $(SRC:%.cc=%.out)

test_comm_info.out:test_comm_info.cpp
	$(CCC) -std=c++11 -I$(PS) $(GTEST_PATH) $(CFLAGS) -o $@ $< $(GTEST_LIB)

$(PROGRAM):$(SRC)
	$(CCC) -std=c++11 -I$(PS) $(CFLAGS) -o $@ $<

clean:
	rm *.s *.o $(PROGRAM)
