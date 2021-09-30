PS_PATH = -I../../../src/

#CXX = mpicxx
CXX = CC
CXXFLAGS = -O3 -march=core-avx2 #-Wall
#CXXFLAGS = -O3 -ffast-math -funroll-loops -Wall 
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DUSEXC30

#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

CXXFLAGS += -DUSEPHANTOMGRAPE
#MULEXP = -DMONOPOLE
MULEXP = -DQUADRUPOLE
MULEXP += -DUSE_AT_PHANTOM

SRC = test_nbody.cpp
#SRC = test_nbody_dev.cc
#SRC = test_sph.cc
#SRC = test_sph_consth.cc
#PROGRAM = $(SRC:%.cc=%.out)
#PROGRAM = test_nbody_flatmpi.out

test_nbody.out:$(SRC) phantomquad_x86.hpp Makefile
	$(CXX) -DUSE_EMA $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

test_nbody_a2a.out:$(SRC)
	$(CXX) -DUSE_EMA -DFAST_ALL_TO_ALL_FOR_K $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm -f *~ test.o* core
distclean:
	make clean
	rm -f *.out

