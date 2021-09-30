PS_PATH = -I../../../src/

CXX = time mpiFCCpx
#CXX = time FCCpx
#CXXFLAGS = -Kfast,openmp -NRtrap
#CXXFLAGS = -Kfast,openmp
CXXFLAGS = -Kfast
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -x32
CXXFLAGS += -Xg
#CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
#CXXFLAGS += -Nfjcex
#CXXFLAGS += -Krestp=all
#CXXFLAGS = -Kfast,openmp,restp=all

#CXX = mpicxx
#CXX = CC
#CCC = g++
#CXXFLAGS = -O3 -Wall -march=core-avx2
#CXXFLAGS = -O3 -ffast-math -funroll-loops -Wall 
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

CXXFLAGS += -DUSEPHANTOMGRAPE
#MULEXP = -DMONOPOLE
MULEXP = -DQUADRUPOLE

SRC = test_nbody.cc
#SRC = test_nbody_dev.cc
#SRC = test_sph.cc
#SRC = test_sph_consth.cc
#PROGRAM = $(SRC:%.cc=%.out)
#PROGRAM = test_nbody_flatmpi.out

test_nbody.out:$(SRC)
	$(CXX) -DUSE_EMA -DFAST_ALL_TO_ALL_FOR_K $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

test_nbody_orgDD.out:$(SRC)
	$(CXX) -DORIGINAL_DD -DUSE_EMA -DFAST_ALL_TO_ALL_FOR_K $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

#test_nbody_a2a.out:$(SRC)
#	$(CXX) -DUSE_EMA -DFAST_ALL_TO_ALL_FOR_K $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm $(PROGRAM) 

