#PS_PATH = -I../../../src/
#PS_PATH = -I../src_sc/
PS_PATH = -I../../../src/

CXX = time mpiFCCpx
#CXX = time FCCpx
#CXXFLAGS = -Kfast,openmp -NRtrap
CXXFLAGS = -Kfast
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
#CXXFLAGS += -Nfjcex
CXXFLAGS += -Krestp=all
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
MULEXP = -DQUADRUPOLE

SRC = nbody_disk_for_sc.cc

nbody_disk_for_sc.out:$(SRC)
	$(CXX) -DUSE_EMA $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_plummer_for_sc.out:$(SRC)
	$(CXX) -DUSE_EMA -DUSE_PLUMMER $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_for_sc_db_a2a.out:$(SRC)
	$(CXX) -DUSE_EMA -DDB_PRINT_ALL_TO_ALL $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_for_sc_snapshot.out:$(SRC)
	$(CXX) -DUSE_EMA -DGET_SNAPSHOT $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_for_sc_k8.out:$(SRC)
	$(CXX) -DUSE_K8 -DUSE_EMA $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_for_sc_kfull.out:$(SRC)
	$(CXX) -DUSE_KFULL -DUSE_EMA $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<


clean:
	rm $(PROGRAM) 

