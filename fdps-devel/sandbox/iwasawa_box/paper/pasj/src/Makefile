#PS_PATH = -I../../../src/
#PS_PATH = -I../src_sc/
PS_PATH = -I../../../../../src/

CXX = time mpiFCCpx
#CXX = time FCCpx
CXXFLAGS = -Kfast
#CXXFLAGS += -NRtrap
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -x32
CXXFLAGS += -Xg
#CXXFLAGS += -Nfjcex
#CXXFLAGS += -Krestp=all
CXXFLAGS += -DFAST_WALK_K
#CXXFLAGS = -Krestp=all
CXXFLAGS += -DUSEPHANTOMGRAPE_K
CXXFLAGS += -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE

#CXX = mpicxx
#CXX = CC
#CCC = g++
#CXXFLAGS = -O3 
#CXXFLAGS += -march=core-avx2
#CXXFLAGS += -ffast-math -funroll-loops
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DUSEPHANTOMGRAPE_X86
#CXXFLAGS += -std=c++11

MULEXP = -DQUADRUPOLE

SRC = nbody_disk.cc class.hpp

#all:nbody_disk.out nbody_disk_simple.out
all:nbody_disk_simple.out nbody_disk_simple_K.out
#all:nbody_disk_simple_K.out
#all:nbody_disk_simple_K2.out
#all:nbody_disk_simple_K3_1.out
#all:nbody_disk_simple_K3_2.out
#all:nbody_disk_simple_K3_4.out
#all:nbody_disk_simple_K3_8.out
#all:nbody_disk_simple.out nbody_disk_simple_K.out nbody_disk_simple_K3_1.out nbody_disk_simple_K3_2.out nbody_disk_simple_K3_4.out nbody_disk_simple_K3_8.out nbody_disk_simple_K3_16.out nbody_disk_simple_K3_32.out nbody_disk_simple_K3_64.out nbody_disk_simple_K3_128.out nbody_disk_simple_K4_1.out nbody_disk_simple_K4_2.out nbody_disk_simple_K4_4.out nbody_disk_simple_K4_8.out nbody_disk_simple_K4_16.out nbody_disk_simple_K4_32.out nbody_disk_simple_K4_64.out nbody_disk_simple_K4_128.out 
#all:nbody_disk_simple_K3_0.out nbody_disk_simple_K4_0.out
#all:nbody_disk_simple_K6_0.out nbody_disk_simple_K6_1.out nbody_disk_simple_K6_2.out nbody_disk_simple_K6_4.out nbody_disk_simple_K6_8.out nbody_disk_simple_K6_16.out nbody_disk_simple_K6_32.out nbody_disk_simple_K6_64.out nbody_disk_simple_K6_128.out

nbody_disk.out:$(SRC)
	$(CXX) $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K2.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K2 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_0.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_0 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_1.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_1 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_2.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_2 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_4.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_4 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_8.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_8 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_16.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_16 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_32.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_32 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_64.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_64 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K3_128.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K3 -DNPTCL_A2A_128 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_0.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_0 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_1.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_1 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_2.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_2 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_4.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_4 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_8.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_8 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_16.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_16 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_32.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_32 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_64.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_64 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K4_128.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K4 -DNPTCL_A2A_128 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<


nbody_disk_simple_K6_0.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_0 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_1.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_1 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_2.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_2 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_4.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_4 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_8.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_8 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_16.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_16 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_32.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_32 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_64.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_64 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

nbody_disk_simple_K6_128.out:$(SRC)
	$(CXX) -DSIMPLE_VERSION -DFAST_ALL_TO_ALL_FOR_K6 -DNPTCL_A2A_128 $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<


clean:
	rm nbody_disk.out

