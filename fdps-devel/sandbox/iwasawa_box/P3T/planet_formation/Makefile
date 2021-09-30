PS_PATH = -I../../../../../fdps/src/
#PS_PATH = -I../src/

#CXX = time mpiFCCpx
#CXXFLAGS = -Kfast
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -x32
#CXXFLAGS += -Xg
#CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
#CXXFLAGS += -DFAST_WALK_K
#CXXFLAGS += -DINTRINSIC_K
#CXXFLAGS += -std=c++11
#CXXFLAGS += -NRtrap
#CXXFLAGS += -Nfjcex
#CXXFLAGS += -Krestp=all

CXX = time mpicxx -Wall
#CXX = time CC
#CCC = g++
CXXFLAGS = -O3 
CXXFLAGS += -Wall
CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops

CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -std=c++11

#CXXFLAGS += -DINTRINSIC_X86
#CXXFLAGS += -DDEBUG_PRINT_PLANET

#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

#CXXFLAGS += -DUSE_STD_SORT

#CXXFLAGS += -DMERGE
#CXXFLAGS += -DFORDEBUG

#CXXFLAGS += -DFREE_MALLOC_TREE 
#CXXFLAGS += -DPARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
#CXXFLAGS += -DPARTICLE_SIMULATOR_EXCHANGE_LET_ALL
#CXXFLAGS += -DKEP_EQ_ACC_10 

SRC = main_kepler.cc hard.hpp class.hpp force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp

#all:p3t_kepler_quad_merge_fix.out
all:p3t_kepler_quad_merge_fix_cont.out
#all:p3t_kepler_quad_merge_phantom_fix.out
#all:p3t_kepler_quad_nomerge_phantom_fix.out
#all:p3t_kepler_quad_nomerge_phantom_cont.out


p3t_kepler_quad_merge_fix.out:$(SRC) Makefile
	$(CXX) -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

p3t_kepler_quad_merge_fix_cont.out:$(SRC) Makefile
	$(CXX) -DUSE_QUAD -DDIV_FIX -DMERGE -DREAD_FILE $(PS_PATH) $(CXXFLAGS) -o $@ $<

p3t_kepler_quad_merge_phantom_fix.out:$(SRC) Makefile
	$(CXX) -DINTRINSIC_X86 -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

p3t_kepler_quad_nomerge_phantom_fix.out:$(SRC) Makefile
	$(CXX) -DINTRINSIC_X86 -DUSE_QUAD -DDIV_FIX $(PS_PATH) $(CXXFLAGS) -o $@ $<

p3t_kepler_quad_nomerge_phantom_cont.out:$(SRC) Makefile
	$(CXX) -DINTRINSIC_X86 -DUSE_QUAD -DREAD_FILE $(PS_PATH) $(CXXFLAGS) -o $@ $<
kepler.out:kepler.cc
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm p3t_kepler.out
cleanall:
	rm p3t_kepler.out *.hpp~ *.cc~ *.h~
