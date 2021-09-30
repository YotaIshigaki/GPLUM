FFTW_LOC  = /opt/local
FFTW_INC  = -I$(FFTW_LOC)/include
FFTW_LIB  = -L$(FFTW_LOC)/lib -lfftw3 
#FFTW_LIB  = -L$(FFTW_LOC)/lib -lfftw3_mpi -lfftw3_omp -lfftw3

FDPS_LOC  = ./
PMM_ROOT  = $(FDPS_LOC)/src/particle_mesh_multipole
DCT_ROOT  = $(FDPS_LOC)/src/direct_sum
EWALD_ROOT = $(FDPS_LOC)/src/ewald
FDPS_INC  = -I$(FDPS_LOC)/src -I$(PMM_ROOT) -I$(DCT_ROOT) -I$(EWALD_ROOT)
FDPS_LIB  = 
LDFLAGS   = $(FDPS_LIB) $(FFTW_LIB)

#CXX = g++-mp-8
CXX = mpicxx-openmpi-gcc8
#CXXFLAGS = -std=c++11 -Wall -Wextra -Wcast-align -Wconversion -O0 -fvar-tracking -ggdb3 -rdynamic
#CXXFLAGS = -std=c++11 -Wall -Wextra -Wcast-align -Wconversion -O0 -fvar-tracking -ggdb3
CXXFLAGS = -std=c++11 -O3 -ffast-math -funroll-loops
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_USE_EPJ_ONLY_TO_EVAL_MM_IN_PMM
CXXFLAGS += $(FDPS_INC) $(FFTW_INC)

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

#use_phantom_grape_x86 = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

ifeq ($(use_phantom_grape_x86),yes)
PG_ROOT = $(FDPS_LOC)/src/phantom_grape_x86/G5/newton/libpg5
CXXFLAGS += -I$(PG_ROOT) -DENABLE_PHANTOM_GRAPE_X86
LDFLAGS += -L$(PG_ROOT) -lpg5
PG_BUILD = cd $(PG_ROOT) && $(MAKE) distclean libpg5.a
PG_CLEAN = cd $(PG_ROOT) && $(MAKE) distclean
else
PG_BUILD =
PG_CLEAN = 
endif

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = pmmm.out

.PHONY:	clean all

all: $(CPPOBJS) $(CPPHDRS)
	$(PG_BUILD)
	$(CXX) $(CXXFLAGS) $(CPPOBJS) -o $(PROGRAM) $(LDFLAGS)

%.o: %.cpp $(CPPHDRS)
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -f $(CPPOBJS)

distclean: clean
	$(PG_CLEAN)
	rm -f $(PROGRAM)
	rm -rf result
	rm -f *.txt

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
