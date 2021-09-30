# (i) PATH setting
FFTW_LOC  = /opt/klocal/fftw-3.3.6pl2-fujitsu
FFTW_INC  = -I$(FFTW_LOC)/include
#FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3 
FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_omp -lfftw3
#FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_mpi -lfftw3_omp -lfftw3

FDPS_LOC  = ./
PMM_ROOT  = $(FDPS_LOC)/src/particle_mesh_multipole
DCT_ROOT  = $(FDPS_LOC)/src/direct_sum
EWALD_ROOT = $(FDPS_LOC)/src/ewald
FDPS_INC  = -I$(FDPS_LOC)/src -I$(PMM_ROOT) -I$(DCT_ROOT) -I$(EWALD_ROOT)
FDPS_LIB  = 
LDFLAGS   = $(FDPS_LIB) $(FFTW_LIB)

# (ii) Compiler and compile options
CXX=time -p FCCpx
#CXX=time -p mpiFCCpx
#CXXFLAGS = -std=c++11 -O0 -g -NRtrap -Nquickdbg=subchk -Nquickdbg=heapchk -Nquickdbg=inf_detail
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DPARTICLE_SIMULATOR_USE_SFINAE_OLD_IMPL
CXXFLAGS += -DPARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
CXXFLAGS += -DPARTICLE_SIMULATOR_USE_ERFC_IN_MATH_H
CXXFLAGS += -DCODE_NUM=2
#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY=2
CXXFLAGS += $(FDPS_INC) $(FFTW_INC)

# (iii) Phantom-GRAPE setting
#use_phantom_grape_x86 = yes
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

# (iv) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/PMMM/develop/debug/
#WORKDIR = /data/ra000008/namekata/PMMM/develop/nproc32/
#WORKDIR = /data/ra000008/namekata/PMMM/develop/nproc128_w_mpifft/
#WORKDIR = /data/ra000008/namekata/PMMM/develop/nproc128_wo_mpifft/
#WORKDIR = /data/ra000008/namekata/PMMM/develop/nproc256_w_mpifft/
#WORKDIR = /data/ra000008/namekata/PMMM/develop/nproc256_wo_mpifft/


TARGET = pmmm.out

all: main.o
	$(PG_BUILD)
	$(CXX) $(CXXFLAGS) $< -o $(TARGET) $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv $(TARGET) $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

main.o: main.cpp user_defined.hpp
	$(CXX) -c main.cpp $(CXXFLAGS)

clean:
	rm -f $(CPPOBJS)

distclean: clean
	$(PG_CLEAN)
	rm -f $(PROGRAM)
	rm -rf result
	rm -f *.txt

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
