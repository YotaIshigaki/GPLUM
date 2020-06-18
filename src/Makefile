PS_PATH = ../FDPS-5.0d2/src/ 
INC = -I$(PS_PATH)

####################
###   Compiler   ###
####################
#CC = clang++ -std=c++11
CC = mpic++ -std=c++11
#CC = time CC

CFLAGS = -O2
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops

########################
###   MPI & OpenMP   ###
########################
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

################
###   SIMD   ###
################
CFLAGS_M += -march=core-avx2
CFLAGS_G += -march=core-avx2
#CFLAGS_G += -mtune=skylake-avx512
#CFLAGS_G += -march=skylake-avx512
#CFLAGS_G += -mavx2
#CFLAGS_G += -mavx512f -mavx512dq

#########################
###   Other Options   ###
#########################
CFLAGS += -DUSE_INDIVIDUAL_CUTOFF
CFLAGS += -DUSE_QUAD
CFLAGS += -DCALC_EPEP_64bit
CFLAGS += -DCALC_EPSP_64bit
CFLAGS += -DUSE_ALLGATHER_EXLET
CFLAGS += -DUSE_RE_SEARCH_NEIGHBOR
#CFLAGS += -DISOTROPIC
#CFLAGS += -DGAS_DRAG
#CFLAGS += -DINDIRECT_TERM
#CFLAGS += -DINTEGRATE_6TH_SUN
#CFLAGS += -DAARSETH
CFLAGS += -DCORRECT_NEIGHBOR
#CFLAGS += -DABSORB_PTCL

CFLAGS += -DCOLLISION
CFLAGS += -DKOMINAMI
#CFLAGS += -DCHAMBERS
CFLAGS += -DMERGE_BINARY

CFLAGS += -DOUTPUT_DETAIL
CFLAGS += -DCALC_WTIME
#CFLAGS += -DTEST_PTCL
CFLAGS += -DMONAR
#CFLAGS += -pg

SRC = main_p3t.cpp
PROGRAM = gplum.out
HEADER = mathfunc.h kepler.h energy.h particle_ep.h particle_fp.h disk.h gravity.h gravity_kernel.hpp collisionA.h collisionB.h hermite.h hard.h read.h time.h func.h

$(PROGRAM): main_p3t.o gravity_kernel.o
	$(CC) -MMD $(INC) $(CFLAGS) $(CFLAGS_M) -o $@ main_p3t.o gravity_kernel.o

main_p3t.o: main_p3t.cpp $(HEADER)
	$(CC) -MMD $(INC) $(CFLAGS) $(CFLAGS_M) -c $<

gravity_kernel.o: gravity_kernel.cpp phantomquad_for_p3t_x86.hpp particle_ep.h 
	$(CC) -MMD $(INC) $(CFLAGS) $(CFLAGS_G) -c $<


all: clean $(PROGRAM)

clean:
	rm -f *.out *.o *.d *~

