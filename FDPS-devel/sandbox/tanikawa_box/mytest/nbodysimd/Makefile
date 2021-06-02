DIR_FDPS = ../../../../src/
PS_PATH = -I$(DIR_FDPS)

#compiler = gcc
compiler = intel
use_mic = yes
use_mpi = yes
#use_intrinsics = yes

ifeq ($(compiler),gcc)
CC = g++
CFLAGS = -O3 -ffast-math -funroll-loops #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
ifeq ($(use_mpi),yes)
CC = mpicxx
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
ifeq ($(use_intrinsics),yes)
CFLAGS  += -DENABLE_INTRINSICS -mavx
endif
endif

ifeq ($(compiler),intel)
CC = icpc
CFLAGS = -O3 #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -openmp
ifeq ($(use_mic),yes)
CFLAGS += -mmic
endif
ifeq ($(use_mpi),yes)
CC = mpiicpc
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
ifeq ($(use_intrinsics),yes)
CFLAGS  += -DENABLE_INTRINSICS
endif
endif

ifeq ($(use_phantom_grape_x86),yes)
PS_PATH += -I$(DIR_FDPS)/phantom_grape_x86/G5/newton/libpg5/
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   += -L$(DIR_FDPS)/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif

PROGRAM = nbody.out

nbody.out: nbody.cpp user-defined.hpp Makefile
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $< $(CLIBS)

asm: nbody.cpp user-defined.hpp Makefile
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -S nbody.cpp $(CLIBS)

clean:	
	rm -f $(PROGRAM)

distclean:	
	rm -f $(PROGRAM) *~
