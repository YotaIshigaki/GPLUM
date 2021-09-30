PS_PATH = ../../../../src
use_phantom_grape_x86 = yes

CC = mpicxx
CFLAGS = -O3 -ffast-math -funroll-loops -Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CFLAGS += -I$(PS_PATH)
ifeq ($(use_phantom_grape_x86),yes)
CFLAGS += -I$(PS_PATH)/phantom_grape_x86/G5/newton/libpg5/
CFLAGS += -DENABLE_PHANTOM_GRAPE_X86
CLIBS  += -L$(PS_PATH)/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif

PROGRAM = nbody.out

nbody.out: nbody.cpp user-defined.hpp
	$(CC) $(CFLAGS) -o nbody.out nbody.cpp $(CLIBS)

clean:	
	rm -f *~ *.o

distclean:	
	rm -f $(PROGRAM)
