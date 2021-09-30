use_mic = yes
#use_mic = no

PS_PATH = -I../../../../src/

CC = mpiicpc
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -openmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DUSE_STD_SORT

HDR = user-defined.hpp vector_knc.hpp
ETC = Makefile

ifeq ($(use_mic),yes)
CFLAGS += -mmic
use_intrinsics = yes
endif
ifeq ($(use_mic),no)
use_phantom_grape_x86 = yes
endif

ifeq ($(use_phantom_grape_x86),yes)
PS_PATH += -I../../src/phantom_grape_x86/G5/newton/libpg5/
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   += -L../../src/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif
ifeq ($(use_intrinsics),yes)
CFLAGS  += -DUSE_INTRINSICS
endif

PROGRAM = nbody.out

nbody.out: nbody.cpp $(HDR) $(ETC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $< $(CLIBS)

asm: nbody.cpp $(HDR) $(ETC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -S nbody.cpp $(CLIBS)

clean:	
	rm -f *~ std*.log

distclean:	
	make clean
	rm -f *.s
	rm -f $(PROGRAM)
