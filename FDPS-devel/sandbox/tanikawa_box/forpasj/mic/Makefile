################################################
use_xc = no
use_coma = yes
################################################

ifeq ($(use_xc),yes)
CCC = CC
CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops -lm
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif

ifeq ($(use_coma),yes)
CCC = mpiicpc
CFLAGS += -mmic
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -openmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif

FDPS_DIR = ../../../../src
PS_PATH  = -I${FDPS_DIR}

OBJ = main.o
HDR = particle.hpp force.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) $(FFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp

clean:
	rm -f *~
	rm -f stderr.log stdout.log
	rm -f core
	rm -f *.s
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
