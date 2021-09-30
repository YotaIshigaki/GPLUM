PS_PATH = -I../src_parallel

#CCC = time mpiFCCpx
#CCC = time FCCpx
#CFLAGS = -Kfast,openmp -NRtrap
#CFLAGS = -Kfast,openmp
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -Xg
#CFLAGS += -Nfjcex
#CFLAGS += -Krestp=all
#CFLAGS = -Kfast,openmp,restp=all

CCC = mpicxx-openmpi-gcc49
#CCC = mpicxx-mpich-gcc49
#CCC = g++

CFLAGS = -O3 -ffast-math -funroll-loops #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#TOWDIM = PARTICLE_SIMULATOR_TOW_DIMENSION

DEFINES = -DLONG_CUTOFF

OBJ = main.o
HDR = myheader.hpp myheader_cutoff.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(DEFINES)

main.o: main.cc $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cc $(DEFINES)

clean:
	rm -f *~
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
