PS_PATH = -I../../../../src

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

CCC = mpicxx

CFLAGS = -O3 -ffast-math -funroll-loops -Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#TOWDIM = PARTICLE_SIMULATOR_TOW_DIMENSION

OBJ = main.o
HDR = myheader.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ)

main.o: main.cc $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cc

clean:
	rm -f *~
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
