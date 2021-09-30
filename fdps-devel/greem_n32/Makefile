PS_PATH  = -I../src
PS_PATH += -I../src/particle_mesh

CCC = mpicxx
#CCC = /nfshome/tanikawa/work3/lib/mpich/bin/mpicxx

CFLAGS = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#TOWDIM = PARTICLE_SIMULATOR_TOW_DIMENSION

INCLUDE_FFTW = -I/nfshome/tanikawa/work3/mylib/include
LIB_FFTW = -L/nfshome/tanikawa/work3/mylib/lib -lfftw3f_mpi -lfftw3f
#INCLUDE_FFTW = -I/nfshome/tanikawa/work3/lib/mpich/include
#LIB_FFTW = -L/nfshome/tanikawa/work3/lib/mpich/lib -lfftw3f_mpi -lfftw3f
LIB_PM   = -L../src/particle_mesh -lpm

CFLAGS += $(INCLUDE_FFTW)
CLIBS += $(LIB_PM) $(LIB_FFTW)

OBJ = main.o
HDR = myheader.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

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
