PS_PATH = -I../../src/ 
FFTW_LIB = -I/work3/tanikawa/mylib/include/ -L/work3/tanikawa/mylib/lib -lfftw3f_mpi -lfftw3f
PM_LIB =  -I../../src/particle_mesh/ -L../../src/particle_mesh -lpm

CC = mpicxx
CFLAGS = -O3 -ffast-math -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

MULEXP = -DMONOPOLE

SRC = test_PM.cc
PROGRAM = $(SRC:%.cc=%.out)

$(PROGRAM):$(SRC) 
#	$(CC) $(MULEXP) $(PS_PATH) $(PM_LIB) $(FFTW_LIB) $(CFLAGS) -o $@ $<
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $< $(PM_LIB) $(FFTW_LIB)

clean:
	rm $(PROGRAM)
