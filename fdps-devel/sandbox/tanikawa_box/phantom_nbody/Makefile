PS_PATH  = -I../../../src/
PS_PATH += -I../../../src/phantom_grape_x86/G5/newton/libpg5/

#CC = g++
CC = mpicxx
CFLAGS = -O3 -ffast-math -funroll-loops -Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION

CLIBS = -L../../../src/phantom_grape_x86/G5/newton/libpg5/ -lpg5

SRC = nbody.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $< $(CLIBS)

clean:
	rm -f $(PROGRAM) *.o *~

