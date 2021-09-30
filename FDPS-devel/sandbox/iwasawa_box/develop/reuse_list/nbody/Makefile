#PS_PATH = ../../../../../src/
PS_PATH = ./src/
INC = -I$(PS_PATH)

#CC = time g++
#CC = time CC
CC = mpicxx
CFLAGS = -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CFLAGS += -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

OBJS = kernel_sunway.o

all:nbody.out

nbody.out:nbody.cpp $(OBJS)
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(CLIBS)

kernel_sunway.o: kernel_sunway.c kernel_sunway.h
	gcc -O3 -c -o $@ kernel_sunway.c -lm

clean:	
	rm -f nbody.out

distclean:	
	rm -f nbody.out *.o *~
