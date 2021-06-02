PS_PATH = ../../../../../src/
INC = -I$(PS_PATH)

CC = time g++-mp-6
CFLAGS = -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

all:a.out

a.out: main.cpp $(OBJS)
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(CLIBS)

run:
	./a.out > 00stdout.log 2>&1

clean:	
	rm -f a.out

distclean:	
	rm -f a.out *.o
