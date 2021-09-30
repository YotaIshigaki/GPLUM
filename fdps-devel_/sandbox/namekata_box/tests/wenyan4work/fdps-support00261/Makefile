PS_PATH = ../../../../../src/
INC = -I$(PS_PATH)

#CC = time g++-mp-8
CC = time g++ -std=c++11
CFLAGS = -O3
CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

TARGET = a.out

all: $(TARGET)

a.out: main.cpp $(OBJS)
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(CLIBS)

clean:	
	rm -f *.o

distclean:
	make clean
	rm -f $(TARGET) *.txt
