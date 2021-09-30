PS_PATH = -I../../../src/

#CC = g++
CC = mpicxx
CFLAGS = -O3
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DCALC_ENG

INC = -I../inc

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

SRC = main.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(PS_PATH) $(CFLAGS) $(INC) -o $@ $<

clean:
	rm -f $(PROGRAM)

distclean: clean
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
