INC = -I ../../../src
CC = g++
#CC = mpic++
CFLAGS = -O3 -ffast-math -funroll-loops #-std=c++11
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
CFLAGS += -DFAST_COMM_FOR_2ND_EXCHANGE
#CFLAGS += -mavx

#use_phantom_grape_x86 = yes

ifeq ($(use_phantom_grape_x86),yes)
INC += -I ../../../src/phantom_grape_x86/G5/newton/libpg5/
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
LIBS   = -L ../../../src/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif

CPPOBJS = main.cpp
CPPHDRS = $(wildcard *.h)
PROGRAM = SPzyH.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(INC)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Building $< ..."
	@$(CC) -c $< $(CFLAGS) $(WARNINGS) $(INC)
	@echo "[$< OK]"

clean:
	-rm *.out *.o *.log


