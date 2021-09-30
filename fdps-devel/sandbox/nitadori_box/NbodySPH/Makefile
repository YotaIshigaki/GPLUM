PS_PATH = -I ../../../src/

CC = mpiFCCpx
CFLAGS = -Kfast -Xg -x32 -Koptmsg=2
#CFLAGS += -Knoswp,nounroll
#CC = g++
#CC = mpic++-mpich-gcc49
#CFLAGS = -Wall -O2 -ffast-math
#CFLAGS += -march=core-avx2
#CFLAGS = -O3 -ffast-math -funroll-loops 
#CFLAGS += -Wall -march=core-avx2
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += #-Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = sph.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CC) -c $< $(CFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

%.s:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CC) -S $< $(CFLAGS) $(WARNINGS) $(PS_PATH) &>> optmsg.log
	@echo "[$< OK]"

clean:
	-rm *.out *.o
