PS_PATH = -I /home/ra000014/a04092/FDPS-master/src/
PS_PATH += -I /home/ra000014/a04092/eigen-3.3.9 #Eigen
PS_PATH += -I ../header

CC = mpiFCCpx
CFLAGS = -O3 -fopenmp -std=c++11
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

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

clean:
	-rm *.out *.o
