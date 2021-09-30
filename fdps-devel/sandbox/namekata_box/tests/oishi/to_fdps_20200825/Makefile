PS_PATH = -I /work/gn99/share/FDPS-master/5.0g/src/
PS_PATH += -I /work/gn99/share/ #Eigen
PS_PATH += -I ./header

#CC = g++
CC = mpiicpc
#CC = mpicxx
CFLAGS = -axMIC-AVX512 -O3 -std=c++11 -ipo -no-prec-div# -parallel# -qopenmp
#CFLAGS = -O3 -ffast-math -funroll-loops -std=c++11# -pg
CFLAGS += -qopt-report -qopt-report-phase=vec
# CFLAGS = -w -Kfast -Kparallel -Xg #-O3
# CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
# CC = mpiFCCpx
#CC = mpiFCC
#CFLAGS =  -Kfast -Kparallel -Xg -O3
#CFLAGS = -O0 -Nquickdbg
# CFLAGS += -Kopenmp
#CFLAGS += -g -O0

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
