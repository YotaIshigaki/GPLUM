#PS_PATH = -I ../../../src/
PS_PATH = -I ../../../../src/

CXX = time mpiFCCpx
CXXFLAGS = -Kfast
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -Xg
CXXFLAGS += -x32
CXXFLAGS += -DFAST_COMM_FOR_2ND_EXCHANGE
CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
CXXFLAGS += -DINTRINSIC_K
CXXFLAGS += -DFAST_WALK_K
#CXXFLAGS += -DPROFILE_K
#CXXFLAGS += -Koptmsg=2
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION

#CXX = time g++
#CXX = time CC
#CXX = time mpic++-mpich-gcc49
#CXX = time mpicxx
#CXXFLAGS = -O3 -Wall -ffast-math -funroll-loops  -march=core-avx2
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += #-Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
#PROGRAM = sph_nx203.out
#PROGRAM = sph_nx406.out
#PROGRAM = sph_nx812.out
#PROGRAM = sph_nx102.out
PROGRAM = sph_nx1394.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CXX) $(CXXFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CXX) -c $< $(CXXFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

#sph.s:	main.cpp $(CPPHDRS)
#	$(CXX) -S -c $(CXXFLAGS) main.cpp -o sph.s $(LIBS) $(PS_PATH)

#clean:
#	rm sph.s

clean:
	-rm *.o
cleanall:
	-rm *.out *.o