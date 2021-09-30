#PS_PATH = -I ../../../src/
PS_PATH = -I ../../../../src/

#CXX = time mpiFCCpx
#CXXFLAGS = -Kfast
#CXXFLAGS += -Koptmsg=2
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -x32
#CXXFLAGS += -Xg
#CXXFLAGS += -DFAST_COMM_FOR_2ND_EXCHANGE
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K

#CXX = g++
#CXX = mpic++-mpich-gcc49
CXX = CC
#CXXFLAGS = -O3 -ffast-math -funroll-loops 
CXXFLAGS += -O3 -march=core-avx2
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += #-Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = sph.out

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
	rm -f *.o test.o* *~

distclean:
	make clean
	rm -f *.out
