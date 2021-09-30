PS_PATH = -I ../../../src/

CXX = g++
CXXFLAGS = -O3 -ffast-math -funroll-loops 
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += #-Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION

#use_phantom_grape_x86 = yes

ifeq ($(use_phantom_grape_x86),yes)
PS_PATH += -I ../../../src/phantom_grape_x86/G5/newton/libpg5/
CXXFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CCCLIBS   += -L ../../src/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = sph.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CXX) $(CXXFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(CXXLIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CXX) -c $< $(CXXFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

clean:
	-rm *.out *.o
	
