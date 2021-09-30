PS_PATH = -I ../../../src/

CXX = time mpiFCCpx
CXXFLAGS = -Kfast
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -DFAST_COMM_FOR_2ND_EXCHANGE

CXXFLAGS += -DWORKING_ON_K

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = sph.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CXX) $(CXXFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"
	cat *.cpp *.h > print.out

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CXX) -c $< $(CXXFLAGS) $(WARNINGS) $(PS_PATH)
	touch print.out
	@echo "[$< OK]"


clean:
	-rm *.out *.o

