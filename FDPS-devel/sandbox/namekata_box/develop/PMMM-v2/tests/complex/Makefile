FFTW_LOC  = /opt/local
FFTW_INC  = -I$(FFTW_LOC)/include
FFTW_LIB  = -L$(FFTW_LOC)/lib -lfftw3 

FDPS_LOC  = ../../
FDPS_INC  = -I$(FDPS_LOC)/src
FDPS_LIB  = 
LDFLAGS   = $(FFTW_LIB)

CXX = g++-mp-8
CXXFLAGS = -Wall -O0
#CXXFLAGS = -O3 -ffast-math -funroll-loops
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += $(FDPS_INC) $(FFTW_INC)

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = a.out

.PHONY:	clean all

all: $(CPPOBJS) $(CPPHDRS)
	$(CXX) $(CXXFLAGS) $(CPPOBJS) -o $(PROGRAM) $(LDFLAGS)

%.o: %.cpp $(CPPHDRS)
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -f $(CPPOBJS)

distclean: clean
	rm -f $(PROGRAM)
