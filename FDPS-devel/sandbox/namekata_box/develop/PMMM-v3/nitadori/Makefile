CXX = g++
CXXFLAGS = -Wall -O2 
#--- Options for PBC ----
#CXXFLAGS += -DCOMPARE_WITH_EWALD_METHOD
#CXXFLAGS += -DEWALD_FILE
CXXFLAGS += -DNON_CHARGE_NEUTRAL
CXXFLAGS += -DCHECK_COSMO_SIM
CXXFLAGS += -DREAD_INITIAL_DATA

LIBFFTW = -I/work4a/local/fftw/include -L/work4a/local/fftw/lib -lfftw3
HEADERS = vector3.h fmm.h ewald.h cell.h particle.h 

all: obc pbc ewald

pbc: pbc.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $< $(LIBFFTW) -o $@

obc: obc.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $< $(LIBFFTW) -o $@

ewald: ewald.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f obc pbc ewald ?diff*.dat

distclean:
	make clean
	rm -f *.txt
