CXX = g++-mp-4.8 -Wall -O2
LIBFFTW = -I /opt/local/include -L /opt/local/lib -lfftw3
HEADERS = vector3.h fmm.h ewald.h cell.h particle.h 

all: obc pbc ewald

pbc: pbc.cpp $(HEADERS)
	$(CXX) -fopenmp $< $(LIBFFTW) -o $@

obc: obc.cpp $(HEADERS)
	$(CXX) -fopenmp $< $(LIBFFTW) -o $@

ewald: ewald.cpp $(HEADERS)
	$(CXX) -fopenmp $< -o $@

clean:
	rm -f obc pbc ewald ?diff*.dat
