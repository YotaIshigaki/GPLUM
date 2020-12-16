#CC = mpicxx
#CCC = mpiFCCpx
CCC = CC
#CFLAGS = -O3 -Wall -DNDEBUG
CFLAGS = -O3
#CFLAGS = -O3 -ffast-math -funroll-loops -Wall -DNDEBUG
#CFLAGS = -O3 -ffast-math -funroll-loops -Wall

all:main.out
#all:main2.out

#OPT = -Kfast,openmp
#OPT = -Kfast,openmp -D__NSIMD__ -D__USING_K__ -NRtrap
#OPT = -Kopenmp -D__NSIMD__ -D__USING_K__ -NRtrap


main.out:main.cc 
	$(CCC) $(CFLAGS) $(OPT) -o $@ $<
main2.out:main2.cc 
	$(CCC) $(CFLAGS) $(OPT) -o $@ $<

clean:
	rm *.s *.o nbody.out 