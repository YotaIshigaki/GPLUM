#CC = mpicxx
CCC = mpiFCCpx
#CCC = CC
#CCC = /opt/FJSVtclang/GM-1.2.0-21/bin/mpiFCCpx

PS_PATH = -I./src/
#PS_PATH = -I../../../../../fdps/src

CFLAGS = -O3
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CFLAGS += -std=c++11
CFLAGS += -Kfast,openmp
#CFLAGS = -O3 -ffast-math -funroll-loops -Wall -DNDEBUG
#CFLAGS = -O3 -ffast-math -funroll-loops -Wall

all:main_group_comm.out
#all:main_group_comm_uniform.out

#OPT = -Kfast
#OPT = -Kfast,openmp
#OPT = -Kfast,openmp -D__NSIMD__ -D__USING_K__ -NRtrap
#OPT = -Kopenmp -D__NSIMD__ -D__USING_K__ -NRtrap

main_group_comm.out:main_group_comm.cc comm.hpp init.hpp
	$(CCC) $(PS_PATH) $(CFLAGS) $(OPT) -o $@ $<

a.out:main_group_comm.s
	$(CCC) -o $@ $<

main_group_comm.s:main_group_comm.cc comm.hpp init.hpp
	$(CCC) -S $(PS_PATH) $(CFLAGS) $(OPT) -o $@ $<

main_group_comm.o:main_group_comm.cc comm.hpp init.hpp
	$(CCC) -c $(PS_PATH) $(CFLAGS) $(OPT) -o $@ $<

test.out:test.cc
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

clean:
	rm *.s *.o *.out 
