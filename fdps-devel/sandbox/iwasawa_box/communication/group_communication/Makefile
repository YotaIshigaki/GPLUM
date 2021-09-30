#CC = mpicxx
CCC = mpiFCCpx
#CCC = CC

#PS_PATH = -I../../../../src/
#CFLAGS = -O3 -Wall -DNDEBUG
#CFLAGS = -O3 -Wall
#CFLAGS = -O3
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS = -O3 -ffast-math -funroll-loops -Wall -DNDEBUG
#CFLAGS = -O3 -ffast-math -funroll-loops -Wall

#all:main_group_comm.out
all:main_gather.out

#all:main3.out
#all:main2.out
#all:main32.out main64.out

#OPT = -Kfast
OPT = -Kfast,openmp
OPT = -Wall
#OPT = -Kfast,openmp -D__NSIMD__ -D__USING_K__ -NRtrap
#OPT = -Kopenmp -D__NSIMD__ -D__USING_K__ -NRtrap

main_group_comm.out:main_group_comm.cc
	$(CCC) $(OPT) -o $@ $<

main_gather.out:main_group_comm.cc
	$(CCC) $(OPT) -DGATHER -o $@ $<

#main3.out:main3.cc comm.hpp
#	$(CCC) $(CFLAGS) $(OPT) -o $@ $<

#main2.out:main2.cc comm.hpp
#	$(CCC) $(CFLAGS) $(OPT) -o $@ $<

#main32.out:main.cc comm.hpp
#	$(CCC) -DBIT32 $(CFLAGS) $(OPT) -o $@ $<

#main64.out:main.cc comm.hpp
#	$(CCC) -DBIT64 $(CFLAGS) $(OPT) -o $@ $<

#main32_omp.out:main.cc comm.hpp
#	$(CCC) -DBIT32 -D__OPENMP__ $(CFLAGS) $(OPT) -o $@ $<

#main64_omp.out:main.cc comm.hpp
#	$(CCC) -DBIT64 -D__OPENMP__ $(CFLAGS) $(OPT) -o $@ $<

clean:
	rm *.s *.o nbody.out 