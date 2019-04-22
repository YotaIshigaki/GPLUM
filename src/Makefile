PS_PATH = ../FDPS-5.0d/src/ 
INC = -I$(PS_PATH)

#CC = clang++ -std=c++11
#CC = time CC
CC = mpic++ -std=c++11
#CC = /usr/local/mpich2-1.4.1p1-gnu64-4.4.7/bin/mpic++ 
CFLAGS = -O3
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DPARTICLE_SIMULATOR_SPMOM_F32

#CFLAGS += -mtune=skylake-avx512
#CFLAGS += -march=skylake-avx512
CFLAGS += -march=core-avx2 
#CFLAGS += -mavx512f -mavx512dq  
#CFLAGS += -mavx2

CFLAGS += -DUSE_QUAD
CFLAGS += -DUSE_INDIVIDUAL_RADII  
#CFLAGS += -DCALC_EP_64bit
CFLAGS += -DCOLLISION
CFLAGS += -DKOMINAMI
#CFLAGS += -DCHAMBERS
#CFLAGS += -DGAS_DRAG
#CFLAGS += -DCALC_WTIME
CFLAGS += -DOUTPUT_DETAIL
#CFLAGS += -DISOTROPIC
CFLAGS += -DMONAR

SRC = main_p3t.cpp
PROGRAM = gplum.out
HEADER = mathfunc.h kepler.h energy.h particle.h disk.h gravity.h gravity_kernel.hpp collisionA.h collisionB.h hermite.h hard.h read.h time.h func.h

$(PROGRAM): main_p3t.o gravity_kernel.o
	$(CC) $(INC) $(CFLAGS) -o $@ main_p3t.o gravity_kernel.o

main_p3t.o: main_p3t.cpp $(HEADER)
	$(CC) $(INC) $(CFLAGS) -c $<

gravity_kernel.o: gravity_kernel.cpp avx.h particle.h
	$(CC) $(INC) $(CFLAGS) -c $<

clean:
	rm -f *.out *.o *~

