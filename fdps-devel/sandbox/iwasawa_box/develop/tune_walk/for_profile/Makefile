#INC_DIR = -I/home/masaki/project/fdps/src/ 
#INC_DIR = -I../../../../../src/
INC_DIR = -I./src/
#INC_DIR += -I/usr/mpi/gcc/openmpi-1.8.4/include/

#DEFS =  -DENABLE_PEZY 
DEFS = -DPARTICLE_SIMULATOR_MPI_PARALLEL
#DEFS += -DMPICH_IGNORE_CXX_SEEK 
#DEFS += -xHOST -O3 -ipo -no-prec-div
#DEFS += -DMPICH_SKIP_MPICXX
#DEFS += -DMONAR
#DEFS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#DEFS += -fopenmp

#CC     = mpiicpc -g
CC     = time CC

CCOPT  = -O3
CCOPT  += -march=core-avx2
CCOPT  += -DFAST_WALK_AVX2
CCOPT  += -Wall

all: nbody.out

nbody.out:nbody.o
	${CC} $< -o $@

nbody.o:nbody.cpp
	${CC} ${DEFS} ${CCOPT} ${INC_DIR} -c $< -o $@

clean:
	rm *.out *.o
