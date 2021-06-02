FDPS_HOME = ./src/
#FDPS_HOME = ./for_profile/src/
INC = -I$(FDPS_HOME)
INCORG = -I./for_profile/src/

#DEFS =  -DENABLE_PEZY 
#DEFS = -DPARTICLE_SIMULATOR_MPI_PARALLEL
DEFS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#DEFS +=  -qopenmp
DEFS += -fopenmp
#DEFS += -DMULTI_WALK
#DEFS += -DMPICH_IGNORE_CXX_SEEK 
#DEFS += -xHOST -O3 -ipo -no-prec-div
#DEFS += -DMPICH_SKIP_MPICXX
#DEFS += -DMONAR


#CC     = icpc -fast
#CC     = mpiicpc -g
#CC     = time CC
CC     = time mpicxx

CCOPT  = -O3
#CCOPT  += -march=core-avx2
#CCOPT  += -DFAST_WALK_AVX2
CCOPT  += -Wall

#all:nbody.simd.out
all:nbody.out

nbody.out:nbody.cpp
	${CC} -march=core-avx2 ${DEFS} ${CCOPT} ${INC} $< -o $@

nbody.simd.out:nbody.cpp
	${CC} -g -march=core-avx2 -DENABLE_AVX2 -DFAST_WALK_AVX2 ${DEFS} ${CCOPT} ${INC} $< -o $@

nbody.simd.s:nbody.cpp
	${CC} -S -march=core-avx2 -DENABLE_AVX2 -DFAST_WALK_AVX2 ${DEFS} ${CCOPT} ${INC} $< -o $@

nbody.org.simd.out:nbody.cpp
	${CC} ${DEFS} -march=core-avx2 -DFAST_WALK_AVX2 ${CCOPT} ${INCORG} $< -o $@

nbody.org.out:nbody.cpp
	${CC} ${DEFS} -march=core-avx2 ${CCOPT} ${INCORG} $< -o $@

#nbody.out:nbody.o
#	${CC} $< -o $@

#nbody.o:nbody.cpp
#	${CC} ${DEFS} ${CCOPT} ${INC} -c $< -o $@

clean:
	rm nbody.out *.o
