#!/bin/bash

#PJM -L "rscunit=ito-a" 
#PJM -L "rscgrp=ito-s-dbg"
#PJM -L "vnode=4"
#PJM -L "vnode-core=36"
#PJM -L "elapse=30:00"
#PJM -j
#PJM -X

module load intel/2019.4

NUM_NODES=$PJM_VNODES
NUM_CORES=36
NUM_PROCS=16
NUM_THREADS=9

export I_MPI_PERHOST=`expr $NUM_CORES / $NUM_THREADS`
export I_MPI_FABRICS=shm:ofi
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_CELL=core

export OMP_NUM_THREADS=$NUM_THREADS
export KMP_STACKSIZE=8m
export KMP_AFFINITY=compact

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

mpiexec.hydra -n $NUM_PROCS ./ring_mono.out -N 5000000 -t 0.5 -T 100.0 -e 0 -n 128 -l 8 -L 1 > N5e6_n128_l8_t0.5_e0_L1.dat

