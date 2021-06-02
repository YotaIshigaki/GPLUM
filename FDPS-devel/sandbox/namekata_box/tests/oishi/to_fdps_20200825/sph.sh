#!/bin/bash

#------ pjsub option -----#
#PJM -S
#PJM -L rscgrp=regular-flat
#PJM -L node=16
#PJM --mpi proc=1024
##PJM --omp thread=16
#PJM -L elapse=5:00:00
#PJM -g po8020


#------ pinning setting -----#
# source /usr/local/bin/hybrid_core_setting.sh 2
source /usr/local/bin/mpi_core_setting.sh
# export KMP_AFFINITY=verbose
export I_MPI_DEBUG=5


#------ program execution ------#
mpiexec.hydra -n ${PJM_MPI_PROC} numactl --preferred=1 ./sph.out
# mpiexec.hydra -n ${PJM_MPI_PROC} ./sph.out
