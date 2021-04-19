#!/bin/bash

#------ pjsub option -----#
#PJM -S
#PJM -L rscgrp=regular-flat
#PJM -L node=8
#PJM --mpi proc=512
##PJM --omp thread=16
#PJM -L elapse=00:10:00
#PJM -g gn99


#------ pinning setting -----#
source /usr/local/bin/mpi_core_setting.sh
# export KMP_AFFINITY=verbose
export I_MPI_DEBUG=5

#------ program execution ------#
 mpiexec.hydra -n ${PJM_MPI_PROC} numactl --preferred=1  ./sph.out


