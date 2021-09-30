#!/bin/bash -x
#PJM -N "nbodysph.cpp"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM -g Q19459
#PJM -L rscunit=gwmpc
#PJM -L rscgrp=batch
#PJM --norestart
#---- for performance measurement (Small, 32 proc. 8 threads)----
#PJM -L node=8
#PJM --mpi proc=32
#PJM -L elapse=24:00:00
#---- for performance measurement (Large, 128 proc. 8 threads)----
##PJM -L node=32
##PJM --mpi proc=128
##PJM -L elapse=0:30:00
#-------------------------------------
##PJM -m b,e
##PJM --mail-list daisuke.namekata@riken.jp
#PJM -s
#PJM --spath "stat-%j.log"
#
module load fftw
export OMP_NUM_THREADS=8
# Perform the job
mpiexec -n 32 ./nbodysph.out
#mpiexec -n 64 ./nbodysph.out
#mpiexec -n 128 ./nbodysph.out
#mpiexec -n 256 ./nbodysph.out
#mpiexec -n 512 ./nbodysph.out

# [Notes]
#   1) In the HOKUSAI Great Wall system, 1 node have 2 CPUs each of which has 16 cores.
#   2) "-j" option merges STDERR to STDOUT.
#   3) "-g" option is used to specify the project ID.
#   4) The system where the program is run is specified by the `-L rscunit=` option.
#      The value `gwmpc` represent the GreatWall system.
#   5) The group of resource is specified by the `-L rscgrp=` option.
#      We have to select a resource group from the following:
#      -----------------------------------------------------------------------------
#         Resource group  | Max. elapsed time | Max. # of nodes  | Max. # of cores
#      -----------------------------------------------------------------------------
#         batch           |   72h             |    16            |    512
#                         |   24h             |   256 (general)  |   8192
#                         |   24h             |    32 (quick)    |   1024
#         gaussian        |   72h             |     1            |     32
#         special         |   48h             |  1080            |  34560
#      -----------------------------------------------------------------------------
#
#   6) "--norestart" option stops the system restarting this job. 
#   7) If "-m" option is given, the system informs the user abount job status
#      by e-mails. The e-mail address should be specified by "--mail-list"
#      option. 
#   8) "-s" option activates the output of statistical information on job. 
#      We can specify the name of output file by "--spath" option.
#   9) With "-Wl,-T" option, we can perfrom IO with little endian.
#
#
