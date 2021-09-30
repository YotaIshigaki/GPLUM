#!/bin/bash -x
#PJM -N "cosmology"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM -g Q19459
#PJM -L rscunit=bwmpc
#PJM -L rscgrp=batch
#PJM --norestart
#---- for test (sequential exec.) ----
##PJM -L vnode=1
##PJM -L vnode-core=1
##PJM -L elapse=3:00:00
#---- for test (parallel exec.,openMP) ----
##PJM -L vnode=1
##PJM -L vnode-core=20
##PJM -L elapse=0:30:00
#---- for test (parallel exec.,flat mpi) ----
##PJM -L vnode=1
##PJM -L vnode-core=40
##PJM -L elapse=0:30:00
#---- for performance measurement (Small, 32 proc. 5 threads)----
#PJM -L vnode=4
#PJM -L vnode-core=40
#PJM -L elapse=0:30:00
#---- for performance measurement (Small, 32 proc. 10 threads)----
##PJM -L vnode=8
##PJM -L vnode-core=40
##PJM -L elapse=0:30:00
#---- for performance measurement (Large, 128 proc. 10 threads)----
##PJM -L vnode=32
##PJM -L vnode-core=40
##PJM -L elapse=0:30:00
#-------------------------------------
##PJM -m b,e
##PJM --mail-list daisuke.namekata@riken.jp
#PJM -s
#PJM --spath "stat-%j.log"
#
module load fftw
export OMP_NUM_THREADS=5
#export OMP_NUM_THREADS=10
#export OMP_NUM_THREADS=20
# Perform the job
#export PARAM_FILE=../../param_files/random.para
export PARAM_FILE=../../param_files/santa_barbara.para
#./cosmology.out ${PARAM_FILE}
mpirun -np 32 -ppn 8 ./cosmology.out ${PARAM_FILE}
#mpirun -np 32 -ppn 4 ./cosmology.out ${PARAM_FILE}
#mpirun -np 64 -ppn 4 ./cosmology.out ${PARAM_FILE}
#mpirun -np 128 -ppn 4 ./cosmology.out ${PARAM_FILE}
#mpirun -np 256 -ppn 4 ./cosmology.out ${PARAM_FILE}
#mpirun -np 1024 -ppn 4 ./cosmology.out ${PARAM_FILE}

# [Notes]
#   1) In the HOKUSAI Big Waterfall system, 1 node have 2 CPUs each of which has 20 cores.
#   2) "-j" option merges STDERR to STDOUT.
#   3) "-g" option is used to specify the project ID.
#   4) The system where the program is run is specified by the `-L rscunit=` option.
#      The value `gwmpc` represent the GreatWall system.
#   5) The group of resource is specified by the `-L rscgrp=` option.
#      We have to select a resource group from the following:
#      -----------------------------------------------------------------------------
#         Resource group  | Max. elapsed time | Max. # of nodes  | Max. # of cores
#      -----------------------------------------------------------------------------
#         batch           |   72h             |    16            |    640
#                         |   24h             |   128 (general)  |   5120
#                         |   24h             |    32 (quick)    |   1280
#         gaussian        |   72h             |     1            |     40
#         special         |   48h             |   840            |  33600
#      -----------------------------------------------------------------------------
#
#   6) "--norestart" option stops the system restarting this job. 
#   7) If "-m" option is given, the system informs the user abount job status
#      by e-mails. The e-mail address should be specified by "--mail-list"
#      option. 
#   8) "-s" option activates the output of statistical information on job. 
#      We can specify the name of output file by "--spath" option.
#   9) With "-Wl,-T" option, we can perfrom IO with little endian.
#   10) "-L vnode=" option specifies the number of nodes.
#   11) "-L vnode-core=" option specifies the number of cores per node used in a job. 
#   12) The configuration of parallelization is specified by the environment variable
#       OMP_NUM_THREADS and mpirun's option "-np X -ppn Y", where X is the total number
#       of MPI processes and Y is the number of MPI processes per nodes.
#
