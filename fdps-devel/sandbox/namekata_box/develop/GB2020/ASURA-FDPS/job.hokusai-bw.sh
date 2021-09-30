#!/bin/bash -x
#PJM -N "check_n_jp_dist"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM -g G19031
#PJM -L rscunit=bwmpc
#PJM -L rscgrp=batch
#PJM --norestart
#---- for test (sequential exec.) ----
##PJM -L vnode=1
##PJM -L vnode-core=1
##PJM -L elapse=1:00:00
#---- for test (parallel exec.,openMP) ----
##PJM -L vnode=1
##PJM -L vnode-core=20
##PJM -L elapse=0:30:00
#---- for test (parallel exec.,flat mpi) ----
#PJM -L vnode=2
#PJM -L vnode-core=40
#PJM -L elapse=00:15:00
##---- for test (parallel exec.,flat mpi) ----
##PJM -L vnode=16
##PJM -L vnode-core=40
##PJM -L elapse=24:00:00
##PJM -L elapse=0:15:00
#---- for performance measurement (Small, 32 proc. 10 threads)----
##PJM -L vnode=8
##PJM -L vnode-core=40
##PJM -L elapse=0:30:00
#---- for performance measurement (Mid., 64 proc. 10 threads)----
##PJM -L vnode=16
##PJM -L vnode-core=40
##PJM -L elapse=00:10:00
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
# Perform the job
# (0) Serial exec.
#./nbodysph.out
# (1) Flat MPI
mpirun -np 40 -ppn 40 ./nbodysph.out
#mpirun -np 80 -ppn 40 ./nbodysph.out
#mpirun -np 160 -ppn 40 ./nbodysph.out
#mpirun -np 320 -ppn 40 ./nbodysph.out
#mpirun -np 640 -ppn 40 ./nbodysph.out
#mpirun -np 640 -ppn 40 ./nbodysph.out -r 14 # check_n_jp_dist
##mpirun -np 640 -ppn 40 ./nbodysph.out -r 50
##mpirun -np 640 -ppn 40 ./nbodysph.out -r 113
# (2) hybrid
#export OMP_NUM_THREADS=10
#mpirun -np 32 -ppn 4 ./nbodysph.out
#mpirun -np 64 -ppn 4 ./nbodysph.out
#mpirun -np 128 -ppn 4 ./nbodysph.out
#mpirun -np 256 -ppn 4 ./nbodysph.out
#mpirun -np 1024 -ppn 4 ./nbodysph.out

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
