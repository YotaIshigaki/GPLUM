#!/bin/bash -x
#PJM -N "prng"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM -g G19031
#PJM -L rscunit=bwmpc
#PJM -L rscgrp=batch
#PJM --norestart
#---- for test (sequential exec.) ----
#PJM -L vnode=1
#PJM -L vnode-core=1
#PJM -L elapse=0:10:00
#-------------------------------------
##PJM -m b,e
##PJM --mail-list daisuke.namekata@riken.jp
#PJM -s
#PJM --spath "stat-%j.log"
#
# Perform the job
./a.out

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
