#!/bin/bash -x
#PJM --name "Making_Tables"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM --rsc-list "node=8"
#PJM --rsc-list "elapse=01:00:00"
#PJM --rsc-list "rscgrp=small"
#PJM --mpi "shape=8, proc=64"
#PJM -m b,e
#PJM --mail-list namekata@ccs.tsukuba.ac.jp
#PJM --stg-transfiles all
#PJM --stgin "./fujitsu.tar.gz ./"
#PJM --stgin "./data.tar.gz ./"
#PJM --stgin "./mktbl.x ./"
#PJM --stgout-dir "./tables ./tables"
#PJM -s
#PJM --spath "stat-%j.log"
#
. /work/system/Env_base
export GSL_LIB=./fujitsu/gsl-1.16/lib
export QD_LIB=./fujitsu/lib
export LD_LIBRARY_PATH=${GSL_LIB}:${QD_LIB}:${LD_LIBRARY_PATH}
#
tar -xvzf fujitsu.tar.gz
tar -xvzf data.tar.gz
mpiexec -n 64 ./mktbl.x -Wl,-T

# [Notes]
#   1) In the K-computer, 1 node have 1 CPU which has 8 cores.
#   2) "-j" option merges STDERR to STDOUT.
#   3) "--norestart" option stops the system restarting this job. 
#   4) If "-m" option is given, the system informs the user abount job status
#      by e-mails. The e-mail address should be specified by "--mail-list"
#      option. 
#   5) "-s" option activate the output of statistical information on job. 
#      We can specify the name of output file by "--spath" option.
#   6) With "-Wl,-T" option, we can perfrom IO with little endian.
