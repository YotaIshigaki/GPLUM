#!/bin/bash -x
#PJM --name "nbody.ftn"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM --rsc-list "node=8x1"
#PJM --rsc-list "elapse=01:00:00"
#PJM --rsc-list "rscgrp=small"
#PJM -m b,e
#PJM --mail-list daisuke.namekata@riken.jp
#PJM --stg-transfiles all
#PJM --stgin "./nbody.out ./"
#PJM -s
#PJM --spath "stat-%j.log"
#
. /work/system/Env_base
ulimit -s unlimited
export OMP_NUM_THREADS=1
export PARALLEL=1
export OMP_STACKSIZE=512M
mkdir -p result
#
./nbody.out -Wl,-T
#mpiexec ./nbody.out -Wl,-T
#mpiexec -n 8 ./nbody.out -debuglib -Wl,-T > 00stdout.log 2>&1
#fipp -C -Icall,hwm -d profile_data -m 40000 mpiexec -n 8 ./nbody.out -Wl,-T
#fapp -C -Impi,hwm -Hevent=Cache,mode=sys,method=normal -d profile_data mpiexec -n 8 ./nbody.out -Wl,-T

# [Notes]
#   1) In the K-computer, 1 node have 1 CPU which has 8 cores.
#   2) "-j" option merges STDERR to STDOUT.
#   3) "--norestart" option stops the system restarting this job. 
#   4) The group of resource is specified by the `--rsc-list "rscgrp="` option.
#      We have to select a resource group from the following:
#      -----------------------------------------------------------
#         Resource group  | Max. # of nodes  | Max. # of cores
#      -----------------------------------------------------------
#         small/interact  |   384            |  3072
#         micro           |   1152           |  9216
#         large           |   36864          |  294912
#         huge            |   82944          |  663552
#      -----------------------------------------------------------
#      Note that the resource group `micro` must be performed in a dedicated
#      directory (/scratch/groupname).
#
#   5) If "-m" option is given, the system informs the user abount job status
#      by e-mails. The e-mail address should be specified by "--mail-list"
#      option. 
#   6) "-s" option activates the output of statistical information on job. 
#      We can specify the name of output file by "--spath" option.
#   7) With "-Wl,-T" option, we can perfrom IO with little endian.
#   8) The profiling data created by fipp can be converted into 
#      a text form by the command:
#      $ fipppx -Icpu,balance,call,hwm -A -d profile_data
#      On the other hand, that created by fapp can be converted by the command:
#      $ fapppx -Impi,hwm -ttext -o profile.txt -A -d profile_data
#
#
