#!/bin/bash -x
#PJM --name "nbody.cpp"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM --rsc-list "node=8x1"
#PJM --rsc-list "elapse=24:00:00"
#PJM --rsc-list "rscgrp=small"
#PJM -m b,e
#PJM --mail-list daisuke.namekata@riken.jp
#PJM --stg-transfiles all
#PJM --stgin "./nbody.out ./"
#PJM --stgout "data.tar.gz ./"
#PJM -s
#PJM --spath "stat-%j.log"
#
. /work/system/Env_base
export OMP_NUM_THREADS=1
export PARALLEL=1
export OMP_STACKSIZE=512M
# Make output directory
DIR=data-k-empty
#DIR=data-k-serial
mkdir -p ${DIR}
# Perform jobs
NMEAS=256
for ((i=0; i < ${NMEAS}; i++)); do
   echo "Start ${i}th measurement."
   id=`printf "%03d" ${i}`
   ./nbody.out -N 1024   > ${DIR}/00stdout-2_10-${id}th.log 2>&1
   ./nbody.out -N 2048   > ${DIR}/00stdout-2_11-${id}th.log 2>&1
   ./nbody.out -N 4096   > ${DIR}/00stdout-2_12-${id}th.log 2>&1
   ./nbody.out -N 8192   > ${DIR}/00stdout-2_13-${id}th.log 2>&1
   ./nbody.out -N 16384  > ${DIR}/00stdout-2_14-${id}th.log 2>&1
   ./nbody.out -N 32768  > ${DIR}/00stdout-2_15-${id}th.log 2>&1
   ./nbody.out -N 65536  > ${DIR}/00stdout-2_16-${id}th.log 2>&1
   ./nbody.out -N 131072 > ${DIR}/00stdout-2_17-${id}th.log 2>&1
   ./nbody.out -N 262144 > ${DIR}/00stdout-2_18-${id}th.log 2>&1 
done
# Compress
tar -cvzf data.tar.gz ${DIR}


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
