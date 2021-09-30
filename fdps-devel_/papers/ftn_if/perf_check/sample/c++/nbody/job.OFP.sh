#!/bin/bash -x
#PJM -g xg17i056 
#PJM --name "nbody.cpp"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM -L rscgrp=regular-flat
#PJM -L node=1
#PJM -L elapse=08:00:00
#PJM -m b,e
#PJM --mail-list daisuke.namekata@riken.jp
##PJM --stg-transfiles all
##PJM --stgin "./nbody.out ./"
##PJM -s
##PJM --spath "stat-%j.log"
#
module load intel
module load impi
ulimit -s unlimited
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=512M
# Make output directory
#DIR=data-ofp-serial
DIR=data-ofp-empty
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
