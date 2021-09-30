#!/bin/bash -x
#PJM -g xg17i056 
#PJM --name "nbody.ftn"
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
mkdir -p result
# Make output directory
DIR=data-ofp-empty
#DIR=data-ofp-serial
mkdir -p ${DIR}
# Perform jobs
NMEAS=256
for ((i=0; i < ${NMEAS}; i++)); do
   echo "Start ${i}th measurement."
   id=`printf "%03d" ${i}`
   echo 1024 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_10-${id}th.log 2>&1
   echo 2048 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_11-${id}th.log 2>&1
   echo 4096 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_12-${id}th.log 2>&1
   echo 8192 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_13-${id}th.log 2>&1
   echo 16384 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_14-${id}th.log 2>&1
   echo 32768 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_15-${id}th.log 2>&1
   echo 65536 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_16-${id}th.log 2>&1
   echo 131072 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_17-${id}th.log 2>&1
   echo 262144 > nptcl.input
   ./nbody.out > ${DIR}/00stdout-2_18-${id}th.log 2>&1
done
