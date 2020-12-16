#!/bin/bash
DIR=data-ajisai-serial
mkdir -p ${DIR}
for ((i=0; i < ${NMEAS}; i++)); do
    echo "Start ${i}th measuturement."
    ./nbody.out -N 1024   > ${DIR}/00stdout-2_10-${i}th.log 2>&1
    ./nbody.out -N 2048   > ${DIR}/00stdout-2_11-${i}th.log 2>&1
    ./nbody.out -N 4096   > ${DIR}/00stdout-2_12-${i}th.log 2>&1 
    ./nbody.out -N 8192   > ${DIR}/00stdout-2_13-${i}th.log 2>&1
    ./nbody.out -N 16384  > ${DIR}/00stdout-2_14-${i}th.log 2>&1
    ./nbody.out -N 32768  > ${DIR}/00stdout-2_15-${i}th.log 2>&1
    ./nbody.out -N 65536  > ${DIR}/00stdout-2_16-${i}th.log 2>&1
    ./nbody.out -N 131072 > ${DIR}/00stdout-2_17-${i}th.log 2>&1
    ./nbody.out -N 262144 > ${DIR}/00stdout-2_18-${i}th.log 2>&1 
done
