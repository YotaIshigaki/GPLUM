#!/bin/bash
DIR=data-ajisai-empty
mkdir -p ${DIR}
#NMEAS=32
NMEAS=128
for ((i=0; i < ${NMEAS}; i++)); do
    echo "Start ${i}th measuturement."
    echo 1024 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_10-${i}th.log 2>&1
    echo 2048 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_11-${i}th.log 2>&1
    echo 4096 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_12-${i}th.log 2>&1
    echo 8192 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_13-${i}th.log 2>&1
    echo 16384 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_14-${i}th.log 2>&1
    echo 32768 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_15-${i}th.log 2>&1
    echo 65536 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_16-${i}th.log 2>&1
    echo 131072 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_17-${i}th.log 2>&1
    echo 262144 > nptcl.input
    ./nbody.out > ${DIR}/00stdout-2_18-${i}th.log 2>&1
done
