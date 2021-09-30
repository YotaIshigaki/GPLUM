#!/bin/bash
NPROC=40000
NPROC_PER_JOB=16
NJOB=$((${NPROC} / ${NPROC_PER_JOB}))

SRC="./submit_job_debug_deep/"
PRG_ORG="./nbody-namekata.out"

for ((i=0; i < ${NJOB}; i++)); do
    rm -rf copy_${i}
    cp -r ${SRC} copy_${i}
    cd copy_$i
    PRG=$(printf "chk-%05d.out" $i)
    mv $(PRG_ORG) $(PRG)
    bsub -I -b -q q_sw_gb -n $(NPROC_PER_JOB) -cgsp 64 -share_size 7000 -host_stack 450 $(PRG) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1 &
    cd ../
done
