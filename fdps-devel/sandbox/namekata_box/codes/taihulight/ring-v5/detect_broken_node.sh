#!/bin/bash
max=10000

SRC="./submit_job_debug_deep/"
PRG="./nbody-namekata.out"

for((i=0; i < $max; i++)); do
rm -rf copy_$i
cp -r $SRC copy_$i
cd copy_$i

bsub -I -b -q q_sw_gb -n 16 -cgsp 64 -share_size 7000 -host_stack 450 $PRG -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1 &
cd ../
done
