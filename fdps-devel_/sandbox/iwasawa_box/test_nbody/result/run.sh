#!/bin/bash -x
#PJM --rsc-list "node=4096"
#PJM --rsc-list "elapse=1:00:00"
#PJM --rsc-list "rscgrp=large"
#PJM -S
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ../test_nbody.out %r:./a.out"
#PJM --stgout "rank=* %r:./stdf.%r ./n512M_np4096_j%j/"
#PJM --stgout "rank=* %r:./out/* ./n512M_np4096_j%j/out%r/"

#
. /work/system/Env_base
#
export OMP_NUM_THREADS=8
export PARALLEL=8

mpiexec /work/system/bin/msh "mkdir ./out"

mpiexec -of-proc stdf ./a.out -N 536870912 -n 256 -T 10.0 -t 0.5 -o ./out/
