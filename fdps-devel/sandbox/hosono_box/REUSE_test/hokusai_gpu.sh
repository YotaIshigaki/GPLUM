#!/bin/sh

#PJM -L rscunit=gwacsg
#PJM -L rscgrp=gpu
#PJM -L vnode=1
#PJM -L vnode-core=12
#PJM -L elapse=00:10:00
#PJM -x gpu_per_vnode=1
#PJM -g Q16324
#PJM -j

export OMP_NUM_THREADS=12
export PARALLEL=12

./sph.out

