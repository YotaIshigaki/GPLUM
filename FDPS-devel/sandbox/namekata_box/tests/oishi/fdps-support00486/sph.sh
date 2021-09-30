#!/bin/sh
#PJM -L  "node=256"
#PJM -L  "rscgrp=eap-small"
#PJM -L  "elapse=02:00:00"
#PJM --mpi "shape=256"
#PJM --mpi "max-proc-per-node=4"
#PJM -s


#------ program execution ------#
 mpiexec ./sph.out


