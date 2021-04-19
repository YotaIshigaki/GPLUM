#!/bin/bash -x
#--------------------------------------------------------------------------------
# [Notes]
#    (1) The site below explains how to describe a job script file in ITO.
#        https://www.cc.kyushu-u.ac.jp/scp/system/ITO/04-1_intel_compiler.html
#    (2) You can find the information on the available resource groups
#        at the follwoing site.
#        https://www.cc.kyushu-u.ac.jp/scp/system/ITO/subsystem/06_limit.html
#--------------------------------------------------------------------------------
#PJM -N "test"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM -L rscunit=ito-a
#PJM -L rscgrp=ito-m-dbg
#PJM --norestart
##PJM -L vnode=1
#PJM -L vnode=4
##PJM -L vnode=16
#PJM -L vnode-core=36
#PJM -L elapse=00:15:00
##PJM -m b,e
##PJM --mail-list daisuke.namekata@riken.jp
#PJM -s
#PJM --spath "stat-%j.log"
#
module load intel/2019.4
module load fftw/3.3.8_intel2019.4
# Set the configuration of parallerization
NUM_NODES=$PJM_VNODES
NUM_CORES=36 # DON'T CHANGE!!
NUM_THREADS=4 # MUST BE BETWEEN 1 TO 36
NUM_PROCS=`expr $NUM_NODES "*" $NUM_CORES / $NUM_THREADS` 
# Set the environment variables for MPI & OpenMP
export I_MPI_PERHOST=`expr $NUM_CORES / $NUM_THREADS`
export I_MPI_FABRICS=shm:ofi
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_CELL=core

export OMP_NUM_THREADS=$NUM_THREADS
export OMP_NESTED=1
export OMP_MAX_ACTIVE_LEVELS=2
export KMP_STACKSIZE=16m
export KMP_AFFINITY=compact

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}
# (Option) Set the environment variables for debugging MPI program
#export I_MPI_DEBUG=1000,pid,host,scope,line,file
#export I_MPI_DEBUG_OUTPUT=mpi_debug_output.txt
#export I_MPI_HYDRA_DEBUG=on
# Perform the job
mpiexec.hydra -n $NUM_PROCS ./test.out

