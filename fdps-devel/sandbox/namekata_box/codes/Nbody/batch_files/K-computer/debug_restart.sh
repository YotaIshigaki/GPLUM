#!/bin/bash
#QSUB -J DEBUG.rst
#QSUB -o 00stdout-%J.out
#QSUB -e 00stdout-%J.err
#QSUB -q gr10339d
#QSUB -ug gr10339
#QSUB -A p=64:t=1:c=1:m=2G
#QSUB -W 0:30
#QSUB -u namekata@ccs.tsukuba.ac.jp
#QSUB -B
#QSUB -N
cd $LS_SUBCWD
export LD_LIBRARY_PATH=/LARGE2/gr10339/b33387/intel/gsl-1.16/lib:/LARGE2/gr10339/b33387/intel/lib:$LD_LIBRARY_PATH
./arrange_calc_info.pl $LSB_JOBID
./make_restart_data.pl $LSB_JOBID
time aprun -n $LSB_PROCS -d $LSB_CPUS -N $LSB_PPN ./rhd2d.x.rst
