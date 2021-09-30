#PBS -N simd_test
#PBS -j oe
#PBS -q large-t
#PBS -l nodes=1
#PBS -l walltime=00:10:00
#PBS -meb
#PBS -M daisuke.namekata@riken.jp
cd $PBS_O_WORKDIR
source /opt/modules/default/init/bash
module swap PrgEnv-cray PrgEnv-gnu
time aprun -n 2 ./main.out >& 00stdout-${PBS_JOBID}.log
