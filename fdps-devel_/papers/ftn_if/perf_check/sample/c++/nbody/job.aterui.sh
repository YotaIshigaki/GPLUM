#PBS -N nbody.cpp
#PBS -j oe
#PBS -q large-b
#PBS -l mppwidth=1
#PBS -l walltime=08:00:00
#PBS -meb
#PBS -M cc67803@gmail.com
cd $PBS_O_WORKDIR
source /opt/modules/default/init/bash
module swap PrgEnv-cray PrgEnv-intel
export MPI_COLL_OPT_ON=1
export MPICH_FAST_MEMCPY=1
# Make output directory
DIR=data-aterui-empty
#DIR=data-aterui-serial
mkdir -p ${DIR}
# Perform jobs
NMEAS=256
for ((i=0; i < ${NMEAS}; i++)); do
   echo "Start ${i}th measurement."
   id=`printf "%03d" ${i}`
   time aprun -n 1 ./nbody.out -N 1024   >& ${DIR}/00stdout-2_10-${id}th.log
   time aprun -n 1 ./nbody.out -N 2048   >& ${DIR}/00stdout-2_11-${id}th.log
   time aprun -n 1 ./nbody.out -N 4096   >& ${DIR}/00stdout-2_12-${id}th.log
   time aprun -n 1 ./nbody.out -N 8192   >& ${DIR}/00stdout-2_13-${id}th.log
   time aprun -n 1 ./nbody.out -N 16384  >& ${DIR}/00stdout-2_14-${id}th.log
   time aprun -n 1 ./nbody.out -N 32768  >& ${DIR}/00stdout-2_15-${id}th.log
   time aprun -n 1 ./nbody.out -N 65536  >& ${DIR}/00stdout-2_16-${id}th.log
   time aprun -n 1 ./nbody.out -N 131072 >& ${DIR}/00stdout-2_17-${id}th.log
   time aprun -n 1 ./nbody.out -N 262144 >& ${DIR}/00stdout-2_18-${id}th.log
done

