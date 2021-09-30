#PBS -N nbody.ftn
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
   echo 1024 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_10-${id}th.log
   echo 2048 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_11-${id}th.log
   echo 4096 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_12-${id}th.log
   echo 8192 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_13-${id}th.log
   echo 16384 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_14-${id}th.log
   echo 32768 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_15-${id}th.log
   echo 65536 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_16-${id}th.log
   echo 131072 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_17-${id}th.log
   echo 262144 > nptcl.input
   time aprun -n 1 ./nbody.out >& ${DIR}/00stdout-2_18-${id}th.log
done
