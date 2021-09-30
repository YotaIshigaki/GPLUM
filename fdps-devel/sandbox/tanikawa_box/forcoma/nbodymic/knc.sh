#!/bin/bash
#SBATCH -J test_native
#SBATCH -p mic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 30
#SBATCH -o stdout.log
#SBATCH -e stderr.log

module load intel intelmpi
export I_MPI_MIC=1
export MIC_PPN=15
cd $SLURM_SUBMIT_DIR

mpirun-mic -m ./nbody.out
#mpirun-mic -m ./nbody.out -N 262144
#mpirun-mic2 -m ./run
