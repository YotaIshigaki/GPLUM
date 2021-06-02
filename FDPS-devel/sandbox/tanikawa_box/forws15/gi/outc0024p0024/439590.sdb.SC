#PBS -N test
#PBS -l mppwidth=24
#PBS -j oe
#PBS -q short-b

NPARALLEL=$NCPUS
NPROCESS=$NCPUS
#NPARALLEL=24
#NPROCESS=24

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir result
cp $0 result

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./sph.out

ncore=`printf "%04d" $NPARALLEL`
nproc=`printf "%04d" $NPROCESS`
mv result outc"$ncore"p"$nproc"
