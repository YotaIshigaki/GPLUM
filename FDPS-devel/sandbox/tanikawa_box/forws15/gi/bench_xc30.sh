#PBS -N test
#PBS -l mppwidth=1536
#PBS -j oe
#PBS -q short-b

NPARALLEL=1536
NPROCESS=$NPARALLEL

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir result
cp $0 result

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./sph.out

ncore=`printf "%04d" $NPARALLEL`
nproc=`printf "%04d" $NPROCESS`
mv result outc"$ncore"p"$nproc"
