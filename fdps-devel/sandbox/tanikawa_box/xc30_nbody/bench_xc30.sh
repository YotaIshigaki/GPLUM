#PBS -N test
#PBS -l mppwidth=24
#PBS -j oe
#PBS -q short-b

NPARALLEL=24
NPROCESS=24
npernode=6291456
ngroup=512

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`
ncore=`printf "%04d" $NPARALLEL`
nproc=`printf "%04d" $NPROCESS`
odir=outc"$ncore"p"$nproc"
rm -rf $odir
mkdir $odir
cp $0 $odir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
echo "nptcl:  $nptcl"   > "$odir"/runparam.log
echo "ngroup: $ngroup" >> "$odir"/runparam.log
aprun -n $NPROCESS -d $OMP_NUM_THREADS ./test_nbody.out -N $nptcl -n $ngroup -T 0.03125 -t 0.4 -o $odir -s 200
