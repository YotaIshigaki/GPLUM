NPROCESS=1
export OMP_NUM_THREADS=1
npernode=6291456

nproc=`printf "%04d" $NPROCESS`
nthrd=`printf "%02d" $OMP_NUM_THREADS`
odir=outp"$nproc"t"$nthrd"
rm -rf $odir
mkdir $odir
cp $0 $odir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
echo "$nptcl" > "$odir"/runparam.log
mpirun -n $NPROCESS ./test_nbody.out -N $nptcl -n 512 -T 0.015625 -t 0.4 -o $odir
#mpirun -n $NPROCESS ./test_nbody.out -N $nptcl -n 384 -T 0.015625 -t 0.4 -o $odir
#mpirun -n $NPROCESS ./test_nbody.out -N $nptcl -n 256 -T 0.015625 -t 0.4 -o $odir
