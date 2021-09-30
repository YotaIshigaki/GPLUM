if test $# -ne 2
then
    echo "$0 <nproc> <ifile>"
    exit
fi

nproc=$1
ifile=$2

mpirun -np $nproc ./run $ifile

awk -f awk/make_dd.awk out/domain_all.txt > out/dd.txt 

awk -f awk/mapDensity.awk xmin=-8 xmax=6 ymin=-8 ymax=6 width=0.0625 out/systeme_00* \
    > out/mapDensity.dat
