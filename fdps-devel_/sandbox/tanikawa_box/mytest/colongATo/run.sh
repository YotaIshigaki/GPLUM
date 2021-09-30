if test $# -ne 2
then
    echo "$0 <nproc> <ifile>"
    exit
fi

nproc=$1
ifile=$2

mpirun-openmpi-gcc49 -np $nproc ./run $ifile
cat dump??.log | sort -n > dump.log
rm -f dump??.log
