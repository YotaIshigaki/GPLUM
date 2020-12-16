if test $# -ne 1
then
    echo "$0 <nproc>"
    exit
fi

nproc=$1

mpirun-openmpi-gcc49 -np $nproc ./run
