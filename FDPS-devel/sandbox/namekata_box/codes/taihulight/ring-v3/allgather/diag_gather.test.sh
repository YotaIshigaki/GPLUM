#!/bin/sh

for nrad in `seq 1 5`
do
  for nphi in `seq 1 10`
  do
    if [ $nrad -lt $nphi ]; then
      np=`expr $nrad \* $nphi`
      echo $nrad $nphi $np
      mpirun-openmpi-gcc5 -n $np ./a.out $nrad $nphi > /dev/null
	fi
  done
done
