#!/bin/sh

st=$1
n_proc=`getconf _NPROCESSORS_ONLN`
n_proc=`expr $n_proc - 10`
n_proc=2

while [ $st -le $2 ]; do
	core=0
	while [ $core -le $n_proc ]; do
		i=`expr $core + $st`
		if [ $i -gt $2 ]; then
			break
		fi
		if [ ! -e $(printf "img/%05d.gif" $i) ]; then
			gnuplot -e "i=${i}" plot.plt &
		fi
		core=`expr $core + 1`
	done
	wait
	st=`expr $st + $n_proc`
done

convert -loop 0 img/*.gif anim.gif
cp anim.gif ~/
