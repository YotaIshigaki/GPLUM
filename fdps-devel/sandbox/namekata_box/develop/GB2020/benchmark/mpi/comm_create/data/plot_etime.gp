set terminal postscript enhanced color eps "Palatino-Bold,24"
set size square
set grid
set title "Total: 65536 processes"
set xlabel "# of processes in a new comm."
set ylabel "Wall clock time [s]"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set key left top
set output "etime_comm_create.eps"
plot "etime.txt" index 0 u 1:2 w linespoints ps 3 pt 4 lw 5 lc 0 t "Single comm", \
     "" index 1 u 1:2 w linespoints ps 3 pt 4 lw 5 lc rgb "#ff0000" t "Multiple comm"

set output "etime_allreduce.eps"
plot "etime.txt" index 2 u 1:2 w linespoints ps 3 pt 4 lw 5 lc 0 t "Allreduce"
