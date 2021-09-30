set terminal postscript enhanced color eps "Palatino-Bold,24"
set size square
set grid
set xlabel "Data size [B]"
set ylabel "Wall clock time [s]"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set output "etime.eps"
plot "etime.txt" u 1:2 w linespoints ps 3 pt 4 lw 5 lc 0 notitle
