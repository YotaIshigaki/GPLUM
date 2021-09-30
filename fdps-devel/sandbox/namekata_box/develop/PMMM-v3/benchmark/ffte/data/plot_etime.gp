set terminal postscript enhanced color eps "Palatino-Bold,24"
set size square
set grid
set xlabel "# of processes"
set ylabel "Wall clock time (fwd+bkw) [s]"
set logscale xy
set yrange [1.0e-4:1.0e2]
set format x "10^{%L}"
set format y "10^{%L}"
set output "etime.eps"
plot "etime.txt" index 0 u 1:2 w linespoints ps 3 pt 4 lw 5 lc 2 t "NC=128^{3}", \
     "" index 1 u 1:2 w linespoints ps 3 pt 4 lw 5 lc 3 t "NC=256^{3}", \
     "" index 2 u 1:2 w linespoints ps 3 pt 4 lw 5 lc 4 t "NC=512^{3}"
