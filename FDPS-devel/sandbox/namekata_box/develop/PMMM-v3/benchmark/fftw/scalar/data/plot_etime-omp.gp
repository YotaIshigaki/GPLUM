set terminal postscript enhanced color eps "Palatino-Bold,24"
set size square
set grid
set xlabel "# of cells"
set ylabel "Wall clock time (fwd+bkw) [s]"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set output "etime-omp-singleFFT.eps"
plot "etime-omp.txt" index 0 u 2:3 w linespoints ps 3 pt 4 lw 5 lc 0 t "2 threads", \
     "" index 0 u 2:4 w linespoints ps 3 pt 4 lw 5 lc 1 t "4 threads",\
     "" index 0 u 2:5 w linespoints ps 3 pt 4 lw 5 lc 2 t "8 threads"

set output "etime-omp-multiFFTs.eps"
plot "etime-omp.txt" index 1 u 2:3 w linespoints ps 3 pt 4 lw 5 lc 0 t "1 threads", \
     "" index 1 u 2:4 w linespoints ps 3 pt 4 lw 5 lc 1 t "2 threads",\
     "" index 1 u 2:5 w linespoints ps 3 pt 4 lw 5 lc 2 t "4 threads", \
     "" index 1 u 2:6 w linespoints ps 3 pt 4 lw 5 lc 3 t "8 threads"
