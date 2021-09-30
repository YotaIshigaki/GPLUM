set terminal postscript enhanced color eps "Palatino-Bold,24"
set size square
set grid
set xlabel "# of particles per process"
set ylabel "Wall clock time [s]"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set output "etime.eps"
f(x)=2.0e-6 * x * log(x)
plot "etime.K.txt" u 1:2 w linespoints ps 3 pt 4 lw 5 lc 0 notitle, \
     f(x) w l lw 3 dt 3 lc rgb "#ff0000" title "{/Symbol \265} N log(N)"
