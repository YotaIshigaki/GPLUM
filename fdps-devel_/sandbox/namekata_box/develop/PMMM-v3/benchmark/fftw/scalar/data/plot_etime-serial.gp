set terminal postscript enhanced color eps "Palatino-Bold,24"
set size square
set grid
set xlabel "# of cells"
set ylabel "Wall clock time (fwd+bkw) [s]"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set output "etime-serial.eps"
f(x)=1.0e-8 * x * log(x)
plot "etime-serial.txt" u 2:3 w linespoints ps 3 pt 4 lw 5 lc 0 notitle, \
     f(x) w l lw 3 dt 3 lc rgb "#ff0000" title "{/Symbol \265} N log(N)"
