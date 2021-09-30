set terminal postscript enhanced color eps
set size square
# [1] W
set xrange [0.0:2.5]
set xlabel 'r'
set ylabel 'W'
set output 'W.eps'
plot 'W_M4.txt' u 1:2 w l lw 2 title 'cubic spline', \
     'W_M5.txt' u 1:2 w l lw 2 title 'quartic spline', \
     'W_M6.txt' u 1:2 w l lw 2 title 'quintic spline', \
     'W_CT.txt' u 1:2 w l lw 2 title 'core triangle'
