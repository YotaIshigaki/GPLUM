#set terminal postscript enhanced color eps "Palatino-Bold,22"
set size square
set grid
set xlabel "# of particles per core"
set ylabel "wall clock time per steps [s]"
set format x "10^{%L}"
set format y "10^{%L}"
set format y2 "10^{%L}"
set yrange  [1.0e-3:1.0e1]
set y2range [1.0e-5:1.0e0]
set ytics nomirror
set y2tics
set logscale xy
set logscale y2
set key right bottom
#set output "nbody_comp.eps"
#### Start multiplot
set multiplot layout 1,3 
#----------------------
# K-computer
#----------------------
set title "K-Computer"
plot "nbody_data_serial.txt" index 2 using 1:2 w lp ps 2 pt 5 lw 2 lc rgb "#ff0000" t "Fortran" axes x1y1, \
     "nbody_data_serial.txt" index 2 using 1:3 w lp ps 2 pt 3 lw 2 lc rgb "#0000ff" t "C++" axes x1y1, \
     "nbody_data_serial.txt" index 2 using 1:(abs(($2)-($3))/(0.5*(($2)+($3)))) w lp ps 1 pt 7 lc rgb "#000000" t "Rel. diff." axes x1y2
#----------------------
# Oakforest-PACS
#----------------------
set title "Oakforest-PACS"
unset ylabel
plot "nbody_data_serial.txt" index 1 using 1:2 w lp ps 2 pt 5 lw 2 lc rgb "#ff0000" t "Fortran" axes x1y1, \
     "nbody_data_serial.txt" index 1 using 1:3 w lp ps 2 pt 3 lw 2 lc rgb "#0000ff" t "C++" axes x1y1, \
     "nbody_data_serial.txt" index 1 using 1:(abs(($2)-($3))/(0.5*(($2)+($3)))) w lp ps 1 pt 7 lc rgb "#000000" t "Rel. diff." axes x1y2
#----------------------
# knl (for reference)
#----------------------
#set title "knl"
#plot "nbody_data_serial.txt" index 0 using 1:2 w lp ps 2 pt 5 lw 2 lc rgb "#ff0000" t "Fortran" axes x1y1, \
#     "nbody_data_serial.txt" index 0 using 1:3 w lp ps 2 pt 3 lw 2 lc rgb "#0000ff" t "C++" axes x1y1, \
#     "nbody_data_serial.txt" index 0 using 1:(abs(($2)-($3))/(0.5*(($2)+($3)))) w lp ps 1 pt 7 lc rgb "#000000" t "Rel. diff." axes x1y2
#----------------------
# ajisai
#----------------------
set title "Ajisai"
set y2label "Relative difference"
plot "nbody_data_serial.txt" index 3 using 1:2 w lp ps 2 pt 5 lw 2 lc rgb "#ff0000" t "Fortran" axes x1y1, \
     "nbody_data_serial.txt" index 3 using 1:3 w lp ps 2 pt 3 lw 2 lc rgb "#0000ff" t "C++" axes x1y1, \
     "nbody_data_serial.txt" index 3 using 1:(abs(($2)-($3))/(0.5*(($2)+($3)))) w lp ps 1 pt 7 lc rgb "#000000" t "Rel. diff." axes x1y2
#----------------------
# Cray XC30 (CfCA)
#----------------------
#set title "Cray XC30 (CfCA)"
#plot "nbody_data_serial.txt" index 4 using 1:2 w lp ps 2 pt 5 lw 2 lc rgb "#ff0000" t "Fortran" axes x1y1, \
#     "nbody_data_serial.txt" index 4 using 1:3 w lp ps 2 pt 3 lw 2 lc rgb "#0000ff" t "C++" axes x1y1, \
#     "nbody_data_serial.txt" index 4 using 1:(abs(($2)-($3))/(0.5*(($2)+($3)))) w lp ps 1 pt 7 lc rgb "#000000" t "Rel. diff." axes x1y2

### End multiplot
unset multiplot
