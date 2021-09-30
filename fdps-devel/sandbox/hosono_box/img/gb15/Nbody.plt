reset

set logscale

#unset key

set term postscript eps enhanced color size 3., 6.
set out "Flops.eps"
set multiplot layout 2, 1
set size square

set format x "10^{%T}"
set format y "10^{%+-T}"
set tics font ", 21pt"

set key bottom right font ", 21pt"



set ylabel "TFlops" font ", 21pt" offset -1, 0
p [][:] "K/SPH_Flops_Wtime.dat" u 1:2 w lp pt 5 ps 1.5 lc rgb "blue" t "SPH K",\
        "XC30/SPH_Flops_Wtime.dat" u 1:2 w lp pt 5 ps 1.5 lc rgb "red" t "SPH XC30",\
        "K/disk_bench_flops.log" u ($1*8):2 w lp pt 7 ps 1.5  lc rgb "blue" t "Nbody K (disk)",\
        "K/pl_bench_flops.log" u ($1*8):2 w lp pt 6 ps 1.5  lc rgb "blue" t "Nbody K (Plummer)",\
        "XC30/Nbody_Flops_Wtime.dat" u 1:2 w lp pt 7 ps 1.5 lc rgb "red" t "Nbody XC30 (Plummer)"
        
             

set xlabel "# of cores" font ", 21pt"
set ylabel "Wallclock time [sec]" offset -1, 0
p [][1:30] "K/SPH_Flops_Wtime.dat" u 1:3 w lp pt 5 ps 1.5 lc rgb "blue" t "SPH K",\
           "XC30/SPH_Flops_Wtime.dat" u 1:3 w lp pt 5 ps 1.5  lc rgb "red" t "SPH XC30",\
           "K/disk2_bench_time.log" u ($1*8):7 w lp pt 7 ps 1.5 lc rgb "blue" t "Nbody K (disk)",\
           "K/pl2_bench_time.log" u ($1*8):7 w lp pt 6 ps 1.5 lc rgb "blue" t "Nbody K (Plummer)",\
           "XC30/Nbody_Flops_Wtime.dat" u 1:3 w lp pt 7 ps 1.5  lc rgb "red" t "Nbody XC30 (Plummer)"

reset
