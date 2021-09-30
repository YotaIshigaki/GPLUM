set terminal postscript enhanced color eps "Palatino-Bold, 24"
set size square
set logscale xy
#=======================
# Plot \rho(r)
#=======================
set xlabel "r"
set ylabel "Density"
set key font "Palatino-Bold, 16"
set key left bottom
set output "density_new_method_v2.eps"
plot "../result/sph00001.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#000030" t "t = 0.112", \
     "../result/sph00002.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#000090" t "t = 0.223", \
     "../result/sph00003.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#000fff" t "t = 0.334", \
     "../result/sph00004.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#0090ff" t "t = 0.445", \
     "../result/sph00005.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#0fffee" t "t = 0.556", \
     "../result/sph00006.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#90ff70" t "t = 0.667", \
     "../result/sph00007.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#ffee00" t "t = 0.778", \
     "../result/sph00008.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#ff7000" t "t = 0.889", \
     "../result/sph00009.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#ee0000" t "t = 0.999"
