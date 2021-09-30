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
set output "density_reference.eps"
plot "../result.ref/sph00001.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#000030" t "t = 0.112", \
     "../result.ref/sph00002.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#000090" t "t = 0.223", \
     "../result.ref/sph00003.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#000fff" t "t = 0.334", \
     "../result.ref/sph00004.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#0090ff" t "t = 0.445", \
     "../result.ref/sph00005.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#0fffee" t "t = 0.556", \
     "../result.ref/sph00006.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#90ff70" t "t = 0.667", \
     "../result.ref/sph00007.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#ffee00" t "t = 0.778", \
     "../result.ref/sph00008.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#ff7000" t "t = 0.889", \
     "../result.ref/sph00009.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):9 w d lc rgb "#ee0000" t "t = 0.999"
#=======================
# Plot dt(r)
#=======================
set xlabel "r"
set ylabel "Timestep"
set key font "Palatino-Bold, 16"
set key left top
set output "timestep_reference.eps"
plot "../result.ref/sph00001.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#000030" t "t = 0.112", \
     "../result.ref/sph00002.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#000090" t "t = 0.223", \
     "../result.ref/sph00003.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#000fff" t "t = 0.334", \
     "../result.ref/sph00004.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#0090ff" t "t = 0.445", \
     "../result.ref/sph00005.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#0fffee" t "t = 0.556", \
     "../result.ref/sph00006.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#90ff70" t "t = 0.667", \
     "../result.ref/sph00007.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#ffee00" t "t = 0.778", \
     "../result.ref/sph00008.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#ff7000" t "t = 0.889", \
     "../result.ref/sph00009.txt" every ::2 u (sqrt($3*$3+$4*$4+$5*$5)):13 w d lc rgb "#ee0000" t "t = 0.999"
