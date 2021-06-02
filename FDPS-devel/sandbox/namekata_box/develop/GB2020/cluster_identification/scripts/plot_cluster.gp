set terminal x11
set xrange [400:500]
splot 'cluster_data' u 1:2:3,'' u 4:5:6,'' u 1:2:3:($4-$1):($5-$2):($6-$3) w vector,'pos_data.txt' w d lc rgb "#000000"
