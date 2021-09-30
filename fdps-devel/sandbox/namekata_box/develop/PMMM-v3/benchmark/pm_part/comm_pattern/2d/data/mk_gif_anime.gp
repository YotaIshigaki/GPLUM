set terminal gif animate delay 5 optimize size 640,480
set output "comm_pattern.gif"
set size square
set key outside
set xrange [-4:12]
set yrange [-4:12]
do for [i=0:63] {
  plot "comm_pattern_2d.txt" u 1:2 w p pt 5 lc rgb "#000000" title "Rank", \
       "comm_pattern_2d.txt" index i u 1:2:3:4 w vector lc rgb "#ff0000" title "direction of MPI\_Send"
}
