set logscale
set xrange[0.95:16.5]
set term pngcairo
set out "time.png"
set ylabel "time [sec]"
set xlabel "# of process"
p "log_double.txt" u 1:(($2 + $3 + $4 + $5) * 3) pt 6 lc rgb "red" t "PEZY-SC (double) x 3",\
  "log_float.txt" u 1:(($2 + $3 + $4 + $5) * 27) pt 4 lc rgb "dark-green" t "PEZY-SC (float) x 27",\
  "log_CPU.txt" u 1:($2 + $3 + $4 + $5) pt 8 lc rgb "blue" t "Xeon E5-2618L (8 cores HT)"#,\
  #"log_PG.txt" u 1:(($2 + $3 + $4 + $5)) pt 2 lc rgb "orange" t "Phantom Grape"
