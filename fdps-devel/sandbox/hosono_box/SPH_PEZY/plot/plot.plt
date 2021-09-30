reset
set size square
set term gif enhanced size 512, 512
set out sprintf("img/%05d.gif", i);
p [-3:3][-3:3]sprintf("< cat ../src/result/%05d_*.dat", i) u 3:4 pt 7 ps 0.5


reset
