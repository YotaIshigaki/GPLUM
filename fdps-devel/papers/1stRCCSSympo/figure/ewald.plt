set term post enh color eps 24
set size square

set xrange [0.0:10]
set yrange [:1.5]

unset xtics
unset ytics

set output "direct.eps"
p 1/x lw 8 not

set output "long.eps"
p erf(x)/x lw 8 not

set output "short.eps"
p erfc(x)/x lw 8 not

set xrange [:2.5]
set output "short_close.eps"
p erfc(x)/x lw 8 not

