# Reset the previous setting if exists
reset
# Plot
proc_num = 0
i_start  = 0
i_end    = 100
i_offset = 1
do for [i = i_start : i_end : i_offset] {
    input = sprintf("snap%05d-proc%05d.dat",i,proc_num)
    set xrange [-6.0e3:6.0e3]
    set yrange [-6.0e3:6.0e3]
    set zrange [-6.0e3:6.0e3]
    splot input u 3:4:5 w p
    pause 0.1
}
# Reset the current setting for the safety
reset
