#!/usr/bin/awk
{
    # In the following, we assume that the input file for which 
    # this awk script is applied is created by the following command:
    # $ grep "s/step" 00stdout-2_10-* | awk '{print $4}'
    x[NR] = $0
}
END {
    if (NR == 0) exit

    max = -3.402823e38
    min =  3.402823e38
    sum = 0.0e0
    for (i in x) {
      if (x[i] > max) max = x[i]
      if (x[i] < min) min = x[i]
      sum += x[i]
    }
    ave = sum / NR

    d2 = 0.0e0
    for (i in x) {
      d2 += (x[i] - ave) ^ 2
    }
    sigma = sqrt(d2 / NR)

    # Output the stat.
    printf "---------------------------------------\n"
    printf "Num.           = %d\n", NR
    printf "Max.           = %.15f\n", max
    printf "Min.           = %.15f\n", min
    printf "Ave.           = %.15f\n", ave
    printf "Std.dev.       = %.15f\n", sigma
    printf "Std.dev.(pct.) = %.15f\n", sigma/ave*1.0e2
    printf "---------------------------------------\n"
}
