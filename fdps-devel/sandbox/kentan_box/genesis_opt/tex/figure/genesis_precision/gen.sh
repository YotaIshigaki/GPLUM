sort -n -k2 genesis_single.dat > tmp1
sort -n -k2 genesis_double.dat > tmp2

paste direct.dat tmp1 tmp2 poly.dat | awk 'BEGIN{printf "# r2 direct single double poly\n"}{printf "%.16e %.16e %.16e %.16e %.16e\n",$3,$4,$9,$14,$18}' > diff.dat
