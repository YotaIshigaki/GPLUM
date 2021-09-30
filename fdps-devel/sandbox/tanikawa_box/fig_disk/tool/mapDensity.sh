if test $# -ne 2
then
    echo "sh $0 <ifile> <ofile>"
    exit
fi

ifile=$1
ofile=$2

awk -f awk/mapDensity.awk xmin=-40 xmax=40 ymin=-40 ymax=40 width=0.125 $ifile \
    > $ofile
