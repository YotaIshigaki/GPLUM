if test $# -ne 1
then
    echo "sh $0 <file>"
    exit
fi

echo "start gnuplot"
gnuplot $1.gpl
echo "end gnuplot"
awk -f $1.awk tmp.eps > tmp2.eps
mv tmp2.eps $1.eps
rm tmp.eps

echo "start convert eps to png"
#convert -density 288 $1.eps $1.png
convert -density 144 $1.eps $1.png
echo "end convert eps to png"

echo "start convert png to eps"
convert $1.png $1.eps
