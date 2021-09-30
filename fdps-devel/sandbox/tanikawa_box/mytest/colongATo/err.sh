if test $# -ne 1
then
    echo "sh $0 <ifile>"
    exit
fi

ifile=$1

devi=6

paste $ifile dump.log \
    | awk '{print $2-$(2+devi), $3-$(3+devi), $4-$(4+devi), $2, $3, $4}' devi=$devi \
    | awk '{acc1=$1*$1+$2*$2+$3*$3;acc2=$4*$4+$5*$5+$6*$6; if(acc1==0.&&acc2==0.){print 0.;}else{print sqrt((acc1)/(acc2))}}' \
    | sort -g \
    | awk '{print $1, NR/8192.}' \
    > acc.log

paste $ifile dump.log \
    | awk '{if($5==0.&&$(5+devi)==0.){print 0.}else{print ($5-$(5+devi))/$(5+devi)}}' devi=$devi \
    | awk '{print sqrt($1*$1)}' \
    | sort -g \
    | awk '{print $1, NR/8192.}' \
    > pot.log

