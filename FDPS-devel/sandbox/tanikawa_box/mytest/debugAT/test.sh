if test $# -ne 1
then
    echo "$0 <nproc>"
    exit
fi

for system in system0 system1
do
    echo "!!!!!!!!!!!!!!!!!! SYSTEM $system TEST !!!!!!!!!!!!!!!!!!"

    xmin=`awk '{if(NR==1) print $1}' out/domain_all.txt`
    ymin=`awk '{if(NR==1) print $2}' out/domain_all.txt`
    zmin=`awk '{if(NR==1) print $3}' out/domain_all.txt`
    xmax=`awk '{if(NR==2) print $1}' out/domain_all.txt`
    ymax=`awk '{if(NR==2) print $2}' out/domain_all.txt`
    zmax=`awk '{if(NR==2) print $3}' out/domain_all.txt`
    nerror=`awk '{if((xmin<=$1&&$1<=xmax)&&(ymin<=$2&&$2<=ymax)&&(zmin<=$3&&$3<=zmax)){;}else{print xmin, xmax, ymin, ymax, zmin, zmax;}}' \
        xmin=$xmin xmax=$xmax ymin=$ymin ymax=$ymax zmin=$zmin zmax=$zmax out/"$system"i_????.txt | wc -l | awk '{print $1}'`
    if test $nerror -eq 0
    then
        echo "test domain is OK."
    else
        echo "test domain is NG."
    fi

    nerror=0
    for iproc in $(seq -f "%04g" 1 1 $nproc)
    do
        nerrort=`diff out/"$system"i_$iproc.txt out/"$system"f_$iproc.txt | wc -l | awk '{print $1}'`
        nerror=`expr $nerror + $nerrort`
    done

    if test $nerror -eq 0
    then
        echo "test recover psys is OK."
    else
        echo "test recover psys is NG."
    fi
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"    
    
done

#cat out/psmp_000[01].txt \
#    | sort -g \
#    > hoge.out
#nerror=`diff out/psmp_all.txt hoge.out | wc -l | awk '{print $1}'`
#if test $nerror -eq 0
#then
#    echo "collecting sample on root process is OK."
#else
#    echo "collecting sample on root process is NG."
#fi

nsample=`cat out/psmp_????.txt | wc -l | awk '{print $1}'`
echo "Number of sample particles is $nsample"

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"    

nerror=0
for iproc in $(seq -f "%04g" 1 1 $nproc)
do
    nerrort=`diff out/posdomain_0000.txt out/posdomain_$iproc.txt | wc -l | awk '{print $1}'`
    nerror=`expr $nerror + $nerrort`
done
if test $nerror -eq 0
then
    echo "test Bcast pos_domain is OK."
else
    echo "test Bcast pos_domain is NG."
fi
