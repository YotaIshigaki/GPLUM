#!/bin/bash
istart=0
iend=110
for ((i=$istart; i < $iend; i++)); do
    n=`printf "%05d" $i`
    diff intlist_else_$n.txt intlist_func_$n.txt
done
