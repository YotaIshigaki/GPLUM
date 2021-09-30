#!/bin/sh
cat le_r_0*txt > le_r.txt
sort -g le_r.txt > le_r_sorted.txt
uniq le_r_sorted.txt le_r_final.txt
