#!/bin/sh
cat mm_r_0*txt > mm_r.txt
sort -g mm_r.txt > mm_r_sorted.txt
uniq mm_r_sorted.txt mm_r_final.txt
