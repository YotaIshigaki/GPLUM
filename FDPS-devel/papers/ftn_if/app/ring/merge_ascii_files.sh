#!/bin/bash
# Set the range of id of snapshot files
i_start=0
i_end=9
i_offset=1
# Merge separated-snapshot files into a single file
str_back='-proc*.dat'
for ((i=${i_start}; i <= ${i_end}; i += ${i_offset}))
do
  str_front=$(printf "snap%05d" ${i})
  files="${str_front}${str_back}"
  cat ${files} > ${str_front}".dat"
done
