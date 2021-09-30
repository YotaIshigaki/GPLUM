#!/bin/bash
# Set parameters
TARGET_DIR=$1
FILE="result.txt"
NPTCL_MAX=18
# Analyze
echo "Start the analysis of ${TARGET_DIR} directory."
cd ${TARGET_DIR}
rm -f ${FILE}
touch ${FILE}
for ((i=10; i<= ${NPTCL_MAX}; i++)); do
    echo "Stat. for 2^{${i}} particles" >> ${FILE}
    grep "s/step" 00stdout-2_${i}-* | awk '{print $3}' | awk -f ../get_stat.awk >> ${FILE}
    echo " " >> ${FILE}
done
