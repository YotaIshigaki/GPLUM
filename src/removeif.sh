#!/bin/bash

function showArgs() {
    for i in $@
    do
        j=`echo ${i} | cut -c 1-2`
        if [ "${j}" = "-D" ]
        then
            echo ${i} | cut -c 3-
        fi
    done
}

function searchArgs() {
    f=-1
    for i in $@
    do
        if [ "${i}" = "$1" ]
        then
            f=`expr ${f} + 1`
        fi
    done
    echo ${f}
}

function multipleArgs() {
    product=1
    for i in $@
    do
        product=`expr ${product} \* ${i}`
    done
    echo ${product}
}

keywords="#ifdef #ifndef #else #endif"
flags=(1 1 1 1 1 1 1 1)
level=0

args=`showArgs $@`

cat $1 | while read line
do
    head=`echo ${line} | cut -d " " -f 1`
    arg=`echo ${line} | cut -d " " -f 2`
    if [ `searchArgs ${head} ${keywords}` -ge 1 ]
    then
        if [ "${head}" = "#ifdef" ]
        then
            level=`expr ${level} + 1`
            if [ `searchArgs ${arg} ${args}` -ge 1 ]
            then
                flags["${level}"]=1
            else
                flags["${level}"]=0
            fi
        elif [ "${head}" = "#ifndef" ]
        then
            level=`expr ${level} + 1`
            if [ `searchArgs ${arg} ${args}` -ge 1 ]
            then
                flags["${level}"]=0
            else
                flags["${level}"]=1
            fi
        elif [ "${head}" = "#else" ]
        then
            if [ ${flags["${level}"]} -gt 0 ]
            then
                flags["${level}"]=0
            else
                flags["${level}"]=1
            fi
        elif [ "${head}" = "#endif" ]
        then
            flags["${level}"]=1
            level=`expr ${level} - 1`
        fi
    else
        if [ `multipleArgs ${flags[@]}` -gt 0 ]
        then
            echo "${line}"
        fi
    fi
done > $2
