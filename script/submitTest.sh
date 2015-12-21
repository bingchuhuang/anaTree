#!/bin/bash
date

if [ $# -lt 2 ]; then
    echo "not enough arguements"
    echo "./submitAll.sh <outName> <mode: 0-mb, 1-ht, 2-mtd>"
	exit
fi

if [ ! -d out ]; then
    mkdir out
fi
if [ ! -d log ]; then
    mkdir log
fi

dir=$(pwd)
ver=$1
echo $dir
if [ ! -d ${dir}/out/out_$1 ]; then
	mkdir ${dir}/out/out_$1
fi
if [ ! -d ${dir}/log/log_$1 ]; then
	mkdir ${dir}/log/log_$1
fi

if [ ! -d submit ]; then
   mkdir submit
fi

cd submit/
 
 
echo $dir

runId=15094070
if [ ! -d ${dir}/submit/$runId ]; then
   mkdir ${dir}/submit/$runId
fi
cp $dir/AuAu200.xml $dir/submit/$runId/.
cd $dir/submit/$runId/
if [ ! -d ${dir}/out/out_$1/$runId ]; then
	mkdir ${dir}/out/out_$1/$runId
fi
if [ ! -d ${dir}/log/log_$1/$runId ]; then
	mkdir ${dir}/log/log_$1/$runId
fi
star-submit-template -template AuAu200.xml -entities path=${dir},ver=$1,mode=$2,run=${runId}
