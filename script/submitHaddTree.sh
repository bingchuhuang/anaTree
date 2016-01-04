#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./submitHaddTree.sh <outName>"
	exit
fi

dir=$(pwd)
echo $dir
if [ ! -d ${dir}/out ]; then
    mkdir ${dir}/out
fi
if [ ! -d ${dir}/log ]; then
    mkdir ${dir}/log
fi


if [ ! -d ${dir}/out/out_$1_anaTree ]; then
	mkdir ${dir}/out/out_$1_anaTree
fi
if [ ! -d ${dir}/log/log_$1_anaTree ]; then
	mkdir ${dir}/log/log_$1_anaTree
fi

echo $dir
cp run.con runHaddAllTree_$1.con
echo "Executable       = haddTreeByRun.sh">>runHaddAllTree_$1.con
for runId in `ls -1 /star/u/huangbc/data01/pico/anaTree/20151125_v4/prod/out/out_$1`
#for runId in `cat runNumber_tobeHadd_$1`
do
      echo $runId
      echo "Arguments        = $1 $runId ">>runHaddAllTree_$1.con
      echo "Output           = log/log_$1_anaTree/run_$runId.out">>runHaddAllTree_$1.con
      echo "Error            = log/log_$1_anaTree/run_$runId.err">>runHaddAllTree_$1.con
      echo "Log              = log/log_$1_anaTree/run_$runId.olog">>runHaddAllTree_$1.con
      echo  "Queue" >>runHaddAllTree_$1.con
      echo  "     " >>runHaddAllTree_$1.con
done
   
condor_submit runHaddAllTree_$1.con
