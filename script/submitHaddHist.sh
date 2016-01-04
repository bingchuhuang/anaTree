#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./submitHaddHist.sh <outName>"
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
if [ ! -d ${dir}/out/out_$1_hist ]; then
	mkdir ${dir}/out/out_$1_hist
fi
if [ ! -d ${dir}/log/log_$1_hist ]; then
	mkdir ${dir}/log/log_$1_hist
fi

echo $dir
cp run.con runHaddAll_$1.con
echo "Executable       = haddHistByRun.sh">>runHaddAll_$1.con
for runId in `cat runNumber_tobeHadd_$1`
do
      echo "Arguments        = $1 $runId ">>runHaddAll_$1.con
      echo "Output           = log/log_$1_hist/run_$runId.out">>runHaddAll_$1.con
      echo "Error            = log/log_$1_hist/run_$runId.err">>runHaddAll_$1.con
      echo "Log              = log/log_$1_hist/run_$runId.olog">>runHaddAll_$1.con
      echo  "Queue" >>runHaddAll_$1.con
      echo  "     " >>runHaddAll_$1.con
done
   
condor_submit runHaddAll_$1.con
