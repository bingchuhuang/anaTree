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
 
if [ -f runNumber_$1_$2 ]; then
   rm runNumber_$1_$2
fi

get_file_list.pl -keys 'runnumber' -cond 'production=P15ic,collision=auau200,trgsetupname=AuAu_200_production_2014||AuAu_200_production_mid_2014||AuAu_200_production_low_2014,filename~st_physics,filetype=daq_reco_Mudst,storage=local' -limit 0 >& runNumber_$1_$2
cp runNumber_$1_$2 ../. 

echo $dir
for runId in `cat runNumber_$1_$2`
do
   if [ ! -d $dir/submit/$runId ]; then
      mkdir $dir/submit/$runId
   fi
   cp $dir/AuAu200.xml $dir/submit/$runId/.
   cd $dir/submit/$runId/
   if [ ! -d ${dir}/out/out_$1/$runId ]; then
   	mkdir ${dir}/out/out_$1/$runId
   fi
   if [ ! -d ${dir}/log/log_$1/$runId ]; then
   	mkdir ${dir}/log/log_$1/$runId
   fi
   star-submit-template -template AuAu200.xml -entities path=${dir},ver=$1,mode=$2,run=$runId
done
