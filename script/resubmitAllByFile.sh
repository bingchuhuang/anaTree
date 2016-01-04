#!/bin/bash
date

if [ $# -lt 2 ]; then
    echo "not enough arguements"
    echo "./resubmitAllByFile.sh <outName> <mode: 0-mb, 1-ht, 2-mtd>"
	exit
fi

dir=$(pwd)
ver=$1
echo $dir
if [ ! -d ${dir}/out/out_$1 ]; then
   echo "no out found! exit!"
   exit
fi

if [ ! -d submit ]; then
   echo "no submit found! exit!"
   exit
fi

if [ ! -d submit/$1 ]; then
   echo "no submit/$1 found! exit!"
   exit
fi

#./findResubmitRun.sh $1

echo $dir
for runId in `cat runNumber_tobeResubmit_$1`
do
   if [ ! -d $dir/submit/$1/$runId ]; then
      mkdir $dir/submit/$1/$runId
   fi
   if [ ! -f $dir/submit/$1/$runId/AuAu200.xml ]; then
      cp $dir/AuAu200.xml $dir/submit/$1/$runId/.
   fi

   cd $dir/submit/$1/$runId/
   if [ ! -d ${dir}/out/out_$1/$runId ]; then
   	mkdir ${dir}/out/out_$1/$runId
   fi
   if [ ! -d ${dir}/log/log_$1/$runId ]; then
   	mkdir ${dir}/log/log_$1/$runId
   fi

   sesName=`basename ${dir}/submit/$1/${runId}/*.session.xml | cut -d . -f1`
   nJobs=`find ${dir}/submit/$1/${runId}/sched${sesName}_*.list | wc -l`
   for (( ijob=0; ijob<$nJobs; ijob++ ))
   do
   str=`cat ${dir}/submit/$1/${runId}/sched${sesName}_${ijob}.list`
   muDst=${str##*/}
   bName=${muDst%%.*}
   anaTree=${dir}/out/out_$1/${runId}/${bName}.anaTree.root
   if [ ! -f $anaTree ]; then
      echo "resubmitting $runIdId job Id $ijob"
      star-submit -r $ijob ${sesName}.session.xml
   fi
   done
done
