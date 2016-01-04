#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./findResubmitFile.sh <outName>"
	exit
fi

#resubmit = no anaTree + evicted 
for run in `cat submit/runNumber_$1_0`
do
    sesName=`basename submit/$1/${run}/*.session.xml | cut -d . -f1`
    nJobs=`find submit/$1/${run}/sched${sesName}_*.list | wc -l`
    for (( ijob=0; ijob<$nJobs; ijob++ ))
    do
    str=`cat submit/$1/${run}/sched${sesName}_${ijob}.list`
    muDst=${str##*/}
    bName=${muDst%%.*}
    anaTree=out/out_$1/${run}/${bName}.anaTree.root
    if [ ! -f $anaTree ]; then
      echo "bName = $bName anaTree = $anaTree"
    fi
    done
done                                    
