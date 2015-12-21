#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./tarSubmitLogs.sh <outName>"
	exit
fi

if [ ! -d submit_backup/$1 ]; then
   mkdir submit_backup/$1
fi

for run in `cat submit/$1/runNumber_$1_0`
do
   if [ ! -d submit/$1/${run} ]; then
      tar -zcvf submit_${run}.tgz submit/$1/${run}
      mv submit_${run}.tgz submit_backup/$1/
      rm -rf submit/$1/${run}
   fi
done
