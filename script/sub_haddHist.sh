#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./sub_haddHist.sh <outName>"
	exit
fi

if [ ! -d log/log_$1_hist ]; then
   mkdir log/log_$1_hist
fi

if [ -f log/log_$1_hist/runTimes ]; then
   ijob=`cat log/log_$1_hist/runTimes`
else
   ijob=0
fi
while [ 1 ]
do
   if [ -f  log/log_$1_hist/log_${ijob} ]; then
      rm log/log_$1_hist/log_${ijob}
   fi
   ./haddHist.sh $1 >& log/log_$1_hist/log_${ijob} 
   let "ijob+=1"
   if [ -f log/log_$1_hist/runTimes ]; then 
      rm log/log_$1_hist/runTimes
   fi
   echo $ijob >& log/log_$1_hist/runTimes
   echo "Done! go to bed!"
   sleep 86400 #24 hours
done
