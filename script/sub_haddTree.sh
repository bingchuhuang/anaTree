#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./sub_haddTree.sh <outName>"
	exit
fi

if [ ! -d log/log_$1_anaTree ]; then
   mkdir log/log_$1_anaTree
fi

if [ -f log/log_$1_anaTree/runTimes ]; then
   ijob=`cat log/log_$1_anaTree/runTimes`
else
   ijob=0
fi
while [ 1 ]
do
   if [ -f  log/log_$1_anaTree/log_${ijob} ]; then
      rm log/log_$1_anaTree/log_${ijob}
   fi
   ./haddTree.sh $1 >& log/log_$1_anaTree/log_${ijob}
   let "ijob+=1"
   if [ -f log/log_$1_anaTree/runTimes ]; then 
      rm log/log_$1_anaTree/runTimes
   fi
   echo $ijob >& log/log_$1_anaTree/runTimes
   echo "Done! go to bed!"
   sleep 28800 #8 hours
done
