#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./findResubmitRun.sh <outName>"
	exit
fi

if [ -f runNumber_tobeResubmit_$1 ]; then
   rm runNumber_tobeResubmit_$1
fi
touch runNumber_tobeResubmit_$1
for run in `cat submit/runNumber_$1_0`
do

   if [ ! -d submit/$1/${run} ]; then
      echo $run >> runNumber_tobeResubmit_$1
   else
      if [ ! -d out/out_$1/${run} ]; then 
         echo $run >> runNumber_tobeResubmit_$1
      else
         sesName=`basename submit/$1/${run}/*.session.xml | cut -d . -f1`
         nMuDst=`cat submit/$1/${run}/sched${sesName}.list | wc -l`
         nAnaTree=`find out/out_$1/${run}/*.anaTree.root | wc -l`
         nMuDstFromCat=`get_file_list.pl -keys 'filename' -cond 'production=P15ic,collision=auau200,'runnumber=${run},'trgsetupname=AuAu_200_production_2014||AuAu_200_production_mid_2014||AuAu_200_production_low_2014,filename~st_physics,filetype=daq_reco_Mudst,storage=local' -limit 0 | wc -l`
         echo "run = ${run} nMuDst = $nMuDst nAnaTree = $nAnaTree nMuDstFromCat = $nMuDstFromCat"
         if [ $nAnaTree -lt $nMuDst ] && [ $nAnaTree -lt $nMuDstFromCat ]; then
            echo $run >> runNumber_tobeResubmit_$1
         fi
      fi
   fi
done                                    
