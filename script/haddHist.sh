#!/bin/bash
date

source $ROOTSYS/bin/thisroot.sh

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./haddHist.sh <outName>"
	exit
fi

if [ ! -d out/out_$1_hist ]; then
   mkdir out/out_$1_hist
fi

for run in `cat runNumberList_run14AuAu200mb`
do
   if [ ! -d out/out_$1/${run} ]; then
      continue
   fi
   ifile=`find out/out_$1/${run} -name "st*${run}*hist.root" | wc -l`
   if [ $ifile -gt 0 ]; then
      echo "run ${run}"
      echo "found $ifile hist files"
      if [ ! -f st_physics_${run}.hist.root ]; then
         if [ ! -f out/out_$1_hist/${run}_hist_unfinish ]; then
            touch out/out_$1_hist/${run}_hist_unfinish
            hadd -f out/out_$1_hist/st_physics_${run}.hist.root `find out/out_$1/${run} -name "st*${run}*hist.root"`
            rm -f out/out_$1_hist/${run}_hist_unfinish
            rm -f out/out_$1/st*${run}*hist.root
         fi
      else
         if [ ! -f out/out_$1_hist/${run}_hist_unfinish ]; then
            touch out/out_$1_hist/${run}_hist_unfinish
            hadd -f out/out_$1_hist/st_physics_${run}.hist_new.root  out/out_$1_hist/st_physics_${run}.hist.root `find out/out_$1/${run}  -name "st*${run}*hist.root"`
            mv out/out_$1_hist/st_physics_${run}.hist_new.root out/out_$1_hist/st_physics_${run}.hist.root
            rm -f out/out_$1_hist/${run}_hist_unfinish
            rm -f out/out_$1/${run}/st*${run}*hist.root
         fi
      fi
   fi

   ifile=`find out/out_$1/${run} -name "st*${run}*qa.root" | wc -l`
   if [ $ifile -gt 0 ]; then
      echo "found $ifile qa files"
      if [ ! -f st_physics_${run}.qa.root ]; then
         if [ ! -f out/out_$1_hist/${run}_qa_unfinish ]; then
            touch out/out_$1_hist/${run}_qa_unfinish
            hadd -f out/out_$1_hist/st_physics_${run}.qa.root `find out/out_$1/${run} -name "st*${run}*qa.root"`
            rm -f out/out_$1_hist/${run}_qa_unfinish
            rm -f out/out_$1/${run}/st*${run}*qa.root
         fi
      else
         if [ ! -f out/out_$1_hist/${run}_qa_unfinish ]; then
            touch out/out_$1_hist/${run}_qa_unfinish
            hadd -f out/out_$1_hist/st_physics_${run}.qa_new.root  out/out_$1_hist/st_physics_${run}.qa.root `find out/out_$1/${run}/ -name "st*${run}*qa.root"`
            mv out/out_$1_hist/st_physics_${run}.qa_new.root out/out_$1_hist/st_physics_${run}.qa.root
            rm -f out/out_$1_hist/${run}_qa_unfinish
            rm -f out/out_$1/${run}/st*${run}*qa.root
         fi
      fi
   fi
done
