#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./compareCode.sh <source dir>"
	exit
fi

srcDir=$1

for dir in `ls -l `
do
   if [ -f $srcDir/$dir ] || [ -d $srcDir/$dir ]; then
      echo "comparing $dir"
      diff $dir $srcDir/$dir
   fi
done
