#!/bin/bash
date

if [ $# -lt 1 ]; then
    echo "not enough arguements"
    echo "./clear.sh <outName>"
	exit
fi

rm -rf submit/*
rm -rf log/log_$1/*                
rm -rf out/out_$1/*
