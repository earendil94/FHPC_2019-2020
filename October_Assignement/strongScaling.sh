#!/bin/bash

ASSIGNMENT_PATH="/home/francesco/Desktop/FHPC_2019-2020/D03/code"
CODE_NAME="pi.x"
N=10000000

cd $ASSIGNMENT_PATH

for proc in 1 2 4 8 16 32 64; do
	./$CODE_NAME $N > ${CODE_NAME}_$proc
	grep "walltime" ${CODE_NAME}_$proc 
done


