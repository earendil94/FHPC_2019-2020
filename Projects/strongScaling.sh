#!/bin/bash

ASSIGNMENT_PATH="/home/fbrand/FHPC_2019-2020/D03/code"
CODE_NAME="mpi_pi.x"
N=10000000
RESULT_FILE="time_data.txt"

cd $ASSIGNMENT_PATH

touch $RESULT_FILE

for proc in 1 2 4; do
	time mpirun -np $proc ./$CODE_NAME $N > ${CODE_NAME}_$proc
	MASTER_TIME=$(grep "walltime" ${CODE_NAME}_$proc | grep "master" | cut -d":" -f2)
	echo "${proc};$MASTER_TIME" >> $RESULT_FILE 
	#rm ${CODE_NAME}_$proc
done

more $RESULT_FILE 
