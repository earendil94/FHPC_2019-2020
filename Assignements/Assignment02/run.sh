#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Usage: ./run.sh source n_dimension n_threads loop_dimension"
    echo "n_dimension: power of 10"
    echo "n_threads: 1-20"
    echo "loop_dimension suggested: 10"
fi

SOURCE=$1
EXE=$( echo $SOURCE | cut -d. -f1)
N_DIMENSION=$2
N_THREADS=$3
LOOP=$4

i=0

if [ $N_THREADS -eq 1 ]
then
    FLAGPS=S
    EXE=${EXE}_SERIAL
else
    FLAGPS=P
    EXE=${EXE}_PARALLEL
    export OMP_NUM_THREADS=${N_THREADS}
fi

./compile.sh $SOURCE $EXE $FLAGPS

while [ $i -lt $LOOP ]
do
    ./${EXE} $N_DIMENSION | grep "wall-clock-time" | cut -d: -f2 >> ${EXE}_${N_DIMENSION}_${N_THREADS}.csv
    i=$(( $i + 1 ))
done


