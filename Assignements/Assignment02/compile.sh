#!/bin/bash

SOURCE=$1
EXE=$2
FLAGPS=$3

#Serial
if [ $FLAGPS == S ]
then
    gcc -o $EXE $SOURCE

else
    gcc -fopenmp -o $EXE $SOURCE
fi
