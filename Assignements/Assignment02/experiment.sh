#!/bin/bash

#I'd like to have a list of numbers to use as a power
PROGRAM=$1
SOURCE=${PROGRAM}.c
EXP=(3 4 5 6 7 8 9 10)
N_THREADS=(1 2 3 4 5 6 7 8)
LOOP=10

n=0
t=0



for n in EXP
do
    N=(( 10 ** n ))
    ./run $SOURCE $N $N_THREADS $LOOP 
    awk 'BEGIN { FS =";"} '
done