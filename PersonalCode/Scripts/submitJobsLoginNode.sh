#!/bin/bash

NODES_NUM=$(grep nodes configSubmitJobLogin.cfg | cut -d"=" -f2)
PROCESS_NUM=$(grep ppn configSubmitJobLogin.cfg | cut -d"=" -f2)

echo $NODES_NUM
echo $PROCESS_NUM

qsub -l nodes=$NODES_NUM:ppn=$PROCESS_NUM  /home/fbrand/FHPC_2019-2020/PersonalCode/Scripts/submitParallel.sh
