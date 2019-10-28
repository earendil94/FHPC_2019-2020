#!/bin/bash

PATH_TO_CODE=$(grep path_to_file /home/fbrand/FHPC_2019-2020/PersonalCode/Scripts/configParallel.cfg | cut -d"=" -f2)
CODE_NAME=$(grep program_name /home/fbrand/FHPC_2019-2020/PersonalCode/Scripts/configParallel.cfg | cut -d"=" -f2)

module load openmpi

cd $PATH_TO_CODE
mpirun -np 4 $CODE_NAME

