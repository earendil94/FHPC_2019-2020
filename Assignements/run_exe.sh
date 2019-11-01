#!/bin/bash
#Please note that it is much better if we just collect all the files we need in our directory
#Pass n as first argument, result file as second argument, program as third argument
run_serial(){

	#echo $1
	#echo $2
	#echo $3
	time ./$3 $1 > $2
}

#Pass p as first argument, n as second argument, result file as third argument, program as fourth argument
run_parallel(){
	time mpirun -np $1 $4 $2 > $3
}

#Pass filename_p as first argument, string to grep as second argument, csv name as third
result_files_to_csv(){

	WALLTIME=$(grep $2 $1 | cut -d":" -f2)
	#USRTIME=$(grep $4 $1 | cut)
	P=$(echo $1 | cut -d"_" -f2)
	N=$(echo $1 | cut -d"_" -f3)
	echo "${P};${OUR_TIME};${N}" >> $3

}

#To use on ulysses
#ASSIGNMENT_PATH="/home/fbrand/FHPC"
#cd $ASSIGNMENT_PATH

CURR_PROGRAM=$1
if [[ $CURR_PROGRAM == *"mpi"* ]]
then
	module load openmpi
fi

#echo $CURR_PROGRAM

P_VALUES=$(grep P_VALUES values.cfg | cut -d"=" -f2)
N_VALUES=$(grep N_VALUES values.cfg | cut -d"=" -f2)
RESULT_FILE=$(echo $CURR_PROGRAM | cut -d"." -f1).res
CSV_FILE=$(echo $CURR_PROGRAM | cut -d"." -f1).csv

#echo $RESULT_FILE
#echo $P_VALUES
#echo $N_VALUES

for n_val in $N_VALUES
do
	for p_val in $P_VALUES
	do
		if [[ $CURR_PROGRAM == *"mpi"* ]]
		then
			run_parallel $p_val $n_val $RESULT_FILE $CURR_PROGRAM
		else
			run_serial $n_val ${RESULT_FILE}_${p_val}_${n_val} $CURR_PROGRAM
			result_files_to_csv ${RESULT_FILE}_${p_val}_${n_val} "walltime|usr" $CSV_FILE
			rm ${RESULT_FILE}_${p_val}_${n_val}
		fi
	done
done
