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
#Our csv will have these columns: Rank,P,Time,N
result_files_to_csv(){

	P=$(echo $1 | cut -d"_" -f2)
	N=$(echo $1 | cut -d"_" -f3)
	#In case we are extracting walltime from serial app
	if (( $P == 1 ))
	then
		WALLTIME=$(grep "walltime" $1 | cut -d":" -f2)
		echo "${P};${P};${WALLTIME};${N}" >> $2
	else
		for ((rank=0;rank<$P;rank++))
		do
			if (( $rank == 0 ))
			then
				WALLTIME=$(grep "walltime" $1 | grep "master" | cut -d":" -f2)
			else
				WALLTIME=$(grep "walltime" $1 | grep -w "processor ${rank}" | cut -d":" -f2)
			fi
			echo "${rank};${P};${WALLTIME};${N}" >> $2
		done
	fi
}

CURR_PROGRAM=$1
SCALING_TEST=$2

if [[ $CURR_PROGRAM == *"mpi"* ]]
then
	module load openmpi
fi

P_VALUES=$(grep P_VALUES_{$SCALING_TEST} values.cfg | cut -d"=" -f2)
N_VALUES=$(grep N_VALUES_{$SCALING_TEST} values.cfg | cut -d"=" -f2)
RESULT_FILE=$(echo $CURR_PROGRAM | cut -d"." -f1)_{$SCALING_TEST}.res
CSV_FILE=$(echo $CURR_PROGRAM | cut -d"." -f1)_{$SCALING_TEST}.csv

#echo $CURR_PROGRAM
if [[ $SCALING_TEST == "STRONG" ]]
then
	for n_val in $N_VALUES
	do
		for p_val in $P_VALUES
		do
			if [[ $CURR_PROGRAM == *"mpi"* ]]
			then
				run_parallel $p_val $n_val ${RESULT_FILE}_${p_val}_${n_val} $CURR_PROGRAM
				result_files_to_csv ${RESULT_FILE}_${p_val}_${n_val} $CSV_FILE
				rm ${RESULT_FILE}_${p_val}_${n_val}
			else
				run_serial $n_val ${RESULT_FILE}_${p_val}_${n_val} $CURR_PROGRAM
				result_files_to_csv ${RESULT_FILE}_${p_val}_${n_val} $CSV_FILE
				rm ${RESULT_FILE}_${p_val}_${n_val}
			fi
		done
	done
else
	for p_val in $P_VALUES

#Use millisecond instead of second
awk 'BEGIN { FS=";" } {printf"%d;%d;%.7f;%d\n", $1, $2, $3*1000, $4}' $CSV_FILE > ${CSV_FILE}_out
#rm $CSV_FILE
#mv ${CSV_FILE}_out $CSV_FILE
