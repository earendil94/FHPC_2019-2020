#!/bin/bash
#Please note that it is much better if we just collect all the files we need in our directory
#Pass n as first argument, result file as second argument, program as third argument
run_serial(){

	#echo $1
	#echo $2
	#echo $3
	#echo "$3 $1"
	{ time ./$3 $1 ; } 1>$2 2> $2
	#time ./$3 $1
}

#Pass p as first argument, n as second argument, result file as third argument, program as fourth argument
run_parallel(){
	{ time mpirun -np $1 $4 $2 ; } 1> $3 2> $3
}

#Pass filename_p as first argument, csv name as second argument
#Our csv will have these columns: Rank,P,Time,N
result_files_to_csv(){

	export P=$(echo $1 | cut -d"_" -f3)
	export N=$(echo $1 | cut -d"_" -f4)
	#In case we are extracting walltime from serial app

	if [[ $CURR_PROGRAM == *"PartialSum"* ]]
	then
		serial_sum_elab $1 $2
	fi

	if (( $P == 1 ))
	then
		WALLTIME=$(grep "walltime" $1 | cut -d":" -f2)
		if [[ $WALLTIME == "" ]]
		then
			USERTIME=$(grep "user" $1 | cut -d"m" -f2 | cut -d"s" -f1)
			echo "${P};${P};${USERTIME};${N}" >> $2
		else
			echo "${P};${P};${WALLTIME};${N}" >> $2
		fi
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

#First argument: file_name; Second argument: csv_name
serial_sum_elab(){
	T_READ=$(grep "T_read" $1 | cut -d":" -f2)
	T_COMM=$(grep "T_comm" $1 | cut -d":" -f2)
	T_COMP=$(grep "T_comp" $1 | cut -d":" -f2)

	echo "T_READ;${P};${T_READ};${N}" >> $2.times
	echo "T_COMM;${P};${T_COMM};${N}" >> $2.times
	echo "T_COMP;${P};${T_COMP};${N}" >> $2.times
}

CURR_PROGRAM=$1
SCALING_TEST=$2

if [[ $2 != "WEAK" ]]
then
	SCALING_TEST=STRONG
fi


if [[ $CURR_PROGRAM == *"mpi"* ]]
then
	module load openmpi
fi

P_VALUES=$(grep P_VALUES_${SCALING_TEST} values.cfg | cut -d"=" -f2)
N_VALUES=$(grep N_VALUES_${SCALING_TEST} values.cfg | cut -d"=" -f2)
RESULT_FILE=$(echo $CURR_PROGRAM | cut -d"." -f1)_${SCALING_TEST}.res
CSV_FILE=$(echo $CURR_PROGRAM | cut -d"." -f1)_${SCALING_TEST}.csv

#echo $CURR_PROGRAM

if [[ $CURR_PROGRAM == *"mpi"* ]]
then
	if [[ $SCALING_TEST == "STRONG" ]]
	then
		for n_val in $N_VALUES
		do
			for p_val in $P_VALUES
			do
				run_parallel $p_val $n_val ${RESULT_FILE}_${p_val}_${n_val} $CURR_PROGRAM
				result_files_to_csv ${RESULT_FILE}_${p_val}_${n_val} $CSV_FILE
				rm ${RESULT_FILE}_${p_val}_${n_val}
			done
		done
	elif [[ $SCALING_TEST == "WEAK" ]]
	then
		for n_ in $N_VALUES
		do
			for p_val in $P_VALUES
			do
				n_val=$(( p_val * $n_ ))
				run_parallel $p_val $n_val ${RESULT_FILE}_${p_val}_${n_val} $CURR_PROGRAM
				result_files_to_csv ${RESULT_FILE}_${p_val}_${n_val} $CSV_FILE
				rm ${RESULT_FILE}_${p_val}_${n_val}
			done
		done
	fi
else
	if [[ $SCALING_TEST == "STRONG" ]]
	then
		for n_val in $N_VALUES
		do
			run_serial $n_val ${RESULT_FILE}_1_${n_val} $CURR_PROGRAM
			result_files_to_csv ${RESULT_FILE}_1_${n_val} $CSV_FILE
			rm ${RESULT_FILE}_1_${n_val}
		done
	elif [[ $SCALING_TEST == "WEAK" ]]
	then
		for n_ in $N_VALUES
		do
			for p_val in $P_VALUES
			do
				n_val=$(( p_val * $n_ ))
				run_serial $n_val ${RESULT_FILE}_1_${n_val} $CURR_PROGRAM
				result_files_to_csv ${RESULT_FILE}_1_${n_val} $CSV_FILE
				rm ${RESULT_FILE}_1_${n_val}
			done
		done
	fi
fi


#Use millisecond instead of second
awk 'BEGIN { FS=";" } {printf"%d;%d;%.7f;%d\n", $1, $2, $3*1000, $4}' $CSV_FILE > ${CSV_FILE}_out
rm $CSV_FILE
mv ${CSV_FILE}_out $CSV_FILE
