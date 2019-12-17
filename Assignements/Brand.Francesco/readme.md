#Assignment 02 
##Francesco Brand 

The files have been organized in three directories with this distinction:

- csv: here we have the csv files containing the data used for the plots shown in the report. The files have been produced on Ulysses cluster
- source: in this directory we have the C source files for the programs used in this assignment. Some we already provided and have been only slightly modified, others have been written from the scratch.
- scripts: here are present the scripts used to run the programs on ulysses.

The user has two options in order to reproduce the results shown here:

1) On ulysses cluster run qsub wrapper after having modified the path in which you cd in the wrapper script itself. The user also wants to adapt the EXE variable to the program name which he/she desires to run. The program name is the source name without the extension .c

2) Manually compile the files on an Ulysses node after having loaded gcc module with the command 
        gcc -std=c99 -fopenmp -o ${EXE}-PARALLEL $SOURCE -lrt
   for the parallel version of the code and
       gcc -std=c99 -o ${EXE}-SERIAL $SOURCE -lrt
   for the non-parallel version.

   Set the desired number of threads to be spawned with export OMP_NUM_THREADS=p (where p can be 1,2,...,20 on Ulysses)
   then run the program with /usr/bin/time -f%e:elapsed ./program N in order to obtain the elapsed time of execution and collect the output in order to produce the csv.
   
If anything is not clear, or for any further question, please contact fbrand@sissa.ità

