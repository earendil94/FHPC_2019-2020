# Usage Instruction

## Compilation

### Mpi code

All the source codes named mpi\* must be compiled in this way:

	- (If on Ulysses): module load openmpi
	- mpicc -o code\_name.x code\_name.c -std=c99

### Non mpi code

All the other source codes not named mpi\* must be compiled with gcc -o code\_name.x code\_name.c -std=c99

## Execution

In the directory provided you will find a general bash script to run the tests needed for the report called run\_exe.sh. Its usage is ./run\_exe name\_of\_the\_program.exe "STRONG/WEAK" ["ELAPSED"]. It can be used for both MPI and non-mpi programs, for weak scalability tests and strong scalability tests.  

The script returns a csv formatted in 4 columns, respectively: RankOfProcessor,NumberOfProcessors,ExecutionTime (ms), N.

As a special case, the mpiPartialSum code also returns an extra csv with $T_{read}, T_{comp}, T_{comm}$ for every processor for every N considered.

## Analysis
The code returns consistent result only for $N \approx 10^9$. If one tries to input values of higher magnitude, the total sum will be affected by overfloating of our numbers.
The plots provided in the report have been produced with an additional jupyter notebook which I will provide in the directory.   
N.B. The csv returned by the serial algorithms sometimes are not sensitive enough to get a meaningful walltime on Ulysses for low values of N, therefore we used for our plot in section 4 with $N=10^6$ the result yielded on a local execution of the program.