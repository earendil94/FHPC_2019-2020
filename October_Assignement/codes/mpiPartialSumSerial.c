#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char **argv){

  MPI_Init(&argc,&argv);

  const int root = 0;

  double start_time, end_time, t_read, t_comm, t_comp, walltime_start, walltime_end;
  int p, rank;
  long long unsigned int i, n;
  //long long unsigned int *partial_array;
  long long unsigned int partial_sum = 0;
  long long unsigned int *partial_sum_ar;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  if( argc != 2){
    fprintf(stderr, "Program should pass the N of the sum to compute \n");
    exit(-1);
  }

  long long unsigned int sum = 0;

  int sub_array_dimension;
  int sub_array_remainder;

  walltime_start = MPI_Wtime();
  if( rank == root){

    //Time def
    start_time = MPI_Wtime();

    n = atoi(argv[1]);

    end_time = MPI_Wtime();
    t_read = end_time - start_time;


    start_time = MPI_Wtime();

    for( i = root+1; i < p; i++){
        MPI_Send(&n, 1, MPI_LONG_LONG, i, i, MPI_COMM_WORLD);
    }

    end_time = MPI_Wtime();
    t_comm = end_time - start_time;

  } else {
    MPI_Recv(&n, 1, MPI_LONG_LONG, root, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  sub_array_dimension = n/p;
  sub_array_remainder = n%p;

  if(rank < sub_array_remainder)
    sub_array_dimension++;


  for(i = 0; i < sub_array_dimension; i++)
    partial_sum += i*p + (rank+1);

  end_time = MPI_Wtime();

  t_comp = end_time - start_time;

  partial_sum_ar = malloc(p*sizeof(long long unsigned int));

  MPI_Ssend(&partial_sum, 1, MPI_UNSIGNED_LONG, root , rank, MPI_COMM_WORLD);


  if( rank == root){

    start_time = MPI_Wtime();

    for(i = 1; i < p; i++)
      MPI_Recv( &partial_sum_ar[i], 1, MPI_UNSIGNED_LONG, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    end_time = MPI_Wtime();

    t_comm = t_comm + end_time - start_time;

    start_time = MPI_Wtime();

    for(i = 1; i<p; i++)
      partial_sum += partial_sum_ar[i];

    end_time = MPI_Wtime();

    t_comp = t_comp + end_time - start_time;
    walltime_end = MPI_Wtime();
    printf("T_read: %.9f\n", t_read);
    printf("T_comp: %.9f\n", t_comp);
    printf("T_comm: %.9f\n", t_comm);
    printf ( "\n # walltime on master processor : %10.8f \n", walltime_end - walltime_start ) ;
    printf("Our partial sum is: %llu\n", partial_sum);
  } else{
    walltime_end = MPI_Wtime();
    printf ( "\n # walltime on processor %i : %10.8f \n",rank, walltime_end - walltime_start ) ;
  }


  MPI_Finalize();

  return 0;

}
