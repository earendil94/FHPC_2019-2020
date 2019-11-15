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
  long long unsigned int *partial_array;
  long long unsigned int partial_sum = 0;
  long long unsigned int total_sum = 0;
  //long long unsigned int *partial_sum_ar;



  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  //printf("Hi, I'm rank %d\n", rank);

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
    //T_read
    n = atoi(argv[1]);

    end_time = MPI_Wtime();
    t_read = end_time - start_time;

  }
  start_time = MPI_Wtime();
  MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
  end_time = MPI_Wtime();
  /*
  else {
    MPI_Recv(&n, 1, MPI_UNSIGNED_LONG, root, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } */

  sub_array_dimension = n/p;
  sub_array_remainder = n%p;
  //printf("sub_array_dimension: %d\n",sub_array_dimension );

  if(rank < sub_array_remainder){
    sub_array_dimension++;
    partial_array = calloc( sub_array_dimension, sizeof(long long unsigned int));
  }else{
    partial_array = calloc( sub_array_dimension, sizeof(long long unsigned int));
  }

  //Populate the array
  for(i = 0; i < sub_array_dimension ; i++){
  	partial_array[i] = i*p + rank+1;
  }

  //Calculation in parallel time
  start_time = MPI_Wtime();

  for(i = 0; i < sub_array_dimension; i++){
    partial_sum += partial_array[i];
    //printf("partial_sum: %d\n", partial_sum );
  }

  end_time = MPI_Wtime();

  t_comp = end_time - start_time;
  //printf("Rank: %d, Our partial sum is: %llu\n", rank, partial_sum);
  //partial_sum_ar = malloc(p*sizeof(long long unsigned int));

  //MPI_Send(&partial_sum, 1, MPI_UNSIGNED_LONG, root , rank, MPI_COMM_WORLD);
  start_time = MPI_Wtime();
  MPI_Reduce(&partial_sum, &total_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, root, MPI_COMM_WORLD );
  end_time = MPI_Wtime();
  t_comm = t_comm + end_time - start_time;

  walltime_end = MPI_Wtime();

  if( rank == root){
    printf( "\n # walltime on processor %i : %10.8f \n",rank, walltime_end - walltime_start ) ;
    printf("Total sum: %llu\n", total_sum);
    printf("T_read: %.9f\n", t_read);
    printf("T_comp: %.9f\n", t_comp);
    printf("T_comm: %.9f\n", t_comm);
  }
  else
    printf( "\n # walltime on processor %i : %10.8f \n",rank, walltime_end - walltime_start ) ;


  MPI_Finalize();

  return 0;

}
