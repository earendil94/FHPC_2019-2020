#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char **argv){

  MPI_Init(&argc,&argv);

  const int root = 0;

  double start_time, end_time;
  int p, rank;
  long long unsigned int i, n;
  long long unsigned int *partial_array;
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

  if( rank == root){

    //Time def
    start_time = MPI_Wtime();
    //T_read
    n = atoi(argv[1]);

    end_time = MPI_Wtime();
    double t_read = end_time - start_time;
    printf("T_read: %.9f\n", t_read);

    start_time = MPI_Wtime();

    for( i = root+1; i < p; i++){
        MPI_Send(&n, 1, MPI_INT, i, i, MPI_COMM_WORLD);
    }

    end_time = MPI_Wtime();
    double t_comm = end_time - start_time;
    printf("T_comm_broadcast: %.9f\n", t_comm);
  } else {
    MPI_Recv(&n, 1, MPI_INT, root, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  sub_array_dimension = n/p;
  sub_array_remainder = n%p;

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
  printf("T_comp is: %.9f\n", end_time - start_time);

  partial_sum_ar = malloc(p*sizeof(long long unsigned int));

  MPI_Send(&partial_sum, 1, MPI_INT, root , rank, MPI_COMM_WORLD);


  if( rank == root){

    start_time = MPI_Wtime();

    for(i = 1; i < p; i++){
      MPI_Recv( &partial_sum_ar[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      partial_sum += partial_sum_ar[i];
    };

    end_time = MPI_Wtime();
    printf("Receving time+calculation time: %.9f\n", end_time - start_time );
    printf("Our partial sum is: %d\n", partial_sum);
  }


  MPI_Finalize();

  return 0;

}
