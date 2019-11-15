#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char **argv){

  if( argc != 2){
    fprintf(stderr, "Program should pass the N of the sum to compute \n");
    exit(-1);
  }

  clock_t start_time, end_time;

  start_time = clock();

  int n = atoi(argv[1]);

  end_time = clock();
  double read_time = (double) (end_time - start_time)/CLOCKS_PER_SEC;

  printf("The code has read n in %f s\n", read_time );

  start_time = clock();
  long long unsigned int sum = 0;

  for(int i = 0; i <= n; i++)
    sum += i;

  printf("The sum from 0 to %d is %lld\n",n, sum );

  end_time = clock();

  double comp_time = (double) (end_time - start_time)/CLOCKS_PER_SEC;

  printf("The code has run in %f s\n", read_time + comp_time );

  return 0;

}
