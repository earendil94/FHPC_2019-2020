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
  double comp_time;
  unsigned int i = 0;
  long long unsigned int sum = 0;

  start_time = clock();

  int n = atoi(argv[1]);


  for(i = 0; i <= n; i++)
    sum += i;

  end_time = clock();


  comp_time = ( (double) (end_time - start_time) )/CLOCKS_PER_SEC ;

  printf("The sum from 0 to %d is %lld\n",n, sum );
  printf("walltime : %10.8f \n",comp_time );

  return 0;

}
