#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>


#define DEFAULT 10

int main(int argc, char **argv){

    int n = DEFAULT;
    int nthreads = 1;
    double *arr;

    //If we actually pass argument take the first argument as N
    if(argc > 1)
        n = atoi(*(argv+1));


    if( (arr = (double *) malloc(n*sizeof(double))) == NULL ){
        printf("Malloc failed, exiting program\n");
        exit(-1);
    }





    //True part of the algorithm
    #if !defined(_OPENMP)

        for(int k = 0; k < n; k++)
              arr[k] = (double)k;


        for(size_t i = 1; i < n; ++i)
            arr[i] = arr[i] + arr[i-1];
    #else

        //Let's start assuming n/p is a perfect integer division
        double *sum;
        int p;
        // if( (sum = malloc(sizeof(double)*p)) == NULL )
        // {
        //     printf("Malloc failed, exiting program\n");
        //     exit(-1);
        // }

        #pragma omp parallel
        {

            //Touch by all
            #pragma omp for 
              for(int register k = 0; k < n; k++)
                arr[k] = (double)k;

            // #pragma omp single 
            // {
            //     printf("Arr:\t");
            //     for(size_t k = 0; k < n; k++){
            //       printf("%f\t",arr[k]);
            //     }
            //     printf("\n\n");
            // }

            #pragma omp single
            {
              p = omp_get_num_threads();
              sum = malloc(sizeof(double)*(p+1));
            }

            sum[0] = 0;
            int t = omp_get_thread_num();         
            
            double tsum = 0;

            #pragma omp for schedule(static)
            for(int i = 0; i < n; i++)
            {
                tsum += arr[i];
                arr[i] = tsum;
            }

            sum[t+1] = tsum;

            #pragma omp barrier

            int offset;

            for(int i = 0; i < t+1; i++)
              offset += sum[i];

            #pragma omp for schedule(static)
            for( int register i = 0; i < n; i++)
            {
                arr[i] += offset;
            }
        } 

    #endif
    
    // printf("Arr:\t");
    // for(size_t k = 0; k < n; k++){
    //      printf("%f\t",arr[k]);
    // }
    // printf("\n\n");

}

