#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>


#if defined(_OPENMP)
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec + \
		     (double)myts.tv_nsec * 1e-9)
#else

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
		   (double)ts.tv_nsec * 1e-9)

#endif

#define DEFAULT 8

int main(int argc, char **argv){

    int n = DEFAULT;
    double *arr;

    //If we actually pass argument take the first argument as N
    if(argc > 1)
        n = atoi(*(argv+1));

    

    if( (arr = (double *) malloc(n*sizeof(double))) == NULL ){
        printf("Malloc failed, exiting program\n");
        exit(-1);
    }

    for(int k = 0; k < n; k++){
        arr[k] = (double)k;
    }

    // printf("Arr:\t");
    // for(size_t k = 0; k < n; k++){
    //      printf("%f\t",arr[k]);
    // }
    // printf("\n\n");
   

    //True part of the algorithm
    #if !defined(_OPENMP)
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
            #pragma omp master
                p = omp_get_num_threads();
           
            int t = omp_get_thread_num();
            int partial_dim = n/p; 

            #pragma omp master
                sum = malloc(sizeof(double)*p);


            sum[t] = arr[partial_dim*t];
            for(size_t i = 1 + partial_dim*t; i < partial_dim*(t+1); i++)
                sum[t] += arr[i];
        }

        for(size_t k = 1; k < p; k++)
            sum[k] += sum[k-1];

        #pragma omp parallel
        {   
            int p = omp_get_num_threads();
            int t = omp_get_thread_num();
            int partial_dim = n/p; 
            // if( t > 0)
            //     arr[partial_dim*(t)] = arr[partial_dim * (t)] + sum[t-1];

            if( t > 0)
                arr[partial_dim * t] += sum[t-1];

            for( size_t i = 1 + partial_dim *t; i < partial_dim * (t+1); i++)
                arr[i] = arr[i] + arr[i-1];
        } 

    #endif


    //#TODO:maybe complete the log(p) version of the algortihm, maybe not
    //Up-Sweep
    // int height = ceil(log2(n));

    // for( int d = 0; d < height; d++)
    // {
    //     int step = (int) pow(2,d+1);
    //     #pragma omp parallel for
    //         for(int i = 0; i < n ; i += step)
    //             arr[i + step - 1] = arr[i + step - 1] + arr[i + step/2 - 1];
    // }

    // //Down-sweep
    // arr[n-1] = 0;

    // for( int d = height - 1; d >= 0; d--)
    // {
    //     int step = (int) pow(2,d+1);
    //     #pragma omp parallel for
    //         for(int i = 0; i < n; i += step)
    //         {
    //             double t = arr[i + step/2 - 1];
    //             arr[i + step/2 - 1] = arr[i + step - 1];
    //             arr[i + step - 1] = arr[i + step - 1] + t;
    //         }
    // }
    
    // printf("\nprefixSumArr:\t");
    // for(size_t k = 0; k < n; ++k){
    //     printf("%f\t",arr[k]);
    // }

}