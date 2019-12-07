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

#define DEFAULT 1000

int main(int argc, char **argv){

    int n = DEFAULT;
    double *arr;

    //If we actually pass argument take the first argument as N
    if(argc > 1)
        n = atoi(*(argv+1));

    

    if( ! (double *) malloc(n*sizeof(double)) ){
        printf("Malloc failed, exiting program\n");
        exit(-1);
    }

    for(int k = 0; k < n; k++){
        arr[k] = (double)k;
    }

    // printf("Arr:\t");
    // for(size_t k = 0; k < n; k++){
    //     printf("%f\t",arr[k]);
    // }

    //Naive mode for serial: on
    //Let's try something fancy for parallel
    #ifndef _OPENMP
        for(size_t i = 1; i < n; ++i){
            arr[i] = arr[i] + arr[i-1];
        }
    #else
        #pragma omp parallel
            //for ... #TODO

    #endif


    
    // printf("\nprefixSumArr:\t");
    // for(size_t k = 0; k < n; ++k){
    //     printf("%f\t",arr[k]);
    // }

}