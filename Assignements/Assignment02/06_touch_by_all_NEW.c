//This code has been specifically designed to be run on Ulysses on 20 threads on 20 different cores.
//It has been designed having in mind the peculiar architecture of ulysses' nodes (2 sockets, 2 numa zones
//10 cores each socket)


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
#include <unistd.h>
#include <sys/syscall.h>
#include <sched.h>
#include <omp.h>


#define N_default 1000

#if defined(_OPENMP)
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +	\
		     (double)myts.tv_nsec * 1e-9)

#else

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)
#endif

#ifdef OUTPUT
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

#define CPU_ID_ENTRY_IN_PROCSTAT 39
#define HOSTNAME_MAX_LENGTH      200

int read_proc__self_stat ( int, int * );
int get_cpu_id           ( void       );



int main( int argc, char **argv )
{

  int     N        = N_default;
  int     nthreads = 1;
  
  struct  timespec ts;
  double *array;
  double *array_1;
  double *array_2;

  /*  -----------------------------------------------------------------------------
   *   initialize 
   *  -----------------------------------------------------------------------------
   */

  // check whether some arg has been passed on
  if ( argc > 1 )
    N = atoi( *(argv+1) );


  
  // just give notice of what will happen and get the number of threads used
#if defined(_OPENMP)  
#pragma omp parallel
  {
      #pragma omp master
      {
        nthreads = omp_get_num_threads();
        PRINTF("omp summation with %d threads\n", nthreads );
      }

      int me = omp_get_thread_num();

      if(me == 0)
        if ( (array_1 = (double*) malloc(N/2 * sizeof(double)) )== NULL )
        {
          printf("Sorry, not enough memory to allocate array1");
          exit(-1);
        }
      
      if(me == 10)
        if ( (array_2 = (double*) malloc(N/2 * sizeof(double)) )== NULL )
        {
          printf("Sorry, not enough memory to allocate array1");
          exit(-1);
        }

      #pragma omp critical
      PRINTF("thread %2d is running on core %2d\n", me, get_cpu_id() );    
  }
#else  
  if ( (array = (double*)malloc( N * sizeof(double) )) == NULL ) {
    printf("I'm sorry, on some thread there is not"
	   "enough memory to host %lu bytes\n",
	   N * sizeof(double) ); return 1; }
#endif

  // initialize the array;
  // each thread is "touching"
  // its own memory as long as
  // the parallel for has the
  // scheduling as the final one
#pragma omp parallel
{
  int t = omp_get_thread_num();

  if( t < 10){
    for ( int ii = t; ii < N/2; ii+=10)
      array_1[ii] = (double) ii;
  }
  if (t >= 10){
      for (int ii = t-10; ii < N/2; ii+=10)
        array_2[ii] = (double) (ii + N/2);
  }

}

  /*  -----------------------------------------------------------------------------
   *   calculate
   *  -----------------------------------------------------------------------------
   */


  double S  = 0;                                       // this will store the summation

  double th_avg_time = 0;                                   // this will be the average thread runtime
  double th_min_time = 1e11;                                // this will be the min thread runtime.
							    // contrasting the average and the min
							    // time taken by the threads, you may
							    // have an idea of the unbalance.
  
  double tstart  = CPU_TIME;

#if !defined(_OPENMP)
  
  for ( int ii = 0; ii < N; ii++ )                          // well, you may notice this implementation
    S += array[ii];                                         // is particularly inefficient anyway

#else

#pragma omp parallel reduction(+:th_avg_time)				\
  reduction(min:th_min_time)                                // in this region there are 2 different
  {                                                         // reductions: the one of runtime, which
    struct  timespec myts;                                  // happens in the whole parallel region;
    double mystart = CPU_TIME_th; 
    int t = omp_get_thread_num();
    double S1 = 0;
    double S2 = 0;

    if(t < 10)       
      for ( int ii = t; ii < N/2; ii+=10 )
        S1 += array_1[ii];
    if (t >= 10)
      for( int ii = t-10; ii < N/2; ii+=10)
        S2 += array_2[ii];

    #pragma omp critical
      S += S1 + S2;


    double mytime = CPU_TIME_th - mystart; 
    th_avg_time += mytime;
    th_min_time  = (mytime < th_min_time)? mytime : th_min_time;
  }

#endif
  
  double tend = CPU_TIME;


  /*  -----------------------------------------------------------------------------
   *   finalize
   *  -----------------------------------------------------------------------------
   */
  //Changing printf format in order to equalize it to the 01 exercise
  printf("Sum is %g\nwall-clock-time:%g\n\n"
	 "<%g> sec of avg thread-time\n"
	 "<%g> sec of min thread-time\n",
	 S, tend - tstart, th_avg_time/nthreads, th_min_time );

#if !defined(_OPENMP)
  
  free( array );
#else
  free( array_1);
  free( array_2);
#endif


return 0;
}






int get_cpu_id( void )
{
#if defined(_GNU_SOURCE)                              // GNU SOURCE ------------
  
  return  sched_getcpu( );

#else

#ifdef SYS_getcpu                                     //     direct sys call ---
  
  int cpuid;
  if ( syscall( SYS_getcpu, &cpuid, NULL, NULL ) == -1 )
    return -1;
  else
    return cpuid;
  
#else      

  unsigned val;
  if ( read_proc__self_stat( CPU_ID_ENTRY_IN_PROCSTAT, &val ) == -1 )
    return -1;

  return (int)val;

#endif                                                // -----------------------
#endif

}



int read_proc__self_stat( int field, int *ret_val )
/*
  Other interesting fields:

  pid      : 0
  father   : 1
  utime    : 13
  cutime   : 14
  nthreads : 18
  rss      : 22
  cpuid    : 39

  read man /proc page for fully detailed infos
 */
{
  // not used, just mnemonic
  // char *table[ 52 ] = { [0]="pid", [1]="father", [13]="utime", [14]="cutime", [18]="nthreads", [22]="rss", [38]="cpuid"};

  *ret_val = 0;

  FILE *file = fopen( "/proc/self/stat", "r" );
  if (file == NULL )
    return -1;

  char   *line = NULL;
  int     ret;
  size_t  len;
  ret = getline( &line, &len, file );
  fclose(file);

  if( ret == -1 )
    return -1;

  char *savetoken = line;
  char *token = strtok_r( line, " ", &savetoken);
  --field;
  do { token = strtok_r( NULL, " ", &savetoken); field--; } while( field );

  *ret_val = atoi(token);

  free(line);

  return 0;
}

