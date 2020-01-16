#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sched.h>
#include <mpi.h>
#include <math.h>

#define G 10e-06
#define epsilon 0.5
#define gamma 0.05

//This represents a particle position, velocity, force (in a 3d vector) and energy
typedef struct { 
    double p[3];
    double v[3];
    double E;
} particle;

particle * readParticles(int rank, int size, int n, MPI_File *ifp, MPI_Datatype *MPI_PARTICLE);
double nBodyEvolution(particle *par, int rank, int size, int displacement, int offset, int n, double m);
void updateParticles(particle *par, int displacement, int offset, double m, double t, double *Fx, double *Fy, double *Fz);
int sgn(double a, double b);

int main(int argc, char **argv){

    //General usage variables
    unsigned short int n_it = 0;

    //To be read from ic file
    int floatPointPrecision;
    int n; //We will use this to sto
    int numOfFiles; //This is the number of files in which we split the I
    double time; //This is the initial condition generation, so we are at time zero


    //MPI variables
    int rank, size;

    //MPI common initial routines
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Our MPI datatype
    int nitems = 3;
    MPI_Datatype types[nitems];
    MPI_Datatype MPI_PARTICLE;
    MPI_Aint offsets[nitems];
    int blocklenghts[nitems];

    types[0] = MPI_DOUBLE;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;

    offsets[0] = offsetof(particle,p);
    offsets[1] = offsetof(particle,v);
    offsets[2] = offsetof(particle,E);

    blocklenghts[0] = 3;
    blocklenghts[1] = 3;
    blocklenghts[2] = 1;

    MPI_Type_create_struct(nitems,blocklenghts,offsets,types,&MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);

    //We need to read these variables in order from our file
    //Files declaration
    MPI_File ifp;
    char inputFileName[] = "initialConditions.ic";
    MPI_File_open(MPI_COMM_WORLD, inputFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &ifp);

    MPI_File_read(ifp, &floatPointPrecision, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read(ifp, &n, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read(ifp, &numOfFiles, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read(ifp, &time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

    #if defined(VERBOSE)
        if(rank == 0){
            printf("floatPointPrecision: %d\n", floatPointPrecision);
            printf("n: %d\n", n);
            printf("numOfFiles: %d\n", numOfFiles);
            printf("time: %f\n", time);
        }
    #endif

    //Initializing the particles from the ic file proper to each process
    const double m = 100./(double) n;
    const int particlePerProcess = n/size + ((n%size) > rank);
    const int quotient = n/size;
    const int remainder = n%size;

    //We first read all the particles
    particle *pars;
    pars = readParticles(rank, size, n, &ifp, &MPI_PARTICLE);

    //MPI_File_close(&ifp);

    #if defined(DEBUG)
        sleep(rank*1);
        printf("\nrank %d i.c.\n", rank);
        for(int q = 0; q < n; ++q){
            printf("Position %d:\t%f; %f; %f\n", q, pars[q].p[0], pars[q].p[1], pars[q].p[2]);
            printf("Velocity %d:\t%f; %f; %f\n", q, pars[q].v[0], pars[q].v[1], pars[q].v[2]);
            printf("Energy %d:\t%f\n", q, pars[q].E);
            printf("\n\n");
        }
    #endif

    //Each process number of particles, so we can use allgatherv
    int *particlesNumAllProcesses;
    int *particlesDisplacement;

    if(( particlesNumAllProcesses = malloc(size*sizeof(int))) == NULL){
        printf("There is not enough space for particlesNumAllProcesses\n");
        exit(-1);
    }

    if(( particlesDisplacement = malloc(size*sizeof(int))) == NULL){
        printf("There is not enough space for particlesNumAllProcesses\n");
        exit(-1);
    }

    for(int register i = 0; i < size; ++i){
        particlesNumAllProcesses[i] = quotient + (i < remainder);
        particlesDisplacement[i] = quotient*i + ((i < remainder) ? i : remainder);
    }

    
    // #if defined(DEBUG)
    //     sleep(rank*1);
    //     printf("Rankd %d All processes number of particles\n", rank);
    //     for(int q = 0; q < size; ++q){
    //         printf("particlesNumAllProcesses[%d]: %d\n", q, particlesNumAllProcesses[q]);
    //         printf("particlesDisplacement[%d]: %d\n", q, particlesDisplacement[q]);
    //     }
    // #endif

    nBodyEvolution(pars, rank, size, particlesDisplacement[rank], particlesNumAllProcesses[rank], n, m);

    particle *newPars;

    if((newPars = malloc(n*sizeof(particle))) == NULL){
        printf("Not enough space to allocate the new particles\n");
        exit(-1);
    }

    //Particles redistribution
    MPI_Allgatherv(pars+particlesDisplacement[rank], particlePerProcess, MPI_PARTICLE, newPars, 
    particlesNumAllProcesses,particlesDisplacement, MPI_PARTICLE, MPI_COMM_WORLD);

    #if defined(DEBUG)
        sleep(rank*1);
        printf("\nrank %d all particles\n", rank);
        for(int q = 0; q < n; ++q){
            printf("Position %d:\t%f; %f; %f\n", q, newPars[q].p[0], newPars[q].p[1], newPars[q].p[2]);
            printf("Velocity %d:\t%f; %f; %f\n", q, newPars[q].v[0], newPars[q].v[1], newPars[q].v[2]);
            printf("Energy %d:\t%f\n", q, newPars[q].E);
            printf("\n\n");
        }
    #endif

    MPI_Finalize();

}

//We read all the particles on every MPI process
particle * readParticles(int rank, int size, int n, MPI_File *ifp, MPI_Datatype *MPI_PARTICLE){

    particle *pars;

    if((pars = malloc(n* sizeof(particle))) == NULL){
        printf("No memory available for particles of rank %d, sorry\n", rank);
        exit(-1);
    }

    //We take into consideration the header that we've already read
    int fileOffset = 3*sizeof(int) + sizeof(double);

    MPI_File_set_view(*ifp, fileOffset, *MPI_PARTICLE, *MPI_PARTICLE, "native", MPI_INFO_NULL);
    MPI_File_read(*ifp, pars, n, *MPI_PARTICLE, MPI_STATUS_IGNORE);

    return pars;
}

//Every MPI process updates its portion of the particles array
double nBodyEvolution(particle *par, int rank, int size, int displacement, int offset, int n, double m){

    //Unreasonably high first value so that we correctly calculate the min
    double tMin = 10;
    double t;

    //We need to save the forces for each particles. Having them in the struct would have meant
    //Useless movement of data (which is a costly operation) between processes, which is something
    //We want to avoid. We don't need the force outside of this function and the function that 
    //is called at the end.
    double *Fx;
    double *Fy;
    double *Fz;

    if((Fx = malloc(offset*sizeof(double))) == NULL){
        printf("Not enough space to allocate Fx\n");
        exit(-1);
    }

    if((Fy = malloc(offset*sizeof(double))) == NULL){
        printf("Not enough space to allocate Fx\n");
        exit(-1);
    }

    if((Fz = malloc(offset*sizeof(double))) == NULL){
        printf("Not enough space to allocate Fx\n");
        exit(-1);
    }
 
    for(int i = displacement; i < (displacement + offset); i++){    
        for( int k = 0; k < n; k++){
            //TODO: change this to k != i ? :
            if( k != i){
                Fx[i-displacement] += G*m*m*sgn(par[i].p[0], par[k].p[0]) / ( (par[i].p[0] - par[k].p[0]) * (par[i].p[0] - par[k].p[0]) + gamma/(double) n);
                Fy[i-displacement] += G*m*m*sgn(par[i].p[1], par[k].p[1]) / ( (par[i].p[1] - par[k].p[1]) * (par[i].p[1] - par[k].p[1]) + gamma/(double) n);
                Fz[i-displacement] += G*m*m*sgn(par[i].p[2], par[k].p[2]) / ( (par[i].p[2] - par[k].p[2]) * (par[i].p[2] - par[k].p[2]) + gamma/(double) n);
            }

        }

        //We have to take the minimum between the maximum of times
        //Each time is calculated as || V || * eps / || a |||
        //TODO:maybe there is space for optimization even here just considering the squares
        t = epsilon * m *
            sqrt(par[i].v[0]*par[i].v[0] + par[i].v[1]*par[i].v[1] + par[i].v[2]*par[i].v[2]) /
            sqrt(Fx[i]*Fx[i] + Fy[i]*Fy[i] + Fz[i]*Fz[i]);

        if( tMin > t)
            tMin = t;

    }

    //Here we need to decide the global time step
    //Each process sends his time to the master that gathers them in an array.
    //The master selects the minimum and broadcast it to the slaves
    double *times;

    if((times = malloc(size*sizeof(double))) == NULL){
        printf("Not enough space to allocate times\n");
        exit(-1);
    }

    MPI_Gather(&tMin, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #if defined(VERBOSE)
        if(rank == 0)
            for(int q = 0; q < size; ++q)
                printf("times[%d]: %f\n", q, times[q]);
    #endif

    double globalMin = 10;

    if(rank == 0){
        for(int register i = 0; i < size; ++i){
            if(globalMin > times[i])
                globalMin = times[i];
        }  
    }

    MPI_Bcast(&globalMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #if defined(VERBOSE)
        printf("Rank %d here, my global min is: %f\n", rank, globalMin);
    #endif

    updateParticles(par, displacement, offset, m, globalMin, Fx, Fy, Fz);

    free(Fx);
    free(Fy);
    free(Fz);

    return tMin;
}

void updateParticles(particle *par, int displacement, int offset, double m, double t, double *Fx, double *Fy, double *Fz){

    //v = a*t + v0
    for( int register i = displacement; i < (displacement+offset); i++){
        par[i].v[0] += Fx[i-displacement]/m*t;
        par[i].v[1] += Fy[i-displacement]/m*t;
        par[i].v[2] += Fz[i-displacement]/m*t; 
    }

    //s  = v*t + s0
    for( int register i = displacement; i < (displacement+offset); i++){
        par[i].p[0] += par[i].v[0]*t;
        par[i].p[1] += par[i].v[1]*t;
        par[i].p[2] += par[i].v[2]*t; 
    }     
}

//Auxiliary function, dunno if we need this
int sgn(double a, double b){
    if(a > b)
        return 1;
    else if( a  < b)
        return -1;    
}