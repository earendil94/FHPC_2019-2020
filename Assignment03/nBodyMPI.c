#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>

#define defaultSeed 100
#define G 10e-06
#define default 10
#define epsilon 0.05

//This represents a particle position, velocity, force (in a 3d vector) and energy
typedef struct { 
    double p[3];
    double v[3];
    double F[3];
    double E;
} particle;

particle * particleInitialize(int n, int size, double m, int rank, double x_border_left, double x_border_right);
particle * readParticles(char *filename, int n, int size, int rank);

//Version 1, to be changed if further indications are given.
//The main idea behind the MPI process is: everybody generates initial condition of the
//particles and stores them, alongside writing them on a file. Then all the other processes
//reads the other particles of the universe and store them in a different array.
int main(int argc, char **argv){

    FILE * ofp;
    char outputFileName[50];
    

    int N, m;

    if(argc < 2){
        printf("Using default number of particles\n\n");
        N = default;
    } else
    {
        N = atoi(argv[1]);
    }

    m = 100/N;

    int rank, size, N_per_proc, remainder;

    //MPI common initial routines
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

        //Defining our MPI_DataType
    int nitems = 4;
    MPI_Datatype types[nitems];
    MPI_Datatype MPI_PARTICLE;
    MPI_Aint offsets[nitems];
    int blocklenghts[nitems];

    types[0] = MPI_DOUBLE;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;
    types[3] = MPI_DOUBLE;

    offsets[0] = offsetof(particle,p);
    offsets[1] = offsetof(particle,v);
    offsets[2] = offsetof(particle,F);
    offsets[3] = offsetof(particle,E);

    blocklenghts[0] = 3;
    blocklenghts[1] = 3;
    blocklenghts[2] = 3;
    blocklenghts[3] = 1;

    MPI_Type_create_struct(nitems,blocklenghts,offsets,types,&MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);

    //Per process parameter initialization
    
    snprintf(outputFileName, sizeof(outputFileName), "Debug_rank_%d.txt", rank);
    double x_border_left = 1. / size * rank;
    double x_border_right = 1. / size * (rank + 1);
    int seed = defaultSeed * (rank+1);
    
    //1. Each process writes n/p particles to the file and stores them    

    particle *ownPars = particleInitialize(N, size, m, rank, x_border_left, x_border_right);
    int ownParsNum = N/size + (rank < N%size);

    sleep(1*rank);
    printf("Rank %d i.c.\n", rank);
    for(int i = 0; i < ownParsNum; i++){
        printf("Position %d:\t%f; %f; %f\n", i, ownPars[i].p[0], ownPars[i].p[1], ownPars[i].p[2]);
        printf("Velocity %d:\t%f; %f; %f\n", i, ownPars[i].v[0], ownPars[i].v[1], ownPars[i].v[2]);
        printf("Force %d:\t%f; %f; %f\n", i, ownPars[i].F[0], ownPars[i].F[1], ownPars[i].F[2]);
        printf("Energy %d:\t%f\n", i, ownPars[i].E);
        printf("\n\n");
    }


    //2. Each process reads the rest of the particles in the file
    particle *otherPars = readParticles("ic.csv", N, size, rank);
    int otherParsNum = N - ownParsNum;

    //DEBUG
    //It seems that we write on the same position of the file and just once
    sleep(1*rank);
    printf("Hey I received these other ones\n");
    for(int i = 0; i < otherParsNum; i++){
        printf("Position %d:\t%f; %f; %f\n", i, otherPars[i].p[0], otherPars[i].p[1], otherPars[i].p[2]);
        printf("Velocity %d:\t%f; %f; %f\n", i, otherPars[i].v[0], otherPars[i].v[1], otherPars[i].v[2]);
        printf("Force %d:\t%f; %f; %f\n", i, otherPars[i].F[0], otherPars[i].F[1], otherPars[i].F[2]);
        printf("Energy %d:\t%f\n", i, otherPars[i].E);
        printf("\n\n");
    }

    MPI_Finalize();
}

//TODO: Still need to implement the random on the proper interval
particle * particleInitialize(int n, int size, double m, int rank, double x_border_left, double x_border_right){

    int seed = rank*10 + 1;
    int procParticles = n/size + (rank < n%size);

    printf("Rank %d here. This is my numberofPars: %d\n", rank, procParticles);

    MPI_File fhw;
    srand(seed);
	MPI_File_open(MPI_COMM_WORLD, "ic.csv", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);
	
    particle *par;
    
    if (( par = malloc(procParticles*sizeof(particle))) == NULL ){
        printf("I am sorry, we ran out of memory\n");
        exit(-1);
    };

    //Defining our MPI_DataType
    int nitems = 4;
    MPI_Datatype types[nitems];
    MPI_Datatype MPI_PARTICLE;
    MPI_Aint offsets[nitems];
    int blocklenghts[nitems];

    types[0] = MPI_DOUBLE;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;
    types[3] = MPI_DOUBLE;

    offsets[0] = offsetof(particle,p);
    offsets[1] = offsetof(particle,v);
    offsets[2] = offsetof(particle,F);
    offsets[3] = offsetof(particle,E);

    blocklenghts[0] = 3;
    blocklenghts[1] = 3;
    blocklenghts[2] = 3;
    blocklenghts[3] = 1;

    MPI_Type_create_struct(nitems,blocklenghts,offsets,types,&MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);

    for(int register i = 0; i < procParticles; i++){
      
        //Random position values in [-1;1]
        par[i].p[0] = seed;
        par[i].p[1] = seed;
        par[i].p[2] = seed;
        par[i].v[0] = seed;
        par[i].v[1] = seed;
        par[i].v[2] = seed;
        par[i].F[0] = seed;
        par[i].F[1] = seed;
        par[i].F[2] = seed;
        par[i].E = seed;
        // par[i].p[0] = drand48()*2 - 1;
        // par[i].p[1] = drand48()*2 - 1;
        // par[i].p[2] = drand48()*2 - 1;

        // //Random velocity values [-0.05; 0.05]
        // par[i].v[0] = (drand48()*2 - 1)/20;
        // par[i].v[1] = (drand48()*2 - 1)/20;
        // par[i].v[2] = (drand48()*2 - 1)/20;

        // par[i].F[0] = 0;
        // par[i].F[1] = 0;
        // par[i].F[2] = 0;

        // par[i].E = 1./2. * m * par[i].v[0] * par[i].v[0] + par[i].v[1] * par[i].v[1] + par[i].v[2] * par[i].v[2];

        
    }

    int offset = rank * (n/size + ((rank - 1) < n%size));

    printf("Rank %d here. This is my offset: %d\n", rank, offset);

    MPI_File_write_at(fhw, offset*sizeof(particle), par, procParticles, MPI_PARTICLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fhw);

    return par;

}


//Needs to read the extra chunck when p > 2
particle * readParticles(char *filename, int n, int size, int rank){

    MPI_File fhw;
    MPI_File_open(MPI_COMM_WORLD, "ic.csv", MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    //The particles we have to read are N - the number of particles already allocated 
    int procParticles = n - (n/size + (rank < n%size));

    particle *pars;
    if (( pars = malloc(procParticles*sizeof(particle))) == NULL ){
        printf("Apologize, no memory left for our particles\n");
        exit(-1);
    }

    //Defining our MPI_DataType
    int nitems = 4;
    MPI_Datatype types[nitems];
    MPI_Datatype MPI_PARTICLE;
    MPI_Aint offsets[nitems];
    int blocklenghts[nitems];

    types[0] = MPI_DOUBLE;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;
    types[3] = MPI_DOUBLE;

    offsets[0] = offsetof(particle,p);
    offsets[1] = offsetof(particle,v);
    offsets[2] = offsetof(particle,F);
    offsets[3] = offsetof(particle,E);

    blocklenghts[0] = 3;
    blocklenghts[1] = 3;
    blocklenghts[2] = 3;
    blocklenghts[3] = 1;

    MPI_Type_create_struct(nitems,blocklenghts,offsets,types,&MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);

    //For now
    if( rank == 0 )
        MPI_File_read_at(fhw, 6*sizeof(particle),pars, 5, MPI_PARTICLE, MPI_STATUS_IGNORE);
    else
        MPI_File_read_at(fhw, 0, pars, 6, MPI_PARTICLE, MPI_STATUS_IGNORE);
    //int offset = (rank+1) * sizeof(MPI_PARTICLE) * n/size;

    const int extraPiece = n/size + 1;
    const int regularPiece = n/size;

    int startRead;
    int stopRead;

    if(rank == 0){
        if( n%size == 0 )
            startRead = regularPiece;
        else
            startRead = extraPiece;

        stopRead = n - startRead;
    }
    else{
        startRead = 0;

        if( rank >= (n%size))
            stopRead = (n%size) * extraPiece + (rank - (n%size)) * regularPiece;
        else
            stopRead = (n%size) * extraPiece;
    }

    printf("Rank %d here, this is where I start: %d\n", rank, startRead);
    printf("Rank %d here, this is where I stop: %d\n", rank, stopRead);

    MPI_File_read_at(fhw, startRead*sizeof(particle),pars, stopRead, MPI_PARTICLE, MPI_STATUS_IGNORE);

    //Probably the best idea to read another chunck is to define a buffer and put our particles on it 
    //if(rank )

    MPI_File_close(&fhw);

    return pars;

}

