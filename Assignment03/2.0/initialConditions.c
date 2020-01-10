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
void particleWrite(int rank, int size, particle *par, int n, int procParticles);
particle *domainDecomposition(int rank, int size, particle *par, int n, int procParticles);

int main(int argc, char **argv){

    FILE * ofp;
    char outputFileName[50];
    
    int N, m;

    if(argc < 2){
        printf("Using default number of particles\n\n");
        N = default;
    } else{
        N = atoi(argv[1]);
    }

    m = 100/N;

    int rank, size, remainder;
    
    //MPI common initial routines
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Define number of variables
    const int ownParsNum = N/size + (rank < N%size);
    const int otherParsNum = N - ownParsNum;

    //Defining our MPI_DataType
    // int nitems = 4;
    // MPI_Datatype types[nitems];
    // MPI_Datatype MPI_PARTICLE;
    // MPI_Aint offsets[nitems];
    // int blocklenghts[nitems];

    // types[0] = MPI_DOUBLE;
    // types[1] = MPI_DOUBLE;
    // types[2] = MPI_DOUBLE;
    // types[3] = MPI_DOUBLE;

    // offsets[0] = offsetof(particle,p);
    // offsets[1] = offsetof(particle,v);
    // offsets[2] = offsetof(particle,F);
    // offsets[3] = offsetof(particle,E);

    // blocklenghts[0] = 3;
    // blocklenghts[1] = 3;
    // blocklenghts[2] = 3;
    // blocklenghts[3] = 1;

    // MPI_Type_create_struct(nitems,blocklenghts,offsets,types,&MPI_PARTICLE);
    // MPI_Type_commit(&MPI_PARTICLE);

    //Per process parameter initialization
    
    snprintf(outputFileName, sizeof(outputFileName), "Debug_rank_%d.txt", rank);
    double x_border_left = 1. / size * rank;
    double x_border_right = 1. / size * (rank + 1);
    int seed = defaultSeed * (rank+1);
    
    //1. Each process writes n/p particles to the file and stores them    

    particle *ownPars = particleInitialize(N, size, m, rank, x_border_left, x_border_right);

    MPI_Finalize();

}


particle * particleInitialize(int n, int size, double m, int rank, double x_border_left, double x_border_right){

    //Let's set a proper seed
    int seed = rank*10 + 1; 
    const int rand_max = RAND_MAX;
    //This is the number of particles that every process has to work on
    int procParticles = n/size + (rank < n%size);

    #if defined(DEBUG)
        printf("Rank %d here. This is my number of particles: %d\n", rank, procParticles);
    #endif

	//Array allocation
    particle *par;
    if (( par = malloc(procParticles*sizeof(particle))) == NULL ){
        printf("I am sorry, we ran out of memory\n");
        exit(-1);
    };

    for(int register i = 0; i < procParticles; i++){
      
        //Random position values in [-1;1]
        // par[i].p[0] = seed;
        // par[i].p[1] = seed;
        // par[i].p[2] = seed;
        // par[i].v[0] = seed;
        // par[i].v[1] = seed;
        // par[i].v[2] = seed;
        // par[i].F[0] = seed;
        // par[i].F[1] = seed;
        // par[i].F[2] = seed;
        // par[i].E = seed;
        //Random position values [0;1]
        par[i].p[0] = (double)rand_r(&seed)/rand_max;
        par[i].p[1] = (double)rand_r(&seed)/rand_max;
        par[i].p[2] = (double)rand_r(&seed)/rand_max;

        //Random velocity values [-0.05; 0.05]
        par[i].v[0] = ((double)rand_r(&seed)*2/rand_max - 1)/20;
        par[i].v[1] = ((double)rand_r(&seed)*2/rand_max - 1)/20;
        par[i].v[2] = ((double)rand_r(&seed)*2/rand_max - 1)/20;

        par[i].F[0] = 0;
        par[i].F[1] = 0;
        par[i].F[2] = 0;

        par[i].E = 1./2. * m * par[i].v[0] * par[i].v[0] + par[i].v[1] * par[i].v[1] + par[i].v[2] * par[i].v[2];

        
    }

    particle *slicePar=domainDecomposition(rank, size, par, n, procParticles);

    return slicePar;

}

int compareParticles(const void *p1, const void *p2){
        particle *par1 = (particle *)p1;
        particle *par2 = (particle *)p2;
        if(par1->p[0] - par2->p[0] > 0)
            return 1;
        else if(par1->p[0] - par2->p[0] < 0)
            return -1;
        else
            return 0;
}

particle *domainDecomposition(int rank, int size, particle *par, int n, int procParticles){

    //MPI custom type declaration
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


    //Let's order the variables
    qsort(par, procParticles,sizeof(particle),compareParticles);

    int *procOffsets;
    if((procOffsets = calloc(size, sizeof(int))) == NULL){
        printf("Not enough memory to allocate the array, sorry!\n");
        exit(-1);
    }

    //Calculate the right amount of particles to send to every process
    int i = 0;
    double slice = 1./size;

    for(int k = 0; k < size; ++k){
        while(par[i].p[0] < (k+1)*slice && i < procParticles){
            procOffsets[k]++;
            ++i;
        }
    }
    
    #if defined(DEBUG)
        for(int k = 0; k < size; ++k)
            printf("In rank %d, procOffset %d is: %d\n", rank, k, procOffsets[k]);
    #endif

    //Let's scatter the procOffsets array
    //In this way, each MPI process knows how many elements to expect from the first, second, ..., p
    //It's hard to use collectives since the offsetsArray is of dimension size instead of size-1
    //DONE:this send and receive is correct
    int *procReceiveOffsets;
    if((procReceiveOffsets = malloc(size*sizeof(int))) == NULL){
        printf("There is no space for the receive offsets array\n");
        exit(-1);
    }

    
    for(int k = 0; k < size; ++k)
        if(rank != k)
            MPI_Send(&procOffsets[k], 1, MPI_INT, k, rank + 50, MPI_COMM_WORLD);

    for(int k = 0; k < size; ++k){
        if(rank != k)
            MPI_Recv(&procReceiveOffsets[k], 1, MPI_INT, k, k+50 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        else
            procReceiveOffsets[k] = procOffsets[k]; 
    }

    #if defined(DEBUG)
        for(int k = 0; k < size; ++k)
            printf("In rank %d, procReceiveOffsets %d is: %d\n", rank, k, procReceiveOffsets[k]);
    #endif

    //These indexes are useful to send and receive the right portion of the array
    i = 0;
    int prevI = 0;

    for(int k = 0; k < size; ++k){

        particle *buffer;

        //I know that mallocing each time is not the most efficient way but no other idea for now
        if((buffer = malloc(procOffsets[k]*sizeof(particle))) == NULL){
            printf("Not enough memory to allocate the buffer, big sorry\n");
            exit(-1);
        }      
        
        //TODO:check indexes because the buffer is actually incorrect when sending to rank 1 from rank 0
        //Send
        while( i-prevI < procOffsets[k]){
            //Copy particles in the buffer
            buffer[i-prevI] = par[i];
            ++i;
        }

        if(k != rank)
            MPI_Send(buffer, procOffsets[k], MPI_PARTICLE, k, rank, MPI_COMM_WORLD);

        #if defined(DEBUG)
            sleep(1*rank);
            printf("Rank %d buffer\n", rank);
            for(int q = 0; q < procOffsets[k]; q++){
                printf("Position %d:\t%f; %f; %f\n", q, buffer[q].p[0], buffer[q].p[1], buffer[q].p[2]);
                printf("Velocity %d:\t%f; %f; %f\n", q, buffer[q].v[0], buffer[q].v[1], buffer[q].v[2]);
                printf("Force %d:\t%f; %f; %f\n", q, buffer[q].F[0], buffer[q].F[1], buffer[q].F[2]);
                printf("Energy %d:\t%f\n", q, buffer[q].E);
                printf("\n\n");
            }
        #endif

        prevI = i;

        free(buffer);
    }

    #if defined(DEBUG)
        printf("Rank %d has finished sending\n\n", rank);
    #endif

    //Index reset
    i = 0;
    prevI = 0;

    //We store all the other particles in a final array
    particle *sliceParticles;
    int sliceParticlesDim = 0;
    //This will be useful to trace where to start feeding our initial array in the final one
    int procDisplacement = 0; 

    //How many particles we are going to receive: reduce sum on the receive offsets array
    for(int k = 0; k < size; ++k)
        sliceParticlesDim += procReceiveOffsets[k];
   
    //How many particles do we need to skip to populate our array?
    for(int k = 0; k < rank; ++k)
        procDisplacement += procOffsets[k];

    if((sliceParticles = malloc(sliceParticlesDim*sizeof(particle))) == NULL){
            printf("Not enough memory to allocate the buffer, big sorry\n");
            exit(-1);
    }

    #if defined(DEBUG)
        printf("Rank %d procDisplacement: %d\n",rank, procDisplacement);
        printf("Rank %d sliceParticlesDim: %d\n",rank, sliceParticlesDim);
    #endif

    //Here we are going to populate the actual slice Array
    for(int k = 0; k < size; ++k){

        particle *recBuffer;

        //I know that mallocing each time is not the most efficient way but no other idea for now
        if((recBuffer = malloc(procReceiveOffsets[k]*sizeof(particle))) == NULL){
            printf("Not enough memory to allocate the buffer, big sorry\n");
            exit(-1);
        }

        //Receive
        if(rank != k){
            MPI_Recv(recBuffer, procReceiveOffsets[k], MPI_PARTICLE, k, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else{
            for(int j = 0; j < procReceiveOffsets[k]; ++j)
                recBuffer[j] = par[procDisplacement+j];
        }

        #if defined(DEBUG)
            sleep(1*rank);
            printf("Rank %d recBuffer\n", rank);
            for(int q = 0; q < procReceiveOffsets[k]; q++){
                printf("Position %d:\t%f; %f; %f\n", q, recBuffer[q].p[0], recBuffer[q].p[1], recBuffer[q].p[2]);
                printf("Velocity %d:\t%f; %f; %f\n", q, recBuffer[q].v[0], recBuffer[q].v[1], recBuffer[q].v[2]);
                printf("Force %d:\t%f; %f; %f\n", q, recBuffer[q].F[0], recBuffer[q].F[1], recBuffer[q].F[2]);
                printf("Energy %d:\t%f\n", q, recBuffer[q].E);
                printf("\n\n");
            }
        #endif

        while( i-prevI < procReceiveOffsets[k]){
            sliceParticles[i] = recBuffer[i-prevI];
            ++i;
        }

        prevI = i;

        free(recBuffer);
        //TODO:I am afraid I was making shallow copies all the time. The last one should not be shallow
        //We need to free some memory I guess
        //free(par) 
    }

    #if defined(DEBUG)
        sleep(1*rank);
        printf("Rank %d i.c.\n", rank);
        for(int i = 0; i < sliceParticlesDim; i++){
            printf("Position %d:\t%f; %f; %f\n", i, sliceParticles[i].p[0], sliceParticles[i].p[1], sliceParticles[i].p[2]);
            printf("Velocity %d:\t%f; %f; %f\n", i, sliceParticles[i].v[0], sliceParticles[i].v[1], sliceParticles[i].v[2]);
            printf("Force %d:\t%f; %f; %f\n", i, sliceParticles[i].F[0], sliceParticles[i].F[1], sliceParticles[i].F[2]);
            printf("Energy %d:\t%f\n", i, sliceParticles[i].E);
            printf("\n\n");
        }
    #endif
    

    //FILE WRITING SECTION
    //Here we will use MPI I/O to write on the same file the slices
    //Each process will occupy the position of the file given by its rank.
    //In order to arrange the offset, anybody needs to know how many particles are preceeding in the file.
    //Therefore we have no choice but to introduce yet another communication

    //Broadcast of the sliceParticlesDim info
    //MPI_Bcast(sliceParticlesDim,1,MPI_INT,rank,MPI_COMM_WORLD);

    int *allSlicesParticles;
    if((allSlicesParticles = malloc(size*sizeof(int))) == NULL){
        printf("No memory left for our communication array, big sorry!\n");
        exit(-1);
    }

    //This collective function come very handy (TODO: modify procOffsets with the same logic)
    MPI_Allgather(&sliceParticlesDim,1,MPI_INT,allSlicesParticles,1,MPI_INT,MPI_COMM_WORLD);

    #if defined(DEBUG)
        printf("rank %d allSlicesParticles:\n", rank);
        for(int k = 0; k < size; ++k)
            printf("\t%d\n",allSlicesParticles[k]);
    #endif

    MPI_File fhw;
	MPI_File_open(MPI_COMM_WORLD, "initialConditions.csv", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);

    int fileOffset = 0;

    for(int k = 0; k < rank; k++)
        fileOffset += allSlicesParticles[k];
 
    MPI_File_write_at(fhw, fileOffset*sizeof(particle), par, sliceParticlesDim, MPI_PARTICLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fhw);

    return sliceParticles;    

}

void particleWrite(int rank, int size, particle *par, int n, int procParticles){

    MPI_File fhw;
	MPI_File_open(MPI_COMM_WORLD, "ic.csv", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);


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

    int offset;

    if(rank < n%size)
        offset = rank * (n/size + 1);
    else
        offset = rank * n/size + n%size;

    #if defined(DEBUG)
        printf("Rank %d here. This is my offset: %d\n", rank, offset);
    #endif

    MPI_File_write_at(fhw, offset*sizeof(particle), par, procParticles, MPI_PARTICLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fhw);

}