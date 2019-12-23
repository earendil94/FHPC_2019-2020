#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define seed 100
#define G 10e-06
#define default 10

//This represents a particle position, velocity, force (in a 3d vector) and energy
typedef struct { 
    double p[3];
    double v[3];
    double F[3];
    double E;
} particle;

void particleInitialize(particle *par, int n, double m);
void readParticles(particle *par, int n);
void updateParticlesForce(particle *par, int n, double m);
void updateParticles(particle *par, int n, double m, double t);
int sgn(double a, double b);

int main(int argc, char **argv){

    int N;

    if(argc < 2){
        printf("Using default number of particles\n");
        N = default;
    } else
    {
        N = atoi(argv[1]);
    }
    
    //The mass of every particle is constant
    double m = 100./(double) N;

    //Particles allocation: maybe this will need to be two arrays for our parallel code
    particle *pars;

    if( (pars = malloc(N * sizeof(particle))) == NULL ){
        printf("Houston, we can't create this many particles!");
        exit(-1);
    }

    particleInitialize(pars, N, m);

    readParticles(pars, N);

    updateParticlesForce(pars, N, m);

    double t = 0.05;

    updateParticles(pars, N, m ,t);

    readParticles(pars,N);

}

void particleInitialize(particle *par, int n, double m){

    srand(seed);

    for( int register i = 0; i < n; i++){
      
        //Random position values
        par[i].p[0] = drand48();
        par[i].p[1] = drand48();
        par[i].p[2] = drand48();

        //Random velocity values (we need to refine the drand48 random generation here)
        par[i].v[0] = drand48();
        par[i].v[1] = drand48();
        par[i].v[2] = drand48();

        par[i].F[0] = 0;
        par[i].F[1] = 0;
        par[i].F[2] = 0;

        par[i].E = 1./2. * m * par[i].v[0] * par[i].v[0] + par[i].v[1] * par[i].v[1] + par[i].v[2] * par[i].v[2];
    }

}

void readParticles(particle *par, int n){

    for(int register i = 0; i < n; i++){
        printf("Particle[%d]:\n"
                "\t x: %f; y: %f; z: %f\n",
                i, par[i].p[0], par[i].p[1], par[i].p[2]);
    }
}


//This needs to update every property of every particle
//First big doubt: do we update the force of each variable AND THEN we make them move?
//Seems reasonable, otherwise we would be implicitly defining an order of variable,
//Which is something that does not actually exist
void updateParticlesForce(particle *par, int n, double m){

    //Still unsure of which is the best way to exclude the i-th component
    for(int register i = 0; i < n; i++){     
        for( int register k = 0; k < n; k++){

            if( k != i){
                par[i].F[0] += G*m*m*sgn(par[i].p[0], par[k].p[0]) / ( (par[i].p[0] - par[k].p[0]) * (par[i].p[0] - par[k].p[0]) );
                par[i].F[1] += G*m*m*sgn(par[i].p[1], par[k].p[1]) / ( (par[i].p[1] - par[k].p[1]) * (par[i].p[1] - par[k].p[1]) );
                par[i].F[2] += G*m*m*sgn(par[i].p[2], par[k].p[2]) / ( (par[i].p[2] - par[k].p[2]) * (par[i].p[2] - par[k].p[2]) );
            }
        }
    }
}

void updateParticles(particle *par, int n, double m, double t){

    //Second doubt: for the way we have defined our deltaV, we are now summing two times our initial v
    //I will suggest a very stupid update like the one that follows

    //v = a*t + v0
    for( int register i = 0; i < n; i++){
        par[i].v[0] += par[i].F[0]/m*t;
        par[i].v[1] += par[i].F[1]/m*t;
        par[i].v[2] += par[i].F[2]/m*t; 
    }

    //s  = v*t + s0
    for( int register i = 0; i < n; i++){
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