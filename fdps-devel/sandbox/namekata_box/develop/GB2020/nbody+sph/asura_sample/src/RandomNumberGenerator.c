#include "config.h"
#include "RandomNumberGenerator.h"
#include <gsl/gsl_rng.h>

static bool FirstCall = true;
void InitializeRandomGenerator(unsigned long int RandomSeed){

    AllocateRandomGenerator();
    //RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RandomGenerator,RandomSeed);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(FirstCall){
            fprintf(stderr,"Random generator name : %s\n",gsl_rng_name(RandomGenerator));
            fprintf(stderr," Seed %ld\n",RandomSeed);
            fprintf(stderr," Max and min value of the random numbers: %ld <-> %ld\n",
                    gsl_rng_min(RandomGenerator),gsl_rng_max(RandomGenerator));
            FirstCall = false;
        }
    }
    
    return;
}

void AllocateRandomGenerator(void){
    RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
    return;
}

void ResetRandomSeedForRandomGenerator(unsigned long int RandomSeed){
    gsl_rng_set(RandomGenerator,RandomSeed);
}
