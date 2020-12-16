#include "config.h"

void UpdateTotalNumber(void){

    unsigned long int Numbers[5],GlobalNumbers[5];

    Numbers[0] = Pall.Ntotal;
    Numbers[1] = Pall.NDM;
    Numbers[2] = Pall.Nhydro;
    Numbers[3] = Pall.Nstars;
    Numbers[4] = Pall.Nsink;

    MPI_Allreduce(Numbers,GlobalNumbers,5,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    Pall.Ntotal_t = GlobalNumbers[0];
    Pall.NDM_t    = GlobalNumbers[1];
    Pall.Nhydro_t = GlobalNumbers[2];
    Pall.Nstars_t = GlobalNumbers[3];
    Pall.Nsink_t  = GlobalNumbers[4];

    return;
}

void UpdateTotalActiveNumber(void){

    unsigned long int Numbers[6],GlobalNumbers[6];

    Numbers[0] = Pall.NActives;
    Numbers[1] = Pall.NActivesDM;
    Numbers[2] = Pall.NActivesHydro;
    Numbers[3] = Pall.NActivesStars;
    Numbers[4] = Pall.NActivesSink;
    Numbers[5] = Pall.NActivesAll;

    MPI_Allreduce(Numbers,GlobalNumbers,6,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    Pall.NActives_t      = GlobalNumbers[0];
    Pall.NActivesDM_t    = GlobalNumbers[1];
    Pall.NActivesHydro_t = GlobalNumbers[2];
    Pall.NActivesStars_t = GlobalNumbers[3];
    Pall.NActivesSink_t  = GlobalNumbers[4];
    Pall.NActivesAll_t   = GlobalNumbers[5];

    return;
}

void UpdateTotalStarNumber(void){

    unsigned long int Numbers[2];

    Numbers[0] = Pall.Ntotal;
    Numbers[1] = Pall.Nstars;

    MPI_Allreduce(MPI_IN_PLACE,Numbers,2,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    Pall.Ntotal_t = Numbers[0];
    Pall.Nstars_t = Numbers[1];

    return;
}
