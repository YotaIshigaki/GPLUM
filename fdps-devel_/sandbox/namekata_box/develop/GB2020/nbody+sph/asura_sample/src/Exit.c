#include "config.h"

void ASURA_Exit(void){

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"ASURA successfully finished this run.\n");
        fprintf(stderr," Data of this run is as follows.\n");
    }

    return ;
}
