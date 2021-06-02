#include "config.h"

/*
 * This function removes particles which satisfy the criteria provided by a
 * function pointer.
 */ 
void ParticleRemoverArbitraryCriteria(bool (*criteria)(const int index)){

    int NLoop = Pall.Ntotal;
    int counter = 0;
    for(int i=0;i<NLoop;i++){
        if((criteria)(i) == true){
            Pbody[i]->Use = OFF;
            Pall.Ntotal --;
            if(Pbody[i]->Type == TypeHydro){
                PbodyHydro(i)->Use = OFF;
                Pall.Nhydro --;
            }else if(Pbody[i]->Type == TypeStar){
                PbodyStar(i)->Use = OFF;
                Pall.Nstars --;
            }else if(Pbody[i]->Type == TypeDM){
                Pall.NDM --;
            }else if(Pbody[i]->Type == TypeSink){
                Pall.Nsink --;
            }
            counter ++;
        }
    }

    int sum_counter = 0;
    MPI_Allreduce(&counter,&sum_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(sum_counter > 0){
        ReConnectPointers();
        UpdateTotalNumber();
    }

    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"%d particles are removed in this simulation.\n",sum_counter);
    }

    return ;
}
