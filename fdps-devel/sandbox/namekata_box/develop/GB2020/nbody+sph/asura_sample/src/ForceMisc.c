#include "config.h"
#include "ForceMisc.h"

void ClearGravitationalForce(void){

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            Pbody[i]->Acc[0] = 
            Pbody[i]->Acc[1] = 
            Pbody[i]->Acc[2] = 
            Pbody[i]->Pot = 0.e0;
            Pbody[i]->InteractionList = 0;
        }
    }
    return;
}

void ClearLocalAccPot(void){

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            Pbody[i]->Acc[0] = 
            Pbody[i]->Acc[1] = 
            Pbody[i]->Acc[2] = 
            Pbody[i]->Pot = 0.e0;
            Pbody[i]->InteractionList = 0;
        }
    }
    return;
}

void ForceEndProcedure(void){

#ifdef COSMOLOGICAL_RUN //{
    if(Pall.OmegaL < TINY)
        return;

    double fact = Pall.OmegaL*SQ(Pall.Hubble);
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            Pbody[i]->Acc[0] += fact*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] += fact*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] += fact*Pbody[i]->PosP[2];
        }
    }
#endif // COSMOLOGICAL_RUN //}
    return;
}

double ReturnInteractionListInDouble(void){

    unsigned long long int SumInteraction = 0;
    for(int i=0;i<Pall.Ntotal;i++)
        SumInteraction += Pbody[i]->InteractionList;
    double dSumInteraction = (double)SumInteraction;
    double GlobaldSumInteraction;

    MPI_Allreduce(&dSumInteraction,&GlobaldSumInteraction,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    return GlobaldSumInteraction;
}

void PrintForceFlops(double Time){

    double dSumInteraction = ReturnInteractionListInDouble();
#define ForceCount  (38)
    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"Total count = %g, %g GFlops\n",
            ForceCount*dSumInteraction,1.e-9*ForceCount*dSumInteraction/Time);
    }

    return ;
}
