#include "config.h"

void InitializeCommunicationOrder(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

#if DECOMPOSITION_TYPE == 0 //{
    int SendRecvRank[NProcs][NProcs];
    for(int i=0;i<NProcs;i++)
        SendRecvRank[0][i] = i;

    int levelmax = (int)log2(NProcs);
    for(int level=0;level<levelmax;level++){
        int remainder,division,header;
        remainder = MyID%(1<<level);
        division = MyID/(1<<level);

        header = (1<<level)*division + Parity(division)*(1<<level);
        for(int k=0;k<NProcs;k++){
            header = k+Parity(k/(1<<level))*(1<<level);
            for(int i=0;i<(1<<level);i++){
                SendRecvRank[(1<<level)+i][k] = SendRecvRank[i][header];
            }
        }
    }

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].SendRank = SendRecvRank[i+1][MyID];
        CommunicationTable[i].RecvRank = SendRecvRank[i+1][MyID];
    }
#elif DECOMPOSITION_TYPE == 1 //}//{
    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].SendRank = (MyID+i+1)%NProcs;
        CommunicationTable[i].RecvRank = (MyID-i-1<0)
                                            ?(MyID-i-1+NProcs)%NProcs
                                            :(MyID-i-1)%NProcs;
        /*
        fprintf(stderr,"[%d] %d %d\n",
                MPIGetMyID(),
                CommunicationTable[i].SendRank,
                CommunicationTable[i].RecvRank);
        */
    }
#endif // DECOMPOSITION_TYPE //}

    return;
}

void InitializeCommunicationTable(void){

    const int Nprocs = MPIGetNumProcs();
    CommunicationTable = malloc(sizeof(struct StructCommunicationTable)*Nprocs);

    InitializeCommunicationOrder();

    return;
}


#ifdef USE_BARYON_COMM //{

void InitializeHydroCommunicationOrder(void){

    if(MPI_HYDRO_COMM_WORLD == MPI_COMM_NULL)
        return ;

    int MyID = MPIGetHydroMyID();
    int NProcs = MPIGetHydroNumProcs();

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].HydroSendRank = (MyID+i+1)%NProcs;
        CommunicationTable[i].HydroRecvRank = (MyID-i-1<0)
                                             ?(MyID-i-1+NProcs)%NProcs
                                             :(MyID-i-1)%NProcs;
    }
    return;
}

void InitializeActiveHydroCommunicationOrder(void){

    if(MPI_ACTIVEHYDRO_COMM_WORLD == MPI_COMM_NULL)
        return ;

    int MyID = MPIGetActiveHydroMyID();
    int NProcs = MPIGetActiveHydroNumProcs();

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].ActiveHydroSendRank = (MyID+i+1)%NProcs;
        CommunicationTable[i].ActiveHydroRecvRank = (MyID-i-1<0)
                                                   ?(MyID-i-1+NProcs)%NProcs
                                                   :(MyID-i-1)%NProcs;
    }

    return;
}

void InitializeBaryonCommunicationOrder(void){

    if(MPI_BARYON_COMM_WORLD == MPI_COMM_NULL)
        return ;

    int MyID = MPIGetBaryonMyID();
    int NProcs = MPIGetBaryonNumProcs();

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].BaryonSendRank = (MyID+i+1)%NProcs;
        CommunicationTable[i].BaryonRecvRank = (MyID-i-1<0)
                                              ?(MyID-i-1+NProcs)%NProcs
                                              :(MyID-i-1)%NProcs;
    }

    return;
}

void InitializeActiveBaryonCommunicationOrder(void){

    if(MPI_ACTIVEBARYON_COMM_WORLD == MPI_COMM_NULL)
        return ;

    int MyID = MPIGetActiveBaryonMyID();
    int NProcs = MPIGetActiveBaryonNumProcs();

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].ActiveBaryonSendRank = (MyID+i+1)%NProcs;
        CommunicationTable[i].ActiveBaryonRecvRank = (MyID-i-1<0)
                                                    ?(MyID-i-1+NProcs)%NProcs
                                                    :(MyID-i-1)%NProcs;
    }

    return;
}

#endif // USE_BARYON_COMM //}
