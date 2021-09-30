#include "config.h"
#include "PreDecomposition.h"
#include "Decomposition.h"
#include "ParallelOperation.h"

static void ParticlesExchangeLimited(void);
static void ParticlesExchangeLimitedNew(void);
static void DataUnPacking(int NContainerRecv[]);
static void DataPacking(const int SendRank, int NContainerSend[]);
static void DomainOutPutKey(void);
static void DomainDecompResult(char *name);

static int DomainKeyMaxAllocated = 0;
static int DomainKeyCurrentSize;
static int *DomainKey;

static void PrepareDomainKey(void){

    if(DomainKeyMaxAllocated < PbodyElementsSize){
        if(DomainKeyMaxAllocated > 0)
            free(DomainKey);
        DomainKeyMaxAllocated = (int)(ForAngelsShare*PbodyElementsSize);
        DomainKey = malloc(sizeof(int)*DomainKeyMaxAllocated);
    }

    return;
}

static void SetDecompositionKey(void){

#if DECOMPOSITION_TYPE == 0 //{
    int Power = MPIGetNumProcsPower();

    for(int i=0;i<PbodyElementsSize;i++){
        if(PbodyElements[i].Use == true){
            DomainKey[i] = 0;
            int index = 0;
            for(int k=0;k<Power;k++){
                //int Axis = k%BisectionDimension;
                int Axis = InfoBiSection[index].Axis;
                if(PbodyElements[i].PosP[Axis] > InfoBiSection[index].Pos){
                    index = InfoBiSection[index].Right;
                    DomainKey[i] += 1<<k;
                }else{
                    index = InfoBiSection[index].Left;
                }
            }
        } else {
            DomainKey[i] = -1;
        }
    }
#elif DECOMPOSITION_TYPE == 1 //}//{
    for(int i=0;i<PbodyElementsSize;i++){
        if(PbodyElements[i].Use == true){
            DomainKey[i] = GetDomainKeyXYZ(PbodyElements[i].PosP);
        }else{
            DomainKey[i] = -1;
        }
    }
#endif // DECOMPOSITION_TYPE //}
    DomainKeyCurrentSize = PbodyElementsSize;
    return;
}

static int DecompositionIterationCount;

void DomainDecompositionOnlyDataExchange(void){

    double TimingResultThisRoutine = GetElapsedTime();

    DecompositionIterationCount = 0;

    if(MPIGetNumProcs() == 1) {
#ifdef USE_BARYON_COMM //{
        DeleteBaryonGroup();
        CreateBaryonGroup();
        DeleteHydroGroup();
        CreateHydroGroup();
#endif // USE_BARYON_COMM //}
        return;
    }


    // dprintlmpi(DecompositionIterationCount);
    // fprintf(stderr,"[%02d] %ld\n",MPIGetMyID(),Pall.Ntotal);
    // fflush(NULL);
    // MPI_Barrier(MPI_COMM_WORLD);

    /*
    char fname[MaxCharactersInLine];
    sprintf(fname,"DomainDecompInit.%02d",MPIGetMyID());
    DomainDecompResult(fname);
    sprintf(fname,"DomainAfterBisection.%02d",MPIGetMyID());
    DomainSamplesOutPutAfterBiSection(fname);
    */

    PrepareDomainKey();

    SetDecompositionKey();

    //ParticlesExchangeLimited();
    ParticlesExchangeLimitedNew();
    //ParticlesExchange();
    //EstimateNextUpdateTime();

#ifdef USE_BARYON_COMM //{
    DeleteBaryonGroup();
    CreateBaryonGroup();
    DeleteHydroGroup();
    CreateHydroGroup();
#endif //USE_BARYON_COMM //}

    TimingResults.DecompositionThisStep += GetElapsedTime()-TimingResultThisRoutine;
    return;
}

void DomainDecomposition(void){

    if(MPIGetNumProcs() == 1) 
        return;

    //if(Pall.DecompositionNextUpdateTime > Pall.TCurrent) 
    //   return;

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Do Domain Decomposition\n");

    //PreDomainDecomposition();
    PreDomainDecomposition(0);
    DomainDecompositionOnlyDataExchange();

    return;
}

StructPbodyptr DecompositionSend;
StructPbodyptr DecompositionRecv;

//#define ExchangeLimit   (200000)
#define ExchangeLimit   (50000)
static void ParticlesExchangeLimited(void){

    if(DecompositionIterationCount > DECOMPOSITION_MAX_ITERATION_TIMES)
        return;

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    int NContainerSend[NProcs][NTypes];
    int NContainerRecv[NProcs][NTypes];
    for(int i=0;i<NProcs;i++)
        for(int k=0;k<NTypes;k++)
            NContainerSend[i][k] = NContainerRecv[i][k] = 0;

    int NotMyID = 0;


    int counter = 0;
    for(int i=0;i<PbodyElementsSize;i++){
        if(DomainKey[i] != -1){ // skip empty structures.
            if(counter < ExchangeLimit){
                NContainerSend[DomainKey[i]][0] ++; // Body
                if(PbodyElements[i].Type == TypeHydro){ // Hydro
                    NContainerSend[DomainKey[i]][TypeHydro] ++;
                }else if(PbodyElements[i].Type == TypeStar){ // Star
                    NContainerSend[DomainKey[i]][TypeStar] ++;
                }else if(PbodyElements[i].Type == TypeSink){ // Sink
                    NContainerSend[DomainKey[i]][TypeSink] ++;
                }
            } else if(DomainKey[i] != MyID){
                NotMyID ++;
            }
            counter ++;
        }
    }


#if 0
//////////////////////////////////////////////////////////////////////

    int ExchangeLimitThisStep = ExchangeLimit/NProcs;

    int counter = 0;
    for(int i=0;i<PbodyElementsSize;i++){
        if(DomainKey[i] != -1){ // skip empty structures.
            if(NContainerSend[DomainKey[i]][0] < ExchangeLimitThisStep){
                NContainerSend[DomainKey[i]][0] ++; // Body
                if(PbodyElements[i].Type == TypeHydro){ // Hydro
                    NContainerSend[DomainKey[i]][TypeHydro] ++;
                }else if(PbodyElements[i].Type == TypeStar){ // Star
                    NContainerSend[DomainKey[i]][TypeStar] ++;
                }else if(PbodyElements[i].Type == TypeSink){ // Sink
                    NContainerSend[DomainKey[i]][TypeSink] ++;
                }
            } else if(DomainKey[i] != MyID){
                NotMyID ++;
            }
            counter ++;
        }
    }
#endif


    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].SendSize 
            = sizeof(StructPbody)*NContainerSend[CommunicationTable[i].SendRank][0]
             +sizeof(StructPhydro)*NContainerSend[CommunicationTable[i].SendRank][TypeHydro]
             +sizeof(StructPstar)*NContainerSend[CommunicationTable[i].SendRank][TypeStar]
             +sizeof(StructPsink)*NContainerSend[CommunicationTable[i].SendRank][TypeSink];
    }


    int NSendAllocated = 0;
    int NRecvAllocated = 0;
    int Stride = (NTypes+1);
    int SendSizes[NProcs*Stride];
    int RecvSizes[NProcs*Stride];
    for(int i=0;i<NProcs-1;i++){
        SendSizes[Stride*CommunicationTable[i].SendRank] = CommunicationTable[i].SendSize;
        for(int k=0;k<NTypes;k++){
            SendSizes[Stride*CommunicationTable[i].SendRank+k+1] =
                NContainerSend[CommunicationTable[i].SendRank][k];
        }
        NSendAllocated = MAX(CommunicationTable[i].SendSize,NSendAllocated);
    }
    MPI_Alltoall(SendSizes,Stride,MPI_INT,RecvSizes,Stride,MPI_INT,MPI_COMM_WORLD);


    int NotMyIDAll = NotMyID;
    MPI_Allreduce(MPI_IN_PLACE,&NotMyIDAll,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    int NImportAll = 0;
    int NContainerRecvAll[NTypes];
    memset(NContainerRecvAll,0,sizeof(int)*NTypes);

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].RecvSize = RecvSizes[Stride*CommunicationTable[i].RecvRank];
        for(int k=0;k<NTypes;k++){
            NContainerRecv[i][k] = RecvSizes[Stride*CommunicationTable[i].RecvRank+k+1];
            NContainerRecvAll[k] = NContainerRecv[i][k];
        }
        NRecvAllocated = MAX(CommunicationTable[i].RecvSize,NRecvAllocated);
        NImportAll += CommunicationTable[i].RecvSize;
    }


    StructPbodyptr DataPackSend[NProcs];
    StructPbodyptr DataPackRecv[NProcs];
    for(int i=0;i<NProcs-1;i++){
        CheckSizeofBufferExportSendIndex(CommunicationTable[i].SendSize,1,i);
        CheckSizeofBufferImportRecvIndex(CommunicationTable[i].RecvSize,1,i);
        DataPackSend[i] = BufferExportSend[i];
        DataPackRecv[i] = BufferImportRecv[i];
    }

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    size_t NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        DecompositionSend = DataPackSend[i];
        DataPacking(CommunicationTable[i].SendRank,
                NContainerSend[CommunicationTable[i].SendRank]);

        if(CommunicationTable[i].SendSize>0){
            MPI_Isend(DataPackSend[i],CommunicationTable[i].SendSize,
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_DECOMPOSITION_SENDDATA+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(CommunicationTable[i].RecvSize>0){
            MPI_Irecv(DataPackRecv[i],CommunicationTable[i].RecvSize,
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_DECOMPOSITION_SENDDATA+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += CommunicationTable[i].RecvSize;
    }

    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

    for(int i=0;i<NProcs-1;i++){
        if(CommunicationTable[i].RecvSize>0){
            DecompositionRecv = DataPackRecv[i];
            DataUnPacking(NContainerRecv[i]);
        }
    }


    ReConnectPointers();


    if(NotMyIDAll != 0){
        DecompositionIterationCount ++;
        PrepareDomainKey();
        SetDecompositionKey();
        ParticlesExchangeLimited();
    }

    return;
}

static void DataPackingNew(const int SendRank, int NContainerSend[restrict], void *DecompositionExportSend){

    StructPbodyptr  PbSend = DecompositionExportSend;
    StructPhydroptr PhSend = (StructPhydroptr)(PbSend+NContainerSend[0]);
    StructPstarptr  PsSend = (StructPstarptr)(PhSend+NContainerSend[TypeHydro]);
    StructPsinkptr  PskSend = (StructPsinkptr)(PsSend+NContainerSend[TypeStar]);

    int nbody = 0;
    int nhydro = 0;
    int nstars = 0;
    int nsinks = 0;

    for(int i=0;i<PbodyElementsSize;i++){
        if(DomainKey[i] == SendRank){
            StructPbodyptr Pb = PbodyElements+i;

            PbSend[nbody] = *Pb;
            Pb->Use = false;
            nbody ++;
            if(Pb->Type == TypeHydro){
                PhSend[nhydro] = *((StructPhydroptr)(Pb->Baryon));
                ((StructPhydroptr)(Pb->Baryon))->Use = false;
                nhydro ++;
                Pall.Nhydro --;
            }else if(Pb->Type == TypeStar){
                PsSend[nstars] = *((StructPstarptr)(Pb->Baryon));
                ((StructPstarptr)(Pb->Baryon))->Use = false;
                nstars ++;
                Pall.Nstars --;
            }else if(Pb->Type == TypeSink){
                PskSend[nsinks] = *((StructPsinkptr)(Pb->Baryon));
                ((StructPsinkptr)(Pb->Baryon))->Use = false;
                nsinks ++;
                Pall.Nsink --;
            }else if(Pb->Type == TypeDM){
                Pall.NDM --;
            }
        }
        if(nbody == NContainerSend[0]){
            Pall.Ntotal -= nbody;
            return;
        }
    }

    return ;
}

static void DataUnPackingNew(int NContainerRecv[restrict], void *DecompositionExportRecv){

    StructPbodyptr  PbRecv = DecompositionExportRecv;
    StructPhydroptr PhRecv = (StructPhydroptr)(PbRecv+NContainerRecv[0]);
    StructPstarptr  PsRecv = (StructPstarptr)(PhRecv+NContainerRecv[TypeHydro]);
    StructPsinkptr  PskRecv = (StructPsinkptr)(PsRecv+NContainerRecv[TypeStar]);

    if(NContainerRecv[0] == 0)
        return;

    // count unuse
    int counter = CountEmptyPbody();
    if(NContainerRecv[0]-counter > 0){
        if(PbodySize == 0){
            GenerateStructPbody(NContainerRecv[0]-counter);
        }else{
            StretchStructPbody(NContainerRecv[0]-counter);
        }
    }

    // Hydro array.
    counter = CountEmptyPhydro();
    if(NContainerRecv[TypeHydro]-counter > 0){
        if(PhydroSize == 0){
            GenerateStructPhydro(NContainerRecv[TypeHydro]-counter);
        }else{
            StretchStructPhydro(NContainerRecv[TypeHydro]-counter);
        }
    }

    // Star array.
    counter = CountEmptyPstar();
    if(NContainerRecv[TypeStar]-counter > 0){
        if(PstarSize == 0){
            GenerateStructPstar(NContainerRecv[TypeStar]-counter);
        }else{
            StretchStructPstar(NContainerRecv[TypeStar]-counter);
        }
    }

    // Sink array.
    counter = CountEmptyPsink();
    if(NContainerRecv[TypeSink]-counter > 0){
        if(PsinkSize == 0){
            GenerateStructPsink(NContainerRecv[TypeSink]-counter);
        }else{
            StretchStructPsink(NContainerRecv[TypeSink]-counter);
        }
    }


    int nbody = 0;
    int nhydro = 0;
    int nstars = 0;
    int nsinks = 0;
    StructPhydroptr Ph = PhydroElements;
    StructPstarptr Ps = PstarElements;
    StructPsinkptr Psk= PsinkElements;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == false){
            StructPbodyptr PbNext = Pb->Next;
            *Pb = PbRecv[nbody];
            Pb->Next = PbNext;
            Pb->Use = true;
            nbody ++;
            Pall.Ntotal ++;

            if(Pb->Type == TypeHydro){
                for(StructPhydroptr Pt=Ph; Pt; Pt = Pt->Next){
                    if(Pt->Use == false){
                        StructPhydroptr PhNext = Pt->Next;
                        *Pt = PhRecv[nhydro];
                        Pt->Next = PhNext;
                        Pt->Use = true;
                        Pb->Baryon = Pt;
                        Pt->Body = Pb;
                        Ph = Pt->Next;
                        break;
                    }
                }
                nhydro ++;
                Pall.Nhydro ++;
            }else if(Pb->Type == TypeStar){
                for(StructPstarptr Pt=Ps; Pt; Pt = Pt->Next){
                    if(Pt->Use == false){
                        StructPstarptr PsNext = Pt->Next;     
                        *Pt = PsRecv[nstars];
                        Pt->Next = PsNext;
                        Pt->Use = true;
                        Pb->Baryon = Pt;
                        Pt->Body = Pb;
                        Ps = Pt->Next;
                        break;
                    }
                }
                nstars ++;
                Pall.Nstars ++;
            }else if(Pb->Type == TypeSink){
                for(StructPsinkptr Pt=Psk; Pt; Pt = Pt->Next){
                    if(Pt->Use == false){
                        StructPsinkptr PskNext = Pt->Next;     
                        *Pt = PskRecv[nsinks];
                        Pt->Next = PskNext;
                        Pt->Use = true;
                        Pb->Baryon = Pt;
                        Pt->Body = Pb;
                        Psk = Pt->Next;
                        break;
                    }
                }
                nsinks ++;
                Pall.Nsink ++;
            } else if (Pb->Type == TypeDM){
                Pall.NDM ++;
            }
        }
        if(nbody == NContainerRecv[0])
            break;
    }

    /*
    fprintf(stderr,"UnPacking [%02d] %d %d %d %d | %d %d %d %d\n",MPIGetMyID(),nbody,nhydro,nstars,nsinks,
            NContainerRecv[0],NContainerRecv[1],NContainerRecv[2],NContainerRecv[3]);

    fflush(NULL);
    if(nhydro != NContainerRecv[TypeHydro]){
        fprintf(stderr,"[%02d] %d %d -- \n",MPIGetMyID(),nhydro,NContainerRecv[TypeHydro]);
        fflush(NULL);
        assert(nhydro == NContainerRecv[TypeHydro]);

    }
    fflush(NULL);
    */
    assert(nbody == NContainerRecv[0]);
    assert(nhydro == NContainerRecv[TypeHydro]);
    assert(nstars == NContainerRecv[TypeStar]);
    assert(nsinks == NContainerRecv[TypeSink]);
    return ;
}

static void ParticlesExchangeLimitedNew(void){

    // dprintlmpi(DecompositionIterationCount);
    if(DecompositionIterationCount > DECOMPOSITION_MAX_ITERATION_TIMES)
        return;

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    int NContainerSend[NProcs][NTypes];
    for(int i=0;i<NProcs;i++)
        for(int k=0;k<NTypes;k++)
            NContainerSend[i][k] = 0;

    // counter for export and import. 
    int NExportThisTime[NProcs][NTypes+1];
    int NImportThisTime[NProcs][NTypes+1];

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    int SendFlag,RecvFlag;

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(NImportThisTime[i],NTypes+1,MPI_INT,
            CommunicationTable[i].RecvRank,TAG_DECOMPOSITION_SENDSIZE+i,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
    }

    int NotMyRank = 0;
    int counter = 0;
    for(int i=0;i<PbodyElementsSize;i++){
        if(DomainKey[i] != -1){ // skip empty structures.
            if((DomainKey[i] != MyID)&&(counter < ExchangeLimit)){
                NContainerSend[DomainKey[i]][0] ++; // Body
                if(PbodyElements[i].Type == TypeHydro){ // Hydro
                    NContainerSend[DomainKey[i]][TypeHydro] ++;
                }else if(PbodyElements[i].Type == TypeStar){ // Star
                    NContainerSend[DomainKey[i]][TypeStar] ++;
                }else if(PbodyElements[i].Type == TypeSink){ // Sink
                    NContainerSend[DomainKey[i]][TypeSink] ++;
                }
                counter ++;
            } else if(DomainKey[i] != MyID){
                NotMyRank ++;
            }
        }
    }

    void *DecompositionExportSend[NProcs];

    for(int i=0;i<NProcs-1;i++){
        CommunicationTable[i].SendSize 
            = sizeof(StructPbody) *NContainerSend[CommunicationTable[i].SendRank][0]
             +sizeof(StructPhydro)*NContainerSend[CommunicationTable[i].SendRank][TypeHydro]
             +sizeof(StructPstar) *NContainerSend[CommunicationTable[i].SendRank][TypeStar]
             +sizeof(StructPsink) *NContainerSend[CommunicationTable[i].SendRank][TypeSink];
        NExportThisTime[i][0] = CommunicationTable[i].SendSize;
        NExportThisTime[i][1] = NContainerSend[CommunicationTable[i].SendRank][0]; // Body or DM
        NExportThisTime[i][TypeHydro+1] = NContainerSend[CommunicationTable[i].SendRank][TypeHydro];
        NExportThisTime[i][TypeStar+1]  = NContainerSend[CommunicationTable[i].SendRank][TypeStar];
        NExportThisTime[i][TypeSink+1]  = NContainerSend[CommunicationTable[i].SendRank][TypeSink];

        MPI_Isend(NExportThisTime[i],NTypes+1,MPI_INT,
                CommunicationTable[i].SendRank,TAG_DECOMPOSITION_SENDSIZE+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);

        CheckSizeofBufferExportSendIndex(NExportThisTime[i][0],1,i);
        DecompositionExportSend[i] = BufferExportSend[i];
    }

    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportAll += NImportThisTime[i][0];
    }

    void *DecompositionExportRecv;
    CheckSizeofBufferExportRecv(NImportAll,1);
    DecompositionExportRecv = BufferExportRecv;

    int counter_recv = 0;
    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NImportThisTime[i][0]>0){
            MPI_Irecv(DecompositionExportRecv+NImportAll,
                NImportThisTime[i][0],MPI_BYTE,CommunicationTable[i].RecvRank,
                    TAG_DECOMPOSITION_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            NImportAll += NImportThisTime[i][0];
            counter_recv ++;
        }
    }

    // int NImport = 0;
    int counter_send = 0;
    for(int i=0;i<NProcs-1;i++){
        // Packing
        DataPackingNew(CommunicationTable[i].SendRank,
                NContainerSend[CommunicationTable[i].SendRank],
                DecompositionExportSend[i]);

        if(NExportThisTime[i][0]>0){
            MPI_Isend(DecompositionExportSend[i],
                NExportThisTime[i][0],MPI_BYTE,CommunicationTable[i].SendRank,
                    TAG_DECOMPOSITION_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
    }

    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);


    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NImportThisTime[i][0]>0){
            DataUnPackingNew(NImportThisTime[i]+1,DecompositionExportRecv+NImportAll);
            NImportAll += NImportThisTime[i][0];
        }
    }

    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);

    ReConnectPointers();

    MPI_Allreduce(MPI_IN_PLACE,&NotMyRank,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(NotMyRank > 0){
        //dprintlmpi(NotMyRank);
        DecompositionIterationCount ++;
        PrepareDomainKey();
        SetDecompositionKey();
        ParticlesExchangeLimitedNew();
    }

    return;
}

static void DataPacking(const int SendRank, int NContainerSend[]){

    StructPbodyptr  PbSend = DecompositionSend;
    StructPhydroptr PhSend = (StructPhydroptr)(PbSend+NContainerSend[0]);
    StructPstarptr  PsSend = (StructPstarptr)(PhSend+NContainerSend[TypeHydro]);
    StructPsinkptr  PskSend = (StructPsinkptr)(PsSend+NContainerSend[TypeStar]);

    int nbody = 0;
    int nhydro = 0;
    int nstars = 0;
    int nsinks = 0;

    int counter = 0;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(DomainKey[counter] == SendRank){
            PbSend[nbody] = *Pb;
            Pb->Use = false;
            nbody ++;
            if(Pb->Type == TypeHydro){
                PhSend[nhydro] = *((StructPhydroptr)(Pb->Baryon));
                ((StructPhydroptr)(Pb->Baryon))->Use = false;
                nhydro ++;
                Pall.Nhydro --;
            }else if(Pb->Type == TypeStar){
                PsSend[nstars] = *((StructPstarptr)(Pb->Baryon));
                ((StructPstarptr)(Pb->Baryon))->Use = false;
                nstars ++;
                Pall.Nstars --;
            }else if(Pb->Type == TypeSink){
                PskSend[nsinks] = *((StructPsinkptr)(Pb->Baryon));
                ((StructPsinkptr)(Pb->Baryon))->Use = false;
                nsinks ++;
                Pall.Nsink --;
            }else if(Pb->Type == TypeDM){
                Pall.NDM --;
            }
        }
        if(nbody == NContainerSend[0]){
            Pall.Ntotal -= nbody;
            return;
        }
        counter ++;
    }
    
}


static void DataUnPacking(int NContainerRecv[]){

    StructPbodyptr  PbRecv = DecompositionRecv;
    StructPhydroptr PhRecv = (StructPhydroptr)(PbRecv+NContainerRecv[0]);
    StructPstarptr  PsRecv = (StructPstarptr)(PhRecv+NContainerRecv[TypeHydro]);
    StructPsinkptr  PskRecv = (StructPsinkptr)(PsRecv+NContainerRecv[TypeStar]);

    if(NContainerRecv[0] == 0)
        return;

    // count unuse
    int counter = 0;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next)
        counter += Pb->Use?0:1;
        //counter += !(Pb->Use);
    if(NContainerRecv[0]-counter > 0){
        if(PbodySize == 0){
            GenerateStructPbody(NContainerRecv[0]-counter);
        }else{
            //AddStructPbody(NContainerRecv[0]-counter);
            StretchStructPbody(NContainerRecv[0]-counter);
        }
    }

    // Hydro array.
    counter = 0;
    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next)
        counter += Ph->Use?0:1;
        //counter += !(Ph->Use);

    if(NContainerRecv[TypeHydro]-counter > 0){
        if(PhydroSize == 0){
            GenerateStructPhydro(NContainerRecv[TypeHydro]-counter);
        }else{
            //AddStructPhydro(NContainerRecv[TypeHydro]-counter);
            StretchStructPhydro(NContainerRecv[TypeHydro]-counter);
        }
    }

    // Star array.
    counter = 0;
    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next)
        counter += Ps->Use?0:1;
        //counter += !(Ps->Use);
    if(NContainerRecv[TypeStar]-counter > 0){
        if(PstarSize == 0){
            GenerateStructPstar(NContainerRecv[TypeStar]-counter);
        }else{
            //AddStructPstar(NContainerRecv[TypeStar]-counter);
            StretchStructPstar(NContainerRecv[TypeStar]-counter);
        }
    }

    // Sink array.
    counter = 0;
    for(StructPsinkptr Psk = PsinkElements; Psk; Psk = Psk->Next)
        counter += Psk->Use?0:1;
        //counter += !(Psk->Use);
    if(NContainerRecv[TypeSink]-counter > 0){
        if(PsinkSize == 0){
            GenerateStructPsink(NContainerRecv[TypeSink]-counter);
        }else{
            //AddStructPstar(NContainerRecv[TypeStar]-counter);
            StretchStructPsink(NContainerRecv[TypeSink]-counter);
        }
    }


    int nbody = 0;
    int nhydro = 0;
    int nstars = 0;
    int nsinks = 0;
    StructPhydroptr Ph = PhydroElements;
    StructPstarptr Ps = PstarElements;
    StructPsinkptr Psk= PsinkElements;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == false){
            StructPbodyptr PbNext = Pb->Next;
            *Pb = PbRecv[nbody];
            Pb->Next = PbNext;
            Pb->Use = true;
            nbody ++;
            Pall.Ntotal ++;

            if(Pb->Type == TypeHydro){
                for(StructPhydroptr Pt=Ph; Pt; Pt = Pt->Next){
                    if(Pt->Use == false){
                        StructPhydroptr PhNext = Pt->Next;
                        *Pt = PhRecv[nhydro];
                        Pt->Next = PhNext;
                        Pt->Use = true;
                        Pb->Baryon = Pt;
                        Pt->Body = Pb;
                        Ph = Pt->Next;
                        break;
                    }
                }
                nhydro ++;
                Pall.Nhydro ++;
            }else if(Pb->Type == TypeStar){
                for(StructPstarptr Pt=Ps; Pt; Pt = Pt->Next){
                    if(Pt->Use == false){
                        StructPstarptr PsNext = Pt->Next;     
                        *Pt = PsRecv[nstars];
                        Pt->Next = PsNext;
                        Pt->Use = true;
                        Pb->Baryon = Pt;
                        Pt->Body = Pb;
                        Ps = Pt->Next;
                        break;
                    }
                }
                nstars ++;
                Pall.Nstars ++;
            }else if(Pb->Type == TypeSink){
                for(StructPsinkptr Pt=Psk; Pt; Pt = Pt->Next){
                    if(Pt->Use == false){
                        StructPsinkptr PskNext = Pt->Next;     
                        *Pt = PskRecv[nsinks];
                        Pt->Next = PskNext;
                        Pt->Use = true;
                        Pb->Baryon = Pt;
                        Pt->Body = Pb;
                        Psk = Pt->Next;
                        break;
                    }
                }
                nsinks ++;
                Pall.Nsink ++;
            } else if (Pb->Type == TypeDM){
                Pall.NDM ++;
            }
        }
        if(nbody == NContainerRecv[0])
            break;
    }
    /*
    fprintf(stderr,"UnPacking [%02d] %d %d %d %d | %d %d %d %d\n",MPIGetMyID(),nbody,nhydro,nstars,nsinks,
            NContainerRecv[0],NContainerRecv[1],NContainerRecv[2],NContainerRecv[3]);

    fflush(NULL);
    if(nhydro != NContainerRecv[TypeHydro]){
        fprintf(stderr,"[%02d] %d %d -- \n",MPIGetMyID(),nhydro,NContainerRecv[TypeHydro]);
        fflush(NULL);
        assert(nhydro == NContainerRecv[TypeHydro]);

    }
    */
    assert(nbody == NContainerRecv[0]);
    assert(nhydro == NContainerRecv[TypeHydro]);
    assert(nstars == NContainerRecv[TypeStar]);
    assert(nsinks == NContainerRecv[TypeSink]);
    return ;
}

static void DomainDecompResult(char *name){

    FILE *fp;

    FileOpen(fp,name,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %d\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],Pbody[i]->Use);
    }
    fclose(fp);
    return;
}

static void DomainOutPutKey(void){

    char name[256];
    FILE *fp;

    sprintf(name,"DomainKey.%02d",MPIGetMyID());
    FileOpen(fp,name,"w");
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next)
        fprintf(fp,"%ld %g %g %g\n",Pb->GlobalID,Pb->Pos[0],Pb->Pos[1],Pb->Pos[2]);
    fclose(fp);
    return;
}
