#include "config.h"
#include "NeighborSearch.h"
#include "IMFParameters.h"
#include "IMF.h"

void AddFeedbackEnergy(const int NLocalTarget);
static void FeedbackAddEnergyEngine(double Pos[restrict], 
        const double Kerneli, const double InvRho, const double Qheat);

static double SNIINumber;
static double SNIIDensity;
static double SNIIEnergy;

void InitializeSNII(void){

    InitializeIMF();
    SNIINumber = IMFSNIINumberPerMass();
    SNIIDensity = IMFSNIINumberPerMassLimitedRange(MSNIIMax,MSNIIMin);
    //SNIIEnergy = SNIIEnergyEfficiency*SNIIEnergyPerNumber*GetUnitEnergy();
    SNIIEnergy = SNIIEnergyEfficiency*SNIIEnergyPerNumber*SNIINumber;

    fprintf(stderr,"SNII number/Msun = %g, SNII density = %g, SNII energy/Msun %g\n",
            SNIINumber,SNIIDensity,SNIIEnergy);
    //fprintf(stderr,"%g,%g,%g\n",Pall.SNIIEnergy*GetUnitEnergy(),
            //Pall.SNIIEnergy,GetUnitEnergy());

    return;
}

#define CandidatesAllActives                    (ON)
#define CandidatesSlabActives                   (OFF)
#define CandidatesHighDensity                   (OFF)
#define CandidatesHighDensityLowTemperature     (OFF)

#define ThresholdDensity        (10000.) // nH/cc
#define ThresholdTemperature    (1000) // Kelvin
#define SlabR                   (51)   // pc 
#define SlabZ                   (4)   // pc 

static int CountCandidates(void){

    int NCandidates = 0;

#if (CandidatesAllActives)
    for(int i=0;i<Pall.Nhydro;i++){ // all actives.
        if(PhydroActive(i)){
            NCandidates ++;
        }
    }
#elif (CandidatesSlabActives)
    for(int i=0;i<Pall.Nhydro;i++){ // only high density particles.
        if(PhydroActive(i)){
            if( (fabs(PhydroPosP(i)[2])<SlabZ)&&(SQ(PhydroPosP(i)[0])+SQ(PhydroPosP(i)[1])< SQ(SlabR)) ){
                NCandidates ++;
            }
        }
    }
#elif (CandidatesHighDensity)
    for(int i=0;i<Pall.Nhydro;i++){ // only high density particles.
        if(PhydroActive(i)){
            if(Pall.ConvertNumberDenstityToCGS*Phydro[i]->Rho > ThresholdDensity){
                NCandidates ++;
            }
        }
    }
#elif (CandidatesHighDensityLowTemperature)
    for(int i=0;i<Pall.Nhydro;i++){ // only high density / low temperature particles.
        if(PhydroActive(i)){
            if( (Pall.ConvertNumberDenstityToCGS*Phydro[i]->Rho > ThresholdDensity)
                    &&(Pall.ConvertUtoT*Phydro[i]->U < ThresholdTemperature)){
                NCandidates ++;
            }
        }
    }
#else
    fprintf(stderr,"Flag error:%s:line %d\n",__FUNCTION__,__LINE__);
    MPI_Abort(MPI_COMM_WORLD,ConstantSNeCandidateError);
#endif


    return NCandidates;
}

static double ConstSNIIResidualMass;
static double ConstSNIINumber;
static double ConstSNIIStandardMass;
/*
*/

void ConstantSNII(const double SFRate){

    double ConvertSFRinSimulationUnit =
        (Pall.UnitTime/YEAR_CGS)/
        (Pall.UnitMass/MSUN_CGS);

    ConstSNIIResidualMass = 
        SFRate*ConvertSFRinSimulationUnit*Pall.TCurrent
        -ConstSNIINumber*ConstSNIIStandardMass;

    if(ConstSNIIResidualMass > ConstSNIIStandardMass){
        int NSNII = 0;
        while(ConstSNIIResidualMass > (NSNII+1)*ConstSNIIStandardMass)
            NSNII ++;

        int NProcs = MPIGetNumProcs();
        if(NSNII > 0){
            // search active particles with some feedback conditions.
            int NCandidates = CountCandidates();
            int NGlobalCandidates;
            MPI_Allreduce(&NCandidates,&NGlobalCandidates,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            if(NGlobalCandidates >= NSNII){
                int GlobalCandidates[NProcs];
                MPI_Gather(&NCandidates,1,MPI_INT,GlobalCandidates,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
                int Distribution[NProcs];
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    int Tail[NProcs];
                    Distribution[0] = 0;
                    Tail[0] = GlobalCandidates[0];
                    for(int i=1;i<NProcs;i++){
                        Distribution[i] = 0;
                        Tail[i] = Tail[i-1] + GlobalCandidates[i];
                    }

                    for(int i=0;i<NGlobalCandidates;i++){
                        int irand = gsl_rng_uniform_int(RandomGenerator,NGlobalCandidates);
                        for(int i=0;i<NProcs;i++){
                            if((irand < Tail[i])&&(Distribution[i] < GlobalCandidates[i]))
                                Distribution[i] ++;
                        }                        
                    }
                }
                int NLocalCandidates;
                MPI_Scatter(Distribution,1,MPI_INT,&NLocalCandidates,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
                fprintf(stderr," [%02d] Number of star forming regions : %d\n",MPIGetMyID(),NLocalCandidates);
                AddFeedbackEnergy(NLocalCandidates);

                ConstSNIINumber += NSNII;

                fprintf(stderr," Current Star formation rate %g Msun/yr | NSNII = %d\n",
                        ConstSNIINumber*(ConstSNIIStandardMass*(Pall.UnitMass/MSUN_CGS))/
                            (Pall.TCurrent*Pall.UnitTime/YEAR_CGS),NSNII);
            }
        }
    }

    return ;
}

#if 0
inline bool OverlapDomainSphereDensity(double Pos[restrict], const double h, const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}


int CheckExportFlagsFeedback(const int Index){

    int NodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;
    unsigned int BitMask = 1<<Index;

    StructNodeptr rnptr = RootHydro.Node;
    StructNodeptr cnptr = rnptr->Children;
    while(cnptr != rnptr){
        if(cnptr->ActiveNumber == 0){
            cnptr = cnptr->Next;
        } else if( !OverlapDomainSphereDensity(cnptr->Pos,cnptr->KernelMax,NodeID) ){ // update
            cnptr = cnptr->Next;
        } else if(cnptr->Children != NULL){
            cnptr = cnptr->Children;
        } else {
            int leaf = TREELEAVES(cnptr);
            for(int k=0;k<cnptr->Number;k++){
                if(PhydroActive(leaf)){
                    if(OverlapDomainSphereDensity(PhydroPosP(leaf),2.0*Phydro[leaf]->Kernel,NodeID)){ // update
                        Phydro[leaf]->ExportFlag |= BitMask;
                        NExport ++;
                    }
                }
                leaf = Phydro[leaf]->NextLeaf;
            }
            cnptr = cnptr->Next;
        }
    }
    return NExport;
}
#endif

void AddFeedbackEnergy(const int NLocalTarget){

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    int *ActiveIndexList;
    ActiveIndexList = malloc(sizeof(int)*NLocalTarget+1);

    for(int i=0;i<NLocalTarget;i++)
        ActiveIndexList[i] = 0;

    int NActives = 0;
    while(NActives != NLocalTarget){
        int leaf = gsl_rng_uniform_int(RandomGenerator,Pall.Nhydro);
        ActiveIndexList[NActives] = leaf;
        int count;
        for(count=0;count<NActives;count++)
            if(ActiveIndexList[count] == ActiveIndexList[NActives])
                break;
        if(count == NActives)
            NActives ++;
    }

    //allocation.
    double *Qheat;
    Qheat = malloc(sizeof(double)*NLocalTarget); 
    double InputEnergyPerMass = SNIIEnergyEfficiency*SNIIEnergyPerNumber*SNIINumber*GetUnitEnergy();
    for(int i=0;i<NLocalTarget;i++){ // Calc feedback energy.
        int leaf = ActiveIndexList[i];
        Qheat[i] = PhydroMass(leaf)*GetUnitMassSolarMass()*InputEnergyPerMass;
    }

struct StructFeedbackExport{
    double    Kernel;  // Kernel size.
    double    InvRho;     // density.
    double    Pos[3];  // Position.
    double    Qheat;   // Energy.
} *FeedbackExportSend[NProcs],*FeedbackExportRecv;
    //int *FeedbackImportSend,*FeedbackImportRecv;

    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    int NExportMax = 0;
    for(int i=0;i<NProcs-1;i++){

        // export flag is necessary!
        NExportThisTime[i] = 0;
        int NodeID = CommunicationTable[i].SendRank;
        for(int j=0;j<NLocalTarget;j++){
            int leaf = ActiveIndexList[j];
            double Dist2 = 0.e0;
            for(int k=0;k<3;k++){
                if(PhydroPosP(leaf)[k] < EdgesForHydro[NodeID].PosMin[k])
                    Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-PhydroPosP(leaf)[k]);
                if(PhydroPosP(leaf)[k] > EdgesForHydro[NodeID].PosMax[k])
                    Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-PhydroPosP(leaf)[k]);
            }
            if(Dist2 < SQ(2.0*Phydro[leaf]->Kernel))
                NExportThisTime[i] ++;
        }
        FeedbackExportSend[i] = malloc(sizeof(struct StructFeedbackExport)*NExportThisTime[i]+1);
        
        int NExport = 0;
        for(int k=0;k<NLocalTarget;k++){
            int leaf = ActiveIndexList[k];
            FeedbackExportSend[i][NExport].InvRho = 1.e0/Phydro[leaf]->Rho;
            FeedbackExportSend[i][NExport].Kernel = Phydro[leaf]->Kernel;
            FeedbackExportSend[i][NExport].Pos[0] = PhydroPosP(leaf)[0];
            FeedbackExportSend[i][NExport].Pos[1] = PhydroPosP(leaf)[1];
            FeedbackExportSend[i][NExport].Pos[2] = PhydroPosP(leaf)[2];
            FeedbackExportSend[i][NExport].Qheat = Qheat[k];
            NExport ++;
        }
        NExportMax = MAX(NExport,NExportMax); 
    }

    for(int i=0;i<NProcs-1;i++)
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_FEEDBACK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_FEEDBACK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);

    int NImportAllocated = 0;
    for(int i=0;i<NProcs-1;i++)
        NImportAllocated += NImportThisTime[i];

    FeedbackExportRecv = malloc(sizeof(struct StructFeedbackExport)*NImportAllocated+1);

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(FeedbackExportSend[i],
            NExportThisTime[i]*sizeof(struct StructFeedbackExport),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_FEEDBACK_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(FeedbackExportRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructFeedbackExport),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FEEDBACK_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);


    // calculation for local
    for(int k=0;k<NLocalTarget;k++){
        int leaf = ActiveIndexList[k];
        FeedbackAddEnergyEngine(PhydroPosP(leaf),Phydro[leaf]->Kernel,1.e0/Phydro[leaf]->Rho,Qheat[k]);
    }


    int NImportAll = NImportAllocated;
    for(int k=0;k<NImportAll;k++){
        FeedbackAddEnergyEngine(FeedbackExportRecv[k].Pos,FeedbackExportRecv[k].Kernel,
                FeedbackExportRecv[k].InvRho,FeedbackExportRecv[k].Qheat);
    }

    free(ActiveIndexList);
    free(Qheat);

    for(int i=0;i<NProcs-1;i++)
        free(FeedbackExportSend[i]);
    free(FeedbackExportRecv);

    return;
}

#if 0
void DecisionAndAdditionFeedbackRadiusAndEnergy(const int NLocalTarget){

    int NProcs = MPIGetNumProcs();
    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
    MPI_Status  mpi_status;

    int *ActiveIndexList;
    ActiveIndexList = malloc(sizeof(int)*NLocalTarget+1);

    for(int i=0;i<NLocalTarget;i++)
        ActiveIndexList[i] = 0;

    int NActives = 0;
    while(NActives != NLocalTarget){
        int leaf = gsl_rng_uniform_int(RandomGenerator,Pall.Nhydro);
        ActiveIndexList[NActives] = leaf;
        int count;
        for(count=0;count<NActives;count++)
            if(ActiveIndexList[count] == ActiveIndexList[NActives])
                break;
        if(count == NActives)
            NActives ++;
    }

    int NLocalActiveLeaves = NActives;
    int NActiveLeaves = NActives; // current active particles.
    bool *LocalActiveFlags;
    LocalActiveFlags = malloc(sizeof(bool)*NLocalActiveLeaves+1);
    for(int i=0;i<NActives;i++)
        LocalActiveFlags[i] = ON;


    // counter for export and import. 
    int BitMask = 0x01; 
    int NExportSizeAllocated[NProcs];
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];
    for(int i=0;i<NProcs;i++)
        NExportSizeAllocated[i] = 0;

struct StructFeedbackExport{
    double    Kernel;  // Kernel size.
    double    Pos[3];  // Position.
} *FeedbackExportSend[NProcs],*FeedbackExportRecv;

    int *FeedbackImportSend,*FeedbackImportRecv;

    int NExportMax = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportSizeAllocated[i] = NExportThisTime[i] = CheckExportFlagsFeedback(i);
        FeedbackExportSend[i] = 
            malloc(sizeof(struct StructFeedbackExport)*NExportSizeAllocated[i]+1);

        int NExport = 0;
        for(int k=0;k<NActives;k++){
            int index = ActiveIndexList[k];
            if(((Phydro[index]->ExportFlag)>>i)&BitMask){ 
                FeedbackExport[i][NExport].Kernel = Phydro[index]->Kernel;
                FeedbackExport[i][NExport].Pos[0] = PhydroPosP(index)[0];
                FeedbackExport[i][NExport].Pos[1] = PhydroPosP(index)[1];
                FeedbackExport[i][NExport].Pos[2] = PhydroPosP(index)[2];
                NExport ++;
            }
        }
        NExportMax = MAX(NExport,NExportMax);
    }
    FeedbackImportRecv = malloc(sizeof(int)*NExportMax+1);

    for(int i=0;i<NProcs-1;i++)
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_FEEDBACK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_FEEDBACK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);

    int NImportAllocated = 0;
    for(int i=0;i<NProcs-1;i++)
        NImportAllocated += NImportThisTime[i];

    FeedbackExportRecv = malloc(sizeof(struct StructFeedbackExport)*(NImportAllocated+1));
    FeedbackImportSend = malloc(sizeof(int)*(NImportAllocated+1));

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(FeedbackExportSend[i],
            NExportThisTime[i]*sizeof(struct StructFeedbackExport),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_FEEDBACK_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(FeedbackExportRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructFeedbackExport),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FEEDBACK_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    
    // start iteration.
    int Niteration = 0;
    do{
        if(Niteration != 0){ 
            for(int i=0;i<NActives;i++){ 
                if(LocalActiveFlags[i]){
                    int index = ActiveIndexList[i];
                    Phydro[index]->Nlist = 0;
                    Phydro[index]->Mnbs = 
                    Phydro[index]->Rho = 
                    Phydro[index]->Div = 
                    Phydro[index]->Rot[0] = 
                    Phydro[index]->Rot[1] = 
                    Phydro[index]->Rot[2] = 
                    Phydro[index]->Omega = 0.e0;
                    Phydro[index]->ExportFlag = OFF;
                }
            }

            int NExportMaxThisTime = 0;
            for(int i=0;i<NProcs-1;i++){
                NExportThisTime[i] = CheckExportFlagsDensityBitOperation(i);
                if(NExportSizeAllocated[i] < NExportThisTime[i]){
                    NExportSizeAllocated[i] = NExportThisTime[i];
                    HydroDensityExportSend[i] = realloc(HydroDensityExportSend[i],
                            sizeof(struct StructHydroDensityExport)*NExportSizeAllocated[i]+1);
                }
                assert(NExportThisTime[i]<=NExportSizeAllocated[i]);

                int NExport = 0;
                for(int k=0;k<NActives;k++){
                    if(LocalActiveFlags[k]){
                        int index = ActiveIndexList[k];
                        if(((Phydro[index]->ExportFlag)>>i)&BitMask){ 
                            E[i][NExport].Kernel = Phydro[index]->Kernel;
                            HydroDensityExportSend[i][NExport].Pos[0] = PhydroPosP(index)[0];
                            HydroDensityExportSend[i][NExport].Pos[1] = PhydroPosP(index)[1];
                            HydroDensityExportSend[i][NExport].Pos[2] = PhydroPosP(index)[2];
                            HydroDensityExportSend[i][NExport].Vel[0] = Phydro[index]->VelP[0];
                            HydroDensityExportSend[i][NExport].Vel[1] = Phydro[index]->VelP[1];
                            HydroDensityExportSend[i][NExport].Vel[2] = Phydro[index]->VelP[2];
                            NExport ++;
                        }
                    }
                }
                NExportThisTime[i] = NExport;

                NExportMaxThisTime = MAX(NExportMaxThisTime,NExport);
            }

            if(NExportMax<NExportMaxThisTime){
                NExportMax = NExportMaxThisTime;
                Delete(HydroDensityImportRecv);
                HydroDensityImportRecv = malloc(sizeof(struct StructHydroDensityImport)*NExportMax+1);
            }

            for(int i=0;i<NProcs-1;i++)
                MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,
                        CommunicationTable[i].SendRank,TAG_SPH_DENSITY_PRECOMM,
                    NImportThisTime+i,1,MPI_INT,
                        CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_PRECOMM,
                            MPI_COMM_WORLD,&mpi_status);

            NImport = 0;
            for(int i=0;i<NProcs-1;i++)
                NImport += NImportThisTime[i];

            if(NImportAllocated<NImport){
                NImportAllocated = NImport;
                Delete(HydroDensityImportSend);
                Delete(HydroDensityExportRecv);
                HydroDensityImportSend = malloc(sizeof(struct StructHydroDensityImport)*NImportAllocated+1);
                HydroDensityExportRecv = malloc(sizeof(struct StructHydroDensityExport)*NImportAllocated+1);
            }

            NImport = 0;
            for(int i=0;i<NProcs-1;i++){
                MPI_Isend(HydroDensityExportSend[i],
                    NExportThisTime[i]*sizeof(struct StructHydroDensityExport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_EXPORT,
                            MPI_COMM_WORLD,mpi_request_Export_Send+i);
                MPI_Irecv(HydroDensityExportRecv+NImport,
                    NImportThisTime[i]*sizeof(struct StructHydroDensityExport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_EXPORT,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+i);
                NImport += NImportThisTime[i];
            }
            double TimeComm = GetElapsedTime();
            MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
            TimingResults.HydroDensityCommThisStep += GetElapsedTime()-TimeComm;
            TimeComm = GetElapsedTime();
            MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
            TimingResults.HydroDensityCommThisStep += GetElapsedTime()-TimeComm;
        }

        for(int i=0;i<NActives;i++){ 
            if(LocalActiveFlags[i]){
                int index = ActiveIndexList[i];

                struct StructHydroDensityExport TempHydroDensityExport;
                TempHydroDensityExport.Kernel = Phydro[index]->Kernel;
                TempHydroDensityExport.Pos[0] = PhydroPosP(index)[0];
                TempHydroDensityExport.Pos[1] = PhydroPosP(index)[1];
                TempHydroDensityExport.Pos[2] = PhydroPosP(index)[2];
                TempHydroDensityExport.Vel[0] = Phydro[index]->VelP[0];
                TempHydroDensityExport.Vel[1] = Phydro[index]->VelP[1];
                TempHydroDensityExport.Vel[2] = Phydro[index]->VelP[2];

                struct StructHydroDensityImport TempHydroDensityImport =
                    ReturnStructureDensityDivRotOmegaEngine_PreFetch(
                        TempHydroDensityExport.Pos,TempHydroDensityExport.Vel,TempHydroDensityExport.Kernel);

                Phydro[index]->Rho += TempHydroDensityImport.Rho;
                Phydro[index]->Div += TempHydroDensityImport.Div;
                Phydro[index]->Rot[0] += TempHydroDensityImport.Rot[0];
                Phydro[index]->Rot[1] += TempHydroDensityImport.Rot[1];
                Phydro[index]->Rot[2] += TempHydroDensityImport.Rot[2];
                Phydro[index]->Omega += TempHydroDensityImport.Omega;
                Phydro[index]->Mnbs += TempHydroDensityImport.Mnbs;
                Phydro[index]->Nlist += TempHydroDensityImport.Nlist;
            }
        }

        int NImportAll = 0;
        for(int i=0;i<NProcs-1;i++)
            NImportAll += NImportThisTime[i];
        assert(NImportAll<=NImportAllocated);

        for(int i=0;i<NImportAll;i++)
            HydroDensityImportSend[i] = ReturnStructureDensityDivRotOmegaEngine_PreFetch(
                HydroDensityExportRecv[i].Pos,HydroDensityExportRecv[i].Vel,
                    HydroDensityExportRecv[i].Kernel);

        NImport = 0;
        for(int i=0;i<NProcs-1;i++){

            MPI_Sendrecv(HydroDensityImportSend+NImport,
                NImportThisTime[i]*sizeof(struct StructHydroDensityImport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_IMPORT,
                        HydroDensityImportRecv,
                NExportThisTime[i]*sizeof(struct StructHydroDensityImport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_IMPORT,
                        MPI_COMM_WORLD,&mpi_status);
            NImport += NImportThisTime[i];

            int NExport = 0;
            for(int k=0;k<NActives;k++){ 
                if(LocalActiveFlags[k]){
                    int index = ActiveIndexList[k];
                    if(((Phydro[index]->ExportFlag)>>i)&BitMask){ 
                        Phydro[index]->Rho += HydroDensityImportRecv[NExport].Rho;
                        Phydro[index]->Div += HydroDensityImportRecv[NExport].Div;
                        Phydro[index]->Rot[0] += HydroDensityImportRecv[NExport].Rot[0];
                        Phydro[index]->Rot[1] += HydroDensityImportRecv[NExport].Rot[1];
                        Phydro[index]->Rot[2] += HydroDensityImportRecv[NExport].Rot[2];
                        Phydro[index]->Omega += HydroDensityImportRecv[NExport].Omega;
                        Phydro[index]->Mnbs += HydroDensityImportRecv[NExport].Mnbs;
                        Phydro[index]->Nlist += HydroDensityImportRecv[NExport].Nlist;
                        NExport ++;
                    }
                }
            }
            assert(NExport == NExportThisTime[i]);
        }

        NLocalActiveLeaves = 0;
        for(int i=0;i<NActives;i++){
            if(LocalActiveFlags[i]){ 
                int index = ActiveIndexList[i];
                int Nlist = GetNeighboringNumber_Rho(index);
                if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                    LocalActiveFlags[i] = OFF;
                }else{
                    KernelUpdatei_Rho(index,Nlist); 
                    NLocalActiveLeaves ++;
                }
            }
        }

        NActiveLeaves = NLocalActiveLeaves;
        for(int i=0;i<NProcs-1;i++){
            int NActiveLeavesComm;
            MPI_Sendrecv(&NLocalActiveLeaves,1,MPI_INT,
                    CommunicationTable[i].SendRank,TAG_SPH_DENSITY_PRECOMM,
                &NActiveLeavesComm,1,MPI_INT,
                    CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_PRECOMM,
                        MPI_COMM_WORLD,&mpi_status);
            NActiveLeaves += NActiveLeavesComm;
        }

        Niteration ++;
        if(Niteration > 10)
            break;

    } while (0<NActiveLeaves);

    for(int i=0;i<NProcs-1;i++){
        Delete(HydroDensityExportSend[i]);
    }
    Delete(HydroDensityExportRecv);
    Delete(HydroDensityImportSend);
    Delete(HydroDensityImportRecv);


    HydroDensityEndProcessing(NActives,ActiveIndexList);


    Delete(LocalActiveFlags);
    Delete(ActiveIndexList);

    return;
}
#endif

//inline double KernelFeedBack(const double r, const double InvKerneli) __attribute__((always_inline));
//inline double KernelFeedBack(const double r, const double InvKerneli){
static inline double __attribute__((always_inline)) KernelFeedBack(const double r, const double InvKerneli){

    double u = r*InvKerneli;
#if (DIMENSION == 1)
    const static double coef1d = 2.0/3.0;
    double coef = coef1d*InvKerneli;
#elif (DIMENSION == 2)   
    const static double coef2d = 10.0/(7.0*M_PI);
    double coef = coef2d*SQ(InvKerneli);
#elif (DIMENSION == 3)
    double coef = IPI*CUBE(InvKerneli);
#endif
    if(u<1.e0){
        return (coef*(1.e0 - 1.5*SQ(u) + 0.75*CUBE(u)));
    } else if (u<2.e0){
        return (coef*(0.25*CUBE(2.e0-u)));
    } else {
        return 0.e0;
    }
}

static void FeedbackAddEnergyEngine(double Pos[restrict], const double Kerneli, 
        const double InvRho, const double Qheat){

    int Neighbors[MaxNeighborSize];

    double TimeComm = GetElapsedTime();
    int Nlist = GetNeighborsLimited(Pos,2.e0*Kerneli,Neighbors);
    TimingResults.FeedbackNeighborSearchThisStep += GetElapsedTime()-TimeComm;

    if(Nlist > MaxNeighborSize){
        fprintf(stderr,"The neighbor list is over flow(Nlist = %d). function:%s,line:%d,file:%s\n",
                Nlist,__FUNCTION__,__LINE__,__FILE__);
        MPI_Abort(MPI_COMM_WORLD,NeighborListOverFlow);
    }

    double InvKerneli = 1.e0/Kerneli;
	for(int k=0;k<Nlist;k++){
        int leaf = Neighbors[k];
        double r = DISTANCE(Pos,PhydroPosP(leaf));
		double w = KernelFeedBack(r,InvKerneli);
        Phydro[leaf]->DQheat += PhydroMass(leaf)*InvRho*Qheat*w;
    }
    return ;
}
