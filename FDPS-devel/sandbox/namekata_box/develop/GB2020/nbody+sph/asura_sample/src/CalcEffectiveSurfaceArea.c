#include "config.h"
#include "CalcEffectiveSurfaceArea.h"
#include "StellarFeedback.h"
#include "NeighborSearch.h"
#include "SizeDetermination.h"
#include "KernelFunctions.h"
#include "PlantHydroTree.h"

#ifdef USE_MOMENTUM_FEEDBACK //{

int NExplosionLog;
//static int NActives = 0;
static int *ActiveIndexList = NULL;
static int ActiveIndexListAllocated = 0;
static int CurrentActiveIndexListSize = 0;
static bool *ExportFlagForGatherScatter = NULL;
static int *ExportFlagForGatherScatterIndexList = NULL;


struct StructEffSAExport{
    int     Leaf;
    double  Kernel;  // Kernel size.
    double  Pos[3];  // Position.
    double  NumberDensity;
    // double  SA;
    // double  WeightCorrection[EffSAVecSize];    // Corrector.
    int     DomainID;
    long int GlobalID;
};

static int NEffSAExportRecv = 0;
static int EffSAExportRecvAllocated = 0;
static struct StructEffSAExport *EffSAExportRecv = NULL;

struct StructEffSAResult{
    int     Leaf;
    int     Nlist;
    double  SA;              // Density.
    double  WeightCorrection[EffSAVecSize];    // Corrector.
    double  CheckWeight;
};

struct StructEffSAInput{
    double Pos[3];
    double Kernel;
    double NumberDensity;
    double SA;
    double WeightCorrection[EffSAVecSize];    // Corrector.
    int    Index;
};


static struct StructEdges *EdgesForEffSA = NULL;
static struct StructEdges *EdgesForEffSAWithExtent = NULL;
static bool FirstCall_CalcActiveDomainEdgesForESA = true;

static void CalcActiveDomainEdgesForESA(struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){ 

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    if(FirstCall_CalcActiveDomainEdgesForESA){
        EdgesForEffSA = malloc(sizeof(struct StructEdges)*NProcs);
        EdgesForEffSAWithExtent = malloc(sizeof(struct StructEdges)*NProcs);
        FirstCall_CalcActiveDomainEdgesForESA = false;
    }

    int NActives = CurrentActiveIndexListSize;

    struct StructEdges TempEdges = {0.e0};
    struct StructEdges TempEdgesWithExtent = {0.e0};

    double max[3] = {0.e0,0.e0,0.e0};
    double min[3] = {0.e0,0.e0,0.e0};
    double maxWithExtent[3] = {0.e0,0.e0,0.e0};
    double minWithExtent[3] = {0.e0,0.e0,0.e0};

    if(NActives > 0){
#if 1
        int head = ActiveIndexList[0];
        for(int k=0;k<3;k++){
            max[k] = ActiveStellarFeedbackParticle[head].Pos[k];
            min[k] = ActiveStellarFeedbackParticle[head].Pos[k];
            maxWithExtent[k] = ActiveStellarFeedbackParticle[head].Pos[k]+2.0*ActiveStellarFeedbackParticle[head].Radius;
            minWithExtent[k] = ActiveStellarFeedbackParticle[head].Pos[k]-2.0*ActiveStellarFeedbackParticle[head].Radius;
        }
        for(int i=1;i<NActives;i++){
            int leaf = ActiveIndexList[i];
            for(int k=0;k<3;k++){
                max[k] = fmax(ActiveStellarFeedbackParticle[leaf].Pos[k],max[k]);
                min[k] = fmin(ActiveStellarFeedbackParticle[leaf].Pos[k],min[k]);
                maxWithExtent[k] = fmax(ActiveStellarFeedbackParticle[leaf].Pos[k]+2.0*ActiveStellarFeedbackParticle[leaf].Radius,maxWithExtent[k]);
                minWithExtent[k] = fmin(ActiveStellarFeedbackParticle[leaf].Pos[k]-2.0*ActiveStellarFeedbackParticle[leaf].Radius,minWithExtent[k]);
            }
        }
#else
        for(int k=0;k<3;k++){
            max[k] = ActiveStellarFeedbackParticle[0].Pos[k];
            min[k] = ActiveStellarFeedbackParticle[0].Pos[k];
            maxWithExtent[k] = ActiveStellarFeedbackParticle[0].Pos[k]+2.0*ActiveStellarFeedbackParticle[0].Radius;
            minWithExtent[k] = ActiveStellarFeedbackParticle[0].Pos[k]-2.0*ActiveStellarFeedbackParticle[0].Radius;
        }

        for(int i=1;i<NExplosionLog;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(ActiveStellarFeedbackParticle[i].Pos[k],max[k]);
                min[k] = fmin(ActiveStellarFeedbackParticle[i].Pos[k],min[k]);
                maxWithExtent[k] = fmax(ActiveStellarFeedbackParticle[i].Pos[k]+2.0*ActiveStellarFeedbackParticle[i].Radius,maxWithExtent[k]);
                minWithExtent[k] = fmin(ActiveStellarFeedbackParticle[i].Pos[k]-2.0*ActiveStellarFeedbackParticle[i].Radius,minWithExtent[k]);
            }
        }

#endif

        for(int k=0;k<3;k++){
            TempEdges.PosMax[k] = max[k];
            TempEdges.PosMin[k] = min[k];
            TempEdgesWithExtent.PosMax[k] = maxWithExtent[k];
            TempEdgesWithExtent.PosMin[k] = minWithExtent[k];
        }
    } else {
        for(int k=0;k<3;k++){
            TempEdges.PosMax[k] = TempEdges.PosMin[k] = HUGE_VAL;
            TempEdgesWithExtent.PosMax[k] = TempEdgesWithExtent.PosMin[k] = HUGE_VAL;
        }
    }

    memset(EdgesForEffSA,0,sizeof(struct StructEdges)*NProcs);
    memset(EdgesForEffSAWithExtent,0,sizeof(struct StructEdges)*NProcs);
    EdgesForEffSA[MyID] = TempEdges;
    EdgesForEffSAWithExtent[MyID] = TempEdgesWithExtent;

    // Exchange 
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];


    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(EdgesForEffSA+CommunicationTable[i].RecvRank,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].RecvRank,TAG_EFFSA_EXPORT,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);

        MPI_Isend(&TempEdges,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].SendRank,TAG_EFFSA_EXPORT,
                MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);
    }

    MPI_Status mpi_status_Export_Recv_WithExtent[NProcs-1];
    MPI_Request mpi_request_Export_Recv_WithExtent[NProcs-1];
    MPI_Status mpi_status_Export_Send_WithExtent[NProcs-1];
    MPI_Request mpi_request_Export_Send_WithExtent[NProcs-1];

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(EdgesForEffSAWithExtent+CommunicationTable[i].RecvRank,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].RecvRank,TAG_EFFSA_EXPORT+1,
                MPI_COMM_WORLD,mpi_request_Export_Recv_WithExtent+i);
        MPI_Test(mpi_request_Export_Recv_WithExtent+i,&RecvFlag,MPI_STATUS_IGNORE);

        MPI_Isend(&TempEdgesWithExtent,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].SendRank,TAG_EFFSA_EXPORT+1,
                MPI_COMM_WORLD,mpi_request_Export_Send_WithExtent+i);
        MPI_Test(mpi_request_Export_Send_WithExtent+i,&SendFlag,MPI_STATUS_IGNORE);
    }

    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

    MPI_Waitall(NProcs-1,mpi_request_Export_Send_WithExtent,mpi_status_Export_Send_WithExtent);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv_WithExtent,mpi_status_Export_Recv_WithExtent);

#if 0
    if(MPIGetMyID() == MPI_ROOT_RANK){

        for(int i=0;i<MPIGetNumProcs();i++){
            fprintf(stderr,"[%02d][%02d] e %g %g %g %g %g %g\n",MPIGetMyID(),i,
                    EdgesForEffSA[i].PosMin[0],
                    EdgesForEffSA[i].PosMax[0],
                    EdgesForEffSA[i].PosMin[1],
                    EdgesForEffSA[i].PosMax[1],
                    EdgesForEffSA[i].PosMin[2],
                    EdgesForEffSA[i].PosMax[2]);
        }
        for(int i=0;i<MPIGetNumProcs();i++){
            fprintf(stderr,"[%02d][%02d] ee %g %g %g %g %g %g\n",MPIGetMyID(),i,
                    EdgesForEffSAWithExtent[i].PosMin[0],
                    EdgesForEffSAWithExtent[i].PosMax[0],
                    EdgesForEffSAWithExtent[i].PosMin[1],
                    EdgesForEffSAWithExtent[i].PosMax[1],
                    EdgesForEffSAWithExtent[i].PosMin[2],
                    EdgesForEffSAWithExtent[i].PosMax[2]);
        }


        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK+1){

        for(int i=0;i<MPIGetNumProcs();i++){
            fprintf(stderr,"[%02d][%02d] e %g %g %g %g %g %g\n",MPIGetMyID(),i,
                    EdgesForEffSA[i].PosMin[0],
                    EdgesForEffSA[i].PosMax[0],
                    EdgesForEffSA[i].PosMin[1],
                    EdgesForEffSA[i].PosMax[1],
                    EdgesForEffSA[i].PosMin[2],
                    EdgesForEffSA[i].PosMax[2]);
        }
        for(int i=0;i<MPIGetNumProcs();i++){
            fprintf(stderr,"[%02d][%02d] ee %g %g %g %g %g %g\n",MPIGetMyID(),i,
                    EdgesForEffSAWithExtent[i].PosMin[0],
                    EdgesForEffSAWithExtent[i].PosMax[0],
                    EdgesForEffSAWithExtent[i].PosMin[1],
                    EdgesForEffSAWithExtent[i].PosMax[1],
                    EdgesForEffSAWithExtent[i].PosMin[2],
                    EdgesForEffSAWithExtent[i].PosMax[2]);
        }


        fflush(NULL);
    }

#endif


    return ;
}

static struct StructEdges *EdgesForHydroEffSA = NULL;
static struct StructEdges *EdgesForHydroEffSAWithExtent = NULL;
static bool FirstCall_CalcHydroDomainEdgesForEffSA = true;

static void CalcHydroDomainEdgesForEffSA(void){ 

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    if(FirstCall_CalcHydroDomainEdgesForEffSA){
        EdgesForHydroEffSA = malloc(sizeof(struct StructEdges)*NProcs);
        EdgesForHydroEffSAWithExtent = malloc(sizeof(struct StructEdges)*NProcs);
        FirstCall_CalcHydroDomainEdgesForEffSA = false;
    }

    struct StructEdges TempEdges = {0.e0};
    struct StructEdges TempEdgesWithExtent = {0.e0};

    double max[3] = {0.e0},min[3] = {0.e0};
    double maxWithExtent[3] = {0.e0},minWithExtent[3] = {0.e0};
    if(Pall.Nhydro > 0){
#if 0
        for(int k=0;k<3;k++){
            max[k] = Phydro[0]->PosP[k];
            min[k] = Phydro[0]->PosP[k];
            maxWithExtent[k] = Phydro[0]->PosP[k]+2.0*Phydro[0]->KernelPred;
            minWithExtent[k] = Phydro[0]->PosP[k]-2.0*Phydro[0]->KernelPred;
        }
        for(int i=1;i<Pall.Nhydro;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(Phydro[i]->PosP[k],max[k]);
                min[k] = fmin(Phydro[i]->PosP[k],min[k]);
                maxWithExtent[k] = fmax(Phydro[i]->PosP[k]+2.0*Phydro[i]->KernelPred,maxWithExtent[k]);
                minWithExtent[k] = fmin(Phydro[i]->PosP[k]-2.0*Phydro[i]->KernelPred,minWithExtent[k]);
            }
        }

        for(int k=0;k<3;k++){
            TempEdges.PosMax[k] = max[k];
            TempEdges.PosMin[k] = min[k];
            TempEdgesWithExtent.PosMax[k] = maxWithExtent[k];
            TempEdgesWithExtent.PosMin[k] = minWithExtent[k];
        }
#else
        for(int k=0;k<3;k++){
            max[k] = NBCache[0].Pos[k];
            min[k] = NBCache[0].Pos[k];
            maxWithExtent[k] = NBCache[0].Pos[k]+2.0*NBCache[0].Kernel;
            minWithExtent[k] = NBCache[0].Pos[k]-2.0*NBCache[0].Kernel;
        }

        int NumberofLeaves = HydroNode[0].NumberofLeaves;
        for(int i=1;i<NumberofLeaves;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(NBCache[i].Pos[k],max[k]);
                min[k] = fmin(NBCache[i].Pos[k],min[k]);
                maxWithExtent[k] = fmax(NBCache[i].Pos[k]+2.0*NBCache[i].Kernel,maxWithExtent[k]);
                minWithExtent[k] = fmin(NBCache[i].Pos[k]-2.0*NBCache[i].Kernel,minWithExtent[k]);
            }
        }
#endif
    } else {
        for(int k=0;k<3;k++){
            TempEdges.PosMax[k] = TempEdges.PosMin[k] = HUGE_VAL;
            TempEdgesWithExtent.PosMax[k] = TempEdgesWithExtent.PosMin[k] = HUGE_VAL;
        }
    }

    memset(EdgesForHydroEffSA,0,sizeof(struct StructEdges)*NProcs);
    memset(EdgesForHydroEffSAWithExtent,0,sizeof(struct StructEdges)*NProcs);
    EdgesForHydroEffSA[MyID] = TempEdges;
    EdgesForHydroEffSAWithExtent[MyID] = TempEdgesWithExtent;

    // Exchange 
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];


    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(EdgesForHydroEffSA+CommunicationTable[i].RecvRank,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].RecvRank,TAG_EFFSA_EXPORT,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);

        MPI_Isend(&TempEdges,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].SendRank,TAG_EFFSA_EXPORT,
                MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);
    }

    MPI_Status mpi_status_Export_Recv_WithExtent[NProcs-1];
    MPI_Request mpi_request_Export_Recv_WithExtent[NProcs-1];
    MPI_Status mpi_status_Export_Send_WithExtent[NProcs-1];
    MPI_Request mpi_request_Export_Send_WithExtent[NProcs-1];

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(EdgesForHydroEffSAWithExtent+CommunicationTable[i].RecvRank,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].RecvRank,TAG_EFFSA_EXPORT+1,
                MPI_COMM_WORLD,mpi_request_Export_Recv_WithExtent+i);
        MPI_Test(mpi_request_Export_Recv_WithExtent+i,&RecvFlag,MPI_STATUS_IGNORE);

        MPI_Isend(&TempEdgesWithExtent,sizeof(struct StructEdges),MPI_BYTE,
            CommunicationTable[i].SendRank,TAG_EFFSA_EXPORT+1,
                MPI_COMM_WORLD,mpi_request_Export_Send_WithExtent+i);
        MPI_Test(mpi_request_Export_Send_WithExtent+i,&SendFlag,MPI_STATUS_IGNORE);
    }

    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

    MPI_Waitall(NProcs-1,mpi_request_Export_Send_WithExtent,mpi_status_Export_Send_WithExtent);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv_WithExtent,mpi_status_Export_Recv_WithExtent);

    return ;
}

StructHydroRoot HydroRootForEffSA;
struct StructHydroNode *HydroNodeForEffSA = NULL;
static struct StructNBCacheForEffSA{
    double Pos[3];            // Position of this leaf.
    double Kernel;            // Kernel size of this leaf.
    double NumberDensity;     // Number density
    int    DomainID;
    int    Interact;
    int    Leaf;
    long int GlobalID;
} *NBCacheForEffSA = NULL, *CachedDataForEffSA = NULL;

//static int NumberofAllocatedLeaves = 0;
static int *sublist = NULL; 
static double DiagHydroForEffSA;

static bool FirstCall_InitializeRootForEffSA = true;
static void InitializeRootForEffSA(void){

    if(!FirstCall_InitializeRootForEffSA) return ;
    FirstCall_InitializeRootForEffSA = false;

    int Nload = FirstAllocationSize;
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)FirstAllocationSize)/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = FirstAllocationSize*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));

    /* HydroRoot */
    HydroRootForEffSA.NumberofAllocatedLeaves = Nload;
	HydroRootForEffSA.Leaves = malloc(sizeof(int)*Nload+1);
    HydroRootForEffSA.NumberofAllocatedNodes = NumberofAllocatedNodes;
    HydroRootForEffSA.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

	HydroRootForEffSA.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForNBS;
    HydroRootForEffSA.MaxLevel = TreeMaxNodeLevel;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        HydroRootForEffSA.WidthFactor[i] = pow(dw,(double)i);
    /* HydroRoot */

    /* HydroNode */
    HydroNodeForEffSA = malloc(sizeof(struct StructHydroNode)*HydroRootForEffSA.NumberofAllocatedNodes);

    struct StructHydroNode HydroNodeTemp;
    memset(&HydroNodeTemp,0,sizeof(struct StructHydroNode));
    HydroNodeTemp.Next = HydroNodeTemp.Parent = HydroNodeTemp.Children = 
    HydroNodeTemp.Sister = NONE;

    for(int i=0;i<HydroRootForEffSA.NumberofAllocatedNodes;i++)
        HydroNodeForEffSA[i] = HydroNodeTemp;
    /* HydroNode */

    /* NBCache */
    NBCacheForEffSA = malloc(sizeof(struct StructNBCacheForEffSA)*Nload);
    /* NBCache */

    sublist = malloc(sizeof(int)*Nload);
    CachedDataForEffSA = malloc(sizeof(struct StructNBCacheForEffSA)*Nload);

	return;
}


static void HydroTreePreprocessingForEffSA(const int NumberofLeaves, struct StructEffSAExport EffSAExportRecv[restrict]){

    if(NumberofLeaves>HydroRootForEffSA.NumberofAllocatedLeaves){
        HydroRootForEffSA.NumberofAllocatedLeaves = (int)(ForAngelsShare*NumberofLeaves);
        HydroRootForEffSA.Leaves = realloc(HydroRootForEffSA.Leaves,sizeof(int)*HydroRootForEffSA.NumberofAllocatedLeaves);
	    sublist = realloc(sublist,sizeof(int)*HydroRootForEffSA.NumberofAllocatedLeaves);
        NBCacheForEffSA = realloc(NBCacheForEffSA,sizeof(struct StructNBCacheForEffSA)*HydroRootForEffSA.NumberofAllocatedLeaves);
        CachedDataForEffSA = realloc(CachedDataForEffSA,sizeof(struct StructNBCacheForEffSA)*HydroRootForEffSA.NumberofAllocatedLeaves);
    }

    for(int i=0;i<NumberofLeaves;i++){
        CachedDataForEffSA[i].Pos[0] = EffSAExportRecv[i].Pos[0];
        CachedDataForEffSA[i].Pos[1] = EffSAExportRecv[i].Pos[1];
        CachedDataForEffSA[i].Pos[2] = EffSAExportRecv[i].Pos[2];
        CachedDataForEffSA[i].Kernel = EffSAExportRecv[i].Kernel;
        CachedDataForEffSA[i].NumberDensity = EffSAExportRecv[i].NumberDensity;
        CachedDataForEffSA[i].DomainID = EffSAExportRecv[i].DomainID;
        CachedDataForEffSA[i].Interact = false;
        CachedDataForEffSA[i].Leaf = i;
        CachedDataForEffSA[i].GlobalID = EffSAExportRecv[i].GlobalID;
    }

    return ;
}

static void MakeHydroRootForEffSA(const int NumberofLeaves, struct StructEffSAExport EffSAExportRecv[restrict]){

	double min[3] = {0.e0,0.e0,0.e0};
    double max[3] = {0.e0,0.e0,0.e0};
    const int RootNodeID = 0;

    if(NumberofLeaves>0){
        min[0] = EffSAExportRecv[0].Pos[0];
        min[1] = EffSAExportRecv[0].Pos[1];
        min[2] = EffSAExportRecv[0].Pos[2];

        max[0] = EffSAExportRecv[0].Pos[0];
        max[1] = EffSAExportRecv[0].Pos[1];
        max[2] = EffSAExportRecv[0].Pos[2];

        for(int i=1;i<NumberofLeaves;i++){
            min[0] = fmin(min[0],EffSAExportRecv[i].Pos[0]);
            min[1] = fmin(min[1],EffSAExportRecv[i].Pos[1]);
            min[2] = fmin(min[2],EffSAExportRecv[i].Pos[2]);

            max[0] = fmax(max[0],EffSAExportRecv[i].Pos[0]);
            max[1] = fmax(max[1],EffSAExportRecv[i].Pos[1]);
            max[2] = fmax(max[2],EffSAExportRecv[i].Pos[2]);
        }
        for(int k=0;k<3;k++){
            HydroRootForEffSA.PosMax[k] = max[k];
            HydroRootForEffSA.PosMin[k] = min[k];
        }
    } else {
        for(int k=0;k<3;k++){
            HydroRootForEffSA.PosMax[k] = HUGE_VAL;
            HydroRootForEffSA.PosMin[k] = HUGE_VAL;
        }
    }


    double WidthMax = 0.e0;
	for(int k=0;k<3;k++){
        HydroNodeForEffSA[RootNodeID].Pos[k] = 0.5*(max[k] + min[k]);
        WidthMax = fmax(WidthMax,max[k]-min[k]);
    }
    HydroRootForEffSA.Width = WidthMax;

    HydroNodeForEffSA[RootNodeID].Next = NONE;
    HydroNodeForEffSA[RootNodeID].Parent = NONE;
    HydroNodeForEffSA[RootNodeID].Sister = NONE;
    HydroNodeForEffSA[RootNodeID].Children = NONE;

    HydroNodeForEffSA[RootNodeID].NumberofLeaves = NumberofLeaves;

    HydroNodeForEffSA[RootNodeID].Level = 0;
    HydroNodeForEffSA[RootNodeID].Leaves = 0;

    for(int i=0;i<NumberofLeaves;i++)
        HydroRootForEffSA.Leaves[i] = i;

    HydroRootForEffSA.NumberofLeaves = NumberofLeaves;

	return ;
}

static inline bool __attribute__((always_inline)) HydroNodeSeparationCriterionForEffSA(const int CurrentNodeID, const int CriticalNumber){

	if( (HydroNodeForEffSA[CurrentNodeID].NumberofLeaves <= CriticalNumber) || HydroNodeForEffSA[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}

static void BuildHydroTreeForEffSA(void){

    int NumberofNodes = 0; 
    int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 

    int NumberofNodeCreationLimit = HydroRootForEffSA.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
	while(1){

		if(HydroNodeSeparationCriterionForEffSA(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(HydroNodeForEffSA[CurrentNodeID].Sister != NONE){
				CurrentNodeID = HydroNodeForEffSA[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(HydroNodeForEffSA[HydroNodeForEffSA[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = HydroNodeForEffSA[HydroNodeForEffSA[NextNodeID].Parent].Sister;
						break;
					} else if(HydroNodeForEffSA[NextNodeID].Parent == RootNodeID){
                        HydroRootForEffSA.CurrentMaxLevel = CurrentMaxLevel;

                        int NumberofLeaves = HydroNodeForEffSA[RootNodeID].NumberofLeaves;
                        for(int k=0;k<NumberofLeaves;k++){
                            int leaf =  HydroRootForEffSA.Leaves[k];
                            NBCacheForEffSA[k] = CachedDataForEffSA[leaf];
                        }
                        HydroRootForEffSA.NumberofNodes = NumberofNodes + 1;
						return;
					}
                    NextNodeID = HydroNodeForEffSA[NextNodeID].Parent;
				}
			}
			continue;
		}


		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

		for(int i=0;i<HydroNodeForEffSA[CurrentNodeID].NumberofLeaves;i++){
            int leaf = HydroRootForEffSA.Leaves[HydroNodeForEffSA[CurrentNodeID].Leaves+i];
			int subindex = 0;

            for(int k=0;k<3;k++)
                if(HydroNodeForEffSA[CurrentNodeID].Pos[k] <= CachedDataForEffSA[leaf].Pos[k])
                    subindex += 1 << k;

			if(subnumber[subindex] == 0){
				subhead[subindex] = subcurrent[subindex] = leaf;
			} else {
				sublist[subcurrent[subindex]] = leaf;
				subcurrent[subindex] = leaf;
			}
			subnumber[subindex] ++;
        }

#if 0
        if(MPIGetMyID()==46){ 
            fprintf(stderr,"[%d] / %d/ %d %d %d %d %d %d %d %d\n",
                    MPIGetMyID(),CurrentNodeID,
                    subnumber[0],subnumber[1],subnumber[2],subnumber[3],
                    subnumber[4],subnumber[5],subnumber[6],subnumber[7]);
            fflush(NULL);
        }
#endif


        ChildNodeID = CurrentNodeID;
		for(int i=0;i<TreeNsub;i++){
			if(subnumber[i] != 0){
				BackwardNodeID = ChildNodeID; 
                // make node
                NumberofNodes ++;
                ChildNodeID = NumberofNodes;
                if(1.1*NumberofNodes >= HydroRootForEffSA.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(ForAngelsShare*NumberofNodes);
                    HydroNodeForEffSA = realloc(HydroNodeForEffSA,sizeof(struct StructHydroNode)*NumberofAllocatedNodes);
                    HydroRootForEffSA.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                HydroNodeForEffSA[ChildNodeID].Next = NONE;
                HydroNodeForEffSA[ChildNodeID].Parent = NONE;
                HydroNodeForEffSA[ChildNodeID].Sister = NONE;
                HydroNodeForEffSA[ChildNodeID].Children = NONE;

                HydroNodeForEffSA[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    HydroNodeForEffSA[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    HydroNodeForEffSA[ChildNodeID].Leaves = HydroNodeForEffSA[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,HydroNodeForEffSA[CurrentNodeID].Level+1);
                } else {
                    HydroNodeForEffSA[BackwardNodeID].Sister = ChildNodeID;
                    HydroNodeForEffSA[ChildNodeID].Leaves = HydroNodeForEffSA[BackwardNodeID].Leaves + HydroNodeForEffSA[BackwardNodeID].NumberofLeaves;
                }

                HydroRootForEffSA.Leaves[HydroNodeForEffSA[ChildNodeID].Leaves] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    HydroRootForEffSA.Leaves[HydroNodeForEffSA[ChildNodeID].Leaves+k] = 
                        sublist[HydroRootForEffSA.Leaves[HydroNodeForEffSA[ChildNodeID].Leaves+k-1]];
                }
                HydroNodeForEffSA[ChildNodeID].NumberofLeaves = subnumber[i];
                HydroNodeForEffSA[ChildNodeID].Level = HydroNodeForEffSA[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    HydroNodeForEffSA[ChildNodeID].Pos[k] = HydroNodeForEffSA[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*HydroRootForEffSA.Width*HydroRootForEffSA.WidthFactor[HydroNodeForEffSA[CurrentNodeID].Level];
            }
		}
        CurrentNodeID = NextNodeID;
	}
}

static int NextHydroNodeForEffSA(const int NodeID){

    int CurrentNodeID = NodeID;

    if(HydroNodeForEffSA[CurrentNodeID].Sister != NONE){
        CurrentNodeID = HydroNodeForEffSA[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(HydroNodeForEffSA[HydroNodeForEffSA[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = HydroNodeForEffSA[HydroNodeForEffSA[NextNodeID].Parent].Sister;
                break;
            } else if(HydroNodeForEffSA[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = HydroNodeForEffSA[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}

static void HydroNodeDataImplantForEffSA(const int CurrentNodeID){

    double Width = DiagHydroForEffSA*HydroRootForEffSA.WidthFactor[HydroNodeForEffSA[CurrentNodeID].Level];
    double DistanceMax = 0.e0;
    double KernelMax = 0.e0;

    if(HydroNodeForEffSA[CurrentNodeID].Children == NONE){
        int NumberofLeaves = HydroNodeForEffSA[CurrentNodeID].NumberofLeaves;
        int header = HydroNodeForEffSA[CurrentNodeID].Leaves;
        for(int k=0;k<NumberofLeaves;k++){
            int leaf = header+k;
            double Distance = DISTANCE(HydroNodeForEffSA[CurrentNodeID].Pos,NBCacheForEffSA[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            KernelMax = fmax(KernelMax,2.e0*NBCacheForEffSA[leaf].Kernel+Distance);
        }
        HydroNodeForEffSA[CurrentNodeID].KernelMax = KernelMax;
        HydroNodeForEffSA[CurrentNodeID].DistanceMax = DistanceMax;
        HydroNodeForEffSA[CurrentNodeID].NumberofActiveLeaves = NumberofLeaves;
    } else {
        int ChildNodeID = HydroNodeForEffSA[CurrentNodeID].Children;
        int NActives = 0;
        do{
            HydroNodeDataImplantForEffSA(ChildNodeID);

            KernelMax = fmax(KernelMax,HydroNodeForEffSA[ChildNodeID].KernelMax);
            DistanceMax = fmax(DistanceMax,HydroNodeForEffSA[ChildNodeID].DistanceMax);
            NActives += HydroNodeForEffSA[ChildNodeID].NumberofActiveLeaves;
            ChildNodeID = HydroNodeForEffSA[ChildNodeID].Sister;
        } while(ChildNodeID != NONE);

        HydroNodeForEffSA[CurrentNodeID].KernelMax = KernelMax + Width;
        HydroNodeForEffSA[CurrentNodeID].DistanceMax = DistanceMax + Width;
        HydroNodeForEffSA[CurrentNodeID].NumberofActiveLeaves = NActives;
    }

    return;
}

static void PlantImportTreeForEffSA(const int NumberofLeaves, struct StructEffSAExport EffSAExportRecv[restrict]){
    
    InitializeRootForEffSA();

    HydroTreePreprocessingForEffSA(NumberofLeaves,EffSAExportRecv);

    MakeHydroRootForEffSA(NumberofLeaves,EffSAExportRecv); 
    if(NumberofLeaves > 0){
        BuildHydroTreeForEffSA();
    }

    for(int i=1;i<HydroRootForEffSA.NumberofNodes;i++){
        HydroNodeForEffSA[i].Next = NextHydroNodeForEffSA(i);
    }
    DiagHydroForEffSA = DISTANCE(HydroNodeForEffSA[0].Pos,HydroNodeForEffSA[HydroNodeForEffSA[0].Children].Pos);
    HydroNodeDataImplantForEffSA(0);

    return ;
}

/*
 * Unlike the other neighbor search routines, this returns an index list of
 * NBCacheForEffSA. Therefore, you need to use the outcome like 
 * int nbindex = list[i];
 * double nbPos[] = {NBCacheForEffSA[nbindex].Pos[0],
 *                   NBCacheForEffSA[nbindex].Pos[1],NBCacheForEffSA[nbindex].Pos[2]};
 */
static int GetNeighborsPairsForESA(double Pos[restrict], const double h, int list[restrict], const int Index){

    const int NProcs = MPIGetNumProcs();
    double hh = h*h;
    int nlist = 0;

#if 1
    const int RootNodeID = 0;
    int CurrentNodeID = HydroNodeForEffSA[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNodeForEffSA[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNodeForEffSA[CurrentNodeID].Pos,Pos);
#endif
        if( (dx2 > SQ(h+HydroNodeForEffSA[CurrentNodeID].DistanceMax))
          &&(dx2 > SQ(HydroNodeForEffSA[CurrentNodeID].KernelMax)) ){
            CurrentNodeID = HydroNodeForEffSA[CurrentNodeID].Next;
        } else if(HydroNodeForEffSA[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNodeForEffSA[CurrentNodeID].Children;
        } else {
			int NumberofLeaves = HydroNodeForEffSA[CurrentNodeID].NumberofLeaves;
            int header = HydroNodeForEffSA[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                // if(HydroRootForEffSA.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCacheForEffSA[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCacheForEffSA[leaf].Pos,Pos);
#endif
                if( (distance2<hh) || (distance2<4.0*SQ(NBCacheForEffSA[leaf].Kernel)) ){
                    //list[nlist] = NBCacheForEffSA[leaf].Leaf;

                    NBCacheForEffSA[leaf].Interact = true;
                    // This is date for ActiveStellarFeedbackParticle[], so
                    // needs to use Index
                    ExportFlagForGatherScatter[NProcs*Index+NBCacheForEffSA[leaf].DomainID] = true;
                    list[nlist] = leaf;
                    nlist ++;
                }
            }
            CurrentNodeID = HydroNodeForEffSA[CurrentNodeID].Next;
        }
    }
#else
    for(int k=0;k<NEffSAExportRecv;k++){
        double distance2 = DISTANCE2(Pos,EffSAExportRecv[k].Pos);
        if( (distance2<hh) || (distance2<4.0*SQ(EffSAExportRecv[k].Kernel)) ){
            //NBCacheForEffSA[leaf].Interact = true;
            ExportFlagForGatherScatter[NProcs*Index+EffSAExportRecv[k].DomainID] = true;
            list[nlist] = k; ////leaf;
            nlist ++;
        }
    }
    return nlist;

    int NumberofLeaves = HydroNodeForEffSA[0].NumberofLeaves;
    int header = HydroNodeForEffSA[0].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
        double distance2 = DISTANCE2(Pos,NBCacheForEffSA[leaf].Pos);
        if( (distance2<hh) || (distance2<4.0*SQ(NBCacheForEffSA[leaf].Kernel)) ){
            NBCacheForEffSA[leaf].Interact = true;
            ExportFlagForGatherScatter[NProcs*Index+NBCacheForEffSA[leaf].DomainID] = true;
            list[nlist] = leaf;
            nlist ++;
        }
    }

#endif
    return nlist;
}

static int GetNeighborsPairsForESADirect(double Pos[restrict], const double h, int list[restrict], const int Index){

    const int NProcs = MPIGetNumProcs();
    double hh = h*h;
    int nlist = 0;

    int NumberofLeaves = HydroNodeForEffSA[0].NumberofLeaves;
    int header = HydroNodeForEffSA[0].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
#ifdef PERIODIC_RUN 
        double distance2 = 0.e0;
        for(int l=0;l<DIMENSION;l++)
            distance2 += SQ(PeriodicDistance(Pos[l],NBCacheForEffSA[leaf].Pos[l],l));
#else
        double distance2 = DISTANCE2(NBCacheForEffSA[leaf].Pos,Pos);
#endif
        if( (distance2<hh) || (distance2<4.0*SQ(NBCacheForEffSA[leaf].Kernel)) ){
            NBCacheForEffSA[leaf].Interact = true;
            ExportFlagForGatherScatter[NProcs*Index+NBCacheForEffSA[leaf].DomainID] = true;
            list[nlist] = leaf;
            nlist ++;
        }
    }

    return nlist;
}

static inline double __attribute__((always_inline)) EffSADomainDistanceSQR(double Pos[restrict], const int DomainID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForEffSA[DomainID].PosMin[k]) 
            Dist2 += SQ(EdgesForEffSA[DomainID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForEffSA[DomainID].PosMax[k])
            Dist2 += SQ(EdgesForEffSA[DomainID].PosMax[k]-Pos[k]);
    }
    return (Dist2);
}

static inline bool __attribute__((always_inline)) EffSAOverlapDomain(double Pos[restrict], const double h, const int DomainID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForEffSA[DomainID].PosMin[k]) 
            Dist2 += SQ(EdgesForEffSA[DomainID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForEffSA[DomainID].PosMax[k])
            Dist2 += SQ(EdgesForEffSA[DomainID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline double __attribute__((always_inline)) EffSADomainDistanceSQRWithExtent(double Pos[restrict], const int DomainID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForEffSAWithExtent[DomainID].PosMin[k]) 
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForEffSAWithExtent[DomainID].PosMax[k])
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMax[k]-Pos[k]);
    }
    return (Dist2);
}

static inline bool __attribute__((always_inline)) EffSAOverlapDomainWithExtent(double Pos[restrict], const double h, const int DomainID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForEffSAWithExtent[DomainID].PosMin[k]) 
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForEffSAWithExtent[DomainID].PosMax[k])
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline bool __attribute__((always_inline)) EffSAOverlapDomainPositionWithExtent(double Pos[restrict], const int DomainID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForEffSAWithExtent[DomainID].PosMin[k]) 
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForEffSAWithExtent[DomainID].PosMax[k])
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMax[k]-Pos[k]);
    }
    return Dist2>0.e0?false:true;
}

static inline bool __attribute__((always_inline)) EffSAOverlapDomainAABBWithExtent(double Pos[restrict], const double h, const int DomainID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForEffSAWithExtent[DomainID].PosMin[k]) 
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForEffSAWithExtent[DomainID].PosMax[k])
            Dist2 += SQ(EdgesForEffSAWithExtent[DomainID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline bool __attribute__((always_inline)) EffSAOverlapDomainToDomain(const struct StructEdges Edge0, const struct StructEdges Edge1){

    bool flag = true;
    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Edge0.PosMax[k] < Edge1.PosMin[k]){
            Dist2 += SQ(Edge1.PosMin[k]-Edge0.PosMax[k]);
            flag = false;
        }
        if(Edge1.PosMax[k] < Edge0.PosMin[k]){
            Dist2 += SQ(Edge0.PosMin[k]-Edge1.PosMax[k]);
            flag = false;
        }
    }
    return flag;
    return Dist2>0.e0?false:true;
}

static bool FirstCall_CheckEffSAExportFlags = true;

/*
 * This function pick up gas particles which will be exported to the other
 * domain.
 */

static inline int __attribute__((always_inline)) CheckEffSAExportFlags(const int Index, bool EffSAExportFlags[restrict]){

    if(Pall.Nhydro == 0)
        return 0;

    const int MyID = MPIGetMyID();
    int DomainID = CommunicationTable[Index].SendRank;

    if(isinf(EdgesForEffSA[DomainID].PosMin[0])){ // This works correctly.
        return 0;
    }

#if 0
    if(!EffSAOverlapDomainToDomain(EdgesForHydroEffSAWithExtent[MyID],EdgesForEffSA[DomainID])){
        return 0;
    }

    if(!EffSAOverlapDomainToDomain(EdgesForHydroEffSA[MyID],EdgesForEffSAWithExtent[DomainID])){
        return 0;
    }
#endif


#if 0
    if(FirstCall_CheckEffSAExportFlags){
        int NExport = 0;
        int NumberofLeaves = HydroNode[0].NumberofLeaves;
        int header = HydroNode[0].Leaves;
        for(int k=0;k<NumberofLeaves;k++){
            int leaf = header+k;
            EffSAExportFlags[NBCache[leaf].Leaf] = true;
            NExport ++;
        }

        //FirstCall_CheckEffSAExportFlags = false;
        return NExport;
    }
#endif

    int NExport = 0;

#if 0
    {

    int NExport = 0;
    int NumberofLeaves = HydroNode[0].NumberofLeaves;
    int header = HydroNode[0].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
        EffSAExportFlags[leaf] = false;

        double dx2    = EffSADomainDistanceSQR(NBCache[leaf].Pos,DomainID);
        double dx2ext = EffSADomainDistanceSQRWithExtent(NBCache[leaf].Pos,DomainID);

        //if((dx2>SQ(2.0*NBCache[leaf].Kernel))&&(dx2ext>0.e0)) continue;
        //if((dx2ext>SQ(2.0*NBCache[leaf].Kernel))&&(dx2ext>0.e0)) continue;
        if(dx2ext>SQ(2.0*NBCache[leaf].Kernel)) continue;

        EffSAExportFlags[leaf] = true;
        NExport ++;
    }

    // fprintf(stderr,"[%02d] %d / %d \n",MPIGetMyID(),NExport,NumberofLeaves);fflush(NULL);
	return NExport;
    }
#endif


#if 1
    {
    const int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
        //double dx2    = EffSADomainDistanceSQR(HydroNode[CurrentNodeID].Pos,DomainID);
        double dx2ext = EffSADomainDistanceSQRWithExtent(HydroNode[CurrentNodeID].Pos,DomainID);

        //if((SQ(HydroNode[CurrentNodeID].KernelMax)   < dx2)&&
           //(SQ(HydroNode[CurrentNodeID].DistanceMax) < dx2ext)){
        if(SQ(HydroNode[CurrentNodeID].KernelMax) < dx2ext){
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else if(HydroNode[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNode[CurrentNodeID].Children;
        } else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                //EffSAExportFlags[NBCache[leaf].Leaf] = true;
                EffSAExportFlags[leaf] = true;
                NExport ++;
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        }

    }
    return NExport;

    }
#endif

	return NExport;
}

static int EffSAExportFlagsMaxAllocated = 0;
static bool *EffSAExportFlags = NULL;

#if 1
static int ImportHydroForEffSA(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){ 

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    // Send/Recv
    MPI_Status  mpi_status;

    if(EffSAExportFlagsMaxAllocated < Pall.Nhydro){
        EffSAExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        EffSAExportFlags = realloc(EffSAExportFlags,sizeof(bool)*EffSAExportFlagsMaxAllocated);
    }
    
    // counter for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    MPI_Status mpi_status_Export_Send_Precomm[NProcs-1];
    MPI_Request mpi_request_Export_Send_Precomm[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    int SendFlag,RecvFlag;

    struct StructEffSAExport *EffSAExportSend[NProcs];

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(NImportThisTime+i,1,MPI_INT,
            CommunicationTable[i].RecvRank,TAG_EFFSA_PRECOMM,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
    }

    for(int i=0;i<NProcs-1;i++){

        memset(EffSAExportFlags,0,Pall.Nhydro*sizeof(bool));
        NExportThisTime[i] = CheckEffSAExportFlags(i,EffSAExportFlags);

        MPI_Isend(NExportThisTime+i,1,MPI_INT,
            CommunicationTable[i].SendRank,TAG_EFFSA_PRECOMM,
                MPI_COMM_WORLD,mpi_request_Export_Send_Precomm+i);
        MPI_Test(mpi_request_Export_Send_Precomm+i,&SendFlag,MPI_STATUS_IGNORE);

        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructEffSAExport),i);
        EffSAExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        if(NExportThisTime[i] > 0){
            int NumberofLeaves = HydroNode[0].NumberofLeaves;
            int header = HydroNode[0].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(EffSAExportFlags[leaf]){
                    int index = NBCache[leaf].Leaf;
                    EffSAExportSend[i][NExport].Kernel = Phydro[index]->KernelPred;
                    // EffSAExportSend[i][NExport].Pos[0] = PhydroPosP(index)[0];
                    // EffSAExportSend[i][NExport].Pos[1] = PhydroPosP(index)[1];
                    // EffSAExportSend[i][NExport].Pos[2] = PhydroPosP(index)[2];
                    EffSAExportSend[i][NExport].Pos[0] = NBCache[leaf].Pos[0];
                    EffSAExportSend[i][NExport].Pos[1] = NBCache[leaf].Pos[1];
                    EffSAExportSend[i][NExport].Pos[2] = NBCache[leaf].Pos[2];
                    EffSAExportSend[i][NExport].NumberDensity = Phydro[index]->NumberDensity;
                    EffSAExportSend[i][NExport].DomainID = MyID;
                    EffSAExportSend[i][NExport].Leaf = index;
                    EffSAExportSend[i][NExport].GlobalID = PhydroBody(index)->GlobalID;
                    NExport ++;
                } 
            }
        }
    }
    MPI_Waitall(NProcs-1,mpi_request_Export_Send_Precomm,mpi_status_Export_Send_Precomm); // NExportThisTime
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv); // NImportThisTime


    int counter_send = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i] > 0){
            MPI_Isend(EffSAExportSend[i],
                NExportThisTime[i]*sizeof(struct StructEffSAExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_EFFSA_EXPORT,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
    }

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportAll += NImportThisTime[i];
    }


    //struct StructHydroDensityExport *HydroDensityExportRecv;
    //struct StructHydroDensityImport *HydroDensityImportSend;
    //CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHydroDensityExport));

    // Prepare buffer
#if 0
    if(EffSAExportRecvAllocated < NImportAll){
        EffSAExportRecvAllocated = ForAngelsShare*NImportAll;
        EffSAExportRecv = realloc(EffSAExportRecv,
                sizeof(struct StructEffSAExport)*EffSAExportRecvAllocated);
    }
#else
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructEffSAExport));
    EffSAExportRecv = BufferExportRecv;
#endif


    int NImport = 0;
    int counter_recv = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NImportThisTime[i]>0){
            MPI_Irecv(EffSAExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructEffSAExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_EFFSA_EXPORT,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }
    NEffSAExportRecv = NImportAll;

    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send); // wait for the jdata send
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv); // wait for the data recive 

    
    // Write imported data
#if 0
    {
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"Imported.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
        FileOpen(fp,fname,"w");
    
        for(int i=0;i<NImport;i++){
            fprintf(fp,"%g %g %g %lg\n",
                    EffSAExportRecv[i].Pos[0],
                    EffSAExportRecv[i].Pos[1],
                    EffSAExportRecv[i].Pos[2],
                    EffSAExportRecv[i].GlobalID);
        }

        fclose(fp);

    }
#endif


    return NImportAll;
}
#else

static int NEffSAExportSendAllocated = 0;
static struct StructEffSAExport *EffSAExportSend;
static int ImportHydroForEffSA(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){ 

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    // Send/Recv
    MPI_Status  mpi_status;

    if(EffSAExportFlagsMaxAllocated < Pall.Nhydro){
        EffSAExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        EffSAExportFlags = realloc(EffSAExportFlags,sizeof(bool)*EffSAExportFlagsMaxAllocated);
    }
    
    // counter for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    MPI_Status mpi_status_Export_Send_Precomm[NProcs-1];
    MPI_Request mpi_request_Export_Send_Precomm[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    int SendFlag,RecvFlag;

    //struct StructEffSAExport *EffSAExportSend[NProcs];
    //struct StructEffSAExport *EffSAExportSend[NProcs];

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(NImportThisTime+i,1,MPI_INT,
            CommunicationTable[i].RecvRank,TAG_EFFSA_PRECOMM,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
    }

    int counter_send = 0;
    for(int i=0;i<NProcs-1;i++){

        memset(EffSAExportFlags,0,Pall.Nhydro*sizeof(bool));
        NExportThisTime[i] = CheckEffSAExportFlags(i,EffSAExportFlags);

        MPI_Isend(NExportThisTime+i,1,MPI_INT,
            CommunicationTable[i].SendRank,TAG_EFFSA_PRECOMM,
                MPI_COMM_WORLD,mpi_request_Export_Send_Precomm+i);
        MPI_Test(mpi_request_Export_Send_Precomm+i,&SendFlag,MPI_STATUS_IGNORE);

        //CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructEffSAExport),i);
        //EffSAExportSend[i] = BufferExportSend[i];
        if(NExportThisTime[i] > NEffSAExportSendAllocated){
            NEffSAExportSendAllocated = MAX(ForAngelsShare*NExportThisTime[i],NAdditionUnit);
            EffSAExportSend = realloc(EffSAExportSend,sizeof(struct StructEffSAExport)*NNEffSAExportSendAllocated);
        }

        int NExport = 0;
        if(NExportThisTime[i] > 0){
#if 0
            for(int k=0;k<Pall.Nhydro;k++){
                if(EffSAExportFlags[k]){ 
                    EffSAExportSend[i][NExport].Kernel = Phydro[k]->KernelPred;
                    EffSAExportSend[i][NExport].Pos[0] = PhydroPosP(k)[0];
                    EffSAExportSend[i][NExport].Pos[1] = PhydroPosP(k)[1];
                    EffSAExportSend[i][NExport].Pos[2] = PhydroPosP(k)[2];
                    EffSAExportSend[i][NExport].NumberDensity = Phydro[k]->NumberDensity;
                    EffSAExportSend[i][NExport].DomainID = MyID;
                    EffSAExportSend[i][NExport].Leaf = k;
                    EffSAExportSend[i][NExport].GlobalID = PhydroBody(k)->GlobalID;
                    NExport ++;
                } 
            }
#else
            int NumberofLeaves = HydroNode[0].NumberofLeaves;
            int header = HydroNode[0].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(EffSAExportFlags[leaf]){
                    int index = NBCache[leaf].Leaf;
                    EffSAExportSend[NExport].Kernel = Phydro[index]->KernelPred;
                    EffSAExportSend[NExport].Pos[0] = PhydroPosP(index)[0];
                    EffSAExportSend[NExport].Pos[1] = PhydroPosP(index)[1];
                    EffSAExportSend[NExport].Pos[2] = PhydroPosP(index)[2];
                    EffSAExportSend[NExport].NumberDensity = Phydro[index]->NumberDensity;
                    EffSAExportSend[NExport].DomainID = MyID;
                    EffSAExportSend[NExport].Leaf = index;
                    EffSAExportSend[NExport].GlobalID = PhydroBody(index)->GlobalID;
                    NExport ++;
                } 
            }
#endif
            MPI_Isend(EffSAExportSend,
                NExportThisTime[i]*sizeof(struct StructEffSAExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_EFFSA_EXPORT,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
    }

    MPI_Waitall(NProcs-1,mpi_request_Export_Send_Precomm,mpi_status_Export_Send_Precomm); // NExportThisTime
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv); // NImportThisTime

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportAll += NImportThisTime[i];
    }

    // Prepare buffer
    if(EffSAExportRecvAllocated < NImportAll){
        EffSAExportRecvAllocated = ForAngelsShare*NImportAll;
        EffSAExportRecv = realloc(EffSAExportRecv,
                sizeof(struct StructEffSAExport)*EffSAExportRecvAllocated);
    }


    int NImport = 0;
    int counter_recv = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NImportThisTime[i]>0){
            MPI_Irecv(EffSAExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructEffSAExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_EFFSA_EXPORT,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }
    NEffSAExportRecv = NImportAll;

    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send); // wait for the jdata send
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv); // wait for the data recive 

    
    // Write imported data
#if 0
    {
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"Imported.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
        FileOpen(fp,fname,"w");
    
        for(int i=0;i<NImport;i++){
            fprintf(fp,"%g %g %g %lg\n",
                    EffSAExportRecv[i].Pos[0],
                    EffSAExportRecv[i].Pos[1],
                    EffSAExportRecv[i].Pos[2],
                    EffSAExportRecv[i].GlobalID);
        }

        fclose(fp);

    }
#endif


    return NImportAll;
}
#endif


/*
 * This is tree for imported hydro data.
 */
static void PlantTreeForEffSA(const int NumberofLeaves, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){ 
    if(NumberofLeaves > 0){
        PlantImportTreeForEffSA(NumberofLeaves,EffSAExportRecv);
    }

    return ;
}


static void PickupTargets(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){

    const int NProcs = MPIGetNumProcs();

    if(ActiveIndexListAllocated < NExplosion){
        ActiveIndexListAllocated = MAX(ForAngelsShare*NExplosion,NAdditionUnit);
        ActiveIndexList = realloc(ActiveIndexList,sizeof(int)*ActiveIndexListAllocated);
        ExportFlagForGatherScatter = realloc(ExportFlagForGatherScatter,sizeof(bool)*(ActiveIndexListAllocated*NProcs));
        ExportFlagForGatherScatterIndexList = realloc(ExportFlagForGatherScatterIndexList,sizeof(int)*ActiveIndexListAllocated);
    }

    memset(ExportFlagForGatherScatter,0,sizeof(bool)*NProcs*ActiveIndexListAllocated);

    int counter = 0;
    for(int i=0;i<NExplosion;i++){
        ExportFlagForGatherScatterIndexList[i] = NONE;
        if((ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_SNII)||
           (ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_SNIa)
#ifdef USE_STELLAR_WIND //{
           || (ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_SW)
#endif // USE_STELLAR_WIND //}
          ){ // type II/Ia
            ActiveIndexList[counter] = i;
            ExportFlagForGatherScatterIndexList[i] = counter;
            counter ++;
        }
    }

    CurrentActiveIndexListSize = counter;

    return ;
}


void CalcEffectiveSurfaceAreaPrev(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){

    NExplosionLog = NExplosion;

    PickupTargets(NExplosion,ActiveStellarFeedbackParticle);

    CalcActiveDomainEdgesForESA(ActiveStellarFeedbackParticle); 
    CalcHydroDomainEdgesForEffSA();

    int NumberofLeaves = ImportHydroForEffSA(NExplosion,ActiveStellarFeedbackParticle); 

    PlantTreeForEffSA(NumberofLeaves,ActiveStellarFeedbackParticle);

    return ;
}

static struct StructEffSAResult ReturnStructureEffSAEngine(struct StructEffSAInput Input){

    double Pos[3];
    Pos[0] = Input.Pos[0];
    Pos[1] = Input.Pos[1];
    Pos[2] = Input.Pos[2];
    double Kerneli = Input.Kernel;
    double InvKerneli = 1.e0/Kerneli;
    double NumberDensityi = Input.NumberDensity;

    int Neighbors[MaxNeighborSize];
    int Nlist = 0;
    if(Pall.Nhydro>0)
        Nlist = GetNeighborsPairsLimitedNBCacheIndex(Pos,2.e0*Kerneli,Neighbors);
        //Nlist = GetNeighborsPairsLimited(Pos,2.e0*Kerneli,Neighbors);

    struct StructEffSAResult TempEffSAResult = {0};
    memset(&TempEffSAResult,0,sizeof(struct StructEffSAResult));

	for(int i=0;i<Nlist;i++){
        int index = Neighbors[i];
        int leaf = NBCache[index].Leaf;
        double xij[3];

#if 0
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],PhydroPosP(leaf)[0],0);
        xij[1] = PeriodicDistance(Pos[1],PhydroPosP(leaf)[1],1);
        xij[2] = PeriodicDistance(Pos[2],PhydroPosP(leaf)[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-PhydroPosP(leaf)[0];
        xij[1] = Pos[1]-PhydroPosP(leaf)[1];
        xij[2] = Pos[2]-PhydroPosP(leaf)[2];
#endif // PERIODIC_RUN //}
#endif
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],NBCache[index].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],NBCache[index].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],NBCache[index].Pos[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-NBCache[index].Pos[0];
        xij[1] = Pos[1]-NBCache[index].Pos[1];
        xij[2] = Pos[2]-NBCache[index].Pos[2];
#endif // PERIODIC_RUN //}

        double r = NORM(xij);
		// double w = SPHKernel(r,InvKerneli);

        double InvKernelj = 1.0/Phydro[leaf]->KernelPred; 
        double dwi = dSPHKernel(r,InvKerneli);
        double dwj = dSPHKernel(r,InvKernelj);

        // Note that
        // dW/dr = dSPHKernel()*r;
#ifdef USE_GRAD_N //{
        double NumberDensityj = Phydro[leaf]->NumberDensity;
#else // USE_GRAD_N //}//{
        double NumberDensityj = Phydro[leaf]->RhoPred/Phydro[leaf]->Mass; //1/Vi
#endif // USE_GRAD_N //}
        double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj));

        double wk = 0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r))));

        double wk_vec[EffSAVecSize] = {0.e0};
        if(xij[0]>0){
            wk_vec[0] = wk*xij[0]/r;
            wk_vec[1] = 0.0;
        }else{
            wk_vec[0] = 0.0;
            wk_vec[1] = wk*xij[0]/r;
        }
        if(xij[1]>0){
            wk_vec[2] = wk*xij[1]/r;
            wk_vec[3] = 0.0;
        }else{
            wk_vec[2] = 0.0;
            wk_vec[3] = wk*xij[1]/r;
        }
        if(xij[2]>0){
            wk_vec[4] = wk*xij[2]/r;
            wk_vec[5] = 0.0;
        } else {
            wk_vec[4] = 0.0;
            wk_vec[5] = wk*xij[2]/r;
        }

        TempEffSAResult.SA += wk;
        for(int k=0;k<6;k++)
            TempEffSAResult.WeightCorrection[k] += wk_vec[k];
    }
    TempEffSAResult.Nlist = Nlist;


    if(NEffSAExportRecv == 0) return TempEffSAResult;

    Nlist = GetNeighborsPairsForESA(Pos,2.e0*Kerneli,Neighbors,Input.Index);

	for(int i=0;i<Nlist;i++){
        int leaf = Neighbors[i];
        double xij[3];

#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],NBCacheForEffSA[leaf].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],NBCacheForEffSA[leaf].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],NBCacheForEffSA[leaf].Pos[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-NBCacheForEffSA[leaf].Pos[0];
        xij[1] = Pos[1]-NBCacheForEffSA[leaf].Pos[1];
        xij[2] = Pos[2]-NBCacheForEffSA[leaf].Pos[2];
#endif // PERIODIC_RUN //}


        double r = NORM(xij);
		// double w = SPHKernel(r,InvKerneli);

        double InvKernelj = 1.0/NBCacheForEffSA[leaf].Kernel; 
        double dwi = dSPHKernel(r,InvKerneli);
        double dwj = dSPHKernel(r,InvKernelj);

        // Note that
        // dW/dr = dSPHKernel()*r;
#ifdef USE_GRAD_N //{
        double NumberDensityj = NBCacheForEffSA[leaf].NumberDensity;
#else // USE_GRAD_N //}//{
        int nbindex = NBCacheForEffSA[leaf].Leaf;
#warning this might not be correct.
        double NumberDensityj = Phydro[nbindex]->RhoPred/Phydro[nbindex]->Mass; //1/Vi
#error This branch is unavilable.
#endif // USE_GRAD_N //}
        double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj));

        double wk = 0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r))));

        double wk_vec[EffSAVecSize] = {0.e0};
        if(xij[0]>0){
            wk_vec[0] = wk*xij[0]/r;
            wk_vec[1] = 0.0;
        }else{
            wk_vec[0] = 0.0;
            wk_vec[1] = wk*xij[0]/r;
        }
        if(xij[1]>0){
            wk_vec[2] = wk*xij[1]/r;
            wk_vec[3] = 0.0;
        }else{
            wk_vec[2] = 0.0;
            wk_vec[3] = wk*xij[1]/r;
        }
        if(xij[2]>0){
            wk_vec[4] = wk*xij[2]/r;
            wk_vec[5] = 0.0;
        } else {
            wk_vec[4] = 0.0;
            wk_vec[5] = wk*xij[2]/r;
        }

        TempEffSAResult.SA += wk;
        for(int k=0;k<6;k++)
            TempEffSAResult.WeightCorrection[k] += wk_vec[k];
    }
    TempEffSAResult.Nlist += Nlist;

    return TempEffSAResult;
}

void CalcEffectiveSurfaceAreaSum(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){

    int NActives = CurrentActiveIndexListSize;
    for(int i=0;i<NActives;i++){ 
        int leaf = ActiveIndexList[i];
        struct StructEffSAResult TempEffSAResult =
            ReturnStructureEffSAEngine((struct StructEffSAInput){
                    .Pos[0]=ActiveStellarFeedbackParticle[leaf].Pos[0],
                    .Pos[1]=ActiveStellarFeedbackParticle[leaf].Pos[1],
                    .Pos[2]=ActiveStellarFeedbackParticle[leaf].Pos[2],
                    .Kernel=ActiveStellarFeedbackParticle[leaf].Radius,
                    .NumberDensity=ActiveStellarFeedbackParticle[leaf].NumberDensity,
                    .Index=leaf});

        ActiveStellarFeedbackParticle[leaf].Nlist = TempEffSAResult.Nlist;
        ActiveStellarFeedbackParticle[leaf].SA = TempEffSAResult.SA;
        ActiveStellarFeedbackParticle[leaf].WeightCorrection[0] = TempEffSAResult.WeightCorrection[0];
        ActiveStellarFeedbackParticle[leaf].WeightCorrection[1] = TempEffSAResult.WeightCorrection[1];
        ActiveStellarFeedbackParticle[leaf].WeightCorrection[2] = TempEffSAResult.WeightCorrection[2];
        ActiveStellarFeedbackParticle[leaf].WeightCorrection[3] = TempEffSAResult.WeightCorrection[3];
        ActiveStellarFeedbackParticle[leaf].WeightCorrection[4] = TempEffSAResult.WeightCorrection[4];
        ActiveStellarFeedbackParticle[leaf].WeightCorrection[5] = TempEffSAResult.WeightCorrection[5];
    } 

    return;

}

static struct StructEffSAResult ReturnStructureEffSAVecEngine(struct StructEffSAInput Input){

    double Pos[3];
    Pos[0] = Input.Pos[0];
    Pos[1] = Input.Pos[1];
    Pos[2] = Input.Pos[2];
    double Kerneli = Input.Kernel;
    double InvKerneli = 1.e0/Kerneli;
    double NumberDensityi = Input.NumberDensity;

    double wk_norm = 1.0/(Input.SA);

    int Neighbors[MaxNeighborSize];
    int Nlist = 0;
    if(Pall.Nhydro>0)
        Nlist = GetNeighborsPairsLimitedNBCacheIndex(Pos,2.e0*Kerneli,Neighbors);


    struct StructEffSAResult TempEffSAResult = {0};

	for(int i=0;i<Nlist;i++){
        int index = Neighbors[i];
        int leaf = NBCache[index].Leaf;
        double xij[3];
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],NBCache[index].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],NBCache[index].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],NBCache[index].Pos[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-NBCache[index].Pos[0];
        xij[1] = Pos[1]-NBCache[index].Pos[1];
        xij[2] = Pos[2]-NBCache[index].Pos[2];
#endif // PERIODIC_RUN //}

        double r = NORM(xij);
		//double w = SPHKernel(r,InvKerneli);

        double InvKernelj = 1.0/Phydro[leaf]->KernelPred; 
        double dwi = dSPHKernel(r,InvKerneli);
        double dwj = dSPHKernel(r,InvKernelj);

        // Note that
        // dW/dr = dSPHKernel()*r;
#ifdef USE_GRAD_N //{
        double NumberDensityj = Phydro[leaf]->NumberDensity;
#else // USE_GRAD_N //}//{
        double NumberDensityj = Phydro[leaf]->RhoPred/Phydro[leaf]->Mass; //1/Vi
#endif // USE_GRAD_N //}
        double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj));
        double wk = 0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r))));

        double wk_vec[EffSAVecSize] = {0.e0};
        if(xij[0]>0){
            wk_vec[0] = wk*xij[0]/r;
            wk_vec[1] = 0.0;
        }else{
            wk_vec[0] = 0.0;
            wk_vec[1] = wk*xij[0]/r;
        }
        if(xij[1]>0){
            wk_vec[2] = wk*xij[1]/r;
            wk_vec[3] = 0.0;
        }else{
            wk_vec[2] = 0.0;
            wk_vec[3] = wk*xij[1]/r;
        }
        if(xij[2]>0){
            wk_vec[4] = wk*xij[2]/r;
            wk_vec[5] = 0.0;
        } else {
            wk_vec[4] = 0.0;
            wk_vec[5] = wk*xij[2]/r;
        }

        wk *= wk_norm;

        //TempEffSAResult.CheckWeight += wk;

        double pnorm = 0.e0;
        double pvec[3] = {0.e0};
        for(int k=0;k<3;k++){
            double q = 0.e0;
            int i1 = 2*k;
            int i2 = i1+1;
            double q_i1 = fabs(Input.WeightCorrection[i1]);
            double q_i2 = fabs(Input.WeightCorrection[i2]);
            if((q_i1>0.e0)&&(q_i2>0.e0)){
                double rr = q_i2/q_i1;
                double rr2 = SQ(rr);
                //if(wk_vec[i1] != 0.e0){
                if(fpclassify(fabs(wk_vec[i1])) != FP_ZERO){
                    q += wk_norm*wk_vec[i1]*sqrt(0.5*(1.e0+rr2));
                }else{
                    q += wk_norm*wk_vec[i2]*sqrt(0.5*(1.e0+1.e0/rr2));
                }
            }else{
                q += wk_norm*(wk_vec[i1]+wk_vec[i2]);
            }
            pvec[k] = -q;
            pnorm += SQ(pvec[k]);
        }
        pnorm = sqrt(pnorm);
        TempEffSAResult.WeightCorrection[6] += pnorm;
    }
    TempEffSAResult.Nlist = Nlist;

    if(NEffSAExportRecv == 0) return TempEffSAResult;

    Nlist = GetNeighborsPairsForESA(Pos,2.e0*Kerneli,Neighbors,Input.Index);

	for(int i=0;i<Nlist;i++){
        int leaf = Neighbors[i];
        double xij[3];
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],NBCacheForEffSA[leaf].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],NBCacheForEffSA[leaf].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],NBCacheForEffSA[leaf].Pos[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-NBCacheForEffSA[leaf].Pos[0];
        xij[1] = Pos[1]-NBCacheForEffSA[leaf].Pos[1];
        xij[2] = Pos[2]-NBCacheForEffSA[leaf].Pos[2];
#endif // PERIODIC_RUN //}

        double r = NORM(xij);

        double InvKernelj = 1.0/NBCacheForEffSA[leaf].Kernel; 
        double dwi = dSPHKernel(r,InvKerneli);
        double dwj = dSPHKernel(r,InvKernelj);

        // Note that
        // dW/dr = dSPHKernel()*r;
#ifdef USE_GRAD_N //{
        double NumberDensityj = NBCacheForEffSA[leaf].NumberDensity;
#else // USE_GRAD_N //}//{
        int nbindex = NBCacheForEffSA[leaf].Leaf;
#warning This part might be incorrect.
        double NumberDensityj = Phydro[nbindex]->RhoPred/Phydro[nbindex]->Mass; //1/Vi
#endif // USE_GRAD_N //}
        double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj));
        double wk = 0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r))));

        double wk_vec[EffSAVecSize] = {0.e0};
        if(xij[0]>0){
            wk_vec[0] = wk*xij[0]/r;
            wk_vec[1] = 0.0;
        }else{
            wk_vec[0] = 0.0;
            wk_vec[1] = wk*xij[0]/r;
        }
        if(xij[1]>0){
            wk_vec[2] = wk*xij[1]/r;
            wk_vec[3] = 0.0;
        }else{
            wk_vec[2] = 0.0;
            wk_vec[3] = wk*xij[1]/r;
        }
        if(xij[2]>0){
            wk_vec[4] = wk*xij[2]/r;
            wk_vec[5] = 0.0;
        } else {
            wk_vec[4] = 0.0;
            wk_vec[5] = wk*xij[2]/r;
        }

        wk *= wk_norm;

        //TempEffSAResult.CheckWeight += wk;

        double pnorm = 0.e0;
        double pvec[3] = {0.e0};
        for(int k=0;k<3;k++){
            double q = 0.e0;
            int i1 = 2*k;
            int i2 = i1+1;
            double q_i1 = fabs(Input.WeightCorrection[i1]);
            double q_i2 = fabs(Input.WeightCorrection[i2]);
            if((q_i1>0.e0)&&(q_i2>0.e0)){
                double rr = q_i2/q_i1;
                double rr2 = SQ(rr);
                //if(wk_vec[i1] != 0.e0){
                if(fpclassify(fabs(wk_vec[i1])) != FP_ZERO){
                    q += wk_norm*wk_vec[i1]*sqrt(0.5*(1.e0+rr2));
                }else{
                    q += wk_norm*wk_vec[i2]*sqrt(0.5*(1.e0+1.e0/rr2));
                }
            }else{
                q += wk_norm*(wk_vec[i1]+wk_vec[i2]);
            }
            pvec[k] = -q;
            pnorm += SQ(pvec[k]);
        }
        pnorm = sqrt(pnorm);
        TempEffSAResult.WeightCorrection[6] += pnorm;
    }
    TempEffSAResult.Nlist += Nlist;

    return TempEffSAResult;
}

void CalcEffectiveSurfaceAreaVec(const int NExplosion, 
        struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){ 

    int NActives = CurrentActiveIndexListSize;
    for(int i=0;i<NActives;i++){
        int leaf = ActiveIndexList[i];

        struct StructEffSAResult TempEffSAResult =
            ReturnStructureEffSAVecEngine((struct StructEffSAInput){
                    .Pos[0]=ActiveStellarFeedbackParticle[leaf].Pos[0],
                    .Pos[1]=ActiveStellarFeedbackParticle[leaf].Pos[1],
                    .Pos[2]=ActiveStellarFeedbackParticle[leaf].Pos[2],
                    .Kernel=ActiveStellarFeedbackParticle[leaf].Radius,
                    .NumberDensity=ActiveStellarFeedbackParticle[leaf].NumberDensity,
                    .SA=ActiveStellarFeedbackParticle[leaf].SA,
                    .WeightCorrection[0]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[0],
                    .WeightCorrection[1]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[1],
                    .WeightCorrection[2]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[2],
                    .WeightCorrection[3]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[3],
                    .WeightCorrection[4]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[4],
                    .WeightCorrection[5]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[5],});

        ActiveStellarFeedbackParticle[leaf].WeightCorrection[6] = TempEffSAResult.WeightCorrection[6];
        ActiveStellarFeedbackParticle[leaf].CheckWeight = TempEffSAResult.CheckWeight;
    } 

    return ;
}

static struct StructEffSAResult ReturnStructureEffSAVecCheckEngine(struct StructEffSAInput Input){

    double Pos[3];
    Pos[0] = Input.Pos[0];
    Pos[1] = Input.Pos[1];
    Pos[2] = Input.Pos[2];
    double Kerneli = Input.Kernel;
    double InvKerneli = 1.e0/Kerneli;
    double NumberDensityi = Input.NumberDensity;

    double wk_norm = 1.0/(Input.SA);
    double pw_norm = 1.0/(Input.WeightCorrection[6]);

    int Neighbors[MaxNeighborSize];
    int Nlist = GetNeighborsPairsLimited(Pos,2.e0*Kerneli,Neighbors);

    struct StructEffSAResult TempEffSAResult = {0};

	for(int i=0;i<Nlist;i++){
        int leaf = Neighbors[i];
        double xij[3];
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],Phydro[leaf]->PosP[0],0);
        xij[1] = PeriodicDistance(Pos[1],Phydro[leaf]->PosP[1],1);
        xij[2] = PeriodicDistance(Pos[2],Phydro[leaf]->PosP[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-Phydro[leaf]->PosP[0];
        xij[1] = Pos[1]-Phydro[leaf]->PosP[1];
        xij[2] = Pos[2]-Phydro[leaf]->PosP[2];
#endif // PERIODIC_RUN //}

        double r = NORM(xij);
		//double w = SPHKernel(r,InvKerneli);

        double InvKernelj = 1.0/Phydro[leaf]->KernelPred; 
        double dwi = dSPHKernel(r,InvKerneli);
        double dwj = dSPHKernel(r,InvKernelj);

        // Note that
        // dW/dr = dSPHKernel()*r;
#ifdef USE_GRAD_N //{
        double NumberDensityj = Phydro[leaf]->NumberDensity;
#else // USE_GRAD_N //}//{
        double NumberDensityj = Phydro[leaf]->RhoPred/Phydro[leaf]->Mass; //1/Vi
#endif // USE_GRAD_N //}
        double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj));
        //double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(Phydro[leaf]->NumberDensity));
        double wk = 0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r))));

        double wk_vec[EffSAVecSize] = {0.e0};
        if(xij[0]>0){
            wk_vec[0] = wk*xij[0]/r;
            wk_vec[1] = 0.0;
        }else{
            wk_vec[0] = 0.0;
            wk_vec[1] = wk*xij[0]/r;
        }
        if(xij[1]>0){
            wk_vec[2] = wk*xij[1]/r;
            wk_vec[3] = 0.0;
        }else{
            wk_vec[2] = 0.0;
            wk_vec[3] = wk*xij[1]/r;
        }
        if(xij[2]>0){
            wk_vec[4] = wk*xij[2]/r;
            wk_vec[5] = 0.0;
        } else {
            wk_vec[4] = 0.0;
            wk_vec[5] = wk*xij[2]/r;
        }

        wk *= wk_norm;

        //TempEffSAResult.CheckWeight += wk;

        double pnorm = 0.e0;
        double pvec[3] = {0.e0};
        for(int k=0;k<3;k++){
            double q = 0.e0;
            int i1 = 2*k;
            int i2 = i1+1;
            double q_i1 = fabs(Input.WeightCorrection[i1]);
            double q_i2 = fabs(Input.WeightCorrection[i2]);
            if((q_i1>0.e0)&&(q_i2>0.e0)){
                double rr = q_i2/q_i1;
                double rr2 = SQ(rr);
                //if(wk_vec[i1] != 0.e0){
                if(fpclassify(fabs(wk_vec[i1])) != FP_ZERO){
                    q += wk_norm*wk_vec[i1]*sqrt(0.5*(1.e0+rr2));
                }else{
                    q += wk_norm*wk_vec[i2]*sqrt(0.5*(1.e0+1.e0/rr2));
                }
            }else{
                q += wk_norm*(wk_vec[i1]+wk_vec[i2]);
            }
            pvec[k] = -q;
            pnorm += SQ(pvec[k]);
        }
        pnorm = sqrt(pnorm);
        TempEffSAResult.WeightCorrection[6] += pnorm;

        pnorm *= pw_norm;
        wk = pnorm;
        TempEffSAResult.CheckWeight += wk;

    }
    TempEffSAResult.Nlist = Nlist;


    Nlist = GetNeighborsPairsForESA(Pos,2.e0*Kerneli,Neighbors,Input.Index);
	for(int i=0;i<Nlist;i++){
        int leaf = Neighbors[i];
        double xij[3];
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],NBCacheForEffSA[leaf].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],NBCacheForEffSA[leaf].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],NBCacheForEffSA[leaf].Pos[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-NBCacheForEffSA[leaf].Pos[0];
        xij[1] = Pos[1]-NBCacheForEffSA[leaf].Pos[1];
        xij[2] = Pos[2]-NBCacheForEffSA[leaf].Pos[2];
#endif // PERIODIC_RUN //}

        double r = NORM(xij);
		//double w = SPHKernel(r,InvKerneli);

        double InvKernelj = 1.0/NBCacheForEffSA[leaf].Kernel; 
        double dwi = dSPHKernel(r,InvKerneli);
        double dwj = dSPHKernel(r,InvKernelj);

        // Note that
        // dW/dr = dSPHKernel()*r;
#ifdef USE_GRAD_N //{
        double NumberDensityj = NBCacheForEffSA[leaf].NumberDensity;
#else // USE_GRAD_N //}//{
        int nbindex = NBCacheForEffSA[leaf].Leaf;
        double NumberDensityj = Phydro[nbindex]->RhoPred/Phydro[nbindex]->Mass; //1/Vi
#endif // USE_GRAD_N //}
        double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj));
        //double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(Phydro[leaf]->NumberDensity));
        double wk = 0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r))));

        double wk_vec[EffSAVecSize] = {0.e0};
        if(xij[0]>0){
            wk_vec[0] = wk*xij[0]/r;
            wk_vec[1] = 0.0;
        }else{
            wk_vec[0] = 0.0;
            wk_vec[1] = wk*xij[0]/r;
        }
        if(xij[1]>0){
            wk_vec[2] = wk*xij[1]/r;
            wk_vec[3] = 0.0;
        }else{
            wk_vec[2] = 0.0;
            wk_vec[3] = wk*xij[1]/r;
        }
        if(xij[2]>0){
            wk_vec[4] = wk*xij[2]/r;
            wk_vec[5] = 0.0;
        } else {
            wk_vec[4] = 0.0;
            wk_vec[5] = wk*xij[2]/r;
        }

        wk *= wk_norm;

        //TempEffSAResult.CheckWeight += wk;

        double pnorm = 0.e0;
        double pvec[3] = {0.e0};
        for(int k=0;k<3;k++){
            double q = 0.e0;
            int i1 = 2*k;
            int i2 = i1+1;
            double q_i1 = fabs(Input.WeightCorrection[i1]);
            double q_i2 = fabs(Input.WeightCorrection[i2]);
            if((q_i1>0.e0)&&(q_i2>0.e0)){
                double rr = q_i2/q_i1;
                double rr2 = SQ(rr);
                //if(wk_vec[i1] != 0.e0){
                if(fpclassify(fabs(wk_vec[i1])) != FP_ZERO){
                    q += wk_norm*wk_vec[i1]*sqrt(0.5*(1.e0+rr2));
                }else{
                    q += wk_norm*wk_vec[i2]*sqrt(0.5*(1.e0+1.e0/rr2));
                }
            }else{
                q += wk_norm*(wk_vec[i1]+wk_vec[i2]);
            }
            pvec[k] = -q;
            pnorm += SQ(pvec[k]);
        }
        pnorm = sqrt(pnorm);
        TempEffSAResult.WeightCorrection[6] += pnorm;

        pnorm *= pw_norm;
        wk = pnorm;
        TempEffSAResult.CheckWeight += wk;

    }
    TempEffSAResult.Nlist += Nlist;

    return TempEffSAResult;
}

void CalcEffectiveSurfaceAreaVecCheck(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle){ 

    int NActives = CurrentActiveIndexListSize;
    for(int i=0;i<NActives;i++){
        int leaf = ActiveIndexList[i];

        struct StructEffSAResult TempEffSAResult =
            ReturnStructureEffSAVecCheckEngine((struct StructEffSAInput){
                    .Pos[0]=ActiveStellarFeedbackParticle[leaf].Pos[0],
                    .Pos[1]=ActiveStellarFeedbackParticle[leaf].Pos[1],
                    .Pos[2]=ActiveStellarFeedbackParticle[leaf].Pos[2],
                    .Kernel=ActiveStellarFeedbackParticle[leaf].Radius,
                    .NumberDensity=ActiveStellarFeedbackParticle[leaf].NumberDensity,
                    .SA=ActiveStellarFeedbackParticle[leaf].SA,
                    .WeightCorrection[0]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[0],
                    .WeightCorrection[1]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[1],
                    .WeightCorrection[2]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[2],
                    .WeightCorrection[3]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[3],
                    .WeightCorrection[4]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[4],
                    .WeightCorrection[5]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[5],
                    .WeightCorrection[6]=ActiveStellarFeedbackParticle[leaf].WeightCorrection[6],});

        ActiveStellarFeedbackParticle[leaf].CheckWeight = TempEffSAResult.CheckWeight;
    } 

#if 0
    const int NProcs = MPIGetNumProcs();

    for(int i=0;i<NProcs-1;i++){
        for(int k=0;k<NExportThisTime[i];k++){
            if(EffSAResult[i][k].Nlist==0) continue;
            int leaf = EffSAResult[i][k].Leaf;
            ActiveStellarFeedbackParticle[leaf].CheckWeight += EffSAResult[i][k].CheckWeight;
        }
    }
#endif

#if 0
    for(int i=0;i<NActives;i++){  // local
        int leaf = ActiveIndexList[i];
        fprintf(stderr,"Leaf %d; 3 Nlist %d, Weight = %g\n",leaf,
                ActiveStellarFeedbackParticle[leaf].Nlist,
                ActiveStellarFeedbackParticle[leaf].CheckWeight);
    }
#endif

    return ;
}

//bool CheckExporFlagGatherScatter(const int Offset){
    //return ExportFlagForGatherScatter[Offset];
//}

bool CheckExporFlagGatherScatter(const int Index, const int SendRank){
    const int NProcs = MPIGetNumProcs();
    //assert(ExportFlagForGatherScatterIndexList[Index]!=NONE);
    // if(ExportFlagForGatherScatterIndexList[Index]==NONE) return false;
    // return ExportFlagForGatherScatter[NProcs*ExportFlagForGatherScatterIndexList[Index]+SendRank];
    //return ExportFlagForGatherScatterIndexList[Index]==NONE?false:ExportFlagForGatherScatter[NProcs*ExportFlagForGatherScatterIndexList[Index]+SendRank];
    return ExportFlagForGatherScatterIndexList[Index]==NONE?false:ExportFlagForGatherScatter[NProcs*Index+SendRank];
}


#endif // USE_MOMENTUM_FEEDBACK //}

