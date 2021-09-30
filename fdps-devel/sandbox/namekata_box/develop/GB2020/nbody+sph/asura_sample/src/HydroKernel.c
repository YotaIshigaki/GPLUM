#include "config.h"
#include "PlantHydroTree.h"
#include "NeighborSearch.h"
#include "HydroDensity.h"
#include "HydroKernel.h"
#include "KernelFunctions.h"

/*! \file HydroKernel.c
 * \brief Kernel size evaluation function.
 */

static int MaxIterationTimes = 20;
static bool OverMaxIterationTimes = false;
static double LocalKernelMax = 0.e0;

#define  USE_KERNEL_LOCAL_UPDATE_KERNELDENSITY 
#define  ADD_PERTURBATION 

#define KernelFactInc   (1.14) // 1.5 ^ (1.3)
#define KernelFactDec   (0.79) // 0.75 ^ (1.3)
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
static double SmoothedMassConversionFactor;
#endif //USE_SMOOTHED_NEIGHBOR_NUMBER

struct StructHydroKernelExport{
    double    Kernel;  // Kernel size.
    double    Pos[3];  // Position.
    int       Leaf;
    bool      ExtraIterationFlag;
#ifdef USE_DEBUG_MODE
     unsigned long int GlobalID;
#endif // USE_DEBUG_MODE
};

struct StructHydroKernelImport{
    double    SmoothedNumber;    // Mass.
    double    Rho;
    int       Nlist;   // Nlist.
    int       Leaf;
    bool      ExtraIterationFlag;
};

struct StructActiveHydroParticle{
    int Index; // NBCache[Index].
    int Nlist; // Number of neighbors.
    double Rho;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    int ExtraIteration;
    bool LocalUpdateFlags; // Flag for the local update.
    double Rvalue;
    double Lvalue;
}; 

static int NContactedDomains;
static int *ContactedDomainID;

static inline void AllocateContactedDomainID(void){
    ContactedDomainID = malloc(sizeof(int)*MPIGetNumProcs());
    return ;
}

static inline bool CheckLocalExternalDomainsContacted(const int MyID, const int ExtID){

    for(int k=0;k<3;k++){
        if((EdgesForHydro[MyID].PosMax[k] < EdgesForHydro[ExtID].PosMin[k])||
           (EdgesForHydro[MyID].PosMin[k] > EdgesForHydro[ExtID].PosMax[k]))  return false;
    }
    return true;

}

/*
 * This function checkes the number of contacted domains by comparing the local
 * domain edge to the external domains. 
 */
static inline void CheckContactedDomain(void){
    int NProcs = MPIGetNumProcs();
    int MyID = MPIGetMyID();
    NContactedDomains = 0;
    for(int i=0;i<NProcs-1;i++){
        int NodeID = CommunicationTable[i].SendRank;
        assert(MPIGetMyID() != NodeID);
        if(CheckLocalExternalDomainsContacted(MyID,NodeID)){
            ContactedDomainID[NContactedDomains] = i;
            NContactedDomains ++;
        }
    }
    return ;
}


static inline bool __attribute__((always_inline)) OverlapDomainKernel(double Pos[restrict], const double h, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline double __attribute__((always_inline)) DomainDistanceSQR(double Pos[restrict], const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2);
}

static inline int __attribute__((always_inline)) CheckHydroKernelExportFlags(const int Index, const int NProcs, 
        bool HydroKernelExportFlags[][NProcs]){

    if(Pall.Nhydro == 0)
        return 0;

    int ExportNodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    while(CurrentNodeID != RootNodeID){
        if(HydroNode[CurrentNodeID].NumberofActiveLeaves == 0){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else 
            if( SQ(HydroNode[CurrentNodeID].KernelMax) < DomainDistanceSQR(HydroNode[CurrentNodeID].Pos,ExportNodeID) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else 
            if(HydroNode[CurrentNodeID].Children != NONE){
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
            int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
                if(NBCache[leaf].Active){
                    if(HydroKernelExportFlags[leaf][NProcs-1]){
                        if(OverlapDomainKernel(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,ExportNodeID)){
                            HydroKernelExportFlags[leaf][Index] = true;
                            NExport ++;
                        }
                    }
                }
            }
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
    }
    
	return NExport;
}

static inline int __attribute__((always_inline)) CheckHydroKernelExportFlagsCheck(const int Index, const int NProcs, bool HydroKernelExportFlags[][NProcs]){

    if(Pall.Nhydro == 0)
        return 0;

    int ExportNodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            if(HydroKernelExportFlags[leaf][NProcs-1]){
                if(OverlapDomainKernel(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,ExportNodeID)){
                    HydroKernelExportFlags[leaf][Index] = true;
                    NExport ++;
                }
            }
        }
    }
    
	return NExport;
}

/*
 * This function return number of particles which should export to the node ID of [Index].
 */
static inline int __attribute__((always_inline)) CheckHydroKernelExportFlagsSequencialModified(const int NodeIndex, const int NProcs, bool HydroKernelExportFlags[restrict], const int NActives, struct StructActiveHydroParticle ActiveHydroParticle[restrict]){

    if(Pall.Nhydro == 0)
        return 0;

    int ExportNodeID = CommunicationTable[NodeIndex].SendRank;

    // Node By Node Comparison
    double BoxCenter[] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]};
    if(!OverlapDomainKernel(BoxCenter,LocalKernelMax,ExportNodeID)){
        return 0;
    }

    int NExport = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(HydroKernelExportFlags[Offset+NProcs-1]){
            int leaf = ActiveHydroParticle[i].Index;
            if(OverlapDomainKernel(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,ExportNodeID)){
                HydroKernelExportFlags[Offset+NodeIndex] = true;
                NExport ++;
            }
        }
    }

    return NExport;
}

static int GetNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict]){

    int nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,2.e0*Kernel,&nlist,Neighbors);
    }while(CurrentNodeID != RootNodeID);

    return nlist;
}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
static int GetSmoothedNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict], double *SmoothedNumber){

    int nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsSmoothedNumberIterativeApproach(CurrentNodeID,
                Pos,2.e0*Kernel,&nlist,Neighbors,SmoothedNumber);
    }while(CurrentNodeID != RootNodeID);

    return nlist;
}
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

static inline void OverwriteNeighborInfo(
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER  //{
        const double Nlist, 
#else 
        const int Nlist, 
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER  //}
        const int leaf){

    int Index = NBCache[leaf].Leaf;
    Phydro[Index]->Nlist = Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
    Phydro[Index]->SmoothedNumber = Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
    Phydro[Index]->Kernel = 
    Phydro[Index]->KernelPred = NBCache[leaf].Kernel;

    return ;
}


static inline bool __attribute__((always_inline)) CheckNeighborNumberAndUpdateKernelModified_i(const int Index, const int NProcs, bool HydroKernelExportFlags_i[restrict], const int leaf, struct StructActiveHydroParticle *ActiveHydroParticle_i){ 

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#ifdef USE_MAXIMUM_KERNEL_SIZE
#ifdef MAXIMUM_KERNEL_SIZE
    double MaximumKernelSize = Pall.AdaptiveSofteningFactor*MAXIMUM_KERNEL_SIZE*KPC_CGS/Pall.UnitLength;
#else
#error Set MAXIMUM_KERNEL_SIZE
#endif
#endif

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = ActiveHydroParticle_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(NBCache[ActiveHydroParticle_i->Index].Kernel);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = ActiveHydroParticle_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        HydroKernelExportFlags_i[NProcs-1] = false;
        OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }else if((Nlist>=(0.8*NBmin))&&(Nlist<=(2*NBmax))
            &&(ActiveHydroParticle_i->Rvalue>0.e0)&&(ActiveHydroParticle_i->Lvalue>0.e0)){
        if(ActiveHydroParticle_i->Rvalue-ActiveHydroParticle_i->Lvalue < 1.e-6*ActiveHydroParticle_i->Lvalue){
            HydroKernelExportFlags_i[NProcs-1] = false;
            OverwriteNeighborInfo(Nlist,leaf);
            return true;
        }
    }
#ifdef USE_MINIMUM_KERNEL_SIZE
    else if (((NBmax)<Nlist)&&(NBCache[leaf].Kernel<=0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor)){
        NBCache[leaf].Kernel = 0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor;
        HydroKernelExportFlags_i[NProcs-1] = false;
        OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }
#endif
#ifdef USE_MAXIMUM_KERNEL_SIZE
    else if (((NBmin)>Nlist)&&(NBCache[leaf].Kernel>MaximumKernelSize)){
        NBCache[leaf].Kernel = MaximumKernelSize;
        HydroKernelExportFlags_i[NProcs-1] = false;
        OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }
#endif
    if(HydroKernelExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            ActiveHydroParticle_i->Lvalue = fmax(ActiveHydroParticle_i->Lvalue,NBCache[leaf].Kernel);
        } else if(Nlist>NBmax){
            if(ActiveHydroParticle_i->Rvalue > 0.e0){
                ActiveHydroParticle_i->Rvalue = fmin(ActiveHydroParticle_i->Rvalue,NBCache[leaf].Kernel);
            }else{
                ActiveHydroParticle_i->Rvalue = NBCache[leaf].Kernel;
            }
        }

        if((ActiveHydroParticle_i->Lvalue>0.e0)&&(ActiveHydroParticle_i->Rvalue>0.e0)){
            NBCache[leaf].Kernel = cbrt(0.5*(CUBE(ActiveHydroParticle_i->Lvalue)+CUBE(ActiveHydroParticle_i->Rvalue)));
        }else{
            if((ActiveHydroParticle_i->Rvalue == 0.e0)&&(ActiveHydroParticle_i->Lvalue > 0.e0)){
                NBCache[leaf].Kernel *= KernelFactInc;
            }else if((ActiveHydroParticle_i->Rvalue > 0.e0)&&(ActiveHydroParticle_i->Lvalue == 0.e0)){
                NBCache[leaf].Kernel *= KernelFactDec;
            }
        }
    }
    
    return false;
}

static inline bool __attribute__((always_inline)) CheckInLocalDomain(double Pos[], double Kernel, const int MyID){
    for(int k=0;k<3;k++){
        if(Pos[k]+2.e0*Kernel > EdgesForHydro[MyID].PosMax[k]) return false;
        if(Pos[k]-2.e0*Kernel < EdgesForHydro[MyID].PosMin[k]) return false;
    }
    return true;
}


static inline void __attribute__((always_inline)) UpdateKernelLocalModified(const int Index, struct StructActiveHydroParticle ActiveHydroParticle[restrict], const int MyID, const int NProcs, bool HydroKernelExportFlags[restrict]){

    int leaf = ActiveHydroParticle[Index].Index;
    if(CheckNeighborNumberAndUpdateKernelModified_i(Index,NProcs,HydroKernelExportFlags+Index*NProcs,leaf,ActiveHydroParticle+Index) == true)
        return;

    int counter = 0;
    do{
        if(!CheckInLocalDomain(NBCache[leaf].Pos,NBCache[leaf].Kernel,MyID)) return;
        for(int i=0;i<NContactedDomains;i++){
            int NodeID = ContactedDomainID[i];
            if(OverlapDomainKernel(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,NodeID)) return;
        }

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveHydroParticle[Index].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        ActiveHydroParticle[Index].Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            GetSmoothedNumberofNeighbors(NBCache[leaf].Pos,NBCache[leaf].Kernel,Neighbors,
                    &(ActiveHydroParticle[Index].SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            GetNeighborNumbers(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        counter ++;
        if(counter > 10) return ;
    }while(CheckNeighborNumberAndUpdateKernelModified_i(Index,NProcs,HydroKernelExportFlags+Index*NProcs,leaf,ActiveHydroParticle+Index) == false);

    return;
}


static inline int __attribute__((always_inline)) CheckNeighborNumberAndUpdateKernelModified(const int NActives, const int NProcs, bool HydroKernelExportFlags[restrict], struct StructActiveHydroParticle ActiveHydroParticle[restrict]){ 

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#ifdef USE_MAXIMUM_KERNEL_SIZE
#ifdef MAXIMUM_KERNEL_SIZE
    double MaximumKernelSize = Pall.AdaptiveSofteningFactor*MAXIMUM_KERNEL_SIZE*KPC_CGS/Pall.UnitLength;
#else
#error Set MAXIMUM_KERNEL_SIZE
#endif
#endif

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(HydroKernelExportFlags[Offset+NProcs-1]){ 
#ifdef USE_KERNEL_LOCAL_UPDATE
        if(ActiveHydroParticle[i].LocalUpdateFlags == false){ 
#endif

        int leaf = ActiveHydroParticle[i].Index;

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            double Nlist = ActiveHydroParticle[i].SmoothedNumber*
                SmoothedMassConversionFactor*CUBE(NBCache[leaf].Kernel);
#else
            int Nlist = ActiveHydroParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}

            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                HydroKernelExportFlags[Offset+NProcs-1] = false;
                OverwriteNeighborInfo(Nlist,leaf);
            }else if((Nlist>=(0.8*NBmin))&&(Nlist<=(2*NBmax))
                    &&(ActiveHydroParticle[i].Rvalue>0.e0)&&(ActiveHydroParticle[i].Lvalue>0.e0)){
                if(ActiveHydroParticle[i].Rvalue-ActiveHydroParticle[i].Lvalue < 1.e-6*ActiveHydroParticle[i].Lvalue){
                    HydroKernelExportFlags[Offset+NProcs-1] = false;
                    OverwriteNeighborInfo(Nlist,leaf);
                }
            }
#ifdef USE_MINIMUM_KERNEL_SIZE
            else if (((NBmax)<Nlist)&&(NBCache[leaf].Kernel<=0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor)){
                NBCache[leaf].Kernel = 0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor;
                HydroKernelExportFlags[Offset+NProcs-1] = false;
                OverwriteNeighborInfo(Nlist,leaf);
            }
#endif
#ifdef USE_MAXIMUM_KERNEL_SIZE
            else if (((NBmin)>Nlist)&&(NBCache[leaf].Kernel>MaximumKernelSize)){
                NBCache[leaf].Kernel = MaximumKernelSize;
                HydroKernelExportFlags[Offset+NProcs-1] = false;
                OverwriteNeighborInfo(Nlist,leaf);
            }
#endif
            if(HydroKernelExportFlags[Offset+NProcs-1]){
                if(Nlist<NBmin){
                    ActiveHydroParticle[i].Lvalue = fmax(ActiveHydroParticle[i].Lvalue,NBCache[leaf].Kernel);
                } else if(Nlist>NBmax){
                    if(ActiveHydroParticle[i].Rvalue > 0.e0){
                        ActiveHydroParticle[i].Rvalue = fmin(ActiveHydroParticle[i].Rvalue,NBCache[leaf].Kernel);
                    }else{
                        ActiveHydroParticle[i].Rvalue = NBCache[leaf].Kernel;
                    }
                }

                if((ActiveHydroParticle[i].Lvalue>0.e0)&&(ActiveHydroParticle[i].Rvalue>0.e0)){
                    NBCache[leaf].Kernel = cbrt(0.5*(CUBE(ActiveHydroParticle[i].Lvalue)+CUBE(ActiveHydroParticle[i].Rvalue)));
                }else{
                    if((ActiveHydroParticle[i].Rvalue == 0.e0)&&(ActiveHydroParticle[i].Lvalue > 0.e0)){
                        NBCache[leaf].Kernel *= KernelFactInc;
                    }else if((ActiveHydroParticle[i].Rvalue > 0.e0)&&(ActiveHydroParticle[i].Lvalue == 0.e0)){
                        NBCache[leaf].Kernel *= KernelFactDec;
                    }
                }
                NLocalActiveLeaves ++;
            }
#ifdef USE_KERNEL_LOCAL_UPDATE
        } else {
            NLocalActiveLeaves ++;
        }
#endif
        }
    }
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(HydroKernelExportFlags[Offset+NProcs-1]){ 
            int leaf = ActiveHydroParticle[i].Index;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            double Nlist = ActiveHydroParticle[i].SmoothedNumber*
                SmoothedMassConversionFactor*CUBE(NBCache[leaf].Kernel);
#else
            int Nlist = ActiveHydroParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}
            OverwriteNeighborInfo(Nlist,leaf);
        }
    }

    return NLocalActiveLeaves;
}

#define _ResetTiming_ 1
static void ResetKernelSizeModified(const int NActives, const int Niteration, const int NProcs, 
            bool HydroKernelExportFlags[restrict], struct StructActiveHydroParticle ActiveHydroParticle[restrict]){ 

    if((Niteration+1)%(_ResetTiming_*MaxIterationTimes) == 0){
        for(int i=0;i<NActives;i++){
            int Offset = i*NProcs;
            if(HydroKernelExportFlags[Offset+NProcs-1]){ 
                int leaf = ActiveHydroParticle[i].Index;
                ActiveHydroParticle[i].Rvalue = ActiveHydroParticle[i].Lvalue = 0.e0;
            }
        }
    }
    return ;
}

static bool first = true;
void CalcKernel(void){

    double TimingResultThisRoutine = GetElapsedTime();

#ifdef EVALUATE_KERNEL_BY_ITERATION
    if(first){
        AllocateContactedDomainID();
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
        SmoothedMassConversionFactor = (4.0*M_PI/3.0)*8.0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
        first = false;
    }

    OverMaxIterationTimes = false;

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int Neighbors[MaxNeighborSize];
    MPI_Status  mpi_status;

    static int HydroKernelExportFlagsMaxAllocated = 0;

    static bool *HydroKernelExportFlags;
    static struct StructActiveHydroParticle *ActiveHydroParticle;

    if(HydroKernelExportFlagsMaxAllocated < MAX(Pall.Nhydro,NAdditionUnit)){
        if(HydroKernelExportFlagsMaxAllocated > 0){
            free(HydroKernelExportFlags);
            free(ActiveHydroParticle);
        }
        HydroKernelExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        HydroKernelExportFlags = malloc(sizeof(bool)*HydroKernelExportFlagsMaxAllocated*NProcs);
        ActiveHydroParticle = malloc(sizeof(struct StructActiveHydroParticle)*HydroKernelExportFlagsMaxAllocated);
    }

    int NActives = 0;
    int RootNodeID = 0; 
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header + i; 
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            HydroKernelExportFlags[NActives*NProcs+NProcs-1] = true;
            ActiveHydroParticle[NActives].Index = leaf;
            ActiveHydroParticle[NActives].Nlist = Phydro[NBCache[leaf].Leaf]->Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
            ActiveHydroParticle[NActives].SmoothedNumber = Phydro[NBCache[leaf].Leaf]->SmoothedNumber;
#endif //ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
            ActiveHydroParticle[NActives].Rvalue =
            ActiveHydroParticle[NActives].Lvalue = 0.e0;

            NActives ++;
        }
    }

    int BitMask = 0x01; 
    int NExportThisTime[NProcs-1];
    int NImportThisTime[NProcs-1];
    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];

    struct StructHydroKernelExport *HydroKernelExportSend[NProcs-1];
    struct StructHydroKernelExport *HydroKernelExportRecv = NULL;
    struct StructHydroKernelImport *HydroKernelImportSend = NULL;
    struct StructHydroKernelImport *HydroKernelImportRecv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Kernel Iteration");
#endif // PRINT_LOG_KERNEL_ITERATION

    int NActiveLeaves;
    MPI_Allreduce(&NActives,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);


    CheckContactedDomain();
    double BoxCenter[] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]};

    int Niteration = 0;
    do{
#ifdef PRINT_LOG_KERNEL_ITERATION
        if(MPIGetMyID()==MPI_ROOT_RANK)
            fprintf(stderr,":[%d] = %d ",Niteration,NActiveLeaves);
#endif // PRINT_LOG_KERNEL_ITERATION

        LocalKernelMax = 0.e0;  
        for(int i=0;i<NActives;i++){ 
            int leaf = ActiveHydroParticle[i].Index;
            int Offset = i*NProcs;
            if(HydroKernelExportFlags[Offset+NProcs-1]){
                ActiveHydroParticle[i].Nlist = 0;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveHydroParticle[i].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
                for(int k=0;k<NProcs-1;k++)
                    HydroKernelExportFlags[Offset+k] = false;
                LocalKernelMax = fmax(LocalKernelMax,
                        DISTANCE(BoxCenter,NBCache[leaf].Pos)+2.0*NBCache[leaf].Kernel);
            }
#ifdef TASK_TEST_HYDRO_QUANTITIES //{
            Phydro[NBCache[leaf].Leaf]->Tag = 0;
#endif // TASK_TEST_HYDRO_QUANTITIES //}
        }

        int NExportMaxThisTime = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime[i] = CheckHydroKernelExportFlagsSequencialModified(i,
                    NProcs,HydroKernelExportFlags,NActives,ActiveHydroParticle);

            CheckSizeofBufferExportSendIndex(NExportThisTime[i],
                    sizeof(struct StructHydroKernelExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],
                    sizeof(struct StructHydroKernelImport),i);
            HydroKernelExportSend[i] = BufferExportSend[i];
            HydroKernelImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            if(NExportThisTime[i] > 0){
            for(int k=0;k<NActives;k++){
                int Offset = k*NProcs;
                if(HydroKernelExportFlags[Offset+NProcs-1]){ 
                    if(HydroKernelExportFlags[Offset+i]&BitMask){ 
                        int leaf = ActiveHydroParticle[k].Index;
                        HydroKernelExportSend[i][NExport].Pos[0] = NBCache[leaf].Pos[0];
                        HydroKernelExportSend[i][NExport].Pos[1] = NBCache[leaf].Pos[1];
                        HydroKernelExportSend[i][NExport].Pos[2] = NBCache[leaf].Pos[2];
                        HydroKernelExportSend[i][NExport].Kernel = NBCache[leaf].Kernel;
                        HydroKernelExportSend[i][NExport].Leaf = k;
#ifdef USE_DEBUG_MODE //{
                        HydroKernelExportSend[i][NExport].GlobalID = PhydroBody(NBCache[leaf].Leaf)->GlobalID;
#endif // USE_DEBUG_MODE //}
                        NExport ++;
                    }
                }
            }
            }
            NExportThisTime[i] = NExport;
            NExportMaxThisTime = MAX(NExportMaxThisTime,NExport);
        }

        int NImportThisTime2[NProcs];
        int NExportThisTime2[NProcs];
        NImportThisTime2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
        }
        MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
        int NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
            NImport += NImportThisTime[i];
        }
        int NImportAll = NImport;

        CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHydroKernelExport));
        CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructHydroKernelImport));
        HydroKernelExportRecv = BufferExportRecv;
        HydroKernelImportSend = BufferImportSend; 

        NImport = 0;

        int counter_send = 0;
        int counter_recv = 0;

        int SendFlag,RecvFlag;
        for(int i=0;i<NProcs-1;i++){
            if(NExportThisTime[i]>0){
                MPI_Isend(HydroKernelExportSend[i],
                    NExportThisTime[i]*sizeof(struct StructHydroKernelExport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NImportThisTime[i]>0){
                MPI_Irecv(HydroKernelExportRecv+NImport,
                    NImportThisTime[i]*sizeof(struct StructHydroKernelExport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
                NImport += NImportThisTime[i];
            }
        }

        double TimeNBS = GetElapsedTime();
        for(int i=0;i<NActives;i++){  // Check local
            int Offset = i*NProcs;
            if(HydroKernelExportFlags[Offset+NProcs-1]){
                int leaf = ActiveHydroParticle[i].Index;
                ActiveHydroParticle[i].Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                    GetSmoothedNumberofNeighbors(NBCache[leaf].Pos,NBCache[leaf].Kernel,
                            Neighbors,&(ActiveHydroParticle[i].SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
                    GetNeighborNumbers(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef USE_KERNEL_LOCAL_UPDATE
                /// Insert Local Update Routine here.
                ActiveHydroParticle[i].LocalUpdateFlags = false;
                int IsLocal = 0;
                for(int k=0;k<NProcs-1;k++){
                    if(HydroKernelExportFlags[Offset+k]&BitMask){
                        IsLocal ++;
                    }
                }

                if(IsLocal == 0){
                    ActiveHydroParticle[i].LocalUpdateFlags = true;
                    UpdateKernelLocalModified(i,ActiveHydroParticle,
                            MyID,NProcs,HydroKernelExportFlags);
                }
#endif // LOCAL_UPDATE
            }
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

        double TimeComm = GetElapsedTime();
        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
        TimingResults.HydroKernelCommThisStep += GetElapsedTime()-TimeComm;


        TimeNBS = GetElapsedTime();
        for(int i=0;i<NImportAll;i++){
            HydroKernelImportSend[i].SmoothedNumber = 0.e0;
            HydroKernelImportSend[i].Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                GetSmoothedNumberofNeighbors(HydroKernelExportRecv[i].Pos,
                    HydroKernelExportRecv[i].Kernel,Neighbors,&(HydroKernelImportSend[i].SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
                GetNeighborNumbers(HydroKernelExportRecv[i].Pos,2.e0*HydroKernelExportRecv[i].Kernel);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
            HydroKernelImportSend[i].Leaf = HydroKernelExportRecv[i].Leaf;
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

        NImportAll = 0;
        int NImportAllNew = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew[i] = 0;
            for(int k=0;k<NImportThisTime[i];k++){
                if(HydroKernelImportSend[NImportAll].Nlist > 0){
                    HydroKernelImportSend[NImportAllNew] = HydroKernelImportSend[NImportAll];
                    NImportThisTimeNew[i] ++;
                    NImportAllNew ++;
                }
                NImportAll ++;
            }
        }

        int NImportThisTimeNew2[NProcs];
        int NExportThisTimeNew2[NProcs];
        NImportThisTimeNew2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew2[CommunicationTable[i].SendRank] = NImportThisTimeNew[i];
        }
        MPI_Alltoall(NImportThisTimeNew2,1,MPI_INT,NExportThisTimeNew2,1,MPI_INT,MPI_COMM_WORLD);
        for(int i=0;i<NProcs-1;i++){
            NExportThisTimeNew[i] = NExportThisTimeNew2[CommunicationTable[i].RecvRank];
        }

        NImport = 0;
        counter_send = counter_recv = 0;
        for(int i=0;i<NProcs-1;i++){
            if(NImportThisTimeNew[i]>0){
                MPI_Isend(HydroKernelImportSend+NImport,
                    NImportThisTimeNew[i]*sizeof(struct StructHydroKernelImport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NExportThisTimeNew[i]>0){
                MPI_Irecv(HydroKernelImportRecv[i],
                    NExportThisTimeNew[i]*sizeof(struct StructHydroKernelImport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
            }
            NImport += NImportThisTimeNew[i];
        }
        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);


        for(int i=0;i<NProcs-1;i++){
            for(int k=0;k<NExportThisTimeNew[i];k++){ 
                int leaf = HydroKernelImportRecv[i][k].Leaf;
                ActiveHydroParticle[leaf].Nlist += HydroKernelImportRecv[i][k].Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
                ActiveHydroParticle[leaf].SmoothedNumber += HydroKernelImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
#ifdef TASK_TEST_HYDRO_QUANTITIES //{
                if(HydroKernelImportRecv[i][k].Nlist > 0)
                    Phydro[NBCache[leaf].Leaf]->Tag = -1;
#endif // TASK_TEST_HYDRO_QUANTITIES //}
            }
        }

#ifdef ADD_PERTURBATION //{
        // assert(Niteration < 1000);
        if(Niteration > 10*MaxIterationTimes) break;
        if(Niteration > MaxIterationTimes)
            OverMaxIterationTimes = true;
        ResetKernelSizeModified(NActives,Niteration,NProcs,HydroKernelExportFlags,ActiveHydroParticle);
#endif // ADD_PERTURBATION //}
        int NLocalActiveLeaves = CheckNeighborNumberAndUpdateKernelModified(NActives,
                NProcs,HydroKernelExportFlags,ActiveHydroParticle);
        MPI_Allreduce(&NLocalActiveLeaves,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        Niteration ++;
#ifndef ADD_PERTURBATION //{
        if(Niteration > MaxIterationTimes)
            break;
#endif // ADD_PERTURBATION //}
    } while (0<NActiveLeaves);

#else // EVALUATE_KERNEL_BY_ITERATION
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->Kernel = KERNEL_FACTOR*
                pow(Phydro[i]->Mass/Phydro[i]->Rho,1.0/((double)DIMENSION));
            Phydro[i]->Rho = 0.e0;
        }
    }
#endif // EVALUATE_KERNEL_BY_ITERATION

    PlantHydroTreeKernelMaxUpdate();
#ifdef EVALUATE_KERNEL_BY_ITERATION
#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"\n");
#else // PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"%d iterations for kernel determination.\n",Niteration);
#endif // PRINT_LOG_KERNEL_ITERATION
#endif // EVALUATE_KERNEL_BY_ITERATION

    TimingResults.HydroKernelThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

void CountNeighborNumber(void){

    int Nmin = Phydro[0]->Nlist;
    int Nmax = Nmin;
    int Nmean = Nmin;
    for(int i=1;i<Pall.Nhydro;i++){
        Nmin = MIN(Nmin,Phydro[i]->Nlist);
        Nmax = MAX(Nmax,Phydro[i]->Nlist);
        Nmean += Phydro[i]->Nlist;
    }
    int GlobalMean,GlobalMin,GlobalMax;
    MPI_Allreduce(&Nmean,&GlobalMean,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&Nmin,&GlobalMin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&Nmax,&GlobalMax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Number of Neighbors : Nmean Nmax Nmin = %ld %d %d %d\n",
            GlobalMean/Pall.Nhydro_t,Nmax,Nmin,Nmean);
    }
    return;
}

void CountDirectNeighborNumber(void){

    int Nmin = Phydro[0]->Nlist;
    int Nmax = Nmin;
    int Nmean = Nmin;
    for(int i=1;i<Pall.Nhydro;i++){
        Nmin = MIN(Nmin,Phydro[i]->Nlist);
        Nmax = MAX(Nmax,Phydro[i]->Nlist);
        Nmean += Phydro[i]->Nlist;
    }
    fprintf(stderr,"Direct Neighbor Number : Nmean Nmax Nmin = %ld %d %d\n",Nmean/Pall.Nhydro,Nmax,Nmin);

    return;
}

#if 0
#define ENTROPY_FORM_INITIAL_DENSITY_EVALUATION (1)

static double __attribute__((always_inline)) CalcRho_i(double Pos[restrict], const double Kernel, const int Nlist, int Neighbors[restrict]){
    double InvKerneli = 1.0/Kernel;
    double Rho = 0.0;
    for(int i=0;i<Nlist;i++){
        int leaf = Neighbors[i];
        double r = DISTANCE(Pos,Phydro[leaf]->PosP);
        double w = SPHKernel(r,InvKerneli);
        Rho += Phydro[leaf]->Mass*w;
    }
    return Rho;
}

static inline bool __attribute__((always_inline)) CheckNeighborNumberAndUpdateKernelDensity_i(const int Index, const int NProcs, bool HydroKernelExportFlags_i[restrict], const int leaf, struct StructActiveHydroParticle *ActiveHydroParticle_i, struct StructBinarySearch *BinarySearch_i, const int Neighbors[restrict]){ 

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#ifdef USE_MAXIMUM_KERNEL_SIZE
#ifdef MAXIMUM_KERNEL_SIZE
    double MaximumKernelSize = Pall.AdaptiveSofteningFactor*MAXIMUM_KERNEL_SIZE*KPC_CGS/Pall.UnitLength;
#else
#error Set MAXIMUM_KERNEL_SIZE
#endif
#endif

    if(ActiveHydroParticle_i->ExtraIteration == 0){
        int index = NBCache[leaf].Leaf;
        Phydro[index]->Rho = Phydro[index]->RhoPred = ActiveHydroParticle_i->Rho;
        HydroKernelExportFlags_i[NProcs-1] = false;
        return true;
    } else {

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        double Nlist = ActiveHydroParticle_i.SmoothedNumber*
            SmoothedMassConversionFactor*CUBE(NBCache[ActiveHydroParticle_i->Index].Kernel);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
        int Nlist = ActiveHydroParticle_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

        if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
            OverwriteNeighborInfo(Nlist,leaf);
            ActiveHydroParticle_i->ExtraIteration --;
        }else if((Nlist<=(NBmax))&&(BinarySearch_i->Rvalue>0.e0)&&(BinarySearch_i->Lvalue>0.e0)){
            if(BinarySearch_i->Rvalue-BinarySearch_i->Lvalue < 1.e-3*BinarySearch_i->Lvalue){
                OverwriteNeighborInfo(Nlist,leaf);
                ActiveHydroParticle_i->ExtraIteration --;
            }
        }
#ifdef USE_MINIMUM_KERNEL_SIZE //{
        else if (((NBmax)<Nlist)&&(NBCache[leaf].Kernel<=0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor)){
            NBCache[leaf].Kernel = 0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor;
            OverwriteNeighborInfo(Nlist,leaf);
            ActiveHydroParticle_i->ExtraIteration --;
        }
#endif // USE_MINIMUM_KERNEL_SIZE //}
#ifdef USE_MAXIMUM_KERNEL_SIZE //{
        else if (((NBmin)>Nlist)&&(NBCache[leaf].Kernel>MaximumKernelSize)){
            NBCache[leaf].Kernel = MaximumKernelSize;
            OverwriteNeighborInfo(Nlist,leaf);
            ActiveHydroParticle_i->ExtraIteration --;
        }
#endif // USE_MAXIMUM_KERNEL_SIZE //}

        if((ActiveHydroParticle_i->ExtraIteration > 0)&&(HydroKernelExportFlags_i[NProcs-1] == true)){
            if(Nlist<NBmin){
                BinarySearch_i->Lvalue = fmax(BinarySearch_i->Lvalue,NBCache[leaf].Kernel);
            } else if(Nlist>NBmax){
                if(BinarySearch_i->Rvalue > 0.e0){
                    BinarySearch_i->Rvalue = fmin(BinarySearch_i->Rvalue,NBCache[leaf].Kernel);
                }else{
                    BinarySearch_i->Rvalue = NBCache[leaf].Kernel;
                }
            }

            if((BinarySearch_i->Lvalue>0.e0)&&(BinarySearch_i->Rvalue>0.e0)){
                NBCache[leaf].Kernel = cbrt(0.5*(CUBE(BinarySearch_i->Lvalue)+CUBE(BinarySearch_i->Rvalue)));
            }else{
                if((BinarySearch_i->Rvalue == 0.e0)&&(BinarySearch_i->Lvalue > 0.e0)){
                    NBCache[leaf].Kernel *= KernelFactInc;
                }else if((BinarySearch_i->Rvalue > 0.e0)&&(BinarySearch_i->Lvalue == 0.e0)){
                    NBCache[leaf].Kernel *= KernelFactDec;
                }
            }
        }
    }

    return false;
}

static inline void UpdateKernelDensityLocal(const int Index, struct StructActiveHydroParticle ActiveHydroParticle[restrict], int Neighbors[restrict], const int MyID, const int NProcs, bool HydroKernelExportFlags[restrict][NProcs], struct StructBinarySearch BinarySearch[restrict]){

    int leaf = ActiveHydroParticle[Index].Index;
    if(CheckNeighborNumberAndUpdateKernelDensity_i(Index,NProcs,HydroKernelExportFlags[leaf],leaf,ActiveHydroParticle+Index,BinarySearch+Index,Neighbors) == true)
        return;

    do{
        if(!CheckInLocalDomain(NBCache[leaf].Pos,NBCache[leaf].Kernel,MyID)) return;
        for(int i=0;i<NContactedDomains;i++){
            int NodeID = ContactedDomainID[i];
            if(OverlapDomainKernel(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,NodeID)) return;
        }

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveHydroParticle[Index].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        ActiveHydroParticle[Index].Nlist = 
        //ActiveNeighborList[Index] = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            GetSmoothedNumberofNeighbors(NBCache[leaf].Pos,NBCache[leaf].Kernel,Neighbors,
                    &(ActiveHydroParticle[Index].SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            GetNumberofNeighbors(NBCache[leaf].Pos,NBCache[leaf].Kernel,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

        if(ActiveHydroParticle[Index].ExtraIteration == 0){
            ActiveHydroParticle[Index].Rho = CalcRho_i(NBCache[leaf].Pos,NBCache[leaf].Kernel,
                    ActiveHydroParticle[Index].Nlist,Neighbors);
        }
    }while(CheckNeighborNumberAndUpdateKernelDensity_i(Index,NProcs,HydroKernelExportFlags[leaf],leaf,ActiveHydroParticle+Index,BinarySearch+Index,Neighbors) == false);

    return;
}


static inline int __attribute__((always_inline)) CheckNeighborNumberAndUpdateKernelDensity(const int NActives, const int NProcs, bool HydroKernelExportFlags[restrict][NProcs], struct StructActiveHydroParticle ActiveHydroParticle[restrict], struct StructBinarySearch BinarySearch[restrict]){ 

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#ifdef USE_MAXIMUM_KERNEL_SIZE
#ifdef MAXIMUM_KERNEL_SIZE
    double MaximumKernelSize = Pall.AdaptiveSofteningFactor*MAXIMUM_KERNEL_SIZE*KPC_CGS/Pall.UnitLength;
#else
#error Set MAXIMUM_KERNEL_SIZE
#endif
#endif

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int leaf = ActiveHydroParticle[i].Index;
        if(HydroKernelExportFlags[leaf][NProcs-1]){ 
#ifdef USE_KERNEL_LOCAL_UPDATE_KERNELDENSITY //{
        if(ActiveHydroParticle[i].LocalUpdateFlags == false){ 
#endif //}

        if(ActiveHydroParticle[i].ExtraIteration == 0){
            Phydro[NBCache[leaf].Leaf]->Rho = Phydro[NBCache[leaf].Leaf]->RhoPred = ActiveHydroParticle[i].Rho;
            HydroKernelExportFlags[leaf][NProcs-1] = false;
        } else {
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            double Nlist = ActiveHydroParticle[i].SmoothedNumber*
                SmoothedMassConversionFactor*CUBE(NBCache[leaf].Kernel);
#else
            int Nlist = ActiveHydroParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}
            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                OverwriteNeighborInfo(Nlist,leaf);
                ActiveHydroParticle[i].ExtraIteration --;
            }else if((Nlist<=(NBmax))&&(BinarySearch[i].Rvalue>0.e0)&&(BinarySearch[i].Lvalue>0.e0)){
                if(BinarySearch[i].Rvalue-BinarySearch[i].Lvalue < 1.e-3*BinarySearch[i].Lvalue){
                    OverwriteNeighborInfo(Nlist,leaf);
                    ActiveHydroParticle[i].ExtraIteration --;
                }
            }
#ifdef USE_MINIMUM_KERNEL_SIZE //{
            else if (((NBmax)<Nlist)&&(NBCache[leaf].Kernel<=0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor)){
                NBCache[leaf].Kernel = 0.5*PhydroBody(NBCache[leaf].Leaf)->Eps*Pall.AdaptiveSofteningFactor;
                OverwriteNeighborInfo(Nlist,leaf);
                ActiveHydroParticle[i].ExtraIteration --;
            }
#endif // USE_MINIMUM_KERNEL_SIZE //}
#ifdef USE_MAXIMUM_KERNEL_SIZE //{
            else if (((NBmin)>Nlist)&&(NBCache[leaf].Kernel>MaximumKernelSize)){
                NBCache[leaf].Kernel = MaximumKernelSize;
                OverwriteNeighborInfo(Nlist,leaf);
                ActiveHydroParticle[i].ExtraIteration --;
            }
#endif // USE_MAXIMUM_KERNEL_SIZE //}

            if((ActiveHydroParticle[i].ExtraIteration != 0)&&(HydroKernelExportFlags[leaf][NProcs-1])){
                if(Nlist<NBmin){
                    BinarySearch[i].Lvalue = fmax(BinarySearch[i].Lvalue,NBCache[leaf].Kernel);
                } else if(Nlist>NBmax){
                    if(BinarySearch[i].Rvalue > 0.e0){
                        BinarySearch[i].Rvalue = fmin(BinarySearch[i].Rvalue,NBCache[leaf].Kernel);
                    }else{
                        BinarySearch[i].Rvalue = NBCache[leaf].Kernel;
                    }
                }

                if((BinarySearch[i].Lvalue>0.e0)&&(BinarySearch[i].Rvalue>0.e0)){
                    NBCache[leaf].Kernel = cbrt(0.5*(CUBE(BinarySearch[i].Lvalue)+CUBE(BinarySearch[i].Rvalue)));
                }else{
                    if((BinarySearch[i].Rvalue == 0.e0)&&(BinarySearch[i].Lvalue > 0.e0)){
                        NBCache[leaf].Kernel *= KernelFactInc;
                    }else if((BinarySearch[i].Rvalue > 0.e0)&&(BinarySearch[i].Lvalue == 0.e0)){
                        NBCache[leaf].Kernel *= KernelFactDec;
                    }
                }

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                double Nlist = ActiveHydroParticle[i].SmoothedNumber*
                    SmoothedMassConversionFactor*CUBE(NBCache[leaf].Kernel);
#else
                int Nlist = ActiveHydroParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}
                OverwriteNeighborInfo(Nlist,leaf);

                NLocalActiveLeaves ++;
            }
        }
#ifdef USE_KERNEL_LOCAL_UPDATE_KERNELDENSITY //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif //}
        }
    }

#if 0
    for(int i=0;i<NActives;i++){
        if(ActiveHydroParticle[i].ExtraIteration != 0){
            int leaf = ActiveHydroParticle[i].Index;
            if(HydroKernelExportFlags[leaf][NProcs-1]){ 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                double Nlist = ActiveHydroParticle[i].SmoothedNumber*
                    SmoothedMassConversionFactor*CUBE(NBCache[leaf].Kernel);
#else
                int Nlist = ActiveHydroParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}
                OverwriteNeighborInfo(Nlist,leaf);
            }
        }
    }
#endif

    return NLocalActiveLeaves;
}



static bool FirstCall = true;
void CalcKernelDensity(void){

    double TimingResultThisRoutine = GetElapsedTime();

#ifdef EVALUATE_KERNEL_BY_ITERATION
    if(FirstCall){
        AllocateContactedDomainID();
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
        SmoothedMassConversionFactor = (4.0*M_PI/3.0)*8.0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
        FirstCall = false;
    }

    int MaxIterationTimes;
    if(Pall.TStepTotal == 0){
        MaxIterationTimes = 100;
    } else {
        MaxIterationTimes = 20;
    }

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int Neighbors[MaxNeighborSize];
    MPI_Status  mpi_status;

    static int HydroKernelExportFlagsMaxAllocated = 0;


    static bool (*HydroKernelExportFlags)[NProcs];
    static struct StructBinarySearch *BinarySearch;
    static struct StructActiveHydroParticle *ActiveHydroParticle;

    if(HydroKernelExportFlagsMaxAllocated < MAX(Pall.Nhydro,NAdditionUnit)){
        if(HydroKernelExportFlagsMaxAllocated > 0){
            free(HydroKernelExportFlags);
            free(BinarySearch);
            free(ActiveHydroParticle);
        }
        HydroKernelExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        HydroKernelExportFlags = malloc(sizeof(bool)*HydroKernelExportFlagsMaxAllocated*NProcs);
        BinarySearch = malloc(sizeof(struct StructBinarySearch)*HydroKernelExportFlagsMaxAllocated);
        ActiveHydroParticle = malloc(sizeof(struct StructActiveHydroParticle)*HydroKernelExportFlagsMaxAllocated);
    }

    /////////////////////////////////////////////////////////////////////

    int NActives = 0;
    int RootNodeID = 0; 
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header + k; 
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){

            HydroKernelExportFlags[leaf][NProcs-1] = true;
            ActiveHydroParticle[NActives].Index = leaf;
            ActiveHydroParticle[NActives].Nlist = Phydro[NBCache[leaf].Leaf]->Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            ActiveHydroParticle[NActives].SmoothedNumber = Phydro[NBCache[leaf].Leaf]->SmoothedNumber;
#endif //ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //}
            ActiveHydroParticle[NActives].ExtraIteration = ENTROPY_FORM_INITIAL_DENSITY_EVALUATION;

            NActives ++;
        }
    }
    for(int i=0;i<NActives;i++)
        BinarySearch[i].Rvalue = BinarySearch[i].Lvalue = 0.e0;

    int BitMask = 0x01; 
    int NExportThisTime[NProcs-1];
    int NImportThisTime[NProcs-1];
    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];

    struct StructHydroKernelExport *HydroKernelExportSend[NProcs-1];
    struct StructHydroKernelExport *HydroKernelExportRecv = NULL;
    struct StructHydroKernelImport *HydroKernelImportSend = NULL;
    struct StructHydroKernelImport *HydroKernelImportRecv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Kernel Iteration");
#endif // PRINT_LOG_KERNEL_ITERATION

    int NActiveLeaves;
    MPI_Allreduce(&NActives,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);


    CheckContactedDomain();

    int Niteration = 0;
    do{
#ifdef PRINT_LOG_KERNEL_ITERATION
        if(MPIGetMyID()==MPI_ROOT_RANK)
            fprintf(stderr,":[%d] = %d ",Niteration,NActiveLeaves);
#endif // PRINT_LOG_KERNEL_ITERATION
        for(int i=0;i<NActives;i++){ 
            int leaf = ActiveHydroParticle[i].Index;
            if(HydroKernelExportFlags[leaf][NProcs-1]){
                ActiveHydroParticle[i].Nlist = 0;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveHydroParticle[i].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
                ActiveHydroParticle[i].Rho = 0.0;
                for(int k=0;k<NProcs-1;k++)
                    HydroKernelExportFlags[leaf][k] = 0;
            }
#ifdef TASK_TEST_HYDRO_QUANTITIES //{
            Phydro[NBCache[leaf].Leaf]->Tag = 0;
#endif // TASK_TEST_HYDRO_QUANTITIES //}
        }
        int NExportMaxThisTime = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime[i] = CheckHydroKernelExportFlagsSequencial(i,
                    NProcs,HydroKernelExportFlags,NActives,ActiveHydroParticle);

            CheckSizeofBufferExportSendIndex(NExportThisTime[i],
                    sizeof(struct StructHydroKernelExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],
                    sizeof(struct StructHydroKernelImport),i);
            HydroKernelExportSend[i] = BufferExportSend[i];
            HydroKernelImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            for(int k=0;k<NActives;k++){
                int leaf = ActiveHydroParticle[k].Index;
                if(HydroKernelExportFlags[leaf][NProcs-1]){ 
                    if(HydroKernelExportFlags[leaf][i]&BitMask){ 
                        HydroKernelExportSend[i][NExport].Pos[0] = NBCache[leaf].Pos[0];
                        HydroKernelExportSend[i][NExport].Pos[1] = NBCache[leaf].Pos[1];
                        HydroKernelExportSend[i][NExport].Pos[2] = NBCache[leaf].Pos[2];
                        HydroKernelExportSend[i][NExport].Kernel = NBCache[leaf].Kernel;
                        if(ActiveHydroParticle[k].ExtraIteration == 0){
                            HydroKernelExportSend[i][NExport].ExtraIterationFlag = true;
                        } else {
                            HydroKernelExportSend[i][NExport].ExtraIterationFlag = false;
                        }
                        HydroKernelExportSend[i][NExport].Leaf = k;
                        NExport ++;
                    }
                }
            }
            NExportThisTime[i] = NExport;
            NExportMaxThisTime = MAX(NExportMaxThisTime,NExport);
        }

        int NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,
                    CommunicationTable[i].SendRank,TAG_SPH_KERNEL_PRECOMM,
                NImportThisTime+i,1,MPI_INT,
                    CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_PRECOMM,
                        MPI_COMM_WORLD,&mpi_status);
            NImport += NImportThisTime[i];
        }
        int NImportAll = NImport;

        CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHydroKernelExport));
        CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructHydroKernelImport));
        HydroKernelExportRecv = BufferExportRecv;
        HydroKernelImportSend = BufferImportSend; 

        NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            MPI_Isend(HydroKernelExportSend[i],
                NExportThisTime[i]*sizeof(struct StructHydroKernelExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+i);
            MPI_Irecv(HydroKernelExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructHydroKernelExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+i);
            NImport += NImportThisTime[i];
        }

        double TimeNBS = GetElapsedTime();
        for(int i=0;i<NActives;i++){  // Check local
            int leaf = ActiveHydroParticle[i].Index;
            if(HydroKernelExportFlags[leaf][NProcs-1]){
                ActiveHydroParticle[i].Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                    GetSmoothedNumberofNeighbors(NBCache[leaf].Pos,NBCache[leaf].Kernel,
                            Neighbors,&(ActiveHydroParticle[i].SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
                    GetNumberofNeighbors(NBCache[leaf].Pos,NBCache[leaf].Kernel,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
                if(ActiveHydroParticle[i].ExtraIteration == 0){
                    ActiveHydroParticle[i].Rho = CalcRho_i(NBCache[leaf].Pos,NBCache[leaf].Kernel,
                            ActiveHydroParticle[i].Nlist,Neighbors);
                }

#ifdef USE_KERNEL_LOCAL_UPDATE_KERNELDENSITY
                /// Insert Local Update Routine here.
                ActiveHydroParticle[i].LocalUpdateFlags = false;
                int IsLocal = 0;
                for(int k=0;k<NProcs-1;k++){
                    if(HydroKernelExportFlags[leaf][k]&BitMask){
                        IsLocal ++;
                    }
                }

                if(IsLocal == 0){
                    ActiveHydroParticle[i].LocalUpdateFlags = true;
                    UpdateKernelDensityLocal(i,ActiveHydroParticle,Neighbors,
                            MyID,NProcs,HydroKernelExportFlags,BinarySearch);
                }
#endif // LOCAL_UPDATE
            }
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

        double TimeComm = GetElapsedTime();
        MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
        MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
        TimingResults.HydroKernelCommThisStep += GetElapsedTime()-TimeComm;


        TimeNBS = GetElapsedTime();
        for(int i=0;i<NImportAll;i++){
            HydroKernelImportSend[i].SmoothedNumber = 0.e0;
            HydroKernelImportSend[i].Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                GetSmoothedNumberofNeighbors(HydroKernelExportRecv[i].Pos,
                    HydroKernelExportRecv[i].Kernel,Neighbors,&(HydroKernelImportSend[i].SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
                GetNumberofNeighbors(HydroKernelExportRecv[i].Pos,
                    HydroKernelExportRecv[i].Kernel,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
            HydroKernelImportSend[i].ExtraIterationFlag = HydroKernelExportRecv[i].ExtraIterationFlag;
            if(HydroKernelExportRecv[i].ExtraIterationFlag == true){
                HydroKernelImportSend[i].Rho = CalcRho_i(HydroKernelExportRecv[i].Pos,HydroKernelExportRecv[i].Kernel,
                        HydroKernelImportSend[i].Nlist,Neighbors);
            }
            HydroKernelImportSend[i].Leaf = HydroKernelExportRecv[i].Leaf;
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

        NImportAll = 0;
        int NImportAllNew = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew[i] = 0;
            for(int k=0;k<NImportThisTime[i];k++){
                if(HydroKernelImportSend[NImportAll].Nlist > 0){
                    HydroKernelImportSend[NImportAllNew] = HydroKernelImportSend[NImportAll];
                    NImportThisTimeNew[i] ++;
                    NImportAllNew ++;
                }
                NImportAll ++;
            }
        }

        for(int i=0;i<NProcs-1;i++){
            MPI_Sendrecv(NImportThisTimeNew+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_PRECOMM,
                NExportThisTimeNew+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_PRECOMM,
                    MPI_COMM_WORLD,&mpi_status);
        }

        NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            MPI_Sendrecv(HydroKernelImportSend+NImport,
                NImportThisTimeNew[i]*sizeof(struct StructHydroKernelImport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_IMPORT,
                        HydroKernelImportRecv[i],
                NExportThisTimeNew[i]*sizeof(struct StructHydroKernelImport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_IMPORT,
                        MPI_COMM_WORLD,&mpi_status);
            NImport += NImportThisTimeNew[i];
        }
        TimeComm = GetElapsedTime();
        MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
        MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
        TimingResults.HydroKernelCommThisStep += GetElapsedTime()-TimeComm;

        for(int i=0;i<NProcs-1;i++){
            for(int k=0;k<NExportThisTimeNew[i];k++){ 
                int leaf = HydroKernelImportRecv[i][k].Leaf;
                ActiveHydroParticle[leaf].Nlist += HydroKernelImportRecv[i][k].Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
                ActiveHydroParticle[leaf].SmoothedNumber += HydroKernelImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
                if(HydroKernelImportRecv[i][k].ExtraIterationFlag == true){
                    ActiveHydroParticle[leaf].Rho += HydroKernelImportRecv[i][k].Rho;
                }
#ifdef TASK_TEST_HYDRO_QUANTITIES //{
                if(HydroKernelImportRecv[i][k].Nlist > 0)
                    Phydro[NBCache[leaf].Leaf]->Tag = -1;
#endif // TASK_TEST_HYDRO_QUANTITIES //}
            }
        }
        
        int NLocalActiveLeaves = CheckNeighborNumberAndUpdateKernelDensity(NActives,
                NProcs,HydroKernelExportFlags,ActiveHydroParticle,BinarySearch);
        MPI_Allreduce(&NLocalActiveLeaves,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        Niteration ++;
        if(Niteration > MaxIterationTimes)
            break;
    } while (0<NActiveLeaves);

#else // EVALUATE_KERNEL_BY_ITERATION
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->Kernel = KERNEL_FACTOR*
                pow(Phydro[i]->Mass/Phydro[i]->Rho,1.0/((double)DIMENSION));
            Phydro[i]->Rho = 0.e0;
        }
    }
#endif // EVALUATE_KERNEL_BY_ITERATION

    PlantHydroTreeKernelMaxUpdate();
#ifdef EVALUATE_KERNEL_BY_ITERATION
#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"\n");
#else // PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"%d iterations for kernel determination.\n",Niteration);
#endif // PRINT_LOG_KERNEL_ITERATION
#endif // EVALUATE_KERNEL_BY_ITERATION

    TimingResults.HydroKernelThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

#endif



/*
static double Coefa,Coefb,Coefc,Coefd;
void InitializeKernelEstimationFucntion(void){
    Coefa = 2.4081; 
    Coefb = -2.11875; 
    Coefc = 1.0/(tanh(Coefa*2+Coefb)-tanh(Coefa*0+Coefb));
    Coefd = -Coefc*tanh(Coefa*0+Coefb);
    fprintf(stderr,"Coef for Kernel Estimation = %g %g %g %g\n",Coefa,Coefb,Coefc,Coefd);
    return;
}

static inline void KernelEstimationFunction(const int index, const int Nlist){

    double s;
    if(Pall.Ns>Nlist){
        s = (double)Nlist/(double)Pall.Ns;
        //eprint(s);
    }else {
        s = (double)Pall.Ns/(double)Nlist;
        //eprint(s);
    }

    double x = Coefc*tanh(Coefa*s+Coefb)+Coefd;
    if(index == 0){
        eprint(s);
        eprint(x);
    }
    if(Pall.Ns>Nlist){
        Phydro[index]->Kernel *= 1.0/x;
    }else {
        Phydro[index]->Kernel *= x;
    }

    return;
}
*/
