#include	"config.h"
#include	"PlantHydroTree.h"

static bool FirstCall = true;
void InitializeRootForHydro(void){

    if(FirstCall == false) return;
    FirstCall = false;

    int NloadHydro =(int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)MAX(Pall.Nhydro,NAdditionUnit))/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = MAX(Pall.Nhydro,NAdditionUnit)*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));

    /* HydroRoot */
    HydroRoot.NumberofAllocatedLeaves = NloadHydro;
	HydroRoot.Leaves = malloc(sizeof(int)*NloadHydro+1);
    HydroRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
    HydroRoot.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

	HydroRoot.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForNBS;
    HydroRoot.MaxLevel = TreeMaxNodeLevel;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        HydroRoot.WidthFactor[i] = pow(dw,(double)i);
    /* HydroRoot */

    /* HydroNode */
    HydroNode = malloc(sizeof(struct StructHydroNode)*HydroRoot.NumberofAllocatedNodes);

    struct StructHydroNode HydroNodeTemp;
    memset(&HydroNodeTemp,0,sizeof(struct StructHydroNode));
    HydroNodeTemp.Next = HydroNodeTemp.Parent = HydroNodeTemp.Children = 
    HydroNodeTemp.Sister = NONE;
    //HydroNodeTemp.Sister = HydroNodeTemp.Traverse = NONE;

    for(int i=0;i<HydroRoot.NumberofAllocatedNodes;i++)
        HydroNode[i] = HydroNodeTemp;
    /* HydroNode */

    /* NBCache */
    NBCache = malloc(sizeof(StructNBCache)*NloadHydro+1);
    /* NBCache */


    int NProcs = MPIGetNumProcs();
    if(EdgesForHydro == NULL)
        EdgesForHydro = malloc(sizeof(struct StructEdges)*NProcs);
    if(EdgesForActiveHydro == NULL)
        EdgesForActiveHydro = malloc(sizeof(struct StructEdges)*NProcs);

	return;
}

static void MakeHydroRoot(void);
static void BuildHydroTree(void);
static void BuildHydroTraversalLink(void);
static int NextHydroNode(const int NodeID);
static int TraversalID[TreeMaxNodeLevel];
static void HydroNodeDataImplant(void);
static void HydroNodeDataImplantNew(const int CurrentNodeID);
//static void HydroNodeDataImplantPartial(void);
static void HydroNodeDataUpdateKernlMaxDistanceMaxNActives(void);
static void HydroNodeDataUpdateKernelMax(void);
static void HydroNodeDataUpdate(void);
static void UpdateEdgesForHydro(void);
static void UpdateEdgesForHydroAsyncronoushComm(const int mode);

static int NumberofAllocatedLeaves = 0;
static int *sublist; 
static double DiagHydro;

/*
 * This structure uses only for the operations in the tree structure
 * construction.
 */
static struct StructCachedData{
    double Pos[3];
    double Kernel;
    bool Active;
} *CachedData;

/*
 * Preprocessing subroutine for the hydro tree construction. In this
 * subroutine, data and cache arrays are prepared.
 */
static void HydroTreePreprocessing(void){

    if(MAX(Pall.Nhydro,NAdditionUnit)>HydroRoot.NumberofAllocatedLeaves){
        HydroRoot.NumberofAllocatedLeaves = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        free(HydroRoot.Leaves);
        free(NBCache);
        HydroRoot.Leaves = malloc(sizeof(int)*HydroRoot.NumberofAllocatedLeaves);
        NBCache = malloc(sizeof(StructNBCache)*HydroRoot.NumberofAllocatedLeaves);
    }

    if(MAX(Pall.Nhydro,NAdditionUnit) >NumberofAllocatedLeaves){
        NumberofAllocatedLeaves = (int)(ForAngelsShare*Pall.Nhydro);
        if(NumberofAllocatedLeaves > 0){
            free(sublist);
            free(CachedData);
        }
	    sublist = malloc(sizeof(int)*NumberofAllocatedLeaves);
        CachedData = malloc(sizeof(struct StructCachedData)*NumberofAllocatedLeaves);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        CachedData[i].Pos[0] = PhydroPosP(i)[0];
        CachedData[i].Pos[1] = PhydroPosP(i)[1];
        CachedData[i].Pos[2] = PhydroPosP(i)[2];
        CachedData[i].Kernel = Phydro[i]->KernelPred;
        CachedData[i].Active = Phydro[i]->Active;
    }
    /*
    int counter = 0; int Nhydro = Pall.Nhydro;
    for(int i=0;i<PhydroElementsSize;i++){
        if((PhydroElements[i].Use)&&(counter < Nhydro)){
            CachedData[counter].Pos[0] = PhydroElements[i].PosP[0];
            CachedData[counter].Pos[1] = PhydroElements[i].PosP[1];
            CachedData[counter].Pos[2] = PhydroElements[i].PosP[2];
            CachedData[counter].Kernel = PhydroElements[i].KernelPred;
            CachedData[counter].Active = PhydroElements[i].Active;
            counter ++;
        }
    }
    */

    return ;
}

/*
 * This subroutine reconnects NBCache[]->Leaf to the renewed index of the hydro
 * particle.
 */
static void ConnectPhydroToLeaves(void){
    //
    for(int i=0;i<Pall.Nhydro;i++){
        //Phydro[HydroRoot.Leaves[i]]->Leaf = i;
        Phydro[NBCache[i].Leaf]->Leaf = i;
    }

    return ;
}

/*
 * This subroutine plants a hydro tree structure. 
 */
void PlantHydroTree(void){

    double TimingResultThisRoutine = GetElapsedTime();

    HydroTreePreprocessing();

    MakeHydroRoot();
    BuildHydroTree();
    ConnectPhydroToLeaves();
    UpdateEdgesForHydroAsyncronoushComm(0);

    for(int i=1;i<HydroRoot.NumberofNodes;i++){
        HydroNode[i].Next = NextHydroNode(i);
    }

    //HydroNodeDataImplant();
    DiagHydro = DISTANCE(HydroNode[0].Pos,HydroNode[HydroNode[0].Children].Pos);
    HydroNodeDataImplantNew(0);
    UpdateEdgesForHydroAsyncronoushComm(1);
    //UpdateEdgesForHydro();

    TimingResults.HydroTreeThisStep += GetElapsedTime()-TimingResultThisRoutine;

	return;
}

/*
 * This subroutine only updates infomation of particles(leafs). The subroutine
 * "PlantHydroTree()" should be called before using this subrountine at least
 * at once.
 */
void PlantHydroTreeUpdate(void){

    double TimingResultThisRoutine = GetElapsedTime();
    //assert(Pall.Nhydro == HydroRoot.NumberofLeaves);
#if 0
    static int RootNodeID = 0; 
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    for(int k=0;k<NumberofLeaves;k++){
        if(HydroRoot.Leaves[k] < 0) continue ;
        int leaf = NBCache[k].Leaf;
        NBCache[k].Kernel = Phydro[leaf]->KernelPred;
        NBCache[k].Pos[0] = Phydro[leaf]->PosP[0];
        NBCache[k].Pos[1] = Phydro[leaf]->PosP[1];
        NBCache[k].Pos[2] = Phydro[leaf]->PosP[2];
        NBCache[k].Active = Phydro[leaf]->Active;
    }
#endif

    for(int i=0;i<Pall.Nhydro;i++){
        int leaf = Phydro[i]->Leaf;
        NBCache[leaf].Kernel = Phydro[i]->KernelPred;
        NBCache[leaf].Pos[0] = Phydro[i]->PosP[0];
        NBCache[leaf].Pos[1] = Phydro[i]->PosP[1];
        NBCache[leaf].Pos[2] = Phydro[i]->PosP[2];
        NBCache[leaf].Active = Phydro[i]->Active;
    }
    UpdateEdgesForHydroAsyncronoushComm(0);

    //HydroNodeDataUpdateKernlMaxDistanceMaxNActives();
    //DiagHydro = DISTANCE(HydroNode[0].Pos,HydroNode[HydroNode[0].Children].Pos);
    HydroNodeDataImplantNew(0);
    UpdateEdgesForHydroAsyncronoushComm(1);
    //UpdateEdgesForHydro();
    TimingResults.HydroTreeThisStep += GetElapsedTime()-TimingResultThisRoutine;

	return;
}

void PlantHydroTreeKernelMaxUpdate(void){

    for(int i=0;i<Pall.Nhydro;i++)
        NBCache[Phydro[i]->Leaf].Kernel = Phydro[i]->KernelPred;
    
    //HydroNodeDataUpdateKernelMax();
    DiagHydro = DISTANCE(HydroNode[0].Pos,HydroNode[HydroNode[0].Children].Pos);
    HydroNodeDataImplantNew(0);

	return;
}

static void MakeHydroRoot(void){

	double min[3],max[3];
    int RootNodeID = 0;

    if(Pall.Nhydro>0){
        for(int k=0;k<3;k++){
            max[k] = CachedData[0].Pos[k];
            min[k] = CachedData[0].Pos[k];
        }
        for(int i=1;i<Pall.Nhydro;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(CachedData[i].Pos[k],max[k]);
                min[k] = fmin(CachedData[i].Pos[k],min[k]);
            }
        }
    }else{ // if no hydro particle case.  // use center of mass
        double mass = 0.e0;
        double COM[3] = {0.e0,0.e0,0.e0};
        for(int i=0;i<Pall.Ntotal;i++){
            for(int k=0;k<3;k++)
                COM[k] = Pbody[i]->Mass*Pbody[i]->PosP[k];
            mass += Pbody[i]->Mass;
        }
        double imass = 1.e0/mass;
        max[0] = min[0] = COM[0]*imass;
        max[1] = min[1] = COM[1]*imass;
        max[2] = min[2] = COM[2]*imass;
    }

	for(int k=0;k<3;k++){
        HydroRoot.PosMax[k] = max[k];
        HydroRoot.PosMin[k] = min[k];
    }

    double WidthMax = 0.e0;
	for(int k=0;k<3;k++){
        HydroNode[RootNodeID].Pos[k] = 0.5*(max[k] + min[k]);
        WidthMax = fmax(WidthMax,max[k]-min[k]);
    }
    HydroRoot.Width = WidthMax;

    HydroNode[RootNodeID].Next = NONE;
    HydroNode[RootNodeID].Parent = NONE;
    HydroNode[RootNodeID].Sister = NONE;
    HydroNode[RootNodeID].Children = NONE;
    // HydroNode[RootNodeID].Traverse = NONE;

    HydroNode[RootNodeID].NumberofLeaves = Pall.Nhydro;

    HydroNode[RootNodeID].Level = 0;
    HydroNode[RootNodeID].Leaves = 0;

    for(int i=0;i<Pall.Nhydro;i++)
        HydroRoot.Leaves[i] = i;

    HydroRoot.NumberofLeaves = Pall.Nhydro;

	return ;
}

static int NextHydroNode(const int NodeID){

    int CurrentNodeID = NodeID;

    if(HydroNode[CurrentNodeID].Sister != NONE){
        CurrentNodeID = HydroNode[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(HydroNode[HydroNode[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = HydroNode[HydroNode[NextNodeID].Parent].Sister;
                break;
            } else if(HydroNode[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = HydroNode[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}


static inline bool HydroNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber) __attribute__((always_inline));
static inline bool HydroNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber){

	if( (HydroNode[CurrentNodeID].NumberofLeaves <= CriticalNumber) || HydroNode[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}

static void BuildHydroTreeForException(void){

    int  NumberofNodes = 1;
    int ChildNodeID = 1;
    int CurrentNodeID = 0;
    int CurrentMaxLevel = 0;
    if(NumberofNodes >= HydroRoot.NumberofAllocatedNodes){
        int NumberofAllocatedNodes = (int)(MAX(ForAngelsShare*NumberofNodes,NAdditionUnit));
        HydroNode = realloc(HydroNode,sizeof(struct StructHydroNode)*NumberofAllocatedNodes);
        HydroRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
        dprintlmpi(NumberofAllocatedNodes);
    }
    HydroNode[ChildNodeID].Next = HydroNode[ChildNodeID].Parent = 
    HydroNode[ChildNodeID].Sister = HydroNode[ChildNodeID].Children = NONE;
    //HydroNode[ChildNodeID].Traverse = NONE;

    HydroNode[ChildNodeID].Parent = CurrentNodeID;

    HydroNode[CurrentNodeID].Children = ChildNodeID;
    int NextNodeID = ChildNodeID;
    HydroNode[ChildNodeID].Leaves = HydroNode[CurrentNodeID].Leaves;
    CurrentMaxLevel = MAX(CurrentMaxLevel,HydroNode[CurrentNodeID].Level+1);

    HydroRoot.Leaves[HydroNode[ChildNodeID].Leaves] = HydroRoot.Leaves[HydroNode[CurrentNodeID].Leaves];
    HydroNode[ChildNodeID].NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
    HydroNode[ChildNodeID].Level = HydroNode[CurrentNodeID].Level+1;

    for(int k=0;k<3;k++)
        HydroNode[ChildNodeID].Pos[k] = HydroNode[CurrentNodeID].Pos[k];

    HydroRoot.CurrentMaxLevel = CurrentMaxLevel;
    HydroRoot.NumberofNodes = NumberofNodes + 1;

    return ;
}

static void BuildHydroTree(void){

    int NumberofNodes = 0; 
    int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 

    int NumberofNodeCreationLimit = HydroRoot.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
    if(Pall.Nhydro == 0){
        BuildHydroTreeForException();
        return ;
    }

	while(1){

		if(HydroNodeSeparationCriterion(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(HydroNode[CurrentNodeID].Sister != NONE){
				CurrentNodeID = HydroNode[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(HydroNode[HydroNode[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = HydroNode[HydroNode[NextNodeID].Parent].Sister;
						break;
					} else if(HydroNode[NextNodeID].Parent == RootNodeID){ // End procedure.
                        HydroRoot.CurrentMaxLevel = CurrentMaxLevel;

                        int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
                        for(int k=0;k<NumberofLeaves;k++){
                            int leaf =  HydroRoot.Leaves[k];
                            NBCache[k].Pos[0] = CachedData[leaf].Pos[0];
                            NBCache[k].Pos[1] = CachedData[leaf].Pos[1];
                            NBCache[k].Pos[2] = CachedData[leaf].Pos[2];
                            NBCache[k].Kernel = CachedData[leaf].Kernel;
                            NBCache[k].Active = CachedData[leaf].Active;
                            NBCache[k].Leaf = leaf;
                        }
                        HydroRoot.NumberofNodes = NumberofNodes + 1;
						return;
					}
                    NextNodeID = HydroNode[NextNodeID].Parent;
				}
			}
			continue;
		}

		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

		for(int i=0;i<HydroNode[CurrentNodeID].NumberofLeaves;i++){
            int leaf = HydroRoot.Leaves[HydroNode[CurrentNodeID].Leaves+i];
			int subindex = 0;

            for(int k=0;k<3;k++)
                if(HydroNode[CurrentNodeID].Pos[k] <= CachedData[leaf].Pos[k])
                    subindex += 1 << k;

			if(subnumber[subindex] == 0){
				subhead[subindex] = subcurrent[subindex] = leaf;
			} else {
				sublist[subcurrent[subindex]] = leaf;
				subcurrent[subindex] = leaf;
			}
			subnumber[subindex] ++;
        }


        ChildNodeID = CurrentNodeID;
		for(int i=0;i<TreeNsub;i++){
			if(subnumber[i] != 0){
				BackwardNodeID = ChildNodeID; 
                // make node
                NumberofNodes ++;
                ChildNodeID = NumberofNodes;
                if(NumberofNodes >= HydroRoot.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(ForAngelsShare*NumberofNodes);
                    HydroNode = realloc(HydroNode,sizeof(struct StructHydroNode)*NumberofAllocatedNodes);
                    HydroRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                HydroNode[ChildNodeID].Next = NONE;
                HydroNode[ChildNodeID].Parent = NONE;
                HydroNode[ChildNodeID].Sister = NONE;
                HydroNode[ChildNodeID].Children = NONE;
                //HydroNode[ChildNodeID].Traverse = NONE;
                // make node

                HydroNode[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    HydroNode[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    HydroNode[ChildNodeID].Leaves = HydroNode[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,HydroNode[CurrentNodeID].Level+1);
                } else {
                    HydroNode[BackwardNodeID].Sister = ChildNodeID;
                    HydroNode[ChildNodeID].Leaves = HydroNode[BackwardNodeID].Leaves + HydroNode[BackwardNodeID].NumberofLeaves;
                }

                HydroRoot.Leaves[HydroNode[ChildNodeID].Leaves] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    HydroRoot.Leaves[HydroNode[ChildNodeID].Leaves+k] = 
                        sublist[HydroRoot.Leaves[HydroNode[ChildNodeID].Leaves+k-1]];
                }
                HydroNode[ChildNodeID].NumberofLeaves = subnumber[i];
                HydroNode[ChildNodeID].Level = HydroNode[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    HydroNode[ChildNodeID].Pos[k] = HydroNode[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*HydroRoot.Width*HydroRoot.WidthFactor[HydroNode[CurrentNodeID].Level];

            }
		}
        CurrentNodeID = NextNodeID;
	}
}

#if 0
static void BuildHydroTraversalLink(void){

    int MarkerID[TreeMaxNodeLevel];
    int RootNodeID = 0;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        TraversalID[i] = NONE;
    TraversalID[RootNodeID] = 0;
    HydroNode[RootNodeID].Traverse = NONE;

    for(int i=1;i<HydroRoot.NumberofNodes;i++){
        HydroNode[i].Traverse = NONE;

        if(TraversalID[HydroNode[i].Level] != NONE){
            HydroNode[MarkerID[HydroNode[i].Level]].Traverse = i;
        } else {
            TraversalID[HydroNode[i].Level] = i;
        }
        MarkerID[HydroNode[i].Level] = i; 

        HydroNode[i].Next = NextHydroNode(i);
    }

	return;
}
#endif

static inline void CalcDistanceMaxKernelMaxActiveNumber(const int CurrentNodeID) __attribute__((always_inline));
static inline void CalcDistanceMaxKernelMaxActiveNumber(const int CurrentNodeID){

    double PosLeaf[3] = {HydroNode[CurrentNodeID].Pos[0],HydroNode[CurrentNodeID].Pos[1],HydroNode[CurrentNodeID].Pos[2]};
    double distancemax = 0.e0;
    double kernelmax = 0.e0;

    int NActives = 0;
    int Number_of_leaf = HydroNode[CurrentNodeID].NumberofLeaves;
    int header = HydroNode[CurrentNodeID].Leaves;
    for(int k=0;k<Number_of_leaf;k++){
        int leaf = header+k;
        if(HydroRoot.Leaves[leaf] < 0) continue;
        double distance = DISTANCE(PosLeaf,NBCache[leaf].Pos);
        distancemax = fmax(distance,distancemax);
        kernelmax = fmax(kernelmax,2.e0*NBCache[leaf].Kernel+distance);
        NActives += NBCache[leaf].Active;
    }

    HydroNode[CurrentNodeID].KernelMax = kernelmax;
    HydroNode[CurrentNodeID].DistanceMax = distancemax;
    HydroNode[CurrentNodeID].NumberofActiveLeaves = NActives;

    return;
}

/*
 * This function implants information of hydro nodes, namely DistanceMax,
 * KernelMax, NActives, and Leaves.  This function should be called after the
 * construction of the hydro tree structure.
 */
static void HydroNodeDataImplantNew(const int CurrentNodeID){

    int NActives = 0;
    double KernelMax = 0.e0;
    double DistanceMax = 0.e0;

    if(HydroNode[CurrentNodeID].Children == NONE){
        int Number_of_leaf = HydroNode[CurrentNodeID].NumberofLeaves;
        int header = HydroNode[CurrentNodeID].Leaves;
        for(int k=0;k<Number_of_leaf;k++){
            int leaf = header+k;
            if(HydroRoot.Leaves[leaf] < 0) continue;
            double Distance = DISTANCE(HydroNode[CurrentNodeID].Pos,NBCache[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            KernelMax = fmax(KernelMax,2.e0*NBCache[leaf].Kernel+Distance);
            NActives += NBCache[leaf].Active;
        }
        HydroNode[CurrentNodeID].KernelMax = KernelMax;
        HydroNode[CurrentNodeID].DistanceMax = DistanceMax;
        HydroNode[CurrentNodeID].NumberofActiveLeaves = NActives;
    } else {
        int ChildNodeID = HydroNode[CurrentNodeID].Children;
        do{
            HydroNodeDataImplantNew(ChildNodeID);

            KernelMax = fmax(KernelMax,HydroNode[ChildNodeID].KernelMax);
            DistanceMax = fmax(DistanceMax,HydroNode[ChildNodeID].DistanceMax);
            NActives += HydroNode[ChildNodeID].NumberofActiveLeaves;
            ChildNodeID = HydroNode[ChildNodeID].Sister;
        } while(ChildNodeID != NONE);
        double Width = DiagHydro*HydroRoot.WidthFactor[HydroNode[CurrentNodeID].Level];
        HydroNode[CurrentNodeID].KernelMax = KernelMax + Width;
        HydroNode[CurrentNodeID].DistanceMax = DistanceMax + Width;
        HydroNode[CurrentNodeID].NumberofActiveLeaves = NActives;
    }

    //HydroNode[CurrentNodeID].Leaves = HydroNode[HydroNode[CurrentNodeID].Children].Leaves;

    return;
}

#if 0
/*
 * This function implants information of hydro nodes, namely DistanceMax,
 * KernelMax, NActives, and Leaves.  This function should be called after the
 * construction of the hydro tree structure.
 */
static void HydroNodeDataImplant(void){

    int RootNodeID = 0;
    double diag = DISTANCE(HydroNode[RootNodeID].Pos,HydroNode[HydroNode[RootNodeID].Children].Pos);

    int ChildNodeID;
    for(int i=HydroRoot.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(HydroNode[CurrentNodeID].Children != NONE){

                ChildNodeID = HydroNode[CurrentNodeID].Children;

                double width = diag*HydroRoot.WidthFactor[i];
                double distancemax = 0.e0;
                double kernelmax = 0.e0;
                int NActives = 0;

                while(ChildNodeID != NONE){
                    kernelmax = fmax(kernelmax,HydroNode[ChildNodeID].KernelMax);
                    distancemax = fmax(distancemax,HydroNode[ChildNodeID].DistanceMax);
                    NActives += HydroNode[ChildNodeID].NumberofActiveLeaves;
                    ChildNodeID = HydroNode[ChildNodeID].Sister;
                }

                HydroNode[CurrentNodeID].KernelMax = kernelmax + width;
                HydroNode[CurrentNodeID].DistanceMax = distancemax + width;
                HydroNode[CurrentNodeID].NumberofActiveLeaves = NActives;
                HydroNode[CurrentNodeID].Leaves = HydroNode[HydroNode[CurrentNodeID].Children].Leaves;
            }else{ // No child node case.
                CalcDistanceMaxKernelMaxActiveNumber(CurrentNodeID); 
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}
#endif

#if 0
/*
 * This function updates information of KernelMax, DistanceMax, and NActives in
 * the hydro tree. The tree structure is kept.
 */
static void HydroNodeDataUpdate(void){

    int RootNodeID = 0;
    double diag = DISTANCE(HydroNode[RootNodeID].Pos,HydroNode[HydroNode[RootNodeID].Children].Pos);

    int ChildNodeID;
    for(int i=HydroRoot.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(HydroNode[CurrentNodeID].Children != NONE){

                ChildNodeID = HydroNode[CurrentNodeID].Children;

                double width = diag*HydroRoot.WidthFactor[i];
                double distancemax = 0.e0;
                double kernelmax = 0.e0;
                int NActives = 0;

                while(ChildNodeID != NONE){
                    kernelmax = fmax(kernelmax,HydroNode[ChildNodeID].KernelMax);
                    distancemax = fmax(distancemax,HydroNode[ChildNodeID].DistanceMax);
                    NActives += HydroNode[ChildNodeID].NumberofActiveLeaves;
                    ChildNodeID = HydroNode[ChildNodeID].Sister;
                }

                HydroNode[CurrentNodeID].KernelMax = kernelmax + width;
                HydroNode[CurrentNodeID].DistanceMax = distancemax + width;
                HydroNode[CurrentNodeID].NumberofActiveLeaves = NActives;
                HydroNode[CurrentNodeID].Leaves = HydroNode[HydroNode[CurrentNodeID].Children].Leaves;
            }else{ // No child node case.
                CalcDistanceMaxKernelMaxActiveNumber(CurrentNodeID); 
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}

/*
 * This function updates information of KernelMax, DistanceMax, and NActives in
 * the hydro tree. The tree structure is kept.
 */
static void HydroNodeDataUpdateKernlMaxDistanceMaxNActives(void){

    int RootNodeID = 0;
    double diag = DISTANCE(HydroNode[RootNodeID].Pos,HydroNode[HydroNode[RootNodeID].Children].Pos);

    for(int i=HydroRoot.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(HydroNode[CurrentNodeID].Children != NONE){

                int ChildNodeID = HydroNode[CurrentNodeID].Children;

                double width = diag*HydroRoot.WidthFactor[i];
                double distancemax = 0.e0;
                double kernelmax = 0.e0;
                int NActives = 0;

                while(ChildNodeID != NONE){
                    kernelmax = fmax(kernelmax,HydroNode[ChildNodeID].KernelMax);
                    distancemax = fmax(distancemax,HydroNode[ChildNodeID].DistanceMax);
                    NActives += HydroNode[ChildNodeID].NumberofActiveLeaves;
                    ChildNodeID = HydroNode[ChildNodeID].Sister;
                }

                HydroNode[CurrentNodeID].KernelMax = kernelmax + width;
                HydroNode[CurrentNodeID].DistanceMax = distancemax + width;
                HydroNode[CurrentNodeID].NumberofActiveLeaves = NActives;
            }else{ // No child node case.
                CalcDistanceMaxKernelMaxActiveNumber(CurrentNodeID); 
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}

/*
 * This function updates information of KernelMax in the hydro tree. The tree
 * structure is kept.
 */
static void HydroNodeDataUpdateKernelMax(void){ // only update kernel size

    int RootNodeID = 0;
    double diag = DISTANCE(HydroNode[RootNodeID].Pos,HydroNode[HydroNode[RootNodeID].Children].Pos);

    int ChildNodeID;
    for(int i=HydroRoot.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(HydroNode[CurrentNodeID].Children != NONE){
                double dist_parent_child_nodes = diag*HydroRoot.WidthFactor[i];
                double kernelmax = 0.e0;

                ChildNodeID = HydroNode[CurrentNodeID].Children;
                while(ChildNodeID != NONE){
                    kernelmax = fmax(kernelmax,HydroNode[ChildNodeID].KernelMax);
                    ChildNodeID = HydroNode[ChildNodeID].Sister;
                }
                HydroNode[CurrentNodeID].KernelMax = kernelmax + dist_parent_child_nodes;
            }else{ // No child node case.
                double kernelmax = 0.e0;

                int Number_of_leaf = HydroNode[CurrentNodeID].NumberofLeaves;
                int header = HydroNode[CurrentNodeID].Leaves;
                for(int k=0;k<Number_of_leaf;k++){ // There in no special treatment even when there is no living leaf in this node.
                    int leaf = header+k;
                    if(HydroRoot.Leaves[leaf] < 0) continue;
                    kernelmax = fmax(kernelmax,2.e0*NBCache[leaf].Kernel+DISTANCE(HydroNode[CurrentNodeID].Pos,NBCache[leaf].Pos));

                }
                HydroNode[CurrentNodeID].KernelMax = kernelmax;
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}
#endif

static void UpdateEdgesForHydro(void){

#ifdef USE_BARYON_COMM //{
    int MyID = MPIGetHydroMyID();
    int NProcs = MPIGetHydroNumProcs();
#else // USE_BARYON_COMM //}//{
    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
#endif // USE_BARYON_COMM //}
    struct StructEdges TempEdges[2];

    double max[3],min[3];
    double maxp[3],minp[3];

    if(Pall.Nhydro == 0){
        // Use domain center.
        double PosLeaf[3] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]};
        //fprintf(stderr,"[%d] %g %g %g\n",MyID,PosLeaf[0],PosLeaf[1],PosLeaf[2]);
        for(int i=0;i<2;i++){
            TempEdges[i].PosMax[0] = PosLeaf[0]; TempEdges[i].PosMin[0] = PosLeaf[0];
            TempEdges[i].PosMax[1] = PosLeaf[1]; TempEdges[i].PosMin[1] = PosLeaf[1];
            TempEdges[i].PosMax[2] = PosLeaf[2]; TempEdges[i].PosMin[2] = PosLeaf[2];
        }
    } else {
        for(int k=0;k<3;k++){
            max[k] = NBCache[0].Pos[k];
            min[k] = NBCache[0].Pos[k];
        }
        for(int i=1;i<Pall.Nhydro;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(NBCache[i].Pos[k],max[k]);
                min[k] = fmin(NBCache[i].Pos[k],min[k]);
            }
        }
        TempEdges[0].PosMax[0] = max[0]; TempEdges[0].PosMin[0] = min[0];
        TempEdges[0].PosMax[1] = max[1]; TempEdges[0].PosMin[1] = min[1];
        TempEdges[0].PosMax[2] = max[2]; TempEdges[0].PosMin[2] = min[2];

        // search for the first active hydro particle
        if(Pall.NActivesHydro > 0){
            int first = true;
            for(int i=0;i<Pall.Nhydro;i++){
                if(NBCache[i].Active){
                    if(first == true){
                        for(int k=0;k<3;k++){
                            max[k] = NBCache[i].Pos[k];
                            min[k] = NBCache[i].Pos[k];
                        }
                        first = false; 
                    } else {
                        for(int k=0;k<3;k++){
                            max[k] = fmax(NBCache[i].Pos[k],max[k]);
                            min[k] = fmin(NBCache[i].Pos[k],min[k]);
                        }

                    }
                }
            }
        } else {
            for(int k=0;k<3;k++)
                max[k] = maxp[k] = min[k] = minp[k] = HydroNode[0].Pos[k];
        }

        TempEdges[1].PosMax[0] = max[0]; TempEdges[1].PosMin[0] = min[0];
        TempEdges[1].PosMax[1] = max[1]; TempEdges[1].PosMin[1] = min[1];
        TempEdges[1].PosMax[2] = max[2]; TempEdges[1].PosMin[2] = min[2];
    }

    // MPI_Status  mpi_status;

    EdgesForHydro[MyID] = TempEdges[0];
    EdgesForActiveHydro[MyID] = TempEdges[1];

#if 0
    for(int i=0;i<NProcs-1;i++){
        struct StructEdges RecvEdges[2];
        MPI_Sendrecv(TempEdges,2*sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PLANTTREE_EXTENSITY,
                     RecvEdges,2*sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY,
                     MPI_COMM_WORLD,&mpi_status);
        EdgesForHydro[CommunicationTable[i].RecvRank] = RecvEdges[0];
        EdgesForActiveHydro[CommunicationTable[i].RecvRank] = RecvEdges[1];
    }

    /*
    for(int i=0;i<MPIGetNumProcs()-1;i++){
        int SendRank = CommunicationTable[i].SendRank;
        fprintf(stderr,"H [%d]->[%d] %g %g %g | %g %g %g \n",MPIGetMyID(),CommunicationTable[i].SendRank,
                EdgesForHydro[SendRank].PosMax[0],
                EdgesForHydro[SendRank].PosMax[1],
                EdgesForHydro[SendRank].PosMax[2],
                EdgesForHydro[SendRank].PosMin[0],
                EdgesForHydro[SendRank].PosMin[1],
                EdgesForHydro[SendRank].PosMin[2]);
        fprintf(stderr,"A [%d]->[%d] %g %g %g | %g %g %g \n",MPIGetMyID(),CommunicationTable[i].SendRank,
                EdgesForActiveHydro[SendRank].PosMax[0],
                EdgesForActiveHydro[SendRank].PosMax[1],
                EdgesForActiveHydro[SendRank].PosMax[2],
                EdgesForActiveHydro[SendRank].PosMin[0],
                EdgesForActiveHydro[SendRank].PosMin[1],
                EdgesForActiveHydro[SendRank].PosMin[2]);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    */

#else

#if 0
    double iRecv[12*NProcs];

    iRecv[12*MyID+0] = TempEdges[0].PosMax[0];
    iRecv[12*MyID+1] = TempEdges[0].PosMax[1];
    iRecv[12*MyID+2] = TempEdges[0].PosMax[2];

    iRecv[12*MyID+3] = TempEdges[0].PosMin[0];
    iRecv[12*MyID+4] = TempEdges[0].PosMin[1];
    iRecv[12*MyID+5] = TempEdges[0].PosMin[2];

    iRecv[12*MyID+6] = TempEdges[1].PosMax[0];
    iRecv[12*MyID+7] = TempEdges[1].PosMax[1];
    iRecv[12*MyID+8] = TempEdges[1].PosMax[2];
   
    iRecv[12*MyID+9]  = TempEdges[1].PosMin[0];
    iRecv[12*MyID+10] = TempEdges[1].PosMin[1];
    iRecv[12*MyID+11] = TempEdges[1].PosMin[2];

    MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,iRecv,12,MPI_DOUBLE,MPI_COMM_WORLD);

    for(int i=0;i<NProcs-1;i++){
        int SendRank = CommunicationTable[i].SendRank;
        EdgesForHydro[SendRank].PosMax[0] = iRecv[SendRank*12+0];
        EdgesForHydro[SendRank].PosMax[1] = iRecv[SendRank*12+1];
        EdgesForHydro[SendRank].PosMax[2] = iRecv[SendRank*12+2];
        EdgesForHydro[SendRank].PosMin[0] = iRecv[SendRank*12+3];
        EdgesForHydro[SendRank].PosMin[1] = iRecv[SendRank*12+4];
        EdgesForHydro[SendRank].PosMin[2] = iRecv[SendRank*12+5];

        EdgesForActiveHydro[SendRank].PosMax[0] = iRecv[SendRank*12+6];
        EdgesForActiveHydro[SendRank].PosMax[1] = iRecv[SendRank*12+7];
        EdgesForActiveHydro[SendRank].PosMax[2] = iRecv[SendRank*12+8];
        EdgesForActiveHydro[SendRank].PosMin[0] = iRecv[SendRank*12+9];
        EdgesForActiveHydro[SendRank].PosMin[1] = iRecv[SendRank*12+10];
        EdgesForActiveHydro[SendRank].PosMin[2] = iRecv[SendRank*12+11];
    }
#else
    double iRecv[12*NProcs];

    iRecv[12*MyID+0] = TempEdges[0].PosMax[0];
    iRecv[12*MyID+1] = TempEdges[0].PosMax[1];
    iRecv[12*MyID+2] = TempEdges[0].PosMax[2];

    iRecv[12*MyID+3] = TempEdges[0].PosMin[0];
    iRecv[12*MyID+4] = TempEdges[0].PosMin[1];
    iRecv[12*MyID+5] = TempEdges[0].PosMin[2];

    iRecv[12*MyID+6] = TempEdges[1].PosMax[0];
    iRecv[12*MyID+7] = TempEdges[1].PosMax[1];
    iRecv[12*MyID+8] = TempEdges[1].PosMax[2];
   
    iRecv[12*MyID+9]  = TempEdges[1].PosMin[0];
    iRecv[12*MyID+10] = TempEdges[1].PosMin[1];
    iRecv[12*MyID+11] = TempEdges[1].PosMin[2];

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    // Isend/Irecv
    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
#ifdef USE_BARYON_COMM //{
        MPI_Isend(iRecv+12*MyID,12,MPI_DOUBLE,
                CommunicationTable[i].HydroSendRank,
                TAG_PLANTTREE_EXTENSITY+i,MPI_HYDRO_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);

        MPI_Irecv(iRecv+12*CommunicationTable[i].HydroRecvRank,12,MPI_DOUBLE,
                CommunicationTable[i].HydroRecvRank,
                TAG_PLANTTREE_EXTENSITY+i,MPI_HYDRO_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
#else // USE_BARYON_COMM //}//{
        MPI_Isend(iRecv+12*MyID,12,MPI_DOUBLE,
                CommunicationTable[i].SendRank,
                TAG_PLANTTREE_EXTENSITY+i,MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);

        MPI_Irecv(iRecv+12*CommunicationTable[i].RecvRank,12,MPI_DOUBLE,
                CommunicationTable[i].RecvRank,
                TAG_PLANTTREE_EXTENSITY+i,MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
#endif // USE_BARYON_COMM //}
    }

    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

    for(int i=0;i<NProcs-1;i++){
#ifdef USE_BARYON_COMM //{
        int SendRank = CommunicationTable[i].HydroSendRank;
#else // USE_BARYON_COMM //}//{
        int SendRank = CommunicationTable[i].SendRank;
#endif // USE_BARYON_COMM //}
        EdgesForHydro[SendRank].PosMax[0] = iRecv[SendRank*12+0];
        EdgesForHydro[SendRank].PosMax[1] = iRecv[SendRank*12+1];
        EdgesForHydro[SendRank].PosMax[2] = iRecv[SendRank*12+2];
        EdgesForHydro[SendRank].PosMin[0] = iRecv[SendRank*12+3];
        EdgesForHydro[SendRank].PosMin[1] = iRecv[SendRank*12+4];
        EdgesForHydro[SendRank].PosMin[2] = iRecv[SendRank*12+5];

        EdgesForActiveHydro[SendRank].PosMax[0] = iRecv[SendRank*12+6];
        EdgesForActiveHydro[SendRank].PosMax[1] = iRecv[SendRank*12+7];
        EdgesForActiveHydro[SendRank].PosMax[2] = iRecv[SendRank*12+8];
        EdgesForActiveHydro[SendRank].PosMin[0] = iRecv[SendRank*12+9];
        EdgesForActiveHydro[SendRank].PosMin[1] = iRecv[SendRank*12+10];
        EdgesForActiveHydro[SendRank].PosMin[2] = iRecv[SendRank*12+11];
    }

#endif
#endif



    return ;
}

static bool firstAC = true;
static double *iRecvAC;
static MPI_Status *mpi_status_Export_SendAC;
static MPI_Request *mpi_request_Export_SendAC;
static MPI_Status *mpi_status_Export_RecvAC;
static MPI_Request *mpi_request_Export_RecvAC;

static void UpdateEdgesForHydroAsyncronoushComm(const int mode){

    if(firstAC == true){
        int MyID = MPIGetMyID();
        int NProcs = MPIGetNumProcs();

        mpi_status_Export_SendAC = malloc(sizeof(MPI_Status)*NProcs);
        mpi_status_Export_RecvAC = malloc(sizeof(MPI_Status)*NProcs);
        mpi_request_Export_SendAC = malloc(sizeof(MPI_Request)*NProcs);
        mpi_request_Export_RecvAC = malloc(sizeof(MPI_Request)*NProcs);
        iRecvAC = malloc(sizeof(double)*12*NProcs);
        firstAC = false;
    }


#ifdef USE_BARYON_COMM //{
    if(MPI_HYDRO_COMM_WORLD == MPI_COMM_NULL){
        return ;
    }

    int MyID = MPIGetHydroMyID();
    int NProcs = MPIGetHydroNumProcs();
#else // USE_BARYON_COMM //}//{
    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
#endif // USE_BARYON_COMM //}


    if(mode == 0){

        struct StructEdges TempEdges[2];
        double max[3],min[3];
        double maxp[3],minp[3];

        if(Pall.Nhydro == 0){
            // Use domain center.
            double PosLeaf[3] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]};
            //fprintf(stderr,"[%d] %g %g %g\n",MyID,PosLeaf[0],PosLeaf[1],PosLeaf[2]);
            for(int i=0;i<2;i++){
                TempEdges[i].PosMax[0] = PosLeaf[0]; TempEdges[i].PosMin[0] = PosLeaf[0];
                TempEdges[i].PosMax[1] = PosLeaf[1]; TempEdges[i].PosMin[1] = PosLeaf[1];
                TempEdges[i].PosMax[2] = PosLeaf[2]; TempEdges[i].PosMin[2] = PosLeaf[2];
            }
        } else {
            for(int k=0;k<3;k++){
                max[k] = NBCache[0].Pos[k];
                min[k] = NBCache[0].Pos[k];
            }
            for(int i=1;i<Pall.Nhydro;i++){
                for(int k=0;k<3;k++){
                    max[k] = fmax(NBCache[i].Pos[k],max[k]);
                    min[k] = fmin(NBCache[i].Pos[k],min[k]);
                }
            }
            TempEdges[0].PosMax[0] = max[0]; TempEdges[0].PosMin[0] = min[0];
            TempEdges[0].PosMax[1] = max[1]; TempEdges[0].PosMin[1] = min[1];
            TempEdges[0].PosMax[2] = max[2]; TempEdges[0].PosMin[2] = min[2];

            // search for the first active hydro particle
            if(Pall.NActivesHydro > 0){
                int first = true;
                for(int i=0;i<Pall.Nhydro;i++){
                    if(NBCache[i].Active){
                        if(first == true){
                            for(int k=0;k<3;k++){
                                max[k] = NBCache[i].Pos[k];
                                min[k] = NBCache[i].Pos[k];
                            }
                            first = false; 
                        } else {
                            for(int k=0;k<3;k++){
                                max[k] = fmax(NBCache[i].Pos[k],max[k]);
                                min[k] = fmin(NBCache[i].Pos[k],min[k]);
                            }

                        }
                    }
                }
            } else {
                for(int k=0;k<3;k++)
                    max[k] = maxp[k] = min[k] = minp[k] = HydroNode[0].Pos[k];
            }

            TempEdges[1].PosMax[0] = max[0]; TempEdges[1].PosMin[0] = min[0];
            TempEdges[1].PosMax[1] = max[1]; TempEdges[1].PosMin[1] = min[1];
            TempEdges[1].PosMax[2] = max[2]; TempEdges[1].PosMin[2] = min[2];
        }

        // MPI_Status  mpi_status;
        EdgesForHydro[MyID] = TempEdges[0];
        EdgesForActiveHydro[MyID] = TempEdges[1];


        iRecvAC[12*MyID+0] = TempEdges[0].PosMax[0];
        iRecvAC[12*MyID+1] = TempEdges[0].PosMax[1];
        iRecvAC[12*MyID+2] = TempEdges[0].PosMax[2];

        iRecvAC[12*MyID+3] = TempEdges[0].PosMin[0];
        iRecvAC[12*MyID+4] = TempEdges[0].PosMin[1];
        iRecvAC[12*MyID+5] = TempEdges[0].PosMin[2];

        iRecvAC[12*MyID+6] = TempEdges[1].PosMax[0];
        iRecvAC[12*MyID+7] = TempEdges[1].PosMax[1];
        iRecvAC[12*MyID+8] = TempEdges[1].PosMax[2];
       
        iRecvAC[12*MyID+9]  = TempEdges[1].PosMin[0];
        iRecvAC[12*MyID+10] = TempEdges[1].PosMin[1];
        iRecvAC[12*MyID+11] = TempEdges[1].PosMin[2];

        // Isend/Irecv
        int SendFlag,RecvFlag;
        for(int i=0;i<NProcs-1;i++){
#ifdef USE_BARYON_COMM //{ 
            MPI_Isend(iRecvAC+12*MyID,12,MPI_DOUBLE,
                    CommunicationTable[i].HydroSendRank,
                    TAG_PLANTTREE_EXTENSITY+i,MPI_HYDRO_COMM_WORLD,mpi_request_Export_SendAC+i);
            MPI_Test(mpi_request_Export_SendAC+i,&SendFlag,MPI_STATUS_IGNORE);

            MPI_Irecv(iRecvAC+12*CommunicationTable[i].HydroRecvRank,12,MPI_DOUBLE,
                    CommunicationTable[i].HydroRecvRank,
                    TAG_PLANTTREE_EXTENSITY+i,MPI_HYDRO_COMM_WORLD,mpi_request_Export_RecvAC+i);
            MPI_Test(mpi_request_Export_RecvAC+i,&RecvFlag,MPI_STATUS_IGNORE);
#else // USE_BARYON_COMM //}//{ 
            MPI_Isend(iRecvAC+12*MyID,12,MPI_DOUBLE,
                    CommunicationTable[i].SendRank,
                    TAG_PLANTTREE_EXTENSITY+i,MPI_COMM_WORLD,mpi_request_Export_SendAC+i);
            MPI_Test(mpi_request_Export_SendAC+i,&SendFlag,MPI_STATUS_IGNORE);

            MPI_Irecv(iRecvAC+12*CommunicationTable[i].RecvRank,12,MPI_DOUBLE,
                    CommunicationTable[i].RecvRank,
                    TAG_PLANTTREE_EXTENSITY+i,MPI_COMM_WORLD,mpi_request_Export_RecvAC+i);
            MPI_Test(mpi_request_Export_RecvAC+i,&RecvFlag,MPI_STATUS_IGNORE);
#endif // USE_BARYON_COMM //}
        }

    } else { 

#if 1
        MPI_Waitall(NProcs-1,mpi_request_Export_SendAC,mpi_status_Export_SendAC);
        MPI_Waitall(NProcs-1,mpi_request_Export_RecvAC,mpi_status_Export_RecvAC);

        for(int i=0;i<NProcs-1;i++){
#ifdef USE_BARYON_COMM //{ 
            int SendRank = CommunicationTable[i].HydroSendRank;
#else // USE_BARYON_COMM //}//{ 
            int SendRank = CommunicationTable[i].SendRank;
#endif // USE_BARYON_COMM //}
            EdgesForHydro[SendRank].PosMax[0] = iRecvAC[SendRank*12+0];
            EdgesForHydro[SendRank].PosMax[1] = iRecvAC[SendRank*12+1];
            EdgesForHydro[SendRank].PosMax[2] = iRecvAC[SendRank*12+2];
            EdgesForHydro[SendRank].PosMin[0] = iRecvAC[SendRank*12+3];
            EdgesForHydro[SendRank].PosMin[1] = iRecvAC[SendRank*12+4];
            EdgesForHydro[SendRank].PosMin[2] = iRecvAC[SendRank*12+5];

            EdgesForActiveHydro[SendRank].PosMax[0] = iRecvAC[SendRank*12+6];
            EdgesForActiveHydro[SendRank].PosMax[1] = iRecvAC[SendRank*12+7];
            EdgesForActiveHydro[SendRank].PosMax[2] = iRecvAC[SendRank*12+8];
            EdgesForActiveHydro[SendRank].PosMin[0] = iRecvAC[SendRank*12+9];
            EdgesForActiveHydro[SendRank].PosMin[1] = iRecvAC[SendRank*12+10];
            EdgesForActiveHydro[SendRank].PosMin[2] = iRecvAC[SendRank*12+11];
        }


#if 0
        for(int i=0;i<NProcs;i++){
            fprintf(stderr,"---- [%d] %g %g / %g %g / %g %g\n",MPIGetMyID(),
                    EdgesForHydro[i].PosMin[0],
                    EdgesForHydro[i].PosMax[0],
                    EdgesForHydro[i].PosMin[1],
                    EdgesForHydro[i].PosMax[1],
                    EdgesForHydro[i].PosMin[2],
                    EdgesForHydro[i].PosMax[2]);
        }
#endif



#else
        // MPI_Waitall(NProcs-1,mpi_request_Export_SendAC,mpi_status_Export_SendAC);
        // MPI_Waitall(NProcs-1,mpi_request_Export_RecvAC,mpi_status_Export_RecvAC);

        int Nrest = NProcs -1;
        int CheckArray[NProcs];
        for(int i=0;i<NProcs-1;i++){
            CheckArray[i] = 0;
        }

        int counter = 0;
        do{
            for(int i=0;i<NProcs-1;i++){
                if(CheckArray[i] == 0){
                    int RecvFlag = 0;
                    MPI_Test(mpi_request_Export_RecvAC+i,&RecvFlag,MPI_STATUS_IGNORE);
                    // MPI_Iprobe(CommunicationTable[i].RecvRank,
                            // TAG_PLANTTREE_EXTENSITY+i,MPI_COMM_WORLD,&RecvFlag,MPI_STATUS_IGNORE);
                    if(RecvFlag != 0){
                        int SendRank = CommunicationTable[i].SendRank;
                        EdgesForHydro[SendRank].PosMax[0] = iRecvAC[SendRank*12+0];
                        EdgesForHydro[SendRank].PosMax[1] = iRecvAC[SendRank*12+1];
                        EdgesForHydro[SendRank].PosMax[2] = iRecvAC[SendRank*12+2];
                        EdgesForHydro[SendRank].PosMin[0] = iRecvAC[SendRank*12+3];
                        EdgesForHydro[SendRank].PosMin[1] = iRecvAC[SendRank*12+4];
                        EdgesForHydro[SendRank].PosMin[2] = iRecvAC[SendRank*12+5];

                        EdgesForActiveHydro[SendRank].PosMax[0] = iRecvAC[SendRank*12+6];
                        EdgesForActiveHydro[SendRank].PosMax[1] = iRecvAC[SendRank*12+7];
                        EdgesForActiveHydro[SendRank].PosMax[2] = iRecvAC[SendRank*12+8];
                        EdgesForActiveHydro[SendRank].PosMin[0] = iRecvAC[SendRank*12+9];
                        EdgesForActiveHydro[SendRank].PosMin[1] = iRecvAC[SendRank*12+10];
                        EdgesForActiveHydro[SendRank].PosMin[2] = iRecvAC[SendRank*12+11];
                        CheckArray[i] = 1;
                        Nrest --;
                    }
                }
            }
            counter ++;
        } while(Nrest != 0);
        MPI_Allreduce(MPI_IN_PLACE,&counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        if(MPIGetMyID() == MPI_ROOT_RANK)
            dprintlmpi(counter);

#endif
    }

    return ;
}

