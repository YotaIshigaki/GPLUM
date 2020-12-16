#include	"config.h"
#include	"PlantStellarTree.h"

static bool StellarTreeInitialized = false;
void InitializeRootForStellar(void){

    if(StellarTreeInitialized == true) return;

    int NloadStellar =(int)(ForAngelsShare*Pall.Nstars);
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)Pall.Nstars)/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = Pall.Nstars*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));

    /* StellarRoot */
    StellarRoot.NumberofAllocatedLeaves = NloadStellar;
	StellarRoot.Leaves = malloc(sizeof(int)*NloadStellar+1);
    StellarRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
    StellarRoot.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

	StellarRoot.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForNBS;
    StellarRoot.MaxLevel = TreeMaxNodeLevel;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        StellarRoot.WidthFactor[i] = pow(dw,(double)i);
    /* StellarRoot */

    /* StellarNode */
    StellarNode = malloc(sizeof(struct StructHydroNode)*StellarRoot.NumberofAllocatedNodes);

    struct StructHydroNode TempStellarNode;
    memset(&TempStellarNode,0,sizeof(struct StructHydroNode));
    TempStellarNode.Next = TempStellarNode.Parent = TempStellarNode.Children = 
    TempStellarNode.Sister = NONE;
    //TempStellarNode.Sister = TempStellarNode.Traverse = NONE;

    for(int i=0;i<StellarRoot.NumberofAllocatedNodes;i++)
        StellarNode[i] = TempStellarNode;
    /* StellarNode */

    /* NBCache */
    StellarNBCache = malloc(sizeof(StructNBCache)*NloadStellar+1);
    /* NBCache */


    int NProcs = MPIGetNumProcs();
    if(EdgesForStars == NULL)
        EdgesForStars = malloc(sizeof(struct StructEdges)*NProcs);

    StellarTreeInitialized = true;

	return;
}

static void MakeStellarRoot(void);
static void BuildStellarTree(void);
static void BuildStellarTraversalLink(void);
static int NextStellarNode(const int NodeID);
static int TraversalID[TreeMaxNodeLevel];
static void StellarNodeDataImplant(void);
static void StellarNodeDataImplantNew(const int CurrentNodeID);
static void StellarNodeDataImplantPartial(void);
static void UpdateEdgesForStars(void);
static void UpdateEdgesForStarsPartial(void);

static int NumberofAllocatedLeaves = 0;
static int *sublist; 
static double DiagStellar;

static struct StructCachedData{
    double Pos[3];
    double Kernel;
    bool Active;
} *CachedData;

static void StellarTreePreprocessing(void){

    if(Pall.Nstars>StellarRoot.NumberofAllocatedLeaves){
        StellarRoot.NumberofAllocatedLeaves = (int)(ForAngelsShare*Pall.Nstars);
        free(StellarRoot.Leaves);
        free(StellarNBCache);
        StellarRoot.Leaves = malloc(sizeof(int)*StellarRoot.NumberofAllocatedLeaves);
        StellarNBCache = malloc(sizeof(StructNBCache)*StellarRoot.NumberofAllocatedLeaves);
    }

    if(NumberofAllocatedLeaves < Pall.Nstars){
        NumberofAllocatedLeaves = (int)(ForAngelsShare*Pall.Nstars);
        if(NumberofAllocatedLeaves > 0){
            free(sublist);
            free(CachedData);
        }
	    sublist = malloc(sizeof(int)*NumberofAllocatedLeaves);
        CachedData = malloc(sizeof(struct StructCachedData)*NumberofAllocatedLeaves);
    }

    for(int i=0;i<Pall.Nstars;i++){
        CachedData[i].Pos[0] = PstarPosP(i)[0];
        CachedData[i].Pos[1] = PstarPosP(i)[1];
        CachedData[i].Pos[2] = PstarPosP(i)[2];
        CachedData[i].Kernel = PstarBody(i)->Eps;
        CachedData[i].Active = PstarActive(i);
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

void PlantStellarTree(void){

    if(StellarTreeInitialized == false)
        InitializeRootForStellar();

    StellarTreePreprocessing();

	MakeStellarRoot();
    BuildStellarTree();
    //BuildStellarTraversalLink();
    //StellarNodeDataImplant();
    for(int i=1;i<StellarRoot.NumberofNodes;i++){
        StellarNode[i].Next = NextStellarNode(i);
    }

    DiagStellar = DISTANCE(StellarNode[0].Pos,StellarNode[StellarNode[0].Children].Pos);
    StellarNodeDataImplantNew(0);
    UpdateEdgesForStars();

	return;
}

void PlantStellarTreePartial(void){

    assert(Pall.Nstars == StellarRoot.NumberofLeaves);

    static int RootNodeID = 0; 
    int Number_of_leaf = StellarNode[RootNodeID].NumberofLeaves;
    for(int k=0;k<Number_of_leaf;k++){
        int leaf =  StellarRoot.Leaves[k];
        StellarNBCache[k].Kernel = PstarBody(leaf)->Eps;
    }

    //StellarNodeDataImplantPartial();
    DiagStellar = DISTANCE(StellarNode[0].Pos,StellarNode[StellarNode[0].Children].Pos);
    StellarNodeDataImplantNew(0);
    UpdateEdgesForStarsPartial();

	return;
}

static void MakeStellarRoot(void){

	double min[3],max[3];
    int RootNodeID = 0;

    if(Pall.Nstars>0){
        for(int k=0;k<3;k++){
            max[k] = CachedData[0].Pos[k];
            min[k] = CachedData[0].Pos[k];
        }
        for(int i=1;i<Pall.Nstars;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(CachedData[i].Pos[k],max[k]);
                min[k] = fmin(CachedData[i].Pos[k],min[k]);
            }
        }
    }else{ // if no hsydro particle case.  // use center of mass
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
        StellarRoot.PosMax[k] = max[k];
        StellarRoot.PosMin[k] = min[k];
    }

    double WidthMax = 0.e0;
	for(int k=0;k<3;k++){
        StellarNode[RootNodeID].Pos[k] = 0.5*(max[k] + min[k]);
        WidthMax = fmax(WidthMax,max[k]-min[k]);
    }
    StellarRoot.Width = WidthMax;

    StellarNode[RootNodeID].Next = NONE;
    StellarNode[RootNodeID].Parent = NONE;
    StellarNode[RootNodeID].Sister = NONE;
    StellarNode[RootNodeID].Children = NONE;
    // StellarNode[RootNodeID].Traverse = NONE;

    StellarNode[RootNodeID].NumberofLeaves = Pall.Nstars;

    StellarNode[RootNodeID].Level = 0;
    StellarNode[RootNodeID].Leaves = 0;

    for(int i=0;i<Pall.Nstars;i++)
        StellarRoot.Leaves[i] = i;

    StellarRoot.NumberofLeaves = Pall.Nstars;

	return ;
}

static int NextStellarNode(const int NodeID){

    int CurrentNodeID = NodeID;

    if(StellarNode[CurrentNodeID].Sister != NONE){
        CurrentNodeID = StellarNode[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(StellarNode[StellarNode[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = StellarNode[StellarNode[NextNodeID].Parent].Sister;
                break;
            } else if(StellarNode[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = StellarNode[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}


static inline bool StellarNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber) __attribute__((always_inline));
static inline bool StellarNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber){

	if( (StellarNode[CurrentNodeID].NumberofLeaves <= CriticalNumber) || StellarNode[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}

static void BuildStellarTree(void){

    int NumberofNodes = 0; 
    int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 

    int NumberofNodeCreationLimit = StellarRoot.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
	while(1){

		if(StellarNodeSeparationCriterion(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(StellarNode[CurrentNodeID].Sister != NONE){
				CurrentNodeID = StellarNode[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(StellarNode[StellarNode[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = StellarNode[StellarNode[NextNodeID].Parent].Sister;
						break;
					} else if(StellarNode[NextNodeID].Parent == RootNodeID){
                        StellarRoot.CurrentMaxLevel = CurrentMaxLevel;

                        int NumberofLeaves = StellarNode[RootNodeID].NumberofLeaves;
                        for(int k=0;k<NumberofLeaves;k++){
                            int leaf =  StellarRoot.Leaves[k];
                            StellarNBCache[k].Pos[0] = CachedData[leaf].Pos[0];
                            StellarNBCache[k].Pos[1] = CachedData[leaf].Pos[1];
                            StellarNBCache[k].Pos[2] = CachedData[leaf].Pos[2];
                            StellarNBCache[k].Kernel = CachedData[leaf].Kernel;
                            StellarNBCache[k].Active = CachedData[leaf].Active;
                            StellarNBCache[k].Leaf = leaf;
                        }
                        StellarRoot.NumberofNodes = NumberofNodes + 1;
						return;
					}
                    NextNodeID = StellarNode[NextNodeID].Parent;
				}
			}
			continue;
		}

		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

		for(int i=0;i<StellarNode[CurrentNodeID].NumberofLeaves;i++){
            int leaf = StellarRoot.Leaves[StellarNode[CurrentNodeID].Leaves+i];
			int subindex = 0;

            for(int k=0;k<3;k++)
                if(StellarNode[CurrentNodeID].Pos[k] <= CachedData[leaf].Pos[k])
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
                if(NumberofNodes >= StellarRoot.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(ForAngelsShare*NumberofNodes);
                    StellarNode = realloc(StellarNode,sizeof(struct StructHydroNode)*NumberofAllocatedNodes);
                    StellarRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                StellarNode[ChildNodeID].Next = NONE;
                StellarNode[ChildNodeID].Parent = NONE;
                StellarNode[ChildNodeID].Sister = NONE;
                StellarNode[ChildNodeID].Children = NONE;
                // StellarNode[ChildNodeID].Traverse = NONE;
                // make node

                StellarNode[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    StellarNode[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    StellarNode[ChildNodeID].Leaves = StellarNode[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,StellarNode[CurrentNodeID].Level+1);
                } else {
                    StellarNode[BackwardNodeID].Sister = ChildNodeID;
                    StellarNode[ChildNodeID].Leaves = StellarNode[BackwardNodeID].Leaves + StellarNode[BackwardNodeID].NumberofLeaves;
                }

                StellarRoot.Leaves[StellarNode[ChildNodeID].Leaves] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    StellarRoot.Leaves[StellarNode[ChildNodeID].Leaves+k] = 
                        sublist[StellarRoot.Leaves[StellarNode[ChildNodeID].Leaves+k-1]];
                }
                StellarNode[ChildNodeID].NumberofLeaves = subnumber[i];
                StellarNode[ChildNodeID].Level = StellarNode[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    StellarNode[ChildNodeID].Pos[k] = StellarNode[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*StellarRoot.Width*StellarRoot.WidthFactor[StellarNode[CurrentNodeID].Level];

                //if(StellarNodeSeparationCriterion(ChildNodeID,NumberofNodeCreationLimit)){
                    //for(int k=0;k<StellarNode[ChildNodeID].NumberofLeaves;k++)
                        //Pstars[StellarRoot.Leaves[StellarNode[ChildNodeID].Leaves+k]]->HostNode = ChildNodeID; 
                //}
            }
		}
        CurrentNodeID = NextNodeID;
	}
}

#if 0
static void BuildStellarTraversalLink(void){

    int MarkerID[TreeMaxNodeLevel];
    int RootNodeID = 0;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        TraversalID[i] = NONE;
    TraversalID[RootNodeID] = 0;
    StellarNode[RootNodeID].Traverse = NONE;

    for(int i=1;i<StellarRoot.NumberofNodes;i++){
        StellarNode[i].Traverse = NONE;

        if(TraversalID[StellarNode[i].Level] != NONE){
            StellarNode[MarkerID[StellarNode[i].Level]].Traverse = i;
        } else {
            TraversalID[StellarNode[i].Level] = i;
        }
        MarkerID[StellarNode[i].Level] = i; 

        StellarNode[i].Next = NextStellarNode(i);
    }

	return;
}
#endif

#if 0
static inline void CalcDistanceMaxKernelMaxActiveNumber(const int NodeID) __attribute__((always_inline));
static inline void CalcDistanceMaxKernelMaxActiveNumber(const int CurrentNodeID){

    double PosLeaf[3] = {StellarNode[CurrentNodeID].Pos[0],StellarNode[CurrentNodeID].Pos[1],StellarNode[CurrentNodeID].Pos[2]};
    double distance,kernel;
    double distancemax = 0.e0;
    double kernelmax = 0.e0;

    int NActives = 0;
    int Number_of_leaf = StellarNode[CurrentNodeID].NumberofLeaves;
    int header = StellarNode[CurrentNodeID].Leaves;
    for(int k=0;k<Number_of_leaf;k++){
        int leaf = header+k;
        distance = DISTANCE(PosLeaf,StellarNBCache[leaf].Pos);
        distancemax = fmax(distance,distancemax);
        kernel = 2.e0*StellarNBCache[leaf].Kernel;
        kernelmax = fmax(kernelmax,kernel+distance);
        NActives += StellarNBCache[leaf].Active;
    }

    StellarNode[CurrentNodeID].KernelMax = kernelmax;
    StellarNode[CurrentNodeID].DistanceMax = distancemax;
    StellarNode[CurrentNodeID].NumberofActiveLeaves = NActives;

    return;
}

static void StellarNodeDataImplant(void){

    int RootNodeID = 0;
    double diag = DISTANCE(StellarNode[RootNodeID].Pos,StellarNode[StellarNode[RootNodeID].Children].Pos);

    int ChildNodeID;
    for(int i=StellarRoot.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(StellarNode[CurrentNodeID].Children != NONE){

                ChildNodeID = StellarNode[CurrentNodeID].Children;

                double width = diag*StellarRoot.WidthFactor[i];
                double distancemax = 0.e0;
                double kernelmax = 0.e0;
                int NActives = 0;

                while(ChildNodeID != NONE){
                    kernelmax = fmax(kernelmax,StellarNode[ChildNodeID].KernelMax);
                    distancemax = fmax(distancemax,StellarNode[ChildNodeID].DistanceMax);
                    NActives += StellarNode[ChildNodeID].NumberofActiveLeaves;
                    ChildNodeID = StellarNode[ChildNodeID].Sister;
                }

                StellarNode[CurrentNodeID].KernelMax = kernelmax + width;
                StellarNode[CurrentNodeID].DistanceMax = distancemax + width;
                StellarNode[CurrentNodeID].NumberofActiveLeaves = NActives;
                StellarNode[CurrentNodeID].Leaves = StellarNode[StellarNode[CurrentNodeID].Children].Leaves;
            }else{ // No child node case.
                CalcDistanceMaxKernelMaxActiveNumber(CurrentNodeID); 
            }
            CurrentNodeID = StellarNode[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}
#endif

static void StellarNodeDataImplantNew(const int CurrentNodeID){

    double Width;
    int NActives = 0;
    double DistanceMax = 0.e0;
    double KernelMax = 0.e0;
    if(StellarNode[CurrentNodeID].Children == NONE){
        int Number_of_leaf = StellarNode[CurrentNodeID].NumberofLeaves;
        int header = StellarNode[CurrentNodeID].Leaves;
        for(int k=0;k<Number_of_leaf;k++){
            int leaf = header+k;
            double Distance = DISTANCE(StellarNode[CurrentNodeID].Pos,StellarNBCache[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            KernelMax = fmax(KernelMax,2.e0*StellarNBCache[leaf].Kernel+Distance);
            NActives += StellarNBCache[leaf].Active;
        }
        Width = 0.e0;
    } else {
        int ChildNodeID = StellarNode[CurrentNodeID].Children;
        do{
            StellarNodeDataImplantNew(ChildNodeID);

            KernelMax = fmax(KernelMax,StellarNode[ChildNodeID].KernelMax);
            DistanceMax = fmax(DistanceMax,StellarNode[ChildNodeID].DistanceMax);
            NActives += StellarNode[ChildNodeID].NumberofActiveLeaves;
            ChildNodeID = StellarNode[ChildNodeID].Sister;
        }while(ChildNodeID != NONE);
        Width = DiagStellar*StellarRoot.WidthFactor[StellarNode[CurrentNodeID].Level];
    }

    StellarNode[CurrentNodeID].KernelMax = KernelMax + Width;
    StellarNode[CurrentNodeID].DistanceMax = DistanceMax + Width;
    StellarNode[CurrentNodeID].NumberofActiveLeaves = NActives;

    return;
}

#if 0
static void StellarNodeDataImplantPartial(void){ // only update kernel size

    int RootNodeID = 0;
    double diag = DISTANCE(StellarNode[RootNodeID].Pos,StellarNode[StellarNode[RootNodeID].Children].Pos);

    int ChildNodeID;
    for(int i=StellarRoot.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(StellarNode[CurrentNodeID].Children != NONE){
                double width = diag*StellarRoot.WidthFactor[i];
                double kernelmax = 0.e0;

                ChildNodeID = StellarNode[CurrentNodeID].Children;
                while(ChildNodeID != NONE){
                    kernelmax = fmax(kernelmax,StellarNode[ChildNodeID].KernelMax);
                    ChildNodeID = StellarNode[ChildNodeID].Sister;
                }
                StellarNode[CurrentNodeID].KernelMax = kernelmax + width;
            }else{ // No child node case.
                double kernelmax = 0.e0;

                int Number_of_leaf = StellarNode[CurrentNodeID].NumberofLeaves;
                int header = StellarNode[CurrentNodeID].Leaves;
                for(int k=0;k<Number_of_leaf;k++){
                    int leaf = header+k;
                    kernelmax = fmax(kernelmax,
                            2.e0*StellarNBCache[leaf].Kernel+DISTANCE(StellarNode[CurrentNodeID].Pos,StellarNBCache[leaf].Pos));
                }
                StellarNode[CurrentNodeID].KernelMax = kernelmax;
            }
            CurrentNodeID = StellarNode[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}
#endif

static void UpdateEdgesForStars(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    struct StructEdges TempEdges;

    double max[3],min[3];
    double maxp[3],minp[3];

    if(Pall.Nstars == 0){
        double PosLeaf[3] = {StellarNode[0].Pos[0],StellarNode[0].Pos[1],StellarNode[0].Pos[2]};
        TempEdges.PosMax[0] = PosLeaf[0]; TempEdges.PosMin[0] = PosLeaf[0];
        TempEdges.PosMax[1] = PosLeaf[1]; TempEdges.PosMin[1] = PosLeaf[1];
        TempEdges.PosMax[2] = PosLeaf[2]; TempEdges.PosMin[2] = PosLeaf[2];
    } else {
        for(int k=0;k<3;k++){
            max[k] = StellarNBCache[0].Pos[k];
            min[k] = StellarNBCache[0].Pos[k];
        }
        for(int i=1;i<Pall.Nstars;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(StellarNBCache[i].Pos[k],max[k]);
                min[k] = fmin(StellarNBCache[i].Pos[k],min[k]);
            }
        }
        TempEdges.PosMax[0] = max[0]; TempEdges.PosMin[0] = min[0];
        TempEdges.PosMax[1] = max[1]; TempEdges.PosMin[1] = min[1];
        TempEdges.PosMax[2] = max[2]; TempEdges.PosMin[2] = min[2];
    }

    MPI_Status  mpi_status;
    
    EdgesForStars[MyID] = TempEdges;

    for(int i=0;i<NProcs-1;i++){
        struct StructEdges RecvEdges;
        MPI_Sendrecv(&TempEdges,sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PLANTTREE_EXTENSITY,
            &RecvEdges,sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY,
            MPI_COMM_WORLD,&mpi_status);
        EdgesForStars[CommunicationTable[i].RecvRank] = RecvEdges;
    }

    return;
}

static void UpdateEdgesForStarsPartial(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    struct StructEdges TempEdges;

    double maxp[3],minp[3];

    if(Pall.Nstars == 0){
        double PosLeaf[3] = {StellarNode[0].Pos[0],StellarNode[0].Pos[1],StellarNode[0].Pos[2]};
        TempEdges.PosMax[0] = PosLeaf[0]; TempEdges.PosMin[0] = PosLeaf[0];
        TempEdges.PosMax[1] = PosLeaf[1]; TempEdges.PosMin[1] = PosLeaf[1];
        TempEdges.PosMax[2] = PosLeaf[2]; TempEdges.PosMin[2] = PosLeaf[2];
    } else {
        for(int k=0;k<3;k++){
            maxp[k] = StellarNBCache[0].Pos[k]+2.0*StellarNBCache[0].Kernel;
            minp[k] = StellarNBCache[0].Pos[k]-2.0*StellarNBCache[0].Kernel;
        }
        for(int i=1;i<Pall.Nstars;i++){
            for(int k=0;k<3;k++){
                maxp[k] = fmax(StellarNBCache[i].Pos[k]+2.0*StellarNBCache[i].Kernel,maxp[k]);
                minp[k] = fmin(StellarNBCache[i].Pos[k]-2.0*StellarNBCache[i].Kernel,minp[k]);
            }
        }
        TempEdges.PosMax[0] = maxp[0]; TempEdges.PosMin[0] = minp[0];
        TempEdges.PosMax[1] = maxp[1]; TempEdges.PosMin[1] = minp[1];
        TempEdges.PosMax[2] = maxp[2]; TempEdges.PosMin[2] = minp[2];
    }

    MPI_Status  mpi_status;

    EdgesForStars[MyID] = TempEdges;

    for(int i=0;i<NProcs-1;i++){
        struct StructEdges RecvEdges;
        MPI_Sendrecv(&TempEdges,sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PLANTTREE_EXTENSITY,
            &RecvEdges,sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY,
            MPI_COMM_WORLD,&mpi_status);
        EdgesForStars[CommunicationTable[i].RecvRank] = RecvEdges;
    }

    return;
}
