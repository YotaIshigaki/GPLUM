#include	"config.h"
#include	"HydroMisc.h"
#include	"PlantHydroTreeImported.h"

void InitializeRootForHydroImported(void){

    int NloadHydro = FirstAllocationSize;
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)FirstAllocationSize)/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = FirstAllocationSize*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));

    /* HydroRoot */
    HydroRootImported.NumberofAllocatedLeaves = NloadHydro;
	HydroRootImported.Leaves = malloc(sizeof(int)*NloadHydro+1);
    HydroRootImported.NumberofAllocatedNodes = NumberofAllocatedNodes;
    HydroRootImported.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

	HydroRootImported.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForNBS;
    HydroRootImported.MaxLevel = TreeMaxNodeLevel;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        HydroRootImported.WidthFactor[i] = pow(dw,(double)i);
    /* HydroRoot */

    /* HydroNode */
    HydroNodeImported = malloc(sizeof(struct StructHydroNode)*HydroRootImported.NumberofAllocatedNodes);

    struct StructHydroNode HydroNodeTemp;
    memset(&HydroNodeTemp,0,sizeof(struct StructHydroNode));
    HydroNodeTemp.Next = HydroNodeTemp.Parent = HydroNodeTemp.Children = 
    HydroNodeTemp.Sister = NONE;
    //HydroNodeTemp.Sister = HydroNodeTemp.Traverse = NONE;

    for(int i=0;i<HydroRootImported.NumberofAllocatedNodes;i++)
        HydroNodeImported[i] = HydroNodeTemp;
    /* HydroNode */

    /* NBCache */
    NBImportedCache = malloc(sizeof(StructNBCache)*NloadHydro+1);
    /* NBCache */

	return;
}

static void MakeHydroRootImported(const int NumberofLeaves);
static void BuildHydroTreeImported(void);
static void BuildHydroTraversalLinkImported(void);
static int NextHydroNodeImported(const int NodeID);
static int TraversalID[TreeMaxNodeLevel];
static void HydroNodeDataImplantImported(void);
static void HydroNodeDataImplantImportedNew(const int CurrentNodeID);

static int NumberofAllocatedLeaves = 0;
static int *sublist; 
static double DiagHydroImported;

static struct StructCachedData{
    double Pos[3];
    double Kernel;
    bool Active;
} *CachedData;

static void HydroTreePreprocessingImported(const int NumberofLeaves, struct StructHydroAccExport HydroAccExportRecv[restrict]){

    if(NumberofLeaves>HydroRootImported.NumberofAllocatedLeaves){
        HydroRootImported.NumberofAllocatedLeaves = (int)(ForAngelsShare*NumberofLeaves);
        free(HydroRootImported.Leaves);
        free(NBImportedCache);
        HydroRootImported.Leaves = malloc(sizeof(int)*HydroRootImported.NumberofAllocatedLeaves);
        NBImportedCache = malloc(sizeof(StructNBCache)*HydroRootImported.NumberofAllocatedLeaves);
    }

    if(NumberofAllocatedLeaves < NumberofLeaves){
        if(NumberofAllocatedLeaves > 0){
            free(sublist);
            free(CachedData);
        }
        NumberofAllocatedLeaves = (int)(ForAngelsShare*NumberofLeaves);
	    sublist = malloc(sizeof(int)*NumberofAllocatedLeaves);
        CachedData = malloc(sizeof(struct StructCachedData)*NumberofAllocatedLeaves);
    }

    for(int i=0;i<NumberofLeaves;i++){
        CachedData[i].Pos[0] = HydroAccExportRecv[i].Pos[0];
        CachedData[i].Pos[1] = HydroAccExportRecv[i].Pos[1];
        CachedData[i].Pos[2] = HydroAccExportRecv[i].Pos[2];
        CachedData[i].Kernel = HydroAccExportRecv[i].Kernel;
        CachedData[i].Active = true;
    }

    return ;
}

void PlantHydroTreeImported(const int NumberofLeaves, struct StructHydroAccExport HydroAccExportRecv[restrict]){

    static bool first = true;
    if(first){
        InitializeRootForHydroImported();
        first = false;
    }

    HydroTreePreprocessingImported(NumberofLeaves,HydroAccExportRecv);

	MakeHydroRootImported(NumberofLeaves);
    BuildHydroTreeImported();
    //BuildHydroTraversalLinkImported();

    //HydroNodeDataImplantImported();
    for(int i=1;i<HydroRootImported.NumberofNodes;i++){
        HydroNodeImported[i].Next = NextHydroNodeImported(i);
    }
    DiagHydroImported = DISTANCE(HydroNodeImported[0].Pos,HydroNodeImported[HydroNodeImported[0].Children].Pos);
    HydroNodeDataImplantImportedNew(0);

	return;
}

static void MakeHydroRootImported(const int NumberofLeaves){

	double min[3],max[3];
    int RootNodeID = 0;

    if(NumberofLeaves>0){
        for(int k=0;k<3;k++){
            max[k] = CachedData[0].Pos[k];
            min[k] = CachedData[0].Pos[k];
        }
        for(int i=1;i<NumberofLeaves;i++){
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
        HydroRootImported.PosMax[k] = max[k];
        HydroRootImported.PosMin[k] = min[k];
    }

    double WidthMax = 0.e0;
	for(int k=0;k<3;k++){
        HydroNodeImported[RootNodeID].Pos[k] = 0.5*(max[k] + min[k]);
        WidthMax = fmax(WidthMax,max[k]-min[k]);
    }
    HydroRootImported.Width = WidthMax;

    HydroNodeImported[RootNodeID].Next = NONE;
    HydroNodeImported[RootNodeID].Parent = NONE;
    HydroNodeImported[RootNodeID].Sister = NONE;
    HydroNodeImported[RootNodeID].Children = NONE;
    // HydroNodeImported[RootNodeID].Traverse = NONE;

    HydroNodeImported[RootNodeID].NumberofLeaves = NumberofLeaves;

    HydroNodeImported[RootNodeID].Level = 0;
    HydroNodeImported[RootNodeID].Leaves = 0;

    for(int i=0;i<NumberofLeaves;i++)
        HydroRootImported.Leaves[i] = i;

    HydroRootImported.NumberofLeaves = NumberofLeaves;

	return ;
}

static int NextHydroNodeImported(const int NodeID){

    int CurrentNodeID = NodeID;

    if(HydroNodeImported[CurrentNodeID].Sister != NONE){
        CurrentNodeID = HydroNodeImported[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(HydroNodeImported[HydroNodeImported[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = HydroNodeImported[HydroNodeImported[NextNodeID].Parent].Sister;
                break;
            } else if(HydroNodeImported[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = HydroNodeImported[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}


static inline bool HydroNodeSeparationCriterionImported(const int CurrentNodeID, const int CriticalNumber) __attribute__((always_inline));
static inline bool HydroNodeSeparationCriterionImported(const int CurrentNodeID, const int CriticalNumber){

	if( (HydroNodeImported[CurrentNodeID].NumberofLeaves <= CriticalNumber) || HydroNodeImported[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}

static void BuildHydroTreeImported(void){

    int NumberofNodes = 0; 
    int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 

    int NumberofNodeCreationLimit = HydroRootImported.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
	while(1){

		if(HydroNodeSeparationCriterionImported(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(HydroNodeImported[CurrentNodeID].Sister != NONE){
				CurrentNodeID = HydroNodeImported[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(HydroNodeImported[HydroNodeImported[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = HydroNodeImported[HydroNodeImported[NextNodeID].Parent].Sister;
						break;
					} else if(HydroNodeImported[NextNodeID].Parent == RootNodeID){
                        HydroRootImported.CurrentMaxLevel = CurrentMaxLevel;

                        int NumberofLeaves = HydroNodeImported[RootNodeID].NumberofLeaves;
                        for(int k=0;k<NumberofLeaves;k++){
                            int leaf =  HydroRootImported.Leaves[k];
                            NBImportedCache[k].Pos[0] = CachedData[leaf].Pos[0];
                            NBImportedCache[k].Pos[1] = CachedData[leaf].Pos[1];
                            NBImportedCache[k].Pos[2] = CachedData[leaf].Pos[2];
                            NBImportedCache[k].Kernel = CachedData[leaf].Kernel;
                            NBImportedCache[k].Active = CachedData[leaf].Active;
                            NBImportedCache[k].Leaf = leaf;
                        }
                        HydroRootImported.NumberofNodes = NumberofNodes + 1;
						return;
					}
                    NextNodeID = HydroNodeImported[NextNodeID].Parent;
				}
			}
			continue;
		}

		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

		for(int i=0;i<HydroNodeImported[CurrentNodeID].NumberofLeaves;i++){
            int leaf = HydroRootImported.Leaves[HydroNodeImported[CurrentNodeID].Leaves+i];
			int subindex = 0;

            for(int k=0;k<3;k++)
                if(HydroNodeImported[CurrentNodeID].Pos[k] <= CachedData[leaf].Pos[k])
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
                if(1.1*NumberofNodes >= HydroRootImported.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(ForAngelsShare*NumberofNodes);
                    HydroNodeImported = realloc(HydroNodeImported,sizeof(struct StructHydroNode)*NumberofAllocatedNodes);
                    HydroRootImported.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                HydroNodeImported[ChildNodeID].Next = NONE;
                HydroNodeImported[ChildNodeID].Parent = NONE;
                HydroNodeImported[ChildNodeID].Sister = NONE;
                HydroNodeImported[ChildNodeID].Children = NONE;
                // HydroNodeImported[ChildNodeID].Traverse = NONE;
                // make node

                HydroNodeImported[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    HydroNodeImported[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    HydroNodeImported[ChildNodeID].Leaves = HydroNodeImported[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,HydroNodeImported[CurrentNodeID].Level+1);
                } else {
                    HydroNodeImported[BackwardNodeID].Sister = ChildNodeID;
                    HydroNodeImported[ChildNodeID].Leaves = HydroNodeImported[BackwardNodeID].Leaves + HydroNodeImported[BackwardNodeID].NumberofLeaves;
                }

                HydroRootImported.Leaves[HydroNodeImported[ChildNodeID].Leaves] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    HydroRootImported.Leaves[HydroNodeImported[ChildNodeID].Leaves+k] = 
                        sublist[HydroRootImported.Leaves[HydroNodeImported[ChildNodeID].Leaves+k-1]];
                }
                HydroNodeImported[ChildNodeID].NumberofLeaves = subnumber[i];
                HydroNodeImported[ChildNodeID].Level = HydroNodeImported[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    HydroNodeImported[ChildNodeID].Pos[k] = HydroNodeImported[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*HydroRootImported.Width*HydroRootImported.WidthFactor[HydroNodeImported[CurrentNodeID].Level];
            }
		}
        CurrentNodeID = NextNodeID;
	}
}

#if 0
static void BuildHydroTraversalLinkImported(void){

    int MarkerID[TreeMaxNodeLevel];
    int RootNodeID = 0;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        TraversalID[i] = NONE;
    TraversalID[RootNodeID] = 0;
    HydroNodeImported[RootNodeID].Traverse = NONE;

    for(int i=1;i<HydroRootImported.NumberofNodes;i++){
        HydroNodeImported[i].Traverse = NONE;

        if(TraversalID[HydroNodeImported[i].Level] != NONE){
            HydroNodeImported[MarkerID[HydroNodeImported[i].Level]].Traverse = i;
        } else {
            TraversalID[HydroNodeImported[i].Level] = i;
        }
        MarkerID[HydroNodeImported[i].Level] = i; 

        HydroNodeImported[i].Next = NextHydroNodeImported(i);
    }

	return;
}
#endif

static inline void CalcDistanceMaxKernelMaxActiveNumberImported(const int CurrentNodeID) __attribute__((always_inline));
static inline void CalcDistanceMaxKernelMaxActiveNumberImported(const int CurrentNodeID){

    double PosLeaf[3] = {HydroNodeImported[CurrentNodeID].Pos[0],
                         HydroNodeImported[CurrentNodeID].Pos[1],
                         HydroNodeImported[CurrentNodeID].Pos[2]};
    double distance,kernel;
    double distancemax = 0.e0;
    double kernelmax = 0.e0;

    //int NActives = 0;
    int NumberofLeaves = HydroNodeImported[CurrentNodeID].NumberofLeaves;
    int header = HydroNodeImported[CurrentNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
        distance = DISTANCE(PosLeaf,NBImportedCache[leaf].Pos);
        distancemax = fmax(distance,distancemax);
        kernel = 2.e0*NBImportedCache[leaf].Kernel;
        kernelmax = fmax(kernelmax,kernel+distance);
        //NActives += NBImportedCache[leaf].Active;
    }

    HydroNodeImported[CurrentNodeID].KernelMax = kernelmax;
    HydroNodeImported[CurrentNodeID].DistanceMax = distancemax;
    //HydroNodeImported[CurrentNodeID].NumberofActiveLeaves = NActives;
    HydroNodeImported[CurrentNodeID].NumberofActiveLeaves = NumberofLeaves;

    return;
}


#if 0
static void HydroNodeDataImplantImported(void){

    int RootNodeID = 0;
    double diag = DISTANCE(HydroNodeImported[RootNodeID].Pos,HydroNodeImported[HydroNodeImported[RootNodeID].Children].Pos);

    int ChildNodeID;
    for(int i=HydroRootImported.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(HydroNodeImported[CurrentNodeID].Children != NONE){

                ChildNodeID = HydroNodeImported[CurrentNodeID].Children;

                double width = diag*HydroRootImported.WidthFactor[i];
                double distancemax = 0.e0;
                double kernelmax = 0.e0;
                int NActives = 0;

                while(ChildNodeID != NONE){
                    kernelmax = fmax(kernelmax,HydroNodeImported[ChildNodeID].KernelMax);
                    distancemax = fmax(distancemax,HydroNodeImported[ChildNodeID].DistanceMax);
                    NActives += HydroNodeImported[ChildNodeID].NumberofActiveLeaves;
                    ChildNodeID = HydroNodeImported[ChildNodeID].Sister;
                }

                HydroNodeImported[CurrentNodeID].KernelMax = kernelmax + width;
                HydroNodeImported[CurrentNodeID].DistanceMax = distancemax + width;
                HydroNodeImported[CurrentNodeID].NumberofActiveLeaves = NActives;
                HydroNodeImported[CurrentNodeID].Leaves = HydroNodeImported[HydroNodeImported[CurrentNodeID].Children].Leaves;
            }else{ // No child node case.
                CalcDistanceMaxKernelMaxActiveNumberImported(CurrentNodeID); 
            }
            CurrentNodeID = HydroNodeImported[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}
#endif

static void HydroNodeDataImplantImportedNew(const int CurrentNodeID){

    double Width = DiagHydroImported*HydroRootImported.WidthFactor[HydroNodeImported[CurrentNodeID].Level];
    double DistanceMax = 0.e0;
    double KernelMax = 0.e0;

    if(HydroNodeImported[CurrentNodeID].Children == NONE){
        int NumberofLeaves = HydroNodeImported[CurrentNodeID].NumberofLeaves;
        int header = HydroNodeImported[CurrentNodeID].Leaves;
        for(int k=0;k<NumberofLeaves;k++){
            int leaf = header+k;
            double Distance = DISTANCE(HydroNodeImported[CurrentNodeID].Pos,NBImportedCache[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            KernelMax = fmax(KernelMax,2.e0*NBImportedCache[leaf].Kernel+Distance);
        }
        HydroNodeImported[CurrentNodeID].KernelMax = KernelMax;
        HydroNodeImported[CurrentNodeID].DistanceMax = DistanceMax;
        HydroNodeImported[CurrentNodeID].NumberofActiveLeaves = NumberofLeaves;
    } else {
        int ChildNodeID = HydroNodeImported[CurrentNodeID].Children;
        int NActives = 0;
        do{
            HydroNodeDataImplantImportedNew(ChildNodeID);

            KernelMax = fmax(KernelMax,HydroNodeImported[ChildNodeID].KernelMax);
            DistanceMax = fmax(DistanceMax,HydroNodeImported[ChildNodeID].DistanceMax);
            NActives += HydroNodeImported[ChildNodeID].NumberofActiveLeaves;
            ChildNodeID = HydroNodeImported[ChildNodeID].Sister;
        } while(ChildNodeID != NONE);

        HydroNodeImported[CurrentNodeID].KernelMax = KernelMax + Width;
        HydroNodeImported[CurrentNodeID].DistanceMax = DistanceMax + Width;
        HydroNodeImported[CurrentNodeID].NumberofActiveLeaves = NActives;
    }

    return;
}
