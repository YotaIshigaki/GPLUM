#include	"config.h"
#include	"PlantGravityTree.h"

void InitializeRootForGravity(void){

    int NloadGravity =(int)(ForAngelsShare*Pall.Ntotal);
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)Pall.Ntotal)/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = (Pall.Ntotal/CUBE(BaseNumberofChildren))*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));
    NloadGravity = MAX(NloadGravity,NAdditionUnit);
    NumberofAllocatedNodes = MAX(NumberofAllocatedNodes,NAdditionUnit);
    // dprintlmpi(NumberofAllocatedNodes);
    // dlprintlmpi(Pall.Ntotal);

    /* GravtityRoot */
    GravityRoot.NumberofAllocatedLeaves = NloadGravity;
	GravityRoot.Leaves = malloc(sizeof(int)*NloadGravity+1);
    GravityRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
    GravityRoot.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

	GravityRoot.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForGrav;
    GravityRoot.MaxLevel = TreeMaxNodeLevel;

    GravityRoot.OpeningAngle = TreeOpeningAngle;
    GravityRoot.NumberofLeavesInGroup = TreeNGroup;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        GravityRoot.WidthFactor[i] = pow(dw,(double)i);
    /* GravtityRoot */

    /* GravityNode */
    GravityNode = malloc(sizeof(struct StructGravityNode)*GravityRoot.NumberofAllocatedNodes);

    struct StructGravityNode TempGravityNode;
    memset(&TempGravityNode,0,sizeof(struct StructGravityNode));
    TempGravityNode.Next = TempGravityNode.Parent = TempGravityNode.Children = 
    TempGravityNode.Sister = NONE;
    //TempGravityNode.Sister = TempGravityNode.Traverse = NONE;

    for(int i=0;i<GravityRoot.NumberofAllocatedNodes;i++)
        GravityNode[i] = TempGravityNode;
    /* GravityNode */

    /* GravityCache */
    GravityCache = malloc(sizeof(struct StructGravityCache)*NloadGravity+1);
    GravityAccPotCache = malloc(sizeof(struct StructGravityAccPotCache)*NloadGravity+1);
    /* GravityCache */
    
    int NProcs = MPIGetNumProcs();
    EdgesForGravity = malloc(sizeof(struct StructEdges)*NProcs);

	return;
}

static void MakeGravityRoot(void);
static void BuildGravityTree(void);
static void BuildGravityTraversalLink(void);
static int NextGravityNode(const int NodeID);
static int TraversalID[TreeMaxNodeLevel+1];
static void GravityNodeDataImplant(void);
//static void GravityNodeDataImplantPartial(void);
static void UpdateDomainEdgesGravityParticles(void);
//static void UpdateDomainEdgesGravityParticlesPartial(void);

static int NumberofAllocatedLeaves = 0;
static int *sublist; 

static struct StructCachedData{
    double Pos[3];
    double Mass;
    double Eps;
    bool Active;
} *CachedData;

static void GravityTreePreprocessing(void){

    if(Pall.Ntotal>GravityRoot.NumberofAllocatedLeaves){
        GravityRoot.NumberofAllocatedLeaves = (int)(MAX(ForAngelsShare*Pall.Ntotal,NAdditionUnit));
        free(GravityRoot.Leaves);
        free(GravityCache);
        free(GravityAccPotCache);

        GravityRoot.Leaves = malloc(sizeof(int)*GravityRoot.NumberofAllocatedLeaves);
        GravityCache = malloc(sizeof(struct StructGravityCache)*GravityRoot.NumberofAllocatedLeaves);
        GravityAccPotCache = malloc(sizeof(struct StructGravityAccPotCache)*GravityRoot.NumberofAllocatedLeaves);
    }

    if(NumberofAllocatedLeaves < Pall.Ntotal){
        if(NumberofAllocatedLeaves > 0){
            free(sublist);
            free(CachedData);
        }
        NumberofAllocatedLeaves = (int)(MAX(ForAngelsShare*Pall.Ntotal,NAdditionUnit));
	    sublist = malloc(sizeof(int)*NumberofAllocatedLeaves);
        CachedData = malloc(sizeof(struct StructCachedData)*NumberofAllocatedLeaves);
    }

    //dprintlmpi(GravityRoot.NumberofAllocatedLeaves);
    //dprintlmpi(NumberofAllocatedLeaves);
    //dprintlmpi(Pall.Ntotal);
#if 0
    fprintf(stderr,"[%04d] FUNCTION : %s, Allocate memory size = %g MB\n",MPIGetMyID(),__FUNCTION__,
            (sizeof(int)*GravityRoot.NumberofAllocatedLeaves+
             sizeof(struct StructGravityCache)*GravityRoot.NumberofAllocatedLeaves+
             sizeof(struct StructGravityAccPotCache)*GravityRoot.NumberofAllocatedLeaves+
             sizeof(int)*NumberofAllocatedLeaves+
             sizeof(struct StructCachedData)*NumberofAllocatedLeaves)/(1024.0*1024.0));
#endif

    double SofteningFactor = Pall.AdaptiveSofteningFactor;

    for(int i=0;i<Pall.Ntotal;i++){
        CachedData[i].Pos[0] = Pbody[i]->PosP[0];
        CachedData[i].Pos[1] = Pbody[i]->PosP[1];
        CachedData[i].Pos[2] = Pbody[i]->PosP[2];
        CachedData[i].Mass   = Pbody[i]->Mass;
        CachedData[i].Eps    = SofteningFactor*Pbody[i]->Eps;
        CachedData[i].Active = Pbody[i]->Active;
    }

    return ;
}

#if 0
static void CalcTreeNodeEPS(const int CurrentNodeID){

    double Mass = 0.e0;
    double Eps2 = 0.e0;
    if(GravityNode[CurrentNodeID].Children == NONE){
        int Number_of_leaf = GravityNode[CurrentNodeID].NumberofLeaves;
        int header = GravityNode[CurrentNodeID].Leaves;
        for(int k=0;k<Number_of_leaf;k++){
            int leaf = header+k;
            Mass += GravityCache[leaf].Mass;
            Eps2 += GravityCache[leaf].Mass*SQ(GravityCache[leaf].Eps);
        }
    } else {
        int ChildNodeID = GravityNode[CurrentNodeID].Children;
        while(ChildNodeID != NONE){
            CalcTreeNodeEPS(ChildNodeID);
            Mass += GravityNode[ChildNodeID].Mass;
            Eps2 += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].Eps2;
            ChildNodeID = GravityNode[ChildNodeID].Sister;
        }
    }
    GravityNode[CurrentNodeID].Eps2 = Eps2/Mass;
    return ;
}
#endif

static double Diag;
static void GravityNodeDataImplantNew(const int CurrentNodeID){
    double Width = Diag*GravityRoot.WidthFactor[GravityNode[CurrentNodeID].Level];

    int NActives = 0;
    double DistanceMax = 0.e0;
    double Mass = 0.e0;
    double COM[3] = {0.e0,0.e0,0.e0};
    double Eps2 = 0.e0;
    double EpsMax,EpsMin;
    if(GravityNode[CurrentNodeID].Children == NONE){
        int Number_of_leaf = GravityNode[CurrentNodeID].NumberofLeaves;
        int header = GravityNode[CurrentNodeID].Leaves;
#ifdef USE_SYMMETRIZED_SOFTENING
        EpsMax = EpsMin = GravityCache[header].Eps;
#endif // USE_SYMMETRIZED_SOFTENING
        for(int k=0;k<Number_of_leaf;k++){
            int leaf = header+k;
            double Distance = DISTANCE(GravityNode[CurrentNodeID].Pos,GravityCache[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            Mass += GravityCache[leaf].Mass;
            COM[0] += GravityCache[leaf].Mass*GravityCache[leaf].Pos[0];
            COM[1] += GravityCache[leaf].Mass*GravityCache[leaf].Pos[1];
            COM[2] += GravityCache[leaf].Mass*GravityCache[leaf].Pos[2];
            NActives += GravityCache[leaf].Active;
#ifdef USE_SYMMETRIZED_SOFTENING
            Eps2 += GravityCache[leaf].Mass*SQ(GravityCache[leaf].Eps);
            EpsMin = fmin(EpsMin,GravityCache[leaf].Eps);
            EpsMax = fmax(EpsMax,GravityCache[leaf].Eps);
#endif // USE_SYMMETRIZED_SOFTENING
        }
        Width = 0.e0;
    } else {
        bool first = true;
        int ChildNodeID = GravityNode[CurrentNodeID].Children;
        do{
            GravityNodeDataImplantNew(ChildNodeID);

            DistanceMax = fmax(DistanceMax,GravityNode[ChildNodeID].DistanceMax);
            Mass += GravityNode[ChildNodeID].Mass;
            COM[0] += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].COM[0];
            COM[1] += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].COM[1];
            COM[2] += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].COM[2];
            NActives += GravityNode[ChildNodeID].NumberofActiveLeaves;
#ifdef USE_SYMMETRIZED_SOFTENING
            Eps2 += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].Eps2;
            if(first == true){
                EpsMin = GravityNode[ChildNodeID].EpsMin;
                EpsMax = GravityNode[ChildNodeID].EpsMax;
                first = false;
            } else {
                EpsMin = fmin(EpsMin,GravityNode[ChildNodeID].EpsMin);
                EpsMax = fmax(EpsMax,GravityNode[ChildNodeID].EpsMax);
            }
#endif // USE_SYMMETRIZED_SOFTENING
            ChildNodeID = GravityNode[ChildNodeID].Sister;
        } while(ChildNodeID != NONE);
    }
    double InvMass = 1.e0/Mass;
    GravityNode[CurrentNodeID].COM[0] = COM[0]*InvMass;
    GravityNode[CurrentNodeID].COM[1] = COM[1]*InvMass;
    GravityNode[CurrentNodeID].COM[2] = COM[2]*InvMass;

    GravityNode[CurrentNodeID].DistanceMax = DistanceMax + Width;
    GravityNode[CurrentNodeID].Mass = Mass;
    GravityNode[CurrentNodeID].NumberofActiveLeaves = NActives;
#ifdef USE_SYMMETRIZED_SOFTENING
    GravityNode[CurrentNodeID].Eps2 = Eps2*InvMass;
    GravityNode[CurrentNodeID].EpsMax = EpsMax;
    GravityNode[CurrentNodeID].EpsMin = EpsMin;
#endif // USE_SYMMETRIZED_SOFTENING
    return ;
}



void PlantGravityTree(void){

    double TimingResultThisRoutine = GetElapsedTime();

    GravityTreePreprocessing();
	MakeGravityRoot();
    if(GravityRoot.NumberofLeaves > 0){
        BuildGravityTree();
        Diag = DISTANCE(GravityNode[0].Pos,GravityNode[GravityNode[0].Children].Pos);
        GravityNodeDataImplantNew(0);

        for(int i=1;i<GravityRoot.NumberofNodes;i++){
            GravityNode[i].Next = NextGravityNode(i);
        }
    }

    // double Com = GetElapsedTime();
    // UpdateDomainEdgesGravityParticles();
    // fprintf(stderr,"Com calculation time = %g[sec]\n",GetElapsedTime()-Com);

    TimingResults.GravityTreeThisStep += GetElapsedTime()-TimingResultThisRoutine;

	return;
}

#if 0
void PlantGravityTreeOld(void){

    double TimingResultThisRoutine = GetElapsedTime();

    double Tpre = GetElapsedTime();
    GravityTreePreprocessing();
    fprintf(stderr,"Tpre calculation time = %g[sec]\n",GetElapsedTime()-Tpre);

    double Tmk = GetElapsedTime();
	MakeGravityRoot();
    fprintf(stderr,"Tmk calculation time = %g[sec]\n",GetElapsedTime()-Tmk);

    double Tbl = GetElapsedTime();
    BuildGravityTree();
    fprintf(stderr,"Tbl calculation time = %g[sec]\n",GetElapsedTime()-Tbl);


    double Ttr = GetElapsedTime();
    BuildGravityTraversalLink();
    fprintf(stderr,"Ttr calculation time = %g[sec]\n",GetElapsedTime()-Ttr);

    double Timp = GetElapsedTime();
    GravityNodeDataImplant();
    fprintf(stderr,"Imp calculation time = %g[sec]\n",GetElapsedTime()-Timp);

    // double Com = GetElapsedTime();
    // UpdateDomainEdgesGravityParticles();
    // fprintf(stderr,"Com calculation time = %g[sec]\n",GetElapsedTime()-Com);

    TimingResults.GravityTreeThisStep += GetElapsedTime()-TimingResultThisRoutine;
    fprintf(stderr,"Tree calculation time = %g[sec]\n",GetElapsedTime()-TimingResultThisRoutine);

	return;
}
#endif

#if 0
void PlantGravityTreePartial(void){

    assert(Pall.Ntotal == GravityRoot.NumberofLeaves);

    static int RootNodeID = 0; 
    int Number_of_leaf = GravityNode[RootNodeID].NumberofLeaves;
    for(int k=0;k<Number_of_leaf;k++){
        int leaf =  GravityRoot.Leaves[k];
        GravityCache[k].Mass = Pbody[leaf]->Mass;
    }

    GravityNodeDataImplantPartial();
    UpdateDomainEdgesGravityParticlesPartial();

	return;
}
#endif

static void MakeGravityRoot(void){

	double min[3],max[3];
    int RootNodeID = 0;
    double SofteningFactor = Pall.AdaptiveSofteningFactor;
    struct StructEdges TempEdges;

    if(Pall.Ntotal > 0){
        for(int k=0;k<3;k++){
            max[k] = CachedData[0].Pos[k];
            min[k] = CachedData[0].Pos[k];
            TempEdges.PosMax[k] = CachedData[0].Pos[k]+CachedData[0].Eps;
            TempEdges.PosMin[k] = CachedData[0].Pos[k]-CachedData[0].Eps;
        }
        for(int i=1;i<Pall.Ntotal;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(CachedData[i].Pos[k],max[k]);
                min[k] = fmin(CachedData[i].Pos[k],min[k]);
                TempEdges.PosMax[k] = fmax(CachedData[i].Pos[k]+CachedData[i].Eps,TempEdges.PosMax[k]);
                TempEdges.PosMin[k] = fmin(CachedData[i].Pos[k]-CachedData[i].Eps,TempEdges.PosMin[k]);
            }
        }
    } else {
        for(int k=0;k<3;k++)
            max[k] = min[k] = TempEdges.PosMax[k] = TempEdges.PosMin[k] = 0.e0;
    }

    int NProcs = MPIGetNumProcs();
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(&TempEdges,sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PLANTTREE_EXTENSITY+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(EdgesForGravity+CommunicationTable[i].RecvRank,
                sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
    }

    /* 
    MPI_Status  mpi_status;
    EdgesForGravity[MPIGetMyID()] = TempEdges;
    for(int i=0;i<MPIGetNumProcs()-1;i++){
        MPI_Sendrecv(&TempEdges,sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PLANTTREE_EXTENSITY,
            EdgesForGravity+CommunicationTable[i].RecvRank,
                    sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY,
            MPI_COMM_WORLD,&mpi_status);
    }
    */

    for(int k=0;k<3;k++){
        GravityRoot.PosMax[k] = max[k];
        GravityRoot.PosMin[k] = min[k];
    }

    double MaxMin[6] = {max[0],max[1],max[2],-min[0],-min[1],-min[2]};
    double GlobalMaxMin[6];
    MPI_Allreduce(MaxMin,GlobalMaxMin,6,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    max[0] = +GlobalMaxMin[0]; max[1] = +GlobalMaxMin[1]; max[2] = +GlobalMaxMin[2];
    min[0] = -GlobalMaxMin[3]; min[1] = -GlobalMaxMin[4]; min[2] = -GlobalMaxMin[5];

    double WidthMax = 0.e0;
    for(int k=0;k<3;k++){
        GravityNode[RootNodeID].Pos[k] = 0.5*(max[k] + min[k]);
        WidthMax = fmax(WidthMax,max[k]-min[k]);
    }
    GravityRoot.Width = WidthMax;

    GravityNode[RootNodeID].Next = NONE;
    GravityNode[RootNodeID].Parent = NONE;
    GravityNode[RootNodeID].Sister = NONE;
    GravityNode[RootNodeID].Children = NONE;
    //GravityNode[RootNodeID].Traverse = NONE;

    GravityNode[RootNodeID].NumberofLeaves = Pall.Ntotal;

    GravityNode[RootNodeID].Level = 0;
    GravityNode[RootNodeID].Leaves = 0;

    for(int i=0;i<Pall.Ntotal;i++)
        GravityRoot.Leaves[i] = i;

    GravityRoot.NumberofLeaves = Pall.Ntotal;

    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);

	return ;
}

static int NextGravityNode(const int NodeID){

    int CurrentNodeID = NodeID;

    if(GravityNode[CurrentNodeID].Sister != NONE){
        CurrentNodeID = GravityNode[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(GravityNode[GravityNode[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = GravityNode[GravityNode[NextNodeID].Parent].Sister;
                break;
            } else if(GravityNode[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = GravityNode[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}


static inline bool GravityNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber) __attribute__((always_inline));
static inline bool GravityNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber){

	if( (GravityNode[CurrentNodeID].NumberofLeaves <= CriticalNumber) || GravityNode[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}

static void BuildGravityTree(void){

    int NumberofNodes = 0; 
    int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 

    int NumberofNodeCreationLimit = GravityRoot.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
	while(1){

		if(GravityNodeSeparationCriterion(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(GravityNode[CurrentNodeID].Sister != NONE){
				CurrentNodeID = GravityNode[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(GravityNode[GravityNode[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = GravityNode[GravityNode[NextNodeID].Parent].Sister;
						break;
					} else if(GravityNode[NextNodeID].Parent == RootNodeID){
                        GravityRoot.CurrentMaxLevel = CurrentMaxLevel;

                        int NumberofLeaves = GravityNode[RootNodeID].NumberofLeaves;
                        for(int k=0;k<NumberofLeaves;k++){
                            int leaf =  GravityRoot.Leaves[k];
                            GravityCache[k].Pos[0] = CachedData[leaf].Pos[0];
                            GravityCache[k].Pos[1] = CachedData[leaf].Pos[1];
                            GravityCache[k].Pos[2] = CachedData[leaf].Pos[2];
                            GravityCache[k].Mass   = CachedData[leaf].Mass;
                            GravityCache[k].Eps    = CachedData[leaf].Eps;
                            GravityCache[k].Active = CachedData[leaf].Active;
                            GravityCache[k].Leaf = leaf;
                        }
                        GravityRoot.NumberofNodes = NumberofNodes + 1;

						return;
					}
                    NextNodeID = GravityNode[NextNodeID].Parent;
				}
			}
			continue;
		}

		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

        int NumberofLeaves = GravityNode[CurrentNodeID].NumberofLeaves;
		for(int i=0;i<NumberofLeaves;i++){
            int leaf = GravityRoot.Leaves[GravityNode[CurrentNodeID].Leaves+i];
			int subindex = 0;

            for(int k=0;k<3;k++)
                if(GravityNode[CurrentNodeID].Pos[k] <= CachedData[leaf].Pos[k])
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
                //if(GravityNode[CurrentNodeID].Level == 20){
                    //fprintf(stderr,"Level = %d, Member = %d\n",GravityNode[CurrentNodeID].Level+1,subnumber[i]);
                //}
				BackwardNodeID = ChildNodeID; 
                // make node
                NumberofNodes ++;
                ChildNodeID = NumberofNodes;
                if(NumberofNodes >= GravityRoot.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(MAX(ForAngelsShare*NumberofNodes,NAdditionUnit));
                    GravityNode = realloc(GravityNode,sizeof(struct StructGravityNode)*NumberofAllocatedNodes);
                    GravityRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                GravityNode[ChildNodeID].Next = NONE;
                GravityNode[ChildNodeID].Parent = NONE;
                GravityNode[ChildNodeID].Sister = NONE;
                GravityNode[ChildNodeID].Children = NONE;
                // GravityNode[ChildNodeID].Traverse = NONE;
                // make node

                GravityNode[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    GravityNode[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    GravityNode[ChildNodeID].Leaves = GravityNode[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,GravityNode[CurrentNodeID].Level+1);
                } else {
                    GravityNode[BackwardNodeID].Sister = ChildNodeID;
                    GravityNode[ChildNodeID].Leaves = 
                        GravityNode[BackwardNodeID].Leaves + GravityNode[BackwardNodeID].NumberofLeaves;
                }

                GravityRoot.Leaves[GravityNode[ChildNodeID].Leaves] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    GravityRoot.Leaves[GravityNode[ChildNodeID].Leaves+k] = 
                        sublist[GravityRoot.Leaves[GravityNode[ChildNodeID].Leaves+k-1]];
                }
                GravityNode[ChildNodeID].NumberofLeaves = subnumber[i];
                GravityNode[ChildNodeID].Level = GravityNode[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    GravityNode[ChildNodeID].Pos[k] = GravityNode[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*GravityRoot.Width*GravityRoot.WidthFactor[GravityNode[CurrentNodeID].Level];
            }
		}
        CurrentNodeID = NextNodeID;

	}
}

#if 0
static void BuildGravityTraversalLink(void){

    int MarkerID[TreeMaxNodeLevel];
    int RootNodeID = 0;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        TraversalID[i] = NONE;
    TraversalID[RootNodeID] = 0;
    GravityNode[RootNodeID].Traverse = NONE;

    for(int i=1;i<GravityRoot.NumberofNodes;i++){
        GravityNode[i].Traverse = NONE;

        if(TraversalID[GravityNode[i].Level] != NONE){
            GravityNode[MarkerID[GravityNode[i].Level]].Traverse = i;
        } else {
            TraversalID[GravityNode[i].Level] = i;
        }
        MarkerID[GravityNode[i].Level] = i; 

        GravityNode[i].Next = NextGravityNode(i);
    }

	return;
}
#endif

static inline void CalcDistanceMaxMassNumberofActiveLeaves(const int NodeID) __attribute__((always_inline));
static inline void CalcDistanceMaxMassNumberofActiveLeaves(const int NodeID){

    double PosLeaf[3] = {GravityNode[NodeID].Pos[0],GravityNode[NodeID].Pos[1],GravityNode[NodeID].Pos[2]};
    double distance;
    double distancemax = 0.e0;
    double mass = 0.e0;
    double COM[3] = {0.e0,0.e0,0.e0};

    int NActives = 0;
    int Number_of_leaf = GravityNode[NodeID].NumberofLeaves;
    int header = GravityNode[NodeID].Leaves;
    for(int k=0;k<Number_of_leaf;k++){
        int leaf = header+k;
        distance = DISTANCE(PosLeaf,GravityCache[leaf].Pos);
        distancemax = fmax(distance,distancemax);
        mass += GravityCache[leaf].Mass;
        COM[0] += GravityCache[leaf].Mass*GravityCache[leaf].Pos[0];
        COM[1] += GravityCache[leaf].Mass*GravityCache[leaf].Pos[1];
        COM[2] += GravityCache[leaf].Mass*GravityCache[leaf].Pos[2];
        NActives += GravityCache[leaf].Active;
    }
    double imass = 1.e0/mass;
    GravityNode[NodeID].COM[0] = COM[0]*imass;
    GravityNode[NodeID].COM[1] = COM[1]*imass;
    GravityNode[NodeID].COM[2] = COM[2]*imass;

    GravityNode[NodeID].Mass = mass;
    GravityNode[NodeID].DistanceMax = distancemax;
    GravityNode[NodeID].NumberofActiveLeaves = NActives;

    return;
}


#if 0
static void GravityNodeDataImplant(void){

    int RootNodeID = 0;
    double diag = DISTANCE(GravityNode[RootNodeID].Pos,GravityNode[GravityNode[RootNodeID].Children].Pos);

    int ChildNodeID;
    for(int i=GravityRoot.CurrentMaxLevel;0<=i;i--){
        int CurrentNodeID = TraversalID[i];
        do{
            if(GravityNode[CurrentNodeID].Children != NONE){

                ChildNodeID = GravityNode[CurrentNodeID].Children;

                double width = diag*GravityRoot.WidthFactor[i];
                double distancemax = 0.e0;
                double mass = 0.e0;
                double COM[3] = {0.e0,0.e0,0.e0};

                int NActives = 0;
                while(ChildNodeID != NONE){
                    distancemax = fmax(distancemax,GravityNode[ChildNodeID].DistanceMax);
                    mass += GravityNode[ChildNodeID].Mass;
                    COM[0] += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].COM[0];
                    COM[1] += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].COM[1];
                    COM[2] += GravityNode[ChildNodeID].Mass*GravityNode[ChildNodeID].COM[2];
                    NActives += GravityNode[ChildNodeID].NumberofActiveLeaves;

                    ChildNodeID = GravityNode[ChildNodeID].Sister;
                }
                double imass = 1.e0/mass;
                GravityNode[CurrentNodeID].COM[0] = COM[0]*imass;
                GravityNode[CurrentNodeID].COM[1] = COM[1]*imass;
                GravityNode[CurrentNodeID].COM[2] = COM[2]*imass;

                GravityNode[CurrentNodeID].DistanceMax = distancemax + width;
                GravityNode[CurrentNodeID].Mass = mass;
                GravityNode[CurrentNodeID].NumberofActiveLeaves = NActives;
                GravityNode[CurrentNodeID].Leaves = GravityNode[GravityNode[CurrentNodeID].Children].Leaves;
            }else{ // No child node case.
                CalcDistanceMaxMassNumberofActiveLeaves(CurrentNodeID); 
            }
            CurrentNodeID = GravityNode[CurrentNodeID].Traverse;
        } while(CurrentNodeID!=NONE);
    }

    return;
}
#endif

static void UpdateDomainEdgesGravityParticles(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double SofteningFactor = Pall.AdaptiveSofteningFactor;
    struct StructEdges TempEdges;

    double max[3],min[3];

    for(int k=0;k<3;k++){
        max[k] = GravityCache[0].Pos[k]+SofteningFactor*Pbody[0]->Eps;
        min[k] = GravityCache[0].Pos[k]-SofteningFactor*Pbody[0]->Eps;
    }
    for(int i=1;i<Pall.Ntotal;i++){
        for(int k=0;k<3;k++){
            max[k] = fmax(GravityCache[i].Pos[k]+SofteningFactor*Pbody[i]->Eps,max[k]);
            min[k] = fmin(GravityCache[i].Pos[k]-SofteningFactor*Pbody[i]->Eps,min[k]);
        }
    }
    TempEdges.PosMax[0] = max[0]; TempEdges.PosMin[0] = min[0];
    TempEdges.PosMax[1] = max[1]; TempEdges.PosMin[1] = min[1];
    TempEdges.PosMax[2] = max[2]; TempEdges.PosMin[2] = min[2];

    MPI_Status  mpi_status;

    EdgesForGravity[MyID] = TempEdges;

    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(&TempEdges,sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PLANTTREE_EXTENSITY,
            EdgesForGravity+CommunicationTable[i].RecvRank,
                    sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY,
            MPI_COMM_WORLD,&mpi_status);
    }

    return;
}
