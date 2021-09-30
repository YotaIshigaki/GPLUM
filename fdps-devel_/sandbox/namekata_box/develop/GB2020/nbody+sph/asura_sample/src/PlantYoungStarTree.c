#include "config.h"
#include "PlantYoungStarTree.h"

static bool FirstCallPlantYoungStarTree = true;
void InitializeRootForYoungStarTree(void){

    int NloadStars =(int)(ForAngelsShare*MAX(Pall.Nstars,NAdditionUnit));
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)Pall.Nstars)/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = (Pall.Nstars/CUBE(BaseNumberofChildren))*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));
    NloadStars = MAX(NloadStars,NAdditionUnit);
    NumberofAllocatedNodes = MAX(NumberofAllocatedNodes,NAdditionUnit);

    /* YoungStarRoot */
    YoungStarRoot.NumberofAllocatedLeaves = NloadStars;
	YoungStarRoot.Leaves = malloc(sizeof(int)*NloadStars+1);
    YoungStarRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
    YoungStarRoot.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

	YoungStarRoot.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForGrav;
    YoungStarRoot.MaxLevel = TreeMaxNodeLevel;

    YoungStarRoot.OpeningAngle = FUVFEEDBACK_THETA;
    YoungStarRoot.NumberofLeavesInGroup = FUVFEEDBACK_NGROUP;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        YoungStarRoot.WidthFactor[i] = pow(dw,(double)i);
    /* YoungStarRoot */

    /* YoungStarode */
    YoungStarNode = malloc(sizeof(struct StructYoungStarNode)*YoungStarRoot.NumberofAllocatedNodes);

    struct StructYoungStarNode TempYoungStarNode;
    memset(&TempYoungStarNode,0,sizeof(struct StructYoungStarNode));
    TempYoungStarNode.Next = TempYoungStarNode.Parent = TempYoungStarNode.Children = 
    TempYoungStarNode.Sister = NONE;

    for(int i=0;i<YoungStarRoot.NumberofAllocatedNodes;i++)
        YoungStarNode[i] = TempYoungStarNode;
    /* YoungStarNode */

    /* YoungStarCache */
    YoungStarCache = malloc(sizeof(struct StructYoungStarCache)*NloadStars+1);
    YoungStarResultCache = malloc(sizeof(struct StructYoungStarResultCache)*NloadStars+1);
    /* YoungStarCache */
    
    int NProcs = MPIGetNumProcs();
    EdgesForYoungStars = malloc(sizeof(struct StructEdges)*NProcs);

	return;
}

static int NumberofAllocatedLeaves = 0;
static int *sublist; 

static struct StructCachedData{
    double Pos[3];
    double LFUV;
    double Eps;
    bool Active;
} *CachedData;

static int YoungStarTreePreprocessing(void){

    if(Pall.Nstars>YoungStarRoot.NumberofAllocatedLeaves){
        YoungStarRoot.NumberofAllocatedLeaves = (int)(MAX(ForAngelsShare*Pall.Nstars,NAdditionUnit));
        free(YoungStarRoot.Leaves);
        free(YoungStarCache);
        free(YoungStarResultCache);

        YoungStarRoot.Leaves = malloc(sizeof(int)*YoungStarRoot.NumberofAllocatedLeaves);
        YoungStarCache = malloc(sizeof(struct StructYoungStarCache)*YoungStarRoot.NumberofAllocatedLeaves);
        YoungStarResultCache = malloc(sizeof(struct StructYoungStarResultCache)*YoungStarRoot.NumberofAllocatedLeaves);
    }

    if(NumberofAllocatedLeaves < Pall.Nstars){
        if(NumberofAllocatedLeaves > 0){
            free(sublist);
            free(CachedData);
        }
        NumberofAllocatedLeaves = (int)(MAX(ForAngelsShare*Pall.Nstars,NAdditionUnit));
	    sublist = malloc(sizeof(int)*NumberofAllocatedLeaves);
        CachedData = malloc(sizeof(struct StructCachedData)*NumberofAllocatedLeaves);
    }

    double SofteningFactor = Pall.AdaptiveSofteningFactor;

    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if((Pall.TCurrent-Pstar[i]->FormationTime)*Pall.UnitTime < YOUNGSTARTREE_AGELIMIT){
            CachedData[counter].Pos[0] = PstarBody(i)->PosP[0];
            CachedData[counter].Pos[1] = PstarBody(i)->PosP[1];
            CachedData[counter].Pos[2] = PstarBody(i)->PosP[2];
            CachedData[counter].LFUV   = Pstar[i]->LFUV;
            CachedData[counter].Eps    = SofteningFactor*PstarBody(i)->Eps;
            CachedData[counter].Active = Pbody[i]->Active;
            counter ++;
        }
    }

    return counter;
}

static void MakeYoungStarRoot(const int NYoungStars){

	double min[3],max[3];
    int RootNodeID = 0;
    double SofteningFactor = Pall.AdaptiveSofteningFactor;
    struct StructEdges TempEdges;

    //dprintlmpi(NYoungStars);

    if(NYoungStars > 0){
        for(int k=0;k<3;k++){
            max[k] = CachedData[0].Pos[k];
            min[k] = CachedData[0].Pos[k];
            TempEdges.PosMax[k] = CachedData[0].Pos[k]+CachedData[0].Eps;
            TempEdges.PosMin[k] = CachedData[0].Pos[k]-CachedData[0].Eps;
        }
        for(int i=1;i<NYoungStars;i++){
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
        MPI_Irecv(EdgesForYoungStars+CommunicationTable[i].RecvRank,
                sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
    }


    for(int k=0;k<3;k++){
        YoungStarRoot.PosMax[k] = max[k];
        YoungStarRoot.PosMin[k] = min[k];
        // fprintf(stderr,"%g %g\n",max[k],min[k]);
    }

    double MaxMin[6] = {max[0],max[1],max[2],-min[0],-min[1],-min[2]};
    double GlobalMaxMin[6];
    MPI_Allreduce(MaxMin,GlobalMaxMin,6,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    max[0] = +GlobalMaxMin[0]; max[1] = +GlobalMaxMin[1]; max[2] = +GlobalMaxMin[2];
    min[0] = -GlobalMaxMin[3]; min[1] = -GlobalMaxMin[4]; min[2] = -GlobalMaxMin[5];

    double WidthMax = 0.e0;
    for(int k=0;k<3;k++){
        YoungStarNode[RootNodeID].Pos[k] = 0.5*(max[k] + min[k]);
        WidthMax = fmax(WidthMax,max[k]-min[k]);
    }
    YoungStarRoot.Width = WidthMax;

    YoungStarNode[RootNodeID].Next = NONE;
    YoungStarNode[RootNodeID].Parent = NONE;
    YoungStarNode[RootNodeID].Sister = NONE;
    YoungStarNode[RootNodeID].Children = NONE;
    YoungStarNode[RootNodeID].NumberofLeaves = NYoungStars;

    YoungStarNode[RootNodeID].Level = 0;
    YoungStarNode[RootNodeID].Leaves = 0;

    for(int i=0;i<NYoungStars;i++)
        YoungStarRoot.Leaves[i] = i;

    YoungStarRoot.NumberofLeaves = NYoungStars;

    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);

	return ;
}


static inline bool __attribute__((always_inline)) YoungStarNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber){

	if( (YoungStarNode[CurrentNodeID].NumberofLeaves <= CriticalNumber) || YoungStarNode[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}


static void BuildYoungStarTree(void){

    int NumberofNodes = 0; 
    int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 

    int NumberofNodeCreationLimit = YoungStarRoot.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
	while(1){
		if(YoungStarNodeSeparationCriterion(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(YoungStarNode[CurrentNodeID].Sister != NONE){
				CurrentNodeID = YoungStarNode[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(YoungStarNode[YoungStarNode[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = YoungStarNode[YoungStarNode[NextNodeID].Parent].Sister;
						break;
					} else if(YoungStarNode[NextNodeID].Parent == RootNodeID){
                        YoungStarRoot.CurrentMaxLevel = CurrentMaxLevel;

                        int NumberofLeaves = YoungStarNode[RootNodeID].NumberofLeaves;
                        for(int k=0;k<NumberofLeaves;k++){
                            int leaf =  YoungStarRoot.Leaves[k];
                            YoungStarCache[k].Pos[0] = CachedData[leaf].Pos[0];
                            YoungStarCache[k].Pos[1] = CachedData[leaf].Pos[1];
                            YoungStarCache[k].Pos[2] = CachedData[leaf].Pos[2];
                            YoungStarCache[k].LFUV   = CachedData[leaf].LFUV;
                            YoungStarCache[k].Eps    = CachedData[leaf].Eps;
                            YoungStarCache[k].Active = CachedData[leaf].Active;
                            YoungStarCache[k].Leaf = leaf;
                        }
                        YoungStarRoot.NumberofNodes = NumberofNodes + 1;

						return;
					}
                    NextNodeID = YoungStarNode[NextNodeID].Parent;
				}
			}
			continue;
		}

		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

        int NumberofLeaves = YoungStarNode[CurrentNodeID].NumberofLeaves;
		for(int i=0;i<NumberofLeaves;i++){
            int leaf = YoungStarRoot.Leaves[YoungStarNode[CurrentNodeID].Leaves+i];
			int subindex = 0;

            for(int k=0;k<3;k++)
                if(YoungStarNode[CurrentNodeID].Pos[k] <= CachedData[leaf].Pos[k])
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
                NumberofNodes ++;
                ChildNodeID = NumberofNodes;
                if(NumberofNodes >= YoungStarRoot.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(MAX(ForAngelsShare*NumberofNodes,NAdditionUnit));
                    YoungStarNode = realloc(YoungStarNode,sizeof(struct StructYoungStarNode)*NumberofAllocatedNodes);
                    YoungStarRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                YoungStarNode[ChildNodeID].Next = NONE;
                YoungStarNode[ChildNodeID].Parent = NONE;
                YoungStarNode[ChildNodeID].Sister = NONE;
                YoungStarNode[ChildNodeID].Children = NONE;

                YoungStarNode[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    YoungStarNode[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    YoungStarNode[ChildNodeID].Leaves = YoungStarNode[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,YoungStarNode[CurrentNodeID].Level+1);
                } else {
                    YoungStarNode[BackwardNodeID].Sister = ChildNodeID;
                    YoungStarNode[ChildNodeID].Leaves = 
                        YoungStarNode[BackwardNodeID].Leaves + YoungStarNode[BackwardNodeID].NumberofLeaves;
                }

                YoungStarRoot.Leaves[YoungStarNode[ChildNodeID].Leaves] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    YoungStarRoot.Leaves[YoungStarNode[ChildNodeID].Leaves+k] = 
                        sublist[YoungStarRoot.Leaves[YoungStarNode[ChildNodeID].Leaves+k-1]];
                }
                YoungStarNode[ChildNodeID].NumberofLeaves = subnumber[i];
                YoungStarNode[ChildNodeID].Level = YoungStarNode[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    YoungStarNode[ChildNodeID].Pos[k] = YoungStarNode[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*YoungStarRoot.Width*YoungStarRoot.WidthFactor[YoungStarNode[CurrentNodeID].Level];
            }
		}
        CurrentNodeID = NextNodeID;

	}
}

//////////////////////////// FRONT LINE ///////////////////////////////////////////

static double Diag;
static void YoungStarNodeDataImplant(const int CurrentNodeID){
    double Width = Diag*YoungStarRoot.WidthFactor[YoungStarNode[CurrentNodeID].Level];

    int NActives = 0;
    double DistanceMax = 0.e0;
    double LFUV = 0.e0;
    double COM[3] = {0.e0,0.e0,0.e0};
    if(YoungStarNode[CurrentNodeID].Children == NONE){
        int Number_of_leaf = YoungStarNode[CurrentNodeID].NumberofLeaves;
        int header = YoungStarNode[CurrentNodeID].Leaves;
        for(int k=0;k<Number_of_leaf;k++){
            int leaf = header+k;
            double Distance = DISTANCE(YoungStarNode[CurrentNodeID].Pos,YoungStarCache[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            LFUV += YoungStarCache[leaf].LFUV;
            COM[0] += YoungStarCache[leaf].LFUV*YoungStarCache[leaf].Pos[0];
            COM[1] += YoungStarCache[leaf].LFUV*YoungStarCache[leaf].Pos[1];
            COM[2] += YoungStarCache[leaf].LFUV*YoungStarCache[leaf].Pos[2];
            NActives += YoungStarCache[leaf].Active;
        }
        Width = 0.e0;
    } else {
        bool first = true;
        int ChildNodeID = YoungStarNode[CurrentNodeID].Children;
        do{
            YoungStarNodeDataImplant(ChildNodeID);

            DistanceMax = fmax(DistanceMax,YoungStarNode[ChildNodeID].DistanceMax);
            LFUV += YoungStarNode[ChildNodeID].LFUV;
            COM[0] += YoungStarNode[ChildNodeID].LFUV*YoungStarNode[ChildNodeID].COM[0];
            COM[1] += YoungStarNode[ChildNodeID].LFUV*YoungStarNode[ChildNodeID].COM[1];
            COM[2] += YoungStarNode[ChildNodeID].LFUV*YoungStarNode[ChildNodeID].COM[2];
            NActives += YoungStarNode[ChildNodeID].NumberofActiveLeaves;
            ChildNodeID = YoungStarNode[ChildNodeID].Sister;
        } while(ChildNodeID != NONE);
    }
    double InvLFUV = 1.e0/LFUV;
    YoungStarNode[CurrentNodeID].COM[0] = COM[0]*InvLFUV;
    YoungStarNode[CurrentNodeID].COM[1] = COM[1]*InvLFUV;
    YoungStarNode[CurrentNodeID].COM[2] = COM[2]*InvLFUV;

    YoungStarNode[CurrentNodeID].DistanceMax = DistanceMax + Width;
    YoungStarNode[CurrentNodeID].LFUV = LFUV;
    YoungStarNode[CurrentNodeID].NumberofActiveLeaves = NActives;
    return ;
}


static int NextYoungStarNode(const int NodeID){

    int CurrentNodeID = NodeID;

    if(YoungStarNode[CurrentNodeID].Sister != NONE){
        CurrentNodeID = YoungStarNode[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(YoungStarNode[YoungStarNode[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = YoungStarNode[YoungStarNode[NextNodeID].Parent].Sister;
                break;
            } else if(YoungStarNode[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = YoungStarNode[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}

void PlantYoungStarTree(void){

    if(FirstCallPlantYoungStarTree){
        InitializeRootForYoungStarTree();
        FirstCallPlantYoungStarTree = false;
    }

    int NYoungStars = YoungStarTreePreprocessing();
    //dprintlmpi(NYoungStars);
	MakeYoungStarRoot(NYoungStars);

    if(YoungStarRoot.NumberofLeaves > 0){
        BuildYoungStarTree();
        Diag = DISTANCE(YoungStarNode[0].Pos,YoungStarNode[YoungStarNode[0].Children].Pos);
        YoungStarNodeDataImplant(0);
        for(int i=1;i<YoungStarRoot.NumberofNodes;i++){
            YoungStarNode[i].Next = NextYoungStarNode(i);
        }
    }

	return;
}
