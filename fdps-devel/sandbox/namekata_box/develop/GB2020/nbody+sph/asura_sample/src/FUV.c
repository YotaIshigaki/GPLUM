#include "config.h"
#include "PlantYoungStarTree.h"
#include "StarFormation.h"

#if 0
#ifdef USE_ASRSIC //{
#include "ASRSIC.h"
# error Incorrect!??????
#undef LOCAL_COMPUTATION
#else
# error Incorrect!
#define LOCAL_COMPUTATION
#   ifdef USE_ASRSIC //{
#   undef LOCAL_COMPUTATION 
#   endif // USE_ASRSIC //}
#endif // USE_ASRSIC //}
#endif

#ifdef USE_ASRSIC //{
#include <ASRSIC.h>
#endif // USE_ASRSIC //}


static int GetNearIntegerUpper(int IntNumber){
    if (IntNumber <= 0)
        return 0;
    IntNumber -= 1;
    IntNumber |= IntNumber >> 1;
    IntNumber |= IntNumber >> 2;
    IntNumber |= IntNumber >> 4;
    IntNumber |= IntNumber >> 8;
    IntNumber |= IntNumber >> 16;
    return IntNumber + 1;
}

static int GetNearIntegerLower(int IntNumber){
    if (IntNumber <= 0)
        return 0;
    IntNumber -= 1;
    IntNumber |= IntNumber >> 1;
    IntNumber |= IntNumber >> 2;
    IntNumber |= IntNumber >> 4;
    IntNumber |= IntNumber >> 8;
    IntNumber |= IntNumber >> 16;
    return (IntNumber+1)>>1;
}


bool CheckFUVStep(void){

#ifdef USE_FUVFEEDBACK_STEP //{

#if 0
    if(CheckLastStepSF() == true){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"Wake up because of SF.\n");
            fflush(NULL);
        }

        // for(int i=0;i<Pall.Nhydro;i++){
            // Phydro[NBCache[i].Leaf]->G0thinNextUpdateTime = Pall.TCurrent;
        // }
        return true;
    }
#endif


#if 0
    // Test part
    // check
    {
        double dtmax = FUVFEEDBACK_G0_EVALUATION_INTERVAL*3;
        int FUVStep = (int)(dtmax/FUVFEEDBACK_G0_EVALUATION_INTERVAL+0.5);
        FUVStep = GetNearIntegerLower(FUVStep);
        fprintf(stderr,"Max check: %g %g, %d\n",dtmax,FUVFEEDBACK_G0_EVALUATION_INTERVAL,FUVStep);
    }

    // check
    {
        double dtmax = FUVFEEDBACK_G0_EVALUATION_INTERVAL*18;
        int FUVStep = (int)(dtmax/FUVFEEDBACK_G0_EVALUATION_INTERVAL+0.5);
        FUVStep = GetNearIntegerLower(FUVStep);
        fprintf(stderr,"Max check: %g %g, %d\n",dtmax,FUVFEEDBACK_G0_EVALUATION_INTERVAL,FUVStep);
    }

    // min check
    {
        double dtmax = FUVFEEDBACK_G0_EVALUATION_INTERVAL*0.5;
        int FUVStep = (int)(dtmax/FUVFEEDBACK_G0_EVALUATION_INTERVAL+0.5);
        FUVStep = GetNearIntegerLower(FUVStep);
        fprintf(stderr,"Min check: %g %g, %d\n",dtmax,FUVFEEDBACK_G0_EVALUATION_INTERVAL,FUVStep);
    }

    // 
    {
        int FUVStep = (int)(FUVFEEDBACK_G0_EVALUATION_INTERVAL/(5*Pall.dtmin*(Pall.UnitTime/YEAR_CGS))+0.5);
        fprintf(stderr,"dtmin: %g %g\n",FUVFEEDBACK_G0_EVALUATION_INTERVAL,Pall.dtmin*(Pall.UnitTime/YEAR_CGS));
        FUVStep = GetNearIntegerLower(FUVStep);
        fprintf(stderr,"%d %d\n",FUVStep, Pall.TStep);
    }
#endif 
    int FUVStep = (int)(fmin(Pall.dtmax*(Pall.UnitTime/YEAR_CGS),FUVFEEDBACK_G0_EVALUATION_INTERVAL)/(Pall.dtmin*(Pall.UnitTime/YEAR_CGS))+0.5);
    FUVStep = MAX(GetNearIntegerLower(FUVStep),1);
    
#if 0
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"dtmin, dt_G0: %g %g\n",Pall.dtmin*(Pall.UnitTime/YEAR_CGS), FUVFEEDBACK_G0_EVALUATION_INTERVAL);
        fprintf(stderr,"FVUStep %d Pall.TStep %d Residual %d\n",FUVStep, Pall.TStep,Pall.TStep%FUVStep);
        //fprintf(stderr,"%d %d\n",Pall.TStep%FUVStep,2%1);
        // exit(1);
    }
#endif 
    if(Pall.TStep%FUVStep == 0){ 
        return true;
    } else {
        return false;
    }

#else // USE_FUVFEEDBACK_STEP //}//{
    return true;
#endif // USE_FUVFEEDBACK_STEP //}
}


StructYoungStarRoot FUVTreeRoot;
struct StructYoungStarNode *FUVTreeNode;

static int NumberofCurrentHydroActive = 0;
static bool *CurrentHydroActive;


struct StructFUVTreeCache{
    double Pos[3];
    double LFUV;
    int  Leaf;
    bool Active;
} *FUVTreeCache;

static int NumberofAllocatedFUVTreeResultCache = 0;
struct StructFUVTreeResultCache{
    double LFUV;
    int  InteractionList; // remove
} *FUVTreeResultCache = NULL;

struct StructFUVLeavesExportImport{
    double Pos[3];
    double LFUV;
#ifdef USE_SYMMETRIZED_SOFTENING
    double Eps;
#endif // USE_SYMMETRIZED_SOFTENING
};

static double FUVTreeSoftening = 1;
static double LFUVNormalization = 1;
static double LFUVInvNormalization = 1;
static bool FirstCallFUVFEEDBACK = true;

static double FUVRmin;
static double FUVRmax;
static double FUVRmax2;

void InitFUV(void){

    FUVTreeSoftening = (10*PC_CGS/Pall.UnitLength);
    LFUVInvNormalization = 1.0/LFUVNormalization;

    NumberofAllocatedFUVTreeResultCache = FirstAllocationSize;
    FUVTreeResultCache =  malloc(sizeof(struct StructFUVTreeResultCache)*NumberofAllocatedFUVTreeResultCache);

    FUVRmin = FUVFEEDBACK_TRANCATION_RADIUS_MIN/Pall.UnitLength;
    FUVRmax = FUVFEEDBACK_TRANCATION_RADIUS_MAX/Pall.UnitLength;
    FUVRmax2 = SQ(FUVRmax);

    return ;
}

static struct StructFUVLeavesExportImport **FUVLeavesExport;
static struct StructFUVLeavesExportImport *FUVLeavesImport;
static int NumberofAllocatedLeavesImport = 0;

void InitializeRootForFUV(void){

    int NloadGravity = FirstAllocationSize;
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)FirstAllocationSize)/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = FirstAllocationSize*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));

    /* GravtityRoot */
    FUVTreeRoot.NumberofAllocatedLeaves = NloadGravity;
    FUVTreeRoot.Leaves = malloc(sizeof(int)*NloadGravity+1);
    FUVTreeRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
    FUVTreeRoot.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

    FUVTreeRoot.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForGrav;
    FUVTreeRoot.MaxLevel = TreeMaxNodeLevel;

    FUVTreeRoot.OpeningAngle = FUVFEEDBACK_THETA;
    FUVTreeRoot.NumberofLeavesInGroup = FUVFEEDBACK_NGROUP;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        FUVTreeRoot.WidthFactor[i] = pow(dw,(double)i);
    /* FUVTreeRoot */

    /* FUVode */
    FUVTreeNode = malloc(sizeof(struct StructYoungStarNode)*FUVTreeRoot.NumberofAllocatedNodes);

    struct StructYoungStarNode FUVTreeNodeTemp;
    memset(&FUVTreeNodeTemp,0,sizeof(struct StructYoungStarNode));
    FUVTreeNodeTemp.Next = FUVTreeNodeTemp.Parent = FUVTreeNodeTemp.Children =
    FUVTreeNodeTemp.Sister = NONE;

    for(int i=0;i<FUVTreeRoot.NumberofAllocatedNodes;i++)
        FUVTreeNode[i] = FUVTreeNodeTemp;
    /* FUVode */

    /* FUVTreeCache */
    FUVTreeCache = malloc(sizeof(struct StructFUVTreeCache)*NloadGravity+1);
    /* FUVTreeCache */

    /* allocate communication buffer for export */
    FUVLeavesExport = malloc(sizeof(struct StructFUVLeavesExportImport*)*MPIGetNumProcs());
    /* allocate communication buffer for export */

    return ;
}


void InitializeFUVTree(void){

#ifdef PHANTOM_FUV //{
	g5_open();
    Npipes = g5_get_number_of_pipelines();

	g5_close();
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr," ==============================\n");
        fprintf(stderr," = Use Phantom GRAPE FUV MODE =\n");
        fprintf(stderr," ==============================\n");
        fprintf(stderr," JMEMSIZE = %d\n",JMEMSIZE);
        fflush(stderr);
    }
#endif //PHANTOM_FUV //}

    return;
}


void UpdateNumberofLeavesInGroupInFUV(const int Ng){
    FUVTreeRoot.NumberofLeavesInGroup = Ng;
    return ;
}


static int FUVMakeExportParticleList(const int Index){

    if(Pall.Nstars == 0){
        return 0;
    }
    if(YoungStarRoot.NumberofLeaves == 0){
        return 0;
    }

#if 0
    CheckSizeofBufferExportSendIndex((YoungStarNode[0].NumberofLeaves+1),
            sizeof(struct StructFUVLeavesExportImport),Index);
    FUVLeavesExport[Index] = BufferExportSend[Index];

    int S = 0;
    int header = YoungStarNode[0].Leaves;
    int NumberofLeaves = YoungStarNode[0].NumberofLeaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
        FUVLeavesExport[Index][S].Pos[0] = YoungStarCache[leaf].Pos[0];
        FUVLeavesExport[Index][S].Pos[1] = YoungStarCache[leaf].Pos[1];
        FUVLeavesExport[Index][S].Pos[2] = YoungStarCache[leaf].Pos[2];
        FUVLeavesExport[Index][S].LFUV   = YoungStarCache[leaf].LFUV;
        S ++;
    }


    return S;
#endif

    int ExportNodeID = CommunicationTable[Index].SendRank;

    double PosMax[3] = {EdgesForActiveHydro[ExportNodeID].PosMax[0],
                        EdgesForActiveHydro[ExportNodeID].PosMax[1],
                        EdgesForActiveHydro[ExportNodeID].PosMax[2]};
    double PosMin[3] = {EdgesForActiveHydro[ExportNodeID].PosMin[0],
                        EdgesForActiveHydro[ExportNodeID].PosMin[1],
                        EdgesForActiveHydro[ExportNodeID].PosMin[2]};

    double Pos[3] = {0.5*(PosMax[0]+PosMin[0]),
                     0.5*(PosMax[1]+PosMin[1]),
                     0.5*(PosMax[2]+PosMin[2])};
    double hWidth[3] = {PosMax[0]-Pos[0],PosMax[1]-Pos[1],PosMax[2]-Pos[2]};
    if(hWidth[0] == 0.e0){ // No target!
        return 0;
    }

#ifdef USE_FUVFEEDBACK_TRANCATION_RADIUS //{
#if 0
    // Distance between the local domain and external domain.
    double NodeToNodeDist[3] = {
    };

    if()
#endif
#endif // USE_FUVFEEDBACK_TRANCATION_RADIUS //}

    int SendThisTime = 0;
    double theta2 = SQ(FUVTreeRoot.OpeningAngle);
    int RootNodeID = 0;
    int CurrentNodeID = YoungStarNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
        double x[3] = {fabs(YoungStarNode[CurrentNodeID].Pos[0]-Pos[0]),
                       fabs(YoungStarNode[CurrentNodeID].Pos[1]-Pos[1]),
                       fabs(YoungStarNode[CurrentNodeID].Pos[2]-Pos[2])};
        double dx2 = 0.e0;
        double hwidth = 0.5*YoungStarRoot.Width*YoungStarRoot.WidthFactor[YoungStarNode[CurrentNodeID].Level];
        // if(x[0] > (hWidth[0]+hwidth)) dx2 += SQ(x[0]-(hWidth[0]+hwidth));
        // if(x[1] > (hWidth[1]+hwidth)) dx2 += SQ(x[1]-(hWidth[1]+hwidth));
        // if(x[2] > (hWidth[2]+hwidth)) dx2 += SQ(x[2]-(hWidth[2]+hwidth));

        dx2 += SQ(fmax(x[0]-(hWidth[0]+hwidth),0.0));
        dx2 += SQ(fmax(x[1]-(hWidth[1]+hwidth),0.0));
        dx2 += SQ(fmax(x[2]-(hWidth[2]+hwidth),0.0));


        if(SQ(2.0*hwidth) < theta2*dx2){
            CheckSizeofBufferExportSendIndex((SendThisTime+1),
                    sizeof(struct StructFUVLeavesExportImport),Index);
            FUVLeavesExport[Index] = BufferExportSend[Index];

            FUVLeavesExport[Index][SendThisTime].Pos[0] = YoungStarNode[CurrentNodeID].COM[0];
            FUVLeavesExport[Index][SendThisTime].Pos[1] = YoungStarNode[CurrentNodeID].COM[1];
            FUVLeavesExport[Index][SendThisTime].Pos[2] = YoungStarNode[CurrentNodeID].COM[2];
            FUVLeavesExport[Index][SendThisTime].LFUV   = YoungStarNode[CurrentNodeID].LFUV;

            SendThisTime ++;
            CurrentNodeID = YoungStarNode[CurrentNodeID].Next;
		} else if (YoungStarNode[CurrentNodeID].Children == NONE) {
            int NumberofLeaves = YoungStarNode[CurrentNodeID].NumberofLeaves;
            CheckSizeofBufferExportSendIndex((SendThisTime+NumberofLeaves),
                    sizeof(struct StructFUVLeavesExportImport),Index);
            FUVLeavesExport[Index] = BufferExportSend[Index];

            int header = YoungStarNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                FUVLeavesExport[Index][SendThisTime].Pos[0] = YoungStarCache[leaf].Pos[0];
                FUVLeavesExport[Index][SendThisTime].Pos[1] = YoungStarCache[leaf].Pos[1];
                FUVLeavesExport[Index][SendThisTime].Pos[2] = YoungStarCache[leaf].Pos[2];
                FUVLeavesExport[Index][SendThisTime].LFUV   = YoungStarCache[leaf].LFUV;
                SendThisTime ++;
            }
			CurrentNodeID = YoungStarNode[CurrentNodeID].Next;
		} else {
			CurrentNodeID = YoungStarNode[CurrentNodeID].Children;
		}
	}

    return SendThisTime;
}


static int NumberofFUVLeavesAllocated = 0;
static int *FUVTreeSublist; 

static void FUVTreePreprocessing(const int NumberofLeavesThisStep){

    if(NumberofLeavesThisStep > FUVTreeRoot.NumberofAllocatedLeaves){
        FUVTreeRoot.NumberofAllocatedLeaves = (int)(ForAngelsShare*NumberofLeavesThisStep);
        free(FUVTreeRoot.Leaves);
        free(FUVTreeCache);
        FUVTreeRoot.Leaves = malloc(sizeof(int)*FUVTreeRoot.NumberofAllocatedLeaves);
        FUVTreeCache = malloc(sizeof(struct StructFUVTreeCache)*FUVTreeRoot.NumberofAllocatedLeaves);
    }

    if(NumberofLeavesThisStep > NumberofFUVLeavesAllocated){
        if(NumberofFUVLeavesAllocated > 0)
            free(FUVTreeSublist);
        NumberofFUVLeavesAllocated = (int)(ForAngelsShare*NumberofLeavesThisStep);
	    FUVTreeSublist = malloc(sizeof(int)*NumberofFUVLeavesAllocated);
    }

    return ;
}


static void MakeFUVTreeRoot(const int NumberofLeavesThisStep){

    int RootNodeID = 0;

	for(int k=0;k<3;k++){
        FUVTreeNode[RootNodeID].Pos[k] = YoungStarNode[RootNodeID].Pos[k];
        FUVTreeRoot.PosMax[k] = YoungStarRoot.PosMax[k];
        FUVTreeRoot.PosMin[k] = YoungStarRoot.PosMin[k];
    }
    FUVTreeRoot.Width = YoungStarRoot.Width;

    FUVTreeNode[RootNodeID].Next = NONE;
    FUVTreeNode[RootNodeID].Parent = NONE;
    FUVTreeNode[RootNodeID].Sister = NONE;
    FUVTreeNode[RootNodeID].Children = NONE;

    FUVTreeNode[RootNodeID].NumberofLeaves = NumberofLeavesThisStep;

    FUVTreeNode[RootNodeID].Level = 0;
    FUVTreeNode[RootNodeID].Leaves = 0;

    for(int i=0;i<NumberofLeavesThisStep;i++)
        FUVTreeRoot.Leaves[i] = i;

    FUVTreeRoot.NumberofLeaves = NumberofLeavesThisStep;

	return ;
}

static inline bool __attribute__((always_inline)) FUVTreeNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber){

	if( (FUVTreeNode[CurrentNodeID].NumberofLeaves <= CriticalNumber) || FUVTreeNode[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}

static void BuildFUVTree(void){


    int NumberofNodes = 0; 

    int NumberofNodeCreationLimit = FUVTreeRoot.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
	while(1){

		if(FUVTreeNodeSeparationCriterion(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(FUVTreeNode[CurrentNodeID].Sister != NONE){
				CurrentNodeID = FUVTreeNode[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(FUVTreeNode[FUVTreeNode[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = FUVTreeNode[FUVTreeNode[NextNodeID].Parent].Sister;
						break;
					} else if(FUVTreeNode[NextNodeID].Parent == RootNodeID){
                        FUVTreeRoot.CurrentMaxLevel = CurrentMaxLevel;
                        FUVTreeRoot.NumberofNodes = NumberofNodes + 1;
						return;
					}
                    NextNodeID = FUVTreeNode[NextNodeID].Parent;
				}
			}
			continue;
		}

        int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 
		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

        int NumberofLeaves = FUVTreeNode[CurrentNodeID].NumberofLeaves;
        int header = FUVTreeNode[CurrentNodeID].Leaves;
		for(int i=0;i<NumberofLeaves;i++){
            int leaf = FUVTreeRoot.Leaves[header+i];

            double Pos[3] = {FUVLeavesImport[leaf].Pos[0],FUVLeavesImport[leaf].Pos[1],FUVLeavesImport[leaf].Pos[2]};

            int subindex0 = ((FUVTreeNode[CurrentNodeID].Pos[0] <= Pos[0])?1:0);
            int subindex1 = ((FUVTreeNode[CurrentNodeID].Pos[1] <= Pos[1])?1:0);
            int subindex2 = ((FUVTreeNode[CurrentNodeID].Pos[2] <= Pos[2])?1:0);
            int subindex = subindex0 | subindex1 << 1 | subindex2 << 2;

			if(subnumber[subindex] > 0){
				FUVTreeSublist[subcurrent[subindex]] = leaf;
			} else {
				subhead[subindex] = leaf;
			}
            subcurrent[subindex] = leaf;
			subnumber[subindex] ++;
        }


        ChildNodeID = CurrentNodeID;
		for(int i=0;i<TreeNsub;i++){
			if(subnumber[i] != 0){
				BackwardNodeID = ChildNodeID; 
                // make node
                NumberofNodes ++;
                ChildNodeID = NumberofNodes;
                if(NumberofNodes >= FUVTreeRoot.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(MAX(ForAngelsShare*NumberofNodes,NAdditionUnit));
                    FUVTreeNode = realloc(FUVTreeNode,sizeof(struct StructYoungStarNode)*NumberofAllocatedNodes);
                    FUVTreeRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                FUVTreeNode[ChildNodeID].Next = NONE;
                FUVTreeNode[ChildNodeID].Parent = NONE;
                FUVTreeNode[ChildNodeID].Sister = NONE;
                FUVTreeNode[ChildNodeID].Children = NONE;

                FUVTreeNode[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    FUVTreeNode[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    FUVTreeNode[ChildNodeID].Leaves = FUVTreeNode[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,FUVTreeNode[CurrentNodeID].Level+1);
                } else {
                    FUVTreeNode[BackwardNodeID].Sister = ChildNodeID;
                    FUVTreeNode[ChildNodeID].Leaves = 
                        FUVTreeNode[BackwardNodeID].Leaves + FUVTreeNode[BackwardNodeID].NumberofLeaves;
                }

                int cheader = FUVTreeNode[ChildNodeID].Leaves;
                FUVTreeRoot.Leaves[cheader] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    FUVTreeRoot.Leaves[cheader+k] = 
                        FUVTreeSublist[FUVTreeRoot.Leaves[cheader+k-1]];
                }
                FUVTreeNode[ChildNodeID].NumberofLeaves = subnumber[i];
                FUVTreeNode[ChildNodeID].Level = FUVTreeNode[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    FUVTreeNode[ChildNodeID].Pos[k] = FUVTreeNode[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*FUVTreeRoot.Width*FUVTreeRoot.WidthFactor[FUVTreeNode[CurrentNodeID].Level];
            }
		}
        CurrentNodeID = NextNodeID;
	}
}


static void CopyFromImportLeavesToFUVTreeCache(void){

    int NumberofLeaves = FUVTreeNode[0].NumberofLeaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf =  FUVTreeRoot.Leaves[k];
        FUVTreeCache[k].Pos[0] = FUVLeavesImport[leaf].Pos[0];
        FUVTreeCache[k].Pos[1] = FUVLeavesImport[leaf].Pos[1];
        FUVTreeCache[k].Pos[2] = FUVLeavesImport[leaf].Pos[2];
        FUVTreeCache[k].LFUV   = FUVLeavesImport[leaf].LFUV;
    }
    return;
}


/////////////////////////////// LET FOR FUV /////////////////////////////// 
static int NextFUVTreeNode(const int NodeID){

    int CurrentNodeID = NodeID;

    if(FUVTreeNode[CurrentNodeID].Sister != NONE){
        CurrentNodeID = FUVTreeNode[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(FUVTreeNode[FUVTreeNode[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = FUVTreeNode[FUVTreeNode[NextNodeID].Parent].Sister;
                break;
            } else if(FUVTreeNode[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = FUVTreeNode[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}

static double Diag;
static void FUVTreeNodeDataImplant(const int CurrentNodeID){
    double Width = Diag*FUVTreeRoot.WidthFactor[FUVTreeNode[CurrentNodeID].Level];

    int NActives = 0;
    double DistanceMax = 0.e0;
    double LFUV = 0.e0;
    double COM[3] = {0.e0,0.e0,0.e0};
    double Eps2 = 0.e0;
    if(FUVTreeNode[CurrentNodeID].Children == NONE){
        int Number_of_leaf = FUVTreeNode[CurrentNodeID].NumberofLeaves;
        int header = FUVTreeNode[CurrentNodeID].Leaves;
        for(int k=0;k<Number_of_leaf;k++){
            int leaf = header+k;
            double Distance = DISTANCE(FUVTreeNode[CurrentNodeID].Pos,FUVTreeCache[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            LFUV += FUVTreeCache[leaf].LFUV;
            COM[0] += FUVTreeCache[leaf].LFUV*FUVTreeCache[leaf].Pos[0];
            COM[1] += FUVTreeCache[leaf].LFUV*FUVTreeCache[leaf].Pos[1];
            COM[2] += FUVTreeCache[leaf].LFUV*FUVTreeCache[leaf].Pos[2];
            NActives += FUVTreeCache[leaf].Active;
        }
        Width = 0.e0;
    } else {
        bool first = true;
        int ChildNodeID = FUVTreeNode[CurrentNodeID].Children;
        while(ChildNodeID != NONE){
            FUVTreeNodeDataImplant(ChildNodeID);

            DistanceMax = fmax(DistanceMax,FUVTreeNode[ChildNodeID].DistanceMax);
            LFUV += FUVTreeNode[ChildNodeID].LFUV;
            COM[0] += FUVTreeNode[ChildNodeID].LFUV*FUVTreeNode[ChildNodeID].COM[0];
            COM[1] += FUVTreeNode[ChildNodeID].LFUV*FUVTreeNode[ChildNodeID].COM[1];
            COM[2] += FUVTreeNode[ChildNodeID].LFUV*FUVTreeNode[ChildNodeID].COM[2];
            NActives += FUVTreeNode[ChildNodeID].NumberofActiveLeaves;
            ChildNodeID = FUVTreeNode[ChildNodeID].Sister;
        }
    }
    double InvLFUV = 1.e0/LFUV;
    FUVTreeNode[CurrentNodeID].COM[0] = COM[0]*InvLFUV;
    FUVTreeNode[CurrentNodeID].COM[1] = COM[1]*InvLFUV;
    FUVTreeNode[CurrentNodeID].COM[2] = COM[2]*InvLFUV;

    FUVTreeNode[CurrentNodeID].DistanceMax = DistanceMax + Width;
    FUVTreeNode[CurrentNodeID].LFUV = LFUV;
    FUVTreeNode[CurrentNodeID].NumberofActiveLeaves = NActives;
    return ;

}

static void PlantFUVTree(const int NumberofLeavesThisStep){

    FUVTreePreprocessing(NumberofLeavesThisStep);

	MakeFUVTreeRoot(NumberofLeavesThisStep);
    BuildFUVTree();
    CopyFromImportLeavesToFUVTreeCache();

    for(int i=1;i<FUVTreeRoot.NumberofNodes;i++){
        FUVTreeNode[i].Next = NextFUVTreeNode(i);
    }
    double diag = DISTANCE(FUVTreeNode[0].Pos,FUVTreeNode[FUVTreeNode[0].Children].Pos);
    FUVTreeNodeDataImplant(0);

	return;
}


/////////////////////////////// LET FOR FUV /////////////////////////////// 
static int NumberofFieldAllocatedFUVFEEDBACK = 0;

static double (*FieldFUVPos)[3]; 
static double *FieldLFUV;

static int GetFieldFromFUVTree(const int CurrentNodeID, const int NField){

    int NumberofField = NField;
    double cwidth = 0.5*GravityRoot.Width*GravityRoot.WidthFactor[HydroNode[CurrentNodeID].Level];
    double theta2 = SQ(FUVTreeRoot.OpeningAngle);

    int RootNodeID = 0;
	int TargetNodeID = FUVTreeNode[RootNodeID].Children;
	while(TargetNodeID != RootNodeID){ 
        double twidth = 0.5*FUVTreeRoot.Width*FUVTreeRoot.WidthFactor[FUVTreeNode[TargetNodeID].Level];
        double sqDist = 0.e0;

        double Dist[3] = {fabs(HydroNode[CurrentNodeID].Pos[0]-FUVTreeNode[TargetNodeID].Pos[0]),
                          fabs(HydroNode[CurrentNodeID].Pos[1]-FUVTreeNode[TargetNodeID].Pos[1]),
                          fabs(HydroNode[CurrentNodeID].Pos[2]-FUVTreeNode[TargetNodeID].Pos[2])};
        double hwidth = cwidth + twidth;
        sqDist += SQ(fmax(Dist[0]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[1]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[2]-hwidth,0.0));

#ifdef USE_FUVFEEDBACK_TRANCATION_RADIUS //{
        if(sqDist >= FUVRmax2){
			TargetNodeID = FUVTreeNode[TargetNodeID].Next;
        } else 
#endif // USE_FUVFEEDBACK_TRANCATION_RADIUS //{
        if(SQ(2.0*twidth)<theta2*sqDist){
            if(NumberofField+1 >= NumberofFieldAllocatedFUVFEEDBACK){
                NumberofFieldAllocatedFUVFEEDBACK = (int)(ForAngelsShare*(NumberofField+1));
                FieldFUVPos = realloc(FieldFUVPos,sizeof(double)*3*NumberofFieldAllocatedFUVFEEDBACK);
                FieldLFUV = realloc(FieldLFUV,sizeof(double)*NumberofFieldAllocatedFUVFEEDBACK);
            }

            FieldFUVPos[NumberofField][0] = FUVTreeNode[TargetNodeID].COM[0]-HydroNode[CurrentNodeID].Pos[0];
            FieldFUVPos[NumberofField][1] = FUVTreeNode[TargetNodeID].COM[1]-HydroNode[CurrentNodeID].Pos[1];
            FieldFUVPos[NumberofField][2] = FUVTreeNode[TargetNodeID].COM[2]-HydroNode[CurrentNodeID].Pos[2];
            FieldLFUV[NumberofField]      = FUVTreeNode[TargetNodeID].LFUV;
    
            NumberofField ++;
			TargetNodeID = FUVTreeNode[TargetNodeID].Next;
		 } else if (FUVTreeNode[TargetNodeID].Children == NONE){

            int NumberofLeaves = FUVTreeNode[TargetNodeID].NumberofLeaves;
            if(NumberofField+NumberofLeaves+1 > NumberofFieldAllocatedFUVFEEDBACK){
                NumberofFieldAllocatedFUVFEEDBACK = (int)(ForAngelsShare*(NumberofField+NumberofLeaves+1));
                FieldFUVPos = realloc(FieldFUVPos,sizeof(double)*3*NumberofFieldAllocatedFUVFEEDBACK);
                FieldLFUV = realloc(FieldLFUV,sizeof(double)*NumberofFieldAllocatedFUVFEEDBACK);
            }
            int header = FUVTreeNode[TargetNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header + k;
                FieldFUVPos[NumberofField][0] = FUVTreeCache[leaf].Pos[0]-HydroNode[CurrentNodeID].Pos[0];
                FieldFUVPos[NumberofField][1] = FUVTreeCache[leaf].Pos[1]-HydroNode[CurrentNodeID].Pos[1];
                FieldFUVPos[NumberofField][2] = FUVTreeCache[leaf].Pos[2]-HydroNode[CurrentNodeID].Pos[2];
                FieldLFUV[NumberofField]      = FUVTreeCache[leaf].LFUV;

                NumberofField ++;
            }
			TargetNodeID = FUVTreeNode[TargetNodeID].Next;
		 } else {
			TargetNodeID = FUVTreeNode[TargetNodeID].Children;
		 }
	}

	return (NumberofField);
}

#ifdef PHANTOM_FUV //{
static void CalculateLFUVEngine(const int CurrentNodeID){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    int LeavesLog[Npipes];

    int NumberofActiveLeaves = HydroNode[CurrentNodeID].NumberofActiveLeaves;
    int header = HydroNode[CurrentNodeID].Leaves;
    int CurrentLeafID = 0;
    for(int i=0;i<NumberofActiveLeaves;i+=Npipes){
		int NActives = MIN(Npipes,NumberofActiveLeaves-i);
        int Count = 0;
        while(Count!=NActives){
            int leaf = header + CurrentLeafID;
            if(NBCache[leaf].Active){
#ifdef USE_SHIFT_GRAVITYFIELD
                Dxi[Count][0] = NBCache[leaf].Pos[0]-HydroNode[CurrentNodeID].Pos[0];
                Dxi[Count][1] = NBCache[leaf].Pos[1]-HydroNode[CurrentNodeID].Pos[1];
                Dxi[Count][2] = NBCache[leaf].Pos[2]-HydroNode[CurrentNodeID].Pos[2];
#else
                Dxi[Count][0] = NBCache[leaf].Pos[0];
                Dxi[Count][1] = NBCache[leaf].Pos[1];
                Dxi[Count][2] = NBCache[leaf].Pos[2];
#endif

                Depsi[Count] = SQ(FUVTreeSoftening);
                LeavesLog[Count] = NBCache[leaf].Leaf;
                Count ++;
            }
            CurrentLeafID ++;
        }
		if(NActives<Npipes){ // fill empty pipelines with dummy data.
			for(int k=NActives;k<Npipes;k++){
				Dxi[k][0] = Dxi[NActives-1][0];
				Dxi[k][1] = Dxi[NActives-1][1];
				Dxi[k][2] = Dxi[NActives-1][2];
				Depsi[k] = Depsi[NActives-1];
			} 
		}
#ifdef PHANTOM_FUV //{

#if 0
#ifdef USE_SYMMETRIZED_SOFTENING
#ifdef HAVE_AVX_PHANTOM_GRAPE //{
        g5_calculate_force_on_x0(Dxi,adummy,phidummy,Npipes,Depsi);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
        g5_set_xi(Npipes,Dxi);
        g5_set_eps(Npipes,Depsi);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else 
        g5_set_eps_to_all(Depsi[0]);
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#endif
#endif

#else // PHANTOM_FUV //}//{


#endif // PHANTOM_FUV //}
        for(int k=0;k<NActives;k++){
			FUVTreeResultCache[LeavesLog[k]].Pot += phidummy[k];
        }
    }

	return;
}
#endif // PHANTOM_FUV //{

#ifdef ASRSIC //{
static void CalculateLFUVEngine(const int CurrentNodeID){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    int LeavesLog[Npipes];

    int NumberofActiveLeaves = HydroNode[CurrentNodeID].NumberofActiveLeaves;
    int header = HydroNode[CurrentNodeID].Leaves;
    int CurrentLeafID = 0;
    for(int i=0;i<NumberofActiveLeaves;i+=Npipes){
		int NActives = MIN(Npipes,NumberofActiveLeaves-i);
        int Count = 0;
        while(Count!=NActives){
            int leaf = header + CurrentLeafID;
            if(NBCache[leaf].Active){
#ifdef USE_SHIFT_GRAVITYFIELD
                Dxi[Count][0] = NBCache[leaf].Pos[0]-HydroNode[CurrentNodeID].Pos[0];
                Dxi[Count][1] = NBCache[leaf].Pos[1]-HydroNode[CurrentNodeID].Pos[1];
                Dxi[Count][2] = NBCache[leaf].Pos[2]-HydroNode[CurrentNodeID].Pos[2];
#else
                Dxi[Count][0] = NBCache[leaf].Pos[0];
                Dxi[Count][1] = NBCache[leaf].Pos[1];
                Dxi[Count][2] = NBCache[leaf].Pos[2];
#endif

                Depsi[Count] = SQ(FUVTreeSoftening);
                LeavesLog[Count] = NBCache[leaf].Leaf;
                Count ++;
            }
            CurrentLeafID ++;
        }
		if(NActives<Npipes){ // fill empty pipelines with dummy data.
			for(int k=NActives;k<Npipes;k++){
				Dxi[k][0] = Dxi[NActives-1][0];
				Dxi[k][1] = Dxi[NActives-1][1];
				Dxi[k][2] = Dxi[NActives-1][2];
				Depsi[k] = Depsi[NActives-1];
			} 
		}
#ifdef PHANTOM_FUV //{

#else // PHANTOM_FUV //}//{


#endif // PHANTOM_FUV //}
        for(int k=0;k<NActives;k++){
			FUVTreeResultCache[LeavesLog[k]].Pot += phidummy[k];
        }
    }

	return;
}
#endif // ASRSIC //{

static double FUVChangeOverFunction(const double r){
    double r_current = fmax(FUVRmin,fmin(FUVRmax,r));
    double t = (FUVRmax-r_current)/(FUVRmax-FUVRmin);
    return -2*CUBE(t)+3*SQ(t);
}

static void AccumulateLFUV(const int CurrentHydroNodeID, const int NumberofField, double Pos[restrict][3], double LFUV[restrict]){

#ifdef USE_ASRSIC //{
    int NumberofLeaves = HydroNode[CurrentHydroNodeID].NumberofLeaves;

    double xi[NumberofLeaves][3];
    double epsi[NumberofLeaves];
    double Flux[NumberofLeaves];

    int header = HydroNode[CurrentHydroNodeID].Leaves;
    int counter = 0;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header+i;
        if(NBCache[leaf].Active){
            xi[counter][0] = NBCache[leaf].Pos[0]-HydroNode[CurrentHydroNodeID].Pos[0];
            xi[counter][1] = NBCache[leaf].Pos[1]-HydroNode[CurrentHydroNodeID].Pos[1];
            xi[counter][2] = NBCache[leaf].Pos[2]-HydroNode[CurrentHydroNodeID].Pos[2];
            epsi[counter] = 0.e0;
            Flux[counter] = 0.e0;
            counter ++;
        }
    }

    ASRSIC_CalcFlux(counter,xi,epsi,Flux,NumberofField,Pos,LFUV);

    counter = 0;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header+i;
        if(NBCache[leaf].Active){
            FUVTreeResultCache[leaf].LFUV += Flux[counter];
            counter ++;
        }
    }
#else // USE_ASRSIC //}//{
    int NumberofLeaves = HydroNode[CurrentHydroNodeID].NumberofLeaves;

    int header = HydroNode[CurrentHydroNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header+i;
        if(NBCache[leaf].Active){
            double xi[3];
            xi[0] = NBCache[leaf].Pos[0]-HydroNode[CurrentHydroNodeID].Pos[0];
            xi[1] = NBCache[leaf].Pos[1]-HydroNode[CurrentHydroNodeID].Pos[1];
            xi[2] = NBCache[leaf].Pos[2]-HydroNode[CurrentHydroNodeID].Pos[2];

            FUVTreeResultCache[leaf].LFUV = 0.e0;
            for(int k=0;k<NumberofField;k++){
                double r2 = DISTANCE2(xi,FieldFUVPos[k]);
                FUVTreeResultCache[leaf].LFUV += FieldLFUV[k]/r2;
            }
        }
    }
#endif // USE_ASRSIC //}

#ifdef PHANTOM_FUV
#error
#endif

    return ;
}

static void WalkLocalTreeAndGetLFUV(const int NImportAll){

    // OK
    if(NumberofFieldAllocatedFUVFEEDBACK == 0){
        NumberofFieldAllocatedFUVFEEDBACK = FirstAllocationSize;
        FieldFUVPos = malloc(sizeof(double)*3*NumberofFieldAllocatedFUVFEEDBACK);
        FieldLFUV = malloc(sizeof(double)*NumberofFieldAllocatedFUVFEEDBACK);
    }

    int RootNodeID = 0;

    /// Walk hydro tree and then check GetFieldFUVTree
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
		if(HydroNode[CurrentNodeID].NumberofActiveLeaves == 0){ 
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}else if( (HydroNode[CurrentNodeID].NumberofLeaves < FUVTreeRoot.NumberofLeavesInGroup)
                || (HydroNode[CurrentNodeID].Children == NONE)){
            int NumberofField = GetFieldFromFUVTree(CurrentNodeID,0);

            AccumulateLFUV(CurrentNodeID,NumberofField,FieldFUVPos,FieldLFUV);

            CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else {
		    CurrentNodeID = HydroNode[CurrentNodeID].Children;
		}
	}
    return;
}


#ifdef FUVFEEDBACK_G0_EVALUATION_INTERVAL //{
static int CheckTargetNumber(void){

    if(NumberofCurrentHydroActive < Pall.Nhydro){
        NumberofCurrentHydroActive = MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit);
        CurrentHydroActive = realloc(CurrentHydroActive,sizeof(bool)*NumberofCurrentHydroActive);
    }


    int NTarget = 0;
    double TimeConvertFactor = YEAR_CGS/Pall.UnitTime;
    for(int i=0;i<Pall.Nhydro;i++){
        CurrentHydroActive[i] = false;
        //if((NBCache[i].Active == true)&&(Pall.TCurrent>=Phydro[NBCache[i].Leaf]->G0thinNextUpdateTime)){
        if(NBCache[i].Active){
            CurrentHydroActive[i] = true;
            Phydro[NBCache[i].Leaf]->G0thinNextUpdateTime = Pall.TCurrent+
                FUVFEEDBACK_G0_EVALUATION_INTERVAL*TimeConvertFactor;
            NTarget ++;
        }
    }

    return NTarget;
}
#endif // FUVFEEDBACK_G0_EVALUATION_INTERVAL //}

#ifdef USE_FUVFEEDBACK_LOCAL_CORRECTION //{
static int GetNeighborsFUV(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = FUVTreeNode[RootNodeID].Children;

	while(CurrentNodeID != RootNodeID){

#ifdef PERIODIC_RUN  //{
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],FUVTreeNode[CurrentNodeID].Pos[k],k));
        }
#else // PERIODIC_RUN //}//{
        double dx2 = DISTANCE2(FUVTreeNode[CurrentNodeID].Pos,Pos);
#endif // PERIODIC_RUN //}

#if 0
#ifdef USE_FUVFEEDBACK_TRANCATION_RADIUS //{
        if(dx2 >= FUVRmax2){
			CurrentNodeID = FUVTreeNode[CurrentNodeID].Next;
        } else 
#endif // USE_FUVFEEDBACK_TRANCATION_RADIUS //{
#endif
        if( dx2 > SQ(h+FUVTreeNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = FUVTreeNode[CurrentNodeID].Next;
		} else if(FUVTreeNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = FUVTreeNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = FUVTreeNode[CurrentNodeID].NumberofLeaves;
            int header = FUVTreeNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
#ifdef PERIODIC_RUN  //{
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],FUVTreeCache[leaf].Pos[l],l));
#else // PERIODIC_RUN //}//{
                double distance2 = DISTANCE2(FUVTreeCache[leaf].Pos,Pos);
#endif // PERIODIC_RUN //}
                if(distance2 < hh){
                    list[nlist] = leaf;
                    nlist ++;
                }
			}
			CurrentNodeID = FUVTreeNode[CurrentNodeID].Next;
		}
	}
	return nlist;

}
#endif // USE_FUVFEEDBACK_LOCAL_CORRECTION //}


static inline double __attribute__((always_inline)) EvaluateColumnDensity(const int index){

    double X = 1.e0 - HeliumAbandance - Phydro[index]->Z;
    double nH = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho*X;

    //return pow(nH,2.0/3.0)*1.4e20;
    //return pow(nH,2.0/3.0)*3e19;
#ifdef USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //{
    return nH*FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE;
#else // USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}//{
    // 2 is a calibration factor.
    return nH*(Phydro[index]->Rho/fabs(2*Phydro[index]->GradRho))*Pall.UnitLength;
#endif // USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}

}

void CalcFUV(void){

#ifndef USE_FUVFEEDBACK //{
    return ;
#else

#ifdef USE_FUVFEEDBACK_CONSTANT_VALUE //{
    return ;
#endif // USE_FUVFEEDBACK_CONSTANT_VALUE //}

#ifdef USE_FUVFEEDBACK_STEP //{ 
    if(CheckFUVStep() == false){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"Skip FUV \n");
        }
        return ;
    }
#endif //USE_FUVFEEDBACK_STEP //}


    if(FirstCallFUVFEEDBACK == true){
        InitFUV();
        InitializeRootForFUV();

        FirstCallFUVFEEDBACK = false;
    }

#ifdef FUVFEEDBACK_G0_EVALUATION_INTERVAL //{
    // If there is no target hydro particle, this routine is skipped.
    int NumberofTargetHydro = CheckTargetNumber();
    MPI_Allreduce(MPI_IN_PLACE,&NumberofTargetHydro,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(NumberofTargetHydro == 0){
        return ;
    }
#endif // FUVFEEDBACK_G0_EVALUATION_INTERVAL //}


    const int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    const double MassConversionFactor = Pall.UnitMass/MSUN_CGS;
    const double TimeConversionFactor = Pall.UnitTime/YEAR_CGS;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            Pstar[i]->LFUV = MassConversionFactor*Pstar[i]->InitialMass*
                ASRFLXGetFUV((Pall.TCurrent-Pstar[i]->FormationTime)*TimeConversionFactor,Pstar[i]->Z); 
        }
    }

    PlantYoungStarTree();
    int GlobalNYongStar;
    MPI_Allreduce(&YoungStarRoot.NumberofLeaves,&GlobalNYongStar,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(GlobalNYongStar== 0){
        return;
    }

    if(NumberofAllocatedFUVTreeResultCache < Pall.Nhydro){
        NumberofAllocatedFUVTreeResultCache = ForAngelsShare*Pall.Nhydro;
        FUVTreeResultCache = realloc(FUVTreeResultCache,
                sizeof(struct StructFUVTreeResultCache)*NumberofAllocatedFUVTreeResultCache);
    }
    for(int i=0;i<Pall.Nhydro;i++){
        FUVTreeResultCache[i].LFUV = 0.e0;
        FUVTreeResultCache[i].InteractionList = 0.0; // Remove
    }

#if 1

    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    int SendFlag,RecvFlag;

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(NImportThisTime+i,1,MPI_INT,
            CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_EXPORT+i,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
    }
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = FUVMakeExportParticleList(i);
        MPI_Isend(NExportThisTime+i,1,MPI_INT,
                CommunicationTable[i].SendRank,TAG_SPH_DENSITY_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);
    }
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportAll += NImportThisTime[i];
    }

    if(NumberofAllocatedLeavesImport < NImportAll){
        NumberofAllocatedLeavesImport = ForAngelsShare*NImportAll;
        FUVLeavesImport = realloc(FUVLeavesImport,NumberofAllocatedLeavesImport*sizeof(struct StructFUVLeavesExportImport));
    }

#else
    //////////////////////////////////
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = FUVMakeExportParticleList(i);
    }

    int NImportAll = 0;
    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MPIGetMyID()] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImportAll += NImportThisTime[i];
    }
    if(NumberofAllocatedLeavesImport < NImportAll){
        NumberofAllocatedLeavesImport = ForAngelsShare*NImportAll;
        FUVLeavesImport = realloc(FUVLeavesImport,NumberofAllocatedLeavesImport*sizeof(struct StructFUVLeavesExportImport));
    }

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
#endif

    int NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    //int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i] > 0){
            MPI_Isend(FUVLeavesExport[i],
                NExportThisTime[i]*sizeof(struct StructFUVLeavesExportImport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_FORCE_PARALLELTREEGRAPE_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i] > 0){
            MPI_Irecv(FUVLeavesImport+NImport,
                NImportThisTime[i]*sizeof(struct StructFUVLeavesExportImport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FORCE_PARALLELTREEGRAPE_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }
    double Tcomm = GetElapsedTime();
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.GravityCommThisStep = GetElapsedTime()-Tcomm;
    //////////////////////////////////


    NImport += YoungStarRoot.NumberofLeaves;
    if(NumberofAllocatedLeavesImport < NImport){
        NumberofAllocatedLeavesImport = ForAngelsShare*NImport;
        FUVLeavesImport = realloc(FUVLeavesImport,NumberofAllocatedLeavesImport*sizeof(struct StructFUVLeavesExportImport));
    }

    int IndexShift = YoungStarRoot.NumberofLeaves;
    for(int i=0;i<NImportAll;i++){ // move
        int BackwardID  = NImportAll-1-i;
        FUVLeavesImport[BackwardID+IndexShift] = FUVLeavesImport[BackwardID];
    }

    // YoungStarCache is prepared in PlantYoungStarTree
    for(int i=0;i<YoungStarRoot.NumberofLeaves;i++){ // copy
        FUVLeavesImport[i].Pos[0] = YoungStarCache[i].Pos[0];
        FUVLeavesImport[i].Pos[1] = YoungStarCache[i].Pos[1];
        FUVLeavesImport[i].Pos[2] = YoungStarCache[i].Pos[2];
        FUVLeavesImport[i].LFUV   = YoungStarCache[i].LFUV;
    }
    NImportAll = NImport;//+YoungStarRoot.NumberofLeaves;

    if(NImportAll>0){
        PlantFUVTree(NImportAll);
        WalkLocalTreeAndGetLFUV(NImportAll);

        // Neighbor search
        static int Neighbors[MaxNeighborSize]; 
        const double factor = 1.0/(1.6e-3*4*M_PI*SQ(Pall.UnitLength));
        for(int i=0;i<Pall.Nhydro;i++){
            if(
#ifdef FUVFEEDBACK_G0_EVALUATION_INTERVAL //{
                    CurrentHydroActive[i]
#else //} FUVFEEDBACK_G0_EVALUATION_INTERVAL //{
                    NBCache[i].Active
#endif //} FUVFEEDBACK_G0_EVALUATION_INTERVAL //{
                    ){
                int leaf = NBCache[i].Leaf;
                Phydro[leaf]->G0thin = factor*FUVTreeResultCache[i].LFUV;
                //Phydro[leaf]->G0thin = FUVTreeResultCache[i].LFUV;

#ifdef USE_FUVFEEDBACK_LOCAL_CORRECTION //{

#define FUV_Epsilon (0.05)
#define FUV_Zsolar (0.0134)

                double X = 1.e0 - HeliumAbandance - Phydro[leaf]->Z;
                double nH = Pall.ConvertNumberDensityToCGS*Phydro[leaf]->Rho*X;

                double Ncol = EvaluateColumnDensity(leaf);

                double Zdep = (Phydro[leaf]->Z/FUV_Zsolar);

#ifdef USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //{
                double Lscale = FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE/Pall.UnitLength;
#else  //USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}//{
                double Lscale = 2*Phydro[leaf]->Kernel;
#endif   //USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}

                int Nlist = GetNeighborsFUV(PhydroBody(leaf)->PosP,Lscale,Neighbors);

                Phydro[leaf]->G0thinLocal = Phydro[leaf]->G0extLocal = 0.e0;
                for(int k=0;k<Nlist;k++){
                    int target_leaf = Neighbors[k];
                    double r2 = DISTANCE2(PhydroBody(leaf)->PosP,FUVTreeCache[target_leaf].Pos);
                    double L = FUVTreeCache[target_leaf].LFUV/r2;

                    Phydro[leaf]->G0thinLocal += L;

#       ifdef USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //{
                    double Exp = exp(-Zdep*FUVFEEDBACK_SIGMA1000*Ncol*(sqrt(r2)/Lscale)); //
#        else // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}//{
                    double Exp = exp(-FUVFEEDBACK_SIGMA1000*Ncol*(sqrt(r2)/Lscale)); //
#       endif // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}

                    // assert(sqrt(r2) < 2*Phydro[leaf]->Kernel);

                    Phydro[leaf]->G0extLocal += L*Exp;
                }
                Phydro[leaf]->G0thinLocal *= factor;
                Phydro[leaf]->G0extLocal *= factor;
#if 0
                if(Nlist > 0){
                    fprintf(stderr,"FUV correction term: %d | %g %g %g\n",leaf,
                            Phydro[leaf]->G0thin,
                            Phydro[leaf]->G0thinLocal,Phydro[leaf]->G0extLocal);
                }
#endif
#endif // USE_FUVFEEDBACK_LOCAL_CORRECTION //}
            }
        }
    }

    return;
#endif //USE_FUVFEEDBACK //}

}
