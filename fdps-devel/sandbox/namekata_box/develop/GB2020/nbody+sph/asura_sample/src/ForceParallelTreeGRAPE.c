#include "config.h"
#include "ForceMisc.h"
#include "PlantGravityTree.h"
#include "GRAPEManager.h"
#if (defined(HAVE_GRAPE7) \
  || defined(HAVE_GRAPE6A) \
  || defined(HAVE_GRAPE5) \
  || defined(HAVE_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE_API2))
#include <gp5util.h>
#else
#include "GRAPEEmulator.h"
#define HAVE_GRAPE_EMULATOR (ON)
#endif

/*
 * It may give wrong results when you turn off the flag MakeGlobalTree.  
 */
#define MakeGlobalTree  (ON)

static int Npipes;
#ifdef HAVE_GRAPE7
static int JMEMSIZE;
#endif

StructGravityRoot LETRoot;
struct StructGravityNode *LETNode;
struct StructGravityCache *LETCache;

struct StructLeavesExportImport{
    double Pos[3];
    double Mass;
#ifdef USE_SYMMETRIZED_SOFTENING //{
    double Eps;
#endif // USE_SYMMETRIZED_SOFTENING //}
};

void CCC(const int ID){
    dbg("[%02d] LETNode[%d].ok = %llu\n",MPIGetMyID(),ID,LETNode[ID].OrderingKey);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    return ;
}

void DDD(const char *func, const int line){

    for(int i=0;i<LETRoot.NumberofAllocatedNodes;i++){
        if(LETNode[i].OrderingKey != 0){
            dbg("[%02d] LETNode[%d].ok = %llu: %s:%d\n",MPIGetMyID(),i,LETNode[i].OrderingKey,func,line);
        }
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    return ;
}

static struct StructLeavesExportImport **LeavesExport;
static struct StructLeavesExportImport *LeavesImport;
static int NumberofAllocatedLeavesImport = 0;

void InitializeRootForLET(void){

    int NloadGravity =(int)(ForAngelsShare*FirstAllocationSize);
    int BaseNumberofChildren = 2;
    int Level = (int)(log((double)FirstAllocationSize)/log((double)CUBE(BaseNumberofChildren)));
    int NumberofAllocatedNodes = FirstAllocationSize*
        (1.e0-pow(1.e0/((double)CUBE(BaseNumberofChildren)),(double)Level))/
            (1.e0-1.e0/((double)CUBE(BaseNumberofChildren)));

    /* GravtityRoot */
    LETRoot.NumberofAllocatedLeaves = NloadGravity;
    LETRoot.Leaves = malloc(sizeof(int)*NloadGravity+1);
    LETRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
    LETRoot.BaseNumberofChildren = BaseNumberofChildren; // the maximum number of children is 2^3.

    LETRoot.NumberofNodeCreationLimit = TreeNodeGenerationLimitNumberForGrav;
    LETRoot.MaxLevel = TreeMaxNodeLevel;

    LETRoot.OpeningAngle = TreeOpeningAngle;
    LETRoot.NumberofLeavesInGroup = TreeNGroup;

    double dw = 1.0/BaseNumberofChildren;
    for(int i=0;i<TreeMaxNodeLevel;i++)
        LETRoot.WidthFactor[i] = pow(dw,(double)i);
    /* GravtityRoot */

    /* GravityNode */
    LETNode = malloc(sizeof(struct StructGravityNode)*LETRoot.NumberofAllocatedNodes);

    struct StructGravityNode LETNodeTemp;
    memset(&LETNodeTemp,0,sizeof(struct StructGravityNode));
    LETNodeTemp.Next = LETNodeTemp.Parent = LETNodeTemp.Children =
    LETNodeTemp.Sister = NONE;
    //LETNodeTemp.Sister = LETNodeTemp.Traverse = NONE;

    for(int i=0;i<LETRoot.NumberofAllocatedNodes;i++)
        LETNode[i] = LETNodeTemp;
    /* GravityNode */

    /* GravityCache */
    LETCache = malloc(sizeof(struct StructGravityCache)*NloadGravity+1);
    /* GravityCache */

    /* allocate communication buffer for export */
    int NProcs = MPIGetNumProcs();
    LeavesExport = malloc(sizeof(struct StructLeavesExportImport*)*NProcs);
    /* allocate communication buffer for export */

    return ;
}


void InitializeParallelTreeGRAPE(void){

#if (defined(HAVE_GRAPE7) \
  || defined(HAVE_GRAPE6A) \
  || defined(HAVE_GRAPE5) \
  || defined(HAVE_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE_API2))
	g5_open();
    Npipes = g5_get_number_of_pipelines();
#ifdef HAVE_GRAPE7
    JMEMSIZE = g5_get_jmemsize();
#endif
	g5_close();
#ifdef HAVE_PHANTOM_GRAPE
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr," =====================\n");
        fprintf(stderr," = Use Phantom GRAPE =\n");
        fprintf(stderr," =====================\n");
        fprintf(stderr," JMEMSIZE = %d\n",JMEMSIZE);
        fflush(stderr);
    }
#endif
#else
	g5_open_emu();
    Npipes = g5_get_number_of_pipelines_emu();
	g5_close_emu();
    InitializeGRAPEEmulator();
#endif

    return;
}

void UpdateNumberofLeavesInGroupInLET(const int Ng){
    LETRoot.NumberofLeavesInGroup = TreeNGroup;
    return ;
}

static void SetGRAPEScaleForParallelTreeGRAPE(void);
static int MakeExportParticleList(const int Index);
static void WalkLocalTreeAndGetAccPot(const int NImportAll);
static void GetAccPotWithoutTreeWalk(const int NImportAll);
static int GetFieldFromLocalTree(const int CurrentNodeID);
static int GetFieldFromLET(const int CurrentNodeID, const int NField);
static void CalculateForceParallelTreeGRAPEEngine(const int CurrentNodeID);
static void PlantLET(const int NumberofLeavesThisStep);

//static double SwitchingDistance;

#if 0
#define NTest 119
void PhantomTest(const int line){

    double x[NTest][3];
    double eps[NTest];
    double eps2[NTest];
    double m[NTest];
    double acc[NTest][3];
    double pot[NTest];

    for(int i=0;i<NTest;i++){
        x[i][0] = Pbody[i]->Pos[0];
        x[i][1] = Pbody[i]->Pos[1];
        x[i][2] = Pbody[i]->Pos[2];
        m[i] = Pbody[i]->Mass;
        eps[i] = Pbody[i]->Eps;
        eps2[i] = SQ(Pbody[i]->Eps);
        // fprintf(stderr,"%d %g %g %g %g %g %g\n",
                // i,x[i][0],x[i][1],x[i][2],m[i],eps[i],eps2[i]);
    }

    g5_open();
#ifdef HAVE_AVX_PHANTOM_GRAPE
    g5_set_xmj0(0,NTest,x,m,eps2);
    g5_set_n(NTest);

    g5_calculate_force_on_x0(x,acc,pot,NTest,eps2);
#else
#endif


    double sum[4] = {0.0,0.0,0.0,0.0};
    for(int i=0;i<NTest;i++){
        sum[0] += acc[i][0];
        sum[1] += acc[i][1];
        sum[2] += acc[i][2];
        sum[3] += pot[i];
        // fprintf(stderr,"A %d %g %g %g %g\n",
                // i,acc[i][0],acc[i][1],acc[i][2],pot[i]);
    }

    fprintf(stderr,"Sum %g %g %g %g %s:%d (%d)\n",
            sum[0],sum[1],sum[2],sum[3],__FUNCTION__,__LINE__,line);

    return ;
}

void SumAcc(const int line){

    double sum[4] = {0.0,0.0,0.0,0.0};
    for(int i=0;i<NTest;i++){
        sum[0] += Pbody[i]->Acc[0];
        sum[1] += Pbody[i]->Acc[1];
        sum[2] += Pbody[i]->Acc[2];
        sum[3] += Pbody[i]->Pot;
        // fprintf(stderr,"A %d %g %g %g %g\n",
                // i,acc[i][0],acc[i][1],acc[i][2],pot[i]);
    }

    fprintf(stderr,"SumR %g %g %g %g %s:%d (%d)\n",
            sum[0],sum[1],sum[2],sum[3],__FUNCTION__,__LINE__,line);

    return ;
}
#endif

void ForceParallelTreeGRAPE(void){

#ifndef GRAVITY_RUN
    return ;
#endif //GRAVITY_RUN

#ifdef TASK_KEPLER
    return ;
#endif //GRAVITY_RUN

    double TimingResultThisRoutine = GetElapsedTime();

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    HoldGRAPE();

    SetGRAPEScaleForParallelTreeGRAPE();
    for(int i=0;i<Pall.Ntotal;i++){
        GravityAccPotCache[i].Acc[0] = 
        GravityAccPotCache[i].Acc[1] = 
        GravityAccPotCache[i].Acc[2] = 
        GravityAccPotCache[i].Pot = 0.e0;
        GravityAccPotCache[i].InteractionList = 0;
    }

    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = MakeExportParticleList(i);
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
        LeavesImport = realloc(LeavesImport,NumberofAllocatedLeavesImport*sizeof(struct StructLeavesExportImport));
    }

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i] > 0){
            MPI_Isend(LeavesExport[i],
                NExportThisTime[i]*sizeof(struct StructLeavesExportImport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_FORCE_PARALLELTREEGRAPE_EXPORT,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i] > 0){
            MPI_Irecv(LeavesImport+NImport,
                NImportThisTime[i]*sizeof(struct StructLeavesExportImport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FORCE_PARALLELTREEGRAPE_EXPORT,
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

    // tree for imported data.
#if MakeGlobalTree
    NImport += Pall.Ntotal;
    if(NumberofAllocatedLeavesImport < NImport){
        NumberofAllocatedLeavesImport = ForAngelsShare*NImport;
        LeavesImport = realloc(LeavesImport,NumberofAllocatedLeavesImport*sizeof(struct StructLeavesExportImport));
    }

    int IndexShift = Pall.Ntotal; 
    for(int i=0;i<NImportAll;i++){ // move
        int BackwardID  = NImportAll-1-i;
        LeavesImport[BackwardID+IndexShift] = LeavesImport[BackwardID];
    }

    for(int i=0;i<Pall.Ntotal;i++){ // copy
        LeavesImport[i].Pos[0] = GravityCache[i].Pos[0];
        LeavesImport[i].Pos[1] = GravityCache[i].Pos[1];
        LeavesImport[i].Pos[2] = GravityCache[i].Pos[2];
        LeavesImport[i].Mass   = GravityCache[i].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
        LeavesImport[i].Eps    = GravityCache[i].Eps;
        // assert(LeavesImport[i].Eps > 0.e0);
#endif
    }
    NImportAll = NImport;
#endif

    if(NImportAll>0)
        PlantLET(NImportAll);

    //double tt = GetElapsedTime();

    WalkLocalTreeAndGetAccPot(NImportAll);

    //double tt2 = GetElapsedTime();

#if 0

    if(MPIGetMyID()==MPI_ROOT_RANK){
        MakeDir("./ppp");
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<MPIGetNumProcs();i++){
        if(i==MPIGetMyID()){
            fprintf(stderr,"%3ld %7d (%7d) %7d %g\n",MPIGetMyID(),
                    Pall.Ntotal,Pall.Nhydro,NImportAll,
                    tt2-tt);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    {
    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./ppp/in.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%1.8g %1.8g %1.8g\n",
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2]);

    }
    fclose(fp);

    Snprintf(fname,"./ppp/out.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<NImportAll;i++){
        fprintf(fp,"%1.8g %1.8g %1.8g\n",
            LeavesImport[i+Pall.Ntotal].Pos[0],
                LeavesImport[i+Pall.Ntotal].Pos[1],
                    LeavesImport[i+Pall.Ntotal].Pos[2]);

    }
    fclose(fp);


    fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    ReleaseGRAPE();

#ifdef USE_SYMMETRIZED_SOFTENING
    double SymmetrizedFactor = 1.0/sqrt(2.0);
#endif
    double AdaptiveSofteningFactor = Pall.AdaptiveSofteningFactor;
    for(int i=0;i<NImportAll;i++){
        if(LETCache[i].Active){
            int leaf = LETCache[i].Leaf;
            Pbody[leaf]->Acc[0] = Pall.GravConst*GravityAccPotCache[leaf].Acc[0];
            Pbody[leaf]->Acc[1] = Pall.GravConst*GravityAccPotCache[leaf].Acc[1];
            Pbody[leaf]->Acc[2] = Pall.GravConst*GravityAccPotCache[leaf].Acc[2];
#ifdef HAVE_GRAPE7
            Pbody[leaf]->Pot = 0.e0;
#else //HAVE_GRAPE7
#ifdef USE_SYMMETRIZED_SOFTENING
            Pbody[leaf]->Pot = GravityAccPotCache[leaf].Pot 
                - SymmetrizedFactor*(Pbody[leaf]->Mass/(AdaptiveSofteningFactor*Pbody[leaf]->Eps));
#else //USE_SYMMETRIZED_SOFTENING
            Pbody[leaf]->Pot = GravityAccPotCache[leaf].Pot 
                - (Pbody[leaf]->Mass/(AdaptiveSofteningFactor*Pbody[leaf]->Eps));
#endif///USE_SYMMETRIZED_SOFTENING
            Pbody[leaf]->Pot = -0.5*Pall.GravConst*Pbody[leaf]->Mass*Pbody[leaf]->Pot;
#endif //HAVE_GRAPE7
            Pbody[leaf]->InteractionList = GravityAccPotCache[leaf].InteractionList;
            //fprintf(stderr,"%g %g %g %g\n",Pbody[leaf]->Acc[0],Pbody[leaf]->Acc[1],Pbody[leaf]->Acc[2],Pbody[leaf]->Pot);
        }
    }

    ForceEndProcedure();

    TimingResults.GravityThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

static void SetGRAPEScaleForParallelTreeGRAPE(void){

#if (defined(HAVE_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE_API2))
    return ;
#endif

#define MinimumMassFactor    (MINIMUM_MASS_FACTOR_FOR_GRAPE)
    const static double InvMinimumMassFactor = 1.e0/MinimumMassFactor;

	double xmin,ymin,zmin,xmax,ymax,zmax;
	double rmin,mmin,size;

	xmin = xmax = GravityCache[0].Pos[0];
	ymin = ymax = GravityCache[0].Pos[1];
	zmin = zmax = GravityCache[0].Pos[2];
	mmin = GravityCache[0].Mass;
    for(int i=0;i<Pall.Ntotal;i++){
		xmin = fmin(xmin,GravityCache[i].Pos[0]); xmax = fmax(xmax,GravityCache[i].Pos[0]);
		ymin = fmin(ymin,GravityCache[i].Pos[1]); ymax = fmax(ymax,GravityCache[i].Pos[1]);
		zmin = fmin(zmin,GravityCache[i].Pos[2]); zmax = fmax(zmax,GravityCache[i].Pos[2]);
		mmin = fmin(mmin,GravityCache[i].Mass);
	}
    assert(mmin>0.e0);

    double PartialMinMaxMass[3],GlobalMinMaxMass[3];
    PartialMinMaxMass[0] = -fmin(fmin(xmin,ymin),zmin);
    PartialMinMaxMass[1] = +fmax(fmax(xmax,ymax),zmax);
    PartialMinMaxMass[2] = -mmin;
    MPI_Allreduce(PartialMinMaxMass,GlobalMinMaxMass,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    double MaximumPos = GlobalMinMaxMass[0];
    double MinimumPos = -GlobalMinMaxMass[1];
    double MinimumMass = -GlobalMinMaxMass[2];

	size = 2.01*fmax(fabs(MaximumPos),fabs(MinimumPos));
	rmin = 1.e-7*size;
	MinimumMass *= InvMinimumMassFactor; 
    mmin = MinimumMass;

    assert(MinimumMass > 0.e0);

    double AdaptiveSofteningFactor = Pall.AdaptiveSofteningFactor;
    for(int i=0;i<Pall.Ntotal;i++){
        if(SQ(AdaptiveSofteningFactor*GravityCache[i].Eps) < 
                SQ(rmin)*(GravityCache[i].Mass*InvMinimumMassFactor)){
            int leaf = GravityCache[i].Leaf;
            fprintf(stderr,"[%02d] The scaling Error is occurred!\n",MPIGetMyID());
            fprintf(stderr,"[%02d] function %s\n",MPIGetMyID(),__FUNCTION__);
            fprintf(stderr,"Index = %d\n",leaf);
            fprintf(stderr,"TYPE = %d\n",Pbody[leaf]->Type);
            fprintf(stderr,"Mass/128 = %g\n",Pbody[leaf]->Mass/128.0);
            fprintf(stderr,"ASF = %g\n",Pall.AdaptiveSofteningFactor);
            fprintf(stderr,"Eps = %g\n",Pbody[leaf]->Eps);
            fprintf(stderr,"ASF*Eps = %g\n",Pall.AdaptiveSofteningFactor*Pbody[leaf]->Eps);
            fprintf(stderr,"mmin = %g\n",mmin);
            fprintf(stderr,"rmin = %g\n",rmin);
            fprintf(stderr,"rmin = %g, sqrt(m[i]/mmin) = %g, rmin*sqrt(mi/mmin) = %g\n",
                    rmin,sqrt(Pbody[leaf]->Mass/mmin),rmin*sqrt(Pbody[leaf]->Mass/mmin));

            fprintf(stderr,"max min= %g, %g\n",MaximumPos,MinimumPos);

            if(Pbody[leaf]->Type == TypeDM){
                dprintlmpi(i);
                StructureReportPbody(leaf);
                StructureReportPbody(GravityCache[i-1].Leaf);
            } else if(Pbody[leaf]->Type == TypeHydro){
                StructureReportPhydro(leaf);
            } else if(Pbody[leaf]->Type == TypeStar){
                StructureReportPstar(leaf);
            } else {
                fprintf(stderr,"Undefined type.\n");
            }

            abort();
            MPI_Finalize();
            exit(GRAPEScaleError);
        }
    }

    //fprintf(stderr,"size = %g, mmin = %g\n",size,mmin);
#ifdef HAVE_GRAPE_EMULATOR //{
	g5_set_range_emu(-size,size,mmin);
#else //HAVE_GRAPE_EMULATOR //}//{
#if (!defined(HAVE_PHANTOM_GRAPE) \
  && !defined(HAVE_AVX_PHANTOM_GRAPE) \
  && !defined(HAVE_AVX_PHANTOM_GRAPE_API2)) //{
	g5_set_range(-size,size,mmin);
#endif //}
#endif //HAVE_GRAPE_EMULATOR //}

	return;
}


static int MakeExportParticleList(const int Index){

    if(Pall.Ntotal == 0){
        return 0;
    }

    int ExportNodeID = CommunicationTable[Index].SendRank;

    double PosMax[3] = {EdgesForGravity[ExportNodeID].PosMax[0],
                        EdgesForGravity[ExportNodeID].PosMax[1],
                        EdgesForGravity[ExportNodeID].PosMax[2]};
    double PosMin[3] = {EdgesForGravity[ExportNodeID].PosMin[0],
                        EdgesForGravity[ExportNodeID].PosMin[1],
                        EdgesForGravity[ExportNodeID].PosMin[2]};

    double Pos[3] = {0.5*(PosMax[0]+PosMin[0]),
                     0.5*(PosMax[1]+PosMin[1]),
                     0.5*(PosMax[2]+PosMin[2])};
    double hWidth[3] = {PosMax[0]-Pos[0],PosMax[1]-Pos[1],PosMax[2]-Pos[2]};
    if(hWidth[0] == 0.e0){
        return 0;
    }
    assert(hWidth[0] > 0.e0);
    assert(hWidth[1] > 0.e0);
    assert(hWidth[2] > 0.e0);

    //fprintf(stderr," hwidth %g %g %g\n",hWidth[0],hWidth[1],hWidth[2]);

    int SendThisTime = 0;
    double theta2 = SQ(GravityRoot.OpeningAngle);
    int RootNodeID = 0;
    int CurrentNodeID = GravityNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
        double x[3] = {fabs(GravityNode[CurrentNodeID].Pos[0]-Pos[0]),
                       fabs(GravityNode[CurrentNodeID].Pos[1]-Pos[1]),
                       fabs(GravityNode[CurrentNodeID].Pos[2]-Pos[2])};
        double dx2 = 0.e0;
        double hwidth = 0.5*GravityRoot.Width*GravityRoot.WidthFactor[GravityNode[CurrentNodeID].Level];
        if(x[0] > (hWidth[0]+hwidth)) dx2 += SQ(x[0]-(hWidth[0]+hwidth));
        if(x[1] > (hWidth[1]+hwidth)) dx2 += SQ(x[1]-(hWidth[1]+hwidth));
        if(x[2] > (hWidth[2]+hwidth)) dx2 += SQ(x[2]-(hWidth[2]+hwidth));


        if(SQ(2.0*hwidth) < theta2*dx2){
            CheckSizeofBufferExportSendIndex((SendThisTime+1),
                    sizeof(struct StructLeavesExportImport),Index);
            LeavesExport[Index] = BufferExportSend[Index];

            LeavesExport[Index][SendThisTime].Pos[0] = GravityNode[CurrentNodeID].COM[0];
            LeavesExport[Index][SendThisTime].Pos[1] = GravityNode[CurrentNodeID].COM[1];
            LeavesExport[Index][SendThisTime].Pos[2] = GravityNode[CurrentNodeID].COM[2];
            LeavesExport[Index][SendThisTime].Mass   = GravityNode[CurrentNodeID].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
            LeavesExport[Index][SendThisTime].Eps    = sqrt(GravityNode[CurrentNodeID].Eps2);
#endif // USE_SYMMETRIZED_SOFTENING

            SendThisTime ++;
            CurrentNodeID = GravityNode[CurrentNodeID].Next;
		} else if (GravityNode[CurrentNodeID].Children == NONE) {
            int NumberofLeaves = GravityNode[CurrentNodeID].NumberofLeaves;
            CheckSizeofBufferExportSendIndex((SendThisTime+NumberofLeaves),
                    sizeof(struct StructLeavesExportImport),Index);
            LeavesExport[Index] = BufferExportSend[Index];

            int header = GravityNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                LeavesExport[Index][SendThisTime].Pos[0] = GravityCache[leaf].Pos[0];
                LeavesExport[Index][SendThisTime].Pos[1] = GravityCache[leaf].Pos[1];
                LeavesExport[Index][SendThisTime].Pos[2] = GravityCache[leaf].Pos[2];
                LeavesExport[Index][SendThisTime].Mass   = GravityCache[leaf].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
                LeavesExport[Index][SendThisTime].Eps    = GravityCache[leaf].Eps;
#endif // USE_SYMMETRIZED_SOFTENING
                SendThisTime ++;
            }
			CurrentNodeID = GravityNode[CurrentNodeID].Next;
		} else {
			CurrentNodeID = GravityNode[CurrentNodeID].Children;
		}
	}

    return SendThisTime;
}

static int NumberofFieldAllocated = 0;

static double (*FieldPos)[3]; 
static double *FieldMass;
static double *FieldEps2;

static void WalkLocalTreeAndGetAccPot(const int NImportAll){

    if(NumberofFieldAllocated == 0){
        NumberofFieldAllocated = FirstAllocationSize;
        FieldPos = malloc(sizeof(double)*3*NumberofFieldAllocated);
        FieldMass = malloc(sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
        FieldEps2 = malloc(sizeof(double)*NumberofFieldAllocated);
#endif
    }

    int counter = 0;

    int RootNodeID = 0;
    int CurrentNodeID = LETNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
		if(LETNode[CurrentNodeID].NumberofActiveLeaves == 0){ 
            CurrentNodeID = LETNode[CurrentNodeID].Next;
		}else if( (LETNode[CurrentNodeID].NumberofLeaves < GravityRoot.NumberofLeavesInGroup)
                || (LETNode[CurrentNodeID].Children == NONE)){
#if MakeGlobalTree
            int NumberofField = GetFieldFromLET(CurrentNodeID,0);
#else
            int NumberofField = GetFieldFromLocalTree(CurrentNodeID);
            if(NImportAll>0)
                NumberofField = GetFieldFromLET(CurrentNodeID,NumberofField);
#endif
            counter ++;

            int NumberofLeaves = LETNode[CurrentNodeID].NumberofLeaves;
            int header = LETNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(LETCache[leaf].Active)
                    GravityAccPotCache[LETCache[leaf].Leaf].InteractionList = NumberofField;
            }

            for(int j=0;j<NumberofField;j+=JMEMSIZE){ 
                int Active_j = MIN(JMEMSIZE,NumberofField-j);
#ifdef USE_SYMMETRIZED_SOFTENING //{
                g5_set_n(Active_j);
#ifdef HAVE_AVX_PHANTOM_GRAPE //}//{
                g5_set_xmj0(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#elif  HAVE_AVX_PHANTOM_GRAPE_API2 //}//{
                g5_set_xmj(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
                g5_set_xmeps2j(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else  // USE_SYMMETRIZED_SOFTENING //} //{
#if defined(HAVE_GRAPE7)
                g5_set_n(Active_j);
                g5_set_jp(0,Active_j,FieldMass+j,FieldPos+j);
#elif (defined(HAVE_GRAPE6A) \
    || defined(HAVE_GRAPE5) \
    || defined(HAVE_PHANTOM_GRAPE) \
    || defined(HAVE_AVX_PHANTOM_GRAPE) \
    || defined(HAVE_AVX_PHANTOM_GRAPE_API2))
                g5_set_n(Active_j);
                g5_set_xmj(0,Active_j,FieldPos+j,FieldMass+j);
#else
                g5_set_n_emu(Active_j);
                g5_set_xmj_emu(0,Active_j,FieldPos+j,FieldMass+j);
#endif
#endif // USE_SYMMETRIZED_SOFTENING //}
                CalculateForceParallelTreeGRAPEEngine(CurrentNodeID);
            }
            CurrentNodeID = LETNode[CurrentNodeID].Next;
		} else {
		    CurrentNodeID = LETNode[CurrentNodeID].Children;
		}
	}

    return;
}

static int GetFieldFromLETForGetAccPotWithoutTreeWalk(double Position[restrict], const double Eps2, const int NField){

    int NumberofField = NField;
    double theta2 = SQ(GravityRoot.OpeningAngle);

    int RootNodeID = 0;
	int TargetNodeID = LETNode[RootNodeID].Children;
	while(TargetNodeID != RootNodeID){ 
        double sqDist = 0.e0;
        double Dist[3] = {fabs(Position[0]-LETNode[TargetNodeID].COM[0]),
                          fabs(Position[1]-LETNode[TargetNodeID].COM[1]),
                          fabs(Position[2]-LETNode[TargetNodeID].COM[2])};
        sqDist = NORM2(Dist);
        double sqDist3D = sqDist;
        sqDist += Eps2;
#ifdef USE_SYMMETRIZED_SOFTENING
        sqDist += LETNode[TargetNodeID].Eps2;
#endif // USE_SYMMETRIZED_SOFTENING

        double width = LETRoot.Width*LETRoot.WidthFactor[LETNode[TargetNodeID].Level];
#ifdef USE_SYMMETRIZED_SOFTENING
        width = sqrt(SQ(width) + SQ(LETNode[TargetNodeID].EpsMax-LETNode[TargetNodeID].EpsMin));
        //width += SQ(LETNode[TargetNodeID].EpsMax-LETNode[TargetNodeID].EpsMin);
#endif // USE_SYMMETRIZED_SOFTENING

#ifdef USE_SYMMETRIZED_SOFTENING
        if( (SQ(width)<theta2*sqDist) && (sqDist3D > TINY)) { // For self-interactions.
#else
        if(SQ(width)<theta2*sqDist){
#endif // USE_SYMMETRIZED_SOFTENING
            if(NumberofField+1 > NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+1));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif // USE_SYMMETRIZED_SOFTENING
            }
#ifdef USE_SHIFT_GRAVITYFIELD
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0]-Position[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1]-Position[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2]-Position[2];
#else 
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2];
#endif // USE_SHIFT_GRAVITYFIELD
            FieldMass[NumberofField]   = LETNode[TargetNodeID].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
            FieldEps2[NumberofField]   = LETNode[TargetNodeID].Eps2;
#endif // USE_SYMMETRIZED_SOFTENING
    
            NumberofField ++;
			TargetNodeID = LETNode[TargetNodeID].Next;
#ifdef USE_SYMMETRIZED_SOFTENING
		 } else if ((LETNode[TargetNodeID].Children == NONE) || (sqDist3D < TINY)) { // For self-interactions.
#else
		 } else if (LETNode[TargetNodeID].Children == NONE) {
#endif // USE_SYMMETRIZED_SOFTENING
            int NumberofLeaves = LETNode[TargetNodeID].NumberofLeaves;
            if(NumberofField+NumberofLeaves > NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+NumberofLeaves));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif
            }
            int header = LETNode[TargetNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header + k;
#ifdef USE_SHIFT_GRAVITYFIELD
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0]-Position[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1]-Position[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2]-Position[2];
#else
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2];
#endif
                FieldMass[NumberofField]   = LETCache[leaf].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2[NumberofField]   = SQ(LETCache[leaf].Eps);
#endif
                NumberofField ++;
            }
			TargetNodeID = LETNode[TargetNodeID].Next;
		 } else {
			TargetNodeID = LETNode[TargetNodeID].Children;
		 }
	}

	return (NumberofField);
}

static void CalculateForceParallelTreeGRAPEEngineForAccPotWithoutTreeWalk(const int CurrentNodeID){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    int LeavesLog[Npipes];

    int NumberofActiveLeaves = LETNode[CurrentNodeID].NumberofActiveLeaves;
    int header = LETNode[CurrentNodeID].Leaves;
    int CurrentLeafID = 0;
    for(int i=0;i<NumberofActiveLeaves;i+=Npipes){
		int NActives = MIN(Npipes,NumberofActiveLeaves-i);
        int Count = 0;
        //int leaf = header+CurrentLeafID;
        while(Count!=NActives){
            int leaf = header + CurrentLeafID;
            if(LETCache[leaf].Active){
#ifdef USE_SHIFT_GRAVITYFIELD
#if MakeGlobalTree
                Dxi[Count][0] = LETCache[leaf].Pos[0]-LETNode[CurrentNodeID].COM[0];
                Dxi[Count][1] = LETCache[leaf].Pos[1]-LETNode[CurrentNodeID].COM[1];
                Dxi[Count][2] = LETCache[leaf].Pos[2]-LETNode[CurrentNodeID].COM[2];
#else
                Dxi[Count][0] = LETCache[leaf].Pos[0]-GravityNode[CurrentNodeID].COM[0];
                Dxi[Count][1] = LETCache[leaf].Pos[1]-GravityNode[CurrentNodeID].COM[1];
                Dxi[Count][2] = LETCache[leaf].Pos[2]-GravityNode[CurrentNodeID].COM[2];
#endif
#else
                Dxi[Count][0] = LETCache[leaf].Pos[0];
                Dxi[Count][1] = LETCache[leaf].Pos[1];
                Dxi[Count][2] = LETCache[leaf].Pos[2];
#endif
#if (defined(HAVE_AVX_PHANTOM_GRAPE)\
  || defined(HAVE_AVX_PHANTOM_GRAPE_API2)) //{
                Depsi[Count] = SQ(LETCache[leaf].Eps);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
                Depsi[Count] = LETCache[leaf].Eps;
#endif // HAVE_AVX_PHANTOM_GRAPE //}
                LeavesLog[Count] = LETCache[leaf].Leaf;
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

#ifdef USE_SYMMETRIZED_SOFTENING
#ifdef HAVE_AVX_PHANTOM_GRAPE //{
        g5_calculate_force_on_x0(Dxi,adummy,phidummy,Npipes,Depsi);
#elif  HAVE_AVX_PHANTOM_GRAPE_API2 //}//{
        g5_calculate_force_on_xe(Dxi,Depsi,adummy,phidummy,Npipes);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
        g5_set_xi(Npipes,Dxi);
        g5_set_eps(Npipes,Depsi);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else 
#if defined(HAVE_GRAPE7)
        g5_set_eps2_to_all(SQ(Depsi[0]));
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE6A)
        g5_set_xepsi(Npipes,Dxi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE5)
        g5_set_ip(Npipes,Dxi,Depsi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif (defined(HAVE_PHANTOM_GRAPE) \
   || defined(HAVE_AVX_PHANTOM_GRAPE) \
   || defined(HAVE_AVX_PHANTOM_GRAPE_API2))
        g5_set_eps_to_all(Depsi[0]);
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#else
        g5_set_ip_emu(Npipes,Dxi,Depsi,Depsi);
        g5_run_emu();
        g5_get_force_emu(Npipes,adummy,phidummy);
#endif
#endif

        // back acc and pot for the temporal data array.
        for(int k=0;k<NActives;k++){
			GravityAccPotCache[LeavesLog[k]].Acc[0] += adummy[k][0];
			GravityAccPotCache[LeavesLog[k]].Acc[1] += adummy[k][1];
			GravityAccPotCache[LeavesLog[k]].Acc[2] += adummy[k][2];
			GravityAccPotCache[LeavesLog[k]].Pot += phidummy[k];
            /*
            if(!isnormal(adummy[k][0])){
                fprintf(stderr,"L[%d] = %d; %d:%d:%d\n",k,LeavesLog[k],i,NActives,NumberofActiveLeaves);
                fflush(NULL);
                assert(isnormal(adummy[k][0]));
            }
            */
        }
    }


	return;
}

extern int MassRatioStart,MassRatioEnd;
static void GetAccPotWithoutTreeWalk(const int NImportAll){

    if(NumberofFieldAllocated == 0){
        NumberofFieldAllocated = FirstAllocationSize;
        FieldPos = malloc(sizeof(double)*3*NumberofFieldAllocated);
        FieldMass = malloc(sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
        FieldEps2 = malloc(sizeof(double)*NumberofFieldAllocated);
#endif
    }

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            int NumberofField = 
                GetFieldFromLETForGetAccPotWithoutTreeWalk(Pbody[i]->PosP,SQ(Pbody[i]->Eps),0);

            //Pbody[i]->Pot = 0.e0;
            //Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;

            for(int k=0;k<NumberofField;k++){
#ifdef USE_SHIFT_GRAVITYFIELD
                double Pos[] = {0.e0,0.e0,0.e0};
#else 
                double Pos[] = {Pbody[i]->PosP[0],Pbody[i]->PosP[1],Pbody[i]->PosP[2]};
#endif
                double X[] = {Pos[0]-FieldPos[k][0],
                              Pos[1]-FieldPos[k][1],
                              Pos[2]-FieldPos[k][2]};

                double Distance2 = NORM2(X);
                Distance2 += SQ(Pbody[i]->Eps);
                Distance2 += FieldEps2[k];
                double Distance = sqrt(Distance2);
                if(NORM2(X) > 0.e0)
                    Pbody[i]->Pot += FieldMass[k]/Distance;
                Pbody[i]->Acc[0] -= (FieldMass[k]/Distance2) * (X[0]/Distance);
                Pbody[i]->Acc[1] -= (FieldMass[k]/Distance2) * (X[1]/Distance);
                Pbody[i]->Acc[2] -= (FieldMass[k]/Distance2) * (X[2]/Distance);

                //fprintf(stderr,"P %g M %g D %g \n",Pbody[i]->Pot,FieldMass[k],Distance);

                /*
                GravityAccPotCache[i].Pot += FieldMass[k]/Distance;
                GravityAccPotCache[i].Acc[0] -= (FieldMass[k]/Distance2) * (X[0]/Distance);
                GravityAccPotCache[i].Acc[1] -= (FieldMass[k]/Distance2) * (X[1]/Distance);
                GravityAccPotCache[i].Acc[2] -= (FieldMass[k]/Distance2) * (X[2]/Distance);
                */
            }
#if 0
            FILE *fp;
            char Fname[MaxCharactersInLine];
            sprintf(Fname,"./Field.%ld.%g.%d.%d",Pall.Ntotal_t,GravityRoot.OpeningAngle,MassRatioStart,MassRatioEnd);
            FileOpen(fp,Fname,"w");
            for(int k=0;k<NumberofField;k++){
                fprintf(fp,"%g %g %g %g %g\n",FieldPos[k][0],FieldPos[k][1],FieldPos[k][2],FieldMass[k],FieldEps2[k]);
            }
            fclose(fp);
#endif

            Pbody[i]->Pot *= -Pall.GravConst;
            Pbody[i]->Acc[0] *= Pall.GravConst;
            Pbody[i]->Acc[1] *= Pall.GravConst;
            Pbody[i]->Acc[2] *= Pall.GravConst;
        }
    }


    return;
}

static int GetFieldFromLocalTree(const int CurrentNodeID){

    int NumberofField = 0;

    double cwidth = 0.5*GravityRoot.Width*GravityRoot.WidthFactor[GravityNode[CurrentNodeID].Level];
    double theta2 = SQ(GravityRoot.OpeningAngle);

    int RootNodeID = 0;
	int TargetNodeID = GravityNode[RootNodeID].Children;
	while(TargetNodeID != RootNodeID){ 
        double twidth = 0.5*GravityRoot.Width*GravityRoot.WidthFactor[GravityNode[TargetNodeID].Level];
        double sqDist = 0.e0;
        double Dist[3] = {fabs(GravityNode[CurrentNodeID].Pos[0]-GravityNode[TargetNodeID].Pos[0]),
                          fabs(GravityNode[CurrentNodeID].Pos[1]-GravityNode[TargetNodeID].Pos[1]),
                          fabs(GravityNode[CurrentNodeID].Pos[2]-GravityNode[TargetNodeID].Pos[2])};
        double hwidth = cwidth + twidth;
        if(Dist[0]>hwidth) sqDist += SQ(Dist[0]-hwidth);
        if(Dist[1]>hwidth) sqDist += SQ(Dist[1]-hwidth);
        if(Dist[2]>hwidth) sqDist += SQ(Dist[2]-hwidth);

       if(SQ(2.0*twidth)<theta2*sqDist){
            if(NumberofField+1 > NumberofFieldAllocated){
                while(NumberofField+1 > NumberofFieldAllocated)
                    NumberofFieldAllocated = (int)(ForAngelsShare*NumberofFieldAllocated);
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
            }
#ifdef USE_SHIFT_GRAVITYFIELD
            FieldPos[NumberofField][0] = GravityNode[TargetNodeID].COM[0]-GravityNode[CurrentNodeID].COM[0];
            FieldPos[NumberofField][1] = GravityNode[TargetNodeID].COM[1]-GravityNode[CurrentNodeID].COM[1];
            FieldPos[NumberofField][2] = GravityNode[TargetNodeID].COM[2]-GravityNode[CurrentNodeID].COM[2];
#else
            FieldPos[NumberofField][0] = GravityNode[TargetNodeID].COM[0];
            FieldPos[NumberofField][1] = GravityNode[TargetNodeID].COM[1];
            FieldPos[NumberofField][2] = GravityNode[TargetNodeID].COM[2];
#endif
            FieldMass[NumberofField]   = GravityNode[TargetNodeID].Mass;

            NumberofField ++;
			TargetNodeID = GravityNode[TargetNodeID].Next;
		} else if (GravityNode[TargetNodeID].Children == NONE) {
            int NumberofLeaves = GravityNode[TargetNodeID].NumberofLeaves;
            if(NumberofField+NumberofLeaves > NumberofFieldAllocated){
                while(NumberofField+NumberofLeaves > NumberofFieldAllocated)
                    NumberofFieldAllocated = (int)(ForAngelsShare*NumberofFieldAllocated);
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
            }
            int header = GravityNode[TargetNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header + k;
#ifdef USE_SHIFT_GRAVITYFIELD
                FieldPos[NumberofField][0] = GravityCache[leaf].Pos[0]-GravityNode[CurrentNodeID].COM[0];
                FieldPos[NumberofField][1] = GravityCache[leaf].Pos[1]-GravityNode[CurrentNodeID].COM[1];
                FieldPos[NumberofField][2] = GravityCache[leaf].Pos[2]-GravityNode[CurrentNodeID].COM[2];
#else
                FieldPos[NumberofField][0] = GravityCache[leaf].Pos[0];
                FieldPos[NumberofField][1] = GravityCache[leaf].Pos[1];
                FieldPos[NumberofField][2] = GravityCache[leaf].Pos[2];
#endif
                FieldMass[NumberofField]   = GravityCache[leaf].Mass;

                NumberofField ++;
            }
			TargetNodeID = GravityNode[TargetNodeID].Next;
		} else {
			TargetNodeID = GravityNode[TargetNodeID].Children;
		}
	}

	 return (NumberofField);
}

#if 0
static int GetFieldFromLET(const int CurrentNodeID, const int NField){

    int NumberofField = NField;
    double cwidth = 0.5*GravityRoot.Width*GravityRoot.WidthFactor[LETNode[CurrentNodeID].Level];
    double theta2 = SQ(GravityRoot.OpeningAngle);
#ifdef USE_SYMMETRIZED_SOFTENING
    double cewidth = 0.5*(LETNode[CurrentNodeID].EpsMax-LETNode[CurrentNodeID].EpsMin);
    double cEpsCenter = 0.5*(LETNode[CurrentNodeID].EpsMax+LETNode[CurrentNodeID].EpsMin);
#endif

    int RootNodeID = 0;
	int TargetNodeID = LETNode[RootNodeID].Children;
	while(TargetNodeID != RootNodeID){ 
        double twidth = 0.5*LETRoot.Width*LETRoot.WidthFactor[LETNode[TargetNodeID].Level];
        double sqDist = 0.e0;
#ifdef USE_SYMMETRIZED_SOFTENING
        double tewidth = 0.5*(LETNode[TargetNodeID].EpsMax-LETNode[TargetNodeID].EpsMin);
        double tEpsCenter = 0.5*(LETNode[TargetNodeID].EpsMax+LETNode[TargetNodeID].EpsMin);
        double Dist[5] = {fabs(LETNode[CurrentNodeID].Pos[0]-LETNode[TargetNodeID].Pos[0]),
                          fabs(LETNode[CurrentNodeID].Pos[1]-LETNode[TargetNodeID].Pos[1]),
                          fabs(LETNode[CurrentNodeID].Pos[2]-LETNode[TargetNodeID].Pos[2]),
                          cEpsCenter,tEpsCenter};
        double hwidth = cwidth + twidth;
        //double hEpsWidth = cewidth + tewidth;
        sqDist += SQ(fmax(Dist[0]-hwidth,0.e0));
        sqDist += SQ(fmax(Dist[1]-hwidth,0.e0));
        sqDist += SQ(fmax(Dist[2]-hwidth,0.e0));
        double sqDist3D = sqDist;
        sqDist += SQ(fmax(Dist[3]-cewidth,0.e0));
        sqDist += SQ(fmax(Dist[4]-tewidth,0.e0));
        // Update twidth
        //twidth = 0.5*sqrt(12*SQ(twidth)+SQ(LETNode[CurrentNodeID].EpsMax-LETNode[CurrentNodeID].EpsMin)); //
        //double twidth2 = 12*SQ(twidth)+SQ(LETNode[CurrentNodeID].EpsMax-LETNode[CurrentNodeID].EpsMin); //
        double twidth2 = 12*SQ(twidth)+SQ(LETNode[TargetNodeID].EpsMax-LETNode[TargetNodeID].EpsMin); //
#else // undef USE_SYMMETRIZED_SOFTENING
        double Dist[3] = {fabs(LETNode[CurrentNodeID].Pos[0]-LETNode[TargetNodeID].Pos[0]),
                          fabs(LETNode[CurrentNodeID].Pos[1]-LETNode[TargetNodeID].Pos[1]),
                          fabs(LETNode[CurrentNodeID].Pos[2]-LETNode[TargetNodeID].Pos[2])};
        double hwidth = cwidth + twidth;
        sqDist += SQ(fmax(Dist[0]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[1]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[2]-hwidth,0.0));
#endif // USE_SYMMETRIZED_SOFTENING

#ifdef USE_SYMMETRIZED_SOFTENING
        if( (twidth2<theta2*sqDist) && (sqDist3D > TINY)) { // For self-interactions.
#else
        if(SQ(2.0*twidth)<theta2*sqDist){
#endif // USE_SYMMETRIZED_SOFTENING
            if(NumberofField+1 > NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+1));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif // USE_SYMMETRIZED_SOFTENING
            }
#ifdef USE_SHIFT_GRAVITYFIELD
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0]-LETNode[CurrentNodeID].COM[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1]-LETNode[CurrentNodeID].COM[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2]-LETNode[CurrentNodeID].COM[2];
#else 
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2];
#endif // USE_SHIFT_GRAVITYFIELD
            FieldMass[NumberofField]   = LETNode[TargetNodeID].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
            FieldEps2[NumberofField]   = LETNode[TargetNodeID].Eps2;
#endif // USE_SYMMETRIZED_SOFTENING
    
            NumberofField ++;
			TargetNodeID = LETNode[TargetNodeID].Next;
#ifdef USE_SYMMETRIZED_SOFTENING
		 } else if ((LETNode[TargetNodeID].Children == NONE) || (sqDist3D < TINY)) { // For self-interactions.
#else
		 } else if (LETNode[TargetNodeID].Children == NONE) {
#endif // USE_SYMMETRIZED_SOFTENING
            int NumberofLeaves = LETNode[TargetNodeID].NumberofLeaves;
            if(NumberofField+NumberofLeaves > NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+NumberofLeaves));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif
            }
            int header = LETNode[TargetNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header + k;
#ifdef USE_SHIFT_GRAVITYFIELD
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0]-LETNode[CurrentNodeID].COM[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1]-LETNode[CurrentNodeID].COM[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2]-LETNode[CurrentNodeID].COM[2];
#else
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2];
#endif
                FieldMass[NumberofField]   = LETCache[leaf].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2[NumberofField]   = SQ(LETCache[leaf].Eps);
#endif
                NumberofField ++;
            }
			TargetNodeID = LETNode[TargetNodeID].Next;
		 } else {
			TargetNodeID = LETNode[TargetNodeID].Children;
		 }
	}

	return (NumberofField);
}
#else  //

#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
extern long int NumberofPPInteractions;
extern long int NumberofPCInteractions;
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR


static int GetFieldFromLET(const int CurrentNodeID, const int NField){

    int NumberofField = NField;
    double cwidth = 0.5*GravityRoot.Width*GravityRoot.WidthFactor[LETNode[CurrentNodeID].Level];
    double theta2 = SQ(GravityRoot.OpeningAngle);

    int RootNodeID = 0;
	int TargetNodeID = LETNode[RootNodeID].Children;
	while(TargetNodeID != RootNodeID){ 
        double twidth = 0.5*LETRoot.Width*LETRoot.WidthFactor[LETNode[TargetNodeID].Level];
        double sqDist = 0.e0;

#if 1
        double Dist[3] = {fabs(LETNode[CurrentNodeID].Pos[0]-LETNode[TargetNodeID].COM[0]),
                          fabs(LETNode[CurrentNodeID].Pos[1]-LETNode[TargetNodeID].COM[1]),
                          fabs(LETNode[CurrentNodeID].Pos[2]-LETNode[TargetNodeID].COM[2])};
        sqDist += SQ(fmax(Dist[0],0.0));
        sqDist += SQ(fmax(Dist[1],0.0));
        sqDist += SQ(fmax(Dist[2],0.0));
#else
        double Dist[3] = {fabs(LETNode[CurrentNodeID].Pos[0]-LETNode[TargetNodeID].Pos[0]),
                          fabs(LETNode[CurrentNodeID].Pos[1]-LETNode[TargetNodeID].Pos[1]),
                          fabs(LETNode[CurrentNodeID].Pos[2]-LETNode[TargetNodeID].Pos[2])};
        double hwidth = cwidth + twidth;
        sqDist += SQ(fmax(Dist[0]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[1]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[2]-hwidth,0.0));
#endif

#ifdef USE_SYMMETRIZED_SOFTENING
        double sqDist3D = sqDist;
        sqDist += SQ(LETNode[CurrentNodeID].EpsMin);
        sqDist += LETNode[TargetNodeID].Eps2;
        double twidth2 = SQ(2.0*twidth)+SQ(LETNode[TargetNodeID].EpsMax)-SQ(LETNode[TargetNodeID].EpsMin); //
        //double twidth2 = SQ(2.0*twidth); //
#endif // USE_SYMMETRIZED_SOFTENING

#ifdef USE_SYMMETRIZED_SOFTENING

        if( (SQ(2*twidth)<theta2*sqDist) && (SQ(LETNode[TargetNodeID].EpsMax-LETNode[TargetNodeID].EpsMin)<GravityRoot.OpeningAngle*sqDist) && (sqDist3D > TINY)) { // For self-interactions.
#else  // USE_SYMMETRIZED_SOFTENING
        if(SQ(2.0*twidth)<theta2*sqDist){
#endif // USE_SYMMETRIZED_SOFTENING
            if(NumberofField+1 >= NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+1));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif // USE_SYMMETRIZED_SOFTENING
            }
#ifdef USE_SHIFT_GRAVITYFIELD
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0]-LETNode[CurrentNodeID].COM[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1]-LETNode[CurrentNodeID].COM[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2]-LETNode[CurrentNodeID].COM[2];
#else 
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2];
#endif // USE_SHIFT_GRAVITYFIELD
            FieldMass[NumberofField]   = LETNode[TargetNodeID].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
            FieldEps2[NumberofField]   = LETNode[TargetNodeID].Eps2;
#endif // USE_SYMMETRIZED_SOFTENING
    
            NumberofField ++;
#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
            if(LETNode[TargetNodeID].NumberofLeaves == 1){
                NumberofPPInteractions ++;
            }else{
                NumberofPCInteractions ++;
            }
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
			TargetNodeID = LETNode[TargetNodeID].Next;
#ifdef USE_SYMMETRIZED_SOFTENING
		 } else if ((LETNode[TargetNodeID].Children == NONE) || (sqDist3D < TINY)) { // For self-interactions.
#else
		 } else if (LETNode[TargetNodeID].Children == NONE) {
#endif // USE_SYMMETRIZED_SOFTENING
            int NumberofLeaves = LETNode[TargetNodeID].NumberofLeaves;
            if(NumberofField+NumberofLeaves+1 > NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+NumberofLeaves+1));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif
            }
            int header = LETNode[TargetNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header + k;
#ifdef USE_SHIFT_GRAVITYFIELD
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0]-LETNode[CurrentNodeID].COM[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1]-LETNode[CurrentNodeID].COM[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2]-LETNode[CurrentNodeID].COM[2];
#else
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2];
#endif
                FieldMass[NumberofField]   = LETCache[leaf].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2[NumberofField]   = SQ(LETCache[leaf].Eps);
#endif
                NumberofField ++;
#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
                NumberofPPInteractions ++;
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR

            }
			TargetNodeID = LETNode[TargetNodeID].Next;
		 } else {
			TargetNodeID = LETNode[TargetNodeID].Children;
		 }
	}

#if 0
    if(NumberofField%2 == 1){
        if(NumberofField == NumberofFieldAllocated){
            NumberofFieldAllocated ++;
            FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
            FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
            FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif 
        }
        FieldPos[NumberofField][0] = 0.e0;
        FieldPos[NumberofField][1] = 0.e0;
        FieldPos[NumberofField][2] = 0.e0;
        FieldMass[NumberofField] = 0.e0;
        FieldEps2[NumberofField] = 1.e0;
        NumberofField++;
    }

#endif


	return (NumberofField);
}
#endif

static void CalculateForceParallelTreeGRAPEEngine(const int CurrentNodeID){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    int LeavesLog[Npipes];

    int NumberofActiveLeaves = LETNode[CurrentNodeID].NumberofActiveLeaves;
    int header = LETNode[CurrentNodeID].Leaves;
    int CurrentLeafID = 0;
    for(int i=0;i<NumberofActiveLeaves;i+=Npipes){
		int NActives = MIN(Npipes,NumberofActiveLeaves-i);
        int Count = 0;
        //int leaf = header+CurrentLeafID;
        while(Count!=NActives){
            int leaf = header + CurrentLeafID;
            if(LETCache[leaf].Active){
#ifdef USE_SHIFT_GRAVITYFIELD
#if MakeGlobalTree
                Dxi[Count][0] = LETCache[leaf].Pos[0]-LETNode[CurrentNodeID].COM[0];
                Dxi[Count][1] = LETCache[leaf].Pos[1]-LETNode[CurrentNodeID].COM[1];
                Dxi[Count][2] = LETCache[leaf].Pos[2]-LETNode[CurrentNodeID].COM[2];
#else
                Dxi[Count][0] = LETCache[leaf].Pos[0]-GravityNode[CurrentNodeID].COM[0];
                Dxi[Count][1] = LETCache[leaf].Pos[1]-GravityNode[CurrentNodeID].COM[1];
                Dxi[Count][2] = LETCache[leaf].Pos[2]-GravityNode[CurrentNodeID].COM[2];
#endif
#else
                Dxi[Count][0] = LETCache[leaf].Pos[0];
                Dxi[Count][1] = LETCache[leaf].Pos[1];
                Dxi[Count][2] = LETCache[leaf].Pos[2];
#endif
#if (defined(HAVE_AVX_PHANTOM_GRAPE)||defined(HAVE_AVX_PHANTOM_GRAPE_API2)) //{
                Depsi[Count] = SQ(LETCache[leaf].Eps);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
                Depsi[Count] = LETCache[leaf].Eps;
#endif // HAVE_AVX_PHANTOM_GRAPE //}
                LeavesLog[Count] = LETCache[leaf].Leaf;
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

#ifdef USE_SYMMETRIZED_SOFTENING
#ifdef HAVE_AVX_PHANTOM_GRAPE //{
        g5_calculate_force_on_x0(Dxi,adummy,phidummy,Npipes,Depsi);
#elif  HAVE_AVX_PHANTOM_GRAPE_API2 //}//{
        g5_calculate_force_on_xe(Dxi,Depsi,adummy,phidummy,Npipes);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
        g5_set_xi(Npipes,Dxi);
        g5_set_eps(Npipes,Depsi);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else 
#if defined(HAVE_GRAPE7)
        g5_set_eps2_to_all(SQ(Depsi[0]));
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE6A)
        g5_set_xepsi(Npipes,Dxi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE5)
        g5_set_ip(Npipes,Dxi,Depsi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_PHANTOM_GRAPE)
        g5_set_eps_to_all(Depsi[0]);
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#else
        g5_set_ip_emu(Npipes,Dxi,Depsi,Depsi);
        g5_run_emu();
        g5_get_force_emu(Npipes,adummy,phidummy);
#endif
#endif

        // back acc and pot for the temporal data array.
        for(int k=0;k<NActives;k++){
			GravityAccPotCache[LeavesLog[k]].Acc[0] += adummy[k][0];
			GravityAccPotCache[LeavesLog[k]].Acc[1] += adummy[k][1];
			GravityAccPotCache[LeavesLog[k]].Acc[2] += adummy[k][2];
			GravityAccPotCache[LeavesLog[k]].Pot += phidummy[k];
        }
    }

	return;
}

///////////////////////////////////////////////////////////////////////////////////////////

static void LETPreprocessing(const int NumberofLeavesThisStep);
static void MakeLETRoot(const int NumberofLeavesThisStep);
static int NextLETNode(const int NodeID);
static void BuildLET(void);
static void CopyFromImportLeavesToLETCache(void);
static void BuildLETTraversalLink(void);
static void LETNodeDataImplant(void);
static void LETNodeDataImplantNew(const int CurrentNodeID);
static int TraversalID[TreeMaxNodeLevel];

static void PlantLET(const int NumberofLeavesThisStep){

    LETPreprocessing(NumberofLeavesThisStep);

	MakeLETRoot(NumberofLeavesThisStep);
    BuildLET();
    CopyFromImportLeavesToLETCache();

    for(int i=1;i<LETRoot.NumberofNodes;i++){
        LETNode[i].Next = NextLETNode(i);
    }
    double diag = DISTANCE(LETNode[0].Pos,LETNode[LETNode[0].Children].Pos);
    LETNodeDataImplantNew(0);

	return;
}

static int NumberofLETLeavesAllocated = 0;
static int *LETSublist; 

static void LETPreprocessing(const int NumberofLeavesThisStep){

    if(NumberofLeavesThisStep > LETRoot.NumberofAllocatedLeaves){
        LETRoot.NumberofAllocatedLeaves = (int)(ForAngelsShare*NumberofLeavesThisStep);
        free(LETRoot.Leaves);
        free(LETCache);
        LETRoot.Leaves = malloc(sizeof(int)*LETRoot.NumberofAllocatedLeaves);
        LETCache = malloc(sizeof(struct StructGravityCache)*LETRoot.NumberofAllocatedLeaves);
    }

    if(NumberofLeavesThisStep > NumberofLETLeavesAllocated){
        if(NumberofLETLeavesAllocated > 0)
            free(LETSublist);
        NumberofLETLeavesAllocated = (int)(ForAngelsShare*NumberofLeavesThisStep);
	    LETSublist = malloc(sizeof(int)*NumberofLETLeavesAllocated);
    }

    return ;
}

static void MakeLETRoot(const int NumberofLeavesThisStep){

    int RootNodeID = 0;

	for(int k=0;k<3;k++){
        LETNode[RootNodeID].Pos[k] = GravityNode[RootNodeID].Pos[k];
        LETRoot.PosMax[k] = GravityRoot.PosMax[k];
        LETRoot.PosMin[k] = GravityRoot.PosMin[k];
    }
    LETRoot.Width = GravityRoot.Width;

    LETNode[RootNodeID].Next = NONE;
    LETNode[RootNodeID].Parent = NONE;
    LETNode[RootNodeID].Sister = NONE;
    LETNode[RootNodeID].Children = NONE;

    LETNode[RootNodeID].NumberofLeaves = NumberofLeavesThisStep;

    LETNode[RootNodeID].Level = 0;
    LETNode[RootNodeID].Leaves = 0;

    for(int i=0;i<NumberofLeavesThisStep;i++)
        LETRoot.Leaves[i] = i;

    LETRoot.NumberofLeaves = NumberofLeavesThisStep;

	return ;
}

static void CopyFromImportLeavesToLETCache(void){

    int NumberofLeaves = LETNode[0].NumberofLeaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf =  LETRoot.Leaves[k];
        LETCache[k].Pos[0] = LeavesImport[leaf].Pos[0];
        LETCache[k].Pos[1] = LeavesImport[leaf].Pos[1];
        LETCache[k].Pos[2] = LeavesImport[leaf].Pos[2];
        LETCache[k].Mass   = LeavesImport[leaf].Mass;

#ifdef USE_SYMMETRIZED_SOFTENING
        LETCache[k].Eps    = LeavesImport[leaf].Eps;
#endif
        if(leaf<Pall.Ntotal){
            LETCache[k].Active = GravityCache[leaf].Active;
            LETCache[k].Leaf   = GravityCache[leaf].Leaf;
#ifndef USE_SYMMETRIZED_SOFTENING
            LETCache[k].Eps    = GravityCache[leaf].Eps;
#endif
        } else {
            LETCache[k].Active = false;
            LETCache[k].Leaf   = NONE;
#ifdef USE_SYMMETRIZED_SOFTENING
            //LETCache[k].Eps = LeavesImport[leaf].Eps;
#endif // USE_SYMMETRIZED_SOFTENING
        }
    }
    return;
}

static int NextLETNode(const int NodeID){

    int CurrentNodeID = NodeID;

    if(LETNode[CurrentNodeID].Sister != NONE){
        CurrentNodeID = LETNode[CurrentNodeID].Sister;
    } else {
        int NextNodeID = CurrentNodeID;
        while(1){
            if(LETNode[LETNode[NextNodeID].Parent].Sister != NONE){
                CurrentNodeID = LETNode[LETNode[NextNodeID].Parent].Sister;
                break;
            } else if(LETNode[NextNodeID].Parent == 0){
                CurrentNodeID = 0;
                break;
            }
            NextNodeID = LETNode[NextNodeID].Parent;
        }
    }
    return CurrentNodeID;
}


static inline bool LETNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber) __attribute__((always_inline));
static inline bool LETNodeSeparationCriterion(const int CurrentNodeID, const int CriticalNumber){

	if( (LETNode[CurrentNodeID].NumberofLeaves <= CriticalNumber) || LETNode[CurrentNodeID].Level+1>=TreeMaxNodeLevel){
        return true;
    } else {
        return false;
    }
}

static void BuildLET(void){

    int NumberofNodes = 0; 

    int NumberofNodeCreationLimit = LETRoot.NumberofNodeCreationLimit;

    int CurrentMaxLevel = 0;
    int RootNodeID = 0; 
    int CurrentNodeID = RootNodeID;
    int ChildNodeID,BackwardNodeID,NextNodeID;
	while(1){

		if(LETNodeSeparationCriterion(CurrentNodeID,NumberofNodeCreationLimit) && (CurrentNodeID != RootNodeID)){
			if(LETNode[CurrentNodeID].Sister != NONE){
				CurrentNodeID = LETNode[CurrentNodeID].Sister;
			}else{
				NextNodeID = CurrentNodeID;
				while(1){
                    if(LETNode[LETNode[NextNodeID].Parent].Sister != NONE){
                        CurrentNodeID = LETNode[LETNode[NextNodeID].Parent].Sister;
						break;
					} else if(LETNode[NextNodeID].Parent == RootNodeID){
                        LETRoot.CurrentMaxLevel = CurrentMaxLevel;
                        LETRoot.NumberofNodes = NumberofNodes + 1;
						return;
					}
                    NextNodeID = LETNode[NextNodeID].Parent;
				}
			}
			continue;
		}

        int subhead[TreeNsub],subcurrent[TreeNsub],subnumber[TreeNsub]; 
		for(int k=0;k<TreeNsub;k++){
			subnumber[k] = 0;
			subhead[k] = subcurrent[k] = NONE;
		}

        int NumberofLeaves = LETNode[CurrentNodeID].NumberofLeaves;
        int header = LETNode[CurrentNodeID].Leaves;
		for(int i=0;i<NumberofLeaves;i++){
            int leaf = LETRoot.Leaves[header+i];

            double Pos[3] = {LeavesImport[leaf].Pos[0],LeavesImport[leaf].Pos[1],LeavesImport[leaf].Pos[2]};

            int subindex0 = ((LETNode[CurrentNodeID].Pos[0] <= Pos[0])?1:0);
            int subindex1 = ((LETNode[CurrentNodeID].Pos[1] <= Pos[1])?1:0);
            int subindex2 = ((LETNode[CurrentNodeID].Pos[2] <= Pos[2])?1:0);
            int subindex = subindex0 | subindex1 << 1 | subindex2 << 2;

			if(subnumber[subindex] > 0){
				LETSublist[subcurrent[subindex]] = leaf;
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
                if(NumberofNodes >= LETRoot.NumberofAllocatedNodes){
                    int NumberofAllocatedNodes = (int)(MAX(ForAngelsShare*NumberofNodes,NAdditionUnit));
                    LETNode = realloc(LETNode,sizeof(struct StructGravityNode)*NumberofAllocatedNodes);
                    LETRoot.NumberofAllocatedNodes = NumberofAllocatedNodes;
                }

                LETNode[ChildNodeID].Next = NONE;
                LETNode[ChildNodeID].Parent = NONE;
                LETNode[ChildNodeID].Sister = NONE;
                LETNode[ChildNodeID].Children = NONE;

                LETNode[ChildNodeID].Parent = CurrentNodeID;

                if(BackwardNodeID == CurrentNodeID){
                    LETNode[CurrentNodeID].Children = ChildNodeID;
					NextNodeID = ChildNodeID;
                    LETNode[ChildNodeID].Leaves = LETNode[CurrentNodeID].Leaves;
                    CurrentMaxLevel = MAX(CurrentMaxLevel,LETNode[CurrentNodeID].Level+1);
                } else {
                    LETNode[BackwardNodeID].Sister = ChildNodeID;
                    LETNode[ChildNodeID].Leaves = 
                        LETNode[BackwardNodeID].Leaves + LETNode[BackwardNodeID].NumberofLeaves;
                }

                int cheader = LETNode[ChildNodeID].Leaves;
                LETRoot.Leaves[cheader] = subhead[i];
                for(int k=1;k<subnumber[i];k++){
                    LETRoot.Leaves[cheader+k] = 
                        LETSublist[LETRoot.Leaves[cheader+k-1]];
                }
                LETNode[ChildNodeID].NumberofLeaves = subnumber[i];
                LETNode[ChildNodeID].Level = LETNode[CurrentNodeID].Level+1;

				for(int k=0;k<3;k++)
                    LETNode[ChildNodeID].Pos[k] = LETNode[CurrentNodeID].Pos[k] +
						+ BitSign((i>>k)&1)*0.25e0*LETRoot.Width*LETRoot.WidthFactor[LETNode[CurrentNodeID].Level];
            }
		}
        CurrentNodeID = NextNodeID;
	}
}


static inline void CalcDistanceMaxMassCOM(const int CurrentNodeID) __attribute__((always_inline));
static inline void CalcDistanceMaxMassCOM(const int CurrentNodeID){

    double CenterofNode[3] = {LETNode[CurrentNodeID].Pos[0],LETNode[CurrentNodeID].Pos[1],LETNode[CurrentNodeID].Pos[2]};
    double distancemax = 0.e0;
    double mass = 0.e0;
    double COM[3] = {0.e0,0.e0,0.e0};
    int NumberofActiveLeaves = 0;

    int NumberofLeaves = LETNode[CurrentNodeID].NumberofLeaves;
    int header = LETNode[CurrentNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
        distancemax = fmax(DISTANCE(CenterofNode,LETCache[leaf].Pos),distancemax);
        mass   += LETCache[leaf].Mass;
        COM[0] += LETCache[leaf].Mass*LETCache[leaf].Pos[0];
        COM[1] += LETCache[leaf].Mass*LETCache[leaf].Pos[1];
        COM[2] += LETCache[leaf].Mass*LETCache[leaf].Pos[2];
        NumberofActiveLeaves += LETCache[leaf].Active;
    }
    double imass = 1.e0/mass;
    LETNode[CurrentNodeID].COM[0] = COM[0]*imass;
    LETNode[CurrentNodeID].COM[1] = COM[1]*imass;
    LETNode[CurrentNodeID].COM[2] = COM[2]*imass;

    LETNode[CurrentNodeID].Mass = mass;
    LETNode[CurrentNodeID].DistanceMax = distancemax;
    LETNode[CurrentNodeID].NumberofActiveLeaves = NumberofActiveLeaves;

    return;
}


static double Diag;
static void LETNodeDataImplantNew(const int CurrentNodeID){
    double Width = Diag*GravityRoot.WidthFactor[LETNode[CurrentNodeID].Level];

    int NActives = 0;
    double DistanceMax = 0.e0;
    double Mass = 0.e0;
    double COM[3] = {0.e0,0.e0,0.e0};
    double Eps2 = 0.e0;
    double EpsMax,EpsMin;
    if(LETNode[CurrentNodeID].Children == NONE){
        int Number_of_leaf = LETNode[CurrentNodeID].NumberofLeaves;
        int header = LETNode[CurrentNodeID].Leaves;
#ifdef USE_SYMMETRIZED_SOFTENING
        EpsMax = EpsMin = LETCache[header].Eps;
#endif // USE_SYMMETRIZED_SOFTENING
        for(int k=0;k<Number_of_leaf;k++){
            int leaf = header+k;
            double Distance = DISTANCE(LETNode[CurrentNodeID].Pos,LETCache[leaf].Pos);
            DistanceMax = fmax(Distance,DistanceMax);
            Mass += LETCache[leaf].Mass;
            COM[0] += LETCache[leaf].Mass*LETCache[leaf].Pos[0];
            COM[1] += LETCache[leaf].Mass*LETCache[leaf].Pos[1];
            COM[2] += LETCache[leaf].Mass*LETCache[leaf].Pos[2];
            NActives += LETCache[leaf].Active;
#ifdef USE_SYMMETRIZED_SOFTENING
            Eps2 += LETCache[leaf].Mass*SQ(LETCache[leaf].Eps);
            EpsMin = fmin(EpsMin,LETCache[leaf].Eps);
            EpsMax = fmax(EpsMax,LETCache[leaf].Eps);
#endif // USE_SYMMETRIZED_SOFTENING
        }
        Width = 0.e0;
    } else {
        bool first = true;
        int ChildNodeID = LETNode[CurrentNodeID].Children;
        while(ChildNodeID != NONE){
            LETNodeDataImplantNew(ChildNodeID);

            DistanceMax = fmax(DistanceMax,LETNode[ChildNodeID].DistanceMax);
            Mass += LETNode[ChildNodeID].Mass;
            COM[0] += LETNode[ChildNodeID].Mass*LETNode[ChildNodeID].COM[0];
            COM[1] += LETNode[ChildNodeID].Mass*LETNode[ChildNodeID].COM[1];
            COM[2] += LETNode[ChildNodeID].Mass*LETNode[ChildNodeID].COM[2];
            NActives += LETNode[ChildNodeID].NumberofActiveLeaves;
#ifdef USE_SYMMETRIZED_SOFTENING
            Eps2 += LETNode[ChildNodeID].Mass*LETNode[ChildNodeID].Eps2;
            if(first == true){
                EpsMin = LETNode[ChildNodeID].EpsMin;
                EpsMax = LETNode[ChildNodeID].EpsMax;
                first = false;
            } else {
                EpsMin = fmin(EpsMin,LETNode[ChildNodeID].EpsMin);
                EpsMax = fmax(EpsMax,LETNode[ChildNodeID].EpsMax);
            }
#endif // USE_SYMMETRIZED_SOFTENING
            ChildNodeID = LETNode[ChildNodeID].Sister;
        }
    }
    double InvMass = 1.e0/Mass;
    LETNode[CurrentNodeID].COM[0] = COM[0]*InvMass;
    LETNode[CurrentNodeID].COM[1] = COM[1]*InvMass;
    LETNode[CurrentNodeID].COM[2] = COM[2]*InvMass;

    LETNode[CurrentNodeID].DistanceMax = DistanceMax + Width;
    LETNode[CurrentNodeID].Mass = Mass;
    LETNode[CurrentNodeID].NumberofActiveLeaves = NActives;
#ifdef USE_SYMMETRIZED_SOFTENING
    LETNode[CurrentNodeID].Eps2 = Eps2*InvMass;
    LETNode[CurrentNodeID].EpsMax = EpsMax;
    LETNode[CurrentNodeID].EpsMin = EpsMin;
#endif // USE_SYMMETRIZED_SOFTENING
    return ;

}

#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
static void WalkLocalTreeAndGetAccPotForInsert(const int NImportAll);
static void CalculateForceParallelTreeGRAPEEngineForInsert(const int CurrentNodeID);
static int GetFieldFromLETForInsert(const int CurrentNodeID, const int NField);

void ForceParallelTreeGRAPEInsert(const int NImport, double Pos[restrict][3], double Mass[restrict], double Eps[restrict]){

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    HoldGRAPE();

    SetGRAPEScaleForParallelTreeGRAPE();
    for(int i=0;i<Pall.Ntotal;i++){
        GravityAccPotCache[i].Acc[0] = 
        GravityAccPotCache[i].Acc[1] = 
        GravityAccPotCache[i].Acc[2] = 
        GravityAccPotCache[i].Pot = 0.e0;
        GravityAccPotCache[i].InteractionList = 0;
    }

    CheckSizeofBufferExportRecv(NImport,sizeof(struct StructLeavesExportImport));
    LeavesImport = BufferExportRecv;

    for(int i=0;i<NImport;i++){ // move
        LeavesImport[i].Pos[0] = Pos[i][0];
        LeavesImport[i].Pos[1] = Pos[i][1];
        LeavesImport[i].Pos[2] = Pos[i][2];
        LeavesImport[i].Mass   = Mass[i];
#ifdef USE_SYMMETRIZED_SOFTENING
        LeavesImport[i].Eps    = Eps[i];
#endif
    }

    if(NImport>0)
        PlantLET(NImport);

    WalkLocalTreeAndGetAccPotForInsert(NImport);

    ReleaseGRAPE();
#ifdef USE_SYMMETRIZED_SOFTENING
    double SymmetrizedFactor = 1.0/sqrt(2.0);
#endif
    double AdaptiveSofteningFactor = Pall.AdaptiveSofteningFactor;
    for(int i=0;i<Pall.Ntotal;i++){
        if(GravityCache[i].Active){
            int leaf = GravityCache[i].Leaf;
            Pbody[leaf]->Acc[0] = Pall.GravConst*GravityAccPotCache[leaf].Acc[0];
            Pbody[leaf]->Acc[1] = Pall.GravConst*GravityAccPotCache[leaf].Acc[1];
            Pbody[leaf]->Acc[2] = Pall.GravConst*GravityAccPotCache[leaf].Acc[2];
#ifdef HAVE_GRAPE7
            Pbody[leaf]->Pot = 0.e0;
#else //HAVE_GRAPE7
#if 0
#ifdef USE_SYMMETRIZED_SOFTENING
            Pbody[leaf]->Pot = GravityAccPotCache[leaf].Pot 
                - SymmetrizedFactor*(Pbody[leaf]->Mass/(AdaptiveSofteningFactor*Pbody[leaf]->Eps));
#else //USE_SYMMETRIZED_SOFTENING
            Pbody[leaf]->Pot = GravityAccPotCache[leaf].Pot 
                - (Pbody[leaf]->Mass/(AdaptiveSofteningFactor*Pbody[leaf]->Eps));
#endif///USE_SYMMETRIZED_SOFTENING
#endif
            Pbody[leaf]->Pot = -0.5*Pall.GravConst*Pbody[leaf]->Mass*GravityAccPotCache[leaf].Pot;
#endif //HAVE_GRAPE7
            Pbody[leaf]->InteractionList = GravityAccPotCache[leaf].InteractionList;
        }
    }

    ForceEndProcedure();

    return;
}

static void WalkLocalTreeAndGetAccPotForInsert(const int NImportAll){

    if(NumberofFieldAllocated == 0){
        NumberofFieldAllocated = FirstAllocationSize;
        FieldPos = malloc(sizeof(double)*3*NumberofFieldAllocated);
        FieldMass = malloc(sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
        FieldEps2 = malloc(sizeof(double)*NumberofFieldAllocated);
#endif
    }

    int RootNodeID = 0;
    int CurrentNodeID = GravityNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
		if(GravityNode[CurrentNodeID].NumberofActiveLeaves == 0){ 
            CurrentNodeID = GravityNode[CurrentNodeID].Next;
		}else if( (GravityNode[CurrentNodeID].NumberofLeaves < GravityRoot.NumberofLeavesInGroup)
                || (GravityNode[CurrentNodeID].Children == NONE)){

            int NumberofField = GetFieldFromLETForInsert(CurrentNodeID,0);
            //dprintlmpi(NumberofField);

            int NumberofLeaves = GravityNode[CurrentNodeID].NumberofLeaves;
            int header = GravityNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(GravityCache[leaf].Active)
                    GravityAccPotCache[GravityCache[leaf].Leaf].InteractionList = NumberofField;
            }

            for(int j=0;j<NumberofField;j+=JMEMSIZE){ 
                int Active_j = MIN(JMEMSIZE,NumberofField-j);
#ifdef USE_SYMMETRIZED_SOFTENING
                g5_set_n(Active_j);
#ifdef HAVE_AVX_PHANTOM_GRAPE //{
                g5_set_xmj0(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#elif  HAVE_AVX_PHANTOM_GRAPE_API2 //{
                g5_set_xmj(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
                g5_set_xmeps2j(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else 
#if defined(HAVE_GRAPE7)
                g5_set_n(Active_j);
                g5_set_jp(0,Active_j,FieldMass+j,FieldPos+j);
#elif (defined(HAVE_GRAPE6A) \
    || defined(HAVE_GRAPE5) \
    || defined(HAVE_PHANTOM_GRAPE) \
    || defined(HAVE_AVX_PHANTOM_GRAPE) \
    || defined(HAVE_AVX_PHANTOM_GRAPE_API2)) 
                g5_set_n(Active_j);
                g5_set_xmj(0,Active_j,FieldPos+j,FieldMass+j);
#else
                g5_set_n_emu(Active_j);
                g5_set_xmj_emu(0,Active_j,FieldPos+j,FieldMass+j);
#endif
#endif
                CalculateForceParallelTreeGRAPEEngineForInsert(CurrentNodeID);
            }
            CurrentNodeID = GravityNode[CurrentNodeID].Next;
		} else {
		    CurrentNodeID = GravityNode[CurrentNodeID].Children;
		}
	}

    return;
}

static int GetFieldFromLETForInsert(const int CurrentNodeID, const int NField){


    int NumberofField = NField;
    double cwidth = 0.5*GravityRoot.Width*GravityRoot.WidthFactor[GravityNode[CurrentNodeID].Level];
    double theta2 = SQ(GravityRoot.OpeningAngle);

    int RootNodeID = 0;
	int TargetNodeID = LETNode[RootNodeID].Children;
	while(TargetNodeID != RootNodeID){ 
        double twidth = 0.5*LETRoot.Width*LETRoot.WidthFactor[LETNode[TargetNodeID].Level];
        double sqDist = 0.e0;

        /*
        double Dist[3] = {fabs(GravityNode[CurrentNodeID].Pos[0]-LETNode[TargetNodeID].Pos[0]),
                          fabs(GravityNode[CurrentNodeID].Pos[1]-LETNode[TargetNodeID].Pos[1]),
                          fabs(GravityNode[CurrentNodeID].Pos[2]-LETNode[TargetNodeID].Pos[2])};
        double hwidth = cwidth + twidth;
        sqDist += SQ(fmax(Dist[0]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[1]-hwidth,0.0));
        sqDist += SQ(fmax(Dist[2]-hwidth,0.0));
        */
        double Dist[3] = {fabs(GravityNode[CurrentNodeID].Pos[0]-LETNode[TargetNodeID].COM[0]),
                          fabs(GravityNode[CurrentNodeID].Pos[1]-LETNode[TargetNodeID].COM[1]),
                          fabs(GravityNode[CurrentNodeID].Pos[2]-LETNode[TargetNodeID].COM[2])};
        sqDist += SQ(fmax(Dist[0],0.0));
        sqDist += SQ(fmax(Dist[1],0.0));
        sqDist += SQ(fmax(Dist[2],0.0));
#ifdef USE_SYMMETRIZED_SOFTENING
        double sqDist3D = sqDist;
        sqDist += SQ(GravityNode[CurrentNodeID].EpsMin);
        sqDist += LETNode[TargetNodeID].Eps2;
        double twidth2 = SQ(2*twidth)+SQ(LETNode[TargetNodeID].EpsMax)-SQ(LETNode[TargetNodeID].EpsMin); //
#endif // USE_SYMMETRIZED_SOFTENING

#ifdef USE_SYMMETRIZED_SOFTENING
        if( (twidth2<theta2*sqDist) && (sqDist3D > TINY)) { // For self-interactions.
        //if( (twidth2<theta2*sqDist) && (twidth2<theta2*sqDist3D)) { // For self-interactions.
        //if( (twidth2<theta2*sqDist)) { // For self-interactions.
#else
        if(SQ(2.0*twidth)<theta2*sqDist){
#endif // USE_SYMMETRIZED_SOFTENING
            if(NumberofField+1 > NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+1));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif // USE_SYMMETRIZED_SOFTENING
            }
#ifdef USE_SHIFT_GRAVITYFIELD
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0]-GravityNode[CurrentNodeID].COM[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1]-GravityNode[CurrentNodeID].COM[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2]-GravityNode[CurrentNodeID].COM[2];
#else 
            FieldPos[NumberofField][0] = LETNode[TargetNodeID].COM[0];
            FieldPos[NumberofField][1] = LETNode[TargetNodeID].COM[1];
            FieldPos[NumberofField][2] = LETNode[TargetNodeID].COM[2];
#endif // USE_SHIFT_GRAVITYFIELD
            FieldMass[NumberofField]   = LETNode[TargetNodeID].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
            FieldEps2[NumberofField]   = LETNode[TargetNodeID].Eps2;
#endif // USE_SYMMETRIZED_SOFTENING
    
            NumberofField ++;
            if(LETNode[TargetNodeID].NumberofLeaves == 1){
                NumberofPPInteractions ++;
            }else{
                NumberofPCInteractions ++;
            }
			TargetNodeID = LETNode[TargetNodeID].Next;

#ifdef USE_SYMMETRIZED_SOFTENING
		 //} else if ((LETNode[TargetNodeID].Children == NONE) || (sqDist3D < TINY)) { // For self-interactions.
		 } else if (LETNode[TargetNodeID].Children == NONE) {
#else
		 } else if (LETNode[TargetNodeID].Children == NONE) {
#endif // USE_SYMMETRIZED_SOFTENING
            int NumberofLeaves = LETNode[TargetNodeID].NumberofLeaves;
            if(NumberofField+NumberofLeaves > NumberofFieldAllocated){
                NumberofFieldAllocated = (int)(ForAngelsShare*(NumberofField+NumberofLeaves));
                FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
                FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif
            }
            int header = LETNode[TargetNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header + k;
#ifdef USE_SHIFT_GRAVITYFIELD
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0]-GravityNode[CurrentNodeID].COM[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1]-GravityNode[CurrentNodeID].COM[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2]-GravityNode[CurrentNodeID].COM[2];
#else
                FieldPos[NumberofField][0] = LETCache[leaf].Pos[0];
                FieldPos[NumberofField][1] = LETCache[leaf].Pos[1];
                FieldPos[NumberofField][2] = LETCache[leaf].Pos[2];
#endif
                FieldMass[NumberofField]   = LETCache[leaf].Mass;
#ifdef USE_SYMMETRIZED_SOFTENING
                FieldEps2[NumberofField]   = SQ(LETCache[leaf].Eps);
#endif
                NumberofPPInteractions ++;
                NumberofField ++;
            }
			TargetNodeID = LETNode[TargetNodeID].Next;
		 } else {
			TargetNodeID = LETNode[TargetNodeID].Children;
		 }
	}

	return (NumberofField);
}

static void CalculateForceParallelTreeGRAPEEngineForInsert(const int CurrentNodeID){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    int LeavesLog[Npipes];

    int NumberofActiveLeaves = GravityNode[CurrentNodeID].NumberofActiveLeaves;
    int header = GravityNode[CurrentNodeID].Leaves;
    int CurrentLeafID = 0;
    for(int i=0;i<NumberofActiveLeaves;i+=Npipes){
		int NActives = MIN(Npipes,NumberofActiveLeaves-i);
        int Count = 0;
        while(Count!=NActives){
            int leaf = header + CurrentLeafID;
            if(LETCache[leaf].Active){
#ifdef USE_SHIFT_GRAVITYFIELD
                Dxi[Count][0] = GravityCache[leaf].Pos[0]-GravityNode[CurrentNodeID].COM[0];
                Dxi[Count][1] = GravityCache[leaf].Pos[1]-GravityNode[CurrentNodeID].COM[1];
                Dxi[Count][2] = GravityCache[leaf].Pos[2]-GravityNode[CurrentNodeID].COM[2];
#else
                Dxi[Count][0] = GravityCache[leaf].Pos[0];
                Dxi[Count][1] = GravityCache[leaf].Pos[1];
                Dxi[Count][2] = GravityCache[leaf].Pos[2];
#endif
#if (defined(HAVE_AVX_PHANTOM_GRAPE)||(defined(HAVE_AVX_PHANTOM_GRAPE_API2)) //{
                Depsi[Count] = SQ(GravityCache[leaf].Eps);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
                Depsi[Count] = GravityCache[leaf].Eps;
#endif // HAVE_AVX_PHANTOM_GRAPE //}
                LeavesLog[Count] = GravityCache[leaf].Leaf;
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

#ifdef USE_SYMMETRIZED_SOFTENING
#ifdef HAVE_AVX_PHANTOM_GRAPE //{
        g5_calculate_force_on_x0(Dxi,adummy,phidummy,Npipes,Depsi);
#elif  HAVE_AVX_PHANTOM_GRAPE_API2 //}//{
        g5_calculate_force_on_xe(Dxi,Depsi,adummy,phidummy,Npipes);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
        g5_set_xi(Npipes,Dxi);
        g5_set_eps(Npipes,Depsi);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else 
#if defined(HAVE_GRAPE7)
        g5_set_eps2_to_all(SQ(Depsi[0]));
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE6A)
        g5_set_xepsi(Npipes,Dxi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE5)
        g5_set_ip(Npipes,Dxi,Depsi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif (defined(HAVE_PHANTOM_GRAPE) \
     ||defined(HAVE_AVX_PHANTOM_GRAPE) \
     ||defined(HAVE_AVX_PHANTOM_GRAPE_API2))
        g5_set_eps_to_all(Depsi[0]);
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#else
        g5_set_ip_emu(Npipes,Dxi,Depsi,Depsi);
        g5_run_emu();
        g5_get_force_emu(Npipes,adummy,phidummy);
#endif
#endif

        // back acc and pot for the temporal data array.
        for(int k=0;k<NActives;k++){
			GravityAccPotCache[LeavesLog[k]].Acc[0] += adummy[k][0];
			GravityAccPotCache[LeavesLog[k]].Acc[1] += adummy[k][1];
			GravityAccPotCache[LeavesLog[k]].Acc[2] += adummy[k][2];
			GravityAccPotCache[LeavesLog[k]].Pot += phidummy[k];
        }
    }

	return;
}

static void CalculateForceEngineForPhantomTest(const int CurrentNodeID){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    int LeavesLog[Npipes];

    int NumberofActiveLeaves = GravityNode[CurrentNodeID].NumberofActiveLeaves;
    int header = GravityNode[CurrentNodeID].Leaves;
    int CurrentLeafID = 0;
    for(int i=0;i<NumberofActiveLeaves;i+=Npipes){
		int NActives = MIN(Npipes,NumberofActiveLeaves-i);
        int Count = 0;
        while(Count!=NActives){
            int leaf = header + CurrentLeafID;
#ifdef USE_SHIFT_GRAVITYFIELD
            Dxi[Count][0] = GravityCache[leaf].Pos[0]-GravityNode[CurrentNodeID].COM[0];
            Dxi[Count][1] = GravityCache[leaf].Pos[1]-GravityNode[CurrentNodeID].COM[1];
            Dxi[Count][2] = GravityCache[leaf].Pos[2]-GravityNode[CurrentNodeID].COM[2];
#else
            Dxi[Count][0] = GravityCache[leaf].Pos[0];
            Dxi[Count][1] = GravityCache[leaf].Pos[1];
            Dxi[Count][2] = GravityCache[leaf].Pos[2];
#endif
#if (defined(HAVE_AVX_PHANTOM_GRAPE)||defined(HAVE_AVX_PHANTOM_GRAPE_API2)) //{
            Depsi[Count] = SQ(GravityCache[leaf].Eps);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
            Depsi[Count] = GravityCache[leaf].Eps;
#endif // HAVE_AVX_PHANTOM_GRAPE //}
            LeavesLog[Count] = GravityCache[leaf].Leaf;
            Count ++;
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

#ifdef USE_SYMMETRIZED_SOFTENING
#ifdef HAVE_AVX_PHANTOM_GRAPE //{
        g5_calculate_force_on_x0(Dxi,adummy,phidummy,Npipes,Depsi);
#elif  HAVE_AVX_PHANTOM_GRAPE_API2 //}//{
        g5_calculate_force_on_xe(Dxi,Depsi,adummy,phidummy,Npipes);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
        g5_set_xi(Npipes,Dxi);
        g5_set_eps(Npipes,Depsi);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else 
#if defined(HAVE_GRAPE7)
        g5_set_eps2_to_all(SQ(Depsi[0]));
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE6A)
        g5_set_xepsi(Npipes,Dxi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif defined(HAVE_GRAPE5)
        g5_set_ip(Npipes,Dxi,Depsi,Depsi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#elif (defined(HAVE_PHANTOM_GRAPE)||defined(HAVE_AVX_PHANTOM_GRAPE))
        g5_set_eps_to_all(Depsi[0]);
        g5_set_xi(Npipes,Dxi);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
#else
        g5_set_ip_emu(Npipes,Dxi,Depsi,Depsi);
        g5_run_emu();
        g5_get_force_emu(Npipes,adummy,phidummy);
#endif
#endif

        // back acc and pot for the temporal data array.
        for(int k=0;k<NActives;k++){
			GravityAccPotCache[LeavesLog[k]].Acc[0] += adummy[k][0];
			GravityAccPotCache[LeavesLog[k]].Acc[1] += adummy[k][1];
			GravityAccPotCache[LeavesLog[k]].Acc[2] += adummy[k][2];
			GravityAccPotCache[LeavesLog[k]].Pot += phidummy[k];
        }
    }

	return;
}

void CalcGravityDirectSymmetrizedPotentialWithPhantom(void){

    if(NumberofFieldAllocated == 0){
        NumberofFieldAllocated = FirstAllocationSize;
        FieldPos = malloc(sizeof(double)*3*NumberofFieldAllocated);
        FieldMass = malloc(sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
        FieldEps2 = malloc(sizeof(double)*NumberofFieldAllocated);
#endif
    }

    if(Pall.Ntotal > NumberofFieldAllocated){
        NumberofFieldAllocated = (int)(ForAngelsShare*(Pall.Ntotal+1));
        FieldPos = realloc(FieldPos,sizeof(double)*3*NumberofFieldAllocated);
        FieldMass = realloc(FieldMass,sizeof(double)*NumberofFieldAllocated);
#ifdef USE_SYMMETRIZED_SOFTENING
        FieldEps2 = realloc(FieldEps2,sizeof(double)*NumberofFieldAllocated);
#endif // USE_SYMMETRIZED_SOFTENING
    }

    for(int i=0;i<Pall.Ntotal;i++){
        FieldPos[i][0] = Pbody[i]->Pos[0];
        FieldPos[i][1] = Pbody[i]->Pos[1];
        FieldPos[i][2] = Pbody[i]->Pos[2];
        FieldMass[i]   = Pbody[i]->Mass;
        FieldEps2[i]   = SQ(Pbody[i]->Eps);
    }


    for(int j=0;j<Pall.Ntotal;j+=JMEMSIZE){ 
        int Active_j = MIN(JMEMSIZE,Pall.Ntotal-j);
        dprintlmpi(Active_j);
#ifdef USE_SYMMETRIZED_SOFTENING
        g5_set_n(Active_j);
#ifdef HAVE_AVX_PHANTOM_GRAPE //{
        g5_set_xmj0(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#else // HAVE_AVX_PHANTOM_GRAPE //}//{
        g5_set_xmeps2j(0,Active_j,FieldPos+j,FieldMass+j,FieldEps2+j);
#endif // HAVE_AVX_PHANTOM_GRAPE //}
#else 
#if defined(HAVE_GRAPE7)
        g5_set_n(Active_j);
        g5_set_jp(0,Active_j,FieldMass+j,FieldPos+j);
#elif (defined(HAVE_GRAPE6A) || defined(HAVE_GRAPE5) || defined(HAVE_PHANTOM_GRAPE) || defined(HAVE_AVX_PHANTOM_GRAPE))
        g5_set_n(Active_j);
        g5_set_xmj(0,Active_j,FieldPos+j,FieldMass+j);
#else
        g5_set_n_emu(Active_j);
        g5_set_xmj_emu(0,Active_j,FieldPos+j,FieldMass+j);
#endif
#endif

        CalculateForceEngineForPhantomTest(0);

    }
}
#endif

