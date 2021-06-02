#include "config.h"
#include "ForceMisc.h"

struct StructForceExport{
    double Pos[3];
    double Eps;
};

struct StructForceImport{
    double Acc[3];
    double Pot;
};


static void CalcGravityPlummer(double x[3], double e, double Acc[3], double *Pot);
static struct StructForceImport CalcGravityTree(struct StructForceExport ForceExport);

static void WriteAccPot(void){

    char name[MaxCharactersInLine];
    FILE *fp;

    sprintf(name,"TreeAccPot.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,name,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
            Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],Pbody[i]->Pot);
    }
    fclose(fp);

    return;
}

void ParallelTreeNew(void){

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    int Actives = 0;
    for(int i=0;i<Pall.Ntotal;i++)
        if(Pbody[i]->Active)
            Actives ++;
    int GlobalActives;
    MPI_Allreduce(&Actives,&GlobalActives,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    int NImportThisTime[NProcs];
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(&Actives,1,MPI_INT,CommunicationTable[i].SendRank,TAG_FORCE_TREE_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_FORCE_TREE_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }
    assert(NImportAll+Actives == GlobalActives);

    struct StructForceExport *ForceExportSendRecv;
    struct StructForceImport *ForceImportRecv;
    ForceImportRecv = malloc(sizeof(struct StructForceImport)*(Actives+1));

    CheckSizeofBufferExportSendIndex(GlobalActives,sizeof(struct StructForceExport),0);
    CheckSizeofBufferExportRecv(Actives,sizeof(struct StructForceImport));
    ForceExportSendRecv = BufferExportSend[0];
    ForceImportRecv = BufferExportRecv;

    Actives = 0;
    int RootNodeID = 0;
    int NumberofLeaves = GravityNode[RootNodeID].NumberofLeaves;
    int header = GravityNode[RootNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header + k;
        if(GravityCache[leaf].Active){
            int index = GravityCache[leaf].Leaf;

            ForceExportSendRecv[Actives].Pos[0] = Pbody[index]->PosP[0];
            ForceExportSendRecv[Actives].Pos[1] = Pbody[index]->PosP[1];
            ForceExportSendRecv[Actives].Pos[2] = Pbody[index]->PosP[2];
            ForceExportSendRecv[Actives].Eps = Pbody[index]->Eps;
            Actives ++;
        }
    }


    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    double TimeComm = GetElapsedTime();
    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(ForceExportSendRecv,
                sizeof(struct StructForceExport)*Actives,
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_FORCE_TREE_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(ForceExportSendRecv+NImportAll+Actives,
                sizeof(struct StructForceExport)*NImportThisTime[i],
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FORCE_TREE_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImportAll += NImportThisTime[i];
    }
    assert(NImportAll+Actives == GlobalActives);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.GravityCommThisStep += GetElapsedTime()-TimeComm;

    TimeComm = GetElapsedTime();
    for(int i=0;i<GlobalActives;i++){
        struct StructForceImport TempForceImport = CalcGravityTree(ForceExportSendRecv[i]);
        ForceExportSendRecv[i].Pos[0] = TempForceImport.Acc[0];
        ForceExportSendRecv[i].Pos[1] = TempForceImport.Acc[1];
        ForceExportSendRecv[i].Pos[2] = TempForceImport.Acc[2];
        ForceExportSendRecv[i].Eps = TempForceImport.Pot;
    }
    TimingResults.GravityThisStep += GetElapsedTime()-TimeComm;

    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(ForceExportSendRecv+NImportAll+Actives,sizeof(struct StructForceExport)*NImportThisTime[i],
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_FORCE_TREE_IMPORT,
            ForceImportRecv,sizeof(struct StructForceImport)*Actives,
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FORCE_TREE_IMPORT,
                    MPI_COMM_WORLD,&mpi_status);
        for(int k=0;k<Actives;k++){
            ForceExportSendRecv[k].Pos[0] += ForceImportRecv[k].Acc[0];
            ForceExportSendRecv[k].Pos[1] += ForceImportRecv[k].Acc[1];
            ForceExportSendRecv[k].Pos[2] += ForceImportRecv[k].Acc[2];
            ForceExportSendRecv[k].Eps += ForceImportRecv[k].Pot;
        }
        NImportAll += NImportThisTime[i];
    }

    Actives = 0;
    RootNodeID = 0;
    NumberofLeaves = GravityNode[RootNodeID].NumberofLeaves;
    header = GravityNode[RootNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header + k;
        if(GravityCache[leaf].Active){
            int index = GravityCache[leaf].Leaf;

            Pbody[index]->Acc[0] = Pall.GravConst*ForceExportSendRecv[Actives].Pos[0];
            Pbody[index]->Acc[1] = Pall.GravConst*ForceExportSendRecv[Actives].Pos[1];
            Pbody[index]->Acc[2] = Pall.GravConst*ForceExportSendRecv[Actives].Pos[2];
            Pbody[index]->Pot = 0.5*Pall.GravConst*ForceExportSendRecv[Actives].Pos[2]*Pbody[leaf]->Mass;

            Actives ++;
        }
    }

    ForceEndProcedure();


    return;
}


static void CalcGravityPlummer(double Pos[3], double e, double Acc[3], double *Pot){

	double r,r2,acc[3],potential,Dist[3],Dist2[3],EPS2,feps;
	double ir,ir3;

    double sqTheta = SQ(GravityRoot.OpeningAngle);

    acc[0] = acc[1] = acc[2] = potential = 0.e0;
	EPS2 = e*e;
    feps = TreeSofteningFactor*e;

    int RootNodeID = 0;
    int CurrentNodeID = GravityNode[RootNodeID].Children;
    while(CurrentNodeID != RootNodeID){
		for(int k=0;k<3;k++){
			Dist[k] = Pos[k] - GravityNode[CurrentNodeID].Pos[k];
			Dist2[k] = SQ(Dist[k]);
		}
		r2 = Dist2[0] + Dist2[1] + Dist2[2];
		r = sqrt(r2);

        double sqWidth = SQ(GravityRoot.Width*GravityRoot.WidthFactor[GravityNode[CurrentNodeID].Level]);

        if((sqWidth<sqTheta*r2)||(GravityNode[CurrentNodeID].Children == NONE)){
			if(GravityNode[CurrentNodeID].Children == NONE){
                int NumberofLeaves = GravityNode[CurrentNodeID].NumberofLeaves;
                int header = GravityNode[CurrentNodeID].Leaves;
                for(int k=0;k<NumberofLeaves;k++){
                    int leaf = header+k;
					Dist[0] = Pos[0] - GravityCache[leaf].Pos[0];
					Dist[1] = Pos[1] - GravityCache[leaf].Pos[1];
					Dist[2] = Pos[2] - GravityCache[leaf].Pos[2];
					double Distance2 = NORM2(Dist);
					r2 = Distance2 + EPS2;
					if( Distance2>EPS2*TINY ){
						r = sqrt(r2);
                        ir = 1.e0/r;
                        ir3 = ir*ir*ir;
						acc[0] -= GravityCache[leaf].Mass*ir3*Dist[0];
						acc[1] -= GravityCache[leaf].Mass*ir3*Dist[1];
						acc[2] -= GravityCache[leaf].Mass*ir3*Dist[2];
						potential -= GravityCache[leaf].Mass*ir;
					}
                }
			} else {
                r2 = DISTANCE2(Pos,GravityNode[CurrentNodeID].COM);
		        r = sqrt(r2+EPS2);  // soften force and potential
				ir = 1.e0/r;
                ir3 = ir*ir*ir;
				acc[0] -= GravityNode[CurrentNodeID].Mass*ir3*Dist[0];
				acc[1] -= GravityNode[CurrentNodeID].Mass*ir3*Dist[1];
				acc[2] -= GravityNode[CurrentNodeID].Mass*ir3*Dist[2];
				potential -= GravityNode[CurrentNodeID].Mass*ir;
			}
            CurrentNodeID = GravityNode[CurrentNodeID].Next;
        } else {
            CurrentNodeID = GravityNode[CurrentNodeID].Children;
        }
    }

    Acc[0] += acc[0];
    Acc[1] += acc[1];
    Acc[2] += acc[2];
    (*Pot) += potential;

	return ;
}

static struct StructForceImport CalcGravityTree(struct StructForceExport ForceExport){

    double Acc[3] = {0.e0,0.e0,0.e0};
    double Potential = 0.e0;
    double Pos[3] = {ForceExport.Pos[0],ForceExport.Pos[1],ForceExport.Pos[2]};
	double EPS2 = SQ(ForceExport.Eps);
	EPS2 *= SQ(Pall.AdaptiveSofteningFactor);

    double sqTheta = SQ(GravityRoot.OpeningAngle);
    int RootNodeID = 0;
    int CurrentNodeID = GravityNode[RootNodeID].Children;
    while(CurrentNodeID != RootNodeID){
        double Dist[3];
        Dist[0] = Pos[0] - GravityNode[CurrentNodeID].Pos[0];
        Dist[1] = Pos[1] - GravityNode[CurrentNodeID].Pos[1];
        Dist[2] = Pos[2] - GravityNode[CurrentNodeID].Pos[2];
		double r2 = NORM2(Dist);
        // or 2xGrvityNode[CurrentNodeID].DistanceMax; ?
        double sqWidth = SQ(GravityRoot.Width*GravityRoot.WidthFactor[GravityNode[CurrentNodeID].Level]);

        if((sqWidth<sqTheta*r2)||(GravityNode[CurrentNodeID].Children == NONE)){
			if(GravityNode[CurrentNodeID].Children == NONE){

                int NumberofLeaves = GravityNode[CurrentNodeID].NumberofLeaves;
                int header = GravityNode[CurrentNodeID].Leaves;
                for(int k=0;k<NumberofLeaves;k++){
                    int leaf = header+k;
					Dist[0] = Pos[0] - GravityCache[leaf].Pos[0];
					Dist[1] = Pos[1] - GravityCache[leaf].Pos[1];
					Dist[2] = Pos[2] - GravityCache[leaf].Pos[2];
					double Distance2 = NORM2(Dist);
					r2 = Distance2 + EPS2;
					if( Distance2>EPS2*TINY ){
						double r = sqrt(r2);
                        double ir = 1.e0/r;
                        double ir3 = ir*ir*ir;
						Acc[0] -= GravityCache[leaf].Mass*ir3*Dist[0];
						Acc[1] -= GravityCache[leaf].Mass*ir3*Dist[1];
						Acc[2] -= GravityCache[leaf].Mass*ir3*Dist[2];
						Potential -= GravityCache[leaf].Mass*ir;
					}
                }
			} else {
                r2 = DISTANCE2(Pos,GravityNode[CurrentNodeID].COM);
		        double r = sqrt(r2+EPS2);  // soften force and potential
				double ir = 1.e0/r;
                double ir2 = ir*ir;
                double ir3 = ir*ir2;
				Acc[0] -= GravityNode[CurrentNodeID].Mass*ir3*Dist[0];
				Acc[1] -= GravityNode[CurrentNodeID].Mass*ir3*Dist[1];
				Acc[2] -= GravityNode[CurrentNodeID].Mass*ir3*Dist[2];
				Potential -= GravityNode[CurrentNodeID].Mass*ir;
			}
            CurrentNodeID = GravityNode[CurrentNodeID].Next;
        } else {
            CurrentNodeID = GravityNode[CurrentNodeID].Children;
        }
    }

	return (struct StructForceImport){{Acc[0],Acc[1],Acc[2]},Potential};
}
