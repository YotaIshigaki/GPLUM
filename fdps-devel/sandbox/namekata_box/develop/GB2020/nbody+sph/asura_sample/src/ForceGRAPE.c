#include "config.h"
#include "ForceMisc.h"
#if ( defined(HAVE_GRAPE7) \
   || defined(HAVE_GRAPE6A) \
   || defined(HAVE_GRAPE5) \
   || defined(HAVE_PHANTOM_GRAPE) \
   || defined(HAVE_AVX_PHANTOM_GRAPE) \
   || defined(HAVE_AVX_PHANTOM_GRAPE_API2))
#include <gp5util.h>
#else
#include "GRAPEEmulator.h"
#define HAVE_GRAPE_EMULATOR
#endif

struct StructForceExport{
    double Pos[3];
    double Mass;
};

static int Npipes;
static int Nboards;
#ifdef HAVE_GRAPE7
static int JMEMSIZE;
#endif

static void CalculateForceGRAPEMyDomain(void);
static void CalculateForceGRAPEAnotherDomain(const int RecvThisTime, struct StructForceExport ForceExportRecv[]);
static void WriteAccPotGrape(void);

void InitializeGRAPE(void){

#if (defined(HAVE_GRAPE7) \
  || defined(HAVE_GRAPE6A) \
  || defined(HAVE_GRAPE5) \
  || defined(HAVE_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE) \
  || defined(HAVE_AVX_PHANTOM_GRAPE_API2))
	g5_open();
    Npipes = g5_get_number_of_pipelines();
	Nboards = 1;
#ifdef HAVE_GRAPE7
    JMEMSIZE = g5_get_jmemsize();
#endif
	g5_close();
#else
	g5_open_emu();
    Npipes = g5_get_number_of_pipelines_emu();
	Nboards = 1;
	g5_close_emu();
    InitializeGRAPEEmulator();
#endif

    return;
}

void SetGRAPE(void){

	double xmin,ymin,zmin,xmax,ymax,zmax,minimum,maximum;
	double rmax,rmin,mmin,size;
    double Wmax,Wmin;

    //sprintlmpi("Start SetGRAPE");

	xmin = Pbody[0]->PosP[0]; xmax = Pbody[0]->PosP[0];
	ymin = Pbody[0]->PosP[1]; ymax = Pbody[0]->PosP[1];
	zmin = Pbody[0]->PosP[2]; zmax = Pbody[0]->PosP[2];
	mmin = Pbody[0]->Mass;
    for(int i=0;i<Pall.Ntotal;i++){
		xmin = MIN(xmin,Pbody[i]->PosP[0]); xmax = MAX(xmax,Pbody[i]->PosP[0]);
		ymin = MIN(ymin,Pbody[i]->PosP[1]); ymax = MAX(ymax,Pbody[i]->PosP[1]);
		zmin = MIN(zmin,Pbody[i]->PosP[2]); zmax = MAX(zmax,Pbody[i]->PosP[2]);
		mmin = MIN(mmin,Pbody[i]->Mass);
	}

	Wmin = MIN(MIN(xmin,ymin),zmin);
	Wmax = MAX(MAX(xmax,ymax),zmax);
    /* Communication and Get Real max and min */
    MPI_Allreduce(&Wmax,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&Wmin,&minimum,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);


	size = 2.01*MAX(ABS(maximum),ABS(minimum));
	rmax = 2.01*MAX(ABS(maximum),ABS(minimum));
	rmin = 2.e-7*rmax;
	mmin *= 1.e-2; 

#ifdef HAVE_GRAPE_EMULATOR 
	g5_set_range_emu(-size,size,mmin);
#else
#if (!defined(HAVE_PHANTOM_GRAPE) \
  && !defined(HAVE_AVX_PHANTOM_GRAPE) \
  && !defined(HAVE_AVX_PHANTOM_GRAPE_API2)) //{
	g5_set_range(-size,size,mmin);
#endif //}
#endif

	return;
}

void ForceGRAPE(void){

    int num_procs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    //sprintl("FORCE GRAPE");
    ClearGravitationalForce();

#ifdef HAVE_GRAPE_EMULATOR
	g5_open_emu();
#else
	g5_open();
#endif
	SetGRAPE();

    //sprintlmpi("Start Mydomain");
    /* Direct GRAPE */
    CalculateForceGRAPEMyDomain();
    //sprintlmpi("End Mydomain");

    struct StructForceExport *ForceExportSend,*ForceExportRecv;

    /* Get Other Domains and Calculate Gravity */
    for(int i=0;i<num_procs-1;i++){
        dprintlmpi(i);
        int SendThisTime,RecvThisTime;
        int SendAlready = 0;
        int RecvAlready = 0;
        int SendData = CommunicationTable[i].SendSize;
        int RecvData = CommunicationTable[i].RecvSize;
        while((SendData != 0)&&(RecvData != 0)){
            SendThisTime = CommunicationTable[i].SendSize-SendAlready;
            CheckSizeofBufferExportSendIndex(SendThisTime,sizeof(struct StructForceExport),0);
            ForceExportSend = BufferExportSend[0];

            for(int k=0;k<SendThisTime;k++){
                ForceExportSend[k].Pos[0] = Pbody[k+SendAlready]->Pos[0];
                ForceExportSend[k].Pos[1] = Pbody[k+SendAlready]->Pos[1];
                ForceExportSend[k].Pos[2] = Pbody[k+SendAlready]->Pos[2];
                ForceExportSend[k].Mass = Pbody[k+SendAlready]->Mass;
            }
            SendAlready += SendThisTime;
            RecvThisTime = CommunicationTable[i].RecvSize-RecvAlready;;
            RecvAlready += RecvThisTime;
            CheckSizeofBufferExportRecv(RecvThisTime,sizeof(struct StructForceExport));
            ForceExportRecv = BufferExportRecv;

            // Communication
            MPI_Sendrecv(ForceExportSend,SendThisTime*sizeof(struct StructForceExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_FORCE_GRAPE_SEND,
                         ForceExportRecv,RecvThisTime*sizeof(struct StructForceExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FORCE_GRAPE_SEND,
                        MPI_COMM_WORLD,&mpi_status);
            
            // Calculation
            CalculateForceGRAPEAnotherDomain(RecvThisTime,ForceExportRecv);

            SendData -= SendThisTime;
            RecvData -= RecvThisTime;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#ifdef HAVE_GRAPE_EMULATOR
    g5_close_emu();
#else
    g5_close();
#endif

    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Acc[0] *= Pall.GravConst;
        Pbody[i]->Acc[1] *= Pall.GravConst;
        Pbody[i]->Acc[2] *= Pall.GravConst;
        Pbody[i]->Pot -= (Pbody[i]->Mass/(Pbody[i]->Eps));
        Pbody[i]->Pot = -0.5*Pall.GravConst*Pbody[i]->Mass*Pbody[i]->Pot;
    }

    ForceEndProcedure();

    return;
}

static void CalculateForceGRAPEMyDomain(void){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    double xj[JMEMSIZE][3],mj[JMEMSIZE];

	for(int k=0;k<Pall.Ntotal;k+=JMEMSIZE){
		int active_j = MIN(JMEMSIZE,Pall.Ntotal-k);

        for(int j=0;j<active_j;j++){
            xj[j][0] = Pbody[j+k]->PosP[0];
            xj[j][1] = Pbody[j+k]->PosP[1];
            xj[j][2] = Pbody[j+k]->PosP[2];
            mj[j] = Pbody[j+k]->Mass;
        }

#if defined(HAVE_GRAPE7)
		g5_set_n(active_j);
		g5_set_jp(0,active_j,mj,xj);
#elif (defined(HAVE_GRAPE6A) || defined(HAVE_GRAPE5) || defined(HAVE_PHANTOM_GRAPE) || defined(HAVE_AVX_PHANTOM_GRAPE))
		g5_set_n(active_j);
		g5_set_xmj(0,active_j,xj,mj);
#else
		g5_set_n_emu(active_j);
		g5_set_xmj_emu(0,active_j,xj,mj);
#endif

		for(int i=0;i<Pall.Ntotal;i+=Npipes){
			int active_i = MIN(Npipes,Pall.Ntotal-i);

			for(int l=0;l<active_i;l++){
                Dxi[l][0] = Pbody[i+l]->PosP[0];
                Dxi[l][1] = Pbody[i+l]->PosP[1];
                Dxi[l][2] = Pbody[i+l]->PosP[2];
                Depsi[l] = Pbody[i+l]->Eps;
            }
			for(int l=active_i;l<Npipes;l++){
                Dxi[l][0] = Dxi[active_i-1][0];
                Dxi[l][1] = Dxi[active_i-1][1];
                Dxi[l][2] = Dxi[active_i-1][2];
                Depsi[l] = Depsi[active_i-1];
            }

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
#elif (defined(HAVE_PHANTOM_GRAPE) || defined(HAVE_AVX_PHANTOM_GRAPE))
            g5_set_eps_to_all(Depsi[0]);
            g5_set_xi(Npipes,Dxi);
            g5_run();
            g5_get_force(Npipes,adummy,phidummy);
#else
            g5_set_ip_emu(Npipes,Dxi,Depsi,Depsi);
            g5_run_emu();
            g5_get_force_emu(Npipes,adummy,phidummy);
#endif

			for(int l=0;l<active_i;l++){
				Pbody[i+l]->Acc[0] += adummy[l][0];
				Pbody[i+l]->Acc[1] += adummy[l][1];
				Pbody[i+l]->Acc[2] += adummy[l][2];
				Pbody[i+l]->Pot += phidummy[l];
			}
		}
	}

    return;
}

static void CalculateForceGRAPEAnotherDomain(const int RecvThisTime, struct StructForceExport ForceExportRecv[]){

	double adummy[Npipes][3],phidummy[Npipes];
	double Dxi[Npipes][3],Depsi[Npipes];
    double xj[JMEMSIZE][3],mj[JMEMSIZE];

	for(int k=0;k<RecvThisTime;k+=JMEMSIZE){
		int active_j = MIN(JMEMSIZE,RecvThisTime-k);

        for(int j=0;j<active_j;j++){
            xj[j][0] = ForceExportRecv[j+k].Pos[0];
            xj[j][1] = ForceExportRecv[j+k].Pos[1];
            xj[j][2] = ForceExportRecv[j+k].Pos[2];
            mj[j] = ForceExportRecv[j+k].Mass;
        }

#if defined(HAVE_GRAPE7)
		g5_set_n(active_j);
		g5_set_jp(0,active_j,mj,xj);
#elif (defined(HAVE_GRAPE6A) || defined(HAVE_GRAPE5) || defined(HAVE_PHANTOM_GRAPE) || defined(HAVE_AVX_PHANTOM_GRAPE))
		g5_set_n(active_j);
		g5_set_xmj(0,active_j,xj,mj);
#else
		g5_set_n_emu(active_j);
		g5_set_xmj_emu(0,active_j,xj,mj);
#endif

		for(int i=0;i<Pall.Ntotal;i+=Npipes){
			int active_i = MIN(Npipes,Pall.Ntotal-i);

			for(int l=0;l<active_i;l++){
                Dxi[l][0] = Pbody[i+l]->PosP[0];
                Dxi[l][1] = Pbody[i+l]->PosP[1];
                Dxi[l][2] = Pbody[i+l]->PosP[2];
                Depsi[l] = Pbody[i+l]->Eps;
            }
			for(int l=active_i;l<Npipes;l++){
                Dxi[l][0] = Dxi[active_i-1][0];
                Dxi[l][1] = Dxi[active_i-1][1];
                Dxi[l][2] = Dxi[active_i-1][2];
                Depsi[l] = Depsi[active_i-1];
            }

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
#elif (defined(HAVE_PHANTOM_GRAPE) || defined(HAVE_AVX_PHANTOM_GRAPE))
            g5_set_eps_to_all(Depsi[0]);
            g5_set_xi(Npipes,Dxi);
            g5_run();
            g5_get_force(Npipes,adummy,phidummy);
#else
            g5_set_ip_emu(Npipes,Dxi,Depsi,Depsi);
            g5_run_emu();
            g5_get_force_emu(Npipes,adummy,phidummy);
#endif

			for(int l=0;l<active_i;l++){
				Pbody[i+l]->Acc[0] += adummy[l][0];
				Pbody[i+l]->Acc[1] += adummy[l][1];
				Pbody[i+l]->Acc[2] += adummy[l][2];
				Pbody[i+l]->Pot += phidummy[l];
			}
		}
	}
    return;
}

static void WriteAccPotGrape(void){

    char name[MaxCharactersInLine];
    FILE *fp;

    sprintf(name,"GrapeAccPot.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,name,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
            Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],Pbody[i]->Pot);
    }
    fclose(fp);

    return;
}
