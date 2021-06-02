#include "config.h"
#include "ForceMisc.h"

struct StructForceExport{
    double Pos[3];
    double Mass;
};

static void CalculateForceFromAnotherDomain(const int RecvThisTime, struct StructForceExport ForceRecv[]);
static void WriteAccPot(void);

void ForceDirect(void){

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    /* Self-Domain Gravity */
    for(int i=0;i<Pall.Ntotal;i++){
        double EPS2 = SQ(Pbody[i]->Eps);
        for(int j=0;j<Pall.Ntotal;j++){
            double xij[3];
            xij[0] = Pbody[i]->Pos[0]-Pbody[j]->Pos[0];
            xij[1] = Pbody[i]->Pos[1]-Pbody[j]->Pos[1];
            xij[2] = Pbody[i]->Pos[2]-Pbody[j]->Pos[2];
            double r2 = NORM2(xij);
            if(r2 < TINY) continue;

            double r = sqrt(r2+EPS2);
            double ir = 1.e0/r; double ir3 = ir*ir*ir;

            Pbody[i]->Acc[0] -= Pbody[j]->Mass*xij[0]*ir3;
            Pbody[i]->Acc[1] -= Pbody[j]->Mass*xij[1]*ir3;
            Pbody[i]->Acc[2] -= Pbody[j]->Mass*xij[2]*ir3;
            Pbody[i]->Pot -= Pbody[j]->Mass*ir;
        }
    }

    struct StructForceExport *ForceExportSend,*ForceExportRecv;

    /* Get Other Domains and Calculate Gravity */
    for(int i=0;i<NProcs-1;i++){
        int SendThisTime,RecvThisTime;
        int SendAlready = 0;
        int RecvAlready = 0;
        int SendData = CommunicationTable[i].SendSize;
        int RecvData = CommunicationTable[i].RecvSize;
        while((SendData != 0)&&(RecvData != 0)){
            // Data Packing 
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
            RecvThisTime = CommunicationTable[i].RecvSize-RecvAlready;
            RecvAlready += RecvThisTime;

            CheckSizeofBufferExportRecv(RecvThisTime,sizeof(struct StructForceExport));
            ForceExportRecv = BufferExportRecv;

            // Communication
            MPI_Sendrecv(ForceExportSend,SendThisTime*sizeof(struct StructForceExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_FORCE_DIRECT_SEND,
                         ForceExportRecv,RecvThisTime*sizeof(struct StructForceExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FORCE_DIRECT_SEND,
                        MPI_COMM_WORLD,&mpi_status);
            
            // Calculation
            CalculateForceFromAnotherDomain(RecvThisTime,ForceExportRecv);

            SendData -= SendThisTime;
            RecvData -= RecvThisTime;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Acc[0] *= Pall.GravConst;
        Pbody[i]->Acc[1] *= Pall.GravConst;
        Pbody[i]->Acc[2] *= Pall.GravConst;
        Pbody[i]->Pot *= 0.5*Pall.GravConst*Pbody[i]->Mass;
    }

    return;
}

static void WriteAccPot(void){

    char name[MaxCharactersInLine];
    FILE *fp;

    sprintf(name,"AccPot.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,name,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
            Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],Pbody[i]->Pot);
    }
    fclose(fp);

    return;
}

static void CalculateForceFromAnotherDomain(const int RecvThisTime, struct StructForceExport ForceRecv[]){

    double r,ir,ir3,xij[3],EPS2;

    for(int i=0;i<Pall.Ntotal;i++){
        EPS2 = SQ(Pbody[i]->Eps);
        for(int j=0;j<RecvThisTime;j++){
            xij[0] = Pbody[i]->Pos[0]-ForceRecv[j].Pos[0];
            xij[1] = Pbody[i]->Pos[1]-ForceRecv[j].Pos[1];
            xij[2] = Pbody[i]->Pos[2]-ForceRecv[j].Pos[2];
            r = NORM2(xij);
            if(r < TINY) continue;

            r += EPS2;
            r = sqrt(r+EPS2);
            ir = 1.e0/r; ir3 = ir*ir*ir;

            Pbody[i]->Acc[0] -= ForceRecv[j].Mass*xij[0]*ir3;
            Pbody[i]->Acc[1] -= ForceRecv[j].Mass*xij[1]*ir3;
            Pbody[i]->Acc[2] -= ForceRecv[j].Mass*xij[2]*ir3;
            Pbody[i]->Pot -= ForceRecv[j].Mass*ir;
        }
    }

    return;
}


void PotentialDirect(void){

    int NumProcs = MPIGetNumProcs();
    double r,r2,xij[3],EPS;
    MPI_Status  mpi_status;

    /* Self-Domain Potential */
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pot = 0.e0;
        for(int j=0;j<Pall.Ntotal;j++){
            xij[0] = Pbody[i]->Pos[0]-Pbody[j]->Pos[0];
            xij[1] = Pbody[i]->Pos[1]-Pbody[j]->Pos[1];
            xij[2] = Pbody[i]->Pos[2]-Pbody[j]->Pos[2];
            r2 = NORM2(xij);
            EPS = 0.5*(Pbody[i]->Eps+Pbody[j]->Eps);
            //if(r2 < SQ(TINY*EPS)) continue;

            r = sqrt(r2+SQ(EPS));
            Pbody[i]->Pot -= Pbody[j]->Mass/r;
        }
    }

    struct StructPotential{
        double Pos[3];
        double Eps;
        double Mass;
    } *PotentialExportSend,*PotentialExportRecv;

    /* Get Other Domains and Calculate Gravity */
    for(int i=0;i<NumProcs-1;i++){
        int SendThisTime,RecvThisTime;
        int SendAlready = 0;
        int RecvAlready = 0;
        int SendData = CommunicationTable[i].SendSize;
        int RecvData = CommunicationTable[i].RecvSize;
        while((SendData != 0)&&(RecvData != 0)){
            SendThisTime = CommunicationTable[i].SendSize-SendAlready;
            CheckSizeofBufferExportSendIndex(SendThisTime,sizeof(struct StructPotential),0);
            PotentialExportSend = BufferExportSend[0];

            int NActives = 0;
            while(NActives != SendThisTime){
                PotentialExportSend[NActives].Pos[0] = Pbody[NActives+SendAlready]->Pos[0];
                PotentialExportSend[NActives].Pos[1] = Pbody[NActives+SendAlready]->Pos[1];
                PotentialExportSend[NActives].Pos[2] = Pbody[NActives+SendAlready]->Pos[2];
                PotentialExportSend[NActives].Mass = Pbody[NActives+SendAlready]->Mass;
                PotentialExportSend[NActives].Eps = Pbody[NActives+SendAlready]->Eps;
                NActives ++;
            }
            RecvThisTime = CommunicationTable[i].RecvSize-RecvAlready;
            CheckSizeofBufferExportRecv(RecvThisTime,sizeof(struct StructPotential));
            PotentialExportRecv = BufferExportRecv;

            // Communication
            MPI_Sendrecv(PotentialExportSend,SendThisTime*sizeof(struct StructPotential),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_POTENTIAL_DIRECT_SENDRECV,
                PotentialExportRecv,RecvThisTime*sizeof(struct StructPotential),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_POTENTIAL_DIRECT_SENDRECV,
                        MPI_COMM_WORLD,&mpi_status);

            // Calculation
            for(int j=0;j<Pall.Ntotal;j++){
                for(int k=0;k<RecvThisTime;k++){
                    xij[0] = Pbody[j]->Pos[0]-PotentialExportRecv[k].Pos[0];
                    xij[1] = Pbody[j]->Pos[1]-PotentialExportRecv[k].Pos[1];
                    xij[2] = Pbody[j]->Pos[2]-PotentialExportRecv[k].Pos[2];
                    r2 = NORM2(xij);
                    EPS = 0.5*(Pbody[j]->Eps+PotentialExportRecv[k].Eps);
                    //if(r2 < SQ(TINY*EPS)) continue;

                    r = sqrt(r2+SQ(EPS));
                    Pbody[j]->Pot -= PotentialExportRecv[k].Mass/r;
                }
            }

            SendData -= SendThisTime;
            RecvData -= RecvThisTime;
            SendAlready += SendThisTime;
            RecvAlready += RecvThisTime;
        }
    }


    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pot += Pbody[i]->Mass/Pbody[i]->Eps;
        Pbody[i]->Pot *= 0.5*Pall.GravConst*Pbody[i]->Mass;
    }

    Delete(PotentialExportSend);
    Delete(PotentialExportRecv);

    return;
}

#ifdef USE_SYMMETRIZED_SOFTENING
#include <gp5util.h>
void CalcSymmetrizedPotential(void){

    g5_open();

    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pot = 0.e0;
    }

    //int Npipes = 4;
    //double FieldPos[Pall.Ntotal][3],FieldMass[Pall.Ntotal],FieldEps2[Pall.Ntotal],FieldEps[Pall.Ntotal];
    double (*FieldPos)[3],*FieldMass,*FieldEps2,*FieldEps;
    
    FieldPos = malloc(sizeof(double)*3*Pall.Ntotal);
    FieldMass = malloc(sizeof(double)*Pall.Ntotal);
    FieldEps = malloc(sizeof(double)*Pall.Ntotal);
    FieldEps2 = malloc(sizeof(double)*Pall.Ntotal);

    for(int i=0;i<Pall.Ntotal;i++){
        FieldPos[i][0] = Pbody[i]->Pos[0];
        FieldPos[i][1] = Pbody[i]->Pos[1];
        FieldPos[i][2] = Pbody[i]->Pos[2];
        FieldMass[i] = Pbody[i]->Mass;
        FieldEps2[i] = SQ(Pbody[i]->Eps);
        FieldEps[i] = Pbody[i]->Eps;
    }
    g5_set_n((int)Pall.Ntotal);
#if defined(HAVE_AVX_PHANTOM_GRAPE) //{
    g5_set_xmj0(0,Pall.Ntotal,FieldPos,FieldMass,FieldEps2);
#elif defined(HAVE_AVX_PHANTOM_GRAPE_API2) //}//{
    g5_set_xmj(0,Pall.Ntotal,FieldPos,FieldMass,FieldEps2);
#else //defined(HAVE_AVX_PHANTOM_GRAPE)//}//{
    g5_set_xmeps2j(0,Pall.Ntotal,FieldPos,FieldMass,FieldEps2);
#endif //}

    //double Acc[Pall.Ntotal][3],Pot[Pall.Ntotal];
    double (*Acc)[3],*Pot;
    Acc = malloc(sizeof(double)*3*Pall.Ntotal);
    Pot = malloc(sizeof(double)*Pall.Ntotal);

#if defined(HAVE_AVX_PHANTOM_GRAPE) //{
    g5_calculate_force_on_x0(FieldPos,Acc,Pot,Pall.Ntotal,FieldEps2);
#elif defined(HAVE_AVX_PHANTOM_GRAPE_API2) //{
    g5_calculate_force_on_xe(FieldPos,FieldEps2,Acc,Pot,Pall.Ntotal);
    //g5_calculate_force_on_x0(xi, ai, pi, NI, epsi2);
#else // defined(HAVE_AVX_PHANTOM_GRAPE) //}//{
    g5_calculate_force_on_xeps(FieldPos,FieldEps,Acc,Pot,Pall.Ntotal);
#endif //defined(HAVE_AVX_PHANTOM_GRAPE) //}
    /*
    double adummy[Npipes][3],phidummy[Npipes];
    double time = GetElapsedTime();
    for(int i=0;i<Pall.Ntotal;i+=Npipes){
        g5_set_xi(Npipes,FieldPos+i);
        g5_set_eps(Npipes,FieldEps+i);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
    }
    */


#ifdef USE_SYMMETRIZED_SOFTENING
    double SymmetrizedFactor = 1.0/sqrt(2.0);
#endif
    for(int i=0;i<Pall.Ntotal;i++){
#ifdef USE_SYMMETRIZED_SOFTENING
        Pbody[i]->Pot = Pot[i] - SymmetrizedFactor*(Pbody[i]->Mass/(Pbody[i]->Eps));
#else
        Pbody[i]->Pot = Pot[i] - (Pbody[i]->Mass/(Pbody[i]->Eps));
#endif
        Pbody[i]->Pot = -0.5*Pall.GravConst*Pbody[i]->Mass*Pbody[i]->Pot;
    }

    free(FieldPos);
    free(FieldMass);
    free(FieldEps);
    free(FieldEps2);

    free(Acc);
    free(Pot);

    g5_close();

    return ;
}
#endif // USE_SYMMETRIZED_SOFTENING
