#include "config.h"
#include "PlantHydroTree.h"
#include "NeighborSearch.h"
#include "IMFParameters.h"
#include "IMF.h"

struct StructDelayedFBExport{
    int       Index;    // Index.
    double    FBRadius; // FBRadius size.
    double    Pos[3];   // Position.
};

struct StructDelayedFBImport{
    int    Index;    // Index.
    int    FBNlist;  // Nlist.
    double FBRho;    // Rho.
#ifdef SET_SNII_TEMPERATURE
    double GasMass;     // Local gas mass
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    unsigned long int  FBDistanceMinID;
    double FBDistanceMin;
#endif
};

struct StructDelayedFBInjection{
    int Mode;        // TypeII? or TypeIa? 
#ifdef MAXIMUM_ENERGY_INPUT
    unsigned long int TargetID;
#endif
    double Radius; // FBRadius size.
    double Rho;    // Rho.
#ifdef SET_SNII_TEMPERATURE
    double GasMass;// Local gas mass.
#endif //SET_SNII_TEMPERATURE
    double Pos[3]; // Position.
    double ZstarII;// Metallicity of this star particle.
    double ZstarIa;// Metallicity of this star particle.
    double Qheat;
    double Mass;
    double Metal;
};

static struct StructDelayedFBImport ReturnStructureDelayedFBImport(double Pos[restrict], const double Kerneli);
static void EnergyMassMetalInjectionEngine(struct StructDelayedFBInjection FBInjection, const int fff);


#define TSNIIstart  (3.3e+6) // yr
#define TSNIIend    (4.5e+7) // yr

double SNIINumber; // per solar mass.
static double SNIIEnergy; // energy per solar mass.
static double SNIIEjector; // Ejector mass per solar mass.
static double SNIIMetal;    // Metal mass per solar mass.
static double iSNIINumber;

static int MaxIndexList = 0;
struct StructDelayedFBActiveStars{
    int Index;
    int Nlist;
    double Radius;
    double Density;
#ifdef SET_SNII_TEMPERATURE
    double GasMass;// Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#if (defined(PRESERVE_SNII_EVENTRATE) || defined(SET_SNII_TEMPERATURE))
    double ReleasedEnergy; 
#endif
    bool ActiveFlag;
    double Rvalue;
    double Lvalue;
} *DelayedFBActiveStars;

static bool first = true;
void InitializeDelayedSNII(void){

    if(!first){
        return ;
    }
    first = false;

    InitializeIMF();
    SNIINumber = IMFSNIINumberPerMass(); 
    SNIIEnergy = SNIIEnergyEfficiency*SNIIEnergyPerNumber*SNIINumber; // unit in ergs.
    SNIIEjector = IMFSNIIEjectaMassPerMass(); // unit in solar mass.
    SNIIMetal = IMFSNIIYieldMassPerMass();  // unit in solar mass.

    iSNIINumber = 1.e0/SNIINumber;

    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"SNII = %g [N/Msun]\n",SNIINumber);
        fprintf(stderr,"SNII energy = %g [erg/Msun]; SNII ENERGY = %g [UnitE/Msun]\n",
                SNIIEnergy,SNIIEnergy*GetUnitEnergy());
        fprintf(stderr,"Ejector Mass = %g [Msun]; Ejector MASS = %g [UnitMass]\n",
                SNIIEjector,SNIIEjector*MSUN_CGS/Pall.UnitMass);
        fprintf(stderr,"Metal Mass = %g [Msun]; Metal MASS = %g [UnitMass]\n",
                SNIIMetal,SNIIMetal*MSUN_CGS/Pall.UnitMass);

        //fprintf(stderr,"SNII = %g [N/Msun], SNII energy = %g [erg/Msun], Ejector Mass = %g [Msun], Metal Mass = %g [Msun]\n",
                //SNIINumber,SNIIEnergy,SNIIEjector,SNIIMetal);
    }

    SNIIEnergy *= GetUnitEnergy();
    //SNIIEjector *= (MSUN_CGS/Pall.UnitMass);
    //SNIIMetal *= (MSUN_CGS/Pall.UnitMass);

    //if(MPIGetMyID()==MPI_ROOT_RANK)
        //fprintf(stderr,"SNII = %g [N/Msun], SNII ENERGY = %g [UnitE/Msun], Ejector MASS = %g [UnitMass], Metal MASS = %g [UnitMass]\n",
                //SNIINumber,SNIIEnergy,SNIIEjector*MSUN_CGS/Pall.UnitMass,SNIIMetal*MSUN_CGS/Pall.UnitMass);

    InitializeLifeTimeLookUpTable();
    if(MPIGetMyID()==MPI_ROOT_RANK){
        gprintlmpi( GetDyingStarMassFromAge(1.e+6, 0.e0) );
        gprintlmpi( GetDyingStarMassFromAge(2.e+6, 0.e0) );
        gprintlmpi( GetDyingStarMassFromAge(3.e+6, 0.e0) );
        gprintlmpi( GetDyingStarMassFromAge(5.e+6, 0.e0) );
        gprintlmpi( GetDyingStarMassFromAge(1.e+7, 0.e0) );
        gprintlmpi( GetDyingStarMassFromAge(3.e+7, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(100, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(80, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(60, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(40, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(20, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(10, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(8, 0.e0) );

        gprintlmpi( GetLifeTimeFromStellarMass(100, 0.e0) );
        gprintlmpi( GetLifeTimeFromStellarMass(100, 0.02) );
        gprintlmpi( GetLifeTimeFromStellarMass(100, 0.1) );
        gprintlmpi( GetLifeTimeFromStellarMass(100, 1.0) );
    }
    fflush(NULL);

#if defined(PRESERVE_SNII_EVENTRATE)
    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"Use PRESERVE_SNII_EVENTRATE option.\n");
    }
    fflush(NULL);
#   ifdef SET_SNII_TEMPERATURE
#   error Exclusive option "SET_SNII_TEMPERATURE" is used.
#   endif
#elif defined(SET_SNII_TEMPERATURE)
    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"Use SET_SNII_TEMPERATURE option.\n");
        fprintf(stderr,"SNII_TEMPERATURE is set to %g [K].\n",SNII_TEMPERATURE);
    }
    fflush(NULL);
#   ifdef PRESERVE_SNII_EVENTRATE
#   error Exclusive option "PRESERVE_SNII_EVENTRATE" is used.
#   endif
#endif
    
    MaxIndexList = NAdditionUnit;
    DelayedFBActiveStars = 
        realloc(DelayedFBActiveStars,sizeof(struct StructDelayedFBActiveStars)*MaxIndexList);

    return;
}


static inline bool CheckCandidateSNII(const int Index) __attribute__((always_inline));
static inline bool CheckCandidateSNII(const int Index){

    if(Pstar[Index]->TypeII == true)
        return false;

    // Get the age of this star particle.
    double Age = (Pall.TCurrent-Pstar[Index]->FormationTime)*Pall.UnitTime/YEAR_CGS;

    // If the age of this star particle is sufficiently long (i.e. Age >
    // Age(MSNIIMin)), SNe should take place in this time-step.
    if(Age > GetLifeTimeFromStellarMass(MSNIIMin,Pstar[Index]->Z))
        return true;

    // If the age of this star is shorter than that with MSNIIMax, this star is
    // unsuitable for SNe.
    if(Age < GetLifeTimeFromStellarMass(MSNIIMax,Pstar[Index]->Z))
        return false;

    // rate and random explosion.
    double Age_old = (Pall.TCurrent-Pstar[Index]->FormationTime-PstarBody(Index)->dt)*Pall.UnitTime/YEAR_CGS;
    if(Age_old < 0.e0){
        fprintf(stderr,"Age_old = %g\n",Age_old);
        return false;
    }

    // Get the mass of the star just dying this epoch.
    double Mstar = GetDyingStarMassFromAge(Age,Pstar[Index]->Z);

    // If Mstar > MSNIIMax, this star is sufficiently young.
    if((Mstar>MSNIIMax)||(Mstar < 0.e0))
        return false;

    double Mstar_old = GetDyingStarMassFromAge(Age_old,Pstar[Index]->Z);
    if((Mstar_old < 0.e0)||(Mstar_old > MSNIIMax))
        Mstar_old = MSNIIMax;
    if(Mstar_old-Mstar <0.e0){
        fprintf(stderr,"Mstar %g, Mstar_old %g\n",Mstar,Mstar_old);
        fprintf(stderr,"Why the mass(Age_old) is larger than the mass(Age)?\n");
        fflush(NULL);
        assert(Mstar_old-Mstar>0.e0);
    }

    double dR = IMFSNIINumberPerMassLimitedRange(Mstar_old,Mstar);
    double R = IMFSNIINumberPerMassLimitedRange(Mstar_old,MSNIIMin);
    //ResetRandomSeedForRandomGenerator(PstarBody(Index)->GlobalID); // FB test
    if(gsl_rng_uniform(RandomGenerator)<dR/R){
#ifdef PRINT_LOG_DELAYED_FEEDBACK
        fprintf(stderr,"FB time has come: %g, n dR %g dR %g |  %g <%g> %g\n",Age,
            dR/R,dR,
            GetLifeTimeFromStellarMass(MSNIIMax,Pstar[Index]->Z)*YEAR_CGS/Pall.UnitTime,Mstar,
            GetLifeTimeFromStellarMass(MSNIIMin,Pstar[Index]->Z)*YEAR_CGS/Pall.UnitTime);
#endif
        return true;
    }

    return false;
}

static int GetSNeNumber(void){

    int NSNe = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarActive(i)){
            if(CheckCandidateSNII(i)){ // FB test
                if(MaxIndexList <= NSNe+1){
                    MaxIndexList += NAdditionUnit;
                    DelayedFBActiveStars = 
                        realloc(DelayedFBActiveStars,sizeof(struct StructDelayedFBActiveStars)*MaxIndexList);
                }
                DelayedFBActiveStars[NSNe].Index = i;
                DelayedFBActiveStars[NSNe].Radius = PstarBody(i)->Eps;
                NSNe ++;
            }
        }
    }

    
    return NSNe;
}

static inline bool OverlapDomainDelayed(double Pos[restrict], 
        const double h, const int NodeID) __attribute__((always_inline));
static inline bool OverlapDomainDelayed(double Pos[restrict], const double h, const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline double DomainDistanceSQR(double Pos[restrict], const int NodeID) __attribute__((always_inline));
static inline double DomainDistanceSQR(double Pos[restrict], const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2);
}


static inline int CheckDelayedExportFlags(const int TargetNodeID, const int NLocalTarget,
            const int NProcs, bool DelayedFBExportFlags[][NProcs-1]) __attribute__((always_inline));
static inline int CheckDelayedExportFlags(const int TargetNodeID, const int NLocalTarget,
            const int NProcs, bool DelayedFBExportFlags[][NProcs-1]){

    int NodeID = CommunicationTable[TargetNodeID].SendRank;

    int NExport = 0;
    for(int i=0;i<NLocalTarget;i++){
        if(DelayedFBActiveStars[i].ActiveFlag){
            int index = DelayedFBActiveStars[i].Index; 
            if(OverlapDomainDelayed(PstarPosP(index),2.0*DelayedFBActiveStars[i].Radius,NodeID)){
                DelayedFBExportFlags[i][TargetNodeID] = ON;
                NExport ++;
            }
        }
    }

	return NExport;
}

static int CheckFBRadiusAndReturnLocalActiceLeaves(const int NLocalTarget, struct StructDelayedFBActiveStars DelayedFBActiveStars[]){

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#define dConverge   (1.e-6)

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NLocalTarget;i++){
        if(DelayedFBActiveStars[i].ActiveFlag){ 
            int Nlist = DelayedFBActiveStars[i].Nlist;
            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                DelayedFBActiveStars[i].ActiveFlag = OFF;
            }else if((DelayedFBActiveStars[i].Rvalue>0.e0)&&(DelayedFBActiveStars[i].Lvalue>0.e0)){
                if(DelayedFBActiveStars[i].Rvalue-DelayedFBActiveStars[i].Lvalue < dConverge*DelayedFBActiveStars[i].Lvalue)
                    DelayedFBActiveStars[i].ActiveFlag = OFF;
            }

            if(DelayedFBActiveStars[i].ActiveFlag){
                if(Nlist<NBmin){
                    DelayedFBActiveStars[i].Lvalue = fmax(DelayedFBActiveStars[i].Lvalue,DelayedFBActiveStars[i].Radius);
                } else if(Nlist>NBmax){
                    if(DelayedFBActiveStars[i].Rvalue > 0.e0){
                        DelayedFBActiveStars[i].Rvalue = fmin(DelayedFBActiveStars[i].Rvalue,DelayedFBActiveStars[i].Radius);
                    }else{
                        DelayedFBActiveStars[i].Rvalue = DelayedFBActiveStars[i].Radius;
                    }
                }

                if((DelayedFBActiveStars[i].Lvalue>0.e0)&&(DelayedFBActiveStars[i].Rvalue>0.e0)){
                    DelayedFBActiveStars[i].Radius = cbrt(0.5*(CUBE(DelayedFBActiveStars[i].Lvalue)+CUBE(DelayedFBActiveStars[i].Rvalue)));
                }else{
                    if((DelayedFBActiveStars[i].Rvalue == 0.e0)&&(DelayedFBActiveStars[i].Lvalue > 0.e0)){
                        DelayedFBActiveStars[i].Radius *= 1.26;
                    }else if((DelayedFBActiveStars[i].Rvalue > 0.e0)&&(DelayedFBActiveStars[i].Lvalue == 0.e0)){
                        DelayedFBActiveStars[i].Radius *= 0.74;
                    }
                }
                NLocalActiveLeaves ++;
            }
        }
    }

    return NLocalActiveLeaves;
}

static void GetFeedbackRadius(const int NLocalTarget, const int NProcs, bool DelayedFBExportFlags[][NProcs-1]){

    MPI_Status  mpi_status;

    int NLocalActiveLeaves = NLocalTarget;
    int NActiveLeaves = NLocalTarget; // Current Active Leaves, decrement in iteration. 

    for(int i=0;i<NLocalTarget;i++){
        DelayedFBActiveStars[i].ActiveFlag = ON;
        DelayedFBActiveStars[i].Rvalue = DelayedFBActiveStars[i].Lvalue = 0.e0;
    }

    // counter for export and import.
    int BitMask = 0x01; 
    int NExportSizeAllocated[NProcs];
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];
    for(int i=0;i<NProcs;i++)
        NExportSizeAllocated[i] = 0;

    struct StructDelayedFBExport *DelayedFBExportSend[NProcs];
    struct StructDelayedFBExport *DelayedFBExportRecv = NULL;
    struct StructDelayedFBImport *DelayedFBImportSend = NULL;
    struct StructDelayedFBImport *DelayedFBImportRecv[NProcs];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    MPI_Status mpi_status_Import_Send[NProcs-1];
    MPI_Request mpi_request_Import_Send[NProcs-1];
    MPI_Status mpi_status_Import_Recv[NProcs-1];
    MPI_Request mpi_request_Import_Recv[NProcs-1];

    // start iteration.
    int Niteration = 0;
    do{
        for(int i=0;i<NLocalTarget;i++){ 
            if(DelayedFBActiveStars[i].ActiveFlag == ON){
                DelayedFBActiveStars[i].Nlist = 0;
                DelayedFBActiveStars[i].Density = 0.e0;
#ifdef SET_SNII_TEMPERATURE
                DelayedFBActiveStars[i].GasMass = 0.e0;
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
                DelayedFBActiveStars[i].DistanceMin = 0.e0;   
                DelayedFBActiveStars[i].DistanceMinGlobalID = NONE;
#endif
                for(int k=0;k<NProcs-1;k++)
                    DelayedFBExportFlags[i][k] = 0;
            }
        }

        int NExportMaxThisTime = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime[i] = CheckDelayedExportFlags(i,NLocalTarget,NProcs,DelayedFBExportFlags);
            CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructDelayedFBExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],sizeof(struct StructDelayedFBImport),i);
            DelayedFBExportSend[i] = BufferExportSend[i];
            DelayedFBImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            for(int k=0;k<NLocalTarget;k++){
                if(DelayedFBActiveStars[k].ActiveFlag){
                    int leaf = DelayedFBActiveStars[k].Index;
                    if(DelayedFBExportFlags[k][i]&BitMask){ 
                        DelayedFBExportSend[i][NExport].Index  = k;
                        DelayedFBExportSend[i][NExport].FBRadius = DelayedFBActiveStars[k].Radius;
                        DelayedFBExportSend[i][NExport].Pos[0] = PstarPosP(leaf)[0];
                        DelayedFBExportSend[i][NExport].Pos[1] = PstarPosP(leaf)[1];
                        DelayedFBExportSend[i][NExport].Pos[2] = PstarPosP(leaf)[2];
                        NExport ++;
                    }
                }
            }
            assert(NExport == NExportThisTime[i]);
            NExportMaxThisTime = MAX(NExportMaxThisTime,NExport);
        }

        int NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,
                    CommunicationTable[i].SendRank,TAG_FEEDBACK_DELAYED_PRECOMM,
                         NImportThisTime+i,1,MPI_INT,
                    CommunicationTable[i].RecvRank,TAG_FEEDBACK_DELAYED_PRECOMM,
                        MPI_COMM_WORLD,&mpi_status);
            NImport += NImportThisTime[i];
        }
        int NImportAll = NImport;

        CheckSizeofBufferExportRecv(NImport,sizeof(struct StructDelayedFBExport));
        CheckSizeofBufferImportSend(NImport,sizeof(struct StructDelayedFBImport));
        DelayedFBExportRecv = BufferExportRecv;
        DelayedFBImportSend = BufferImportSend;

        NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            MPI_Isend(DelayedFBExportSend[i],
                NExportThisTime[i]*sizeof(struct StructDelayedFBExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_FEEDBACK_DELAYED_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+i);
            MPI_Irecv(DelayedFBExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructDelayedFBExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FEEDBACK_DELAYED_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+i);
            NImport += NImportThisTime[i];
        }

        // search local
        for(int i=0;i<NLocalTarget;i++){ 
            if(DelayedFBActiveStars[i].ActiveFlag == ON){
                int leaf = DelayedFBActiveStars[i].Index;
                struct StructDelayedFBImport TemporalData = 
                    ReturnStructureDelayedFBImport(PstarPosP(leaf),DelayedFBActiveStars[i].Radius);
                DelayedFBActiveStars[i].Nlist = TemporalData.FBNlist;
                DelayedFBActiveStars[i].Density = TemporalData.FBRho;
#ifdef SET_SNII_TEMPERATURE
                DelayedFBActiveStars[i].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
                if(TemporalData.FBNlist > 0){
                    DelayedFBActiveStars[i].DistanceMin = TemporalData.FBDistanceMin;
                    DelayedFBActiveStars[i].DistanceMinGlobalID = TemporalData.FBDistanceMinID;
                }
#endif
            }
        }
        double TimeComm = GetElapsedTime();
        MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
        MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
        TimingResults.FeedbackCommThisStep += GetElapsedTime()-TimeComm;


        for(int i=0;i<NImportAll;i++){
            struct StructDelayedFBImport TemporalData = 
                ReturnStructureDelayedFBImport(DelayedFBExportRecv[i].Pos,
                        DelayedFBExportRecv[i].FBRadius);
            DelayedFBImportSend[i].Index = DelayedFBExportRecv[i].Index;
            DelayedFBImportSend[i].FBNlist = TemporalData.FBNlist;
            DelayedFBImportSend[i].FBRho = TemporalData.FBRho;
#ifdef SET_SNII_TEMPERATURE
            DelayedFBImportSend[i].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
            if(TemporalData.FBNlist > 0){
                DelayedFBImportSend[i].FBDistanceMinID = TemporalData.FBDistanceMinID;
                DelayedFBImportSend[i].FBDistanceMin = TemporalData.FBDistanceMin;
            } else {
                DelayedFBImportSend[i].FBRho = 0.e0;
                DelayedFBImportSend[i].FBDistanceMinID = NONE;
                DelayedFBImportSend[i].FBDistanceMin = 0.e0;
            }
#endif
        }

        NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            MPI_Isend(DelayedFBImportSend+NImport,
                NImportThisTime[i]*sizeof(struct StructDelayedFBImport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_FEEDBACK_DELAYED_IMPORT+i,
                        MPI_COMM_WORLD,mpi_request_Import_Send+i);
            MPI_Irecv(DelayedFBImportRecv[i],
                NExportThisTime[i]*sizeof(struct StructDelayedFBImport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FEEDBACK_DELAYED_IMPORT+i,
                        MPI_COMM_WORLD,mpi_request_Import_Recv+i);
            NImport += NImportThisTime[i];
        }
        TimeComm = GetElapsedTime();
        MPI_Waitall(NProcs-1,mpi_request_Import_Recv,mpi_status_Import_Recv);
        MPI_Waitall(NProcs-1,mpi_request_Import_Send,mpi_status_Import_Send);
        TimingResults.FeedbackCommThisStep += GetElapsedTime()-TimeComm;

        for(int i=0;i<NProcs-1;i++){
            for(int k=0;k<NExportThisTime[i];k++){ 
                if(DelayedFBImportRecv[i][k].FBNlist>0){
                    int leaf = DelayedFBImportRecv[i][k].Index;
#ifdef MAXIMUM_ENERGY_INPUT
                    if(DelayedFBActiveStars[leaf].Nlist == 0){
                        DelayedFBActiveStars[leaf].DistanceMin = DelayedFBImportRecv[i][k].FBDistanceMin;
                        DelayedFBActiveStars[leaf].DistanceMinGlobalID = DelayedFBImportRecv[i][k].FBDistanceMinID;
                    } else {
                        if(DelayedFBActiveStars[leaf].DistanceMin>DelayedFBImportRecv[i][k].FBDistanceMin){
                            DelayedFBActiveStars[leaf].DistanceMin = DelayedFBImportRecv[i][k].FBDistanceMin;
                            DelayedFBActiveStars[leaf].DistanceMinGlobalID = DelayedFBImportRecv[i][k].FBDistanceMinID;
                        }
                    }
#endif
                    DelayedFBActiveStars[leaf].Nlist += DelayedFBImportRecv[i][k].FBNlist;
                    DelayedFBActiveStars[leaf].Density += DelayedFBImportRecv[i][k].FBRho;
#ifdef SET_SNII_TEMPERATURE
                    DelayedFBActiveStars[leaf].GasMass += DelayedFBImportRecv[i][k].GasMass;
#endif //SET_SNII_TEMPERATURE
                }
            }
        }


        NLocalActiveLeaves = CheckFBRadiusAndReturnLocalActiceLeaves(NLocalTarget,DelayedFBActiveStars);
        MPI_Allreduce(&NLocalActiveLeaves,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        Niteration ++;
    } while (0<NActiveLeaves);

    return;
}

static inline double KernelDelayedFB(const double r, const double InvKerneli) __attribute__((always_inline));
static inline double KernelDelayedFB(const double r, const double InvKerneli){

    double u = r*InvKerneli;
    double coef = IPI*CUBE(InvKerneli);
    if(u<1.e0){
        return (coef*(1.e0 - 1.5*SQ(u) + 0.75*CUBE(u)));
    } else if (u<2.e0){
        return (coef*(0.25*CUBE(2.e0-u)));
    } else {
        return 0.e0;
    }
}

static int Neighbors[MaxNeighborSize];
static struct StructDelayedFBImport ReturnStructureDelayedFBImport(double Pos[restrict], const double Kerneli){

    double InvKerneli = 1.e0/Kerneli;
    //struct StructDelayedFBImport TempDelayedFBImport = InitializeDelayedFBImport;
    struct StructDelayedFBImport TempDelayedFBImport;
    memset(&TempDelayedFBImport,0,sizeof(struct StructDelayedFBImport));

    double TimeComm = GetElapsedTime();
    int Nlist = GetNeighborsLimited(Pos,2.e0*Kerneli,Neighbors);
    //int Nlist2 = GetNeighborsDirect(Pos,2.e0*Kerneli,Neighbors);
    //assert(Nlist == Nlist);
    TimingResults.FeedbackNeighborSearchThisStep += GetElapsedTime()-TimeComm;

    if(Nlist > MaxNeighborSize){
        fprintf(stderr,"The neighbor list is over flow(Nlist = %d). function:%s,line:%d,file:%s\n",
                Nlist,__FUNCTION__,__LINE__,__FILE__);
        MPI_Finalize();
        exit(NeighborListOverFlow);
    } else if (Nlist == 0){
        return TempDelayedFBImport;
    }

    struct StructDelayedFBPreFetch{
        double Mass;
        double Kernel;
        double Pos[3];
#ifdef MAXIMUM_ENERGY_INPUT
        unsigned long int GlobalID;
#endif
    } DelayedFBPreFetch[Nlist];

	for(int k=0;k<Nlist;k++){
		int leaf = Neighbors[k];
        DelayedFBPreFetch[k].Mass = PhydroMass(leaf);
        DelayedFBPreFetch[k].Kernel = Phydro[leaf]->KernelPred;
        DelayedFBPreFetch[k].Pos[0] = PhydroPosP(leaf)[0];
        DelayedFBPreFetch[k].Pos[1] = PhydroPosP(leaf)[1];
        DelayedFBPreFetch[k].Pos[2] = PhydroPosP(leaf)[2];
#ifdef MAXIMUM_ENERGY_INPUT
        DelayedFBPreFetch[k].GlobalID = PhydroBody(leaf)->GlobalID;
#endif
    }

	for(int k=0;k<Nlist;k++){
        double xij[3];
#ifdef PERIODIC_RUN
        xij[0] = PeriodicDistance(Pos[0],DelayedFBPreFetch[k].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],DelayedFBPreFetch[k].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],DelayedFBPreFetch[k].Pos[2],2);
#else // PERIODIC_RUN
        xij[0] = Pos[0]-DelayedFBPreFetch[k].Pos[0];
        xij[1] = Pos[1]-DelayedFBPreFetch[k].Pos[1];
        xij[2] = Pos[2]-DelayedFBPreFetch[k].Pos[2];
#endif // PERIODIC_RUN

        double r = NORM(xij); 
		double w = KernelDelayedFB(r,InvKerneli);

        TempDelayedFBImport.FBRho += DelayedFBPreFetch[k].Mass*w;
#ifdef SET_SNII_TEMPERATURE
        TempDelayedFBImport.GasMass += DelayedFBPreFetch[k].Mass;
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
        if(k==0){
            TempDelayedFBImport.FBDistanceMinID = DelayedFBPreFetch[k].GlobalID;       
            TempDelayedFBImport.FBDistanceMin = r;       
        } else if (TempDelayedFBImport.FBDistanceMin > r){
            TempDelayedFBImport.FBDistanceMinID = DelayedFBPreFetch[k].GlobalID;       
            TempDelayedFBImport.FBDistanceMin = r;       
        }
#endif
    }
    TempDelayedFBImport.FBNlist = Nlist;

    return TempDelayedFBImport;
}


static void EnergyMassMetalInjectionEngine(struct StructDelayedFBInjection FBInjection, const int fff){

    double InvKerneli = 1.e0/FBInjection.Radius;

    double TimeComm = GetElapsedTime();
    int Nlist = GetNeighborsLimited(FBInjection.Pos,2.e0*FBInjection.Radius,Neighbors);
    TimingResults.FeedbackNeighborSearchThisStep += GetElapsedTime()-TimeComm;

    if(Nlist > MaxNeighborSize){
        fprintf(stderr,"The neighbor list is over flow(Nlist = %d). function:%s,line:%d,file:%s\n",
                Nlist,__FUNCTION__,__LINE__,__FILE__);
        MPI_Finalize();
        exit(NeighborListOverFlow);
    }

    double iRhoi = 1.e0/FBInjection.Rho;
	for(int k=0;k<Nlist;k++){
		int leaf = Neighbors[k];
        double xij[3],Posj[3];

        Posj[0] = PhydroPosP(leaf)[0];
        Posj[1] = PhydroPosP(leaf)[1];
        Posj[2] = PhydroPosP(leaf)[2];
#ifdef PERIODIC_RUN 
        xij[0] = PeriodicDistance(FBInjection.Pos[0],Posj[0],0);
        xij[1] = PeriodicDistance(FBInjection.Pos[1],Posj[1],1);
        xij[2] = PeriodicDistance(FBInjection.Pos[2],Posj[2],2);
#else // PERIODIC_RUN 
        xij[0] = FBInjection.Pos[0]-Posj[0];
        xij[1] = FBInjection.Pos[1]-Posj[1];
        xij[2] = FBInjection.Pos[2]-Posj[2];
#endif // PERIODIC_RUN

        double r = NORM(xij); 
		double w = KernelDelayedFB(r,InvKerneli);

        double Weight = PhydroMass(leaf)*iRhoi*w;
#ifdef MAXIMUM_ENERGY_INPUT
        if(PhydroBody(leaf)->GlobalID == FBInjection.TargetID){
            Phydro[leaf]->DQheat += FBInjection.Qheat;
        } 
#else 
        Phydro[leaf]->DQheat += FBInjection.Qheat*Weight;
#endif

        double dM = FBInjection.Mass*Weight;
        if(FBInjection.Mode == 0){
            Phydro[leaf]->dZII += FBInjection.Metal*Weight+FBInjection.ZstarII*dM;
            Phydro[leaf]->dZIa += FBInjection.ZstarIa*dM;
        } else {
            Phydro[leaf]->dZII += FBInjection.Metal*Weight+FBInjection.ZstarIa*dM;
            Phydro[leaf]->dZIa += FBInjection.ZstarII*dM;
        }
        Phydro[leaf]->dMass += dM;
    }

    return ;
}

#if (defined(PRESERVE_SNII_EVENTRATE) || defined(SET_SNII_TEMPERATURE))
static void EnergyMassMetalInjection(const int NLocalTarget, const int NProcs, bool DelayedFBExportFlags[][NProcs-1]){

    MPI_Status  mpi_status;

    // counter for export and import. 
    int BitMask = 0x01; 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    for(int i=0;i<NLocalTarget;i++){
        int leaf = DelayedFBActiveStars[i].Index;
        if(Pstar[leaf]->TypeII == OFF){ // If SNII case...
            DelayedFBActiveStars[i].ReleasedEnergy = 0.e0;
            double MassInMsun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS;
#if defined(PRESERVE_SNII_EVENTRATE)
            if(Pstar[leaf]->TypeIIProb){
                double prob = SNIINumber*MassInMsun;
                if(prob >= 1.0){
                    DelayedFBActiveStars[i].ReleasedEnergy = SNIIEnergy*MassInMsun;
                } else {
                    DelayedFBActiveStars[i].ReleasedEnergy 
                        = SNIIEnergyEfficiency*SNIIEnergyPerNumber*GetUnitEnergy();
                }
            } else {
                DelayedFBActiveStars[i].ReleasedEnergy = 0.e0;
            }
#if 0
            double prob = SNIINumber*MassInMsun;
            if(prob >= 1.0){
                DelayedFBActiveStars[i].ReleasedEnergy = SNIIEnergy*MassInMsun;
            }else{
                if(prob > gsl_rng_uniform(RandomGenerator)){
                    DelayedFBActiveStars[i].ReleasedEnergy 
                        = SNIIEnergyEfficiency*SNIIEnergyPerNumber*GetUnitEnergy();
                }else{
                    DelayedFBActiveStars[i].ReleasedEnergy = 0.e0;
                }
            }
            fprintf(stderr,"p = %g, %g; %g > %g\n",
                prob,SNIIEnergy*MassInMsun,
                    SNIIEnergyEfficiency*SNIIEnergyPerNumber*GetUnitEnergy(),0.e0); // FB test
#endif
#elif defined(SET_SNII_TEMPERATURE)
            double Usn_in_sim_unit = (SNIIEnergy*MassInMsun)/(DelayedFBActiveStars[i].GasMass);
            double Tsn = Pall.ConvertUtoT*Usn_in_sim_unit;
            double prob = Tsn/SNII_TEMPERATURE;
            if(prob >= 1.0){
                DelayedFBActiveStars[i].ReleasedEnergy = SNIIEnergy*MassInMsun;
            }else{
                if(prob > gsl_rng_uniform(RandomGenerator)){
                    DelayedFBActiveStars[i].ReleasedEnergy 
                        = Pall.ConvertTtoU*SNII_TEMPERATURE*DelayedFBActiveStars[i].GasMass;
                }else{
                    DelayedFBActiveStars[i].ReleasedEnergy = 0.e0;
                }
            }
            fprintf(stderr,"p = %g, Tsn = %g [k], E = %g\n",
                    prob,Tsn,DelayedFBActiveStars[i].ReleasedEnergy);
#endif // defined(PRESERVE_SNII_EVENTRATE)
        }
    }


    struct StructDelayedFBInjection *DelayedFBInjectionSend[NProcs],*DelayedFBInjectionRecv;

    for(int i=0;i<NProcs-1;i++){ // Prepare data for export 
        int NExport = 0;
        for(int k=0;k<NLocalTarget;k++){
            if(DelayedFBExportFlags[k][i]&BitMask)
                NExport ++;
        }
        NExportThisTime[i] = NExport;
        DelayedFBInjectionSend[i] = 
            malloc(sizeof(struct StructDelayedFBInjection)*(NExportThisTime[i]+1));

        NExport = 0;
        for(int k=0;k<NLocalTarget;k++){
            int leaf = DelayedFBActiveStars[k].Index;
            if(DelayedFBExportFlags[k][i]&BitMask){
                DelayedFBInjectionSend[i][NExport].Radius = DelayedFBActiveStars[k].Radius;
                DelayedFBInjectionSend[i][NExport].Rho = DelayedFBActiveStars[k].Density;
                DelayedFBInjectionSend[i][NExport].Pos[0] = PstarPosP(leaf)[0];
                DelayedFBInjectionSend[i][NExport].Pos[1] = PstarPosP(leaf)[1];
                DelayedFBInjectionSend[i][NExport].Pos[2] = PstarPosP(leaf)[2];
#ifdef MAXIMUM_ENERGY_INPUT
                DelayedFBInjectionSend[i][NExport].TargetID = DelayedFBActiveStars[k].DistanceMinGlobalID;
                assert(DelayedFBActiveStars[k].DistanceMinGlobalID != NONE);
#endif
                if(Pstar[leaf]->TypeII == OFF){ // If SNII case...
                    DelayedFBInjectionSend[i][NExport].Qheat = DelayedFBActiveStars[k].ReleasedEnergy;
                    //DelayedFBInjectionSend[i][NExport].Qheat = SNIIEnergy*MassInMsun; //
                    DelayedFBInjectionSend[i][NExport].Mass = SNIIEjector*Pstar[leaf]->InitialMass; 
                    DelayedFBInjectionSend[i][NExport].Metal = SNIIMetal*Pstar[leaf]->InitialMass;
                    DelayedFBInjectionSend[i][NExport].ZstarII = Pstar[leaf]->ZII;
                    DelayedFBInjectionSend[i][NExport].ZstarIa = Pstar[leaf]->ZIa;
                    DelayedFBInjectionSend[i][NExport].Mode = 0;
                } else if(Pstar[leaf]->TypeIa == OFF){  // If SNIa case...
                    double MassInMsun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS;
                    DelayedFBInjectionSend[i][NExport].Qheat = SNIIEnergy*MassInMsun; 
                    DelayedFBInjectionSend[i][NExport].Mass = SNIIEjector*Pstar[leaf]->InitialMass;  
                    DelayedFBInjectionSend[i][NExport].Metal = SNIIMetal*Pstar[leaf]->InitialMass; 
                    DelayedFBInjectionSend[i][NExport].ZstarII = Pstar[leaf]->ZII;
                    DelayedFBInjectionSend[i][NExport].ZstarIa = Pstar[leaf]->ZIa;
                    DelayedFBInjectionSend[i][NExport].Mode = 1;
                }
                NExport ++;
            }
        }
        assert(NExport == NExportThisTime[i]);
    }

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_FEEDBACK_DELAYED_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_FEEDBACK_DELAYED_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }

    DelayedFBInjectionRecv = malloc(sizeof(struct StructDelayedFBInjection)*(NImportAll+1));

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(DelayedFBInjectionSend[i],
            NExportThisTime[i]*sizeof(struct StructDelayedFBInjection),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_FEEDBACK_DELAYED_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(DelayedFBInjectionRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructDelayedFBInjection),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FEEDBACK_DELAYED_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);


    for(int i=0;i<NLocalTarget;i++){  // Feedback from stars in the local domain. 
        int leaf = DelayedFBActiveStars[i].Index;
        struct StructDelayedFBInjection FBInjection;

        FBInjection.Pos[0] = PstarPosP(leaf)[0];
        FBInjection.Pos[1] = PstarPosP(leaf)[1];
        FBInjection.Pos[2] = PstarPosP(leaf)[2];
        FBInjection.Rho = DelayedFBActiveStars[i].Density;
        FBInjection.Radius = DelayedFBActiveStars[i].Radius;
#ifdef MAXIMUM_ENERGY_INPUT
        FBInjection.TargetID = DelayedFBActiveStars[i].DistanceMinGlobalID;
#endif
        if(Pstar[leaf]->TypeII == OFF){ // If SNII case...
            FBInjection.Qheat = DelayedFBActiveStars[i].ReleasedEnergy;
            //FBInjection.Qheat = SNIIEnergy*MassInMsun;
            FBInjection.Mass = SNIIEjector*Pstar[leaf]->InitialMass; 
            FBInjection.Metal = SNIIMetal*Pstar[leaf]->InitialMass;
            FBInjection.ZstarII = Pstar[leaf]->ZII;
            FBInjection.ZstarIa = Pstar[leaf]->ZIa;
            FBInjection.Mode = 0;
        } else if(Pstar[leaf]->TypeIa == OFF){  // If SNIa case...
            double MassInMsun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS;
            FBInjection.Qheat = SNIIEnergy*MassInMsun;
            FBInjection.Mass = SNIIEjector*Pstar[leaf]->InitialMass;  
            FBInjection.Metal = SNIIMetal*Pstar[leaf]->InitialMass; 
            FBInjection.ZstarII = Pstar[leaf]->ZII;
            FBInjection.ZstarIa = Pstar[leaf]->ZIa;
            FBInjection.Mode = 1;
        }
        EnergyMassMetalInjectionEngine(FBInjection,0);
    }

    double TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.FeedbackCommThisStep += GetElapsedTime()-TimeComm;


    for(int i=0;i<NImportAll;i++) // Feedback from stars in external domains.
        EnergyMassMetalInjectionEngine(DelayedFBInjectionRecv[i],1);

    for(int i=0;i<NLocalTarget;i++){ 
        int leaf = DelayedFBActiveStars[i].Index;
        if(Pstar[leaf]->TypeII == OFF){ // If SNII case...
            PstarMass(leaf) -= SNIIEjector*Pstar[leaf]->InitialMass; 
            Pstar[leaf]->Mass -= SNIIEjector*Pstar[leaf]->InitialMass; 
            if( PstarMass(leaf)<=0.e0){
                fprintf(stderr,"to small mass: M %g IM %g SNIIE %g\n",
                        PstarMass(leaf),Pstar[leaf]->InitialMass,SNIIEjector);
                fprintf(stderr,"Mold %g IdM %g\n",
                        PstarMass(leaf)+SNIIEjector*Pstar[leaf]->InitialMass,
                        SNIIEjector*Pstar[leaf]->InitialMass);
            }
            assert(PstarMass(leaf)>0.e0);
            Pstar[leaf]->TypeII = ON;
        } else if(Pstar[leaf]->TypeIa == OFF){  // If SNIa case...
            PstarMass(leaf) -= SNIIEjector*Pstar[leaf]->InitialMass;  
            assert(PstarMass(leaf)>0.e0);
            Pstar[leaf]->TypeIa = ON;
        }
    }

    for(int i=0;i<NProcs-1;i++)
        free(DelayedFBInjectionSend[i]);
    free(DelayedFBInjectionRecv);

    return;
}


#else // (defined(PRESERVE_SNII_EVENTRATE) || defined(SET_SNII_TEMPERATURE))
static void EnergyMassMetalInjection(const int NLocalTarget, const int NProcs, bool DelayedFBExportFlags[][NProcs-1]){

    MPI_Status  mpi_status;

    // counter for export and import. 
    int BitMask = 0x01; 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    struct StructDelayedFBInjection *DelayedFBInjectionSend[NProcs],*DelayedFBInjectionRecv;

    for(int i=0;i<NProcs-1;i++){
        int NExport = 0;
        for(int k=0;k<NLocalTarget;k++){
            if(DelayedFBExportFlags[k][i]&BitMask)
                NExport ++;
        }
        NExportThisTime[i] = NExport;
        DelayedFBInjectionSend[i] = 
            malloc(sizeof(struct StructDelayedFBInjection)*(NExportThisTime[i]+1));

        NExport = 0;
        for(int k=0;k<NLocalTarget;k++){
            int leaf = DelayedFBActiveStars[k].Index;
            if(DelayedFBExportFlags[k][i]&BitMask){
                DelayedFBInjectionSend[i][NExport].Radius = DelayedFBActiveStars[k].Radius;
                DelayedFBInjectionSend[i][NExport].Rho = DelayedFBActiveStars[k].Density;
                DelayedFBInjectionSend[i][NExport].Pos[0] = PstarPosP(leaf)[0];
                DelayedFBInjectionSend[i][NExport].Pos[1] = PstarPosP(leaf)[1];
                DelayedFBInjectionSend[i][NExport].Pos[2] = PstarPosP(leaf)[2];
#ifdef MAXIMUM_ENERGY_INPUT
                DelayedFBInjectionSend[i][NExport].TargetID = DelayedFBActiveStars[k].DistanceMinGlobalID;
                assert(DelayedFBActiveStars[k].DistanceMinGlobalID != NONE);
#endif
                double MassInMsun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS;
                if(Pstar[leaf]->TypeII == OFF){ // If SNII case...
                    DelayedFBInjectionSend[i][NExport].Qheat = SNIIEnergy*MassInMsun;
                    DelayedFBInjectionSend[i][NExport].Mass = SNIIEjector*Pstar[leaf]->InitialMass; 
                    DelayedFBInjectionSend[i][NExport].Metal = SNIIMetal*Pstar[leaf]->InitialMass;
                    DelayedFBInjectionSend[i][NExport].ZstarII = Pstar[leaf]->ZII;
                    DelayedFBInjectionSend[i][NExport].ZstarIa = Pstar[leaf]->ZIa;
                    DelayedFBInjectionSend[i][NExport].Mode = 0;
                } else if(Pstar[leaf]->TypeIa == OFF){  // If SNIa case...
                    DelayedFBInjectionSend[i][NExport].Qheat = SNIIEnergy*MassInMsun; 
                    DelayedFBInjectionSend[i][NExport].Mass = SNIIEjector*Pstar[leaf]->InitialMass;  
                    DelayedFBInjectionSend[i][NExport].Metal = SNIIMetal*Pstar[leaf]->InitialMass; 
                    DelayedFBInjectionSend[i][NExport].ZstarII = Pstar[leaf]->ZII;
                    DelayedFBInjectionSend[i][NExport].ZstarIa = Pstar[leaf]->ZIa;
                    DelayedFBInjectionSend[i][NExport].Mode = 1;
                }
                NExport ++;
            }
        }
        assert(NExport == NExportThisTime[i]);
    }

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_FEEDBACK_DELAYED_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_FEEDBACK_DELAYED_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }

    DelayedFBInjectionRecv = malloc(sizeof(struct StructDelayedFBInjection)*(NImportAll+1));

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(DelayedFBInjectionSend[i],
            NExportThisTime[i]*sizeof(struct StructDelayedFBInjection),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_FEEDBACK_DELAYED_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(DelayedFBInjectionRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructDelayedFBInjection),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_FEEDBACK_DELAYED_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);

    // ok.
    for(int i=0;i<NLocalTarget;i++){ 
        int leaf = DelayedFBActiveStars[i].Index;
        struct StructDelayedFBInjection FBInjection;

        FBInjection.Pos[0] = PstarPosP(leaf)[0];
        FBInjection.Pos[1] = PstarPosP(leaf)[1];
        FBInjection.Pos[2] = PstarPosP(leaf)[2];
        FBInjection.Rho = DelayedFBActiveStars[i].Density;
        FBInjection.Radius = DelayedFBActiveStars[i].Radius;
#ifdef MAXIMUM_ENERGY_INPUT
        FBInjection.TargetID = DelayedFBActiveStars[i].DistanceMinGlobalID;
#endif
        double MassInMsun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS;
        if(Pstar[leaf]->TypeII == OFF){ // If SNII case...
            FBInjection.Qheat = SNIIEnergy*MassInMsun;
            FBInjection.Mass = SNIIEjector*Pstar[leaf]->InitialMass; 
            FBInjection.Metal = SNIIMetal*Pstar[leaf]->InitialMass;
            FBInjection.ZstarII = Pstar[leaf]->ZII;
            FBInjection.ZstarIa = Pstar[leaf]->ZIa;
            FBInjection.Mode = 0;
            /*
            fprintf(stderr,"N = %g, Q, Mass, Metal = %g %g %g\n",
                MassInMsun,FBInjection.Qheat,FBInjection.Mass,FBInjection.Metal);
            fprintf(stderr,"Q, Mass, Metal = %g ergs, %g Msun, %g Msun\n",
                FBInjection.Qheat/GetUnitEnergy(),FBInjection.Mass*Pall.UnitMass/MSUN_CGS,
                    FBInjection.Metal*Pall.UnitMass/MSUN_CGS);
            */
        } else if(Pstar[leaf]->TypeIa == OFF){  // If SNIa case...
            FBInjection.Qheat = SNIIEnergy*MassInMsun;
            FBInjection.Mass = SNIIEjector*Pstar[leaf]->InitialMass;  
            FBInjection.Metal = SNIIMetal*Pstar[leaf]->InitialMass; 
            FBInjection.ZstarII = Pstar[leaf]->ZII;
            FBInjection.ZstarIa = Pstar[leaf]->ZIa;
            FBInjection.Mode = 1;
        }
        EnergyMassMetalInjectionEngine(FBInjection,0);
    }

    double TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.FeedbackCommThisStep += GetElapsedTime()-TimeComm;


    for(int i=0;i<NImportAll;i++)
        EnergyMassMetalInjectionEngine(DelayedFBInjectionRecv[i],1);

    for(int i=0;i<NLocalTarget;i++){ 
        int leaf = DelayedFBActiveStars[i].Index;
        if(Pstar[leaf]->TypeII == OFF){ // If SNII case...
            PstarMass(leaf) -= SNIIEjector*Pstar[leaf]->InitialMass; 
            Pstar[leaf]->Mass -= SNIIEjector*Pstar[leaf]->InitialMass; 
            if( PstarMass(leaf)<=0.e0){
                fprintf(stderr,"to small mass: M %g IM %g SNIIE %g\n",
                        PstarMass(leaf),Pstar[leaf]->InitialMass,SNIIEjector);
                fprintf(stderr,"Mold %g IdM %g\n",
                        PstarMass(leaf)+SNIIEjector*Pstar[leaf]->InitialMass,
                        SNIIEjector*Pstar[leaf]->InitialMass);
            }
            assert(PstarMass(leaf)>0.e0);
            Pstar[leaf]->TypeII = ON;
        } else if(Pstar[leaf]->TypeIa == OFF){  // If SNIa case...
            PstarMass(leaf) -= SNIIEjector*Pstar[leaf]->InitialMass;  
            assert(PstarMass(leaf)>0.e0);
            Pstar[leaf]->TypeIa = ON;
        }
    }

    for(int i=0;i<NProcs-1;i++)
        free(DelayedFBInjectionSend[i]);
    free(DelayedFBInjectionRecv);

    return;
}
#endif // (defined(PRESERVE_SNII_EVENTRATE) || defined(SET_SNII_TEMPERATURE))


static int DelayedFBExportFlagsMaxAllocated = 0;
void DelayedSNe(void){

    double TimingResultThisRoutine = GetElapsedTime();

#ifdef DELAYED_FEEDBACK
    int NProcs = MPIGetNumProcs();

    int NLocalTarget = GetSNeNumber();
    int NGlobalTarget;
    MPI_Allreduce(&NLocalTarget,&NGlobalTarget,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(NGlobalTarget==0)
        return;

    static bool (*DelayedFBExportFlags)[NProcs-1];
    if(DelayedFBExportFlagsMaxAllocated < NLocalTarget){
        if(DelayedFBExportFlagsMaxAllocated > 0)
            free(DelayedFBExportFlags);
        DelayedFBExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*NLocalTarget,NAdditionUnit));
        DelayedFBExportFlags = malloc(sizeof(bool)*DelayedFBExportFlagsMaxAllocated*(NProcs-1));
    }


    if(Pall.NActivesHydro_t==0){
        TimingResults.HydroTreeThisStep = GetElapsedTime();
        PlantHydroTreeUpdate();
        TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;
    }

    // get nearest hydro particle's kernel size
    GetFeedbackRadius(NLocalTarget,NProcs,DelayedFBExportFlags);

    // clear dMass;
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->dZII = 0.e0;
        Phydro[i]->dZIa = 0.e0;
        Phydro[i]->dMass = 0.e0;
    }

    // add mass, metal, and energy.
    EnergyMassMetalInjection(NLocalTarget,NProcs,DelayedFBExportFlags);

    // add dMass;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->dZII > 0.e0)
            Phydro[i]->ZII = (Phydro[i]->ZII*PhydroMass(i)+Phydro[i]->dZII)/(PhydroMass(i)+Phydro[i]->dMass);
        if(Phydro[i]->dZIa > 0.e0)
            Phydro[i]->ZIa = (Phydro[i]->ZIa*PhydroMass(i)+Phydro[i]->dZIa)/(PhydroMass(i)+Phydro[i]->dMass);
        if(Phydro[i]->dMass > 0.e0){
            Phydro[i]->Mass += Phydro[i]->dMass;
            PhydroBody(i)->Mass = Phydro[i]->Mass;
        }
        Phydro[i]->Z = Phydro[i]->ZII;
    }
#ifdef UPDATE_TIMESTEP_IN_HEATEDREGION
    const static double dtFact = TFactorCourant*2.0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(!Phydro[i]->Active){
        if(Phydro[i]->DQheat > 0.e0){
	        double U_heat = Phydro[i]->DQheat/PhydroMass(i);
            double cs_heat = sqrt(Pall.GGm1*U_heat);
            double dt_heat = dtFact*Phydro[i]->KernelPred/cs_heat;
            int k_heat = (int)(log2(dt_heat/Pall.dtmin));
            Phydro[i]->k_hydro_localmin = MIN(k_heat-MAX_K_LOCAL,Phydro[i]->k_hydro_localmin);
            // fprintf(stderr,"[%02d] dt [%g->%g] shrinks takes place due to heating for particle %ld\n",
                    // MPIGetMyID(),Phydro[i]->dt_hydro,dt_heat,PhydroBody(i)->GlobalID);

            dt_heat = Pall.dtmin*exp2((double)k_heat);
            int step = (int)(Pall.EraLocal/dt_heat);
            double NextUpdateEra;
            do{
                NextUpdateEra = step*dt_heat;
                step ++;
            }while(NextUpdateEra <= Pall.EraLocal);
            Phydro[i]->NextUpdateEra = NextUpdateEra;
        }
        }
    }
#endif 

#endif 

    TimingResults.FeedbackThisStep = GetElapsedTime()-TimingResultThisRoutine;

    return;
}

