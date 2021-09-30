#include "config.h"
#include "NeighborSearch.h"
#include "HydroAcc.h"
#include "PlantStellarTree.h"
#include "SizeDetermination.h"


#define ShowLog OFF

/*
 * This file involves functions concerning sink particles. There are two modes.
 * The first one only consider the sinking length. In this mode, all particles
 * in the sinking length are sinked by the sink particles.  Adaptive and
 * non-adaptive sinking lenghts are available. The second one consider both the
 * sinking length and the binding condition. The scan of the binding condition
 * is only done in the case within the sinking length. This mode also adopts
 * adaptive and non-adaptive sinking lengths.
 */

// Proto-type declarations.
static void MakeSinkParticles(void);
static void SinkHydroAccretion(void);
static void SinkSinkMerging(void);
static bool CheckSinkConditions(const int index);
static void ReConnectHydroAndSinkPointers(void);
static void CheckPermissionForMakeSinkParticle(void);
static void CheckSinkFormationCondition(void);
static void CalcSinkPressureForces(void);
// Proto-type declarations.
static int NumberofAbsorbedParticles = 0;
static bool *AllowSinkConversion; // This is the acceptance/regection flag of the gas->sink conversion.

struct StructSinkExport{ 
    // This structure includes information relating to a sink particle.
    double  AccretionRadius;  // Sink radius.
    double  Pos[3];  // Positions.
    double  Vel[3];  // Volocities.
    double  Mass;    // Mass.

#ifdef USE_CELIB
    double  Elements[CELibYield_Number]; //H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, Ni
#endif // USE_CELIB

    double  dt_localmin; // the local minimum dt.
    double  Eps; // the softening length.
    unsigned long int GlobalID; // GlobalID
    int     Index;   // Index of the sink particle.
};

struct StructSinkCheckExport{ 
    // This structure includes information relating to a sink particle.
    double  AccretionRadius;  // Sink radius.
    double  Pos[3];  // Positions.
};

struct StructSinkImport{ 
    // This structure includes information relating to a victime who will be sank.
    double  Pos[3];  // Positions.
    double  Vel[3];  // Velocities.
    double  Mass;    // Mass.

#ifdef USE_CELIB
    double  Elements[CELibYield_Number]; //H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, Ni
#endif // USE_CELIB

    int     TargetSinkIndex; // Index of the target sink particle.
    int     Type;   // Type of the victime (Gas/Star).
};

struct StructHydroSinkFlags{
    int TargetSinkID;   // The target sink particle ID.
    int NodeID;         // The node ID of the target sink particle.
    double Distance;    // Distance to the target sink particle.
};

struct StructSinkSinkFlags{
    int TargetSinkID;   // The target sink particle ID.
    int NodeID;         // The node ID of the target sink particle.
    double Distance;    // Distance to the target sink particle.
};

struct StructSinkSinkFlags2{
    int TargetSinkID;   // The target sink particle ID.
    int NodeID;         // The node ID of the target sink particle.
    double TotalEnergy; // The total energy of pair particles.
};

/*
 * In this function, the parameters for sink particle operation are converted
 * into the simulation unit.
 */
static double SinkGasAbsorptionThresholdDensity;
static void InitSinkParticles(void){

    Pall.SinkThresholdDensity = SINK_THRESHOLD_DENSITY/Pall.ConvertDensityToCGS;
    if(Pall.Nhydro>0){
        AllowSinkConversion = malloc(sizeof(bool)*MAX(Pall.Nhydro,NAdditionUnit));
        NumberofAbsorbedParticles = MAX(Pall.Nhydro,NAdditionUnit);
    }
    SinkGasAbsorptionThresholdDensity = SINK_GASABSORPTION_THRESHOLD_DENSITY/Pall.ConvertDensityToCGS;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The threshold density for sink particle formation is %g g/cm^3 and %g in simulation unit.\n",
                SINK_THRESHOLD_DENSITY,Pall.SinkThresholdDensity);
        //fprintf(stderr,"%g %g \n",SINK_THRESHOLD_DENSITY, SINK_THRESHOLD_DENSITY*CUBE(Pall.UnitLength)/Pall.UnitMass);
        //fflush(NULL);
        //exit(1);
    }
    
    return;
}

/*
 * This inline function returns the distance from the particle to the external
 * domain.
 */
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

/*
 * This inline function returns the judgement of the overlap.
 */
static inline bool CheckDomainOverlapForActiveHydroEdges(double Pos[restrict], const double h, const int NodeID) __attribute__((always_inline));
static inline bool CheckDomainOverlapForActiveHydroEdges(double Pos[restrict], const double h, const int DomainID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForActiveHydro[DomainID].PosMin[k]) 
            Dist2 += SQ(EdgesForActiveHydro[DomainID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForActiveHydro[DomainID].PosMax[k])
            Dist2 += SQ(EdgesForActiveHydro[DomainID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

/*
 * This inline function returns the judged result of the overlap for the
 * "Active Sink Edges ".
 */
static inline bool CheckDomainOverlapForActiveSinkEdges(double Pos[restrict], const double h, const int NodeID) __attribute__((always_inline));
static inline bool CheckDomainOverlapForActiveSinkEdges(double Pos[restrict], const double h, const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForActiveSink[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForActiveSink[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForActiveSink[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForActiveSink[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

/*
 * This function checks and returns number of sink particles which have to
 * export. The current implementation is quite simple. It may be desired for us
 * to make further efficient check algorithm.
 */
static inline int CheckSinkExportFlags(const int Index, const int NProcs, bool SinkExportFlags[][NProcs-1]) __attribute__((always_inline));
static inline int CheckSinkExportFlags(const int Index, const int NProcs, bool SinkExportFlags[][NProcs-1]){

    ///// Check must be necessary.

    int DomainID = CommunicationTable[Index].SendRank;
    int NExport = 0;
    for(int i=0;i<Pall.Nsink;i++){
        if(PsinkBody(i)->Active){
            if(CheckDomainOverlapForActiveHydroEdges(Psink[i]->PosP,Psink[i]->AccretionRadius,DomainID)){
                SinkExportFlags[i][Index] = true;
                NExport ++;
            }
        }
    }

    return NExport;
}

/*
 * This function checks and returns number of sink particles which have to
 * export. The check is done for the "Active Sink Edges".  The current
 * implementation is quite simple. It may be desired for us to make further
 * efficient check algorithm.
 */
static inline int CheckSinkExportFlagsForActiveEdges(const int Index, const int NProcs, bool SinkExportFlags[][NProcs-1]) __attribute__((always_inline));
static inline int CheckSinkExportFlagsForActiveEdges(const int Index, const int NProcs, bool SinkExportFlags[][NProcs-1]){

    int NodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;

    for(int i=0;i<Pall.Nsink;i++){
        if(PsinkBody(i)->Active){
            //if(OverlapDomainDensity(Psink[i]->PosP,Psink[i]->Radius,NodeID)){
            if(CheckDomainOverlapForActiveHydroEdges(Psink[i]->PosP,Psink[i]->AccretionRadius,NodeID)){
                SinkExportFlags[i][Index] = ON;
                NExport ++;
            }
        }
    }

    return NExport;
}

/*
 * This function checks the bound condistion between two particles (sink and
 * hydro particles). If the bound condition is satisfied, the fucntion returns
 * `true', otherwise the function returns `false'.
 */
static inline bool __attribute__((always_inline)) SinkHydroBoundCheck(const double Mass_sink, const double Pos_sink[], const double Vel_sink[],
        const double Mass_hydro, const double Pos_hydro[], const double Vel_hydro[], const double U){

    double distance = DISTANCE(Pos_sink,Pos_hydro);
    double velocity2 = DISTANCE2(Vel_sink,Vel_hydro);
    double mu = Mass_sink*Mass_hydro/(Mass_sink+Mass_hydro);

    //double E = 0.5*mu*velocity2 + U*Mass_hydro - Pall.GravConst*SQ(mu)/distance;
    //double E = 0.5*mu*velocity2 - Pall.GravConst*SQ(mu)/distance;
    double E = 0.5*mu*velocity2 - Pall.GravConst*Mass_sink*Mass_hydro/distance;

    //return true ;
    if(E < 0){
        //fprintf(stderr,"Bound check clear %g %g %g %g\n",E,0.5*mu*velocity2,U*Mass_hydro,-Pall.GravConst*SQ(mu)/distance);
        return true;
    } else {
        //fprintf(stderr,"Bound check failure %g %g %g %g\n",E,0.5*mu*velocity2,U*Mass_hydro,-Pall.GravConst*SQ(mu)/distance);
        return false;
    }
}

/*
 * This function checks the bound condistion between two particles (sink and
 * star particles). If the bound condition is satisfied, the fucntion returns
 * `true', otherwise the function returns `false'.
 */
static inline bool __attribute__((always_inline)) SinkStarBoundCheck(const double Mass_sink, const double Pos_sink[], const double Vel_sink[],
        const double Mass_star, const double Pos_star[], const double Vel_star[]){

    double distance = DISTANCE(Pos_sink,Pos_star);
    double velocity2 = DISTANCE2(Vel_sink,Vel_star);
    double mu = Mass_sink*Mass_star/(Mass_sink+Mass_star);

    double E = 0.5*mu*velocity2 - Pall.GravConst*Mass_sink*Mass_star/distance;

    //return true ;
    if(E < 0){
        return true;
    } else {
        return false;
    }
}

static inline bool SinkHydroAngularMomentumCheck(const double AccretionRadius, const double Mass_sink, const double Pos_sink[], const double Vel_sink[], const double Eps_sink, const double Mass_hydro, const double Pos_hydro[], const double Vel_hydro[], const double Eps_hydro) __attribute__((always_inline));
static inline bool SinkHydroAngularMomentumCheck(const double AccretionRadius, 
        const double Mass_sink, const double Pos_sink[], const double Vel_sink[], const double Eps_sink,
        const double Mass_hydro, const double Pos_hydro[], const double Vel_hydro[], const double Eps_hydro){

    /*
     * Lacc = r_acc*v(r_acc) 
     *      = r_acc*sqrt(G*M/r), where r^2 = r_acc^2 + Eps2;
     *      = sqrt(r_acc^2*G*M/sqrt(r^2)
     */
    double Racc2 = SQ(AccretionRadius);
#ifdef USE_SYMMETRIZED_SOFTENING
    double Lacc = sqrt(Racc2*Pall.GravConst*Mass_sink/sqrt(Racc2+SQ(Eps_sink)+SQ(Eps_hydro)));
#else  // USE_SYMMETRIZED_SOFTENING
    double Lacc = sqrt(Racc2*Pall.GravConst*Mass_sink/sqrt(Racc2+SQ(Eps_sink)));
#endif // USE_SYMMETRIZED_SOFTENING
    
    double l[3] = {(Pos_hydro[1]-Pos_sink[1])*(Vel_hydro[2]-Vel_sink[2]),
                   (Pos_hydro[2]-Pos_sink[2])*(Vel_hydro[0]-Vel_sink[0]),
                   (Pos_hydro[0]-Pos_sink[0])*(Vel_hydro[1]-Vel_sink[1])};
    double Lhydro = NORM(l);
    
    if(Lhydro < Lacc){
        return true;
    } else {
        return false;
    }
}

static inline bool SinkStarAngularMomentumCheck(const double AccretionRadius, const double Mass_sink, const double Pos_sink[], const double Vel_sink[], const double Eps_sink, const double Mass_star, const double Pos_star[], const double Vel_star[], const double Eps_star) __attribute__((always_inline));
static inline bool SinkStarAngularMomentumCheck(const double AccretionRadius, 
        const double Mass_sink, const double Pos_sink[], const double Vel_sink[], const double Eps_sink,
        const double Mass_star, const double Pos_star[], const double Vel_star[], const double Eps_star){

    /*
     * Lacc = r_acc*v(r_acc) 
     *      = r_acc*sqrt(G*M/r), where r^2 = r_acc^2 + Eps2;
     *      = sqrt(r_acc^2*G*M/sqrt(r^2)
     */
    double Racc2 = SQ(AccretionRadius);
#ifdef USE_SYMMETRIZED_SOFTENING
    double Lacc = sqrt(Racc2*Pall.GravConst*Mass_sink/sqrt(Racc2+SQ(Eps_sink)+SQ(Eps_star)));
#else  // USE_SYMMETRIZED_SOFTENING
    double Lacc = sqrt(Racc2*Pall.GravConst*Mass_sink/sqrt(Racc2+SQ(Eps_sink)));
#endif // USE_SYMMETRIZED_SOFTENING
    
    double l[3] = {(Pos_star[1]-Pos_sink[1])*(Vel_star[2]-Vel_sink[2]),
                   (Pos_star[2]-Pos_sink[2])*(Vel_star[0]-Vel_sink[0]),
                   (Pos_star[0]-Pos_sink[0])*(Vel_star[1]-Vel_sink[1])};
    double Lhydro = NORM(l);
    
    if(Lhydro < Lacc){
        return true;
    } else {
        return false;
    }
}

#ifdef USE_SINK_HILL_CONDITION //{
#define X_HILL (4.0)
static bool CheckSinkHillCondition(const double Rho, const double Pos[restrict], const double Acc[restrict]){

    for(int i=0;i<Pall.Nsink;i++){
        double dx[] = {Pos[0]-PsinkBody(i)->PosP[0],
                       Pos[1]-PsinkBody(i)->PosP[1],
                       Pos[2]-PsinkBody(i)->PosP[2]};
        double da[] = {Acc[0]-PsinkBody(i)->Acc[0],
                       Acc[1]-PsinkBody(i)->Acc[1],
                       Acc[2]-PsinkBody(i)->Acc[2]};

        double Rho_Hill = 3.0*X_HILL*(-DOT_PRODUCT(dx,da))/(4.0*M_PI*Pall.GravConst*NORM2(dx));

        if(Rho < Rho_Hill){
            return false;
        }
    }

    return true;
}
#endif // USE_SINK_NOLOCALSINK_CONDITION //}

static bool __attribute__((always_inline)) SinkSinkBoundCheck(const double Mass_sink1, const double Pos_sink1[], const double Vel_sink1[],
        const double Mass_sink2, const double Pos_sink2[], const double Vel_sink2[]){

    double distance = DISTANCE(Pos_sink1,Pos_sink2);
    double velocity2 = DISTANCE2(Vel_sink1,Vel_sink2);
    double mu = Mass_sink1*Mass_sink2/(Mass_sink1+Mass_sink2);

    double E = 0.5*mu*velocity2 - Pall.GravConst*SQ(mu)/distance;
    if(E < 0){
        fprintf(stderr,"Bound check clear %g %g %g\n",E,0.5*mu*velocity2,-Pall.GravConst*SQ(mu)/distance);
        return true;
    } else {
        fprintf(stderr,"Bound check failure %g %g %g\n",E,0.5*mu*velocity2,-Pall.GravConst*SQ(mu)/distance);
        return false;
    }

    // return true ;
    /*
    if(E < 0){
        fprintf(stderr,"Bound check clear %g %g %g\n",E,0.5*mu*velocity2,-Pall.GravConst*SQ(mu)/distance);
        return true;
    } else {
        fprintf(stderr,"Bound check failure %g %g %g\n",E,0.5*mu*velocity2,-Pall.GravConst*SQ(mu)/distance);
        return false;
    }
    */
}

/*
 * Operations relating particle absorption are described in this function.
 */
static inline void __attribute__((always_inline)) ParticleAbsorption(const int leaf_sink, const double Pos[restrict], const double Vel[restrict], const double Mass, double Elements[], const int Type){

    // Calculate center of mass (COM) and velocities of COM.
    double COM[3] = {PsinkMass(leaf_sink)*PsinkPos(leaf_sink)[0] + Mass*Pos[0],
                     PsinkMass(leaf_sink)*PsinkPos(leaf_sink)[1] + Mass*Pos[1],
                     PsinkMass(leaf_sink)*PsinkPos(leaf_sink)[2] + Mass*Pos[2]};

    double dt_half = 0.5*PsinkBody(leaf_sink)->dt;
    double VCOM[3] = {
        PsinkMass(leaf_sink)*PsinkVel(leaf_sink)[0] + Mass*Vel[0],
        PsinkMass(leaf_sink)*PsinkVel(leaf_sink)[1] + Mass*Vel[1],
        PsinkMass(leaf_sink)*PsinkVel(leaf_sink)[2] + Mass*Vel[2]};
    double mu = PsinkMass(leaf_sink)+Mass;
    double imu = 1.e0/mu;
    for(int k=0;k<3;k++){
        COM[k] *= imu;
        VCOM[k] *= imu;
    }

    // Update values of the sink particle.
    for(int k=0;k<3;k++){
        PsinkBody(leaf_sink)->Pos[k] = PsinkBody(leaf_sink)->PosP[k] = 
        Psink[leaf_sink]->PosP[k] = COM[k];
        PsinkBody(leaf_sink)->Vel[k] = 
        Psink[leaf_sink]->VelP[k] = VCOM[k];
    }
    PsinkMass(leaf_sink) = mu;
    Psink[leaf_sink]->NumberofAbsorbedParticles ++;

    if(Type == TypeHydro){
        Psink[leaf_sink]->AccretionMassGas += Mass;
    } else if(Type == TypeStar){
        Psink[leaf_sink]->AccretionMassStar += Mass;
    }
    Psink[leaf_sink]->AccretionMass += Mass;

    // Update metallicity
    
#ifdef USE_CELIB //{
    if(Type == TypeHydro){
        for(int k=0;k<CELibYield_Number;k++){
            Psink[leaf_sink]->Elements[k] += Elements[k];
        }
    } else if(Type == TypeStar){
        for(int k=0;k<CELibYield_Number;k++){
            Psink[leaf_sink]->Elements[k] += Elements[k];
        }
    }
#else

#endif // USE_CELIB //}

    return ;
}


FILE *fp_sinklog;
static bool LogSinkParticlesFirstFlag = true;
static void LogSinkParticles(void){

    if(LogSinkParticlesFirstFlag == true){
        if(MPIGetMyID() == MPI_ROOT_RANK)
            FileOpen(fp_sinklog,"./SinkLog.data","w");
        LogSinkParticlesFirstFlag = false;
    }
    double mass[2];
    mass[0] = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        mass[0] += Phydro[i]->Mass;
    }
    mass[1] = 0.e0;
    for(int i=0;i<Pall.Nsink;i++){
        mass[1] += PsinkMass(i);
    }
    double GlobalMass[2];
    MPI_Allreduce(mass,GlobalMass,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // allreduce masses.
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(fp_sinklog,"%g %ld %ld %ld %g %g %g\n",
            Pall.TCurrent,Pall.Nhydro_t,Pall.Nsink_t,Pall.Nhydro_t+Pall.Nsink_t,
            GlobalMass[0],GlobalMass[1],GlobalMass[0]+GlobalMass[1]);

    return ;
}

/*
 * This function generates new sink particles. We do not care the internal
 * energy of hydro particles. 
 */
static void MakeSinkParticlesNew(void){

    int NewSinkParticles = 0; 
    for(int i=0;i<Pall.Nhydro;i++){ // Search all hydro particles.
        if(AllowSinkConversion[i] == false) continue;

        StructPsinkptr Psk = ReturnEmptySinkStructurePointer();

        // Setup sink partcile.
        Psk->Use = ON;
        Psk->ParentGlobalID = PhydroBody(i)->GlobalID;
        Psk->Body = PhydroBody(i);
        Psk->Body->Baryon = (void*)Psk;
        Psk->Body->Type = TypeSink;
        Psk->FormationTime = Pall.TCurrent;
        Psk->dt_localmin = PhydroBody(i)->dt;

        Psk->PosP[0] = Phydro[i]->PosP[0];
        Psk->PosP[1] = Phydro[i]->PosP[1];
        Psk->PosP[2] = Phydro[i]->PosP[2];

        Psk->VelP[0] = Phydro[i]->VelP[0];
        Psk->VelP[1] = Phydro[i]->VelP[1];
        Psk->VelP[2] = Phydro[i]->VelP[2];

        // follow Z
        Psk->Z   = Phydro[i]->Z;
        Psk->ZII = Phydro[i]->ZII;
        Psk->ZIa = Phydro[i]->ZIa;

#ifdef USE_CELIB //{
        for(int k=0;k<CELibYield_Number;k++){
            Psk->Elements[k] = Phydro[i]->Elements[k];
        }
#endif // USE_CELIB //}

        Psk->AccretionRadius = SINKHYDRO_ACCRETION_RADIUS/Pall.UnitLength;
        Psk->MergingDistance = SINKSINK_MERGING_DISTANCE/Pall.UnitLength;

        HydroRoot.Leaves[Phydro[i]->Leaf] *= -1;

        //Remove the hydro particle.
        Phydro[i]->Use = OFF;

        NewSinkParticles ++;
    }

    int GlobalNewSinkParticles;
    MPI_Allreduce(&NewSinkParticles,&GlobalNewSinkParticles,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(GlobalNewSinkParticles > 0){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            if(GlobalNewSinkParticles == 1){
                fprintf(stderr,"<')_+++< %d new sink particle is born in this step\n",GlobalNewSinkParticles);
            } else {
                fprintf(stderr,"<')_+++< %d new sink particles are born in this step\n",GlobalNewSinkParticles);
            }
        }
    }


    /*
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    */

    if(GlobalNewSinkParticles > 0){
        Pall.Nhydro -= NewSinkParticles;
        Pall.Nhydro_t -= GlobalNewSinkParticles;
        Pall.Nsink += NewSinkParticles;
        Pall.Nsink_t += GlobalNewSinkParticles;

        ReConnectPointers(); // Is this necessary?
        //ReConnectPointers();
        UpdateTotalNumber();
        //UpdateTotalActiveNumber();

        for(int i=0;i<HydroRoot.NumberofLeaves;i++)
            HydroRoot.Leaves[i] = NONE;

        for(int i=0;i<Pall.Nhydro;i++){
            int index = Phydro[i]->Leaf;
            NBCache[index].Leaf = i;
            HydroRoot.Leaves[index] = i;
        }
        //UpdateTotalActiveNumber();
    }

    /*
    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"New/GNew %d, %d, Nhydro/Nsink %ld, %ld\n",
            NewSinkParticles,GlobalNewSinkParticles,Pall.Nhydro,Pall.Nsink);
        fprintf(stderr,"Pall.Ntotal_t,Nhydro_t,Nsink_t = %ld, %ld, %ld\n",
            Pall.Ntotal_t,Pall.Nhydro_t,Pall.Nsink_t);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    */

    return ;
}

static void UpdateSinkEdges(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    struct StructEdges TempEdges[2];
    int RootNodeID = 0;

    if(EdgesForSink == NULL)        
        EdgesForSink = malloc(sizeof(struct StructEdges)*NProcs);
    if(EdgesForActiveSink == NULL)  
        EdgesForActiveSink = malloc(sizeof(struct StructEdges)*NProcs);


    if(Pall.Nsink == 0){
        double PosLeaf[3] = {HydroNode[RootNodeID].Pos[0],HydroNode[RootNodeID].Pos[1],HydroNode[RootNodeID].Pos[2]};
        for(int i=0;i<2;i++){
            TempEdges[i].PosMax[0] = PosLeaf[0]; TempEdges[i].PosMin[0] = PosLeaf[0];
            TempEdges[i].PosMax[1] = PosLeaf[1]; TempEdges[i].PosMin[1] = PosLeaf[1];
            TempEdges[i].PosMax[2] = PosLeaf[2]; TempEdges[i].PosMin[2] = PosLeaf[2];
        }
    } else {
        double max[3],min[3];
        for(int k=0;k<3;k++){
            max[k] = Psink[0]->PosP[k];
            min[k] = Psink[0]->PosP[k];
        }
        for(int i=1;i<Pall.Nsink;i++){
            for(int k=0;k<3;k++){
                max[k] = fmax(Psink[i]->PosP[k],max[k]);
                min[k] = fmin(Psink[i]->PosP[k],min[k]);
            }
        }
        TempEdges[0].PosMax[0] = max[0]; TempEdges[0].PosMin[0] = min[0];
        TempEdges[0].PosMax[1] = max[1]; TempEdges[0].PosMin[1] = min[1];
        TempEdges[0].PosMax[2] = max[2]; TempEdges[0].PosMin[2] = min[2];

        // search for the first active hydro particle
        if(Pall.NActivesSink > 0){
            int first = true;
            for(int i=0;i<Pall.Nsink;i++){
                if(PsinkActive(i)){
                    if(first == true){
                        for(int k=0;k<3;k++){
                            max[k] = Psink[i]->PosP[k];
                            min[k] = Psink[i]->PosP[k];
                        }
                        first = false; 
                    } else {
                        for(int k=0;k<3;k++){
                            max[k] = fmax(Psink[i]->PosP[k],max[k]);
                            min[k] = fmin(Psink[i]->PosP[k],min[k]);
                        }
                    }
                }
            }
        } else {
            for(int k=0;k<3;k++)
                max[k] = min[k] = HydroNode[RootNodeID].Pos[k];
        }

        TempEdges[1].PosMax[0] = max[0]; TempEdges[1].PosMin[0] = min[0];
        TempEdges[1].PosMax[1] = max[1]; TempEdges[1].PosMin[1] = min[1];
        TempEdges[1].PosMax[2] = max[2]; TempEdges[1].PosMin[2] = min[2];
    }

    MPI_Status  mpi_status;

    EdgesForSink[MyID] = TempEdges[0];
    EdgesForActiveSink[MyID] = TempEdges[1];

    for(int i=0;i<NProcs-1;i++){
        struct StructEdges RecvEdges[2];
        MPI_Sendrecv(TempEdges,2*sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_PLANTTREE_EXTENSITY,
                     RecvEdges,2*sizeof(struct StructEdges),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_PLANTTREE_EXTENSITY,
                     MPI_COMM_WORLD,&mpi_status);
        EdgesForSink[CommunicationTable[i].RecvRank] = RecvEdges[0];
        EdgesForActiveSink[CommunicationTable[i].RecvRank] = RecvEdges[1];
    }

#if 0
    for(int i=0;i<NProcs;i++){
        fprintf(stderr,"[%02d->%02d] |%d| %g %g %g %g %g %g\n",
                MPIGetMyID(),CommunicationTable[i].SendRank,i,
                EdgesForActiveSink[i].PosMax[0],EdgesForActiveSink[i].PosMax[1],EdgesForActiveSink[i].PosMax[2],
                EdgesForActiveSink[i].PosMin[0],EdgesForActiveSink[i].PosMin[1],EdgesForActiveSink[i].PosMin[2]);
        fprintf(stderr,"[%02d->%02d] |%d| %g %g %g %g %g %g\n",
                MPIGetMyID(),CommunicationTable[i].SendRank,i,
                EdgesForSink[i].PosMax[0],EdgesForSink[i].PosMax[1],EdgesForSink[i].PosMax[2],
                EdgesForSink[i].PosMin[0],EdgesForSink[i].PosMin[1],EdgesForSink[i].PosMin[2]);
    }
#endif

    return ;
}

/*
 * This function handles the sink particles. 
 */
static bool InitSinkParticlesFlag = false;
void SinkParticles(void){

#if 0
    // 
    int counter = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->RhoPred > Pall.SinkThresholdDensity){
            counter ++;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE,&counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(MPIGetMyID() == 0)
        fprintf(stderr,"Target count: %d %d \n",counter,Pall.Nhydro_t);
    fflush(NULL);

#endif


#ifndef USE_SINK_PARTICLE
    return ;
#endif
    if(InitSinkParticlesFlag == false){
        InitSinkParticles();
        //fprintf(stderr,"SINK_RADIUS = %g\n",SINK_RADIUS);
        InitSinkParticlesFlag = true;
    }
    if(Pall.Nhydro>0){
        if(Pall.Nhydro > NumberofAbsorbedParticles){
            NumberofAbsorbedParticles = ForAngelsShare*MAX(Pall.Nhydro,NAdditionUnit);
            AllowSinkConversion = realloc(AllowSinkConversion,sizeof(bool)*NumberofAbsorbedParticles);
        }
    }

    for(int i=0;i<Pall.Nhydro;i++){
        AllowSinkConversion[i] = true;
    }


    // update active sink edges.
    UpdateSinkEdges();
    if(Pall.Nsink_t > 0){ // Pall.NActiveSink_t;
        SinkHydroAccretion();
    }
    SinkSinkMerging();


    // Make new function
#ifdef MAKE_SINK_PARTICLE
    //CheckPermissionForMakeSinkParticle();
    // MakeSinkParticles();
    CheckSinkFormationCondition();
    MakeSinkParticlesNew();
#endif

    return ;
}


/*
 * This function checks and returns number of sink particles which should be
 * exported. This function is only for the acception/regection check of sink
 * formation for dense gas particle.
 */
static inline int CheckSinkCheckExportFlags(const int Index, const int NProcs, bool SinkCheckExportFlags[][NProcs-1]) __attribute__((always_inline));
static inline int CheckSinkCheckExportFlags(const int Index, const int NProcs, bool SinkCheckExportFlags[][NProcs-1]){

    int NodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;
    for(int i=0;i<Pall.Nsink;i++){
        SinkCheckExportFlags[i][Index] = OFF;
        if(CheckDomainOverlapForActiveHydroEdges(Psink[i]->PosP,Psink[i]->AccretionRadius,NodeID)){
            SinkCheckExportFlags[i][Index] = ON;
            NExport ++;
        }
    }

    return NExport;
}


static void CheckSinkFormationHillCondition(void){

    int NProcs = MPIGetNumProcs();

    struct StructSinkForHill{
        double Pos[3];
        double Acc[3];
    };
    //SinkForHill[Pall.Nsink_t];

    // counter for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    int SendFlag,RecvFlag;

    struct StructSinkForHill *SinkForHillExportSend[NProcs];

    struct StructSinkForHill *SinkForHillExportRecv;
    CheckSizeofBufferExportRecv(Pall.Nsink_t,sizeof(struct StructSinkForHill));
    SinkForHillExportRecv = BufferExportRecv;

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(NImportThisTime+i,1,MPI_INT,
            CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);

        NExportThisTime[i] = Pall.Nsink;
        MPI_Isend(NExportThisTime+i,1,MPI_INT,
                CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);

        CheckSizeofBufferExportSendIndex(NExportThisTime[i],
                sizeof(struct StructSinkForHill),i);
        SinkForHillExportSend[i] = BufferExportSend[i];

        for(int k=0;k<Pall.Nsink;k++){
            SinkForHillExportSend[i][k].Pos[0] = PsinkBody(k)->Pos[0];
            SinkForHillExportSend[i][k].Pos[1] = PsinkBody(k)->Pos[1];
            SinkForHillExportSend[i][k].Pos[2] = PsinkBody(k)->Pos[2];
            SinkForHillExportSend[i][k].Acc[0] = PsinkBody(k)->Acc[0];
            SinkForHillExportSend[i][k].Acc[1] = PsinkBody(k)->Acc[1];
            SinkForHillExportSend[i][k].Acc[2] = PsinkBody(k)->Acc[2];
        }
    }

    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

    int NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(SinkForHillExportSend[i],
                NExportThisTime[i]*sizeof(struct StructSinkForHill),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(SinkForHillExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructSinkForHill),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }

#define Xhill (4.0)
    const double Factor = 3*Xhill/(4*M_PI*Pall.GravConst);

    for(int i=0;i<Pall.Nhydro;i++){
        if(AllowSinkConversion[i]){
            // check hill consition
            double Pos[3] = {PhydroBody(i)->PosP[0],PhydroBody(i)->PosP[1],PhydroBody(i)->PosP[2]};
#if 1
            double Acc[3] = {PhydroBody(i)->Acc[0]+Phydro[i]->HydroAcc[0],
                             PhydroBody(i)->Acc[1]+Phydro[i]->HydroAcc[1],
                             PhydroBody(i)->Acc[2]+Phydro[i]->HydroAcc[2]};
#else
            double Acc[3] = {PhydroBody(i)->Acc[0],
                             PhydroBody(i)->Acc[1],
                             PhydroBody(i)->Acc[2]};
#endif
            double rho = Phydro[i]->RhoPred;
            int counter = 0;
            for(int k=0;k<NImport;k++){
                double dr[] = {Pos[0]-SinkForHillExportRecv[k].Pos[0],
                               Pos[1]-SinkForHillExportRecv[k].Pos[1],
                               Pos[2]-SinkForHillExportRecv[k].Pos[2]};
                double da[] = {Acc[0]-SinkForHillExportRecv[k].Acc[0],
                               Acc[1]-SinkForHillExportRecv[k].Acc[1],
                               Acc[2]-SinkForHillExportRecv[k].Acc[2]};

                double rho_hill = -Factor*DOT_PRODUCT(dr,da)/NORM2(dr);
                if(rho > rho_hill){
                    counter ++;
                } else { 
                }
            }
            if(counter != NImport){
                //fprintf(stderr,"Hill condition does not satisfy for %d %d, %g < %g\n",i,k,rho,rho_hill);
                fprintf(stderr,"Hill condition does not satisfy\n");
                AllowSinkConversion[i] = false;
            } else {
                //fprintf(stderr,"Hill condition satisfy for %d %d, %g > %g\n",i,k,rho,rho_hill);
                fprintf(stderr,"Hill condition satisfies\n");
            }
        }
    }

    return ;
}

/*
 * This function checks the sink formation conditions proposed by Federrath et
 * al. ApJ 713 269-290 (2010). In particular, this function checks their condisions a,b,c,d
 * Condition a: the converging flow.
 * Condition b: the local gravitational potential minimum.
 * Condition c: the Jeans-unstable.
 * Condition d: the negative total energy.
 * Extra condision is that the density is higher than a threshold density.
 */
static void CheckSinkFormationCondition(void){

    int counter = 0;
    int Target = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        //AllowSinkConversion[i] = true;

        if(!Phydro[i]->Active){ // Check active
            AllowSinkConversion[i] = false;
            continue;
        }

        counter ++; // Note that the counter is incremented 
                    //  before it used in CalcSizeGetHydroInfo_i.

        if(!PhydroBody(i)->Active){ // Check active
            AllowSinkConversion[i] = false;
            continue;
        }

        if(Phydro[i]->Rho <= Pall.SinkThresholdDensity){ // !condition 0
            AllowSinkConversion[i] = false;
            continue;
        }

#if 0
        if(Phydro[i]->Div > 0.e0){ // !condition a
            AllowSinkConversion[i] = false;
            continue;
        }
#endif 

        
        double PotentialMin;
        double DensityMax;
        double VCOM[3];
        bool NoLocalSink;
        //CalcSizeGetHydroInfo_i(i,&PotentialMin,VCOM);
        //CalcSizeGetHydroInfo_i(counter-1,&PotentialMin,VCOM,&NoLocalSink);
        CalcSizeGetHydroInfo_i(counter-1,&PotentialMin,&DensityMax,VCOM,&NoLocalSink);
#if 1
#if 0
#ifdef USE_SMOOTHED_POTENTIAL //{
        double Ep = Phydro[i]->Pot;
#else // USE_SMOOTHED_POTENTIAL //}//{
        double Ep = PhydroBody(i)->Pot/PhydroBody(i)->Mass;
#endif // USE_SMOOTHED_POTENTIAL //}

        if(Ep != PotentialMin){ // !condition b
            AllowSinkConversion[i] = false;
            // fprintf(stderr,"Potential minimum condition\n");
            continue;
        }
#else
        if(PotentialMin == 0){ // !condition b
            AllowSinkConversion[i] = false;
            // fprintf(stderr,"Potential minimum condition\n");
            continue;
        }
#endif
#else
        if(DensityMax == 0){ // !condition b
            AllowSinkConversion[i] = false;
            // fprintf(stderr,"Potential minimum condition\n");
            continue;
        }

#endif
        if(DensityMax == 0){ // !condition b
            AllowSinkConversion[i] = false;
            // fprintf(stderr,"Potential minimum condition\n");
            continue;
        }


        //////////
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
        if((Pall.Nsink_t>0)&&(!NoLocalSink)){
            AllowSinkConversion[i] = false;
            // fprintf(stderr,"Nosink condition\n");
            // dprintlmpi(300);
            continue;
        }
#endif // USE_SINK_NOLOCALSINK_CONDITION //}

#if 0

        double tsound2 = SQ(2.0*Phydro[i]->Kernel)/(Pall.GGm1*Phydro[i]->U);
        double tdyn2 = 1.0/(4.0*Pall.GravConst*Phydro[i]->Rho);
        if(tsound2 < tdyn2 ){ // !condition c
            AllowSinkConversion[i] = false;
            continue;
        }
#endif

#if 0
        double Eth = PhydroBody(i)->Mass*Phydro[i]->U;
        double Ek = 0.5*PhydroBody(i)->Mass*DISTANCE2(PhydroBody(i)->Vel,VCOM);

        if(Eth + Ek + Ep >= 0.e0){ // !condition d
            AllowSinkConversion[i] = false;
            fprintf(stderr,"Negative energy condition\n");
            continue;
        }
#endif 
        Target ++;
    }

    MPI_Allreduce(MPI_IN_PLACE,&Target,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    
    if(Target > 0){
        //Check Hill condition
        CheckSinkFormationHillCondition();
    }

    return ;
}


#if 1
/*
 * In this function, all gas particles which are active & denser than the
 * threshold density are checked whether they are out of sink radii of sink
 * particles.
 */
static void CheckPermissionForMakeSinkParticle(void){

    // get boundary sink particles.
    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;
    
    // Initialize data arrays.
    static int SinkCheckExportFlagsMaxAllocated = 0;
    static bool (*SinkCheckExportFlags)[NProcs-1];
    if(SinkCheckExportFlagsMaxAllocated < Pall.Nsink){
        if(SinkCheckExportFlagsMaxAllocated > 0){
            free(SinkCheckExportFlags);
        }
        SinkCheckExportFlagsMaxAllocated = (int)(ForAngelsShare*MAX(Pall.Nsink,NAdditionUnit));
        SinkCheckExportFlags = malloc(sizeof(bool)*SinkCheckExportFlagsMaxAllocated*(NProcs-1)+1);
    }
    // Counters for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    // Export structure.
    struct StructSinkCheckExport *SinkCheckExportSend[NProcs];

    // Pick up sink particles who need to be exported.
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckSinkCheckExportFlags(i,NProcs,SinkCheckExportFlags);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructSinkCheckExport),i);
        SinkCheckExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        for(int k=0;k<Pall.Nsink;k++){
            if(SinkCheckExportFlags[k][i]){
                SinkCheckExportSend[i][NExport].AccretionRadius = Psink[k]->AccretionRadius;
                SinkCheckExportSend[i][NExport].Pos[0] = PsinkPosP(k)[0];
                SinkCheckExportSend[i][NExport].Pos[1] = PsinkPosP(k)[1];
                SinkCheckExportSend[i][NExport].Pos[2] = PsinkPosP(k)[2];
                NExport ++;
            }
        }
        //fprintf(stderr,"[%02d] %d %d\n",MPIGetMyID(),NExport,NExportThisTime[i]);
        //fflush(NULL);
        assert(NExport == NExportThisTime[i]);
    }

    // Exchange the number of sink particles who will be exported.
#if 0
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }
#else
    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MPIGetMyID()] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImport += NImportThisTime[i];
    }
    int NImportAll = NImport;
#endif

    // Prepare buffers.
    struct StructSinkCheckExport *SinkCheckExportRecv;
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructSinkCheckExport));
    SinkCheckExportRecv = BufferExportRecv;

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    // Start asynchronous data communication.
#if 0
    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkCheckExportSend[i],
            NExportThisTime[i]*sizeof(struct StructSinkCheckExport),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(SinkCheckExportRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructSinkCheckExport),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

#else

    int counter_send = 0;
    int counter_recv = 0;

    int SendFlag,RecvFlag;
    NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(SinkCheckExportSend[i],
                NExportThisTime[i]*sizeof(struct StructSinkCheckExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(SinkCheckExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructSinkCheckExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
            NImport += NImportThisTime[i];
        }
    }
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
#endif

    // Check accretion condition for hydros.
    // if the gas density is sufficiently high, the check will be done.
    for(int i=0;i<Pall.Nhydro;i++){
        AllowSinkConversion[i] = true;
        if((Phydro[i]->Active) && (Phydro[i]->Rho > Pall.SinkThresholdDensity)){
            double Pos[] = {Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2]};
            for(int k=0;k<Pall.Nsink;k++){ // for local data.
                double Distance2 = DISTANCE2(Pos,PsinkBody(k)->PosP);
                if(SQ(Psink[k]->AccretionRadius) > Distance2)
                    AllowSinkConversion[k] = false;
            }

            for(int k=0;k<NImportAll;k++){ // for imported data.
                double Distance2 = DISTANCE2(Pos,SinkCheckExportRecv[k].Pos);
                if(SQ(SinkCheckExportRecv[k].AccretionRadius) > Distance2)
                    AllowSinkConversion[k] = false;
            }
        }
    }

    return ;
}
#else
static void CheckPermissionForMakeSinkParticle(void){

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;
    
    static int SinkCheckExportFlagsMaxAllocated = 0;
    static bool (*SinkCheckExportFlags)[NProcs-1];
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    struct StructSinkCheckExport *SinkCheckExportSend[NProcs];
    int NImportAll = 0;
    struct StructSinkCheckExport *SinkCheckExportRecv;
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    int NImport = 0;

    if(SinkCheckExportFlagsMaxAllocated < Pall.Nsink){
        if(SinkCheckExportFlagsMaxAllocated > 0){
            free(SinkCheckExportFlags);
        }
        SinkCheckExportFlagsMaxAllocated = (int)(ForAngelsShare*MAX(Pall.Nsink,NAdditionUnit));
        SinkCheckExportFlags = malloc(sizeof(bool)*SinkCheckExportFlagsMaxAllocated*(NProcs-1)+1);
    }


    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckSinkCheckExportFlags(i,NProcs,SinkCheckExportFlags);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructSinkCheckExport),i);
        SinkCheckExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        for(int k=0;k<Pall.Nsink;k++){
            if(SinkCheckExportFlags[k][i]){
                SinkCheckExportSend[i][NExport].AccretionRadius = Psink[k]->AccretionRadius;
                SinkCheckExportSend[i][NExport].Pos[0] = PsinkPosP(k)[0];
                SinkCheckExportSend[i][NExport].Pos[1] = PsinkPosP(k)[1];
                SinkCheckExportSend[i][NExport].Pos[2] = PsinkPosP(k)[2];
                NExport ++;
            }
        }
    }

#if 0
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }
#else
    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MPIGetMyID()] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImport += NImportThisTime[i];
    }
    NImportAll = NImport;
#endif

    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructSinkCheckExport));
    SinkCheckExportRecv = BufferExportRecv;
mtrace();
for(int i=0;i<1000000;i++){
    dprintlmpi(1);

#if 0
    NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkCheckExportSend[i],
            NExportThisTime[i]*sizeof(struct StructSinkCheckExport),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(SinkCheckExportRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructSinkCheckExport),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
#else

    int counter_send = 0;
    int counter_recv = 0;

    int SendFlag,RecvFlag;
    NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(SinkCheckExportSend[i],
                NExportThisTime[i]*sizeof(struct StructSinkCheckExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(SinkCheckExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructSinkCheckExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
            NImport += NImportThisTime[i];
        }
    }
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
#endif
}
muntrace();
MPI_Finalize();
exit(1);

    // Check accretion condition for hydros.
    // if the gas density is sufficiently high, the check will be done.
    for(int i=0;i<Pall.Nhydro;i++){
        AllowSinkConversion[i] = true;
        if((Phydro[i]->Active) && (Phydro[i]->Rho > Pall.SinkThresholdDensity)){
            double Pos[] = {Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2]};
            for(int k=0;k<Pall.Nsink;k++){ // for local data.
                double Distance2 = DISTANCE2(Pos,PsinkBody(k)->PosP);
                if(SQ(Psink[k]->AccretionRadius) > Distance2)
                    AllowSinkConversion[k] = false;
            }

            for(int k=0;k<NImportAll;k++){ // for imported data.
                double Distance2 = DISTANCE2(Pos,SinkCheckExportRecv[k].Pos);
                if(SQ(SinkCheckExportRecv[k].AccretionRadius) > Distance2)
                    AllowSinkConversion[k] = false;
            }
        }
    }

    return ;
}
#endif

void SinkHydroAccretionExported(void){
    SinkHydroAccretion();
    return;
}

static void SinkHydroMergingCheckForLocalData(const int NActives, int ActiveIndexList[restrict],
        struct StructHydroSinkFlags HydroSinkFlags[restrict]){

    // Check Sink conditions for neighbors. This loop is only for the local domain.

    int ThisNode = NONE;
    double TimeStepLimiterFactor = (double)(1<<SINK_TIMESTEP_MAX_K_LOCAL);

    int Neighbors[MaxNeighborSize];
    for(int i=0;i<NActives;i++){  // For local sinks.
        int leaf = ActiveIndexList[i];

        int RootNodeID = 0;
        int CurrentNodeID = HydroNode[RootNodeID].Children;
        do {
            int nlist = 0;
            CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,PsinkPos(leaf),
                    Psink[leaf]->AccretionRadius,&nlist,Neighbors);

            for(int k=0;k<nlist;k++){ // Check Bound conditions.
                int leaf_j = Neighbors[k];
                AllowSinkConversion[leaf_j] = false;
                if((PhydroBody(leaf_j)->Active)&&(Phydro[leaf_j]->Active)){
#if 1
                    if((Phydro[leaf_j]->Rho > SinkGasAbsorptionThresholdDensity)&&
                            SinkHydroBoundCheck(PsinkMass(leaf),PsinkPos(leaf),PsinkVel(leaf),
                        PhydroMass(leaf_j),PhydroPos(leaf_j),PhydroVel(leaf_j),Phydro[leaf_j]->U)&&
                       SinkHydroAngularMomentumCheck(Psink[leaf]->AccretionRadius,
                                PsinkMass(leaf),PsinkPos(leaf),PsinkVel(leaf),PsinkBody(leaf)->Eps,
                                PhydroMass(leaf_j),PhydroPos(leaf_j),PhydroVel(leaf_j),PhydroBody(leaf_j)->Eps)) { // E and L check.
#else
                    if((Phydro[leaf_j]->Rho > SinkGasAbsorptionThresholdDensity)&&
                            SinkHydroBoundCheck(PsinkMass(leaf),PsinkPos(leaf),PsinkVel(leaf),
                        PhydroMass(leaf_j),PhydroPos(leaf_j),PhydroVel(leaf_j),Phydro[leaf_j]->U)){
#endif
                        if(HydroSinkFlags[leaf_j].TargetSinkID == NONE){ // First!
                            HydroSinkFlags[leaf_j].TargetSinkID = leaf;
                            HydroSinkFlags[leaf_j].Distance = DISTANCE(PsinkPos(leaf),PhydroPos(leaf_j));
                        } else { // Else!
                            double distance = DISTANCE(PsinkPos(leaf),PhydroPos(leaf_j));
                            if(distance <= HydroSinkFlags[leaf_j].Distance){
                                HydroSinkFlags[leaf_j].TargetSinkID = leaf;
                                HydroSinkFlags[leaf_j].Distance = distance;
                            }
                        }
                        HydroSinkFlags[leaf_j].NodeID = ThisNode;
                    }
                    Psink[leaf]->dt_localmin = fmin(Psink[leaf]->dt_localmin,
                            TimeStepLimiterFactor*Phydro[leaf_j]->dt_hydro);
                }
            }
        }while(CurrentNodeID != RootNodeID);
    }
    return ;
}

static void SinkStarMergingCheckForLocalData(const int NActives, int ActiveIndexList[restrict],
        struct StructHydroSinkFlags StarSinkFlags[restrict]){

    // Check Sink conditions for neighbors. This loop is only for the local domain.

    int ThisNode = NONE;
    double TimeStepLimiterFactor = (double)(1<<SINK_TIMESTEP_MAX_K_LOCAL);

    int Neighbors[MaxNeighborSize];
    for(int i=0;i<NActives;i++){  // For local sinks.
        int leaf = ActiveIndexList[i];

        int RootNodeID = 0;
        int CurrentNodeID = StellarNode[RootNodeID].Children;
        do {
            int nlist = 0;
            CurrentNodeID = GetNeighborsIterativeApproachFromStellarTree(CurrentNodeID,PsinkPos(leaf),
                    Psink[leaf]->AccretionRadius,&nlist,Neighbors);

            for(int k=0;k<nlist;k++){ // Check Bound conditions.
                int leaf_j = Neighbors[k];
                if(PstarBody(leaf_j)->Active){
                    if(SinkStarBoundCheck(PsinkMass(leaf),PsinkPos(leaf),PsinkVel(leaf),
                                PstarMass(leaf_j),PstarPos(leaf_j),PstarVel(leaf_j))&&
                       SinkStarAngularMomentumCheck(Psink[leaf]->AccretionRadius,
                                PsinkMass(leaf),PsinkPos(leaf),PsinkVel(leaf),PsinkBody(leaf)->Eps,
                                PstarMass(leaf_j),PstarPos(leaf_j),PstarVel(leaf_j),PstarBody(leaf_j)->Eps)) { // E and L check.
                        if(StarSinkFlags[leaf_j].TargetSinkID == NONE){ // First!
                            StarSinkFlags[leaf_j].TargetSinkID = leaf;
                            StarSinkFlags[leaf_j].Distance = DISTANCE(PsinkPos(leaf),PstarPos(leaf_j));
                        } else { // Else!
                            double distance = DISTANCE(PsinkPos(leaf),PstarPos(leaf_j));
                            if(distance <= StarSinkFlags[leaf_j].Distance){
                                StarSinkFlags[leaf_j].TargetSinkID = leaf;
                                StarSinkFlags[leaf_j].Distance = distance;
                            }
                        }
                        StarSinkFlags[leaf_j].NodeID = ThisNode;
                    }
                    Psink[leaf]->dt_localmin = fmin(Psink[leaf]->dt_localmin,
                            TimeStepLimiterFactor*PstarBody(leaf_j)->dt);
                }
            }
        }while(CurrentNodeID != RootNodeID);
    }
    return ;
}


static void SinkHydroMergingCheckForImportedData(const int NProcs, int NImportThisTime[restrict],
        struct StructSinkExport SinkExportRecv[restrict], struct StructHydroSinkFlags HydroSinkFlags[restrict]){

    int NImport = 0;
    int Neighbors[MaxNeighborSize];
    double TimeStepLimiterFactor = (double)(1<<SINK_TIMESTEP_MAX_K_LOCAL);
    for(int i=0;i<NProcs-1;i++){ // For imported sinks.
        for(int k=0;k<NImportThisTime[i];k++){

            double Pos[3] = {SinkExportRecv[NImport].Pos[0],
                             SinkExportRecv[NImport].Pos[1],
                             SinkExportRecv[NImport].Pos[2]};
            double Vel[3] = {SinkExportRecv[NImport].Vel[0],
                             SinkExportRecv[NImport].Vel[1],
                             SinkExportRecv[NImport].Vel[2]};
            double AccretionRadius = SinkExportRecv[NImport].AccretionRadius;
            double Mass = SinkExportRecv[NImport].Mass;
            double Eps = SinkExportRecv[NImport].Eps;
                
            int Nlist = GetNeighborsLimited(Pos,AccretionRadius,Neighbors);
            for(int l=0;l<Nlist;l++){ // Check Bound conditions.
                int leaf_j = Neighbors[l];
                AllowSinkConversion[leaf_j] = false;
                if((PhydroBody(leaf_j)->Active)&&(Phydro[leaf_j]->Active)){
#if 1
                    if((Phydro[leaf_j]->Rho > SinkGasAbsorptionThresholdDensity)&&
                        SinkHydroBoundCheck(Mass,Pos,Vel,
                            PhydroMass(leaf_j),PhydroPos(leaf_j),PhydroVel(leaf_j),Phydro[leaf_j]->U) &&
                        SinkHydroAngularMomentumCheck(AccretionRadius,Mass,Pos,Vel,Eps,
                                PhydroMass(leaf_j),PhydroPos(leaf_j),PhydroVel(leaf_j),PhydroBody(leaf_j)->Eps)) { // E and L check.
#else
                    if((Phydro[leaf_j]->Rho > SinkGasAbsorptionThresholdDensity)&&
                        SinkHydroBoundCheck(Mass,Pos,Vel,
                            PhydroMass(leaf_j),PhydroPos(leaf_j),PhydroVel(leaf_j),Phydro[leaf_j]->U)) { // E and L check.

#endif
                        if(HydroSinkFlags[leaf_j].TargetSinkID == NONE){ // First!
                            HydroSinkFlags[leaf_j].TargetSinkID = SinkExportRecv[NImport].Index;
                            HydroSinkFlags[leaf_j].Distance = DISTANCE(Pos,PhydroPos(leaf_j));
                            HydroSinkFlags[leaf_j].NodeID = i;
                        } else { // Else!
                            double distance = DISTANCE(Pos,PhydroPosP(leaf_j));
                            if(distance < HydroSinkFlags[leaf_j].Distance){
                                HydroSinkFlags[leaf_j].TargetSinkID = SinkExportRecv[NImport].Index;
                                HydroSinkFlags[leaf_j].Distance = DISTANCE(Pos,PhydroPos(leaf_j));
                                HydroSinkFlags[leaf_j].NodeID = i;
                            }
                        }
                    }
                    SinkExportRecv[NImport].dt_localmin = fmin(SinkExportRecv[NImport].dt_localmin,
                            TimeStepLimiterFactor*Phydro[leaf_j]->dt_hydro);
                }
            }
            NImport ++;
        }
    }

    return ;
}

static void SinkStarMergingCheckForImportedData(const int NProcs, int NImportThisTime[restrict],
        struct StructSinkExport SinkExportRecv[restrict], struct StructHydroSinkFlags StarSinkFlags[restrict]){

    int NImport = 0;
    int Neighbors[MaxNeighborSize];
    double TimeStepLimiterFactor = (double)(1<<SINK_TIMESTEP_MAX_K_LOCAL);
    for(int i=0;i<NProcs-1;i++){ // For imported sinks.
        for(int k=0;k<NImportThisTime[i];k++){

            double Pos[3] = {SinkExportRecv[NImport].Pos[0],
                             SinkExportRecv[NImport].Pos[1],
                             SinkExportRecv[NImport].Pos[2]};
            double Vel[3] = {SinkExportRecv[NImport].Vel[0],
                             SinkExportRecv[NImport].Vel[1],
                             SinkExportRecv[NImport].Vel[2]};
            double AccretionRadius = SinkExportRecv[NImport].AccretionRadius;
            double Mass = SinkExportRecv[NImport].Mass;
            double Eps = SinkExportRecv[NImport].Eps;
                
            int Nlist = GetNeighborsLimitedFromStellarTree(Pos,AccretionRadius,Neighbors);

            for(int l=0;l<Nlist;l++){ // Check Bound conditions.
                int leaf_j = Neighbors[l];
                if(PstarBody(leaf_j)->Active){
                    if(SinkStarBoundCheck(Mass,Pos,Vel,
                        PstarMass(leaf_j),PstarPos(leaf_j),PstarVel(leaf_j)) &&
                       SinkStarAngularMomentumCheck(AccretionRadius,Mass,Pos,Vel,Eps,
                                PstarMass(leaf_j),PstarPos(leaf_j),PstarVel(leaf_j),PstarBody(leaf_j)->Eps)) { // E and L check.
                        if(StarSinkFlags[leaf_j].TargetSinkID == NONE){ // First!
                            StarSinkFlags[leaf_j].TargetSinkID = SinkExportRecv[NImport].Index;
                            StarSinkFlags[leaf_j].Distance = DISTANCE(Pos,PstarPos(leaf_j));
                            StarSinkFlags[leaf_j].NodeID = i;
                        } else { // Else!
                            double distance = DISTANCE(Pos,PhydroPosP(leaf_j));
                            if(distance < StarSinkFlags[leaf_j].Distance){
                                StarSinkFlags[leaf_j].TargetSinkID = SinkExportRecv[NImport].Index;
                                StarSinkFlags[leaf_j].Distance = DISTANCE(Pos,PstarPos(leaf_j));
                                StarSinkFlags[leaf_j].NodeID = i;
                            }
                        }
                    }
                    SinkExportRecv[NImport].dt_localmin = fmin(SinkExportRecv[NImport].dt_localmin,
                            TimeStepLimiterFactor*PstarBody(leaf_j)->dt);
                }
            }
            NImport ++;
        }
    }

    return ;
}

static void UpdateSinkLocalTimeStep(const int NExportThisTime[restrict], const int NImportThisTime[restrict], 
        const int NProcs, struct StructSinkExport *SinkExportSend[NProcs],
            struct StructSinkExport SinkExportRecv[restrict]){

#ifdef USE_SINK_TIMESTEP_LIMITER

    double *SinkLocalTimeStepSend;
    double *SinkLocalTimeStepRecv[NProcs];
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++)
        NImportAll += NImportThisTime[i];
    SinkLocalTimeStepSend = malloc(sizeof(double)*NImportAll+1);
    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        SinkLocalTimeStepRecv[i] = malloc(sizeof(double)*NExportThisTime[i]+1);
        for(int k=0;k<NImportThisTime[i];k++){
            SinkLocalTimeStepSend[NImportAll] = SinkExportRecv[NImportAll].dt_localmin;
            NImportAll ++;
        }
    }

    MPI_Status mpi_status_Export_Send_LocalTimeStep[NProcs-1];
    MPI_Request mpi_request_Export_Send_LocalTimeStep[NProcs-1];
    MPI_Status mpi_status_Export_Recv_LocalTimeStep[NProcs-1];
    MPI_Request mpi_request_Export_Recv_LocalTimeStep[NProcs-1];

    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkLocalTimeStepSend+NImport,NImportThisTime[i],
                MPI_DOUBLE,CommunicationTable[i].RecvRank,TAG_SINK_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send_LocalTimeStep+i);
        MPI_Irecv(SinkLocalTimeStepRecv[i],NExportThisTime[i],
                MPI_DOUBLE,CommunicationTable[i].SendRank,TAG_SINK_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv_LocalTimeStep+i);
        NImport += NImportThisTime[i];
    }
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv_LocalTimeStep,mpi_status_Export_Recv_LocalTimeStep);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send_LocalTimeStep,mpi_status_Export_Send_LocalTimeStep);
    assert(NImport == NImportAll);

    for(int i=0;i<NProcs-1;i++){
        for(int k=0;k<NExportThisTime[i];k++){
            int leaf = (*(SinkExportSend+i))[k].Index;
            Psink[leaf]->dt_localmin = fmin(Psink[leaf]->dt_localmin,SinkLocalTimeStepRecv[i][k]);
        }
    }

    free(SinkLocalTimeStepSend);
    for(int i=0;i<NProcs-1;i++)
        free(SinkLocalTimeStepRecv[i]);
#endif

    return ;
}

/*
 * This function handles the gas accretion onto the sink particles.
 */
static void SinkHydroAccretion(void){

    if((Pall.NActivesSink_t == 0)||(Pall.Nhydro_t == 0)) 
        return;

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    int ThisNode = NONE;

    // Initialize data arrays.
    static int SinkExportFlagsMaxAllocated = 0;
    static bool (*SinkExportFlags)[NProcs-1];
    static int *ActiveIndexList;
    if(SinkExportFlagsMaxAllocated < Pall.Nsink){
        if(SinkExportFlagsMaxAllocated > 0){
            free(ActiveIndexList);
            free(SinkExportFlags);
        }
        SinkExportFlagsMaxAllocated = (int)(ForAngelsShare*Pall.Nsink);
        SinkExportFlags = malloc(sizeof(bool)*SinkExportFlagsMaxAllocated*(NProcs-1)+1);
        ActiveIndexList = malloc(sizeof(int)*SinkExportFlagsMaxAllocated);
    }

    // Find active sink particles and clear corresponding export_flags.
    int NActives = 0;
    for(int i=0;i<Pall.Nsink;i++){
        for(int k=0;k<NProcs-1;k++)
            SinkExportFlags[i][k] = false;

        if(PsinkActive(i)){
            ActiveIndexList[NActives] = i;
            NActives ++;
        }
    }

    // Counters for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    // Export structure.
    struct StructSinkExport *SinkExportSend[NProcs];



#if 0
    // Pick up sink particles who need to be exported.
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckSinkExportFlags(i,NProcs,SinkExportFlags);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructSinkExport),i);
        SinkExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        for(int k=0;k<NActives;k++){
            int leaf = ActiveIndexList[k];
            if(SinkExportFlags[leaf][i]){
                SinkExportSend[i][NExport].Index  = leaf;
                SinkExportSend[i][NExport].AccretionRadius = Psink[leaf]->AccretionRadius;
                SinkExportSend[i][NExport].Mass   = PsinkMass(leaf);
                SinkExportSend[i][NExport].Pos[0] = PsinkPosP(leaf)[0];
                SinkExportSend[i][NExport].Pos[1] = PsinkPosP(leaf)[1];
                SinkExportSend[i][NExport].Pos[2] = PsinkPosP(leaf)[2];

                SinkExportSend[i][NExport].Vel[0] = PsinkBody(leaf)->Vel[0];
                SinkExportSend[i][NExport].Vel[1] = PsinkBody(leaf)->Vel[1];
                SinkExportSend[i][NExport].Vel[2] = PsinkBody(leaf)->Vel[2];
                SinkExportSend[i][NExport].dt_localmin = Psink[leaf]->dt_localmin;
                SinkExportSend[i][NExport].Eps = PsinkBody(leaf)->Eps;
                NExport ++;
            }
        }
        assert(NExport == NExportThisTime[i]);
    }

    // Exchange the number of sink particles who will be exported.
#if 0
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }
#else
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
#endif


#else 

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    int SendFlag,RecvFlag;

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(NImportThisTime+i,1,MPI_INT,
            CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM+i,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
    }

    // Pick up sink particles who need to be exported.
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckSinkExportFlags(i,NProcs,SinkExportFlags);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructSinkExport),i);
        SinkExportSend[i] = BufferExportSend[i];

        MPI_Isend(NExportThisTime+i,1,MPI_INT,
                CommunicationTable[i].SendRank,TAG_SINK_PRECOMM+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);

        int NExport = 0;
        for(int k=0;k<NActives;k++){
            int leaf = ActiveIndexList[k];
            if(SinkExportFlags[leaf][i]){
                SinkExportSend[i][NExport].Index  = leaf;
                SinkExportSend[i][NExport].AccretionRadius = Psink[leaf]->AccretionRadius;
                SinkExportSend[i][NExport].Mass   = PsinkMass(leaf);
                SinkExportSend[i][NExport].Pos[0] = PsinkPosP(leaf)[0];
                SinkExportSend[i][NExport].Pos[1] = PsinkPosP(leaf)[1];
                SinkExportSend[i][NExport].Pos[2] = PsinkPosP(leaf)[2];

                SinkExportSend[i][NExport].Vel[0] = PsinkBody(leaf)->Vel[0];
                SinkExportSend[i][NExport].Vel[1] = PsinkBody(leaf)->Vel[1];
                SinkExportSend[i][NExport].Vel[2] = PsinkBody(leaf)->Vel[2];
                SinkExportSend[i][NExport].dt_localmin = Psink[leaf]->dt_localmin;
                SinkExportSend[i][NExport].Eps = PsinkBody(leaf)->Eps;
                NExport ++;
            }
        }
        assert(NExport == NExportThisTime[i]);
    }

    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportAll += NImportThisTime[i];
    }

#endif


    // Prepare buffers.
    struct StructSinkExport *SinkExportRecv;
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructSinkExport));
    SinkExportRecv = BufferExportRecv;

    // MPI_Status mpi_status_Export_Send[NProcs-1];
    // MPI_Request mpi_request_Export_Send[NProcs-1];
    // MPI_Status mpi_status_Export_Recv[NProcs-1];
    // MPI_Request mpi_request_Export_Recv[NProcs-1];


    // Start asynchronous data communication.
    int NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(SinkExportSend[i],
                NExportThisTime[i]*sizeof(struct StructSinkExport),MPI_BYTE,
                    CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(SinkExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructSinkExport),MPI_BYTE,
                    CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);

    
#ifdef SINK_STAR_PARTICLES
    // Make stellar tree!
    PlantStellarTree();
#endif // SINK_STAR_PARTICLES

    // Prepare sink flag arrays for hydro particles.
    struct StructHydroSinkFlags *HydroSinkFlags;
    HydroSinkFlags = malloc(sizeof(struct StructHydroSinkFlags)*Pall.Nhydro);
    for(int i=0;i<Pall.Nhydro;i++)
        HydroSinkFlags[i].TargetSinkID = ThisNode;

#ifdef SINK_STAR_PARTICLES
    // Prepare sink flag arrays for star particles.
    struct StructHydroSinkFlags *StarSinkFlags;
    StarSinkFlags = malloc(sizeof(struct StructHydroSinkFlags)*Pall.Nstars);
    for(int i=0;i<Pall.Nstars;i++)
        StarSinkFlags[i].TargetSinkID = ThisNode;
#endif // SINK_STAR_PARTICLES


    SinkHydroMergingCheckForLocalData(NActives,ActiveIndexList,HydroSinkFlags); 
#ifdef SINK_STAR_PARTICLES
    SinkStarMergingCheckForLocalData(NActives,ActiveIndexList,HydroSinkFlags); 
#endif // SINK_STAR_PARTICLES

    // MPI wait
    double TimeComm = GetElapsedTime();
    // MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    // MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.SinkCommThisStep += GetElapsedTime()-TimeComm;

    SinkHydroMergingCheckForImportedData(NProcs,NImportThisTime,SinkExportRecv,HydroSinkFlags);
#ifdef SINK_STAR_PARTICLES
    SinkStarMergingCheckForImportedData(NProcs,NImportThisTime,SinkExportRecv,StarSinkFlags);
#endif // SINK_STAR_PARTICLES

    UpdateSinkLocalTimeStep(NExportThisTime,NImportThisTime,NProcs,SinkExportSend,SinkExportRecv);


    // Counters for export and import. 
    int NExportSinkThisTime[NProcs];
    int NImportSinkThisTime[NProcs];
    memset(NExportSinkThisTime,0,sizeof(int)*NProcs);
    memset(NImportSinkThisTime,0,sizeof(int)*NProcs);

    int NLocalVictims = 0;
    int NExportSinkThisTimeAll = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroSinkFlags[i].TargetSinkID != NONE){
            if(HydroSinkFlags[i].NodeID == ThisNode){
                NLocalVictims ++;
            } else {
                NExportSinkThisTime[HydroSinkFlags[i].NodeID] ++;
                NExportSinkThisTimeAll ++;
            }
        }
    }
#ifdef SINK_STAR_PARTICLES
    int NLocalVictimStars = 0;
    int NExportSinkThisTimeAllStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(StarSinkFlags[i].TargetSinkID != NONE){
            if(StarSinkFlags[i].NodeID == ThisNode){
                NLocalVictims ++;
                NLocalVictimStars ++;
            } else {
                NExportSinkThisTime[StarSinkFlags[i].NodeID] ++;
                NExportSinkThisTimeAll ++;
                NExportSinkThisTimeAllStars ++;
            }
        }
    }
#endif //SINK_STAR_PARTICLES


    struct StructSinkImport *SinkImportSend;
    CheckSizeofBufferImportSend(NExportSinkThisTimeAll,sizeof(struct StructSinkImport));
    SinkImportSend = BufferImportSend;

    // prepare headers and counters.
    int Headers[MPIGetNumProcs()];
    int Counters[MPIGetNumProcs()];
    for(int i=0;i<MPIGetNumProcs();i++) Counters[i] = 0;
    Headers[0] = 0;
    for(int i=1;i<MPIGetNumProcs();i++)
        Headers[i] = Headers[i-1]+ NExportSinkThisTime[i-1];

    int check_vanished_export = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroSinkFlags[i].TargetSinkID != NONE){
            if(HydroSinkFlags[i].NodeID != ThisNode){
                int Index = Headers[HydroSinkFlags[i].NodeID] + Counters[HydroSinkFlags[i].NodeID];

                SinkImportSend[Index].Pos[0] = PhydroPos(i)[0];
                SinkImportSend[Index].Pos[1] = PhydroPos(i)[1];
                SinkImportSend[Index].Pos[2] = PhydroPos(i)[2];

                SinkImportSend[Index].Vel[0] = PhydroVel(i)[0];
                SinkImportSend[Index].Vel[1] = PhydroVel(i)[1];
                SinkImportSend[Index].Vel[2] = PhydroVel(i)[2];

#ifdef USE_CELIB //{
                for(int k=0;k<CELibYield_Number;k++){
                    SinkImportSend[Index].Elements[k] = Phydro[i]->Elements[k];
                }
#endif // USE_CELIB //}
                SinkImportSend[Index].Mass = PhydroMass(i);
                SinkImportSend[Index].TargetSinkIndex = HydroSinkFlags[i].TargetSinkID;
                SinkImportSend[Index].Type = PhydroBody(i)->Type;
                Phydro[i]->Use = false;
                PhydroBody(i)->Use = false;
                check_vanished_export ++;

                Counters[HydroSinkFlags[i].NodeID] ++;
            }
        }
    }
#ifdef SINK_STAR_PARTICLES
    for(int i=0;i<Pall.Nstars;i++){
        if(StarSinkFlags[i].TargetSinkID != NONE){
            if(StarSinkFlags[i].NodeID != ThisNode){
                int Index = Headers[StarSinkFlags[i].NodeID] + Counters[StarSinkFlags[i].NodeID];

                SinkImportSend[Index].Pos[0] = PstarPos(i)[0];
                SinkImportSend[Index].Pos[1] = PstarPos(i)[1];
                SinkImportSend[Index].Pos[2] = PstarPos(i)[2];

                SinkImportSend[Index].Vel[0] = PstarVel(i)[0];
                SinkImportSend[Index].Vel[1] = PstarVel(i)[1];
                SinkImportSend[Index].Vel[2] = PstarVel(i)[2];

#ifdef USE_CELIB //{
                for(int k=0;k<CELibYield_Number;k++){
                    SinkImportSend[Index].Elements[k] = Pstar[i]->Elements[k];
                }
#endif // USE_CELIB //}

                SinkImportSend[Index].Mass = PstarMass(i);
                SinkImportSend[Index].TargetSinkIndex = StarSinkFlags[i].TargetSinkID;
                SinkImportSend[Index].Type = PstarBody(i)->Type;
                Pstar[i]->Use = false;
                PstarBody(i)->Use = false;
                check_vanished_export ++;

                Counters[StarSinkFlags[i].NodeID] ++;
            }
        }
    }
#endif // SINK_STAR_PARTICLES

    ////////////////////////////////////////////////////////
    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportSinkThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
            NImportSinkThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportSinkThisTime[i];
    }

    struct StructSinkImport *SinkImportRecv[NProcs];
    CheckSizeofBufferImportRecv(NImportSinkThisTime,sizeof(struct StructSinkImport));
    for(int i=0;i<NProcs-1;i++)
        SinkImportRecv[i] = BufferImportRecv[i];

    NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkImportSend+NImport,
            NExportSinkThisTime[i]*sizeof(struct StructSinkImport),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SINK_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(SinkImportRecv[i],
            NImportSinkThisTime[i]*sizeof(struct StructSinkImport),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_SINK_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NExportSinkThisTime[i];
    }


    ///////////////////////////////////////////
    // Sink operations for the local domain. //
    ///////////////////////////////////////////

#ifdef USE_CELIB //{  
#define InputLocalHydroElements(_x) (Phydro[_x]->Elements)
#else // USE_CELIB //}//{ 
#define InputLocalHydroElements(_x) (NULL)
#endif // USE_CELIB //}

    int check_vanished_local = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroSinkFlags[i].TargetSinkID != NONE){
            if(HydroSinkFlags[i].NodeID == ThisNode){
                ParticleAbsorption(HydroSinkFlags[i].TargetSinkID,
                        PhydroPos(i),PhydroVel(i),PhydroMass(i),InputLocalHydroElements(i),TypeHydro);
                Phydro[i]->Use = false;
                PhydroBody(i)->Use = false;
                check_vanished_local ++;
            }
        }
    }


#ifdef SINK_STAR_PARTICLES
    for(int i=0;i<Pall.Nstars;i++){
        if(StarSinkFlags[i].TargetSinkID != NONE){
            if(StarSinkFlags[i].NodeID == ThisNode){
                ParticleAbsorption(StarSinkFlags[i].TargetSinkID,
                        PstarPos(i),PstarVel(i),PstarMass(i),TypeStar);
                Pstar[i]->Use = false;
                PstarBody(i)->Use = false;
                check_vanished_local ++;
            }
        }
    }
#endif //SINK_STAR_PARTICLES


    TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.SinkCommThisStep += GetElapsedTime()-TimeComm;

#ifdef USE_CELIB //{  
#define InputImportHydroElements(_x,_y) (SinkImportRecv[_x][_y].Elements)
#else // USE_CELIB //}//{ 
#define InputImportHydroElements(_x,_y) (NULL)
#endif // USE_CELIB //}


    ////////////////////////////////////////
    // Sink operations for imported data. //
    ////////////////////////////////////////

#if 0

    for(int i=0;i<NProcs-1;i++){
        if(NImportSinkThisTime[i]>0) dprintlmpi(NImportSinkThisTime[i]);
        for(int k=0;k<NImportSinkThisTime[i];k++){
            dprintlmpi(SinkImportRecv[i][k].TargetSinkIndex);
            gprintlmpi(SinkImportRecv[i][k].Pos[0]);
            gprintlmpi(SinkImportRecv[i][k].Mass);
            ParticleAbsorption(SinkImportRecv[i][k].TargetSinkIndex,
                    SinkImportRecv[i][k].Pos,SinkImportRecv[i][k].Vel,
                    SinkImportRecv[i][k].Mass,InputImportHydroElements(i,k),SinkImportRecv[i][k].Type);
        }
    }
#else
    if(Pall.Nsink > 0){
        for(int i=0;i<NProcs-1;i++){
            for(int k=0;k<NImportSinkThisTime[i];k++){
                ParticleAbsorption(SinkImportRecv[i][k].TargetSinkIndex,
                        SinkImportRecv[i][k].Pos,SinkImportRecv[i][k].Vel,
                        SinkImportRecv[i][k].Mass,InputImportHydroElements(i,k),SinkImportRecv[i][k].Type);
            }
        }
    }

#endif

    int TotalVictims = NLocalVictims + NExportSinkThisTimeAll;


    ////////////////////// HERE IS THE FRONTLINE!! /////////////////////// 
    int GlobalTotalVictims;
    MPI_Allreduce(&TotalVictims,&GlobalTotalVictims,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        dprintlmpi(GlobalTotalVictims);
    }
    if(GlobalTotalVictims > 0){
#if ShowLog //{
        if(NLocalVictims>0){
            fprintf(stderr,"[%02d] Local and importerd victims = %d and %d %d:%s:%s\n",
                    MPIGetMyID(),NLocalVictims,NExportSinkThisTimeAll,
                    __LINE__,__FUNCTION__,__FILE__);
            fflush(stderr);
        }
#endif // ShowLog //}

#ifdef SINK_STAR_PARTICLES
        int VictimsStars = NLocalVictimStars+NExportSinkThisTimeAllStars;
        int VictimsHydro = TotalVictims-VictimsStars;
        Pall.Nhydro -= VictimsHydro;
        Pall.Nstars -= VictimsStars;
        Pall.Ntotal -= TotalVictims;
        fprintf(stderr,"Eating [%02d] Vt %d, Vh %d, Vs %d\n",MPIGetMyID(),
                TotalVictims,VictimsHydro,VictimsStars);

        int LogHydro,LogStars;
        LogHydro = VictimsHydro;
        LogStars = VictimsStars;
        int GlobalLogHydro,GlobalLogStars;
        MPI_Allreduce(&LogHydro,&GlobalLogHydro,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&LogStars,&GlobalLogStars,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        //if(MPIGetMyID() == MPI_ROOT_RANK){
        if(Pall.Nsink > 0){
            FILE *fp_;
            FileOpen(fp_,"./sink.dat","a");
            fprintf(fp_,"%g %d %d %g\n",Pall.TCurrent,
                    GlobalLogHydro,GlobalLogStars,PsinkBody(0)->Mass);
            fflush(fp_);
            fclose(fp_);
        }


#else  // SINK_STAR_PARTICLES
        Pall.Nhydro -= TotalVictims;
        Pall.Ntotal -= TotalVictims;
#endif //SINK_STAR_PARTICLES

        ReConnectPointers();
        UpdateTotalNumber();
        UpdateTotalActiveNumber(); // What is this?
    }


    free(HydroSinkFlags);
#ifdef SINK_STAR_PARTICLES
    free(StarSinkFlags);
#endif // SINK_STAR_PARTICLES

    return ;
}


#ifdef USE_SINK_PRESSURE_FORCE
/*
 * This function calculates the pressure forces from sinks to gas particles.
 */

struct StructSinkPosSizeExport{ 
    int Index; // ActiveIndex.
    double  AccretionRadius;  // Sink radius.
    double  Pos[3];  // Positions.
};

struct StructSinkPressureForces{
    double    Mass;    // Mass.
    double    Kernel;  // Kernel size.
    double    Pos[3];  // Position.
    double    Vel[3];  // Volocity.
    double    Rho;     // Density.
    int       Nlist;
};

static void CalcSinkPressureForces(void){

    if((Pall.NActivesSink_t == 0)||(Pall.Nhydro_t == 0)) 
        return;

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    int ThisNode = NONE;

    // Initialize data arrays.
    static int SinkExportFlagsMaxAllocated = 0;
    static bool (*SinkExportFlags)[NProcs-1];
    static int *ActiveIndexList;
    if(SinkExportFlagsMaxAllocated < Pall.Nsink){
        if(SinkExportFlagsMaxAllocated > 0){
            free(ActiveIndexList);
            free(SinkExportFlags);
        }
        SinkExportFlagsMaxAllocated = (int)(ForAngelsShare*Pall.Nsink);
        SinkExportFlags = malloc(sizeof(bool)*SinkExportFlagsMaxAllocated*(NProcs-1)+1);
        ActiveIndexList = malloc(sizeof(int)*SinkExportFlagsMaxAllocated);
    }

    // Find active sink particles and clear corresponding export_flags.
    int NActives = 0;
    for(int i=0;i<Pall.Nsink;i++){
        if(PsinkActive(i)){
            ActiveIndexList[NActives] = i;
            for(int k=0;k<NProcs-1;k++)
                SinkExportFlags[i][k] = 0;
            NActives ++;
        }
    }

    // Counters for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    // Export structure.
    struct StructSinkPressureForces *SinkPressureForcesSend[NProcs];

    // Pick up sink particles who need to be exported.
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckSinkExportFlags(i,NProcs,SinkExportFlags);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructSinkPressureForces),i);
        SinkPressureForcesSend[i] = BufferExportSend[i];

        int NExport = 0;
        for(int k=0;k<NActives;k++){
            int leaf = ActiveIndexList[k];
            if(SinkExportFlags[leaf][i]){
                SinkPressureForcesSend[i][NExport].Mass = PsinkMass(leaf);
                SinkPressureForcesSend[i][NExport].Kernel = SINK_PRESSURE_FORCE_SCALEFACTOR*Psink[leaf]->AccretionRadius;
                SinkPressureForcesSend[i][NExport].Pos[0] = PsinkPosP(leaf)[0];
                SinkPressureForcesSend[i][NExport].Pos[1] = PsinkPosP(leaf)[1];
                SinkPressureForcesSend[i][NExport].Pos[2] = PsinkPosP(leaf)[2];

                SinkPressureForcesSend[i][NExport].Vel[0] = PsinkVel(leaf)[0];
                SinkPressureForcesSend[i][NExport].Vel[1] = PsinkVel(leaf)[1];
                SinkPressureForcesSend[i][NExport].Vel[2] = PsinkVel(leaf)[2];
                NExport ++;
            }
        }
        assert(NExport == NExportThisTime[i]);
    }

    // Exchange the number of sink particles who will be exported.
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }

    // Prepare buffers.
    struct StructSinkPressureForces *SinkPressureForcesRecv;
    CheckSizeofBufferExportRecv(NActives+NImportAll,sizeof(struct StructSinkPressureForces));
    SinkPressureForcesRecv = BufferExportRecv;

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    // Start asynchronous data communication.
    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkPressureForcesSend[i],
            NExportThisTime[i]*sizeof(struct StructSinkPressureForces),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(SinkPressureForcesRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructSinkPressureForces),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);

    // Copy the local data.
    for(int i=0;i<NActives;i++){
        int leaf = ActiveIndexList[i];
        SinkPressureForcesRecv[i+NImportAll].Mass = PsinkMass(leaf);
        SinkPressureForcesRecv[i+NImportAll].Kernel = SINK_PRESSURE_FORCE_SCALEFACTOR*Psink[leaf]->AccretionRadius;

        SinkPressureForcesRecv[i+NImportAll].Pos[0] = PsinkPosP(leaf)[0];
        SinkPressureForcesRecv[i+NImportAll].Pos[1] = PsinkPosP(leaf)[1];
        SinkPressureForcesRecv[i+NImportAll].Pos[2] = PsinkPosP(leaf)[2];

        SinkPressureForcesRecv[i+NImportAll].Vel[0] = PsinkVel(leaf)[0];
        SinkPressureForcesRecv[i+NImportAll].Vel[1] = PsinkVel(leaf)[1];
        SinkPressureForcesRecv[i+NImportAll].Vel[2] = PsinkVel(leaf)[2];
    }

    double *ArrayDensityMax;
    ArrayDensityMax = malloc(sizeof(double)*(NActives+NImportAll));
    
    // Neighbor search and check maximum density.
    int Neighbors[MaxNeighborSize];
    for(int i=0;i<NActives;i++){
        int Nlist = GetNeighborsLimited(SinkPressureForcesRecv[i+NImportAll].Pos,SinkPressureForcesRecv[i+NImportAll].Kernel,Neighbors);
        SinkPressureForcesRecv[i+NImportAll].Nlist = Nlist;
        ArrayDensityMax[i+NImportAll] = 0.e0;
        for(int k=0;k<Nlist;k++){
            if(Phydro[Neighbors[k]]->Active)
                ArrayDensityMax[i+NImportAll] = fmax(ArrayDensityMax[i+NImportAll],Phydro[Neighbors[k]]->RhoPred);
                //ArrayDensityMax[i+NImportAll] = fmax(ArrayDensityMax[i+NImportAll],NORM(Phydro[Neighbors[k]]->HydroAcc));
        }
    }

    // MPI wait
    double TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.SinkCommThisStep += GetElapsedTime()-TimeComm;

    for(int i=0;i<NImportAll;i++){
        int Nlist = GetNeighborsLimited(SinkPressureForcesRecv[i].Pos,SinkPressureForcesRecv[i].Kernel,Neighbors);
        SinkPressureForcesRecv[i].Nlist = Nlist;
        ArrayDensityMax[i] = 0.e0;
        for(int k=0;k<Nlist;k++){
            if(Phydro[Neighbors[k]]->Active)
                ArrayDensityMax[i] = fmax(ArrayDensityMax[i],Phydro[Neighbors[k]]->RhoPred);
                //ArrayDensityMax[i] = fmax(ArrayDensityMax[i],NORM(Phydro[Neighbors[k]]->HydroAcc));
        }
    }

    // All reduce!
    double *GlobalArrayDensityMax;
    GlobalArrayDensityMax = malloc(sizeof(double)*(NActives+NImportAll));
    MPI_Allreduce(ArrayDensityMax,GlobalArrayDensityMax,(NActives+NImportAll),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    for(int i=0;i<(NActives+NImportAll);i++){
        SinkPressureForcesRecv[i].Rho = GlobalArrayDensityMax[i];
        //SinkPressureForcesRecv[i].Rho = GlobalArrayDensityMax[i]/3.0;
        SinkPressureForcesRecv[i].Mass = (4.0*M_PI/3.0)*GlobalArrayDensityMax[i]*CUBE(2*SinkPressureForcesRecv[i].Kernel);
        // CallHydroAcc
        if(SinkPressureForcesRecv[i].Nlist > 0)
            CalcSinkPressureForcesToHydro(
                    SinkPressureForcesRecv[i].Pos,
                    SinkPressureForcesRecv[i].Vel,
                    SinkPressureForcesRecv[i].Mass,
                    SinkPressureForcesRecv[i].Kernel,
                    SinkPressureForcesRecv[i].Rho);
    }

    // Is density max surely required?
    
    // >>>>>>>>>>>>>>>>>>> OK.
    



    free(ArrayDensityMax);
    free(GlobalArrayDensityMax);
    //free(SinkPressureForcesRecv);
    //free(HydroSinkFlags);

    return ;
}
#endif //USE_SINK_PRESSURE_FORCE

#if 0
/*
 * This function should be called when the run without self/external-gravity to
 * determine the sink time-step.
 */
static void SinkTimeStepsForNonGravityRuns(void){

    if(Pall.Nsink_t == 0) return;

    if(InitSinkParticlesFlag == false){
        InitSinkParticles();
        InitSinkParticlesFlag = true;
    }


    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    int ThisNode = NONE;

    // Initialize data arrays.
    static int SinkExportFlagsMaxAllocated = 0;
    static bool (*SinkExportFlags)[NProcs-1];
    static int *ActiveIndexList;
    if(SinkExportFlagsMaxAllocated < Pall.Nsink){
        if(SinkExportFlagsMaxAllocated > 0){
            free(ActiveIndexList);
            free(SinkExportFlags);
        }
        SinkExportFlagsMaxAllocated = (int)(ForAngelsShare*Pall.Nsink);
        SinkExportFlags = malloc(sizeof(bool)*SinkExportFlagsMaxAllocated*(NProcs-1)+1);
        ActiveIndexList = malloc(sizeof(int)*SinkExportFlagsMaxAllocated);
    }
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> O.K.

    // Find active sink particles and clear corresponding export_flags.
    int NActives = 0;
    for(int i=0;i<Pall.Nsink;i++){
        if(PsinkActive(i)){
            ActiveIndexList[NActives] = i;
            for(int k=0;k<NProcs-1;k++)
                SinkExportFlags[i][k] = 0;
            NActives ++;
        }
    }

    // Counters for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    // Export structure.
    struct StructSinkExport *SinkExportSend[NProcs];

    // Pick up sink particles who need to be exported.
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckSinkExportFlags(i,NProcs,SinkExportFlags);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructSinkExport),i);
        SinkExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        for(int k=0;k<NActives;k++){
            int leaf = ActiveIndexList[k];
            if(SinkExportFlags[leaf][i]){
                SinkExportSend[i][NExport].Index  = leaf;
                SinkExportSend[i][NExport].AccretionRadius = Psink[leaf]->AccretionRadius;
                SinkExportSend[i][NExport].Mass   = PsinkMass(leaf);
                SinkExportSend[i][NExport].Pos[0] = PsinkPosP(leaf)[0];
                SinkExportSend[i][NExport].Pos[1] = PsinkPosP(leaf)[1];
                SinkExportSend[i][NExport].Pos[2] = PsinkPosP(leaf)[2];

                SinkExportSend[i][NExport].Vel[0] = PsinkBody(leaf)->Vel[0];
                SinkExportSend[i][NExport].Vel[1] = PsinkBody(leaf)->Vel[1];
                SinkExportSend[i][NExport].Vel[2] = PsinkBody(leaf)->Vel[2];
                SinkExportSend[i][NExport].dt_localmin = Psink[leaf]->dt_localmin;
                NExport ++;
            }
        }
        assert(NExport == NExportThisTime[i]);
    }

    // Exchange the number of sink particles who will be exported.
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }

    // Prepare buffers.
    struct StructSinkExport *SinkExportRecv;
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructSinkExport));
    SinkExportRecv = BufferExportRecv;

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    // Start asynchronous data communication.
    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkExportSend[i],
            NExportThisTime[i]*sizeof(struct StructSinkExport),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(SinkExportRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructSinkExport),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);


    // Prepare sink flag arrays for hydro particles.
    struct StructHydroSinkFlags *HydroSinkFlags;
    HydroSinkFlags = malloc(sizeof(struct StructHydroSinkFlags)*Pall.Nhydro);
    for(int i=0;i<Pall.Nhydro;i++)
        HydroSinkFlags[i].TargetSinkID = ThisNode;

    SinkHydroMergingCheckForLocalData(NActives,ActiveIndexList,HydroSinkFlags); 

    // MPI wait
    double TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.SinkCommThisStep += GetElapsedTime()-TimeComm;

    SinkHydroMergingCheckForImportedData(NProcs,NImportThisTime,SinkExportRecv,HydroSinkFlags);

    UpdateSinkLocalTimeStep(NExportThisTime,NImportThisTime,NProcs,SinkExportSend,SinkExportRecv);

    // Counters for export and import. 
    int NExportSinkThisTime[NProcs];
    int NImportSinkThisTime[NProcs];
    memset(NExportSinkThisTime,0,sizeof(int)*NProcs);
    memset(NImportSinkThisTime,0,sizeof(int)*NProcs);

    int NLocalVictims = 0;
    int NExportSinkThisTimeAll = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroSinkFlags[i].TargetSinkID != NONE){
            if(HydroSinkFlags[i].NodeID == ThisNode){
                NLocalVictims ++;
            } else {
                NExportSinkThisTime[HydroSinkFlags[i].NodeID] ++;
                NExportSinkThisTimeAll ++;
            }
        }
    }
    // dprintlmpi(NLocalVictims);
    // dprintlmpi(NExportSinkThisTimeAll);

    struct StructSinkImport *SinkImportSend;
    CheckSizeofBufferImportSend(NExportSinkThisTimeAll,sizeof(struct StructSinkImport));
    SinkImportSend = BufferImportSend;

    // prepare headers and counters.
    int Headers[MPIGetNumProcs()];
    int Counters[MPIGetNumProcs()];
    for(int i=0;i<MPIGetNumProcs();i++) Counters[i] = 0;
    Headers[0] = 0;
    for(int i=1;i<MPIGetNumProcs();i++)
        Headers[i] = Headers[i-1]+ NExportSinkThisTime[i-1];

    int check_vanished_export = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroSinkFlags[i].TargetSinkID != NONE){
            if(HydroSinkFlags[i].NodeID != ThisNode){
                int Index = Headers[HydroSinkFlags[i].NodeID] + Counters[HydroSinkFlags[i].NodeID];

                SinkImportSend[Index].Pos[0] = PhydroPos(i)[0];
                SinkImportSend[Index].Pos[1] = PhydroPos(i)[1];
                SinkImportSend[Index].Pos[2] = PhydroPos(i)[2];

                SinkImportSend[Index].Vel[0] = PhydroBody(i)->Vel[0];
                SinkImportSend[Index].Vel[1] = PhydroBody(i)->Vel[1];
                SinkImportSend[Index].Vel[2] = PhydroBody(i)->Vel[2];

                SinkImportSend[Index].Mass = PhydroMass(i);
                SinkImportSend[Index].TargetSinkIndex =  HydroSinkFlags[i].TargetSinkID;
                Phydro[i]->Use = OFF;
                PhydroBody(i)->Use = OFF;
                check_vanished_export ++;

                Counters[HydroSinkFlags[i].NodeID] ++;
            }
        }
    }

    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportSinkThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
            NImportSinkThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportSinkThisTime[i];
    }

    struct StructSinkImport *SinkImportRecv[NProcs];
    CheckSizeofBufferImportRecv(NImportSinkThisTime,sizeof(struct StructSinkImport));
    for(int i=0;i<NProcs-1;i++)
        SinkImportRecv[i] = BufferImportRecv[i];

    NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkImportSend+NImport,
            NExportSinkThisTime[i]*sizeof(struct StructSinkImport),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_SINK_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(SinkImportRecv[i],
            NImportSinkThisTime[i]*sizeof(struct StructSinkImport),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SINK_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NExportSinkThisTime[i];
    }
    //assert(NImport == NImportAll);


    ///////////////////////////////////////////
    // Sink operations for the local domain. //
    ///////////////////////////////////////////
    int check_vanished_local = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroSinkFlags[i].TargetSinkID != NONE){
            if(HydroSinkFlags[i].NodeID == ThisNode){
                ParticleAbsorption(HydroSinkFlags[i].TargetSinkID,
                        PhydroPos(i),PhydroVel(i),PhydroMass(i));
                Phydro[i]->Use = OFF;
                PhydroBody(i)->Use = OFF;
                check_vanished_local ++;
            }
        }
    }

    TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.SinkCommThisStep += GetElapsedTime()-TimeComm;

    ////////////////////////////////////////
    // Sink operations for imported data. //
    ////////////////////////////////////////
    for(int i=0;i<NProcs-1;i++){
        if(NImportSinkThisTime[i]>0){
            for(int k=0;k<NImportSinkThisTime[i];k++){
                ParticleAbsorption(SinkImportRecv[i][k].TargetSinkIndex,
                    SinkImportRecv[i][k].Pos,SinkImportRecv[i][k].Vel,SinkImportRecv[i][k].Mass);
            }
        }
    }

    int TotalVictims = NLocalVictims + NExportSinkThisTimeAll;

    int GlobalTotalVictims;
    MPI_Allreduce(&TotalVictims,&GlobalTotalVictims,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(GlobalTotalVictims > 0){
#if ShowLog //{
        fprintf(stderr,"[%02d] Local and importerd victims = %d and %d\n",
                MPIGetMyID(),NLocalVictims,NExportSinkThisTimeAll);
        fflush(stderr);
#endif // ShowLog //}
        Pall.Nhydro -= TotalVictims;
        Pall.Ntotal -= TotalVictims;
        ReConnectPointers();
        UpdateTotalNumber();
        UpdateTotalActiveNumber(); // What is this?
    }


    free(HydroSinkFlags);

    return ;
}
#endif

void SinkSinkMerging(void);
void SinkSinkMergingExported(void){
    SinkSinkMerging();
    return;
}

/* Friend-of-Friend */
int SinkFOFGroups;
int SinkFOFSize;
int SinkFOFCatalogSize;

struct StructSinkFOF{
    int Next;
    int Head;
    int Tail;
} *SinkFOF;

struct StructSinkFOFCatalog{
    int Number;
    int Head;
    // double Mass;
    // double Pos[3];
    // double Vel[3];
} *SinkFOFCatalog;


static void AllocateStructSinkFOF(const int Number){

    SinkFOF = malloc(MAX(sizeof(struct StructSinkFOF)*Number,1));
    SinkFOFSize = Number;

    return;
}

static void ReleaseStructSinkFOF(void){
    free(SinkFOF);
    return ;
}

static void StretchSinkFOFStructure(const int Number){

    SinkFOF = realloc(SinkFOF,MAX(sizeof(struct StructSinkFOF)*Number,1));
    SinkFOFSize = Number;

    return ;
}

static void InitializeStructSinkFOF(void){

    for(int i=0;i<SinkFOFSize;i++){
        SinkFOF[i].Next = NONE;
        SinkFOF[i].Head = SinkFOF[i].Tail = i;
    }

    return ;
}

static void ReleaseStructSinkFOFCatalog(void){

    if(SinkFOFCatalogSize > 0)
        free(SinkFOFCatalog);
    SinkFOFCatalogSize = 0;
    return ;
}

static void AllocateStructSinkFOFCatalog(const int NFOFGroups){

    if(SinkFOFCatalogSize > 0)
        ReleaseStructSinkFOFCatalog();

    SinkFOFCatalog = malloc(sizeof(struct StructSinkFOFCatalog)*NFOFGroups);
    if(SinkFOFCatalog == NULL){
        fprintf(stderr,"Allcoation of FOFCatalog is a failture.\n");
        fprintf(stderr,"In line %d, the function name is %s, the file name is %s\n",
            __LINE__,__FUNCTION__,__FILE__);
        exit(EXIT_FAILURE);
    }
    SinkFOFCatalogSize = NFOFGroups;
    return;
}


/* Friend-of-Friend */

struct StructParticleDataForActiveSinkFOF{
    double Pos[3];
    double Vel[3];
    double Mass;
    double MergingDistance;
    double dt_localmin;
    double AccretionMass;
    double AccretionMassGas;
    double AccretionMassStar;
    double AccretionMassToBH;
    unsigned long int GlobalID;
    int Leaf;
    bool Local;
};

static double ReturnTotalEnergySinkSinkBoundCheck( 
        struct StructParticleDataForActiveSinkFOF P1,
        struct StructParticleDataForActiveSinkFOF P2,
        const double distance){

    double velocity2 = DISTANCE2(P1.Vel,P2.Vel);
    double mu = P1.Mass*P2.Mass/(P1.Mass+P2.Mass);
    double E = 0.5*mu*velocity2 - Pall.GravConst*SQ(mu)/distance;

    /*
    if(E < 0){
        // fprintf(stderr,"Bound check clear %g %g %g\n",E,0.5*mu*velocity2,-Pall.GravConst*SQ(mu)/distance);
    } else {
        // fprintf(stderr,"Bound check failure %g %g %g\n",E,0.5*mu*velocity2,-Pall.GravConst*SQ(mu)/distance);
    }
    */
    return E;
}


static void SinkLocalFOFLink(struct StructParticleDataForActiveSinkFOF ParticleDataForActiveFOF[restrict]){

#if ShowLog //{
    fprintf(stderr,"Start %s\n",__FUNCTION__);
    fprintf(stderr,"SinkFOFSize %d\n",SinkFOFSize);
    fflush(NULL);
#endif // ShowLog //}

    for(int i=0;i<SinkFOFSize;i++){
        double Pos[3] = {ParticleDataForActiveFOF[i].Pos[0],
                         ParticleDataForActiveFOF[i].Pos[1],
                         ParticleDataForActiveFOF[i].Pos[2]};
        double Radiusi2 = SQ(ParticleDataForActiveFOF[i].MergingDistance);

        for(int k=0;k<SinkFOFSize;k++){
            if(i == k) continue;

            double distance2 = DISTANCE2(ParticleDataForActiveFOF[i].Pos,ParticleDataForActiveFOF[k].Pos);
            //if( (distance2 > SQ(ParticleDataForActiveFOF[i].Radius))&&
            if( (distance2 > Radiusi2)&&
                (distance2 > SQ(ParticleDataForActiveFOF[k].MergingDistance)) ) continue;

 
#ifdef USE_SINKSINK_MERGING_BOUND_CONDITION
            // check bound condition.
            double TotalEnergy = ReturnTotalEnergySinkSinkBoundCheck(
                    ParticleDataForActiveFOF[i],ParticleDataForActiveFOF[k],
                    sqrt(distance2));

            if(TotalEnergy > 0.e0) continue;
#endif //USE_SINKSINK_MERGING_BOUND_CONDITION

            if(SinkFOF[k].Head != SinkFOF[i].Head){
                int new_head = SinkFOF[i].Head;
                int new_tail = SinkFOF[k].Tail;
                int prev_tail = SinkFOF[i].Tail;

                SinkFOF[prev_tail].Next = SinkFOF[k].Head; // some problems in this line.

                int current_leaf = SinkFOF[i].Head;
                do{
                    SinkFOF[current_leaf].Head = new_head;
                    SinkFOF[current_leaf].Tail = new_tail;
                    current_leaf = SinkFOF[current_leaf].Next;
                } while(current_leaf != NONE);
            }
        }
    }
#if ShowLog //{
    fprintf(stderr,"End %s\n",__FUNCTION__);
    fflush(NULL);
#endif // ShowLog //}

    return ;
}

static void AdditionalSinkFOFLink(const int StartID, struct StructParticleDataForActiveSinkFOF ParticleDataForActiveFOF[restrict]){

    for(int i=StartID;i<SinkFOFSize;i++){
        for(int k=0;k<SinkFOFSize;k++){
            if(i == k) continue;

            double distance2 = DISTANCE2(ParticleDataForActiveFOF[i].Pos,ParticleDataForActiveFOF[k].Pos);
            if(distance2 > 
                (SQ(ParticleDataForActiveFOF[i].MergingDistance) + SQ(ParticleDataForActiveFOF[k].MergingDistance)+
                    2*ParticleDataForActiveFOF[i].MergingDistance*ParticleDataForActiveFOF[k].MergingDistance)) continue;

            // check bound condition.
            double TotalEnergy = ReturnTotalEnergySinkSinkBoundCheck(
                    ParticleDataForActiveFOF[i],ParticleDataForActiveFOF[k],
                    sqrt(distance2));
            if(TotalEnergy > 0.e0) continue;

            if(SinkFOF[k].Head != SinkFOF[i].Head){
                int new_head = SinkFOF[i].Head;
                int new_tail = SinkFOF[k].Tail;
                int prev_tail = SinkFOF[i].Tail;

                SinkFOF[prev_tail].Next = SinkFOF[k].Head; // some problems in this line.

                int current_leaf = SinkFOF[i].Head;
                do{
                    SinkFOF[current_leaf].Head = new_head;
                    SinkFOF[current_leaf].Tail = new_tail;
                    current_leaf = SinkFOF[current_leaf].Next;
                } while(current_leaf != NONE);
            }
        }
    }

    return ;
}

static void MakeSinkFOFCatalog(void){

    SinkFOFGroups = 0;
    for(int i=0;i<SinkFOFSize;i++){
        if(i == SinkFOF[i].Head){
            SinkFOFGroups ++;
        }
    }
    if(SinkFOFGroups == 0) return;

    ReleaseStructSinkFOFCatalog();
    AllocateStructSinkFOFCatalog(SinkFOFGroups);
#if ShowLog //{
    fprintf(stderr,"[%02d] SinkFOFSize = %d:%s:%d\n",MPIGetMyID(),SinkFOFSize,__FUNCTION__,__LINE__);
    fflush(NULL);
#endif // ShowLog //}

    int counter = 0;
    for(int i=0;i<SinkFOFSize;i++){
        if(i == SinkFOF[i].Head){
            SinkFOFCatalog[counter].Head = SinkFOF[i].Head;
            SinkFOFCatalog[counter].Number = 0;
            int current_leaf = SinkFOF[i].Head;
            do{
                SinkFOFCatalog[counter].Number ++;
                current_leaf = SinkFOF[current_leaf].Next;
            } while(current_leaf != NONE);
#if ShowLog //{
            fprintf(stderr,"[%02d] !Number of FOF  = %d:%s:%d\n",
                    MPIGetMyID(),SinkFOFCatalog[counter].Number,__FUNCTION__,__LINE__);
            fflush(NULL);
#endif // ShowLog //}
            counter ++;
        }
    }

#if ShowLog //{
    fprintf(stderr,"[%02d] Number of FOF[0]  = %d:%s:%d\n",MPIGetMyID(),SinkFOFCatalog[0].Number,__FUNCTION__,__LINE__);
    fprintf(stderr,"[%02d] Number of FOF groups = %d:%s:%d\n",MPIGetMyID(),SinkFOFGroups,__FUNCTION__,__LINE__);
    fflush(NULL);
#endif // ShowLog //}

    return ;
}

/*
 * This function checks and returns number of sink particles which have to
 * export. The check is done for every FOF clump to the "Active Sink Edges". 
 */
static inline int __attribute__((always_inline)) CheckSinkSinkExportFlagsForActiveEdges(const int Index, const int NProcs, bool SinkExportFlags[][NProcs-1], struct StructParticleDataForActiveSinkFOF ParticleDataForActiveFOF[restrict]){

    int NodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;

#if ShowLog //{
    fprintf(stderr,"Number of groups = %d %s:%d\n",SinkFOFGroups,__FUNCTION__,__LINE__);
    fflush(NULL);
#endif // ShowLog //}

    for(int i=0;i<SinkFOFGroups;i++){
        bool export_flag = false;
        int current_leaf = SinkFOFCatalog[i].Head;
        do{
            if(CheckDomainOverlapForActiveSinkEdges(ParticleDataForActiveFOF[current_leaf].Pos,
                        ParticleDataForActiveFOF[current_leaf].MergingDistance,NodeID)){
                export_flag = true;
            }
            current_leaf = SinkFOF[current_leaf].Next;
        } while(current_leaf != NONE);

        //if(export_flag == true){
            current_leaf = SinkFOFCatalog[i].Head;
            do{
                SinkExportFlags[ParticleDataForActiveFOF[current_leaf].Leaf][Index] = ON;
                current_leaf = SinkFOF[current_leaf].Next;
            } while(current_leaf != NONE);
            NExport += SinkFOFCatalog[i].Number;
        //}
    }

    return NExport;
}


static int SinkSinkMergingEngine(struct StructParticleDataForActiveSinkFOF ParticleDataForActiveFOF[restrict]){

    // search all fof clumps.

    int check_vanished_local = 0;
    for(int i=0;i<SinkFOFGroups;i++){
        if(SinkFOFCatalog[i].Number > 1){

            int current_leaf = SinkFOFCatalog[i].Head;   
            int GlobalID_minimum = ParticleDataForActiveFOF[current_leaf].GlobalID;
            int Leaf_minimum = ParticleDataForActiveFOF[current_leaf].Leaf;
            int domain_local = ParticleDataForActiveFOF[current_leaf].Local;
            
            double PxMax = ParticleDataForActiveFOF[current_leaf].Mass*
                            ParticleDataForActiveFOF[current_leaf].Pos[0];
            do{
                if(PxMax < ParticleDataForActiveFOF[current_leaf].Mass*
                            ParticleDataForActiveFOF[current_leaf].Pos[0]){

                    PxMax = ParticleDataForActiveFOF[current_leaf].Mass*
                            ParticleDataForActiveFOF[current_leaf].Pos[0];
                    Leaf_minimum = ParticleDataForActiveFOF[current_leaf].Leaf;
                    domain_local = ParticleDataForActiveFOF[current_leaf].Local;
                }
                current_leaf = SinkFOF[current_leaf].Next;
            }while(current_leaf != NONE);

#if 0
            do{
                if(GlobalID_minimum > ParticleDataForActiveFOF[current_leaf].GlobalID){
                    GlobalID_minimum = ParticleDataForActiveFOF[current_leaf].GlobalID;
                    Leaf_minimum = ParticleDataForActiveFOF[current_leaf].Leaf;
                    domain_local = ParticleDataForActiveFOF[current_leaf].Local;
                }
                current_leaf = SinkFOF[current_leaf].Next;
            }while(current_leaf != NONE);

            if(ShowLog){
                fprintf(stderr,"[%02d] Number of Members[%d] = %d, %d\n",
                        MPIGetMyID(),i,SinkFOFCatalog[i].Number,GlobalID_minimum);
                fflush(NULL);
            }
#endif


            if(domain_local == true){
                // Sink all particles.

                int current_leaf = SinkFOFCatalog[i].Head;   
                double COM[3] = {0.e0,0.e0,0.e0};
                double VCOM[3] = {0.e0,0.e0,0.e0};
                double CMASS = 0.e0;
                double dt_localmin = ParticleDataForActiveFOF[current_leaf].dt_localmin;
                double AccretionMass = 0.e0;
                double AccretionMassGas = 0.e0;
                double AccretionMassStar = 0.e0;
                double AccretionMassToBH = 0.e0;
                do{
                    CMASS += ParticleDataForActiveFOF[current_leaf].Mass;
                    COM[0] += ParticleDataForActiveFOF[current_leaf].Mass*ParticleDataForActiveFOF[current_leaf].Pos[0];
                    COM[1] += ParticleDataForActiveFOF[current_leaf].Mass*ParticleDataForActiveFOF[current_leaf].Pos[1];
                    COM[2] += ParticleDataForActiveFOF[current_leaf].Mass*ParticleDataForActiveFOF[current_leaf].Pos[2];

                    VCOM[0] += ParticleDataForActiveFOF[current_leaf].Mass*ParticleDataForActiveFOF[current_leaf].Vel[0];
                    VCOM[1] += ParticleDataForActiveFOF[current_leaf].Mass*ParticleDataForActiveFOF[current_leaf].Vel[1];
                    VCOM[2] += ParticleDataForActiveFOF[current_leaf].Mass*ParticleDataForActiveFOF[current_leaf].Vel[2];
                    dt_localmin = fmin(dt_localmin,ParticleDataForActiveFOF[current_leaf].dt_localmin);
                    AccretionMass += ParticleDataForActiveFOF[current_leaf].AccretionMass;
                    AccretionMassGas += ParticleDataForActiveFOF[current_leaf].AccretionMassGas;
                    AccretionMassStar += ParticleDataForActiveFOF[current_leaf].AccretionMassStar;
                    AccretionMassToBH += ParticleDataForActiveFOF[current_leaf].AccretionMassToBH;
#if 0
                    if((ParticleDataForActiveFOF[current_leaf].Local == true)&&
                            (ParticleDataForActiveFOF[current_leaf].GlobalID != GlobalID_minimum))
#else
                    //if((ParticleDataForActiveFOF[current_leaf].Local == true)&&(counter > 0))
                    double Px = ParticleDataForActiveFOF[current_leaf].Mass*
                                ParticleDataForActiveFOF[current_leaf].Pos[0];
                    if((ParticleDataForActiveFOF[current_leaf].Local == true)&&(PxMax != Px))
#endif
                    {

                        /*
                        fprintf(stderr,"%d(%ld) -> %d(%ld), dt %g -> %g | %g\n",
                                ParticleDataForActiveFOF[current_leaf].Local,
                                ParticleDataForActiveFOF[current_leaf].GlobalID,
                                Leaf_minimum,
                                PsinkBody(Leaf_minimum)->GlobalID,
                                PsinkBody(ParticleDataForActiveFOF[current_leaf].Local)->dt,
                                PsinkBody(Leaf_minimum)->dt,
                                Psink[Leaf_minimum]->dt_localmin);
                        */
                        int index = ParticleDataForActiveFOF[current_leaf].Leaf;
                        Psink[index]->Use = false;
                        PsinkBody(index)->Use = false;
                        fprintf(stderr,"Vanish![%02d] %d %d:%s:%s\n",MPIGetMyID(),check_vanished_local,
                                __LINE__,__FUNCTION__,__FILE__);
                        check_vanished_local ++;
                    }
                    current_leaf = SinkFOF[current_leaf].Next;
                }while(current_leaf != NONE);

                Psink[Leaf_minimum]->dt_localmin = fmin(Psink[Leaf_minimum]->dt_localmin,dt_localmin);

                double iCMASS = 1.e0/CMASS;       
                for(int k=0;k<3;k++){
                    PsinkBody(Leaf_minimum)->Pos[k] = PsinkBody(Leaf_minimum)->PosP[k] = 
                        Psink[Leaf_minimum]->PosP[k] = COM[k]*iCMASS;
                    PsinkBody(Leaf_minimum)->Vel[k] = 
                    Psink[Leaf_minimum]->VelP[k] = VCOM[k]*iCMASS;
                }
                PsinkMass(Leaf_minimum) = CMASS;
                Psink[Leaf_minimum]->AccretionMass = AccretionMass;
                Psink[Leaf_minimum]->AccretionMassGas = AccretionMassGas;
                Psink[Leaf_minimum]->AccretionMassStar = AccretionMassStar;
                Psink[Leaf_minimum]->AccretionMassToBH = AccretionMassToBH;

                Psink[Leaf_minimum]->NumberofAbsorbedParticles += SinkFOFCatalog[i].Number-1;
            } else {
                int current_leaf = SinkFOFCatalog[i].Head;   
                do{
                    if(ParticleDataForActiveFOF[current_leaf].Local == true){
                        int index = ParticleDataForActiveFOF[current_leaf].Leaf;
                        Psink[index]->Use = false;
                        PsinkBody(index)->Use = false;

                        fprintf(stderr,"Vanish![%02d] %d %d:%s:%s\n",MPIGetMyID(),check_vanished_local,
                                __LINE__,__FUNCTION__,__FILE__);
                        check_vanished_local ++;
                    }
                    current_leaf = SinkFOF[current_leaf].Next;
                }while(current_leaf != NONE);
            }
        }
    }

#if ShowLog //{
    fprintf(stderr,"Number of Groups = %d\n",SinkFOFGroups);
    fprintf(stderr,"Number of vanished local Sinks = %d, Pall.Nsinks = %ld\n",check_vanished_local,Pall.Nsink);
    fflush(NULL);
#endif // ShowLog //}

    return check_vanished_local;

}


void SinkSinkMerging(void){

    if(Pall.NActivesSink_t < 2) 
        return;
#if USE_SINK_TIMESTEP_LIMITER
    static const double InvSinkTimeStepLimiterFactor = 1.0/((double)(1<<SINK_TIMESTEP_MAX_K_LOCAL));
#endif

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    // Initialize data arrays.
    static int SinkExportFlagsMaxAllocated = 0;
    static bool (*SinkExportFlags)[NProcs-1];
    static int *ActiveIndexList;
    if(SinkExportFlagsMaxAllocated < Pall.Nsink){
        if(SinkExportFlagsMaxAllocated > 0){
            free(ActiveIndexList);
            free(SinkExportFlags);
        }
        SinkExportFlagsMaxAllocated = MAX((int)(ForAngelsShare*Pall.Nsink),NAdditionUnit);
        SinkExportFlags = malloc(sizeof(bool)*SinkExportFlagsMaxAllocated*(NProcs-1)+1);
        ActiveIndexList = malloc(sizeof(int)*SinkExportFlagsMaxAllocated);
    }

    // Find active sink particles and clear corresponding export_flags.
    int NActives = 0;
    for(int i=0;i<Pall.Nsink;i++){
        if(PsinkActive(i)){
            ActiveIndexList[NActives] = i;
            for(int k=0;k<NProcs-1;k++)
                SinkExportFlags[i][k] = 0;
            NActives ++;
        }
    }

    // Do local FOF.
    // FOF Start!
    static struct StructParticleDataForActiveSinkFOF *ParticleDataForActiveFOF;

    ParticleDataForActiveFOF = malloc(sizeof(struct StructParticleDataForActiveSinkFOF)*(Pall.NActivesSink)+1);
    for(int i=0;i<Pall.NActivesSink;i++){
        int leaf = ActiveIndexList[i];
        // Do second kick before input the temporarly array.
        ParticleDataForActiveFOF[i].Pos[0] = PsinkBody(leaf)->Pos[0];
        ParticleDataForActiveFOF[i].Pos[1] = PsinkBody(leaf)->Pos[1];
        ParticleDataForActiveFOF[i].Pos[2] = PsinkBody(leaf)->Pos[2];

        double dt_half = 0.5*PsinkBody(leaf)->dt;
        ParticleDataForActiveFOF[i].Vel[0] = PsinkBody(leaf)->Vel[0];
        ParticleDataForActiveFOF[i].Vel[1] = PsinkBody(leaf)->Vel[1];
        ParticleDataForActiveFOF[i].Vel[2] = PsinkBody(leaf)->Vel[2];

        ParticleDataForActiveFOF[i].Mass = PsinkBody(leaf)->Mass;
        ParticleDataForActiveFOF[i].MergingDistance = Psink[leaf]->MergingDistance;
#if USE_SINK_TIMESTEP_LIMITER
        ParticleDataForActiveFOF[i].dt_localmin = fmin(Psink[leaf]->dt_localmin,InvSinkTimeStepLimiterFactor*PsinkBody(leaf)->dt);
#else
        ParticleDataForActiveFOF[i].dt_localmin = Psink[leaf]->dt_localmin;
#endif
        ParticleDataForActiveFOF[i].AccretionMass = Psink[leaf]->AccretionMass;
        ParticleDataForActiveFOF[i].AccretionMassGas = Psink[leaf]->AccretionMassGas;
        ParticleDataForActiveFOF[i].AccretionMassStar = Psink[leaf]->AccretionMassStar;
        ParticleDataForActiveFOF[i].AccretionMassToBH = Psink[leaf]->AccretionMassToBH;

        ParticleDataForActiveFOF[i].GlobalID = PsinkBody(leaf)->GlobalID;
        ParticleDataForActiveFOF[i].Leaf = leaf;
        ParticleDataForActiveFOF[i].Local = true;
    }

    AllocateStructSinkFOF(Pall.NActivesSink);
    InitializeStructSinkFOF();

    // Do FOF for active sinks.
    SinkLocalFOFLink(ParticleDataForActiveFOF);
    MakeSinkFOFCatalog();

    // Counters for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    // Data structure for Export.
    struct StructParticleDataForActiveSinkFOF *SinkSinkExportSend[NProcs];

    // Pick up sink particles who need to be exported.
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckSinkSinkExportFlagsForActiveEdges(i,NProcs,SinkExportFlags,
                ParticleDataForActiveFOF);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructParticleDataForActiveSinkFOF),i);
        SinkSinkExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        for(int k=0;k<NActives;k++){
            int leaf = ActiveIndexList[k];
            if(SinkExportFlags[leaf][i]){
                SinkSinkExportSend[i][NExport] = ParticleDataForActiveFOF[k];
                SinkSinkExportSend[i][NExport].Local = false;
                NExport ++;
            }
        }
        assert(NExport == NExportThisTime[i]);
    }

    // Exchange the number of sink particles who will be exported.
    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SINK_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SINK_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }

    // realloc //
    ParticleDataForActiveFOF = realloc(ParticleDataForActiveFOF,
            MAX(sizeof(struct StructParticleDataForActiveSinkFOF)*(Pall.NActivesSink+NImportAll),1));

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    // Start asynchronous data communication.
    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(SinkSinkExportSend[i],
            NExportThisTime[i]*sizeof(struct StructParticleDataForActiveSinkFOF),MPI_BYTE,
                CommunicationTable[i].SendRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);

        MPI_Irecv(ParticleDataForActiveFOF+Pall.NActivesSink+NImport,
            NImportThisTime[i]*sizeof(struct StructParticleDataForActiveSinkFOF),MPI_BYTE,
                CommunicationTable[i].RecvRank,TAG_SINK_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);

    // MPI wait
    double TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.SinkCommThisStep += GetElapsedTime()-TimeComm;


    // Strach FOF structures.
    StretchSinkFOFStructure(Pall.NActivesSink+NImportAll);
    InitializeStructSinkFOF();


    // Additional FOF.
    //AdditionalSinkFOFLink(Pall.NActivesSink,ParticleDataForActiveFOF);
    SinkLocalFOFLink(ParticleDataForActiveFOF);
    MakeSinkFOFCatalog();

    // check FOF clumps; if the minimum GlobalID is in this node, sink alls.
    int TotalVictims = SinkSinkMergingEngine(ParticleDataForActiveFOF); // calc
    int GlobalTotalVictims;
    MPI_Allreduce(&TotalVictims,&GlobalTotalVictims,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(GlobalTotalVictims > 0){
#if ShowLog //{
        fprintf(stderr,"[%02d] Local and importerd victims = %d and %d:%s%d\n",
                MPIGetMyID(),TotalVictims,GlobalTotalVictims,__FUNCTION__,__LINE__);
        fprintf(stderr,"[%02d] Number of Active Sinks = %ld -> %ld:%s%d\n",
                MPIGetMyID(),Pall.NActivesSink,Pall.NActivesSink-TotalVictims,__FUNCTION__,__LINE__);
        fflush(stderr);
#endif // ShowLog //}
        Pall.Nsink -= TotalVictims;
        Pall.Ntotal -= TotalVictims;
        ReConnectPointers();
        UpdateTotalNumber();
        UpdateTotalActiveNumber(); // What is this?
    }

    free(ParticleDataForActiveFOF);
    ReleaseStructSinkFOF();
    ReleaseStructSinkFOFCatalog();

    return ;
}


/*
 * This function generates new sink particles. We do not care the internal
 * energy of hydro particles.
 */
static void MakeSinkParticles(void){

//#warning This function is under tests.

    int NewSinkParticles = 0; 
    for(int i=0;i<Pall.Nhydro;i++){ // Search all hydro particles.
        //if(PhydroActive(i)){ //
        if(AllowSinkConversion[i] == false) continue;
        if((Phydro[i]->Active == true)&&(PhydroBody(i)->Active == true)){ //
            // if the particle i satisfies the threshold condition(s), convert the
            // particle into a sink particle.
            if(CheckSinkConditions(i)){
                // Allocate a sink particle in the sink array.  If the array does
                // not have any free element, the array should be strach.
                StructPsinkptr Psk = ReturnEmptySinkStructurePointer();

                // Setup sink partcile.
                Psk->Use = ON;
                Psk->ParentGlobalID = PhydroBody(i)->GlobalID;
                Psk->Body = PhydroBody(i);
                Psk->Body->Baryon = (void*)Psk;
                Psk->Body->Type = TypeSink;
                Psk->FormationTime = Pall.TCurrent;
                Psk->dt_localmin = PhydroBody(i)->dt;

                Psk->PosP[0] = Phydro[i]->PosP[0];
                Psk->PosP[1] = Phydro[i]->PosP[1];
                Psk->PosP[2] = Phydro[i]->PosP[2];

                Psk->VelP[0] = Phydro[i]->VelP[0];
                Psk->VelP[1] = Phydro[i]->VelP[1];
                Psk->VelP[2] = Phydro[i]->VelP[2];

                // follow Z
                Psk->Z   = Phydro[i]->Z;
                Psk->ZII = Phydro[i]->ZII;
                Psk->ZIa = Phydro[i]->ZIa;

#ifdef USE_CELIB //{
                for(int k=0;k<CELibYield_Number;k++){
                    Psk->Elements[k] = Phydro[i]->Elements[k];
                }
#endif // USE_CELIB //}


                Psk->AccretionRadius = SINKHYDRO_ACCRETION_RADIUS/Pall.UnitLength;
                Psk->MergingDistance = SINKSINK_MERGING_DISTANCE/Pall.UnitLength;

                HydroRoot.Leaves[Phydro[i]->Leaf] *= -1;

                //Remove the hydro particle.
                Phydro[i]->Use = OFF;

                NewSinkParticles ++;
            }
        }
    }
    int GlobalNewSinkParticles;
    MPI_Allreduce(&NewSinkParticles,&GlobalNewSinkParticles,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(GlobalNewSinkParticles > 0){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"%d new sink particles are born in this step\n",GlobalNewSinkParticles);
        }
    }
    if(GlobalNewSinkParticles > 0){
        Pall.Nhydro -= NewSinkParticles;
        Pall.Nhydro_t -= GlobalNewSinkParticles;
        Pall.Nsink += NewSinkParticles;
        Pall.Nsink_t += GlobalNewSinkParticles;

        ReConnectPointers(); // Is this necessary?

        for(int i=0;i<HydroRoot.NumberofLeaves;i++)
            HydroRoot.Leaves[i] = NONE;

        for(int i=0;i<Pall.Nhydro;i++){
            int index = Phydro[i]->Leaf;
            NBCache[index].Leaf = i;
            HydroRoot.Leaves[index] = i;
        }
        //UpdateTotalActiveNumber();
    }
    /*
    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"New/GNew %d, %d, Nhydro/Nsink %ld, %ld\n",
            NewSinkParticles,GlobalNewSinkParticles,Pall.Nhydro,Pall.Nsink);
        fprintf(stderr,"Pall.Ntotal_t,Nhydro_t,Nsink_t = %ld, %ld, %ld\n",
            Pall.Ntotal_t,Pall.Nhydro_t,Pall.Nsink_t);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    */

    return ;
}

/*
 * When all of sink particle formation criteria are satisfied, this function
 * returns true in boolean. Otherwise this function returns false.
 */
static bool CheckSinkConditions(const int index){
    int flag = false;

    // When this hydro particle satisfies all of conditions in below,
    // the hydro particle becomes a sink particle.
    if(Phydro[index]->Rho > Pall.SinkThresholdDensity){
        flag = true;
    }

    return flag;
}


static void ReConnectHydroAndSinkPointers(void){

    /* For Body */
    if(Pall.Ntotal > PbodySize){
        PbodySize = (int)(ForAngelsShare*Pall.Ntotal);
        free(Pbody);
        Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
        //Pbody = realloc(Pbody,PbodySize*sizeof(StructPbodyptr));
    }
    int counter = 0;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == ON){
            Pbody[counter] = Pb;
            counter ++;
        }
    }
    if(counter != Pall.Ntotal){
        MPI_Finalize();
        fprintf(stderr,"Element Number of Body is not correct! (StarFormation.c)\n");
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nhydro = %ld, Pall.Ntotal_t = %ld\n",
                MPIGetMyID(),counter,Pall.Ntotal,Pall.Ntotal_t);
        exit(StarFormationElementNumberError);
    }


    /* For Hydro */
    if(Pall.Nhydro > PhydroSize){
        PhydroSize = (int)(ForAngelsShare*Pall.Nhydro);
        free(Phydro);
        Phydro = malloc(PhydroSize*sizeof(StructPhydroptr));
        //Phydro = realloc(Phydro,PhydroSize*sizeof(StructPhydroptr));
        //dprintlmpi(PhydroSize);
    }
    counter = 0;
    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next){
        if(Ph->Use == ON){
            Phydro[counter] = Ph;
            counter ++;
        }
    }
    if(counter != Pall.Nhydro){
        MPI_Finalize();
        fprintf(stderr,"Element Number of Hydro is not correct! MyID = %d (StarForamtion.c)\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nhydro = %ld, Pall.Nhydro_t = %ld, Pall.Ntotal_t = %ld\n",
                MPIGetMyID(),counter,Pall.Nhydro,Pall.Nhydro_t,Pall.Ntotal_t);
        exit(StarFormationElementNumberError);
    }

    /* For Stars */
    if(Pall.Nstars > PstarSize){
        PstarSize = (int)(ForAngelsShare*Pall.Nstars);
        free(Pstar);
        Pstar = malloc(PstarSize*sizeof(StructPstarptr));
        //Pstar = realloc(Pstar,PstarSize*sizeof(StructPstarptr));
        //dprintlmpi(PstarSize);
    }
    counter = 0;
    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next){
        if(Ps->Use == ON){
            Pstar[counter] = Ps;
            counter ++;
        }
    }
    if(counter != Pall.Nstars){
        MPI_Finalize();
        fprintf(stderr,"Element Number of Star is not correct! MyID = %d (StarForamtion.c)\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nstars = %ld, Pall.Nstar_t = %ld, Pall.Ntotal_t = %ld\n",
                MPIGetMyID(),counter,Pall.Nstars, Pall.Nstars_t,Pall.Ntotal_t);
        exit(StarFormationElementNumberError);
    }

    return;
}
