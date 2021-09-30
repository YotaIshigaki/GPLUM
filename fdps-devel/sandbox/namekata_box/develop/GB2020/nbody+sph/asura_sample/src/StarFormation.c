#include "config.h"
#include "PlantHydroTree.h"
#include "StructureOperation.h"
#include "ParallelOperation.h"
#ifdef USE_CELIB
#include "StellarFeedback.h"
#endif // USE_CELIB

extern double SNIINumber; // per solar mass.

static bool ConvertGasParticleIntoStarParticle(const int index);
static bool ConvertGasParticleIntoDeadStarParticle(const int index);
static bool SpawnStarParticle(const int index);
#ifdef USE_INSTANTANEOUS_STARFORMATION //{
static bool InstantaneousSFCondition(const int index);
#endif // USE_INSTANTANEOUS_STARFORMATION //}
static bool SFCondition(const int index);
static bool ForcibleSFCondition(const int index);
static void ReConnectBaryonPointers(void);

static double TimePrev = 0.e0;
static double MassPrev = 0.e0;
FILE *FpSFR;
char FileStarFormationRate[MaxCharactersInLine];
int TempNstars,TempNtotal,TempNhydro;

static inline __attribute__((always_inline)) void *SendBuf(const bool target, void *sendbuf){
    return target
        ? MPI_IN_PLACE
        : sendbuf;
}

#ifdef USE_STARFORMATION_BIGSTEP //{
extern int Kcurrent;

static inline __attribute__((always_inline)) double GetSFBigStepFactor(const int index){

    if(Phydro[index]->k_hydro - Kcurrent >= STARFORMATION_BIGSTEP_K_FACTOR){
        return 1.0;
    } else {
        int kdiff = Phydro[index]->k_hydro - Kcurrent; // 0,1,... STARFORMATION_BIGSTEP_K_FACTOR -1;
        int k_step = STARFORMATION_BIGSTEP_K_FACTOR - kdiff;
        return (double)(1<<k_step);
    }

    return ;
}
#endif // USE_STARFORMATION_BIGSTEP //}

void InitializeStarFormation(void){

    TimePrev = Pall.TCurrent;
    MassPrev = 0.e0;
    for(int i=0;i<Pall.Nstars;i++)
        MassPrev += PstarMass(i);

    MPI_Reduce(SendBuf(MPIGetMyID()==MPI_ROOT_RANK,&MassPrev),&MassPrev,1,MPI_DOUBLE,MPI_SUM,MPI_ROOT_RANK,MPI_COMM_WORLD);

    if(MPIGetMyID()==MPI_ROOT_RANK){
        sprintf(FileStarFormationRate,"./StarFormationRate.data");
        if(Pall.RunStatus == NewSimulation){
            FileOpen(FpSFR,FileStarFormationRate,"w");
        } else {
            FileOpen(FpSFR,FileStarFormationRate,"a");
        }
    }

    return;
}

void LogStarFormationRate(void){

#ifdef STARFORMATION //{
    double dt = Pall.TCurrent - TimePrev; 


#define Tave 3.e+7 // 30 Myr average.
    double TaveInUnit = Tave*YEAR_CGS/Pall.UnitTime;
    double MassCurrent[2] = {0.e0,0.e0};
    for(int i=0;i<Pall.Nstars;i++){
        MassCurrent[0] += PstarMass(i);
        if(TaveInUnit > Pall.TCurrent - Pstar[i]->FormationTime){
            MassCurrent[1] += PstarMass(i);
        }
    }

    MPI_Reduce(SendBuf(MPIGetMyID()==MPI_ROOT_RANK,MassCurrent),MassCurrent,2,MPI_DOUBLE,MPI_SUM,MPI_ROOT_RANK,MPI_COMM_WORLD);

    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(FpSFR,"%g %g %g %ld %g %g\n",
            Pall.UnitTime*Pall.TCurrent/YEAR_CGS,Pall.Redshift,Pall.UnitMass*(MassCurrent[0]-MassPrev)/MSUN_CGS,
                Pall.Nstars_t,Pall.UnitMass*(MassCurrent[0]-MassPrev)/MSUN_CGS/(Pall.UnitTime*dt/YEAR_CGS),
                (Pall.UnitMass*MassCurrent[1]/MSUN_CGS)/(Tave));
    }

    TimePrev = Pall.TCurrent;
    MassPrev = MassCurrent[0];

#endif //STARFORMATION //}

    return;
}

void StarFormation(void){

    // ResetLastStepSF();

#ifdef STARFORMATION
    if(Pall.NActivesHydro_t <= 1)
        return;

#ifdef USE_STARFORMATION_BIGSTEP //{
    if(Pall.kmax >= STARFORMATION_BIGSTEP_K_FACTOR){
        if((Pall.TStep)%(1<<STARFORMATION_BIGSTEP_K_FACTOR)!=0){
            if(MPIGetMyID() == MPI_ROOT_RANK){
                fprintf(stderr,"Skip SF\n");
                fflush(NULL);
            }
            return ;
        }
    }
#endif // USE_STARFORMATION_BIGSTEP //}



    double TimingResultThisRoutine = GetElapsedTime();

    TempNstars = TempNhydro = TempNtotal = 0;

    int CountNewStar = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroActive(i)){
#ifdef USE_INSTANTANEOUS_STARFORMATION //{
            if(InstantaneousSFCondition(i)){
                if(ConvertGasParticleIntoStarParticle(i) == true){
                    CountNewStar ++;
                }
            } else 
#endif // USE_INSTANTANEOUS_STARFORMATION //}
            if(SFCondition(i)){
#if (UseSFModelConvert)
                if(ConvertGasParticleIntoStarParticle(i) == true){
                    CountNewStar ++;
                }
#endif

#if (UseSFModelSpawn) //{
                if(PhydroMass(i) > 1.5*Phydro[i]->SpawnMass){
                    if(SpawnStarParticle(i) == true){
                        CountNewStar ++;
                    }
                }else{
                    if(ConvertGasParticleIntoStarParticle(i) == true){
                        CountNewStar ++;
                    }
                }
#endif // UseSFModelSpawn //}
            }
#ifdef COSMOLOGICAL_RUN //{
#ifdef USE_FORCIBLE_STARFORMATION //{
            else if(ForcibleSFCondition(i)){
                if(ConvertGasParticleIntoDeadStarParticle(i) == true){
                    CountNewStar ++;
                }
            }
#endif // USE_FORCIBLE_STARFORMATION //{
#endif // COSMOLOGICAL_RUN //}
        }
    }

    Pall.Nhydro += TempNhydro;
    Pall.NActivesHydro += TempNhydro;
    Pall.Nstars += TempNstars;
    Pall.NActivesStars += TempNstars;
    Pall.Ntotal += TempNtotal;
    Pall.NActives += TempNtotal;

    int Numbers[4] = {CountNewStar,Pall.Ntotal,Pall.Nstars,Pall.Nhydro};
    MPI_Allreduce(MPI_IN_PLACE,Numbers,4,MPI_INT,MPI_SUM,MPI_COMM_WORLD);


    int GlobalCountNewStar = Numbers[0];
    Pall.Ntotal_t = Numbers[1];
    Pall.Nstars_t = Numbers[2];
    Pall.Nhydro_t = Numbers[3];

#if 0
    if(GlobalCountNewStar > 0){
        TurnonLastStepSF();
    }
#endif

    if(CountNewStar > 0){
        ReConnectPointers();
        for(int i=0;i<HydroRoot.NumberofLeaves;i++)
            HydroRoot.Leaves[i] = NONE;

        for(int i=0;i<Pall.Nhydro;i++){
            int index = Phydro[i]->Leaf;
            NBCache[index].Leaf = i;
            HydroRoot.Leaves[index] = i;
        }
    }



    TimingResults.StarformationThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
#endif

}

#ifdef USE_INSTANTANEOUS_STARFORMATION //{
static bool InstantaneousSFCondition(const int index){

#if 0
#ifdef USE_MULTIPHASE_MODEL //{
    if(Phydro[index]->MultiphaseFlag)
        return false;
#endif // USE_MULTIPHASE_MODEL //}
#endif

    const double SFDensityCritIns = SFCONDITION_DENSITY_CRITERION
        *INSTANTANEOUS_STARFORMATION_DENSITY_THRESHOLD_FACTOR;

    if(0.76*Pall.ConvertNumberDensityToCGS*Phydro[index]->RhoPred > SFDensityCritIns){
        return true;
    } else {
        return false;
    }
}
#endif // USE_INSTANTANEOUS_STARFORMATION //}


static bool SFCondition(const int index){

#if 0
#ifdef USE_INSTANTANEOUS_STARFORMATION //{
    const double SFDensityCritIns = SFCONDITION_DENSITY_CRITERION
        *INSTANTANEOUS_STARFORMATION_DENSITY_THRESHOLD_FACTOR;
    if(0.76*Pall.ConvertNumberDensityToCGS*Phydro[index]->RhoPred > SFDensityCritIns)
        return true;
#endif // USE_INSTANTANEOUS_STARFORMATION //}
#endif


#ifdef USE_MULTIPHASE_MODEL //{
    if(Phydro[index]->MultiphaseFlag)
        return false;
#endif // USE_MULTIPHASE_MODEL //}


#if (CosmologicalRun) //{
    if(Phydro[index]->RhoPred < 200.0*Pall.RHOcZ)
        return false;
#endif // CosmologicalRun //}

#ifdef SFCONDITION_DENSITY //{
    if(0.76*Pall.ConvertNumberDensityToCGS*Phydro[index]->RhoPred < SFDensityCrit)
        return false;
#endif // SFCONDITION_DENSITY //}

#ifdef SFCONDITION_TEMPERATURE //{
    if((Pall.ConvertUtoT*Phydro[index]->UPred > SFTemperatureCrit))
        return false; 
#endif // SFCONDITION_TEMPERATURE //} 

#ifdef SFCONDITION_CONVERGING_FLOW //{
    if(Phydro[index]->Div > 0.0)
        return false;
#endif // SFCONDITION_CONVERGING_FLOW //}

// Note that this branch has no meaning now, since the star formation routine is
// evaluated after the cooling function where the DQheat is used and flushed.
#ifdef SFCONDITION_HAVE_SNE_HEAT //{
    if(Phydro[index]->DQheat > 0.0)
        return false;
#endif // SFCONDITION_HAVE_SNE_HEAT //}

#ifdef USE_WARMSTART //{
    double p = (-0.5*SQ((Pall.TCurrent*Pall.UnitTime/MEGAYEAR_CGS)/(100.0)));
    if(gsl_rng_uniform(RandomGenerator) < p){
        Phydro[index]->U = 
        Phydro[index]->UPred = 
        Pall.ConvertTtoU*5.e+7;
        return false;
    }
#endif //WARMSTART //}


    return true;
}


static bool ForcibleSFCondition(const int index){

    if(Phydro[index]->RhoPred < USE_FORCIBLE_STARFORMATION_DENSITY_CRITERION_FACTOR*Pall.RHOcZ){
        return true;
    }
    return false;

}

#define StarFormationTimeScale  (1.e+9*YEAR_CGS) // 1 Gyr in sec
static bool ConvertGasParticleIntoStarParticle(const int index){

#ifdef USE_STARFORMATION_BIGSTEP //{
    double SFBigStepFactor = GetSFBigStepFactor(index);
#if (UseConstantSFrate)
    double InvTsf = Pall.UnitTime/StarFormationTimeScale;
    double p = 1.0-exp(-SFBigStepFactor*Phydro[index]->dt_hydro*InvTsf);
#else
    double InvTdyn = sqrt(4.0*PI*Pall.GravConst*Phydro[index]->Rho);
    double p = 1.0-exp(-SFeff*SFBigStepFactor*Phydro[index]->dt_hydro*InvTdyn);
#endif
#else  // USE_STARFORMATION_BIGSTEP //}//{
#if (UseConstantSFrate)
    double InvTsf = Pall.UnitTime/StarFormationTimeScale;
    double p = 1.0-exp(-Phydro[index]->dt_hydro*InvTsf);
#else
    double InvTdyn = sqrt(4.0*PI*Pall.GravConst*Phydro[index]->Rho);
    double p = 1.0-exp(-SFeff*Phydro[index]->dt_hydro*InvTdyn);
#endif
#endif // USE_STARFORMATION_BIGSTEP //}
    
    if(gsl_rng_uniform(RandomGenerator)<p){ // born a star particle

        StructPstarptr Ps = ReturnEmptyStarStructurePointer();

        Ps->Use = ON;
        Ps->IMFTYPE = IMFTYPE_SP;
        Ps->Mass = PhydroMass(index);
        Ps->InitialMass = PhydroMass(index);
        Ps->FormationTime = Pall.TCurrent;
        Ps->Z = Phydro[index]->Z;
        Ps->ZII = Phydro[index]->ZII;
        Ps->ZIa = Phydro[index]->ZIa;
        Ps->TypeII = false;
        Ps->TypeIa = false;

        Ps->TempForm = Pall.ConvertUtoT*Phydro[index]->U;
        Ps->RhoForm = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho;

#if (UseSFModelSpawn)
        Ps->NthChildren = Phydro[index]->SpawnTimes;
#else
        Ps->NthChildren = 0;
#endif
        Ps->ParentGlobalID = PhydroBody(index)->GlobalID;

        Ps->Body = PhydroBody(index);
        Ps->Body->Baryon = (void*)Ps;

        // Energy loss
        Pall.CoolingEnergyLoss += 
            PhydroBody(index)->Mass*
                (Phydro[index]->U+0.5*Phydro[index]->dt_hydro*Phydro[index]->Du);

        Ps->Body->Type = TypeStar;
        
        // Add hydro kick.
        double dt_half_hydro = 0.5*Phydro[index]->dt_hydro;
        PhydroBody(index)->Vel[0] = PhydroBody(index)->Velh[0]+dt_half_hydro*Phydro[index]->HydroAcc[0];
        PhydroBody(index)->Vel[1] = PhydroBody(index)->Velh[1]+dt_half_hydro*Phydro[index]->HydroAcc[1];
        PhydroBody(index)->Vel[2] = PhydroBody(index)->Velh[2]+dt_half_hydro*Phydro[index]->HydroAcc[2];
        
        HydroRoot.Leaves[Phydro[index]->Leaf] *= -1;

        Phydro[index]->Use = OFF;

#ifdef USE_CELIB
        Ps->IMFTYPE = StellarFeedbackGetIMFType();
        Ps->SNIICount = 0;
        Ps->SNIaCount = 0;
        Ps->EventTimeSNII = Ps->FormationTime 
            +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                .R = gsl_rng_uniform(RandomGenerator),
                .InitialMass_in_Msun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS,
                .Metallicity = Ps->Z,
                .Count = 0,
#ifdef USE_CELIB_SNII_INDIVIDUAL //{
                .Mode = CELibSNIIRateModelID_Individual,
                .Nassociation = CHEMICALEVOLUTION_SNII_EVENTNUMBER,
#else // USE_CELIB_SNII_INDIVIDUAL //}//{
                .Mode = CELibSNIIRateModelID_Once,
#endif // USE_CELIB_SNII_INDIVIDUAL //}
                },CELibFeedbackType_SNII)
                *YEAR_CGS/Pall.UnitTime;
#ifdef USE_CELIB_AGB //{
        Ps->AGBCount = 0;
        Ps->EventTimeAGB = Ps->FormationTime 
            +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                .R = gsl_rng_uniform(RandomGenerator),
                .InitialMass_in_Msun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS,
                .Metallicity = Ps->Z,
                .Count = 0,
                },CELibFeedbackType_AGB)
                *YEAR_CGS/Pall.UnitTime;

#endif //USE_CELIB_AGB //}

#ifdef USE_CELIB_NSM //{
        Ps->NSMCount = 0;
        Ps->EventTimeNSM = Ps->FormationTime
            +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                .R = gsl_rng_uniform(RandomGenerator),
                .InitialMass_in_Msun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS,
                .Metallicity = Ps->Z,
                .Count = 0,
                },CELibFeedbackType_NSM)
                *YEAR_CGS/Pall.UnitTime;

#endif //USE_CELIB_NSM //}

        for(int k=0;k<CELibYield_Number;k++){
            Ps->Elements[k] = Phydro[index]->Elements[k];
        }
#endif // USE_CELIB

#ifdef USE_FUVFEEDBACK //{
        Ps->LFUV = (Pall.UnitMass/MSUN_CGS)*ASRFLXGetFUV(0.0,Ps->Z); 
#endif // USE_FUVFEEDBACK //}

#if defined(PRESERVE_SNII_EVENTRATE) //{
        double MassInMsun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS;
        double prob = SNIINumber*MassInMsun;
        if(prob >= 1.0){
            Ps->TypeIIProb = true;
        } else {
            if(prob > gsl_rng_uniform(RandomGenerator)){
                Ps->TypeIIProb = true;
            }else{
                Ps->TypeIIProb = false;
            }
        }
#endif // PRESERVE_SNII_EVENTRATE //}

        TempNhydro --;
        TempNstars ++;
        return true;
    } else {
        return false;
    }
}


static bool ConvertGasParticleIntoDeadStarParticle(const int index){

    StructPstarptr Ps = ReturnEmptyStarStructurePointer();

    Ps->Use = ON;
    Ps->IMFTYPE = IMFTYPE_SP;
    Ps->Mass = PhydroMass(index);
    Ps->InitialMass = PhydroMass(index);
    Ps->FormationTime = Pall.TCurrent;
    Ps->Z = Phydro[index]->Z;
    Ps->ZII = Phydro[index]->ZII;
    Ps->ZIa = Phydro[index]->ZIa;
    Ps->TypeII = true;
    Ps->TypeIa = false;

    Ps->TempForm = Pall.ConvertUtoT*Phydro[index]->U;
    Ps->RhoForm = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho;

#if (UseSFModelSpawn)
    Ps->NthChildren = Phydro[index]->SpawnTimes;
#else
    Ps->NthChildren = 0;
#endif
    Ps->ParentGlobalID = PhydroBody(index)->GlobalID;

    Ps->Body = PhydroBody(index);
    Ps->Body->Baryon = (void*)Ps;

    // Energy loss
    Pall.CoolingEnergyLoss += 
        PhydroBody(index)->Mass*
            (Phydro[index]->U+0.5*Phydro[index]->dt_hydro*Phydro[index]->Du);

    Ps->Body->Type = TypeStar;
    
    // Add hydro kick.
    double dt_half_hydro = 0.5*Phydro[index]->dt_hydro;
    PhydroBody(index)->Vel[0] = PhydroBody(index)->Velh[0]+dt_half_hydro*Phydro[index]->HydroAcc[0];
    PhydroBody(index)->Vel[1] = PhydroBody(index)->Velh[1]+dt_half_hydro*Phydro[index]->HydroAcc[1];
    PhydroBody(index)->Vel[2] = PhydroBody(index)->Velh[2]+dt_half_hydro*Phydro[index]->HydroAcc[2];
    
    HydroRoot.Leaves[Phydro[index]->Leaf] *= -1;

    Phydro[index]->Use = OFF;

#ifdef USE_CELIB
    Ps->IMFTYPE = StellarFeedbackGetIMFType();
    Ps->SNIICount = NONE;
    Ps->EventTimeSNII = 10*Pall.TEnd;
    Ps->SNIaCount = 0;
    Ps->EventTimeSNIa = 10*Pall.TEnd;
#ifdef USE_CELIB_AGB //{
    Ps->AGBCount = 0;
    Ps->EventTimeAGB = 10*Pall.TEnd;
#endif //USE_CELIB_AGB //}
#ifdef USE_CELIB_NSM //{
    Ps->NSMCount = 0;
    Ps->EventTimeNSM = 10*Pall.TEnd;
#endif //USE_CELIB_NSM //}

    for(int k=0;k<CELibYield_Number;k++){
        Ps->Elements[k] = Phydro[index]->Elements[k];
    }
#endif // USE_CELIB

#ifdef USE_FUVFEEDBACK //{
        Ps->LFUV = 1.e-30*(Pall.UnitMass/MSUN_CGS)*ASRFLXGetFUV(0.0,Ps->Z); 
#endif // USE_FUVFEEDBACK //}

#if defined(PRESERVE_SNII_EVENTRATE)
    Ps->TypeIIProb = false;
#endif

    TempNhydro --;
    TempNstars ++;

    return true;
}

#if (UseSFModelSpawn) //{
static bool SpawnStarParticle(const int index){

#ifdef USE_STARFORMATION_BIGSTEP //{
    double SFBigStepFactor = GetSFBigStepFactor(index);
#if (UseConstantSFrate) //{
    double InvTsf = Pall.UnitTime/StarFormationTimeScale;
    double p = 1.0-exp(-SFBigStepFactor*Phydro[index]->dt_hydro*InvTsf);
#else // UseConstantSFrate //}//{
    double InvTdyn = sqrt(4.0*PI*Pall.GravConst*Phydro[index]->Rho);
    double p = 1.0-exp(-SFeff*SFBigStepFactor*Phydro[index]->dt_hydro*InvTdyn);
#endif // UseConstantSFrate //}
#else  // USE_STARFORMATION_BIGSTEP //}//{
#if (UseConstantSFrate) //{
    double InvTsf = Pall.UnitTime/StarFormationTimeScale;
    double p = 1.0-exp(-Phydro[index]->dt_hydro*InvTsf);
#else // UseConstantSFrate //}//{
    double InvTdyn = sqrt(4.0*PI*Pall.GravConst*Phydro[index]->Rho);
    double p = 1.0-exp(-SFeff*Phydro[index]->dt_hydro*InvTdyn);
#endif // UseConstantSFrate //}
#endif // USE_STARFORMATION_BIGSTEP //}

    p *= PhydroBody(index)->Mass/Phydro[index]->SpawnMass;

    double R = gsl_rng_uniform(RandomGenerator);
    if(R<p){ // born a star particle
        StructPstarptr Ps = ReturnEmptyStarStructurePointer();
        StructPbodyptr Pb = ReturnEmptyBodyStructurePointer();

        StructPbodyptr PbNext = Pb->Next;
        *Pb = *(PhydroBody(index));
        Pb->Next = PbNext;
        Pb->Baryon = (void *)Ps;
        Pb->Type = TypeStar;

        Ps->Use = ON;
        Ps->IMFTYPE = IMFTYPE_SP;
        Ps->Mass = Phydro[index]->SpawnMass;
        Ps->InitialMass = Phydro[index]->SpawnMass;
        Ps->FormationTime = Pall.TCurrent;
        Ps->Z = Phydro[index]->Z;
        Ps->ZII = Phydro[index]->ZII;
        Ps->ZIa = Phydro[index]->ZIa;
        Ps->TypeII = false;
        Ps->TypeIa = false;
        Ps->NthChildren = Phydro[index]->SpawnTimes;
        Ps->ParentGlobalID = PhydroBody(index)->GlobalID;

        Ps->TempForm = Pall.ConvertUtoT*Phydro[index]->U;
        Ps->RhoForm = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho;

        Ps->Body = Pb;
        Pb->Mass = Phydro[index]->SpawnMass;

        // Energy loss
        Pall.CoolingEnergyLoss += 
            Ps->Mass*(Phydro[index]->U+0.5*Phydro[index]->dt_hydro*Phydro[index]->Du);

#ifdef USE_CELIB //{
        double OriginalMass = Phydro[index]->Mass;
#endif // USE_CELIB //}
        Phydro[index]->Mass -= Phydro[index]->SpawnMass;
        PhydroBody(index)->Mass = Phydro[index]->Mass;
        Phydro[index]->SpawnTimes ++;

        // Add hydro kick.
        double dt_half_hydro = 0.5*Phydro[index]->dt_hydro;
        PhydroBody(index)->Vel[0] = PhydroBody(index)->Velh[0]+dt_half_hydro*Phydro[index]->HydroAcc[0];
        PhydroBody(index)->Vel[1] = PhydroBody(index)->Velh[1]+dt_half_hydro*Phydro[index]->HydroAcc[1];
        PhydroBody(index)->Vel[2] = PhydroBody(index)->Velh[2]+dt_half_hydro*Phydro[index]->HydroAcc[2];

#ifdef USE_CELIB
        Ps->IMFTYPE = StellarFeedbackGetIMFType();
        Ps->SNIICount = 0;
        Ps->SNIaCount = 0;
        Ps->EventTimeSNII = Ps->FormationTime 
                        +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                            .R = gsl_rng_uniform(RandomGenerator),
                            .InitialMass_in_Msun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS,
                            .Metallicity = Ps->Z,
                            .Count = 0,
#ifdef USE_CELIB_SNII_INDIVIDUAL //{
                            .Mode = CELibSNIIRateModelID_Individual,
                            .Nassociation = CHEMICALEVOLUTION_SNII_EVENTNUMBER,
#else // USE_CELIB_SNII_INDIVIDUAL //}//{
                            .Mode = CELibSNIIRateModelID_Once,
#endif // USE_CELIB_SNII_INDIVIDUAL //}
                            },CELibFeedbackType_SNII)
                            *YEAR_CGS/Pall.UnitTime;
#ifdef USE_CELIB_AGB //{
        Ps->AGBCount = 0;
        Ps->EventTimeAGB = Ps->FormationTime 
                        +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                            .R = gsl_rng_uniform(RandomGenerator),
                            .InitialMass_in_Msun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS,
                            .Metallicity = Ps->Z,
                            .Count = 0,
                            },CELibFeedbackType_AGB)
                            *YEAR_CGS/Pall.UnitTime;
#endif //USE_CELIB_AGB //}

#ifdef USE_CELIB_NSM //{
        Ps->NSMCount = 0;
        Ps->EventTimeNSM = Ps->FormationTime
            +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                .R = gsl_rng_uniform(RandomGenerator),
                .InitialMass_in_Msun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS,
                .Metallicity = Ps->Z,
                .Count = 0,
                },CELibFeedbackType_NSM)
                *YEAR_CGS/Pall.UnitTime;

#endif //USE_CELIB_NSM //}

        double fs = Pb->Mass/OriginalMass;
        for(int k=0;k<CELibYield_Number;k++){
            Ps->Elements[k] = fs*Phydro[index]->Elements[k];
            Phydro[index]->Elements[k] -= Ps->Elements[k];
        }
#endif // USE_CELIB

#ifdef USE_FUVFEEDBACK //{
        Ps->LFUV = (Pall.UnitMass/MSUN_CGS)*ASRFLXGetFUV(0.0,Ps->Z); 
#endif // USE_FUVFEEDBACK //}

#if defined(PRESERVE_SNII_EVENTRATE) //{
        double MassInMsun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS;
        double prob = SNIINumber*MassInMsun;
        if(prob >= 1.0){
            Ps->TypeIIProb = true;
        } else {
            if(prob > gsl_rng_uniform(RandomGenerator)){
                Ps->TypeIIProb = true;
            }else{
                Ps->TypeIIProb = false;
            }
        }
#endif // PRESERVE_SNII_EVENTRATE //}

        TempNstars ++;
        TempNtotal ++;
        return true;
    } else {
        return false;
    }
}
#endif //UseSFModelSpawn //}

static void ReConnectBaryonPointers(void){

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

double CheckTotalBaryonMass(void){ 

    double Mass = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++)
        Mass += PhydroMass(i);
    for(int i=0;i<Pall.Nstars;i++)
        Mass += PstarMass(i);
    double GlobalMass;
    MPI_Allreduce(&Mass,&GlobalMass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return GlobalMass;
}

#if (UseSFModelSpawn)
/*
 *This function is used for ArtificialStarFormation().
 *
 */
int ArtificialTempNstars,ArtificialTempNtotal,ArtificialTempNhydro;
static void ArtificialSpawnStarParticle(const int index, const double FormationTime){

    StructPstarptr Ps = ReturnEmptyStarStructurePointer();
    StructPbodyptr Pb = ReturnEmptyBodyStructurePointer();

    StructPbodyptr PbNext = Pb->Next;
    *Pb = *(PhydroBody(index));
    Pb->Next = PbNext;
    Pb->Baryon = (void *)Ps;
    Pb->Type = TypeStar;

    Ps->Use = ON;
    Ps->IMFTYPE = IMFTYPE_SP;
    Ps->Mass = Phydro[index]->SpawnMass;
    Ps->InitialMass = Phydro[index]->SpawnMass;
    Ps->FormationTime = FormationTime;
    Ps->Z = Phydro[index]->Z;
    Ps->ZII = Phydro[index]->ZII;
    Ps->ZIa = Phydro[index]->ZIa;
    Ps->TypeII = false;
    Ps->TypeIa = false;
    Ps->NthChildren = Phydro[index]->SpawnTimes;
    Ps->ParentGlobalID = PhydroBody(index)->GlobalID;

    Ps->TempForm = Pall.ConvertUtoT*Phydro[index]->U;
    Ps->RhoForm = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho;

    Ps->Body = Pb;
    Pb->Mass = Phydro[index]->SpawnMass;

    Pall.CoolingEnergyLoss +=
        Ps->Mass*(Phydro[index]->U+0.5*Phydro[index]->dt_hydro*Phydro[index]->Du);

    Phydro[index]->Mass -= Phydro[index]->SpawnMass;
    PhydroBody(index)->Mass = Phydro[index]->Mass;
    Phydro[index]->SpawnTimes ++;

    double dt_half_hydro = 0.5*Phydro[index]->dt_hydro;
    PhydroBody(index)->Vel[0] = PhydroBody(index)->Velh[0]+dt_half_hydro*Phydro[index]->HydroAcc[0];
    PhydroBody(index)->Vel[1] = PhydroBody(index)->Velh[1]+dt_half_hydro*Phydro[index]->HydroAcc[1];
    PhydroBody(index)->Vel[2] = PhydroBody(index)->Velh[2]+dt_half_hydro*Phydro[index]->HydroAcc[2];

#if defined(PRESERVE_SNII_EVENTRATE)
    double MassInMsun = Ps->InitialMass*Pall.UnitMass/MSUN_CGS;
    double prob = SNIINumber*MassInMsun;
    if(prob >= 1.0){
        Ps->TypeIIProb = true;
    } else {
        double p = gsl_rng_uniform(RandomGenerator);
        //if(prob > gsl_rng_uniform(RandomGenerator)){
        // fprintf(stderr,"%g %g\n",p,prob);
        if(prob > p){
            Ps->TypeIIProb = true;
        }else{
            Ps->TypeIIProb = false;
        }
    }
#endif

    ArtificialTempNstars ++;
    ArtificialTempNtotal ++;

    return;
}
#endif


/*
 * This function convert a part of hydro particles into star particles
 * following input SFR and duration time.
 * The unit of `SFR' is Msun/yr, and that of `Duration' is yr.
 */
void ArtificialStarFormation(const double SFR, const double Duration){

    if(Pall.RunStatus != NewSimulation) return;

    double DurationTime = Duration*YEAR_CGS/Pall.UnitTime;
    double NewStarsMass = (SFR*MSUN_CGS/YEAR_CGS)*Duration*YEAR_CGS/Pall.UnitMass;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"SFR,DurationTime %g %g\n",SFR,Duration);
        fprintf(stderr,"DurationTime in simulation unit %g\n",DurationTime);
        fprintf(stderr,"Stellar Mass %g\n",NewStarsMass);
    }

    if(Pall.Nhydro == 0)
        abort();
    MPI_Barrier(MPI_COMM_WORLD);

    double stellar_mass = Phydro[0]->SpawnMass;
    double global_stellar_mass;
    MPI_Allreduce(&stellar_mass,&global_stellar_mass,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    if(global_stellar_mass == 0)
        abort();
    MPI_Barrier(MPI_COMM_WORLD);
    stellar_mass = global_stellar_mass;

    //how many stars should we form?
    int Nstars = NewStarsMass/stellar_mass;
    int stride = Pall.Nhydro_t/Nstars;
    if(stride < 1){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"Star formation is difficult.\n");
            fprintf(stderr,"The stride is %d.\n",stride);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The stride is %d.\n",stride);
    }

    ArtificialTempNstars = ArtificialTempNhydro = ArtificialTempNtotal = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroBody(i)->GlobalID%stride == 0){
        //if(PhydroBody(i)->GlobalID%stride == (int)(gsl_rng_uniform(RandomGenerator)+0.5)){
        //if((PhydroBody(i)->GlobalID+(int)(4*gsl_rng_uniform(RandomGenerator)+0.5))%stride == 0){
            ResetRandomSeedForRandomGenerator(PhydroBody(i)->GlobalID);
            double FormationTime = DurationTime*gsl_rng_uniform(RandomGenerator);
            ArtificialSpawnStarParticle(i,-FormationTime);
            //dprintlmpi(i);
        }
    }
    dprintlmpi(ArtificialTempNtotal);
    dprintlmpi(ArtificialTempNhydro);
    dprintlmpi(ArtificialTempNstars);

    Pall.Nhydro += ArtificialTempNhydro;
    Pall.NActivesHydro += ArtificialTempNhydro;
    Pall.Nstars += ArtificialTempNstars;
    Pall.NActivesStars += ArtificialTempNstars;
    Pall.Ntotal += ArtificialTempNtotal;
    Pall.NActives += ArtificialTempNtotal;

    ReConnectPointers();
    UpdateTotalNumber();
    ResetRandomSeedForRandomGenerator(1977+MPIGetMyID());

    return ;
}   
