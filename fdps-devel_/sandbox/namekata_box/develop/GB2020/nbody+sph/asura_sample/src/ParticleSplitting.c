#include "config.h"
#include "StructureOperation.h"

#ifdef USE_PARTICLE_SPLITTING //{

static bool SplitCondition(const int Index){
    
    if(Phydro[Index]->SplitGeneration >= MAX_PARTICLE_SPLIT_GENERATION)
        return false;

    double cs = sqrt(Pall.GGm1*Phydro[Index]->UPred);
    double LJeans = cs/sqrt(Pall.GravConst*Phydro[Index]->RhoPred);

    const double factor = M_PI/6.0;
    double MJeans = factor*Phydro[Index]->RhoPred*CUBE(LJeans);

    if(MJeans < PARTICLE_SPLITTING_MASSFACTOR*Phydro[Index]->Mass){
        return true;
    } else {
        return false;
    }
}

#ifdef USE_AGRESSIVE_PARTICLE_SPLITTING_STRATEGY //{
enum SplitColor{
    Split,
    SplitInAdvance,
    NoSplit,
    SplitTypes,
};

static int SplitConditionForAgressiveParticleSplittingStrategy(const int Index){

    if(Phydro[Index]->SplitGeneration >= MAX_PARTICLE_SPLIT_GENERATION)
        return NoSplit;
    
    double cs = sqrt(Pall.GGm1*Phydro[Index]->UPred);
    double LJeans = cs/sqrt(Pall.GravConst*Phydro[Index]->RhoPred);

    const double factor = M_PI/6.0;
    double MJeans = factor*Phydro[Index]->RhoPred*CUBE(LJeans);

    if(MJeans < PARTICLE_SPLITTING_MASSFACTOR*Phydro[Index]->Mass){
        return Split;
    } else if(AGRESSIVE_PARTICLE_SPLITTING_THRESHOLD*MJeans < PARTICLE_SPLITTING_MASSFACTOR*Phydro[Index]->Mass){
        return SplitInAdvance;
    } else {
        return NoSplit;
    }
}
#endif // USE_AGRESSIVE_PARTICLE_SPLITTING_STRATEGY //}

static int CheckSplitCondition(void){

    int counter = 0;
#ifdef USE_AGRESSIVE_PARTICLE_SPLITTING_STRATEGY //{
    int SplitFlags[Pall.Nhydro];
    int Colors[SplitTypes] = {0};
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            SplitFlags[i] = SplitConditionForAgressiveParticleSplittingStrategy(i);
            Colors[SplitFlags[i]] ++;
            Phydro[i]->SplitFlag = false;
        }
    }
    if(Colors[Split] > 0){
        for(int i=0;i<Pall.Nhydro;i++){
            if((Phydro[i]->Active)&&(SplitFlags[i]!=NoSplit)){
                Phydro[i]->SplitFlag = true;
            }
        }
    }
    counter = Colors[0] + Colors[1];
#else // USE_AGRESSIVE_PARTICLE_SPLITTING_STRATEGY //}//{
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->SplitFlag = SplitCondition(i);
            if(Phydro[i]->SplitFlag)
                counter ++;
        }
    }
#endif // USE_AGRESSIVE_PARTICLE_SPLITTING_STRATEGY //}

    return counter;
}


#if (PARTICLE_SPLITTING_TYPE == 0) //{
double x[] = {0.0, 1.0, 0.5, -0.5, -1.0, -0.5, 0.5, -0.5, 0.5, 1.76763e-16, -0.5, 0.5, 1.76763e-16};
double y[] = {0.0, 0.0, 0.866025, 0.866025, 1.22465e-16, -0.866025, -0.866025, -0.288675, -0.288675, 0.57735, -0.288675, -0.288675, 0.57735};
double z[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.816497, 0.816497, 0.816497, -0.816497, -0.816497, -0.816497};


static void DaughterParticlePositions(double Pos[][3]){

    double alpha = 2*M_PI*(gsl_rng_uniform(RandomGenerator)-0.5);
    double beta  = M_PI*gsl_rng_uniform(RandomGenerator);
    double gamma = 2*M_PI*(gsl_rng_uniform(RandomGenerator)-0.5);

// Z1X2Z3 transformation.

    double c1 = cos(alpha);
    double s1 = sin(alpha);
    double c2 = cos(beta);
    double s2 = sin(beta);
    double c3 = cos(gamma);
    double s3 = sin(gamma);


    for(int i=0;i<13;i++){
        Pos[i][0] = (c1*c3-c2*s1*s3)*x[i] + (-c1*s3-c2*c3*s1)*y[i] + s1*s2*z[i];
        Pos[i][1] = (c3*s1+c1*c2*s3)*x[i] + (c1*c2*c3-s1*s3) *y[i] - c1*s2*z[i];
        Pos[i][2] = s2*s3           *x[i] + c3*s2            *y[i] + c2   *z[i];

    }

    return; 
}

static void GenChildParticles(const int LocalSplittingTarget){

    const double Inv13 = 1.0/13.0;
    const double SqrtInv13 = sqrt(Inv13);
    const double CubeInv13 = pow(13.0,-1.0/3.0);

    int counter = 0;
    int PhCounter = 0;
    int PbCounter = 0;
    if(LocalSplittingTarget>0){
        for(int i=0;i<Pall.Nhydro;i++){
            if((Phydro[i]->Active)&&(Phydro[i]->SplitFlag)){

                double dr[13][3];
                DaughterParticlePositions(dr);

                double hchild = CubeInv13*Phydro[i]->Kernel;
                double Lscale = 1.5*hchild;

                // Get 12 daughters
                for(int k=0;k<12;k++){
                    //StructPhydroptr Ph = ReturnEmptyHydroStructurePointer();
                    StructPhydroptr Ph = ReturnEmptyHydroStructurePointerWithCounter(&PhCounter);
                    //StructPbodyptr  Pb = ReturnEmptyBodyStructurePointer();
                    StructPbodyptr  Pb = ReturnEmptyBodyStructurePointerWithCounter(&PbCounter);

                    StructPbodyptr  PbNext = Pb->Next;
                    StructPhydroptr PhNext = Ph->Next;
                    *Ph = *(Phydro[i]);
                    *Pb = *(PhydroBody(i));
                    Pb->Next = PbNext;
                    Ph->Next = PhNext;
                    Pb->Baryon = (void *)Ph;
                    Ph->Body = (void *)Pb;
                    
                    Pb->Mass *= Inv13;
                    Ph->Mass *= Inv13;
#ifdef USE_MULTIPHASE_MODEL //{
                    Ph->Mhot *= Inv13;
#endif //USE_MULTIPHASE_MODEL //}
                    for(int l=0;l<CELibYield_Number;l++){
                        Ph->Elements[l] *= Inv13;
                    }

                    // Position updates
                    Pb->PosP[0] += dr[k+1][0]*Lscale;
                    Pb->PosP[1] += dr[k+1][1]*Lscale;
                    Pb->PosP[2] += dr[k+1][2]*Lscale;
                    Ph->PosP[0] = Pb->Pos[0] = Pb->PosP[0];
                    Ph->PosP[1] = Pb->Pos[1] = Pb->PosP[1];
                    Ph->PosP[2] = Pb->Pos[2] = Pb->PosP[2];

                    // Pb->Eps *= SqrtInv13;
                    Ph->KernelPred = Ph->Kernel = hchild;

                    Ph->SplitGeneration ++;
                    Ph->SplitNthChildren = k+1;
                    Ph->SplitFlag = false; 
                    if(Ph->SplitGeneration == 1){
                        Ph->AncestorGlobalID = PhydroBody(i)->GlobalID;
                    } else {
                        Ph->AncestorGlobalID = Phydro[i]->AncestorGlobalID;
                    }
                    Ph->Body->GlobalID = PowderSnowIDGenerator();
                }

                Phydro[i]->Mass *= Inv13;
                PhydroBody(i)->Mass *= Inv13;
#ifdef USE_MULTIPHASE_MODEL //{
                Phydro[i]->Mhot *= Inv13;
#endif //USE_MULTIPHASE_MODEL //}

                for(int l=0;l<CELibYield_Number;l++){
                    Phydro[i]->Elements[l] *= Inv13;
                }
                // PhydroBody(i)->Eps *= SqrtInv13;
                Phydro[i]->KernelPred = Phydro[i]->Kernel = hchild;

                Phydro[i]->SplitGeneration ++;
                Phydro[i]->SplitNthChildren = 0;
                Phydro[i]->SplitFlag = false; 

                // Keep global ID
                if(Phydro[i]->SplitGeneration == 1){
                    Phydro[i]->AncestorGlobalID = Phydro[i]->GlobalID;
                } else {
                    Phydro[i]->AncestorGlobalID = Phydro[i]->AncestorGlobalID;
                }

                counter += 12;
            }
        }
    }

    Pall.Nhydro        += counter;
    Pall.NActivesHydro += counter;
    Pall.Ntotal        += counter;
    Pall.NActives      += counter;

    UpdateTotalNumber();
    ReConnectPointers();

    return ;
}
#elif (PARTICLE_SPLITTING_TYPE == 1) //}//{
double x[] = {-1, 1,-1, 1,-1, 1,-1, 1};
double y[] = {-1,-1, 1, 1,-1,-1, 1, 1};
double z[] = {-1,-1,-1,-1, 1, 1, 1, 1};

static void DaughterParticlePositions(double Pos[][3]){

    double alpha = 2*M_PI*(gsl_rng_uniform(RandomGenerator)-0.5);
    double beta  = M_PI*gsl_rng_uniform(RandomGenerator);
    double gamma = 2*M_PI*(gsl_rng_uniform(RandomGenerator)-0.5);

// Z1X2Z3 transformation.

    double c1 = cos(alpha);
    double s1 = sin(alpha);
    double c2 = cos(beta);
    double s2 = sin(beta);
    double c3 = cos(gamma);
    double s3 = sin(gamma);


    for(int i=0;i<8;i++){
        Pos[i][0] = (c1*c3-c2*s1*s3)*x[i] + (-c1*s3-c2*c3*s1)*y[i] + s1*s2*z[i];
        Pos[i][1] = (c3*s1+c1*c2*s3)*x[i] + (c1*c2*c3-s1*s3) *y[i] - c1*s2*z[i];
        Pos[i][2] = s2*s3           *x[i] + c3*s2            *y[i] + c2   *z[i];

    }

    return; 
}

static void GenChildParticles(const int LocalSplittingTarget){

    const double Inv8 = 1.0/8.0;
    const double SqrtInv8 = sqrt(Inv8);
    const double CubeInv8 = pow(8.0,-1.0/3.0);
    // const double offset_scale = 2.0/8.0;
    const double offset_scale = 1.0/cbrt(Pall.Ns);

    int counter = 0;
    int PhCounter = 0;
    int PbCounter = 0;
    if(LocalSplittingTarget>0){
        for(int i=0;i<Pall.Nhydro;i++){
            if((Phydro[i]->Active)&&(Phydro[i]->SplitFlag)){

                double dr[8][3];
                DaughterParticlePositions(dr);

                double hchild = Phydro[i]->Kernel/2.0;
                double Lscale = Phydro[i]->Kernel*offset_scale;

                // double hchild = Phydro[i]->Kernel/2.0;
                // double Lscale = 1.5*hchild;

                // Get 7 daughters
                for(int k=0;k<7;k++){
                    StructPhydroptr Ph = ReturnEmptyHydroStructurePointerWithCounter(&PhCounter);
                    StructPbodyptr  Pb = ReturnEmptyBodyStructurePointerWithCounter(&PbCounter);

                    StructPbodyptr  PbNext = Pb->Next;
                    StructPhydroptr PhNext = Ph->Next;
                    *Ph = *(Phydro[i]);
                    *Pb = *(PhydroBody(i));
                    Pb->Next = PbNext;
                    Ph->Next = PhNext;
                    Pb->Baryon = (void *)Ph;
                    Ph->Body = (void *)Pb;
                    
                    Pb->Mass *= Inv8;
                    Ph->Mass *= Inv8;
#ifdef USE_MULTIPHASE_MODEL //{
                    Ph->Mhot *= Inv8;
#endif //USE_MULTIPHASE_MODEL //}
                    for(int l=0;l<CELibYield_Number;l++){
                        Ph->Elements[l] *= Inv8;
                    }

                    // Position updates
                    Pb->PosP[0] += dr[k+1][0]*Lscale;
                    Pb->PosP[1] += dr[k+1][1]*Lscale;
                    Pb->PosP[2] += dr[k+1][2]*Lscale;
                    Ph->PosP[0] = Pb->Pos[0] = Pb->PosP[0];
                    Ph->PosP[1] = Pb->Pos[1] = Pb->PosP[1];
                    Ph->PosP[2] = Pb->Pos[2] = Pb->PosP[2];

                    Ph->KernelPred = Ph->Kernel = hchild;

                    Ph->SplitGeneration ++;
                    Ph->SplitNthChildren = k+1;
                    Ph->SplitFlag = false; 
                    if(Ph->SplitGeneration == 1){
                        Ph->AncestorGlobalID = PhydroBody(i)->GlobalID;
                    } else {
                        Ph->AncestorGlobalID = Phydro[i]->AncestorGlobalID;
                    }
                    Ph->Body->GlobalID = PowderSnowIDGenerator();
                }

                Phydro[i]->Mass *= Inv8;
                PhydroBody(i)->Mass *= Inv8;
#ifdef USE_MULTIPHASE_MODEL //{
                Phydro[i]->Mhot *= Inv8;
#endif //USE_MULTIPHASE_MODEL //}
                for(int l=0;l<CELibYield_Number;l++){
                    Phydro[i]->Elements[l] *= Inv8;
                }

                // Position updates
                PhydroBody(i)->PosP[0] += dr[0][0]*Lscale;
                PhydroBody(i)->PosP[1] += dr[0][1]*Lscale;
                PhydroBody(i)->PosP[2] += dr[0][2]*Lscale;
                Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0] = PhydroBody(i)->PosP[0];
                Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1] = PhydroBody(i)->PosP[1];
                Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2] = PhydroBody(i)->PosP[2];


                Phydro[i]->KernelPred = Phydro[i]->Kernel = hchild;

                Phydro[i]->SplitGeneration ++;
                Phydro[i]->SplitNthChildren = 0;
                Phydro[i]->SplitFlag = false; 

                // Keep global ID
                if(Phydro[i]->SplitGeneration == 1){
                    Phydro[i]->AncestorGlobalID = PhydroBody(i)->GlobalID;
                } else {
                    Phydro[i]->AncestorGlobalID = Phydro[i]->AncestorGlobalID;
                }

                counter += 7;
            }
        }
    }

    Pall.Nhydro        += counter;
    Pall.NActivesHydro += counter;
    Pall.Ntotal        += counter;
    Pall.NActives      += counter;

    UpdateTotalNumber();
    ReConnectPointers();

    return ;
}

#elif (PARTICLE_SPLITTING_TYPE == 2) //}//{
double x[] = {0,-1, 1,-1, 1,-1, 1,-1, 1};
double y[] = {0,-1,-1, 1, 1,-1,-1, 1, 1};
double z[] = {0,-1,-1,-1,-1, 1, 1, 1, 1};

static void DaughterParticlePositions(double Pos[][3]){

    double alpha = 2*M_PI*(gsl_rng_uniform(RandomGenerator)-0.5);
    double beta  = M_PI*gsl_rng_uniform(RandomGenerator);
    double gamma = 2*M_PI*(gsl_rng_uniform(RandomGenerator)-0.5);

// Z1X2Z3 transformation.

    double c1 = cos(alpha);
    double s1 = sin(alpha);
    double c2 = cos(beta);
    double s2 = sin(beta);
    double c3 = cos(gamma);
    double s3 = sin(gamma);


    for(int i=0;i<9;i++){
        Pos[i][0] = (c1*c3-c2*s1*s3)*x[i] + (-c1*s3-c2*c3*s1)*y[i] + s1*s2*z[i];
        Pos[i][1] = (c3*s1+c1*c2*s3)*x[i] + (c1*c2*c3-s1*s3) *y[i] - c1*s2*z[i];
        Pos[i][2] = s2*s3           *x[i] + c3*s2            *y[i] + c2   *z[i];
    }

    return; 
}

static void GenChildParticles(const int LocalSplittingTarget){

    const double Inv9 = 1.0/9.0;
    const double SqrtInv9 = sqrt(Inv9);
    const double CubeInv9 = pow(9.0,-1.0/3.0);
    const double offset_scale = 1.0/9.0;

    int counter = 0;
    int PhCounter = 0;
    int PbCounter = 0;
    if(LocalSplittingTarget>0){
        for(int i=0;i<Pall.Nhydro;i++){
            if((Phydro[i]->Active)&&(Phydro[i]->SplitFlag)){

                double dr[9][3];
                DaughterParticlePositions(dr);

                double hchild = Phydro[i]->Kernel/2.0;
                double Lscale = Phydro[i]->Kernel*offset_scale;

                // Get 8 daughters
                for(int k=0;k<8;k++){
                    StructPhydroptr Ph = ReturnEmptyHydroStructurePointerWithCounter(&PhCounter);
                    StructPbodyptr  Pb = ReturnEmptyBodyStructurePointerWithCounter(&PbCounter);

                    StructPbodyptr  PbNext = Pb->Next;
                    StructPhydroptr PhNext = Ph->Next;
                    *Ph = *(Phydro[i]);
                    *Pb = *(PhydroBody(i));
                    Pb->Next = PbNext;
                    Ph->Next = PhNext;
                    Pb->Baryon = (void *)Ph;
                    Ph->Body = (void *)Pb;
                    
                    Pb->Mass *= Inv9;
                    Ph->Mass *= Inv9;
#ifdef USE_MULTIPHASE_MODEL //{
                    Ph->Mhot *= Inv9;
#endif //USE_MULTIPHASE_MODEL //}
                    for(int l=0;l<CELibYield_Number;l++){
                        Ph->Elements[l] *= Inv9;
                    }

                    // Position updates
                    Pb->PosP[0] += dr[k+1][0]*Lscale;
                    Pb->PosP[1] += dr[k+1][1]*Lscale;
                    Pb->PosP[2] += dr[k+1][2]*Lscale;
                    Ph->PosP[0] = Pb->Pos[0] = Pb->PosP[0];
                    Ph->PosP[1] = Pb->Pos[1] = Pb->PosP[1];
                    Ph->PosP[2] = Pb->Pos[2] = Pb->PosP[2];

                    Ph->KernelPred = Ph->Kernel = hchild;

                    Ph->SplitGeneration ++;
                    Ph->SplitNthChildren = k+1;
                    Ph->SplitFlag = false; 

                    if(Ph->SplitGeneration == 1){
                        Ph->AncestorGlobalID = PhydroBody(i)->GlobalID;
                    } else {
                        Ph->AncestorGlobalID = Phydro[i]->AncestorGlobalID;
                    }
                    Ph->Body->GlobalID = PowderSnowIDGenerator();
                }

                Phydro[i]->Mass *= Inv9;
                PhydroBody(i)->Mass *= Inv9;
#ifdef USE_MULTIPHASE_MODEL //{
                Phydro[i]->Mhot *= Inv9;
#endif //USE_MULTIPHASE_MODEL //}
                for(int l=0;l<CELibYield_Number;l++){
                    Phydro[i]->Elements[l] *= Inv9;
                }

                // Position updates
                PhydroBody(i)->PosP[0] += dr[0][0]*Lscale;
                PhydroBody(i)->PosP[1] += dr[0][1]*Lscale;
                PhydroBody(i)->PosP[2] += dr[0][2]*Lscale;
                Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0] = PhydroBody(i)->PosP[0];
                Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1] = PhydroBody(i)->PosP[1];
                Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2] = PhydroBody(i)->PosP[2];


                Phydro[i]->KernelPred = Phydro[i]->Kernel = hchild;

                Phydro[i]->SplitGeneration ++;
                Phydro[i]->SplitNthChildren = 0;
                Phydro[i]->SplitFlag = false; 

                // Keep global ID
                if(Phydro[i]->SplitGeneration == 1){
                    Phydro[i]->AncestorGlobalID = Phydro[i]->GlobalID;
                } else {
                    Phydro[i]->AncestorGlobalID = Phydro[i]->AncestorGlobalID;
                }

                counter += 8;
            }
        }
    }

    Pall.Nhydro        += counter;
    Pall.NActivesHydro += counter;
    Pall.Ntotal        += counter;
    Pall.NActives      += counter;

    UpdateTotalNumber();
    ReConnectPointers();

    return ;
}

#endif // (PARTICLE_SPLITTING_TYPE == ?) //}


static void StretchDataArrayes(const int LocalSplittingTarget){

#if (PARTICLE_SPLITTING_TYPE == 0) //{
    const int BaseNumber = 12;
#elif (PARTICLE_SPLITTING_TYPE == 1) //}//{
    const int BaseNumber = 7;
#elif (PARTICLE_SPLITTING_TYPE == 2) //}//{
    const int BaseNumber = 8;
#endif // (PARTICLE_SPLITTING_TYPE == ?) //}

    int EmptyPbody = CountEmptyPbody();
    if(EmptyPbody < BaseNumber*LocalSplittingTarget){
        StretchStructPbody(BaseNumber*LocalSplittingTarget-EmptyPbody);
    }

    int EmptyPhydro = CountEmptyPhydro();
    if(EmptyPhydro < BaseNumber*LocalSplittingTarget){
        StretchStructPhydro(BaseNumber*LocalSplittingTarget-EmptyPhydro);
    }

    ReConnectPointers();

    return ;
}

static bool ChechSplitStep(void){

    if(Pall.TStep%PARTICLE_SPLITTING_INTERAVAL == 0){
        return true;
    } else if(Pall.NActivesHydro_t > 0.1*Pall.Nhydro_t){
        return true;
    } else {
        return false;
    }
}
#endif // USE_PARTICLE_SPLITTING //}

static bool FirstCall = true;

void ParticleSplitting(void){

#ifdef USE_PARTICLE_SPLITTING //{
    if(FirstCall){
        PowderSnowSetShardID(MPIGetMyID());
        FirstCall = false;
    }
    if(ChechSplitStep()){
        return;
    }

    int SplittingTarget = CheckSplitCondition();

    int LocalSplittingTarget = SplittingTarget;
    MPI_Allreduce(MPI_IN_PLACE,&SplittingTarget,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(SplittingTarget == 0) return ;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Particle Splitting! %d particles \n",SplittingTarget);
    }

    StretchDataArrayes(LocalSplittingTarget);

    GenChildParticles(LocalSplittingTarget);
#endif // USE_PARTICLE_SPLITTING //}

    return ;
}
