#include "config.h"
#include "PreDecomposition.h"
#include "Decomposition.h"
#include "ForceMisc.h"
#include "ForceGRAPE.h"
#include "ForceParallelTreeGRAPE.h"
#include "ForceFromExternalPotentials.h"
#include "PlantGravityTree.h"
#include "PlantHydroTree.h"
#include "Integral.h"
#include "TimeStep.h"
#include "HydroDensity.h"
#include "HydroAcc.h"
#include "HydroMisc.h"
#include "HydroKernel.h"
#include "HydroExtraOperations.h"
#include "TimeStep.h"
#include "Cooling.h"
#include "Heating.h"
#include "StarFormation.h"
#include "Delayed.h"
#include "StellarFeedback.h"
#include "HIIregion.h"
#include "SetUpTestRun.h"
#include "SinkParticle.h"
#include "SizeDetermination.h"
#include "ThermalConductivity.h"
#include "Logs.h"
#include "RunLogs.h"
#include "FUV.h"
#ifdef USE_RADIATION_PRESSURE //{
#include "RadiationPressure.h"
#endif // USE_RADIATION_PRESSURE //}
#include "ParticleSplitting.h"
#ifdef USE_STELLAR_WIND //{
#include "StellarWind.h"
#endif // USE_STELLAR_WIND //}

//#include <mcheck.h>


static void WriteData(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./phase.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(fabs(Phydro[i]->PosP[2]*Pall.UnitLength/KPC_CGS)<0.1)
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                PhydroBody(i)->GlobalID,
                Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho,
                Pall.ConvertUtoT*Phydro[i]->UPred,
                Phydro[i]->G0thin,Phydro[i]->GradRho,
                Phydro[i]->Z,Phydro[i]->GradN,
                //Phydro[i]->Gradh,Phydro[i]->GradN,
                Phydro[i]->fij,
                Phydro[i]->dt_hydro,
                Phydro[i]->dt_hydro*Pall.UnitTime/YEAR_CGS,
                PhydroBody(i)->dt*Pall.UnitTime/YEAR_CGS,
                PhydroBody(i)->PosP[0]*Pall.UnitLength/KPC_CGS,
                PhydroBody(i)->PosP[1]*Pall.UnitLength/KPC_CGS,
                PhydroBody(i)->PosP[2]*Pall.UnitLength/KPC_CGS,
                Phydro[i]->VelP[0],
                Phydro[i]->VelP[1],
                Phydro[i]->VelP[2]
                );
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./phase.??? | sort -n > ./phase.data");
        fflush(NULL);
    }
    fflush(NULL);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("rm -rf ./phase.???");
        fflush(NULL);
    }

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ;
}

static void WriteS(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./phase.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nstars;i++){
        fprintf(fp,"%ld %g %g %g %g %g\n",
                PstarBody(i)->GlobalID,
                PstarBody(i)->PosP[0]*Pall.UnitLength/KPC_CGS,
                PstarBody(i)->PosP[1]*Pall.UnitLength/KPC_CGS,
                PstarBody(i)->PosP[2]*Pall.UnitLength/KPC_CGS,
                Pstar[i]->Z,
                Pstar[i]->FormationTime*Pall.UnitTime/YEAR_CGS);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./phase.??? | sort -n > ./phase.data");
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        //system("rm -rf ./phase.???");
        fflush(NULL);
    }

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Finalize();

    return ;
}

static void WriteHydro(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./hydro.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active)
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                PhydroBody(i)->GlobalID,
                Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho, //2 
                Pall.ConvertUtoT*Phydro[i]->UPred,             //3
                //Phydro[i]->GradRho, 
                Phydro[i]->Gradh,                              //4
                Phydro[i]->GradN,                              //5
                Phydro[i]->fij,                                //6
                Phydro[i]->dt_hydro,                           //7
                Phydro[i]->dt_hydro*Pall.UnitTime/YEAR_CGS,    //8
                Phydro[i]->Alpha,                              //9
                Phydro[i]->DQheat,                             //10
                PhydroBody(i)->dt*Pall.UnitTime/YEAR_CGS,      //11
                PhydroBody(i)->PosP[0]*Pall.UnitLength/KPC_CGS,// 12
                PhydroBody(i)->PosP[1]*Pall.UnitLength/KPC_CGS,
                PhydroBody(i)->PosP[2]*Pall.UnitLength/KPC_CGS,
                PhydroBody(i)->Vel[0]*Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS, //15
                PhydroBody(i)->Vel[1]*Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS,
                PhydroBody(i)->Vel[2]*Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS, 
                Phydro[i]->HydroAcc[0], //16
                Phydro[i]->HydroAcc[1], 
                Phydro[i]->HydroAcc[2]
                );
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"cat ./hydro.??? | sort -n > ./hydro.%05d.data",Pall.TStepTotal);
        system(command);
        fflush(NULL);
    }
    fflush(NULL);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("rm -rf ./hydro.???");
        fflush(NULL);
    }

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();

    return ;
}

static void WriteTimeStepData(void);
static void WriteTimeStepDataWithType(char base[], const int Type);
static void WritePhaseData(void);

static int DecompFrequencyCounter = 0;

static void checker(const int GID, const int LINE){

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
        if(Pbody[i]->GlobalID == GID){
            fprintf(stderr,"-%d-P %g %g %g\n",LINE,Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2]); 
            fprintf(stderr,"-%d-V %g %g %g\n",LINE,Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]); 
            fprintf(stderr,"-%d-H %g %g %g\n",LINE,Pbody[i]->Velh[0],Pbody[i]->Velh[1],Pbody[i]->Velh[2]); 
            fprintf(stderr,"-%d-A %g %g %g\n",LINE,Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]); 
            fprintf(stderr,"-%d-HA %g %g %g\n",LINE, PbodyHydro(i)->HydroAcc[0], PbodyHydro(i)->HydroAcc[1], PbodyHydro(i)->HydroAcc[2]);
            fprintf(stderr,"-%d-Al,dt %g %g\n",LINE, PbodyHydro(i)->Alpha, PbodyHydro(i)->dt_hydro);

        }
        }
    }
    fflush(NULL);

    return ;
}

static void checker_u(const int GID, const int LINE){

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
        if(Pbody[i]->GlobalID == GID){
            fprintf(stderr,"-%d-P %g %g %g\n",LINE,Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2]); 
            fprintf(stderr,"-%d-V %g %g %g\n",LINE,Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]); 
            fprintf(stderr,"-%d-H %g %g %g\n",LINE,Pbody[i]->Velh[0],Pbody[i]->Velh[1],Pbody[i]->Velh[2]); 
            fprintf(stderr,"-%d-A %g %g %g\n",LINE,Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]); 
            fprintf(stderr,"-%d-HA %g %g %g\n",LINE, PbodyHydro(i)->HydroAcc[0], PbodyHydro(i)->HydroAcc[1], PbodyHydro(i)->HydroAcc[2]);
            fprintf(stderr,"-%d-UT %g %g\n",LINE, PbodyHydro(i)->UPred, Pall.ConvertUtoT*PbodyHydro(i)->UPred);
            fprintf(stderr,"-%d-Al,dt %g %g\n",LINE, PbodyHydro(i)->Alpha, PbodyHydro(i)->dt_hydro);

        }
        }
    }
    fflush(NULL);

    return ;
}

static void checker2(const int LINE){

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroBody(i)->GlobalID == 609620){
            CheckHydroStructures(i);
        }

        // check velocity

        double v = NORM(PhydroBody(i)->Vel)*(Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
        double vh = NORM(Phydro[i]->VelP)*(Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
        if(v > 1000){
            if(v > 10000){
                fprintf(stderr,"Vel[%ld] = %g km/s , %g km/s high speed! %d\n",PhydroBody(i)->GlobalID,v,vh,LINE);
            } else {
                fprintf(stderr,"Vel[%ld] = %g km/s , %g km/s \n",PhydroBody(i)->GlobalID,v,vh);
            }
        }


        const double TimeFactor = Pall.UnitTime/YEAR_CGS;
        const double LengthFactor = Pall.UnitLength/KPC_CGS;
        const double MassFactor = Pall.UnitMass/MSUN_CGS;
        const double VelFactor = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
        if(Phydro[i]->dt_hydro*TimeFactor < 4000){
            fprintf(stderr,"Small time step! GID %ld\n",PhydroBody(i)->GlobalID);

            fprintf(stderr,"%ld %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroBody(i)->dt*TimeFactor,
                    Phydro[i]->dt_hydro*TimeFactor,
                    Phydro[i]->Kernel*LengthFactor,
                    Phydro[i]->Vsig*VelFactor,
                    Phydro[i]->Kernel/Phydro[i]->Vsig*TimeFactor,
                    PhydroBody(i)->Mass*MassFactor,
                    PhydroBody(i)->Eps*LengthFactor);
        }
    }
    fflush(NULL);

    return ;
}

static void checker3(const int LINE){

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroBody(i)->GlobalID == 31168){
            CheckHydroStructures(i);

            // check velocity

            double v = NORM(PhydroBody(i)->Vel)*(Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
            double vh = NORM(Phydro[i]->VelP)*(Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
            if(v > 1000){
                if(v > 10000){
                    fprintf(stderr,"Vel[%ld] c3= %g km/s , %g km/s high speed! %d\n",PhydroBody(i)->GlobalID,v,vh,LINE);
                } else {
                    fprintf(stderr,"Vel[%ld] c3= %g km/s , %g km/s %d \n",PhydroBody(i)->GlobalID,v,vh,LINE);
                }
            }
        }

    }
    fflush(NULL);

    return ;
}

static void PrintMultiphase(const int Tag){

    for(int i=0;i<Pall.Nhydro;i++){
        if((Phydro[i]->Active)&&(Phydro[i]->MultiphaseFlag)){
            //if(PhydroBody(i)->GlobalID%1000 == 0){
            if(PhydroBody(i)->GlobalID == 347){
                fprintf(stderr,"M%d %ld %g %g %g %g %g %g %g %g %g\n",Tag,
                        PhydroBody(i)->GlobalID,
                        Pall.TCurrent*Pall.UnitTime/YEAR_CGS,
                        Phydro[i]->Mhot,Pall.ConvertUtoT*Phydro[i]->Uhot,Phydro[i]->Rhohot,
                        Phydro[i]->Mhot*Phydro[i]->Uhot,
                        Phydro[i]->Mass,Pall.ConvertUtoT*Phydro[i]->U,Phydro[i]->Rho,
                        Phydro[i]->Mass*Phydro[i]->U);
            }
        }
    }
    fflush(NULL);

    return ;
}

static void PrintMultiphaseEcheck(const int Tag){

    for(int i=0;i<Pall.Nhydro;i++){
        if((Phydro[i]->Active)&&(Phydro[i]->MultiphaseFlag)){
            double E = Phydro[i]->Mass*Phydro[i]->U;
            double Ehot = Phydro[i]->Mhot*Phydro[i]->Uhot;
            if(E<= Ehot){
                fprintf(stderr,"H%d %ld %g %g %g %g %g %g %g %g %g\n",Tag,
                        PhydroBody(i)->GlobalID,
                        Pall.TCurrent*Pall.UnitTime/YEAR_CGS,
                        Phydro[i]->Mhot,Pall.ConvertUtoT*Phydro[i]->Uhot,Phydro[i]->Rhohot,
                        Phydro[i]->Mhot*Phydro[i]->Uhot,
                        Phydro[i]->Mass,Pall.ConvertUtoT*Phydro[i]->U,Phydro[i]->Rho,
                        Phydro[i]->Mass*Phydro[i]->U);
            }
        }
    }
    fflush(NULL);

    return ;
}

#if 1
void Run(void){

#if 0
    //WriteTimeStepData();
    WriteTimeStepDataWithType("dt.h",TypeHydro);
    WriteTimeStepDataWithType("dt.s",TypeStar);
    WriteTimeStepDataWithType("dt.d",TypeDM);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
#endif


    while(Pall.TCurrent < Pall.TEnd){

        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"%u : ",Pall.TStepTotal);
            fflush(NULL);
        }

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        UpdateCosmologicalParameters();

        if(Pall.TCurrent >= Pall.Era){
            //PreDomainDecomposition();
            PreDomainDecomposition(0);
            DomainDecompositionOnlyDataExchange();

            SortAllStructures();

            FirstTimeStep();
            BuildHierarchicalTimeStep();
        } else {
            if(DecompFrequencyCounter > DECOMPOSITION_FREQUENCY){
                if (10*Pall.NActives_t > Pall.Ntotal_t){
                    PreDomainDecomposition(1);
                    DomainDecompositionOnlyDataExchange();
                    DecompFrequencyCounter = 0;
                }
            } else {
                DecompFrequencyCounter ++;
            }
            BuildNewTimeStep();
        }

        RaiseActiveFlags();

        Kick1Drift(); 
        BuildPredictors();


#ifdef GRAVITY_RUN //{
        if(Pall.NActives_t>0){ // Gravity
            PlantGravityTree();

            ClearGravitationalForce();
            ForceParallelTreeGRAPE(); // Split this function? insert PlantHydroTree and ClearHydroData?
            ForceFromExternalPotentials();
        } else {
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(stderr,"Skip Gravity\n");
        }
#endif // GRAVITY_RUN //}

#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            PlantHydroTree();

            ClearHydroData();
#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
            CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
            CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
            //PlantHydroTree(); // Mistake?
            CalcDensityDivRot(); 
            CalcDuDtAcc();
        }
#endif // HYDRO_RUN //}

        Kick2();

#ifdef HYDRO_RUN //{

        if((Pall.NActivesHydro_t==0)&&(Pall.NActivesStars_t+Pall.NActivesSink_t > 0)){
            PlantHydroTree();
#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
            CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
            CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
        }

#ifdef USE_CELIB //{
        //StellarFeedback();
        StellarFeedbackGatherScatterESA();
#else // USE_CELIB //}//{
        DelayedSNe();
#endif  // USE_CELIB //}

        if(Pall.NActivesHydro_t>0){ // Hydro

#ifdef USE_FUVFEEDBACK //{
            CalcFUV();
#endif // USE_FUVFEEDBACK //}
    // Update Particle predictors.
            // AdjustTimeStepAfterFeedback();

            CalcCooling();
            HIIregions();
            CalcDuDtAccEnergyDensityForCorrection();
            StarFormation();
            SinkParticles();
        }
#endif // HYDRO_RUN //}

        AddDragTerm();

        UpdateGravityKickFlag();

        // Calc New Time Step or Out Put Logs
        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();

        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            OutPutAllParticlesInASCIIFormat();

            if(Pall.Nhydro_t>0){
                CountNeighborNumber();
            }

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            // DataDump();
            FileOutPutConstantInterval();

            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            LogTotalMass();
            ReportAllocatedMemorySizes();
            EndLogs();
            fflush(NULL);
            MPI_Barrier(MPI_COMM_WORLD);

            // DataFullDump();

            // WriteData();

            ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }

        DataFullDump();
#ifdef USE_ON_THE_FLY_4D2U //{
#ifdef USE_ON_THE_FLY_4D2U_EXTRA_INFO //{
        Write4D2UExtraInfo();
#endif // USE_ON_THE_FLY_4D2U_EXTRA_INFO //}
#endif //USE_ON_THE_FLY_4D2U //}
    }

    return ;
}
#else 
//if(MPIGetMyID() == MPI_ROOT_RANK){ fprintf(stderr,"%s:%d\n",__func__,__LINE__); fflush(NULL); } MPI_Barrier(MPI_COMM_WORLD);
enum {
    TimeCosmo,
    TimePreDecomp,
    TimeDecomp,
    TimeSort,
    TimeFirstT,
    TimeBuildT,
    TimeNewT,
    TimeRaise,
    TimePredict,
    TimeKick1,
    TimeDrift,
    TimeKick2,
    TimeGravityAll,
    TimeGravKick,
    TimeGravity,
    TimeGravityTree,
    TimeGravityExt,
    TimeHydroAll,
    TimeHydroKernel,
    TimeHydroDensity,
    TimeHydroAcc,
    TimeHydroAccCor,
    TimeHydroTree,
    TimeOthersAll,
    TimeCooling,
    TimeSF,
    TimeFB,
    TimeHII,
    TimeSink,
    TimeIO,
    TimeAll,
    NTime,
};


#define Start(_x) { \
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime();\
    }

#define End(_x) {\
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime()-Times[_x];\
    }

#define Add(_x) {\
    TimesAll[_x] += Times[_x];\
    }


#define ShowResults {\
    MPI_Barrier(MPI_COMM_WORLD); \
    if(MPIGetMyID() == MPI_ROOT_RANK){ \
        fprintf(stderr,"All %g\n",Times[TimeAll]); \
        fprintf(stderr,"Cosm %g, PreDec %g, Dec %g, Sort %g\n",Times[TimeCosmo],Times[TimePreDecomp],Times[TimeDecomp],Times[TimeSort]); \
        fprintf(stderr,"First %g, Build %g, New %g, Rais %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"Kick1 %g, Kick2 %g, GravKick %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"GravTree %g, Grav %g, GravExt %g\n",Times[TimeGravityTree],Times[TimeGravity],Times[TimeGravityExt]); \
        fprintf(stderr,"HTree %g, HK %g, HD %g, HA %g, HAC %g\n",Times[TimeHydroTree],Times[TimeHydroKernel],Times[TimeHydroDensity],Times[TimeHydroAcc],Times[TimeHydroAccCor]); \
        fprintf(stderr,"Cooling %g, SF %g, FB %g, HII %g, Sink %g IO %g\n",Times[TimeCooling],Times[TimeSF],Times[TimeFB],Times[TimeHII],Times[TimeSink],Times[TimeIO]); \
        fflush(stderr); \
    }\
    }

void Run(void){

    // Start
    // End(TimeGravityTree)

    double Times[NTime];
    double TimesAll[NTime];
    for(int i=0;i<NTime;i++)
        TimesAll[i] = 0.e0;


    int counter = 0;
    while(Pall.TCurrent < Pall.TEnd){


        for(int i=0;i<NTime;i++)
            Times[i] = 0.e0;

        Start(TimeAll)

        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"%ld : ",Pall.TStepTotal);
            fflush(NULL);
        }

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();
        
        Start(TimeCosmo)
        UpdateCosmologicalParameters();
        End(TimeCosmo)

        if(Pall.TCurrent >= Pall.Era){
            Start(TimePreDecomp)
            //PreDomainDecomposition();
            PreDomainDecomposition(0);
            End(TimePreDecomp)

            Start(TimeDecomp)
            DomainDecompositionOnlyDataExchange();
            End(TimeDecomp)

            Start(TimeSort)
            SortAllStructures();
            End(TimeSort)

            Start(TimeFirstT)
            FirstTimeStep();
            End(TimeFirstT)

            Start(TimeBuildT)
            BuildHierarchicalTimeStep();
            End(TimeBuildT)

            HydroRoot.LifeGauge = 0;
        } else {

            if (10*Pall.NActives_t > Pall.Ntotal_t){
                Start(TimePreDecomp)
                //PreDomainDecomposition();
                PreDomainDecomposition(1);
                End(TimePreDecomp)

                Start(TimeDecomp)
                DomainDecompositionOnlyDataExchange();
                End(TimeDecomp)
                HydroRoot.LifeGauge = 0;
            }

            Start(TimeNewT)
            BuildNewTimeStep();
            End(TimeNewT)
        }

        Start(TimeRaise)
        RaiseActiveFlags();
        End(TimeRaise)

        Start(TimeKick1)
        Kick1Drift(); 
        End(TimeKick1)

        Start(TimePredict)
        BuildPredictors();
        End(TimePredict)

#ifdef GRAVITY_RUN //{
        Start(TimeGravityAll)
        if(Pall.NActives_t>0){ // Gravity

            Start(TimeGravityTree)
            PlantGravityTree();
            End(TimeGravityTree)

            ClearGravitationalForce();

            Start(TimeGravity)
            ForceParallelTreeGRAPE();
            End(TimeGravity)

            Start(TimeGravityExt)
            ForceFromExternalPotentials();
            End(TimeGravityExt)

        } else {
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(stderr,"Skip Gravity\n");
        }
        End(TimeGravityAll)
        Add(TimeGravityAll)
#endif // GRAVITY_RUN //}


#ifdef HYDRO_RUN //{
        Start(TimeHydroAll)
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeHydroTree)
            PlantHydroTree();
            End(TimeHydroTree)

            ClearHydroData();
            Start(TimeHydroKernel)
#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
            CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
            CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
            End(TimeHydroKernel)

            Start(TimeHydroDensity)
            CalcDensityDivRot(); 
            End(TimeHydroDensity)

            Start(TimeHydroAcc)
            CalcDuDtAcc();
            End(TimeHydroAcc)
        }
        End(TimeHydroAll)
        Add(TimeHydroAll)
#endif // HYDRO_RUN //}

        Start(TimeKick2)
        Kick2();
        End(TimeKick2)

        Start(TimeOthersAll)
#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeFB)
#ifdef USE_CELIB //{
            StellarFeedback();
#else
            DelayedSNe();
#endif  // USE_CELIB //}
            End(TimeFB)

            Start(TimeCooling)
            CalcCooling();
            End(TimeCooling)

            Start(TimeHII)
            HIIregions();
            End(TimeHII)

            Start(TimeHydroAccCor)
            CalcDuDtAccEnergyDensityForCorrection();
            End(TimeHydroAccCor)

            Start(TimeSF)
            StarFormation();
            End(TimeSF)

            Start(TimeSink)
            SinkParticles();
            End(TimeSink)
        }
#endif // HYDRO_RUN //}
        End(TimeOthersAll)
        Add(TimeOthersAll)

        Start(TimeGravKick)
        UpdateGravityKickFlag();
        End(TimeGravKick)

        /// Post output
        // WriteCurrentActiveParticles(Pall.TStepTotal,"Post");

        // Calc New Time Step or Out Put Logs
        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        LogsThisTimeStep();

        Start(TimeIO)
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            // OutPutAllParticlesInASCIIFormat();

            if(Pall.Nhydro_t>0){
                CountNeighborNumber();
            }

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            // DataDump();
            // FileOutPutConstantInterval();

            // LogOutPutEnergyMomentumAngularMomentum();
            // LogOutPutElapsedTime();
            // LogTotalMass();
            // ReportAllocatedMemorySizes();
            // EndLogs();
            // fflush(NULL);
            // MPI_Barrier(MPI_COMM_WORLD);

            // DataFullDump();

            ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }

        // DataFullDump();

        End(TimeIO)
        End(TimeAll)
        Add(TimeAll)

        ShowResults;

        if(counter > 30){

            if(MPIGetMyID() == MPI_ROOT_RANK){
                FILE *fp;
                char fname[MaxCharactersInLine];
                Snprintf(fname,"Bench.%04d",MPIGetNumProcs());
                FileOpen(fp,fname,"w");
                fprintf(fp,"%04d %g %g %g %g %g\n",MPIGetNumProcs(),TimesAll[TimeAll],
                TimesAll[TimeGravityAll]+TimesAll[TimeHydroAll]+TimesAll[TimeOthersAll],
                TimesAll[TimeGravityAll],TimesAll[TimeHydroAll],TimesAll[TimeOthersAll]);
                fflush(NULL);
            }


            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        counter ++;
    }

    return ;
}


#endif

#if 0
enum {
    TimeCosmo,
    TimePreDecomp,
    TimeDecomp,
    TimeSort,
    TimeFirstT,
    TimeBuildT,
    TimeNewT,
    TimeRaise,
    TimePredict,
    TimeKick1,
    TimeDrift,
    TimeKick2,
    TimeGravKick,
    TimeGravity,
    TimeGravityTree,
    TimeGravityExt,
    TimeHydroKernel,
    TimeHydroDensity,
    TimeHydroAcc,
    TimeHydroAccCor,
    TimeHydroTree,
    TimeCooling,
    TimeSF,
    TimeFB,
    TimeHII,
    TimeSink,
    TimeIO,
    TimeAll,
    NTime,
};

#define Start(_x) { \
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime();\
    }

#define End(_x) {\
    MPI_Barrier(MPI_COMM_WORLD); \
    Times[_x] = GetElapsedTime()-Times[_x];\
    }

#define ShowResults {\
    MPI_Barrier(MPI_COMM_WORLD); \
    if(MPIGetMyID() == MPI_ROOT_RANK){ \
        fprintf(stderr,"All %g\n",Times[TimeAll]); \
        fprintf(stderr,"Cosm %g, PreDec %g, Dec %g, Sort %g\n",Times[TimeCosmo],Times[TimePreDecomp],Times[TimeDecomp],Times[TimeSort]); \
        fprintf(stderr,"First %g, Build %g, New %g, Rais %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"Kick1 %g, Kick2 %g, GravKick %g, Pred %g\n",Times[TimeFirstT],Times[TimeBuildT],Times[TimeNewT],Times[TimeRaise]); \
        fprintf(stderr,"GravTree %g, Grav %g, GravExt %g\n",Times[TimeGravityTree],Times[TimeGravity],Times[TimeGravityExt]); \
        fprintf(stderr,"HTree %g, HK %g, HD %g, HA %g, HAC %g\n",Times[TimeHydroTree],Times[TimeHydroKernel],Times[TimeHydroDensity],Times[TimeHydroAcc],Times[TimeHydroAccCor]); \
        fprintf(stderr,"Cooling %g, SF %g, FB %g, HII %g, Sink %g IO %g\n",Times[TimeCooling],Times[TimeSF],Times[TimeFB],Times[TimeHII],Times[TimeSink],Times[TimeIO]); \
        fflush(stderr); \
    }\
    }

void Run(void){

    // Start
    // End(TimeGravityTree)

    double Times[NTime];
    double ___t;

    int counter = 0;
    while(Pall.TCurrent < Pall.TEnd){

        for(int i=0;i<NTime;i++)
            Times[i] = 0.e0;

        Start(TimeAll)

        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"%ld : ",Pall.TStepTotal);
            fflush(NULL);
        }

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();
        
        Start(TimeCosmo)
        UpdateCosmologicalParameters();
        End(TimeCosmo)

        if(Pall.TCurrent >= Pall.Era){
            Start(TimePreDecomp)
            PreDomainDecomposition();
            End(TimePreDecomp)

            Start(TimeDecomp)
            DomainDecompositionOnlyDataExchange();
            End(TimeDecomp)

            Start(TimeSort)
            SortAllStructures();
            End(TimeSort)

            Start(TimeFirstT)
            FirstTimeStep();
            End(TimeFirstT)

            Start(TimeBuildT)
            BuildHierarchicalTimeStep();
            End(TimeBuildT)

            HydroRoot.LifeGauge = 0;
        } else {

            if (10*Pall.NActives_t > Pall.Ntotal_t){
                Start(TimePreDecomp)
                PreDomainDecomposition();
                End(TimePreDecomp)

                Start(TimeDecomp)
                DomainDecompositionOnlyDataExchange();
                End(TimeDecomp)
                HydroRoot.LifeGauge = 0;
            }

            Start(TimeNewT)
            BuildNewTimeStep();
            End(TimeNewT)
        }

        Start(TimeRaise)
        RaiseActiveFlags();
        End(TimeRaise)

        Start(TimeKick1)
        Kick1Drift(); 
        End(TimeKick1)

        Start(TimePredict)
        BuildPredictors();
        End(TimePredict)

        /// Pre output
        // WriteCurrentActiveParticles(Pall.TStepTotal,"Pre");


#ifdef GRAVITY_RUN //{
        if(Pall.NActives_t>0){ // Gravity


            Start(TimeGravityTree)
            PlantGravityTree();
            End(TimeGravityTree)

            ClearGravitationalForce();

            Start(TimeGravity)
            ForceParallelTreeGRAPE();
            End(TimeGravity)

            Start(TimeGravityExt)
            ForceFromExternalPotentials();
            End(TimeGravityExt)

        } else {
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(stderr,"Skip Gravity\n");
        }
#endif // GRAVITY_RUN //}


#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeHydroTree)
            PlantHydroTree();
            End(TimeHydroTree)

            ClearHydroData();
            Start(TimeHydroKernel)
            CalcKernel();
            End(TimeHydroKernel)

            Start(TimeHydroDensity)
            CalcDensityDivRot(); 
            End(TimeHydroDensity)

            Start(TimeHydroAcc)
            CalcDuDtAcc();
            End(TimeHydroAcc)
        }
#endif // HYDRO_RUN //}

        Start(TimeKick2)
        Kick2();
        End(TimeKick2)

#ifdef HYDRO_RUN //{
        if(Pall.NActivesHydro_t>0){ // Hydro
            Start(TimeFB)
#ifdef USE_CELIB //{
            StellarFeedback();
#else
            DelayedSNe();
#endif  // USE_CELIB //}
            End(TimeFB)

            Start(TimeCooling)
            CalcCooling();
            End(TimeCooling)

            Start(TimeHII)
            HIIregions();
            End(TimeHII)

            Start(TimeHydroAccCor)
            CalcDuDtAccEnergyDensityForCorrection();
            End(TimeHydroAccCor)

            Start(TimeSF)
            StarFormation();
            End(TimeSF)

            Start(TimeSink)
            SinkParticles();
            End(TimeSink)
        }
#endif // HYDRO_RUN //}

        Start(TimeGravKick)
        UpdateGravityKickFlag();
        End(TimeGravKick)

        /// Post output
        // WriteCurrentActiveParticles(Pall.TStepTotal,"Post");

        // Calc New Time Step or Out Put Logs
        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        //LogsThisTimeStep();

        Start(TimeIO)
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            OutPutAllParticlesInASCIIFormat();

            if(Pall.Nhydro_t>0){
                CountNeighborNumber();
            }

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            // DataDump();
            FileOutPutConstantInterval();

            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            LogTotalMass();
            ReportAllocatedMemorySizes();
            EndLogs();
            fflush(NULL);
            MPI_Barrier(MPI_COMM_WORLD);

            // DataFullDump();

            ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }

        DataFullDump();

        End(TimeIO)
        End(TimeAll)

        ShowResults;

        if(counter > 10){
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        counter ++;
    }

    return ;
}
#endif 

#ifdef USE_SPSPH //{
static void ResetInitialZ(void){
    
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->Zw = Phydro[i]->ZwPred = PhydroBody(i)->Mass/Phydro[i]->Rho;
        Phydro[i]->PseudoDensity = Phydro[i]->PseudoDensityPred = 1.e0;
    }

    CalcDensityDivRot();

    return ;
}
#endif // USE_SPSPH //}

#if VISCOSITY_TYPE == 1 //{
void InitCopyHydroPredictors(void){

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->RhoPred = Phydro[i]->Rho;
        Phydro[i]->KernelPred = Phydro[i]->Kernel;
#ifdef USE_DISPH //{
        Phydro[i]->EnergyDensityPred = Phydro[i]->EnergyDensity;
#endif //USE_DISPH //}
#ifdef USE_SPSPH //{
        Phydro[i]->PseudoDensityPred;
        Phydro[i]->ZwPred = Phydro[i]->Zw;
#endif // USE_SPSPH //}
    }
    return ;
}
#endif // VISCOSITY_TYPE == 1 //}

//int NNN;
void InitializeRun(void){

    //Kick1Drift();
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    InitializeDecomposition();
    DomainDecomposition();

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

#ifdef GRAVITY_RUN
    PlantGravityTree();
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();
    ForceFromExternalPotentials();
#endif // GRAVITY_RUN

    if(Pall.Nhydro_t > 0){
        InitializeRootForHydro();
#ifdef COOLING_RUN //{
        InitializeCoolingTable();
#endif // COOLING_RUN //}
        InitializeFarUltraVioletField();
#ifdef USE_THERMAL_CONDUCTIVITY //{ 
        InitThermalConductivity();
#endif // USE_THERMAL_CONDUCTIVITY //}

        ASRFLXRunParameters.IMFType = CHEMICALEVOLUTION_IMFTYPE;
// #ifdef USE_FUVFEEDBACK //{
#if CHEMICALEVOLUTION_POPIII_IMF == 1 //{
        ASRFLXRunParameters.UsePopIII = true;
        // Use default values for PopIII
#endif // CHEMICALEVOLUTION_POPIII_IMF //}
        ASRFLXInitFUV();
// #endif // USE_FUVFEEDBACK //}
#ifdef USE_HIIREGION_MODEL //{
        ASRFLXInitNLY();
#endif // USE_HIIREGION_MODEL //}
#ifdef USE_RADIATION_PRESSURE //{
        InitRadiationPressure();
        ASRFLXInitBol();
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
        InitStellarWind();
        //ASRFLXInitBol();
#endif // USE_STELLAR_WIND //}

        InitializeStarFormation();
#ifdef USE_CELIB
        InitializeStellarFeedback();
#endif
#ifdef DELAYED_FEEDBACK
#error do not use this routine. (Sendrecv is not fit to the new domain decomposition)
        InitializeDelayedSNII();
#endif

        ClearHydroData();
        PlantHydroTree();

#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
#if 0
#if VISCOSITY_TYPE == 1 //{  
        CalcKernel();
        CalcDensityDivRot();
        InitCopyHydroPredictors();
#endif // VISCOSITY_TYPE == 1 //}
#endif 
        CalcSize();

#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
        CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}

        //CalcKernel();
        CalcDensityDivRot();
#ifdef USE_SPSPH //{
        ResetInitialZ();
#endif // USE_SPSPH //}
        CalcDuDtAcc();
    }

#ifdef USE_SYMMETRIZED_SOFTENING
    //CalcSymmetrizedPotential();
#endif
    LogOutPutEnergyMomentumAngularMomentum();

    FirstTimeStep();
    if(Pall.Nhydro_t > 0){
        PlantHydroTree();
#ifdef HYDRO_TIMESTEP_LIMITER 
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;
        }
        BuildHierarchicalTimeStep();
        CalcDensityDivRot();
        for(int i=0;i<Pall.Nhydro;i++){
            //Phydro[i]->k_hydro = MIN(Phydro[i]->k_hydro,Phydro[i]->k_hydro_localmin+MAX_K_LOCAL);
            Phydro[i]->k_hydro = MIN(Phydro[i]->k_hydro,Phydro[i]->k_hydro_localmin);
            Phydro[i]->dt_hydro = Pall.dtmin*exp2(Phydro[i]->k_hydro);
#ifndef USE_FAST_SCHEME 
            PhydroBody(i)->k = Phydro[i]->k_hydro;
            PhydroBody(i)->dt = Phydro[i]->dt_hydro;
#else
            if(Phydro[i]->k_hydro > PhydroBody(i)->k){
                Phydro[i]->k_hydro = PhydroBody(i)->k;
                Phydro[i]->dt_hydro = PhydroBody(i)->dt;
            }
#endif
        }
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;
        }
#endif
    }


    UpdateGravityKickFlag();


#ifdef USE_ON_THE_FLY_4D2U //{
    InitWrite4D2U();
#endif //USE_ON_THE_FLY_4D2U //}

    FileOutPutConstantInterval();
    // Pall.OutPutFileNumber ++;

#ifdef WRITE_KERNEL_SHAPE //{
    WriteKernelShape();
#endif //WRITE_KERNEL_SHAPE //}


    return ;
}

void RestartRun(void){

    // Decomposition 
    InitializeDecomposition();
    // PreDomainDecomposition();

    // First force / hydro calculation.
#ifdef GRAVITY_RUN
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

#endif // GRAVITY_RUN

    if(Pall.Nhydro_t > 0){
        InitializeRootForHydro();
        InitializeCoolingTable();
        InitializeFarUltraVioletField();
        InitializeStarFormation();
#ifdef USE_CELIB
        InitializeStellarFeedback();
#else
        InitializeDelayedSNII();
#endif
        ASRFLXRunParameters.IMFType = CHEMICALEVOLUTION_IMFTYPE;
#ifdef USE_FUVFEEDBACK //{
#if CHEMICALEVOLUTION_POPIII_IMF == 1 //{
        ASRFLXRunParameters.UsePopIII = true;
        // Use default values for PopIII
#endif // CHEMICALEVOLUTION_POPIII_IMF //}
        ASRFLXInitFUV();
#endif // USE_FUVFEEDBACK //}
#ifdef USE_HIIREGION_MODEL //{
        ASRFLXRunParameters.IMFType = CHEMICALEVOLUTION_IMFTYPE;
        ASRFLXInitNLY();
#endif // USE_HIIREGION_MODEL //}
#ifdef USE_RADIATION_PRESSURE //{
        InitRadiationPressure();
        ASRFLXInitBol();
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
        InitStellarWind();
        //ASRFLXInitBol();
#endif // USE_STELLAR_WIND //}
    }

    LogOutPutEnergyMomentumAngularMomentum();
    UpdateGravityKickFlag();

    ImposeBoundaryCondition(BOUNDARY_CONDITION_TYPE);

#ifdef USE_ON_THE_FLY_4D2U //{
    InitWrite4D2U();
    ResetWrite4D2UExtraCounter();
#endif //USE_ON_THE_FLY_4D2U //}

    return ;
}


static void WriteTimeStepData(void){

    const double TimeFactor = Pall.UnitTime/YEAR_CGS;
    const double LengthFactor = Pall.UnitLength/KPC_CGS;
    const double MassFactor = Pall.UnitMass/MSUN_CGS;
    const double VelFactor = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./dt.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){

        double ParticleEpsSize = Pbody[i]->Eps*Pall.AdaptiveSofteningFactor;
        double Acc[3] = {Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]};
        double normAcc = NORM(Acc);
        double DiffAcc[3] = {Acc[0]-Pbody[i]->AccOld[0],
                             Acc[1]-Pbody[i]->AccOld[1],
                             Acc[2]-Pbody[i]->AccOld[2]};
        double normDiffAcc = NORM(DiffAcc);
        double dt_diffa = TFactorDiffAcc*normAcc/normDiffAcc*Pbody[i]->dt;
        double dt_a = TFactorAcc*sqrt(ParticleEpsSize/normAcc);

        double dt_h = 0.01/Pall.HubbleZ;


        if(Pbody[i]->Type == TypeHydro){
            fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->dt*TimeFactor,
                    dt_diffa*TimeFactor,
                    dt_a*TimeFactor,
                    dt_h*TimeFactor,
                    Pbody[i]->Mass*MassFactor,
                    Pbody[i]->Eps*LengthFactor,
                    ParticleEpsSize*LengthFactor,
                    PbodyHydro(i)->dt_hydro*TimeFactor,
                    PbodyHydro(i)->Kernel*LengthFactor,
                    PbodyHydro(i)->Vsig*VelFactor,
                    PbodyHydro(i)->Kernel/PbodyHydro(i)->Vsig*TimeFactor
                    );
        } else {
            fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->dt*TimeFactor,
                    dt_diffa*TimeFactor,
                    dt_a*TimeFactor,
                    dt_h*TimeFactor,
                    Pbody[i]->Mass*MassFactor,
                    Pbody[i]->Eps*LengthFactor,
                    ParticleEpsSize*LengthFactor
                    );
        }
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./dt.??? | sort -n > ./dt.data");
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("rm -rf ./dt.???");
        fflush(NULL);
    }

    return ;
}

static void WriteTimeStepDataWithType(char base[], const int Type){

    const double TimeFactor = Pall.UnitTime/YEAR_CGS;
    const double LengthFactor = Pall.UnitLength/KPC_CGS;
    const double MassFactor = Pall.UnitMass/MSUN_CGS;
    const double VelFactor = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./%s.%03d",base,MPIGetMyID());
    FileOpen(fp,fname,"w");

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == Type){

            double ParticleEpsSize = Pbody[i]->Eps*Pall.AdaptiveSofteningFactor;
            double Acc[3] = {Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]};
            double normAcc = NORM(Acc);
            double dt_a = TFactorAcc*sqrt(ParticleEpsSize/normAcc);

            if(Pbody[i]->Type == TypeHydro){
                fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                        Pbody[i]->dt*TimeFactor,
                        dt_a*TimeFactor,
                        PbodyHydro(i)->dt_hydro*TimeFactor,
                        PbodyHydro(i)->Kernel*LengthFactor,
                        PbodyHydro(i)->Vsig*VelFactor,
                        PbodyHydro(i)->Kernel/PbodyHydro(i)->Vsig*TimeFactor,
                        Pbody[i]->Mass*MassFactor,
                        Pbody[i]->Eps*LengthFactor,
                        ParticleEpsSize*LengthFactor
                        );
            } else if(Pbody[i]->Type == TypeStar){
#ifdef USE_CELIB //{
                const double FactorMass = Pall.UnitMass/MSUN_CGS;
                const double FactorTime = YEAR_CGS/Pall.UnitTime;

                double dt_fb = fmax(CHEMICALEVOLUTION_SNII_TIMEINTERVAL_MIN*FactorTime,
                    2.0*
                    (CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = FactorMass*PbodyStar(i)->InitialMass,
                        .Metallicity = PbodyStar(i)->Z,
                        .Count = PbodyStar(i)->SNIICount,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII)-
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = FactorMass*PbodyStar(i)->InitialMass,
                        .Metallicity = PbodyStar(i)->Z,
                        .Count = PbodyStar(i)->SNIICount,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII))
                        *FactorTime);


                if(PbodyStar(i)->EventTimeSNII > Pall.TEnd){

                }

                fprintf(fp,"%ld %g %g %g %g %g %d %g %d %g %g %g %g %g\n",Pbody[i]->GlobalID,
                        Pbody[i]->dt*TimeFactor,
                        dt_a*TimeFactor,
                        //PbodyStar(i)->dt_fb*TimeFactor,
                        dt_fb*TimeFactor, //4

                        CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = FactorMass*PbodyStar(i)->InitialMass,
                        .Metallicity = PbodyStar(i)->Z,
                        .Count = PbodyStar(i)->SNIICount,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),
                        CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = FactorMass*PbodyStar(i)->InitialMass,
                        .Metallicity = PbodyStar(i)->Z,
                        .Count = PbodyStar(i)->SNIICount,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),
                        PbodyStar(i)->SNIICount, //7
                        PbodyStar(i)->EventTimeSNII*TimeFactor, //8
                        PbodyStar(i)->SNIaCount, //9
                        PbodyStar(i)->EventTimeSNIa*TimeFactor,
                        Pbody[i]->Mass*MassFactor, //11 
                        PbodyStar(i)->InitialMass*MassFactor, //12
                        Pbody[i]->Eps*LengthFactor,
                        ParticleEpsSize*LengthFactor
                        );
#endif //USE_CELIB //}
            } else {
                fprintf(fp,"%ld %g %g %g %g %g\n",Pbody[i]->GlobalID,
                        Pbody[i]->dt*TimeFactor,
                        dt_a*TimeFactor,
                        Pbody[i]->Mass*MassFactor,
                        Pbody[i]->Eps*LengthFactor,
                        ParticleEpsSize*LengthFactor
                        );
            }
        }
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"cat ./%s.??? | sort -n > ./%s.data",base,base);
        //system("cat ./dt.??? | sort -n > ./dt.data");
        system(command);
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"rm -rf ./%s.???",base);
        system(command);
        fflush(NULL);
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(TypeStar == Type){
            double Mass = 4.e5;
            int counter = Mass*0.005*10;

            FILE *fp;
            char fname[] = "./sne.data";
            FileOpen(fp,fname,"w");

            for(int i=0;i<counter;i++){

                fprintf(fp,"%d %g %g %g\n",i,
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                    .R = 0.e0,
                    .InitialMass_in_Msun = Mass,
                    //.Metallicity = 0.013,
                    .Metallicity = 0.0,
                    .Count = i,
                    .Mode = CELibSNIIRateModelID_Individual,
                    },CELibFeedbackType_SNII),
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                    .R = 1.e0,
                    .InitialMass_in_Msun = Mass,
                    //.Metallicity = 0.013,
                    .Metallicity = 0.0,
                    .Count = i,
                    .Mode = CELibSNIIRateModelID_Individual,
                    },CELibFeedbackType_SNII),

                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                    .R = 1.e0,
                    .InitialMass_in_Msun = Mass,
                    //.Metallicity = 0.013,
                    .Metallicity = 0.0,
                    .Count = i,
                    .Mode = CELibSNIIRateModelID_Individual,
                    },CELibFeedbackType_SNII)-
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                    .R = 0.e0,
                    .InitialMass_in_Msun = Mass,
                    .Metallicity = 0.0,
                    .Count = i,
                    .Mode = CELibSNIIRateModelID_Individual,
                    },CELibFeedbackType_SNII)

                );
            }
            fclose(fp);
        }
    }

    return ;
}


static void WritePhaseData(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./phase.%03d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                PhydroBody(i)->GlobalID,
                Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho,
                Pall.ConvertUtoT*Phydro[i]->UPred,
                Phydro[i]->KernelPred,Phydro[i]->Vsig,
                Phydro[i]->Gradh,Phydro[i]->GradN,
                Phydro[i]->fij,
                Phydro[i]->dt_hydro,
                Phydro[i]->dt_hydro*Pall.UnitTime/YEAR_CGS,
                PhydroBody(i)->dt*Pall.UnitTime/YEAR_CGS,
                PhydroBody(i)->PosP[0]*Pall.UnitLength/KPC_CGS,
                PhydroBody(i)->PosP[1]*Pall.UnitLength/KPC_CGS,
                PhydroBody(i)->PosP[2]*Pall.UnitLength/KPC_CGS,
                Phydro[i]->VelP[0],
                Phydro[i]->VelP[1],
                Phydro[i]->VelP[2]
                );
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./phase.??? | sort -n > ./phase.data");
        fflush(NULL);
    }
    fflush(NULL);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("rm -rf ./phase.???");
        fflush(NULL);
    }

    return ;
}
