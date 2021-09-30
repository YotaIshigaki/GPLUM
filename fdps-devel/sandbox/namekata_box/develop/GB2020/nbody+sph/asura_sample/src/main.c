#include "config.h"
#include "Run.h"
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
#include "NeighborSearch.h"
#include "TimeStep.h"
#include "Cooling.h"
#include "Heating.h"
#include "StarFormation.h"
#include "Delayed.h"
#include "SetUpTestRun.h"
#include "TestRuns.h"
#include "Read.h"
#include "CommunicationTable.h"
#include "Logs.h"
#include "FOF.h"
#include "RunLogs.h"
#include "Exit.h"
#include "ThermalConductivity.h"

#define Iter 3

int main_IOTest(const int argc, char *argv[]);
int main_Orbit(void);
int main_3DShockTube(void);
int main_SelfSimilarCooling(void);
int main_IsothermalSphericalCollapse(void);
int main_BreakingDam(void);
int main_ColdCollapseTest_ActiveParticles(void);
int main_NavarroWhiteTest(void);
int main_WadaNorman(const int argc, char *argv[]);
int main_TreeTest(void);
int main_ForceHydroTest(void);
int main_CosmologicalRun(const int argc, char *argv[]);
int main_NeighborSearchTest(void);
int main_ReadAndPlantTreeTest(const int argc, char *argv[]);
int main_BenchMarkTest(void);
int main_BenchMarkTestForForce(void);
int main_HydroCheckDirectAndTree(const int argc, char *argv[]);
int main_NeighborSearchTest2(void);

#ifdef TASK_MERGER
int main_IsolateDisk(const int argc, char *argv[]);
#endif //TASK_MERGER
#ifdef TASK_DICE_RUN
int main_DiceRun(const int argc, char *argv[]);
#endif //TASK_DICE_RUN
#ifdef TASK_MW
int main_MilkywayModel(const int argc, char *argv[]);
#endif //TASK_MW
#ifdef TASK_NFW
int main_NFWrun(const int argc, char *argv[]);
#endif
#ifdef TASK_CLOUD_EXPLOSION
int main_CloudExplosion(const int argc, char *argv[]);
int main_MakeInitialConditionForCloudExplosion(const int argc, char *argv[]);
#endif //TASK_CLOUD_EXPLOSION
#ifdef TASK_3D_COLLAPSE
int main_3DCollapse(const int argc, char *argv[]);
#endif //TASK_3D_COLLAPSE
#ifdef TASK_BLAST_WAVE
int main_BlastWave(void);
#endif //TASK_BLAST_WAVE
#ifdef TASK_SINUSOIDAL_WAVE
int main_SinusoidalWave(void);
#endif //TASK_SINUSOIDAL_WAVE
#ifdef TASK_M2_COLLAPSE
int main_M2SphericalCollapse(const int argc, char *argv[]);
#endif //TASK_B_COLLAPSELAST_WAVE
#ifdef TASK_TURBULENCE
int main_turbulence(const int argc, char *argv[]);
#endif //TASK_TURBULENCE
#ifdef TASK_ROTATINGDISK_WITH_SINK
int main_rotatingdisk_with_sink(const int argc, char *argv[]);
#endif //TASK_ROTATINGDISK_WITH_SINK
#ifdef TASK_AGNTORUS
int main_AGNTorus(const int argc, char *argv[]);
#endif //TASK_AGNTORUS
#ifdef TASK_COLD_COLLAPSE
int main_ColdCollapseTest(const int argc, char *argv[]);
#endif //TASK_COLD_COLLAPSE
#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
int main_TestSymmetrizedPotentialError(const int argc, char *argv[]);
int main_TestSymmetrizedPotentialError2(const int argc, char *argv[]);
int main_TestSymmetrizedPotentialError_ColdCollapse(const int argc, char *argv[]);
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
#ifdef TASK_1D_SHOCKE_TUBE
int main_ShockTube(void);
#endif //TASK_1D_SHOCKE_TUBE
#ifdef TASK_HYDROSTATIC
int main_HydroStatic(const int argc, char *argv[]);
#endif // TASK_HYDROSTATIC
#ifdef TASK_KELVINHELMHOLTZ_INSTABILITY
int main_KHInstability(const int argc, char *argv[]);
#endif //TASK_KELVINHELMHOLTZ_INSTABILITY
#ifdef TASK_GALACTIC_CENTER
int main_GalacticCenter(const int argc, char *argv[]);
#endif //TASK_GALACTIC_CENTER
#ifdef TASK_SANTABARBARA
int main_SantaBarbara(const int argc, char *argv[]);
#endif //TASK_SANTABARBARA
#ifdef TASK_TEST_NEIGHBORSEARCH
int main_Test_NeighborSearch(const int argc, char *argv[]);
#endif // TASK_TEST_NEIGHBORSEARCH
#ifdef TASK_1D_TWOFLUIDS
int main_TwoFluids(void);
#endif //TASK_TWOFLUIDS
#ifdef TASK_KEPLER
int main_Kepler(void);
#endif //TASK_KEPLER
#ifdef TASK_GALAXY_FORMATION //{
int main_GalaxyFormation(const int argc, char *argv[]);
#endif //TASK_GALAXY_FORMATION //}

int main(int argc, char *argv[]){

    InitializeMPIEnv(&argc,argv);
    GetRunStatus(argc,argv);
    atexit(ASURA_Exit);

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"This binary was compiled at %s on %s\n",__TIME__,__DATE__);

    fprintf(stderr,"Process %d on %s\n",MPIGetMyID(),MPIGetProcessorName());

    WriteMakeflags();

#ifdef TASK_CLOUD_EXPLOSION
    //return main_MakeInitialConditionForCloudExplosion(argc,argv);
    return main_CloudExplosion(argc,argv);
#endif
#ifdef TASK_BLAST_WAVE
    return main_BlastWave();
#endif

#ifdef TASK_MERGER
    return main_IsolateDisk(argc,argv);
#endif
#ifdef TASK_DICE_RUN
    return main_DiceRun(argc,argv);
#endif //TASK_DICE_RUN
#ifdef TASK_MW
    return main_MilkywayModel(argc,argv);
#endif

#ifdef TASK_3D_COLLAPSE
    return main_3DCollapse(argc,argv);
#endif
#ifdef TASK_NFW
    return main_NFWrun(argc,argv);
#endif
#ifdef TASK_1D_SHOCKE_TUBE
    return main_ShockTube();
#endif
#ifdef TASK_M2_COLLAPSE
    return main_M2SphericalCollapse(argc,argv);
#endif

#ifdef TASK_SINUSOIDAL_WAVE
    return main_SinusoidalWave();
#endif //TASK_SINUSOIDAL_WAVE

#ifdef TASK_TURBULENCE
    return main_turbulence(argc,argv);
#endif
#ifdef TASK_ROTATINGDISK_WITH_SINK
    return main_rotatingdisk_with_sink(argc,argv);
#endif //TASK_ROTATINGDISK_WITH_SINK
#ifdef TASK_AGNTORUS
    return main_AGNTorus(argc,argv);
#endif
#ifdef TASK_COLD_COLLAPSE
    return main_ColdCollapseTest(argc,argv);
#endif // TASK_COLD_COLLAPSE

#ifdef TASK_TEST_STROMGRENSPHERE
    return main_TestStromgrenSphere(argc,argv);
#endif //TASK_TEST_STROMGRENSPHERE

#ifdef TASK_TEST_RADIATION_PRESSURE
    return main_TestRadiationPressure(argc,argv);
#endif //TASK_TEST_RADIATION_PRESSURE

#ifdef TASK_HYDROSTATIC
    return main_HydroStatic(argc,argv);
#endif // TASK_HYDROSTATIC

#ifdef TASK_KELVINHELMHOLTZ_INSTABILITY
    return main_KHInstability(argc,argv);
#endif //TASK_KELVINHELMHOLTZ_INSTABILITY

#ifdef TASK_TEST_FUVFEEDBACK
    return main_TestFUVFeedback(argc,argv);
#endif //TASK_TEST_FUVFEEDBACK



#ifdef TASK_GALACTIC_CENTER
    return main_GalacticCenter(argc,argv);
#endif //TASK_GALACTIC_CENTER

#ifdef TASK_SANTABARBARA
    return main_SantaBarbara(argc,argv);
#endif //TASK_SANTABARBARA

#ifdef TASK_1D_TWOFLUIDS //{
    return main_TwoFluids();
#endif // TASK_1D_TWOFLUIDS //}

#ifdef TASK_KEPLER
    return main_Kepler();
#endif //TASK_KEPLER

#ifdef TASK_GALAXY_FORMATION
    return main_GalaxyFormation(argc,argv);
#endif //TASK_GALAXY_FORMATION

#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
    //return main_TestSymmetrizedPotentialError(argc,argv);
    return main_TestSymmetrizedPotentialError2(argc,argv);
    //return main_TestSymmetrizedPotentialError_ColdCollapse(argc,argv);
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR

#ifdef TASK_TEST_NEIGHBORSEARCH
    return main_Test_NeighborSearch(argc,argv);
#endif

#ifdef TASK_TEST_HYDRO_QUANTITIES //{
    return main_Test_HydroQuantities(argc,argv);
#endif //TASK_TEST_HYDRO_QUANTITIES //}

#ifdef TASK_TEST_EQUILIBRIUM_TEMPERATURE //{
    return main_Test_EquilibriumTemperature(argc,argv);
#endif //TASK_TEST_EQUILIBRIUM_TEMPERATURE //}

#ifdef TASK_TEST_STELLARFEEDBACK // TASK_TEST_STELLARFEEDBACK //{
    return main_StellarFeedbackTest(argc,argv);
#endif // TASK_TEST_STELLARFEEDBACK //}

#ifdef TASK_TEST_MOMENTUMFEEDBACK // TASK_TEST_MOMENTUMFEEDBACK //{
    return main_MomentumFeedbackTest(argc,argv);
#endif // TASK_TEST_MOMENTUMFEEDBACK //}

#ifdef TASK_TEST_1D_THERMAL_CONDUCTIVITY //{
    return main_Test_1DThermalConductivity(argc,argv);
#endif // TASK_TEST_1D_THERMAL_CONDUCTIVITY //{

#ifdef TASK_TEST_3D_THERMAL_CONDUCTIVITY //{
    return main_Test_3DThermalConductivity(argc,argv);
#endif // TASK_TEST_1D_THERMAL_CONDUCTIVITY //{

#ifdef TASK_TEST_GRAVITYTREE
    return main_GravityTreeTest(argc,argv);
#endif
#ifdef TASK_TEST_HYDROTREE
    //return main_HydroTreeTest(argc,argv);
    //return main_HydroTreeIntegralTest(argc,argv);
    //return main_HydroTreeRobustTest(argc,argv);
#endif
#ifdef TASK_TEST_SINKPARTICLE
    //return main_Test_Sink(argc,argv);
    return main_Test_SinkParticleRun(argc,argv);
#endif


#ifdef TASK_TEST_DOMAINDECOMPOSITION //{
    return main_TestDomainDecomposition(argc,argv);
#endif // TASK_TEST_DOMAINDECOMPOSITION //}

#ifdef TASK_TEST_PARTICLE_SPLITTING //{
    return main_Test_ParticleSplitting(argc,argv);
#endif // TASK_TEST_PARTICLE_SPLITTING //}

    //return main_NeighborSearchTest2();
    // return main_HydroCheckDirectAndTree(argc,argv);

    //return main_IOTest(argc,argv);
    //return main_CosmologicalRun(argc,argv);
    //return main_IsothermalSphericalCollapse();
    //return main_M2SphericalCollapse();
    //return main_BreakingDam();
    //return main_NavarroWhiteTest();
    //return main_WadaNorman(argc,argv);

    //return main_SelfSimilarCooling();
    //return main_3DShockTube();
    //return main_Orbit();
    //return main_BenchMarkTest();
    //return main_BenchMarkTestForForce();
    //return main_FOFTest(argc,argv);


    //return main_ColdCollapseTest_ActiveParticles();
    //return main_NeighborSearchTest();
    //return main_ForceHydroTest();
    //return main_TreeTest();
    //return main_ReadAndPlantTreeTest(argc,argv);

    FinalizeMPIEnv();

    return (EXIT_SUCCESS);
}

/////////////////////////////////////
//                                 //
// Main routines for product runs. //
//                                 //
/////////////////////////////////////


#ifdef TASK_COLD_COLLAPSE
int main_ColdCollapseTest(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    //InitColdCollapseTest(1024,0.1,ON);
    InitColdCollapseTestMixed(1024,0.1,1.0/1.0,ON);
   
    fprintf(stderr,"%g %g\n",Pbody[0]->Mass,Pbody[0]->Eps);
    fprintf(stderr,"%g %g\n",Pbody[1]->Mass,Pbody[1]->Eps);

    //
    InitColdCollapseTestMixed(1024,0.1,1.0/8.0,ON);
   
    fprintf(stderr,"%g %g\n",Pbody[0]->Mass,Pbody[0]->Eps);
    fprintf(stderr,"%g %g\n",Pbody[1]->Mass,Pbody[1]->Eps);
    InitColdCollapseTestMixed(1024,0.1,1.0/64.0,ON);
   
    fprintf(stderr,"%g %g\n",Pbody[0]->Mass,Pbody[0]->Eps);
    fprintf(stderr,"%g %g\n",Pbody[1]->Mass,Pbody[1]->Eps);
    exit(1);

    if (Pall.RunStatus == RestartSimulation){
        fprintf(stderr,"This task does not support the restart run!\n");
        exit(EXIT_FAILURE);
    }

    InitLogFiles();
    InitializeRun();


#if 0
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();
#include <gp5util.h>

    int Npipes = 4;
    double FieldPos[Pall.Ntotal_t][3],FieldMass[Pall.Ntotal_t],FieldEps2[Pall.Ntotal_t],FieldEps[Pall.Ntotal_t];
    for(int i=0;i<Pall.Ntotal;i++){
        FieldPos[i][0] = Pbody[i]->Pos[0];
        FieldPos[i][1] = Pbody[i]->Pos[1];
        FieldPos[i][2] = Pbody[i]->Pos[2];
        FieldMass[i] = Pbody[i]->Mass;
        FieldEps2[i] = SQ(Pbody[i]->Eps);
        FieldEps[i] = Pbody[i]->Eps;
    }
    g5_set_n(Pall.Ntotal_t);
    g5_set_xmeps2j(0,Pall.Ntotal,FieldPos,FieldMass,FieldEps2);

    double adummy[Npipes][3],phidummy[Npipes];
    double time = GetElapsedTime();
    for(int k=0;k<100;k++){
    for(int i=0;i<Pall.Ntotal;i+=Npipes){
        g5_set_xi(Npipes,FieldPos+i);
        g5_set_eps(Npipes,FieldEps+i);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
    }
    }
    fprintf(stderr,"%g \n",GetElapsedTime()-time);

    time = GetElapsedTime();
    for(int k=0;k<100;k++){
    for(int i=0;i<Pall.Ntotal;i+=Npipes){
        g5_set_xi(Npipes,FieldPos+i);
        g5_set_eps(Npipes,FieldEps+i);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
    }
    }
    fprintf(stderr,"%g \n",GetElapsedTime()-time);

    fflush(NULL);

    exit(1);
#endif
    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    CalcSymmetrizedPotential();
    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif

#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
static int Nstart,Nend;
static double Rsphere,Distance;
static double ThetaStart,ThetaEnd,ThetaStep;
int MassRatioStart,MassRatioEnd,MassRatioStep;
static double Ndeg;
static int LogTheta;
static void ReadParametersForSymmetrizedPotentialError(void){

    FILE *fp;
    FileOpen(fp,"./PotentialTest.param","r");
    fscanf(fp,"%d",&Nstart);
    fscanf(fp,"%d",&Nend);
    fprintf(stderr,"Nstart and Nend = %d %d\n",Nstart,Nend);

    fscanf(fp,"%le",&Rsphere);
    fscanf(fp,"%le",&Distance);
    fprintf(stderr,"Rsphere and Distance = %g %g\n",Rsphere,Distance);

    fscanf(fp,"%d",&LogTheta);
    fprintf(stderr,"LogTheta = %d\n",LogTheta);

    fscanf(fp,"%le",&ThetaStart);
    fscanf(fp,"%le",&ThetaEnd);
    fscanf(fp,"%le",&ThetaStep);
    fprintf(stderr,"ThetaStep, ThetaEnd, and ThetaStep = %g %g %g\n",ThetaStart,ThetaEnd,ThetaStep);

    fscanf(fp,"%d",&MassRatioStart);
    fscanf(fp,"%d",&MassRatioEnd);
    fscanf(fp,"%d",&MassRatioStep);
    fprintf(stderr,"MassRatioStart, MassRatioEnd, MassRatioStep  = %d %d %d\n",MassRatioStart,MassRatioEnd,MassRatioStep);

    fscanf(fp,"%le",&Ndeg);
    fprintf(stderr,"Ndeg = %g\n",Ndeg);

    fclose(fp);

    return;
}

static double ReturnPotForSymmetrizedPotentailError(void){

    double Pos[] = {Pbody[Pall.Ntotal-1]->Pos[0],Pbody[Pall.Ntotal-1]->Pos[1],Pbody[Pall.Ntotal-1]->Pos[2]};
    double Eps2 = SQ(Pbody[Pall.Ntotal-1]->Eps);
    double pot = 0.e0;
    for(int i=0;i<Pall.Ntotal-1;i++){
        double Distance = sqrt(DISTANCE2(Pos,Pbody[i]->Pos)+Eps2+SQ(Pbody[i]->Eps));
        pot += Pbody[i]->Mass/Distance;
        //fprintf(stderr,"[T] P %g M %g D %g \n",pot,Pbody[i]->Mass,Distance);
    }

    return -Pall.GravConst*pot;
}

static void CalcPotAccForSymmetrizedPotentialErrorEstimation(void){

    for(int i=0;i<Pall.Ntotal;i++){
        double Pos[] = {Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2]};
        double Eps2 = SQ(Pbody[i]->Eps);
        
        Pbody[i]->Pot = Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
        for(int k=0;k<Pall.Ntotal-1;k++){
            if(k == i) continue;
            double dx[] = {Pos[0]-Pbody[k]->Pos[0],Pos[1]-Pbody[k]->Pos[1],Pos[2]-Pbody[k]->Pos[2]};

            double Distance2 = NORM2(dx)+Eps2+SQ(Pbody[k]->Eps);
            double Distance = sqrt(Distance2);
            Pbody[i]->Pot += -Pall.GravConst*Pbody[k]->Mass/Distance;
            Pbody[i]->Acc[0] += -Pall.GravConst*Pbody[k]->Mass/Distance2 * (dx[0]/Distance);
            Pbody[i]->Acc[1] += -Pall.GravConst*Pbody[k]->Mass/Distance2 * (dx[1]/Distance);
            Pbody[i]->Acc[2] += -Pall.GravConst*Pbody[k]->Mass/Distance2 * (dx[2]/Distance);
        }
        Pbody[i]->Pot *= 0.5*Pbody[i]->Mass;
    }
    
    return ;
}

static void MakeSamplingPoints(double PosOriginal[][3], const bool Random){

    if(Random == true){
        int counter = 0;
        while(counter<Ndeg){
            double Pos[] ={2.0*gsl_rng_uniform(RandomGenerator)-1.0,
                           2.0*gsl_rng_uniform(RandomGenerator)-1.0,
                           2.0*gsl_rng_uniform(RandomGenerator)-1.0};
            if(NORM(Pos)<1.0){
                PosOriginal[counter][0] = Rsphere*Pos[0];
                PosOriginal[counter][1] = Rsphere*Pos[1];
                PosOriginal[counter][2] = Rsphere*Pos[2];
                counter ++;
            }
        } 
    }else{
        double dT = 2*M_PI/Ndeg;
        double Pos[] = {Distance,0.e0};
        for(int i=0;i<Ndeg;i++){
            double T = i*dT;
            PosOriginal[i][0] = cos(T)*Pos[0]-sin(T)*Pos[1];
            PosOriginal[i][1] = sin(T)*Pos[0]+cos(T)*Pos[1];
            PosOriginal[i][2] = 0.e0;
        }
    }

    return ;
}

int main_TestSymmetrizedPotentialError(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    ReadParametersForSymmetrizedPotentialError();
    InitLogFiles();

    int NTheta = 0; 
    for(double Theta=ThetaStart;Theta<ThetaEnd;Theta+=ThetaStep){
        NTheta ++;
    }   

    InitSphericalParticleDistributionWithTwoTypes(Nstart,1.0,Rsphere,Distance);

    double (*PosOriginal)[3];
    double *Mass;
    double *Eps;
    PosOriginal = malloc(sizeof(double)*Ndeg*3);
    Eps = malloc(sizeof(double)*Ndeg);
    MakeSamplingPoints(PosOriginal,false);

    InitializeDecomposition();
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();
    FreeStructPbody();

    fprintf(stderr,"Nstart and Nend = %d %d\n",Nstart,Nend);
    for(int l=Nstart;l<Nend;l=2*l){
        for(int j=MassRatioStart;j<=MassRatioEnd;j*=MassRatioStep){
            // Make Initial Condition.
            InitSphericalParticleDistributionWithTwoTypes(l,1.0/(double)j,Rsphere,Distance);
            fprintf(stderr,"## %g %g\n",Pbody[Pall.Ntotal-1]->Pos[0],Pbody[Pall.Ntotal-1]->Pos[1]);
            for(int p=0;p<Ndeg;p++)
                Eps[p] = Pbody[Pall.Ntotal_t-1]->Eps;

            //Kick1Drift();
            BuildPredictors(); // Pos -> PosP/ Vel -> VelP
            //DomainDecomposition();
            
            double dT = 2*M_PI/Ndeg;
            
            // Calc true potential.
            double PotLog_Direct[NTheta];
            int counter_d = 0;
            for(double Theta=ThetaStart;Theta<ThetaEnd;Theta+=ThetaStep){ //
                ClearGravitationalForce();
                //PotLog_Direct[counter_d] = ReturnPotForSymmetrizedPotentailError();
                double Pos[] = {Pbody[Pall.Ntotal-1]->Pos[0],Pbody[Pall.Ntotal-1]->Pos[1],Pbody[Pall.Ntotal-1]->Pos[2]};
                double PosP[] = {Pbody[Pall.Ntotal-1]->PosP[0],Pbody[Pall.Ntotal-1]->PosP[1],Pbody[Pall.Ntotal-1]->PosP[2]};
                PotLog_Direct[counter_d] = 0.e0;
#if 0

                for(int p=0;p<Ndeg;p++){
                    double T = p*dT;
                    Pbody[Pall.Ntotal-1]->Pos[0] = cos(T)*Pos[0]-sin(T)*Pos[1];
                    Pbody[Pall.Ntotal-1]->Pos[1] = sin(T)*Pos[0]+cos(T)*Pos[1];
                    
                    Pbody[Pall.Ntotal-1]->PosP[0] = cos(T)*PosP[0]-sin(T)*PosP[1];
                    Pbody[Pall.Ntotal-1]->PosP[1] = sin(T)*PosP[0]+cos(T)*PosP[1];
                    PotLog_Direct[counter_d] += ReturnPotForSymmetrizedPotentailError();
                    fprintf(stderr,"%g %g | %g %g || %g %g\n",
                        Pbody[Pall.Ntotal-1]->Pos[0],Pbody[Pall.Ntotal-1]->Pos[1],
                        PosOriginal[p][0],PosOriginal[p][1],Pos[0],Pos[1]);
                }
#else
                for(int p=0;p<Ndeg;p++){
                    Pbody[Pall.Ntotal-1]->Pos[0] = PosOriginal[p][0];
                    Pbody[Pall.Ntotal-1]->Pos[1] = PosOriginal[p][1];
                    Pbody[Pall.Ntotal-1]->PosP[0] = PosOriginal[p][0];
                    Pbody[Pall.Ntotal-1]->PosP[1] = PosOriginal[p][1];
                    Pbody[Pall.Ntotal-1]->Eps = Eps[p];
                    PotLog_Direct[counter_d] += ReturnPotForSymmetrizedPotentailError();
                    fprintf(stderr,"POTENTIAL TRUE %1.15g, (%1.15g,%1.15g;%1.15g)\n",ReturnPotForSymmetrizedPotentailError(),
                            Pbody[Pall.Ntotal-1]->PosP[0],Pbody[Pall.Ntotal-1]->PosP[1],Pbody[Pall.Ntotal-1]->Eps = Eps[p]);
                }
#endif
                Pbody[Pall.Ntotal-1]->PosP[0] = Pbody[Pall.Ntotal-1]->Pos[0] = Pos[0];
                Pbody[Pall.Ntotal-1]->PosP[1] = Pbody[Pall.Ntotal-1]->Pos[1] = Pos[1];

                PotLog_Direct[counter_d] /= Ndeg;

                counter_d ++;
            }
            

            // First force / hydro calculation.
            PlantGravityTree();

            fprintf(stderr,"Number of particles = %d, N = %d, Mass ratio = %g\n",
                    Pall.Ntotal_t,l,1.0/(double)j);

            double PotLog_Sym[NTheta];
            int counter = 0;
            for(double Theta=ThetaStart;Theta<ThetaEnd;Theta+=ThetaStep){ //

                GravityRoot.OpeningAngle = Theta;

                // Calc force from Two trees.
                ClearGravitationalForce();
                ////  ForceParallelTreeGRAPE();

                //fprintf(stderr,"Potential[Theta = %1.2g] = %1.15g\n",Theta,Pbody[Pall.Ntotal_t-1]->Pot);

                ////  PotLog_Sym[counter] = Pbody[Pall.Ntotal_t-1]->Pot;
                //PotLog_Sym[counter] = ReturnPotForSymmetrizedPotentailError();
                
                double Pos[] = {Pbody[Pall.Ntotal-1]->Pos[0],Pbody[Pall.Ntotal-1]->Pos[1],Pbody[Pall.Ntotal-1]->Pos[2]};
                double PosP[] = {Pbody[Pall.Ntotal-1]->PosP[0],Pbody[Pall.Ntotal-1]->PosP[1],Pbody[Pall.Ntotal-1]->PosP[2]};
                PotLog_Sym[counter] = 0.e0;
#if 0
                for(int p=0;p<Ndeg;p++){
                    double T = p*dT;
                    Pbody[Pall.Ntotal-1]->Pos[0] = cos(T)*Pos[0]-sin(T)*Pos[1];
                    Pbody[Pall.Ntotal-1]->Pos[1] = sin(T)*Pos[0]+cos(T)*Pos[1];
                    
                    Pbody[Pall.Ntotal-1]->PosP[0] = cos(T)*PosP[0]-sin(T)*PosP[1];
                    Pbody[Pall.Ntotal-1]->PosP[1] = sin(T)*PosP[0]+cos(T)*PosP[1];
                    ClearGravitationalForce();
                    ForceParallelTreeGRAPE();
                    PotLog_Sym[counter] += Pbody[Pall.Ntotal_t-1]->Pot;
                }
#else
                for(int p=0;p<Ndeg;p++){
                    Pbody[Pall.Ntotal-1]->Pos[0] = PosOriginal[p][0];
                    Pbody[Pall.Ntotal-1]->Pos[1] = PosOriginal[p][1];
                    Pbody[Pall.Ntotal-1]->PosP[0] = PosOriginal[p][0];
                    Pbody[Pall.Ntotal-1]->PosP[1] = PosOriginal[p][1];
                    Pbody[Pall.Ntotal-1]->Eps = Eps[p];
                    ClearGravitationalForce();
                    PlantGravityTree();
                    ForceParallelTreeGRAPE();
                    PotLog_Sym[counter] += Pbody[Pall.Ntotal_t-1]->Pot;
                    fprintf(stderr,"POTENTIAL TG %1.15g, (%1.15g,%1.15g;%1.15g)\n",Pbody[Pall.Ntotal_t-1]->Pot,
                            Pbody[Pall.Ntotal-1]->PosP[0],Pbody[Pall.Ntotal-1]->PosP[1],Pbody[Pall.Ntotal-1]->Eps = Eps[p]);
                }
#endif
                Pbody[Pall.Ntotal-1]->PosP[0] = Pbody[Pall.Ntotal-1]->Pos[0] = Pos[0];
                Pbody[Pall.Ntotal-1]->PosP[1] = Pbody[Pall.Ntotal-1]->Pos[1] = Pos[1];

                PotLog_Sym[counter] /= Ndeg;
                counter ++;
            }

            {
                FILE *fp;
                FileOpen(fp,"./all","w");
                for(int k=0;k<Pall.Ntotal_t;k++){
                    fprintf(fp,"%g %g %g %g %g\n",
                        Pbody[k]->Pos[0],Pbody[k]->Pos[1],Pbody[k]->Pos[2],Pbody[k]->Mass,Pbody[k]->Eps);
                }
                fclose(fp);
            }


            int Nparticles_Sym = l>>1;
            dprintlmpi(Nparticles_Sym);
            double (*Pos_1_Sym)[3];
            double *Mass_1_Sym;
            double *Eps_1_Sym;
            double (*Pos_2_Sym)[3];
            double *Mass_2_Sym;
            double *Eps_2_Sym;

#if 1
            Pos_1_Sym = malloc(sizeof(double)*3*Nparticles_Sym);
            Pos_2_Sym = malloc(sizeof(double)*3*Nparticles_Sym);
            Mass_1_Sym = malloc(sizeof(double)*Nparticles_Sym);
            Mass_2_Sym = malloc(sizeof(double)*Nparticles_Sym);
            Eps_1_Sym = malloc(sizeof(double)*Nparticles_Sym);
            Eps_2_Sym = malloc(sizeof(double)*Nparticles_Sym);
#endif

#if 1
            for(int k=0;k<l/2;k++){
                int index = 2*k;
                Pos_1_Sym[k][0] = Pbody[index]->Pos[0];
                Pos_1_Sym[k][1] = Pbody[index]->Pos[1];
                Pos_1_Sym[k][2] = Pbody[index]->Pos[2];
                Mass_1_Sym[k]   = Pbody[index]->Mass;
                Eps_1_Sym[k]    = Pbody[index]->Eps;
            }
            for(int k=0;k<l/2;k++){
                int index = 2*k+1;
                Pos_2_Sym[k][0] = Pbody[index]->Pos[0];
                Pos_2_Sym[k][1] = Pbody[index]->Pos[1];
                Pos_2_Sym[k][2] = Pbody[index]->Pos[2];
                Mass_2_Sym[k]   = Pbody[index]->Mass;
                Eps_2_Sym[k]    = Pbody[index]->Eps;
            }
#endif

            double Mass_p = Pbody[Pall.Ntotal_t-1]->Mass;
            double Eps_p = Pbody[Pall.Ntotal_t-1]->Eps;

            FreeStructPbody();

            //////////////// Copy Particle 1
            RestoreSphericalParticleDistribution(Nparticles_Sym,Pos_1_Sym,Mass_1_Sym,Eps_1_Sym,Mass_p,Eps_p,Rsphere,Distance);
            BuildPredictors(); 
            PlantGravityTree();
            fprintf(stderr,"Pall.Ntotal_t = %ld\n",Pall.Ntotal_t);
            double PotLog_flat[NTheta];
            counter = 0;
            for(double Theta=ThetaStart;Theta<ThetaEnd;Theta+=ThetaStep){ 
                GravityRoot.OpeningAngle = Theta;
                ///  ClearGravitationalForce();
                ///  ForceParallelTreeGRAPE();
                //fprintf(stderr,"Potential[Theta = %1.2g] = %1.15g\n",Theta,Pbody[Pall.Ntotal_t-1]->Pot);
                ///  PotLog_flat[counter] = Pbody[Pall.Ntotal_t-1]->Pot;
                //PotLog_flat[counter] = ReturnPotForSymmetrizedPotentailError();
                
                double Pos[] = {Pbody[Pall.Ntotal-1]->Pos[0],Pbody[Pall.Ntotal-1]->Pos[1],Pbody[Pall.Ntotal-1]->Pos[2]};
                double PosP[] = {Pbody[Pall.Ntotal-1]->PosP[0],Pbody[Pall.Ntotal-1]->PosP[1],Pbody[Pall.Ntotal-1]->PosP[2]};
                PotLog_flat[counter] = 0.e0;
                for(int p=0;p<Ndeg;p++){
                    double T = p*dT;
                    Pbody[Pall.Ntotal-1]->Pos[0] = cos(T)*Pos[0]-sin(T)*Pos[1];
                    Pbody[Pall.Ntotal-1]->Pos[1] = sin(T)*Pos[0]+cos(T)*Pos[1];
                    
                    Pbody[Pall.Ntotal-1]->PosP[0] = cos(T)*PosP[0]-sin(T)*PosP[1];
                    Pbody[Pall.Ntotal-1]->PosP[1] = sin(T)*PosP[0]+cos(T)*PosP[1];
                    ClearGravitationalForce();
                    PlantGravityTree();
                    ForceParallelTreeGRAPE();
                    PotLog_flat[counter] += Pbody[Pall.Ntotal_t-1]->Pot;
                }
                Pbody[Pall.Ntotal-1]->PosP[0] = Pbody[Pall.Ntotal-1]->Pos[0] = Pos[0];
                Pbody[Pall.Ntotal-1]->PosP[1] = Pbody[Pall.Ntotal-1]->Pos[1] = Pos[1];

                //PotLog_flat[counter] /= Ndeg;
                counter ++;
            }
            {
                FILE *fp;
                FileOpen(fp,"./p1","w");
                for(int k=0;k<Pall.Ntotal_t;k++){
                    fprintf(fp,"%g %g %g %g %g\n",
                        Pbody[k]->Pos[0],Pbody[k]->Pos[1],Pbody[k]->Pos[2],Pbody[k]->Mass,Pbody[k]->Eps);
                }
                fclose(fp);
            }
            FreeStructPbody();

            //////////////// Copy Particle 2
            RestoreSphericalParticleDistribution(Nparticles_Sym,Pos_2_Sym,Mass_2_Sym,Eps_2_Sym,Mass_p,Eps_p,Rsphere,Distance);
            BuildPredictors(); 
            PlantGravityTree();
            fprintf(stderr,"Pall.Ntotal_t = %ld\n",Pall.Ntotal_t);

            counter = 0;
            for(double Theta=ThetaStart;Theta<ThetaEnd;Theta+=ThetaStep){ 
                GravityRoot.OpeningAngle = Theta;
                //  ClearGravitationalForce();
                //  ForceParallelTreeGRAPE();
                //fprintf(stderr,"Potential[Theta = %1.2g] = %1.15g\n",Theta,Pbody[Pall.Ntotal_t-1]->Pot);
                //  PotLog_flat[counter] += Pbody[Pall.Ntotal_t-1]->Pot;
                //PotLog_flat[counter] += ReturnPotForSymmetrizedPotentailError();
                
                double Pos[] = {Pbody[Pall.Ntotal-1]->Pos[0],Pbody[Pall.Ntotal-1]->Pos[1],Pbody[Pall.Ntotal-1]->Pos[2]};
                double PosP[] = {Pbody[Pall.Ntotal-1]->PosP[0],Pbody[Pall.Ntotal-1]->PosP[1],Pbody[Pall.Ntotal-1]->PosP[2]};
                for(int p=0;p<Ndeg;p++){
                    double T = p*dT;
                    Pbody[Pall.Ntotal-1]->Pos[0] = cos(T)*Pos[0]-sin(T)*Pos[1];
                    Pbody[Pall.Ntotal-1]->Pos[1] = sin(T)*Pos[0]+cos(T)*Pos[1];
                    
                    Pbody[Pall.Ntotal-1]->PosP[0] = cos(T)*PosP[0]-sin(T)*PosP[1];
                    Pbody[Pall.Ntotal-1]->PosP[1] = sin(T)*PosP[0]+cos(T)*PosP[1];
                    ClearGravitationalForce();
                    PlantGravityTree();
                    ForceParallelTreeGRAPE();
                    PotLog_flat[counter] += Pbody[Pall.Ntotal_t-1]->Pot;
                    fprintf(stderr,"POTENTIAL TGb %1.15g, (%1.15g,%1.15g;%1.15g)\n",Pbody[Pall.Ntotal_t-1]->Pot,
                            Pbody[Pall.Ntotal-1]->PosP[0],Pbody[Pall.Ntotal-1]->PosP[1],Pbody[Pall.Ntotal-1]->Eps = Eps[p]);
                }
                Pbody[Pall.Ntotal-1]->PosP[0] = Pbody[Pall.Ntotal-1]->Pos[0] = Pos[0];
                Pbody[Pall.Ntotal-1]->PosP[1] = Pbody[Pall.Ntotal-1]->Pos[1] = Pos[1];

                PotLog_flat[counter] /= Ndeg;
                counter ++;
            }
            {
                FILE *fp;
                FileOpen(fp,"./p2","w");
                for(int k=0;k<Pall.Ntotal_t;k++){
                    fprintf(fp,"%g %g %g %g %g\n",
                        Pbody[k]->Pos[0],Pbody[k]->Pos[1],Pbody[k]->Pos[2],Pbody[k]->Mass,Pbody[k]->Eps);
                }
                fclose(fp);
            }
            FreeStructPbody();

            counter = 0;
            for(double Theta=ThetaStart;Theta<ThetaEnd;Theta+=ThetaStep){ 
                fprintf(stderr,"|(Pot_sym - Pot_Dir)|, |(Pot_flat - Pot_Dir)| (Theta = %g) = %1.15g %1.15g %1.15g %1.15g %1.15g\n",Theta,
                    //fabs(PotLog_flat[counter]-PotLog_Sym[counter]),
#if 0
                    fabs((SQ(PotLog_Sym[counter])-SQ(PotLog_Direct[counter]))/(PotLog_Sym[counter]+PotLog_Direct[counter])),
                    fabs((SQ(PotLog_flat[counter])-SQ(PotLog_Direct[counter]))/(PotLog_flat[counter]+PotLog_Direct[counter])),
#else
                    fabs(PotLog_Sym[counter]-PotLog_Direct[counter]),fabs(PotLog_flat[counter]-PotLog_Direct[counter]),
#endif
                    PotLog_flat[counter],PotLog_Sym[counter],
                    PotLog_Direct[counter]);
                //fprintf(stderr,"|(Pot_sym - Pot_t)/Pot_t| (Theta = %g) = %1.15g, %g %g\n",Theta,
                    //fabs(PotLog_flat[counter]-PotLog_Sym[counter])/PotLog_Sym[counter],
                    //PotLog_flat[counter],PotLog_Sym[counter]);
                counter ++;
            }

#if 1
            free(Pos_1_Sym);
            free(Pos_2_Sym);
            free(Mass_1_Sym);
            free(Mass_2_Sym);
            free(Eps_1_Sym);
            free(Eps_2_Sym);
#endif
        }
    }

    CloseLogFiles();

#if 0

    //InitColdCollapseTest(1024,0.1,ON);
    InitColdCollapseTestMixed(1024,0.1,1.0/1.0,ON);
   
    fprintf(stderr,"%g %g\n",Pbody[0]->Mass,Pbody[0]->Eps);
    fprintf(stderr,"%g %g\n",Pbody[1]->Mass,Pbody[1]->Eps);

    //
    InitColdCollapseTestMixed(1024,0.1,1.0/8.0,ON);
   
    fprintf(stderr,"%g %g\n",Pbody[0]->Mass,Pbody[0]->Eps);
    fprintf(stderr,"%g %g\n",Pbody[1]->Mass,Pbody[1]->Eps);
    InitColdCollapseTestMixed(1024,0.1,1.0/64.0,ON);
   
    fprintf(stderr,"%g %g\n",Pbody[0]->Mass,Pbody[0]->Eps);
    fprintf(stderr,"%g %g\n",Pbody[1]->Mass,Pbody[1]->Eps);
    exit(1);

    if (Pall.RunStatus == RestartSimulation){
        fprintf(stderr,"This task does not support the restart run!\n");
        exit(EXIT_FAILURE);
    }

    InitLogFiles();
    InitializeRun();

#endif

#if 0
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();
#include <gp5util.h>

    int Npipes = 4;
    double FieldPos[Pall.Ntotal_t][3],FieldMass[Pall.Ntotal_t],FieldEps2[Pall.Ntotal_t],FieldEps[Pall.Ntotal_t];
    for(int i=0;i<Pall.Ntotal;i++){
        FieldPos[i][0] = Pbody[i]->Pos[0];
        FieldPos[i][1] = Pbody[i]->Pos[1];
        FieldPos[i][2] = Pbody[i]->Pos[2];
        FieldMass[i] = Pbody[i]->Mass;
        FieldEps2[i] = SQ(Pbody[i]->Eps);
        FieldEps[i] = Pbody[i]->Eps;
    }
    g5_set_n(Pall.Ntotal_t);
    g5_set_xmeps2j(0,Pall.Ntotal,FieldPos,FieldMass,FieldEps2);

    double adummy[Npipes][3],phidummy[Npipes];
    double time = GetElapsedTime();
    for(int k=0;k<100;k++){
    for(int i=0;i<Pall.Ntotal;i+=Npipes){
        g5_set_xi(Npipes,FieldPos+i);
        g5_set_eps(Npipes,FieldEps+i);
        g5_run_symmetrized_softening();
        g5_get_force(Npipes,adummy,phidummy);
    }
    }
    fprintf(stderr,"%g \n",GetElapsedTime()-time);

    time = GetElapsedTime();
    for(int k=0;k<100;k++){
    for(int i=0;i<Pall.Ntotal;i+=Npipes){
        g5_set_xi(Npipes,FieldPos+i);
        g5_set_eps(Npipes,FieldEps+i);
        g5_run();
        g5_get_force(Npipes,adummy,phidummy);
    }
    }
    fprintf(stderr,"%g \n",GetElapsedTime()-time);

    fflush(NULL);

    exit(1);
#endif
    //Run();

    //if(MPIGetMyID() == MPI_ROOT_RANK)
        //fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    //CalcSymmetrizedPotential();
    //LogOutPutEnergyMomentumAngularMomentum();
    //CloseLogFiles();

    return EXIT_SUCCESS;
}


long int NumberofPCInteractions;
long int NumberofPPInteractions;

int main_TestSymmetrizedPotentialError2(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    ReadParametersForSymmetrizedPotentialError();
    InitLogFiles();

    int NTheta = 0; 
    if(LogTheta == 1){
        for(double Theta=ThetaEnd;Theta>ThetaStart;Theta*=ThetaStep){
            NTheta ++;
        }
        NTheta ++;
    } else {
        for(double Theta=ThetaStart;Theta<ThetaEnd;Theta+=ThetaStep){
            NTheta ++;
        }
    }

    double ThetaArray[NTheta];
    if(LogTheta == 1){
        int counter = 0;
        for(double Theta=ThetaEnd;Theta>ThetaStart;Theta*=ThetaStep){
            ThetaArray[NTheta-1-counter] = Theta;
            gprintlmpi(ThetaArray[NTheta-1-counter]);
            counter ++;
        }
        ThetaArray[0] = 0.e0;
        gprintlmpi(ThetaArray[0]);
    } else {
        for(int i=0;i<NTheta;i++){
            ThetaArray[i] = ThetaStart+ThetaStep*i;
            gprintlmpi(ThetaArray[i]);
        }
    }
    //exit(1);

    double *Pot_true[NTheta];
    double *Pot_TG[NTheta];
    double *Pot_TGb1[NTheta];
    double *Pot_TGb2[NTheta];
    for(int i=0;i<NTheta;i++){
        Pot_true[i] = malloc(sizeof(double)*Nstart);
        Pot_TG[i] = malloc(sizeof(double)*Nstart);
        Pot_TGb1[i] = malloc(sizeof(double)*Nstart);
        Pot_TGb2[i] = malloc(sizeof(double)*Nstart);
    }
    double (*Acc_true[NTheta])[3];
    double (*Acc_TG[NTheta])[3];
    double (*Acc_TGb1[NTheta])[3];
    double (*Acc_TGb2[NTheta])[3];
    for(int i=0;i<NTheta;i++){
        Acc_true[i] = malloc(sizeof(double)*3*Nstart);
        Acc_TG[i] = malloc(sizeof(double)*3*Nstart);
        Acc_TGb1[i] = malloc(sizeof(double)*3*Nstart);
        Acc_TGb2[i] = malloc(sizeof(double)*3*Nstart);
    }
    double PPInteractions[NTheta];
    double PCInteractions[NTheta];
    double PPInteractionsSep[NTheta];
    double PCInteractionsSep[NTheta];

    FILE *fp_log;
    FileOpen(fp_log,"./Accuracy.dat","w");

    FILE *fp_melog;
    FileOpen(fp_melog,"./MassEps.dat","w");


    static bool first_error_estimation = true;
    for(int j=MassRatioStart;j<=MassRatioEnd;j*=MassRatioStep){
        dprintlmpi(j);
        InitSphericalParticleDistributionForErrorEstimation(Nstart,1.0/(double)j,Rsphere);

        fprintf(fp_melog,"%d %d %1.15g %1.15g %1.15g %1.15g\n",1,j,
                Pbody[0]->Mass,Pbody[0]->Eps,Pbody[1]->Mass,Pbody[1]->Eps);


        for(int i=0;i<Pall.Ntotal;i++)
            //Pbody[i]->Eps *= 10;

        if(first_error_estimation == true){
            InitializeDecomposition();
            InitializeRootForGravity();
            InitializeRootForLET();
            InitializeParallelTreeGRAPE();
            first_error_estimation = false;
        }

        BuildPredictors();

        PlantGravityTree();
        /// Are there necessary routines?
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
        ClearGravitationalForce();
        //CalcPotAccForSymmetrizedPotentialErrorEstimation();
#if 0
        HoldGRAPE();
        CalcGravityDirectSymmetrizedPotentialWithPhantom();
        ReleaseGRAPE();

        for(int k=0;k<NTheta;k++){
            for(int i=0;i<Pall.Ntotal;i++){
                Pot_true[k][i] = Pbody[i]->Pot; 
                Acc_true[k][i][0] = Pbody[i]->Acc[0]; 
                Acc_true[k][i][1] = Pbody[i]->Acc[1]; 
                Acc_true[k][i][2] = Pbody[i]->Acc[2]; 
            }
        }
#endif

        //PlantGravityTree();
        //GravityRoot.NumberofLeavesInGroup = TreeNGroup;
        //UpdateNumberofLeavesInGroupInLET(TreeNGroup);
        for(int k=0;k<NTheta;k++){
            GravityRoot.OpeningAngle = ThetaArray[k];
            ClearGravitationalForce();
            NumberofPCInteractions = NumberofPPInteractions = 0;
            ForceParallelTreeGRAPE();
            for(int i=0;i<Pall.Ntotal;i++){
                Pot_TG[k][i] = Pbody[i]->Pot; 
                Acc_TG[k][i][0] = Pbody[i]->Acc[0]; 
                Acc_TG[k][i][1] = Pbody[i]->Acc[1]; 
                Acc_TG[k][i][2] = Pbody[i]->Acc[2]; 
            }
            PPInteractions[k] = NumberofPPInteractions;
            PCInteractions[k] = NumberofPCInteractions;
        }

        //FreeStructPbody();
        //continue;

        int Nparticles_Sym = Nstart>>1;
        double (*Pos_1_Sym)[3];
        double *Mass_1_Sym;
        double *Eps_1_Sym;
        double (*Pos_2_Sym)[3];
        double *Mass_2_Sym;
        double *Eps_2_Sym;

        Pos_1_Sym = malloc(sizeof(double)*3*Nparticles_Sym);
        Pos_2_Sym = malloc(sizeof(double)*3*Nparticles_Sym);
        Mass_1_Sym = malloc(sizeof(double)*Nparticles_Sym);
        Mass_2_Sym = malloc(sizeof(double)*Nparticles_Sym);
        Eps_1_Sym = malloc(sizeof(double)*Nparticles_Sym);
        Eps_2_Sym = malloc(sizeof(double)*Nparticles_Sym);

        for(int k=0;k<Nstart/2;k++){
            int index = 2*k;
            Pos_1_Sym[k][0] = Pbody[index]->Pos[0];
            Pos_1_Sym[k][1] = Pbody[index]->Pos[1];
            Pos_1_Sym[k][2] = Pbody[index]->Pos[2];
            Mass_1_Sym[k]   = Pbody[index]->Mass;
            Eps_1_Sym[k]    = Pbody[index]->Eps;
        }
        for(int k=0;k<Nstart/2;k++){
            int index = 2*k+1;
            Pos_2_Sym[k][0] = Pbody[index]->Pos[0];
            Pos_2_Sym[k][1] = Pbody[index]->Pos[1];
            Pos_2_Sym[k][2] = Pbody[index]->Pos[2];
            Mass_2_Sym[k]   = Pbody[index]->Mass;
            Eps_2_Sym[k]    = Pbody[index]->Eps;
        }

        //////////////// Copy Particle 1
        RestoreSphericalParticleDistributionErrorEstimation(Nparticles_Sym,Pos_1_Sym,Mass_1_Sym,Eps_1_Sym);
        BuildPredictors();
        PlantGravityTree();
        dprintlmpi(GravityRoot.NumberofNodes);

        //GravityRoot.NumberofLeavesInGroup = TreeNGroup/2;
        //UpdateNumberofLeavesInGroupInLET(TreeNGroup/2);
        for(int k=0;k<NTheta;k++){
            GravityRoot.OpeningAngle = ThetaArray[k];
            ClearGravitationalForce();
            NumberofPCInteractions = NumberofPPInteractions = 0;
            ForceParallelTreeGRAPE();
            for(int i=0;i<Nparticles_Sym;i++){
                Pot_TGb1[k][i] = Pbody[i]->Pot; 
                Acc_TGb1[k][i][0] = Pbody[i]->Acc[0]; 
                Acc_TGb1[k][i][1] = Pbody[i]->Acc[1]; 
                Acc_TGb1[k][i][2] = Pbody[i]->Acc[2]; 
            }
            ClearGravitationalForce();
            ForceParallelTreeGRAPEInsert(Nparticles_Sym,Pos_2_Sym,Mass_2_Sym,Eps_2_Sym);
            for(int i=0;i<Nparticles_Sym;i++){
                Pot_TGb1[k][i] += Pbody[i]->Pot; 
                Acc_TGb1[k][i][0] += Pbody[i]->Acc[0]; 
                Acc_TGb1[k][i][1] += Pbody[i]->Acc[1]; 
                Acc_TGb1[k][i][2] += Pbody[i]->Acc[2]; 
            }
            PPInteractionsSep[k] = NumberofPPInteractions;
            PCInteractionsSep[k] = NumberofPCInteractions;
        }
        FreeStructPbody();


        //////////////// Copy Particle 2
        RestoreSphericalParticleDistributionErrorEstimation(Nparticles_Sym,Pos_2_Sym,Mass_2_Sym,Eps_2_Sym);
        BuildPredictors();
        PlantGravityTree();
        for(int k=0;k<NTheta;k++){
            GravityRoot.OpeningAngle = ThetaArray[k];
            ClearGravitationalForce();
            NumberofPCInteractions = NumberofPPInteractions = 0;
            ForceParallelTreeGRAPE();
            for(int i=0;i<Nparticles_Sym;i++){
                Pot_TGb2[k][i] = Pbody[i]->Pot; 
                Acc_TGb2[k][i][0] = Pbody[i]->Acc[0]; 
                Acc_TGb2[k][i][1] = Pbody[i]->Acc[1]; 
                Acc_TGb2[k][i][2] = Pbody[i]->Acc[2]; 
            }
            ClearGravitationalForce();
            ForceParallelTreeGRAPEInsert(Nparticles_Sym,Pos_1_Sym,Mass_1_Sym,Eps_1_Sym);
            for(int i=0;i<Nparticles_Sym;i++){
                Pot_TGb2[k][i] += Pbody[i]->Pot; 
                Acc_TGb2[k][i][0] += Pbody[i]->Acc[0]; 
                Acc_TGb2[k][i][1] += Pbody[i]->Acc[1]; 
                Acc_TGb2[k][i][2] += Pbody[i]->Acc[2]; 
            }
            PPInteractionsSep[k] += NumberofPPInteractions;
            PCInteractionsSep[k] += NumberofPCInteractions;
        }
        FreeStructPbody();


        free(Pos_1_Sym);
        free(Pos_2_Sym);
        free(Mass_1_Sym);
        free(Mass_2_Sym);
        free(Eps_1_Sym);
        free(Eps_2_Sym);

        FILE *fp_log2;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"./Accuracy.%d.%d.dat",j,(int)(10*Rsphere));
        fprintf(stderr,"---- %s\n",Fname);
        FileOpen(fp_log2,Fname,"w");

        FILE *fp_log3;
        Snprintf(Fname,"./Interactions.%d.%d.dat",j,(int)(10*Rsphere));
        fprintf(stderr,"---- %s\n",Fname);
        FileOpen(fp_log3,Fname,"w");

        //fprintf(fp_log,"Mass Ratio = %d:%d\n",1,j);
        for(int k=0;k<NTheta;k++){
            double PotSub = 0.e0;
            double AccSub = 0.e0;
            for(int i=0;i<Nstart;i++){
                PotSub += fabs((Pot_TG[0][i] - Pot_TG[k][i])/Pot_TG[0][i]);
                //PotSub += fabs((Pot_true[k][i] - Pot_TG[k][i])/Pot_true[k][i]);
                AccSub += fabs((
                        fabs(Acc_TG[0][i][0] - Acc_TG[k][i][0])+
                        fabs(Acc_TG[0][i][1] - Acc_TG[k][i][1])+
                        fabs(Acc_TG[0][i][2] - Acc_TG[k][i][2]))
                        /NORM(Acc_TG[0][i]));
            }
            double PotSub_b = 0.e0;
            double AccSub_b = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++){
                //PotSub_b += fabs((Pot_TG[0][i*2] - Pot_TGb1[k][i])/Pot_TG[0][i*2]);
                PotSub_b += fabs((Pot_TGb1[0][i] - Pot_TGb1[k][i])/Pot_TGb1[0][i]);
                AccSub_b += fabs((
                        fabs(Acc_TGb1[0][i][0] - Acc_TGb1[k][i][0])+
                        fabs(Acc_TGb1[0][i][1] - Acc_TGb1[k][i][1])+
                        fabs(Acc_TGb1[0][i][2] - Acc_TGb1[k][i][2]))
                        /NORM(Acc_TGb1[0][i]));
            }
            for(int i=0;i<Pall.Ntotal;i++){
                //PotSub_b += fabs((Pot_TG[0][i*2+1] - Pot_TGb2[k][i])/Pot_TG[0][i*2+1]);
                PotSub_b += fabs((Pot_TGb2[0][i] - Pot_TGb2[k][i])/Pot_TGb2[0][i]);
                AccSub_b += fabs((
                        fabs(Acc_TGb2[0][i][0] - Acc_TGb2[k][i][0])+
                        fabs(Acc_TGb2[0][i][1] - Acc_TGb2[k][i][1])+
                        fabs(Acc_TGb2[0][i][2] - Acc_TGb2[k][i][2]))
                        /NORM(Acc_TGb2[0][i]));
            }
            PotSub /= (double) Nstart;
            PotSub_b /= (double) Nstart;
            //fprintf(stderr,"comp %1.15g %1.15g\n", Pot_true[k][0], Pot_TG[k][0]);

            fprintf(stderr,"|(Pot_sym - Pot_Dir)|, |(Pot_flat - Pot_Dir)| (Theta = %g) = %1.15g %1.15g\n",
                    ThetaArray[k],PotSub,PotSub_b);

            AccSub /= (double) Nstart;
            AccSub_b /= (double) Nstart;

            fprintf(stderr,"|(Acc_sym - Acc_Dir)|, |(Acc_flat - Acc_Dir)| (Theta = %g) = %1.15g %1.15g\n",
                    ThetaArray[k],AccSub,AccSub_b);

            fprintf(fp_log,"%g %g %g %g %g\n",ThetaArray[k],PotSub,PotSub_b,AccSub,AccSub_b);
            fprintf(fp_log2,"%g %g %g %g %g\n",ThetaArray[k],PotSub,PotSub_b,AccSub,AccSub_b);
            fprintf(fp_log3,"%g %g %g %g %g\n",ThetaArray[k],
                    PPInteractions[k]/(double)Nstart,PCInteractions[k]/(double)Nstart,
                    PPInteractionsSep[k]/(double)Nstart,PCInteractionsSep[k]/(double)Nstart);
        }
        fclose(fp_log2);
        fclose(fp_log3);
    }
    fclose(fp_log);
    fclose(fp_melog);


    CloseLogFiles();

    return EXIT_SUCCESS;
}

static void InitializeRunPotentialError(void){

    PlantGravityTree();
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();
    ForceFromExternalPotentials();

    CalcSymmetrizedPotential();

    FirstTimeStep();
    UpdateGravityKickFlag();
    FileOutPutConstantInterval();
}

int main_TestSymmetrizedPotentialError_ColdCollapse(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    ReadParametersForSymmetrizedPotentialError();

    InitColdCollapseTestMixed(Nstart,0.1,1.0/1.0,ON);
    Pall.dt_const = exp2(-3.0);
    InitLogFiles();
    InitializeRun();
    FreeStructPbody();

    GravityRoot.OpeningAngle = 0.5;

    GravityRoot.NumberofLeavesInGroup = 128;
    for(int j=MassRatioStart;j<=MassRatioEnd;j*=MassRatioStep){
        FILE *fp;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"./dedt.1.%02d.%03d.dat",j,128);
        FileOpen(fp,Fname,"w");
        for(int k=-3;k>-11;k--){
            InitColdCollapseTestMixed(Nstart,0.1,1.0/(double)j,ON);
            Pall.dt_const = exp2((double)k);
            gprintlmpi(Pall.dt_const);

            InitializeRunPotentialError();
            CalcSymmetrizedPotential();
            double E_init = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++) 
                E_init += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel) + Pbody[i]->Pot;

            Run();

            CalcSymmetrizedPotential();
            double E_end = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++) 
                E_end += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel) + Pbody[i]->Pot;
    
            fprintf(fp,"%1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g\n",
                    Pall.dt_const,fabs((E_init-E_end)/E_init),E_init,E_end,
                    Pbody[0]->Mass,Pbody[0]->Eps,Pbody[1]->Mass,Pbody[1]->Eps,cbrt(1.0/Nstart));

            FreeStructPbody();
        }
        fclose(fp);
    }


    GravityRoot.NumberofLeavesInGroup = 8;
    for(int j=MassRatioStart;j<=MassRatioEnd;j*=MassRatioStep){
        FILE *fp;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"./dedt.1.%02d.%03d.dat",j,8);
        FileOpen(fp,Fname,"w");
        for(int k=-3;k>-11;k--){
            InitColdCollapseTestMixed(Nstart,0.1,1.0/(double)j,ON);
            Pall.dt_const = exp2((double)k);
            gprintlmpi(Pall.dt_const);

            InitializeRunPotentialError();
            CalcSymmetrizedPotential();
            double E_init = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++) 
                E_init += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel) + Pbody[i]->Pot;

            Run();

            CalcSymmetrizedPotential();
            double E_end = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++) 
                E_end += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel) + Pbody[i]->Pot;

            fprintf(fp,"%1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g\n",
                    Pall.dt_const,fabs((E_init-E_end)/E_init),E_init,E_end,
                    Pbody[0]->Mass,Pbody[0]->Eps,Pbody[1]->Mass,Pbody[1]->Eps,cbrt(1.0/Nstart));

            FreeStructPbody();
        }
        fclose(fp);
    }


    GravityRoot.OpeningAngle = 0.e0;
    for(int j=MassRatioStart;j<=MassRatioEnd;j*=MassRatioStep){
        FILE *fp;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"./dedt.direct.1.%02d.dat",j);
        FileOpen(fp,Fname,"w");
        for(int k=-3;k>-11;k--){
            InitColdCollapseTestMixed(Nstart,0.1,1.0/(double)j,ON);
            Pall.dt_const = exp2((double)k);
            gprintlmpi(Pall.dt_const);

            InitializeRunPotentialError();
            CalcSymmetrizedPotential();
            double E_init = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++) 
                E_init += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel) + Pbody[i]->Pot;

            Run();

            CalcSymmetrizedPotential();
            double E_end = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++) 
                E_end += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel) + Pbody[i]->Pot;
    
            fprintf(fp,"%1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g %1.15g\n",
                    Pall.dt_const,fabs((E_init-E_end)/E_init),E_init,E_end,
                    Pbody[0]->Mass,Pbody[0]->Eps,Pbody[1]->Mass,Pbody[1]->Eps,cbrt(1.0/Nstart));

            FreeStructPbody();
        }
        fclose(fp);
    }


    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR

int main_ColdCollapseTest_ActiveParticles(void){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    //InitColdCollapseTest(1024<<10,0.1,OFF);
    //InitSphereWithMeanParticleDistance(1024<<13);
    //InitColdCollapseTest(4*1024,0.1);
    //ReadCosmologicalData();
    //InitUniformSphereTest(88);
    //InitUniformSphereTest(64);
    //InitUniformSphereTest(128);
    //InitUniformSphereTest(168);
    //InitUniformSphereTest(256);
    //InitRandomBoxTest(6250000);
    //InitRandomBoxTest(4000000*MPIGetNumProcs());
    //InitRandomBoxTest(4000000);
    //InitRandomBoxTest(1000000);

    ReadOldCosmologicalData(10);

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    double TimeTree = GetElapsedTime();
    for(int i=0;i<Iter;i++)
        PlantGravityTree();
    TimeTree = GetElapsedTime()-TimeTree;

    long long int count_f = 0;
    for(int i=0;i<Pall.Ntotal;i++)
        count_f += Pbody[i]->InteractionList;
    fprintf(stderr,"[%02d] =F= Just after the decomposition = %llu, %lld\n",
            MPIGetMyID(),count_f/Pall.Ntotal,count_f);

    GravityRoot.OpeningAngle = 0.75;
    for(int i=0;i<Pall.Ntotal;i++)
        Pbody[i]->Active = ON;
    Pall.NActives_t = Pall.Ntotal_t;

    PlantGravityTree();
    ClearLocalAccPot();
    ForceParallelTreeGRAPE();

    for(int i=0;i<3;i++){
        PlantGravityTree();
        ClearLocalAccPot();
        ForceParallelTreeGRAPE();
        DomainDecomposition();
    }

#if 0
    {
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"TestActives.%08ld.%02d",Pall.Ntotal_t,MPIGetNumProcs());
    if(MPIGetMyID()==0)
        FileOpen(fp,fname,"w");
    MPI_Barrier(MPI_COMM_WORLD);
    int Stride = 1;
    //while(Stride < Pall.Ntotal_t){
    while(Stride < 2){
        double Time,TimeGadget,TimeLET;
        double CheckAccPot[4];
        double CheckAccPotGlobal[4];

        for(int i=0;i<Pall.Ntotal;i++)
            Pbody[i]->Active = OFF;
        int NActives = 0;
        for(int i=0;i<Pall.Ntotal;i++){
            if(Pbody[i]->GlobalID%Stride == 0){
                Pbody[i]->Active = ON;
                NActives ++;
            }
        }
        int GlobalActives;
        MPI_Allreduce(&NActives,&GlobalActives,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        Pall.NActives_t = GlobalActives;

        //PlantTreePartialUpdate(TreeForGravFine);
        //PlantTreeFullUpdate(TreeForGravFine);

        Time = GetElapsedTime();
        for(int i=0;i<Iter;i++){
            ClearLocalAccPot();
            ForceTreeGRAPE_PreCommunication();
        }
        TimeGadget = GetElapsedTime()-Time;
        CheckAccPot[0] = CheckAccPot[1] = CheckAccPot[2] = CheckAccPot[3] = 0.e0;
        for(int i=0;i<Pall.Ntotal;i++){
            if(Pbody[i]->Active){
                CheckAccPot[0] += Pbody[i]->Mass*Pbody[i]->Acc[0];
                CheckAccPot[1] += Pbody[i]->Mass*Pbody[i]->Acc[1];
                CheckAccPot[2] += Pbody[i]->Mass*Pbody[i]->Acc[2];
                CheckAccPot[3] += Pbody[i]->Pot;
            }
        }
        MPI_Allreduce(CheckAccPot,CheckAccPotGlobal,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(MPIGetMyID()==0)
            fprintf(stderr,"[Type Gadget] Check = %g, %g, %g, %g\n",
                CheckAccPotGlobal[0],CheckAccPotGlobal[1],
                    CheckAccPotGlobal[2],CheckAccPotGlobal[3]);

        Time = GetElapsedTime();
        for(int i=0;i<Iter;i++){
            ClearLocalAccPot();
            ForceForestTreeGRAPE();
        }
        TimeLET = GetElapsedTime()-Time;
        CheckAccPot[0] = CheckAccPot[1] = CheckAccPot[2] = CheckAccPot[3] = 0.e0;
        for(int i=0;i<Pall.Ntotal;i++){
            if(Pbody[i]->Active){
                CheckAccPot[0] += Pbody[i]->Mass*Pbody[i]->Acc[0];
                CheckAccPot[1] += Pbody[i]->Mass*Pbody[i]->Acc[1];
                CheckAccPot[2] += Pbody[i]->Mass*Pbody[i]->Acc[2];
                CheckAccPot[3] += Pbody[i]->Pot;
            }
        }
        MPI_Allreduce(CheckAccPot,CheckAccPotGlobal,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(MPIGetMyID()==0)
            fprintf(stderr,"[Type LET] Check = %g, %g, %g, %g\n",
                CheckAccPotGlobal[0],CheckAccPotGlobal[1],
                    CheckAccPotGlobal[2],CheckAccPotGlobal[3]);
        if(MPIGetMyID()==0)
            fprintf(fp,"%ld %ld %g %g %g\n",
                Pall.NActives_t,Pall.Ntotal_t,
                    TimeGadget/(double)Iter,TimeLET/(double)Iter,TimeTree/(double)Iter);

        Stride *= 2;
    }

    if(MPIGetMyID()==0)
        fclose(fp);
    }
#endif
#if 1
    //double Theta[] = {0.1,0.25,0.50,0.75,1.0};
    double Theta[] = {1.0,0.75,0.50,0.25,0.10};
    int Ntheta = sizeof(Theta)/sizeof(double);
    dprintlmpi(Ntheta);
    int CurrentTheta = 0;
    FILE *fp;

    Ntheta = 1;
    Theta[0] = 0.75;

    while(CurrentTheta<Ntheta){
        GravityRoot.OpeningAngle = Theta[CurrentTheta];

        MPI_Barrier(MPI_COMM_WORLD);
        double TimeTree = GetElapsedTime();
        for(int i=0;i<Iter;i++)
            PlantGravityTree();
        TimeTree = GetElapsedTime()-TimeTree;
        //TimeTree *= Iter; 

        char fname[MaxCharactersInLine];
        sprintf(fname,"TestActives.%08ld.%02d.T%03d",Pall.Ntotal_t,MPIGetNumProcs(),(int)(100*Theta[CurrentTheta]));
        if(MPIGetMyID()==0){
            //dprintlmpi(MPIGetMyID());
            FileOpen(fp,fname,"w");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int Stride = 1;
        //while(Stride < Pall.Ntotal_t){
        while(Stride < 2){
            double Time,TimeGadget,TimeLET;
            double CheckAccPot[4];
            double CheckAccPotGlobal[4];

            for(int i=0;i<Pall.Ntotal;i++)
                Pbody[i]->Active = OFF;
            int NActives = 0;
            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->GlobalID%Stride == 0){
                    Pbody[i]->Active = ON;
                    NActives ++;
                }
            }
            Pall.NActives = NActives;
            int GlobalActives;
            MPI_Allreduce(&NActives,&GlobalActives,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            Pall.NActives_t = GlobalActives;

            PlantGravityTree();

            // set 1
            Time = GetElapsedTime();
            for(int i=0;i<Iter;i++){
                ClearLocalAccPot();
                ForceParallelTreeGRAPE();
                double sum = 0.e0;
                for(int i=0;i<Pall.Ntotal;i++)
                    sum += (double)(Pbody[i]->InteractionList);
                    //fprintf(stderr,"%d = %d\n", Pbody[i]->InternalLeaves , Pbody[i]->ExternalLeaves);
                double wsum;
                MPI_Allreduce(&sum,&wsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                fprintf(stderr,"[%02d] Sendrecv external = %g, %ld\n",MPIGetMyID(),sum/Pall.Ntotal,Pall.Ntotal);
                if(MPIGetMyID()==MPI_ROOT_RANK)
                    fprintf(stderr,"sum_external = %g\n",wsum/Pall.Ntotal_t);
            }
            TimeGadget = GetElapsedTime()-Time;
            CheckAccPot[0] = CheckAccPot[1] = CheckAccPot[2] = CheckAccPot[3] = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Active){
                    CheckAccPot[0] += Pbody[i]->Mass*Pbody[i]->Acc[0];
                    CheckAccPot[1] += Pbody[i]->Mass*Pbody[i]->Acc[1];
                    CheckAccPot[2] += Pbody[i]->Mass*Pbody[i]->Acc[2];
                    CheckAccPot[3] += Pbody[i]->Pot;
                }
            }
            MPI_Allreduce(CheckAccPot,CheckAccPotGlobal,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            if(MPIGetMyID()==MPI_ROOT_RANK)
                fprintf(stderr,"[Type LET Sendrecv] Check = %g, %g, %g, %g\n",
                    CheckAccPotGlobal[0],CheckAccPotGlobal[1],
                        CheckAccPotGlobal[2],CheckAccPotGlobal[3]);
#if 0
            // set 2
            Time = GetElapsedTime();
            for(int i=0;i<Iter;i++){
                ClearLocalAccPot();
                //for(int i=0;i<Pall.Ntotal;i++)
                    //Pbody[i]->InternalLeaves = Pbody[i]->ExternalLeaves = 0;
                //ForceForestTreeGRAPE_IsendIrecv();
                ForceForestTreeGRAPE_Sendrecv_SingleTree();

                double sum = 0.e0;
                for(int i=0;i<Pall.Ntotal;i++)
                    sum += (double)(Pbody[i]->InteractionList);
                    //fprintf(stderr,"%d = %d\n", Pbody[i]->InternalLeaves , Pbody[i]->ExternalLeaves);
                double wsum;
                MPI_Allreduce(&sum,&wsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                fprintf(stderr,"[%02d] SendrecvSingleTree external = %g\n",MPIGetMyID(),sum/Pall.Ntotal);
                if(MPIGetMyID()==MPI_ROOT_RANK)
                    fprintf(stderr,"sum_all = %g\n",wsum/Pall.Ntotal_t);
            }
            TimeLET = GetElapsedTime()-Time;
            CheckAccPot[0] = CheckAccPot[1] = CheckAccPot[2] = CheckAccPot[3] = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Active){
                    CheckAccPot[0] += Pbody[i]->Mass*Pbody[i]->Acc[0];
                    CheckAccPot[1] += Pbody[i]->Mass*Pbody[i]->Acc[1];
                    CheckAccPot[2] += Pbody[i]->Mass*Pbody[i]->Acc[2];
                    CheckAccPot[3] += Pbody[i]->Pot;
                }
            }
            MPI_Allreduce(CheckAccPot,CheckAccPotGlobal,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            if(MPIGetMyID()==MPI_ROOT_RANK)
                fprintf(stderr,"[Type LET SendrecvSingleTree] Check = %g, %g, %g, %g\n",
                    CheckAccPotGlobal[0],CheckAccPotGlobal[1],
                        CheckAccPotGlobal[2],CheckAccPotGlobal[3]);
#endif
            // set 3
#if 0
            Time = GetElapsedTime();
            for(int i=0;i<Iter;i++){
                ClearLocalAccPot();
                ForceForestTreeGRAPE_IsendIrecv();

                double sum = 0.e0;
                for(int i=0;i<Pall.Ntotal;i++)
                    sum += (double)(Pbody[i]->InternalLeaves + Pbody[i]->ExternalLeaves);
                double wsum;
                MPI_Allreduce(&sum,&wsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                fprintf(stderr,"sum_all = %g\n",wsum/Pall.Ntotal_t);
            }
            double TimeLETSendrecv = GetElapsedTime()-Time;
            CheckAccPot[0] = CheckAccPot[1] = CheckAccPot[2] = CheckAccPot[3] = 0.e0;
            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Active){
                    CheckAccPot[0] += Pbody[i]->Mass*Pbody[i]->Acc[0];
                    CheckAccPot[1] += Pbody[i]->Mass*Pbody[i]->Acc[1];
                    CheckAccPot[2] += Pbody[i]->Mass*Pbody[i]->Acc[2];
                    CheckAccPot[3] += Pbody[i]->Pot;
                }
            }
            MPI_Allreduce(CheckAccPot,CheckAccPotGlobal,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            if(MPIGetMyID()==0)
                fprintf(stderr,"[Type LET Sendrecv] Check = %g, %g, %g, %g\n",
                    CheckAccPotGlobal[0],CheckAccPotGlobal[1],
                        CheckAccPotGlobal[2],CheckAccPotGlobal[3]);
            if(MPIGetMyID()==0)
                fprintf(fp,"%ld %ld %g %g %g %g\n",
                    Pall.NActives_t,Pall.Ntotal_t,
                        TimeGadget/(double)Iter,TimeLET/(double)Iter,TimeLET/(double)Iter,TimeTree/(double)Iter);
#endif

            Stride *= 2;
        }

        if(MPIGetMyID()==0)
            fclose(fp);

        CurrentTheta ++;
    }
#endif

    return EXIT_SUCCESS;
}

int main_NavarroWhiteTest(void){

#if 0
    int NProcs,MyID,NameLen;
    char ProcessorName[MPI_MAX_PROCESSOR_NAME];

    // Initialize MPI.
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Get_processor_name(ProcessorName,&NameLen);

    MPISetMyID(MyID);
    MPISetNumProcs(NProcs);
    MPISetNumProcsPower(NProcs);
    MPISetNumGrapes(MIN(MPIGetNumProcs(),4));
    MPISetNameLen(NameLen);
    MPISetProcessorName(ProcessorName);

    fprintf(stderr,"Process %d on %s\n",MPIGetMyID(),MPIGetProcessorName());
#endif

    // if NX = 11, totally 484 particles generated
    // if NX = 16, totally 1736 particles generated
    // if NX = 22, totally 4776 particles generated
    // if NX = 40, totally 30976 particles generated
    // if NX = 45, totally 44413 particles generated
    // if NX = 50, totally 61432 particles generated
    // if NX = 60, totally 107392 particles generated
    // if NX = 80, totally 257776 particles generated
    // if NX = 100, totally 507376 particles generated
    // if NX = 128, totally 1072184 particles generated

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    //InitializeNavarroWhiteTest(40);
    //InitializeNavarroWhiteTest(22);
    //InitializeNavarroWhiteTest(11);
    InitializeNavarroWhiteTest(15);

    if(Pall.RunStatus == NewSimulation)
        BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    if(Pall.RunStatus == NewSimulation)
        DomainDecomposition();

    // First force / hydro calculation.

    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    if(Pall.RunStatus == NewSimulation){
        PlantGravityTree();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
    }

    InitializeRootForHydro();
    InitializeCoolingTable();

    if(Pall.RunStatus == NewSimulation){
        PlantHydroTree();
        ClearHydroData();
        CalcKernel();
        CalcDensityDivRot();
        CalcDuDtAcc();
    }

    InitializeStarFormation();

    LogOutPutEnergyMomentumAngularMomentum();
    double Einit = CalcTotalEnergyForColdCollapse();

    while(Pall.TCurrent < Pall.TEnd){
        //fprintf(stderr,"%d step\n",Pall.TStepTotal);

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        TimingResults.DecompositionThisStep = GetElapsedTime();
        DomainDecompositionOnlyDataExchange();
        TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;


        if(Pall.TCurrent >= Pall.Era){
            double Tdecomposition = GetElapsedTime();
            PreDomainDecomposition(0);
            TimingResults.DecompositionThisStep += GetElapsedTime()-Tdecomposition;
            BuildHierarchicalTimeStep();
        }

        TimingResults.IntegralThisStep = GetElapsedTime();
        RaiseActiveFlags();

        // Kick Drift and Predictors
        Kick1Drift();
        BuildPredictors();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep;

        // Calc Force
        TimingResults.GravityTreeThisStep = GetElapsedTime();
        PlantGravityTree();
        TimingResults.GravityTreeThisStep = GetElapsedTime()-TimingResults.GravityTreeThisStep;

        TimingResults.GravityThisStep = GetElapsedTime();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
        TimingResults.GravityThisStep = GetElapsedTime()-TimingResults.GravityThisStep;

        if(Pall.NActivesHydro_t>0){
            TimingResults.HydroTreeThisStep = GetElapsedTime();
            PlantHydroTree();
            TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;

            // Hydro
            ClearHydroData();
            TimingResults.HydroKernelThisStep = GetElapsedTime();
            CalcKernel();
            TimingResults.HydroKernelThisStep = GetElapsedTime()-TimingResults.HydroKernelThisStep;

            TimingResults.HydroDensityThisStep = GetElapsedTime();
            CalcDensityDivRot();
            TimingResults.HydroDensityThisStep = GetElapsedTime()-TimingResults.HydroDensityThisStep;

            TimingResults.HydroAccThisStep = GetElapsedTime();
            CalcDuDtAcc();
            TimingResults.HydroAccThisStep = GetElapsedTime()-TimingResults.HydroAccThisStep;

            // StarFormation 
            StarFormation();

            TimingResults.CoolingThisStep = GetElapsedTime();
            if(Pall.TCurrent>0.1){
                CalcCooling();
            }
            TimingResults.CoolingThisStep = GetElapsedTime()-TimingResults.CoolingThisStep;
        }

        // Kick
        double TimeLog = TimingResults.IntegralThisStep;
        TimingResults.IntegralThisStep = GetElapsedTime();
        Kick2();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep + TimeLog;

        // Calc New Time Step or Out Put Logs
        if(Pall.EraStart + Pall.EraLocal < Pall.Era){
            BuildNewTimeStep();
        }

        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            char fname_gas[MaxCharactersInLine],
                    fname_star[MaxCharactersInLine],
                    fname_dm[MaxCharactersInLine];
            sprintf(fname_gas,"%s.Hydro.%02d.%02d",Pall.ASCIIFileName,MPIGetMyID(),MPIGetNumProcs());
            sprintf(fname_star,"%s.Star.%02d.%02d",Pall.ASCIIFileName,MPIGetMyID(),MPIGetNumProcs());
            sprintf(fname_dm,"%s.DM.%02d.%02d",Pall.ASCIIFileName,MPIGetMyID(),MPIGetNumProcs());
            OutPutNavarroWhite(fname_gas,fname_star,fname_dm);

            CountNeighborNumber();
            fprintf(stderr,"Pall.Nstars, Pall.Nhydro = %ld %ld\n",Pall.Nstars,Pall.Nhydro);

            FileOutPutConstantInterval();

            EndLogs();
        }
        // Add This Step -> Total

        Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
        // Decomposition if necessary.

        Pall.TStepTotal ++;
    }
    dprintlmpi(Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();

    double Eend = CalcTotalEnergyForColdCollapse();

    char fname[MaxCharactersInLine];
    sprintf(fname,"ResultNavarroWhite.%02d.%02d",MPIGetMyID(),MPIGetNumProcs());
    //OutPutNavarroWhite(fname);

    CloseLogFiles();

    return EXIT_SUCCESS;
}

#ifdef TASK_NFW
int main_NFWrun(const int argc, char *argv[]){

#ifdef COSMOLOGICAL_RUN //{
    fprintf(stderr,"COSMOLOGICAL_RUN is unnecessary!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif // COSMOLOGICAL_RUN //}

    InitThermalConductivity();
    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    //int NParticles = 250000;
    //int M200 = 10;
    if(Pall.RunStatus == NewSimulation){
        ReadNFWInitialCondition();
        //CheckKennicutLaw(0);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        ReadNFWTdyn();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    // WriteTimeSteps();

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif

#ifdef TASK_MW
int main_MilkywayModel(const int argc, char *argv[]){

#ifdef COSMOLOGICAL_RUN //{
    fprintf(stderr,"COSMOLOGICAL_RUN is unnecessary!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif // COSMOLOGICAL_RUN //}

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    int NParticles = 1000000;
    //int NParticles = 3000000;
    if(Pall.RunStatus == NewSimulation){
        //Tdyn = ReadNFWInitialCondition(NParticles,M200);
        InitializeParticleDistributionForMilkyWay(NParticles,2.0);
        //InitializeParticleDistributionForMilkyWayWithInvRProf(NParticles,1.0);
        //InitializeParticleDistributionForMilkyWayWithHaloGas(NParticles);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif //TASK_MW


#ifdef TASK_MERGER //{
int main_IsolateDisk(const int argc, char *argv[]){

#if (CosmologicalRun == ON)
    fprintf(stderr,"CosmologicalRun is unnecessary!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        //ReadBinaryGalactICSFlag("IsolateDisk.init");
        ReadBinaryGalactICSFlag("MergingGalaxies.init");
        //ReadMergerRestart("MergingRestart.dat");
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();


    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif //}


#ifdef TASK_DICE_RUN //{
int main_DiceRun(const int argc, char *argv[]){

#ifdef COSMOLOGICAL_RUN 
    fprintf(stderr,"CosmologicalRun is unnecessary!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif //COSMOLOGICAL_RUN

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        //ReadDiceASCII("dice.asura");
        ReadDiceBinary("dice.asura");
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        // Pall.TEnd = 2*GIGAYEAR_CGS/Pall.UnitTime;
        // Pall.OutPutInterval = Pall.TEnd/400.0; 
        Pall.RunStatus = RestartSimulation;

    } else {
        exit(RunStatusError);
    }
    InitLogFiles();


    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_DICE_RUN //}

#ifdef TASK_MAKE_SMOOTHDISK //{
int main_SmoothIsolateDisk(const int argc, char *argv[]){

#if (CosmologicalRun == ON)
    fprintf(stderr,"CosmologicalRun is unnecessary!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();

    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        // ReadMergerRestart("MergingRestart.dat");
        ReadMergerRestart("MergingRestart.dat");
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_MAKE_SMOOTHDISK //}

#ifdef TASK_CLOUD_EXPLOSION
int main_CloudExplosion(const int argc, char *argv[]){

#if (CosmologicalRun == ON)
    fprintf(stderr,"CosmologicalRun is unnecessary!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        // //ReadBinaryGalactICSFlag("MergingGalaxies.init");
        // Init3DCollapseTest(22);
        ReadInitCloudExplosion();
        //Pall.Ns = 64;
        //Pall.Npm = 2;
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    fprintf(stderr,"Ntotal_t = %ld, NDM_t = %ld\n",Pall.Ntotal_t,Pall.NDM_t);
    fprintf(stderr,"Ntotal = %ld, NDM = %ld\n",Pall.Ntotal,Pall.NDM);

    //double Tdyn = 50*MEGAYEAR_CGS/Pall.UnitTime;
    //double Tdyn = 0.0;

    if(Pall.RunStatus == NewSimulation)
        BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();

    if(Pall.RunStatus == NewSimulation)
        DomainDecomposition();
    else
        PreDomainDecomposition(0);

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    if(Pall.RunStatus == NewSimulation){
        PlantGravityTree();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
    }

    // Check potential energy, chenge internal energy and add energy on the center.


    if(Pall.Nhydro_t > 0){
        InitializeRootForHydro();
        InitializeCoolingTable();
        InitializeFarUltraVioletField();

        if(Pall.RunStatus == NewSimulation){
            ClearHydroData();
            PlantHydroTree();
            CalcKernel();
            CalcDensityDivRot();
            CalcDuDtAcc();
        }
        //InitializeStarFormation();
        //InitializeDelayedSNII();
    }

    LogOutPutEnergyMomentumAngularMomentum();

    if(Pall.RunStatus == NewSimulation){
        FirstTimeStep();
#ifdef HYDRO_TIMESTEP_LIMITER 
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;
        }
        BuildHierarchicalTimeStep();
        CalcDensityDivRot();
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro = MIN(Phydro[i]->k_hydro,Phydro[i]->k_hydro_localmin+MAX_K_LOCAL);
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

    if(Pall.RunStatus == NewSimulation)
        FileOutPutConstantInterval();

    FILE *fp_log;
    FileOpen(fp_log,"dt.log","w");

    int flag = 0;
    int count_hydro = 0; 
    int count_gravity = 0; 
    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        if(Pall.TCurrent >= Pall.Era){
            TimingResults.DecompositionThisStep = GetElapsedTime();
            PreDomainDecomposition(0);
            DomainDecompositionOnlyDataExchange();
            TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;

            SortAllStructures();

            FirstTimeStep();
            BuildHierarchicalTimeStep();
        } else {
            if (10*Pall.NActives_t > Pall.Ntotal_t){
                TimingResults.DecompositionThisStep = GetElapsedTime();
                DomainDecompositionOnlyDataExchange();
                TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;
            }
            BuildNewTimeStep();
        }

        TimingResults.IntegralThisStep = GetElapsedTime();
        RaiseActiveFlags();

        // Kick Drift and Predictors
        Kick1Drift();
        BuildPredictors();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep;

        if(Pall.NActives_t>0){ // Calc Force
            TimingResults.GravityTreeThisStep = GetElapsedTime();
            PlantGravityTree();
            TimingResults.GravityTreeThisStep = GetElapsedTime()-TimingResults.GravityTreeThisStep;

            TimingResults.GravityThisStep = GetElapsedTime();
            ClearGravitationalForce();
            ForceParallelTreeGRAPE();
            TimingResults.GravityThisStep = GetElapsedTime()-TimingResults.GravityThisStep;
            count_gravity ++;
        } else {
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(stderr,"Skip Gravity\n");
        }

        if(Pall.NActivesHydro_t>0){
            TimingResults.HydroTreeThisStep = GetElapsedTime();
            PlantHydroTree();
            TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;

            // Hydro
            ClearHydroData();
            TimingResults.HydroKernelThisStep = GetElapsedTime();
            CalcKernel();
            TimingResults.HydroKernelThisStep = GetElapsedTime()-TimingResults.HydroKernelThisStep;

            TimingResults.HydroDensityThisStep = GetElapsedTime();
            CalcDensityDivRot();
            TimingResults.HydroDensityThisStep = GetElapsedTime()-TimingResults.HydroDensityThisStep;

            TimingResults.HydroAccThisStep = GetElapsedTime();
            CalcDuDtAcc();
            TimingResults.HydroAccThisStep = GetElapsedTime()-TimingResults.HydroAccThisStep;

            // StarFormation 
#if 0
            TimingResults.StarformationThisStep = GetElapsedTime();
            if(Pall.TCurrent>Tdyn)
                StarFormation();
            TimingResults.StarformationThisStep = GetElapsedTime()-TimingResults.StarformationThisStep;
#endif
            count_hydro ++;
        }

#if 0
        // Delayed SNe 
        if(Pall.NActivesStars_t>0){
            TimingResults.FeedbackThisStep = GetElapsedTime();
            DelayedSNe();
            TimingResults.FeedbackThisStep = GetElapsedTime()-TimingResults.FeedbackThisStep;
        }
#endif
#if 0
        if(Pall.NActivesHydro_t>0){
            TimingResults.CoolingThisStep = GetElapsedTime();
            if(Pall.TCurrent>Tdyn){
                CalcCooling();
            }
            TimingResults.CoolingThisStep = GetElapsedTime()-TimingResults.CoolingThisStep;
        }
#endif

        // Kick
        double TimeLog = TimingResults.IntegralThisStep;
        TimingResults.IntegralThisStep = GetElapsedTime();
        Kick2();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep + TimeLog;

        UpdateGravityKickFlag();

        // Calc New Time Step or Out Put Logs
        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            OutPutAllParticlesInASCIIFormat();

            if(Pall.Nhydro_t>0){
                CountNeighborNumber();
            }

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            FileOutPutConstantInterval();
            DataDump();

            ReportAllocatedMemorySizes();
            EndLogs();
            fflush(NULL);
            MPI_Barrier(MPI_COMM_WORLD);

        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"End Steps\n");
    }
    fclose(fp_log);

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Number of gravity/hydro update : %d/%d\n",count_gravity,count_hydro);

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}

static double SPHKernel(const double r, const double InvKerneli){

    double u = r*InvKerneli;
    const static double coef3d = M_1_PI;
    double coef = 1.e0;

    if(u<1.e0){
        return (coef*(1.e0 - (1.5)*SQ(u) + (0.75)*CUBE(u)));
    } else if (u<2.e0){
        return (coef*((0.25)*CUBE(2.e0-(u))));
    } else {
        return (0.e0);
    }
}

int main_MakeInitialConditionForCloudExplosion(const int argc, char *argv[]){

#if (CosmologicalRun == ON)
    fprintf(stderr,"CosmologicalRun is unnecessary!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    //Init3DCollapseTest(22);
    Init3DCollapseTest(40);

    //InitLogFiles();

    fprintf(stderr,"Ntotal_t = %ld, NDM_t = %ld\n",Pall.Ntotal_t,Pall.NDM_t);
    fprintf(stderr,"Ntotal = %ld, NDM = %ld\n",Pall.Ntotal,Pall.NDM);

    BuildPredictors();

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    PlantGravityTree();
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();

    InitializeRootForHydro();
    //InitializeCoolingTable();
    //InitializeFarUltraVioletField();

    ClearHydroData();
    PlantHydroTree();
    CalcKernel();
    CalcDensityDivRot();
    //CalcDuDtAcc();

#define EnergyFraction  (0.1)

    // get the total potentail energy.
    double TotalPot = 0.e0;
    for(int i=0;i<Pall.Ntotal;i++){
        TotalPot += Pbody[i]->Pot;
        PbodyHydro(i)->U *= 0.001;
    }
    fprintf(stderr,"Total Potential Energy = %g\n",TotalPot);

    // Search a density peak and add thermal energy around the peak.
    int index = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Rho > Phydro[index]->Rho)
            index = i;
    }
    // then index is the id of the gas with density peak.
    dprintlmpi(index);
    fprintf(stderr,"Center = %g %g %g\n",
            PhydroBody(index)->Pos[0],PhydroBody(index)->Pos[1],PhydroBody(index)->Pos[2]);

    int Neighbors[MaxNeighborSize];
    int Nlist = GetNeighborsLimited(PhydroBody(index)->Pos,2.0*Phydro[index]->Kernel,Neighbors);

    double wt = 0;
    for(int leaf=0;leaf<Nlist;leaf++){
        double r = DISTANCE(PhydroBody(index)->Pos,PhydroBody(Neighbors[leaf])->Pos);
        double w = SPHKernel(r,1.e0/Phydro[index]->Kernel);
        wt += w;
    }
    double iwt = 1.e0/wt;
    eprintlmpi(iwt);

    double Uinit = (EnergyFraction*fabs(TotalPot))/PhydroBody(index)->Mass;
    double totale = 0; 
    for(int leaf=0;leaf<Nlist;leaf++){
        double r = DISTANCE(PhydroBody(index)->Pos,PhydroBody(Neighbors[leaf])->Pos);
        double w = SPHKernel(r,1.e0/Phydro[index]->Kernel);
        PbodyHydro(Neighbors[leaf])->U += Uinit*w*iwt;
        totale += Uinit*w*iwt;
    }
    eprintlmpi(totale);
    eprintlmpi(totale*PhydroBody(index)->Mass);

    // write
    // file 名は、CE.Nparticles.fraction
    //          CE.10000.010 //1
    //                   100 //10
    //                   001 //0.1

    char fname[MaxCharactersInLine];
    FILE *fp;
    // ASCII
    Snprintf(fname,"CE.%05d.%03d",(int)(Pall.Ntotal),(int)(10*EnergyFraction));
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%g %g %g %g %g %g %g %g %g\n",
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                PhydroBody(i)->Mass,Phydro[i]->Kernel,Phydro[i]->U);
    }
    fclose(fp);
    // BIN
    Snprintf(fname,"CE.%05d.%03d.bin",(int)Pall.Ntotal,(int)(10*EnergyFraction));
    FileOpen(fp,fname,"wb");
    struct StructWriteCE{
        float Pos[3];
        float Vel[3];
        float Mass,Kernel,U;
    } WriteCE;
    int Number = Pall.Ntotal;
    fwrite(&Number,sizeof(int),1,fp);
    for(int i=0;i<Pall.Nhydro;i++){
        WriteCE.Pos[0] = PhydroBody(i)->Pos[0];
        WriteCE.Pos[1] = PhydroBody(i)->Pos[1];
        WriteCE.Pos[2] = PhydroBody(i)->Pos[2];
        WriteCE.Vel[0] = PhydroBody(i)->Vel[0];
        WriteCE.Vel[1] = PhydroBody(i)->Vel[1];
        WriteCE.Vel[2] = PhydroBody(i)->Vel[2];
        WriteCE.Mass = PhydroBody(i)->Mass;
        WriteCE.Kernel = Phydro[i]->Kernel;
        WriteCE.U = Phydro[i]->U;
        fwrite(&WriteCE,sizeof(struct StructWriteCE),1,fp);
    }
    fclose(fp);



    double Ek=0.e0,Ep=0.e0,Ethermal=0.e0;
    for(int i=0;i<Pall.Ntotal;i++){
        Ek += Pbody[i]->Mass*NORM2(Pbody[i]->Vel);
        Ep += Pbody[i]->Pot;
        Ethermal += Pbody[i]->Mass*PbodyHydroU(i);
    }
    fprintf(stderr,"Et,Ek,Ep,Eth = %g %g %g %g\n",
            Ek+Ep+Ethermal,Ek,Ep,Ethermal); 

    return ;
}
#endif

int main_IsothermalSphericalCollapse(void){

#if 0
    int NProcs,MyID,NameLen;
    char ProcessorName[MPI_MAX_PROCESSOR_NAME];

    // Initialize MPI.
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Get_processor_name(ProcessorName,&NameLen);

    MPISetMyID(MyID);
    MPISetNumprocs(NProcs);
    MPISetNumProcsPower(NProcs);
    MPISetNumgrapes(MIN(MPIGetNumProcs(),4));
    MPISetNamelen(NameLen);
    MPISetProcessorName(ProcessorName);

    fprintf(stderr,"Process %d on %s\n",MPIGetMyID(),MPIGetProcessorName());
#endif

#ifndef ISOTHERMALEOS_RUN
    fprintf(stderr,"ISOTHERMALEOS_RUN flag is necessary!!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return EXIT_SUCCESS;
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    //InitIsothermalSphericalCollapseTest(20);
    InitStandardIsothermalTestCase(40);

    if(Pall.RunStatus == NewSimulation)
        BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    if(Pall.RunStatus == NewSimulation)
        DomainDecomposition();

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    if(Pall.RunStatus == NewSimulation){
        PlantGravityTree();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
    }

    InitializeRootForHydro();
    if(Pall.RunStatus == NewSimulation){
        ClearHydroData();
        PlantHydroTree();
        CalcKernel();
        CalcDensityDivRot();
        CalcDuDtAcc();
    }

    FILE *fdp;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        FileOpen(fdp,"DensityPeak.dat","w");

    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        TimingResults.DecompositionThisStep = GetElapsedTime();
        DomainDecompositionOnlyDataExchange();
        TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;

        if(Pall.TCurrent >= Pall.Era){
            double Tdecomposition = GetElapsedTime();
            PreDomainDecomposition(0);
            TimingResults.DecompositionThisStep += GetElapsedTime()-Tdecomposition;
            BuildHierarchicalTimeStep();
        }

        TimingResults.IntegralThisStep = GetElapsedTime();
        RaiseActiveFlags();

        // Kick Drift and Predictors
        Kick1Drift();
        BuildPredictors();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep;

        // Calc Force
        TimingResults.GravityTreeThisStep = GetElapsedTime();
        PlantGravityTree();
        TimingResults.GravityTreeThisStep = GetElapsedTime()-TimingResults.GravityTreeThisStep;

        TimingResults.GravityThisStep = GetElapsedTime();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
        TimingResults.GravityThisStep = GetElapsedTime()-TimingResults.GravityThisStep;

        if(Pall.NActivesHydro_t>0){ // Hydro Part
            TimingResults.HydroTreeThisStep = GetElapsedTime();
            PlantHydroTree();
            TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;

            ClearHydroData();

            TimingResults.HydroKernelThisStep = GetElapsedTime();
            CalcKernel();
            TimingResults.HydroKernelThisStep = GetElapsedTime()-TimingResults.HydroKernelThisStep;

            TimingResults.HydroDensityThisStep = GetElapsedTime();
            CalcDensityDivRot();
            TimingResults.HydroDensityThisStep = GetElapsedTime()-TimingResults.HydroDensityThisStep;

            TimingResults.HydroAccThisStep = GetElapsedTime();
            CalcDuDtAcc();
            TimingResults.HydroAccThisStep = GetElapsedTime()-TimingResults.HydroAccThisStep;

        }

        // Kick
        double TimeLog = TimingResults.IntegralThisStep;
        TimingResults.IntegralThisStep = GetElapsedTime();
        Kick2();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep + TimeLog;

        // Calc New Time Step or Out Put Logs
        if(Pall.EraStart + Pall.EraLocal < Pall.Era){
            BuildNewTimeStep();
        }

        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            OutPutASCIIDATA();
            FileOutPutConstantInterval();
            EndLogs();


            /* Search and write density peak. */
            double Peak = 0.e0;
            for(int i=0;i<Pall.Nhydro;i++){
                if(Peak < Phydro[i]->Rho){
                    Peak = Phydro[i]->Rho;
                } 
            }

            double GlobalPeak;
            MPI_Allreduce(&Peak,&GlobalPeak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(fdp,"%g %g\n",Pall.TCurrent,
                    (Pall.UnitMass/CUBE(Pall.UnitLength))*GlobalPeak);
            fflush(NULL);
        }

        FileOutputIsothermalSphericalCollapse();

        Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
        // Decomposition if necessary.

        Pall.TStepTotal ++;
    }
    dprintlmpi(Pall.TStepTotal);

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fclose(fdp);

#if 0
    FILE *fp;
    char fname_r[MaxCharactersInLine];
    sprintf(fname_r,"./data/Isothermal.Result.%02d",MPIGetMyID());
    FileOpen(fp,fname_r,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }   
    fclose(fp);
#endif

#if 0
    char fname[MaxCharactersInLine];
    sprintf(fname,"./data/Result3D.%02d",MPIGetMyID());
    OutPut3DCollapseTest(fname);
#endif

    return EXIT_SUCCESS;
}

#ifdef TASK_M2_COLLAPSE
int main_M2SphericalCollapse(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        //InitM2SphericalCollapseTest(50);
        InitM2SphericalCollapseTest(200);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_M2_COLLAPSE

#ifdef TASK_TURBULENCE
int main_turbulence(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        InitTurbulentCloud("./TC.init",1,0.5,0.1);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_TURBULENCE

#ifdef TASK_ROTATINGDISK_WITH_SINK
int main_rotatingdisk_with_sink(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        InitRotatingDiskWithSink(5000);
        //InitTurbulentCloud();
        //InitM2SphericalCollapseTest(40);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif //TASK_ROTATINGDISK_WITH_SINK

#ifdef TASK_AGNTORUS

int NGasParticles_AGN;
double Redge_AGN_CGS;
double MassBH_AGN_CGS;
double MassGas_AGN_CGS;
int PRESNFLAG;
double PRESN_SNR; // in Msun/year.
double PRESN_DURATION; // in year.
int DISTORTION_FLAG;

static void ReadAGNParams(char *fname){

    if(CheckFile(fname)){
        FILE *fp;
        FileOpen(fp,fname,"r");
        fscanf(fp,"%d",&NGasParticles_AGN);
        fscanf(fp,"%le",&Redge_AGN_CGS);
        fscanf(fp,"%le",&MassGas_AGN_CGS);
        fscanf(fp,"%le",&MassBH_AGN_CGS);
        fscanf(fp,"%d",&PRESNFLAG);
        fscanf(fp,"%le",&PRESN_SNR);
        fscanf(fp,"%le",&PRESN_DURATION);
        fscanf(fp,"%d",&DISTORTION_FLAG);
        fclose(fp);

        Redge_AGN_CGS *= PC_CGS;
        MassBH_AGN_CGS *= MSUN_CGS;
        MassGas_AGN_CGS *= MSUN_CGS;

        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"Parameters obtained in %s.\n",fname);
    } else {
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"Since there is no parameter file, fidutial values are adopted.\n");

        NGasParticles_AGN = 100000;
        Redge_AGN_CGS = 32.0*PC_CGS;
        MassGas_AGN_CGS = 6.e+6*MSUN_CGS;
        MassBH_AGN_CGS = 1.3e+7*MSUN_CGS;
        PRESNFLAG = 1;
        PRESN_SNR = 1.0;
        PRESN_DURATION = 3.e+6;
        DISTORTION_FLAG = 0;
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"NGasParticles_AGN = %d\n",NGasParticles_AGN);
        fprintf(stderr,"Redge = %g [pc]\n",Redge_AGN_CGS/PC_CGS);
        fprintf(stderr,"MGas = %g [Msun]\n",MassGas_AGN_CGS/MSUN_CGS);
        fprintf(stderr,"MBH = %g [Msun]\n",MassBH_AGN_CGS/MSUN_CGS);
        fprintf(stderr,"SN flag = %d\n",PRESNFLAG);
        fprintf(stderr,"SNR  = %g\n",PRESN_SNR);
        fprintf(stderr,"SNDuration = %g\n",PRESN_DURATION);
    }
    return ;
}

int main_AGNTorus(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    ReadAGNParams("./AGNParams.txt");

    if(Pall.RunStatus == NewSimulation){
        //InitAGNTorus(250000);
        InitAGNTorus(NGasParticles_AGN,DISTORTION_FLAG);
        //OutPutAllParticlesInASCIIFormat();
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();
    InitializeDelayedSNII();
    if(PRESNFLAG == 1){
        ArtificialStarFormation(PRESN_SNR,PRESN_DURATION);
        //ArtificialStarFormation(1.,3.e+6);
        //ArtificialStarFormation(0.1,3.e+6);
    }

    /*
    fflush(NULL);
    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++)
        if(Pstar[i]->TypeIIProb)
            counter ++;
    dprintlmpi(counter);
    */

#if 0
    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->TypeIIProb)
            counter ++;
    }
    dlprintlmpi(Pall.Nstars);
    dprintlmpi(counter);
    int wcounter = 0;
    MPI_Allreduce(&counter,&wcounter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    dprintlmpi(wcounter);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    exit(1);
#endif

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_AGNTORUS

#ifdef TASK_GALACTIC_CENTER
int GCCloudType;
double GCCloudPower;
int GCCloudN;
double GCCloudR;
double GCCloudM;
double GCMBH;
double GCEPSBH;
double GCEATHALOGAS;
double GCR_peri;
double GCDay_peri;
double GCEccentricity;
double GCTorbit;
double GCDay_start;
double GCDay_end;
int GCNorbit;
char GCOutputDir[MaxCharactersInLine];
// halo gas parameter
int GCHaloType; // Type 0 simple hydrostatic
int GCHaloN;
int GCHaloRotation; // if 1 halo has angular momentum.
                    // This flag works only when GCHaloType == 0
double GCHaloRotationFactor; // if GCHaloRotation == 1, halo rotates with ThisParam X keplar.
int GCHaloSpinVectorChange; // if 1, change the spin vector's direction.
double GCHaloSpinTheta;
double GCHaloSpinPhi;
double GCHaloEtaHot;
double GCHaloTruncR;
// halo gas parameter

static void ReadGCParams(char fname[]){

    if(CheckFile(fname)){
        FILE *fp;
        FileOpen(fp,fname,"r");
        fscanf(fp,"%s",GCOutputDir);
        fscanf(fp,"%d",&GCCloudType);
        fscanf(fp,"%le",&GCCloudPower);
        fscanf(fp,"%d",&GCCloudN);
        fscanf(fp,"%le",&GCCloudR);
        fscanf(fp,"%le",&GCCloudM);

        fscanf(fp,"%le",&GCMBH);
        fscanf(fp,"%le",&GCEPSBH);
        fscanf(fp,"%d",&GCEATHALOGAS);
        fscanf(fp,"%le",&GCR_peri);
        fscanf(fp,"%le",&GCDay_peri);
        fscanf(fp,"%le",&GCEccentricity);
        fscanf(fp,"%le",&GCTorbit);
        fscanf(fp,"%le",&GCDay_start);
        fscanf(fp,"%le",&GCDay_end);

        fscanf(fp,"%d",&GCNorbit);

        fscanf(fp,"%d",&GCHaloType);
        fscanf(fp,"%d",&GCHaloN);
        fscanf(fp,"%d",&GCHaloRotation);
        fscanf(fp,"%le",&GCHaloRotationFactor);
        fscanf(fp,"%d",&GCHaloSpinVectorChange);
        fscanf(fp,"%le",&GCHaloSpinTheta);
        fscanf(fp,"%le",&GCHaloSpinPhi);
        fscanf(fp,"%le",&GCHaloEtaHot);
        fscanf(fp,"%le",&GCHaloTruncR);

        GCCloudR *= AU_CGS;
        GCCloudM *= MEARTH_CGS;
        GCMBH *= MSUN_CGS;
        GCEPSBH *= AU_CGS;

        fclose(fp);
    } else {
        // Use fiducial parameters.
        sprintf(GCOutputDir,"./data");
        GCCloudType = 0; // 0=uniform,
        GCCloudPower = -2; //
        GCCloudN = 10000;// number of particles.
        GCCloudR = 125*AU_CGS;// in AU
        GCCloudM = 3*MEARTH_CGS;// in g
        GCMBH  =  4.31e+6*MSUN_CGS;//  # in g. SrgA* mass.
        GCEPSBH  = 10*AU_CGS;//  # in g. SrgA* mass.
        GCEATHALOGAS  = 1;//  # 1 BH eats halo gas
        GCR_peri = 4e+15;//    # cm.
        GCDay_peri = 2013.5;// # year.
        GCEccentricity = 0.938;// #
        GCNorbit = 1000;//
        GCDay_start = -10.0;// # year
        GCDay_end   = +10.0;// # year
        GCTorbit = 137.77;//  # orbital period in year. 

        GCHaloType = 0; // Halo type
        GCHaloN = 10000; // Halo particle number
        GCHaloRotation = 0;     // Halo rotation parameter.
        GCHaloRotationFactor = 0.0;     // Halo rotation parameter.
        GCHaloSpinVectorChange = 0;     // Halo rotation parameter.
        GCHaloSpinTheta = 0.0;  // Halo spin parameter: theta.
        GCHaloSpinPhi = 0.0;    // Halo spin parameter: phi.
        GCHaloEtaHot = 0.5; // Fiducial value, from Burkert et al. 2012, ApJ, 750, 58
        GCHaloTruncR = 1.26e+17; // Truncation radius of halo. 
                                 // Fiducial value is apocenter distance.
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(!CheckDir(GCOutputDir)){
            MakeDir(GCOutputDir);
        }
        // Write Parameters.
        FILE *fp;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"%s/Param.Log",GCOutputDir);
        FileOpen(fp,Fname,"w");

        fprintf(fp,"%s #GCOutputDir\n",GCOutputDir);
        fprintf(fp,"%d #GCCloudType\n",GCCloudType);
        fprintf(fp,"%le #GCCloudPower\n",GCCloudPower);
        fprintf(fp,"%d #GCCloudN\n",GCCloudN);
        fprintf(fp,"%le #GCCloudR\n",GCCloudR);
        fprintf(fp,"%le #GCCloudM\n",GCCloudM);

        fprintf(fp,"%le #GCMBH\n",GCMBH);
        fprintf(fp,"%le #GCEPSBH\n",GCEPSBH);
        fprintf(fp,"%d #GCEATHALOGAS\n",GCEATHALOGAS);
        fprintf(fp,"%le #GCR_peri\n",GCR_peri);
        fprintf(fp,"%le #GCDay_peri\n",GCDay_peri);
        fprintf(fp,"%le #GCEccentricity\n",GCEccentricity);
        fprintf(fp,"%le #GCTorbit\n",GCTorbit);
        fprintf(fp,"%le #GCDay_start\n",GCDay_start);
        fprintf(fp,"%le #GCDay_end\n",GCDay_end);

        fprintf(fp,"%d #GCNorbit\n",GCNorbit);

        fprintf(fp,"%d #GCHaloType\n",GCHaloType);
        fprintf(fp,"%d #GCHaloN\n",GCHaloN);
        fprintf(fp,"%d #GCHaloRotation\n",GCHaloRotation);
        fprintf(fp,"%le #GCHaloRotationFactor\n",GCHaloRotationFactor);
        fprintf(fp,"%d #GCHaloSpinVectorChange\n",GCHaloSpinVectorChange);
        fprintf(fp,"%le #GCHaloSpinTheta\n",GCHaloSpinTheta);
        fprintf(fp,"%le #GCHaloSpinPhi\n",GCHaloSpinPhi);
        fprintf(fp,"%le #GCHaloEtaHot\n",GCHaloEtaHot);
        fprintf(fp,"%le #GCHaloTruncR\n",GCHaloTruncR);
        fclose(fp);
    }

    return ;
}

static double ReturnGCSemiMajorAxis(void){
    return GCR_peri/(1.0-GCEccentricity);
}
static double ReturnGCSemiMinorAxis(void){
        return GCR_peri/(1.0+GCEccentricity);
} 
static double ReturnGCOrbitalPeriod(void){
    double a_in_cgs = ReturnGCSemiMajorAxis();
    double M_in_cgs = GCMBH;
    return 2.0*M_PI*sqrt(CUBE(a_in_cgs)/(GRAVITY_CONSTANT_CGS*M_in_cgs));
}
static double ReturnGCTFromEta(double Eta, double Tr){
    return (Tr/(2.0*M_PI))*(Eta-GCEccentricity*sin(Eta));
}
static double ReturnGCPhiFromEta(const double Eta){
    return 2.0*atan(sqrt((1.0+GCEccentricity)/(1.0-GCEccentricity)) * tan(0.5*Eta) );
}
static double ReturnGCEtaFromT(double Tstart){

    double Eta_down = -M_PI;
    double Eta_up = M_PI;
    double Tr = ReturnGCOrbitalPeriod();

    int counter =0;
    do{
        double Eta_Center = 0.5*(Eta_up+Eta_down);
        if(ReturnGCTFromEta(Eta_Center,Tr) < Tstart){
            Eta_down = Eta_Center;
        } else if(ReturnGCTFromEta(Eta_Center,Tr) > Tstart){
            Eta_up = Eta_Center;
        }
        // fprintf(stderr,"Eta iteration: %g %g %g | %g %g\n",Eta_Center,Tstart,ReturnTFromEta(Eta_Center,Tr),Eta_down,Eta_up);

        counter ++;
    } while (fabs(Eta_down-Eta_up)>1.e-6);

    return Eta_down;
}


static void CalcGCParmetersOnOrbit(void){

    double a_in_cgs =  ReturnGCSemiMajorAxis();
    double M_in_cgs = GCMBH;
    double Tr = ReturnGCOrbitalPeriod();

    fprintf(stderr,"Tr = %g [sec], %g [yr]\n",Tr,Tr/YEAR_CGS);

    FILE *fp;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"%s/Orbit.dat",GCOutputDir);
    FileOpen(fp,Fname,"w");

    double deta = M_PI/(double)GCNorbit;
    for(int i=0;i<GCNorbit;i++){
        double eta = deta * i;
        double r = a_in_cgs*(1.0-GCEccentricity*cos(eta));
        double T = ReturnGCTFromEta(eta,Tr);
        double v = a_in_cgs*GCEccentricity*sin(eta)*2.0*M_PI/(1.0-GCEccentricity*cos(eta))/Tr;
        double Phi_Now = ReturnGCPhiFromEta(eta);
        fprintf(fp,"%g %g %g %g %g %g\n",eta,+T,r*cos(Phi_Now),-r*sin(Phi_Now),-v*cos(Phi_Now),-v*sin(Phi_Now));
    }
    for(int i=0;i<GCNorbit;i++){
        double eta = deta * i;
        double r = a_in_cgs*(1.0-GCEccentricity*cos(eta));
        double T = ReturnGCTFromEta(eta,Tr);
        double v = a_in_cgs*GCEccentricity*sin(eta)*2.0*M_PI/(1.0-GCEccentricity*cos(eta))/Tr;
        double Phi_Now = ReturnGCPhiFromEta(eta);
        fprintf(fp,"%g %g %g %g %g %g\n",eta,+T,r*cos(Phi_Now),-r*sin(Phi_Now),-v*cos(Phi_Now),-v*sin(Phi_Now));
    }
    fclose(fp);

    return ;
}

static int GCTStep = 0;
static double ReturnGCAcc(const double Pos[], const int Dim){
    double dist = NORM(Pos);
    return -(GRAVITY_CONSTANT_CGS*GCMBH/CUBE(dist))*Pos[Dim];
}

static void CalcGCAcc(const double Pos[], double Acc[]){

    for(int k=0;k<3;k++){
        Acc[k] = ReturnGCAcc(Pos,k);
    }
    return ;
}

static void GCKick(const double dt, double Vel[], const double Acc[]){

    for(int k=0;k<3;k++){
        Vel[k] += Acc[k]*dt;
    }
    return;
}
static void GCDrift(const double dt, double Pos[], const double Vel[]){

    for(int k=0;k<3;k++){
        Pos[k] += Vel[k]*dt;
    }
    return ;
}

static void InitParticles(double Eta_start, double Pos[], 
        double Vel[], double Acc[]){

    // Shift Particles, Add velocities.
    // r = a(1-Eccentricity*cos(eta))
    double a_in_cgs = ReturnGCSemiMajorAxis();
    double b_in_cgs = ReturnGCSemiMinorAxis();
    double Phi_start = ReturnGCPhiFromEta(Eta_start);
    double R = a_in_cgs*(1.0-SQ(GCEccentricity))
        /(1.0+GCEccentricity*cos(Phi_start));
    double Tr = ReturnGCOrbitalPeriod();
    double V = a_in_cgs*GCEccentricity*sin(Eta_start)
        *2.0*M_PI/(1.0-GCEccentricity*cos(Eta_start))/Tr;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        gprintlmpi((Tr/(2.0*M_PI))*(Eta_start-GCEccentricity*sin(Eta_start))/YEAR_CGS);
        gprintlmpi(sqrt(CUBE(a_in_cgs)/(GRAVITY_CONSTANT_CGS*GCMBH))
                *(Eta_start-GCEccentricity*sin(Eta_start))/YEAR_CGS);
        gprintlmpi(Tr/YEAR_CGS);
    }

    double p = a_in_cgs*(1.0-SQ(GCEccentricity));
    double mu = GRAVITY_CONSTANT_CGS*GCMBH;
    double vx = - sqrt(mu/p) * sin(Phi_start);
    double vy = + sqrt(mu/p) * (GCEccentricity+cos(Phi_start));
    if(MPIGetMyID() == MPI_ROOT_RANK)
        gprintlmpi(2*M_PI*sqrt(CUBE(a_in_cgs)/mu)/YEAR_CGS);

    Pos[0] = R*cos(Phi_start);
    Pos[1] = R*sin(Phi_start);
    Pos[2] = 0.e0;

    Vel[0] = vx;
    Vel[1] = vy;
    Vel[2] = 0.e0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"%s/Orbit.dat",GCOutputDir);

        FILE *fp;
        FileOpen(fp,Fname,"w");

        double dphi = M_PI/500;
        for(int i=0;i<500;i++){
            double _R = a_in_cgs*(1.0-SQ(GCEccentricity))
                /(1.0+GCEccentricity*cos(dphi*i));
            double x = _R*cos(dphi*i);
            double y = _R*sin(dphi*i);

            fprintf(fp,"%g %g %g\n",x,y,0.0);
        }
        dphi *= -1.0;
        for(int i=0;i<500;i++){
            double phi = M_PI-dphi*i;
            double _R = a_in_cgs*(1.0-SQ(GCEccentricity))
                /(1.0+GCEccentricity*cos(phi));
            double x = _R*cos(phi);
            double y = _R*sin(phi);

            fprintf(fp,"%g %g %g\n",x,y,0.0);
        }
        fclose(fp);


        Snprintf(Fname,"%s/PeriApoCenters.dat",GCOutputDir);
        FileOpen(fp,Fname,"w");
        fprintf(fp,"%g %g %g\n",a_in_cgs*(1.0-SQ(GCEccentricity))
                /(1.0+GCEccentricity),0.0,0.0);
        fprintf(fp,"%g %g %g\n",-a_in_cgs*(1.0-SQ(GCEccentricity))
                /(1.0-GCEccentricity),0.0,0.0);
        fprintf(fp,"%g %g %g\n",R*cos(Phi_start),R*sin(Phi_start),0.0);
        fprintf(fp,"%g %g %g\n",a_in_cgs*(cos(Eta_start)-GCEccentricity),
                b_in_cgs*sin(Eta_start),0.0);
        fprintf(fp,"Velocities of the pericenter and apocenter\n");
        fprintf(fp,"%g, %g [cm/s]\n",
                sqrt(mu/p) * (GCEccentricity+cos(0.0)),
                sqrt(mu/p) * (GCEccentricity+cos(M_PI)));
        fprintf(fp,"%g, %g [km/s]\n",
                1.e-5*sqrt(mu/p) * (GCEccentricity+cos(0.0)),
                1.e-5*sqrt(mu/p) * (GCEccentricity+cos(M_PI)));
        fclose(fp);
    }

    return ;
}

static void IntegrateGCOrbit(double InitPos[], double InitVel[]){

    double Tstart = GCDay_start*YEAR_CGS; // in sec.

    double Eta_start = ReturnGCEtaFromT(Tstart);       
    if(MPIGetMyID() == MPI_ROOT_RANK){
        gprintl(Tstart);
        gprintl(Tstart/YEAR_CGS);
        gprintl(Eta_start);
    }

    double Pos[3];
    double Vel[3];
    double Acc[3];
    InitParticles(Eta_start,Pos,Vel,Acc);
    CalcGCAcc(Pos,Acc);

    for(int k=0;k<3;k++){
        InitPos[k] = Pos[k]; // in cm.
        InitVel[k] = Vel[k]; // in cm/s.
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"%g %g\n",Pos[0],Pos[1]);
        fprintf(stderr,"%g %g\n",InitPos[0],InitPos[1]);
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"%s/Integral.log",GCOutputDir);
        FileOpen(fp,Fname,"w");

        double t = 0.0;
        double t_end = (GCDay_end-GCDay_start)*YEAR_CGS;
        double dt = t_end/(double)GCNorbit;
        do{
            GCKick(0.5*dt,Vel,Acc);
            GCDrift(dt,Pos,Vel);
            CalcGCAcc(Pos,Acc);
            GCKick(0.5*dt,Vel,Acc);

            fprintf(fp,"%g %g %g %g %g %g %g %g %g %g\n",t,
                Pos[0],Pos[1],Pos[2],Vel[0],Vel[1],Vel[2],Acc[0],Acc[1],Acc[2]);

            t += dt;
            GCTStep ++;
        } while(t<t_end);

        fclose(fp);

    fprintf(stderr,"Time integration was finished. t = %g [yr] -> %g [yr] with %d steps\n",
            GCDay_start+GCDay_peri,GCDay_end+GCDay_peri,GCNorbit);
    }

    return ;
}

int main_GalacticCenter(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    ReadGCParams("./GCParam.txt");
    double InitPos[3],InitVel[3];
    IntegrateGCOrbit(InitPos,InitVel);

    // exit(1);
    if(Pall.RunStatus == NewSimulation){
        InitGCCloud(InitPos,InitVel);
        // exit(1);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();
    InitializeDelayedSNII();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_GALACTIC_CENTER


#ifdef TASK_SANTABARBARA
int main_SantaBarbara(const int argc, char *argv[]){

    // Test Cosmological Time 
    InitializeCommunicationTable();
    InitializeCommunicationBuffers();

    // Read Cosmological Initial Conditions
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        ReadSantaBarbaraInitialCondition(1);
        //ReadCosmologicalData();
        //ReadMultiMassNbodyInitialCondition();
        //ReadMultiMassCosmologicalInitialCondition("TheSmallBox.Body.dat","TheSmallBox.Boundary.dat");
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    char fname[MaxCharactersInLine];
    sprintf(fname,"%s.final.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    ParallelWriteAllData(fname);

    return EXIT_SUCCESS;
}
#endif // TASK_SANTABARBARA


#ifdef TASK_GALAXY_FORMATION //{
int main_GalaxyFormation(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        ReadMusic("music.asura",1);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();
    UpdateCosmologicalParameters();


    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_GALAXY_FORMATION //}

int main_BreakingDam(void){

#if 0
    int NProcs,MyID,NameLen;
    char ProcessorName[MPI_MAX_PROCESSOR_NAME];

    // Initialize MPI.
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Get_processor_name(ProcessorName,&NameLen);

    MPISetMyID(MyID);
    MPISetNumprocs(NProcs);
    MPISetNumProcsPower(NProcs);
    MPISetNumgrapes(MIN(MPIGetNumProcs(),4));
    MPISetNamelen(NameLen);
    MPISetProcessorName(ProcessorName);

    fprintf(stderr,"Process %d on %s\n",MPIGetMyID(),MPIGetProcessorName());
#endif

    fprintf(stderr,"Must not use the self-gravity routine.\n");

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();

#define MultiplicationConstantForDam    (1) 
    InitBreakingDam(MultiplicationConstantForDam);

    if(Pall.RunStatus == NewSimulation)
        BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    if(Pall.RunStatus == NewSimulation)
        DomainDecomposition();

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    if(Pall.RunStatus == NewSimulation){
        PlantGravityTree();
        ClearGravitationalForce();
        GravitationalForceFromEarth();
    }

    InitializeRootForHydro();
    if(Pall.RunStatus == NewSimulation){
        ClearHydroData();
        PlantHydroTree();
        CalcKernel();
        CalcDensityDivRot();
        CalcDuDtAcc();
    }

    FILE *fdp;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        FileOpen(fdp,"HightandSurgeFront.dat","w");

    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        TimingResults.DecompositionThisStep = GetElapsedTime();
        DomainDecompositionOnlyDataExchange();
        TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;

        if(Pall.TCurrent >= Pall.Era){
            double Tdecomposition = GetElapsedTime();
            PreDomainDecomposition(0);
            TimingResults.DecompositionThisStep += GetElapsedTime()-Tdecomposition;
            BuildHierarchicalTimeStep();
        }

        TimingResults.IntegralThisStep = GetElapsedTime();
        RaiseActiveFlags();

        // Kick Drift and Predictors
        Kick1Drift();
        BuildPredictors();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep;

        // Calc Force
        TimingResults.GravityTreeThisStep = GetElapsedTime();
        PlantGravityTree();
        TimingResults.GravityTreeThisStep = GetElapsedTime()-TimingResults.GravityTreeThisStep;

        TimingResults.GravityThisStep = GetElapsedTime();
        ClearGravitationalForce();
        GravitationalForceFromEarth();
        TimingResults.GravityThisStep = GetElapsedTime()-TimingResults.GravityThisStep;

        if(Pall.NActivesHydro_t>0){ // Hydro Part
            TimingResults.HydroTreeThisStep = GetElapsedTime();
            PlantHydroTree();
            TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;

            ClearHydroData();
            TimingResults.HydroKernelThisStep = GetElapsedTime();
            CalcKernel();
            TimingResults.HydroKernelThisStep = GetElapsedTime()-TimingResults.HydroKernelThisStep;

            TimingResults.HydroDensityThisStep = GetElapsedTime();
            CalcDensityDivRot();
            TimingResults.HydroDensityThisStep = GetElapsedTime()-TimingResults.HydroDensityThisStep;

            TimingResults.HydroAccThisStep = GetElapsedTime();
            CalcDuDtAcc();
            TimingResults.HydroAccThisStep = GetElapsedTime()-TimingResults.HydroAccThisStep;

        }

        // add anti force from boundary regions
        AntiForceFromBounaryRegions(MultiplicationConstantForDam);

        // clear acc and hydro acc for bounary particles.
#if 1
        for(int i=0;i<Pall.Nhydro;i++){
            if(PhydroBody(i)->GlobalID == -1){
                Phydro[i]->HydroAcc[0] = 
                Phydro[i]->HydroAcc[1] =
                Phydro[i]->HydroAcc[2] = 
                PhydroAcc(i)[0] = 
                PhydroAcc(i)[1] =
                PhydroAcc(i)[2] = 0.e0;
            }
        }
#endif
        // Kick
        double TimeLog = TimingResults.IntegralThisStep;
        TimingResults.IntegralThisStep = GetElapsedTime();
        Kick2();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep + TimeLog;

        // Calc New Time Step or Out Put Logs
        if(Pall.EraStart + Pall.EraLocal < Pall.Era){
            BuildNewTimeStep();
        }

        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            OutPutASCIIDATA();
            FileOutPutConstantInterval();
            EndLogs();

            /* Search and write density peak. */
            double Peak = 0.e0;
            for(int i=0;i<Pall.Nhydro;i++){
                if(Peak < Phydro[i]->Rho){
                    Peak = Phydro[i]->Rho;
                } 
            }

            double GlobalPeak;
            MPI_Allreduce(&Peak,&GlobalPeak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            if(MPIGetMyID() == MPI_ROOT_RANK)
                fprintf(fdp,"%g %g\n",Pall.UnitTime*Pall.TCurrent/1.0774e+12,
                    (Pall.UnitMass/CUBE(Pall.UnitLength))*GlobalPeak);
            fflush(NULL);
        }

        //FileIO3DCollapse();

        Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
        // Decomposition if necessary.

        Pall.TStepTotal ++;
    }
    dprintlmpi(Pall.TStepTotal);

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fclose(fdp);

#if 1
    FILE *fp;
    char fname_r[MaxCharactersInLine];
    sprintf(fname_r,"./data/.Result.%02d",MPIGetMyID());
    FileOpen(fp,fname_r,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif

#if 1
    char fname[MaxCharactersInLine];
    sprintf(fname,"./data/Result3D.%02d",MPIGetMyID());
    OutPut3DCollapseTest(fname);
#endif

    return EXIT_SUCCESS;
}


int main_CosmologicalRun(const int argc, char *argv[]){

#if 0
    int NProcs,MyID,NameLen;
    char ProcessorName[MPI_MAX_PROCESSOR_NAME];

    // Initialize MPI.
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Get_processor_name(ProcessorName,&NameLen);

    MPISetMyID(MyID);
    MPISetNumProcs(NProcs);
    MPISetNumProcsPower(NProcs);
    MPISetNumGrapes(MIN(MPIGetNumProcs(),4));
    MPISetNameLen(NameLen);
    MPISetProcessorName(ProcessorName);

    fprintf(stderr,"Process %d on %s\n",MPIGetMyID(),MPIGetProcessorName());
#endif

    // Test Cosmological Time 
    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    // Read Cosmological Initial Conditions
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        //ReadCosmologicalData();
        ReadCosmologicalHalfwayData("convert006");
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation)
        BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    if(Pall.RunStatus == NewSimulation)
        DomainDecomposition();
    else
        PreDomainDecomposition(0);

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    if(Pall.RunStatus == NewSimulation){
        PlantGravityTree();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
    }

    InitializeRootForHydro();
    InitializeCoolingTable();
    //InitializeMinimumTemperature(100.0,PhydroMass(0));

    if(Pall.RunStatus == NewSimulation){
        ClearHydroData();
        PlantHydroTree();
        CalcKernel();
        CalcDensityDivRot();
        CalcDuDtAcc();
    }

    InitializeStarFormation();
    InitializeDelayedSNII();

    if(Pall.RunStatus == RestartSimulation)
        BuildHierarchicalTimeStep();

    LogOutPutEnergyMomentumAngularMomentum();

    int TheLastASCIIOutputTimeStep = 0;
    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        /* Cosmological Paramters in this step */
        Pall.Redshift = CalcCurrentTtoZ();
        UpdateAdaptiveSofteningFactor();
        Pall.OverDensityZ = CalcOverDensityZ();
        Pall.RHOcZ = CalcVirializedDensityZ();
        Pall.HubbleZ = CalcHubbleParameterZ();
        /* Cosmological Paramters in this step */

        if(Pall.TCurrent >= Pall.Era){
            TimingResults.DecompositionThisStep = GetElapsedTime();
            PreDomainDecomposition(0);
            DomainDecompositionOnlyDataExchange();
            TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;
            BuildHierarchicalTimeStep();

            SortAllStructures();
        } else if (10*Pall.NActives_t > Pall.Ntotal_t){
            TimingResults.DecompositionThisStep = GetElapsedTime();
            DomainDecompositionOnlyDataExchange();
            TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;
        }


        TimingResults.IntegralThisStep = GetElapsedTime();
        RaiseActiveFlags();

        // Kick Drift and Predictors
        Kick1Drift();
        BuildPredictors();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep;

        // Calc Force
        TimingResults.GravityTreeThisStep = GetElapsedTime();
        PlantGravityTree();
        TimingResults.GravityTreeThisStep = GetElapsedTime()-TimingResults.GravityTreeThisStep;

        TimingResults.GravityThisStep = GetElapsedTime();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
        TimingResults.GravityThisStep = GetElapsedTime()-TimingResults.GravityThisStep;

        if(Pall.NActivesHydro_t>0){
            TimingResults.HydroTreeThisStep = GetElapsedTime();
            PlantHydroTree();
            TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;

            // Hydro
            ClearHydroData();
            TimingResults.HydroKernelThisStep = GetElapsedTime();
            CalcKernel();
            TimingResults.HydroKernelThisStep = GetElapsedTime()-TimingResults.HydroKernelThisStep;

            TimingResults.HydroDensityThisStep = GetElapsedTime();
            CalcDensityDivRot();
            TimingResults.HydroDensityThisStep = GetElapsedTime()-TimingResults.HydroDensityThisStep;

            TimingResults.HydroAccThisStep = GetElapsedTime();
            CalcDuDtAcc();
            TimingResults.HydroAccThisStep = GetElapsedTime()-TimingResults.HydroAccThisStep;

            // StarFormation 
            StarFormation();
        }

        if(Pall.NActivesStars_t>0)
            DelayedSNe();

        if(Pall.NActivesHydro_t>0){
            TimingResults.CoolingThisStep = GetElapsedTime();
            CalcCooling();
            TimingResults.CoolingThisStep = GetElapsedTime()-TimingResults.CoolingThisStep;
        }

        // Kick
        double TimeLog = TimingResults.IntegralThisStep;
        TimingResults.IntegralThisStep = GetElapsedTime();
        Kick2();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep + TimeLog;

        // Calc New Time Step or Out Put Logs
        if(Pall.EraStart + Pall.EraLocal < Pall.Era){
            BuildNewTimeStep();
        }

        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
#ifdef PRINT_LOG_THIS_STEP
        LogsThisTimeStep();
#endif
        CountDirectNeighborNumber();
        CountNeighborNumber();
        EndLogs();

        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();
            LogStarFormationRate();

            if( Pall.TStepTotal - TheLastASCIIOutputTimeStep > 10 ){
                OutPutAllParticlesInASCIIFormat();
                TheLastASCIIOutputTimeStep = Pall.TStepTotal;
            }

            CountNeighborNumber();
            fprintf(stderr,"Pall.Nstars, Pall.Nhydro = %ld %ld\n",Pall.Nstars,Pall.Nhydro);

            FileOutPutConstantInterval();
            EndLogs();

            BinaryDump();

            fflush(NULL);
        }
        // Add This Step -> Total
        Pall.TCurrent = Pall.EraStart + Pall.EraLocal;

        Pall.TStepTotal ++;
    }
    dprintlmpi(Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();

    CloseLogFiles();

    return EXIT_SUCCESS;
}

int main_NeighborSearchTest(void){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    //ReadCosmologicalData();

    //InitUniformSphereTest(128);
    //InitUniformSphereTest(100);
    InitUniformSphereTest(60);
    //InitUniformSphereTest(40);
    //InitUniformSphereTest(22);

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    // First force / hydro calculation.
    ClearHydroData();
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    InitializeRootForHydro();
    PlantHydroTree();

    ClearHydroData();
    CalcKernel();

    int Neighbors[MaxNeighborSize];
    int Ntotal;
    // Gather
    double tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=0;i<Pall.Nhydro;i++){
            Ntotal += GetNeighborsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-G- [%02d] Nmean = %ld, Ntotal = %d, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,GetElapsedTime()-tstart_each);
    }
    double tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-G- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);


    // Gather and Scatter
    tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=0;i<Pall.Nhydro;i++){
            Ntotal += GetNeighborsPairsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-GS- [%02d] Nmean = %ld, Ntotal = %d, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,GetElapsedTime()-tstart_each);
    }
    tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-GS- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);

    return EXIT_SUCCESS;
}

int main_ReadAndPlantTreeTest(const int argc, char *argv[]){

    if(MPIGetNumProcs()!= 1){
        fprintf(stderr,"NProcs is not 1 but %d\n",MPIGetNumProcs());
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if(argc != 2){
        fprintf(stderr,"Usage : One input file name is necessary!\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    FILE *fp;
    FileOpen(fp,argv[1],"r");
    int count = 0;
    while(fscanf(fp,"%*e %*e %*e %*e\n") != EOF)
        count ++;
    fclose(fp);
    dprintlmpi(count);

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));

    Pall.Ntotal = Pall.Ntotal_t = count;
    int AllocationSize = count; 
    PbodySize = count;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
    }


    double mass = 1.0/(double)count;
    double eps = 0.1*1.e-3;
    double Vel_init = 1.e-8;

    FileOpen(fp,argv[1],"r");
    for(int i=0;i<count;i++){
        fscanf(fp,"%le %le %le %*e",&(Pbody[i]->Pos[0]),&(Pbody[i]->Pos[1]),&(Pbody[i]->Pos[2]));
        Pbody[i]->Mass = mass;
        Pbody[i]->Eps = eps;
        Pbody[i]->Vel[0] = Pbody[i]->Vel[1] = Pbody[i]->Vel[2] = Vel_init;
        Pbody[i]->Use = ON;
        Pbody[i]->Active = ON;
        Pbody[i]->Type = TypeDM;
        Pbody[i]->GlobalID = i;
    }
    fclose(fp);

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // PlantTree
    InitializeRootForGravity();
    sprint("Start PlantTree");
    double Time = GetElapsedTime();
    PlantGravityTree();
    Time = GetElapsedTime()-Time;
    fprintf(stderr,"The time of make tree = %g\n",Time);
    fprintf(stderr," :The loaded file name is %s\n",argv[1]);

    // GetForce
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();


    Pall.AdaptiveSofteningFactor = 1.e0;
    sprint("Start Force Calculation");
    Time = GetElapsedTime();
    ForceParallelTreeGRAPE();
    Time = GetElapsedTime()-Time;
    fprintf(stderr,"The force calculcation time = %g\n",Time);

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;
    Pall.GravConst  = 1.e0;


    return EXIT_SUCCESS;
}

static void CountActiveParticles(void){

    Pall.NActives = Pall.NActivesDM = Pall.NActivesHydro = Pall.NActivesStars = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            Pall.NActives ++;
            if(Pbody[i]->Type == TypeHydro){
                Pall.NActivesHydro ++;
                PbodyHydro(i)->Active = ON;
            } else if(Pbody[i]->Type == TypeStar){
                Pall.NActivesStars ++;
            } else if(Pbody[i]->Type == TypeDM){
                Pall.NActivesDM ++;
            }
        }
    }

    UpdateTotalActiveNumber();

    return ;
}

int main_BenchMarkTest(void){

    //zdprintlmpi(sizeof(int));
    //zdprintlmpi(sizeof(long int));
    //zdprintlmpi(sizeof(long long int));

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();

    //InitRandomBoxTestHydro(1000000);
    ReadOldCosmologicalDataFull(48);
    CountActiveParticles();

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    SortAllStructures();

#define Test (MIN(1000,Pall.Nhydro))
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // New version from here.

    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"\n");
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);

    /*
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->Kernel = 
        Phydro[i]->KernelPred = 
        10*PhydroBody(i)->Eps;
    }
    */

    InitializeRootForHydro();
    double tnbtree = GetElapsedTime();
    PlantHydroTree();
    tnbtree = GetElapsedTime()-tnbtree;
    fprintf(stderr,"[%d] New Plant Tree %g [sec]\n",MPIGetMyID(),tnbtree);

    int Nlist = 0;
    int Neighbors[MaxNeighborSize];
    tnbtree = GetElapsedTime();
    for(int i=0;i<Test;i++){
        Nlist += GetNeighborsLimited(PhydroPosP(i),2.e0*Phydro[i]->KernelPred,Neighbors);
        //fprintf(stderr,"[%d] ",Nlist);
        //for(int k=0;k<Nlist;k++)
            //NIndex += Neighbors[k];
    }
    fprintf(stderr,"\n");
    tnbtree = GetElapsedTime()-tnbtree;
    fprintf(stderr,"[%d] New NeighborSearchNew %g [sec]\n",MPIGetMyID(),tnbtree);
    fprintf(stderr,"[%d] New Sum of Neighbors = %d\n",MPIGetMyID(),Nlist);
    //fprintf(stderr,"[%d] Sum of Neighbor Indexes = %d\n",MPIGetMyID(),NIndex);

    ClearHydroData();

    double tkenrel_s = GetElapsedTime();
    CalcKernel();
    fprintf(stderr,"[%d] New End Kernel %g [sec]\n",MPIGetMyID(),GetElapsedTime()-tkenrel_s);

    int total_neighbor = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        total_neighbor += Phydro[i]->Nlist;
    }
    int global_total_neighbor = 0;
    MPI_Allreduce(&total_neighbor,&global_total_neighbor,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"New total number of neighbors = %d\n",global_total_neighbor);


    /////////////////////////////////////////////////////////////////
    // select actives
#if 0
    for(int i=0;i<Pall.Ntotal>>10;i++){
        int leaf = gsl_rng_uniform_int(RandomGenerator,Pall.Ntotal);
        Pbody[leaf]->Active = false;
        PbodyHydro(leaf)->Active = false;
    }
    PlantHydroTree();
#endif
    /////////////////////////////////////////////////////////////////


    TimingResults.HydroNeighborSearchThisStep = 0.e0;
    MPI_Barrier(MPI_COMM_WORLD);
    double trho_f = GetElapsedTime();
    CalcDensityDivRot();
    MPI_Barrier(MPI_COMM_WORLD);
    trho_f = GetElapsedTime()-trho_f;
    fprintf(stderr,"[%d] New End Hydro First loop %g [sec], NeighborSearch = %g [sec], Total = %g [sec]\n",
            MPIGetMyID(),trho_f-TimingResults.HydroNeighborSearchThisStep,
                TimingResults.HydroNeighborSearchThisStep,trho_f);

    // insert dummy data into non-active particles.
    for(int i=0;i<Pall.Ntotal;i++){
        if(!(Pbody[i]->Active)){
            PbodyHydro(i)->Rho = 0.125;
            PbodyHydro(i)->RhoPred = 0.125;
        }
    }   

    TimingResults.CoolingThisStep = 
    TimingResults.IntegralThisStep = 
    TimingResults.HydroNeighborSearchThisStep = 0.e0;

    MPI_Barrier(MPI_COMM_WORLD);
    double tacc_f = GetElapsedTime();
    CalcDuDtAcc();
    MPI_Barrier(MPI_COMM_WORLD);
    tacc_f = GetElapsedTime()-tacc_f;

    fprintf(stderr,"[%d] New End Hydro Second loop %g [sec], NeighborSearch = %g [sec], Total = %g [sec]\n",
            MPIGetMyID(),tacc_f-TimingResults.HydroNeighborSearchThisStep,
                TimingResults.HydroNeighborSearchThisStep,tacc_f);

    double rho_even = 0.e0;
    double rho_odd = 0.e0;
    double acc_even = 0.e0;
    double acc_odd = 0.e0;

    for(int i=0;i<Pall.Nhydro;i+=2){
        if(Phydro[i]->Active){
            rho_odd += Phydro[i]->Rho;
            acc_odd += NORM(Phydro[i]->HydroAcc);
        }
    }
    for(int i=1;i<Pall.Nhydro;i+=2){
        if(Phydro[i]->Active){
            rho_even += Phydro[i]->Rho;
            acc_even += NORM(Phydro[i]->HydroAcc);
        }
    }
    double LocalPack[4],GlobalPack[4];
    LocalPack[0] = rho_odd;
    LocalPack[1] = rho_even;
    LocalPack[2] = acc_odd;
    LocalPack[3] = acc_even;
    MPI_Allreduce(LocalPack,GlobalPack,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"HydroCheck %g %g %g %g\n",GlobalPack[0],GlobalPack[1],GlobalPack[2],GlobalPack[3]);

#if 0 // output hydro
    {
    FILE *fp_check_hydro;
    char fname_check_hydro[MaxCharactersInLine];
    sprintf(fname_check_hydro,"CheckHydroNew.%08ld.%02d.%02d",Pall.Nhydro_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_hydro,fname_check_hydro,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active)
            fprintf(fp_check_hydro,"%ld %g %g %g %d %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                Phydro[i]->Nlist,Phydro[i]->Kernel,
                Phydro[i]->Rho,Phydro[i]->Div,NORM(Phydro[i]->Rot),
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
    }
    fclose(fp_check_hydro);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID()==MPI_ROOT_RANK){
        char fname[MaxCharactersInLine];
        sprintf(fname,"cat CheckHydroNew.%08ld.%02d.?? | sort -n > CheckHydroNew.%08ld.%02d",
                Pall.Nhydro_t,MPIGetNumProcs(),Pall.Nhydro_t,MPIGetNumProcs());
        system(fname);
        sprintf(fname,"rm -rf CheckHydroNew.%08ld.%02d.??",Pall.Nhydro_t,MPIGetNumProcs());
        system(fname);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    ReportAllocatedMemorySizes();
    fflush(NULL);

    int l = 1;
    int NActives = 0;
    do{
        for(int i=0;i<Pall.Nhydro;i++){
            PhydroBody(i)->Active = Phydro[i]->Active = false;
        }

        for(int i=0;i<Pall.Nhydro;i++){
            if(PhydroBody(i)->GlobalID%l == 0){
                PhydroBody(i)->Active = Phydro[i]->Active = true;
            }
        }
        CountActiveParticles();
        NActives = Pall.NActivesHydro_t;

        ClearHydroData();
        double ttree_f = GetElapsedTime();
        PlantHydroTree();
        MPI_Barrier(MPI_COMM_WORLD);
        ttree_f = GetElapsedTime()-ttree_f;

        double th_f = GetElapsedTime();
        CalcKernel();
        MPI_Barrier(MPI_COMM_WORLD);
        th_f = GetElapsedTime()-th_f;
        
        double trho_f = GetElapsedTime();
        CalcDensityDivRot();
        MPI_Barrier(MPI_COMM_WORLD);
        trho_f = GetElapsedTime()-trho_f;

        double tacc_f = GetElapsedTime();
        CalcDuDtAcc();
        MPI_Barrier(MPI_COMM_WORLD);
        tacc_f = GetElapsedTime()-tacc_f;

        long int Nlist = 0;
        for(int i=0;i<Pall.Nhydro;i++){
            if(Phydro[i]->Active)
                Nlist += Phydro[i]->Nlist;
        }
        long int GNlist;
        MPI_Allreduce(&Nlist,&GNlist,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"Actives = %d, GNlist = %ld, time = %g / %g / %g / %g\n",
                    NActives,GNlist,ttree_f,th_f,trho_f,tacc_f);
        l = l<<1;
    }while(NActives > 10);


    return EXIT_SUCCESS;
}

int main_BenchMarkTestForForce(void){

    zdprintlmpi(sizeof(int));
    zdprintlmpi(sizeof(long int));
    zdprintlmpi(sizeof(long long int));
    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
#if 0
#ifdef MAKE_INIT 
    //InitRandomBoxTestHydro(100000);
    ReadOldCosmologicalDataFull(48);
#else
    //ReadCosmologicalDataForBenchParallel();
    //ReadKeisokuBenchParallel();
    ReadKeisokuCheckParallel(100000);
#endif
#endif
    InitRandomBoxTestHydro(100000);
    //ReadOldCosmologicalDataFull(48);

    // select actives
#if 0
    for(int i=0;i<Pall.Ntotal>>2;i++){
        int leaf = gsl_rng_uniform_int(RandomGenerator,Pall.Ntotal);
        Pbody[leaf]->Active = false;
    }
#endif

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    SortAllStructures();

    InitializeRootForGravity();
    InitializeRootForLET(); 
    InitializeParallelTreeGRAPE(); 


    ///////////////////////////////////////////////////////////////////////// 
    // test for Parallel tree.  
    /////////////////////////////////////////////////////////////////////////

#define OpeningAngleForTest (TreeOpeningAngle)

    ClearGravitationalForce();
    MPI_Barrier(MPI_COMM_WORLD);


    ClearGravitationalForce();
    MPI_Barrier(MPI_COMM_WORLD);

    double TimeGravPtree = GetElapsedTime();
    //ParallelTreeNew();
    TimeGravPtree = GetElapsedTime()-TimeGravPtree;
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr,"Force Tree [%d] %g [sec]\n",MPIGetMyID(),TimeGravPtree);
    PrintForceFlops(TimeGravPtree);
     
#if 0
    {
    FILE *fp_check_grav;
    char fname_check_grav[MaxCharactersInLine];
    sprintf(fname_check_grav,"CheckGravParallelTreeNew.%03d.%08ld.%02d.%02d",
            (int)(100*RootGrav.OpeningAngle),Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_grav,fname_check_grav,"w");
    for(int i=0;i<Pall.Ntotal;i++)
        fprintf(fp_check_grav,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],Pbody[i]->Pot,
                NORM(Pbody[i]->Pos),NORM(Pbody[i]->Acc));
    fclose(fp_check_grav);
    }
#endif

    ///////////////////////////////////////////////////////////////////////// 
    // test for Parallel Tree+GRAPE.  
    /////////////////////////////////////////////////////////////////////////

    ClearGravitationalForce();
    MPI_Barrier(MPI_COMM_WORLD);

    double TimeGravPTG = GetElapsedTime();
    ForceParallelTreeGRAPE();
    MPI_Barrier(MPI_COMM_WORLD);
    TimeGravPTG = GetElapsedTime()-TimeGravPTG;
    fprintf(stderr,"Force ParallelTreeGRAPE New [%d] %g [sec]\n",MPIGetMyID(),TimeGravPTG);
    PrintForceFlops(TimeGravPTG);

#if 1
    {
    FILE *fp_check_grav;
    char fname_check_grav[MaxCharactersInLine];
    sprintf(fname_check_grav,"CheckGravParallelTreeGRAPENew.%03d.%08ld.%02d.%02d",
            (int)(100*GravityRoot.OpeningAngle),Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_grav,fname_check_grav,"w");
    for(int i=0;i<Pall.Ntotal;i++)
        if(Pbody[i]->Active)
            fprintf(fp_check_grav,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],Pbody[i]->Pot,
                NORM(Pbody[i]->Pos),NORM(Pbody[i]->Acc));
    fclose(fp_check_grav);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    Snprintf(fname_check_grav,
        "cat ./CheckGravParallelTreeGRAPENew.%03d.%08ld.%02d.?? | sort -n > ./CheckGravParallelTreeGRAPENew.%03d.%08ld.%02d",
            (int)(100*GravityRoot.OpeningAngle),Pall.Ntotal_t,MPIGetNumProcs(),
                (int)(100*GravityRoot.OpeningAngle),Pall.Ntotal_t,MPIGetNumProcs()); 
    system(fname_check_grav);
    Snprintf(fname_check_grav,
        "rm -rf ./CheckGravParallelTreeGRAPENew.%03d.%08ld.%02d.??",
            (int)(100*GravityRoot.OpeningAngle),Pall.Ntotal_t,MPIGetNumProcs()); 
    system(fname_check_grav);
    }
#endif

    return EXIT_SUCCESS;
}
//////////////////////////////////
//                              //
// Main routines for test runs. //
//                              //
//////////////////////////////////

#ifdef TASK_3D_COLLAPSE //{

static void WriteASCIIData3DCollapseTest(char fname[]){

    FILE *fp;
    char fname_r[MaxCharactersInLine];
    sprintf(fname_r,"%s.%03d",fname,MPIGetMyID());
    FileOpen(fp,fname_r,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight = Phydro[i]->Mass*Phydro[i]->U;
#elif defined(USE_SPSPH) //}//{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight = Phydro[i]->Zw;
#else // USE_SPSPH //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight = Phydro[i]->Mass;
#endif // USE_DISPH //}
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i),Phydro[i]->Nlist,
                Smoothed,Weight);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"cat %s.??? | sort -n > %s",fname,fname);
        fprintf(stderr,"%s\n",command);
        system(command);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"rm %s.???",fname);
        fprintf(stderr,"%s\n",command);
        system(command);
    }


    return ;
}

int main_3DCollapse(const int argc, char *argv[]){

#ifdef ISOTHERMAL_EOS_RUN
    fprintf(stderr,"The flag  ISOTHERMAL_REOS_RUN is not allowed in this run!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){
        //Init3DCollapseTest(22);
        Init3DCollapseTest(40);
        //Init3DCollapseTest(30);
        //Init3DCollapseTest(60);
        //Init3DCollapseTest(80);
        //Init3DCollapseTest(100);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    WriteASCIIData3DCollapseTest("3Dcollapse.check");

    Run();

    EndLogs();

    WriteASCIIData3DCollapseTest("3Dcollapse.Result");
#if 1
    //char fname[MaxCharactersInLine];
    //sprintf(fname,"Result3D.%02d",MPIGetMyID());
    OutPut3DCollapseTest("Result3D");
#endif

    ReleaseAllBuffers();

    return EXIT_SUCCESS;
}

#endif //TASK_3D_COLLAPSE //}

#ifdef TASK_BLAST_WAVE //{

struct StructDt{
    int index;
    double dt;
};

static int DtCmp(const void *x, const void *y){
    const struct StructDt *pointer1 = x;
    const struct StructDt *pointer2 = y;
    if( pointer1->dt > pointer2->dt)
        return 1;
    else if( pointer1->dt < pointer2->dt)
        return -1;
    else
        return 0;
}

int main_BlastWave(void){

#ifndef PERIODIC_RUN
    fprintf(stderr,"Periodic is not defined!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    if(Pall.RunStatus == NewSimulation){
        ReadBlastWave();
        //InitBlastWave(64);
        //InitBlastWave(16);
        //InitBlastWave(32);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    char fname[MaxCharactersInLine];
    sprintf(fname,"Result3D.%02d",MPIGetMyID());
    OutPut3DCollapseTest(fname);
    OutPutZeroPlane(fname);

    return EXIT_SUCCESS;

#if 0

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP
    InitializeDecomposition();
    DomainDecomposition();

    ClearGravitationalForce();
    InitializeRootForHydro();
    ClearHydroData();
    PlantHydroTree();
    CalcKernel();
    CalcDensityDivRot();
    CalcDuDtAcc();

    LogOutPutEnergyMomentumAngularMomentum();
    if(Pall.RunStatus == NewSimulation){
        FirstTimeStep();
#ifdef HYDRO_TIMESTEP_LIMITER 
        for(int i=0;i<Pall.Nhydro;i++){
            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;
        }
#endif // HYDRO_TIMESTEP_LIMITER
        //Pall.OutPutFileNumber = 10;
        //BuildHierarchicalTimeStep();

        FILE *fp;
        FileOpen(fp,"dt.log","w");
        for(int i=0;i<Pall.Nhydro;i++){
            //fprintf(fp,"%g %g %g %g %g\n",NORM(PhydroBody(i)->Pos),Phydro[i]->dt_hydro,Phydro[i]->Kernel,Phydro[i]->U,sqrt(Pall.GGm1*Phydro[i]->U));
            /*
            fprintf(fp,"%g %g %g %g %g %g %g %g %g %g\n",NORM(PhydroBody(i)->Pos),Phydro[i]->dt_hydro,
                    Phydro[i]->Kernel,NORM(Phydro[i]->HydroAcc),
                    TFactorCourant*sqrt(2.0*Phydro[i]->Kernel/NORM(Phydro[i]->HydroAcc)),
                    2.e0*Phydro[i]->Kernel,Phydro[i]->Vsig,
                    TFactorCourant*2.e0*Phydro[i]->Kernel/Phydro[i]->Vsig,
                    Phydro[i]->U,sqrt(Pall.GGm1*Phydro[i]->U));
            */
            fprintf(fp,"%g %g %g %g\n",Phydro[i]->dt_hydro,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2]);
        }
        fclose(fp);
    }

    // reset time-steps
    if(Pall.RunStatus == NewSimulation){
        //sort dt 
        struct StructDt Dt[Pall.Nhydro];

        for(int i=0;i<Pall.Nhydro;i++){
            Dt[i].index = i;
            Dt[i].dt = Phydro[i]->dt_hydro;
        }
        qsort(Dt,Pall.Nhydro,sizeof(struct StructDt),(int(*)(const void*, const void*))DtCmp);

        //search muximum
        fprintf(stderr,"Min and Max = %g and %g\n",Dt[0].dt,Dt[Pall.Nhydro-1].dt);

        // overwrite dts of ambient particles
        for(int i=32;i<Pall.Nhydro;i++){
            int index = Dt[i].index;
            Phydro[index]->dt_hydro = Dt[Pall.Nhydro-1].dt;
            PhydroBody(index)->dt = Dt[Pall.Nhydro-1].dt;
        }
        FILE *fp;
        FileOpen(fp,"dt_new.log","w");
        for(int i=0;i<Pall.Nhydro;i++){
            fprintf(fp,"%g %g %g %g %g %g %g %g %g %g\n",NORM(PhydroBody(i)->Pos),Phydro[i]->dt_hydro,
                    Phydro[i]->Kernel,NORM(Phydro[i]->HydroAcc),
                    TFactorCourant*sqrt(2.0*Phydro[i]->Kernel/NORM(Phydro[i]->HydroAcc)),
                    2.e0*Phydro[i]->Kernel,Phydro[i]->Vsig,
                    TFactorCourant*2.e0*Phydro[i]->Kernel/Phydro[i]->Vsig,
                    Phydro[i]->U,sqrt(Pall.GGm1*Phydro[i]->U));
        }
        fclose(fp);
    }

    UpdateGravityKickFlag();
    if(Pall.RunStatus == NewSimulation)
        FileOutPutConstantInterval();
    FileIOBlastWave();

    Run();

    EndLogs();

#if 1
    FILE *fp;
    char fname_r[MaxCharactersInLine];
    sprintf(fname_r,"BlastWave.Result.%02d",MPIGetMyID());
    FileOpen(fp,fname_r,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i),Phydro[i]->Nlist);
    }
    fclose(fp);
#endif

#if 1
    char fname[MaxCharactersInLine];
    sprintf(fname,"Result3D.%02d",MPIGetMyID());
    OutPut3DCollapseTest(fname);
    OutPutZeroPlane(fname);
#endif
#endif 

    return EXIT_SUCCESS;
}
#endif //TASK_BLAST_WAVE //}


#ifdef TASK_SINUSOIDAL_WAVE //{
static void WriteSinusoidalWave(char *fname){

    FILE *fp;
    char fname_r[MaxCharactersInLine];
    sprintf(fname_r,"%s.%02d",fname,MPIGetMyID());
    FileOpen(fp,fname_r,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroVel(i)[0],
                Phydro[i]->Kernel,Phydro[i]->U,
                Phydro[i]->Alpha,Phydro[i]->DAlpha,
                Phydro[i]->Nlist);
    }
    fclose(fp);

    return ;
}

int main_SinusoidalWave(void){

#ifndef PERIODIC_RUN
    fprintf(stderr,"Periodic is not defined!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    InitSinusoidalWave();
    InitLogFiles();

    InitializeRun();

    WriteSinusoidalWave("Sinusoidal.Init");
    Run();
    EndLogs();
    WriteSinusoidalWave("Sinusoidal.Result");

    return EXIT_SUCCESS;
}
#endif //TASK_SINUSOIDAL_WAVE //}


int main_SelfSimilarCooling(void){

#ifdef ISOTHERMAL_EOS_RUN
    fprintf(stderr,"The flag  ISOTHERMAL_REOS_RUN is not allowed in this run!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    InitSelfSimilarCooling(36);

    ClearHydroData();
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    PlantGravityTree();
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();

    InitializeRootForHydro();
    PlantHydroTree();

    CalcKernel();
    CalcDensityDivRot();
    CalcDuDtAcc();

    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        TimingResults.DecompositionThisStep = GetElapsedTime();
        DomainDecompositionOnlyDataExchange();
        TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;

        if(Pall.TCurrent >= Pall.Era){
            DomainDecomposition();
            BuildHierarchicalTimeStep();
        }

        TimingResults.IntegralThisStep = GetElapsedTime();
        RaiseActiveFlags();

        // Kick Drift and Predictors
        Kick1Drift();
        BuildPredictors();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep;

        // Calc Force
        TimingResults.GravityTreeThisStep = GetElapsedTime();
        PlantGravityTree();
        TimingResults.GravityTreeThisStep = GetElapsedTime()-TimingResults.GravityTreeThisStep;

        TimingResults.GravityThisStep = GetElapsedTime();
        ClearGravitationalForce();
        ForceParallelTreeGRAPE();
        TimingResults.GravityThisStep = GetElapsedTime()-TimingResults.GravityThisStep;

        // Hydro Part
        TimingResults.HydroTreeThisStep = GetElapsedTime();
        PlantHydroTree();
        TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;

        ClearHydroData();
        TimingResults.HydroKernelThisStep = GetElapsedTime();
        CalcKernel();
        TimingResults.HydroKernelThisStep = GetElapsedTime()-TimingResults.HydroKernelThisStep;

        TimingResults.HydroDensityThisStep = GetElapsedTime();
        CalcDensityDivRot();
        TimingResults.HydroDensityThisStep = GetElapsedTime()-TimingResults.HydroDensityThisStep;

        TimingResults.HydroAccThisStep = GetElapsedTime();
        CalcDuDtAcc();
        TimingResults.HydroAccThisStep = GetElapsedTime()-TimingResults.HydroAccThisStep;

        // Kick
        double TimeLog = TimingResults.IntegralThisStep;
        TimingResults.IntegralThisStep = GetElapsedTime();
        Kick2();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep + TimeLog;

        // Calc New Time Step or Out Put Logs
        if(Pall.EraStart + Pall.EraLocal < Pall.Era){
            BuildNewTimeStep();
        }

        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        LogsThisTimeStep();
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();

            CountNeighborNumber();
            fprintf(stderr,"Pall.Nstars, Pall.Nhydro = %ld %ld\n",Pall.Nstars,Pall.Nhydro);
        }

        FileIOSelfSimilar();

        Pall.TActiveParticles += Pall.NActives_t;
        Pall.TActiveParticlesDM += Pall.NActivesDM_t;
        Pall.TActiveParticlesHydro += Pall.NActivesHydro_t;
        Pall.TActiveParticlesStar += Pall.NActivesStars_t;

        Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
        // Decomposition if necessary.

        Pall.TStepTotal ++;
    }
    dprintlmpi(Pall.TStepTotal);

    EndLogs();

#if 1
    FILE *fp;
    char fname_r[MaxCharactersInLine];
    sprintf(fname_r,"Damp.SelfSimilarCooling.Result.%02d",MPIGetMyID());
    FileOpen(fp,fname_r,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i),Phydro[i]->Nlist);
    }
    fclose(fp);
#endif

#if 1
    char fname[MaxCharactersInLine];
    sprintf(fname,"SelfSimilarCooling.Result.%02d",MPIGetMyID());
    OutPut3DCollapseTest(fname);
#endif

    return EXIT_SUCCESS;
}

int main_3DShockTube(void){

#if 0
    int NProcs,MyID,NameLen;
    char ProcessorName[MPI_MAX_PROCESSOR_NAME];

    // Initialize MPI.
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Get_processor_name(ProcessorName,&NameLen);

    MPISetMyID(MyID);
    MPISetNumprocs(NProcs);
    MPISetNumProcsPower(NProcs);
    MPISetNumgrapes(MIN(MPIGetNumProcs(),4));
    MPISetNamelen(NameLen);
    MPISetProcessorName(ProcessorName);

    fprintf(stderr,"Process %d on %s\n",MPIGetMyID(),MPIGetProcessorName());
#endif

#ifdef ISOTHERMAL_EOS_RUN
    fprintf(stderr,"The flag ISOTHERMAEOS_RUN is not allowed in this run!!\n");
    MPI_Finalize();
    exit(EXIT_SUCCESS);
#endif

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    //Read3DShockTube(5000);
    //Read3DShockTube(5000,"ASURA_3Dshocktube_5000.dat");
    Read3DShockTube(50000,"ASURA_3Dshocktube_50000.dat");
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    InitializeRootForHydro();

    InitializeRootForHydro();
    PlantHydroTree();

    ClearHydroData();

    CalcKernel();
    CalcDensityDivRot();
    //GlassCondition3DShockTube();

    BuildPredictors();
#if 0
    {
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"Init3Dshocktube.Adaptive.%02d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->Rho,PhydroMass(i));
    }
    fclose(fp);
    exit(0);
    }
#endif
    CalcDuDtAcc();

    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        TimingResults.TotalThisStep = GetElapsedTime();

        PeriodicWrapping();

#if 0
        TimingResults.DecompositionThisStep = GetElapsedTime();
        DomainDecompositionOnlyDataExchange();
        TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;

        if(Pall.TCurrent >= Pall.Era){
            DomainDecomposition();
            BuildHierarchicalTimeStep();
        }
#endif

        if(Pall.TCurrent >= Pall.Era){
            TimingResults.DecompositionThisStep = GetElapsedTime();
            PreDomainDecomposition(0);
            DomainDecompositionOnlyDataExchange();
            TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;

            //SortAllStructures();

            BuildHierarchicalTimeStep();
        } else if (10*Pall.NActives_t > Pall.Ntotal_t){
            TimingResults.DecompositionThisStep = GetElapsedTime();
            DomainDecompositionOnlyDataExchange();
            TimingResults.DecompositionThisStep = GetElapsedTime()-TimingResults.DecompositionThisStep;
        }


        TimingResults.IntegralThisStep = GetElapsedTime();
        RaiseActiveFlags();

        // Kick Drift and Predictors
        Kick1Drift();
        BuildPredictors();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep;

        // Hydro Part
        TimingResults.HydroTreeThisStep = GetElapsedTime();
        PlantHydroTree();
        TimingResults.HydroTreeThisStep = GetElapsedTime()-TimingResults.HydroTreeThisStep;

        ClearHydroData();
        TimingResults.HydroDensityThisStep = GetElapsedTime();
        CalcKernel();
        CalcDensityDivRot();
        TimingResults.HydroDensityThisStep = GetElapsedTime()-TimingResults.HydroDensityThisStep;

        TimingResults.HydroAccThisStep = GetElapsedTime();
        CalcDuDtAcc();
        TimingResults.HydroAccThisStep = GetElapsedTime()-TimingResults.HydroAccThisStep;

        // Kick
        double TimeLog = TimingResults.IntegralThisStep;
        TimingResults.IntegralThisStep = GetElapsedTime();
        Kick2();
        TimingResults.IntegralThisStep = GetElapsedTime()-TimingResults.IntegralThisStep + TimeLog;

        // Calc New Time Step or Out Put Logs
        if(Pall.EraStart + Pall.EraLocal < Pall.Era){
            BuildNewTimeStep();
        }

        TimingResults.TotalThisStep = GetElapsedTime()-TimingResults.TotalThisStep;
        UpdateTimeLogs();
        LogsThisTimeStep();
        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            LogOutPutEnergyMomentumAngularMomentum();
            LogOutPutElapsedTime();

            CountNeighborNumber();
            fprintf(stderr,"Pall.Nstars, Pall.Nhydro = %ld %ld\n",Pall.Nstars,Pall.Nhydro);

            PeriodicWrapping();
            Pall.OutPutFileNumber = 0;
            //FileOutPutConstantInterval();
            //OutPutZeroPlane("Current.dat");
            //OutPut3DCollapseTest("Radial.Current.dat");
        }
        //FileIOBlastWave();

        Pall.TActiveParticles += Pall.NActives_t;
        Pall.TActiveParticlesDM += Pall.NActivesDM_t;
        Pall.TActiveParticlesHydro += Pall.NActivesHydro_t;
        Pall.TActiveParticlesStar += Pall.NActivesStars_t;

        Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
        // Decomposition if necessary.

        Pall.TStepTotal ++;
    }
    dprintlmpi(Pall.TStepTotal);

    EndLogs();

#if 1
    FILE *fp;
    char fname_r[MaxCharactersInLine];
    sprintf(fname_r,"3DShockTube.Result.%02d",MPIGetMyID());
    FileOpen(fp,fname_r,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],Phydro[i]->Rho,
                Phydro[i]->Kernel,Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U,
                Phydro[i]->U,PhydroMass(i),Phydro[i]->Nlist);
    }
    fclose(fp);
#endif

#if 1
    char fname[MaxCharactersInLine];
    sprintf(fname,"Result3D.%02d",MPIGetMyID());
    OutPut3DCollapseTest(fname);
    OutPutZeroPlane(fname);
#endif

    return EXIT_SUCCESS;
}


#ifdef TASK_1D_SHOCKE_TUBE
int main_ShockTube(void){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
   
    if(Pall.RunStatus == NewSimulation){
        //InitShockTube(1280);
        InitSodShockTube(1000);
        //InitShockTubeHK(400);
        //InitShockTube123Problem(1000);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }


    // char fname_s[MaxCharactersInLine];
    // sprintf(fname_s,"ShockTube.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    // OutPutShockTube(fname_s);
    char fname[MaxCharactersInLine];
    sprintf(fname,"%sInit.%02d.%d.data",
            Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    OutPutShockTube(fname);

    Run();

    sprintf(fname,"%sResult.%02d.%02d.data",
            Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    OutPutShockTube(fname);

    return EXIT_SUCCESS;
}
#endif



#ifdef TASK_HYDROSTATIC
static void ParticleDistributions(char basename[]){

    // Write final state.
    FILE *fp0;
    FILE *fp1;
    char fname0[MaxCharactersInLine],fname1[MaxCharactersInLine];
    Snprintf(fname0,"./data/hs.%s.0.dat",basename);
    Snprintf(fname1,"./data/hs.%s.1.dat",basename);
    FileOpen(fp0,fname0,"w");
    FileOpen(fp1,fname1,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight= Phydro[i]->Mass*Phydro[i]->U;
#elif defined(USE_SPSPH) //}//{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight= Phydro[i]->Zw;
#else  //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight= Phydro[i]->Mass;
#endif // USE_DISPH //}
        if(Phydro[i]->Tag == 0){
            fprintf(fp0,"%ld %g %g %g %g %d\n",PhydroBody(i)->GlobalID,PhydroPos(i)[0],PhydroPos(i)[1],
                    Smoothed,Weight,Phydro[i]->Tag);
        } else {
            fprintf(fp1,"%ld %g %g %g %g %d\n",PhydroBody(i)->GlobalID,PhydroPos(i)[0],PhydroPos(i)[1],
                    Smoothed,Weight,Phydro[i]->Tag);
        }
    }
    fflush(fp0);
    fflush(fp1);
    fclose(fp0);
    fclose(fp1);

    return ;
}

int main_HydroStatic(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
   
    if(Pall.RunStatus == NewSimulation){
        //InitHydroStatic(64,0,4); // FlagEqualMass, DensityRatio
        InitHydroStatic(64,0,4); // FlagEqualMass, DensityRatio
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();
    OutPutAllParticlesInASCIIFormat();
    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    ParticleDistributions("Init");

    Run();

    ParticleDistributions("Final");
    fflush(NULL);

    return EXIT_SUCCESS;
}
#endif //TASK_HYDROSTATIC

#ifdef TASK_KELVINHELMHOLTZ_INSTABILITY
static void ParticleDistributionsForKH(char basename[]){

    // Write final state.
    FILE *fp0;
    FILE *fp1;
    char fname0[MaxCharactersInLine],fname1[MaxCharactersInLine];
    Snprintf(fname0,"./data/%s.0.dat",basename);
    Snprintf(fname1,"./data/%s.1.dat",basename);
    FileOpen(fp0,fname0,"w");
    FileOpen(fp1,fname1,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight= Phydro[i]->Mass*Phydro[i]->U;
        double DWeight= Phydro[i]->Du;
#elif defined(USE_SPSPH) //}//{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight= Phydro[i]->Zw;
        double DWeight= Phydro[i]->DZw;
#else  //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight= Phydro[i]->Mass;
        double DWeight= 1.0;
#endif // USE_DISPH //}
        if(Phydro[i]->Tag == 0){
            fprintf(fp0,"%ld %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,PhydroPos(i)[0],PhydroPos(i)[1],
                    Smoothed,Weight,DWeight,Phydro[i]->Tag);
        } else {
            fprintf(fp1,"%ld %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,PhydroPos(i)[0],PhydroPos(i)[1],
                    Smoothed,Weight,DWeight,Phydro[i]->Tag);
        }
    }
    fflush(fp0);
    fflush(fp1);
    fclose(fp0);
    fclose(fp1);

    return ;
}
int main_KHInstability(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
   
    if(Pall.RunStatus == NewSimulation){
        //InitKH(256,1,0,2.0,8.0); // NGrid, mode, seed, density contrast, Tend.
        InitKH(256,0,0,2.0,8.0); // NGrid, mode, seed, density contrast, Tend.
        //InitKH(64,0,0,2.0,8.0); // NGrid, mode, seed, density contrast, Tend.
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();
    OutPutAllParticlesInASCIIFormat();
    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }


    ParticleDistributionsForKH("kh.Init");

    Run();

    ParticleDistributionsForKH("kh.Final");

    return EXIT_SUCCESS;
}
#endif //TASK_KELVINHELMHOLTZ_INSTABILITY

#ifdef TASK_1D_TWOFLUIDS //{
int main_TwoFluids(void){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
   
    if(Pall.RunStatus == NewSimulation){
        Init1DTwoFluids(512);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    // char fname[MaxCharactersInLine];
    // sprintf(fname,"%sInit.%02d.%d.data",
            // Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    // OutPutShockTube(fname);

    Run();

    // sprintf(fname,"%sResult.%02d.%02d.data",
            // Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    // OutPutShockTube(fname);

    return EXIT_SUCCESS;
}
#endif //TASK_TWOFLUIDS //}

#ifdef TASK_KEPLER //{
int main_Kepler(void){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
   
    if(Pall.RunStatus == NewSimulation){
        InitKepler(128);
        //InitKepler(129);
    } else if (Pall.RunStatus == RestartSimulation){
        ParallelReadAllData();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();

    return EXIT_SUCCESS;
}
#endif //TASK_KEPLER //}

///////////////////////////////////////
//                                   //
// Main routines for low-level runs. //
//                                   //
///////////////////////////////////////

int main_IOTest(const int argc, char *argv[]){

    dprintlmpi(argc);
    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    InitializeDecomposition();
    // check restart file name
    GetRunStatus(argc,argv);

    if(Pall.RunStatus == NewSimulation){
        //InitSphereWithMeanParticleDistance(1024<<7);
        InitializeNavarroWhiteTest(20);
        
        PreDomainDecomposition(0);
        DomainDecomposition();

        strcpy(Pall.BaseFileName,"FileWriteTest");
        WriteAllData();
        char fname_s[MaxCharactersInLine];
        sprintf(fname_s,"InitSphere.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
        OutPutColdCollapse(fname_s);

    }else if (Pall.RunStatus == RestartSimulation){
        ReadAllData();

        PreDomainDecomposition(0);
        DomainDecomposition();

        char fname_s[MaxCharactersInLine];
        sprintf(fname_s,"ReloadSphere.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
        OutPutColdCollapse(fname_s);
    }

    return EXIT_SUCCESS;
}

void StaticPotential(void){
#define Mh  (1.0)

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            double Dist = NORM(Pbody[i]->PosP);
            double DistCube = CUBE(Dist);
            Pbody[i]->Acc[0] = -(Pall.GravConst*Mh/DistCube)*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] = -(Pall.GravConst*Mh/DistCube)*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] = -(Pall.GravConst*Mh/DistCube)*Pbody[i]->PosP[2];
        }
    }

    return;
}

#define Ecc     0.9
#define Aaxis   1.0
void InitElliptical(void){

    FILE *fp;
    char name[MaxCharactersInLine];

    ParticleDistributionGenetator(1,MPIGetMyID());

    //Pbody[0]->Pos[0] = -1.e-1;
    Pbody[0]->Pos[0] = 1.9e0;
    Pbody[0]->Pos[1] = 0.e0;
    Pbody[0]->Pos[2] = 0.e0;
    Pbody[0]->Vel[0] = 0.e0;
    //Pbody[0]->Vel[1] = sqrt((1.e0+Ecc)/(1.e0-Ecc));
    Pbody[0]->Vel[1] = sqrt((1.e0-Ecc)/(1.e0+Ecc));
    Pbody[0]->Vel[2] = 0.e0;
    Pbody[0]->Eps = 1.e-1;

    Pbody[0]->PosP[0] = Pbody[0]->Pos[0];
    Pbody[0]->PosP[1] = Pbody[0]->Pos[1];
    Pbody[0]->PosP[2] = Pbody[0]->Pos[2];

    /*
    Pbody[1]->Pos[1] = -1.e-1;
    Pbody[1]->Pos[0] = 0.e0;
    Pbody[1]->Pos[2] = 0.e0;
    Pbody[1]->Vel[1] = 0.e0;
    Pbody[1]->Vel[0] = sqrt((1.e0+Ecc)/(1.e0-Ecc));
    Pbody[1]->Vel[2] = 0.e0;
    Pbody[1]->Eps = 1.e-1;

    Pbody[1]->PosP[1] = -1.e-1;
    Pbody[1]->PosP[0] = 0.e0;
    Pbody[1]->PosP[2] = 0.e0;
    */
#if 0
    for( int n=0; n<N; n++ ){
        M[n].p  = Vector2( 0 , sqrt((1+0.9)/(1-0.9)) );
        M[n].q  = Vector2( -0.1 , 0 );
    }
#endif

    Pall.RunStatus = NewSimulation;
    Pall.GravConst  = 1.e0;

    for( int n=0; n<Pall.Ntotal; n++ ){
        sprintf(name,"evolution.%02d",n);
        fp =fopen(name,"w");
        fclose(fp);
    }

  return;
}

void InitEllipticalMany(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    FILE *fp;
    char name[MaxCharactersInLine];

    ParticleDistributionGenetator(Number,MPIGetMyID());

    double ang = 2.0*PI/(NProcs*Number);
    double x = -0.1;
    double y = 0.e0;
    double vx = 0.e0;
    double vy = sqrt((1.e0+Ecc)/(1.e0-Ecc));
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pos[0] = x*cos(ang*(i*NProcs+MyID))-y*sin(ang*(i*NProcs+MyID));
        Pbody[i]->Pos[1] = x*sin(ang*(i*NProcs+MyID))+y*cos(ang*(i*NProcs+MyID));
        Pbody[i]->Pos[2] = 0.e0;
        Pbody[i]->Vel[0] = vx*cos(ang*(i*NProcs+MyID))-vy*sin(ang*(i*NProcs+MyID));
        Pbody[i]->Vel[1] = vx*sin(ang*(i*NProcs+MyID))+vy*cos(ang*(i*NProcs+MyID));
        Pbody[i]->Vel[2] = 0.e0;

        Pbody[i]->PosP[0] = Pbody[i]->Pos[0];
        Pbody[i]->PosP[1] = Pbody[i]->Pos[1];
        Pbody[i]->PosP[2] = Pbody[i]->Pos[2];

        Pbody[i]->Eps = 1.e-1;
    }

    dlprint(Pall.Ntotal);

    sprintf(name,"evolution.%02d",MyID);
    fp =fopen(name,"w");
    fclose(fp);

    return;
}

void OutPutEllipticalMany(const int Number){

    FILE *fp;
    char name[MaxCharactersInLine];

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->GlobalID%Number == 0){
            sprintf(name,"evolution.%02ld",Pbody[i]->GlobalID/Number);
            fp =fopen(name,"a");
            fprintf(fp,"%d %g %g %g %g %g %g %g\n",MPIGetMyID(),Pall.TCurrent+Pbody[i]->dt,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
            fclose(fp);
        }
    }
    return;
}

void OutPutElliptical(void){

    FILE *fp;
    char name[MaxCharactersInLine];

    for( int n=0; n<Pall.Ntotal; n++ ){
        sprintf(name,"evolution.%02d",n);
        fp =fopen(name,"a");
        fprintf(fp,"%g %g %g %g %g\n",Pall.TCurrent+Pbody[n]->dt,
                Pbody[n]->Pos[0],Pbody[n]->Pos[1],Pbody[n]->Vel[0],Pbody[n]->Vel[1]);
        fclose(fp);
    }
    return;
}

void InitElliptical2Particles(void){

    FILE *fp;
    char name[MaxCharactersInLine];

    ParticleDistributionGenetator(2,MPIGetMyID());

    double ang = 0.0;
    //double x = 1.9;
    double x = -0.1;
    double y = 0.e0;
    double vx = 0.e0;
    double vy = sqrt((1.e0+Ecc)/(1.e0-Ecc));
    //double vy = sqrt((1.e0-Ecc)/(1.e0+Ecc));

    Pbody[0]->Pos[0] = x*cos(ang)-y*sin(ang);
    Pbody[0]->Pos[1] = x*sin(ang)+y*cos(ang);
    Pbody[0]->Pos[2] = 0.e0;
    Pbody[0]->Vel[0] = vx*cos(ang)-vy*sin(ang);
    Pbody[0]->Vel[1] = vx*sin(ang)+vy*cos(ang);
    Pbody[0]->Vel[2] = 0.e0;

    Pbody[0]->PosP[0] = Pbody[0]->Pos[0];
    Pbody[0]->PosP[1] = Pbody[0]->Pos[1];
    Pbody[0]->PosP[2] = Pbody[0]->Pos[2];

    Pbody[0]->Eps = 1.e-1;

    ang = PI;
    x = 1.9;
    y = 0.e0;
    vx = 0.e0;
    vy = sqrt((1.e0-Ecc)/(1.e0+Ecc));

    Pbody[1]->Pos[0] = x*cos(ang)-y*sin(ang);
    Pbody[1]->Pos[1] = x*sin(ang)+y*cos(ang);
    Pbody[1]->Pos[2] = 0.e0;
    Pbody[1]->Vel[0] = vx*cos(ang)-vy*sin(ang);
    Pbody[1]->Vel[1] = vx*sin(ang)+vy*cos(ang);
    Pbody[1]->Vel[2] = 0.e0;

    Pbody[1]->PosP[0] = Pbody[1]->Pos[0];
    Pbody[1]->PosP[1] = Pbody[1]->Pos[1];
    Pbody[1]->PosP[2] = Pbody[1]->Pos[2];

    Pbody[1]->Eps = 1.e-1;

    dlprint(Pall.Ntotal);

    for(int i=0;i<Pall.Ntotal;i++){
        sprintf(name,"evolution.%02d",i);
        fp =fopen(name,"w");
        fclose(fp);
    }

    return;
}

void OutPutElliptical2Particles(void){

    FILE *fp;
    char name[MaxCharactersInLine];

    for(int i=0;i<Pall.Ntotal;i++){
        sprintf(name,"evolution.%02d",i);
        fp =fopen(name,"a");
        fprintf(fp,"%d %g %g %g %g %g %g %g\n",MPIGetMyID(),Pall.TCurrent+Pbody[i]->dt,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
        fclose(fp);
    }
    return;
}

int main_Orbit(void){

#if 0
    int NProcs,MyID,NameLen;
    char ProcessorName[MPI_MAX_PROCESSOR_NAME];

    // Initialize MPI.
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Get_processor_name(ProcessorName,&NameLen);

    MPISetMyID(MyID);
    MPISetNumprocs(NProcs);
    MPISetNumProcsPower(NProcs);
    MPISetNumgrapes(MIN(MPIGetNumProcs(),4));
    MPISetNamelen(NameLen);
    MPISetProcessorName(ProcessorName);

    fprintf(stderr,"Process %d on %s\n",MPIGetMyID(),MPIGetProcessorName());
#endif
    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    // Load Particle data.
    InitElliptical();
#define _NTEST   1
    //InitEllipticalMany(_NTEST);
    //InitElliptical2Particles();

    // First force / hydro calculation.
    StaticPotential();

    Pall.TEnd = 2*PI*200;
    //Pall.TEnd =2*PI*20;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    int step = 0;
    while(Pall.TCurrent < Pall.TEnd){

        if(Pall.TCurrent >= Pall.Era)
            BuildHierarchicalTimeStep();

        RaiseActiveFlags();

        Kick1Drift();
        BuildPredictors();

        // Calc Force
        StaticPotential();

        Kick2();

        // Calc new step
        if(Pall.EraStart + Pall.EraLocal < Pall.Era)
            BuildNewTimeStep();

        Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
        //eprintlmpi(Pall.TCurrent);

        // Decomposition if necessary.


        // Data Output.
        OutPutElliptical();
        //OutPutEllipticalMany(_NTEST);
        //OutPutElliptical2Particles();

        step ++;
        Pall.TStepTotal ++;
    }
    dprint(step);
    dprint(step/200);

    return EXIT_SUCCESS;
}

void AdjustActiveNumber(const int nstride){
    // nstride = 1 means full active case.
    // nstride = N means active particles are every N. 

    for(int i=0;i<Pall.Ntotal;i++)
        Pbody[i]->Active = OFF;
    Pall.NActives = 0;
    for(int i=0;i<Pall.Ntotal;i+=nstride){
        Pbody[i]->Active = ON;
        Pall.NActives ++;
    }

    unsigned long int GlobalActives[2],PartialActives[2];
    PartialActives[0] = Pall.NActives;
    PartialActives[1] = Pall.NActivesHydro;
    MPI_Allreduce(PartialActives,GlobalActives,2,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    Pall.NActives_t = GlobalActives[0];
    Pall.NActivesHydro_t = GlobalActives[1];

    fprintf(stderr,"NActive particles = %ld, Active fraction = %g\n",Pall.NActives,(double)Pall.NActives/(double)Pall.Ntotal);

    return;
}


int main_HydroCheckDirectAndTree(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();

    GetRunStatus(argc,argv);
    if(Pall.RunStatus == NewSimulation){ // Make an uniform sphere.
        InitUniformSphereTest(8);
    } else if (Pall.RunStatus == RestartSimulation){ // Read restart files.
        //ParallelReadAllData();
        ReadParallelDataOnSingleNode();
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }

    for(int i=0;i<Pall.Ntotal;i++){
        //if(Pbody[i]->Type == TypeHydro)
            //PbodyHydro(i)->PosP[0] = Pbody[i]->PosP[0];

        Pbody[i]->Pos[0] = Pbody[i]->Pos[1];
    }

    //if(Pall.RunStatus == NewSimulation)
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    // First force / hydro calculation.
    ClearHydroData();
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    InitializeRootForHydro();

    PlantHydroTree();
    ClearHydroData();

    // count Numbers.
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int LNHydro[NProcs];
    int WNHydro[NProcs];

    for(int i=0;i<NProcs;i++)
        LNHydro[i] = 0;
    LNHydro[MyID] = Pall.Nhydro;
    MPI_Allreduce(LNHydro,WNHydro,NProcs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    int *NodeArray;
    NodeArray = malloc(sizeof(int)*Pall.Nhydro_t);
    int count = 0;
    for(int i=0;i<NProcs;i++){
        for(int k=0;k<WNHydro[i];k++){
            NodeArray[count] = i;
            count ++;
        }
    }


#if 0
    dlprintlmpi(Pall.Nhydro_t);
    count = 0;
    for(int i=0;i<Pall.Nhydro_t;i++){
        if(MyID == MPI_ROOT_RANK)
            dprintlmpi(i);

        for(int k=0;k<Pall.Nhydro;k++){
            Phydro[k]->Active = PhydroBody(k)->Active = OFF;
        }
        if(NodeArray[i] == MyID){
            Phydro[count]->Active = 
            PhydroBody(count)->Active = ON;
            Pall.NActivesHydro = 1;
            count ++;
        } else {
            Pall.NActivesHydro = 0;
        }
        Pall.NActivesHydro_t = 1;

        PlantHydroTree();
        ClearHydroData();
        //CalcKernel();
        CalcDensityDivRot();
        for(int k=0;k<Pall.Nhydro;k++){
            if(PhydroBody(k)->Active){
                fprintf(stderr,"T NodeID = %d, ID = %d, GlobalID = %ld |+",
                        MyID,k,PhydroBody(k)->GlobalID);
                fprintf(stderr,"%d | %g %g %g %g \n",
                        Phydro[k]->Nlist,Phydro[k]->Rho,Phydro[k]->Kernel,Phydro[k]->Div,NORM(Phydro[k]->Rot));

            }
        }
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);

        ClearHydroData();
        CalcDensityDivRotDirect();
        for(int k=0;k<Pall.Nhydro;k++){
            if(PhydroBody(k)->Active){
                fprintf(stderr,"D NodeID = %d, ID = %d, GlobalID = %ld |+",
                        MyID,k,PhydroBody(k)->GlobalID);
                fprintf(stderr,"%d | %g %g %g %g \n",
                        Phydro[k]->Nlist,Phydro[k]->Rho,Phydro[k]->Kernel,Phydro[k]->Div,NORM(Phydro[k]->Rot));

            }
        }
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);

        for(int l=0;l<Pall.Nhydro;l++){
            if(PhydroActive(l)){
                Phydro[i]->Du = 
                Phydro[i]->HydroAcc[0] = 
                Phydro[i]->HydroAcc[1] = 
                Phydro[i]->HydroAcc[2] = 0.e0;
            }
        }

        CalcDuDtAcc();
        for(int k=0;k<Pall.Nhydro;k++){
            if(PhydroBody(k)->Active){
                fprintf(stderr,"T2 NodeID = %d, ID = %d, GlobalID = %ld |+",
                        MyID,k,PhydroBody(k)->GlobalID);
                fprintf(stderr,"%d | %g %g %g \n",Phydro[k]->Nlist,
                        Phydro[k]->Du,Phydro[k]->Vsig,NORM(Phydro[k]->HydroAcc));

            }
        }
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);

        for(int l=0;l<Pall.Nhydro;l++){
            if(PhydroActive(l)){
                Phydro[i]->Du = 
                Phydro[i]->HydroAcc[0] = 
                Phydro[i]->HydroAcc[1] = 
                Phydro[i]->HydroAcc[2] = 0.e0;
            }
        }

        CalcDuDtAccDirect();
        for(int k=0;k<Pall.Nhydro;k++){
            if(PhydroBody(k)->Active){
                fprintf(stderr,"D2 NodeID = %d, ID = %d, GlobalID = %ld |+",
                        MyID,k,PhydroBody(k)->GlobalID);
                fprintf(stderr,"%d | %g %g %g \n",Phydro[k]->Nlist,
                        Phydro[k]->Du,Phydro[k]->Vsig,NORM(Phydro[k]->HydroAcc));

            }
        }
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
#if 0
        for(int k=0;k<Pall.Nhydro;k++){
            if(PhydroBody(k)->Active){
                fprintf(stderr,"NodeID = %d, ID = %d, GlobalID = %ld\n",
                        MyID,k,PhydroBody(k)->GlobalID);
                fprintf(stderr,"%d | %g %g %g %g | %g %g %g \n",
                        Phydro[k]->Nlist,
                        Phydro[k]->Rho,Phydro[k]->Kernel,Phydro[k]->Div,NORM(Phydro[k]->Rot),
                        Phydro[k]->Du,Phydro[k]->U,Phydro[k]->Vsig,NORM(Phydro[k]->HydroAcc));

            }
        }
#endif
    }
#endif


#if 0
    // Set Initial kernel size.
#define InitKernelSize  (0.1)
    for(int i=0;i<Pall.Nhydro;i++)
        Phydro[i]->Kernel = Phydro[i]->KernelPred = 0.e0;

    // Tree get kernel size
    CalcKernel();


    // copy the kernel data


    // Set Initial kernel size.
    for(int i=0;i<Pall.Nhydro;i++)
        Phydro[i]->Kernel = Phydro[i]->KernelPred = 0.e0;
    // Direct get kernel size

#endif

#if 0
    fprintf(stderr,"Kernel size = %g %g\n",Phydro[0]->Kernel,Phydro[0]->KernelPred);
    fprintf(stderr,"Pos = %g %g %g\n",PhydroPos(0)[0],PhydroPos(0)[1],PhydroPos(0)[2]);

    int Neighbors[MaxNeighborSize];
    int Ntotal;
    // Gather
    double tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=0;i<Pall.Nhydro;i++){
            Ntotal += GetNeighborsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            if(i == 0){
                int count = 0;
                for(int l=0;l<Ntotal;l++){
                    int index =  Neighbors[l];
                    if(DISTANCE2(PhydroPosP(0),PhydroPosP(index))<SQ(2.e0*Phydro[i]->Kernel))
                        count++;
                    fprintf(stderr,"%d ",index);
                }
                fprintf(stderr,"\n");
                fprintf(stderr,"NB = %d %d, Kernel size = %g %g\n",
                        Ntotal,count,Phydro[0]->Kernel,Phydro[0]->KernelPred);
            }
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-G- [%02d] Nmean = %ld, Ntotal = %d, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,GetElapsedTime()-tstart_each);
    }
    double tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-G- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);

    // Direct Gather
    tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=0;i<Pall.Nhydro;i++){
            //Ntotal += GetNeighborsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            Ntotal += GetNeighborsDirect(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            //int GetOnlyNumberofNeighbors(double Pos[restrict], const doubleh);
        }
        fprintf(stderr,"-GD- [%02d] Nmean = %ld, Ntotal = %d, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,GetElapsedTime()-tstart_each);
    }
    tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-GD- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);


    // Gather and Scatter
    tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=0;i<Pall.Nhydro;i++){
            Ntotal += GetNeighborsPairsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-GS- [%02d] Nmean = %ld, Ntotal = %d, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,GetElapsedTime()-tstart_each);
    }
    tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-GS- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);

    // Direct Gather and Scatter
    tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=0;i<Pall.Nhydro;i++){
            //Ntotal += GetNeighborsPairsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            Ntotal += GetNeighborsPairsDirect(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-GSD- [%02d] Nmean = %ld, Ntotal = %d, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,GetElapsedTime()-tstart_each);
    }
    tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-GSD- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);
#endif

    return EXIT_SUCCESS;
}

int main_NeighborSearchTest2(void){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    //ReadCosmologicalData();

    //InitUniformSphereTest(60);
    //InitUniformSphereTest(40);
    //InitUniformSphereTest(22);
    InitUniformSphereTest(16);
    //InitUniformSphereTest(8);

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    // First force / hydro calculation.
    ClearHydroData();
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    InitializeRootForHydro();
    PlantHydroTree();

    ClearHydroData();
    CalcKernel();

    for(int i=0;i<Pall.Nhydro;i++){ // add fluctuations.
        Phydro[i]->Kernel = 
        Phydro[i]->KernelPred = Phydro[i]->Kernel*(1.0-0.5*gsl_ran_gaussian(RandomGenerator,1.0));
        //Phydro[i]->KernelPred = Phydro[i]->Kernel*(1.0-0.5*gsl_rng_uniform(RandomGenerator));
    }
    PlantHydroTree();

    fprintf(stderr,"Kernel size = %g %g\n",Phydro[0]->Kernel,Phydro[0]->KernelPred);
    fprintf(stderr,"Pos = %g %g %g\n",PhydroPos(0)[0],PhydroPos(0)[1],PhydroPos(0)[2]);

    int Neighbors[MaxNeighborSize];
    int Ntotal;
    unsigned long int NeighborTotal;
    // Gather
    double tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=NeighborTotal=0;i<Pall.Nhydro;i++){
            int Nlist = GetNeighborsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            Ntotal += Nlist;
            for(int l=0;l<Nlist;l++){
                NeighborTotal += Neighbors[l];
            }

#if 0
            if(i == 0){
                int count = 0;
                for(int l=0;l<Ntotal;l++){
                    int index =  Neighbors[l];
                    if(DISTANCE2(PhydroPosP(0),PhydroPosP(index))<SQ(2.e0*Phydro[i]->Kernel))
                        count++;
                    fprintf(stderr,"%d ",index);
                }
                fprintf(stderr,"\n");
                fprintf(stderr,"NB = %d %d, Kernel size = %g %g\n",
                        Ntotal,count,Phydro[0]->Kernel,Phydro[0]->KernelPred);
            }
#endif
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-G- [%02d] Nmean = %ld, Ntotal = %d, NNT = %ld, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,NeighborTotal,GetElapsedTime()-tstart_each);
    }
    double tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-G- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);

    // Direct Gather
    tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=NeighborTotal=0;i<Pall.Nhydro;i++){
            //Ntotal += GetNeighborsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            int Nlist = GetNeighborsDirect(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            Ntotal += Nlist;
            for(int l=0;l<Nlist;l++){
                NeighborTotal += Neighbors[l];
            }
//int GetOnlyNumberofNeighbors(double Pos[restrict], const double h){
        }
        fprintf(stderr,"-GD- [%02d] Nmean = %ld, Ntotal = %d, NNT = %ld, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,NeighborTotal,GetElapsedTime()-tstart_each);
    }
    tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-GD- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);


    // Gather and Scatter
    tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=NeighborTotal=0;i<Pall.Nhydro;i++){
            int Nlist = GetNeighborsPairsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            Ntotal += Nlist;
            for(int l=0;l<Nlist;l++){
                NeighborTotal += Neighbors[l];
            }
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-GS- [%02d] Nmean = %ld, Ntotal = %d, NNT = %ld, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,NeighborTotal,GetElapsedTime()-tstart_each);
    }
    tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-GS- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);

    // Direct Gather and Scatter
    tstart = GetElapsedTime();
    for(int k=0;k<Iter;k++){
        double tstart_each = GetElapsedTime();
        for(int i=Ntotal=NeighborTotal=0;i<Pall.Nhydro;i++){
            //Ntotal += GetNeighborsPairsLimited(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            int Nlist = GetNeighborsPairsDirect(PhydroPosP(i),2.e0*Phydro[i]->Kernel,Neighbors);
            Ntotal += Nlist;
            for(int l=0;l<Nlist;l++){
                NeighborTotal += Neighbors[l];
            }
        }
        //PlantTreeFullUpdate(TreeForHydro);
        fprintf(stderr,"-GSD- [%02d] Nmean = %ld, Ntotal = %d, NNT = %ld, Time = %g\n",k,
                Ntotal/Pall.Nhydro,Ntotal,NeighborTotal,GetElapsedTime()-tstart_each);
    }
    tend = GetElapsedTime()-tstart;
    fprintf(stderr,"-GSD- [all] Nmean = %ld, Ntotal = %d, Time = %g\n",Ntotal/Pall.Nhydro,Ntotal,tend);

    return EXIT_SUCCESS;
}

int main_ForceHydroTest(void){

    // if NX = 11, totally 484 particles generated
    // if NX = 16, totally 1736 particles generated
    // if NX = 22, totally 4776 particles generated
    // if NX = 40, totally 30976 particles generated
    // if NX = 45, totally 44413 particles generated
    // if NX = 50, totally 61432 particles generated
    // if NX = 60, totally 107392 particles generated
    // if NX = 80, totally 257776 particles generated
    // if NX = 100, totally 507376 particles generated
    // if NX = 128, totally 1072184 particles generated

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();

    //InitColdCollapseTest(40*40*40,0.1,OFF);
    //InitRandomSphereTest(128*128*128);
    //InitUniformSphereTest(256);
    //InitUniformSphereTest(128);
    //InitUniformSphereTest(22);
    //InitUniformSphereTest(10);
    //Init3DCollapseTest(128);
    //Init3DCollapseTest(40);
    //Init3DCollapseTest(22);
    //InitRandomBoxTest(100000);

    //InitRandomBoxTestHydro(10000);

    ReadOldCosmologicalDataFull(10);

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    // Decomposition 
    InitializeDecomposition();
    DomainDecomposition();

    // First force / hydro calculation.
    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

#define IT  1
    PlantGravityTree();

// Check Acc and Pot.
#if 0
    FILE *fp_check;
    char fname_check[MaxCharactersInLine];

    sprintf(fname_check,"CheckGrav.%08d.%02d.%02d",Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check,fname_check,"w");
    for(int i=0;i<Pall.Ntotal;i++)
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = Pbody[i]->Pot = 0.e0;
    ForceTreeGRAPE();
    for(int i=0;i<Pall.Ntotal;i++)
        fprintf(fp_check,"%ld %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],Pbody[i]->Pot);
    fclose(fp_check);

    sprintf(fname_check,"CheckGrav_.%08d.%02d.%02d",Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check,fname_check,"w");
    for(int i=0;i<Pall.Ntotal;i++)
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = Pbody[i]->Pot = 0.e0;
    ForceTreeGRAPE_();
    for(int i=0;i<Pall.Ntotal;i++)
        fprintf(fp_check,"%ld %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],Pbody[i]->Pot);
    fclose(fp_check);

    MPI_Finalize();
    exit(-1098);
#endif


    InitializeRootForHydro();
    PlantHydroTree();

    ClearHydroData();
    CalcKernel();

// Check Rho Div Rot.
#if 0
    FILE *fp_check_hydro_rho;
    char fname_check_hydro_rho[MaxCharactersInLine];

    sprintf(fname_check_hydro_rho,"CheckHydroRho_.%08ld.%02d.%02d",Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_hydro_rho,fname_check_hydro_rho,"w");
    ClearHydroData();
    CalcDensityDivRotOmega_();
    for(int i=0;i<Pall.Nhydro;i++)
        fprintf(fp_check_hydro_rho,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Phydro[i]->Rho,Phydro[i]->Div,
                Phydro[i]->Rot[0],Phydro[i]->Rot[1],Phydro[i]->Rot[2],Phydro[i]->Omega);
    fclose(fp_check_hydro_rho);

    sprintf(fname_check_hydro_rho,"CheckHydroRho_ExportFlag_IsendIrecv.%08ld.%02d.%02d",Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_hydro_rho,fname_check_hydro_rho,"w");
    ClearHydroData();
    CalcDensityDivRotOmega_ExportFlag_IsendIrecv_WithRenew_DataUpdate();
    for(int i=0;i<Pall.Nhydro;i++)
        fprintf(fp_check_hydro_rho,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Phydro[i]->Rho,Phydro[i]->Div,
                Phydro[i]->Rot[0],Phydro[i]->Rot[1],Phydro[i]->Rot[2],Phydro[i]->Omega);
    fclose(fp_check_hydro_rho);

    /*
    sprintf(fname_check_hydro_rho,"CheckHydroRho_IsendIrecvPair.%08ld.%02d.%02d",Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_hydro_rho,fname_check_hydro_rho,"w");
    ClearHydroData();
    CalcDensityDivRotOmega_OmitSmallActiveJ(100);
    for(int i=0;i<Pall.Nhydro;i++)
        fprintf(fp_check_hydro_rho,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Phydro[i]->Rho,Phydro[i]->Div,
                Phydro[i]->Rot[0],Phydro[i]->Rot[1],Phydro[i]->Rot[2],Phydro[i]->Omega);
    fclose(fp_check_hydro_rho);
    */
    MPI_Finalize();
    exit(-1098);
#endif

//Rho Div Rot Omega Timing Result
#if 1
    {
    FILE *fp_check_hydro_rho_time;
    char fname_check_hydro_rho_time[MaxCharactersInLine];

    if(MPIGetMyID() == 0){
        sprintf(fname_check_hydro_rho_time,"TimeCheckHydroRho.%08ld.%02d",Pall.Ntotal_t,MPIGetNumProcs()); 
        FileOpen(fp_check_hydro_rho_time,fname_check_hydro_rho_time,"w");
    }

    TimingResults.HydroNeighborSearchThisStep = 0.e0;
    double TimeComm = GetElapsedTime();
    for(int i=0;i<IT;i++){
        ClearHydroData();
        CalcDensityDivRot();
    }
    if(MPIGetMyID() == 0){
        fprintf(fp_check_hydro_rho_time,"%ld %g %g\n",
                Pall.Ntotal_t,(GetElapsedTime()-TimeComm)/IT,TimingResults.HydroNeighborSearchThisStep/IT);
        fprintf(stderr,"1 : %ld %g %g\n",
                Pall.Ntotal_t,(GetElapsedTime()-TimeComm)/IT,TimingResults.HydroNeighborSearchThisStep/IT);
    }

    TimingResults.HydroNeighborSearchThisStep = 0.e0; 
    TimeComm = GetElapsedTime();
    for(int i=0;i<IT;i++){
        ClearHydroData();
        CalcDensityDivRot();
    }
    if(MPIGetMyID() == 0){
        fprintf(fp_check_hydro_rho_time,"%ld %g %g\n",
                Pall.Ntotal_t,(GetElapsedTime()-TimeComm)/IT,TimingResults.HydroNeighborSearchThisStep/IT);
        fprintf(stderr,"2 : %ld %g %g\n",
                Pall.Ntotal_t,(GetElapsedTime()-TimeComm)/IT,TimingResults.HydroNeighborSearchThisStep/IT);
    }
    if(MPIGetMyID() == 0)
        fclose(fp_check_hydro_rho_time);
    }
#endif


    ClearHydroData();
    CalcDensityDivRot();

// Check DuDt Acc.
#if 0
    FILE *fp_check_hydro_acc;
    char fname_check_hydro_acc[MaxCharactersInLine];

    sprintf(fname_check_hydro_acc,"CheckDuAcc.%08ld.%02d.%02d",Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_hydro_acc,fname_check_hydro_acc,"w");

    sprint("start CalcDuDtAcc_");
    CalcDuDtAcc_();
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp_check_hydro_acc,"%ld %g %g %g %g\n",Pbody[i]->GlobalID,
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
        Phydro[i]->Du = Phydro[i]->HydroAcc[0] = Phydro[i]->HydroAcc[1] = Phydro[i]->HydroAcc[2] = 0.e0;
    }
    fclose(fp_check_hydro_acc);

    sprintf(fname_check_hydro_acc,"CheckDuAcc_ExportFlag_IsendIrecv.%08ld.%02d.%02d",
        Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_hydro_acc,fname_check_hydro_acc,"w");

    sprint("start CalcDuDtAcc_ExportFlag_IsendIrecv");
    CalcDuDtAcc_ExportFlag_IsendIrecv();
    for(int i=0;i<Pall.Nhydro;i++)
        fprintf(fp_check_hydro_acc,"%ld %g %g %g %g\n",Pbody[i]->GlobalID,
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
    fclose(fp_check_hydro_acc);

    sprintf(fname_check_hydro_acc,"CheckDuAcc_ExportFlag_IsendIrecv_Group.%08ld.%02d.%02d",
        Pall.Ntotal_t,MPIGetNumProcs(),MPIGetMyID()); 
    FileOpen(fp_check_hydro_acc,fname_check_hydro_acc,"w");

    sprint("start CalcDuDtAcc_ExportFlag_IsendIrecv_Group");
    CalcDuDtAcc_ExportFlag_IsendIrecv_Group();
    for(int i=0;i<Pall.Nhydro;i++)
        fprintf(fp_check_hydro_acc,"%ld %g %g %g %g\n",Pbody[i]->GlobalID,
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
    fclose(fp_check_hydro_acc);


    MPI_Finalize();
    exit(-1098);
#endif

//Du Acc Timing Result
#if 1
    {
    FILE *fp_check_hydro_acc_time;
    char fname_check_hydro_acc_time[MaxCharactersInLine];

    if(MPIGetMyID() == 0){
        sprintf(fname_check_hydro_acc_time,"TimeCheckHydroAcc.%08ld.%02d",Pall.Ntotal_t,MPIGetNumProcs()); 
        FileOpen(fp_check_hydro_acc_time,fname_check_hydro_acc_time,"w");
    }

    PlantHydroTree();

    //
    TimingResults.HydroNeighborSearchThisStep = 0.e0;
    double TimeComm = GetElapsedTime();
    for(int i=0;i<IT;i++){
        CalcDuDtAcc();
    }
    if(MPIGetMyID() == 0)
        fprintf(fp_check_hydro_acc_time,"%ld %g %g\n",
            Pall.Ntotal_t,(GetElapsedTime()-TimeComm)/IT,TimingResults.HydroNeighborSearchThisStep/IT);

    //
    TimingResults.HydroNeighborSearchThisStep = 0.e0;
    TimeComm = GetElapsedTime();
    for(int i=0;i<IT;i++){
        CalcDuDtAcc();
    }
    if(MPIGetMyID() == 0)
        fprintf(fp_check_hydro_acc_time,"%ld %g %g\n",
            Pall.Ntotal_t,(GetElapsedTime()-TimeComm)/IT,TimingResults.HydroNeighborSearchThisStep/IT);

    if(MPIGetMyID() == 0)
        fclose(fp_check_hydro_acc_time);
    }
#endif

    return EXIT_SUCCESS;
}

int main_TreeTest(void){

    // if NX = 11, totally 484 particles generated
    // if NX = 16, totally 1736 particles generated
    // if NX = 22, totally 4776 particles generated
    // if NX = 40, totally 30976 particles generated
    // if NX = 45, totally 44413 particles generated
    // if NX = 50, totally 61432 particles generated
    // if NX = 60, totally 107392 particles generated
    // if NX = 80, totally 257776 particles generated
    // if NX = 100, totally 507376 particles generated
    // if NX = 128, totally 1072184 particles generated

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    //InitializeTreeTest(16);
    //InitColdCollapseTest(2*2*2*2*2*2*2*10*1024,0.1);
    //InitTreeTest(4*4*4*4*4*1024);
    InitTreeTest(4*1024);


    // First force / hydro calculation.
    BuildPredictors(); // Pos -> PosP/ Vel -> VelP


    // Decomposition 
    InitializeDecomposition();
    PreDomainDecomposition(0);
    DomainDecomposition();
    UpdateTotalActiveNumber();

    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();

    PlantGravityTree();
    for(int i = 0;i<10;i++){
        ClearTimingLogsThisStep();
        TimingResults.GravityThisStep = GetElapsedTime();
        ForceParallelTreeGRAPE();
        TimingResults.GravityThisStep = GetElapsedTime()-TimingResults.GravityThisStep;
        UpdateTimeLogs();
    }

    if(MPIGetMyID() == 0){
        FILE *fp;
        char fname[MaxCharactersInLine];
        sprintf(fname,"Force.log.%02d",MPIGetNumProcs());
        FileOpen(fp,fname,"a");
        fprintf(fp,"%ld %e %e\n",Pall.Ntotal_t,(TimingResults.Gravity-TimingResults.GravityComm)/10.0,TimingResults.GravityComm/10.0);
        fclose(fp);
    }

    return EXIT_SUCCESS;
}

