#include "config.h"
#include "Run.h"
#include "PreDecomposition.h"
#include "Decomposition.h"
#include "ForceMisc.h"
#include "ForceGRAPE.h"
#include "ForceParallelTreeGRAPE.h"
#include "PlantGravityTree.h"
#include "PlantHydroTree.h"
#include "Integral.h"
#include "TimeStep.h"
#include "HydroDensity.h"
#include "HydroAcc.h"
#include "HydroMisc.h"
#include "HydroKernel.h"
#include "SizeDetermination.h"
#include "NeighborSearch.h"
#include "TimeStep.h"
#include "Cooling.h"
#include "Heating.h"
#include "StarFormation.h"
#include "Delayed.h"
#include "HIIregion.h"
#include "SinkParticle.h"
#include "SetUpTestRun.h"
#include "Read.h"
#include "CommunicationTable.h"
#include "Logs.h"
#include "FOF.h"
#include "RunLogs.h"
#include "Cooling.h"
#include "StellarFeedback.h"
#include "StellarWind.h"
#include "FUV.h"
#include "ParticleSplitting.h"



#ifdef TASK_TEST_HYDRO_QUANTITIES //{


#   ifdef USE_NEIGHBOR_LIST //{

static int NeighborIndexCmp(const void *_x, const void *_y){
    return *(int*)_x - *(int*)_y;
}


static void CompareNeighbors(const int Index){

    for(int k=0;k<MPIGetNumProcs();k++){
        if(k==MPIGetMyID()){
            struct StructGetLocalNeighborList NB;
            NB = GetLocalNeighrborList(Index);
            fprintf(stderr,"L [%d]:%d: ",k,NB.Nlist); fflush(NULL);
            for(int i=0;i<NB.Nlist;i++)
                fprintf(stderr,"%d ",NB.Neighbors[i]);
            fprintf(stderr,"\n");
            fflush(NULL);

            
            int Target = 0;
            int RootNodeID = 0;
            int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
            int header = HydroNode[RootNodeID].Leaves;
            for(int l=0;l<NumberofLeaves;l++){
                int leaf = header + l;
                if(HydroRoot.Leaves[leaf] < 0) continue;
                if(NBCache[leaf].Active){
                    if(Target == Index) break;
                    Target ++;

                }
            }

            int leaf = NBCache[Target].Leaf;
            int Neighbors[MaxNeighborSize];
            int Nlist = GetNeighbors(PhydroBody(leaf)->PosP,
                    2*Phydro[leaf]->KernelPred,Neighbors);
            fprintf(stderr,"D [%d]:%d: ",
                    k,Nlist); fflush(NULL);
            for(int i=0;i<Nlist;i++)
                fprintf(stderr,"%d ",Neighbors[i]);
            fprintf(stderr,"\n");
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    return ;
}


static void CompareNlistSum(void){

    int Nlist = 0;
    int Neighbors[MaxNeighborSize];

    double Tnb = GetElapsedTime();
    int NActives = 0;
    long int SumNeighbors = 0; 
    int RootNodeID = 0;
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int l=0;l<NumberofLeaves;l++){
        int leaf = header + l;
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            int CurrentNlist = GetNeighbors(NBCache[leaf].Pos,2*NBCache[leaf].Kernel,Neighbors);
            Nlist += CurrentNlist;
            for(int k=0;k<CurrentNlist;k++){
                SumNeighbors += Neighbors[k];
            }
            NActives ++;

        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&Nlist,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Time[Old] = %g[sec], %d, %ld\n",GetElapsedTime()-Tnb,Nlist,SumNeighbors);
        fflush(NULL);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    SumNeighbors = 0;

    Nlist = 0; 
    Tnb = GetElapsedTime();
    int counter = 0;
    int counter_out = 0;
    for(int i=0;i<Pall.NActivesHydro;i++){
        struct StructGetLocalNeighborList NB;
        NB = GetLocalNeighrborList(i);
        if(NB.Nlist>0){
            Nlist += NB.Nlist;
            for(int k=0;k<NB.Nlist;k++){
                SumNeighbors += NB.Neighbors[k];
            }
            counter ++;
        } else {
            int CurrentNlist = GetNeighbors(NBCache[i].Pos,2*NBCache[i].Kernel,Neighbors);
            Nlist += CurrentNlist;
            for(int k=0;k<CurrentNlist;k++){
                SumNeighbors += Neighbors[k];
            }
            counter_out ++;
        }
    }
    // NB.Nlist != GetNeighbors();
    MPI_Allreduce(MPI_IN_PLACE,&Nlist,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Time[NBList] = %g[sec], %d, %ld\n",GetElapsedTime()-Tnb,Nlist,SumNeighbors);
        fflush(NULL);
    }
    fprintf(stderr,"[%d] counter %d, %d\n",MPIGetMyID(),counter,counter_out);

    MPI_Barrier(MPI_COMM_WORLD);
    return ;
}
#   endif // USE_NEIGHBOR_LIST //}

static void CompareNlistSumOriginal(void){

    int Nlist = 0;
    int Neighbors[MaxNeighborSize];

    double Tnb = GetElapsedTime();
    int NActives = 0;
    long int SumNeighbors = 0; 
    int RootNodeID = 0;
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int l=0;l<NumberofLeaves;l++){
        int leaf = header + l;
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            int CurrentNlist = GetNeighbors(NBCache[leaf].Pos,2*NBCache[leaf].Kernel,Neighbors);
            Nlist += CurrentNlist;
            for(int k=0;k<CurrentNlist;k++){
                SumNeighbors += Neighbors[k];
            }
            NActives ++;

        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&Nlist,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Time[Old] = %g[sec], %d, %ld\n",GetElapsedTime()-Tnb,Nlist,SumNeighbors);
        fflush(NULL);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    return ;
}

static void InitTestHydroQuantities(const int number, const int mode);
static void InitTestHydroQuantitiesRand(const int number);
int main_Test_HydroQuantities(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    if(Pall.RunStatus == NewSimulation){
        //InitTestHydroQuantities(20,0);
        //InitTestHydroQuantitiesRand(100000);
        InitTestHydroQuantitiesRand(50000);
    } else {
        fprintf(stderr,"You cannot use this function with the restart mode\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    InitLogFiles();
    ActivateAllparticles();

    InitializeRun();
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"***** Computation time = %g [sec]\n",
                TimingResults.HydroKernelThisStep);
        fflush(NULL);
    }
    ActivateAllparticles();

    PlantHydroTree();
    TimingResults.HydroKernelThisStep = 0.e0;

    CalcSize();

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"***** Computation time = %g [sec]\n",
                TimingResults.HydroKernelThisStep);
        fflush(NULL);
    }


    ClearHydroData();
    double trho = GetElapsedTime();
    CalcDensityDivRot();
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"***** Computation time/rho = %g [sec]\n",
                GetElapsedTime()-trho);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double tacc = GetElapsedTime();
    CalcDuDtAcc();
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"***** Computation time/acc = %g [sec]\n",
                GetElapsedTime()-tacc);
        fflush(NULL);
    }

    // CheckData //
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"./HydroQs.Result.%02d.%03d",NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        double  SmoothedNumber = Phydro[i]->SmoothedNumber;
#else  // USE_SMOOTHED_NEIGHBOR_NUMBER //}//{
        double  SmoothedNumber = Phydro[i]->Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef USE_GRAD_H //{
        double Gradh = Phydro[i]->Gradh;
#else // USE_GRAD_H //}//{
        double Gradh = 5.0;
#endif // USE_GRAD_H //}
#ifdef USE_GRAD_N //{
        double GradN = Phydro[i]->GradN;
#else // USE_GRAD_N //}//{
        double GradN = 6.0;
#endif // USE_GRAD_N //}

#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight = Phydro[i]->Mass*Phydro[i]->U;
        double Pressure = Pall.Gm1*Phydro[i]->EnergyDensity;
#elif defined(USE_SPSPH) //{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight = Phydro[i]->Zw;
        double Pressure = Pall.Gm1*Phydro[i]->Mass*(Phydro[i]->PseudoDensity/Phydro[i]->Zw)*Phydro[i]->U;
#else // USE_DISPH //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight = Phydro[i]->Mass;
        double Pressure = Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U;
#endif // USE_DISPH //}

        // double PressurePower = 2.0;
        double Entropy = 3.0;


#if 1
        fprintf(fp,"%ld %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->Tag,Phydro[i]->Nlist,SmoothedNumber,Phydro[i]->Mass,
                PhydroBody(i)->PosP[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->U,Smoothed,Weight,
                Pressure,Entropy,Gradh,GradN,
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
#else
#if 0
    int header = HydroNode[RootNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header + i;
        //if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            ActiveIndexList[NActives] = NBCache[leaf].Leaf;
            NActives ++;
        }
    }
#endif
        fprintf(fp,"%ld %g\n",PhydroBody(NBCache[i].Leaf)->GlobalID,NBCache[i].Kernel);

#if 0
        fprintf(fp,"%ld %d %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->Nlist,NORM(PhydroPos(i)),
                Phydro[i]->Rho,Phydro[i]->Kernel,Smoothed,Weight);
                // Pressure,Entropy,Gradh,GradN,
                //PressurePower,Entropy,Phydro[i]->EntropyPred,Phydro[i]->EntropyPredPower,
                // Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
#endif
#endif
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
#ifdef USE_DISPH //{
        sprintf(fname,"cat ./HydroQs.Result.%02d.* | sort -n > ./DHydroQs.Result.%02d",NProcs,NProcs);
#elif defined(USE_SPSPH) //}//{
        sprintf(fname,"cat ./HydroQs.Result.%02d.* | sort -n > ./YHydroQs.Result.%02d",NProcs,NProcs);
#else //}//{
        sprintf(fname,"cat ./HydroQs.Result.%02d.* | sort -n > ./SHydroQs.Result.%02d",NProcs,NProcs);
#endif // USE_DISPH //}
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./HydroQs.Result.%02d.*",NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    

    //////////////////////////////////////////////
    // Check Correction routine.

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->HydroAcc[0] =
        Phydro[i]->HydroAcc[1] =
        Phydro[i]->HydroAcc[2] =
        Phydro[i]->Du = 0.e0;
    }

    CalcDuDtAccEnergyDensityForCorrection();
    // CalcDensityDivRot();
    // CalcDuDtAcc();

    sprintf(fname,"./HydroQs.Result.Corr.%02d.%03d",NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        double  SmoothedNumber = Phydro[i]->SmoothedNumber;
#else  // USE_SMOOTHED_NEIGHBOR_NUMBER //}//{
        double  SmoothedNumber = Phydro[i]->Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef USE_GRAD_H //{
        double Gradh = Phydro[i]->Gradh;
#else // USE_GRAD_H //}//{
        double Gradh = 5.0;
#endif // USE_GRAD_H //}
#ifdef USE_GRAD_N //{
        double GradN = Phydro[i]->GradN;
#else // USE_GRAD_N //}//{
        double GradN = 6.0;
#endif // USE_GRAD_N //}

#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight = Phydro[i]->Mass*Phydro[i]->U;
        double Pressure = Pall.Gm1*Phydro[i]->EnergyDensity;
#elif defined(USE_SPSPH) //{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight = Phydro[i]->Zw;
        double Pressure = Pall.Gm1*Phydro[i]->Mass*(Phydro[i]->PseudoDensity/Phydro[i]->Zw)*Phydro[i]->U;
#else // USE_DISPH //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight = Phydro[i]->Mass;
        double Pressure = Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U;
#endif // USE_DISPH //}


        //double PressurePower = 2.0;
        double Entropy = 3.0;

        fprintf(fp,"%ld %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->Tag,Phydro[i]->Nlist,SmoothedNumber,Phydro[i]->Mass,
                //NORM(PhydroPos(i)),
                PhydroBody(i)->PosP[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                // PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->U,Smoothed,Weight,
                Pressure,Entropy,Gradh,GradN,
                //PressurePower,Entropy,Phydro[i]->EntropyPred,Phydro[i]->EntropyPredPower,
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
#ifdef USE_DISPH //{
        sprintf(fname,"cat ./HydroQs.Result.Corr.%02d.* | sort -n > ./DHydroQs.Result.Corr.%02d",NProcs,NProcs);
#elif defined(USE_SPSPH) //}//{
        sprintf(fname,"cat ./HydroQs.Result.Corr.%02d.* | sort -n > ./YHydroQs.Result.Corr.%02d",NProcs,NProcs);
#else //}//{
        sprintf(fname,"cat ./HydroQs.Result.Corr.%02d.* | sort -n > ./SHydroQs.Result.Corr.%02d",NProcs,NProcs);
#endif // USE_DISPH //}
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./HydroQs.Result.Corr.%02d.*",NProcs);
        system(fname);
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    // Check Correction routine.
    //////////////////////////////////////////////

#ifdef USE_NEIGHBOR_LIST //{
    //////////////////////////////////////////////
    // Check neighbor list.

    /// checker
    CompareNeighbors(0);
    CompareNeighbors(1);
    CompareNeighbors(3);
    CompareNeighbors(Pall.Nhydro-1);

    CompareNlistSum();

    // Check neighbor list.
    //////////////////////////////////////////////
#else  // USE_NEIGHBOR_LIST //}//{
    CompareNlistSumOriginal();
#endif // USE_NEIGHBOR_LIST //}


    fprintf(stderr,"Finish test!\n");
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    return EXIT_SUCCESS;
}

static void InitTestHydroQuantities(const int number, const int mode){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    InitializeRandomGenerator(1977);

    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    if(mode != 0){
                        Pos[0] /= r;
                        Pos[1] /= r;
                        Pos[2] /= r;
                        /*** strech grid(r_new = r_old**3) ***/
                        Pos[0] *= (r)*sqrt(r);
                        Pos[1] *= (r)*sqrt(r);
                        Pos[2] *= (r)*sqrt(r);
                    }
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
    SetViscosityParameters(0.1,1.0,1.0,0.1);

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    double Uinit = 0.05;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;

	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    if(mode != 0){
                        Pos[0] /= r;
                        Pos[1] /= r;
                        Pos[2] /= r;
				
                        /*** strech grid(r_new = r_old**3) ***/
                        Pos[0] *= (r)*sqrt(r);
                        Pos[1] *= (r)*sqrt(r);
                        Pos[2] *= (r)*sqrt(r);
                    }
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
                double Vel[3] = {
                    2.0*gsl_rng_uniform(RandomGenerator)-1.0,
                    2.0*gsl_rng_uniform(RandomGenerator)-1.0,
                    2.0*gsl_rng_uniform(RandomGenerator)-1.0};
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        /*
                        Pbody[mycount]->Vel[0] = TINY;
                        Pbody[mycount]->Vel[1] = TINY;
                        Pbody[mycount]->Vel[2] = TINY;
                        */
                        Pbody[mycount]->Vel[0] = Vel[0];
                        Pbody[mycount]->Vel[1] = Vel[1];
                        Pbody[mycount]->Vel[2] = Vel[2];

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        if(mode == 0){
                            PbodyHydro(mycount)->Rho = 3.0/(4*M_PI);
                        } else {
                            PbodyHydro(mycount)->Rho = (1.0/(2.0*M_PI))/(NORM(Pos)+eps);
                        }
                        PbodyHydro(mycount)->U = Uinit;
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    ActivateAllparticles();


    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);
    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();
    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"./HydroQs.Init.%02d.%03d",NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %d %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->Tag,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"cat ./HydroQs.Init.%02d.* | sort -n > ./HydroQs.Init.%02d",NProcs,NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./HydroQs.Init.%02d.*",NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 8;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 1.0;
    Pall.ViscousL = 0.1;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/10.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/HQ.ASCII");
    strcpy(Pall.BaseFileName,"./data/HQ");
    strcpy(Pall.RestartFileName,"./data/HQ.dump");


    return;
}

static void InitTestHydroQuantitiesRand(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    InitializeRandomGenerator(1977);

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
    SetViscosityParameters(0.1,1.0,1.0,0.1);


    int count = 0;
    int mycount = 0;
    for(int i=0;i<number;i++){
        if(count%NProcs == MyID){
            mycount ++;
        }
        count ++;
    }
    dprintlmpi(count);
    dprintlmpi(mycount);


    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)number;
    double eps = (cbrt((double)1736/(double)number));
    double Uinit = 0.05;

	count = 0; mycount = 0;
    for(int i=0;i<number;i++){
        double Pos[3] = {
            2.0*gsl_rng_uniform(RandomGenerator)-1.0,
            2.0*gsl_rng_uniform(RandomGenerator)-1.0,
            2.0*gsl_rng_uniform(RandomGenerator)-1.0};
        double Vel[3] = {
            2.0*gsl_rng_uniform(RandomGenerator)-1.0,
            2.0*gsl_rng_uniform(RandomGenerator)-1.0,
            2.0*gsl_rng_uniform(RandomGenerator)-1.0};

        if(count%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeHydro;
            Pbody[mycount]->GlobalID = count;

            Pbody[mycount]->Pos[0] = Pos[0];
            Pbody[mycount]->Pos[1] = Pos[1];
            Pbody[mycount]->Pos[2] = Pos[2];

            Pbody[mycount]->Vel[0] = Vel[0];
            Pbody[mycount]->Vel[1] = Vel[1];
            Pbody[mycount]->Vel[2] = Vel[2];

            Pbody[mycount]->Mass = mass;
            Pbody[mycount]->Eps = eps;

            PbodyHydro(mycount)->Use = ON;
            PbodyHydro(mycount)->Kernel = 0.1*eps;
            PbodyHydro(mycount)->Rho = 3.0/(4*M_PI);
            PbodyHydro(mycount)->U = Uinit;
            mycount ++;
        }
        count ++;
    }


    ActivateAllparticles();


    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);
    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();
    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"./HydroQs.Init.%02d.%03d",NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %d %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->Tag,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"cat ./HydroQs.Init.%02d.* | sort -n > ./HydroQs.Init.%02d",NProcs,NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./HydroQs.Init.%02d.*",NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 8;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 1.0;
    Pall.ViscousL = 0.1;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/10.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/HQ.ASCII");
    strcpy(Pall.BaseFileName,"./data/HQ");
    strcpy(Pall.RestartFileName,"./data/HQ.dump");


    return;
}
#endif // TASK_TEST_HYDRO_QUANTITIES //{

#ifdef TASK_TEST_GRAVITYTREE //{
static void WriteGravitationalForce(const int Tag);
static void InitTestGravityTree(const int number, const double Power);
int main_GravityTreeTest(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    // Make initial condition.
    InitTestGravityTree(128,1.5);
    //InitTestGravityTree(40,1.5);
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of particles is %ld\n",Pall.Ntotal_t);

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP
    InitializeDecomposition();
    DomainDecomposition();

    InitializeRootForGravity();
    InitializeRootForLET();
    InitializeParallelTreeGRAPE();


    PlantGravityTree();
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();

    // Write force.
    WriteGravitationalForce(1);

#if 0
    PlantGravityTreeOld();
    ClearGravitationalForce();
    ForceParallelTreeGRAPE();

    // Write force.
    WriteGravitationalForce(0);
#endif

    return EXIT_SUCCESS;
}

static void WriteGravitationalForce(const int Tag){

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"Force.%02d.%02d.%02d",Tag,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],
                Pbody[i]->Mass,Pbody[i]->Eps);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char conv[MaxCharactersInLine];
        Snprintf(conv,"cat ./Force.%02d.%02d.?? | sort -n > ./AllForce.%02d.%02d",
                Tag,MPIGetNumProcs(),Tag,MPIGetNumProcs());
        system(conv);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}

static void InitTestGravityTree(const int number, const double Power){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    InitializeRandomGenerator(1977+MPIGetMyID());
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                double Pos[] = {-1.e0+i*dx,-1.e0+j*dx,-1.e0+k*dx};

				double r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
                    
                    /*** strech grid(r_new = r_old**3) ***/
                    double r_fact = pow(r,Power);
                    Pos[0] *= r_fact;
                    Pos[1] *= r_fact;
                    Pos[2] *= r_fact;
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.NDM = Pall.Ntotal = mycount;
    Pall.NDM_t = Pall.Ntotal_t = count;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);

    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    //double eps = 0.01*(cbrt((double)1736/(double)count));
    //eps = 0.01;

    // make Hydro sphere. +5
	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                double Pos[] = {-1.e0+i*dx,-1.e0+j*dx,-1.e0+k*dx};

				double r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    double r_fact = pow(r,Power);
                    Pos[0] *= r_fact;
                    Pos[1] *= r_fact;
                    Pos[2] *= r_fact;
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeDM;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = TINY;
                        Pbody[mycount]->Vel[1] = TINY;
                        Pbody[mycount]->Vel[2] = TINY;

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    ReConnectPointers();
    UpdateTotalNumber();

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"GravityTreeTest.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],Pbody[i]->Mass);
    }
    fclose(fp);
#endif
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    return;
}

#endif //TASK_TEST_GRAVITYTREE //}

#ifdef TASK_TEST_HYDROTREE //{
static void InitTestHydroTree(const int number);
//static void WriteDataTestHydroTree(void);
int main_HydroTreeTest(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    if(Pall.RunStatus == NewSimulation){
        InitTestHydroTree(20);
    } else {
        fprintf(stderr,"You cannot use this function with the restart mode\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP
    InitializeDecomposition();
    DomainDecomposition();

    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.Nhydro);
    dlprintlmpi(Pall.NDM);

    InitializeRootForHydro();
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->Active = ON;
    }
    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesHydro_t = Pall.Nhydro_t;


    ClearHydroData();

    PlantHydroTree();
    CalcKernel();
    CalcDensityDivRot();
    CalcDuDtAcc();

    // Write data.
    for(int k=0;k<MPIGetNumProcs();k++){
        if(k == MPIGetMyID()){
            double Temperature;
            FILE *fp_hydro,*fp_star,*fp_dm;
            char fname_hydro[MaxCharactersInLine],fname_star[MaxCharactersInLine],fname_dm[MaxCharactersInLine];
            sprintf(fname_hydro,"%s.Hydro.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());
            sprintf(fname_dm,"%s.DM.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());

            if(Pall.Nhydro_t>0){
                FileOpen(fp_hydro,fname_hydro,"w");
            }
            if(Pall.NDM_t>0){
                FileOpen(fp_dm,fname_dm,"w");
            }

            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Type == TypeHydro){
                    if(Pall.Nhydro_t>0)
                    fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,PbodyHydro(i)->Rho,
                            PbodyHydroU(i),PbodyHydroKernel(i),
                            PbodyHydroAcc(i)[0],PbodyHydroAcc(i)[1],PbodyHydroAcc(i)[2],
                            PbodyHydro(i)->Nlist);
                } else if(Pbody[i]->Type == TypeDM){
                    if(Pall.NDM_t>0)
                    fprintf(fp_dm,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass);
                }
            }

            if(Pall.Nhydro_t>0){
                fclose(fp_hydro);
            }
            if(Pall.NDM_t>0){
                fclose(fp_dm);
            }
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    fprintf(stderr,"Finish test!\n");
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    return EXIT_SUCCESS;
}

int main_HydroTreeIntegralTest(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    if(Pall.RunStatus == NewSimulation){
        InitTestHydroTree(20);
    } else {
        fprintf(stderr,"You cannot use this function with the restart mode\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    InitLogFiles();

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    Run();
    OutPutAllParticlesInASCIIFormat();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    MPI_Barrier(MPI_COMM_WORLD);
    return EXIT_SUCCESS;
}

static void WriteHydroDataForHydroTest(const int AdditionalFlag){

    for(int k=0;k<MPIGetNumProcs();k++){
        if(k == MPIGetMyID()){
            FILE *fp_hydro;
            char fname_hydro[MaxCharactersInLine];
            sprintf(fname_hydro,"%s.Hydro.%02d.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),AdditionalFlag,MPIGetMyID());
            FileOpen(fp_hydro,fname_hydro,"w");

            for(int i=0;i<Pall.Nhydro;i++){
                    //fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",Pbody[i]->GlobalID,
                            //Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            //Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            //Pbody[i]->Mass,PbodyHydro(i)->Rho,
                            //PbodyHydroU(i),PbodyHydroKernel(i),
                            //PbodyHydroAcc(i)[0],PbodyHydroAcc(i)[1],PbodyHydroAcc(i)[2],
                            //PbodyHydro(i)->Nlist);
                    //fprintf(fp_hydro,"%ld %d %1.15le %1.15le %1.15le\n",PhydroBody(i)->GlobalID,
                    fprintf(fp_hydro,"%ld %d %1.7le %1.7le %1.7le\n",PhydroBody(i)->GlobalID,
                            Phydro[i]->Nlist,Phydro[i]->Rho,Phydro[i]->Div,Phydro[i]->Kernel);
                            //Phydro[i]->Nlist,Phydro[i]->Rho,Phydro[i]->Du,NORM(Phydro[i]->HydroAcc));
                            //PbodyHydro(i)->Nlist,PbodyHydro(i)->Rho,PbodyHydroU(i),PbodyHydroKernel(i));

            }

            fclose(fp_hydro);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return ;

}

static void InitHydroTreeRobustTest(const int number, const int power){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

    //InitializeRandomGenerator(1977+MPIGetMyID());
    InitializeRandomGenerator(1977);

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx; Pos[1] = -1.e0+j*dx; Pos[2] = -1.e0+k*dx;
				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    double ir = 1.0/r;
                    double sqrtr = sqrt(r);
                    Pos[0] *= ir; Pos[1] *= ir; Pos[2] *= ir;
                    /*** strech grid(r_new = r_old**3) ***/
                    double r_new = pow(sqrtr,(double)power);
                    Pos[0] *= r_new; Pos[1] *= r_new; Pos[2] *= r_new;

                    //Pos[0] *= r*sqrtr; Pos[1] *= r*sqrtr; Pos[2] *= r*sqrtr;
                } else {
                    Pos[0] = Pos[1] = Pos[2] = 0.e0;
                }
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    double Uinit = 0.05;
    double InitVelocity = 1.0;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

                double rand[3] = {gsl_rng_uniform(RandomGenerator),
                                  gsl_rng_uniform(RandomGenerator),
                                  gsl_rng_uniform(RandomGenerator)};

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    double ir = 1.0/r;
                    double sqrtr = sqrt(r);
                    Pos[0] *= ir; Pos[1] *= ir; Pos[2] *= ir;
                    /*** strech grid(r_new = r_old**3) ***/
                    double r_new = pow(sqrtr,(double)power);
                    Pos[0] *= r_new; Pos[1] *= r_new; Pos[2] *= r_new;
                    //Pos[0] *= r*sqrtr; Pos[1] *= r*sqrtr; Pos[2] *= r*sqrtr;
                } else {
                    Pos[0] = Pos[1] = Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = InitVelocity*(2.0*rand[0]-1.0);
                        Pbody[mycount]->Vel[1] = InitVelocity*(2.0*rand[1]-1.0);
                        Pbody[mycount]->Vel[2] = InitVelocity*(2.0*rand[2]-1.0);
                        //fprintf(stderr,"%g %g %g\n",Pbody[mycount]->Vel[0],Pbody[mycount]->Vel[1],Pbody[mycount]->Vel[2]);

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        PbodyHydro(mycount)->U = Uinit;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    ReConnectPointers();
    UpdateTotalNumber();

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"3Dcollapse.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 100.0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 3.0;
    //Pall.TEnd = 10.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/10.0; 
    //Pall.OutPutInterval = Pall.TEnd/1.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/3D.ASCII");
    strcpy(Pall.BaseFileName,"./data/3D");
    strcpy(Pall.RestartFileName,"./data/3D.dump");

    return;
}

/*
 * This routine checkes the robustness of the hydro tree.  First this makes a
 * hydro particle distribution and calc hydrodynamical values. Then vanishes
 * several hydro particles.  After that check the hydrodynamical values again
 * and compare the values with the hydrodynamical values before the particle
 * vanishement.
 */
int main_HydroTreeRobustTest(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    if(Pall.RunStatus == NewSimulation){
        //InitTestHydroTree(20);
        //Init3DCollapseTest(12);
        InitHydroTreeRobustTest(12,3);
    } else {
        fprintf(stderr,"You cannot use this function with the restart mode\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    InitLogFiles();

    /*
    */
    for(int i=0;i<Pall.Ntotal;i++) Pbody[i]->Active = ON;
    for(int i=0;i<Pall.Nhydro;i++) Phydro[i]->Active = ON;
    Pall.NActives = Pall.Ntotal;        Pall.NActives_t = Pall.Ntotal_t;
    Pall.NActivesHydro = Pall.Nhydro;   Pall.NActivesHydro_t = Pall.Nhydro_t;

    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }
    int nactives = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->GlobalID%2 == 0){
            Pbody[i]->Active = ON;
            PbodyHydro(i)->Active = ON;
            nactives++;
        } else {
            Pbody[i]->Active = OFF;
            PbodyHydro(i)->Active = OFF;
        }
    }
    Pall.NActives =  Pall.NActives_t = Pall.NActivesHydro = Pall.NActivesHydro_t = nactives;
    dprintlmpi(nactives);

    for(int i=0;i<Pall.Nhydro_t;i++)
        Phydro[i]->Nlist = 0;


    int nactives_t = 0; int nactives_h = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active == ON)      nactives_t ++;
        if(PbodyHydro(i)->Active == ON) nactives_h ++;
    }
    Pall.NActives =  Pall.NActives_t = nactives_t;
    Pall.NActivesHydro = Pall.NActivesHydro_t = nactives_h;

    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active == ON){
            Phydro[i]->Nlist = 0;
        }
    }

    gprintlmpi(Phydro[0]->VelP[0]);
    gprintlmpi(Phydro[0]->PosP[0]);
    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroBody(i)->GlobalID == 116){
            fprintf(stderr,"116 %1.15g %1.15g %1.15g %1.15g\n",
                    Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2],Phydro[i]->KernelPred);
        }
        if(PhydroBody(i)->GlobalID == 286){
            fprintf(stderr,"286 %1.15g %1.15g %1.15g %1.15g\n",
                    Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2],Phydro[i]->KernelPred);
        }
    }

    PlantHydroTree();
    ClearHydroData();
    CalcDensityDivRot();
    CalcDuDtAcc();
    WriteHydroDataForHydroTest(0);
    char conv[MaxCharactersInLine];
    Snprintf(conv,"cat ./data/3D.ASCII.Hydro.%02d.00.0? | sort -n > ./check/3D.ASCII.Hydro.%02d.00",MPIGetNumProcs(),MPIGetNumProcs());
    system(conv);


    // search maximum density.
    double rho = 0.0;
    for(int i=0;i<Pall.Nhydro;i++){
        rho = fmax(rho,Phydro[i]->Rho);
    }
    double wrho;
    MPI_Allreduce(&rho,&wrho,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    wrho *=0.19;
    int count = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Rho > wrho)
            count++;
    }
    int wcount;
    MPI_Allreduce(&count,&wcount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    dprintlmpi(wcount);
    
    // Sink several particles.
    Pall.SinkThresholdDensity = wrho;
    SinkParticles();

    for(int i=0;i<HydroRoot.NumberofLeaves;i++){
        if(HydroRoot.Leaves[i] < 0){
            dlprintlmpi(PhydroBody(-1*HydroRoot.Leaves[i])->GlobalID);
            dprintlmpi(-1*HydroRoot.Leaves[i]);
        }
    }
    for(int i=0;i<Pall.Nsink;i++){
        dlprintlmpi(PsinkBody(i)->GlobalID);
    }
    nactives_t = 0; nactives_h = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active == ON)      nactives_t ++;
        if(PbodyHydro(i)->Active == ON) nactives_h ++;
    }
    Pall.NActives =  Pall.NActives_t = nactives_t;
    Pall.NActivesHydro = Pall.NActivesHydro_t = nactives_h;

    ClearHydroData();
    PlantHydroTreeUpdate();
    CalcDensityDivRot();
    CalcDuDtAcc();
    WriteHydroDataForHydroTest(1);
    Snprintf(conv,"cat ./data/3D.ASCII.Hydro.%02d.01.0? | sort -n > ./check/3D.ASCII.Hydro.%02d.01",MPIGetNumProcs(),MPIGetNumProcs());
    system(conv);


    // Add 0.5*Kernel dispersion for particle positions.
#if 1
    for(int i=0;i<Pall.Nhydro;i++){
        /*
        InitializeRandomGenerator(PhydroBody(i)->GlobalID);
        PhydroBody(i)->Pos[0] = PhydroBody(i)->PosP[0] = Phydro[i]->PosP[0] = 
            + gsl_ran_gaussian(RandomGenerator,Phydro[i]->Kernel);
            //+ gsl_ran_gaussian(RandomGenerator,Phydro[i]->Kernel);
        PhydroBody(i)->Pos[1] = PhydroBody(i)->PosP[1] = Phydro[i]->PosP[1] = 
            + gsl_ran_gaussian(RandomGenerator,Phydro[i]->Kernel);
        PhydroBody(i)->Pos[2] = PhydroBody(i)->PosP[2] = Phydro[i]->PosP[2] = 
            + gsl_ran_gaussian(RandomGenerator,Phydro[i]->Kernel);
        */
        PhydroBody(i)->Pos[0] = PhydroBody(i)->PosP[0] = Phydro[i]->PosP[0] += 1.0;
        PhydroBody(i)->Pos[1] = PhydroBody(i)->PosP[1] = Phydro[i]->PosP[1] += 1.0;
        PhydroBody(i)->Pos[2] = PhydroBody(i)->PosP[2] = Phydro[i]->PosP[2] += 1.0;
    }
    for(int i=0;i<Pall.Nhydro;i++)
        if(PhydroBody(i)->GlobalID == 0)
            fprintf(stderr,"======= %g %g %g\n",Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2]);

    ClearHydroData();
    PlantHydroTreeUpdate();
    CalcKernel();
    CalcDensityDivRot();
    CalcDuDtAcc();
    WriteHydroDataForHydroTest(3);
    Snprintf(conv,"cat ./data/3D.ASCII.Hydro.%02d.03.0? | sort -n > ./check/3D.ASCII.Hydro.%02d.03",MPIGetNumProcs(),MPIGetNumProcs());
    system(conv);

    for(int i=0;i<Pall.Nhydro;i++)
        if(PhydroBody(i)->GlobalID == 0)
            fprintf(stderr,"======= %g %g %g\n",Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2]);

#endif

    ClearHydroData();
    PlantHydroTree();
    CalcDensityDivRot();
    CalcDuDtAcc();
    WriteHydroDataForHydroTest(2);
    Snprintf(conv,"cat ./data/3D.ASCII.Hydro.%02d.02.0? | sort -n > ./check/3D.ASCII.Hydro.%02d.02",MPIGetNumProcs(),MPIGetNumProcs());
    system(conv);

    fflush(NULL);

    CloseLogFiles();

    MPI_Barrier(MPI_COMM_WORLD);
    return EXIT_SUCCESS;
}

#define Offset +5
static void InitTestHydroTree(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
                    
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = mycount*3;
    Pall.Ntotal_t = count*3;
    Pall.Nhydro = mycount;
    Pall.Nhydro_t = count;
    Pall.NDM = mycount*2;
    Pall.NDM_t = count*2;

    int AllocationSize = mycount; 
    GenerateStructPbody(3*AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    //double eps = 0.01*(cbrt((double)1736/(double)count));
    //eps = 0.01;
    double Uinit = 0.05;


//////////////////////////////////// OK 

    // make Hydro sphere. +5
	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0]+Offset;
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = TINY;
                        Pbody[mycount]->Vel[1] = TINY;
                        Pbody[mycount]->Vel[2] = TINY;

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        PbodyHydro(mycount)->U = Uinit;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    // make Nbody sphere. -5
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount+Pall.Nhydro]->Active = ON;
                        Pbody[mycount+Pall.Nhydro]->Use = ON;
                        Pbody[mycount+Pall.Nhydro]->Type = TypeDM;
                        Pbody[mycount+Pall.Nhydro]->GlobalID = count+Pall.Nhydro;

                        Pbody[mycount+Pall.Nhydro]->Pos[0] = Pos[0]-Offset;
                        Pbody[mycount+Pall.Nhydro]->Pos[1] = Pos[1];
                        Pbody[mycount+Pall.Nhydro]->Pos[2] = Pos[2];

                        Pbody[mycount+Pall.Nhydro]->Vel[0] = TINY;
                        Pbody[mycount+Pall.Nhydro]->Vel[1] = TINY;
                        Pbody[mycount+Pall.Nhydro]->Vel[2] = TINY;

                        Pbody[mycount+Pall.Nhydro]->Mass = mass;
                        Pbody[mycount+Pall.Nhydro]->Eps = eps;
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    // make Nbody sphere. -5
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount+Pall.Nhydro*2]->Active = ON;
                        Pbody[mycount+Pall.Nhydro*2]->Use = ON;
                        Pbody[mycount+Pall.Nhydro*2]->Type = TypeDM;
                        Pbody[mycount+Pall.Nhydro*2]->GlobalID = count+Pall.Nhydro*2;

                        Pbody[mycount+Pall.Nhydro*2]->Pos[0] = Pos[0]-Offset;
                        Pbody[mycount+Pall.Nhydro*2]->Pos[1] = Pos[1];
                        Pbody[mycount+Pall.Nhydro*2]->Pos[2] = Pos[2];

                        Pbody[mycount+Pall.Nhydro*2]->Vel[0] = TINY;
                        Pbody[mycount+Pall.Nhydro*2]->Vel[1] = TINY;
                        Pbody[mycount+Pall.Nhydro*2]->Vel[2] = TINY;

                        Pbody[mycount+Pall.Nhydro*2]->Mass = mass;
                        Pbody[mycount+Pall.Nhydro*2]->Eps = eps;
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    ReConnectPointers();
    UpdateTotalNumber();

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"HydtoTreeTest.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 100.0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/10.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/HT.ASCII");
    strcpy(Pall.BaseFileName,"./data/HT");
    strcpy(Pall.RestartFileName,"./data/HT.dump");
    InitializeRandomGenerator(1977+MPIGetMyID());

    return;
}

#endif // TASK_TEST_HYDROTREE //}

#ifdef TASK_TEST_NEIGHBORSEARCH //{
static void Write_Test_NeighborSearch(char KeyWord[]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    ////////////////////////////////////////////////////////////////////////////////////
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"./%s.%02d.%03d",KeyWord,NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                Phydro[i]->Kernel);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"cat ./%s.%02d.* | sort -n > ./%s.%02d",
                KeyWord,NProcs,KeyWord,NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./%s.%02d.*",KeyWord,NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ////////////////////////////////////////////////////////////////////////////////////

}
    
static void InitTestNeighborSearch(const int number);
int main_Test_NeighborSearch(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    if(Pall.RunStatus == NewSimulation){
        InitTestNeighborSearch(80);
    } else {
        fprintf(stderr,"You cannot use this function with the restart mode\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP
    InitializeDecomposition();
    DomainDecomposition();

    InitializeRootForHydro();
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->Active = ON;
    }
    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesHydro_t = Pall.Nhydro_t;

    Write_Test_NeighborSearch("NeighborTest.Init");

    ClearHydroData();

    PlantHydroTree();
    double t = GetElapsedTime();
    CalcKernel();
    CalcKernel();
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Computation time = %g [sec]\n",GetElapsedTime()-t);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    Write_Test_NeighborSearch("NeighborTest.Fin");


    int Nlist = 0; 
    int Neighbors[MaxNeighborSize];
    double Tnb = GetElapsedTime();
    for(int i=0;i<Pall.Nhydro;i++){
        Nlist += GetNeighbors(PhydroBody(i)->Pos,2*Phydro[i]->Kernel,Neighbors);
        //Nlist += ReturnNeighborNumber(0,PhydroBody(i)->Pos,2*Phydro[i]->Kernel);
    }
    MPI_Allreduce(MPI_IN_PLACE,&Nlist,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Time[Old] = %g[sec], %d\n",GetElapsedTime()-Tnb,Nlist);
        fflush(NULL);
    }

    Nlist = 0; 
    Tnb = GetElapsedTime();
    for(int i=0;i<Pall.Nhydro;i++){
        Nlist += ReturnNeighborNumber(0,PhydroBody(i)->Pos,2*Phydro[i]->Kernel);
        //Nlist += GetNeighbors(PhydroBody(i)->Pos,2*Phydro[i]->Kernel,Neighbors);
    }
    MPI_Allreduce(MPI_IN_PLACE,&Nlist,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Time[New] = %g[sec], %d\n",GetElapsedTime()-Tnb,Nlist);
        fflush(NULL);
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Finish test!\n");
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return EXIT_SUCCESS;
}

static void InitTestNeighborSearch(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
                    
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    double Uinit = 0.05;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if((count%NProcs) == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

#if 0
                        int sign;
                        if((count%2)==0) sign = +1.0;
                        else             sign = -1.0;
#else
                        int sign = 1.0;
#endif

                        Pbody[mycount]->Vel[0] = +sign*Pos[0];
                        Pbody[mycount]->Vel[1] = -sign*Pos[1];
                        Pbody[mycount]->Vel[2] = +sign*Pos[2];

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        //PbodyHydro(mycount)->Kernel = 0.0361373;
                        PbodyHydro(mycount)->U = Uinit;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}


    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);
    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();
    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"HydtoTreeTest.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 100.0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/10.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/NT.ASCII");
    strcpy(Pall.BaseFileName,"./data/NT");
    strcpy(Pall.RestartFileName,"./data/NT.dump");
    InitializeRandomGenerator(1977+MPIGetMyID());

    return;
}

#endif //TASK_TEST_NEIGHBORSEARCH //}



#if 0
#ifdef TASK_TEST_HYDRO_QUANTITIES //{
static void InitTestHydroQuantities(const int number, const int mode);
int main_Test_HydroQuantities(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    if(Pall.RunStatus == NewSimulation){
        InitTestHydroQuantities(20,0);
    } else {
        fprintf(stderr,"You cannot use this function with the restart mode\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    InitLogFiles();

    InitializeRun();

    // CheckData //
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"./HydroQs.Result.%02d.%03d",NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        double  SmoothedNumber = Phydro[i]->SmoothedNumber;
#else  // USE_SMOOTHED_NEIGHBOR_NUMBER //}//{
        double  SmoothedNumber = Phydro[i]->Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef USE_GRAD_H //{
        double Gradh = Phydro[i]->Gradh;
#else // USE_GRAD_H //}//{
        double Gradh = 5.0;
#endif // USE_GRAD_H //}
#ifdef USE_GRAD_N //{
        double GradN = Phydro[i]->GradN;
#else // USE_GRAD_N //}//{
        double GradN = 6.0;
#endif // USE_GRAD_N //}

#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight = Phydro[i]->Mass*Phydro[i]->U;
        double Pressure = Pall.Gm1*Phydro[i]->EnergyDensity;
#elif defined(USE_SPSPH) //{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight = Phydro[i]->Zw;
        double Pressure = Pall.Gm1*Phydro[i]->Mass*(Phydro[i]->PseudoDensity/Phydro[i]->Zw)*Phydro[i]->U;
#else // USE_DISPH //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight = Phydro[i]->Mass;
        double Pressure = Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U;
#endif // USE_DISPH //}

        // double PressurePower = 2.0;
        double Entropy = 3.0;


        fprintf(fp,"%ld %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->Tag,Phydro[i]->Nlist,SmoothedNumber,Phydro[i]->Mass,
                //NORM(PhydroPos(i)),
                PhydroBody(i)->PosP[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                // PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->U,Smoothed,Weight,
                Pressure,Entropy,Gradh,GradN,
                //PressurePower,Entropy,Phydro[i]->EntropyPred,Phydro[i]->EntropyPredPower,
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
#ifdef USE_DISPH //{
        sprintf(fname,"cat ./HydroQs.Result.%02d.* | sort -n > ./DHydroQs.Result.%02d",NProcs,NProcs);
#elif defined(USE_SPSPH) //}//{
        sprintf(fname,"cat ./HydroQs.Result.%02d.* | sort -n > ./YHydroQs.Result.%02d",NProcs,NProcs);
#else //}//{
        sprintf(fname,"cat ./HydroQs.Result.%02d.* | sort -n > ./SHydroQs.Result.%02d",NProcs,NProcs);
#endif // USE_DISPH //}
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./HydroQs.Result.%02d.*",NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    

    //////////////////////////////////////////////
    // Check Correction routine.


    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->HydroAcc[0] =
        Phydro[i]->HydroAcc[1] =
        Phydro[i]->HydroAcc[2] =
        Phydro[i]->Du = 0.e0;
    }

    CalcDuDtAccEnergyDensityForCorrection();
    // CalcDensityDivRot();
    // CalcDuDtAcc();

    sprintf(fname,"./HydroQs.Result.Corr.%02d.%03d",NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        double  SmoothedNumber = Phydro[i]->SmoothedNumber;
#else  // USE_SMOOTHED_NEIGHBOR_NUMBER //}//{
        double  SmoothedNumber = Phydro[i]->Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef USE_GRAD_H //{
        double Gradh = Phydro[i]->Gradh;
#else // USE_GRAD_H //}//{
        double Gradh = 5.0;
#endif // USE_GRAD_H //}
#ifdef USE_GRAD_N //{
        double GradN = Phydro[i]->GradN;
#else // USE_GRAD_N //}//{
        double GradN = 6.0;
#endif // USE_GRAD_N //}

#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight = Phydro[i]->Mass*Phydro[i]->U;
        double Pressure = Pall.Gm1*Phydro[i]->EnergyDensity;
#elif defined(USE_SPSPH) //{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight = Phydro[i]->Zw;
        double Pressure = Pall.Gm1*Phydro[i]->Mass*(Phydro[i]->PseudoDensity/Phydro[i]->Zw)*Phydro[i]->U;
#else // USE_DISPH //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight = Phydro[i]->Mass;
        double Pressure = Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U;
#endif // USE_DISPH //}


        //double PressurePower = 2.0;
        double Entropy = 3.0;


        fprintf(fp,"%ld %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->Tag,Phydro[i]->Nlist,SmoothedNumber,Phydro[i]->Mass,
                //NORM(PhydroPos(i)),
                PhydroBody(i)->PosP[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                // PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->U,Smoothed,Weight,
                Pressure,Entropy,Gradh,GradN,
                //PressurePower,Entropy,Phydro[i]->EntropyPred,Phydro[i]->EntropyPredPower,
                Phydro[i]->Du,Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
#ifdef USE_DISPH //{
        sprintf(fname,"cat ./HydroQs.Result.Corr.%02d.* | sort -n > ./DHydroQs.Result.Corr.%02d",NProcs,NProcs);
#elif defined(USE_SPSPH) //}//{
        sprintf(fname,"cat ./HydroQs.Result.Corr.%02d.* | sort -n > ./YHydroQs.Result.Corr.%02d",NProcs,NProcs);
#else //}//{
        sprintf(fname,"cat ./HydroQs.Result.Corr.%02d.* | sort -n > ./SHydroQs.Result.Corr.%02d",NProcs,NProcs);
#endif // USE_DISPH //}
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./HydroQs.Result.Corr.%02d.*",NProcs);
        system(fname);
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    // Check Correction routine.
    //////////////////////////////////////////////


    fprintf(stderr,"Finish test!\n");
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    return EXIT_SUCCESS;
}

static void InitTestHydroQuantities(const int number, const int mode){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    if(mode != 0){
                        Pos[0] /= r;
                        Pos[1] /= r;
                        Pos[2] /= r;
                        /*** strech grid(r_new = r_old**3) ***/
                        Pos[0] *= (r)*sqrt(r);
                        Pos[1] *= (r)*sqrt(r);
                        Pos[2] *= (r)*sqrt(r);
                    }
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
    SetViscosityParameters(0.1,1.0,1.0,0.1);

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    double Uinit = 0.05;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    if(mode != 0){
                        Pos[0] /= r;
                        Pos[1] /= r;
                        Pos[2] /= r;
				
                        /*** strech grid(r_new = r_old**3) ***/
                        Pos[0] *= (r)*sqrt(r);
                        Pos[1] *= (r)*sqrt(r);
                        Pos[2] *= (r)*sqrt(r);
                    }
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = TINY;
                        Pbody[mycount]->Vel[1] = TINY;
                        Pbody[mycount]->Vel[2] = TINY;

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        if(mode == 0){
                            PbodyHydro(mycount)->Rho = 3.0/(4*M_PI);
                        } else {
                            PbodyHydro(mycount)->Rho = (1.0/(2.0*M_PI))/(NORM(Pos)+eps);
                        }
                        PbodyHydro(mycount)->U = Uinit;
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    ActivateAllparticles();

#if 0
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        PhydroBody(i)->Vel[0] = Phydro[i]->VelP[0] = PhydroBody(i)->Pos[0];
        PhydroBody(i)->Vel[1] = Phydro[i]->VelP[1] = PhydroBody(i)->Pos[1];
        PhydroBody(i)->Vel[2] = Phydro[i]->VelP[2] = PhydroBody(i)->Pos[2];

        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }
#endif

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);
    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();
    fprintf(stderr,"Pall. NActivesHydro = %d\n",Pall.NActivesHydro);

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"./HydroQs.Init.%02d.%03d",NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"cat ./HydroQs.Init.%02d.* | sort -n > ./HydroQs.Init.%02d",NProcs,NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./HydroQs.Init.%02d.*",NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 1.0;
    Pall.ViscousL = 0.1;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/10.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/HQ.ASCII");
    strcpy(Pall.BaseFileName,"./data/HQ");
    strcpy(Pall.RestartFileName,"./data/HQ.dump");
    InitializeRandomGenerator(1977+MPIGetMyID());


    return;
}
#endif // TASK_TEST_HYDRO_QUANTITIES //{
#endif

#ifdef TASK_TEST_SINKPARTICLE //{
#ifndef GRAVITY_RUN
#   error This test run requires GRVITY_RUN flag.
#endif
#ifndef HYDRO_RUN
#   error This test run requires HYDRO_RUN flag.
#endif
#ifndef ISOTHERMAL_EOS_RUN
#   error This test run requires ISOTHERMAL_EOS_RUN flag.
#endif
#ifndef USE_SINK_PARTICLE 
#   error This test run requires USE_SINK_PARTICLE flag.
#endif
#ifdef STARFORMATION
#   error Star formation routines are unnecessary in this run.
#endif
#ifdef DELAYED_FEEDBACK
#   error Feedback routines are unnecessary in this run. 
#endif


static void SetParticleDistribution(const int number);
static void ResetSinkThresholdDensity(void);

int main_Test_SinkParticleRun(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    SetParticleDistribution(18);

    // Run several steps.
    if(Pall.RunStatus == NewSimulation){
        InitializeRun();
    } else if (Pall.RunStatus == RestartSimulation){
        RestartRun();
    }

    // Set new threshold density.
    ResetSinkThresholdDensity();

    Run();

    // Check data consistency.

    return EXIT_SUCCESS;
}

static void ResetSinkThresholdDensity(void){

#if 0
    // Cloud mass
    double mass = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        mass += PhydroMass(i);
    }
    double wmass;
    MPI_Allreduce(&mass,&wmass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


    // Cloud size
    double max_distance = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        double r = NORM(PhydroPos(i));
        if(max_distance < r){
            max_distance = r;
        }
    }
    double wmax_distance;
    MPI_Allreduce(&max_distance,&wmax_distance,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    double mean_density = wmass/((4.0*M_PI/3.0)*CUBE(wmax_distance));

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"mean density = %g\n",mean_density);
        fprintf(stderr,"mean density = %g\n",Pall.ConvertNumberDenstityToCGS*mean_density);
    }
    // reset the sink density.
    Pall.SinkThresholdDensity = 1.1*mean_density;
#else
    // Search the maximum density.
    double max_density = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(max_density < Phydro[i]->Rho)
            max_density = Phydro[i]->Rho;
    }
    double GlobalMaxDensity;
    MPI_Allreduce(&max_density,&GlobalMaxDensity,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    // Reset the sink density.
    Pall.SinkThresholdDensity = 1000.0*GlobalMaxDensity;
    gprintlmpi(Pall.SinkThresholdDensity);
#endif


    return ;
}

static void SetParticleDistribution(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    // Set unit
    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDenstityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDenstityToCGS = GetUnitNumberDensityCGS();

    eprintlmpi(Pall.GravConst);
    eprintlmpi(Pall.DegreeOfFreedom);
    eprintlmpi(Pall.HeliumWeight);
    eprintlmpi(Pall.MeanMolecularWeight);
    eprintlmpi(Pall.ConvertUtoT);
    eprintlmpi(Pall.ConvertTtoU);

    eprintlmpi(Pall.ConvertDenstityToCGS);
    eprintlmpi(Pall.ConvertNumberDenstityToCGS);

    eprintlmpi(GetUnitVel());

    // Set Parameters
    // Rsphere, MCloud, cs, and dOmega/dt.
    double Rsphere = 5e+16/Pall.UnitLength; // 0.02pc
    double MCloud = 10.e0; // 1 Msun
    double cs = 0.1666*VELOCITY_KMS_CGS*GetUnitVel(); // 0.17 km/s
    //double omega = 7.2e-13*Pall.UnitTime; // 7.2x10^{-13} rad/s
    //double omega = 7.2e-13*GetUnitVel();
    double omega = 7.2e-13*Pall.UnitTime;

    eprintlmpi(cs);

    double Nhalf_real = (double)NX;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=-NX;i<NX;i++){
		for(int j=-NY;j<NY;j++){
			for(int k=-NZ;k<NZ;k++){
				Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
				Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
				Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

				if(Rsphere>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}


    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = MCloud/(double)count;
    //double eps = 0.1*(cbrt((double)1736/(double)count));
    double eps = 1.e-7*PC_CGS/Pall.UnitLength;
    //double Uinit = SQ(cs)/(5.0/3.0 * (5.0/3.0-1.0));

    double rho_crit = 7.e-24*SQ(Pall.Nhydro_t);
    //eps = cbrt((6*50*mass/(M_PI*rho_crit*GetUnitDensity())))/5.6;
    eps = cbrt((6*50*mass/(M_PI*rho_crit*GetUnitDensity())))/4.0;
    fprintf(stderr,"esp = %g [pc], critical density = %g [g/cm^3]\n",eps,rho_crit);
    //exit(0);

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;

    double v_tangential,vx,vy;
    double Rxy;
    double phi;
	for(int i=-NX;i<NX;i++){
		for(int j=-NY;j<NY;j++){
			for(int k=-NZ;k<NZ;k++){
				Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
				Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
				Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

				if(Rsphere>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Rxy = sqrt(SQ(Pos[0])+SQ(Pos[1]));

                        if(Rxy < TINY*Rsphere)
                            exit(444);

                        if(Pos[1] > 0.e0)
                            phi = acos(Pos[0]/Rxy);
                        else
                            phi = 2.0*PI-acos(Pos[0]/Rxy);

                        Pbody[mycount]->Mass = mass*(1.0+0.1*cos(2.0*phi));

                        //v_tangential = (Rxy * KM_PER_PC) * omega_in_rad_per_sec;
                        v_tangential = Rxy * omega;
                        vx=  -v_tangential * Pos[1]/Rxy;
                        vy=   v_tangential * Pos[0]/Rxy;

                        Pbody[mycount]->Vel[0] = vx;
                        Pbody[mycount]->Vel[1] = vy;
                        Pbody[mycount]->Vel[2] = 0;

                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        //PbodyHydro(mycount)->Kernel = 0.1*dx;
                        PbodyHydro(mycount)->Kernel = Rsphere/NX;
                        //eprintlmpi(PbodyHydro(mycount)->Kernel);

                        //cout<<"Internal energy..."<<endl;
                        // P = (gamma-1) rho u = cs_2 rho
                        // thus 2/3 u = kT/mp = cs_2
                        // thus  cs^2 = 3/2 kT/mp = u ; internal energy per mass

                        //for(i=0;i<N_gas;i++)
                            //buf1[i]=pow(cs_in_cm_per_sec, 2) * (UnitMass/UnitEnergy);

                        //PbodyHydro(mycount)->U = Uinit;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }
    ReConnectPointers();
    UpdateTotalNumber();

    // Check Total Mass

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"M2Collapse.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.CS = cs;

    //Pall.TEnd = 3.0; // ~1.3 dynamical time
    Pall.TEnd = 100.3*1.0774e+12/Pall.UnitTime;
    //Pall.TEnd = 3*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
            Pall.TEnd,Pall.TEnd*Pall.UnitTime);

    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 1.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    Pall.AdaptiveSofteningFactor = 1.e0;

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 100.0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    strcpy(Pall.ASCIIFileName,"./data/Sink.ASCII");
    strcpy(Pall.BaseFileName,"./data/Sink");
    strcpy(Pall.RestartFileName,"./data/Sink.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.05*1.e+12/Pall.UnitTime;

    //Pall.SinkThresholdDensity = 100000*Phydro[0]->Mass/CUBE(Phydro[0]->Kernel);

    OutPutASCIIDATA();

    FileOutPutConstantInterval();

    InitLogFiles();

    fflush(NULL);

    return;
}


static double MeanDistance;
static void MakeGridData(const int BaseNumber){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = BaseNumber;

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    // Set unit
    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDenstityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDenstityToCGS = GetUnitNumberDensityCGS();

    eprintlmpi(Pall.GravConst);
    eprintlmpi(Pall.DegreeOfFreedom);
    eprintlmpi(Pall.HeliumWeight);
    eprintlmpi(Pall.MeanMolecularWeight);
    eprintlmpi(Pall.ConvertUtoT);
    eprintlmpi(Pall.ConvertTtoU);

    eprintlmpi(Pall.ConvertDenstityToCGS);
    eprintlmpi(Pall.ConvertNumberDenstityToCGS);

    eprintlmpi(GetUnitVel());

    // Set Parameters
    // Rsphere, MCloud, cs, and dOmega/dt.
    double Rsphere = 5e+16/Pall.UnitLength; // 0.02pc
    double MCloud = 10.e0; // 1 Msun
    double cs = 0.1666*VELOCITY_KMS_CGS*GetUnitVel(); // 0.17 km/s
    //double omega = 7.2e-13*Pall.UnitTime; // 7.2x10^{-13} rad/s
    //double omega = 7.2e-13*GetUnitVel();
    double omega = 7.2e-13*Pall.UnitTime;

    eprintlmpi(cs);

    ///////////////// THIS IS THE FRONTLINE ////////////////////

    double Nhalf_real = (double)NX;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=-NX;i<NX;i++){
		for(int j=-NY;j<NY;j++){
			for(int k=-NZ;k<NZ;k++){
				Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
				Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
				Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

                if(count%NProcs == MyID){
                    mycount ++;
                }
                count ++;
			}
		}
	}
    MeanDistance = Rsphere/Nhalf_real;


    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = MCloud/(double)count;
    double eps = 1.e-7*PC_CGS/Pall.UnitLength;

    double rho_crit = 7.e-24*SQ(Pall.Nhydro_t);
    eps = cbrt((6*50*mass/(M_PI*rho_crit*GetUnitDensity())))/4.0;
    fprintf(stderr,"esp = %g [pc], critical density = %g [g/cm^3]\n",eps,rho_crit);

	dx = 2.e0/(double)(NX-1);
	count = mycount = 0;

    double v_tangential,vx,vy;
    double Rxy;
    double phi;
	for(int i=-NX;i<NX;i++){
		for(int j=-NY;j<NY;j++){
			for(int k=-NZ;k<NZ;k++){
				Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
				Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
				Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

                if(count%NProcs == MyID){
                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;
                    Pbody[mycount]->Type = TypeHydro;
                    Pbody[mycount]->GlobalID = count;

                    Pbody[mycount]->Pos[0] = Pos[0];
                    Pbody[mycount]->Pos[1] = Pos[1];
                    Pbody[mycount]->Pos[2] = Pos[2];

                    Rxy = sqrt(SQ(Pos[0])+SQ(Pos[1]));

                    if(Rxy < TINY*Rsphere)
                        exit(444);

                    if(Pos[1] > 0.e0)
                        phi = acos(Pos[0]/Rxy);
                    else
                        phi = 2.0*PI-acos(Pos[0]/Rxy);

                    Pbody[mycount]->Mass = mass*(1.0+0.1*cos(2.0*phi));

                    //v_tangential = (Rxy * KM_PER_PC) * omega_in_rad_per_sec;
                    v_tangential = Rxy * omega;
                    vx=  -v_tangential * Pos[1]/Rxy;
                    vy=   v_tangential * Pos[0]/Rxy;

                    Pbody[mycount]->Vel[0] = vx;
                    Pbody[mycount]->Vel[1] = vy;
                    Pbody[mycount]->Vel[2] = 0;

                    Pbody[mycount]->Eps = eps;

                    PbodyHydro(mycount)->Use = ON;
                    PbodyHydro(mycount)->Kernel = Rsphere/NX;

                    mycount ++;
                }
                count ++;
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }
    ReConnectPointers();
    UpdateTotalNumber();


#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"SinkGrid.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.CS = cs;

    Pall.TEnd = 10.3*1.0774e+12/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
            Pall.TEnd,Pall.TEnd*Pall.UnitTime);

    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 1.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    Pall.AdaptiveSofteningFactor = 1.e0;

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 100.0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    strcpy(Pall.ASCIIFileName,"./data/Sink.ASCII");
    strcpy(Pall.BaseFileName,"./data/Sink");
    strcpy(Pall.RestartFileName,"./data/Sink.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.05*1.e+12/Pall.UnitTime;

    Pall.SinkThresholdDensity = 100*Phydro[0]->Mass/Phydro[0]->Kernel;

    OutPutASCIIDATA();

    FileOutPutConstantInterval();
    InitLogFiles();
    fflush(NULL);
    return ;
}

struct StructSinkTarget{
    int index;
    double distance;
};

static int SinkDistCmp(const void *_x, const void *_y){
    const struct StructSinkTarget *pointer1 = (struct StructSinkTarget *)_x;
    const struct StructSinkTarget *pointer2 = (struct StructSinkTarget *)_y;
    if(pointer1->distance > pointer2->distance){
        return 1;
    } else if(pointer1->distance < pointer2->distance){
        return -1;
    } else {
        return 0;
    }
}

static void ConvertSinks(const int SinkNumber, const int BaseNumber){

#if 0
    struct StructSinkTarget *SinkTarget;
    SinkTarget = malloc(sizeof(struct StructSinkTarget)*Pall.Nhydro);
    for(int i=0;i<Pall.Nhydro;i++){
        SinkTarget[i].index = i;
        SinkTarget[i].distance = NORM(PhydroBody(i)->Pos);
    }
    qsort(SinkTarget,Pall.Nhydro,sizeof(struct StructSinkTarget),(int(*)(const void*, const void*))SinkDistCmp);

    // dx ?
    double dx = 2.e0/(double)(BaseNumber-1);

    int NewSinkParticles = 0;
    for(int k=0;k<SinkNumber;k++){
        int i = SinkTarget[k].index;
        fprintf(stderr,"Sink target = %ld\n",(PhydroBody(SinkTarget[k].index)->GlobalID));
        // do sink convert!
        StructPsinkptr Psk = ReturnEmptySinkStructurePointer();

        // Setup sink partcile.
        Psk->Use = ON;
        Psk->ParentGlobalID = PhydroBody(i)->GlobalID;
        Psk->Body = PhydroBody(i);
        Psk->Body->Baryon = (void*)Psk;
        Psk->Body->Type = TypeSink;
        Psk->FormationTime = Pall.TCurrent;

        Psk->NumberofAbsorbedParticles = 0;
        Psk->PosP[0] = Phydro[i]->PosP[0];
        Psk->PosP[1] = Phydro[i]->PosP[1];
        Psk->PosP[2] = Phydro[i]->PosP[2];

        Psk->VelP[0] = Phydro[i]->VelP[0];
        Psk->VelP[1] = Phydro[i]->VelP[1];
        Psk->VelP[2] = Phydro[i]->VelP[2];
        Psk->Radius = SINK_RADIUS;

        double dt_half_hydro = 0.5*Phydro[i]->dt_hydro;
        PhydroBody(i)->Vel[0] = PhydroBody(i)->Velh[0]+dt_half_hydro*Phydro[i]->HydroAcc[0];
        PhydroBody(i)->Vel[1] = PhydroBody(i)->Velh[1]+dt_half_hydro*Phydro[i]->HydroAcc[1];
        PhydroBody(i)->Vel[2] = PhydroBody(i)->Velh[2]+dt_half_hydro*Phydro[i]->HydroAcc[2];

        HydroRoot.Leaves[Phydro[i]->Leaf] *= -1;
        Phydro[i]->Use = OFF;

        NewSinkParticles ++;
    }
#elif 0
    int NewSinkParticles = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroBody(i)->GlobalID != 31179) continue;

        // do sink convert!
        StructPsinkptr Psk = ReturnEmptySinkStructurePointer();

        // Setup sink partcile.
        Psk->Use = ON;
        Psk->ParentGlobalID = PhydroBody(i)->GlobalID;
        Psk->Body = PhydroBody(i);
        Psk->Body->Baryon = (void*)Psk;
        Psk->Body->Type = TypeSink;
        Psk->FormationTime = Pall.TCurrent;

        Psk->NumberofAbsorbedParticles = 0;
        Psk->PosP[0] = Phydro[i]->PosP[0];
        Psk->PosP[1] = Phydro[i]->PosP[1];
        Psk->PosP[2] = Phydro[i]->PosP[2];

        // Psk->VelP[0] = Phydro[i]->VelP[0];
        // Psk->VelP[1] = Phydro[i]->VelP[1];
        // Psk->VelP[2] = Phydro[i]->VelP[2];
        Psk->Radius = SINK_RADIUS;

        double dt_half_hydro = 0.5*Phydro[i]->dt_hydro;
        Psk->VelP[0] = PhydroBody(i)->Vel[0] = PhydroBody(i)->Velh[0]+dt_half_hydro*Phydro[i]->HydroAcc[0];
        Psk->VelP[1] = PhydroBody(i)->Vel[1] = PhydroBody(i)->Velh[1]+dt_half_hydro*Phydro[i]->HydroAcc[1];
        Psk->VelP[2] = PhydroBody(i)->Vel[2] = PhydroBody(i)->Velh[2]+dt_half_hydro*Phydro[i]->HydroAcc[2];

        HydroRoot.Leaves[Phydro[i]->Leaf] *= -1;
        Phydro[i]->Use = OFF;

        NewSinkParticles ++;
    }
#else
    int NewSinkParticles = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if((PhydroBody(i)->GlobalID %1829 == 0)||((PhydroBody(i)->GlobalID-1) %1829 == 0)){

            // do sink convert!
            StructPsinkptr Psk = ReturnEmptySinkStructurePointer();

            // Setup sink partcile.
            Psk->Use = ON;
            Psk->ParentGlobalID = PhydroBody(i)->GlobalID;
            Psk->Body = PhydroBody(i);
            Psk->Body->Baryon = (void*)Psk;
            Psk->Body->Type = TypeSink;
            Psk->FormationTime = Pall.TCurrent;

            Psk->NumberofAbsorbedParticles = 0;
            Psk->PosP[0] = Phydro[i]->PosP[0];
            Psk->PosP[1] = Phydro[i]->PosP[1];
            Psk->PosP[2] = Phydro[i]->PosP[2];

            Psk->VelP[0] = Phydro[i]->VelP[0];
            Psk->VelP[1] = Phydro[i]->VelP[1];
            Psk->VelP[2] = Phydro[i]->VelP[2];
            Psk->Radius = SINK_RADIUS;

            double dt_half_hydro = 0.5*Phydro[i]->dt_hydro;
            Psk->PosP[0] = PhydroBody(i)->Vel[0] = PhydroBody(i)->Velh[0]+dt_half_hydro*Phydro[i]->HydroAcc[0];
            Psk->PosP[1] = PhydroBody(i)->Vel[1] = PhydroBody(i)->Velh[1]+dt_half_hydro*Phydro[i]->HydroAcc[1];
            Psk->PosP[2] = PhydroBody(i)->Vel[2] = PhydroBody(i)->Velh[2]+dt_half_hydro*Phydro[i]->HydroAcc[2];

            HydroRoot.Leaves[Phydro[i]->Leaf] *= -1;
            Phydro[i]->Use = OFF;

            NewSinkParticles ++;
        }
    }
#endif

    int GlobalNewSinkParticles;
    MPI_Allreduce(&NewSinkParticles,&GlobalNewSinkParticles,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
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
    Pall.NActivesSink = Pall.Nsink; 
    Pall.NActivesSink_t = Pall.Nsink_t; 

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The number of new sink particles = %d\n",GlobalNewSinkParticles);
        fflush(stderr);
    }

    return ;
}

static void OutputParticleDistribution(const int flag){

    FILE *fp;
    char fname[MaxCharactersInLine];
    double Mass = 0.e0;
    for(int k=0;k<MPIGetNumProcs();k++){
        if(MPIGetMyID() == k){
            Snprintf(fname,"./data/H.%02d.%02d.%02d",flag,MPIGetNumProcs(),MPIGetMyID());
            FileOpen(fp,fname,"w");
            for(int i=0;i<Pall.Nhydro;i++){
                fprintf(fp,"%ld %g %g %g %g\n",PhydroBody(i)->GlobalID,
                        Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2],PhydroBody(i)->Mass);
                Mass += PhydroBody(i)->Mass;
            }
            fclose(fp);

            Snprintf(fname,"./data/Hs.%02d.%02d.%02d",flag,MPIGetNumProcs(),MPIGetMyID());
            FileOpen(fp,fname,"w");
            for(int i=0;i<Pall.Nhydro;i++){
                if(5*MeanDistance > fabs(Phydro[i]->PosP[2]))
                    fprintf(fp,"%ld %g %g %g %g\n",PhydroBody(i)->GlobalID,
                            Phydro[i]->PosP[0],Phydro[i]->PosP[1],Phydro[i]->PosP[2],PhydroBody(i)->Mass);
            }
            fclose(fp);


            Snprintf(fname,"./data/S.%02d.%02d.%02d",flag,MPIGetNumProcs(),MPIGetMyID());
            FileOpen(fp,fname,"w");
            for(int i=0;i<Pall.Nsink;i++){
                fprintf(fp,"%ld %g %g %g %g %g %d\n",PsinkBody(i)->GlobalID,
                        Psink[i]->PosP[0],Psink[i]->PosP[1],Psink[i]->PosP[2],PsinkBody(i)->Mass,
                        Psink[i]->Radius,Psink[i]->NumberofAbsorbedParticles);
                Mass += PsinkBody(i)->Mass;
            }
            fclose(fp);

            //fflush(fp);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    double GlobalMass;
    MPI_Allreduce(&Mass,&GlobalMass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"The total mass is %g\n",GlobalMass);
        fflush(stderr);
    }

    // Center of mass, vel.
    Mass = 0.e0;
    double COM[3] = {0.e0,0.e0,0.e0};
    double VCOM[3] = {0.e0,0.e0,0.e0};
    for(int i=0;i<Pall.Ntotal;i++){
        Mass += Pbody[i]->Mass;
        COM[0] += Pbody[i]->Mass*Pbody[i]->Pos[0];
        COM[1] += Pbody[i]->Mass*Pbody[i]->Pos[1];
        COM[2] += Pbody[i]->Mass*Pbody[i]->Pos[2];

        VCOM[0] += Pbody[i]->Mass*Pbody[i]->Vel[0];
        VCOM[1] += Pbody[i]->Mass*Pbody[i]->Vel[1];
        VCOM[2] += Pbody[i]->Mass*Pbody[i]->Vel[2];
    }
    double GMass = 0.e0;
    double GCOM[3] = {0.e0,0.e0,0.e0};
    double GVCOM[3] = {0.e0,0.e0,0.e0};
    MPI_Allreduce(&Mass,&GMass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(COM,GCOM,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(VCOM,GVCOM,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for(int k=0;k<3;k++){
        GCOM[k] /= GMass;
        GVCOM[k] /= GMass;
    }
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Mass = %g\n",GMass);
        fprintf(stderr,"COM = %g %g %g\n",GCOM[0],GCOM[1],GCOM[2]);
        fprintf(stderr,"VCOM = %g %g %g\n",GVCOM[0],GVCOM[1],GVCOM[2]);
    }

    Mass = 0.e0;
    COM[0] = COM[1] = COM[2] = 0.e0;
    VCOM[0] = VCOM[1] = VCOM[2] = 0.e0;
    for(int i=0;i<Pall.Nsink;i++){
        Mass += PsinkBody(i)->Mass;
        COM[0] += PsinkBody(i)->Mass*PsinkBody(i)->Pos[0];
        COM[1] += PsinkBody(i)->Mass*PsinkBody(i)->Pos[1];
        COM[2] += PsinkBody(i)->Mass*PsinkBody(i)->Pos[2];

        VCOM[0] += PsinkBody(i)->Mass*PsinkBody(i)->Vel[0];
        VCOM[1] += PsinkBody(i)->Mass*PsinkBody(i)->Vel[1];
        VCOM[2] += PsinkBody(i)->Mass*PsinkBody(i)->Vel[2];
    }
    for(int k=0;k<3;k++){
        COM[k] /= Mass;
        VCOM[k] /= Mass;
    }
    fprintf(stderr,"Sink Mass = %g\n",Mass);
    fprintf(stderr,"Sink COM = %g %g %g\n",COM[0],COM[1],COM[2]);
    fprintf(stderr,"Sink VCOM = %g %g %g\n",VCOM[0],VCOM[1],VCOM[2]);

    return ;
}

int main_Test_Sink(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv); 

    // make uniform grid.
#define BaseNumberForGrid (20)
    MakeGridData(BaseNumberForGrid);

    BuildPredictors(); // Pos -> PosP/ Vel -> VelP
    InitializeDecomposition();
    DomainDecomposition();

    InitializeRootForHydro();
    ClearHydroData();
    PlantHydroTree();
    CalcKernel();

    ConvertSinks(100,BaseNumberForGrid);

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Nhydro,Nsink = %d,%d\n",Pall.Nhydro_t,Pall.Nsink_t);
    
    for(int i=0;i<Pall.Nsink;i++)
        //Psink[i]->Radius = 0.1*MeanDistance;
        Psink[i]->Radius = 5.1*MeanDistance;
        //Psink[i]->Radius = 17.1*MeanDistance;
        //Psink[i]->Radius = 1.1*MeanDistance;
        //Psink[i]->Radius = 100.1*MeanDistance;

    // Check mass.
    OutputParticleDistribution(0);

    SinkHydroAccretionExported();
    //SinkSinkMergingExported();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Nhydro,Nsink = %d,%d\n",Pall.Nhydro_t,Pall.Nsink_t);

    // Check mass.
    OutputParticleDistribution(1);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    exit(1);
    return ;
}
#endif //TASK_TEST_SINKPARTICLE //}

#ifdef TASK_TEST_FOF //{
int main_FOFTest(const int argc, char *argv[]){

    double LinkingLength;

    InitializeCommunicationTable(); 
    InitializeCommunicationBuffers(); 

    int TotalFileNumber = GetRunStatusForSingleNodeAnalysis(argc,argv);
    TotalFileNumber = 4;
    if(Pall.RunStatus == NewSimulation){
        InitUniformSphereTest(20);
        //LinkingLength = 0.45; 
        LinkingLength = 0.1; 
    } else if (Pall.RunStatus == RestartSimulation){
        //ReadParallelDataOnSingleNode(TotalFileNumber); // interface changed.
        FILE *fp;
        FileOpen(fp,"grid.dx.astart","r");
        fscanf(fp,"%*d %*d %*d %le %*g",&LinkingLength);
        fclose(fp);
        //LinkingLength = 0.45; 
        LinkingLength *= 0.2/(1.0+Pall.Redshift); 
        Pall.RunStatus = RestartSimulation;
    } else {
        exit(RunStatusError);
    }

    BuildPredictors();

    InitializeRootForGravity();
    PlantGravityTree();

    AllocateStructFOF();
    SetFOFLinkingLength(LinkingLength);

    FOFHaloFinder();
    ListUpFOFCatalog();

    WriteFOFCatalog();

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_FOF //}

#ifdef TASK_TEST_STELLARFEEDBACK // TASK_TEST_STELLARFEEDBACK //{
#include <CELib.h>
/*
 * Test routines for stellar feedback with CELib library.
 * Test1: Check delayed time distribution.
 * Test2: Check feedback radius.
 * Test3: Check yield, energy distribution.
 * Test4: Check SNeII Individual mode.
 * Test5: Check SNeII Individual mode with association.
 */

static void StellarFeedbackTest1(void){

    if(MPIGetMyID() != MPI_ROOT_RANK) 
        return ;

    fprintf(stderr,"####################\n");
    fprintf(stderr,"##  Start Test 1  ##\n");
    fprintf(stderr,"####################\n");


    int Number = 10000;
    InitializeRandomGenerator(1977+MPIGetMyID());

    // Test Type II DTD.
    FILE *fp;
    FileOpen(fp,"./DTD2.dat","w");

    for(int i=0;i<Number;i++){
        double R = gsl_rng_uniform(RandomGenerator);
        fprintf(fp,"%g %g\n",R,CELibGetSNIIExplosionTime(R,0.02));
    }
    fclose(fp);


    // Test Type Ia DTD.
    FileOpen(fp,"./DTDIa.dat","w");
    for(int i=0;i<Number;i++){
        double R = gsl_rng_uniform(RandomGenerator);
        fprintf(fp,"%g %g %g\n",R,
            CELibGetSNIaExplosionTime(R,0.02,1000.0,0),
            CELibGetSNIaExplosionTime(R,0.02,1000.0,1));
    }
    fclose(fp);

    return ;
}

void MakeParticleDistributionForStellarFeedbackTest2(const double nH, const double Radius, const int Nhydro, const int Nstars, double StarPos[restrict][3], const double Mstar_in_SolarMass){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    MakeDir("./data");

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();


    /* Number of particles N is modified to N' = 2^p, where p = (int)log2(N) */
    int m = (int)log2(Nhydro);
    int Nhydro_p2 = 1<<m;
    dprintlmpi(Nhydro_p2);

    double ParticleMass = (4*M_PI/3.0)*(nH*PROTON_MASS_CGS/0.76)*CUBE(Radius*PC_CGS)
                            /MSUN_CGS/(double)Nhydro_p2;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"ParticleMass = %g, Total Mass = %g\n",ParticleMass,ParticleMass*Nhydro_p2);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    } else {
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = 0;
        Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    }
    fprintf(stderr,"[%02d] %ld %ld %ld %ld %ld | %ld %ld %ld %ld %ld\n",MPIGetMyID(),
            Pall.Nhydro,Pall.Nstars,Pall.Nsink,Pall.NDM,Pall.Ntotal,
            Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t,Pall.NDM_t,Pall.Ntotal_t);


    int AllocationSize = Pall.Ntotal; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        GenerateStructPstar(Nstars);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Nstars;i++){
            Pbody[Pall.Nhydro+i]->Baryon = (void *)(Pstar[i]);
            Pstar[i]->Body = Pbody[Pall.Nhydro+i];
        }
    }

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    double eps = 0.1*PC_CGS/Pall.UnitLength;
    double Uinit = 100*Pall.ConvertTtoU;
    double cs = sqrt(Pall.GGm1*Uinit);

#if 1
    int count = 0;
    int count_passall = 0;
    while(count_passall < Pall.Nhydro_t){
        double Pos[] ={Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
        double r = NORM(Pos);
        if(r<Radius){
            // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
            if(count_passall%NProcs == MyID){
                Pbody[count]->Active = ON;
                Pbody[count]->Use = ON;
                Pbody[count]->Type = TypeHydro;
                Pbody[count]->GlobalID = count_passall;

                Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
                Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
                Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];

                Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 1.e-30;
                Pbody[count]->Acc[0] = Pbody[count]->Acc[1] = Pbody[count]->Acc[2] = 0.e0;

                Pbody[count]->Eps = eps;
                Pbody[count]->Mass = ParticleMass;
                
                PbodyHydro(count)->Mass = Pbody[count]->Mass = ParticleMass;
                CELibSetPrimordialMetallicity(ParticleMass,PbodyHydro(count)->Elements);

                double MassLightElements = Phydro[count]->Elements[CELibYield_H]+Phydro[count]->Elements[CELibYield_He];
                double Z = (ParticleMass-MassLightElements)/ParticleMass;

                PbodyHydro(count)->Z = PbodyHydro(count)->ZII = PbodyHydro(count)->ZIa = Z;
#if (UseSFModelSpawn) 
                PbodyHydro(count)->SpawnMass = ParticleMass/((double)MaxSpawnTimes);
#endif

#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
                PbodyHydro(count)->G0 = 10;
                PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

                PbodyHydro(count)->DQheat = 0.e0;
                PbodyHydro(count)->Use = ON;
                PbodyHydro(count)->Kernel = eps;
                PbodyHydro(count)->U = Uinit;

                count ++;
            }
            count_passall ++;
        }
    }


#else
    int count = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        double R = Radius*gsl_rng_uniform(RandomGenerator);
        double costh = 2.0*gsl_rng_uniform(RandomGenerator)-1;
        double phi = 360*gsl_rng_uniform(RandomGenerator);
        double phi_rad = (M_PI/180.0)*phi;
        if(i%NProcs == MyID){
            // fprintf(stderr,"%d\n",i);
            double Pos[3] = {R*sqrt(1-SQ(costh))*cos(phi_rad),
                            R*sqrt(1-SQ(costh))*sin(phi_rad),
                            R*costh};
            
            Pbody[count]->Active = ON;
            Pbody[count]->Use = ON;
            Pbody[count]->Type = TypeHydro;
            Pbody[count]->GlobalID = i;

            Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
            Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
            Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];

            Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 1.e-30;
            Pbody[count]->Acc[0] = Pbody[count]->Acc[1] = Pbody[count]->Acc[2] = 0.e0;

            Pbody[count]->Eps = eps;
            Pbody[count]->Mass = ParticleMass;

            PbodyHydro(count)->Mass = Pbody[count]->Mass = ParticleMass;
            CELibSetPrimordialMetallicity(ParticleMass,PbodyHydro(count)->Elements);

            double MassLightElements = Phydro[count]->Elements[CELibYield_H]+Phydro[count]->Elements[CELibYield_He];
            double Z = (ParticleMass-MassLightElements)/ParticleMass;

            PbodyHydro(count)->Z = PbodyHydro(count)->ZII = PbodyHydro(count)->ZIa = Z;
            // PbodyHydro(count)->Z = 0.02;
            // PbodyHydro(count)->ZII = 0.02;
            // PbodyHydro(count)->ZIa = 0.02;
#if (UseSFModelSpawn) 
            PbodyHydro(count)->SpawnMass = ParticleMass/((double)MaxSpawnTimes);
#endif

#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
            PbodyHydro(count)->G0 = 10;
            PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

            PbodyHydro(count)->DQheat = 0.e0;
            PbodyHydro(count)->Use = ON;
            PbodyHydro(count)->Kernel = eps;
            PbodyHydro(count)->U = Uinit;

            count ++;
        }
    }
#endif

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }
    

// star particles.
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Pall.Nstars;i++){
#if 0
            int index = Pall.Nhydro+i;
            Pbody[index]->Active = ON;
            Pbody[index]->Use = ON;
            Pbody[index]->Type = TypeStar;
            Pbody[index]->GlobalID = Pall.Nhydro_t+i;

            Pbody[index]->Pos[0] = StarPos[i][0];
            Pbody[index]->Pos[1] = StarPos[i][1];
            Pbody[index]->Pos[2] = StarPos[i][2];
            Pbody[index]->Vel[0] = Pbody[index]->Vel[1] = Pbody[index]->Vel[2] = 0.e0;
            Pbody[index]->Acc[0] = Pbody[index]->Acc[1] = Pbody[index]->Acc[2] = 0.e0;
            Pbody[index]->Eps = 10*eps;
            gprintlmpi(eps);
            gprintlmpi(Pbody[index]->Eps);
            Pbody[index]->Mass = Mstar_in_SolarMass;

            PbodyStar(index)->Use = ON;
            PbodyStar(index)->ParentGlobalID = Pbody[index]->GlobalID;
            PbodyStar(index)->FormationTime = 0.e0;
            PbodyStar(index)->TypeIIProb = true; 
            PbodyStar(index)->InitialMass = Pbody[index]->Mass; 
            fprintf(stderr,"Stellar Mass is %g [Msun]\n",PbodyStar(index)->InitialMass*Pall.UnitMass/MSUN_CGS);
#else
            PstarBody(i)->Active = ON;
            PstarBody(i)->Use = ON;
            PstarBody(i)->Type = TypeStar;
            dlprintlmpi(PstarBody(i)->GlobalID);
            PstarBody(i)->GlobalID = Pall.Nhydro_t+i;
            dlprintlmpi(PstarBody(i)->GlobalID);

            PstarBody(i)->Pos[0] = StarPos[i][0];
            PstarBody(i)->Pos[1] = StarPos[i][1];
            PstarBody(i)->Pos[2] = StarPos[i][2];
            PstarBody(i)->Vel[0] = PstarBody(i)->Vel[1] = PstarBody(i)->Vel[2] = 0.e0;
            PstarBody(i)->Acc[0] = PstarBody(i)->Acc[1] = PstarBody(i)->Acc[2] = 0.e0;
            PstarBody(i)->Eps = 10*eps;
            gprintlmpi(eps);
            gprintlmpi(PstarBody(i)->Eps);
            PstarBody(i)->Mass = Mstar_in_SolarMass;

            Pstar[i]->Use = ON;
            Pstar[i]->ParentGlobalID = PstarBody(i)->GlobalID;
            Pstar[i]->FormationTime = 0.e0;
            // Pstar[i]->TypeIIProb = true; 
            Pstar[i]->InitialMass = PstarBody(i)->Mass; 
            fprintf(stderr,"Stellar Mass is %g [Msun]\n",Pstar[i]->InitialMass*Pall.UnitMass/MSUN_CGS);
#endif

        }
    }


    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./data/init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);
    fflush(fp_init);

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char ComName[MaxCharactersInLine];
        Snprintf(ComName,"cat ./data/init.%02d.* | sort -n > ./data/init.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        system(ComName);
    }

#if 1
    int Nbin = 30;
    double Rad[Nbin];
    double Mass[Nbin];
    int Pnum[Nbin];
    double dr = Radius/(double)Nbin;
    // gprintlmpi(dr);
    for(int i=0;i<Nbin;i++){
        Rad[i] = (i+0.5)*dr;
        Mass[i] = 0.e0;
        Pnum[i] = 0;
    }

    for(int i=0;i<Pall.Nhydro;i++){
        double r = NORM(PhydroBody(i)->Pos);
        int int_r = r/dr;
        if(int_r < 0) continue;
        if(int_r >= Nbin) continue;

        Mass[int_r] += PhydroBody(i)->Mass;
        Pnum[int_r] ++;
    }
    double GlobalMass[Nbin];
    MPI_Allreduce(Mass,GlobalMass,Nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    int GlobalPnum[Nbin];
    MPI_Allreduce(Pnum,GlobalPnum,Nbin,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,"./data/rad.dist","w");

        double m = 0.e0; 
        // output mean (number) density.
        for(int i=0;i<Nbin;i++){
            m += Mass[i];
            double density = m*Pall.UnitMass/((4*M_PI/3.0)*CUBE(dr*(i+1)*PC_CGS));
            double Ndensity = 0.76*density/PROTON_MASS_CGS;
            // fprintf(stderr,"Mass %g, R %g, density %g, nH %g\n",
                    // m,dr*(i+1),density,Ndensity);
            // fprintf(fp,"%g %g %d %g %g\n",Rad[i],GlobalMass[i],GlobalPnum[i],density,Ndensity);
        }
        fclose(fp);
    }
#endif

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;
   

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 100.0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 100*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
    }

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.0001*Pall.TEnd;

    strcpy(Pall.ASCIIFileName,"./data/StellarFBTest.ASCII");
    strcpy(Pall.BaseFileName,"./data/StellarFBTest");
    strcpy(Pall.RestartFileName,"./data/StellarFBTest.dump");

    return;
}

static void StellarFeedbackTest2(void){

    fprintf(stderr,"####################\n");
    fprintf(stderr,"##  Start Test 2  ##\n");
    fprintf(stderr,"####################\n");


    // make particle distribution.
    int Nstar = 100;
    int IndexList[Nstar];
    int Nnblist[Nstar];
    int CheckSum[Nstar];
    int Type[Nstar];
    double StarPos[Nstar][3];
    double Radius[Nstar];

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Nstar = 2;
        StarPos[0][0] = StarPos[0][1] = StarPos[0][2] = 0.0;
        StarPos[1][0] = StarPos[1][1] = StarPos[1][2] = +5.0;
    } else {
        Nstar = 0;
    }
    assert(Nstar<100);

    MakeParticleDistributionForStellarFeedbackTest2(100,10,110000,Nstar,StarPos,1000.0);

    for(int i=0;i<Nstar;i++){
        dlprintlmpi(PstarBody(i)->GlobalID);
        dlprintlmpi(Pstar[i]->ParentGlobalID);
        gprintlmpi(PstarBody(i)->Eps);
    }

    InitializeDecomposition();
    DomainDecomposition();

    for(int i=0;i<Pall.Nstars;i++){
        dlprintlmpi(PstarBody(i)->GlobalID);
        dlprintlmpi(Pstar[i]->ParentGlobalID);
        gprintlmpi(PstarBody(i)->Eps);
    }

    // Init and plant tree
    InitializeRootForHydro();
    PlantHydroTree();

    // Type II feedback
    // raise feedback flag
    int count_activestars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        Pstar[i]->EventTimeSNII = 0.e0;
        Pstar[i]->EventTimeSNIa = 0.e0;
        gprintlmpi(PstarBody(i)->Eps);
        IndexList[i] = i;
        Type[i] = CELibFeedbackType_SNII;
        count_activestars ++;
    }
    dprintlmpi(count_activestars);

    // search feedback radii
    StellarFeedbackTestCalcFeedbackRadius(count_activestars,IndexList,Type,Radius,Nnblist,CheckSum);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    int Neighbors[MaxNeighborSize];
    // check size and neighbors
    for(int i=0;i<count_activestars;i++){
        int Index = IndexList[i];
        int Nlist = GetNeighborsLimited(PstarBody(i)->Pos,2.0*Radius[i],Neighbors);
        fprintf(stderr,"-II- ID:%d, Nlist_local = %d, Nlist = %d, CheckSum = %d, Radius = %g, Pos = %g %g %g\n",
                i,Nlist,Nnblist[i],CheckSum[i],Radius[i],
                PstarBody(i)->Pos[0],PstarBody(i)->Pos[1],PstarBody(i)->Pos[2]);
    }


    // Type Ia feedback
    // raise feedback flag
    count_activestars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        Pstar[i]->EventTimeSNII = 0.e0;
        Pstar[i]->EventTimeSNIa = 0.e0;
        IndexList[i] = i;
        Type[i] = CELibFeedbackType_SNIa;
        count_activestars ++;
    }

    // search feedback radii
    StellarFeedbackTestCalcFeedbackRadius(count_activestars,IndexList,Type,Radius,Nnblist,CheckSum);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    // check size and neighbors
    for(int i=0;i<count_activestars;i++){
        int Index = IndexList[i];
        int Nlist = GetNeighborsLimited(PstarBody(i)->Pos,2.0*Radius[i],Neighbors);
        fprintf(stderr,"-Ia- ID:%d, Nlist_local = %d, Nlist = %d, CheckSum = %d, Radius = %g, Pos = %g %g %g\n",
                i,Nlist,Nnblist[i],CheckSum[i],Radius[i],
                PstarBody(i)->Pos[0],PstarBody(i)->Pos[1],PstarBody(i)->Pos[2]);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}

static void StellarFeedbackTestWriteData(char *fname){

    // write hydro info: pos and mass
    FILE *fp;
    char FileName[MaxCharactersInLine];
    Snprintf(FileName,"./data/%s.%03d.%03d",fname,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,FileName,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                //      p  p  p  m  m  H  He C  N  O  Ne Mg Si S  Ca Fe Ni Z
                PhydroBody(i)->GlobalID,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                PhydroBody(i)->Mass,Phydro[i]->Mass,
                Phydro[i]->Elements[CELibYield_H],Phydro[i]->Elements[CELibYield_He],
                Phydro[i]->Elements[CELibYield_C],Phydro[i]->Elements[CELibYield_N],
                Phydro[i]->Elements[CELibYield_O],Phydro[i]->Elements[CELibYield_Ne],
                Phydro[i]->Elements[CELibYield_Mg],Phydro[i]->Elements[CELibYield_Si],
                Phydro[i]->Elements[CELibYield_S],Phydro[i]->Elements[CELibYield_Ca],
                Phydro[i]->Elements[CELibYield_Fe],Phydro[i]->Elements[CELibYield_Ni],
                Phydro[i]->Z);

    }
    fflush(fp);
    fclose(fp);

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char ComName[MaxCharactersInLine];
        Snprintf(ComName,"cat ./data/%s.%03d.* | sort -n > ./data/%s.%03d",
                fname,MPIGetNumProcs(),
                fname,MPIGetNumProcs());
        system(ComName);
    }



    return ;
}

/*
 * Show the total system mass.
 */
static void StellarFeedbackTestCheckMassConservation(void){

    double Mass = 0.e0;
    for(int i=0;i<Pall.Ntotal;i++){
        Mass += Pbody[i]->Mass;
    }
    double GlobalMass;

    MPI_Allreduce(&Mass,&GlobalMass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total system mass is %1.15g\n",GlobalMass);

    return ;
}


static void StellarFeedbackTest3(void){

    fprintf(stderr,"####################\n");
    fprintf(stderr,"##  Start Test 3  ##\n");
    fprintf(stderr,"####################\n");

    // make particle distribution.
    int Nstar = 100;
    int IndexList[Nstar];
    int Nnblist[Nstar];
    int CheckSum[Nstar];
    int Type[Nstar];
    double StarPos[Nstar][3];
    double Radius[Nstar];

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Nstar = 2;
        StarPos[0][0] = StarPos[0][1] = StarPos[0][2] = 0.0;
        StarPos[1][0] = StarPos[1][1] = StarPos[1][2] = +5.0;
    } else {
        Nstar = 0;
    }
    assert(Nstar<100);

    MakeParticleDistributionForStellarFeedbackTest2(100,10,11000,Nstar,StarPos,1000.0);

    StellarFeedbackTestCheckMassConservation();

    InitializeDecomposition();
    DomainDecomposition();

    // Init and plant tree
    InitializeRootForGravity();
    PlantGravityTree();
    InitializeRootForHydro();
    PlantHydroTree();

    // Type II feedback
    // raise feedback flag
    int count_activestars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        Pstar[i]->EventTimeSNII = 0.e0;
        IndexList[i] = i;
        Type[i] = CELibFeedbackType_SNII;
        count_activestars ++;
    }
    dprintlmpi(count_activestars);

    // search feedback radii
    StellarFeedbackTestCalcFeedbackRadius(count_activestars,IndexList,Type,Radius,Nnblist,CheckSum);

    // release energy and elements.
    StellarFeedbackTestReleaseEnergyHeavyElements(count_activestars,IndexList,Type);

    // Check mass conservation.
    StellarFeedbackTestCheckMassConservation();

    // write pos, mass, and elements.
    StellarFeedbackTestWriteData("TypeII");

    // Type Ia feedback
    // raise feedback flag
    count_activestars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        Pstar[i]->EventTimeSNII = 0.e0;
        IndexList[i] = i;
        Type[i] = CELibFeedbackType_SNIa;
        count_activestars ++;
    }

    // search feedback radii
    StellarFeedbackTestCalcFeedbackRadius(count_activestars,IndexList,Type,Radius,Nnblist,CheckSum);

    // release energy and elements.
    StellarFeedbackTestReleaseEnergyHeavyElements(count_activestars,IndexList,Type);

    // Check mass conservation.
    StellarFeedbackTestCheckMassConservation();

    // write pos, mass, and elements.
    StellarFeedbackTestWriteData("TypeIa");
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}

static void StellarFeedbackTest4(void){

    if(MPIGetMyID() != MPI_ROOT_RANK) 
        return ;


    fprintf(stderr,"####################\n");
    fprintf(stderr,"##  Start Test 4  ##\n");
    fprintf(stderr,"####################\n");

    double Mass = 4.e5;
    int counter = Mass*0.005*10;
    double Z[] = {0.e0,0.013};
    char fname[][MaxCharactersInLine] = {"./SNIIEventInterval.Z0.data","./SNIIEventInterval.Zsun.data"};

    for(int k=0;k<2;k++){

        FILE *fp;
        //char fname[] = "./SNIIEventInterval.Z0.data";
        FileOpen(fp,fname[k],"w");
        for(int i=0;i<counter;i++){

            fprintf(fp,"%d %g %g %g\n",i,
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),

                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII)-
                CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII)
            );
        }
        fclose(fp);
    }

    // 

    return ;
}

static void StellarFeedbackTest5(void){

    if(MPIGetMyID() != MPI_ROOT_RANK) 
        return ;

    fprintf(stderr,"####################\n");
    fprintf(stderr,"##  Start Test 5  ##\n");
    fprintf(stderr,"####################\n");

    double Mass = 4.e5;
    int counter = Mass*0.005*10;
    double Z[] = {0.e0,0.013};
    char fname[][MaxCharactersInLine] = {"./SNIIEventIntervalNa.Z0.data","./SNIIEventIntervalNa.Zsun.data"};

    int step = 10;

    for(int k=0;k<2;k++){

        FILE *fp;
        //char fname[] = "./SNIIEventInterval.Z0.data";
        FileOpen(fp,fname[k],"w");
        for(int i=0;i<counter;i+=step){

            fprintf(fp,"%d %g %g %g\n",i,
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),

                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII)-
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII)
            );
        }
        fclose(fp);
    }


    {
    double Mass = 100;
    int counter = MAX(1,Mass*0.005*10);
    double Z[] = {0.e0,0.013};
    char fname[][MaxCharactersInLine] = {"./SNIIEventIntervalNa.Z0.Ext.data","./SNIIEventIntervalNa.Zsun.Ext.data"};

    int step = 10;

    for(int k=0;k<2;k++){

        FILE *fp;
        FileOpen(fp,fname[k],"w");
        for(int i=0;i<counter;i+=step){

            fprintf(fp,"%d %g %g %g\n",i,
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII),

                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII)-
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = Mass,
                        .Metallicity = Z[k],
                        .Count = i,
                        .Nassociation = step,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII)
            );
        }
        fclose(fp);
    }
    }

    // 

    return ;
}

int main_StellarFeedbackTest(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    //CELibSetRunParameterTestMode(true);
    CELibSetRunParameterTestMode(false);

    CELibSetRunParameterIMFType(CHEMICALEVOLUTION_IMFTYPE);
    CELibSetRunParameterSNIIRange(CHEMICALEVOLUTION_SNII_UPPERMASS,CHEMICALEVOLUTION_SNII_LOWERMASS);
    CELibSetRunParameterSNIaType(CHEMICALEVOLUTION_SNIa_TYPE);
#if CHEMICALEVOLUTION_SNIa_TYPE == 0
    CELibSetRunParameterSNIaRange(CHEMICALEVOLUTION_SNIa_UPPERMASS,CHEMICALEVOLUTION_SNIa_LOWERMASS);
#endif // CHEMICALEVOLUTION_SNIa_TYPE

    CELibInit();

    if(MPIGetMyID() == MPI_ROOT_RANK){
        CELibShowCurrentRunParameters();
    }

    // Test 1
    StellarFeedbackTest1();

    // Test 2
    //StellarFeedbackTest2();

    // Test 3
    // StellarFeedbackTest3();

    // Test 4
    StellarFeedbackTest4();

    // Test 5
    StellarFeedbackTest5();

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_STELLARFEEDBACK //}


#ifdef TASK_TEST_MOMENTUMFEEDBACK //{
void MakeParticleDistributionForMomentumFeedbackTest(const double nH, const double Radius, const int Nhydro, const int Nstars, double StarPos[restrict][3], const double Mstar_in_SolarMass){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    MakeDir("./data");

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();


    /* Number of particles N is modified to N' = 2^p, where p = (int)log2(N) */
    int m = (int)log2(Nhydro);
    int Nhydro_p2 = 1<<m;
    dprintlmpi(Nhydro_p2);

    double ParticleMass = (4*M_PI/3.0)*(nH*PROTON_MASS_CGS/0.76)*CUBE(Radius*PC_CGS)
                            /MSUN_CGS/(double)Nhydro_p2;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"ParticleMass = %g, Total Mass = %g\n",ParticleMass,ParticleMass*Nhydro_p2);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    } else {
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = 0;
        Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    }
    fprintf(stderr,"[%02d] %ld %ld %ld %ld %ld | %ld %ld %ld %ld %ld\n",MPIGetMyID(),
            Pall.Nhydro,Pall.Nstars,Pall.Nsink,Pall.NDM,Pall.Ntotal,
            Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t,Pall.NDM_t,Pall.Ntotal_t);


    int AllocationSize = Pall.Ntotal; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        GenerateStructPstar(Nstars);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Nstars;i++){
            Pbody[Pall.Nhydro+i]->Baryon = (void *)(Pstar[i]);
            Pstar[i]->Body = Pbody[Pall.Nhydro+i];
        }
    }

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    double eps = 0.1*PC_CGS/Pall.UnitLength;
    double Uinit = 100*Pall.ConvertTtoU;
    double cs = sqrt(Pall.GGm1*Uinit);

    int count = 0;
    int count_passall = 0;
    while(count_passall < Pall.Nhydro_t){
        double Pos[] ={Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
        double r = NORM(Pos);
        if(r<Radius){
            // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
            if(count_passall%NProcs == MyID){
                Pbody[count]->Active = ON;
                Pbody[count]->Use = ON;
                Pbody[count]->Type = TypeHydro;
                Pbody[count]->GlobalID = count_passall;

                Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
                Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
                Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];

                Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 1.e-30;
                Pbody[count]->Acc[0] = Pbody[count]->Acc[1] = Pbody[count]->Acc[2] = 0.e0;

                Pbody[count]->Eps = eps;
                Pbody[count]->Mass = ParticleMass;
                
                PbodyHydro(count)->Mass = Pbody[count]->Mass = ParticleMass;
                CELibSetSolarMetallicity(ParticleMass,PbodyHydro(count)->Elements);
                double MassLightElements = Phydro[count]->Elements[CELibYield_H]+Phydro[count]->Elements[CELibYield_He];
                double Z = (ParticleMass-MassLightElements)/ParticleMass;
                PbodyHydro(count)->Z = PbodyHydro(count)->ZII = PbodyHydro(count)->ZIa = Z;

#if (UseSFModelSpawn) 
                PbodyHydro(count)->SpawnMass = ParticleMass/((double)MaxSpawnTimes);
#endif

#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
                PbodyHydro(count)->G0 = 10;
                PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

                PbodyHydro(count)->DQheat = 0.e0;
                PbodyHydro(count)->Use = ON;
                PbodyHydro(count)->Kernel = eps;
                PbodyHydro(count)->U = Uinit;

                count ++;
            }
            count_passall ++;
        }
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }
    

// star particles.
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Pall.Nstars;i++){
            PstarBody(i)->Active = ON;
            PstarBody(i)->Use = ON;
            PstarBody(i)->Type = TypeStar;
            //dlprintlmpi(PstarBody(i)->GlobalID);
            PstarBody(i)->GlobalID = Pall.Nhydro_t+i;
            dlprintlmpi(PstarBody(i)->GlobalID);

            PstarBody(i)->Pos[0] = StarPos[i][0];
            PstarBody(i)->Pos[1] = StarPos[i][1];
            PstarBody(i)->Pos[2] = StarPos[i][2];
            PstarBody(i)->Vel[0] = PstarBody(i)->Vel[1] = PstarBody(i)->Vel[2] = 0.e0;
            PstarBody(i)->Acc[0] = PstarBody(i)->Acc[1] = PstarBody(i)->Acc[2] = 0.e0;
            PstarBody(i)->Eps = 10*eps;
            PstarBody(i)->Mass = Mstar_in_SolarMass;

            CELibSetSolarMetallicity(PstarBody(i)->Mass,Pstar[i]->Elements);
            Pstar[i]->Z = (Pstar[i]->Mass
                    -Pstar[i]->Elements[CELibYield_H]-Pstar[i]->Elements[CELibYield_He])/PstarBody(i)->Mass;
            Pstar[i]->ZII = Pstar[i]->Z;
            Pstar[i]->ZIa = Pstar[i]->Z;

            Pstar[i]->Use = ON;
            Pstar[i]->ParentGlobalID = PstarBody(i)->GlobalID;
            Pstar[i]->FormationTime = 0.e0;
            Pstar[i]->InitialMass = PstarBody(i)->Mass; 
            fprintf(stderr,"Stellar Mass is %g [Msun]\n",Pstar[i]->InitialMass*Pall.UnitMass/MSUN_CGS);
        }
    }


    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./data/init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);
    fflush(fp_init);

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char ComName[MaxCharactersInLine];
        Snprintf(ComName,"cat ./data/init.%02d.* | sort -n > ./data/init.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        system(ComName);
    }

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    fprintf(stderr,"[%02d] %ld %ld %ld %ld %ld | %ld %ld %ld %ld %ld\n",MPIGetMyID(),
            Pall.Nhydro,Pall.Nstars,Pall.Nsink,Pall.NDM,Pall.Ntotal,
            Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t,Pall.NDM_t,Pall.Ntotal_t);


    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 128;
    Pall.Npm =  8;
    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;
   

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.0;
    Pall.ViscousS = 100.0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    Pall.TEnd = 1*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
    }

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.0001*Pall.TEnd;

    strcpy(Pall.ASCIIFileName,"./data/MomentumFBTest.ASCII");
    strcpy(Pall.BaseFileName,"./data/MomentumFBTest");
    strcpy(Pall.RestartFileName,"./data/MomentumFBTest.dump");

    return;
}

static void PlotData(void){

    for(int i=0;i<Pall.Nstars;i++){
        fprintf(stderr,"S: %g %g %g | %g %g\n",
                PstarBody(i)->Pos[0],PstarBody(i)->Pos[1],PstarBody(i)->Pos[2],
                PstarBody(i)->Mass,
                Pstar[i]->LFUV);
    }

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"tmpfb.%03d.%03d",MPIGetNumProcs(),MPIGetMyID());
    //FileOpen(fp,"fb.dat","w");
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%g %g %g %g %g %g %g %g\n",
                PhydroBody(i)->PosP[0],
                PhydroBody(i)->PosP[1],
                PhydroBody(i)->PosP[2],
                PhydroBody(i)->Vel[0],
                PhydroBody(i)->Vel[1],
                PhydroBody(i)->Vel[2],
                NORM(PhydroBody(i)->PosP),
                DOT_PRODUCT(PhydroBody(i)->PosP,PhydroBody(i)->Vel)/NORM(PhydroBody(i)->PosP)
                );
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"cat tmpfb.%03d.??? | sort -n > fb.%03d.dat",
                MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",command);
        system(command);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"rm tmpfb.%03d.???",MPIGetNumProcs());
        fprintf(stderr,"%s\n",command);
        system(command);
    }

    return ;
}

int main_MomentumFeedbackTest(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    MakeDir("Log");

    // set init
    // make particle distribution.
    int Nstar = 100;
    int IndexList[Nstar];
    int Nnblist[Nstar];
    int CheckSum[Nstar];
    double StarPos[Nstar][3];
    double Radius[Nstar];

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Nstar = 1;
        StarPos[0][0] = StarPos[0][1] = StarPos[0][2] = 0.0;
    } else {
        Nstar = 0;
    }
    assert(Nstar<100);

    MakeParticleDistributionForMomentumFeedbackTest(1,10,10000,1,StarPos,1000);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pstar[0]->EventTimeSNII = 0.e0;
    }

    InitializeDecomposition();
    DomainDecomposition();
    BuildPredictors();

    InitializeRootForHydro();
    PlantHydroTree();

    InitStellarWind();
    InitializeStarFormation();
    InitializeStellarFeedback();
    ClearHydroData();
    PlantHydroTree();

    Pall.TCurrent = 1*YEAR_CGS/Pall.UnitTime;
    CalcSize();

    for(int i=0;i<Pall.Nhydro;i++){
        double R = NORM(PhydroBody(i)->Pos);
        // if(R > 9){
            // Phydro[i]->Kernel = 
            // Phydro[i]->KernelPred *= 10; 
        // }
#if 0
        if((R<3)&&(R>2)){
            Phydro[i]->Kernel = 
            Phydro[i]->KernelPred *= 100; 
        }
#endif
        Phydro[i]->Kernel = 
        Phydro[i]->KernelPred *= 5*gsl_rng_uniform(RandomGenerator)*(Phydro[i]->Kernel); 
    }
    PlantHydroTree();


    CalcDensityDivRot();

    //StellarFeedback();
    StellarFeedbackGatherScatterESA();

    // plot particle data
    PlotData();

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_MOMENTUMFEEDBACK //}

#ifdef TASK_TEST_EQUILIBRIUM_TEMPERATURE //{
#define EQ_INITIAL_TEMPERATURE  (1.e4)

void InitializeParticleDistributionForEquilibriumTemperatureTest(const int Number, 
        const double EndTime /* in cgs */, const int OutputNumber){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MPIGetMyID());


    Pall.UnitLength = KPC_CGS;
    Pall.UnitTime = KPC_CGS/VELOCITY_KMS_CGS;
    Pall.UnitMass = 1.e+10*MSUN_CGS;

    Pall.TCMB = CMB_TEMPERATURE;
    Pall.hubble = 0.7;
    Pall.Hubble = Pall.hubble*100*(Pall.UnitLength/MPC_CGS)*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    //Pall.OmegaM = 1.0;
    //Pall.OmegaL = 0.0; 

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 10.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    int mycount_body = 0;
    for(int k=0;k<Number;k+=NProcs)
        mycount_body ++;

    Pall.Ntotal = Pall.Nhydro = mycount_body;
    Pall.Ntotal_t = Pall.Nhydro_t = Number;

    PbodySize = Pall.Ntotal;
    PhydroSize = Pall.Nhydro;
    PbodyElementsSize = Pall.Ntotal;
    PhydroElementsSize = Pall.Nhydro;

    dlprintlmpi(Pall.NDM);
    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.Nhydro);

    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<PhydroSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double Rhomin = -5;
    double Rhomax = +5;
    double dRho = (Rhomax-Rhomin)/Pall.Ntotal;

    ////////////////////////////////////////////////////////////////
    for(int i=0;i<Pall.Ntotal;i++){

        double r = gsl_rng_uniform(RandomGenerator);
        double theta = 2.0*M_PI*gsl_rng_uniform(RandomGenerator);
        double Pos[3];
        Pos[0] = r*cos(theta);
        Pos[1] = r*sin(theta);
        Pos[2] = gsl_ran_gaussian(RandomGenerator,1.0);

        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = (i*NProcs+MyID);

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = Pos[2];

        Pbody[i]->PosP[0] = Pbody[i]->Pos[0];
        Pbody[i]->PosP[1] = Pbody[i]->Pos[1];
        Pbody[i]->PosP[2] = Pbody[i]->Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = 1.0;
        Pbody[i]->Eps = 1.0;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = 1.0;
        PbodyHydro(i)->UPred = PbodyHydro(i)->U = Pall.ConvertTtoU*EQ_INITIAL_TEMPERATURE;

        PbodyHydro(i)->Z   = 0.02;
        PbodyHydro(i)->ZII = 0.02;
        PbodyHydro(i)->ZIa = 0.00;
        PbodyHydro(i)->RhoPred = PbodyHydro(i)->Rho = pow(10.0,Rhomin + i*dRho)/Pall.ConvertNumberDensityToCGS;

        PbodyHydro(i)->Du = 0.e0;
        Pbody[i]->Acc[0] = PbodyHydro(i)->HydroAcc[0] =
        Pbody[i]->Acc[1] = PbodyHydro(i)->HydroAcc[1] =
        Pbody[i]->Acc[2] = PbodyHydro(i)->HydroAcc[2] = 0.e0;
    }

    ClearGravitationalForce();

    //AcctivateAllparticles();
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
#ifdef USE_VARIABLE_ALPHA //{
        Phydro[i]->Alpha = Pall.ViscousAlphaMax;
#endif // USE_VARIABLE_ALPHA //}
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Active = ON;
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
    }

    // count all active particles.
    int NActives = 0;
    int NActivesDM = 0;
    int NActivesHydro = 0;
    int NActivesStars = 0;
    int NActivesSink = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
            NActivesHydro ++;
        } else if(Pbody[i]->Type == TypeStar){
            NActivesStars ++;
        } else if(Pbody[i]->Type == TypeSink){
            NActivesSink ++;
        } else if(Pbody[i]->Type == TypeDM){
            NActivesDM ++;
        }
        NActives ++;
    }

    Pall.NActivesAll = Pall.NActives = NActives;
    Pall.NActivesDM = NActivesDM;
    Pall.NActivesHydro = NActivesHydro;
    Pall.NActivesStars = NActivesStars;
    Pall.NActivesSink = NActivesSink;

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    // hydro parameters

    Pall.TEnd = EndTime/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = EndTime/Pall.UnitTime/(double)OutputNumber; 

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/EQ.ASCII");
    strcpy(Pall.BaseFileName,"./data/EQ");
    strcpy(Pall.RestartFileName,"./data/EQ.dump");

    return ;
}


static void WritePhaseDiagram(const int Index){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"Phase.%04d.%02d.%03d",Index,NProcs,MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g\n",PhydroBody(i)->GlobalID,
                Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho,Pall.ConvertUtoT*Phydro[i]->U);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"cat ./Phase.%04d.%02d.* | sort -n > ./Phase.%04d.%02d",
                Index,NProcs,Index,NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        sprintf(fname,"rm -rf ./Phase.%04d.%02d.*",Index,NProcs);
        system(fname);
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}





int main_Test_EquilibriumTemperature(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);
    int NParticles = 1000;

    //InitializeParticleDistributionForEquilibriumTemperatureTest(NParticles,10*MEGAYEAR_CGS,100);
    InitializeParticleDistributionForEquilibriumTemperatureTest(NParticles,0.1*MEGAYEAR_CGS,10);
    InitializeCoolingTable();


    double dt = 1000*YEAR_CGS;
    int steps = 10*MEGAYEAR_CGS/dt;

    dt /= Pall.UnitTime;

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->dt_hydro = dt;
    }

    for(int i=0;i<steps;i++){
        dprintlmpi(i);
        CalcCooling();
        FileOutPutConstantInterval();
        if(i%1000==0)
            WritePhaseDiagram(i);
        Pall.TCurrent += dt;
    }


    return EXIT_SUCCESS;
}
#endif // TASK_TEST_EQUILIBRIUM_TEMPERATURE //}


#ifdef TASK_TEST_REINICKE_MEYER_TER_VEHN //{
#error

#define _RMTVSmoothedNumber_  (32)
void AddSmoothedEnergyRMTV(const int PeakIndex, const int number){

    if(MPIGetNumProcs()>1)
        exit(-1);
    
    struct StructDistance *Distances; 
    Distances = malloc(sizeof(struct StructDistance)*Pall.Nhydro);
    for(int i=0;i<Pall.Nhydro;i++){
        Distances[i].dist = DISTANCE(Phydro[i]->PosP,Phydro[PeakIndex]->PosP);
        Distances[i].index = i;
    }

    qsort(Distances,Pall.Nhydro,sizeof(struct StructDistance),(int(*)(const void*, const void*))DistCmp);

    // Peak first Pall.Ns particles.
    int Neighbors[_RMTVSmoothedNumber_];
    for(int i=0;i<_RMTVSmoothedNumber_;i++){
        Neighbors[i] = Distances[i].index;
    }
    double kernel = 0.5*DISTANCE(Phydro[Neighbors[0]]->PosP,Phydro[Neighbors[_RMTVSmoothedNumber_-1]]->PosP);

    double wt = 0;
    for(int i=0;i<_RMTVSmoothedNumber_;i++){
        double r = DISTANCE(Phydro[PeakIndex]->PosP,Phydro[Neighbors[i]]->PosP);
        double w = BlastWaveKernel(r,1.e0/kernel);
        wt += w;
    }

    FILE *fp;
    FileOpen(fp,"heated.dat","w");
    fprintf(fp,"%d\n",_RMTVSmoothedNumber_);

    double iwt = 1.e0/wt;
    eprintlmpi(iwt);
    double Uinit = 1.0/Phydro[PeakIndex]->Mass;
    double totale = 0.0;
    for(int i=0;i<_RMTVSmoothedNumber_;i++){
        double r = DISTANCE(Phydro[PeakIndex]->PosP,Phydro[Neighbors[i]]->PosP);
        //double w = BlastWaveKernel(r,1.e0/Phydro[PeakIndex]->Kernel);
        double w = BlastWaveKernel(r,1.e0/kernel);
        fprintf(stderr,"Add energy[%d] %g -> %g, dU = %g, r %g r/h %g \n",i,
                //Phydro[Neighbors[i]]->U,Phydro[Neighbors[i]]->U+Uinit*w*iwt,Uinit*w*iwt,r,r/Phydro[PeakIndex]->Kernel);
                Phydro[Neighbors[i]]->U,Phydro[Neighbors[i]]->U+Uinit*w*iwt,Uinit*w*iwt,r,r/kernel);
        Phydro[Neighbors[i]]->U += Uinit*w*iwt;
        fprintf(fp,"%ld %g %g %g\n",PhydroBody(Neighbors[i])->GlobalID,
                Phydro[Neighbors[i]]->PosP[0],
                Phydro[Neighbors[i]]->PosP[1],
                Phydro[Neighbors[i]]->PosP[2]);
        totale += Uinit*w*iwt;
    }
    gprintlmpi(totale);
    fclose(fp);
    fflush(fp);

#if 1 //cheker
    FILE *fh;
    char fname[MaxCharactersInLine];
    sprintf(fname,"heated.dat");
    FileOpen(fh,fname,"r");
    int Nsmooth;
    fscanf(fh,"%d",&Nsmooth);
    dprintlmpi(Nsmooth);
    int GID[Nsmooth];
    for(int i=0;i<Nsmooth;i++){
        fscanf(fh,"%d %*g %*g %*g",GID+i);
        fprintf(stderr,"GDI = %d\n",GID[i]);
    }
    fclose(fh);

    sprintf(fname,"heated_check.dat");
    FileOpen(fh,fname,"w");
    for(int k=0;k<Nsmooth;k++){
        for(int i=0;i<Pall.Nhydro;i++){
            if(GID[k] == PhydroBody(i)->GlobalID){
                fprintf(fh,"%ld %g %g %g\n",PhydroBody(i)->GlobalID,
                        Phydro[i]->PosP[0],
                        Phydro[i]->PosP[1],
                        Phydro[i]->PosP[2]);
            }
        }
    }
    fclose(fh);
#endif
    free(Distances);

    return;
}

static void InitReinickemeyerTerVehn(const int Ngrid){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = Ngrid;

    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                if(count%NProcs == MyID){
                    mycount ++;
                }
                count ++;
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));

    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1.0;

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 2.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 1.0;
    Pall.BoxCenter[0] = Pall.BoxCenter[1] = Pall.BoxCenter[2] = 0.e0;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    ///////////////////////////////////////////////////////////////////

	double dx = 1.e0/((double)NX);
	double Lhalf = 0.5;
    double mass = 1.e0/(double)count;
    double Uinit = 1.e-10/mass;
    int PeakIndex = 0;

	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                double Pos[3];

				Pos[0] = i*dx-Lhalf;
				Pos[1] = j*dx-Lhalf;
				Pos[2] = k*dx-Lhalf;

                if(count%NProcs == MyID){
                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;
                    Pbody[mycount]->Type = TypeHydro;
                    Pbody[mycount]->GlobalID = count;

                    Pbody[mycount]->Pos[0] = Pos[0];
                    Pbody[mycount]->Pos[1] = Pos[1];
                    Pbody[mycount]->Pos[2] = Pos[2];

                    Pbody[mycount]->Vel[0] = TINY;
                    Pbody[mycount]->Vel[1] = 0.e0;
                    Pbody[mycount]->Vel[2] = 0.e0;

                    Pbody[mycount]->Mass = mass;
                    Pbody[mycount]->Eps = 1.0;

                    PbodyHydro(mycount)->Use = ON;
                    PbodyHydro(mycount)->Kernel = 2.0*dx;
                    PbodyHydro(mycount)->U = Uinit;
                    if( (i==NX/2)&&(j==NY/2)&&(k==NZ/2) ){
                        PeakIndex = mycount;
                        //PbodyHydro(mycount)->U = 1.e0/mass;
                    }

                    mycount ++;
                }
                count ++;
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }


    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    //Pall.Ns = 64;
    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = 0.041;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(3.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/BlastWave.ASCII");
    strcpy(Pall.BaseFileName,"./data/BlastWave");
    strcpy(Pall.RestartFileName,"./data/BlastWave.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 1.01;

    InitializeRandomGenerator(1977);


#if 0
    Phydro[PeakIndex]->U = 1.0/PhydroMass(PeakIndex);
#else
    //AddSmoothedEnergyBlastWave(PeakIndex);
    AddSmoothedEnergyBlastWaveDirectSearch(PeakIndex, number);
#endif


    return ;
}

int main_Test_ReinickeMeyerTerVehn(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    int Number = 32;
    InitReinickemeyerTerVehn(Ngrid);

    InitLogFiles();
    InitializeRun();

    Run();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_REINICKE_MEYER_TER_VEHN //}


#if (defined(TASK_TEST_1D_THERMAL_CONDUCTIVITY)||defined(TASK_TEST_3D_THERMAL_CONDUCTIVITY)) //{

static void WriteThermalConductivitySolution(const double T_JUMP){

    int counter = 0;
    for(double T=0;T<=Pall.TEnd;T+=Pall.OutPutInterval){
    //for(int i=0;i<(int)(Pall.TEnd/Pall.OutPutInterval+0.5);i++){

        FILE *fp;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"./TCSolution.%03d.dat",counter);
        FileOpen(fp,Fname,"w");

#define Nbin (200)
        double dx = 1.0/Nbin;
        for(int k=0;k<Nbin;k++){
            double pos = -0.5+dx*(k+0.5);
            double Kappa = 1.0/(Phydro[0]->Rho*Pall.ConvertTtoU);
            double X = pos/2.0/sqrt(Kappa*T);
            double U = (50.0/2.0)*erfc(X)+(50.0/T_JUMP/2.0)*erfc(-X);
            fprintf(fp,"%g %g\n",pos+0.5,U);
        }

        fclose(fp);
        counter ++;
    }

    return ;
}
#endif // (defined(TASK_TEST_1D_THERMAL_CONDUCTIVITY)||defined(TASK_TEST_3D_THERMAL_CONDUCTIVITY)) //}


#ifdef TASK_TEST_1D_THERMAL_CONDUCTIVITY //{
static void FlushAccVel(void){

    for(int i=0;i<Pall.Ntotal;i++){
        PhydroBody(i)->Acc[0] = PhydroBody(i)->Acc[1] = PhydroBody(i)->Acc[2] =
        PhydroBody(i)->Vel[0] = PhydroBody(i)->Vel[1] = PhydroBody(i)->Vel[2] =
        PhydroBody(i)->Velh[0] = PhydroBody(i)->Velh[1] = PhydroBody(i)->Velh[2] =
        Phydro[i]->HydroAcc[0] = Phydro[i]->HydroAcc[1] = Phydro[i]->HydroAcc[2] = 0.e0;
    }
    return ;
}

void RunFor1DConductivity(void){

    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        
        if(Pall.TCurrent >= Pall.Era){
            FirstTimeStep();
            BuildHierarchicalTimeStep();
        } else {
            if (10*Pall.NActives_t > Pall.Ntotal_t){
                PreDomainDecomposition();
                DomainDecompositionOnlyDataExchange();
            }
            BuildNewTimeStep();
        }


    ////////////////////////////// check timestep
#if 0
#include "ThermalConductivity.h"
        FILE *fp;
        FileOpen(fp,"./dt.dat","w");
        for(int i=0;i<Pall.Nhydro;i++){
            double dt = Phydro[i]->dt_hydro;
            double Kappai = CalcThermalConductivity(i);
            double dt_diff = SQ(2.0*Phydro[i]->Kernel)/(2*Kappai);
            fprintf(fp,"%d %g %g %g\n",i,dt,dt_diff,Kappai);
        }
        fclose(fp);


        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        exit(-1);

#endif

        RaiseActiveFlags();

        FlushAccVel();
        Kick1Drift(); 
        BuildPredictors();

        if(Pall.NActivesHydro_t>0){ // Hydro
            PlantHydroTree();
            ClearHydroData();
            CalcKernel();
            CalcDensityDivRot(); 
            CalcDuDtAcc();
        }
        // Flush Accerelation

        FlushAccVel();
        Kick2();

        CalcDuDtAccEnergyDensityForCorrection();

        UpdateGravityKickFlag();


        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            OutPutAllParticlesInASCIIFormat();

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            if(Pall.TCurrent >= Pall.OutPutFileNumber*Pall.OutPutInterval){
                FILE *fp;
                char fname[MaxCharactersInLine];
                Snprintf(fname,"./thermal.%03d.dat",Pall.OutPutFileNumber);
                FileOpen(fp,fname,"w");
                for(int i=0;i<Pall.Nhydro;i++){
                    fprintf(fp,"%d %g %g\n",i,PhydroBody(i)->Pos[0],Phydro[i]->U);
                }
                fclose(fp);
                fflush(NULL);
            }

            FileOutPutConstantInterval();
        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }


        DataFullDump();
    }

    return ;
}


int main_Test_1DThermalConductivity(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

#define T_JUMP (100)
    Init1DThermalConductivity(1000,T_JUMP);

    InitLogFiles();
    InitializeRun();

    WriteThermalConductivitySolution(T_JUMP);
    RunFor1DConductivity();
#undef T_JUMP 

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif //TASK_TEST_1D_THERMAL_CONDUCTIVITY //}

#ifdef TASK_TEST_3D_THERMAL_CONDUCTIVITY //{
static void FlushAccVel(void){

    for(int i=0;i<Pall.Ntotal;i++){
        PhydroBody(i)->Acc[0] = PhydroBody(i)->Acc[1] = PhydroBody(i)->Acc[2] =
        PhydroBody(i)->Vel[0] = PhydroBody(i)->Vel[1] = PhydroBody(i)->Vel[2] =
        PhydroBody(i)->Velh[0] = PhydroBody(i)->Velh[1] = PhydroBody(i)->Velh[2] =
        Phydro[i]->HydroAcc[0] = Phydro[i]->HydroAcc[1] = Phydro[i]->HydroAcc[2] = 0.e0;
    }
    return ;
}

void RunFor3DConductivity(void){

    while(Pall.TCurrent < Pall.TEnd){

        ClearTimingLogsThisStep();
        
        if(Pall.TCurrent >= Pall.Era){
            FirstTimeStep();
            BuildHierarchicalTimeStep();
        } else {
            if (10*Pall.NActives_t > Pall.Ntotal_t){
                PreDomainDecomposition();
                DomainDecompositionOnlyDataExchange();
            }
            BuildNewTimeStep();
        }


    ////////////////////////////// check timestep
#if 0
#include "ThermalConductivity.h"
        FILE *fp;
        FileOpen(fp,"./dt.dat","w");
        for(int i=0;i<Pall.Nhydro;i++){
            double dt = Phydro[i]->dt_hydro;
            double Kappai = CalcThermalConductivity(i);
            double dt_diff = SQ(2.0*Phydro[i]->Kernel)/(2*Kappai);
            fprintf(fp,"%d %g %g %g\n",i,dt,dt_diff,Kappai);
        }
        fclose(fp);


        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        exit(-1);

#endif


#if 0
        int MyID = MPIGetMyID();
    FILE *fp;
    char fcloud[MaxCharactersInLine];
    Snprintf(fcloud,"_tc.%03d",MyID);
    FileOpen(fp,fcloud,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%d %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->U);
    }
    fclose(fp);
    fflush(NULL);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat _tc.* > tc.dat");
        fflush(NULL);
        system("rm -rf ./_tc.*");
        fflush(NULL);
    }
    fflush(NULL);

        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        exit(-1);
#endif

        RaiseActiveFlags();

        FlushAccVel();
        Kick1Drift(); 
        BuildPredictors();

        if(Pall.NActivesHydro_t>0){ // Hydro
            PlantHydroTree();
            ClearHydroData();
            CalcKernel();
            CalcDensityDivRot(); 
            CalcDuDtAcc();
        }
        // Flush Accerelation

        FlushAccVel();
        Kick2();

        CalcDuDtAccEnergyDensityForCorrection();

        UpdateGravityKickFlag();


        if(Pall.EraStart + Pall.EraLocal >= Pall.Era){
            if(Pall.Nhydro_t>0)
                LogStarFormationRate();
            OutPutAllParticlesInASCIIFormat();

            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;

            if(Pall.TCurrent >= Pall.OutPutFileNumber*Pall.OutPutInterval){
                FILE *fp;
                char fname[MaxCharactersInLine];
                Snprintf(fname,"./thermal3d.%03d.dat",Pall.OutPutFileNumber);
                FileOpen(fp,fname,"w");
                for(int i=0;i<Pall.Nhydro;i++){
                    fprintf(fp,"%d %g %g %g %g\n",i,
                            PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                            Phydro[i]->U);
                            //Pall.ConvertUtoT*Phydro[i]->U);
                }
                fclose(fp);
                fflush(NULL);
            }

            FileOutPutConstantInterval();
            LogOutPutEnergyMomentumAngularMomentum();
        }else{ 
            Pall.TCurrent = Pall.EraStart + Pall.EraLocal;
            Pall.TStepTotal ++;
        }


        DataFullDump();
    }

    return ;
}

int main_Test_3DThermalConductivity(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

#define T_JUMP (100)
    Init3DThermalConductivity(20,0,T_JUMP);

    InitLogFiles();
    InitializeRun();


    WriteThermalConductivitySolution(T_JUMP);
    RunFor3DConductivity();

#undef  T_JUMP 

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif //TASK_TEST_3D_THERMAL_CONDUCTIVITY //}


#ifdef TASK_TEST_FUVFEEDBACK //{
static int Ngas;
static int Nstar;
static int ActiveStep;
static double Radius;
static double Mass;
static double Z;
static double AgeMax;
static double PowerLawIndex;
static double dt;
static void ReadTestFUVFeedbackParams(void){

    if(CheckFile("./fuvparam.txt")){
        FILE *fp;
        FileOpen(fp,"./fuvparam.txt","r");
        fscanf(fp,"%d",&Ngas);
        fscanf(fp,"%d",&Nstar);
        fscanf(fp,"%d",&ActiveStep);
        fscanf(fp,"%le",&Radius);
        fscanf(fp,"%le",&PowerLawIndex);
        fscanf(fp,"%le",&Mass);
        fscanf(fp,"%le",&Z);
        fscanf(fp,"%le",&AgeMax);
        fscanf(fp,"%le",&dt);
        fclose(fp);
    } else {
        Ngas = 10000;
        Nstar = 10000;
        ActiveStep = 1;
        Radius = 10;   // kpc
        PowerLawIndex = 0.0;
        Mass = 1.e11;  // Msun
        Z = 0.013;
        AgeMax = 10*GIGAYEAR_CGS/YEAR_CGS;
        dt = 0.1*MEGAYEAR_CGS/YEAR_CGS;
    }

    return ;
}

static void ShowTestFUVFeedbackParams(void){

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Ngas = %d\n",Ngas);
        fprintf(stderr,"Nstar = %d\n",Nstar);
        fprintf(stderr,"Radius = %g [kpc]\n",Radius);
        fprintf(stderr,"PowerLawIndex = %g\n",PowerLawIndex);
        fprintf(stderr,"Mass = %g [Msun]\n",Mass);
        fprintf(stderr,"Z = %g\n",Z);
        fprintf(stderr,"AgeMax = %g [yr]\n",AgeMax);
    }

    return ;
}

static bool TestFUVFeedbackFirstSetViscosityParameters = true;

static void TestFUVFeedbackSetConstantViscosityParameters(const double Alpha){

    Pall.HydroAlpha = Alpha;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;

    TestFUVFeedbackFirstSetViscosityParameters = false;

    return;
}

static void TestFUVFeedbackSetVariableViscosityParameters(const double AlphaMin, const double AlphaMax, const double SourceFactor, const double InvDumpingScale){

    Pall.ViscousAlphaMin = AlphaMin; // 0.1
    Pall.ViscousAlphaMax = AlphaMax; // 
    Pall.ViscousS = SourceFactor;    // 1.0
    Pall.ViscousL = InvDumpingScale; // 0.1

    TestFUVFeedbackFirstSetViscosityParameters = false;

    return;
}

static void TestFUVFeedbackInsertData(const int TargetIndex, float Pos[], float Vel[], float Mass, float Rho, float U, const double Age, const double Metal, const int Type, const int GlobalID, const double Eps, const double R){

    Pbody[TargetIndex]->GlobalID = GlobalID;
    Pbody[TargetIndex]->Pos[0] = Pos[0];
    Pbody[TargetIndex]->Pos[1] = Pos[1];
    Pbody[TargetIndex]->Pos[2] = Pos[2];

    Pbody[TargetIndex]->Vel[0] = Vel[0];
    Pbody[TargetIndex]->Vel[1] = Vel[1];
    Pbody[TargetIndex]->Vel[2] = Vel[2];
    
    Pbody[TargetIndex]->Mass = Mass;

    Pbody[TargetIndex]->Active = ON;
    Pbody[TargetIndex]->Use = ON;
    Pbody[TargetIndex]->Type = Type;
    Pbody[TargetIndex]->Eps = Eps;

    if(Type == TypeHydro){
        PbodyHydro(TargetIndex)->Use = ON;
        PbodyHydro(TargetIndex)->Rho = Rho;
        PbodyHydro(TargetIndex)->U = U;
        PbodyHydro(TargetIndex)->Mass = Mass;
        PbodyHydro(TargetIndex)->Kernel = 2.0*Eps;
#if (UseSFModelSpawn) 
        PbodyHydro(TargetIndex)->SpawnMass = Pbody[TargetIndex]->Mass/(double)MaxSpawnTimes;
#endif
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
        PbodyHydro(TargetIndex)->G0 = 1;
        PbodyHydro(TargetIndex)->fH2 = 0.1;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

#ifdef USE_CELIB //{
        CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements,Metal);
        double MassLightElements = PbodyHydro(TargetIndex)->Elements[CELibYield_H]
                                  +PbodyHydro(TargetIndex)->Elements[CELibYield_He];
        double Z = (Pbody[TargetIndex]->Mass-MassLightElements)/Pbody[TargetIndex]->Mass;
        PbodyHydro(TargetIndex)->Z = PbodyHydro(TargetIndex)->ZII = PbodyHydro(TargetIndex)->ZIa = Z;
#else 
        PbodyHydro(TargetIndex)->Z   = 0.02;
        PbodyHydro(TargetIndex)->ZII = 0.02;
        PbodyHydro(TargetIndex)->ZIa = 0.00;
#endif // USE_CELIB //}
    } else if(Type == TypeStar){
#include "StellarFeedback.h"
        PbodyStar(TargetIndex)->Use = ON;
        PbodyStar(TargetIndex)->IMFTYPE = CHEMICALEVOLUTION_IMFTYPE;
        PbodyStar(TargetIndex)->Mass = Pbody[TargetIndex]->Mass;
        PbodyStar(TargetIndex)->InitialMass = Pbody[TargetIndex]->Mass;
        PbodyStar(TargetIndex)->FormationTime = Age;
        CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyStar(TargetIndex)->Elements,Metal);

        PbodyStar(TargetIndex)->Z = (PbodyStar(TargetIndex)->Mass-PbodyStar(TargetIndex)->Elements[CELibYield_H]-PbodyStar(TargetIndex)->Elements[CELibYield_He])/PbodyStar(TargetIndex)->Mass;
        PbodyStar(TargetIndex)->ZII = PbodyStar(TargetIndex)->Z;
        PbodyStar(TargetIndex)->ZIa = 0;

        if(fabs(Age)*Pall.UnitTime < 40.0*MEGAYEAR_CGS){
            PbodyStar(TargetIndex)->TypeII = false;
        } else {
            PbodyStar(TargetIndex)->TypeII = true;
        }
        PbodyStar(TargetIndex)->TypeIa = false;
        PbodyStar(TargetIndex)->NthChildren = 0;
        PbodyStar(TargetIndex)->ParentGlobalID = GlobalID;

#ifdef USE_CELIB
        PbodyStar(TargetIndex)->IMFTYPE = StellarFeedbackGetIMFType();
        PbodyStar(TargetIndex)->SNIaCount = -1;
        if(PbodyStar(TargetIndex)->TypeII == true){
            PbodyStar(TargetIndex)->EventTimeSNII = 100*GIGAYEAR_CGS/Pall.UnitTime;
            PbodyStar(TargetIndex)->SNIaCount = 100000;
        } else {
            PbodyStar(TargetIndex)->EventTimeSNII = PbodyStar(TargetIndex)->FormationTime 
                            +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                                .R = R,
                                .InitialMass_in_Msun = PbodyStar(TargetIndex)->InitialMass*Pall.UnitMass/MSUN_CGS,
                                .Metallicity = PbodyStar(TargetIndex)->Z,
                                .Count = 0,
                                },CELibFeedbackType_SNII)
                                *YEAR_CGS/Pall.UnitTime;
        }
#ifdef USE_CELIB_AGB //{
        PbodyStar(TargetIndex)->EventTimeAGB = 100*GIGAYEAR_CGS/Pall.UnitTime;
        PbodyStar(TargetIndex)->AGBCount = 100000000;
#endif //USE_CELIB_AGB //}
#endif //USE_CELIB //}

    }

    return ;
}

void TestFUVFeedbackActivateAllparticles(void){

    if(TestFUVFeedbackFirstSetViscosityParameters == true){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"Need to call SetVariableViscosityParameters before calling %s.\n",__FUNCTION__);
            fflush(NULL);
        }
        assert(TestFUVFeedbackFirstSetViscosityParameters == false);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
#ifdef USE_VARIABLE_ALPHA //{
        Phydro[i]->Alpha = Pall.ViscousAlphaMax;
#endif // USE_VARIABLE_ALPHA //}
        Phydro[i]->Active = PhydroBody(i)->Active;
#ifdef USE_SPSPH //{
        // m/rho = Zw/y
        // y = Zw*rho/m
        // Z = m*y/rho
        // if(Phydro[i]->Rho > 0.e0){
            // Phydro[i]->PseudoDensity = Phydro[i]->Rho;
            // Phydro[i]->Zw = Phydro[i]->Mass;
        // } else {
            //Phydro[i]->Zw = 1.0;
            //Phydro[i]->PseudoDensity = Phydro[i]->Zw*Phydro[i]->Rho/Phydro[i]->Mass;
        // }
        if(Phydro[i]->Rho > 0.e0){
            Phydro[i]->PseudoDensity = 1.e0;
            Phydro[i]->Zw = Phydro[i]->PseudoDensity*Phydro[i]->Mass/Phydro[i]->Rho;
        } else {
            Phydro[i]->Zw = 1.0;
            Phydro[i]->PseudoDensity = Phydro[i]->Zw*Phydro[i]->Rho/Phydro[i]->Mass;
        }
#ifdef TASK_KELVINHELMHOLTZ_INSTABILITY //{
        Phydro[i]->PseudoDensity = 1.e0;
        Phydro[i]->Zw = Phydro[i]->PseudoDensity*Phydro[i]->Mass/Phydro[i]->Rho;
#endif // TASK_HYDROSTATIC //}

#endif //USE_SPSPH //}
    }

    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Active = ON;
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
    }

    // count all active particles.
    int NActives = 0;
    int NActivesDM = 0;
    int NActivesHydro = 0;
    int NActivesStars = 0;
    int NActivesSink = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
            NActivesHydro ++;
        } else if(Pbody[i]->Type == TypeStar){
            NActivesStars ++;
        } else if(Pbody[i]->Type == TypeSink){
            NActivesSink ++;
        } else if(Pbody[i]->Type == TypeDM){
            NActivesDM ++;
        }
        NActives ++;
    }

    Pall.NActivesAll = Pall.NActives = NActives;
    Pall.NActivesDM = NActivesDM;
    Pall.NActivesHydro = NActivesHydro;
    Pall.NActivesStars = NActivesStars;
    Pall.NActivesSink = NActivesSink;

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    return ;
}

static void InitTestFUVFeedback(const int Ngas, const int Nstar, const double Radius_kpc, const double Mass_Msun, const double Z, const double AgeMax_year, const double dt){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977);

    Pall.UnitMass = MSUN_CGS;
    Pall.UnitLength = KPC_CGS;
    Pall.UnitTime = YEAR_CGS;
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    TestFUVFeedbackSetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    TestFUVFeedbackSetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    int mycount_gas = 0;
    for(int i=0;i<Ngas;i++){
        if(i%NProcs == MyID){
            mycount_gas ++;
        }
    }
    int mycount_star = 0;
    for(int i=0;i<Nstar;i++){
        if(i%NProcs == MyID){
            mycount_star ++;
        }
    }
    dprintlmpi(mycount_gas);
    dprintlmpi(mycount_star);

    Pall.Ntotal = mycount_gas+mycount_star;
    Pall.Nhydro = mycount_gas;
    Pall.NDM = 0;
    Pall.Nstars = mycount_star;

    Pall.Ntotal_t = Ngas+Nstar;
    Pall.Nhydro_t = Ngas;
    Pall.NDM_t = 0;
    Pall.Nstars_t = Nstar;

    int AllocationSize = mycount_gas+mycount_star;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    GenerateStructPstar(Pall.Nstars);

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
        Pbody[i]->Type = TypeHydro;
    }
    for(int i=0;i<Pall.Nstars;i++){
        Pbody[i+Pall.Nhydro]->Baryon = (void *)(Pstar[i]);
        Pstar[i]->Body = Pbody[i+Pall.Nhydro];
        Pbody[i+Pall.Nhydro]->Type = TypeStar;
    }

    double m = Mass_Msun/(Ngas+Nstar);
    double V = 1.0/CUBE(2.0*Radius_kpc);
    double rho = Mass_Msun/V;
    double u = 1.e+3*Pall.ConvertTtoU;
    double eps = 1.0*PC_CGS/Pall.UnitLength;

    fprintf(stderr,"Note: Initialization part has some strange behavior...\n");fflush(NULL);

    mycount_gas = 0;
    for(int i=0;i<Ngas;i++){
#if 0
        float Pos[] = {Radius_kpc*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius_kpc*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius_kpc*(2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        float Vel[] = {(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        double r = NORM(Pos);
        if(r > Radius_kpc){
            i--;
            continue;
        }

        if(PowerLawIndex > 0.0){
            Pos[0] /= r*Radius_kpc;
            Pos[1] /= r*Radius_kpc;
            Pos[2] /= r*Radius_kpc;
            r /= Radius_kpc;
            Pos[0] *= (r)*sqrt(r);
            Pos[1] *= (r)*sqrt(r);
            Pos[2] *= (r)*sqrt(r);
            Pos[0] *= Radius_kpc;
            Pos[1] *= Radius_kpc;
            Pos[2] *= Radius_kpc;
        }
#else
        float Pos[] = {(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        float Vel[] = {(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        double r = NORM(Pos);
        if(r > 1.0){
            i--;
            continue;
        }

        if(PowerLawIndex > 0.0){
            Pos[0] /= r;
            Pos[1] /= r;
            Pos[2] /= r;
            Pos[0] *= (r)*sqrt(r);
            Pos[1] *= (r)*sqrt(r);
            Pos[2] *= (r)*sqrt(r);
        }
        Pos[0] *= Radius_kpc;
        Pos[1] *= Radius_kpc;
        Pos[2] *= Radius_kpc;
#endif


        if(i%NProcs == MyID){
            TestFUVFeedbackInsertData(mycount_gas,Pos,Vel,m,rho,u,0,Z,TypeHydro,i,eps,0.0);
            PbodyHydro(mycount_gas)->Kernel = 0.01;
            //PbodyHydro(mycount_gas)->dt_hydro = 0.1*MEGAYEAR_CGS/Pall.UnitTime;
            PbodyHydro(mycount_gas)->dt_hydro = dt*YEAR_CGS/Pall.UnitTime;
            mycount_gas ++;
        }
    }
    gprintlmpi(Phydro[0]->Kernel);
    mycount_star = 0;
    for(int i=0;i<Nstar;i++){
#if 0 
        float Pos[] = {Radius_kpc*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius_kpc*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius_kpc*(2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        float Vel[] = {(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        double r = NORM(Pos);
        if(r > Radius_kpc){
            i--;
            continue;
        }

        if(PowerLawIndex > 0.0){
            Pos[0] /= r*Radius_kpc;
            Pos[1] /= r*Radius_kpc;
            Pos[2] /= r*Radius_kpc;
            r /= Radius_kpc;
            Pos[0] *= (r)*sqrt(r);
            Pos[1] *= (r)*sqrt(r);
            Pos[2] *= (r)*sqrt(r);
            Pos[0] *= Radius_kpc;
            Pos[1] *= Radius_kpc;
            Pos[2] *= Radius_kpc;
        }
#else
        float Pos[] = {(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        float Vel[] = {(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        double r = NORM(Pos);
        if(r > 1.0){
            i--;
            continue;
        }

        if(PowerLawIndex > 0.0){
            Pos[0] /= r;
            Pos[1] /= r;
            Pos[2] /= r;
            Pos[0] *= (r)*sqrt(r);
            Pos[1] *= (r)*sqrt(r);
            Pos[2] *= (r)*sqrt(r);
        }
        Pos[0] *= Radius_kpc;
        Pos[1] *= Radius_kpc;
        Pos[2] *= Radius_kpc;
#endif

        double Age = -AgeMax_year*gsl_rng_uniform(RandomGenerator)*Pall.UnitTime/YEAR_CGS;
        double R = gsl_rng_uniform(RandomGenerator);
        if(i%NProcs == MyID){
            TestFUVFeedbackInsertData(mycount_star+mycount_gas,Pos,Vel,m,rho,u,Age,Z,TypeStar,i+Ngas,eps,R);

            Pstar[mycount_star]->LFUV = (Pall.UnitMass/MSUN_CGS)*PstarBody(mycount_star)->Mass*
                ASRFLXGetFUV(-Pstar[mycount_star]->FormationTime*Pall.UnitTime/YEAR_CGS,Z); 
            mycount_star ++;
        }
    }

    TestFUVFeedbackActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.hubble = 0.01*70;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

    Pall.TEnd = 5*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd; 

    if(MPIGetMyID() == MPI_ROOT_RANK)
        MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/FUVFB.ASCII");
    strcpy(Pall.BaseFileName,"./data/FUVFB");
    strcpy(Pall.RestartFileName,"./data/FUVFB.dump");


    return ;
}

static void WriteDistributionFUVFeedback(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"FUVFeedback.Init.G.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho,
                Phydro[i]->GradRho,Phydro[i]->G0thin,Pall.ConvertUtoT*Phydro[i]->UPred);
    }
    fclose(fp);
    fflush(NULL);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        sprintf(command,"cat FUVFeedback.Init.G.* | sort -n > FUVFeedback.G.%03d.dat",NProcs);
        //system("cat FUVFeedback.Init.G.* | sort -n > FUVFeedback.G.dat");
        system(command);
        fflush(NULL);
        system("rm -rf ./FUVFeedback.Init.G.*");
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    sprintf(fname,"FUVFeedback.Init.S.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nstars;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PstarBody(i)->GlobalID,
                PstarPos(i)[0],PstarPos(i)[1],PstarPos(i)[2],
                PstarVel(i)[0],PstarVel(i)[1],PstarVel(i)[2],
                Pstar[i]->FormationTime,Pstar[i]->LFUV,Pstar[i]->InitialMass);
    }
    fclose(fp);
    fflush(NULL);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        sprintf(command,"cat FUVFeedback.Init.S.* | sort -n > FUVFeedback.S.%03d.dat",NProcs);
        system(command);
        //system("cat FUVFeedback.Init.S.* | sort -n > FUVFeedback.S.dat");
        fflush(NULL);
        system("rm -rf ./FUVFeedback.Init.S.*");
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}

static void WriteDirectFUVFeedback(void){

    const double factor = 1.0/(1.6e-3*4*M_PI*SQ(Pall.UnitLength));
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->G0thin = 0.e0;
        for(int k=0;k<Pall.Nstars;k++){
            if((Pall.TCurrent-Pstar[k]->FormationTime)*Pall.UnitTime < YOUNGSTARTREE_AGELIMIT){
                double Dist2 = DISTANCE2(PhydroPosP(i),PstarPosP(k));
                //Phydro[i]->LFUV += Dist2;
                Phydro[i]->G0thin += factor*Pstar[k]->LFUV/Dist2;
            }
        }
    }


    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"FUVFeedback.Init.D.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        // fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                // PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                // PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                // Phydro[i]->Kernel,Phydro[i]->GradRho,Phydro[i]->G0thin);
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho,
                Phydro[i]->GradRho,Phydro[i]->G0thin,Pall.ConvertUtoT*Phydro[i]->UPred);
    }
    fclose(fp);
    fflush(NULL);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        sprintf(command,"cat FUVFeedback.Init.D.* | sort -n > FUVFeedback.D.%03d.dat",NProcs);
        system(command);
        fflush(NULL);
        system("rm -rf ./FUVFeedback.Init.D.*");
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}

static void CountLFUV(void){

    double LFUV_Gas = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        LFUV_Gas += Phydro[i]->G0thin;
    }

    double LFUV = 0;
    for(int i=0;i<Pall.Nstars;i++){
        LFUV += Pstar[i]->LFUV;
    }
    fprintf(stderr,"Total LFUV_Gas = %g [erg/s/cm^2], LFUV_Star %g [erg/s]\n",LFUV_Gas,LFUV);

    return ;
}

static void ClearLFUV(void){

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->G0thin = 0.e0;
    }

    return ;
}

static void SelectActives(void){

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroBody(i)->GlobalID%ActiveStep == 0){
            Phydro[i]->Active = true;
            PhydroBody(i)->Active = true;
        } else {
            Phydro[i]->Active = false;
            PhydroBody(i)->Active = false;
        }
    }

    return ;
}

int main_TestFUVFeedback(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    CELibSetRunParameterIntegrationSteps(100);
    InitializeStellarFeedback();
    //ASRFLXRunParameters.IMFWeighted = true;
    ASRFLXRunParameters.IMFType = CHEMICALEVOLUTION_IMFTYPE;
    ASRFLXInitFUV();

    ReadTestFUVFeedbackParams();

    InitTestFUVFeedback(Ngas,Nstar,Radius,Mass,Z,AgeMax,dt);


    BuildPredictors(); // Pos -> PosP/ Vel -> VelP

    SelectActives();

    InitializeDecomposition();
    DomainDecomposition();

    InitializeRootForHydro();
    PlantHydroTree();
    ClearHydroData();
    //CalcKernel();
    CalcSize();
    CalcDensityDivRot();


    WriteDirectFUVFeedback();
    CountLFUV();

    ClearLFUV();
    MPI_Barrier(MPI_COMM_WORLD);
    double t = GetElapsedTime();
    CalcFUV();
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK) gprintlmpi(GetElapsedTime()-t);
    

    InitializeCoolingTable();
    InitializeFarUltraVioletField();
    CalcCooling();

    WriteDistributionFUVFeedback();
    CountLFUV();

    ShowTestFUVFeedbackParams();

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_FUVFEEDBACK //}


#ifdef TASK_TEST_STROMGRENSPHERE
double nH; // in number/cc.
double Radius; // in pc.
int Nhydro; 
int Nstars;
double (*Pos)[3];
double *LyAlpha;
bool *HIIFlag; 

void ReadTestStromgrenSphereParams(char *fname){

    if(CheckFile(fname)){
        FILE *fp;
        FileOpen(fp,fname,"r");
        fscanf(fp,"%le",&nH);
        fscanf(fp,"%le",&Radius);
        fscanf(fp,"%d",&Nhydro);
        fscanf(fp,"%d",&Nstars);
        Pos = realloc(Pos,sizeof(double)*3*(Nstars));
        LyAlpha = realloc(LyAlpha,sizeof(double)*Nstars);
        for(int i=0;i<Nstars;i++){
            fscanf(fp,"%le %le %le %le",Pos[i],Pos[i]+1,Pos[i]+2,LyAlpha+i);
        }
        fclose(fp);
    }

    return ;
}

static double ReturnStromgrenRadius(const double nH, const double LyAlpha){
    const double Alpha = 2.6e-13; // cm^3 s^-1.
    double r3 = (3.0/(4.0*M_PI))*LyAlpha/(nH*nH*Alpha);
    return cbrt(r3);
}

static void WriteStromgrenSphere(const double Z_in_PC){

    FILE *fp;
    char fname[MaxCharactersInLine];

    // write on
    Snprintf(fname,"./HII.all.on.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(HIIFlag[i])
            fprintf(fp,"%ld %g %g %g\n",PhydroBody(i)->GlobalID,PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./HII.all.on.%02d.* | sort -n > HII.all.on.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }

    // write off
    Snprintf(fname,"./HII.all.off.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(!HIIFlag[i])
            fprintf(fp,"%ld %g %g %g\n",PhydroBody(i)->GlobalID,PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./HII.all.off.%02d.* | sort -n > HII.all.off.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }

    
    // slice on
    Snprintf(fname,"./HII.slice.on.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(fabs(PhydroBody(i)->Pos[2])*Pall.UnitLength/PC_CGS<Z_in_PC)
        if(HIIFlag[i])
            fprintf(fp,"%ld %g %g %g\n",PhydroBody(i)->GlobalID,PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./HII.slice.on.%02d.* | sort -n > HII.slice.on.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // slice off
    Snprintf(fname,"./HII.slice.off.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(fabs(PhydroBody(i)->Pos[2])*Pall.UnitLength/PC_CGS<Z_in_PC)
        if(!HIIFlag[i])
            fprintf(fp,"%ld %g %g %g\n",PhydroBody(i)->GlobalID,PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./HII.slice.off.%02d.* | sort -n > HII.slice.off.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }

    return ;
}

int main_TestStromgrenSphere(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    nH = 1.0; // in number/cc.
    Radius = 10; // in pc.
    Nhydro = 100000; 
    Nstars = 1;
    Pos = malloc(sizeof(double)*3*Nstars);
    LyAlpha = malloc(sizeof(double)*Nstars);
    Pos[0][0] = Pos[0][1] = Pos[0][2] = LyAlpha[0] = 0.e0;

    ReadTestStromgrenSphereParams("StromgrenSphereParams.txt");

    fprintf(stderr,"%g %g %d %d\n",nH,Radius,Nhydro,Nstars);
    for(int i=0;i<Nstars;i++)
        fprintf(stderr,"%g %g %g %g\n",Pos[i][0],Pos[i][1],Pos[i][2],LyAlpha[i]);
    fflush(NULL);

    fprintf(stderr,"Expected Stromgren radii are...\n");
    for(int i=0;i<Nstars;i++){
        fprintf(stderr," %g [pc]\n",ReturnStromgrenRadius(nH,LyAlpha[i])/PC_CGS);
    }

    InitTestStromgrenSphere(nH,Radius,Nhydro,Nstars,Pos,LyAlpha);

    BuildPredictors();
    InitializeDecomposition();
    DomainDecomposition();
    InitializeRootForHydro();
    InitializeCoolingTable();
    InitializeStarFormation();
    InitializeDelayedSNII();
    InitializeHIIregions();

    ClearHydroData();
    BuildPredictors();
    PlantHydroTree();

    // Make pho predictor
#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
    CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
    CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
    CalcDensityDivRot();
    BuildPredictors();
    // Make pho predictor

#ifdef EVALUATE_SIZES_ALL_TOGETHER //{
    CalcSize();
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{
    CalcKernel();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}
    CalcDensityDivRot();

    InitializeRootForGravity();
    PlantGravityTree();

    HIIFlag = malloc(sizeof(bool)*Pall.Nhydro);
    HIIregions();

    for(int i=0;i<Pall.Nstars;i++){
        fprintf(stderr,"Stromgren Radius = %g [pc]\n",Pstar[i]->StromgrenRadius);
    }

    WriteStromgrenSphere(1.0);

    return ;
}
#endif //TASK_TEST_STROMGRENSPHERE

#ifdef TASK_TEST_RADIATION_PRESSURE //{
#include "RadiationPressure.h"

double nH; // in number/cc.
double Radius; // in pc.
int Nhydro; 
int Nstars;
double (*Pos)[3];

void ReadTestRadiationPressureParams(char *fname){

    if(CheckFile(fname)){
        FILE *fp;
        FileOpen(fp,fname,"r");
        fscanf(fp,"%le",&nH);
        fscanf(fp,"%le",&Radius);
        fscanf(fp,"%d",&Nhydro);
        fscanf(fp,"%d",&Nstars);
        Pos = realloc(Pos,sizeof(double)*3*(Nstars));
        for(int i=0;i<Nstars;i++){
            fscanf(fp,"%le %le %le",Pos[i],Pos[i]+1,Pos[i]+2);
        }
        fclose(fp);
    }

    return ;
}

static void WriteRadiationPressure(const int Z_in_PC){

    FILE *fp;
    char fname[MaxCharactersInLine];

    // write all
    Snprintf(fname,"./RP.all.full.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2]);
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./RP.all.full.%02d.* | sort -n > RP.all.full.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }


    // write on
    Snprintf(fname,"./RP.all.on.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(NORM(PhydroBody(i)->Vel) > 10*TINY){
            fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2]);
        }
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./RP.all.on.%02d.* | sort -n > RP.all.on.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }

    // write off
    Snprintf(fname,"./RP.all.off.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(NORM(PhydroBody(i)->Vel)<10*TINY){
            fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2]);
        }
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./RP.all.off.%02d.* | sort -n > RP.all.off.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }

    
    // slice on
    Snprintf(fname,"./RP.slice.on.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(fabs(PhydroBody(i)->Pos[2])*Pall.UnitLength/PC_CGS<Z_in_PC){
            if(NORM(PhydroBody(i)->Vel)>10*TINY){
                fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                        PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                        PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2]);
            }
        }
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./RP.slice.on.%02d.* | sort -n > RP.slice.on.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // slice off
    Snprintf(fname,"./RP.slice.off.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(fabs(PhydroBody(i)->Pos[2])*Pall.UnitLength/PC_CGS<Z_in_PC){
            if(NORM(PhydroBody(i)->Vel)<10*TINY){
            fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    Phydro[i]->MomentumRP[0],Phydro[i]->MomentumRP[1],Phydro[i]->MomentumRP[2]);
            }
        }
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./RP.slice.off.%02d.* | sort -n > RP.slice.off.%02d",MPIGetNumProcs(),MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
    }

    return ;
}

void InitTestRadiationPressure(const double nH, const double Radius, const int Nhydro, const int Nstars, double Pos[restrict][3]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
    SetViscosityParameters(0.1,1.0,1.0,0.1);

    // 
    int m = (int)log2(Nhydro);
    int Nhydro_p2 = 1<<m;
    dprintlmpi(Nhydro_p2);

    double ParticleMass = (4*M_PI/3.0)*(nH*PROTON_MASS_CGS/0.76)*CUBE(Radius*PC_CGS)
                            /MSUN_CGS/(double)Nhydro_p2;
    fprintf(stderr,"ParticleMass = %g, Total Mass = %g\n",ParticleMass,ParticleMass*Nhydro_p2);


    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    } else {
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = 0;
        Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    }
    fprintf(stderr,"[%02d] %ld %ld %ld %ld %ld | %ld %ld %ld %ld %ld\n",MPIGetMyID(),
            Pall.Nhydro,Pall.Nstars,Pall.Nsink,Pall.NDM,Pall.Ntotal,
            Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t,Pall.NDM_t,Pall.Ntotal_t);


    int AllocationSize = Pall.Ntotal; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        GenerateStructPstar(Nstars);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Nstars;i++){
            Pbody[Pall.Nhydro+i]->Baryon = (void *)(Pstar[i]);
            Pstar[i]->Body = Pbody[Pall.Nhydro+i];
        }
    }

    double eps = 0.1*PC_CGS/Pall.UnitLength;
    double Uinit = 100*Pall.ConvertTtoU;
    double cs = sqrt(Pall.GGm1*Uinit);

    int count = 0;
    int count_passall = 0;
    while(count_passall < Pall.Nhydro_t){
        double Pos[] ={Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
        double r = NORM(Pos);
        if(r<Radius){
            // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
            if(count_passall%NProcs == MyID){
                Pbody[count]->Active = ON;
                Pbody[count]->Use = ON;
                Pbody[count]->Type = TypeHydro;
                Pbody[count]->GlobalID = count_passall;

                Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
                Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
                Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];

                Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 1.e-30;
                Pbody[count]->Acc[0] = Pbody[count]->Acc[1] = Pbody[count]->Acc[2] = 0.e0;

                Pbody[count]->Eps = eps;
                Pbody[count]->Mass = ParticleMass;

                PbodyHydro(count)->Mass = Pbody[count]->Mass = ParticleMass;
                PbodyHydro(count)->Z = 0.02;
                PbodyHydro(count)->ZII = 0.02;
                PbodyHydro(count)->ZIa = 0.02;
#if (UseSFModelSpawn) 
                PbodyHydro(count)->SpawnMass = ParticleMass/((double)MaxSpawnTimes);
#endif

#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
                PbodyHydro(count)->G0 = 10;
                PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

                PbodyHydro(count)->DQheat = 0.e0;
                PbodyHydro(count)->Use = ON;
                PbodyHydro(count)->Kernel = eps;
                PbodyHydro(count)->U = Uinit;

#ifdef USE_CELIB //{
                CELibSetSolarMetallicity(Pbody[count]->Mass,PbodyHydro(count)->Elements);

                double MassLightElements = PbodyHydro(count)->Elements[CELibYield_H]
                    +PbodyHydro(count)->Elements[CELibYield_He];
                double Z = (Pbody[count]->Mass-MassLightElements)/Pbody[count]->Mass;
                PbodyHydro(count)->Z = PbodyHydro(count)->ZII = PbodyHydro(count)->ZIa = Z;
#else 
                PbodyHydro(count)->Z   = 0.02;
                PbodyHydro(count)->ZII = 0.02;
                PbodyHydro(count)->ZIa = 0.00;
#endif // USE_CELIB //}

                count ++;
            }
            count_passall ++;
        }
    }


// star particles.
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Pall.Nstars;i++){
            int index = Pall.Nhydro+i;
            Pbody[index]->Active = ON;
            Pbody[index]->Use = ON;
            Pbody[index]->Type = TypeStar;
            Pbody[index]->GlobalID = Pall.Nhydro_t+i;

            Pbody[index]->Pos[0] = Pos[i][0];
            Pbody[index]->Pos[1] = Pos[i][1];
            Pbody[index]->Pos[2] = Pos[i][2];
            Pbody[index]->Vel[0] = Pbody[index]->Vel[1] = Pbody[index]->Vel[2] = 0.e0;
            Pbody[index]->Acc[0] = Pbody[index]->Acc[1] = Pbody[index]->Acc[2] = 0.e0;
            Pbody[index]->Eps = eps;
            Pbody[index]->Mass = 100;

            PbodyStar(index)->Use = ON;
            PbodyStar(index)->ParentGlobalID = Pbody[index]->GlobalID;
            PbodyStar(index)->FormationTime = 0.e0;
#ifdef PRESERVE_SNII_EVENTRATE //{
            PbodyStar(index)->TypeIIProb = true; 
#endif // PRESERVE_SNII_EVENTRATE //}
            PbodyStar(index)->InitialMass = Pbody[index]->Mass; 
            fprintf(stderr,"Stellar Mass is %g [Msun]\n",PbodyStar(index)->InitialMass*Pall.UnitMass/MSUN_CGS);

#ifdef USE_CELIB
            PbodyStar(index)->IMFTYPE = StellarFeedbackGetIMFType();
#if 0
            PbodyStar(index)->SNIaCount = -1;
            if(PbodyStar(index)->TypeII == true){
                PbodyStar(index)->EventTime = 100*GIGAYEAR_CGS/Pall.UnitTime;
                PbodyStar(index)->SNIaCount = 100000;
            } else {
                PbodyStar(index)->EventTime = PbodyStar(index)->FormationTime 
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                            .R = gsl_rng_uniform(RandomGenerator),
                            .InitialMass_in_Msun = PbodyStar(index)->InitialMass*Pall.UnitMass/MSUN_CGS,
                            .Metallicity = PbodyStar(index)->Z,
                            .Count = 0,
                            },CELibFeedbackType_SNII)
                *YEAR_CGS/Pall.UnitTime;
            }
#endif 
#endif //USE_CELIB_AGB //}
        }
    }

    ActivateAllparticles();


    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);

#if 1
    int Nbin = 30;
    double Rad[Nbin];
    double Mass[Nbin];
    int Pnum[Nbin];
    double dr = Radius/(double)Nbin;
    gprintlmpi(dr);
    for(int i=0;i<Nbin;i++){
        Rad[i] = (i+0.5)*dr;
        Mass[i] = 0.e0;
        Pnum[i] = 0;
    }

    for(int i=0;i<Pall.Nhydro;i++){
        double r = NORM(PhydroBody(i)->Pos);
        int int_r = r/dr;
        if(int_r < 0) continue;
        if(int_r >= Nbin) continue;

        Mass[int_r] += PhydroBody(i)->Mass;
        Pnum[int_r] ++;
    }
    double GlobalMass[Nbin];
    MPI_Allreduce(Mass,GlobalMass,Nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    int GlobalPnum[Nbin];
    MPI_Allreduce(Pnum,GlobalPnum,Nbin,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,"./rad.dist","w");

        double m = 0.e0; 
        // output mean (number) density.
        for(int i=0;i<Nbin;i++){
            m += Mass[i];
            double density = m*Pall.UnitMass/((4*M_PI/3.0)*CUBE(dr*(i+1)*PC_CGS));
            double Ndensity = 0.76*density/PROTON_MASS_CGS;
            // fprintf(stderr,"Mass %g, R %g, density %g, nH %g\n",
                    // m,dr*(i+1),density,Ndensity);
            fprintf(fp,"%g %g %d %g %g\n",Rad[i],GlobalMass[i],GlobalPnum[i],density,Ndensity);
        }
        fclose(fp);
    }
#endif

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;

    SetViscosityParameters(0.1,1.0,1.0,0.1);

    Pall.TEnd = 100*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
    }


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.0001*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/RP.ASCII");
    strcpy(Pall.BaseFileName,"./data/RP");
    strcpy(Pall.RestartFileName,"./data/RP.dump");

    return;
}

int main_TestRadiationPressure(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    nH = 1.0; // in number/cc.
    Radius = 10; // in pc.
    Nhydro = 1024; 
    Nstars = 1;
    Pos = malloc(sizeof(double)*3*Nstars);
    Pos[0][0] = Pos[0][1] = Pos[0][2] = 0.e0;

    ReadTestRadiationPressureParams("RadiationPressureParams.txt");

    fprintf(stderr,"%g %g %d %d\n",nH,Radius,Nhydro,Nstars);
    fflush(NULL);

    InitTestRadiationPressure(nH,Radius,Nhydro,Nstars,Pos);

    BuildPredictors();
    InitializeDecomposition();
    DomainDecomposition();
    InitializeRootForHydro();
    InitializeCoolingTable();
    InitializeStarFormation();
    InitializeStellarFeedback();

    ASRFLXRunParameters.IMFType = CHEMICALEVOLUTION_IMFTYPE;
    ASRFLXRunParameters.UsePopIII = true;
    ASRFLXInitFUV();
    ASRFLXInitNLY();
    ASRFLXInitBol();
    InitRadiationPressure();

    ClearHydroData();
    BuildPredictors();
    PlantHydroTree();

    // Make pho predictor
    CalcSize();
    CalcDensityDivRot();
    BuildPredictors();
    // Make pho predictor

    CalcSize();
    CalcDensityDivRot();

    for(int i=0;i<Pall.Ntotal;i++)
        Pbody[i]->dt = 1*MEGAYEAR_CGS/Pall.UnitTime;
    StellarFeedback();

    for(int i=0;i<Pall.Nstars;i++){
        fprintf(stderr,"RP Radius = %g [pc]\n",Pstar[i]->RadiusRP);
    }

    WriteRadiationPressure(1.0);

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_RADIATION_PRESSURE //}

#ifdef TASK_TEST_EXACTCOOLING //{
#include "CloudyCooling.h"

static struct StructECParams{
    double InitialTemperature; // Temperature [K]
    double NumberDensity; //Number density
    double Redshift;//Redshift
    double Metallicity;//Metallicity
    double G0;//Metallicity
    double TimeStep; // Year
    double TimeStart;// Year
    double TimeEnd;  // Year
} EXParams;

static void ReadExactCoolingParams(char fname[]){

    if(CheckFile(fname)){
        FILE *fp;
        FileOpen(fp,fname,"r");
        fscanf(fp,"%le",&EXParams.InitialTemperature);
        fscanf(fp,"%le",&EXParams.NumberDensity);
        fscanf(fp,"%le",&EXParams.Redshift);
        fscanf(fp,"%le",&EXParams.Metallicity);
        fscanf(fp,"%le",&EXParams.G0);
        fscanf(fp,"%le",&EXParams.TimeStep);
        fscanf(fp,"%le",&EXParams.TimeStart);
        fscanf(fp,"%le",&EXParams.TimeEnd);
        fclose(fp);
    } else {
        // Use fiducial parameters.
        EXParams.InitialTemperature = 1.e6;
        EXParams.NumberDensity   = 0.1;
        EXParams.Redshift  = 0;
        EXParams.Metallicity  = 0.0013;
        EXParams.G0  = 1.0;
        EXParams.TimeStep  = 0.1*MEGAYEAR_CGS/YEAR_CGS;
        EXParams.TimeStart = 0;
        EXParams.TimeEnd   = 30*MEGAYEAR_CGS/YEAR_CGS;
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        char Fname[MaxCharactersInLine];
        FileOpen(fp,"Params.Log","w");

        fprintf(fp,"%g #InitialTemperature\n",EXParams.InitialTemperature);
        fprintf(fp,"%g #NumberDensity\n",EXParams.NumberDensity);
        fprintf(fp,"%g #Redshift\n",EXParams.Redshift);
        fprintf(fp,"%g #Metallicity\n",EXParams.Metallicity);
        fprintf(fp,"%g #G0\n",EXParams.G0);
        fprintf(fp,"%g #TimeStep\n",EXParams.TimeStep);
        fprintf(fp,"%g #TimeStart\n",EXParams.TimeStart);
        fprintf(fp,"%g #TimeEnd\n",EXParams.TimeEnd);

        fclose(fp);
    }

    return ;
}

static void InitTestExactCooling(void){

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = YEAR_CGS;
    Pall.UnitMass = MSUN_CGS;
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();
    return ;
}


int main_TestExactCooling(const int argc, char *argv[]){

    ReadExactCoolingParams("ExactCooling.params");

    InitTestExactCooling();

    InitializeCoolingTable();
    CloudyDataInterpolationRedshift(EXParams.Redshift);

    double Gamma_Heating = 1.e-24*(0.05)*EXParams.G0; // [erg/s]

    double LognH = log10(EXParams.NumberDensity);
    double nH = (EXParams.NumberDensity);
    double LogT = log10(EXParams.InitialTemperature);
    const double factor = CUBE(Pall.UnitTime)/SQ(Pall.UnitLength); // [erg/g/s], [cm^2/s^3]
    double Uinit = (3.0/2.0)*BOLTZMANN_CONSTANT_CGS*EXParams.InitialTemperature/(0.6*PROTON_MASS_CGS);

    FILE *fp;
    FileOpen(fp,"ExactCoolingTest.dat","w");

    double t = EXParams.TimeStart+EXParams.TimeStep;
    do{
    
        double dudt = ReturnCloudyCoolingHeatingValueExactSolverWithFUVHeating(LognH,LogT,
                        EXParams.Metallicity,t*Pall.UnitTime,Gamma_Heating/nH);
        double Unew = Uinit + dudt*(t*YEAR_CGS);
        double Tnew = Unew*(2.0/3.0)*(0.6*PROTON_MASS_CGS)/BOLTZMANN_CONSTANT_CGS;

        fprintf(stderr,"%g %g\n",t,Tnew);
        fprintf(fp,"%g %g\n",t,Tnew);
        t += EXParams.TimeStep;
    }while(t<EXParams.TimeEnd);
    fclose(fp);

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_EXACTCOOLING //}


#ifdef TASK_TEST_DOMAINDECOMPOSITION //{
static void InsertDataTestDomainDecomposition(const int TargetIndex, float Pos[], float Vel[], float Mass, float Rho, float U,
        const int Type, const int GlobalID, const float Eps){

    Pbody[TargetIndex]->GlobalID = GlobalID;
    Pbody[TargetIndex]->Pos[0] = Pos[0];
    Pbody[TargetIndex]->Pos[1] = Pos[1];
    Pbody[TargetIndex]->Pos[2] = Pos[2];

    Pbody[TargetIndex]->Vel[0] = Vel[0];
    Pbody[TargetIndex]->Vel[1] = Vel[1];
    Pbody[TargetIndex]->Vel[2] = Vel[2];
    
    Pbody[TargetIndex]->Mass = Mass;

    Pbody[TargetIndex]->Active = ON;
    Pbody[TargetIndex]->Use = ON;
    Pbody[TargetIndex]->Type = Type;
    Pbody[TargetIndex]->Eps = Eps;

    if(Type == TypeHydro){
        PbodyHydro(TargetIndex)->Use = ON;
        PbodyHydro(TargetIndex)->Rho = Rho;
        PbodyHydro(TargetIndex)->U = U;
        PbodyHydro(TargetIndex)->Kernel = 2.0*Eps;
    }
    return ;
}

static void InsertDataTestDomainDecompositionDouble(const int TargetIndex, double Pos[], double Vel[], double Mass, double Rho, double U,
        const int Type, const int GlobalID, double Eps){

    Pbody[TargetIndex]->GlobalID = GlobalID;
    Pbody[TargetIndex]->Pos[0] = Pos[0];
    Pbody[TargetIndex]->Pos[1] = Pos[1];
    Pbody[TargetIndex]->Pos[2] = Pos[2];

    Pbody[TargetIndex]->Vel[0] = Vel[0];
    Pbody[TargetIndex]->Vel[1] = Vel[1];
    Pbody[TargetIndex]->Vel[2] = Vel[2];
    
    Pbody[TargetIndex]->Mass = Mass;

    Pbody[TargetIndex]->Active = ON;
    Pbody[TargetIndex]->Use = ON;
    Pbody[TargetIndex]->Type = Type;
    Pbody[TargetIndex]->Eps = Eps;

    if(Type == TypeHydro){
        PbodyHydro(TargetIndex)->Use = ON;
        PbodyHydro(TargetIndex)->Rho = Rho;
        PbodyHydro(TargetIndex)->U = U;
        PbodyHydro(TargetIndex)->Kernel = 2.0*Eps;
    }
    return ;
}


static void InitDomainDecomposition(const int TotalNumber, const int Dimension){

    if(Dimension < 2){ 
        fprintf(stderr,"Dimension is less than 2.\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    } else if(Dimension > 3){
        fprintf(stderr,"Dimension is larger than 3.\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    
    memset(&Pall,0,sizeof(struct StructPall));
    //InitializeRandomGenerator(1977+MPIGetMyID());
    InitializeRandomGenerator(1977);

    Pall.UnitLength = 1;
    Pall.UnitTime = 1;
    Pall.UnitMass = 1;
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
    SetViscosityParameters(0.1,1.0,1.0,0.1);

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();


    int mycount = 0;
    for(int i=0;i<TotalNumber;i++){
        if(i%NProcs == MyID){
            mycount ++;
        }
    }
    Pall.Ntotal = mycount;
    Pall.Nhydro = 0;
    Pall.NDM = mycount;
    Pall.Nstars = 0;
    Pall.Ntotal_t = TotalNumber;

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
        Pbody[i]->Type = TypeHydro;
    }

    float mass = 1.0/Pall.Ntotal_t;
    float Eps = 0.01/cbrt(Pall.Ntotal_t);
    float Vel[] = {0.e0,0.e0,0.e0};
    int counter = 0;

    Pall.NDM = 0;
    for(int i=0;i<Pall.Ntotal_t;i++){
        float Pos[3] = {2.0*gsl_rng_uniform(RandomGenerator)-1.0,
                        2.0*gsl_rng_uniform(RandomGenerator)-1.0,
                        2.0*gsl_rng_uniform(RandomGenerator)-1.0};
        if(Dimension == 2){
            Pos[2] = 0.0;
        }
        if(i%NProcs == MyID){
            InsertDataTestDomainDecomposition(counter,Pos,Vel,mass,0.e0,0.e0,TypeDM,i,Eps);
            Pall.NDM ++;
            counter ++;
        }
    }
    // dlprintlmpi(Pall.NDM);

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"_Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],Pbody[i]->Mass);
    }
    fclose(fp);
    fflush(NULL);

    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK){
        char Command[MaxCharactersInLine];
        Snprintf(Command,"cat ./_Init.* | sort -n > Init.%02d",MPIGetNumProcs());
        system(Command);
        //system("cat _Init.* | sort -n > Init.dat");
        fflush(NULL);
        system("rm -rf ./_Init.*");
        fflush(NULL);
    }
#endif


    MPI_Barrier(MPI_COMM_WORLD);

    ActivateAllparticles();

    // dlprintlmpi(Pall.Ntotal);
    // dlprintlmpi(Pall.NDM);

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.TEnd = 1.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/200;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(!CheckDir("./data"))
            MakeDir("./data");
    }
    strcpy(Pall.ASCIIFileName,"./data/MB.ASCII");
    strcpy(Pall.BaseFileName,"./data/MB");
    strcpy(Pall.RestartFileName,"./data/MB.dump");

    // dlprintlmpi(Pall.Ntotal);
    // dlprintlmpi(Pall.NDM);

    return;

}

static void TestDomainDecompositionWriteDecompInfo(void){

    MakeDir("Decomp");

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./Decomp/Decomp.DM.%03d.%03d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%g %g %g %d\n",Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],Pbody[i]->Type);
    }
    fclose(fp);

    return ;
}

int main_TestDomainDecomposition(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    InitDomainDecomposition(10000,3);
    InitLogFiles();
    UpdateCosmologicalParameters();

    BuildPredictors();
    InitializeDecomposition();
    DomainDecomposition();

    TestDomainDecompositionWriteDecompInfo();

////////////////////////////////////////////////////////////////

    // PreDomainDecomposition(DomainUpdate);
    // DomainDecompositionOnlyDataExchange();

////////////////////////////////////////////////////////////////


    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Total number of time step: %d\n",Pall.TStepTotal);

    LogOutPutEnergyMomentumAngularMomentum();
    CloseLogFiles();

    return EXIT_SUCCESS;
}
#endif // TASK_TEST_DOMAINDECOMPOSITION //}

#ifdef TASK_TEST_PARTICLE_SPLITTING //{
double CloudMass; // in number/cc.
double CloudRadius; // in pc.
double CloudPower; // in pc.
double CloudTemperature; // in pc.
int    CloudParticles; 

void ReadParticleSplittingTestParams(char *fname){

    if(CheckFile(fname)){
        FILE *fp;
        FileOpen(fp,fname,"r");
        fscanf(fp,"%d",&CloudParticles);
        fscanf(fp,"%le",&CloudMass);
        fscanf(fp,"%le",&CloudRadius);
        fscanf(fp,"%le",&CloudTemperature);
        fscanf(fp,"%le",&CloudPower);
        fclose(fp);
    }

    return ;
}


static void WriteParticleSplitting(const int Counter, const double dz){

    FILE *fp;
    char fname[MaxCharactersInLine];

    if(MPIGetMyID() == MPI_ROOT_RANK)
        MakeDir("data");
    MPI_Barrier(MPI_COMM_WORLD);

    // write on
    Snprintf(fname,"./data/SP.all.%03d.%02d.%02d",Counter,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %d %d %g %g %g %g\n",PhydroBody(i)->GlobalID,
                Phydro[i]->SplitGeneration,Phydro[i]->SplitNthChildren,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                Phydro[i]->Rho*Pall.UnitMass/CUBE(Pall.UnitLength));
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./data/SP.all.%03d.%02d.* | sort -n > ./data/SP.all.%03d.%02d",
                Counter,MPIGetNumProcs(),
                Counter,MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
        Snprintf(fname,"rm  ./data/SP.all.%03d.%02d.*",Counter,MPIGetNumProcs());
        system(fname);
    }

    
    // slice on
    Snprintf(fname,"./data/SP.slice.%03d.%02d.%02d",Counter,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(fabs(PhydroBody(i)->Pos[2])*Pall.UnitLength/PC_CGS<dz)
            fprintf(fp,"%ld %d %d %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    Phydro[i]->SplitGeneration,Phydro[i]->SplitNthChildren,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    Phydro[i]->Rho*Pall.UnitMass/CUBE(Pall.UnitLength));
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Snprintf(fname,"cat ./data/SP.slice.%03d.%02d.* | sort -n > ./data/SP.slice.%03d.%02d",
                Counter,MPIGetNumProcs(),
                Counter,MPIGetNumProcs());
        fprintf(stderr,"%s\n",fname);
        system(fname);
        Snprintf(fname,"rm  ./data/SP.slice.%03d.%02d.*",Counter,MPIGetNumProcs());
        system(fname);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    
    return ;
}

static void InitTestParticleSplitting(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
    SetViscosityParameters(0.1,1.0,1.0,0.1);

    //////////////////////////////////////////////////////

    double ParticleMass = CloudMass/(double)CloudParticles;


    int mycount = 0;
    for(int i=0;i<CloudParticles;i++){
        if(i%NProcs == MyID){
            mycount ++;
        }
    }
    dprintlmpi(mycount);

    Pall.Nhydro = mycount;
    Pall.Nhydro_t = CloudParticles;
    Pall.Nstars = Pall.Nstars_t = 0;
    Pall.Nsink = Pall.Nsink_t = 0;
    Pall.Ntotal = Pall.Nhydro;
    Pall.Ntotal_t = CloudParticles;

    fprintf(stderr,"[%02d] %ld %ld %ld %ld %ld | %ld %ld %ld %ld %ld\n",MPIGetMyID(),
            Pall.Nhydro,Pall.Nstars,Pall.Nsink,Pall.NDM,Pall.Ntotal,
            Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t,Pall.NDM_t,Pall.Ntotal_t);


    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double eps = 0.01*PC_CGS/Pall.UnitLength;
    double Uinit = CloudTemperature*Pall.ConvertTtoU;
    double cs = sqrt(Pall.GGm1*Uinit);

    int count = 0;
    int count_passall = 0;
    while(count_passall < Pall.Nhydro_t){
        double Pos[] ={(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       (2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
        double r = NORM(Pos);
        if(r<1.0){
            // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
            if(count_passall%NProcs == MyID){
                Pbody[count]->Active = ON;
                Pbody[count]->Use = ON;
                Pbody[count]->Type = TypeHydro;
                Pbody[count]->GlobalID = count_passall;

                Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
                Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
                Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];

                Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 1.e-30;
                Pbody[count]->Acc[0] = Pbody[count]->Acc[1] = Pbody[count]->Acc[2] = 0.e0;

                Pbody[count]->Eps = eps;
                Pbody[count]->Mass = ParticleMass;

                PbodyHydro(count)->Mass = Pbody[count]->Mass = ParticleMass;
                PbodyHydro(count)->Z = 0.02;
                PbodyHydro(count)->ZII = 0.02;
                PbodyHydro(count)->ZIa = 0.02;
                PbodyHydro(count)->SpawnMass = ParticleMass/((double)MaxSpawnTimes);

#warning CELIB

                PbodyHydro(count)->DQheat = 0.e0;
                PbodyHydro(count)->Use = ON;
                PbodyHydro(count)->Kernel = eps;
                PbodyHydro(count)->U = Uinit;

                count ++;
            }
            count_passall ++;
        }
    }

    // Stretch positions to be fit the power law index.


    if(CloudParticles > 0){ // change profile.
// rho=const, dv = r^2*dr
//  M \propto rho*dv = r^3
// rho = r^n, dv = r^2*dr
//  M \propto rho*dv = r^{n+3}

        for(int i=0;i<Pall.Nhydro;i++){

            double Pos[] ={PhydroBody(i)->Pos[0],
                           PhydroBody(i)->Pos[1],
                           PhydroBody(i)->Pos[2]};

            double r = NORM(Pos);

            if(r>TINY){
                Pos[0] /= r;
                Pos[1] /= r;
                Pos[2] /= r;
            }
            /*
            Pos[0] *= (r)*sqrt(r);
            Pos[1] *= (r)*sqrt(r);
            Pos[2] *= (r)*sqrt(r);
            */
            Pos[0] *= SQ(r)*sqrt(r);
            Pos[1] *= SQ(r)*sqrt(r);
            Pos[2] *= SQ(r)*sqrt(r);


            PhydroBody(i)->Pos[0] = PhydroBody(i)->PosP[0] = Phydro[i]->PosP[0] = Pos[0];
            PhydroBody(i)->Pos[1] = PhydroBody(i)->PosP[1] = Phydro[i]->PosP[1] = Pos[1];
            PhydroBody(i)->Pos[2] = PhydroBody(i)->PosP[2] = Phydro[i]->PosP[2] = Pos[2];
        }
    }

    ActivateAllparticles();


    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);

#if 0
    int Nbin = 30;
    double Rad[Nbin];
    double Mass[Nbin];
    int Pnum[Nbin];
    double dr = Radius/(double)Nbin;
    gprintlmpi(dr);
    for(int i=0;i<Nbin;i++){
        Rad[i] = (i+0.5)*dr;
        Mass[i] = 0.e0;
        Pnum[i] = 0;
    }

    for(int i=0;i<Pall.Nhydro;i++){
        double r = NORM(PhydroBody(i)->Pos);
        int int_r = r/dr;
        if(int_r < 0) continue;
        if(int_r >= Nbin) continue;

        Mass[int_r] += PhydroBody(i)->Mass;
        Pnum[int_r] ++;
    }
    double GlobalMass[Nbin];
    MPI_Allreduce(Mass,GlobalMass,Nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    int GlobalPnum[Nbin];
    MPI_Allreduce(Pnum,GlobalPnum,Nbin,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,"./rad.dist","w");

        double m = 0.e0; 
        // output mean (number) density.
        for(int i=0;i<Nbin;i++){
            m += Mass[i];
            double density = m*Pall.UnitMass/((4*M_PI/3.0)*CUBE(dr*(i+1)*PC_CGS));
            double Ndensity = 0.76*density/PROTON_MASS_CGS;
            // fprintf(stderr,"Mass %g, R %g, density %g, nH %g\n",
                    // m,dr*(i+1),density,Ndensity);
            fprintf(fp,"%g %g %d %g %g\n",Rad[i],GlobalMass[i],GlobalPnum[i],density,Ndensity);
        }
        fclose(fp);
    }
#endif

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.TEnd = 1*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
    }

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.01*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/Sp.ASCII");
    strcpy(Pall.BaseFileName,"./data/Sp");
    strcpy(Pall.RestartFileName,"./data/Sp.dump");

    return;
}

static void UpParticleMass(void){
    
    double UpFactor = 3;
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->Mass *= UpFactor;
        PhydroBody(i)->Mass *= UpFactor;
        for(int k=0;k<CELibYield_Number;k++){
            Phydro[i]->Elements[k] *= UpFactor;
        }
    }
    return ;
}

int main_Test_ParticleSplitting(const int argc, char *argv[]){

    InitializeCommunicationTable();
    InitializeCommunicationBuffers();
    GetRunStatus(argc,argv);

    CloudMass   = 1.0; // in Msun.
    CloudRadius = 1.0;      // in pc.
    CloudPower  = -3.0;      // 
    CloudTemperature  = 1.e+4;      // 
    CloudParticles = 10000; // 

    ReadParticleSplittingTestParams("ParticleSplittingTestParams.txt");

    InitTestParticleSplitting();

    BuildPredictors();
    InitializeDecomposition();
    DomainDecomposition();
    InitializeRootForHydro();

    // loop changing gas temperature and write data

    for(int i=0;i<10;i++){
        // up particle mass;
        UpParticleMass();


        ClearHydroData();
        BuildPredictors();
        PlantHydroTree();

        CalcSize();
        CalcDensityDivRot();
        ParticleSplitting();

        WriteParticleSplitting(i,0.1);
    }

    return EXIT_SUCCESS;
}
#endif //TASK_TEST_PARTICLE_SPLITTING
