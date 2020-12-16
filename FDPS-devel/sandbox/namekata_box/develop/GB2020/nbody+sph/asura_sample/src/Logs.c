#include "config.h"

void InitLogFiles(void){

    if(MPIGetMyID()==MPI_ROOT_RANK){

        MakeDir("./Log");

        sprintf(FnameEnergy,"./Log/Energy.data");
        sprintf(FnameMomentum,"./Log/Momentum.data");
        sprintf(FnameAMomentum,"./Log/AMomentum.data");
        sprintf(FnameElapsedTime,"./Log/ElapsedTime.data");
        sprintf(FnameStepLog,"./Log/StepLog.data");
        if(Pall.RunStatus == NewSimulation){
            FileOpen(FpEnergy,FnameEnergy,"w");
            FileOpen(FpMomentum,FnameMomentum,"w");
            FileOpen(FpAMomentum,FnameAMomentum,"w");
            FileOpen(FpElapsedTime,FnameElapsedTime,"w");
        }else if(Pall.RunStatus == RestartSimulation){
            FileOpen(FpEnergy,FnameEnergy,"a");
            FileOpen(FpMomentum,FnameMomentum,"a");
            FileOpen(FpAMomentum,FnameAMomentum,"a");
            FileOpen(FpElapsedTime,FnameElapsedTime,"a");
        }
        FileOpen(FpStepLog,FnameStepLog,"a");

    }

    return;
}

void CloseLogFiles(void){

    if(MPIGetMyID()==MPI_ROOT_RANK){
        fclose(FpEnergy);
        fclose(FpMomentum);
        fclose(FpAMomentum);
        fclose(FpElapsedTime);
        fclose(FpStepLog);
    }

    return;
}

static inline __attribute__((always_inline)) void *SendBuf(const bool target, void *sendbuf){
    return target
        ? MPI_IN_PLACE
        : sendbuf;
}

void LogOutPutEnergyMomentumAngularMomentum(void){

    struct LogPhysicalQuantities{
        double Ek,Ep,Ethermal,Eloss;    // Ek,Ep,Ethermal,Ecoolingloss
        double Momentum[3];     // Mx,My,Mz, not specific
        double AMomentum[3];    // AMx,AMy,AMz, not specific
        double Mass;            // Total mass
    } PhysicalQuantities = {0.e0};

    // calc local
    for(int i=0;i<Pall.Ntotal;i++){
        PhysicalQuantities.Ek += Pbody[i]->Mass*NORM2(Pbody[i]->Vel);
        PhysicalQuantities.Ep += Pbody[i]->Pot;
        if(Pbody[i]->Type == TypeHydro){
            PhysicalQuantities.Ethermal += Pbody[i]->Mass*PbodyHydroU(i);
        }

        PhysicalQuantities.Momentum[0] += Pbody[i]->Mass*Pbody[i]->Vel[0];
        PhysicalQuantities.Momentum[1] += Pbody[i]->Mass*Pbody[i]->Vel[1];
        PhysicalQuantities.Momentum[2] += Pbody[i]->Mass*Pbody[i]->Vel[2];

        PhysicalQuantities.AMomentum[0] += 
            Pbody[i]->Mass*(Pbody[i]->Pos[1]*Pbody[i]->Vel[2]-Pbody[i]->Pos[2]*Pbody[i]->Vel[1]);
        PhysicalQuantities.AMomentum[1] += 
            Pbody[i]->Mass*(Pbody[i]->Pos[2]*Pbody[i]->Vel[0]-Pbody[i]->Pos[0]*Pbody[i]->Vel[2]);
        PhysicalQuantities.AMomentum[2] += 
            Pbody[i]->Mass*(Pbody[i]->Pos[0]*Pbody[i]->Vel[1]-Pbody[i]->Pos[1]*Pbody[i]->Vel[0]);

        PhysicalQuantities.Mass += Pbody[i]->Mass;
    }
    PhysicalQuantities.Eloss = Pall.CoolingEnergyLoss;

    const int MyID = MPIGetMyID();
    MPI_Reduce(SendBuf(MyID==MPI_ROOT_RANK,&PhysicalQuantities),&PhysicalQuantities,11,MPI_DOUBLE,MPI_SUM,MPI_ROOT_RANK,MPI_COMM_WORLD);

    // IO
    if(MyID == MPI_ROOT_RANK){
        fprintf(FpEnergy,"%e %e %e %e %e %e %e\n",Pall.TCurrent,Pall.Redshift,
                0.5*PhysicalQuantities.Ek,
                PhysicalQuantities.Ep,
                PhysicalQuantities.Ethermal,
                PhysicalQuantities.Eloss,
                0.5*PhysicalQuantities.Ek+PhysicalQuantities.Ep
                 +PhysicalQuantities.Ethermal+PhysicalQuantities.Eloss);

        fprintf(FpMomentum,"%e %e %e %e %e %e %e %e\n",Pall.TCurrent,Pall.Redshift,
                PhysicalQuantities.Momentum[0],PhysicalQuantities.Momentum[1],
                PhysicalQuantities.Momentum[2],
                PhysicalQuantities.Momentum[0]/PhysicalQuantities.Mass,
                    PhysicalQuantities.Momentum[1]/PhysicalQuantities.Mass,
                        PhysicalQuantities.Momentum[2]/PhysicalQuantities.Mass);

        fprintf(FpAMomentum,"%e %e %e %e %e %e %e %e\n",Pall.TCurrent,Pall.Redshift,
                PhysicalQuantities.AMomentum[0],PhysicalQuantities.AMomentum[1],
                PhysicalQuantities.AMomentum[2],
                PhysicalQuantities.AMomentum[0]/PhysicalQuantities.Mass,
                    PhysicalQuantities.AMomentum[1]/PhysicalQuantities.Mass,
                        PhysicalQuantities.AMomentum[2]/PhysicalQuantities.Mass);
    }

    fflush(NULL);
    return;
}

void LogOutPutEnergyMomentumAngularMomentumPredictors(void){

    struct LogPhysicalQuantities{
        double Ek,Ep,Ethermal,Eloss;    // Ek,Ep,Ethermal,Ecoolingloss
        double Momentum[3];     // Mx,My,Mz, not specific
        double AMomentum[3];    // AMx,AMy,AMz, not specific
        double Mass;            // Total mass
    } PhysicalQuantities = {0.e0};

    // calc local
    for(int i=0;i<Pall.Ntotal;i++){
        PhysicalQuantities.Ek += Pbody[i]->Mass*NORM2(Pbody[i]->Vel);
        PhysicalQuantities.Ep += Pbody[i]->Pot;
        if(Pbody[i]->Type == TypeHydro)
            PhysicalQuantities.Ethermal += Pbody[i]->Mass*PbodyHydroUPred(i);

        PhysicalQuantities.Momentum[0] += Pbody[i]->Mass*Pbody[i]->Vel[0];
        PhysicalQuantities.Momentum[1] += Pbody[i]->Mass*Pbody[i]->Vel[1];
        PhysicalQuantities.Momentum[2] += Pbody[i]->Mass*Pbody[i]->Vel[2];

        PhysicalQuantities.AMomentum[0] += 
            Pbody[i]->Mass*(Pbody[i]->PosP[1]*Pbody[i]->Vel[2]-Pbody[i]->PosP[2]*Pbody[i]->Vel[1]);
        PhysicalQuantities.AMomentum[1] += 
            Pbody[i]->Mass*(Pbody[i]->PosP[2]*Pbody[i]->Vel[0]-Pbody[i]->PosP[0]*Pbody[i]->Vel[2]);
        PhysicalQuantities.AMomentum[2] += 
            Pbody[i]->Mass*(Pbody[i]->PosP[0]*Pbody[i]->Vel[1]-Pbody[i]->PosP[1]*Pbody[i]->Vel[0]);

        PhysicalQuantities.Mass += Pbody[i]->Mass;
    }
    PhysicalQuantities.Eloss = Pall.CoolingEnergyLoss;

    const int MyID = MPIGetMyID();
    MPI_Reduce(SendBuf(MyID==MPI_ROOT_RANK,&PhysicalQuantities),&PhysicalQuantities,11,MPI_DOUBLE,MPI_SUM,MPI_ROOT_RANK,MPI_COMM_WORLD);

    // IO
    if(MyID == MPI_ROOT_RANK){
        fprintf(FpEnergy,"%e %e %e %e %e %e %e\n",Pall.TCurrent,Pall.Redshift,
                0.5*PhysicalQuantities.Ek,
                PhysicalQuantities.Ep,
                PhysicalQuantities.Ethermal,
                PhysicalQuantities.Eloss,
                0.5*PhysicalQuantities.Ek+PhysicalQuantities.Ep
                 +PhysicalQuantities.Ethermal+PhysicalQuantities.Eloss);

        fprintf(FpMomentum,"%e %e %e %e %e %e %e %e\n",Pall.TCurrent,Pall.Redshift,
                PhysicalQuantities.Momentum[0],PhysicalQuantities.Momentum[1],
                PhysicalQuantities.Momentum[2],
                PhysicalQuantities.Momentum[0]/PhysicalQuantities.Mass,
                    PhysicalQuantities.Momentum[1]/PhysicalQuantities.Mass,
                        PhysicalQuantities.Momentum[2]/PhysicalQuantities.Mass);

        fprintf(FpAMomentum,"%e %e %e %e %e %e %e %e\n",Pall.TCurrent,Pall.Redshift,
                PhysicalQuantities.AMomentum[0],PhysicalQuantities.AMomentum[1],
                PhysicalQuantities.AMomentum[2],
                PhysicalQuantities.AMomentum[0]/PhysicalQuantities.Mass,
                    PhysicalQuantities.AMomentum[1]/PhysicalQuantities.Mass,
                        PhysicalQuantities.AMomentum[2]/PhysicalQuantities.Mass);
    }


    fflush(NULL);
    return;
}

void LogOutPutElapsedTime(void){

    // IO
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(FpElapsedTime,"%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            Pall.TStepTotal,Pall.TCurrent,Pall.Redshift,
            TimingResults.Total,                                            // <4>
            TimingResults.Decomposition,                                    // <5>
            TimingResults.Gravity,                                          // <6>
            TimingResults.GravityTree,                                      // <7>
            TimingResults.GravityComm,    // A part of TimeGravity          // <8>
            TimingResults.Hydro,                                            // <9>
            TimingResults.HydroNeighborSearch,                              // <10>
            TimingResults.HydroTree,                                        // <11>
            TimingResults.HydroComm,                                        // <12>
            TimingResults.HydroKernel,                                      // <13>
            TimingResults.HydroKernelNeighborSearch,                        // <14>
            TimingResults.HydroKernelComm,                                  // <15>
            TimingResults.HydroDensity,                                     // <16>
            TimingResults.HydroDensityNeighborSearch,                       // <17>
            TimingResults.HydroDensityComm,                                 // <18>
            TimingResults.HydroAcc,                                         // <19>
            TimingResults.HydroAccNeighborSearch,                           // <20>
            TimingResults.HydroAccComm,                                     // <21>
            TimingResults.Integral,                                         // <22>
            TimingResults.Cooling,                                          // <23>
            TimingResults.GravImbarance,                                    // <24>
            TimingResults.HydroImbarance,                                   // <25>
            TimingResults.Starformation);                                   // <26>

    return;
}

void EndLogs(void){

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"\n");
        fprintf(stderr,"------------------------------\n");
        //fprintf(stderr,"This run is finished.\n");
        fprintf(stderr,"Since the start of this run...\n");
        fprintf(stderr,"Number of steps = %d\n",Pall.TStepTotal);
#ifdef USE_SINK_PARTICLE //{
        fprintf(stderr,"Number of particles (total,dm,hydro,stars,sinks) = (%ld,%ld,%ld,%ld,%ld)\n",
                Pall.Ntotal_t,Pall.NDM_t,Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t);
        fprintf(stderr,"Update Number (total,dm,hydro,stars,sinks) = (%ld,%ld,%ld,%ld,%ld)\n",
                Pall.TActiveParticles,Pall.TActiveParticlesDM,
                Pall.TActiveParticlesHydro,Pall.TActiveParticlesStar,
                Pall.TActiveParticlesSink);
#else  // USE_SINK_PARTICLE //}//{
        fprintf(stderr,"Number of particles (total,dm,hydro,stars) = (%ld,%ld,%ld,%ld)\n",
                Pall.Ntotal_t,Pall.NDM_t,Pall.Nhydro_t,Pall.Nstars_t);
        fprintf(stderr,"Update Number (total,dm,hydro,stars) = (%ld,%ld,%ld,%ld)\n",
                Pall.TActiveParticles,Pall.TActiveParticlesDM,
                Pall.TActiveParticlesHydro,Pall.TActiveParticlesStar);
#endif // USE_SINK_PARTICLE //}

        fprintf(stderr,"\n");
        fprintf(stderr,"Total Elapsed time: %g [sec]\n",TimingResults.Total);
        fprintf(stderr,"\n");
        fprintf(stderr,"-Total Gravity: %g [sec]\n",TimingResults.Gravity+TimingResults.GravityTree);
        if(Pall.Nhydro_t>0)
            fprintf(stderr,"-Total Hydro: %g [sec]\n",TimingResults.HydroKernel+TimingResults.HydroDensity+
                TimingResults.HydroAcc+TimingResults.Feedback+TimingResults.Sink+TimingResults.HydroTree+
                TimingResults.Cooling+TimingResults.HydroNeighborSearch);
        fprintf(stderr,"-Total Others: %g [sec]\n",TimingResults.Decomposition+
                TimingResults.Integral+TimingResults.KeyGeneration
                +TimingResults.SortStructures+TimingResults.TimeStep);
        fprintf(stderr,"\n");
        fprintf(stderr,"=Gravity= GetForce: %g [sec], Tree: %g [sec], Comm: %g [sec]",
                TimingResults.Gravity,TimingResults.GravityTree,TimingResults.GravityComm);
        if(Pall.Nhydro_t>0){
            fprintf(stderr,"\n");
            fprintf(stderr,"=Hydro= Interaction: %g [sec], Tree: %g [sec], NBS: %g [sec], Comm = %g [sec]\n",
                    TimingResults.Hydro,TimingResults.HydroTree,
                        TimingResults.HydroNeighborSearch,TimingResults.HydroComm);
            fprintf(stderr,"  =HydroKernel= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    TimingResults.HydroKernel,
                    TimingResults.HydroKernelNeighborSearch,TimingResults.HydroKernelComm);
            fprintf(stderr,"  =HydroDensity= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    TimingResults.HydroDensity,
                        TimingResults.HydroDensityNeighborSearch,TimingResults.HydroDensityComm);
            fprintf(stderr,"  =HydroAcc= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    TimingResults.HydroAcc,
                        TimingResults.HydroAccNeighborSearch,TimingResults.HydroAccComm);
#ifdef COOLING_RUN //{
            fprintf(stderr,"  =Cooling= %g [sec]\n",TimingResults.Cooling);
#endif // COOLING_RUN //}
#ifdef STARFORMATION //{
        fprintf(stderr,"  =Starformation= : %g [sec]\n",TimingResults.Starformation);
#endif // STARFORMATION //}
#ifdef DELAYED_FEEDBACK //}
        fprintf(stderr,"  =Feedback= Interaction: %g [sec], NBS: %g [sec], Comm = %g [sec]\n",
            TimingResults.Feedback,TimingResults.FeedbackNeighborSearch,TimingResults.FeedbackComm);
#endif // DELAYED_FEEDBACK //}
#ifdef USE_SINK_PARTICLE //{
        fprintf(stderr,"  =Sink particles= Sink formation: %g [sec], Accretion: %g [sec], Comm = %g [sec]\n",
            TimingResults.SinkFormation,TimingResults.SinkAccretion,TimingResults.SinkComm);
#endif // USE_SINK_PARTICLE //}
        }
        fprintf(stderr,"\n");
        fprintf(stderr,"Decomposition: %g [sec]\n",TimingResults.Decomposition);
        fprintf(stderr,"Time Integral : %g [sec]\n",TimingResults.Integral);
        fprintf(stderr,"Key generation : %g [sec]\n",TimingResults.KeyGeneration);
        fprintf(stderr,"Particle sort : %g [sec]\n",TimingResults.SortStructures);
        fprintf(stderr,"Time Step : %g [sec]\n",TimingResults.TimeStep);
        fprintf(stderr,"\n");
        //fprintf(stderr,"=Imbalances=\n");
        //fprintf(stderr,"Gravity Imbalance = %g\n",TimingResults.GravImbarance);
        //if(Pall.Nhydro_t)
            //fprintf(stderr,"Hydro Imbalance = %g\n",TimingResults.HydroImbarance);
        fprintf(stderr,"------------------------------\n");
    }

    return;
}

void LogsThisTimeStep(void){

    //int NProcs = MPIGetNumProcs();
    //double GravityImbarance;
    //double HydroImbarance;
    //MPI_Reduce(&TimingResults.GravImbaranceThisStep,&GravityImbarance,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    //MPI_Reduce(&TimingResults.HydroImbaranceThisStep,&HydroImbarance,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    //TimingResults.GravImbaranceThisStep /= NProcs;
    //TimingResults.HydroImbaranceThisStep /= NProcs;
#if 0

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"\n");
        fprintf(stderr,"------------------------------\n");
        fprintf(stderr,"Timing results of this step.\n");
        fprintf(stderr,"Current step number = %d\n",Pall.TStepTotal);
        fprintf(stderr,"Number of particles (total,dm,hydro,stars) = (%ld,%ld,%ld,%ld)\n",
                Pall.NActives_t,Pall.NActivesDM_t,Pall.NActivesHydro_t,Pall.NActivesStars_t);
        fprintf(stderr,"Update Number until now (total,dm,hydro,stars) = (%ld,%ld,%ld,%ld)\n",
                Pall.TActiveParticles,Pall.TActiveParticlesDM,
                Pall.TActiveParticlesHydro,Pall.TActiveParticlesStar);

        fprintf(stderr,"\n");
        fprintf(stderr,"Total Elapsed time: %g [sec]\n",TimingResults.TotalThisStep);
        fprintf(stderr,"-Total Gravity: %g [sec]\n",TimingResults.GravityThisStep+TimingResults.GravityTreeThisStep);
        if(Pall.Nhydro_t>0)
            fprintf(stderr,"-Total Hydro: %g [sec]\n",TimingResults.HydroKernelThisStep+TimingResults.HydroDensityThisStep+
                TimingResults.HydroAccThisStep+TimingResults.FeedbackThisStep+
                TimingResults.HydroTreeThisStep+TimingResults.CoolingThisStep+TimingResults.HydroNeighborSearchThisStep);
        fprintf(stderr,"-Total Others: %g [sec]\n",TimingResults.DecompositionThisStep+
                TimingResults.IntegralThisStep+TimingResults.KeyGenerationThisStep
                +TimingResults.SortStructuresThisStep+TimingResults.TimeStepThisStep);
        fprintf(stderr,"\n");
        fprintf(stderr,"=Gravity= GetForce: %g [sec], Tree: %g [sec], Comm: %g [sec]",
                TimingResults.GravityThisStep-TimingResults.GravityCommThisStep,
                TimingResults.GravityTreeThisStep,TimingResults.GravityCommThisStep);

        if(Pall.Nhydro_t>0){
            fprintf(stderr,"\n");
            double HydroKernel = TimingResults.HydroKernelThisStep
                                -TimingResults.HydroKernelNeighborSearchThisStep
                                -TimingResults.HydroKernelCommThisStep;
            double HydroDensity = TimingResults.HydroDensityThisStep
                                 -TimingResults.HydroDensityNeighborSearchThisStep
                                 -TimingResults.HydroDensityCommThisStep;
            double HydroAcc = TimingResults.HydroAccThisStep
                             -TimingResults.HydroAccNeighborSearchThisStep
                             -TimingResults.HydroAccCommThisStep;
            double nbs = TimingResults.HydroKernelNeighborSearchThisStep
                        +TimingResults.HydroDensityNeighborSearchThisStep
                        +TimingResults.HydroAccNeighborSearchThisStep;
            double comm = TimingResults.HydroKernelCommThisStep
                        +TimingResults.HydroDensityCommThisStep
                        +TimingResults.HydroAccCommThisStep;
            fprintf(stderr,"=Hydro= Interaction: %g [sec], Tree: %g [sec], NBS: %g [sec], Comm = %g [sec]\n",
                    HydroKernel+HydroDensity+HydroAcc,TimingResults.HydroTreeThisStep,nbs,comm);
            fprintf(stderr,"  =HydroKernel= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    HydroKernel,TimingResults.HydroKernelNeighborSearchThisStep,TimingResults.HydroKernelCommThisStep);
            fprintf(stderr,"  =HydroDensity= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    HydroDensity,TimingResults.HydroDensityNeighborSearchThisStep,TimingResults.HydroDensityCommThisStep);
            fprintf(stderr,"  =HydroAcc= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    HydroAcc,TimingResults.HydroAccNeighborSearchThisStep,TimingResults.HydroAccCommThisStep);
#ifdef COOLING_RUN
            fprintf(stderr,"  =Cooling= %g [sec]\n",TimingResults.CoolingThisStep);
#endif
        }
#ifdef STARFORMATION
        fprintf(stderr,"  =Starformation= : %g [sec]\n",TimingResults.StarformationThisStep);
#endif
#ifdef DELAYED_FEEDBACK
        fprintf(stderr,"  =Feedback= Interaction: %g [sec], NBS: %g [sec], Comm = %g [sec]\n",
            TimingResults.FeedbackThisStep-TimingResults.FeedbackNeighborSearchThisStep
                                          -TimingResults.FeedbackCommThisStep,
                TimingResults.FeedbackNeighborSearchThisStep,TimingResults.FeedbackCommThisStep);
#endif
        fprintf(stderr,"\n");
        fprintf(stderr,"Decomposition: %g [sec]\n",TimingResults.DecompositionThisStep);
        fprintf(stderr,"Time Integral: %g [sec]\n",TimingResults.IntegralThisStep);
        fprintf(stderr,"Key generation: %g [sec]\n",TimingResults.KeyGenerationThisStep);
        fprintf(stderr,"Particle sort: %g [sec]\n",TimingResults.SortStructuresThisStep);
        fprintf(stderr,"Time Step: %g [sec]\n",TimingResults.TimeStepThisStep);
        //fprintf(stderr,"\n");
        //fprintf(stderr,"=Imbalances=\n");
        //fprintf(stderr,"Gravity Imbalance = %g\n",TimingResults.GravImbaranceThisStep);
        //if(Pall.Nhydro_t)
            //fprintf(stderr,"Hydro Imbalance = %g\n",TimingResults.HydroImbaranceThisStep);
        fprintf(stderr,"------------------------------\n");
        fprintf(stderr,"\n");
    }
#else

    for(int i=0;i<MPIGetNumProcs();i++){

    if(MPIGetMyID() == i){
        fprintf(stderr,"[%d]\n",MPIGetMyID());
        fprintf(stderr,"------------------------------\n");
        fprintf(stderr,"Timing results of this step.\n");
        fprintf(stderr,"Current step number = %d\n",Pall.TStepTotal);
        fprintf(stderr,"Number of particles (total,dm,hydro,stars) = (%ld,%ld,%ld,%ld)\n",
                Pall.NActives_t,Pall.NActivesDM_t,Pall.NActivesHydro_t,Pall.NActivesStars_t);
        fprintf(stderr,"Update Number until now (total,dm,hydro,stars) = (%ld,%ld,%ld,%ld)\n",
                Pall.TActiveParticles,Pall.TActiveParticlesDM,
                Pall.TActiveParticlesHydro,Pall.TActiveParticlesStar);

        fprintf(stderr,"\n");
        fprintf(stderr,"Total Elapsed time: %g [sec]\n",TimingResults.TotalThisStep);
        fprintf(stderr,"-Total Gravity: %g [sec]\n",TimingResults.GravityThisStep+TimingResults.GravityTreeThisStep);
        if(Pall.Nhydro_t>0)
            fprintf(stderr,"-Total Hydro: %g [sec]\n",TimingResults.HydroKernelThisStep+TimingResults.HydroDensityThisStep+
                TimingResults.HydroAccThisStep+TimingResults.FeedbackThisStep+
                TimingResults.HydroTreeThisStep+TimingResults.CoolingThisStep+TimingResults.HydroNeighborSearchThisStep);
        fprintf(stderr,"-Total Others: %g [sec]\n",TimingResults.DecompositionThisStep+
                TimingResults.IntegralThisStep+TimingResults.KeyGenerationThisStep
                +TimingResults.SortStructuresThisStep+TimingResults.TimeStepThisStep);
        fprintf(stderr,"\n");
        fprintf(stderr,"=Gravity= GetForce: %g [sec], Tree: %g [sec], Comm: %g [sec]",
                TimingResults.GravityThisStep-TimingResults.GravityCommThisStep,
                TimingResults.GravityTreeThisStep,TimingResults.GravityCommThisStep);

        if(Pall.Nhydro_t>0){
            fprintf(stderr,"\n");
            double HydroKernel = TimingResults.HydroKernelThisStep
                                -TimingResults.HydroKernelNeighborSearchThisStep
                                -TimingResults.HydroKernelCommThisStep;
            double HydroDensity = TimingResults.HydroDensityThisStep
                                 -TimingResults.HydroDensityNeighborSearchThisStep
                                 -TimingResults.HydroDensityCommThisStep;
            double HydroAcc = TimingResults.HydroAccThisStep
                             -TimingResults.HydroAccNeighborSearchThisStep
                             -TimingResults.HydroAccCommThisStep;
            double nbs = TimingResults.HydroKernelNeighborSearchThisStep
                        +TimingResults.HydroDensityNeighborSearchThisStep
                        +TimingResults.HydroAccNeighborSearchThisStep;
            double comm = TimingResults.HydroKernelCommThisStep
                        +TimingResults.HydroDensityCommThisStep
                        +TimingResults.HydroAccCommThisStep;
            fprintf(stderr,"=Hydro= Interaction: %g [sec], Tree: %g [sec], NBS: %g [sec], Comm = %g [sec]\n",
                    HydroKernel+HydroDensity+HydroAcc,TimingResults.HydroTreeThisStep,nbs,comm);
            fprintf(stderr,"  =HydroKernel= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    HydroKernel,TimingResults.HydroKernelNeighborSearchThisStep,TimingResults.HydroKernelCommThisStep);
            fprintf(stderr,"  =HydroDensity= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    HydroDensity,TimingResults.HydroDensityNeighborSearchThisStep,TimingResults.HydroDensityCommThisStep);
            fprintf(stderr,"  =HydroAcc= Interaction: %g [sec], NBS: %g [sec], Comm: %g [sec]\n",
                    HydroAcc,TimingResults.HydroAccNeighborSearchThisStep,TimingResults.HydroAccCommThisStep);
#ifdef COOLING_RUN
            fprintf(stderr,"  =Cooling= %g [sec]\n",TimingResults.CoolingThisStep);
#endif
        }
#ifdef STARFORMATION
        fprintf(stderr,"  =Starformation= : %g [sec]\n",TimingResults.StarformationThisStep);
#endif
#ifdef DELAYED_FEEDBACK
        fprintf(stderr,"  =Feedback= Interaction: %g [sec], NBS: %g [sec], Comm = %g [sec]\n",
            TimingResults.FeedbackThisStep-TimingResults.FeedbackNeighborSearchThisStep
                                          -TimingResults.FeedbackCommThisStep,
                TimingResults.FeedbackNeighborSearchThisStep,TimingResults.FeedbackCommThisStep);
#endif
        fprintf(stderr,"\n");
        fprintf(stderr,"Decomposition: %g [sec]\n",TimingResults.DecompositionThisStep);
        fprintf(stderr,"Time Integral: %g [sec]\n",TimingResults.IntegralThisStep);
        fprintf(stderr,"Key generation: %g [sec]\n",TimingResults.KeyGenerationThisStep);
        fprintf(stderr,"Particle sort: %g [sec]\n",TimingResults.SortStructuresThisStep);
        fprintf(stderr,"Time Step: %g [sec]\n",TimingResults.TimeStepThisStep);
        //fprintf(stderr,"\n");
        //fprintf(stderr,"=Imbalances=\n");
        //fprintf(stderr,"Gravity Imbalance = %g\n",TimingResults.GravImbaranceThisStep);
        //if(Pall.Nhydro_t)
            //fprintf(stderr,"Hydro Imbalance = %g\n",TimingResults.HydroImbaranceThisStep);
        fprintf(stderr,"------------------------------\n");
        fprintf(stderr,"\n");
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    return;
}

void ClearTimingLogsThisStep(void){

    TimingResults.TotalThisStep = 0.e0;

    TimingResults.DecompositionThisStep = 0.e0;

    TimingResults.GravityThisStep = 0.e0;
    TimingResults.GravityTreeThisStep = 0.e0;
    TimingResults.GravityCommThisStep = 0.e0;

    TimingResults.HydroThisStep = 0.e0;
    TimingResults.HydroNeighborSearchThisStep = 0.e0;
    TimingResults.HydroTreeThisStep = 0.e0;
    TimingResults.HydroCommThisStep = 0.e0;

    TimingResults.HydroKernelThisStep = 0.e0;
    TimingResults.HydroKernelNeighborSearchThisStep = 0.e0;
    TimingResults.HydroKernelCommThisStep = 0.e0;

    TimingResults.HydroDensityThisStep = 0.e0;
    TimingResults.HydroDensityNeighborSearchThisStep = 0.e0;
    TimingResults.HydroDensityCommThisStep = 0.e0;

    TimingResults.HydroAccThisStep = 0.e0;
    TimingResults.HydroAccNeighborSearchThisStep = 0.e0;
    TimingResults.HydroAccCommThisStep = 0.e0;
   
    TimingResults.CoolingThisStep = 0.e0;

    // feedback
    TimingResults.FeedbackThisStep = 0.e0;
    TimingResults.FeedbackNeighborSearchThisStep = 0.e0;
    TimingResults.FeedbackCommThisStep = 0.e0;

    TimingResults.StarformationThisStep = 0.e0;
    
    // others
    TimingResults.IntegralThisStep = 0.e0;
    TimingResults.SortStructuresThisStep = 0.e0;
    TimingResults.KeyGenerationThisStep = 0.e0;
    TimingResults.TimeStepThisStep = 0.e0;

    return;
}

void UpdateTimeLogs(void){

    TimingResults.Total += TimingResults.TotalThisStep;

    TimingResults.Decomposition += TimingResults.DecompositionThisStep;   

    // gravity
    TimingResults.Gravity += TimingResults.GravityThisStep;
    TimingResults.GravityTree += TimingResults.GravityTreeThisStep;
    TimingResults.GravityComm += TimingResults.GravityCommThisStep;

    // hydro
    double HydroKernel = (TimingResults.HydroKernelThisStep
                         -TimingResults.HydroKernelNeighborSearchThisStep
                         -TimingResults.HydroKernelCommThisStep);
    TimingResults.HydroKernel += HydroKernel;
    TimingResults.HydroKernelNeighborSearch += TimingResults.HydroKernelNeighborSearchThisStep;
    TimingResults.HydroKernelComm += TimingResults.HydroKernelCommThisStep;

    double HydroDensity = (TimingResults.HydroDensityThisStep
                          -TimingResults.HydroDensityNeighborSearchThisStep
                          -TimingResults.HydroDensityCommThisStep);
    TimingResults.HydroDensity += HydroDensity;
    TimingResults.HydroDensityNeighborSearch += TimingResults.HydroDensityNeighborSearchThisStep;
    TimingResults.HydroDensityComm += TimingResults.HydroDensityCommThisStep;

    double HydroAcc = (TimingResults.HydroAccThisStep
                      -TimingResults.HydroAccNeighborSearchThisStep
                      -TimingResults.HydroAccCommThisStep);
    TimingResults.HydroAcc += HydroAcc;
    TimingResults.HydroAccNeighborSearch += TimingResults.HydroAccNeighborSearchThisStep;
    TimingResults.HydroAccComm += TimingResults.HydroAccCommThisStep;


    TimingResults.HydroThisStep = HydroKernel+HydroDensity+HydroAcc;
    TimingResults.HydroNeighborSearchThisStep = (TimingResults.HydroKernelNeighborSearchThisStep
                                                +TimingResults.HydroDensityNeighborSearchThisStep
                                                +TimingResults.HydroAccNeighborSearchThisStep);
    TimingResults.HydroCommThisStep = (TimingResults.HydroKernelCommThisStep
                                       +TimingResults.HydroDensityCommThisStep
                                       +TimingResults.HydroAccCommThisStep);
    TimingResults.Hydro += TimingResults.HydroThisStep;
    TimingResults.HydroNeighborSearch += TimingResults.HydroNeighborSearchThisStep;
    TimingResults.HydroComm += TimingResults.HydroCommThisStep;
    TimingResults.HydroTree += TimingResults.HydroTreeThisStep;
    
    TimingResults.Cooling += TimingResults.CoolingThisStep;

    // Star formation
    TimingResults.Starformation += TimingResults.StarformationThisStep;

    // feedback
    TimingResults.Feedback += TimingResults.FeedbackThisStep
                             -TimingResults.FeedbackNeighborSearchThisStep
                             -TimingResults.FeedbackCommThisStep;
    TimingResults.FeedbackNeighborSearch += TimingResults.FeedbackNeighborSearchThisStep;
    TimingResults.FeedbackComm += TimingResults.FeedbackCommThisStep;
    
    // Sink particles
    TimingResults.Sink += TimingResults.SinkThisStep;
    TimingResults.SinkFormation += TimingResults.SinkFormationThisStep;
    TimingResults.SinkAccretion += TimingResults.SinkAccretionThisStep;
    TimingResults.SinkComm += TimingResults.SinkCommThisStep;

    // others
    TimingResults.Integral += TimingResults.IntegralThisStep;
    TimingResults.SortStructures += TimingResults.SortStructuresThisStep;
    TimingResults.KeyGeneration += TimingResults.KeyGenerationThisStep;
    TimingResults.TimeStep += TimingResults.TimeStepThisStep;

    // imbalance
    TimingResults.GravImbarance += TimingResults.GravImbaranceThisStep;
    TimingResults.HydroImbarance += TimingResults.HydroImbaranceThisStep;


    // particles
    //Pall.TActiveParticlesAll += Pall.NActivesAll_t; 
    Pall.TActiveParticles += Pall.NActives_t; 
    Pall.TActiveParticlesDM += Pall.NActivesDM_t; 
    Pall.TActiveParticlesHydro += Pall.NActivesHydro_t; 
    Pall.TActiveParticlesStar += Pall.NActivesStars_t; 
    Pall.TActiveParticlesSink += Pall.NActivesSink_t; 

    return;
}


void LogTotalMass(void){

    double Mass[NTypes],GlobalMass[NTypes];

    for(int i=0;i<NTypes;i++)
        Mass[i] = 0.e0;

    //
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Mass < 0.e0){
            fprintf(stderr,"[%02d] leaf %d, Type = %d Mass = %g\n",
                    MPIGetMyID(),i,PstarBody(i)->Type,PstarBody(i)->Mass);
            StructureReportPstar(i);
            fflush(NULL);
            assert(PstarBody(i)->Mass > 0.e0);
        }
    }
    
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
            Mass[TypeHydro] += Pbody[i]->Mass;
        } else if(Pbody[i]->Type == TypeStar){
            Mass[TypeStar] += Pbody[i]->Mass;
        } else if(Pbody[i]->Type == TypeSink){
            Mass[TypeSink] += Pbody[i]->Mass;
        } else if(Pbody[i]->Type == TypeDM){
            Mass[TypeDM] += Pbody[i]->Mass;
        }
        if(Pbody[i]->Mass < 0.e0){
            fprintf(stderr,"[%02d] leaf %d, Type = %d Mass = %g\n",
                    MPIGetMyID(),i,Pbody[i]->Type,Pbody[i]->Mass);
            StructureReportPbody(i);
            fflush(NULL);
            assert(Pbody[i]->Mass > 0.e0);
        }
    }

    MPI_Allreduce(Mass,GlobalMass,NTypes,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        double GlobalMassTotal = 0.e0;
        for(int i=0;i<NTypes;i++)
            GlobalMassTotal += GlobalMass[i];

        FILE *fp;
        FileOpen(fp,"./data/Mass.dat","a");
        fprintf(fp,"%g %1.15g ",Pall.TCurrent,GlobalMassTotal);
        for(int i=0;i<NTypes;i++)
            fprintf(fp,"%1.15g ",GlobalMass[i]);
        fprintf(fp,"\n");
        fclose(fp);
    }

    return ;
}
