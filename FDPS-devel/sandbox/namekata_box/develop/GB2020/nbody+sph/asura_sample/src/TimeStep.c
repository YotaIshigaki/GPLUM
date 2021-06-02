#include "config.h"
#include "ThermalConductivity.h"
#ifdef USE_RADIATION_PRESSURE //{
#include "RadiationPressure.h"
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
#include "StellarWind.h"
#endif // USE_STELLAR_WIND //}

//#define USE_TIME_DERIVATIVE_OF_U

#ifdef USE_ACC_DOT_TIMESTEP //{
#define TimeStepDiffAcc (ON)
#else // USE_ACC_DOT_TIMESTEP //} //{
#define TimeStepDiffAcc (OFF)
#endif // USE_ACC_DOT_TIMESTEP //} 


static inline double __attribute__((always_inline)) CFLTimeStepi(const int index){
    /* Springel 2005, Eq.16                */
    /* This version uses a signal velocity */
    /* in order to estimate the hydro time */
    /* step.                               */
    //return 2.e0*PbodyHydroKernelPred(index)/PbodyHydroVsig(index);
// #if VISCOSITY_TYPE == 0 //{
#ifdef USE_TIMESTEP_PREDICTOR_CHANGE
    return fmin(2.e0*PbodyHydroKernel(index)/PbodyHydroVsig(index),
                1.3863/fabs(PbodyHydro(index)->Div)); // 2*log(2) = 1.3863; log(2.0)=0.69315
#else 
    return 2.e0*PbodyHydroKernel(index)/PbodyHydroVsig(index);
#endif
// #elif VISCOSITY_TYPE == 1 //}//{
// 
    // double h = PbodyHydro(index)->KernelPred;
    // double cs = sqrt(Pall.GGm1*PbodyHydro(index)->UPred);
    // double DivV = PbodyHydro(index)->Div;
    // double Alpha = PbodyHydro(index)->Alpha;
// 
    // return h/(h*fabs(DivV)+cs+1.2*(Alpha*cs+2.0*Alpha*h*fabs(fmin(DivV,0.0))));
// 
// #endif // VISCOSITY_TYPE == 1 //}//{
}

#if 1
#ifdef USE_INDIVIDUAL_TIMESTEP
static inline void __attribute__((always_inline)) TimeStepi(const int index){

#ifdef GRAVITY_RUN
    double ParticleEpsSize = Pbody[index]->Eps*Pall.AdaptiveSofteningFactor;
    double dt = TFactorAcc*sqrt(ParticleEpsSize/NORM(Pbody[index]->Acc));

#if TimeStepDiffAcc
    if(Pbody[index]->Type == TypeHydro){
        double Acc[3] = {Pbody[index]->Acc[0],Pbody[index]->Acc[1],Pbody[index]->Acc[2]};
        double normAcc = NORM(Acc);
        double DiffAcc[3] = {Acc[0]-Pbody[index]->AccOld[0],
                             Acc[1]-Pbody[index]->AccOld[1],
                             Acc[2]-Pbody[index]->AccOld[2]};
        double normDiffAcc = NORM(DiffAcc);
        double normVel = NORM(Pbody[index]->Vel);

        if((normDiffAcc > TINY)&&(Pbody[index]->dt>0.e0))
            dt = fmin(dt,TFactorDiffAcc*normAcc/normDiffAcc*Pbody[index]->dt);
    }
#endif

#ifdef USE_RADIATION_PRESSURE //{
    if(Pbody[index]->Type == TypeStar){
        //if(Pall.TCurrent-PbodyStar(index)->FormationTime < GetRadiationPressureSourceLifeTime()){
        if(CheckRadiationPressureSourceThroughPbody(index)){
            dt = fmin(dt,GetRadiationPressureMaxTimestep());
        }
    }
#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_FEEDBACK_TIMESTEP //{
    if(Pbody[index]->Type == TypeStar){
        dt = PbodyStar(index)->dt_fb>0.e0
            ?fmin(dt,PbodyStar(index)->dt_fb)
            :dt;
        PbodyStar(index)->dt_fb = 0.e0;
    }
#endif // USE_FEEDBACK_TIMESTEP //}


#ifdef USE_STELLAR_WIND //{
#ifdef USE_STELLAR_WIND_TIME_STEP //{
    if(Pbody[index]->Type == TypeStar){
        if(CheckStellarWindSourceThroughPbody(index)){
            dt = fmin(dt,GetStellarWindMaxTimestep());
            //dt = fmin(dt,1.e+5*YEAR_CGS/Pall.UnitTime);
        }
    }
#endif // USE_STELLAR_WIND_TIME_STEP //}
#endif // USE_RADIATION_PRESSURE //}

#ifdef COSMOLOGICAL_RUN
    //if(Pall.hubble > TINY)
        //dt = fmin(dt,0.01/Pall.HubbleZ);
    dt = fmin(dt,0.01/Pall.HubbleZ);
    Pbody[index]->dt = dt;
#endif 
    Pbody[index]->dt = dt;
#endif

#ifdef USE_SINK_TIMESTEP_LIMITER
    static const double SinkTimeStepLimiterFactor = (double)(1<<SINK_TIMESTEP_MAX_K_LOCAL); 
    double dt_sink;
    if(Pbody[index]->Type == TypeSink){
        if(PbodySink(index)->dt_localmin > TINY){
#ifdef GRAVITY_RUN
            dt_sink = fmin(dt,SinkTimeStepLimiterFactor*PbodySink(index)->dt_localmin);
#else
            dt_sink = SinkTimeStepLimiterFactor*PbodySink(index)->dt_localmin;
#endif
        }
        Pbody[index]->dt = dt_sink;
        if(Pbody[index]->dt == 0.e0)
            Pbody[index]->dt = Pall.OutPutInterval;
    }
#endif

    if(Pbody[index]->Type == TypeHydro){
        // double dt_hydro = TFactorCourant*fmin(CFLTimeStepi(index),sqrt(2.0*PbodyHydro(index)->Kernel/NORM(PbodyHydroAcc(index))));
        double dt_hydro = TFactorCourant*CFLTimeStepi(index);

#ifdef USE_TIME_DERIVATIVE_OF_U //{
        dt_hydro = fmin(dt_hydro,TFactorCourant*fabs(PbodyHydro(index)->U/PbodyHydro(index)->Du));
#endif //USE_TIME_DERIVATIVE_OF_U //}

#ifdef USE_THERMAL_CONDUCTIVITY // {
#ifdef USE_THERMAL_CONDUCTIVITY_TIMESTEP //{
        double Temperaturei = Pall.ConvertUtoT*PbodyHydro(index)->UPred;
        double Kappai = CalcThermalConductivityByTemperature(Temperaturei);
        dt_hydro = fmin(dt_hydro,THERMAL_CONDUCTIVITY_TIMESTEP_FACTOR*SQ(2.0*PbodyHydro(index)->Kernel)/(2.0*Kappai));
        // if(Pbody[index]->PosP[0] < 0.53)
            // dt_hydro /= 3.0;
        dt_hydro = dt_hydro*(gsl_rng_uniform(RandomGenerator)+0.5);
#endif //USE_THERMAL_CONDUCTIVITY_TIMESTEP //}
#endif //USE_THERMAL_CONDUCTIVITY //}

#ifdef GRAVITY_RUN
        PbodyHydro(index)->dt_hydro = fmin(dt_hydro,Pbody[index]->dt);
#else // GRAVITY_RUN
        Pbody[index]->dt = PbodyHydro(index)->dt_hydro = dt_hydro;
#endif // GRAVITY_RUN

#ifdef GRAVITY_RUN //{ 
#if MAX_K_GRAVITY >= 0 //{
        const static int fact = 1 << MAX_K_GRAVITY;
        Pbody[index]->dt = fmin(Pbody[index]->dt,fact*PbodyHydro(index)->dt_hydro);
#endif // MAX_K_GRAVITY //}

#ifndef USE_FAST_SCHEME //{
        Pbody[index]->dt = PbodyHydro(index)->dt_hydro;
#endif // ifndef USE_FAST_SCHEME //}
#endif // GRAVITY_RUN //}

    }

    return;
}
#else //Do not use USE_VARIABLE_TIMESTEP
inline void __attribute__((always_inline)) TimeStepi(const int index){

    if(Pbody[index]->Type == TypeHydro)
        PbodyHydro(index)->dt_hydro = Pall.dt_const;
    Pbody[index]->dt = Pall.dt_const;
    return;
}
#endif


#else
inline void __attribute__((always_inline)) TimeStepi(const int index){

    double ParticleEpsSize = Pbody[index]->Eps*Pall.AdaptiveSofteningFactor;

    double Vel[3] = {Pbody[index]->Vel[0],Pbody[index]->Vel[1],Pbody[index]->Vel[2]};

    double Acc[3] = {Pbody[index]->Acc[0]+PbodyHydroAcc(index)[0],
                     Pbody[index]->Acc[1]+PbodyHydroAcc(index)[1],
                     Pbody[index]->Acc[2]+PbodyHydroAcc(index)[2]};


    if(Pbody[index]->Type == TypeHydro){
        double dt_hydro = TFactorCourant*CFLTimeStepi(index);
        //double cs = sqrt(Pall.GGm1*PbodyHydroU(index));
        //double dt = 0.3*fmin(sqrt(ParticleEpsSize/NORM(Acc)),PbodyHydroKernel(index)/(NORM(Vel)+cs));
        double dt = TFactorCourant*sqrt(2.0*PbodyHydroKernel(index)->Kernel/NORM(PbodyHydroAcc(index)));
        Pbody[index]->dt = PbodyHydro(index)->dt_hydro = fmin(dt_hydro,dt);
        //Pbody[index]->dt = PbodyHydro(index)->dt_hydro = dt_hydro;
    } else {
        double dt = 0.3*fmin(sqrt(ParticleEpsSize/NORM(Acc)),PbodyHydroKernel(index)/NORM(Vel));
        Pbody[index]->dt = dt;
    }
}

#endif

static inline double __attribute__((always_inline)) CFLTimeStepHydroi(const int index){
    /* Springel 2005, Eq.16                */
    /* This version uses a signal velocity */
    /* in order to estimate the hydro time */
    /* step.                               */

// #if VISCOSITY_TYPE == 0 //}//{
#ifdef USE_TIMESTEP_PREDICTOR_CHANGE
    return fmin(2.e0*Phydro[index]->Kernel/Phydro[index]->Vsig,
                1.3863/fabs(Phydro[index]->Div)); // 2*log(2) = 1.3863; log(2.0)=0.69315
#else 
    /*
    fprintf(stderr,"%d %ld %g %g %g [km/s]\n",index,PhydroBody(index)->GlobalID,
            2.e0*Phydro[index]->Kernel/Phydro[index]->Vsig*Pall.UnitTime/YEAR_CGS,
             Phydro[index]->Kernel*Pall.UnitLength/PC_CGS,
             Phydro[index]->Vsig*(Pall.UnitLength)/(Pall.UnitTime)/1.e5);
             */
    return 2.e0*Phydro[index]->Kernel/Phydro[index]->Vsig;
#endif
// #elif VISCOSITY_TYPE == 1 //}//{
    // double h = Phydro[index]->KernelPred;
    // double cs = sqrt(Pall.GGm1*Phydro[index]->UPred);
    // double DivV = Phydro[index]->Div;
    // double Alpha = Phydro[index]->Alpha;
// 
    // return h/(h*fabs(DivV)+cs+1.2*(Alpha*cs+2.0*Alpha*h*fabs(fmin(DivV,0.0))));
// #endif // VISCOSITY_TYPE == 1 //}//{
}

static inline void __attribute__((always_inline)) TimeStepHydroi(const int index){

    //double ParticleEpsSize = PhydroBody(index)->Eps*Pall.AdaptiveSofteningFactor;

    double dt_hydro = TFactorCourant*CFLTimeStepHydroi(index);

#ifdef USE_THERMAL_CONDUCTIVITY //{
#ifdef USE_THERMAL_CONDUCTIVITY_TIMESTEP //{
    double Kappai = CalcThermalConductivity(index);
    dt_hydro = fmin(dt_hydro,THERMAL_CONDUCTIVITY_TIMESTEP_FACTOR*SQ(2.0*Phydro[index]->Kernel)/(2.0*Kappai));
#endif //USE_THERMAL_CONDUCTIVITY_TIMESTEP //}
#endif // USE_THERMAL_CONDUCTIVITY //}

#ifdef USE_FAST_SCHEME 
    Phydro[index]->dt_hydro = dt_hydro;
#else
    PhydroBody(index)->dt = Phydro[index]->dt_hydro = dt_hydro;
#endif // USE_FAST_SCHEME

    return ;
}

static inline void __attribute__((always_inline)) TimeStepGravi(const int index){
#ifdef USE_VARIABLE_TIMESTEP
    double ParticleEpsSize = Pbody[index]->Eps*Pall.AdaptiveSofteningFactor;
    double dt = TFactorAcc*sqrt(ParticleEpsSize/NORM(Pbody[index]->Acc));
    //dt = fmin(dt,TFactorVel*ParticleEpsSize/NORM(Pbody[index]->Vel));

#if TimeStepDiffAcc
    if(Pbody[index]->Type == TypeHydro){
        double Acc[3] = {Pbody[index]->Acc[0],Pbody[index]->Acc[1],Pbody[index]->Acc[2]};
        double normAcc = NORM(Acc);
        double DiffAcc[3] = {Acc[0]-Pbody[index]->AccOld[0],
                             Acc[1]-Pbody[index]->AccOld[1],
                             Acc[2]-Pbody[index]->AccOld[2]};
        double normDiffAcc = NORM(DiffAcc);
        double normVel = NORM(Pbody[index]->Vel);

        if(normDiffAcc > TINY){
            //dt = TFactor*fmax(sqrt(normVel/normDiffAcc),normAcc/normDiffAcc);
            //dt = fmin(dt,TFactorDiffAcc*normAcc/normDiffAcc);
            // a/a_dot, where a_dot = DiffAcc/dt
            /*
            fprintf(stderr,"%ld dt %g dt_n %g normA %g normDA %g dt_old %g | %d\n",Pbody[index]->GlobalID,
                    dt,TFactorDiffAcc*normAcc/normDiffAcc*Pbody[index]->dt,
                    normAcc,normDiffAcc,Pbody[index]->dt,Pbody[index]->Type);
            fprintf(stderr,"%g %g %g | %g %g %g | %g %g %g\n",
                    Acc[0], Acc[1], Acc[2],
                    Pbody[index]->AccOld[0],
                    Pbody[index]->AccOld[1],
                    Pbody[index]->AccOld[2],
                    Acc[0]-Pbody[index]->AccOld[0],
                    Acc[1]-Pbody[index]->AccOld[1],
                    Acc[2]-Pbody[index]->AccOld[2]);
            fflush(NULL);
            */
            dt = fmin(dt,TFactorDiffAcc*normAcc/normDiffAcc*Pbody[index]->dt);
        }
    }
#endif

#ifdef USE_RADIATION_PRESSURE //{
    if(Pbody[index]->Type == TypeStar){
        //if(Pall.TCurrent-PbodyStar(index)->FormationTime < GetRadiationPressureSourceLifeTime()){
        if(CheckRadiationPressureSourceThroughPbody(index)){
            dt = fmin(dt,GetRadiationPressureMaxTimestep());
        }
    }
#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_FEEDBACK_TIMESTEP //{
    if(Pbody[index]->Type == TypeStar){
        dt = PbodyStar(index)->dt_fb>0.e0
            ?fmin(dt,PbodyStar(index)->dt_fb)
            :dt;
        PbodyStar(index)->dt_fb = 0.e0;
    }
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_STELLAR_WIND //{
#ifdef USE_STELLAR_WIND_TIME_STEP //{
    if(Pbody[index]->Type == TypeStar){
        if(CheckStellarWindSourceThroughPbody(index)){
            dt = fmin(dt,GetStellarWindMaxTimestep());
            //dt = fmin(dt,1.e+5*YEAR_CGS/Pall.UnitTime);
        }
    }
#endif // USE_STELLAR_WIND_TIME_STEP //}
#endif // USE_STELLAR_WIND //}

#ifdef COSMOLOGICAL_RUN //{
    if(Pall.hubble > TINY)
        dt = fmin(dt,0.01/Pall.HubbleZ);
#endif // COSMOLOGICAL_RUN //}

    Pbody[index]->dt = dt;
#else
    Pbody[index]->dt = Pall.dt_const;
#endif

    return;
}


void FirstTimeStep(void){

    double TimeStepThisStep = GetElapsedTime();

    for(int i=0;i<Pall.Ntotal;i++){
        TimeStepi(i);
        if(Pbody[i]->Type == TypeHydro)
            Pbody[i]->dt = fmin(Pbody[i]->dt,PbodyHydro(i)->dt_hydro);
    }

    TimingResults.TimeStepThisStep += GetElapsedTime()-TimeStepThisStep;
    return ;
}

void BuildHierarchicalTimeStep(void){

    double TimeStepThisStep = GetElapsedTime();

    double MinTime,MaxTime;
    double MinMax[2],GlobalMinMax[2];

    MinTime = MaxTime = Pbody[0]->dt;
    for(int i=1;i<Pall.Ntotal;i++){
        MinTime = fmin(MinTime,Pbody[i]->dt);
        MaxTime = fmax(MaxTime,Pbody[i]->dt);
    }
    for(int i=0;i<Pall.Nhydro;i++){
        MinTime = fmin(MinTime,Phydro[i]->dt_hydro);
        MaxTime = fmax(MaxTime,Phydro[i]->dt_hydro);
    }
    MinMax[0] = -MinTime;
    MinMax[1] = MaxTime;

    MPI_Allreduce(MinMax,GlobalMinMax,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    Pall.dtmin = -GlobalMinMax[0];
    Pall.dtmax = GlobalMinMax[1];

    Pall.kmax = (int)(log2(Pall.dtmax/Pall.dtmin));

#ifdef COSMOLOGICAL_RUN //{
    double Toffset = CalcZtoT(Pall.InitialRedshift);
    if(Pall.TCurrent-Toffset+Pall.dtmax > Pall.OutPutFileNumber*Pall.OutPutInterval){
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"TCurrent %g, Pall.dtmax %g, OutPutFileNumber %d, OutPutInterval %g\n",
                    Pall.TCurrent,Pall.dtmax,Pall.OutPutFileNumber,Pall.OutPutInterval);
        Pall.dtmax = Pall.OutPutFileNumber*Pall.OutPutInterval-Pall.TCurrent+Toffset;
        if(Pall.dtmax < Pall.dtmin){
            Pall.dtmin = Pall.dtmax;
            Pall.kmax = 0;
        }else{
            Pall.kmax = (int)(log2(Pall.dtmax/Pall.dtmin))+1;
            Pall.dtmin = Pall.dtmax*exp2(-(double)Pall.kmax);
        }
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"TCurrent %g, Pall.dtmin %g, Pall.kmax %d, Toffset = %g\n",
                    Pall.TCurrent,Pall.dtmin,Pall.kmax,Toffset);
    }
#else // COSMOLOGICAL_RUN //}//{
    if(Pall.TCurrent+Pall.dtmax > Pall.OutPutFileNumber*Pall.OutPutInterval){
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"TCurrent %g, Pall.dtmax %g, OutPutFileNumber %d, OutPutInterval %g\n",
                    Pall.TCurrent,Pall.dtmax,Pall.OutPutFileNumber,Pall.OutPutInterval);
        Pall.dtmax = Pall.OutPutFileNumber*Pall.OutPutInterval-Pall.TCurrent;
        if(Pall.dtmax < Pall.dtmin){
            Pall.dtmin = Pall.dtmax;
            Pall.kmax = 0;
        }else{
            Pall.kmax = (int)(log2(Pall.dtmax/Pall.dtmin))+1;
            Pall.dtmin = Pall.dtmax*exp2(-(double)Pall.kmax);
        }
    }
#endif // COSMOLOGICAL_RUN //}

#ifdef OUTPUT_CONSTANT_INTERVAL //{
#   ifdef COSMOLOGICAL_RUN //{
    if(Pall.TCurrent-Toffset+Pall.dtmin > Pall.OutPutFileNumber*Pall.OutPutInterval){
        Pall.dtmin = Pall.OutPutFileNumber*Pall.OutPutInterval-Pall.TCurrent+Toffset;
        Pall.dtmax = Pall.dtmin;
        Pall.kmax = 0;
    }
#   else // COSMOLOGICAL_RUN //}//{
    if(Pall.TCurrent+Pall.dtmin > Pall.OutPutFileNumber*Pall.OutPutInterval){
        Pall.dtmin = Pall.OutPutFileNumber*Pall.OutPutInterval-Pall.TCurrent;
        Pall.dtmax = Pall.dtmin;
        Pall.kmax = 0;
    }
#   endif // COSMOLOGICAL_RUN //}
#endif // OUTPUT_CONSTANT_INTERVAL //}

#ifdef DTMAX_RESCALE_FACTOR //{
    if((Pall.dtmax > 0.e0) && (Pall.kmax > DTMAX_RESCALE_LIMIT_KMAX)){
        Pall.dtmax *= DTMAX_RESCALE_FACTOR;
        if(Pall.dtmax < Pall.dtmin){
            Pall.dtmax = Pall.dtmin;
            Pall.kmax = 0;
        } else {
            Pall.kmax = (int)(log2(Pall.dtmax/Pall.dtmin));
            Pall.dtmin = Pall.dtmax*exp2(-(double)Pall.kmax);
        }
    }
#endif // DTMAX_RESCALE_FACTOR //}

    if(Pall.TCurrent+Pall.dtmax > Pall.TEnd){
        Pall.dtmax = Pall.TEnd-Pall.TCurrent;
        if(Pall.dtmax < Pall.dtmin){
            Pall.dtmin = Pall.dtmax;
            Pall.kmax = 0;
        }else{
            Pall.kmax = (int)(log2(Pall.dtmax/Pall.dtmin))+1;
            Pall.dtmin = Pall.dtmax*exp2(-(double)Pall.kmax);
        }
    }

    if(Pall.TCurrent+Pall.dtmin > Pall.TEnd){
        Pall.dtmin = Pall.TEnd-Pall.TCurrent;
        Pall.dtmax = Pall.dtmin;
        Pall.kmax = 0;
    }

    if(Pall.kmax > MaximumTimeHierarchy)
        Pall.kmax = MaximumTimeHierarchy;

    Pall.dtmax = Pall.dtmin*exp2((double)Pall.kmax);

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->dt < Pall.dtmax){
            Pbody[i]->k = (int)(log2(Pbody[i]->dt/Pall.dtmin));
            Pbody[i]->dt = Pall.dtmin*pow(2.e0,(double)Pbody[i]->k);
        } else {
            Pbody[i]->k = Pall.kmax;
            Pbody[i]->dt = Pall.dtmax;
        }
        Pbody[i]->EraLocal = 0.e0;
        if(Pbody[i]->Type == TypeHydro){
            if(PbodyHydro(i)->dt_hydro < Pall.dtmax){
                PbodyHydro(i)->k_hydro = (int)(log2(PbodyHydro(i)->dt_hydro/Pall.dtmin));
#ifdef HYDRO_TIMESTEP_LIMITER
                //if(PbodyHydro(i)->k_hydro-PbodyHydro(i)->k_hydro_localmin > MAX_K_LOCAL){
                    //PbodyHydro(i)->k_hydro = PbodyHydro(i)->k_hydro_localmin+MAX_K_LOCAL;
                if(PbodyHydro(i)->k_hydro > PbodyHydro(i)->k_hydro_localmin){
                    PbodyHydro(i)->k_hydro = PbodyHydro(i)->k_hydro_localmin;
                }
#endif
                PbodyHydro(i)->dt_hydro = Pall.dtmin*pow(2.e0,(double)PbodyHydro(i)->k_hydro);
#ifndef GRAVITY_RUN
                Pbody[i]->k = PbodyHydro(i)->k_hydro;
                Pbody[i]->dt = PbodyHydro(i)->dt_hydro;
#endif
            } else {
                PbodyHydro(i)->k_hydro = Pall.kmax;
                PbodyHydro(i)->dt_hydro = Pall.dtmax;
            }
            PbodyHydro(i)->EraLocal_hydro = 0.e0;
        }
    }

    Pall.Era = Pall.TCurrent+Pall.dtmax;
    Pall.EraStart = Pall.TCurrent; 
    Pall.EraLocal = 0.e0;
    Pall.EraLocalEnd = Pall.dtmax;

    TimingResults.TimeStepThisStep += GetElapsedTime()-TimeStepThisStep;

    return;
}


#ifdef USE_FAST_SCHEME // define USE_FAST_SCHEME
void BuildNewTimeStep(void){

    double TimeStepThisStep = GetElapsedTime();
#ifdef USE_SINK_TIMESTEP_LIMITER
    static const double SinkTimeStepLimiterFactor = (double)(1<<SINK_TIMESTEP_MAX_K_LOCAL); 
#endif

    int step;
    double dtold,dtnew;

    for(int i=0;i<Pall.Ntotal;i++){ // Gravity part.
        if(Pbody[i]->Active){

            dtold = Pbody[i]->dt;
            step = (int)((Pbody[i]->EraLocal/dtold)+0.5);

            Pbody[i]->EraLocal = step*dtold;
            TimeStepGravi(i);

#ifdef USE_SINK_TIMESTEP_LIMITER
            if(Pbody[i]->Type == TypeSink){
                if(PbodySink(i)->dt_localmin > TINY){
                    Pbody[i]->dt = fmin(Pbody[i]->dt,SinkTimeStepLimiterFactor*PbodySink(i)->dt_localmin);
                }
            }
#endif // USE_SINK_TIMESTEP_LIMITER

            if(Pbody[i]->dt < dtold){
                do{
                    Pbody[i]->k--;
                    dtnew = Pall.dtmin*exp2((double)(Pbody[i]->k));
                }while(Pbody[i]->dt < dtnew);
                Pbody[i]->dt = dtnew;
            } else if((step%2 == 0)&&(Pbody[i]->dt > 2.e0*dtold)){
                Pbody[i]->k ++;
                Pbody[i]->dt = Pall.dtmin*exp2((double)(Pbody[i]->k));
            } else {
                Pbody[i]->dt = dtold;
            }
        }
    }


    double dt_localmin;
    for(int i=0;i<Pall.Nhydro;i++){ // Hydro part.
        if(Phydro[i]->Active){
            double EL = Phydro[i]->EraLocal_hydro;
#ifdef HYDRO_TIMESTEP_LIMITER
            if(Phydro[i]->HydroTimeStepLimiterFlag){
                int oldk = Phydro[i]->k_hydro;
                int olddt = Phydro[i]->dt_hydro;
                //Phydro[i]->k_hydro =  Phydro[i]->k_hydro_localmin_old+MAX_K_LOCAL;
                Phydro[i]->k_hydro =  Phydro[i]->k_hydro_localmin_old;

                assert(Phydro[i]->k_hydro_localmin_old > -100);

                dt_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro));
                step = (int)((Pall.EraLocal/dt_localmin)+0.5);
                Phydro[i]->EraLocal_hydro = step*dt_localmin;
                dtold = dt_localmin;
            } else 
#endif // HYDRO_TIMESTEP_LIMITER
            {
                dtold = Phydro[i]->dt_hydro;
                step = (int)((Phydro[i]->EraLocal_hydro/dtold)+0.5);
                Phydro[i]->EraLocal_hydro = step*dtold;
            }

            TimeStepHydroi(i);


#ifdef HYDRO_TIMESTEP_LIMITER
            // Phydro[i]->dt_hydro = fmin(Phydro[i]->dt_hydro,Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin+MAX_K_LOCAL)));
            // Phydro[i]->k_hydro = MIN(Phydro[i]->k_hydro,Phydro[i]->k_hydro_localmin+MAX_K_LOCAL);
            Phydro[i]->dt_hydro = fmin(Phydro[i]->dt_hydro,Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin)));
            Phydro[i]->k_hydro = MIN(Phydro[i]->k_hydro,Phydro[i]->k_hydro_localmin);
#endif // HYDRO_TIMESTEP_LIMITER
            Phydro[i]->dt_hydro = fmin(PhydroBody(i)->dt,Phydro[i]->dt_hydro);


            if(Phydro[i]->dt_hydro < dtold ){ // Shrink time step
                do{
                    Phydro[i]->k_hydro--;
                    dtnew = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro));
                }while(Phydro[i]->dt_hydro < dtnew);
                Phydro[i]->dt_hydro = dtnew;
                // fprintf(stderr,"shrink dt: %ld %d %g < %g\n",
                        // PhydroBody(i)->GlobalID,
                        // step,Phydro[i]->dt_hydro,dtold);
            //} else if((step % 2 == 0)&&(Phydro[i]->dt_hydro > 2.e0*dtold)&&(Phydro[i]->Du*2.e0*dtold)){ // Double time step
            } else if((step % 2 == 0)&&(Phydro[i]->dt_hydro > 2.e0*dtold)
                    // &&(Phydro[i]->U+Phydro[i]->Du*2.e0*dtold>0.e0)
                    ){ // Double time step
#ifdef HYDRO_TIMESTEP_LIMITER  //{
                if(Phydro[i]->HydroTimeStepLimiterFlag){
                    // fprintf(stderr,"test: %ld %d %g > %g | %g > %g\n",
                            // PhydroBody(i)->GlobalID,
                            // step,Phydro[i]->dt_hydro,2*dtold,
                            // dt_localmin,2*dtold);
#if 1
                    if(dt_localmin > 2*dtold){
                        //Phydro[i]->k_hydro = (Phydro[i]->k_hydro_localmin_old+MAX_K_LOCAL+1);
                        Phydro[i]->k_hydro = (Phydro[i]->k_hydro_localmin_old+1);
                        Phydro[i]->dt_hydro = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro));
                    } else {
                        Phydro[i]->k_hydro = (int)log2(dtold/Pall.dtmin);
                        Phydro[i]->dt_hydro = dtold;
                    }
#else
                    Phydro[i]->k_hydro = (int)log2(dtold/Pall.dtmin);
                    Phydro[i]->dt_hydro = dtold;
#endif
                } else 
#endif //HYDRO_TIMESTEP_LIMITER //}
                {
                    // fprintf(stderr,"test2: %ld %d %g > %g\n",
                            // PhydroBody(i)->GlobalID,
                            // step,Phydro[i]->dt_hydro,2*dtold);
                    Phydro[i]->k_hydro ++;
                    Phydro[i]->dt_hydro = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro));
                }
            } else { // Keep time step
                Phydro[i]->dt_hydro = dtold;
            }
#ifdef HYDRO_TIMESTEP_LIMITER //{
            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;
#endif // HYDRO_TIMESTEP_LIMITER //}

#if MAX_K_GRAVITY >= 0 //{
            const static int fact = 1 << MAX_K_GRAVITY;
            if(PhydroBody(i)->Active){
                if( PhydroBody(i)->dt > fact*Phydro[i]->dt_hydro ){
                    PhydroBody(i)->k = Phydro[i]->k_hydro+MAX_K_GRAVITY;
                    PhydroBody(i)->dt = Pall.dtmin*exp2((double)(PhydroBody(i)->k));
                }
            }
#endif // MAX_K_GRAVITY //}
        }
    }

    TimingResults.TimeStepThisStep += GetElapsedTime()-TimeStepThisStep;

    return;
}

#else // ndef USE_FAST_SCHEME
void BuildNewTimeStep(void){

    double TimeStepThisStep = GetElapsedTime();

    int step;
    double dtold,dtnew;

    double dt_localmin;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
#ifdef HYDRO_TIMESTEP_LIMITER
            if((Pbody[i]->Type == TypeHydro)&&(PbodyHydro(i)->HydroTimeStepLimiterFlag)){
                Pbody[i]->k = PbodyHydro(i)->k_hydro = 
                    //PbodyHydro(i)->k_hydro_localmin_old+MAX_K_LOCAL;
                    PbodyHydro(i)->k_hydro_localmin_old;
                dt_localmin = Pall.dtmin*exp2((double)(PbodyHydro(i)->k_hydro));
                step = (int)((Pall.EraLocal/dt_localmin)+0.5);
                Pbody[i]->EraLocal = PbodyHydro(i)->EraLocal_hydro = step*dt_localmin;
                dtold = dt_localmin;
            } else 
#endif
            {
                dtold = Pbody[i]->dt;
                step = (int)((Pbody[i]->EraLocal/dtold)+0.5);
                Pbody[i]->EraLocal = step*dtold;
            }

            TimeStepi(i);
            if(Pbody[i]->Type == TypeHydro){
#ifdef HYDRO_TIMESTEP_LIMITER
                //PbodyHydro(i)->dt_hydro = Pbody[i]->dt
                    //= fmin(Pbody[i]->dt,Pall.dtmin*exp2((double)(PbodyHydro(i)->k_hydro_localmin+MAX_K_LOCAL)));
                //Pbody[i]->k = MIN(Pbody[i]->k,PbodyHydro(i)->k_hydro_localmin+MAX_K_LOCAL);
                PbodyHydro(i)->dt_hydro = Pbody[i]->dt
                    = fmin(Pbody[i]->dt,Pall.dtmin*exp2((double)(PbodyHydro(i)->k_hydro_localmin)));
                Pbody[i]->k = MIN(Pbody[i]->k,PbodyHydro(i)->k_hydro_localmin);
#endif
#ifdef GRAVITY_RUN
                Pbody[i]->dt = fmin(Pbody[i]->dt,PbodyHydro(i)->dt_hydro);
#endif
            }

            if(Pbody[i]->dt < dtold ){
                do{
                    Pbody[i]->k--;
                    dtnew = Pall.dtmin*exp2((double)(Pbody[i]->k));
                }while(Pbody[i]->dt < dtnew);
                Pbody[i]->dt = dtnew;
            } else if((step % 2 == 0)&&(Pbody[i]->dt > 2.e0*dtold)){
#ifdef HYDRO_TIMESTEP_LIMITER
                if((Pbody[i]->Type == TypeHydro)&&(PbodyHydro(i)->HydroTimeStepLimiterFlag)){
                    if(dt_localmin > 2*dtold){
                    //if((dt_localmin > 2*dtold)&&(PbodyHydro(i)->U+2*dtold*PbodyHydro(i)->Du > 0.e0)){
                        Pbody[i]->k = PbodyHydro(i)->k_hydro 
                            //= (PbodyHydro(i)->k_hydro_localmin_old+MAX_K_LOCAL+1);
                            = (PbodyHydro(i)->k_hydro_localmin_old+1);
                        Pbody[i]->dt = PbodyHydro(i)->dt_hydro 
                            = Pall.dtmin*exp2((double)(PbodyHydro(i)->k_hydro));
                    }else{
                        Pbody[i]->k = PbodyHydro(i)->k_hydro = (int)log2(dtold/Pall.dtmin);
                        Pbody[i]->dt = dtold;
                    }
                } else 
#endif
                {
                    Pbody[i]->k ++;
                    Pbody[i]->dt = Pall.dtmin*exp2((double)(Pbody[i]->k));
                }
            } else {
                Pbody[i]->dt = dtold;
            }
#ifdef HYDRO_TIMESTEP_LIMITER
            if(Pbody[i]->Type == TypeHydro){
                PbodyHydro(i)->k_hydro_localmin = MaximumTimeHierarchy;
                PbodyHydro(i)->NextUpdateEra = Pall.TEnd;
            }
#endif

            if(Pbody[i]->Type == TypeHydro){
                PbodyHydro(i)->dt_hydro = Pbody[i]->dt;
                PbodyHydro(i)->k_hydro = Pbody[i]->k;
                PbodyHydro(i)->EraLocal_hydro = Pbody[i]->EraLocal;
            }
        }
    }
    TimingResults.TimeStepThisStep += GetElapsedTime()-TimeStepThisStep;

    return;
}

#endif

int Kcurrent;

void RaiseActiveFlags(void){

    double TimingResultThisRoutine = GetElapsedTime();

    double dtnow = Pbody[0]->dt;
    double EraLocalnow = Pbody[0]->EraLocal + Pbody[0]->dt;

    if(Pall.Ntotal == 0){
        dtnow = EraLocalnow = Pall.TEnd;
    }

    for(int i=1;i<Pall.Ntotal;i++){
        dtnow = fmin(dtnow,Pbody[i]->dt);
        EraLocalnow = fmin(EraLocalnow,Pbody[i]->EraLocal+Pbody[i]->dt);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        dtnow = fmin(dtnow,Phydro[i]->dt_hydro);
        EraLocalnow = fmin(EraLocalnow,Phydro[i]->EraLocal_hydro+Phydro[i]->dt_hydro);
    }

    double TimeMin[2],GlobalTimeMin[2];
    TimeMin[0] = dtnow;
    TimeMin[1] = EraLocalnow;

    MPI_Allreduce(TimeMin,GlobalTimeMin,2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    dtnow = GlobalTimeMin[0];
    EraLocalnow = GlobalTimeMin[1];

    Pall.TStep = (int)(EraLocalnow/dtnow + 0.5e0);
    int TStepEnd = (int)(Pall.dtmax/dtnow+0.5);

    double kdummy= log2(dtnow/Pall.dtmin);
    int ktilde = (dtnow >= Pall.dtmin) ? (int)(kdummy + 0.5e0) : (int)(kdummy - 0.5e0);

    Kcurrent = ktilde;

    Pall.dtnow = Pall.dtmin*exp2((double)ktilde);
    /*
    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"#TIME STEP# dtnow = %g, EraLocal = %g, NextEra = %g, TCurrent = %g\n",
                dtnow,EraLocalnow,Pall.EraLocal+dtnow,Pall.TCurrent);
        fprintf(stderr,"#TIME STEP# Pall.dtnow = %g, Pall.dtmin = %g, ktilde = %d\n",Pall.dtnow,Pall.dtmin,ktilde);
        fflush(NULL);
    }
    */

    double EraMin = Pall.EraLocal+0.9*(Pall.dtnow);
    double EraMax = Pall.EraLocal+1.1*(Pall.dtnow);
    Pall.NActivesAll = 0;
    int count_shortcut = 0;
    for(int i=0;i<Pall.Ntotal;i++){ // Count number of active particles.
        double Timei = Pbody[i]->EraLocal+Pbody[i]->dt;
        if( (EraMin<=Timei)&&(Timei<=EraMax) ){
            Pall.NActivesAll ++;
            continue;
        }
        if(Pbody[i]->Type == TypeHydro){
#ifdef HYDRO_TIMESTEP_LIMITER //{
            double Timei_hydro = PbodyHydro(i)->EraLocal_hydro+PbodyHydro(i)->dt_hydro;
            if( (EraMin<=Timei_hydro)&&(Timei_hydro<=EraMax) ){
                Pall.NActivesAll ++;
            } else if(PbodyHydro(i)->NextUpdateEra <= EraMax){
                count_shortcut ++;
                Pall.NActivesAll ++;
            }
#else // HYDRO_TIMESTEP_LIMITER //}//{
            double Timei_hydro = PbodyHydro(i)->EraLocal_hydro+PbodyHydro(i)->dt_hydro;
            if( (EraMin<=Timei_hydro)&&(Timei_hydro<=EraMax) ){
                Pall.NActivesAll ++;
            }
#endif // HYDRO_TIMESTEP_LIMITER //}
        }
    }

    // fprintf(stderr,"This step:  %g <-%g-> %g, dtnow = %g, %s:%d\n",
            // EraMin,Pall.EraLocal,EraMax,Pall.dtnow,__FUNCTION__,__LINE__);
    // fflush(NULL);

    Pall.NActivesHydro = 0;
    for(int i=0;i<Pall.Nhydro;i++){ // Count number of active hydro particles.
        Phydro[i]->HydroTimeStepLimiterFlag = false;
        double Timei = Phydro[i]->EraLocal_hydro+Phydro[i]->dt_hydro;

        if( (EraMin<=Timei) && (Timei<=EraMax) ){
            Phydro[i]->EraLocal_hydro += Phydro[i]->dt_hydro;
            Phydro[i]->Active = ON;
            Pall.NActivesHydro ++;
#ifdef HYDRO_TIMESTEP_LIMITER
            if (Phydro[i]->NextUpdateEra <= EraMax){
                Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
                Phydro[i]->NextUpdateEra = Pall.TEnd;
            }
#endif
        }
#ifdef HYDRO_TIMESTEP_LIMITER
        else  if (Phydro[i]->NextUpdateEra <= EraMax){
            if (Phydro[i]->NextUpdateEra <= EraMin){
                fprintf(stderr,"The problem [%ld] NextEra %g, EraLocal %g, dt %g, %g <-%g-> %g\n",
                    PhydroBody(i)->GlobalID,Phydro[i]->NextUpdateEra,Phydro[i]->EraLocal_hydro,
                        Phydro[i]->dt_hydro,EraMin,Pall.EraLocal,EraMax);
            }
            if(Phydro[i]->dt_hydro < Pall.EraLocal + Pall.dtnow - Phydro[i]->EraLocal_hydro){
                fprintf(stderr,"Some error [%d], %g and %g \n",
                        i,Phydro[i]->dt_hydro,Pall.EraLocal + Pall.dtnow - Phydro[i]->EraLocal_hydro);
                fprintf(stderr,"hydro dt %g and Era %g \n",Phydro[i]->dt_hydro,Phydro[i]->EraLocal_hydro);
                fprintf(stderr,"all   dt %g and Era %g \n",Pall.dtnow,Pall.EraLocal);
                MPI_Abort(MPI_COMM_WORLD,-1);
            }

            Phydro[i]->dt_hydro = (Pall.EraLocal + Pall.dtnow) - Phydro[i]->EraLocal_hydro;
            //Phydro[i]->EraLocal_hydro += Phydro[i]->dt_hydro;
            Phydro[i]->EraLocal_hydro = Phydro[i]->NextUpdateEra;
            //Phydro[i]->dt_hydro_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin+MAX_K_LOCAL));
            Phydro[i]->dt_hydro_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin));
            Phydro[i]->Active = ON;

            //Phydro[i]->k_hydro = Phydro[i]->k_hydro_localmin+MAX_K_LOCAL;
            Phydro[i]->k_hydro_localmin_old = Phydro[i]->k_hydro_localmin;
            Phydro[i]->HydroTimeStepLimiterFlag = true;
            Pall.NActivesHydro ++;

            Phydro[i]->k_hydro_localmin = MaximumTimeHierarchy;
            Phydro[i]->NextUpdateEra = Pall.TEnd;

#ifndef USE_FAST_SCHEME //{
            PhydroBody(i)->dt = Phydro[i]->dt_hydro;
            PhydroBody(i)->k = Phydro[i]->k_hydro_localmin;
            PhydroBody(i)->EraLocal = Phydro[i]->EraLocal_hydro-Phydro[i]->dt_hydro;
#endif // USE_FAST_SCHEME //}
        }
#endif
        else{
            Phydro[i]->Active = OFF;
        }
    }

    Pall.NActives = Pall.NActivesDM = Pall.NActivesStars = Pall.NActivesSink = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        double Timei = Pbody[i]->EraLocal+Pbody[i]->dt;
        if((EraMin<=Timei)&&(Timei<=EraMax)){
            Pbody[i]->EraLocal += Pbody[i]->dt;
            Pbody[i]->Active = ON;
            Pall.NActives ++;
            if(Pbody[i]->Type == TypeStar){
                Pall.NActivesStars ++;
            }else if(Pbody[i]->Type == TypeSink){
                Pall.NActivesSink ++;
            }else if(Pbody[i]->Type == TypeDM){
                Pall.NActivesDM ++;
            }
        }else{
            Pbody[i]->Active = OFF;
        }
    }

#ifdef USE_SINK_TIMESTEP_LIMITER
    for(int i=0;i<Pall.Nsink;i++){
        if(PsinkBody(i)->Active)
            Psink[i]->dt_localmin = PsinkBody(i)->dt;
            //PbodySink(i)->dt_localmin = Pbody[i]->dt;
    }
#endif

    UpdateTotalActiveNumber();

#ifndef USE_FAST_SCHEME
    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(Pall.NActivesHydro != Pall.NActives)
            fprintf(stderr,"Nhydro %ld != Ngrav %ld\n",Pall.NActivesHydro,Pall.NActives);
    }
#endif
#if 1
    if(Pall.NActivesAll_t == 0){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"[%d] Pall.NActivesAll = %ld, Pall.NActivesAll_t = %ld\n",
                    MPIGetMyID(),Pall.NActivesAll,Pall.NActivesAll_t);
            fprintf(stderr,"[%d] Pall.NActives = %ld, Pall.NActives_t = %ld\n",
                    MPIGetMyID(),Pall.NActives,Pall.NActives_t);
            fprintf(stderr,"[%d] Pall.NActivesDM = %ld, Pall.NActivesDM_t = %ld\n",
                    MPIGetMyID(),Pall.NActivesDM,Pall.NActivesDM_t);
            fprintf(stderr,"[%d] Pall.NActivesHydro = %ld, Pall.NActivesHydro_t = %ld\n",
                    MPIGetMyID(),Pall.NActivesHydro,Pall.NActivesHydro_t);
            fprintf(stderr,"[%d] Pall.NActivesStars = %ld, Pall.NActivesStars_t = %ld\n",
                    MPIGetMyID(),Pall.NActivesStars,Pall.NActivesStars_t);
            fprintf(stderr,"[%d] Pall.NActivesSink = %ld, Pall.NActivesSink_t = %ld\n",
                    MPIGetMyID(),Pall.NActivesSink,Pall.NActivesSink_t);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        assert(Pall.NActivesAll_t > 0);
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
#endif 

#if 1
    if(Pall.dtnow > 0){
        MPI_Barrier(MPI_COMM_WORLD);
        //dprintlmpi(TStepEnd);
        //fflush(NULL);
        assert(TStepEnd < 1.e10);
        //MPI_Abort(MPI_COMM_WORLD,-1);
    }
#endif


    if(MPIGetMyID() == MPI_ROOT_RANK){
#ifdef COSMOLOGICAL_RUN
        fprintf(stderr,"Step %d/%d %.3e Myr, %.3e Myr, t = %.3e Gyr, z = %2.5g, NActives(h,s,sk,d,t) = (%ld,%ld,%ld,%ld,%ld)\n",
            Pall.TStep,TStepEnd,Pall.dtnow*Pall.UnitTime/MEGAYEAR_CGS,
                Pall.EraLocal*Pall.UnitTime/MEGAYEAR_CGS,
                Pall.TCurrent*Pall.UnitTime/GIGAYEAR_CGS,Pall.Redshift,
                    Pall.NActivesHydro_t,Pall.NActivesStars_t,
                    Pall.NActivesSink_t,Pall.NActivesDM_t,Pall.NActivesAll_t);
        fprintf(FpStepLog,"Total Step %u Local Step %d/%d %.3e Myr, %.3e Myr, t = %.3e Gyr, z = %2.5g, NActives(h,s,sk,d,t) = (%ld,%ld,%ld,%ld,%ld)\n",
                Pall.TStepTotal,Pall.TStep,TStepEnd,Pall.dtnow*Pall.UnitTime/MEGAYEAR_CGS,
                Pall.EraLocal*Pall.UnitTime/MEGAYEAR_CGS,
                Pall.TCurrent*Pall.UnitTime/GIGAYEAR_CGS,Pall.Redshift,
                    Pall.NActivesHydro_t,Pall.NActivesStars_t,
                    Pall.NActivesSink_t,Pall.NActivesDM_t,Pall.NActivesAll_t);
#else
#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
        fprintf(stderr,"Step %d/%d %.3e yr, %.3e yr, t = %.3e yr, NActives(h,s,sk,d,t) = (%ld,%ld,%ld,%ld,%ld)\n",
            Pall.TStep,TStepEnd,Pall.dtnow*Pall.UnitTime/YEAR_CGS,
                Pall.EraLocal*Pall.UnitTime/YEAR_CGS,
                Pall.TCurrent*Pall.UnitTime/YEAR_CGS,
                    Pall.NActivesHydro_t,Pall.NActivesStars_t,
                    Pall.NActivesSink_t,Pall.NActivesDM_t,Pall.NActivesAll_t);
        fprintf(FpStepLog,"Total Step %u Local Step %d/%d %.3e yr, %.3e yr, t = %.3e yr, NActives(h,s,sk,d,t) = (%ld,%ld,%ld,%ld,%ld)\n",
                Pall.TStepTotal,Pall.TStep,TStepEnd,Pall.dtnow*Pall.UnitTime/YEAR_CGS,
                Pall.EraLocal*Pall.UnitTime/YEAR_CGS,
                Pall.TCurrent*Pall.UnitTime/YEAR_CGS,
                    Pall.NActivesHydro_t,Pall.NActivesStars_t,
                    Pall.NActivesSink_t,Pall.NActivesDM_t,Pall.NActivesAll_t);
#else
        fprintf(stderr,"Step %d/%d %.3e, %.3e, t = %.3e -> %.3e, NActives(h,s,sk,d,t) = (%ld,%ld,%ld,%ld,%ld)\n",
            Pall.TStep,TStepEnd,Pall.dtnow*Pall.UnitTime,Pall.EraLocal*Pall.UnitTime,
                Pall.TCurrent*Pall.UnitTime,(Pall.TCurrent+dtnow)*Pall.UnitTime,
                    Pall.NActivesHydro_t,Pall.NActivesStars_t,
                        Pall.NActivesSink_t,Pall.NActivesDM_t,Pall.NActivesAll_t);
        fprintf(FpStepLog,"Total Step %u Local Step %d/%d %.3e, %.3e, t = %.3e -> %.3e, NActives(h,s,sk,d,t) = (%ld,%ld,%ld,%ld,%ld)\n",
                Pall.TStepTotal,Pall.TStep,TStepEnd,Pall.dtnow*Pall.UnitTime,Pall.EraLocal*Pall.UnitTime,
                Pall.TCurrent*Pall.UnitTime,(Pall.TCurrent+dtnow)*Pall.UnitTime,
                    Pall.NActivesHydro_t,Pall.NActivesStars_t,
                        Pall.NActivesSink_t,Pall.NActivesDM_t,Pall.NActivesAll_t);
#endif
#endif
    }


#define TINY_TIMESTEP (1.e-15)
    if(Pall.dtnow < TINY_TIMESTEP*Pall.TEnd){
        dtnow = Pbody[0]->dt;
        EraLocalnow = Pbody[0]->EraLocal + Pbody[0]->dt;
        for(int i=1;i<Pall.Ntotal;i++){
            dtnow = fmin(dtnow,Pbody[i]->dt);
            EraLocalnow = fmin(EraLocalnow,Pbody[i]->EraLocal+Pbody[i]->dt);
        }
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"Pall.dtnow < %g*Pall.TEnd\n",TINY_TIMESTEP);
            fprintf(stderr,"%d | dtnow %g EraLocal =%g\n",MPIGetMyID(),dtnow,EraLocalnow);

            fprintf(stderr,"%d | Step dtnow dtnow*step %d %.3e %.3e, TCurrent = %.3e, NActives = %ld\n",
                    MPIGetMyID(), Pall.TStep,Pall.dtnow,Pall.EraLocal,Pall.TCurrent,Pall.NActives_t);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        assert(Pall.dtnow > TINY_TIMESTEP*Pall.TEnd);

        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    Pall.EraLocal = Pall.TStep*Pall.dtnow;
    //if(MPIGetMyID() == MPI_ROOT_RANK){
        //fprintf(stderr,"EraLocal = %g, dtnow = %g\n",Pall.EraLocal,Pall.dtnow);
    //}

    TimingResults.IntegralThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

void AdjustTimeStepAfterFeedback(void){

#ifdef USE_HYDRO_TIMESTEP_LIMITER_FOR_FB //{

    int LimiterDiff = MAX_K_LOCAL-MAX_K_LOCAL_FB;

    // int counter = 0;
    double EraMinimum = Pall.EraLocal+1.1*Pall.dtnow;
    for(int i=0;i<Pall.Nhydro;i++){
        if((Phydro[i]->Active == false)&&(Phydro[i]->DQheat > 0.e0)){
            int new_k_hydro_localmin = Phydro[i]->k_hydro_localmin>>(LimiterDiff);
            if(Phydro[i]->k_hydro - new_k_hydro_localmin > MAX_K_LOCAL_FB){
                Phydro[i]->k_hydro_localmin = new_k_hydro_localmin;

                double dt_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin+MAX_K_LOCAL_FB));
                int step = (int)(Pall.EraLocal/dt_localmin);
                double NextUpdateEra;
                do{
                    NextUpdateEra = step*dt_localmin;
                    step ++;
                }while(NextUpdateEra < EraMinimum);
                Phydro[i]->NextUpdateEra = NextUpdateEra;
                // counter ++;
            }
        }
    }

    // MPI_Allreduce(MPI_IN_PLACE,&counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    // if(MPIGetMyID() == MPI_ROOT_RANK){
        // dprintlmpi(counter);
    // }
#endif // USE_HYDRO_TIMESTEP_LIMITER_FOR_FB //}

    return ;
}


