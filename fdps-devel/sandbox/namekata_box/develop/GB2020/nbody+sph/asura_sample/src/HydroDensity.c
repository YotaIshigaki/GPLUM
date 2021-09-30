#include "config.h"
#include "PlantHydroTree.h"
#include "HydroMisc.h"
#include "NeighborSearch.h"
#include "SizeDetermination.h"
#include "KernelFunctions.h"

#if defined(BAROTROPIC_EOS_RUN)
static inline double __attribute__((always_inline)) ReturnSoundSpeedForBarotoropicRun(const double rho){
    /*
     * P = K \rho^\eta, where K = cs(T=10K)^2
     * Cs(\rho)^2 = K \eta \rho^(\eta-1)
     */
    double K = SQ(SOUND_VELOCITY_FOR_BAROTROPIC_RUN*Pall.UnitTime/Pall.UnitLength);
    if(rho*Pall.ConvertDensityToCGS < CRITICAL_DENSITY_FOR_BAROTROPIC_RUN){
        return (K);
    } else {
        return (1.4*K*pow(rho,0.4));
    }
}
#endif

static void HydroDensityEndProcessing(const int NActives, const int ActiveIndexList[restrict]){ 

    for(int k=0;k<NActives;k++){
        int leaf = ActiveIndexList[k];
        Phydro[leaf]->RhoPred = Phydro[leaf]->Rho;

#ifdef USE_DISPH //{
        Phydro[leaf]->EnergyDensityPred = Phydro[leaf]->EnergyDensity;
        double SmoothedQuantity = Phydro[leaf]->EnergyDensity;
#elif defined(USE_SPSPH)
        Phydro[leaf]->PseudoDensityPred = Phydro[leaf]->PseudoDensity;
        double SmoothedQuantity = Phydro[leaf]->PseudoDensity;
#else // USE_DISPH //}//{
        double SmoothedQuantity = Phydro[leaf]->Rho;
#endif // USE_DISPH //}

        double InvValue = 1.e0/SmoothedQuantity;
        Phydro[leaf]->Div *= InvValue;
        Phydro[leaf]->Rot[0] *= InvValue;
        Phydro[leaf]->Rot[1] *= InvValue;
        Phydro[leaf]->Rot[2] *= InvValue;

#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
        Phydro[leaf]->GradN = -DIMENSION*Phydro[leaf]->NumberDensity/Phydro[leaf]->GradN;
#   ifdef USE_GRAD_H_LIMITER //{
        if(Phydro[leaf]->GradN > GRAD_H_MAX){
            Phydro[leaf]->GradN = GRAD_H_MAX;
        }
        if(Phydro[leaf]->GradN < GRAD_H_MIN){
            Phydro[leaf]->GradN = GRAD_H_MIN;
        }
#   endif // USE_GRAD_H_LIMITER //}
        Phydro[leaf]->fij = -(SmoothedQuantity+Phydro[leaf]->Gradh/DIMENSION)*Phydro[leaf]->GradN/Phydro[leaf]->NumberDensity;


#endif // USE_GRAD_N //}
        Phydro[leaf]->Gradh = -DIMENSION*SmoothedQuantity/Phydro[leaf]->Gradh;
#ifdef USE_GRAD_H_LIMITER
#if 0
        double sign = Sign(Phydro[leaf]->Gradh); 
        double Gradh = fmin(GRAD_H_MAX,fabs(Phydro[leaf]->Gradh));
        Phydro[leaf]->Gradh = sign*fmax(GRAD_H_MIN,Gradh);
#else
        if(Phydro[leaf]->Gradh > GRAD_H_MAX){
            Phydro[leaf]->Gradh = GRAD_H_MAX;
        }
        if(Phydro[leaf]->Gradh < GRAD_H_MIN){
            Phydro[leaf]->Gradh = GRAD_H_MIN;
        }
#endif
#endif

#endif // USE_GRAD_H //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
        Phydro[leaf]->DZdiff = 4.0*Phydro[leaf]->Rho*Phydro[leaf]->Kernel*sqrt(Phydro[leaf]->DZdiff)/(double)Phydro[leaf]->Nlist;
#elif DIFFUSION_TYPE==1 //}//{
        Phydro[leaf]->Sxy *= InvValue;  
        Phydro[leaf]->Sxz *= InvValue; 
        Phydro[leaf]->Syx *= InvValue; 
        Phydro[leaf]->Syz *= InvValue; 
        Phydro[leaf]->Szx *= InvValue; 
        Phydro[leaf]->Szy *= InvValue; 
#endif  // DIFFUSION_TYPE //}
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}

#if 0
#ifdef USE_DISPH //{

        Phydro[leaf]->EnergyDensityPred = Phydro[leaf]->EnergyDensity;
        double InvValue = 1.e0/Phydro[leaf]->EnergyDensity;
        Phydro[leaf]->Div *= InvValue;
        Phydro[leaf]->Rot[0] *= InvValue;
        Phydro[leaf]->Rot[1] *= InvValue;
        Phydro[leaf]->Rot[2] *= InvValue;
#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
        Phydro[leaf]->GradN = -DIMENSION*Phydro[leaf]->NumberDensity/Phydro[leaf]->GradN;
        Phydro[leaf]->fij = -(Phydro[leaf]->EnergyDensity+Phydro[leaf]->Gradh/DIMENSION)*Phydro[leaf]->GradN/Phydro[leaf]->NumberDensity;
#endif // USE_GRAD_N //}
        Phydro[leaf]->Gradh = -DIMENSION*Phydro[leaf]->EnergyDensity/Phydro[leaf]->Gradh;
#endif // USE_GRAD_H //}

#elif defined(USE_SPSPH)

        Phydro[leaf]->PseudoDensityPred = Phydro[leaf]->PseudoDensity;
        double InvValue = 1.e0/Phydro[leaf]->PseudoDensity;
        Phydro[leaf]->Div *= InvValue;
        Phydro[leaf]->Rot[0] *= InvValue;
        Phydro[leaf]->Rot[1] *= InvValue;
        Phydro[leaf]->Rot[2] *= InvValue;
#error Evaluate diffusion coeff here
#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
        Phydro[leaf]->GradN = -DIMENSION*Phydro[leaf]->NumberDensity/Phydro[leaf]->GradN;
        Phydro[leaf]->fij = -(Phydro[leaf]->PseudoDensity+Phydro[leaf]->Gradh/DIMENSION)*Phydro[leaf]->GradN/Phydro[leaf]->NumberDensity;
#endif // USE_GRAD_N //}
        Phydro[leaf]->Gradh = -DIMENSION*Phydro[leaf]->PseudoDensity/Phydro[leaf]->Gradh;
#endif // USE_GRAD_H //}

#else // USE_DISPH //}//{
        double InvValue = 1.e0/Phydro[leaf]->Rho;
        Phydro[leaf]->Div *= InvValue;
        Phydro[leaf]->Rot[0] *= InvValue;
        Phydro[leaf]->Rot[1] *= InvValue;
        Phydro[leaf]->Rot[2] *= InvValue;
#ifdef USE_GRAD_H //{
#ifdef USE_GRAD_N //{
        Phydro[leaf]->GradN = -DIMENSION*Phydro[leaf]->NumberDensity/Phydro[leaf]->GradN;
        Phydro[leaf]->fij = -(Phydro[leaf]->Rho+Phydro[leaf]->Gradh/DIMENSION)*Phydro[leaf]->GradN/Phydro[leaf]->NumberDensity;
#endif // USE_GRAD_N //}
        Phydro[leaf]->Gradh = -DIMENSION*Phydro[leaf]->Rho/Phydro[leaf]->Gradh;
#endif // USE_GRAD_H //}
#endif // USE_DISPH //}
#endif

#if defined(ISOTHERMAL_EOS_RUN)
        double cs = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
        double cs = ReturnSoundSpeedForBarotoropicRun(Phydro[leaf]->RhoPred);
#else
        double cs = sqrt(Pall.GGm1*Phydro[leaf]->UPred);
#endif
        double div = fabs(Phydro[leaf]->Div);
        Phydro[leaf]->F = (div/(div+NORM(Phydro[leaf]->Rot)+0.0001*cs*Phydro[leaf]->KernelPred));
        if(Phydro[leaf]->F > 1.0)
            Phydro[leaf]->F = 1.0;

#ifdef USE_SPSPH //{
        double dt_CFL = 2.0*Phydro[leaf]->KernelPred/cs;
        Phydro[leaf]->Ddif = SPSPH_DIFFUSION_FACTOR*SQ(2.0*Phydro[leaf]->KernelPred)/dt_CFL;
#endif //USE_SPSPH //}

#ifdef USE_VARIABLE_ALPHA //{
        double InvTau = Pall.ViscousL*cs/(2.0*Phydro[leaf]->KernelPred);

        // Phydro[leaf]->DAlpha = -(Phydro[leaf]->Alpha-Pall.ViscousAlphaMin)*InvTau
                        // +Pall.ViscousS*Phydro[leaf]->F*fmax(0.e0,-Phydro[leaf]->Div)
                            // *(Pall.ViscousAlphaMax-Phydro[leaf]->Alpha);
                            //*CUBE(Pall.ViscousAlphaMax-Phydro[leaf]->Alpha);
        Phydro[leaf]->DAlpha = -(Phydro[leaf]->Alpha-Pall.ViscousAlphaMin)*InvTau
                        +fmax(0.e0,-Phydro[leaf]->Div*(Pall.ViscousAlphaMax-Phydro[leaf]->Alpha));
                           //;
#endif // USE_VARIABLE_ALPHA //}
    }

#ifdef HYDRO_TIMESTEP_LIMITER //{
    double EraMinimum = Pall.EraLocal+0.1*Pall.dtnow;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active == false){
            //if(Phydro[i]->k_hydro - Phydro[i]->k_hydro_localmin > MAX_K_LOCAL){
            if(Phydro[i]->k_hydro > Phydro[i]->k_hydro_localmin){
                //double dt_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin+MAX_K_LOCAL));
                double dt_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin));
                int step = (int)(Pall.EraLocal/dt_localmin);
                double NextUpdateEra;
                do{
                    NextUpdateEra = step*dt_localmin;
                    step ++;
                }while(NextUpdateEra < EraMinimum);
                Phydro[i]->NextUpdateEra = NextUpdateEra;
            }
        }
    }
#endif // HYDRO_TIMESTEP_LIMITER //}

    return;
}

static inline bool __attribute__((always_inline)) OverlapDomainDensity(double Pos[restrict], const double h, const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline double __attribute__((always_inline)) DomainDistanceSQR(double Pos[restrict], const int NodeID){

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2);
}


//static inline int __attribute__((always_inline)) CheckHydroDensityExportFlags(const int Index, const int NProcs, bool HydroDensityExportFlags[][NProcs-1]){
static inline int __attribute__((always_inline)) CheckHydroDensityExportFlags(const int Index, bool HydroDensityExportFlags[restrict]){

    if(Pall.Nhydro == 0)
        return 0;

    int DomainID = CommunicationTable[Index].SendRank;
    int NExport = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
        if(HydroNode[CurrentNodeID].NumberofActiveLeaves == 0){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else if( SQ(HydroNode[CurrentNodeID].KernelMax) < DomainDistanceSQR(HydroNode[CurrentNodeID].Pos,DomainID) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
                if(NBCache[leaf].Active){
                    if(OverlapDomainDensity(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,DomainID)){
                        HydroDensityExportFlags[NBCache[leaf].Leaf] = true;
                        //HydroDensityExportFlags[NBCache[leaf].Leaf*NProcs+Index] = true;
                        //HydroDensityExportFlags[NBCache[leaf].Leaf][Index] = ON;
                        NExport ++;
                    }
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}
	return NExport;
}

static struct StructHydroDensityImport InitializeHydroDensityImport = {
    .Rho = 0.e0,
    .Div = 0.e0,
    .Rot = {0.e0,0.e0,0.e0},
#ifdef USE_DISPH
    .EnergyDensity = 0.e0,
#endif // USE_DISPH
#ifdef USE_GRAD_H //{
    .Gradh = 0.e0,
#ifdef USE_GRAD_N //{
    .NumberDensity = 0.e0,
    .GradN = 0.e0,
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
    .DZdiff = 0.e0,
#elif DIFFUSION_TYPE==1 //}/{
    .Sxy = 0.e0,
    .Sxz = 0.e0,
    .Syx = 0.e0,
    .Syz = 0.e0,
    .Szx = 0.e0,
    .Szy = 0.e0,
#endif  // DIFFUSION_TYPE //}
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
#ifdef USE_FUVFEEDBACK //{
    .GradRho = 0.e0,
#endif // USE_FUVFEEDBACK //}
    .Nlist = 0,
    .Leaf = 0,
};


struct StructHydroDensityInput{
    int index;
    int nbindex; 
};


/*
 * This function retruns pieces of density, div, rot and grad-h term.
 * It also returns pieces of Entropy and Pressure^1/Gamma if flags are on.
 * The first argument controls the behavior of this function. If this argument
 * is >= 0, this function returns physical quantities for the local particle of
 * which index is "index". If this argument is == NONE (=-1), this function
 * evaluates physical quantities for the imported particle data.
 */
static struct StructHydroDensityImport ReturnStructureDensityDivRotEngine(
        const struct StructHydroDensityInput Input, 
            struct StructHydroDensityExport HydroDensityExportRecv, 
                bool HydroInteractionFlags[restrict]){

    double Pos[3];
    double Vel[3];
    double Kerneli;
    double B[3][3];
    int k_hydro;

    if(Input.index == NONE){
        Pos[0] = HydroDensityExportRecv.Pos[0];
        Pos[1] = HydroDensityExportRecv.Pos[1];
        Pos[2] = HydroDensityExportRecv.Pos[2];
        Vel[0] = HydroDensityExportRecv.Vel[0];
        Vel[1] = HydroDensityExportRecv.Vel[1];
        Vel[2] = HydroDensityExportRecv.Vel[2];
        Kerneli = HydroDensityExportRecv.Kernel;
#if VISCOSITY_TYPE == 1 //{
        B[0][0] = HydroDensityExportRecv.B[0][0];
        B[0][1] = HydroDensityExportRecv.B[0][1];
        B[0][2] = HydroDensityExportRecv.B[0][2];
        B[1][0] = HydroDensityExportRecv.B[1][0];
        B[1][1] = HydroDensityExportRecv.B[1][1];
        B[1][2] = HydroDensityExportRecv.B[1][2];
        B[2][0] = HydroDensityExportRecv.B[2][0];
        B[2][1] = HydroDensityExportRecv.B[2][1];
        B[2][2] = HydroDensityExportRecv.B[2][2];
#endif // VISCOSITY_TYPE == 1 //}
#ifdef HYDRO_TIMESTEP_LIMITER
        k_hydro = HydroDensityExportRecv.k_hydro;
#endif // HYDRO_TIMESTEP_LIMITER
    } else {
        Pos[0] = PhydroPosP(Input.index)[0];
        Pos[1] = PhydroPosP(Input.index)[1];
        Pos[2] = PhydroPosP(Input.index)[2];
        Vel[0] = Phydro[Input.index]->VelP[0];
        Vel[1] = Phydro[Input.index]->VelP[1];
        Vel[2] = Phydro[Input.index]->VelP[2];
        Kerneli = Phydro[Input.index]->KernelPred;
#if VISCOSITY_TYPE == 1 //{
        B[0][0] = Phydro[Input.index]->Bxx;
        B[0][1] = Phydro[Input.index]->Bxy;
        B[0][2] = Phydro[Input.index]->Bxz;
        B[1][0] = Phydro[Input.index]->Byx;
        B[1][1] = Phydro[Input.index]->Byy;
        B[1][2] = Phydro[Input.index]->Byz;
        B[2][0] = Phydro[Input.index]->Bzx;
        B[2][1] = Phydro[Input.index]->Bzy;
        B[2][2] = Phydro[Input.index]->Bzz;
#endif // VISCOSITY_TYPE == 1 //}
#ifdef HYDRO_TIMESTEP_LIMITER
        k_hydro = Phydro[Input.index]->k_hydro;
#endif // HYDRO_TIMESTEP_LIMITER
    }

    double InvKerneli = 1.e0/Kerneli;
    int Neighbors[MaxNeighborSize];
    struct StructHydroDensityImport TempHydroDensityImport = InitializeHydroDensityImport;

    double TimeComm = GetElapsedTime();
#ifdef USE_NEIGHBOR_LIST //{
    int Nlist;
    if(Input.index == NONE){
        Nlist = GetNeighborsLimited(Pos,2.e0*Kerneli,Neighbors);
    } else {
        struct StructGetLocalNeighborList NB = GetLocalNeighrborList(Input.nbindex);
        if(NB.Nlist>0){
            Nlist = NB.Nlist;
            for(int k=0;k<Nlist;k++){
                Neighbors[k] = NB.Neighbors[k];
            }
        } else {
            Nlist = GetNeighborsLimited(Pos,2.e0*Kerneli,Neighbors);
        }
    }
#else //USE_NEIGHBOR_LIST //}//{
    int Nlist = GetNeighborsLimited(Pos,2.e0*Kerneli,Neighbors);
#endif //}
    TimingResults.HydroDensityNeighborSearchThisStep += GetElapsedTime()-TimeComm;

    struct StructHydroDensityPreFetch{
        double Mass;
        double Kernel;
        double Pos[3];
        double Vel[3];
        double Weight; 
#if VISCOSITY_TYPE == 1 //{
        double Volume; 
        double B[3][3];
#endif // VISCOSITY_TYPE == 1 //}
#ifdef USE_SMOOTHED_POTENTIAL //{
        double Pot;
#endif // USE_SMOOTHED_POTENTIAL //}
    } HydroDensityPreFetch[Nlist];

	for(int k=0;k<Nlist;k++){
		int leaf = Neighbors[k];
        HydroDensityPreFetch[k].Mass = PhydroMass(leaf);
        HydroDensityPreFetch[k].Kernel = Phydro[leaf]->KernelPred;
        HydroDensityPreFetch[k].Pos[0] = PhydroPosP(leaf)[0];
        HydroDensityPreFetch[k].Pos[1] = PhydroPosP(leaf)[1];
        HydroDensityPreFetch[k].Pos[2] = PhydroPosP(leaf)[2];
        HydroDensityPreFetch[k].Vel[0] = Phydro[leaf]->VelP[0];
        HydroDensityPreFetch[k].Vel[1] = Phydro[leaf]->VelP[1];
        HydroDensityPreFetch[k].Vel[2] = Phydro[leaf]->VelP[2];
#ifdef USE_DISPH //{
        HydroDensityPreFetch[k].Weight = Phydro[leaf]->Mass*Phydro[leaf]->UPred;
#elif defined(USE_SPSPH) // USE_DISPH//USE_SPSPH //}//{
        HydroDensityPreFetch[k].Weight = Phydro[leaf]->ZwPred;
#else // USE_SPSPH //}//{
        HydroDensityPreFetch[k].Weight = Phydro[leaf]->Mass;
#endif // USE_SPSPH //}
#if VISCOSITY_TYPE == 1 //{
#ifdef USE_DISPH //{
        HydroDensityPreFetch[k].Volume = Phydro[leaf]->Mass*Phydro[leaf]->UPred/Phydro[leaf]->EnergyDensityPred;
#elif defined(USE_SPSPH) // USE_DISPH//USE_SPSPH //}//{
        HydroDensityPreFetch[k].Volume = Phydro[leaf]->ZwPred/Phydro[leaf]->PseudoDensityPred;
#else // USE_SPSPH //}//{
        HydroDensityPreFetch[k].Volume = Phydro[leaf]->Mass/Phydro[leaf]->Rho;
#endif // USE_SPSPH //}
#endif // VISCOSITY_TYPE == 1 //}
//#ifdef HYDRO_TIMESTEP_LIMITER //{
        //Phydro[leaf]->k_hydro_localmin = MIN(Phydro[leaf]->k_hydro_localmin,k_hydro);
//#endif // HYDRO_TIMESTEP_LIMITER //}
#ifdef USE_SMOOTHED_POTENTIAL //{
        HydroDensityPreFetch[k].Pot = PhydroBody(leaf)->Pot/PhydroBody(leaf)->Mass;
#endif // USE_SMOOTHED_POTENTIAL //}
        if(Input.index == NONE)
            HydroInteractionFlags[leaf] = true;
    }

	for(int k=0;k<Nlist;k++){
        double xij[3],vij[3];
#ifdef PERIODIC_RUN //{
        xij[0] = PeriodicDistance(Pos[0],HydroDensityPreFetch[k].Pos[0],0);
        xij[1] = PeriodicDistance(Pos[1],HydroDensityPreFetch[k].Pos[1],1);
        xij[2] = PeriodicDistance(Pos[2],HydroDensityPreFetch[k].Pos[2],2);
#else // PERIODIC_RUN //}//{
        xij[0] = Pos[0]-HydroDensityPreFetch[k].Pos[0];
        xij[1] = Pos[1]-HydroDensityPreFetch[k].Pos[1];
        xij[2] = Pos[2]-HydroDensityPreFetch[k].Pos[2];
#endif // PERIODIC_RUN //}

        vij[0] = Vel[0]-HydroDensityPreFetch[k].Vel[0];
        vij[1] = Vel[1]-HydroDensityPreFetch[k].Vel[1];
        vij[2] = Vel[2]-HydroDensityPreFetch[k].Vel[2];

        double vx = DOT_PRODUCT(xij,vij);
        double r2 = SQ(xij[0])+SQ(xij[1])+SQ(xij[2]);
        double r = sqrt(r2); 
		double w = SPHKernel(r,InvKerneli);

        TempHydroDensityImport.Rho += HydroDensityPreFetch[k].Mass*w;

#ifdef USE_DISPH //{
        TempHydroDensityImport.EnergyDensity += HydroDensityPreFetch[k].Weight*w;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
        TempHydroDensityImport.PseudoDensity += HydroDensityPreFetch[k].Weight*w;
#endif // USE_SPSPH //}

        double dw = dSPHKernel(r,InvKerneli);
        double Weightdw = HydroDensityPreFetch[k].Weight*dw;
#ifdef USE_GRAD_H //{
        TempHydroDensityImport.Gradh += Weightdw*r2;
#ifdef USE_GRAD_N //{
        TempHydroDensityImport.NumberDensity += w;
        TempHydroDensityImport.GradN += r2*dw;
#endif //USE_GRAD_N //}
#endif //USE_GRAD_H //}

#if VISCOSITY_TYPE == 0 //{
        TempHydroDensityImport.Div -= Weightdw*vx;
        TempHydroDensityImport.Rot[0] -= Weightdw*(vij[2]*xij[1]-vij[1]*xij[2]);
        TempHydroDensityImport.Rot[1] -= Weightdw*(vij[0]*xij[2]-vij[2]*xij[0]);
        TempHydroDensityImport.Rot[2] -= Weightdw*(vij[1]*xij[0]-vij[0]*xij[1]);
#elif VISCOSITY_TYPE == 1 //{
        double psi_ji[3] = {
            B[0][0]*xij[0]+B[0][1]*xij[1]+B[0][2]*xij[2],
            B[1][0]*xij[0]+B[1][1]*xij[1]+B[1][2]*xij[2],
            B[2][0]*xij[0]+B[2][1]*xij[1]+B[2][2]*xij[2]};

        psi_ji[0] *= w*HydroDensityPreFetch[k].Volume;
        psi_ji[1] *= w*HydroDensityPreFetch[k].Volume;
        psi_ji[2] *= w*HydroDensityPreFetch[k].Volume;

        TempHydroDensityImport.Div += DOT_PRODUCT(vij, psi_ji);

        TempHydroDensityImport.Rot[0] += vij[1]*psi_ji[2]-vij[2]*psi_ji[1];
        TempHydroDensityImport.Rot[1] += vij[2]*psi_ji[0]-vij[0]*psi_ji[2];
        TempHydroDensityImport.Rot[2] += vij[0]*psi_ji[1]-vij[1]*psi_ji[0];
#endif // VISCOSITY_TYPE == 1 //}

#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE == 0 //{
        TempHydroDensityImport.DZdiff += NORM2(vij);
#elif DIFFUSION_TYPE == 1 //} //{
        TempHydroDensityImport.Sxy += Weightdw*vij[1]*xij[0];
        TempHydroDensityImport.Sxz += Weightdw*vij[2]*xij[0];
        TempHydroDensityImport.Syx += Weightdw*vij[2]*xij[1];
        TempHydroDensityImport.Syz += Weightdw*vij[2]*xij[1];
        TempHydroDensityImport.Szx += Weightdw*vij[0]*xij[2];
        TempHydroDensityImport.Szy += Weightdw*vij[1]*xij[2];
#endif //DIFFUSION_TYPE //} 
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}

#ifdef USE_FUVFEEDBACK //{
        TempHydroDensityImport.GradRho += HydroDensityPreFetch[k].Mass*r*dw;
#endif // USE_FUVFEEDBACK //}
#ifdef USE_SMOOTHED_POTENTIAL //{
        TempHydroDensityImport.Pot += HydroDensityPreFetch[k].Pot*HydroDensityPreFetch[k].Volume*w;
#endif // USE_SMOOTHED_POTENTIAL //}
#ifdef HYDRO_TIMESTEP_LIMITER //{
        int k_shift = r*InvKerneli>1.0?1:0;
        //int k_shift = 0;
        Phydro[Neighbors[k]]->k_hydro_localmin = MIN(Phydro[Neighbors[k]]->k_hydro_localmin,k_hydro+k_shift+MAX_K_LOCAL);
#endif // HYDRO_TIMESTEP_LIMITER //}

    }
    TempHydroDensityImport.Nlist += Nlist;

    return TempHydroDensityImport;
}

void CalcDensityDivRot(void){

    double TimingResultThisRoutine = GetElapsedTime();

    int NProcs = MPIGetNumProcs();
    int MyID = MPIGetMyID();


    MPI_Status  mpi_status;

    static int HydroDensityExportFlagsMaxAllocated = 0;
    static bool *HydroDensityExportFlags = NULL;
    static int *ActiveIndexList = NULL;
    if(HydroDensityExportFlagsMaxAllocated < MAX(Pall.Nhydro,NAdditionUnit)){
        if(HydroDensityExportFlagsMaxAllocated > 0){
            free(ActiveIndexList);
            free(HydroDensityExportFlags);
        }
        HydroDensityExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*Pall.Nhydro,NAdditionUnit));
        HydroDensityExportFlags = malloc(sizeof(bool)*HydroDensityExportFlagsMaxAllocated);
        ActiveIndexList = malloc(sizeof(int)*HydroDensityExportFlagsMaxAllocated);
    }

    int NActives = 0;
    int RootNodeID = 0;
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header + i;
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            ActiveIndexList[NActives] = NBCache[leaf].Leaf;
            // for(int k=0;k<NProcs-1;k++)
              //   HydroDensityExportFlags[NBCache[leaf].Leaf*NProcs+k] = false;
            NActives ++;
        }
    }

    // counter for export and import. 
    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];
    int SendFlag,RecvFlag;

    struct StructHydroDensityExport *HydroDensityExportSend[NProcs];

    for(int i=0;i<NProcs-1;i++){
        MPI_Irecv(NImportThisTime+i,1,MPI_INT,
            CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_EXPORT+i,
                MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        MPI_Test(mpi_request_Export_Recv+i,&RecvFlag,MPI_STATUS_IGNORE);
    }


    int NExportMax = 0;
    for(int i=0;i<NProcs-1;i++){
        //memset(HydroDensityExportFlags,0,Pall.Nhydro*sizeof(bool));
        for(int k=0;k<Pall.Nhydro;k++) HydroDensityExportFlags[k] = false;
        NExportThisTime[i] = CheckHydroDensityExportFlags(i,HydroDensityExportFlags);
        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructHydroDensityExport),i);
        HydroDensityExportSend[i] = BufferExportSend[i];

        MPI_Isend(NExportThisTime+i,1,MPI_INT,
                CommunicationTable[i].SendRank,TAG_SPH_DENSITY_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Test(mpi_request_Export_Send+i,&SendFlag,MPI_STATUS_IGNORE);


        int NExport = 0;
        if(NExportThisTime[i] > 0){
        for(int k=0;k<NActives;k++){
            int leaf = ActiveIndexList[k];
            //if(HydroDensityExportFlags[leaf][i]){ 
            //if(HydroDensityExportFlags[leaf*NProcs+i]){ 
            if(HydroDensityExportFlags[leaf]){ 
                HydroDensityExportSend[i][NExport].Kernel = Phydro[leaf]->KernelPred;
                HydroDensityExportSend[i][NExport].Pos[0] = PhydroPosP(leaf)[0];
                HydroDensityExportSend[i][NExport].Pos[1] = PhydroPosP(leaf)[1];
                HydroDensityExportSend[i][NExport].Pos[2] = PhydroPosP(leaf)[2];
                HydroDensityExportSend[i][NExport].Vel[0] = Phydro[leaf]->VelP[0];
                HydroDensityExportSend[i][NExport].Vel[1] = Phydro[leaf]->VelP[1];
                HydroDensityExportSend[i][NExport].Vel[2] = Phydro[leaf]->VelP[2];
                HydroDensityExportSend[i][NExport].Leaf = leaf;
#if VISCOSITY_TYPE == 1 //{
                HydroDensityExportSend[i][NExport].B[0][0] = Phydro[leaf]->Bxx;
                HydroDensityExportSend[i][NExport].B[0][1] = Phydro[leaf]->Bxy;
                HydroDensityExportSend[i][NExport].B[0][2] = Phydro[leaf]->Bxz;
                HydroDensityExportSend[i][NExport].B[1][0] = Phydro[leaf]->Byx;
                HydroDensityExportSend[i][NExport].B[1][1] = Phydro[leaf]->Byy;
                HydroDensityExportSend[i][NExport].B[1][2] = Phydro[leaf]->Byz;
                HydroDensityExportSend[i][NExport].B[2][0] = Phydro[leaf]->Bzx;
                HydroDensityExportSend[i][NExport].B[2][1] = Phydro[leaf]->Bzy;
                HydroDensityExportSend[i][NExport].B[2][2] = Phydro[leaf]->Bzz;
#endif // VISCOSITY_TYPE == 1 //}
#ifdef HYDRO_TIMESTEP_LIMITER //{
                HydroDensityExportSend[i][NExport].k_hydro = Phydro[leaf]->k_hydro;
#endif // HYDRO_TIMESTEP_LIMITER //}
#ifdef USE_DEBUG_MODE //{
                HydroDensityExportSend[i][NExport].GlobalID = PhydroBody(leaf)->GlobalID;
#endif // USE_DEBUG_MODE //}
#ifdef USE_SMOOTHED_POTENTIAL //{
                HydroDensityExportSend[i][NExport].Pot = PhydroBody(leaf)->Pot/PhydroBody(leaf)->Mass;
#endif // USE_SMOOTHED_POTENTIAL //{
                NExport ++;
            }
        }
        }
        // assert(NExport == NExportThisTime[i]);
        NExportMax = MAX(NExport,NExportMax);
    }

    double TimeComm1 = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.HydroAccCommThisStep += GetElapsedTime()-TimeComm1;

    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportAll += NImportThisTime[i];
    }

    struct StructHydroDensityExport *HydroDensityExportRecv;
    struct StructHydroDensityImport *HydroDensityImportSend;
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHydroDensityExport));
    CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructHydroDensityImport));
    HydroDensityExportRecv = BufferExportRecv;
    HydroDensityImportSend = BufferImportSend;

    int NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(HydroDensityExportSend[i],
                NExportThisTime[i]*sizeof(struct StructHydroDensityExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(HydroDensityExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructHydroDensityExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }

    for(int i=0;i<NActives;i++){  // local
        int leaf = ActiveIndexList[i];

        struct StructHydroDensityImport TempHydroDensityImport =
            ReturnStructureDensityDivRotEngine((struct StructHydroDensityInput){.index=leaf,.nbindex=i},(struct StructHydroDensityExport){0.0},&(bool){true});
        Phydro[leaf]->Nlist = TempHydroDensityImport.Nlist;
        Phydro[leaf]->Rho = TempHydroDensityImport.Rho;
        Phydro[leaf]->Div = TempHydroDensityImport.Div;
        Phydro[leaf]->Rot[0] = TempHydroDensityImport.Rot[0];
        Phydro[leaf]->Rot[1] = TempHydroDensityImport.Rot[1];
        Phydro[leaf]->Rot[2] = TempHydroDensityImport.Rot[2];
#ifdef USE_DISPH //{
        Phydro[leaf]->EnergyDensity = TempHydroDensityImport.EnergyDensity;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
        Phydro[leaf]->PseudoDensity = TempHydroDensityImport.PseudoDensity;
#endif // USE_SPSPH //}

#ifdef USE_GRAD_H //{
        Phydro[leaf]->Gradh = TempHydroDensityImport.Gradh;
#ifdef USE_GRAD_N //{
        Phydro[leaf]->NumberDensity = TempHydroDensityImport.NumberDensity;
        Phydro[leaf]->GradN = TempHydroDensityImport.GradN;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}

#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
        Phydro[leaf]->DZdiff = TempHydroDensityImport.DZdiff; 
#elif DIFFUSION_TYPE==1 //} //{
        Phydro[leaf]->Sxy = TempHydroDensityImport.Sxy; 
        Phydro[leaf]->Sxz = TempHydroDensityImport.Sxz; 
        Phydro[leaf]->Syx = TempHydroDensityImport.Syx; 
        Phydro[leaf]->Syz = TempHydroDensityImport.Syz; 
        Phydro[leaf]->Szx = TempHydroDensityImport.Szx; 
        Phydro[leaf]->Szy = TempHydroDensityImport.Szy; 
#endif  // DIFFUSION_TYPE //}
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
#ifdef USE_FUVFEEDBACK //{
        Phydro[leaf]->GradRho = TempHydroDensityImport.GradRho; 
#endif // USE_FUVFEEDBACK //}
#ifdef USE_SMOOTHED_POTENTIAL //{
        Phydro[leaf]->Pot = TempHydroDensityImport.Pot; 
#endif // USE_SMOOTHED_POTENTIAL //}
    } 
    

    double TimeComm = GetElapsedTime();
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.HydroDensityCommThisStep += GetElapsedTime()-TimeComm;

    // allocated buffer hydro interaction flags
    CheckSizeofBufferHydroInteractionFlags(Pall.Nhydro);  
    bool *HydroInteractionFlags;
    HydroInteractionFlags = (void *)BufferHydroInteractionFlags;

    int HydroInteractionFlagsLoopSize = Pall.Nhydro*(NProcs-1);
    for(int i=0;i<HydroInteractionFlagsLoopSize;i++){
        HydroInteractionFlags[i] = false;
    }

    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        for(int k=0;k<NImportThisTime[i];k++){
            HydroDensityImportSend[NImportAll] = 
                ReturnStructureDensityDivRotEngine(
                    (struct StructHydroDensityInput){.index=NONE},
                        HydroDensityExportRecv[NImportAll],HydroInteractionFlags+Pall.Nhydro*i);
            HydroDensityImportSend[NImportAll].Leaf = HydroDensityExportRecv[NImportAll].Leaf;

            NImportAll ++;
        }
    }

    struct StructHydroDensityImport *HydroDensityImportRecv[NProcs];
    CheckSizeofBufferImportRecv(NExportThisTime,sizeof(struct StructHydroDensityImport));
    for(int i=0;i<NProcs-1;i++)
        HydroDensityImportRecv[i] = BufferImportRecv[i];

    NImport = 0;
    counter_send = counter_recv = 0;
    for(int i=0;i<NProcs-1;i++){
        if(NImportThisTime[i]>0){
            MPI_Isend(HydroDensityImportSend+NImport,
                NImportThisTime[i]*sizeof(struct StructHydroDensityImport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_IMPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
        }
        if(NExportThisTime[i]>0){
            MPI_Irecv(HydroDensityImportRecv[i],
                NExportThisTime[i]*sizeof(struct StructHydroDensityImport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_IMPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
        }
        NImport += NImportThisTime[i];
    }

    TimeComm = GetElapsedTime();
    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
    TimingResults.HydroDensityCommThisStep += GetElapsedTime()-TimeComm;

    for(int i=0;i<NProcs-1;i++){
        for(int k=0;k<NExportThisTime[i];k++){
            if(HydroDensityImportRecv[i][k].Nlist==0) continue;
            int leaf = HydroDensityImportRecv[i][k].Leaf;
            Phydro[leaf]->Nlist += HydroDensityImportRecv[i][k].Nlist;
            Phydro[leaf]->Rho += HydroDensityImportRecv[i][k].Rho;
            Phydro[leaf]->Div += HydroDensityImportRecv[i][k].Div;
            Phydro[leaf]->Rot[0] += HydroDensityImportRecv[i][k].Rot[0];
            Phydro[leaf]->Rot[1] += HydroDensityImportRecv[i][k].Rot[1];
            Phydro[leaf]->Rot[2] += HydroDensityImportRecv[i][k].Rot[2];
#ifdef USE_DISPH //{
            Phydro[leaf]->EnergyDensity += HydroDensityImportRecv[i][k].EnergyDensity;
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
            Phydro[leaf]->PseudoDensity += HydroDensityImportRecv[i][k].PseudoDensity;
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
            Phydro[leaf]->Gradh += HydroDensityImportRecv[i][k].Gradh;
#ifdef USE_GRAD_N //{
            Phydro[leaf]->NumberDensity += HydroDensityImportRecv[i][k].NumberDensity;
            Phydro[leaf]->GradN += HydroDensityImportRecv[i][k].GradN;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}

#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
            Phydro[leaf]->DZdiff += HydroDensityImportRecv[i][k].DZdiff; 
#elif DIFFUSION_TYPE==1 //}//{
            Phydro[leaf]->Sxy += HydroDensityImportRecv[i][k].Sxy; 
            Phydro[leaf]->Sxz += HydroDensityImportRecv[i][k].Sxz; 
            Phydro[leaf]->Syx += HydroDensityImportRecv[i][k].Syx; 
            Phydro[leaf]->Syz += HydroDensityImportRecv[i][k].Syz; 
            Phydro[leaf]->Szx += HydroDensityImportRecv[i][k].Szx; 
            Phydro[leaf]->Szy += HydroDensityImportRecv[i][k].Szy; 
#endif  // DIFFUSION_TYPE //}
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
#ifdef USE_FUVFEEDBACK //{
           Phydro[leaf]->GradRho += HydroDensityImportRecv[i][k].GradRho; 
#endif // USE_FUVFEEDBACK //}
#ifdef USE_SMOOTHED_POTENTIAL //{
           Phydro[leaf]->Pot += HydroDensityImportRecv[i][k].Pot; 
#endif // USE_SMOOTHED_POTENTIAL //}
        }
    }

    HydroDensityEndProcessing(NActives,ActiveIndexList);

    TimingResults.HydroDensityThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}

struct StructBinarySearch {
    double Rvalue;
    double Lvalue;
};

static inline int CheckNeighborNumberAndUpdateKernelRho(const int NActives, const int NProcs, 
        bool HydroDensityExportFlags[restrict][NProcs], int ActiveIndexList[restrict], 
            struct StructBinarySearch BinarySearch[restrict])  __attribute__((always_inline));
static inline int CheckNeighborNumberAndUpdateKernelRho(const int NActives, const int NProcs, 
        bool HydroDensityExportFlags[restrict][NProcs], int ActiveIndexList[restrict], 
            struct StructBinarySearch BinarySearch[restrict]){

#define KernelFactInc   (1.14) // 1.5 ^ (1.3)
#define KernelFactDec   (0.79) // 0.75 ^ (1.3)

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#ifdef USE_MAXIMUM_KERNEL_SIZE
#ifdef MAXIMUM_KERNEL_SIZE
    double MaximumKernelSize = Pall.AdaptiveSofteningFactor*MAXIMUM_KERNEL_SIZE*KPC_CGS/Pall.UnitLength;
#else
#error Set MAXIMUM_KERNEL_SIZE
#endif
#endif

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int leaf = NBCache[ActiveIndexList[i]].Leaf;
        if(HydroDensityExportFlags[leaf][NProcs-1]){ 
            int Nlist = Phydro[leaf]->Nlist;
            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                HydroDensityExportFlags[leaf][NProcs-1] = OFF;
            }else if((BinarySearch[i].Rvalue>0.e0)&&(BinarySearch[i].Lvalue>0.e0)){
                if(BinarySearch[i].Rvalue-BinarySearch[i].Lvalue < 1.e-3*BinarySearch[i].Lvalue)
                    HydroDensityExportFlags[leaf][NProcs-1] = OFF;
            }
#ifdef USE_MINIMUM_KERNEL_SIZE
            else if (((NBmax)<Nlist)&&(Phydro[leaf]->KernelPred<=0.5*PhydroBody(leaf)->Eps*Pall.AdaptiveSofteningFactor)){
                Phydro[leaf]->KernelPred = 0.5*PhydroBody(leaf)->Eps*Pall.AdaptiveSofteningFactor;
                HydroDensityExportFlags[leaf][NProcs-1] = OFF;
            }
#endif
#ifdef USE_MAXIMUM_KERNEL_SIZE
            else if (((NBmin)>Nlist)&&(Phydro[leaf]->KernelPred>MaximumKernelSize)){
                Phydro[leaf]->Kernel =
                Phydro[leaf]->KernelPred = MaximumKernelSize;
                HydroDensityExportFlags[leaf][NProcs-1] = OFF;
            }
#endif
            if(HydroDensityExportFlags[leaf][NProcs-1]){
                if(Nlist<NBmin){
                    BinarySearch[i].Lvalue = fmax(BinarySearch[i].Lvalue,Phydro[leaf]->KernelPred);
                } else if(Nlist>NBmax){
                    if(BinarySearch[i].Rvalue > 0.e0){
                        BinarySearch[i].Rvalue = fmin(BinarySearch[i].Rvalue,Phydro[leaf]->KernelPred);
                    }else{
                        BinarySearch[i].Rvalue = Phydro[leaf]->KernelPred;
                    }
                }

                if((BinarySearch[i].Lvalue>0.e0)&&(BinarySearch[i].Rvalue>0.e0)){
                    Phydro[leaf]->Kernel =
                    Phydro[leaf]->KernelPred = cbrt(0.5*(CUBE(BinarySearch[i].Lvalue)+CUBE(BinarySearch[i].Rvalue)));
                }else{
                    if((BinarySearch[i].Rvalue == 0.e0)&&(BinarySearch[i].Lvalue > 0.e0)){
                        Phydro[leaf]->Kernel =
                        Phydro[leaf]->KernelPred *= KernelFactInc;
                    }else if((BinarySearch[i].Rvalue > 0.e0)&&(BinarySearch[i].Lvalue == 0.e0)){
                        Phydro[leaf]->Kernel =
                        Phydro[leaf]->KernelPred *= KernelFactDec;
                    }
                }
                NLocalActiveLeaves ++;
            }
        }
    }

    return NLocalActiveLeaves;
}

static inline int CheckHydroKernelDensityExportFlags(const int Index, const int NProcs, 
        bool HydroKernelDensityExportFlags[][NProcs]) __attribute__((always_inline));
static inline int CheckHydroKernelDensityExportFlags(const int Index, const int NProcs, 
        bool HydroKernelDensityExportFlags[][NProcs]){ 

    int NodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    while(CurrentNodeID != RootNodeID){
        if(HydroNode[CurrentNodeID].NumberofActiveLeaves == 0){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else if( SQ(HydroNode[CurrentNodeID].KernelMax) < DomainDistanceSQR(HydroNode[CurrentNodeID].Pos,NodeID) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
            int Number_of_leaf = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<Number_of_leaf;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
                if(NBCache[leaf].Active){
                    if(HydroKernelDensityExportFlags[NBCache[leaf].Leaf][NProcs-1]){
                        if(OverlapDomainDensity(NBCache[leaf].Pos,2.0*NBCache[leaf].Kernel,NodeID)){
                            HydroKernelDensityExportFlags[NBCache[leaf].Leaf][Index] = ON;
                            NExport ++;
                        }
                    }
                }
            }
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
    }
    
	return NExport;
}


static struct StructHydroDensityImport ReturnStructureDensityDivRotEngineDirect(
        double Pos[restrict], double Vel[restrict], const double Kerneli){

    double InvKerneli = 1.e0/Kerneli;
    struct StructHydroDensityImport TempHydroDensityImport = InitializeHydroDensityImport;

	for(int k=0;k<Pall.Nhydro;k++){
		int leaf = k;
        double xij[3],vij[3];
        xij[0] = Pos[0]-PhydroPosP(leaf)[0];
        xij[1] = Pos[1]-PhydroPosP(leaf)[1];
        xij[2] = Pos[2]-PhydroPosP(leaf)[2];

        vij[0] = Vel[0]-Phydro[leaf]->VelP[0];
        vij[1] = Vel[1]-Phydro[leaf]->VelP[1];
        vij[2] = Vel[2]-Phydro[leaf]->VelP[2];

        double vx = DOT_PRODUCT(xij,vij);
        double r2 = SQ(xij[0])+SQ(xij[1])+SQ(xij[2]);
        double r = sqrt(r2); 

	    double u = r*InvKerneli;
        if(u > 2.0) continue;
        TempHydroDensityImport.Nlist ++;

		double w = SPHKernel(r,InvKerneli);
        TempHydroDensityImport.Rho += PhydroMass(leaf)*w;
        if(r>0.e0){
		    double dw = dSPHKernel(r,InvKerneli);
            double mdwj = PhydroMass(leaf)*dw;
            TempHydroDensityImport.Div -= mdwj*vx;
            TempHydroDensityImport.Rot[0] -= mdwj*(vij[2]*xij[1]-vij[1]*xij[2]);
            TempHydroDensityImport.Rot[1] -= mdwj*(vij[0]*xij[2]-vij[2]*xij[0]);
            TempHydroDensityImport.Rot[2] -= mdwj*(vij[1]*xij[0]-vij[0]*xij[1]);
#ifdef USE_FUVFEEDBACK //{
            TempHydroDensityImport.GradRho += PhydroMass(leaf)*dw;
#endif // USE_FUVFEEDBACK //}
        }
    }

    return TempHydroDensityImport;
}

void CalcDensityDivRotDirect(void){

    int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;


    // count Actives
    int Actives = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active)
            Actives ++;
    }
    int WActives;
    MPI_Allreduce(&Actives,&WActives,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];
    struct StructHydroDensityExport *HydroDensityExportSend[NProcs];

    for(int i=0;i<NProcs-1;i++){ // put data into output
        CheckSizeofBufferExportSendIndex(Actives,sizeof(struct StructHydroDensityExport),i);
        HydroDensityExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        for(int k=0;k<Pall.Nhydro;k++){
            if(Phydro[k]->Active == true){
                int leaf = k;
                HydroDensityExportSend[i][NExport].Kernel = Phydro[leaf]->KernelPred;
                HydroDensityExportSend[i][NExport].Pos[0] = PhydroPosP(leaf)[0];
                HydroDensityExportSend[i][NExport].Pos[1] = PhydroPosP(leaf)[1];
                HydroDensityExportSend[i][NExport].Pos[2] = PhydroPosP(leaf)[2];
                HydroDensityExportSend[i][NExport].Vel[0] = Phydro[leaf]->VelP[0];
                HydroDensityExportSend[i][NExport].Vel[1] = Phydro[leaf]->VelP[1];
                HydroDensityExportSend[i][NExport].Vel[2] = Phydro[leaf]->VelP[2];
                HydroDensityExportSend[i][NExport].Leaf = leaf;
                NExport ++;
            }
        }
        NExportThisTime[i] = NExport;
    }


    int NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){ 
        MPI_Sendrecv(NExportThisTime+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_PRECOMM,
            NImportThisTime+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
        NImportAll += NImportThisTime[i];
    }

    struct StructHydroDensityExport *HydroDensityExportRecv;
    struct StructHydroDensityImport *HydroDensityImportSend;
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHydroDensityExport));
    CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructHydroDensityImport));
    HydroDensityExportRecv = BufferExportRecv;
    HydroDensityImportSend = BufferImportSend;

    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){ 
        MPI_Isend(HydroDensityExportSend[i],
            NExportThisTime[i]*sizeof(struct StructHydroDensityExport),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(HydroDensityExportRecv+NImport,
            NImportThisTime[i]*sizeof(struct StructHydroDensityExport),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_EXPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTime[i];
    }
    assert(NImport == NImportAll);

    for(int k=0;k<Pall.Nhydro;k++){ 
        if(Phydro[k]->Active == true){
            int leaf = k;
            struct StructHydroDensityImport TempHydroDensityImport =
                ReturnStructureDensityDivRotEngineDirect(
                    PhydroPosP(leaf),Phydro[leaf]->VelP,Phydro[leaf]->KernelPred);

            Phydro[leaf]->Nlist = TempHydroDensityImport.Nlist;
            Phydro[leaf]->Rho = TempHydroDensityImport.Rho;
            Phydro[leaf]->Div = TempHydroDensityImport.Div;
            Phydro[leaf]->Rot[0] = TempHydroDensityImport.Rot[0];
            Phydro[leaf]->Rot[1] = TempHydroDensityImport.Rot[1];
            Phydro[leaf]->Rot[2] = TempHydroDensityImport.Rot[2];
#ifdef USE_FUVFEEDBACK //{
            Phydro[leaf]->GradRho = TempHydroDensityImport.GradRho;
#endif // USE_FUVFEEDBACK //}
        }
    }

    double TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.HydroDensityCommThisStep += GetElapsedTime()-TimeComm;

    CheckSizeofBufferHydroInteractionFlags(Pall.Nhydro);  
    bool (*HydroInteractionFlags)[NProcs-1];
    HydroInteractionFlags = (void *)BufferHydroInteractionFlags;
    for(int i=0;i<Pall.Nhydro;i++){
        for(int k=0;k<NProcs-1;k++){
            HydroInteractionFlags[i][k] = false;
        }
    }

    NImportAll = 0;
    for(int i=0;i<NProcs-1;i++){
        for(int k=0;k<NImportThisTime[i];k++){
            HydroDensityImportSend[NImportAll] = 
            //ReturnStructureDensityDivRotEngineWithGatherFlags(
                //HydroDensityExportRecv[NImportAll].Pos,HydroDensityExportRecv[NImportAll].Vel,
                    //HydroDensityExportRecv[NImportAll].Kernel,i,NProcs,HydroInteractionFlags);
            ReturnStructureDensityDivRotEngineDirect(
                HydroDensityExportRecv[NImportAll].Pos,HydroDensityExportRecv[NImportAll].Vel,
                    HydroDensityExportRecv[NImportAll].Kernel);
            HydroDensityImportSend[NImportAll].Leaf = HydroDensityExportRecv[NImportAll].Leaf;
            NImportAll ++;
        }
    }


    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];
    NImportAll = 0;
    int NImportAllNew = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportThisTimeNew[i] = 0;
        for(int k=0;k<NImportThisTime[i];k++){
            if(HydroDensityImportSend[NImportAll].Nlist > 0){
                HydroDensityImportSend[NImportAllNew] = HydroDensityImportSend[NImportAll];
                NImportThisTimeNew[i] ++;
                NImportAllNew ++;
            }
            NImportAll ++;
        }
    }

    for(int i=0;i<NProcs-1;i++){
        MPI_Sendrecv(NImportThisTimeNew+i,1,MPI_INT,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_PRECOMM,
            NExportThisTimeNew+i,1,MPI_INT,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_PRECOMM,
                MPI_COMM_WORLD,&mpi_status);
    }

    struct StructHydroDensityImport *HydroDensityImportRecv[NProcs];
    CheckSizeofBufferImportRecv(NExportThisTimeNew,sizeof(struct StructHydroDensityImport));
    for(int i=0;i<NProcs-1;i++)
        HydroDensityImportRecv[i] = BufferImportRecv[i];

    NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        MPI_Isend(HydroDensityImportSend+NImport,
            NImportThisTimeNew[i]*sizeof(struct StructHydroDensityImport),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_DENSITY_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Send+i);
        MPI_Irecv(HydroDensityImportRecv[i],
            NExportThisTimeNew[i]*sizeof(struct StructHydroDensityImport),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_DENSITY_IMPORT+i,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+i);
        NImport += NImportThisTimeNew[i];
    }

    TimeComm = GetElapsedTime();
    MPI_Waitall(NProcs-1,mpi_request_Export_Recv,mpi_status_Export_Recv);
    MPI_Waitall(NProcs-1,mpi_request_Export_Send,mpi_status_Export_Send);
    TimingResults.HydroDensityCommThisStep += GetElapsedTime()-TimeComm;

    for(int i=0;i<NProcs-1;i++){
        for(int k=0;k<NExportThisTimeNew[i];k++){
            int leaf = HydroDensityImportRecv[i][k].Leaf;
            Phydro[leaf]->Nlist += HydroDensityImportRecv[i][k].Nlist;
            Phydro[leaf]->Rho += HydroDensityImportRecv[i][k].Rho;
            Phydro[leaf]->Div += HydroDensityImportRecv[i][k].Div;
            Phydro[leaf]->Rot[0] += HydroDensityImportRecv[i][k].Rot[0];
            Phydro[leaf]->Rot[1] += HydroDensityImportRecv[i][k].Rot[1];
            Phydro[leaf]->Rot[2] += HydroDensityImportRecv[i][k].Rot[2];
#ifdef USE_FUVFEEDBACK //{
            Phydro[leaf]->GradRho += HydroDensityImportRecv[i][k].GradRho;
#endif // USE_FUVFEEDBACK //}
        }
    }
    for(int k=0;k<Pall.Nhydro;k++){
        if(Phydro[k]->Active == true){
            int leaf = k;
            double iRho = 1.e0/Phydro[leaf]->Rho;
            Phydro[leaf]->RhoPred = Phydro[leaf]->Rho;
            Phydro[leaf]->Div *= iRho;
            Phydro[leaf]->Rot[0] *= iRho;
            Phydro[leaf]->Rot[1] *= iRho;
            Phydro[leaf]->Rot[2] *= iRho;

#if defined(ISOTHERMAL_EOS_RUN)
            double cs = Pall.CS;
#elif defined(BAROTROPIC_EOS_RUN)
            double cs = ReturnSoundSpeedForBarotoropicRun(Phydro[leaf]->RhoPred);
#else
            double cs = sqrt(Pall.GGm1*Phydro[leaf]->UPred);
#endif
            double div = fabs(Phydro[leaf]->Div);
            Phydro[leaf]->F = (div/(div+NORM(Phydro[leaf]->Rot)+0.0001*cs*Phydro[leaf]->KernelPred));

#ifdef USE_VARIABLE_ALPHA
            double InvTau = Pall.ViscousL*cs/Phydro[leaf]->KernelPred;
            // Phydro[leaf]->DAlpha = -(Phydro[leaf]->Alpha-Pall.ViscousAlphaMin)*InvTau
                            // +Pall.ViscousS*Phydro[leaf]->F*fmax(0.e0,-Phydro[leaf]->Div);
            Phydro[leaf]->DAlpha = -(Phydro[leaf]->Alpha-Pall.ViscousAlphaMin)*InvTau
                            +fmax(0.e0,-Phydro[leaf]->Div*(Pall.ViscousAlphaMax-Phydro[leaf]->Alpha));
#endif // USE_VARIABLE_ALPHA
        }
    }

    return;
}
