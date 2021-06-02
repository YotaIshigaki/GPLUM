/* HydroMisc.c */

struct StructHydroDensityExport{
    double    Kernel;  // Kernel size.
    double    Pos[3];  // Position.
    double    Vel[3];  // Volocity.
    int       Leaf;
#if VISCOSITY_TYPE == 1 //{
    double    B[3][3];
#endif // VISCOSITY_TYPE == 1 //}
#ifdef HYDRO_TIMESTEP_LIMITER
    int       k_hydro;
#endif // HYDRO_TIMESTEP_LIMITER
#ifdef USE_DEBUG_MODE
     unsigned long int GlobalID;
#endif // USE_DEBUG_MODE
#ifdef USE_SMOOTHED_POTENTIAL //{
     double   Pot;
#endif // USE_SMOOTHED_POTENTIAL //}
};

struct StructHydroDensityImport{
    double    Rho;     // Density.
    double    Div;     // Divergense of velocity.
    double    Rot[3];  // Rotation of velocity.
#ifdef USE_DISPH //{
    double    EnergyDensity; // Energy Density
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
    double    PseudoDensity; // Pseudo Density
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
    double    Gradh;   // Gradh
#ifdef USE_GRAD_N //{
    double NumberDensity;
    double GradN;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
#if DIFFUSION_TYPE==0 //{
    double DZdiff;
#elif DIFFUSION_TYPE==1 //}//{
    double Sxy;
    double Sxz;
    double Syx;
    double Syz;
    double Szx;
    double Szy;
#endif //}
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
#ifdef USE_FUVFEEDBACK //{
    double GradRho;
#endif // USE_FUVFEEDBACK //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    double   Pot;
#endif // USE_SMOOTHED_POTENTIAL //}
    int       Nlist;   // Nlist.
    int       Leaf;
};

struct StructHydroAccExport{
    double    Mass;    // Mass.
    double    Kernel;  // Kernel size.
    double    Pos[3];  // Position.
    double    Vel[3];  // Volocity.
    double    U;       // Thermal energy.
    double    Rho;     // Density.
    double    F;       // div/(div+rot).
#ifdef USE_DISPH //{
    double    EnergyDensity; // Energy Density
#endif // USE_DISPH //}
#ifdef USE_SPSPH //{
    double    PseudoDensity; // Pseudo-Density
    double    Zw;            // Weight for SPSPH
    double    Ddif;          // Diffusion coef.
#endif // USE_SPSPH //}
#ifdef USE_GRAD_H //{
    double    Gradh;   // Gradh
#ifdef USE_GRAD_N //{
    double    fij;
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
#ifdef USE_VARIABLE_ALPHA
    double    Alpha;   // Alpha parameter for artificial viscosity.
#endif // USE_VARIABLE_ALPHA
#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
    bool Active;
    double dt;
#endif //(defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}
#ifdef USE_METAL_DIFFUSION //{
#ifdef USE_CELIB //{
#if DIFFUSION_TYPE==0 //{
    double DZdiff;
#elif DIFFUSION_TYPE==1 //}/{
    double Sxy;
    double Sxz;
    double Syx;
    double Syz;
    double Szx;
    double Szy;
#endif  // DIFFUSION_TYPE //}
    double Elements[CELibYield_Number-1]; //He, C, N, O, Ne, Mg, Si, S, Ca, Fe, Ni, (Eu)
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
#if VISCOSITY_TYPE == 1 //{
    double Div;
    double Rot[3];
#endif // VISCOSITY_TYPE == 1 //{
#ifdef USE_MOMENTUM_FEEDBACK //{
    double Z;
#endif // USE_MOMENTUM_FEEDBACK //}

    int       Leaf;
};

struct StructHydroAccImport{
    double    Du;     // du/dt.
    double    HydroAcc[3];  // dp/dr.
    double    Vsig;   // Velocity dispresion in the local scale.
#ifdef USE_SPSPH //{
    double    DZw;
#endif // USE_SPSPH //}
#if (defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //{
    bool Active;
    double dt;
#endif //(defined(USE_SPSPH) || defined(USE_THERMAL_CONDUCTIVITY) || defined(USE_METAL_DIFFUSION)) //}
#ifdef USE_CELIB //{
#ifdef USE_METAL_DIFFUSION //{
    double dZ[CELibYield_Number-1]; //He, C, N, O, Ne, Mg, Si, S, Ca, Fe, Ni, (Eu)
#endif // USE_METAL_DIFFUSION //}
#endif // USE_CELIB //}
    int       Nlist;
    int       Leaf;
};


void ClearHydroData(void);
void InitializeHydroInteractionFlags(void);
void CheckSizeofBufferHydroInteractionFlags(const int NumberofBufferHydroInteractionFlags);
void ReleaseBufferHydroInteractionFlags(void);
void WriteKernelShape(void);


