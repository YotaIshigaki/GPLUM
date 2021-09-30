#include "config.h"

void UpdateCosmologicalParameters(void){

#ifdef COSMOLOGICAL_RUN //{
    Pall.Redshift = CalcCurrentTtoZ();
    UpdateAdaptiveSofteningFactor();
    Pall.OverDensityZ = CalcOverDensityZ();
    //Pall.RHOcZ = CalcVirializedDensityZ();
    Pall.RHOcZ = CalcCriticalDensityZ();
    Pall.HubbleZ = CalcHubbleParameterZ();
#endif // COSMOLOGICAL_RUN //}
    return ;
}

void UpdateAdaptiveSofteningFactor(void){

    double fact = (1.e0+Pall.FrozenRedshift);

    if(Pall.Redshift < Pall.FrozenRedshift){
        Pall.AdaptiveSofteningFactor = 1.0/(1.e0+Pall.FrozenRedshift);
    } else {
        Pall.AdaptiveSofteningFactor = 1.0/(1.e0+Pall.Redshift);
    }
    Pall.AdaptiveSofteningFactor *= fact;

    return;
}

void OldUpdateAdaptiveSofteningFactor(void){

    if(Pall.Redshift < Pall.FrozenRedshift){
        Pall.AdaptiveSofteningFactor = 1.e0;
    } else {
        Pall.AdaptiveSofteningFactor = (Pall.FrozenRedshift+1.0)/(1.e0+Pall.Redshift);
    }

    return;
}

double CalcCriticalDensityZ0(void){
    const static double Factor = 3.0/(8.0*M_PI);
    return (Factor*SQ(CalcHubbleParameterZ0())/Pall.GravConst);
}

double CalcHubbleParameterZ0(void){
    //return (Pall.Hubble*sqrt(Pall.OmegaM*CUBE(1.e0+Pall.Redshift) + Pall.OmegaL));
    return (Pall.Hubble*sqrt(Pall.OmegaM*CUBE(1.e0+0.0) + Pall.OmegaL));
}

double CalcCriticalDensityZ(void){
    const static double Factor = 3.0/(8.0*M_PI);
    return (Factor*SQ(CalcHubbleParameterZ())/Pall.GravConst);
}

double CalcHubbleParameterZ(void){
    return (Pall.Hubble*sqrt(Pall.OmegaM*CUBE(1.e0+Pall.Redshift) + Pall.OmegaL));
}

double CalcOverDensityZ(void){

    /* References                            */
    /* Bryan and Norman ApJ 495 1998,        */
    /* Somerville and Primack MNRAS 310 1087.*/

    double OD,x,Oz,Ez;
    Ez = Pall.OmegaM*CUBE(1.e0+Pall.Redshift) + Pall.OmegaL;
    Oz = Pall.OmegaM*CUBE(1.e0+Pall.Redshift)/Ez;
    x = Oz - 1.e0;
    OD = 18.e0*SQ(PI) + 82.e0*x - 39.e0*SQ(x);
    return (OD);
}

double CalcVirializedDensityZ(void){ // OverDensity*CriticalDensity

    /* References                            */
    /* Bryan and Norman ApJ 495 1998,        */
    /* Somerville and Primack MNRAS 310 1087.*/

    double OD,RHOc;
    OD = CalcOverDensityZ();
    RHOc = CalcCriticalDensityZ();
    return (OD*RHOc);
}

void CalcUniformSphereForce(void){

    double fact = Pall.OmegaL*SQ(Pall.Hubble);
    //fprintf(stderr,"fact %g ,OL %g, H %g\n",fact,Pall.OmegaL,Pall.Hubble);

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            double crit = CalcCriticalDensityZ();
            double R = NORM(Pbody[i]->PosP);
            double acc = -Pall.GravConst*(4.0*M_PI/3.0)*crit*R;
            //double acc = -Pall.GravConst*(4.0*M_PI/3.0)*crit*(SQ(R)+SQ(Pbody[i]->Eps));
            double z = Pbody[i]->PosP[2]; 
            double r = sqrt(SQ(Pbody[i]->PosP[0])+SQ(Pbody[i]->PosP[1])); 
            double sintheta = sqrt(1-SQ(z/R));

            if(r>TINY){
                Pbody[i]->Acc[0] = acc*sintheta*Pbody[i]->PosP[0]/r;
                Pbody[i]->Acc[1] = acc*sintheta*Pbody[i]->PosP[1]/r;
                Pbody[i]->Acc[2] = acc*z/R;
            }
            if(Pall.OmegaL > TINY){
                Pbody[i]->Acc[0] += fact*Pbody[i]->PosP[0];
                Pbody[i]->Acc[1] += fact*Pbody[i]->PosP[1];
                Pbody[i]->Acc[2] += fact*Pbody[i]->PosP[2];

            }
        }
    }

    return ;
}

