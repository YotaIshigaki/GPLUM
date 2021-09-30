#include "config.h"

double GetMeanMolecularWeight(double MetallicityWeight){
    return (Pall.MeanMolecularWeight);
}

//[dyne m^2/kg^2] = [cm^3/g/s^2]
double GetUnitGravitationalConstant(void){
    return (GRAVITY_CONSTANT_CGS*(Pall.UnitMass*SQ(Pall.UnitTime))/CUBE(Pall.UnitLength));
    //return (GRAVITY_CONSTANT_CGS/(CUBE(Pall.UnitLength)/(Pall.UnitMass*SQ(Pall.UnitTime))));
}

double GetUnitAu(void){
    return (AU_CGS/Pall.UnitLength);
}

double GetUnitKpc(void){
    return (KPC_CGS/Pall.UnitLength);
}

double GetUnitMpc(void){
    return (MPC_CGS/Pall.UnitLength);
}

double GetUnitVel(void){
    return (Pall.UnitTime/Pall.UnitLength);
}

double GetUnitLightVel(void){
    return (LIGHT_C_CGS/(Pall.UnitLength/Pall.UnitTime));
}

double GetUnitHubble(void){
    return (HUBBLE_CGS/(Pall.UnitLength/(1.e+5*Pall.UnitTime)));
}

double GetUnitEnergy(void){
    return (SQ(Pall.UnitTime)/(Pall.UnitMass*SQ(Pall.UnitLength)));
}

double GetUnitSpecificEnergy(void){
    return (SQ(Pall.UnitTime)/SQ(Pall.UnitLength));
}

double GetUnitCoolingRate(void){
    return (CUBE(Pall.UnitTime)*Pall.UnitLength/Pall.UnitMass);
}

double GetUnitEnergyPerUnitTime(void){
    return (CUBE(Pall.UnitTime)/(Pall.UnitMass*SQ(Pall.UnitLength)));
}

double GetUnitInnerEnergy(void){
    return (Pall.DegreeOfFreedom*(BOLTZMANN_CONSTANT_CGS/(GetMeanMolecularWeight(0.e0)*PROTON_MASS_CGS))*GetUnitSpecificEnergy());
}

double GetUnitConversionFactorTemperatureToInnerEnergy(void){
    return (GetUnitInnerEnergy());
}

double GetUnitConversionFactorInnerEnergyToTemperature(void){
    return (1.0/GetUnitInnerEnergy());
}

double GetUnitDensity(void){
    return (CUBE(Pall.UnitLength)/(Pall.UnitMass));
}

double GetUnitDensityCGS(void){
    return (Pall.UnitMass/CUBE(Pall.UnitLength));
}

double GetUnitNumberDensityCGS(void){
    return (1.0/(NEUTRON_MASS_CGS*GetUnitDensity()));
}

double GetUnitMassSolarMass(void){
    return (Pall.UnitMass/MSUN_CGS);
}
