#ifndef UNITPARM_H
#define UNITPARM_H

#include "Common.h"

namespace UnitParameter
{

const double Avogdro = 6.0221367e+23;     // avogdro's number (mol^-1)
const double Boltzmann = 8.617080363e-5;  // Boltzmann's number (eV K^-1)
const double CarbonWeight = 12.0107;      // atomic weight of C

const double unitMass = CarbonWeight/Avogdro*1e-3;           // unit mass (kg) 
const double unitLength = 1e-10;                             // unit length (m) 
const double unitEnergy_eV = 14.39965173;                    // unit energy (eV) 
const double unitEnergy_J  = unitEnergy_eV*1.60219e-19;      // unit energy (J) 
const double unitEnergy_kcal_mol = unitEnergy_eV*2.30492e+1; // unit energy (kcal/mol) 
const double unitTime = sqrt(unitMass / unitEnergy_J) * unitLength;    // unit time (s) 
const double kJmol2kcalmol = 0.2388886898;
const double eV2kcalmol = 2.30491e+1;
const double unitPressure_Pa = UnitParameter::unitEnergy_J/(UnitParameter::unitLength*UnitParameter::unitLength*UnitParameter::unitLength);

namespace {
inline double normalizeLength(const double length){ return length/unitLength;}
inline double normalizeVelocity(const double velocity){ return velocity/(unitLength/unitTime);}
inline double denormalizeVelocity(const double velocity){ return velocity*(unitLength/unitTime);}
inline double normalizeTime(const double time){ return time/unitTime;}
inline double normalizeEnergy_eV(const double energy){ return energy/unitEnergy_eV;}
inline double normalizeEnergy_kJmol(const double energy){ 
  return energy*kJmol2kcalmol/unitEnergy_kcal_mol;
}
inline double normalizeEnergy_kcal_mol(const double energy){ 
  return energy/2.30492e+1/unitEnergy_eV;
}
inline double normalizeMassAtomicWeight(const double mass){ return mass/CarbonWeight;}
inline double normalizeTemperature(const double temp){ return temp*Boltzmann/unitEnergy_eV;}
inline double realizeTemperature(const double temp){ return temp*unitEnergy_eV/Boltzmann;}
inline double realizeEnergy_kcal_mol(const double erg){ return erg*unitEnergy_kcal_mol;}
//inline double multiUnitMass(const double mass) {return mass*unitMass;}
}

}
#endif
