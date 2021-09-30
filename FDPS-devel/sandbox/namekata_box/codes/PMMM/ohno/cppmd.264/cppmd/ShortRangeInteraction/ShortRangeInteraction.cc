#include "ShortRangeInteraction.h"

LJMixparameterArray ShortRange::ljmixparameters;
long ShortRange::sizeofljmixparameters;
double ShortRange::alpha=0.8;
double ShortRange::LJCutoffEnergyCorrection;
double ShortRange::LJCutoffEnergyCorrectionV;
double ShortRange::TotalLJCutoffEnergyCorrection;
#ifdef CHECK_ENERGY
double ShortRange::short_coulomb;
double ShortRange::short_lj;
#endif
double ShortRange::r1co0, ShortRange::r1co1, ShortRange::r1co2;
double ShortRange::r6co0, ShortRange::r6co1, ShortRange::r6co2;
double ShortRange::r12co0, ShortRange::r12co1, ShortRange::r12co2;
double ShortRange::r3co1, ShortRange::r3co2;
double ShortRange::r8co1, ShortRange::r8co2;
double ShortRange::r14co1, ShortRange::r14co2;


void ShortRange::set_Shift_Function(double cutoff2)
{
  double cutoff;
  
  cutoff = sqrt(cutoff2);

  double cutoff3 = cutoff2*cutoff;
  double cutoff4 = cutoff2*cutoff2;
  double cutoff5 = cutoff2*cutoff3;
  double cutoff10 = cutoff5*cutoff5;

  r3co1  =  5.0/cutoff4;
  r3co2  =  4.0/cutoff5;
  r8co1  = 10.0/(cutoff4*cutoff5);
  r8co2  =  9.0/(cutoff10);
  r14co1 = 16.0/(cutoff10*cutoff5);
  r14co2 = 15.0/(cutoff10*cutoff4*cutoff2);

  r1co1  =  r3co1/3.0*1.0;
  r1co2  =  r3co2/4.0*1.0;
  r6co1  =  r8co1/3.0*6.0;
  r6co2  =  r8co2/4.0*6.0;
  r12co1 = r14co1/3.0*12.0;
  r12co2 = r14co2/4.0*12.0;

  r1co0  = 1.0/cutoff            + r1co1*(cutoff3)- r1co2*(cutoff4);
  r6co0  = 1.0/(cutoff2*cutoff4 )+ r6co1*(cutoff3)- r6co2*(cutoff4);
  r12co0 = 1.0/(cutoff2*cutoff10)+r12co1*(cutoff3)-r12co2*(cutoff4);

}

size_t setljmixparameters(LJMixparameterArray& ljmt)
{
  ShortRange::ljmixparameters = ljmt;
  return ljmt.size();
}


/// TODO : ljcec and number_of_lj_particle depend atomtype



double calcLJCutoffEnergyCorrection(double cutoff,int ati, int atj)
{
  ShortRange::LJCutoffEnergyCorrection = 8./9.*M_PI*(ShortRange::ljmixparameters[ati][atj].potFirst*pow(cutoff,-9)
					 - 3.*ShortRange::ljmixparameters[ati][atj].potSecond*pow(cutoff,-3));
  return ShortRange::LJCutoffEnergyCorrection;
}


void
setLJCutoffEnergyCorrectionParmeter(double ljcec, int num_lj)
{
  ShortRange::LJCutoffEnergyCorrectionV = ljcec*num_lj*num_lj;
}

double
calcTotalLJCutoffEnergyCorrection(double volume)
{
  ShortRange::TotalLJCutoffEnergyCorrection = ShortRange::LJCutoffEnergyCorrectionV/volume;
  return ShortRange::TotalLJCutoffEnergyCorrection;
}
