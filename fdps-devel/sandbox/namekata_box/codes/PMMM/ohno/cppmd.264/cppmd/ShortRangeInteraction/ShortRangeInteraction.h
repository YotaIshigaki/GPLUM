#ifndef SHORTRANGEINTERACTION_H
#define SHORTRANGEINTERACTION_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H

#include <iostream>

#include "LJParameter.h"


#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
#include "EwaldRealInterpolation.h"
#define EWALD_REAL_FORCE(x) ewaldRealInterpolation.ewaldRealForceTbl(x)
#define EWALD_REAL_POT(x) ewaldRealInterpolation.ewaldRealPotTbl(x)
#else  // CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
#define EWALD_REAL_FORCE(x) (M_2_SQRTPI*x*exp(-x*x) + erfc(x))
#define EWALD_REAL_POT(x) erfc(x)
#endif  // CPPMD_ENABLE_TABLE_FOR_EWALD_REAL

#include "ParticleInfo.h"
#include "SetPairList.h"

//! common data and indentifier for ShortRange
namespace ShortRange {
  extern LJMixparameterArray ljmixparameters;
  extern long sizeofljmixparameters;
  extern double alpha;     // alpha for Ewald
  // TODO : It may be depend pair of atomtype, must be array
  extern double LJCutoffEnergyCorrection;
  extern double LJCutoffEnergyCorrectionV; //! LJ energey cutoff correction per  inv volume
  extern double TotalLJCutoffEnergyCorrection;
#ifdef CHECK_ENERGY
  extern double short_coulomb;
  extern double short_lj;
#endif
  enum CoulombType{
    NoCoulomb,
    OriginalCoulomb,
    ForEwaldReal,
    IPSCoulomb,
    MWolf,
    ZeroDipole
  };
  enum ActiveInteraction{
    ActiveLJ,
    ActiveCoulomb
  };
  enum InteractionType{
    Coulomb,
    EwaldReal,
    LJ,
    LJCoulomb,
    LJEwaldReal
  };
//! coefficient for GROMACS Shift Function
  extern double r1co0, r1co1, r1co2;
  extern double r6co0, r6co1, r6co2;
  extern double r12co0, r12co1, r12co2;
  extern double r3co1, r3co2;
  extern double r8co1, r8co2;
  extern double r14co1, r14co2;

void set_Shift_Function(double cutoff2);

}

size_t setljmixparameters(LJMixparameterArray& ljmt);

double calcLJCutoffEnergyCorrection(double cutoff,int ati, int atj);

void
setLJCutoffEnergyCorrectionParmeter(double ljcec, int num_lj);

double
calcTotalLJCutoffEnergyCorrection(double volume);

template<ShortRange::InteractionType T, typename  TPI, typename TPJ>
  inline void Interaction(double r2, TPI pi, TPJ pj,
                          double& energy, double& dp)
{
}

template<> inline void Interaction<ShortRange::Coulomb>(double r2, Particle pi, Particle pj,
 double& energy, double& dp)
{
  double _r = 1.0/sqrt(r2);
  double _r3 = _r*_r*_r;
#ifdef CHECK_ENERGY
  double sp = pi.charge*pj.charge*_r;
  ShortRange::short_coulomb += sp;
  energy += sp;
#else
  energy    += pi.charge*pj.charge*_r;
#endif
  dp  = pi.charge*pj.charge*_r3;
}

template<> inline void Interaction<ShortRange::EwaldReal>(double r2, Particle pi, Particle pj,
 double& energy, double& dp)
{
  double r = sqrt(r2);
  double _r = 1.0/r;
  double _r3 = _r*_r*_r;
#ifdef CHECK_ENERGY
  double sp = pi.charge*pj.charge*_r*EWALD_REAL_POT(ShortRange::alpha*r);
  ShortRange::short_coulomb += sp;
  energy += sp;
#else
  energy += pi.charge*pj.charge*_r*EWALD_REAL_POT(ShortRange::alpha*r);
#endif
  dp      = pi.charge*pj.charge*_r3*EWALD_REAL_FORCE(ShortRange::alpha*r);
}

template<> inline void Interaction<ShortRange::LJ>(double r2, Particle pi, Particle pj,
 double& energy, double& dp)
{
  int atomi = pi.atomtype;
  int atomj = pj.atomtype;
  double _r2 = 1.0/r2;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
#ifdef CHECK_ENERGY
  double sp = ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
    - ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6;
  ShortRange::short_lj += sp;
  energy += sp;
#else
  energy   += ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
    - ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6;
#endif
  dp = ShortRange::ljmixparameters[atomi][atomj].forceFirst*_r14 
    - ShortRange::ljmixparameters[atomi][atomj].forceSecond*_r8;
}

template<> inline void Interaction<ShortRange::LJCoulomb>(double r2, Particle pi, Particle pj,
 double& energy, double& dp)
{
  int atomi = pi.atomtype;
  int atomj = pj.atomtype;
  double _r = 1.0/sqrt(r2);
  double _r2 = _r*_r;
  double _r3 = _r2*_r;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
#ifdef CHECK_ENERGY
  double scp = pi.charge*pj.charge*_r;
  double slj = ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
    - ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6;
  ShortRange::short_coulomb += scp;
  ShortRange::short_lj += slj;
  energy += slj + scp;
#else
  energy   += ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
    - ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6
    + pi.charge*pj.charge*_r;
#endif
  dp = ShortRange::ljmixparameters[atomi][atomj].forceFirst*_r14 
    - ShortRange::ljmixparameters[atomi][atomj].forceSecond*_r8
    + pi.charge*pj.charge*_r3;
}

template<> inline void Interaction<ShortRange::LJEwaldReal>(double r2, Particle pi, Particle pj,
 double& energy, double& dp)
{
  int atomi = pi.atomtype;
  int atomj = pj.atomtype;
  double  r = sqrt(r2);
  double _r = 1.0/r;
  double _r2 = _r*_r;
  double _r6 = _r2*_r2*_r2;
  double _r12 = _r6*_r6;
#ifdef CHECK_ENERGY
  double scp = pi.charge*pj.charge*_r*EWALD_REAL_POT(ShortRange::alpha*r);
  double slj = ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
    - ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6;
  ShortRange::short_coulomb += scp;
  ShortRange::short_lj += slj;
  energy += slj + scp;
#else
  energy   +=    ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
    -    ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6
    +    pi.charge*pj.charge*_r*EWALD_REAL_POT(ShortRange::alpha*r);
#endif
  dp = (  ShortRange::ljmixparameters[atomi][atomj].forceFirst*_r12 
          - ShortRange::ljmixparameters[atomi][atomj].forceSecond*_r6
          + pi.charge*pj.charge*_r*EWALD_REAL_FORCE(ShortRange::alpha*r))*_r2;
}

inline void Interaction(const double r2, 
                        const double chargei,
                        const double chargej,
                        const Atomtype atomi, 
                        const Atomtype atomj,
                        double& energy, double& dp)
{
  using namespace ShortRange;
  double _r = 1.0/sqrt(r2);
  double _r2 = _r*_r;
  double _r3 = _r2*_r;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
  LJMixParameter &lj = ShortRange::ljmixparameters[atomi][atomj];
#ifdef CHECK_ENERGY
  double scp = chargei*chargej*(_r);
  double slj =  lj.potFirst*(_r12) - lj.potSecond*(_r6);
  short_coulomb += scp;
  short_lj += slj;
  energy += slj + scp;
#else
  energy   += lj.potFirst*(_r12)
      - lj.potSecond*(_r6)
      + chargei*chargej*(_r);
#endif
  dp = lj.forceFirst*(_r14) 
    - lj.forceSecond*(_r8)
    + chargei*chargej*(_r3);
}
inline void Interaction(const double r2, 
                        const double chargei,
                        const double chargej,
                        const Atomtype atomi, 
                        const Atomtype atomj,
                        double& ljenergy, double& ljdp,
                        double& energy, double& dp)
{
  using namespace ShortRange;
  double _r = 1.0/sqrt(r2);
  double _r2 = _r*_r;
  double _r3 = _r2*_r;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
  LJMixParameter &lj = ShortRange::ljmixparameters[atomi][atomj];
#ifdef CHECK_ENERGY
  double scp = chargei*chargej*(_r);
  double slj =  lj.potFirst*(_r12) - lj.potSecond*(_r6);
  short_coulomb += scp;
  short_lj += slj;
  ljenergy += slj;
  energy += scp;
#else
  ljenergy   += lj.potFirst*(_r12)
      - lj.potSecond*(_r6);
  energy +=  chargei*chargej*(_r);
#endif
  ljdp = lj.forceFirst*(_r14) 
      - lj.forceSecond*(_r8);
  dp = chargei*chargej*(_r3);
}
//! See GROMACS Manual
inline void Interaction_LJShiftCoulombShift(const double r2, 
                                            const double chargei,
                                            const double chargej,
                                            const Atomtype atomi, 
                                            const Atomtype atomj,
                                            double& energy, double& dp)
{
  using namespace ShortRange;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double r3 = r2*r;
  double _r2 = _r*_r;
  double _r3 = _r2*_r;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
  LJMixParameter &lj = ShortRange::ljmixparameters[atomi][atomj];
#ifdef CHECK_ENERGY
  double scp = chargei*chargej*(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ));
  double slj = lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
  short_coulomb += scp;
  short_lj += slj;
  energy += scp + slj;
#else
  energy   += lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
      + chargei*chargej*(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ));
#endif
  dp = lj.forceFirst*(_r14-r*(r14co1-r14co2*r)) 
    - lj.forceSecond*(_r8 -r*(r8co1 -r8co2 *r))
    + chargei*chargej*(_r3 -r*(r3co1 -r3co2 *r));
}
inline void Interaction_LJShiftCoulombShift(const double r2, 
                                            const double chargei,
                                            const double chargej,
                                            const Atomtype atomi, 
                                            const Atomtype atomj,
                                            double& ljenergy, double& ljdp,
                                            double& energy, double& dp)
{
  using namespace ShortRange;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double r3 = r2*r;
  double _r2 = _r*_r;
  double _r3 = _r2*_r;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
  LJMixParameter &lj = ShortRange::ljmixparameters[atomi][atomj];
#ifdef CHECK_ENERGY
  double scp = chargei*chargej*(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ));
  double slj = lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
  short_coulomb += scp;
  short_lj += slj;
  ljenergy += slj;
  energy += scp;
#else
  ljenergy   += lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
  energy +=  chargei*chargej*(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ));
#endif
  ljdp = lj.forceFirst*(_r14-r*(r14co1-r14co2*r)) 
      - lj.forceSecond*(_r8 -r*(r8co1 -r8co2 *r));
  dp = chargei*chargej*(_r3 -r*(r3co1 -r3co2 *r));
}
inline void Interaction_LJShiftCoulomb(const double r2, 
				       const double chargei,
				       const double chargej,
				       const Atomtype atomi, 
				       const Atomtype atomj,
				       double& energy, double& dp)
{
  using namespace ShortRange;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double r3 = r2*r;
  double _r2 = _r*_r;
  double _r3 = _r2*_r;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
  LJMixParameter &lj = ShortRange::ljmixparameters[atomi][atomj];
#ifdef CHECK_ENERGY
  double scp = chargei*chargej*(_r);
  double slj = lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
  short_coulomb += scp;
  short_lj += slj;
  energy += scp + slj;
#else
  energy   += lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
      + chargei*chargej*(_r);
#endif
  dp = lj.forceFirst*(_r14-r*(r14co1-r14co2*r)) 
    - lj.forceSecond*(_r8 -r*(r8co1 -r8co2 *r))
    + chargei*chargej*(_r3);
}
inline void Interaction_LJShiftCoulomb(const double r2, 
				       const double chargei,
				       const double chargej,
				       const Atomtype atomi, 
				       const Atomtype atomj,
				       double& ljenergy, double& ljdp,
				       double& energy, double& dp)
{
  using namespace ShortRange;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double r3 = r2*r;
  double _r2 = _r*_r;
  double _r3 = _r2*_r;
  double _r6 = _r2*_r2*_r2;
  double _r8 = _r6*_r2;
  double _r12 = _r6*_r6;
  double _r14 = _r8*_r6;
  LJMixParameter &lj = ShortRange::ljmixparameters[atomi][atomj];
#ifdef CHECK_ENERGY
  double scp = chargei*chargej*(_r);
  double slj = lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
  short_coulomb += scp;
  short_lj += slj;
  ljenergy += slj;
  energy += scp;
#else
  ljenergy   += lj.potFirst*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - lj.potSecond*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
  energy +=  chargei*chargej*(_r);
#endif
  ljdp = lj.forceFirst*(_r14-r*(r14co1-r14co2*r)) 
      - lj.forceSecond*(_r8 -r*(r8co1 -r8co2 *r));
  dp = chargei*chargej*(_r3);
}


//! coulomb force and potential
class Coulomb {
public:
  template<typename TPI, typename TPJ>
  inline
  void Interaction(double r2, TPI pi, TPJ pj,
                   double& energy, double& dp)
  {
    double _r = 1.0/sqrt(r2);
    double _r3 = _r*_r*_r;
    energy    += pi.charge*pj.charge*_r;
    dp  = pi.charge*pj.charge*_r3;

  }
};
//! Ewald Real part force and potential 
class EwaldReal {
public:
  template<typename TPI, typename TPJ>
  inline
  void Interaction(double r2, TPI pi, TPJ pj,
                   double& energy, double& dp)
  {
    double r = sqrt(r2);
    double _r = 1.0/r;
    double _r3 = _r*_r*_r;
    energy += pi.charge*pj.charge*_r*EWALD_REAL_POT(ShortRange::alpha*r);
    dp      = pi.charge*pj.charge*_r3*EWALD_REAL_FORCE(ShortRange::alpha*r);
  }
};

//! Lennard Jones force and potential
class LJ {
public:
  template<typename TPI, typename TPJ>
  inline
  void Interaction(double r2, TPI pi, TPJ pj,
                   double& energy, double& dp)
  {
    int atomi = pi.atomtype;
    int atomj = pj.atomtype;
    double _r2 = 1.0/r2;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
    energy   += ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
      - ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6;
    dp = ShortRange::ljmixparameters[atomi][atomj].forceFirst*_r14 
      - ShortRange::ljmixparameters[atomi][atomj].forceSecond*_r8;
  }
};
//! Lennard Jones and Coulomb force and potential
class LJCoulomb {
public:
  template<typename TPI, typename TPJ>
  inline
  void Interaction(double r2, TPI pi, TPJ pj,
                   double& energy, double& dp)
  {
    int atomi = pi.atomtype;
    int atomj = pj.atomtype;
    double _r = 1.0/sqrt(r2);
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
    energy   += ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
      - ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6
      + pi.charge*pj.charge*_r;
    dp = ShortRange::ljmixparameters[atomi][atomj].forceFirst*_r14 
      - ShortRange::ljmixparameters[atomi][atomj].forceSecond*_r8
      + pi.charge*pj.charge*_r3;
  }
};
//! Lennard Jones and Ewald Real part force and potential
class LJEwaldReal {
public:
  template<typename TPI, typename TPJ>
  inline
  void Interaction(double r2, TPI pi, TPJ pj,
                   double& energy, double& dp)
  {
    int atomi = pi.atomtype;
    int atomj = pj.atomtype;
    double  r = sqrt(r2);
    double _r = 1.0/r;
    double _r2 = _r*_r;
    double _r6 = _r2*_r2*_r2;
    double _r12 = _r6*_r6;
    energy   +=    ShortRange::ljmixparameters[atomi][atomj].potFirst*_r12
      -    ShortRange::ljmixparameters[atomi][atomj].potSecond*_r6
      +    pi.charge*pj.charge*_r*EWALD_REAL_POT(ShortRange::alpha*r);
    dp = (  ShortRange::ljmixparameters[atomi][atomj].forceFirst*_r12 
            - ShortRange::ljmixparameters[atomi][atomj].forceSecond*_r6
            + pi.charge*pj.charge*_r*EWALD_REAL_FORCE(ShortRange::alpha*r))*_r2;
  }
};

//! Short range interaction loops (ver. PairRangeList)
template<class T, typename TPI = ParticleParameter1, typename TPJ=ParticleParameter1>
class ShortRangeInteraction {
public:
  T interaction;

  ShortRangeInteraction() : interaction() {}

  void loopiself(std::vector<TPI>& particlei, 
                 double &energy) {
    for(long i=0;i<particlei.size()-1;i++){
      Force forcei(0.0,0.0,0.0);
      double energyi=0.0;
      for(long j=i+1;j<particlei.size();j++){
        double energyij = 0.0;
        double dp;
        Position d = particlei[i].position - particlei[j].position;
        double r2 = d.norm2();
        interaction.Interaction(r2,particlei[i],particlei[j],energyij,dp);
        Force forceij = d*dp;
        energyi += energyij;
        particlei[j].force -= forceij;
        forcei += forceij;
      }
      particlei[i].force += forcei;
      energy += energyi;
    }
  }
  void loopiself(std::vector<TPI>& particlei, 
                 double &energy,
                 double cutoff2) {
    for(long i=0;i<particlei.size()-1;i++){
      Force forcei(0.0,0.0,0.0);
      double energyi=0.0;
      for(long j=i+1;j<particlei.size();j++){
        double energyij = 0.0;
        double dp;
        Position d = particlei[i].position - particlei[j].position;
        double r2 = d.norm2();
        if(r2<cutoff2){
          interaction.Interaction(r2,particlei[i],particlei[j],energyij,dp);
          Force forceij = d*dp;
          energyi += energyij;
          particlei[j].force -= forceij;
          forcei += forceij;
        }
      }
      particlei[i].force += forcei;
      energy += energyi;
    }
  }

  void loopij(std::vector<TPI>& particlei, 
              std::vector<TPJ>& particlej,
              double &energy) {
    for(long i=0;i<particlei.size();i++){
      Force forcei(0.0,0.0,0.0);
      double energyi=0.0;
      for(long j=0;j<particlej.size();j++){
        double energyij = 0.0;
        double dp;
        Position d = particlei[i].position - particlej[j].position;
        double r2 = d.norm2();
        interaction.Interaction(r2,particlei[i],particlej[j],energyij,dp);
        Force forceij = d*dp;
        energyi += energyij;
        particlej[j].force -= forceij;
        forcei += forceij;
      }
      particlei[i].force += forcei;
      energy += energyi;
    }
  }
  void loopij(std::vector<TPI>& particlei, 
              std::vector<TPJ>& particlej,
              double &energy,
              double cutoff2) {
    for(long i=0;i<particlei.size();i++){
      Force forcei(0.0,0.0,0.0);
      double energyi=0.0;
      for(long j=0;j<particlej.size();j++){
        double energyij = 0.0;
        double dp;
        Position d = particlei[i].position - particlej[j].position;
        double r2 = d.norm2();
        if(r2<cutoff2){
          interaction.Interaction(r2,particlei[i],particlej[j],energyij,dp);
          Force forceij = d*dp;
          energyi += energyij;
          particlej[j].force -= forceij;
          forcei += forceij;
        }
      }
      particlei[i].force += forcei;
      energy += energyi;
    }
  }
  
  void loopjset(const PairRangeList& jsets, 
                 std::vector<TPI>& particlei, 
                 std::vector<TPJ>& particlej,
                 double &energy) {
    for(long i=0;i<particlei.size();i++){
      Force forcei(0.0,0.0,0.0);
      double energyi=0.0;
      for(long jset=0;jset<jsets.size();jset++){
        for(long j=jsets[jset].begin;j<jsets[jset].end;j++){
          double energyij = 0.0;
          double dp;
          Position d = particlei[i].position - particlej[j].position;
          double r2 = d.norm2();
          interaction.Interaction(r2,particlei[i],particlej[j],energyij,dp);
          Force forceij = d*dp;
          energyi += energyij;
          particlej[j].force -= forceij;
          forcei += forceij;
        }
      }
      particlei[i].force += forcei;
      energy += energyi;
    }
  }

  void loopjset(const PairRangeList& jsets, 
                std::vector<TPI>& particlei, 
                std::vector<TPJ>& particlej,
                double &energy,
                double cutoff2) {
    for(long i=0;i<particlei.size();i++){
      Force forcei(0.0,0.0,0.0);
      double energyi=0.0;
      for(long jset=0;jset<jsets.size();jset++){
        for(long j=jsets[jset].begin;j<jsets[jset].end;j++){
          double energyij = 0.0;
          double dp;
          Position d = particlei[i].position - particlej[j].position;
          double r2 = d.norm2();
          if(r2<cutoff2){
            interaction.Interaction(r2,particlei[i],particlej[j],energyij,dp);
            Force forceij = d*dp;
            energyi += energyij;
            particlej[j].force -= forceij;
            forcei += forceij;
          }
        }
      }
      particlei[i].force += forcei;
      energy += energyi;
    }
  }

  void loopijsets(const PairRangeList& isets,
                  const std::vector<PairRangeList>& jsets, 
                  std::vector<TPI>& particlei, 
                  std::vector<TPJ>& particlej,
                  double &energy) {
    for(long iset=0;iset<isets.size();iset++){
      for(long i=isets[iset].begin;i<isets[iset].end;i++){
        Force forcei(0.0,0.0,0.0);
        double energyi=0.0;
        for(long jset=0;jset<jsets[iset].size();jset++){
          for(long j=jsets[iset][jset].begin;j<jsets[iset][jset].end;j++){
            double energyij = 0.0;
            double dp;
            Position d = particlei[i].position - particlej[j].position;
            double r2 = d.norm2();
            interaction.Interaction(r2,particlei[i],particlej[j],energyij,dp);
            Force forceij = d*dp;
            energyi += energyij;
            particlej[j].force -= forceij;
            forcei += forceij;
          }
        }
        particlei[i].force += forcei;
        energy += energyi;
      }
    }
  }

  void loopijsets(const PairRangeList& isets,
                  const std::vector<PairRangeList>& jsets, 
                  std::vector<TPI>& particlei, 
                  std::vector<TPJ>& particlej,
                  double &energy,
                  double cutoff2) {
    for(long iset=0;iset<isets.size();iset++){
      for(long i=isets[iset].begin;i<isets[iset].end;i++){
        Force forcei(0.0,0.0,0.0);
        double energyi=0.0;
        for(long jset=0;jset<jsets[iset].size();jset++){
          for(long j=jsets[iset][jset].begin;j<jsets[iset][jset].end;j++){
            double energyij = 0.0;
            double dp;
            Position d = particlei[i].position - particlej[j].position;
            double r2 = d.norm2();
            if(r2<cutoff2){
              interaction.Interaction(r2,particlei[i],particlej[j],energyij,dp);
              Force forceij = d*dp;
              energyi += energyij;
              particlej[j].force -= forceij;
              forcei += forceij;
            }
          }
        }
        particlei[i].force += forcei;
        energy += energyi;
      }
    }
  }
};

//! Short range interaction loops for multiple type (ver. PairRangeList)
template<class TCL, class TLJ, 
         typename TPI, typename TPJ>
class ShortRangeInteractions {
public:
  ShortRangeInteraction<LJ> lj;
  ShortRangeInteraction<TCL,TPI,TPJ> coulomb;
  ShortRangeInteraction<TLJ,TPI,TPJ> ljc;
  
  ShortRangeInteractions() : coulomb(), lj() {}
  
  void loopijsets(const PairRangeList& ljcisets,
                  const std::vector<PairRangeList>& ljcjsets,
                  const PairRangeList& coulombisets,
                  const std::vector<PairRangeList>& coulombjsets, 
                  std::vector<TPI>& particlei, 
                  std::vector<TPJ>& particlej,
                  double &energy)
  {
    ljc.loopijsets(ljcisets,ljcjsets,particlei,particlej,energy);
    coulomb.loopijsets(coulombisets,ljcjsets,particlei,particlej,energy);
    coulomb.loopijsets(ljcisets,coulombjsets,particlei,particlej,energy);
    coulomb.loopijsets(coulombisets,coulombjsets,particlei,particlej,energy);
  }

  void loopijsets(PairRangeList& ljcisets,
                  std::vector<PairRangeList>& ljcjsets,
                  PairRangeList& coulombisets,
                  std::vector<PairRangeList>& coulombjsets, 
                  std::vector<TPI>& particlei, 
                  std::vector<TPJ>& particlej,
                  double &energy,
                  double cutoff2)
  {
    ljc.loopijsets(ljcisets,ljcjsets,particlei,particlej,energy,cutoff2);
    coulomb.loopijsets(coulombisets,ljcjsets,particlei,particlej,energy,cutoff2);
    coulomb.loopijsets(ljcisets,coulombjsets,particlei,particlej,energy,cutoff2);
    coulomb.loopijsets(coulombisets,coulombjsets,particlei,particlej,energy,cutoff2);
  }
  
};

//! Short range interaction loops for LJ, LJEwaldReal and EwaldReal with ParticleArray
typedef ShortRangeInteractions<EwaldReal, LJEwaldReal, Particle, Particle> EwaldAndLJEInteractions;

//! Short range interaction loops for LJ, LJCoulomb and Coulomb with Particlearray
typedef ShortRangeInteractions<Coulomb, LJ, Particle, Particle> CoulombAndLJInteractions;

//! Short range interaction loops for LJ, LJEwaldReal and EwaldReal with ParticleParameters1
typedef ShortRangeInteractions<EwaldReal, LJEwaldReal,ParticleParameter1,ParticleParameter1> EwaldAndLJEInteractions1;

//! Short range interaction loops for LJ, LJCoulomb and Coulomb with ParticleParameters1
typedef ShortRangeInteractions<Coulomb, LJ,ParticleParameter1,ParticleParameter1> CoulombAndLJInteractions1;

template<ShortRange::CoulombType T>
inline double Potential_Form_Full(double r2, double cutoff2)
{
  double _r = 1.0/sqrt(r2);
  return _r;
}

template<ShortRange::CoulombType T>
inline double Force_Form_Full(double r2, double cutoff2)
{
  double _r = 1.0/sqrt(r2);
  return _r*_r*_r;
}

template<ShortRange::CoulombType T>
inline double ForcePotential_Form_Full(double r2, double cutoff2, double& e)
{
  double _r = 1.0/sqrt(r2);
  e = _r;
  return _r*_r*_r;
}

template<> inline double Potential_Form_Full<ShortRange::ForEwaldReal>(double r2, double cutoff2)
{
  double r = sqrt(r2);
  double _r = 1.0/r;
  double ep;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
  ep = EWALD_REAL_POT(ShortRange::alpha*r);
#else
  ep = erfc(ShortRange::alpha*r);
#endif
  double correct_ewaldpot = _r*ep;
#ifdef DUMP_ENERGY_CONTENT
  printf("correct ewald %e : r2 = %e, r = %e, alpha = %e, alpha*r = %e, ewald pot %e, erfc(ar) %e\n",correct_ewaldpot,r2,r,ShortRange::alpha,ShortRange::alpha*r,ep,erfc(ShortRange::alpha*r));
#endif
  return correct_ewaldpot;
}

template<> inline double Force_Form_Full<ShortRange::ForEwaldReal>(double r2, double cutoff2)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double _r3 = _r*_r*_r;
  const double alphar = ShortRange::alpha*r;
  double ee;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
  ee = EWALD_REAL_FORCE(alphar);
#else
  ee = M_2_SQRTPI*alphar*exp(-alphar*alphar) + erfc(alphar);
#endif
  de = _r3*ee;

  return de;
}

template<> inline double ForcePotential_Form_Full<ShortRange::ForEwaldReal>(double r2, double cutoff2, double& e)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double _r3 = _r*_r*_r;
  const double alphar = ShortRange::alpha*r;
  double ep,dep;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
  ewaldrealforceandpottbl_(&ep,&dep,&alphar);
#else
  ep = erfc(alphar);
  dep = M_2_SQRTPI*alphar*exp(-alphar*alphar) + ep;
#endif
  e = _r*ep;
  de = _r3*dep;

  return de;
}


template<> inline double Potential_Form_Full<ShortRange::ZeroDipole>(double r2, double cutoff2)
{
  double irc = 1.0/sqrt(cutoff2);
  double e;
  double _r = 1.0/sqrt(r2);
#ifdef ZERODIPOLE0
  e = _r+0.5*irc*(r2*irc*irc-3.0);
#else
  double r = r2*_r;
  const double alpharc2 = ShortRange::alpha*ShortRange::alpha*cutoff2;
  const double alpharc = ShortRange::alpha*cutoff2*irc;
  const double ec2 = 0.5*M_2_SQRTPI*ShortRange::alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)*irc;
  const double el = (ec1+ec2)*irc*irc;
  const double ec = -3.0*ec1 - ec2;
  const double ee = erfc(ShortRange::alpha*r);
  e = _r*ee + el*r2 + ec;
#endif
  return e;
}

template<> inline double Force_Form_Full<ShortRange::ZeroDipole>(double r2, double cutoff2)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
  double _r = 1.0/sqrt(r2);
  double _r3 = _r*_r*_r;
#ifdef ZERODIPOLE0
  double ec1 = 0.5*irc;
  double el = ec1*irc*irc;
  de = _r3-2.0*el;
#else
  double r = r2*_r;
  const double alpharc2 = ShortRange::alpha*ShortRange::alpha*cutoff2;
  const double alpharc = ShortRange::alpha*cutoff2*irc;
  const double alphar = ShortRange::alpha*r;
  const double ec2 = 0.5*M_2_SQRTPI*ShortRange::alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)*irc;
  const double el = (ec1+ec2)*irc*irc;
  const double ee = erfc(alphar);
  de = _r3*(M_2_SQRTPI*alphar*exp(-alphar*alphar)+ee) - 2.0*el;
#endif
  return de;
}

template<> inline double ForcePotential_Form_Full<ShortRange::ZeroDipole>(double r2, double cutoff2, double& e)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
  double _r = 1.0/sqrt(r2);
  double _r3 = _r*_r*_r;
#ifdef ZERODIPOLE0
  double ec1 = 0.5*irc;
  double el = ec1*irc*irc;
  e = _r + el*r2 - 3.0*ec1;
  de = _r3 -2.0*el;
#else
  double r = r2*_r;
  const double alpharc2 = ShortRange::alpha*ShortRange::alpha*cutoff2;
  const double alpharc = ShortRange::alpha*cutoff2*irc;
  const double alphar = ShortRange::alpha*r;
  const double ec2 = 0.5*M_2_SQRTPI*ShortRange::alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)*irc;
  const double el = (ec1+ec2)*irc*irc;
  const double ec = -3.0*ec1 - ec2;
  const double ee = erfc(alphar);
  e = _r*ee + el*r2 + ec;
  de = _r3*(M_2_SQRTPI*alphar*exp(-alphar*alphar)+ee) - 2.0*el;
#endif
  return de;
}

template<ShortRange::CoulombType T>
double energy_correction_fixed_water_full(double co, double ch, double cutoff2)
{
  double angleO = 104.520 * M_PI / 180.0;
  double lengthOH = 0.9572;
  double lengthHH = 2.0 * lengthOH * sin(0.5 * angleO);
  
  double eoh = co*ch*Potential_Form_Full<T>(lengthOH*lengthOH,cutoff2);
  double ehh = ch*ch*Potential_Form_Full<T>(lengthHH*lengthHH,cutoff2);
  return eoh*2.0 + ehh;
}


template<class PA, ShortRange::CoulombType T>
double energy_correction_excluded_water_full(PA& particle_array,
					const WaterList& waterlist,
					double cutoff2)
{
  double e=0.0;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
    Position ro = getpos(particle_array,o_index);
    Position rh1 = getpos(particle_array,h1_index);
    Position rh2 = getpos(particle_array,h2_index);
    double co = getcharge(particle_array,o_index);
    double ch1 = getcharge(particle_array,h1_index);
    double ch2 = getcharge(particle_array,h2_index);
    double r2oh1 = (ro-rh1).norm2();
    double r2oh2 = (ro-rh2).norm2();
    double r2h1h2 = (rh1-rh2).norm2();

//    e += co*(ch1*Potential_Form<T>(r2oh1,cutoff2) + ch2*Potential_Form<T>(r2oh2,cutoff2)) + ch1*ch2*Potential_Form<T>(r2h1h2,cutoff2);
    double eoh1 = co*ch1*Potential_Form_Full<T>(r2oh1,cutoff2);
    double eoh2 = co*ch2*Potential_Form_Full<T>(r2oh2,cutoff2);
    double eh1h2 = ch1*ch2*Potential_Form_Full<T>(r2h1h2,cutoff2);
    e += eoh1 + eoh2 + eh1h2;
#ifdef DUMP_ENERGY_CONTENT
    printf("correct water energy %e : OH1 %e , OH2 %e , H1H2 %e\n",eoh1 + eoh2 + eh1h2,eoh1,eoh2,eh1h2);
#endif
  }
  return e;
}

template<class PA, ShortRange::CoulombType T>
void force_correction_excluded_water_full(PA& particle_array,
				     const WaterList& waterlist,
				     double cutoff2)
{
  double e=0.0;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
    Position ro = getpos(particle_array,o_index);
    Position rh1 = getpos(particle_array,h1_index);
    Position rh2 = getpos(particle_array,h2_index);
    double co = getcharge(particle_array,o_index);
    double ch1 = getcharge(particle_array,h1_index);
    double ch2 = getcharge(particle_array,h2_index);
    double r2oh1 = (ro-rh1).norm2();
    double r2oh2 = (ro-rh2).norm2();
    double r2h1h2 = (rh1-rh2).norm2();

//    e += co*(ch1*Potential_Form<T>(r2oh1,cutoff2) + ch2*Potential_Form<T>(r2oh2,cutoff2)) + ch1*ch2*Potential_Form<T>(r2h1h2,cutoff2);
    double deoh1 = co*ch1*Force_Form_Full<T>(r2oh1,cutoff2);
    double deoh2 = co*ch2*Force_Form_Full<T>(r2oh2,cutoff2);
    double deh1h2 = ch1*ch2*Force_Form_Full<T>(r2h1h2,cutoff2);
    getforce(particle_array,o_index) += (deoh1*(ro-rh1)+deoh2*(ro-rh2));
    getforce(particle_array,h1_index) += (-deoh1*(ro-rh1)+deh1h2*(rh1-rh2));
    getforce(particle_array,h2_index) += (-deoh2*(ro-rh2)-deh1h2*(rh1-rh2));
#ifdef DUMP_ENERGY_CONTENT
    printf("correct water denergy  OH1 %e , OH2 %e , H1H2 %e\n",deoh1,deoh2,deh1h2);
#endif
  }
}
				

template<class PA, ShortRange::CoulombType T>
double forceenergy_correction_excluded_water_full(PA& particle_array,
					     const WaterList& waterlist,
					     double cutoff2)
{
  double e=0.0;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
    Position ro = getpos(particle_array,o_index);
    Position rh1 = getpos(particle_array,h1_index);
    Position rh2 = getpos(particle_array,h2_index);
    double co = getcharge(particle_array,o_index);
    double ch1 = getcharge(particle_array,h1_index);
    double ch2 = getcharge(particle_array,h2_index);
    double r2oh1 = (ro-rh1).norm2();
    double r2oh2 = (ro-rh2).norm2();
    double r2h1h2 = (rh1-rh2).norm2();

//    e += co*(ch1*Potential_Form<T>(r2oh1,cutoff2) + ch2*Potential_Form<T>(r2oh2,cutoff2)) + ch1*ch2*Potential_Form<T>(r2h1h2,cutoff2);
    double eoh1, eoh2, eh1h2;
    double deoh1 = co*ch1*ForcePotential_Form_Full<T>(r2oh1,cutoff2,eoh1);
    double deoh2 = co*ch2*ForcePotential_Form_Full<T>(r2oh2,cutoff2,eoh2);
    double deh1h2 = ch1*ch2*ForcePotential_Form_Full<T>(r2h1h2,cutoff2,eh1h2);
    getforce(particle_array,o_index) += (deoh1*(ro-rh1)+deoh2*(ro-rh2));
    getforce(particle_array,h1_index) += (-deoh1*(ro-rh1)+deh1h2*(rh1-rh2));
    getforce(particle_array,h2_index) += (-deoh2*(ro-rh2)-deh1h2*(rh1-rh2));
    e += co*ch1*eoh1 + co*ch2*eoh2 + ch1*ch2*eh1h2;
#ifdef DUMP_ENERGY_CONTENT
    printf("correct water energy %e : OH1 %e , OH2 %e , H1H2 %e\n",eoh1 + eoh2 + eh1h2,eoh1,eoh2,eh1h2);
    printf("correct water denergy  OH1 %e , OH2 %e , H1H2 %e\n",eoh1,eoh2,eh1h2);
#endif
  }
  return e;
}

template<ShortRange::CoulombType T>
inline double Potential_Form(double r2, double cutoff2)
{
}

template<ShortRange::CoulombType T>
inline double Force_Form(double r2, double cutoff2)
{
}

template<ShortRange::CoulombType T>
inline double ForcePotential_Form(double r2, double cutoff2, double& e)
{
}

template<> inline double Potential_Form<ShortRange::ForEwaldReal>(double r2, double cutoff2)
{
  double r = sqrt(r2);
  double _r = 1.0/r;
  double ep;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
  ep = EWALD_REAL_POT(ShortRange::alpha*r)-1.0;
#else
  ep = -erf(ShortRange::alpha*r);
#endif
  double correct_ewaldpot = _r*ep;
#ifdef DUMP_ENERGY_CONTENT
  printf("correct ewald %e : r2 = %e, r = %e, alpha = %e, alpha*r = %e, ewald pot %e, erf(ar) %e\n",correct_ewaldpot,r2,r,ShortRange::alpha,ShortRange::alpha*r,ep,erf(ShortRange::alpha*r));
#endif
  return correct_ewaldpot;
}

template<> inline double Force_Form<ShortRange::ForEwaldReal>(double r2, double cutoff2)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double _r3 = _r*_r*_r;
  const double alphar = ShortRange::alpha*r;
  double ee;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
  ee = EWALD_REAL_FORCE(alphar)-1.0;
#else
  ee = M_2_SQRTPI*alphar*exp(-alphar*alphar) - erf(alphar);
#endif
  de = _r3*ee;

  return de;
}

template<> inline double ForcePotential_Form<ShortRange::ForEwaldReal>(double r2, double cutoff2, double& e)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double _r3 = _r*_r*_r;
  const double alphar = ShortRange::alpha*r;
  double ep,dep;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
  ewaldrealforceandpottbl_(&ep,&dep,&alphar);
  ep -= 1.0;
  dep -= 1.0;
#else
  ep = -erf(alphar);
  dep = M_2_SQRTPI*alphar*exp(-alphar*alphar) + ep;
#endif
  e = _r*ep;
  de = _r3*dep;

  return de;
}


template<> inline double Potential_Form<ShortRange::ZeroDipole>(double r2, double cutoff2)
{
  double irc = 1.0/sqrt(cutoff2);
  double e;
#ifdef ZERODIPOLE0
  e = 0.5*irc*(r2*irc*irc-3.0);
#else
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  const double alpharc2 = ShortRange::alpha*ShortRange::alpha*cutoff2;
  const double alpharc = ShortRange::alpha*cutoff2*irc;
  const double ec2 = 0.5*M_2_SQRTPI*ShortRange::alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)*irc;
  const double el = (ec1+ec2)*irc*irc;
  const double ec = -3.0*ec1 - ec2;
  const double ee = -erf(ShortRange::alpha*r);
  e = _r*ee + el*r2 + ec;
#endif
  return e;
}

template<> inline double Force_Form<ShortRange::ZeroDipole>(double r2, double cutoff2)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
#ifdef ZERODIPOLE0
  double ec1 = 0.5*irc;
  double el = ec1*irc*irc;
  de = -2.0*el;
#else
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double _r3 = _r*_r*_r;
  const double alpharc2 = ShortRange::alpha*ShortRange::alpha*cutoff2;
  const double alpharc = ShortRange::alpha*cutoff2*irc;
  const double alphar = ShortRange::alpha*r;
  const double ec2 = 0.5*M_2_SQRTPI*ShortRange::alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)*irc;
  const double el = (ec1+ec2)*irc*irc;
  const double ee = -erf(alphar);
  de = _r3*(M_2_SQRTPI*alphar*exp(-alphar*alphar)+ee) - 2.0*el;
#endif
  return de;
}

template<> inline double ForcePotential_Form<ShortRange::ZeroDipole>(double r2, double cutoff2, double& e)
{
  double irc = 1.0/sqrt(cutoff2);
  double de;
#ifdef ZERODIPOLE0
  double ec1 = 0.5*irc;
  double el = ec1*irc*irc;
  e = el*r2 - 3.0*ec1;
  de = -2.0*el;
#else
  double _r = 1.0/sqrt(r2);
  double r = r2*_r;
  double _r3 = _r*_r*_r;
  const double alpharc2 = ShortRange::alpha*ShortRange::alpha*cutoff2;
  const double alpharc = ShortRange::alpha*cutoff2*irc;
  const double alphar = ShortRange::alpha*r;
  const double ec2 = 0.5*M_2_SQRTPI*ShortRange::alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)*irc;
  const double el = (ec1+ec2)*irc*irc;
  const double ec = -3.0*ec1 - ec2;
  const double ee = -erf(alphar);
  e = _r*ee + el*r2 + ec;
  de = _r3*(M_2_SQRTPI*alphar*exp(-alphar*alphar)+ee) - 2.0*el;
#endif
  return de;
}


template<class PA, ShortRange::CoulombType T>
double energy_correction_excluded_water(PA& particle_array,
					const WaterList& waterlist,
					double cutoff2)
{
  double e=0.0;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
    Position ro = getpos(particle_array,o_index);
    Position rh1 = getpos(particle_array,h1_index);
    Position rh2 = getpos(particle_array,h2_index);
    double co = getcharge(particle_array,o_index);
    double ch1 = getcharge(particle_array,h1_index);
    double ch2 = getcharge(particle_array,h2_index);
    double r2oh1 = (ro-rh1).norm2();
    double r2oh2 = (ro-rh2).norm2();
    double r2h1h2 = (rh1-rh2).norm2();

//    e += co*(ch1*Potential_Form<T>(r2oh1,cutoff2) + ch2*Potential_Form<T>(r2oh2,cutoff2)) + ch1*ch2*Potential_Form<T>(r2h1h2,cutoff2);
    double eoh1 = co*ch1*Potential_Form<T>(r2oh1,cutoff2);
    double eoh2 = co*ch2*Potential_Form<T>(r2oh2,cutoff2);
    double eh1h2 = ch1*ch2*Potential_Form<T>(r2h1h2,cutoff2);
    e += eoh1 + eoh2 + eh1h2;
#ifdef DUMP_ENERGY_CONTENT
    printf("correct water energy %e : OH1 %e , OH2 %e , H1H2 %e\n",eoh1 + eoh2 + eh1h2,eoh1,eoh2,eh1h2);
#endif
  }
  return e;
}

template<class PA, ShortRange::CoulombType T>
void force_correction_excluded_water(PA& particle_array,
				     const WaterList& waterlist,
				     double cutoff2)
{
  double e=0.0;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
    Position ro = getpos(particle_array,o_index);
    Position rh1 = getpos(particle_array,h1_index);
    Position rh2 = getpos(particle_array,h2_index);
    double co = getcharge(particle_array,o_index);
    double ch1 = getcharge(particle_array,h1_index);
    double ch2 = getcharge(particle_array,h2_index);
    double r2oh1 = (ro-rh1).norm2();
    double r2oh2 = (ro-rh2).norm2();
    double r2h1h2 = (rh1-rh2).norm2();

//    e += co*(ch1*Potential_Form<T>(r2oh1,cutoff2) + ch2*Potential_Form<T>(r2oh2,cutoff2)) + ch1*ch2*Potential_Form<T>(r2h1h2,cutoff2);
    double deoh1 = co*ch1*Force_Form<T>(r2oh1,cutoff2);
    double deoh2 = co*ch2*Force_Form<T>(r2oh2,cutoff2);
    double deh1h2 = ch1*ch2*Force_Form<T>(r2h1h2,cutoff2);
    getforce(particle_array,o_index) += (deoh1*(ro-rh1)+deoh2*(ro-rh2));
    getforce(particle_array,h1_index) += (-deoh1*(ro-rh1)+deh1h2*(rh1-rh2));
    getforce(particle_array,h2_index) += (-deoh2*(ro-rh2)-deh1h2*(rh1-rh2));
#ifdef DUMP_ENERGY_CONTENT
    printf("correct water denergy  OH1 %e , OH2 %e , H1H2 %e\n",deoh1,deoh2,deh1h2);
#endif
  }
}
				

template<class PA, ShortRange::CoulombType T>
double forceenergy_correction_excluded_water(PA& particle_array,
					     const WaterList& waterlist,
					     double cutoff2)
{
  double e=0.0;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
    Position ro = getpos(particle_array,o_index);
    Position rh1 = getpos(particle_array,h1_index);
    Position rh2 = getpos(particle_array,h2_index);
    double co = getcharge(particle_array,o_index);
    double ch1 = getcharge(particle_array,h1_index);
    double ch2 = getcharge(particle_array,h2_index);
    double r2oh1 = (ro-rh1).norm2();
    double r2oh2 = (ro-rh2).norm2();
    double r2h1h2 = (rh1-rh2).norm2();

//    e += co*(ch1*Potential_Form<T>(r2oh1,cutoff2) + ch2*Potential_Form<T>(r2oh2,cutoff2)) + ch1*ch2*Potential_Form<T>(r2h1h2,cutoff2);
    double eoh1, eoh2, eh1h2;
    double deoh1 = co*ch1*ForcePotential_Form<T>(r2oh1,cutoff2,eoh1);
    double deoh2 = co*ch2*ForcePotential_Form<T>(r2oh2,cutoff2,eoh2);
    double deh1h2 = ch1*ch2*ForcePotential_Form<T>(r2h1h2,cutoff2,eh1h2);
    getforce(particle_array,o_index) += (deoh1*(ro-rh1)+deoh2*(ro-rh2));
    getforce(particle_array,h1_index) += (-deoh1*(ro-rh1)+deh1h2*(rh1-rh2));
    getforce(particle_array,h2_index) += (-deoh2*(ro-rh2)-deh1h2*(rh1-rh2));
    e += co*ch1*eoh1 + co*ch2*eoh2 + ch1*ch2*eh1h2;
#ifdef DUMP_ENERGY_CONTENT
    printf("correct water energy %e : OH1 %e , OH2 %e , H1H2 %e\n",eoh1 + eoh2 + eh1h2,eoh1,eoh2,eh1h2);
    printf("correct water denergy  OH1 %e , OH2 %e , H1H2 %e\n",eoh1,eoh2,eh1h2);
#endif
  }
  return e;
}
				

#endif

