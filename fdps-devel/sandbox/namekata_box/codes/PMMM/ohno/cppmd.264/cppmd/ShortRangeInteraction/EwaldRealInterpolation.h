#ifndef EWALD_REAL_INTERPOLATION_H
#define EWALD_REAL_INTERPOLATION_H

#include <iostream>
#include <cmath>

extern struct EwaldRealInterpolation ewaldRealInterpolation;
struct ForceAndPot {
  double force;
  double pot;
};

struct EwaldRealInterpolation {
  static const double ewald_real_force_table[320];
  static const double ewald_real_pot_table[320];

  inline double ewaldRealForceTbl(double x)
  {
    int ix = (int)(x*4.0);
    const double *c=ewald_real_force_table+ix*5;
    double dx;

    // Note: the following line may be unnecessary for normal situations.
    if ((x<0.0)||(x>=16.0)) {
//      std::cerr << "EwaldRealInterpolation::ewaldRealForceTbl Range Error" << x << '\n';
      return 0.0;
    }
    dx = (x*4.0-ix)*2.0-1.0;
    return c[0]+dx*(c[1]+dx*(c[2]+dx*(c[3]+dx*c[4])));
  }

  inline double ewaldRealPotTbl(double x)
  {
    int ix = int(x*4.0);
    const double *c=ewald_real_pot_table+ix*5;
    double dx;

    // Note: the following line may be unnecessary for normal situations.
    if ((x<0.0)||(x>=16.0)) return 0.0;

    dx = (x*4.0-ix)*2.0-1.0;
    return c[0]+dx*(c[1]+dx*(c[2]+dx*(c[3]+dx*c[4])));
  }

  inline ForceAndPot ewaldRealForceAndPotTbl(double x)
  {
    ForceAndPot ret={0.0,0.0};
    int ix = int(x*4.0);
    const double *c=ewald_real_force_table+ix*5;
    const double *pc=ewald_real_pot_table+ix*5;

    // Note: the following line may be unnecessary for normal situations.
    if ((x<0.0)||(x>=16.0)) return ret;

    double dx = (x*4.0-ix)*2.0-1.0;
    ret.force = c[0]+dx*(c[1]+dx*(c[2]+dx*(c[3]+dx*c[4])));
    ret.pot = pc[0]+dx*(pc[1]+dx*(pc[2]+dx*(pc[3]+dx*pc[4])));
    return ret;
  }
};

// EwaldRealInterpolation.ewaldRealForceAndPotTbl for Fortran kernel
extern "C"
inline
void ewaldrealforceandpottbl_(double *pot, double *dp, const double *x)
{
  
  int ix = int((*x)*4.0);
  const double *c=ewaldRealInterpolation.ewald_real_force_table+ix*5;
  const double *pc=ewaldRealInterpolation.ewald_real_pot_table+ix*5;

  // Note: the following line may be unnecessary for normal situations.
  //  if ((x<0.0)||(x>=16.0)) return ret;

  double dx = ((*x)*4.0-ix)*2.0-1.0;
  *dp = c[0]+dx*(c[1]+dx*(c[2]+dx*(c[3]+dx*c[4])));
  *pot = pc[0]+dx*(pc[1]+dx*(pc[2]+dx*(pc[3]+dx*pc[4])));
};

#endif
