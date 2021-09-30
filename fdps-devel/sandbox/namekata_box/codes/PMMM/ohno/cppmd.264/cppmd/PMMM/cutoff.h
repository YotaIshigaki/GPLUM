#pragma once

#include <cmath>
#include <complex>
#include <vector>

struct CutFunc{
  typedef std::complex<double> cplx_t;
  int lmax;
  std::vector<double> gammmainv;
  std::vector<double> rcut;
  std::vector<cplx_t> kcut;
  

CutFunc(int p2=10) : lmax(p2) {
    gammmainv.resize(lmax+1);
    rcut.resize(lmax+1);
    kcut.resize(lmax+1);
    double val = sqrt(4.0*atan(1.0));
    for(int ell=0; ell<=lmax; ell++){
      gammmainv[ell] = 1.0 / val;
      val *= (0.5 + ell);
    }
  }

    void eval_rcut(const double x){
      const double xexp = exp(-x);
      double xpow = sqrt(x);
      double val  = sqrt(4.0*atan(1.0)) * erfc(xpow);
      for(int ell=0; ell<=lmax; ell++){
	rcut[ell] = val * gammmainv[ell];
	val = (0.5 + ell)*val + xpow * xexp;
	xpow *= x;
      }
    }

  void eval_kcut(const double kn2, const double kappa){
    const double pi = 4.0 * atan(1.0);
    const double gauss = exp(-kn2 * ((pi/kappa)*(pi/kappa)));

    const cplx_t ratio(0.0, -pi * kn2);
    cplx_t coef = 1.0 / sqrt(pi * kn2);
    for(int ell=0; ell<=lmax; ell++){
      kcut[ell] = coef * (gauss * gammmainv[ell]);
      coef *= ratio;
    }
  }
};
