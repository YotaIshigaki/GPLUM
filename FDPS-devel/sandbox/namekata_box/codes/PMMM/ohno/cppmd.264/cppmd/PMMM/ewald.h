#pragma once

#include <vector>
#include "Common.h"



typedef std::complex<double> complex;

static inline complex expi(double theta){
	return complex(cos(theta), sin(theta));
}

struct Waves{
  std::vector<complex> wx;
  std::vector<complex> wy;
  std::vector<complex> wz;
  int mmax;

Waves(const Position &pos, int _mmax=5) : mmax(_mmax), wx(2*_mmax+1), wy(2*_mmax+1), wz(2*_mmax+1) {
  const double twopi = 8.0 * atan(1.0);
  const complex ex = expi(twopi * pos.x);
  const complex ey = expi(twopi * pos.y);
  const complex ez = expi(twopi * pos.z);

  complex emx(1.0, 0.0);
  complex emy(1.0, 0.0);
  complex emz(1.0, 0.0);
		
  for(int m=0; m<=mmax; m++){
    wx[mmax + m] = emx;
    wx[mmax - m] = conj(emx);
    emx *= ex;
    wy[mmax + m] = emy;
    wy[mmax - m] = conj(emy);
    emy *= ey;
    wz[mmax + m] = emz;
    wz[mmax - m] = conj(emz);
    emz *= ez;
  }
}
};

struct EwaldSum{
  int nx;
  int ny;
  int nz;
  std::vector<complex> sum;
  int mmax;
  
EwaldSum(int _mmax=5) : mmax(_mmax), nx(2*_mmax+1), ny(2*_mmax+1), nz(2*_mmax+1) {
  sum.resize(nx*ny*nz);
}
  void clear(){
    for(int i=0; i<nx*ny*nz; i++){
      sum[i] = complex(0.0, 0.0);
    }
  }

  void assign(const ParticlePosCharge &p){
    const double q = p.charge;
    const Position  r = p.position;
    Waves  w(r, mmax);
    for(int iz=0;  iz<nz; iz++){
      for(int iy=0;  iy<ny; iy++){
	for(int ix=0;  ix<nx; ix++){
	  sum[(iz*ny+iy)*nx+ix] += ((q * w.wz[iz]) * w.wy[iy]) * w.wx[ix];
	}
      }
    }
  }

  void filter(const double alpha){
    const double pi = 4.0 * atan(1.0);
    const double coef = (pi/alpha) * (pi/alpha);
    for(int iz=0;  iz<nz; iz++){
      for(int iy=0;  iy<ny; iy++){
	for(int ix=0;  ix<nx; ix++){
	  const int mx = ix - mmax;
	  const int my = iy - mmax;
	  const int mz = iz - mmax;
	  const double m2 = double(mx*mx + my*my + mz*mz);
	  const double factor = exp(-coef * m2) / m2;
	  sum[(iz*ny+iy)*nx+ix] *= factor;
	}
      }
    }
    sum[(mmax*ny+mmax)*nx+mmax] = 0.0;
  }

  void phi_and_grad(ParticlePosCharge &p, Force &acc_direct, double &phi_direct) const
  {
    const Position  r = p.position;
    Waves w(r,mmax);
    double phi (0.0);
    Position  grad(0.0);
    for(int iz=0;  iz<nz; iz++){
      for(int iy=0;  iy<ny; iy++){
	for(int ix=0;  ix<nx; ix++){
	  const int mx = ix - mmax;
	  const int my = iy - mmax;
	  const int mz = iz - mmax;
	  const Position mvec = Position(double(mx), double(my), double(mz));
	  const complex cs = (w.wz[iz] * w.wy[iy]) * w.wx[ix];
	  const double c = real(cs);
	  const double s = imag(cs);
	  const double csum = real(sum[(iz*ny+iy)*nx+ix]);
	  const double ssum = imag(sum[(iz*ny+iy)*nx+ix]);

	  phi  += c*csum + s*ssum;
	  grad -= (s*csum - c*ssum) * mvec;
	}
      }
    }
    const double pi = 4.0 * atan(1.0);
    phi_direct += (1.0/pi) * phi;
    acc_direct += -2.0 * grad;
  }
};

void eval_k_space(
		  const int    N,
		  const double alpha,
		  PosChargeArray    & ptcl,
		  ForceArray &acc,
		  std::vector<double> &phi,
		  const int mmax=5)
{
  static EwaldSum sum(mmax);

  sum.clear();
  for(int i=0; i<N; i++){
    sum.assign(ptcl[i]);
  }
  sum.filter(alpha);
  for(int i=0; i<N; i++){
    sum.phi_and_grad(ptcl[i], acc[i], phi[i]);
  }
}

inline void cutfunc(
		    const double r,
		    const double alpha,
		    const double rsq,
		    const double asq,
		    double &pcut,
		    double &fcut)
{
  const double c = 2.0 / sqrt(4.0 * atan(1.0)); // 2/sqrt(pi)  M_2_SQRTPI
  const double tmp = erfc(alpha * r);
  pcut = tmp;
  fcut = tmp + c * (alpha * r) * exp(-asq * rsq);
}

static inline Position minimum_image(const Position &inp){
  return Position(
		  inp.x - round(inp.x),
		  inp.y - round(inp.y),
		  inp.z - round(inp.z));
}

void eval_r_space(
		  const int    N,
		  const double alpha,
		  const double msum,
		  PosChargeArray  &ptcl,
		  ForceArray &acc_direct,
		  std::vector<double> &phi_direct,
		  const int nmir=3)
{
  const double ainv2 = (4.0 * atan(1.0)) / (alpha * alpha);
#pragma omp parallel for
  for(int i=0; i<N; i++){
    if(N>1000){
      if(i%10==0){
	printf("eval_r_space %d in %d\n",i,N);
      }
    }
    double phi(0.0);
    Position  grad(0.0);
    for(int mz=-nmir; mz<=nmir; mz++){
      for(int my=-nmir; my<=nmir; my++){
	for(int mx=-nmir; mx<=nmir; mx++){
	  const Position off = Position(double(mx), double(my), double(mz));
	  double ptmp(0.0);
	  Position  gtmp(0.0);
	  for(int j=0; j<N; j++){
	    Position dr = ptcl[j].position - ptcl[i].position;
	    dr = minimum_image(dr) + off;
	    
	    const double r2 = dr * dr;
	    if(0.0 == r2) continue;
	    const double r  = sqrt(r2);
	    const double rinv = 1.0 / r;
	    const double qri  = ptcl[j].charge * rinv;
	    const double qri3 = qri * (rinv * rinv);
	    
	    double pcut, fcut;
	    cutfunc(r, alpha, r2, alpha*alpha, pcut, fcut);
	    
	    ptmp += pcut * qri;
	    gtmp -= (fcut * qri3) * dr;
	  } // for(j)
	  phi  += ptmp;
	  grad += gtmp;
	}
      }
    }
    // self energy
    phi -= alpha * (2.0/sqrt(4.0 * atan(1.0))) * ptcl[i].charge;
    phi -= msum * ainv2; // charged system correction
    phi_direct[i] += phi;
    acc_direct[i] += grad;
    //    printf("%e %e %e %e\n",phi,grad.x,grad.y,grad.z);
  } // for(i)
}

void brute_nbody(
		 const int    N,
		 PosChargeArray  &ptcl,
		 ForceArray &acc_direct,
		 std::vector<double> &phi_direct,
		 double &e0,
		 const int nmir=3)
{
  double ee0 = 0.0;
#pragma omp parallel for reduction(+:ee0)
  for(int i=0; i<N; i++){
    double phi(0.0);
    Position  grad(0.0);
    for(int mz=-nmir; mz<=nmir; mz++){
      for(int my=-nmir; my<=nmir; my++){
	for(int mx=-nmir; mx<=nmir; mx++){
	  const Position off = Position(double(mx), double(my), double(mz));
	  double ptmp(0.0);
	  Position  gtmp(0.0);
	  for(int j=0; j<N; j++){
	    Position dr = ptcl[j].position - ptcl[i].position;
	    // minimum image
	    dr.x -= round(dr.x);
	    dr.y -= round(dr.y);
	    dr.z -= round(dr.z);

	    dr += off;
	    
	    const double r2 = dr * dr;
	    if(0.0 == r2) continue;
	    const double r  = sqrt(r2);
	    const double rinv = 1.0 / r;
	    const double qri  = ptcl[j].charge * rinv;
	    const double qri3 = qri * (rinv * rinv);
	    
	    ptmp += qri;
	    gtmp -= qri3 * dr;
	  } // for(j)
	  if((mx==0)&&(my==0)&&(mz==0))ee0+=ptmp*ptcl[i].charge;
	  phi  += ptmp;
	  grad += gtmp;
	}
      }
    }
    phi_direct[i] += phi;
    acc_direct[i] += grad;
  } // for(i)
  e0 = ee0;
}

/*
  f = q/|dr|^3 * dr = F /rs^2
  p = q/|dr|        = P /rs

  dR = dr/rs
  F = q/|dR|^3 = q/|dr/rs|^3 * dr/rs = q/|dr|^3 * dr * rs^2 = f * rs^2
  P = q/|dR|   = q/|dr/rs|           = q/|dr| *rs           = p * rs



  f = (erfc(alpha * dr) + M_2_SQRTPI*(alpha * dr) * exp(-alpha^2 * dr^2)) * q/|dr|^3 * dr;
  p = (erfc(alpha * dr)) * q/|dr|

  P = (erfc(Alpha*dR))*q/|dR| = (erfc(alpha*rs*dr/rs))*q/|dr/rs| = erfc(alpha*dr)*q/|dr| * rs = P*rs

 */
