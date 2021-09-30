#pragma once
#include "particle.h"

typedef std::complex<double> complex;

static inline complex expi(double theta){
	return complex(cos(theta), sin(theta));
}

template <int MMAX>
struct Waves{
	complex wx[2*MMAX+1];
	complex wy[2*MMAX+1];
	complex wz[2*MMAX+1];

	Waves(const dvec3 &pos){
		const double twopi = 8.0 * atan(1.0);
		const complex ex = expi(twopi * pos.x);
		const complex ey = expi(twopi * pos.y);
		const complex ez = expi(twopi * pos.z);

		complex emx(1.0, 0.0);
		complex emy(1.0, 0.0);
		complex emz(1.0, 0.0);
		
		for(int m=0; m<=MMAX; m++){
			wx[MMAX + m] = emx;
			wx[MMAX - m] = conj(emx);
			emx *= ex;
			wy[MMAX + m] = emy;
			wy[MMAX - m] = conj(emy);
			emy *= ey;
			wz[MMAX + m] = emz;
			wz[MMAX - m] = conj(emz);
			emz *= ez;
		}
	}
};

template <int MMAX>
struct EwaldSum{
	enum{
		NX = 2*MMAX+1,
		NY = 2*MMAX+1,
		NZ = 2*MMAX+1,
	};
	complex sum[NZ][NY][NX];

	void clear(){
		for(int i=0; i<NX*NY*NZ; i++){
			sum[0][0][i] = complex(0.0, 0.0);
		}
	}

	void assign(const Particle &p){
		const double q = p.mass;
		const dvec3  r = p.pos;
		Waves<MMAX> w(r);
		for(int iz=0;  iz<NZ; iz++){
			for(int iy=0;  iy<NY; iy++){
				for(int ix=0;  ix<NX; ix++){
					sum[iz][iy][ix] += ((q * w.wz[iz]) * w.wy[iy]) * w.wx[ix];
				}
			}
		}
	}

	void filter(const double alpha){
		const double pi = 4.0 * atan(1.0);
		const double coef = (pi/alpha) * (pi/alpha);
		for(int iz=0;  iz<NZ; iz++){
			for(int iy=0;  iy<NY; iy++){
				for(int ix=0;  ix<NX; ix++){
					const int mx = ix - MMAX;
					const int my = iy - MMAX;
					const int mz = iz - MMAX;
					const double m2 = double(mx*mx + my*my + mz*mz);
					const double factor = exp(-coef * m2) / m2;
					sum[iz][iy][ix] *= factor;
				}
			}
		}
		sum[MMAX][MMAX][MMAX] = 0.0;
	}

	void phi_and_grad(Particle &p) const
	{
		const dvec3  r = p.pos;
		Waves<MMAX> w(r);
		double phi (0.0);
		dvec3  grad(0.0);
		for(int iz=0;  iz<NZ; iz++){
			for(int iy=0;  iy<NY; iy++){
				for(int ix=0;  ix<NX; ix++){
					const int mx = ix - MMAX;
					const int my = iy - MMAX;
					const int mz = iz - MMAX;
					const dvec3 mvec = dvec3(double(mx), double(my), double(mz));
					const complex cs = (w.wz[iz] * w.wy[iy]) * w.wx[ix];
					const double c = real(cs);
					const double s = imag(cs);
					const double csum = real(sum[iz][iy][ix]);
					const double ssum = imag(sum[iz][iy][ix]);

					phi  += c*csum + s*ssum;
					grad += (s*csum - c*ssum) * mvec;
				}
			}
		}
		const double pi = 4.0 * atan(1.0);
		p.phi_direct += (1.0/pi) * phi;
		p.acc_direct += -2.0 * grad;
	}
};

template <int MMAX>
void eval_k_space(
		const int    N,
		const double alpha,
		Particle     ptcl[])
{
	static EwaldSum<MMAX> sum;

	sum.clear();
	for(int i=0; i<N; i++){
		sum.assign(ptcl[i]);
	}
	sum.filter(alpha);
	for(int i=0; i<N; i++){
		sum.phi_and_grad(ptcl[i]);
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
	const double c = 2.0 / sqrt(4.0 * atan(1.0));
	const double tmp = erfc(alpha * r);
	pcut = tmp;
	fcut = tmp + c * (alpha * r) * exp(-asq * rsq);
}

template <int NMIR>
void eval_r_space(
		const int    N,
		const double alpha,
		const double msum,
		Particle     ptcl[])
{
	const double ainv2 = (4.0 * atan(1.0)) / (alpha * alpha);
#pragma omp parallel for
	for(int i=0; i<N; i++){
		double phi(0.0);
		dvec3  grad(0.0);
		for(int mz=-NMIR; mz<=NMIR; mz++){
			for(int my=-NMIR; my<=NMIR; my++){
				for(int mx=-NMIR; mx<=NMIR; mx++){
					const dvec3 off = dvec3(double(mx), double(my), double(mz));
					double ptmp(0.0);
					dvec3  gtmp(0.0);
					for(int j=0; j<N; j++){
						dvec3 dr = ptcl[j].pos - ptcl[i].pos;
						dr = minimum_image(dr) + off;

						const double r2 = dr * dr;
						if(0.0 == r2) continue;
						const double r  = sqrt(r2);
						const double rinv = 1.0 / r;
						const double qri  = ptcl[j].mass * rinv;
						const double qri3 = qri * (rinv * rinv);

						double pcut, fcut;
						cutfunc(r, alpha, r2, alpha*alpha, pcut, fcut);

						ptmp += pcut * qri;
						gtmp += (fcut * qri3) * dr;
					} // for(j)
					phi  += ptmp;
					grad += gtmp;
				}
			}
		}
		// self energy
		phi -= alpha * (2.0/sqrt(4.0 * atan(1.0))) * ptcl[i].mass;
		phi -= msum * ainv2; // charged system correction
		ptcl[i].phi_direct += phi;
		ptcl[i].acc_direct += grad;
	} // for(i)
}

template <int NMIR>
void brute_nbody(
		const int    N,
		Particle     ptcl[])
{
#pragma omp parallel for
	for(int i=0; i<N; i++){
		double phi(0.0);
		dvec3  grad(0.0);
		for(int mz=-NMIR; mz<=NMIR; mz++){
			for(int my=-NMIR; my<=NMIR; my++){
				for(int mx=-NMIR; mx<=NMIR; mx++){
					const dvec3 off = dvec3(double(mx), double(my), double(mz));
					double ptmp(0.0);
					dvec3  gtmp(0.0);
					for(int j=0; j<N; j++){
						dvec3 dr = ptcl[j].pos - ptcl[i].pos;
						// minimum image
						dr.x -= round(dr.x);
						dr.y -= round(dr.y);
						dr.z -= round(dr.z);

						dr += off;

						const double r2 = dr * dr;
						if(0.0 == r2) continue;
						const double r  = sqrt(r2);
						const double rinv = 1.0 / r;
						const double qri  = ptcl[j].mass * rinv;
						const double qri3 = qri * (rinv * rinv);

						ptmp += qri;
						gtmp += qri3 * dr;
					} // for(j)
					phi  += ptmp;
					grad += gtmp;
				}
			}
		}
		ptcl[i].phi_direct += phi;
		ptcl[i].acc_direct += grad;
	} // for(i)
}

