#include <cmath>
#include <complex>

template <int LMAX>
struct CutFunc{
	typedef std::complex<double> cplx_t;
	double gammmainv[1+LMAX];
	double rcut     [1+LMAX];
	cplx_t kcut     [1+LMAX];

	CutFunc(){
		double val = sqrt(4.0*atan(1.0));
		for(int ell=0; ell<=LMAX; ell++){
			gammmainv[ell] = 1.0 / val;
			val *= (0.5 + ell);
		}
	}

	void eval_rcut(const double x){
		const double xexp = exp(-x);
		double xpow = sqrt(x);
		double val  = sqrt(4.0*atan(1.0)) * erfc(xpow);
		for(int ell=0; ell<=LMAX; ell++){
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
		for(int ell=0; ell<=LMAX; ell++){
			kcut[ell] = coef * (gauss * gammmainv[ell]);
			coef *= ratio;
		}
	}
};
