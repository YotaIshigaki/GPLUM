#pragma once

// PMMM for periodic boundary condition
#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
/*
#include "vector3.h"

#include "cell.h"
#include "ewald.h"
*/

#include "fmm.h"
#include "cutoff.h"

template<int p, int NX, int NY, int NZ, int NMAX, int MMAX, int ICUT>
struct GreenFunction_PBC{
	typedef double               real_t;
	typedef std::complex<real_t> cplx_t;
	// typedef double _Complex      cplx_t;
	enum{
		LEN  = lmbuf<p>::length,
		LEN2 = lmbuf<2*p>::length,
	};

	real_t gf_r[NZ][NY][NX]    [LEN2];
	cplx_t gf_k[NZ][NY][1+NX/2][LEN2];


	void gen_gf_k(){
		static real_t rbuf[NZ][NY][NX];
		static cplx_t kbuf[NZ][NY][1+NX/2];
		fftw_plan plan_fwd = 
			fftw_plan_dft_r2c_3d(
				NZ, NY, NX, 
				(double       *)(rbuf),
				(fftw_complex *)(kbuf),
				FFTW_ESTIMATE);
		for(int lm=0; lm<LEN2; lm++){
			int i, j, k;
			for(k=0; k<NZ; k++) for(int j=0; j<NY; j++) for(i=0; i<NX; i++)
			{
				rbuf[k][j][i] = gf_r[k][j][i] [lm];
			}
			// CALL FFTW
			fftw_execute(plan_fwd);

			for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
			{
				gf_k[k][j][i] [lm] = kbuf[k][j][i];
			}
		}
		fftw_destroy_plan(plan_fwd);
	}

	typedef MultipoleMoment <p, cplx_t> mm_t;
	typedef LocalExpansion  <p, cplx_t> le_t;

	void transform(
			const mm_t mm_k[NZ][NY][1+NX/2],
			      le_t le_k[NZ][NY][1+NX/2]) const
	{
		int i, j, k;
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
		{
			typedef Slm<2*p, cplx_t> slm_t;
			((slm_t *)(gf_k[k][j][i]))
				-> template transform_M2L<p, p, false>(
						mm_k[k][j][i], le_k[k][j][i]);
		}
	}

  //  template <int NMAX, int MMAX, int ICUT>
	  void gen_gf_r(const double alpha, const double cell_length)
  {
		{
			Slm<2*p, real_t> dry;
			dry.init();
		}
#pragma omp parallel for
		for(int k=0; k<NZ; k++) for(int j=0; j<NY; j++) for(int i=0; i<NX; i++)
		{
			CutFunc<2*p> cf;
			// real-space sum
			Slm<2*p, real_t> rsum;
			rsum.clear();

			const int kk = (k>=NZ/2) ? k - NZ : k;
			const int jj = (j>=NY/2) ? j - NY : j;
			const int ii = (i>=NX/2) ? i - NX : i;
			int nx, ny, nz;
			for(nz=-NMAX; nz<=NMAX; nz++) for(ny=-NMAX; ny<=NMAX; ny++) for(nx=-NMAX; nx<=NMAX; nx++)
			{
				const int kkk = kk + nz*NZ;
				const int jjj = jj + ny*NY;
				const int iii = ii + nx*NX;
				if( 0 == (iii|jjj|kkk) ) continue;
				const double dx = cell_length * double(iii);
				const double dy = cell_length * double(jjj);
				const double dz = cell_length * double(kkk);
				const double dr2 = dx*dx + dy*dy + dz*dz;
				cf.eval_rcut(dr2 * (alpha*alpha));
				Slm<2*p, real_t> slm;
				slm.eval_opt(-dx, dy, dz); // eval S_l^{-m}
				// near cell correction
				const bool near = (abs(kkk)<=ICUT && abs(jjj)<=ICUT && abs(iii)<=ICUT);
				if(near){
					for(int l=0; l<=2*p; l++){
						cf.rcut[l] -= 1.0;
					}
				}
				for(int l=0; l<=2*p; l++){
					for(int m=0; m<=l; m++){
						const cplx_t val = cf.rcut[l] * slm.val_at(l, m);
						rsum.accum_at(l, m, val);
					}
				}
			}

			// wave-space sum
			Slm<2*p, real_t> ksum;
			ksum.clear();

			int mx, my, mz;
			for(mz=-MMAX; mz<=MMAX; mz++) for(my=-MMAX; my<=MMAX; my++) for(mx=-MMAX; mx<=MMAX; mx++)
			{
				if(0 == (mx|my|mz)) continue;
				const double dx = cell_length * double(i);
				const double dy = cell_length * double(j);
				const double dz = cell_length * double(k);

				const double theta = (+8.0 * atan(1.0)) * (dx*mx + dy*my + dz*mz);
				const cplx_t phase(cos(theta), sin(theta));

				Slm<2*p, real_t> slm;
				slm.eval_opt(-double(mx), double(my), double(mz));

				const double m2 = mx*mx + my*my + mz*mz;
				cf.eval_kcut(m2, alpha);

				for(int l=0; l<=2*p; l++){
					for(int m=0; m<=l; m++){
						const cplx_t val = (cf.kcut[l] * phase) * slm.val_at(l, m);
						ksum.accum_at(l, m, val);
					}
				}
			}
			// store sum
#ifdef IGNORE_RSPACE
			rsum.clear();
#endif
#ifdef IGNORE_KSPACE
			ksum.clear();
#endif
			for(int lm=0; lm<LEN2; lm++){
				gf_r[k][j][i][lm] = rsum.buf[lm] + ksum.buf[lm];
			}
#if 0
			if(0 == (i|j|k)){
				rsum.show();
				ksum.show();
			};
#endif
		} // for(k, j, i)
  }

};
