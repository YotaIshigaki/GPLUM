#pragma once

#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
#include "vector3.h"

#include "cell.h"
#include "GreenFunction_PBC.h"

template<int p, int NX, int NY, int NZ, int NMAX, int MMAX, int ICUT>
struct M2L_convolution_PBC
{
  enum{
    LEN  = lmbuf<p>::length,
  };

  typedef typename GreenFunction_PBC<p, NX, NY, NZ, NMAX, MMAX, ICUT>::real_t real_t;
  typedef typename GreenFunction_PBC<p, NX, NY, NZ, NMAX, MMAX, ICUT>::cplx_t cplx_t;

  std::vector<int> mm_list;
  // Multipole Moments
  cplx_t mm_k[NZ][NY][1+NX/2][LEN];
  // Local Expansions
  cplx_t le_k[NZ][NY][1+NX/2][LEN];
  // FFT buffer
  real_t rbuf[NZ][NY][NX];
  cplx_t kbuf[NZ][NY][1+NX/2];

  void set_mm_list(std::vector<int> _mm_list)
  {
    mm_list.resize(_mm_list.size());
    std::copy(_mm_list.begin(),_mm_list.end(),mm_list.begin());
  }

  void forward(MultipoleMoment<p> mm[NZ][NY][NX])
  {


	fftw_plan plan_fwd = 
		fftw_plan_dft_r2c_3d(
			NZ, NY, NX, 
			(double       *)(rbuf),
			(fftw_complex *)(kbuf),
			FFTW_ESTIMATE);

	int i, j, k;
	// forward multipole
	for(int li=0; li<mm_list.size(); li++){
	  int lm = mm_list[li];
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
		{
			rbuf[k][j][i] = mm[k][j][i].buf[lm];
		}

		fftw_execute(plan_fwd);

		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
		{
			mm_k[k][j][i][lm] = kbuf[k][j][i];
		}
	}

	fftw_destroy_plan(plan_fwd);
  }

  void transform(
	    const GreenFunction_PBC<p, NX, NY, NZ, NMAX, MMAX, ICUT> &gf
	    ){
	// M2L transformation
	typedef MultipoleMoment<p, cplx_t> (*mm_t)[NY][1+NX/2];
	typedef LocalExpansion <p, cplx_t> (*le_t)[NY][1+NX/2];
	gf.transform((mm_t)mm_k, (le_t)le_k);
  }

  void backward(LocalExpansion <p> le[NZ][NY][NX])
  {
	fftw_plan plan_bkw  =
		fftw_plan_dft_c2r_3d(
			NZ, NY, NX, 
			(fftw_complex *)(kbuf),
			(double       *)(rbuf),
			FFTW_ESTIMATE);
	int i, j, k;
	// backward local expansion
	for(int li=0; li<mm_list.size(); li++){
	  int lm = mm_list[li];
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
		{
			kbuf[k][j][i] = le_k[k][j][i][lm];
		}

		fftw_execute(plan_bkw);

		const double norm = 1.0 / (NX*NY*NZ);
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
		{
			le[k][j][i].buf[lm] = norm * rbuf[k][j][i];
		}

	}
	fftw_destroy_plan(plan_bkw);
}

};
