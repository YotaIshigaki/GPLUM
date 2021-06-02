#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
#include "Common.h"

//#include "cell.h"
#include "GreenFunction_PBC.h"

struct M2L_convolution_PBC
{
  int mm_length;
  int num_cell_real;
  int num_cell_cplx;
  int cell_size[3];

  typedef  GreenFunction_PBC::real_t real_t;
  typedef  GreenFunction_PBC::cplx_t cplx_t;

  int m_begin;
  int m_end;
  int num_local_m;

  // Multipole Moments
  std::vector< MultipoleMoment<cplx_t> > mm_k;
  // Local Expansions
  std::vector< LocalExpansion<cplx_t> > le_k;
  // FFT buffer
  std::vector<real_t> rbuf;
  std::vector<cplx_t> kbuf;

  void set_size(int _cell_size[3], int p){
    cell_size[0] = _cell_size[0];
    cell_size[1] = _cell_size[1];
    cell_size[2] = _cell_size[2];
    mm_length = (p+1)*(p+1);
    num_cell_real = cell_size[2]*cell_size[1]*cell_size[0];
    num_cell_cplx = cell_size[2]*cell_size[1]*(1+cell_size[0]/2);
    mm_k.resize(num_cell_cplx,p);
    //    for(int m=0;m<num_cell_cplx;m++) mm_k[m].set_order(p);
    le_k.resize(num_cell_cplx,p);
    //    for(int m=0;m<num_cell_cplx;m++) le_k[m].set_order(p);
    rbuf.resize(num_cell_real);
    kbuf.resize(num_cell_cplx);
  }

  void set_mm_range(int begin, int end)
  {
    m_begin = begin;
    m_end = end;
    num_local_m = m_end-m_begin;
  }

  void forward(std::vector< std::vector<real_t> > &mm)
  {
    static bool initcall = true;
    static fftw_plan plan_fwd;
    if(initcall){
      initcall = false;
      plan_fwd = fftw_plan_dft_r2c_3d(
			   cell_size[2], cell_size[1], cell_size[0],
			   (double       *)(&(rbuf[0])),
			   (fftw_complex *)(&(kbuf[0])),
			   FFTW_ESTIMATE);
    }
    int i;
    // forward multipole
    for(int lm=0; lm<num_local_m; lm++){
      for(i=0;i<num_cell_real;i++){
	rbuf[i] = mm[lm][i];
      }
      fftw_execute(plan_fwd);
      for(i=0;i<num_cell_cplx;i++){
	mm_k[i].buf[m_begin+lm] = kbuf[i];
      }
    }
    //    fftw_destroy_plan(plan_fwd);
  }

  void transform(
		 const GreenFunction_PBC &gf
		 ){
	// M2L transformation
	typedef MultipoleMoment<cplx_t> *mm_t;
	typedef LocalExpansion <cplx_t> *le_t;
	gf.transform((mm_t)(&(mm_k[0])), (le_t)(&(le_k[0])),
		     m_begin, m_end);
  }

  void backward(std::vector< std::vector<real_t> > &le)
  {
    static bool initcall = true;
    static fftw_plan plan_bkw;
    if(initcall){
      initcall = false;
      plan_bkw = fftw_plan_dft_c2r_3d(
			 cell_size[2], cell_size[1], cell_size[0],
			 (fftw_complex *)(&(kbuf[0])),
			 (double       *)(&(rbuf[0])),
			 FFTW_ESTIMATE);
    }
    int i;
    // backward local expansion
    for(int lm=0; lm<num_local_m; lm++){
      for(i=0; i<num_cell_cplx; i++)
	{
	  kbuf[i] = le_k[i].buf[m_begin+lm];
	}
      fftw_execute(plan_bkw);
      const double norm = 1.0 / (num_cell_real);
      for(i=0; i<num_cell_real; i++)
	{
	  le[lm][i] = norm * rbuf[i];
	}
    }
    //    fftw_destroy_plan(plan_bkw);
  }
  
};
