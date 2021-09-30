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

#include "pmmm_fmm.h"
#include "cutoff.h"

struct GreenFunction_PBC{
  typedef double               real_t;
  typedef std::complex<real_t> cplx_t;
  // typedef double _Complex      cplx_t;
  int order;
  int length;
  int length2;
  int num_cell_real;
  int num_cell_cplx;
  int cell_size[3];

  typedef std::vector<real_t> gf_r_t;
  typedef lmbuf<cplx_t> gf_k_t;
  std::vector<gf_r_t> gf_r;
  std::vector<gf_k_t> gf_k;

  void set_size(int _cell_size[3], int p){
    order = p;
    cell_size[0] = _cell_size[0];
    cell_size[1] = _cell_size[1];
    cell_size[2] = _cell_size[2];
    num_cell_real = cell_size[2]*cell_size[1]*cell_size[0];
    num_cell_cplx = cell_size[2]*cell_size[1]*(1+cell_size[0]/2);
    length = (p+1)*(p+1);
    length2 = (2*p+1)*(2*p+1);
    gf_r.resize(num_cell_real);
    for(int m=0;m<num_cell_real;m++)gf_r[m].resize(length2);
    gf_k.resize(num_cell_cplx);
    for(int m=0;m<num_cell_cplx;m++)gf_k[m].set_order(2*p);
  }

  void gen_gf_k(){
    std::vector<real_t> rbuf(num_cell_real);
    std::vector<cplx_t> kbuf(num_cell_cplx);
    
    fftw_plan plan_fwd = 
      fftw_plan_dft_r2c_3d(
			   cell_size[2], cell_size[1], cell_size[0],
			   (double       *)(&(rbuf[0])),
			   (fftw_complex *)(&(kbuf[0])),
			   FFTW_ESTIMATE);
    for(int lm=0; lm<length2; lm++){
      int i;
      for(i=0; i<num_cell_real; i++)
	{
	  rbuf[i] = gf_r[i][lm];
	}
      // CALL FFTW
      fftw_execute(plan_fwd);

      for(i=0; i<num_cell_cplx; i++)
	{
	  gf_k[i].buf[lm] = kbuf[i];
	}
    }
    fftw_destroy_plan(plan_fwd);
  }

  typedef MultipoleMoment <cplx_t> *mm_t;
  typedef LocalExpansion  <cplx_t> *le_t;

  void transform(
		 const mm_t mm_k,
		 le_t le_k, int m_begin, int m_end) const
  {
    int i;
#pragma omp parallel for
    for(i=0; i<num_cell_cplx; i++)
      {
	//	printf("transform %d\n",i);
	typedef Slm<cplx_t> slm_t;
	((slm_t *)(&(gf_k[i])))
	  ->  transform_M2L<false>(
				   mm_k[i], le_k[i], m_begin, m_end);
      }
  }

  //  template <int NMAX, int MMAX, int ICUT>
  void gen_gf_r(const double alpha, const double cell_length, 
		const int nmax, const int mmax, const int icut)
  {
    {
      Slm<real_t> dry(2*order);
      dry.init();
    }
#pragma omp parallel for
    for(int k=0; k<cell_size[2]; k++) 
      for(int j=0; j<cell_size[1]; j++) 
	for(int i=0; i<cell_size[0]; i++)
	  {
	    CutFunc cf(2*order);
	    // real-space sum
	    Slm<real_t> rsum(2*order);
	    rsum.clear();

	    const int kk = (k>=cell_size[2]/2) ? k - cell_size[2] : k;
	    const int jj = (j>=cell_size[1]/2) ? j - cell_size[1] : j;
	    const int ii = (i>=cell_size[0]/2) ? i - cell_size[0] : i;
	    int nx, ny, nz;
	    for(nz=-nmax; nz<=nmax; nz++) 
	      for(ny=-nmax; ny<=nmax; ny++) 
		for(nx=-nmax; nx<=nmax; nx++)
		  {
		    const int kkk = kk + nz*cell_size[2];
		    const int jjj = jj + ny*cell_size[1];
		    const int iii = ii + nx*cell_size[0];
		    if( 0 == (iii|jjj|kkk) ) continue;
		    const double dx = cell_length * double(iii);
		    const double dy = cell_length * double(jjj);
		    const double dz = cell_length * double(kkk);
		    const double dr2 = dx*dx + dy*dy + dz*dz;
		    cf.eval_rcut(dr2 * (alpha*alpha));
		    Slm<real_t> slm(2*order);
		    slm.eval_opt(-dx, dy, dz); // eval S_l^{-m}
		    // near cell correction
		    const bool near = (abs(kkk)<=icut && abs(jjj)<=icut && abs(iii)<=icut);
		    if(near){
		      for(int l=0; l<=2*order; l++){
			cf.rcut[l] -= 1.0;
		      }
		    }
		    for(int l=0; l<=2*order; l++){
		      for(int m=0; m<=l; m++){
			const cplx_t val = cf.rcut[l] * slm.val_at(l, m);
			rsum.accum_at(l, m, val);
		      }
		    }
		  }
	    
	// wave-space sum
	    Slm<real_t> ksum(2*order);
	    ksum.clear();

	    int mx, my, mz;
	    for(mz=-mmax; mz<=mmax; mz++)
	      for(my=-mmax; my<=mmax; my++)
		for(mx=-mmax; mx<=mmax; mx++)
		  {
		    if(0 == (mx|my|mz)) continue;
		    const double dx = cell_length * double(i);
		    const double dy = cell_length * double(j);
		    const double dz = cell_length * double(k);
		    
		    const double theta = (+8.0 * atan(1.0)) * (dx*mx + dy*my + dz*mz);
		    const cplx_t phase(cos(theta), sin(theta));
		
		    Slm<real_t> slm(2*order);
		    slm.eval_opt(-double(mx), double(my), double(mz));

		    const double m2 = mx*mx + my*my + mz*mz;
		    cf.eval_kcut(m2, alpha);

		    for(int l=0; l<=2*order; l++){
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
	    for(int lm=0; lm<length2; lm++){
	      gf_r[i+cell_size[0]*(j+cell_size[1]*k)][lm] = rsum.buf[lm] + ksum.buf[lm];
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
