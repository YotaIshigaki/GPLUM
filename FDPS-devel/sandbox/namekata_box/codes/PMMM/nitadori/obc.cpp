#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
#include "vector3.h"

#include "fmm.h"

typedef vector3<int> ivec3;

inline ivec3 cell_nearest(const dvec3 &pos, const double d){
	return ivec3(pos / d);
}
inline dvec3 cell_pos(const ivec3 &idx, const double d){
	return d * (dvec3(0.5) + dvec3(idx));
}

#if 0
static double wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + 1.e-6 * (double)tv.tv_usec;
}
#endif

#include "cell.h"

template<int p, int NX, int NY, int NZ>
struct GreenFunction_OBC{
	typedef double               real_t;
	// typedef std::complex<real_t> cplx_t;
	typedef double _Complex      cplx_t;
	enum{
		LEN  = lmbuf<p>::length,
		LEN2 = lmbuf<2*p>::length,
	};

	real_t gf_r[2*NZ][2*NY][2*NX][LEN2];
	cplx_t gf_k[2*NZ][2*NY][1+NX][LEN2];

	void gen_gf_r(const int icut, const double cell_length){
		int i, j, k;
		for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<2*NX; i++)
		{
			if(k==NZ || j==NY || i==NX){
				for(int lm=0; lm<LEN2; lm++){
					gf_r[k][j][i] [lm] = 0.0;
				}
				continue;
			}
			const int kk = (k>NZ) ? k - 2*NZ : k;
			const int jj = (j>NY) ? j - 2*NY : j;
			const int ii = (i>NX) ? i - 2*NX : i;
			if(abs(kk)<=icut && abs(jj)<=icut && abs(ii)<=icut){
				for(int lm=0; lm<LEN2; lm++){
					gf_r[k][j][i] [lm] = 0.0;
				}
				continue;
			}
			const double dx = cell_length * double(ii);
			const double dy = cell_length * double(jj);
			const double dz = cell_length * double(kk);
			Slm<2*p, real_t> slm;
			slm.eval_opt(-dx, dy, dz); // eval S_l^{-m}
			for(int lm=0; lm<LEN2; lm++){
				gf_r[k][j][i] [lm] = slm.buf[lm];
			}
		}
#if 0
		for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<2*NX; i++)
		{
			printf("(%d, %d, %d) : %e\n", i, j, k, gf_r[k][j][i][0]);
		}
#endif
	}

	void gen_gf_k(){
		static real_t rbuf[2*NZ][2*NY][2*NX];
		static cplx_t kbuf[2*NZ][2*NY][1+NX];
		fftw_plan plan_fwd = 
			fftw_plan_dft_r2c_3d(
				2*NZ, 2*NY, 2*NX, 
				(double       *)(rbuf),
				(fftw_complex *)(kbuf),
				FFTW_ESTIMATE);
		for(int lm=0; lm<LEN2; lm++){
			int i, j, k;
			for(k=0; k<2*NZ; k++) for(int j=0; j<2*NY; j++) for(i=0; i<2*NX; i++)
			{
				rbuf[k][j][i] = gf_r[k][j][i] [lm];
			}
			// CALL FFTW
			fftw_execute(plan_fwd);

			for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<1+NX; i++)
			{
				gf_k[k][j][i] [lm] = kbuf[k][j][i];
			}
		}
	}

	typedef MultipoleMoment <p, cplx_t> mm_t;
	typedef LocalExpansion  <p, cplx_t> le_t;

	void transform(
			const mm_t mm_k[2*NZ][2*NY][1+NX],
			      le_t le_k[2*NZ][2*NY][1+NX]) const
	{
		for(int k=0; k<2*NZ; k++){
			for(int j=0; j<2*NY; j++){
				for(int i=0; i<1+NX; i++){
					typedef Slm<2*p, cplx_t> slm_t;
					((slm_t *)(gf_k[k][j][i]))
						-> template transform_M2L<p, p, false>(
								mm_k[k][j][i], le_k[k][j][i]);
				}
			}
		}
	}
};

static void PP_interact_inner(Cell &ci, const Cell &cj){
	const int ni = ci.plist.size();
	const int nj = cj.plist.size();
	for(int i=0; i<ni; i++){
		Particle &pi = *ci.plist[i];
		for(int j=0; j<nj; j++){
			const Particle &pj = *cj.plist[j];
			if(&pi == &pj) continue;

			const dvec3  dr  = pj.pos - pi.pos;
			const double r2  = dr*dr;
			const double ri2 = 1.0 / r2;
			const double ri  = sqrt(ri2);
			const double ri3 = ri * ri2;
			pi.phi_direct += pj.mass * ri;
			pi.acc_direct += (pj.mass * ri3) * dr;
		}
	}
}

template <int PFMM, int ICUT, int NX, int NY, int NZ>
void PP_interact_OBC(Cell_FMM<PFMM> cell[NZ][NY][NX])
{
	int i, j, k;
	for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
	{
		int ii, jj, kk;
		for(kk=k-ICUT; kk<=k+ICUT; kk++) for(jj=j-ICUT; jj<=j+ICUT; jj++) for(ii=i-ICUT; ii<=i+ICUT; ii++)
		{
			if(kk < 0 || kk >= NZ) continue;
			if(jj < 0 || jj >= NY) continue;
			if(ii < 0 || ii >= NX) continue;
			PP_interact_inner(cell[k][j][i], cell[kk][jj][ii]);
		}
	}
}

template<int p, int NX, int NY, int NZ>
void M2L_convolution_OBC(
		const GreenFunction_OBC<p, NX, NY, NZ> &gf,
		Cell_FMM<p> cell[NZ][NY][NX])
{
	typedef typename GreenFunction_OBC<p, NX, NY, NZ>::real_t real_t;
	typedef typename GreenFunction_OBC<p, NX, NY, NZ>::cplx_t cplx_t;
	enum{
		LEN  = lmbuf<p>::length,
	};

	// Multipole Moments
	// static real_t mm_r[2*NZ][2*NY][2*NX][LEN];
	static cplx_t mm_k[2*NZ][2*NY][1+NX][LEN];
	// Local Expansions
	// static real_t le_r[2*NZ][2*NY][2*NX][LEN];
	static cplx_t le_k[2*NZ][2*NY][1+NX][LEN];
	// FFT buffer
	static real_t rbuf[2*NZ][2*NY][2*NX];
	static cplx_t kbuf[2*NZ][2*NY][1+NX];

	fftw_plan plan_fwd = 
		fftw_plan_dft_r2c_3d(
			2*NZ, 2*NY, 2*NX, 
			(double       *)(rbuf),
			(fftw_complex *)(kbuf),
			FFTW_ESTIMATE);
	fftw_plan plan_bkw  =
		fftw_plan_dft_c2r_3d(
			2*NZ, 2*NY, 2*NX, 
			(fftw_complex *)(kbuf),
			(double       *)(rbuf),
			FFTW_ESTIMATE);

	int i, j, k;
	// clear rbuf
	for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<2*NX; i++)
	{
		rbuf[k][j][i] = 0.0;
	}
	// forward multipole
	for(int lm=0; lm<LEN; lm++){
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
		{
			rbuf[k][j][i] = cell[k][j][i].mm.buf[lm];
		}

		fftw_execute(plan_fwd);

		for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<1+NX; i++)
		{
			mm_k[k][j][i][lm] = kbuf[k][j][i];
		}
	}
	// M2L transformation
	typedef MultipoleMoment<p, cplx_t> (*mmarray)[2*NY][1+NX];
	typedef LocalExpansion <p, cplx_t> (*learray)[2*NY][1+NX];
	gf.transform((mmarray)mm_k, (learray)le_k);

	// backward local expansion
	for(int lm=0; lm<LEN; lm++){
		for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<1+NX; i++)
		{
			kbuf[k][j][i] = le_k[k][j][i][lm];
		}

		fftw_execute(plan_bkw);

		const double norm = 1.0 / (8*NX*NY*NZ);
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
		{
			cell[k][j][i].le.buf[lm] = norm * rbuf[k][j][i];
		}

	}
	fftw_destroy_plan(plan_fwd);
	fftw_destroy_plan(plan_bkw);
}

static void print_err(
		const std::vector<double> &err,
		const char * name,
		const int icut,
		const int p)
{
	static char fname[256];
	sprintf(fname, "%s.c%dp%d.dat", name, icut, p);
	FILE *fp = fopen(fname, "w");
	assert(fp);
	const int len = err.size();
	for(int i=0; i<len; i++){
		fprintf(fp, "%e %e\n", double(i)/len, err[i]);
	}
	fclose(fp);
}

int main(){
	enum{
		NP   = 128,
		NC   = 8,
		NC3  = NC*NC*NC,
		PFMM = 7,
		ICUT = 2,
	};
	typedef Cell_FMM<PFMM>	Cell_t;

	static Particle ptcl[NP];
	static Cell_t   cell[NC][NC][NC];
	const double clen = 1.0 / NC;

	for(int k=0; k<NC; k++){
		for(int j=0; j<NC; j++){
			for(int i=0; i<NC; i++){
				cell[k][j][i].set(ivec3(i,j,k), clen);
			}
		}
	}

	Particle::gen_rand_dist(NP, ptcl);
	double msum = 0.0;
	for(int i=0; i<NP; i++){
		msum += ptcl[i].mass;
	}
	// printf("msum = %e\n", msum);
	for(int i=0; i<NP; i++){
		const ivec3 idx = cell_nearest(ptcl[i].pos, clen);
		assert(0 <= idx.x && idx.x < NC);
		assert(0 <= idx.y && idx.y < NC);
		assert(0 <= idx.z && idx.z < NC);
		cell[idx.z][idx.y][idx.x].plist.push_back(&ptcl[i]);
	}

	for(int k=0; k<NC; k++) for(int j=0; j<NC; j++) for(int i=0; i<NC; i++){
		cell[k][j][i].sanity_check();
	}

	puts("Eval PP");
	PP_interact_OBC<PFMM, ICUT, NC, NC, NC>  (cell);

	puts("Gen Green");
	static GreenFunction_OBC<PFMM, NC, NC, NC> gf;
	gf.gen_gf_r(ICUT, 1./NC);
	gf.gen_gf_k();

	puts("Eval PM");
	Cell_t *cell1d = cell[0][0];
	for(int i=0; i<NC3; i++){
		cell1d[i].do_P2M();
	}

	M2L_convolution_OBC<PFMM, NC, NC, NC> (gf, cell);

	for(int i=0; i<NC3; i++){
		cell1d[i].do_L2P();
	}

	dvec3 fpp(0.0), fpm(0.0);
	for(int i=0; i<NP; i++){
		fpp += ptcl[i].mass * ptcl[i].acc_direct;
		fpm += ptcl[i].mass * ptcl[i].acc_app;
	}
	printf("PP ftot : (%e, %e, %e)\n", fpp.x, fpp.y, fpp.z);
	printf("PM ftot : (%e, %e, %e)\n", fpm.x, fpm.y, fpm.z);

	for(int i=0; i<NP; i++){
		ptcl[i].move_accp();
	}

	puts("eval dirct force");
#pragma omp parallel for
	for(int i=0; i<NP; i++){
		Particle &pi = ptcl[i];
		for(int j=0; j<NP; j++){
			const Particle pj = ptcl[j];
			if(j == i) continue;
			const dvec3  dr  = pj.pos - pi.pos;
			const double r2  = dr*dr;
			const double ri2 = 1.0 / r2;
			const double ri  = sqrt(ri2);
			const double ri3 = ri * ri2;
			pi.phi_direct += pj.mass * ri;
			pi.acc_direct += (pj.mass * ri3) * dr;
		}
	}

#if 1
	std::vector<double> err(NP);
	for(int i=0; i<NP; i++) err[i] = ptcl[i].adiff_rel();
	std::sort(err.begin(), err.end());
	print_err(err, "adiffr", ICUT, PFMM);

	for(int i=0; i<NP; i++) err[i] = ptcl[i].pdiff_rel();
	std::sort(err.begin(), err.end());
	print_err(err, "pdiffr", ICUT, PFMM);
#endif

	return 0;
}
