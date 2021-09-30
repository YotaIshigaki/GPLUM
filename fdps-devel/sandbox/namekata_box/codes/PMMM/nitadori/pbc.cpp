// PMMM for periodic boundary condition
#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
#include "vector3.h"

#define EWALD_FILE
// #define NON_CHARGE_NEUTRAL

typedef vector3<int> ivec3;

inline ivec3 cell_nearest(const dvec3 &pos, const double d){
	return ivec3(pos / d);
}
inline dvec3 cell_pos(const ivec3 &idx, const double d){
	return d * (dvec3(0.5) + dvec3(idx));
}

inline dvec3 minimum_image(const dvec3 &inp){
	return dvec3(
			inp.x - round(inp.x),
			inp.y - round(inp.y),
			inp.z - round(inp.z));
}

#include "cell.h"
#include "cutoff.h"
#include "ewald.h"

#if 0
static double wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + 1.e-6 * (double)tv.tv_usec;
}
#endif

template<int p, int NX, int NY, int NZ>
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

	template <int NMAX, int MMAX, int ICUT>
	void gen_gf_r(const double alpha, const double cell_length){
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
};

static void PP_interact_inner(Cell &ci, const Cell &cj){
	const int ni = ci.plist.size();
	const int nj = cj.plist.size();
	for(int i=0; i<ni; i++){
		Particle &pi = *ci.plist[i];
		for(int j=0; j<nj; j++){
			const Particle &pj = *cj.plist[j];
			if(&pi == &pj) continue;

			const dvec3  dr  = minimum_image(pj.pos - pi.pos);
			const double r2  = dr*dr;
			const double ri2 = 1.0 / r2;
			const double ri  = sqrt(ri2);
			const double ri3 = ri * ri2;
			pi.phi_direct += pj.mass * ri;
			pi.acc_direct += (pj.mass * ri3) * dr;
			// puts("evaluated PP");
		}
	}
}

template <int PFMM, int ICUT, int NX, int NY, int NZ>
void PP_interact_PBC(Cell_FMM<PFMM> cell[NZ][NY][NX])
{
	int i, j, k;
	for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
	{
		int ii, jj, kk;
		for(kk=k-ICUT; kk<=k+ICUT; kk++) for(jj=j-ICUT; jj<=j+ICUT; jj++) for(ii=i-ICUT; ii<=i+ICUT; ii++)
		{
			const int iii = (ii+NX)%NX;
			const int jjj = (jj+NY)%NY;
			const int kkk = (kk+NX)%NZ;
			PP_interact_inner(cell[k][j][i], cell[kkk][jjj][iii]);
		}
	}
}

template<int p, int NX, int NY, int NZ>
void M2L_convolution_PBC(
		const GreenFunction_PBC<p, NX, NY, NZ> &gf,
		Cell_FMM<p> cell[NZ][NY][NX])
{
	typedef typename GreenFunction_PBC<p, NX, NY, NZ>::real_t real_t;
	typedef typename GreenFunction_PBC<p, NX, NY, NZ>::cplx_t cplx_t;
	enum{
		LEN  = lmbuf<p>::length,
	};

	// Multipole Moments
	static cplx_t mm_k[NZ][NY][1+NX/2][LEN];
	// Local Expansions
	static cplx_t le_k[NZ][NY][1+NX/2][LEN];
	// FFT buffer
	static real_t rbuf[NZ][NY][NX];
	static cplx_t kbuf[NZ][NY][1+NX/2];

	fftw_plan plan_fwd = 
		fftw_plan_dft_r2c_3d(
			NZ, NY, NX, 
			(double       *)(rbuf),
			(fftw_complex *)(kbuf),
			FFTW_ESTIMATE);
	fftw_plan plan_bkw  =
		fftw_plan_dft_c2r_3d(
			NZ, NY, NX, 
			(fftw_complex *)(kbuf),
			(double       *)(rbuf),
			FFTW_ESTIMATE);

	int i, j, k;
	// forward multipole
	for(int lm=0; lm<LEN; lm++){
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
		{
			rbuf[k][j][i] = cell[k][j][i].mm.buf[lm];
		}

		fftw_execute(plan_fwd);

		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
		{
			mm_k[k][j][i][lm] = kbuf[k][j][i];
		}
	}
	// M2L transformation
	typedef MultipoleMoment<p, cplx_t> (*mmarray_t)[NY][1+NX/2];
	typedef LocalExpansion <p, cplx_t> (*learray_t)[NY][1+NX/2];
	gf.transform((mmarray_t)mm_k, (learray_t)le_k);

	// backward local expansion
	for(int lm=0; lm<LEN; lm++){
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
		{
			kbuf[k][j][i] = le_k[k][j][i][lm];
		}

		fftw_execute(plan_bkw);

		const double norm = 1.0 / (NX*NY*NZ);
		for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
		{
			cell[k][j][i].le.buf[lm] = norm * rbuf[k][j][i];
		}

	}
	fftw_destroy_plan(plan_fwd);
	fftw_destroy_plan(plan_bkw);
}

#if 1
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
#endif

int main(){
	enum{
		NP   = 128,
		NC   = 8,
		NC3  = NC*NC*NC,
		PFMM = 5,
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
#ifdef EWALD_FILE
	{
		FILE *fp = fopen("qpos.dat", "r");
		assert(fp);
		int n;
		fscanf(fp, "%d", &n);
		assert(NP == n);
		for(int i=0; i<NP; i++){
			fscanf(fp, "%lf %lf %lf %lf",
					&ptcl[i].mass,
					&ptcl[i].pos.x,
					&ptcl[i].pos.y,
					&ptcl[i].pos.z);
		}
		fclose(fp);
	}
#else
	Particle::gen_rand_dist(NP, ptcl);
	double msum = 0.0;
	for(int i=0; i<NP; i++){
		msum += ptcl[i].mass;
	}
#endif
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
	PP_interact_PBC<PFMM, ICUT, NC, NC, NC>  (cell);

	puts("Gen Green");
	static GreenFunction_PBC<PFMM, NC, NC, NC> gf;
#if 1
	const double alpha = 2.4;
	gf.gen_gf_r<3, 5, ICUT>(alpha, 1./NC);
#else
	const double alpha = 1.5;
	gf.gen_gf_r<4, 3, ICUT>(alpha, 1./NC);
#endif
	gf.gen_gf_k();

	puts("Eval PM");
	Cell_t *cell1d = cell[0][0];
	for(int i=0; i<NC3; i++){
		cell1d[i].do_P2M();
	}

	M2L_convolution_PBC<PFMM, NC, NC, NC> (gf, cell);

	// Dipole correction
#if 1
	dvec3 dipole(0.0);
	double quad0 = 0.0;
	for(int i=0; i<NC3; i++){
		dipole.x += cell1d[i].mm.buf[3];
		dipole.y += cell1d[i].mm.buf[1];
		dipole.z += cell1d[i].mm.buf[2];
		quad0 += cell1d[i].dispersion();
#endif
	}
	const double pi = 4.0 * atan(1.0);
	dipole *= (4./3.) * pi;
	printf("quad : %e\n", quad0);
	for(int i=0; i<NC3; i++){
		cell1d[i].le.buf[3] += 2.0 * dipole.x;
		cell1d[i].le.buf[1] -= 2.0 * dipole.y;
		cell1d[i].le.buf[2] += 1.0 * dipole.z;
		cell1d[i].le.buf[0] += ((2./3.) * pi) * quad0;
		// self energy correction
		cell1d[i].le.buf[0] -= 
			alpha * (2.0/sqrt(pi)) * cell1d[i].mm.buf[0];
	}

	for(int i=0; i<NC3; i++){
		cell1d[i].do_L2P();
#ifdef NON_CHARGE_NEUTRAL
		cell1d[i].do_L2P_corr(msum, alpha);
#endif
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
#if 0
	for(int i=0; i<NP; i++){
		const dvec3 acc = ptcl[i].acc_app;
		printf("%24.16e, (%24.16e, %24.16e, %24.16e)\n", 
				ptcl[i].phi_app, acc.x, acc.y, acc.z);
	}
#endif

	puts("Ewald sum");
	// Direct Ewald
#ifdef EWALD_FILE
	{
		FILE *fp = fopen("ewald.dat", "r");
		assert(fp);
		int n;
		fscanf(fp, "%d", &n);
		assert(NP == n);
		for(int i=0; i<NP; i++){
			fscanf(fp, "%lf %lf %lf %lf",
					&ptcl[i].phi_direct,
					&ptcl[i].acc_direct.x,
					&ptcl[i].acc_direct.y,
					&ptcl[i].acc_direct.z);
		}
		fclose(fp);
	}
#else
	const double alpha_ewald = 2.4;
	eval_k_space<5>(NP, alpha_ewald, ptcl);
	eval_r_space<3>(NP, alpha_ewald, msum, ptcl);
#endif
#if 0
	for(int i=0; i<NP; i++){
		const dvec3 acc = ptcl[i].acc_direct;
		printf("%24.16e, (%24.16e, %24.16e, %24.16e)\n", 
				ptcl[i].phi_direct, acc.x, acc.y, acc.z);
	}
#endif
	double en_app=0.0, en_dir=0.0;
	for(int i=0; i<NP; i++){
		en_app += 0.5 * ptcl[i].mass * ptcl[i].phi_app;
		en_dir += 0.5 * ptcl[i].mass * ptcl[i].phi_direct;
	}
	printf("energy : %24.16e, %24.16e\n", en_app, en_dir);

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

