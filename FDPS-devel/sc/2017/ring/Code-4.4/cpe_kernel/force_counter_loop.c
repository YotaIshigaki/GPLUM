#include <slave.h>
#include <simd.h>

typedef struct Work{
	doublev4 dx;    //   0
	doublev4 dy;    //  32
	doublev4 dz;    //  64
	doublev4 rdotv; //  96
	doublev4 r2;
	doublev4 rinv;  // 160
} Work, *pWork;

typedef double REAL;
typedef struct{
  REAL r[4]; //x,y,z,m
  REAL v[4]; //vx,vy,vz,id
} EpiMM;

typedef EpiMM EpjMM;

typedef struct{
  REAL ax[4],ay[4],az[4],pot[4];
} ForceMM4;

typedef struct{
  REAL ax,ay,az,pot;
} ForceMM;

typedef struct{
  doublev4 m_r_coll_inv_qv4;
  doublev4 k_v_miv4;
  doublev4 eta_v_miv4;   
  doublev4 r_collv4;   
} Force_coff;

static inline void dv4_trans(doublev4* a, doublev4* b, doublev4* c, doublev4* d){
  const doublev4 a1 = simd_vshff(*b, *a, 0x44); // b1, b0, a1, a0 
  const doublev4 b1 = simd_vshff(*b, *a, 0xee); // b3, b2, a3, a2 
  const doublev4 a2 = simd_vshff(a1, a1, 0xd8); // b1, a1, b0, a0
  const doublev4 b2 = simd_vshff(b1, b1, 0xd8); // b3, a3, b2, a2

  const doublev4 c1 = simd_vshff(*d, *c, 0x44); // d1, d0, c1, c0 
  const doublev4 d1 = simd_vshff(*d, *c, 0xee); // d3, d2, c3, c2 
  const doublev4 c2 = simd_vshff(c1, c1, 0xd8); // d1, c1, d0, c0
  const doublev4 d2 = simd_vshff(d1, d1, 0xd8); // d3, c3, d2, c2
  
  *a = simd_vshff(c2, a2, 0x44); // d0, c0, b0, a0
  *b = simd_vshff(c2, a2, 0xee); // d1, c1, b1, a1
  *c = simd_vshff(d2, b2, 0x44); // d2, c2, b2, a2
  *d = simd_vshff(d2, b2, 0xee); // d3, c3, b3, a3
}

void force_counter_loop(
		const int         n_epi, 
        const EpiMM      *epi,     // double r[4], v[4]
		const EpjMM      *epj,     // double r[4], v[4]
		ForceMM          *force,   // double ax[4], ay[4], az[4], pot[4];
		const Force_coff *pcoff,   // doublev4 m_r_coll_inv_qv4, k_v_miv4, eta_v_miv4;   
		Work             *work     // doublev4 dx, dy, dz, rdotv, r2, rinv
		)
{
	int i, j;
	const doublev4 m_r_coll_inv_qv4 = pcoff->m_r_coll_inv_qv4;
	const doublev4 k_v_miv4         = pcoff->k_v_miv4; 
	const doublev4 eta_v_miv4       = pcoff->eta_v_miv4;   
	const doublev4 r_collv4         = pcoff->r_collv4;   

    //	const double mi[4] = {epi[0].r[3], epj[1].r[3], epj[2].r[3], epj[3].r[3],};
    const double zero[4]={0};
    for(j=0; j<4; j++){

      doublev4 fsx= *(doublev4 *)zero;
      doublev4 fsy= *(doublev4 *)zero;
      doublev4 fsz= *(doublev4 *)zero;
      doublev4 fsp= *(doublev4 *)zero;

      for(i=0; i<n_epi; i+=4){
            const doublev4 mi = *(doublev4 *)epi[i+3].r;
      
			// const double mj = epj[j].r[3];
          //			const double mj = mjs[j];

			const doublev4 dx    = work[i+j].dx;
			const doublev4 dy    = work[i+j].dy;
			const doublev4 dz    = work[i+j].dz;
			const doublev4 rdotv = work[i+j].rdotv;
			// const doublev4 r_sq  = work[i+j].r2;
			const doublev4 r_inv = work[i+j].rinv;

			const doublev4 r_inv_sq = r_inv * r_inv;
			const doublev4 xij      = 1.0 - r_collv4 * r_inv;
			const doublev4 rvji     = rdotv * r_inv_sq; 

			const doublev4 fsd = m_r_coll_inv_qv4 + k_v_miv4 * xij + eta_v_miv4 * rvji;
			const doublev4 pij = mi * r_inv;
			const doublev4 mri3 = - simd_vsellt(xij, fsd, r_inv_sq * pij);

			fsx -= mri3 * dx;
			fsy -= mri3 * dy;
			fsz -= mri3 * dz;
			fsp -= pij;
      }
      dv4_trans(&fsx, &fsy, &fsz, &fsp);
      *(doublev4 *)&force[j] += fsx + fsy + fsz + fsp;
	}
}
