#include"slave.h"
#include"simd.h"
#include"dma.h"

//#include<math.h>
//#include<stdio.h>
//#include<stdlib.h>
//#include<string.h>
#include"cpe_prof.h"
#include"cpe_func.h"

void *ldm_malloc(size_t);
void ldm_free(void *, size_t);
//void *memset(void *b, int c, size_t len);
//void *memcpy(void *b, void* c, size_t len);
void abort();

#define ID_CHECK 63
#define REP4(A) A,A,A,A

#define HUGE 3.40282347e+38F

typedef struct{
  doublev4 m_r_coll_inv_qv4; // mj/r_coll^3
  doublev4 k_v_miv4;         // kappa/mi
  doublev4 eta_v_miv4;       // eta/mi
  doublev4 r_collv4;   
} Force_coff;

#ifdef SPLIT
typedef struct Work{
	doublev4 dx;
	doublev4 dy;
	doublev4 dz;
	doublev4 rdotv;
	doublev4 r2;
	doublev4 rinv;
} Work, *pWork;
#endif

//#define simd_vseleqw(__A,__B,__C) (intv8)__builtin_sw64_vseleqw(__A,__B,__C)
//__thread_local dma_desc epj_input;
//__thread_local dma_desc pj_index;
//__thread_local dma_desc spj_input;
//__thread_local dma_desc sat_input;
//__thread_local dma_desc force_sat_input;
//__thread_local dma_desc force_sat_output;
// 
//__thread_local volatile unsigned int par_reply, epi_reply;
//__thread_local volatile unsigned int epj_reply, spj_reply;
//__thread_local volatile unsigned int index_reply;
//__thread_local volatile unsigned int sat_reply, fin_reply, fout_reply;
//__thread_local Force_Kernel_Pars pars;
//__thread_local pWork work;
//__thread_local ForceMM* force_sat;
//__thread_local EpiMM*   sat_lm;


__thread_local unsigned int my_id, my_row, my_col;
//__thread_local const intv8 ilowmask={0,0xffffffff,0,0xffffffff,0,0xffffffff,0,0xffffffff};
//__thread_local const intv8 ihighmask={0xffffffff,0,0xffffffff,0,0xffffffff,0,0xffffffff,0};
//__thread_local const intv8 ilowone ={0,1,0,1,0,1,0,1};
//__thread_local const intv8 i0v8={0,0,0,0,0,0,0,0};
//__thread_local const doublev4 twov4 = {2.0,2.0,2.0,2.0};
__thread_local volatile doublev4 d1v4=1.0;
__thread_local volatile doublev4 d0v4=0.0;
__thread_local volatile doublev4 dhv4=0.5;
__thread_local volatile doublev4 d8v4=8.0;
//__thread_local doublev4 dlowone;
//#if !defined(SPLIT) && !defined(FULL_ASM)
__thread_local volatile doublev4 dhuge=HUGE;
#if defined(SPLIT)
__thread_local const uint256 glb_tiny = { REP4(1L << 55) };
#endif


#ifdef PROFILE
__thread_local unsigned long p_count, sp_count, t_count, tr_count, kd_count, init_count, mem_count, jlst_count, dmaj_count, dmas_count, geti_count, dmaw_count, st_count, pl_count, ss_count,gets_count, puti_count, puts_count, rfst_count, t2_count, final_sync_count, final_max_count, final_min_count, final_tot_count;
#ifdef SUNWAY_FORCE_KERNEL_CHECK
__thread_local unsigned long ch_count;
#endif
#endif

/*static inline void rsqrtd(const doublev4 *x, doublev4 *yout){
	doublev4 xh = 0.5 * *x;
	uint256 ivec = *(uint256 *)x;
	// ivec = 0x5fe6eb50c7b537a9 - (ivec >> 1);
	uint256 sr = ivec >> 1;
	uint256 magic = 0x5fe6eb50c7b537a9;
	uint256 mask  = 0x7fffffffffffffff;
	ivec = simd_vsubl(magic, sr&mask);

	doublev4 y = *(doublev4 *)&ivec;

	y += y * (0.5 - xh*(y*y));
	y += y * (0.5 - xh*(y*y));
    y += y * (0.5 - xh*(y*y));
    y += y * (0.5 - xh*(y*y));

	*yout = y;
	return;
    }*/


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

static inline void epi_trans(const int ni, EpiMM* epi) {
  int i;
  for(i=0; i<ni; i+=4) {
    doublev4 p0,p1,p2,p3;

    // x,y,z,id
    simd_load(p0, (double*)epi[i  ].r);
    simd_load(p1, (double*)epi[i+1].r);
    simd_load(p2, (double*)epi[i+2].r);
    simd_load(p3, (double*)epi[i+3].r);
    dv4_trans(&p0,&p1,&p2,&p3);
    simd_store(p0, (double*)epi[i  ].r);
    simd_store(p1, (double*)epi[i+1].r);
    simd_store(p2, (double*)epi[i+2].r);
    simd_store(p3, (double*)epi[i+3].r);
    
    // vx,vy,vz,pot
    simd_load(p0, (double*)epi[i  ].v);
    simd_load(p1, (double*)epi[i+1].v);
    simd_load(p2, (double*)epi[i+2].v);
    simd_load(p3, (double*)epi[i+3].v);
    dv4_trans(&p0,&p1,&p2,&p3);
    simd_store(p0, (double*)epi[i  ].v);
    simd_store(p1, (double*)epi[i+1].v);
    simd_store(p2, (double*)epi[i+2].v);
    simd_store(p3, (double*)epi[i+3].v);
  }
}

static inline void reduce(doublev4 *d) {
  doublev4 t;
  unsigned int dest;
  if((my_row&1)==0) {
    dest = my_row+1;
    REG_PUTC(*d,dest);
  }
  else {
    REG_GETC(t);  
    *d += t;

    if((my_row&2)==0) {
      dest = my_row+2;
      REG_PUTC(*d,dest);
    }
    else {
      REG_GETC(t);  
      *d += t;
      
      if(my_row==3) {
        dest = 7;
        REG_PUTC(*d,dest);
      }
      else {
        REG_GETC(t);  
        *d += t;

        if((my_col&1)==0) {
          dest = my_col+1;
          REG_PUTR(*d,dest);
        }
        else {
          REG_GETR(t);  
          *d += t;

          if((my_col&2)==0) {
            dest = my_col+2;
            REG_PUTR(*d,dest);
          }
          else {
            REG_GETR(t);  
            *d += t;

            if(my_col==3) {
               dest = 7;
               REG_PUTR(*d,dest);
            }
            else {
              REG_GETR(t);  
              *d += t;
            }
          }
        }
      }
    }
  }
}

static inline void reduce_ulong(unsigned long *input) {
  doublev4 dd,tt;
  unsigned long *d=(unsigned long*)&dd;
  unsigned long *t=(unsigned long*)&tt;
  *d = *input;
  unsigned int dest;
  if((my_row&1)==0) {
    dest = my_row+1;
    REG_PUTC(*d,dest);
  }
  else {
    REG_GETC(tt);  
    *d += *t;

    if((my_row&2)==0) {
      dest = my_row+2;
      REG_PUTC(dd,dest);
    }
    else {
      REG_GETC(tt);  
      *d += *t;
      
      if(my_row==3) {
        dest = 7;
        REG_PUTC(dd,dest);
      }
      else {
        REG_GETC(tt);  
        *d += *t;

        if((my_col&1)==0) {
          dest = my_col+1;
          REG_PUTR(dd,dest);
        }
        else {
          REG_GETR(tt);  
          *d += *t;

          if((my_col&2)==0) {
            dest = my_col+2;
            REG_PUTR(dd,dest);
          }
          else {
            REG_GETR(tt);  
            *d += *t;

            if(my_col==3) {
               dest = 7;
               REG_PUTR(dd,dest);
            }
            else {
              REG_GETR(tt);  
              *d += *t;
            }
          }
        }
      }
    }
  }
  *input = *d;
}

static inline void min_ulong(unsigned long *input) {
  doublev4 dd,tt;
  unsigned long *d=(unsigned long*)&dd;
  unsigned long *t=(unsigned long*)&tt;
  *d = *input;
  unsigned int dest;
  if((my_row&1)==0) {
    dest = my_row+1;
    REG_PUTC(dd,dest);
  }
  else {
    REG_GETC(tt);  
    *d = (*d<*t)?*d:*t;

    if((my_row&2)==0) {
      dest = my_row+2;
      REG_PUTC(dd,dest);
    }
    else {
      REG_GETC(tt);  
      *d = (*d<*t)?*d:*t;
      
      if(my_row==3) {
        dest = 7;
        REG_PUTC(dd,dest);
      }
      else {
        REG_GETC(tt);  
        *d = (*d<*t)?*d:*t;

        if((my_col&1)==0) {
          dest = my_col+1;
          REG_PUTR(dd,dest);
        }
        else {
          REG_GETR(tt);  
          *d = (*d<*t)?*d:*t;

          if((my_col&2)==0) {
            dest = my_col+2;
            REG_PUTR(dd,dest);
          }
          else {
            REG_GETR(tt);  
            *d = (*d<*t)?*d:*t;

            if(my_col==3) {
               dest = 7;
               REG_PUTR(dd,dest);
            }
            else {
              REG_GETR(tt);  
              *d = (*d<*t)?*d:*t;
            }
          }
        }
      }
    }
  }
//  if (my_id==63) REG_PUTC(dd,0);
//  else if(my_id==7) {
//    REG_GETC(dd);
//    REG_PUTR(dd,0);
//  }
//  else if(my_id==0) REG_GETR(dd);
  *input = *d;
}

static inline void max_ulong(unsigned long *input) {
  doublev4 dd,tt;
  unsigned long *d=(unsigned long*)&dd;
  unsigned long *t=(unsigned long*)&tt;
  *d = *input;
  unsigned int dest;
  if((my_row&1)==0) {
    dest = my_row+1;
    REG_PUTC(dd,dest);
  }
  else {
    REG_GETC(tt);  
    *d = (*d>*t)?*d:*t;

    if((my_row&2)==0) {
      dest = my_row+2;
      REG_PUTC(dd,dest);
    }
    else {
      REG_GETC(tt);  
      *d = (*d>*t)?*d:*t;
      
      if(my_row==3) {
        dest = 7;
        REG_PUTC(dd,dest);
      }
      else {
        REG_GETC(tt);  
        *d = (*d>*t)?*d:*t;

        if((my_col&1)==0) {
          dest = my_col+1;
          REG_PUTR(dd,dest);
        }
        else {
          REG_GETR(tt);  
          *d = (*d>*t)?*d:*t;

          if((my_col&2)==0) {
            dest = my_col+2;
            REG_PUTR(dd,dest);
          }
          else {
            REG_GETR(tt);  
            *d = (*d>*t)?*d:*t;

            if(my_col==3) {
               dest = 7;
               REG_PUTR(dd,dest);
            }
            else {
              REG_GETR(tt);  
              *d = (*d>*t)?*d:*t;
            }
          }
        }
      }
    }
  }
//  if (my_id==63) REG_PUTC(dd,0);
//  else if(my_id==7) {
//    REG_GETC(dd);
//    REG_PUTR(dd,0);
//  }
//  else if(my_id==0) REG_GETR(dd);
  *input = *d;
}

// boardcasting (frozen if row_cut and col_cut are incosistent with the using CPE core number)
// n: data size in unit of double
static inline void cpe_broadcast(const int n, double* data) {
  int i;
  doublev4 ttmp;
  if(my_id==0) {
    for (i=0;i<n;i+=4) {
      simd_load(ttmp,&data[i]);
      REG_PUTR(ttmp,8);
      REG_PUTC(ttmp,8);
    }
  }
  else {
    if(my_row==0) {
      for (i=0;i<n;i+=4) {
        REG_GETR(ttmp);
        REG_PUTC(ttmp,8);
        simd_store(ttmp,&data[i]);
      }
    }
    else {
      for (i=0;i<n;i+=4) {
        REG_GETC(ttmp);
        simd_store(ttmp,&data[i]);
      }
    }
  }
}

#ifdef SPLIT
// Longer inner-loop (512) version
static inline void CalcGrav_p_split2(
		const int         n_epi, 
		const EpiMM      *epi,     // double r[4], v[4]
		const EpjMM      *epj,     // double r[4], v[4]
		ForceMM4         *force,   // double ax[4], ay[4], az[4], pot[4];
		const Force_coff *pcoff,   // doublev4 m_r_coll_inv_qv4, k_v_miv4, eta_v_miv4;   
		Work             *work     // doublev4 dx, dy, dz, rdotv, r2, rinv
		)
{
#ifdef DEEP_PROFILE
	long clk0 = rtc_();
#endif

	posvel_loop(n_epi, epi, epj, work);

#ifdef DEEP_PROFILE
	long clk1 = rtc_();
#endif

	long nloop = 4*(n_epi/4 + (n_epi%4 ? 1 : 0) );
	rsqrt_loop(nloop, work);

#ifdef DEEP_PROFILE
	long clk2 = rtc_();
#endif

	force_loop(n_epi, epj, force, pcoff, work);

#ifdef DEEP_PROFILE
	long clk3 = rtc_();

	long clk01 = clk1 - clk0;
	long clk12 = clk2 - clk1;
	long clk23 = clk3 - clk2;

	if(ID_CHECK == my_id){
		long nloop = 4*(n_epi/4 + (n_epi%4 ? 1 : 0) );
		printf("n_epi=%d, nloop=%ld\n", n_epi, nloop);
		printf("loop1: %f cycle/loop\n", (double)clk01 / nloop);
		printf("loop2: %f cycle/loop\n", (double)clk12 / nloop);
		printf("loop3: %f cycle/loop\n", (double)clk23 / nloop);
	}
#endif

}

// Loop split version
static inline void CalcGrav_p_split(
		const int         n_epi, 
		const EpiMM      *epi,     // double r[4], v[4]
		const EpjMM      *epj,     // double r[4], v[4]
		ForceMM4         *force,   // double ax[4], ay[4], az[4], pot[4];
		const Force_coff *pcoff,   // doublev4 m_r_coll_inv_qv4, k_v_miv4, eta_v_miv4;   
		Work             *work     // doublev4 dx, dy, dz, rdotv, r2, rinv
		)
{
	long clk01 = 0;
	long clk12 = 0;
	long clk23 = 0;

	int j;
	for(j=0; j<4; j++){
		doublev4 rj = *(doublev4 *)epj[j].r;
		doublev4 vj = *(doublev4 *)epj[j].v;
		const doublev4 xj = simd_vshff(rj, rj, 0x00);
		const doublev4 yj = simd_vshff(rj, rj, 0x55);
		const doublev4 zj = simd_vshff(rj, rj, 0xaa);
		const doublev4 mj = simd_vshff(rj, rj, 0xff);

		const doublev4 vxj = simd_vshff(vj, vj, 0x00);
		const doublev4 vyj = simd_vshff(vj, vj, 0x55);
		const doublev4 vzj = simd_vshff(vj, vj, 0xaa);

		int i;
		const doublev4 tiny = *(doublev4 *)&glb_tiny;

		long clk0 = rtc_();

		for(i=0; i<n_epi; i+=4){
			const doublev4 xi = *(doublev4 *)epi[i+0].r;
			const doublev4 yi = *(doublev4 *)epi[i+1].r;
			const doublev4 zi = *(doublev4 *)epi[i+2].r;

			const doublev4 vxi = *(doublev4 *)epi[i+0].v;
			const doublev4 vyi = *(doublev4 *)epi[i+1].v;
			const doublev4 vzi = *(doublev4 *)epi[i+2].v;

			const doublev4 dx = xi - xj; 
			const doublev4 dy = yi - yj; 
			const doublev4 dz = zi - zj; 

			const doublev4 dvx = vxi - vxj; 
			const doublev4 dvy = vyi - vyj; 
			const doublev4 dvz = vzi - vzj; 

			const doublev4 r2    = tiny + dx*dx + dy*dy + dz*dz;
			const doublev4 rdotv = dx*dvx + dy*dvy + dz*dvz; 

			work[i/4].dx    = dx;
			work[i/4].dy    = dy;
			work[i/4].dz    = dz;
			work[i/4].rdotv = rdotv;
			work[i/4].r2    = r2;
		}

		long clk1 = rtc_();

#if 0
		for(i=0; i<n_epi; i+=4){
			const doublev4 r2 = work[i/4].r2;
			const doublev4 rinv = 1.0 / simd_vsqrtd(r2);
			work[i/4].rinv = rinv;
		}
#else
		long nloop = n_epi/4 + (n_epi%4 ? 1 : 0);
		rsqrt_loop(nloop, work);
#endif

		long clk2 = rtc_();

		for(i=0; i<n_epi; i+=4){
			const doublev4 dx    = work[i/4].dx;
			const doublev4 dy    = work[i/4].dy;
			const doublev4 dz    = work[i/4].dz;
			const doublev4 rdotv = work[i/4].rdotv;
			// const doublev4 r_sq  = work[i/4].r2;
			const doublev4 r_inv = work[i/4].rinv;

			const doublev4 pij = mj * r_inv;
			const doublev4 r_inv_sq = r_inv * r_inv;
			const doublev4 xij      = d1v4 - pcoff->r_collv4*r_inv;
			const doublev4 rvji     = rdotv * r_inv_sq; 

			const doublev4 fsd = pcoff->m_r_coll_inv_qv4 + pcoff->k_v_miv4 * xij + pcoff->eta_v_miv4 * rvji;
			const doublev4 mri3 = simd_vsellt(xij, fsd, r_inv_sq * pij);

			doublev4 fx = *(doublev4 *)force[i/4].ax;
			doublev4 fy = *(doublev4 *)force[i/4].ay;
			doublev4 fz = *(doublev4 *)force[i/4].az;
			doublev4 fp = *(doublev4 *)force[i/4].pot;

			fx -= mri3 * dx;
			fy -= mri3 * dy;
			fz -= mri3 * dz;
			fp -= pij;

			*(doublev4 *)force[i/4].ax  = fx;
			*(doublev4 *)force[i/4].ay  = fy;
			*(doublev4 *)force[i/4].az  = fz;
			*(doublev4 *)force[i/4].pot = fp;
		}

		long clk3 = rtc_();

		clk01 += clk1 - clk0;
		clk12 += clk2 - clk1;
		clk23 += clk3 - clk2;
	} // for(j)

#ifdef DEEP_PROFILE
	if(ID_CHECK == my_id){
		long nloop = 4 * ( n_epi/4 + (n_epi%4 ? 1 : 0) );
		printf("n_epi=%d, nloop=%ld\n", n_epi, nloop);
		printf("loop1: %f cycle/loop\n", (double)clk01 / nloop);
		printf("loop2: %f cycle/loop\n", (double)clk12 / nloop);
		printf("loop3: %f cycle/loop\n", (double)clk23 / nloop);
	}
#endif
}

#elif !defined(FULL_ASM)
void CalcGrav_p(const int n_epi, const EpiMM* epi, const EpjMM* epj, ForceMM4* force, const Force_coff* pcoff){
  int i;
//  doublev4 r_coll_sq4,eps2v4;
//  simd_load(r_coll_sq4, r_coll_sq);
//  simd_load(eps2v4    , eps2     );

  doublev4 pj1,pj2,pj3,pj4;
//! x,y,z,m
  simd_load(pj1,(double*)epj[0].r);
  simd_load(pj2,(double*)epj[1].r);
  simd_load(pj3,(double*)epj[2].r);
  simd_load(pj4,(double*)epj[3].r);

#ifdef DEEP_DEBUG
  if(my_id==ID_CHECK) {
    printf("load epj[0].r=%e\n",epj[3].r[0]);
    simd_print_doublev4(pj1);
    simd_print_doublev4(pj2);
    simd_print_doublev4(pj3);
    simd_print_doublev4(pj4);
  }
#endif

  const doublev4 xj1 = simd_vshff(pj1,pj1,0x00);
  const doublev4 yj1 = simd_vshff(pj1,pj1,0x55);
  const doublev4 zj1 = simd_vshff(pj1,pj1,0xaa);
  const doublev4 mj1 = simd_vshff(pj1,pj1,0xff);

  const doublev4 xj2 = simd_vshff(pj2,pj2,0x00);
  const doublev4 yj2 = simd_vshff(pj2,pj2,0x55);
  const doublev4 zj2 = simd_vshff(pj2,pj2,0xaa);
  const doublev4 mj2 = simd_vshff(pj2,pj2,0xff);

  const doublev4 xj3 = simd_vshff(pj3,pj3,0x00);
  const doublev4 yj3 = simd_vshff(pj3,pj3,0x55);
  const doublev4 zj3 = simd_vshff(pj3,pj3,0xaa);
  const doublev4 mj3 = simd_vshff(pj3,pj3,0xff);

  const doublev4 xj4 = simd_vshff(pj4,pj4,0x00);
  const doublev4 yj4 = simd_vshff(pj4,pj4,0x55);
  const doublev4 zj4 = simd_vshff(pj4,pj4,0xaa);
  const doublev4 mj4 = simd_vshff(pj4,pj4,0xff);

  //! vx,vy,vz,id
  simd_load(pj1,(double*)epj[0].v);
  simd_load(pj2,(double*)epj[1].v);
  simd_load(pj3,(double*)epj[2].v);
  simd_load(pj4,(double*)epj[3].v);

  const doublev4 vxj1 = simd_vshff(pj1,pj1,0x00);
  const doublev4 vyj1 = simd_vshff(pj1,pj1,0x55);
  const doublev4 vzj1 = simd_vshff(pj1,pj1,0xaa);
  const doublev4 idj1 = simd_vshff(pj1,pj1,0xff);

  const doublev4 vxj2 = simd_vshff(pj2,pj2,0x00);
  const doublev4 vyj2 = simd_vshff(pj2,pj2,0x55);
  const doublev4 vzj2 = simd_vshff(pj2,pj2,0xaa);
  const doublev4 idj2 = simd_vshff(pj2,pj2,0xff);

  const doublev4 vxj3 = simd_vshff(pj3,pj3,0x00);
  const doublev4 vyj3 = simd_vshff(pj3,pj3,0x55);
  const doublev4 vzj3 = simd_vshff(pj3,pj3,0xaa);
  const doublev4 idj3 = simd_vshff(pj3,pj3,0xff);

  const doublev4 vxj4 = simd_vshff(pj4,pj4,0x00);
  const doublev4 vyj4 = simd_vshff(pj4,pj4,0x55);
  const doublev4 vzj4 = simd_vshff(pj4,pj4,0xaa);
  const doublev4 idj4 = simd_vshff(pj4,pj4,0xff);
  
  for(i=0; i<n_epi; i+=4){

    doublev4 xi,yi,zi,mi;
    simd_load(xi, (double*)epi[i  ].r);
    simd_load(yi, (double*)epi[i+1].r);
    simd_load(zi, (double*)epi[i+2].r);
    simd_load(mi, (double*)epi[i+3].r);

    doublev4 vxi,vyi,vzi,idi;
    simd_load(vxi,(double*)epi[i  ].v);
    simd_load(vyi,(double*)epi[i+1].v);
    simd_load(vzi,(double*)epi[i+2].v);
    simd_load(idi,(double*)epi[i+3].v);
    
    //    idi &= ilowmask;
    //const doublev4 idid = *(doublev4*) &idi;

#ifdef DEEP_DEBUG
    if(my_id==ID_CHECK) {
      simd_print_doublev4(xi);
      printf("x=%e %e %e %e\n",
             epi[i].r[0],
             epi[i].r[1],
             epi[i].r[2],
             epi[i].r[3]);
    }
#endif

    doublev4 fx,fy,fz,fp;
    int ii=i>>2;
    simd_load(fx, force[ii].ax);
    simd_load(fy, force[ii].ay);
    simd_load(fz, force[ii].az);
    simd_load(fp, force[ii].pot);

    //doublev4 frmin;
    //simd_load(frmin,force[i].r_ngb_sq);

    //intv8 fn;
    //simd_load(fn,force[i].icoll[0]);


    // const doublev4 dx1 = xi - simd_vshff(xj,xj,0x00);
    // const doublev4 dx2 = xi - simd_vshff(xj,xj,0x55);
    // const doublev4 dx3 = xi - simd_vshff(xj,xj,0xaa);
    // const doublev4 dx4 = xi - simd_vshff(xj,xj,0xff);
    //  
    // const doublev4 dy1 = yi - simd_vshff(yj,yj,0x00);
    // const doublev4 dy2 = yi - simd_vshff(yj,yj,0x55);
    // const doublev4 dy3 = yi - simd_vshff(yj,yj,0xaa);
    // const doublev4 dy4 = yi - simd_vshff(yj,yj,0xff);
    //  
    // const doublev4 dz1 = zi - simd_vshff(zj,zj,0x00);
    // const doublev4 dz2 = zi - simd_vshff(zj,zj,0x55);
    // const doublev4 dz3 = zi - simd_vshff(zj,zj,0xaa);
    // const doublev4 dz4 = zi - simd_vshff(zj,zj,0xff);
    //     
    // const doublev4 mj1 = simd_vshff(mj,mj,0x00);
    // const doublev4 mj2 = simd_vshff(mj,mj,0x55);
    // const doublev4 mj3 = simd_vshff(mj,mj,0xaa);
    // const doublev4 mj4 = simd_vshff(mj,mj,0xff);

    const doublev4 dx1 = xi - xj1; 
    const doublev4 dx2 = xi - xj2; 
    const doublev4 dx3 = xi - xj3; 
    const doublev4 dx4 = xi - xj4; 
                              
    const doublev4 dy1 = yi - yj1; 
    const doublev4 dy2 = yi - yj2; 
    const doublev4 dy3 = yi - yj3; 
    const doublev4 dy4 = yi - yj4; 
                              
    const doublev4 dz1 = zi - zj1; 
    const doublev4 dz2 = zi - zj2; 
    const doublev4 dz3 = zi - zj3; 
    const doublev4 dz4 = zi - zj4; 
    
//    const doublev4 r_sq_org1 = *eps2v4 + dx1*dx1 + dy1*dy1 + dz1*dz1;
//    const doublev4 r_sq_org2 = *eps2v4 + dx2*dx2 + dy2*dy2 + dz2*dz2;
//    const doublev4 r_sq_org3 = *eps2v4 + dx3*dx3 + dy3*dy3 + dz3*dz3;
//    const doublev4 r_sq_org4 = *eps2v4 + dx4*dx4 + dy4*dy4 + dz4*dz4;

    const doublev4 r_sq_org1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
    const doublev4 r_sq_org2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    const doublev4 r_sq_org3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
    const doublev4 r_sq_org4 = dx4*dx4 + dy4*dy4 + dz4*dz4;

    
    // intv8 id_adrj;
    // simd_load(id_adrj,(int*)epj[j].id[0]);
    // const intv8 idj = id_adrj & ilowmask;
    // const doublev4 idjd = *(doublev4*)&idj;

    //#ifdef DEEP_DEBUG
    //  if(my_id==ID_CHECK) simd_print_intv8(idj);
    //#endif
        
    // const doublev4 idj1 = simd_vshff(idjd,idjd,0x00);
    // const doublev4 idj2 = simd_vshff(idjd,idjd,0x55);
    // const doublev4 idj3 = simd_vshff(idjd,idjd,0xaa);
    // const doublev4 idj4 = simd_vshff(idjd,idjd,0xff);
    
    const doublev4 mask1 = simd_vfcmpeq(idi,idj1);
    const doublev4 mask2 = simd_vfcmpeq(idi,idj2);
    const doublev4 mask3 = simd_vfcmpeq(idi,idj3);
    const doublev4 mask4 = simd_vfcmpeq(idi,idj4);

    const doublev4 r_sq1 = simd_vseleq(mask1, r_sq_org1, dhuge);
    const doublev4 r_sq2 = simd_vseleq(mask2, r_sq_org2, dhuge);
    const doublev4 r_sq3 = simd_vseleq(mask3, r_sq_org3, dhuge); 
    const doublev4 r_sq4 = simd_vseleq(mask4, r_sq_org4, dhuge);

#ifdef DEEP_DEBUG
    if(my_id==ID_CHECK) {
      printf("idi:\n");
      simd_print_intv8(*(intv8*)&idi);
      printf("idj:\n");
      simd_print_intv8(*(intv8*)&idj1);
      simd_print_intv8(*(intv8*)&idj2);
      simd_print_intv8(*(intv8*)&idj3);
      simd_print_intv8(*(intv8*)&idj4);
      printf("mask:\n");
      simd_print_doublev4(mask1);
      simd_print_doublev4(mask2);
      simd_print_doublev4(mask3);
      simd_print_doublev4(mask4);
    }
#endif        

    // const doublev4 r_sq1 = eps2v4 + r_sq_real1;
    // const doublev4 r_sq2 = eps2v4 + r_sq_real2;
    // const doublev4 r_sq3 = eps2v4 + r_sq_real3;
    // const doublev4 r_sq4 = eps2v4 + r_sq_real4;

//      doublev4 r_inv1, r_inv2, r_inv3, r_inv4;
//      rsqrtd(&r_sq1, &r_inv1);
//      rsqrtd(&r_sq2, &r_inv2);
//      rsqrtd(&r_sq3, &r_inv3);
//      rsqrtd(&r_sq4, &r_inv4);

    const doublev4 r_inv1 = d1v4 / simd_vsqrtd(r_sq1);
    const doublev4 r_inv2 = d1v4 / simd_vsqrtd(r_sq2);
    const doublev4 r_inv3 = d1v4 / simd_vsqrtd(r_sq3);
    const doublev4 r_inv4 = d1v4 / simd_vsqrtd(r_sq4);

    const doublev4 r_inv_sq1 = r_inv1 * r_inv1;
    const doublev4 r_inv_sq2 = r_inv2 * r_inv2;
    const doublev4 r_inv_sq3 = r_inv3 * r_inv3;
    const doublev4 r_inv_sq4 = r_inv4 * r_inv4;

//      const doublev4 r_inv1 = 0.5 / simd_vsqrtd(r_sq1) * mask1;
//      const doublev4 r_inv2 = 0.5 / simd_vsqrtd(r_sq2) * mask2;
//      const doublev4 r_inv3 = 0.5 / simd_vsqrtd(r_sq3) * mask3;
//      const doublev4 r_inv4 = 0.5 / simd_vsqrtd(r_sq4) * mask4;

    const doublev4 rvxji1 = (vxi - vxj1) * dx1;
    const doublev4 rvxji2 = (vxi - vxj2) * dx2;
    const doublev4 rvxji3 = (vxi - vxj3) * dx3;
    const doublev4 rvxji4 = (vxi - vxj4) * dx4;

    const doublev4 rvyji1 = (vyi - vyj1) * dy1;
    const doublev4 rvyji2 = (vyi - vyj2) * dy2;
    const doublev4 rvyji3 = (vyi - vyj3) * dy3;
    const doublev4 rvyji4 = (vyi - vyj4) * dy4;

    const doublev4 rvzji1 = (vzi - vzj1) * dz1;
    const doublev4 rvzji2 = (vzi - vzj2) * dz2;
    const doublev4 rvzji3 = (vzi - vzj3) * dz3;
    const doublev4 rvzji4 = (vzi - vzj4) * dz4;

    const doublev4 rvji1 = (rvxji1 + rvyji1 + rvzji1) * r_inv_sq1; 
    const doublev4 rvji2 = (rvxji2 + rvyji2 + rvzji2) * r_inv_sq2; 
    const doublev4 rvji3 = (rvxji3 + rvyji3 + rvzji3) * r_inv_sq3; 
    const doublev4 rvji4 = (rvxji4 + rvyji4 + rvzji4) * r_inv_sq4; 

    const doublev4 pij1 = mj1 * r_inv1;
    const doublev4 pij2 = mj2 * r_inv2;
    const doublev4 pij3 = mj3 * r_inv3;
    const doublev4 pij4 = mj4 * r_inv4;

    doublev4 mri3_1 = r_inv_sq1*pij1;
    doublev4 mri3_2 = r_inv_sq2*pij2;
    doublev4 mri3_3 = r_inv_sq3*pij3;
    doublev4 mri3_4 = r_inv_sq4*pij4;

    const doublev4 xij1 = d1v4 - pcoff->r_collv4*r_inv1;
    const doublev4 xij2 = d1v4 - pcoff->r_collv4*r_inv2;
    const doublev4 xij3 = d1v4 - pcoff->r_collv4*r_inv3;
    const doublev4 xij4 = d1v4 - pcoff->r_collv4*r_inv4;
    
    const doublev4 fsd1 = pcoff->k_v_miv4 * xij1 + pcoff->eta_v_miv4 * rvji1 + pcoff->m_r_coll_inv_qv4;
    const doublev4 fsd2 = pcoff->k_v_miv4 * xij2 + pcoff->eta_v_miv4 * rvji2 + pcoff->m_r_coll_inv_qv4;
    const doublev4 fsd3 = pcoff->k_v_miv4 * xij3 + pcoff->eta_v_miv4 * rvji3 + pcoff->m_r_coll_inv_qv4;
    const doublev4 fsd4 = pcoff->k_v_miv4 * xij4 + pcoff->eta_v_miv4 * rvji4 + pcoff->m_r_coll_inv_qv4;

    mri3_1 = simd_vsellt(xij1, fsd1, mri3_1);
    mri3_2 = simd_vsellt(xij2, fsd2, mri3_2);
    mri3_3 = simd_vsellt(xij3, fsd3, mri3_3);
    mri3_4 = simd_vsellt(xij4, fsd4, mri3_4);

    
//    if(my_id==ID_CHECK) {
//      printf("\n");
//      simd_print_doublev4(xij1);
//      simd_print_doublev4(xij2);
//      simd_print_doublev4(xij3);
//      simd_print_doublev4(xij4);
//    }
        
    const doublev4 mdx1 = mri3_1 * dx1;
    const doublev4 mdx2 = mri3_2 * dx2;
    const doublev4 mdx3 = mri3_3 * dx3;
    const doublev4 mdx4 = mri3_4 * dx4;

    const doublev4 mdy1 = mri3_1 * dy1;
    const doublev4 mdy2 = mri3_2 * dy2;
    const doublev4 mdy3 = mri3_3 * dy3;
    const doublev4 mdy4 = mri3_4 * dy4;
                                  
    const doublev4 mdz1 = mri3_1 * dz1;
    const doublev4 mdz2 = mri3_2 * dz2;
    const doublev4 mdz3 = mri3_3 * dz3;
    const doublev4 mdz4 = mri3_4 * dz4;

    //      dv4_trans(&mdx, &mdy, &mdz, &pij);
            
    fx -= mdx1 + mdx2 + mdx3 + mdx4;
    fy -= mdy1 + mdy2 + mdy3 + mdy4;
    fz -= mdz1 + mdz2 + mdz3 + mdz4;
    fp -= pij1 + pij2 + pij3 + pij4;

    //        f4 -= mdx + mdy + mdz + pij;

    //        simd_store(f4,&force[i]);

//    const doublev4 rcmp1 = simd_vfcmplt(r_sq1, r_coll_sqv4);
//    const doublev4 rcmp2 = simd_vfcmplt(r_sq2, r_coll_sqv4);
//    const doublev4 rcmp3 = simd_vfcmplt(r_sq3, r_coll_sqv4);
//    const doublev4 rcmp4 = simd_vfcmplt(r_sq4, r_coll_sqv4);

//    intv8 ncoll1, ncoll2, ncoll3, ncoll4;
// 
//    asm volatile("vseleq %1,%2,%3,%0" : "=r"(ncoll1) : "r"(rcmp1), "r"(i0v8), "r"(ilowone));
//    asm volatile("vseleq %1,%2,%3,%0" : "=r"(ncoll2) : "r"(rcmp2), "r"(i0v8), "r"(ilowone));
//    asm volatile("vseleq %1,%2,%3,%0" : "=r"(ncoll3) : "r"(rcmp3), "r"(i0v8), "r"(ilowone));
//    asm volatile("vseleq %1,%2,%3,%0" : "=r"(ncoll4) : "r"(rcmp4), "r"(i0v8), "r"(ilowone));
//      const doublev4 dncoll1 = simd_vseleq(rcmp1, d0v4, dlowone);
//      const doublev4 dncoll2 = simd_vseleq(rcmp2, d0v4, dlowone);
//      const doublev4 dncoll3 = simd_vseleq(rcmp3, d0v4, dlowone);
//      const doublev4 dncoll4 = simd_vseleq(rcmp4, d0v4, dlowone);

//      const intv8 ncoll1 = *(intv8*)&dncoll1;
//      const intv8 ncoll2 = *(intv8*)&dncoll2;
//      const intv8 ncoll3 = *(intv8*)&dncoll3;
//      const intv8 ncoll4 = *(intv8*)&dncoll4;
//#ifdef DEEP_DEBUG
//      if(my_id==ID_CHECK) {
//        simd_print_doublev4(r_sq_real1);
//        simd_print_doublev4(r_coll_sqv4);
//        simd_print_doublev4(rcmp1);
//        simd_print_intv8(ncoll1);
//        simd_print_intv8(ncoll2);
//        simd_print_intv8(ncoll3);
//        simd_print_intv8(ncoll4);
//      }
//#endif
      
//      intv8 fnc  = fn & ilowmask;
//      fnc += ncoll1 + ncoll2 + ncoll3 + ncoll4;
// 
//#ifdef DEEP_DEBUG
//      if(my_id==ID_CHECK) simd_print_intv8(fnc);
//#endif

      //      asm volatile("vshfw %1,%1,%2,%0" : "=r"(ida1) : "r"(id_adrj), "r"(0x00));
//      const doublev4 ida1= simd_vshff(*(doublev4*)&id_adrj,*(doublev4*)&id_adrj,0x00);
//      const doublev4 ida2= simd_vshff(*(doublev4*)&id_adrj,*(doublev4*)&id_adrj,0x55);
//      const doublev4 ida3= simd_vshff(*(doublev4*)&id_adrj,*(doublev4*)&id_adrj,0xaa);
//      const doublev4 ida4= simd_vshff(*(doublev4*)&id_adrj,*(doublev4*)&id_adrj,0xff);
      
//      const intv8 ida1= simd_vshff(id_adrj,id_adrj,0x00000000);
//      const intv8 ida2= simd_vshff(id_adrj,id_adrj,0x22222222);
//      const intv8 ida3= simd_vshff(id_adrj,id_adrj,0x44444444);
//      const intv8 ida4= simd_vshff(id_adrj,id_adrj,0x66666666);
// 
//      const doublev4 rcmpm1= simd_vfcmplt(r_sq_real1, frmin);
//      asm volatile("vseleq %1,%2,%3,%0" : "=r"(fn) : "r"(rcmpm1), "r"(fn), "r"(ida1));
//      // fid   = simd_vseleq(rcmpm1, fid,   *(doublev4*)&ida1);
//      frmin = simd_vseleq(rcmpm1, frmin, r_sq_real1);
// 
//#ifdef DEEP_DEBUG
//      if(my_id==ID_CHECK) {
//        simd_print_doublev4(rcmpm1);
//        simd_print_intv8(ida1);
//        simd_print_intv8(fn);
//      }
//#endif
// 
//      const doublev4 rcmpm2= simd_vfcmplt(r_sq_real2, frmin);
//      asm volatile("vseleq %1,%2,%3,%0" : "=r"(fn) : "r"(rcmpm2), "r"(fn), "r"(ida2));
//      //      const doublev4 fid  = simd_vseleq(rcmpm2, *(doublev4*)&fn,   *(doublev4*)&ida2);
//      //      fn = *(intv8*)&fid;
//      frmin = simd_vseleq(rcmpm2, frmin, r_sq_real2);
// 
//#ifdef DEEP_DEBUG
//      if(my_id==ID_CHECK) {
//        simd_print_doublev4(rcmpm2);
//        simd_print_intv8(ida2);
//        simd_print_intv8(fn);
//      }
//#endif
// 
//      const doublev4 rcmpm3= simd_vfcmplt(r_sq_real3, frmin);
//      asm volatile("vseleq %1,%2,%3,%0" : "=r"(fn) : "r"(rcmpm3), "r"(fn), "r"(ida3));
//      // fid   = simd_vseleq(rcmpm3, fid,   *(doublev4*)&ida3);
//      frmin = simd_vseleq(rcmpm3, frmin, r_sq_real3);
// 
//#ifdef DEEP_DEBUG
//      if(my_id==ID_CHECK) {
//        simd_print_doublev4(rcmpm3);
//        simd_print_intv8(ida3);
//        simd_print_intv8(fn);
//      }      
//#endif
//      
//      const doublev4 rcmpm4= simd_vfcmplt(r_sq_real4, frmin);
//      asm volatile("vseleq %1,%2,%3,%0" : "=r"(fn) : "r"(rcmpm4), "r"(fn), "r"(ida4));
//      // fid   = simd_vseleq(rcmpm4, fid,   *(doublev4*)&ida4); 
//      frmin = simd_vseleq(rcmpm4, frmin, r_sq_real4);
// 
//#ifdef DEEP_DEBUG
//      if(my_id==ID_CHECK) {
//        simd_print_doublev4(rcmpm4);
//        simd_print_intv8(ida4);
//        simd_print_intv8(fn);
//      }      
//#endif
// 
//      fn = (fn & ihighmask) + fnc;
// 
//#ifdef DEEP_DEBUG
//      if(my_id==ID_CHECK) simd_print_intv8(fn);
//#endif
// 
//    }
    simd_store(fx,force[ii].ax);
    simd_store(fy,force[ii].ay);
    simd_store(fz,force[ii].az);
    simd_store(fp,force[ii].pot);

#ifdef DEEP_DEBUG
    if(my_id==ID_CHECK) {
      printf("force[%d].ax=%e %e %e %e\n",ii,force[ii].ax[0],force[ii].ax[1],force[ii].ax[2],force[ii].ax[3]);
      printf("force[%d].ay=%e %e %e %e\n",ii,force[ii].ay[0],force[ii].ay[1],force[ii].ay[2],force[ii].ay[3]);
      printf("force[%d].az=%e %e %e %e\n",ii,force[ii].az[0],force[ii].az[1],force[ii].az[2],force[ii].az[3]);
      printf("force[%d].pot=%e %e %e %e\n",ii,force[ii].pot[0],force[ii].pot[1],force[ii].pot[2],force[ii].pot[3]);
    }    
#endif

    // simd_store(frmin,force[i].r_ngb_sq);
    // simd_store(fn,&(force[i].icoll[0][0]));

//#ifdef DEEP_DEBUG
//    if(my_id==ID_CHECK) {
//      simd_print_intv8(fn);
//      printf("icoll[%d]=%d %d %d %d %d %d %d %d\n",i,
//           force[i].icoll[0][0],force[i].icoll[0][1],
//           force[i].icoll[1][0],force[i].icoll[1][1],
//           force[i].icoll[2][0],force[i].icoll[2][1],
//           force[i].icoll[3][0],force[i].icoll[3][1]);
//    }    
//#endif
  }
}
#endif


#ifdef SPLIT
static inline CalcGrav_ps_split(const int n_epi, const EpiMM* epi, const EpiMM* sat, ForceMM4* forcei, ForceMM* forces, const Force_coff* ps_coff, const Force_coff* sp_coff, Work* work){

#ifdef DEEP_PROFILE
	long clk0 = rtc_();
#endif

	posvel_loop(n_epi, epi, sat, work);

#ifdef DEEP_PROFILE
	long clk1 = rtc_();
#endif

	long nloop = 4*(n_epi/4 + (n_epi%4 ? 1 : 0) );
	rsqrt_loop(nloop, work);

#ifdef DEEP_PROFILE
	long clk2 = rtc_();
#endif

	force_loop(n_epi, sat, forcei, ps_coff, work);

#ifdef DEEP_PROFILE
	long clk3 = rtc_();
#endif

    force_counter_loop(n_epi, epi, sat, forces, sp_coff, work);

#ifdef DEEP_PROFILE
	long clk4 = rtc_();
#endif

#ifdef DEEP_PROFILE
	long clk01 = clk1 - clk0;
	long clk12 = clk2 - clk1;
	long clk23 = clk3 - clk2;
	long clk34 = clk4 - clk3;

	if(ID_CHECK == my_id){
		long nloop = 4*(n_epi/4 + (n_epi%4 ? 1 : 0) );
		printf("n_epi=%d, nloop=%ld\n", n_epi, nloop);
		printf("loop1: %f cycle/loop\n", (double)clk01 / nloop);
		printf("loop2: %f cycle/loop\n", (double)clk12 / nloop);
		printf("loop3: %f cycle/loop\n", (double)clk23 / nloop);
		printf("loop4: %f cycle/loop\n", (double)clk34 / nloop);
	}
#endif
    
}

#elif defined(FULL_ASM)
void collect_sat_force(ForceMM* forces, const doublev4* mp_over_ms) {
    doublev4 fsat[16];
    read_satellite_force(fsat);
    int i,i4;
    for(i=0;i<4;i++) {
      i4=i*4;
      dv4_trans(&fsat[i4],&fsat[1+i4],&fsat[2+i4],&fsat[3+i4]);
      *(doublev4*)&forces[i] += (fsat[i4]+fsat[1+i4]+fsat[2+i4]+fsat[3+i4])* *mp_over_ms;
    }
    return;
}

#else

__attribute__((noinline))
void CalcGrav_ps(const int n_epi, const EpiMM* epi, const EpiMM* sat, ForceMM4* forcei, ForceMM* forces, const Force_coff* ps_coff, const Force_coff* sp_coff){
  int i,j;

  for(i=0; i<n_epi; i+=4){

    doublev4 xi,yi,zi,mi;
    simd_load(xi, (double*)epi[i  ].r);
    simd_load(yi, (double*)epi[i+1].r);
    simd_load(zi, (double*)epi[i+2].r);
    simd_load(mi, (double*)epi[i+3].r);

    doublev4 vxi,vyi,vzi;
    simd_load(vxi,(double*)epi[i  ].v);
    simd_load(vyi,(double*)epi[i+1].v);
    simd_load(vzi,(double*)epi[i+2].v);
    
    //    idi &= ilowmask;
    //const doublev4 idid = *(doublev4*) &idi;

#ifdef DEEP_DEBUG
    if(my_id==ID_CHECK) {
      simd_print_doublev4(xi);
      printf("x=%e %e %e %e\n",
             epi[i].r[0],
             epi[i].r[1],
             epi[i].r[2],
             epi[i].r[3]);
    }
#endif

    doublev4 fx,fy,fz,fp;
    int ii=i>>2;
    simd_load(fx, forcei[ii].ax);
    simd_load(fy, forcei[ii].ay);
    simd_load(fz, forcei[ii].az);
    simd_load(fp, forcei[ii].pot);

    //    for(j=0; j<n_sat; j+=4){
    j=0;

      doublev4 pj1,pj2,pj3,pj4;
      //! x,y,z,m
      simd_load(pj1,(double*)sat[j  ].r);
      simd_load(pj2,(double*)sat[j+1].r);
      simd_load(pj3,(double*)sat[j+2].r);
      simd_load(pj4,(double*)sat[j+3].r);

#ifdef DEEP_DEBUG
      if(my_id==ID_CHECK) {
        printf("load sat[0].r=%e\n",sat[3].r[0]);
        simd_print_doublev4(pj1);
        simd_print_doublev4(pj2);
        simd_print_doublev4(pj3);
        simd_print_doublev4(pj4);
      }
#endif

      const doublev4 xj1 = simd_vshff(pj1,pj1,0x00);
      const doublev4 yj1 = simd_vshff(pj1,pj1,0x55);
      const doublev4 zj1 = simd_vshff(pj1,pj1,0xaa);
      const doublev4 mj1 = simd_vshff(pj1,pj1,0xff);

      const doublev4 xj2 = simd_vshff(pj2,pj2,0x00);
      const doublev4 yj2 = simd_vshff(pj2,pj2,0x55);
      const doublev4 zj2 = simd_vshff(pj2,pj2,0xaa);
      const doublev4 mj2 = simd_vshff(pj2,pj2,0xff);

      const doublev4 xj3 = simd_vshff(pj3,pj3,0x00);
      const doublev4 yj3 = simd_vshff(pj3,pj3,0x55);
      const doublev4 zj3 = simd_vshff(pj3,pj3,0xaa);
      const doublev4 mj3 = simd_vshff(pj3,pj3,0xff);

      const doublev4 xj4 = simd_vshff(pj4,pj4,0x00);
      const doublev4 yj4 = simd_vshff(pj4,pj4,0x55);
      const doublev4 zj4 = simd_vshff(pj4,pj4,0xaa);
      const doublev4 mj4 = simd_vshff(pj4,pj4,0xff);

      //! vx,vy,vz,id
      simd_load(pj1,(double*)sat[j  ].v);
      simd_load(pj2,(double*)sat[j+1].v);
      simd_load(pj3,(double*)sat[j+2].v);
      simd_load(pj4,(double*)sat[j+3].v);

      const doublev4 vxj1 = simd_vshff(pj1,pj1,0x00);
      const doublev4 vyj1 = simd_vshff(pj1,pj1,0x55);
      const doublev4 vzj1 = simd_vshff(pj1,pj1,0xaa);

      const doublev4 vxj2 = simd_vshff(pj2,pj2,0x00);
      const doublev4 vyj2 = simd_vshff(pj2,pj2,0x55);
      const doublev4 vzj2 = simd_vshff(pj2,pj2,0xaa);

      const doublev4 vxj3 = simd_vshff(pj3,pj3,0x00);
      const doublev4 vyj3 = simd_vshff(pj3,pj3,0x55);
      const doublev4 vzj3 = simd_vshff(pj3,pj3,0xaa);

      const doublev4 vxj4 = simd_vshff(pj4,pj4,0x00);
      const doublev4 vyj4 = simd_vshff(pj4,pj4,0x55);
      const doublev4 vzj4 = simd_vshff(pj4,pj4,0xaa);
  
      const doublev4 dx1 = xi - xj1; 
      const doublev4 dx2 = xi - xj2; 
      const doublev4 dx3 = xi - xj3; 
      const doublev4 dx4 = xi - xj4; 
                              
      const doublev4 dy1 = yi - yj1; 
      const doublev4 dy2 = yi - yj2; 
      const doublev4 dy3 = yi - yj3; 
      const doublev4 dy4 = yi - yj4; 
                              
      const doublev4 dz1 = zi - zj1; 
      const doublev4 dz2 = zi - zj2; 
      const doublev4 dz3 = zi - zj3; 
      const doublev4 dz4 = zi - zj4; 
    
      const doublev4 r_sq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
      const doublev4 r_sq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
      const doublev4 r_sq3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
      const doublev4 r_sq4 = dx4*dx4 + dy4*dy4 + dz4*dz4;

      const doublev4 r_inv1 = d1v4 / simd_vsqrtd(r_sq1);
      const doublev4 r_inv2 = d1v4 / simd_vsqrtd(r_sq2);
      const doublev4 r_inv3 = d1v4 / simd_vsqrtd(r_sq3);
      const doublev4 r_inv4 = d1v4 / simd_vsqrtd(r_sq4);

      const doublev4 r_inv_sq1 = r_inv1 * r_inv1;
      const doublev4 r_inv_sq2 = r_inv2 * r_inv2;
      const doublev4 r_inv_sq3 = r_inv3 * r_inv3;
      const doublev4 r_inv_sq4 = r_inv4 * r_inv4;

      const doublev4 rvxji1 = (vxi - vxj1) * dx1;
      const doublev4 rvxji2 = (vxi - vxj2) * dx2;
      const doublev4 rvxji3 = (vxi - vxj3) * dx3;
      const doublev4 rvxji4 = (vxi - vxj4) * dx4;

      const doublev4 rvyji1 = (vyi - vyj1) * dy1;
      const doublev4 rvyji2 = (vyi - vyj2) * dy2;
      const doublev4 rvyji3 = (vyi - vyj3) * dy3;
      const doublev4 rvyji4 = (vyi - vyj4) * dy4;

      const doublev4 rvzji1 = (vzi - vzj1) * dz1;
      const doublev4 rvzji2 = (vzi - vzj2) * dz2;
      const doublev4 rvzji3 = (vzi - vzj3) * dz3;
      const doublev4 rvzji4 = (vzi - vzj4) * dz4;

      const doublev4 rvji1 = (rvxji1 + rvyji1 + rvzji1) * r_inv_sq1; 
      const doublev4 rvji2 = (rvxji2 + rvyji2 + rvzji2) * r_inv_sq2; 
      const doublev4 rvji3 = (rvxji3 + rvyji3 + rvzji3) * r_inv_sq3; 
      const doublev4 rvji4 = (rvxji4 + rvyji4 + rvzji4) * r_inv_sq4; 

      // particle
      doublev4 pij1 = mj1 * r_inv1;
      doublev4 pij2 = mj2 * r_inv2;
      doublev4 pij3 = mj3 * r_inv3;
      doublev4 pij4 = mj4 * r_inv4;

      doublev4 mri3_1 = r_inv_sq1*pij1;
      doublev4 mri3_2 = r_inv_sq2*pij2;
      doublev4 mri3_3 = r_inv_sq3*pij3;
      doublev4 mri3_4 = r_inv_sq4*pij4;

      doublev4 xij1 = d1v4 - ps_coff->r_collv4*r_inv1;
      doublev4 xij2 = d1v4 - ps_coff->r_collv4*r_inv2;
      doublev4 xij3 = d1v4 - ps_coff->r_collv4*r_inv3;
      doublev4 xij4 = d1v4 - ps_coff->r_collv4*r_inv4;

      doublev4 fsd1 = ps_coff->k_v_miv4 * xij1 + ps_coff->eta_v_miv4 * rvji1 + ps_coff->m_r_coll_inv_qv4;
      doublev4 fsd2 = ps_coff->k_v_miv4 * xij2 + ps_coff->eta_v_miv4 * rvji2 + ps_coff->m_r_coll_inv_qv4;
      doublev4 fsd3 = ps_coff->k_v_miv4 * xij3 + ps_coff->eta_v_miv4 * rvji3 + ps_coff->m_r_coll_inv_qv4;
      doublev4 fsd4 = ps_coff->k_v_miv4 * xij4 + ps_coff->eta_v_miv4 * rvji4 + ps_coff->m_r_coll_inv_qv4;

      mri3_1 = simd_vsellt(xij1, fsd1, mri3_1);
      mri3_2 = simd_vsellt(xij2, fsd2, mri3_2);
      mri3_3 = simd_vsellt(xij3, fsd3, mri3_3);
      mri3_4 = simd_vsellt(xij4, fsd4, mri3_4);

      doublev4 mdx1 = mri3_1 * dx1;
      doublev4 mdx2 = mri3_2 * dx2;
      doublev4 mdx3 = mri3_3 * dx3;
      doublev4 mdx4 = mri3_4 * dx4;
      
      doublev4 mdy1 = mri3_1 * dy1;
      doublev4 mdy2 = mri3_2 * dy2;
      doublev4 mdy3 = mri3_3 * dy3;
      doublev4 mdy4 = mri3_4 * dy4;
                            
      doublev4 mdz1 = mri3_1 * dz1;
      doublev4 mdz2 = mri3_2 * dz2;
      doublev4 mdz3 = mri3_3 * dz3;
      doublev4 mdz4 = mri3_4 * dz4;

      fx -= mdx1 + mdx2 + mdx3 + mdx4;
      fy -= mdy1 + mdy2 + mdy3 + mdy4;
      fz -= mdz1 + mdz2 + mdz3 + mdz4;
      fp -= pij1 + pij2 + pij3 + pij4;

      //satelite      
      pij1 = mi * r_inv1;
      pij2 = mi * r_inv2;
      pij3 = mi * r_inv3;
      pij4 = mi * r_inv4;

      mri3_1 = r_inv_sq1*pij1;
      mri3_2 = r_inv_sq2*pij2;
      mri3_3 = r_inv_sq3*pij3;
      mri3_4 = r_inv_sq4*pij4;

      xij1 = d1v4 - sp_coff->r_collv4*r_inv1;
      xij2 = d1v4 - sp_coff->r_collv4*r_inv2;
      xij3 = d1v4 - sp_coff->r_collv4*r_inv3;
      xij4 = d1v4 - sp_coff->r_collv4*r_inv4;

      fsd1 = sp_coff->k_v_miv4 * xij1 + sp_coff->eta_v_miv4 * rvji1 + sp_coff->m_r_coll_inv_qv4;
      fsd2 = sp_coff->k_v_miv4 * xij2 + sp_coff->eta_v_miv4 * rvji2 + sp_coff->m_r_coll_inv_qv4;
      fsd3 = sp_coff->k_v_miv4 * xij3 + sp_coff->eta_v_miv4 * rvji3 + sp_coff->m_r_coll_inv_qv4;
      fsd4 = sp_coff->k_v_miv4 * xij4 + sp_coff->eta_v_miv4 * rvji4 + sp_coff->m_r_coll_inv_qv4;
      
      mri3_1 = -simd_vsellt(xij1, fsd1, mri3_1);
      mri3_2 = -simd_vsellt(xij2, fsd2, mri3_2);
      mri3_3 = -simd_vsellt(xij3, fsd3, mri3_3);
      mri3_4 = -simd_vsellt(xij4, fsd4, mri3_4);

      mdx1 = mri3_1 * dx1;
      mdx2 = mri3_2 * dx2;
      mdx3 = mri3_3 * dx3;
      mdx4 = mri3_4 * dx4;

      mdy1 = mri3_1 * dy1;
      mdy2 = mri3_2 * dy2;
      mdy3 = mri3_3 * dy3;
      mdy4 = mri3_4 * dy4;
                   
      mdz1 = mri3_1 * dz1;
      mdz2 = mri3_2 * dz2;
      mdz3 = mri3_3 * dz3;
      mdz4 = mri3_4 * dz4;
      
      doublev4 fs1, fs2, fs3, fs4;

      simd_load(fs1, (double*)&forces[j  ]);
      simd_load(fs2, (double*)&forces[j+1]);
      simd_load(fs3, (double*)&forces[j+2]);
      simd_load(fs4, (double*)&forces[j+3]);

      dv4_trans(&mdx1, &mdy1, &mdz1, &pij1);
      dv4_trans(&mdx2, &mdy2, &mdz2, &pij2);
      dv4_trans(&mdx3, &mdy3, &mdz3, &pij3);
      dv4_trans(&mdx4, &mdy4, &mdz4, &pij4);

      fs1 -= mdx1 + mdy1 + mdz1 + pij1;
      fs2 -= mdx2 + mdy2 + mdz2 + pij2;
      fs3 -= mdx3 + mdy3 + mdz3 + pij3;
      fs4 -= mdx4 + mdy4 + mdz4 + pij4;

      simd_store(fs1, (double*)&forces[j  ]);
      simd_store(fs2, (double*)&forces[j+1]);
      simd_store(fs3, (double*)&forces[j+2]);
      simd_store(fs4, (double*)&forces[j+3]);
#ifdef DEEP_DEBUG
      if(my_id==ID_CHECK) {
        printf("forces[%d]=%e %e %e %e\n", j  ,forces[j  ].ax, forces[j  ].ay, forces[j  ].az, forces[j  ].pot);
        printf("forces[%d]=%e %e %e %e\n", j+1,forces[j+1].ax, forces[j+1].ay, forces[j+1].az, forces[j+1].pot);
        printf("forces[%d]=%e %e %e %e\n", j+2,forces[j+2].ax, forces[j+2].ay, forces[j+2].az, forces[j+2].pot);
        printf("forces[%d]=%e %e %e %e\n", j+3,forces[j+3].ax, forces[j+3].ay, forces[j+3].az, forces[j+3].pot);
      }    
#endif
      //    }

    simd_store(fx,forcei[ii].ax);
    simd_store(fy,forcei[ii].ay);
    simd_store(fz,forcei[ii].az);
    simd_store(fp,forcei[ii].pot);

#ifdef DEEP_DEBUG
    if(my_id==ID_CHECK) {
      printf("forcei[%d].ax=%e %e %e %e\n",ii,forcei[ii].ax[0],forcei[ii].ax[1],forcei[ii].ax[2],forcei[ii].ax[3]);
      printf("forcei[%d].ay=%e %e %e %e\n",ii,forcei[ii].ay[0],forcei[ii].ay[1],forcei[ii].ay[2],forcei[ii].ay[3]);
      printf("forcei[%d].az=%e %e %e %e\n",ii,forcei[ii].az[0],forcei[ii].az[1],forcei[ii].az[2],forcei[ii].az[3]);
      printf("forcei[%d].pot=%e %e %e %e\n",ii,forcei[ii].pot[0],forcei[ii].pot[1],forcei[ii].pot[2],forcei[ii].pot[3]);
    }    
#endif

  }
}
#endif

void ReduceForceSat(const int n_sat, ForceMM* forces) {
    int i;
    for (i=0; i<n_sat; i+=4){
        //* Load
        doublev4 f1,f2,f3,f4;
        simd_load(f1, (double*)&forces[i  ]);
        simd_load(f2, (double*)&forces[i+1]);
        simd_load(f3, (double*)&forces[i+2]);
        simd_load(f4, (double*)&forces[i+3]);
        //* Reduction
        sync_array_();
        reduce(&f1);
        sync_array_();
        reduce(&f2);
        sync_array_();
        reduce(&f3);
        sync_array_();
        reduce(&f4);
        //* Check 
#ifdef SUNWAY_FORCE_KERNEL_DEBUG
        if (my_id==ID_CHECK) {
           printf("s[%d]-s[%d], reduced force:\n",i,j);
           simd_print_doublev4(f1);
           simd_print_doublev4(f2);
           simd_print_doublev4(f3);
           simd_print_doublev4(f4);
        }
#endif
        simd_store(f1,(double*)&forces[i  ]);
        simd_store(f2,(double*)&forces[i+1]);
        simd_store(f3,(double*)&forces[i+2]);
        simd_store(f4,(double*)&forces[i+3]);
    }
}

__attribute__((noinline))
void CalcGrav_ss(const int n_sat, EpiMM* sat, const int nj, const EpiMM* satj, ForceMM* forces, const Force_coff* ss_coff, const doublev4* dt, const doublev4* dtv){
  int i,j;
  const doublev4 vmask={1.0,1.0,1.0,0.0};

  for(i=0; i<n_sat; i++){

    doublev4 ri,vi;
    simd_load(ri, (double*)sat[i].r);
    simd_load(vi, (double*)sat[i].v);

    const doublev4 xi = simd_vshff(ri,ri,0x00);
    const doublev4 yi = simd_vshff(ri,ri,0x55);
    const doublev4 zi = simd_vshff(ri,ri,0xaa);
    // const doublev4 mi = simd_vshff(ri,ri,0xff);

    const doublev4 vxi = simd_vshff(vi,vi,0x00);
    const doublev4 vyi = simd_vshff(vi,vi,0x55);
    const doublev4 vzi = simd_vshff(vi,vi,0xaa);
    const doublev4 idi = simd_vshff(vi,vi,0xff);

    doublev4 fi;
    simd_load(fi, (double*)&forces[i]);

#ifdef SAT_DEBUG
    if(my_id==ID_CHECK) {
      printf("reading sat[%d] r,v,f:\n",i);
      simd_print_doublev4(ri);
      simd_print_doublev4(vi);
      simd_print_doublev4(fi);
    }
    
#endif

    doublev4 mdx = d0v4;
    doublev4 mdy = d0v4;
    doublev4 mdz = d0v4;
    doublev4 pot = d0v4;

    for(j=0; j<nj; j+=4) {
      doublev4 xj,yj,zj,mj;
      //! x,y,z,m
      simd_load(xj,(double*)satj[j  ].r);
      simd_load(yj,(double*)satj[j+1].r);
      simd_load(zj,(double*)satj[j+2].r);
      simd_load(mj,(double*)satj[j+3].r);

#ifdef SAT_DEBUG
      if(my_id==ID_CHECK) {
        printf("load satj[%d-%d]\n",j,j+3);
        simd_print_doublev4(xj);
        simd_print_doublev4(yj);
        simd_print_doublev4(zj);
        simd_print_doublev4(mj);
      }
#endif

      //! vx,vy,vz,id
      doublev4 vxj,vyj,vzj,idj;
      simd_load(vxj,(double*)satj[j  ].v);
      simd_load(vyj,(double*)satj[j+1].v);
      simd_load(vzj,(double*)satj[j+2].v);
      simd_load(idj,(double*)satj[j+3].v);

#ifdef SAT_DEBUG
      if(my_id==ID_CHECK) {
        simd_print_doublev4(vxj);
        simd_print_doublev4(vyj);
        simd_print_doublev4(vzj);
        simd_print_doublev4(idj);
      }
#endif
      
      const doublev4 dx = xi - xj; 
      const doublev4 dy = yi - yj; 
      const doublev4 dz = zi - zj; 
    
      const doublev4 r_sq_org = dx*dx + dy*dy + dz*dz;
      const doublev4 mask = simd_vfcmpeq(idi,idj);

#ifdef SAT_DEBUG
      if(my_id==ID_CHECK) {
        printf("mask i=%d, j=%d\n",i,j);
        simd_print_doublev4(idi);
        simd_print_doublev4(idj);
        simd_print_doublev4(mask);
      }
#endif

      const doublev4 r_sq = simd_vseleq(mask, r_sq_org, dhuge);
      const doublev4 r_inv = d1v4 / simd_vsqrtd(r_sq);
      const doublev4 r_inv_sq = r_inv * r_inv;

      const doublev4 rvxji = (vxi - vxj) * dx;
      const doublev4 rvyji = (vyi - vyj) * dy;
      const doublev4 rvzji = (vzi - vzj) * dz;
      const doublev4 rvji = (rvxji + rvyji + rvzji) * r_inv_sq; 

      const doublev4 pij = mj * r_inv;
      doublev4 mri3 = r_inv_sq*pij;

      const doublev4 xij = d1v4 - ss_coff->r_collv4*r_inv;
      const doublev4 fsd = ss_coff->k_v_miv4 * xij + ss_coff->eta_v_miv4 * rvji + ss_coff->m_r_coll_inv_qv4;

      mri3 = simd_vsellt(xij, fsd, mri3);

      mdx += mri3 * dx;
      mdy += mri3 * dy;
      mdz += mri3 * dz;
      pot += pij;
    }
    
    dv4_trans(&mdx,&mdy,&mdz,&pot);
    fi -= mdx + mdy + mdz + pot;

#ifdef DEEP_DEBUG
    if(my_id==ID_CHECK) {
      printf("s[%d]-s[%d-%d], force:\n",i,j,j+3);
      simd_print_doublev4(mdx);
      simd_print_doublev4(mdy);
      simd_print_doublev4(mdz);
      simd_print_doublev4(pot);
    }
#endif

#ifdef SUNWAY_FORCE_KERNEL_CHECK
    simd_store(fi,(double*)&forces[i]);
#else 
    double id_back=sat[i].v[3];
    vi += fi* *dtv;
    simd_store(vi,sat[i].v);
    sat[i].v[3] = id_back;
    ri += vi* *dt * vmask;
    simd_store(ri,sat[i].r);
#endif
  }
}

#if !defined(SPLIT) && !defined(FULL_ASM)
__attribute__((noinline))
void CalcGrav_sp(const int n_epi, const EpiMM* epi, const SpjMM* spj, ForceMM4* force) {
  int i;

  doublev4 pj1,pj2,pj3,pj4;
  //! x,y,z,m
  simd_load(pj1,(double*)spj[0].r);
  simd_load(pj2,(double*)spj[1].r);
  simd_load(pj3,(double*)spj[2].r);
  simd_load(pj4,(double*)spj[3].r);

  const doublev4 xj1 = simd_vshff(pj1,pj1,0x00);
  const doublev4 yj1 = simd_vshff(pj1,pj1,0x55);
  const doublev4 zj1 = simd_vshff(pj1,pj1,0xaa);
  const doublev4 mj1 = simd_vshff(pj1,pj1,0xff);

  const doublev4 xj2 = simd_vshff(pj2,pj2,0x00);
  const doublev4 yj2 = simd_vshff(pj2,pj2,0x55);
  const doublev4 zj2 = simd_vshff(pj2,pj2,0xaa);
  const doublev4 mj2 = simd_vshff(pj2,pj2,0xff);

  const doublev4 xj3 = simd_vshff(pj3,pj3,0x00);
  const doublev4 yj3 = simd_vshff(pj3,pj3,0x55);
  const doublev4 zj3 = simd_vshff(pj3,pj3,0xaa);
  const doublev4 mj3 = simd_vshff(pj3,pj3,0xff);

  const doublev4 xj4 = simd_vshff(pj4,pj4,0x00);
  const doublev4 yj4 = simd_vshff(pj4,pj4,0x55);
  const doublev4 zj4 = simd_vshff(pj4,pj4,0xaa);
  const doublev4 mj4 = simd_vshff(pj4,pj4,0xff);

  for(i=0; i<n_epi; i+=4){
    doublev4 xi,yi,zi;
    simd_load(xi, (double*)epi[i  ].r);
    simd_load(yi, (double*)epi[i+1].r);
    simd_load(zi, (double*)epi[i+2].r);
      
    doublev4 fx,fy,fz,fp;
    int ii=i>>2;
    simd_load(fx,force[ii].ax);
    simd_load(fy,force[ii].ay);
    simd_load(fz,force[ii].az);
    simd_load(fp,force[ii].pot);
        
    const doublev4 dx1 = xi - xj1; 
    const doublev4 dx2 = xi - xj2; 
    const doublev4 dx3 = xi - xj3; 
    const doublev4 dx4 = xi - xj4; 
                              
    const doublev4 dy1 = yi - yj1; 
    const doublev4 dy2 = yi - yj2; 
    const doublev4 dy3 = yi - yj3; 
    const doublev4 dy4 = yi - yj4; 
                              
    const doublev4 dz1 = zi - zj1; 
    const doublev4 dz2 = zi - zj2; 
    const doublev4 dz3 = zi - zj3; 
    const doublev4 dz4 = zi - zj4; 
    
    const doublev4 r_sq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
    const doublev4 r_sq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    const doublev4 r_sq3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
    const doublev4 r_sq4 = dx4*dx4 + dy4*dy4 + dz4*dz4;

    const doublev4 r_inv1 = 1.0 / simd_vsqrtd(r_sq1);
    const doublev4 r_inv2 = 1.0 / simd_vsqrtd(r_sq2);
    const doublev4 r_inv3 = 1.0 / simd_vsqrtd(r_sq3);
    const doublev4 r_inv4 = 1.0 / simd_vsqrtd(r_sq4);

    const doublev4 pij1 = mj1 * r_inv1;
    const doublev4 pij2 = mj2 * r_inv2;
    const doublev4 pij3 = mj3 * r_inv3;
    const doublev4 pij4 = mj4 * r_inv4;
        
    const doublev4 mri3_1 = r_inv1*r_inv1*pij1;
    const doublev4 mri3_2 = r_inv2*r_inv2*pij2;
    const doublev4 mri3_3 = r_inv3*r_inv3*pij3;
    const doublev4 mri3_4 = r_inv4*r_inv4*pij4;
        
    const doublev4 mdx1 = mri3_1 * dx1;
    const doublev4 mdx2 = mri3_2 * dx2;
    const doublev4 mdx3 = mri3_3 * dx3;
    const doublev4 mdx4 = mri3_4 * dx4;

    const doublev4 mdy1 = mri3_1 * dy1;
    const doublev4 mdy2 = mri3_2 * dy2;
    const doublev4 mdy3 = mri3_3 * dy3;
    const doublev4 mdy4 = mri3_4 * dy4;

    const doublev4 mdz1 = mri3_1 * dz1;
    const doublev4 mdz2 = mri3_2 * dz2;
    const doublev4 mdz3 = mri3_3 * dz3;
    const doublev4 mdz4 = mri3_4 * dz4;
#ifdef DEEP_DEBUG
    if(my_id==ID_CHECK) {
      printf("spj force[%d].fx:\n",ii);
      simd_print_doublev4(mdx1);
      simd_print_doublev4(mdx2);
      simd_print_doublev4(mdx3);
      simd_print_doublev4(mdx4);
    }
#endif
    fx -= mdx1 + mdx2 + mdx3 + mdx4;
    fy -= mdy1 + mdy2 + mdy3 + mdy4;
    fz -= mdz1 + mdz2 + mdz3 + mdz4;
    fp -= pij1 + pij2 + pij3 + pij4;

    simd_store(fx,force[ii].ax);
    simd_store(fy,force[ii].ay);
    simd_store(fz,force[ii].az);
    simd_store(fp,force[ii].pot);
  }

  return;
}
#endif

__attribute__((noinline))
void CalcGrav_planet(const int n_epi, const EpiMM* epi, const doublev4 *mp, ForceMM4* force) {
  int i;

  for(i=0; i<n_epi; i+=4){
    doublev4 xi,yi,zi;
    simd_load(xi, (double*)epi[i  ].r);
    simd_load(yi, (double*)epi[i+1].r);
    simd_load(zi, (double*)epi[i+2].r);
      
    doublev4 fx,fy,fz,fp;
    int ii=i>>2;
    simd_load(fx,force[ii].ax);
    simd_load(fy,force[ii].ay);
    simd_load(fz,force[ii].az);
    simd_load(fp,force[ii].pot);
        
    const doublev4 r_sq = xi*xi + yi*yi + zi*zi;
    const doublev4 r_inv = 1.0 / simd_vsqrtd(r_sq);
    const doublev4 pij = *mp * r_inv;
    const doublev4 mri3 = r_inv*r_inv*pij;
        
    fx -= mri3 * xi;
    fy -= mri3 * yi;
    fz -= mri3 * zi;
    fp -= pij;

    simd_store(fx,force[ii].ax);
    simd_store(fy,force[ii].ay);
    simd_store(fz,force[ii].az);
    simd_store(fp,force[ii].pot);
  }

  return;
}

static inline void Kick(const doublev4* dt_lm, const int ni, EpiMM* epi, ForceMM4* force) {
  int i;
  for (i=0; i<ni; i+=4) {
    doublev4 vx,vy,vz,fx,fy,fz;
    simd_load(vx, (double*)epi[i  ].v);
    simd_load(vy, (double*)epi[i+1].v);
    simd_load(vz, (double*)epi[i+2].v);
    int ii=i>>2;
    simd_load(fx, (double*)force[ii].ax);
    simd_load(fy, (double*)force[ii].ay);
    simd_load(fz, (double*)force[ii].az);

    // kick
    vx += fx * *dt_lm;
    vy += fy * *dt_lm;
    vz += fz * *dt_lm;

    simd_store(vx, (double*)epi[i  ].v);
    simd_store(vy, (double*)epi[i+1].v);
    simd_store(vz, (double*)epi[i+2].v);
  }
}

static inline void Drift(const doublev4* dt_lm, const int ni, EpiMM* epi) {
  int i;
  for (i=0; i<ni; i+=4) {
    doublev4 x,y,z,vx,vy,vz;
    simd_load(x,  (double*)epi[i  ].r);
    simd_load(y,  (double*)epi[i+1].r);
    simd_load(z,  (double*)epi[i+2].r);
    simd_load(vx, (double*)epi[i  ].v);
    simd_load(vy, (double*)epi[i+1].v);
    simd_load(vz, (double*)epi[i+2].v);

    // drift
    x += vx * *dt_lm;
    y += vy * *dt_lm;
    z += vz * *dt_lm;

    simd_store(x,  (double*)epi[i  ].r);
    simd_store(y,  (double*)epi[i+1].r);
    simd_store(z,  (double*)epi[i+2].r);
  }
}

static inline void CalcEnergy(const int ni, EpiMM* epi, ForceMM4* force, double* ekin, double* pot) {
  int i;
  int nloop=ni>>2;
  int noff=ni-(nloop<<2);
  doublev4 ekin4=0;
  doublev4 pot4=0;
  doublev4 vx,vy,vz,m,p;
  for (i=0; i<nloop; i++) {
    int ii=i<<2;
    simd_load(vx, (double*)epi[ii  ].v);
    simd_load(vy, (double*)epi[ii+1].v);
    simd_load(vz, (double*)epi[ii+2].v);
    simd_load(m,  (double*)epi[ii+3].v);
    simd_load(p,  (double*)force[i].pot);

    const doublev4 dv2 = vx*vx + vy*vy + vz*vz;
    ekin4 += m*dv2;
    pot4  += m*p;
  }
  double ek_[4],p_[4];
  simd_store(ekin4,ek_);
  simd_store(pot4, p_ );
  *ekin += ek_[0] + ek_[1] + ek_[2] + ek_[3];
  *pot  +=  p_[0] +  p_[1] +  p_[2] +  p_[3];
  
  if (noff>0) {
    int k = nloop+1;
    int kk = k<<2;
    simd_load(vx, (double*)epi[kk  ].v);
    simd_load(vy, (double*)epi[kk+1].v);
    simd_load(vz, (double*)epi[kk+2].v);
    simd_load(m,  (double*)epi[kk+3].v);
    simd_load(p,  (double*)force[k].pot);

    const doublev4 dv2 = vx*vx + vy*vy + vz*vz;
    ekin4 = m*dv2;
    pot4  = m*p;

    simd_store(ekin4,ek_);
    simd_store(pot4, p_ );
    
    for(i=0; i<noff; i++) {
      *ekin += ek_[i];
      *pot  += p_[i];
    }
  }
}

void ForceKernelSunWay1st(void* pin){
#ifdef PROFILE
  p_count    = 0;
  sp_count   = 0;
  t_count    = 0;
  kd_count   = 0;
  tr_count   = 0;
  init_count = 0;
  mem_count  = 0;
  jlst_count = 0;
  dmaj_count = 0;
  dmas_count = 0;
  geti_count = 0;
  dmaw_count = 0;
  st_count   = 0;
  pl_count   = 0;
  ss_count   = 0;
  gets_count = 0;
  puti_count = 0;
  puts_count = 0;
  rfst_count = 0;
  unsigned long ttmp;
#ifdef SUNWAY_FORCE_KERNEL_CHECK
  ch_count=0;
#endif
  t_count = rtc_();
#endif
  sync_array_();
#ifdef PROFILE
  unsigned long init_sync_count = rtc_()-t_count;
  //if(athread_get_id(-1)==ID_CHECK) printf("int_sync_count position=%lu\n",(long)&init_sync_count);
  ttmp = rtc_();
#endif
  int i, j, k, l;
  for(i=0;i<5;i++) NOP();

  my_id = athread_get_id(-1);
  for(i=0;i<my_id*500;i++) NOP(); // delay each CPE with different cycles to avoid DMA competition.

  NOP();
#ifdef PROFILE
  unsigned long init_delay_count = rtc_()-ttmp;
  ttmp = rtc_();
#endif
  get_col_id_(&my_col);
  get_row_id_(&my_row);

#ifdef FULL_ASM
  int occupy_size = 3232;
  unsigned long* occupy_space = (unsigned long*)(long)ldm_malloc(occupy_size);
  //  for(i=0;i<occupy_size;i++) occupy_space[i]=0;
#ifdef M_CHECK
  if(my_id==ID_CHECK) printf("Allocate occupy region from %lu to %lu\n",(long)occupy_space,(long)occupy_space+occupy_size);
#endif
#endif  

  Force_Kernel_Pars pars;
  volatile unsigned int par_reply = 0;
  athread_get(PE_MODE,pin,&pars,sizeof(Force_Kernel_Pars),(void*)&par_reply,0,0,0);
  while(par_reply!=1);
  //  const double r_coll_lm = pars.r_coll;
  //const double k_lm = pars.kappa;
  //  const double eta_lm = pars.eta;
  //const double mp_lm = pars.m_particle;
  //const double ms_lm = pars.m_satelite;
  //  const double eps2_lm = pars.eps2;
  //  const double dt_lm = pars.dt;
  //  const double dt_lmh = 0.5*dt_lm;
#ifdef PROFILE
  unsigned long init_par_count = rtc_()-ttmp;
  ttmp = rtc_();
#endif

  const doublev4 dt_lm4 = simd_set_doublev4(REP4(pars.dt));
  const doublev4 dt_lmh4= dhv4*dt_lm4;
  //  const doublev4 eps2v4     = simd_set_doublev4(REP4(eps2_lm));
  const doublev4 r_epv4 = simd_set_doublev4(REP4(pars.r_ep));
  const doublev4 r_satv4 = simd_set_doublev4(REP4(pars.r_sat));
  //const doublev4 r_collv4 = simd_set_doublev4(REP4(pars.r_coll));
  //const doublev4 r_coll_inv_qv4 = d1v4/(r_collv4*r_collv4*r_collv4);
  const doublev4 kv4 = simd_set_doublev4(REP4(pars.kappa));
  const doublev4 etav4 = simd_set_doublev4(REP4(pars.eta));

  const doublev4 mpv4 = simd_set_doublev4(REP4(pars.m_particle));
  const doublev4 mp_invv4 = d1v4/mpv4;
  const doublev4 msv4 = simd_set_doublev4(REP4(pars.m_satelite));
  const doublev4 ms_invv4 = d1v4/msv4;

  const doublev4 mplanetv4 = simd_set_doublev4(REP4(pars.m_planet));

  // i particle j particle
  Force_coff pp_coff;
  //  pp_coff.mr_collv4        = r_collv4      *mp_invv4;
  //pp_coff.m_r_coll_inv_qv4 = r_coll_inv_qv4*mpv4;
  pp_coff.r_collv4         = (r_epv4 + r_epv4); 
  pp_coff.m_r_coll_inv_qv4 = mpv4/(pp_coff.r_collv4*pp_coff.r_collv4*pp_coff.r_collv4); 
  //pp_coff.k_v_miv4         = kv4           *mp_invv4;
  pp_coff.k_v_miv4         = (kv4 * dhv4 * mpv4)          *mp_invv4; 
  //pp_coff.eta_v_miv4       = etav4         *mp_invv4;
  pp_coff.eta_v_miv4       = (etav4 * dhv4 * mpv4)        *mp_invv4;
  //pp_coff.r_collv4         = r_collv4;

  const doublev4 r_collspv4 = r_epv4 + r_satv4;
  const doublev4 uspv4 = mpv4 * msv4/(mpv4 + msv4);
  const doublev4 r_coll_inv_spv4 = d1v4/(r_collspv4*r_collspv4*r_collspv4);

#ifdef FULL_ASM
  const doublev4 mp_over_ms = - mpv4 / msv4;
#else
  // i satelite j particle
  Force_coff sp_coff;
  //  sp_coff.mr_collv4        = pp_coff.mr_collv4;
  //sp_coff.m_r_coll_inv_qv4 = pp_coff.m_r_coll_inv_qv4;
  sp_coff.r_collv4         = r_collspv4;
  sp_coff.m_r_coll_inv_qv4 = r_coll_inv_spv4    *mpv4;
  //sp_coff.k_v_miv4         = kv4           *ms_invv4;
  sp_coff.k_v_miv4         = kv4   *uspv4       *ms_invv4;
  //sp_coff.eta_v_miv4       = etav4         *ms_invv4;
  sp_coff.eta_v_miv4       = etav4 *uspv4       *ms_invv4;
#endif

  // i particle j satelite
  Force_coff ps_coff;
  //  ps_coff.mr_collv4        = r_collv4      *ms_invv4;
  ps_coff.r_collv4         = r_collspv4;
  //ps_coff.m_r_coll_inv_qv4 = r_coll_inv_qv4*msv4;
  ps_coff.m_r_coll_inv_qv4 = r_coll_inv_spv4    *msv4;
  //ps_coff.k_v_miv4         = pp_coff.k_v_miv4;
  ps_coff.k_v_miv4         = kv4   *uspv4       *mp_invv4;
  //ps_coff.eta_v_miv4       = pp_coff.eta_v_miv4;
  ps_coff.eta_v_miv4       = etav4 *uspv4       *mp_invv4;


  //  int energy_flag_lm = pars.energy_flag;
  double ekin_lm = 0;
  double pot_lm = 0;

//  double dt_lm4[4]     = {REP4(dt_lm)};
//  double dt_lmh[4]     = {REP4(dt_lmh)};
//  double r_coll_sq4[4] = {REP4(r_coll_sq_lm)};
//  double eps2v4[4]     = {REP4(eps2_lm)};

//  int n_walk_lm = pars.n_walk;
  int n_walk_lm = pars.n_walk_cpe[my_id];
  //  int n_walk_block=n_walk_lm+/64;
//  int col_cut=n_walk_remain%8;
//  int row_cut=(n_walk_remain+7)/8;


#ifdef DEBUG
  if(my_id==ID_CHECK) printf("pid=%d col=%d row=%d CPE starting, nwalk=%d\n",my_id,my_col,my_row,n_walk_lm);
#endif


#ifdef DEBUG
  if(my_id==ID_CHECK)  printf("Memory allocation\nAllocate work buffer\n");
#endif
  
#ifdef PROFILE
  unsigned long init_const_count = rtc_() - ttmp;
  ttmp = rtc_();
#endif

#ifdef SPLIT
  const size_t work_off  = 4;
  const size_t work_size = (NI_LIMIT + 2*work_off) * sizeof(Work);
  pWork work = (pWork)ldm_malloc(work_size) + work_off;

#ifdef M_CHECK
  if(my_id==ID_CHECK){
    printf("pWork work = [%ld:%ld]\n", 
           (long)work, (long)work + work_size - 1 );
  }
#endif                  
#endif

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Allocating epi_lm, size=%d\n",NI_LIMIT);
#endif

#ifndef SPLIT
  const size_t epi_off  = 12;
  EpiMM* epi_lm=(EpiMM*)(long)ldm_malloc((NI_LIMIT+epi_off)*EPI_SIZE);
//  for(i=0;i<epi_off;i++) {
//    simd_store(d0v4,epi_lm[i].r);
//    simd_store(d0v4,epi_lm[i].v);
//  }
  epi_lm += epi_off;
#else  
  //! i particle 
  EpiMM* epi_lm=(EpiMM*)(long)ldm_malloc(NI_LIMIT*EPI_SIZE);
#endif

#ifdef M_CHECK
  if(my_id==ID_CHECK) printf("EPI_ALLOC = %d, position =%ld\n",NI_LIMIT*EPI_SIZE,(long)epi_lm);
#endif

  //! each force structure packs four particle togethers
  ForceMM4* force_lm=(ForceMM4*)(long)ldm_malloc((2+NI_LIMIT/4)*FORCE_SIZE);
  force_lm +=2;

#ifdef M_CHECK              
  if(my_id==ID_CHECK) printf("FORCE_ALLOC = %d, position=%ld\n",(2+NI_LIMIT/4)*FORCE_SIZE,(long)force_lm);
#endif

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Allocating epj_lm and spj_lm\n");
#endif

  //! prepare two blocks of j particles for switching,
  //! j particle spaces are set as same size as i particle
  EpjMM* epj_lm[2];
  epj_lm[0]=(EpjMM*)(long)ldm_malloc(4*EPJ_SIZE);  
  epj_lm[1]=(EpjMM*)(long)ldm_malloc(4*EPJ_SIZE);

#ifdef M_CHECK          
  if(my_id==ID_CHECK) printf("EPJ_ALLOC = %d, position=%ld %ld\n",8*EPJ_SIZE,(long)epj_lm[0],(long)epj_lm[1]);
#endif                  

  //! two groups of sp particles
  SpjMM* spj_lm[2];
  spj_lm[0]=(SpjMM*)(long)ldm_malloc(4*SPJ_SIZE);
  spj_lm[1]=(SpjMM*)(long)ldm_malloc(4*SPJ_SIZE);

#ifdef M_CHECK              
  if(my_id==ID_CHECK) printf("SPJ_ALLOC = %d, position=%ld %ld\n",8*SPJ_SIZE,(long)spj_lm[0],(long)spj_lm[1]);
#endif

  //! two groups of sp particles
  EpiMM* sat_lm[2];
  sat_lm[0]=(EpiMM*)(long)ldm_malloc(4*EPI_SIZE);
  sat_lm[1]=(EpiMM*)(long)ldm_malloc(4*EPI_SIZE);

#ifdef M_CHECK              
  if(my_id==ID_CHECK) printf("SAT_ALLOC = %d, position=%ld %ld\n",8*EPI_SIZE,(long)sat_lm[0],(long)sat_lm[1]);
#endif

  ForceMM* force_sat[2];
  force_sat[0]=(ForceMM*)(long)ldm_malloc(4*FORCEMM_SIZE);
  force_sat[1]=(ForceMM*)(long)ldm_malloc(4*FORCEMM_SIZE);

#ifdef M_CHECK
  if(my_id==ID_CHECK) printf("FORCE_SAT_ALLOC = %d, position=%ld %ld\n",8*FORCEMM_SIZE,(long)force_sat[0],(long)force_sat[1]);
#endif

  //! left bytes (limit - epi - force - epj/spj)
#ifdef SPLIT
  int LDM_REMAIN = LDM_B_LIMIT - NI_LIMIT*(EPI_SIZE+FORCEMM_SIZE)-8*(EPI_SIZE+EPJ_SIZE+SPJ_SIZE+FORCEMM_SIZE)-2*FORCE_SIZE - work_size - 2*FORCEMM_SIZE;
#else
  int LDM_REMAIN = LDM_B_LIMIT - NI_LIMIT*(EPI_SIZE+FORCEMM_SIZE)-8*(EPI_SIZE+EPJ_SIZE+SPJ_SIZE+FORCEMM_SIZE)-2*FORCE_SIZE - epi_off*EPI_SIZE - 2*FORCEMM_SIZE;
#endif

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("LDM remaining size=%d\n",LDM_REMAIN);
#endif

  //! Epj/Spj list 
  int NLST_LIMIT= (LDM_REMAIN/sizeof(int))>>2<<2;   
  
#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Address list allocation, size=%d\n",NLST_LIMIT);
#endif

  int* jlst=(int*)(long)ldm_malloc(NLST_LIMIT*sizeof(int));

#ifdef M_CHECK          
  if(my_id==ID_CHECK) printf("JLST_ALLOC = %lu, position=%ld\n",NLST_LIMIT*sizeof(int),(long)jlst);
#endif                  

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Initializing satellite number\n");
#endif

//satellite number
  // EpiMM* sat_lm;
  // ForceMM* force_sat;
  int n_sat4 = (pars.n_sat+3)>>2<<2;
  int sat_lm_size = n_sat4*sizeof(EpiMM);
  //  if (n_sat4+4>NSAT_LIMIT) {
  //    printf("Error: Satelite number too large (%d), should be less than %d!\n",n_sat4+4,NSAT_LIMIT);
  //    abort();
  //  }
  /*
  sat_lm = (EpiMM*) work;
  force_sat = (ForceMM*) ((long)sat_lm+sat_lm_size);
      
#ifdef M_CHECK
  if(my_id==ID_CHECK) printf("Pointer position, sat_lm=%ld, force_sat=%ld\n",(long)sat_lm,(long)force_sat);
#endif

#ifdef DEBUG
  int force_sat_size = n_sat4*sizeof(ForceMM);
  if(sat_lm_size+force_sat_size>work_size*sizeof(Work)) {
    printf("Error: satellite memory alloc (%d) overflow, maximum size: %d\n",sat_lm_size+force_sat_size,work_size*sizeof(Work));
    abort();
  }
#endif
  //force_sat = (ForceMM*)(long)ldm_malloc(n_sat4*FORCEMM_SIZE);
  */

#ifdef PROFILE
  unsigned long init_alloc_count = rtc_()-ttmp;
  ttmp = rtc_();
#endif
  
  //! for j particle
  dma_desc epj_input;
  volatile unsigned int epj_reply;
  dma_set_op(   &epj_input, DMA_GET);
  dma_set_mode( &epj_input, PE_MODE);
  dma_set_reply(&epj_input, &epj_reply);
  dma_set_size( &epj_input, EPJ_SIZE);  // real j particle size

  //! for sp particle
  dma_desc spj_input;
  volatile unsigned int spj_reply;
  dma_set_op(   &spj_input, DMA_GET);
  dma_set_mode( &spj_input, PE_MODE);
  dma_set_reply(&spj_input, &spj_reply);
  dma_set_size( &spj_input, SPJ_SIZE); 

  //! for index      
  dma_desc pj_index;
  volatile unsigned int index_reply;
  dma_set_op(   &pj_index, DMA_GET);
  dma_set_mode( &pj_index, PE_MODE);
  dma_set_reply(&pj_index, &index_reply);

  //! for sat particles
  dma_desc sat_input;
  volatile unsigned int sat_reply;
  dma_set_op(   &sat_input, DMA_GET);
  dma_set_mode( &sat_input, PE_MODE);
  dma_set_reply(&sat_input, &sat_reply);
  dma_set_size( &sat_input, 4*EPI_SIZE);

  dma_desc force_sat_input;
  volatile unsigned int fin_reply;
  dma_set_op(   &force_sat_input, DMA_GET);
  dma_set_mode( &force_sat_input, PE_MODE);
  dma_set_reply(&force_sat_input, &fin_reply);
  dma_set_size( &force_sat_input, 4*FORCEMM_SIZE);

  dma_desc force_sat_output;
  volatile unsigned int fout_reply;
  dma_set_op(   &force_sat_output, DMA_PUT);
  dma_set_mode( &force_sat_output, PE_MODE);
  dma_set_reply(&force_sat_output, &fout_reply);
  dma_set_size( &force_sat_output, 4*FORCEMM_SIZE);

  
#ifdef PROFILE
  unsigned long init_dma_count = rtc_() - ttmp;
  //  unsigned long I_init_if_count=0, I_init_get_n=0, I_init_I_alloc=0;
#endif

  int first_loop_flag=1;

  int nw_shift = pars.n_disp_walk[my_id];
  // n_walk loop
  for(i=0; i<n_walk_lm; i++){
#ifdef PROFILE
    unsigned long tistart = rtc_();
#endif
//    int ki=i+my_id;
//    if(ki>=n_walk_lm) {
////      int n_epi_block = (n_epi+NI_LIMIT)/NI_LIMIT;
////      for (j=0; j<n_epi_block; j++) {
////        sync_array_();
////        cpe_broadcast(pars.n_sat*8,sat_lm[0].r);
//        if(first_loop_flag) {
// 
//          EpiMM* sat_lm_all = (EpiMM*) work;
//          ForceMM* force_sat_all = (ForceMM*) ((long)sat_lm_all+sat_lm_size);
//          int kk;
//          for(kk=0;kk<n_sat4;kk++) {
//            simd_store(d0v4,(double*)&force_sat_all[kk].ax);
//          }
//          //   memset(force_sat_all,0,n_sat4*FORCEMM_SIZE);
//          first_loop_flag=0;
////        
//          volatile unsigned int fsat_reply=0;
//          athread_put(PE_MODE, force_sat_all, &pars.force_sat[my_id*n_sat4], n_sat4*FORCEMM_SIZE, (void*)&fsat_reply, 0,0);
//          while(fsat_reply<1);
////          sat_reply = 0;
////          athread_get(PE_MODE, pars.sat, sat_lm, pars.n_sat*EPI_SIZE, (void*)&sat_reply, 0,0,0);
////          int n_sat_off = n_sat4-pars.n_sat;
////          if (n_sat_off>0) memset(&sat_lm[pars.n_sat],0,n_sat_off*EPI_SIZE);
////          while(sat_reply<1);
////        }
////        else cpe_broadcast(pars.n_sat*4,(double*)force_sat);
//      }
//        
// 
//      break;
//    }
//#ifdef PROFILE
//    ttmp = rtc_();
//    I_init_if_count += ttmp-tistart;
//#endif
    
    int ki = pars.adr_n_walk[nw_shift+i];
#ifdef BW_DEBUG
    if(my_id==ID_CHECK) printf("CPE_ID=%d staring i_walk=%d\n",my_id,ki);
#endif

    int epi_head = pars.n_disp_epi[ki];
    int epj_head = pars.n_disp_epj[ki];
    int spj_head = pars.n_disp_spj[ki];
    int n_epi = pars.n_epi[ki];
    int n_epj = pars.n_epj[ki];
    int n_spj = pars.n_spj[ki];

    int n_epi_block = n_epi/NI_LIMIT;
    int n_epi_remain = n_epi%NI_LIMIT;
    int ni_shift=0;
    if (n_epi_remain) n_epi_block++;
    else n_epi_remain = NI_LIMIT;
    //int alloc_size=NI_LIMIT;
    //if (n_epi_block==1) alloc_size=(n_epi_remain+3)>>2<<2;
    //if (n_epi_block==1) alloc_size=n_epi_remain;

    //! j particle index list size limit
    int n_epj_block = n_epj/NLST_LIMIT;     //! index list block number, one more for planet
    int n_epj_remain = n_epj%NLST_LIMIT;    //! last index list size, one more fore planet
    if (n_epj_remain) n_epj_block++;
    else n_epj_remain = NLST_LIMIT;

    //! sp particle index list size limit (same as j particle)
    int n_spj_block = n_spj/NLST_LIMIT;    //! index list block number
    int n_spj_remain = n_spj%NLST_LIMIT;
    if (n_spj_remain) n_spj_block++;
    else n_spj_remain = NLST_LIMIT;

#ifdef DEBUG
    if(my_id==ID_CHECK) {
      printf("Current number ni=%d nj=%d ns=%d nsat=%d, staring index istart=%d jstart=%d sstart=%d\n",n_epi,n_epj,n_spj,pars.n_sat,epi_head,epj_head,spj_head);
      printf("blocksize(i,j,sp)=%d %d %d\n",n_epi_block,n_epj_block,n_spj_block);
    }
#endif
#ifdef PROFILE
    unsigned long tiend = rtc_();
    init_count += tiend-tistart;
#endif

    //------------------------EpI-----------------------------//
    for (j=0; j<n_epi_block; j++) {
      //! current number of i particle
      int ni_local = (j+1==n_epi_block)?n_epi_remain:NI_LIMIT;
      int ni_local4 = (ni_local+3)>>2<<2;
#ifdef DEBUG
      if(my_id==ID_CHECK) {
        printf("I block loop index=%d,  ni_local=%d\n",j,ni_local);
        printf("Athread get I particles, staring index=%d size=%d\n",epi_head+ni_shift, ni_local*EPI_SIZE);
      }
#endif
#ifdef PROFILE
      unsigned long tjstart = rtc_();
#endif
#ifdef USE_DMA
      volatile unsigned int epi_reply=0;
      athread_get(PE_MODE, &pars.epi[epi_head+ni_shift], &epi_lm[0], ni_local*EPI_SIZE, (void*)&epi_reply,0,0,0);
#ifdef NO_DOUBLE_BUFFER
      while (epi_reply!=1);
#endif
#endif
#ifdef PROFILE
      unsigned long tjend = rtc_();
      geti_count += tjend-tjstart;
#endif

#ifdef DEBUG
      if(my_id==ID_CHECK) printf("Clearing force_lm memory");
#endif
#ifdef PROFILE
      tjstart = rtc_();
#endif
      //! initialize force
      //      memset(force_lm,0,ni_local4*FORCEMM_SIZE);
      int ni_force=(ni_local4>>2);
      for(k=0; k<ni_force; k++) {
        simd_store(d0v4,force_lm[k].ax);
        simd_store(d0v4,force_lm[k].ay);
        simd_store(d0v4,force_lm[k].az);
        simd_store(d0v4,force_lm[k].pot);
      }
      //! initialize empty i mass to zero
      int ni_off = ni_local4 - ni_local;
      if(ni_off>0) {
        for (k=ni_local; k<ni_local4; k++) 
          simd_store(d0v4,epi_lm[k].r);
      }
#ifdef PROFILE
      tjend = rtc_();
      mem_count += tjend- tjstart;
#endif

#ifdef DEBUG
      if(my_id==ID_CHECK) printf("Wait for athread_get i particle finished\n");
#endif

#ifdef PROFILE
      tjstart = rtc_();
#endif
#ifdef USE_DMA
      while(epi_reply!=1);
#endif
#ifdef PROFILE
      tjend = rtc_();
      geti_count += tjend-tjstart;
#endif

#ifdef DEBUG
      if(my_id==ID_CHECK) {
         printf("I particle data checking, epi_lm[0].x=%e\n",epi_lm[0].r[0]);
         printf("epi_lm transforming\n");
      }
#endif
      //! transformation
#ifdef PROFILE
      tjstart=rtc_();
#endif
      epi_trans(ni_local, epi_lm);
#ifdef PROFILE
      tjend=rtc_();
      tr_count += tjend-tjstart;
#endif

      //-----------------EpJ--------------------//

      //! index memory allocation
      //int nj_alloc_size=NLST_LIMIT;
      //if (n_epj_block==1) {
      //  nj_alloc_size=(n_epj_remain+3)>>2<<2;
      //  if (n_spj_block==1&&n_spj_remain>n_epj_remain) nj_alloc_size=(n_spj_remain+3)>>2<<2;
      //}

      // double buffering switcher
      int sw=0;

      //! loop j index block
      int nj_shift = 0;
      for (k=0; k<n_epj_block; k++) {
        
        int nj_local = NLST_LIMIT;
        int nj_local4 = nj_local;
        int last_loop_flag=(k+1==n_epj_block);

        if (last_loop_flag) {
          nj_local = n_epj_remain; 
          nj_local4 = (nj_local+3)>>2<<2;
          nj_local--; // remove planet 
          //! make sure jlst are filled
          for (l=nj_local; l<nj_local4; l++) jlst[l]=0;
        }
        
#ifdef DEBUG
        if(my_id==ID_CHECK) {
          printf("J block loop index=%d-%d nj_local=%d\n",my_id,k,nj_local, epj_head);
          printf("Dma transferring of j index, staring index=%d, size=%d\n",epj_head+nj_shift,nj_local*sizeof(int));
        }
#endif

        //! get list
#ifdef PROFILE
        unsigned long tkstart = rtc_();
#endif
        index_reply = 0;
        if(nj_local>0) {
          dma_set_size( &pj_index, nj_local*sizeof(int));
          dma(pj_index, (long)&pars.adr_epj[epj_head+nj_shift], (long)jlst);
          dma_wait(&index_reply,1);
        }
#ifdef PROFILE
        unsigned long tkend = rtc_();
        jlst_count += tkend-tkstart;
#endif

#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Staring J particle loop, jlst[0-3]=%d %d %d %d; jlst_mpe[%d]=%d %d %d %d\n",jlst[0],jlst[1],jlst[2],jlst[3],epj_head+nj_shift,pars.adr_epj[epj_head+nj_shift],pars.adr_epj[epj_head+nj_shift+1],pars.adr_epj[epj_head+nj_shift+2],pars.adr_epj[epj_head+nj_shift+3]);
#endif
        // get first 4 j particles

#ifdef PROFILE
        tkstart = rtc_();
#endif
#ifdef USE_DMA
        epj_reply = 0;
        dma(epj_input,(long)&pars.epj[jlst[0]],(long)&epj_lm[sw][0]);
        dma(epj_input,(long)&pars.epj[jlst[1]],(long)&epj_lm[sw][1]);
        dma(epj_input,(long)&pars.epj[jlst[2]],(long)&epj_lm[sw][2]);
        dma(epj_input,(long)&pars.epj[jlst[3]],(long)&epj_lm[sw][3]);
#endif
#ifdef PROFILE
        tkend = rtc_();
        dmaj_count += tkend-tkstart;
#endif

        sw = 1-sw;
#ifdef FULL_ASM
        p_init_ldm_constants((void*)&pp_coff);
#endif

        
#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Dma get first 4 jparticles, data checking epj_mm[3].x=%e, epj_lm[3].x=%e\n",pars.epj[jlst[3]].r[0],epj_lm[1-sw][3].r[0]);
#endif

        // loop j particles
        for (l=4; l<nj_local4; l+=4) {

#ifdef PROFILE
          unsigned long tlstart = rtc_();
          unsigned long pstart = tlstart;
#endif
#ifdef USE_DMA
          while(epj_reply<l);
          //          dma_wait(&ep_reply,l);
#endif
#ifdef PROFILE
          unsigned long pend = rtc_();
#endif
#ifdef USE_DMA
          dma(epj_input,(long)&pars.epj[jlst[l]],  (long)&epj_lm[sw][0]);
          dma(epj_input,(long)&pars.epj[jlst[l+1]],(long)&epj_lm[sw][1]);
          dma(epj_input,(long)&pars.epj[jlst[l+2]],(long)&epj_lm[sw][2]);
          dma(epj_input,(long)&pars.epj[jlst[l+3]],(long)&epj_lm[sw][3]);
#endif
          sw = 1-sw;
#ifdef USE_DMA
#ifdef NO_DOUBLE_BUFFER
          while (epj_reply<l+4);
#endif
#endif
#ifdef PROFILE
          unsigned long tlend = rtc_();
          dmaj_count += tlend - tlstart;
          dmaw_count += pend- pstart;
#endif

          //#ifdef DEBUG
          //          if(my_id==ID_CHECK) {
          //            printf("inner loop l=%d\n",l);
          //            printf("dma get jparticles:\n");
//          int ik;
//          for (ik=0;ik<4;ik++) 
//            if(((EpjMM*)(EP_POINTER[1]))[jlst[l-4+ik]].r[0]!=epj_lm[sw][ik].r[0]) 
//              printf("epmpe[%d].x=%e, epj[%d][%d].x=%e\n",l-4,((EpjMM*)(EP_POINTER[1]))[jlst[l-4]].r[0],sw,0,epj_lm[sw][0].r[0]);
//            printf("epmpe[%d].x=%e, epj[%d][%d].x=%e\n",l-3,((EpjMM*)(EP_POINTER[1]))[jlst[l-3]].r[0],sw,1,epj_lm[sw][1].r[0]);
//            printf("epmpe[%d].x=%e, epj[%d][%d].x=%e\n",l-2,((EpjMM*)(EP_POINTER[1]))[jlst[l-2]].r[0],sw,2,epj_lm[sw][2].r[0]);
//            printf("epmpe[%d].x=%e, epj[%d][%d].x=%e\n",l-1,((EpjMM*)(EP_POINTER[1]))[jlst[l-1]].r[0],sw,3,epj_lm[sw][3].r[0]);
            //          }
          //#endif
          
#ifdef PROFILE
          tlstart=rtc_();
#endif
#ifdef SPLIT
          CalcGrav_p_split2(ni_local, epi_lm, epj_lm[sw], force_lm, &pp_coff, work);
#elif defined(FULL_ASM)
          CalcGrav_p(ni_local, epi_lm, epj_lm[sw], force_lm);
#else
          CalcGrav_p(ni_local, epi_lm, epj_lm[sw], force_lm, &pp_coff);
#endif
#ifdef PROFILE
          tlend=rtc_();
          p_count +=tlend-tlstart;
#endif

#ifdef DEBUG
          if(my_id==ID_CHECK) printf("J particle inner loop=%d finished, data checking: force_lm[0].ax=%e\n",l,force_lm[0].ax[0]);
#endif
        }
        sw = 1-sw;
#ifdef USE_DMA
        while(epj_reply<nj_local4);
#endif
#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Begining last j particle group, nj_local4=%d\n",nj_local4);
#endif
        int nj_edge=4-(nj_local4 - nj_local);
        //        if(nj_edge<4) memset(&epj_lm[sw][nj_edge],0,(4-nj_edge)*EPJ_SIZE);
//        if(last_loop_flag) {
//          epj_lm[sw][nj_edge].r[3] = pars.m_planet;
//          epj_lm[sw][nj_edge].v[3] = -1.0;
//        }
        for (l=nj_edge; l<4; l++) {
          simd_store(dhuge,epj_lm[sw][l].r);
          //epj_lm[sw][l].r[0]=HUGE;
          //epj_lm[sw][l].r[1]=HUGE;          
          //epj_lm[sw][l].r[2]=HUGE;
          epj_lm[sw][l].r[3]=0;
        }
#ifdef PROFILE
        tkstart=rtc_();
#endif
#ifdef SPLIT
        CalcGrav_p_split2(ni_local, epi_lm, epj_lm[sw], force_lm, &pp_coff, work);
#elif defined(FULL_ASM)
        CalcGrav_p(ni_local, epi_lm, epj_lm[sw], force_lm);
#else
        CalcGrav_p(ni_local, epi_lm, epj_lm[sw], force_lm, &pp_coff);
#endif
#ifdef PROFILE
        tkend=rtc_();
        p_count +=tkend-tkstart;
#endif
        nj_shift += nj_local;
      }

      //--------------SpJ------------------//

      //nj_alloc_size=NLST_LIMIT;
      //if (n_spj_block==1) nj_alloc_size=(n_spj_remain+3)>>2<<2;
      //if (n_spj_block==1) nj_alloc_size=n_spj_remain;

      sw = 0;
#if defined(SPLIT) || defined(FULL_ASM)
      sp_init_ldm_constants();
#endif      

      int ns_shift=0;
      for (k=0; k<n_spj_block; k++) {

        int ns_local = NLST_LIMIT;
        int ns_local4 = ns_local;

        if (k+1==n_spj_block) {
          ns_local = n_spj_remain;
          ns_local4 = (ns_local+3)>>2<<2;
          //! make sure jlst are filled
          for (l=ns_local; l<ns_local4; l++) jlst[l]=0;
        }

        //! get list
#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Dma transferring of sp index, staring index=%d, size=%d\n",spj_head+ns_shift,ns_local*sizeof(int));
#endif  

#ifdef PROFILE
        unsigned long tkstart = rtc_();
#endif
        index_reply = 0;
        dma_set_size( &pj_index, ns_local*sizeof(int));
        dma(pj_index, (long)&pars.adr_spj[spj_head+ns_shift], (long)jlst);
        dma_wait(&index_reply,1);
#ifdef PROFILE
        unsigned long tkend = rtc_();
        jlst_count += tkend- tkstart;
#endif

#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Staring sp particle loop, jlst[0-3]=%d %d %d %d\n",jlst[0],jlst[1],jlst[2],jlst[3]);
#endif

        
#ifdef PROFILE
        tkstart = rtc_();
#endif
#ifdef USE_DMA
        spj_reply = 0;
        dma(spj_input,(long)&pars.spj[jlst[0]],(long)&(spj_lm[sw][0]));
        dma(spj_input,(long)&pars.spj[jlst[1]],(long)&(spj_lm[sw][1]));
        dma(spj_input,(long)&pars.spj[jlst[2]],(long)&(spj_lm[sw][2]));
        dma(spj_input,(long)&pars.spj[jlst[3]],(long)&(spj_lm[sw][3]));
#endif
        sw=1-sw;
#ifdef PROFILE
        tkend = rtc_();
        dmas_count += tkend - tkstart;
#endif

#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Dma first get 4 sp particles, data checking: spj_mm[4].x=%e, spj_lm[4].x=%e\n",pars.spj[jlst[3]].r[0],spj_lm[1-sw][3].r[0]);
#endif
        // loop j particles
        for (l=4; l<ns_local; l+=4) {

#ifdef PROFILE
          unsigned long tlstart = rtc_();
#endif
#ifdef USE_DMA
          while(spj_reply<l);
          dma(spj_input,(long)&pars.spj[jlst[l]],(long)&spj_lm[sw][0]);
          dma(spj_input,(long)&pars.spj[jlst[l+1]],(long)&spj_lm[sw][1]);
          dma(spj_input,(long)&pars.spj[jlst[l+2]],(long)&spj_lm[sw][2]);
          dma(spj_input,(long)&pars.spj[jlst[l+3]],(long)&spj_lm[sw][3]);
#endif
          sw=1-sw;
#ifdef USE_DMA
#ifdef NO_DOUBLE_BUFFER
          while (spj_reply<l+4);
#endif
#endif
#ifdef PROFILE
          unsigned long tlend = rtc_();
          dmas_count += tlend - tlstart;
#endif

//#ifdef DEBUG
//          if(my_id==ID_CHECK) printf("inner loop l=%d\n",l);
//          if(my_id==ID_CHECK) printf("dma first get 4 sp particles, sp_p[4].x=%e, spj[4].x=%e\n",pars.spj[jlst[3]].r[0],spj_lm[sw][3].r[0]);
//#endif

#ifdef PROFILE
          tlstart=rtc_();
#endif
#if defined(SPLIT) || defined(FULL_ASM)
          CalcGrav_sp_loop(ni_local, epi_lm, spj_lm[sw], force_lm);
#else
          CalcGrav_sp(ni_local, epi_lm, spj_lm[sw], force_lm);
#endif
#ifdef PROFILE
          tlend=rtc_();
          sp_count +=tlend-tlstart;
#endif
#ifdef DEBUG
          if(my_id==ID_CHECK) printf("SP particle inner loop=%d finished, data checking: force_lm[0].ax=%e\n",l,force_lm[0].ax[0]);
#endif
        }
#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Begining last sp particle group, ns_local=%d ns_local4=%d\n",ns_local,ns_local4);
#endif        
        
        sw=1-sw;
#ifdef USE_DMA
        dma_wait(&spj_reply,ns_local4);
#endif
        int ns_edge=4-(ns_local4 - ns_local);
        for (l=ns_edge; l<4; l++) {
          simd_store(dhuge,spj_lm[sw][l].r);
          spj_lm[sw][l].r[3]=0.0;
        }
        // if(ns_edge<4) 
          //memset(&spj_lm[sw][ns_edge],0,(4-ns_edge)*SPJ_SIZE);

#ifdef PROFILE
        tkstart=rtc_();
#endif
#if defined(SPLIT) || defined(FULL_ASM)
        CalcGrav_sp_loop(ni_local, epi_lm, spj_lm[sw], force_lm);
#else
        CalcGrav_sp(ni_local, epi_lm, spj_lm[sw], force_lm);
#endif
#ifdef PROFILE
        tkend=rtc_();
        sp_count += tkend-tkstart;
#endif
        ns_shift += ns_local;
      }

      // satellite 
#ifdef DEBUG
      if(my_id==ID_CHECK) printf("Starting satellites - i particles forces\n");
#endif

#ifdef PROFILE
      tjstart = rtc_();
#endif
#ifdef USE_DMA
      sat_reply = 0;
      fout_reply = 1;
      dma(sat_input, (long)&pars.sat[0],(long)&sat_lm[sw][0]);
      if(!first_loop_flag) {
        fin_reply = 0;
        dma(force_sat_input, (long)&pars.force_sat[my_id*n_sat4],(long)&force_sat[sw][0]);
      }
      else {
        for(l=0; l<4; l++) simd_store(d0v4,(double*)&force_sat[sw][l].ax);
        //        memset(&force_sat[sw][0],0,4*FORCEMM_SIZE);
      }
#endif

#ifdef PROFILE
      tjend = rtc_();
      gets_count += tjend-tjstart;
#endif

      sw=1-sw;
#ifdef FULL_ASM
      p_init_ldm_constants((void*)&ps_coff);
#endif

#ifdef DEBUG
      if(my_id==ID_CHECK) printf("Dma get first 4 sat particles, data checking sat_mm[3].x=%e, sat_lm[3].x=%e\n",pars.sat[3].r[0],sat_lm[1-sw][3].r[0]);
#endif

      for (l=4; l<n_sat4; l+=4) {
#ifdef PROFILE
        unsigned long tlstart = rtc_();
#endif
        volatile unsigned int lk = (l>>2);
#ifdef USE_DMA
        while(sat_reply<lk);
        dma(sat_input, (long)&pars.sat[l],(long)&sat_lm[sw][0]);
        while(fout_reply<lk);
        if(!first_loop_flag) {
          while(fin_reply<lk);
          dma(force_sat_input, (long)&pars.force_sat[my_id*n_sat4+l],(long)&force_sat[sw][0]);
#ifdef NO_DOUBLE_BUFFER
          while (fin_reply<lk+1);
#endif
        }
        else {
          int lkk;
          for(lkk=0; lkk<4; lkk++) 
            simd_store(d0v4,(double*)&force_sat[sw][lkk].ax);
          //  memset(&force_sat[sw][0],0,4*FORCEMM_SIZE);
        }
#endif
        sw=1-sw;
#ifdef USE_DMA
#ifdef NO_DOUBLE_BUFFER
        while (sat_reply<lk+1);
#endif
#endif
#ifdef PROFILE
        unsigned long tlend = rtc_();
        gets_count += tlend - tlstart;
#endif


#ifdef DEBUG
        if(my_id==ID_CHECK) printf("Inner loop=%d, force calculation, data checking sat_mm[%d].x=%e, sat_lm[3].x=%e\n",lk,l-4,pars.sat[l-4].r[0],sat_lm[sw][0].r[0]);
#endif
        //! for satelite force
#ifdef PROFILE
        tlstart=rtc_();
#endif
#ifdef SPLIT
        CalcGrav_ps_split(ni_local, epi_lm, sat_lm[sw], force_lm, force_sat[sw], &ps_coff, &sp_coff, work);
#elif defined(FULL_ASM)
        clear_satellite_force();
        CalcGrav_ps(ni_local, epi_lm, sat_lm[sw], force_lm);
        collect_sat_force(force_sat[sw],&mp_over_ms);
#else
        CalcGrav_ps(ni_local, epi_lm, sat_lm[sw], force_lm, force_sat[sw], &ps_coff, &sp_coff);
#endif
#ifdef PROFILE
        tlend=rtc_();
        st_count +=tlend-tlstart;
#endif

#ifdef PROFILE
        tlstart = rtc_();
#endif
#ifdef USE_DMA
        while(fout_reply<lk);
        dma(force_sat_output, (long)&pars.force_sat[my_id*n_sat4+l-4],(long)&force_sat[sw][0]);
#ifdef NO_DOUBLE_BUFFER
        while (fout_reply<lk+1);
#endif
#endif
#ifdef PROFILE
        tlend = rtc_();
        puts_count += tlend -tlstart;
#endif

#ifdef DEBUG
        if(my_id==ID_CHECK) printf("SAT particle inner loop=%d finished, data checking: force_lm[0].ax=%e force_sat[0].ax=%e\n",l,force_lm[0].ax[0],force_sat[sw][0].ax);
#endif
      }
#ifdef DEBUG
      if(my_id==ID_CHECK) printf("Begining last sat particle group, n_sat4=%d n_sat=%d\n",n_sat4,pars.n_sat);
#endif        

#ifdef PROFILE
      tjstart=rtc_();
#endif
      sw=1-sw;
      int n_sat_d4 = (n_sat4>>2);
#ifdef USE_DMA
      dma_wait(&sat_reply,n_sat_d4);
#endif
      int nsat_edge = 4-(n_sat4-pars.n_sat);
      for(l=nsat_edge; l<4; l++) {
        simd_store(dhuge,sat_lm[sw][l].r);
        sat_lm[sw][l].r[3]=0.0;
        // simd_store(dhuge,&sat_lm[sw][l].v);
      }
      // if(nsat_edge<4)
        // memset(&sat_lm[sw][nsat_edge],0,(4-nsat_edge)*EPI_SIZE);

#ifdef USE_DMA
      if(!first_loop_flag) dma_wait(&fin_reply,n_sat_d4);
      else first_loop_flag=0;
#endif
#ifdef PROFILE
      tjend=rtc_();
      gets_count += tjend - tjstart;
#endif      
      
      //! for satelite force
#ifdef PROFILE
      tjstart=rtc_();
#endif
#ifdef SPLIT
      CalcGrav_ps_split(ni_local, epi_lm, sat_lm[sw], force_lm, force_sat[sw], &ps_coff, &sp_coff, work);
#elif defined(FULL_ASM)
      clear_satellite_force();
      CalcGrav_ps(ni_local, epi_lm, sat_lm[sw], force_lm);
      collect_sat_force(force_sat[sw],&mp_over_ms);
#else
      CalcGrav_ps(ni_local, epi_lm, sat_lm[sw], force_lm, force_sat[sw], &ps_coff, &sp_coff);
#endif
#ifdef PROFILE
      tjend=rtc_();
      st_count +=tjend-tjstart;
#endif

#ifdef PROFILE
      tjstart=rtc_();
#endif
#ifdef USE_DMA
      while(fout_reply<n_sat_d4); 
      //if(my_id==ID_CHECK)  printf("fout_reply=%d n_sat_d4=%d\n",fout_reply,n_sat_d4);
      dma(force_sat_output, (long)&pars.force_sat[(my_id+1)*n_sat4-4],(long)&force_sat[sw][0]);
      while(fout_reply<=n_sat_d4);
#endif
#ifdef PROFILE
      tjend =rtc_();
      puts_count +=tjend -tjstart;
#endif
      
      //! for planet force
#ifdef PROFILE
      tjstart=rtc_();
#endif
      CalcGrav_planet(ni_local,epi_lm,&mplanetv4,force_lm);
#ifdef PROFILE
      tjend=rtc_();
      pl_count +=tjend-tjstart;
#endif

//#ifdef DEBUG
//      if(my_id==ID_CHECK) {
//        for(k=0;k<n_sat4;k++) {
//          //          printf("sat[%d].r,v=%e %e %e %e %e %e %e %e\n",k,sat_lm[k].r[0],sat_lm[k].r[1],sat_lm[k].r[2],sat_lm[k].r[3],sat_lm[k].v[0],sat_lm[k].v[1],sat_lm[k].v[2],sat_lm[k].v[3]);
//          printf("force_sat[%d]=%e %e %e %e\n",k,force_sat[k].ax,force_sat[k].ay,force_sat[k].az,force_sat[k].pot);
//        }
//      }
//#endif
      
#ifdef SUNWAY_FORCE_KERNEL_CHECK
#ifdef PROFILE
      tjstart = rtc_();
#endif
      int nforce=(ni_local4>>2);
      for (k=0;k<nforce;k++) {
        doublev4 fx,fy,fz,fp;
        simd_load(fx, force_lm[k].ax);
        simd_load(fy, force_lm[k].ay);
        simd_load(fz, force_lm[k].az);
        simd_load(fp, force_lm[k].pot);
        
        dv4_trans(&fx,&fy,&fz,&fp);
        simd_store(fx, force_lm[k].ax);
        simd_store(fy, force_lm[k].ay);
        simd_store(fz, force_lm[k].az);
        simd_store(fp, force_lm[k].pot);
        
      }
#ifdef DEBUG
      if(my_id==ID_CHECK) {
        printf("Sending i particle force to MPE\nForce check before sending:\n");
        for(l=0;l<4;l++) printf("force_lm[%d].ax=%e\n",0,force_lm[0].ax[l]);
      }
#endif
      volatile unsigned int fput_reply = 0;
      athread_put(PE_MODE, &force_lm[0], &pars.forcei[epi_head+ni_shift], ni_local*FORCEMM_SIZE, (void*)&fput_reply,0,0);
      while(fput_reply!=1);
#ifdef DEBUG
      if(my_id==ID_CHECK) printf("Athread_put i particle force, start index=%d, size=%d\n",epi_head+ni_shift,ni_local*FORCEMM_SIZE);
#endif
#ifdef PROFILE
      tjend = rtc_();
      ch_count += tjend-tjstart;
#endif
#else

      //! drift and kick
#ifdef PROFILE
      tjstart=rtc_();
#endif
      //* [Original implementation by Long]
      if (pars.energy_flag) {
        if (!pars.first_flag) Kick (&dt_lmh4, ni_local, epi_lm, force_lm);
        CalcEnergy(ni_local, epi_lm, force_lm, &ekin_lm, &pot_lm);
        pars.ekin = ekin_lm;
        pars.pot = pot_lm;
        Kick (&dt_lmh4, ni_local, epi_lm, force_lm);
      }
      else if (pars.first_flag) Kick(&dt_lmh4, ni_local, epi_lm, force_lm);
      else Kick(&dt_lm4, ni_local, epi_lm, force_lm);
      if (pars.last_flag) Drift(&dt_lmh4, ni_local, epi_lm);
      else Drift(&dt_lm4, ni_local, epi_lm);
#ifdef PROFILE
      tjend=rtc_();
      kd_count +=tjend-tjstart;
#endif

      //! transformation
#ifdef PROFILE
      tjstart=rtc_();
#endif
      epi_trans(ni_local, epi_lm);
#ifdef PROFILE
      tjend=rtc_();
      tr_count +=tjend-tjstart;
      tjstart=tjend;
#endif

#ifdef USE_DMA
      volatile unsigned int fput_reply = 0;
      // if(my_id==ID_CHECK) printf("pid=%d force[0]_lm=%e\n",my_id,force_lm[0].ax);
      athread_put(PE_MODE, &epi_lm[0], &pars.epi[epi_head+ni_shift], ni_local*EPI_SIZE, (void*)&fput_reply,0,0);
      while(fput_reply!=1);
#endif
#endif
#ifdef PROFILE
      tjend =rtc_();
      puti_count +=tjend -tjstart;
#endif

      /*
      // sending satellite
      if(i+64<n_walk_lm) {
#ifdef PROFILE
        tjstart = rtc_();
#endif
#ifdef DEBUG
        if(my_id==ID_CHECK) {
          printf("Sending satellite force to MPE, offset=%d, checking:\n",my_id*pars.n_sat);
          for(l=0;l<4;l++) printf("force_sat[%d].ax=%e\n",1,force_sat[l].ax);
        }
#endif
        ep_reply=0;
        athread_put(PE_MODE, force_sat, &pars.force_sat[my_id*pars.n_sat], pars.n_sat*FORCEMM_SIZE, (void*)&ep_reply, 0,0);
        while(ep_reply<1);
#ifdef PROFILE
        tjend = rtc_();
        puts_count += tjend -tjstart;
#endif
      }
      */
      
      ni_shift += ni_local;
      // if(my_id==ID_CHECK) printf("pid=%d fin i=%d ni_shift=%d\n",my_id,i,ni_shift);
    }

    // if(my_id==ID_CHECK) printf("pid=%d free size=%d",my_id,alloc_size);
  }

#ifdef PROFILE
  ttmp = rtc_();
  sync_array_();
  final_sync_count = rtc_()-ttmp;
  final_max_count = final_sync_count;
  final_min_count = final_sync_count;
  final_tot_count = final_sync_count;
  max_ulong(&final_max_count);
  sync_array_();
  min_ulong(&final_min_count);
  sync_array_();
  reduce_ulong(&final_tot_count);
#endif

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Free memory allocation\n");
#endif

#ifdef FULL_ASM
  ldm_free(occupy_space,occupy_size);
#endif

#ifdef SPLIT
  ldm_free(work-work_off, work_size);
  ldm_free(epi_lm,NI_LIMIT*EPI_SIZE);
#else
  ldm_free(epi_lm-epi_off,(NI_LIMIT+epi_off)*EPI_SIZE);
#endif

  ldm_free(force_lm-2,(2+NI_LIMIT/4)*FORCE_SIZE);

  ldm_free(epj_lm[0],4*EPJ_SIZE);
  ldm_free(epj_lm[1],4*EPJ_SIZE);

  ldm_free(spj_lm[0],4*SPJ_SIZE);
  ldm_free(spj_lm[1],4*SPJ_SIZE);

  ldm_free(sat_lm[0],4*EPI_SIZE);
  ldm_free(sat_lm[1],4*EPI_SIZE);
    
  ldm_free(force_sat[0],4*FORCEMM_SIZE);
  ldm_free(force_sat[1],4*FORCEMM_SIZE);

  ldm_free(jlst,NLST_LIMIT*sizeof(int));

#ifdef M_CHECK          
    if(my_id==ID_CHECK) {
#ifdef FULL_ASM
      printf("occupy free = %d\n",occupy_size);
#endif
#ifdef SPLIT
      printf("WORK_FREE = %d\n",work_size);
#endif
      printf("EPI_FREE = %d\n",NI_LIMIT*EPI_SIZE);
      printf("FORCE_FREE = %d\n",(2+NI_LIMIT/4)*FORCE_SIZE);
      printf("EPJ_FREE = %d\n",8*EPJ_SIZE);
      printf("SPJ_FREE = %d\n",8*SPJ_SIZE);
      printf("FORCE_SAT_FREE = %d\n",8*FORCEMM_SIZE);
      printf("JLST_FREE = %lu\n",NLST_LIMIT*sizeof(int));
    }
#endif                      

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Reducing satellite force from i particles\n");
#endif

  if(n_walk_lm==0) {
#ifdef SPLIT
    EpiMM* sat_lm_all = (EpiMM*) work;
    ForceMM* force_sat_all = (ForceMM*) ((long)sat_lm_all+sat_lm_size);
#else
    EpiMM* sat_lm_all = (EpiMM*)(long)ldm_malloc(sat_lm_size);
    ForceMM* force_sat_all = (ForceMM*)(long)ldm_malloc(n_sat4*FORCEMM_SIZE);
#endif
    int kk;
//#ifdef M_CHECK
//    if(my_id==ID_CHECK) printf("force_sat_all position = %lu, size=%d\n",(long)force_sat_all,n_sat4*FORCEMM_SIZE);
//#endif

    for(kk=0;kk<n_sat4;kk++) {
      simd_store(d0v4,(double*)&force_sat_all[kk].ax);
    }
    //   memset(force_sat_all,0,n_sat4*FORCEMM_SIZE);
    //   first_loop_flag=0;
    //        
    volatile unsigned int fsat_reply=0;
    athread_put(PE_MODE, force_sat_all, &pars.force_sat[my_id*n_sat4], n_sat4*FORCEMM_SIZE, (void*)&fsat_reply, 0,0);
    while(fsat_reply<1);

#ifndef SPLIT
    ldm_free(sat_lm_all, sat_lm_size);
    ldm_free(force_sat_all, n_sat4*FORCEMM_SIZE);
#endif
  }

  ForceMM* force_sat_all = (ForceMM*)(long)ldm_malloc(n_sat4*FORCEMM_SIZE);

#ifdef M_CHECK          
  if(my_id==ID_CHECK){
    printf("force_sat_all allocation. Position= [%ld:%ld]\n", (long)force_sat_all, (long)(force_sat_all+n_sat4));
  }
#endif                      

#ifdef PROFILE
  unsigned long tstart = rtc_();
#endif
//  ForceMM* force_sat_all = (ForceMM*) ((long)work+sat_lm_size);
  volatile unsigned int fsat_reply=0;
  athread_get(PE_MODE, &pars.force_sat[my_id*n_sat4], force_sat_all, pars.n_sat*FORCEMM_SIZE, (void*)&fsat_reply, 0,0,0);
  while(fsat_reply<1);
#ifdef PROFILE
  unsigned long tend = rtc_();
  gets_count += tend-tstart;
#endif

#ifdef PROFILE
  tstart = rtc_();
#endif
  ReduceForceSat(pars.n_sat, force_sat_all);
#ifdef PROFILE
  tend = rtc_();
  rfst_count += tend-tstart;
#endif

  if(my_id==63) {
#ifdef PROFILE
    tstart = rtc_();
#endif
#ifdef DEBUG
    if(my_id==ID_CHECK) printf("Sending satellite force from i particles to MPE\n");
#endif
#ifdef USE_DMA
    fsat_reply=0;
    athread_put(PE_MODE, force_sat_all, pars.force_sat, pars.n_sat*FORCEMM_SIZE, (void*)&fsat_reply, 0,0);
    while(fsat_reply!=1);
#endif
#ifdef PROFILE
    tend = rtc_();
    puts_count += tend - tstart;
#endif
  }
  ldm_free(force_sat_all,n_sat4*FORCEMM_SIZE);
#ifdef M_CHECK
  if(my_id==ID_CHECK) printf("force_sat_all_FREE = %d\n",n_sat4*FORCEMM_SIZE);
#endif

  sync_array_();
#ifdef PROFILE
  t_count = rtc_() - t_count;
  
#ifdef PROF_PRINT
  if(my_id==ID_CHECK) {
    printf("Initialization kernel cycles:\n");
    printf("Init-sync : %ld\n",init_sync_count);
    printf("Init-delay: %ld\n",init_delay_count);
    printf("Init-par:   %ld\n",init_par_count);
    printf("Init-const: %ld\n",init_const_count);
    printf("Init-alloc:  %ld\n",init_alloc_count);    
    printf("Init-dma:   %ld\n",init_dma_count);    
//    printf("I-Init-if:     %ld\n",I_init_if_count);
//    printf("I-Init-get_n:  %ld\n",I_init_get_n);
//    printf("I-Init-Ialloc: %ld\n",I_init_I_alloc);
  }
#endif
#endif
}


void ForceKernelSunWay2nd(void* pin) {
#ifdef PROFILE
  unsigned long ststart=rtc_();
#endif

  Force_Kernel_Pars pars;
  volatile unsigned int par_reply = 0;
  athread_get(PE_MODE,pin,&pars,sizeof(Force_Kernel_Pars),(void*)&par_reply,0,0,0);
  while(par_reply!=1);

  const doublev4 r_satv4 = simd_set_doublev4(REP4(pars.r_sat));
  //const doublev4 r_collv4 = simd_set_doublev4(REP4(pars.r_coll));
  //const doublev4 r_coll_inv_qv4 = d1v4/(r_collv4*r_collv4*r_collv4);
  const doublev4 kv4 = simd_set_doublev4(REP4(pars.kappa));
  const doublev4 etav4 = simd_set_doublev4(REP4(pars.eta));
  const doublev4 msv4 = simd_set_doublev4(REP4(pars.m_satelite));
  const doublev4 ms_invv4 = d1v4/msv4;
  
  // i satelite j satelite
  Force_coff ss_coff;
  //  ss_coff.mr_collv4        = r_collv4      *ms_invv4;
  ss_coff.r_collv4         = r_satv4 + r_satv4;
  //ss_coff.m_r_coll_inv_qv4 = r_coll_inv_qv4*msv4;
  ss_coff.m_r_coll_inv_qv4 = msv4/(ss_coff.r_collv4*ss_coff.r_collv4*ss_coff.r_collv4);
  //ss_coff.k_v_miv4         = kv4    *ms_invv4;
  ss_coff.k_v_miv4         = (kv4 * dhv4 * msv4)   *ms_invv4;
  //ss_coff.eta_v_miv4       = etav4  *ms_invv4;
  ss_coff.eta_v_miv4       = (etav4 * dhv4 * msv4) *ms_invv4;
  
  
  //satelite self-force
  int n_sat_part = pars.n_sat/64;
  int n_sat_off = pars.n_sat%64;
  int n_sat_start = my_id*n_sat_part;
  int n_sat = n_sat_part;
  if (my_id<n_sat_off) {
    n_sat_start +=my_id;
    n_sat++;
  }
  else n_sat_start += n_sat_off;

  int n_sat_t4 = ((pars.n_sat+4)>>2<<2);
  EpiMM* sat_lm=(EpiMM*)(long)ldm_malloc(n_sat_t4*EPI_SIZE);
  ForceMM* force_sat=(ForceMM*)(long)ldm_malloc(n_sat*FORCEMM_SIZE);
  EpiMM* sati=(EpiMM*)(long)ldm_malloc(n_sat*EPI_SIZE);

#ifdef M_CHECK
  if(my_id==ID_CHECK) printf("Allocate sat_lm, size=%d, position=%ld\n",n_sat*EPI_SIZE,(long)sat_lm);
  if(my_id==ID_CHECK) printf("Allocate force_sat, size=%d, position=%ld\n",n_sat*FORCEMM_SIZE,(long)force_sat);
  if(my_id==ID_CHECK) printf("Allocate sati, size=%d, position=%ld\n",n_sat*EPI_SIZE,(long)sati);
#endif

  volatile unsigned int sat_reply = 0;
  athread_get(PE_MODE, pars.sat, sat_lm, pars.n_sat*EPI_SIZE, (void*)&sat_reply, 0,0,0);
  
  //int n_sat_toff = n_sat_t4-pars.n_sat;
  //if (n_sat_toff>0) {
  int k;
  for(k=pars.n_sat+1;k<n_sat_t4;k++) {
    simd_store(dhuge,sat_lm[k].r);
    sat_lm[k].r[3]=0.0;
      //      simd_store(dhuge,sat_lm[k].v);
  }
    // memset(&sat_lm[pars.n_sat],0,n_sat_toff*EPI_SIZE);
  //  }
  simd_store(d0v4,sat_lm[pars.n_sat].r);
  simd_store(d0v4,sat_lm[pars.n_sat].v);
  sat_lm[pars.n_sat].r[3]=pars.m_planet;
  while(sat_reply<1);
  
#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Athread_get satellite particles, n_sat=%d\ndata check, sat_lm[0]=%e %e %e %e %e %e %e %e\n",pars.n_sat,sat_lm[0].r[0],sat_lm[0].r[1],sat_lm[0].r[2],sat_lm[0].r[3],sat_lm[0].v[0],sat_lm[0].v[1],sat_lm[0].v[2],sat_lm[0].v[3]);
#endif


#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Athread_get satellite forces, n_sat=%d, start index=%d\n",n_sat,n_sat_start);
#endif

  volatile unsigned int fin_reply=0;
  athread_get(PE_MODE, &pars.force_sat[n_sat_start], force_sat, n_sat*FORCEMM_SIZE, (void*)&fin_reply, 0,0,0);
  while(fin_reply<1);

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Memory copy sati and Transform satj\n");
  if(my_id==ID_CHECK) printf("dat check: sat_lm[0]=%e %e %e %e %e %e %e %e\n",sat_lm[0].r[0],sat_lm[0].r[1],sat_lm[0].r[2],sat_lm[0].r[3],sat_lm[0].v[0],sat_lm[0].v[1],sat_lm[0].v[2],sat_lm[0].v[3]);
#endif

  int i;
  for(i=0;i<n_sat;i++) {
    *(doublev4*)sati[i].r = *(doublev4*)sat_lm[n_sat_start+i].r;
    *(doublev4*)sati[i].v = *(doublev4*)sat_lm[n_sat_start+i].v;
  }
  //  memcpy(sati,&sat_lm[n_sat_start],n_sat*EPI_SIZE);
  epi_trans(pars.n_sat+1, sat_lm); // one more for planet
  
#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Satellite - satellite force, data check: sati.r=%e %e %e %e, force_sat[start].ax=%e\n",sati[0].r[0],sati[0].r[1],sati[0].r[2],sati[0].r[3],force_sat[0].ax);
#endif

  doublev4 dtsat = simd_set_doublev4(REP4(pars.dt));
  doublev4 dtvsat = dtsat;
  if(pars.first_flag) dtvsat *= dhv4;
  if(pars.last_flag)  dtsat  *= dhv4;
#ifdef PROFILE
  unsigned long pstart=rtc_();
#endif
  CalcGrav_ss(n_sat, sati, pars.n_sat+1, sat_lm, force_sat, &ss_coff,&dtsat,&dtvsat);
#ifdef PROFILE
  unsigned long pend = rtc_();
  ss_count +=pend-pstart;
#endif

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Calculate satellite force from planet\n");
#endif

#ifdef SUNWAY_FORCE_KERNEL_CHECK
#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Sending satellite force to MPE, data checking. force[0].ax=%e\n",force_sat[n_sat_start].ax);
#endif

  volatile unsigned int fout_reply=0;
  athread_put(PE_MODE, force_sat, &pars.force_sat[n_sat_start], n_sat*FORCEMM_SIZE, (void*)&fout_reply, 0,0);
  while(fout_reply!=1);
#else

#ifdef DEBUG
  if(my_id==ID_CHECK) printf("Sending satellite particle to MPE, data checking. sat[0].x=%e\n",sati[0].r[0]);
#endif
  sat_reply=0;
  athread_put(PE_MODE, sati, &pars.sat[n_sat_start], n_sat*EPI_SIZE, (void*)&sat_reply, 0,0);
  while(sat_reply!=1);  

#endif

  ldm_free(sat_lm, n_sat_t4*EPI_SIZE);
  ldm_free(sati, n_sat*EPI_SIZE);
  ldm_free(force_sat, n_sat*FORCEMM_SIZE);

#ifdef M_CHECK
  if(my_id==ID_CHECK){
      printf("SAT_FREE = %d\n",n_sat_t4*EPI_SIZE);
      printf("SATI_FREE = %d\n",n_sat*EPI_SIZE);
      printf("FORCE_SAT_FREE = %d\n",n_sat*FORCEMM_SIZE);
  }
#endif

#ifdef PROFILE
  sync_array_();
  unsigned long stend=rtc_();  
  t2_count = stend- ststart;
  t_count += t2_count;

  unsigned long dma_total = geti_count + jlst_count + dmaj_count + dmas_count + gets_count + puti_count + puts_count + rfst_count;
  unsigned long dma_max_total = dma_total;
  sync_array_();
  max_ulong(&dma_max_total);
  unsigned long dma_min_total = dma_total;
  sync_array_();
  min_ulong(&dma_min_total);
  unsigned long dma_tot_total = dma_total;
  sync_array_();
  reduce_ulong(&dma_tot_total);

#ifdef PROF_PRINT
  if(my_id==ID_CHECK) {
    // double ptime =  (double)p_count/1.5e9;
    // double sptime = (double)sp_count/1.5e9;
    // double ttime =  (double)t_count/1.5e9;
    // double kdtime = (double)kd_count/1.5e9;
    // double trtime = (double)tr_count/1.5e9;
    printf("PID=%d CPE profiling:\n\ttime[sec]\tcycles\n",my_id);
    printf("I-J    \t%e\t%lu\n", (double)   p_count/1.45e9,   p_count); // I-J force
    printf("I-SP   \t%e\t%lu\n", (double)  sp_count/1.45e9,  sp_count); // I-SP force
    printf("I-SAT  \t%e\t%lu\n", (double)  st_count/1.45e9,  st_count); // I-SAT force
    printf("I-PT   \t%e\t%lu\n", (double)  pl_count/1.45e9,  pl_count); // I-Planet force
    printf("SAT-SAT\t%e\t%lu\n", (double)  ss_count/1.45e9,  ss_count); // I-SAT force
    printf("I-K/D  \t%e\t%lu\n", (double)  kd_count/1.45e9,  kd_count); // I kick/drift
    printf("I-Trans\t%e\t%lu\n", (double)  tr_count/1.45e9,  tr_count); // I transformation
    printf("I-init \t%e\t%lu\n", (double)init_count/1.45e9,init_count); // before I block loop
    printf("IF-init\t%e\t%lu\n", (double) mem_count/1.45e9, mem_count); // I force initial
    printf("GET-i  \t%e\t%lu\n", (double)geti_count/1.45e9,geti_count); // DMA i
    printf("GET-lst\t%e\t%lu\n", (double)jlst_count/1.45e9,jlst_count); // DMA Jlst
    printf("GET-j  \t%e\t%lu\n", (double)dmaj_count/1.45e9,dmaj_count); // DMA J
    printf("GET-jw \t%e\t%lu\n", (double)dmaw_count/1.45e9,dmaw_count); // DMA J wait
    printf("GET-sp \t%e\t%lu\n", (double)dmas_count/1.45e9,dmas_count); // DMA sp
    printf("GET-sat\t%e\t%lu\n", (double)gets_count/1.45e9,gets_count); // DMA satellite
    printf("PUT-i  \t%e\t%lu\n", (double)puti_count/1.45e9,puti_count); // DMA put i
    printf("PUT-fst\t%e\t%lu\n", (double)puts_count/1.45e9,puts_count); // DMA put sat force
    printf("RD-fst \t%e\t%lu\n", (double)rfst_count/1.45e9,rfst_count); // reduce sat force
    
    printf("Tot    \t%e\t%lu\n", (double)   t_count/1.45e9,   t_count);
    printf("Tot    \t%e\t%lu\n", (double)  t2_count/1.45e9,  t2_count);
    printf("DMA    \t%e\t%lu\n", (double) dma_total/1.45e9, dma_total);
    printf("DMAmin \t%e\t%lu\n", (double) dma_min_total/1.45e9, dma_min_total);
    printf("DMAmax \t%e\t%lu\n", (double) dma_max_total/1.45e9, dma_max_total);
    printf("DMAtot \t%e\t%lu\n", (double) dma_tot_total/1.45e9, dma_tot_total);

#ifdef SUNWAY_FORCE_KERNEL_CHECK
    printf("CHF \t%e\t%lu\n", (double)  ch_count/1.45e9,  ch_count);
#endif

    printf("sync_wait_cpe[%d]\t%e\t%lu\n", my_id, (double)final_sync_count/1.45e9, final_sync_count);
    printf("sync_wait_max    \t%e\t%lu\n", (double)final_max_count/1.45e9, final_max_count);
    printf("sync_wait_min    \t%e\t%lu\n", (double)final_min_count/1.45e9, final_min_count);
    printf("sync_wait_total  \t%e\t%lu\n", (double)final_tot_count/1.45e9, final_tot_count);

  }
#endif
#endif
}
