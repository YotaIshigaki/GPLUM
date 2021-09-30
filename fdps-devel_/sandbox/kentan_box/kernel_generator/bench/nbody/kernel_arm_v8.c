#include <stdio.h>
#include <math.h>

#include "user_defined_class.h"
#include "constant.h"

#include <arm_sve.h>
#define CHECK_FORCE
//#define CHECK_FORCE_DETAIL
//#define DETAILED_PROFILE

#ifdef PROFILE
#include "fj_tool/fipp.h"
#endif

#if 0
void transpose4x4_f32(svfloat32_t* x,
		      svfloat32_t* y,
		      svfloat32_t* z,
		      svfloat32_t* w){
  svfloat32_t xy0 = svtrn1_f32(*x,*y);
  svfloat32_t xy1 = svtrn2_f32(*x,*y);
  svfloat32_t zw0 = svtrn1_f32(*z,*w);
  svfloat32_t zw1 = svtrn2_f32(*z,*w);
  *x = svtrn1_f64(xy0,zw0);
  *y = svtrn1_f64(xy1,zw1);
  *z = svtrn2_f64(xy0,zw0);
  *w = svtrn2_f64(xy1,zw1);
}

void gather_st4_f32(svfloat32_t* x,
		    svfloat32_t* y,
		    svfloat32_t* z,
		    svfloat32_t* w){
  const unsigned int tmp[16] = {0,4,8,12,
				1,5,9,13,
				2,6,10,14,
				3,7,11,15};
  const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
  *x = svtbl_f32(*x,index);
  *y = svtbl_f32(*y,index);
  *z = svtbl_f32(*z,index);
  *w = svtbl_f32(*w,index);
  transpose4x4_f32(x,y,z,w);
}
#endif

#define DEBUG_PRINT(x,I,J) float tmp[16];svst1_f32(pg,tmp,x);if(((I)/16)==i && (jjj)==(J)) printf("%d %d: %e\n",I,J,tmp[I%16]);
#define DEBUG_PRINT2(x) float tmp[16];svst1_f32(pg,tmp,x);for(int ii=0;ii<16;ii++) printf("%d %d: %e\n",i+ii,j+jj,tmp[ii]);

//#define NJ_MAX 1024
//#define RECALC_DX

// L1$ size = 64 kB = 16 k single words = 1 k loop for 1 word

void CalcForceLongEpEp(const EPI* __restrict__ epi, const int ni,
		       const EPJ* __restrict__ epj, const int nj,
		       Force* __restrict__ force){
#ifdef DETAILED_PROFILE
  double t0 = 0.0, t1 = 0.0;
#endif

#ifdef PROFILE
  fipp_start();
#endif

  svbool_t pg = svptrue_b32();
  const float r_out2 = r_out*r_out;
  //printf("--- optimized epep kernel ---\n");

  for(int i=0; i<ni; i+=16){
    const svfloat32x3_t xi = svld3_f32(pg,(float*)&epi[i]);
    //svfloat32x4_t ai = {svdup_n_f32(0.f),svdup_n_f32(1.f),svdup_n_f32(2.f),svdup_n_f32(3.f)};
    svfloat32_t ax = svdup_n_f32(0.f);
    svfloat32_t ay = svdup_n_f32(0.f);
    svfloat32_t az = svdup_n_f32(0.f);
    svfloat32_t ap = svdup_n_f32(0.f);
    int j = 0;
    // main j loop
    for(; j<nj; j+=NJ_EPEP){
#ifdef LOOP_DIVISION
      float dx_tmp[16*NJ_EPEP];
      float dy_tmp[16*NJ_EPEP];
      float dz_tmp[16*NJ_EPEP];
      float r2_tmp[16*NJ_EPEP];
#ifdef DETAILED_PROFILE
      double tmp = get_dtime();
#endif
#endif
      for(int jj=0;jj<NJ_EPEP;jj++){
	const int jjj = j + jj;
	svfloat32_t dx = svdup_n_f32(epj[jjj].x);
	svfloat32_t dy = svdup_n_f32(epj[jjj].y);
	svfloat32_t dz = svdup_n_f32(epj[jjj].z);
	dx = svsubr_f32_x(pg,xi.v0,dx);
	dy = svsubr_f32_x(pg,xi.v1,dy);
	dz = svsubr_f32_x(pg,xi.v2,dz);
	svfloat32_t r2 = svmad_n_f32_x(pg,dx,dx,eps2);
	r2 = svmad_f32_x(pg,dy,dy,r2);
	r2 = svmad_f32_x(pg,dz,dz,r2);// 3 sub, 3 mad
#ifdef LOOP_DIVISION
	svst1_f32(pg,dx_tmp+16*jj,dx);
	svst1_f32(pg,dy_tmp+16*jj,dy);
	svst1_f32(pg,dz_tmp+16*jj,dz);
	svst1_f32(pg,r2_tmp+16*jj,r2);
      }
#ifdef DETAILED_PROFILE
      t0 += get_dtime() - tmp;
      tmp = get_dtime();
#endif
      for(int jj=0;jj<NJ_EPEP;jj++){
	const int jjj = j + jj;
	svfloat32_t r2 = svld1_f32(pg,r2_tmp+16*jj);
	svfloat32_t dx = svld1_f32(pg,dx_tmp+16*jj);
	svfloat32_t dy = svld1_f32(pg,dy_tmp+16*jj);
	svfloat32_t dz = svld1_f32(pg,dz_tmp+16*jj);
#endif
	r2 = svmax_n_f32_x(pg,r2,r_out2);
	svfloat32_t rinv  = svrsqrte_f32(r2);
	svfloat32_t h = svmul_f32_x(pg,r2,rinv);
	h = svmsb_n_f32_x(pg,h,rinv,1.f);
	svfloat32_t poly = svmad_n_f32_x(pg,h,svdup_f32(0.375f),0.5f);
	poly = svmul_f32_x(pg,poly,h);
	rinv = svmad_f32_x(pg,rinv,poly,rinv); // 2 mul, 3 fma, 1 sqrt, 1 max
	svfloat32_t mri3 = svdup_n_f32(epj[jjj].m);
	mri3 = svmul_f32_x(pg,mri3,rinv);
	//ap = svsub_f32_x(pg,ap,mri3);
	mri3 = svmul_f32_x(pg,mri3,rinv);
	mri3 = svmul_f32_x(pg,mri3,rinv);
	ax = svmad_f32_x(pg,mri3,dx,ax);
	ay = svmad_f32_x(pg,mri3,dy,ay);
	az = svmad_f32_x(pg,mri3,dz,az); // 3 mul, 3 fma
	// total 8 nonfma, 9 fma, 2 math func
      }
#ifdef LOOP_DIVISION
#ifdef DETAILED_PROFILE
      t1 += get_dtime() - tmp;
      tmp = get_dtime();
#endif
#endif
    }
    // tail
    for(; j<nj; j++){
      svfloat32_t dx = svsubr_n_f32_x(pg,xi.v0,epj[j].x);
      svfloat32_t dy = svsubr_n_f32_x(pg,xi.v1,epj[j].y);
      svfloat32_t dz = svsubr_n_f32_x(pg,xi.v2,epj[j].z);
      svfloat32_t r2 = svmad_n_f32_x(pg,dx,dx,eps2);
      r2 = svmad_f32_x(pg,dy,dy,r2);
      r2 = svmad_f32_x(pg,dz,dz,r2);
      r2 = svmax_n_f32_x(pg,r2,r_out2);
      svfloat32_t rinv  = svrsqrte_f32(r2);
      svfloat32_t h = svmul_f32_x(pg,r2,rinv);
      h = svmsb_n_f32_x(pg,h,rinv,1.f);
      svfloat32_t poly = svmad_n_f32_x(pg,h,svdup_f32(0.375f),0.5f);
      poly = svmul_f32_x(pg,poly,h);
      rinv = svmad_f32_x(pg,rinv,poly,rinv); // 2 mul, 3 mad
      svfloat32_t mri3 = svdup_n_f32(epj[j].m);
      mri3 = svmul_f32_x(pg,mri3,rinv);
      ap = svsub_f32_x(pg,ap,mri3);
      mri3 = svmul_f32_x(pg,mri3,rinv);
      mri3 = svmul_f32_x(pg,mri3,rinv);
      ax = svmad_f32_x(pg,mri3,dx,ax);
      ay = svmad_f32_x(pg,mri3,dy,ay);
      az = svmad_f32_x(pg,mri3,dz,az);
    }
#if 0
    svfloat32x4_t ai = {ax,ay,az,ap};
    svst4_f32(pg,(float*)&force[i],ai);
#else
    //gather_st4_f32(&ax,&ay,&az,&ap);
    float tmpx[16],tmpy[16],tmpz[16],tmpp[16];
    svst1_f32(pg,tmpx,ax);
    svst1_f32(pg,tmpy,ay);
    svst1_f32(pg,tmpz,az);
    svst1_f32(pg,tmpp,ap);
    for(int j=0;j<16;j++){
      force[i+j].x += tmpx[j];
      force[i+j].y += tmpy[j];
      force[i+j].z += tmpz[j];
      force[i+j].p += tmpp[j];
    }
#endif
  }
#ifdef PROFILE
  fipp_stop();
#endif

#ifdef DETAILED_PROFILE
  printf("performance of loop0 is %lf GFlops (%lf sec)\n",9.5*ni*nj/(t0*1e9),t0);
  printf("performance of loop1 is %lf GFlops (%lf sec)\n",17.5*ni*nj/(t1*1e9),t1);
#endif
}

//#define NO_LOOP_DIVISION
//__attribute__ ((noinline))
void CalcForceLongEpSp(const EPI* epi,const int ni,
		       const SPJ* spj,const int nj,
		       Force* force){
#ifdef PROFILE
  fipp_start();
#endif

  svbool_t pg = svptrue_b32();
  const float r_out2 = r_out*r_out;
  for(int i=0; i<ni; i+=16){
    const svfloat32x3_t xi = svld3_f32(pg,(float*)&epi[i]);

    //svfloat32x4_t ai = {svdup_n_f32(0.f),svdup_n_f32(0.f),svdup_n_f32(0.f),svdup_n_f32(0.f)};
    svfloat32_t ax = svdup_n_f32(0.f);
    svfloat32_t ay = svdup_n_f32(0.f);
    svfloat32_t az = svdup_n_f32(0.f);
    svfloat32_t ap = svdup_n_f32(0.f);
#ifdef LOOP_DIVISION
    float dx_tmp[16*NJ_EPSP],dy_tmp[16*NJ_EPSP],dz_tmp[16*NJ_EPSP],r2_tmp[16*NJ_EPSP];
#endif
    int j = 0;
    for(; j<nj; j+=NJ_EPSP){
      for(int jj=0;jj<NJ_EPSP;jj++){
	const int jjj = j + jj;
	svfloat32_t dx = svsub_n_f32_x(pg,xi.v0,spj[jjj].x);
	svfloat32_t dy = svsub_n_f32_x(pg,xi.v1,spj[jjj].y);
	svfloat32_t dz = svsub_n_f32_x(pg,xi.v2,spj[jjj].z);
	svfloat32_t r2  = svmad_n_f32_x(pg,dx,dx,eps2);
	r2 = svmad_f32_x(pg,dy,dy,r2);
	r2 = svmad_f32_x(pg,dz,dz,r2);
	r2 = svmax_n_f32_x(pg,r2,r_out2); // 3 sub, 4 fma, 1 max
#ifdef LOOP_DIVISION
	svst1_f32(pg,dx_tmp+16*jj,dx);
	svst1_f32(pg,dy_tmp+16*jj,dy);
	svst1_f32(pg,dz_tmp+16*jj,dz);
	svst1_f32(pg,r2_tmp+16*jj,r2);
      }
      float ri2_tmp[16*NJ_EPSP],ri3_tmp[16*NJ_EPSP];
      for(int jj=0;jj<NJ_EPSP;jj++){
	svfloat32_t r2 =  svld1_f32(pg,r2_tmp+16*jj);
#endif
	svfloat32_t rinv  = svrsqrte_f32(r2);
	svfloat32_t h = svmul_f32_x(pg,r2,rinv);
	h = svmsb_n_f32_x(pg,h,rinv,1.f);
	svfloat32_t poly = svmad_n_f32_x(pg,h,svdup_f32(0.375f),0.5f);
	poly = svmul_f32_x(pg,poly,h);
	rinv = svmad_f32_x(pg,rinv,poly,rinv);
	svfloat32_t ri2 = svmul_f32_x(pg,rinv,rinv);
	svfloat32_t ri3 = svmul_f32_x(pg,rinv,ri2);  // 4 mul, 3 fma, 1 sqrt
#ifdef LOOP_DIVISION
	svst1_f32(pg,ri2_tmp+16*jj,ri2);
	svst1_f32(pg,ri3_tmp+16*jj,ri3);
      }
      float fx_tmp[16*NJ_EPSP],fy_tmp[16*NJ_EPSP],fz_tmp[16*NJ_EPSP];
      for(int jj=0;jj<NJ_EPSP;jj++){
	const int jjj = j+jj;
	svfloat32_t dx  = svld1_f32(pg, dx_tmp+16*jj);
	svfloat32_t dy  = svld1_f32(pg, dy_tmp+16*jj);
	svfloat32_t dz  = svld1_f32(pg, dz_tmp+16*jj);
#endif
	svfloat32_t fx;
	fx = svmul_n_f32_x(pg,dx,spj[jjj].xx);
	fx = svmla_n_f32_x(pg,fx,dy,spj[jjj].xy);
	fx = svmla_n_f32_x(pg,fx,dz,spj[jjj].zx);
	svfloat32_t fy;
	fy = svmul_n_f32_x(pg,dy,spj[jjj].yy);
	fy = svmla_n_f32_x(pg,fy,dx,spj[jjj].xy);
	fy = svmla_n_f32_x(pg,fy,dz,spj[jjj].yz);
	svfloat32_t fz;
	fz = svmul_n_f32_x(pg,dz,spj[jjj].zz);
	fz = svmla_n_f32_x(pg,fz,dx,spj[jjj].zx);
	fz = svmla_n_f32_x(pg,fz,dy,spj[jjj].yz); // 3 mul, 6 fma
	//az = svmsb_f32_x(pg,ri5,fz,az);
#ifdef LOOP_DIVISION
	svst1_f32(pg,fx_tmp+16*jj,fx);
	svst1_f32(pg,fy_tmp+16*jj,fy);
	svst1_f32(pg,fz_tmp+16*jj,fz);
      }
      for(int jj=0;jj<NJ_EPSP;jj++){
	const int jjj=j+jj;
	svfloat32_t dx = svld1_f32(pg,dx_tmp+16*jj);
	svfloat32_t dy = svld1_f32(pg,dy_tmp+16*jj);
	svfloat32_t dz = svld1_f32(pg,dz_tmp+16*jj);
	svfloat32_t ri2 = svld1_f32(pg,ri2_tmp+16*jj);
	svfloat32_t ri3 = svld1_f32(pg,ri3_tmp+16*jj);
	svfloat32_t qrx = svld1_f32(pg,fx_tmp+16*jj);
	svfloat32_t qry = svld1_f32(pg,fy_tmp+16*jj);
	svfloat32_t qrz = svld1_f32(pg,fz_tmp+16*jj);
#endif
	svfloat32_t meff3 = svmad_n_f32_x(pg,qrx,dx,spj[jjj].tr);
	meff3 = svmad_f32_x(pg,qry,dy,meff3);
	meff3 = svmad_f32_x(pg,qrz,dz,meff3);
	meff3 = svmul_f32_x(pg,meff3,ri2);
	meff3 = svmul_n_f32_x(pg,meff3,2.5f); // 2 mul, 3 fma

	svfloat32_t mj = svdup_n_f32(spj[j+jj].m);
	svfloat32_t fx = svmsb_f32_x(pg,meff3,dx,qrx);
	svfloat32_t mjx = svmul_f32_x(pg,mj,dx);
	fx = svmsb_f32_x(pg,fx,ri2,mjx); // 1 mul, 2 fma

	svfloat32_t fy = svmsb_f32_x(pg,meff3,dy,qry);
	svfloat32_t mjy = svmul_f32_x(pg,mj,dy);
	fy = svmsb_f32_x(pg,fy,ri2,mjy); // 1 mul, 2 fma

	svfloat32_t fz = svmsb_f32_x(pg,meff3,dz,qrz);
	svfloat32_t mjz = svmul_f32_x(pg,mj,dz);
	fz = svmsb_f32_x(pg,fz,ri2,mjz); // 1 mul, 2 fma
	//svfloat32_t meff  = svmla_n_f32_x(pg,mj,rqr_ri4,rqr_ri4,0.5f);
	ax = svmad_f32_x(pg,ri3,fx,ax);
	ay = svmad_f32_x(pg,ri3,fy,ay);
	az = svmad_f32_x(pg,ri3,fz,az); // 3 fma
	//ai.v3 = svmsb_f32_x(pg,meff,rinv,pot);
	// total 15 nonfma, 25 fma, 2 math func(= 1 nonfma)
      }
    }
    #if 1
    for(;j<nj;j++){
      svfloat32_t xj = svdup_n_f32(spj[j].x);
      svfloat32_t yj = svdup_n_f32(spj[j].y);
      svfloat32_t zj = svdup_n_f32(spj[j].z);

      svfloat32_t dx = svsub_f32_x(pg,xi.v0,xj);
      svfloat32_t dy = svsub_f32_x(pg,xi.v1,yj);
      svfloat32_t dz = svsub_f32_x(pg,xi.v2,zj);
      svfloat32_t r2  = svmad_n_f32_x(pg,dx,dx,eps2);
      r2 = svmad_f32_x(pg,dy,dy,r2);
      r2 = svmad_f32_x(pg,dz,dz,r2);
      r2 = svmax_n_f32_x(pg,r2,r_out2); // 3 sub, 4 fma, 1 max
      svfloat32_t rinv  = svrsqrte_f32(r2);
      svfloat32_t h = svmul_f32_x(pg,r2,rinv);
      h = svmsb_n_f32_x(pg,h,rinv,1.f);
      svfloat32_t poly = svmad_n_f32_x(pg,h,svdup_f32(0.375f),0.5f);
      poly = svmul_f32_x(pg,poly,h);
      rinv = svmad_f32_x(pg,rinv,poly,rinv);
      svfloat32_t ri2 = svmul_f32_x(pg,rinv,rinv);
      svfloat32_t ri3 = svmul_f32_x(pg,rinv,ri2);  // 4 mul, 3 fma, 1 sqrt
      svfloat32_t fx;
      fx = svmul_n_f32_x(pg,dx,spj[j].xx);
      fx = svmla_n_f32_x(pg,fx,dy,spj[j].xy);
      fx = svmla_n_f32_x(pg,fx,dz,spj[j].zx);
      svfloat32_t fy;
      fy = svmul_n_f32_x(pg,dy,spj[j].yy);
      fy = svmla_n_f32_x(pg,fy,dx,spj[j].xy);
      fy = svmla_n_f32_x(pg,fy,dz,spj[j].yz);
      svfloat32_t fz;
      fz = svmul_n_f32_x(pg,dz,spj[j].zz);
      fz = svmla_n_f32_x(pg,fz,dx,spj[j].zx);
      fz = svmla_n_f32_x(pg,fz,dy,spj[j].yz); // 3 mul, 6 fma
      svfloat32_t meff3 = svmad_n_f32_x(pg,fx,dx,spj[j].tr);
      meff3 = svmad_f32_x(pg,fy,dy,meff3);
      meff3 = svmad_f32_x(pg,fz,dz,meff3);
      meff3 = svmul_f32_x(pg,meff3,ri2);
      meff3 = svmul_n_f32_x(pg,meff3,2.5f); // 2 mul, 3 fma

      svfloat32_t mj = svdup_n_f32(spj[j].m);
      fx = svmsb_f32_x(pg,meff3,dx,fx);
      svfloat32_t mjx = svmul_f32_x(pg,mj,dx);
      fx = svmsb_f32_x(pg,fx,ri2,mjx); // 1 mul, 2 fma

      fy = svmsb_f32_x(pg,meff3,dy,fy);
      svfloat32_t mjy = svmul_f32_x(pg,mj,dy);
      fy = svmsb_f32_x(pg,fy,ri2,mjy); // 1 mul, 2 fma

      fz = svmsb_f32_x(pg,meff3,dz,fz);
      svfloat32_t mjz = svmul_f32_x(pg,mj,dz);
      fz = svmsb_f32_x(pg,fz,ri2,mjz); // 1 mul, 2 fma
      //svfloat32_t meff  = svmla_n_f32_x(pg,mj,rqr_ri4,rqr_ri4,0.5f);
      ax = svmad_f32_x(pg,ri3,fx,ax);
      ay = svmad_f32_x(pg,ri3,fy,ay);
      az = svmad_f32_x(pg,ri3,fz,az); // 3 fma
      //ai.v3 = svmsb_f32_x(pg,meff,rinv,pot);
      // total 15 nonfma, 25 fma, 2 math func(= 1 nonfma)
    }
    #endif
    float tmpx[16],tmpy[16],tmpz[16],tmpp[16];
    svst1_f32(pg,tmpx,ax);
    svst1_f32(pg,tmpy,ay);
    svst1_f32(pg,tmpz,az);
    svst1_f32(pg,tmpp,ap);
    for(int jj=0;jj<16;jj++){
      //printf("%d: %e %e %e\n",i,tmpx[jj],tmpy[jj],tmpz[jj]);
      force[i+jj].x += tmpx[jj];
      force[i+jj].y += tmpy[jj];
      force[i+jj].z += tmpz[jj];
      force[i+jj].p += tmpp[jj];
    }
    //svst4_f32(pg,(float*)&force[i],ai);
  }
#ifdef PROFILE
  fipp_stop();
#endif
}


