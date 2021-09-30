#include "user_defined_class.h"
#include "constant.h"
#include <arm_sve.h>

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

#define DEBUG_PRINT(x) float tmp[16];svst1_f32(pg,tmp,x);for(int ii=0;ii<16;ii++) printf("%d %d: %e\n",i+ii,j,tmp[ii]);

#define NJ_MAX 1024
void CalcForceLongEpEp(const EPI* epi,const int ni,
		       const EPJ* epj,const int nj,
		       Force* force){
  svbool_t pg = svptrue_b32();
  const float r_out2 = r_out*r_out;
  printf("--- original epep kernel ---\n");
  for(int i=0; i<ni; i+=16){
    const svfloat32x3_t xi = svld3_f32(pg,(float*)&epi[i]);
    //svfloat32x4_t ai = {svdup_n_f32(0.f),svdup_n_f32(1.f),svdup_n_f32(2.f),svdup_n_f32(3.f)};
    svfloat32_t ax = svdup_n_f32(0.f);
    svfloat32_t ay = svdup_n_f32(0.f);
    svfloat32_t az = svdup_n_f32(0.f);
    svfloat32_t ap = svdup_n_f32(0.f);
    for(int j=0; j<nj; j++){
      svfloat32_t dx = svdup_n_f32(epj[j].x);
      svfloat32_t dy = svdup_n_f32(epj[j].y);
      svfloat32_t dz = svdup_n_f32(epj[j].z);
      dx = svsub_f32_x(pg,dx,xi.v0);
      dy = svsub_f32_x(pg,dy,xi.v1);
      dz = svsub_f32_x(pg,dz,xi.v2);
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
      force[i+j].x = tmpx[j];
      force[i+j].y = tmpy[j];
      force[i+j].z = tmpz[j];
      force[i+j].p = tmpp[j];
    }
#endif
  }
}

//#define NO_LOOP_DIVISION
//__attribute__ ((noinline))
void CalcForceLongEpSp(const EPI* epi,const int ni,
		       const SPJ* spj,const int nj,
		       Force* force){
  svbool_t pg = svptrue_b32();
  const float r_out2 = r_out*r_out;
  for(int i=0; i<ni; i+=8){
    const svfloat32x3_t xi = svld3_f32(pg,(float*)&epi[i]);

    //svfloat32x4_t ai = {svdup_n_f32(0.f),svdup_n_f32(0.f),svdup_n_f32(0.f),svdup_n_f32(0.f)};
    svfloat32_t ax = svdup_n_f32(0.f);
    svfloat32_t ay = svdup_n_f32(0.f);
    svfloat32_t az = svdup_n_f32(0.f);
    svfloat32_t ap = svdup_n_f32(0.f);
    for(int j=0; j<nj; j+=2){
      svfloat32_t xj = svdup_n_f32(spj[j].x);
      svfloat32_t yj = svdup_n_f32(spj[j].y);
      svfloat32_t zj = svdup_n_f32(spj[j].z);

      svfloat32_t dx = svsub_f32_x(pg,xj,xi.v0);
      svfloat32_t dy = svsub_f32_x(pg,yj,xi.v1);
      svfloat32_t dz = svsub_f32_x(pg,zj,xi.v2);

      svfloat32_t r2  = svmad_n_f32_x(pg,dx,dx,eps2);
      r2 = svmad_f32_x(pg,dy,dy,r2);
      r2 = svmad_f32_x(pg,dz,dz,r2);
      r2 = svmax_n_f32_x(pg,r2,r_out2);
      svfloat32_t rinv  = svrsqrte_f32(r2);
      svfloat32_t h = svmul_f32_x(pg,r2,rinv);
      h = svnmsb_n_f32_x(pg,h,rinv,1.f);
      svfloat32_t poly = svmad_n_f32_x(pg,h,svdup_f32(3.f/8.f),0.5f);
      poly = svmul_f32_x(pg,poly,h);
      rinv = svmad_f32_x(pg,rinv,poly,rinv); // 2 mul, 3 mad
      svfloat32_t qxx = svdup_n_f32(spj[j].xx);
      svfloat32_t qyy = svdup_n_f32(spj[j].yy);
      svfloat32_t qzz = svdup_n_f32(spj[j].zz);

      svfloat32_t qxy = svdup_n_f32(spj[j].xy);
      svfloat32_t qyz = svdup_n_f32(spj[j].yz);
      svfloat32_t qzx = svdup_n_f32(spj[j].zx);

      svfloat32_t qr_x,qr_y,qr_z;
      qr_x = svmul_f32_x(pg,qxx,dx);
      qr_x = svmad_f32_x(pg,qxy,dy,qr_x);
      qr_x = svmad_f32_x(pg,qzx,dz,qr_x);
      qr_y = svmul_f32_x(pg,qyy,dy);
      qr_y = svmad_f32_x(pg,qxy,dx,qr_y);
      qr_y = svmad_f32_x(pg,qyz,dz,qr_y);
      qr_z = svmul_f32_x(pg,qzz,dz);
      qr_z = svmad_f32_x(pg,qzx,dx,qr_z);
      qr_z = svmad_f32_x(pg,qyz,dy,qr_z);
      svfloat32_t ri2 = svmul_f32_x(pg,rinv,rinv);
      svfloat32_t ri3 = svmul_f32_x(pg,rinv,ri2);
      svfloat32_t ri4 = svmul_f32_x(pg,ri2,ri2);
      svfloat32_t ri5 = svmul_f32_x(pg,ri2,ri3);
      ax = svmsb_f32_x(pg,ri5,qr_x,ax);
      ay = svmsb_f32_x(pg,ri5,qr_y,ay);
      az = svmsb_f32_x(pg,ri5,qr_z,az);

      svfloat32_t mtr = svdup_n_f32(spj[j].tr);
      svfloat32_t rqr,rqr_ri4;
      rqr = svmad_f32_x(pg,qr_x,dx,mtr);
      rqr = svmad_f32_x(pg,qr_y,dy,rqr);
      rqr = svmad_f32_x(pg,qr_z,dz,rqr);
      rqr_ri4 = svmul_f32_x(pg,rqr,ri4);
      svfloat32_t mj = svdup_n_f32(spj[j].m);
      //svfloat32_t meff  = svmla_n_f32_x(pg,mj,rqr_ri4,rqr_ri4,0.5f);
      svfloat32_t meff3 = svmla_n_f32_x(pg,mj,rqr_ri4,2.5f);
      meff3 = svmul_f32_x(pg,meff3,ri3);
      ax = svmad_f32_x(pg,meff3,dx,ax);
      ay = svmad_f32_x(pg,meff3,dy,ay);
      az = svmad_f32_x(pg,meff3,dz,az);
      //ai.v3 = svmsb_f32_x(pg,meff,rinv,pot);
    }
    float tmpx[16],tmpy[16],tmpz[16],tmpp[16];
    svst1_f32(pg,tmpx,ax);
    svst1_f32(pg,tmpy,ay);
    svst1_f32(pg,tmpz,az);
    svst1_f32(pg,tmpp,ap);
    for(int j=0;j<16;j++){
      force[i+j].x = tmpx[j];
      force[i+j].y = tmpy[j];
      force[i+j].z = tmpz[j];
      force[i+j].p = tmpp[j];
    }
    //svst4_f32(pg,(float*)&force[i],ai);
  }
}
