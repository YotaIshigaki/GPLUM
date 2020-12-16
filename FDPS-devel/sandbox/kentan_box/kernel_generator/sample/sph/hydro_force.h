#include <arm_sve.h>
template<typename Tepi,typename Tepj,typename Tforce>
struct Kernel{
PS::F32 gamma;
PS::F32 cfl;
Kernel(PS::F32 gamma,PS::F32 cfl):gamma(gamma),cfl(cfl){}
void operator()(const Tepi* epi,const int ni,const Tepj* epj,const int nj,Tforce *force){
svbool_t pg0 = svptrue_b32();
svfloat32x3_t rij;
svfloat32_t r2;
svfloat32_t r_inv;
svfloat32_t r;
svfloat32x3_t vij;
svfloat32_t dvdr;
svfloat32_t wij;
svfloat32_t v_sig;
svfloat32_t aij;
svfloat32_t AV;
svfloat32_t __fkg_tmp0;
svfloat32_t __fkg_tmp1;
svfloat32_t __fkg_tmp2;
svfloat32_t __fkg_tmp3;
svfloat32_t __fkg_tmp4;
svfloat32_t __fkg_tmp5;
svfloat32_t __fkg_tmp6;
svfloat32_t __fkg_tmp7;
svfloat32_t __fkg_tmp8;
svfloat32_t __fkg_tmp20;
svfloat32x3_t gradWi;
svfloat32_t __fkg_tmp10;
svfloat32_t __fkg_tmp11;
svfloat32_t __fkg_tmp12;
svfloat32_t __fkg_tmp13;
svfloat32_t __fkg_tmp14;
svfloat32_t __fkg_tmp15;
svfloat32_t __fkg_tmp16;
svfloat32_t __fkg_tmp17;
svfloat32_t __fkg_tmp18;
svfloat32_t __fkg_tmp21;
svfloat32x3_t gradWj;
svfloat32x3_t gradWij;
svfloat32_t g1u;
svfloat32_t fij;
svfloat32_t fji;
svfloat32_t Pi_inv;
svfloat32_t Pj_inv;
svfloat32_t g1u2;
svfloat32_t tmp0;
svfloat32_t tmp1;
svfloat32_t __fkg_tmp22;
PS::F32vec xiloc_tmp[ni];
PS::F32vec xjloc_tmp[nj];
for(int i=0;i<ni;i++){
xiloc_tmp[i].x = (epi[i].pos.x-epi[0].pos.x);
}
for(int i=0;i<ni;i++){
xiloc_tmp[i].y = (epi[i].pos.y-epi[0].pos.y);
}
for(int i=0;i<ni;i++){
xiloc_tmp[i].z = (epi[i].pos.z-epi[0].pos.z);
}
for(int i=0;i<nj;i++){
xjloc_tmp[i].x = (epj[i].pos.x-epi[0].pos.x);
}
for(int i=0;i<nj;i++){
xjloc_tmp[i].y = (epj[i].pos.y-epi[0].pos.y);
}
for(int i=0;i<nj;i++){
xjloc_tmp[i].z = (epj[i].pos.z-epi[0].pos.z);
}
for(int i=0;i<((ni+15)/16)*16;i+=16){
pg0 = svwhilelt_b32_s32(i,ni);
svfloat32x3_t vi;
uint32_t index_viv0[16] = {6,26,46,66,86,106,126,146,166,186,206,226,246,266,286,306};
svuint32_t vindex_viv0 = svld1_u32(svptrue_b32(),index_viv0);
vi.v0 = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_viv0);
uint32_t index_viv1[16] = {7,27,47,67,87,107,127,147,167,187,207,227,247,267,287,307};
svuint32_t vindex_viv1 = svld1_u32(svptrue_b32(),index_viv1);
vi.v1 = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_viv1);
uint32_t index_viv2[16] = {8,28,48,68,88,108,128,148,168,188,208,228,248,268,288,308};
svuint32_t vindex_viv2 = svld1_u32(svptrue_b32(),index_viv2);
vi.v2 = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_viv2);
svfloat32_t mi;
uint32_t index_mi[16] = {9,29,49,69,89,109,129,149,169,189,209,229,249,269,289,309};
svuint32_t vindex_mi = svld1_u32(svptrue_b32(),index_mi);
mi = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_mi);
svfloat32_t ui;
uint32_t index_ui[16] = {10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310};
svuint32_t vindex_ui = svld1_u32(svptrue_b32(),index_ui);
ui = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_ui);
svfloat32_t hi;
uint32_t index_hi[16] = {11,31,51,71,91,111,131,151,171,191,211,231,251,271,291,311};
svuint32_t vindex_hi = svld1_u32(svptrue_b32(),index_hi);
hi = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_hi);
svfloat32_t rhoi;
uint32_t index_rhoi[16] = {12,32,52,72,92,112,132,152,172,192,212,232,252,272,292,312};
svuint32_t vindex_rhoi = svld1_u32(svptrue_b32(),index_rhoi);
rhoi = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_rhoi);
svfloat32_t Pi;
uint32_t index_Pi[16] = {13,33,53,73,93,113,133,153,173,193,213,233,253,273,293,313};
svuint32_t vindex_Pi = svld1_u32(svptrue_b32(),index_Pi);
Pi = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_Pi);
svfloat32_t fi;
uint32_t index_fi[16] = {14,34,54,74,94,114,134,154,174,194,214,234,254,274,294,314};
svuint32_t vindex_fi = svld1_u32(svptrue_b32(),index_fi);
fi = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_fi);
svfloat32_t ci;
uint32_t index_ci[16] = {15,35,55,75,95,115,135,155,175,195,215,235,255,275,295,315};
svuint32_t vindex_ci = svld1_u32(svptrue_b32(),index_ci);
ci = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_ci);
svfloat32_t BalSWi;
uint32_t index_BalSWi[16] = {16,36,56,76,96,116,136,156,176,196,216,236,256,276,296,316};
svuint32_t vindex_BalSWi = svld1_u32(svptrue_b32(),index_BalSWi);
BalSWi = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_BalSWi);
svfloat32_t ai;
uint32_t index_ai[16] = {17,37,57,77,97,117,137,157,177,197,217,237,257,277,297,317};
svuint32_t vindex_ai = svld1_u32(svptrue_b32(),index_ai);
ai = svld1_gather_u32index_f32(pg0,(float*)&epi[i],vindex_ai);
svfloat32x3_t xiloc;
xiloc = svld3_f32(pg0,(float*)(xiloc_tmp+i));
svfloat32x3_t acc;
uint32_t index_accv0[16] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
svuint32_t vindex_accv0 = svld1_u32(svptrue_b32(),index_accv0);
acc.v0 = svld1_gather_u32index_f32(pg0,(float*)&force[i],vindex_accv0);
uint32_t index_accv1[16] = {1,6,11,16,21,26,31,36,41,46,51,56,61,66,71,76};
svuint32_t vindex_accv1 = svld1_u32(svptrue_b32(),index_accv1);
acc.v1 = svld1_gather_u32index_f32(pg0,(float*)&force[i],vindex_accv1);
uint32_t index_accv2[16] = {2,7,12,17,22,27,32,37,42,47,52,57,62,67,72,77};
svuint32_t vindex_accv2 = svld1_u32(svptrue_b32(),index_accv2);
acc.v2 = svld1_gather_u32index_f32(pg0,(float*)&force[i],vindex_accv2);
svfloat32_t eng;
uint32_t index_eng[16] = {3,8,13,18,23,28,33,38,43,48,53,58,63,68,73,78};
svuint32_t vindex_eng = svld1_u32(svptrue_b32(),index_eng);
eng = svld1_gather_u32index_f32(pg0,(float*)&force[i],vindex_eng);
svfloat32_t dt;
uint32_t index_dt[16] = {4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79};
svuint32_t vindex_dt = svld1_u32(svptrue_b32(),index_dt);
dt = svld1_gather_u32index_f32(pg0,(float*)&force[i],vindex_dt);
#pragma loop loop_fission
#pragma loop loop_fission_stripmining L1
for(int j=0;j<nj;j++){
svfloat32x3_t vj;
vj.v0 = svdup_f32(epj[j].vel.x);
vj.v1 = svdup_f32(epj[j].vel.y);
vj.v2 = svdup_f32(epj[j].vel.z);
svfloat32_t mj;
mj = svdup_f32(epj[j].mass);
svfloat32_t uj;
uj = svdup_f32(epj[j].eng);
svfloat32_t hj;
hj = svdup_f32(epj[j].h);
svfloat32_t rhoj;
rhoj = svdup_f32(epj[j].dens);
svfloat32_t Pj;
Pj = svdup_f32(epj[j].pres);
svfloat32_t fj;
fj = svdup_f32(epj[j].gradh);
svfloat32_t cj;
cj = svdup_f32(epj[j].snds);
svfloat32_t BalSWj;
BalSWj = svdup_f32(epj[j].BalSW);
svfloat32_t aj;
aj = svdup_f32(epj[j].alpha);
svfloat32x3_t xjloc;
xjloc.v0 = svdup_f32(xjloc_tmp[j].x);
xjloc.v1 = svdup_f32(xjloc_tmp[j].y);
xjloc.v2 = svdup_f32(xjloc_tmp[j].z);
rij.v0 = svsub_f32_z(pg0,xiloc.v0,xjloc.v0);
rij.v1 = svsub_f32_z(pg0,xiloc.v1,xjloc.v1);
rij.v2 = svsub_f32_z(pg0,xiloc.v2,xjloc.v2);
r2 = madd(pg0,rij.v2,rij.v2,madd(pg0,rij.v0,rij.v0,svmul_f32_z(pg0,rij.v1,rij.v1)));
svbool_t pg1;
pg1=svcmpgt_f32(pg0,r2,svdup_n_f32(0.0f));
r_inv = rsqrt(pg1,r2);

r = svmul_f32_z(pg0,r2,r_inv);
#pragma statement fission_point
vij.v0 = svsub_f32_z(pg0,vi.v0,vj.v0);
vij.v1 = svsub_f32_z(pg0,vi.v1,vj.v1);
vij.v2 = svsub_f32_z(pg0,vi.v2,vj.v2);
dvdr = madd(pg0,rij.v2,vij.v2,madd(pg0,rij.v0,vij.v0,svmul_f32_z(pg0,rij.v1,vij.v1)));
pg1=svcmpge_f32(pg0,dvdr,svdup_n_f32(0.0f));
wij = svdup_n_f32(0.0f);
pg1 = svnot_b_z(pg0,svcmpge_f32(pg1,dvdr,svdup_n_f32(0.0f)));
wij = svmul_f32_z(pg1,dvdr,r_inv);

#pragma statement fission_point
v_sig = nmadd(pg0,svdup_n_f32(3.0f),wij,svneg_f32_z(pg0,svadd_f32_z(pg0,ci,cj)));
dt = max(pg0,dt,v_sig);
aij = svmul_f32_z(pg0,svdup_n_f32(0.5f),svadd_f32_z(pg0,ai,aj));
AV = svmul_f32_z(pg0,svmul_f32_z(pg0,svmul_f32_z(pg0,svmul_f32_z(pg0,svmul_f32_z(pg0,svneg_f32_z(pg0,svdup_n_f32(0.5f)),aij),v_sig),wij),inv(pg0,svadd_f32_z(pg0,rhoi,rhoj))),svadd_f32_z(pg0,BalSWi,BalSWj));
#pragma statement fission_point
__fkg_tmp0 = r;
__fkg_tmp1 = hi;
__fkg_tmp2 = svmul_f32_z(pg0,__fkg_tmp0,inv(pg0,__fkg_tmp1));
__fkg_tmp3 = max(pg0,svdup_n_f32(0.0f),svsub_f32_z(pg0,svdup_n_f32(1.0f),__fkg_tmp2));
__fkg_tmp4 = svmul_f32_z(pg0,__fkg_tmp1,__fkg_tmp1);
__fkg_tmp5 = svmul_f32_z(pg0,svmul_f32_z(pg0,__fkg_tmp4,__fkg_tmp4),__fkg_tmp1);
__fkg_tmp6 = svmul_f32_z(pg0,svdup_n_f32(1155.0f),inv(pg0,svmul_f32_z(pg0,svdup_n_f32(12.5663706144f),__fkg_tmp5)));
__fkg_tmp7 = svmul_f32_z(pg0,__fkg_tmp3,__fkg_tmp3);
__fkg_tmp8 = svmul_f32_z(pg0,svmul_f32_z(pg0,__fkg_tmp7,__fkg_tmp7),__fkg_tmp3);
__fkg_tmp20 = svmul_f32_z(pg0,svmul_f32_z(pg0,svneg_f32_z(pg0,__fkg_tmp6),__fkg_tmp8),madd(pg0,svdup_n_f32(5.0f),__fkg_tmp2,svdup_n_f32(1.0f)));
gradWi.v0 = svmul_f32_z(pg0,__fkg_tmp20,rij.v0);
gradWi.v1 = svmul_f32_z(pg0,__fkg_tmp20,rij.v1);
gradWi.v2 = svmul_f32_z(pg0,__fkg_tmp20,rij.v2);
#pragma statement fission_point
__fkg_tmp10 = r;
__fkg_tmp11 = hj;
__fkg_tmp12 = svmul_f32_z(pg0,__fkg_tmp10,inv(pg0,__fkg_tmp11));
__fkg_tmp13 = max(pg0,svdup_n_f32(0.0f),svsub_f32_z(pg0,svdup_n_f32(1.0f),__fkg_tmp12));
__fkg_tmp14 = svmul_f32_z(pg0,__fkg_tmp11,__fkg_tmp11);
__fkg_tmp15 = svmul_f32_z(pg0,svmul_f32_z(pg0,__fkg_tmp14,__fkg_tmp14),__fkg_tmp11);
__fkg_tmp16 = svmul_f32_z(pg0,svdup_n_f32(1155.0f),inv(pg0,svmul_f32_z(pg0,svdup_n_f32(12.5663706144f),__fkg_tmp15)));
__fkg_tmp17 = svmul_f32_z(pg0,__fkg_tmp13,__fkg_tmp13);
__fkg_tmp18 = svmul_f32_z(pg0,svmul_f32_z(pg0,__fkg_tmp17,__fkg_tmp17),__fkg_tmp13);
__fkg_tmp21 = svmul_f32_z(pg0,svmul_f32_z(pg0,svneg_f32_z(pg0,__fkg_tmp16),__fkg_tmp18),madd(pg0,svdup_n_f32(5.0f),__fkg_tmp12,svdup_n_f32(1.0f)));
gradWj.v0 = svmul_f32_z(pg0,__fkg_tmp21,rij.v0);
gradWj.v1 = svmul_f32_z(pg0,__fkg_tmp21,rij.v1);
gradWj.v2 = svmul_f32_z(pg0,__fkg_tmp21,rij.v2);
#pragma statement fission_point
gradWij.v0 = svmul_f32_z(pg0,svdup_n_f32(0.5f),svadd_f32_z(pg0,gradWi.v0,gradWj.v0));
gradWij.v1 = svmul_f32_z(pg0,svdup_n_f32(0.5f),svadd_f32_z(pg0,gradWi.v1,gradWj.v1));
gradWij.v2 = svmul_f32_z(pg0,svdup_n_f32(0.5f),svadd_f32_z(pg0,gradWi.v2,gradWj.v2));
g1u = svsub_f32_z(pg0,svdup_n_f32(gamma),svdup_n_f32(1.0f));
fij = nmadd(pg0,fi,inv(pg0,svmul_f32_z(pg0,svmul_f32_z(pg0,g1u,mj),uj)),svneg_f32_z(pg0,svdup_n_f32(1.0f)));
fji = nmadd(pg0,fj,inv(pg0,svmul_f32_z(pg0,svmul_f32_z(pg0,g1u,mi),ui)),svneg_f32_z(pg0,svdup_n_f32(1.0f)));
#pragma statement fission_point
Pi_inv = inv(pg0,Pi);
Pj_inv = inv(pg0,Pj);
#pragma statement fission_point
g1u2 = svmul_f32_z(pg0,g1u,g1u);
tmp0 = svmul_f32_z(pg0,svmul_f32_z(pg0,svmul_f32_z(pg0,g1u2,mj),ui),uj);
tmp1 = svmul_f32_z(pg0,mj,AV);
#pragma statement fission_point
__fkg_tmp22 = svmul_f32_z(pg0,mj,AV);
acc.v0 = nmadd(pg0,__fkg_tmp22,gradWij.v0,madd(pg0,tmp0,madd(pg0,svmul_f32_z(pg0,fji,Pj_inv),gradWj.v0,svmul_f32_z(pg0,svmul_f32_z(pg0,fij,Pi_inv),gradWi.v0)),svneg_f32_z(pg0,acc.v0)));
acc.v1 = nmadd(pg0,__fkg_tmp22,gradWij.v1,madd(pg0,tmp0,madd(pg0,svmul_f32_z(pg0,fji,Pj_inv),gradWj.v1,svmul_f32_z(pg0,svmul_f32_z(pg0,fij,Pi_inv),gradWi.v1)),svneg_f32_z(pg0,acc.v1)));
acc.v2 = nmadd(pg0,__fkg_tmp22,gradWij.v2,madd(pg0,tmp0,madd(pg0,svmul_f32_z(pg0,fji,Pj_inv),gradWj.v2,svmul_f32_z(pg0,svmul_f32_z(pg0,fij,Pi_inv),gradWi.v2)),svneg_f32_z(pg0,acc.v2)));
eng = madd(pg0,svmul_f32_z(pg0,svdup_n_f32(0.5f),tmp1),madd(pg0,gradWij.v2,vij.v2,madd(pg0,gradWij.v0,vij.v0,svmul_f32_z(pg0,gradWij.v1,vij.v1))),madd(pg0,svmul_f32_z(pg0,svmul_f32_z(pg0,tmp0,fij),madd(pg0,gradWi.v2,vij.v2,madd(pg0,gradWi.v0,vij.v0,svmul_f32_z(pg0,gradWi.v1,vij.v1)))),Pi_inv,eng));
}
svst1_scatter_u32index_f32(pg0,(float*)&force[i],vindex_accv0,acc.v0);
svst1_scatter_u32index_f32(pg0,(float*)&force[i],vindex_accv1,acc.v1);
svst1_scatter_u32index_f32(pg0,(float*)&force[i],vindex_accv2,acc.v2);
svst1_scatter_u32index_f32(pg0,(float*)&force[i],vindex_eng,eng);
svst1_scatter_u32index_f32(pg0,(float*)&force[i],vindex_dt,dt);
}
}
svfloat32_t rsqrt(svbool_t pg,svfloat32_t op){
svfloat32_t rinv = svrsqrte_f32(op);
svfloat32_t h = svmul_f32_z(pg,op,rinv);
h = svmsb_n_f32_z(pg,h,rinv,1.f);
svfloat32_t poly = svmad_n_f32_z(pg,h,svdup_f32(0.375f),0.5f);
poly = svmul_f32_z(pg,poly,h);
rinv = svmad_f32_z(pg,rinv,poly,rinv);
return rinv;
}
svfloat32_t sqrt(svbool_t pg,svfloat32_t op){ return svsqrt_f32_z(pg,op); }
svfloat32_t inv(svbool_t pg,svfloat32_t op){
svfloat32_t x1 = svrecpe_f32(op);
svfloat32_t x2 = svmsb_n_f32_z(pg,op,x1,2.f);
x2 = svmul_f32_z(pg,x2,x1);
svfloat32_t ret = svmsb_n_f32_z(pg,op,x2,2.f);
ret = svmul_f32_z(pg,ret,x2);
return ret;
}
svfloat32_t  madd(svbool_t pg,svfloat32_t a,svfloat32_t b,svfloat32_t c){ return  svmad_f32_z(pg,a,b,c); }
svfloat32_t  msub(svbool_t pg,svfloat32_t a,svfloat32_t b,svfloat32_t c){ return  svmsb_f32_z(pg,a,b,c); }
svfloat32_t nmadd(svbool_t pg,svfloat32_t a,svfloat32_t b,svfloat32_t c){ return svnmad_f32_z(pg,a,b,c); }
svfloat32_t nmsub(svbool_t pg,svfloat32_t a,svfloat32_t b,svfloat32_t c){ return svnmsb_f32_z(pg,a,b,c); }
svfloat32_t max(svbool_t pg,svfloat32_t a,svfloat32_t b){ return svmax_f32_z(pg,a,b);}
svfloat32_t min(svbool_t pg,svfloat32_t a,svfloat32_t b){ return svmin_f32_z(pg,a,b);}
void transpose4x4(svfloat32x4_t& v){
const unsigned int tmp[16] = { 0, 2, 1, 3, 4, 6, 5, 7, 8,10, 9,11,12,14,13,15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
v.v0 = svtbl_f32(v.v0,index);
v.v1 = svtbl_f32(v.v1,index);
v.v2 = svtbl_f32(v.v2,index);
v.v3 = svtbl_f32(v.v3,index);
svfloat64_t xy0 = svreinterpret_f64_f32(svtrn1_f32(v.v0,v.v1));
svfloat64_t xy1 = svreinterpret_f64_f32(svtrn2_f32(v.v0,v.v1));
svfloat64_t zw0 = svreinterpret_f64_f32(svtrn1_f32(v.v2,v.v3));
svfloat64_t zw1 = svreinterpret_f64_f32(svtrn2_f32(v.v2,v.v3));
v.v0 = svreinterpret_f32_f64(svtrn1_f64(xy0,zw0));
v.v1 = svreinterpret_f32_f64(svtrn2_f64(xy0,zw0));
v.v2 = svreinterpret_f32_f64(svtrn1_f64(xy1,zw1));
v.v3 = svreinterpret_f32_f64(svtrn2_f64(xy1,zw1));
}
void gather8(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
const unsigned int tmp[16] = {  0,  8,  1,  9,  2, 10,  3, 11, 4, 12,  5, 13,  6, 14,  7, 15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
svfloat64_t a = svreinterpret_f64_f32(svtbl_f32(v0,index));
svfloat64_t b = svreinterpret_f64_f32(svtbl_f32(v1,index));
svfloat64_t c = svreinterpret_f64_f32(svtbl_f32(v2,index));
svfloat64_t d = svreinterpret_f64_f32(svtbl_f32(v3,index));
svfloat64_t e = svreinterpret_f64_f32(svtbl_f32(v4,index));
svfloat64_t f = svreinterpret_f64_f32(svtbl_f32(v5,index));
svfloat64_t g = svreinterpret_f64_f32(svtbl_f32(v6,index));
svfloat64_t h = svreinterpret_f64_f32(svtbl_f32(v7,index));
svfloat64_t ae0 = svzip1_f64(a,e);
svfloat64_t ae1 = svzip2_f64(a,e);
svfloat64_t bf0 = svzip1_f64(b,f);
svfloat64_t bf1 = svzip2_f64(b,f);
svfloat64_t cg0 = svzip1_f64(c,g);
svfloat64_t cg1 = svzip2_f64(c,g);
svfloat64_t dh0 = svzip1_f64(d,h);
svfloat64_t dh1 = svzip2_f64(d,h);
svfloat64_t aceg0 = svzip1_f64(ae0,cg0);
svfloat64_t aceg1 = svzip2_f64(ae0,cg0);
svfloat64_t aceg2 = svzip1_f64(ae1,cg1);
svfloat64_t aceg3 = svzip2_f64(ae1,cg1);
svfloat64_t bdfh0 = svzip1_f64(bf0,dh0);
svfloat64_t bdfh1 = svzip2_f64(bf0,dh0);
svfloat64_t bdfh2 = svzip1_f64(bf1,dh1);
svfloat64_t bdfh3 = svzip2_f64(bf1,dh1);
v0 = svreinterpret_f32_f64(svzip1_f64(aceg0,bdfh0));
v1 = svreinterpret_f32_f64(svzip2_f64(aceg0,bdfh0));
v2 = svreinterpret_f32_f64(svzip1_f64(aceg1,bdfh1));
v3 = svreinterpret_f32_f64(svzip2_f64(aceg1,bdfh1));
v4 = svreinterpret_f32_f64(svzip1_f64(aceg2,bdfh2));
v5 = svreinterpret_f32_f64(svzip2_f64(aceg2,bdfh2));
v6 = svreinterpret_f32_f64(svzip1_f64(aceg3,bdfh3));
v7 = svreinterpret_f32_f64(svzip2_f64(aceg3,bdfh3));
}
void gather5(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
const unsigned int tmp[16] = {0,1,2,3,4,10,11,12,5,6,7,8,9,13,14,15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
gather8(v0,v1,v2,v3,v4,v5,v6,v7);
}
void gather6(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
const unsigned int tmp[16] = {0,1,2,3,4,5,12,13,6,7,8,9,10,11,14,15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
gather8(v0,v1,v2,v3,v4,v5,v6,v7);
}
void gather7(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
const unsigned int tmp[16] = {0,1,2,3,4,5,6,14,7,8,9,10,11,12,13,15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
gather8(v0,v1,v2,v3,v4,v5,v6,v7);
}
void scatter8(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
svfloat32_t ae0 = svzip1_f32(v0,v4);
svfloat32_t ae1 = svzip2_f32(v0,v4);
svfloat32_t bf0 = svzip1_f32(v1,v5);
svfloat32_t bf1 = svzip2_f32(v1,v5);
svfloat32_t cg0 = svzip1_f32(v2,v6);
svfloat32_t cg1 = svzip2_f32(v2,v6);
svfloat32_t dh0 = svzip1_f32(v3,v7);
svfloat32_t dh1 = svzip2_f32(v3,v7);
svfloat32_t aceg0 = svzip1_f32(ae0,cg0);
svfloat32_t aceg1 = svzip2_f32(ae0,cg0);
svfloat32_t aceg2 = svzip1_f32(ae1,cg1);
svfloat32_t aceg3 = svzip2_f32(ae1,cg1);
svfloat32_t bdfh0 = svzip1_f32(bf0,dh0);
svfloat32_t bdfh1 = svzip2_f32(bf0,dh0);
svfloat32_t bdfh2 = svzip1_f32(bf1,dh1);
svfloat32_t bdfh3 = svzip2_f32(bf1,dh1);
v0 = svzip1_f32(aceg0,bdfh0);
v1 = svzip2_f32(aceg0,bdfh0);
v2 = svzip1_f32(aceg1,bdfh1);
v3 = svzip2_f32(aceg1,bdfh1);
v4 = svzip1_f32(aceg2,bdfh2);
v5 = svzip2_f32(aceg2,bdfh2);
v6 = svzip1_f32(aceg3,bdfh3);
v7 = svzip2_f32(aceg3,bdfh3);
}
void scatter5(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
scatter8(v0,v1,v2,v3,v4,v5,v6,v7);
const unsigned int tmp[16] = {0,1,2,3,4,8,9,10,11,12,5,6,7,13,14,15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
v0 = svtbl_f32(v0,index);
v1 = svtbl_f32(v1,index);
v2 = svtbl_f32(v2,index);
v3 = svtbl_f32(v3,index);
v4 = svtbl_f32(v4,index);
v5 = svtbl_f32(v5,index);
v6 = svtbl_f32(v6,index);
v7 = svtbl_f32(v7,index);
}
void scatter6(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
scatter8(v0,v1,v2,v3,v4,v5,v6,v7);
const unsigned int tmp[16] = {0,1,2,3,4,5,8,9,10,11,12,13,6,7,14,15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
v0 = svtbl_f32(v0,index);
v1 = svtbl_f32(v1,index);
v2 = svtbl_f32(v2,index);
v3 = svtbl_f32(v3,index);
v4 = svtbl_f32(v4,index);
v5 = svtbl_f32(v5,index);
v6 = svtbl_f32(v6,index);
v7 = svtbl_f32(v7,index);
}
void scatter7(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7){
scatter8(v0,v1,v2,v3,v4,v5,v6,v7);
const unsigned int tmp[16] = {0,1,2,3,4,5,6,8,9,10,11,12,13,14,7,15};
const svuint32_t index = svld1_u32(svptrue_b32(),tmp);
v0 = svtbl_f32(v0,index);
v1 = svtbl_f32(v1,index);
v2 = svtbl_f32(v2,index);
v3 = svtbl_f32(v3,index);
v4 = svtbl_f32(v4,index);
v5 = svtbl_f32(v5,index);
v6 = svtbl_f32(v6,index);
}
void transpose16x16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){
svfloat32_t ai0 = svzip1_f32(v0,v8);
svfloat32_t ai1 = svzip2_f32(v0,v8);
svfloat32_t bj0 = svzip1_f32(v1,v9);
svfloat32_t bj1 = svzip2_f32(v1,v9);
svfloat32_t ck0 = svzip1_f32(v2,v10);
svfloat32_t ck1 = svzip2_f32(v2,v10);
svfloat32_t dl0 = svzip1_f32(v3,v11);
svfloat32_t dl1 = svzip2_f32(v3,v11);
svfloat32_t em0 = svzip1_f32(v4,v12);
svfloat32_t em1 = svzip2_f32(v4,v12);
svfloat32_t fn0 = svzip1_f32(v5,v13);
svfloat32_t fn1 = svzip2_f32(v5,v13);
svfloat32_t go0 = svzip1_f32(v6,v14);
svfloat32_t go1 = svzip2_f32(v6,v14);
svfloat32_t hp0 = svzip1_f32(v7,v15);
svfloat32_t hp1 = svzip2_f32(v7,v15);
svfloat32_t aeim0 = svzip1_f32(ai0,em0);
svfloat32_t aeim1 = svzip2_f32(ai0,em0);
svfloat32_t aeim2 = svzip1_f32(ai1,em1);
svfloat32_t aeim3 = svzip2_f32(ai1,em1);
svfloat32_t bfjn0 = svzip1_f32(bj0,fn0);
svfloat32_t bfjn1 = svzip2_f32(bj0,fn0);
svfloat32_t bfjn2 = svzip1_f32(bj1,fn1);
svfloat32_t bfjn3 = svzip2_f32(bj1,fn1);
svfloat32_t cgko0 = svzip1_f32(ck0,go0);
svfloat32_t cgko1 = svzip2_f32(ck0,go0);
svfloat32_t cgko2 = svzip1_f32(ck1,go1);
svfloat32_t cgko3 = svzip2_f32(ck1,go1);
svfloat32_t dhlp0 = svzip1_f32(dl0,hp0);
svfloat32_t dhlp1 = svzip2_f32(dl0,hp0);
svfloat32_t dhlp2 = svzip1_f32(dl1,hp1);
svfloat32_t dhlp3 = svzip2_f32(dl1,hp1);
svfloat32_t acegikmo0 = svzip1_f32(aeim0,cgko0);
svfloat32_t acegikmo1 = svzip2_f32(aeim0,cgko0);
svfloat32_t acegikmo2 = svzip1_f32(aeim1,cgko1);
svfloat32_t acegikmo3 = svzip2_f32(aeim1,cgko1);
svfloat32_t acegikmo4 = svzip1_f32(aeim2,cgko2);
svfloat32_t acegikmo5 = svzip2_f32(aeim2,cgko2);
svfloat32_t acegikmo6 = svzip1_f32(aeim3,cgko3);
svfloat32_t acegikmo7 = svzip2_f32(aeim3,cgko3);
svfloat32_t bdfhjlnp0 = svzip1_f32(bfjn0,dhlp0);
svfloat32_t bdfhjlnp1 = svzip2_f32(bfjn0,dhlp0);
svfloat32_t bdfhjlnp2 = svzip1_f32(bfjn1,dhlp1);
svfloat32_t bdfhjlnp3 = svzip2_f32(bfjn1,dhlp1);
svfloat32_t bdfhjlnp4 = svzip1_f32(bfjn2,dhlp2);
svfloat32_t bdfhjlnp5 = svzip2_f32(bfjn2,dhlp2);
svfloat32_t bdfhjlnp6 = svzip1_f32(bfjn3,dhlp3);
svfloat32_t bdfhjlnp7 = svzip2_f32(bfjn3,dhlp3);
v0  = svzip1_f32(acegikmo0,bdfhjlnp0);
v1  = svzip2_f32(acegikmo0,bdfhjlnp0);
v2  = svzip1_f32(acegikmo1,bdfhjlnp1);
v3  = svzip2_f32(acegikmo1,bdfhjlnp1);
v4  = svzip1_f32(acegikmo2,bdfhjlnp2);
v5  = svzip2_f32(acegikmo2,bdfhjlnp2);
v6  = svzip1_f32(acegikmo3,bdfhjlnp3);
v7  = svzip2_f32(acegikmo3,bdfhjlnp3);
v8  = svzip1_f32(acegikmo4,bdfhjlnp4);
v9  = svzip2_f32(acegikmo4,bdfhjlnp4);
v10 = svzip1_f32(acegikmo5,bdfhjlnp5);
v11 = svzip2_f32(acegikmo5,bdfhjlnp5);
v12 = svzip1_f32(acegikmo6,bdfhjlnp6);
v13 = svzip2_f32(acegikmo6,bdfhjlnp6);
v14 = svzip1_f32(acegikmo7,bdfhjlnp7);
v15 = svzip2_f32(acegikmo7,bdfhjlnp7);
}
void gather16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){ transpose16x16(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15); };
void scatter16(svfloat32_t& v0,svfloat32_t& v1,svfloat32_t& v2,svfloat32_t& v3,svfloat32_t& v4,svfloat32_t& v5,svfloat32_t& v6,svfloat32_t& v7,svfloat32_t& v8,svfloat32_t& v9,svfloat32_t& v10,svfloat32_t& v11,svfloat32_t& v12,svfloat32_t& v13,svfloat32_t& v14,svfloat32_t& v15){ transpose16x16(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15); };
};
