#include<pikg_vector.hpp>
#include<cmath>
#include<limits>
#include<chrono>

#include <pikg_avx2.hpp>
struct CalcForceLongEPSP{
PIKG::F32 eps2;
CalcForceLongEPSP(){}
CalcForceLongEPSP(PIKG::F32 eps2):eps2(eps2){}
void initialize(PIKG::F32 eps2_){
eps2 = eps2_;
}
int kernel_id = 0;
void operator()(const EPI_t* __restrict__ epi,const int ni,const SPJ_t* __restrict__ epj,const int nj,Force_t* __restrict__ force,const int kernel_select = 1){
static_assert(sizeof(EPI_t) == 48,"check consistency of EPI member variable definition between PIKG source and original source");
static_assert(sizeof(SPJ_t) == 80,"check consistency of EPJ member variable definition between PIKG source and original source");
static_assert(sizeof(Force_t) == 32,"check consistency of FORCE member variable definition between PIKG source and original source");
if(kernel_select>=0) kernel_id = kernel_select;
if(kernel_id == 0){
std::cout << "ni: " << ni << " nj:" << nj << std::endl;
Force_t* force_tmp = new Force_t[ni];
std::chrono::system_clock::time_point  start, end;
double min_time = std::numeric_limits<double>::max();
{ // test Kernel_I8_J1
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I8_J1(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 1: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 1;
}
}
{ // test Kernel_I1_J8
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I1_J8(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 2: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 2;
}
}
delete[] force_tmp;
} // if(kernel_id == 0)
if(kernel_id == 1) Kernel_I8_J1(epi,ni,epj,nj,force);
if(kernel_id == 2) Kernel_I1_J8(epi,ni,epj,nj,force);
} // operator() definition 
void Kernel_I8_J1(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const SPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xxloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zxloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (epi[i].pos.x-epi[0].pos.x);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (epi[i].pos.y-epi[0].pos.y);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (epi[i].pos.z-epi[0].pos.z);
} // loop of i
for(j = 0;j < nj;++j){
xjloc_tmp_x[j] = (epj[j].pos.x-epi[0].pos.x);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_y[j] = (epj[j].pos.y-epi[0].pos.y);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_z[j] = (epj[j].pos.z-epi[0].pos.z);
} // loop of j
for(j = 0;j < nj;++j){
mjloc_tmp[j] = epj[j].mass;
} // loop of j
for(j = 0;j < nj;++j){
qj_xxloc_tmp[j] = epj[j].quad.xx;
} // loop of j
for(j = 0;j < nj;++j){
qj_yyloc_tmp[j] = epj[j].quad.yy;
} // loop of j
for(j = 0;j < nj;++j){
qj_zzloc_tmp[j] = epj[j].quad.zz;
} // loop of j
for(j = 0;j < nj;++j){
qj_xyloc_tmp[j] = epj[j].quad.xy;
} // loop of j
for(j = 0;j < nj;++j){
qj_yzloc_tmp[j] = epj[j].quad.yz;
} // loop of j
for(j = 0;j < nj;++j){
qj_zxloc_tmp[j] = epj[j].quad.xz;
} // loop of j
for(i = 0;i < (ni/8)*8;i += 8){
__m256x3 xiloc;

xiloc.v0 = _mm256_loadu_ps(((float*)&xiloc_tmp_x[i+0]));
xiloc.v1 = _mm256_loadu_ps(((float*)&xiloc_tmp_y[i+0]));
xiloc.v2 = _mm256_loadu_ps(((float*)&xiloc_tmp_z[i+0]));
__m256x3 af;

af.v0 = _mm256_set1_ps(0.0f);
af.v1 = _mm256_set1_ps(0.0f);
af.v2 = _mm256_set1_ps(0.0f);
__m256 pf;

pf = _mm256_set1_ps(0.0f);
for(j = 0;j < (nj/1)*1;++j){
__m256 mjloc;

mjloc = _mm256_set1_ps(mjloc_tmp[j+0]);
__m256 qj_xxloc;

qj_xxloc = _mm256_set1_ps(qj_xxloc_tmp[j+0]);
__m256 qj_xyloc;

qj_xyloc = _mm256_set1_ps(qj_xyloc_tmp[j+0]);
__m256 qj_yyloc;

qj_yyloc = _mm256_set1_ps(qj_yyloc_tmp[j+0]);
__m256 qj_yzloc;

qj_yzloc = _mm256_set1_ps(qj_yzloc_tmp[j+0]);
__m256 qj_zxloc;

qj_zxloc = _mm256_set1_ps(qj_zxloc_tmp[j+0]);
__m256 qj_zzloc;

qj_zzloc = _mm256_set1_ps(qj_zzloc_tmp[j+0]);
__m256x3 xjloc;

xjloc.v0 = _mm256_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm256_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm256_set1_ps(xjloc_tmp_z[j+0]);
__m256x3 rij;

__m256 __fkg_tmp1;

__m256 __fkg_tmp0;

__m256 r2;

__m256 r_inv;

__m256 __fkg_tmp2;

__m256 tmp;

__m256 __fkg_tmp3;

__m256 r2_inv;

__m256 r3_inv;

__m256 r4_inv;

__m256 r5_inv;

__m256 __fkg_tmp4;

__m256 tr;

__m256 qxx;

__m256 qyy;

__m256 qzz;

__m256 qxy;

__m256 qyz;

__m256 qzx;

__m256 __fkg_tmp5;

__m256 mtr;

__m256 __fkg_tmp7;

__m256 __fkg_tmp6;

__m256x3 qr;

__m256 __fkg_tmp9;

__m256 __fkg_tmp8;

__m256 __fkg_tmp11;

__m256 __fkg_tmp10;

__m256 __fkg_tmp13;

__m256 __fkg_tmp12;

__m256 rqr;

__m256 rqr_r4_inv;

__m256 meff;

__m256 __fkg_tmp14;

__m256 meff3;

__m256 __fkg_tmp15;

__m256 __fkg_tmp16;

__m256 __fkg_tmp17;

rij.v0 = _mm256_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm256_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm256_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp1 = _mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2));
__fkg_tmp0 = _mm256_fmadd_ps(rij.v1,rij.v1,__fkg_tmp1);
r2 = _mm256_fmadd_ps(rij.v2,rij.v2,__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = _mm256_mul_ps(r_inv,r_inv);
tmp = _mm256_fnmadd_ps(r2,__fkg_tmp2,_mm256_set1_ps(3.0f));
__fkg_tmp3 = _mm256_mul_ps(tmp,_mm256_set1_ps(0.5f));
r_inv = _mm256_mul_ps(r_inv,__fkg_tmp3);
r2_inv = _mm256_mul_ps(r_inv,r_inv);
r3_inv = _mm256_mul_ps(r2_inv,r_inv);
r4_inv = _mm256_mul_ps(r2_inv,r2_inv);
r5_inv = _mm256_mul_ps(r2_inv,r3_inv);
__fkg_tmp4 = _mm256_add_ps(qj_xxloc,qj_yyloc);
tr = _mm256_add_ps(__fkg_tmp4,qj_zzloc);
qxx = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_xxloc,tr);
qyy = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_yyloc,tr);
qzz = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_zzloc,tr);
qxy = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_xyloc);
qyz = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_yzloc);
qzx = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_zxloc);
__fkg_tmp5 = _mm256_mul_ps(_mm256_set1_ps(eps2),tr);
mtr = _mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),__fkg_tmp5);
__fkg_tmp7 = _mm256_mul_ps(qxx,rij.v0);
__fkg_tmp6 = _mm256_fmadd_ps(qxy,rij.v1,__fkg_tmp7);
qr.v0 = _mm256_fmadd_ps(qzx,rij.v2,__fkg_tmp6);
__fkg_tmp9 = _mm256_mul_ps(qyy,rij.v1);
__fkg_tmp8 = _mm256_fmadd_ps(qyz,rij.v2,__fkg_tmp9);
qr.v1 = _mm256_fmadd_ps(qxy,rij.v0,__fkg_tmp8);
__fkg_tmp11 = _mm256_mul_ps(qzz,rij.v2);
__fkg_tmp10 = _mm256_fmadd_ps(qzx,rij.v0,__fkg_tmp11);
qr.v2 = _mm256_fmadd_ps(qyz,rij.v1,__fkg_tmp10);
__fkg_tmp13 = _mm256_fmadd_ps(qr.v0,rij.v0,mtr);
__fkg_tmp12 = _mm256_fmadd_ps(qr.v1,rij.v1,__fkg_tmp13);
rqr = _mm256_fmadd_ps(qr.v2,rij.v2,__fkg_tmp12);
rqr_r4_inv = _mm256_mul_ps(rqr,r4_inv);
meff = _mm256_fmadd_ps(_mm256_set1_ps(0.5f),rqr_r4_inv,mjloc);
__fkg_tmp14 = _mm256_fmadd_ps(_mm256_set1_ps(2.5f),rqr_r4_inv,mjloc);
meff3 = _mm256_mul_ps(__fkg_tmp14,r3_inv);
pf = _mm256_fnmadd_ps(meff,r_inv,pf);
__fkg_tmp15 = _mm256_fnmadd_ps(r5_inv,qr.v0,af.v0);
af.v0 = _mm256_fmadd_ps(meff3,rij.v0,__fkg_tmp15);
__fkg_tmp16 = _mm256_fnmadd_ps(r5_inv,qr.v1,af.v1);
af.v1 = _mm256_fmadd_ps(meff3,rij.v1,__fkg_tmp16);
__fkg_tmp17 = _mm256_fnmadd_ps(r5_inv,qr.v2,af.v2);
af.v2 = _mm256_fmadd_ps(meff3,rij.v2,__fkg_tmp17);
} // loop of j

{
__m256 __fkg_tmp_accum;
alignas(32) int index_gather_load0[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load0 = _mm256_load_si256((const __m256i*)index_gather_load0);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.x),vindex_gather_load0,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,af.v0);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc.x)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc.x)[8] = __fkg_store_tmp[1];
((float*)&force[i+0].acc.x)[16] = __fkg_store_tmp[2];
((float*)&force[i+0].acc.x)[24] = __fkg_store_tmp[3];
((float*)&force[i+0].acc.x)[32] = __fkg_store_tmp[4];
((float*)&force[i+0].acc.x)[40] = __fkg_store_tmp[5];
((float*)&force[i+0].acc.x)[48] = __fkg_store_tmp[6];
((float*)&force[i+0].acc.x)[56] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
alignas(32) int index_gather_load1[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load1 = _mm256_load_si256((const __m256i*)index_gather_load1);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.y),vindex_gather_load1,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,af.v1);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc.y)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc.y)[8] = __fkg_store_tmp[1];
((float*)&force[i+0].acc.y)[16] = __fkg_store_tmp[2];
((float*)&force[i+0].acc.y)[24] = __fkg_store_tmp[3];
((float*)&force[i+0].acc.y)[32] = __fkg_store_tmp[4];
((float*)&force[i+0].acc.y)[40] = __fkg_store_tmp[5];
((float*)&force[i+0].acc.y)[48] = __fkg_store_tmp[6];
((float*)&force[i+0].acc.y)[56] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
alignas(32) int index_gather_load2[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load2 = _mm256_load_si256((const __m256i*)index_gather_load2);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.z),vindex_gather_load2,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,af.v2);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc.z)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc.z)[8] = __fkg_store_tmp[1];
((float*)&force[i+0].acc.z)[16] = __fkg_store_tmp[2];
((float*)&force[i+0].acc.z)[24] = __fkg_store_tmp[3];
((float*)&force[i+0].acc.z)[32] = __fkg_store_tmp[4];
((float*)&force[i+0].acc.z)[40] = __fkg_store_tmp[5];
((float*)&force[i+0].acc.z)[48] = __fkg_store_tmp[6];
((float*)&force[i+0].acc.z)[56] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
alignas(32) int index_gather_load3[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load3 = _mm256_load_si256((const __m256i*)index_gather_load3);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].phi),vindex_gather_load3,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,pf);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].phi)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].phi)[8] = __fkg_store_tmp[1];
((float*)&force[i+0].phi)[16] = __fkg_store_tmp[2];
((float*)&force[i+0].phi)[24] = __fkg_store_tmp[3];
((float*)&force[i+0].phi)[32] = __fkg_store_tmp[4];
((float*)&force[i+0].phi)[40] = __fkg_store_tmp[5];
((float*)&force[i+0].phi)[48] = __fkg_store_tmp[6];
((float*)&force[i+0].phi)[56] = __fkg_store_tmp[7];
}
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec af;

af.x = 0.0f;
af.y = 0.0f;
af.z = 0.0f;
PIKG::F32 pf;

pf = 0.0f;
for(j = 0;j < nj;++j){
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::F32 qj_xxloc;

qj_xxloc = qj_xxloc_tmp[j+0];
PIKG::F32 qj_xyloc;

qj_xyloc = qj_xyloc_tmp[j+0];
PIKG::F32 qj_yyloc;

qj_yyloc = qj_yyloc_tmp[j+0];
PIKG::F32 qj_yzloc;

qj_yzloc = qj_yzloc_tmp[j+0];
PIKG::F32 qj_zxloc;

qj_zxloc = qj_zxloc_tmp[j+0];
PIKG::F32 qj_zzloc;

qj_zzloc = qj_zzloc_tmp[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32vec rij;

PIKG::F32 __fkg_tmp1;

PIKG::F32 __fkg_tmp0;

PIKG::F32 r2;

PIKG::F32 r_inv;

PIKG::F32 __fkg_tmp2;

PIKG::F32 tmp;

PIKG::F32 __fkg_tmp3;

PIKG::F32 r2_inv;

PIKG::F32 r3_inv;

PIKG::F32 r4_inv;

PIKG::F32 r5_inv;

PIKG::F32 __fkg_tmp4;

PIKG::F32 tr;

PIKG::F32 qxx;

PIKG::F32 qyy;

PIKG::F32 qzz;

PIKG::F32 qxy;

PIKG::F32 qyz;

PIKG::F32 qzx;

PIKG::F32 __fkg_tmp5;

PIKG::F32 mtr;

PIKG::F32 __fkg_tmp7;

PIKG::F32 __fkg_tmp6;

PIKG::F32vec qr;

PIKG::F32 __fkg_tmp9;

PIKG::F32 __fkg_tmp8;

PIKG::F32 __fkg_tmp11;

PIKG::F32 __fkg_tmp10;

PIKG::F32 __fkg_tmp13;

PIKG::F32 __fkg_tmp12;

PIKG::F32 rqr;

PIKG::F32 rqr_r4_inv;

PIKG::F32 meff;

PIKG::F32 __fkg_tmp14;

PIKG::F32 meff3;

PIKG::F32 __fkg_tmp15;

PIKG::F32 __fkg_tmp16;

PIKG::F32 __fkg_tmp17;

rij.x = (xjloc.x-xiloc.x);
rij.y = (xjloc.y-xiloc.y);
rij.z = (xjloc.z-xiloc.z);
__fkg_tmp1 = (rij.x*rij.x+eps2);
__fkg_tmp0 = (rij.y*rij.y+__fkg_tmp1);
r2 = (rij.z*rij.z+__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = (r_inv*r_inv);
tmp = (3.0f - r2*__fkg_tmp2);
__fkg_tmp3 = (tmp*0.5f);
r_inv = (r_inv*__fkg_tmp3);
r2_inv = (r_inv*r_inv);
r3_inv = (r2_inv*r_inv);
r4_inv = (r2_inv*r2_inv);
r5_inv = (r2_inv*r3_inv);
__fkg_tmp4 = (qj_xxloc+qj_yyloc);
tr = (__fkg_tmp4+qj_zzloc);
qxx = (3.0f*qj_xxloc-tr);
qyy = (3.0f*qj_yyloc-tr);
qzz = (3.0f*qj_zzloc-tr);
qxy = (3.0f*qj_xyloc);
qyz = (3.0f*qj_yzloc);
qzx = (3.0f*qj_zxloc);
__fkg_tmp5 = (eps2*tr);
mtr = -(__fkg_tmp5);
__fkg_tmp7 = (qxx*rij.x);
__fkg_tmp6 = (qxy*rij.y+__fkg_tmp7);
qr.x = (qzx*rij.z+__fkg_tmp6);
__fkg_tmp9 = (qyy*rij.y);
__fkg_tmp8 = (qyz*rij.z+__fkg_tmp9);
qr.y = (qxy*rij.x+__fkg_tmp8);
__fkg_tmp11 = (qzz*rij.z);
__fkg_tmp10 = (qzx*rij.x+__fkg_tmp11);
qr.z = (qyz*rij.y+__fkg_tmp10);
__fkg_tmp13 = (qr.x*rij.x+mtr);
__fkg_tmp12 = (qr.y*rij.y+__fkg_tmp13);
rqr = (qr.z*rij.z+__fkg_tmp12);
rqr_r4_inv = (rqr*r4_inv);
meff = (0.5f*rqr_r4_inv+mjloc);
__fkg_tmp14 = (2.5f*rqr_r4_inv+mjloc);
meff3 = (__fkg_tmp14*r3_inv);
pf = (pf - meff*r_inv);
__fkg_tmp15 = (af.x - r5_inv*qr.x);
af.x = (meff3*rij.x+__fkg_tmp15);
__fkg_tmp16 = (af.y - r5_inv*qr.y);
af.y = (meff3*rij.y+__fkg_tmp16);
__fkg_tmp17 = (af.z - r5_inv*qr.z);
af.z = (meff3*rij.z+__fkg_tmp17);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+af.x);
force[i+0].acc.y = (force[i+0].acc.y+af.y);
force[i+0].acc.z = (force[i+0].acc.z+af.z);
force[i+0].phi = (force[i+0].phi+pf);
} // loop of i
} // end loop of reference 
} // Kernel_I8_J1 definition 
void Kernel_I1_J8(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const SPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xxloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_xyloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_yzloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) qj_zxloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (epi[i].pos.x-epi[0].pos.x);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (epi[i].pos.y-epi[0].pos.y);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (epi[i].pos.z-epi[0].pos.z);
} // loop of i
for(j = 0;j < nj;++j){
xjloc_tmp_x[j] = (epj[j].pos.x-epi[0].pos.x);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_y[j] = (epj[j].pos.y-epi[0].pos.y);
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_z[j] = (epj[j].pos.z-epi[0].pos.z);
} // loop of j
for(j = 0;j < nj;++j){
mjloc_tmp[j] = epj[j].mass;
} // loop of j
for(j = 0;j < nj;++j){
qj_xxloc_tmp[j] = epj[j].quad.xx;
} // loop of j
for(j = 0;j < nj;++j){
qj_yyloc_tmp[j] = epj[j].quad.yy;
} // loop of j
for(j = 0;j < nj;++j){
qj_zzloc_tmp[j] = epj[j].quad.zz;
} // loop of j
for(j = 0;j < nj;++j){
qj_xyloc_tmp[j] = epj[j].quad.xy;
} // loop of j
for(j = 0;j < nj;++j){
qj_yzloc_tmp[j] = epj[j].quad.yz;
} // loop of j
for(j = 0;j < nj;++j){
qj_zxloc_tmp[j] = epj[j].quad.xz;
} // loop of j
for(i = 0;i < (ni/1)*1;++i){
__m256x3 xiloc;

xiloc.v0 = _mm256_set1_ps(xiloc_tmp_x[i+0]);
xiloc.v1 = _mm256_set1_ps(xiloc_tmp_y[i+0]);
xiloc.v2 = _mm256_set1_ps(xiloc_tmp_z[i+0]);
__m256x3 af;

af.v0 = _mm256_set1_ps(0.0f);
af.v1 = _mm256_set1_ps(0.0f);
af.v2 = _mm256_set1_ps(0.0f);
__m256 pf;

pf = _mm256_set1_ps(0.0f);
for(j = 0;j < (nj/8)*8;j += 8){
__m256 mjloc;

mjloc = _mm256_loadu_ps(((float*)&mjloc_tmp[j+0]));
__m256 qj_xxloc;

qj_xxloc = _mm256_loadu_ps(((float*)&qj_xxloc_tmp[j+0]));
__m256 qj_xyloc;

qj_xyloc = _mm256_loadu_ps(((float*)&qj_xyloc_tmp[j+0]));
__m256 qj_yyloc;

qj_yyloc = _mm256_loadu_ps(((float*)&qj_yyloc_tmp[j+0]));
__m256 qj_yzloc;

qj_yzloc = _mm256_loadu_ps(((float*)&qj_yzloc_tmp[j+0]));
__m256 qj_zxloc;

qj_zxloc = _mm256_loadu_ps(((float*)&qj_zxloc_tmp[j+0]));
__m256 qj_zzloc;

qj_zzloc = _mm256_loadu_ps(((float*)&qj_zzloc_tmp[j+0]));
__m256x3 xjloc;

xjloc.v0 = _mm256_loadu_ps(((float*)&xjloc_tmp_x[j+0]));
xjloc.v1 = _mm256_loadu_ps(((float*)&xjloc_tmp_y[j+0]));
xjloc.v2 = _mm256_loadu_ps(((float*)&xjloc_tmp_z[j+0]));
__m256x3 rij;

__m256 __fkg_tmp1;

__m256 __fkg_tmp0;

__m256 r2;

__m256 r_inv;

__m256 __fkg_tmp2;

__m256 tmp;

__m256 __fkg_tmp3;

__m256 r2_inv;

__m256 r3_inv;

__m256 r4_inv;

__m256 r5_inv;

__m256 __fkg_tmp4;

__m256 tr;

__m256 qxx;

__m256 qyy;

__m256 qzz;

__m256 qxy;

__m256 qyz;

__m256 qzx;

__m256 __fkg_tmp5;

__m256 mtr;

__m256 __fkg_tmp7;

__m256 __fkg_tmp6;

__m256x3 qr;

__m256 __fkg_tmp9;

__m256 __fkg_tmp8;

__m256 __fkg_tmp11;

__m256 __fkg_tmp10;

__m256 __fkg_tmp13;

__m256 __fkg_tmp12;

__m256 rqr;

__m256 rqr_r4_inv;

__m256 meff;

__m256 __fkg_tmp14;

__m256 meff3;

__m256 __fkg_tmp15;

__m256 __fkg_tmp16;

__m256 __fkg_tmp17;

rij.v0 = _mm256_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm256_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm256_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp1 = _mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2));
__fkg_tmp0 = _mm256_fmadd_ps(rij.v1,rij.v1,__fkg_tmp1);
r2 = _mm256_fmadd_ps(rij.v2,rij.v2,__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = _mm256_mul_ps(r_inv,r_inv);
tmp = _mm256_fnmadd_ps(r2,__fkg_tmp2,_mm256_set1_ps(3.0f));
__fkg_tmp3 = _mm256_mul_ps(tmp,_mm256_set1_ps(0.5f));
r_inv = _mm256_mul_ps(r_inv,__fkg_tmp3);
r2_inv = _mm256_mul_ps(r_inv,r_inv);
r3_inv = _mm256_mul_ps(r2_inv,r_inv);
r4_inv = _mm256_mul_ps(r2_inv,r2_inv);
r5_inv = _mm256_mul_ps(r2_inv,r3_inv);
__fkg_tmp4 = _mm256_add_ps(qj_xxloc,qj_yyloc);
tr = _mm256_add_ps(__fkg_tmp4,qj_zzloc);
qxx = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_xxloc,tr);
qyy = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_yyloc,tr);
qzz = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_zzloc,tr);
qxy = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_xyloc);
qyz = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_yzloc);
qzx = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_zxloc);
__fkg_tmp5 = _mm256_mul_ps(_mm256_set1_ps(eps2),tr);
mtr = _mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),__fkg_tmp5);
__fkg_tmp7 = _mm256_mul_ps(qxx,rij.v0);
__fkg_tmp6 = _mm256_fmadd_ps(qxy,rij.v1,__fkg_tmp7);
qr.v0 = _mm256_fmadd_ps(qzx,rij.v2,__fkg_tmp6);
__fkg_tmp9 = _mm256_mul_ps(qyy,rij.v1);
__fkg_tmp8 = _mm256_fmadd_ps(qyz,rij.v2,__fkg_tmp9);
qr.v1 = _mm256_fmadd_ps(qxy,rij.v0,__fkg_tmp8);
__fkg_tmp11 = _mm256_mul_ps(qzz,rij.v2);
__fkg_tmp10 = _mm256_fmadd_ps(qzx,rij.v0,__fkg_tmp11);
qr.v2 = _mm256_fmadd_ps(qyz,rij.v1,__fkg_tmp10);
__fkg_tmp13 = _mm256_fmadd_ps(qr.v0,rij.v0,mtr);
__fkg_tmp12 = _mm256_fmadd_ps(qr.v1,rij.v1,__fkg_tmp13);
rqr = _mm256_fmadd_ps(qr.v2,rij.v2,__fkg_tmp12);
rqr_r4_inv = _mm256_mul_ps(rqr,r4_inv);
meff = _mm256_fmadd_ps(_mm256_set1_ps(0.5f),rqr_r4_inv,mjloc);
__fkg_tmp14 = _mm256_fmadd_ps(_mm256_set1_ps(2.5f),rqr_r4_inv,mjloc);
meff3 = _mm256_mul_ps(__fkg_tmp14,r3_inv);
pf = _mm256_fnmadd_ps(meff,r_inv,pf);
__fkg_tmp15 = _mm256_fnmadd_ps(r5_inv,qr.v0,af.v0);
af.v0 = _mm256_fmadd_ps(meff3,rij.v0,__fkg_tmp15);
__fkg_tmp16 = _mm256_fnmadd_ps(r5_inv,qr.v1,af.v1);
af.v1 = _mm256_fmadd_ps(meff3,rij.v1,__fkg_tmp16);
__fkg_tmp17 = _mm256_fnmadd_ps(r5_inv,qr.v2,af.v2);
af.v2 = _mm256_fmadd_ps(meff3,rij.v2,__fkg_tmp17);
} // loop of j

if(j<nj){ // tail j loop
__m256x3 __fkg_tmp18;

__fkg_tmp18.v0 = af.v0;
__fkg_tmp18.v1 = af.v1;
__fkg_tmp18.v2 = af.v2;
__m256 __fkg_tmp19;

__fkg_tmp19 = pf;
for(;j < nj;++j){
__m256 mjloc;

mjloc = _mm256_set1_ps(mjloc_tmp[j+0]);
__m256 qj_xxloc;

qj_xxloc = _mm256_set1_ps(qj_xxloc_tmp[j+0]);
__m256 qj_xyloc;

qj_xyloc = _mm256_set1_ps(qj_xyloc_tmp[j+0]);
__m256 qj_yyloc;

qj_yyloc = _mm256_set1_ps(qj_yyloc_tmp[j+0]);
__m256 qj_yzloc;

qj_yzloc = _mm256_set1_ps(qj_yzloc_tmp[j+0]);
__m256 qj_zxloc;

qj_zxloc = _mm256_set1_ps(qj_zxloc_tmp[j+0]);
__m256 qj_zzloc;

qj_zzloc = _mm256_set1_ps(qj_zzloc_tmp[j+0]);
__m256x3 xjloc;

xjloc.v0 = _mm256_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm256_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm256_set1_ps(xjloc_tmp_z[j+0]);
__m256x3 rij;

__m256 __fkg_tmp1;

__m256 __fkg_tmp0;

__m256 r2;

__m256 r_inv;

__m256 __fkg_tmp2;

__m256 tmp;

__m256 __fkg_tmp3;

__m256 r2_inv;

__m256 r3_inv;

__m256 r4_inv;

__m256 r5_inv;

__m256 __fkg_tmp4;

__m256 tr;

__m256 qxx;

__m256 qyy;

__m256 qzz;

__m256 qxy;

__m256 qyz;

__m256 qzx;

__m256 __fkg_tmp5;

__m256 mtr;

__m256 __fkg_tmp7;

__m256 __fkg_tmp6;

__m256x3 qr;

__m256 __fkg_tmp9;

__m256 __fkg_tmp8;

__m256 __fkg_tmp11;

__m256 __fkg_tmp10;

__m256 __fkg_tmp13;

__m256 __fkg_tmp12;

__m256 rqr;

__m256 rqr_r4_inv;

__m256 meff;

__m256 __fkg_tmp14;

__m256 meff3;

__m256 __fkg_tmp15;

__m256 __fkg_tmp16;

__m256 __fkg_tmp17;

rij.v0 = _mm256_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm256_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm256_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp1 = _mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2));
__fkg_tmp0 = _mm256_fmadd_ps(rij.v1,rij.v1,__fkg_tmp1);
r2 = _mm256_fmadd_ps(rij.v2,rij.v2,__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = _mm256_mul_ps(r_inv,r_inv);
tmp = _mm256_fnmadd_ps(r2,__fkg_tmp2,_mm256_set1_ps(3.0f));
__fkg_tmp3 = _mm256_mul_ps(tmp,_mm256_set1_ps(0.5f));
r_inv = _mm256_mul_ps(r_inv,__fkg_tmp3);
r2_inv = _mm256_mul_ps(r_inv,r_inv);
r3_inv = _mm256_mul_ps(r2_inv,r_inv);
r4_inv = _mm256_mul_ps(r2_inv,r2_inv);
r5_inv = _mm256_mul_ps(r2_inv,r3_inv);
__fkg_tmp4 = _mm256_add_ps(qj_xxloc,qj_yyloc);
tr = _mm256_add_ps(__fkg_tmp4,qj_zzloc);
qxx = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_xxloc,tr);
qyy = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_yyloc,tr);
qzz = _mm256_fmsub_ps(_mm256_set1_ps(3.0f),qj_zzloc,tr);
qxy = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_xyloc);
qyz = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_yzloc);
qzx = _mm256_mul_ps(_mm256_set1_ps(3.0f),qj_zxloc);
__fkg_tmp5 = _mm256_mul_ps(_mm256_set1_ps(eps2),tr);
mtr = _mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),__fkg_tmp5);
__fkg_tmp7 = _mm256_mul_ps(qxx,rij.v0);
__fkg_tmp6 = _mm256_fmadd_ps(qxy,rij.v1,__fkg_tmp7);
qr.v0 = _mm256_fmadd_ps(qzx,rij.v2,__fkg_tmp6);
__fkg_tmp9 = _mm256_mul_ps(qyy,rij.v1);
__fkg_tmp8 = _mm256_fmadd_ps(qyz,rij.v2,__fkg_tmp9);
qr.v1 = _mm256_fmadd_ps(qxy,rij.v0,__fkg_tmp8);
__fkg_tmp11 = _mm256_mul_ps(qzz,rij.v2);
__fkg_tmp10 = _mm256_fmadd_ps(qzx,rij.v0,__fkg_tmp11);
qr.v2 = _mm256_fmadd_ps(qyz,rij.v1,__fkg_tmp10);
__fkg_tmp13 = _mm256_fmadd_ps(qr.v0,rij.v0,mtr);
__fkg_tmp12 = _mm256_fmadd_ps(qr.v1,rij.v1,__fkg_tmp13);
rqr = _mm256_fmadd_ps(qr.v2,rij.v2,__fkg_tmp12);
rqr_r4_inv = _mm256_mul_ps(rqr,r4_inv);
meff = _mm256_fmadd_ps(_mm256_set1_ps(0.5f),rqr_r4_inv,mjloc);
__fkg_tmp14 = _mm256_fmadd_ps(_mm256_set1_ps(2.5f),rqr_r4_inv,mjloc);
meff3 = _mm256_mul_ps(__fkg_tmp14,r3_inv);
pf = _mm256_fnmadd_ps(meff,r_inv,pf);
__fkg_tmp15 = _mm256_fnmadd_ps(r5_inv,qr.v0,af.v0);
af.v0 = _mm256_fmadd_ps(meff3,rij.v0,__fkg_tmp15);
__fkg_tmp16 = _mm256_fnmadd_ps(r5_inv,qr.v1,af.v1);
af.v1 = _mm256_fmadd_ps(meff3,rij.v1,__fkg_tmp16);
__fkg_tmp17 = _mm256_fnmadd_ps(r5_inv,qr.v2,af.v2);
af.v2 = _mm256_fmadd_ps(meff3,rij.v2,__fkg_tmp17);
} // loop of j
af.v0 = _mm256_blend_ps(__fkg_tmp18.v0,af.v0,0b00000001);
af.v1 = _mm256_blend_ps(__fkg_tmp18.v1,af.v1,0b00000001);
af.v2 = _mm256_blend_ps(__fkg_tmp18.v2,af.v2,0b00000001);
pf = _mm256_blend_ps(__fkg_tmp19,pf,0b00000001);
} // if of j tail loop

{
af.v0 = _mm256_add_ps(af.v0,_mm256_shuffle_ps(af.v0,af.v0,0xb1));
af.v0 = _mm256_add_ps(af.v0,_mm256_shuffle_ps(af.v0,af.v0,0xee));
af.v0 = _mm256_add_ps(af.v0,_mm256_castps128_ps256(_mm256_extractf128_ps(af.v0,1)));
((float*)&force[i+0].acc.x)[0] = (((float*)&force[i+0].acc.x)[0]+af.v0[0]);
}

{
af.v1 = _mm256_add_ps(af.v1,_mm256_shuffle_ps(af.v1,af.v1,0xb1));
af.v1 = _mm256_add_ps(af.v1,_mm256_shuffle_ps(af.v1,af.v1,0xee));
af.v1 = _mm256_add_ps(af.v1,_mm256_castps128_ps256(_mm256_extractf128_ps(af.v1,1)));
((float*)&force[i+0].acc.y)[0] = (((float*)&force[i+0].acc.y)[0]+af.v1[0]);
}

{
af.v2 = _mm256_add_ps(af.v2,_mm256_shuffle_ps(af.v2,af.v2,0xb1));
af.v2 = _mm256_add_ps(af.v2,_mm256_shuffle_ps(af.v2,af.v2,0xee));
af.v2 = _mm256_add_ps(af.v2,_mm256_castps128_ps256(_mm256_extractf128_ps(af.v2,1)));
((float*)&force[i+0].acc.z)[0] = (((float*)&force[i+0].acc.z)[0]+af.v2[0]);
}

{
pf = _mm256_add_ps(pf,_mm256_shuffle_ps(pf,pf,0xb1));
pf = _mm256_add_ps(pf,_mm256_shuffle_ps(pf,pf,0xee));
pf = _mm256_add_ps(pf,_mm256_castps128_ps256(_mm256_extractf128_ps(pf,1)));
((float*)&force[i+0].phi)[0] = (((float*)&force[i+0].phi)[0]+pf[0]);
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec af;

af.x = 0.0f;
af.y = 0.0f;
af.z = 0.0f;
PIKG::F32 pf;

pf = 0.0f;
for(j = 0;j < nj;++j){
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::F32 qj_xxloc;

qj_xxloc = qj_xxloc_tmp[j+0];
PIKG::F32 qj_xyloc;

qj_xyloc = qj_xyloc_tmp[j+0];
PIKG::F32 qj_yyloc;

qj_yyloc = qj_yyloc_tmp[j+0];
PIKG::F32 qj_yzloc;

qj_yzloc = qj_yzloc_tmp[j+0];
PIKG::F32 qj_zxloc;

qj_zxloc = qj_zxloc_tmp[j+0];
PIKG::F32 qj_zzloc;

qj_zzloc = qj_zzloc_tmp[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32vec rij;

PIKG::F32 __fkg_tmp1;

PIKG::F32 __fkg_tmp0;

PIKG::F32 r2;

PIKG::F32 r_inv;

PIKG::F32 __fkg_tmp2;

PIKG::F32 tmp;

PIKG::F32 __fkg_tmp3;

PIKG::F32 r2_inv;

PIKG::F32 r3_inv;

PIKG::F32 r4_inv;

PIKG::F32 r5_inv;

PIKG::F32 __fkg_tmp4;

PIKG::F32 tr;

PIKG::F32 qxx;

PIKG::F32 qyy;

PIKG::F32 qzz;

PIKG::F32 qxy;

PIKG::F32 qyz;

PIKG::F32 qzx;

PIKG::F32 __fkg_tmp5;

PIKG::F32 mtr;

PIKG::F32 __fkg_tmp7;

PIKG::F32 __fkg_tmp6;

PIKG::F32vec qr;

PIKG::F32 __fkg_tmp9;

PIKG::F32 __fkg_tmp8;

PIKG::F32 __fkg_tmp11;

PIKG::F32 __fkg_tmp10;

PIKG::F32 __fkg_tmp13;

PIKG::F32 __fkg_tmp12;

PIKG::F32 rqr;

PIKG::F32 rqr_r4_inv;

PIKG::F32 meff;

PIKG::F32 __fkg_tmp14;

PIKG::F32 meff3;

PIKG::F32 __fkg_tmp15;

PIKG::F32 __fkg_tmp16;

PIKG::F32 __fkg_tmp17;

rij.x = (xjloc.x-xiloc.x);
rij.y = (xjloc.y-xiloc.y);
rij.z = (xjloc.z-xiloc.z);
__fkg_tmp1 = (rij.x*rij.x+eps2);
__fkg_tmp0 = (rij.y*rij.y+__fkg_tmp1);
r2 = (rij.z*rij.z+__fkg_tmp0);
r_inv = rsqrt(r2);
__fkg_tmp2 = (r_inv*r_inv);
tmp = (3.0f - r2*__fkg_tmp2);
__fkg_tmp3 = (tmp*0.5f);
r_inv = (r_inv*__fkg_tmp3);
r2_inv = (r_inv*r_inv);
r3_inv = (r2_inv*r_inv);
r4_inv = (r2_inv*r2_inv);
r5_inv = (r2_inv*r3_inv);
__fkg_tmp4 = (qj_xxloc+qj_yyloc);
tr = (__fkg_tmp4+qj_zzloc);
qxx = (3.0f*qj_xxloc-tr);
qyy = (3.0f*qj_yyloc-tr);
qzz = (3.0f*qj_zzloc-tr);
qxy = (3.0f*qj_xyloc);
qyz = (3.0f*qj_yzloc);
qzx = (3.0f*qj_zxloc);
__fkg_tmp5 = (eps2*tr);
mtr = -(__fkg_tmp5);
__fkg_tmp7 = (qxx*rij.x);
__fkg_tmp6 = (qxy*rij.y+__fkg_tmp7);
qr.x = (qzx*rij.z+__fkg_tmp6);
__fkg_tmp9 = (qyy*rij.y);
__fkg_tmp8 = (qyz*rij.z+__fkg_tmp9);
qr.y = (qxy*rij.x+__fkg_tmp8);
__fkg_tmp11 = (qzz*rij.z);
__fkg_tmp10 = (qzx*rij.x+__fkg_tmp11);
qr.z = (qyz*rij.y+__fkg_tmp10);
__fkg_tmp13 = (qr.x*rij.x+mtr);
__fkg_tmp12 = (qr.y*rij.y+__fkg_tmp13);
rqr = (qr.z*rij.z+__fkg_tmp12);
rqr_r4_inv = (rqr*r4_inv);
meff = (0.5f*rqr_r4_inv+mjloc);
__fkg_tmp14 = (2.5f*rqr_r4_inv+mjloc);
meff3 = (__fkg_tmp14*r3_inv);
pf = (pf - meff*r_inv);
__fkg_tmp15 = (af.x - r5_inv*qr.x);
af.x = (meff3*rij.x+__fkg_tmp15);
__fkg_tmp16 = (af.y - r5_inv*qr.y);
af.y = (meff3*rij.y+__fkg_tmp16);
__fkg_tmp17 = (af.z - r5_inv*qr.z);
af.z = (meff3*rij.z+__fkg_tmp17);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+af.x);
force[i+0].acc.y = (force[i+0].acc.y+af.y);
force[i+0].acc.z = (force[i+0].acc.z+af.z);
force[i+0].phi = (force[i+0].phi+pf);
} // loop of i
} // end loop of reference 
} // Kernel_I1_J8 definition 
PIKG::F64 rsqrt(PIKG::F64 op){ return 1.0/std::sqrt(op); }
PIKG::F64 sqrt(PIKG::F64 op){ return std::sqrt(op); }
PIKG::F64 inv(PIKG::F64 op){ return 1.0/op; }
PIKG::F64 max(PIKG::F64 a,PIKG::F64 b){ return std::max(a,b);}
PIKG::F64 min(PIKG::F64 a,PIKG::F64 b){ return std::min(a,b);}
PIKG::F32 rsqrt(PIKG::F32 op){ return 1.f/std::sqrt(op); }
PIKG::F32 sqrt(PIKG::F32 op){ return std::sqrt(op); }
PIKG::F32 inv(PIKG::F32 op){ return 1.f/op; }
PIKG::S64 max(PIKG::S64 a,PIKG::S64 b){ return std::max(a,b);}
PIKG::S64 min(PIKG::S64 a,PIKG::S64 b){ return std::min(a,b);}
PIKG::S32 max(PIKG::S32 a,PIKG::S32 b){ return std::max(a,b);}
PIKG::S32 min(PIKG::S32 a,PIKG::S32 b){ return std::min(a,b);}
PIKG::F64 table(PIKG::F64 tab[],PIKG::S64 i){ return tab[i]; }
PIKG::F32 table(PIKG::F32 tab[],PIKG::S32 i){ return tab[i]; }
PIKG::F64 to_float(PIKG::U64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::U32 op){return (PIKG::F32)op;}
PIKG::F64 to_float(PIKG::S64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::S32 op){return (PIKG::F32)op;}
PIKG::S64   to_int(PIKG::F64 op){return (PIKG::S64)op;}
PIKG::S32   to_int(PIKG::F32 op){return (PIKG::S32)op;}
PIKG::U64  to_uint(PIKG::F64 op){return (PIKG::U64)op;}
PIKG::U32  to_uint(PIKG::F32 op){return (PIKG::U32)op;}
template<typename T> PIKG::F64 to_f64(const T& op){return (PIKG::F64)op;}
template<typename T> PIKG::F32 to_f32(const T& op){return (PIKG::F32)op;}
template<typename T> PIKG::S64 to_s64(const T& op){return (PIKG::S64)op;}
template<typename T> PIKG::S32 to_s32(const T& op){return (PIKG::S32)op;}
template<typename T> PIKG::U64 to_u64(const T& op){return (PIKG::U64)op;}
template<typename T> PIKG::U32 to_u32(const T& op){return (PIKG::U32)op;}
__m256 rsqrt(__m256 op){
  return _mm256_rsqrt_ps(op);
}
__m256 sqrt(__m256 op){ return _mm256_sqrt_ps(op); }
__m256 inv(__m256 op){
return _mm256_rcp_ps(op);
}
__m256d rsqrt(__m256d op){
  __m256d y = _mm256_castsi256_pd(_mm256_sub_epi64(_mm256_set1_epi64x(0x5fe6eb50c7b537a9LL),_mm256_srlv_epi64(_mm256_castpd_si256(op),_mm256_set1_epi64x(1))));
  __m256d h = _mm256_mul_pd(op,y);
  h = _mm256_fnmadd_pd(h,y,_mm256_set1_pd(1.0));
  __m256d poly = _mm256_fmadd_pd(h,_mm256_set1_pd(0.375),_mm256_set1_pd(0.5));
  poly = _mm256_mul_pd(poly,h);
  y = _mm256_fmadd_pd(y,poly,y);
  return y;
}__m256d sqrt(__m256d op){
  return _mm256_sqrt_pd(op);
}__m256d inv(__m256d op){
  __m256 x = _mm256_castps128_ps256(_mm256_cvtpd_ps(op));
  x = inv(x);
  return _mm256_cvtps_pd(_mm256_castps256_ps128(x));
}__m256d max(__m256d a,__m256d b){ return _mm256_max_pd(a,b);}
__m256d min(__m256d a,__m256d b){ return _mm256_min_pd(a,b);}
__m256  max(__m256  a,__m256  b){ return _mm256_max_ps(a,b);}
__m256  min(__m256  a,__m256  b){ return _mm256_min_ps(a,b);}
__m256i max(__m256i a,__m256i b){ return _mm256_max_epi32(a,b);}
__m256i min(__m256i a,__m256i b){ return _mm256_min_epi32(a,b);}
__m256d table(__m256d tab,__m256i index){ return _mm256_permutexvar_pd(index,tab);}
__m256  table(__m256  tab,__m256i index){ return _mm256_permutexvar_ps(index,tab);}
__m256  to_float(__m256i op){ return _mm256_cvtepi32_ps(op);}
__m256i  to_int(__m256  op){ return _mm256_cvtps_epi32(op);}
void transpose2x2_pd(__m256d& a, __m256d& b){
  __m256d tmp = _mm256_unpacklo_pd(a,b);
  b = _mm256_unpackhi_pd(a,b);
  a = tmp;
}
void transpose4x4_pd(__m256d& a, __m256d& b,__m256d& c,__m256d& d){
  __m256d tmp0 = _mm256_unpacklo_pd(a, b);
  __m256d tmp1 = _mm256_unpackhi_pd(a, b);
  __m256d tmp2 = _mm256_unpacklo_pd(c, d);
  __m256d tmp3 = _mm256_unpackhi_pd(c, d);
  a = _mm256_permute2f128_pd(tmp0, tmp2, 0|(2<<4));
  b = _mm256_permute2f128_pd(tmp1, tmp3, 0|(2<<4));
  c = _mm256_permute2f128_pd(tmp0, tmp2, 1|(3<<4));
  d = _mm256_permute2f128_pd(tmp1, tmp3, 1|(3<<4));
}
void unpack2x2_pd(__m256d& a,__m256d& b){
  transpose2x2_pd(a,b);
}
void pack2x2_pd(__m256d& a,__m256d& b){
  transpose2x2_pd(a,b);
}
void unpack4x4_pd(__m256d& a,__m256d& b,__m256d& c,__m256d& d){
  transpose4x4_pd(a,b,c,d);
}
void pack4x4_pd(__m256d& a,__m256d& b,__m256d& c,__m256d& d){
  transpose4x4_pd(a,b,c,d);
}
void unpack2x2_ps(__m256& a, __m256& b){
  __m256 tmp = _mm256_shuffle_ps(a,b,0xdd);
  b = _mm256_shuffle_ps(a,b,0x88);
  a = tmp;
}
void pack2x2_ps(__m256& a, __m256& b){
  __m256 tmp = _mm256_unpackhi_ps(a,b);
  b = _mm256_unpacklo_ps(a,b);
  a = tmp;
}
void gather4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  __m256 src0 = _mm256_permute2f128_ps(a,c,0x20); // x0 y0 z0 w0 x4 y4 z4 w4
  __m256 src1 = _mm256_permute2f128_ps(a,c,0x31); // x1 y1 z1 w1 x5 y5 z5 w5
  __m256 src2 = _mm256_permute2f128_ps(b,d,0x20);
  __m256 src3 = _mm256_permute2f128_ps(b,d,0x31);
  __m256d tmp0 = _mm256_castps_pd(_mm256_unpacklo_ps(src0, src1)); // x0 x2 y0 y2 x1 y1 x3 y3
  __m256d tmp1 = _mm256_castps_pd(_mm256_unpackhi_ps(src0, src1)); // z0 z2 w0 w2 z1 w1 z3 w3
  __m256d tmp2 = _mm256_castps_pd(_mm256_unpacklo_ps(src2, src3)); // x4 x6 y4 y6 x5 x7 y5 y7
  __m256d tmp3 = _mm256_castps_pd(_mm256_unpackhi_ps(src2, src3)); // z4 z6 w4 w6 z5 z7 w5 w7
  a = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp0, tmp2)); // x0 x2 x4 x6 x1 x3 x5 x7
  b = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp0, tmp2)); // y0 y2 y4 y6 y1 y3 y5 y7
  c = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp1, tmp3)); // z0 z2 z4 z6 z1 z3 z5 z7
  d = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp1, tmp3)); // w0 w2 w4 w6 w1 w3 w5 w7
}
void scatter4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  __m256d tmp0 = _mm256_castps_pd(_mm256_unpacklo_ps(a, b)); // x0 x2 y0 y2 x1 y1 x3 y3
  __m256d tmp1 = _mm256_castps_pd(_mm256_unpackhi_ps(a, b)); // z0 z2 w0 w2 z1 w1 z3 w3
  __m256d tmp2 = _mm256_castps_pd(_mm256_unpacklo_ps(c, d)); // x4 x6 y4 y6 x5 x7 y5 y7
  __m256d tmp3 = _mm256_castps_pd(_mm256_unpackhi_ps(c, d)); // z4 z6 w4 w6 z5 z7 w5 w7
  __m256 dst0 = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp0, tmp2)); // x0 x2 x4 x6 x1 x3 x5 x7
  __m256 dst1 = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp0, tmp2)); // y0 y2 y4 y6 y1 y3 y5 y7
  __m256 dst2 = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp1, tmp3)); // z0 z2 z4 z6 z1 z3 z5 z7
  __m256 dst3 = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp1, tmp3)); // w0 w2 w4 w6 w1 w3 w5 w7
  a = _mm256_permute2f128_ps(dst0,dst1,0x20);
  b = _mm256_permute2f128_ps(dst2,dst3,0x20);
  c = _mm256_permute2f128_ps(dst0,dst1,0x31);
  d = _mm256_permute2f128_ps(dst2,dst3,0x31);
}
void unpack4x4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  gather4_ps(a,b,c,d);
}
void pack4x4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  scatter4_ps(a,b,c,d);
}
void unpack4x4_ps(__m256x4& v){
  unpack4x4_ps(v.v0,v.v1,v.v2,v.v3);
}
void pack4x4_ps(__m256x4& v){
  unpack4x4_ps(v.v0,v.v1,v.v2,v.v3);
}
void unpack4x4_pd(__m256dx4& v){
  unpack4x4_pd(v.v0,v.v1,v.v2,v.v3);
}
void pack4x4_pd(__m256dx4& v){
  unpack4x4_pd(v.v0,v.v1,v.v2,v.v3);
}
void print_ps(const __m256 v){ for(int i=0;i<8;i++) printf(" %f",v[i]); printf("\n");}
void print_ps(const __m256x4 v){
  print_ps(v.v0);
  print_ps(v.v1);
  print_ps(v.v2);
  print_ps(v.v3);
}
};// kernel functor definition 
