#include<pikg_vector.hpp>
#include<cmath>
#include<limits>
#include<chrono>

struct CalcForceLongEPEP{
PIKG::F32 eps2;
CalcForceLongEPEP(){}
CalcForceLongEPEP(PIKG::F32 eps2):eps2(eps2){}
void initialize(PIKG::F32 eps2_){
eps2 = eps2_;
}
int kernel_id = 0;
void operator()(const EPI_t* __restrict__ epi,const int ni,const EPJ_t* __restrict__ epj,const int nj,Force_t* __restrict__ force,const int kernel_select = 1){
static_assert(sizeof(EPI_t) == 48,"check consistency of EPI member variable definition between PIKG source and original source");
static_assert(sizeof(EPJ_t) == 112,"check consistency of EPJ member variable definition between PIKG source and original source");
static_assert(sizeof(Force_t) == 32,"check consistency of FORCE member variable definition between PIKG source and original source");
if(kernel_select>=0) kernel_id = kernel_select;
if(kernel_id == 0){
std::cout << "ni: " << ni << " nj:" << nj << std::endl;
Force_t* force_tmp = new Force_t[ni];
std::chrono::system_clock::time_point  start, end;
double min_time = std::numeric_limits<double>::max();
{ // test Kernel_I1_J1
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I1_J1(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 1: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 1;
}
}
delete[] force_tmp;
} // if(kernel_id == 0)
if(kernel_id == 1) Kernel_I1_J1(epi,ni,epj,nj,force);
} // operator() definition 
void Kernel_I1_J1(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const EPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) rsearchiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) routiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) rsearchjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) routjloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (epi[i].pos.x-epi[0].pos.x);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (epi[i].pos.y-epi[0].pos.y);
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (epi[i].pos.z-epi[0].pos.z);
} // loop of i
for(i = 0;i < ni;++i){
rsearchiloc_tmp[i] = epi[i].r_search;
} // loop of i
for(i = 0;i < ni;++i){
routiloc_tmp[i] = epi[i].r_out;
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
rsearchjloc_tmp[j] = epj[j].r_search;
} // loop of j
for(j = 0;j < nj;++j){
routjloc_tmp[j] = epj[j].r_out;
} // loop of j
for(i = 0;i < ni;++i){
PIKG::S32 idi;

idi = epi[i+0].id_local;
PIKG::S32 ranki;

ranki = epi[i+0].myrank;
PIKG::F32 routiloc;

routiloc = routiloc_tmp[i+0];
PIKG::F32 rsearchiloc;

rsearchiloc = rsearchiloc_tmp[i+0];
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec af;

af.x = 0.0f;
af.y = 0.0f;
af.z = 0.0f;
PIKG::S32 id_maxf;

id_maxf = std::numeric_limits<int32_t>::lowest();
PIKG::S32 id_minf;

id_minf = std::numeric_limits<int32_t>::max();
PIKG::S32 ngbf;

ngbf = 0;
PIKG::F32 pf;

pf = 0.0f;
PIKG::S32 rankf;

rankf = 0;
for(j = 0;j < nj;++j){
PIKG::S32 idj;

idj = epj[j].id_local;
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::S32 rankj;

rankj = epj[j].myrank;
PIKG::F32 routjloc;

routjloc = routjloc_tmp[j+0];
PIKG::F32 rsearchjloc;

rsearchjloc = rsearchjloc_tmp[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32 rout;

PIKG::F32 rsearch;

PIKG::F32 rout2;

PIKG::F32 __fkg_tmp4;

PIKG::F32 rsearch2;

PIKG::F32vec rij;

PIKG::F32 __fkg_tmp6;

PIKG::F32 __fkg_tmp5;

PIKG::F32 r2_real;

PIKG::F32 r2;

PIKG::S32 rank_sub;

PIKG::S32 rank_sub2;

PIKG::S32 __fkg_tmp0;

PIKG::S32 __fkg_tmp1;

PIKG::S32 __fkg_tmp2;

PIKG::S32 __fkg_tmp3;

PIKG::F32 r_inv;

PIKG::F32 __fkg_tmp7;

PIKG::F32 tmp;

PIKG::F32 __fkg_tmp8;

PIKG::F32 r2_inv;

PIKG::F32 mr_inv;

PIKG::F32 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = (rout*rout);
__fkg_tmp4 = (rsearch*rsearch);
rsearch2 = (__fkg_tmp4*1.0201f);
rij.x = (xjloc.x-xiloc.x);
rij.y = (xjloc.y-xiloc.y);
rij.z = (xjloc.z-xiloc.z);
__fkg_tmp6 = (rij.x*rij.x+eps2);
__fkg_tmp5 = (rij.y*rij.y+__fkg_tmp6);
r2_real = (rij.z*rij.z+__fkg_tmp5);
r2 = max(r2_real,rout2);
rank_sub = (ranki-rankj);
rank_sub2 = (rank_sub*rank_sub);
if((r2_real<rsearch2)){
if(((idi!=idj)||(ranki!=rankj))){
__fkg_tmp0 = (ngbf+1);
__fkg_tmp1 = (rankf+rank_sub2);
__fkg_tmp2 = max(id_maxf,idj);
__fkg_tmp3 = min(id_minf,idj);
ngbf = __fkg_tmp0;
rankf = __fkg_tmp1;
id_maxf = __fkg_tmp2;
id_minf = __fkg_tmp3;
}
}
r_inv = rsqrt(r2);
__fkg_tmp7 = (r_inv*r_inv);
tmp = (3.0f - r2*__fkg_tmp7);
__fkg_tmp8 = (tmp*0.5f);
r_inv = (r_inv*__fkg_tmp8);
r2_inv = (r_inv*r_inv);
mr_inv = (mjloc*r_inv);
mr3_inv = (r2_inv*mr_inv);
af.x = (mr3_inv*rij.x+af.x);
af.y = (mr3_inv*rij.y+af.y);
af.z = (mr3_inv*rij.z+af.z);
pf = (pf-mr_inv);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+af.x);
force[i+0].acc.y = (force[i+0].acc.y+af.y);
force[i+0].acc.z = (force[i+0].acc.z+af.z);
force[i+0].neighbor.id_max = max(id_maxf,force[i+0].neighbor.id_max);
force[i+0].neighbor.id_min = min(id_minf,force[i+0].neighbor.id_min);
force[i+0].neighbor.number = (force[i+0].neighbor.number+ngbf);
force[i+0].phi = (force[i+0].phi+pf);
force[i+0].neighbor.rank = (force[i+0].neighbor.rank+rankf);
} // loop of i
} // Kernel_I1_J1 definition 
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
};// kernel functor definition 
