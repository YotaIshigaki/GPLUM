#include<pikg_vector.hpp>
#include<cmath>
#include<limits>
#include<chrono>

#include <pikg_avx2.hpp>
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
void Kernel_I8_J1(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const EPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
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
for(i = 0;i < (ni/8)*8;i += 8){
__m256i idi;

alignas(32) int index_gather_load0[8] = {0,12,24,36,48,60,72,84};
__m256i vindex_gather_load0 = _mm256_load_si256((const __m256i*)index_gather_load0);
idi = _mm256_i32gather_epi32(((int*)&epi[i+0].id_local),vindex_gather_load0,4);
__m256i ranki;

alignas(32) int index_gather_load1[8] = {0,12,24,36,48,60,72,84};
__m256i vindex_gather_load1 = _mm256_load_si256((const __m256i*)index_gather_load1);
ranki = _mm256_i32gather_epi32(((int*)&epi[i+0].myrank),vindex_gather_load1,4);
__m256 routiloc;

routiloc = _mm256_loadu_ps(((float*)&routiloc_tmp[i+0]));
__m256 rsearchiloc;

rsearchiloc = _mm256_loadu_ps(((float*)&rsearchiloc_tmp[i+0]));
__m256x3 xiloc;

xiloc.v0 = _mm256_loadu_ps(((float*)&xiloc_tmp_x[i+0]));
xiloc.v1 = _mm256_loadu_ps(((float*)&xiloc_tmp_y[i+0]));
xiloc.v2 = _mm256_loadu_ps(((float*)&xiloc_tmp_z[i+0]));
__m256x3 af;

af.v0 = _mm256_set1_ps(0.0f);
af.v1 = _mm256_set1_ps(0.0f);
af.v2 = _mm256_set1_ps(0.0f);
__m256i id_maxf;

id_maxf = _mm256_set1_epi32(std::numeric_limits<int32_t>::lowest());
__m256i id_minf;

id_minf = _mm256_set1_epi32(std::numeric_limits<int32_t>::max());
__m256i ngbf;

ngbf = _mm256_set1_epi32(0);
__m256 pf;

pf = _mm256_set1_ps(0.0f);
__m256i rankf;

rankf = _mm256_set1_epi32(0);
for(j = 0;j < (nj/1)*1;++j){
__m256i idj;

idj = _mm256_set1_epi32(epj[j].id_local);

__m256 mjloc;

mjloc = _mm256_set1_ps(mjloc_tmp[j+0]);
__m256i rankj;

rankj = _mm256_set1_epi32(epj[j].myrank);

__m256 routjloc;

routjloc = _mm256_set1_ps(routjloc_tmp[j+0]);
__m256 rsearchjloc;

rsearchjloc = _mm256_set1_ps(rsearchjloc_tmp[j+0]);
__m256x3 xjloc;

xjloc.v0 = _mm256_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm256_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm256_set1_ps(xjloc_tmp_z[j+0]);
__m256 rout;

__m256 rsearch;

__m256 rout2;

__m256 __fkg_tmp4;

__m256 rsearch2;

__m256x3 rij;

__m256 __fkg_tmp6;

__m256 __fkg_tmp5;

__m256 r2_real;

__m256 r2;

__m256i rank_sub;

__m256i rank_sub2;

__m256i __fkg_tmp0;

__m256i __fkg_tmp1;

__m256i __fkg_tmp2;

__m256i __fkg_tmp3;

__m256 r_inv;

__m256 __fkg_tmp7;

__m256 tmp;

__m256 __fkg_tmp8;

__m256 r2_inv;

__m256 mr_inv;

__m256 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = _mm256_mul_ps(rout,rout);
__fkg_tmp4 = _mm256_mul_ps(rsearch,rsearch);
rsearch2 = _mm256_mul_ps(__fkg_tmp4,_mm256_set1_ps(1.0201f));
rij.v0 = _mm256_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm256_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm256_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp6 = _mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2));
__fkg_tmp5 = _mm256_fmadd_ps(rij.v1,rij.v1,__fkg_tmp6);
r2_real = _mm256_fmadd_ps(rij.v2,rij.v2,__fkg_tmp5);
r2 = max(r2_real,rout2);
rank_sub = _mm256_sub_epi32(ranki,rankj);
rank_sub2 = _mm256_mul_epi32(rank_sub,rank_sub);
{
__m256 pg1;
__m256 pg0;
pg1 = _mm256_cmp_ps(r2_real,rsearch2,_CMP_LT_OS);
pg0 = pg1;

{
__m256 pg3;
__m256 pg2;
pg3 = _mm256_or_ps(_mm256_cmp_ps(_mm256_castsi256_ps(idi),_mm256_castsi256_ps(idj),_CMP_NEQ_OQ),_mm256_cmp_ps(_mm256_castsi256_ps(ranki),_mm256_castsi256_ps(rankj),_CMP_NEQ_OQ));
pg2 = pg3;
pg3 = _mm256_and_ps(pg3,pg1);

__fkg_tmp0 = _mm256_add_epi32(ngbf,_mm256_set1_epi32(1));
__fkg_tmp1 = _mm256_add_epi32(rankf,rank_sub2);
__fkg_tmp2 = max(id_maxf,idj);
__fkg_tmp3 = min(id_minf,idj);
ngbf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(ngbf),_mm256_castsi256_ps(__fkg_tmp0),pg3));;
rankf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(rankf),_mm256_castsi256_ps(__fkg_tmp1),pg3));;
id_maxf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(id_maxf),_mm256_castsi256_ps(__fkg_tmp2),pg3));;
id_minf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(id_minf),_mm256_castsi256_ps(__fkg_tmp3),pg3));;
}

}

r_inv = rsqrt(r2);
__fkg_tmp7 = _mm256_mul_ps(r_inv,r_inv);
tmp = _mm256_fnmadd_ps(r2,__fkg_tmp7,_mm256_set1_ps(3.0f));
__fkg_tmp8 = _mm256_mul_ps(tmp,_mm256_set1_ps(0.5f));
r_inv = _mm256_mul_ps(r_inv,__fkg_tmp8);
r2_inv = _mm256_mul_ps(r_inv,r_inv);
mr_inv = _mm256_mul_ps(mjloc,r_inv);
mr3_inv = _mm256_mul_ps(r2_inv,mr_inv);
af.v0 = _mm256_fmadd_ps(mr3_inv,rij.v0,af.v0);
af.v1 = _mm256_fmadd_ps(mr3_inv,rij.v1,af.v1);
af.v2 = _mm256_fmadd_ps(mr3_inv,rij.v2,af.v2);
pf = _mm256_sub_ps(pf,mr_inv);
} // loop of j

{
__m256 __fkg_tmp_accum;
alignas(32) int index_gather_load2[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load2 = _mm256_load_si256((const __m256i*)index_gather_load2);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.x),vindex_gather_load2,4);
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
alignas(32) int index_gather_load3[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load3 = _mm256_load_si256((const __m256i*)index_gather_load3);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.y),vindex_gather_load3,4);
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
alignas(32) int index_gather_load4[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load4 = _mm256_load_si256((const __m256i*)index_gather_load4);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.z),vindex_gather_load4,4);
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
__m256i __fkg_tmp_accum;
alignas(32) int index_gather_load5[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load5 = _mm256_load_si256((const __m256i*)index_gather_load5);
__fkg_tmp_accum = _mm256_i32gather_epi32(((int*)&force[i+0].neighbor.id_max),vindex_gather_load5,4);
__fkg_tmp_accum = max(__fkg_tmp_accum,id_maxf);
{
PIKG::S32 __fkg_store_tmp[8];
_mm256_storeu_si256((__m256i*)__fkg_store_tmp,__fkg_tmp_accum);
((int*)&force[i+0].neighbor.id_max)[0] = __fkg_store_tmp[0];
((int*)&force[i+0].neighbor.id_max)[8] = __fkg_store_tmp[1];
((int*)&force[i+0].neighbor.id_max)[16] = __fkg_store_tmp[2];
((int*)&force[i+0].neighbor.id_max)[24] = __fkg_store_tmp[3];
((int*)&force[i+0].neighbor.id_max)[32] = __fkg_store_tmp[4];
((int*)&force[i+0].neighbor.id_max)[40] = __fkg_store_tmp[5];
((int*)&force[i+0].neighbor.id_max)[48] = __fkg_store_tmp[6];
((int*)&force[i+0].neighbor.id_max)[56] = __fkg_store_tmp[7];
}
}

{
__m256i __fkg_tmp_accum;
alignas(32) int index_gather_load6[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load6 = _mm256_load_si256((const __m256i*)index_gather_load6);
__fkg_tmp_accum = _mm256_i32gather_epi32(((int*)&force[i+0].neighbor.id_min),vindex_gather_load6,4);
__fkg_tmp_accum = min(__fkg_tmp_accum,id_minf);
{
PIKG::S32 __fkg_store_tmp[8];
_mm256_storeu_si256((__m256i*)__fkg_store_tmp,__fkg_tmp_accum);
((int*)&force[i+0].neighbor.id_min)[0] = __fkg_store_tmp[0];
((int*)&force[i+0].neighbor.id_min)[8] = __fkg_store_tmp[1];
((int*)&force[i+0].neighbor.id_min)[16] = __fkg_store_tmp[2];
((int*)&force[i+0].neighbor.id_min)[24] = __fkg_store_tmp[3];
((int*)&force[i+0].neighbor.id_min)[32] = __fkg_store_tmp[4];
((int*)&force[i+0].neighbor.id_min)[40] = __fkg_store_tmp[5];
((int*)&force[i+0].neighbor.id_min)[48] = __fkg_store_tmp[6];
((int*)&force[i+0].neighbor.id_min)[56] = __fkg_store_tmp[7];
}
}

{
__m256i __fkg_tmp_accum;
alignas(32) int index_gather_load7[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load7 = _mm256_load_si256((const __m256i*)index_gather_load7);
__fkg_tmp_accum = _mm256_i32gather_epi32(((int*)&force[i+0].neighbor.number),vindex_gather_load7,4);
__fkg_tmp_accum = _mm256_add_epi32(__fkg_tmp_accum,ngbf);
{
PIKG::S32 __fkg_store_tmp[8];
_mm256_storeu_si256((__m256i*)__fkg_store_tmp,__fkg_tmp_accum);
((int*)&force[i+0].neighbor.number)[0] = __fkg_store_tmp[0];
((int*)&force[i+0].neighbor.number)[8] = __fkg_store_tmp[1];
((int*)&force[i+0].neighbor.number)[16] = __fkg_store_tmp[2];
((int*)&force[i+0].neighbor.number)[24] = __fkg_store_tmp[3];
((int*)&force[i+0].neighbor.number)[32] = __fkg_store_tmp[4];
((int*)&force[i+0].neighbor.number)[40] = __fkg_store_tmp[5];
((int*)&force[i+0].neighbor.number)[48] = __fkg_store_tmp[6];
((int*)&force[i+0].neighbor.number)[56] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
alignas(32) int index_gather_load8[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load8 = _mm256_load_si256((const __m256i*)index_gather_load8);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].phi),vindex_gather_load8,4);
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

{
__m256i __fkg_tmp_accum;
alignas(32) int index_gather_load9[8] = {0,8,16,24,32,40,48,56};
__m256i vindex_gather_load9 = _mm256_load_si256((const __m256i*)index_gather_load9);
__fkg_tmp_accum = _mm256_i32gather_epi32(((int*)&force[i+0].neighbor.rank),vindex_gather_load9,4);
__fkg_tmp_accum = _mm256_add_epi32(__fkg_tmp_accum,rankf);
{
PIKG::S32 __fkg_store_tmp[8];
_mm256_storeu_si256((__m256i*)__fkg_store_tmp,__fkg_tmp_accum);
((int*)&force[i+0].neighbor.rank)[0] = __fkg_store_tmp[0];
((int*)&force[i+0].neighbor.rank)[8] = __fkg_store_tmp[1];
((int*)&force[i+0].neighbor.rank)[16] = __fkg_store_tmp[2];
((int*)&force[i+0].neighbor.rank)[24] = __fkg_store_tmp[3];
((int*)&force[i+0].neighbor.rank)[32] = __fkg_store_tmp[4];
((int*)&force[i+0].neighbor.rank)[40] = __fkg_store_tmp[5];
((int*)&force[i+0].neighbor.rank)[48] = __fkg_store_tmp[6];
((int*)&force[i+0].neighbor.rank)[56] = __fkg_store_tmp[7];
}
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
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
} // end loop of reference 
} // Kernel_I8_J1 definition 
void Kernel_I1_J8(const EPI_t* __restrict__ epi,const PIKG::S32 ni,const EPJ_t* __restrict__ epj,const PIKG::S32 nj,Force_t* __restrict__ force){
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
for(i = 0;i < (ni/1)*1;++i){
__m256i idi;

idi = _mm256_set1_epi32(epi[i+0].id_local);

__m256i ranki;

ranki = _mm256_set1_epi32(epi[i+0].myrank);

__m256 routiloc;

routiloc = _mm256_set1_ps(routiloc_tmp[i+0]);
__m256 rsearchiloc;

rsearchiloc = _mm256_set1_ps(rsearchiloc_tmp[i+0]);
__m256x3 xiloc;

xiloc.v0 = _mm256_set1_ps(xiloc_tmp_x[i+0]);
xiloc.v1 = _mm256_set1_ps(xiloc_tmp_y[i+0]);
xiloc.v2 = _mm256_set1_ps(xiloc_tmp_z[i+0]);
__m256x3 af;

af.v0 = _mm256_set1_ps(0.0f);
af.v1 = _mm256_set1_ps(0.0f);
af.v2 = _mm256_set1_ps(0.0f);
__m256i id_maxf;

id_maxf = _mm256_set1_epi32(std::numeric_limits<int32_t>::lowest());
__m256i id_minf;

id_minf = _mm256_set1_epi32(std::numeric_limits<int32_t>::max());
__m256i ngbf;

ngbf = _mm256_set1_epi32(0);
__m256 pf;

pf = _mm256_set1_ps(0.0f);
__m256i rankf;

rankf = _mm256_set1_epi32(0);
for(j = 0;j < (nj/8)*8;j += 8){
__m256i idj;

alignas(32) int index_gather_load10[8] = {0,28,56,84,112,140,168,196};
__m256i vindex_gather_load10 = _mm256_load_si256((const __m256i*)index_gather_load10);
idj = _mm256_i32gather_epi32(((int*)&epj[j].id_local),vindex_gather_load10,4);
__m256 mjloc;

mjloc = _mm256_loadu_ps(((float*)&mjloc_tmp[j+0]));
__m256i rankj;

alignas(32) int index_gather_load11[8] = {0,28,56,84,112,140,168,196};
__m256i vindex_gather_load11 = _mm256_load_si256((const __m256i*)index_gather_load11);
rankj = _mm256_i32gather_epi32(((int*)&epj[j].myrank),vindex_gather_load11,4);
__m256 routjloc;

routjloc = _mm256_loadu_ps(((float*)&routjloc_tmp[j+0]));
__m256 rsearchjloc;

rsearchjloc = _mm256_loadu_ps(((float*)&rsearchjloc_tmp[j+0]));
__m256x3 xjloc;

xjloc.v0 = _mm256_loadu_ps(((float*)&xjloc_tmp_x[j+0]));
xjloc.v1 = _mm256_loadu_ps(((float*)&xjloc_tmp_y[j+0]));
xjloc.v2 = _mm256_loadu_ps(((float*)&xjloc_tmp_z[j+0]));
__m256 rout;

__m256 rsearch;

__m256 rout2;

__m256 __fkg_tmp4;

__m256 rsearch2;

__m256x3 rij;

__m256 __fkg_tmp6;

__m256 __fkg_tmp5;

__m256 r2_real;

__m256 r2;

__m256i rank_sub;

__m256i rank_sub2;

__m256i __fkg_tmp0;

__m256i __fkg_tmp1;

__m256i __fkg_tmp2;

__m256i __fkg_tmp3;

__m256 r_inv;

__m256 __fkg_tmp7;

__m256 tmp;

__m256 __fkg_tmp8;

__m256 r2_inv;

__m256 mr_inv;

__m256 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = _mm256_mul_ps(rout,rout);
__fkg_tmp4 = _mm256_mul_ps(rsearch,rsearch);
rsearch2 = _mm256_mul_ps(__fkg_tmp4,_mm256_set1_ps(1.0201f));
rij.v0 = _mm256_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm256_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm256_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp6 = _mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2));
__fkg_tmp5 = _mm256_fmadd_ps(rij.v1,rij.v1,__fkg_tmp6);
r2_real = _mm256_fmadd_ps(rij.v2,rij.v2,__fkg_tmp5);
r2 = max(r2_real,rout2);
rank_sub = _mm256_sub_epi32(ranki,rankj);
rank_sub2 = _mm256_mul_epi32(rank_sub,rank_sub);
{
__m256 pg1;
__m256 pg0;
pg1 = _mm256_cmp_ps(r2_real,rsearch2,_CMP_LT_OS);
pg0 = pg1;

{
__m256 pg3;
__m256 pg2;
pg3 = _mm256_or_ps(_mm256_cmp_ps(_mm256_castsi256_ps(idi),_mm256_castsi256_ps(idj),_CMP_NEQ_OQ),_mm256_cmp_ps(_mm256_castsi256_ps(ranki),_mm256_castsi256_ps(rankj),_CMP_NEQ_OQ));
pg2 = pg3;
pg3 = _mm256_and_ps(pg3,pg1);

__fkg_tmp0 = _mm256_add_epi32(ngbf,_mm256_set1_epi32(1));
__fkg_tmp1 = _mm256_add_epi32(rankf,rank_sub2);
__fkg_tmp2 = max(id_maxf,idj);
__fkg_tmp3 = min(id_minf,idj);
ngbf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(ngbf),_mm256_castsi256_ps(__fkg_tmp0),pg3));;
rankf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(rankf),_mm256_castsi256_ps(__fkg_tmp1),pg3));;
id_maxf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(id_maxf),_mm256_castsi256_ps(__fkg_tmp2),pg3));;
id_minf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(id_minf),_mm256_castsi256_ps(__fkg_tmp3),pg3));;
}

}

r_inv = rsqrt(r2);
__fkg_tmp7 = _mm256_mul_ps(r_inv,r_inv);
tmp = _mm256_fnmadd_ps(r2,__fkg_tmp7,_mm256_set1_ps(3.0f));
__fkg_tmp8 = _mm256_mul_ps(tmp,_mm256_set1_ps(0.5f));
r_inv = _mm256_mul_ps(r_inv,__fkg_tmp8);
r2_inv = _mm256_mul_ps(r_inv,r_inv);
mr_inv = _mm256_mul_ps(mjloc,r_inv);
mr3_inv = _mm256_mul_ps(r2_inv,mr_inv);
af.v0 = _mm256_fmadd_ps(mr3_inv,rij.v0,af.v0);
af.v1 = _mm256_fmadd_ps(mr3_inv,rij.v1,af.v1);
af.v2 = _mm256_fmadd_ps(mr3_inv,rij.v2,af.v2);
pf = _mm256_sub_ps(pf,mr_inv);
} // loop of j

if(j<nj){ // tail j loop
__m256x3 __fkg_tmp9;

__fkg_tmp9.v0 = af.v0;
__fkg_tmp9.v1 = af.v1;
__fkg_tmp9.v2 = af.v2;
__m256i __fkg_tmp10;

__fkg_tmp10 = id_maxf;
__m256i __fkg_tmp11;

__fkg_tmp11 = id_minf;
__m256i __fkg_tmp12;

__fkg_tmp12 = ngbf;
__m256 __fkg_tmp13;

__fkg_tmp13 = pf;
__m256i __fkg_tmp14;

__fkg_tmp14 = rankf;
for(;j < nj;++j){
__m256i idj;

idj = _mm256_set1_epi32(epj[j].id_local);

__m256 mjloc;

mjloc = _mm256_set1_ps(mjloc_tmp[j+0]);
__m256i rankj;

rankj = _mm256_set1_epi32(epj[j].myrank);

__m256 routjloc;

routjloc = _mm256_set1_ps(routjloc_tmp[j+0]);
__m256 rsearchjloc;

rsearchjloc = _mm256_set1_ps(rsearchjloc_tmp[j+0]);
__m256x3 xjloc;

xjloc.v0 = _mm256_set1_ps(xjloc_tmp_x[j+0]);
xjloc.v1 = _mm256_set1_ps(xjloc_tmp_y[j+0]);
xjloc.v2 = _mm256_set1_ps(xjloc_tmp_z[j+0]);
__m256 rout;

__m256 rsearch;

__m256 rout2;

__m256 __fkg_tmp4;

__m256 rsearch2;

__m256x3 rij;

__m256 __fkg_tmp6;

__m256 __fkg_tmp5;

__m256 r2_real;

__m256 r2;

__m256i rank_sub;

__m256i rank_sub2;

__m256i __fkg_tmp0;

__m256i __fkg_tmp1;

__m256i __fkg_tmp2;

__m256i __fkg_tmp3;

__m256 r_inv;

__m256 __fkg_tmp7;

__m256 tmp;

__m256 __fkg_tmp8;

__m256 r2_inv;

__m256 mr_inv;

__m256 mr3_inv;

rout = max(routjloc,routiloc);
rsearch = max(rsearchjloc,rsearchiloc);
rout2 = _mm256_mul_ps(rout,rout);
__fkg_tmp4 = _mm256_mul_ps(rsearch,rsearch);
rsearch2 = _mm256_mul_ps(__fkg_tmp4,_mm256_set1_ps(1.0201f));
rij.v0 = _mm256_sub_ps(xjloc.v0,xiloc.v0);
rij.v1 = _mm256_sub_ps(xjloc.v1,xiloc.v1);
rij.v2 = _mm256_sub_ps(xjloc.v2,xiloc.v2);
__fkg_tmp6 = _mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2));
__fkg_tmp5 = _mm256_fmadd_ps(rij.v1,rij.v1,__fkg_tmp6);
r2_real = _mm256_fmadd_ps(rij.v2,rij.v2,__fkg_tmp5);
r2 = max(r2_real,rout2);
rank_sub = _mm256_sub_epi32(ranki,rankj);
rank_sub2 = _mm256_mul_epi32(rank_sub,rank_sub);
{
__m256 pg1;
__m256 pg0;
pg1 = _mm256_cmp_ps(r2_real,rsearch2,_CMP_LT_OS);
pg0 = pg1;

{
__m256 pg3;
__m256 pg2;
pg3 = _mm256_or_ps(_mm256_cmp_ps(_mm256_castsi256_ps(idi),_mm256_castsi256_ps(idj),_CMP_NEQ_OQ),_mm256_cmp_ps(_mm256_castsi256_ps(ranki),_mm256_castsi256_ps(rankj),_CMP_NEQ_OQ));
pg2 = pg3;
pg3 = _mm256_and_ps(pg3,pg1);

__fkg_tmp0 = _mm256_add_epi32(ngbf,_mm256_set1_epi32(1));
__fkg_tmp1 = _mm256_add_epi32(rankf,rank_sub2);
__fkg_tmp2 = max(id_maxf,idj);
__fkg_tmp3 = min(id_minf,idj);
ngbf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(ngbf),_mm256_castsi256_ps(__fkg_tmp0),pg3));;
rankf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(rankf),_mm256_castsi256_ps(__fkg_tmp1),pg3));;
id_maxf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(id_maxf),_mm256_castsi256_ps(__fkg_tmp2),pg3));;
id_minf = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(id_minf),_mm256_castsi256_ps(__fkg_tmp3),pg3));;
}

}

r_inv = rsqrt(r2);
__fkg_tmp7 = _mm256_mul_ps(r_inv,r_inv);
tmp = _mm256_fnmadd_ps(r2,__fkg_tmp7,_mm256_set1_ps(3.0f));
__fkg_tmp8 = _mm256_mul_ps(tmp,_mm256_set1_ps(0.5f));
r_inv = _mm256_mul_ps(r_inv,__fkg_tmp8);
r2_inv = _mm256_mul_ps(r_inv,r_inv);
mr_inv = _mm256_mul_ps(mjloc,r_inv);
mr3_inv = _mm256_mul_ps(r2_inv,mr_inv);
af.v0 = _mm256_fmadd_ps(mr3_inv,rij.v0,af.v0);
af.v1 = _mm256_fmadd_ps(mr3_inv,rij.v1,af.v1);
af.v2 = _mm256_fmadd_ps(mr3_inv,rij.v2,af.v2);
pf = _mm256_sub_ps(pf,mr_inv);
} // loop of j
af.v0 = _mm256_blend_ps(__fkg_tmp9.v0,af.v0,0b00000001);
af.v1 = _mm256_blend_ps(__fkg_tmp9.v1,af.v1,0b00000001);
af.v2 = _mm256_blend_ps(__fkg_tmp9.v2,af.v2,0b00000001);
id_maxf = _mm256_castps_si256(_mm256_blend_ps(_mm256_castsi256_ps(__fkg_tmp10),_mm256_castsi256_ps(id_maxf),0b00000001));
id_minf = _mm256_castps_si256(_mm256_blend_ps(_mm256_castsi256_ps(__fkg_tmp11),_mm256_castsi256_ps(id_minf),0b00000001));
ngbf = _mm256_castps_si256(_mm256_blend_ps(_mm256_castsi256_ps(__fkg_tmp12),_mm256_castsi256_ps(ngbf),0b00000001));
pf = _mm256_blend_ps(__fkg_tmp13,pf,0b00000001);
rankf = _mm256_castps_si256(_mm256_blend_ps(_mm256_castsi256_ps(__fkg_tmp14),_mm256_castsi256_ps(rankf),0b00000001));
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
id_maxf = _mm256_max_epi32(id_maxf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(id_maxf),_mm256_castsi256_ps(id_maxf),0xb1)));
id_maxf = _mm256_max_epi32(id_maxf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(id_maxf),_mm256_castsi256_ps(id_maxf),0xee)));
id_maxf = _mm256_max_epi32(id_maxf,_mm256_castps_si256(_mm256_castps128_ps256(_mm256_extractf128_ps(_mm256_castsi256_ps(id_maxf),1))));
((int*)&force[i+0].neighbor.id_max)[0] = max(((int*)&force[i+0].neighbor.id_max)[0],(PIKG::S32)id_maxf[0]);}

{
id_minf = _mm256_min_epi32(id_minf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(id_minf),_mm256_castsi256_ps(id_minf),0xb1)));
id_minf = _mm256_min_epi32(id_minf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(id_minf),_mm256_castsi256_ps(id_minf),0xee)));
id_minf = _mm256_min_epi32(id_minf,_mm256_castps_si256(_mm256_castps128_ps256(_mm256_extractf128_ps(_mm256_castsi256_ps(id_minf),1))));
((int*)&force[i+0].neighbor.id_min)[0] = min(((int*)&force[i+0].neighbor.id_min)[0],(PIKG::S32)id_minf[0]);}

{
ngbf = _mm256_add_epi32(ngbf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(ngbf),_mm256_castsi256_ps(ngbf),0xb1)));
ngbf = _mm256_add_epi32(ngbf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(ngbf),_mm256_castsi256_ps(ngbf),0xee)));
ngbf = _mm256_add_epi32(ngbf,_mm256_castps_si256(_mm256_castps128_ps256(_mm256_extractf128_ps(_mm256_castsi256_ps(ngbf),1))));
((int*)&force[i+0].neighbor.number)[0] = (((int*)&force[i+0].neighbor.number)[0]+ngbf[0]);
}

{
pf = _mm256_add_ps(pf,_mm256_shuffle_ps(pf,pf,0xb1));
pf = _mm256_add_ps(pf,_mm256_shuffle_ps(pf,pf,0xee));
pf = _mm256_add_ps(pf,_mm256_castps128_ps256(_mm256_extractf128_ps(pf,1)));
((float*)&force[i+0].phi)[0] = (((float*)&force[i+0].phi)[0]+pf[0]);
}

{
rankf = _mm256_add_epi32(rankf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(rankf),_mm256_castsi256_ps(rankf),0xb1)));
rankf = _mm256_add_epi32(rankf,_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(rankf),_mm256_castsi256_ps(rankf),0xee)));
rankf = _mm256_add_epi32(rankf,_mm256_castps_si256(_mm256_castps128_ps256(_mm256_extractf128_ps(_mm256_castsi256_ps(rankf),1))));
((int*)&force[i+0].neighbor.rank)[0] = (((int*)&force[i+0].neighbor.rank)[0]+rankf[0]);
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
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
