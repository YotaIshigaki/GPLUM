#include "user_defined_class.h"
struct CalcForceLongEpEp{
CalcForceLongEpEp(){}
void operator()(const EPI* epi,const int ni,const EPJ* epj,const int nj,Force *force){
PS::F32vec rij;
PS::F32 r2;
PS::F32 rinv;
PS::F32 mrinv;
PS::F32 __fkg_tmp0;
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
for(int i=0;i<ni;i++){
PS::F32 eps2i = epi[i].eps2;
PS::F32vec xiloc = xiloc_tmp[i];
PS::F32vec f = force[i].acc;
PS::F32 phi = force[i].pot;
for(int j=0;j<nj;j++){
PS::F32 mj = epj[j].mass;
PS::F32 eps2j = epj[j].eps2;
PS::F32vec xjloc = xjloc_tmp[j];
rij.x = (xiloc.x-xjloc.x);
rij.y = (xiloc.y-xjloc.y);
rij.z = (xiloc.z-xjloc.z);
r2 = ((((eps2i+eps2j) + rij.x*rij.x) + rij.y*rij.y) + rij.z*rij.z);
#pragma statement fission_point
rinv = rsqrt<PS::F32,PS::F32>(r2);
#pragma statement fission_point
mrinv = (mj*rinv);
__fkg_tmp0 = ((mrinv*rinv)*rinv);
f.x = (f.x - __fkg_tmp0*rij.x);
f.y = (f.y - __fkg_tmp0*rij.y);
f.z = (f.z - __fkg_tmp0*rij.z);
phi = (phi-mrinv);
}
force[i].acc = f;
force[i].pot = phi;
}
}
template<typename Tret,typename Top>
Tret rsqrt(Top op){ return (Tret)1.0/std::sqrt(op); }
template<typename Tret,typename Top>
Tret sqrt(Top op){ return std::sqrt(op); }
template<typename Tret,typename Top>
Tret inv(Top op){ return 1.0/op; }
template<typename Tret,typename Ta,typename Tb,typename Tc>
Tret madd(Ta a,Tb b,Tc c){ return a*b+c; }
template<typename Tret,typename Ta,typename Tb,typename Tc>
Tret msub(Ta a,Tb b,Tc c){ return c - a*b; }
template<typename Tret,typename Ta,typename Tb,typename Tc>
Tret nmadd(Ta a,Tb b,Tc c){ return -(a*b+c); }
template<typename Tret,typename Ta,typename Tb,typename Tc>
Tret nmsub(Ta a,Tb b,Tc c){ return a*b-c; }
template<typename Tret,typename Ta,typename Tb>
Tret max(Ta a,Tb b){ return std::max(a,b);}template<typename Tret,typename Ta,typename Tb>
Tret min(Ta a,Tb b){ return std::min(a,b);}};
