/* ideal generated hpp file */
#pragma once

#include<particle_simulator.hpp>

struct Force{
  PS::F64vec f;
  PS::F64 u;
  void clear(){
    f = 0.0;
    u = 0.0;
  }
};

struct FP{
  PS::F64vec r;
  PS::F64vec v;
  PS::F64vec f;
  PS::F64    u;
  PS::F64 sigma;
  PS::F64 epsilon;
  void copyFromForce(const Force& _f){
    f = _f.f;
    u = _f.u;
  }
  PS::F64vec getPos() const {
    return r;
  }
  void setPos(const PS::F64vec& _r){
    r = _r;
  }
};

struct EPI{
  PS::F64vec r;
  PS::F64 sigma;
  PS::F64 epsilon;

  void copyFromFP(const FP& fp){
    r       = fp.r;
    sigma   = fp.sigma;
    epsilon = fp.epsilon
  }
  PS::F64vec getPos(){return r;}
  void setPos(const PS::F64vec r){
    this->r = r;
  }
  PS::F64 getRSearch(){
    return 4.0;
  }
};

struct EPJ{
  PS::F64vec r;
  PS::F64 sigma;
  PS::F64 epsilon;
  void copyFromFP(const FP& fp){
    r       = fp.r;
    sigma   = fp.sigma;
    epsilon = fp.epsilon
  }
  PS::F64vec getPos(){return r;}
  void setPos(const PS::F64vec r){
    this->r = r;
  }
  PS::F64 getRSearch(){
    return 4.0;
  }
};
/*
PS::F64 eps = 1e-32;

.
.
.

calcForceAllAnd...(DispatchKernel(eps)),RetrieveKernel(),...);
*/

template <typename Tpsys,typename Ttree>
void calcForceLJ(Tpsys& psys,Ttree& tree,..., F32 eps){
  #ifdef ENABLE_GPU
  tree.calcForceAllAndWriteBackMultiWalkIndex
    (DispatchKernelLJ(eps),
     RetrieveKernelLJ(),
     ...);
  #else
  tree.calcForceAllAndWriteBack
    (KernelLJ(eps), n_group_limit, ..., psys, dinfo);
  tree.calcForceAllAndWriteBack
    (KernelLJ(eps),psys,dinfo,reuse);
  #endif
}
template <typename Tpsys,typename Ttree>
void calcForceSPH(Tpsys& psys,Ttree& tree){
  #ifdef GPU
  tree.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernel(eps),
					 RetrieveKernel(),
					 ...);
  #else
  tree.calcForceAllAndWriteBack(Kernel(eps), ...);
  #endif
}
 

struct Kernel{
  PS::F32 eps;
  Kernel(eps):eps(eps){};
  void operator()(const EPI* epi,
		  const int  nepi,
		  const EPJ* epj,
		  const int  nepj,
		  Force*     force){
    for(int i=0;i<nepi;i++){
      const PS::F64vec r_i       = epi[i].r;
      const PS::F64    sigma_i   = epi[i].sigma;
      const PS::F64    epsilon_i = epi[i].epsilon;
      for(int j=0;j<nepj;j++){
	const PS::F64vec r_j       = epj[j].r;
	const PS::F64    sigma_j   = epj[j].sigma;
	const PS::F64    epsilon_j = epj[j].epsilon;

	const PS::F64vec dr = r_i - r_j;
	const PS::F64    r2 = dr*dr;
	if(r2 != 0.0 && r2 < 16.0){
	  const PS::F64 epsilon_ij = sqrt(epsilon_i * epsilon_j);
	  const PS::F64 sigma_ij = 0.5*(sigma_i + sigma_j);
	  const PS::F64 tmp_ij = sigma_ij*sigma_ij*sigma_ij*sigma_ij*sigma_ij*sigma_ij / (r2*r2*r2);
	  force[i].f += (24*epsilon_ij*tmp_ij*(2.0*tmp_ij - 1.0))*dr;
	  force[i].u += 4*epsilon_ij*tmp_ij*(tmp_ij - 1.0);
	}
      }
    }
  }
}
