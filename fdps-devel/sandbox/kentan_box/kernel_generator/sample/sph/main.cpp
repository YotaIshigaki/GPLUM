#include <iostream>
#include <fstream>

#include <string>
#include <cassert>

#include <particle_simulator.hpp>

#include "hydro_force.h"

class EP_hydro{
public:
  PS::F64vec pos;
  PS::F32vec vel;
  PS::F32    mass;
  PS::F32    eng;
  PS::F32    h;
  PS::F32    dens;
  PS::F32    pres;
  PS::F32    gradh;
  PS::F32    snds;
  PS::F32    BalSW;
  PS::F32    alpha;
  PS::S32    dummy;

  friend std::istream& operator>>(std::istream& is,EP_hydro& ep){
    is >> ep.dummy
       >> ep.pos.x
       >> ep.pos.y
       >> ep.pos.z
       >> ep.vel.x
       >> ep.vel.y
       >> ep.vel.z
       >> ep.mass
       >> ep.eng
       >> ep.h
       >> ep.dens
       >> ep.pres
       >> ep.gradh
       >> ep.snds
       >> ep.BalSW
       >> ep.alpha;
    return is;
  }
  friend std::ostream& operator<<(std::ostream& os,EP_hydro& ep){
    os << ep.dummy  << " "
       << ep.pos.x  << " "
       << ep.pos.y  << " "
       << ep.pos.z  << " "
       << ep.vel.x  << " "
       << ep.vel.y  << " "
       << ep.vel.z  << " "
       << ep.mass   << " "
       << ep.eng    << " "
       << ep.h      << " "
       << ep.dens   << " "
       << ep.pres   << " "
       << ep.gradh  << " "
       << ep.snds   << " "
       << ep.BalSW  << " "
       << ep.alpha; 
    return os;
  }
};

class Force_hydro{
public:
  PS::F32vec acc;
  PS::F32 eng_dot;
  PS::F32 dt;
  void clear(){acc=0.f;eng_dot = dt = 0.f;}
  friend std::ostream& operator<<(std::ostream& os,const Force_hydro f){
    os << f.acc << " "
       << f.eng_dot << " "
       << f.dt;
    return os;
  }
};

//#include "kernel.h"
#define SQ(x)   ((x)*(x))
#define PWR5(x) ((x)*(x)*(x)*(x)*(x))

PS::F32vec gradW(const PS::F32vec dr, const PS::F32 h) {
    // Spatial gradiaent of Wendland C4 kernel
    const PS::F32 r = std::sqrt(dr * dr);
    const PS::F32 u = r/h;
    const PS::F32 p1u = std::max(0.0, 1.0 - u);
    const PS::F32 coeff = 1155.0 / (4.0 * M_PI * PWR5(h));
    return - coeff * PWR5(p1u) * (1.0 + 5.0 * u) * dr;

}

class CalcHydroForce{
public:
  PS::F32 gamma;
  PS::F32 CFL;
  CalcHydroForce(PS::F32 gamma,PS::F32 CFL):gamma(gamma),CFL(CFL){}
  void operator () (const EP_hydro * ep_i,
		    const PS::S32 n_ip,
		    const EP_hydro * ep_j,
		    const PS::S32 n_jp,
		    Force_hydro * force) {
    PS::F32vec xi[n_ip];
    for(int i=0;i<n_ip;i++){
      xi[i] = ep_i[i].pos - ep_i[0].pos;
    }
    PS::F32vec xj[n_jp];
    for(int j=0;j<n_jp;j++){
      xj[j] = ep_j[j].pos - ep_i[0].pos;
    }
    for (PS::S32 i = 0; i < n_ip; i++){
      const PS::F32vec pos_i = xi[i];
      const PS::F32vec vel_i = ep_i[i].vel;
      const PS::F32 mi = ep_i[i].mass;
      const PS::F32 hi = ep_i[i].h;
      const PS::F32 rhoi = ep_i[i].dens;
      const PS::F32 ui = ep_i[i].eng;
      const PS::F32 Pi = ep_i[i].pres;
      const PS::F32 fi = ep_i[i].gradh;
      const PS::F32 ci = ep_i[i].snds;
      const PS::F32 BalSWi = ep_i[i].BalSW;
      const PS::F32 ai = ep_i[i].alpha;
      PS::F32 v_sig_max = 0.0;
      for (PS::S32 j = 0; j < n_jp; j++){
	const PS::F32vec dr = pos_i - xj[j];
	const PS::F32vec dv = vel_i - ep_j[j].vel;
	const PS::F32 w_ij = (dv * dr < 0) ? dv * dr / std::sqrt(dr * dr) : 0;
	const PS::F32 v_sig = ci + ep_j[j].snds - 3.0 * w_ij;
	v_sig_max = std::max(v_sig_max, v_sig);
	const PS::F32 alpha_ij = 0.5 * (ai + ep_j[j].alpha);
	const PS::F32 AV = - 0.5 * alpha_ij * v_sig * w_ij
	  / (0.5 * (rhoi + ep_j[j].dens))
	  * 0.5 * (BalSWi + ep_j[j].BalSW);
	const PS::F32vec gradW_i  = gradW(dr, hi);
	const PS::F32vec gradW_j  = gradW(dr, ep_j[j].h);
	const PS::F32vec gradW_ij = 0.5 * (gradW_i + gradW_j);
	const PS::F32 mj = ep_j[j].mass;
	const PS::F32 uj = ep_j[j].eng;
	const PS::F32 Pj = ep_j[j].pres;
	const PS::F32 fj = ep_j[j].gradh;
	const PS::F32 fij = 1.0 - fi/((gamma - 1.0) * mj * uj);
	const PS::F32 fji = 1.0 - fj/((gamma - 1.0) * mi * ui);
	force[i].acc -= SQ(gamma - 1.0) * mj * ui * uj
	  * ( fij * gradW_i / Pi
	      + fji * gradW_j / Pj)
	  + mj * AV * gradW_ij;
	force[i].eng_dot += SQ(gamma - 1.0) * mj * ui * uj * fij * gradW_i * dv / Pi
	  + mj * 0.5 * AV * gradW_ij * dv;
      }
      force[i].dt = CFL * 2.0 * ep_i[i].h / v_sig_max;
    }
  }
};


int main(){
  int nepi,nepj;
  const std::string epi_in = "epi.in";
  const std::string epj_in = "epj.in";
  std::ifstream ifs_epi(epi_in); assert(!ifs_epi.fail());
  ifs_epi >> nepi;
  EP_hydro *epi = (EP_hydro*)malloc(nepi*sizeof(EP_hydro));
  for(int i=0;i<nepi;i++) ifs_epi >> epi[i];

  std::ifstream ifs_epj(epj_in); assert(!ifs_epj.fail());
  ifs_epj >> nepj;
  EP_hydro *epj = (EP_hydro*)malloc(nepj*sizeof(EP_hydro));
  for(int i=0;i<nepj;i++) ifs_epj >> epj[i];
#if 0
  std::cout << nepi << std::endl;
  for(int i=0;i<nepi;i++) std::cout << epi[i] << std::endl;
  std::cout << nepj << std::endl;
  for(int i=0;i<nepj;i++) std::cout << epj[i] << std::endl;
#endif

  Force_hydro *ref = (Force_hydro*)malloc(nepi*sizeof(Force_hydro));
  for(int i=0;i<nepi;i++) ref[i].clear();
  Force_hydro *tmp = (Force_hydro*)malloc(nepi*sizeof(Force_hydro));
  for(int i=0;i<nepi;i++) tmp[i].clear();

  const PS::F32 gamma = 5.f/3.f;
  const PS::F32 CFL = 0.3f;
  CalcHydroForce(gamma,CFL)(epi,nepi,epj,nepj,ref);
  Kernel<EP_hydro,EP_hydro,Force_hydro> kernel(gamma,CFL);
  kernel(epi,nepi,epj,nepj,tmp);
  printf("size EP:    %d\n",sizeof(EP_hydro));
  printf("size FORCE: %d\n",sizeof(Force_hydro));
  for(int i=0;i<nepi;i++){
    tmp[i].dt = CFL * 2.0 * epi[i].h / tmp[i].dt;
    tmp[i].acc     -= ref[i].acc;
    tmp[i].eng_dot -= ref[i].eng_dot;
    tmp[i].dt      -= ref[i].dt;
    std::cout << ref[i] << " " << tmp[i] << std::endl;
  }
  return 0;
}
