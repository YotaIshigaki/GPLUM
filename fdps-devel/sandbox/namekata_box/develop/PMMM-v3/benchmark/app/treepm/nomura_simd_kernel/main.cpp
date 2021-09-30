//
// main.cpp
//
// nbody simulation code using P^3M method
//
// Kentaro NOMURA 2018/02/25
//
// Known problem:
// nothing so far

#include<iostream>
#include<particle_simulator.hpp>
#include<particle_mesh.hpp>
#include<param_fdps.h>

#include <cstdio>
#include <string>
#include <vector>

#include "user_defined_class.h"
#include "kernel.hpp"


template <class Tp>
void setup_random_particle(Tp & ptcl,
                           const PS::S32 npart_total,
                           const PS::F32 omegam = 0.27)
{
    PS::S32 rank = PS::Comm::getRank();
    PS::S64 npart_local = (rank == 0) ? npart_total : 0;

    ptcl.setNumberOfParticleLocal(npart_local);

    PS::MT::init_genrand(0);

    for (PS::S32 i=0;i<npart_local;i++) {
        ptcl[i].charge = 3.0*omegam/(8.0*M_PI*(PS::F32)npart_total);
        ptcl[i].pos[0] = PS::MT::genrand_res53();
        ptcl[i].pos[1] = PS::MT::genrand_res53();
        ptcl[i].pos[2] = PS::MT::genrand_res53();
        ptcl[i].vel[0] = 0.0;
        ptcl[i].vel[1] = 0.0;
        ptcl[i].vel[2] = 0.0;
    }
}

int main(int argc, char *argv[]){
  PS::Initialize(argc, argv);

  PS::S32 nmol = 8192;
  PS::F64 theta = 0.0;

  PS::S32 n_group_limit = 256;
  PS::S32 n_leaf_limit = 32;

  int c;
  while((c=getopt(argc,argv,"N:t:g:l:h")) != -1){
    switch(c){
    case 'N':
      nmol = atoi(optarg);
      if(PS::Comm::getRank()==0)std::cout<<"nmol="<<nmol<<std::endl;
      break;
    case 't':
      theta = atof(optarg);
      if(PS::Comm::getRank()==0)std::cout<<"theta="<<theta<<std::endl;
      break;
    case 'g':
      n_group_limit = atoi(optarg);
      if(PS::Comm::getRank()==0)std::cout<<"n_group_limit="<<n_group_limit<<std::endl;
      break;
    case 'l':
      n_leaf_limit = atoi(optarg);
      if(PS::Comm::getRank()==0)std::cout<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
      break;
    case 'h':
      if(PS::Comm::getRank()==0){
        std::cout<<"N: number of particle (default: 1024)"<<std::endl;
        std::cout<<"t: theta for tree calculation (default: 0.0)"<<std::endl;
        std::cout<<"g: n_group_limit (default: 64)"<<std::endl;
        std::cout<<"l: n_leaf_limit (default: 32)"<<std::endl;
      }
    default:
      PS::Abort();
    }
  }

  PS::ParticleSystem<FP>   system_nbody;
  system_nbody.initialize();
  const PS::F64 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
                         PS::F64vec(1.0,1.0,1.0));
  setup_random_particle(system_nbody,nmol);
  system_nbody.adjustPositionIntoRootDomain(dinfo);

  dinfo.decomposeDomainAll(system_nbody);
  system_nbody.exchangeParticle(dinfo);

  CalcForceAll p3m_tree;
  p3m_tree.initialize(system_nbody.getNumberOfParticleLocal(),theta,n_leaf_limit,n_group_limit);

  const PS::F64 rcut = RCUT;
  PS::F32vec *wAVX,*woAVX;
  wAVX  = new PS::F32vec[system_nbody.getNumberOfParticleLocal()];
  woAVX = new PS::F32vec[system_nbody.getNumberOfParticleLocal()];
  p3m_tree.clearTimeProfile();
  p3m_tree.calcForceAllAndWriteBack(PP::CalcForceEpEp<true>(rcut),
                                    PP::CalcForceEpSp<true>(rcut),
                                    system_nbody,dinfo);
  PS::TimeProfile et_wAVX = p3m_tree.getTimeProfile();
  for(int i=0;i<system_nbody.getNumberOfParticleLocal();i++) wAVX[i] = system_nbody[i].acc;
  p3m_tree.clearTimeProfile();
  p3m_tree.calcForceAllAndWriteBack(PP::CalcForceEpEp<false>(rcut),
                                    PP::CalcForceEpSp<false>(rcut),
                                    system_nbody,dinfo);
  PS::TimeProfile et_woAVX = p3m_tree.getTimeProfile();
  for(int i=0;i<system_nbody.getNumberOfParticleLocal();i++) woAVX[i] = system_nbody[i].acc;

  if (PS::Comm::getRank() == 0) {
    std::cout << "time (w/ AVX)  = " << et_wAVX.calc_force_only << " [s]" << std::endl;
    std::cout << "time (w/o AVX) = " << et_woAVX.calc_force_only << " [s]" << std::endl;
  }

  PS::F64vec error = 0.0;
  for(int i=0;i<system_nbody.getNumberOfParticleLocal();i++){
    PS::F64vec diff;
    diff.x = std::abs(wAVX[i].x - woAVX[i].x);
    diff.y = std::abs(wAVX[i].y - woAVX[i].y);
    diff.z = std::abs(wAVX[i].z - woAVX[i].z);
    error.x = diff.x > error.x ? diff.x : error.x;
    error.y = diff.y > error.y ? diff.y : error.y;
    error.z = diff.z > error.z ? diff.z : error.z;
  }
  error = PS::Comm::getMaxValue(error);
  if(PS::Comm::getRank() == 0) printf("max_error: %e %e %e\n",error.x,error.y,error.z);

  delete[] wAVX;
  delete[] woAVX;

  PS::Finalize();
  return 0;
}
