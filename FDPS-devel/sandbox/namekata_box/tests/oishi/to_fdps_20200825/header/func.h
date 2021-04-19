#pragma once

#include "Eigen/Core"
#include "Eigen/Eigen"

void makeOutputDirectory(char *dir_name){
  struct stat st;
  PS::S32 ret;
  if(PS::Comm::getRank()==0){
    if(stat(dir_name,&st)!=0){
      ret=mkdir(dir_name,0777);
    }else{
      ret =0;//the directory named dir_name already exists.
    }
  }
  PS::Comm::broadcast(&ret,1);
  if(ret==0){
    if(PS::Comm::getRank()==0)
      fprintf(stderr,"Directory \"%s\" is successfully made.\n",dir_name);
    }else{
      if(PS::Comm::getRank()==0)
      fprintf(stderr,"Directory %s fails to be made./n",dir_name);
      PS::Abort();
    }
}

void InitialKick(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
  const PS::F64 dt_half = 0.5 * dt;
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++i){
    if( sph_system[i].itype > S_W ) continue;
    // initial kick
    sph_system[i].dens_half = sph_system[i].dens + dt_half * sph_system[i].dens_rate;
    sph_system[i].vel_half  = sph_system[i].vel + dt_half * sph_system[i].acc;
    sph_system[i].sig_half  = sph_system[i].sig + product_02_3d( dt_half, sph_system[i].sig_rate);
    // full drift
    sph_system[i].pos      += dt * sph_system[i].vel_half;
    // predict
    sph_system[i].dens     += dt * sph_system[i].dens_rate;
    sph_system[i].vel      += dt * sph_system[i].acc;
    sph_system[i].sig      += product_02_3d( dt, sph_system[i].sig_rate);
    // sph_system[i].ReturnMapping(dt);
  }
}

void EigenSolver(PS::ParticleSystem<FP> &sph_system){
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> ES;
  Eigen::Matrix3d sigma, Rp, N, R;
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++i){
    const PS::F64 densi_sq = sph_system[i].dens * sph_system[i].dens;
    Rp = Eigen::Matrix3d::Zero();
    sigma <<  sph_system[i].sig.xx, sph_system[i].sig.xy, sph_system[i].sig.xz,
              sph_system[i].sig.xy, sph_system[i].sig.yy, sph_system[i].sig.yz,
              sph_system[i].sig.xz, sph_system[i].sig.yz, sph_system[i].sig.zz;
    ES.compute(sigma);
    PS::F64 ps_x = ES.eigenvalues()(0), ps_y = ES.eigenvalues()(1), ps_z = ES.eigenvalues()(2);
    N = ES.eigenvectors();
    if (ps_x > 0.0) Rp(0,0) = - eps_AS * ps_x / densi_sq;
    if (ps_y > 0.0) Rp(1,1) = - eps_AS * ps_y / densi_sq;
    if (ps_z > 0.0) Rp(2,2) = - eps_AS * ps_z / densi_sq;
    R = N * Rp * N.transpose();
    // sph_system[i].R_tensor( R(0,0), R(1,1), R(2,2), R(0,1), R(0,2), R(1,2) );
    sph_system[i].R_tensor.xx = R(0,0);
    sph_system[i].R_tensor.yy = R(1,1);
    sph_system[i].R_tensor.zz = R(2,2);
    sph_system[i].R_tensor.xy = R(0,1);
    sph_system[i].R_tensor.xz = R(0,2);
    sph_system[i].R_tensor.yz = R(1,2);
  }
}


void FinalKick(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
  const PS::F64 dt_half = 0.5 * dt;
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++i){
    if( sph_system[i].itype > S_W ) continue;
    // calculate stress rate (Drucker-Prager)
    sph_system[i].CalcStressRateDP_3d(dt);
    // final kick
    sph_system[i].dens = sph_system[i].dens_half + dt_half * sph_system[i].dens_rate;
    sph_system[i].vel  = sph_system[i].vel_half + dt_half * sph_system[i].acc;
    sph_system[i].eps += product_02_3d( dt, sph_system[i].eps_rate);
    sph_system[i].sig  = sph_system[i].sig_half + product_02_3d( dt_half, sph_system[i].sig_rate);
    sph_system[i].ReturnMapping(dt);
    // sph_system[i].pres = sph_system[i].eps.xy;
  }
}

void GetPRank(PS::ParticleSystem<FP>& sph_system){
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++i){
    sph_system[i].prank = PS::Comm::getRank();
  }
}

void CalcSum(const PS::ParticleSystem<FP>& sph_system, const PS::F64 time){
  // ene zero clear
  PS::F64 Hene  = zero;//potential energy
  PS::F64 Kene  = zero;//kinetic energy
  PS::F64 HKene = zero;//potential + kinetic energy
  PS::F64 Tene  = zero;//mechanical energy TOTAL!
  PS::F64 Eene  = zero;//elastic energy
  PS::F64 Pene  = zero;//plastic energy
  PS::F64 I1    = zero;//total I1
  PS::F64 XX    = zero;
  PS::F64 YY    = zero;
  PS::F64 XY    = zero;
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++i){
    if (sph_system[i].itype > S_W) continue;// calc energy only soil particle
    Hene += sph_system[i].mass * gravity * sph_system[i].pos.y;
    // PS::F64 dv2 = sph_system[i].vel * sph_system[i].vel;
    Kene += 0.5 * sph_system[i].mass * (sph_system[i].vel * sph_system[i].vel);
    Eene += dx2 * 0.5 * (sph_system[i].sig.xx * sph_system[i].epse.xx + sph_system[i].sig.yy * sph_system[i].epse.yy
    + sph_system[i].sig.zz * sph_system[i].epse.zz + 2.0 * (sph_system[i].sig.xy * sph_system[i].epse.xy 
    + sph_system[i].sig.yz * sph_system[i].epse.yz + sph_system[i].sig.xz * sph_system[i].epse.xz) );
    Pene += dx2 * (sph_system[i].sig.xx * sph_system[i].epsp.xx + sph_system[i].sig.yy * sph_system[i].epsp.yy
    + sph_system[i].sig.zz * sph_system[i].epsp.zz + 2.0 * (sph_system[i].sig.xy * sph_system[i].epsp.xy
    + sph_system[i].sig.yz * sph_system[i].epsp.yz + sph_system[i].sig.xz * sph_system[i].epsp.xz) );
    // XX += sqrt(sph_system[i].sig.xx*sph_system[i].sig.xx) * sqrt(sph_system[i].epsp.xx*sph_system[i].epsp.xx);
    // YY += sqrt(sph_system[i].sig.yy*sph_system[i].sig.yy) * sqrt(sph_system[i].epsp.yy*sph_system[i].epsp.yy);
    // XY += sqrt(sph_system[i].sig.xy*sph_system[i].sig.xy) * sqrt(sph_system[i].epsp.xy*sph_system[i].epsp.xy);
    // I1 += sph_system[i].sig.xx + sph_system[i].sig.yy;
  }
  // HKene = Hene + Kene;
  // Tene = HKene + Eene + Pene;
  Hene  = PS::Comm::getSum(Hene);
  Kene  = PS::Comm::getSum(Kene);
  Eene  = PS::Comm::getSum(Eene);
  Pene  = PS::Comm::getSum(Pene);
  HKene = Hene + Kene;
  Tene  = HKene + Eene + Pene;
  // I1 = PS::Comm::getSum(I1);
  // HKene = PS::Comm::getSum(HKene);
  // Tene = PS::Comm::getSum(Tene);
  // XX = PS::Comm::getSum(XX);
  // YY = PS::Comm::getSum(YY);
  // XY = PS::Comm::getSum(XY);
  if (PS::Comm::getRank() == 0){
    //std::cout<<"potential = "<<Hene<<" ; kinetic = "<<Kene<<" ; total = "<<Tene<<std::endl;
    ene <<time <<"," <<Hene <<"," <<Kene << "," << Eene <<"," <<Pene << "," <<Tene << "," <<HKene <<"\n";
  }
}
  
