#include <math.h>

#include "user_defined_class.h"
#include "constant.h"

void CalcForceLongEpEpOrig(const EPI * ep_i,
			   const int n_ip,
			   const EPJ * ep_j,
			   const int n_jp,
			   Force * force){
  PS::F32vec xiloc[n_ip],xjloc[n_jp];
  for(int i=0;i<n_ip;i++) xiloc[i] = ep_i[i].pos - ep_i[0].pos;
  for(int j=0;j<n_jp;j++) xjloc[j] = ep_j[j].pos - ep_i[0].pos;

  for(int i=0; i<n_ip; i++){
    const PS::F32vec pos_i = xiloc[i];
    const PS::F32    eps2i = ep_i[i].eps2;
    PS::F32vec acc_i = 0.f;
    PS::F32    pot_i = 0.f;
    for(int j=0; j<n_jp; j++){
      PS::F32vec rij = xjloc[j] - pos_i;
      const float r2 = rij*rij + eps2i + ep_j[j].eps2;
      const float r_inv = 1.f/sqrt(r2);
      const float m_r = ep_j[j].mass * r_inv;
      const float m_r3 = m_r * r_inv * r_inv;
      acc_i += m_r3 * rij;
      pot_i -= m_r;
    }
    force[i].acc += acc_i;
    force[i].pot += pot_i;
  }
};

void CalcForceLongEpSpOrig(const EPI * ep_i,
			   const int n_ip,
			   const SPJ * sp_j,
			   const int n_jp,
			   Force * force){
  for(int i=0; i<n_ip; i++){
    PS::F64vec ri = ep_i[i].pos;
    PS::F64    ei = ep_i[i].eps2;
    PS::F64vec fi;
    PS::F64    pi;
    fi = 0.0;
    pi = 0.0;
    for(int j=0; j<n_jp; j++){
      PS::F64vec rj = sp_j[j].pos;
      PS::F64 qxx,qyy,qzz,qxy,qyz,qzx,mj,ej,tr;
      mj  = sp_j[j].mass;
      ej  = sp_j[j].eps2;      
      qxx = sp_j[j].xx;
      qyy = sp_j[j].yy;
      qzz = sp_j[j].zz;
      qxy = sp_j[j].xy;
      qyz = sp_j[j].yz;
      qzx = sp_j[j].zx;
      tr  = qxx + qyy + qzz;
      const PS::F64 dx = rj.x - ri.x;
      const PS::F64 dy = rj.y - ri.y;
      const PS::F64 dz = rj.z - ri.z; // 3
      const PS::F64 r2 = ei + ej + dx*dx + (dy*dy + (dz*dz + eps2)); //6
      const PS::F64 qrx = qxx*dx + qxy*dy + qzx*dz;
      const PS::F64 qry = qyy*dy + qyz*dz + qxy*dx;
      const PS::F64 qrz = qzz*dz + qzx*dx + qyz*dy; // 15
      PS::F64 qrr    = qrx * dx + qry * dy + qrz * dz; // 5
      PS::F64 r_inv  = 1./sqrt(r2); // 4
      PS::F64 r2_inv = r_inv * r_inv;
      PS::F64 r3_inv = r2_inv * r_inv;
      PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
      PS::F64 qrr_r5 = r5_inv * qrr;
      PS::F64 qrr_r7 = r2_inv * qrr_r5;
      PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
      PS::F64 B = -2.0*r5_inv; //12
      fi.x -= A*dx + B*qrx;
      fi.y -= A*dy + B*qry;
      fi.z -= A*dz + B*qrz; // 12
      pi -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
      // total 57 flop
    }
    force[i].acc.x += fi.x;
    force[i].acc.y += fi.y;
    force[i].acc.z += fi.z;
    force[i].pot   += pi;
  }
}
