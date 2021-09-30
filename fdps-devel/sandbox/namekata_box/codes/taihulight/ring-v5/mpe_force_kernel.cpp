#include"mpe_func.h"
#include <cmath>

void calcForcePlanet(const EpiMM epi[], const int n_epi, const double mp, ForceMM force[]) {
  for (int i=0; i<n_epi; i++) {
    const double x=epi[i].x;
    const double y=epi[i].y;
    const double z=epi[i].z;
    const double r_sq = x*x+y*y+z*z;
    const double r_inv = 1.0 / sqrt(r_sq);
    const double pij = mp * r_inv;
    const double mri3 = r_inv*r_inv*pij;
    force[i].ax -= mri3 * x;
    force[i].ay -= mri3 * y;
    force[i].az -= mri3 * z;
    force[i].pot -= pij;
  }
}

void calcDirectGrav(const EpiMM epi[], const int n_epi, 
                    const EpjMM epj[], const int n_epj, const int adr_epj[],
                    ForceMM  force[], double r_epi, double r_epj, double kappa, double eta, const bool adr_flag){
  int j;
  const double r_coll = r_epi + r_epj;
  const double r_coll_inv_q3 = 1.0/(r_coll*r_coll*r_coll);
    for(int i=0; i<n_epi; i++){
      const int id_i = epi[i].id;
      //      force[i].ax = force[i].ay = force[i].az = force[i].pot = 0.0;
      //        force[i].n_coll = force[i].adr_ngb = 0;
      //        double r_min_sq = 9999999.9;
      const double mi = epi[i].mass;
      for(int jk=0; jk<n_epj; jk++){
        j = adr_flag?adr_epj[jk]:jk;
        if(id_i == epj[j].id) continue;
        const double dx = epi[i].x - epj[j].x;
        const double dy = epi[i].y - epj[j].y;
        const double dz = epi[i].z - epj[j].z;
        const double dvx = epi[i].vx - epj[j].vx;
        const double dvy = epi[i].vy - epj[j].vy;
        const double dvz = epi[i].vz - epj[j].vz;
        const double mj = epj[j].mass;
        const double r_sq = dx*dx + dy*dy + dz*dz;
        const double r_inv = 1.0 / sqrt(r_sq);
        const double r_inv_sq = r_inv*r_inv;
        const double xij  = 1.0 - r_coll*r_inv;
        const double rv = (dx*dvx + dy*dvy + dz*dvz)*r_inv_sq;
        const double pij = mj * r_inv;
        double mri3;
        if (xij<0) {
          const double uij = mi*mj/(mi+mj);
          mri3 = mj*r_coll_inv_q3 + kappa*uij/mi * xij + eta*uij/mi * rv;
          //          std::cerr<<"coll ["<<i<<"] - ["<<j<<"] r_inv="<<r_inv<<std::endl;
        }
        else mri3 = r_inv_sq*pij;
        force[i].ax -= mri3 * dx;
        force[i].ay -= mri3 * dy;
        force[i].az -= mri3 * dz;
        force[i].pot -= pij;
//           if(r_sq_real < r_coll_sq_mm){
//                force[i].n_coll++;
//                if(r_sq_real < r_min_sq){
//                    r_min_sq = r_sq;
//                    force[i].adr_ngb = epj[j].adr;
//                }
//            }
//            force[i].r_ngb_sq = r_min_sq;            
      }
    }
}

void calcSPGrav(const EpiMM epi[], const int n_epi,
                const SpjMM spj[], const int n_spj, const int adr_spj[], 
                ForceMM  force[], const bool adr_flag){
    int j;
    for(int i=0; i<n_epi; i++){
      //      force[i].ax = force[i].ay = force[i].az = force[i].pot = 0.0;
      const double mi = epi[i].mass;
      for(int jk=0; jk<n_spj; jk++){
        j = adr_flag?adr_spj[jk]:jk;
        const double dx = epi[i].x - spj[j].x;
        const double dy = epi[i].y - spj[j].y;
        const double dz = epi[i].z - spj[j].z;
        const double mj = spj[j].mass;
        const double r_sq = dx*dx + dy*dy + dz*dz;
        const double r_inv = 1.0 / sqrt(r_sq);
        const double r_inv_sq = r_inv*r_inv;
        const double pij = mj * r_inv;
        double mri3 = r_inv_sq*pij;
        force[i].ax -= mri3 * dx;
        force[i].ay -= mri3 * dy;
        force[i].az -= mri3 * dz;
        force[i].pot -= pij;
      }
    }
}
