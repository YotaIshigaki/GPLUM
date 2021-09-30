#include<math.h>
#include<stdio.h>
#include"kernel_sunway.h"

int     disp_lm[N_WALK_LIMIT+2][3];
EpiMM   epi_lm[NI_LIMIT];
EpjMM   epj_lm[NJ_LIMIT];
SpjMM   spj_lm[NJ_LIMIT];
ForceMM force_lm[NI_LIMIT];
double r_coll_sq_lm;

void CalcGrav(const int n_epi, const EpiMM * epi, 
              const int n_epj, const EpjMM * epj,
              const int n_spj, const SpjMM * spj,
              ForceMM * force, const double eps2){
    int i, j;
    double mass = 0.0;
    double mass_epj = 0.0;
    double mass_spj = 0.0;
    for(i=0; i<n_epi; i++){
        mass = 0.0;
        mass_epj = 0.0;
        mass_spj = 0.0;
        const double xi = epi[i].x;
        const double yi = epi[i].y;
        const double zi = epi[i].z;
        //double r_min_sq = r_coll_sq_lm;
        double r_min_sq = 9999999.9;
        const int id_i = epi[i].id;
        for(j=0; j<n_epj; j++){
            if(id_i == epj[j].id) continue;
            const double dx = xi - epj[j].x;
            const double dy = yi - epj[j].y;
            const double dz = zi - epj[j].z;
            const double mj = epj[j].mass;
            const double r_sq_real = dx*dx + dy*dy + dz*dz;
            const double r_sq = eps2 + r_sq_real;
            //const double r_sq = eps2 + dx*dx + dy*dy + dz*dz;
            const double r_inv = 1.0 / sqrt(r_sq);
            const double pij = mj * r_inv;
            const double mri3 = r_inv*r_inv*pij;
            force[i].ax -= mri3 * dx;
            force[i].ay -= mri3 * dy;
            force[i].az -= mri3 * dz;
            force[i].pot -= pij;
            mass_epj += mj;
            if(r_sq_real < r_coll_sq_lm){
                force[i].n_coll++;
                //printf("r_min_sq=%e, r_sq=%e, epj[j].id=%d \n", r_min_sq, r_sq, epj[j].id);
                if(r_sq_real < r_min_sq){
                    //printf("r_min_sq=%e, r_sq=%e, epj[j].id=%d \n", r_min_sq, r_sq, epj[j].id);
                    r_min_sq = r_sq;
                    force[i].adr_ngb = epj[j].adr;
                }
            }
            force[i].r_ngb_sq = r_min_sq;
        }
        for(j=0; j<n_spj; j++){
            const double dx = xi - spj[j].x;
            const double dy = yi - spj[j].y;
            const double dz = zi - spj[j].z;
            const double mj = spj[j].mass;
            const double r_sq = eps2 + dx*dx + dy*dy + dz*dz;
            const double r_inv = 1.0 / sqrt(r_sq);
            const double pij = mj * r_inv;
            const double mri3 = r_inv*r_inv*pij;
            force[i].ax -= mri3 * dx;
            force[i].ay -= mri3 * dy;
            force[i].az -= mri3 * dz;
            force[i].pot -= pij;
            mass_spj += mj;
        }
        mass = mass_epj + mass_spj;
        //printf("pot(0)= %8.5e \n", force[i].pot);
    }
    //printf("i =%d  n_epj= %d n_spj= %d m_epj= %e m_spj= %e m= %e \n", i, n_epj, n_spj, mass_epj, mass_spj, mass);
}

//double POT = 0.0;
void DispatchKernelSunWay(const int n_walk, 
                          const double eps2){
    int i, j;
    r_coll_sq_lm = r_coll_sq_mm;
    //printf("r_coll_sq_lm=%e \n", r_coll_sq_lm);
    for(i=0; i<n_walk; i++){
        int epi_head = disp_mm[i][0];
        int epj_head = disp_mm[i][1];
        int spj_head = disp_mm[i][2];
        int n_epi = disp_mm[i+1][0] - epi_head;
        int n_epj = disp_mm[i+1][1] - epj_head;
        int n_spj = disp_mm[i+1][2] - spj_head;
        for(j=0; j<n_epi; j++) epi_lm[epi_head+j] = epi_mm[epi_head+j];
        for(j=0; j<n_epj; j++) epj_lm[epj_head+j] = epj_mm[epj_head+j];
        for(j=0; j<n_spj; j++) spj_lm[spj_head+j] = spj_mm[spj_head+j];
        for(j=0; j<n_epi; j++){
            const int adr_j = epi_head+j;
            force_lm[adr_j].ax = force_lm[adr_j].ay = force_lm[adr_j].az = force_lm[adr_j].pot = 0.0;
            force_lm[adr_j].n_coll = 0;
        }
        CalcGrav(n_epi, epi_lm+epi_head, n_epj, epj_lm+epj_head, n_spj, spj_lm+spj_head, force_lm+epi_head, eps2);
    }
}


void RetrieveKernelSunWay(const int n_epi){
    int i;
    for(i=0; i<n_epi; i++){
        force_mm[i].ax  = force_lm[i].ax;
        force_mm[i].ay  = force_lm[i].ay;
        force_mm[i].az  = force_lm[i].az;
        force_mm[i].pot = force_lm[i].pot;
        force_mm[i].n_coll  = force_lm[i].n_coll;
        force_mm[i].adr_ngb = force_lm[i].adr_ngb;
        force_mm[i].r_ngb_sq = force_lm[i].r_ngb_sq;
    }
}

