#include "user_defined.h"
#include "parameter.h"

#if 1
void calc_BS_ep_ep(Full_particle *ep_i,
                   int n_ip,
                   Full_particle *ep_j,
                   int n_jp,
                   Full_particle *f)
{
    int i, j;

    for (i=0; i<n_ip; i++) {
        double tmpUx = 0.0;
        double tmpUy = 0.0;
        Full_particle *pi = ep_i + i;
        const double posXi = pi->pos.x;
        const double posYi = pi->pos.y;
        for (j = 0; j < n_jp; j++) {
            Full_particle *pj = ep_j + j;
            // real
            const double posXj = pj->pos.x;
            const double posYj = pj->pos.y;
            const double rj2 = posXj * posXj + posYj * posYj;
            const double posXij = posXi - posXj;
            const double posYij = posYi - posYj;
            const double rij2 = posXij * posXij + posYij * posYij;
            double invrij2;
            if (rij2 > 0.0) invrij2 = 1.0/rij2;
            else invrij2 = 0.0;
            const double wzj = pj->wz;
            // pseudo
            const double pposXj = posXj * R0 * R0 / rj2;
            const double pposYj = posYj * R0 * R0 / rj2;
            const double pposXij = posXi - pposXj;
            const double pposYij = posYi - pposYj;
            const double prij2 = pposXij * pposXij + pposYij * pposYij;
            const double invprij2 = 1.0/prij2;
            tmpUx -= wzj * posYij * invrij2  - wzj * pposYij * invprij2;
            tmpUy += wzj * posXij * invrij2  - wzj * pposXij * invprij2;
        }
        Full_particle *pfi = f + i;
        pfi->vel.x += tmpUx;
        pfi->vel.y += tmpUy;
    }
}

void calc_BS_ep_sp(Full_particle *ep_i,
                   int n_ip,
                   fdps_spj_monopole *ep_j,
                   int n_jp,
                   Full_particle *f)
{
    int i, j;

    for (i=0; i<n_ip; i++) {
        double tmpUx = 0.0;
        double tmpUy = 0.0;
        Full_particle *pi = ep_i + i;
        const double posXi = pi->pos.x;
        const double posYi = pi->pos.y;
        for (j = 0; j < n_jp; j++) {
            fdps_spj_monopole *pj = ep_j + j;
            // real
            const double posXj = pj->pos.x;
            const double posYj = pj->pos.y;
            const double rj2 = posXj * posXj + posYj * posYj;
            const double posXij = posXi - posXj;
            const double posYij = posYi - posYj;
            const double rij2 = posXij * posXij + posYij * posYij;
            double invrij2;
            if (rij2 > 0.0) invrij2 = 1.0 / rij2;
            else invrij2 = 0.0;
            const double wzj = pj->mass;
            // pseudo
            const double pposXj = posXj * R0 * R0 / rj2;
            const double pposYj = posYj * R0 * R0 / rj2;
            const double pposXij = posXi - pposXj;
            const double pposYij = posYi - pposYj;
            const double prij2 = pposXij * pposXij + pposYij * pposYij;
            const double pinvrij2 = 1.0 / prij2;
            tmpUx -= wzj * posYij * invrij2 - wzj * pposYij * pinvrij2;
            tmpUy += wzj * posXij * invrij2 - wzj * pposXij * pinvrij2;
        }
        Full_particle *pfi = f + i;
        pfi->vel.x += tmpUx;
        pfi->vel.y += tmpUy;
    }
}
#else
void calc_BS_ep_ep(Full_particle *ep_i,
                        int n_ip,
                        Full_particle *ep_j,
                        int n_jp,
                        Full_particle *f)
{
    int i, j;

    for (i=0; i<n_ip; i++) {
        double tmpUx = 0.0;
        double tmpUy = 0.0;
        Full_particle *pi = ep_i + i;
        const double posXi = pi->pos.x;
        const double posYi = pi->pos.y;
        // j < i
        for (j = 0; j < i; j++) {
            Full_particle *pj = ep_j + j;
            // real
            const double posXj = pj->pos.x;
            const double posYj = pj->pos.y;
            const double rj2 = posXj * posXj + posYj * posYj;
            const double posXij = posXi - posXj;
            const double posYij = posYi - posYj;
            const double rij2 = posXij * posXij + posYij * posYij;
            const double wzj = pj->wz;
            // pseudo
            const double pposXj = posXj * R0 * R0 / rj2;
            const double pposYj = posYj * R0 * R0 / rj2;
            const double pposXij = posXi - pposXj;
            const double pposYij = posYi - pposYj;
            const double prij2 = pposXij * pposXij + pposYij * pposYij;
            tmpUx -= wzj * posYij / rij2 - wzj * pposYij / prij2;
            tmpUy += wzj * posXij / rij2 - wzj * pposXij / prij2;
        }
        // i=j pseudo only
        Full_particle *pj = ep_j + i;
        const double posXj = pj->pos.x;
        const double posYj = pj->pos.y;
        const double rj2 = posXj * posXj + posYj * posYj;
        const double wzj = pj->wz;
        const double pposXj = posXj * R0 * R0 / rj2;
        const double pposYj = posYj * R0 * R0 / rj2;
        const double pposXij = posXi - pposXj;
        const double pposYij = posYi - pposYj;
        const double prij2 = pposXij * pposXij + pposYij * pposYij;
        tmpUx -= - wzj * pposYij / prij2;
        tmpUy += - wzj * pposXij / prij2;
        // j > i
        for (j = i + 1; j < n_jp; j++) {
            Full_particle *pj = ep_j + j;
            // real
            const double posXj = pj->pos.x;
            const double posYj = pj->pos.y;
            const double rj2 = posXj * posXj + posYj * posYj;
            const double posXij = posXi - posXj;
            const double posYij = posYi - posYj;
            const double rij2 = posXij * posXij + posYij * posYij;
            const double wzj = pj->wz;
            // pseudo
            const double pposXj = posXj * R0 * R0 / rj2;
            const double pposYj = posYj * R0 * R0 / rj2;
            const double pposXij = posXi - pposXj;
            const double pposYij = posYi - pposYj;
            const double prij2 = pposXij * pposXij + pposYij * pposYij;
            tmpUx -= wzj * posYij / rij2 - wzj * pposYij / prij2;
            tmpUy += wzj * posXij / rij2 - wzj * pposXij / prij2;
        }
        Full_particle *pfi = f + i;
        pfi->vel.x = tmpUx;
        pfi->vel.y = tmpUy;
    }
}

void calc_BS_ep_sp(Full_particle *ep_i,
                                int n_ip,
                                fdps_spj_monopole *ep_j,
                                int n_jp,
                                Full_particle *f)
{
    int i, j;

    for (i=0; i<n_ip; i++) {
        double tmpUx = 0.0;
        double tmpUy = 0.0;
        Full_particle *pi = ep_i + i;
        const double posXi = pi->pos.x;
        const double posYi = pi->pos.y;
        printf("ok3 at outer.\n");
        // j < i
        for (j = 0; j < i; j++) {
            fdps_spj_monopole *pj = ep_j + j;
            // real
            const double posXj = pj->pos.x;
            const double posYj = pj->pos.y;
            const double rj2 = posXj * posXj + posYj * posYj;
            const double posXij = posXi - posXj;
            const double posYij = posYi - posYj;
            const double rij2 = posXij * posXij + posYij * posYij;
            const double invrij2 = 1.0 / rij2;
            const double wzj = pj->mass;
            // pseudo
            const double pposXj = posXj * R0 * R0 / rj2;
            const double pposYj = posYj * R0 * R0 / rj2;
            const double pposXij = posXi - pposXj;
            const double pposYij = posYi - pposYj;
            const double prij2 = pposXij * pposXij + pposYij * pposYij;
            const double pinvrij2 = 1.0 / prij2;
            tmpUx -= wzj * posYij * invrij2 - wzj * pposYij * pinvrij2;
            tmpUy += wzj * posXij * invrij2 - wzj * pposXij * pinvrij2;
        }
        // i=j pseudo only
        fdps_spj_monopole *pj = ep_j + i;
        const double posXj = pj->pos.x;
        const double posYj = pj->pos.y;
        const double rj2 = posXj * posXj + posYj * posYj;
        const double wzj = pj->mass;
        const double pposXj = posXj * R0 * R0 / rj2;
        const double pposYj = posYj * R0 * R0 / rj2;
        const double pposXij = posXi - pposXj;
        const double pposYij = posYi - pposYj;
        const double prij2 = pposXij * pposXij + pposYij * pposYij;
        const double pinvrij2 = 1.0 / prij2;
        tmpUx -= - wzj * pposYij * pinvrij2;
        tmpUy += - wzj * pposXij * pinvrij2;
        // j > i
        for (j = i + 1; j < n_jp; j++) {
            fdps_spj_monopole *pj = ep_j + j;
            // real
            const double posXj = pj->pos.x;
            const double posYj = pj->pos.y;
            const double rj2 = posXj * posXj + posYj * posYj;
            const double posXij = posXi - posXj;
            const double posYij = posYi - posYj;
            const double rij2 = posXij * posXij + posYij * posYij;
            const double invrij2 = 1.0 / rij2;
            const double wzj = pj->mass;
            // pseudo
            const double pposXj = posXj * R0 * R0 / rj2;
            const double pposYj = posYj * R0 * R0 / rj2;
            const double pposXij = posXi - pposXj;
            const double pposYij = posYi - pposYj;
            const double prij2 = pposXij * pposXij + pposYij * pposYij;
            const double pinvrij2 = 1.0 / prij2;
            tmpUx -= wzj * posYij * invrij2 - wzj * pposYij * pinvrij2;
            tmpUy += wzj * posXij * invrij2 - wzj * pposXij * pinvrij2;
        }
        Full_particle *pfi = f + i;
        pfi->vel.x = tmpUx;
        pfi->vel.y = tmpUy;
    }
}
#endif
