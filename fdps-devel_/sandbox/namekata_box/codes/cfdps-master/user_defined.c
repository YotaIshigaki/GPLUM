//user_defined.c
#include "math.h"
#include "./FDPS_vector.h"
#include "FDPS_types.h"
#include "user_defined.h"

Cvec_Float64 Cvec_Float64_zero(){
    Cvec_Float64 zero;
    zero.x = zero.y = zero.z = 0;
    return zero;
}

void  calc_gravity_c_epep(PFull_particle ep_i,
			   int n_ip,
			   PFull_particle ep_j,
			   int n_jp,
			   PFull_particle f)
{
    int i, j;
    for(i=0;i<n_ip;i++){

	PFull_particle pi = ep_i + i;
	double eps2 = pi->eps*pi->eps;
	double xi = pi->pos.x;
	double yi = pi->pos.y;
	double zi = pi->pos.z;
	double ax, ay, az, pot;
	ax = ay = az = pot = 0;
	for(j=0;j<n_jp;j++){
	    PFull_particle pj = ep_j + j;
	    double dx = xi - pj->pos.x;
	    double dy = yi - pj->pos.y;
	    double dz = zi - pj->pos.z;
            double r2 = dx*dx+dy*dy+dz*dz+eps2;
	    double rinv = 1.0/sqrt(r2);
	    double mrinv = pj->mass* rinv;
	    double mr3inv = mrinv*rinv*rinv;
	    ax -= dx*mr3inv;
	    ay -= dy*mr3inv;
	    az -= dz*mr3inv;
	    pot = pot - mrinv;
	}
	PFull_particle  pfi = f+i;
	pfi->pot += pot;
	pfi->acc.x += ax;
	pfi->acc.y += ay;
	pfi->acc.z += az;
    }
}

void  calc_gravity_c_epsp(PFull_particle ep_i,
			   int n_ip,
			   PSPJMonopole ep_j,
			   int n_jp,
			   PFull_particle f)
{
    int i, j;
    for(i=0;i<n_ip;i++){

	PFull_particle pi = ep_i + i;
	double eps2 = pi->eps*pi->eps;
	double xi = pi->pos.x;
	double yi = pi->pos.y;
	double zi = pi->pos.z;
	double ax, ay, az, pot;
	ax = ay = az = pot = 0;
	for(j=0;j<n_jp;j++){
	    PSPJMonopole pj = ep_j + j;
	    double dx = xi - pj->pos.x;
	    double dy = yi - pj->pos.y;
	    double dz = zi - pj->pos.z;
            double r2 = dx*dx+dy*dy+dz*dz+eps2;
	    double rinv = 1.0/sqrt(r2);
	    double mrinv = pj->mass* rinv;
	    double mr3inv = mrinv*rinv*rinv;
	    ax -= dx*mr3inv;
	    ay -= dy*mr3inv;
	    az -= dz*mr3inv;
	    pot = pot - mrinv;
	}
	PFull_particle  pfi = f+i;
	pfi->pot += pot;
	pfi->acc.x += ax;
	pfi->acc.y += ay;
	pfi->acc.z += az;
    }
}
