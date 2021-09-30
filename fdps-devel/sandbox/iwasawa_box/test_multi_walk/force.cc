//#include <iostream>
//#include <particle_simulator.hpp>
//#include "class.hpp"
#include "force.hpp"

static const int g_n_walk_limit = 128;
static const int g_ni_limit = 1024;
static Force g_force[g_n_walk_limit][g_ni_limit];

#if 1

int DispatchKernelWithSP(PS::S32    nwalk,
                         const Epi     ** epi,
                         const PS::S32  * n_epi,
                         const Epj     ** epj,
                         const PS::S32  * n_epj,
                         const PS::SPJMonopole ** spj,
                         const PS::S32  * n_spj){
    PS::F32 eps2 = Epi::eps * Epi::eps;
    //std::cerr<<"nwalk="<<nwalk<<std::endl;
    for(int iw=0; iw<nwalk; iw++){
	for(int ip=0; ip<n_epi[iw]; ip++){
	    g_force[iw][ip].acc = 0.0;
	    g_force[iw][ip].pot = 0.0;
	    //std::cerr<<"n_epj[iw]="<<n_epj[iw]<<std::endl;
	    //std::cerr<<"epi[iw][ip].id="<<epi[iw][ip].id<<std::endl;
	    for(int jp=0; jp<n_epj[iw]; jp++){
		//std::cerr<<"ip="<<ip<<" jp="<<jp<<std::endl;
		//std::cerr<<"epj[iw][jp].id="<<epj[iw][jp].id<<std::endl;
		if( epi[iw][ip].id == epj[iw][jp].id) continue;
		const PS::F64 mj = epj[iw][jp].mass;
		PS::F64vec rij = epi[iw][ip].pos - epj[iw][jp].pos;
		PS::F64 r2 = rij*rij + eps2;
		PS::F64 inv_r = 1.0 / sqrt(r2);
		PS::F64 inv_r3 = inv_r * inv_r * inv_r;
		g_force[iw][ip].acc -= mj * inv_r3 *rij;
		g_force[iw][ip].pot -= mj * inv_r;
	    }
	    //std::cerr<<"n_spj[iw]="<<n_spj[iw]<<std::endl;
	    for(int jp=0; jp<n_spj[iw]; jp++){
		//std::cerr<<"jp="<<jp<<std::endl;
		const PS::F64 mj = spj[iw][jp].getCharge();
		PS::F64vec rij = epi[iw][ip].pos - spj[iw][jp].getPos();
		PS::F64 r2 = rij*rij + eps2;
		PS::F64 inv_r = 1.0 / sqrt(r2);
		PS::F64 inv_r3 = inv_r * inv_r * inv_r;
		g_force[iw][ip].acc -= mj * inv_r3 *rij;
		g_force[iw][ip].pot -= mj * inv_r;
	    }
	}
    }
    return 0;
}

int RetrieveKernel(const PS::S32    nwalk,
                   const PS::S32  * ni,
                   Force         ** force){
    for(int iw=0; iw<nwalk; iw++){
	for(int ip=0; ip<ni[iw]; ip++){
	    force[iw][ip] = g_force[iw][ip];
	}
    }
    return 0;
}

#else
int DispatchKernel::operator () (PS::S32    nwalk,
				 const Epi     ** epi,
				 const PS::S32  * ni,
				 const Epj     ** epj,
				 const PS::S32  * nj){
    PS::F32 eps2 = Epi::eps * Epi::eps;
    for(int iw=0; iw<nwalk; iw++){
	for(int ip=0; ip<ni[iw]; ip++){
	    g_force[iw][ip].acc = 0.0;
	    g_force[iw][ip].pot = 0.0;
	    for(int jp=0; jp<nj[iw]; jp++){
		const PS::F64 mj = epj[iw][jp].mass;
		PS::F64vec rij = epi[iw][ip].pos - epj[iw][jp].pos;
		PS::F64 r2 = rij*rij + eps2;
		PS::F64 inv_r = 1.0 / sqrt(r2);
		PS::F64 inv_r3 = inv_r * inv_r * inv_r;
		g_force[iw][ip].acc -= mj * inv_r3 *rij;
		g_force[iw][ip].pot -= mj * inv_r;
	    }
	}
    }
    return 0;
}

int DispatchKernelWithSP::operator () (PS::S32    nwalk,
				       const Epi     ** epi,
				       const PS::S32  * n_epi,
				       const Epj     ** epj,
				       const PS::S32  * n_epj,
				       const PS::SPJMonopole ** spj,
				       const PS::S32  * n_spj){
    PS::F32 eps2 = Epi::eps * Epi::eps;
    //std::cerr<<"nwalk="<<nwalk<<std::endl;
    for(int iw=0; iw<nwalk; iw++){
	for(int ip=0; ip<n_epi[iw]; ip++){
	    g_force[iw][ip].acc = 0.0;
	    g_force[iw][ip].pot = 0.0;
	    //std::cerr<<"n_epj[iw]="<<n_epj[iw]<<std::endl;
	    //std::cerr<<"epi[iw][ip].id="<<epi[iw][ip].id<<std::endl;
	    for(int jp=0; jp<n_epj[iw]; jp++){
		//std::cerr<<"ip="<<ip<<" jp="<<jp<<std::endl;
		//std::cerr<<"epj[iw][jp].id="<<epj[iw][jp].id<<std::endl;
		if( epi[iw][ip].id == epj[iw][jp].id) continue;
		const PS::F64 mj = epj[iw][jp].mass;
		PS::F64vec rij = epi[iw][ip].pos - epj[iw][jp].pos;
		PS::F64 r2 = rij*rij + eps2;
		PS::F64 inv_r = 1.0 / sqrt(r2);
		PS::F64 inv_r3 = inv_r * inv_r * inv_r;
		g_force[iw][ip].acc -= mj * inv_r3 *rij;
		g_force[iw][ip].pot -= mj * inv_r;
	    }
	    //std::cerr<<"n_spj[iw]="<<n_spj[iw]<<std::endl;
	    for(int jp=0; jp<n_spj[iw]; jp++){
		//std::cerr<<"jp="<<jp<<std::endl;
		const PS::F64 mj = spj[iw][jp].getCharge();
		PS::F64vec rij = epi[iw][ip].pos - spj[iw][jp].getPos();
		PS::F64 r2 = rij*rij + eps2;
		PS::F64 inv_r = 1.0 / sqrt(r2);
		PS::F64 inv_r3 = inv_r * inv_r * inv_r;
		g_force[iw][ip].acc -= mj * inv_r3 *rij;
		g_force[iw][ip].pot -= mj * inv_r;
	    }
	}
    }
    return 0;
}

int RetrieveKernel::operator () (const PS::S32    nwalk,
				 const PS::S32  * ni,
				 Force         ** force){
    for(int iw=0; iw<nwalk; iw++){
	for(int ip=0; ip<ni[iw]; ip++){
	    force[iw][ip] = g_force[iw][ip];
	}
    }
    return 0;
}

#endif
