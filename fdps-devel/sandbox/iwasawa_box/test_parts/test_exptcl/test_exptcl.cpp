#include<iostream>
#include<fstream>
#include<unistd.h>
#include<random>
#include<particle_simulator.hpp>

class FP{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec getPos() const { return pos; }
    PS::F64vec getCharge() const { return mass; }
};

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
    /*
    const PS::S32 nx = 7;
    const PS::S32 ny = 13;
    const PS::S32 nz = 2;
    PS::S32 n_proc = nx*ny*nz;
    PS::S32vec n_domain(nx, ny, nz);
    for(auto i=0; i<n_proc; i++){
	PS::S32vec rank_1d = PS::GetRank1D(n_domain, i, 3);
	std::cerr<<"i= "<<i<<" rank_1d= "<<rank_1d<<std::endl;
    }
    */
    
    PS::ParticleSystem<FP> ptcl;
    ptcl.initialize();
    PS::S32 n_loc = 128;
    PS::F64vec len_1d(2.0, 2.0, 2.0);
    PS::S32 n_domain[] = {8, 8, 8};
    ptcl.setNumberOfParticleLocal(n_loc);
    PS::MTTS mt;
    for(auto i=0; i<n_loc; i++){
	ptcl[i].pos.x = (mt.genrand_res53() - 0.5) * len_1d.x;
	ptcl[i].pos.y = (mt.genrand_res53() - 0.5) * len_1d.y;
	ptcl[i].pos.z = (mt.genrand_res53() - 0.5) * len_1d.z;
    }
    PS::S32 n_proc = n_domain[0] * n_domain[1] * n_domain[2];
    PS::F64ort * domain = new PS::F64ort[n_proc];
    for(auto i=0; i<n_proc; i++){
	const auto rank_1d = PS::GetRankVec(n_domain, i, 3);
	for(auto k=0; k<3; k++){
	    domain[i].low_[k]  = -len_1d[k]*0.5 + (len_1d[k]/n_domain[k])*rank_1d[k];
	    domain[i].high_[k] = -len_1d[k]*0.5 + (len_1d[k]/n_domain[k])*(rank_1d[k]+1);
	}
	//std::cerr<<"i= "<<i<<" domain[i]= "<<domain[i]<<std::endl;
    }
    
    PS::S32vec rank_new_vec, my_dnp;
    PS::F64vec my_dis_domain;
    PS::F64vec len_peri = {0.0, 0.0, 0.0};
    for(auto i=0; i<n_loc; i++){
	PS::S32 my_rank_glb = mt.genrand_int31() % n_proc;
	std::cerr<<" ************ "<<std::endl;
	std::cerr<<"i= "<<i<<" my_rank_glb= "<<my_rank_glb<<" ptcl[i].pos= "<<ptcl[i].pos<<std::endl;
	for(auto k=0; k<3; k++){
	    CalcRank1D(rank_new_vec[k], my_dnp[k], my_dis_domain[k], ptcl[i].pos, my_rank_glb, n_domain, domain, len_peri[k], k);
	}
	std::cerr<<"rank_new_vec= "<<rank_new_vec<<std::endl;
	std::cerr<<"my_dnp= "<<my_dnp<<std::endl;
	std::cerr<<"my_dis_domain= "<<my_dis_domain<<std::endl;
    }




    
    PS::Finalize();
    return 0;
}
