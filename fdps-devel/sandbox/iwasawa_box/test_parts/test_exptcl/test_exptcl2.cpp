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
    const auto n_proc  = PS::Comm::getNumberOfProc();
    const auto my_rank = PS::Comm::getRank();
    PS::ParticleSystem<FP> ptcl;
    ptcl.initialize();
    PS::S32 n_loc = 128;
    PS::F64vec len_1d(2.0, 2.0, 2.0);
    ptcl.setNumberOfParticleLocal(n_loc);
    PS::MTTS mt;
    //mt.init_genrand(my_rank);
    for(auto i=0; i<n_loc; i++){
	ptcl[i].pos.x = (mt.genrand_res53() - 0.5) * len_1d.x;
	ptcl[i].pos.y = (mt.genrand_res53() - 0.5) * len_1d.y;
	ptcl[i].pos.z = (mt.genrand_res53() - 0.5) * len_1d.z;
    }
    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(-len_1d*0.5, len_1d*0.5);
    dinfo.decomposeDomainAll(ptcl);
    /*
    if(my_rank == 0){
	for(auto i=0; i<n_proc; i++){
	    std::cerr<<"i= "<<i<<" dinfo.getPosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
	}
    }
    */
    ptcl.exchangeParticle(dinfo);
    n_loc = ptcl.getNumberOfParticleLocal();
    for(auto i=0; i<n_loc; i++){
	if(dinfo.getPosDomain(my_rank).notContained(ptcl[i].pos)){
	    std::cerr<<"my_rank= "<<my_rank
		     <<" dinfo.getPosDomain(my_rank)= "<<dinfo.getPosDomain(my_rank)
		     <<" ptcl[i].pos= "<<ptcl[i].pos
		     <<std::endl;
	}
    }
    
    PS::Finalize();
    return 0;
}
