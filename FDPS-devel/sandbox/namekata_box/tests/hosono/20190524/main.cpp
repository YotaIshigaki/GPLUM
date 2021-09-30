#include <particle_simulator.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>


//Density summation
class Force{
	public:
	void clear(){
	}
};

class RealPtcl{
	public:
	PS::F64 mass;
	PS::F64vec pos;
	PS::F64 rad;

	//Constructor
	RealPtcl(){
	}
	//Copy functions
	void copyFromForce(const Force& dens){
	}
	//Give necessary values to FDPS
	PS::F64 getCharge() const{
		return this->mass;
	}
	PS::F64vec getPos() const{
		return this->pos;
	}
	void setPos(const PS::F64vec& pos){
		this->pos = pos;
	}
};

class EPI{
	public:
	PS::F64vec pos;
	PS::F64    mass;
	PS::F64    rad;
	void copyFromFP(const RealPtcl& rp){
		this->pos  = rp.pos;
		this->mass = rp.mass;
		this->rad = rp.rad;
	}
	PS::F64vec getPos() const{
		return this->pos;
	}
	PS::F64 getRSearch() const{
		return this->rad;
	}
	void setPos(const PS::F64vec& pos){
		this->pos = pos;
	}
};

class EPJ{
    public:
	PS::F64    mass;
	PS::F64vec pos;
	PS::F64    rad;
	void copyFromFP(const RealPtcl& rp){
		this->mass = rp.mass;
		this->pos  = rp.pos;
		this->rad  = rp.rad;
	}
	PS::F64vec getPos() const{
		return this->pos;
	}
	void setPos(const PS::F64vec& pos){
		this->pos = pos;
	}
	PS::F64 getRSearch() const{
		return this->rad;
	}
};

class ShortRangeForce{
	public:
	void operator () (const EPI* const ep_i, const PS::S32 Nip, const EPJ* const ep_j, const PS::S32 Njp, Force* const dens){
	}
};

int main(int argc, char* argv[]){
	PS::Initialize(argc, argv);
	PS::ParticleSystem<RealPtcl> ptcl;
	ptcl.initialize();
	PS::DomainInfo dinfo;
	dinfo.initialize();
	ptcl.setAverageTargetNumberOfSampleParticlePerProcess(200);

	dinfo.decomposeDomainAll(ptcl);
	ptcl.exchangeParticle(dinfo);
	PS::TreeForForceShort<Force, EPI, EPJ>::Gather short_tree;

	short_tree.initialize(ptcl.getNumberOfParticleLocal());
	short_tree.calcForceAllAndWriteBack(ShortRangeForce(), ptcl, dinfo);

	PS::Finalize();
	return 0;
}

