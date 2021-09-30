#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include "particle_simulator.hpp"
#include "particle_mesh.hpp"
//#include<particle_simulator.hpp>
//#include<particle_mesh.hpp>

/*
class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64 pot;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

class FP{
public:
    PS::S64 id;
    PS::F64 charge;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc_direct;
    PS::F64vec acc_pm;
    PS::F64 pot_direct;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const ForceGrav & force){
        acc_direct = force.acc;
        pot_direct = force.pot;
    }
    void writeAscii(FILE* fp) const{
	fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                this->id, this->charge, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z);
    }

    void readAscii(FILE* fp){
	fscanf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->charge, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z);
    }
    void copyFromForceParticleMesh( const PS::F32vec & force){
	acc_pm = force;
    }
    PS::F64 getChargeParticleMesh() const {
	return charge;
    }
};
*/

int main(int argc, char *argv[]){
    PS::Initialize(argc, argv);
    /*
    PS::DomainInfo dinfo;
    PS::ParticleSystem<FP> system;
    dinfo.initialize();
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(0.0, 1.0);
    dinfo.decomposeDomainAll(system);
    */
    PS::PM::ParticleMesh pm;
    //pm.calcForceAllAndWriteBack(system, dinfo);
    
    //dinfo.collectSampleParticle(system);

    PS::Finalize();
    return 0;
}
