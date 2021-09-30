/* C++ header */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
/* FDPS header */
#include <particle_simulator.hpp>

//Density summation
class Force{
public:
    void clear() {}
};

class FP {
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64 rad;

    void copyFromForce(const Force& f){ }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64vec getPos() const{
        return this->pos;
    }
    void setPos(const PS::F64vec& pos){
        this->pos = pos;
    }
    void writeRestartData(FILE *fp) const { }
};

class EPI {
public:
    PS::F64vec pos;
    PS::F64    mass;
    PS::F64    rad;
    void copyFromFP(const FP& fp){
        this->pos  = fp.pos;
        this->mass = fp.mass;
        this->rad = fp.rad;
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
    void copyFromFP(const FP& fp){
        this->mass = fp.mass;
        this->pos  = fp.pos;
        this->rad  = fp.rad;
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

class CalcForce {
public:
    void operator () (const EPI * const ep_i,
                      const PS::S32 n_ip,
                      const EPJ * const ep_j,
                      const PS::S32 n_jp,
                      Force * const f){
    }
};

int main (int argc, char* argv[]) {
    PS::Initialize(argc, argv);
    PS::ParticleSystem<FP> psys;
    psys.initialize();
    
    char basename[256], format[256];
    sprintf(format,"%%s_np%%05d_r%%05d_rst.dat");
    sprintf(basename, "result/nbody_sn%05d_", 0); 
    psys.writeParticleBinary(basename, format, &FP::writeRestartData);

    PS::Finalize();
    return 0;
}

