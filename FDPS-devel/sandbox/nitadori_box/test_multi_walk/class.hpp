#pragma once
#include <particle_simulator.hpp>

//#define N_THREAD_GPU 128

class Force{
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
    PS::S32 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const Force & force){
        acc = force.acc;
        pot = force.pot;
    }
    /*
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
                this->id, this->mass,
		this->pos.x, this->pos.y, this->pos.z,
		this->vel.x, this->vel.y, this->vel.z);
    }

    void readAscii(FILE* fp){
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
	       &this->pos.x, &this->pos.y, &this->pos.z,
	       &this->vel.x, &this->vel.y, &this->vel.z);
    }
    */
};

class Epi{
public:
    PS::S32 id;
    PS::F64 mass;
    PS::F64vec pos;
    //PS::F64vec vel;
    static PS::F64 eps;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FP & fp){ 
        pos = fp.pos;
        id = fp.id;
    }
};


class Epj{
public:
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64vec vel;
    PS::S32 id, pad;
    void copyFromFP(const FP & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getCharge() const { return mass; }
    
};

