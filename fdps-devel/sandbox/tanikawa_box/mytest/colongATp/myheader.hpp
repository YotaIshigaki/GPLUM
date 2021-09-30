#include "particle_simulator.hpp"

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
    PS::F64vec pos;
    PS::F64vec getPos() const {
        return this->pos;
    }
    void copyFromForce(const Force & force){
        this->acc = force.acc;
        this->pot = force.pot;
    }
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
        fout<<"acc="<<acc<<std::endl;
        fout<<"pot="<<pot<<std::endl;
    }
};

class EPI{
public:
    PS::F64vec pos;
    PS::S64 id;
    static PS::F64 eps;
    void copyFromFP(const FP & fp){ 
        pos = fp.pos;
        id  = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
};

PS::F64 EPI::eps = 1.0/32.0;

class EPJ{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    void copyFromFP(const FP & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
        id   = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

class Moment{
public:
    PS::F64vec pos;
    PS::F64 mass;
    void init(){
        pos  = 0.0;
        mass = 0.0;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    void accumulateAtLeaf(const EPJ & epj){
        mass += epj.mass;
        pos  += epj.mass * epj.pos;
    }
    void set(){
        pos = pos / mass;
    }
    void accumulate(const Moment & mom){
        mass += mom.mass;
        pos  += mom.mass * mom.pos;
    }
    // if you want to use
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

class SPJ{
public:
    void copyFromMoment(const Moment & mom){
        mass = mom.mass;
        pos  = mom.pos;
    }
    void clear(){
        mass = 0.0;
        pos  = 0.0;
    }
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec getPos() const {
        return this->pos;
    }
    Moment convertToMoment() {
        Moment mom;
        mom.mass = mass;
        mom.pos  = pos;
        return mom;
    }
};


struct calcForceEpEp {
    void operator () (const EPI *epi,
                      const PS::S32 nip,
                      const EPJ *epj,
                      const PS::S32 njp,
                      Force *force) {

        PS::F64 eps2 = epi[0].eps * epi[0].eps;
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S64 idi     = epi[i].id;
            PS::F64vec posi = epi[i].pos;
            PS::F64vec acci = 0.;
            PS::F64 poti    = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S64 idj     = epj[j].id;
                PS::F64 mj      = epj[j].mass;
                PS::F64vec posj = epj[j].pos;
                PS::F64vec dx   = posj - posi;

                PS::F64 r2   = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + eps2;
                PS::F64 rinv = ((idj != idi) ? 1. / sqrt(r2) : 0.);
                PS::F64 pot  = mj * rinv;
                poti -= pot;

                PS::F64 fij = pot * rinv * rinv;
                acci += fij * dx;
            }
            force[i].acc += acci;
            force[i].pot += poti;
        }
    }
};

struct calcForceEpSp {
    void operator () (const EPI *epi,
                      const PS::S32 nip,
                      const SPJ *spj,
                      const PS::S32 njp,
                      Force *force) {

        PS::F64 eps2 = epi[0].eps * epi[0].eps;
        for(PS::S32 i = 0; i < nip; i++) {
            PS::F64vec posi = epi[i].pos;
            PS::F64vec acci = 0.;
            PS::F64 poti    = 0.;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::F64 mj      = spj[j].mass;
                PS::F64vec posj = spj[j].pos;
                PS::F64vec dx   = posj - posi;

                PS::F64 r2   = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + eps2;
                PS::F64 rinv = 1. / sqrt(r2);
                PS::F64 pot  = mj * rinv;
                poti -= pot;

                PS::F64 fij = pot * rinv * rinv;
                acci += fij * dx;
            }
            force[i].acc += acci;
            force[i].pot += poti;
        }
    }
};

