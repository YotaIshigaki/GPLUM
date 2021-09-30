#include "particle_simulator.hpp"

class Force{
public:
    PS::F64vec acc;
    PS::F64 pot;
    PS::S64 nj;
    PS::S64 njreal;
    void clear(){
        acc = 0.0;
        pot = 0.0;
        nj  = 0;
        njreal = 0;
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
        this->nj  = force.nj;
        this->njreal = force.njreal;
    }
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::S64 nj;
    PS::S64 njreal;
    static PS::F64 cutoff;
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
        fout<<"acc="<<acc<<std::endl;
        fout<<"pot="<<pot<<std::endl;
    }
};

PS::F64 FP::cutoff = 1.0 / 8.0;

class EPI{
public:
    PS::F64vec pos;
    PS::S64 id;
    PS::S64 nj;
    PS::S64 njreal;
    static PS::F64 eps;
    void copyFromFP(const FP & fp){ 
        pos = fp.pos;
        id  = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
};

PS::F64 EPI::eps    = 1.0/32.0;

class EPJ{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    static PS::F64 cutoff;
    void copyFromFP(const FP & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
        id   = fp.id;
        cutoff = fp.cutoff;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
    PS::F64 getRSearch() {
        return this->cutoff;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

PS::F64 EPJ::cutoff = 0.;

struct calcForceEpEp {
    void operator () (const EPI *epi,
                      const PS::S32 nip,
                      const EPJ *epj,
                      const PS::S32 njp,
                      Force *force) {

        PS::F64 eps2 = epi[0].eps * epi[0].eps;
        PS::F64 cutoff2 = eps2 + epj[0].cutoff * epj[0].cutoff;
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S64 idi     = epi[i].id;
            PS::F64vec posi = epi[i].pos;
            PS::F64vec acci = 0.;
            PS::F64 poti    = 0.;
            PS::S64 njreal  = 0;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S64 idj     = epj[j].id;
                PS::F64 mj      = epj[j].mass;
                PS::F64vec posj = epj[j].pos;
                PS::F64vec dx   = posj - posi;

                PS::F64 r2   = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + eps2;
                PS::F64 rinv = ((idj != idi) ? 1. / sqrt(r2) : 0.);
                rinv         = ((r2 < cutoff2) ? rinv : 0.);
                njreal      += ((r2 < cutoff2) ? 1 : 0);
                PS::F64 pot  = mj * rinv;
                poti -= pot;

                PS::F64 fij = pot * rinv * rinv;
                acci += fij * dx;
            }
            force[i].acc += acci;
            force[i].pot += poti;
            force[i].nj  += njp;
            force[i].njreal += njreal;

        }

    }
};

struct calcForceEpSp {
    void operator () (const EPI *epi,
                      const PS::S32 nip,
                      const PS::SPJMonoPolePeriodic *spj,
                      const PS::S32 njp,
                      Force *force) {

        PS::F64 eps2 = epi[0].eps * epi[0].eps;
        PS::F64 cutoff2 = eps2 + FP::cutoff * FP::cutoff;
        for(PS::S32 i = 0; i < nip; i++) {
            PS::F64vec posi = epi[i].pos;
            PS::F64vec acci = 0.;
            PS::F64 poti    = 0.;
            PS::S64 njreal  = 0;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::F64 mj      = spj[j].mass;
                PS::F64vec posj = spj[j].pos;
                PS::F64vec dx   = posj - posi;

                PS::F64 r2   = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + eps2;
                PS::F64 rinv = 1. / sqrt(r2);
                rinv         = ((r2 < cutoff2) ? rinv : 0.);
                njreal      += ((r2 < cutoff2) ? 1 : 0);
                PS::F64 pot  = mj * rinv;
                poti -= pot;

                PS::F64 fij = pot * rinv * rinv;
                acci += fij * dx;
            }
            force[i].acc += acci;
            force[i].pot += poti;
            force[i].nj  += njp;
            force[i].njreal += njreal;
        }
    }
};
