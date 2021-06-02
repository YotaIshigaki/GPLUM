#pragma once

class FP_nbody {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;    

    static PS::F64 eps;

    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64 getChargePMMM() const {
        return this->mass;
    }
    void copyFromFP(const FP_nbody & fp){ 
        this->mass = fp.mass;
        this->pos  = fp.pos;
    }
    void copyFromForce(const FP_nbody & force) {
        this->acc = force.acc;
        this->pot = force.pot;
    }
    void copyFromForcePMMM(const PS::F64vec & acc,
                           const PS::F64 & pot) {
        this->acc = acc;
        this->pot = pot;
    }
    void clear() {
        this->acc = 0.0;
        this->pot = 0.0;
    }
    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
    }
};
