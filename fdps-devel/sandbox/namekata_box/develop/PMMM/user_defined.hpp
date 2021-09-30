#pragma once

class FP_nbody {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;
    PS::F64vec acc_pm;
    PS::F64    pot_pm;

    static PS::F64 eps;

    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F64vec & pos) {
        this->pos = pos;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::S64 getId() const {
        return this->id;
    }
    void copyFromFP(const FP_nbody & fp){ 
        this->id   = fp.id;
        this->mass = fp.mass;
        this->pos  = fp.pos;
    }
    void copyFromForce(const FP_nbody & force) {
        this->acc = force.acc;
        this->pot = force.pot;
    }
    void copyFromForcePMM(const PS::F64vec & acc, const PS::F64 & pot) {
        this->acc_pm = acc;
        this->pot_pm = pot;
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

template <class Tptcl>
void CalcForce(const FP_nbody * ep_i,
               const PS::S32 n_ip,
               const Tptcl * ep_j,
               const PS::S32 n_jp,
               FP_nbody * force) {
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec rij    = xi - ep_j[j].getPos();
            PS::F64    r3_inv = rij * rij;
            if (r3_inv == 0.0) continue;
            PS::F64    r_inv  = 1.0/sqrt(r3_inv);
            r3_inv  = r_inv * r_inv;
            r_inv  *= ep_j[j].getCharge();
            r3_inv *= r_inv;
            ai     += r3_inv * rij; // Coulomb's force 
            poti   += r_inv;        // Coulomb's potential
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}
