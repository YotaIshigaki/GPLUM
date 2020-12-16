#pragma once
#include <cstdint>

class Force_clus_id {
public:
    void clear() { }
};

class FP_nbody {
public:
    PS::S64    id;
    int16_t    type;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    dt;
    PS::F64    h;

    PS::S64 getId() const {
        return id;
    }

    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromForce(const Force_clus_id & f) { }

};

class EP_clus_id {
public:
    PS::S64 id;
    int16_t type;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64 dt;
    PS::F64 h;

    void copyFromFP(const FP_nbody & fp) {
        this->id   = fp.id;
        this->type = fp.type;
        this->mass = fp.mass;
        this->pos  = fp.pos;
        this->dt   = fp.dt;
        this->h    = fp.h;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getRSearch() const {
        return this->h;
    }
    PS::F64 getTimestep() const {
        return this->dt;
    }
    void setPos(const PS::F64vec & pos) {
        this->pos = pos;
    }
};

void CalcForceEmpty (const EP_clus_id * ep_i,
                     const PS::S32 n_ip,
                     const EP_clus_id * ep_j,
                     const PS::S32 n_jp,
                     Force_clus_id * force) {
}
