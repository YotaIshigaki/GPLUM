class Force;

class FullParticle {
public:
    PS::F64    mass;
    PS::F64    size;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;

    PS::F64vec getPos() const {
        return this->pos;
    }

    void setPos(PS::F64vec pos_new) {
        this->pos = pos_new;
    }
    
    PS::F64 getCharge() const {
        return this->mass;
    }

    PS::F64 getRSearch() const {
        return this->size;
    }

    void copyFromForce(const Force & force);

    void copyFromFP(const FullParticle & fp);

    void integrateOrbit(PS::F64 dt) {
        this->pos += this->vel * dt;
    }

};

typedef FullParticle FP;
