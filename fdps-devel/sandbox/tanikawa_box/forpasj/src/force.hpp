class Force {
public:
    PS::F64vec acc;
    void clear() {
        acc = 0.;
    }
};

void FP::copyFromForce(const Force & force) {
    this->acc = force.acc;
}

void FP::copyFromFP(const FP & fp) {
    this->mass = fp.mass;
    this->size = fp.size;
    this->pos  = fp.pos;
}

struct calcForce {

    void operator () (const FP * epi,
                      const PS::S32 nip,
                      const FP * epj,
                      const PS::S32 njp,
                      Force * force) {
        for(PS::S32 i = 0; i < nip; i++) {
            force[i].acc = epi[i].pos;
        }
    }

};

