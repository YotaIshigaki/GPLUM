class FPGrav{
public:
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;    
    PS::F64vec vel2;

    static PS::F64 eps;

    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromFP(const FPGrav & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
    }

    void copyFromForce(const FPGrav & force) {
        acc = force.acc;
        pot = force.pot;
    }    

    void clear() {
        acc = 0.0;
        pot = 0.0;
    }

    void predict(PS::F32 dt) {
        pos  = pos  +       vel * dt + 0.5 * acc * dt * dt;
        vel2 = vel  + 0.5 * acc * dt;
    }

    void correct(PS::F32 dt) {
        vel  = vel2 + 0.5 * acc * dt;
    }

	void writeAscii(FILE* fp) const {
		fprintf(fp, "%+e %+e %+e %+e %+e %+e %+e\n", 
                this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
	}

};

PS::F64 FPGrav::eps = 1.0 / 32.0;

template <class TParticleJ>
struct CalcGravity{

    void operator () (const FPGrav * iptcl,
                      const PS::S32 ni,
                      const TParticleJ * jptcl,
                      const PS::S32 nj,
                      FPGrav * force) {

        const PS::S32 nipipe = ni;
        const PS::S32 njpipe = nj;
        PS::F64 (*xi)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
        PS::F64 (*ai)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
        PS::F64  *pi     = (PS::F64  *    )malloc(sizeof(PS::F64) * nipipe);
        PS::F64 (*xj)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * njpipe * PS::DIMENSION);
        PS::F64  *mj     = (PS::F64  *    )malloc(sizeof(PS::F64) * njpipe);

        for(PS::S32 i = 0; i < ni; i++) {
            xi[i][0] = iptcl[i].pos[0];
            xi[i][1] = iptcl[i].pos[1];
            xi[i][2] = iptcl[i].pos[2];
            ai[i][0] = 0.0;
            ai[i][1] = 0.0;
            ai[i][2] = 0.0;
            pi[i]    = 0.0;
        }

        for(PS::S32 j = 0; j < nj; j++) {
            xj[j][0] = jptcl[j].pos[0];
            xj[j][1] = jptcl[j].pos[1];
            xj[j][2] = jptcl[j].pos[2];
            mj[j]    = jptcl[j].mass;
        }

        PS::S32 devid = omp_get_thread_num();
        g5_set_xmjMC(devid, 0, nj, xj, mj);
        g5_set_nMC(devid, nj);
        g5_calculate_force_on_xMC(devid, xi, ai, pi, ni);

        for(PS::S32 i = 0; i < ni; i++) {
            force[i].acc[0] += ai[i][0];
            force[i].acc[1] += ai[i][1];
            force[i].acc[2] += ai[i][2];
            force[i].pot    -= pi[i];
        }

        free(xi);
        free(ai);
        free(pi);
        free(xj);
        free(mj);
          
    }
};

