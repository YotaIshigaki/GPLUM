class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp) {
		fscanf(fp, "%lf\n", &time);
		fscanf(fp, "%lld\n", &n_body);
		return n_body;
    }
    void writeAscii(FILE* fp) const {
	fprintf(fp, "%e\n", time);
	fprintf(fp, "%lld\n", n_body);
    }
};

class FPGrav{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;    

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

	void writeAscii(FILE* fp) const {
		fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
        /*
        fprintf(fp, "%5lld %+.16e %+.16e %+.16e %+.16e\n",
                this->id, this->acc[0], this->acc[1], this->acc[2], this->pot);
        */
	}

	void readAscii(FILE* fp) {
		fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
	}

};

PS::F64 FPGrav::eps = 1.0/32.0;

#ifdef ENABLE_PHANTOM_GRAPE_X86


template <class TParticleJ>
void CalcGravity(const FPGrav * iptcl,
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
        xi[i][0] = iptcl[i].getPos()[0];
        xi[i][1] = iptcl[i].getPos()[1];
        xi[i][2] = iptcl[i].getPos()[2];
        ai[i][0] = 0.0;
        ai[i][1] = 0.0;
        ai[i][2] = 0.0;
        pi[i]    = 0.0;
    }
    for(PS::S32 j = 0; j < nj; j++) {
        xj[j][0] = jptcl[j].getPos()[0];
        xj[j][1] = jptcl[j].getPos()[1];
        xj[j][2] = jptcl[j].getPos()[2];
        mj[j]    = jptcl[j].getCharge();
	xj[j][0] = jptcl[j].pos[0];
        xj[j][1] = jptcl[j].pos[1];
        xj[j][2] = jptcl[j].pos[2];
        mj[j]    = jptcl[j].mass;
    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
    PS::S32 devid = omp_get_thread_num();
#else
    PS::S32 devid = 0;
#endif
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

#elif defined ENABLE_INTRINSICS

#include "vector_x86.hpp"

template <class TParticleJ>
void CalcGravity(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {    
    PS::S32 nvl = v8sf::getVectorLength();
    v8sf eps2(FPGrav::eps * FPGrav::eps);
    for(PS::S32 i = 0; i < n_ip; i += nvl) {
        v8sf pxi(ep_i[i  ].pos[0], ep_i[i+1].pos[0], ep_i[i+2].pos[0], ep_i[i+3].pos[0],
                 ep_i[i+4].pos[0], ep_i[i+5].pos[0], ep_i[i+6].pos[0], ep_i[i+7].pos[0]);
        v8sf pyi(ep_i[i  ].pos[1], ep_i[i+1].pos[1], ep_i[i+2].pos[1], ep_i[i+3].pos[1],
                 ep_i[i+4].pos[1], ep_i[i+5].pos[1], ep_i[i+6].pos[1], ep_i[i+7].pos[1]);
        v8sf pzi(ep_i[i  ].pos[2], ep_i[i+1].pos[2], ep_i[i+2].pos[2], ep_i[i+3].pos[2],
                 ep_i[i+4].pos[2], ep_i[i+5].pos[2], ep_i[i+6].pos[2], ep_i[i+7].pos[2]);
        v8sf axi(0.0);
        v8sf ayi(0.0);
        v8sf azi(0.0);        
        v8sf pot(0.0);
        for(PS::S32 j = 0; j < n_jp; j++){
            v8sf dx = pxi - ep_j[j].pos[0];
            v8sf dy = pyi - ep_j[j].pos[1];
            v8sf dz = pzi - ep_j[j].pos[2];

            v8sf r2 = eps2 + dx * dx + dy * dy + dz * dz;

            v8sf ri1 = v8sf::rsqrt_0th(r2);
            v8sf mr1 = ri1 * ep_j[j].mass;
            pot -= mr1;

            v8sf mr3 = mr1 * ri1 * ri1;
            axi -= mr3 * dx;
            ayi -= mr3 * dy;
            azi -= mr3 * dz;
        }

        PS::F32 buf0[nvl], buf1[nvl], buf2[nvl], buf3[nvl];
        axi.store(buf0);
        ayi.store(buf1);
        azi.store(buf2);
        pot.store(buf3);

        PS::S32 nii = (n_ip - i > nvl) ? nvl : n_ip - i;
        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf0[ii];
            force[i+ii].acc[1] += buf1[ii];
            force[i+ii].acc[2] += buf2[ii];
            force[i+ii].pot    += buf3[ii];
        }        

    }

}

#else

template <class TParticleJ>
void CalcGravity(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec rij    = xi - ep_j[j].getPos();
            PS::F64    r3_inv = rij * rij + eps2;
            PS::F64    r_inv  = 1.0/sqrt(r3_inv);
            r3_inv  = r_inv * r_inv;
            r_inv  *= ep_j[j].getCharge();
            r3_inv *= r_inv;
            ai     -= r3_inv * rij;
            poti   -= r_inv;
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}

#endif
