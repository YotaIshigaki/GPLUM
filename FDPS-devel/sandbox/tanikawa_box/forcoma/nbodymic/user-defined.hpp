PS::F64 getWallclockTime();
static PS::F64 wt_kernel[320];


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
        /*
		fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
        */
		fprintf(fp, "hoge %8d %+e %+e %+e %+e %+e %+e %+e %+e\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->acc.x, this->acc.y, this->acc.z,
                this->pot);
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

#elif defined USE_INTRINSICS

#if 0

template <class TParticleJ>
void CalcGravity(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {
    const PS::S32 nvector = 16;

    v16sf ep2(FPGrav::eps * FPGrav::eps);
    for(PS::S32 i = 0; i < n_ip; i += nvector){
        v16sf pxi(ep_i[i+ 0].pos[0], ep_i[i+ 1].pos[0], ep_i[i+ 2].pos[0], ep_i[i+ 3].pos[0],
                  ep_i[i+ 4].pos[0], ep_i[i+ 5].pos[0], ep_i[i+ 6].pos[0], ep_i[i+ 7].pos[0],
                  ep_i[i+ 8].pos[0], ep_i[i+ 9].pos[0], ep_i[i+10].pos[0], ep_i[i+11].pos[0],
                  ep_i[i+12].pos[0], ep_i[i+13].pos[0], ep_i[i+14].pos[0], ep_i[i+15].pos[0]);
        v16sf pyi(ep_i[i+ 0].pos[1], ep_i[i+ 1].pos[1], ep_i[i+ 2].pos[1], ep_i[i+ 3].pos[1],
                  ep_i[i+ 4].pos[1], ep_i[i+ 5].pos[1], ep_i[i+ 6].pos[1], ep_i[i+ 7].pos[1],
                  ep_i[i+ 8].pos[1], ep_i[i+ 9].pos[1], ep_i[i+10].pos[1], ep_i[i+11].pos[1],
                  ep_i[i+12].pos[1], ep_i[i+13].pos[1], ep_i[i+14].pos[1], ep_i[i+15].pos[1]);
        v16sf pzi(ep_i[i+ 0].pos[2], ep_i[i+ 1].pos[2], ep_i[i+ 2].pos[2], ep_i[i+ 3].pos[2],
                  ep_i[i+ 4].pos[2], ep_i[i+ 5].pos[2], ep_i[i+ 6].pos[2], ep_i[i+ 7].pos[2],
                  ep_i[i+ 8].pos[2], ep_i[i+ 9].pos[2], ep_i[i+10].pos[2], ep_i[i+11].pos[2],
                  ep_i[i+12].pos[2], ep_i[i+13].pos[2], ep_i[i+14].pos[2], ep_i[i+15].pos[2]);
        v16sf axi(0.), ayi(0.), azi(0.), pti(0.);

        static __thread PS::F32 tjp[nvector];
        tjp[0] = ep_j[0].pos[0];
        tjp[1] = ep_j[0].pos[1];
        tjp[2] = ep_j[0].pos[2];
        tjp[3] = ep_j[0].mass;
        v16sf jp;
        jp.load(tjp);
        for(PS::S32 j = 0; j < n_jp; j++) {
            jp = v16sf::permute(jp, (_MM_PERM_ENUM)0b00000000);
            v16sf dx  = pxi - v16sf::swizzle(jp, _MM_SWIZ_REG_AAAA);
            v16sf dy  = pyi - v16sf::swizzle(jp, _MM_SWIZ_REG_BBBB);
            v16sf dz  = pzi - v16sf::swizzle(jp, _MM_SWIZ_REG_CCCC);            
            v16sf mj  = v16sf::swizzle(jp, _MM_SWIZ_REG_DDDD);

            tjp[0] = ep_j[j+1].pos[0];
            tjp[1] = ep_j[j+1].pos[1];
            tjp[2] = ep_j[j+1].pos[2];
            tjp[3] = ep_j[j+1].mass;
            jp.load(tjp);

            v16sf r2  = ep2 + dx * dx + dy * dy + dz * dz;
            v16sf ri  = v16sf::rsqrt_0th(r2);
            v16sf mri = ri * mj;
            v16sf mr3 = mri * ri * ri;
            pti -= mri;
            axi -= mr3 * dx;
            ayi -= mr3 * dy;
            azi -= mr3 * dz;                        
        }

        PS::F32 buf_pt[nvector], buf_ax[nvector], buf_ay[nvector], buf_az[nvector];
        pti.store(buf_pt);
        axi.store(buf_ax);
        ayi.store(buf_ay);
        azi.store(buf_az);
        const PS::S32 nii = std::min((n_ip - i), nvector);
        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].pot    += buf_pt[ii];
            force[i+ii].acc[0] += buf_ax[ii];
            force[i+ii].acc[1] += buf_ay[ii];
            force[i+ii].acc[2] += buf_az[ii];
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

    PS::F64 ts = getWallclockTime();

    const PS::S32 njmax   = 65536;
    const PS::S32 nvector = 16;
    const PS::S32 nipar   = 4;
    const PS::S32 njpar   = 4;
    v16sf * jptcl;
    PS::S32 ret = posix_memalign((void **)&jptcl, 1024, sizeof(v16sf)*n_jp);
    //PS::S32 ret = posix_memalign((void **)&jptcl, 1024, sizeof(v16sf)*njmax);
    
    //assert(n_jp <= njmax);
    for(PS::S32 j = 0, jarray = 0; j < n_jp; j += njpar, jarray++) {
        static __thread PS::F32 jp[nvector];
        PS::S32 njj = std::min(n_jp - j, njpar);
        for(PS::S32 jj = 0; jj < njj; jj++) {
            jp[4*jj+0] = ep_j[j+jj].pos[0];
            jp[4*jj+1] = ep_j[j+jj].pos[1];
            jp[4*jj+2] = ep_j[j+jj].pos[2];
            jp[4*jj+3] = ep_j[j+jj].mass;
        }
        for(PS::S32 jj = njj; jj < njpar; jj++) {
            jp[4*jj+0] = 0.;
            jp[4*jj+1] = 0.;
            jp[4*jj+2] = 0.;
            jp[4*jj+3] = 0.;
        }
        jptcl[jarray].load(jp);
        //jptcl[jarray] = _mm512_load_ps(jp);
    }

    v16sf ep2(FPGrav::eps * FPGrav::eps);
    for(PS::S32 i = 0; i < n_ip; i += nipar){

        v16sf pxi(ep_i[i+0].pos[0], ep_i[i+1].pos[0], ep_i[i+2].pos[0], ep_i[i+3].pos[0],
                  ep_i[i+0].pos[0], ep_i[i+1].pos[0], ep_i[i+2].pos[0], ep_i[i+3].pos[0],
                  ep_i[i+0].pos[0], ep_i[i+1].pos[0], ep_i[i+2].pos[0], ep_i[i+3].pos[0],
                  ep_i[i+0].pos[0], ep_i[i+1].pos[0], ep_i[i+2].pos[0], ep_i[i+3].pos[0]);
        v16sf pyi(ep_i[i+0].pos[1], ep_i[i+1].pos[1], ep_i[i+2].pos[1], ep_i[i+3].pos[1],
                  ep_i[i+0].pos[1], ep_i[i+1].pos[1], ep_i[i+2].pos[1], ep_i[i+3].pos[1],
                  ep_i[i+0].pos[1], ep_i[i+1].pos[1], ep_i[i+2].pos[1], ep_i[i+3].pos[1],
                  ep_i[i+0].pos[1], ep_i[i+1].pos[1], ep_i[i+2].pos[1], ep_i[i+3].pos[1]);
        v16sf pzi(ep_i[i+0].pos[2], ep_i[i+1].pos[2], ep_i[i+2].pos[2], ep_i[i+3].pos[2],
                  ep_i[i+0].pos[2], ep_i[i+1].pos[2], ep_i[i+2].pos[2], ep_i[i+3].pos[2],
                  ep_i[i+0].pos[2], ep_i[i+1].pos[2], ep_i[i+2].pos[2], ep_i[i+3].pos[2],
                  ep_i[i+0].pos[2], ep_i[i+1].pos[2], ep_i[i+2].pos[2], ep_i[i+3].pos[2]);

        v16sf axi(0.), ayi(0.), azi(0.), pti(0.);

        PS::S32 jarray = 0;
        //v16sf jp = _mm512_load_ps(&jptcl[jarray]);
        v16sf jp = jptcl[jarray];
        jarray++;
        for(PS::S32 j = 0; j < n_jp; j += njpar) {
            v16sf dx, dy, dz, mj;
            v16sf r2, ri, mri, mr3;

            dx  = pxi - v16sf::swizzle(jp, _MM_SWIZ_REG_AAAA);
            dy  = pyi - v16sf::swizzle(jp, _MM_SWIZ_REG_BBBB);
            dz  = pzi - v16sf::swizzle(jp, _MM_SWIZ_REG_CCCC);
            mj  = v16sf::swizzle(jp, _MM_SWIZ_REG_DDDD);
            //jp  = _mm512_load_ps(&jptcl[jarray]);
            jp  = jptcl[jarray];
            jarray++;

            r2  = ep2 + dx * dx + dy * dy + dz * dz;
            ri  = v16sf::rsqrt_0th(r2);
            mri = ri * mj;
            mr3 = mri * ri * ri;
            pti -= mri;
            axi -= mr3 * dx;
            ayi -= mr3 * dy;
            azi -= mr3 * dz;                                    
        }

        const __mmask16 imm[nipar] = {0x1111, 0x2222, 0x4444, 0x8888};
        const PS::S32 nii = std::min((n_ip - i), nipar);
        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].pot    += v16sf::mask_reduce_add(imm[ii], pti);
            force[i+ii].acc[0] += v16sf::mask_reduce_add(imm[ii], axi);
            force[i+ii].acc[1] += v16sf::mask_reduce_add(imm[ii], ayi);
            force[i+ii].acc[2] += v16sf::mask_reduce_add(imm[ii], azi);
        }

    }

    free(jptcl);

    wt_kernel[PS::Comm::getThreadNum()] += (getWallclockTime() - ts);

}

#endif

#else

template <class TParticleJ>
void CalcGravity(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
#pragma simd
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
