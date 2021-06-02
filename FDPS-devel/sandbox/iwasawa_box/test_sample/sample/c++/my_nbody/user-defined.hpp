#pragma once
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

class ForceGrav0{
public:
    PS::F64vec acc;
    PS::F64    pot;
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
};

class FPGrav0{
public:
    //PS::S64    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    mass;
    PS::F64    pot;    
    static PS::F64 eps;
    static PS::F64 dt_tree;
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    void copyFromFP(const FPGrav0 & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
        //id   = fp.id;
    }
    void copyFromForce(const ForceGrav0 & force) {
        acc = force.acc;
        pot = force.pot;
    }
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
    void writeAscii(FILE* fp) const {
        /*
        fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
        */
    }
    void readAscii(FILE* fp) {
        /*
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
        */
    }
    void integrate(){
        vel += 0.5*dt_tree*acc;
        pos += dt_tree*vel;
    }
};

class EPIGrav0{
public:
    PS::F64vec pos;
    //PS::F64    mass;
    static PS::F64 eps;
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {}
    /*
    PS::F64 getCharge() const {
        return mass;
    }
    */
    void copyFromFP(const FPGrav0 & fp){ 
        //mass = fp.mass;
        pos  = fp.pos;
    }
};

class EPJGrav0{
public:
    PS::F64vec pos;
    PS::F64    mass;
    static PS::F64 eps;
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    void copyFromFP(const FPGrav0 & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
    }
};



/////////////////////////
// OPTIMIZED GPU VERSION
class ForceGrav{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    void clear() {}
    /*
    PS::F64vec acc;
    PS::F64    pot;    
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
    */
};

class FPGrav{
public:
    //PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    //PS::F64vec acc;
    //PS::F64    pot;    
    static PS::F64 eps;

    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromForce(const ForceGrav & force) {
        pos = force.pos;
        vel = force.vel;
        //acc = force.acc;
        //pot = force.pot;
    }

    /*
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
    */

    void writeAscii(FILE* fp) const {
        /*
        fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
        */
    }


    void readAscii(FILE* fp) {
        /*
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
        */
    }

};

class EPIGrav{
public:
    PS::F64vec pos;
    //PS::F64    mass;
    PS::F64vec vel;
    static PS::F64 eps;

    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {}
    /*
    PS::F64 getCharge() const {
        return mass;
    }
    */
    void copyFromFP(const FPGrav & fp){ 
        //mass = fp.mass;
        pos  = fp.pos;
        vel = fp.vel;
    }
};

class EPJGrav{
public:
    PS::F64vec pos;
    PS::F64    mass;
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
};


#ifdef ENABLE_PHANTOM_GRAPE_X86


template <class TParticleI, class TParticleJ, class Tforce>
void CalcGravity(const TParticleI * iptcl,
                 const PS::S32 ni,
                 const TParticleJ * jptcl,
                 const PS::S32 nj,
                 Tforce * force) {
    const PS::S32 nipipe = ni;
    const PS::S32 njpipe = nj;
    PS::F64 (*xi)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64 (*ai)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64  *pi     = (PS::F64  *    )malloc(sizeof(PS::F64) * nipipe);
    PS::F64 (*xj)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * njpipe * PS::DIMENSION);
    PS::F64  *mj     = (PS::F64  *    )malloc(sizeof(PS::F64) * njpipe);
    for(PS::S32 i = 0; i < ni; i++) {
        xi[i][0] = iptcl[i].getPos().x;
        xi[i][1] = iptcl[i].getPos().y;
        xi[i][2] = iptcl[i].getPos().z;
        ai[i][0] = 0.0;
        ai[i][1] = 0.0;
        ai[i][2] = 0.0;
        pi[i]    = 0.0;
    }
    for(PS::S32 j = 0; j < nj; j++) {
        xj[j][0] = jptcl[j].getPos().x;
        xj[j][1] = jptcl[j].getPos().y;
        xj[j][2] = jptcl[j].getPos().z;
        mj[j]    = jptcl[j].getCharge();
        xj[j][0] = jptcl[j].pos.x;
        xj[j][1] = jptcl[j].pos.y;
        xj[j][2] = jptcl[j].pos.z;
        mj[j]    = jptcl[j].mass;
    }
    PS::S32 devid = PS::Comm::getThreadNum();
    g5_set_xmjMC(devid, 0, nj, xj, mj);
    g5_set_nMC(devid, nj);
    g5_calculate_force_on_xMC(devid, xi, ai, pi, ni);
    for(PS::S32 i = 0; i < ni; i++) {
        force[i].acc.x += ai[i][0];
        force[i].acc.y += ai[i][1];
        force[i].acc.z += ai[i][2];
        force[i].pot   -= pi[i];
    }
    free(xi);
    free(ai);
    free(pi);
    free(xj);
    free(mj);
}

#else

template <class TParticleI, class TParticleJ, class Tforce>
void CalcGravity(const TParticleI * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 Tforce * force) {
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        //if(ep_i[i].id == 0) std::cerr<<"xi= "<<xi<<std::endl;
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




#ifdef OPTIMIZED_VERSION
    typedef EPIGrav    EPI;
    typedef EPJGrav    EPJ;
    typedef FPGrav     FP;
    typedef ForceGrav  Force;
#else
    typedef EPIGrav0    EPI;
    typedef EPJGrav0    EPJ;
    typedef FPGrav0     FP;
    typedef ForceGrav0  Force;
#endif
