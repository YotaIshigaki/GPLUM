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

class Planet{
public:
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64vec vel;
    PS::S64 id;
    static PS::F64 r_coll;
    // [Notes]
    // (1) The order of members must be consistent with that of
    //     Long's force kernel, namely, EpiMM type.
    // (2) In the original ring code, r_coll was a normal member.
    //     But, it is now changed to a static member because
    //     Long's force kernel assumes that satellite type consists
    //     of pos,mass,vel,id only.
    Planet(){
        pos = 0.0;
        mass = 0.0;
        vel = 0.0;
    } 
    Planet(const PS::F64 _mass, const PS::F64vec & _pos, 
           const PS::F64vec & _vel): 
        mass(_mass), pos(_pos), vel(_vel){}
    void clear(){}
};
typedef Planet Satellite; // Planet::r_coll is identical to Satellite::r_coll

class SatelliteForce{
public:
    PS::F64vec acc;
    PS::F64    pot;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

//////////////// for satellite particles /////////////
class NgbSat{
public:
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    NgbSat():mass(0.0), pos(0.0), vel(0.0){}
    void clear(){
        pos = vel = 0.0;
        mass = 0.0;
    }
};

class ForceSat{
public:
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64 n_coll;
    ForceSat():acc(0.0), pot(0.0), n_coll(0.0){}
    void clear(){
        acc = 0.0;
        pot = n_coll = 0.0;
    }
};

class RNgbSqRank{
public:
    double r_ngb_sq;
    int rank;
};

//////////////// for satellite particles /////////////

#ifdef SUNWAY_FORCE_KERNEL
class Force{
public:
    void clear() {}
};

class FPGrav{
public:
    PS::F64vec pos;
    PS::F64    mass;
    PS::F64vec vel;
    PS::S64    id;
    //PS::F64vec acc;
    //PS::F64    pot;
    //PS::F64 r_init;
    PS::F64vec getPos() const {
        return pos;
    }
    void setPos(const PS::F64vec & pos_new){
        pos = pos_new;
    }    
    void copyFromForce(const Force & force) { }    
    /*
    void copyFromForce(const Force & force) {
        acc    = force.acc;
        pot    = force.pot;
    }
    */
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

#else
class Force{
public:
    PS::S64    id;
    PS::F64vec acc;
    PS::F64    pot;
    //Ngb        ngb;
    PS::F64 r_ngb_sq;
    PS::S64    n_coll;
    PS::S32 adr_ngb; // point to EPJ_POINTER in defined in force_sunway.hpp    
    Force(): id(-1), acc(0.0), pot(0.0), n_coll(0), r_ngb_sq(PS::LARGE_FLOAT), adr_ngb(-1){}
    void clear() {
        acc = 0.0;
        pot = 0.0;
        n_coll = 0;
        //ngb.clear();
        r_ngb_sq = PS::LARGE_FLOAT;
        adr_ngb = -1;        
    }
};

class FPGrav{
public:
    PS::S64    id;
    //PS::S64    id_ngb;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;
    //PS::S32    rank_org;
    //PS::S32    n_coll;
    //PS::F64    r_ngb_sq;
    //static PS::F64 r_coll; // physical radius of particle
    PS::F64vec getPos() const {
        return pos;
    }
    void setPos(const PS::F64vec & pos_new){
        pos = pos_new;
    }
    void copyFromForce(const Force & force) {
        //id_ngb = force.ngb.id;
        acc    = force.acc;
        pot    = force.pot;
        //n_coll = force.n_coll;
        //r_ngb_sq = force.r_ngb_sq;
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

#endif

class Epi{
public:
    PS::F64vec pos;
    PS::F64    mass;
    PS::F64vec vel;
    PS::S64    id;
    static PS::F64 eps;
    PS::F64vec getPos() const {
        return pos;
    }
    void copyFromFP(const FPGrav & fp){ 
        pos  = fp.pos;
        mass = fp.mass;
        vel  = fp.vel;
        id   = fp.id;
    }
};

class Epj{
public:
    PS::F64vec pos;
    PS::F64    mass;
    PS::F64vec vel;
    PS::S64    id;
#ifdef PHI_R_TREE
    //PS::F64    pos_phi;
    //PS::F64    pos_r;
#endif
    static PS::F64 r_search;
    static PS::F64 r_coll; // physical radius of particle
    PS::F64vec getPos() const {
        return pos;
    }
#ifdef PHI_R_TREE
    // for debug, actually not used
    PS::F64vec getPosForTree() const {
        const PS::F64 MY_PI = 4.0*atan(1.0);
        PS::F64 pos_phi_tmp = atan2(pos.y, pos.x);
        if(pos_phi_tmp < 0.0){
            pos_phi_tmp += 2.0*MY_PI;
        }
        assert(pos_phi_tmp >= 0.0 && pos_phi_tmp < 2.0*MY_PI);
        PS::F64 pos_r_tmp   = sqrt(pos.x*pos.x + pos.y*pos.y);
        return PS::F64vec(pos_phi_tmp, pos_r_tmp, pos.z);
    }
#endif
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64 getRSearch() const {
        return r_search;
    }
    template<class Tepi>
    void copyFromEPI(const Tepi & epi){
        pos = epi.pos;
        mass = epi.mass;
        id = epi.id;
    }
    void copyFromFP(const FPGrav & fp){ 
        pos      = fp.pos;
        mass     = fp.mass;
        vel      = fp.vel;
        id       = fp.id;
        //rank_org = fp.rank_org;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass= "<<mass<<std::endl;
        fout<<"pos= "<<pos<<std::endl;
#ifdef PHI_R_TREE
        //fout<<"pos_phi= "<<pos_phi<<" pos_r= "<<pos_r<<std::endl;
#endif
    }
};

//* Moment classes for Sunway Taihulight application
class MomentMonopoleScatterSW{
public:
    PS::F32 mass;
    PS::F32vec pos;
    PS::F64ort vertex_out_;
    PS::F64ort vertex_in_;
#ifdef PHI_R_TREE
    PS::F32 pos_phi;
    PS::F32 pos_r;
#endif
    MomentMonopoleScatterSW(){
        mass = 0.0;
        pos = 0.0;
        vertex_out_.init();
        vertex_in_.init();
#ifdef PHI_R_TREE
        pos_phi = pos_r = 0.0;
#endif
    }
#ifdef PHI_R_TREE
    MomentMonopoleScatterSW(const PS::F32 m, const PS::F32vec & p,
                            const PS::F32 p_phi, const PS::F32 p_r){
        mass = m;
        pos = p;
        pos_phi = p_phi;
        pos_r = p_r;
    }
#else
    MomentMonopoleScatterSW(const PS::F32 m, const PS::F32vec & p){ 
        mass = m;
        pos = p;
    }
#endif
    PS::F64ort getVertexOut() const { return vertex_out_; }
    PS::F64ort getVertexIn() const { return vertex_in_; }
    void init(){
        mass = 0.0;
        pos = 0.0;
        vertex_out_.init();
        vertex_in_.init();
#ifdef PHI_R_TREE
        pos_phi = pos_r = 0.0;
#endif
    }
    PS::F32vec getPos() const {
        return pos;
    }
#ifdef PHI_R_TREE
    PS::F32vec getPosForTree() const {
        return PS::F64vec(pos_phi, pos_r, pos.z);
    }
#endif
    PS::F32 getCharge() const {
        return mass;
    }
    
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        // for debug
        this->mass += epj.getCharge();
        this->pos += epj.getCharge() * epj.getPos();
#ifdef PHI_R_TREE
        (this->vertex_out_).merge(epj.getPosForTree(), epj.getRSearch());
        (this->vertex_in_).merge(epj.getPosForTree());
        this->pos_phi += epj.getCharge() * epj.getPosForTree().x;
        this->pos_r   += epj.getCharge() * epj.getPosForTree().y;
#else
        (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        (this->vertex_in_).merge(epj.getPos());
#endif
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void set(){
        pos /=  mass;
#ifdef PHI_R_TREE
        pos_phi /= mass;
        pos_r   /= mass;
#endif
    }
    void accumulate(const MomentMonopoleScatterSW & mom){
        this->mass += mom.mass;
        this->pos += mom.mass * mom.pos;
        (this->vertex_out_).merge(mom.vertex_out_);
        (this->vertex_in_).merge(mom.vertex_in_);        
#ifdef PHI_R_TREE
        this->pos_phi += mom.mass * mom.pos_phi;
        this->pos_r   += mom.mass * mom.pos_r;
#endif
    }
    void accumulate2(const MomentMonopoleScatterSW & mom){}
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
        fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
        fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
        fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
    }
};

//* SuperParticleJ classes for Sunway Taihulight application
class SPJMonopoleScatterSW{
public:
    PS::F32vec pos;
    PS::F32 mass;
#ifdef PHI_R_TREE
    PS::F32 pos_phi;
    PS::F32 pos_r;
    //PS::F32 pad[2];
#endif
    // The order of members must be consistent with that of
    // Long's force kernel.
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        PS::F32 mass_tmp = mom.mass;
        PS::F32vec pos_tmp = mom.pos;
        this->mass = mass_tmp;
        this->pos = pos_tmp;
#ifdef PHI_R_TREE
        PS::F32 pos_phi_tmp = mom.pos_phi;
        PS::F32 pos_r_tmp   = mom.pos_r;
        this->pos_phi = pos_phi_tmp;
        this->pos_r   = pos_r_tmp;
#endif
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
#ifdef PHI_R_TREE
        pos_phi = pos_r = 0.0;
#endif
    }
    PS::F32 getCharge() const {
        return mass;
    }
    PS::F32vec getPos() const {
        return pos;
    }
#ifdef PHI_R_TREE
    PS::F32vec getPosForTree() const {
        return PS::F32vec(pos_phi, pos_r, pos.z);
    }
#endif
    MomentMonopoleScatterSW convertToMoment() const {
#ifdef PHI_R_TREE
        return MomentMonopoleScatterSW(mass, pos, pos_phi, pos_r);
#else
        return MomentMonopoleScatterSW(mass, pos);
#endif
    }
    void dump(std::ostream & fout=std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
#ifdef PHI_R_TREE
        fout<<"pos_phi= "<<pos_phi<<" pos_r= "<<pos_r<<std::endl;
#endif
    }
};

typedef Epi EPI;
typedef Epj EPJ;
typedef Force FORCE;
//typedef PS::TreeForForceLong<FORCE, EPI, EPJ>::MonopoleWithScatterSearch TreeType;
typedef PS::TreeForForce
<PS::SEARCH_MODE_LONG_SCATTER,
 FORCE, EPI, EPJ,
 MomentMonopoleScatterSW,
 MomentMonopoleScatterSW,
 SPJMonopoleScatterSW> TreeType;

