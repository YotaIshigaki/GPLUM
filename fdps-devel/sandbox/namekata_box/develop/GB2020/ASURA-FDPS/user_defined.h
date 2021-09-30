#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
/* CELib header */
#include "CELib.h"

enum class ParticleType {Gas, Star, DarkMatter};
enum class FunctionUseType {Main, PostProcess};

/* Definitions of user-defined classes */
//** File header class (used for IO)
class FileHeader {
public:
    PS::F64 time;
    PS::S64 numPtcl;
    PS::S32 readAscii(FILE* fp);
    void writeAscii(FILE* fp) const;
    PS::S32 readBinary(FILE* fp);
    void writeBinary(FILE* fp) const;
};

class TestHeader {
public:
    PS::F64 time;
    PS::S32 readAscii(FILE*fp);
    void writeAscii(FILE* fp) const;
};

//** Moment class
class Mom_grav {
public:
    PS::F64 mass;
    PS::F64 eps2;
    PS::F64vec pos;
    Mom_grav() {
        this->mass = 0.0;
        this->eps2 = 0.0;
        this->pos = 0.0;
    }
    Mom_grav(const PS::F64 mass,
             const PS::F64 eps2,
             const PS::F64vec & pos){
        this->mass = mass;
        this->eps2 = eps2;
        this->pos = pos;
    }
    void init() {
        this->mass = 0.0;
        this->eps2 = 0.0;
        this->pos = 0.0;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getEps2() const {
        return this->eps2;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        this->mass += epj.mass;
        this->eps2 += epj.mass * (epj.eps * epj.eps);
        this->pos += epj.mass * epj.pos;
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void set(){
        this->eps2 /= this->mass;
        this->pos /= this->mass;
    }
    void accumulate(const Mom_grav & mom){
        this->mass += mom.mass;
        this->eps2 += mom.mass * mom.eps2;
        this->pos += mom.mass * mom.pos;
    }
    void accumulate2(const Mom_grav & mom){}
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

//** SPJ class
class SPJ_grav {
public:
    PS::F64 mass;
    PS::F64 eps2;
    PS::F64vec pos;
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        this->mass = mom.mass;
        this->eps2 = mom.eps2;
        this->pos = mom.pos;
    }
    void clear(){
        this->mass = 0.0;
        this->eps2 = 0.0;
        this->pos = 0.0;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64 getEps2() const {
        return this->eps2;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F64vec & pos_new) {
        this->pos = pos_new;
    }
    Mom_grav convertToMoment() const {
        return Mom_grav(this->mass, this->eps2, this->pos);
    }
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const {
    fout << msg_id
         << " (SPJ_grav)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " mass = " << this->mass
         << " eps2 = " << this->eps2
         << " pos = " << this->pos
         << " @" << caller_name
         << std::endl;
    }
};

//** Force class for gravity calculation
class Force_grav {
public:
    PS::F64vec acc;
    PS::F64 pot; 
    void clear();
};

//** Force classes for SPH calculation
class Force_knl_sz {
public:
    PS::S32 flag;
    PS::F64 rad;
    void clear();
};


class Force_dens{
public:
    PS::F64 dens;
    PS::F64 pres;
    PS::F64 gradh;
    PS::F64 divv;
    PS::F64vec rotv;
    void clear();
};

class Force_hydro{
public:
    PS::F64vec acc;
    PS::F64 eng_dot;
    PS::F64 dt;
    void clear();
};

//** Full Particle Classes
class FP_dm{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;

    /* Member functions required by FDPS */
    PS::S64 getId() const;
    PS::F64vec getPos() const;
    PS::F64 getCharge() const;
    void setPos(const PS::F64vec& pos);
    void copyFromForce(const Force_grav & f);
    void writeAscii(FILE* fp) const;
    void readAscii(FILE* fp);
    void writeRestartData(FILE* fp) const;
    void readRestartData(FILE* fp);
    /* Other member functions */
    void clearGravitationalForce();
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

class FP_star{
private:
    static const PS::S32 nbit = 2; // the number of bits used to represent cid
public:
    PS::S32    flag;
    PS::S64    id;
    PS::F64    mass0; // initial mass
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;

    PS::F64    mabn[CELibYield_Number]; // mass abundance

    PS::F64    t_form; // formation time
    PS::U64    FBcnt;  // feedback counter (we assign 16 bits for each FB type)
    PS::F64    t_SNII; // event time of SNe II
    PS::F64    t_SNIa; // event time of SNe Ia;
    PS::F64    t_AGB;  // event time of AGB
    PS::F64    t_NSM;  // event time of neutron star merger
    PS::F64    dens;   // actually the sum of weights
    PS::F64    FBrad;  // feedback radius
   
    FP_star();

    /* Member functions required by FDPS */
    PS::S64 getId() const;
    PS::F64vec getPos() const;
    PS::F64 getCharge() const;
    PS::F64 getRSearch() const;
    void setPos(const PS::F64vec& pos);
    void copyFromForce(const Force_grav & f);
    void copyFromForce(const Force_knl_sz & f);
    void copyFromForce(const Force_dens & f);
    void writeAscii(FILE* fp) const;
    void readAscii(FILE* fp);
    void writeRestartData(FILE* fp) const;
    void readRestartData(FILE* fp);
    /* Other member functions */
    void setId(const PS::S64 cid, const PS::S64 pid);
    PS::S64 getPid() const;
    PS::S64 getCid() const;
    void copyAbundanceFrom(const PS::F64 mabn[]);
    PS::F64 getMetallicity() const;
    PS::U32 getSNIICount() const;
    void setSNIICount(const PS::U64 cnt);
    PS::U32 getSNIaCount() const;
    void setSNIaCount(const PS::U64 cnt);
    PS::U32 getAGBCount() const;
    void setAGBCount(const PS::U64 cnt);
    PS::U32 getNSMCount() const;
    void setNSMCount(const PS::U64 cnt);
    bool feedbackAsSNII(const PS::F64 t, const PS::F64 dt) const;
    bool feedbackAsSNIa(const PS::F64 t, const PS::F64 dt) const;
    bool feedbackAsAGB(const PS::F64 t, const PS::F64 dt) const;
    bool feedbackAsNSM(const PS::F64 t, const PS::F64 dt) const;
    void clearGravitationalForce();
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

class FP_gas {
public:
    PS::S32    idx; // array index in ParticleSystem
    PS::S64    id;
    PS::F64    mass0; // initial mass
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc_grav; // gravitational acceleration
    PS::F64    pot_grav; // gravitational potential
    PS::F64vec acc_hydro; // acceleration due to pressure-gradient
    PS::S32    flag;
    PS::F64    dens; // mass density
    PS::F64    eng; // specific internal energy
    PS::F64    pres; // pressure
    PS::F64    h; // smoothing length
    PS::F64    h_prev; // smoothing length at the previous step
    PS::F64    h_dot_prev; // time rate of change of smoothing length at the previous step
    PS::F64    gradh; // grad-h term
    PS::F64    divv; // divergence of velocity
    PS::F64vec rotv; // rotation of velocity
    PS::F64    BalSW; // Balsara switch
    PS::F64    alpha; // \alpha parameter in artificial viscosity
    PS::F64    snds; // sound speed
    PS::F64    eng_dot; // time rate of change of `eng`
    PS::F64    dt; // hydrodynamic time step for this particle
    PS::F64vec vel_half;
    PS::F64    eng_half;

    PS::F64    mabn[CELibYield_Number]; // mass abundance
    PS::S32    n_stars; // # of stars spawned from this particle

    FP_gas(); 

    /* Member functions required by FDPS */
    PS::S64 getId() const;
    PS::F64 getCharge() const;
    PS::F64vec getPos() const;
    PS::F64 getRSearch() const;
    void setPos(const PS::F64vec& pos);
    void copyFromForce(const Force_grav & f);
    void copyFromForce(const Force_knl_sz & f);
    void copyFromForce(const Force_dens & f);
    void copyFromForce(const Force_hydro & f);
    void writeAscii(FILE* fp) const;
    void readAscii(FILE* fp);
    void writeRestartData(FILE* fp) const;
    void readRestartData(FILE* fp);
    void writeGlassData(FILE* fp) const;
    void readGlassData(FILE* fp);
    /* Other member functions */
    PS::F64 getMass() const;
    PS::F64 getKernelSupportRadius() const;
    PS::F64 getMassDensity() const;
    PS::F64 getTemperature() const;
    PS::F64 getInternalEnergyDensity() const;
    PS::F64 getMetallicity() const;
    void calcSoundSpeed();
    void applyTemperatureLimits();
    void clearGravitationalForce();
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

//** Essential Particle Class
class EP_grav {
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64 eps;

    /* Member functions required by FDPS */
    PS::S64 getId() const;
    PS::F64 getCharge() const;
    PS::F64 getEps2() const;
    PS::F64vec getPos() const;
    void copyFromFP(const FP_dm& fp);
    void copyFromFP(const FP_star& fp);
    void copyFromFP(const FP_gas& fp);
    void writeAscii(std::ostream & fout) const;
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};


class EPI_knl_sz {
public:
    ParticleType type;
    PS::S32 flag;
    PS::S64 id;
    PS::F64vec pos;
    PS::F64 rad;

    PS::S64 getId() const;
    PS::F64vec getPos() const;
    PS::F64 getRSearch() const;
    PS::F64 getRPhysical() const;
    void setPos(const PS::F64vec& pos);
    void copyFromFP(const FP_gas& fp);
    void copyFromFP(const FP_star& fp);
    void writeAscii(std::ostream & fout) const;
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

class EPJ_knl_sz {
public:
    ParticleType type;
    PS::F64vec pos;

    PS::F64vec getPos() const;
    void setPos(const PS::F64vec& pos);
    void copyFromFP(const FP_gas& fp);
    void copyFromFP(const FP_star& fp);
    void writeAscii(std::ostream & fout) const;
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

class EP_hydro {
public:
    ParticleType type; 
    PS::S32    rank; // MPI rank
    PS::S32    idx;  // array index in ParticleSystem
    PS::S64    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    mass;
    PS::F64    eng;
    PS::F64    h;
    PS::F64    dens;
    PS::F64    pres;
    PS::F64    gradh;
    PS::F64    snds;
    PS::F64    BalSW;
    PS::F64    alpha;

    PS::S64 getId() const;
    PS::F64vec getPos() const;
    PS::F64 getRSearch() const;
    PS::F64 getRPhysical() const;
    void setPos(const PS::F64vec& pos);
    void copyFromFP(const FP_gas& fp);
    void copyFromFP(const FP_star& fp);
    void writeAscii(std::ostream & fout) const;
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

/* Interaction functions */
#if defined(ENABLE_PHANTOM_GRAPE_X86)
template <class TParticleJ>
class CalcGravity {
public:
    void operator () (const EP_grav * iptcl,
                      const PS::S32 ni,
                      const TParticleJ * jptcl,
                      const PS::S32 nj,
                      Force_grav * force) {
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
        }
        PS::S32 devid = PS::Comm::getThreadNum();
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
#else
template <class TParticleJ>
class CalcGravity {
public:
    void operator () (const EP_grav * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * ep_j,
                      const PS::S32 n_jp,
                      Force_grav * force) {
        for (PS::S32 i = 0; i < n_ip; i++){
            const PS::F64vec xi = ep_i[i].getPos();
            const PS::F64 eps2i = ep_i[i].getEps2();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j = 0; j < n_jp; j++){
                const PS::F64vec rij = xi - ep_j[j].getPos();
                const PS::F64 eps2j = ep_j[j].getEps2();
                PS::F64 r3_inv = rij * rij + eps2i + eps2j;
                PS::F64 r_inv  = 1.0/std::sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].getCharge();
                r3_inv *= r_inv;
                ai   -= r3_inv * rij;
                poti -= r_inv;

#if defined(ENABLE_NAN_CHECK)
                //if (ai.isnan() ||
                //    ai.isinf() ||
                //    std::isnan(poti) ||
                //    std::isinf(poti)) {
                //    ep_i[i].dump("[nan(i)]", i, __func__, dbg_utils::fout);
                //    ep_j[j].dump("[nan(j)]", j, __func__, dbg_utils::fout);
                //    assert(false);
                //}
#endif
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#endif

extern bool output_flag;
extern PS::F64 et_rij_calc;
extern PS::F64 et_bisec;

PS::F64 getPredictedKernelSize(const PS::S32 N_ngb_cur,
                               const PS::F64 h_cur);

class CalcKernelSize {
public:
    void operator () (const EPI_knl_sz * ep_i,
                      const PS::S32 n_ip,
                      const EPJ_knl_sz * ep_j,
                      const PS::S32 n_jp,
                      Force_knl_sz * force);
};


class CalcDensityAndPressure {
public:
    void operator () (const EP_hydro * ep_i,
                      const PS::S32 n_ip,
                      const EP_hydro * ep_j,
                      const PS::S32 n_jp,
                      Force_dens * force);
};

class CalcHydroForce {
public:
    void operator () (const EP_hydro * ep_i,
                      const PS::S32 n_ip,
                      const EP_hydro * ep_j,
                      const PS::S32 n_jp,
                      Force_hydro * force);
};
