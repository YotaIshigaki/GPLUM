#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
/* CELib header */
#include "CELib.h"

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

//** Force class for gravity calculation
class Force_grav {
public:
    PS::F64vec acc;
    PS::F64 pot; 
    void clear();
};

//** Force classes for SPH calculation
class Force_dens{
public:
    PS::S32 flag;
    PS::F64 dens;
    PS::F64 smth;
    PS::F64 gradh;
    PS::F64 divv;
    PS::F64vec rotv;
    void clear();
};
class Force_hydro{
public:
    PS::F64vec acc;
    PS::F64 eng_dot;
    PS::F64 ent_dot;
    PS::F64 dt;
    void clear();
};

//** Full Particle Classes
class FP_nbody{
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
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

class FP_star{
public:
    PS::S32    flag;
    PS::S64    pid; // parent id
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
    void copyFromForce(const Force_dens & f);
    void writeAscii(FILE* fp) const;
    void readAscii(FILE* fp);
    void writeRestartData(FILE* fp) const;
    void readRestartData(FILE* fp);
    /* Other member functions */
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
    PS::F64    ent; // entropy
    PS::F64    pres; // pressure
    PS::F64    smth; // smoothing length
    PS::F64    gradh; // grad-h term
    PS::F64    divv; // divergence of velocity
    PS::F64vec rotv; // rotation of velocity
    PS::F64    BalSW; // Balsara switch
    PS::F64    snds; // sound speed
    PS::F64    eng_dot; // time rate of change of `eng`
    PS::F64    ent_dot; // time rate of change of `ent`
    PS::F64    dt; // hydrodynamic time step for this particle
    PS::F64vec vel_half;
    PS::F64    eng_half;
    PS::F64    ent_half;

    PS::F64    mabn[CELibYield_Number]; // mass abundance
    PS::S32    n_stars; // # of stars spawned from this particle

    FP_gas(); 

    /* Member functions required by FDPS */
    PS::S64 getId() const;
    PS::F64 getCharge() const;
    PS::F64vec getPos() const;
    PS::F64 getRSearch() const;
    void setPos(const PS::F64vec& pos);
    void copyFromForce(const Force_grav& f);
    void copyFromForce(const Force_dens& f);
    void copyFromForce(const Force_hydro& f);
    void writeAscii(FILE* fp) const;
    void readAscii(FILE* fp);
    void writeRestartData(FILE* fp) const;
    void readRestartData(FILE* fp);
    /* Other member functions */
    PS::F64 getMass() const;
    PS::F64 getKernelSupportRadius() const;
    PS::F64 getMassDensity() const;
    PS::F64 getTemperature() const;
    PS::F64 getInternalEnergyDensity() const;
    PS::F64 getMetallicity() const;
    void setEntropy();
    void setPressure();
    void setPressureFromSpecificInternalEnergy(const PS::F64 eng_new);
    void setPressureFromInternalEnergyDensity(const PS::F64 U);
    void dump(const std::string msg_id,
              const PS::S32 arr_idx,
              const std::string caller_name,
              std::ostream & fout=std::cout) const;
};

//** Essential Particle Class
class EP_grav {
public:
    static PS::F64 eps;

    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;

    /* Member functions required by FDPS */
    PS::S64 getId() const;
    PS::F64 getCharge() const;
    PS::F64vec getPos() const;
    void copyFromFP(const FP_nbody& fp);
    void copyFromFP(const FP_star& fp);
    void copyFromFP(const FP_gas& fp);
};

class EP_hydro {
public:
    PS::S32    type; // particle type
    PS::S32    rank; // MPI rank
    PS::S32    idx;  // array index in ParticleSystem
    PS::S64    id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    mass;
    PS::F64    smth;
    PS::F64    dens;
    PS::F64    pres;
    PS::F64    gradh;
    PS::F64    snds;
    PS::F64    BalSW;

    PS::S64 getId() const;
    PS::F64vec getPos() const;
    PS::F64 getRSearch() const;
    void setPos(const PS::F64vec& pos);
    void copyFromFP(const FP_gas& fp);
    void copyFromFP(const FP_star& fp);
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
        const PS::F64 eps2 = SQ(EP_grav::eps);
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
};
#endif

class CalcDensity {
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
