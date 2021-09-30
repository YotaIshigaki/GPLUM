#pragma once
#include <particle_simulator.hpp>
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
#include "run_param.hpp"
#include "cosmology.hpp"

#define TINY (1.0e-30)

class Force_grav {
public:
    PS::F64vec acc;
    PS::F64    pot;
    PS::F64vec acc_pm;
    PS::F64    pot_pm;

    void clear() {
        this->acc = 0.0;
        this->pot = 0.0;
    }

    void clearPMM() {
        this->acc_pm = 0.0;
        this->pot_pm = 0.0;
    }

    void accumulateForcePMM(const PS::F64vec & acc, const PS::F64 & pot) {
        this->acc_pm += acc;
        this->pot_pm += pot;
    }
};


class FP_grav {
private:
    template<class T>
    T reverseEndian(T value){
        char * first = reinterpret_cast<char*>(&value);
        char * last = first + sizeof(T);
        std::reverse(first, last);
        return value;
    }
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec acc_pm;
    PS::F64 pot_pm;
    
    static PS::F64 eps;
    static PS::F64 H0;
    static PS::F64 Lbnd;
   
    PS::S64 getId() const {
        return this->id;
    } 
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F64vec & pos) {
        this->pos = pos;
    }
    void copyFromForce(const Force_grav & f) {
        this->acc = f.acc;
        this->pot = f.pot;
    }
    void copyFromForcePMM(const Force_grav & f) {
        this->acc_pm = f.acc_pm;
        this->pot_pm = f.pot_pm;
    }
    void clearEwald() {
        this->acc = 0.0;
        this->pot = 0.0;
    }
    void accumulateForceEwald(const PS::F64vec & acc, const PS::F64 & pot) {
        this->acc += acc;
        this->pot += pot;
    }
    void calcTotalForce() {
        this->acc += this->acc_pm;
        this->pot -= this->pot_pm;
    }

    void writeParticleBinary(FILE *fp) {
        PS::F32 x = pos[0];
        PS::F32 y = pos[1];
        PS::F32 z = pos[2];
        PS::F32 vx = vel[0];
        PS::F32 vy = vel[1];
        PS::F32 vz = vel[2];    
        PS::S32 i = id;
        PS::F32 m = mass;
        fwrite(&x,  sizeof(PS::F32),1,fp);
        fwrite(&vx, sizeof(PS::F32),1,fp);
        fwrite(&y,  sizeof(PS::F32),1,fp);
        fwrite(&vy, sizeof(PS::F32),1,fp);
        fwrite(&z,  sizeof(PS::F32),1,fp);
        fwrite(&vz, sizeof(PS::F32),1,fp);    
        fwrite(&m,  sizeof(PS::F32),1,fp);
        fwrite(&i,  sizeof(PS::S32),1,fp);
    }

    // for API of FDPS 
    // in snapshot, L unit is Mpc/h, M unit is Msun, v unit is km/s
    void readBinary(FILE *fp){
        static PS::S32 ONE = 1;
        static bool is_little_endian = *reinterpret_cast<char*>(&ONE) == ONE;
        static const PS::F64 Mpc_m = 3.08567e22; // unit is m
        static const PS::F64 Mpc_km = 3.08567e19; // unit is km    
        static const PS::F64 Msun_kg = 1.9884e30; // unit is kg
        static const PS::F64 G = 6.67428e-11; // m^3*kg^-1*s^-2
        static const PS::F64 Cl = 1.0 / FP_grav::Lbnd;
        static const PS::F64 Cv = 1.0 / (FP_grav::Lbnd * FP_grav::H0);
        static const PS::F64 Cm = 1.0 / (pow(Mpc_m*FP_grav::Lbnd, 3.0) / pow(Mpc_km/FP_grav::H0, 2.0) / G / Msun_kg);
        PS::F32 x, y, z, vx, vy, vz, m;
        PS::S32 i;
        fread(&x,  4, 1, fp);
        fread(&vx, 4, 1, fp);
        fread(&y,  4, 1, fp);
        fread(&vy, 4, 1, fp);
        fread(&z,  4, 1, fp);
        fread(&vz, 4, 1, fp);
        fread(&m,  4, 1, fp);
        fread(&i,  4, 1, fp);
        if (is_little_endian) {
            pos.x = x * Cl;
            pos.y = y * Cl;
            pos.z = z * Cl;
            vel.x = vx * Cv;
            vel.y = vy * Cv;
            vel.z = vz * Cv;
            mass = m * Cm;
            id = i;
        }
        else{
            pos.x = reverseEndian(x) * Cl;
            pos.y = reverseEndian(y) * Cl;
            pos.z = reverseEndian(z) * Cl;
            vel.x = reverseEndian(vx) * Cv;
            vel.y = reverseEndian(vy) * Cv;
            vel.z = reverseEndian(vz) * Cv;
            mass = reverseEndian(m) * Cm;
            id = reverseEndian(i);
        }
    }

    // for API of FDPS 
    void writeBinary(FILE *fp){
        static const PS::F64 Mpc_m = 3.08567e22; // unit is m
        static const PS::F64 Mpc_km = 3.08567e19; // unit is km    
        static const PS::F64 Msun_kg = 1.9884e30; // unit is kg
        static const PS::F64 G = 6.67428e-11; // m^3*kg^-1*s^-2
        static const PS::F64 Cl = FP_grav::Lbnd;
        static const PS::F64 Cv = (FP_grav::Lbnd * FP_grav::H0);
        static const PS::F64 Cm = (pow(Mpc_m*FP_grav::Lbnd, 3.0) / pow(Mpc_km/FP_grav::H0, 2.0) / G / Msun_kg);    
        PS::F32vec x = pos * Cl;
        PS::F32vec v = vel * Cv;
        PS::F32 m = mass * Cm;
        PS::S32 i = id;
        fwrite(&x.x, sizeof(PS::F32), 1, fp);
        fwrite(&v.x, sizeof(PS::F32), 1, fp);
        fwrite(&x.y, sizeof(PS::F32), 1, fp);
        fwrite(&v.y, sizeof(PS::F32), 1, fp);
        fwrite(&x.z, sizeof(PS::F32), 1, fp);
        fwrite(&v.z, sizeof(PS::F32), 1, fp);
        fwrite(&m,   sizeof(PS::F32), 1, fp);
        fwrite(&i,   sizeof(PS::S32), 1, fp);
    }

    PS::F64 calcDtime(run_param &this_run) {
        PS::F64 dtime_v, dtime_a, dtime;
        PS::F64 vnorm, anorm;
        vnorm = sqrt(SQR(this->vel))+TINY;
        anorm = sqrt(SQR(this->acc))+TINY;
        dtime_v = this->eps/vnorm;
        dtime_a = sqrt(this->eps/anorm)*CUBE(this_run.anow);
        dtime = fmin(0.5*dtime_v, dtime_a);
        return dtime;
    }
};
PS::F64 FP_grav::eps;
PS::F64 FP_grav::H0;
PS::F64 FP_grav::Lbnd;

class EPI_grav {
public:
    PS::S64    id;
    PS::F64vec pos;

    PS::S64 getId() const {
        return this->id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void copyFromFP(const FP_grav & fp) {
        this->id = fp.id;
        this->pos = fp.pos;
    }
};

class EPJ_grav {
public:
    PS::S64    id;
    PS::F64vec pos;
    PS::F64    mass;

    PS::S64 getId() const {
        return this->id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    void copyFromFP(const FP_grav & fp) {
        this->id = fp.id;
        this->mass = fp.mass;
        this->pos = fp.pos;
    }
    void setPos(const PS::F64vec pos_new) {
        this->pos = pos_new;
    }
};


#ifdef ENABLE_PHANTOM_GRAPE_X86
template <class TParticleJ>
class CalcGravity {
public:
    void operator () (const EPI_grav * iptcl,
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
            xj[j][0] = jptcl[j].pos[0];
            xj[j][1] = jptcl[j].pos[1];
            xj[j][2] = jptcl[j].pos[2];
            mj[j]    = jptcl[j].mass;
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
    void operator () (const EPI_grav * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * ep_j,
                      const PS::S32 n_jp,
                      Force_grav * force) {
        const PS::F64 eps2 = FP_grav::eps * FP_grav::eps;
        for (PS::S32 i = 0; i < n_ip; i++) {
            const PS::F64vec xi = ep_i[i].getPos();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for (PS::S32 j = 0; j < n_jp; j++) {
                const PS::F64 mj = ep_j[j].getCharge();
                const PS::F64vec dr = xi - ep_j[j].getPos();
                PS::F64 r3_inv = dr * dr + eps2;
#if defined(CHECK_DIVIDE_BY_ZERO_ERROR_IN_FORCE_KERNEL)
                if (r3_inv == 0.0) continue;
#endif
                const PS::F64 r_inv = 1.0/std::sqrt(r3_inv);
                r3_inv = r_inv * r_inv * r_inv;
                ai -= mj * r3_inv * dr;
                poti -= mj * r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#endif
