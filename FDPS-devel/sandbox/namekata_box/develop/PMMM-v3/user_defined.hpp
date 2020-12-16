#pragma once

class FP_nbody {
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64    mtot; // total mass of interacting particles
    PS::F64vec acc;
    PS::F64    pot;
    PS::F64vec acc_pm;
    PS::F64    pot_pm;
    PS::F64vec acc_exact;
    PS::F64    pot_exact;

    static PS::F64 eps;

    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F64vec & pos) {
        this->pos = pos;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::S64 getId() const {
        return this->id;
    }
    void copyFromFP(const FP_nbody & fp){ 
        this->id   = fp.id;
        this->mass = fp.mass;
        this->pos  = fp.pos;
    }
    void copyFromForce(const FP_nbody & force) {
        this->mtot = force.mtot;
        this->acc  = force.acc;
        this->pot  = force.pot;
    }
    void copyFromForcePMM(const FP_nbody & force) {
        this->acc_pm = force.acc_pm;
        this->pot_pm = force.pot_pm;
    }
    void accumulateForcePMM(const PS::F64vec & acc, const PS::F64 & pot) {
        this->acc_pm += acc;
        this->pot_pm += pot;
    }
    void accumulateForceDirect(const PS::F64vec & acc, const PS::F64 & pot) {
        this->acc_exact += acc;
        this->pot_exact += pot;
    }
    void accumulateForceEwald(const PS::F64vec & acc, const PS::F64 & pot) {
        this->acc_exact += acc;
        this->pot_exact += pot;
    }
    void clear() {
        this->mtot = 0.0;
        this->acc = 0.0;
        this->pot = 0.0;
    }
    void clearPMM() {
        this->acc_pm = 0.0;
        this->pot_pm = 0.0;
    }
    void clearDirect() {
        this->acc_exact = 0.0;
        this->pot_exact = 0.0;
    }
    void clearEwald() {
        this->acc_exact = 0.0;
        this->pot_exact = 0.0;
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

template <class T>
static inline T tp_sign(T val){
    return (val>0)-(val<0);
}

template <class Tptcl>
void CalcForceEpEp(const FP_nbody * ep_i,
                   const PS::S32 n_ip,
                   const Tptcl * ep_j,
                   const PS::S32 n_jp,
                   FP_nbody * force) {
    for (PS::S32 i = 0; i < n_ip; i++){
        const PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++){
            const PS::F64vec dr = ep_j[j].getPos() - xi; // rji 
            PS::F64 r3_inv = dr * dr;
            if (r3_inv == 0.0) continue;
            PS::F64 r_inv = 1.0/sqrt(r3_inv);
            r3_inv  = r_inv * r_inv;
            r_inv  *= ep_j[j].getCharge();
            r3_inv *= r_inv;
            ai     += r3_inv * dr; // gradient of potential
            poti   += r_inv; // potential
            // Here we assumed that the potential has the form of \phi(r) = m/r,
            // instead of \phi(r) = - m/r. 
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}

template <class Tptcl>
void CalcForceEpSp(const FP_nbody * ep_i,
                   const PS::S32 n_ip,
                   const Tptcl * ep_j,
                   const PS::S32 n_jp,
                   FP_nbody * force) {
    for (PS::S32 i = 0; i < n_ip; i++){
        const PS::F64vec xi = ep_i[i].getPos();
        PS::F64 mtoti = 0.0;
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++){
            const PS::F64 mj = ep_j[j].getCharge();
            const PS::F64vec di = PS::GetMyDipole(ep_j[j]);
            const PS::F64mat quad = PS::GetMyQuadrupole(ep_j[j]);
            const PS::F64vec dr = ep_j[j].getPos() - xi; // rji (not rij) to calculate gradient
            PS::F64 r3_inv = dr * dr;
            const PS::F64 r = std::sqrt(r3_inv);
            const PS::F64 r_inv  = 1.0/r;
            r3_inv = r_inv * r_inv * r_inv;
            mtoti  += mj; // Total charge of interacting particles
            ai     += mj * r3_inv * dr; // gradient of potential
            poti   += mj * r_inv;       // Coulomb's potential
#if 1
            //----------------------
            // Dipole corrections
            //----------------------
            const PS::F64 r2_inv = r_inv * r_inv;
            const PS::F64vec n = (-dr) * r_inv;
            ai   += (di - 3.0 * (di * n) * n) * r3_inv;
            poti += r2_inv * (di * n); 
#endif
#if 1
            //--------------------------
            // Quadrupole corrections
            //--------------------------
            const PS::F64 r4_inv = r2_inv * r2_inv;
            const PS::F64 tr = quad.xx + quad.yy + quad.zz;
            const PS::F64 Qxx = quad.xx - tr/3.0;
            const PS::F64 Qxy = quad.xy;
            const PS::F64 Qxz = quad.xz;
            const PS::F64 Qyy = quad.yy - tr/3.0;
            const PS::F64 Qyz = quad.yz;
            const PS::F64 Qzz = quad.zz - tr/3.0;
            const PS::F64 Q  = ( Qxx * n.x * n.x
                               + 2.0 * Qxy * n.x * n.y 
                               + 2.0 * Qxz * n.x * n.z 
                               + Qyy * n.y * n.y
                               + 2.0 * Qyz * n.y * n.z
                               + Qzz * n.z * n.z);
            const PS::F64vec ex = PS::F64vec(1,0,0);
            const PS::F64vec ey = PS::F64vec(0,1,0);
            const PS::F64vec ez = PS::F64vec(0,0,1);
            const PS::F64vec r_gradQ = ( Qxx * (2.0 * n.x * ex 
                                               - 2.0 * n.x * n.x * n)
                                       + 2.0 * Qxy * (n.x * ey + n.y * ex
                                                     - 2.0 * n.x * n.y * n)
                                       + 2.0 * Qxz * (n.x * ez + n.z * ex
                                                     - 2.0 * n.x * n.z * n)
                                       + Qyy * (2.0 * n.y * ey
                                               - 2.0 * n.y * n.y * n)
                                       + 2.0 * Qyz * (n.y * ez + n.z * ey
                                                     - 2.0 * n.y * n.z * n)
                                       + Qzz * (2.0 * n.z * ez
                                               -2.0 * n.z * n.z * n)
                                       );
            ai   += 1.5 * (r_gradQ - 3.0 * Q * n) * r4_inv;
            poti += 1.5 * Q * r3_inv;
#endif
#if 0
            //------------------
            // Output for debug
            //------------------
            if (ep_i[i].id == 2270) {
                const PS::S64 id = PS::GetMyId(ep_j[j]);
                const PS::F64 id_f = id;
                const PS::F64vec pos = ep_j[j].getPos();
                //std::cout << id_f << " " << pos << " " << mj << std::endl;
                std::cout << pos << "   " << mj << std::endl;
            }
#endif
        }
        force[i].mtot += mtoti;
        force[i].acc += ai;
        force[i].pot += poti;
    }
}

#ifdef USE_FDPS_PMM_EXPERIMENTAL_FEATURE
template <int p>
void CalcForceEpSp(const FP_nbody * ep_i,
                   const PS::S32 n_ip,
                   const PS::SPJMultipoleGeometricCenterPMMM<p> * ep_j,
                   const PS::S32 n_jp,
                   FP_nbody * force) {
    for (PS::S32 i = 0; i < n_ip; i++) {
        const PS::F64vec xi = ep_i[i].getPos();
        PS::F64 mtoti = 0.0;
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++) {
            const PS::F64 mj = ep_j[j].getCharge();
            const PS::F64vec xj = ep_j[j].getPos();
            PS::LocalExpansion0<p> le;
            le.assign_from_MM(ep_j[j].mm, xi, xj); // M2L
            PS::F64 phitmp {0.0};
            PS::F64vec gradtmp(0.0);
            le.read_phi_and_grad(phitmp, gradtmp);
            mtoti  += mj; // Total charge of interacting particles
            ai     += gradtmp; // gradient of potential
            poti   += phitmp;  // Coulomb's potential
        }
        force[i].mtot += mtoti;
        force[i].acc += ai;
        force[i].pot += poti;
    }
}
#endif // USE_FDPS_PMM_EXPERIMENTAL_FEATURE
