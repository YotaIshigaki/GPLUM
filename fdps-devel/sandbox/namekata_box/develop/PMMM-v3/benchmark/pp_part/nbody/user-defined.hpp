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


void CalcGravityEpEp(const FPGrav * ep_i,
                     const PS::S32 n_ip,
                     const FPGrav * ep_j,
                     const PS::S32 n_jp,
                     FPGrav * force) {
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec dr = ep_j[j].getPos() - xi; // rji 
            PS::F64 r3_inv = dr * dr;
            if (r3_inv == 0.0) continue;
            PS::F64 r_inv  = 1.0/sqrt(r3_inv);
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

void CalcGravityEpSp(const FPGrav * ep_i,
                     const PS::S32 n_ip,
                     const PS::SPJQuadrupoleGeometricCenter * ep_j,
                     const PS::S32 n_jp,
                     FPGrav * force) {
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            const PS::F64 mj = ep_j[j].getCharge();
            const PS::F64vec di = ep_j[j].dipole;
            const PS::F64mat quad = ep_j[j].quadrupole;
            PS::F64vec dr = ep_j[j].getPos() - xi; // rji
            PS::F64 r3_inv = dr * dr;
            PS::F64 r = std::sqrt(r3_inv);
            PS::F64 r_inv  = 1.0/r;
            r3_inv  = r_inv * r_inv * r_inv;
            ai     += mj * r3_inv * dr;
            poti   += mj * r_inv;
            //---------------------
            // Dipole terms
            //---------------------
            const PS::F64 r2_inv = r_inv * r_inv;
            const PS::F64vec n = (-dr) * r_inv;
            ai   += (di - 3.0 * (di * n) * n) * r3_inv;
            poti += r2_inv * (di * n); 
            //---------------------
            // Quadrupole terms
            //---------------------
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
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}
