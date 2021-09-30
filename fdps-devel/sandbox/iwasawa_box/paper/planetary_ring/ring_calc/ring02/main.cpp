#include<iostream>
#include<cstdio>
#include<unistd.h>
#include<particle_simulator.hpp>
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif

//#define ADD_RNADOM_SHIFT

PS::CommInfo COMM_INFO;

PS::F64 DT_MERGE;
PS::F64 VEL_SHEAR;

PS::S32 TARGET_ADR = 21;

void PrintHelp() {
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

template<typename Tpsys>
void Rotate(Tpsys & psys,
            const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
    for(auto i=0; i<n; i++){
        const auto pos = psys[i].pos;
        const auto vel = psys[i].vel;
        const auto r = sqrt(pos*pos);
        const auto mm = sqrt(1.0/(r*r*r)); // mean mortion
        const auto theta = mm*dt;
        const auto cth = std::cos(theta);
        const auto sth = std::sin(theta);
        const auto l = pos ^ vel;
        const auto n = l / sqrt(l*l);

        const auto ex_pos = pos / sqrt(pos*pos);
        const auto ey_pos = n ^ ex_pos;
        const auto pos_new = sqrt(pos*pos)*(cth*ex_pos + sth*ey_pos);

        const auto ex_vel = vel / sqrt(vel*vel);
        const auto ey_vel = n ^ ex_vel;
        const auto vel_new = sqrt(vel*vel)*(cth*ex_vel + sth*ey_vel);
        const auto l_new = pos_new ^ vel_new;

        psys[i].pos = pos_new;
        psys[i].vel = vel_new;
        
        //if(PS::Comm::getRank() == 0){
        //    std::cerr<<"pos= "<<pos<<" pos_new= "<<pos_new
        //             <<" l= "<<l<<" l_new= "<<l_new
        //             <<std::endl;
        //}
    }
}

PS::F64vec ConvertCar2Cyl(const PS::F64vec & pos){
    const auto pos_r   = sqrt(pos.x*pos.x + pos.y*pos.y);
    const auto pos_phi = atan2(pos.y, pos.x);
    return PS::F64vec(pos_phi, pos_r, pos.z);
}

PS::F64vec ConvertCyl2Car(const PS::F64vec & pos){
    const auto cth = cos(pos.x);
    const auto sth = sin(pos.x);
    const auto r = pos.y;
    const auto pos_x = r*cth;
    const auto pos_y = r*sth;
    return PS::F64vec(pos_x, pos_y, pos.z);
}

class Force_t{
public:
    PS::F64vec acc; // total
    PS::F64    pot;
    PS::F64vec acc_dash;
    //PS::F64    cm_mass;
    //PS::F64vec cm_pos;
    void clear(){
        acc = 0.0;
        pot = 0.0;
        acc_dash = 0.0;
        //cm_mass = 0.0;
        //cm_pos = 0.0;
    }
};

class FP_t{
public:
    PS::F64vec pos; // cartesian
    PS::F64    mass;
    PS::F64vec vel;
    PS::F64vec vel_full; // for calculating disipation energy
    PS::S64    id;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64 r_coll;
    PS::F64vec acc_dash;
    static inline PS::F64 r_search;
    static inline PS::F64 eps;
    static inline PS::F64 kappa;
    static inline PS::F64 eta;
    PS::F64vec getPos() const {
        return ConvertCar2Cyl(pos);
    }
    PS::F64vec getPosCar() const {
        return pos;
    }
    void setPos(const PS::F64vec & pos_new){
        const auto pos_new_car = ConvertCyl2Car(pos_new);
        pos = pos_new_car;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64 getRSearch() const {
        return r_search;
    }
    void copyFromFP(const FP_t & fp){ 
        pos  = fp.pos;
        mass = fp.mass;
        vel  = fp.vel;
        id   = fp.id;
        acc  = fp.acc;
        pot  = fp.pot;
        r_coll = fp.r_coll;
        acc_dash  = fp.acc_dash;
        vel_full  = fp.vel_full;
    }
    void copyFromForce(const Force_t & force) {
        acc = force.acc;
        pot = force.pot;
        acc_dash = force.acc_dash;
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

using EPI_t = FP_t;
using EPJ_t = FP_t;

void DivideNProc(PS::S32 & nx,
                 PS::S32 & ny,
                 PS::S32 & nz,
                 const PS::S32 n_proc,
                 const PS::F64 delta_ax){
    nz = 1;
    ny = 1;
    nx = n_proc / ny;
    if(n_proc == 1) return;
    const PS::F64 dx = 2.0*M_PI   / nx;
    double dy = delta_ax / ny;
    double ratio = (dx < dy) ? dx / dy : dy / dx;
    PS::S32 ny_tmp = ny;
    PS::S32 nx_tmp = nx;
    double dx_tmp = dx;
    double dy_tmp = dy;
    double ratio_tmp = ratio;
    do{
        ny = ny_tmp;
        nx = nx_tmp;
        ratio = ratio_tmp;
        ny_tmp += 1;
        while( n_proc % ny_tmp != 0) ny_tmp++;
        nx_tmp = n_proc / ny_tmp;
        dx_tmp = 2.0*M_PI   / nx_tmp;
        dy_tmp = delta_ax / ny_tmp;
        ratio_tmp = (dx_tmp < dy_tmp) ? dx_tmp / dy_tmp : dy_tmp / dx_tmp;
    }while( fabs(ratio_tmp-1.0) < fabs(ratio-1.0));
    //PS::Comm::barrier();
    COMM_INFO.barrier();
    if (COMM_INFO.getRank() == 0){
        std::cerr<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
    }
    assert(n_proc == nx*ny*nz);
}

void GetPosDomainCyl(const PS::F64 delta_ax,
                     PS::F64ort & pos_domain,
                     const PS::S32 nx,
                     const PS::S32 ny){
    constexpr PS::F64 len_x = 2.0 * M_PI;
    //PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 my_rank = COMM_INFO.getRank();
    PS::S32 rank_x = my_rank / ny;
    PS::S32 rank_y = my_rank % ny;
    PS::F64 dx = len_x / nx;
    PS::F64 dy = delta_ax / ny;
    PS::F64 dy_offset = 1.0-delta_ax*0.5;
    PS::F64 x_offset = -M_PI;
    pos_domain.low_.x = dx*rank_x + x_offset;
    pos_domain.low_.y = dy*rank_y + dy_offset;
    pos_domain.low_.z = -M_PI;
    pos_domain.high_.x = dx*(rank_x+1) + x_offset;
    pos_domain.high_.y = dy*(rank_y+1) + dy_offset;
    pos_domain.high_.z = M_PI;
}

PS::F64vec GetVel(const PS::F64vec & pos){
    const PS::F64 dr = sqrt(pos*pos);
    const PS::F64 v_kep = sqrt(1.0/(dr));
    const PS::F64 theta = atan2(pos.y, pos.x) + M_PI*0.5;
    const PS::F64 vx = v_kep * cos(theta);
    const PS::F64 vy = v_kep * sin(theta);
    const PS::F64 vz = 0.0;
    return PS::F64vec(vx, vy, vz);
}

template<class Tpsys>
void SetParticleKeplerDiskCyl(Tpsys & psys,
                              const PS::S64 n_glb,
                              const PS::F64 ax_in,
                              const PS::F64 ax_out,
                              const PS::F64 dens, // surface density at ax_in
                              const PS::F64ort box, // domain
                              const PS::F64 r_phy, // radius of particle
                              const bool layer=true){
    PS::MTTS mt;
    //mt.init_genrand( PS::Comm::getRank() );
    mt.init_genrand( COMM_INFO.getRank() );
    const PS::F64 area = M_PI*(ax_out*ax_out-ax_in*ax_in);

    PS::F64 dS = (area*3.0) / n_glb;
    if(!layer){
        dS = area / n_glb;
    }

    const PS::F64 mass = (dens*area) / n_glb; // particle mass
    const PS::F64 dl = sqrt(dS);
    const PS::F64 dz = dl;
    const PS::F64 ax_cen = (ax_in+ax_out) * 0.5;
    const PS::S64 n_theta = (2.0*M_PI*ax_cen / dl);
    const PS::F64 dtheta = 2.0*M_PI / n_theta;
    const PS::F64 dax = ax_out - ax_in;
    assert(dax > dl);
    const PS::S64 n_r = (dax / dl);
    const PS::F64 dr = dax / n_r;

    const PS::S64 id_x_head = box.low_.x / dtheta;
    const PS::S64 id_x_tail = box.high_.x / dtheta;
    const PS::S64 id_y_head = box.low_.y / dr;
    const PS::S64 id_y_tail = box.high_.y / dr;
    const PS::F64 offset_theta = dtheta*(sqrt(2.0)-1.0)*0.5;

    const PS::F64 eps = (dl > 2.0*r_phy) ? (0.5*(dl-2.0*r_phy))*0.9 : 0.0;
    
    COMM_INFO.barrier();
    if(COMM_INFO.getRank()==0){
        std::cerr<<"dax= "<<dax
                 <<" dS = "<<dS
                 <<" dl = "<<dl
                 <<" n_r= "<<n_r
                 <<" dr= "<<dr
                 <<" eps= "<<eps
                 <<std::endl;
        std::cerr<<"n_theta= "<<n_theta
                 <<" dtheta= "<<dtheta
                 <<" offset_theta= "<<offset_theta
                 <<std::endl;
    }
    COMM_INFO.barrier();
    std::vector<PS::F64vec> pos;
    for(PS::S64 i=id_x_head-3; i<=id_x_tail+3; i++){
        PS::F64 pos_x = ((PS::F64)i)*dtheta + offset_theta;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-3; j<=id_y_tail+3; j++){
            PS::F64 pos_y = j*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
            PS::F64 eps_x, eps_y, eps_z;
            eps_x = eps_y = eps_z = 0.0;
#ifdef ADD_RNADOM_SHIFT
            eps_x = eps*(mt.genrand_res53()-0.5)*2.0;
            eps_y = eps*(mt.genrand_res53()-0.5)*2.0;
            eps_z = eps*(mt.genrand_res53()-0.5)*2.0;
#endif
            pos.push_back(PS::F64vec(pos_x+eps_x, pos_y+eps_y, eps_z));
        }
    }

    if(layer == true){
    
    for(PS::S64 i=id_x_head-3; i<=id_x_tail+3; i++){
        PS::F64 pos_x = ((PS::F64)i)*dtheta+0.5*dtheta + offset_theta;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-3; j<=id_y_tail+3; j++){
            PS::F64 pos_y = j*dr+0.5*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
            PS::F64 eps_x, eps_y, eps_z;
            eps_x = eps_y = eps_z = 0.0;
#ifdef ADD_RNADOM_SHIFT
            eps_x = eps*(mt.genrand_res53()-0.5)*2.0;
            eps_y = eps*(mt.genrand_res53()-0.5)*2.0;
            eps_z = eps*(mt.genrand_res53()-0.5)*2.0;
#endif
            pos.push_back(PS::F64vec(pos_x+eps_x, pos_y+eps_y, -0.5*dz+eps_z));
        }
    }
    for(PS::S64 i=id_x_head-3; i<=id_x_tail+3; i++){
        PS::F64 pos_x = ((PS::F64)i)*dtheta+0.5*dtheta + offset_theta;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-3; j<=id_y_tail+3; j++){
            PS::F64 pos_y = j*dr+0.5*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
            PS::F64 eps_x, eps_y, eps_z;
            eps_x = eps_y = eps_z = 0.0;
#ifdef ADD_RNADOM_SHIFT
            eps_x = eps*(mt.genrand_res53()-0.5)*2.0;
            eps_y = eps*(mt.genrand_res53()-0.5)*2.0;
            eps_z = eps*(mt.genrand_res53()-0.5)*2.0;
#endif
            pos.push_back(PS::F64vec(pos_x+eps_x, pos_y+eps_y, 0.5*dz+eps_z));
        }
    }
    }
    PS::S64 n_loc = pos.size();
    psys.setNumberOfParticleLocal(n_loc);
    for(PS::S64 i=0; i<n_loc; i++){
        PS::F64 pos_x = pos[i].y*cos(pos[i].x);
        PS::F64 pos_y = pos[i].y*sin(pos[i].x);
        PS::F64 pos_z = pos[i].z;
        psys[i].pos.x  = pos_x;
        psys[i].pos.y  = pos_y;
        psys[i].pos.z  = pos_z;
        psys[i].vel  = GetVel(psys[i].pos);
        psys[i].mass = mass;
        psys[i].r_coll = r_phy;
    }
}

#ifdef ENABLE_PHANTOM_GRAPE_X86
template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceEp{
    void operator () (const Tpi * iptcl,
                      const PS::S32 ni,
                      const Tpj * jptcl,
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
        g5_set_xmj(0, nj, xj, mj);
        g5_set_n(nj);
        g5_calculate_force_on_x(xi, ai, pi, ni);
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
#elif defined(GRAV_ONLY)
template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceEp{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FP_t::eps*FP_t::eps;
        for(auto i=0; i<ni; i++){
            if(pi[i].id == TARGET_ADR){
                std::cerr<<"ep: nj= "<<nj<<std::endl;
            }
            const PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64vec ai_dash = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                const auto r_coll    = (pi[i].r_coll + pj[j].r_coll);
                const auto r_coll_sq = r_coll*r_coll;
                PS::F64vec rij    = xi - pj[j].getPosCar();
                if(pi[i].id == pj[j].id) continue;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64 r_inv  = 1.0/sqrt(r2);
                PS::F64 r2_inv  = r_inv * r_inv;
                PS::F64 pot = r_inv * pj[j].getCharge();
                ai     -= pot * r2_inv * rij;
                poti   -= 0.5 * pot;
            }
            force[i].acc += ai;
            force[i].acc_dash += ai_dash;
            force[i].pot += poti;
        }
    }
};
#else
template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceEp{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FP_t::eps*FP_t::eps;
        const auto kappa = FP_t::kappa;
        const auto eta   = FP_t::eta;
        for(auto i=0; i<ni; i++){
            const PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64vec ai_dash = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                const auto r_coll    = (pi[i].r_coll + pj[j].r_coll);
                const auto r_coll_sq = r_coll*r_coll;
                PS::F64vec rij       = xi - pj[j].getPosCar();
                if(pi[i].id == pj[j].id) continue;
                PS::F64 r2_real = rij * rij + eps2;
                PS::F64 r2      = std::max(r2_real, r_coll_sq);
                PS::F64 r_inv   = 1.0/sqrt(r2);
                PS::F64 r2_inv  = r_inv * r_inv;
                PS::F64 pot = r_inv * pj[j].getCharge();
                if(r_coll_sq > r2_real){
                    ai     -= pj[j].getCharge() / (r_coll_sq*r_coll) * rij;
                    PS::F64 pot_offset = -1.5/r_coll;
                    poti   += 0.5*pj[j].getCharge()*(0.5*r2_real/(r_coll_sq*r_coll) + pot_offset);
                }
                else{
                    ai     -= pot * r2_inv * rij;
                    poti   -= 0.5 * pot;
                }
                //poti   -= 0.5 * pj[j].getCharge() * sqrt(r2_real) * r2_inv;
                if(r_coll_sq > r2_real){
                    PS::F64 m_r = pj[j].mass / (pi[i].mass+pj[j].mass);
                    PS::F64 r   = sqrt(r2_real);
                    PS::F64 dr  = r_coll-r ;
                    ai += kappa * m_r * dr/r * rij;
                    poti += 0.5*kappa*m_r*dr*dr * 0.5;
                    PS::F64vec vij = pi[i].vel_full - pj[j].vel_full;
                    PS::F64 rv = rij*vij;
                    //PS::F64vec a_eta = eta * m_r * rv * r2_inv * rij;
                    PS::F64vec a_eta = eta * m_r * rv / r2_real * rij;
                    ai_dash += a_eta;
                    ai += a_eta;
                }
            }
            force[i].acc += ai;
            force[i].acc_dash += ai_dash;
            force[i].pot += poti;
        }
    }
};
#endif


template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpMono{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
#if 1
        const auto eps2 = FP_t::eps*FP_t::eps;
        for(auto i=0; i<ni; i++){
            const auto xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            if(pi[i].id == TARGET_ADR){
                std::cerr<<"sp: nj= "<<nj<<std::endl;
            }            
            for(auto j=0; j<nj; j++){
                PS::F64vec rij    = xi - pj[j].getPosCar();
                PS::F64    r3_inv = rij * rij + eps2;
                PS::F64    r_inv  = 1.0/sqrt(r3_inv);
                r3_inv  = r_inv * r_inv;
                r_inv  *= pj[j].getCharge();
                r3_inv *= r_inv;
                ai     -= r3_inv * rij;
                poti   -= 0.5*r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
#endif
    }
};

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpQuad{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FP_t::eps*FP_t::eps;
        for(auto i=0; i<ni; i++){
            PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                PS::F64 mj = pj[j].getCharge();
                PS::F64vec xj= pj[j].getPosCar();
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = pj[j].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[i].acc += ai;
            force[i].pot += 0.5*poti;
        }
    }
};


class MyMomentMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64ort boundary;
    MyMomentMonopole() : mass(0.0), pos(PS::F64vec(0.0)), pos_car(PS::F64vec(0.0)), boundary(PS::F64ort(PS::F64vec(-100.0), PS::F64vec(100.0))){}
    MyMomentMonopole(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & p_car, const PS::F64ort & b) : mass(m), pos(p), pos_car(p_car), boundary(b){}
    void init(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        mass += epj.getCharge();
        pos += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
        boundary.merge(epj.getPos());
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
    }
    void accumulate(const MyMomentMonopole & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
        boundary.merge(mom.boundary);
    }
    void accumulate2(const MyMomentMonopole & mom){}
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"pos_car="<<pos_car<<std::endl;
    }
};

class MySPJMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64ort boundary;
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        this->mass     = mom.mass;
        this->pos      = mom.pos;
        this->pos_car  = mom.pos_car;
        this->boundary = mom.boundary;
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void setPos(const PS::F64vec & pos_new) {
        pos = pos_new;
    }
    MyMomentMonopole convertToMoment() const {
        return MyMomentMonopole(mass, pos, pos_car, boundary);
    }
};

class MyMomentQuadrupole{
public:
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64mat quad;
    PS::F64vec pos_car;
    void init(){
        pos = 0.0;
        mass = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(const PS::F64 m, const PS::F64vec & p, const PS::F64mat & q, const PS::F64vec & p_car){
        mass = m;
        pos = p;
        quad = q;
        pos_car = p_car;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        mass += epj.getCharge();
        pos  += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){
        PS::F64 ctmp = epj.getCharge();
        PS::F64vec ptmp = epj.getPosCar() - this->pos_car;
        PS::F64 cx = ctmp * ptmp.x;
        PS::F64 cy = ctmp * ptmp.y;
        PS::F64 cz = ctmp * ptmp.z;
        this->quad.xx += cx * ptmp.x;
        this->quad.yy += cy * ptmp.y;
        this->quad.zz += cz * ptmp.z;
        this->quad.xy += cx * ptmp.y;
        this->quad.xz += cx * ptmp.z;
        this->quad.yz += cy * ptmp.z;
    }
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
    }
    void accumulate(const MyMomentQuadrupole & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
    }
    void accumulate2(const MyMomentQuadrupole & mom){
        PS::F64 mtmp = mom.mass;
        PS::F64vec ptmp = mom.pos_car - this->pos_car;
        PS::F64 cx = mtmp * ptmp.x;
        PS::F64 cy = mtmp * ptmp.y;
        PS::F64 cz = mtmp * ptmp.z;
        this->quad.xx += cx * ptmp.x + mom.quad.xx;
        this->quad.yy += cy * ptmp.y + mom.quad.yy;
        this->quad.zz += cz * ptmp.z + mom.quad.zz;
        this->quad.xy += cx * ptmp.y + mom.quad.xy;
        this->quad.xz += cx * ptmp.z + mom.quad.xz;
        this->quad.yz += cy * ptmp.z + mom.quad.yz;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass= "<<mass<<std::endl;
        fout<<"pos= "<<pos<<std::endl;
        fout<<"quad= "<<quad<<std::endl;
    }
};

class MySPJQuadrupole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64mat quad;
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void setPos(const PS::F64vec & pos_new) {
        pos = pos_new;
    }
    void copyFromMoment(const MyMomentQuadrupole & mom){
        mass = mom.mass;
        pos = mom.pos;
        quad = mom.quad;
        pos_car = mom.pos_car;
    }
    MyMomentQuadrupole convertToMoment() const {
        return MyMomentQuadrupole(mass, pos, quad, pos_car);
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
};

FP_t PLANET;

#ifdef QUAD
using SPJ_t    = MySPJQuadrupole;
using Moment_t = MyMomentQuadrupole;
using CalcForceSp = CalcForceSpQuad<EPI_t, SPJ_t, Force_t>;
#else
using SPJ_t    = MySPJMonopole;
using Moment_t = MyMomentMonopole;
    #ifdef ENABLE_PHANTOM_GRAPE_X86
using CalcForceSp = CalcForceEp<EPI_t, SPJ_t, Force_t>;
    #else
using CalcForceSp = CalcForceSpMono<EPI_t, SPJ_t, Force_t>;
    #endif
#endif

using MY_SEARCH_MODE = PS::SEARCH_MODE_LONG_SCATTER;
using Tree_t = PS::TreeForForce<MY_SEARCH_MODE, Force_t, EPI_t, EPJ_t, Moment_t, Moment_t, SPJ_t>;
//using Tree_t = PS::TreeForForce<PS::SEARCH_MODE_LONG_SCATTER, Force_t, EPI_t, EPJ_t, Moment_t, Moment_t, SPJ_t>;
//using Tree_t = PS::TreeForForce<PS::SEARCH_MODE_LONG_SYMMETRY, Force_t, EPI_t, EPJ_t, Moment_t, Moment_t, SPJ_t>;

template<typename Ttree, typename Tepj, typename Tspj>
inline void MyMakeInteractionList(Ttree & tree, const PS::DomainInfo & dinfo){
    
    COMM_INFO.barrier();
    if(COMM_INFO.getRank()==0) std::cerr<<"check d1"<<std::endl;
    
    using namespace PS;
    const S32 n_thread = Comm::getNumberOfThread();
    std::vector<ReallocatableArray<S32>> adr_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_spj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_ipg_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_disp_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_disp_spj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));

    //PS::Comm::barrier();
    //if(PS::Comm::getRank()==0) std::cerr<<"check d2"<<std::endl;
    
    const S32 n_ipg = tree.ipg_.size();
    //tree.n_walk_local_ += n_ipg;
    const S32 adr_tc = 0;
    const S32 adr_tree_sp_first = tree.spj_sorted_.size() - tree.tc_glb_.size();
    tree.interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
    tree.interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
    tree.interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
    tree.interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
    const auto len_peri = dinfo.getPosRootDomain().getFullLength();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
    {
        const S32 ith = Comm::getThreadNum();
        S32 n_ep_cum_prev = 0;
        S32 n_sp_cum_prev = 0;
        adr_epj_tmp[ith].clearSize();
        adr_spj_tmp[ith].clearSize();
        adr_ipg_tmp[ith].clearSize();
        n_disp_epj_tmp[ith].clearSize();
        n_disp_spj_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
        for(S32 i=0; i<n_ipg; i++){
            adr_ipg_tmp[ith].push_back(i);
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
            n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
            //TargetBox<SEARCH_MODE_LONG_SCATTER> target_box;
            //TargetBox<SEARCH_MODE_LONG_SYMMETRY> target_box;
            TargetBox<MY_SEARCH_MODE> target_box;
            target_box.set(tree.ipg_[i]);
            MakeListUsingTreeRecursiveTop
                //<SEARCH_MODE_LONG_SCATTER, TreeCell<Moment_t, Geometry<typename SEARCH_MODE_LONG_SCATTER::tree_cell_loc_geometry_type>>,
                //<SEARCH_MODE_LONG_SYMMETRY, TreeCell<Moment_t, Geometry<typename SEARCH_MODE_LONG_SYMMETRY::tree_cell_loc_geometry_type>>,
                <MY_SEARCH_MODE, TreeCell<Moment_t, Geometry<typename MY_SEARCH_MODE::tree_cell_loc_geometry_type>>,
                 TreeParticle, Tepj, Tspj,
                 TagWalkModeNearestImage, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl>
                (tree.tc_glb_,  adr_tc, tree.tp_glb_,
                 tree.epj_sorted_, adr_epj_tmp[ith],
                 tree.spj_sorted_, adr_spj_tmp[ith],
                 target_box,
                 tree.n_leaf_limit_,
                 adr_tree_sp_first, len_peri, tree.theta_);
            tree.interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
            tree.interaction_list_.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;
            n_ep_cum_prev = adr_epj_tmp[ith].size();
            n_sp_cum_prev = adr_spj_tmp[ith].size();
        }
        n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
        n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
    } // end of OMP

    //PS::Comm::barrier();
    //if(PS::Comm::getRank()==0) std::cerr<<"check d3"<<std::endl;
        
    tree.interaction_list_.n_disp_ep_[0] = 0;
    tree.interaction_list_.n_disp_sp_[0] = 0;
    for(S32 i=0; i<n_ipg; i++){
        tree.interaction_list_.n_disp_ep_[i+1] = tree.interaction_list_.n_disp_ep_[i] + tree.interaction_list_.n_ep_[i];
        tree.interaction_list_.n_disp_sp_[i+1] = tree.interaction_list_.n_disp_sp_[i] + tree.interaction_list_.n_sp_[i];
    }
    tree.interaction_list_.adr_ep_.resizeNoInitialize( tree.interaction_list_.n_disp_ep_[n_ipg] );
    tree.interaction_list_.adr_sp_.resizeNoInitialize( tree.interaction_list_.n_disp_sp_[n_ipg] );

    //PS::Comm::barrier();
    //if(PS::Comm::getRank()==0) std::cerr<<"check d4"<<std::endl;
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
    for(S32 i=0; i<n_thread; i++){
        for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
            const S32 adr_ipg = adr_ipg_tmp[i][j];
            S32 adr_ep = tree.interaction_list_.n_disp_ep_[adr_ipg];
            const S32 k_ep_h = n_disp_epj_tmp[i][j];
            const S32 k_ep_e = n_disp_epj_tmp[i][j+1];
            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                tree.interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
            }
            S32 adr_sp = tree.interaction_list_.n_disp_sp_[adr_ipg];
            const S32 k_sp_h = n_disp_spj_tmp[i][j];
            const S32 k_sp_e = n_disp_spj_tmp[i][j+1];
            for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                tree.interaction_list_.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
            }
        }
    }
    //PS::Comm::barrier();
    //if(PS::Comm::getRank()==0) std::cerr<<"check d5"<<std::endl;
}

template<typename Tfunc_ep_ep, typename Tfunc_ep_sp, typename Tpsys, typename Ttree>
inline void MyCalcForceAllAndWriteBack(Ttree & tree,
                                Tfunc_ep_ep pfunc_ep_ep,
                                Tfunc_ep_sp pfunc_ep_sp,
                                Tpsys & psys,
                                PS::DomainInfo & dinfo,
                                const bool clear_force=true,
                                const PS::INTERACTION_LIST_MODE list_mode = PS::MAKE_LIST_FOR_REUSE,
                                const bool flag_serialize=false){
    std::cerr<<std::setprecision(15);
    tree.setParticleLocalTree(psys, true);
    if(list_mode == PS::MAKE_LIST || list_mode == PS::MAKE_LIST_FOR_REUSE){
        PS::F64 shift_z_loc = 0.0;
        for(auto i=0; i<psys.getNumberOfParticleLocal(); i++){
            shift_z_loc = std::max(fabs(shift_z_loc), psys[i].getPos().z);
        }
        //PS::F64 shift_z_glb = PS::Comm::getMaxValue(shift_z_loc);
        PS::F64 shift_z_glb = COMM_INFO.getMaxValue(shift_z_loc);
        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"shift_z_glb= "<<shift_z_glb<<std::endl;
        PS::F64vec cen = PS::F64vec(0.0, 0.0, shift_z_glb*2.01);
        PS::F64    len = 2.0*M_PI;
        tree.setRootCell(len, cen);
        tree.mortonSortLocalTreeOnly();
        tree.linkCellLocalTreeOnly();
        tree.calcMomentLocalTreeOnly();
        
#ifdef EXPAND_BOX
        for(auto i=0; i<tree.tc_loc_.size(); i++){
            tree.tc_loc_[i].geo_.size_ += FP_t::r_search*2.0;
            tree.tc_loc_[i].geo_.vertex_in_.low_  -= PS::F64vec(FP_t::r_search);
            tree.tc_loc_[i].geo_.vertex_in_.high_  += PS::F64vec(FP_t::r_search);
            tree.tc_loc_[i].geo_.vertex_out_.low_  -= PS::F64vec(FP_t::r_search);
            tree.tc_loc_[i].geo_.vertex_out_.high_  += PS::F64vec(FP_t::r_search);
        }
#endif
        /*
        PS::F64 cm_mass_loc = 0.0;
        PS::F64vec cm_pos_cyl_loc = 0.0;
        PS::F64vec cm_pos_car_loc = 0.0;
        for(auto i=0; i<psys.getNumberOfParticleLocal(); i++){
            cm_mass_loc += psys[i].mass;
            cm_pos_cyl_loc  += psys[i].mass*psys[i].getPos();
            cm_pos_car_loc  += psys[i].mass*psys[i].getPosCar();
        }
        cm_pos_cyl_loc /= cm_mass_loc;
        */
        
        //std::cerr<<"cm_mass_loc= "<<cm_mass_loc
        //         <<" cm_pos_cyl_loc= "<<cm_pos_cyl_loc
        //         <<std::endl;
        //std::cerr<<"tc_loc_[0].mom_.mass= "<<tree.tc_loc_[0].mom_.mass
        //         <<" tc_loc_[0].mom_.pos= "<<tree.tc_loc_[0].mom_.pos
        //         <<std::endl;

        //PS::F64 cm_mass_dir    = PS::Comm::getSum(cm_mass_loc);
        //PS::F64vec cm_pos_car_dir = PS::Comm::getSum(cm_pos_car_loc);
        //if(PS::Comm::getRank()==0){
        //    std::cerr<<"cm_mass_dir= "<<cm_mass_dir
        //             <<" cm_pos_car_dir= "<<cm_pos_car_dir
        //             <<std::endl;
        //}

        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"check a"<<std::endl;
        
        tree.template exchangeLocalEssentialTreeLong<PS::TagWalkModeNearestImage>(dinfo, false);

        /*
        if(PS::Comm::getRank()==0){
            for(auto i=0; i<PS::Comm::getNumberOfProc(); i++){
                std::cerr<<"i= "<<i
                         <<" tree.comm_table_.n_ep_send_[i]= "<<tree.comm_table_.n_ep_send_[i]
                         <<" tree.comm_table_.n_sp_send_[i]= "<<tree.comm_table_.n_sp_send_[i]
                         <<" tree.comm_table_.n_ep_recv_[i]= "<<tree.comm_table_.n_ep_recv_[i]
                         <<" tree.comm_table_.n_sp_recv_[i]= "<<tree.comm_table_.n_sp_recv_[i]
                         <<std::endl;
            }
        }
        */
        
        /*
        PS::F64 cm_mass_glb = 0.0;
        for(auto i=0; i<tree.epj_org_.size(); i++){
            cm_mass_glb += tree.epj_org_[i].mass;
        }
        for(auto i=0; i<tree.spj_org_.size(); i++){
            cm_mass_glb += tree.spj_org_[i].mass;
        }
        std::cerr<<"my_rank= "<<PS::Comm::getRank()
                 <<" cm_mass_glb= "<<cm_mass_glb
                 <<" cm_mass_dir= "<<cm_mass_dir
                 <<std::endl;
        */

        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"check b"<<std::endl;
        
        tree.setLocalEssentialTreeToGlobalTree(false);

        tree.mortonSortGlobalTreeOnly();
        tree.linkCellGlobalTreeOnly();
        tree.calcMomentGlobalTreeOnly();

#ifdef EXPAND_BOX
        for(auto i=0; i<tree.tc_loc_.size(); i++){
            //tree.tc_glb_[i].geo_.size_ += FP_t::r_search*2.0;
            //tree.tc_glb_[i].geo_.vertex_in_.low_  -= PS::F64vec(FP_t::r_search);
            //tree.tc_glb_[i].geo_.vertex_in_.high_  += PS::F64vec(FP_t::r_search);
            //tree.tc_glb_[i].geo_.vertex_out_.low_  -= PS::F64vec(FP_t::r_search);
            //tree.tc_glb_[i].geo_.vertex_out_.high_  += PS::F64vec(FP_t::r_search);
        }
#endif
        
        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"check c"<<std::endl;
        
        //if(PS::Comm::getRank()==0){        
        //    std::cerr<<"tc_glb_[0].mom_.mass= "<<tree.tc_glb_[0].mom_.mass
        //             <<" tc_glb_[0].mom_.pos= "<<tree.tc_glb_[0].mom_.pos
        //             <<std::endl;
        //}

        tree.addMomentAsSp();
        tree.makeIPGroup();
        //tree.n_walk_local_ += tree.ipg_.size();
#ifdef EXPAND_BOX
        for(auto i=0; i<tree.ipg_.size(); i++){
            tree.ipg_[i].vertex_in_.low_  -= PS::F64vec(FP_t::r_search);
            tree.ipg_[i].vertex_in_.high_ += PS::F64vec(FP_t::r_search);
            //tree.ipg_[i].vertex_out_.low_  -= PS::F64vec(FP_t::r_search);
            //tree.ipg_[i].vertex_out_.high_ += PS::F64vec(FP_t::r_search);
        }
#endif
        
        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"check d"<<std::endl;
        
        MyMakeInteractionList<Tree_t, EPJ_t, SPJ_t>(tree, dinfo);

        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"check e"<<std::endl;
        
        bool clear_force = true;
        tree.calcForceNoWalk(pfunc_ep_ep, pfunc_ep_sp, clear_force);

        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"check f"<<std::endl;
        
        //PS::Finalize();
        //exit(1);
        
        //if(PS::Comm::getRank()==0){
        //    for(auto i=0; i<psys.getNumberOfParticleLocal(); i++){
        //        std::cerr<<"tree.force_org_[i].cm_mass= "<<tree.force_org_[i].cm_mass
        //                 <<" cm_pos= "<<tree.force_org_[i].cm_pos<<std::endl;
        //    }
        //}
    }
    else if(list_mode == PS::REUSE_LIST){
        tree.mortonSortLocalTreeOnly(true);
        tree.calcMomentLocalTreeOnly();
        tree.template exchangeLocalEssentialTreeLong<PS::TagWalkModeNearestImage>(dinfo, true);
        tree.setLocalEssentialTreeToGlobalTree(true);
        tree.mortonSortGlobalTreeOnly(true);
        tree.calcMomentGlobalTreeOnly();
        //AddMomentAsSpImpl(typename PS::SEARCH_MODE_LONG_SCATTER::force_type(), tree.tc_glb_, tree.spj_sorted_.size(), tree.spj_sorted_);
        //AddMomentAsSpImpl(typename PS::SEARCH_MODE_LONG_SYMMETRY::force_type(), tree.tc_glb_, tree.spj_sorted_.size(), tree.spj_sorted_);
        AddMomentAsSpImpl(typename MY_SEARCH_MODE::force_type(), tree.tc_glb_, tree.spj_sorted_.size(), tree.spj_sorted_);
        tree.n_walk_local_ += tree.ipg_.size();
        tree.calcForceNoWalk(pfunc_ep_ep, pfunc_ep_sp, clear_force);
    }
    else{
        PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
        std::cerr<<"INTERACTION_LIST_MODE: "<<list_mode<<std::endl;
        PS::Abort(-1);
    }
    tree.freeObjectFromMemoryPool();
    tree.writeBack(psys);
}

template<typename Tpsys>
void CalcForceFromPlanet(Tpsys & psys, const FP_t & pla){
    const auto n = psys.getNumberOfParticleLocal();
    for(auto i=0; i<n; i++){
        const auto rij = psys[i].pos - pla.pos;
        const auto r_sq = rij*rij;
        const auto r_inv = 1.0 / sqrt(r_sq);
        const auto pot   = pla.mass * r_inv;
        psys[i].acc -= pot * r_inv * r_inv * rij;
        psys[i].pot -= pot;
    }
}

template<typename Tpi>
void Kick(PS::ReallocatableArray<Tpi> & pi,
          const PS::F64 dt,
          const PS::S32 n){
    for(auto i=0; i<n; i++){
        pi[i].vel += pi[i].acc*dt;
        pi[i].vel_full = pi[i].vel + pi[i].acc*dt;
    }
}
template<typename Tpi>
void Drift(PS::ReallocatableArray<Tpi> & pi,
           const PS::F64 dt,
           const PS::S32 n){
    for(auto i=0; i<n; i++){
        pi[i].pos += pi[i].vel*dt;
    }
}
template<typename Tpi, typename Tpj>
void Copy(const PS::ReallocatableArray<Tpi> & pi,
          PS::ReallocatableArray<Tpj> & pj,
          const PS::S32 n){
    for(auto i=0; i<n; i++){
        pj[i] = pi[i];
    }
}
template<typename Tpsys>
    void Kick(Tpsys & psys,
          const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
    for(auto i=0; i<n; i++){
        psys[i].vel += psys[i].acc*dt;
        psys[i].vel_full = psys[i].vel + psys[i].acc*dt;
    }
}
template<typename Tpsys>
void Drift(Tpsys & psys,
           const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
    for(auto i=0; i<n; i++){
        psys[i].pos += psys[i].vel*dt;
    }
}

class Energy{
public:
    PS::F64 tot;
    PS::F64 pot;
    PS::F64 kin;
    PS::F64 eng_disp;
    Energy(){
        pot = kin = tot = eng_disp = 0.0;
    }
    template<typename Tpsys>
    void calc(const Tpsys & psys){
        const auto n = psys.getNumberOfParticleLocal();
        PS::F64 pot_loc = 0.0;
        PS::F64 kin_loc = 0.0;
        for(auto i=0; i<n; i++){
            kin_loc += 0.5*psys[i].mass*psys[i].vel*psys[i].vel;
            pot_loc += psys[i].mass*psys[i].pot;
        }
        pot = COMM_INFO.getSum(pot_loc);
        kin = COMM_INFO.getSum(kin_loc);
        tot = pot + kin;
    }
    template<typename Tptcl, typename Tforce>
    void calc(const Tptcl  * psys,
              const Tforce * force,
              const FP_t & planet,
              const int n){
        PS::F64 pot_loc = 0.0;
        PS::F64 kin_loc = 0.0;
        for(auto i=0; i<n; i++){
            kin_loc += 0.5*psys[i].mass*psys[i].vel*psys[i].vel;
            pot_loc += psys[i].mass*force[i].pot;
            PS::F64vec rij = planet.pos - psys[i].pos;
            pot_loc -= psys[i].mass*planet.mass / sqrt(rij*rij);
        }
        pot = COMM_INFO.getSum(pot_loc);
        kin = COMM_INFO.getSum(kin_loc);
        tot = pot + kin;
    }
    void dump(std::ostream & fout = std::cerr){
        fout<<"tot= "<<tot<<" pot= "<<pot<<" kin= "<<kin<<std::endl;
    }
};

template<typename Tpsys>
void SetID(Tpsys & psys){
    PS::S32 n = psys.getNumberOfParticleLocal();
    PS::S32 n_proc = COMM_INFO.getNumberOfProc();
    PS::S32 my_rank = COMM_INFO.getRank();
    PS::ReallocatableArray<PS::S32> n_ar;
    n_ar.resizeNoInitialize(n_proc);
    //PS::Comm::allGather(&n, 1, n_ar.getPointer());
    COMM_INFO.allGather(&n, 1, n_ar.getPointer());
    PS::S32 offset = 0;
    for(auto i=0; i<my_rank; i++){
        offset += n_ar[i];
    }
    for(auto i=0; i<n; i++){
        psys[i].id = i + offset;
    }
    //std::cerr<<"psys[0].id= "<<psys[0].id<<std::endl;
}

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    //PS::Initialize(argc, argv, 1500000000);
    PS::Initialize(argc, argv, (long long int)(3.0e10));

    //COMM_INFO = PS::Comm::split(PS::Comm::getRank()%2, PS::Comm::getRank());
    
    const auto my_rank = COMM_INFO.getRank();
    const auto n_proc  = COMM_INFO.getNumberOfProc();
    PS::S64 n_loop_merge = 1;

    PS::F64 time_sys = 0.0;
    PS::F64 dt = 1.0 / 256.0;
    PS::S32 n_smp = 100;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F64 theta = 0.5;
    PS::F64 ax_cen = 1.0;
    //PS::F64 delta_ax = 0.5;
    PS::F64 delta_ax = 1e-3;
    PS::S64 n_glb = 1024;
    PS::S64 n_loc = 0;
    PS::F64 time_end = 1.0;


    PS::S32 c;
    while((c=getopt(argc,argv,"t:a:d:T:l:L:n:N:h")) != -1){
        switch(c){
        case 'a':
            delta_ax = atof(optarg);
            break;
        case 'd':
            dt = 1.0/atoi(optarg);
            break;
        case 't':
            theta = atof(optarg);
            break;
        case 'T':
            time_end = atof(optarg);
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            break;
        case 'L':
            n_loop_merge = atoi(optarg);
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            break;
        case 'N':
            n_glb = atol(optarg) * COMM_INFO.getNumberOfProc();
            break;
        case 'h':
            if(PS::Comm::getRank() == 0)
                PrintHelp();
            PS::Finalize();
            return 0;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
                PrintHelp();
            }
            PS::Abort();
        }
    }
    PS::F64 ax_in  = ax_cen - 0.5*delta_ax;
    PS::F64 ax_out = ax_cen + 0.5*delta_ax;
    PLANET.mass = 1.0;
    PLANET.pos = PLANET.vel = 0.0;
    DT_MERGE = dt*n_loop_merge;
    VEL_SHEAR = sqrt(PLANET.mass/ax_in) - sqrt(PLANET.mass/ax_out);
    

    PS::F64 tau = 1.0;
    PS::F64 ratio_r_phy_r_hill = 1.0;
    PS::F64 r_phy = sqrt(tau*(ax_out*ax_out - ax_in*ax_in) / n_glb);
#ifdef GRAV_ONLY
    FP_t::eps = r_phy;
#else
    FP_t::eps = 0.0;
#endif
    PS::F64 r_hill = r_phy / ratio_r_phy_r_hill;
    PS::F64 r_search = 6.0 * r_hill;
    //PS::F64 r_search = 12.0 * r_hill;
    //PS::F64 r_search = 24.0 * r_hill;
    FP_t::r_search = r_search;
    
    PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5;
    //PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5 * 0.01;
    
    PS::F64 dens = mass_dust * n_glb / (M_PI*(ax_out*ax_out - ax_in*ax_in));

    PS::F64 t_dur = dt * 64;
    //PS::F64 t_dur = dt * 128;
    PS::F64 e_refl = 0.5;
    PS::F64 ln_e_refl = std::log(e_refl);
    PS::F64 kappa     = std::pow((2.0*M_PI/t_dur),2.0); // k^{'}
    //PS::F64 eta       = 4.0*M_PI/(t_dur*std::sqrt(1.0+std::pow((M_PI/ln_e_refl),2.0))); // \eta^{'}
    PS::F64 eta       = 4.0*M_PI*ln_e_refl/(t_dur*std::sqrt(M_PI*M_PI+ln_e_refl)); // \eta^{'}
    FP_t::kappa = kappa;
    FP_t::eta   = eta;
    if(my_rank==0){
        std::cerr<<"kappa= "<<kappa<<" eta= "<<eta<<" r_phy= "<<r_phy<<std::endl;
    }

#if 0
    PS::F64 eng_disp_tmp = 0.0;
    PS::ReallocatableArray<FP_t> pi;
    pi.resizeNoInitialize(2);
    PS::ReallocatableArray<Force_t> fi;
    fi.resizeNoInitialize(2);
    for(auto i=0; i<2; i++) fi[i].clear();
    pi[0].mass = pi[1].mass = mass_dust;
    pi[0].pos = PS::F64vec(-r_phy*1.01, 0.0, 0.0);
    //pi[0].pos = PS::F64vec(-r_phy*0.99, 0.0, 0.0);
    pi[1].pos = -pi[0].pos;
    //pi[0].vel = PS::F64vec(sqrt(pi[0].mass/(2.0*r_phy))*0.1, 0.0, 0.0);
    pi[0].vel = PS::F64vec(r_phy, 0.0, 0.0);
    pi[1].vel = -pi[0].vel;
    pi[0].vel_full = pi[0].vel;
    pi[1].vel_full = pi[1].vel;
    pi[0].r_coll = pi[1].r_coll = r_phy;

    pi[0].id = 0;
    pi[1].id = 1;
    PS::ReallocatableArray<FP_t> pj;
    pj.resizeNoInitialize(2); 
    pj[0] = pi[0];
    pj[1] = pi[1];
    
    CalcForceEp<FP_t, FP_t, Force_t> func;
    for(int i=0; i<2; i++){fi[i].clear();}
    func(pi.getPointer(), 2, pj.getPointer(), 2, fi.getPointer());
    for(int i=0; i<2; i++){pi[i].copyFromForce(fi[i]);}
    
    PS::F64 e_kin0 = 0;
    PS::F64 e_pot0 = 0;
    for(int i=0; i<2; i++){
        e_kin0 += 0.5*pi[i].mass*pi[i].vel*pi[i].vel;
        e_pot0 += pi[i].mass*pi[i].pot;
    }
    PS::F64 e_tot0 = e_kin0 + e_pot0;
    if(my_rank==0){
        std::cout<<"e_kin0= "<<e_kin0
                 <<" e_pot0= "<<e_pot0
                 <<" e_tot0= "<<e_tot0
                 <<" pi[0].pos= "<<pi[0].pos.x
                 <<" vel= "<<pi[0].vel.x
                 <<" acc= "<<pi[0].acc.x
                 <<" acc_dash= "<<pi[0].acc_dash.x
                 <<std::endl;
    }
    //dt *= 0.125;
    dt *= 0.0625;
    dt *= 0.5;
    for(int loop=0; loop<200; loop++){
        PS::F64 va0 = 0.0;
        for(int i=0; i<2; i++){
            va0 += pi[i].mass*pi[i].acc_dash*pi[i].vel;
        }
        Kick(pi, 0.5*dt, 2);
        Drift(pi, dt, 2);
        Copy(pi, pj, 2);
        for(auto i=0; i<2; i++) fi[i].clear();
        func(pi.getPointer(), 2, pj.getPointer(), 2, fi.getPointer());
        for(int i=0; i<2; i++){pi[i].copyFromForce(fi[i]);}
        PS::F64 va1 = 0.0;
        Kick(pi, 0.5*dt, 2);
        for(int i=0; i<2; i++){
            va1 += pi[i].mass*pi[i].acc_dash*pi[i].vel;
        }
        eng_disp_tmp += (va0+va1)*dt*0.5;
        PS::F64 e_kin1 = 0;
        PS::F64 e_pot1 = 0;
        for(int i=0; i<2; i++){
            e_kin1 += 0.5*pi[i].mass*pi[i].vel*pi[i].vel;
            e_pot1 += pi[i].mass*pi[i].pot;
        }
        PS::F64 e_tot1 = e_kin1 + e_pot1;
        if(my_rank==0){
            std::cout<<"loop"<<loop
                //<<" e_kin1= "<<e_kin1
                //<<" e_pot1= "<<e_pot1
                     <<" e_tot1= "<<e_tot1
                //<<" de= "<<(e_tot1-e_tot0)/e_tot0
                     <<" de-e_disp= "<<(e_tot1-e_tot0-eng_disp_tmp)/e_tot0
                //<<" de+eng_disp= "<<(e_tot1-e_tot0+eng_disp)/e_tot0
                     <<" pi[0].pos= "<<pi[0].pos.x
                     <<" vel= "<<pi[0].vel.x
                     <<" acc= "<<pi[0].acc.x
                     <<" acc_dash= "<<pi[0].acc_dash.x
                     <<" eng_disp_tmp= "<<eng_disp_tmp
                     <<" va0= "<<va0
                     <<" va1= "<<va1
                     <<std::endl;
        }
    }
    PS::Finalize();
    return 0;
#endif
    
    PS::S32 nx, ny, nz;
    DivideNProc(nx, ny, nz, n_proc, delta_ax);

    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.setCommInfo(COMM_INFO);
    dinfo.setDomain(nx, ny, nz);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setPosRootDomain(PS::F64vec(-M_PI, -M_PI, -M_PI),
                           PS::F64vec( M_PI,  M_PI,  M_PI),
                           7u);

    PS::ParticleSystem<FP_t> system;
    system.initialize();
    system.setCommInfo(COMM_INFO);
    system.setAverageTargetNumberOfSampleParticlePerProcess(n_smp);
    PS::F64ort pos_domain;
    GetPosDomainCyl(delta_ax, pos_domain, nx, ny);
    SetParticleKeplerDiskCyl(system, n_glb, ax_in, ax_out, dens, pos_domain, r_phy, true);
    //SetParticleKeplerDiskCyl(system, n_glb, ax_in, ax_out, dens, pos_domain, r_phy, false);
    n_loc = system.getNumberOfParticleLocal();
    n_glb = system.getNumberOfParticleGlobal();
    std::cerr<<"my_rank= "<<my_rank<<" n_loc= "<<n_loc<<" n_glb= "<<n_glb
             <<" mass_dust= "<<mass_dust<<" mass_dust*n_glb= "<<mass_dust*n_glb
             <<std::endl;
    
    for(auto i=0; i<n_loc; i++){
        assert(dinfo.getPosRootDomain().contained(system[i].getPos()));
    }
    dinfo.decomposeDomainAll(system);
    
    if(my_rank==0){
        std::cerr<<"dinfo.getPosRootDomain()= "<<dinfo.getPosRootDomain()<<std::endl;
        for(auto i=0; i<n_proc; i++){
            std::cerr<<"i= "<<i<<" dinfo.getPosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }
    
    system.adjustPositionIntoRootDomain(dinfo);
    
    n_loc = system.getNumberOfParticleLocal();
    
    system.exchangeParticle(dinfo);
    n_loc = system.getNumberOfParticleLocal();
    SetID(system);


#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FP_t::eps);
#endif
    
    Tree_t tree;
    tree.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    tree.setCommInfo(COMM_INFO);
    MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);


#ifdef FORCE_CHECK
    PS::ReallocatableArray<Force_t> force_mono(n_loc, n_loc, PS::MemoryAllocMode::Default);
    MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSpMono<EPI_t, SPJ_t, Force_t>(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
    for(auto i=0; i<n_loc; i++){
        force_mono[i].acc = system[i].acc;
        force_mono[i].pot = system[i].pot;
    }

    PS::ReallocatableArray<Force_t> force_quad(n_loc, n_loc, PS::MemoryAllocMode::Default);
    MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSpQuad<EPI_t, SPJ_t, Force_t>(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
    for(auto i=0; i<n_loc; i++){
        force_quad[i].acc = system[i].acc;
        force_quad[i].pot = system[i].pot;
    }
    PS::ReallocatableArray<Force_t> force_dir(n_loc, n_loc, PS::MemoryAllocMode::Default);
    system.calcForceDirectParallel(CalcForceEp<EPI_t, EPJ_t, Force_t>(), force_dir.getPointer(), dinfo, true);
    
    if(PS::Comm::getRank()==0){
        for(auto i=0; i<n_loc; i++){
            std::cout<<"force_dir[i].acc= "<<force_dir[i].acc
                     <<"force_mono[i].acc= "<<force_mono[i].acc<<" ( "<<force_mono[i].acc-force_dir[i].acc
                     <<" ) force_quad[i].acc= "<<force_quad[i].acc<<" ( "<<force_quad[i].acc-force_dir[i].acc<<" ) "
                     <<std::endl;            
        }
    }
    PS::Finalize();
    return 0;
#endif
    
    CalcForceFromPlanet(system, PLANET);

    Energy eng_init;
    Force_t * force_dir = new Force_t[n_loc];
    //system.calcForceDirectParallel(CalcForceEp<EPI_t, EPJ_t, Force_t>(), &force_dir[0], dinfo, true);
    //eng_init.calc(&system[0], force_dir, PLANET, n_loc);
    eng_init.calc(system);

    std::cerr<<"my_rank= "<<my_rank<<" tree.getUsedMemorySize()= "<<tree.getUsedMemorySize()<<std::endl;
    std::cerr<<" tree.dumpMemSizeUsed():"<<std::endl;
    tree.dumpMemSizeUsed(std::cerr);
    std::cerr<<" PS::MemeoryPool::dump():"<<std::endl;
    PS::MemoryPool::dump(std::cerr);


    //Rotate(system, dt*n_loop_merge/2);
    //exit(1);
    
    PS::F64 * force_diff_00 = new PS::F64[n_loc];
    PS::F64 * force_diff_31 = new PS::F64[n_loc];
    PS::F64 * force_diff_63 = new PS::F64[n_loc];
    PS::S64 n_loop = 0;
    PS::F64 eng_disp = 0.0;
    /*
    std::ofstream * fout_err = new std::ofstream[n_loc];
    for(auto i=0; i<n_loc; i++){
        std::string name = "./result/err" + std::to_string(i);
        fout_err[i].open(name);
    }
    */
    while(time_sys <= time_end){
        COMM_INFO.barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"n_loop= "<<n_loop<<std::endl;

        n_loc = system.getNumberOfParticleLocal();
        PS::F64 va0_loc = 0.0;
        for(int i=0; i<n_loc; i++){
            va0_loc += system[i].mass*system[i].acc_dash*system[i].vel;
        }

        Kick(system, 0.5*dt);
        Drift(system, dt);

        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"after drift"<<std::endl;
        
        if(n_loop % n_loop_merge == 0){
#ifdef ROTATE
            Rotate(system, dt*n_loop_merge*0.5);
#endif
            dinfo.decomposeDomainAll(system);
            system.adjustPositionIntoRootDomain(dinfo);
            system.exchangeParticle(dinfo);
            n_loc = system.getNumberOfParticleLocal();
            n_glb = system.getNumberOfParticleGlobal();
        }

        //system.adjustPositionIntoRootDomain(dinfo);
        
        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"after exchange"<<std::endl;
        
        /*
        for(auto i=0; i<n_loc; i++){
            if(!dinfo.getPosDomain(my_rank).contained(system[i].getPos())){
                std::cerr<<"dinfo.getPosDomain(my_rank)= "<<dinfo.getPosDomain(my_rank)
                         <<" system[i].getPos()= "<<system[i].getPos()
                         <<std::endl;
            }
            assert(dinfo.getPosDomain(my_rank).contained(system[i].getPos()));
        }
        */
        tree.clearCounterAll();
        if(n_loop % n_loop_merge == 0){
#ifdef ROTATE
            MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
            Rotate(system, -dt*n_loop_merge*0.5);
            MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::REUSE_LIST);
#else
            MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
#endif
        }
        else{
            MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::REUSE_LIST);
        }


        //if( n_loop==0  || n_loop==31 || n_loop==63){
        if(0){
            //Force_t * force_dir = new Force_t[n_loc];
            PS::F64 * force_diff = new PS::F64[n_loc];
            system.calcForceDirectParallel(CalcForceEp<EPI_t, EPJ_t, Force_t>(), &force_dir[0], dinfo, true);
            PS::F64 diff_max = 0.0;
            PS::S32 adr_max = 0;
            for(auto i=0; i<n_loc; i++){
                force_diff[i] = sqrt((system[i].acc - force_dir[i].acc)*(system[i].acc - force_dir[i].acc));
                PS::F64 a2 = (force_dir[i].acc - system[i].acc)*(force_dir[i].acc - system[i].acc);
                if(a2 > diff_max){
                    diff_max = a2;
                    adr_max = i;
                }
            }
            //for(auto i=0; i<n_loc; i++)
            //    fout_err[i]<<time_sys<<"  "<<force_diff[i]<<"  "<<system[i].getPos()<<std::endl;
            
            std::ofstream fout;
            if(n_loop==0){
                fout.open("./file00.dat");
                for(auto i=0; i<n_loc; i++)
                    force_diff_00[i] = force_diff[i];
            }
            else if(n_loop==31){
                fout.open("./file31.dat");
                for(auto i=0; i<n_loc; i++)
                    force_diff_31[i] = force_diff[i];
            }
            else if(n_loop==63){
                fout.open("./file63.dat");
                for(auto i=0; i<n_loc; i++)
                    force_diff_63[i] = force_diff[i];
            }
            std::sort(force_diff, force_diff+n_loc);
            for(auto i=0; i<n_loc; i++){
                fout<<(double)(i+1)/n_loc<<"  "<<force_diff[i]<<std::endl;
            }
            fout.close();
            std::cerr<<"adr_max= "<<adr_max
                     <<" system[adr_max].acc= "<<system[adr_max].acc
                     <<" force_dir[adr_max].acc= "<<force_dir[adr_max].acc
                     <<std::endl;
            std::cerr<<"diff= "<<system[adr_max].acc-force_dir[adr_max].acc
                     <<" |a_diff|= "<<sqrt(diff_max)
                     <<std::endl;
            //delete [] force_dir;
            delete [] force_diff;
        }
        /*
        if(n_loop==63){
            std::ofstream fout;
            fout.open("./diff_acc.dat");
            for(auto i=0; i<n_loc; i++){
                fout<<i<<"  "<<force_diff_00[i]<<"  "<<force_diff_31[i]<<"  "<<force_diff_63[i]<<"  "<<system[i].pos<<std::endl;
            }
            //return 0;
        }
        */
        
        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0) std::cerr<<"after force"<<std::endl;


        
        PS::MemoryPool::checkEmpty();
        CalcForceFromPlanet(system, PLANET);
    
        Kick(system, 0.5*dt);

        n_loc = system.getNumberOfParticleLocal();
        PS::F64 va1_loc = 0.0;
        for(int i=0; i<n_loc; i++){
            va1_loc += system[i].mass*system[i].acc_dash*system[i].vel;
        }
        eng_disp += (va0_loc+va1_loc)*dt*0.5;
        PS::F64 eng_disp_glb = COMM_INFO.getSum(eng_disp);
        Energy eng_now;
        //system.calcForceDirectParallel(CalcForceEp<EPI_t, EPJ_t, Force_t>(), &force_dir[0], dinfo, true);
        //eng_now.calc(&system[0], force_dir, PLANET, n_loc);
        eng_now.calc(system);
        
        auto n_walk = tree.getNumberOfWalkGlobal();
        auto n_int_ep_ep = tree.getNumberOfInteractionEPEPGlobal();
        auto n_int_ep_sp = tree.getNumberOfInteractionEPSPGlobal();
        if(PS::Comm::getRank()==1){
            eng_init.dump(std::cout);
            eng_now.dump(std::cout);
            std::cout<<"eng_disp_glb= "<<eng_disp_glb<<std::endl;
            std::cout<<"time_sys= "<<time_sys<<" n_loop= "<<n_loop<<" (eng_now.tot-eng_init.tot)/eng_init.tot= "<<(eng_now.tot-eng_init.tot)/eng_init.tot
                     <<" (eng_now.tot-eng_init.tot-eng_disp)/eng_init.tot= "<<(eng_now.tot-eng_init.tot-eng_disp)/eng_init.tot
                     <<std::endl;
            std::cout<<"n_int_ep_ep= "<<n_int_ep_ep
                     <<" n_int_ep_sp= "<<n_int_ep_sp
                     <<" <n_epi>= "<<((PS::F64)n_glb) / n_walk
                     <<" <n_epj>= "<<((PS::F64)n_int_ep_ep) / n_glb
                     <<" <n_spj>= "<<((PS::F64)n_int_ep_sp) / n_glb
                     <<std::endl;
        }
        tree.clearCounterAll();
        PS::S64 size_used_loc = tree.getUsedMemorySize();
        PS::S64 size_used_glb = COMM_INFO.getMaxValue(size_used_loc);
        if(size_used_loc == size_used_glb){
            std::cerr<<"my_rank= "<<my_rank<<" tree.getUsedMemorySize()= "<<tree.getUsedMemorySize()<<std::endl;
            std::cerr<<" tree.dumpMemSizeUsed():"<<std::endl;
            tree.dumpMemSizeUsed(std::cerr);
        }
        
        time_sys += dt;
        n_loop++;
    }

#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif
    PS::Finalize();
    return 0;
}
