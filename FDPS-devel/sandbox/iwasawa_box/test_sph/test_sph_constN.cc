#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>

//#define ONE_DIM_SIM
//#define TC_CUSP_KERNEL
//#define TWO_DIM_SIM
#define THREE_DIM_SIM
//#define DEBUG_PRINT
//#define DEBUG

class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp){
		fscanf(fp, "%lf\n", &time);
		fscanf(fp, "%lld\n", &n_body);
		return n_body;
    }
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%lld\n", n_body);
	}
};

///////////
/// SPH ///
PS::F64 CubicSpline(const PS::F64 r_sq,
                    const PS::F64 h_inv){
    PS::F64 xi = sqrt(r_sq) * h_inv;
    PS::F64 xi10 = (1.0-xi > 0.0) ? 1.0-xi : 0.0;
    PS::F64 xi05 = (0.5-xi > 0.0) ? 0.5-xi : 0.0;
    return xi10*xi10*xi10 - 4.0*xi05*xi05*xi05;
}

PS::F64vec CubicSpline_ri(const PS::F64vec & rij,
                          const PS::F64 h_inv){
    PS::F64 r_sq = rij * rij;
    PS::F64 r = sqrt(r_sq);
    PS::F64 xi = r * h_inv;
#ifdef TC_CUSP_KERNEL
    xi = std::max(xi, 1.0/3.0);
#endif
    PS::F64 xi10 = (1.0-xi > 0.0) ? 1.0-xi : 0.0;
    PS::F64 xi05 = (0.5-xi > 0.0) ? 0.5-xi : 0.0;
    PS::F64 C = (-3.0*xi10*xi10 + 12.0*xi05*xi05) * h_inv;
    return C * rij / r;
}

class ResultDens{
public:
    PS::F64 dens;
    PS::F64 divv;
    PS::F64 r_kernel;
    PS::F64 r_search;
    PS::S32 n_ngb;
    void clear(){
        dens = divv = r_kernel = r_search = 0.0;
        n_ngb = 0;
    }
};

class ResultForce{
public:
    PS::F64 eng_dot;
    PS::F64 dt;
    PS::F64vec acc;
    void clear(){
        acc = 0.0;
        eng_dot = dt = 0.0;
    }
};

class FullPtcl{
public:
    void writeAscii(FILE* fp) const{
        fprintf(fp, 
                "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, 
                this->vel.x, this->vel.y, this->vel.z, 
                this->acc.x, this->acc.y, this->acc.z,
                this->dens, this->pres, this->eng,
                this->eng_dot, this->vel_sound);
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F64vec & pos_new ) {
        pos = pos_new;
    }
    PS::F64 getRSearch() const {
        return this->r_search;
    }
    void copyFromForce(const ResultDens & rd){
        dens = rd.dens;
        divv = rd.divv;
        r_kernel = rd.r_kernel;
        r_search = rd.r_search;
        n_ngb = rd.n_ngb;
    }
    void copyFromForce(const ResultForce & rf){
        acc = rf.acc;
        eng_dot = rf.eng_dot;
        dt = rf.dt;
    }

    PS::F64 mass;
    PS::F64vec pos;
    PS::F64 r_kernel; // H
    PS::F64 r_search;
    PS::F64 dens;
    PS::F64 divv;
    PS::F64 pres;
    PS::F64 vel_sound;
    PS::F64vec acc;
    PS::F64 eng_dot;
    PS::S64 id;
    PS::F64vec vel;
    PS::F64 eng;
    PS::F64vec vel_half;
    PS::F64 eng_half;
    PS::S64 n_ngb;
    PS::F64 dt;
};

class EPIDens{
public:
#ifdef ONE_DIM_SIM
    static const size_t n_neighbour_crit = 8;
    //static const size_t n_neighbour_crit = 12;
#elif defined(TWO_DIM_SIM)
    static const size_t n_neighbour_crit = 32;
#elif defined(THREE_DIM_SIM)
    static const size_t n_neighbour_crit = 60;
#endif
    PS::S64 id;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 r_search;
    PS::F64 getRSearch() const { return r_search; }
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FullPtcl & fp){
        id = fp.id;
        pos = fp.getPos();
        vel = fp.vel;
        r_search = fp.r_search;
    }
};

class EPJDens{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64 r_search;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 getCharge() const { return mass; }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec pos_new){ pos = pos_new; }
    PS::F64 getRSearch() const { return r_search; }
    void copyFromFP(const FullPtcl & fp){
        id = fp.id;
        mass = fp.getCharge();
        r_search = fp.getRSearch();
        pos = fp.getPos();
        vel = fp.vel;
    }
};

struct CalcDensEpEp{
    void operator () (const EPIDens * ep_i,
                      const PS::S32 n_ip,
                      const EPJDens * ep_j,
                      const PS::S32 n_jp,
                      ResultDens * dens){
#ifdef ONE_DIM_SIM
        static const PS::F64 Cnorm = 8.0 / 3.0; // 1D
#elif defined(TWO_DIM_SIM)
        static const PS::F64 Cnorm = 80.0 / (7.0*M_PI);
#elif defined(THREE_DIM_SIM)
        static PS::F64 Cnorm = 16.0 / M_PI; // 3D
#endif
        //static const PS::F64 Cnorm = 8.0 / 3.0 * h_inv; // 1D
        //static PS::F64 Cnorm = 80.0/(7.0*M_PI) * h_inv * h_inv; // 2D
        //static PS::F64 Cnorm = 16.0/M_PI * h_inv * h_inv * h_inv; // 3D

        std::vector < std::pair < PS::F64, PS::S64 > > rsq_id;
        rsq_id.reserve(EPIDens::n_neighbour_crit*3+100);
        for(PS::S32 i=0; i<n_ip; i++){
            if (ep_i[i].getRSearch() == 0.0) continue;
            dens[i].clear();
            rsq_id.clear();
            const PS::F64 r_search_sq = ep_i[i].getRSearch() * ep_i[i].getRSearch();

#ifdef DEBUG
            std::vector<int> j_id;
            j_id.reserve(20);
            for(PS::S32 j=0; j<n_jp; j++) j_id.push_back(ep_j[j].id);
            std::sort(j_id.begin(), j_id.end());
            for(int jj=0; jj<n_jp-1; jj++){
                if(j_id[jj]+1 != j_id[jj+1]){
                    std::cerr<<"jj="<<jj<<std::endl;
                    std::cerr<<"ep_i[i].id="<<ep_i[i].id<<std::endl;
                    std::cerr<<"j_id[jj]="<<j_id[jj]<<std::endl;
                    std::cerr<<"j_id[jj+1]="<<j_id[jj+1]<<std::endl;
                    std::cerr<<std::endl;
                }
            }
#endif

            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec rij = ep_i[i].getPos() - ep_j[j].getPos();
                const PS::F64 r_sq = rij * rij;
                if(r_sq < r_search_sq) rsq_id.push_back( std::make_pair(r_sq, j) );
            }


            const size_t ncrit = EPIDens::n_neighbour_crit;
            dens[i].n_ngb = rsq_id.size();
            if( rsq_id.size() <=  ncrit ){
                // search radius was too small !!
                dens[i].r_search = ep_i[i].getRSearch() * cbrt(2.0);
                dens[i].r_kernel = -1.0;
                continue;
            }
            else{
                // DON'T need to search next
                dens[i].r_search = 0.0;
                std::sort(rsq_id.begin(), rsq_id.end());
                //dens[i].r_kernel = sqrt( rsq_id[ncrit-1].first );
                dens[i].r_kernel = sqrt( rsq_id[ncrit-1].first ) * 1.2;
                PS::F64 dens_tmp = 0.0;
                PS::F64 h_inv = 1.0/dens[i].r_kernel;
                for(PS::S32 k=0; k<PS::S32(ncrit); k++){
                    const PS::S32 id = rsq_id[k].second;
                    const PS::F64 mj = ep_j[id].getCharge();
                    //dens_tmp += mj * CubicSpline(rsq_id[id].first, h_inv);
                    dens_tmp += mj * CubicSpline(rsq_id[k].first, h_inv);
                }
#ifdef ONE_DIM_SIM
                dens[i].dens = Cnorm * h_inv * dens_tmp;
#elif  defined(TWO_DIM_SIM)
                dens[i].dens = Cnorm * h_inv * h_inv * dens_tmp;
#elif  defined(THREE_DIM_SIM)
                dens[i].dens = Cnorm * h_inv * h_inv * h_inv * dens_tmp;
#endif
                PS::F64 divv_tmp = 0.0;
                for(PS::S32 k=0; k<PS::S32(ncrit); k++){
                    const PS::S32 id = rsq_id[k].second;
                    const PS::F64 mj = ep_j[id].getCharge();
                    if(ep_i[i].id == ep_j[id].id) continue;
                    const PS::F64vec rij = ep_i[i].pos - ep_j[id].pos;
                    const PS::F64vec vij = ep_i[i].vel - ep_j[id].vel;
                    divv_tmp -= mj * vij * CubicSpline_ri(rij, h_inv);
                }
#ifdef ONE_DIM_SIM
                dens[i].divv = divv_tmp * Cnorm * h_inv / dens[i].dens;
#elif  defined(TWO_DIM_SIM)
                dens[i].divv = divv_tmp * Cnorm * h_inv * h_inv / dens[i].dens;
#elif  defined(THREE_DIM_SIM)
                dens[i].divv = divv_tmp * Cnorm * h_inv * h_inv * h_inv / dens[i].dens;
#endif
            }
        }
    }
};


template<class Tpsys>
void CalcPressureAndSoundVelocity(Tpsys & system, const PS::F32 gamma){
    const PS::F32 C = gamma - 1.0;
    const PS::S32 n_ptcl = system.getNumberOfParticleLocal();
    for(size_t i=0; i<n_ptcl; i++){
        system[i].pres = C * system[i].dens * system[i].eng;
        system[i].vel_sound = sqrt(gamma * system[i].pres / system[i].dens);
    }
}

template<class Tpsys>
void CalcEnergyAndSoundVelocity(Tpsys & system, const PS::F32 gamma){
    const PS::F32 C = gamma - 1.0;
    const PS::S32 n_ptcl = system.getNumberOfParticleLocal();
    for(size_t i=0; i<n_ptcl; i++){
        system[i].eng = system[i].pres / (C * system[i].dens);
        system[i].vel_sound = sqrt(gamma * system[i].pres / system[i].dens);
    }
}


class EPIForce{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 divv;
    PS::F64 vel_sound;
    PS::F64 dens;
    PS::F64 pres;
    PS::F64 r_kernel;
    PS::S64 id;
    PS::F64 getRSearch() const {return r_kernel;}
    void copyFromFP(const FullPtcl & fp){
        pos = fp.pos;
        vel = fp.vel;
        divv = fp.divv;
        vel_sound = fp.vel_sound;
        dens = fp.dens;
        pres = fp.pres;
        r_kernel = fp.r_kernel;
        id = fp.id;
    }
    PS::F64vec getPos() const{
        return pos;
    }
};

class EPJForce{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 divv;
    PS::F64 vel_sound;
    PS::F64 pres;
    PS::F64 dens;
    PS::F64 r_kernel;
    PS::F64 mass;
    PS::S64 id;
    PS::F64vec getPos() const {return pos;}
    void setPos (const PS::F64vec & pos_new) {pos = pos_new;}
    PS::F64 getRSearch() const {return r_kernel;}
    void copyFromFP(const FullPtcl & fp){
        pos = fp.pos;
        vel = fp.vel;
        divv = fp.divv;
        vel_sound = fp.vel_sound;
        pres = fp.pres;
        dens = fp.dens;
        r_kernel = fp.r_kernel;
        mass = fp.mass;
        id = fp.id;
    }
};

struct CalcForceEpEp{
    void operator () (const EPIForce * ep_i,
                      const PS::S32   n_ip,
                      const EPJForce * ep_j,
                      const PS::S32   n_jp,
                      ResultForce    * force){
        static const PS::F64 alpha = 1.0;
        //static const PS::F64 alpha = 2.0;
        //static const PS::F64 alpha = 3.0;
        //static const PS::F64 alpha = 5.0;
        static const PS::F64 beta = 2.0 * alpha;
        static const PS::F64 Ccfl = 0.1;
#ifdef ONE_DIM_SIM
        static const PS::F64 Cnorm = 8.0 / 3.0; // 1D
#elif defined(TWO_DIM_SIM)
        static const PS::F64 Cnorm = 80.0 / (7.0*M_PI); // 2D
#elif defined(THREE_DIM_SIM)
        static const PS::F64 Cnorm = 16.0 / M_PI; // 3D
#endif
        //static const PS::F64 Cnorm = 8.0 / 3.0 * h_inv; // 1D
        //static PS::F64 Cnorm = 80.0/(7.0*M_PI) * h_inv * h_inv; // 2D
        //static PS::F64 Cnorm = 16.0/M_PI * h_inv * h_inv * h_inv; // 3D
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64 cs_i = ep_i[i].vel_sound;
            const PS::F64 dens_i = ep_i[i].dens;
            const PS::F64 p_dens2_i = ep_i[i].pres / (dens_i * dens_i);
            //std::cerr<<"ep_i[i].r_kernel="<<ep_i[i].r_kernel<<std::endl;
            const PS::F64 h_inv_i = 1.0 / ep_i[i].r_kernel;
            PS::F64vec acc = 0.0;
            PS::F64 eng_dot = 0.0;
            PS::F64 mu_ij_max = 0.0;
            PS::F64 vel_sig_max = 0.1;
            for(PS::S32 j=0; j<n_jp; j++){
                //std::cerr<<"ep_i[i].id="<<ep_i[i].id<<" ep_j[j].id="<<ep_j[j].id<<std::endl;
                if(ep_i[i].id == ep_j[j].id) continue;
                const PS::F64 mj = ep_j[j].mass;
                const PS::F64vec r_ij = ep_i[i].pos - ep_j[j].pos;
                const PS::F64vec v_ij = ep_i[i].vel - ep_j[j].vel;
                const PS::F64 r_sq = r_ij * r_ij;
                const PS::F64 r_crit = std::max(ep_i[i].r_kernel, ep_j[j].r_kernel);
                //std::cerr<<"r_crit*r_crit="<<r_crit*r_crit<<std::endl;
                //std::cerr<<"r_sq="<<r_sq<<std::endl;
                if(r_sq > r_crit * r_crit) continue;
                //std::cerr<<"ep_j[j].r_kernel="<<ep_j[j].r_kernel<<std::endl;
                const PS::F64 h_inv_j = 1.0 / ep_j[j].r_kernel;
                const PS::F64 h_ij = 0.5 * (ep_i[i].r_kernel + ep_j[j].r_kernel);
                const PS::F64 cs_j = ep_j[j].vel_sound;
                const PS::F64 dens_j = ep_j[j].dens;
                const PS::F64 p_dens2_j = ep_j[j].pres / (dens_j * dens_j);
                const PS::F64 rijvij = r_ij * v_ij;
                const PS::F64 dens_ij = 0.5 * (dens_i + dens_j);
                PS::F64 mu_ij = rijvij * h_ij / (r_ij*r_ij + 0.01*h_ij*h_ij);
                //PS::F64 mu_ij = rijvij * h_ij / (r_ij*r_ij + 0.0001*h_ij*h_ij);
#if 0
                PS::F64 pi_ij = (-alpha * 0.5*(cs_i + cs_j) * mu_ij + beta * mu_ij * mu_ij)/ dens_ij;
#else
                //const PS::F64 w_ij =  (rijvij < 0.0) ? rijvij * r_inv : 0.0;
                //std::cerr<<"sqrt(r_sq)="<<sqrt(r_sq)<<std::endl;
                const PS::F64 w_ij =  rijvij / sqrt(r_sq);
                PS::F64 vel_sig = 0.5 * (cs_i + cs_j - 3.0*w_ij);
                PS::F64 pi_ij = -alpha * vel_sig * w_ij / dens_ij;
                vel_sig_max = std::max(vel_sig, vel_sig_max);
#endif
                pi_ij = ( rijvij <= 0.0 ) ? pi_ij : 0.0;
                const PS::F64 h_inv_ij = 1.0 / h_ij;
#ifdef ONE_DIM_SIM
                const PS::F64vec grad_w_i = Cnorm  * h_inv_i * CubicSpline_ri(r_ij, h_inv_i);
                const PS::F64vec grad_w_j = Cnorm  * h_inv_j * CubicSpline_ri(r_ij, h_inv_j);
                //const PS::F64vec grad_w_ij = Cnorm * h_inv_ij * CubicSpline_ri(r_ij, h_inv_ij);
                const PS::F64vec grad_w_ij =  0.5 * (grad_w_i +  grad_w_j);
#elif defined(TWO_DIM_SIM)
                const PS::F64vec grad_w_i = Cnorm  * h_inv_i  * h_inv_i  * CubicSpline_ri(r_ij, h_inv_i);
                const PS::F64vec grad_w_j = Cnorm  * h_inv_j  * h_inv_j  * CubicSpline_ri(r_ij, h_inv_j);
                const PS::F64vec grad_w_ij = Cnorm * h_inv_ij * h_inv_ij * CubicSpline_ri(r_ij, h_inv_ij);
#elif defined(THREE_DIM_SIM)
                const PS::F64vec grad_w_i = Cnorm  * h_inv_i  * h_inv_i  * h_inv_i  * CubicSpline_ri(r_ij, h_inv_i);
                const PS::F64vec grad_w_j = Cnorm  * h_inv_j  * h_inv_j  * h_inv_j  * CubicSpline_ri(r_ij, h_inv_j);
                const PS::F64vec grad_w_ij = Cnorm * h_inv_ij * h_inv_ij * h_inv_ij * CubicSpline_ri(r_ij, h_inv_ij);
                //std::cerr<<"grad_w_ij="<<grad_w_ij<<std::endl;
#endif
                acc -= mj * (p_dens2_i*grad_w_i + p_dens2_j*grad_w_j + pi_ij * grad_w_ij);
                eng_dot += mj * (p_dens2_i * v_ij * grad_w_i + 0.5 * pi_ij * v_ij * grad_w_ij);
                mu_ij_max = std::max(mu_ij, mu_ij_max);
            }
            force[i].acc = acc;
            force[i].eng_dot = eng_dot;
#if 0
            const PS::F64 dt_cv = ep_i[i].r_kernel / (cs_i + 0.6*(cs_i + 2.0*mu_ij_max));
            const PS::F64 tmp = 1.0 / sqrt(acc*acc);
            const PS::F64 dt_f = sqrt( ep_i[i].r_kernel * tmp );
            force[i].dt = Ccfl * std::min(dt_cv, dt_f);
#else
            force[i].dt = Ccfl * ep_i[i].r_kernel / vel_sig_max;
            //if(vel_sig_max >= 0.0) force[i].dt = Ccfl * ep_i[i].r_kernel / vel_sig_max;
            //else force[i].dt = 0.1;
#endif
        }
    }
};


template<class Tpsys>
void Kick1( Tpsys & system,
            const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        system[i].vel_half  = system[i].vel + system[i].acc * dt;
        system[i].eng_half  = system[i].eng + system[i].eng_dot * dt;
    }
}

template<class Tpsys>
void Drift( Tpsys & system,
            const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        system[i].pos  += system[i].vel_half * dt; // corrected value
        system[i].vel  += system[i].acc * dt;           // predicted value
        system[i].eng  += system[i].eng_dot * dt;        // predicted value
    }
}

template<class Tpsys>
void Kick2( Tpsys & system,
            const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        system[i].vel  = system[i].vel_half + system[i].acc * dt;
        system[i].eng  = system[i].eng_half + system[i].eng_dot * dt;
    }
}

template<class Tpsys>
PS::F64 CalcDt(Tpsys & system){
    //std::cerr<<"PS::Comm::getRank()="<<PS::Comm::getRank()<<std::endl;
    PS::S32 n = system.getNumberOfParticleLocal();
    //std::cerr<<"n="<<n<<std::endl;
    PS::F64 dt_min_loc = system[0].dt;
    //std::cerr<<"dt_min_loc="<<dt_min_loc<<std::endl;
    for(PS::S32 i=1; i<n; i++){
        dt_min_loc = std::min(dt_min_loc, system[i].dt);
    }
    //std::cerr<<"dt_min_loc2="<<dt_min_loc<<std::endl;
    return PS::Comm::getMinValue(dt_min_loc);
}

#ifdef ONE_DIM_SIM
template<class Tpsys>
void WriteSPHFormat(const Tpsys & psys,
                    const PS::F32 time_sys,
                    const PS::S32 snp_id,
                    const char * dir_name){
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::S32 n_glb = 0;
    FullPtcl * fp;
    PS::AllGatherParticle(fp, n_glb, &psys[0], n_loc);
    if(PS::Comm::getRank () == 0){
        const PS::S32 STRINGSIZE = 1024;
        char sout[STRINGSIZE];
        sprintf(sout,"%s/snap%5d.dat", dir_name, snp_id);
        for(int i=0;i<STRINGSIZE;i++)if(sout[i]==' ')sout[i]='0';
        std::ofstream foutput;
        foutput.open(sout);
        foutput<<std::setprecision(15);
        foutput<<n_glb<<std::endl;
        foutput<<"3"<<std::endl;
        foutput<<time_sys<<std::endl;
        for(PS::S32 i=0; i<n_glb; i++){
            foutput<<fp[i].pos.x<<"   "<<fp[i].dens<<"   "<<fp[i].pres<<"   "<<fp[i].vel.x<<"   "
                   <<fp[i].eng<<"   "<<fp[i].eng_dot<<"   "<<fp[i].vel_sound<<"   "<<fp[i].acc.x
                   <<std::endl;
        }
        foutput.close();
    }
}

#else
template<class Tpsys>
void WriteSPHFormat(const Tpsys & psys,
                    const PS::F32 time_sys,
                    const PS::S32 snp_id,
                    const char * dir_name){
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::S32 n_glb = 0;
    FullPtcl * fp;
    PS::AllGatherParticle(fp, n_glb, &psys[0], n_loc);
    if(PS::Comm::getRank () == 0){
        const PS::S32 STRINGSIZE = 1024;
        char sout[STRINGSIZE];
        sprintf(sout,"%s/snap%5d.dat", dir_name, snp_id);
        for(int i=0;i<STRINGSIZE;i++)if(sout[i]==' ')sout[i]='0';
        std::ofstream foutput;
        foutput.open(sout);
        foutput<<std::setprecision(15);
        foutput<<n_glb<<std::endl;
        foutput<<"3"<<std::endl;
        foutput<<time_sys<<std::endl;
        for(PS::S32 i=0; i<n_glb; i++){
            foutput<<fp[i].pos<<"   "<<fp[i].vel<<"   "<<fp[i].acc<<"   "<<fp[i].dens<<"   "<<fp[i].pres<<"   "
                   <<fp[i].eng<<"   "<<fp[i].eng_dot<<"   "<<fp[i].vel_sound<<std::endl;
        }
        foutput.close();
    }
}
#endif

template<class Tpsys>
void SetParticlesShockTube(Tpsys & psys,
                           const PS::S32 n_glb,                           
                           PS::S32 & n_loc,
                           PS::F64 & time_sys){
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    time_sys = 0.0;
    PS::S32 n_left = n_glb*0.8;
    PS::F64 m_eq = 1.0 / n_glb;
    PS::F64 rho_left = 1.0;
    PS::F64 dx_left = m_eq / rho_left;
    PS::F64 dx_right = dx_left*4;
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        psys[i].mass = 1.0 / n_glb;
        psys[i].pos.y = psys[i].pos.z = 0.0;
        psys[i].vel = 0.0;
        //psys[i].r_search = dx_right * 3.0;
        psys[i].r_search = dx_right * 1.0;  // too short
        psys[i].id = i_h + i;
        if( (i_h+i) < n_left ){
            psys[i].pos.x = (i_h+i) * dx_left;
            psys[i].eng = 2.5;
        }
        else{
            psys[i].pos.x = (n_left-1)*dx_left + (i_h+i-(n_left-1)) * dx_right;
            psys[i].eng = 1.795;
        }
    }
}

template<class Tpsys>
void SetParticlesStrongShock(Tpsys & psys,
                             const PS::S32 n_glb,
                             PS::S32 & n_loc,
                             PS::F64 & time_sys){
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc;
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    time_sys = 0.0;
    PS::F64 L = 2.0;
    PS::F64 rho = 1.0;
    PS::F64 gamma = 1.4;
    PS::F64 pres_left = 1000.0;
    PS::F64 pres_right = 0.01;
    PS::F64 dl = L / n_glb;
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        psys[i].mass = rho * L / n_glb;
        psys[i].pos.x = (i_h+i) * dl - 1.0;
        psys[i].pos.y = psys[i].pos.z = 0.0;
        psys[i].vel = 0.0;
        psys[i].r_search = dl * 3.0;
        psys[i].id = i_h + i;
        if( psys[i].pos.x < 0.0 ){
            psys[i].pres = pres_left;
            psys[i].eng = pres_left / ((gamma - 1.0) * rho);
        }
        else{
            psys[i].pres = pres_right;
            psys[i].eng = pres_right / ((gamma - 1.0) * rho);
        }
    }
}


template<class Tpsys>
void SetParticlesContact2D(Tpsys & psys,
                           const PS::S32 n_glb,                           
                           PS::S32 & n_loc,
                           PS::F64 & time_sys){
    const PS::S32 n_one_dim = sqrt( (double)n_glb * 1.00001);
    assert(n_one_dim*n_one_dim == n_glb);
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::F64 dl = 1.0 / n_one_dim;
    const PS::F64 pres = 2.5;
    const PS::F64 gamma = 5.0 / 3.0;
    const PS::F64 rho_in = 4.0;
    const PS::F64 rho_out = 1.0;
    const PS::F64 low = 0.25;
    const PS::F64 high = 0.75;
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    time_sys = 0.0;
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        PS::S32 ix = (i_h+i) / n_one_dim;
        PS::S32 iy = (i_h+i) % n_one_dim;
        psys[i].id = i_h + i;
        psys[i].mass = 1.0 / n_glb;
        psys[i].vel = 0.0;
        psys[i].r_search = dl * 2.0;
        psys[i].pos.x = ix * dl;
        psys[i].pos.y = iy * dl;
        psys[i].pos.z = 0.0;
        if( low<=psys[i].pos.x && psys[i].pos.x<= high
            && low<=psys[i].pos.y && psys[i].pos.y<= high){
            psys[i].mass = rho_in / n_glb;
            psys[i].eng = pres / ((gamma - 1.0) * rho_in);
        }
        else{
            psys[i].mass = rho_out / n_glb;
            psys[i].eng = pres / ((gamma - 1.0) * rho_out);
        }
    }
}


template<class Tpsys>
void SetParticlesContact3D(Tpsys & psys,
                           const PS::S32 n_glb,                           
                           PS::S32 & n_loc,
                           PS::F64 & time_sys){
    const PS::S32 n_one_dim = cbrt( (double)n_glb * 1.00001);
    assert(n_one_dim*n_one_dim*n_one_dim == n_glb);
    const PS::S32 n_two_dim = n_one_dim * n_one_dim;
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::F64 dl = 1.0 / n_one_dim;
    const PS::F64 pres = 2.5;
    const PS::F64 gamma = 5.0 / 3.0;
    const PS::F64 rho_in = 4.0;
    const PS::F64 rho_out = 1.0;
    const PS::F64 low = 0.25;
    const PS::F64 high = 0.75;
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    time_sys = 0.0;
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        PS::S32 ix = (i_h+i) / n_two_dim;
        PS::S32 iy = ((i_h+i) % n_two_dim) / n_one_dim;
        PS::S32 iz = ((i_h+i) % n_two_dim) % n_one_dim;
        psys[i].id = i_h + i;
        psys[i].mass = 1.0 / n_glb;
        psys[i].vel = 0.0;
        psys[i].r_search = dl * 2.0;
        psys[i].pos.x = ix * dl;
        psys[i].pos.y = iy * dl;
        psys[i].pos.z = iz * dl;
        if( low<=psys[i].pos.x && psys[i].pos.x<= high
            && low<=psys[i].pos.y && psys[i].pos.y<= high 
            && low<=psys[i].pos.z && psys[i].pos.z<= high ){
            psys[i].mass = rho_in / n_glb;
            psys[i].eng = pres / ((gamma - 1.0) * rho_in);
        }
        else{
            psys[i].mass = rho_out / n_glb;
            psys[i].eng = pres / ((gamma - 1.0) * rho_out);
        }
    }
}

template<class Tpsys>
void SetParticlesUniform3D(Tpsys & psys,
                           const PS::S32 n_glb,                           
                           PS::S32 & n_loc,
                           PS::F64 & time_sys){
    const PS::S32 n_one_dim = cbrt( (double)n_glb * 1.00001);
    assert(n_one_dim*n_one_dim*n_one_dim == n_glb);
    const PS::S32 n_two_dim = n_one_dim * n_one_dim;
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::F64 rho = 1.0;
    const PS::F64 L = 1.0;
    const PS::F64 M = L * L * L * rho;
    const PS::F64 dl = L / n_one_dim;
    const PS::F64 pres = 2.5;
    const PS::F64 gamma = 5.0 / 3.0;

    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    time_sys = 0.0;
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        PS::S32 ix = (i_h+i) / n_two_dim;
        PS::S32 iy = ((i_h+i) % n_two_dim) / n_one_dim;
        PS::S32 iz = ((i_h+i) % n_two_dim) % n_one_dim;
        psys[i].id = i_h + i;
        psys[i].mass = M / n_glb;
        psys[i].vel = 0.0;
        psys[i].r_search = dl * 2.0;
        psys[i].pos.x = ix * dl;
        psys[i].pos.y = iy * dl;
        psys[i].pos.z = iz * dl;
        psys[i].mass = rho / n_glb;
        psys[i].eng = pres / ((gamma - 1.0) * rho);
    }
}


template<class Tpsys>
void SetParticlesUniformRandom3D(Tpsys & psys,
				 const PS::S32 n_glb,                           
				 PS::S32 & n_loc,
				 PS::F64 & time_sys){
    const PS::S32 n_one_dim = cbrt( (double)n_glb * 1.00001);
    assert(n_one_dim*n_one_dim*n_one_dim == n_glb);
    //const PS::S32 n_two_dim = n_one_dim * n_one_dim;
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::F64 rho = 1.0;
    const PS::F64 L = 1.0;
    const PS::F64 M = L * L * L * rho;
    const PS::F64 dl = L / n_one_dim;
    const PS::F64 pres = 2.5;
    const PS::F64 gamma = 5.0 / 3.0;

    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    time_sys = 0.0;
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    PS::MTTS mt;
    mt.init_genrand(PS::Comm::getRank());
    for(size_t i=0; i<n_loc; i++){
        psys[i].id = i_h + i;
        psys[i].mass = M / n_glb;
        psys[i].vel = 0.0;
        psys[i].r_search = dl * 2.0;
        psys[i].pos.x = L * mt.genrand_res53();
        psys[i].pos.y = L * mt.genrand_res53();
        psys[i].pos.z = L * mt.genrand_res53();
        psys[i].mass = rho / n_glb;
        psys[i].eng = pres / ((gamma - 1.0) * rho);
    }
}



int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);
    PS::F64 theta = 0.5;
    const PS::S32 n_leaf_limit = 4;
    PS::S32 n_group_limit = 8;
    char dir_name[1024];
    PS::S32 n_glb = 100;
    int c;
    while((c=getopt(argc,argv,"o:n:N:h")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            std::cerr<<"dir_name: "<<dir_name<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'N':
            n_glb = atoi(optarg);
            std::cerr<<"n_glb="<<n_glb<<std::endl;
            break;
        case 'h':
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 8.0)"<<std::endl;
            std::cerr<<"N: number of particles"<<std::endl;
            return 0;
        }
    }
    std::ofstream fout_tcal;
#if 0
    if(PS::Comm::getRank() == 0){
#else
    if(1){
#endif
        char sout_tcal[1024];
        sprintf(sout_tcal, "%s/t-tcal_%05d_%05d.dat", 
		dir_name, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
        std::cerr<<sout_tcal<<std::endl;
        fout_tcal.open(sout_tcal);
    }

    const PS::F64 gamma = 1.4;
#ifdef ONE_DIM_SIM
    const PS::F64 time_end = 0.1;
//#elif  defined(THREE_DIM_SIM)
#else
    const PS::F64 time_end = 10.0;
#endif

    PS::ParticleSystem<FullPtcl> system_sph;
    system_sph.initialize();
    PS::S32 n_loc = 0;
    PS::F64 time_sys = 0.0;
#ifdef ONE_DIM_SIM
    //SetParticlesShockTube(system_sph, n_glb, n_loc, time_sys);
    SetParticlesStrongShock(system_sph, n_glb, n_loc, time_sys);
    for(int i=0; i<10; i++){
        std::cerr<<"system_sph[i].pos="<<system_sph[i].pos<<std::endl;
    }
#elif  defined(TWO_DIM_SIM)
    SetParticlesContact2D(system_sph, n_glb, n_loc, time_sys);
#elif  defined(THREE_DIM_SIM)
    //SetParticlesContact(system_sph, n_glb, n_loc, time_sys);
    //SetParticlesUniform3D(system_sph, n_glb, n_loc, time_sys);
    SetParticlesUniformRandom3D(system_sph, n_glb, n_loc, time_sys);
#endif
    PS::DomainInfo dinfo;
    dinfo.initialize();

#ifdef ONE_DIM_SIM
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
#elif  defined(TWO_DIM_SIM)
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
    dinfo.setPosRootDomain(0.0, 1.0);
#elif  defined(THREE_DIM_SIM)
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(0.0, 1.0);
#endif

    //dinfo.collectSampleParticle(system_sph, system_sph.getNumberOfParticleLocal());
    dinfo.collectSampleParticle(system_sph);
    dinfo.decomposeDomain();
    system_sph.exchangeParticle(dinfo);
    
    n_loc = system_sph.getNumberOfParticleLocal();

    PS::TreeForForceShort<ResultDens, EPIDens, EPJDens>::Symmetry tree_dens;
    //PS::TreeForForceShort<ResultDens, EPIDens, EPJDens>::Gather tree_dens;
    tree_dens.initialize(n_glb, theta, n_leaf_limit, n_group_limit);

    for(bool repeat = true; repeat;){
        bool repeat_loc = false;
        bool clear_force = false;
        tree_dens.calcForceAllAndWriteBack(CalcDensEpEp(), system_sph, dinfo, clear_force);
        //tree_dens.calcForceAllAndWriteBackWithCheck(CalcDensEpEp(), system_sph, dinfo, clear_force);
        for(PS::S32 i=0; i<n_loc; i++){
            if(system_sph[i].r_search > 0.0){
                repeat_loc = true;
            }
        }
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
    }
    CalcPressureAndSoundVelocity(system_sph, gamma);
    //CalcEnergyAndSoundVelocity(system_sph, gamma);

    PS::TreeForForceShort<ResultForce, EPIForce, EPJForce>::Symmetry tree_force;
    tree_force.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    tree_force.calcForceAllAndWriteBack(CalcForceEpEp(), system_sph, dinfo);
    //tree_force.calcForceAllAndWriteBackWithCheck(CalcForceEpEp(), system_sph, dinfo);

    PS::F64 Tloop = 0.0;
    PS::S32 snp_id = 0;
    PS::S64 n_loop = 0;
    PS::Timer timer;
    timer.initialize(fout_tcal);
    while(time_sys < time_end){

#ifdef DEBUG_PRINT
        PS::S32 check_id = 0;
        std::cerr<<"check: "<<check_id<<std::endl;
        std::cerr<<"n_loop="<<n_loop<<std::endl;
        check_id++;
#endif

        if(n_loop % 30 == 0){
            timer.initialize(fout_tcal);
        }
        timer.reset();
        timer.start();

/*
        if(n_loop % 100 == 0){
#if 0
            WriteSPHFormat(system_sph, time_sys, snp_id++, dir_name);
#else
            FileHeader header;
            header.time = time_sys;
            header.n_body = system_sph.getNumberOfParticleLocal();
			char filename[256];
			sprintf(filename, "%s/snap%05d.dat", dir_name, snp_id++);
			system_sph.writeParticleAscii(filename, "%s_%05d_%05d.dat", header);
#endif
        }
*/

#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("WriteFile");
        PS::F64 dt = CalcDt(system_sph);

#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("CalcDt");
        time_sys += dt;
        Kick1(system_sph, dt*0.5);

#ifdef DEBUG
        std::vector< std::pair<PS::S64, PS::F64> > pos_id_old;
        pos_id_old.reserve(1000);
        for(int i=0; i<n_loc; i++)  pos_id_old.push_back( std::make_pair(system_sph[i].id, system_sph[i].pos.x) );
        std::sort(pos_id_old.begin(), pos_id_old.end());
#endif
        Drift(system_sph, dt);
#ifdef DEBUG
        std::vector< std::pair<PS::S64, PS::F64> > pos_id_new;
        pos_id_new.reserve(1000);
        for(int i=0; i<n_loc; i++)  pos_id_new.push_back( std::make_pair(system_sph[i].id, system_sph[i].pos.x) );
        std::sort(pos_id_new.begin(), pos_id_new.end());
        for(int i=0; i<pos_id_new.size()-1; i++){
            if( pos_id_new[i].second > pos_id_new[i+1].second ){
                std::cerr<<"i="<<i<<std::endl;
                std::cerr<<"pos_id_new[i].second="<<pos_id_new[i].second<<std::endl;
                std::cerr<<"pos_id_new[i+1].second="<<pos_id_new[i+1].second<<std::endl;
                std::cerr<<"pos_id_old[i].second="<<pos_id_old[i].second<<std::endl;
                std::cerr<<"pos_id_old[i+1].second="<<pos_id_old[i+1].second<<std::endl;
            }
        }
#endif

#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("kick+drift");

        system_sph.adjustPositionIntoRootDomain(dinfo);

#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("adjustPos");

        dinfo.collectSampleParticle(system_sph);

#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("collect");

        dinfo.decomposeDomain();

#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("decompose");

        system_sph.exchangeParticle(dinfo);
	
#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("exchangeParticle");

        n_loc = system_sph.getNumberOfParticleLocal();
#pragma omp parallel for
        for(PS::S32 i=0; i<n_loc; i++){
            system_sph[i].r_search = system_sph[i].r_kernel * 1.2;
        }

        //tree_dens.calcForceAllAndWriteBack(CalcDensEpEp(), system_sph, dinfo);
        Tloop = PS::GetWtime();
        PS::S32 n_repeat = 0;
        for(bool repeat = true, clear_force = true; repeat;){
            n_repeat++;
            bool repeat_loc = false;
            //tree_dens.calcForceAllAndWriteBack(CalcDensEpEp(), system_sph, dinfo, clear_force);
            //tree_dens.calcForceAllAndWriteBackWithCheck(CalcDensEpEp(), system_sph, dinfo, clear_force);
            tree_dens.calcForceAllAndWriteBackWithTimer(CalcDensEpEp(), system_sph, dinfo, timer, clear_force);
            clear_force = false;
            for(PS::S32 i=0; i<n_loc; i++){
                if(system_sph[i].r_search > 0.0){
                    repeat_loc = true;
                }
            }
            repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
            timer.restart("checkRepeatCondition");
        }	



        CalcPressureAndSoundVelocity(system_sph, gamma);
	
#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("calcPAndC");

        //tree_force.calcForceAllAndWriteBack(CalcForceEpEp(), system_sph, dinfo);
        //tree_force.calcForceAllAndWriteBackWithCheck(CalcForceEpEp(), system_sph, dinfo);
        tree_force.calcForceAllAndWriteBackWithTimer(CalcForceEpEp(), system_sph, dinfo, timer);

        Tloop = PS::GetWtime() - Tloop;

        Kick2(system_sph, dt*0.5);

#ifdef DEBUG_PRINT
        std::cerr<<"check: "<<check_id<<std::endl;
        check_id++;
#endif
        timer.restart("kick2");

        std::cout<<"time_sys= "<<time_sys<<" n_loop= "<<n_loop<<std::endl;
        fout_tcal<<"time_sys= "<<time_sys<<" n_loop= "<<n_loop<<std::endl;
        fout_tcal<<"system_sph.getNumberOfParticleLocal()= "<<system_sph.getNumberOfParticleLocal()<<std::endl;
        fout_tcal<<"tree_dens.getMemSizeUsed()= "<<tree_dens.getMemSizeUsed()*1e-9<<" [Gbyte]";
        fout_tcal<<" tree_force.getMemSizeUsed()= "<<tree_force.getMemSizeUsed()*1e-9<<" [Gbyte]";
        fout_tcal<<" system_sph.getMemSizeUsed()= "<<system_sph.getMemSizeUsed()*1e-9<<" [Gbyte]"<<std::endl;
        fout_tcal<<"n_repeat= "<<n_repeat<<std::endl;
        tree_dens.dump_calc_cost(PS::Comm::getMaxValue(Tloop), fout_tcal);
        tree_force.dump_calc_cost(PS::Comm::getMaxValue(Tloop), fout_tcal);
        timer.dump(fout_tcal);
        fout_tcal<<std::endl;
        n_loop++;
    }

    PS::Finalize();

    return 0;
}
