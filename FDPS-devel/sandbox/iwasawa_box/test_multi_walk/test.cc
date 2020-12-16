#include<iostream>
#include<unistd.h>
#include<particle_simulator.hpp>
#include"class.hpp"
#include"force.hpp"

struct ClacForceDirect{
    void operator () (const Epi     * epi,
                      const PS::S32  n_epi,
                      const Epj     * epj,
                      const PS::S32  n_epj,
                      Force * force){
        PS::F32 eps2 = Epi::eps * Epi::eps;
        for(int ip=0; ip<n_epi; ip++){
            force[ip].acc = 0.0;
            force[ip].pot = 0.0;
            for(int jp=0; jp<n_epj; jp++){
                if( epi[ip].id == epj[jp].id) continue;
                const PS::F64 mj = epj[jp].mass;
                PS::F64vec rij = epi[ip].pos - epj[jp].pos;
                PS::F64 r2 = rij*rij + eps2;
                PS::F64 inv_r = 1.0 / sqrt(r2);
                PS::F64 inv_r3 = inv_r * inv_r * inv_r;
                force[ip].acc -= mj * inv_r3 *rij;
                force[ip].pot -= mj * inv_r;
            }
        }
    }
};


void MakePlummerModel(const double mass_glb,
                      const long long int n_glb,
                      const long long int n_loc,
                      double *& mass,
                      PS::F64vec *& pos,
                      PS::F64vec *& vel,
                      const double eng = -0.25,
                      const int seed = 0){

    assert(eng < 0.0);
    static const double PI = atan(1.0) * 4.0;
    const double r_cutoff = 22.8 / (-3.0 * PI * mass_glb * mass_glb / (64.0 * -0.25)); // 22.8 is cutoff in Nbody units
    //const double r_cutoff = 22.8 * 0.25;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];

    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        double r_tmp = 9999.9;
        while(r_tmp > r_cutoff){ 
            double m_tmp = mt.genrand_res53();
            r_tmp = 1.0 / sqrt( pow(m_tmp, (-2.0/3.0)) - 1.0);
        }
        double phi = 2.0 * PI * mt.genrand_res53();
        double cth = 2.0 * (mt.genrand_real2() - 0.5);
        double sth = sqrt(1.0 - cth*cth);
        pos[i][0] = r_tmp * sth * cos(phi);
        pos[i][1] = r_tmp * sth * sin(phi);
        pos[i][2] = r_tmp * cth;
        while(1){
            const double v_max = 0.1;
            const double v_try = mt.genrand_res53();
            const double v_crit = v_max * mt.genrand_res53();
            if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                const double ve = sqrt(2.0) * pow( (r_tmp*r_tmp + 1.0), -0.25 );
                phi = 2.0 * PI * mt.genrand_res53();
                cth = 2.0 * (mt.genrand_res53() - 0.5);
                sth = sqrt(1.0 - cth*cth);
                vel[i][0] = ve * v_try * sth * cos(phi);
                vel[i][1] = ve * v_try * sth * sin(phi);
                vel[i][2] = ve * v_try * cth;
                break;
            }
        }
    }

    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    double  cm_mass = 0.0;
    for(long long int i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(long long int i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(long long int i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }

    double r_max_sq = -1.0;
    for(int i=0; i<n_loc; i++){
        if(r_max_sq < pos[i] * pos[i]){
            r_max_sq = pos[i] * pos[i];
        }
    }
    std::cout<<"r_max= "<<sqrt(r_max_sq)<<std::endl;
}


template<class Tpsys>
void SetParticlesPlummer(Tpsys & psys,
                         const PS::S64 n_glb,
                         PS::S32 & n_loc,  
                         PS::F32 & t_sys){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(long long int i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

PS::F64 Epi::eps = 1.0/1024.0;
int main(int argc, char *argv[]){

    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    
    PS::Initialize(argc, argv);
    std::cout<<"PS::Comm::getNumberOfThread="<<PS::Comm::getNumberOfThread()<<std::endl;
    std::cout<<"omp_get_num_threads()="<<omp_get_num_threads()<<std::endl;
    std::cout<<"omp_get_max_threads()="<<omp_get_max_threads()<<std::endl;

    PS::F32 theta = 0.0;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 n_smp_ave = 30;
    PS::F32 time_end = 10.0;
    char dir_name[1024];
    long long int n_tot = 16384;
    int c;
    while((c=getopt(argc,argv,"o:t:T:n:N:s:l:h")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 't':
            theta = atof(optarg);
            std::cerr<<"theta="<<theta<<std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr<<"time_end="<<time_end<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'N':
            n_tot = atol(optarg);
            std::cerr<<"n_tot="<<n_tot<<std::endl;
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'h':
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 10.0)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"N: n_tot (dafult: 16384)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 30)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
	    PS::Finalize();
            return 0;
        }
    }

    PS::ParticleSystem<FP> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    PS::S32 n_grav_loc;
    PS::F32 time_sys;
    PS::S64 n_grav_glb = n_tot;
    SetParticlesPlummer(system_grav, n_grav_glb, n_grav_loc, time_sys);
    
    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);

    PS::S64 n_walk_limit = 13;
    PS::TreeForForceLong<Force, Epi, Epj>::Monopole tree_grav;
    tree_grav.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP, RetrieveKernel,
						system_grav, dinfo, n_walk_limit, true);


    Force * fd = new Force[system_grav.getNumberOfParticleLocal()];
    tree_grav.calcForceDirect(ClacForceDirect(), fd, dinfo);
    for(int i=0; i<10; i++){
	std::cerr<<"i="<<i<<" id="<<system_grav[i].id<<" pos="<<system_grav[i].pos<<std::endl;
        std::cerr<<"acc=      "<<system_grav[i].acc<<" pot=      "<<system_grav[i].pot<<std::endl;
        std::cerr<<"fd[i].acc="<<fd[i].acc         <<" fd[i].pot="<<fd[i].pot<<std::endl;
    }
    double Epot_tr_loc = 0.0;
    double Epot_di_loc = 0.0;
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        Epot_tr_loc += system_grav[i].mass * system_grav[i].pot;
        Epot_di_loc += system_grav[i].mass * fd[i].pot;
    }
    double Epot_tr_glb = 0.5 * PS::Comm::getSum(Epot_tr_loc);
    double Epot_di_glb = 0.5 * PS::Comm::getSum(Epot_di_loc);
    std::cerr<<"Epot_tr_glb="<<Epot_tr_glb<<" Epot_di_glb="<<Epot_di_glb<<std::endl;

    PS::Finalize();
    
    return 0;
}

