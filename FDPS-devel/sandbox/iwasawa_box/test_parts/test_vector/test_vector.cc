#include<iostream>
#include<fstream>
#include<unistd.h>
#include<random>
#include<particle_simulator.hpp>

class FPGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    template<class Tforce>
    void copyFromForce(const Tforce & force){
        acc = force.acc;
        pot = force.pot;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z);
    }
    void readAscii(FILE* fp){
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z);
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
    for(int i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(int i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(int i=0; i<n_loc; i++){
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
    for(int i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

int main(int argc, char *argv[]){
    
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

    PS::F64vec x;
    PS::F64vec y;
    x[0] = 0.0;
    x[1] = 1.0;
    x[2] = 2.0;
    std::cerr<<"x[-1]="<<x[-1]<<std::endl;    
    x[100] = 3.0;
    std::cerr<<"x[0]="<<x[0]<<std::endl;
    std::cerr<<"y[0]="<<y[0]<<std::endl;
    /*
    long long int n_tot = 16384;
    int c;
    while((c=getopt(argc,argv,"N:h")) != -1){
        switch(c){
        case 'N':
            n_tot = atol(optarg);
            std::cerr<<"n_tot="<<n_tot<<std::endl;
            break;
        case 'h':
            std::cerr<<"N: n_tot (dafult: 16384)"<<std::endl;
	    PS::Finalize();
            return 0;
        }
    }

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_grav_loc;
    PS::F32 time_sys;
    PS::S64 n_grav_glb = n_tot;
    SetParticlesPlummer(system_grav, n_grav_glb, n_grav_loc, time_sys);

    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.collectSampleParticle(system_grav, true);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_grav_loc = system_grav.getNumberOfParticleLocal();

    FPGrav fp_add[4];
    fp_add[0].id = n_grav_loc;
    fp_add[1].id = n_grav_loc + 1;
    fp_add[2].id = n_grav_loc + 2;
    fp_add[3].id = n_grav_loc + 3;
    system_grav.addOneParticle(fp_add[0]);
    system_grav.addOneParticle(fp_add[1]);
    system_grav.addOneParticle(fp_add[2]);
    system_grav.addOneParticle(fp_add[3]);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
    std::cerr<<"after adding"<<std::endl;
    std::cerr<<"n_grav_loc="<<n_grav_loc<<std::endl;
    const PS::S32 n_remove = n_grav_loc / 2;
    
    PS::S32 * idx_remove = new PS::S32[n_remove];
    PS::S32 seed = 0;
    std::mt19937 mt(seed);
    std::uniform_int_distribution<PS::S32> rand(0, n_grav_loc-1);
    std::cerr<<"remove index"<<std::endl;
    std::cerr<<"n_remove="<<n_remove<<std::endl;
    for(PS::S32 i=0; i<n_remove; i++){
	idx_remove[i] = rand(mt);
	std::cerr<<"idx_remove[i]="<<idx_remove[i]<<std::endl;
    }
    system_grav.removeParticle(idx_remove, n_remove);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
    std::cerr<<"after removing"<<std::endl;
    std::cerr<<"n_grav_loc="<<n_grav_loc<<std::endl;
    
    ///////////
    // CHECK
    for(PS::S32 i=0; i<n_grav_loc; i++){
	for(PS::S32 j=0; j<n_remove; j++){
	    assert(system_grav[i].id != idx_remove[j]);
	}
    }
    std::sort(idx_remove, idx_remove+n_remove);
    PS::S32 n_remove_unique = std::unique(idx_remove, idx_remove+n_remove) - idx_remove;
    std::cerr<<"n_remove="<<n_remove<<std::endl;
    std::cerr<<"n_remove_unique="<<n_remove_unique<<std::endl;
    for(PS::S32 i=0; i<n_remove_unique; i++){
	//std::cerr<<"i= "<<i<<" idx_remove[i]= "<<idx_remove[i]<<std::endl;
    }
    PS::S32 * idx_check = new PS::S32[n_grav_loc+n_remove_unique];
    for(PS::S32 i=0; i<n_grav_loc; i++){
	idx_check[i] = system_grav[i].id;
    }
    PS::S32 n_cnt = 0;
    for(PS::S32 i=n_grav_loc; i<n_grav_loc+n_remove_unique; i++){
	idx_check[i] = idx_remove[n_cnt];
	n_cnt++;
    }
    std::sort(idx_check, idx_check+(n_grav_loc+n_remove_unique));
    bool flag_err = false;
    for(PS::S32 i=0; i<n_grav_loc+n_remove_unique; i++){
	//std::cerr<<"i="<<i<<" idx_check[i]="<<idx_check[i]<<std::endl;
	//assert( i==idx_check[i] );
	//assert( (idx_check[i]-idx_check[i-1]) >=0 );
	if( i != idx_check[i] ){
	    std::cerr<<"i="<<i<<" idx_check[i]="<<idx_check[i]<<std::endl;
	    std::cerr<<"FAIL"<<std::endl;
	    flag_err = true;
	    break;
	}
    }
    if(!flag_err) std::cerr<<"PASS"<<std::endl;
    */
    
    PS::Finalize();
    return 0;
}
