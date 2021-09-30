#include"force_sunway_impl.hpp"
#include"kepler.hpp"

// f(x) = x/sigma^2*exp(-x^2/(2sigma^2))
// F(x) = 1.0 - exp(-x^2/(2sigma^2))
// x = sqrt(-2*sigma^2*ln(1-F))
// 2sigma^2 = <e^2> or <i^2>
// <> means R.M.S.
double RayleighDistribution(const double sigma){
    static PS::MTTS mt;
    static bool first = true;
    if(first){
        mt.init_genrand( PS::Comm::getRank() );
        first = false;
    }
    double F = mt.genrand_res53();
    return sqrt( -2.0*sigma*sigma*log(1.0-F));
}

// p(a)da = C*a^(p+1)*da
// P(a) = 2*C*(a^0.5 - a_in^0.5)
double HayashiDistributionWithIceLine(const double a_in,
                                      const double a_out,
                                      const double a_ice, 
                                      const double f_ice=1.0,
                                      const double p_sigma=-1.5){
    const PS::F64 p_mass = 2.0 + p_sigma;
    static PS::MTTS mt;
    static bool first = true;
    if(first){
        mt.init_genrand( PS::Comm::getRank() );
        first = false;
    }
    double P = mt.genrand_res53();
    double ret = 0.0;
    double C = 0.0;
    double P_ice = 0.0;
    if(a_ice <= a_in || a_ice >= a_out){
        C = p_mass / (pow(a_out, p_mass) - pow(a_in, p_mass));
        ret = pow(P*p_mass/C + pow(a_in, p_mass), 1.0/p_mass);
    }
    else{
        C = p_mass / ( f_ice*(pow(a_out, p_mass)-pow(a_ice, p_mass)) + pow(a_ice, p_mass) - pow(a_in, p_mass) );
        P_ice = C/p_mass*(pow(a_ice, p_mass) - pow(a_in, p_mass));
        if(P < P_ice){
            ret = pow(P*p_mass/C + pow(a_in, p_mass), 1.0/p_mass);
        }
        else{
            ret = pow( ((P*p_mass/C + pow(a_in, p_mass) - pow(a_ice, p_mass))/f_ice + pow(a_ice, p_mass)), 1.0/p_mass); // beyond the ice line
        }
    }
    return ret;
}

void MakeKeplerDisk(PS::F64 & mass_planet_glb,
                    PS::F64 *& mass,
                    PS::F64vec *& pos,
                    PS::F64vec *& vel,
                    const long long int n_glb,
                    const long long int n_loc,
                    const double a_in, // [AU]
                    const double a_out, // [AU]
                    const double e_rms, // normalized
                    const double i_rms, // normalized
                    const double dens = 10.0, // [g/cm^2]
                    const double mass_sun = 1.0, //[m_sun]
                    const double a_ice = 0.0,
                    const double f_ice = 1.0,
                    const double power = -1.5,
                    const int seed = 0){
    static const double PI = atan(1.0) * 4.0;
    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank()+seed*PS::Comm::getNumberOfProc() );
    mass_planet_glb = 2.0 * PI * dens / (2.0+power) * ( pow(a_out, 2.0+power) - pow(a_in, 2.0+power) );
    const double m_planet = mass_planet_glb / n_glb;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];
    const double h = pow(2.0*m_planet / (3.0*mass_sun), 1.0/3.0);
    PS::F64 e_ave = 0.0;
    PS::F64 i_ave = 0.0;    
    const PS::F64 e_sigma = sqrt(0.5*e_rms*e_rms); // this is right procedure
    const PS::F64 i_sigma = sqrt(0.5*i_rms*i_rms);
    for(long long int i=0; i<n_loc; i++){
        mass[i] = m_planet;
        double ax = HayashiDistributionWithIceLine(a_in, a_out, a_ice, f_ice, power);
        double ecc_tmp;
        do{
            ecc_tmp = RayleighDistribution(e_sigma) * h;
        }while(ecc_tmp >= 1.0);
        double ecc = ecc_tmp;
        double inc_tmp;
        do{
            inc_tmp = RayleighDistribution(i_sigma) * h;
        }while(inc_tmp >= 0.5*PI || inc_tmp <= -0.5*PI);
        double inc = inc_tmp;
        PS::F64vec pos_dummy, vel_dummy;
        double omg = 2.0 * PI * mt.genrand_res53();
        double OMG = 2.0 * PI * mt.genrand_res53();
        double l = 2.0 * PI * mt.genrand_res53(); // mean anomayl
        double u = solve_keplereq(l, ecc); // eccentric anomayl
        OrbParam2PosVel(pos_dummy, pos[i],   vel_dummy, vel[i],
                        mass_sun,  m_planet, ax,        ecc,
                        inc,       OMG,      omg,       u);
        //std::cout<<"pos[i]= "<<pos[i]<<" ax= "<<ax<<" ecc= "<<ecc<<" inc= "<<inc<<std::endl;
        e_ave += ecc*ecc;
        i_ave += inc*inc;        
    }
    PS::F64 e_ave_glb = sqrt(PS::Comm::getSum(e_ave)/n_glb);
    PS::F64 i_ave_glb = sqrt(PS::Comm::getSum(i_ave)/n_glb);    
    if(PS::Comm::getRank() == 0){
        std::cerr<<"e_ave_glb="<<e_ave_glb<<std::endl;
        std::cerr<<"e_ave_hill="<<e_ave_glb/h<<std::endl;
        std::cerr<<"e_rms="<<e_rms<<std::endl;
        std::cerr<<"i_ave_glb="<<i_ave_glb<<std::endl;
        std::cerr<<"i_ave_hill="<<i_ave_glb/h<<std::endl;
        std::cerr<<"i_rms="<<i_rms<<std::endl;        
    }
}


template<class Tpsys>
void SetParticleKeplerDisk(Tpsys & psys,
                           const PS::S64 n_glb,
                           const PS::F64 ax_in,
                           const PS::F64 ax_out,
                           const PS::F64 ecc_rms, // normalized
                           const PS::F64 inc_rms, // normalized
                           const PS::F64 dens = 10.0, // disk density
                           const PS::F64 mass_sun = 1.0,
                           const double a_ice = 0.0,
                           const double f_ice = 1.0,
                           const double power = -1.5,
                           const PS::S32 seed = 0){
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    //auto n_loc = n_glb;
    PS::S64 n_loc = n_glb;
#if 1 // for debug
    if(my_rank == 0){
        n_loc = n_glb;
    }
    else{
        n_loc = 0;
    }
    psys.setNumberOfParticleLocal(n_loc);
#else
    // original
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
#endif
    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    PS::F64 mass_planet_glb;
    MakeKeplerDisk(mass_planet_glb, mass, pos, vel, n_glb, n_loc,
                   ax_in, ax_out, ecc_rms, inc_rms, dens, mass_sun, a_ice, f_ice, power, seed);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(PS::S32 i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}


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
    for(size_t i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(size_t i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(size_t i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }

    double r_max_sq = -1.0;
    for(int i=0; i<n_loc; i++){
        if(r_max_sq < pos[i] * pos[i]){
            r_max_sq = pos[i] * pos[i];
        }
    }
    //std::cout<<"r_max= "<<sqrt(r_max_sq)<<std::endl;
}


template<class Tpsys>
void SetParticlePlummer(Tpsys & psys,
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
    for(size_t i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

template<class Tpsys>
void Kick(Tpsys & system,
          const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].vel  += system[i].acc * dt;
    }
}
