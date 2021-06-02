#pragma once

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
    /*
    double ret = 0.0;
    do{
	ret = sqrt( -2.0*sigma*sigma*log(1.0-F));
    }while(ret >= 1.0);
    return ret;    
    */
    return sqrt( -2.0*sigma*sigma*log(1.0-F));
}

// p(a)da = C*a^-0.5*da
// P(a) = 2*C*(a^0.5 - a_in^0.5)
// a = (P/(2C) + a_in^0.5)^2
double HayashiDistribution(const double a_in, const double a_out){
    static PS::MTTS mt;
    static bool first = true;
    if(first){
	mt.init_genrand( PS::Comm::getRank() );
	first = false;
    }
    const double C = 0.5 / (sqrt(a_out) - sqrt(a_in));
    double P = mt.genrand_res53();
    double ret = P/(2*C) + sqrt(a_in);
    ret *= ret;
    return ret;
}


double AU2CM(const double x){
    return x * 1.49597871e13;
}

double CM2AU(const double x){
    return x / 1.49597871e13;
}

// [M] = 1 [g]
// [L] = 1 [cm]
// [T] = 3.871e3[sec]
// [V] = 2.5833118e-4 [cm/sec] = 2.5833118e-9 [km/sec]
// 1yr = 8.146732e3[T]
// or
// G = 1 [L]^3/([M][T]^2)
// [L] = 1[AU] = 1.495978707e8[km]
// [M] = 1[Msun] = 1.989e30[kg]
// [T] = 5.02198050479e6 [sec] = 0.15924595715 [yr]
// [V] = 29.7886203575 [km/sec]
void MakeKeplerDisk(PS::F64 & mass_planet_glb,
                    PS::F64 *& mass,
                    PS::F64vec *& pos,
                    PS::F64vec *& vel,
                    const long long int n_glb,
                    const long long int n_loc,
                    const double a_in, // [AU]
                    const double a_out, // [AU]
                    const double e_sigma, // normalized
                    const double i_sigma, // normalized
                    const double dens=10.0, // [g/cm^2]
                    const double mass_sun = 1.0, //[m_sun]
		    const int seed = 0
    ){
    static const double mass_sun_gram = 1.989e33; //[g]
    PS::MTTS mt;
    //mt.init_genrand( PS::Comm::getRank() );
    mt.init_genrand( PS::Comm::getRank()+seed*PS::Comm::getNumberOfProc() );
    static const double PI = atan(1.0) * 4.0;
    const double AU = AU2CM(1.0);
    mass_planet_glb = 4.0 * PI * dens * ( sqrt(a_out) - sqrt(a_in) ) * AU * AU / mass_sun_gram; // [Msun]
    const double m_planet = mass_planet_glb / n_glb;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];
    const double h = pow(2.0*m_planet / (3.0*mass_sun), 1.0/3.0);
    //std::cerr<<"h="<<h<<std::endl;
    //std::cerr<<"e_sigma="<<e_sigma<<std::endl;
    //std::cerr<<"i_sigma="<<i_sigma<<std::endl;
    PS::F64 e_ave = 0.0;
    for(long long int i=0; i<n_loc; i++){
        mass[i] = m_planet;
        double ax = HayashiDistribution(a_in, a_out);
        double ecc = RayleighDistribution(e_sigma) * h;
        double inc = RayleighDistribution(i_sigma) * h;
        //double ecc = RayleighDistribution( sqrt(e_sigma*e_sigma/(2.0-0.5*PI)) ) * h;
        //double inc = RayleighDistribution( sqrt(i_sigma*i_sigma/(2.0-0.5*PI)) ) * h;
        //std::cerr<<"mass[i]="<<mass[i]<<" h="<<h<<" inc="<<inc<<" ecc="<<ecc<<std::endl;
        PS::F64vec pos_dummy, vel_dummy;
        double omg = 2.0 * PI * mt.genrand_res53();
        double OMG = 2.0 * PI * mt.genrand_res53();
        double l = 2.0 * PI * mt.genrand_res53(); // mean anomayl
        double u = solve_keplereq(l, ecc); // eccentric anomayl
        OrbParam2PosVel(pos_dummy, pos[i],   vel_dummy, vel[i],
                        mass_sun,  m_planet, ax,        ecc,
                        inc,       OMG,      omg,       u);
        e_ave += ecc;
    }
    //std::cerr<<"e_ave="<<e_ave/n_loc<<std::endl;
    //std::cerr<<"e_ave_hill="<<e_ave/n_loc/h<<std::endl;
}

template<class Tpsys>
void SetParticleKeplerDisk(Tpsys & psys,
                           const PS::S64 n_glb,
                           PS::S32 & n_loc,
                           PS::F64 & t_sys,
                           const PS::F64 ax_in, // [AU]
                           const PS::F64 ax_out, // [AU]
                           const PS::F64 ecc_sigma, // normalized
                           const PS::F64 inc_sigma, // normalized
                           const PS::F64 dens = 10.0, // [g/cm^2]
                           const PS::F64 mass_sun = 1.0, //[m_sun]
			   const PS::S32 seed = 0
    ){ 
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
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
    t_sys = 0.0;

    PS::F64 mass_planet_glb;
    MakeKeplerDisk(mass_planet_glb, mass, pos, vel, n_glb, n_loc,
                   ax_in, ax_out, ecc_sigma, inc_sigma, dens, mass_sun, seed);
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
    for(PS::S32 i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(PS::S32 i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(PS::S32 i=0; i<n_loc; i++){
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
void SetParticlePlummer(Tpsys & psys,
                        const PS::S64 n_glb,
                        PS::S32 & n_loc,  
                        PS::F64 & t_sys){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
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
    t_sys = 0.0;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
#ifdef BINARY
    for(PS::S32 i=0; i<n_loc; i+=2){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
	psys[i+1].mass = mass[i];
        psys[i+1].pos = pos[i] + PS::F64vec(0.02, 0.02, 0.02);
        psys[i+1].vel = vel[i];
        psys[i+1].id = i_h + i+1;
    }
#else
    for(PS::S32 i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
#endif
    delete [] mass;
    delete [] pos;
    delete [] vel;

}

template<class Tpsys>
void WriteFile(const Tpsys & psys,
               const char * file_name,
               const PS::F64 time_sys){
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S64 n_glb = psys.getNumberOfParticleGlobal();
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    std::ofstream fout;
    if(my_rank == 0){
        fout.open(file_name);
        fout<<std::setprecision(15);
        fout<<n_glb<<std::endl;
        fout<<time_sys<<std::endl;
    }
    const PS::S32 n_tmp = 10;
    FPSoft * ptcl_loc = new FPSoft[n_tmp];
    FPSoft * ptcl_glb = new FPSoft[n_tmp];
    PS::S32 * n_ptcl_array = new PS::S32[n_proc];
    PS::S32 * n_ptcl_disp = new PS::S32[n_proc+1];
    for(PS::S32 i=0; i<(n_glb-1)/n_tmp+1; i++){
	PS::S32 i_head = i * n_tmp;
	PS::S32 i_tail = (i+1) * n_tmp;
	PS::S32 n_cnt = 0;
	for(PS::S32 j=0; j<n_loc; j++){
	    if( i_head<=psys[j].id && psys[j].id<i_tail){
		ptcl_loc[n_cnt].id = psys[j].id;
		ptcl_loc[n_cnt].mass = psys[j].mass;
		ptcl_loc[n_cnt].pos = psys[j].pos;
		ptcl_loc[n_cnt].vel = psys[j].vel;
		n_cnt++;
	    }
	}
	PS::Comm::allGather(&n_cnt, 1, n_ptcl_array);
	n_ptcl_disp[0] = 0;
	for(PS::S32 j=0; j<n_proc; j++){
	    n_ptcl_disp[j+1] = n_ptcl_disp[j] + n_ptcl_array[j];
	}
	PS::Comm::gatherV(ptcl_loc, n_cnt, ptcl_glb, n_ptcl_array, n_ptcl_disp);
	if(my_rank == 0){
	    const PS::S32 n = n_ptcl_disp[n_proc];
	    //std::sort(ptcl_glb, ptcl_glb+n, SortID());
	    std::sort(ptcl_glb, ptcl_glb+n,
		      [](const FPSoft & left, const FPSoft & right)->bool{return left.id < right.id;}
		      );
	    for(PS::S32 j=0; j<n; j++){
		fout<<ptcl_glb[j].id<<"   "<<ptcl_glb[j].mass<<"   "<<ptcl_glb[j].pos<<"   "<<ptcl_glb[j].vel<<std::endl;
		//fout<<ptcl_glb[j].mass<<"   "<<ptcl_glb[j].pos<<"   "<<ptcl_glb[j].vel<<std::endl;
	    }
	}
    }
    fout.close();
}
