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
                    const int seed = 0,
                    const PS::F64 phi_min = 0.0,
                    PS::F64 phi_max = 0.0){
    static const double PI = atan(1.0) * 4.0;
    if(phi_max == 0.0) phi_max = 2.0 * PI;
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
        //double l = 2.0 * PI * mt.genrand_res53(); // mean anomayl
        double l = ((phi_max-phi_min)*mt.genrand_res53()) + phi_min; // mean anomayl
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
                           const PS::S32 seed = 0,
                           const PS::F64 phi_min = 0.0,
                           PS::F64 phi_max = 0.0){
    if(phi_max == 0.0) phi_max = 2.0 * atan(1.0) * 4.0;
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    //auto n_loc = n_glb;
    PS::S64 n_loc = n_glb;
#if 0 // for debug
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
                   ax_in, ax_out, ecc_rms, inc_rms, dens, mass_sun, a_ice, f_ice, power, seed,
                   phi_min, phi_max);
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

PS::F64vec getVel(const PS::F64vec & pos){
    const double PI = atan(1.0) * 4.0;
    const PS::F64 dr = sqrt(pos*pos);
    const PS::F64 v_kep = sqrt(1.0/(dr));
    const PS::F64 theta = atan2(pos.y, pos.x) + PI*0.5;
    const PS::F64 vx = v_kep * cos(theta);
    const PS::F64 vy = v_kep * sin(theta);
    const PS::F64 vz = 0.0;
    return PS::F64vec(vx, vy, vz);
}

template<class Tpsys>
void SetParticleKeplerDisk2(Tpsys & psys,
                            const PS::S64 n_glb,
                            const PS::F64 ax_in,
                            const PS::F64 ax_out,
                            const PS::F64 dens, // surface density
                            const PS::F64 dz, // hight thickeness is 2*dz
                            const PS::F64ort box) // domain
{
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok1" << std::endl;
    const PS::F64ort boundary(PS::F64vec(-2.0, -2.0, -2.0), PS::F64vec(2.0, 2.0, 2.0));
    const PS::F64 pos_x_low = (box.low_.x > boundary.low_.x) ? box.low_.x : boundary.low_.x;
    const PS::F64 pos_x_high = (box.high_.x < boundary.high_.x) ? box.high_.x : boundary.high_.x;
    const PS::F64 pos_y_low = (box.low_.y > boundary.low_.y) ? box.low_.y : boundary.low_.y;
    const PS::F64 pos_y_high = (box.high_.y < boundary.high_.y) ? box.high_.y : boundary.high_.y;
    
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const double PI = atan(1.0) * 4.0;
    const PS::F64 area = PI*(ax_out*ax_out-ax_in*ax_in);
    const PS::F64 mass = (dens*area) / n_glb; // particle mass
    const PS::F64 dx = sqrt(area / (n_glb / 3)); // grid size
    const PS::F64 dy = dx;
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok1-a" << std::endl;
    //if (PS::Comm::getRank() == 0) {
    //    std::cout << "dx = " << dx << std::endl;
    //    std::cout << "dy = " << dy << std::endl;
    //    std::cout << "box.low_.x = " << box.low_.x << std::endl;
    //    std::cout << "box.high_.x = " << box.high_.x << std::endl;
    //    std::cout << "box.low_.y = " << box.low_.y << std::endl;
    //    std::cout << "box.high_.y = " << box.high_.y << std::endl;
    //}
    //PS::Comm::barrier(); 
    //PS::S64 x_low  = (PS::S64)(box.low_.x/dx) - 2;
    //PS::S64 x_high = (PS::S64)(box.high_.x/dx) + 2;
    //PS::S64 y_low  = (PS::S64)(box.low_.y/dy) - 2;
    //PS::S64 y_high = (PS::S64)(box.high_.y/dy) + 2;
    PS::S64 x_low  = (PS::S64)(pos_x_low/dx) - 2;
    PS::S64 x_high = (PS::S64)(pos_x_high/dx) + 2;
    PS::S64 y_low  = (PS::S64)(pos_y_low/dy) - 2;
    PS::S64 y_high = (PS::S64)(pos_y_high/dy) + 2;

    
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok1-b" << std::endl;
    PS::ReallocatableArray<PS::F64vec> pos;
    pos.clearSize();
    pos.reserve(n_glb/PS::Comm::getNumberOfProc());
    const PS::F64 ax_out_sq = ax_out*ax_out;
    const PS::F64 ax_in_sq = ax_in*ax_in;
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok2" << std::endl;
    //if (PS::Comm::getRank() == 0) {
    ////if (PS::Comm::getRank() == 5) {
    //     std::cout << "x_low  = " << x_low << std::endl;
    //     std::cout << "x_high = " << x_high << std::endl;
    //     std::cout << "y_low  = " << y_low << std::endl;
    //     std::cout << "y_high  = " << y_high << std::endl;
    //}
    for(PS::S64 ix=x_low; ix<=x_high; ix++){
        const PS::F64 x = ix*dx;
        const PS::F64 x2 = x+0.5*dx;
        for(PS::S64 iy=y_low; iy<=y_high; iy++){
            const PS::F64 y = iy*dy;
            const PS::F64 y2 = y+0.5*dy;
            const PS::F64vec pos_cen(x, y, 0.0);
            const PS::F64vec pos_up(x2, y2, dz);
            const PS::F64vec pos_dw(x2, y2, -dz);
            //std::cerr<<"x= "<<x<<" y= "<<y<<" r= "<<sqrt(x*x+y*y)<<std::endl;
            //std::cerr<<"ax_in= "<<ax_in<<" ax_out= "<<ax_out<<std::endl;
            if( (x*x+y*y) >= ax_in_sq && (x*x+y*y) < ax_out_sq && box.contained(pos_cen)){
                pos.push_back(pos_cen);
            }
            if( (x2*x2+y2*y2) >= ax_in_sq && (x2*x2+y2*y2) < ax_out_sq && box.contained(pos_up)){
                pos.push_back(pos_up);
            }
            if( (x2*x2+y2*y2) >= ax_in_sq && (x2*x2+y2*y2) < ax_out_sq && box.contained(pos_dw)){
                pos.push_back(pos_dw);
            }            
        }
    }
    //std::cout << "myrank = " << PS::Comm::getRank() << std::endl;
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok3" << std::endl;
    PS::S32 n_loc = pos.size();
    psys.setNumberOfParticleLocal(n_loc);
    for(PS::S32 i=0; i<n_loc; i++){
        psys[i].pos  = pos[i];
        psys[i].vel  = getVel(psys[i].pos);
        psys[i].mass = mass;
    }
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok4" << std::endl;
}

template<class Tpsys>
void SetParticleKeplerDisk3(Tpsys & psys,
                            const PS::S64 n_glb,
                            const PS::F64 ax_in,
                            const PS::F64 ax_out,
                            const PS::F64 dens, // surface density at ax_in
                            const PS::F64ort box, // domain
                            const bool layer=true)
{
    assert(PS::Comm::getNumberOfProc()%4 == 0);
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok1" << std::endl;
    const double PI = atan(1.0) * 4.0;
    const PS::F64 area = PI*(ax_out*ax_out-ax_in*ax_in);
    PS::F64 dS = (area*3.0) / n_glb;
    if(!layer){
        dS = area / n_glb;
    }
    //std::cerr<<"area= "<<area<<" dS= "<<dS<<std::endl;
    const PS::F64 mass = (dens*area) / n_glb; // particle mass
    const PS::F64 dr = sqrt(dS);
    const PS::F64 ax_cen = (ax_in+ax_out) * 0.5;
    const PS::F64 dtheta = sqrt(dS) / ax_cen;

    //std::cerr<<"dr= "<<dr<<std::endl;
    //std::cerr<<"dtheta= "<<dtheta<<std::endl;
    //std::cerr<<"ax_cen*dtheta= "<<ax_cen*dtheta<<std::endl;

    const PS::F64ort boundary(PS::F64vec(-ax_out*1.00001, -ax_out*1.00001, -1000000), PS::F64vec(ax_out*1.00001, ax_out*1.00001, 1000000.0));
    PS::F64 pos_x_low = (box.low_.x > boundary.low_.x) ? box.low_.x : boundary.low_.x;
    PS::F64 pos_x_high = (box.high_.x < boundary.high_.x) ? box.high_.x : boundary.high_.x;
    PS::F64 pos_y_low = (box.low_.y > boundary.low_.y) ? box.low_.y : boundary.low_.y;
    PS::F64 pos_y_high = (box.high_.y < boundary.high_.y) ? box.high_.y : boundary.high_.y;

    /*
    std::cerr<<"pos_x_low= "<<pos_x_low<<" pos_x_high= "<<pos_x_high
             <<" pos_y_low= "<<pos_y_low<<" pos_y_high= "<<pos_y_high
             <<std::endl;
    */

    
    PS::F64 r_max = -99999999.9;
    PS::F64 r_min =  99999999.9;
    PS::F64 theta_max = 0.0;
    PS::F64 theta_min = 0.0;
    PS::F64 dis[4];
    dis[0] = sqrt(pos_x_low*pos_x_low + pos_y_low*pos_y_low);
    dis[1] = sqrt(pos_x_high*pos_x_high + pos_y_low*pos_y_low);
    dis[2] = sqrt(pos_x_high*pos_x_high + pos_y_high*pos_y_high);
    dis[3] = sqrt(pos_x_low*pos_x_low + pos_y_high*pos_y_high);    
    for(PS::S32 i=0; i<4; i++){
        if(dis[i] > r_max){
            r_max = dis[i];
        }
        if(dis[i] < r_min){
            r_min = dis[i];
        }
    }
    //r_max = sqrt(std::min(r_max*r_max, ax_out*ax_out));
    //r_min = sqrt(std::max(r_min*r_min, ax_in*ax_in));
    r_max = std::min(r_max, ax_out);
    r_min = std::max(r_min, ax_in);    

    //std::cerr<<"r_max= "<<r_max<<" r_min= "<<r_min<<std::endl;

    PS::S64 ir_min = r_min / dr - 2;
    PS::S64 ir_max = r_max / dr + 2;

    //std::cerr<<"ir_min= "<<ir_min<<" ir_max= "<<ir_max<<std::endl;
    PS::ReallocatableArray<PS::F64vec> pos;
    pos.reserve(1000);
    for(PS::S64 ir=ir_min; ir<ir_max; ir++){
        PS::F64 r_tmp = dr*ir;
        //std::cerr<<"r_tmp= "<<r_tmp<<std::endl;
        if(r_tmp <= ax_in ||  r_tmp >= ax_out) continue;
        //if(r_tmp*r_tmp <= ax_in*ax_in ||  r_tmp*r_tmp >= ax_out*ax_out) continue;
        //std::cerr<<"ir= "<<ir<<std::endl;
        //std::cerr<<"r_tmp= "<<r_tmp<<" pos_x_high= "<<pos_x_high<<" r_tmp= "<<r_tmp<<std::endl;

        PS::F64 cos_th0 = pos_x_high / r_tmp;
        PS::F64 cos_th1 = pos_x_low / r_tmp;
        PS::F64 sin_th2 = pos_y_high / r_tmp;
        PS::F64 sin_th3 = pos_y_low / r_tmp;
        /*
        std::cerr<<"cos_th0= "<<cos_th0
                 <<" cos_th1= "<<cos_th1
                 <<" sin_th2= "<<sin_th2
                 <<" sin_th3= "<<sin_th3
                 <<std::endl;
        */
        cos_th0 = std::min(cos_th0, 1.0);
        cos_th0 = std::max(cos_th0, -1.0);

        cos_th1 = std::min(cos_th1, 1.0);
        cos_th1 = std::max(cos_th1, -1.0);

        sin_th2 = std::min(sin_th2, 1.0);
        sin_th2 = std::max(sin_th2, -1.0);

        sin_th3 = std::min(sin_th3, 1.0);
        sin_th3 = std::max(sin_th3, -1.0);
        /*
        std::cerr<<"cos_th0= "<<cos_th0
                 <<" cos_th1= "<<cos_th1
                 <<" sin_th2= "<<sin_th2
                 <<" sin_th3= "<<sin_th3
                 <<std::endl;
        */
        PS::F64 th_max, th_min;
        PS::F64 th_min_c, th_max_c, th_min_s, th_max_s;        
        if( sin_th2 > 0.0 || sin_th3 > 0.0){
            if(cos_th0 > 0.0 || cos_th1 > 0.0){
                //std::cerr<<"1st"<<std::endl;
                PS::F64 th0 = acos(cos_th0);
                PS::F64 th1 = acos(cos_th1);
                PS::F64 th2 = asin(sin_th2);
                PS::F64 th3 = asin(sin_th3);
                //std::cerr<<"th0= "<<th0<<" th1= "<<th1<<" th2= "<<th2<<" th3= "<<th3<<std::endl;
                if(th0 > th1){
                    th_max_c = th0;
                    th_min_c = th1;
                }
                else{
                    th_min_c = th0;
                    th_max_c = th1;
                }
                if(th2 > th3){
                    th_max_s = th2;
                    th_min_s = th3;
                }
                else{
                    th_min_s = th2;
                    th_max_s = th3;
                }
            }
            else if(cos_th0 < 0.0 || cos_th1 < 0.0){
                // 2nd
                //std::cerr<<"2nd"<<std::endl;
                PS::F64 th0 = acos(cos_th0);
                PS::F64 th1 = acos(cos_th1);
                PS::F64 th2 = PI - asin(sin_th2);
                PS::F64 th3 = PI - asin(sin_th3);
                //std::cerr<<"th0= "<<th0<<" th1= "<<th1<<" th2= "<<th2<<" th3= "<<th3<<std::endl;
                if(th0 > th1){
                    th_max_c = th0;
                    th_min_c = th1;
                }
                else{
                    th_min_c = th0;
                    th_max_c = th1;                    
                }
                if(th2 > th3){
                    th_max_s = th2;
                    th_min_s = th3;
                }
                else{
                    th_min_s = th2;
                    th_max_s = th3;
                }
            }
            else{
                std::cerr<<"ERROR0"<<std::endl;
                exit(1);
            }
            th_max = std::min(th_max_s, th_max_c);
            th_min = std::max(th_min_s, th_min_c);            
        }
        else if(sin_th2 < 0.0 || sin_th3 < 0.0){
            if(cos_th0 < 0.0 || cos_th1 < 0.0){
                //std::cerr<<"3rd"<<std::endl;
                // 3rd
                //PS::F64 th0 = 2.0*PI - acos(cos_th0);
                //PS::F64 th1 = 2.0*PI - acos(cos_th1);
                //PS::F64 th2 = PI - asin(sin_th2); 
                //PS::F64 th3 = PI - asin(sin_th3);
                PS::F64 th0 = 2.0*PI - acos(cos_th0);
                PS::F64 th1 = 2.0*PI - acos(cos_th1);
                PS::F64 th2 = PI - asin(sin_th2); 
                PS::F64 th3 = PI - asin(sin_th3);
                //std::cerr<<"th0= "<<th0<<" th1= "<<th1<<" th2= "<<th2<<" th3= "<<th3<<std::endl;
                if(th0 > th1){
                    th_max_c = th0;
                    th_min_c = th1;
                }
                else{
                    th_min_c = th0;
                    th_max_c = th1;                    
                }
                if(th2 > th3){
                    th_max_s = th2;
                    th_min_s = th3;
                }
                else{
                    th_min_s = th2;
                    th_max_s = th3;
                }
            }
            else if(cos_th0 > 0.0 || cos_th1 > 0.0){
                // 4th
                //std::cerr<<"4th"<<std::endl;
                PS::F64 th0 = 2.0*PI - acos(cos_th0);
                PS::F64 th1 = 2.0*PI - acos(cos_th1);
                PS::F64 th2 = 2.0*PI + asin(sin_th2);
                PS::F64 th3 = 2.0*PI + asin(sin_th3);
                if(th0 > th1){
                    th_max_c = th0;
                    th_min_c = th1;
                }
                else{
                    th_min_c = th0;
                    th_max_c = th1;                    
                }
                if(th2 > th3){
                    th_max_s = th2;
                    th_min_s = th3;
                }
                else{
                    th_min_s = th2;
                    th_max_s = th3;
                }
            }
            else{
                std::cerr<<"ERROR1"<<std::endl;
                exit(1);
            }
            th_max = std::min(th_max_s, th_max_c);
            th_min = std::max(th_min_s, th_min_c);
        }
        else{
            std::cerr<<"ERROR2"<<std::endl;
            exit(1);
        }
        
        //std::cerr<<"th_min= "<<th_min<<" th_max= "<<th_max<<std::endl;
        PS::S64 it_min = (th_min / dtheta) - 2;
        PS::S64 it_max = (th_max / dtheta) + 2;
        //std::cerr<<"it_min= "<<it_min<<" it_max= "<<it_max<<std::endl;
        for(PS::S64 it=it_min; it<it_max; it++){
            const PS::F64 theta_tmp = dtheta*it;
            const PS::F64 x_tmp = r_tmp*cos(theta_tmp);
            const PS::F64 y_tmp = r_tmp*sin(theta_tmp);
            const PS::F64vec pos_cen(x_tmp, y_tmp, 0.0);
            //std::cerr<<"sqrt(pos_cen*pos_cen)= "<<sqrt(pos_cen*pos_cen)<<std::endl;
            //if( box.contained(pos_cen))  pos.push_back(pos_cen);
            if( r_tmp >= ax_in && r_tmp < ax_out && box.contained(pos_cen))  pos.push_back(pos_cen);
            if(layer){
                const PS::F64 x_tmp2 = (r_tmp+dr*0.5)*cos(theta_tmp+0.5*dtheta);
                const PS::F64 y_tmp2 = (r_tmp+dr*0.5)*sin(theta_tmp+0.5*dtheta);
                const PS::F64vec pos_up(x_tmp2, y_tmp2, dr*0.5);
                const PS::F64vec pos_dw(x_tmp2, y_tmp2, -dr*0.5);
                if( r_tmp+dr*0.5 >= ax_in && r_tmp+dr*0.5 < ax_out && box.contained(pos_up))  pos.push_back(pos_up);
                if( r_tmp+dr*0.5 >= ax_in && r_tmp+dr*0.5 < ax_out && box.contained(pos_dw))  pos.push_back(pos_dw);                
                //if( box.contained(pos_up))  pos.push_back(pos_up);
                //if( box.contained(pos_dw))  pos.push_back(pos_dw);
            }
        }
    }
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok2" << std::endl;
    PS::S64 n_loc = pos.size();
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "n_loc = " << n_loc << std::endl;
    psys.setNumberOfParticleLocal(n_loc);
    for(PS::S64 i=0; i<n_loc; i++){
        psys[i].pos  = pos[i];
        psys[i].vel  = getVel(psys[i].pos);
        psys[i].mass = mass;
    }
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "ok3" << std::endl;
}


template<class Tpsys, class Tsat>
void ChangePtclToSatellite(Tpsys & sys,
                           Tsat & sat,
                           const PS::S32 n_sat){
    //const PS::F64 mass_factor = 100.0;
    const PS::F64 mass_factor = 1.0;
    sat.resizeNoInitialize(n_sat);
    Satellite::r_coll = Epj::r_coll * pow(mass_factor, 1.0/3.0); // added by DN
    if(PS::Comm::getRank() == 0){
        PS::S32 * adr_remove = new PS::S32[n_sat];
        /*
        for(PS::S32 i=0; i<n_sat; i++){
            adr_remove[i] = i;
            SATELLITE[i].mass = sys[i].mass * mass_factor;
#if 0
            SATELLITE[i].pos  = sys[i].pos*1e12;
            SATELLITE[i].vel  = sys[i].vel*1e6;
#else
            SATELLITE[i].pos  = sys[i].pos;
            SATELLITE[i].vel  = sys[i].vel;
#endif
            SATELLITE[i].id   = i; // added by DN
            //SATELLITE[i].r_coll = Epj::r_coll * pow(mass_factor, 1.0/3.0);
        }
        */
        PS::S32 j = 0;
        PS::S32 n_loc = sys.getNumberOfParticleLocal();
        const PS::S32 stride = n_loc / n_sat;
        for(PS::S32 i=0, j=0; i<n_sat; i++, j+=stride){
            adr_remove[i] = i;
            SATELLITE[i].mass = sys[i].mass * mass_factor;
#if 0
            SATELLITE[i].pos  = sys[i].pos*1e12;
            SATELLITE[i].vel  = sys[i].vel*1e6;
#else
            SATELLITE[i].pos  = sys[i].pos;
            SATELLITE[i].vel  = sys[i].vel;
#endif
            SATELLITE[i].id   = i; // added by DN
            //SATELLITE[i].r_coll = Epj::r_coll * pow(mass_factor, 1.0/3.0);
        }
        
        sys.removeParticle(adr_remove, n_sat);
        delete [] adr_remove;
    }
    PS::Comm::broadcast(SATELLITE.getPointer(), n_sat);
}

template<class Tpsys>
void CollisionTest_IC(Tpsys & psys, 
                      const PS::S64 n_glb, 
                      const PS::F64 mass,
                      const PS::F64 ax,
                      const PS::F64 dz, 
                      const PS::S32 mode) {
   assert(n_glb == 2);
   PS::S32 myrank = PS::Comm::getRank();

   switch (mode) {
       case 0:
           if (myrank == 0) {
               psys.setNumberOfParticleLocal(n_glb);
               //* 1st particle
               psys[0].pos.x = - 0.5*dz;
               psys[0].pos.y = 0.0;
               psys[0].pos.z = 0.0;
               psys[0].vel.x = 0.0; 
               psys[0].vel.y = 0.0;
               psys[0].vel.z = 0.0;
               //* 2nd particle
               psys[1].pos.x =   0.5*dz;
               psys[1].pos.y = 0.0;
               psys[1].pos.z = 0.0;
               psys[1].vel.x = - dz;
               psys[1].vel.y = 0.0;
               psys[1].vel.z = 0.0;
               // Here, we assume that the velocity of head-on collision
               // is the order of shear velocity at r = 1, where the angular
               // velocity, \Omega, is 1 in our unit. The shear velocity is
               // roughly \Omega * dz.
               //* Set mass,vel,id
               for (PS::S32 i=0; i<n_glb; i++) {
                  psys[i].mass = mass;
                  psys[i].id   = i;
               }
           } else {
               psys.setNumberOfParticleLocal(0);
           }
           break;
       case 1:
           if (myrank == 0) {
               psys.setNumberOfParticleLocal(n_glb);
               //* 1st particle
               psys[0].pos.x = ax;
               psys[0].pos.y = 0.0;
               psys[0].pos.z = - dz;
               //* 2nd particle
               psys[1].pos.x = ax;
               psys[1].pos.y = 0.0;
               psys[1].pos.z = dz;
               //* Set mass,vel,id
               for (PS::S32 i=0; i<n_glb; i++) {
                  psys[i].mass = mass;
                  psys[i].vel  = getVel(psys[i].pos);
                  psys[i].id   = i;
               }
           } else {
               psys.setNumberOfParticleLocal(0);
           }
           break;
       default:
           std::cout << "Invalid mode specified!!" << std::endl;
           PS::Finalize();
           std::exit(1);
   }
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
#if 1
    if(my_rank == 0) n_loc = n_glb;
    else n_loc = 0;
#else
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
#endif
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
void SetParticleUniformSphere(Tpsys & psys,
                              const PS::S64 n_glb,
                              PS::S32 & n_loc,  
                              PS::F32 & t_sys){
    //* Compute # of local particles (n_loc)
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
#if 1
    if(my_rank == 0) n_loc = n_glb;
    else n_loc = 0;
#else
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
#endif
    psys.setNumberOfParticleLocal(n_loc);
    //* Memory allocation and initialize t_sys
    PS::F64    * mass = new PS::F64[n_loc];
    PS::F64vec * pos  = new PS::F64vec[n_loc];
    PS::F64vec * vel  = new PS::F64vec[n_loc];
    t_sys = 0.0;
    //* Make uniform sphere
    const PS::F64 mass_glb = 1.0e-3; // much smaller than the planes mass (=1).
    const PS::F64 radius = 1.0;
    PS::MTTS mt;
    mt.init_genrand(0);
    for (PS::S32 i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        do {
            pos[i][0] = (2. * mt.genrand_res53() - 1.) * radius;
            pos[i][1] = (2. * mt.genrand_res53() - 1.) * radius;
            pos[i][2] = (2. * mt.genrand_res53() - 1.) * radius;
        } while (pos[i] * pos[i] >= radius * radius);
        vel[i][0] = 0.0;
        vel[i][1] = 0.0;
        vel[i][2] = 0.0;
    }
    PS::F64vec cm_pos  = 0.0;
    PS::F64vec cm_vel  = 0.0;
    PS::F64    cm_mass = 0.0;
    for(PS::S32 i = 0; i < n_loc; i++){
        cm_pos  += mass[i] * pos[i];
        cm_vel  += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(PS::S32 i = 0; i < n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }
    //* Set to ParticleSystem                   
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        psys[i].pos = pos[i];
        psys[i].mass = mass[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

template<class Tsys, class Tsat>
void RemovePartcile(Tsys & sys, const Tsat & satellite){
    std::vector<PS::S32> adr_remove;
    adr_remove.clear();
    const PS::S32 n_loc = sys.getNumberOfParticleLocal();
    const PS::S32 n_sat = satellite.size();
    for(PS::S32 i=0; i<n_loc; i++){
        for(PS::S32 j=0; j<n_sat; j++){
            PS::F64vec rij = sys[i].pos - satellite[j].pos;
            if(rij*rij < satellite[j].r_coll*satellite[j].r_coll){
                adr_remove.push_back(i);
            }
        }
    }
    sys.removeParticle(&adr_remove[0], adr_remove.size());
    std::cerr<<"adr_remove.size()= "<<adr_remove.size()<<std::endl;
}

/*
template<class Tsys, class Tsat>
void RemovePartcile(Tsys & sys, const Tsat * satellite, const PS::S32 n_sat){
    std::vector<PS::S32> adr_remove;
    adr_remove.clear();
    const PS::S32 n_loc = sys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        for(PS::S32 j=0; j<n_sat; j++){
            PS::F64vec rij = sys[i].pos - satellite[j].pos;
            if(rij*rij < satellite[j].r_coll*satellite[j].r_coll){
                adr_remove.push_back(i);
            }
        }
    }
    sys.removeParticle(&adr_remove[0], adr_remove.size());
    std::cerr<<"adr_remove.size()= "<<adr_remove.size()<<std::endl;
}
*/

/*
template<class Tsys, class Tsat>
void RemovePartcile(Tsys & sys){
    std::vector<PS::S32> adr_remove;
    adr_remove.clear();
    const PS::S32 n_loc = sys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        if(sys[i].n_coll > )
        for(PS::S32 j=0; j<n_sat; j++){
            PS::F64vec rij = sys[i].pos - satellite[j].pos;
            if(rij*rij < satellite[j].r_coll*satellite[j].r_coll){
                adr_remove.push_back(i);
            }
        }
    }
    sys.removeParticle(&adr_remove[0], adr_remove.size());
    std::cerr<<"adr_remove.size()= "<<adr_remove.size()<<std::endl;
}
*/

#if 0
void MakeSimpleRing(PS::F64vec & pos,
                    PS::F64vec & vel,
                    const PS::F64    & mass_pla,
                    const PS::F64vec & pos_pla,
                    const PS::F64vec & vel_pla,
                    const double ax_out,
                    const double ax_in,
                    const double hight,
                    const double tau){
    static PS::MTTS mt;
    static bool first = true;
    static const PS::F64 PI = atan(1.0) * 4.0;
    if(kfirst){
        mt.init_genrand(PS::Comm::getRank());
        first = false;
    }
    double phi_tmp = mt.genrand_res53()*2.0*PI;
    double z_tmp = (mt.genrand_res53()-0.5)*2.0*hight;
    pos.z = z_tmp;
    pos.x = cos(phi_tmp);
    pos.y = sin(phi_tmp);
}

template<class Tpsys>
void SetParticleSimpleRing(Tpsys & psys,
                           const PS::S64 n_glb,
                           PS::S32 & n_loc,  
                           PS::F32 & t_sys){
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
#if 1
    // for debug
    if(my_rank == 0) n_loc = n_glb;
    else n_loc = 0;
#else
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
#endif
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
#endif

template<class Tpsys>
void SetParticleKeplerDiskCyl(Tpsys & psys,
                              const PS::S64 n_glb,
                              const PS::F64 ax_in,
                              const PS::F64 ax_out,
                              const PS::F64 dens, // surface density at ax_in
                              const PS::F64ort box, // domain
                              const bool layer=true){
    const double PI = atan(1.0) * 4.0;
    const PS::F64 area = PI*(ax_out*ax_out-ax_in*ax_in);
    PS::F64 dS = (area*3.0) / n_glb;
    if(!layer){
        dS = area / n_glb;
    }
    const PS::F64 mass = (dens*area) / n_glb; // particle mass
    const PS::F64 dr = sqrt(dS);
    const PS::F64 ax_cen = (ax_in+ax_out) * 0.5;
    const PS::F64 dtheta = sqrt(dS) / ax_cen;
    const PS::S64 id_x_head = box.low_.x / dr;
    const PS::S64 id_x_tail = box.high_.x / dr;
    const PS::S64 id_y_head = box.low_.y / dr;
    const PS::S64 id_y_tail = box.high_.y / dr;
    PS::ReallocatableArray<PS::F64vec> pos;
    pos.clearSize();
    for(PS::S64 i=id_x_head-2; i<=id_x_tail+2; i++){
        PS::F64 pos_x = i*dr;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-2; j<=id_y_tail+2; j++){
            PS::F64 pos_y = j*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
            pos.push_back(PS::F64vec(pos_x, pos_y, 0.0));
        }
    }
    for(PS::S64 i=id_x_head-2; i<=id_x_tail+2; i++){
        PS::F64 pos_x = i*dr+0.5*dr;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-2; j<=id_y_tail+2; j++){
            PS::F64 pos_y = j*dr+0.5*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
            pos.push_back(PS::F64vec(pos_x, pos_y, -0.5*dr));
            pos.push_back(PS::F64vec(pos_x, pos_y, 0.5*dr));
        }
    }
    PS::S64 n_loc = pos.size();
    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0) std::cerr<<"OK 01 @SetParticleKeplerDiskCyl"<<std::endl;
    psys.setNumberOfParticleLocal(n_loc);
    for(PS::S64 i=0; i<n_loc; i++){
        PS::F64 pos_x = pos[i].y*cos(pos[i].x);
        PS::F64 pos_y = pos[i].y*sin(pos[i].x);
        PS::F64 pos_z = pos[i].z;
        psys[i].pos.x  = pos_x;
        psys[i].pos.y  = pos_y;
        psys[i].pos.z  = pos_z;
        psys[i].vel  = getVel(psys[i].pos);
        psys[i].mass = mass;
    }
    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0) std::cerr<<"OK 02 @SetParticleKeplerDiskCyl"<<std::endl;
}

template<class Tpsys>
void SetParticleKeplerDiskCyl2(Tpsys & psys,
                               const PS::S64 n_glb,
                               const PS::F64 ax_in,
                               const PS::F64 ax_out,
                               const PS::F64 dens, // surface density at ax_in
                               const PS::F64ort box, // domain
                               const PS::F64 r_phy, // radius of particle
                               const bool layer=true){
    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    
    const double PI = atan(1.0) * 4.0;
    const PS::F64 area = PI*(ax_out*ax_out-ax_in*ax_in);
    PS::F64 dS = (area*3.0) / n_glb;
    if(!layer){
        dS = area / n_glb;
    }
    const PS::F64 mass = (dens*area) / n_glb; // particle mass
    const PS::F64 dl = sqrt(dS);
    const PS::F64 dz = dl;
    const PS::F64 ax_cen = (ax_in+ax_out) * 0.5;
    const PS::S64 n_theta = (2.0*PI*ax_cen / dl);
    const PS::F64 dtheta = 2.0*PI / n_theta;
    const PS::F64 dax = ax_out - ax_in;
    const PS::S64 n_r = (dax / dl);
    const PS::F64 dr = dax / n_r;

    const PS::S64 id_x_head = box.low_.x / dtheta;
    const PS::S64 id_x_tail = box.high_.x / dtheta;
    const PS::S64 id_y_head = box.low_.y / dr;
    const PS::S64 id_y_tail = box.high_.y / dr;
    const PS::F64 offset_theta = dtheta*(sqrt(2.0)-1.0)*0.5;

    const PS::F64 eps = (dl > 2.0*r_phy) ? (0.5*(dl-2.0*r_phy))*0.9 : 0.0;
    
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"n_theta= "<<n_theta
                 <<" dtheta= "<<dtheta
                 <<" offset_theta= "<<offset_theta
                 <<std::endl;
    }
    PS::Comm::barrier();
    
    PS::TempArray<PS::F64vec> pos;
    pos.reserve(pos.getLeftSizeOfMemoryPool()-100);
    for(PS::S64 i=id_x_head-3; i<=id_x_tail+3; i++){
        PS::F64 pos_x = ((PS::F64)i)*dtheta + offset_theta;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-3; j<=id_y_tail+3; j++){
            PS::F64 pos_y = j*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
#if 1
            pos.pushBackNoCheck(PS::F64vec(pos_x, pos_y, 0.0));
#else
            PS::F64 eps_x = eps*(mt.genrand_res53()-0.5)*2.0;
            PS::F64 eps_y = eps*(mt.genrand_res53()-0.5)*2.0;
            PS::F64 eps_z = eps*(mt.genrand_res53()-0.5)*2.0;
            pos.pushBackNoCheck(PS::F64vec(pos_x+eps_x, pos_y+eps_y, eps_z));
#endif
        }
    }
    for(PS::S64 i=id_x_head-3; i<=id_x_tail+3; i++){
        PS::F64 pos_x = ((PS::F64)i)*dtheta+0.5*dtheta + offset_theta;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-3; j<=id_y_tail+3; j++){
            PS::F64 pos_y = j*dr+0.5*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
#if 1
            pos.pushBackNoCheck(PS::F64vec(pos_x, pos_y, -0.5*dz));
#else
            PS::F64 eps_x = eps*(mt.genrand_res53()-0.5)*2.0;
            PS::F64 eps_y = eps*(mt.genrand_res53()-0.5)*2.0;
            PS::F64 eps_z = eps*(mt.genrand_res53()-0.5)*2.0;
            pos.pushBackNoCheck(PS::F64vec(pos_x+eps_x, pos_y+eps_y, -0.5*dz+eps_z));
#endif
        }
    }
    for(PS::S64 i=id_x_head-3; i<=id_x_tail+3; i++){
        PS::F64 pos_x = ((PS::F64)i)*dtheta+0.5*dtheta + offset_theta;
        if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
        for(PS::S64 j=id_y_head-3; j<=id_y_tail+3; j++){
            PS::F64 pos_y = j*dr+0.5*dr;
            if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
#if 1
            pos.pushBackNoCheck(PS::F64vec(pos_x, pos_y, 0.5*dz));
#else
            PS::F64 eps_x = eps*(mt.genrand_res53()-0.5)*2.0;
            PS::F64 eps_y = eps*(mt.genrand_res53()-0.5)*2.0;
            PS::F64 eps_z = eps*(mt.genrand_res53()-0.5)*2.0;
            pos.pushBackNoCheck(PS::F64vec(pos_x+eps_x, pos_y+eps_y, 0.5*dz+eps_z));
#endif
        }
    }
    PS::S64 n_loc = pos.size();
    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0) std::cerr<<"OK 01 @SetParticleKeplerDiskCyl2"<<std::endl;
    psys.setNumberOfParticleLocal(n_loc);

    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0) std::cerr<<"OK 02 @SetParticleKeplerDiskCyl2"<<std::endl;
    
    for(PS::S64 i=0; i<n_loc; i++){
        PS::F64 pos_x = pos[i].y*cos(pos[i].x);
        PS::F64 pos_y = pos[i].y*sin(pos[i].x);
        PS::F64 pos_z = pos[i].z;
        psys[i].pos.x  = pos_x;
        psys[i].pos.y  = pos_y;
        psys[i].pos.z  = pos_z;
        psys[i].vel  = getVel(psys[i].pos);
        psys[i].mass = mass;
    }

    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0) std::cerr<<"OK 03 @SetParticleKeplerDiskCyl2"<<std::endl;
    
    pos.free();
    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0) std::cerr<<"OK 04 @SetParticleKeplerDiskCyl2"<<std::endl;
}

