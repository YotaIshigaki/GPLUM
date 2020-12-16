#include<iostream>
#include<cstdio>
#include<unistd.h>
#include<random>
#include<particle_simulator.hpp>
#include"user_defined.hpp"

constexpr PS::F64 MY_PI = 3.14159265358979323846;

FP_t PLANET;

void PrintHelp() {
    std::cerr<<"a: ring width (default: 1e-3)"<<std::endl;
    std::cerr<<"d: inverse of dt (default: 256)"<<std::endl;
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"L: # of steps sharing the same list (default: 1)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: # of particles per process (default: 1024)"<<std::endl;
    std::cerr<<"e: exchange LET mode"<<std::endl;
    std::cerr<<"   0: PS::EXCHANGE_LET_A2A(default)"<<std::endl;
    std::cerr<<"   1: PS::EXCHANGE_LET_P2P_EXACT"<<std::endl;
    std::cerr<<"   2: PS::EXCHANGE_LET_P2P_FAST"<<std::endl;
    std::cerr<<"E: restitution coefficient (default: 0.5)"<<std::endl;
    std::cerr<<"D: inverse of duration time (default: 4)"<<std::endl;
    std::cerr<<"O: optical depth (default: 1)"<<std::endl;
    std::cerr<<"R: physical radius normalized hill radius (default: 1)"<<std::endl;
    std::cerr<<"m: mid point tree method (default: off)"<<std::endl;
    std::cerr<<"r: rotating frame (default: off)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

inline PS::F64vec ConvertCar2Cyl(const PS::F64vec & pos){
    const auto pos_r   = sqrt(pos.x*pos.x + pos.y*pos.y);
    const auto pos_phi = atan2(pos.y, pos.x);
    return PS::F64vec(pos_phi, pos_r, pos.z);
}

inline PS::F64vec ConvertCyl2Car(const PS::F64vec & pos){
    const auto cth = cos(pos.x);
    const auto sth = sin(pos.x);
    const auto r = pos.y;
    const auto pos_x = r*cth;
    const auto pos_y = r*sth;
    return PS::F64vec(pos_x, pos_y, pos.z);
}

template<typename Tpsys>
void UpdateCyl(Tpsys & psys){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
      psys[i].pos_cyl = ConvertCar2Cyl(psys[i].pos_car);
    }
}


template<typename Tpsys>
void Rotate(Tpsys & psys, const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        const auto pos = psys[i].pos_car;
        const auto vel = psys[i].vel;
        const auto r = sqrt(pos*pos);
        const auto mm = sqrt(PLANET.mass/(r*r*r)); // mean mortion
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
        //const auto l_new = pos_new ^ vel_new;

        psys[i].pos_car = pos_new;
	psys[i].pos_cyl = ConvertCar2Cyl(pos_new);
        psys[i].vel     = vel_new;
    }
}

template<typename Tsys>
void RigidRotation(Tsys & sys, const PS::F64 dt){
    const auto n = sys.getNumberOfParticleLocal();
    const auto ax = 1.0;
    const auto mm = sqrt(PLANET.mass/(ax*ax*ax)); // mean mortion
    const auto theta = mm*dt;
    const auto cth = std::cos(theta);
    const auto sth = std::sin(theta);    
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        const auto x_new  = cth*sys[i].pos_car.x - sth*sys[i].pos_car.y;
        const auto y_new  = sth*sys[i].pos_car.x + cth*sys[i].pos_car.y;
        const auto vx_new = cth*sys[i].vel.x - sth*sys[i].vel.y;
        const auto vy_new = sth*sys[i].vel.x + cth*sys[i].vel.y;
        sys[i].pos_car = PS::F64vec(x_new, y_new, sys[i].pos_car.z);
	sys[i].pos_cyl = ConvertCar2Cyl(sys[i].pos_car);
        sys[i].vel     = PS::F64vec(vx_new, vy_new, sys[i].vel.z);
    }
}



void DivideNProc(PS::S32 & nx,
                 PS::S32 & ny,
                 PS::S32 & nz,
                 const PS::S32 n_proc,
                 const PS::F64 delta_ax){
    nz = 1;
    ny = 1;
    nx = n_proc / ny;
    if(n_proc == 1) return;
    const PS::F64 dx = 2.0*MY_PI   / nx;
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
        dx_tmp = 2.0*MY_PI   / nx_tmp;
        dy_tmp = delta_ax / ny_tmp;
        ratio_tmp = (dx_tmp < dy_tmp) ? dx_tmp / dy_tmp : dy_tmp / dx_tmp;
    }while( fabs(ratio_tmp-1.0) < fabs(ratio-1.0));
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
    }
    assert(n_proc == nx*ny*nz);
}

PS::F64ort GetPosDomainCyl(const PS::F64 delta_ax,
			   const PS::S32 nx,
			   const PS::S32 ny){
    constexpr PS::F64 len_x = 2.0 * MY_PI;
    PS::F64ort pos_domain;
    const auto my_rank = PS::Comm::getRank();
    const auto rank_x = my_rank / ny;
    const auto rank_y = my_rank % ny;
    const auto dx = len_x / nx;
    const auto dy = delta_ax / ny;
    const auto dy_offset = 1.0-delta_ax*0.5;
    const auto x_offset = -MY_PI;
    pos_domain.low_.x = dx*rank_x + x_offset;
    pos_domain.low_.y = dy*rank_y + dy_offset;
    pos_domain.low_.z = -MY_PI;
    pos_domain.high_.x = dx*(rank_x+1) + x_offset;
    pos_domain.high_.y = dy*(rank_y+1) + dy_offset;
    pos_domain.high_.z = MY_PI;
    return pos_domain;
}

PS::F64vec GetVel(const PS::F64vec & pos){
    const PS::F64 dr = sqrt(pos*pos);
    const PS::F64 v_kep = sqrt(1.0/(dr));
    const PS::F64 theta = atan2(pos.y, pos.x) + MY_PI*0.5;
    const PS::F64 vx = v_kep * cos(theta);
    const PS::F64 vy = v_kep * sin(theta);
    const PS::F64 vz = 0.0;
    return PS::F64vec(vx, vy, vz);
}

#ifdef QUAD
using SPJ_t    = MySPJQuadrupole;
using Moment_t = MyMomentQuadrupole;
using CalcForceSp = CalcForceSpQuad<EPI_t, SPJ_t, Force_t>;
#else
using SPJ_t    = MySPJMonopole;
using Moment_t = MyMomentMonopole;
using CalcForceSp = CalcForceSpMono<EPI_t, SPJ_t, Force_t>;
#endif

using MY_SEARCH_MODE = PS::SEARCH_MODE_LONG_SCATTER;
using Tree_t = PS::TreeForForce<MY_SEARCH_MODE, Force_t, EPI_t, EPJ_t, Moment_t, Moment_t, SPJ_t, PS::CALC_DISTANCE_TYPE_NEAREST_X>;

template<typename Tpsys>
void CalcForceFromPlanet(Tpsys & psys, const FP_t & pla){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        const auto rij = psys[i].pos_car - pla.pos_car;
        const auto r_sq = rij*rij;
        const auto r_inv = 1.0 / sqrt(r_sq);
        const auto pot   = pla.mass * r_inv;
        psys[i].acc -= pot * r_inv * r_inv * rij;
        psys[i].pot -= pot;
    }
}

template<typename Tpsys>
void Kick(Tpsys & psys, const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        psys[i].vel += psys[i].acc*dt;
        psys[i].vel_full = psys[i].vel + psys[i].acc*dt;
    }
}
template<typename Tpsys>
void Drift(Tpsys & psys,
           const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        psys[i].pos_car += psys[i].vel*dt;
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
        pot = PS::Comm::getSum(pot_loc);
        kin = PS::Comm::getSum(kin_loc);
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
            PS::F64vec rij = planet.pos_car - psys[i].pos_car;
            pot_loc -= psys[i].mass*planet.mass / sqrt(rij*rij);
        }
        pot = PS::Comm::getSum(pot_loc);
        kin = PS::Comm::getSum(kin_loc);
        tot = pot + kin;
    }
    void dump(std::ostream & fout = std::cerr){
        fout<<"tot= "<<tot<<" pot= "<<pot<<" kin= "<<kin<<std::endl;
    }
};

template<typename Tpsys>
void SetID(Tpsys & psys){
    PS::S32 n = psys.getNumberOfParticleLocal();
    PS::S32 n_proc  = PS::Comm::getNumberOfProc();
    PS::S32 my_rank = PS::Comm::getRank();
    PS::ReallocatableArray<PS::S32> n_ar;
    n_ar.resizeNoInitialize(n_proc);
    PS::Comm::allGather(&n, 1, n_ar.getPointer());
    PS::S32 offset = 0;
    for(auto i=0; i<my_rank; i++){
        offset += n_ar[i];
    }
    for(auto i=0; i<n; i++){
        psys[i].id = i + offset;
    }
}


class DiskInfo{
    PS::F64 ax_;
    PS::F64 delta_ax_;
    PS::F64 r_hill_;
    PS::F64 r_phy_;
    PS::F64 dens_;
    PS::F64 m_ptcl_;
    PS::F64 kappa_;
    PS::F64 eta_;
    PS::F64 e_refl_;
    PS::F64 t_dur_;
    PS::F64 tau_;
    PS::S64 n_glb_;
    PS::F64 getKappa(){
	return std::pow((2.0*MY_PI/t_dur_), 2.0); // k^{'};
    }
    PS::F64 getEta(){
	const auto ln_e_refl = std::log(e_refl_);
	return 4.0*MY_PI*ln_e_refl/(t_dur_*std::sqrt(MY_PI*MY_PI+ln_e_refl)); // \eta^{'}
    }
    void setParticlesOnLayer(std::vector<PS::F64vec> & pos, const PS::S64 id_x_head, const PS::S64 id_x_tail, const PS::S64 id_y_head, const PS::S64 id_y_tail, const PS::F64 dtheta, const PS::F64 dr, const PS::F64ort box,
			     const PS::F64 offset_theta = 0.0, const PS::F64 offset_r = 0.0, const PS::F64 offset_z = 0.0, const PS::F64 eps = 0.0){
	static std::mt19937 mt(PS::Comm::getRank());
	static std::uniform_real_distribution<double> dist(0.0,1.0);
	PS::F64 pos_z = offset_z;
	for(auto i=id_x_head-3; i<=id_x_tail+3; i++){
	    PS::F64 pos_x = ((PS::F64)i)*dtheta + offset_theta;
	    if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
	    for(auto j=id_y_head-3; j<=id_y_tail+3; j++){
		PS::F64 pos_y = j*dr + offset_r;
		if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
		PS::F64 eps_x = eps*(dist(mt)-0.5)*2.0;
		PS::F64 eps_y = eps*(dist(mt)-0.5)*2.0;
		PS::F64 eps_z = eps*(dist(mt)-0.5)*2.0;
		pos.push_back(PS::F64vec(pos_x+eps_x, pos_y+eps_y, pos_z+eps_z));
	    }
	}
    }
public:
    DiskInfo(const PS::F64 delta_ax, const PS::F64 e_refl, const PS::F64 t_dur, const PS::F64 tau, const PS::F64 rphy_over_rhill, const PS::S64 n_glb){
	ax_ = 1.0;
	delta_ax_ = delta_ax;
	e_refl_ = e_refl;
	t_dur_  = t_dur;
	tau_    = tau;
	n_glb_ = n_glb;
	PS::F64 ax_in  = ax_ - 0.5*delta_ax_;
	PS::F64 ax_out = ax_ + 0.5*delta_ax_;
	r_phy_  = sqrt(tau_*(ax_out*ax_out - ax_in*ax_in) / n_glb_);
	r_hill_ = r_phy_ / rphy_over_rhill;
	m_ptcl_ = (r_hill_/((ax_out+ax_in)*0.5))*(r_hill_/((ax_out+ax_in)*0.5))*(r_hill_/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5;
	dens_ = m_ptcl_ * n_glb / (MY_PI*(ax_out*ax_out - ax_in*ax_in));
	kappa_ = getKappa();
	eta_ = getEta();
    }
    
    template<class Tpsys>
    void setParticles(Tpsys & psys, const PS::F64ort box, const bool layer=true, const bool random_shift=true){
	const auto ax_in  = ax_ - 0.5*delta_ax_;
	const auto ax_out = ax_ + 0.5*delta_ax_;
	const auto area = MY_PI*(ax_out*ax_out-ax_in*ax_in);
	PS::F64	dS = area / n_glb_;
	if(layer){
	    dS *= 3.0;
	}
	const auto dl = sqrt(dS);
	assert(delta_ax_ > dl);
	const auto dz = dl;
	const PS::S64 n_theta = (2.0*MY_PI*ax_ / dl);
	const auto dtheta = 2.0*MY_PI / n_theta;
	const PS::S64 n_r = (delta_ax_ / dl);
	const auto dr = delta_ax_ / n_r;

	const PS::S64 id_x_head = box.low_.x / dtheta;
	const PS::S64 id_x_tail = box.high_.x / dtheta;
	const PS::S64 id_y_head = box.low_.y / dr;
	const PS::S64 id_y_tail = box.high_.y / dr;

	PS::F64 eps = 0.0;
	if(random_shift){
	    eps = (dl > 2.0*r_phy_) ? (0.5*(dl-2.0*r_phy_))*0.9 : 0.0;
	}
    
	PS::Comm::barrier();
	if(PS::Comm::getRank()==0){
	    std::cerr<<"delta_ax_= "<<delta_ax_
		     <<" dS = "<<dS
		     <<" dl = "<<dl
		     <<" n_r= "<<n_r
		     <<" dr= "<<dr
		     <<" eps= "<<eps
		     <<std::endl;
	    std::cerr<<"n_theta= "<<n_theta
		     <<" dtheta= "<<dtheta
		     <<std::endl;
	}
	PS::Comm::barrier();
	std::vector<PS::F64vec> pos;
	PS::F64 offset_theta = dtheta*(sqrt(2.0)-1.0)*0.5;
	PS::F64 offset_r = 0.0;
	PS::F64 offset_z = 0.0;
	setParticlesOnLayer(pos, id_x_head, id_x_tail, id_y_head, id_y_tail, dtheta, dr, box, offset_theta, offset_r, offset_z, eps);
	if(layer == true){
	    offset_theta = 0.5*dtheta + dtheta*(sqrt(2.0)-1.0)*0.5;
	    offset_r = 0.5*dr;
	    offset_z = -0.5*dz;
	    setParticlesOnLayer(pos, id_x_head, id_x_tail, id_y_head, id_y_tail, dtheta, dr, box, offset_theta, offset_r, offset_z, eps);
	    offset_z = 0.5*dz;
	    setParticlesOnLayer(pos, id_x_head, id_x_tail, id_y_head, id_y_tail, dtheta, dr, box, offset_theta, offset_r, offset_z, eps);

	}
	
	PS::S64 n_loc = pos.size();
	psys.setNumberOfParticleLocal(n_loc);
	for(PS::S64 i=0; i<n_loc; i++){
	    psys[i].pos_cyl = pos[i];
	    psys[i].pos_car = ConvertCyl2Car(pos[i]);
	    psys[i].vel     = GetVel(psys[i].pos_car);
	    psys[i].mass    = m_ptcl_;
	    psys[i].r_coll  = r_phy_;
	    psys[i].r_search = 6.0*r_hill_;
	    psys[i].eps = 0.0;
	    psys[i].kappa = kappa_;
	    psys[i].eta = eta_;
	}
    }
};


int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    
    PS::Initialize(argc, argv);
    const auto my_rank = PS::Comm::getRank();
    const auto n_proc  = PS::Comm::getNumberOfProc();
    PS::S64 n_loop_merge = 1;

    PLANET.mass = 1.0;
    PLANET.pos_car = PLANET.vel = 0.0;
    
    PS::F64 time_sys = 0.0;
    PS::F64 dt = 1.0 / 256.0;
    PS::S32 n_smp = 100;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F64 theta = 0.5;
    //PS::F64 ax_cen = 1.0;
    PS::F64 delta_ax = 1e-3;
    PS::S64 n_glb = 1024;
    PS::S64 n_loc = 0;
    PS::F64 time_end = 1.0;
    PS::F64 e_refl = 0.5;    
    PS::F64 t_dur = 1.0 / 4.0;
    PS::F64 tau = 1.0;
    PS::F64 rphy_over_rhill = 1.0;

    auto ex_let_mode = PS::EXCHANGE_LET_A2A;
    bool flag_mtm = false; // mid point tree method
    bool flag_rot = false; // on the rotating reference frame

    PS::S32 id_snp = 0;
    PS::F64 dt_snp = 1.0;
    
    PS::S32 c;
    while((c=getopt(argc,argv,"t:a:d:T:l:L:n:N:e:E:D:O:R:mrh")) != -1){
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
	    n_glb = atol(optarg) * PS::Comm::getNumberOfProc();
            break;
	case 'e':
	    if(atoi(optarg) == 0){
		ex_let_mode = PS::EXCHANGE_LET_A2A;
	    }
	    else if(atoi(optarg) == 1){
		ex_let_mode = PS::EXCHANGE_LET_P2P_EXACT;
	    }
	    else if(atoi(optarg) == 2){
		ex_let_mode = PS::EXCHANGE_LET_P2P_FAST;
	    }
	    break;
        case 'E':
	    e_refl = atof(optarg);
            break;
        case 'D':
	    t_dur = 1.0 / atoi(optarg);
            break;
        case 'O':
	    tau = atof(optarg);
            break;
        case 'R':
	    rphy_over_rhill = atof(optarg);
            break;
        case 'm':
	    flag_mtm = true;
            break;
        case 'r':
	    flag_rot = true;
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

    PS::S32 nx, ny, nz;
    DivideNProc(nx, ny, nz, n_proc, delta_ax);
    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.setDomain(nx, ny, nz);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setPosRootDomainX(-MY_PI, MY_PI);
    PS::ParticleSystem<FP_t> system;
    system.initialize();
    system.setAverageTargetNumberOfSampleParticlePerProcess(n_smp);
    const auto pos_domain = GetPosDomainCyl(delta_ax, nx, ny);
    if(my_rank==0){
	std::cerr<<"my_rank= "<<my_rank<<" pos_domain= "<<pos_domain<<std::endl;
    }
    DiskInfo disk_info(delta_ax, e_refl, t_dur, tau, rphy_over_rhill, n_glb);
    disk_info.setParticles(system, pos_domain);
    
    n_loc = system.getNumberOfParticleLocal();
    n_glb = system.getNumberOfParticleGlobal();
    if(my_rank==0){
	std::cerr<<"my_rank= "<<my_rank<<" n_loc= "<<n_loc<<" n_glb= "<<n_glb<<std::endl;
    }
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

    PS::F64 mass_loc = 0.0;
    for(auto i=0; i<system.getNumberOfParticleLocal(); i++){
      mass_loc += system[i].mass;
    }
    PS::F64 mass_glb = PS::Comm::getSum(mass_loc);
    if(my_rank==0){
      std::cerr<<"mass_glb= "<<mass_glb<<std::endl;
    }
    
    Tree_t tree;
    tree.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    tree.setExchagneLETMode(ex_let_mode);
    tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);

#if defined(FORCE_CHECK)
    PS::ReallocatableArray<Force_t> force_mono(n_loc, n_loc, PS::MemoryAllocMode::Default);
    MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSpMono<EPI_t, SPJ_t, Force_t>(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
    for(auto i=0; i<n_loc; i++){
        force_mono[i].acc = system[i].acc;
        force_mono[i].pot = system[i].pot;
    }

    /*
    PS::ReallocatableArray<Force_t> force_quad(n_loc, n_loc, PS::MemoryAllocMode::Default);
    MyCalcForceAllAndWriteBack(tree, CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSpQuad<EPI_t, SPJ_t, Force_t>(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
    for(auto i=0; i<n_loc; i++){
        force_quad[i].acc = system[i].acc;
        force_quad[i].pot = system[i].pot;
    }
    */
    
    PS::ReallocatableArray<Force_t> force_dir(n_loc, n_loc, PS::MemoryAllocMode::Default);
    system.calcForceDirectParallel(CalcForceEp<EPI_t, EPJ_t, Force_t>(), force_dir.getPointer(), dinfo, true);
    
    if(PS::Comm::getRank()==0){
        for(auto i=0; i<n_loc; i++){
            std::cout<<"force_dir[i].acc= "<<force_dir[i].acc
                     <<"force_mono[i].acc= "<<force_mono[i].acc<<" ( "<<force_mono[i].acc-force_dir[i].acc
	      //<<" ) force_quad[i].acc= "<<force_quad[i].acc<<" ( "<<force_quad[i].acc-force_dir[i].acc<<" ) "
                     <<std::endl;            
        }
    }
    PS::Finalize();
    return 0;
#endif

    CalcForceFromPlanet(system, PLANET);

    Energy eng_init;
    eng_init.calc(system);

    PS::S64 n_loop = 0;
    PS::S64 n_loop_prev = 0;
    PS::F64 eng_disp = 0.0;
    while(time_sys <= time_end){

        PS::Comm::barrier();

        n_loc = system.getNumberOfParticleLocal();
        PS::F64 va0_loc = 0.0;
        for(int i=0; i<n_loc; i++){
            va0_loc += system[i].mass*system[i].acc_dash*system[i].vel;
        }

        Kick(system, 0.5*dt);
        Drift(system, dt);
        time_sys += dt;
        n_loop++;
        if(n_loop % n_loop_merge == 0){
	    if(my_rank==0){
	        std::cerr<<"n_loop= "<<n_loop<<" n_loop_merge= "<<n_loop_merge<<std::endl;
	    }
            UpdateCyl(system);
	    if(flag_rot){
		RigidRotation(system, -dt*(n_loop-n_loop_prev) );
		n_loop_prev = n_loop;
	    }
	    
	    if(flag_mtm){
		Rotate(system, dt*n_loop_merge*0.5);
	    }
            dinfo.decomposeDomainAll(system);
            system.adjustPositionIntoRootDomain(dinfo);
            system.exchangeParticle(dinfo);
            n_loc = system.getNumberOfParticleLocal();
            n_glb = system.getNumberOfParticleGlobal();
        }
        tree.clearCounterAll();
	
	if(flag_mtm){	
	    if(n_loop % n_loop_merge == 0){
		tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
		Rotate(system, -dt*n_loop_merge*0.5);
	    }
	    tree.clearCounterAll();
	    tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::REUSE_LIST);
	} else {
	    const auto int_mode = (n_loop % n_loop_merge == 0) ? PS::MAKE_LIST_FOR_REUSE : PS::REUSE_LIST;
	    tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, int_mode);
	}
	
        PS::MemoryPool::checkEmpty();
	
        CalcForceFromPlanet(system, PLANET);
    
        Kick(system, 0.5*dt);

        n_loc = system.getNumberOfParticleLocal();
        PS::F64 va1_loc = 0.0;
        for(int i=0; i<n_loc; i++){
            va1_loc += system[i].mass*system[i].acc_dash*system[i].vel;
        }
        eng_disp += (va0_loc+va1_loc)*dt*0.5;
        PS::F64 eng_disp_glb = PS::Comm::getSum(eng_disp);
        Energy eng_now;
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
        PS::S64 size_used_glb = PS::Comm::getMaxValue(size_used_loc);
        if(size_used_loc == size_used_glb){
            std::cerr<<"my_rank= "<<my_rank<<" tree.getUsedMemorySize()= "<<tree.getUsedMemorySize()<<std::endl;
            std::cerr<<" tree.dumpMemSizeUsed():"<<std::endl;
            tree.dumpMemSizeUsed(std::cerr);
        }

	/*
	if(time_sys >= time_snp){
	    FileHeader file_header;
	    char * file_base;
	    system.writeParticleBinary(file_prefix, "%s_%05d_%05d_", file_header)
	    time_snp += time_sys + dt_snp;
	    id_snp++;
	}
	*/

    }
    PS::Finalize();
    return 0;
}
