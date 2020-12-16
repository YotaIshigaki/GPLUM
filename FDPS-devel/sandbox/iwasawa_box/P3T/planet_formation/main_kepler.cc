#include<iostream>
#include<fstream>
#include<iomanip>
#include<unistd.h>
#include<unordered_map>
#include<particle_simulator.hpp>
#include"class.hpp"
#include"hard.hpp"
#include"kepler.hpp"
#include"io.hpp"
#include"profile.hpp"

void CalculateBoundaryOfDomain(const PS::S32 &np,
			       const PS::F64vec pos_sample[],
			       const PS::S32 cid,
			       const PS::S32 &istart,
			       const PS::S32 &iend,
			       const PS::F64ort pos_root_domain,
			       PS::F64 & xlow,
			       PS::F64 & xhigh) {
    if(istart == 0) {
	xlow  = pos_root_domain.low_[cid];
    } 
    else {
	xlow  = 0.5 * (pos_sample[istart-1][cid] + pos_sample[istart][cid]);
    }
    if(iend == np - 1) {
	xhigh = pos_root_domain.high_[cid];
    } 
    else {
	xhigh = 0.5 * (pos_sample[iend][cid] + pos_sample[iend+1][cid]);
    }
}

template<class Tpsys>
PS::F64 GetRootFullLenght(const Tpsys & psys, const PS::F64vec & cen){
    PS::S64 nloc = psys.getNumberOfParticleLocal();
    PS::F64 len_loc_max = 0.0;
    for(PS::S32 i=0; i<nloc; i++){
	PS::F64vec dr = psys[i].pos - cen;
	for(PS::S32 k=0; k<3; k++){
	    if(len_loc_max < dr[k]) len_loc_max = dr[k];
	}
    }
    return 2.1*fabs(PS::Comm::getMaxValue(len_loc_max));
}

template<class T>
void Print(const T str, std::ostream & fout){
#ifdef DEBUG_PRINT_PLANET
    fout<<str<<std::endl;
#endif
}

class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    FileHeader(){
        n_body = 0;
        time = 0.0;
    }
    FileHeader(const PS::S64 n, const PS::F64 t){
        n_body = n;
        time = t;
    }
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lld\t%lf\n", &n_body, &time);
	std::cout<<"n_body="<<n_body<<" time="<<time<<std::endl;
        return n_body;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%lf\n", n_body, time);
    }
};

template<class Tpsys>
PS::F64 GetMassMax(const Tpsys & system, const PS::S64 n){
    PS::F64 m_max_loc = -1.0;
    for(PS::S64 i=0; i<n; i++){
        if(m_max_loc < system[i].mass) m_max_loc = system[i].mass;
    }
    return PS::Comm::getMaxValue(m_max_loc);
}

template<class Tpsys>
void CalcEng(const Tpsys & psys, Energy & eng_glb, bool clear=true){
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot += 0.5 * psys[i].mass * psys[i].pot_tot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
    }
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.tot = eng_glb.kin + eng_glb.pot;
}

template<class Tpsys>
void CalcEngKepler(const Tpsys & psys, 
		   Energy & eng_glb, 
		   PS::F64 m_sun, 
		   bool clear=true){
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot_planet += 0.5 * psys[i].mass * psys[i].pot_tot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
        eng_loc.pot -= psys[i].mass * m_sun / sqrt(psys[i].pos * psys[i].pos);
    }
    eng_loc.pot += eng_loc.pot_planet;
    eng_glb.pot_planet += PS::Comm::getSum(eng_loc.pot_planet);
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.tot = eng_glb.kin + eng_glb.pot;
}

template<class Tpsys>
void CalcEngKepler(const Tpsys & psys, 
		   const PS::F64 eng_disp_loc, 
                   Energy & eng_glb, 
		   PS::F64 m_sun, 
		   bool clear=true){
    PS::F64 eng_disp_cum_glb = eng_glb.disp;
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot_planet += 0.5 * psys[i].mass * psys[i].pot_tot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
        eng_loc.pot -= psys[i].mass * m_sun / sqrt(psys[i].pos * psys[i].pos);
    }
    eng_loc.pot += eng_loc.pot_planet;
    eng_loc.disp = eng_disp_loc;
    eng_glb.pot_planet += PS::Comm::getSum(eng_loc.pot_planet);
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.disp += PS::Comm::getSum(eng_loc.disp) + eng_disp_cum_glb;
    eng_glb.tot = eng_glb.kin + eng_glb.pot + eng_glb.disp;
}


template<class Tpsys, class Ttree>
void Kick(Tpsys & system,
          const Ttree & tree,
          const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
	system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys, class Ttree>
void Drift(Tpsys & system,
           const Ttree & tree,
           const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
        if(tree.getForce(i).n_ngb <= 0){
            system[i].pos  += system[i].vel * dt;
        }
    }
}

#if 1
template<class Tpsys, class Ttree>
void DriftKepler(Tpsys & system,
                 const Ttree & tree,
                 const PS::F64 dt,
                 const PS::F64 mass_sun = 1.0){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        if(tree.getForce(i).n_ngb <= 0){
	    PS::F64vec pos0 = 0.0;
	    PS::F64vec vel0 = 0.0;
	    const PS::F64 mass1 = 0.0;
	    PS::F64vec pos1 = system[i].pos;
	    PS::F64vec vel1 = system[i].vel;
	    DriveKepler(mass_sun, mass1, pos0, pos1, vel0, vel1, dt);
	    system[i].pos = pos1;;
	    system[i].vel = vel1;
        }
    }
}

template<class Tpsys, class Ttree>
void DriftKeplerDebug(Tpsys & system,
		      const Ttree & tree,
		      const PS::F64 dt,
		      const PS::F64 mass_sun,
		      std::ostream & fout){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        if(tree.getForce(i).n_ngb <= 0){
	    fout<<"i0= "<<i<<std::endl;
	    PS::F64vec pos0 = 0.0;
	    PS::F64vec vel0 = 0.0;
	    const PS::F64 mass1 = 0.0;
	    PS::F64vec & pos1 = system[i].pos;
	    PS::F64vec & vel1 = system[i].vel;
	    fout<<"i1= "<<i<<std::endl;
	    DriveKepler(mass_sun, mass1, pos0, pos1, vel0, vel1, dt);
	    fout<<"i2= "<<i<<std::endl;
        }
    }
}
#else

template<class Tpsys, class Ttree>
void DriftKepler(Tpsys & system,
                 const Ttree & tree,
                 const PS::F64 dt,
                 const PS::F64 mass_sun = 1.0){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
	static __thread PS::F64vec pos0, vel0, pos1, vel1;
	static __thread PS::F64 mass1;
        if(tree.getForce(i).n_ngb <= 0){
	    pos0 = 0.0;
	    vel0 = 0.0;
	    mass1 = 0.0;
	    pos1 = system[i].pos;
	    vel1 = system[i].vel;
	    DriveKepler(mass_sun, mass1, pos0, pos1, vel0, vel1, dt);
	    system[i].pos = pos1;
	    system[i].vel = vel1;
        }
    }
}

#endif

template<class Tpsys>
void setRoutRinRhillVdisp(const Tpsys & system_soft,
                          const PS::F64vec pos_sun,
                          const PS::F64vec vel_sun,
                          const PS::F64 mass_sun,
                          const PS::F64 ratio_r_cut,
                          const PS::F64 ratio_r_search,
                          PS::F64 & r_out, 
                          PS::F64 & r_in, 
                          PS::F64 & r_hill_max_glb,
                          PS::F64 & vel_disp,
			  PS::F64 & ecc_rms,
			  PS::F64 & inc_rms,
                          PS::F64 & mass_planet_max_glb,
                          PS::F64 & mass_planet_tot_glb){
    r_out = r_in = r_hill_max_glb = vel_disp = mass_planet_max_glb = mass_planet_tot_glb = 0.0;
    const PS::S32 n_loc = system_soft.getNumberOfParticleLocal();
    const PS::S32 n_glb = system_soft.getNumberOfParticleGlobal();
    PS::F64 mass_planet_max_loc = 0.0;
    PS::F64 mass_planet_tot_loc = 0.0;
    PS::F64 r_hill_max_loc = 0.0;
    PS::F64 v_kep_max_loc = 0.0;
    PS::F64 ecc_sq_tot_loc = 0.0;
    PS::F64 inc_sq_tot_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 ax, ecc, inc;
        PosVel2AxEccInc(ax, ecc, inc,
                        pos_sun, system_soft[i].pos,
                        vel_sun, system_soft[i].vel,
                        mass_sun,        PS::F64(0.0));
        ecc_sq_tot_loc += ecc * ecc;
        inc_sq_tot_loc += inc * inc;
        PS::F64 v_kep = mass_sun / ax;
        PS::F64 r_hill = ax*ax*ax*(2.0*system_soft[i].mass) / (3.0*mass_sun);
        if(r_hill_max_loc < r_hill) r_hill_max_loc = r_hill;
        if(v_kep_max_loc < v_kep) v_kep_max_loc = v_kep;
        if(mass_planet_max_loc < system_soft[i].mass) mass_planet_max_loc = system_soft[i].mass;
        mass_planet_tot_loc += system_soft[i].mass;
    }
    mass_planet_max_glb = PS::Comm::getMaxValue(mass_planet_max_loc);
    mass_planet_tot_glb = PS::Comm::getSum(mass_planet_tot_loc);
    r_hill_max_loc = cbrt(r_hill_max_loc);
    r_hill_max_glb = PS::Comm::getMaxValue(r_hill_max_loc);
    v_kep_max_loc = sqrt(v_kep_max_loc);
    PS::F64 v_kep_max_glb = PS::Comm::getMaxValue(v_kep_max_loc);
    ecc_rms = sqrt( PS::Comm::getSum(ecc_sq_tot_loc) / n_glb );
    inc_rms = sqrt( PS::Comm::getSum(inc_sq_tot_loc) / n_glb );
    vel_disp = v_kep_max_glb * sqrt( ecc_rms*ecc_rms + inc_rms*inc_rms);
    r_out = r_hill_max_glb * ratio_r_search;
    r_in = r_out * ratio_r_cut;
}

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

    //const PS::F64 flops_per_proc = 1.28e11; // for K computer (8thread)
    const PS::F64 flops_per_core = 2.6e9*8*2*2;
    const PS::F64 n_op_ep_ep = 30.0;
    const PS::F64 n_op_ep_sp = 60.0;
    //const PS::F64 n_op_ep_ep = 24.0*2.0;
    //const PS::F64 n_op_ep_sp = 24.0;

    //PS::F64 Tbegin = PS::GetWtime();
    PS::F64 r_merge_factor = 1.0;
    PS::F64 ratio_r_cut = 0.1;
    PS::F64 time_sys = 0.0;
    PS::F64 theta = 0.4;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 n_smp_ave = 100;
    PS::F64 time_end = 1000.0 * 8.0; // 1000[yr] * 8[T/yr]
    PS::F64 r_cut_factor = 3.0;
    //PS::F64 dt_soft = 1.0; // about 1/8 [yr]
    PS::F64 dt_soft = 1.0/16.0; // about 1/128 [yr]
    PS::F64 eta = 0.1;
    char dir_name[1024];
    PS::S64 n_glb = 16384;
    PS::F64 ax_in = 0.98; //[AU]
    PS::F64 ax_out = 1.02; //[AU]
    PS::F64 ecc_sigma_hill = 2.0;
    PS::F64 inc_sigma_hill = 1.0;
    PS::F64 dens = 10.0; //[g/cm^2]
    PS::F64 mass_sun = 1.0; //[Msun]
    PS::F64vec pos_sun = 0.0;
    PS::F64vec vel_sun = 0.0;
    PS::F64 dt_snp = 100.0;
    PS::F64 search_factor = 3.0;

    int c;
#ifdef READ_FILE
    char sinput[2048];
    while((c=getopt(argc,argv,"i:o:d:D:E:t:r:T:n:s:S:l:R:r:h")) != -1){
        switch(c){
        case 'i':
            sprintf(sinput, optarg);
            std::cerr<<"sinput="<<sinput<<std::endl;
            break;
        case 'o':
            sprintf(dir_name, optarg);
            std::cerr<<"dir_name="<<dir_name<<std::endl;
            break;
        case 'd':
            dt_soft = 1.0 / atof(optarg);
            std::cerr<<"dt_soft="<<dt_soft<<std::endl;
            break;
        case 'D':
            dt_snp = atof(optarg);
            std::cerr<<"dt_snp="<<dt_snp<<std::endl;
            break;
        case 'E':
            eta = atof(optarg);
            std::cerr<<"eta="<<eta<<std::endl;
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
        case 's':
            n_smp_ave = atoi(optarg);
            std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'S':
            search_factor = atoi(optarg);
            std::cerr<<"search_factor="<<search_factor<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'r':
            r_merge_factor = atof(optarg);
            std::cerr<<"r_merge_factor="<<r_merge_factor<<std::endl;
            break;
        case 'R':
            r_cut_factor = atof(optarg);
            std::cerr<<"r_cut_factor="<<r_cut_factor<<std::endl;
            break;
        case 'h':
            std::cerr<<"i: input_file"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"d: inv_dt (dafult 16 ~ 1/128yr  )"<<std::endl;
            std::cerr<<"D: dt_snp (dafult 100 ~ 12yr  )"<<std::endl;
            std::cerr<<"E: eta (dafult 0.1)"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 8000 ~ 1yr)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"N: n_glb (dafult: 16384)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 100)"<<std::endl;
            std::cerr<<"S: search_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            std::cerr<<"r: r_merge_factor (dafult: 1.0)"<<std::endl;
            std::cerr<<"R: r_cut_factor (dafult: 3.0)"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }
#else //READ_FILE
    PS::S32 seed = 0;
    while((c=getopt(argc,argv,"a:A:i:e:o:d:D:E:t:r:T:n:N:s:S:l:R:r:x:h")) != -1){
        switch(c){
        case 'a':
            ax_in = atof(optarg);
	    std::cerr<<"ax_in="<<ax_in<<std::endl;
            break;
        case 'A':
            ax_out = atof(optarg);
            std::cerr<<"ax_out="<<ax_out<<std::endl;
            break;
        case 'i':
            inc_sigma_hill = atof(optarg);
            std::cerr<<"inc_sigma_hill="<<inc_sigma_hill<<std::endl;
            break;
        case 'e':
            ecc_sigma_hill = atof(optarg);
            std::cerr<<"ecc_sigma_hill="<<ecc_sigma_hill<<std::endl;
            break;
        case 'o':
            sprintf(dir_name,optarg);
            std::cerr<<"dir_name="<<dir_name<<std::endl;
            break;
        case 'd':
            dt_soft = 1.0 / atof(optarg);
            std::cerr<<"dt_soft="<<dt_soft<<std::endl;
            break;
        case 'D':
            dt_snp = atof(optarg);
            std::cerr<<"dt_snp="<<dt_snp<<std::endl;
            break;
        case 'E':
            eta = atof(optarg);
            std::cerr<<"eta="<<eta<<std::endl;
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
            n_glb = atol(optarg);
            std::cerr<<"n_glb="<<n_glb<<std::endl;
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'S':
            search_factor = atoi(optarg);
            std::cerr<<"search_factor="<<search_factor<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'r':
            r_merge_factor = atof(optarg);
            std::cerr<<"r_merge_factor="<<r_merge_factor<<std::endl;
            break;
        case 'R':
            r_cut_factor = atof(optarg);
            std::cerr<<"r_cut_factor="<<r_cut_factor<<std::endl;
            break;
        case 'x':
            seed = atoi(optarg);
            std::cerr<<"seed="<<seed<<std::endl;
            break;
        case 'h':
            std::cerr<<"a: ax_in [AU] (dafult 0.98AU)"<<std::endl;
            std::cerr<<"A: ax_out [AU] (dafult 1.02AU)"<<std::endl;
            std::cerr<<"i: inc_sigma_hill (dafult 1.0)"<<std::endl;
            std::cerr<<"e: ecc_sigma_hill (dafult 2.0)"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"d: inv_dt (dafult 16 ~ 1/128yr  )"<<std::endl;
            std::cerr<<"D: dt_snp (dafult 100 ~ 12yr  )"<<std::endl;
            std::cerr<<"E: eta (dafult 0.1)"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 8000 ~ 1yr)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"N: n_glb (dafult: 16384)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 100)"<<std::endl;
            std::cerr<<"S: search_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            std::cerr<<"r: r_merge_factor (dafult: 1.0)"<<std::endl;
            std::cerr<<"R: r_cut_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"x: seed (dafult: 0)"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }
#endif//READ_FILE

    PTCLHard::r_factor = r_merge_factor;
    std::ofstream fout_debug;
    std::ofstream fout_domain;
#ifdef DEBUG_PRINT_PLANET
    char sout_debug[1024];
    sprintf(sout_debug, "%s/debug_%05d.dat", dir_name, PS::Comm::getRank());
    fout_debug.open(sout_debug);
#endif
    std::ofstream fout_tcal;
    char sout_tcal[1024];
    sprintf(sout_tcal, "%s/tcal_%05d.dat", dir_name, PS::Comm::getRank());
    fout_tcal.open(sout_tcal);

    std::ofstream fout_merge;
    std::ofstream fout_diag;
    //std::ofstream fout_eng;
    std::ofstream fout_log;
#if 0
    if(PS::Comm::getRank() == 0){
        char sout_merge[1024];
        sprintf(sout_merge, "%s/merge.dat", dir_name);
        fout_merge.open(sout_merge);
    }
#else
    char sout_merge[1024];
    sprintf(sout_merge, "%s/merge_%05d.dat", dir_name, PS::Comm::getRank());
    fout_merge.open(sout_merge);
    fout_merge<<std::setprecision(15);

    char sout_debug[1024];
    sprintf(sout_debug, "%s/debug_%05d.dat", dir_name, PS::Comm::getRank());
    fout_debug.open(sout_debug);
    fout_debug<<std::setprecision(15);

    char sout_domain[1024];
    sprintf(sout_domain, "%s/domain_%05d.dat", dir_name, PS::Comm::getRank());
    fout_domain.open(sout_domain);
    fout_domain<<std::setprecision(15);

    if(PS::Comm::getRank() == 0){
        char sout_diag[1024];
        sprintf(sout_diag, "%s/diag.dat", dir_name);
        fout_diag.open(sout_diag);
        fout_diag<<std::setprecision(15);
        char sout_log[1024];
        sprintf(sout_log, "%s/log.dat", dir_name);
        fout_log.open(sout_log);
        fout_log<<std::setprecision(15);
    }
#endif
    
    PS::ParticleSystem<FPSoft> system_soft;
    system_soft.initialize();
    system_soft.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    PS::S32 n_loc;
#ifdef READ_FILE
    FileHeader file_header_read;
    std::cerr<<"reading file"<<std::endl;
    system_soft.readParticleAscii(sinput, file_header_read);
    std::cerr<<"finish reading file"<<std::endl;
    time_sys = file_header_read.time;
    std::cout<<"time_sys="<<time_sys<<std::endl;
    //std::cerr<<"time_sys="<<time_sys<<std::endl;
    n_loc = system_soft.getNumberOfParticleLocal();
    std::cerr<<"n_loc="<<n_loc<<std::endl;
#else //READ_FILE
    std::cerr<<"SetParticleKeplerDisk Beg"<<std::endl;
    SetParticleKeplerDisk(system_soft, n_glb, n_loc, time_sys,
			  ax_in, ax_out, ecc_sigma_hill, inc_sigma_hill, dens, mass_sun, seed);
    std::cerr<<"SetParticleKeplerDisk Fin"<<std::endl;
#endif
#ifdef REV_VEL
    for(PS::S32 i=0; i<n_loc; i++){
	system_soft[i].vel *= -1.0;
    }
#endif
    PS::Comm::barrier();

    PS::F64 r_out, r_in, vel_disp, r_hill_max_glb, mass_planet_max_glb, mass_planet_tot_glb;
    PS::F64 ecc_rms, inc_rms;
    setRoutRinRhillVdisp(system_soft,           pos_sun,
                         vel_sun,               mass_sun,
                         ratio_r_cut,           r_cut_factor,
                         r_out,                 r_in, 
                         r_hill_max_glb,        vel_disp,
			 ecc_rms,               inc_rms,
                         mass_planet_max_glb,   mass_planet_tot_glb);

#ifdef MERGE
    EPISoft::eps = 0.0; // eps should be zero between the sun and planets, otherwise kepler solver doesn't work well
    //EPISoft::eps = 4.0*vel_disp*dt_soft;
#else //MERGE
    static const PS::F64 rho_ave = 3.0 * ( (1.49597871e13*1.49597871e13*1.49597871e13) / 1.989e33); // [Msun/AU^3]
    static const PS::F64 PI = 4.0*atan(1.0);
    static const PS::F64 C = 3.0/(4.0*PI*rho_ave);
    EPISoft::eps = cbrt(C*mass_planet_tot_glb/system_soft.getNumberOfParticleGlobal()) * r_merge_factor;
    //EPISoft::eps = cbrt(C*system_soft[0].mass) * r_merge_factor;
#endif
    std::cerr<<"EPISoft::eps="<<EPISoft::eps<<std::endl;

    //EPJSoft::r_search = r_out + 3.0*vel_disp*dt_soft;
    EPJSoft::r_search = r_out + search_factor*vel_disp*dt_soft;
    if(r_out <= EPISoft::eps){
        EPJSoft::r_search = 0.0;
    }

    PS::F64 unit_m = 1.989e30; //[kg]
    PS::F64 unit_v = 29.7886203575; //[km/sec]
    PS::F64 unit_t = 0.15924595715; //[yr]

    if(PS::Comm::getRank() == 0){
#ifdef READ_FILE
#else //READ_FILE
        fout_log<<"ax_in= "<<ax_in<<std::endl;
        fout_log<<"ax_out= "<<ax_out<<std::endl;
        fout_log<<"inc_sigma_hill= "<<inc_sigma_hill<<std::endl;
        fout_log<<"ecc_sigma_hill= "<<ecc_sigma_hill<<std::endl;
#endif //READ_FILE
	fout_log<<"r_out= "<<r_out<<std::endl;
	fout_log<<"r_in= "<<r_in<<std::endl;
	fout_log<<"vel_disp= "<<vel_disp<<std::endl;;
	fout_log<<"r_hill_max_glb= "<<r_hill_max_glb<<std::endl;
	fout_log<<"mass_planet_max_glb= "<<mass_planet_max_glb<<std::endl;
	fout_log<<"mass_planet_tot_glb= "<<mass_planet_tot_glb<<std::endl;
	fout_log<<"ecc_rms= "<<ecc_rms<<std::endl;
	fout_log<<"inc_rms= "<<inc_rms<<std::endl;
        fout_log<<"dir_name= "<<dir_name<<std::endl;
        fout_log<<"dt_soft= "<<dt_soft<<std::endl;
        fout_log<<"dt_snp= "<<dt_snp<<std::endl;
        fout_log<<"eta= "<<eta<<std::endl;
        fout_log<<"theta= "<<theta<<std::endl;
        fout_log<<"time_end= "<<time_end<<std::endl;
        fout_log<<"n_group_limit= "<<n_group_limit<<std::endl;
        fout_log<<"n_glb= "<<n_glb<<std::endl;
        fout_log<<"n_smp_ave= "<<n_smp_ave<<std::endl;
        fout_log<<"n_leaf_limit= "<<n_leaf_limit<<std::endl;
        fout_log<<"r_cut_factor= "<<r_cut_factor<<std::endl;
        fout_log<<"PTCLHard::r_factor= "<<PTCLHard::r_factor<<std::endl;
        fout_log<<"EPISoft::eps= "<<EPISoft::eps<<std::endl;
        fout_log<<"EPJSoft::r_search= "<<EPJSoft::r_search<<std::endl;
        fout_log<<"unit_m= "<<unit_m<<"[kg]"<<std::endl;
        fout_log<<"unit_v= "<<unit_v<<"[km/sec]"<<std::endl;
        fout_log<<"unit_t= "<<unit_t<<"[yr]"<<std::endl;
    }    


    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
    if(PS::Comm::getRank() == 0){
        fout_log<<"nx,ny,nz= "<<nx<<" "<<ny<<" "<<nz<<std::endl;
    }
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
    dinfo.collectSampleParticle(system_soft, true);
    dinfo.decomposeDomain();

#ifdef DIV_FIX
    if(n_proc%4 == 0){
	int n_loc_tmp = system_soft.getNumberOfParticleLocal();
	int * n_recv_tmp = new int[n_proc];
	int * n_recv_disp_tmp = new int[n_proc+1];
	PS::F64vec * pos_loc_tmp = new PS::F64vec[n_loc_tmp];
	for(int i=0; i<n_loc_tmp; i++){
	    pos_loc_tmp[i].z = system_soft[i].pos.z;
	    if(system_soft[i].pos.x * system_soft[i].pos.y > 0.0){
		pos_loc_tmp[i].x = fabs(system_soft[i].pos.x);
		pos_loc_tmp[i].y = fabs(system_soft[i].pos.y);
	    }
	    else{
		pos_loc_tmp[i].x = fabs(system_soft[i].pos.y);
		pos_loc_tmp[i].y = fabs(system_soft[i].pos.x);
	    }
	}
	PS::Comm::allGather(&n_loc_tmp, 1, n_recv_tmp);
	n_recv_disp_tmp[0] = 0;
	for(int i=0; i<n_proc; i++){
	    n_recv_disp_tmp[i+1] = n_recv_disp_tmp[i] + n_recv_tmp[i];
	}
	int n_glb_tmp = n_recv_disp_tmp[n_proc];
	PS::F64vec * pos_glb_tmp = new PS::F64vec[n_glb_tmp];
	PS::Comm::allGatherV(pos_loc_tmp, n_loc_tmp, pos_glb_tmp, n_recv_tmp, n_recv_disp_tmp);
	PS::S32 nx_half = nx / 2;
	PS::S32 ny_half = ny / 2;
	PS::S32 n_proc_quat = n_proc / 4;
	PS::F64ort * pos_domain_tmp = new PS::F64ort[n_proc];
	if(PS::Comm::getRank() == 0){
	    PS::S32 * istart = new PS::S32[n_proc_quat];
	    PS::S32 * iend   = new PS::S32[n_proc_quat];
	    PS::F64ort pos_root_domain_tmp = PS::F64ort( PS::F64vec(0.0, 0.0, -PS::LARGE_FLOAT), PS::F64vec(PS::LARGE_FLOAT, PS::LARGE_FLOAT, PS::LARGE_FLOAT));
	    std::sort(pos_glb_tmp, pos_glb_tmp+n_glb_tmp, 
		      [](const PS::F64vec & l, const PS::F64vec & r)->bool{return l.x < r.x;}
		      );
	    for(PS::S32 i = 0; i < n_proc_quat; i++) {
		istart[i] = ((PS::S64)(i) * (PS::S64)(n_glb_tmp)) / (PS::S64)(n_proc_quat);
		if(i > 0) iend[i-1] = istart[i] - 1;
	    }
	    iend[n_proc_quat-1] = n_glb_tmp - 1;
	    for(PS::S32 ix = 0; ix<nx_half; ix++) {
		PS::S32 ix0 =  ix      * ny_half * nz;
		PS::S32 ix1 = (ix + 1) * ny_half * nz;
		PS::F64 x0, x1;
		CalculateBoundaryOfDomain(n_glb_tmp, pos_glb_tmp, 0, istart[ix0], iend[ix1-1], pos_root_domain_tmp, x0, x1);
		for(PS::S32 i=0; i<ny_half*nz; i++) {
		    PS::S32 offset = (nx_half+ix)*ny;
		    pos_domain_tmp[offset+i].low_[0]  = x0;
		    pos_domain_tmp[offset+i].high_[0] = x1;
		    pos_domain_tmp[offset+ny_half+i].low_[0]  = x0;
		    pos_domain_tmp[offset+ny_half+i].high_[0] = x1;

		    offset = (nx_half-1-ix)*ny;
		    pos_domain_tmp[offset+i].low_[0]  = -x1;
		    pos_domain_tmp[offset+i].high_[0] = -x0;
		    pos_domain_tmp[offset+ny_half+i].low_[0]  = -x1;
		    pos_domain_tmp[offset+ny_half+i].high_[0] = -x0;
		}
	    }

	    for(PS::S32 ix = 0; ix<nx_half; ix++) {
		PS::S32 ix0 =  ix      * ny_half * nz;
		PS::S32 ix1 = (ix + 1) * ny_half * nz;
		std::sort(pos_glb_tmp+istart[ix0], pos_glb_tmp+iend[ix1-1]+1,
			  [](const PS::F64vec & l, const PS::F64vec & r)->bool{return l.y < r.y;}
			  );
		PS::S32 n_tmp_y = iend[ix1-1] - istart[ix0] + 1;
		for(PS::S32 iy = 0; iy<ny_half; iy++) {
		    PS::S32 iy0 = ix0 +  iy      * nz;
		    PS::S32 iy1 = ix0 + (iy + 1) * nz;
		    PS::F64 y0, y1;
		    CalculateBoundaryOfDomain(n_tmp_y, pos_glb_tmp+istart[ix0], 1, istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], pos_root_domain_tmp, y0, y1);
		    for(PS::S32 i=0; i<nz; i++) {
			PS::S32 offset = (nx_half+ix)*ny + ny_half + iy;
			pos_domain_tmp[offset+i].low_[1]  = y0;
			pos_domain_tmp[offset+i].high_[1] = y1;
			offset = (nx_half-ix-1)*ny + ny_half + iy;
			pos_domain_tmp[offset+i].low_[1]  = y0;
			pos_domain_tmp[offset+i].high_[1] = y1;
			offset = (nx_half+ix)*ny + ny_half-iy-1;
			pos_domain_tmp[offset+i].low_[1]  = -y1;
			pos_domain_tmp[offset+i].high_[1] = -y0;
			offset = (nx_half-ix-1)*ny + ny_half-iy-1;
			pos_domain_tmp[offset+i].low_[1]  = -y1;
			pos_domain_tmp[offset+i].high_[1] = -y0;
		    }
		}
	    }
	    delete [] istart;
	    delete [] iend;
	}
	PS::Comm::broadcast(pos_domain_tmp, n_proc);
	for(PS::S32 i=0; i<n_proc; i++){
	    pos_domain_tmp[i].low_.z  = -PS::LARGE_FLOAT;
	    pos_domain_tmp[i].high_.z =  PS::LARGE_FLOAT;
	    dinfo.setPosDomain(i, pos_domain_tmp[i]);
	}
	delete [] n_recv_tmp;
	delete [] n_recv_disp_tmp;
	delete [] pos_loc_tmp;
	delete [] pos_glb_tmp;
	delete [] pos_domain_tmp;
    }
    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0){
	for(PS::S32 i=0; i<n_proc; i++){
	    std::cout<<std::setprecision(15)<<i<<"   "<<dinfo.getPosDomain(i)<<std::endl;
	}
    }
#endif // DIV_FIX

    system_soft.exchangeParticle(dinfo); 
    n_loc = system_soft.getNumberOfParticleLocal();

#ifdef USE_QUAD
    PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithScatterSearch tree_soft;
#else
    PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithScatterSearch tree_soft;
#endif

    tree_soft.initialize(n_glb, theta, n_leaf_limit, n_group_limit);

    const PS::F64vec root_cen(0.0);
    PS::F64 root_len = GetRootFullLenght(system_soft, root_cen);
    tree_soft.setParticleLocalTree(system_soft);
    tree_soft.setRootCell(root_len, root_cen);
    tree_soft.mortonSortLocalTreeOnly();
    tree_soft.linkCellLocalTreeOnly();
    tree_soft.calcMomentLocalTreeOnly();
    tree_soft.exchangeLocalEssentialTree(dinfo);
    tree_soft.setLocalEssentialTreeToGlobalTree();
    tree_soft.mortonSortGlobalTreeOnly();
    tree_soft.linkCellGlobalTreeOnly();
    tree_soft.calcMomentGlobalTreeOnly();
    tree_soft.makeIPGroup();
#ifdef USE_QUAD
    tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSPQuad(), true);
    //tree_soft.calcForceAllAndWriteBack(CalcForceEPEP(), CalcForceEPSPQuad(), system_soft, dinfo);
#else
    tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSP(), true);
    //tree_soft.calcForceAllAndWriteBack(CalcForceEPEP(), CalcForceEPSP(), system_soft, dinfo);
#endif
    n_loc = system_soft.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
	system_soft[i].copyFromForce(tree_soft.getForce(i));
    }
    /*
    for(int i=0; i<system_soft.getNumberOfParticleLocal(); i++){
	system_soft[i].acc_direct  = system_soft[i].acc;
	system_soft[i].pot_tot_direct = system_soft[i].pot_tot;
	system_soft[i].n_ep_direct  = system_soft[i].n_ep;
	system_soft[i].n_sp_direct  = system_soft[i].n_sp;
    }
    tree_soft.calcForceDirectAndWriteBack(CalcForceEPEP(), system_soft, dinfo);
    for(int i=0; i<system_soft.getNumberOfParticleLocal(); i++){
	std::swap(system_soft[i].acc_direct,  system_soft[i].acc);
	std::swap(system_soft[i].pot_tot_direct,  system_soft[i].pot_tot);
	std::swap(system_soft[i].n_ep_direct, system_soft[i].n_ep);
	std::swap(system_soft[i].n_sp_direct, system_soft[i].n_sp);
    }
    */

    HardSystem system_hard;
    system_hard.setSun(mass_sun, pos_sun, vel_sun);

    //const PS::F64 dt_limit_hard = dt_soft * 0.25;
    //const PS::F64 dt_limit_hard = dt_soft * 0.125;
    const PS::F64 dt_limit_hard = dt_soft * 0.0625;
    const PS::F64 eta_s = eta * 0.01;
    ///////////////
    // hard part

    for(int i=0; i<system_soft.getNumberOfParticleLocal(); i++){
	system_soft[i].acc_pla = system_soft[i].acc;
    }
    system_hard.setUp(tree_soft, system_soft, system_soft.getNumberOfParticleLocal(), r_out, r_in, time_sys);

    Energy eng_init;
    system_hard.clearEngDisp();
    CalcEngKepler(system_soft, system_hard.eng_disp_, eng_init, PS::F64(mass_sun), true);
    eng_init.dump(std::cerr);
    Energy eng = eng_init;

    PS::S64 snp_id = 0; 
    PS::S64 n_loop = 0;
    PS::S64 n_loop_offset = 0;
    Wtime::cum_offset = PS::GetWtime();
    Wtime::interval_offset = PS::GetWtime();
    Wtime::clear();
    system_hard.clear_counter();
    while(time_sys < time_end){
	Print("check a", fout_debug);
	Print(n_loop, fout_debug);

	if(1){


	    if(fmod(time_sys, dt_snp) == 0.0){
		PS::Comm::barrier();
		Wtime::take_snp_offset = PS::GetWtime();
		char file_snp[1024];
		sprintf(file_snp, "%s/snap%05d.dat", dir_name, (int)snp_id++);
		FileHeader header(system_soft.getNumberOfParticleGlobal(), time_sys);
		system_soft.writeParticleAscii(file_snp, header);
		PS::Comm::barrier();
		Wtime::take_snp += PS::GetWtime() - Wtime::take_snp_offset;
	    }

	    if(n_loop != 0 && fmod(time_sys, 10.0) == 0.0){
		PS::Comm::barrier();
		Wtime::interval = PS::GetWtime() - Wtime::interval_offset;
		Wtime::cum = PS::GetWtime() - Wtime::cum_offset;
		Wtime::dump(fout_tcal, time_sys, n_loop-n_loop_offset);
		system_hard.dump_counter(fout_tcal);
		fout_tcal<<"tree_soft.getNumberOfLETEPSend1stLocal()= "<<tree_soft.getNumberOfLETEPSend1stLocal()
			 <<" tree_soft.getNumberOfLETSPSend1stLocal()= "<<tree_soft.getNumberOfLETSPSend1stLocal()
			 <<" tree_soft.getNumberOfLETEPSend1stGlobal()= "<<tree_soft.getNumberOfLETEPSend1stGlobal()
			 <<" tree_soft.getNumberOfLETSPSend1stGlobal()= "<<tree_soft.getNumberOfLETSPSend1stGlobal()<<std::endl;
		fout_tcal<<"tree_soft.getNumberOfLETEPRecv1stLocal()= "<<tree_soft.getNumberOfLETEPRecv1stLocal()
			 <<" tree_soft.getNumberOfLETSPRecv1stLocal()= "<<tree_soft.getNumberOfLETSPRecv1stLocal()
			 <<" tree_soft.getNumberOfLETEPRecv1stGlobal()= "<<tree_soft.getNumberOfLETEPRecv1stGlobal()
			 <<" tree_soft.getNumberOfLETSPRecv1stGlobal()= "<<tree_soft.getNumberOfLETSPRecv1stGlobal()<<std::endl;
		fout_tcal<<" system_soft.getNumberOfParticleLocal()= "<<system_soft.getNumberOfParticleLocal()<<std::endl;
		fout_tcal<<"time_sys= "<<time_sys<<" dinfo.getPosDomain(PS::Comm::getRank())="<<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
		fout_tcal<<"-----------------------"<<std::endl;
		fout_tcal<<std::endl;
		fout_debug<<"time_sys= "<<time_sys<<" dinfo.getPosDomain(PS::Comm::getRank())="<<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
		system_hard.clear_counter();

		CalcEngKepler(system_soft, system_hard.eng_disp_, eng, PS::F64(mass_sun), true);
		system_hard.clearEngDisp();
		PS::S64 n_body_tmp = system_soft.getNumberOfParticleGlobal();
		if(PS::Comm::getRank() == 0){
		    fout_diag<<"time_sys= "<<time_sys
			     <<" n_body= "<<n_body_tmp
			     <<" mass_planet_max= "<<mass_planet_max_glb
			     <<" r_hill_max= "<<r_hill_max_glb
			     <<" ecc_rms= "<<ecc_rms
			     <<" inc_rms= "<<inc_rms
			     <<" (eng.tot-eng_init.tot)/eng_init.tot= "<<(eng.tot-eng_init.tot)/eng_init.tot
			     <<" (eng.tot-eng_init.tot)/eng_init.pot_planet= "<<(eng.tot-eng_init.tot)/eng_init.pot_planet
			     <<" eng.tot= "<<eng.tot
			     <<" eng.pot= "<<eng.pot
			     <<" eng.kin= "<<eng.kin
			     <<" eng.pot_planet= "<<eng.pot_planet
			     <<" eng.disp= "<<eng.disp
			     <<" eng_init.tot= "<<eng_init.tot
			     <<" eng_init.pot= "<<eng_init.pot
			     <<" eng_init.kin= "<<eng_init.kin
			     <<" eng_init.pot_planet= "<<eng_init.pot_planet
			     <<" eng_init.disp= "<<eng_init.disp<<std::endl;
		}
		n_loop_offset = n_loop;
		Wtime::clear();
		Wtime::interval_offset = PS::GetWtime();
            }
        }
	Print("check b", fout_debug);

	Wtime::soft_offset = PS::GetWtime();
        /////////////
        // 1st KICK
        Kick(system_soft, tree_soft, dt_soft*0.5);
	Print("check c", fout_debug);
        system_hard.copyVelSoftToHard(system_soft);
	Print("check d", fout_debug);
        // 1st KICK
        /////////////
	Wtime::soft += PS::GetWtime() - Wtime::soft_offset;

        //////////////
        // HARD PART
	Wtime::hard_offset = PS::GetWtime();
        //////////////
        // DRIFT 1body
	Wtime::hard_1body_offset = PS::GetWtime();
        DriftKepler(system_soft, tree_soft, dt_soft);
	Wtime::hard_1body += PS::GetWtime() - Wtime::hard_1body_offset;
	Print("check e", fout_debug);
        // DRIFT 1body
        //////////////


        ////////
        // DRIFT 2body
	Wtime::hard_2body_offset = PS::GetWtime();
        system_hard.selectIsolatedSystem();
	Print("check f", fout_debug);
        system_hard.evolveIsolatedSystem(system_soft, r_out, r_in, mass_sun, pos_sun, vel_sun,
                                         eta_s, eta,  dt_soft,  dt_limit_hard);
	Print("check g", fout_debug);
	Wtime::hard_2body += PS::GetWtime() - Wtime::hard_2body_offset;
        // DRIFT 2body
        ////////

        /////////
        // GATHER HARD PTCL
	Wtime::hard_gather_data_offset = PS::GetWtime(); 
        system_hard.gatherData();
	Wtime::hard_gather_data += PS::GetWtime() - Wtime::hard_gather_data_offset;
        // GATHER HARD PTCL
        /////////

	Print("check h", fout_debug);
        ////////
        // DRIFT multi body
	Wtime::hard_multi_offset = PS::GetWtime();
        if(PS::Comm::getRank() == 0){
            if( system_hard.ptcl_multi_glb_.size() > 0){
		system_hard.evolve_multi_sirial(dt_soft, dt_limit_hard, r_out, r_in, eta_s, eta);
            }
        }
	Wtime::hard_multi += PS::GetWtime() - Wtime::hard_multi_offset;
        // DRIFT multi body
        ////////

	Print("check i", fout_debug);
        /////////
        // SCATTER HARD PTCL
	Wtime::hard_scatter_data_offset = PS::GetWtime(); 
        system_hard.scatterData();
	Wtime::hard_scatter_data += PS::GetWtime() - Wtime::hard_scatter_data_offset;
        // SCATTER HARD PTCL
        /////////
	Print("check j", fout_debug);

        /////////
        // COPY HARD PTCL
	Wtime::hard_copy_h2s_offset = PS::GetWtime();
        system_hard.copyPtclHardToSoft(system_soft);
	Wtime::hard_copy_h2s += PS::GetWtime() - Wtime::hard_copy_h2s_offset;
        // COPY HARD PTCL
        /////////

	Print("check k", fout_debug);
	
        /////////
        // MERGE HARD PTCL
	Wtime::hard_merge_offset = PS::GetWtime();
#ifdef MERGE
        const PS::S32 n_loc_old = system_soft.getNumberOfParticleLocal();
        PS::S32 n_loc_new = 0;
	Print("check l", fout_debug);
        for(PS::S32 ip=0; ip<n_loc_old; ip++){
            if(system_soft[ip].mass > 0.0){ // not merged
                system_soft[n_loc_new] = system_soft[ip];
                n_loc_new++;
            }
        }
	Print("check m", fout_debug);
        if(system_hard.merge_history_.size() > 0){
            fout_merge<<"system_hard.merge_history_.size()="<<system_hard.merge_history_.size()<<std::endl;
            fout_merge<<"n_loc_new="<<n_loc_new<<" n_loc_old="<<n_loc_old<<std::endl;
            for(size_t ip=0; ip<system_hard.merge_history_.size(); ip++){
                system_hard.merge_history_[ip].dump(fout_merge);
            }
            fout_merge<<std::endl;
            system_hard.merge_history_.clear(); // clear history
        }
	Print("check n", fout_debug);
        PS::S32 n_glb_old = PS::Comm::getSum(n_loc_old);
        PS::S32 n_glb_new = PS::Comm::getSum(n_loc_new);
        if(n_glb_new != n_glb_old){
            mass_planet_max_glb = GetMassMax(system_soft, n_loc_new);
            if(PS::Comm::getRank() == 0){
                std::cout<<time_sys<<"   "<<n_glb_new<<std::endl;
                //fout_diag<<time_sys<<"   "<<n_glb_new<<"   "<<mass_planet_max_glb<<std::endl;
            }
        }
	Print("check o", fout_debug);
        system_soft.setNumberOfParticleLocal(n_loc_new);
#endif // MERGE
	Wtime::hard_merge += PS::GetWtime() - Wtime::hard_merge_offset;
	system_hard.accumulate_counter();
	Wtime::hard += PS::GetWtime() - Wtime::hard_offset;
        // HARD PART
        //////////////


        //////////////
        // SOFT PART
	Wtime::soft_offset = PS::GetWtime();
#ifndef DIV_FIX
	if( n_loop < 12 || n_loop % 16 == 0){
	    if(PS::Comm::getRank() == 0){
		fout_domain<<time_sys<<"   "<<tree_soft.pos_root_cell_;
		for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++){
		    fout_domain<<"   "<<dinfo.getPosDomain(i);
		}
		fout_domain<<std::endl;
	    }
	    dinfo.collectSampleParticle(system_soft);
	    dinfo.decomposeDomain();
	}
#endif
        //////////////////////
        // EVALUATE SOFT FORCE
	Print("check s", fout_debug);
        system_soft.exchangeParticle(dinfo);
	Print("check t", fout_debug);
        if(fmod(time_sys, 1000.0) == 0.0){
            setRoutRinRhillVdisp(system_soft,           pos_sun,
                                 vel_sun,               mass_sun,
                                 ratio_r_cut,           r_cut_factor,
                                 r_out,                 r_in, 
                                 r_hill_max_glb,        vel_disp,
				 ecc_rms,               inc_rms,
                                 mass_planet_max_glb,   mass_planet_tot_glb);
            EPJSoft::r_search = r_out + search_factor*vel_disp*dt_soft;
            if(r_out <= EPISoft::eps){
                EPJSoft::r_search = 0.0;
            }
        }
	Print("check u", fout_debug);
	Wtime::soft_force_offset = PS::GetWtime();
#ifdef FREE_MALLOC_TREE
	tree_soft.freeMem();
	tree_soft.reallocMem();
#endif //FREE_MALLOC_TREE

	root_len = GetRootFullLenght(system_soft, root_cen);
	tree_soft.setParticleLocalTree(system_soft);
	tree_soft.setRootCell(root_len, root_cen);
	tree_soft.mortonSortLocalTreeOnly();
	tree_soft.linkCellLocalTreeOnly();
	tree_soft.calcMomentLocalTreeOnly();
	tree_soft.exchangeLocalEssentialTree(dinfo);
	tree_soft.setLocalEssentialTreeToGlobalTree();
	tree_soft.mortonSortGlobalTreeOnly();
	tree_soft.linkCellGlobalTreeOnly();
	tree_soft.calcMomentGlobalTreeOnly();
	tree_soft.makeIPGroup();
#ifdef USE_QUAD
	tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSPQuad(), true);
	//tree_soft.calcForceAllAndWriteBack(CalcForceEPEP(), CalcForceEPSPQuad(), system_soft, dinfo);
#else
	tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSP(), true);
	//tree_soft.calcForceAllAndWriteBack(CalcForceEPEP(), CalcForceEPSP(), system_soft, dinfo);
#endif
	n_loc = system_soft.getNumberOfParticleLocal();
#pragma omp parallel for
	for(PS::S32 i=0; i<n_loc; i++){
	    system_soft[i].copyFromForce(tree_soft.getForce(i));
	}

	Wtime::soft_force += PS::GetWtime() - Wtime::soft_force_offset;
	Print("check v", fout_debug);
	for(int i=0; i<system_soft.getNumberOfParticleLocal(); i++){
	    system_soft[i].acc_pla = system_soft[i].acc;
	}
        system_hard.setUp(tree_soft, system_soft, system_soft.getNumberOfParticleLocal(), r_out, r_in, time_sys);

	Print("check w", fout_debug);
        // EVALUATE SOFT FORCE
        //////////////////////

        //////////////
        // second KICK
        Kick(system_soft, tree_soft, dt_soft*0.5);
        time_sys += dt_soft;
        // second KICK
        //////////////
	Wtime::soft += PS::GetWtime() - Wtime::soft_offset;
        // SOFT PART
        //////////////

        n_loop++;
	Print("check x", fout_debug);
    }

    std::cerr<<"system_hard.merge_history_.size()="<<system_hard.merge_history_.size()<<std::endl;
    std::cerr<<"system_soft.getNumberOfParticleGlobal()="<<system_soft.getNumberOfParticleGlobal()<<std::endl;
    for(size_t ip=0; ip<system_hard.merge_history_.size(); ip++){
        system_hard.merge_history_[ip].dump();
    }

    PS::Finalize();
    return 0;
}

