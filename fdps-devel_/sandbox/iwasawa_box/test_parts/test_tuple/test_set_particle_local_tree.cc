#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<tuple>
#include<particle_simulator.hpp>



class Force{
public:
    PS::F64vec acc;
    PS::F64    pot;
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
};
class FPGrav0{
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
    void copyFromFP(const FPGrav0 & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
    }
    void copyFromForce(const Force & force) {
        acc = force.acc;
        pot = force.pot;
    }
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
};

class FPGrav1{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;
    PS::F64    x, y, z;
    static PS::F64 eps;
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    void copyFromFP(const FPGrav1 & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
    }
    void copyFromForce(const Force & force) {
        acc = force.acc;
        pot = force.pot;
    }
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
};
class EPGrav{
public:
    PS::F64    mass;
    PS::F64vec pos;
    static PS::F64 eps;
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    void copyFromFP(const FPGrav0 & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
    }
    void copyFromFP(const FPGrav1 & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
    }
};

PS::F64 FPGrav0::eps = 1.0/32.0;
PS::F64 FPGrav1::eps = 1.0/32.0;
PS::F64 EPGrav::eps = 1.0/32.0;

template <class TParticleJ>
void CalcGravity(const EPGrav * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 Force * force) {
    PS::F64 eps2 = EPGrav::eps * EPGrav::eps;
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec rij    = xi - ep_j[j].getPos();
            PS::F64    r3_inv = rij * rij + eps2;
            PS::F64    r_inv  = 1.0/sqrt(r3_inv);
            r3_inv  = r_inv * r_inv;
            r_inv  *= ep_j[j].getCharge();
            r3_inv *= r_inv;
            ai     -= r3_inv * rij;
            poti   -= r_inv;
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}

void makeColdUniformSphere(const PS::F64 mass_glb,
                           const PS::S64 n_glb,
                           const PS::S64 n_loc,
                           PS::F64 *& mass,
                           PS::F64vec *& pos,
                           PS::F64vec *& vel,
                           const PS::F64 eng = -0.25,
                           const PS::S32 seed = 0) {
    
    assert(eng < 0.0);
    {
        PS::MTTS mt;
        mt.init_genrand(seed);
        for(PS::S32 i = 0; i < n_loc; i++){
            mass[i] = mass_glb / n_glb;
            const PS::F64 radius = 3.0;
            do {
                pos[i][0] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][1] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][2] = (2. * mt.genrand_res53() - 1.) * radius;
            }while(pos[i] * pos[i] >= radius * radius);
            vel[i][0] = 0.0;
            vel[i][1] = 0.0;
            vel[i][2] = 0.0;
        }
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
}

template<class Tpsys>
void setParticlesColdUniformSphere(Tpsys & psys,
                                   const PS::S32 n_glb,
                                   PS::S32 & n_loc,
                                   const PS::S32 seed) {

    n_loc = n_glb; 
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64    * mass = new PS::F64[n_loc];
    PS::F64vec * pos  = new PS::F64vec[n_loc];
    PS::F64vec * vel  = new PS::F64vec[n_loc];
    const PS::F64 m_tot = 1.0;
    const PS::F64 eng   = -0.25;
    makeColdUniformSphere(m_tot, n_glb, n_loc, mass, pos, vel, eng, seed);
    for(PS::S32 i = 0; i < n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos  = pos[i];
        psys[i].vel  = vel[i];
        psys[i].id   = i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

template<class Tpsys>
void kick(Tpsys & system,
          const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].pos  += system[i].vel * dt;
    }
}

template<class Tpsys>
void calcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
    if(clear){
        etot = ekin = epot = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * (system[i].pot + system[i].mass / FPGrav0::eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);    
}

void printHelp() {
    std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
    std::cerr<<"s: time_step (default: 1.0 / 128.0)"<<std::endl;
    std::cerr<<"d: dt_diag (default: 1.0 / 8.0)"<<std::endl;
    std::cerr<<"D: dt_snap (default: 1.0)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

void makeOutputDirectory(char * dir_name) {
    struct stat st;
    PS::S32 ret;
    if (PS::Comm::getRank() == 0) {
        if (stat(dir_name, &st) != 0) {
            ret = mkdir(dir_name, 0777);
        } else {
            ret = 0; // the directory named dir_name already exists.
        }
    } 
    PS::Comm::broadcast(&ret, 1);
    if (ret == 0) {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
        PS::Abort();
    }
}



int main(int argc, char *argv[]) {

    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F32 time_end = 10.0;
    PS::F32 dt = 1.0 / 128.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_tot = 1024;
    PS::S32 c;
    sprintf(dir_name,"./result");
    opterr = 0;
    while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
        switch(c){
        //case 'o':
        //    sprintf(dir_name,optarg);
        //    break;
        case 't':
            theta = atof(optarg);
            std::cerr << "theta =" << theta << std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr << "time_end = " << time_end << std::endl;
            break;
        case 's':
            dt = atof(optarg);
            std::cerr << "time_step = " << dt << std::endl;
            break;
        case 'd':
            dt_diag = atof(optarg);
            std::cerr << "dt_diag = " << dt_diag << std::endl;
            break;
        case 'D':
            dt_snap = atof(optarg);
            std::cerr << "dt_snap = " << dt_snap << std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr << "n_group_limit = " << n_group_limit << std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            std::cerr << "n_tot = " << n_tot << std::endl;
            break;
        case 'h':
            if(PS::Comm::getRank() == 0) {
                printHelp();
            }
            PS::Finalize();
            return 0;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
                printHelp();
            }
            PS::Abort();
        }
    }

    makeOutputDirectory(dir_name);

    std::ofstream fout_eng;

    if(PS::Comm::getRank() == 0) {
        char sout_de[1024];
	sprintf(sout_de, "%s/t-de.dat", dir_name);
        fout_eng.open(sout_de);
        fprintf(stdout, "This is a sample program of N-body simulation on FDPS!\n");
        fprintf(stdout, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stdout, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }

    PS::ParticleSystem<FPGrav0> system_grav0;
    PS::ParticleSystem<FPGrav1> system_grav1;
    PS::ParticleSystem<FPGrav1> system_grav2;
    system_grav0.initialize();
    system_grav1.initialize();
    system_grav2.initialize();
    PS::S32 n_loc0    = 0;
    PS::S32 n_loc1    = 0;
    PS::S32 n_loc2    = 0;
    PS::F32 time_sys = 0.0;
    if(PS::Comm::getRank() == 0) {
        setParticlesColdUniformSphere(system_grav0, n_tot, n_loc0, 0);
        setParticlesColdUniformSphere(system_grav1, n_tot, n_loc1, 1);
        setParticlesColdUniformSphere(system_grav2, n_tot, n_loc2, 2);
    } else {
        system_grav0.setNumberOfParticleLocal(n_loc0);
        system_grav1.setNumberOfParticleLocal(n_loc1);
        system_grav2.setNumberOfParticleLocal(n_loc2);
    }
    const PS::F32 wgh = 1.23;
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
#if 0
    //dinfo.collectSampleParticle(system_grav0, true, wgh);
    //dinfo.collectSampleParticle(system_grav1, false, wgh);
    //dinfo.collectSampleParticle(system_grav0, true);
    //dinfo.collectSampleParticle(system_grav1, false);
    //dinfo.decomposeDomain();
#elif 0
    dinfo.decomposeDomainAll(system_grav0, wgh);
#else
    dinfo.decomposeDomainAll(system_grav0, wgh);
    auto sys_tpl = std::tie(system_grav1, system_grav2);
    dinfo.decomposeDomainAll(sys_tpl, wgh);
#endif
    if(PS::Comm::getRank()==0){
        for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++){
            std::cerr<<"dinfo.getPosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }
    PS::TreeForForceLong<Force, EPGrav, EPGrav>::Monopole tree_grav;
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);

#if 0
    bool clear = true;
    //tree_grav.setParticleLocalTree(system_grav0, system_grav1);
    //tree_grav.calcForceAll(CalcGravity<EPGrav>,
    //                       CalcGravity<PS::SPJMonopole>,
    //                       system_grav0,
    //                       system_grav1,
    //                       dinfo,
    //                       true);
    tree_grav.calcForceAll(CalcGravity<EPGrav>,
                           CalcGravity<PS::SPJMonopole>,
                           system_grav0,
                           system_grav1,
                           dinfo);
#else
    //tree_grav.setParticleLocalTree(sys_tpl);
    //tree_grav.calcForceAll(CalcGravity<EPGrav>,
    //                       CalcGravity<PS::SPJMonopole>,
    //                       sys_tpl,
    //                       dinfo);

    tree_grav.calcForceAllAndWriteBack(CalcGravity<EPGrav>,
                                       system_grav0,
                                       dinfo);
    
    tree_grav.calcForceAllAndWriteBack(CalcGravity<EPGrav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       system_grav0,
                                       dinfo);
    
    tree_grav.calcForceAllAndWriteBack(CalcGravity<EPGrav>,
                                       sys_tpl,
                                       dinfo);
    tree_grav.calcForceAllAndWriteBack(CalcGravity<EPGrav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       sys_tpl,
                                       dinfo);
#endif

    
    /*
    tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       system_grav,
                                       dinfo);
    */
    PS::Finalize();
    return 0;
}
