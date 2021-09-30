#include<iostream>
#include<fstream>
#include<unistd.h>
#include<random>
#include<particle_simulator.hpp>

typedef PS::F64 REAL;
typedef PS::F64vec REALVEC;

struct Force{
    PS::S32 n_ngb;
    REAL mass;
    REAL pot;
    void clear() {
        n_ngb = 0.0;
        mass = 0.0;
        pot = 0.0;
    }
};

class FP{
public:
    PS::S64 id;
    REAL mass;
    REAL mass_tmp;
    REALVEC pos;
    REALVEC vel;
    PS::S32 n_ngb;
    REAL r_search;
    REALVEC getPos() const { return pos; }
    void setPos(const REALVEC & p) { pos = p; }
    REAL getRSearch() const {
	return r_search;
    }
    void copyFromForce(const Force & force){
        n_ngb = force.n_ngb;
        mass_tmp = force.mass;
    }
};

class EP{
public:
    PS::S64 id;
    REAL mass;
    REALVEC pos;
    REAL r_search;
    PS::S64 getId() const {
        return id;
    }
    REAL getCharge() const {
        return mass;
    }
    REAL getRSearch() const {
	return r_search;
    }
    REALVEC getPos() const {
        return pos;
    }
    void setPos(const REALVEC & _pos){
	pos = _pos;
    }    
    void copyFromFP(const FP & fp){
        id = fp.id;
        pos  = fp.pos;
        mass  = fp.mass;
        r_search = fp.r_search;
    }
};

struct SearchNeighborSymmetry{
    template<class Tepi, class Tepj>
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Force * force){
	for(PS::S32 i=0; i<n_ip; i++){
	    const REALVEC xi = ep_i[i].getPos();
	    const REAL r_search_i_sq = ep_i[i].getRSearch() * ep_i[i].getRSearch();
	    for(PS::S32 j=0; j<n_jp; j++){
		const REALVEC xj = ep_j[j].getPos();
		const REAL r_search_sq = (ep_j[j].getRSearch() * ep_j[j].getRSearch() > r_search_i_sq) ? (ep_j[j].getRSearch()*ep_j[j].getRSearch()) : r_search_i_sq;
		REALVEC rij = xi - xj;
		if(rij*rij <= r_search_sq){
		    force[i].n_ngb++;
                    force[i].mass += ep_j[j].mass;
		}
	    }
	}
    }
};

void MakeUniformQubeModel(const double mass_glb,
                          const long long int n_glb,
                          const long long int n_loc,
                          double *& mass,
                          PS::F64vec *& pos,
                          PS::F64vec *& vel,
                          PS::F64 full_len = 1.0,
                          const PS::F64vec offset = PS::F64vec(0.0, 0.0, 0.0),
                          const int seed = 0){

    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];
    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        pos[i].x = full_len * mt.genrand_real2() - offset.x;
        pos[i].y = full_len * mt.genrand_real2() - offset.y;
        pos[i].z = full_len * mt.genrand_real2() - offset.z;
        vel[i].x = (mt.genrand_real2() - 0.5);
        vel[i].y = (mt.genrand_real2() - 0.5);
        vel[i].z = (mt.genrand_real2() - 0.5);
    }
}

template<class Tpsys>
void SetParticlesUniformQube(Tpsys & psys,
                             const PS::S64 n_glb,
                             PS::S32 & n_loc,  
                             REAL & t_sys,
                             const PS::F64 full_len = 1.0,
                             const PS::F64vec offset = PS::F64vec(0.0, 0.0, 0.0),
                             const bool random=true){
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc;
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    
    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;
    const PS::F64 m_tot = 1.0;
    const PS::S32 seed = 0;
    if(random){
        MakeUniformQubeModel(m_tot, n_glb, n_loc, mass, pos, vel, full_len, seed);
    }
    else{
        PS::S32 n_glb_1d = cbrt(n_glb);
        assert(n_glb == (n_glb_1d*n_glb_1d*n_glb_1d) );
        PS::F64 dx = full_len / n_glb_1d;
        mass = new double[n_loc];
        pos = new PS::F64vec[n_loc];
        vel = new PS::F64vec[n_loc];
        PS::S64 i_e = i_h + n_loc;
        PS::S32 n_cnt = 0;
        for(PS::S32 i=i_h; i<i_e; i++, n_cnt++){
            PS::S32 ix = i / (n_glb_1d*n_glb_1d);
            PS::S32 iy = (i/n_glb_1d) % n_glb_1d;
            PS::S32 iz = i%n_glb_1d;
            pos[n_cnt].x = ix*dx - offset.x;
            pos[n_cnt].y = iy*dx - offset.y;
            pos[n_cnt].z = iz*dx - offset.z;
        }
        
    }
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

void printHelp() {
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);

    REAL theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S64 n_tot = 1024;
    PS::S32 c;
    while((c=getopt(argc,argv,"l:n:N:h")) != -1){
        switch(c){
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr << "n_leaf_limit= " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr << "n_group_limit= " << n_group_limit << std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            std::cerr << "n_tot= " << n_tot << std::endl;
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

    PS::F64 full_len = 1.0;
    PS::F64vec offset = PS::F64vec(0.0);
    PS::ParticleSystem<FP> system;
    system.initialize();
    system.setAverageTargetNumberOfSampleParticlePerProcess(200);
    PS::S32 n_loc    = 0;
    REAL time_sys = 0.0;
    SetParticlesUniformQube(system, n_tot, n_loc, time_sys, full_len, offset, true); // random
    /*
    if(PS::Comm::getRank() == 1){
        for(PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++){
            std::cerr<<"system[i].pos= "<<system[i].pos<<std::endl;
        }
    }
    */
    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.collectSampleParticle(system);
    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);
    n_loc = system.getNumberOfParticleLocal();
    std::mt19937_64 gen(PS::Comm::getRank());
    std::uniform_real_distribution<PS::F64> dist(1.0, 3.0);
    for(PS::S32 i=0; i<n_loc; i++){
	system[i].r_search = full_len / cbrt((PS::F64)(n_tot)) * dist(gen);
        system[i].r_search *= 1.1;
        system[i].mass = 1.0 / n_tot;
    }

    PS::TreeForForceShort<Force, EP, EP>::Symmetry tree_symmetry;
    tree_symmetry.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    tree_symmetry.calcForceAll(SearchNeighborSymmetry(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
    /*
    if(PS::Comm::getRank()==0){
        for(PS::S32 i=0; i<n_tot; i += (n_tot/50)){
            EP * epj = tree_symmetry.getEpjFromId(i);
            std::cerr<<"i= "<<i;
            if(epj == NULL) std::cerr<<" it's null"<<std::endl;
            else std::cerr<<" id= "<<epj->id<<std::endl;
        }
    }
    */
    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(PS::S32 n_loop=0; n_loop<5; n_loop++){
        for(PS::S32 i=0; i<n_loc; i++){
            system[i].pos.x = full_len * mt.genrand_real2();
            system[i].pos.y = full_len * mt.genrand_real2();
            system[i].pos.z = full_len * mt.genrand_real2();
        }
        dinfo.decomposeDomainAll(system);
        system.exchangeParticle(dinfo);
        tree_symmetry.calcForceAll(SearchNeighborSymmetry(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
        PS::S32 n_cnt_epj = 0;
        for(PS::S32 i=0; i<n_tot; i++){
            EP * epj = tree_symmetry.getEpjFromId(i);
            if(epj != NULL){ 
                n_cnt_epj++;
                assert(i==epj->id);
            }
        }
        assert(n_cnt_epj == tree_symmetry.getNumberOfEpjSorted());
        if(PS::Comm::getRank()==0){
            std::cerr<<"n_cnt_epj= "<<n_cnt_epj
                     <<" tree_symmetry.getNumberOfEpjSorted()= "
                     <<tree_symmetry.getNumberOfEpjSorted()
                     <<std::endl;
        }
        /*
        if(PS::Comm::getRank()==0){
            std::cerr<<"****** n_loop= "<<n_loop<<std::endl;
            for(PS::S32 i=0; i<n_tot; i += (n_tot/50)){
                EP * epj = tree_symmetry.getEpjFromId(i);
                std::cerr<<"i= "<<i;
                if(epj == NULL) std::cerr<<" it's null"<<std::endl;
                else std::cerr<<" id= "<<epj->id<<" pos= "<<epj->pos<<std::endl;
            }
        }
        */
    }

    PS::Finalize();
    return 0;
}
