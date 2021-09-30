#pragma once

template<class Tpsys, class Tsat, class Tpla>
void CalcForceFromPlanet(Tpsys & system,
                         Tsat satellite[],
                         const PS::S32 n_sat,
                         const Tpla & planet){
    const PS::S32 n = system.getNumberOfParticleLocal();
    const PS::F64 mj = planet.mass;
    const PS::F64vec rj = planet.pos;
#pragma omp parallel
    for(PS::S32 i=0; i<n; i++) {
        const PS::F64vec rij = system[i].pos - rj;
        const PS::F64 r_sq = rij * rij;
        const PS::F64 r_inv = 1.0 / sqrt(r_sq);
        const PS::F64 pij = mj * r_inv;
        const PS::F64 mri3 = pij * r_inv * r_inv;
        system[i].acc -= mri3 * rij;
        system[i].pot -= pij * 2.0; // factor two is just a trick to calculate energy
    }
#pragma omp parallel
    for(PS::S32 i=0; i<n_sat; i++) {
        const PS::F64vec rij = satellite[i].pos - rj;
        const PS::F64 r_sq = rij * rij;
        const PS::F64 r_inv = 1.0 / sqrt(r_sq);
        const PS::F64 pij = mj * r_inv;
        const PS::F64 mri3 = pij * r_inv * r_inv;
        satellite[i].acc -= mri3 * rij;
        satellite[i].pot -= pij * 2.0; // factor two is just a trick to calculate energy
    }
}

class NgbSat{
public:
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    NgbSat():mass(0.0), pos(0.0), vel(0.0){}
    void clear(){
        pos = vel = 0.0;
        mass = 0.0;
    }
};

template<class Tpsys, class Tptcl, class Tforce>
void CalcForceBetweenDustAndSatellite(const Tpsys & system,
                                      Tptcl satellite[],
                                      Tforce force_dust[],
                                      const PS::S32 n_sat){
    const PS::S32 n = system.getNumberOfParticleLocal();
    struct ForceSat{
        PS::F64vec acc;
        PS::F64 pot;
        PS::F64 n_coll;
        ForceSat():acc(0.0), pot(0.0), n_coll(0.0){}
        void clear(){
            acc = 0.0;
            pot = n_coll = 0.0;
        }
    };

    ForceSat * force_sat_loc = new ForceSat[n_sat];
    ForceSat * force_sat_glb = new ForceSat[n_sat];
    NgbSat *   ngb_sat_loc = new NgbSat[n_sat];
    NgbSat *   ngb_sat_glb = new NgbSat[n_sat];
    PS::F64 * r_ngb_sq_sat = new PS::F64[n_sat];
    for(int i=0; i<n_sat; i++){
        r_ngb_sq_sat[i] = PS::LARGE_FLOAT;
        force_sat_loc[i].clear();
        force_sat_glb[i].clear();
        ngb_sat_loc[i].clear();
        ngb_sat_glb[i].clear();
    }
#pragma omp parallel
    for(PS::S32 i=0; i<n; i++) {
        const PS::F64vec ri = system[i].pos;
        const PS::F64 r_coll_i = system[i].r_coll;
        force_dust[i].clear();
        force_dust[i].r_ngb_sq = PS::LARGE_FLOAT;
        for(PS::S32 j=0; j<n_sat; j++){
            const PS::F64 r_coll_sq = (satellite[j].r_coll+r_coll_i)*(satellite[j].r_coll+r_coll_i);
            const PS::F64 mj = satellite[j].mass;
            const PS::F64vec rij = ri - satellite[j].pos;
            const PS::F64 r_sq = rij * rij;
            const PS::F64 r_inv = 1.0 / sqrt(r_sq);
            const PS::F64 pij = mj * r_inv;
            const PS::F64 mri3 = pij * r_inv * r_inv;
            const PS::F64vec aij = mri3 * rij;
            force_dust[i].acc -= aij;
            force_dust[i].pot -= pij;
            force_sat_loc[j].acc += aij;
            force_sat_loc[j].pot -= pij;
            if(r_coll_sq < rij*rij){
                force_dust[i].n_coll++;
                force_sat_loc[j].n_coll += 1.0+1e-10;
                if(rij*rij < force_dust[i].r_ngb_sq){
                    // for dust
                    force_dust[i].r_ngb_sq = rij*rij;
                    //force_dust[i].ngb.mass = satellite[j].mass;
                    //force_dust[i].ngb.pos = satellite[j].pos;
                    //force_dust[i].ngb.vel = satellite[j].vel;
                }
                if(rij*rij < r_ngb_sq_sat[j]){
                    r_ngb_sq_sat[j] = rij*rij;
                    ngb_sat_loc[j].mass = system[i].mass;
                    ngb_sat_loc[j].pos = system[i].pos;
                    ngb_sat_loc[j].vel = system[i].vel;
                }
            }
        }
    }
    MPI_Allreduce( (double*)force_sat_loc,  (double*)force_sat_glb, 5*n_sat, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    class RNgbSqRank{
    public:
        double r_ngb_sq;
        int rank;
    };
    RNgbSqRank * r_ngb_sq_rank_in  = new RNgbSqRank[n_sat];
    RNgbSqRank * r_ngb_sq_rank_out = new RNgbSqRank[n_sat];
    for(PS::S32 i=0; i<n_sat; i++){
        r_ngb_sq_rank_in[i].r_ngb_sq = r_ngb_sq_sat[i];
        r_ngb_sq_rank_in[i].rank = PS::Comm::getRank();
    }

    for(PS::S32 i=0; i<n_sat; i++){
        MPI_Allreduce(r_ngb_sq_rank_in+i, r_ngb_sq_rank_out+i, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        MPI_Bcast(ngb_sat_loc+i, 1, PS::GetDataType<NgbSat>(), r_ngb_sq_rank_out[i].rank, MPI_COMM_WORLD);
        //std::cerr<<"ngb_sat_loc[i].pos= "<<ngb_sat_loc[i].pos<<std::endl;
    }
    for(PS::S32 i=0; i<n_sat; i++){
        satellite[i].acc = force_sat_glb[i].acc;
        satellite[i].pot = force_sat_glb[i].pot;
    }
    delete [] r_ngb_sq_rank_in;
    delete [] r_ngb_sq_rank_out;
    delete [] ngb_sat_loc;
    delete [] ngb_sat_glb;
    delete [] force_sat_loc;
    delete [] force_sat_glb;
}

template<class Tptcl>
void kick(Tptcl ptcl[],
          const PS::S32 n,
          const PS::F64 dt) {
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++) {
        ptcl[i].vel  += ptcl[i].acc * dt;
    }
}

template<class Tptcl>
void drift(Tptcl ptcl[],
           const PS::S32 n,
           const PS::F64 dt) {
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++) {
        ptcl[i].pos  += ptcl[i].vel * dt;
    }
}

template<class Tpsys>
void kick(Tpsys & system,
          const PS::F64 dt) {
    const PS::S32 n = system.getNumberOfParticleLocal();
    kick(&system[0], n, dt);
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    const PS::S32 n = system.getNumberOfParticleLocal();
    drift(&system[0], n, dt);
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
        epot_loc += system[i].mass * system[i].pot;
        //epot_loc += system[i].mass * (system[i].pot + system[i].mass / Epi::eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
#else
    etot = etot_loc;
    epot = epot_loc;
    ekin = ekin_loc;
#endif
}

template<class Tpsys, class Tptcl>
void calcEnergy(const Tpsys & system,
                const Tptcl satellite,
                const PS::S32 n_sat,
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
        epot_loc += system[i].mass * system[i].pot;
        //epot_loc += system[i].mass * (system[i].pot + system[i].mass / Epi::eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
    for(PS::S32 i=0; i<n_sat; i++){
        ekin += 0.5 * satellite[i].mass * satellite[i].vel * satellite[i].vel;
        epot += 0.5 * satellite[i].mass * satellite[i].pot;
    }
    etot = ekin + epot;
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
    if(stat(dir_name, &st) != 0) {
        PS::S32 ret_loc = 0;
        PS::S32 ret     = 0;
        if(PS::Comm::getRank() == 0)
            ret_loc = mkdir(dir_name, 0777);
        PS::Comm::broadcast(&ret_loc, ret);
        if(ret == 0) {
            if(PS::Comm::getRank() == 0)
                fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
        } else {
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
            PS::Abort();
        }
    }
}
