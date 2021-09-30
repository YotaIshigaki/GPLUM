#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>

#ifdef USEPHANTOMGRAPE
//#include"phantomgrape.hpp"
#include"phantomquad.hpp"
//#include"phantomquad_x86.hpp"
#endif

class FileHeader{
public:
    PS::S64 n_tot_glb;
    PS::S64 n_disk_glb;
    PS::S64 n_bulge_glb;
    PS::S64 n_dark_glb;
    PS::S64 n_tot_loc;
    PS::F64 time;
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%e\n", time);
        fprintf(fp, "%lld\n", n_tot_glb);
        fprintf(fp, "%lld\n", n_disk_glb);
        fprintf(fp, "%lld\n", n_bulge_glb);
        fprintf(fp, "%lld\n", n_dark_glb);
        fprintf(fp, "%lld\n", n_tot_loc);
    }
};


class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64 pot;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

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
    void copyFromForce(const ForceGrav & force){
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

class EPIGrav{
public:
    PS::S64 id;
    PS::F64vec pos;
    static PS::F64 eps;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FPGrav & fp){ 
        pos = fp.pos;
        id = fp.id;
    }
};

//PS::F64 EPIGrav::eps = 1.0/32.0;
PS::F64 EPIGrav::eps = 1e-3; // 1PC

class EPJGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    void copyFromFP(const FPGrav & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getCharge() const { return mass; }
};

#ifdef USEPHANTOMGRAPE
struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
        }
        for(PS::S32 i=0; i<n_jp; i++){
            const PS::F64 m_j = ep_j[i].getCharge();
            const PS::F64vec pos_j = ep_j[i].getPos();
            pg.set_epj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j);
        }
        pg.run_epj(n_ip, n_jp);
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64 * p = &(force[i].pot);
            PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
            pg.get_accp_one(i, a[0], a[1], a[2], *p);
        }
    }
};

#else
struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){

        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S64 idi = ep_i[i].id;
            for(PS::S32 j=0; j<n_jp; j++){
                if( idi == ep_j[j].id ) continue;
                PS::F64vec rij = xi - ep_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#endif // USEPHANTOMGRAPE



#ifdef USEPHANTOMGRAPE
#ifdef MONOPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){


        if(n_ip > PhantomGRAPE::PG_NIMAX){
            std::cout<<"n_ip(ep-sp)="<<n_ip<<std::endl;
            std::cout<<"PhantomGRAPE::PG_NIMAX(ep-sp)="<<PhantomGRAPE::PG_NIMAX<<std::endl;
        }
        if(n_jp > PhantomGRAPE::PG_NJMAX){
            std::cout<<"n_jp(ep-sp)="<<n_jp<<std::endl;
            std::cout<<"PhantomGRAPE::PG_NJMAX(ep-sp)="<<PhantomGRAPE::PG_NJMAX<<std::endl;
        }

        double xi[PhantomGRAPE::PG_NIMAX][3];
        double mxj[PhantomGRAPE::PG_NJMAX][4];
        double ai[PhantomGRAPE::PG_NIMAX][3];
        double pi[PhantomGRAPE::PG_NIMAX];
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            xi[i][0] = ep_i[i].pos.x;
            xi[i][1] = ep_i[i].pos.y;
            xi[i][2] = ep_i[i].pos.z;
        }
        for(PS::S32 j=0; j<n_jp; j++){
            mxj[j][0] = sp_j[j].getPos().x;
            mxj[j][1] = sp_j[j].getPos().y;
            mxj[j][2] = sp_j[j].getPos().z;
            mxj[j][3] = sp_j[j].getCharge();
        }
        PhantomGRAPE pg;
        pg.set_eps2(eps2);
        pg.set_xj(n_jp, mxj);
        pg.set_xi(n_ip, xi);
        pg.run(n_ip, n_jp);
        pg.get_ai(n_ip, ai, pi);

        for(PS::S32 i=0; i<n_ip; i++){
            force[i].acc.x += ai[i][0];
            force[i].acc.y += ai[i][1];
            force[i].acc.z += ai[i][2];
            force[i].pot += pi[i];
        }

    }
};
#elif QUADRUPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJQuadrupole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){

        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
        }

        for(PS::S32 j=0; j<n_jp; j++){
            const PS::F64 m_j = sp_j[j].getCharge();
            const PS::F64vec pos_j = sp_j[j].getPos();
            const PS::F64mat q = sp_j[j].quad;
            pg.set_spj_one(j, pos_j.x, pos_j.y, pos_j.z, m_j,
                           q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
        }
        pg.run_spj(n_ip, n_jp);
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64 * p = &(force[i].pot);
            PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
            pg.accum_accp_one(i, a[0], a[1], a[2], *p);
        }
    }
};
#endif // MONOPOLE QUADRUPOLE
#else // USEPHANTOMGRAPE

#ifdef MONOPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#elif QUADRUPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJQuadrupole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){



        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].mass;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};
#endif // MONOPOLE
#endif // USEPHANTOMGRAPE

template<class Tpsys>
void WriteNemoAscii(const Tpsys & psys,
                    const PS::F32 time_sys,
                    const PS::S32 snp_id,
                    const char * dir_name){
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::S64 n_glb = 0;
    FPGrav * fp;
    PS::AllGatherParticle(fp, n_glb, &psys[0], n_loc);
    if(PS::Comm::getRank () == 0){
        const PS::S32 STRINGSIZE = 1024;
        char sout[STRINGSIZE];
        sprintf(sout,"%s/snap%5d.dat", dir_name, snp_id);
        for(int i=0;i<STRINGSIZE;i++)if(sout[i]==' ')sout[i]='0';
        std::ofstream foutput;
        foutput.open(sout);
        foutput<<std::setprecision(15);
        foutput<<n_glb<<std::endl;
        foutput<<"3"<<std::endl;
        foutput<<time_sys<<std::endl;
        for(PS::S64 i=0; i<n_glb; i++) foutput<<fp[i].mass<<std::endl;
        for(PS::S64 i=0; i<n_glb; i++) foutput<<fp[i].pos<<std::endl;
        for(PS::S64 i=0; i<n_glb; i++) foutput<<fp[i].vel<<std::endl;
        foutput.close();
    }
    delete [] fp;
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
    std::cout<<"r_max= "<<sqrt(r_max_sq)<<std::endl;
}


template<class Tpsys>
void SetParticlesPlummer(Tpsys & psys,
                         const PS::S64 n_glb,
                         PS::S64 & n_loc,  
                         PS::F64 & t_sys){

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

template<class Tpsys>
void Drift(Tpsys & system,
           const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].pos  += system[i].vel * dt;
    }
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
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
    for(PS::S32 i=0; i<nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * system[i].pot;
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
}

void ReadSnapShot(PS::ParticleSystem<FPGrav> & sys, 
                  const PS::S64 & n_factor,
                  PS::F64 & time_sys,
                  PS::S64 & n_tot_glb,
                  PS::S64 & n_tot_loc,
                  PS::S64 & n_disk_glb,
                  PS::S64 & n_bulge_glb,
                  PS::S64 & n_dark_glb,
                  const char * input_dir_name, 
                  const PS::S64 n_snp_init,
                  const PS::S64 n_snp_init_limit=1024){

    std::cout<<"n_snp_init="<<n_snp_init<<std::endl;

    const PS::S32 n_proc_glb = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank_glb = PS::Comm::getRank();
    assert( n_proc_glb >= n_snp_init);
    assert(n_proc_glb % n_snp_init == 0);
    MPI_Comm comm_sub;
    int color = my_rank_glb % n_snp_init;

    std::cout<<"color"<<std::endl;

    MPI_Comm_split(MPI_COMM_WORLD, color, my_rank_glb, &comm_sub);
    int my_rank_sub, n_proc_sub;
    MPI_Comm_rank(comm_sub, &my_rank_sub);
    MPI_Comm_size(comm_sub, &n_proc_sub);
    std::cout<<"n_proc_sub="<<n_proc_sub<<" my_rank_sub="<<my_rank_sub<<std::endl;


    char sinput[1024];
    std::ifstream fin;

    PS::S64 n_tot_glb_snp, n_disk_glb_snp, n_bulge_glb_snp, n_dark_glb_snp;
    PS::S64 n_tot_loc_org = 0;
    if( PS::Comm::getRank() < n_snp_init ){
        PS::S64 XXX;
        sprintf(sinput, "%s/disk_%d_%d.dat", input_dir_name, n_snp_init_limit, PS::Comm::getRank());
        std::cout<<"sinput: "<<sinput<<std::endl;
        fin.open(sinput);
        fin>>n_tot_glb_snp;
        fin>>n_disk_glb_snp;
        fin>>n_bulge_glb_snp;
        fin>>n_dark_glb_snp;
        fin>>XXX;
        n_tot_loc_org += XXX;
        fin>>XXX;
        n_tot_loc_org += XXX;
        fin>>XXX;
        n_tot_loc_org += XXX;
        fin>>time_sys;
    }
    else{
        n_tot_loc_org = 0;
    }
    PS::Comm::broadcast(&n_tot_glb_snp, 1, 0);
    PS::Comm::broadcast(&n_disk_glb_snp, 1, 0);
    PS::Comm::broadcast(&n_bulge_glb_snp, 1, 0);
    PS::Comm::broadcast(&n_dark_glb_snp, 1, 0);
    PS::Comm::broadcast(&time_sys, 1, 0);
    MPI_Bcast(&n_tot_loc_org, 1, PS::GetDataType<PS::S64>(), 0, comm_sub);
    //PS::S64 n_tot_glb_org = PS::Comm::getSum(n_tot_loc_org); // at first evaluate n_tot_glb_org befor broadcast n_tot_loc_org
    //n_tot_glb = n_tot_glb_org * n_factor;
    //std::cerr<<"time_sys="<<time_sys<<std::endl;
    //std::cout<<"n_tot_loc_org="<<n_tot_loc_org<<std::endl;


    PS::S64 * id = new PS::S64[n_tot_loc_org];
    PS::F64 * mass = new PS::F64[n_tot_loc_org];
    PS::F64vec * pos = new PS::F64vec[n_tot_loc_org];
    PS::F64vec * vel = new PS::F64vec[n_tot_loc_org];
    if( PS::Comm::getRank() < n_snp_init ){
        for(long i=0; i<n_tot_loc_org; i++){
            fin >> id[i] >> mass[i] >> pos[i] >> vel[i];
        }
        fin.close();
    }
    MPI_Bcast(id, n_tot_loc_org, PS::GetDataType<PS::S64>(), 0, comm_sub);
    MPI_Bcast(mass, n_tot_loc_org, PS::GetDataType<PS::F64>(), 0, comm_sub);
    MPI_Bcast(pos, n_tot_loc_org, PS::GetDataType<PS::F64vec>(), 0, comm_sub);
    MPI_Bcast(vel, n_tot_loc_org, PS::GetDataType<PS::F64vec>(), 0, comm_sub);

    //std::cout<<"finish of coping of snap shot. each proc has n_tot_loc_org ptcls"<<std::endl;

    n_tot_loc = (n_tot_loc_org * n_factor) / n_proc_sub;
    PS::S64 id_head = n_tot_loc * my_rank_sub;
    if(my_rank_sub < (n_tot_loc_org * n_factor) % n_proc_sub ){
        n_tot_loc++;
        id_head += my_rank_sub;
    }
    else{
        id_head += (n_tot_loc_org * n_factor) % n_proc_sub;
    }
    //std::cout<<"my_rank_sub="<<my_rank_sub<<" n_proc_sub="<<n_proc_sub<<std::endl;
    //std::cout<<"id_head="<<id_head<<" n_tot_loc="<<n_tot_loc<<std::endl;

    sys.initialize();
    sys.createParticle(n_tot_loc*3+100);
    sys.setNumberOfParticleLocal(n_tot_loc);
    n_tot_glb = PS::Comm::getSum(n_tot_loc);
    PS::F64 m_factor = (PS::F64)n_tot_glb_snp / (PS::F64)n_tot_glb;
    std::cout<<"m_factor="<<m_factor<<std::endl;
    PS::MTTS mt;
    mt.init_genrand(my_rank_glb);
    PS::S64 cnt = 0;
    for(PS::S64 i=id_head; i<id_head+n_tot_loc; i++){
        PS::S64 i_tmp = i % n_tot_loc_org;
        PS::F64 theta = 2.0 * M_PI * mt.genrand_real2();
        PS::F64 cth = cos(theta);
        PS::F64 sth = sin(theta);
        //PS::F64 cth = 1;
        //PS::F64 sth = 0;
        sys[cnt].id = id[i_tmp];
        sys[cnt].mass = mass[i_tmp] * m_factor;
        sys[cnt].pos.x = cth * pos[i_tmp].x - sth * pos[i_tmp].y;
        sys[cnt].pos.y = sth * pos[i_tmp].x + cth * pos[i_tmp].y;
        sys[cnt].vel.x = cth * vel[i_tmp].x - sth * vel[i_tmp].y;
        sys[cnt].vel.y = sth * vel[i_tmp].x + cth * vel[i_tmp].y;
        sys[cnt].pos.z = pos[i_tmp].z;
        sys[cnt].vel.z = vel[i_tmp].z;
        cnt++;
    }

    assert(cnt == n_tot_loc);

    PS::S64 n_disk_loc = 0;
    PS::S64 n_bulge_loc = 0;
    PS::S64 n_dark_loc = 0;
    for(PS::S64 i=0; i<n_tot_loc; i++){
        if( sys[i].id < n_disk_glb_snp) n_disk_loc++;
        else if( sys[i].id < n_disk_glb_snp + n_bulge_glb_snp) n_bulge_loc++;
        else n_dark_loc++;
    }

    n_disk_glb = PS::Comm::getSum(n_disk_loc);
    n_bulge_glb = PS::Comm::getSum(n_bulge_loc);
    n_dark_glb = PS::Comm::getSum(n_dark_loc);

    std::cout<<"n_disk_glb="<<n_disk_glb<<std::endl;
    std::cout<<"n_bulge_glb="<<n_bulge_glb<<std::endl;
    std::cout<<"n_dark_glb="<<n_dark_glb<<std::endl;

    assert( n_disk_glb+n_bulge_glb+n_dark_glb == n_tot_glb);

    PS::S64 * n_disk_array = new PS::S64[n_proc_glb];
    PS::S64 * n_bulge_array = new PS::S64[n_proc_glb];
    PS::S64 * n_dark_array = new PS::S64[n_proc_glb];

    PS::Comm::allGather(&n_disk_loc, 1, n_disk_array);
    PS::Comm::allGather(&n_bulge_loc, 1, n_bulge_array);
    PS::Comm::allGather(&n_dark_loc, 1, n_dark_array);

    PS::S64 * n_disk_array_disp = new PS::S64[n_proc_glb+1];
    PS::S64 * n_bulge_array_disp = new PS::S64[n_proc_glb+1];
    PS::S64 * n_dark_array_disp = new PS::S64[n_proc_glb+1];

    n_disk_array_disp[0] = n_bulge_array_disp[0] = n_dark_array_disp[0] = 0;
    for(PS::S64 i=0; i<n_proc_glb; i++){
        n_disk_array_disp[i+1] = n_disk_array_disp[i] + n_disk_array[i];
        n_bulge_array_disp[i+1] = n_bulge_array_disp[i] + n_bulge_array[i];
        n_dark_array_disp[i+1] = n_dark_array_disp[i] + n_dark_array[i];
    }

    std::cout<<"n_disk_array_disp[n_proc_glb]="<<n_disk_array_disp[n_proc_glb]<<std::endl;
    std::cout<<"n_bulge_array_disp[n_proc_glb]="<<n_bulge_array_disp[n_proc_glb]<<std::endl;
    std::cout<<"n_dark_array_disp[n_proc_glb]="<<n_dark_array_disp[n_proc_glb]<<std::endl;

    PS::S64 id_disk = n_disk_array_disp[my_rank_glb];
    PS::S64 id_bulge = n_disk_array_disp[n_proc_glb] + n_bulge_array_disp[my_rank_glb];
    PS::S64 id_dark = n_disk_array_disp[n_proc_glb] + n_bulge_array_disp[n_proc_glb] + n_dark_array_disp[my_rank_glb];

    std::cout<<"id_disk="<<id_disk<<std::endl;
    std::cout<<"id_bulge="<<id_bulge<<std::endl;
    std::cout<<"id_dark="<<id_dark<<std::endl;

    for(PS::S64 i=0; i<n_tot_loc; i++){
        if( sys[i].id < n_disk_glb_snp){
            sys[i].id = id_disk;
            id_disk++;
        }
        else if( sys[i].id < n_disk_glb_snp + n_bulge_glb_snp){
            sys[i].id = id_bulge;
            id_bulge++;
        }
        else{
            sys[i].id = id_dark;
            id_dark++;
        }
    }

    std::cout<<"id_disk2="<<id_disk<<std::endl;
    std::cout<<"id_bulge2="<<id_bulge<<std::endl;
    std::cout<<"id_dark2="<<id_dark<<std::endl;

    PS::S64 cnt_disk = 0;
    PS::S64 cnt_bulge = 0;
    PS::S64 cnt_dark = 0;
    for(PS::S64 i=0; i<n_tot_loc; i++){
        if( sys[i].id < n_disk_glb) cnt_disk++;
        else if( sys[i].id < n_disk_glb + n_bulge_glb) cnt_bulge++;
        else cnt_dark++;
    }

    std::cout<<"cnt_disk="<<cnt_disk<<std::endl;
    std::cout<<"cnt_bulge="<<cnt_bulge<<std::endl;
    std::cout<<"cnt_dark="<<cnt_dark<<std::endl;

    assert(cnt_disk == n_disk_loc);
    assert(cnt_bulge == n_bulge_loc);
    assert(cnt_dark == n_dark_loc);

    PS::F64 mass_CM_loc = 0.0;
    PS::F64vec pos_CM_loc = 0.0;
    PS::F64vec vel_CM_loc = 0.0;

    for(PS::S64 i=0; i<n_tot_loc; i++){
        mass_CM_loc += sys[i].mass;
        pos_CM_loc += sys[i].mass * sys[i].pos;
        vel_CM_loc += sys[i].mass * sys[i].vel;
    }

    std::cerr<<"mass_CM_loc="<<mass_CM_loc<<std::endl;
    std::cerr<<"pos_CM_loc="<<pos_CM_loc<<std::endl;
    std::cerr<<"vel_CM_loc="<<vel_CM_loc<<std::endl;

    PS::F64 mass_CM_glb = PS::Comm::getSum(mass_CM_loc);
    PS::F64vec pos_CM_glb = PS::Comm::getSum(pos_CM_loc);
    PS::F64vec vel_CM_glb = PS::Comm::getSum(vel_CM_loc);

    pos_CM_glb /= mass_CM_glb;
    vel_CM_glb /= mass_CM_glb;

    std::cerr<<"mass_CM_glb="<<mass_CM_glb<<std::endl;
    std::cerr<<"pos_CM_glb="<<pos_CM_glb<<std::endl;
    std::cerr<<"vel_CM_glb="<<vel_CM_glb<<std::endl;

    for(PS::S64 i=0; i<n_tot_loc; i++){
        sys[i].pos -= pos_CM_glb;
        sys[i].vel -= vel_CM_glb;
    }

#if 0
    mass_CM_loc = 0.0;
    pos_CM_loc = 0.0;
    vel_CM_loc = 0.0;
    for(PS::S64 i=0; i<n_tot_loc; i++){
        mass_CM_loc += sys[i].mass;
        pos_CM_loc += sys[i].mass * sys[i].pos;
        vel_CM_loc += sys[i].mass * sys[i].vel;
    }
    mass_CM_glb = PS::Comm::getSum(mass_CM_loc);
    pos_CM_glb = PS::Comm::getSum(pos_CM_loc);
    vel_CM_glb = PS::Comm::getSum(vel_CM_loc);
    pos_CM_glb /= mass_CM_glb;
    vel_CM_glb /= mass_CM_glb;
    std::cerr<<"pos_CM_glb="<<pos_CM_glb<<std::endl;
    std::cerr<<"vel_CM_glb="<<vel_CM_glb<<std::endl;
#endif

    std::cout<<"before delete"<<std::endl;

    delete [] n_disk_array;
    delete [] n_bulge_array;
    delete [] n_dark_array;
    delete [] n_disk_array_disp;
    delete [] n_bulge_array_disp;
    delete [] n_dark_array_disp;

    delete [] id;
    delete [] mass;
    delete [] pos;
    delete [] vel;

    std::cout<<"after delete"<<std::endl;

}

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);

    PS::F64 Tbegin = PS::GetWtime();
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 n_smp_ave = 30;
    PS::F32 time_end = 10.0;
    PS::S64 factor = 1;
    PS::S64 n_snp_init = 1;
    char input_dir_name[1024];
    char output_dir_name[1024];
    int c;
    while((c=getopt(argc,argv,"i:o:t:T:n:N:f:s:l:h")) != -1){
        switch(c){
        case 'i':
            sprintf(input_dir_name, optarg);
            break;
        case 'o':
            sprintf(output_dir_name, optarg);
            break;
        case 't':
            theta = atof(optarg);
            //std::cerr<<"theta="<<theta<<std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            //std::cerr<<"time_end="<<time_end<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            //std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'N':
            n_snp_init = atoi(optarg);
            //std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'f':
            factor = atol(optarg);
            //std::cerr<<"factor="<<factor<<std::endl;
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            //std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            //std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'h':
            std::cerr<<"i: input file name (nemo ascii)"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 10.0)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"N: n_snp_init (# of initial snap shot used. dafult: 1)"<<std::endl;
            std::cerr<<"f: n_factor (must >= 1)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 30)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            return 0;
        }
    }
    std::ofstream fout_eng;
    std::ofstream fout_tcal;
    std::ofstream fout_log;

#if 0
    if(PS::Comm::getRank() == 0){
#else
    if(1){
#endif
        char sout_de[1024];
        char sout_tcal[1024];
        char sout_log[1024];
        sprintf(sout_de, "%s/t-de_%05d_%05d.dat", output_dir_name, PS::Comm::getNumberOfProc(), PS::Comm::getRank() );
        sprintf(sout_tcal, "%s/t-tcal_%05d_%05d.dat", output_dir_name, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
        sprintf(sout_log, "%s/log_%05d_%05d.dat", output_dir_name, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
        if(PS::Comm::getRank() == 0){
            std::cerr<<sout_de<<std::endl;
            std::cerr<<sout_tcal<<std::endl;
            std::cerr<<sout_log<<std::endl;
        }
        fout_eng.open(sout_de);
        fout_tcal.open(sout_tcal);
        fout_log.open(sout_log);
    }
    fout_log<<"Comm::getNumberOfProc()="<<PS::Comm::getNumberOfProc()<<std::endl;
    fout_log<<"Comm::getNumberOfThread()="<<PS::Comm::getNumberOfThread()<<std::endl;
    fout_log<<"output_dir_name="<<output_dir_name<<std::endl;
    fout_log<<"input_dir_name="<<input_dir_name<<std::endl;
    fout_log<<"n_snp_init="<<n_snp_init<<std::endl;
    fout_log<<"theta="<<theta<<std::endl;
    fout_log<<"time_end="<<time_end<<std::endl;
    fout_log<<"n_group_limit="<<n_group_limit<<std::endl;
    fout_log<<"n_smp_ave="<<n_smp_ave<<std::endl;
    fout_log<<"n_leaf_limit="<<n_leaf_limit<<std::endl;

    PS::ParticleSystem<FPGrav> system_grav;

    PS::F64 time_sys = 0.0;
    PS::S64 n_tot_glb = 0;
    PS::S64 n_tot_loc = 0;
#ifdef USE_PLUMMER
    n_tot_glb = 1048576 * PS::Comm::getNumberOfProc();
    system_grav.initialize();
    SetParticlesPlummer(system_grav, n_tot_glb, n_tot_loc, time_sys);
#else
    PS::S64 n_disk_glb = 0;
    PS::S64 n_bulge_glb = 0;
    PS::S64 n_dark_glb = 0;
    PS::S64 n_snp_init_limit = 1024;
    ReadSnapShot(system_grav, factor, time_sys, n_tot_glb, n_tot_loc, n_disk_glb, n_bulge_glb, n_dark_glb, input_dir_name, n_snp_init, n_snp_init_limit);
#endif
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);

    fout_log<<"finish set particle: time="<<PS::GetWtime() - Tbegin<<std::endl;

    fout_log<<"n_tot_glb="<<n_tot_glb<<std::endl;
    fout_log<<"n_tot_loc="<<n_tot_loc<<std::endl;
#ifndef USE_PLUMMER
    fout_log<<"n_disk_glb="<<n_disk_glb<<std::endl;
    fout_log<<"n_bulge_glb="<<n_bulge_glb<<std::endl;
    fout_log<<"n_dark_glb="<<n_dark_glb<<std::endl;
#endif

#if 1

    const PS::F32 coef_ema = 0.4;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    PS::S32 nx = 48;
    PS::S32 ny = 54;
    PS::S32 nz = 32;
#ifdef USE_KFULL
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
#elif defined(USE_K8)
    dinfo.setNumberOfDomainMultiDimension(nx/2, ny/2, nz/2);
#elif defined(USE_K64)
    dinfo.setNumberOfDomainMultiDimension(nx/4, ny/4, nz/4);
#endif
    fout_log<<"finish dinof.initialize="<<PS::GetWtime() - Tbegin<<std::endl;

    dinfo.collectSampleParticle(system_grav, true);

    fout_log<<"finish dinof.collect="<<PS::GetWtime() - Tbegin<<std::endl;

    dinfo.decomposeDomain();
    fout_log<<"dinfo.getPosDomain(PS::Comm::getRank())="<<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    fout_log<<"finish dinof.decompose="<<PS::GetWtime() - Tbegin<<std::endl;


    system_grav.exchangeParticle(dinfo);
    //system_grav.exchangeParticleSafely(dinfo);
    //system_grav.exchangeParticleSafeMode(dinfo);

    fout_log<<"finish dinof.exchangeparticle="<<PS::GetWtime() - Tbegin<<std::endl;

    PS::S64 n_grav_loc = system_grav.getNumberOfParticleLocal();
    PS::S64 n_grav_glb = PS::Comm::getSum(n_grav_loc);
#ifdef MONOPOLE
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Monopole tree_grav;
#elif QUADRUPOLE
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Quadrupole tree_grav;
#elif MONOPOLEGEOMETRICCENTER
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::MonopoleGeometricCenter tree_grav;
#elif DIPOLEGEOMETRICCENTER
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::DipoleGeometricCenter tree_grav;
#elif QUADRUPOLEGEOMETRICCENTER
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::QuadrupoleGeometricCenter tree_grav;
#endif

    tree_grav.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    
    fout_log<<"finish tree_grav.initialize="<<PS::GetWtime() - Tbegin<<std::endl;

    fout_log<<"tree_grav.getMemSizeUsed()= "<<tree_grav.getMemSizeUsed()*1e-9<<" [Gbyte]";
    fout_log<<" system_grav.getMemSizeUsed()= "<<system_grav.getMemSizeUsed()*1e-9<<" [Gbyte]"<<std::endl;

    //tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);
    tree_grav.calcForceAllAndWriteBackWithCheck(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true, fout_log);

    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    CalcEnergy(system_grav, Etot0, Ekin0, Epot0);

    fout_log<<"Epot0= "<<Epot0<<" Ekin0= "<<Ekin0<<" Etot0= "<<Etot0<<std::endl;

    const PS::F32 dt = 1.0/128.0; // 7.66e4 yr

    Kick(system_grav, dt*0.5);

    PS::F64 Tloop = 0.0;

    PS::S32 snp_id = 0;
    PS::S64 n_loop = 0;

    PS::Timer timer;
    while(time_sys < time_end){
        //if( fmod(time_sys, 0.5) == 0.0){
        if( n_loop % 8 == 0){
            timer.initialize(fout_tcal);
        }
        timer.reset();
        timer.start();

#ifdef GET_SNAPSHOT
        if( fmod(time_sys, 1.0) == 0.0){
            FileHeader header;
            header.time = time_sys;
            header.n_tot_glb = n_tot_glb;
            header.n_disk_glb = n_disk_glb;
            header.n_bulge_glb = n_bulge_glb;
            header.n_dark_glb = n_dark_glb;
            header.n_tot_loc = system_grav.getNumberOfParticleLocal();
            char filename[256];
            sprintf(filename, "%s/snap_%d", output_dir_name, snp_id++);
            system_grav.writeParticleAscii(filename, "%s_%d_%d.dat", header);
        }
#endif
        timer.restart("WriteNemoAscii");


        time_sys += dt;
        Drift(system_grav, dt);

        timer.restart("Drift");
        if( n_loop < 4){
            dinfo.collectSampleParticle(system_grav);
            timer.restart("collect");
            dinfo.decomposeDomain();
        }
        if( n_loop > 0 && (n_loop < 12 || n_loop % 8 == 0) ){
            dinfo.collectSampleParticle(system_grav, true, Tloop);
            timer.restart("collect");
            //dinfo.decomposeDomain();
            dinfo.decomposeDomainMultiStep();
        }
        else{
            timer.restart("collect");
        }
/*            
        //if( fmod(time_sys, 1.0/16.0) == 0.0 && n_loop > 1){
        if( n_loop % 8 == 0 && n_loop > 1){
#ifdef USE_EMA 
            dinfo.collectSampleParticle(system_grav, true, Tloop);
#else
            dinfo.collectSampleParticle(system_grav);
#endif
            timer.restart("collect");
            dinfo.decomposeDomain();
        }
        else{
            timer.restart("collect");
        }
*/
        timer.restart("decompose");
            
        system_grav.exchangeParticle(dinfo);

        timer.restart("exchangeParticle");

        Tloop = PS::GetWtime();
                
        //tree_grav.calcForceAllAndWriteBackWithTimer(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, timer, true);
        tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true);
        
        Tloop = PS::GetWtime() - Tloop;
        
        Kick(system_grav, dt*0.5);
            
        timer.stop("Kick");
        
#ifdef GET_SNAPSHOT
        if( fmod(time_sys, 1.0) == 0.0){
            fout_tcal<<"time_sys= "<<time_sys<<" n_loop="<<n_loop<<std::endl;
            fout_tcal<<"tree_grav.getMemSizeUsed()= "<<tree_grav.getMemSizeUsed()*1e-9<<" [Gbyte]";
            fout_tcal<<" system_grav.getMemSizeUsed()= "<<system_grav.getMemSizeUsed()*1e-9<<" [Gbyte]"<<std::endl;
            //tree_grav.dump_calc_cost(PS::Comm::getMaxValue(Tloop), fout_tcal);
            fout_tcal<<"Tloop= "<<Tloop<<" Ttot="<<PS::GetWtime()-Tbegin<<std::endl;
            timer.dump(fout_tcal);
            fout_tcal<<std::endl;

            CalcEnergy(system_grav, Etot1, Ekin1, Epot1);
            if(PS::Comm::getRank() == 0){
                fout_eng<<time_sys<<"   "<<(Etot1-Etot0)/Etot0<<std::endl;
            }
        }
#else
        fout_tcal<<"time_sys= "<<time_sys<<" n_loop="<<n_loop<<std::endl;
        fout_tcal<<"tree_grav.getMemSizeUsed()= "<<tree_grav.getMemSizeUsed()*1e-9<<" [Gbyte]";
        fout_tcal<<" system_grav.getMemSizeUsed()= "<<system_grav.getMemSizeUsed()*1e-9<<" [Gbyte]"<<std::endl;
        //tree_grav.dump_calc_cost(PS::Comm::getMaxValue(Tloop), fout_tcal);
        fout_tcal<<"Tloop= "<<Tloop<<" Ttot="<<PS::GetWtime()-Tbegin<<std::endl;
        timer.dump(fout_tcal);
        fout_tcal<<std::endl;

        CalcEnergy(system_grav, Etot1, Ekin1, Epot1);
        if(PS::Comm::getRank() == 0){
            fout_eng<<time_sys<<"   "<<(Etot1-Etot0)/Etot0<<std::endl;
        }
#endif
        Kick(system_grav, dt*0.5);
        n_loop++;
    }

#endif

    PS::Finalize();
    return 0;
}
