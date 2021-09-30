#pragma once

class Collision0{
 public:
    PS::F64 time;

    PS::S64 id_imp;
    PS::S64 id_tar;
    PS::S32 id_c_imp;
    PS::S32 id_c_tar;

    PS::F64vec pos_imp;
    PS::F64vec pos_tar;
    PS::F64vec pos_g;
    PS::F64vec vel_imp;
    PS::F64vec vel_tar;
    PS::F64vec vel_g;
    PS::F64 col_angle;
    
    PS::F64vec pos_imp_new;
    PS::F64vec pos_tar_new;
    PS::F64vec pos_g_new;
    PS::F64vec vel_imp_new;
    PS::F64vec vel_tar_new;
    PS::F64vec vel_g_new;

    PS::F64 mass_imp;
    PS::F64 mass_tar;
    PS::F64 mass_frag;

    PS::F64 r_planet_imp;
    PS::F64 r_planet_tar;
    PS::F64 f_imp;
    PS::F64 f_tar;
    PS::F64 r_planet_imp_new;
    PS::F64 r_planet_tar_new;
    PS::F64 f_imp_new;
    PS::F64 f_tar_new;
#ifdef MERGE_BINARY
    bool imp_is_binary;
    bool tar_is_binary;
    bool imp_will_be_binary;
    bool tar_will_be_binary;
#endif
    
    PS::F64 edisp;
    PS::F64 edisp_d;

    PS::S64 id_frag;
    PS::S32 id_c_frag;
    PS::S32 n_frag;

    //PS::S32 HitAndRun;
    PS::S32 flag_merge;

    //static PS::F64 f;
    static PS::F64 m_min;

    PS::S32 getNumberOfFragment() const { return n_frag; }
    PS::S64 getFragmentID() const { return id_frag; }
    PS::S32 getFragmentIDCluster() const { return id_c_frag; }
    PS::F64 getEnergyDissipation() const { return edisp; }
    PS::F64 getHardEnergyDissipation() const { return edisp_d; }

#ifdef MERGE_BINARY
    template <class Tp>
    PS::S32 mergeOutcome(std::vector<Tp> & pfrag) {
        PS::S32 n = collisionOutcome(pfrag);
        flag_merge = 2;
        return n;
    }
#endif
    
    template <class Tpsys>
    void inputPair(Tpsys & pp,
                   std::multimap<PS::S32,PS::S32> & merge_list,
                   std::pair<PS::S32,PS::S32> col_pair);
    template <class Tp>
    PS::S32 collisionOutcome(std::vector<Tp> & pfrag){
        pfrag.clear();
        mass_frag = 0.;

        pos_imp_new = pos_tar_new = (mass_imp*pos_imp + mass_tar*pos_tar)/(mass_imp+mass_tar);
        vel_imp_new = vel_tar_new = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);
        
        n_frag = id_c_frag = 0;
        id_frag = 0;
        //HitAndRun = 0;
        flag_merge = 1;
        return n_frag;
    }
    template <class Tpsys, class Tp>
    void setParticle(Tpsys & pp,
                     std::vector<Tp> & pfrag,
                     std::multimap<PS::S32,PS::S32> & merge_list,
                     PS::S64 & id_next);
    
    template <class Tpsys>
    PS::F64 calcEnergyDissipation(Tpsys & pp,
                                  std::multimap<PS::S32,PS::S32> & merge_list);
    template <class Tpsys>
    void setNeighbors(Tpsys & pp);
    
    void setNewFragmentID( std::map<PS::S64, PS::S64> id_old2new ){
        if ( id_imp < 0 ) id_imp = id_old2new.at(id_imp);
        if ( id_tar < 0 ) id_tar = id_old2new.at(id_tar);
        if ( id_frag < 0 ) id_frag = id_old2new.at(id_frag);
    }

    void write2File(std::ofstream & fp) const {
        PS::F64vec ximp = pos_imp - pos_tar;
        PS::F64vec vimp = vel_imp - vel_tar;
        PS::F64vec dpos_g = pos_g_new - pos_g;
        PS::F64vec dvel_g = vel_g_new - vel_g;
        PS::S32 Flag = 0;
        Flag |= ((PS::S32)flag_merge)<<0;
#ifdef MERGE_BINARY
        Flag |= ((PS::S32)imp_is_binary)<<2;
        Flag |= ((PS::S32)tar_is_binary)<<3;
#endif
        
        fp << std::fixed<<std::setprecision(8)<< this->time << "\t"
           << this->id_imp << "\t" << this->id_tar  << "\t"
           << this->n_frag << "\t" << this->id_frag << "\t"
           << std::scientific<<std::setprecision(15)
           << this->mass_imp  << "\t" << this->mass_tar  << "\t" << this->mass_frag << "\t"
           << sqrt(ximp*ximp) << "\t" << sqrt(vimp*vimp) << "\t"
           << this->col_angle << "\t" << Flag << "\t"
           << sqrt((dpos_g*dpos_g)/(pos_g*pos_g)) << "\t" << sqrt((dvel_g*dvel_g)/(vel_g*vel_g)) 
           << std::endl;
    }

    static void readParameter(std::string name,
                              std::string value){}
    static void showParameter(std::ostream & fout = std::cout) {}
};

PS::F64 Collision0::m_min = 9.426627927538057e-12;


template <class Tpsys>
inline void Collision0::inputPair(Tpsys & pp,
                                  std::multimap<PS::S32,PS::S32> & merge_list,
                                  std::pair<PS::S32,PS::S32> col_pair){
    id_c_imp = col_pair.first;
    id_c_tar = col_pair.second;
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::pair<iterator, iterator> imp_range = merge_list.equal_range(id_c_imp);
    std::pair<iterator, iterator> tar_range = merge_list.equal_range(id_c_tar);
        
    id_imp = pp[id_c_imp].id;
    id_tar = pp[id_c_tar].id;
    //time = pp[id_c_imp].time;
    //assert ( time == pp[id_c_tar].time );
    time = pp[id_c_imp].getTime();
    assert ( time == pp[id_c_tar].getTime() );
        
    pos_imp = pp[id_c_imp].pos;
    pos_tar = pp[id_c_tar].pos;
    vel_imp = pp[id_c_imp].vel;
    vel_tar = pp[id_c_tar].vel;
    PS::F64vec ximp = pos_imp - pos_tar;
    PS::F64vec vimp = vel_imp - vel_tar;
    PS::F64 costheta = std::max(std::min(-(ximp*vimp)/(sqrt(ximp*ximp)*sqrt(vimp*vimp)), 1.), -1.);
    col_angle = std::acos(costheta);
 
    mass_imp = pp[id_c_imp].mass;
    mass_tar = pp[id_c_tar].mass;
    if ( pp[id_c_imp].isMerged ) {
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            mass_imp += pp[id_i].mass;
            assert( time == pp[id_i].getTime() );
        }
    }
    if ( pp[id_c_tar].isMerged ) {
        for (iterator it = tar_range.first; it != tar_range.second; ++it){
            PS::S32 id_j = it->second;
            mass_tar += pp[id_j].mass;
            assert( time == pp[id_j].getTime() );
        }
    }
    r_planet_imp = pp[id_c_imp].r_planet;
    r_planet_tar = pp[id_c_tar].r_planet;
    r_planet_imp_new = -1.;
    r_planet_tar_new = -1.;
    f_imp = f_imp_new = pp[id_c_imp].f;
    f_tar = f_tar_new = pp[id_c_tar].f;
#ifdef MERGE_BINARY
    imp_is_binary = imp_will_be_binary = pp[id_c_imp].isBinary;
    tar_is_binary = tar_will_be_binary = pp[id_c_tar].isBinary;
#endif

    pos_g = (mass_imp*pos_imp + mass_tar*pos_tar)/(mass_imp+mass_tar);
    vel_g = (mass_imp*vel_imp + mass_tar*vel_tar)/(mass_imp+mass_tar);
 
    std::cerr << std::fixed << std::setprecision(8)
              << "Time: " << time << "\tImpactor ID: " << id_imp << "\tTarget ID: " << id_tar << std::endl
              << std::scientific << std::setprecision(15)
              << "Imp_pos_vel:     " << pos_imp.x << "\t" << pos_imp.y << "\t" << pos_imp.z << "\t"
              << vel_imp.x << "\t" << vel_imp.y << "\t" << vel_imp.z << std::endl
              << "Tar_pos_vel:     " << pos_tar.x << "\t" << pos_tar.y << "\t" << pos_tar.z << "\t"
              << vel_tar.x << "\t" << vel_tar.y << "\t" << vel_tar.z << std::endl;
}

template <class Tpsys, class Tp>
inline void Collision0::setParticle(Tpsys & pp,
                                    std::vector<Tp> & pfrag,
                                    std::multimap<PS::S32,PS::S32> & merge_list,
                                    PS::S64 & id_next)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::pair<iterator, iterator> imp_range = merge_list.equal_range(id_c_imp);
    std::pair<iterator, iterator> tar_range = merge_list.equal_range(id_c_tar);

    PS::F64 mass_rem = mass_imp + mass_tar - mass_frag;
    PS::F64vec masspos = 0.;
    PS::F64vec massvel = 0.;
    PS::F64 mass_tot = 0.;
    
    ///////////////////
    /*   Accretion   */
    ///////////////////
    //Mass & Radius
    //if ( HitAndRun ) {
    if ( !flag_merge ) {
        assert ( pp[id_c_tar].r_planet == r_planet_tar );
        if ( r_planet_imp_new < 0. ) {
            pp[id_c_imp].r_planet *= pow((mass_imp - mass_frag)/mass_imp, 1./3.);
            //pp[id_c_imp].setRPlanet(mass_imp - mass_frag);
        } else {
            pp[id_c_imp].r_planet = r_planet_imp_new;
        }
    } else {
        pp[id_c_imp].r_planet = 0.;
        if ( r_planet_tar_new < 0. ) {
            pp[id_c_tar].r_planet *= pow(mass_rem/mass_tar, 1./3.);
            //pp[id_c_tar].setRPlanet(mass_rem);
        } else {
            pp[id_c_tar].r_planet = r_planet_tar_new;
        }
    }
    pp[id_c_imp].mass -= mass_frag;
    
    pp[id_c_imp].f = f_imp_new;
    pp[id_c_tar].f = f_tar_new;
#ifdef MERGE_BINARY
    pp[id_c_imp].isBinary = imp_will_be_binary;
    pp[id_c_tar].isBinary = tar_will_be_binary;
#endif
    
    // ID
    //if ( !HitAndRun ) {
    if ( flag_merge ) {
        pp[id_c_imp].id = pp[id_c_tar].id;
        pp[id_c_imp].isDead = true;
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            pp[id_i].id = pp[id_c_tar].id;
            pp[id_i].isDead = true;
        }
        pp[id_c_tar].isMerged = true;
    }
    
    //Add Fragments
    id_c_frag = n_frag ? pp.size() : 0;
    id_frag   = n_frag ? ( - (128*pp[id_c_imp].id_cluster+id_next+1) ) : 0;
    for ( PS::S32 i=0; i<n_frag; i++ ){

        std::cerr << std::scientific << std::setprecision(15)
                  << "frag_mass: " << pfrag[i].mass << " frag_pos_vel: "
                  << "\t" << pfrag[i].pos.x << "\t" << pfrag[i].pos.y << "\t" << pfrag[i].pos.z
                  << "\t" << pfrag[i].vel.x << "\t" << pfrag[i].vel.y << "\t" << pfrag[i].vel.z << std::endl;

        masspos += pfrag[i].mass * pfrag[i].pos;
        massvel += pfrag[i].mass * pfrag[i].vel;
        mass_tot += pfrag[i].mass;

        if ( pfrag[i].r_planet <= 0 ) pfrag[i].setRPlanet();
        if ( pfrag[i].f == 0 ) pfrag[i].f = pp[id_c_imp].f;

#ifdef USE_INDIVIDUAL_CUTOFF
        pfrag[i].r_out     = pp[id_c_imp].r_out;
        pfrag[i].r_out_inv = pp[id_c_imp].r_out_inv;
        pfrag[i].r_search  = pp[id_c_imp].r_search;
        pfrag[i].v_disp    = pp[id_c_imp].v_disp;
#endif
        pfrag[i].time   = pp[id_c_imp].time;
        pfrag[i].time_c = pp[id_c_imp].time_c;
        pfrag[i].dt     = pp[id_c_imp].dt;
        pfrag[i].phi    = pp[id_c_imp].phi;
        pfrag[i].phi_d  = pp[id_c_imp].phi_d;
        pfrag[i].phi_s  = pp[id_c_imp].phi_s;

        pfrag[i].id         = id_frag - i;
        pfrag[i].id_cluster = pp[id_c_imp].id_cluster;
        pfrag[i].inDomain = true;
        pfrag[i].isDead   = false;
        pfrag[i].isMerged = false;
#ifdef MERGE_BINARY
        pfrag[i].isBinary = false;
#endif
        id_next ++;
        
        pp.push_back(pfrag[i]);
    }
    if ( n_frag ) assert ( (PS::S32)(pp.size()) == id_c_frag + n_frag );

    //Pos & Vel
    pp[id_c_imp].pos = pos_imp_new;
    pp[id_c_imp].vel = vel_imp_new;
    for (iterator it = imp_range.first; it != imp_range.second; ++it){
        PS::S32 id_i = it->second;
        pp[id_i].pos = pos_imp_new;
        pp[id_i].vel = vel_imp_new;
    }    
    pp[id_c_tar].pos = pos_tar_new;
    pp[id_c_tar].vel = vel_tar_new;
    for (iterator it = tar_range.first; it != tar_range.second; ++it){
        PS::S32 id_j = it->second;
        pp[id_j].pos = pos_tar_new;
        pp[id_j].vel = vel_tar_new;
    }
    //if ( !HitAndRun ) {
    if ( flag_merge ) {
            assert( pos_imp_new == pos_tar_new );
            assert( vel_imp_new == vel_tar_new );
    }

    masspos += (mass_imp-mass_frag)*pos_imp_new + mass_tar*pos_tar_new;
    massvel += (mass_imp-mass_frag)*vel_imp_new + mass_tar*vel_tar_new;
    mass_tot += mass_imp-mass_frag + mass_tar;
    pos_g_new = masspos / mass_tot;
    vel_g_new = massvel / mass_tot;

    std::cerr << std::scientific << std::setprecision(15)
              << "Imp_pos_vel_new: " << pos_imp_new.x << "\t" << pos_imp_new.y << "\t" << pos_imp_new.z << "\t"
              << vel_imp_new.x << "\t" << vel_imp_new.y << "\t" << vel_imp_new.z << std::endl
              << "Tar_pos_vel_new: " << pos_tar_new.x << "\t" << pos_tar_new.y << "\t" << pos_tar_new.z << "\t"
              << vel_tar_new.x << "\t" << vel_tar_new.y << "\t" << vel_tar_new.z << std::endl;

}

template <class Tpsys>
inline PS::F64 Collision0::calcEnergyDissipation(Tpsys & pp,
                                                 std::multimap<PS::S32,PS::S32> & merge_list)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::pair<iterator, iterator> imp_range = merge_list.equal_range(id_c_imp);
    std::pair<iterator, iterator> tar_range = merge_list.equal_range(id_c_tar);
    PS::F64 mass_i = pp[id_c_imp].mass + mass_frag;
    
    const PS::F64 eps2  = FP_t::eps2;
    const PS::F64 m_sun = FP_t::m_sun;
    
    ///////////////////////////
    /*   Energy Dissipation  */
    ///////////////////////////
    //Kinetic energy
    PS::F64 e_kin = 0.;
#if 1
    PS::F64vec vel_rel = vel_imp - vel_tar;   
    e_kin -= 0.5 * mass_imp*mass_tar/(mass_imp+mass_tar) * vel_rel*vel_rel;

    vel_rel = vel_imp_new - vel_g;
    e_kin += 0.5 * (mass_imp-mass_frag)* vel_rel*vel_rel;
    vel_rel = vel_tar_new - vel_g;
    e_kin += 0.5 * mass_tar* vel_rel*vel_rel;
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        vel_rel = pp[id_f].vel - vel_g;
        e_kin += 0.5 * pp[id_f].mass * vel_rel*vel_rel;
    }
#else
    e_kin -= mass_imp * vel_imp*vel_imp + mass_tar * vel_tar*vel_tar;

    e_kin += (mass_imp-mass_frag) * vel_imp_new*vel_imp_new + mass_tar * vel_tar_new*vel_tar_new;
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        e_kin += pp[id_f].mass *pp[id_f].vel*pp[id_f].vel;
    }
    e_kin *= 0.5;
#endif
    
    //Interactive energy
    PS::F64 e_int = 0.;
    PS::F64 e_int_d = 0.;
    PS::F64 dphi = 0;
    
    PS::F64vec dr = pos_imp - pos_tar;
    PS::F64 dr2 = dr*dr + eps2;
    PS::F64 rinv = sqrt(1./dr2);
    PS::F64 rinv_new;
    dphi = mass_i * rinv;
    e_int   += pp[id_c_tar].mass * dphi;
    e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
        * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out_inv, pp[id_c_tar].r_out_inv));
#else
        * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
    for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
        PS::S32 id_j = it2->second;
        e_int   += pp[id_j].mass * dphi;
        e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out_inv, pp[id_j].r_out_inv));
#else
            * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
    }
        
    for (iterator it = imp_range.first; it != imp_range.second; ++it){
        PS::S32 id_i = it->second;
        dphi = pp[id_i].mass * rinv;
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
            * (1.-cutoff_W2(dr2, pp[id_i].r_out_inv, pp[id_c_tar].r_out_inv));
#else
            * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
            PS::S32 id_j = it2->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_i].r_out_inv, pp[id_j].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }
    }

    //if ( HitAndRun ){
    if ( !flag_merge ) {
        dr = pos_imp_new - pos_tar_new;
        dr2 = dr*dr + eps2;
        rinv = sqrt(1./dr2);
        dphi = -pp[id_c_imp].mass * rinv;
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out_inv, pp[id_c_tar].r_out_inv));
#else
            * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
            PS::S32 id_j = it2->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out_inv, pp[id_j].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }
        
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            dphi = -pp[id_i].mass * rinv;
            e_int   += pp[id_c_tar].mass * dphi;
            e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_i].r_out_inv, pp[id_c_tar].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
            for (iterator it2 = tar_range.first; it2 != tar_range.second; ++it2){
                PS::S32 id_j = it2->second;
                e_int   += pp[id_j].mass * dphi;
                e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                    * (1.-cutoff_W2(dr2, pp[id_i].r_out_inv, pp[id_j].r_out_inv));
#else
                    * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
            }
        }
    } else {
        assert ( pos_imp_new == pos_tar_new );
    }
    
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f0 = id_c_frag + i;
        dr = pp[id_f0].pos - pp[id_c_imp].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = -pp[id_f0].mass * rinv_new;
        e_int   += pp[id_c_imp].mass * dphi;
        e_int_d += pp[id_c_imp].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out_inv, pp[id_f0].r_out_inv));
#else
            * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            e_int   += pp[id_i].mass * dphi;
            e_int_d += pp[id_i].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_i].r_out_inv, pp[id_f0].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }

        dr = pp[id_f0].pos - pp[id_c_tar].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = -pp[id_f0].mass * rinv_new;
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
            * (1.-cutoff_W2(dr2, pp[id_c_tar].r_out_inv, pp[id_f0].r_out_inv));
#else
            * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        for (iterator it = tar_range.first; it != tar_range.second; ++it){
            PS::S32 id_j = it->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_j].r_out_inv, pp[id_f0].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }      
        for ( PS::S32 j=0; j<i; j++ ){
            PS::S32 id_f1 = id_c_frag + j;
            dr = pp[id_f0].pos - pp[id_f1].pos;
            dr2 = dr*dr + eps2;
            rinv_new = sqrt(1./dr2);
            dphi = -pp[id_f0].mass * pp[id_f1].mass * rinv_new;
            e_int   += dphi;
            e_int_d += dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_f0].r_out_inv, pp[id_f1].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }
    }
    
#if 1
    for ( PS::S32 i=0; i<pp[id_c_imp].neighbor; i++ ){
        PS::S32 id_nei = pp[id_c_imp].n_hard_list.at(i);
        if ( pp[id_nei].id == pp[id_c_imp].id || pp[id_nei].id == pp[id_c_tar].id ) continue;
        PS::F64 m_nei = pp[id_nei].mass;
        
        dr = pos_imp - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv = sqrt(1./dr2);
        dr = pp[id_c_imp].pos - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = m_nei * ( rinv - rinv_new );
        e_int   += m_nei * ( mass_i * rinv - pp[id_c_imp].mass * rinv_new );
        e_int_d += m_nei * ( mass_i * rinv - pp[id_c_imp].mass * rinv_new )
#ifdef USE_INDIVIDUAL_CUTOFF
            * (1.-cutoff_W2(dr2, pp[id_c_imp].r_out_inv, pp[id_nei].r_out_inv));
#else
            * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        for (iterator it = imp_range.first; it != imp_range.second; ++it){
            PS::S32 id_i = it->second;
            e_int   += pp[id_i].mass * dphi;
            e_int_d += pp[id_i].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_i].r_out_inv, pp[id_nei].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }
        for ( PS::S32 j=0; j<n_frag; j++ ){
            PS::S32 id_f = id_c_frag + j;
            dr = pp[id_f].pos - pp[id_nei].pos;
            dr2 = dr*dr + eps2;
            rinv_new = sqrt(1./dr2);
            dphi = -m_nei * pp[id_f].mass * rinv_new;
            e_int   += dphi;
            e_int_d += dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_f].r_out_inv, pp[id_nei].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }
    }
    for ( PS::S32 i=0; i<pp[id_c_tar].neighbor; i++ ){
        PS::S32 id_nei = pp[id_c_tar].n_hard_list.at(i);
        if ( pp[id_nei].id == pp[id_c_imp].id || pp[id_nei].id == pp[id_c_tar].id ) continue;
        PS::F64 m_nei = pp[id_nei].mass;
        
        dr = pos_tar - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv = sqrt(1./dr2);
        dr = pp[id_c_tar].pos - pp[id_nei].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        dphi = m_nei * ( rinv - rinv_new );
        e_int   += pp[id_c_tar].mass * dphi;
        e_int_d += pp[id_c_tar].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
            * (1.-cutoff_W2(dr2, pp[id_c_tar].r_out_inv, pp[id_nei].r_out_inv));
#else
            * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        for (iterator it = tar_range.first; it != tar_range.second; ++it){
            PS::S32 id_j = it->second;
            e_int   += pp[id_j].mass * dphi;
            e_int_d += pp[id_j].mass * dphi
#ifdef USE_INDIVIDUAL_CUTOFF
                * (1.-cutoff_W2(dr2, pp[id_j].r_out_inv, pp[id_nei].r_out_inv));
#else
                * (1.-cutoff_W2(dr2, FP_t::r_out_inv));
#endif
        }
    }
#endif
    
    // Gravitational energy
    PS::F64 e_sun = 0.;
    e_sun +=  mass_imp/sqrt(pos_imp*pos_imp+eps2)
        + mass_tar/sqrt(pos_tar*pos_tar+eps2);

    dr = pp[id_c_tar].pos;
    dr2 = dr*dr + eps2;
    rinv_new = sqrt(1./dr2);
    e_sun -= mass_tar * rinv_new;
    
    dr = pp[id_c_imp].pos;
    dr2 = dr*dr + eps2;
    rinv_new = sqrt(1./dr2);
    e_sun -= (mass_imp-mass_frag) * rinv_new;

    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        dr = pp[id_f].pos;
        dr2 = dr*dr + eps2;
        rinv_new = sqrt(1./dr2);
        e_sun -= pp[id_f].mass * rinv_new;
    }
    e_sun *= m_sun;
    
    edisp = e_kin + e_int + e_sun;
    edisp_d = e_kin + e_int_d + e_sun;
    //PRC(e_kin); PRC(e_int);PRC(e_int_d); PRC(e_sun); PRL(edisp);
    std::cerr << std::scientific << std::setprecision(15)
              << "DEkin: " << e_kin << "\tDEint: " << e_int << "\tDEint_hard: " << e_int_d
              << "\tDEsun: " << e_sun << "\tDEtot: " << edisp << std::endl;

    return edisp;
}

template <class Tpsys>
inline void Collision0::setNeighbors(Tpsys & pp)
{   
    ///////////////////////
    /*   Set Neighbors   */
    ///////////////////////
    for ( PS::S32 i=0; i<n_frag; i++ ){
        PS::S32 id_f = id_c_frag + i;
        pp[id_f].neighbor = pp[id_c_imp].neighbor + n_frag;
        pp[id_f].n_hard_list.clear();
        pp[id_f].n_hard_list.push_back(id_c_imp);
        for ( PS::S32 j=0; j<pp[id_c_imp].neighbor; j++ ){
            PS::S32 id_nei = pp[id_c_imp].n_hard_list.at(j);
            pp[id_f].n_hard_list.push_back(id_nei);
            pp[id_nei].n_hard_list.push_back(id_f);
            pp[id_nei].neighbor ++;
        }
        for ( PS::S32 j=0; j<n_frag; j++ ){
            if ( i != j ) pp[id_f].n_hard_list.push_back(id_c_frag + j);
        }
        assert ( pp[id_f].neighbor == (PS::S32)(pp[id_f].n_hard_list.size()) );
    }

    pp[id_c_imp].neighbor += n_frag;
    for ( PS::S32 i=0; i<n_frag; i++ ){
        pp[id_c_imp].n_hard_list.push_back(id_c_frag + i);
    }
    assert ( pp[id_c_imp].neighbor == (PS::S32)(pp[id_c_imp].n_hard_list.size()) );
    for ( PS::S32 i=0; i<pp[id_c_imp].neighbor; i++ ){
        PS::S32 id_nei = pp[id_c_imp].n_hard_list.at(i);
        assert ( pp[id_nei].neighbor == (PS::S32)(pp[id_nei].n_hard_list.size()) );
    }
}

template <class Tp>
void setFragmentCircle(std::vector<Tp> & pfrag,
                       PS::F64vec & masspos,
                       PS::F64vec & massvel,
                       std::vector<PS::F64> mass,
                       PS::F64vec pos_g,
                       PS::F64vec vel_g,
                       PS::F64 x_frag,
                       PS::F64 v_frag,
                       PS::F64vec vec1,
                       PS::F64vec vec2)
{
    PS::S32 n_frag = mass.size();
    PS::F64vec e0 = vec1;
    PS::F64vec e1 = vec2 - ((vec2*vec1)/(vec1*vec1))*vec1;
    e0 = e0 / sqrt(e0*e0);
    e1 = e1 / sqrt(e1*e1);

    pfrag.clear();
    PS::F64 theta0 = 2. * M_PI * drand48();
    for ( PS::S32 i=0; i<n_frag; i++ ){
        Tp new_frag;
        new_frag.mass = mass[i];

        PS::F64 theta = theta0 + 2.*M_PI*i/n_frag;
        new_frag.pos = pos_g + x_frag * ( cos(theta)*e0 + sin(theta)*e1 );
        new_frag.vel = vel_g + v_frag * ( cos(theta)*e0 + sin(theta)*e1 );

        masspos += new_frag.mass * new_frag.pos;
        massvel += new_frag.mass * new_frag.vel;
        
        pfrag.push_back(new_frag);
    }

}
