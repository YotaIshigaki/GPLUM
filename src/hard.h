#pragma once

class ParticleCluster{
public:
    static std::vector<FPH_t> ptcl_hard;   
    std::vector<PS::S32> id_list;

    static void resize_hard(PS::S32 n) { ptcl_hard.resize(n); }
    static void reserve_hard(PS::S32 n) { ptcl_hard.reserve(n); }
    static void clear_hard() { ptcl_hard.clear(); }
    
    FPH_t & operator [](PS::S32 i) { return ptcl_hard[id_list[i]]; }
    FPH_t const & operator [](PS::S32 i) const { return ptcl_hard[id_list[i]]; }
    FPH_t & at(PS::S32 i) { return ptcl_hard.at(id_list.at(i)); }
    FPH_t const & at(PS::S32 i) const { return ptcl_hard.at(id_list.at(i)); }

    void clear() { id_list.clear(); }
    
    void push_back(const FPH_t & pi) {
        PS::S32 i = 0;
#pragma omp critical
        {
            i = ptcl_hard.size();
            ptcl_hard.push_back(pi);
        }
        id_list.push_back(i);
    }
    void push_back(FPH_t && pi) {
        PS::S32 i = 0;
#pragma omp critical
        {
            i = ptcl_hard.size();
            ptcl_hard.push_back(pi);
        }
        id_list.push_back(i);
    }
    void resize(const PS::U32 n) {
        PS::S32 i = 0;
        PS::S32 m = n - id_list.size();
        if ( m > 0 ) {
#pragma omp critical
            {
                i = ptcl_hard.size();
                ptcl_hard.resize(i + m);
            }
        }
        PS::S32 j = id_list.size();
        id_list.resize(n);
        for (PS::U32 k=j; k<n; k++) id_list.at(k) = i + k;
    }
    void reserve (const PS::U32 n) { id_list.reserve(n); }

    PS::U32 size() const { return id_list.size(); }
    PS::U32 capacity() const { return id_list.capacity(); }
};

std::vector<FPH_t> ParticleCluster::ptcl_hard;

class ClusterSystem : public std::vector<ParticleCluster> {
public:
    void clear() { 
        std::vector<ParticleCluster>::clear();
        ParticleCluster::clear_hard();
    }

    void reserve_hard(const PS::U32 n) { ParticleCluster::reserve_hard(n); }
};

class HardSystem{
public:
    std::vector<PS::S32> list_iso;
    std::vector<std::vector<std::pair<bool,PS::S32> > > list_multi;
    //std::vector<std::vector<FPH_t> > ptcl_multi;
    ClusterSystem ptcl_multi;
    std::map<PS::S32,PS::S32> mp_cluster;  // cluster ID -> cluster adress in HardSystem

    std::vector<Collision> collision_list;
    std::vector<std::pair<PS::S32,PS::S32> > frag_list;

    PS::S32 n_col;
    PS::S32 n_frag;
    PS::F64 edisp;
    PS::F64 edisp_d;

    //static PS::F64 f;

    std::vector<PS::S32>   n_col_list;
    std::vector<PS::S32>   n_frag_list;
    std::vector<Collision> col_list_tot;
    std::vector<PS::S64>   id_frag_list;
    std::vector<PS::S64>   id_frag_loc;
    std::vector<PS::S32>   col_recv;
    std::vector<PS::S32>   frag_recv;
    

    PS::S32 getNumberOfClusterLocal() const { return ptcl_multi.size(); }
    PS::S32 getNumberOfClusterGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_cluster_loc = ptcl_multi.size();
        return PS::Comm::getSum(n_cluster_loc);
#else
        return ptcl_multi.size();
#endif
    }
    PS::S32 getNumberOfIsolatedParticleLocal() const { return list_iso.size(); }
    PS::S32 getNumberOfIsolatedParticleGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_ptcl_loc = list_iso.size();
        return PS::Comm::getSum(n_ptcl_loc);
#else
        return list_iso.size();
#endif
    }
    PS::S32 getNumberOfParticleLocal(){
        PS::S32 n = 0;
        PS::S32 size = ptcl_multi.size();
        for ( PS::S32 i=0; i<size; i++ ) n += ptcl_multi[i].size();
        n += list_iso.size();
        return n;
    }
    PS::S32 getNumberOfParticleGlobal() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_ptcl_loc = getNumberOfParticleLocal();
        return PS::Comm::getSum(n_ptcl_loc);
#else
        return getNumberOfParticleLocal();
#endif
    }
    PS::S32 getNumberOfParticleInLargestClusterLocal(){
        PS::S32 n = 0;
        PS::S32 size = ptcl_multi.size();
        for ( PS::S32 i=0; i<size; i++ )
            if ( n < size ) n = ptcl_multi[i].size();
        return n;
    }
    PS::S32 getNumberOfParticleInLargestClusterGlobal() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_loc = getNumberOfParticleInLargestClusterLocal();
        return PS::Comm::getMaxValue(n_loc);
#else
        return getNumberOfParticleInLargestClusterLocal();
#endif
    }
    
    PS::S32 getNumberOfFragmentLocal() const { return n_frag; }
    PS::S32 getNumberOfFragmentGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_frag_ = n_frag;
        return PS::Comm::getSum(n_frag_);
#else
        return n_frag;
#endif
    }
    PS::S32 getNumberOfCollisionLocal() const { return n_col; }
    PS::S32 getNumberOfCollisionGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::S32 n_col_ = n_col;
        return PS::Comm::getSum(n_col_);
#else
        return n_col;
#endif
    }
    PS::F64 getEnergyDissipationLocal() const { return edisp; }
    PS::F64 getEnergyDissipationGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::F64 edisp_ = edisp;
        return PS::Comm::getSum(edisp_);
#else
        return edisp;
#endif
    }
    PS::F64 getHardEnergyDissipationLocal() const { return edisp_d; }
    PS::F64 getHardEnergyDissipationGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        PS::F64 edisp_d_ = edisp_d;
        return PS::Comm::getSum(edisp_d_);
#else
        return edisp_d;
#endif
    }

    void reserve_hard(const PS::U32 n) { ptcl_multi.reserve_hard(n); }
    
    void showParticleID() const {
        PS::S32 size = ptcl_multi.size();
        for ( PS::S32 i=0; i<size; i++ ){
            bool flag = ptcl_multi.at(i).at(0).inDomain;
            std::cout << "Rank " << PS::Comm::getRank()
                      << " Cluster " << ptcl_multi.at(i).at(0).id_cluster
                      << " (" << flag << "): ";
            PS::S32 sizei = ptcl_multi.at(i).size();
            for ( PS::S32 j=0; j<sizei; j++ ){
                std::cout << " " << ptcl_multi.at(i).at(j).id;
                assert( flag == ptcl_multi.at(i).at(j).inDomain );
            }
            std::cout << std::endl;
        }
    }

    void clear(){
        list_iso.clear();
        list_multi.clear();        
        ptcl_multi.clear();
        mp_cluster.clear();
        collision_list.clear();
        frag_list.clear();

        n_col = n_frag = 0;
        edisp = edisp_d = 0.;
    }

    template <class Tpsys, class Tpsys2>
    PS::S32 makeList(Tpsys & pp,
                     Tpsys2 & ex_pp);

    template <class Tpsys, class Tpsys2, class NL, class NL2>
    PS::S32 timeIntegrate(Tpsys & pp,
                          Tpsys2 & ex_pp,
                          NL & NList,
                          NL2 & ex_NList,
                          const PS::S32 istep);

    static void rewriteFragmentID(std::vector<PS::S64> & id_frag_list,
                                  std::vector<Collision> & col_list,
                                  PS::S32 n_col_tot,
                                  PS::S32 n_frag_tot,
                                  PS::S64 & id_next);
    template <class Tpsys, class Tpsys2>
    PS::S32 addFragment2ParticleSystem(Tpsys & pp,
                                       Tpsys2 & ex_pp,
                                       PS::S64 & id_next,
                                       std::ofstream & fp); 
};


template <class Tpsys, class Tpsys2>
inline PS::S32 HardSystem::makeList(Tpsys & pp,
                                    Tpsys2 & ex_pp)
{
    const PS::S32 n_pp    = pp.getNumberOfParticleLocal();
    const PS::S32 n_ex_pp = ex_pp.getNumberOfParticleLocal();
    PS::S32 n_ptcl_loc = 0;
    PS::S32 n_cluster_loc = 0;
    PS::U32 tmp = 0;
        
    for ( PS::S32 i=0; i<n_pp; i++ ){
        if ( pp[i].neighbor ){
            if ( pp[i].isSent ) continue;
            auto itr = mp_cluster.find(pp[i].id_cluster);                    
            if ( itr == mp_cluster.end() ){
                mp_cluster[pp[i].id_cluster] = tmp;
                list_multi.push_back(std::vector<std::pair<bool,PS::S32> >{std::make_pair(true, i)});
                tmp ++;
                assert( list_multi.size() == tmp );
            } else {
                assert( mp_cluster.at(pp[i].id_cluster) == itr->second );
                list_multi.at(itr->second).push_back(std::make_pair(true, i));
            }
            n_ptcl_loc ++;
        } else {
            list_iso.push_back(i);
            n_ptcl_loc ++;
        }
    }
    
    for ( PS::S32 i=0; i<n_ex_pp; i++ ){
        assert( ex_pp[i].neighbor > 0 );
        assert( ex_pp[i].isSent );
        assert( !ex_pp[i].inDomain );
        auto itr = mp_cluster.find(ex_pp[i].id_cluster);                    
        if ( itr == mp_cluster.end() ){
            mp_cluster[ex_pp[i].id_cluster] = tmp;
            list_multi.push_back(std::vector<std::pair<bool,PS::S32> >{std::make_pair(false, i)});
            tmp ++;
            assert( list_multi.size() == tmp );
        } else {
            assert( mp_cluster.at(ex_pp[i].id_cluster) == itr->second );
            list_multi.at(itr->second).push_back(std::make_pair(false, i));
        }
        n_ptcl_loc ++;
    }
    n_cluster_loc = tmp;
    assert( list_multi.size() == tmp );
    ptcl_multi.resize(n_cluster_loc);

    return n_ptcl_loc;
}

#if 0
template <class Tpsys, class Tpsys2, class NL, class NL2>
inline PS::S32 HardSystem::timeIntegrate(Tpsys & pp,
                                         Tpsys2 & ex_pp,
                                         NL & NList,
                                         NL2 & ex_NList,
                                         const PS::S32 istep)
{
    PS::S32 n_all = list_multi.size() + list_iso.size();
    PS::S32 n_ptcl_loc = 0;
    
#pragma omp parallel for reduction(+:n_ptcl_loc) schedule (dynamic)
    for ( PS::S32 ii=0; ii<n_all; ii++ ){
        PS::S32 size = list_multi.size();
        if ( ii<size ){
            PS::S32 i = ii;
            //for ( PS::S32 i=0; i<size; i++ ){
            PS::S32 n_p = list_multi.at(i).size();
            PS::S32 id_cluster = 0;
            std::map<PS::S64, PS::S32> id_map;
            ptcl_multi[i].clear();
            id_map.clear();
            
            // Add Particle To Hard System
            ptcl_multi[i].reserve(n_p);
            for ( PS::S32 j=0; j<n_p; j++ ){
                std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
                PS::S32 id_loc = adr.second;
                if ( adr.first ) {
                    ptcl_multi[i].push_back(FPH_t(pp[id_loc]));
                    ptcl_multi[i][j].copyList(NList[id_loc]);
                } else {
                    ptcl_multi[i].push_back(FPH_t(ex_pp[id_loc]));
                    ptcl_multi[i][j].copyList(ex_NList[id_loc]);
                }
                if ( j==0 ) id_cluster = ptcl_multi[i][j].id_cluster;
                assert ( ptcl_multi[i][j].id_cluster == id_cluster );
                id_map[ptcl_multi[i][j].id] = j;
            }

            // Make Neighbor List
#ifdef TEST_PTCL
            for ( PS::S32 j=0; j<n_p; j++ ) ptcl_multi[i][j].makeHardList(id_map, ptcl_multi[i]);
#else
            for ( PS::S32 j=0; j<n_p; j++ ) ptcl_multi[i][j].makeHardList(id_map);
#endif

            PS::S32 n_col_tmp  = 0;
            PS::S32 n_frag_tmp = 0;
            PS::F64 edisp_tmp   = 0.;
            PS::F64 edisp_d_tmp = 0.;
            timeIntegrate_multi(ptcl_multi[i], 0., FP_t::dt_tree,
                                n_col_tmp, n_frag_tmp, edisp_tmp, edisp_d_tmp, collision_list_tmp);
            if ( n_col_tmp > 0 ){
#pragma omp critical
                {
                    n_col  += n_col_tmp;
                    n_frag += n_frag_tmp;
                    edisp   += edisp_tmp;
                    edisp_d += edisp_d_tmp;
                    for ( PS::S32 j=0; j<n_col_tmp; j++ ) 
                        collision_list.push_back(collision_list_tmp.at(j));
                    for ( PS::S32 j=0; j<n_frag_tmp; j++ ) 
                        frag_list.push_back(std::make_pair(i, ptcl_multi[i].size()-n_frag_tmp+j));
                }
            }

            PS::S32 sizei = ptcl_multi[i].size();
            for ( PS::S32 j=0; j<sizei; j++ ) {
                ptcl_multi[i][j].n_cluster = ptcl_multi[i].size();
                ptcl_multi[i][j].resetTime();
            }
                
            for ( PS::S32 j=0; j<sizei-n_frag_tmp; j++ ){
                std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
                PS::S32 id_loc = adr.second;
                if ( adr.first ) {
                    if ( !ptcl_multi[i][j].isDead ) assert ( pp[id_loc].id == ptcl_multi[i][j].id );
                    pp[id_loc] = FP_t(ptcl_multi[i][j]);
                    pp[id_loc].neighbor = ptcl_multi[i][j].n_list.size();
                } else {
                    if ( !ptcl_multi[i][j].isDead ) assert ( ex_pp[id_loc].id == ptcl_multi[i][j].id );
                    ex_pp[id_loc] = FP_t(ptcl_multi[i][j]);
                    ex_pp[id_loc].neighbor = ptcl_multi[i][j].n_list.size();
                }
                n_ptcl_loc ++;
            }
            
            
        } else {
            PS::S32 i = ii - list_multi.size();
            //#pragma omp parallel for reduction(+:n_ptcl_loc)
            //for ( PS::S32 i=0; i<list_iso.size(); i++ ){
#ifndef WITHOUT_SUN
            if ( pp[list_iso[i]].getEccentricity() < 0.8 && FP_t::eps2 == 0. ){
                timeIntegrateKepler_isolated(pp[list_iso[i]], (istep-1)*FP_t::dt_tree, istep*FP_t::dt_tree);
            } else {
                FPH_t pi = FPH_t(pp[list_iso[i]]);
                timeIntegrate_isolated(pi, 0., FP_t::dt_tree);
                pi.resetTime();
                pp[list_iso[i]] = FP_t(pi);
            }
#else
            freeMotion(pp[list_iso[i]], (istep-1)*FP_t::dt_tree, istep*FP_t::dt_tree);
#endif
            n_ptcl_loc ++;

            pp[list_iso[i]].n_cluster = 1;
        }
    }
    
    return n_ptcl_loc;
}
#else
template <class Tpsys, class Tpsys2, class NL, class NL2>
inline PS::S32 HardSystem::timeIntegrate(Tpsys & pp,
                                         Tpsys2 & ex_pp,
                                         NL & NList,
                                         NL2 & ex_NList,
                                         const PS::S32 istep)
{
    PS::S32 n_ptcl_loc = 0;
    const PS::S32 large_cluster_size = 64;

    PS::S32 size = list_multi.size();
    std::vector<PS::S32> small_cluster, large_cluster;
    small_cluster.clear();
    large_cluster.clear();
    for ( PS::S32 i=0; i<size; i++ ){
        if ( list_multi.at(i).size() < large_cluster_size ){
            small_cluster.push_back(i);
        } else {
            large_cluster.push_back(i);
        }
    }
    PS::S32 n_small = small_cluster.size();
    PS::S32 n_large = large_cluster.size();
    PS::S32 n_all = n_small + list_iso.size();

    ///////////////////////////////////
    ///   Integrate Large Cluster   ///
    ///////////////////////////////////
    for ( PS::S32 ii=0; ii<n_large; ii++ ){
        PS::S32 i = large_cluster.at(ii);
        PS::S32 n_p = list_multi.at(i).size();
        PS::S32 id_cluster = 0;
        std::map<PS::S64, PS::S32> id_map;
        ptcl_multi[i].clear();
        id_map.clear();
        assert( large_cluster_size <= n_p );
        
        // Add Particle To Hard System
        ptcl_multi[i].reserve(n_p);
        for ( PS::S32 j=0; j<n_p; j++ ){
            std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
            PS::S32 id_loc = adr.second;
            if ( adr.first ) {
                ptcl_multi[i].push_back(FPH_t(pp[id_loc]));
                ptcl_multi[i][j].copyList(NList[id_loc]);
            } else {
                ptcl_multi[i].push_back(FPH_t(ex_pp[id_loc]));
                ptcl_multi[i][j].copyList(ex_NList[id_loc]);
            }
            if ( j==0 ) id_cluster = ptcl_multi[i][j].id_cluster;
            assert ( ptcl_multi[i][j].id_cluster == id_cluster );
            id_map[ptcl_multi[i][j].id] = j;
        }

        // Make Neighbor List
#ifdef TEST_PTCL
#pragma omp parallel for
        for ( PS::S32 j=0; j<n_p; j++ ) ptcl_multi[i][j].makeHardList(id_map, ptcl_multi[i]);
#else
#pragma omp parallel for
        for ( PS::S32 j=0; j<n_p; j++ ) ptcl_multi[i][j].makeHardList(id_map);
#endif
        
        PS::S32 n_col_tmp  = 0;
        PS::S32 n_frag_tmp = 0;
        PS::F64 edisp_tmp   = 0.;
        PS::F64 edisp_d_tmp = 0.;
        std::vector<Collision> collision_list_tmp;
        timeIntegrate_multi_omp(ptcl_multi[i], 0., FP_t::dt_tree,
                                n_col_tmp, n_frag_tmp, edisp_tmp, edisp_d_tmp, collision_list_tmp);
        if ( n_col_tmp > 0 ){
            n_col  += n_col_tmp;
            n_frag += n_frag_tmp;
            edisp   += edisp_tmp;
            edisp_d += edisp_d_tmp;
            for ( PS::S32 j=0; j<n_col_tmp; j++ )
                collision_list.push_back(collision_list_tmp.at(j));
            for ( PS::S32 j=0; j<n_frag_tmp; j++ ) 
                frag_list.push_back(std::make_pair(i, ptcl_multi[i].size()-n_frag_tmp+j));
        }
        
        PS::S32 sizei = ptcl_multi[i].size();
#pragma omp parallel for
        for ( PS::S32 j=0; j<sizei; j++ ) {
            ptcl_multi[i][j].n_cluster = ptcl_multi[i].size();
            ptcl_multi[i][j].resetTime();
        }

#pragma omp parallel for reduction(+:n_ptcl_loc)
        for ( PS::S32 j=0; j<sizei-n_frag_tmp; j++ ){
            std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
            PS::S32 id_loc = adr.second;
            if ( adr.first ) {
                if ( !ptcl_multi[i][j].isDead ) assert ( pp[id_loc].id == ptcl_multi[i][j].id );
                pp[id_loc] = FP_t(ptcl_multi[i][j]);
                pp[id_loc].neighbor = ptcl_multi[i][j].n_list.size();
            } else {
                if ( !ptcl_multi[i][j].isDead ) assert ( ex_pp[id_loc].id == ptcl_multi[i][j].id );
                ex_pp[id_loc] = FP_t(ptcl_multi[i][j]);
                ex_pp[id_loc].neighbor = ptcl_multi[i][j].n_list.size();
            }
            n_ptcl_loc ++;
        }
    }

    //////////////////////////////////////////////////////
    ///   Integrate Small Cluster & Isolated Particle  ///
    //////////////////////////////////////////////////////
#pragma omp parallel for reduction(+:n_ptcl_loc) schedule (dynamic)
    for ( PS::S32 ii=0; ii<n_all; ii++ ){
        if ( ii<n_small ){
            PS::S32 i = small_cluster.at(ii);
            PS::S32 n_p = list_multi.at(i).size();
            PS::S32 id_cluster = 0;
            std::map<PS::S64, PS::S32> id_map;
            ptcl_multi[i].clear();
            id_map.clear();
            assert( large_cluster_size > n_p );
            
            // Add Particle To Hard System
            ptcl_multi[i].reserve(n_p);
            for ( PS::S32 j=0; j<n_p; j++ ){
                std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
                PS::S32 id_loc = adr.second;
                if ( adr.first ) {
                    ptcl_multi[i].push_back(FPH_t(pp[id_loc]));
                    ptcl_multi[i][j].copyList(NList[id_loc]);
                } else {
                    ptcl_multi[i].push_back(FPH_t(ex_pp[id_loc]));
                    ptcl_multi[i][j].copyList(ex_NList[id_loc]);
                }
                if ( j==0 ) id_cluster = ptcl_multi[i][j].id_cluster;
                assert ( ptcl_multi[i][j].id_cluster == id_cluster );
                id_map[ptcl_multi[i][j].id] = j;
            }

            // Make Neighbor List
#ifdef TEST_PTCL
            for ( PS::S32 j=0; j<n_p; j++ ) ptcl_multi[i][j].makeHardList(id_map, ptcl_multi[i]);
#else
            for ( PS::S32 j=0; j<n_p; j++ ) ptcl_multi[i][j].makeHardList(id_map);
#endif

            PS::S32 n_col_tmp  = 0;
            PS::S32 n_frag_tmp = 0;
            PS::F64 edisp_tmp   = 0.;
            PS::F64 edisp_d_tmp = 0.;
            std::vector<Collision> collision_list_tmp;
            //timeIntegrate_multi(ptcl_multi[i], (istep-1)*FP_t::dt_tree, istep*FP_t::dt_tree, f,
            //                    n_col_tmp, n_frag_tmp, edisp_tmp, edisp_d_tmp, collision_list_tmp);
            timeIntegrate_multi(ptcl_multi[i], 0., FP_t::dt_tree,
                                n_col_tmp, n_frag_tmp, edisp_tmp, edisp_d_tmp, collision_list_tmp);
            if ( n_col_tmp > 0 ){
#pragma omp critical
                {
                    n_col  += n_col_tmp;
                    n_frag += n_frag_tmp;
                    edisp   += edisp_tmp;
                    edisp_d += edisp_d_tmp;
                    for ( PS::S32 j=0; j<n_col_tmp; j++ )
                        collision_list.push_back(collision_list_tmp.at(j));
                    for ( PS::S32 j=0; j<n_frag_tmp; j++ )
                        frag_list.push_back(std::make_pair(i, ptcl_multi[i].size()-n_frag_tmp+j));
                }
            }

            PS::S32 sizei = ptcl_multi[i].size();
            for ( PS::S32 j=0; j<sizei; j++ ) {
                ptcl_multi[i][j].n_cluster = ptcl_multi[i].size();
                ptcl_multi[i][j].resetTime();
            }
                
            for ( PS::S32 j=0; j<sizei-n_frag_tmp; j++ ){
                std::pair<bool,PS::S32> adr = list_multi.at(i).at(j);
                PS::S32 id_loc = adr.second;
                if ( adr.first ) {
                    if ( !ptcl_multi[i][j].isDead ) assert ( pp[id_loc].id == ptcl_multi[i][j].id );
                    pp[id_loc] = FP_t(ptcl_multi[i][j]);
                    pp[id_loc].neighbor = ptcl_multi[i][j].n_list.size();
                } else {
                    if ( !ptcl_multi[i][j].isDead ) assert ( ex_pp[id_loc].id == ptcl_multi[i][j].id );
                    ex_pp[id_loc] = FP_t(ptcl_multi[i][j]);
                    ex_pp[id_loc].neighbor = ptcl_multi[i][j].n_list.size();
                }
                n_ptcl_loc ++;
            }
            
            
        } else {
            PS::S32 i = ii - n_small;
#ifndef WITHOUT_SUN
            if ( pp[list_iso[i]].getEccentricity() < 0.8 && FP_t::eps2 == 0. ){
                timeIntegrateKepler_isolated(pp[list_iso[i]], (istep-1)*FP_t::dt_tree, istep*FP_t::dt_tree);
            } else {
                FPH_t pi = FPH_t(pp[list_iso[i]]);
                timeIntegrate_isolated(pi, 0., FP_t::dt_tree);
                //timeIntegrate_isolated(pi, (istep-1)*FP_t::dt_tree, istep*FP_t::dt_tree);
                pi.resetTime();
                pp[list_iso[i]] = FP_t(pi);
            }
#else
            freeMotion(pp[list_iso[i]], (istep-1)*FP_t::dt_tree, istep*FP_t::dt_tree);
#endif
            n_ptcl_loc ++;

            pp[list_iso[i]].n_cluster = 1;
        }
    }
    
    return n_ptcl_loc;
}
#endif

inline void HardSystem::rewriteFragmentID(std::vector<PS::S64> & id_frag_list,
                                          std::vector<Collision> & col_list,
                                          PS::S32 n_col_tot,
                                          PS::S32 n_frag_tot,
                                          PS::S64 & id_next)
{
    std::map<PS::S64, PS::S64> id_old2new;
    id_old2new.clear();
        
    for ( PS::S32 i=0; i<n_col_tot; i++ ){
        PS::S32 n_fragi  = col_list[i].getNumberOfFragment();
        PS::S64 id_fragi = col_list[i].getFragmentID();
        if ( n_fragi == 0 ) continue;
        for ( PS::S32 j=0; j<n_fragi; j++ ){
            id_old2new[id_fragi - j] = id_next;
            id_next ++;
        }
    }
    
    for ( PS::S32 i=0; i<n_frag_tot; i++ )
        if ( id_frag_list[i] < 0 ) id_frag_list[i] = id_old2new.at(id_frag_list[i]);
    for ( PS::S32 i=0; i<n_col_tot; i++ )
        col_list[i].setNewFragmentID(id_old2new);
}

template <class Tpsys, class Tpsys2>
inline PS::S32 HardSystem::addFragment2ParticleSystem(Tpsys & pp,
                                                      Tpsys2 & ex_pp,
                                                      PS::S64 & id_next,
                                                      std::ofstream & fp)
{
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    //PS::S32 * n_col_list     = nullptr;
    //PS::S32 * n_frag_list    = nullptr;
    //Collision * col_list_tot = nullptr;
    //PS::S64 * id_frag_list   = nullptr;
    //PS::S64 * id_frag_loc    = nullptr;
    //PS::S32 * col_recv  = nullptr;
    //PS::S32 * frag_recv = nullptr;

    //id_frag_loc = new PS::S64[n_frag];
    id_frag_loc.resize(n_frag);
    for ( PS::S32 i=0; i<n_frag; i++ ){
        std::pair<PS::S32, PS::S32> adr = frag_list.at(i);
        id_frag_loc[i] = ptcl_multi[adr.first][adr.second].id;
    }
    if ( PS::Comm::getRank() == 0 ){
        //n_col_list  = new PS::S32[n_proc];
        //n_frag_list = new PS::S32[n_proc];
        //col_recv    = new PS::S32[n_proc];
        //frag_recv   = new PS::S32[n_proc];
        n_col_list.resize(n_proc);
        n_frag_list.resize(n_proc);
        col_recv.resize(n_proc);
        frag_recv.resize(n_proc);
    }
    // Send Number of Collision & Fragments
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Gather(&n_col,      1, PS::GetDataType(n_col),
               &n_col_list[0],  1, PS::GetDataType(n_col_list[0]),  0, MPI_COMM_WORLD);
    MPI_Gather(&n_frag,     1, PS::GetDataType(n_frag),
               &n_frag_list[0], 1, PS::GetDataType(n_frag_list[0]), 0, MPI_COMM_WORLD);
#else
    n_col_list[0]  = n_col;
    n_frag_list[0] = n_frag;
#endif
    //PS::Comm::gather(&n_col, 1, n_col_list);
    //PS::Comm::gather(&n_frag, 1, n_frag_list);

    PS::S32 n_col_tot = 0;
    PS::S32 n_frag_tot = 0;
    if ( PS::Comm::getRank() == 0 ){
        PS::S32 tmp_col = 0;
        PS::S32 tmp_frag = 0;
        for ( PS::S32 i=0; i<n_proc; i++ ){
            col_recv[i]  = tmp_col;
            frag_recv[i] = tmp_frag;
            tmp_col  += n_col_list[i];
            tmp_frag += n_frag_list[i];
        }
        //col_recv[n_proc]  = tmp_col;
        //frag_recv[n_proc] = tmp_frag;
        n_col_tot = tmp_col;
        n_frag_tot = tmp_frag;
        //col_list_tot = new Collision[n_col_tot];
        //id_frag_list = new PS::S64[n_frag_tot];
        col_list_tot.resize(n_col_tot);
        id_frag_list.resize(n_frag_tot);
    }
    // Send Collision Information & Fragments ID
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Gatherv(&collision_list[0],  n_col,                         PS::GetDataType(collision_list[0]),
                &col_list_tot[0],   &n_col_list[0],  &col_recv[0],  PS::GetDataType(col_list_tot[0]),   0, MPI_COMM_WORLD);
    MPI_Gatherv(&id_frag_loc[0],     n_frag,                        PS::GetDataType(id_frag_loc[0]),
                &id_frag_list[0],   &n_frag_list[0], &frag_recv[0], PS::GetDataType(id_frag_list[0]),   0, MPI_COMM_WORLD);
#else
    for(PS::S32 i=0; i<n_col; i++)  col_list_tot[i] = collision_list[i];
    for(PS::S32 i=0; i<n_frag; i++) id_frag_list[i] = id_frag_loc[i];
#endif
    //PS::Comm::gatherV(&collision_list[0], n_col,  col_list_tot, n_col_list,  col_recv);
    //PS::Comm::gatherV(id_frag_loc,        n_frag, id_frag_list, n_frag_list, frag_recv);
    
    // Rewrite Fragments ID
    if ( PS::Comm::getRank() == 0 ){
        rewriteFragmentID(id_frag_list, col_list_tot, n_col_tot, n_frag_tot, id_next);
    }
    
    // Return Fragments ID
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Scatterv(&id_frag_list[0], &n_frag_list[0], &frag_recv[0], PS::GetDataType(id_frag_list[0]),
                 &id_frag_loc[0],   n_frag,                        PS::GetDataType(id_frag_loc[0]),  0, MPI_COMM_WORLD);
#else
    for(int i=0; i<n_frag_list[0]; i++) id_frag_loc[i] = id_frag_list[i];
#endif
    //PS::Comm::scatterV(id_frag_list, n_frag_list, frag_recv, id_frag_loc, n_frag);
    assert( (PS::S32)(frag_list.size()) == n_frag );
    for ( PS::S32 i=0; i<n_frag; i++ ){
        std::pair<PS::S32, PS::S32> adr = frag_list.at(i);
        
        if ( ptcl_multi[adr.first][adr.second].isMerged ) {
            PS::S32 sizef = list_multi.at(adr.first).size();
            for ( PS::S32 j=0; j<sizef; j++ ) {
                if ( ptcl_multi[adr.first][j].id == ptcl_multi[adr.first][adr.second].id &&
                     j != adr.second ) {
                    std::pair<bool,PS::S32> adr_j = list_multi.at(adr.first).at(j);
                    if ( adr_j.first ) {
                        pp[adr_j.second].id = id_frag_loc[i];
                        assert ( pp[adr_j.second].isDead );
                    } else {
                        ex_pp[adr_j.second].id = id_frag_loc[i];
                        assert ( ex_pp[adr_j.second].isDead );
                    }
                    ptcl_multi[adr.first][j].id = id_frag_loc[i];
                    assert ( ptcl_multi[adr.first][j].isDead );
                }
            }
        }
        
        ptcl_multi[adr.first][adr.second].id = id_frag_loc[i];
        pp.addOneParticle(ptcl_multi[adr.first][adr.second]);
        assert( ptcl_multi[adr.first][adr.second].time_c == 0. );
    }
    PS::Comm::broadcast(&id_next, 1);

    if ( PS::Comm::getRank() == 0 ){
        // Output Collision Information
        for ( PS::S32 i=0; i<n_col_tot; i++ ) col_list_tot[i].write2File(fp);
            
        //delete [] n_col_list;
        //delete [] n_frag_list;
        //delete [] col_list_tot;
        //delete [] id_frag_list;
        //delete [] col_recv;
        //delete [] frag_recv;
    }
    //delete [] id_frag_loc;

    return n_frag_tot;
}
