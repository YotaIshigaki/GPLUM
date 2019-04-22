#pragma once


template <class Tpsys>
PS::F64 makeActiveList(Tpsys & pp,
                       std::vector<PS::S32> & active_list)
{
    PS::F64 t_s = 1.e10;
    active_list.clear();
    PS::S32 psize = pp.size();
    for(PS::S32 i=0; i<psize; i++){
        if ( pp[i].isDead ) continue;
        PS::F64 t_i = pp[i].time + pp[i].dt;
        if ( t_s > t_i ) {
            t_s = t_i;
            active_list.clear();
            active_list.push_back(i);
        } else if ( t_s == t_i ) {
            active_list.push_back(i);
        }
    }

    return t_s;
}
template <class Tpsys>
PS::F64 getSystemTime(Tpsys & pp)
{
    PS::F64 t_s = 1.e10;
    PS::S32 psize = pp.size();
    for(PS::S32 i=0; i<psize; i++){
        if ( pp[i].isDead ) continue;
        if ( t_s > pp[i].time ) t_s = pp[i].time;
    }
    return t_s;
}
template <class Tpsys>
void mergeAccJerk(Tpsys & pp,
                  std::multimap<PS::S32,PS::S32> & merge_list,
                  PS::S32 i,
                  std::pair<std::multimap<PS::S32,PS::S32>::iterator,std::multimap<PS::S32,PS::S32>::iterator> & i_range)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    assert ( pp[i].isMerged );
    assert ( i_range.first->first == i );
    PS::F64vec acc = pp[i].mass * pp[i].acc_d;
    PS::F64vec jerk = pp[i].mass * pp[i].jerk;
    PS::F64 mass = pp[i].mass;
    for (iterator it = i_range.first; it != i_range.second; ++it){
        PS::S32 j_id = it->second;
        acc += pp[j_id].mass * pp[j_id].acc_d;
        jerk += pp[j_id].mass * pp[j_id].jerk;
        mass += pp[j_id].mass;
        assert ( pp[i].id == pp[j_id].id );
        assert( pp[j_id].isDead );
    }
    acc /= mass;
    jerk /= mass;
    pp[i].acc_d = acc;
    pp[i].jerk = jerk;
    for (iterator it = i_range.first; it != i_range.second; ++it){
        PS::S32 j_id = it->second;
        pp[j_id].acc_d = acc;
        pp[j_id].jerk = jerk;
    }  
}
template <class Tpsys>
void mergeAccJerk(Tpsys & pp,
                  std::multimap<PS::S32,PS::S32> & merge_list,
                  PS::S32 i)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::pair<iterator,iterator> i_range = merge_list.equal_range(i);
    mergeAccJerk(pp, merge_list, i, i_range);
}
template <class Tpsys>
void mergeAccJerk(Tpsys & pp,
                  std::multimap<PS::S32,PS::S32> & merge_list)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    iterator it = merge_list.begin();
    while ( it != merge_list.end() ){
        PS::S32 i = it->first;
        iterator it2 = merge_list.upper_bound(i);
        std::pair<iterator,iterator> i_range = std::make_pair(it,it2);
        mergeAccJerk(pp, merge_list, i, i_range);
        it = it2;
    }
}
template <class Tpsys>
bool collisionDetermination(Tpsys & pp,
                            std::pair<PS::S32,PS::S32> & col_pair,
                            const PS::F64 f)
{
    PS::F64 R_min = 1.;
    bool flag_col = false;

    //if ( col_pair.first > -1 && col_pair.second > -1 ){
    //    PS::F64vec dr = pp[col_pair.first].xp - pp[col_pair.second].xp;
    //    PS::F64 r1 = sqrt(dr*dr);
    //    PS::F64 r2 = pp[col_pair.first].r_planet + pp[col_pair.second].r_planet;
    //    PS::F64 R = r1 / ( f * r2 );
    //    if ( R > 1. ) col_pair = std::make_pair(-1,-1);
    //}

    PS::S32 psize = pp.size();
    for(PS::S32 i=0; i<psize; i++){
        if ( pp[i].isDead ) continue;
        for(PS::S32 j=0; j<pp[i].neighbor; j++){
            PS::S32 pj_id = pp[i].n_hard_list.at(j);
            if ( pp[pj_id].isDead ) continue;
            PS::F64vec dr = pp[i].xp - pp[pj_id].xp;
            //PS::F64vec dv = pp[i].vp - pp[pj_id].vp;
            PS::F64 r1 = sqrt(dr*dr);
            PS::F64 r2 = pp[i].r_planet + pp[pj_id].r_planet;
            PS::F64 R = r1 / ( f * r2 );
            if ( R < R_min ){
                //if ( ( col_pair == std::make_pair(i, pj_id) || col_pair == std::make_pair(pj_id, i) )
                //     || dr*dv > 0. ) continue;
                R_min = R;
                if ( pp[i].mass < pp[pj_id].mass ){
                    col_pair.first = i;
                    col_pair.second = pj_id;
                } else {
                    col_pair.first = pj_id;
                    col_pair.second = i;
                }
                flag_col = true;
            }
        }
    }
    return flag_col;
}
template <class Tpsys>
PS::F64 calcEnergyCluster(Tpsys & pp,
                          PS::F64 & ekin,
                          PS::F64 & ephi_d,
                          PS::F64 & ephi_s)
{
    ekin = ephi_d = ephi_s = 0.;
    
    for(PS::S32 i=0; i<pp.size(); i++){
        ekin += pp[i].mass * pp[i].vel * pp[i].vel;
        ephi_d += pp[i].mass * pp[i].phi_d;
        ephi_s += pp[i].mass * pp[i].phi_s;
    }
    ekin *= 0.5;
    ephi_d *= 0.5;

    return ekin + ephi_d + ephi_s;
}
template <class Tpsys>
PS::F64 calcEnergyCluster(Tpsys & pp)
{
    PS::F64 ekin = 0.;
    PS::F64 ephi_d = 0.;
    PS::F64 ephi_s = 0.;
    return calcEnergyCluster(pp, ekin, ephi_d, ephi_s);
}
template <class Tpsys>
void timeIntegrate_multi(Tpsys & pp,
                         PS::F64 time_start,
                         PS::F64 time_end,
                         const PS::F64 f,
                         PS::S32 & n_col,
                         PS::S32 & n_frag,
                         PS::F64 & edisp,
                         PS::F64 & edisp_d,
                         std::vector<Collision> & collision_list)
{
    using iterator = std::multimap<PS::S32,PS::S32>::iterator;
    std::vector<PS::S32> active_list;
    std::pair<PS::S32,PS::S32> col_pair;
    std::multimap<PS::S32,PS::S32> merge_list;
    active_list.clear();
    merge_list.clear();
    PS::F64 time = time_start;
    PS::F64 time_s = 0.;
    PS::S32 loop = 0;
    PS::S32 id_next = 0;
    bool flag_col = false;
    n_col = n_frag = 0;
    edisp = edisp_d = 0.;
    collision_list.clear();
    col_pair = std::make_pair(-1,-1);

    PS::S32 a_id = 0;
    PS::S32 j_id = 0;
    PS::F64 t_p = 0.;
    PS::F64 t_c = 0.;

    assert( pp.size() > 0 );
#ifdef FORDEBUG
    for(PS::S32 i=0; i<pp.size(); i++) calcGravity(pp[i], pp);
    PS::F64 e0 = calcEnergyCluster(pp);
#endif

    PS::S32 psize = pp.size();
    PS::S32 asize = 0;
    for(PS::S32 i=0; i<psize; i++){
        calcJerk(pp[i], pp);
        pp[i].calcDeltatInitial();
        pp[i].isDead = pp[i].isMerged = false;
        assert ( pp[i].time == time );
    }
    
    while ( time < time_end ) {
        
        time_s = makeActiveList(pp, active_list);

        // Predict
        psize = pp.size();
        for ( PS::S32 i=0; i<psize; i++ ) {
            pp[i].xp = 0.0;
            pp[i].vp = 0.0;
            t_p = time_s - pp[i].time;
            pp[i].predict(t_p);
        }

#ifdef COLLISION
        // Collison?
        flag_col = collisionDetermination(pp, col_pair, f);
        if ( flag_col ){
            active_list.clear();
            for(PS::S32 i=0; i<psize; i++){
                pp[i].dt = time_s - pp[i].time;
                if ( pp[i].isDead ) continue;
                active_list.push_back(i);
            }
        }
#endif
            
        // Correct
        asize = active_list.size();
        for(PS::S32 i=0; i<asize; i++){
            a_id = active_list.at(i);
            
            std::pair<iterator, iterator> range;
            calcGravity_p(pp[a_id], pp);    
            if ( pp[a_id].isMerged ){
                range = merge_list.equal_range(a_id);
                for (iterator it = range.first; it != range.second; ++it){
                    j_id = it->second;
                    assert( pp[j_id].isDead );
                    assert( pp[j_id].time == pp[a_id].time );
                    calcGravity_p(pp[j_id], pp);
                }
                mergeAccJerk(pp, merge_list, a_id);
            }

            t_c = time_s - pp[a_id].time;
            assert( t_c == pp[a_id].dt );
            ////////////////
            //  iteration
            for (PS::S32 ite=0; ite<2; ite++){
                pp[a_id].correct(t_c);
                if ( pp[a_id].isMerged ){
                    for (iterator it = range.first; it != range.second; ++it){
                        j_id = it->second;
                        pp[j_id].correct(t_c);
                    }
                }
                
                calcGravity_c(pp[a_id], pp);
                if ( pp[a_id].isMerged ){
                    for (iterator it = range.first; it != range.second; ++it){
                        j_id = it->second;
                        calcGravity_c(pp[j_id], pp);
                    }
                    mergeAccJerk(pp, merge_list, a_id);
                }
            }
            //  iteration
            ////////////////
            
            pp[a_id].time += t_c;
            assert( pp[a_id].time == time_s );
            if ( pp[a_id].isMerged ){
                for (iterator it = range.first; it != range.second; ++it){
                    j_id = it->second;
                    pp[j_id].time += t_c;
                    assert( pp[j_id].time == time_s );
                    assert ( pp[a_id].time == pp[j_id].time );
                }
            }
            
            if ( !flag_col ){
                //if ( pp[a_id].time < time_end )
                pp[a_id].calcDeltat();
                if ( pp[a_id].isMerged ){
                    for (iterator it = range.first; it != range.second; ++it){
                        j_id = it->second;
                        pp[j_id].dt = pp[a_id].dt;
                    }
                }
            }
        }
        
#ifdef COLLISION
        if ( flag_col ){
#ifdef FORDEBUG
            psize = pp.size();
            for(PS::S32 i=0; i<psize; i++) calcGravity(pp[i], pp);
            mergeAccJerk(pp, merge_list);
            PS::F64 ekin0, ephi_d0, ephi_s0;
            PS::F64 e2 = calcEnergyCluster(pp, ekin0, ephi_d0, ephi_s0);
            std::cerr << pp.size() << " " << (e2-e0-edisp_d)/e0 << " " << n_col << std::endl;
#endif
            ///////////////////
            //   Collision   //
            ///////////////////
            Collision col;
            std::vector<FPHard> pfrag;

            col.inputPair(pp, merge_list, col_pair);
            n_frag += col.collisionOutcome(pfrag);
            col.setParticle(pp, pfrag, merge_list, id_next);
            edisp += col.calcEnergyDissipation(pp, merge_list);
            edisp_d += col.getHardEnergyDissipation();
            col.setNeighbors(pp);

            if ( !col.HitAndRun ){
                merge_list.insert(std::make_pair(col_pair.second,col_pair.first));
                if ( pp[col_pair.first].isMerged ) {
                    std::pair<iterator, iterator> range = merge_list.equal_range(col_pair.first);
                    for (iterator it = range.first; it != range.second; ++it){
                        merge_list.insert(std::make_pair(col_pair.second,it->second));
                    }
                    merge_list.erase(col_pair.first);
                    pp[col_pair.first].isMerged = false;
                }
            }
            //for(iterator it = merge_list.begin(); it != merge_list.end() ; ++it) {
            //    PRC(it->first);PRC(it->second);PRC(pp[it->first].id);PRL(pp[it->second].id);
            //}
            
            collision_list.push_back(col);
            n_col ++;

            psize = pp.size();
            for(PS::S32 i=0; i<psize; i++) calcGravity(pp[i], pp);
            mergeAccJerk(pp, merge_list);
#ifdef FORDEBUG
            PS::F64 ekin1, ephi_d1, ephi_s1;
            PS::F64 e3 = calcEnergyCluster(pp, ekin1, ephi_d1, ephi_s1);
            std::cerr << pp.size() << " " << (e3-e0-edisp_d)/e0 << " " << n_col << std::endl;
            PRC(ekin1-ekin0);PRC(ephi_d1-ephi_d0);PRL(ephi_s1-ephi_s0);
#endif

            for(PS::S32 i=0; i<psize; i++) {
                if ( !pp[i].isDead ) pp[i].calcDeltatInitial();
                if ( pp[i].isMerged ){
                    std::pair<iterator, iterator> range = merge_list.equal_range(i);
                    for (iterator it = range.first; it != range.second; ++it){
                        j_id = it->second;
                        pp[j_id].dt = pp[i].dt;
                    }
                }
            }
        }
#endif
        time = getSystemTime(pp);
        loop ++;
    }

    psize = pp.size();
    for(PS::S32 i=0; i<psize; i++){
        calcGravity(pp[i], pp);
        assert ( pp[i].time == time_end );
    }
    if ( collision_list.size() > 0 ) mergeAccJerk(pp, merge_list);
#ifdef FORDEBUG
    //PS::F64 e1 = calcEnergyCluster(pp);
    //std::cerr << loop << " " << (e1-e0-edisp_d)/e0 << " " << n_col << std::endl;
#endif
}

template <class Tp>
void timeIntegrate_isolated(Tp & pi,
                         PS::F64 time_start,
                         PS::F64 time_end)
{
    pi.jerk = 0.;
    calcStarJerk(pi);
    pi.calcDeltatInitial();
    assert ( pi.time == time_start );
    
    while( pi.time < time_end ){
        pi.predict(pi.dt);
                    
        pi.acc_d = 0.0;
        pi.jerk = 0.0;
        calcStarGravity_p(pi);

        ////////////////
        //  iteration
        for ( PS::S32 ite=0; ite<2; ite++ ){
            pi.correct(pi.dt);
            pi.acc_d = 0.0;
            pi.jerk = 0.0; 
            calcStarGravity_c(pi);
        }
        //  iteration
        ////////////////
        
        pi.time += pi.dt;
        pi.calcDeltat();
    }
    pi.phi_d = 0.;
    pi.phi_s = 0.;
    pi.acc_d = 0.0;
    pi.jerk = 0.0;
    calcStarGravity(pi);
    assert ( pi.time ==  time_end );
}

template <class Tp>
void timeIntegrateKepler_isolated(Tp & pi,
                                  PS::F64 time_start,
                                  PS::F64 time_end)
{
    const PS::F64 m_sun = FPGrav::m_sun;
    PS::F64 ax, ecc, n;
    PS::F64 u, l;
    PS::F64vec P, Q;
    assert ( pi.time == time_start );
    posVel2OrbitalElement(pi.pos, pi.vel, m_sun, ax, ecc, n, u, P, Q);
    l = KeplerEq(u, ecc);
    l += n * (time_end - time_start);
    u = solveKeplerEq(l, ecc);
    orbitalElement2PosVel(pi.pos, pi.vel, m_sun, ax, ecc, n, u, P, Q);
    pi.time += (time_end - time_start);
    
    pi.phi_d = 0.;
    pi.phi_s = 0.;
    pi.acc_d = 0.;
    pi.jerk = 0.;
    calcStarGravity(pi);
    pi.calcDeltatInitial();
    assert ( pi.time == time_end );
}
