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
    PS::F64vec acc_d  = pp[i].mass * pp[i].acc_d;
    PS::F64vec jerk_d = pp[i].mass * pp[i].jerk_d;
    PS::F64    mass   = pp[i].mass;
    for (iterator it = i_range.first; it != i_range.second; ++it){
        PS::S32 j_id = it->second;
        acc_d  += pp[j_id].mass * pp[j_id].acc_d;
        jerk_d += pp[j_id].mass * pp[j_id].jerk_d;
        mass   += pp[j_id].mass;
        assert ( pp[i].id == pp[j_id].id );
        assert ( pp[i].acc_s  == pp[j_id].acc_s );
        assert ( pp[i].jerk_s == pp[j_id].jerk_s );
        assert ( pp[j_id].isDead );
    }
    PS::F64 mass_inv = 1./mass;
    acc_d  *= mass_inv;
    jerk_d *= mass_inv;
    pp[i].acc_d  = acc_d;
    pp[i].jerk_d = jerk_d;
#ifdef INTEGRATE_6TH_SUN
    pp[i].setAcc_();
#endif
    for (iterator it = i_range.first; it != i_range.second; ++it){
        PS::S32 j_id = it->second;
        pp[j_id].acc_d  = acc_d;
        pp[j_id].jerk_d = jerk_d;
#ifdef INTEGRATE_6TH_SUN
        pp[j_id].setAcc_();
#endif
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
PS::S32 collisionDetermination(Tpsys & pp,
                               std::pair<PS::S32,PS::S32> & col_pair)
{
    PS::F64 R_min = 1.;
#ifdef MERGE_BINARY
    PS::F64 E_min = 0.;
#endif
    PS::S32 flag_col = 0;
    
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
            PS::F64 r2 = pp[i].f*pp[i].r_planet + pp[pj_id].f*pp[pj_id].r_planet;
            PS::F64 R = r1 / r2;
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
                flag_col = 1;
            }

#ifdef MERGE_BINARY
            //Determination of Binary Partcles
            if ( !flag_col ){
                PS::F64    mi  = pp[i].mass;
                PS::F64    mj  = pp[pj_id].mass;
                PS::F64vec ri  = pp[i].pos;
                PS::F64vec rj  = pp[pj_id].pos;
                PS::F64vec r_c = (mi*ri + mj*rj) / (mi + mj);
                PS::F64    ax  = sqrt(r_c * r_c);
                PS::F64    r_H = pow((mi+mj)/(3.*FP_t::m_sun), 1./3.) * ax;
 
                PS::F64 R_merge = std::min(pp[i].R_merge, pp[pj_id].R_merge);  
                R = r1 / ( R_merge * r_H );
                if ( R < 1. ){
                    PS::F64vec vi  = pp[i].vel;
                    PS::F64vec vj  = pp[pj_id].vel;
                    PS::F64vec v_c = (mi*vi + mj*vj) / (mi + mj);
                    
                    PS::F64vec ex  = r_c / ax;
                    PS::F64vec ez  = 0.;
                    ez.x = v_c.y*ex.z - v_c.z*ex.y;
                    ez.y = v_c.z*ex.x - v_c.x*ex.z;
                    ez.z = v_c.x*ex.y - v_c.y*ex.x;
                    ez = ez / sqrt(ez*ez) ;
                    
                    PS::F64vec dr = ri - rj;
                    PS::F64vec dv = vi - vj;
                    PS::F64 dr_x = dr * ex;
                    PS::F64 dr_z = dr * ez;
                    
                    PS::F64 Omega = sqrt(FP_t::m_sun / (ax*ax*ax));
                    
                    PS::F64 E_J = 0.5*dv*dv
                        + Omega*Omega*(-1.5*dr_x*dr_x + 0.5*dr_z*dr_z +4.5*r_H*r_H)
                        - (mi + mj) / sqrt(dr*dr);
                    
                    if ( E_J < E_min ){
                        E_min = E_J;
                        if ( pp[i].mass < pp[pj_id].mass ){
                            col_pair.first = i;
                            col_pair.second = pj_id;
                        } else {
                            col_pair.first = pj_id;
                            col_pair.second = i;
                        }
                        flag_col = 2;
                    }
                }
            }
#endif
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
    PS::S64 id_next = 0;
    PS::S32 flag_col = 0;
    n_col = n_frag = 0;
    edisp = edisp_d = 0.;
    collision_list.clear();
    col_pair = std::make_pair(-1,-1);

    PS::S32 a_id = 0;
    PS::S32 j_id = 0;
    PS::F64 t_p = 0.;
    PS::F64 t_c = 0.;

    assert( pp.size() > 0 );
    PS::S32 psize = pp.size();
    PS::S32 asize = 0;
    
    for(PS::S32 i=0; i<psize; i++) calcGravity(pp[i], pp);
#ifdef FORDEBUG
    PS::F64 e0 = calcEnergyCluster(pp);
#endif

    for(PS::S32 i=0; i<psize; i++){
        //calcJerk(pp[i], pp);
        //pp[i].calcDeltatInitial();
        pp[i].isDead = pp[i].isMerged = false;
        assert ( pp[i].time == time );
        assert ( pp[i].dt != 0 );
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
        flag_col = collisionDetermination(pp, col_pair);
        
        if ( flag_col ){
            active_list.clear();
            for(PS::S32 i=0; i<psize; i++){
                pp[i].dt = time_s - pp[i].time;
                if ( pp[i].isDead ) continue;
                active_list.push_back(i);
            }
        }
#endif //COLLISION
            
        // Correct
        asize = active_list.size();
        for(PS::S32 i=0; i<asize; i++){
            a_id = active_list.at(i);
            
            std::pair<iterator, iterator> range;
            calcGravity_p(pp[a_id], pp);
#ifdef INTEGRATE_6TH_SUN
            pp[a_id].setAcc_();
#endif
            if ( pp[a_id].isMerged ){
                range = merge_list.equal_range(a_id);
                for (iterator it = range.first; it != range.second; ++it){
                    j_id = it->second;
                    assert( pp[j_id].isDead );
                    assert( pp[j_id].time == pp[a_id].time );
                    calcGravity_p(pp[j_id], pp);
#ifdef INTEGRATE_6TH_SUN
                    pp[j_id].setAcc_();
#endif
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
#ifdef INTEGRATE_6TH_SUN
                pp[a_id].setAcc_();
#endif
                if ( pp[a_id].isMerged ){
                    for (iterator it = range.first; it != range.second; ++it){
                        j_id = it->second;
                        calcGravity_c(pp[j_id], pp);
#ifdef INTEGRATE_6TH_SUN
                        pp[j_id].setAcc_();
#endif
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
            for(PS::S32 i=0; i<psize; i++) {
                calcGravity(pp[i], pp);
#ifdef INTEGRATE_6TH_SUN
                pp[i].setAcc_();
#endif
            }
            mergeAccJerk(pp, merge_list);
            PS::F64 ekin0, ephi_d0, ephi_s0;
            PS::F64 e2 = calcEnergyCluster(pp, ekin0, ephi_d0, ephi_s0);
            std::cerr << pp.size() << " " << (e2-e0-edisp_d)/e0 << " " << n_col << std::endl;
#endif //FORDEBUG
            ///////////////////
            //   Collision   //
            ///////////////////
            Collision col;
            std::vector<FPHard> pfrag;

            col.inputPair(pp, merge_list, col_pair);
#ifdef MERGE_BINARY
            if ( flag_col == 2 ) {
                n_frag += col.mergeOutcome(pfrag);
            } else {
                n_frag += col.collisionOutcome(pfrag);
            }
#else
            n_frag += col.collisionOutcome(pfrag);
#endif
            col.setParticle(pp, pfrag, merge_list, id_next);
            edisp += col.calcEnergyDissipation(pp, merge_list);
            edisp_d += col.getHardEnergyDissipation();
            col.setNeighbors(pp);

            if ( col.flag_merge ){
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
#ifndef INTEGRATE_6TH_SUN
            for(PS::S32 i=0; i<psize; i++) calcGravity(pp[i], pp);
#else
            for(PS::S32 i=0; i<psize; i++) calcAccJerk(pp[i], pp);
            for(PS::S32 i=0; i<psize; i++) {
                pp[i].setAcc_();
                calcStarSnap(pp[i]);
            }
#endif
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
#endif //COLLISION
        
        time = getSystemTime(pp);
        loop ++;
    }

    psize = pp.size();
    for(PS::S32 i=0; i<psize; i++){
        calcGravity(pp[i], pp);
#ifdef INTEGRATE_6TH_SUN
        pp[i].setAcc_();
#endif
        assert ( pp[i].time == time_end );
    }
    if ( collision_list.size() > 0 ) mergeAccJerk(pp, merge_list);
#ifdef FORDEBUG
    //PS::F64 e1 = calcEnergyCluster(pp);
    //std::cerr << loop << " " << (e1-e0-edisp_d)/e0 << " " << n_col << std::endl;
#endif
}

template <class Tpsys>
void timeIntegrate_multi_omp(Tpsys & pp,
                             PS::F64 time_start,
                             PS::F64 time_end,
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
    PS::S64 id_next = 0;
    PS::S32 flag_col = 0;
    n_col = n_frag = 0;
    edisp = edisp_d = 0.;
    collision_list.clear();
    col_pair = std::make_pair(-1,-1);

    PS::S32 a_id = 0;
    PS::S32 j_id = 0;
    PS::F64 t_p = 0.;
    PS::F64 t_c = 0.;

    assert( pp.size() > 0 );
    PS::S32 psize = pp.size();
    PS::S32 asize = 0;
    
    for(PS::S32 i=0; i<psize; i++) calcGravity(pp[i], pp);
#ifdef FORDEBUG
    PS::F64 e0 = calcEnergyCluster(pp);
#endif

    for(PS::S32 i=0; i<psize; i++){
        //calcJerk(pp[i], pp);
        //pp[i].calcDeltatInitial();
        pp[i].isDead = pp[i].isMerged = false;
        assert ( pp[i].time == time );
        assert ( pp[i].dt != 0. );
    }
    
    while ( time < time_end ) {
        
        time_s = makeActiveList(pp, active_list);

        // Predict
        psize = pp.size();
#pragma omp parallel for private (t_p)
        for ( PS::S32 i=0; i<psize; i++ ) {
            pp[i].xp = 0.0;
            pp[i].vp = 0.0;
            t_p = time_s - pp[i].time;
            pp[i].predict(t_p);
        }

#ifdef COLLISION
        // Collison?
        flag_col = collisionDetermination(pp, col_pair);
        
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
#pragma omp parallel for private (a_id, j_id, t_c)
        for(PS::S32 i=0; i<asize; i++){
            a_id = active_list.at(i);
            
            std::pair<iterator, iterator> range;
            calcGravity_p(pp[a_id], pp);
#ifdef INTEGRATE_6TH_SUN
                pp[a_id].setAcc_();
#endif
            if ( pp[a_id].isMerged ){
                range = merge_list.equal_range(a_id);
                for (iterator it = range.first; it != range.second; ++it){
                    j_id = it->second;
                    assert( pp[j_id].isDead );
                    assert( pp[j_id].time == pp[a_id].time );
                    calcGravity_p(pp[j_id], pp);
#ifdef INTEGRATE_6TH_SUN
                    pp[j_id].setAcc_();
#endif
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
#ifdef INTEGRATE_6TH_SUN
                pp[a_id].setAcc_();
#endif
                if ( pp[a_id].isMerged ){
                    for (iterator it = range.first; it != range.second; ++it){
                        j_id = it->second;
                        calcGravity_c(pp[j_id], pp);
#ifdef INTEGRATE_6TH_SUN
                        pp[j_id].setAcc_();
#endif
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
#ifdef MERGE_BINARY
            if ( flag_col == 2 ) {
                n_frag += col.mergeOutcome(pfrag);
            } else {
                n_frag += col.collisionOutcome(pfrag);
            }
#else
            n_frag += col.collisionOutcome(pfrag);
#endif
            col.setParticle(pp, pfrag, merge_list, id_next);
            edisp += col.calcEnergyDissipation(pp, merge_list);
            edisp_d += col.getHardEnergyDissipation();
            col.setNeighbors(pp);

            if ( col.flag_merge ){
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
#ifndef INTEGRATE_6TH_SUN
            for(PS::S32 i=0; i<psize; i++) calcGravity(pp[i], pp);
#else
            for(PS::S32 i=0; i<psize; i++) calcAccJerk(pp[i], pp);
            for(PS::S32 i=0; i<psize; i++) {
                pp[i].setAcc_();
                calcStarSnap(pp[i]);
            }
#endif
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
#ifdef INTEGRATE_6TH_SUN
        pp[i].setAcc_();
#endif
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
    //pi.jerk_d = 0.;
    //calcStarJerk(pi);
    //pi.calcDeltatInitial();
    assert ( pi.time == time_start );
    assert ( pi.dt != 0. );

    calcStarGravity(pi);
    
    while( pi.time < time_end ){
        pi.predict(pi.dt);

        //pi.acc_s  = 0.;
        pi.acc_d  = 0.;
        //pi.jerk_s = 0.;
        pi.jerk_d = 0.;
        
        calcStarGravity_p(pi);
#ifdef INTEGRATE_6TH_SUN
        pi.setAcc_();
#endif

        ////////////////
        //  iteration
        for ( PS::S32 ite=0; ite<2; ite++ ){
            pi.correct(pi.dt);
            calcStarGravity_c(pi);
#ifdef INTEGRATE_6TH_SUN
            pi.setAcc_();
#endif
        }
        //  iteration
        ////////////////
        
        pi.time += pi.dt;
        pi.calcDeltat();
    }
    pi.phi_s = 0.;
    pi.phi_d = 0.;
    //#ifndef INTEGRATE_6TH_SUN
    calcStarGravity(pi);
    //#else
    //calcStarAccJerk(pi);
    //pi.setAcc_();
    //calcStarSnap(pi);
    //#endif
    assert ( pi.time ==  time_end );
}

template <class Tp>
void timeIntegrateKepler_isolated(Tp & pi,
                                  PS::F64 time_start,
                                  PS::F64 time_end)
{
    const PS::F64 m_sun = FP_t::m_sun;
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
    
    pi.phi_d  = 0.;
    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
#ifndef INTEGRATE_6TH_SUN
    calcStarGravity(pi);
#else
    calcStarAccJerk(pi);
    pi.setAcc_();
    calcStarSnap(pi);
#endif
    pi.calcDeltatInitial();
    assert ( pi.time == time_end );
}

template <class Tp>
void freeMotion(Tp & pi,
                PS::F64 time_start,
                PS::F64 time_end)
{
    pi.pos  += pi.vel * (time_end - time_start);
    pi.time += (time_end - time_start);
    
    pi.phi_d  = 0.;
    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
#ifndef INTEGRATE_6TH_SUN
    calcStarGravity(pi);
#else
    calcStarAccJerk(pi);
    pi.setAcc_();
    calcStarSnap(pi);
#endif
    pi.calcDeltatInitial();
    assert ( pi.time == time_end );
}
