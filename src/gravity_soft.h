#pragma once

#include "cutfunc.h"

template <class Tpsys, class Tptree>
void correctForceLong(Tpsys & pp,
                      Tptree & tree_grav,
                      NeighborList & NList,
                      PS::S32 & n_ngb_tot,
                      PS::S32 & n_with_ngb)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::F64 eps2 = FP_t::eps2;
    NList.initializeList(pp);

#ifdef CORRECT_NEIGHBOR
    PS::F64    phi_d[n_loc];
    PS::F64vec acc_d[n_loc];
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        phi_d[i] = pp[i].phi_d;
        acc_d[i] = pp[i].acc_d;
    }
#endif

    n_ngb_tot    = 0;
    n_with_ngb = 0;
#pragma omp parallel for reduction (+:n_ngb_tot, n_with_ngb)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec acci = 0.;
        PS::F64    phii = 0.;
        PS::F64    acc0i = 0.;
        PS::S64    j_id = 0;
        PS::S32    j_id_local = 0;
        PS::S32    j_rank = 0;
        PS::S32    neighbor = pp[i].neighbor - 1;
        PS::F64vec posi = pp[i].pos;
        pp[i].neighbor = 0.;
        pp[i].id_cluster = pp[i].id;
#if defined(CORRECT_NEIGHBOR) || defined(CHECK_NEIGHBOR)
        PS::F64vec acc_di = 0.;
        PS::F64    phi_di = 0.;
#endif

        if ( neighbor == 0 ) {
#ifdef CHECK_NEIGHBOR
            PS::S32 nei_loc = 0;
            EPJ_t* next = NULL;
            nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
            //assert( nei_loc == 1 );
	    for (PS::S32 j=0; j<nei_loc; j++) {
                PS::F64vec posj = (next+j)->pos;
                PS::F64vec dr   = posj - posi;
                PS::F64 dr2 = dr * dr;
                dr2 += eps2;
                PS::F64 rij  = sqrt(dr2);
                PS::F64 r_search = std::max(pp[i].r_search, (next+j)->r_search);
                if ( rij < r_search && pp[i].id != (next+j)->id ) {
                    std::cerr << "Particle " << pp[i].id << " has the neighbor " << (next+j)->id << std::endl;
	            std::cerr << "Search radii: " << pp[i].r_search << " " << (next+j)->r_search << std::endl;
                    std::cerr << "Distance: " << rij << std::endl;
                }
            }
#endif
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
#endif
            phii += pp[i].mass * r_out_inv;
            
        } else {
            assert ( neighbor > 0 );
            j_id = pp[i].id_neighbor;
            auto iter = NList.id_map.find(j_id);
            bool isMerged = (iter == NList.id_map.end()) ? true : pp[iter->second].isMerged;
            
            if ( !( neighbor > 1 || isMerged ) ) {
#ifdef CHECK_NEIGHBOR
                PS::S32 nei_loc = 0;
                EPJ_t* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                //assert( nei_loc < 3 );
                for (PS::S32 j=0; j<nei_loc; j++) {
                    PS::F64vec posj = (next+j)->pos;
                    PS::F64vec dr   = posj - posi;
                    PS::F64 dr2 = dr * dr;
                    dr2 += eps2;
                    PS::F64 rij  = sqrt(dr2);
                    PS::F64 r_search = std::max(pp[i].r_search, (next+j)->r_search);
                    if ( rij < r_search && pp[i].id != (next+j)->id && j_id != (next+j)->id ) {
                        std::cerr << "Particle " << pp[i].id << " has the neighbor " << (next+j)->id << std::endl;
                        std::cerr << "Search radii: " << pp[i].r_search << " " << (next+j)->r_search << std::endl;
                        std::cerr << "Distance: " << rij << std::endl;
                    }
                }
#endif
                j_id_local = iter->second;
                j_rank     = pp[j_id_local].myrank;
                assert( pp[j_id_local].id == j_id );
                assert( j_id != pp[i].id );
                //assert( !pp[j_id_local].isDead );

#ifdef USE_INDIVIDUAL_CUTOFF
                PS::F64 r_out_inv = pp[i].r_out_inv;
                PS::F64 r_out     = pp[i].r_out;
#else
                PS::F64 r_out_inv = FP_t::r_out_inv;
                PS::F64 r_out     = FP_t::r_out;
#endif
                phii += pp[i].mass * r_out_inv;
                
                PS::F64 massj  = pp[j_id_local].mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                r_out_inv = std::min(pp[i].r_out_inv, pp[j_id_local].r_out_inv);
                r_out     = std::max(pp[i].r_out,     pp[j_id_local].r_out);
                PS::F64 r_search = std::max(pp[i].r_search, pp[j_id_local].r_search);
#else
                PS::F64 r_search = FP_t::r_search;
#endif
                
                PS::F64vec posj = pp[j_id_local].pos;
                PS::F64vec dr   = posj - posi;
                PS::F64 dr2 = dr * dr;
                assert( dr2 != 0.0 );
                assert( j_id > -1 );
                dr2 += eps2;
                PS::F64 rij   = sqrt(dr2);
                
#ifdef USE_RE_SEARCH_NEIGHBOR
                PS::F64vec dv      = pp[j_id_local].vel   - pp[i].vel;
                //PS::F64vec da      = (pp[j_id_local].acc_s - pp[i].acc_s) + (pp[j_id_local].acc_d - pp[i].acc_d);
                PS::F64vec da      = pp[j_id_local].acc_d - pp[i].acc_d;
                PS::F64    drdv    = dr * dv;
                PS::F64    dv2     = dv * dv;
                PS::F64    da2     = da * da;
                PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FP_t::dt_tree);
                PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                //assert( dr2 >= dr2_min );
                
                PS::F64    r_crit   = FP_t::R_search2 * r_out;
                PS::F64    v_crit_a = FP_t::R_search3*0.5*FP_t::dt_tree;
                if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 ) {
                    //if ( dr2_min < r_crit * r_crit ) {
#endif
                    if ( rij < r_search ) {
                        NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                        acc0i += r_out * r_out / massj;
                    }
#ifdef USE_RE_SEARCH_NEIGHBOR
                }
#endif
                
                if ( rij < r_out ) {    
                    PS::F64 rinv  = 1. / rij;
                    PS::F64 r2inv = rinv * rinv;
                    PS::F64 r3inv = rinv * r2inv;
                    
                    PS::F64 W  = cutoff_W(rij, r_out_inv);
                    PS::F64 K  = cutoff_K(rij, r_out_inv);
                    PS::F64 r_min = std::min(rinv, r_out_inv);
                    
#if defined(CORRECT_NEIGHBOR) || defined(CHECK_NEIGHBOR)
                    phi_di -= massj * rinv * (1.-W);
                    acc_di += massj * r3inv * (1.-K) * dr;
#endif
                    phii  -= massj * ( rinv * W  - r_min );
                    acci  += massj * ( r3inv * K - r_min * r_min * r_min ) * dr;
                }
                
            } else { // neighbor > 2 || j_id is not in map
                PS::S32 nei_loc = 0;
                EPJ_t* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                
                for(PS::S32 j=0; j<nei_loc; j++){
                    j_id       = (next+j)->id;
                    j_id_local = (next+j)->id_local;
                    j_rank     = (next+j)->myrank;
                    //j_id_local = ( pp[i].myrank == j_rank ) ? NList.id_map.at(j_id) : -1;
                
                    PS::F64 massj  = (next+j)->mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                    //PS::F64 r_out_inv = std::min(pp[i].r_out_inv, (next+j)->r_out_inv);
                    PS::F64 r_out     = std::max(pp[i].r_out,     (next+j)->r_out);
                    PS::F64 r_out_inv = 1. / r_out;
                    PS::F64 r_search  = std::max(pp[i].r_search,  (next+j)->r_search);
#else
                    PS::F64 r_out_inv = FP_t::r_out_inv;
                    PS::F64 r_out     = FP_t::r_out;
                    PS::F64 r_search  = FP_t::r_search;
#endif
                    
                    if ( j_id == pp[i].id ) {
                        phii += massj * r_out_inv;
                        continue;
                    }
            
                    PS::F64vec posj = (next+j)->pos;
                    PS::F64vec dr   = posj - posi;
                    PS::F64 dr2 = dr * dr;              
                    assert( dr2 != 0.0 );
                    assert( j_id > -1 );
                    dr2 += eps2;
                    PS::F64 rij   = sqrt(dr2);

#ifdef USE_RE_SEARCH_NEIGHBOR
                    PS::F64vec dv      = (next+j)->vel   - pp[i].vel;
                    //PS::F64vec da      = ((next+j)->acc_s - pp[i].acc_s) + ((next+j)->acc_d - pp[i].acc_d);
                    PS::F64vec da      = (next+j)->acc_d - pp[i].acc_d;
                    PS::F64    drdv    = dr * dv;
                    PS::F64    dv2     = dv * dv;
                    PS::F64    da2     = da * da;
                    PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FP_t::dt_tree);
                    PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                    //assert( dr2 >= dr2_min );

                    PS::F64    r_crit   = FP_t::R_search2 * r_out;
                    PS::F64    v_crit_a = FP_t::R_search3*0.5*FP_t::dt_tree;
                    if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 ){
                        //if ( dr2_min < r_crit * r_crit ){
#endif
                        if ( rij < r_search ) {
                            NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                            acc0i += r_out * r_out / massj;
                        }
#ifdef USE_RE_SEARCH_NEIGHBOR
                    }
#endif
                    
                    if ( rij < r_out ) {
                        PS::F64 rinv  = 1. / rij;
                        PS::F64 r2inv = rinv * rinv;
                        PS::F64 r3inv = rinv * r2inv;
                        
                        PS::F64 W  = cutoff_W(rij, r_out_inv);
                        PS::F64 K  = cutoff_K(rij, r_out_inv);
                        PS::F64 r_min = std::min(rinv, r_out_inv);
                        
#if defined(CORRECT_NEIGHBOR) || defined(CHECK_NEIGHBOR)
                        phi_di -= massj * rinv * (1.-W);
                        acc_di += massj * r3inv * (1.-K) * dr;
#endif
                        phii  -= massj * ( rinv * W  - r_min );
                        acci  += massj * ( r3inv * K - r_min * r_min * r_min ) * dr;
                    }
                }
            }
        }

        if ( pp[i].neighbor ){
#pragma omp critical
            {
                NList.with_neighbor_list.push_back(i);
            }
        }

        n_ngb_tot += pp[i].neighbor;
        n_with_ngb ++;

        pp[i].acc  += acci;
        pp[i].phi  += phii;
        pp[i].acc0  = ( acc0i > 0. ) ? pp[i].neighbor/acc0i : 0.;
#ifdef CORRECT_NEIGHBOR
        pp[i].phi_correct = phi_di - phi_d[i];
        pp[i].acc_correct = acc_di - acc_d[i];
        //pp[i].acc_d  += acc_di;
#endif
        
#ifdef CHECK_NEIGHBOR
        PS::F64  dphi_d = pp[i].phi_d - phi_di;
        if ( abs((pp[i].phi_d - phi_di)/std::max(pp[i].phi_d, phi_di)) > 1.e-14 ) {
            std::cout << pp[i].id << "\t" << pp[i].phi_d << "\t"
                      << phi_di << "\t" << dphi_d/pp[i].phi_s << std::endl;

        }
#endif
    }
    //NList.checkNeighbor(pp);
}

template <class Tpsys, class Tptree>
void correctForceLongInitial(Tpsys & pp,
                             Tptree & tree_grav,
                             NeighborList & NList,
                             PS::S32 & n_ngb_tot,
                             PS::S32 & n_with_ngb)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::F64 eps2 = FP_t::eps2;
    NList.initializeList(pp);

    PS::F64vec acc_s[n_loc];
    PS::F64vec acc_d[n_loc];
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        acc_s[i] = pp[i].acc_s;
        acc_d[i] = pp[i].acc_d;
    }

    n_ngb_tot    = 0;
    n_with_ngb = 0;
#pragma omp parallel for reduction (+:n_ngb_tot, n_with_ngb)
    for(PS::S32 i=0; i<n_loc; i++){
        //pp[i].acc_s  = 0.;
        //pp[i].acc_d  = 0.;
        //pp[i].jerk_s = 0.;
        //pp[i].jerk_d = 0.;
        //pp[i].phi_s  = 0.;
        //pp[i].phi_d  = 0.; 

#ifndef INTEGRATE_6TH_SUN
        calcStarGravity(pp[i]);
#else
        calcStarAccJerk(pp[i]);
#endif

        PS::F64vec acci   = 0.;
        PS::F64vec acc_di = 0.;
        PS::F64vec jerki  = 0.;
        PS::F64    acc0i  = 0.;
        PS::F64    phii   = 0.;
        PS::F64    phi_di = 0.;
        PS::S64    j_id = 0;
        PS::S32    j_id_local = 0;
        PS::S32    j_rank = 0;
        PS::S32    neighbor = pp[i].neighbor - 1;
        PS::F64vec posi = pp[i].pos;
        pp[i].neighbor = 0.;
        pp[i].id_cluster = pp[i].id;
        

        if ( neighbor == 0 ) {      
#ifdef CHECK_NEIGHBOR
            PS::S32 nei_loc = 0;
            EPJ_t* next = NULL;
            nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
            //assert( nei_loc == 1 );
            for (PS::S32 j=0; j<nei_loc; j++) {
                PS::F64vec posj = (next+j)->pos;
                PS::F64vec dr   = posj - posi;
                PS::F64 dr2 = dr * dr;
                dr2 += eps2;
                PS::F64 rij  = sqrt(dr2);
                PS::F64 r_search = std::max(pp[i].r_search, (next+j)->r_search);
                if ( rij < r_search && pp[i].id != (next+j)->id ) {
                    std::cerr << "Particle " << pp[i].id << " has the neighbor " << (next+j)->id << std::endl;
                    std::cerr << "Search radii: " << pp[i].r_search << " " << (next+j)->r_search << std::endl;
                    std::cerr << "Distance: " << rij << std::endl;
                }
            }
#endif
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
#endif
            phii += pp[i].mass * r_out_inv;
            
        } else {
            assert ( neighbor > 0 );
            j_id = pp[i].id_neighbor;
            auto iter = NList.id_map.find(j_id);
            bool isMerged = (iter == NList.id_map.end()) ? true : pp[iter->second].isMerged;
            
            if ( !( neighbor > 1 || isMerged ) ) {          
#ifdef CHECK_NEIGHBOR
                PS::S32 nei_loc = 0;
                EPJ_t* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                //assert( nei_loc < 3 );
                for (PS::S32 j=0; j<nei_loc; j++) {
                    PS::F64vec posj = (next+j)->pos;
                    PS::F64vec dr   = posj - posi;
                    PS::F64 dr2 = dr * dr;
                    dr2 += eps2;
                    PS::F64 rij  = sqrt(dr2);
                    PS::F64 r_search = std::max(pp[i].r_search, (next+j)->r_search);
                    if ( rij < r_search && pp[i].id != (next+j)->id && j_id != (next+j)->id ) {
                        std::cerr << "Particle " << pp[i].id << " has the neighbor " << (next+j)->id << std::endl;
                        std::cerr << "Search radii: " << pp[i].r_search << " " << (next+j)->r_search << std::endl;
                        std::cerr << "Distance: " << rij << std::endl;
                    }
                }
#endif
                j_id_local = iter->second;
                j_rank     = pp[j_id_local].myrank;
                assert( pp[j_id_local].id == j_id );
                assert( j_id != pp[i].id );
                //assert( !pp[j_id_local].isDead );

#ifdef USE_INDIVIDUAL_CUTOFF
                PS::F64 r_out_inv = pp[i].r_out_inv;
                PS::F64 r_out     = pp[i].r_out;
#else
                PS::F64 r_out_inv = FP_t::r_out_inv;
                PS::F64 r_out     = FP_t::r_out;
#endif
                phii += pp[i].mass * r_out_inv;
                
                PS::F64 massj  = pp[j_id_local].mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                r_out_inv = std::min(pp[i].r_out_inv, pp[j_id_local].r_out_inv);
                r_out     = std::max(pp[i].r_out,     pp[j_id_local].r_out);
                PS::F64 r_search = std::max(pp[i].r_search, pp[j_id_local].r_search);
#else
                PS::F64 r_search = FP_t::r_search;
#endif
                
                PS::F64vec posj = pp[j_id_local].pos;
                PS::F64vec dr   = posj - posi;
                PS::F64 dr2 = dr * dr;
                assert( dr2 != 0.0 );
                assert( j_id > -1 );
                dr2 += eps2;
                PS::F64 rij   = sqrt(dr2);
                
                PS::F64vec dv      = pp[j_id_local].vel   - pp[i].vel;
                PS::F64    drdv    = dr * dv;
#ifdef USE_RE_SEARCH_NEIGHBOR
                //PS::F64vec da      = (acc_s[j_id_local] - acc_s[i]) + (acc_d[j_id_local] - acc_d[i]);
                PS::F64vec da      = acc_d[j_id_local] - acc_d[i];
                PS::F64    dv2     = dv * dv;
                PS::F64    da2     = da * da;
                PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FP_t::dt_tree);
                PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                //assert( dr2 >= dr2_min );
                
                PS::F64    r_crit   = FP_t::R_search2 * r_out;
                PS::F64    v_crit_a = FP_t::R_search3*0.5*FP_t::dt_tree;
                if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 || da2 == 0. ){
                //if ( dr2_min < r_crit * r_crit ){
#endif
                    if ( rij < r_search ) {
                        NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                        acc0i += r_out * r_out / massj;
                    }
#ifdef USE_RE_SEARCH_NEIGHBOR
                }
#endif

                if ( rij < r_out ) {
                    PS::F64 rinv  = 1. / rij;
                    PS::F64 r2inv = rinv * rinv;
                    PS::F64 r3inv = rinv * r2inv;
                    //PS::F64 r5inv = r3inv * r2inv;

                    PS::F64 alpha = drdv * r2inv;
                    
                    PS::F64 W    = cutoff_W(rij, r_out_inv);
                    PS::F64 K    = cutoff_K(rij, r_out_inv);
                    PS::F64 dKdt = cutoff_dKdt(rij, r_out_inv, alpha);
                    PS::F64 alpha_c = alpha*(1.-K);
                    PS::F64 r_min = std::min(rinv, r_out_inv);
                    
                    phii   -= massj * ( rinv * W - r_min );
                    phi_di -= massj * rinv * (1.-W);
                    acci   += massj * ( r3inv * K - r_min * r_min * r_min ) * dr;
                    acc_di += massj * r3inv * (1.-K) * dr;
                    jerki  += massj * r3inv * ( (1.-K) * dv - (3.*alpha_c + dKdt) * dr );
                }

            } else { // neighbor > 2 || j_id is not in map
                PS::S32 nei_loc = 0;
                EPJ_t* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                
                for(PS::S32 j=0; j<nei_loc; j++){
                    j_id       = (next+j)->id;
                    j_id_local = (next+j)->id_local;
                    j_rank     = (next+j)->myrank;
                    //j_id_local = ( pp[i].myrank == j_rank ) ? NList.id_map.at(j_id) : -1;
                
                    PS::F64 massj  = (next+j)->mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                    //PS::F64 r_out_inv = std::min(pp[i].r_out_inv, (next+j)->r_out_inv);
                    PS::F64 r_out     = std::max(pp[i].r_out,     (next+j)->r_out);
                    PS::F64 r_out_inv = 1. / r_out;
                    PS::F64 r_search  = std::max(pp[i].r_search,  (next+j)->r_search);
#else
                    PS::F64 r_out_inv = FP_t::r_out_inv;
                    PS::F64 r_out     = FP_t::r_out;
                    PS::F64 r_search  = FP_t::r_search;
#endif

                    if ( j_id == pp[i].id ) {
                        phii += massj * r_out_inv;
                        continue;
                    }
            
                    PS::F64vec posj = (next+j)->pos;
                    PS::F64vec dr   = posj - posi;
                    PS::F64 dr2 = dr * dr;              
                    assert( dr2 != 0.0 );
                    assert( j_id > -1 );
                    dr2 += eps2;
                    PS::F64 rij   = sqrt(dr2);
                    
                    PS::F64vec dv      = (next+j)->vel   - pp[i].vel;
                    PS::F64    drdv    = dr * dv;
#ifdef USE_RE_SEARCH_NEIGHBOR
                    //PS::F64vec da      = ((next+j)->acc_s - acc_s[i]) + ((next+j)->acc_d - acc_d[i]);
                    PS::F64vec da      = (next+j)->acc_d - acc_d[i];
                    PS::F64    dv2     = dv * dv;
                    PS::F64    da2     = da * da;
                    PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FP_t::dt_tree);
                    PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                    //assert( dr2 >= dr2_min );
                    
                    PS::F64    r_crit   = FP_t::R_search2 * r_out;
                    PS::F64    v_crit_a = FP_t::R_search3*0.5*FP_t::dt_tree;
                    if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 || da2 == 0. ){
                    //if ( dr2_min < r_crit * r_crit ){
#endif
                        if ( rij < r_search ) {
                            NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                            acc0i += r_out * r_out / massj;
                        }
#ifdef USE_RE_SEARCH_NEIGHBOR
                    }
#endif
                    
                    if ( rij < r_out ) {
                        PS::F64 rinv  = 1. / rij;
                        PS::F64 r2inv = rinv * rinv;
                        PS::F64 r3inv = rinv * r2inv;
                        //PS::F64 r5inv = r3inv * r2inv;

                        PS::F64 alpha = drdv * r2inv;
                    
                        PS::F64 W    = cutoff_W(rij, r_out_inv);
                        PS::F64 K    = cutoff_K(rij, r_out_inv);
                        PS::F64 dKdt = cutoff_dKdt(rij, r_out_inv, alpha);
                        PS::F64 alpha_c = alpha*(1.-K);
                        PS::F64 r_min = std::min(rinv, r_out_inv);
                    
                        phii   -= massj * ( rinv * W - r_min );
                        phi_di -= massj * rinv * (1.-W);
                        acci   += massj * ( r3inv * K - r_min * r_min * r_min ) * dr;
                        acc_di += massj * r3inv * (1.-K) * dr;
                        jerki  += massj * r3inv * ( (1.-K) * dv - (3.*alpha_c + dKdt) * dr );
                    }
                }
            }
        }

        if ( pp[i].neighbor ){
#pragma omp critical
            {
                NList.with_neighbor_list.push_back(i);
            }
        }
        
        n_ngb_tot += pp[i].neighbor;
        n_with_ngb ++;
        
        pp[i].acc   += acci;
        pp[i].phi   += phii;
        pp[i].acc_d  = acc_di;
        pp[i].phi_d  = phi_di;
        pp[i].jerk_d = jerki;
        pp[i].acc0   = ( acc0i > 0. ) ? pp[i].neighbor/acc0i : 0.;

#ifndef INTEGRATE_6TH_SUN
        pp[i].calcDeltatInitial();
#endif
    }
#ifdef INTEGRATE_6TH_SUN
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].setAcc_();
        calcStarSnap(pp[i]);
        pp[i].calcDeltatInitial();
    }
#endif
    //NList.checkNeighbor(pp);
}


#ifdef INDIRECT_TERM
template <class Tpsys>
void calcIndirectTerm(Tpsys & pp)
{
    PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::F64 eps2 = FP_t::eps2;
    
    PS::F64vec acc_loc = 0;
    PS::F64vec pos_grav_loc = 0;
    PS::F64vec vel_grav_loc = 0;
    PS::F64    mass_tot_loc = 0;
    PS::F64    mass_tot_glb = 0;
#pragma omp parallel for reduction (+:acc_loc, pos_grav_loc, vel_grav_loc, mass_tot_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        
        PS::F64vec posi  = pp[i].pos;
        PS::F64vec veli  = pp[i].vel;
        PS::F64    massi = pp[i].mass;
        
        PS::F64vec dr    = - posi;
        PS::F64    r2inv = 1.0 / (dr * dr + eps2);
        PS::F64    rinv  = sqrt(r2inv);
        PS::F64    r3inv = rinv * r2inv;
        
        acc_loc += massi * r3inv * dr;
        
        pos_grav_loc += massi * posi;
        vel_grav_loc += massi * veli;
        mass_tot_loc += massi;
    }
    
    FP_t::acc_indirect = PS::Comm::getSum(acc_loc);
    
    FP_t::mass_tot = mass_tot_glb = PS::Comm::getSum(mass_tot_loc) + FP_t::m_sun;
    FP_t::pos_g = PS::Comm::getSum(pos_grav_loc) / mass_tot_glb;
    FP_t::vel_g = PS::Comm::getSum(vel_grav_loc) / mass_tot_glb;
}

template <class Tpsys>
PS::F64 calcIndirectEnergy(Tpsys & pp)
{
    PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::F64 eps2 = FP_t::eps2;
    
    PS::F64vec vel_grav_loc = 0;
    PS::F64vec vel_grav_glb = 0;
    PS::F64    mass_tot_loc = 0;
    PS::F64    mass_tot_glb = 0;
#pragma omp parallel for reduction (+:pos_grav_loc, vel_grav_loc, mass_tot_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec posi  = pp[i].pos;
        PS::F64vec veli  = pp[i].vel;
        PS::F64    massi = pp[i].mass;
        pos_grav_loc += massi * posi;
        vel_grav_loc += massi * veli;
        mass_tot_loc += massi;
    }
    mass_tot_glb = PS::Comm::getSum(mass_tot_loc) + FP_t::m_sun;
    vel_grav_glb = PS::Comm::getSum(vel_grav_loc) / mass_tot_glb;
    
    return -0.5*(mass_tot_glb + FP_t::m_sun) * vel_grav_glb*vel_grav_glb;
}

template <class Tpsys>
PS::F64 getIndirectEnergy(Tpsys & pp)
{
    return -0.5 * (FP_t::mass_tot + FP_t::m_sun) * FP_t::vel_g * FP_t::vel_g;
}
#endif


 
