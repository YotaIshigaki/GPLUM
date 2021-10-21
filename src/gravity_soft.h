#pragma once

#include "cutfunc.h"

template <class Tpsys>
void velKick(Tpsys & pp){
    const PS::S32 n = pp.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        pp[i].velKick();
    }
}

#ifdef CORRECT_NEIGHBOR
template <class Tpsys>
void velKick2nd(Tpsys & pp){
    const PS::S32 n = pp.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        pp[i].velKick2nd();
    }
}
#endif


template<class Tp0, class Tp1> inline PS::F64 getROutInvMin(Tp0 & p0, Tp1 & p1) { return std::min(p0.r_out_inv, p1.r_out_inv); }
template<> inline PS::F64 getROutInvMin<FP_t, EPJ_t>(FP_t & p0, EPJ_t & p1) { return std::min(p0.r_out_inv, 1./p1.r_out); }
template<> inline PS::F64 getROutInvMin<EPJ_t, FP_t>(EPJ_t & p0, FP_t & p1) { return std::min(1./p0.r_out, p1.r_out_inv); }
template<> inline PS::F64 getROutInvMin<EPJ_t, EPJ_t>(EPJ_t & p0, EPJ_t & p1) { return std::min(1./p0.r_out, 1./p1.r_out); }
template<class Tp0, class Tp1> inline PS::F64 getROutInverseMin(Tp0 & p0, Tp1 & p1){ return getROutInvMin<Tp0, Tp1>(p0, p1); }

template <class Tpsys, class Tp>
void correctForceBetween2Particles(Tpsys &   pp,
                                   PS::S32 & i,
                                   Tp &      epj,
                                   PS::F64 &    phii,
                                   PS::F64vec & acci,
                                   PS::F64 &    acc0i,
                                   NeighborList & NList)
{
    const PS::F64 eps2 = FP_t::eps2;
    
    PS::S64 j_id       = epj.id;
    PS::S32 j_id_local = epj.id_local;
    PS::S32 j_rank     = epj.myrank;
    //if (j_rank == pp[i].myrank && !pp[i].isMerged) assert( NList.id_map.at(j_id) == j_id_local );

    
    PS::F64 massj  = epj.mass;
#ifdef USE_INDIVIDUAL_CUTOFF
    PS::F64 r_out     = std::max(pp[i].r_out,     epj.r_out);
    //PS::F64 r_out_inv = std::min(pp[i].r_out_inv, epj.r_out_inv);
    PS::F64 r_out_inv = getROutInverseMin(pp[i], epj);
    PS::F64 r_search  = std::max(pp[i].r_search,  epj.r_search);
#else
    PS::F64 r_out     = FP_t::r_out;
    PS::F64 r_out_inv = FP_t::r_out_inv;
    PS::F64 r_search  = FP_t::r_search;
#endif

    if ( j_id == pp[i].id ) {
        phii += massj * r_out_inv;
        return;
    }

    //PS::F64vec posj = epj.pos;
    PS::F64vec dr  = epj.pos - pp[i].pos;
    PS::F64    dr2 = dr * dr;
    assert( dr2 != 0.0 );
    assert( j_id > -1 );
    dr2 += eps2;
    PS::F64 rij = sqrt(dr2);
    assert( rij < r_search * 1.2 );
    
#ifdef USE_RE_SEARCH_NEIGHBOR
    PS::F64vec dv      = epj.vel   - pp[i].vel;
    PS::F64vec da      = epj.acc_d - pp[i].acc_d;
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
        
        phii  -= massj * ( rinv * W  - r_min );
        acci  += massj * ( r3inv * K - r_min * r_min * r_min ) * dr;
    }
}

template <class Tpsys, class Tp>
void correctForceBetween2ParticlesInitial(Tpsys &   pp,
                                          PS::S32 & i,
                                          Tp &      epj,
                                          PS::F64 &    phii,
                                          PS::F64 &    phi_di,
                                          PS::F64vec & acci,
                                          PS::F64vec & acc_di,
                                          PS::F64vec & jerki,
                                          PS::F64 &    acc0i,
                                          NeighborList & NList)
{
    const PS::F64 eps2 = FP_t::eps2;
    
    PS::S64 j_id       = epj.id;
    PS::S32 j_id_local = epj.id_local;
    PS::S32 j_rank     = epj.myrank;
    //if (j_rank == pp[i].myrank && !pp[i].isMerged) assert( NList.id_map.at(j_id) == j_id_local );
                
    PS::F64 massj  = epj.mass;
#ifdef USE_INDIVIDUAL_CUTOFF
    PS::F64 r_out     = std::max(pp[i].r_out,     epj.r_out);
    //PS::F64 r_out_inv = std::min(pp[i].r_out_inv, epj.r_out_inv);
    PS::F64 r_out_inv = getROutInverseMin(pp[i], epj);
    PS::F64 r_search  = std::max(pp[i].r_search,  epj.r_search);  
#else
    PS::F64 r_out     = FP_t::r_out;
    PS::F64 r_out_inv = FP_t::r_out_inv;
    PS::F64 r_search  = FP_t::r_search;
#endif

    if ( j_id == pp[i].id ) {
        phii += massj * r_out_inv;
        return;
    }
            
    PS::F64vec dr  = epj.pos - pp[i].pos;
    PS::F64    dr2 = dr * dr;              
    assert( dr2 != 0.0 );
    assert( j_id > -1 );
    dr2 += eps2;
    PS::F64 rij   = sqrt(dr2);
    assert( rij < r_search * 1.2 );
                    
    PS::F64vec dv      = epj.vel - pp[i].vel;
    PS::F64    drdv    = dr * dv;
#ifdef USE_RE_SEARCH_NEIGHBOR
    PS::F64vec da      = epj.acc_d - pp[i].acc_d;
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


template <class Tpsys, class Tptree>
void correctForceLong(Tpsys & pp,
                      Tptree & tree_grav,
                      NeighborList & NList,
                      PS::S32 & n_ngb_tot,
                      PS::S32 & n_with_ngb)
{
    const PS::S32 n_loc  = pp.getNumberOfParticleLocal();
    const PS::F64 eps2 = FP_t::eps2;
	//auto time0 = PS::GetWtime();    
    NList.initializeList(pp);
    //auto time1 = PS::GetWtime();
    
    //NList.makeTemporaryNeighborList(pp, tree_grav);
    //auto time2 = PS::GetWtime();    
    //NList.exchangeNeighborInfo(pp);
    //auto time3 = PS::GetWtime();

    n_ngb_tot    = 0;
    n_with_ngb = 0;
#pragma omp parallel for reduction (+:n_ngb_tot, n_with_ngb)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec acci = 0.;
        PS::F64    phii = 0.;
        PS::F64    acc0i = 0.;
        PS::S32    neighbor = pp[i].neighbor.number;
        pp[i].neighbor.number = 0.;
        pp[i].id_cluster = pp[i].id;

        if ( neighbor == 0 ) {
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
#endif
            phii += pp[i].mass * r_out_inv;
            
            //} else if ( neighbor == 1 ) {
        } else {
            assert ( neighbor > 0 );
            //EPNgb epj;
            
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
#endif
            phii += pp[i].mass * r_out_inv;
            
            
            if ( neighbor > 2 || ! pp[i].neighbor.inDomain() ) {
                EPJ_t* next = NULL;
                neighbor = tree_grav.getNeighborListOneParticle(pp[i], next);
                
                for ( PS::S32 j=0; j<neighbor; j++ ) {
                    if ( next[j].id ==  pp[i].id ) continue;
                    correctForceBetween2Particles(pp, i, next[j], phii, acci, acc0i, NList);
                }
            } else {
                for ( PS::S32 j=0; j<neighbor; j++ ) {
                    PS::S32 id_loc = pp[i].neighbor.getId(j);
                    assert( id_loc != i );
                    correctForceBetween2Particles(pp, i, pp[id_loc], phii, acci, acc0i, NList);
                }
            }

            /*
            for ( PS::S32 j=0; j<neighbor; j++ ) {
                PS::S32 rank   = NList.n_list_tmp[i][j].rank;
                PS::S32 id_loc = NList.n_list_tmp[i][j].id_local;
                assert( rank > -1 );
                assert( id_loc > -1 );
                
                if ( pp[i].myrank == rank ){
                    epj.copyFromFP(pp[id_loc]);
                } else {
                    epj = NList.getNeighborInfo(rank, id_loc);
                    //epj.dump();
                }
                assert(epj.myrank == rank);
                assert(epj.id_local == id_loc);
                
            
                //PS::S32 j = 0;
                //if ( pp[i].neighbor.getRank(j) == pp[i].myrank ){
                //epj.setNeighbor(pp[i].neighbor.getId(j));
                //epj.copyNeighbor(pp);
                //} else {
                //epj = NList.getNeighborInfo(pp[i].neighbor.getRank(j), pp[i].neighbor.getId(j));
                //}
            
                correctForceBetween2Particles(pp, i, epj, phii, acci, acc0i, NList);
            
            //} else { // neighbor > 2 
            //PS::S32 n_ngb = 0;
            //EPJ_t* next = NULL;
            //n_ngb = tree_grav.getNeighborListOneParticle(pp[i], next);
                
            //for(PS::S32 j=0; j<n_ngb; j++){
            //    correctForceBetween2Particles(pp, i, next[j], phii, acci, acc0i, NList);
            //}
            
            }
            */
        }

        if ( pp[i].neighbor.number ){
#pragma omp critical
            {
                NList.with_neighbor_list.push_back(i);
            }

            n_ngb_tot += pp[i].neighbor.number;
            n_with_ngb ++;
        }

        pp[i].acc  += acci;
        pp[i].phi  += phii;
        pp[i].acc0  = ( acc0i > 0. ) ? pp[i].neighbor.number/acc0i : 0.;
    }
    //NList.checkNeighbor(pp);
    //auto time4 = PS::GetWtime();
    //std::cerr << time1-time0 << " " << time2-time1 << " "
    //          << time3-time2 << " " << time4-time3 << std::endl;
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

    //NList.makeTemporaryNeighborList(pp, tree_grav);
    //NList.exchangeNeighborInfo(pp);

    PS::F64vec acc_d[n_loc];
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        acc_d[i] = pp[i].acc_d;
    }

    n_ngb_tot    = 0;
    n_with_ngb = 0;
#pragma omp parallel for reduction (+:n_ngb_tot, n_with_ngb)
    for(PS::S32 i=0; i<n_loc; i++){

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
        PS::S32    neighbor = pp[i].neighbor.number;
        pp[i].neighbor.number = 0.;
        pp[i].id_cluster = pp[i].id;
        

        if ( neighbor == 0 ) {      
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
#endif
            phii += pp[i].mass * r_out_inv;
            
            //} else if ( neighbor == 1 ) {
        } else {
            assert ( neighbor > 0 );
            //EPNgb epj;
            
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
#endif
            phii += pp[i].mass * r_out_inv;
            
            
            if ( neighbor > 2 || ! pp[i].neighbor.inDomain() ) {
                EPJ_t* next = NULL;
                neighbor = tree_grav.getNeighborListOneParticle(pp[i], next);
                
                for ( PS::S32 j=0; j<neighbor; j++ ) {
                    if ( next[j].id ==  pp[i].id ) continue;
                    correctForceBetween2ParticlesInitial(pp, i, next[j], phii, phi_di, acci, acc_di, jerki, acc0i, NList);
                }
            } else {
                for ( PS::S32 j=0; j<neighbor; j++ ) {
                    PS::S32 id_loc = pp[i].neighbor.getId(j);
                    assert( id_loc != i );
                    correctForceBetween2ParticlesInitial(pp, i, pp[id_loc], phii, phi_di, acci, acc_di, jerki, acc0i, NList);
                }
            }
            /*
            PS::U32 n_size = neighbor;
            for ( PS::S32 j=0; j<n_size; j++ ) {
                PS::S32 rank   = NList.n_list_tmp[i][j].rank;
                PS::S32 id_loc = NList.n_list_tmp[i][j].id_local;
                assert( rank > -1 );
                assert( id_loc > -1 );
                
                if ( pp[i].myrank == rank ){
                    epj.copyFromFP(pp[id_loc]);
                } else {
                    epj = NList.getNeighborInfo(rank, id_loc);
                }
                assert(epj.myrank == rank);
                assert(epj.id_local == id_loc);
                
                //PS::S32 j = 0;
                //if ( pp[i].neighbor.getRank(j) == pp[i].myrank ){
                //epj.setNeighbor(pp[i].neighbor.getId(j));
                //epj.copyNeighbor(pp);
                //} else {
                //epj = NList.getNeighborInfo(pp[i].neighbor.getRank(j), pp[i].neighbor.getId(j));
                //}
            
                correctForceBetween2ParticlesInitial(pp, i, epj, phii, phi_di, acci, acc_di, jerki, acc0i, NList);

            //} else { // neighbor > 2 
            //PS::S32 n_ngb = 0;
            //EPJ_t* next = NULL;
            //n_ngb = tree_grav.getNeighborListOneParticle(pp[i], next);
            //for(PS::S32 j=0; j<n_ngb; j++){
            //    correctForceBetween2ParticlesInitial(pp, i, next[j], phii, phi_di, acci, acc_di, jerki, acc0i, NList);
            //}
            }
            */
        }

        if ( pp[i].neighbor.number ){
#pragma omp critical
            {
                NList.with_neighbor_list.push_back(i);
            }
        
            n_ngb_tot += pp[i].neighbor.number;
            n_with_ngb ++;
        }
        
        pp[i].acc   += acci;
        pp[i].phi   += phii;
        acc_d[i]     = acc_di;
        pp[i].phi_d  = phi_di;
        pp[i].jerk_d = jerki;
        pp[i].acc0   = ( acc0i > 0. ) ? pp[i].neighbor.number/acc0i : 0.;

#ifndef INTEGRATE_6TH_SUN
        pp[i].calcDeltatInitial();
#endif
    }

#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].acc_d =  acc_d[i];
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


 
