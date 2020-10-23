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

template <class Tpsys>
void calcStarGravity(Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;
    
   PS::F64vec posi = pp.pos;
    PS::F64vec veli = pp.vel;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acci = pp.acc_;
#endif
    
    PS::F64vec dr = - posi;
    PS::F64vec dv = - veli;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec da = - acci;
#endif
    PS::F64    r2inv = 1. / (dr*dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    //PS::F64    r5inv = r3inv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;
#ifdef INTEGRATE_6TH_SUN
    PS::F64    beta  = (dv*dv + dr*da) * r2inv - 5. * alpha*alpha;
#endif

    pp.phi_s  = -m_sun * rinv;
    pp.acc_s  =  mj_rij3 * dr;
    pp.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
#ifdef INTEGRATE_6TH_SUN
    pp.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
#endif
}

template <class Tpsys>
void calcStarGravity_p(Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.xp;
    PS::F64vec veli = pp.vp;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acci = pp.ap;
#endif

    PS::F64vec dr = - posi;
    PS::F64vec dv = - veli;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec da = - acci;
#endif
    PS::F64    r2inv = 1. / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    //PS::F64    r5inv = r3inv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;
#ifdef INTEGRATE_6TH_SUN
    PS::F64    beta  = (dv*dv + dr*da) * r2inv - 5. * alpha*alpha;
#endif
    
    //pp.phi_s  = -m_sun * rinv;
    pp.acc_s  =  mj_rij3 * dr;
    pp.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
#ifdef INTEGRATE_6TH_SUN
    pp.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
#endif
}

template <class Tpsys>
void calcStarGravity_c(Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.pos;
    PS::F64vec veli = pp.vel;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acci = pp.acc_;
#endif

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec da = - acci;
#endif
    PS::F64    r2inv = 1. / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    //PS::F64    r5inv = r3inv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;
#ifdef INTEGRATE_6TH_SUN
    PS::F64    beta  = (dv*dv + dr*da) * r2inv - 5. * alpha*alpha;
#endif

    //pp.phi_s  = -m_sun * rinv;
    pp.acc_s  =  mj_rij3 * dr;
    pp.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
#ifdef INTEGRATE_6TH_SUN
    pp.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
#endif
}

template <class Tpsys>
void calcStarAccJerk(Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;
    
    PS::F64vec posi = pp.pos;
    PS::F64vec veli = pp.vel;
    
    PS::F64vec dr = - posi;
    PS::F64vec dv = - veli;
    PS::F64    r2inv = 1. / (dr*dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;

    pp.phi_s  = -m_sun * rinv;
    pp.acc_s  =  mj_rij3 * dr;
    pp.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
}

#ifdef INTEGRATE_6TH_SUN
template <class Tpsys>
void calcStarSnap(Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;
    
    PS::F64vec posi = pp.pos;
    PS::F64vec veli = pp.vel;
    PS::F64vec acci = pp.acc_;
    
    PS::F64vec dr = - posi;
    PS::F64vec dv = - veli;
    PS::F64vec da = - acci;
    PS::F64    r2inv = 1. / (dr*dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;
    PS::F64    beta  = (dv*dv + dr*da) * r2inv - 5. * alpha*alpha;

    pp.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
}
#endif

#if 0
template <class Tpsys>
void calcStarAcc(Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.pos;

    PS::F64vec dr  = - posi;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;

    pp.phi_s = -m_sun * rinv;
    pp.acc_s =  m_sun * r3inv * dr;
}
#endif

template <class Tpsys>
void calcStarJerk(Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    const PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.pos;
    PS::F64vec veli = pp.vel;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    //PS::F64    r5inv = r3inv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;
    
    pp.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
}

template <class Tp, class Tpsys>
void calcGravity(Tp & pi,
                 Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    
    //#ifndef INTEGRATE_6TH_SUN
    calcStarGravity(pi);
    //#else
    //calcStarAccJerk(pi);
    //pi.setAcc_();
    //calcStarSnap(pi);
    //#endif

    pi.phi_d  = 0.;
    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    
#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++)
#else
    for(PS::S32 j=0; j<pp.size(); j++)
#endif
        {
#ifndef FORDEBUG
            pj_id = pi.n_hard_list.at(j);
#else
            pj_id = j;
#endif  
            
            if ( pi.id == pp[pj_id].id ) continue;
        
            PS::F64vec xj = pp[pj_id].pos;
            PS::F64vec dr  = xj - xi;
            PS::F64 dr2 = dr * dr;
            assert( dr2 != 0.0 );
            dr2 += eps2;
        
            PS::F64vec vj   = pp[pj_id].vel;
            PS::F64vec dv   = vj - vi;
            PS::F64    massj = pp[pj_id].mass;
        
            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv  * rinv;
            PS::F64 r3inv = r2inv * rinv;
            //PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
            PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif

            PS::F64 mj_rij3 = massj * r3inv;
            PS::F64 alpha = (dr*dv) * r2inv;

            PS::F64 _W    = 1.-cutoff_W(rij, r_out_inv);
            PS::F64 _K    = 1.-cutoff_K(rij, r_out_inv);
            PS::F64 dKdt = cutoff_dKdt(rij, r_out_inv, alpha);
            PS::F64 alpha_c = alpha*_K;

            pi.phi_d  -= massj * rinv * _W;
            pi.acc_d  += mj_rij3 *   _K * dr;
            pi.jerk_d += mj_rij3 * ( _K * dv - (3.*alpha_c + dKdt) * dr );
        }
}
    
template <class Tp, class Tpsys>
void calcGravity_p(Tp & pi,
                   Tpsys & pp)
{
    const PS::F64 eps2  = EPGrav::eps2;
    
    calcStarGravity_p(pi);

    //pi.phi_d  = 0.;
    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    
    PS::S32 pj_id = 0;
    PS::F64vec xpi = pi.xp;
    PS::F64vec vpi = pi.vp;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++)
#else
    for(PS::S32 j=0; j<pp.size(); j++)
#endif
        {
#ifndef FORDEBUG
            pj_id = pi.n_hard_list.at(j);
#else
            pj_id = j;
#endif
            if ( pi.id == pp[pj_id].id ) continue;

            PS::F64vec xpj = pp[pj_id].xp;
            PS::F64vec dr  = xpj - xpi;
            PS::F64 dr2 = dr * dr;
            assert( dr2 != 0.0 );
            dr2 += eps2;

            PS::F64vec vpj   = pp[pj_id].vp;
            PS::F64vec dv    = vpj - vpi;
            PS::F64    massj = pp[pj_id].mass;

            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv  * rinv;
            PS::F64 r3inv = r2inv * rinv;
            //PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
            PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif
            
            PS::F64 mj_rij3 = massj * r3inv;
            PS::F64 alpha = (dr*dv) * r2inv;
            
            //PS::F64 _W    = 1.-cutoff_W(rij, r_out_inv);
            PS::F64 _K    = 1.-cutoff_K(rij, r_out_inv);
            PS::F64 dKdt = cutoff_dKdt(rij, r_out_inv, alpha);
            PS::F64 alpha_c = alpha*_K;

            //pi.phi_d  -= massj * rinv * _W;
            pi.acc_d  += mj_rij3 *   _K * dr;
            pi.jerk_d += mj_rij3 * ( _K * dv - (3.*alpha_c + dKdt) * dr );
        }
}

template <class Tp, class Tpsys>
void calcGravity_c(Tp & pi,
                   Tpsys & pp)
{
    //assert( pi.neighbor != 0 );
    
    const PS::F64 eps2  = EPGrav::eps2;
    
    calcStarGravity_c(pi);
    
    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++)
#else
    for(PS::S32 j=0; j<pp.size(); j++)
#endif
        {
#ifndef FORDEBUG
            pj_id = pi.n_hard_list.at(j);
#else
            pj_id = j;
#endif
            if ( pi.id == pp[pj_id].id ) continue;
        
            PS::F64vec xpj = pp[pj_id].xp;
            PS::F64vec dr  = xpj - xi;
            PS::F64 dr2 = dr * dr;
            assert( dr2 != 0.0 );
            dr2 += eps2;
        
            PS::F64vec vpj   = pp[pj_id].vp;
            PS::F64vec dv    = vpj - vi;
            PS::F64    massj = pp[pj_id].mass;

            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv  * rinv;
            PS::F64 r3inv = r2inv * rinv;
            //PS::F64 r5inv = r3inv * r2inv;
        
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
            PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif

            PS::F64 mj_rij3 = massj * r3inv;
            PS::F64 alpha = (dr*dv) * r2inv;
            
            //PS::F64 _W   = 1.-cutoff_W(rij, r_out_inv);
            PS::F64 _K   = 1.-cutoff_K(rij, r_out_inv);
            PS::F64 dKdt = cutoff_dKdt(rij, r_out_inv, alpha);
            PS::F64 alpha_c = alpha*_K;

            //pi.phi_d  -= massj * rinv * _W;
            pi.acc_d  += mj_rij3 *   _K * dr;
            pi.jerk_d += mj_rij3 * ( _K * dv - (3.*alpha_c + dKdt) * dr );
        }
}

template <class Tp, class Tpsys>
void calcJerk(Tp & pi,
              Tpsys & pp)
{   
    const PS::F64 eps2  = EPGrav::eps2;
    
    calcStarJerk(pi);

    pi.jerk  = 0.;
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++)
#else
    for(PS::S32 j=0; j<pp.size(); j++)
#endif
        {
#ifndef FORDEBUG
            pj_id = pi.n_hard_list.at(j);
#else
            pj_id = j;
#endif
            if ( pi.id == pp[pj_id].id ) continue;
        
            PS::F64vec xj = pp[pj_id].pos;
            PS::F64vec dr  = xj - xi;
            PS::F64 dr2 = dr * dr;
            assert( dr2 != 0.0 );
            dr2 += eps2;
        
            PS::F64vec vj   = pp[pj_id].vel;
            PS::F64vec dv    = vj - vi;
            PS::F64    massj = pp[pj_id].mass;

            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv  * rinv;
            PS::F64 r3inv = r2inv * rinv;
            //PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
            PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif

            PS::F64 mj_rij3 = massj * r3inv;
            PS::F64 alpha = (dr*dv) * r2inv;
            
            //PS::F64 _W    = 1.-cutoff_W(rij, r_out_inv);
            PS::F64 _K    = 1.-cutoff_K(rij, r_out_inv);
            PS::F64 dKdt = cutoff_dKdt(rij, r_out_inv, alpha);
            PS::F64 alpha_c = alpha*_K;

            pi.jerk_d += mj_rij3 * ( _K * dv - (3.*alpha_c + dKdt) * dr );
        }
}

template <class Tp, class Tpsys>
void calcAccJerk(Tp & pi,
                 Tpsys & pp)
{
    //assert( pi.neighbor != 0 );
    
    const PS::F64 eps2  = EPGrav::eps2;
    
    calcStarAccJerk(pi);

    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    pi.phi_d  = 0.;
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    
#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++)
#else
    for(PS::S32 j=0; j<pp.size(); j++)
#endif
        {
#ifndef FORDEBUG
            pj_id = pi.n_hard_list.at(j);
#else
            pj_id = j;
#endif  
            
            if ( pi.id == pp[pj_id].id ) continue;
            
            PS::F64vec xj = pp[pj_id].pos;
            PS::F64vec dr  = xj - xi;
            PS::F64 dr2 = dr * dr;
            assert( dr2 != 0.0 );
            dr2 += eps2;
        
            PS::F64vec vj   = pp[pj_id].vel;
            PS::F64vec dv    = vj - vi;
            PS::F64    massj = pp[pj_id].mass;
        
            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv  * rinv;
            PS::F64 r3inv = r2inv * rinv;
            //PS::F64 r5inv = r3inv * r2inv;
        
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
            PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif

            PS::F64 mj_rij3 = massj * r3inv;
            PS::F64 alpha = (dr*dv) * r2inv;
            
            PS::F64 _W    = 1.-cutoff_W(rij, r_out_inv);
            PS::F64 _K    = 1.-cutoff_K(rij, r_out_inv);
            PS::F64 dKdt = cutoff_dKdt(rij, r_out_inv, alpha);
            PS::F64 alpha_c = alpha*_K;

            pi.phi_d  -= massj * rinv * _W;
            pi.acc_d  += mj_rij3 *   _K * dr;
            pi.jerk_d += mj_rij3 * ( _K * dv - (3.*alpha_c + dKdt) * dr );
        }
}


template <class Tpsys, class Tptree>
void correctForceLong(Tpsys & pp,
                      Tptree & tree_grav,
                      NeighborList & NList,
                      PS::S32 & nei_dist,
                      PS::S32 & nei_tot_loc)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::F64 eps2 = EPGrav::eps2;
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

    nei_dist    = 0;
    nei_tot_loc = 0;
#pragma omp parallel for reduction (+:nei_dist, nei_tot_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec acci = 0.;
        PS::F64    phii = 0.;
        PS::F64    acc0i = 0.;
        PS::S32    j_id = 0;
        PS::S32    j_id_local = 0;
        PS::S32    j_rank = 0;
        PS::S32    neighbor = pp[i].neighbor - 1;
        PS::F64vec posi = pp[i].getPos();
        pp[i].neighbor = 0.;
        pp[i].id_cluster = pp[i].id;
#if defined(CORRECT_NEIGHBOR) || defined(CHECK_NEIGHBOR)
        PS::F64vec acc_di = 0.;
        PS::F64    phi_di = 0.;
#endif

        if ( neighbor == 0 ) {
#ifdef CHECK_NEIGHBOR
            PS::S32 nei_loc = 0;
            EPGrav* next = NULL;
            nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
            assert( nei_loc == 1 );
#endif
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = EPGrav::r_out_inv;
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
                EPGrav* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                assert( nei_loc < 3 );
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
                PS::F64 r_out_inv = EPGrav::r_out_inv;
                PS::F64 r_out     = EPGrav::r_out;
#endif
                phii += pp[i].mass * r_out_inv;
                
                PS::F64 massj  = pp[j_id_local].mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                r_out_inv = std::min(pp[i].r_out_inv, pp[j_id_local].r_out_inv);
                r_out     = std::max(pp[i].r_out,     pp[j_id_local].r_out);
                PS::F64 r_search = std::max(pp[i].r_search, pp[j_id_local].r_search);
#else
                PS::F64 r_search = EPGrav::r_search;
#endif
                
                PS::F64vec posj = pp[j_id_local].getPos();
                PS::F64vec dr   = posj - posi;
                PS::F64 dr2 = dr * dr;
                assert( dr2 != 0.0 );
                assert( j_id > -1 );
                dr2 += eps2;
                PS::F64 rij   = sqrt(dr2);
                
#ifdef USE_RE_SEARCH_NEIGHBOR
                PS::F64vec dv      = pp[j_id_local].vel   - pp[i].vel;
                PS::F64vec da      = (pp[j_id_local].acc_s - pp[i].acc_s) + (pp[j_id_local].acc_d - pp[i].acc_d);
                PS::F64    drdv    = dr * dv;
                PS::F64    dv2     = dv * dv;
                PS::F64    da2     = da * da;
                PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FPGrav::dt_tree);
                PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                //assert( dr2 >= dr2_min );
                
                PS::F64    r_crit   = FPGrav::R_search2 * r_out;
                PS::F64    v_crit_a = FPGrav::R_search3*0.5*FPGrav::dt_tree;
                if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 ) {
                    //if ( dr2_min < r_crit * r_crit ) {
#endif
                    if ( rij < r_search ) {
#ifdef TEST_PTCL
                        if ( r_out > 0. ){
#endif
                            NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                            acc0i += r_out * r_out / massj;
#ifdef TEST_PTCL
                        }
#endif
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
                EPGrav* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                
                for(PS::S32 j=0; j<nei_loc; j++){
                    j_id       = (next+j)->id;
                    j_id_local = (next+j)->id_local;
                    j_rank     = (next+j)->myrank;
                
                    PS::F64 massj  = (next+j)->mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                    PS::F64 r_out_inv = std::min(pp[i].r_out_inv, (next+j)->r_out_inv);
                    PS::F64 r_out     = std::max(pp[i].r_out,     (next+j)->r_out);
                    PS::F64 r_search  = std::max(pp[i].r_search,  (next+j)->r_search);
#else
                    PS::F64 r_out_inv = EPGrav::r_out_inv;
                    PS::F64 r_out     = EPGrav::r_out;
                    PS::F64 r_search  = EPGrav::r_search;
#endif
                    
                    if ( j_id == pp[i].id ) {
                        phii += massj * r_out_inv;
                        continue;
                    }
            
                    PS::F64vec posj = (next+j)->getPos();
                    PS::F64vec dr   = posj - posi;
                    PS::F64 dr2 = dr * dr;              
                    assert( dr2 != 0.0 );
                    assert( j_id > -1 );
                    dr2 += eps2;
                    PS::F64 rij   = sqrt(dr2);

#ifdef USE_RE_SEARCH_NEIGHBOR
                    PS::F64vec dv      = (next+j)->vel   - pp[i].vel;
                    PS::F64vec da      = ((next+j)->acc_s - pp[i].acc_s) + ((next+j)->acc_d - pp[i].acc_d);
                    PS::F64    drdv    = dr * dv;
                    PS::F64    dv2     = dv * dv;
                    PS::F64    da2     = da * da;
                    PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FPGrav::dt_tree);
                    PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                    //assert( dr2 >= dr2_min );

                    PS::F64    r_crit   = FPGrav::R_search2 * r_out;
                    PS::F64    v_crit_a = FPGrav::R_search3*0.5*FPGrav::dt_tree;
                    if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 ){
                        //if ( dr2_min < r_crit * r_crit ){
#endif
                        if ( rij < r_search ) {
#ifdef TEST_PTCL
                            if ( r_out > 0. ){
#endif
                                NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                                acc0i += r_out * r_out / massj;
#ifdef TEST_PTCL
                            }
#endif
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

        nei_dist += pp[i].neighbor;
        nei_tot_loc ++;

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
                             PS::S32 & nei_dist,
                             PS::S32 & nei_tot_loc)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::F64 eps2 = EPGrav::eps2;
    NList.initializeList(pp);

    PS::F64vec acc_s[n_loc];
    PS::F64vec acc_d[n_loc];
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        acc_s[i] = pp[i].acc_s;
        acc_d[i] = pp[i].acc_d;
    }

    nei_dist    = 0;
    nei_tot_loc = 0;
#pragma omp parallel for reduction (+:nei_dist, nei_tot_loc)
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
        PS::S32    j_id = 0;
        PS::S32    j_id_local = 0;
        PS::S32    j_rank = 0;
        PS::S32    neighbor = pp[i].neighbor - 1;
        PS::F64vec posi = pp[i].getPos();
        pp[i].neighbor = 0.;
        pp[i].id_cluster = pp[i].id;
        

        if ( neighbor == 0 ) {      
#ifdef CHECK_NEIGHBOR
            PS::S32 nei_loc = 0;
            EPGrav* next = NULL;
            nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
            assert( nei_loc == 1 );
#endif
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = pp[i].r_out_inv;
#else
            PS::F64 r_out_inv = EPGrav::r_out_inv;
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
                EPGrav* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                assert( nei_loc < 3 );
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
                PS::F64 r_out_inv = EPGrav::r_out_inv;
                PS::F64 r_out     = EPGrav::r_out;
#endif
                phii += pp[i].mass * r_out_inv;
                
                PS::F64 massj  = pp[j_id_local].mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                r_out_inv = std::min(pp[i].r_out_inv, pp[j_id_local].r_out_inv);
                r_out     = std::max(pp[i].r_out,     pp[j_id_local].r_out);
                PS::F64 r_search = std::max(pp[i].r_search, pp[j_id_local].r_search);
#else
                PS::F64 r_search = EPGrav::r_search;
#endif
                
                PS::F64vec posj = pp[j_id_local].getPos();
                PS::F64vec dr   = posj - posi;
                PS::F64 dr2 = dr * dr;
                assert( dr2 != 0.0 );
                assert( j_id > -1 );
                dr2 += eps2;
                PS::F64 rij   = sqrt(dr2);
                
                PS::F64vec dv      = pp[j_id_local].vel   - pp[i].vel;
                PS::F64    drdv    = dr * dv;
#ifdef USE_RE_SEARCH_NEIGHBOR
                PS::F64vec da      = (acc_s[j_id_local] - acc_s[i]) + (acc_d[j_id_local] - acc_d[i]);
                PS::F64    dv2     = dv * dv;
                PS::F64    da2     = da * da;
                PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FPGrav::dt_tree);
                PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                //assert( dr2 >= dr2_min );
                
                PS::F64    r_crit   = FPGrav::R_search2 * r_out;
                PS::F64    v_crit_a = FPGrav::R_search3*0.5*FPGrav::dt_tree;
                if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 || da2 == 0. ){
                //if ( dr2_min < r_crit * r_crit ){
#endif
                    if ( rij < r_search ) {
#ifdef TEST_PTCL
                        if ( r_out > 0. ){
#endif
                            NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                            acc0i += r_out * r_out / massj;
#ifdef TEST_PTCL
                        }
#endif
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
                EPGrav* next = NULL;
                nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
                
                for(PS::S32 j=0; j<nei_loc; j++){
                    j_id       = (next+j)->id;
                    j_id_local = (next+j)->id_local;
                    j_rank     = (next+j)->myrank;
                
                    PS::F64 massj  = (next+j)->mass;
#ifdef USE_INDIVIDUAL_CUTOFF
                    PS::F64 r_out_inv = std::min(pp[i].r_out_inv, (next+j)->r_out_inv);
                    PS::F64 r_out     = std::max(pp[i].r_out,     (next+j)->r_out);
                    PS::F64 r_search  = std::max(pp[i].r_search,  (next+j)->r_search);
#else
                    PS::F64 r_out_inv = EPGrav::r_out_inv;
                    PS::F64 r_out     = EPGrav::r_out;
                    PS::F64 r_search  = EPGrav::r_search;
#endif

                    if ( j_id == pp[i].id ) {
                        phii += massj * r_out_inv;
                        continue;
                    }
            
                    PS::F64vec posj = (next+j)->getPos();
                    PS::F64vec dr   = posj - posi;
                    PS::F64 dr2 = dr * dr;              
                    assert( dr2 != 0.0 );
                    assert( j_id > -1 );
                    dr2 += eps2;
                    PS::F64 rij   = sqrt(dr2);
                    
                    PS::F64vec dv      = (next+j)->vel   - pp[i].vel;
                    PS::F64    drdv    = dr * dv;
#ifdef USE_RE_SEARCH_NEIGHBOR
                    PS::F64vec da      = ((next+j)->acc_s - acc_s[i]) + ((next+j)->acc_d - acc_d[i]);
                    PS::F64    dv2     = dv * dv;
                    PS::F64    da2     = da * da;
                    PS::F64    t_min   = std::min(std::max(-drdv/sqrt(dv2), 0.), FPGrav::dt_tree);
                    PS::F64    dr2_min = std::min(dr2, dr2 + 2.*drdv*t_min + dv2*t_min*t_min);
                    //assert( dr2 >= dr2_min );
                    
                    PS::F64    r_crit   = FPGrav::R_search2 * r_out;
                    PS::F64    v_crit_a = FPGrav::R_search3*0.5*FPGrav::dt_tree;
                    if ( dr2_min < r_crit * r_crit || dv2 < v_crit_a * v_crit_a * da2 || da2 == 0. ){
                    //if ( dr2_min < r_crit * r_crit ){
#endif
                        if ( rij < r_search ) {
#ifdef TEST_PTCL
                            if ( r_out > 0. ){
#endif
                                NList.addNeighbor(pp, i, j_id, j_rank, j_id_local);
                                acc0i += r_out * r_out / massj;
#ifdef TEST_PTCL
                            }
#endif
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
        
        nei_dist += pp[i].neighbor;
        nei_tot_loc ++;
        
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
    PS::F64 eps2 = EPGrav::eps2;
    
    PS::F64vec acc_loc = 0;
    PS::F64vec pos_grav_loc = 0;
    PS::F64vec vel_grav_loc = 0;
    PS::F64    mass_tot_loc = 0;
    PS::F64    mass_tot_glb = 0;
#pragma omp parallel for reduction (+:acc_loc, pos_grav_loc, vel_grav_loc, mass_tot_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        
        PS::F64vec posi  = pp[i].getPos();
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
    
    FPGrav::acc_indirect = PS::Comm::getSum(acc_loc);
    
    FPGrav::mass_tot = mass_tot_glb = PS::Comm::getSum(mass_tot_loc) + FPGrav::m_sun;
    FPGrav::pos_g = PS::Comm::getSum(pos_grav_loc) / mass_tot_glb;
    FPGrav::vel_g = PS::Comm::getSum(vel_grav_loc) / mass_tot_glb;
}

template <class Tpsys>
PS::F64 calcIndirectEnergy(Tpsys & pp)
{
    PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::F64 eps2 = EPGrav::eps2;
    
    PS::F64vec vel_grav_loc = 0;
    PS::F64vec vel_grav_glb = 0;
    PS::F64    mass_tot_loc = 0;
    PS::F64    mass_tot_glb = 0;
#pragma omp parallel for reduction (+:pos_grav_loc, vel_grav_loc, mass_tot_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec posi  = pp[i].getPos();
        PS::F64vec veli  = pp[i].vel;
        PS::F64    massi = pp[i].mass;
        pos_grav_loc += massi * posi;
        vel_grav_loc += massi * veli;
        mass_tot_loc += massi;
    }
    mass_tot_glb = PS::Comm::getSum(mass_tot_loc) + FPGrav::m_sun;
    vel_grav_glb = PS::Comm::getSum(vel_grav_loc) / mass_tot_glb;
    
    return -0.5*(mass_tot_glb + FPGrav::m_sun) * vel_grav_glb*vel_grav_glb;
}

template <class Tpsys>
PS::F64 getIndirectEnergy(Tpsys & pp)
{
    return -0.5 * (FPGrav::mass_tot + FPGrav::m_sun) * FPGrav::vel_g * FPGrav::vel_g;
}
#endif


 
