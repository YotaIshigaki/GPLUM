#pragma once

#include "cutfunc.h"

template <class Tpsys>
void velKick(Tpsys & pp){
    PS::S32 n = pp.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
        pp[i].velKick();
    }
}

template <class Tpsys>
void calcStarGravity(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;
    
    PS::F64vec posi = pp.getPos();
    PS::F64vec veli = pp.vel;
    
    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    //PS::F64    rinv  = rsqrt(dr * dr + eps2);
    //PS::F64    r2inv = rinv * rinv;
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.phi_s  = - m_sun * rinv;
    pp.acc_d += m_sun * r3inv * dr;
    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);

#ifdef CHECK_NEIGHBOR
    pp.true_neighbor = 0;
#endif
}

template <class Tpsys>
void calcStarGravity_p(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.xp;
    PS::F64vec veli = pp.vp;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    //PS::F64    rinv  = rsqrt(dr * dr + eps2);
    //PS::F64    r2inv = rinv * rinv;
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.acc_d += m_sun * r3inv * dr;
    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);
}

template <class Tpsys>
void calcStarGravity_c(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.getPos();
    PS::F64vec veli = pp.vel;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    //PS::F64    rinv  = rsqrt(dr * dr + eps2);
    //PS::F64    r2inv = rinv * rinv;
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.acc_d += m_sun * r3inv * dr;
    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);
}

#if 0
template <class Tpsys>
void calcStarAcc(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.getPos();

    PS::F64vec dr  = - posi;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    //PS::F64    rinv  = rsqrt(dr * dr + eps2);
    PS::F64    r3inv = rinv * r2inv;

    pp.phi_s  = - m_sun * rinv;
    pp.acc_d += m_sun * r3inv * dr;
}
#endif

template <class Tpsys>
void calcStarJerk(Tpsys & pp)
{
    PS::F64 eps2 = EPGrav::eps2;
    PS::F64 m_sun = FPGrav::m_sun;

    PS::F64vec posi = pp.getPos();
    PS::F64vec veli = pp.vel;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    //PS::F64    rinv  = rsqrt(dr * dr + eps2);
    //PS::F64    r2inv = rinv * rinv;
    PS::F64    r3inv = rinv * r2inv;
    PS::F64    r5inv = r3inv * r2inv;

    pp.jerk  += m_sun * (r3inv*dv -3.0*r5inv*(dr*dv)*dr);
}

template <class Tp, class Tpsys>
void calcGravity(Tp & pi,
                 Tpsys & pp)
{
    assert( pi.neighbor != 0 );
    
    PS::F64 eps2  = EPGrav::eps2;
    
    pi.phi_d = 0.0;
    pi.phi_s = 0.0;
    pi.acc_d = 0.0;
    pi.jerk  = 0.0;
    
    calcStarGravity(pi);
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    PS::F64vec acci = 0.0;
    PS::F64    phii = 0.0;
    PS::F64vec jerki= 0.0;
    
#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
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
        PS::F64    massj = pp[pj_id].getCharge();
        
        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 rinv  = 1. / rij;
        PS::F64 r2inv = rinv  * rinv;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;
        
#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
        PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif
        PS::F64 W  = cutoff_W(rij, r_out_inv);
        PS::F64 K  = cutoff_K(rij, r_out_inv);
        PS::F64 dK = cutoff_dK(rij, rinv, drdv, r_out_inv);
        
        phii  -= massj * rinv * (1.-W);
        acci  += massj * r3inv * dr * (1.-K);
        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );

#ifdef CHECK_NEIGHBOR
#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out = std::max(pi.r_out, pp[pj_id].r_out);
#else
        PS::F64 r_out = FPGrav::r_out;
#endif
        if ( rij < r_out ) pi.true_neighbor ++;
#endif
    }
    pi.acc_d += acci;
    pi.phi_d += phii;
    pi.jerk  += jerki;
}

template <class Tp, class Tpsys>
void calcGravity_p(Tp & pi,
                   Tpsys & pp)
{
    assert( pi.neighbor != 0 );

    PS::F64 eps2  = EPGrav::eps2;
  
    pi.acc_d = 0.0;
    pi.jerk  = 0.0;
    calcStarGravity_p(pi);
        
    PS::S32 pj_id = 0;
    PS::F64vec xpi = pi.xp;
    PS::F64vec vpi = pi.vp;
    PS::F64vec acci = 0.0;
    PS::F64vec jerki= 0.0;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
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
        PS::F64    massj = pp[pj_id].getCharge();

        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 rinv  = 1. / rij;
        PS::F64 r2inv = rinv  * rinv;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
        PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif
        PS::F64 K  = cutoff_K(rij, r_out_inv);
        PS::F64 dK = cutoff_dK(rij, rinv, drdv, r_out_inv);
        
        acci  += massj * r3inv * dr * (1.-K);
        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
    }
    pi.acc_d += acci;
    pi.jerk  += jerki;
}

template <class Tp, class Tpsys>
void calcGravity_c(Tp & pi,
                   Tpsys & pp)
{
    assert( pi.neighbor != 0 );
    
    PS::F64 eps2  = EPGrav::eps2;
    pi.acc_d = 0.0;
    pi.jerk  = 0.0;
    
    calcStarGravity_c(pi);
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    PS::F64vec acci = 0.0;
    PS::F64vec jerki= 0.0;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
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
        PS::F64    massj = pp[pj_id].getCharge();

        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 rinv  = 1. / rij;
        PS::F64 r2inv = rinv  * rinv;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;
        
#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
        PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif
        PS::F64 K  = cutoff_K(rij, r_out_inv);
        PS::F64 dK = cutoff_dK(rij, rinv, drdv, r_out_inv);
        
        acci  += massj * r3inv * dr * (1.-K);
        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
    }
    pi.acc_d += acci;
    pi.jerk  += jerki;
}

template <class Tp, class Tpsys>
void calcJerk(Tp & pi,
              Tpsys & pp)
{   
    PS::F64 eps2  = EPGrav::eps2;
    pi.jerk  = 0.0;
    
    calcStarJerk(pi);
    
    PS::S32 pj_id = 0;
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;
    PS::F64vec jerki= 0.0;

#ifndef FORDEBUG
    for(PS::S32 j=0; j<pi.neighbor; j++){
        pj_id = pi.n_hard_list.at(j);
#else
    for(PS::S32 j=0; j<pp.size(); j++){
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
        PS::F64    massj = pp[pj_id].getCharge();

        PS::F64 rij   = sqrt(dr2);
        PS::F64 drdv  = dr*dv;
        PS::F64 rinv  = 1. / rij;
        PS::F64 r2inv = rinv  * rinv;
        PS::F64 r3inv = r2inv * rinv;
        PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_RADII
        PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[pj_id].r_out_inv);
#else
        PS::F64 r_out_inv = FPGrav::r_out_inv;
#endif
        PS::F64 K  = cutoff_K(rij, r_out_inv);
        PS::F64 dK = cutoff_dK(rij, rinv, drdv, r_out_inv);

        jerki += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
    }
    pi.jerk  += jerki;
}



template <class Tpsys, class Tptree>
void correctForceLong(Tpsys & pp,
                      Tptree & tree_grav,
                      NeighborList & NList,
                      PS::S32 & nei_dist,
                      PS::S32 & nei_tot_loc)
{
    PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::F64 eps2 = EPGrav::eps2;
    NList.initializeList(pp);

    PS::S32 nei_dist_tmp    = 0;
    PS::S32 nei_tot_loc_tmp = 0;
#pragma omp parallel for reduction (+:nei_tot_loc_tmp, nei_dist_tmp)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::S32 j_id    = 0;
        PS::S32 j_rank  = 0;
        PS::S32 nei_loc = 0;
        EPGrav* next = NULL;
        nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
        pp[i].neighbor = 0;
        //pp[i].inDomain = true;
        pp[i].id_cluster = pp[i].id;

        PS::F64vec posi = pp[i].getPos();
        PS::F64vec acci= 0.0;
        PS::F64 phii = 0.0;
#ifdef CHECK_NEIGHBOR
        PS::S32 true_neighbor = 0;
#endif

        for(PS::S32 j=0; j<nei_loc; j++){
            j_id   = (next+j)->id;
            j_rank = (next+j)->myrank;

            PS::F64 massj  = (next+j)->getCharge();
#ifdef USE_INDIVIDUAL_RADII
            PS::F64 r_out_inv = std::min(pp[i].r_out_inv, (next+j)->r_out_inv);
#else
            PS::F64 r_out_inv = EPGrav::r_out_inv;
#endif
            
            if ( pp[i].id == j_id ) {
#ifdef CUTOFF_0
                phii  += massj * r_out_inv;
#endif
                continue;
            }
            
            PS::F64vec posj = (next+j)->getPos();
            PS::F64vec dr   = posj - posi;
            PS::F64 dr2 = dr * dr;
            
            assert( dr2 != 0.0 );
            assert( j_id > -1 );
            dr2 += eps2;

            NList.addNeighbor(pp, i, j_id, j_rank);
            
            
            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv * rinv;
            PS::F64 r3inv = rinv * r2inv;
            
            PS::F64 W  = cutoff_W(rij, r_out_inv);
            PS::F64 K  = cutoff_K(rij, r_out_inv);

#ifdef CUTOFF_0
            phii  -= massj * ( rinv * W  - r_out_inv );
            acci  += massj * ( r3inv * K - r_out_inv * r_out_inv * r_out_inv ) * dr;
#else
            PS::F64 K0  = std::min( std::max( EPGrav::g2_1_inv*(EPGrav::g2-dr2*r_out_inv*r_out_inv), 0.), 1.);

            phii  -= massj * rinv * ( W  - K0 );
            acci  += massj * r3inv * ( K - K0 ) * dr;
#endif

#ifdef CHECK_NEIGHBOR
#ifdef USE_INDIVIDUAL_RADII
            PS::F64 r_out = std::max(pp[i].r_out, (next+j)->r_out);
#else
            PS::F64 r_out = EPGrav::r_out;
#endif
            if ( rij < r_out ) true_neighbor ++;
#endif
        }
        
        nei_tot_loc_tmp += pp[i].neighbor;
        nei_dist_tmp ++;

        pp[i].acc  += acci;
        pp[i].phi  += phii;

#ifdef CHECK_NEIGHBOR
        if ( pp[i].true_neighbor < true_neighbor ) {
            std::cerr << "Unknown Neighbor is Exist on Particle " << pp[i].id
                      << " at Time "  << pp[i].time << std::endl;
        }
#endif
    }

    nei_dist    = nei_tot_loc_tmp;
    nei_tot_loc = nei_dist_tmp;
}
 
template <class Tpsys, class Tptree>
void correctForceLongInitial(Tpsys & pp,
                             Tptree & tree_grav,
                             NeighborList & NList,
                             PS::S32 & nei_dist,
                             PS::S32 & nei_tot_loc)
{
    PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::F64 eps2 = EPGrav::eps2;
    NList.initializeList(pp);

    PS::S32 nei_dist_tmp    = 0;
    PS::S32 nei_tot_loc_tmp = 0;
#pragma omp parallel for reduction (+:nei_tot_loc_tmp, nei_dist_tmp)
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].acc_d = 0.0;
        pp[i].jerk = 0.0;
        pp[i].phi_d = 0.0;
        pp[i].phi_s = 0.0;

        calcStarGravity(pp[i]);
        
        PS::S32 j_id = 0;
        PS::S32 j_rank = 0;
        PS::S32 nei_loc = 0;
        EPGrav* next = NULL;//Pointer Of Neighbor List
        nei_loc = tree_grav.getNeighborListOneParticle(pp[i], next);
        pp[i].neighbor = 0;
        //pp[i].inDomain = true;
        pp[i].id_cluster = pp[i].id;

        PS::F64vec posi = pp[i].getPos();
        PS::F64vec veli = pp[i].vel;
        PS::F64vec acci = 0.0;
        PS::F64vec acc_di = 0.0;
        PS::F64vec jerki  = 0.0;
        PS::F64 phii   = 0.0;
        PS::F64 phi_di = 0.0;

        for(PS::S32 j=0; j<nei_loc; j++){
            j_id   = (next+j)->id;
            j_rank = (next+j)->myrank;

            PS::F64 massj  = (next+j)->getCharge();
#ifdef USE_INDIVIDUAL_RADII
            PS::F64 r_out_inv = std::min(pp[i].r_out_inv, (next+j)->r_out_inv);
#else
            PS::F64 r_out_inv = EPGrav::r_out_inv;
#endif
            
            if ( pp[i].id == j_id ) {
#ifdef CUTOFF_0
                phii  += massj * r_out_inv;
#endif
                continue;
            }
            
            PS::F64vec posj   = (next+j)->getPos();
            PS::F64vec dr  = posj - posi;
            PS::F64 dr2 = dr * dr;
            
            assert( dr2 != 0.0 );
            assert( j_id > -1 );
            dr2 += eps2;
            
            NList.addNeighbor(pp, i, j_id, j_rank);
            
            PS::F64vec velj   = (next+j)->vel;
            PS::F64vec dv  = velj - veli;
            PS::F64 drdv  = dr*dv;

            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv * rinv;
            PS::F64 r3inv = rinv * r2inv;
            PS::F64 r5inv = r3inv * r2inv;
            
            PS::F64 W  = cutoff_W(rij, r_out_inv);
            PS::F64 K  = cutoff_K(rij, r_out_inv);
            PS::F64 dK = cutoff_dK(rij, rinv, drdv, r_out_inv);

#ifdef CUTOFF_0
            phii   -= massj * ( rinv * W - r_out_inv );
            phi_di -= massj * rinv * (1.-W);
            acci   += massj * ( r3inv * K - r_out_inv * r_out_inv * r_out_inv ) * dr;
            acc_di += massj * r3inv * dr * (1.-K);
            jerki  += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
#else
            PS::F64 K0  = std::min( std::max( EPGrav::g2_1_inv*(EPGrav::g2-dr2*r_out_inv*r_out_inv), 0.), 1.);
            
            phii   -= massj * rinv * ( W  - K0 );
            phi_di -= massj * rinv * (1.-W);
            acci   += massj * r3inv * ( K - K0 ) * dr;
            acc_di += massj * r3inv * dr * (1.-K);
            jerki  += massj * ( (r3inv*dv -3.*r5inv*(drdv)*dr)*(1.-K) - r3inv*dr*dK );
#endif
        }
        
        nei_tot_loc_tmp += pp[i].neighbor;
        nei_dist_tmp ++;
        
        pp[i].acc  += acci;
        pp[i].phi  += phii;
        pp[i].acc_d += acc_di;
        pp[i].phi_d += phi_di;
        pp[i].jerk += jerki;

        pp[i].calcDeltatInitial();
    }

    nei_dist    = nei_tot_loc_tmp;
    nei_tot_loc = nei_dist_tmp;
}
