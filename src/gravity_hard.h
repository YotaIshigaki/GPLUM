#pragma once

#include "cutfunc.h"

template <class Tpsys>
void calcStarGravity(Tpsys & pi)
{
    const PS::F64 eps2  = FP_t::eps2_sun;
    const PS::F64 m_sun = FP_t::m_sun;
    
   PS::F64vec posi = pi.pos;
    PS::F64vec veli = pi.vel;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acci = pi.acc_;
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

    pi.phi_s  = -m_sun * rinv;
    pi.acc_s  =  mj_rij3 * dr;
    pi.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
#ifdef INTEGRATE_6TH_SUN
    pi.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
#endif
}

template <class Tpsys>
void calcStarGravity_p(Tpsys & pi)
{
    const PS::F64 eps2  = FP_t::eps2_sun;
    const PS::F64 m_sun = FP_t::m_sun;

    PS::F64vec posi = pi.xp;
    PS::F64vec veli = pi.vp;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acci = pi.ap;
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
    
    //pi.phi_s  = -m_sun * rinv;
    pi.acc_s  =  mj_rij3 * dr;
    pi.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
#ifdef INTEGRATE_6TH_SUN
    pi.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
#endif
}

template <class Tpsys>
void calcStarGravity_c(Tpsys & pi)
{
    const PS::F64 eps2  = FP_t::eps2_sun;
    const PS::F64 m_sun = FP_t::m_sun;

    PS::F64vec posi = pi.pos;
    PS::F64vec veli = pi.vel;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acci = pi.acc_;
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

    //pi.phi_s  = -m_sun * rinv;
    pi.acc_s  =  mj_rij3 * dr;
    pi.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
#ifdef INTEGRATE_6TH_SUN
    pi.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
#endif
}

template <class Tpsys>
void calcStarAccJerk(Tpsys & pi)
{
    const PS::F64 eps2  = FP_t::eps2_sun;
    const PS::F64 m_sun = FP_t::m_sun;
    
    PS::F64vec posi = pi.pos;
    PS::F64vec veli = pi.vel;
    
    PS::F64vec dr = - posi;
    PS::F64vec dv = - veli;
    PS::F64    r2inv = 1. / (dr*dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;

    pi.phi_s  = -m_sun * rinv;
    pi.acc_s  =  mj_rij3 * dr;
    pi.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
}

#ifdef INTEGRATE_6TH_SUN
template <class Tpsys>
void calcStarSnap(Tpsys & pi)
{
    const PS::F64 eps2  = FP_t::eps2_sun;
    const PS::F64 m_sun = FP_t::m_sun;
    
    PS::F64vec posi = pi.pos;
    PS::F64vec veli = pi.vel;
    PS::F64vec acci = pi.acc_;
    
    PS::F64vec dr = - posi;
    PS::F64vec dv = - veli;
    PS::F64vec da = - acci;
    PS::F64    r2inv = 1. / (dr*dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;
    PS::F64    beta  = (dv*dv + dr*da) * r2inv - 5. * alpha*alpha;

    pi.snap_s =  mj_rij3 * (da - 6.*alpha * dv - 3.*beta * dr);
}
#endif

#if 0
template <class Tpsys>
void calcStarAcc(Tpsys & pi)
{
    const PS::F64 eps2  = FP_t::eps2_sun;
    const PS::F64 m_sun = FP_t::m_sun;

    PS::F64vec posi = pi.pos;

    PS::F64vec dr  = - posi;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;

    pi.phi_s = -m_sun * rinv;
    pi.acc_s =  m_sun * r3inv * dr;
}
#endif

template <class Tpsys>
void calcStarJerk(Tpsys & pi)
{
    const PS::F64 eps2  = FP_t::eps2_sun;
    const PS::F64 m_sun = FP_t::m_sun;

    PS::F64vec posi = pi.pos;
    PS::F64vec veli = pi.vel;

    PS::F64vec dr  = - posi;
    PS::F64vec dv  = - veli;
    PS::F64    r2inv = 1.0 / (dr * dr + eps2);
    PS::F64    rinv  = sqrt(r2inv);
    PS::F64    r3inv = rinv * r2inv;
    //PS::F64    r5inv = r3inv * r2inv;

    PS::F64    mj_rij3 = m_sun * r3inv;
    PS::F64    alpha = (dr*dv) * r2inv;
    
    pi.jerk_s =  mj_rij3 * (dv - 3.*alpha * dr);
}

template <class Tp, class Tpsys>
void calcGravity(Tp & pi,
                 Tpsys & pp)
{
    const PS::F64 eps2  = FP_t::eps2;
    
    calcStarGravity(pi);

    pi.phi_d  = 0.;
    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;

    const PS::S32 np = pp.size();
    for(PS::S32 j=0; j<np; j++) {
            if ( pi.id == pp[j].id ) continue;
        
            PS::F64vec xj = pp[j].pos;
            PS::F64vec dr  = xj - xi;
            PS::F64 dr2 = dr * dr;
            assert( dr2 != 0.0 );
            dr2 += eps2;
        
            PS::F64vec vj   = pp[j].vel;
            PS::F64vec dv   = vj - vi;
            PS::F64    massj = pp[j].mass;
        
            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv  * rinv;
            PS::F64 r3inv = r2inv * rinv;
            //PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[j].r_out_inv);
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
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
    const PS::F64 eps2  = FP_t::eps2;
    
    calcStarGravity_p(pi);

    //pi.phi_d  = 0.;
    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    
    PS::F64vec xpi = pi.xp;
    PS::F64vec vpi = pi.vp;

    const PS::S32 np = pp.size();
    for(PS::S32 j=0; j<np; j++) {
            if ( pi.id == pp[j].id ) continue;

            PS::F64vec xpj = pp[j].xp;
            PS::F64vec dr  = xpj - xpi;
            PS::F64 dr2 = dr * dr;
            assert( dr2 != 0.0 );
            dr2 += eps2;

            PS::F64vec vpj   = pp[j].vp;
            PS::F64vec dv    = vpj - vpi;
            PS::F64    massj = pp[j].mass;

            PS::F64 rij   = sqrt(dr2);
            PS::F64 rinv  = 1. / rij;
            PS::F64 r2inv = rinv  * rinv;
            PS::F64 r3inv = r2inv * rinv;
            //PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[j].r_out_inv);
#else
            PS::F64 r_out_inv = FP_t::r_out_inv;
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
    
    const PS::F64 eps2  = FP_t::eps2;
    
    calcStarGravity_c(pi);

    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;

    const PS::S32 np = pp.size();
    for(PS::S32 j=0; j<np; j++) {
        if ( pi.id == pp[j].id ) continue;
        
        PS::F64vec xpj = pp[j].xp;
        PS::F64vec dr  = xpj - xi;
        PS::F64 dr2 = dr * dr;
        assert( dr2 != 0.0 );
        dr2 += eps2;
        
        PS::F64vec vpj   = pp[j].vp;
        PS::F64vec dv    = vpj - vi;
        PS::F64    massj = pp[j].mass;
        
        PS::F64 rij   = sqrt(dr2);
        PS::F64 rinv  = 1. / rij;
        PS::F64 r2inv = rinv  * rinv;
        PS::F64 r3inv = r2inv * rinv;
        //PS::F64 r5inv = r3inv * r2inv;
        
#ifdef USE_INDIVIDUAL_CUTOFF
        PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[j].r_out_inv);
#else
        PS::F64 r_out_inv = FP_t::r_out_inv;
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
    const PS::F64 eps2  = FP_t::eps2;
    
    calcStarJerk(pi);

    pi.jerk  = 0.;
    
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;

    const PS::S32 np = pp.size();
    for(PS::S32 j=0; j<np; j++){
        if ( pi.id == pp[j].id ) continue;
        
        PS::F64vec xj = pp[j].pos;
        PS::F64vec dr  = xj - xi;
        PS::F64 dr2 = dr * dr;
        assert( dr2 != 0.0 );
        dr2 += eps2;
        
        PS::F64vec vj   = pp[j].vel;
        PS::F64vec dv    = vj - vi;
        PS::F64    massj = pp[j].mass;

        PS::F64 rij   = sqrt(dr2);
        PS::F64 rinv  = 1. / rij;
        PS::F64 r2inv = rinv  * rinv;
        PS::F64 r3inv = r2inv * rinv;
        //PS::F64 r5inv = r3inv * r2inv;

#ifdef USE_INDIVIDUAL_CUTOFF
        PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[j].r_out_inv);
#else
        PS::F64 r_out_inv = FP_t::r_out_inv;
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
    
    const PS::F64 eps2  = FP_t::eps2;
    
    calcStarAccJerk(pi);

    pi.acc_d  = 0.;
    pi.jerk_d = 0.;
    pi.phi_d  = 0.;
    
    PS::F64vec xi = pi.pos;
    PS::F64vec vi = pi.vel;

    const PS::S32 np = pp.size();
    for(PS::S32 j=0; j<np; j++) {
        if ( pi.id == pp[j].id ) continue;
            
        PS::F64vec xj = pp[j].pos;
        PS::F64vec dr  = xj - xi;
        PS::F64 dr2 = dr * dr;
        assert( dr2 != 0.0 );
        dr2 += eps2;
        
        PS::F64vec vj   = pp[j].vel;
        PS::F64vec dv    = vj - vi;
        PS::F64    massj = pp[j].mass;
        
        PS::F64 rij   = sqrt(dr2);
        PS::F64 rinv  = 1. / rij;
        PS::F64 r2inv = rinv  * rinv;
        PS::F64 r3inv = r2inv * rinv;
        //PS::F64 r5inv = r3inv * r2inv;
        
#ifdef USE_INDIVIDUAL_CUTOFF
        PS::F64 r_out_inv = std::min(pi.r_out_inv, pp[j].r_out_inv);
#else
        PS::F64 r_out_inv = FP_t::r_out_inv;
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


 
