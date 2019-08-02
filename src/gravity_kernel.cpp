#include <particle_simulator.hpp>

#ifdef CALC_EP_64bit
#define RSQRT_NR_EPJ_X4
#else
#define RSQRT_NR_EPJ_X2
#endif
#define RSQRT_NR_SPJ_X2
#define EP_OFFSET

#if defined(__AVX2__) || defined(__AVX512DQ__)
#include "phantomquad_for_p3t_x86.hpp"
#endif
#include "particle.h"

PS::S32 getAVXVersion(char * avx_var){
#ifdef __AVX512DQ__
    sprintf(avx_var, "AVX512DQ"); return 1;
#elif defined(__AVX2__)
    sprintf(avx_var, "AVX2"); return 1;
#else
    return 0;
#endif
}

PS::S32 getNIMAX(){
#if defined(__AVX2__) || defined(__AVX512DQ__)
    return PhantomGrapeQuad::NIMAX;
#else
    return std::numeric_limits<PS::S32>::max();
#endif
}
PS::S32 getNJMAX(){
#if defined(__AVX2__) || defined(__AVX512DQ__)
    return PhantomGrapeQuad::NJMAX;
#else
    return std::numeric_limits<PS::S32>::max();
#endif
}


#if defined(__AVX2__) || defined(__AVX512DQ__)

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPEP {
    //void CalcForceLongEPEP
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * ep_j,
                      const PS::S32 n_jp,
                      TForce * force){
        const PS::F64 eps2  = (*ep_i).getEps2();
        const PS::F64 gamma = (*ep_i).getGamma();

        //std::cout << "ni=" << n_ip << " nj=" << n_jp << std::endl;
        
        PhantomGrapeQuad pg;
        //static __thread PhantomGrapeQuad pg;
        //static thread_local PhantomGrapeQuad pg;
    
        pg.set_eps2(eps2);
#ifndef USE_INDIVIDUAL_CUTOFF
        const PS::F64 r_out = (*ep_i).getROut();
        pg.set_cutoff(r_out, gamma*r_out);
#endif
    
#ifdef EP_OFFSET
        const PS::F64vec pos_ofs = ep_i[0].getPos();
#endif

        for(PS::S32 i=0; i<n_ip; i++){
#ifdef EP_OFFSET
            const PS::F64vec pos_i = ep_i[i].getPos()-pos_ofs;
#else
            const PS::F64vec pos_i = ep_i[i].getPos();
#endif //EP_OFFSET
#ifdef USE_INDIVIDUAL_CUTOFF
            const PS::F64 r_outi = ep_i[i].getROut();
#ifdef CALC_EPEP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z, r_outi, r_outi*gamma);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, r_outi, r_outi*gamma);
#endif
#else //USE_INDIVIDUAL_CUTOFF
#ifdef CALC_EPEP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif
#endif //USE_INDIVIDUAL_CUTOFF
        }

        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it =ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){

                //const PS::F64 m_j = ep_j[i].getCharge();
                const PS::F64 m_j = ep_j[i].mass;
#ifdef EP_OFFSET
                const PS::F64vec pos_j = ep_j[i].getPos()-pos_ofs;
#else
                const PS::F64vec pos_j = ep_j[i].getPos();
#endif //EP_OFFSET
#ifdef USE_INDIVIDUAL_CUTOFF
                const PS::F64 r_outj = ep_j[i].getROut();
#ifdef CALC_EPEP_64bit
                pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, r_outj, r_outj*gamma);
#else
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, r_outj, r_outj*gamma);
#endif
#else //USE_INDIVIDUAL_CUTOFF
#ifdef CALC_EPEP_64bit
                pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#endif
#endif //USE_INDIVIDUAL_CUTOFF
            }
#ifdef CALC_EPEP_64bit
            pg.run_epj_for_p3t_with_linear_cutoff_d(n_ip, n_jp_tmp);
#else
            pg.run_epj_for_p3t_with_linear_cutoff(n_ip, n_jp_tmp);
#endif

            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].phi);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
#ifdef CALC_EPEP_64bit
                pg.accum_accp_one_d(i, a[0], a[1], a[2], *p);
#else
                pg.accum_accp_one(i, a[0], a[1], a[2], *p);
#endif
            }
        }
    }
};

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPSPQuad {
    //void CalcForceLongEPSPQuad
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * sp_j,
                      const PS::S32 n_jp,
                      TForce * force){
        const PS::F64 eps2  = (*ep_i).getEps2();
        const PS::F64 gamma = (*ep_i).getGamma();

        //std::cout << "ni=" << n_ip << " nj=" << n_jp << std::endl;
        
        //#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
        //#else
        //static __thread PhantomGrapeQuad pg;
        //static thread_local PhantomGrapeQuad pg;
        //#endif
    
#ifdef EP_OFFSET
        const PS::F64vec pos_ofs = ep_i[0].getPos();
#endif
    
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
#ifdef EP_OFFSET
            const PS::F64vec pos_i = ep_i[i].getPos()-pos_ofs;
#else        
            const PS::F64vec pos_i = ep_i[i].getPos();
#endif //EP_OFFSET
#ifdef USE_INDIVIDUAL_CUTOFF
            const PS::F64 r_outi = ep_i[i].getROut();
#ifdef CALC_EPEP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z, r_outi, r_outi*gamma);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, r_outi, r_outi*gamma);
#endif
#else //USE_INDIVIDUAL_CUTOFF
#ifdef CALC_EPEP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif
#endif //USE_INDIVIDUAL_CUTOFF
        }

        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it =ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
            
                const PS::F64 m_j = sp_j[i].getCharge();
#ifdef EP_OFFSET
                const PS::F64vec pos_j = sp_j[i].getPos()-pos_ofs;
#else
                const PS::F64vec pos_j = sp_j[i].getPos();
#endif //EP_OFFSET
                const PS::F64mat q = sp_j[i].quad;
#ifdef CALC_EPSP_64bit
                pg.set_spj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j,
                                 q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
#else
                pg.set_spj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j,
                               q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
#endif
            }
#ifdef CALC_EPSP_64bit
            pg.run_spj_d(n_ip, n_jp_tmp);
#else
            pg.run_spj(n_ip, n_jp_tmp);
#endif
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].phi);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
#ifdef CALC_EPSP_64bit
                pg.accum_accp_one_d(i, a[0], a[1], a[2], *p);
#else
                pg.accum_accp_one(i, a[0], a[1], a[2], *p);
#endif
                
            }
        }
    }
};

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPSPMono {
    //void CalcForceLongEPSPMono
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * sp_j,
                      const PS::S32 n_jp,
                      TForce * force){
        const PS::F64 eps2  = (*ep_i).getEps2();
        const PS::F64 gamma = (*ep_i).getGamma();

        //std::cout << "ni=" << n_ip << " nj=" << n_jp << std::endl;

        //#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
        //#else
        //static __thread PhantomGrapeQuad pg;
        //static thread_local PhantomGrapeQuad pg;
        //#endif
    
#ifdef EP_OFFSET
        const PS::F64vec pos_ofs = ep_i[0].getPos();
#endif
 
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
#ifdef EP_OFFSET
            const PS::F64vec pos_i = ep_i[i].getPos()-pos_ofs;
#else
            const PS::F64vec pos_i = ep_i[i].getPos();
#endif //EP_OFFSET
#ifdef USE_INDIVIDUAL_CUTOFF
            const PS::F64 r_outi = ep_i[i].getROut();
#ifdef CALC_EPEP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z, r_outi, r_outi*gamma);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z, r_outi, r_outi*gamma);
#endif
#else //USE_INDIVIDUAL_CUTOFF
#ifdef CALC_EPEP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif
#endif //USE_INDIVIDUAL_CUTOFF
        }

        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
            const PS::S32 it =ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                          
                const PS::F64 m_j = sp_j[i].getCharge();
#ifdef EP_OFFSET
                const PS::F64vec pos_j = sp_j[i].getPos()-pos_ofs;
#else
                const PS::F64vec pos_j = sp_j[i].getPos();
#endif //EP_OFFSET
#ifdef USE_INDIVIDUAL_CUTOFF
#ifdef CALC_EPSP_64bit
                pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, 0., 0.);
#else
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j, 0., 0.);
#endif
#else //USE_INDIVIDUAL_CUTOFF
#ifdef CALC_EPSP_64bit
                pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
                pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#endif
#endif //USE_INDIVIDUAL_CUTOFF
            }
#ifdef CALC_EPSP_64bit
            pg.run_epj_d(n_ip, n_jp_tmp);
#else
            pg.run_epj(n_ip, n_jp_tmp);
#endif
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].phi);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
#ifdef CALC_EPSP_64bit
                pg.accum_accp_one_d(i, a[0], a[1], a[2], *p);
#else
                pg.accum_accp_one(i, a[0], a[1], a[2], *p);
#endif
            }
        }
    }
};

#else //SIMD instruction not used

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPEP {
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * ep_j,
                      const PS::S32 n_jp,
                      TForce * force){
        const PS::F64 eps2  = (*ep_i).getEps2();
        const PS::F64 gamma = (*ep_i).getGamma();  
#ifndef USE_INDIVIDUAL_CUTOFF
        const PS::F64 r_out  = (*ep_i).getROut();
        const PS::F64 r_out2 = r_out * r_out;
#endif
        
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
#ifdef USE_INDIVIDUAL_CUTOFF
            const PS::F64 r_outi = ep_i[i].getROut();
#endif           
            PS::F64vec acc_i = 0.0;
            PS::F64    pot_i = 0.0;
            
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec rij = ep_j[j].getPos() - pos_i;
#ifdef USE_INDIVIDUAL_CUTOFF
                const PS::F64 r_out = std::max(r_outi, ep_j[j].getROut());
                const PS::F64 r_out2 = r_out * r_out;
#endif
                const PS::F64 r2_real = rij * rij + eps2;
                const PS::F64 r2      = std::max(r2_real, r_out2);
                
                const PS::F64 r_inv = 1./sqrt(r2);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                
                acc_i += m_r3 * rij;
                pot_i -= m_r;
            }
            force[i].acc += acc_i;
            force[i].phi += pot_i;
        }
    }
};

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPSPQuad {
    //void CalcForceLongEPSPQuad
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * sp_j,
                      const PS::S32 n_jp,
                      TForce * force){
        const PS::F64 eps2  = (*ep_i).getEps2();
        
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec pos_i = ep_i[ip].getPos();
            PS::F64vec acc_i = 0.;
            PS::F64    pot_i = 0.;
            
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64    mass_j = sp_j[jp].mass;
                PS::F64vec pos_j = sp_j[jp].getPos();
                PS::F64vec rij= pos_i - pos_j;
                
                PS::F64    r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64    tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr    = qr * rij;
                PS::F64 r_inv  = 1./sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                
                PS::F64 A = mass_j*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                
                acc_i -= A*rij + B*qr;
                pot_i -= mass_j*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += acc_i;
            force[ip].phi += pot_i;
        }
    }
};

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPSPMono {
    //void CalcForceLongEPSPMono
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * sp_j,
                      const PS::S32 n_jp,
                      TForce * force){
        const PS::F64 eps2  = (*ep_i).getEps2();
        const PS::F64 gamma = (*ep_i).getGamma();  
        
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();          
            PS::F64vec acc_i = 0.0;
            PS::F64    pot_i = 0.0;
            
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec rij = sp_j[j].getPos() - pos_i;
                const PS::F64    r2 = rij * rij + eps2;
                
                const PS::F64 r_inv = 1./sqrt(r2);
                const PS::F64 m_r = sp_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                
                acc_i += m_r3 * rij;
                pot_i -= m_r;
            }
            force[i].acc += acc_i;
            force[i].phi += pot_i;
        }
    }
};

#endif

template void CalcForceLongEPEP<EPGrav, EPGrav, ForceGrav>::operator ()
    (const EPGrav * ep_i, const PS::S32 n_ip, const EPGrav * ep_j, const PS::S32 n_jp, ForceGrav * force);
//void CalcForceLongEPEP<EPGrav, EPGrav, ForceGrav>

template void CalcForceLongEPSPQuad<EPGrav, PS::SPJQuadrupoleScatter, ForceGrav>::operator ()
    (const EPGrav * ep_i, const PS::S32 n_ip, const PS::SPJQuadrupoleScatter * ep_j, const PS::S32 n_jp, ForceGrav * force);
    //void CalcForceLongEPSPQuad<EPGrav, PS::SPJQuadrupoleScatter, ForceGrav>

template void CalcForceLongEPSPMono<EPGrav, PS::SPJMonopoleScatter, ForceGrav>::operator ()
    (const EPGrav * ep_i, const PS::S32 n_ip, const PS::SPJMonopoleScatter * ep_j, const PS::S32 n_jp, ForceGrav * force);
    //void CalcForceLongEPSPMono<EPGrav, PS::SPJMonopoleScatter, ForceGrav>

template void CalcForceLongEPSPQuad<EPGrav, PS::SPJQuadrupoleInAndOut, ForceGrav>::operator ()
    (const EPGrav * ep_i, const PS::S32 n_ip, const PS::SPJQuadrupoleInAndOut * ep_j, const PS::S32 n_jp, ForceGrav * force);
    //void CalcForceLongEPSPQuad<EPGrav, PS::SPJQuadrupoleInAndOut, ForceGrav>

template void CalcForceLongEPSPMono<EPGrav, PS::SPJMonopoleInAndOut, ForceGrav>::operator ()
    (const EPGrav * ep_i, const PS::S32 n_ip, const PS::SPJMonopoleInAndOut * ep_j, const PS::S32 n_jp, ForceGrav * force);
    //void CalcForceLongEPSPMono<EPGrav, PS::SPJMonopoleInAndOut, ForceGrav>
