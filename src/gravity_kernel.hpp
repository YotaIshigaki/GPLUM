#pragma once

PS::S32 getAVXVersion(char * avx_var);

PS::S32 getNIMAX();
PS::S32 getNJMAX();

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPEP {
    //void CalcForceLongEPEP
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * ep_j,
                      const PS::S32 n_jp,
                      TForce * force);
};

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPSPQuad {
    //void CalcForceLongEPSPQuad
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * sp_j,
                      const PS::S32 n_jp,
                      TForce * force);
};

template <class TParticleI, class TParticleJ, class TForce>
struct CalcForceLongEPSPMono {
    //void CalcForceLongEPSPMono
    void operator () (const TParticleI * ep_i,
                      const PS::S32 n_ip,
                      const TParticleJ * sp_j,
                      const PS::S32 n_jp,
                      TForce * force);
};
