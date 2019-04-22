#pragma once

template <class TParticleI, class TParticleJ, class TForce>
void CalcForceLongEP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force);
template <class TParticleI, class TParticleJ, class TForce>
void CalcForceLongSP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force);
