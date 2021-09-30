/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "debug_utilities.h"
#include "user_defined.h"
#include "leapfrog.h"

/* Leapfrog integrators */
void InitialKick(PS::ParticleSystem<FP_nbody> & psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += 0.5 * dt * psys[i].acc;
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].vel.isnan() ||
            psys[i].vel.isinf() ||
            psys[i].acc.isnan() ||
            psys[i].acc.isinf()) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}
void InitialKick(PS::ParticleSystem<FP_star> & psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += 0.5 * dt * psys[i].acc;
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].vel.isnan() ||
            psys[i].vel.isinf() ||
            psys[i].acc.isnan() ||
            psys[i].acc.isinf()) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}
void InitialKick(PS::ParticleSystem<FP_gas> & psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel_half = psys[i].vel + 0.5 * dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng_half = psys[i].eng + 0.5 * dt * psys[i].eng_dot;
        psys[i].ent_half = psys[i].ent + 0.5 * dt * psys[i].ent_dot;
#endif
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].vel.isnan() || 
            psys[i].vel.isinf() ||
            psys[i].vel_half.isnan() ||
            psys[i].vel_half.isinf() ||
            psys[i].acc_grav.isnan() ||
            psys[i].acc_grav.isinf() ||
            psys[i].acc_hydro.isnan() ||
            psys[i].acc_hydro.isinf() ||
            std::isnan(psys[i].eng) ||
            std::isinf(psys[i].eng) ||
            (psys[i].eng <= 0.0)    ||
            std::isnan(psys[i].ent) ||
            std::isinf(psys[i].ent) ||
            (psys[i].ent <= 0.0)    ||
            std::isnan(psys[i].eng_half) ||
            std::isinf(psys[i].eng_half) ||
            (psys[i].eng_half <= 0.0) ||
            std::isnan(psys[i].ent_half) ||
            std::isinf(psys[i].ent_half) ||
            (psys[i].ent_half <= 0.0) ||
            std::isnan(psys[i].eng_dot) ||
            std::isinf(psys[i].eng_dot) ||
            std::isnan(psys[i].ent_dot) ||
            std::isinf(psys[i].ent_dot) ||
            (psys[i].ent_dot < 0.0)) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}


void FullDrift(PS::ParticleSystem<FP_nbody>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].pos += dt * psys[i].vel;
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].pos.isnan() ||
            psys[i].pos.isinf() ||
            psys[i].vel.isnan() ||
            psys[i].vel.isinf()) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}
void FullDrift(PS::ParticleSystem<FP_star>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].pos += dt * psys[i].vel;
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].pos.isnan() ||
            psys[i].pos.isinf() ||
            psys[i].vel.isnan() ||
            psys[i].vel.isinf()) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}
void FullDrift(PS::ParticleSystem<FP_gas>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].pos += dt * psys[i].vel_half;
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].pos.isnan() ||
            psys[i].pos.isinf() ||
            psys[i].vel_half.isnan() ||
            psys[i].vel_half.isinf()) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}


void Predict(PS::ParticleSystem<FP_gas>& psys, const PS::F64 dt){
    for(PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng += dt * psys[i].eng_dot;
        psys[i].ent += dt * psys[i].ent_dot;
#endif
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].vel.isnan() ||
            psys[i].vel.isinf() ||
            psys[i].acc_grav.isnan() ||
            psys[i].acc_grav.isinf() ||
            psys[i].acc_hydro.isnan() ||
            psys[i].acc_hydro.isinf() ||
            std::isnan(psys[i].eng) ||
            std::isinf(psys[i].eng) ||
            (psys[i].eng <= 0.0) ||
            std::isnan(psys[i].ent) ||
            std::isinf(psys[i].ent) ||
            (psys[i].ent <= 0.0) ||
            std::isnan(psys[i].eng_dot) ||
            std::isinf(psys[i].eng_dot) ||
            std::isnan(psys[i].ent_dot) ||
            std::isinf(psys[i].ent_dot) ||
            (psys[i].ent_dot < 0.0)) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}


void FinalKick(PS::ParticleSystem<FP_nbody>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += 0.5 * dt * psys[i].acc;
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].vel.isnan() ||
            psys[i].vel.isinf() ||
            psys[i].acc.isnan() ||
            psys[i].acc.isinf()) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}
void FinalKick(PS::ParticleSystem<FP_star>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += 0.5 * dt * psys[i].acc;
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].vel.isnan() ||
            psys[i].vel.isinf() ||
            psys[i].acc.isnan() ||
            psys[i].acc.isinf()) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}
void FinalKick(PS::ParticleSystem<FP_gas>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel = psys[i].vel_half + 0.5 * dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng = psys[i].eng_half + 0.5 * dt * psys[i].eng_dot;
        psys[i].ent = psys[i].ent_half + 0.5 * dt * psys[i].ent_dot;
#endif
        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (psys[i].vel.isnan() || 
            psys[i].vel.isinf() ||
            psys[i].vel_half.isnan() ||
            psys[i].vel_half.isinf() ||
            psys[i].acc_grav.isnan() ||
            psys[i].acc_grav.isinf() ||
            psys[i].acc_hydro.isnan() ||
            psys[i].acc_hydro.isinf() ||
            std::isnan(psys[i].eng) ||
            std::isinf(psys[i].eng) ||
            (psys[i].eng <= 0.0) ||
            std::isnan(psys[i].ent) ||
            std::isinf(psys[i].ent) ||
            (psys[i].ent <= 0.0) ||
            std::isnan(psys[i].eng_half) ||
            std::isinf(psys[i].eng_half) ||
            (psys[i].eng_half <= 0.0) ||
            std::isnan(psys[i].ent_half) ||
            std::isinf(psys[i].ent_half) ||
            (psys[i].ent_half <= 0.0) ||
            std::isnan(psys[i].eng_dot) ||
            std::isinf(psys[i].eng_dot) ||
            std::isnan(psys[i].ent_dot) ||
            std::isinf(psys[i].ent_dot) ||
            (psys[i].ent_dot < 0.0)) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif
    }
}
