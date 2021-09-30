
/* Leapfrog integrators */
void InitialKick(PS::ParticleSystem<FP_sph> & psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel_half = psys[i].vel + 0.5 * dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng_half = psys[i].eng + 0.5 * dt * psys[i].eng_dot;
        psys[i].ent_half = psys[i].ent + 0.5 * dt * psys[i].ent_dot;
#endif
    }
}
void InitialKick(PS::ParticleSystem<FP_sph> & psys,
                 const PS::F64 dt[]){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        const PS::S32 l = psys[i].level;
        if (l == 0) {
            psys[i].vel_half = psys[i].vel_save + 0.5 * dt[l] * (psys[i].acc_grav_save
                                                                +psys[i].acc_hydro_save);
#if !defined(ISOTHERMAL_EOS)
            psys[i].eng_half = psys[i].eng_save + 0.5 * dt[l] * psys[i].eng_dot_save;
            psys[i].ent_half = psys[i].ent_save + 0.5 * dt[l] * psys[i].ent_dot_save;
#endif
            // where dt[l] is either the time step at level 0 or
            // the time from the time when the previous time integration
            // was performed to the time when prediction is performed.
        } else if (l == 1) {
            psys[i].vel_half = psys[i].vel + 0.5 * dt[l] * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
            psys[i].eng_half = psys[i].eng + 0.5 * dt[l] * psys[i].eng_dot;
            psys[i].ent_half = psys[i].ent + 0.5 * dt[l] * psys[i].ent_dot;
#endif
            // where dt[l] must be the time step at level 1.
        }
    }
}


void FullDrift(PS::ParticleSystem<FP_sph>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].pos += dt * psys[i].vel_half;
    }
}
void FullDrift(PS::ParticleSystem<FP_sph>& psys,
               const PS::F64 dt[]){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        const PS::S32 l = psys[i].level;
        if (l == 0) {
            psys[i].pos = psys[i].pos_save + dt[l] * psys[i].vel_half;
            // where dt[l] is either the time step at level 0 or
            // the time from the time when the previous time integration
            // was performed to the time when prediction is performed.
        } else if (l == 1) {
            psys[i].pos += dt[l] * psys[i].vel_half;
            // where dt[l] must be the time step at level 1.
        }
    }
}


void Predict(PS::ParticleSystem<FP_sph>& psys, const PS::F64 dt){
    for(PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel += dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng += dt * psys[i].eng_dot;
        psys[i].ent += dt * psys[i].ent_dot;
#endif
    }
}
void Predict(PS::ParticleSystem<FP_sph>& psys,
             const PS::F64 dt[]){
    for(PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        const PS::S32 l = psys[i].level;
        if (l == 0) {
            psys[i].vel = psys[i].vel_save + dt[l] * (psys[i].acc_grav_save
                                                     +psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
            psys[i].eng = psys[i].eng_save + dt[l] * psys[i].eng_dot_save;
            psys[i].ent = psys[i].ent_save + dt[l] * psys[i].ent_dot_save;
#endif
        } else if (l == 1) {
            psys[i].vel += dt[l] * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
            psys[i].eng += dt[l] * psys[i].eng_dot;
            psys[i].ent += dt[l] * psys[i].ent_dot;
#endif
        }
    }
}


void FinalKick(PS::ParticleSystem<FP_sph>& psys, const PS::F64 dt){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        psys[i].vel = psys[i].vel_half + 0.5 * dt * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng = psys[i].eng_half + 0.5 * dt * psys[i].eng_dot;
        psys[i].ent = psys[i].ent_half + 0.5 * dt * psys[i].ent_dot;
#endif
    }
#if defined(DUMP_VELOCITY_OF_SPH_PARTICLE)
    const PS::F64 coeff = 0.1;
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].vel *= std::exp(- coeff * (CFL_hydro/0.1) * psys[i].snds * dt / psys[i].smth);
    }
#endif
}
void FinalKick(PS::ParticleSystem<FP_sph>& psys,
               const PS::F64 dt[]){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        const PS::S32 l = psys[i].level;
        psys[i].vel = psys[i].vel_half + 0.5 * dt[l] * (psys[i].acc_grav + psys[i].acc_hydro);
#if !defined(ISOTHERMAL_EOS)
        psys[i].eng = psys[i].eng_half + 0.5 * dt[l] * psys[i].eng_dot;
        psys[i].ent = psys[i].ent_half + 0.5 * dt[l] * psys[i].ent_dot;
#endif
    }
}
