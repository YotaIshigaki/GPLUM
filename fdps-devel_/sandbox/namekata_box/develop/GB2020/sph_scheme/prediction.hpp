void predictDensityEtc(PS::ParticleSystem<FP_sph> & psys,
                       const PS::F64 dt[]){
    // This function overwrite dens, smth, gradh, divv, rotv by
    // predicted values.
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        const PS::S32 l = psys[i].level;
        if (l == 0) {
            //std::cout << "id = " << psys[i].id
            //          << " r = " << std::sqrt(psys[i].pos * psys[i].pos)
            //          << std::endl;
            //if (psys[i].id == 42799) std::cout << "[bef.] dens = " << psys[i].dens << std::endl;
            psys[i].dens = psys[i].dens_save + dt[l] * psys[i].dens_dot_save;
            //if (psys[i].id == 42799) std::cout << "[aft.] dens = " << psys[i].dens << std::endl;
#if defined(USE_BALSARA_SWITCH)
            psys[i].divv = psys[i].divv_save + dt[l] * psys[i].divv_dot_save;
            psys[i].rotv = psys[i].rotv_save + dt[l] * psys[i].rotv_dot_save;
#endif
            // Predict the smoothing length using the constraint
            psys[i].smth = std::cbrt( (3.0 * mass_avg * N_neighbor)
                                    / (4.0 * math_const::pi * psys[i].dens));

            // Predict the grad-h term
            const PS::F64 drhodh = psys[i].drhodh_save + dt[l] * psys[i].drhodh_dot_save;
            psys[i].gradh = 1.0 / (1.0 + (psys[i].smth * drhodh) / (3.0 * psys[i].dens));

            // where dt[l] is either the time step at level 0 or
            // the time from the time when the previous time integration
            // was performed to the time when prediction is performed.
        }
    }
}

void saveDataForPrediction(PS::ParticleSystem<FP_sph> & psys) {
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++){
        const PS::S32 l = psys[i].level;
        if (l == 0) psys[i].save();
    }
}
