/*
template <class T>
PS::F32ort getParticleDomain(const PS::ParticleSystem<T> & system)
{
    PS::F32ort domain;
    PS::S32 ndim = PS::DIMENSION;
    PS::S32 nloc = system.getNumberOfParticleLocal();

    for(PS::S32 k = 0; k < ndim; k++) {
        domain.low_[k]  = FLT_MAX;
        domain.high_[k] = - FLT_MAX;
    }

    for(PS::S32 i = 0; i < nloc; i++) {
        for(PS::S32 k = 0; k < ndim; k++) {
        if(system[i].pos[k] < domain.low_[k])
            domain.low_[k]  = system[i].pos[k];
        if(system[i].pos[k] > domain.high_[k])
            domain.high_[k] = system[i].pos[k];
        }
    }

    return domain;
}

PS::S32 getUniformDistributionFromArg1ToArg2(PS::S32 arg1,
                                             PS::S32 arg2)
{
    PS::S32 random_number;

    random_number = (PS::S32)((arg2 - arg1 + 1) * frand()) + arg1; // Function frand should be MT.

    return random_number;
}

template <class T>
void getSampleParticle(const PS::ParticleSystem<T> & system,
                       const PS::S32 target_number_of_sample_particle,
                       PS::S32 & nsample,
                       PS::F32vec pos_sample[])
{
    PS::S32 nglb = 0;
    PS::S32 nloc = system.getNumberOfParticleLocal();    

    MPI::COMM_WORLD.Allreduce(&nloc, &nglb, 1, MPI::INT, MPI::SUM);

    nsample = (PS::S32)(nloc * target_number_of_sample_particle / (PS::F32)nglb);
    nsample = (nsample < nloc) ? nsample : nloc;

    PS::S32 *record = new PS::S32[nsample];

    for(PS::S32 i = 0; i < nsample; i++) {
        PS::S32 j = getUniformDistributionFromArg1ToArg2(i, nloc-1);
        T hold   = system[j];
        system[j] = system[i];
        system[i] = hold;        
        record[i] = j;
    }

    for(PS::S32 i = 0; i < nsample; i++) {
        pos_sample[i] = system[i].pos;
    }

    for(PS::S32 i = nsample - 1; i >= 0; i--) {
        PS::S32 j = record[i];
        T hold   = system[j];
        system[j] = system[i];
        system[i] = hold;
    }

    delete record;

    return;
}
*/
