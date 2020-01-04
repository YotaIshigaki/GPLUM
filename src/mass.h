#pragma once

class MassGrowth{
public:
    static PS::F64 tau_mass;
    static PS::F64 mass_crit;
    static PS::F64 mass_max;
    std::vector<PS::S32> id_list;

    MassGrowth(){
        id_list.clear();
    }

    template <class Tpsys>
    PS::S32 makeList(Tpsys & pp,
                     PS::F64 time){
        const PS::S32 n_loc  = pp.getNumberOfParticleLocal();
        
        id_list.clear();
        for ( PS::S32 i=0; i<n_loc; i++ ){
            if ( mass_crit < pp[i].mass && mass_max > pp[i].mass ) {
                id_list.push_back(i);
            }
        }
    }

    template <class Tpsys>
    PS::S32 calcMassGrowth(Tpsys & pp,
                           PS::F64 time,
                           PS::F64 dt_tree){
        const PS::S32 n_loc  = pp.getNumberOfParticleLocal();
        const PS::S32 n_list = id_list.size();
        
#pragma omp parallel for
        for ( PS::S32 i=0; i<n_loc; i++ ){
            for ( PS::S32 j=0; j<n_list; j++ ){
                if ( pp[i] == id_list.at(j) ){
                    
                    if ( tau_mass != 0. )
                        pp[i].mass *= exp(dt_tree / tau_mass);
                    
                }
            }
        }
        return n_list;
    }
    
};

PS::F64 MassGrowth::tau_mass = 0.;
