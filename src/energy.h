#pragma once

class Energy{
public:
    PS::F64 etot;
    PS::F64 ekin;
    PS::F64 ephi_sun;
    PS::F64 ephi_planet;
    PS::F64 ephi;
    PS::F64 ephi_d;
    PS::F64 edisp;

    Energy(){ etot = ekin = ephi_sun = ephi_planet = ephi = ephi_d = edisp = 0.; }
    
    Energy(PS::F64 ekin0,
           PS::F64 ephi_sun0,
           PS::F64 ephi0,
           PS::F64 ephi_d0,
           PS::F64 edisp0){
        ekin = ekin0;
        ephi_sun = ephi_sun0;
        ephi = ephi0;
        ephi_d = ephi_d0;
        edisp = edisp0;
        ephi_planet = ephi + ephi_d;
        etot = ekin + ephi_sun + ephi_planet;
    }

    template<class Tpsys>
    void calcEnergy(const Tpsys & pp,
                    const bool clear=true){
        if ( clear ) etot = ekin = ephi_sun = ephi_planet = ephi = ephi_d = 0.0;

        PS::F64 ekin_loc = 0.0;
        PS::F64 ephi_sun_loc = 0.0;
        PS::F64 ephi_loc = 0.0;
        PS::F64 ephi_d_loc = 0.0;
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < n_loc; i++){
            ekin_loc     += pp[i].mass * pp[i].vel * pp[i].vel;
            ephi_sun_loc += pp[i].mass * pp[i].phi_s;
            ephi_loc     += pp[i].mass * pp[i].phi;
            ephi_d_loc   += pp[i].mass * pp[i].phi_d;   
        }
        ekin_loc *= 0.5;
        ephi_loc *= 0.5;
        ephi_d_loc *= 0.5;
        
        ekin     += PS::Comm::getSum(ekin_loc);
        ephi_sun += PS::Comm::getSum(ephi_sun_loc);
        ephi     += PS::Comm::getSum(ephi_loc);
        ephi_d   += PS::Comm::getSum(ephi_d_loc);
        ephi_planet =  ephi + ephi_d;
        etot = ekin + ephi_sun + ephi_planet;
    }

    PS::F64 calcEnergyError(const Energy e_init){
        return (etot - e_init.etot - edisp)/e_init.etot;
    }
};

class FileHeader{
public:
    PS::S32 n_body;
    PS::F64 time;
    //PS::F64 etot0;
    //PS::F64 etot1;
    //PS::F64 edisp;
    Energy e_init;
    Energy e_now;
    PS::S32 readAscii(FILE * fp) {
        fscanf(fp, "%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &time, &n_body,
               &e_init.etot, &e_init.ekin, &e_init.ephi_sun, &e_init.ephi_planet, &e_init.edisp,
               &e_now.etot,  &e_now.ekin,  &e_now.ephi_sun,  &e_now.ephi_planet,  &e_now.edisp);
		return n_body;
    }
    void writeAscii(FILE* fp) const {
        fprintf(fp,
                "%g\t%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\n",
                time, n_body,
                e_init.etot, e_init.ekin, e_init.ephi_sun, e_init.ephi_planet, e_init.edisp,
                e_now.etot,  e_now.ekin,  e_now.ephi_sun,  e_now.ephi_planet,  e_now.edisp);
    }

    FileHeader(){ n_body = 0; time = 0.; }
    FileHeader(PS::S32 n_body0,
               PS::F64 time0,
               Energy e_init0,
               Energy e_now0){
        n_body = n_body0;
        time = time0;
        e_init = e_init0;
        e_now = e_now0;
    }
};

#ifdef OUTPUT_DETAIL
template<class Tpsys>
void calcKineticEnergy(const Tpsys & pp,
                       PS::F64 & ekin)
{
    ekin = 0.;   
    PS::F64 ekin_loc = 0.0;    
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    
    for(PS::S32 i = 0; i < n_loc; i++){
        ekin_loc += pp[i].mass * pp[i].vel * pp[i].vel;
    }
    ekin_loc *= 0.5;
    
    ekin = PS::Comm::getSum(ekin_loc);
}
#endif


