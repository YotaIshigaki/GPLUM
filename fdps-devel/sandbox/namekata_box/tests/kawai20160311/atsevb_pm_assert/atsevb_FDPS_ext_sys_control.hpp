//***************************************************************************************
//  This is external system control routine.
//    This code is used by "atsevb_main.cpp"
//***************************************************************************************
#pragma once

namespace EXT_SYS_CONTROL {
    
    class EnergyBuf {
      private:
        //--- singeton rule
        EnergyBuf() = default;
        ~EnergyBuf() = default;
        
      
        //--- total of system
        PS::F64 eng_tot;
        PS::F64 eng_kin;
        PS::F64 eng_pot;
        
        //--- LJ
        PS::F64 pot_LJ;
        //--- coulomb
        PS::F64 pot_coulomb;
        //--- intra
        PS::F64 pot_intra;
        
        //--- reference value
        PS::F64 ref_eng_kin;
        PS::F64 ref_eng_pot;
        PS::F64 ref_eng_tot;
        PS::F64 ref_pot_LJ;
        PS::F64 ref_pot_coulomb;
        PS::F64 ref_pot_intra;
        
      public:
        //--- unable object copy
        EnergyBuf(const EnergyBuf&) = delete;
        EnergyBuf& operator=(const EnergyBuf&) = delete;
        EnergyBuf(EnergyBuf&&) = delete;
        EnergyBuf& operator=(EnergyBuf&&) = delete;
        
        //--- call singleton object
        static EnergyBuf* getInstance(void) {
            static EnergyBuf singleton;
            return &singleton;
        }
        
        //--- access
        void clear(){
            this->eng_tot = 0.0;
            this->eng_kin = 0.0;
            this->eng_pot = 0.0;
            this->pot_LJ      = 0.0;
            this->pot_coulomb = 0.0;
            this->pot_intra   = 0.0;
        }
        void setEngKin(const PS::F64 kin){     this->eng_kin = kin; }
        void setPotLJ(const PS::F64 pot){      this->pot_LJ = pot; }
        void setPotCoulomb(const PS::F64 pot){ this->pot_coulomb = pot; }
        void setPotIntra(const PS::F64 pot){   this->pot_intra = pot; }
        void calcTotal(){
            this->eng_pot =   this->pot_LJ
                            + this->pot_coulomb
                            + this->pot_intra;
            this->eng_tot = this->eng_kin + this->eng_pot;
        }
        void setRef(){
            this->ref_eng_tot = this->eng_tot;
            this->ref_eng_kin = this->eng_kin;
            this->ref_eng_pot = this->eng_pot;
            this->ref_pot_LJ      = this->pot_LJ;
            this->ref_pot_coulomb = this->pot_coulomb;
            this->ref_pot_intra   = this->pot_intra;
        }
        PS::F64 getEng(){    return this->eng_tot - ref_eng_tot; }
        PS::F64 getEngKin(){ return this->eng_kin - ref_eng_kin; }
        PS::F64 getEngPot(){ return this->eng_pot - ref_eng_pot; }
        PS::F64 getPotLJ(){      return this->pot_LJ - ref_pot_LJ; }
        PS::F64 getPotCoulomb(){ return this->pot_coulomb - ref_pot_coulomb; }
        PS::F64 getPotIntra(){   return this->pot_intra - ref_pot_intra; }
        
        PS::F64 getR_error(){
            return (this->eng_tot - this->ref_eng_tot)/(this->ref_eng_tot);
        }
    };
}

template<class Tpsys>
void CalcVirial(const Tpsys & system,
                PS::F64vec & virt_tot){
    PS::F64vec virt_loc = 0.0;
               virt_tot = 0.0;
    
    //--- sumation in local
    PS::S32 n_loc = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        virt_loc += system[i].getVirt();
    }
    
    //--- calculate total sumation
    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
      virt_tot = PS::Comm::getSum(virt_loc);
    #else
      virt_tot = virt_loc;
    #endif
}

template<class Tpsys>
void calcEnergy(const Tpsys & system){
    
    EXT_SYS_CONTROL::EnergyBuf *energy_buf = EXT_SYS_CONTROL::EnergyBuf::getInstance();
    energy_buf->clear();
    
    PS::F64 eng_kin_loc = 0.0;
    
    PS::F64 pot_LJ_loc      = 0.0;
    PS::F64 pot_coulomb_loc = 0.0;
    PS::F64 pot_intra_loc   = 0.0;
    
    //--- sumation in local
    PS::S32 n_loc = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        //--- kinetic energy
        PS::F64vec vel = system[i].getVel();
        eng_kin_loc += system[i].getMass()*(vel*vel);
        
        //--- potential energy
        PS::F64vec virt = system[i].getVirtCoulomb();
        pot_LJ_loc      += system[i].getPotLJ();
        pot_coulomb_loc += (virt[0] + virt[1] + virt[2]);  // using coulomb virial
        pot_intra_loc   += system[i].getPotIntra();

        //--- error check
        if (std::isnan(vel.x) || std::isnan(vel.y) || std::isnan(vel.z) ||
            std::isnan(virt.x) || std::isnan(virt.y) || std::isnan(virt.z) ||
            std::isnan(pot_LJ_loc) || std::isnan(pot_coulomb_loc) ||
            std::isnan(pot_intra_loc)) {
            std::cout << "i           = " << i << std::endl;
            std::cout << "vel.x       = " << vel.x << std::endl;
            std::cout << "vel.y       = " << vel.y << std::endl;
            std::cout << "vel.z       = " << vel.z << std::endl;
            std::cout << "virt.x      = " << virt.x << std::endl;
            std::cout << "virt.y      = " << virt.y << std::endl;
            std::cout << "virt.z      = " << virt.z << std::endl;
            std::cout << "pot (LJ)    = " << pot_LJ_loc << std::endl;
            std::cout << "pot (clmb)  = " << pot_coulomb_loc << std::endl;
            std::cout << "pot (intra) = " << pot_intra_loc << std::endl;
            MPI::COMM_WORLD.Abort(9);
            std::exit(1);
         }
    }
    
    eng_kin_loc = 0.5*eng_kin_loc;
    
    //--- calculate total sumation
    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
      energy_buf->setEngKin(     PS::Comm::getSum(eng_kin_loc) );
      energy_buf->setPotLJ(      PS::Comm::getSum(pot_LJ_loc) );
      energy_buf->setPotCoulomb( PS::Comm::getSum(pot_coulomb_loc) );
      energy_buf->setPotIntra(   PS::Comm::getSum(pot_intra_loc) );
    #else
      energy_buf->setEngKin(     eng_kin_loc );
      energy_buf->setPotLJ(      pot_LJ_loc );
      energy_buf->setPotCoulomb( pot_coulomb_loc );
      energy_buf->setPotIntra(   pot_intra_loc );
    #endif
    
    energy_buf->calcTotal();
}
