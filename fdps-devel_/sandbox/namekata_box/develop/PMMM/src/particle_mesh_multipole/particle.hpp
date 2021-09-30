#pragma once
#include "../ps_defs.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        class Particle {
        public:
            F64vec pos;
            F64 charge;
            F64vec acc;
            F64 pot;

            F64vec getPos() const { return this->pos; }
            void setPos(const F64vec & pos) { this->pos = pos; }
            F64 getCharge() const { return this->charge; }
            F64vec getAcc() const { return this->acc; }
            F64 getPot() const {return this->pot; }
            void clear() {
                this->acc = 0.0;
                this->pot = 0.0;
            }
        };

    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
