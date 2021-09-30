#pragma once
#include <cstdlib>
#include "../ps_defs.hpp"

namespace ParticleSimulator {

    struct Particle{
        F64vec  pos;
        F64     mass;
        F64vec  acc_direct;
        F64     phi_direct;
        F64vec  acc_app;
        F64     phi_app;
    
        void clear_pot(){
            acc_direct = F64vec(0.0);
            phi_direct = 0.0;
            acc_app    = F64vec(0.0);
            phi_app    = 0.0;
        }
        void move_accp(){
            acc_app += acc_direct;
            phi_app += phi_direct;
            acc_direct = 0.0;
            phi_direct = 0.0;
        }
        F64vec avecdiff_rel() const {
            const F64 acc_direct_abs = std::sqrt(acc_direct * acc_direct);
            return (acc_app - acc_direct) / acc_direct_abs;
        }
        double adiff_rel() const {
            const F64vec acc = acc_app - acc_direct;
            const F64 acc_abs = std::sqrt(acc * acc);
            const F64 acc_direct_abs = std::sqrt(acc_direct * acc_direct);
            return acc_abs / acc_direct_abs;
        }
        double adiff_abs() const {
            const F64vec acc = acc_app - acc_direct;
            return std::sqrt(acc * acc);
        }
        double aabs() const {
            return std::sqrt(acc_direct * acc_direct);
        }
        double pdiff_rel() const {
            return fabs((phi_app - phi_direct) / phi_direct);
        }
        double pdiff_abs() const {
            return fabs(phi_app - phi_direct);
        }
        double pabs() const {
            return fabs(phi_direct);
        }
    
        static F64vec rand_vec(){
            return F64vec(drand48(), drand48(), drand48());
        }
    
        static void gen_rand_dist(
                const int NP,
                Particle ptcl[],
                const long seed = 19810614)
        {
            srand48(seed);
    
            for(int i=0; i<NP; i++){
                ptcl[i].pos = rand_vec();
                ptcl[i].clear_pot();
            }
            double msum = 0.0;
            for(int i=0; i<NP; i++){
                msum += 
                    (ptcl[i].mass  = drand48() * (1.0/NP));
            }
#ifndef NON_CHARGE_NEUTRAL
            for(int i=0; i<NP; i++){
                ptcl[i].mass -= msum / NP;
            }
#endif
        }
    };

} // END of namespace of ParticleSimulator
