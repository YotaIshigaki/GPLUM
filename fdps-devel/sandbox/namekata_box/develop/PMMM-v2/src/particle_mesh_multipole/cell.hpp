#pragma once
#include "../ps_defs.hpp"
#include "fmm.hpp"

namespace ParticleSimulator {

    namespace ParticleMeshMultipole {

        template <class real_t, class cplx_t,
                  class Tforce, class Tepi, class Tepj, class Tspj>
        class Cell_FMM {
        public:
            F64vec center;

            bool epj_is_shared;
            std::vector<S32> n_epj;
            std::vector<Tepj * > epj_first; 
            std::vector<S32> n_spj;
            std::vector<Tspj * > spj_first;
            // The information of j particles used to calculate multipole moments
            // of this PM cell is managed by lists of pairs of the begining
            // address of particle array and the number of particles.
            // epj_first[][] and spj_first[][] point to epj_sorted_[] and
            // spj_sorted_[] in TreeForForce class, respectively.

            S32 n_epi; 
            Tepi * epi_first; 
            Tforce * force_first; 
            // The information of i particles is managed by pairs of the
            // begining address and the number of particles.
            // epi_first and force_first are the pointers to epi_sorted_[]
            // and force_sorted_[] in TreeForForce class, respectively.
            
            MultipoleMoment<real_t, cplx_t> mm;
            LocalExpansion<real_t, cplx_t>  le;
    
            Cell_FMM() {
                center = F64vec(0.0);
                
                epj_is_shared = false;
                n_epj.clear();
                epj_first.clear();
                n_spj.clear();
                spj_first.clear();

                n_epi = 0;
                epi_first = nullptr;
                force_first = nullptr;
            }
    
            void init(const S32 p) {
                mm.alloc(p);
                le.alloc(p);
            }
    
            void setPos(const F64vec & pos) {
                center = pos;
            }
    
            inline void setEpiToCell(const S32 n_epi, Tepi * epi_first) {
                this->n_epi = n_epi;
                this->epi_first = epi_first;
            }
            inline void setEpjToCell(const bool epj_is_shared,
                                     const std::vector<S32> & n_epj,
                                     const std::vector<Tepj * > & epj_first) {
                this->epj_is_shared = epj_is_shared;
                this->n_epj = n_epj;
                this->epj_first = epj_first;
            }
            inline void setSpjToCell(const std::vector<S32> & n_spj,
                                     const std::vector<Tspj * > & spj_first) {
                this->n_spj = n_spj;
                this->spj_first = spj_first;
            }
            inline void setForceToCell(Tforce * force_first) {
                this->force_first = force_first;
            }

            F64 getTotalCharge(const F64ort & pos_my_domain) const {
                F64 msum = 0;
                if (!epj_is_shared) {
                    for (S32 n=0; n < n_epj.size(); n++) {
                        for (S32 k=0; k < n_epj[n]; k++) {
                            msum += epj_first[n][k].getCharge();
                        }
                    }
                    for (S32 n=0; n < n_spj.size(); n++) {
                        for (S32 k=0; k < n_spj[n]; k++) {
                            msum += spj_first[n][k].getCharge();
                        }
                    }
                } else {
                    for (S32 n=0; n < n_epj.size(); n++) {
                        for (S32 k=0; k < n_epj[n]; k++) {
                            const F64vec pos = epj_first[n][k].getPos();
                            if (pos_my_domain.contains(pos)) {
                                msum += epj_first[n][k].getCharge();
                            }
                        }
                    }
                }
                return msum;
            }
        
            void do_P2M(const F64ort & pos_my_domain){
                if (!epj_is_shared) {
                    for (S32 n=0; n < n_epj.size(); n++) {
                        for (S32 k=0; k < n_epj[n]; k++){
                            mm.assign_particle(center,
                                               epj_first[n][k].getPos(),
                                               epj_first[n][k].getCharge());
                        }
                    }
                    for (S32 n=0; n < n_spj.size(); n++) {
                        for (S32 k=0; k < n_spj[n]; k++){
                            mm.assign_particle(center,
                                               spj_first[n][k].getPos(),
                                               spj_first[n][k].getCharge());
                        }
                    }
                } else {
                    for (S32 n=0; n < n_epj.size(); n++) {
                        for (S32 k=0; k < n_epj[n]; k++){
                            const F64vec pos = epj_first[n][k].getPos();
                            if (pos_my_domain.contains(pos)) {
                                mm.assign_particle(center, pos,
                                                   epj_first[n][k].getCharge());
                            }
                        }
                    }
                }
            }
            void do_L2P(const bool clear = true){
                for (S32 k=0; k < n_epi; k++){
                    LocalExpansion<real_t, cplx_t> le1;
                    le1.alloc(1);
                    le1.assign_from_LE(le, epi_first[k].getPos(), center);
                    F64 pot, ax, ay, az;
                    le1.read_phi_and_grad(pot, ax, ay, az);
                    if (clear) force_first[k].clearPMM();
                    force_first[k].accumulateForcePMM(F64vec(ax, ay, az), pot);
                }
            }
            void do_L2P_corr(const F64 msum, const F64 alpha){
                const F64 pi = 4.0 * atan(1.0);
                const F64 ainv2 = pi / (alpha * alpha);
                const F64 mc = ((2./3.) * pi) * msum;
                for (S32 k=0; k < n_epi; k++){
                    const F64vec dr = epi_first[k].getPos() - center;
                    const F64vec acc_corr = (2.0 * mc) * dr;
                    const F64 pot_corr = mc * (dr*dr) - (msum*ainv2);
                    force_first[k].accumulateForcePMM(acc_corr, pot_corr);
                }
            }
            
            F64 dispersion(const F64ort & pos_my_domain) const{
                F64 ret = 0.0;
                if (!epj_is_shared) {
                    for (S32 n=0; n < n_epj.size(); n++) {
                        for (S32 k=0; k < n_epj[n]; k++){
                            const F64vec dr = epj_first[n][k].getPos() - center;
                            ret += epj_first[n][k].getCharge() * (dr*dr);
                        }
                    }
                    for (S32 n=0; n < n_spj.size(); n++) {
                        for (S32 k=0; k < n_spj[n]; k++){
                            const F64vec dr = spj_first[n][k].getPos() - center;
                            ret += spj_first[n][k].getCharge() * (dr*dr);
                        }
                    }
                } else {
                    for (S32 n=0; n < n_epj.size(); n++) {
                        for (S32 k=0; k < n_epj[n]; k++){
                            const F64vec pos = epj_first[n][k].getPos();
                            if (pos_my_domain.contains(pos)) {
                                const F64vec dr = pos - center;
                                ret += epj_first[n][k].getCharge() * (dr*dr);
                            }
                        }
                    }
                }
                return ret;
            }
        };

    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
