#pragma once
#include "../ps_defs.hpp"

namespace ParticleSimulator {

    namespace ParticleMeshMultipole {

        template <class real_t, class cplx_t, class Tepi, class Tforce>
        class Cell_FMM {
        public:
            S32 idx;
            F64vec center;

            S32 n_epi; 
            Tepi * epi_first; 
            Tforce * force_first; 
            // The information of i particles is managed by pairs of the
            // begining address and the number of particles.
            // epi_first and force_first are the pointers to epi_sorted_[]
            // and force_sorted_[] in TreeForForce class, respectively.
 
            bool is_mm_defined;           
            MultipoleMoment<real_t, cplx_t> mm;
            LocalExpansion<real_t, cplx_t>  le;
    
            Cell_FMM() {
                center = F64vec(0.0);

                n_epi = 0;
                epi_first = nullptr;
                force_first = nullptr;
    
                is_mm_defined = {false};
            }
    
            void init(const S32 p) {
                mm.alloc(p);
                le.alloc(p);
            }

            void setIdx(const S32 idx) {
                this->idx = idx;
            }

            void setPos(const F64vec & pos) {
                center = pos;
            }

            inline void setIParticleInfoToCell(const S32 n_epi,
                                               Tepi * epi_first,
                                               Tforce * force_first) {
                this->n_epi = n_epi;
                this->epi_first = epi_first;
                this->force_first = force_first;
            }

            void setMultipoleMoment(MultipoleMoment<real_t, cplx_t> & _mm) {
                assert(_mm.buf.size() == mm.buf.size());
                is_mm_defined = true;
                for (S32 i=0; i<mm.buf.size(); i++) mm.buf[i] = _mm.buf[i];
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
        };

    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
