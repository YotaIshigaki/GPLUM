#pragma once
#include "../ps_defs.hpp"
#include "particle.hpp"
#include "fmm.hpp"

namespace ParticleSimulator {

    namespace ParticleMeshMultipole {

        template <typename real_t=double, typename cplx_t = std::complex<real_t> >
        class Cell_FMM {
        public:
            F64vec center;
            std::vector<Particle *> plist;
            typedef std::vector<Particle *>::iterator piter;
            typedef std::vector<Particle *>::const_iterator cpiter;
            
            MultipoleMoment<real_t, cplx_t> mm;
            LocalExpansion<real_t, cplx_t>  le;
    
            Cell_FMM() {}
    
            void init(const S32 p) {
                mm.alloc(p);
                le.alloc(p);
            }
    
            void setPos(const F64vec & pos) {
                center = pos;
            }
    
            inline void setParticleToCell(Particle * ptr_to_ptcl) {
                this->plist.push_back(ptr_to_ptcl);
            }
        
            void do_P2M(){
                const S32 nk = plist.size();
                for (S32 k=0; k<nk; k++){
                    mm.assign_particle(center, plist[k]->getPos(), plist[k]->getCharge());
                }
            }
            void do_L2P(){
                const S32 nk = plist.size();
                for (S32 k=0; k < nk; k++){
                    LocalExpansion<real_t, cplx_t> le1;
                    le1.alloc(1);
                    le1.assign_from_LE(le, plist[k]->getPos(), center);
                    F64 phi, ax, ay, az;
                    le1.read_phi_and_grad(phi, ax, ay, az);
                    plist[k]->clear();
                    plist[k]->acc += F64vec(ax, ay, az);
                    plist[k]->pot += phi;
                }
            }
            void do_L2P_corr(const F64 msum, const F64 alpha){
                const F64 pi = 4.0 * atan(1.0);
                const F64 ainv2 = pi / (alpha * alpha);
                const F64 mc = ((2./3.) * pi) * msum;
                const S32 nk = plist.size();
                for (S32 k=0; k < nk; k++){
                    const F64vec dr = plist[k]->getPos() - center;
                    plist[k]->acc += (2.0 * mc) * dr;
                    plist[k]->pot += mc * (dr*dr);
                    plist[k]->pot -= (msum*ainv2);
                }
            }
            
            F64 dispersion() const{
                const S32 nk = plist.size();
                F64 ret = 0.0;
                for (S32 k=0; k < nk; k++){
                    const F64vec dr = plist[k]->getPos() - center;
                    ret += plist[k]->getCharge() * (dr*dr);
                }
                return ret;
            }
        };

    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
