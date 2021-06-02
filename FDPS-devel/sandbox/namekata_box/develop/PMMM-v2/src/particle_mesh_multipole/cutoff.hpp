#pragma once
#include <cmath>
#include <complex>
#include <vector>
#include "../ps_defs.hpp"

namespace ParticleSimulator {

    namespace ParticleMeshMultipole {

        template <typename real_t=double, typename cplx_t=std::complex<real_t>>
        struct CutFunc{
            S32 LMAX;
            std::vector<real_t> gammainv;
            std::vector<real_t> rcut;
            std::vector<cplx_t> kcut;
    
            CutFunc() {}
    
            void init(S32 _LMAX){
                LMAX = _LMAX;
                gammainv.reserve(1+LMAX);
                rcut.reserve(1+LMAX);
                kcut.reserve(1+LMAX);
                F64 val = sqrt(4.0*atan(1.0));
                for(S32 ell=0; ell<=LMAX; ell++){
                    gammainv[ell] = 1.0 / val;
                    val *= (0.5 + ell);
                }
            }
    
            void eval_rcut(const F64 x){
                const F64 xexp = exp(-x);
                F64 xpow = sqrt(x);
                F64 val  = sqrt(4.0*atan(1.0)) * erfc(xpow);
                for(S32 ell=0; ell<=LMAX; ell++){
                    rcut[ell] = val * gammainv[ell];
                    val = (0.5 + ell)*val + xpow * xexp;
                    xpow *= x;
                }
            }
    
            void eval_kcut(const F64 kn2, const F64 kappa){
                const F64 pi = 4.0 * atan(1.0);
                const F64 gauss = exp(-kn2 * ((pi/kappa)*(pi/kappa)));
    
                const cplx_t ratio(0.0, -pi * kn2);
                cplx_t coef = 1.0 / sqrt(pi * kn2);
                for(S32 ell=0; ell<=LMAX; ell++){
                    kcut[ell] = coef * (gauss * gammainv[ell]);
                    coef *= ratio;
                }
            }

            void freeMem() {
                // This function explicitly free memories.
                std::vector<real_t>().swap(gammainv);
                std::vector<real_t>().swap(rcut);
                std::vector<cplx_t>().swap(kcut);
            }
        };

    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
