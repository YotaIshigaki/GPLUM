#pragma once
#include <vector>
#include "../ps_defs.hpp"
#include "fmm.hpp"
#include "particle.hpp"

namespace ParticleSimulator {

    inline S32vec cell_nearest(const F64vec &pos, const F64 d){
    	return S32vec(pos / d);
    }
    inline F64vec cell_pos(const S32vec &idx, const F64 d){
    	return d * (F64vec(0.5) + F64vec(idx));
    }

    inline F64vec minimum_image(const F64vec &inp){
        // This function is used in both pbc.cpp and ewald.hpp.
        // This should be replaced by a function that supports
        // arbitrary box sizes (bx, by, bz).
    	return F64vec(
    			inp.x - round(inp.x),
    			inp.y - round(inp.y),
    			inp.z - round(inp.z));
    }

    struct Cell{
        F64vec  center;
        F64     length;
        std::vector<Particle *> plist;
        typedef std::vector<Particle *>::iterator piter;
        typedef std::vector<Particle *>::const_iterator cpiter;
    
        void set(const S32vec &idx, const F64 d){
            center = cell_pos(idx, d);
            length = d;
            plist.clear();
        }
    
        void sanity_check() const {
            for(cpiter it = plist.begin(); it != plist.end(); ++it){
                const F64vec dr = (*it)->pos - center;
                assert(fabs(dr.x) <= 0.5*length);
                assert(fabs(dr.y) <= 0.5*length);
                assert(fabs(dr.z) <= 0.5*length);
            }
        }
    };
    
    template <int p>
    struct Cell_FMM : public Cell {
        MultipoleMoment<p> mm;
        LocalExpansion <p> le;
    
        void do_P2M(){
            const S32 nk = plist.size();
            for(S32 k=0; k<nk; k++){
                // const F64vec dr = center - plist[k]->pos;
                mm.assign_particle(center, plist[k]->pos, plist[k]->mass);
            }
        }
        void do_L2P(){
            const S32 nk = plist.size();
            for(S32 k=0; k<nk; k++){
                // const F64vec dr = plist[k]->pos - center;
                LocalExpansion<1> le1;
                le1.assign_from_LE(le, plist[k]->pos, center);
                F64 phi, ax, ay, az;
                le1.read_phi_and_grad(phi, ax, ay, az);
                plist[k]->acc_app += F64vec(ax, ay, az);
                plist[k]->phi_app += phi;
            }
        }
        void do_L2P_corr(const F64 msum, const F64 alpha){
            const F64 pi = 4.0 * atan(1.0);
            const F64 ainv2 = pi / (alpha * alpha);
            const F64 mc = ((2./3.) * pi) * msum;
            const S32 nk = plist.size();
            for(S32 k=0; k<nk; k++){
                const F64vec dr = plist[k]->pos - center;
                plist[k]->acc_app += (2.0 * mc) * dr;
                plist[k]->phi_app += mc * (dr*dr);
                plist[k]->phi_app -= (msum*ainv2);
            }
        }
        
        double dispersion() const{
            const S32 nk = plist.size();
            F64 ret = 0.0;
            for(S32 k=0; k<nk; k++){
                const F64vec dr = plist[k]->pos - center;
                ret += plist[k]->mass * (dr*dr);
            }
            return ret;
        }
    };

} // END of namespace of ParticleSimulator
