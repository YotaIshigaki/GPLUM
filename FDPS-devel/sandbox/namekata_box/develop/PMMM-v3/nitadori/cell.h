#pragma once
#include <vector>
#include "vector3.h"
#include "fmm.h"
#include "particle.h"

typedef vector3<int> ivec3;

struct Cell{
    dvec3  center;
    double length;
    std::vector<Particle *> plist;
    typedef std::vector<Particle *>::iterator piter;
    typedef std::vector<Particle *>::const_iterator cpiter;

    void set(const ivec3 &idx, const double d){
        center = cell_pos(idx, d);
        length = d;
        plist.clear();
    }

    void sanity_check() const {
        for(cpiter it = plist.begin(); it != plist.end(); ++it){
            const dvec3 dr = (*it)->pos - center;
            assert(fabs(dr.x) <= 0.5*length);
            assert(fabs(dr.y) <= 0.5*length);
            assert(fabs(dr.z) <= 0.5*length);
        }
    }

    void dump_plist(std::ostream & fout=std::cout) const {
        for(cpiter it = plist.begin(); it != plist.end(); ++it){
            fout << (*it)->pos.x << " "
                 << (*it)->pos.y << " "
                 << (*it)->pos.z << " "
                 << (*it)->mass << std::endl;
        }
    }

};

template <int p>
struct Cell_FMM : public Cell {
    MultipoleMoment<p> mm;
    LocalExpansion <p> le;

    void do_P2M(){
        const int nk = plist.size();
        for(int k=0; k<nk; k++){
            // const dvec3 dr = center - plist[k]->pos;
            mm.assign_particle(center, plist[k]->pos, plist[k]->mass);
        }
    }
    void do_L2P(){
        const int nk = plist.size();
        for(int k=0; k<nk; k++){
            // const dvec3 dr = plist[k]->pos - center;
            LocalExpansion<1> le1;
            // ---- debug ----
            //const double x = plist[k]->pos.x;
            //const double y = plist[k]->pos.y;
            //const double z = plist[k]->pos.z;
            bool debug_flag {false};
            //if (0.0959307 < x && x < 0.0959309) debug_flag = true;
            // ---- debug ----
            le1.assign_from_LE(le, plist[k]->pos, center, debug_flag);
            double phi, ax, ay, az;
            le1.read_phi_and_grad(phi, ax, ay, az);
            // ---- debug ----
            //if (0.0959307 < x && x < 0.0959309) {
            //    std::cout << "center = " << center.x << "   "
            //              << center.y << "   " << center.z << std::endl;
            //    std::cout << "--------------------" << std::endl;
            //    std::cout << "Content of le:" << std::endl;
            //    for (int lm = 0; lm < le.size(); lm++) {
            //        std::cout << "le.buf[" << lm << "] = " << le.buf[lm] << std::endl;
            //    }
            //    std::cout << "--------------------" << std::endl;
            //    std::cout << "Content of le1:" << std::endl;
            //    for (int lm = 0; lm < le1.size(); lm++) {
            //        std::cout << "le1.buf[" << lm << "] = " << le1.buf[lm] << std::endl;
            //    }
            //    std::cout << "--------------------" << std::endl;
            //    std::cout << "pos = " << x << "   " << y << "   " << z  << std::endl;
            //    std::cout << "acc = " << ax << "   " << ay << "   " << az << std::endl;
            //    std::cout << "--------------------" << std::endl;
            //}
            // ---- debug ----
            plist[k]->acc_app += dvec3(ax, ay, az);
            plist[k]->phi_app += phi;
        }
    }
    void do_L2P_corr(const double msum, const double alpha){
        const double pi = 4.0 * atan(1.0);
        const double ainv2 = pi / (alpha * alpha);
        const double mc = ((2./3.) * pi) * msum;
        const int nk = plist.size();
        for(int k=0; k<nk; k++){
            const dvec3 dr = plist[k]->pos - center;
            plist[k]->acc_app += (2.0 * mc) * dr;
            plist[k]->phi_app += mc * (dr*dr);
            plist[k]->phi_app -= (msum*ainv2);
        }
    }

    double msum() const{
        const int nk = plist.size();
        double ret = 0.0;
        for(int k=0; k<nk; k++){
            ret += plist[k]->mass;
        }
        return ret;
    }
    
    double dispersion() const{
        const int nk = plist.size();
        double ret = 0.0;
        for(int k=0; k<nk; k++){
            const dvec3 dr = plist[k]->pos - center;
            ret += plist[k]->mass * (dr*dr);
        }
        return ret;
    }
};
