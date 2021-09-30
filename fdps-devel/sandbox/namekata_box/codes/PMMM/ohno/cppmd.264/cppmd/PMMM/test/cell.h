#pragma once
#include <algorithm>
#include <vector>
#include "vector3.h"
#include "fmm.h"
#include "particle.h"

typedef vector3<int> ivec3;

inline ivec3 cell_nearest(const dvec3 &pos, const double d){
	return ivec3(pos / d);
}

inline dvec3 minimum_image(const dvec3 &inp){
	return dvec3(
			inp.x - round(inp.x),
			inp.y - round(inp.y),
			inp.z - round(inp.z));
}

inline dvec3 cell_pos(const ivec3 &idx, const double d){
        return d * (dvec3(0.5) + dvec3(idx));
}

struct CellRange{
  int min[3];
  int max[3];
  int size;

  void recalc_size(){
    size = (max[0]-min[0])*(max[1]-min[1])*(max[2]-min[2]);
  }

  void reorder(){
    if(min[0]>max[0])std::swap(min[0],max[0]);
    if(min[1]>max[1])std::swap(min[1],max[1]);
    if(min[2]>max[2])std::swap(min[2],max[2]);
  }

  void set_min(const int m[3]){
    min[0] = m[0]; min[1] = m[1]; min[2] = m[2];
    reorder();
    recalc_size();
  }

  void set_max(const int m[3]){
    max[0] = m[0]; max[1] = m[1]; max[2] = m[2];
    reorder();
    recalc_size();
  }

  void set_min_max(const int _min[3], const int _max[3]){
    min[0] = _min[0]; min[1] = _min[1]; min[2] = _min[2];
    max[0] = _max[0]; max[1] = _max[1]; max[2] = _max[2];
    reorder();
    recalc_size();
  }

};

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
			le1.assign_from_LE(le, plist[k]->pos, center);
			double phi, ax, ay, az;
			le1.read_phi_and_grad(phi, ax, ay, az);
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
