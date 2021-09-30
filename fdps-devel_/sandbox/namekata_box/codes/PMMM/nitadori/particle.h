#pragma once
#include <cstdlib>
#include "vector3.h"

struct Particle{
	dvec3  pos;
	double mass;
	dvec3  acc_direct;
	double phi_direct;
	dvec3  acc_app;
	double phi_app;

	void clear_pot(){
		acc_direct = dvec3(0.0);
		phi_direct = 0.0;
		acc_app    = dvec3(0.0);
		phi_app    = 0.0;
	}
	void move_accp(){
		acc_app += acc_direct;
		phi_app += phi_direct;
		acc_direct = 0.0;
		phi_direct = 0.0;
	}
	dvec3 avecdiff_rel() const {
		return (acc_app - acc_direct) / acc_direct.abs();
	}
	double adiff_rel() const {
		return (acc_app - acc_direct).abs() / acc_direct.abs();
	}
	double adiff_abs() const {
		return (acc_app - acc_direct).abs();
	}
	double aabs() const {
		return acc_direct.abs();
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

	static dvec3 rand_vec(){
		return dvec3(drand48(), drand48(), drand48());
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
