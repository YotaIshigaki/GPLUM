// To make a file output of Ewald sum
#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>

#include "vector3.h"
inline dvec3 minimum_image(const dvec3 &inp){
	return dvec3(
			inp.x - round(inp.x),
			inp.y - round(inp.y),
			inp.z - round(inp.z));
}
// #define NON_CHARGE_NEUTRAL
#include "ewald.h"

int main(){
	enum{
		NP   = 128,
	};

	static Particle ptcl[NP];
	Particle::gen_rand_dist(NP, ptcl);
	double msum = 0.0;
	for(int i=0; i<NP; i++){
		msum += ptcl[i].mass;
	}

	puts("Ewald sum");

	const double alpha= 2.4;
	eval_k_space<5>(NP, alpha, ptcl);
	eval_r_space<3>(NP, alpha, msum, ptcl);

#if 0
	for(int i=0; i<NP; i++){
		const dvec3 acc = ptcl[i].acc_direct;
		printf("%24.16e, (%24.16e, %24.16e, %24.16e)\n", 
				ptcl[i].phi_direct, acc.x, acc.y, acc.z);
	}
#endif
	FILE *fp;
	assert(fp = fopen("qpos.dat", "w"));
	fprintf(fp, "%d\n", NP);
	for(int i=0; i<NP; i++){
		fprintf(fp, "%A  %A %A %A\n",
				ptcl[i].mass,
				ptcl[i].pos.x,
				ptcl[i].pos.y,
				ptcl[i].pos.z);
	}
	fclose(fp);

	assert(fp = fopen("ewald.dat", "w"));
	fprintf(fp, "%d\n", NP);
	for(int i=0; i<NP; i++){
		fprintf(fp, "%A  %A %A %A\n",
				ptcl[i].phi_direct,
				ptcl[i].acc_direct.x,
				ptcl[i].acc_direct.y,
				ptcl[i].acc_direct.z);
	}
	fclose(fp);

	return 0;
};
