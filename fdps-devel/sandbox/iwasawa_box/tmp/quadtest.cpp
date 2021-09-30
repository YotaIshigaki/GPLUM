#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "phantomquad.hpp"

struct Particle{
	double x, y, z, m;
	double q_xx, q_yy, q_zz;
	double q_xy, q_yz, q_zx;

	void set_rand(double mscale, double qscale){
		x = drand48() - 0.5;
		y = drand48() - 0.5;
		z = drand48() - 0.5;
		m = mscale * drand48();

		double dx = qscale * (drand48() - 0.5);
		double dy = qscale * (drand48() - 0.5);
		double dz = qscale * (drand48() - 0.5);

		q_xx = dx*dx;
		q_yy = dy*dy;
		q_zz = dz*dz;
		q_xy = dx*dy;
		q_yz = dy*dz;
		q_zx = dz*dx;
	}

	static void grav_mono(
			const Particle &pi, 
			const Particle &pj, 
			const double eps2,
			double accp[4])
	{
		double dx = pj.x - pi.x;
		double dy = pj.y - pi.y;
		double dz = pj.z - pi.z;
		double r2 = eps2 + dx*dx + dy*dy + dz*dz;
		double ri1 = 1.0 / sqrt(r2);
		double ri2 = ri1 * ri1;
		double mri1 = pj.m * ri1;
		double mri3 = mri1 * ri2;
		accp[0] += mri3 * dx;
		accp[1] += mri3 * dy;
		accp[2] += mri3 * dz;
		accp[3] -= mri1;
	}
	static void grav_quad(
			const Particle &pi, 
			const Particle &pj, 
			const double eps2,
			double accp[4])
	{
		double dx = pj.x - pi.x;
		double dy = pj.y - pi.y;
		double dz = pj.z - pi.z;
		double r2 = eps2 + dx*dx + dy*dy + dz*dz;

		double ri1 = 1.0 / sqrt(r2);
		double ri2 = ri1 * ri1;
		double ri3 = ri1 * ri2;
		double ri4 = ri2 * ri2;
		double ri5 = ri2 * ri3;

		double tr = pj.q_xx + pj.q_yy + pj.q_zz;
		double q_xx = 3.0 * pj.q_xx - tr;
		double q_yy = 3.0 * pj.q_yy - tr;
		double q_zz = 3.0 * pj.q_zz - tr;
		double q_xy = 3.0 * pj.q_xy;
		double q_yz = 3.0 * pj.q_yz;
		double q_zx = 3.0 * pj.q_zx;

		double qr_x = q_xx*dx + q_xy*dy + q_zx*dz;
		double qr_y = q_xy*dx + q_yy*dy + q_yz*dz;
		double qr_z = q_zx*dx + q_yz*dy + q_zz*dz;

		double rqr = -eps2*tr + qr_x*dx + qr_y*dy + qr_z*dz;

		double meff05     = (pj.m + 0.5*(rqr*ri4));
		double meff25_ri3 = (pj.m + 2.5*(rqr*ri4)) *ri3;

		accp[0] += meff25_ri3*dx - ri5*qr_x;
		accp[1] += meff25_ri3*dy - ri5*qr_y;
		accp[2] += meff25_ri3*dz - ri5*qr_z;
		// if(r2 > eps2)
		accp[3] -= meff05 * ri1;
	}
};

void dump_accp(
		const int ni,
		const double accp[][4],
		const char *fname)
{
	FILE *fp = fopen(fname, "w");
	assert(fp);
	for(int i=0; i<ni; i++){
		fprintf(fp, "%4d : %+20.10e, %+20.10e, %+20.10e : %+20.10e\n",
				i, accp[i][0], accp[i][1], accp[i][2], accp[i][3]);
	}
	fclose(fp);
}

int main(){
	enum{
		NI = 10,
		NJ = 100,
	};

	const double eps2 = 1.0e-2;

	static Particle ptcl[NJ];
	for(int i=0; i<NJ; i++) ptcl[i].set_rand(1./NJ, 1.0/NJ);

	static double accp_mono_ref[NI][4];
	static double accp_mono_pgk[NI][4];
	static double accp_quad_ref[NI][4];
	static double accp_quad_pgk[NI][4];

	for(int i=0; i<NI; i++){
		for(int k=0; k<4; k++) {
			accp_mono_ref[i][k] = 0.0;
			accp_mono_pgk[i][k] = 0.0;
			accp_quad_ref[i][k] = 0.0;
			accp_quad_pgk[i][k] = 0.0;
		}
	}

	for(int i=0; i<NI; i++){
		for(int j=0; j<NJ; j++) {
			Particle::grav_mono(ptcl[i], ptcl[j], eps2, accp_mono_ref[i]);
		}
	}
	dump_accp(NI, accp_mono_ref, "accp_mono_ref");

	static PhantomGrapeQuad pg;
	pg.set_eps2(eps2);
	for(int i=0; i<NI; i++){
		const Particle &p = ptcl[i];
		pg.set_xi_one(i, p.x, p.y, p.z);
	}
	for(int j=0; j<NJ; j++) {
		const Particle &p = ptcl[j];
		pg.set_epj_one(j, p.x, p.y, p.z, p.m);
	}
	pg.run_epj(NI, NJ);
	for(int i=0; i<NI; i++){
		double *a = accp_mono_pgk[i];
		pg.get_accp_one(i, a[0], a[1], a[2], a[3]);
	}
	dump_accp(NI, accp_mono_pgk, "accp_mono_pgk");

	for(int i=0; i<NI; i++){
		for(int j=0; j<NJ; j++) {
			Particle::grav_quad(ptcl[i], ptcl[j], eps2, accp_quad_ref[i]);
		}
	}
	dump_accp(NI, accp_quad_ref, "accp_quad_ref");

	for(int j=0; j<NJ; j++) {
		const Particle &p = ptcl[j];
		pg.set_spj_one(j, p.x, p.y, p.z, p.m, 
				p.q_xx, p.q_yy, p.q_zz, p.q_xy, p.q_yz, p.q_zx);
	}
	pg.run_spj(NI, NJ);
	for(int i=0; i<NI; i++){
		double *a = accp_quad_pgk[i];
		pg.get_accp_one(i, a[0], a[1], a[2], a[3]);
	}
	dump_accp(NI, accp_quad_pgk, "accp_quad_pgk");

	return 0;
}
