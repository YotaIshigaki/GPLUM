#pragma once

#include <cassert>
#include "v2r8.h"

#if 1
#include "mathfunc.h"
#include "kernel.h"
#include "EoS.h"
#include "class.h"
#include "param.h"
#else
#include "header.h"
#endif

namespace kernel{
	static const double pi = 4.0 * atan(1.0);
	static const double supportRadius = 3.5;
	static const double normalzeFactor = 1365./(64.*pi);
	// static const double normalzeFactor = 6.788953041263661;

	template <typename real>
	inline real pow4(const real x){
		real x2 = x*x;
		return x2*x2;
	}
	template <typename real>
	inline real pow7(const real x){
		real x2 = x*x;
		real x4 = x2*x2;
		return x*x2*x4;
	}
	template <typename real>
	inline real pow8(const real x){
		real x2 = x*x;
		real x4 = x2*x2;
		return x4*x4;
	}

	__attribute__ ((always_inline))
	inline v2r8 W(const v2r8 s){
		v2r8 plus =  v2r8(1.0) - s;
		plus = v2r8::max(plus, v2r8(0.0));
		v2r8 r_value = pow8(plus);
		r_value *= (v2r8(1.0) + s*(v2r8(8.0) + s*(v2r8(25.0) + s*(v2r8(32.0)))));
		return r_value;
	}
	__attribute__ ((always_inline))
	inline v2r8 gradW(const v2r8 s, const v2r8 rinv){
		v2r8 plus =  v2r8(1.0) - s;
		plus = v2r8::max(plus, v2r8(0.0));
		v2r8 r_value = pow7(plus);
		r_value *= (plus * (v2r8(8.0) + s*(v2r8(50.0) + s*(v2r8(96.0)))))
		           - (v2r8(8.0) + s*(v2r8(64.0) + s*(v2r8(200.0) + s*(v2r8(256.0)))));
		return (r_value * rinv);
	}
}

struct CalcDensity{
	__attribute__ ((noinline))
	void operator () (
			const EPI::Dens * const __restrict ep_i, 
			const PS::S32                      Nip, 
			const EPJ::Dens * const __restrict ep_j, 
			const PS::S32                      Njp, 
			RESULT::Dens * const __restrict    dens)
	{
		for(int i=0; i<Nip; i+=2){
			const v2r8 xi(ep_i[i+0].pos.x, ep_i[i+1].pos.x);
			const v2r8 yi(ep_i[i+0].pos.y, ep_i[i+1].pos.y);
			const v2r8 zi(ep_i[i+0].pos.z, ep_i[i+1].pos.z);
			const v2r8 hi(ep_i[i+0].smth,  ep_i[i+1].smth);
			const v2r8 Hinv = (hi * v2r8(kernel::supportRadius) ).rcpa_x8();
			const v2r8 eps2(1.e-16);

			v2r8 rho(0.0);

			for(int j=0; j<Njp; j++){
				const double *ptr = (const double *)(ep_j + j);
				const v2r8 jp0 = v2r8::load(ptr + 0);
				const v2r8 jp1 = v2r8::load(ptr + 2);
				const v2r8_bcl mj(jp0);
				const v2r8_bch xj(jp0);
				const v2r8_bcl yj(jp1);
				const v2r8_bch zj(jp1);

				const v2r8 dx = xj - xi;
				const v2r8 dy = yj - yi;
				const v2r8 dz = zj - zi;

				const v2r8 r2   = ((eps2 + dx*dx) + dy*dy) + dz*dz;
				const v2r8 rinv = r2.rsqrta_x8();
				const v2r8 s    = (r2 * rinv) * Hinv;
				const v2r8 ker  = kernel::W(s);

				rho = mj.madd(ker, rho);
			}

			rho *= (Hinv*Hinv)*(Hinv*v2r8(kernel::normalzeFactor));
			rho.storel(&dens[i+0].dens);
			dens[i+0].smth = PARAM::SMTH * pow(ep_i[i+0].mass / dens[i+0].dens, 1.0/(PS::F64)(PARAM::Dim));
			if(i+1 < Nip){
				rho.storeh(&dens[i+1].dens);
				dens[i+1].smth = PARAM::SMTH * pow(ep_i[i+1].mass / dens[i+1].dens, 1.0/(PS::F64)(PARAM::Dim));
			}
		}
	}
};

struct CalcDerivative{
	__attribute__ ((noinline))
	void operator () (
			const EPI::Drvt * const __restrict ep_i, 
			const PS::S32                      Nip,
		   	const EPJ::Drvt * const __restrict ep_j,
		   	const PS::S32                      Njp,
		   	RESULT::Drvt * const               drvt)
	{
		for(int i=0; i<Nip; i+=2){
			const v2r8 xi (ep_i[i+0].pos.x, ep_i[i+1].pos.x);
			const v2r8 yi (ep_i[i+0].pos.y, ep_i[i+1].pos.y);
			const v2r8 zi (ep_i[i+0].pos.z, ep_i[i+1].pos.z);
			const v2r8 vxi(ep_i[i+0].vel.x, ep_i[i+1].vel.x);
			const v2r8 vyi(ep_i[i+0].vel.y, ep_i[i+1].vel.y);
			const v2r8 vzi(ep_i[i+0].vel.z, ep_i[i+1].vel.z);
			const v2r8 hi (ep_i[i+0].smth,  ep_i[i+1].smth);
			const v2r8 Hinv = (hi * v2r8(kernel::supportRadius) ).rcpa_x8();
			const v2r8 eps2(1.e-16);

			v2r8 div(0.0);
			v2r8 rot_x(0.0);
			v2r8 rot_y(0.0);
			v2r8 rot_z(0.0);

			for(int j=0; j<Njp; j++){
				const double *ptr = (const double *)(ep_j + j);
				const v2r8 jp0 = v2r8::load(ptr + 0);
				const v2r8 jp1 = v2r8::load(ptr + 2);
				const v2r8 jp2 = v2r8::load(ptr + 4);
				const v2r8 jp3 = v2r8::load(ptr + 6);

				const v2r8_bcl mj(jp0);
				const v2r8_bch xj(jp0);
				const v2r8_bcl yj(jp1);
				const v2r8_bch zj(jp1);
				const v2r8_bcl vxj(jp2);
				const v2r8_bch vyj(jp2);
				const v2r8_bcl vzj(jp3);

				const v2r8 dx = xj - xi;
				const v2r8 dy = yj - yi;
				const v2r8 dz = zj - zi;
				const v2r8 dvx = vxj - vxi;
				const v2r8 dvy = vyj - vyi;
				const v2r8 dvz = vzj - vzi;

				const v2r8 r2   = ((eps2 + dx*dx) + dy*dy) + dz*dz;
				const v2r8 rinv = r2.rsqrta_x8();
				const v2r8 s    = (r2 * rinv) * Hinv;
				const v2r8 ker  = mj * kernel::gradW(s, rinv);

				div   -= ker * (dvx*dx + dvy*dy + dvz*dz);
				rot_x -= ker * (dvy*dz - dvz*dy);
				rot_y -= ker * (dvz*dx - dvx*dz);
				rot_z -= ker * (dvx*dy - dvy*dx);
			}
			const v2r8 rhoinv  = v2r8(ep_i[i+0].dens,  ep_i[i+1].dens).rcpa_x8();
			const v2r8 c = rhoinv * v2r8(kernel::normalzeFactor) * kernel::pow4(Hinv);
			div   *= c;
			rot_x *= c;
			rot_y *= c;
			rot_z *= c;

			div  .storel(&drvt[i+0].div_v);
			rot_x.storel(&drvt[i+0].rot_v.x);
			rot_y.storel(&drvt[i+0].rot_v.y);
			rot_z.storel(&drvt[i+0].rot_v.z);

			if(i+1 < Nip){
				div  .storeh(&drvt[i+1].div_v);
				rot_x.storeh(&drvt[i+1].rot_v.x);
				rot_y.storeh(&drvt[i+1].rot_v.y);
				rot_z.storeh(&drvt[i+1].rot_v.z);
			}
		}
	}
};

struct CalcHydroForce{
	__attribute__ ((noinline))
	void operator () (
			const EPI::Hydro * const __restrict  ep_i,
		   	const PS::S32                        Nip,
		   	const EPJ::Hydro * const __restrict  ep_j,
		   	const PS::S32                        Njp,
		   	RESULT::Hydro * const                hydro)
	{
#if 0
		enum {NJMAX = 65536, };
		assert(Njp <= NJMAX);
		struct{
			PS::F64 Hinv;
			PS::F64 P_ov_rho2;
		} jcache[NJMAX];
		for(int j=0; j<Njp ; j++){
			const EPJ::Hydro& jth = ep_j[j];
			jcache[j].Hinv = 1.0 / (kernel_t::supportRadius() * jth.smth);
			jcache[j].P_ov_rho2 = jth.pres / (jth.dens * jth.dens);
		}
#else
		for(int j=0; j<Njp ; j++){
			double *ptr = (double *)(ep_j + j);
			const EPJ::Hydro& jth = ep_j[j];
			const double Hinv = 1.0 / (kernel_t::supportRadius() * jth.smth);
			const double P_ov_rho2 = jth.pres / (jth.dens * jth.dens);
			ptr[8] = Hinv;
			ptr[9] = P_ov_rho2;
		}
#endif
		for(int i=0; i<Nip; i+=2){
			const v2r8 xi (ep_i[i+0].pos.x, ep_i[i+1].pos.x);
			const v2r8 yi (ep_i[i+0].pos.y, ep_i[i+1].pos.y);
			const v2r8 zi (ep_i[i+0].pos.z, ep_i[i+1].pos.z);
			const v2r8 vxi(ep_i[i+0].vel.x, ep_i[i+1].vel.x);
			const v2r8 vyi(ep_i[i+0].vel.y, ep_i[i+1].vel.y);
			const v2r8 vzi(ep_i[i+0].vel.z, ep_i[i+1].vel.z);

			const v2r8 idens(ep_i[i+0].dens, ep_i[i+1].dens);
			const v2r8 ipres(ep_i[i+0].pres, ep_i[i+1].pres);
			const v2r8 ismth(ep_i[i+0].smth, ep_i[i+1].smth);
			const v2r8 isnds(ep_i[i+0].snds, ep_i[i+1].snds);
			const v2r8 iBal (ep_i[i+0].Bal,  ep_i[i+1].Bal);


			const v2r8 hi (ep_i[i+0].smth,  ep_i[i+1].smth);
			const v2r8 iHinv  = (hi * v2r8(kernel::supportRadius) ).rcpa_x8();
			const v2r8 iHinv4 = kernel::pow4(iHinv);
			const v2r8 iP_ov_rho2 = ipres * (idens * idens).rcpa_x8();
			const v2r8 eps2(1.e-16);

			v2r8 v_sig_max(0.0);
			v2r8 acc_x(0.0), acc_y(0.0), acc_z(0.0), edot(0.0);

			for(int j=0; j<Njp; j++){
				const double *ptr = (const double *)(ep_j + j);
				const v2r8 jp0 = v2r8::load(ptr + 0); // x, y
				const v2r8 jp1 = v2r8::load(ptr + 2); // z, vx
				const v2r8 jp2 = v2r8::load(ptr + 4); // vy, vz
				const v2r8 jp3 = v2r8::load(ptr + 6); // dens, mass
				// const v2r8 jp4 = v2r8::load(ptr + 8); // smth, pres
				const v2r8 jp5 = v2r8::load(ptr + 10); // snds, Bal
#if 0
				const v2r8 jp6 = v2r8::load((double *)&jcache[j]);
#else
				const v2r8 jp6 = v2r8::load(ptr + 8);
#endif

				const v2r8_bcl xj (jp0);
				const v2r8_bch yj (jp0);
				const v2r8_bcl zj (jp1);
				const v2r8_bch vxj(jp1);
				const v2r8_bcl vyj(jp2);
				const v2r8_bch vzj(jp2);

				const v2r8_bcl jdens(jp3);
				const v2r8_bch jmass(jp3);
				// jsmth, jpres not used


				const v2r8_bcl jsnds(jp5);
				const v2r8_bch jBal (jp5);

				const v2r8_bcl jHinv     (jp6);
				const v2r8_bch jP_ov_rho2(jp6);

				const v2r8 dx = xj - xi;
				const v2r8 dy = yj - yi;
				const v2r8 dz = zj - zi;
				const v2r8 dvx = vxj - vxi;
				const v2r8 dvy = vyj - vyi;
				const v2r8 dvz = vzj - vzi;

				const v2r8 r2   = ((eps2 + dx*dx) + dy*dy) + dz*dz;
				const v2r8 rinv = r2.rsqrta_x8();
				const v2r8 r1   = rinv * r2;
				const v2r8 si   = r1 * iHinv;
				const v2r8 sj   = r1 * jHinv;

				const v2r8 dvdr = dvx*dx + dvy*dy + dvz*dz;
				const v2r8 w_ij = v2r8::min(v2r8(0.0), rinv*dvdr);
				const v2r8 v_sig = (isnds + jsnds) - v2r8(3.0)*w_ij;
				v_sig_max = v2r8::max(v_sig_max, v_sig);

				#if 1 //With Balsara Switch
				// const v2r8 AV = v2r8(-0.5) * ((v_sig * w_ij) * (idens + jdens).rcpa_x8() * (iBal + jBal));
				const v2r8 AV = ((v_sig * w_ij) * (idens + jdens).rcpa_x8() * (iBal + jBal));
				#else
				// const v2r8 AV = (-v_sig * w_ij) * (idens + jdens).rcpa_x8;
				const v2r8 AV = v2r8(-2.0) * (v_sig * w_ij) * (idens + jdens).rcpa_x8;
				#endif

				const v2r8 gradWs = jmass * 
					(kernel::gradW(si, rinv)*iHinv4 + kernel::gradW(sj, rinv)*kernel::pow4(v2r8(jHinv)));

				// const v2r8 cacc = (gradWs) * ((iP_ov_rho2 + jP_ov_rho2) + AV);
				const v2r8 cacc = (gradWs) * ((iP_ov_rho2 + jP_ov_rho2) - v2r8(0.5)*AV);
				acc_x += cacc * dx;
				acc_y += cacc * dy;
				acc_z += cacc * dz;
				edot  += (gradWs) * (iP_ov_rho2 - v2r8(0.25)*AV) * dvdr;
			}
			const v2r8 norm(0.5 * kernel::normalzeFactor);
			acc_x *= norm;
			acc_y *= norm;
			acc_z *= norm;
			edot  *= norm;
			const v2r8 dt = v2r8(2.0 * PARAM::C_CFL) * ismth * v_sig_max.rcpa_x8();

			acc_x.storel(&hydro[i+0].acc.x);
			acc_y.storel(&hydro[i+0].acc.y); 
			acc_z.storel(&hydro[i+0].acc.z);
			edot .storel(&hydro[i+0].eng_dot);
			dt   .storel(&hydro[i+0].dt  );
			if(i+1 < Nip){
				acc_x.storeh(&hydro[i+1].acc.x);
				acc_y.storeh(&hydro[i+1].acc.y);
				acc_z.storeh(&hydro[i+1].acc.z);
				edot .storeh(&hydro[i+1].eng_dot);
				dt   .storeh(&hydro[i+1].dt);
			}
		}
	}
};

template <class TPtclJ> 
struct CalcGravityForce{
	void operator () (const EPI::Grav* const __restrict ep_i, const PS::S32 Nip, const TPtclJ* const __restrict ep_j, const PS::S32 Njp, RESULT::Grav* const grav){
		for(PS::S32 i = 0; i < Nip ; ++ i){
			grav[i].clear();
			const EPI::Grav& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const TPtclJ& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64 dr2 = dr * dr;
				const PS::F64 dr_inv = 1.0 / sqrt(dr2 + ith.getEps2());
				const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
				grav[i].acc -= m_dr3_inv * dr;
				grav[i].pot -= jth.mass * dr_inv;
			}
		}
	}
};
