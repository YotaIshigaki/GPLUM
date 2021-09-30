#pragma once

#include "phantomquad.hpp"

namespace kernel{
	static const double pi = 4.0 * atan(1.0);
	static const double supportRadius = 2.5;
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
	} // 2sub 5add 12mul 1max ()
}

class CalcDensity{
	kernel_t kernel;
	public:
	//__attribute__ ((noinline))
	void operator () (const EPI::Dens* const __restrict ep_i,
                      const PS::S32 Nip,
                      const EPJ::Dens* const __restrict ep_j,
                      const PS::S32 Njp,
                      RESULT::Dens* const __restrict dens){
		const v2r8 dens_min(5.0);//Genda+ 2012
		for(int i=0; i<Nip; i+=2){
			const v2r8 xi(ep_i[i+0].pos.x, ep_i[i+1].pos.x);
			const v2r8 yi(ep_i[i+0].pos.y, ep_i[i+1].pos.y);
			const v2r8 zi(ep_i[i+0].pos.z, ep_i[i+1].pos.z);
			v2r8 hi(ep_i[i+0].smth,  ep_i[i+1].smth);
			v2r8 Hinv = (hi * v2r8(kernel::supportRadius) ).rcpa_x8();
			const v2r8 eps2(1.e-16);

            dens[i+0].itr = false;
            if(i+1 < Nip){
                dens[i+1].itr = false;
            }
            for(int loop =0; loop<3; loop++){
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
                    const v2r8 dz = zj - zi; // 3sub (3OP)

                    const v2r8 r2   = ((eps2 + dx*dx) + dy*dy) + dz*dz; // 3sub 3mul (9OP)
                    const v2r8 rinv = r2.rsqrta_x8(); // 6add 10mul 1rsqrt (26OP)
                    const v2r8 s    = (r2 * rinv) * Hinv; // 2mul (28OP)
                    const v2r8 ker  = kernel::W(s); // 1sub 3add 7mul 1max (40OP)

                    rho = mj.madd(ker, rho); // 1add 1mul (42OP)
                } // 42OP * 3 = 123OP
                 

                rho *= (Hinv*Hinv)*(Hinv*v2r8(kernel::normalzeFactor));
				rho = v2r8::max(rho, dens_min);

                rho.storel(&dens[i+0].dens);
                dens[i+0].smth = PARAM::SMTH * pow(ep_i[i+0].mass / dens[i+0].dens, 1.0/(PS::F64)(PARAM::Dim));
                dens[i+0].itr = (dens[i+0].smth*kernel_t::supportRadius() < ep_i[i+0].Rsearch) ? dens[i+0].itr : true;
                if(i+1 < Nip){
                    rho.storeh(&dens[i+1].dens);
                    dens[i+1].smth = PARAM::SMTH * pow(ep_i[i+1].mass / dens[i+1].dens, 1.0/(PS::F64)(PARAM::Dim));
                    dens[i+1].itr = (dens[i+0].smth*kernel_t::supportRadius() < ep_i[i+1].Rsearch) ? dens[i+1].itr : true;
                }
                hi = v2r8(dens[i+0].smth, dens[i+1].smth);
                //hi = __builtin_fj_set_v2r8(dens[i+0].smth, dens[i+1].smth);
                Hinv = (hi * v2r8(kernel::supportRadius) ).rcpa_x8();
                //rho = __builtin_fj_set_v2r8(0.0, 0.0);
                rho = v2r8(0.0, 0.0);
            }

// calc drvt
            const v2r8 vxi(ep_i[i+0].vel.x, ep_i[i+1].vel.x);
            const v2r8 vyi(ep_i[i+0].vel.y, ep_i[i+1].vel.y);
            const v2r8 vzi(ep_i[i+0].vel.z, ep_i[i+1].vel.z);
            //hi = __builtin_fj_set_v2r8(dens[i+0].smth, dens[i+1].smth);
			//Hinv = (hi * v2r8(kernel::supportRadius) ).rcpa_x8();
            
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
                const v2r8_bcl vzj(jp3); // 

                const v2r8 dx = xj - xi;
                const v2r8 dy = yj - yi;
                const v2r8 dz = zj - zi;
                const v2r8 dvx = vxj - vxi;
                const v2r8 dvy = vyj - vyi;
                const v2r8 dvz = vzj - vzi; // 6sub (6op)

                const v2r8 r2   = ((eps2 + dx*dx) + dy*dy) + dz*dz; // 3add 3mul (9sub 3mul 12op)
                const v2r8 rinv = r2.rsqrta_x8(); // 3add 3sub 10mul 1rsqrt (15sub 13mul 1rsqrt 29op)
                const v2r8 s    = (r2 * rinv) * Hinv; // 0add 2mul (15sub 15mul 1rsqrt 31op)
                const v2r8 ker  = mj * kernel::gradW(s, rinv); // 2sub 5add 13mul 1max (22sub 28mul 1rsqrt 1max 47op)

                div   -= ker * (dvx*dx + dvy*dy + dvz*dz); 
                rot_x -= ker * (dvy*dz - dvz*dy);
                rot_y -= ker * (dvz*dx - dvx*dz);
                rot_z -= ker * (dvx*dy - dvy*dx); // 7sub 2 add 13mul (31sub 41mul 1rsqrt 1max 74op)
            } // total 58 + 129 = 187op
            
            const v2r8 rhoinv  = v2r8(dens[i+0].dens,  dens[i+1].dens).rcpa_x8();
            const v2r8 c = rhoinv * v2r8(kernel::normalzeFactor) * kernel::pow4(Hinv);
            div   *= c;
            rot_x *= c;
            rot_y *= c;
            rot_z *= c;

            div  .storel(&dens[i+0].div_v);
            rot_x.storel(&dens[i+0].rot_v.x);
            rot_y.storel(&dens[i+0].rot_v.y);
            rot_z.storel(&dens[i+0].rot_v.z);
            
            if(i+1 < Nip){
                div  .storeh(&dens[i+1].div_v);
                rot_x.storeh(&dens[i+1].rot_v.x);
                rot_y.storeh(&dens[i+1].rot_v.y);
                rot_z.storeh(&dens[i+1].rot_v.z);
            }
        }
    } // end of class

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
				const v2r8 dvz = vzj - vzi; // 6sub (6op)

				const v2r8 r2   = ((eps2 + dx*dx) + dy*dy) + dz*dz; // 3add 3mul (9add 3mul 12op)
				const v2r8 rinv = r2.rsqrta_x8(); // 6add 10mul 1rsqrt (15add 13mul 1rsqrt 29op)
				const v2r8 r1   = rinv * r2; // 1mul (15add 14mul 1rsqrt 30op)
				const v2r8 si   = r1 * iHinv; // 1mul (15add 15mul 1rsqrt 31op)
				const v2r8 sj   = r1 * jHinv; // 1mul (15add 16mul 1rsqrt 32op)

				const v2r8 dvdr = dvx*dx + dvy*dy + dvz*dz; // 2add 3mul (17add 19mul 1rsqrt 37op)
				const v2r8 w_ij = v2r8::min(v2r8(0.0), rinv*dvdr); // 1mul 1min (17add 20mul 1rsqrt 1min 39op)
				const v2r8 v_sig = isnds + jsnds - v2r8(3.0)*w_ij; // 1sub 1 add 1 mul(19add 21mul 1rsqrt 1min 42op)
				v_sig_max = v2r8::max(v_sig_max, v_sig); // 1max (19add 21mul 1rsqrt 2min 43op)

				#if 1 //With Balsara Switch
				// const v2r8 AV = v2r8(-0.5) * ((v_sig * w_ij) * (idens + jdens).rcpa_x8() * (iBal + jBal));
				const v2r8 AV = ((v_sig * w_ij) * (idens + jdens).rcpa_x8() * (iBal + jBal)); 
                // (rcpa_x8: 1rcpa 1sub 3add 7mul (12op))
                // 1rcpa 1sub 5add 10mul (25add 31mul 2rsqrt 2min 50op)
				#else
				// const v2r8 AV = (-v_sig * w_ij) * (idens + jdens).rcpa_x8;
				const v2r8 AV = v2r8(-2.0) * (v_sig * w_ij) * (idens + jdens).rcpa_x8;
				#endif

				const v2r8 gradWs = jmass * 
					(kernel::gradW(si, rinv)*iHinv4 + kernel::gradW(sj, rinv)*kernel::pow4(v2r8(jHinv)));  
                // (gradW: 2sub 5add 12mul 1max (20op))
                // 4sub 11add 27mul 2max (40add 58mul 2rsqrt 4min 104op)

				// const v2r8 cacc = (gradWs) * ((iP_ov_rho2 + jP_ov_rho2) + AV);
				const v2r8 cacc = (gradWs) * ((iP_ov_rho2 + jP_ov_rho2) - v2r8(0.5)*AV); // 1sub 1add 2mul (42add 60mul 2rsqrt 4min 108op)
				acc_x += cacc * dx;
				acc_y += cacc * dy;
				acc_z += cacc * dz;
				edot  += (gradWs) * (iP_ov_rho2 - v2r8(0.25)*AV) * dvdr; // 1sub 4add 6mul (47add 66mul 2rsqrt 4min 119op)
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


template <class TPtclJ> class CalcGravityForce{
	public:
	void operator () (const EPI::Grav* const __restrict ep_i, const PS::S32 Nip, const TPtclJ* const __restrict ep_j, const PS::S32 Njp, RESULT::Grav* const grav){
		//const PS::F64 eps2 = (0.01 * 6400e+3) * (0.01 * 6400e+3);
		const PS::F64 eps2 = ep_i[0].getEps2();
		PhantomGrapeQuad pg;
		pg.set_eps2(eps2);
		for(PS::S32 i = 0 ; i < Nip ; ++ i){
			const PS::F64vec pos_i = ep_i[i].getPos();
			pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
		}
		for(PS::S32 i = 0 ; i < Njp ; ++ i){
			const PS::F64 m_j = ep_j[i].getCharge();
			const PS::F64vec pos_j = ep_j[i].getPos();
			pg.set_epj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j);
		}
		pg.run_epj(Nip, Njp);
		for(PS::S32 i = 0 ; i < Nip ; ++ i){
			PS::F64 * p = &(grav[i].pot);
			PS::F64 * a = (PS::F64 * )(&grav[i].acc[0]);
			pg.get_accp_one(i, a[0], a[1], a[2], *p);
		}
	}
};

