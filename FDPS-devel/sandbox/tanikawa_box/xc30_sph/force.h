#pragma once

//#include "phantomquad.hpp"
#include "phantomquad_x86.hpp"

class CalcDensity{
	kernel_t kernel;
	public:
	//__attribute__ ((noinline))
	void operator () (const EPI::Dens* const __restrict ep_i, const PS::S32 Nip, const EPJ::Dens* const __restrict ep_j, const PS::S32 Njp, RESULT::Dens* const __restrict dens){
        //std::cout<<"calc dens"<<std::endl;
        for(PS::S32 i = 0 ; i < Nip ; ++ i){
            const EPI::Dens& ith = ep_i[i];
            const PS::F64 r_search_i = ep_i[i].getRSearch();
            //const PS::F64 Hinv = 1.0 / (kernel_t::supportRadius() * ith.smth);
            PS::F64 Hinv = 1.0 / (kernel_t::supportRadius() * ith.smth);
            bool itr = false;
            PS::F64 smth = 0.0;
            PS::F64 rhoi = 0.0;
            for(PS::S32 j = 0 ; j < Njp ; ++ j){
                const EPJ::Dens& jth = ep_j[j];
                const PS::F64vec dr = jth.pos - ith.pos;
                const PS::F64 dr1  = sqrt(dr * dr);
                const PS::F64 s    = dr1 * Hinv;
                rhoi += jth.mass * kernel_t::W(s, Hinv);
            }
            rhoi = std::max(rhoi, 5.0);
            smth = PARAM::SMTH * pow(ith.mass / rhoi, 1.0/(PS::F64)(PARAM::Dim));
            itr = (smth*kernel_t::supportRadius() < r_search_i) ? itr : true;

            rhoi = 0.0;
            Hinv = 1.0 / (kernel_t::supportRadius() * smth);
            for(PS::S32 j = 0 ; j < Njp ; ++ j){
                const EPJ::Dens& jth = ep_j[j];
                const PS::F64vec dr = jth.pos - ith.pos;
                const PS::F64 dr1  = sqrt(dr * dr);
                const PS::F64 s    = dr1 * Hinv;
                rhoi += jth.mass * kernel_t::W(s, Hinv);
            }
            rhoi = std::max(rhoi, 5.0);
            smth = PARAM::SMTH * pow(ith.mass / rhoi, 1.0/(PS::F64)(PARAM::Dim));
            itr = (smth*kernel_t::supportRadius() < r_search_i) ? itr : true;

            rhoi = 0.0;
            Hinv = 1.0 / (kernel_t::supportRadius() * smth);
            for(PS::S32 j = 0 ; j < Njp ; ++ j){
                const EPJ::Dens& jth = ep_j[j];
                const PS::F64vec dr = jth.pos - ith.pos;
                const PS::F64 dr1  = sqrt(dr * dr);
                const PS::F64 s    = dr1 * Hinv;
                rhoi += jth.mass * kernel_t::W(s, Hinv);
            }
            rhoi = std::max(rhoi, 5.0);
            smth = PARAM::SMTH * pow(ith.mass / rhoi, 1.0/(PS::F64)(PARAM::Dim));
            itr = (smth*kernel_t::supportRadius() < r_search_i) ? itr : true;

            dens[i].dens = rhoi;
            dens[i].smth = smth;
            dens[i].itr = itr;

			//Hinv = 1.0 / (kernel_t::supportRadius() * ith.smth);
            Hinv = 1.0 / (kernel_t::supportRadius() * smth);
			PS::F64    div(0.0);
			PS::F64vec rot(0.0);
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const EPJ::Dens& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;
				//const PS::F64 dr2  = kernel_t::safe_norm2(dr, ith.smth);
                const PS::F64 dr2  = kernel_t::safe_norm2(dr, smth);
				const PS::F64 rinv  = 1.0 / sqrt(dr2);
				const PS::F64 dr1   = dr2 * rinv;
				const PS::F64 s     = dr1 * Hinv;
				const PS::F64 gWs = kernel_t::gradWs(s, Hinv, rinv);
				div -= (jth.mass * gWs) * (dv * dr);
				rot -= (jth.mass * gWs) * (dv ^ dr);
			}

			//const PS::F64 rhoinv = 1.0 / ith.dens;
            const PS::F64 rhoinv = 1.0 / rhoi;
			dens[i].div_v = rhoinv * div;
			dens[i].rot_v = rhoinv * rot;

/*
            for(PS::S32 j = 0 ; j < Njp ; ++ j){
                const EPJ::Dens& jth = ep_j[j];
                const PS::F64vec dr = jth.pos - ith.pos;
                //const PS::F64 dr1  = sqrt(dr * dr);
                //const PS::F64 s    = dr1 * Hinv;
				const PS::F64 dr2  = kernel_t::safe_norm2(dr, ith.smth);
				const PS::F64 rinv  = 1.0 / sqrt(dr2);
				const PS::F64 dr1   = dr2 * rinv;
				const PS::F64 s     = dr1 * Hinv;
                rhoi += jth.mass * kernel_t::W(s, Hinv);
				const PS::F64vec dv = ith.vel - jth.vel;
				const PS::F64 gWs = kernel_t::gradWs(s, Hinv, rinv);
				div -= (jth.mass * gWs) * (dv * dr);
				rot -= (jth.mass * gWs) * (dv ^ dr);
            }
            rhoi = std::max(rhoi, 5.0);
            smth = PARAM::SMTH * pow(ith.mass / rhoi, 1.0/(PS::F64)(PARAM::Dim));
            itr = (smth*kernel_t::supportRadius() < r_search_i) ? itr : true;
			const PS::F64 rhoinv = 1.0 / ith.dens;

            dens[i].dens = rhoi;
            dens[i].smth = smth;
            dens[i].itr = itr;
			dens[i].div_v = rhoinv * div;
			dens[i].rot_v = rhoinv * rot;
*/
        }
	}
};


/*
class CalcDerivative{
	kernel_t kernel;
	public:
	//__attribute__ ((noinline))
	void operator () (const EPI::Drvt* const __restrict ep_i, const PS::S32 Nip, const EPJ::Drvt* const __restrict ep_j, const PS::S32 Njp, RESULT::Drvt* const drvt){
		for(PS::S32 i = 0; i < Nip ; ++ i){
#if 0
			drvt[i].clear();
			const EPI::Drvt& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const EPJ::Drvt& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;
				const PS::F64vec gW = kernel.gradW(dr, ith.smth);
				drvt[i].div_v += - jth.mass * (dv * gW);
				drvt[i].rot_v += - jth.mass * (dv ^ gW);
			}
			const PS::F64 rhoinv = 1.0 / ith.dens;
			drvt[i].div_v *= rhoinv;
			drvt[i].rot_v *= rhoinv;
#else
			const EPI::Drvt& ith = ep_i[i];
			const PS::F64 Hinv = 1.0 / (kernel_t::supportRadius() * ith.smth);
			// drvt[i].clear();
			PS::F64    div(0.0);
			PS::F64vec rot(0.0);
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const EPJ::Drvt& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;
				const PS::F64 dr2  = kernel_t::safe_norm2(dr, ith.smth);
				const PS::F64 rinv  = 1.0 / sqrt(dr2);
				// const PS::F64 rinv2 = rinv * rinv;
				const PS::F64 dr1   = dr2 * rinv;
				const PS::F64 s     = dr1 * Hinv;

#if 0
				// const PS::F64vec gW = kernel.gradW(dr, ith.smth);
				const PS::F64vec gW = kernel.gradWs(s, Hinv, rinv) * dr;
				div -= jth.mass * (dv * gW);
				rot -= jth.mass * (dv ^ gW);
#else
				const PS::F64 gWs = kernel_t::gradWs(s, Hinv, rinv);
				div -= (jth.mass * gWs) * (dv * dr);
				rot -= (jth.mass * gWs) * (dv ^ dr);
#endif
			}
			const PS::F64 rhoinv = 1.0 / ith.dens;
			drvt[i].div_v = rhoinv * div;
			drvt[i].rot_v = rhoinv * rot;
#endif
		}
	}
};
*/


class CalcHydroForce{
	const kernel_t kernel;
	public:
	//__attribute__ ((noinline))
	void operator () (const EPI::Hydro* const  __restrict ep_i, const PS::S32 Nip, const EPJ::Hydro* const __restrict ep_j, const PS::S32 Njp, RESULT::Hydro* const hydro){
#if 0
		for(PS::S32 i = 0; i < Nip ; ++ i){
			hydro[i].clear();
			PS::F64 v_sig_max = 0.0;
			const EPI::Hydro& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const EPJ::Hydro& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;

				PS::F64 w_ij = (dv * dr) / sqrt(dr * dr);
				w_ij = w_ij < 0.0 ? w_ij : 0.0;

				const PS::F64 v_sig = ith.snds + jth.snds - 3.0 * w_ij;
				v_sig_max = v_sig_max >  v_sig ? v_sig_max : v_sig;
				#if 1 //With Balsara Switch
				const PS::F64 AV = -0.5 * v_sig * w_ij / (ith.dens + jth.dens) * (ith.Bal + jth.Bal);
				//if(Bal_i > 0.1 && ep_i[i].pos.x < 0.6)std::cout << ep_i[i].pos.x << ", " << Bal_i << std::endl;
				#else
				const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ith.dens + jth.dens));
				#endif
				const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ith.smth) + kernel.gradW(dr, jth.smth));
				const PS::F64 iP_ov_rho2 = ith.pres / (ith.dens * ith.dens); // j-loop invariant
				const PS::F64 jP_ov_rho2 = jth.pres / (jth.dens * jth.dens);
				hydro[i].acc     -= (jth.mass * (iP_ov_rho2 + jP_ov_rho2 + AV)) * gradW;
				hydro[i].eng_dot += jth.mass * (iP_ov_rho2 + 0.5 * AV) * (dv * gradW);
			}
			hydro[i].dt = PARAM::C_CFL * 2.0 * ith.smth / v_sig_max;
		}
#else
		//enum {NJMAX = 65536, };
        enum {NJMAX = 131072, };
		assert(Njp <= NJMAX);
		struct{
			PS::F64 Hinv;
			PS::F64 P_ov_rho2;
		} jcache[NJMAX];
		for(PS::S32 j = 0; j < Njp ; ++ j){
			const EPJ::Hydro& jth = ep_j[j];
			jcache[j].Hinv = 1.0 / (kernel_t::supportRadius() * jth.smth);
			jcache[j].P_ov_rho2 = jth.pres / (jth.dens * jth.dens);
		}

		for(PS::S32 i = 0; i < Nip ; ++ i){
			// hydro[i].clear();
			PS::F64vec acc    = 0.0;
			PS::F64 edot      = 0.0;
			PS::F64 v_sig_max = 0.0;

			const EPI::Hydro& ith = ep_i[i];
			const PS::F64 iP_ov_rho2 = ith.pres / (ith.dens * ith.dens); // j-loop invariant
			const PS::F64 iHinv = 1.0 / (kernel_t::supportRadius() * ith.smth);
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const EPJ::Hydro& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;

				const PS::F64 dr2  = kernel_t::safe_norm2(dr, ith.smth);
				const PS::F64 rinv  = 1.0 / sqrt(dr2);
				// const PS::F64 rinv2 = rinv * rinv;
				const PS::F64 dr1   = dr2 * rinv;
				const PS::F64 si    = dr1 * iHinv;
				const PS::F64 jHinv = jcache[j].Hinv;
				const PS::F64 sj    = dr1 * jHinv;
				const PS::F64 dvdr  = dv * dr;

				PS::F64 w_ij = dvdr * rinv;
				w_ij = w_ij < 0.0 ? w_ij : 0.0;

				const PS::F64 v_sig = ith.snds + jth.snds - 3.0 * w_ij;
				v_sig_max = v_sig_max >  v_sig ? v_sig_max : v_sig;


				#if 1 //With Balsara Switch
				const PS::F64 AV = -0.5 * v_sig * w_ij / (ith.dens + jth.dens) * (ith.Bal + jth.Bal);
				#else
				const PS::F64 AV = -v_sig * w_ij / (ith.dens + jth.dens);
				#endif

				// const PS::F64vec gradW =  0.5 *  (kernel_t::gradW(si, iHinv, dr, rinv2) + kernel_t::gradW(sj, jHinv, dr, rinv2));
				const PS::F64 gradWs = (kernel_t::gradWs(si, iHinv, rinv) + kernel_t::gradWs(sj, jHinv, rinv));
				// const PS::F64 iP_ov_rho2 = ith.pres / (ith.dens * ith.dens); // j-loop invariant
				const PS::F64 jP_ov_rho2 = jcache[j].P_ov_rho2;
				acc  -= ((jth.mass * gradWs) * (iP_ov_rho2 + jP_ov_rho2 + AV)) * dr;
				edot +=  (jth.mass * gradWs) * (iP_ov_rho2 + 0.5 * AV) * dvdr;
			}
			hydro[i].acc     = 0.5 * acc;
			hydro[i].eng_dot = 0.5 * edot;
			hydro[i].dt      = PARAM::C_CFL * 2.0 * ith.smth / v_sig_max;
		}
#endif
	}
};

template <class TPtclJ> class CalcGravityForce{
public:
	void operator () (const EPI::Grav* const __restrict ep_i, const PS::S32 Nip, const TPtclJ* const __restrict ep_j, const PS::S32 Njp, RESULT::Grav* const grav){
#if 0
		for(PS::S32 i = 0; i < Nip ; ++ i){
			grav[i].clear();
			const EPI::Grav& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const TPtclJ& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64 dr2 = dr * dr;
				const PS::F64 dr_inv = 1.0 / sqrt(dr2 + ith.getEps2());
				const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
#if 0
				grav[i].acc += - m_dr3_inv * dr;
				grav[i].pot += - jth.mass * dr_inv;
#else
				grav[i].acc -= m_dr3_inv * dr;
				grav[i].pot -= jth.mass * dr_inv;
#endif
			}
		}
#else
        const PS::F64 eps2 = (0.01 * 6400e+3) * (0.01 * 6400e+3);
        static __thread PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        pg.check_range(Nip, Njp);
        for(PS::S32 i=0; i<Nip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
        }
        for(PS::S32 i=0; i<Njp; i++){
            const PS::F64 m_j = ep_j[i].getCharge();
            const PS::F64vec pos_j = ep_j[i].getPos();
            pg.set_epj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j);
        }
        pg.run_epj(Nip, Njp);
        for(PS::S32 i=0; i<Nip; i++){
            PS::F64 * p = &(grav[i].pot);
            PS::F64 * a = (PS::F64 * )(&grav[i].acc[0]);
            pg.get_accp_one(i, a[0], a[1], a[2], *p);
        }

		#endif
	}
};
