#pragma once
#include "math.h"
#include "param.h"
#include "EoS.h"
#ifdef __AVX__
#include "avx.h"
#include <x86intrin.h>
#endif

PS::F32 get_dtime(void);

//Wendland C6
struct kernel_t{
	kernel_t(){
	}
	//W
	PS::F32 W(const PS::F32vec dr, const PS::F32 h) const{
		const PS::F32 H = supportRadius() * h;
		const PS::F32 s = sqrt(dr * dr) / H;
		PS::F32 r_value;
		r_value = (1.0f + s * (8.0f + s * (25.0f + s * (32.0f)))) * math::pow8(math::plus(1.0f - s));
		r_value *= (1365.f/64.f) / (H * H * H * math::pi);
		return r_value;
	}
	//gradW
	PS::F32vec gradW(const PS::F32vec dr, const PS::F32 h) const{
		const PS::F32 H = supportRadius() * h;
		const PS::F32 s = sqrt(dr * dr) / H;
		PS::F32 r_value;
		r_value = math::pow7(math::plus(1.0f - s)) * (math::plus(1.0f - s) * (8.0f + s * (50.0f + s * (96.0f))) - 8.0f * (1.0f + s * (8.0f + s * (25.0f + s * (32.0f)))));
		r_value *= (1365.f/64.f) / (H * H * H * math::pi);
		return dr * r_value / (sqrt(dr * dr) * H  + 0.01f * h);
	}
	#ifdef __AVX__
	AVXf W(const AVXf r, const AVXf h) const{
		const AVXf H = AVXf(supportRadius()) * h;
		const AVXf s = r / H;
		AVXf r_value;
		r_value = (AVXf(1.0f) + s * (AVXf(8.0f) + s * (AVXf(25.0f) + s * AVXf(32.0f)))) * (AVXf(1.0f - s).max(0.0f)).pow(8);
		r_value = r_value * AVXf(1365.f/64.f) / (H * H * H * AVXf(math::pi));
		return r_value;
	}
	AVXf abs_gradW(const AVXf r, const AVXf h) const{
		const AVXf H = AVXf(supportRadius()) * h;
		const AVXf s = r / H;
		AVXf r_value;
		r_value = (AVXf(1.0f - s).max(0.0f)).pow(7) * AVXf(1.0f - s).max(0.0f) * (AVXf(8.0f) + s * (AVXf(50.0f) + s * AVXf(96.0f))) - AVXf(8.0f) * (AVXf(1.0f) + s * (AVXf(8.0f) + s * (AVXf(25.0f) + s * (AVXf(32.0f)))));
		r_value = r_value * AVXf(1365.f/64.f) / (H * H * H * AVXf(math::pi));
		return r_value / H;
	}
	#endif
	#ifdef ENABLE_KNL
	__m512 W(const __m512 r, const __m512 h) const{
		const __m512 H = _mm512_mul_ps(h, _mm512_set1_ps(supportRadius()));
		const __m512 s = _mm512_div_ps(s, H);
		__m512 r_value = _mm512_set1_ps(32.0f);
		r_value = _mm512_fmadd_ps(r_value, s, _mm512_set1_ps(25.0f)); /* 25.0f + s * 32.0f */;
		r_value = _mm512_fmadd_ps(r_value, s, _mm512_set1_ps(8.0f));  /*  8.0f + s * r     */;
		r_value = _mm512_fmadd_ps(r_value, s, _mm512_set1_ps(1.0f));  /*  1.0f + s * r     */;
		__m512 tmp = _mm512_pow_ps(_mm512_max_ps(_mm512_sub_ps(_mm512_set1_ps(1.0), s), _mm512_set1_ps(0.0f)), _mm512_set1_ps(8.0f));
		return _mm512_mul_ps(tmp, r_value);
	}
	#endif
	static PS::F32 supportRadius(){
		return 3.5f;
	}
};

namespace RESULT{
	//Density summation
	class Dens{
		public:
		PS::F32 dens;
		void clear(){
			dens = 0;
		}
	};
	//for Balsara switch
	class Drvt{
		public:
		PS::F32 div_v;
		PS::F32vec rot_v;
		PS::F32 grad_smth;
		void clear(){
			div_v = 0.0;
			grad_smth = 0.0;
			rot_v = 0.0;
		}
	};
	//Hydro force
	class Hydr{
		public:
		PS::F32vec acc;
		PS::F32 eng_dot;
		PS::F32 dt;
		void clear(){
			acc = 0;
			eng_dot = 0;
			dt = 1.0e+30;
		}
	};
	//Self gravity
	class Grav{
		public:
		PS::F32vec acc;
		PS::F32 pot;
		PS::F32 dt;
		void clear(){
			acc = 0.0;
			pot = 0.0;
			dt = 1.0e+30;
		}
	};
}

class RealPtcl{
	public:
	PS::F32 mass;
	PS::F32vec pos, vel, acc, grav;
	PS::F32 dens;//DENSity
	PS::F32 eng; //ENerGy
	PS::F32 pres;//PRESsure
	PS::F32 smth;//SMooTHing length
	PS::F32 snds; //SouND Speed
	PS::F32 div_v;
	PS::F32vec rot_v;
	PS::F32 pot;
	PS::F32 Bal; //Balsala switch
	PS::F32 grad_smth;

	PS::F32 eng_dot;
	PS::F32vec vel_half;
	PS::F32 eng_half;
	PS::F32 dt;
	PS::S64 id;
	PS::S64 mat;//material
	const EoS::EoS_t<PS::F32>* EoS;

	//Copy functions
	void copyFromForce(const RESULT::Dens& dens){
		this->dens = dens.dens;
	}
	void copyFromForce(const RESULT::Drvt& drvt){
		this->div_v = drvt.div_v;
		this->rot_v = drvt.rot_v;
		this->Bal = std::abs(drvt.div_v) / (std::abs(drvt.div_v) + sqrt(drvt.rot_v * drvt.rot_v) + 1.0e-4 * this->snds / this->smth); //Balsala switch
		this->grad_smth = drvt.grad_smth;
	}
	void copyFromForce(const RESULT::Hydr& force){
		this->acc     = force.acc;
		this->eng_dot = force.eng_dot;
		this->dt      = force.dt;
	}
	void copyFromForce(const RESULT::Grav& force){
		this->grav = force.acc;
		this->pot  = force.pot;
		//not copy dt
	}
	//Give necessary values to FDPS
	PS::F32 getCharge() const{
		return this->mass;
	}
	PS::F32vec getPos() const{
		return this->pos;
	}
	PS::F32 getRSearch() const{
		return kernel_t::supportRadius() * this->smth;
	}
	void setPos(const PS::F32vec& pos){
		this->pos = pos;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", id, mat, mass, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, dens, eng, pres, pot);
	}
	void readAscii(FILE* fp){
		fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &mat, &mass, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z, &dens, &eng, &pres, &pot);
	}
	void setPressure(void){
		pres = EoS->Pressure(dens, eng);
		snds = EoS->SoundSpeed(dens, eng);
	}
	void initialize(void){
		smth = PARAM::SMTH * powf(mass / dens, 1.0/PARAM::Dim);
	}
};

namespace EPI{
	class Dens{
		public:
		PS::F32vec pos;
		PS::F32    mass;
		PS::F32    smth;
		PS::S64    id;
		void copyFromFP(const RealPtcl& rp){
			this->pos  = rp.pos;
			this->mass = rp.mass;
			this->smth = rp.smth;
			this->id   = rp.id;
		}
		PS::F32vec getPos() const{
			return this->pos;
		}
		PS::F32 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F32vec& pos){
			this->pos = pos;
		}
	};
	class Drvt{
		public:
		PS::F32vec pos;
		PS::F32vec vel;
		PS::F32    smth;
		PS::F32    dens;
		void copyFromFP(const RealPtcl& rp){
			pos  = rp.pos;
			vel  = rp.vel;
			dens = rp.dens;
			smth = rp.smth;
		}
		PS::F32vec getPos() const{
			return this->pos;
		}
		PS::F32 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F32vec& pos){
			this->pos = pos;
		}
	};
	class Hydr{
		public:
		PS::F32vec pos;
		PS::F32vec vel;
		PS::F32    smth;
		PS::F32    dens;
		PS::F32    pres;
		PS::F32    snds;
		PS::F32    grad_smth;
		PS::F32    Bal;
		PS::S64    id;///DEBUG
		void copyFromFP(const RealPtcl& rp){
			this->pos   = rp.pos;
			this->vel   = rp.vel;
			this->smth  = rp.smth;
			this->dens  = rp.dens;
			this->pres  = rp.pres;
			this->snds  = rp.snds;
			this->grad_smth = rp.grad_smth;
			this->Bal   = rp.Bal;
			this->id    = rp.id;///DEBUG
		}
		PS::F32vec getPos() const{
			return this->pos;
		}
		PS::F32 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F32vec& pos){
			this->pos = pos;
		}
	};
	class Grav{
		public:
		PS::F32vec pos;
		PS::F32 eps2;
		PS::S64 id;
		PS::F32vec getPos() const{
			return this->pos;
		}
		PS::F32 getEps2(void) const{
			return eps2;
		}
		void copyFromFP(const RealPtcl& rp){
			pos  = rp.pos;
			id   = rp.id;
			eps2 = rp.smth * rp.smth * 1.0e-4;
		}
	};
}

namespace EPJ{
	class Dens{
	public:
		PS::F32    mass;
		PS::F32vec pos;
		PS::F32    smth;
		void copyFromFP(const RealPtcl& rp){
			this->mass = rp.mass;
			this->pos  = rp.pos;
			this->smth = rp.smth;
		}
		PS::F32vec getPos() const{
			return this->pos;
		}
		void setPos(const PS::F32vec& pos){
			this->pos = pos;
		}
		PS::F32 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Drvt{
		public:
		PS::F32    mass;
		PS::F32vec pos;
		PS::F32vec vel;
		PS::F32    smth;
		void copyFromFP(const RealPtcl& rp){
			this->mass = rp.mass;
			this->pos  = rp.pos;
			this->vel  = rp.vel;
			this->smth = rp.smth;
		}
		PS::F32vec getPos() const{
			return this->pos;
		}
		void setPos(const PS::F32vec& pos){
			this->pos = pos;
		}
		PS::F32 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Hydr{
		public:
		PS::F32vec pos;
		PS::F32vec vel;
		PS::F32    dens;
		PS::F32    mass;
		PS::F32    smth;
		PS::F32    pres;
		PS::F32    grad_smth;
		PS::F32    snds;
		PS::F32    Bal;
		PS::S64    id;///DEBUG
		void copyFromFP(const RealPtcl& rp){
			this->pos   = rp.pos;
			this->vel   = rp.vel;
			this->dens  = rp.dens;
			this->pres  = rp.pres;
			this->smth  = rp.smth;
			this->mass  = rp.mass;
			this->snds  = rp.snds;
			this->grad_smth = rp.grad_smth;
			this->Bal = rp.Bal;
			this->id    = rp.id;///DEBUG
		}
		PS::F32vec getPos() const{
			return this->pos;
		}
		PS::F32 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F32vec& pos){
			this->pos = pos;
		}
	};
	class Grav{
		public:
		PS::F32vec pos;
		PS::F32    mass;
		PS::S64    id;
		PS::F32vec getPos() const{
			return this->pos;
		}
		PS::F32 getCharge(void) const{
			return this->mass;
		}
		void copyFromFP(const RealPtcl& rp){
			this->mass = rp.mass;
			this->pos  = rp.pos;
			this->id   = rp.id;
		}
	};
}

template <class Ptcl> class Problem{
	Problem(){
	}
	public:
};

class CalcDensity{
	kernel_t kernel;
	public:
	void operator () (const EPI::Dens* const ep_i, const PS::S32 Nip, const EPJ::Dens* const ep_j, const PS::S32 Njp, RESULT::Dens* const dens){
		#ifdef __AVX__
			#ifdef ENABLE_KNL
				#warning KNL MODE
				const int excess = Nip % 16;
				for(PS::S32 i = 0 ; i + 15 < Nip ; i += 16){
					__m512 dens = _mm512_set1_ps(0.0f);
					__m512i idx = _mm512_set_epi32(60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0);
					__m512  irx = _mm512_i32gather_ps(idx, &ep_i[i].pos.x, 4);
					__m512  iry = _mm512_i32gather_ps(idx, &ep_i[i].pos.y, 4);
					__m512  irz = _mm512_i32gather_ps(idx, &ep_i[i].pos.z, 4);
					__m512  ih  = _mm512_i32gather_ps(idx, &ep_i[i].smth , 4);
					for(PS::S32 j = 0 ; j < Njp ; ++ j){
						__m512 drx = _mm512_sub_ps(irx, _mm512_set1_ps(ep_j[j].pos.x));
						__m512 dry = _mm512_sub_ps(iry, _mm512_set1_ps(ep_j[j].pos.y));
						__m512 drz = _mm512_sub_ps(irz, _mm512_set1_ps(ep_j[j].pos.z));
						__m512 jm  = _mm512_set1_ps(ep_j[j].mass);
						//__m512 r   = _mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(drx, drx), _mm512_mul_ps(dry, dry)), _mm512_mul_ps(drz, drz));
						__m512 r   = _mm512_sqrt_ps(_mm512_fmadd_ps(drx, drx, _mm512_fmadd_ps(dry, dry, _mm512_mul_ps(drz, drz))));
						dens = _mm512_fmadd_ps(kernel.W(r, ih), jm, dens);
					}
				}
			#else
				const int excess = Nip % 8;
				for(PS::S32 i = 0 ; i + 7 < Nip ; i += 8){
					AVXf idens = 0.0f;
					const AVXf irx = v8sf{ep_i[i+0].pos.x, ep_i[i+1].pos.x, ep_i[i+2].pos.x, ep_i[i+3].pos.x, ep_i[i+4].pos.x, ep_i[i+5].pos.x, ep_i[i+6].pos.x, ep_i[i+7].pos.x};
					const AVXf iry = v8sf{ep_i[i+0].pos.y, ep_i[i+1].pos.y, ep_i[i+2].pos.y, ep_i[i+3].pos.y, ep_i[i+4].pos.y, ep_i[i+5].pos.y, ep_i[i+6].pos.y, ep_i[i+7].pos.y};
					const AVXf irz = v8sf{ep_i[i+0].pos.z, ep_i[i+1].pos.z, ep_i[i+2].pos.z, ep_i[i+3].pos.z, ep_i[i+4].pos.z, ep_i[i+5].pos.z, ep_i[i+6].pos.z, ep_i[i+7].pos.z};
					const AVXf ih  = v8sf{ep_i[i+0].smth, ep_i[i+1].smth, ep_i[i+2].smth, ep_i[i+3].smth, ep_i[i+4].smth, ep_i[i+5].smth, ep_i[i+6].smth, ep_i[i+7].smth};
					for(PS::S32 j = 0 ; j < Njp ; ++ j){
						const EPJ::Dens& jth = ep_j[j];
						const AVXf jrx = jth.pos.x;
						const AVXf jry = jth.pos.y;
						const AVXf jrz = jth.pos.z;
						const AVXf jm  = jth.mass;
						const AVXf drx = jrx - irx;
						const AVXf dry = jry - iry;
						const AVXf drz = jrz - irz;
						const AVXf r = (drx * drx + dry * dry + drz * drz).sqrt();
						idens += jm * kernel.W(r, ih);
					}
					for(int s = 0 ; s < 8 ; ++ s) dens[i + s].dens = idens[s];
				}
				for(PS::S32 i = Nip - excess ; i < Nip ; ++ i){
					const EPI::Dens& ith = ep_i[i];
					for(PS::S32 j = 0 ; j < Njp ; ++ j){
						const EPJ::Dens& jth = ep_j[j];
						const PS::F32vec dr = jth.pos - ith.pos;
						dens[i].dens += jth.mass * kernel.W(dr, ith.smth);
					}
				}
			#endif
		#else
			#pragma simd
			for(PS::S32 i = 0 ; i < Nip ; ++ i){
				const EPI::Dens& ith = ep_i[i];
				for(PS::S32 j = 0 ; j < Njp ; ++ j){
					const EPJ::Dens& jth = ep_j[j];
					const PS::F32vec dr = jth.pos - ith.pos;
					dens[i].dens += jth.mass * kernel.W(dr, ith.smth);
				}
			}
		#endif
	}
};

class CalcHydroForce{
	const kernel_t kernel;
	public:
	void operator () (const EPI::Hydr* const ep_i, const PS::S32 Nip, const EPJ::Hydr* const ep_j, const PS::S32 Njp, RESULT::Hydr* const hydro){
		#ifdef __AVX__
			#ifdef ENABLE_KNL
				#warning KNL MODE
				__m512 ax;
			#else
				const int excess = Nip % 8;
				for(PS::S32 i = 0 ; i + 7 < Nip ; i += 8){
					AVXf v_sig_max = 0.0f;
					AVXf accx = 0.0f;
					AVXf accy = 0.0f;
					AVXf accz = 0.0f;
					AVXf eng_dot = 0.0f;

					const AVXf irx = v8sf{ep_i[i+0].pos.x, ep_i[i+1].pos.x, ep_i[i+2].pos.x, ep_i[i+3].pos.x, ep_i[i+4].pos.x, ep_i[i+5].pos.x, ep_i[i+6].pos.x, ep_i[i+7].pos.x};
					const AVXf iry = v8sf{ep_i[i+0].pos.y, ep_i[i+1].pos.y, ep_i[i+2].pos.y, ep_i[i+3].pos.y, ep_i[i+4].pos.y, ep_i[i+5].pos.y, ep_i[i+6].pos.y, ep_i[i+7].pos.y};
					const AVXf irz = v8sf{ep_i[i+0].pos.z, ep_i[i+1].pos.z, ep_i[i+2].pos.z, ep_i[i+3].pos.z, ep_i[i+4].pos.z, ep_i[i+5].pos.z, ep_i[i+6].pos.z, ep_i[i+7].pos.z};
					const AVXf ivx = v8sf{ep_i[i+0].vel.x, ep_i[i+1].vel.x, ep_i[i+2].vel.x, ep_i[i+3].vel.x, ep_i[i+4].vel.x, ep_i[i+5].vel.x, ep_i[i+6].vel.x, ep_i[i+7].vel.x};
					const AVXf ivy = v8sf{ep_i[i+0].vel.y, ep_i[i+1].vel.y, ep_i[i+2].vel.y, ep_i[i+3].vel.y, ep_i[i+4].vel.y, ep_i[i+5].vel.y, ep_i[i+6].vel.y, ep_i[i+7].vel.y};
					const AVXf ivz = v8sf{ep_i[i+0].vel.z, ep_i[i+1].vel.z, ep_i[i+2].vel.z, ep_i[i+3].vel.z, ep_i[i+4].vel.z, ep_i[i+5].vel.z, ep_i[i+6].vel.z, ep_i[i+7].vel.z};
					const AVXf ih  = v8sf{ep_i[i+0].smth, ep_i[i+1].smth, ep_i[i+2].smth, ep_i[i+3].smth, ep_i[i+4].smth, ep_i[i+5].smth, ep_i[i+6].smth, ep_i[i+7].smth};
					const AVXf idens = v8sf{ep_i[i+0].dens, ep_i[i+1].dens, ep_i[i+2].dens, ep_i[i+3].dens, ep_i[i+4].dens, ep_i[i+5].dens, ep_i[i+6].dens, ep_i[i+7].dens};
					const AVXf ipres = v8sf{ep_i[i+0].pres, ep_i[i+1].pres, ep_i[i+2].pres, ep_i[i+3].pres, ep_i[i+4].pres, ep_i[i+5].pres, ep_i[i+6].pres, ep_i[i+7].pres};
					const AVXf isnds = v8sf{ep_i[i+0].snds, ep_i[i+1].snds, ep_i[i+2].snds, ep_i[i+3].snds, ep_i[i+4].snds, ep_i[i+5].snds, ep_i[i+6].snds, ep_i[i+7].snds};
					const AVXf iBal  = v8sf{ep_i[i+0].Bal, ep_i[i+1].Bal, ep_i[i+2].Bal, ep_i[i+3].Bal, ep_i[i+4].Bal, ep_i[i+5].Bal, ep_i[i+6].Bal, ep_i[i+7].Bal};
					for(PS::S32 j = 0; j < Njp ; ++ j){
						const EPJ::Hydr& jth = ep_j[j];
						const AVXf jrx = jth.pos.x;
						const AVXf jry = jth.pos.y;
						const AVXf jrz = jth.pos.z;
						const AVXf jvx = jth.vel.x;
						const AVXf jvy = jth.vel.y;
						const AVXf jvz = jth.vel.z;
						const AVXf jm  = jth.mass;
						const AVXf jh  = jth.smth;
						const AVXf jdens = jth.dens;
						const AVXf jpres = jth.pres;
						const AVXf jsnds = jth.snds;
						const AVXf jBal  = jth.Bal;

						const AVXf drx = jrx - irx;
						const AVXf dry = jry - iry;
						const AVXf drz = jrz - irz;
						const AVXf dvx = jvx - ivx;
						const AVXf dvy = jvy - ivy;
						const AVXf dvz = jvz - ivz;
						const AVXf r    = (drx * drx + dry * dry + drz * drz).sqrt();
						const AVXf drdv = (drx * dvx + dry * dvy + drz * dvz);
						const AVXf w_ij = __builtin_ia32_minps256((drdv / r).v, v8sf{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f});
						const AVXf v_sig = isnds + jsnds - 3.0f * w_ij;
						v_sig_max = __builtin_ia32_maxps256(v_sig_max.v, v_sig.v);
						const AVXf AV = - v_sig * w_ij / (idens + jdens) * 0.5f * (iBal + jBal);
						const AVXf gradW = 0.5f * (kernel.abs_gradW(r, ih) + kernel.abs_gradW(r, jh));

						accx += - jm * (ipres / (idens * idens) + jpres / (jdens * jdens) + AV) * gradW * drx / r;
						accy += - jm * (ipres / (idens * idens) + jpres / (jdens * jdens) + AV) * gradW * dry / r;
						accz += - jm * (ipres / (idens * idens) + jpres / (jdens * jdens) + AV) * gradW * drz / r;
						eng_dot += jm * (ipres / (idens * idens) + 0.5f * AV) * drdv * gradW / r;

					}
					for(int s = 0 ; s < 8 ; ++ s){
						hydro[i + s].acc.x = accx[s];
						hydro[i + s].acc.y = accy[s];
						hydro[i + s].acc.z = accz[s];
						hydro[i + s].eng_dot = eng_dot[s];
						hydro[i + s].dt = PARAM::C_CFL * 2.0f * ih[s] / v_sig_max[s];
					}
				}
				for(PS::S32 i = Nip - excess ; i < Nip ; ++ i){
					PS::F32 v_sig_max = 0.0f;
					const EPI::Hydr& ith = ep_i[i];
					for(PS::S32 j = 0; j < Njp ; ++ j){
						const EPJ::Hydr& jth = ep_j[j];
						const PS::F32vec dr = ith.pos - jth.pos;
						const PS::F32vec dv = ith.vel - jth.vel;
						const PS::F32 w_ij = (dv * dr < 0.f) ? dv * dr / sqrt(dr * dr) : 0.f;
						const PS::F32 v_sig = ith.snds + jth.snds - 3.0f * w_ij;
						v_sig_max = std::max(v_sig_max, v_sig);
						const PS::F32 AV = - 0.5f * v_sig * w_ij / (0.5f * (ith.dens + jth.dens)) * 0.5f * (ith.Bal + jth.Bal);
						const PS::F32vec gradW = 0.5f * (kernel.gradW(dr, ith.smth) + kernel.gradW(dr, jth.smth));
						hydro[i].acc     -= jth.mass * (ith.pres / (ith.dens * ith.dens) + jth.pres / (jth.dens * jth.dens) + AV) * gradW;
						hydro[i].eng_dot += jth.mass * (ith.pres / (ith.dens * ith.dens) + 0.5f * AV) * dv * gradW;
					}
					hydro[i].dt = PARAM::C_CFL * 2.0f * ith.smth / v_sig_max;
				}
			#endif
		#else
		#pragma simd
			for(PS::S32 i = 0; i < Nip ; ++ i){
				PS::F32 v_sig_max = 0.0f;
				const EPI::Hydr& ith = ep_i[i];
				for(PS::S32 j = 0; j < Njp ; ++ j){
					const EPJ::Hydr& jth = ep_j[j];
					const PS::F32vec dr = ith.pos - jth.pos;
					const PS::F32vec dv = ith.vel - jth.vel;
					const PS::F32 w_ij = (dv * dr < 0.f) ? dv * dr / sqrt(dr * dr) : 0.f;
					const PS::F32 v_sig = ith.snds + jth.snds - 3.0f * w_ij;
					v_sig_max = std::max(v_sig_max, v_sig);
					const PS::F32 AV = - 0.5f * v_sig * w_ij / (0.5f * (ith.dens + jth.dens)) * 0.5f * (ith.Bal + jth.Bal);
					const PS::F32vec gradW = 0.5f * (kernel.gradW(dr, ith.smth) + kernel.gradW(dr, jth.smth));
					hydro[i].acc     -= jth.mass * (ith.pres / (ith.dens * ith.dens) + jth.pres / (jth.dens * jth.dens) + AV) * gradW;
					hydro[i].eng_dot += jth.mass * (ith.pres / (ith.dens * ith.dens) + 0.5f * AV) * dv * gradW;
				}
				hydro[i].dt = PARAM::C_CFL * 2.0f * ith.smth / v_sig_max;
			}
		#endif
	}
};

class CalcDerivative{
	kernel_t kernel;
	public:
	void operator () (const EPI::Drvt* ep_i, const PS::S32 Nip, const EPJ::Drvt* ep_j, const PS::S32 Njp, RESULT::Drvt* const drvt){
		#ifdef __AVX__
			#ifdef ENABLE_KNL
				#warning KNL MODE
				__m512 ax;
			#else
				const int excess = Nip % 8;
				for(PS::S32 i = 0 ; i + 7 < Nip ; i += 8){
					AVXf idiv_v  = 0.0f;
					AVXf irot_vx = 0.0f;
					AVXf irot_vy = 0.0f;
					AVXf irot_vz = 0.0f;
					const AVXf irx = v8sf{ep_i[i+0].pos.x, ep_i[i+1].pos.x, ep_i[i+2].pos.x, ep_i[i+3].pos.x, ep_i[i+4].pos.x, ep_i[i+5].pos.x, ep_i[i+6].pos.x, ep_i[i+7].pos.x};
					const AVXf iry = v8sf{ep_i[i+0].pos.y, ep_i[i+1].pos.y, ep_i[i+2].pos.y, ep_i[i+3].pos.y, ep_i[i+4].pos.y, ep_i[i+5].pos.y, ep_i[i+6].pos.y, ep_i[i+7].pos.y};
					const AVXf irz = v8sf{ep_i[i+0].pos.z, ep_i[i+1].pos.z, ep_i[i+2].pos.z, ep_i[i+3].pos.z, ep_i[i+4].pos.z, ep_i[i+5].pos.z, ep_i[i+6].pos.z, ep_i[i+7].pos.z};
					const AVXf ivx = v8sf{ep_i[i+0].vel.x, ep_i[i+1].vel.x, ep_i[i+2].vel.x, ep_i[i+3].vel.x, ep_i[i+4].vel.x, ep_i[i+5].vel.x, ep_i[i+6].vel.x, ep_i[i+7].vel.x};
					const AVXf ivy = v8sf{ep_i[i+0].vel.y, ep_i[i+1].vel.y, ep_i[i+2].vel.y, ep_i[i+3].vel.y, ep_i[i+4].vel.y, ep_i[i+5].vel.y, ep_i[i+6].vel.y, ep_i[i+7].vel.y};
					const AVXf ivz = v8sf{ep_i[i+0].vel.z, ep_i[i+1].vel.z, ep_i[i+2].vel.z, ep_i[i+3].vel.z, ep_i[i+4].vel.z, ep_i[i+5].vel.z, ep_i[i+6].vel.z, ep_i[i+7].vel.z};
					const AVXf ih  = v8sf{ep_i[i+0].smth, ep_i[i+1].smth, ep_i[i+2].smth, ep_i[i+3].smth, ep_i[i+4].smth, ep_i[i+5].smth, ep_i[i+6].smth, ep_i[i+7].smth};
					const AVXf idens = v8sf{ep_i[i+0].dens, ep_i[i+1].dens, ep_i[i+2].dens, ep_i[i+3].dens, ep_i[i+4].dens, ep_i[i+5].dens, ep_i[i+6].dens, ep_i[i+7].dens};
					for(PS::S32 j = 0 ; j < Njp ; ++ j){
						const EPJ::Drvt& jth = ep_j[j];
						const AVXf jrx = jth.pos.x;
						const AVXf jry = jth.pos.y;
						const AVXf jrz = jth.pos.z;
						const AVXf jvx = jth.vel.x;
						const AVXf jvy = jth.vel.y;
						const AVXf jvz = jth.vel.z;
						const AVXf jm  = jth.mass;
						const AVXf drx = jrx - irx;
						const AVXf dry = jry - iry;
						const AVXf drz = jrz - irz;
						const AVXf dvx = jvx - ivx;
						const AVXf dvy = jvy - ivy;
						const AVXf dvz = jvz - ivz;
						const AVXf r    = (drx * drx + dry * dry + drz * drz).sqrt();
						const AVXf drdv = (drx * dvx + dry * dvy + drz * dvz);
						idiv_v  += - jm * drdv / r * kernel.abs_gradW(r, ih);
						irot_vx += - jm * (dry * dvz - drz * dvy) / r * kernel.abs_gradW(r, ih);
						irot_vy += - jm * (drz * dvx - drx * dvz) / r * kernel.abs_gradW(r, ih);
						irot_vz += - jm * (drx * dvy - dry * dvx) / r * kernel.abs_gradW(r, ih);
					}
					idiv_v = idiv_v / idens;
					irot_vx = irot_vx / idens;
					irot_vy = irot_vy / idens;
					irot_vz = irot_vz / idens;
					for(int s = 0 ; s < 8 ; ++ s){
						drvt[i + s].div_v = idiv_v[s];
						drvt[i + s].rot_v.x = irot_vx[s];
						drvt[i + s].rot_v.y = irot_vy[s];
						drvt[i + s].rot_v.z = irot_vz[s];
					}
				}
				for(PS::S32 i = Nip - excess ; i < Nip ; ++ i){
					const EPI::Drvt& ith = ep_i[i];
					for(PS::S32 j = 0; j < Njp ; ++ j){
						const EPJ::Drvt& jth = ep_j[j];
						const PS::F32vec dr = ith.pos - jth.pos;
						const PS::F32vec dv = ith.vel - jth.vel;
						drvt[i].div_v += - jth.mass * dv * kernel.gradW(dr, ith.smth);
						drvt[i].rot_v += - jth.mass * dv ^ kernel.gradW(dr, ith.smth);
					}
					drvt[i].div_v /= ith.dens;
					drvt[i].rot_v /= ith.dens;
				}
			#endif
		#else
			#pragma simd
			for(PS::S32 i = 0; i < Nip ; ++ i){
				const EPI::Drvt& ith = ep_i[i];
				for(PS::S32 j = 0; j < Njp ; ++ j){
					const EPJ::Drvt& jth = ep_j[j];
					const PS::F32vec dr = ith.pos - jth.pos;
					const PS::F32vec dv = ith.vel - jth.vel;
					drvt[i].div_v += - jth.mass * dv * kernel.gradW(dr, ith.smth);
					drvt[i].rot_v += - jth.mass * dv ^ kernel.gradW(dr, ith.smth);
				}
				drvt[i].div_v /= ith.dens;
				drvt[i].rot_v /= ith.dens;
			}
		#endif
	}
};

#ifdef ENABLE_PHANTOM_GRAPE_X86
template <class TPtclJ> class CalcGravityForce{
	static const PS::F32 G = 6.67e-11;
	public:
	void operator () (const EPI::Grav * iptcl, const PS::S32 ni, const TPtclJ * jptcl, const PS::S32 nj, RESULT::Grav * force){
		const PS::S32 nipipe = ni;
		const PS::S32 njpipe = nj;
		PS::F32 (*xi)[3] = (PS::F32 (*)[3])malloc(sizeof(PS::F32) * nipipe * PS::DIMENSION);
		PS::F32 (*ai)[3] = (PS::F32 (*)[3])malloc(sizeof(PS::F32) * nipipe * PS::DIMENSION);
		PS::F32  *pi     = (PS::F32  *    )malloc(sizeof(PS::F32) * nipipe);
		PS::F32 (*xj)[3] = (PS::F32 (*)[3])malloc(sizeof(PS::F32) * njpipe * PS::DIMENSION);
		PS::F32  *mj     = (PS::F32  *    )malloc(sizeof(PS::F32) * njpipe);
		for(PS::S32 i = 0; i < ni; i++) {
			xi[i][0] = iptcl[i].getPos()[0];
			xi[i][1] = iptcl[i].getPos()[1];
			xi[i][2] = iptcl[i].getPos()[2];
			ai[i][0] = 0.0;
			ai[i][1] = 0.0;
			ai[i][2] = 0.0;
			pi[i]    = 0.0;
		}
		for(PS::S32 j = 0; j < nj; j++) {
			xj[j][0] = jptcl[j].getPos()[0];
			xj[j][1] = jptcl[j].getPos()[1];
			xj[j][2] = jptcl[j].getPos()[2];
			mj[j]    = jptcl[j].getCharge();
			xj[j][0] = jptcl[j].pos[0];
			xj[j][1] = jptcl[j].pos[1];
			xj[j][2] = jptcl[j].pos[2];
			mj[j]    = jptcl[j].mass;
		}
		PS::S32 devid = PS::Comm::getThreadNum();
		g5_set_xmjMC(devid, 0, nj, xj, mj);
		g5_set_nMC(devid, nj);
		g5_calculate_force_on_xMC(devid, xi, ai, pi, ni);
		for(PS::S32 i = 0; i < ni; i++) {
			force[i].acc[0] += ai[i][0] * G;
			force[i].acc[1] += ai[i][1] * G;
			force[i].acc[2] += ai[i][2] * G;
			force[i].pot    -= pi[i] * G;
		}
		free(xi);
		free(ai);
		free(pi);
		free(xj);
		free(mj);
	}
};
#else
template <class TPtclJ> class CalcGravityForce{
	static const PS::F32 G = 6.67e-11;
	public:
	void operator () (const EPI::Grav* const __restrict ep_i, const PS::S32 Nip, const TPtclJ* const __restrict ep_j, const PS::S32 Njp, RESULT::Grav* const grav){
		#pragma simd
		for(PS::S32 i = 0; i < Nip ; ++ i){
			const EPI::Grav& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const TPtclJ& jth = ep_j[j];
				const PS::F32vec dr = ith.pos - jth.pos;
				const PS::F32 dr2 = dr * dr;
				const PS::F32 dr_inv = 1.0 / sqrt(dr2 + ith.getEps2());
				const PS::F32 m_dr3_inv = jth.mass * math::pow3(dr_inv);
				grav[i].acc -= m_dr3_inv * dr;
				grav[i].pot -= jth.mass * dr_inv;
			}
			grav[i].acc *= G;
			grav[i].pot *= G;
		}
	}
};
#endif


