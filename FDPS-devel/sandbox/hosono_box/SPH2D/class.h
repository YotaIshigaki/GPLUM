#pragma once

enum TYPE{
	HYDRO,
	FREEZE,
};

struct system_t{
	PS::F64 dt, time;
	PS::S64 step;
	system_t() : step(0), time(0.0){
	}
};

class FileHeader{
public:
	int Nbody;
	double time;
	int readAscii(FILE* fp){
		fscanf(fp, "%lf\n", &time);
		fscanf(fp, "%d\n", &Nbody);
		return Nbody;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%d\n", Nbody);
	}
};

namespace STD{
	namespace RESULT{
		//Density summation
		class Dens{
			public:
			PS::F64 dens;
			PS::F64 smth;
			void clear(){
				dens = smth = 0;
			}
		};
		//for Balsara switch
		class Drvt{
			public:
			PS::F64 div_v;
			PS::F64vec rot_v;
			PS::F64 grad_smth;
			void clear(){
				div_v = 0.0;
				grad_smth = 0.0;
				rot_v = 0.0;
			}
		};
		//Hydro force
		class Hydro{
			public:
			PS::F64vec acc;
			PS::F64 eng_dot;
			PS::F64 dt;
			void clear(){
				acc = 0;
				eng_dot = 0;
				dt = 1.0e+30;
			}
		};
		//Self gravity
		class Grav{
			public:
			PS::F64vec acc;
			PS::F64    pot;
			PS::F64    dt;
			void clear(){
				acc = 0.0;
				pot = 0.0;
				dt = 1.0e+30;
			}
		};
	}

	class RealPtcl{
		public:
		PS::F64 mass;
		PS::F64vec pos, vel, acc;
		PS::F64 dens;//DENSity
		PS::F64 eng; //ENerGy
		PS::F64 pres;//PRESsure
		PS::F64 smth;//SMooTHing length
		PS::F64 snds; //SouND Speed
		PS::F64 div_v;
		PS::F64vec rot_v;
		PS::F64 Bal; //Balsala switch
		PS::F64 AVa; //Time dependent AV_alpha
		PS::F64 grad_smth;

		PS::F64 eng_dot;
		PS::F64vec vel_half;
		PS::F64 eng_half;
		PS::F64 dt;
		PS::S64 id, tag;
		const EoS::EoS_t<PS::F64>* EoS;

		PS::F64vec grav;
		PS::F64    pot;

		TYPE type;
		//Constructor
		RealPtcl(){
			type = HYDRO;
			AVa = 1.0;
			Bal = 1.0;
		}
		//Copy functions
		void copyFromForce(const RESULT::Dens& dens){
			this->dens = dens.dens;
			this->smth = dens.smth;
		}
		void copyFromForce(const RESULT::Drvt& drvt){
			this->div_v = drvt.div_v;
			this->rot_v = drvt.rot_v;
			if(PARAM::FLAG_B95 == true){
				this->Bal = std::abs(drvt.div_v) / (std::abs(drvt.div_v) + sqrt(drvt.rot_v * drvt.rot_v) + 1.0e-4 * this->snds / this->smth); //Balsala switch
			}else{
				this->Bal = 1.0;
			}
			this->grad_smth = drvt.grad_smth;
		}
		void copyFromForce(const RESULT::Hydro& force){
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
		PS::F64 getCharge() const{
			return this->mass;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		PS::F64 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
		void writeAscii(FILE* fp) const{
			#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
			fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",  id,  tag,  mass,  pos.x,  pos.y,  0.0  ,  vel.x,  vel.y,  0.0  ,  dens,  eng,  pres, rot_v.x);
			#else
			fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",  id,  tag,  mass,  pos.x,  pos.y,  pos.z,  vel.x,  vel.y,  vel.z,  dens,  eng,  pres, pot);
			#endif
		}
		void readAscii(FILE* fp){
			#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
			double dummy1, dummy2;
			fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &tag, &mass, &pos.x, &pos.y, &dummy1, &vel.x, &vel.y, &dummy2, &dens, &eng, &pres, &pot);
			#else
			fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &tag, &mass, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z, &dens, &eng, &pres, &pot);
			#endif
		}
		void setPressure(const EoS::EoS_t<PS::F64>* const _EoS){
			EoS = _EoS;
		}
		void initialize(){
			smth = PARAM::SMTH * pow(mass / dens, 1.0/(PS::F64)(PARAM::Dim));
			grad_smth = 1.0;
		}
		void initialKick(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			vel_half = vel + 0.5 * dt_glb * (acc + grav);
			eng_half = eng + 0.5 * dt_glb * eng_dot;
		}
		void fullDrift(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			pos += dt_glb * vel_half;
			if(PARAM::FLAG_R00 == true){
				AVa += ((2.0 - AVa) * std::max(- div_v, 0.0) - (AVa - 0.1) / (smth / snds)) * dt_glb;
			}
		}
		void predict(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			vel += dt_glb * (acc + grav);
			eng += dt_glb * eng_dot;
		}
		void finalKick(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			vel = vel_half + 0.5 * dt_glb * (acc + grav);
			eng = eng_half + 0.5 * dt_glb * eng_dot;
		}
	};

	namespace EPI{
		class Dens{
			public:
			PS::F64vec pos;
			PS::F64    mass;
			PS::F64    smth;
			PS::S64    id;
			void copyFromFP(const RealPtcl& rp){
				this->pos  = rp.pos;
				this->mass = rp.mass;
				this->smth = rp.smth;
				this->id   = rp.id;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
		};

		class Drvt{
			public:
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    smth;
			PS::F64    dens;
			void copyFromFP(const RealPtcl& rp){
				pos  = rp.pos;
				vel  = rp.vel;
				dens = rp.dens;
				smth = rp.smth;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
		};

		class Hydro{
			public:
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    smth;
			PS::F64    dens;
			PS::F64    pres;
			PS::F64    snds;
			PS::F64    grad_smth;
			PS::F64    Bal;
			PS::F64    AVa;
			PS::S64    id;///DEBUG
			void copyFromFP(const RealPtcl& rp){
				this->pos       = rp.pos;
				this->vel       = rp.vel;
				this->smth      = rp.smth;
				this->dens      = rp.dens;
				this->pres      = rp.pres;
				this->snds      = rp.snds;
				this->grad_smth = rp.grad_smth;
				this->Bal       = rp.Bal;
				this->AVa       = rp.AVa;
				this->id        = rp.id;///DEBUG
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
		};
		class Grav{
			public:
			PS::F64vec pos;
			PS::F64    eps2;
			PS::S64    id;
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getEps2(void) const{
				//return (1.0e-4 * 6400.0e+3) * (1.0e-4 * 6400.0e+3);//GI Unit
				return eps2;
			}
			void copyFromFP(const RealPtcl& rp){
				pos = rp.pos;
				id  = rp.id;
				eps2 = 1.0e-2 * rp.smth * rp.smth;
				//eps2 = rp.smth * rp.smth;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
		};
	}

	namespace EPJ{
		class Dens{
		public:
			PS::F64    mass;
			PS::F64vec pos;
			PS::F64    smth;
			void copyFromFP(const RealPtcl& rp){
				this->mass = rp.mass;
				this->pos  = rp.pos;
				this->smth = rp.smth;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
		};
		class Drvt{
			public:
			PS::F64    mass;
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    smth;
			void copyFromFP(const RealPtcl& rp){
				this->mass = rp.mass;
				this->pos  = rp.pos;
				this->vel  = rp.vel;
				this->smth = rp.smth;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
		};

		class Hydro{
			public:
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    dens;
			PS::F64    mass;
			PS::F64    smth;
			PS::F64    pres;
			PS::F64    grad_smth;
			PS::F64    snds;
			PS::F64    Bal;
			PS::F64    AVa;
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
				this->AVa   = rp.AVa;
				this->id    = rp.id;///DEBUG
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
		};

		class Grav{
			public:
			PS::F64vec pos;
			PS::F64    mass;
			PS::S64    id;
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getCharge(void) const{
				return this->mass;
			}
			void copyFromFP(const RealPtcl& rp){
				this->mass = rp.mass;
				this->pos  = rp.pos;
				this->id   = rp.id;
			}
		};
	}
}

namespace DI{
	namespace RESULT{
		//Density summation
		class Dens{
			public:
			PS::F64 dens_smth;
			PS::F64 pres_smth;
			PS::F64 smth;
			void clear(){
				dens_smth = pres_smth = smth = 0;
			}
		};
		//for Balsara switch
		class Drvt{
			public:
			PS::F64 div_v;
			PS::F64vec rot_v;
			PS::F64 grad_smth;
			void clear(){
				div_v = 0.0;
				grad_smth = 0.0;
				rot_v = 0.0;
			}
		};
		//Hydro force
		class Hydro{
			public:
			PS::F64vec acc;
			PS::F64 eng_dot;
			PS::F64 dt;
			void clear(){
				acc = 0;
				eng_dot = 0;
				dt = 1.0e+30;
			}
		};
		//Self gravity
		class Grav{
			public:
			PS::F64vec acc;
			PS::F64    pot;
			PS::F64    dt;
			void clear(){
				acc = 0.0;
				pot = 0.0;
				dt = 1.0e+30;
			}
		};
	}

	class RealPtcl{
		public:
		PS::F64 mass;
		PS::F64vec pos, vel, acc;
		PS::F64 dens, dens_smth;//DENSity
		PS::F64 eng; //ENerGy
		PS::F64 pres_smth;//PRESsure
		PS::F64 smth;//SMooTHing length
		PS::F64 snds; //SouND Speed
		PS::F64 div_v;
		PS::F64vec rot_v;
		PS::F64 Bal; //Balsala switch
		PS::F64 AVa; //Time dependent AV_alpha
		PS::F64 grad_smth;
		PS::F64 pV;//Y in Hosono+ (2013)

		PS::F64 eng_dot;
		PS::F64vec vel_half;
		PS::F64 eng_half;
		PS::F64 dt;

		PS::S64 id, tag;
		const EoS::EoS_t<PS::F64>* EoS;

		PS::F64vec grav;
		PS::F64    pot;

		TYPE type;
		//Constructor
		RealPtcl(){
			//type = HYDRO;
			AVa = 1.0;
			Bal = 1.0;
		}
		//Copy functions
		void copyFromForce(const RESULT::Dens& dens){
			this->dens_smth = dens.dens_smth;
			this->pres_smth = dens.pres_smth;
			this->smth = dens.smth;
		}
		void copyFromForce(const RESULT::Drvt& drvt){
			this->div_v = drvt.div_v;
			this->rot_v = drvt.rot_v;
			if(PARAM::FLAG_B95 == true){
				this->Bal = std::abs(drvt.div_v) / (std::abs(drvt.div_v) + sqrt(drvt.rot_v * drvt.rot_v) + 1.0e-4 * this->snds / this->smth); //Balsala switch
			}else{
				this->Bal = 1.0;
			}
			this->grad_smth = drvt.grad_smth;
		}
		void copyFromForce(const RESULT::Hydro& force){
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
		PS::F64 getCharge() const{
			return this->mass;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		PS::F64 getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
		void writeAscii(FILE* fp) const{
			#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
			fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",  id,  tag,  mass,  pos.x,  pos.y,  0.0  ,  vel.x,  vel.y,  0.0  ,  dens,  eng, pow(pres_smth, 1.0 / PARAM::DISPH_POWER), pot);
			#else
			fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",  id,  tag,  mass,  pos.x,  pos.y,  pos.z,  vel.x,  vel.y,  vel.z,  dens,  eng, pow(pres_smth, 1.0 / PARAM::DISPH_POWER), pot);
			#endif
		}
		void readAscii(FILE* fp){
			#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
			double dummy1, dummy2;
			fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &tag, &mass, &pos.x, &pos.y, &dummy1, &vel.x, &vel.y, &dummy2, &dens, &eng, &pres_smth, &pot);
			#else
			fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &tag, &mass, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z, &dens, &eng, &pres_smth, &pot);
			#endif
			pres_smth = pow(pres_smth, PARAM::DISPH_POWER);
		}
		void setPressure(const EoS::EoS_t<PS::F64>* const _EoS){
			EoS = _EoS;
		}
		void initialize(){
			smth = PARAM::SMTH * pow(mass / dens, 1.0/(PS::F64)(PARAM::Dim));
			pV   = pow(EoS->Pressure(dens, eng), PARAM::DISPH_POWER) * mass / dens;
			grad_smth = 1.0;
		}
		void initialKick(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			vel_half = vel + 0.5 * dt_glb * (acc + grav);
			eng_half = eng + 0.5 * dt_glb * eng_dot;
		}
		void fullDrift(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			pos += dt_glb * vel_half;
			dens += - dens * dt_glb * div_v;
			if(PARAM::FLAG_R00 == true){
				AVa += ((2.0 - AVa) * std::max(- div_v, 0.0) - (AVa - 0.1) / (smth / snds)) * dt_glb;
			}
		}
		void predict(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			vel += dt_glb * (acc + grav);
			eng += dt_glb * eng_dot;
			pV   = pow(EoS->Pressure(dens, eng), PARAM::DISPH_POWER) * mass / dens;
		}
		void finalKick(const PS::F64 dt_glb){
			//if(type == FREEZE) return;
			vel = vel_half + 0.5 * dt_glb * (acc + grav);
			eng = eng_half + 0.5 * dt_glb * eng_dot;
		}
	};
	namespace EPI{
		class Dens{
			public:
			PS::F64vec pos;
			PS::F64    mass;
			PS::F64    smth;
			PS::S64    id;
			void copyFromFP(const RealPtcl& rp){
				this->pos  = rp.pos;
				this->mass = rp.mass;
				this->smth = rp.smth;
				this->id   = rp.id;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
		};
		class Drvt{
			public:
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    smth;
			PS::F64    dens_smth;
			PS::F64    pres_smth;
			//DEBUG
			PS::S64    id;
			PS::F64    snds;
			void copyFromFP(const RealPtcl& rp){
				pos  = rp.pos;
				vel  = rp.vel;
				smth = rp.smth;
				dens_smth = rp.dens_smth;
				pres_smth = rp.pres_smth;
				id   = rp.id;
				snds = rp.snds;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
		};

		class Hydro{
			public:
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    smth;
			PS::F64    dens_smth;
			PS::F64    pres_smth;
			PS::F64    mass;
			PS::F64    pV;
			PS::F64    snds;
			PS::F64    grad_smth;
			PS::F64    Bal;
			PS::F64    AVa;
			PS::S64    id;///DEBUG
			void copyFromFP(const RealPtcl& rp){
				this->pos   = rp.pos;
				this->vel   = rp.vel;
				this->smth  = rp.smth;
				this->dens_smth = rp.dens_smth;
				this->pres_smth = rp.pres_smth;
				this->mass  = rp.mass;
				this->pV    = rp.pV;
				this->snds  = rp.snds;
				this->grad_smth = rp.grad_smth;
				this->Bal   = rp.Bal;
				this->AVa   = rp.AVa;
				this->id    = rp.id;///DEBUG
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
		};

		class Grav{
			public:
			PS::F64vec pos;
			PS::F64    eps2;
			PS::S64    id;
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getEps2(void) const{
				//return (1.0e-4 * 6400.0e+3) * (1.0e-4 * 6400.0e+3);//GI Unit
				return eps2;
			}
			void copyFromFP(const RealPtcl& rp){
				pos = rp.pos;
				id  = rp.id;
				//eps2 = 1.0e-2 * rp.smth * rp.smth;
				eps2 = rp.smth * rp.smth;
			}
		};
	}

	namespace EPJ{
		class Dens{
		public:
			PS::F64    mass;
			PS::F64vec pos;
			PS::F64    dens;
			PS::F64    eng;
			PS::F64    smth;
			PS::F64    pV;
			void copyFromFP(const RealPtcl& rp){
				this->mass = rp.mass;
				this->pos  = rp.pos;
				this->dens = rp.dens;
				this->eng  = rp.eng;
				this->smth = rp.smth;
				this->pV   = rp.pV  ;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
		};
		class Drvt{
			public:
			PS::F64    pV ;
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    smth;
			PS::F64    mass;
			void copyFromFP(const RealPtcl& rp){
				this->pV   = rp.pV;
				this->pos  = rp.pos;
				this->vel  = rp.vel;
				this->smth = rp.smth;
				this->mass = rp.mass;
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
		};

		class Hydro{
			public:
			PS::F64vec pos;
			PS::F64vec vel;
			PS::F64    dens_smth;
			PS::F64    mass;
			PS::F64    smth;
			PS::F64    pres_smth;
			PS::F64    grad_smth;
			PS::F64    pV;
			PS::F64    snds;
			PS::F64    Bal;
			PS::F64    AVa;
			PS::S64    id;///DEBUG
			void copyFromFP(const RealPtcl& rp){
				this->pos       = rp.pos;
				this->vel       = rp.vel;
				this->dens_smth = rp.dens_smth;
				this->pres_smth = rp.pres_smth;
				this->smth      = rp.smth;
				this->mass      = rp.mass;
				this->snds      = rp.snds;
				this->grad_smth = rp.grad_smth;
				this->pV        = rp.pV;
				this->Bal       = rp.Bal;
				this->AVa       = rp.AVa;
				this->id        = rp.id;///DEBUG
			}
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getRSearch() const{
				return kernel_t::supportRadius() * this->smth;
			}
			void setPos(const PS::F64vec& pos){
				this->pos = pos;
			}
		};

		class Grav{
			public:
			PS::F64vec pos;
			PS::F64    mass;
			PS::S64    id;
			PS::F64vec getPos() const{
				return this->pos;
			}
			PS::F64 getCharge(void) const{
				return this->mass;
			}
			void copyFromFP(const RealPtcl& rp){
				this->mass = rp.mass;
				this->pos  = rp.pos;
				this->id   = rp.id;
			}
		};
	}
}

template <class Ptcl> class Problem{
	Problem(){
	}
	public:
	static void setupIC(PS::ParticleSystem<Ptcl>&, system_t&, PS::DomainInfo&){
	}
	static void setEoS(PS::ParticleSystem<Ptcl>&){
	}
	static void setDomain(PS::DomainInfo&){
	}
	static void addExternalForce(PS::ParticleSystem<Ptcl>&, system_t&){
		std::cout << "No Ext. Force" << std::endl;
	}
	static void postTimestepProcess(PS::ParticleSystem<Ptcl>&, system_t&){
	}
};

