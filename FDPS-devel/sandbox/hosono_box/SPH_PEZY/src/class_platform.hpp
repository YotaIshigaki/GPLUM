#pragma once
#include "math.h"
#include "param.h"

struct system_t{
	double dt, time, end_time, Grav;
	PS::S32 step;
	system_t(){
		step = 0;
		time = 0.0;
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
		fprintf(fp, "%lf\n", time);
		fprintf(fp, "%d\n", Nbody);
	}
};


namespace RESULT{
	//Density summation
	class Dens{
		public:
		double dens;
		void clear(){
			dens = 0;
		}
	};
	//for Balsara switch
	class Drvt{
		public:
		double div_v;
		PS::F64vec rot_v;
		double grad_smth;
		void clear(){
			div_v = 0.0;
			grad_smth = 0.0;
			rot_v = 0.0;
		}
	};
	//Hydro force
	class Hydr{
		public:
		PS::F64vec acc;
		double eng_dot;
		double dt;
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
		double pot;
		double dt;
		void clear(){
			acc = 0.0;
			pot = 0.0;
			dt = 1.0e+30;
		}
	};
}

class RealPtcl{
	public:
	double mass;
	PS::F64vec pos, vel, acc, grav;
	double dens;//DENSity
	double eng; //ENerGy
	double pres;//PRESsure
	double smth;//SMooTHing length
	double snds; //SouND Speed
	double div_v;
	PS::F64vec rot_v;
	double pot;
	double Bal; //Balsala switch
	double grad_smth;

	double eng_dot;
	PS::F64vec vel_half;
	double eng_half;
	double dt;
	PS::S64 id;
	PS::S64 mat;//material
	const EoS::EoS_t<PS::F64>* EoS;

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
	double getCharge() const{
		return this->mass;
	}
	PS::F64vec getPos() const{
		return this->pos;
	}
	double getRSearch() const{
		return kernel_t::supportRadius() * this->smth;
	}
	void setPos(const PS::F64vec& pos){
		this->pos = pos;
	}
	void writeAscii(FILE* fp) const{
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", id, mat, mass, pos.x, pos.y, 0.0  , vel.x, vel.y, 0.0  , dens, eng, pres, pot);
		#else
		fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", id, mat, mass, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, dens, eng, pres, pot);
		#endif
		//fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", id, mat, mass, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, dens, eng, pres, pot, acc.x, div_v, rot_v.x, Bal);
	}
	void readAscii(FILE* fp){
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &mat, &mass, &pos.x, &pos.y,   NULL, &vel.x, &vel.y,   NULL, &dens, &eng, &pres, &pot);
		#else
		fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &mat, &mass, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z, &dens, &eng, &pres, &pot);
		#endif
		//fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &mat, &mass, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z, &dens, &eng, &pres, &pot, &acc.x, &acc.y, &acc.z, &eng_dot);
	}
	void setPressure(void){
		pres = EoS->Pressure(dens, eng);
		snds = EoS->SoundSpeed(dens, eng);
	}
	void initialize(void){
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		smth = PARAM::SMTH * sqrt(mass / dens);
		#else
		smth = PARAM::SMTH * cbrt(mass / dens);
		#endif
	}
};

namespace EPI{
	class Dens{
		public:
		PS::F64vec pos;
		double    mass;
		double    smth;
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
		double getRSearch() const{
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
		double    smth;
		double    dens;
		void copyFromFP(const RealPtcl& rp){
			pos  = rp.pos;
			vel  = rp.vel;
			dens = rp.dens;
			smth = rp.smth;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		double getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Hydr{
		public:
		PS::F64vec pos;
		PS::F64vec vel;
		double    smth;
		double    dens;
		double    pres;
		double    snds;
		double    grad_smth;
		double    Bal;
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
		PS::F64vec getPos() const{
			return this->pos;
		}
		double getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Grav{
		public:
		PS::F64vec pos;
		double eps2;
		PS::S64 id;
		PS::F64vec getPos() const{
			return this->pos;
		}
		double getEps2(void) const{
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
		double    mass;
		PS::F64vec pos;
		double    smth;
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
		double getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Drvt{
		public:
		double    mass;
		PS::F64vec pos;
		PS::F64vec vel;
		double    smth;
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
		double getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Hydr{
		public:
		PS::F64vec pos;
		PS::F64vec vel;
		double    dens;
		double    mass;
		double    smth;
		double    pres;
		double    grad_smth;
		double    snds;
		double    Bal;
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
		PS::F64vec getPos() const{
			return this->pos;
		}
		double getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
	};
	class Grav{
		public:
		PS::F64vec pos;
		double    mass;
		PS::S64    id;
		PS::F64vec getPos() const{
			return this->pos;
		}
		double getCharge(void) const{
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
	static void setupIC(PS::ParticleSystem<Ptcl>&, system_t&, PS::DomainInfo&){
	}
	static void setDomain(PS::DomainInfo&){
	}
	static void addExternalForce(PS::ParticleSystem<Ptcl>&, system_t&){
		std::cout << "No Ext. Force" << std::endl;
	}
	static void postTimestepProcess(PS::ParticleSystem<Ptcl>&, system_t&){
	}
};



