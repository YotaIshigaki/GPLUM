#pragma once

struct system_t{
	PS::F64 dt, time, end_time;
	PS::F64 Grav;
	system_t(){
		Grav = 0.0;
	}
};

class FileHeader{
public:
	int Nbody;
	double time;
	int readAscii(FILE* fp){
		fscanf(fp, "%e\n", &time);
		fscanf(fp, "%d\n", &Nbody);
		return Nbody;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%d\n", Nbody);
	}
};

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
		void clear(){
			div_v = 0.0;
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
	PS::F64 Rsearch;
	PS::F64 snds;//SouND Speed
	PS::F64 div_v;
	PS::F64vec rot_v;
	PS::F64 Bal; //Balsala switch

	PS::F64 eng_dot;
	PS::F64vec vel_half;
	PS::F64 eng_half;
	PS::F64 dt;
	PS::S64 id;
	PS::S64 mat;//material
	PS::S32 tag;

	PS::F64vec grav;
	PS::F64    pot;
	//Copy functions
	void copyFromForce(const RESULT::Dens& dens){
		this->dens = dens.dens;
		this->smth = dens.smth;
	}
	void copyFromForce(const RESULT::Drvt& drvt){
		this->div_v = drvt.div_v;
		this->rot_v = drvt.rot_v;
		this->Bal = abs(drvt.div_v) / (abs(drvt.div_v) + sqrt(drvt.rot_v * drvt.rot_v) + 1.0e-4 * this->snds / this->smth); //Balsala switch
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
		#ifdef FAST
		return this->Rsearch;
		#else
		return kernel_t::supportRadius() * this->smth;
		#endif
	}
	void setPos(const PS::F64vec& pos){
		this->pos = pos;
	}
/*
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",  this->id,  this->mat,  this->mass,  this->pos.x,  this->pos.y,  this->pos.z,  this->vel.x,  this->vel.y,  this->vel.z,  this->dens,  this->eng,  this->pres);
	}
*/

	void writeAscii(FILE* fp) const{
		fprintf(fp, "%ld\t%ld\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",  
                this->id,  this->mat,  this->mass,  this->pos.x,  this->pos.y,  this->pos.z,  this->vel.x,  this->vel.y,  this->vel.z,  this->dens,  this->eng,  this->pres,  this->smth, 
                this->acc.x, this->acc.y, this->acc.z, this->grav.x, this->grav.y, this->grav.z, this->eng_dot, this->dt, this->div_v, this->rot_v.x, this->rot_v.y, this->rot_v.z, this->snds);
	}

	void readAscii(FILE* fp){
		fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &this->id, &this->mat, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z, &this->dens, &this->eng, &this->pres);
	}
	void setPressure(const EoS::EoS_t<PS::F64>* const EoS){
		pres = EoS->Pressure(dens, eng);
		snds = EoS->SoundSpeed(dens, eng);
	}
};

namespace EPI{
	class Dens{
		public:
		PS::F64vec pos;
		PS::F64    mass;
		PS::F64    smth;
		PS::F64    Rsearch;
		PS::S64    id;
		void copyFromFP(const RealPtcl& rp){
			this->pos  = rp.pos;
			this->mass = rp.mass;
			this->smth = rp.smth;
			this->Rsearch = rp.Rsearch;
			this->id   = rp.id;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		PS::F64 getRSearch() const{
			#ifdef FAST
			return this->Rsearch;
			#else
			return kernel_t::supportRadius() * this->smth;
			#endif
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
	};
	class Hydro{
		public:
		PS::F64vec pos;
		PS::F64vec vel;
		PS::F64    smth;
		PS::F64    dens;
		PS::F64    pres;
		PS::F64    snds;
		PS::F64    Bal;
		PS::S64    id;///DEBUG
		void copyFromFP(const RealPtcl& rp){
			this->pos   = rp.pos;
			this->vel   = rp.vel;
			this->smth  = rp.smth;
			this->dens  = rp.dens;
			this->pres  = rp.pres;
			this->snds  = rp.snds;
			this->Bal   = rp.Bal;
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
			return (1.0e-2 * 6400.0e+3) * (1.0e-2 * 6400.0e+3);//GI Unit
		}
		void copyFromFP(const RealPtcl& rp){
			pos = rp.pos;
			id  = rp.id;
		}
	};
}

namespace EPJ{
	class Dens{
	public:
		PS::F64    mass;
		PS::F64vec pos;
		PS::F64    smth;
		PS::F64    Rsearch;
		PS::S64    id;
		void copyFromFP(const RealPtcl& rp){
			this->mass = rp.mass;
			this->pos  = rp.pos;
			this->smth = rp.smth;
			this->Rsearch = rp.Rsearch;
			this->id = rp.id;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
		PS::F64 getRSearch() const{
			#ifdef FAST
			return this->Rsearch;
			#else
			return kernel_t::supportRadius() * this->smth;
			#endif
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
		PS::F64    pres;
		PS::F64    smth;
		PS::F64    mass;
		PS::F64    snds;
		PS::F64    Bal;
		PS::S64    id;///DEBUG
		void copyFromFP(const RealPtcl& rp){
			this->pos   = rp.pos;
			this->vel   = rp.vel;
			this->dens  = rp.dens;
			this->pres  = rp.pres;
			this->smth  = rp.smth;
			this->mass  = rp.mass;
			this->snds  = rp.snds;
			this->Bal = rp.Bal;
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

