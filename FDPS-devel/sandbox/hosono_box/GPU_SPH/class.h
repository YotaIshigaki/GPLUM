#pragma once

struct system_t{
	float dt, time, end_time;
	PS::S32 step;
	system_t(){
		step = 0;
		time = 0.0;
	}
};

class FileHeader{
public:
	int Nbody;
	float time;
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

//Wendland C6
struct kernel_t{
	kernel_t(){
	}
	//W
	float W(const PS::F64vec dr, const PS::F64 h) const{
		const float H = supportRadius() * h;
		const float s = sqrt(dr * dr) / H;
		float r_value;
		r_value = (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))) * math::pow8(math::plus(1.0 - s));
		r_value *= (1365./64.) / (H * H * H * math::pi);
		return r_value;
	}
	//gradW
	PS::F32vec gradW(const PS::F64vec dr, const PS::F64 h) const{
		const float H = supportRadius() * h;
		const float s = sqrt(dr * dr) / H;
		float r_value;
		r_value = math::pow7(math::plus(1.0 - s)) * (math::plus(1.0 - s) * (8.0 + s * (50.0 + s * (96.0))) - 8.0 * (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))));
		r_value *= (1365./64.) / (H * H * H * math::pi);
		return dr * r_value / (sqrt(dr * dr) * H  + 1.0e-6 * h);
	}
	static float supportRadius(){
		return 3.5;
	}
};


namespace RESULT{
	//Density summation
	class Dens{
		public:
		float dens;
		float smth;
		void clear(){
			dens = smth = 0;
		}
	};
	//for Balsara switch
	class Drvt{
		public:
		float div_v;
		PS::F32vec rot_v;
		float grad_smth;
		void clear(){
			div_v = 0.0;
			grad_smth = 0.0;
			rot_v = 0.0;
		}
	};
	//Hydro force
	class Hydro{
		public:
		PS::F32vec acc;
		float eng_dot;
		float dt;
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
		float pot;
		float dt;
		void clear(){
			acc = 0.0;
			pot = 0.0;
			dt = 1.0e+30;
		}
	};
}

class RealPtcl{
	public:
	float mass;
	PS::F32vec pos, vel, acc;
	float dens;//DENSity
	float eng; //ENerGy
	float pres;//PRESsure
	float smth;//SMooTHing length
	float snds; //SouND Speed
	float div_v;
	PS::F32vec rot_v;
	float Bal; //Balsala switch
	float grad_smth;

	float eng_dot;
	PS::F32vec vel_half;
	float eng_half;
	float dt;
	PS::S64 id;
	PS::S64 mat;//material
	PS::S32 tag;

	//Copy functions
	void copyFromForce(const RESULT::Dens& dens){
		this->dens = dens.dens;
		this->smth = dens.smth;
	}
	void copyFromForce(const RESULT::Drvt& drvt){
		this->div_v = drvt.div_v;
		this->rot_v = drvt.rot_v;
		this->Bal = abs(drvt.div_v) / (abs(drvt.div_v) + sqrt(drvt.rot_v * drvt.rot_v) + 1.0e-4 * this->snds / this->smth); //Balsala switch
		this->grad_smth = drvt.grad_smth;
	}
	void copyFromForce(const RESULT::Hydro& force){
		this->acc     = force.acc;
		this->eng_dot = force.eng_dot;
		this->dt      = force.dt;
	}
	void copyFromForce(const RESULT::Grav& force){
		//not copy dt
	}
	//Give necessary values to FDPS
	float getCharge() const{
		return this->mass;
	}
	PS::F32vec getPos() const{
		return this->pos;
	}
	float getRSearch() const{
		return kernel_t::supportRadius() * this->smth;
	}
	void setPos(const PS::F32vec& pos){
		this->pos = pos;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",  id,  mat,  mass,  pos.x,  pos.y,  pos.z,  vel.x,  vel.y,  vel.z,  dens,  eng,  pres);
		//fprintf(fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",  id,  mat,  mass,  pos.x,  pos.y,  pos.z,  vel.x,  vel.y,  vel.z,  acc.x,  acc.y,  acc.z);
	}
	void readAscii(FILE* fp){
		fscanf (fp, "%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &id, &mat, &mass, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z, &dens, &eng, &pres);
	}
	/*
	void setPressure(const EoS::EoS_t<float>* const EoS){
		pres = EoS->Pressure(dens, eng);
		snds = EoS->SoundSpeed(dens, eng);
	}
	*/
	void setPressure(void){
		const float hcr = 1.4;
		pres = (hcr - 1.0) * dens * eng;
		pres = sqrt((hcr - 1.0) * hcr * eng);
	}
};

namespace EPI{
	class Dens{
		public:
		PS::F32vec pos;
		float    mass;
		float    smth;
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
		float getRSearch() const{
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
		float    smth;
		float    dens;
		void copyFromFP(const RealPtcl& rp){
			pos  = rp.pos;
			vel  = rp.vel;
			dens = rp.dens;
			smth = rp.smth;
		}
		PS::F32vec getPos() const{
			return this->pos;
		}
		float getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Hydro{
		public:
		PS::F32vec pos;
		PS::F32vec vel;
		float    smth;
		float    dens;
		float    pres;
		float    snds;
		float    grad_smth;
		float    Bal;
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
		float getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Grav{
		public:
		PS::F32vec pos;
		float eps2;
		PS::S64 id;
		PS::F32vec getPos() const{
			return this->pos;
		}
		float getEps2(void) const{
			return 1.0e-4;
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
		float    mass;
		PS::F32vec pos;
		float    smth;
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
		float getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Drvt{
		public:
		float    mass;
		PS::F32vec pos;
		PS::F32vec vel;
		float    smth;
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
		float getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
	};
	class Hydro{
		public:
		PS::F32vec pos;
		PS::F32vec vel;
		float    dens;
#if 0
		float    pres;
		float    smth;
		float    mass;
#else
		float    mass;
		float    smth;
		float    pres;
#endif
		float    grad_smth;
		float    snds;
		float    Bal;
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
		float getRSearch() const{
			return kernel_t::supportRadius() * this->smth;
		}
		void setPos(const PS::F32vec& pos){
			this->pos = pos;
		}
	};
	class Grav{
		public:
		PS::F32vec pos;
		float    mass;
		PS::S64    id;
		PS::F32vec getPos() const{
			return this->pos;
		}
		float getCharge(void) const{
			return this->mass;
		}
		void copyFromFP(const RealPtcl& rp){
			this->mass = rp.mass;
			this->pos  = rp.pos;
			this->id   = rp.id;
		}
	};
}

