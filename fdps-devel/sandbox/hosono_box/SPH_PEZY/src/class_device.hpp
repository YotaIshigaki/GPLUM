typedef float real;

namespace Dens{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real mass;
		real smth;
		int id_walk;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real mass;
		real smth;
	};
	struct ForceDev{
		real dens;
	};
}

namespace Drvt{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real dens;
		real smth;
		int id_walk;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real mass;
		real smth;
	};
	struct ForceDev{
		real div_v;
		real rot_vx;
		real rot_vy;
		real rot_vz;
	};
}

namespace Hydr{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real dens;
		real pres;
		real snds;
		real smth;
		real Bal;
		int  id_walk;
		real grad_smth;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real vx;
		real vy;
		real vz;
		real dens;
		real pres;
		real snds;
		real mass;
		real smth;
		real Bal;
		real grad_smth;
	};
	struct ForceDev{
		real ax;
		real ay;
		real az;
		real eng_dot;
		real dt;
	};	
}

namespace Grav{
	struct EpiDev{
		real rx;
		real ry;
		real rz;
		real eps2;
		int id_walk;
	};

	struct EpjDev{
		real rx;
		real ry;
		real rz;
		real mass;
	};
	struct ForceDev{
		real ax;
		real ay;
		real az;
		real pot;
	};	
}


struct kernel_t{
	const real pi = M_PI;
	real plus(const real x) const{
		return (x > 0) ? x : 0;
	}
	real pow8(const real x) const{
		const real x2 = x  * x ;
		const real x4 = x2 * x2;
		return x4 * x4;
	}
	real pow7(const real x) const{
		const real x2 = x  * x ;
		const real x4 = x2 * x2;
		return x4 * x2 * x;
	}
	kernel_t(){
	}
	//W
	real W(const real dr, const real h) const{
		const real H = supportRadius() * h;
		const real s = dr / H;
		real r_value;
		r_value = (1.0f + s * (8.0f + s * (25.0f + s * (32.0f)))) * pow8(plus(1.0f - s));
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		r_value *= (78.f/7.f) / (H * H * pi);
		#else
		r_value *= (1365.f/64.f) / (H * H * H * pi);
		#endif
		return r_value;
	}
	//gradW
	real gradW(const real dr, const real h) const{
		const real H = supportRadius() * h;
		const real s = dr / H;
		real r_value;
		r_value = pow7(plus(1.0f - s)) * (plus(1.0f - s) * (8.0f + s * (50.0f + s * (96.0f))) - 8.0f * (1.0f + s * (8.0f + s * (25.0f + s * (32.0f)))));
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		r_value *= (78.f/7.f) / (H * H * pi);
		#else
		r_value *= (1365.f/64.f) / (H * H * H * pi);
		#endif
		return r_value / (H);
	}
	static real supportRadius(){
		return 3.5f;
	}
};

