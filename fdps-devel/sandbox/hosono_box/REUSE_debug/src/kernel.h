#pragma once

//Wendland C6
struct WendlandC6{
	WendlandC6(){}
	//W
	PS::F64 W(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))) * math::pow8(math::plus(1.0 - s));
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		r_value *= (78./7.) / (H * H * math::pi);
		#else
		r_value *= (1365./64.) / (H * H * H * math::pi);
		#endif
		return r_value;
	}
	//gradW
	PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = math::pow7(math::plus(1.0 - s)) * (math::plus(1.0 - s) * (8.0 + s * (50.0 + s * (96.0))) - 8.0 * (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))));
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		r_value *= (78./7.) / (H * H * math::pi);
		#else
		r_value *= (1365./64.) / (H * H * H * math::pi);
		#endif
		return dr * r_value / (sqrt(dr * dr) * H  + 1.0e-6 * h);
	}
	static PS::F64 supportRadius(){
		return 3.5;
	}
};

typedef WendlandC6 kernel_t;

