#pragma once


//Wendland C6
struct kernel_t{
	kernel_t(){
	}
	//W
	PS::F64 W(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
#if 0
		r_value = math::pow8(math::plus(1.0 - s)) * (1.0 + 8.0 * s + 25.0 * s * s + 32.0 * s * s * s);
		r_value *= (1365./64.) / (H * H * H * math::pi);
#else
		r_value = math::pow8(math::plus(1.0 - s)) * (1.0 + s*(8.0 + s*(25.0 + s*(32.0))));
		r_value *= (1365./(64.*math::pi)) / (H * H * H);
#endif
		return r_value;
	}
	//gradW
	PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
#if 0
		r_value = math::pow7(math::plus(1.0 - s)) * (math::plus(1.0 - s) * (8.0 + 50.0 * s + 96.0 * s * s) - 8.0 * (1.0 + 8.0 * s + 25.0 * s * s + 32.0 * s * s * s));
		r_value *= (1365./64.) / (H * H * H * math::pi);
#else
		const PS::F64 plus = math::plus(1.0 - s);
		r_value = math::pow7(plus) * (
			(plus * (8.0 + s*(50.0 + s*(96.0))))
		   	- (8.0 + s*(64.0 + s*(200.0 + s*(256.0)))));
		r_value *= (1365./(64.*math::pi)) / (H * H * H);
#endif
		return dr * r_value / (sqrt(dr * dr) * H  + 1.0e-6 * h);
	}
	static PS::F64 supportRadius(){
		return 2.5;
	}

	static PS::F64 W(const PS::F64 s, const PS::F64 Hinv) {
		PS::F64 r_value = math::pow8(math::plus(1.0 - s));
		r_value *= (1.0 + s*(8.0 + s*(25.0 + s*(32.0))));
		r_value *= (1365./(64.*math::pi)) * (Hinv * Hinv * Hinv);

		return r_value;
	}
	static PS::F64vec gradW(const PS::F64 s, const PS::F64 Hinv, const PS::F64vec dr, const PS::F64 rinv) { 
		return gradWs(s, Hinv, rinv) * dr;
	}
	// multiply dr later
	static PS::F64 gradWs(const PS::F64 s, const PS::F64 Hinv, const PS::F64 rinv) { 
		const PS::F64 plus = math::plus(1.0 - s);
		PS::F64 r_value = math::pow7(plus);
		r_value *= (plus * (8.0 + s*(50.0 + s*(96.0))))
		           - (8.0 + s*(64.0 + s*(200.0 + s*(256.0))));
		r_value *= (1365./(64.*math::pi)) * (Hinv * Hinv * Hinv * Hinv);
		return (r_value * rinv);
	}
	static PS::F64 safe_norm2(const PS::F64vec dr, const PS::F64 h) {
		//const PS::F64 eps2 = (1.e-16)*(h*h);
		const PS::F64 eps2 = (1.e-6)*(h*h);
		return (eps2 + (dr*dr));
	}
};

