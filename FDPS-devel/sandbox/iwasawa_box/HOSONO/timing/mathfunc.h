#pragma once

namespace math{
	const PS::F64 pi = atan(1.0) * 4.0;
	const PS::F64 NaN = + 0.0 / 0.0;
	const PS::F64 VERY_LARGE_VALUE = 1.0e+30;
	template <typename type> inline type abs(type arg){
		return (arg > 0) ? arg : - arg;
	}
	template <typename type> inline type plus(type arg){
		return (arg > 0) ? arg : 0;
	}
	template <typename type> inline type sign(type arg){
		return (arg > 0) ? 1.0 : - 1.0;
	}
	///////////
	template <typename type> inline type power(type base, int exp){
		return pow(base, exp);
	}
	template <typename type> inline type root(type base, int idx){
		if(idx == 1){
			return base;
		}else if(idx == 2){
			return sqrt(base);
		}else if(idx == 3){
			return cbrt(base);
		}else{
			return NaN;
		}
	}
	template <typename type> inline type power(type base, long long int exp){
		return pow(base, exp);
	}
	template <typename type> inline type root(type base, long long int idx){
		if(idx == 1){
			return base;
		}else if(idx == 2){
			return sqrt(base);
		}else if(idx == 3){
			return cbrt(base);
		}else{
			return NaN;
		}
	}
	///////////
	template <typename type> inline type pow2(type arg){
		return arg * arg;
	}
	template <typename type> inline type pow3(type arg){
		return arg * arg * arg;
	}
	template <typename type> inline type pow4(type arg){
		const type arg2 = arg * arg;
		return arg2 * arg2;
	}
	template <typename type> inline type pow5(type arg){
		const type arg2 = arg * arg;
		return arg2 * arg2 * arg;
	}
	template <typename type> inline type pow6(type arg){
		const type arg3 = arg * arg * arg;
		return arg3 * arg3;
	}
	template <typename type> inline type pow7(type arg){
		const type arg3 = arg * arg * arg;
		return arg3 * arg3 * arg;
	}
	template <typename type> inline type pow8(type arg){
		const type arg2 = arg * arg;
		const type arg4 = arg2 * arg2;
		return arg4 * arg4;
	}
	template <typename type> inline type min(type a, type b){
		return (a < b) ? a : b;
	}
	template <typename type> inline type max(type a, type b){
		return (a > b) ? a : b;
	}
}
