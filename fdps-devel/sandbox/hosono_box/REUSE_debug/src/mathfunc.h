#pragma once

namespace math{
	const PS::F64 pi = atan(1.0) * 4.0;
	const PS::F64 NaN = + 0.0 / 0.0;
	const PS::F64 VERY_LARGE_VALUE = 1.0e+30;
	template <typename type> inline type plus(const type arg){
		return (arg > 0) ? arg : 0;
	}
	template <typename type> inline type sign(const type arg){
		return (arg > 0) ? 1.0 : - 1.0;
	}
	template <typename type> inline type pow2(const type arg){
		return arg * arg;
	}
	template <typename type> inline type pow3(const type arg){
		return arg * arg * arg;
	}
	template <typename type> inline type pow4(const type arg){
		const type arg2 = arg * arg;
		return arg2 * arg2;
	}
	template <typename type> inline type pow5(const type arg){
		const type arg2 = arg * arg;
		return arg2 * arg2 * arg;
	}
	template <typename type> inline type pow6(const type arg){
		const type arg3 = arg * arg * arg;
		return arg3 * arg3;
	}
	template <typename type> inline type pow7(const type arg){
		const type arg3 = arg * arg * arg;
		return arg3 * arg3 * arg;
	}
	template <typename type> inline type pow8(const type arg){
		const type arg2 = arg * arg;
		const type arg4 = arg2 * arg2;
		return arg4 * arg4;
	}
}

