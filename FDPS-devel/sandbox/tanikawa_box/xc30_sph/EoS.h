#pragma once

namespace EoS{
	enum{
		Monoatomic,
		Diatomic,
		Granite,
	};
	//////////////////
	//abstract class
	//////////////////
	template <typename type> class EoS_t{
		public:
		EoS_t(){
			return ;
		}
		~EoS_t(){
			return ;
		}
		virtual type Pressure  (type dens, type eng) const = 0;
		virtual type SoundSpeed(type dens, type eng) const = 0;
	};
	//////////////////
	//EoSs
	//////////////////
	template <typename type> class IdealGas : public EoS_t<type>{
		const type hcr;//heat capacity ratio;
		public:
		IdealGas(type _hcr) : hcr(_hcr){
		}
		inline type Pressure(type dens, type eng) const{
			return math::max((hcr - 1.0) * dens * eng, 0.0);
		}
		inline type SoundSpeed(type dens, type eng) const{
			return sqrt(math::abs(hcr * (hcr - 1.0) * eng));
		}
		inline type HeatCapacityRatio() const{
			return hcr;
		}
	};
	template <typename type> class Tillotson : public EoS_t<type>{
		type rho0, a, b, A, B, u0, alpha, beta, uiv, ucv;
		inline type P_co(type dens, type eng) const{
			type eta = dens / rho0;
			type mu  = eta - 1.0;
			return (a + b / (eng / u0 / eta / eta + 1.0)) * dens * eng + A * mu + B * mu * mu;
		}
		inline type P_ex(type dens, type eng) const{
			type eta = dens / rho0;
			type mu  = eta - 1.0;
			return a * dens * eng + (b * dens * eng / (eng / u0 / eta / eta + 1.0) + A * mu * exp(- alpha * (1.0 / eta - 1.0))) * exp(- beta * (1.0 / eta - 1.0) * (1.0 / eta - 1.0));
		}
		inline type dPdrho(type rho, type u) const{
			type drho = 0.0001;
			return (Pressure(rho + drho, u) - Pressure(rho - drho, u)) / (2.0 * drho);
		}
		inline type dPdu(type rho, type u) const{
			type du = 0.0001;
			return (Pressure(rho, u + du) - Pressure(rho, u - du)) / (2.0 * du);
		}
		public:
		Tillotson(type a_rho0, type a_u0, type a_uiv, type a_ucv, type a_A, type a_B, type a_a, type a_b, type a_alpha, type a_beta){
			//in MKS unit...
			rho0  = a_rho0; // kg/m^3
			u0    = a_u0;   // J/kg
			uiv   = a_uiv;  // J/kg
			ucv   = a_ucv;  // J/kg
			A     = a_A;    // Pa
			B     = a_B;    // Pa
			a     = a_a;    // dimension-less
			b     = a_b;    //
			alpha = a_alpha;//
			beta  = a_beta; //
		}
		inline type Pressure(type dens, type eng) const{
			const double p_min = 1.0e+8;
			if(dens >= rho0 || eng < uiv){
				return std::max(P_co(dens, eng), p_min);
			}else if(dens < rho0 && eng > ucv){
				return std::max(P_ex(dens, eng), p_min);
			}else{
				return std::max(((eng - uiv) * P_ex(dens, eng) + (ucv - eng) * P_co(dens, eng)) / (ucv - uiv), p_min);
			}
		}
		inline type SoundSpeed(type dens, type eng) const{
			//std::cout << sqrt(std::max(Pressure(dens, eng) / (dens * dens) * dPdu(dens, eng) + dPdrho(dens, eng), 1.0e-16)) << std::endl;
			return sqrt(std::max(Pressure(dens, eng) / (dens * dens) * dPdu(dens, eng) + dPdrho(dens, eng), 1.0e-16));
		}
	};
}
