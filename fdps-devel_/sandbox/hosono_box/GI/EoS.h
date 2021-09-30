#pragma once

namespace EoS{
	enum{
		Monoatomic,
		Diatomic  ,
		Granite   ,
		Iron      ,
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
		virtual type Pressure  (const type dens, const type eng) const = 0;
		virtual type SoundSpeed(const type dens, const type eng) const = 0;
	};
	//////////////////
	//EoSs
	//////////////////
	template <typename type> class IdealGas : public EoS_t<type>{
		const type hcr;//heat capacity ratio;
		public:
		IdealGas(const type _hcr) : hcr(_hcr){
		}
		inline type Pressure(const type dens, const type eng) const{
			return math::max((hcr - 1.0) * dens * eng, 0.0);
		}
		inline type SoundSpeed(const type dens, const type eng) const{
			return sqrt(math::abs(hcr * (hcr - 1.0) * eng));
		}
		inline type HeatCapacityRatio() const{
			return hcr;
		}
	};
	template <typename type> class Tillotson : public EoS_t<type>{
		type rho0, a, b, A, B, u0, alpha, beta, uiv, ucv;
		inline type P_co(const type dens, const type eng) const{
			const type eta = dens / rho0;
			const type mu  = eta - 1.0;
			return (a + b / (eng / u0 / eta / eta + 1.0)) * dens * eng + A * mu + B * mu * mu;
		}
		inline type P_ex(const type dens, const type eng) const{
			const type eta = dens / rho0;
			const type mu  = eta - 1.0;
			return a * dens * eng + (b * dens * eng / (eng / u0 / eta / eta + 1.0) + A * mu * exp(- alpha * (1.0 / eta - 1.0))) * exp(- beta * (1.0 / eta - 1.0) * (1.0 / eta - 1.0));
		}
		inline type dPdrho(const type rho, const type u) const{
			const type drho = 0.0001;
			return (Pressure(rho + drho, u) - Pressure(rho - drho, u)) / (2.0 * drho);
		}
		inline type dPdu(const type rho, const type u) const{
			const type du = 0.0001;
			return (Pressure(rho, u + du) - Pressure(rho, u - du)) / (2.0 * du);
		}
		public:
		Tillotson(const type a_rho0, const type a_u0, const type a_uiv, const type a_ucv, const type a_A, const type a_B, const type a_a, const type a_b, const type a_alpha, const type a_beta){
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
		inline type Pressure(const type dens, const type eng) const{
			const double p_min = 1.0e+7;
			#ifdef INITIAL_CONDITION
			#warning this is for IC!
			//if IC
			if(dens >= rho0 || eng < uiv){
				return std::max(P_co(dens, eng), p_min);
			}else if(dens < rho0 && eng > ucv){
				return std::max(P_ex(dens, eng), p_min);
			}else{
				return std::max(((eng - uiv) * P_ex(dens, eng) + (ucv - eng) * P_co(dens, eng)) / (ucv - uiv), p_min);
			}
			#else
			//if run
			double p;
			if(dens >= rho0 || eng < uiv){
				p = P_co(dens, eng);
				//if(dens <= 0.9 * rho0) return 1.0e-16;
			}else if(dens < rho0 && eng > ucv){
				p = P_ex(dens, eng);
			}else{
				p = ((eng - uiv) * P_ex(dens, eng) + (ucv - eng) * P_co(dens, eng)) / (ucv - uiv);
			}
			return (p > p_min) ? p : p_min;
			#endif
		}
		inline type SoundSpeed(const type dens, const type eng) const{
			#if 0
			const type drho = 0.0001;
			if(dens >= rho0 || eng < uiv){
				if(dens - drho <= 0.9 * rho0) return sqrt(A / dens);
			}
			#endif
			return sqrt(std::max(Pressure(dens, eng) / (dens * dens) * dPdu(dens, eng) + dPdrho(dens, eng), 1.0e-16));
		}
	};
}

