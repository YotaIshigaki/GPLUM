#pragma once

namespace EoS{
	//////////////////
	//abstract class
	//////////////////
	template <typename type> class EoS_t{
		public:
		EoS_t(){
			return ;
		}
		/* virtual */ ~EoS_t(){
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
			return (hcr - 1.0) * dens * eng;
		}
		inline type SoundSpeed(const type dens, const type eng) const{
			return sqrt(hcr * (hcr - 1.0) * eng);
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
			const type p_min = 1.0e+7;
			if(dens >= rho0 || eng < uiv){
				return std::max(P_co(dens, eng), p_min);
			}else if(dens < rho0 && eng > ucv){
				return std::max(P_ex(dens, eng), p_min);
			}else{
				return std::max(((eng - uiv) * P_ex(dens, eng) + (ucv - eng) * P_co(dens, eng)) / (ucv - uiv), p_min);
			}
		}
		inline type SoundSpeed(const type dens, const type eng) const{
			return sqrt(std::max(Pressure(dens, eng) / (dens * dens) * dPdu(dens, eng) + dPdrho(dens, eng), 0.0) + 1.0e-16);
		}
	};
}

static const EoS::IdealGas<PS::F64>  Monoatomic(5./3.);
static const EoS::IdealGas<PS::F64>  Diatomic  (1.4);
static const EoS::Tillotson<PS::F64> Granite   (2680.0, 16.0e+6, 3.5e+6, 18.00e+6,  18.0e+9,  18.0e+9, 0.5, 1.3, 5.0, 5.0);
static const EoS::Tillotson<PS::F64> Iron      (7800.0,  9.5e+6, 2.4e+6 , 8.67e+6, 128.0e+9, 105.0e+9, 0.5, 1.5, 5.0, 5.0);
static const EoS::Tillotson<PS::F64> Water     ( 998.0,  7.0e+6, 0.419e+6,  2.69e+6,   2.18e+9,  /*13.25e+9*/2.18e+9, 0.7, 0.15, 10.0, 5.0);
static const EoS::Tillotson<PS::F64> WetTuff   (1970.0, 11.0e+6, 3.200e+6, 16.00e+6,  10.00e+9,   6.00e+9, 0.5, 1.30,  5.0, 5.0);

