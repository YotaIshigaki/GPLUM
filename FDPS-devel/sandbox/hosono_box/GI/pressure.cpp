#include "header.h"

static const EoS::IdealGas<PS::F64>  Monoatomic(5./3.);
static const EoS::IdealGas<PS::F64>  Diatomic  (1.4);
static const EoS::Tillotson<PS::F64> Granite(2680.0, 16.0e+6 , 3.5e+6 , 18.00e+6,  18.0e+9,  18.0e+9, 0.5, 1.3, 5.0, 5.0);
static const EoS::Tillotson<PS::F64> Iron   (7800.0,  9.5e+6 , 2.4e+6 ,  8.67e+6, 128.0e+9, 105.0e+9, 0.5, 1.5, 5.0, 5.0);

void CalcPressure(PS::ParticleSystem<RealPtcl>& sph_system){
	std::vector<const EoS::EoS_t<PS::F64>*> EoSs;
	EoSs.push_back(&Monoatomic);
	EoSs.push_back(&Diatomic);
	EoSs.push_back(&Granite);
	EoSs.push_back(&Iron);
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].setPressure(EoSs[sph_system[i].mat]);
	}
}

