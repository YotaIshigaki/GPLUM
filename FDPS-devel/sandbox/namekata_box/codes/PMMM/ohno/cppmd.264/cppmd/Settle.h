#ifndef SETTLE_H
#define SETTLE_H

#include <vector>
#include <cmath>
#include "SpaceVector.h"
#include "ParticleInfo.h"
#include "CovalentBondInfo.h"

class Settle {
public:
	double dt, dtInv;
	double massO, massH, massOH, massH2O;
	double massOinv, massHinv, massH2Oinv;
	double massOinvHinv, massHinvHinv;
	double angleO;
	double lengthOH, lengthHH;
	double length2OH, length2HH;

	double model_HX, model_HY, model_OY; //r_c, r_b, r_a in original paper

	double tolerance;
	double tolerance2PosOH, tolerance2PosHH;
	double toleranceVelOH,	toleranceVelHH;
	
	Settle();

	virtual ~Settle() {}

  template<class PA>
  void constrain_position(PA& particle_array,
                          const PosVelArray& particle_array_prev,
                          const WaterList& waterlist, 
                          const TypeRange& type_range,
                          const double dt);
  template<class PA>
  void constrain_position(PA& particle_array,
                          const PosVelArray& particle_array_prev,
                          const WaterList& waterlist, 
                          const TypeRange& type_range,
                          const double dt,
			  ForceArray& water_settle_force);
  template<class PA>
  void constrain_velocity(PA& particle_array,
			  const WaterList& waterlist, 
			  const TypeRange& type_range,
			  const double dt,
			  ForceArray& water_settle_force);
  template<class PA>
  void constrain_velocity(PA& particle_array,
			  const WaterList& waterlist, 
			  const TypeRange& type_range,
			  const double dt,
			  double &virial);
  void constrain_velocity(CombinedParticleArray& particle_array,
                          const WaterList& waterlist, 
                          const TypeRange& type_range);
  void constrain_position(ParticleArray& particlearray,
                          const ParticleArray& particle_array_prev,
                          const WaterList& water_list, 
                          const TypeRange& type_range,
                          double dt);
  void constrain_velocity(ParticleArray& particlearray, 
                          const WaterList& water_list, 
                          const TypeRange& type_range);

  void rattle_position(ParticleArray& particlearray,
                       const ParticleArray& particle_array_prev,
                       const WaterList& water_list, 
                       const TypeRange& type_range,
                       double dt);
  void rattle_velocity(ParticleArray& particlearray, 
                       const WaterList& water_list, 
                       const TypeRange& type_range,
                       double dt);

  template<class PA,class PVA>
  void rattle_position(PA& particlearray,
                       const PVA& particle_array_prev,
                       const ShakeList& shake_list, 
                       const CovalentBondParameterList* param_list, 
                       std::vector<CovalentBondInfo::BondList>& bondlistarray,
                       double dt,
                       int shake_max_iterate,
                       double shake_tolerance,
		       ForceArray& rattle_force);
  template<class PA,class PVA>
  void rattle_position(PA& particlearray,
                       const PVA& particle_array_prev,
                       const ShakeList& shake_list, 
                       const CovalentBondParameterList* param_list, 
                       std::vector<CovalentBondInfo::BondList>& bondlistarray,
                       double dt,
                       int shake_max_iterate,
                       double shake_tolerance);
  template<class PA>
  void rattle_velocity(PA& particlearray, 
                       const ShakeList& shake_list, 
                       const CovalentBondParameterList* param_list, 
                       std::vector<CovalentBondInfo::BondList>& bondlistarray,
                       double dt,
                       int shake_max_iterate,
                       double shake_tolerance);
  template<class PA>
  void rattle_velocity(PA& particle_array, 
		       const ShakeList& shakelist, 
		       const CovalentBondParameterList* param_list, 
		       std::vector<CovalentBondInfo::BondList>& bondlistarray,
		       double dt,
		       int shake_max_iterate, double shake_tolerance,
		       ForceArray& rattle_force);
  template<class PA>
  void rattle_velocity(PA& particle_array, 
		       const ShakeList& shakelist, 
		       const CovalentBondParameterList* param_list, 
		       std::vector<CovalentBondInfo::BondList>& bondlistarray,
		       double dt,
		       int shake_max_iterate, double shake_tolerance,
		       double& rattle_virial);
};

#endif
