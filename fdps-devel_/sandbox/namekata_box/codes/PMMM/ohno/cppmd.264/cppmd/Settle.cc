#include "Settle.h"

#define DEBUG_SHAKE 0

Settle::Settle()
{
  massO = 16.00/12.01;
  massH = 1.008/12.01;
  angleO = 104.520 * M_PI / 180.0;
  lengthOH = 0.9572;

  massOH = massO + massH;
  massH2O = massO + massH + massH;
  massOinv = 1.0 / massO;
  massHinv = 1.0 / massH;
  massH2Oinv = 1.0 / massH2O;

  massOinvHinv = massOinv + massHinv;
  massHinvHinv = massHinv + massHinv;
        
  lengthHH = 2.0 * lengthOH * sin(0.5 * angleO);
  length2OH = lengthOH * lengthOH;
  length2HH = lengthHH * lengthHH;
 
  model_HX = lengthOH * sin(0.5 * angleO);
  model_OY = lengthOH * cos(0.5 * angleO) * 2.0 * massH / massH2O;
  model_HY = lengthOH * cos(0.5 * angleO) * massO / massH2O;

  tolerance = 1e-10;

  tolerance2PosOH = 2.0 * tolerance * length2OH;
  tolerance2PosHH = 2.0 * tolerance * length2HH;
        
  //     toleranceVelOH = tolerance * length2OH / dt;
  //     toleranceVelHH = tolerance * length2HH / dt;

}

template<class PA>
void Settle::constrain_position(PA& particle_array,
                                const PosVelArray& particle_array_prev,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
                                const double dt)
{
  using namespace std;
  double dtInv = 1.0 / dt;
  
  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
		
    // _crd0 : coordinates in world (global) coordinate system
    // _crd1 or none : coordinates in local coordinate system based on the molecule position

    //vector from the center of mass to new (unconstrained) {O,H1,H2} in old coordinate system
    Position com = ((massO * getpos(particle_array,o_index))
                    +(massH * getpos(particle_array,h1_index))
                    +(massH * getpos(particle_array,h2_index))) * massH2Oinv;
    Position o_crd0(getpos(particle_array,o_index) - com);
    Position h1_crd0(getpos(particle_array,h1_index) - com);
    Position h2_crd0(getpos(particle_array,h2_index) - com);

    //vector from the 'com' to prev {O,H1,H2} in old coordinate system
    Position o_prev_crd0(getpos(particle_array_prev,o_index) - com);
    Position h1_prev_crd0(getpos(particle_array_prev,h1_index) - com);
    Position h2_prev_crd0(getpos(particle_array_prev,h2_index) - com);

    // unit vectors of new axes in old coordinate system
    Position ez_new_crd0((h1_prev_crd0 - o_prev_crd0)%(h2_prev_crd0 - o_prev_crd0) );
    ez_new_crd0 /= ez_new_crd0.norm();
    Position ex_new_crd0(o_crd0%ez_new_crd0);
    ex_new_crd0 /= ex_new_crd0.norm();
    Position ey_new_crd0(ez_new_crd0%ex_new_crd0);
    ey_new_crd0 /= ey_new_crd0.norm();
    // unit vectors of old axes in new coordinate system
    Position ex_old_crd1(ex_new_crd0.x, ey_new_crd0.x, ez_new_crd0.x);
    Position ey_old_crd1(ex_new_crd0.y, ey_new_crd0.y, ez_new_crd0.y);
    Position ez_old_crd1(ex_new_crd0.z, ey_new_crd0.z, ez_new_crd0.z);

    // positions in new coordinate system ('crd1' are omitted)
    Position o( ex_new_crd0*o_crd0,	 ey_new_crd0*o_crd0,	ez_new_crd0*o_crd0);
    Position h1(ex_new_crd0*h1_crd0, ey_new_crd0*h1_crd0, ez_new_crd0*h1_crd0);
    Position h2(ex_new_crd0*h2_crd0, ey_new_crd0*h2_crd0, ez_new_crd0*h2_crd0);
    // prev positions in new coordinate system ('crd1' are omitted)
    Position o_prev( ex_new_crd0*o_prev_crd0,	 ey_new_crd0*o_prev_crd0,	 ez_new_crd0*o_prev_crd0);
    Position h1_prev(ex_new_crd0*h1_prev_crd0, ey_new_crd0*h1_prev_crd0, ez_new_crd0*h1_prev_crd0);
    Position h2_prev(ex_new_crd0*h2_prev_crd0, ey_new_crd0*h2_prev_crd0, ez_new_crd0*h2_prev_crd0);

    // first and second rotations (psi and phi)
    double sin_phi = (o.z / model_OY);
    double cos_phi = sqrt(1 - sin_phi * sin_phi);
    double sin_psi = (h1.z - h2.z)/(2 * model_HX * cos_phi);
    double cos_psi = sqrt(1 - sin_psi * sin_psi);
    Position o_2(0.0,
                 model_OY * cos_phi,
                 model_OY * sin_phi);
    Position h1_2(-model_HX * cos_psi,
                  -model_HY * cos_phi - model_HX * sin_psi * sin_phi,
                  -model_HY * sin_phi + model_HX * sin_psi * cos_phi);
    Position h2_2(model_HX * cos_psi,
                  -model_HY * cos_phi + model_HX * sin_psi * sin_phi,
                  -model_HY * sin_phi - model_HX * sin_psi * cos_phi);

    // third rotation (theta)
    double alpha(h1_2.x * (h1_prev.x - h2_prev.x)
                 + h1_2.y * (h1_prev.y - o_prev.y)
                 + h2_2.y * (h2_prev.y - o_prev.y) );
    double beta( h1_2.x	* (h2_prev.y - h1_prev.y)
                 + h1_2.y * (h1_prev.x - o_prev.x)
                 + h2_2.y * (h2_prev.x - o_prev.x) );
    double gamma(h1.y *	(h1_prev.x - o_prev.x)
                 - h1.x * (h1_prev.y - o_prev.y)
                 + h2.y * (h2_prev.x - o_prev.x)
                 - h2.x * (h2_prev.y - o_prev.y) );
    double sin_theta = (alpha * gamma - beta * sqrt(alpha*alpha + beta*beta - gamma*gamma))
        / (alpha*alpha + beta*beta);
    double cos_theta = sqrt(1 - sin_theta*sin_theta);
    
    Position o_3( o_2.x * cos_theta - o_2.y * sin_theta,
                  o_2.x * sin_theta + o_2.y * cos_theta,
                  o_2.z);
    Position h1_3( h1_2.x * cos_theta - h1_2.y * sin_theta,
                   h1_2.x * sin_theta + h1_2.y * cos_theta,
                   h1_2.z);
    Position h2_3( h2_2.x * cos_theta - h2_2.y * sin_theta,
                   h2_2.x * sin_theta + h2_2.y * cos_theta,
                   h2_2.z);

    // new (constrained) positions in old coordinate system
    Position o_3_crd0( ex_old_crd1*o_3,	 ey_old_crd1*o_3,	 ez_old_crd1*o_3);
    Position h1_3_crd0(ex_old_crd1*h1_3, ey_old_crd1*h1_3, ez_old_crd1*h1_3);
    Position h2_3_crd0(ex_old_crd1*h2_3, ey_old_crd1*h2_3, ez_old_crd1*h2_3);

    // update position
    getpos(particle_array,o_index) = o_3_crd0 + com;
    getpos(particle_array,h1_index) = h1_3_crd0 + com;
    getpos(particle_array,h2_index) = h2_3_crd0 + com;

    // update velocity
    getvelocity(particle_array,o_index) += dtInv * (o_3_crd0 - o_crd0);
    getvelocity(particle_array,h1_index) += dtInv * (h1_3_crd0 - h1_crd0);
    getvelocity(particle_array,h2_index) += dtInv * (h2_3_crd0 - h2_crd0);

  }

}
template
void Settle::constrain_position(CombinedParticleArray& particle_array,
                                const PosVelArray& particle_array_prev,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
                                const double dt);
template
void Settle::constrain_position(ParticleArray& particle_array,
                                const PosVelArray& particle_array_prev,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
                                const double dt);
template<class PA>
void Settle::constrain_position(PA& particle_array,
                                const PosVelArray& particle_array_prev,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
                                const double dt,
				ForceArray& water_settle_force)
{
  using namespace std;
  double dtInv = 1.0 / dt;
  double dt2Inv2 = 2.0*dtInv*dtInv;
  
  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
		
    // _crd0 : coordinates in world (global) coordinate system
    // _crd1 or none : coordinates in local coordinate system based on the molecule position

    //vector from the center of mass to new (unconstrained) {O,H1,H2} in old coordinate system
    Position com = ((massO * getpos(particle_array,o_index))
                    +(massH * getpos(particle_array,h1_index))
                    +(massH * getpos(particle_array,h2_index))) * massH2Oinv;
    Position o_crd0(getpos(particle_array,o_index) - com);
    Position h1_crd0(getpos(particle_array,h1_index) - com);
    Position h2_crd0(getpos(particle_array,h2_index) - com);

    //vector from the 'com' to prev {O,H1,H2} in old coordinate system
    Position o_prev_crd0(getpos(particle_array_prev,o_index) - com);
    Position h1_prev_crd0(getpos(particle_array_prev,h1_index) - com);
    Position h2_prev_crd0(getpos(particle_array_prev,h2_index) - com);

    // unit vectors of new axes in old coordinate system
    Position ez_new_crd0((h1_prev_crd0 - o_prev_crd0)%(h2_prev_crd0 - o_prev_crd0) );
    ez_new_crd0 /= ez_new_crd0.norm();
    Position ex_new_crd0(o_crd0%ez_new_crd0);
    ex_new_crd0 /= ex_new_crd0.norm();
    Position ey_new_crd0(ez_new_crd0%ex_new_crd0);
    ey_new_crd0 /= ey_new_crd0.norm();
    // unit vectors of old axes in new coordinate system
    Position ex_old_crd1(ex_new_crd0.x, ey_new_crd0.x, ez_new_crd0.x);
    Position ey_old_crd1(ex_new_crd0.y, ey_new_crd0.y, ez_new_crd0.y);
    Position ez_old_crd1(ex_new_crd0.z, ey_new_crd0.z, ez_new_crd0.z);

    // positions in new coordinate system ('crd1' are omitted)
    Position o( ex_new_crd0*o_crd0,	 ey_new_crd0*o_crd0,	ez_new_crd0*o_crd0);
    Position h1(ex_new_crd0*h1_crd0, ey_new_crd0*h1_crd0, ez_new_crd0*h1_crd0);
    Position h2(ex_new_crd0*h2_crd0, ey_new_crd0*h2_crd0, ez_new_crd0*h2_crd0);
    // prev positions in new coordinate system ('crd1' are omitted)
    Position o_prev( ex_new_crd0*o_prev_crd0,	 ey_new_crd0*o_prev_crd0,	 ez_new_crd0*o_prev_crd0);
    Position h1_prev(ex_new_crd0*h1_prev_crd0, ey_new_crd0*h1_prev_crd0, ez_new_crd0*h1_prev_crd0);
    Position h2_prev(ex_new_crd0*h2_prev_crd0, ey_new_crd0*h2_prev_crd0, ez_new_crd0*h2_prev_crd0);

    // first and second rotations (psi and phi)
    double sin_phi = (o.z / model_OY);
    double cos_phi = sqrt(1 - sin_phi * sin_phi);
    double sin_psi = (h1.z - h2.z)/(2 * model_HX * cos_phi);
    double cos_psi = sqrt(1 - sin_psi * sin_psi);
    Position o_2(0.0,
                 model_OY * cos_phi,
                 model_OY * sin_phi);
    Position h1_2(-model_HX * cos_psi,
                  -model_HY * cos_phi - model_HX * sin_psi * sin_phi,
                  -model_HY * sin_phi + model_HX * sin_psi * cos_phi);
    Position h2_2(model_HX * cos_psi,
                  -model_HY * cos_phi + model_HX * sin_psi * sin_phi,
                  -model_HY * sin_phi - model_HX * sin_psi * cos_phi);

    // third rotation (theta)
    double alpha(h1_2.x * (h1_prev.x - h2_prev.x)
                 + h1_2.y * (h1_prev.y - o_prev.y)
                 + h2_2.y * (h2_prev.y - o_prev.y) );
    double beta( h1_2.x	* (h2_prev.y - h1_prev.y)
                 + h1_2.y * (h1_prev.x - o_prev.x)
                 + h2_2.y * (h2_prev.x - o_prev.x) );
    double gamma(h1.y *	(h1_prev.x - o_prev.x)
                 - h1.x * (h1_prev.y - o_prev.y)
                 + h2.y * (h2_prev.x - o_prev.x)
                 - h2.x * (h2_prev.y - o_prev.y) );
    double sin_theta = (alpha * gamma - beta * sqrt(alpha*alpha + beta*beta - gamma*gamma))
        / (alpha*alpha + beta*beta);
    double cos_theta = sqrt(1 - sin_theta*sin_theta);
    
    Position o_3( o_2.x * cos_theta - o_2.y * sin_theta,
                  o_2.x * sin_theta + o_2.y * cos_theta,
                  o_2.z);
    Position h1_3( h1_2.x * cos_theta - h1_2.y * sin_theta,
                   h1_2.x * sin_theta + h1_2.y * cos_theta,
                   h1_2.z);
    Position h2_3( h2_2.x * cos_theta - h2_2.y * sin_theta,
                   h2_2.x * sin_theta + h2_2.y * cos_theta,
                   h2_2.z);

    // new (constrained) positions in old coordinate system
    Position o_3_crd0( ex_old_crd1*o_3,	 ey_old_crd1*o_3,	 ez_old_crd1*o_3);
    Position h1_3_crd0(ex_old_crd1*h1_3, ey_old_crd1*h1_3, ez_old_crd1*h1_3);
    Position h2_3_crd0(ex_old_crd1*h2_3, ey_old_crd1*h2_3, ez_old_crd1*h2_3);

    // update position
    getpos(particle_array,o_index) = o_3_crd0 + com;
    getpos(particle_array,h1_index) = h1_3_crd0 + com;
    getpos(particle_array,h2_index) = h2_3_crd0 + com;

    // update velocity
    getvelocity(particle_array,o_index) += dtInv * (o_3_crd0 - o_crd0);
    getvelocity(particle_array,h1_index) += dtInv * (h1_3_crd0 - h1_crd0);
    getvelocity(particle_array,h2_index) += dtInv * (h2_3_crd0 - h2_crd0);

    water_settle_force[o_index] += dt2Inv2 * (o_3_crd0 - o_crd0)*getmass(particle_array,o_index);
    water_settle_force[h1_index] += dt2Inv2 * (h1_3_crd0 - h1_crd0)*getmass(particle_array,h1_index);
    water_settle_force[h2_index] += dt2Inv2 * (h2_3_crd0 - h2_crd0)*getmass(particle_array,h2_index);
  }

}
template
void Settle::constrain_position(CombinedParticleArray& particle_array,
                                const PosVelArray& particle_array_prev,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
                                const double dt,
				ForceArray& water_settle_force);
template<class PA>
void Settle::constrain_velocity(PA& particle_array,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
				const double dt,
				ForceArray& water_settle_force)
{
  double dtInv = 1.0 / dt;
  //	for(int i=type_range.ljcoulomb.begin;i<type_range.ljcoulomb.end;i++) {
  //		int id_wh = type_range.coulomb.begin + (i - type_range.ljcoulomb.begin) * 2;
  //	 for(int iw=0; iw<water_list.size(); ++iw){
  //		 int o_index = water_list[iw].o;
  //		 int h1_index = water_list[iw].h1;
  //		 int h2_index = water_list[iw].h2;

  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  //	for(std::map<int, WHPair>::const_iterator it=water_list.begin(); it != water_list.end(); ++it){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;

    // relative velocity
    Velocity vel_oh1(getvelocity(particle_array,h1_index) - getvelocity(particle_array,o_index));
    Velocity vel_oh2(getvelocity(particle_array,h2_index) - getvelocity(particle_array,o_index));
    Velocity vel_h1h2(getvelocity(particle_array,h2_index) - getvelocity(particle_array,h1_index));
	
    Position unit_vec_oh1(getpos(particle_array,h1_index) - getpos(particle_array,o_index));
    unit_vec_oh1 /= unit_vec_oh1.norm();
    Position unit_vec_oh2(getpos(particle_array,h2_index) - getpos(particle_array,o_index));
    unit_vec_oh2 /= unit_vec_oh2.norm();
    Position unit_vec_h1h2(getpos(particle_array,h2_index) - getpos(particle_array,h1_index));
    unit_vec_h1h2 /= unit_vec_h1h2.norm();
	
    double cos_O	= (unit_vec_oh1 * unit_vec_oh2);
    double cos_H1 = -(unit_vec_oh1 * unit_vec_h1h2);
    double cos_H2 = (unit_vec_oh2 * unit_vec_h1h2);
	
    double vel0_oh1 = unit_vec_oh1 * vel_oh1;
    double vel0_oh2 = unit_vec_oh2 * vel_oh2;
    double vel0_h1h2 = unit_vec_h1h2 * vel_h1h2;

    // determinant
    double det =
        ( 2 * massOH * massOH
          + 2 * massO * massH * cos_O * cos_H1 * cos_H2
          - 2 * massH * massH * cos_O * cos_O
          - massO * massOH * ( (cos_H1 * cos_H1) + (cos_H2 * cos_H2) ))
        * 0.5 * massHinv;
    double det_inv = 1.0 / det;

    double tau_OH1 = 
        ( (vel0_oh1	 * (2*massOH - massO * cos_H2 * cos_H2))
          +(vel0_oh2 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_h1h2 * (massH * cos_O * cos_H2 - massOH * cos_H1))
          ) * massO * det_inv;
    double tau_OH2 =
        ( (vel0_oh1 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_oh2 * (2 * massOH - massO * cos_H1 * cos_H1))
          +(vel0_h1h2 * (massH * cos_O * cos_H1 - massOH * cos_H2))
          ) * massO * det_inv;
    double tau_H1H2 =
        ( (vel0_oh1 * massO * (massH * cos_O * cos_H2 - massOH * cos_H1))
          +(vel0_oh2 * massO * (massH * cos_O * cos_H1 - massOH * cos_H2))
          +(vel0_h1h2 * (massOH * massOH - massH * massH * cos_O * cos_O))
          ) * det_inv;
	
    //update
    getvelocity(particle_array,o_index)	+= 0.5 * massOinv * ( tau_OH1 * unit_vec_oh1 + tau_OH2 * unit_vec_oh2);
    getvelocity(particle_array,h1_index) += 0.5 * massHinv * (-tau_OH1 * unit_vec_oh1 + tau_H1H2 * unit_vec_h1h2);
    getvelocity(particle_array,h2_index) += 0.5 * massHinv * (-tau_OH2 * unit_vec_oh2 - tau_H1H2 * unit_vec_h1h2);
    
    water_settle_force[o_index] = dtInv*( tau_OH1 * unit_vec_oh1 + tau_OH2 * unit_vec_oh2);
    water_settle_force[h1_index] = dtInv*(-tau_OH1 * unit_vec_oh1 + tau_H1H2 * unit_vec_h1h2);
    water_settle_force[h2_index] = dtInv*(-tau_OH2 * unit_vec_oh2 - tau_H1H2 * unit_vec_h1h2);
  }
}
template
void Settle::constrain_velocity(CombinedParticleArray& particle_array,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
				const double dt,
				ForceArray& water_settle_force);

template<class PA>
void Settle::constrain_velocity(PA& particle_array,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
				const double dt,
				double& virial)
{
  double dtInv = 1.0 / dt;
  double v = virial;
  //	for(int i=type_range.ljcoulomb.begin;i<type_range.ljcoulomb.end;i++) {
  //		int id_wh = type_range.coulomb.begin + (i - type_range.ljcoulomb.begin) * 2;
  //	 for(int iw=0; iw<water_list.size(); ++iw){
  //		 int o_index = water_list[iw].o;
  //		 int h1_index = water_list[iw].h1;
  //		 int h2_index = water_list[iw].h2;

  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  //	for(std::map<int, WHPair>::const_iterator it=water_list.begin(); it != water_list.end(); ++it){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;

    // relative velocity
    Velocity vel_oh1(getvelocity(particle_array,h1_index) - getvelocity(particle_array,o_index));
    Velocity vel_oh2(getvelocity(particle_array,h2_index) - getvelocity(particle_array,o_index));
    Velocity vel_h1h2(getvelocity(particle_array,h2_index) - getvelocity(particle_array,h1_index));
	
    Position unit_vec_oh1(getpos(particle_array,h1_index) - getpos(particle_array,o_index));
    double norm_oh1 = unit_vec_oh1.norm();
    unit_vec_oh1 /= norm_oh1;
    Position unit_vec_oh2(getpos(particle_array,h2_index) - getpos(particle_array,o_index));
    double norm_oh2 = unit_vec_oh2.norm();
    unit_vec_oh2 /= norm_oh2;
    Position unit_vec_h1h2(getpos(particle_array,h2_index) - getpos(particle_array,h1_index));
    double norm_h1h2 = unit_vec_h1h2.norm();
    unit_vec_h1h2 /= norm_h1h2;
	
    double cos_O	= (unit_vec_oh1 * unit_vec_oh2);
    double cos_H1 = -(unit_vec_oh1 * unit_vec_h1h2);
    double cos_H2 = (unit_vec_oh2 * unit_vec_h1h2);
	
    double vel0_oh1 = unit_vec_oh1 * vel_oh1;
    double vel0_oh2 = unit_vec_oh2 * vel_oh2;
    double vel0_h1h2 = unit_vec_h1h2 * vel_h1h2;

    // determinant
    double det =
        ( 2 * massOH * massOH
          + 2 * massO * massH * cos_O * cos_H1 * cos_H2
          - 2 * massH * massH * cos_O * cos_O
          - massO * massOH * ( (cos_H1 * cos_H1) + (cos_H2 * cos_H2) ))
        * 0.5 * massHinv;
    double det_inv = 1.0 / det;

    double tau_OH1 = 
        ( (vel0_oh1	 * (2*massOH - massO * cos_H2 * cos_H2))
          +(vel0_oh2 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_h1h2 * (massH * cos_O * cos_H2 - massOH * cos_H1))
          ) * massO * det_inv;
    double tau_OH2 =
        ( (vel0_oh1 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_oh2 * (2 * massOH - massO * cos_H1 * cos_H1))
          +(vel0_h1h2 * (massH * cos_O * cos_H1 - massOH * cos_H2))
          ) * massO * det_inv;
    double tau_H1H2 =
        ( (vel0_oh1 * massO * (massH * cos_O * cos_H2 - massOH * cos_H1))
          +(vel0_oh2 * massO * (massH * cos_O * cos_H1 - massOH * cos_H2))
          +(vel0_h1h2 * (massOH * massOH - massH * massH * cos_O * cos_O))
          ) * det_inv;
	
    //update
    getvelocity(particle_array,o_index)	+= 0.5 * massOinv * ( tau_OH1 * unit_vec_oh1 + tau_OH2 * unit_vec_oh2);
    getvelocity(particle_array,h1_index) += 0.5 * massHinv * (-tau_OH1 * unit_vec_oh1 + tau_H1H2 * unit_vec_h1h2);
    getvelocity(particle_array,h2_index) += 0.5 * massHinv * (-tau_OH2 * unit_vec_oh2 - tau_H1H2 * unit_vec_h1h2);
    /*
    water_settle_force[o_index] = dtInv*( tau_OH1 * unit_vec_oh1 + tau_OH2 * unit_vec_oh2);
    water_settle_force[h1_index] = dtInv*(-tau_OH1 * unit_vec_oh1 + tau_H1H2 * unit_vec_h1h2);
    water_settle_force[h2_index] = dtInv*(-tau_OH2 * unit_vec_oh2 - tau_H1H2 * unit_vec_h1h2);

    force
    OH1   dtInv*tau_OH1*unit_vec_oh1
    OH2   dtInv*tau_OH2*unit_vec_oh2
    H1O  -dtInv*tau_OH1*unit_vec_oh1
    H1H2  dtInv*tau_H1H2*unit_vec_h1h2
    H2O  -dtInv*tau_OH2*unit_vec_oh2
    H2H1 -dtInv*tau_H1H2*unit_vec_h1h2

    virial
    0.5*( vec_h1o*force_OH1+vec_h2o*force_OH2
         +vec_oh1*force_H1O+vec_h2h1*forceH1H2
         +vec_oh2*force_H2O+vec_h1h2*forceH2H1)
    = 0.5*( -vec_oh1*force_OH1-vec_oh2*force_OH2
         -vec_oh1*force_OH1-vec_h1h2*forceH1H2
         -vec_oh2*force_OH2-vec_h1h2*forceH1H2)
    = - vec_oh1*force_OH1 - vec_oh2*force_OH2 - vec_h1h2*forceH1H2
    = - dtInv*(vec_oh1*tau_OH1*unit_vec_oh1 + vec_oh2*tau_OH2*unit_vec_oh2 + vec_h1h2*tau_H1H2*unit_vec_h1h2)
    = - dtInv*(tau_OH1*vec_oh1*vec_oh1/norm(vec_oh1) + ...)
    = - dtInv*(tau_OH1*norm(vec_oh1) + tau_OH2*norm(vec_oh2) + tau_H1H2*norm(vec_h1h2)
    */
    v -= dtInv*(tau_OH1*norm_oh1 + tau_OH2*norm_oh2 + tau_H1H2*norm_h1h2);
  }
  virial = v;
}
template
void Settle::constrain_velocity(CombinedParticleArray& particle_array,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
				const double dt,
				double& virial);

void Settle::constrain_velocity(CombinedParticleArray& particle_array,
                                const WaterList& waterlist, 
                                const TypeRange& type_range)
{
  //	for(int i=type_range.ljcoulomb.begin;i<type_range.ljcoulomb.end;i++) {
  //		int id_wh = type_range.coulomb.begin + (i - type_range.ljcoulomb.begin) * 2;
  //	 for(int iw=0; iw<water_list.size(); ++iw){
  //		 int o_index = water_list[iw].o;
  //		 int h1_index = water_list[iw].h1;
  //		 int h2_index = water_list[iw].h2;

  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  //	for(std::map<int, WHPair>::const_iterator it=water_list.begin(); it != water_list.end(); ++it){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;

    // relative velocity
    Velocity vel_oh1(particle_array.parameters[h1_index].velocity - particle_array.parameters[o_index].velocity);
    Velocity vel_oh2(particle_array.parameters[h2_index].velocity - particle_array.parameters[o_index].velocity);
    Velocity vel_h1h2(particle_array.parameters[h2_index].velocity - particle_array.parameters[h1_index].velocity);
	
    Position unit_vec_oh1(particle_array.poscharge[h1_index].position - particle_array.poscharge[o_index].position);
    unit_vec_oh1 /= unit_vec_oh1.norm();
    Position unit_vec_oh2(particle_array.poscharge[h2_index].position - particle_array.poscharge[o_index].position);
    unit_vec_oh2 /= unit_vec_oh2.norm();
    Position unit_vec_h1h2(particle_array.poscharge[h2_index].position - particle_array.poscharge[h1_index].position);
    unit_vec_h1h2 /= unit_vec_h1h2.norm();
	
    double cos_O	= (unit_vec_oh1 * unit_vec_oh2);
    double cos_H1 = -(unit_vec_oh1 * unit_vec_h1h2);
    double cos_H2 = (unit_vec_oh2 * unit_vec_h1h2);
	
    double vel0_oh1 = unit_vec_oh1 * vel_oh1;
    double vel0_oh2 = unit_vec_oh2 * vel_oh2;
    double vel0_h1h2 = unit_vec_h1h2 * vel_h1h2;

    // determinant
    double det =
        ( 2 * massOH * massOH
          + 2 * massO * massH * cos_O * cos_H1 * cos_H2
          - 2 * massH * massH * cos_O * cos_O
          - massO * massOH * ( (cos_H1 * cos_H1) + (cos_H2 * cos_H2) ))
        * 0.5 * massHinv;
    double det_inv = 1.0 / det;

    double tau_OH1 = 
        ( (vel0_oh1	 * (2*massOH - massO * cos_H2 * cos_H2))
          +(vel0_oh2 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_h1h2 * (massH * cos_O * cos_H2 - massOH * cos_H1))
          ) * massO * det_inv;
    double tau_OH2 =
        ( (vel0_oh1 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_oh2 * (2 * massOH - massO * cos_H1 * cos_H1))
          +(vel0_h1h2 * (massH * cos_O * cos_H1 - massOH * cos_H2))
          ) * massO * det_inv;
    double tau_H1H2 =
        ( (vel0_oh1 * massO * (massH * cos_O * cos_H2 - massOH * cos_H1))
          +(vel0_oh2 * massO * (massH * cos_O * cos_H1 - massOH * cos_H2))
          +(vel0_h1h2 * (massOH * massOH - massH * massH * cos_O * cos_O))
          ) * det_inv;
	
    //update
    particle_array.parameters[o_index].velocity	+= 0.5 * massOinv * ( tau_OH1 * unit_vec_oh1 + tau_OH2 * unit_vec_oh2);
    particle_array.parameters[h1_index].velocity += 0.5 * massHinv * (-tau_OH1 * unit_vec_oh1 + tau_H1H2 * unit_vec_h1h2);
    particle_array.parameters[h2_index].velocity += 0.5 * massHinv * (-tau_OH2 * unit_vec_oh2 - tau_H1H2 * unit_vec_h1h2);
  }
}

void Settle::constrain_position(ParticleArray& particle_array,
                                const ParticleArray& particle_array_prev,
                                const WaterList& waterlist, 
                                const TypeRange& type_range,
                                double dt)
{
  using namespace std;
  double dtInv = 1.0 / dt;
  
  // DEBUG
  //  double accel_max=0.0;

  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
		
    // _crd0 : coordinates in world (global) coordinate system
    // _crd1 or none : coordinates in local coordinate system based on the molecule position

    //vector from the center of mass to new (unconstrained) {O,H1,H2} in old coordinate system
    Position com = ((massO * particle_array[o_index].position)
                    +(massH * particle_array[h1_index].position)
                    +(massH * particle_array[h2_index].position)) * massH2Oinv;
    Position o_crd0(particle_array[o_index].position - com);
    Position h1_crd0(particle_array[h1_index].position - com);
    Position h2_crd0(particle_array[h2_index].position - com);

    //vector from the 'com' to prev {O,H1,H2} in old coordinate system
    Position o_prev_crd0(particle_array_prev[o_index].position - com);
    Position h1_prev_crd0(particle_array_prev[h1_index].position - com);
    Position h2_prev_crd0(particle_array_prev[h2_index].position - com);

    // unit vectors of new axes in old coordinate system
    Position ez_new_crd0((h1_prev_crd0 - o_prev_crd0)%(h2_prev_crd0 - o_prev_crd0) );
    ez_new_crd0 /= ez_new_crd0.norm();
    Position ex_new_crd0(o_crd0%ez_new_crd0);
    ex_new_crd0 /= ex_new_crd0.norm();
    Position ey_new_crd0(ez_new_crd0%ex_new_crd0);
    ey_new_crd0 /= ey_new_crd0.norm();
    // unit vectors of old axes in new coordinate system
    Position ex_old_crd1(ex_new_crd0.x, ey_new_crd0.x, ez_new_crd0.x);
    Position ey_old_crd1(ex_new_crd0.y, ey_new_crd0.y, ez_new_crd0.y);
    Position ez_old_crd1(ex_new_crd0.z, ey_new_crd0.z, ez_new_crd0.z);

    // positions in new coordinate system ('crd1' are omitted)
    Position o( ex_new_crd0*o_crd0,	 ey_new_crd0*o_crd0,	ez_new_crd0*o_crd0);
    Position h1(ex_new_crd0*h1_crd0, ey_new_crd0*h1_crd0, ez_new_crd0*h1_crd0);
    Position h2(ex_new_crd0*h2_crd0, ey_new_crd0*h2_crd0, ez_new_crd0*h2_crd0);
    // prev positions in new coordinate system ('crd1' are omitted)
    Position o_prev( ex_new_crd0*o_prev_crd0,	 ey_new_crd0*o_prev_crd0,	 ez_new_crd0*o_prev_crd0);
    Position h1_prev(ex_new_crd0*h1_prev_crd0, ey_new_crd0*h1_prev_crd0, ez_new_crd0*h1_prev_crd0);
    Position h2_prev(ex_new_crd0*h2_prev_crd0, ey_new_crd0*h2_prev_crd0, ez_new_crd0*h2_prev_crd0);

    // first and second rotations (psi and phi)
    double sin_phi = (o.z / model_OY);
    double cos_phi = sqrt(1 - sin_phi * sin_phi);
    double sin_psi = (h1.z - h2.z)/(2 * model_HX * cos_phi);
    double cos_psi = sqrt(1 - sin_psi * sin_psi);
    Position o_2(0.0,
                 model_OY * cos_phi,
                 model_OY * sin_phi);
    Position h1_2(-model_HX * cos_psi,
                  -model_HY * cos_phi - model_HX * sin_psi * sin_phi,
                  -model_HY * sin_phi + model_HX * sin_psi * cos_phi);
    Position h2_2(model_HX * cos_psi,
                  -model_HY * cos_phi + model_HX * sin_psi * sin_phi,
                  -model_HY * sin_phi - model_HX * sin_psi * cos_phi);

    // third rotation (theta)
    double alpha(h1_2.x * (h1_prev.x - h2_prev.x)
                 + h1_2.y * (h1_prev.y - o_prev.y)
                 + h2_2.y * (h2_prev.y - o_prev.y) );
    double beta( h1_2.x	* (h2_prev.y - h1_prev.y)
                 + h1_2.y * (h1_prev.x - o_prev.x)
                 + h2_2.y * (h2_prev.x - o_prev.x) );
    double gamma(h1.y *	(h1_prev.x - o_prev.x)
                 - h1.x * (h1_prev.y - o_prev.y)
                 + h2.y * (h2_prev.x - o_prev.x)
                 - h2.x * (h2_prev.y - o_prev.y) );
    double sin_theta = (alpha * gamma - beta * sqrt(alpha*alpha + beta*beta - gamma*gamma))
        / (alpha*alpha + beta*beta);
    double cos_theta = sqrt(1 - sin_theta*sin_theta);
    
    Position o_3( o_2.x * cos_theta - o_2.y * sin_theta,
                  o_2.x * sin_theta + o_2.y * cos_theta,
                  o_2.z);
    Position h1_3( h1_2.x * cos_theta - h1_2.y * sin_theta,
                   h1_2.x * sin_theta + h1_2.y * cos_theta,
                   h1_2.z);
    Position h2_3( h2_2.x * cos_theta - h2_2.y * sin_theta,
                   h2_2.x * sin_theta + h2_2.y * cos_theta,
                   h2_2.z);

    // new (constrained) positions in old coordinate system
    Position o_3_crd0( ex_old_crd1*o_3,	 ey_old_crd1*o_3,	 ez_old_crd1*o_3);
    Position h1_3_crd0(ex_old_crd1*h1_3, ey_old_crd1*h1_3, ez_old_crd1*h1_3);
    Position h2_3_crd0(ex_old_crd1*h2_3, ey_old_crd1*h2_3, ez_old_crd1*h2_3);

    // update position
    particle_array[o_index].position = o_3_crd0 + com;
    particle_array[h1_index].position = h1_3_crd0 + com;
    particle_array[h2_index].position = h2_3_crd0 + com;

    // update velocity
    particle_array[o_index].velocity += dtInv * (o_3_crd0 - o_crd0);
    particle_array[h1_index].velocity += dtInv * (h1_3_crd0 - h1_crd0);
    particle_array[h2_index].velocity += dtInv * (h2_3_crd0 - h2_crd0);

    //DEBUG
    /*
    {
      Position dv = (o_3_crd0 - o_crd0);
      accel_max = std::max(accel_max,dv.norm2());
      dv = (h1_3_crd0 - h1_crd0);
      accel_max = std::max(accel_max,dv.norm2());
      dv = (h2_3_crd0 - h2_crd0);
      accel_max = std::max(accel_max,dv.norm2());
    }
    */
  }
  //DEBUG
  //  std::cout << " max settle accel " << accel_max << std::endl;
}

void Settle::constrain_velocity(ParticleArray& particle_array,
                                const WaterList& waterlist, 
                                const TypeRange& type_range)
{
  //	for(int i=type_range.ljcoulomb.begin;i<type_range.ljcoulomb.end;i++) {
  //		int id_wh = type_range.coulomb.begin + (i - type_range.ljcoulomb.begin) * 2;
  //	 for(int iw=0; iw<water_list.size(); ++iw){
  //		 int o_index = water_list[iw].o;
  //		 int h1_index = water_list[iw].h1;
  //		 int h2_index = water_list[iw].h2;

  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  //	for(std::map<int, WHPair>::const_iterator it=water_list.begin(); it != water_list.end(); ++it){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;

    // relative velocity
    Velocity vel_oh1(particle_array[h1_index].velocity - particle_array[o_index].velocity);
    Velocity vel_oh2(particle_array[h2_index].velocity - particle_array[o_index].velocity);
    Velocity vel_h1h2(particle_array[h2_index].velocity - particle_array[h1_index].velocity);
	
    Position unit_vec_oh1(particle_array[h1_index].position - particle_array[o_index].position);
    unit_vec_oh1 /= unit_vec_oh1.norm();
    Position unit_vec_oh2(particle_array[h2_index].position - particle_array[o_index].position);
    unit_vec_oh2 /= unit_vec_oh2.norm();
    Position unit_vec_h1h2(particle_array[h2_index].position - particle_array[h1_index].position);
    unit_vec_h1h2 /= unit_vec_h1h2.norm();
	
    double cos_O	= (unit_vec_oh1 * unit_vec_oh2);
    double cos_H1 = -(unit_vec_oh1 * unit_vec_h1h2);
    double cos_H2 = (unit_vec_oh2 * unit_vec_h1h2);
	
    double vel0_oh1 = unit_vec_oh1 * vel_oh1;
    double vel0_oh2 = unit_vec_oh2 * vel_oh2;
    double vel0_h1h2 = unit_vec_h1h2 * vel_h1h2;

    // determinant
    double det =
        ( 2 * massOH * massOH
          + 2 * massO * massH * cos_O * cos_H1 * cos_H2
          - 2 * massH * massH * cos_O * cos_O
          - massO * massOH * ( (cos_H1 * cos_H1) + (cos_H2 * cos_H2) ))
        * 0.5 * massHinv;
    double det_inv = 1.0 / det;

    double tau_OH1 = 
        ( (vel0_oh1	 * (2*massOH - massO * cos_H2 * cos_H2))
          +(vel0_oh2 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_h1h2 * (massH * cos_O * cos_H2 - massOH * cos_H1))
          ) * massO * det_inv;
    double tau_OH2 =
        ( (vel0_oh1 * (massO * cos_H1 * cos_H2 - 2 * massH * cos_O))
          +(vel0_oh2 * (2 * massOH - massO * cos_H1 * cos_H1))
          +(vel0_h1h2 * (massH * cos_O * cos_H1 - massOH * cos_H2))
          ) * massO * det_inv;
    double tau_H1H2 =
        ( (vel0_oh1 * massO * (massH * cos_O * cos_H2 - massOH * cos_H1))
          +(vel0_oh2 * massO * (massH * cos_O * cos_H1 - massOH * cos_H2))
          +(vel0_h1h2 * (massOH * massOH - massH * massH * cos_O * cos_O))
          ) * det_inv;
	
    //update
    //    std::cout << " settle_v O H1 H2 " << particle_array[o_index].velocity << particle_array[h1_index].velocity << particle_array[h2_index].velocity << std::endl;
    particle_array[o_index].velocity	+= 0.5 * massOinv * ( tau_OH1 * unit_vec_oh1 + tau_OH2 * unit_vec_oh2);
    particle_array[h1_index].velocity += 0.5 * massHinv * (-tau_OH1 * unit_vec_oh1 + tau_H1H2 * unit_vec_h1h2);
    particle_array[h2_index].velocity += 0.5 * massHinv * (-tau_OH2 * unit_vec_oh2 - tau_H1H2 * unit_vec_h1h2);
    //    std::cout << " settle_v O H1 H2 " << particle_array[o_index].velocity << particle_array[h1_index].velocity << particle_array[h2_index].velocity << std::endl;
  }
}

void Settle::rattle_position(ParticleArray& particle_array,
                             const ParticleArray& particle_array_prev,
                             const WaterList& waterlist, 
                             const TypeRange& type_range,
                             double dt)
{
  int maxIterate = 3000;

  //	 for(int i=type_range.ljcoulomb.begin;i<type_range.ljcoulomb.end;i++) {
  //		 int id_wh = type_range.coulomb.begin + (i - type_range.ljcoulomb.begin) * 2;
  //	 for(int iw=0; iw<water_list.size(); ++iw){
  //		 int o_index = water_list[iw].o;
  //		 int h1_index = water_list[iw].h1;
  //		 int h2_index = water_list[iw].h2;
  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;

    Position oh1_prev(particle_array_prev[h1_index].position - particle_array_prev[o_index].position);
    Position oh2_prev(particle_array_prev[h2_index].position - particle_array_prev[o_index].position);
    Position h1h2_prev(particle_array_prev[h2_index].position - particle_array_prev[h1_index].position);

    bool is_modified = true;

    for (int numIterate=0; is_modified && numIterate< maxIterate; ++numIterate) {
      is_modified = false;

      Position oh1(particle_array[h1_index].position - particle_array[o_index].position);
      double diff2_oh1 = oh1.norm2() - length2OH;
      if(std::abs(diff2_oh1) > tolerance2PosOH){
        double g = diff2_oh1 / (2.0 * dt * massOinvHinv * (oh1 * oh1_prev));
        particle_array[o_index].position	+= g * dt * massOinv * oh1_prev;
        particle_array[h1_index].position -= g * dt * massHinv * oh1_prev;
        particle_array[o_index].velocity	+= g * massOinv * oh1_prev;
        particle_array[h1_index].velocity -= g * massHinv * oh1_prev;
        is_modified = true;
      }

      Position oh2(particle_array[h2_index].position - particle_array[o_index].position);
      double diff2_oh2 = oh2.norm2() - length2OH;
      if(std::abs(diff2_oh2) > tolerance2PosOH){
        double g = diff2_oh2 / (2.0 * dt * massOinvHinv * (oh2 * oh2_prev));
        particle_array[o_index].position	+= g * dt * massOinv * oh2_prev;
        particle_array[h2_index].position -= g * dt * massHinv * oh2_prev;
        particle_array[o_index].velocity	+= g * massOinv * oh2_prev;
        particle_array[h2_index].velocity -= g * massHinv * oh2_prev;
        is_modified = true;
      }

      Position h1h2(particle_array[h2_index].position - particle_array[h1_index].position);
      double diff2_h1h2 = h1h2.norm2() - length2HH;
      if(std::abs(diff2_h1h2) > tolerance2PosHH){
        double g = diff2_h1h2 / (2.0 * dt * massHinvHinv * (h1h2 * h1h2_prev));
        particle_array[h1_index].position		+= g * dt * massHinv * h1h2_prev;
        particle_array[h2_index].position -= g * dt * massHinv * h1h2_prev;
        particle_array[h1_index].velocity		+= g * massHinv * h1h2_prev;
        particle_array[h2_index].velocity -= g * massHinv * h1h2_prev;
        is_modified = true;
      }
    }
  }			 
}


void Settle::rattle_velocity(ParticleArray& particle_array, 
                             const WaterList& waterlist, 
                             const TypeRange& type_range, double dt)
{
  toleranceVelOH = tolerance * length2OH / dt;
  toleranceVelHH = tolerance * length2HH / dt;

  int maxIterate = 3000;

  //	 for(int i=type_range.ljcoulomb.begin;i<type_range.ljcoulomb.end;i++) {
  //		 int h1_index = type_range.coulomb.begin + (i - type_range.ljcoulomb.begin) * 2;
  //	 for(int iw=0; iw<water_list.size(); ++iw){
  //		 int o_index = water_list[iw].o;
  //		 int h1_index = water_list[iw].h1;
  //		 int h2_index = water_list[iw].h2;
  // 	map<int,WHPair> wlist = waterlist.water_list;
  // 	for(map<int,WHPair>::iterator it=wlist.begin(); it != wlist.end();it++ ){
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end();it++ ){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;

    bool is_modified = true;
    for (int numIterate=0; is_modified && numIterate< maxIterate; ++numIterate) {
      is_modified = false;
			
      Position oh1_pos(particle_array[h1_index].position - particle_array[o_index].position);
      Velocity oh1_vel(particle_array[h1_index].velocity - particle_array[o_index].velocity);
      double oh1_pv = oh1_pos * oh1_vel;
      if (std::abs(oh1_pv) > toleranceVelOH){
        double k = oh1_pv / (massOinvHinv * oh1_pos.norm2());
        particle_array[o_index].velocity	+= k * massOinv * oh1_pos;
        particle_array[h1_index].velocity -= k * massHinv * oh1_pos;
        is_modified = true;
      }
			
      Position oh2_pos(particle_array[h2_index].position - particle_array[o_index].position);
      Velocity oh2_vel(particle_array[h2_index].velocity - particle_array[o_index].velocity);
      double oh2_pv = oh2_pos * oh2_vel;
      if (std::abs(oh2_pv) > toleranceVelOH){
        double k = oh2_pv / (massOinvHinv * oh2_pos.norm2());
        particle_array[o_index].velocity	+= k * massOinv * oh2_pos;
        particle_array[h2_index].velocity -= k * massHinv * oh2_pos;
        is_modified = true;
      }
			
      Position h1h2_pos(particle_array[h2_index].position - particle_array[h1_index].position);
      Velocity h1h2_vel(particle_array[h2_index].velocity - particle_array[h1_index].velocity);
      double h1h2_pv = h1h2_pos * h1h2_vel;
      if (std::abs(h1h2_pv) > toleranceVelHH){  /// fix toleranceVelOH to toleranceVelHH ?
        double k = h1h2_pv / (massHinvHinv * h1h2_pos.norm2());
        particle_array[h1_index].velocity		+= k * massHinv * h1h2_pos;
        particle_array[h2_index].velocity -= k * massHinv * h1h2_pos;
        is_modified = true;
      }
    }
  }
}

template<class PA, class PVA>
void Settle::rattle_position(PA& particle_array,
                             const PVA& particle_array_prev,
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance)
{
  double shake_tolerance2 = 2.0 * shake_tolerance;
  double dt2= 2.0 * dt;

#if DEBUG_SHAKE_COUNT
  int num_shake = 0;
  int num_shake_bond = 0;
  int num_search_length = 0;
  int num_iter=0;
#endif

  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int ha_index = it->first;
    double massHAinv = getinvmass(particle_array,ha_index);
    double dt_massHAinv = dt * massHAinv;
    Position ha_pos(getpos(particle_array,ha_index));
    Velocity ha_vel(getvelocity(particle_array,ha_index));

    int nh1 = it->second.nh1;
#if DEBUG_SHAKE_COUNT
    num_shake++;
#endif

    double length2HAH1[MAX_HBOND];
    double tolerance2HAH1[MAX_HBOND];
    
    Position z0[MAX_HBOND];
    Position za[MAX_HBOND];
    Position z1[MAX_HBOND];
    Position zza[MAX_HBOND];
    Position zz1[MAX_HBOND];
    Position h1_pos[MAX_HBOND];
    Velocity h1_vel[MAX_HBOND];

    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      int h1_index = it->second.h1[n1];

#if DEBUG_SHAKE_COUNT
      num_shake_bond++;
#endif
      double lengthHAH1 = 1.0;    // default
      bool get_lengthHAH1 = false;
      int bondtype = it->second.bondtype[n1];
      if(bondtype==-1){
	std::cout <<  "Waning rattle_position : search bondtype for index(" <<ha_index << "," << h1_index << ")" << std::endl;
	for(size_t s=0; !get_lengthHAH1 && s<bondlistarray.size(); s++){
	  for(size_t sb=0; !get_lengthHAH1 && sb<bondlistarray[s].BondArray.size(); sb++){
	    CovalentBondInfo::Bond bond = bondlistarray[s].BondArray[sb];
	    
#if DEBUG_SHAKE_COUNT
	    num_search_length++;
#endif
	    if(!get_lengthHAH1 &&
	       ((getatomid(particle_array,ha_index) == bond.id_of_atom[0] && getatomid(particle_array,h1_index) == bond.id_of_atom[1])
		|| (getatomid(particle_array,ha_index) == bond.id_of_atom[1] && getatomid(particle_array,h1_index) == bond.id_of_atom[0]))){
	      //	      lengthHAH1 = param_list->bond[bond.typeofbond].equilibrium_length;
	      bondtype = bond.typeofbond;
	      get_lengthHAH1 = true;
	      break;
	    }
	  }
	}
      }else{
	get_lengthHAH1 = true;
      }

      if(!get_lengthHAH1){
        std::cout << "ERROR: in Settle::rattle_position : not hit bond" << std::endl;
        continue;
      }else{
	lengthHAH1 = param_list->bond[bondtype].equilibrium_length;
      }

#if DEBUG_SHAKE
      std::cout << "in Settle::rattle_position lengthHAH1=" << lengthHAH1 << " ha_index=" << ha_index << " h1_index=" << h1_index << std::endl;
#endif

      length2HAH1[n1] = lengthHAH1 * lengthHAH1;
      tolerance2HAH1[n1] = shake_tolerance2 * length2HAH1[n1];
      double massH1inv = getinvmass(particle_array,h1_index);
      double dt_massH1inv = dt * massH1inv;
      double massHAinvH1inv = massHAinv + massH1inv;
      
      h1_pos[n1] = getpos(particle_array,h1_index);
      h1_vel[n1] = getvelocity(particle_array,h1_index);
      Position hah1_prev(getpos(particle_array_prev,h1_index) - getpos(particle_array_prev,ha_index));
      z0[n1] = dt2 * massHAinvH1inv * hah1_prev;
      za[n1] = dt_massHAinv * hah1_prev;
      z1[n1] = dt_massH1inv * hah1_prev;
      zza[n1] = massHAinv * hah1_prev;
      zz1[n1] = massH1inv * hah1_prev;
    }
    bool is_modified = true;
    for (int numIterate=0; is_modified && numIterate< shake_max_iterate; ++numIterate) {
      is_modified = false;

#if DEBUG_SHAKE_COUNT
      num_iter++;
#endif
      for(int n1=0; n1<nh1; n1++){
	Position hah1(h1_pos[n1] - ha_pos);
	double diff2_hah1 = hah1.norm2() - length2HAH1[n1];

#if DEBUG_SHAKE
	std::cout << "n=" << numIterate << " diff2_hah1=" << diff2_hah1 << " tolerance=" << tolerance << std::endl;
#endif

	if(std::abs(diff2_hah1) > tolerance2HAH1[n1]){
	  double g = diff2_hah1 / (hah1 * z0[n1]);
	  ha_pos += g * za[n1];
	  h1_pos[n1] -= g * z1[n1];
	  ha_vel += g * zza[n1];
	  h1_vel[n1] -= g * zz1[n1];
	  is_modified = true;
	}
      }
    }
    getpos(particle_array,ha_index) = ha_pos;
    getvelocity(particle_array,ha_index) = ha_vel;
    for(int n1=0; n1<nh1; n1++){
      int h1_index = it->second.h1[n1];
      getpos(particle_array,h1_index) = h1_pos[n1];
      getvelocity(particle_array,h1_index) = h1_vel[n1];
    }
  }			 
#if DEBUG_SHAKE||DEBUG_SHAKE_COUNT
  std::cout << "in Settle::rattle_position num_shake=" << num_shake << " num_shake_bond=" <<  num_shake_bond << " num_search_length=" << num_search_length << " num_iter=" << num_iter << " average_iter" << double(num_iter)/num_shake << std::endl;
#endif
}
template
void Settle::rattle_position(CombinedParticleArray& particle_array,
                             const PosVelArray& particle_array_prev,
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance);
template
void Settle::rattle_position(ParticleArray& particle_array,
                             const PosVelArray& particle_array_prev,
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance);

template<class PA, class PVA>
void Settle::rattle_position(PA& particle_array,
                             const PVA& particle_array_prev,
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     ForceArray& rattle_force)
{
  double Inv2dt = 2.0/dt;
  double shake_tolerance2 = 2.0 * shake_tolerance;
  double dt2= 2.0 * dt;

#if DEBUG_SHAKE_COUNT
  int num_shake = 0;
  int num_shake_bond = 0;
  int num_search_length = 0;
  int num_iter=0;
#endif

  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int ha_index = it->first;
    double massHAinv = getinvmass(particle_array,ha_index);
    double dt_massHAinv = dt * massHAinv;
    Position ha_pos(getpos(particle_array,ha_index));
    Velocity ha_vel(getvelocity(particle_array,ha_index));

    int nh1 = it->second.nh1;
#if DEBUG_SHAKE_COUNT
    num_shake++;
#endif

    double length2HAH1[MAX_HBOND];
    double tolerance2HAH1[MAX_HBOND];
    
    Position z0[MAX_HBOND];
    Position za[MAX_HBOND];
    Position z1[MAX_HBOND];
    Position zza[MAX_HBOND];
    Position zz1[MAX_HBOND];
    Position h1_pos[MAX_HBOND];
    Velocity h1_vel[MAX_HBOND];

    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      int h1_index = it->second.h1[n1];

#if DEBUG_SHAKE_COUNT
      num_shake_bond++;
#endif
      double lengthHAH1 = 1.0;    // default
      bool get_lengthHAH1 = false;
      int bondtype = it->second.bondtype[n1];
      if(bondtype==-1){
	std::cout <<  "Waning rattle_position : search bondtype for index(" <<ha_index << "," << h1_index << ")" << std::endl;
	for(size_t s=0; !get_lengthHAH1 && s<bondlistarray.size(); s++){
	  for(size_t sb=0; !get_lengthHAH1 && sb<bondlistarray[s].BondArray.size(); sb++){
	    CovalentBondInfo::Bond bond = bondlistarray[s].BondArray[sb];
	    
#if DEBUG_SHAKE_COUNT
	    num_search_length++;
#endif
	    if(!get_lengthHAH1 &&
	       ((getatomid(particle_array,ha_index) == bond.id_of_atom[0] && getatomid(particle_array,h1_index) == bond.id_of_atom[1])
		|| (getatomid(particle_array,ha_index) == bond.id_of_atom[1] && getatomid(particle_array,h1_index) == bond.id_of_atom[0]))){
	      //	      lengthHAH1 = param_list->bond[bond.typeofbond].equilibrium_length;
	      bondtype = bond.typeofbond;
	      get_lengthHAH1 = true;
	      break;
	    }
	  }
	}
      }else{
	get_lengthHAH1 = true;
      }

      if(!get_lengthHAH1){
        std::cout << "ERROR: in Settle::rattle_position : not hit bond" << std::endl;
        continue;
      }else{
	lengthHAH1 = param_list->bond[bondtype].equilibrium_length;
      }

#if DEBUG_SHAKE
      std::cout << "in Settle::rattle_position lengthHAH1=" << lengthHAH1 << " ha_index=" << ha_index << " h1_index=" << h1_index << std::endl;
#endif

      length2HAH1[n1] = lengthHAH1 * lengthHAH1;
      tolerance2HAH1[n1] = shake_tolerance2 * length2HAH1[n1];
      double massH1inv = getinvmass(particle_array,h1_index);
      double dt_massH1inv = dt * massH1inv;
      double massHAinvH1inv = massHAinv + massH1inv;
      
      h1_pos[n1] = getpos(particle_array,h1_index);
      h1_vel[n1] = getvelocity(particle_array,h1_index);
      Position hah1_prev(getpos(particle_array_prev,h1_index) - getpos(particle_array_prev,ha_index));
      z0[n1] = dt2 * massHAinvH1inv * hah1_prev;
      za[n1] = dt_massHAinv * hah1_prev;
      z1[n1] = dt_massH1inv * hah1_prev;
      zza[n1] = massHAinv * hah1_prev;
      zz1[n1] = massH1inv * hah1_prev;
    }
    bool is_modified = true;
    for (int numIterate=0; is_modified && numIterate< shake_max_iterate; ++numIterate) {
      is_modified = false;

#if DEBUG_SHAKE_COUNT
      num_iter++;
#endif
      for(int n1=0; n1<nh1; n1++){
	Position hah1(h1_pos[n1] - ha_pos);
	double diff2_hah1 = hah1.norm2() - length2HAH1[n1];

#if DEBUG_SHAKE
	std::cout << "n=" << numIterate << " diff2_hah1=" << diff2_hah1 << " tolerance=" << tolerance << std::endl;
#endif

	if(std::abs(diff2_hah1) > tolerance2HAH1[n1]){
	  double g = diff2_hah1 / (hah1 * z0[n1]);
	  ha_pos += g * za[n1];
	  h1_pos[n1] -= g * z1[n1];
	  ha_vel += g * zza[n1];
	  h1_vel[n1] -= g * zz1[n1];
	  is_modified = true;
	}
      }
    }
    getpos(particle_array,ha_index) = ha_pos;
    rattle_force[ha_index] += (ha_vel - getvelocity(particle_array,ha_index))*Inv2dt*getmass(particle_array,ha_index);
    getvelocity(particle_array,ha_index) = ha_vel;
    for(int n1=0; n1<nh1; n1++){
      int h1_index = it->second.h1[n1];
      getpos(particle_array,h1_index) = h1_pos[n1];
      rattle_force[h1_index] += (h1_vel[n1] - getvelocity(particle_array,h1_index))*Inv2dt*getmass(particle_array,h1_index);
      getvelocity(particle_array,h1_index) = h1_vel[n1];
    }
  }			 
#if DEBUG_SHAKE||DEBUG_SHAKE_COUNT
  std::cout << "in Settle::rattle_position num_shake=" << num_shake << " num_shake_bond=" <<  num_shake_bond << " num_search_length=" << num_search_length << " num_iter=" << num_iter << " average_iter" << double(num_iter)/num_shake << std::endl;
#endif
}
template
void Settle::rattle_position(CombinedParticleArray& particle_array,
                             const PosVelArray& particle_array_prev,
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     ForceArray& rattle_force);
template
void Settle::rattle_position(ParticleArray& particle_array,
                             const PosVelArray& particle_array_prev,
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     ForceArray& rattle_force);

template<class PA>
void Settle::rattle_velocity(PA& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance)
{
  double tolerance_per_dt = shake_tolerance / dt;

  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int ha_index = it->first;
    double massHAinv = getinvmass(particle_array,ha_index);
    Velocity ha_vel(getvelocity(particle_array,ha_index));

    int nh1 = it->second.nh1;

    double tolerance2HAH1[MAX_HBOND];
    double massH1inv[MAX_HBOND];
    Position hah1_pos[MAX_HBOND];
    Velocity h1_vel[MAX_HBOND];
    double massHAinvH1invHAH1posnorm2inv[MAX_HBOND];

    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      int h1_index = it->second.h1[n1];

      double lengthHAH1 = 1.0;    // default
      bool get_lengthHAH1 = false;
      int bondtype = it->second.bondtype[n1];
      if(bondtype==-1){
	std::cout <<  "Waning rattle_position : search bondtype for index(" <<ha_index << "," << h1_index << ")" << std::endl;
	for(size_t s=0; !get_lengthHAH1 && s<bondlistarray.size(); s++){
	  for(size_t sb=0; !get_lengthHAH1 && sb<bondlistarray[s].BondArray.size(); sb++){
	    CovalentBondInfo::Bond bond = bondlistarray[s].BondArray[sb];

	    if(!get_lengthHAH1 &&
	       ((getatomid(particle_array,ha_index) == bond.id_of_atom[0] && getatomid(particle_array,h1_index) == bond.id_of_atom[1])
		|| (getatomid(particle_array,ha_index) == bond.id_of_atom[1] && getatomid(particle_array,h1_index) == bond.id_of_atom[0]))){
	      //	      lengthHAH1 = param_list->bond[bond.typeofbond].equilibrium_length;
	      bondtype = bond.typeofbond;
	      get_lengthHAH1 = true;
	      break;
	    }
	  }
	}
      }else{
	get_lengthHAH1 = true;
      }
      if(!get_lengthHAH1){
        std::cout << "ERROR: in Settle::rattle_velocity : not hit bond" << std::endl;
        continue;
      }else{
	lengthHAH1 = param_list->bond[bondtype].equilibrium_length;
      }

#if DEBUG_SHAKE
      std::cout << "in Settle::rattle_velocity lengthHAH1=" << lengthHAH1 << " ha_index=" << ha_index << " h1_index=" << h1_index << std::endl;
#endif

      double length2HAH1 = lengthHAH1 * lengthHAH1;
      tolerance2HAH1[n1] = tolerance_per_dt * length2HAH1;
      massH1inv[n1] = getinvmass(particle_array,h1_index);

      hah1_pos[n1] = (getpos(particle_array,h1_index) - getpos(particle_array,ha_index));
      h1_vel[n1] = getvelocity(particle_array,h1_index);
      double massHAinvH1inv = massHAinv + massH1inv[n1];
      massHAinvH1invHAH1posnorm2inv[n1] = 1.0/(massHAinvH1inv*hah1_pos[n1].norm2());
    }

    bool is_modified = true;
    for (int numIterate=0; is_modified && numIterate< shake_max_iterate; ++numIterate) {
      is_modified = false;

      for(int n1=0; n1<nh1; n1++){
        Velocity hah1_vel(h1_vel[n1] - ha_vel);
        double hah1_pv = hah1_pos[n1] * hah1_vel;

#if DEBUG_SHAKE
        std::cout << "n=" << numIterate << " hah1_pv=" << hah1_pv << " tolerance=" << tolerance << std::endl;
#endif

        if (std::abs(hah1_pv) > tolerance){
          double k = hah1_pv * massHAinvH1invHAH1posnorm2inv[n1];
          ha_vel += k * massHAinv * hah1_pos[n1];
          h1_vel[n1] -= k * massH1inv[n1] * hah1_pos[n1];
          is_modified = true;
        }
      }
    }
    getvelocity(particle_array,ha_index) = ha_vel;
    for(int n1=0; n1<nh1; n1++){
      int h1_index = it->second.h1[n1];
      getvelocity(particle_array,h1_index) = h1_vel[n1];
    }
  }
}
template
void Settle::rattle_velocity(CombinedParticleArray& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance);
template
void Settle::rattle_velocity(ParticleArray& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance);

template<class PA>
void Settle::rattle_velocity(PA& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     ForceArray& rattle_force)
{
  double Inv2dt = 2.0/dt;
  double tolerance_per_dt = shake_tolerance / dt;

  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int ha_index = it->first;
    double massHAinv = getinvmass(particle_array,ha_index);
    Velocity ha_vel(getvelocity(particle_array,ha_index));

    int nh1 = it->second.nh1;

    double tolerance2HAH1[MAX_HBOND];
    double massH1inv[MAX_HBOND];
    Position hah1_pos[MAX_HBOND];
    Velocity h1_vel[MAX_HBOND];
    double massHAinvH1invHAH1posnorm2inv[MAX_HBOND];

    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      int h1_index = it->second.h1[n1];

      double lengthHAH1 = 1.0;    // default
      bool get_lengthHAH1 = false;
      int bondtype = it->second.bondtype[n1];
      if(bondtype==-1){
	std::cout <<  "Waning rattle_position : search bondtype for index(" <<ha_index << "," << h1_index << ")" << std::endl;
	for(size_t s=0; !get_lengthHAH1 && s<bondlistarray.size(); s++){
	  for(size_t sb=0; !get_lengthHAH1 && sb<bondlistarray[s].BondArray.size(); sb++){
	    CovalentBondInfo::Bond bond = bondlistarray[s].BondArray[sb];

	    if(!get_lengthHAH1 &&
	       ((getatomid(particle_array,ha_index) == bond.id_of_atom[0] && getatomid(particle_array,h1_index) == bond.id_of_atom[1])
		|| (getatomid(particle_array,ha_index) == bond.id_of_atom[1] && getatomid(particle_array,h1_index) == bond.id_of_atom[0]))){
	      //	      lengthHAH1 = param_list->bond[bond.typeofbond].equilibrium_length;
	      bondtype = bond.typeofbond;
	      get_lengthHAH1 = true;
	      break;
	    }
	  }
	}
      }else{
	get_lengthHAH1 = true;
      }
      if(!get_lengthHAH1){
        std::cout << "ERROR: in Settle::rattle_velocity : not hit bond" << std::endl;
        continue;
      }else{
	lengthHAH1 = param_list->bond[bondtype].equilibrium_length;
      }

#if DEBUG_SHAKE
      std::cout << "in Settle::rattle_velocity lengthHAH1=" << lengthHAH1 << " ha_index=" << ha_index << " h1_index=" << h1_index << std::endl;
#endif

      double length2HAH1 = lengthHAH1 * lengthHAH1;
      tolerance2HAH1[n1] = tolerance_per_dt * length2HAH1;
      massH1inv[n1] = getinvmass(particle_array,h1_index);

      hah1_pos[n1] = (getpos(particle_array,h1_index) - getpos(particle_array,ha_index));
      h1_vel[n1] = getvelocity(particle_array,h1_index);
      double massHAinvH1inv = massHAinv + massH1inv[n1];
      massHAinvH1invHAH1posnorm2inv[n1] = 1.0/(massHAinvH1inv*hah1_pos[n1].norm2());
    }

    bool is_modified = true;
    for (int numIterate=0; is_modified && numIterate< shake_max_iterate; ++numIterate) {
      is_modified = false;

      for(int n1=0; n1<nh1; n1++){
        Velocity hah1_vel(h1_vel[n1] - ha_vel);
        double hah1_pv = hah1_pos[n1] * hah1_vel;

#if DEBUG_SHAKE
        std::cout << "n=" << numIterate << " hah1_pv=" << hah1_pv << " tolerance=" << tolerance << std::endl;
#endif

        if (std::abs(hah1_pv) > tolerance){
          double k = hah1_pv * massHAinvH1invHAH1posnorm2inv[n1];
          ha_vel += k * massHAinv * hah1_pos[n1];
          h1_vel[n1] -= k * massH1inv[n1] * hah1_pos[n1];
          is_modified = true;
        }
      }
    }
    rattle_force[ha_index] = (ha_vel - getvelocity(particle_array,ha_index))*Inv2dt*getmass(particle_array,ha_index);
    getvelocity(particle_array,ha_index) = ha_vel;
    for(int n1=0; n1<nh1; n1++){
      int h1_index = it->second.h1[n1];
      rattle_force[h1_index] = (h1_vel[n1] - getvelocity(particle_array,h1_index))*Inv2dt*getmass(particle_array,h1_index);
      getvelocity(particle_array,h1_index) = h1_vel[n1];
    }
  }
}
template
void Settle::rattle_velocity(CombinedParticleArray& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     ForceArray& rattle_force);
template
void Settle::rattle_velocity(ParticleArray& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     ForceArray& rattle_force);

template<class PA>
void Settle::rattle_velocity(PA& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     double& rattle_virial)
{
  double Inv2dt = 2.0/dt;
  double tolerance_per_dt = shake_tolerance / dt;
  double virial=rattle_virial;

  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int ha_index = it->first;
    double massHAinv = getinvmass(particle_array,ha_index);
    Velocity ha_vel(getvelocity(particle_array,ha_index));

    int nh1 = it->second.nh1;

    double tolerance2HAH1[MAX_HBOND];
    double massH1inv[MAX_HBOND];
    Position hah1_pos[MAX_HBOND];
    Velocity h1_vel[MAX_HBOND];
    Velocity hah1_deltavel[MAX_HBOND];
    double massHAinvH1invHAH1posnorm2inv[MAX_HBOND];

    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      int h1_index = it->second.h1[n1];
      hah1_deltavel[n1] = 0.0;

      double lengthHAH1 = 1.0;    // default
      bool get_lengthHAH1 = false;
      int bondtype = it->second.bondtype[n1];
      if(bondtype==-1){
	std::cout <<  "Waning rattle_position : search bondtype for index(" <<ha_index << "," << h1_index << ")" << std::endl;
	for(size_t s=0; !get_lengthHAH1 && s<bondlistarray.size(); s++){
	  for(size_t sb=0; !get_lengthHAH1 && sb<bondlistarray[s].BondArray.size(); sb++){
	    CovalentBondInfo::Bond bond = bondlistarray[s].BondArray[sb];

	    if(!get_lengthHAH1 &&
	       ((getatomid(particle_array,ha_index) == bond.id_of_atom[0] && getatomid(particle_array,h1_index) == bond.id_of_atom[1])
		|| (getatomid(particle_array,ha_index) == bond.id_of_atom[1] && getatomid(particle_array,h1_index) == bond.id_of_atom[0]))){
	      //	      lengthHAH1 = param_list->bond[bond.typeofbond].equilibrium_length;
	      bondtype = bond.typeofbond;
	      get_lengthHAH1 = true;
	      break;
	    }
	  }
	}
      }else{
	get_lengthHAH1 = true;
      }
      if(!get_lengthHAH1){
        std::cout << "ERROR: in Settle::rattle_velocity : not hit bond" << std::endl;
        continue;
      }else{
	lengthHAH1 = param_list->bond[bondtype].equilibrium_length;
      }

#if DEBUG_SHAKE
      std::cout << "in Settle::rattle_velocity lengthHAH1=" << lengthHAH1 << " ha_index=" << ha_index << " h1_index=" << h1_index << std::endl;
#endif

      double length2HAH1 = lengthHAH1 * lengthHAH1;
      tolerance2HAH1[n1] = tolerance_per_dt * length2HAH1;
      massH1inv[n1] = getinvmass(particle_array,h1_index);

      hah1_pos[n1] = (getpos(particle_array,h1_index) - getpos(particle_array,ha_index));
      h1_vel[n1] = getvelocity(particle_array,h1_index);
      double massHAinvH1inv = massHAinv + massH1inv[n1];
      massHAinvH1invHAH1posnorm2inv[n1] = 1.0/(massHAinvH1inv*hah1_pos[n1].norm2());
    }

    bool is_modified = true;
    for (int numIterate=0; is_modified && numIterate< shake_max_iterate; ++numIterate) {
      is_modified = false;

      for(int n1=0; n1<nh1; n1++){
        Velocity hah1_vel(h1_vel[n1] - ha_vel);
        double hah1_pv = hah1_pos[n1] * hah1_vel;

#if DEBUG_SHAKE
        std::cout << "n=" << numIterate << " hah1_pv=" << hah1_pv << " tolerance=" << tolerance << std::endl;
#endif

        if (std::abs(hah1_pv) > tolerance){
          double k = hah1_pv * massHAinvH1invHAH1posnorm2inv[n1];
          ha_vel += k * massHAinv * hah1_pos[n1];
	  hah1_deltavel[n1] += k * hah1_pos[n1];
          h1_vel[n1] -= k * massH1inv[n1] * hah1_pos[n1];
          is_modified = true;
        }
      }
    }
    /*
    rattle_force[ha_index] = (ha_vel - getvelocity(particle_array,ha_index))*Inv2dt*getmass(particle_array,ha_index);
    getvelocity(particle_array,ha_index) = ha_vel;
    for(int n1=0; n1<nh1; n1++){
      int h1_index = it->second.h1[n1];
      rattle_force[h1_index] = (h1_vel[n1] - getvelocity(particle_array,h1_index))*Inv2dt*getmass(particle_array,h1_index);
      getvelocity(particle_array,h1_index) = h1_vel[n1];
    }
    */
    getvelocity(particle_array,ha_index) = ha_vel;
    for(int n1=0; n1<nh1; n1++){
      getvelocity(particle_array,it->second.h1[n1]) = h1_vel[n1];
      virial -= (hah1_deltavel[n1]*Inv2dt*hah1_pos[n1]);
    }
  }
  rattle_virial = virial;
}
template
void Settle::rattle_velocity(CombinedParticleArray& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     double& rattle_virial);
template
void Settle::rattle_velocity(ParticleArray& particle_array, 
                             const ShakeList& shakelist, 
                             const CovalentBondParameterList* param_list, 
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             double dt,
                             int shake_max_iterate, double shake_tolerance,
			     double& rattle_virial);
