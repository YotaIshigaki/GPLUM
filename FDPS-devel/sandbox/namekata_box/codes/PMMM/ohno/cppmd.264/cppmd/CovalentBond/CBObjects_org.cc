
#include "CBInterface.h"
#include "CBObjects.h"
#include "CBInterfaceImpl.h"

using namespace std;
using namespace CBModule;


namespace CBModule{
template <> CBInterface::ParticleLocation CBInterface::getParticleLocation(Context pContext,
                                                const AtomID ai) {
#ifdef HAVE_MAP_AT
  return particleMap.at(ai);
#else
  ParticleMap::iterator pf = particleMap.find(ai);
  if(pf!=particleMap.end()){
    return pf->second;
  }else{
    std::cout << " not found atomid " << ai << std::endl;
    return NULL;
  }
#endif
}

template <> CBInterface::ParticleIndex CBInterface::getParticleLocation(Context pContext,
                                                const AtomID ai) {
#ifdef HAVE_MAP_AT
    return particleIndexMap.at(ai);
#else
    ParticleIndexMap::iterator pf = particleIndexMap.find(ai);
    if(pf!=particleIndexMap.end()){
      return pf->second;
    }else{
      std::cout << " not found atomid " << ai << std::endl;
      ParticleIndex p0 = {0,0};
      return p0;
    }
#endif
}

template<>
void CBInterface::calcBond<CBInterface::ParticleIndex>(const BondArray& bond, double& energy, double& shortenergy, double& shortvirial) 
{
  //    printf("old calcBond\n");
  CBContext context(bond, energy, shortenergy, shortvirial);
  cbObjects.calcCBBond<CBInterface::ParticleIndex>(&context);
 }

template<>
void CBInterface::calcAngle<CBInterface::ParticleIndex>(const AngleArray& angle, double& energy, double& shortenergy, double& shortvirial) {
  //    printf("old calcAngle\n");
  CBContext context(angle, energy, shortenergy, shortvirial);
  cbObjects.calcCBAngle<CBInterface::ParticleIndex>(&context);
}


template<>
void CBInterface::calcTorsion<CBInterface::ParticleIndex>(const TorsionArray& torsion, double& energy, double& shortenergy, double& shortvirial) {
  //  printf("old calcTorsion\n");
  CBContext context(torsion, energy, shortenergy, shortvirial);
  cbObjects.calcCBTorsion<CBInterface::ParticleIndex>(&context);
}

template<>
void CBInterface::calcImproper<CBInterface::ParticleIndex>(const ImproperArray& improper, double& energy, double& shortenergy, double& shortvirial) {
  //  printf("old calcImproper\n");
  CBContext context(improper, energy, shortenergy, shortvirial);
  cbObjects.calcCBImproper<CBInterface::ParticleIndex>(&context);
}
template<>
bool CBInterface::isCalc14Interaction(Context pContext, const Torsion& p)
{
  return p.calc14interaction;
}

template<>
bool CBInterface::isCalc14Interaction(Context pContext, const Improper& p)
{
  return false;
}

}

template<typename PR>
void CBObjects::calcCBBond(CBInterface::Context pContext)
{
  Position d;
  Force force;
  double r,p,dp;
  // DEBUG LJC cancel  int cc=0;
  //  std::cout << "cbInterface at CBObjects::calcCBBond" << cbInterface << std::endl;
  /*
  if(cancel_short){
    std::cout << "CBObjects::calcCBBond with cancel_short" << std::endl;
  }else{
    std::cout << "CBObjects::calcCBBond without cancel_short" << std::endl;
  }
  */
  for (CBInterface::CBBondIterator it =
       cbInterface->begin<CBInterface::Bond>(pContext);
       it != cbInterface->end<CBInterface::Bond>(pContext);++it) {
    const CBInterface::Bond& pBond = *it;
    PR p0 =
      cbInterface->getParticleLocation<PR>(pContext, pBond.id_of_atom[0]);
    PR p1 =
      cbInterface->getParticleLocation<PR>(pContext, pBond.id_of_atom[1]);
    cbInterface->calcDistance(pContext, p0, p1, d);
    if (!cbInterface->isFixed(pBond)) {
      r = sqrt(d.norm2());
      if(r>9.0){
        std::cout << r << " is large bond(Bond) length (" << pBond.id_of_atom[0] << "," << pBond.id_of_atom[1] << ")" << std::endl;
      }
      cbInterface->calcInteraction(pBond, r, p, dp);
      cbInterface->addPotentialEnergy<CBInterface::Bond>(pContext, p);
      force = static_cast<Force>(d*dp);
      cbInterface->addForce(pContext, p0, force);
      cbInterface->subForce(pContext, p1, force);
      /*
      if(cancel_short_virial){
	cbInterface->addshortVirial(pContext,-r*r*dp);
      }
      */
    }
    //! cancel 1-2 LJC
    if(cancel_short){
      CBInterface::ForceLocation f0 =
          cbInterface->getForceLocation(pContext, pBond, 0);
      CBInterface::ForceLocation f1 =
          cbInterface->getForceLocation(pContext, pBond, 1);
      subtractCBParticleForce<PR>(pContext, p0, p1, f0, f1, d);
      /* DEBUG LJC cancel
      cbInterface->cc.push_back(cc);
      cc++;
      cbInterface->ct.push_back(0);
      */
    }
  }
}
template
void CBObjects::calcCBBond<CBInterface::ParticleLocation>(CBInterface::Context pContext);
template
void CBObjects::calcCBBond<CBInterface::ParticleIndex>(CBInterface::Context pContext);

template <typename PR>
void CBObjects::calcCBAngle(CBInterface::Context pContext)
{
  Position d01,d21,d02;
  Force force01,force21;
  double ir01,ir21,p,dp;
  double theta;
// DEBUG LJC cancel    int cc = 0;

  for (CBInterface::CBAngleIterator it =
       cbInterface->begin<CBInterface::Angle>(pContext);
       it != cbInterface->end<CBInterface::Angle>(pContext);++it) {
    const CBInterface::Angle& pAngle = *it;
    PR p0 =
      cbInterface->getParticleLocation<PR>(pContext, pAngle.id_of_atom[0]);
    PR p1 =
      cbInterface->getParticleLocation<PR>(pContext, pAngle.id_of_atom[1]);
    PR p2 =
      cbInterface->getParticleLocation<PR>(pContext, pAngle.id_of_atom[2]);
    cbInterface->calcDistance(pContext, p0, p1, d01);
    cbInterface->calcDistance(pContext, p2, p1, d21);
    if(d01.norm2()>81.0){
      std::cout << d01.norm2() << " is large bond(Angle d01) length (" << pAngle.id_of_atom[0] << "," << pAngle.id_of_atom[1] << ")" << std::endl;
    }
    if(d01.norm2()>81.0){
      std::cout << d21.norm2() << " is large bond(Angle d21) length (" << pAngle.id_of_atom[2] << "," << pAngle.id_of_atom[1] << ")" << std::endl;
    }
    ir01 = 1.0/sqrt(d01.norm2());
    ir21 = 1.0/sqrt(d21.norm2());
    theta = acos((d01*d21)*ir01*ir21);
    cbInterface->calcInteraction(pAngle, theta, p, dp);
    cbInterface->addPotentialEnergy<CBInterface::Angle>(pContext, p);
    double f = dp/sin(theta);
    double cos_theta = cos(theta);
    force01 = f*(d21*ir21 - d01*cos_theta*ir01)*ir01;
    force21 = f*(d01*ir01 - d21*cos_theta*ir21)*ir21;
    cbInterface->addForce(pContext, p0, force01);
    cbInterface->addForce(pContext, p2, force21);
    cbInterface->subForce(pContext, p1, force01+force21);
    if(cancel_short){
      CBInterface::ForceLocation f0 =
          cbInterface->getForceLocation(pContext, pAngle, 0);
      CBInterface::ForceLocation f2 =
          cbInterface->getForceLocation(pContext, pAngle, 2);
      cbInterface->calcDistance(pContext, p0, p2, d02);
      //! cancel 1-3 LJC
      subtractCBParticleForce(pContext, p0, p2, f0, f2, d02);
      /* DEBUG LJC cancel  
      cbInterface->cc.push_back(cc);
      cc++;
      cbInterface->ct.push_back(1);
      */
    }
  }
}
template
void CBObjects::calcCBAngle<CBInterface::ParticleLocation>(CBInterface::Context pContext);
template
void CBObjects::calcCBAngle<CBInterface::ParticleIndex>(CBInterface::Context pContext);

template <typename PR>
void CBObjects::calcCBTorsion(CBInterface::Context pContext)
{
  calcCBTorsionImproper<PR, CBInterface::Torsion,
                        CBInterface::CBTorsionIterator>(pContext);
}
template
void CBObjects::calcCBTorsion<CBInterface::ParticleLocation>(CBInterface::Context pContext);
template
void CBObjects::calcCBTorsion<CBInterface::ParticleIndex>(CBInterface::Context pContext);

template <typename PR>
void CBObjects::calcCBImproper(CBInterface::Context pContext)
{
  calcCBTorsionImproper<PR, CBInterface::Improper,
                        CBInterface::CBImproperIterator>(pContext);
}
template
void CBObjects::calcCBImproper<CBInterface::ParticleLocation>(CBInterface::Context pContext);
template
void CBObjects::calcCBImproper<CBInterface::ParticleIndex>(CBInterface::Context pContext);



template <typename PR, typename T, typename ITERATOR>
void CBObjects::calcCBTorsionImproper(CBInterface::Context pContext)
{
  Position rkj,rkl,rik,rij,rjl;
  Position nv1,nv2,nv;
  double nv1INorm,nv2INorm;
  Force force0,force1;
  double p,dp;
  // DEBUG LJC cancel  int cc=0;

  for (ITERATOR it = cbInterface->begin<T>(pContext);
       it != cbInterface->end<T>(pContext);++it) {
    const T& pTorsion = *it;
    PR pi =
      cbInterface->getParticleLocation<PR>(pContext, pTorsion.id_of_atom[0]);
    PR pj =
      cbInterface->getParticleLocation<PR>(pContext, pTorsion.id_of_atom[1]);
    PR pk =
      cbInterface->getParticleLocation<PR>(pContext, pTorsion.id_of_atom[2]);
    PR pl =
      cbInterface->getParticleLocation<PR>(pContext, pTorsion.id_of_atom[3]);
    cbInterface->calcDistance(pContext, pk, pj, rkj);
    cbInterface->calcDistance(pContext, pk, pl, rkl);
    cbInterface->calcDistance(pContext, pi, pk, rik);
    cbInterface->calcDistance(pContext, pi, pj, rij);
    cbInterface->calcDistance(pContext, pj, pl, rjl);
    int flag=0;
    if(rkj.norm2()>81.0){
      std::cout << rkj.norm2() << " is large bond(Torsion rkj) length (" << pTorsion.id_of_atom[2] << "," << pTorsion.id_of_atom[1] << ")" << std::endl;
      flag++;
    }
    if(rkl.norm2()>81.0){
      std::cout << rkl.norm2() << " is large bond(Torsion rkl) length (" << pTorsion.id_of_atom[2] << "," << pTorsion.id_of_atom[3] << ")" << std::endl;
      flag++;
    }
    if(rik.norm2()>81.0){
      std::cout << rik.norm2() << " is large bond(Torsion rik) length (" << pTorsion.id_of_atom[0] << "," << pTorsion.id_of_atom[2] << ")" << std::endl;
      flag++;
    }
    if(rij.norm2()>81.0){
      std::cout << rij.norm2() << " is large bond(Torsion rij) length (" << pTorsion.id_of_atom[0] << "," << pTorsion.id_of_atom[1] << ")" << std::endl;
      flag++;
    }
    if(rjl.norm2()>81.0){
      std::cout << rjl.norm2() << " is large bond(Torsion rjl) length (" << pTorsion.id_of_atom[1] << "," << pTorsion.id_of_atom[3] << ")" << std::endl;
      flag++;
    }
    if(flag>0){
      std::cout << " at Torsion/Improper (" << pTorsion.id_of_atom[0] << "," << pTorsion.id_of_atom[1] << "," << pTorsion.id_of_atom[2] << "," << pTorsion.id_of_atom[3] << ")" << std::endl;
    }
    nv1 = (rij % rkj);
    nv2 = (rkj % rkl);
    nv1INorm = 1.0/nv1.norm();
    nv2INorm = 1.0/nv2.norm();
    nv = nv1 % nv2;
    double sign = nv*rkj;
    double dot = nv1*nv2*nv1INorm*nv2INorm;
    if(dot > 1.0){
      dot = 1.0;
    } else if(dot < -1.0){
      dot = -1.0;
    }
    double theta  = acos(dot)*(sign > 0 ? 1:-1);
    double cos_theta = cos(theta);
    for(int i = 0;i < cbInterface->getParameterNum(pTorsion);++i){
      cbInterface->calcInteraction(pTorsion, i, theta, p, dp);
      cbInterface->addPotentialEnergy<T>(pContext, p);
      force0 = (nv2*nv2INorm - cos_theta*nv1*nv1INorm)*nv1INorm;
      force1 = (nv1*nv1INorm - cos_theta*nv2*nv2INorm)*nv2INorm;
      cbInterface->addForce(pContext, pi,dp*force0%rkj);
      cbInterface->addForce(pContext, pj,dp*(force0%rik - force1%rkl));
      cbInterface->addForce(pContext, pk,dp*(force1%rjl - force0%rij));
      cbInterface->addForce(pContext, pl,dp*force1%rkj);
    }
    if(cancel_short){
      CBInterface::ForceLocation fi =
          cbInterface->getForceLocation(pContext, pTorsion, 0);
      CBInterface::ForceLocation fl =
          cbInterface->getForceLocation(pContext, pTorsion, 3);
    //! cancel 1-4 LJC
      if(cbInterface->isCalc14Interaction(pContext, pTorsion)){
        Position dil;
        cbInterface->calcDistance(pContext, pi, pl, dil);
        calcCB14Interaction(pContext, pi, pl, fi, fl, dil);
	/* DEBUG LJC cancel
	cbInterface->cc.push_back(cc);
	cc++;
	cbInterface->ct.push_back(2);
	*/
        subtractCBParticleForce(pContext, pi, pl, fi, fl, dil);
	/* DEBUG LJC cancel
	cbInterface->cc.push_back(cc);
	cc++;
	cbInterface->ct.push_back(2);
	*/
      }
    }
  }
}

template<typename PR>
void CBObjects::calcCB14Interaction(CBInterface::Context pContext,
                                    PR pi,
                                    PR pl, 
                                    CBInterface::ForceLocation fi, 
                                    CBInterface::ForceLocation fl, 
				    const Position& d)
{
  SpaceVector<double> force;
  double r2;
  double shortRangePotential;
  double longRangePotential;

  r2 = d.norm2();
  if (cbInterface->inCutoffSphere(pContext, r2)) {
    cbInterface->calcInteraction(pContext,
                                 pi, pl,
		                 scnbInv, sceeInv,
		                 d, r2,
		                 shortRangePotential,
		                 longRangePotential,
		                 force); 
    /* DEBUG LJC cancel
    tmpshort -= shortRangePotential;
    cbInterface->c0.push_back(pi->atomid);
    cbInterface->c1.push_back(pl->atomid);
    cbInterface->ce.push_back(-shortRangePotential);
    */
    cbInterface->addShortPotentialEnergy(pContext, shortRangePotential);
    cbInterface->addLongPotentialEnergy(pContext, longRangePotential);
    /* no separate short force, for Single Time Step
    cbInterface->addForce(pContext, pi,force);
    cbInterface->subForce(pContext, pl,force);
    */
    cbInterface->addshortForce(pContext, fi,force);
    cbInterface->subshortForce(pContext, fl,force);
  }
  /* DEBUG LJC cancel
  else{
    cbInterface->c0_off.push_back(pi->atomid);
    cbInterface->c1_off.push_back(pl->atomid);
  }
  */
}

//! cancel LJC force
template<typename PR>
void CBObjects::subtractCBParticleForce(CBInterface::Context pContext,
                                        PR pi,
                                        PR pj,
                                        CBInterface::ForceLocation fi,
                                        CBInterface::ForceLocation fj,
					const Position& d)
{
  SpaceVector<double> force;
  double r2;
  double shortRangePotential;
  double longRangePotential;
  double virial;

  r2 = d.norm2();
  if (cbInterface->inCutoffSphere(pContext, r2)) {
    if(cancel_short_virial){
      cbInterface->calcInteraction(pContext,
				   pi, pj,
				   d,r2,
				   shortRangePotential,
				   longRangePotential,
				   force, virial);
    }else{
      cbInterface->calcInteraction(pContext,
				   pi, pj,
				   d,r2,
				   shortRangePotential,
				   longRangePotential,
				   force);
    }
    /* DEBUG LJC cancel
    tmpshort += shortRangePotential;
    cbInterface->c0.push_back(pi->atomid);
    cbInterface->c1.push_back(pj->atomid);
    cbInterface->ce.push_back(shortRangePotential);
    */
    cbInterface->addShortPotentialEnergy(pContext, -shortRangePotential);
    cbInterface->addLongPotentialEnergy(pContext, -longRangePotential);
    /* no separate short force, for Single Time Step
    cbInterface->subForce(pContext, pi, force);
    cbInterface->addForce(pContext, pj, force);
    */
    cbInterface->subshortForce(pContext, fi, force);
    cbInterface->addshortForce(pContext, fj, force);
    if(cancel_short_virial){
      cbInterface->addshortVirial(pContext,-virial);
    }
  }
  /* DEBUG LJC cancel
  else{
    cbInterface->c0_off.push_back(pi->atomid);
    cbInterface->c1_off.push_back(pj->atomid);
  }
  */
}

/*
  cancel short

  Covalent Bond force includes no-bonded force(Coulomb and LJ).
  The no-bonded force must be excluded or canceled.

  potential form
                  Original-Coulomb    Ewald
  no-bonded pair    q_i q_j / r_ij     q_i q_j erfc(beta r_ijn) / r_ijn
  cancel          - q_i q_j / r_ij    -q_i q_j /r_ijn
  exculded            0               -q_i q_j erf(beta r_ijn) / r_ijn

  where   n is periodic shift
  erfc(X) = 1 - erf(X)
  no-bonded Ewald  q_i q_j erfc(beta r_ijn) / r_ijn
                 = q_i q_j / r_ijn - q_i q_j erf(beta r_ijn) / r_ijn
  first term is icluded CB, must be canceled when bonded pair is not excluded
  secondary term caused from Ewald sum, must be calculated for bonded pair when exluded

  if(nb_correction){
    if(excluded){
      if(ewald){
        E -= qi*qj*erf(beta*rij)/rij
	Fij[axis] += qi*qj/rij^3 * rij[axis] * (2*beta*rij/sqrt(M_PI)*exp(-(beta
*rij)^2) - erf(beta*rij))
// or Fij[axis] += ewaldforce(rij) - qi*qj/rij^3 * rij[axis]
        Virial += qi*qj * (2*beta/sqrt(M_PI)*exp(-(beta*rij)^2) - erf(beta*rij)/rij)
      }
    }else{   // cancel
      E -= qi*qj/rij
      fij[axis] -= qi*qj/rij^3 * rij[axis]
      Virial -= qi*qj/rij
    }
  }
 */
