
#include "CBInterface.h"
#include "CBObjects.h"
#include "CBInterfaceImpl.h"

using namespace std;
using namespace CBModule;


namespace CBModule{
CBInterface::ParticleIndex CBInterface::getParticleLocation(CBInterface::ParticleIndexMap& pmap, const int i)
{
#ifdef HAVE_MAP_AT
  return pmap.at(i);
#else
  CBInterface::ParticleIndexMap::iterator pf = pmap.find(i);
  if(pf!=pmap.end()){
    return pf->second;
  }else{
    std::cout << " not found atomid " << i << std::endl;
    ParticleIndex p0;
    noneParticleLocation(p0);
    return p0;
  }
#endif
}

template<typename T>
void CBInterface::convert_torsion_atomid_to_index(const CBInterface::ParticleIndexMap& pmap,
						  const std::vector<T>& torsion, 
						  std::vector< std::vector<CBInterface::ParticleIndex> >& torsion_index)
{
  torsion_index.resize(torsion.size());
  for(int i=0; i<torsion.size(); i++){
    torsion_index[i].resize(4);
    for(int j=0;j<4;j++){
      int aid = torsion[i].id_of_atom[j];
      CBInterface::ParticleIndexMap::const_iterator pf = pmap.find(aid);
      if(pf!=pmap.end()){
	torsion_index[i][j] = pf->second;
      }else{
	std::cout << " not found atomid " << aid << std::endl;
	noneParticleLocation(torsion_index[i][j]);
      }
    }
  }
}

void CBInterface::convert_torsion_improper_atomid_to_index(const TorsionArray& torsion, const ImproperArray& improper)
{
  convert_torsion_atomid_to_index(particleIndexMap,torsion,torsion_index);
  convert_torsion_atomid_to_index(particleIndexMap,improper,improper_index);
}

//! cancel LJC force
template<>
void CBInterface::subtractCBParticleForce(double& energy, double& shortenergy, double& shortvirial,
					  CBInterface::ParticleIndex pi,
					  CBInterface::ParticleIndex pj,
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
  if (inCutoffSphere(r2)) {
    if(cancel_short_virial){
      calcInteraction(pi, pj,
		      d,r2,
		      shortRangePotential,
		      longRangePotential,
		      force, virial);
    }else{
      calcInteraction(pi, pj,
		      d,r2,
		      shortRangePotential,
		      longRangePotential,
		      force);
    }
    shortenergy -= shortRangePotential;
    // energy -= longRangePotential;
    (*fi) -= force;
    (*fj) += force;
    if(cancel_short_virial){
      shortvirial -= virial;
    }
  }
}

void CBInterface::clear_tmpforce()
{
#pragma omp parallel
  {
    int t = 0;
#ifdef _OPENMP
    t = omp_get_thread_num();
#endif
    //   printf("clear tmpforce[%d]\n",t);
    double (*tf0)[3] = (double (*)[3])(&(tmpforce[t][0][0].x));
    int nt = tmpforce[t][0].size();
#pragma loop norecurrence
#pragma loop noalias
    for(size_t i=0;i<nt;i++){
      tf0[i][0] = 0.0; //      tmpforce[t][0][i].x = 0.0;
      tf0[i][1] = 0.0; //      tmpforce[t][0][i].y = 0.0;
      tf0[i][2] = 0.0; //      tmpforce[t][0][i].z = 0.0;
    }
    double (*tf1)[3] = (double (*)[3])(&(tmpforce[t][1][0].x));
    nt = tmpforce[t][1].size();
#pragma loop norecurrence
    for(size_t i=0;i<nt;i++){
      tf1[i][0] = 0.0;
      tf1[i][1] = 0.0;
      tf1[i][2] = 0.0;
    }
  }
}

void CBInterface::merge_tmpforce()
{
  //  printf("CB reduct %d thread tmpforce\n",num_threads);
  int i;
#ifdef K_SIMD
  double (*forceaca)[3] = (double (*)[3] )(&((*forcea[0])[0].x));
  int nf = tmpforce[0][0].size();
  for(int t=0;t<num_threads;t++){
    const double (*tf)[3] = (double (*)[3])(&(tmpforce[t][0][0].x));
#pragma omp parallel for private(i)
#pragma loop norecurrence
    for(i=0;i<nf;i++){
      forceaca[i][0] += tf[i][0];
      forceaca[i][1] += tf[i][1];
      forceaca[i][2] += tf[i][2];
    }
  }
  double (*forcea1)[3]  = (double (*)[3] )(&((*forcea[1])[0].x));
  nf = tmpforce[0][1].size();
  for(int t=0;t<num_threads;t++){
    const double (*tf)[3] = (double (*)[3])(&(tmpforce[t][1][0].x));
#pragma omp parallel for private(i)
    for(i=0;i<nf;i++){
      forcea1[i][0] += tf[i][0];
      forcea1[i][1] += tf[i][1];
      forcea1[i][2] += tf[i][2];
    }
  }
#else
  for(int t=0;t<num_threads;t++){
#pragma omp parallel for
#pragma loop norecurrence
    for(i=0;i<tmpforce[0][0].size();i++){
      (*forcea[0])[i].x += tmpforce[t][0][i].x;
      (*forcea[0])[i].y += tmpforce[t][0][i].y;
      (*forcea[0])[i].z += tmpforce[t][0][i].z;
    }
  }
#pragma omp parallel for
  for(i=0;i<tmpforce[0][1].size();i++){
    for(int t=0;t<num_threads;t++){
      (*forcea[1])[i] += tmpforce[t][1][i];
    }
  }
#endif
}

template<>
void CBInterface::calcBond<CBInterface::ParticleIndex>(const BondArray& bond, double& energy, double& shortenergy, double& shortvirial) 
{
  double *cbenergy = new double[num_threads];
  //    printf("flat calcBond\n");
#define CBOMP_BOND
# ifdef CBOMP_BOND
  //    clear_tmpforce();
# endif


# ifdef CBOMP_BOND
#pragma omp parallel 
# endif
  {
    int i;
    int thread_id=0;
# ifdef CBOMP_BOND
#ifdef _OPENMP
    thread_id=omp_get_thread_num();
#endif
# endif
    cbenergy[thread_id] = 0.0;
# ifdef CBOMP_BOND
#define THREAD_FORCE_POINTER
#  ifdef THREAD_FORCE_POINTER
    double (*tf[2])[3];
    tf[0] = (double (*)[3])(&(tmpforce[thread_id][0][0].x));
    tf[1] = (double (*)[3])(&(tmpforce[thread_id][1][0].x));
#  endif
# endif
# ifdef CBOMP_BOND
#pragma omp for 
# endif
#pragma loop noalias
    for (i=0; i<bond.size();i++) {
      const CBInterface::Bond& pBond = bond[i];
      Position d;
      Force force;
      double r,p,dp;
# ifdef TUNE_CBMAP
      CBInterface::ParticleIndex p0 = bondIndexArray[i].index_of_atom[0];
      CBInterface::ParticleIndex p1 = bondIndexArray[i].index_of_atom[1];
# else  // TUNE_CBMAP
      CBInterface::ParticleIndex p0 = getParticleLocation(particleIndexMap,pBond.id_of_atom[0]);
      CBInterface::ParticleIndex p1 = getParticleLocation(particleIndexMap,pBond.id_of_atom[1]);
# endif  // TUNE_CBMAP
      d = (*poscharge[p0.array])[p0.index].position - (*poscharge[p1.array])[p1.index].position;
      if (!pBond.shake) {
	r = sqrt(d.norm2());
#ifndef NDEBUG
	if(r>9.0){
	  std::cout << r << " is large bond(Bond) length (" << pBond.id_of_atom[0] << "," << pBond.id_of_atom[1] << ")" << std::endl;
	}
#endif
	{
	  const BondParameter& param = pParameterList->bond[pBond.typeofbond]; 
	  double db = r - param.equilibrium_length;
	  p = param.force_constant * db * db;
	  dp = -2.0 * param.force_constant * db / r;
	}
	cbenergy[thread_id] += p;
	force = static_cast<Force>(d*dp);
# ifdef CBOMP_BOND
#  ifdef THREAD_FORCE_POINTER
	tf[p0.array][p0.index][0] += force.x;
	tf[p0.array][p0.index][1] += force.y;
	tf[p0.array][p0.index][2] += force.z;
	tf[p1.array][p1.index][0] -= force.x;
	tf[p1.array][p1.index][1] -= force.y;
	tf[p1.array][p1.index][2] -= force.z;
#  else
	tmpforce[thread_id][p0.array][p0.index] += force;
	tmpforce[thread_id][p1.array][p1.index] -= force;
#  endif
# else
	(*forcea[p0.array])[p0.index] += force;
	(*forcea[p1.array])[p1.index] -= force;
# endif
	/*
	  if(cancel_short_virial){
	  cbInterface->addshortVirial(pContext,-r*r*dp);
	  }
	*/
      }
    }
  }
# ifdef CBOMP_BOND
  //  merge_tmpforce();
  /*
  {
    //    printf("CB reduct %d thread tmpforce\n",num_threads);
    int i;
#pragma omp parallel for
    for(i=0;i<tmpforce[0][0].size();i++){
      for(int t=0;t<num_threads;t++){
	(*forcea[0])[i] += tmpforce[t][0][i];
      }
    }
#pragma omp parallel for
    for(i=0;i<tmpforce[0][1].size();i++){
      for(int t=0;t<num_threads;t++){
	(*forcea[1])[i] += tmpforce[t][1][i];
      }
    }
  }
  */
# endif
  for(int t=0;t<num_threads;t++){
    energy += cbenergy[t];
  }

  if(cancel_short){
    int i;
    for (i=0; i<bond.size();i++) {
      const CBInterface::Bond& pBond = bond[i];
      Position d;
      Force force;
# ifdef TUNE_CBMAP
      CBInterface::ParticleIndex p0 = bondIndexArray[i].index_of_atom[0];
      CBInterface::ParticleIndex p1 = bondIndexArray[i].index_of_atom[1];
# else  // TUNE_CBMAP
      CBInterface::ParticleIndex p0 = getParticleLocation(particleIndexMap,pBond.id_of_atom[0]);
      CBInterface::ParticleIndex p1 = getParticleLocation(particleIndexMap,pBond.id_of_atom[1]);
# endif  // TUNE_CBMAP
      d = (*poscharge[p0.array])[p0.index].position - (*poscharge[p1.array])[p1.index].position;
      //! cancel 1-2 LJC
      {
	CBInterface::ForceLocation f0 = getForceLocation(pBond.id_of_atom[0]);
	CBInterface::ForceLocation f1 = getForceLocation(pBond.id_of_atom[1]);
	{
	  SpaceVector<double> force;
	  double r2;
	  double shortRangePotential;
	  double longRangePotential;
	  double virial;

	  r2 = d.norm2();
	  if (inCutoffSphere(r2)) {
	    if(cancel_short_virial){
	      calcInteraction(p0, p1,
			      d,r2,
			      shortRangePotential,
			      longRangePotential,
			      force, virial);
	    }else{
	      calcInteraction(p0, p1,
			      d,r2,
			      shortRangePotential,
			      longRangePotential,
			      force);
	    }
	    shortenergy -= shortRangePotential;
	    //	  cbInterface->addLongPotentialEnergy(pContext, -longRangePotential);
	    (*f0) -= force;
	    (*f1) += force;
	    if(cancel_short_virial){
	      shortvirial -= virial;
	    }
	  }
	}
      } // cancel short
    }
  }
 }

template<>
void CBInterface::calcAngle<CBInterface::ParticleIndex>(const AngleArray& angle, double& energy, double& shortenergy, double& shortvirial) {

  double *cbenergy = new double[num_threads];
#define CBOMP_ANGLE
# ifdef CBOMP_ANGLE
  //      clear_tmpforce();
# endif

# ifdef CBOMP_ANGLE
#pragma omp parallel 
# endif  
  {
    int i;

    int thread_id=0;
# ifdef CBOMP_ANGLE
#ifdef _OPENMP
    thread_id=omp_get_thread_num();
#endif
# endif
    cbenergy[thread_id] = 0.0;
# ifdef CBOMP_ANGLE
#  ifdef THREAD_FORCE_POINTER
    double (*tf[2])[3];
    tf[0] = (double (*)[3])(&(tmpforce[thread_id][0][0].x));
    tf[1] = (double (*)[3])(&(tmpforce[thread_id][1][0].x));
#  endif
# endif
#pragma loop norecurrence tf
# ifdef CBOMP_ANGLE
#pragma omp for 
# endif
#pragma loop noalias
    for(i=0;i<angle.size();i++){
      const CBInterface::Angle& pAngle = angle[i];
      Position d01,d21,d02;
      Force force01,force21;
      double ir01,ir21,p,dp;
      double theta;
# ifdef TUNE_CBMAP
      CBInterface::ParticleIndex p0 = angleIndexArray[i].index_of_atom[0];
      CBInterface::ParticleIndex p1 = angleIndexArray[i].index_of_atom[1];
      CBInterface::ParticleIndex p2 = angleIndexArray[i].index_of_atom[2];
# else  // TUNE_CBMAP
      CBInterface::ParticleIndex p0 = getParticleLocation(particleIndexMap,pAngle.id_of_atom[0]);
      CBInterface::ParticleIndex p1 = getParticleLocation(particleIndexMap,pAngle.id_of_atom[1]);
      CBInterface::ParticleIndex p2 = getParticleLocation(particleIndexMap,pAngle.id_of_atom[2]);
# endif  // TUNE_CBMAP
      d01 = (*poscharge[p0.array])[p0.index].position - (*poscharge[p1.array])[p1.index].position;
      d21 = (*poscharge[p2.array])[p2.index].position - (*poscharge[p1.array])[p1.index].position;
#ifndef NDEBUG 
      if(d01.norm2()>81.0){
	std::cout << d01.norm2() << " is large bond(Angle d01) length (" << pAngle.id_of_atom[0] << "," << pAngle.id_of_atom[1] << ")" << std::endl;
      }
      if(d21.norm2()>81.0){
	std::cout << d21.norm2() << " is large bond(Angle d21) length (" << pAngle.id_of_atom[2] << "," << pAngle.id_of_atom[1] << ")" << std::endl;
      }
#endif
      ir01 = 1.0/sqrt(d01.norm2());
      ir21 = 1.0/sqrt(d21.norm2());
#ifdef REDUCE_COS
      double cos_theta = (d01*d21)*ir01*ir21;
      theta = acos(cos_theta);
      double recp_sin_theta = 1.0/sqrt(1.0-cos_theta*cos_theta);
#else
      theta = acos((d01*d21)*ir01*ir21);
      double cos_theta = cos(theta);
#endif
      {
	const AngleParameter& param = pParameterList->angle[pAngle.typeofangle];
	double d = theta - param.equilibrium_angle;
	p = param.force_constant * d * d;
	dp = 2.0 * param.force_constant * d;
      }
      cbenergy[thread_id] += p;
#ifdef REDUCE_COS
      double f = dp*recp_sin_theta;
#else
      double f = dp/sin(theta);
#endif
      force01 = f*(d21*ir21 - d01*cos_theta*ir01)*ir01;
      force21 = f*(d01*ir01 - d21*cos_theta*ir21)*ir21;
# ifdef CBOMP_ANGLE
#  ifdef THREAD_FORCE_POINTER
      Force force012 = force01+force21;
      tf[p0.array][p0.index][0] += force01.x;
      tf[p0.array][p0.index][1] += force01.y;
      tf[p0.array][p0.index][2] += force01.z;
      tf[p2.array][p2.index][0] += force21.x;
      tf[p2.array][p2.index][1] += force21.y;
      tf[p2.array][p2.index][2] += force21.z;
      tf[p1.array][p1.index][0] -= force012.x;
      tf[p1.array][p1.index][1] -= force012.y;
      tf[p1.array][p1.index][2] -= force012.z;
#  else
      tmpforce[thread_id][p0.array][p0.index] += force01;
      tmpforce[thread_id][p2.array][p2.index] += force21;
      tmpforce[thread_id][p1.array][p1.index] -= (force01+force21);
#  endif
# else
      (*forcea[p0.array])[p0.index] += force01;
      (*forcea[p2.array])[p2.index] += force21;
      (*forcea[p1.array])[p1.index] -= (force01+force21);
# endif
    }
  }
# ifdef CBOMP_ANGLE

  //  merge_tmpforce();
  /*
  {
    //    printf("CB reduct %d thread tmpforce\n",num_threads);
    int i;
    for(int t=0;t<num_threads;t++){
#pragma omp parallel for
#pragma loop norecurrence
      for(i=0;i<tmpforce[0][0].size();i++){
	(*forcea[0])[i].x += tmpforce[t][0][i].x;
	(*forcea[0])[i].y += tmpforce[t][0][i].y;
	(*forcea[0])[i].z += tmpforce[t][0][i].z;
      }
    }
#pragma omp parallel for
    for(i=0;i<tmpforce[0][1].size();i++){
      for(int t=0;t<num_threads;t++){
	(*forcea[1])[i] += tmpforce[t][1][i];
      }
    }
  }
  */
# endif
  for(int t=0;t<num_threads;t++){
    energy += cbenergy[t];
  }

  if(cancel_short){
    for(int i=0 ;i<angle.size();i++){
      Position d01,d21,d02;
      Force force01,force21;
      const CBInterface::Angle& pAngle = angle[i];
# ifdef TUNE_CBMAP
      CBInterface::ParticleIndex p0 = angleIndexArray[i].index_of_atom[0];
      CBInterface::ParticleIndex p2 = angleIndexArray[i].index_of_atom[2];
# else  // TUNE_CBMAP
      CBInterface::ParticleIndex p0 = getParticleLocation(particleIndexMap,pAngle.id_of_atom[0]);
      CBInterface::ParticleIndex p2 = getParticleLocation(particleIndexMap,pAngle.id_of_atom[2]);
# endif  // TUNE_CBMAP
      {
	CBInterface::ForceLocation f0 = getForceLocation(pAngle.id_of_atom[0]);
	CBInterface::ForceLocation f2 = getForceLocation(pAngle.id_of_atom[2]);
	d02 = (*poscharge[p0.array])[p0.index].position - (*poscharge[p2.array])[p2.index].position;
      //! cancel 1-3 LJC
	subtractCBParticleForce(energy, shortenergy, shortvirial, p0, p2, f0, f2, d02);
      }
    }

  }
}

#ifdef TUNE_CBMAP
template
void CBInterface::calcCBTorsionImproper<CBInterface::Torsion, CBInterface::CBTorsionIterator>(const CBInterface::TorsionArray & torsion, double& energy, double& shortenergy, double& shortvirial, const std::vector<PI4>& indexArray);
#else  // TUNE_CBMAP
template
void CBInterface::calcCBTorsionImproper<CBInterface::Torsion, CBInterface::CBTorsionIterator>(const CBInterface::TorsionArray & torsion, double& energy, double& shortenergy, double& shortvirial);
#endif  // TUNE_CBMAP

#ifdef TUNE_CBMAP
template<>
void CBInterface::calcTorsion<CBInterface::ParticleIndex>(const TorsionArray& torsion, double& energy, double& shortenergy, double& shortvirial, const std::vector<PI4>& torsionIndexArray) {
  //  printf("flat calcTorsion\n");
  calcCBTorsionImproper<Torsion,CBTorsionIterator>(torsion, energy, shortenergy, shortvirial, torsionIndexArray);
}
#else  // TUNE_CBMAP
template<>
void CBInterface::calcTorsion<CBInterface::ParticleIndex>(const TorsionArray& torsion, double& energy, double& shortenergy, double& shortvirial) {
  //  printf("flat calcTorsion\n");
  calcCBTorsionImproper<Torsion,CBTorsionIterator>(torsion, energy, shortenergy, shortvirial);
}
#endif  // TUNE_CBMAP

#ifdef TUNE_CBMAP
template<>
void CBInterface::calcImproper<CBInterface::ParticleIndex>(const ImproperArray& improper, double& energy, double& shortenergy, double& shortvirial, const std::vector<PI4>& improperIndexArray) {
  //  printf("flat calcImproper\n");
  calcCBTorsionImproper<Improper,CBImproperIterator>(improper, energy, shortenergy, shortvirial, improperIndexArray);
}
#else  // TUNE_CBMAP
template<>
void CBInterface::calcImproper<CBInterface::ParticleIndex>(const ImproperArray& improper, double& energy, double& shortenergy, double& shortvirial) {
  //  printf("flat calcImproper\n");
  calcCBTorsionImproper<Improper,CBImproperIterator>(improper, energy, shortenergy, shortvirial);
}
#endif  // TUNE_CBMAP

template<>
void CBInterface::calcInteraction(const Bond& p,
                                  const double r,
                                  double& potential,
                                  double& force)
{
  if (!pParameterList) {
    throw std::runtime_error(errorPos("CovalentBondParameterList is not set"));
  }
  const BondParameter& param = pParameterList->bond[p.typeofbond];
  double d = r - param.equilibrium_length;
  potential = param.force_constant * d * d;
  force = -2.0 * param.force_constant * d / r;
}

template<>
void CBInterface::calcInteraction(const Angle& p,
                                  const double r,
                                  double& potential,
                                  double& force)
{
  if (!pParameterList) {
    throw std::runtime_error(errorPos("CovalentBondParameterList is not set"));
  }
  const AngleParameter& param = pParameterList->angle[p.typeofangle];
  double d = r - param.equilibrium_angle;
  potential = param.force_constant * d * d;
  force = 2.0 * param.force_constant * d;
}


template <typename TP>
void CBInterface::calcInteractionTorsion(const TP& param,
                                         const double theta,
                                         double& potential,
                                         double& force)
{
  potential = param.force_constant *
              (1 + cos(theta*param.periodicity - param.phase));
  force = param.force_constant * param.periodicity * sin(theta*param.periodicity - param.phase)/sin(theta);
}

template<>
void CBInterface::calcInteraction(const Torsion& p,
                                  int iparam,
                                  const double theta,
                                  double& potential,
                                  double& force)
{
  if (!pParameterList) {
    throw std::runtime_error(errorPos("CovalentBondParameterList is not set"));
  }
  const TorsionParameter& param = pParameterList->torsion[p.typeoftorsion];
  calcInteractionTorsion(param, theta, potential, force);
}

template<>
void CBInterface::calcInteraction(const Improper& p,
                                  int iparam,
                                  const double theta,
                                  double& potential,
                                  double& force)
{
  if (!pParameterList) {
    throw std::runtime_error(errorPos("CovalentBondParameterList is not set"));
  }
  const ImproperParameter& param = pParameterList->improper[p.typeofimproper];
  calcInteractionTorsion(param, theta, potential, force);
}

#ifdef BOND_COUNT
static int bc[5] = {0,0,0,0,0};
#endif

template <typename TP>
void CBInterface::calcInteractionTorsion(const TP& param,
                                         const double cos_theta,
					 const double sign,
                                         double& potential,
                                         double& force)
{
  double cos_theta_2 = cos_theta*cos_theta;
  double sin_theta_2 = 1.0-cos_theta*cos_theta;
  double cos_p = cos(param.phase);
  //  double sin_p = sin(param.phase);  = 0,  phase = 0 or pi
  double cos_nt, sin_nt_sin_ti;
  if (param.periodicity == 1) {
    cos_nt = cos_theta;
    sin_nt_sin_ti = 1.0;
#ifdef BOND_COUNT
    bc[1]++;
#endif
  } else if (param.periodicity == 2) {
    cos_nt = cos_theta_2 - sin_theta_2;
    sin_nt_sin_ti = 2.0*cos_theta;
#ifdef BOND_COUNT
    bc[2]++;
#endif
  } else if (param.periodicity == 3) {
    cos_nt = cos_theta*(4.0*cos_theta_2 - 3.0);
    sin_nt_sin_ti = (3.0-4.0*sin_theta_2);
#ifdef BOND_COUNT
    bc[3]++;
#endif
  } else if (param.periodicity == 4) {
    cos_nt = (cos_theta_2 - sin_theta_2)*(cos_theta_2 - sin_theta_2)
      - 4.0*sin_theta_2*cos_theta_2;
    sin_nt_sin_ti = 4.0*cos_theta*(cos_theta_2 - sin_theta_2);
#ifdef BOND_COUNT
    bc[4]++;
#endif
  } else {
    double theta  = acos(cos_theta)*(sign > 0.0 ? 1.0:-1.0);
    cos_nt = cos(param.periodicity*theta);
    sin_nt_sin_ti = sin(param.periodicity*theta)/sin(theta);
#ifdef BOND_COUNT
    bc[5]++;
#endif
  }
  force = param.force_constant * param.periodicity * sin_nt_sin_ti*cos_p; // sin_p = 0
  potential = param.force_constant*(1.0 + cos_nt*cos_p);  // sin_p = 0
}

template<>
void CBInterface::calcInteraction(const Torsion& p,
                                  int iparam,
				  const double cos_theta,
				  const double sign,
                                  double& potential,
                                  double& force)
{
  if (!pParameterList) {
    throw std::runtime_error(errorPos("CovalentBondParameterList is not set"));
  }
  const TorsionParameter& param = pParameterList->torsion[p.typeoftorsion];
  calcInteractionTorsion(param, cos_theta, sign, potential, force);
}

template<>
void CBInterface::calcInteraction(const Improper& p,
                                  int iparam,
				  const double cos_theta,
				  const double sign,
                                  double& potential,
                                  double& force)
{
  if (!pParameterList) {
    throw std::runtime_error(errorPos("CovalentBondParameterList is not set"));
  }
  const ImproperParameter& param = pParameterList->improper[p.typeofimproper];
  calcInteractionTorsion(param, cos_theta, sign, potential, force);
}

template<>
bool CBInterface::isCalc14Interaction(const Torsion& p)
{
  return p.calc14interaction;
}

template<>
bool CBInterface::isCalc14Interaction(const Improper& p)
{
  return false;
}

}

#ifdef TUNE_CBMAP
template <typename T, typename TI>
void CBInterface::calcCBTorsionImproper(const std::vector<T>& torsion, double& energy, double& shortenergy, double& shortvirial, const std::vector<PI4>& indexArray)
#else  // TUNE_CBMAP
template <typename T, typename TI>
void CBInterface::calcCBTorsionImproper(const std::vector<T>& torsion, double& energy, double& shortenergy, double& shortvirial)
#endif  // TUNE_CBMAP
{

  double *cbenergy = new double[num_threads];
#define CBOMP_TORSION
# ifdef CBOMP_TORSION
  //     clear_tmpforce();
# endif

# ifdef CBOMP_TORSION
#pragma omp parallel 
# endif
  {
    int i;
    int thread_id=0;
# ifdef CBOMP_TORSION
#ifdef _OPENMP
    thread_id=omp_get_thread_num();
#endif
# endif
    cbenergy[thread_id] = 0.0;
# ifdef CBOMP_TORSION
#  ifdef THREAD_FORCE_POINTER
    double (*tf[2])[3];
    tf[0] = (double (*)[3])(&(tmpforce[thread_id][0][0].x));
    tf[1] = (double (*)[3])(&(tmpforce[thread_id][1][0].x));
#  endif
# endif
    //    printf("calcCBTorsionImproper thread %d\n", thread_id);
    /*
#pragma loop norecurrence tf
    */
# ifdef CBOMP_TORSION
#pragma omp for
# endif
#pragma loop noalias
    for(i=0; i<torsion.size(); i++){
      Position nv1,nv2,nv;
      double nv1INorm,nv2INorm;
      Force force0,force1;
      double p,dp;

      const T& pTorsion = torsion[i];
# ifdef TUNE_CBMAP
      CBInterface::ParticleIndex pi = indexArray[i].index_of_atom[0];
      CBInterface::ParticleIndex pj = indexArray[i].index_of_atom[1];
      CBInterface::ParticleIndex pk = indexArray[i].index_of_atom[2];
      CBInterface::ParticleIndex pl = indexArray[i].index_of_atom[3];
# else  // TUNE_CBMAP
      CBInterface::ParticleIndex pi = getParticleLocation(particleIndexMap,pTorsion.id_of_atom[0]);
      CBInterface::ParticleIndex pj = getParticleLocation(particleIndexMap,pTorsion.id_of_atom[1]);
      CBInterface::ParticleIndex pk = getParticleLocation(particleIndexMap,pTorsion.id_of_atom[2]);
      CBInterface::ParticleIndex pl = getParticleLocation(particleIndexMap,pTorsion.id_of_atom[3]);
# endif  // TUNE_CBMAP
      const Position rkj = (*poscharge[pk.array])[pk.index].position - (*poscharge[pj.array])[pj.index].position;
      const Position rkl = (*poscharge[pk.array])[pk.index].position - (*poscharge[pl.array])[pl.index].position;
      const Position rik = (*poscharge[pi.array])[pi.index].position - (*poscharge[pk.array])[pk.index].position;
      const Position rij = (*poscharge[pi.array])[pi.index].position - (*poscharge[pj.array])[pj.index].position;
      const Position rjl = (*poscharge[pj.array])[pj.index].position - (*poscharge[pl.array])[pl.index].position;
#ifndef NDEBUG
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
#endif
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
#ifdef REDUCE_COS
      double cos_theta = dot;
#else
      double theta  = acos(dot)*(sign > 0.0 ? 1.0:-1.0);
      double cos_theta = cos(theta);
#endif
#ifdef MULTI_PARAM_TORSION
      for(int pti = 0;pti < getParameterNum(pTorsion);++pti)
#else
	int pti = 0;
#endif
      {
#ifdef REDUCE_COS
	calcInteraction(pTorsion, pti, cos_theta, sign, p, dp);
#else
	calcInteraction(pTorsion, pti, theta, p, dp); 
#endif
	cbenergy[thread_id] += p;
	force0 = (nv2*nv2INorm - cos_theta*nv1*nv1INorm)*nv1INorm;
	force1 = (nv1*nv1INorm - cos_theta*nv2*nv2INorm)*nv2INorm;
# ifdef CBOMP_TORSION
#  ifdef THREAD_FORCE_POINTER
	Force fi = dp*force0%rkj;
	Force fj = dp*(force0%rik - force1%rkl);
	Force fk = dp*(force1%rjl - force0%rij);
	Force fl = dp*force1%rkj;
	tf[pi.array][pi.index][0] += fi.x;
	tf[pi.array][pi.index][1] += fi.y;
	tf[pi.array][pi.index][2] += fi.z;
	tf[pj.array][pj.index][0] += fj.x;
	tf[pj.array][pj.index][1] += fj.y;
	tf[pj.array][pj.index][2] += fj.z;
	tf[pk.array][pk.index][0] += fk.x;
	tf[pk.array][pk.index][1] += fk.y;
	tf[pk.array][pk.index][2] += fk.z;
	tf[pl.array][pl.index][0] += fl.x;
	tf[pl.array][pl.index][1] += fl.y;
	tf[pl.array][pl.index][2] += fl.z;
#  else
	tmpforce[thread_id][pi.array][pi.index] += dp*force0%rkj;
	tmpforce[thread_id][pj.array][pj.index] += dp*(force0%rik - force1%rkl);
	tmpforce[thread_id][pk.array][pk.index] += dp*(force1%rjl - force0%rij);
	tmpforce[thread_id][pl.array][pl.index] += dp*force1%rkj;
#  endif
# else
	(*forcea[pi.array])[pi.index] += dp*force0%rkj;
	(*forcea[pj.array])[pj.index] += dp*(force0%rik - force1%rkl);
	(*forcea[pk.array])[pk.index] += dp*(force1%rjl - force0%rij);
	(*forcea[pl.array])[pl.index] += dp*force1%rkj;
#endif
      }

    }
  }
#ifdef BOND_COUNT
  printf("Num torsion %d %d %d %d %d\n",bc[1],bc[2],bc[3],bc[4],bc[0]);
#endif
# ifdef CBOMP_TORSION
  //  merge_tmpforce();
  /*
  {
    //    printf("CB reduct %d thread tmpforce\n",num_threads);
    int i;
    for(int t=0;t<num_threads;t++){
#pragma omp parallel for
#pragma loop norecurrence
      for(i=0;i<tmpforce[0][0].size();i++){
	(*forcea[0])[i].x += tmpforce[t][0][i].x;
	(*forcea[0])[i].y += tmpforce[t][0][i].y;
	(*forcea[0])[i].z += tmpforce[t][0][i].z;
      }
    }
#pragma omp parallel for
    for(i=0;i<tmpforce[0][1].size();i++){
      for(int t=0;t<num_threads;t++){
	(*forcea[1])[i] += tmpforce[t][1][i];
      }
    }
  }
  */
# endif
  for(int t=0;t<num_threads;t++){
    energy += cbenergy[t];
  }

  if(cancel_short){
    //  for (ITERATOR it = cbInterface->begin<T>(pContext);     it != cbInterface->end<T>(pContext);++it) {
    for(int i=0; i<torsion.size(); i++){
      const T& pTorsion = torsion[i];
# ifdef TUNE_CBMAP
      CBInterface::ParticleIndex pi = indexArray[i].index_of_atom[0];
      CBInterface::ParticleIndex pl = indexArray[i].index_of_atom[3];
# else  // TUNE_CBMAP
      CBInterface::ParticleIndex pi =
	//      cbInterface->getParticleLocation<PR>(pContext, pTorsion.id_of_atom[0]);
	getParticleLocation(particleIndexMap,pTorsion.id_of_atom[0]);
      CBInterface::ParticleIndex pl =
	//      cbInterface->getParticleLocation<PR>(pContext, pTorsion.id_of_atom[3]);
	getParticleLocation(particleIndexMap,pTorsion.id_of_atom[3]);
# endif  // TUNE_CBMAP

      {
	CBInterface::ForceLocation fi =
          getForceLocation(pTorsion.id_of_atom[0]); // cbInterface->getForceLocation(pContext, pTorsion, 0);
	CBInterface::ForceLocation fl =
          getForceLocation(pTorsion.id_of_atom[3]); // cbInterface->getForceLocation(pContext, pTorsion, 3);
	//! cancel 1-4 LJC
	if(isCalc14Interaction(pTorsion)){
	  Position dil;
	  dil = (*poscharge[pi.array])[pi.index].position - (*poscharge[pl.array])[pl.index].position; //cbInterface->calcDistance(pContext, pi, pl, dil);
	  //	  calcCB14Interaction(pContext, pi, pl, fi, fl, dil);
	  {
	    SpaceVector<double> force;
	    double r2;
	    double shortRangePotential;
	    double longRangePotential;

	    r2 = dil.norm2();
	    if (inCutoffSphere(r2)) {
	      calcInteraction(pi, pl,
			      scnbInv, sceeInv,
			      dil, r2,
			      shortRangePotential,
			      longRangePotential,
			      force); 
	      shortenergy += shortRangePotential; //cbInterface->addShortPotentialEnergy(pContext, shortRangePotential);
#ifdef DUMP_ENERGY_CONTENT
	      cb14energy += shortRangePotential;
#endif
	      //energy += longRangePotential;  //cbInterface->addLongPotentialEnergy(pContext, longRangePotential);
	      (*fi) += force; // cbInterface->addshortForce(pContext, fi,force);
	      (*fl) -= force; // cbInterface->subshortForce(pContext, fl,force);
	    }
	  }
	  subtractCBParticleForce(energy, shortenergy, shortvirial, pi, pl, fi, fl, dil); //subtractCBParticleForce(pContext, pi, pl, fi, fl, dil);
	}
      }
    }
  }
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
