#include "EwaldInterface.h"
#include "EwaldInterfaceImpl.h"
#include "EwaldMethod.h"

using namespace std;
using namespace EwaldModule;

void EwaldMethod::setSide(double side)
{
  EwaldBase::setSide(side);
  initializeDFT();
}
void EwaldMethod::initializeDFT()
{
  if (kvec.size() > 0) return;
  setWaveVector(ewaldInterface->getKCutoff());
  kstart = 0;
  kend = kvec.size();
  if (getSurfaceDipole()) {
    dftr.resize(kvec.size()+2);
  }
  else {
    dftr.resize(kvec.size());
  }
  ewaldInterface->setDFT(kvec, dftr);
}
void EwaldMethod::setWaveVector(const double kCutoff)
{
  int c;
  int totalNum;
  c = static_cast<int>(ceil(kCutoff)+.5);
  totalNum = (c*2+1)*(c*2+1)*(c*2+1);

  vector<WaveVector> tmp;
  tmp.reserve(totalNum);
  for(int i = -c; i < c+1; ++i){
    for(int j = -c; j < c+1; ++j){
      for(int k = -c; k < c+1; ++k){
	WaveVector wv;
	wv.norm2 = i*i+j*j+k*k;
	wv.v[0] = i;
	wv.v[1] = j;
	wv.v[2] = k;
	tmp.push_back(wv);
      }
    }
  }

  int vecMax2 = static_cast<int>(kCutoff*kCutoff);
  int kvecNum = 0;
  for(int v2 = 1; v2 < vecMax2+1; ++v2){
    for(int i = totalNum/2+1; i < totalNum; ++i){
      if(tmp[i].norm2 == v2){
        ++kvecNum;
      }
    }
  }
  //  cout << kvecNum << endl;

  kvec.reserve(kvecNum);

  for(int v2 = 1; v2 < vecMax2+1; ++v2){
    for(int i = totalNum/2+1; i < totalNum; ++i){
      if(tmp[i].norm2 == v2){
	kvec.push_back(tmp[i]);
      }
    }
  }
  /*
  cout << kvec.size() << " " << kvec.capacity() << endl;
  for(int i = 0; i < kvec.size(); ++i){
    cout << i << " " << kvec[i].norm2 << " " << kvec[i].v[0] << " " << kvec[i].v[1] << " " << kvec[i].v[2] << endl;
  }
  */
}
void EwaldMethod::addPotentialEnergy(EwaldInterface::Context pContext,
                                            double& potentialEnergy)
{
  const double volInv = 1.0/ewaldInterface->getVolume(pContext);
  const double pipi = M_PI*M_PI;
  const double alpha2Inv = 1.0/(getAlpha()*getAlpha());
  double pot = 0.0;
  for(int k = kstart;k < kend ; ++k){
    SpaceVector<double> wv(ewaldInterface->getWaveVector(pContext, kvec[k].v));
    double norm2 = wv.norm2();
    double factorPot = M_2_PI*volInv/norm2*exp(-pipi*alpha2Inv*norm2);
    pot += factorPot*(dftr[k].c*dftr[k].c + dftr[k].s*dftr[k].s);
    dftr[k].c *= factorPot;
    dftr[k].s *= factorPot;
  }
  potentialEnergy += 0.5 * pot * getShare();
  parallel->aggregate(potentialEnergy);
}

void EwaldMethod::initializeKLoopRange()
{
  ewaldInterface->getStartEnd(kstart, kend, kvec.size());
}

void EwaldMethod::clearDFT()
{
  parallel->deliverDataComplete();
  for(vector<DFTResult>::iterator it = dftr.begin(); it != dftr.end(); ++it){
    it->c = 0;
    it->s = 0;
  }
}

void EwaldMethod::completeDFT(EwaldInterface::Context pContext)
{
  parallel->calcDFT(pContext, this);
  if (getSurfaceDipole()) {
    SpaceVector<double>& dipoleMoment = getDipoleMoment();
    dftr[kvec.size()].c = dipoleMoment.x;
    dftr[kvec.size()].s = dipoleMoment.y;
    dftr[kvec.size()+1].c = dipoleMoment.z;
  }
  parallel->aggregateDFTResult(dftr);
  if (getSurfaceDipole()) {
    SpaceVector<double>& dipoleMoment = getDipoleMoment();
    dipoleMoment.x = dftr[kvec.size()].c;
    dipoleMoment.y = dftr[kvec.size()].s;
    dipoleMoment.z = dftr[kvec.size()+1].c;
  }
}

double EwaldMethod::calculateEwaldEnergy(EwaldInterface::Context pContext,
					 double &virial)
{
  double potEwald = 0.0;
  initializeDFT();
  addPotentialEnergy(pContext, potEwald);
  return potEwald;
}

void EwaldMethod::calculateDFT(EwaldInterface::Context pContext,
                                      const std::vector<Position>& cd,
				      const std::vector<double>& charge)
{
  const double twopi = 2.0*M_PI;
  for(int k = kstart;k < kend; ++k){
    SpaceVector<double> wv(ewaldInterface->getWaveVector(pContext, kvec[k].v));
    //      cout << k << " " << kvec[k].v << " " << wv << endl;
    for(vector<Position>::size_type i = 0; i < cd.size(); ++i){
      double inner = cd[i]*wv;
      //	cout << " " << i << " " << inner << " " << cd[i] << endl;
      dftr[k].c += cos(twopi*inner)*charge[i];
      dftr[k].s += sin(twopi*inner)*charge[i];
    }
  }
}

void EwaldMethod::calculateEwaldForce(EwaldInterface::Context pContext,
                                             const std::vector<Position>& cd,
				             const std::vector<double>& charge,
				             std::vector<Position>& fc)
{
  const double twopi = 2.0*M_PI;
  for(int k = kstart;k < kend; ++k){
    SpaceVector<double> wv(ewaldInterface->getWaveVector(pContext,
                                                           kvec[k].v));
    for(vector<Position>::size_type i = 0; i < cd.size(); ++i){
      double inner = cd[i]*wv;
      fc[i] += twopi*
          (dftr[k].c*sin(twopi*inner)-dftr[k].s*cos(twopi*inner))*charge[i]*wv;
    }
  }
}

void EwaldMethod::calculate(EwaldInterface_::Context pContext,
                                   const std::vector<Position>& cd,
                                   const std::vector<double>& charge,
                                   std::vector<Force>& fc,
                                   double& potentialEnergy)
{
  prepare(pContext, cd, charge);
  clearDipoleMoment(ewaldInterface->getVolume(pContext));
  prepareDFT(cd, charge);
  calculateDipoleMoment(cd, charge);
  clearDFT();
  calculateDFT(pContext, cd, charge);
  completeDFT(pContext);
  potentialEnergy += calculateEwaldEnergy(pContext);
  calculateEwaldForce(pContext, cd, charge, fc);
  calculateDipoleForce(fc, charge);
  complete(pContext);
  potentialEnergy += getSelfEnergy();
  potentialEnergy += getDipoleEnergy();
  completeForce(fc);
}
