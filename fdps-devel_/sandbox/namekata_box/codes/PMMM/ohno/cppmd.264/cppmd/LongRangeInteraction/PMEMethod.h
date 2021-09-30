#ifndef PMEMETHOD_H
#define PMEMETHOD_H

#include <vector>
#include <LongRangeParameter.h>
#include "PMEInterfaceFwd.h"
#include "FFT3D.h"
#include "Lagrange.h"
#include "BSpline.h"
#include "PMEArray.h"
#include "EwaldBase.h"

namespace PMEModule {

class ParallelPMEMethod {
public:
  typedef ParallelPMEMethod* Ptr;

  ParallelPMEMethod(): cdPME(), chargePME(), fcPME(), itfcPME() {
    //    std::cout << "ParallelPMEMethod()" << std::endl;
  }
  virtual ~ParallelPMEMethod() {}

  virtual void setFFT(const FFT3D::Ptr& _fft) {}
  virtual void clear() {
    cdPME.clear();
    chargePME.clear();
  }
  virtual void setPositionsAndCharges(const std::vector<Position>& cd,
				      const std::vector<double>& charge) {
    assert(cd.size() == charge.size());
    cdPME.insert(cdPME.end(), cd.begin(), cd.end());
    chargePME.insert(chargePME.end(), charge.begin(), charge.end());
  }
  virtual void communicatePositionsAndCharges() {}
  virtual void communicateForces() {}
  virtual void communicateForces(SpaceVector<double>& dipoleMoment) {}
  virtual void setCommunicator() {}
  virtual void setSide(const SpaceVector<double>& side) {
    //    std::cout << "setSide" << std::endl;
  }
  virtual void setAxis(const SpaceVector<std::size_t>& axis) {}
  virtual void setLog(std::ostream& logfs) {}

  const std::vector<Position>& getPositions() {
    return cdPME;
  }
  const std::vector<double>& getCharges() {
    return chargePME;
  }
  std::vector<Force>& getForce() {
    return fcPME;
  }
  void clearForce() {
    fcPME.assign(cdPME.size(), Force());
    itfcPME = fcPME.begin();
  }
  void addForce(std::vector<Force>& fc) {
    for (std::vector<Force>::iterator it=fc.begin();it != fc.end();
      ++it,++itfcPME) {
      (*it) += (*itfcPME);
    }
  }
  int getIt() {
    return itfcPME-fcPME.begin();
  }
  virtual std::ostream& getLog() { return std::cout; }
protected:
  std::vector<Position> cdPME;
  std::vector<double> chargePME;
  std::vector<Force> fcPME;
  std::vector<Force>::const_iterator itfcPME;
};

ParallelPMEMethod::Ptr createParallelPMEMethod(PMEInterface* pmeInterface,
                                               int order);

typedef EwaldModule::EwaldBase EwaldBase;
typedef EwaldModule::SurfaceDipole SurfaceDipole;

class PMEMethod : public EwaldBase {
public:
  /*  use PMEType in LongRangeParameter.h
  typedef enum PMEType {
    PME,
    SmoothPME
  } PMEType;
  */
  PMEMethod(PMEInterface* _pmeInterface, const double cutoff,
            const int _order,
            const PMEType _pmeType=SmoothPME, const double alpha=0.0,
            const double kCutoff=0.0, const int surfaceDipole=0)
    : EwaldBase(cutoff, alpha, kCutoff, surfaceDipole),
      pmeInterface(_pmeInterface),
      pmeType(_pmeType), order(_order),
      nGrid(), indexes(), lagX(), lagY(), lagZ(), bspX(), bspY(), bspZ(),
      fft(), pmeArray(), axis(0,1,2),
      funcX(), funcY(), funcZ(), dfuncX(), dfuncY(), dfuncZ(),
      start(), end(), prepared(false),
      parallel(createParallelPMEMethod(pmeInterface, order))
  {
    initialize();
    num_threads=1;
  }
  virtual ~PMEMethod() {
    delete(pmeArray);
    delete(fft);
    delete(parallel);
    delete(bspX);
    delete(bspY);
    delete(bspZ);
    delete(lagX);
    delete(lagY);
    delete(lagZ);
  }

  void prepare(PMEInterface_::Context pContext,
               const std::vector<Position>& cd,
               const std::vector<double>& charge) {
    if (!prepared) {
      initialize(pContext);
      calculateSelfEnergy(charge);
      prepared = true;
    }
  }
  void calculate(PMEInterface_::Context pContext,
                 const std::vector<Position>& cd,
                 const std::vector<double>& charge,
                 std::vector<Force>& fc,
                 double& potentialEnergy,
		 double& virial);
  
  void initialize_t();
  void initialize(PMEInterface_::Context pContext);
  void clearChargeDistribution(PMEInterface_::Context pContext);
  void calculateChargeDistribution(const std::vector<Position>& cd,
                                   const std::vector<double>& charge);
  void completeChargeDistribution();
  double calculateEwaldEnergy(PMEInterface_::Context pContext, 
			      double &virial);
  void calculateEwaldForce();
  void addEwaldForce(std::vector<Force>& fc);
  void setAxisOrder(int ix, int iy, int iz) {
    axis[0] = ix;
    axis[1] = iy;
    axis[2] = iz;
  }
  void clear();
  void setLog(std::ostream& logfs) { parallel->setLog(logfs); }
  void pme_time_dump();

  void setCommunicator(){
    parallel->setCommunicator();
  }

private:

  std::size_t* getIndexPointer(int iaxis, double u, int& sinx, int& einx) {
    int k = 2*nGrid[iaxis]-int(floor(u));
    sinx = start[iaxis][k];
    einx = end[iaxis][k];
    return &indexes[iaxis][k+sinx];
  }
  static int minByAbs(int a, int b) {
    return (abs(a)<abs(b)?(a):(b));
  }
  double getScaledFractionalCoordinate(const Position& pos, int iaxis) const {
    return pos[axis[iaxis]]*gridInvs[iaxis];
  }
  void initialize();
  void calcb2(std::vector<double>& bv, int length, BSpline::Ptr& bsp);
  void calcInterpolation_nodiff(double ux, double uy, double uz);
  void calcInterpolation_nodiff_t(int t, double ux, double uy, double uz);
  void calcInterpolation(double ux, double uy, double uz);
  void calcInterpolation_t(int t, double ux, double uy, double uz);
  void addForce(Force &f, const Force& dEdu) const {
     f[axis[0]] -= dEdu[0]*gridInvs[0];
     f[axis[1]] -= dEdu[1]*gridInvs[1];
     f[axis[2]] -= dEdu[2]*gridInvs[2];
  }
  void calcQ(const std::vector<Position>& cd,const std::vector<double>& charge);
#if 0
  template<class PA>
    void calcQ(const PA& particlearray,
	       const std::vector<ParticleRange>& particlerange);
#endif
  void calcForce(std::vector<Force>& fc, const std::vector<Position>& cd,
	         const std::vector<double>& charge);
  void initializePMEArray() {
    if (pmeType == SmoothPME) {
      std::vector<double> b2X, b2Y, b2Z;
      calcb2(b2X, nGrid.getX(), bspX);
      calcb2(b2Y, nGrid.getY(), bspY);
      calcb2(b2Z, nGrid.getZ(), bspZ);
      pmeArray = PMEArray::Ptr(new PMEArray(fft, getAlpha(), axis,
					     b2X, b2Y, b2Z));
    }
    else {
      pmeArray = PMEArray::Ptr(new PMEArray(fft, getAlpha(), axis));
    }
  }
  
  PMEInterface* pmeInterface;
  const PMEType pmeType;
  SpaceVector<double> gridInvs;
  int order;
  SpaceVector<std::size_t> nGrid;
  std::vector<std::size_t> indexes[3];
  Lagrange::Ptr lagX, lagY, lagZ;
  Lagrange::Ptr *tlagX, *tlagY, *tlagZ;
  BSpline::Ptr bspX, bspY, bspZ;
  BSpline::Ptr *tbspX, *tbspY, *tbspZ;
  FFT3D::Ptr fft;
  PMEArray::Ptr pmeArray;
  SpaceVector<std::size_t> axis;
  double *funcX, *funcY, *funcZ;
  double *dfuncX, *dfuncY, *dfuncZ;
  double **tfuncX, **tfuncY, **tfuncZ;
  double **tdfuncX, **tdfuncY, **tdfuncZ;
  std::vector<std::size_t> start[3];
  std::vector<std::size_t> end[3];
  bool prepared;
  ParallelPMEMethod::Ptr parallel;

  int num_threads;

  int nx,ny,nz,nxyz;
  int qsx,qsy,qsz;
  double **q;
  FFT3D::Real3D_Ptr *tQ;
};
}
#endif
