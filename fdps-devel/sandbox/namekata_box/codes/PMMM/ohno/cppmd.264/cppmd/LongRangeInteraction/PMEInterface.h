#ifndef PMEINTERFACE_H
#define PMEINTERFACE_H

#include <mpi.h>
#include "EwaldInterfaceBase.h"
#include "PMEMethod.h"

namespace PMEModule {

typedef EwaldModule::EwaldInterfaceBase<PMEMethod> PMEInterfaceBase;

class PMEInterface : public PMEInterfaceBase {
public:
  typedef PMEInterface_::Context Context;
  typedef PMEInterfaceBase::BaseContext PMEContext;
  /* use PMEType in LongRangeParameter.h
  typedef PMEMethod::PMEType PMEType;
  static const PMEType PME = PMEMethod::PME;
  static const PMEType SmoothPME = PMEMethod::SmoothPME;
  */

  PMEInterface(int unitid, const LongRangeParameter& param, MPI_Comm lcomm=MPI_COMM_WORLD) 
      : PMEInterfaceBase(param),
#ifdef USE_MPI
        myworld(lcomm),
#endif
      inPMEWorld(true) {
    //    std::cout << "construct PMEInterface" << std::endl;
#ifdef USE_MPI
    if(myworld==MPI_COMM_NULL){
      inPMEWorld=false;
      //      std::cout << " this node out PMEWorld" << std::endl;
      return;
    }else{
      //      std::cout << " this node in PMEWorld" << std::endl;
    }
#endif

    int longrank=0;
#ifdef USE_MPI
    MPI_Comm_rank(myworld,&longrank);
#endif
    if (longrank == 0) {
      if (param.pmeType == SmoothPME) {
        std::cout << "SmoothPME ";
      }else{
        std::cout << "PME ";
      }
      std::cout << "order " << param.order <<" Grid " << param.grid_num << std::endl;
    }

    method = new PMEMethod(this, param.cutoff, param.order,
                           param.pmeType, param.alpha,
                           param.kCutoff, param.surfaceDipole);
    if (longrank != 0) {
      method->setShare(0.0);
    }
    isEwaldBaseNode = inPMEWorld;
    PMEContext context(param);
    method->setCommunicator();
#ifdef USE_MPI
    if(myworld!=MPI_COMM_NULL){
      MPI_Barrier(myworld);
    }
#endif
    fft = FFT3D::createFFT3D(param.grid_num[0],
                             param.grid_num[1],
                             param.grid_num[2]);
    if(DebugLog::verbose>1)std::cout << " createdFFT3D " <<std::endl;
    method->initialize(&context);
    if(DebugLog::verbose>1)std::cout << " PMEMethod initialize " <<std::endl;
  }
  ~PMEInterface() {
    if (method) {
      delete(method);
      method = 0;
    }
  }

  void initialize(){
    method->initialize_t();
  }

  /* PMEMethod -> PMEInterface methods */
  bool isPMENode() { return inPMEWorld; }
  FFT3D::Ptr getFFT(Context pContext) { return fft; }
#ifdef USE_MPI
  MPI_Comm getPMECommunicator() { return myworld; }
#endif
  /* end of PMEMethod -> PMEInterface methods */
private:
  void chargeAssignMain(Context pContext) {
    method->clearChargeDistribution(pContext);
    method->clearDipoleMoment(getVolume(pContext));
    method->calculateChargeDistribution(cd, charge);
    method->calculateDipoleMoment(cd, charge);
    method->completeChargeDistribution();
  }
  void backInterpolateMain(Context pContext) {
    method->calculateEwaldForce();
    method->addEwaldForce(fc);
    method->calculateDipoleForce(fc, charge);
  }
  void copyFromGrid(const GridData& data) {
    FFT3D::Real3D& r3d = fft->getReal();
    int xs = fft->getRealStart(0);
    int xe = fft->getRealEnd(0);
    int ys = fft->getRealStart(1);
    int ye = fft->getRealEnd(1);
    int zs = fft->getRealStart(2);
    int ze = fft->getRealEnd(2);
#ifndef NDEBUG
    int nx = fft->getSize(0);
#endif  // NDEBUG
    int ny = fft->getSize(1);
    int nz = fft->getSize(2);
    int dx = ny-ye+ys;
    assert(static_cast<size_t>(nx*ny*nz) == data.size);
    const double* p = data.gridvalue+((xs*ny)+ys)*nz;
    for (int ix=xs;ix<xe;++ix,p+=dx) {
      for (int iy=ys;iy<ye;++iy,p+=nz) {
        std::copy(p, p+ze-zs, &r3d[ix][iy][zs]);
      }
    }
  }
  void copyToGrid(GridData& data) {
    FFT3D::Real3D& r3d = fft->getReal();
    int xs = fft->getRealStart(0);
    int xe = fft->getRealEnd(0);
    int ys = fft->getRealStart(1);
    int ye = fft->getRealEnd(1);
    int zs = fft->getRealStart(2);
    int ze = fft->getRealEnd(2);
#ifndef NDEBUG
    int nx = fft->getSize(0);
#endif  // NDEBUG
    int ny = fft->getSize(1);
    int nz = fft->getSize(2);
    int dx = ny-ye+ys;
    assert(static_cast<size_t>(nx*ny*nz) == data.size);
    double* p = data.gridvalue+((xs*ny)+ys)*nz;
    for (int ix=xs;ix<xe;++ix,p+=dx) {
      for (int iy=ys;iy<ye;++iy,p+=nz) {
        double* pstart = &r3d[ix][iy][zs];
        std::copy(pstart, pstart+ze-zs, p);
      }
    }
  }
  
  FFT3D::Ptr fft;
#ifdef USE_MPI
  MPI_Comm myworld;
#endif
  bool inPMEWorld;
  bool selfEnergyCalculated;
};
}
#endif
