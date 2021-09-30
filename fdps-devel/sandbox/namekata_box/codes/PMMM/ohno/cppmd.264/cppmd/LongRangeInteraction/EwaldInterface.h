#ifndef EWALDINTERFACE_H
#define EWALDINTERFACE_H

#include "EwaldInterfaceBase.h"
#include "EwaldMethod.h"

namespace EwaldModule {

typedef EwaldModule::EwaldInterfaceBase<EwaldMethod> EwaldOrgInterfaceBase;

class EwaldInterface : public EwaldOrgInterfaceBase {
public:
  EwaldInterface(int unitid, const LongRangeParameter& param, MPI_Comm lcomm=MPI_COMM_WORLD)
    : EwaldOrgInterfaceBase(param),
      nGrid(param.grid_num),
      kCutoff(static_cast<double>((nGrid.getComponentMin()-1)/2)),
      pkvec(), pdftr(), lazygrid(), 
#ifdef USE_MPI
      myworld(lcomm),
#endif
      inEwaldWorld(true) {
    //    std::cout << "construct PMEInterface" << std::endl;
#ifdef USE_MPI
    if(myworld==MPI_COMM_NULL){
      inEwaldWorld=false;
      std::cout << " this node out EwaldWorld" << std::endl;
      return;
    }else{
      std::cout << " this node in EwaldWorld" << std::endl;
    }
#endif

    int longrank=0;
#ifdef USE_MPI
    MPI_Comm_rank(myworld,&longrank);
#endif

    if (longrank == 0) {
      std::cout << "Grid " << param.grid_num << std::endl;
    }
    method = new EwaldMethod(this,
                             param.cutoff, param.alpha,
                             param.kCutoff, param.surfaceDipole);
    if (longrank != 0) {
      method->setShare(0.0);
    }
    isEwaldBaseNode = inEwaldWorld;
    method->initializeDFT();
    method->initializeKLoopRange();
  }
  
  void setSide(double side){
    method->setSide(side);
  }
  /* EwaldMethod -> EwaldInterface */
  double getKCutoff() { return kCutoff; }
  void setDFT(std::vector<WaveVector>& kvec, std::vector<DFTResult>& dftr) {
    pkvec = &kvec;
    pdftr = &dftr;
    if (lazygrid) {
      assert(pkvec != 0 && pdftr != 0);
      copyFromGrid(*lazygrid);
      lazygrid = 0;
    }
  }
  void getStartEnd(int& kstart, int &kend, int ksize) {}
#ifdef USE_MPI
  bool isEwaldNode() { return inEwaldWorld; }
  bool isWaveDevided() { return false; }
  bool isWaveParent() { return true; }
  int waveParentRank() { return 0; }
  MPI_Comm dftrCommunicator() { return myworld; }
  MPI_Comm waveCommunicator() { return MPI_COMM_NULL; }
#endif
  /* end of EwaldMethod -> EwaldInterface */

private:
  void chargeAssignMain(Context pContext) {
    method->clearDipoleMoment(getVolume(pContext));
    method->prepareDFT(cd, charge);
    method->calculateDipoleMoment(cd, charge);
    method->clearDFT();
    method->calculateDFT(pContext, cd, charge);
    method->completeDFT(pContext);
  }
  void backInterpolateMain(Context pContext) {
    method->calculateEwaldForce(pContext, cd, charge, fc);
    method->calculateDipoleForce(fc, charge);
    method->complete(pContext);
    method->completeForce(fc);
  }
  void copyToGrid(GridData& data) {
    assert(pkvec != 0 && pdftr != 0);
    const std::vector<WaveVector>& kvec = *pkvec;
    const std::vector<DFTResult>& dftr = *pdftr;
    for(std::vector<WaveVector>::size_type k = 0;k < kvec.size(); ++k) {
      int ix = (kvec[k].v[0]+nGrid[0])%nGrid[0];
      int iy = (kvec[k].v[1]+nGrid[1])%nGrid[1];
      int iz = (kvec[k].v[2]+nGrid[2])%nGrid[2];
      data.gridvalue[(ix*nGrid[1]+iy)*nGrid[2]+iz] = dftr[k].c;
      //std::cout << "C " << k << " " << kvec[k].v << " " << dftr[k].c << " " << dftr[k].s << " " << ix << " " << iy << " " << iz << std::endl;
      ix = (nGrid[0]-ix)%nGrid[0];
      iy = (nGrid[1]-iy)%nGrid[1];
      iz = (nGrid[2]-iz)%nGrid[2];
      data.gridvalue[(ix*nGrid[1]+iy)*nGrid[2]+iz] = dftr[k].s;
      //std::cout << "S " << k << " " << kvec[k].v << " " << dftr[k].c << " " << dftr[k].s << " " << ix << " " << iy << " " << iz << std::endl;
    }
  }
  void copyFromGrid(const GridData& data) {
    if (pkvec == 0 && pkvec == 0) {
      lazygrid = &data;
      return;
    }
    const std::vector<WaveVector>& kvec = *pkvec;
    std::vector<DFTResult>& dftr = *pdftr;
    for(std::vector<WaveVector>::size_type k = 0;k < kvec.size(); ++k) {
      int ix = (kvec[k].v[0]+nGrid[0])%nGrid[0];
      int iy = (kvec[k].v[1]+nGrid[1])%nGrid[1];
      int iz = (kvec[k].v[2]+nGrid[2])%nGrid[2];
      dftr[k].c = data.gridvalue[(ix*nGrid[1]+iy)*nGrid[2]+iz];
      //std::cout << "CF " << k << " " << kvec[k].v << " " << dftr[k].c << " " << dftr[k].s << " " << ix << " " << iy << " " << iz << std::endl;
      ix = (nGrid[0]-ix)%nGrid[0];
      iy = (nGrid[1]-iy)%nGrid[1];
      iz = (nGrid[2]-iz)%nGrid[2];
      dftr[k].s = data.gridvalue[(ix*nGrid[1]+iy)*nGrid[2]+iz];
      //std::cout << "SF " << k << " " << kvec[k].v << " " << dftr[k].c << " " << dftr[k].s << " " << ix << " " << iy << " " << iz << std::endl;
    }
  }

  SpaceVector<int> nGrid;
  double kCutoff;
  std::vector<WaveVector>* pkvec;
  std::vector<DFTResult>* pdftr;
  const GridData* lazygrid;
#ifdef USE_MPI
  MPI_Comm myworld;
#endif
  bool inEwaldWorld;
};
}
#endif
