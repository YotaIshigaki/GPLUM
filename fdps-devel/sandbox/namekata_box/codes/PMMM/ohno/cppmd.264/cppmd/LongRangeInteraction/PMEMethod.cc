#include "PMEInterface.h"
#include "PMEInterfaceImpl.h"
#include "PMEMethod.h"
#include "ErrorPos.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace PMEModule;

void PMEMethod::initialize()
{
  if (order <= 0) {
    throw runtime_error(errorPos("PME order must be positive"));
  }
  if (order % 2 > 0) {
    throw runtime_error(errorPos("PME order must be even"));
  }
  if (pmeType == SmoothPME) {
    bspX = BSpline::Ptr(new BSpline(order));
    bspY = BSpline::Ptr(new BSpline(order));
    bspZ = BSpline::Ptr(new BSpline(order));
    funcX = bspX->getFunction();
    funcY = bspY->getFunction();
    funcZ = bspZ->getFunction();
    dfuncX = bspX->getDerivative();
    dfuncY = bspY->getDerivative();
    dfuncZ = bspZ->getDerivative();
  }
  else {
    int p = order/2;
    lagX = Lagrange::Ptr(new Lagrange(p));
    lagY = Lagrange::Ptr(new Lagrange(p));
    lagZ = Lagrange::Ptr(new Lagrange(p));
    funcX = lagX->getFunction();
    funcY = lagY->getFunction();
    funcZ = lagZ->getFunction();
    dfuncX = lagX->getDerivative();
    dfuncY = lagY->getDerivative();
    dfuncZ = lagZ->getDerivative();
  }
}

void PMEMethod::initialize_t()
{
  //  std::cout << "PME buffers for thread " << std::endl;
#if 1
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
  //  std::cout << "PME Thread max " << num_threads << std::endl;
#else
  num_threads = 1;
#endif
  tfuncX = new double*[num_threads];
  tfuncY = new double*[num_threads];
  tfuncZ = new double*[num_threads];
  tdfuncX = new double*[num_threads];
  tdfuncY = new double*[num_threads];
  tdfuncZ = new double*[num_threads];
  if (pmeType == SmoothPME) {
    tbspX = new BSpline::Ptr[num_threads];
    tbspY = new BSpline::Ptr[num_threads];
    tbspZ = new BSpline::Ptr[num_threads];
    for(int t=0;t<num_threads;t++){
      tbspX[t] =  BSpline::Ptr(new BSpline(order));
      tbspY[t] =  BSpline::Ptr(new BSpline(order));
      tbspZ[t] =  BSpline::Ptr(new BSpline(order));
      tfuncX[t] = tbspX[t]->getFunction();
      tfuncY[t] = tbspY[t]->getFunction();
      tfuncZ[t] = tbspZ[t]->getFunction();
      tdfuncX[t] = tbspX[t]->getDerivative();
      tdfuncY[t] = tbspY[t]->getDerivative();
      tdfuncZ[t] = tbspZ[t]->getDerivative();
    }
  }else{
    tlagX = new Lagrange::Ptr[num_threads];
    tlagY = new Lagrange::Ptr[num_threads];
    tlagZ = new Lagrange::Ptr[num_threads];
    for(int t=0;t<num_threads;t++){
      tlagX[t] =  Lagrange::Ptr(new Lagrange(order));
      tlagY[t] =  Lagrange::Ptr(new Lagrange(order));
      tlagZ[t] =  Lagrange::Ptr(new Lagrange(order));
      tfuncX[t] = tlagX[t]->getFunction();
      tfuncY[t] = tlagY[t]->getFunction();
      tfuncZ[t] = tlagZ[t]->getFunction();
      tdfuncX[t] = tlagX[t]->getDerivative();
      tdfuncY[t] = tlagY[t]->getDerivative();
      tdfuncZ[t] = tlagZ[t]->getDerivative();
    }
  }
  q = new double*[num_threads];
  tQ = new FFT3D::Real3D_Ptr[num_threads];
  fft->get_r3d_n(nx,ny,nz);
  nxyz = nx*ny*nz;
  qsx = fft->getRealStart(0);
  qsy = fft->getRealStart(1);
  qsz = fft->getRealStart(2);
  FFT3D::Dims base = {{qsx,qsy,qsz}};
  for(int t=0;t<num_threads;t++){
    q[t] = new double[nxyz*2];
//    printf("q[%d] %016lx +%08x %d\n",t,q[t],nxyz,nxyz);
    tQ[t] = FFT3D::Real3D_Ptr(FFT3D::R3D::createArray(q[t],nx,ny,nz));
    tQ[t]->reindex(base);
  }
  //  std::cout << "PME buffers for thread Done " << std::endl;
#else
  //  std::cout << "PME buffers for thread Done (Debug)" << std::endl;
#endif
}

//  Ulrich Essmann 1995  (4.4) (4.8)
void PMEMethod::calcb2(vector<double>& bv, int length, 
                       BSpline::Ptr& bsp)
{
  const complex<double> factor(0.0, 2.0*M_PI/length);
  bsp->calculate(0.0);
  const double* const M = bsp->getFunction();
  for (int m=0;m<length;++m) {
    complex<double> s = 0.0;
    for (int k=0;k<order-1;++k) {
      s += M[k+1]*exp(factor*double(m*k));
    }
    bv.push_back(norm(exp(factor*double((order-1)*m))/s));
  }
}

void PMEMethod::calcInterpolation_nodiff(double ux, double uy, double uz)
{
  if (pmeType == SmoothPME) {
    bspX->calculate_nodiff(ux);
    bspY->calculate_nodiff(uy);
    bspZ->calculate_nodiff(uz);
  }
  else {
    lagX->calculate_nodiff(ux);
    lagY->calculate_nodiff(uy);
    lagZ->calculate_nodiff(uz);
  }
}

void PMEMethod::calcInterpolation_nodiff_t(int t, double ux, double uy, double uz)
{
  if (pmeType == SmoothPME) {
    tbspX[t]->calculate_nodiff(ux);
    tbspY[t]->calculate_nodiff(uy);
    tbspZ[t]->calculate_nodiff(uz);
  }
  else {
    tlagX[t]->calculate_nodiff(ux);
    tlagY[t]->calculate_nodiff(uy);
    tlagZ[t]->calculate_nodiff(uz);
  }
}

void PMEMethod::calcInterpolation(double ux, double uy, double uz)
{
  if (pmeType == SmoothPME) {
    bspX->calculate(ux);
    bspY->calculate(uy);
    bspZ->calculate(uz);
  }
  else {
    lagX->calculate(ux);
    lagY->calculate(uy);
    lagZ->calculate(uz);
  }
}

void PMEMethod::calcInterpolation_t(int t, double ux, double uy, double uz)
{
  if (pmeType == SmoothPME) {
    tbspX[t]->calculate(ux);
    tbspY[t]->calculate(uy);
    tbspZ[t]->calculate(uz);
  }
  else {
    tlagX[t]->calculate(ux);
    tlagY[t]->calculate(uy);
    tlagZ[t]->calculate(uz);
  }
}

void PMEMethod::initialize(PMEInterface::Context pContext)
{
  //  std::cout << "PMEMethod::initialize" << std::endl;
  parallel->setCommunicator();
  fft = pmeInterface->getFFT(pContext);
  SpaceVector<double> side;
  for (int i=0;i<SpaceVector<std::size_t>::Dim;++i) {
    side[i] = pmeInterface->getSide(pContext)[axis[i]];
    nGrid[i] = fft->getSize(i);
    if(DebugLog::verbose>1){
      std::cout << "nGrid["<<i<<"] " << nGrid[i] << std::endl;
    }
    gridInvs[i] = nGrid[i]/side[i];
    int indexLength = 3*nGrid[i]+order;
    indexes[i].clear();
    int p = order/2;
    int negativeBlock = nGrid[i]*(1+(indexLength-p)/nGrid[i]);
    for (int j=0;j<indexLength;++j) {

      indexes[i].push_back((negativeBlock+p-j)%nGrid[i]);
    }
  }
  parallel->setSide(pmeInterface->getSide(pContext));
  parallel->setAxis(axis);
  if(DebugLog::verbose>2){
    if (pmeType == SmoothPME) {
      cout << "SmoothPME " << order << " " << nGrid << endl;
    }
    else {
      cout << "PME " << order << " " << nGrid << endl;
    }
  }
  parallel->setFFT(fft);
  initializePMEArray();
  pmeArray->calculate(pmeInterface, pContext);

  for (int i=0;i<SpaceVector<std::size_t>::Dim;++i) {
    start[i].clear();
    end[i].clear();
    for (std::vector<std::size_t>::size_type j=0;
         j<indexes[i].size()-order;++j) {
      int sinx = 0;
      int einx = 0;
      for (int k=0;k<order;++k) {
        if (indexes[i][j+k] >=
            static_cast<std::size_t>(fft->getRealStart(i)) &&
            indexes[i][j+k] < static_cast<std::size_t>(fft->getRealEnd(i))) {
          sinx = k;
          break;
        }
      }
      for (int k=sinx;k<order;++k) {
        if (indexes[i][j+k] >=
            static_cast<std::size_t>(fft->getRealStart(i)) &&
            indexes[i][j+k] < static_cast<std::size_t>(fft->getRealEnd(i))) {
          einx = k+1;
        }
        else {
          break;
        }
      }
      for (int k=einx;k<order;++k) {
        if (indexes[i][j+k] >=
            static_cast<std::size_t>(fft->getRealStart(i)) &&
            indexes[i][j+k] < static_cast<std::size_t>(fft->getRealEnd(i))) {
          throw runtime_error(errorPos("PMEMethod pme_order grid_length unbalance. use smaller pme_order or smaller grid_length."));
        }
      }
      start[i].push_back(sinx);
      end[i].push_back(einx);
    }
  }
}

/// Ulrich Essmann 1995  (4.1) - (4.6)
void PMEMethod::calcQ(const vector<Position>& cd,
                      const vector<double>& charge)
{
  //  std::cout << "calcQ " << cd.size() << std::endl;
//  std::cout << fft->getRealStart(0) << ":" << fft->getRealEnd(0) << std::endl;
#ifdef K_SIMD 
  int cdsize = cd.size();
  const double (*cdc)[3] = (double (*)[3])(&(cd[0].x));
  const double *chargec = (double *)(&(charge[0]));
  const double *gridInvsc = (double *)(&(gridInvs[0]));
  const int ax0 = axis[0];
  const int ax1 = axis[1];
  const int ax2 = axis[2];
# ifndef CALCQ_OPENMP
  FFT3D::Real3D& Q = fft->getReal();
  const double* const Mx = funcX;
  const double* const My = funcY;
  const double* const Mz = funcZ;
# else
#  ifdef _OPENMP
# pragma omp parallel
#  endif
  {
#  ifdef _OPENMP
    const int t = omp_get_thread_num();
#  else
    const int t = 0;
#  endif
    FFT3D::Real3D& Q = *(tQ[t]);
    const double* const Mx = tfuncX[t];
    const double* const My = tfuncY[t];
    const double* const Mz = tfuncZ[t];
    memset(&(Q[qsx][qsy][qsz]),0,nxyz*sizeof(double));
#  ifdef _OPENMP
#pragma omp for schedule(runtime)
#  endif
# endif
  for (int i=0;i<cdsize;i++) {
    //    const double ux = getScaledFractionalCoordinate(cd[i], 0);
    //    const double uy = getScaledFractionalCoordinate(cd[i], 1);
    //    const double uz = getScaledFractionalCoordinate(cd[i], 2);
    double ux = cdc[i][ax0]*gridInvsc[0];
    double uy = cdc[i][ax1]*gridInvsc[1];
    double uz = cdc[i][ax2]*gridInvsc[2];
    int sx, ex;
    std::size_t* px0 = getIndexPointer(0, ux, sx, ex);
    int sy, ey;
    std::size_t* py0 = getIndexPointer(1, uy, sy, ey);
    int sz, ez;
    std::size_t* pz0 = getIndexPointer(2, uz, sz, ez);
# ifdef CALCQ_OPENMP
    calcInterpolation_nodiff_t(t,ux, uy, uz);
# else
    calcInterpolation_nodiff(ux, uy, uz);
# endif
    std::size_t* px = px0;
    double ci=charge[i];
//printf("Qrange %d:%d %d:%d %d:%d\n",px0[0],px0[ex-sx-1],py0[0],py0[ey-sy-1],pz0[0],pz0[ez-sz-1]);
    for (int jx=sx;jx<ex;++jx,++px) {
      std::size_t* py = py0;
      for (int jy=sy;jy<ey;++jy,++py) {
	double *Qc = &(Q[*px][*py][qsz]);
#pragma loop norecurrence
	for (int j=0;j<ez-sz;j++) {
    std::size_t pzj = pz0[j]-qsz;
# if 1
	  Qc[pzj] += ci*Mx[jx]*My[jy]*Mz[sz+j];
# else
	  Q[*px][*py][pz0[j]] += ci*Mx[jx]*My[jy]*Mz[sz+j];
# endif
	}
      }
    }
  }
# ifdef CALCQ_OPENMP
  }
//  printf("reduce Q %d %d %d %d\n",nx,ny,nz,nxyz);
  const int nt = omp_get_max_threads();
  FFT3D::Real3D& Q = fft->getReal();
  double *Qc = &(Q[qsx][qsy][qsz]);
  double **tQc;
//  printf("Qc %016lx\n",Qc);
  tQc = new double*[nt];
  for(int t=0;t<nt;t++){
    tQc[t] = &((*tQ[t])[qsx][qsy][qsz]);
//    printf("tQc[%d] %016lx\n",t,tQc[t]);
  }
  for(int t=0;t<nt;t++){
#pragma omp parallel for
#pragma loop norecurrence
    for(int i=0;i<nxyz;i++){
      Qc[i] += tQc[t][i];
    }
  }
//  printf("reduce Q %016lx done\n",Qc);
# endif
#else
  FFT3D::Real3D& Q = fft->getReal();
  const double* const Mx = funcX;    // equal to &(bspX->f[0])
  const double* const My = funcY;    // equal to &(bspY->f[0])
  const double* const Mz = funcZ;    // equal to &(bspZ->f[0])
  for (vector<Position>::size_type i=0;i<cd.size();++i) {
    const double ux = getScaledFractionalCoordinate(cd[i], 0);
    const double uy = getScaledFractionalCoordinate(cd[i], 1);
    const double uz = getScaledFractionalCoordinate(cd[i], 2);
    int sx, ex;
    std::size_t* px0 = getIndexPointer(0, ux, sx, ex);
    int sy, ey;
    std::size_t* py0 = getIndexPointer(1, uy, sy, ey);
    int sz, ez;
    std::size_t* pz0 = getIndexPointer(2, uz, sz, ez);
    calcInterpolation_nodiff(ux, uy, uz);   //  Ulrich Essmann 1995  (4.1)
    std::size_t* px = px0;
    for (int jx=sx;jx<ex;++jx,++px) {
      std::size_t* py = py0;
      for (int jy=sy;jy<ey;++jy,++py) {
        std::size_t* pz = pz0;
        for (int jz=sz;jz<ez;++jz,++pz) {
          Q[*px][*py][*pz] += charge[i]*Mx[jx]*My[jy]*Mz[jz];  //  Ulrich Essmann 1995  (4.6)
        }
      }
    }
  }
#endif
}

#if 0
template<class PA>
void PMEMethod::calcQ(const PA& particlearray,
		      const std::vector<ParticleRange>& particlerange)
{
  //  std::cout << "calcQ " << cd.size() << std::endl;
  FFT3D::Real3D& Q = fft->getReal();
  const double* const Mx = funcX;
  const double* const My = funcY;
  const double* const Mz = funcZ;
  int pr;
  for(pr=0;pr<particlerange.size();pr++){
    int i;
    int si=particlerange[pr].begin;
    int ei=particlerange[pr].end;
    for(i=si;i<ei;i++){
      const double ux = getScaledFractionalCoordinate(cd[i], 0);
      const double uy = getScaledFractionalCoordinate(cd[i], 1);
      const double uz = getScaledFractionalCoordinate(cd[i], 2);
      int sx, ex;
      std::size_t* px0 = getIndexPointer(0, ux, sx, ex);
      int sy, ey;
      std::size_t* py0 = getIndexPointer(1, uy, sy, ey);
      int sz, ez;
      std::size_t* pz0 = getIndexPointer(2, uz, sz, ez);
      calcInterpolation_nodiff(ux, uy, uz);
      std::size_t* px = px0;
      for (int jx=sx;jx<ex;++jx,++px) {
        std::size_t* py = py0;
	for (int jy=sy;jy<ey;++jy,++py) {
    std::size_t* pz = pz0;
	  for (int jz=sz;jz<ez;++jz,++pz) {
	    Q[*px][*py][*pz] += charge[i]*Mx[jx]*My[jy]*Mz[jz];
	  }
	}
      }
    }
  }
}
#endif

void PMEMethod::calcForce(vector<Force>& fc, const vector<Position>& cd,
                          const vector<double>& charge)
{
  FFT3D::Real3D& Q = fft->getReal();
#ifdef K_SIMD 
  int cdsize = cd.size();
  double (*fcc)[3] = (double (*)[3])(&(fc[0].x));
  const double (*cdc)[3] = (double (*)[3])(&(cd[0].x));
  const double *chargec = (double *)(&(charge[0]));
  const double *gridInvsc = (double *)(&(gridInvs[0]));
  const int ax0 = axis[0];
  const int ax1 = axis[1];
  const int ax2 = axis[2];
#ifdef _OPENMP
#pragma omp parallel 
#endif
  {
#ifdef _OPENMP
    const int t = omp_get_thread_num();
#else
    const int t = 0;
#endif
    const double* const Mx = tfuncX[t];
    const double* const My = tfuncY[t];
    const double* const Mz = tfuncZ[t];
    const double* const dMx = tdfuncX[t];
    const double* const dMy = tdfuncY[t];
    const double* const dMz = tdfuncZ[t];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int i=0;i<cdsize;++i) {
    //    const double ux = getScaledFractionalCoordinate(cd[i], 0);
    //    const double uy = getScaledFractionalCoordinate(cd[i], 1);
    //    const double uz = getScaledFractionalCoordinate(cd[i], 2);
    const double ux = cdc[i][ax0]*gridInvsc[0];
    const double uy = cdc[i][ax1]*gridInvsc[1];
    const double uz = cdc[i][ax2]*gridInvsc[2];
    int sx, ex;
    std::size_t* px0 = getIndexPointer(0, ux, sx, ex);
    int sy, ey;
    std::size_t* py0 = getIndexPointer(1, uy, sy, ey);
    int sz, ez;
    std::size_t* pz0 = getIndexPointer(2, uz, sz, ez);
    calcInterpolation_t(t,ux, uy, uz);
    double dEdu[3] = {0.0,0.0,0.0};
#if 0
    std::size_t* px = px0;
    for (int jx=sx;jx<ex;++jx,++px) {
      std::size_t* py = py0;
      for (int jy=sy;jy<ey;++jy,++py) {
        std::size_t* pz = pz0;
        for (int jz=sz;jz<ez;++jz,++pz) {
          dEdu[0] += chargec[i]*dMx[jx]* My[jy]* Mz[jz]*Q[*px][*py][*pz];
          dEdu[1] += chargec[i]* Mx[jx]*dMy[jy]* Mz[jz]*Q[*px][*py][*pz];
          dEdu[2] += chargec[i]* Mx[jx]* My[jy]*dMz[jz]*Q[*px][*py][*pz];
        }
      }
    }
#else
    std::size_t* px = px0;
    for (int jx=sx;jx<ex;++jx,++px) {
      std::size_t* py = py0;
      for (int jy=sy;jy<ey;++jy,++py) {
	//        std::size_t* pz = pz0;
	double *Qc = &(Q[*px][*py][qsz]);
        for (int j=0;j<ez-sz;j++) {
          std::size_t pzj = pz0[j]-qsz;
          dEdu[0] += chargec[i]*dMx[jx]* My[jy]* Mz[sz+j]*Qc[pzj];
          dEdu[1] += chargec[i]* Mx[jx]*dMy[jy]* Mz[sz+j]*Qc[pzj];
          dEdu[2] += chargec[i]* Mx[jx]* My[jy]*dMz[sz+j]*Qc[pzj];
        }
      }
    }
#endif
    //    addForce(fc[i], dEdu);
    fcc[i][ax0] -= dEdu[0]*gridInvsc[0];
    fcc[i][ax1] -= dEdu[1]*gridInvsc[1];
    fcc[i][ax2] -= dEdu[2]*gridInvsc[2];
  }
  }
#else
  const double* const Mx = funcX;
  const double* const My = funcY;
  const double* const Mz = funcZ;
  const double* const dMx = dfuncX;
  const double* const dMy = dfuncY;
  const double* const dMz = dfuncZ;
  for (vector<Position>::size_type i=0;i<cd.size();++i) {
    const double ux = getScaledFractionalCoordinate(cd[i], 0);
    const double uy = getScaledFractionalCoordinate(cd[i], 1);
    const double uz = getScaledFractionalCoordinate(cd[i], 2);
    int sx, ex;
    std::size_t* px0 = getIndexPointer(0, ux, sx, ex);
    int sy, ey;
    std::size_t* py0 = getIndexPointer(1, uy, sy, ey);
    int sz, ez;
    std::size_t* pz0 = getIndexPointer(2, uz, sz, ez);
    calcInterpolation(ux, uy, uz);
    Force dEdu;
    std::size_t* px = px0;
    for (int jx=sx;jx<ex;++jx,++px) {    //  Ulrich Essmann 1995  (4.9)
      std::size_t* py = py0;
      for (int jy=sy;jy<ey;++jy,++py) {
        std::size_t* pz = pz0;
        for (int jz=sz;jz<ez;++jz,++pz) { 
          dEdu.x += charge[i]*dMx[jx]* My[jy]* Mz[jz]*Q[*px][*py][*pz];
          dEdu.y += charge[i]* Mx[jx]*dMy[jy]* Mz[jz]*Q[*px][*py][*pz];
          dEdu.z += charge[i]* Mx[jx]* My[jy]*dMz[jz]*Q[*px][*py][*pz];
        }
      }
    }
    addForce(fc[i], dEdu);
  }
#endif

}

void PMEMethod::clearChargeDistribution(PMEInterface::Context pContext)
{
  fft->assign(0.0);
  parallel->clear();
  SpaceVector<double> side;
  for (int i=0;i<SpaceVector<std::size_t>::Dim;++i) {
    side[i] = pmeInterface->getSide(pContext)[axis[i]];
    gridInvs[i] = nGrid[i]/side[i];
  }
  parallel->setSide(pmeInterface->getSide(pContext));
}

void PMEMethod::calculateChargeDistribution(const vector<Position>& cd,
                                            const vector<double>& charge)
{
  parallel->setPositionsAndCharges(cd, charge);
}

void PMEMethod::completeChargeDistribution()
{
  parallel->communicatePositionsAndCharges();
  calcQ(parallel->getPositions(), parallel->getCharges());
}

double PMEMethod::calculateEwaldEnergy(PMEInterface::Context pContext,
				       double &virial)
{
  return fft->convolution(pmeArray->calculate(pmeInterface, pContext),
			  pmeArray->getvirialarray(),virial);
}

void PMEMethod::calculateEwaldForce()
{
  parallel->clearForce();
  calcForce(parallel->getForce(),
            parallel->getPositions(), 
            parallel->getCharges());
  if (getSurfaceDipole()) {
    parallel->communicateForces(getDipoleMoment());
  }
  else {
    parallel->communicateForces();
  }
}

void PMEMethod::addEwaldForce(vector<Force>& fc)
{
  parallel->addForce(fc);
}

void PMEMethod::clear()
{
  delete(pmeArray);
  pmeArray = 0;
  delete(fft);
  fft = 0;
}

void PMEMethod::calculate(PMEInterface::Context pContext,
                          const std::vector<Position>& cd,
                          const std::vector<double>& charge,
                          std::vector<Force>& fc,
                          double& potentialEnergy,
			  double& virial)
{
  prepare(pContext, cd, charge);
  clearChargeDistribution(pContext);
  clearDipoleMoment(pmeInterface->getVolume(pContext));
  calculateChargeDistribution(cd, charge);
  calculateDipoleMoment(cd, charge);
  completeChargeDistribution();
  potentialEnergy += calculateEwaldEnergy(pContext,virial);
  calculateEwaldForce();
  addEwaldForce(fc);
  calculateDipoleForce(fc, charge);
  potentialEnergy += getSelfEnergy();
  potentialEnergy += getDipoleEnergy();
}
