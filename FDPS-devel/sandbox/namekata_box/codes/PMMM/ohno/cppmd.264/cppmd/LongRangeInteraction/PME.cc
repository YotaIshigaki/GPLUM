/*
#define HAVE_BOOST_ARRAY_HPP
#define HAVE_BOOST_MULTI_ARRAY_HPP
#define USE_MPI
*/

#include "UseMPI.h"
#include "Common.h"
#include "LongRangeParameter.h"
#include "SurfaceDipole.h"
#ifdef CPPMD_ENABLE_SIMPLE_FFT
#include "FFT.h"
typedef FFTModule::FFT3D FFT3D;
#else
#include "MPIFFT3D.h"
#endif
#include <bitset>
#include "BSpline.h"
#include "Lagrange.h"
#include "PMEArray.h"
#include "ErrorPos.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _OPENMP
# ifdef PME_OPENMP
#  undefine PME_OPENMP
# endif
#endif

typedef EwaldModule::SurfaceDipole SurfaceDipole;

namespace PMEModule {
  typedef struct {
    int start;
    int end;
  } StartEnd;
  typedef struct {
    int border;
    int nodeid;
    int width;
  } BorderNode;
  typedef enum {
    KEEP,
    SEND,
    KEEP_SEND_NUM
  } KeepSendFlag;
  typedef std::bitset<KEEP_SEND_NUM> KeepSend;
  typedef struct {
    int bitIndex;
    int destNode;
    bool sendFirst;
    bool keepAll;
    bool sendOnly;
    bool recvOnly;
    int tag;
    std::vector<KeepSend> keepSends;
  } SendPlan;
  typedef union {
    struct {
      int bitIndex;
      int gridIndex;
      int particleIndex;
    } bg;
    double data;
  } BitGrid;
  struct lessBorderNode {
    bool operator()(BorderNode q1, BorderNode q2) const
    {
      return
          (q1.width > q2.width)
          || (((q1.width == q2.width) && (q1.border < q2.border))
              || ((q1.width == q2.width) && (q1.border == q2.border)
                  && (q1.nodeid < q2.nodeid)
                 )
             );
    }
  };

  static const int PME_TAG = 100;
  static const int PME_BACK_TAG = 150;

//PMEMethod* method;                     //  EwaldInterfaceBase.h : 181
std::vector<Position> cd;              //  EwaldInterfaceBase.h : 182
std::vector<double> charge;            //  EwaldInterfaceBase.h : 183
std::vector<Force> fc;                 //  EwaldInterfaceBase.h : 184
bool selfEnergyCalculated;             //  EwaldInterfaceBase.h : 185
bool isEwaldBaseNode;                  //  EwaldInterfaceBase.h : 186, PMEInterface.h : 155

MPI_Comm myworld;                      //  PMEInterface.h : 152, ParallelPMEMethod.cc : 128

bool inPMEWorld;                       //  PMEInterface.h : 154

LongRangeParameter param;
// from PMEMethod : EwaldBase
double cutoff;     // EwaldBase.h : 60
double alpha;       // EwaldBase.h : 61
double kCutoff;     // EwaldBase.h : 62
int calcSurfaceDipole; // EwaldBase.h : 63
SurfaceDipole surface; // EwaldBase.h : 64
double selfEnergy;  // EwaldBase.h : 65
double shareValue;  // EwaldBase.h : 66
double erfcCutValue;  // EwaldBase.h : 67
double ewaldExpCutValue;  // EwaldBase.h : 68

std::vector<Position> cdPME;   // PMEMethod.h : 69
std::vector<double> chargePME; // PMEMethod.h : 70
std::vector<Force> fcPME; // PMEMethod.h : 71
std::vector<Force>::const_iterator itfcPME; // PMEMethod.h : 72

int myrank;      // ParallelPMEMethod.cc : 129
int nodeNum;      // ParallelPMEMethod.cc : 130
int gridNum;      // ParallelPMEMethod.cc : 131
std::vector<StartEnd> sevec;      // ParallelPMEMethod.cc : 132
SpaceVector<double> side;      // ParallelPMEMethod.cc : 133
std::vector<BorderNode> borderNodes;     // ParallelPMEMethod.cc : 134
int bitIndex;     // ParallelPMEMethod.cc : 135
std::size_t particleIndex;     // ParallelPMEMethod.cc : 136
std::size_t particleNum;     // ParallelPMEMethod.cc : 137
std::vector<SendPlan> sendPlans;     // ParallelPMEMethod.cc : 138
std::vector<SendPlan> sendBackPlans;     // ParallelPMEMethod.cc : 139
int overlap;     // ParallelPMEMethod.cc : 140
int divideAxis;     // ParallelPMEMethod.cc : 141
std::vector<Position>cdSend;     // ParallelPMEMethod.cc : 142
std::vector<double>chargeSend;     // ParallelPMEMethod.cc : 143
std::vector<BitGrid>bitGridSend;     // ParallelPMEMethod.cc : 145
std::vector<Position>cdKeep;     // ParallelPMEMethod.cc : 145
std::vector<double>chargeKeep;     // ParallelPMEMethod.cc : 146
std::vector<BitGrid>bitGridKeep;     // ParallelPMEMethod.cc : 147
std::vector<BitGrid>bitGridPME;     // ParallelPMEMethod.cc : 148

std::vector<Force>fcKeep;     // ParallelPMEMethod.cc : 151
std::vector<Force>fcSend;     // ParallelPMEMethod.cc : 152
bool withDipoleMoment;     // ParallelPMEMethod.cc : 153
SpaceVector<double> dipoleMoment;     // ParallelPMEMethod.cc : 154
SpaceVector<double> recvDipoleMoment;     // ParallelPMEMethod.cc : 155
#ifdef CPPMD_ENABLE_SIMPLE_FFT
  FFTModule::Communicator comms;     // ParallelPMEMethod.cc : 156
#else
Communicator comms;     // ParallelPMEMethod.cc : 156
#endif
  //NullStream logfs;     // ParallelPMEMethod.cc : 157

  //PMEInterface* pmeInterface;  // ParallelPMEMethod.cc : 127, PMEMethod.h : 200
PMEType pmeType;  // PMEMethod.h 201
SpaceVector<double> gridInvs;  // PMEMethod.h 202
int order; // PMEMethod.h 203
SpaceVector<std::size_t> nGrid; // PMEMethod.h 204
std::vector<std::size_t> indexes[3]; // PMEMethod.h 205
Lagrange::Ptr lagX, lagY, lagZ; // PMEMethod.h 206
Lagrange::Ptr *tlagX, *tlagY, *tlagZ;    /// for OpenMP
BSpline::Ptr bspX, bspY, bspZ; // PMEMethod.h 208
BSpline::Ptr *tbspX, *tbspY, *tbspZ;    /// for OpenMP
FFT3D::Ptr fft; // PMEMethod.h 210
PMEArray::Ptr pmeArray;   // PMEMethod.h 211
SpaceVector<std::size_t> axis;   // PMEMethod.h 212
double *funcX, *funcY, *funcZ;   // PMEMethod.h 213
double *dfuncX, *dfuncY, *dfuncZ;   // PMEMethod.h 214
double **tfuncX, **tfuncY, **tfuncZ;
double **tdfuncX, **tdfuncY, **tdfuncZ;
std::vector<std::size_t> start[3]; // PMEMethod.h 217
std::vector<std::size_t> end[3]; // PMEMethod.h 218
bool prepared;  // PMEMethod.h : 219
  //ParallelPMEMethod::Ptr parallel;  // PMEMethod.h : 220

int num_threads;   // PMEMethod.h : 222


int nx,ny,nz,nxyz;
int qsx,qsy,qsz;
double **q;
FFT3D::Real3D_Ptr *tQ;


// EwaldInterfaceBase.h : 27 - 31
inline void clear_cd_charge()
{
  cd.clear();
  charge.clear();
}

/*
// EwaldInterfaceBase.h : 174 - 179
// EwaldBase.h : 28 - 33
*/
void calculateSelfEnergy()
{
  if (!selfEnergyCalculated) {
    for (std::vector<double>::const_iterator it = charge.begin(); it != charge.end();++it) {
      selfEnergy += (-(*it)*(*it)*alpha*M_2_SQRTPI*.5);
    }
    selfEnergyCalculated = true;
#ifdef DUMP_ENERGY_CONTENT
    printf("selfenergy at long %f with %d atom with factor %f alpha %f\n",selfEnergy,charge.size(),alpha*M_2_SQRTPI*.5,alpha);
#endif
  }
}


// PMEMethod::clearChargeDistribution PMEMethod.cc : 457 - 467
void clearChargeDistribution(LongRangeParameter param)
{
  fft->assign(0.0);  // FFT3D.h : 104
  cdPME.clear();
  chargePME.clear();
  side = param.boxSize;
  for (int i=0;i<SpaceVector<std::size_t>::Dim;++i) {
    gridInvs[i] = nGrid[i]/side[axis[i]];
  }
}


// EwaldBase::clearDipoleMoment EwaldBase.h : 34 - 38
void clearDipoleMoment(double volume) {
  if (calcSurfaceDipole) {
    surface.initialize(volume);
  }
}

/*
  @param[in] cd positions of particles
  @param[out] charge chareges of particles
  \sa{cdPME} positions of particles
  \sa{chargePME} charge chareges of particles
 */
// PMEMethod::calculateChargeDistribution PMEMethod.cc : 469 - 473
void calculateChargeDistribution(const std::vector<Position>& cd,
				 const std::vector<double>& charge)
{
  /*
  parallel->setPositionsAndCharges(cd, charge); // PMEMethod.h : 29
  */
  assert(cd.size() == charge.size());
  cdPME.insert(cdPME.end(), cd.begin(), cd.end());
  chargePME.insert(chargePME.end(), charge.begin(), charge.end());
}

// EwaldBase::calculateDipoleMoment EwaldBase.h : 39 - 45
void calculateDipoleMoment(const std::vector<Position>& cd,
			   const std::vector<double>& charge)
{
  if (calcSurfaceDipole) {
    surface.addDipoleMoment(cd, charge);
  }
}


// PMEMethod.h 165 - 167
double getScaledFractionalCoordinate(const Position& pos, int iaxis) {
  return pos[axis[iaxis]]*gridInvs[iaxis];
}

/*
  @param[in] iaxis id of axis
  @param[in] u position
  @param[out] sinx start index
  @param[out] einx end index
  @retrun pointer for indexes of grid point for position u
  \sa{nGrid} number of grid point
  \sa{start} start index for axis iaxis, grid point k
  \sa{end} end index for axis iaxis, grid point k
  \sa{indexes} indexes of grid point for axis iaxis
 */
// PMEMethod.h 156 - 161
std::size_t* getIndexPointer(int iaxis, double u, int& sinx, int& einx) {
  int k = 2*nGrid[iaxis]-int(floor(u));
  sinx = start[iaxis][k];
  einx = end[iaxis][k];
  return &indexes[iaxis][k+sinx];
}

// PMEMethod.cc 113 - 125
/*
  @param[in] ux x position of particle
  @param[in] uy y position of particle
  @param[in] uz z position of particle
  \sa{bspX} pointer of B-spline instance for x
  \sa{bspY} pointer of B-spline instance for y
  \sa{bspZ} pointer of B-spline instance for z
  \sa{pmeType} indicate type of PME such as SmoothPME 
 */
void calcInterpolation_nodiff(double ux, double uy, double uz)
{
  if (pmeType == SmoothPME) {
    bspX->calculate_nodiff(ux);    // BSpline.h : 22
    bspY->calculate_nodiff(uy);
    bspZ->calculate_nodiff(uz);
  }
  else {
    lagX->calculate_nodiff(ux);    // Lagrange.h : 17
    lagY->calculate_nodiff(uy);
    lagZ->calculate_nodiff(uz);
  }
}


// PMEMethod.h : 229 - 297
/*
  @param[in] cd positions of particles
  @param[in] charge charges of particles
  \sa{gridInvs} scale factor, convert position to grid
  \sa{fancX}, \sa{fancY}, \sa{fancZ} assign functions (pointer of matrix in BSpline)
 */
void calcQ(const std::vector<Position>& cd, 
	   const std::vector<double>& charge)
{
  //  std::cout << "calcQ " << cd.size() << std::endl;
  FFT3D::Real3D& Q = fft->getReal();
  const double* const Mx = funcX;    // equal to &(bspX->f[0])
  const double* const My = funcY;    // equal to &(bspY->f[0])
  const double* const Mz = funcZ;    // equal to &(bspZ->f[0])
#ifdef K_SIMD 
  int cdsize = cd.size();
  const double (*cdc)[3] = (double (*)[3])(&(cd[0].x));
  const double *chargec = (double *)(&(charge[0]));
  const double *gridInvsc = (double *)(&(gridInvs[0]));
  const int ax0 = axis[0];
  const int ax1 = axis[1];
  const int ax2 = axis[2];
  for (int i=0;i<cdsize;++i) {
    const double ux = cdc[i][ax0]*gridInvsc[0];
    const double uy = cdc[i][ax1]*gridInvsc[1];
    const double uz = cdc[i][ax2]*gridInvsc[2];
    int sx, ex;
    std::size_t* px0 = getIndexPointer(0, ux, sx, ex); // PMEMethod.h 156
    int sy, ey;
    std::size_t* py0 = getIndexPointer(1, uy, sy, ey); // PMEMethod.h 156
    int sz, ez;
    std::size_t* pz0 = getIndexPointer(2, uz, sz, ez); // PMEMethod.h 156
    calcInterpolation_nodiff(ux, uy, uz);   // calculate bsp{X,Y,Z}->fc
    std::size_t* px = px0;
    double ci=charge[i];
    for (int jx=sx;jx<ex;++jx,++px) {
      std::size_t* py = py0;
      for (int jy=sy;jy<ey;++jy,++py) {
	double *Qc = &(Q[*px][*py][0]);
	for (int j=0;j<ez-sz;j++) {
	  std::size_t pzj = pz0[j];
	  Qc[pzj] += ci*Mx[jx]*My[jy]*Mz[sz+j];
	}
      }
    }
  }
#else
  for (std::vector<Position>::size_type i=0;i<cd.size();++i) {
    const double ux = getScaledFractionalCoordinate(cd[i], 0); // PMEMethod.h 165
    const double uy = getScaledFractionalCoordinate(cd[i], 1);
    const double uz = getScaledFractionalCoordinate(cd[i], 2);
    int sx, ex;
    std::size_t* px0 = getIndexPointer(0, ux, sx, ex); // PMEMethod.h 156
    int sy, ey;
    std::size_t* py0 = getIndexPointer(1, uy, sy, ey);
    int sz, ez;
    std::size_t* pz0 = getIndexPointer(2, uz, sz, ez);
    calcInterpolation_nodiff(ux, uy, uz);         // PMEMethod.cc : 113
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
#endif
}


// PMEMethod::completeChargeDistribution PMEMethod.cc : 475 - 479
void completeChargeDistribution()
{
  calcQ(cdPME,chargePME);
}

void chargeAssignMain()
{
  clearChargeDistribution(param); 
  clearDipoleMoment(param.boxSize[0]*param.boxSize[1]*param.boxSize[2]); 
  calculateChargeDistribution(cd, charge); 
  calculateDipoleMoment(cd, charge);  
  completeChargeDistribution();           
}

// EwaldInterfaceBase.h : 127 - 137
template<class PA>
void convert(PA& particlearray,
	     const std::vector<ParticleRange>& selfrange) {
  for (std::vector<ParticleRange>::size_type it=0;it < selfrange.size(); ++it) {
    for (int i = selfrange[it].begin; i < selfrange[it].end; ++i) {
      cd.push_back(getpos(particlearray,i));
      charge.push_back(getcharge(particlearray,i));
    }
  }
}

// PMELongRangeInteraction.cc : 24
  /* 
     typedef EwaldModule::EwaldInterfaceBase<PMEMethod> PMEInterfaceBase;  // PMEInterface.h : 10
     class PMEInterface : public PMEInterfaceBase  {...};  //  PMEInterface.h : 12 - 156
  EwaldInterfaceBase.h : 48
  */
template<class PA, class GPA>
void assign(PA& particlearray,
	    const std::vector<ParticleRange>& selfrange,
	    GPA& ghost,
	    const std::vector<ParticleRange>& ghostrange,
	    //	    GridData& gridcharge,
	    const std::vector<ParticleRange>& self_selfenergy_range,
	    const std::vector<ParticleRange>& ghost_selfenergy_range)
{
  //  std::cout << "assign " << std::endl;
  {
    if (!isEwaldBaseNode) return;
    //  std::cout << "isEwaldBaseNode " << std::endl;
    if (!selfEnergyCalculated) {
      clear_cd_charge();
      convert(particlearray, self_selfenergy_range);
      convert(ghost, ghost_selfenergy_range);
      calculateSelfEnergy();
    }
    //  std::cout << " clear_cd_charge " << std::endl;
    clear_cd_charge();
    //  std::cout << " convert self " << std::endl;
    convert(particlearray, selfrange);
    //  std::cout << " convert ghost " << std::endl;
    convert(ghost, ghostrange);
    //    std::cout << " chargeAssignMain " << std::endl;
    chargeAssignMain();         //    chargeAssignMain(&context);
  }
}

void solvePoisson(double& energy)
{
  if (!isEwaldBaseNode) {
    energy = 0.0;
  }else{
    energy = fft->convolution(pmeArray->calculate(param.boxSize));
  }
}

void solvePoisson(double& energy, double& virial)
{
  if (!isEwaldBaseNode) {
    energy = 0.0;
  }else{
    energy = fft->convolution(pmeArray->calculate(param.boxSize),pmeArray->getvirialarray(),virial);
  }
}


// PMEMethod.h  54 - 57
void clearForce()
{
  fcPME.assign(cdPME.size(), Force());
  itfcPME = fcPME.begin();
}

// PMEMethod.cc 127 - 139
void calcInterpolation(double ux, double uy, double uz)
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


// PMEMethod.cc 141 - 153 
// for OpenMP
void calcInterpolation_t(int t, double ux, double uy, double uz)
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

void addForce(Force &f, const Force& dEdu) 
{
  f[axis[0]] -= dEdu[0]*gridInvs[0];
  f[axis[1]] -= dEdu[1]*gridInvs[1];
  f[axis[2]] -= dEdu[2]*gridInvs[2];
}

// PMEMethod.cc 339 - 455
void calcForce(std::vector<Force>& fc, const std::vector<Position>& cd,
	       const std::vector<double>& charge)
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
#ifdef PME_OPENMP
#pragma omp parallel 
#endif
  {
#ifdef PME_OPENMP
    const int t = omp_get_thread_num();
    const double* const Mx = tfuncX[t];
    const double* const My = tfuncY[t];
    const double* const Mz = tfuncZ[t];
    const double* const dMx = tdfuncX[t];
    const double* const dMy = tdfuncY[t];
    const double* const dMz = tdfuncZ[t];
#else
    const double* const Mx = funcX;
    const double* const My = funcY;
    const double* const Mz = funcZ;
    const double* const dMx = dfuncX;
    const double* const dMy = dfuncY;
    const double* const dMz = dfuncZ;
#endif
#ifdef PME_OPENMP
#pragma omp for schedule(static)
#endif
    for (int i=0;i<cdsize;++i) {
      const double ux = cdc[i][ax0]*gridInvsc[0];
      const double uy = cdc[i][ax1]*gridInvsc[1];
      const double uz = cdc[i][ax2]*gridInvsc[2];
      int sx, ex;
      std::size_t* px0 = getIndexPointer(0, ux, sx, ex);
      int sy, ey;
      std::size_t* py0 = getIndexPointer(1, uy, sy, ey);
      int sz, ez;
      std::size_t* pz0 = getIndexPointer(2, uz, sz, ez);
#ifdef PME_OPENMP
      calcInterpolation_t(t,ux, uy, uz);
#else
      calcInterpolation(ux, uy, uz);
#endif
      double dEdu[3] = {0.0,0.0,0.0};
      std::size_t* px = px0;
      for (int jx=sx;jx<ex;++jx,++px) {
	std::size_t* py = py0;
	for (int jy=sy;jy<ey;++jy,++py) {
	  //        std::size_t* pz = pz0;
	  double *Qc = &(Q[*px][*py][0]);
	  for (int j=0;j<ez-sz;j++) {
	    std::size_t pzj = pz0[j];
	    dEdu[0] += chargec[i]*dMx[jx]* My[jy]* Mz[sz+j]*Qc[pzj];
	    dEdu[1] += chargec[i]* Mx[jx]*dMy[jy]* Mz[sz+j]*Qc[pzj];
	    dEdu[2] += chargec[i]* Mx[jx]* My[jy]*dMz[sz+j]*Qc[pzj];
	  }
	}
      }
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
  for (std::vector<Position>::size_type i=0;i<cd.size();++i) {
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
    for (int jx=sx;jx<ex;++jx,++px) {
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


// PMEMethod.cc 486 - 498
void calculateEwaldForce()
{
  clearForce();
  calcForce(fcPME, cdPME, chargePME);
}

// PMEMethod.h 58 - 63
void addForce(std::vector<Force>& fc) 
{
  for (std::vector<Force>::iterator it=fc.begin();it != fc.end(); ++it,++itfcPME) {
    (*it) += (*itfcPME);
  }
}

// PMEMethod.cc 500 - 504
void addEwaldForce(std::vector<Force>& fc)
{
  addForce(fc); //  parallel->addForce(fc);
}

// EwaldBace.h 46 - 52
void calculateDipoleForce(std::vector<Force>& fc, const std::vector<double>& charge)
{
  if (calcSurfaceDipole) {
    surface.addDipoleForce(fc, charge);
  }
}

// PMEInterface.h 99 - 103
void backInterpolateMain()
{
  calculateEwaldForce();
  addEwaldForce(fc);
  calculateDipoleForce(fc, charge);
}



// EwaldInterfaceBase.h 150 - 169
template<class PA, class GPA>
void convertForce(PA& particlearray,
		  const std::vector<ParticleRange>& selfrange,
		  GPA& ghost,
		  const std::vector<ParticleRange>& ghostrange)
{
  std::vector<Force>::const_iterator itfc = fc.begin();
  for (std::vector<ParticleRange>::size_type it=0;it < selfrange.size();
       ++it) {
    for (int i = selfrange[it].begin; i < selfrange[it].end; ++i,++itfc) {
      getforce(particlearray,i) += *itfc;
    }
  }
  for (std::vector<ParticleRange>::size_type it=0;it < ghostrange.size();
       ++it) {
    for (int i = ghostrange[it].begin; i < ghostrange[it].end; ++i,++itfc) {
      getforce(ghost,i) += *itfc;
    }
  }
}

/*
// PMELongRangeInteraction.cc 63 - 75
// EwaldInterfaceBase.h 91 - 103
*/
template<class PA, class GPA>
void backinterpolate(PA& particlearray,
		     const std::vector<ParticleRange>& selfrange,
		     GPA& ghost,
		     const std::vector<ParticleRange>& ghostrange)
{
  if (!isEwaldBaseNode) return;
  fc.assign(cd.size(), Force());
  backInterpolateMain();
  convertForce(particlearray, selfrange, ghost, ghostrange);
}

// PMEMethod.h 36
void communicateForces() {}
// PMEMethod.h 37
void communicateForces(SpaceVector<double>& dipoleMoment) {}
// EwaldBase.h 27
int getSurfaceDipole() { return calcSurfaceDipole;}
// PMEMethod.h 45 - 47
const std::vector<Position>& getPositions() 
{
  return cdPME;
}
// PMEMethod.h 48 - 50
const std::vector<double>& getCharges()
{
  return chargePME;
}
// PMEMethod.h 51 - 53
std::vector<Force>& getForce() 
{
  return fcPME;
}



// PMEMethod.cc 12 - 43
void initialize_spread()
{
  if (order <= 0) {
    throw std::runtime_error(errorPos("PME order must be positive"));
  }
  if (order % 2 > 0) {
    throw std::runtime_error(errorPos("PME order must be even"));
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

void initialize_spread_t()
{
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
}


// EwaldBase.h 56
void setShare(double _shareValue) { shareValue = _shareValue; }

/*
// PMEMethod.h 150 - 152
// ParallelPMEMethod.cc 174 - 189
*/
void setCommunicator()
{
  comms.communicator = myworld;
  FFT3D::setCommunicator(&comms);
  if(myworld==MPI_COMM_NULL){
    myrank=0;
    nodeNum=0;
  }else{
    MPI_Comm_rank(myworld, &myrank);
    MPI_Comm_size(myworld, &nodeNum);
  }
  if(DebugLog::verbose>1){
    std::cout << " FFT Comm size rank " << nodeNum << " " << myrank << std::endl;
  }
}

// ParallelPMEMethod.cc  244 - 252
void setAxis(const SpaceVector<std::size_t>& axis)
{
  for (int i=0;i<SpaceVector<std::size_t>::Dim;++i) {
    if (axis[i] == 0) {
      divideAxis = i;
      break;
    }
  }
}

// PMEInterface.h  85
bool isPMENode() { return inPMEWorld; }


// ParallelPMEMethod.cc 94 - 96
int getRankFromBitIndex(int i)
{
  return borderNodes[i].nodeid;
}

// ParallelPMEMethod.cc 97 - 101
void initKeepSends(std::vector<KeepSend>& keepSends) 
{
  KeepSend ks;
  ks.reset();
  keepSends.assign(gridNum, ks);
}

// arallelPMEMethod.cc 102 - 112
void setKeepSends(std::vector<KeepSend>& keepSends, int start, int end, KeepSendFlag flag) 
{
  for (int i=start;i<end;++i) {
    if (borderNodes[i].width > 0) {
      for (int grid=borderNodes[i].border-overlap;
	   grid<borderNodes[i+1].border+overlap-1;++grid) {
	keepSends[(grid+gridNum)%gridNum].set(flag);
      }
    }
  }
}

// ParallelPMEMethod.cc  254 - 313
void calcSendPlan()
{
  int bitMax = 1;
  while (2*bitMax < nodeNum) bitMax *= 2;
  for (int bit=bitMax,tag=PME_TAG;bit > 0;bit /= 2,++tag) {
    SendPlan plan;
    plan.sendFirst = ((bitIndex & bit) != 0);
    plan.bitIndex = bitIndex ^ bit; /* xor */
    plan.tag = tag;
    int base = (bitIndex/bit)*bit;
    int basebit = base+bit;
    int planbase = (plan.bitIndex/bit)*bit;
    int planbasebit = planbase+bit;
    if (plan.bitIndex < nodeNum) {
      plan.keepAll = false;
      plan.sendOnly = false;
      plan.recvOnly = false;
      plan.destNode = getRankFromBitIndex(plan.bitIndex);
      initKeepSends(plan.keepSends);
      int sendStart = planbase;
      int sendEnd = std::min(planbasebit, nodeNum);
      setKeepSends(plan.keepSends, sendStart, sendEnd, SEND);
      int keepStart = base;
      int keepEnd = std::min(basebit, nodeNum);
      setKeepSends(plan.keepSends, keepStart, keepEnd, KEEP);
    }
    else if (nodeNum <= basebit) {
      plan.keepAll = true;
      plan.recvOnly = false;
    }
    else {
      plan.keepAll = false;
      plan.sendOnly = true;
      plan.recvOnly = false;
      plan.bitIndex = basebit + (plan.bitIndex-nodeNum)%(nodeNum-basebit);
      plan.destNode = getRankFromBitIndex(plan.bitIndex);
      initKeepSends(plan.keepSends);
      int sendStart = planbase;
      int sendEnd = std::min(planbasebit, nodeNum);
      setKeepSends(plan.keepSends, sendStart, sendEnd, SEND);
      int keepStart = base;
      int keepEnd = std::min(basebit, nodeNum);
      setKeepSends(plan.keepSends, keepStart, keepEnd, KEEP);
    }
    sendPlans.push_back(plan);

    if (bitIndex >= planbasebit) {
      plan.keepAll = true;
      plan.sendOnly = false;
      plan.recvOnly = true;
      for (plan.bitIndex = nodeNum+bitIndex-2*planbasebit;
           plan.bitIndex < planbasebit;plan.bitIndex += nodeNum-planbasebit) {
        if (plan.bitIndex+bit < nodeNum) continue;
        plan.destNode = getRankFromBitIndex(plan.bitIndex);
        initKeepSends(plan.keepSends);
        sendPlans.push_back(plan);
      }
    }
  }
}
// ParallelPMEMethod.cc  315 - 367
void calcSendBackPlan()
{
  int bitMax = 1;
  while (2*bitMax < nodeNum) bitMax *= 2;
  for (int bit=1,tag=PME_BACK_TAG,mask=1;bit <= bitMax;bit *= 2,++tag,
       mask |= bit) {
    SendPlan plan;
    plan.sendFirst = ((bitIndex & bit) != 0);
    plan.bitIndex = bitIndex ^ bit; /* xor */
    plan.tag = tag;
    plan.recvOnly = false;
    plan.sendOnly = false;
    for (int destBitIndex=0;destBitIndex<nodeNum;++destBitIndex) {
      KeepSend ks;
      ks.reset();
      if ((bitIndex & bit) == (destBitIndex & bit)) {
        ks.set(KEEP);
      }
      else {
        ks.set(SEND);
      }
      plan.keepSends.push_back(ks);
    }
    if (plan.bitIndex < nodeNum) {
      plan.keepAll = false;
      plan.destNode = getRankFromBitIndex(plan.bitIndex);
      sendBackPlans.push_back(plan);
    }
    for (int srcBitIndex=0;srcBitIndex<nodeNum;++srcBitIndex) {
      int destBitIndex = srcBitIndex ^ bit; /* xor */
      if (destBitIndex >= nodeNum) {
        destBitIndex &= mask;
        assert(srcBitIndex != destBitIndex);
        if (destBitIndex < nodeNum) {
          if (srcBitIndex == bitIndex) {
            plan.bitIndex = destBitIndex;
            plan.destNode = getRankFromBitIndex(plan.bitIndex);
            plan.keepAll = false;
            plan.sendOnly = true;
            sendBackPlans.push_back(plan);
          }
          else if (destBitIndex == bitIndex) {
            plan.bitIndex = srcBitIndex;
            plan.destNode = getRankFromBitIndex(plan.bitIndex);
            plan.keepAll = true;
            plan.recvOnly = true;
            sendBackPlans.push_back(plan);
          }
        }
      }
    }
  }
}

//  ParallelPMEMethod.cc 191 - 237
void setFFT(const FFT3D::Ptr& fft) 
{
  if (!isPMENode()) return;
  if (fft->getRealDivideType() != FFT3D::DIVIDE_X) {
    throw(std::runtime_error(errorPos("divide type not support")));
  }
  sevec.resize(nodeNum);
  StartEnd se;
  se.start = fft->getRealStart(0);
  se.end = fft->getRealEnd(0);
  MPI_Allgather(&se, 2, MPI_INTEGER, &sevec[0], 2, MPI_INTEGER, myworld);
  int endmax = 0;
  for (int i=0;i<nodeNum;++i) {
    endmax = std::max(endmax, sevec[i].end);
  }
  gridNum = fft->getSize(0);
  if (gridNum != endmax) {
    throw(std::runtime_error(errorPos("fft->getSize(0) illegal")));
  }
  borderNodes.clear();
  BorderNode bn;
  for (int i=0;i<nodeNum;++i) {
    bn.border = sevec[i].start;
    bn.nodeid = i;
    bn.width = sevec[i].end-sevec[i].start;
    if (bn.width > 0) {
      bn.width = 1;
    }
    else {
      bn.border = gridNum;
    }
    borderNodes.push_back(bn);
  }
  std::sort(borderNodes.begin(), borderNodes.end(), lessBorderNode());
  bn.border = gridNum;
  bn.nodeid = nodeNum;
  bn.width = 1;
  borderNodes.push_back(bn);
  for (int i=0;i<=nodeNum;++i) {
    if (borderNodes[i].nodeid == myrank) {
      bitIndex = i;
      break;
    }
  }
  calcSendPlan();
  calcSendBackPlan();
}

// PMEMethod.cc 98 - 111
void calcb2(std::vector<double>& bv, int length, BSpline::Ptr& bsp)
{
  const std::complex<double> factor(0.0, 2.0*M_PI/length);
  bsp->calculate(0.0);
  const double* const M = bsp->getFunction();
  for (int m=0;m<length;++m) {
    std::complex<double> s = 0.0;
    for (int k=0;k<order-1;++k) {
      s += M[k+1]*exp(factor*double(m*k));
    }
    bv.push_back(norm(exp(factor*double((order-1)*m))/s));
  }
}

//  PMEMethod.h 186 - 198
void initializePMEArray()
{
  if (pmeType == SmoothPME) {
    std::vector<double> b2X, b2Y, b2Z;
    calcb2(b2X, nGrid.getX(), bspX);
    calcb2(b2Y, nGrid.getY(), bspY);
    calcb2(b2Z, nGrid.getZ(), bspZ);
    pmeArray = PMEArray::Ptr(new PMEArray(fft, alpha, axis, b2X, b2Y, b2Z));
  }
  else {
    pmeArray = PMEArray::Ptr(new PMEArray(fft, alpha, axis));
  }
}


//// start,end
// PMEMethod.cc 155-227
void initialize_method()
{
  for (int i=0;i<SpaceVector<std::size_t>::Dim;++i) {
    //    side[i] = pmeInterface->getSide(pContext)[axis[i]];
    side[i] = param.boxSize[axis[i]];
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
  //  parallel->setSide(pmeInterface->getSide(pContext));
  setAxis(axis);
  if(DebugLog::verbose>2){
    if (pmeType == SmoothPME) {
      std::cout << "SmoothPME " << order << " " << nGrid << std::endl;
    }
    else {
      std::cout << "PME " << order << " " << nGrid << std::endl;
    }
  }
  setFFT(fft);
  initializePMEArray();
  //  pmeArray->calculate(pmeInterface, pContext);
  pmeArray->calculate(param.boxSize);

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
          throw std::runtime_error(errorPos("PMEMethod pme_order grid_length unbalance. use smaller pme_order or smaller grid_length."));
        }
      }
      start[i].push_back(sinx);
      end[i].push_back(einx);
    }
  }
}

// constractor
void constructPME(const LongRangeParameter& _param, MPI_Comm lcomm=MPI_COMM_WORLD)
{
  param = _param;
  selfEnergyCalculated = false;
  isEwaldBaseNode = true;
  /*    
// PMEInterface.h 22 - 72
  */
  myworld = lcomm;
  inPMEWorld = true;

  {
    if(myworld==MPI_COMM_NULL){
      inPMEWorld=false;
      return;
    }else{
      // debug message
    }

    int longrank=0;
    MPI_Comm_rank(myworld,&longrank);
    if (longrank == 0) 
      {
      if (param.pmeType == SmoothPME) {
	std::cout << "SmoothPME ";
      }else{
	std::cout << "PME ";
      }
      std::cout << "order " << param.order <<" Grid " << param.grid_num << std::endl;
    }

    {
    /*
// EwaldBase.h 10 - 20
*/
      {
	cutoff = param.cutoff;
	alpha = param.alpha;
	kCutoff = param.kCutoff;
	calcSurfaceDipole = param.surfaceDipole;
	// surface()    // construct surface
	selfEnergy = 0.0;
	shareValue = 1.0;
	erfcCutValue = 1.0e-6;
	ewaldExpCutValue = 5.0e-9;
      }
      /*
// PMEMethod.h 89 - 104
      */
      {
	pmeType = param.pmeType;
	order = param.order;
	/* initialized at initialize()
	   lagX(), lagY(), lagZ(), bspX(), bspY(), bspZ(),
	   funcX(), funcY(), funcZ(), dfuncX(), dfuncY(), dfuncZ(),
	 */
	/* initialized at initialize(&pcontext)
	   nGrid(), indexes(), pmeArray(), start(), end(), 
	 */
	axis[0]=0;axis[1]=1;axis[2]=2;
	prepared = false;
	// construct MPIPMEMethod parallel
	/*
//  ParallelPMEMethod.cc 161 - 171
// PMEMethod.h 19 - 21
// ParallelPMEMethod.cc 18 - 29
	 */
	{
	  overlap = order/2;
	  withDipoleMoment = false;
	}
	initialize_spread();
      }
    }
    if (longrank != 0) {
      setShare(0.0);
    }
    isEwaldBaseNode = inPMEWorld;
    setCommunicator();
    if(myworld!=MPI_COMM_NULL){
      MPI_Barrier(myworld);
    }
    fft = FFT3D::createFFT3D(param.grid_num[0],
			     param.grid_num[1],
			     param.grid_num[2]);
    if(DebugLog::verbose>1)std::cout << " createdFFT3D " <<std::endl;
    initialize_method();  // method->initialize(&context);
    if(DebugLog::verbose>1)std::cout << " PMEMethod initialize " <<std::endl;
  }
}

// PMEInterface.h 88
MPI_Comm getPMECommunicator() { return myworld; }

// PMEInterface.h 86
FFT3D::Ptr getFFT() { return fft; }


double getSelfEnergy() { return selfEnergy; }
double getDipoleEnergy() { return shareValue*surface.dipoleEnergy(); }


}
//// nGrid 64 , node 0 [0:15]
/// indexes
/*
     0  1  2  3  4  5 ...  50 51 52 53 54 55 56 57 58 59 60 61 62 63 63 65 66 67
     2  1  0 63 62 61 ...  16 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1  0 63


 */


#include "PME.h"

template<class PA, class GPA>
void PMEChargeAssign::assign(PA& particlearray,
                             const std::vector<ParticleRange>& selfrange,
                             GPA& ghost,
                             const std::vector<ParticleRange>& ghostrange,
                             GridData& gridcharge,
                             const std::vector<ParticleRange>& self_selfenergy_range,
                             const std::vector<ParticleRange>& ghost_selfenergy_range)
{
  PMEModule::assign(particlearray,selfrange,ghost,ghostrange,self_selfenergy_range,ghost_selfenergy_range);
}
#ifdef OLDPARTICLE
template
void PMEChargeAssign::assign(ParticleArray& particlearray,
                             const std::vector<ParticleRange>& selfrange,
                             ParticleArray& ghost,
                             const std::vector<ParticleRange>& ghostrange,
                             GridData& gridcharge,
                             const std::vector<ParticleRange>& self_selfenergy_range,
                             const std::vector<ParticleRange>& ghost_selfenergy_range);
#else
template
void PMEChargeAssign::assign(CombinedParticleArray& particlearray,
                             const std::vector<ParticleRange>& selfrange,
                             GhostParticleArray& ghost,
                             const std::vector<ParticleRange>& ghostrange,
                             GridData& gridcharge,
                             const std::vector<ParticleRange>& self_selfenergy_range,
                             const std::vector<ParticleRange>& ghost_selfenergy_range);
#endif

template<class PA, class GPA>
void PMEChargeAssign::backinterpolate
                      (PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential)
{
  PMEModule::backinterpolate(particlearray, selfrange, ghost, ghostrange);
}
#ifdef OLDPARTICLE
template
void PMEChargeAssign::backinterpolate
                      (ParticleArray& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       ParticleArray& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
#else
template
void PMEChargeAssign::backinterpolate
                      (CombinedParticleArray& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GhostParticleArray& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
#endif

template<class PA, class GPA>
void PMEChargeAssign::addenergy(PA& particlearray,
                                const std::vector<ParticleRange>& selfrange,
                                GPA& ghost,
                                const std::vector<ParticleRange>& ghostrange,
                                double& energy)
{
  energy += PMEModule::getSelfEnergy();
  energy += PMEModule::getDipoleEnergy();
}
#ifdef OLDPARTICLE
template
void PMEChargeAssign::addenergy(ParticleArray& particlearray,
                                const std::vector<ParticleRange>& selfrange,
                                ParticleArray& ghost,
                                const std::vector<ParticleRange>& ghostrange,
                                double& energy);
#else
template
void PMEChargeAssign::addenergy(CombinedParticleArray& particlearray,
                                const std::vector<ParticleRange>& selfrange,
                                GhostParticleArray& ghost,
                                const std::vector<ParticleRange>& ghostrange,
                                double& energy);
#endif

PMEPoissonSolver::PMEPoissonSolver(int unitid,
                                   const LongRangeParameter& _param,
                                   PMEModuleInterface* pmemi) 
{
}

void PMEPoissonSolver::solvePoisson(GridData& gridcharge,
                                    GridData& gridpotential, double& energy)
{
  //  std::cout << "PMEChargeAssign::solvePoisson" << std::endl;
  PMEModule::solvePoisson(energy);
}

void PMEPoissonSolver::solvePoisson(GridData& gridcharge,
                                    GridData& gridpotential, double& energy, double &virial)
{
  //  std::cout << "PMEChargeAssign::solvePoisson" << std::endl;
  PMEModule::solvePoisson(energy,virial);
}


template<>
void PMELongRangeInteraction::initialize()
{
  PMEModule::constructPME(chargeassign.pmemoduleinterface.param, chargeassign.pmemoduleinterface.myworld);
  PMEModule::initialize_spread_t();
}
