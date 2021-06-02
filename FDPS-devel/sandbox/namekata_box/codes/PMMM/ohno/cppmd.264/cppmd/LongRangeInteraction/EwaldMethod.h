#ifndef BODYEWALDORGMETHOD_H
#define BODYEWALDORGMETHOD_H

#include <vector>
#include <cassert>
#include "EwaldInterfaceFwd.h"
#include "SurfaceDipole.h"
#include "EwaldBase.h"

namespace EwaldModule {

typedef struct {
  double s;
  double c;
} DFTResult;
typedef struct {
  int norm2;
  SpaceVector<int> v;
} WaveVector;

class EwaldMethod;

class ParallelEwald {
public:
  typedef ParallelEwald* Ptr;
  
  ParallelEwald() {}
  virtual ~ParallelEwald() {}

  virtual void aggregateDFTResult(
    std::vector<DFTResult>& dftr) {}
  virtual void deliverData(const std::vector<Position>& cd,
		   const std::vector<double>& charge) {}
  virtual void deliverDataComplete() {}
  virtual void aggregateForce(std::vector<Force>& fc) {}
  virtual void aggregate(double& potentialEnergy) {}
  virtual void calcDFT(EwaldInterface_::Context pContext,
                       EwaldMethod* p) {}
  virtual void calcIDFT(EwaldInterface_::Context pContext,
                        EwaldMethod* p) {}
};

ParallelEwald::Ptr createParallelEwald
                            (EwaldInterface* ewaldInterface);

class EwaldMethod : public EwaldBase {
public:
  EwaldMethod(EwaldInterface* _ewaldInterface,
                     const double cutoff,
                     const double alpha=0.0, const double kCutoff=0.0,
                     const int surfaceDipole=0)
    : EwaldBase(cutoff, alpha, kCutoff, surfaceDipole),
      ewaldInterface(_ewaldInterface),
      dftr(), kvec(), kstart(0), kend(0), prepared(false),
      parallel(createParallelEwald(_ewaldInterface)) {}

  virtual ~EwaldMethod() {
    delete(parallel);
  }

  void prepare(EwaldInterface_::Context pContext,
	       const std::vector<Position>& cd,
	       const std::vector<double>& charge) {
    if (!prepared) {
      initializeDFT();
      initializeKLoopRange();
      calculateSelfEnergy(charge);
      prepared = true;
    }
  }
  void calculate(EwaldInterface_::Context pContext,
	         const std::vector<Position>& cd,
	         const std::vector<double>& charge,
	         std::vector<Force>& fc,
                 double& potentialEnergy,
		 double& virial);

  void setSide(double side);
  void initializeDFT();
  void initializeKLoopRange();
  void prepareDFT(const std::vector<Position>& cd,
		  const std::vector<double>& charge) {
    parallel->deliverData(cd, charge);
  }
  void clearDFT();
  void calculateDFT(EwaldInterface_::Context pContext,
		    const std::vector<Position>& cd,
		    const std::vector<double>& charge);
  void completeDFT(EwaldInterface_::Context pContext);
  double calculateEwaldEnergy(EwaldInterface_::Context pContext,
			      double &virial);
  void calculateEwaldForce(EwaldInterface_::Context pContext,
		           const std::vector<Position>& cd,
		           const std::vector<double>& charge,
		           std::vector<Force>& fc);
  void complete(EwaldInterface_::Context pContext) {
    parallel->calcIDFT(pContext, this);
  }
  void completeForce(std::vector<Force>& fc) {
    parallel->aggregateForce(fc);
  }
private:
  void setWaveVector(const double kCutoff);
  void addPotentialEnergy(EwaldInterface_::Context pContext,
                          double& potentialEnergy);

  EwaldInterface* ewaldInterface;
  std::vector<DFTResult> dftr;
  std::vector<WaveVector> kvec;
  int kstart, kend;
  bool prepared;
  ParallelEwald::Ptr parallel;
};
}
#endif
