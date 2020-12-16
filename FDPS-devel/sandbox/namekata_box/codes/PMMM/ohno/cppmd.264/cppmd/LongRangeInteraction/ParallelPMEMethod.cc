#include <mpi.h>
#include "MPIFFT3D.h"
#include "PMEInterface.h"
#include "PMEMethod.h"
#include <bitset>
#include "NullStream.h"
#include "ErrorPos.h"

using namespace PMEModule;

#define PACKSEND 1


#ifdef USE_MPI
/* for DIVIDE_X code */
class MPIPMEMethod: public ParallelPMEMethod {
 public:
  MPIPMEMethod(PMEInterface* _pmeInterface, int order):
      ParallelPMEMethod(),
      pmeInterface(_pmeInterface),
      myworld(), myrank(), nodeNum(), gridNum(),
      sevec(), side(), borderNodes(), bitIndex(),
      particleIndex(), particleNum(),
      overlap(order/2), divideAxis(),
      cdSend(), chargeSend(), bitGridSend(),
      cdKeep(), chargeKeep(), bitGridKeep(), bitGridPME(),
      fcKeep(), fcSend(),
      withDipoleMoment(false), dipoleMoment(), recvDipoleMoment(),
      logfs() {}

  void setCommunicator();
  void setFFT(const FFT3D::Ptr& fft);
  void setSide(const SpaceVector<double>& side);
  void setAxis(const SpaceVector<std::size_t>& axis);
#ifdef USE_MODULE_LS_COMM
  void clear();
  void setPositionsAndCharges(const std::vector<Position>& cd,
                              const std::vector<double>& charge);
  void communicatePositionsAndCharges();
  void communicateForces();
  void communicateForces(SpaceVector<double>& _dipoleMoment);
#endif
  void setLog(std::ostream& _logfs) {
    logfs.rdbuf(_logfs.rdbuf());
  }
  std::ostream& getLog() { return logfs; }
  double share();
 private:
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
  int getRankFromBitIndex(int i) {
    return borderNodes[i].nodeid;
  }
  void initKeepSends(std::vector<KeepSend>& keepSends) {
    KeepSend ks;
    ks.reset();
    keepSends.assign(gridNum, ks);
  }
  void setKeepSends(std::vector<KeepSend>& keepSends, int start, int end,
                    KeepSendFlag flag) {
    for (int i=start;i<end;++i) {
      if (borderNodes[i].width > 0) {
        for (int grid=borderNodes[i].border-overlap;
             grid<borderNodes[i+1].border+overlap-1;++grid) {
          keepSends[(grid+gridNum)%gridNum].set(flag);
        }
      }
    }
  }
  void calcSendPlan();
  void calcSendBackPlan();
  void sendPositionsAndCharges(int destNode, int tag);
  void recvPositionsAndCharges(int destNode, int tag);
  void sendForces(int destNode, int tag);
  void recvForces(int destNode, int tag);
#ifdef USE_MODULE_LS_COMM
  void setPositionsAndCharges(int planIndex);
#endif
  void setForces(int planIndex);

  static const int PME_TAG = 100;
  static const int PME_BACK_TAG = 150;

  PMEInterface* pmeInterface;
  MPI_Comm myworld;
  int myrank;
  int nodeNum;
  int gridNum;
  std::vector<StartEnd> sevec;
  SpaceVector<double> side;
  std::vector<BorderNode> borderNodes;
  int bitIndex;
  std::size_t particleIndex;
  std::size_t particleNum;
  std::vector<SendPlan> sendPlans;
  std::vector<SendPlan> sendBackPlans;
  int overlap;
  int divideAxis;
  std::vector<Position>cdSend;
  std::vector<double>chargeSend;
  std::vector<BitGrid>bitGridSend;
  std::vector<Position>cdKeep;
  std::vector<double>chargeKeep;
  std::vector<BitGrid>bitGridKeep;
  std::vector<BitGrid>bitGridPME;
  std::vector<double>sendBuf;
  std::vector<double>recvBuf;
  std::vector<Force>fcKeep;
  std::vector<Force>fcSend;
  bool withDipoleMoment;
  SpaceVector<double> dipoleMoment;
  SpaceVector<double> recvDipoleMoment;
  Communicator comms;
  NullStream logfs;
};
#endif

ParallelPMEMethod::Ptr PMEModule::createParallelPMEMethod
(PMEInterface* pmeInterface, int order)
{
#ifdef USE_MPI
  if(DebugLog::verbose>1)std::cout << "PMEModule::createParallelPMEMethod use MPI" << std::endl;
  return ParallelPMEMethod::Ptr(new MPIPMEMethod(pmeInterface, order));
#else
  if(DebugLog::verbose>1)std::cout << "PMEModule::createParallelPMEMethod no MPI" << std::endl;
  return ParallelPMEMethod::Ptr(new ParallelPMEMethod());
#endif
}

#ifdef USE_MPI
void MPIPMEMethod::setCommunicator()
{
  myworld = pmeInterface->getPMECommunicator();
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

void MPIPMEMethod::setFFT(const FFT3D::Ptr& fft)
{
  if (!pmeInterface->isPMENode()) return;
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

void MPIPMEMethod::setSide(const SpaceVector<double>& _side)
{
  side = _side;
}

void MPIPMEMethod::setAxis(const SpaceVector<std::size_t>& axis)
{
  for (int i=0;i<SpaceVector<std::size_t>::Dim;++i) {
    if (axis[i] == 0) {
      divideAxis = i;
      break;
    }
  }
}

void MPIPMEMethod::calcSendPlan()
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

void MPIPMEMethod::calcSendBackPlan()
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

#ifdef USE_MODULE_LS_COMM
void MPIPMEMethod::clear()
{
  cdSend.clear();
  chargeSend.clear();
  bitGridSend.clear();
  cdKeep.clear();
  chargeKeep.clear();
  bitGridKeep.clear();
  particleIndex = 0;
}

void MPIPMEMethod::setPositionsAndCharges(const std::vector<Position>& cd,
                                          const std::vector<double>& charge)
{
  if (!pmeInterface->isPMENode()) return;
  assert(cd.size() == charge.size());
  assert(sendPlans.size() > 0);
  double factor = gridNum/side[divideAxis];
  SendPlan& plan = sendPlans[0];
  for (int i=0;i<cd.size();++i,++particleIndex) {
    BitGrid bg;
    bg.bg.bitIndex = bitIndex;
    bg.bg.gridIndex = (int(floor(factor*cd[i][divideAxis]))+gridNum) % gridNum;
    bg.bg.particleIndex = particleIndex;
    KeepSend ks = plan.keepSends[bg.bg.gridIndex];
    if (plan.keepAll || ks[KEEP]) {
      cdKeep.push_back(cd[i]);
      chargeKeep.push_back(charge[i]);
      bitGridKeep.push_back(bg);
    }
    if (!plan.keepAll && ks[SEND]) {
      cdSend.push_back(cd[i]);
      chargeSend.push_back(charge[i]);
      bitGridSend.push_back(bg);
    }
  }
}

void MPIPMEMethod::setPositionsAndCharges(int planIndex)
{
  if (!pmeInterface->isPMENode()) return;
  assert(cdPME.size() == chargePME.size());
  assert(cdPME.size() == bitGridPME.size());
  assert(sendPlans.size() > 0);
  SendPlan& plan = sendPlans[planIndex];
  for (int i=0;i<cdPME.size();++i) {
    KeepSend ks = plan.keepSends[bitGridPME[i].bg.gridIndex];
    if (plan.keepAll || ks[KEEP]) {
      cdKeep.push_back(cdPME[i]);
      chargeKeep.push_back(chargePME[i]);
      bitGridKeep.push_back(bitGridPME[i]);
    }
    if (!plan.keepAll && ks[SEND]) {
      cdSend.push_back(cdPME[i]);
      chargeSend.push_back(chargePME[i]);
      bitGridSend.push_back(bitGridPME[i]);
    }
    if (!plan.keepAll && ks.none()) {
      throw(std::runtime_error(errorPos("logic error")));
    }
  }
}
#endif

void MPIPMEMethod::setForces(int planIndex)
{
  if (!pmeInterface->isPMENode()) return;
  assert(fcPME.size() == bitGridPME.size());
  fcKeep.clear();
  fcSend.clear();
  bitGridKeep.clear();
  bitGridSend.clear();
  SendPlan& plan = sendBackPlans[planIndex];
  for (std::vector<Force>::size_type i=0;i<fcPME.size();++i) {
    KeepSend ks = plan.keepSends[bitGridPME[i].bg.bitIndex];
    if (ks[SEND]) {
      fcSend.push_back(fcPME[i]);
      bitGridSend.push_back(bitGridPME[i]);
    } 
    if (ks[KEEP]) {
      fcKeep.push_back(fcPME[i]);
      bitGridKeep.push_back(bitGridPME[i]);
    }
  }
}

#ifdef USE_MODULE_LS_COMM
void MPIPMEMethod::communicatePositionsAndCharges()
{
  std::cout << "communicatePositionsAndCharges ";
  if (!pmeInterface->isPMENode()){
    std::cout << "not isPMENode()" << std::endl;
    return;
  }
  std::cout << "isPMENode()" << std::endl;
  particleNum = particleIndex;
  bool first = true;
  for (int planIndex=0;planIndex < sendPlans.size();++planIndex) {
    if (first) {
      first = false;
    }
    else {
      clear();
      setPositionsAndCharges(planIndex);
    }
    SendPlan& plan = sendPlans[planIndex];
    cdPME.swap(cdKeep);
    chargePME.swap(chargeKeep);
    bitGridPME.swap(bitGridKeep);
    if (!plan.keepAll) {
      if (plan.recvOnly) {
        recvPositionsAndCharges(plan.destNode, plan.tag);
      }
      else if (plan.sendOnly) {
        sendPositionsAndCharges(plan.destNode, plan.tag);
      }
      else if (plan.sendFirst) {
        sendPositionsAndCharges(plan.destNode, plan.tag);
        recvPositionsAndCharges(plan.destNode, plan.tag);
      }
      else {
        recvPositionsAndCharges(plan.destNode, plan.tag);
        sendPositionsAndCharges(plan.destNode, plan.tag);
      }
    }
    else {
      if (plan.recvOnly) {
        recvPositionsAndCharges(plan.destNode, plan.tag);
      }
    }
  }
  if (cdPME.size() > 0) {
    double minZ = cdPME[0][2];
    double maxZ = cdPME[0][2];
    for (int i=0;i<cdPME.size();++i) {
      minZ = std::min(cdPME[i][2], minZ);
      maxZ = std::max(cdPME[i][2], maxZ);
    }
  }
}
#endif

void MPIPMEMethod::sendPositionsAndCharges(int destNode, int tag)
{
  std::cout << " MPIPMEMethod::sendPositionsAndCharges" << std::endl;
#if PACKSEND
  int cdsize = cdSend.size()*(Position::Dim);
  int bgsize = (sizeof(BitGrid)*bitGridSend.size()+sizeof(double)-1)
      /sizeof(double);
  int datasize = cdsize + chargeSend.size() + bgsize;

  sendBuf.resize(datasize);
  std::vector<double>::iterator it = sendBuf.begin();
  it = std::copy(&(cdSend[0][0]), &(cdSend[0][0])+cdsize, it);
  it = std::copy(chargeSend.begin(), chargeSend.end(), it);
  it = std::copy(&(bitGridSend[0].data), (&(bitGridSend[0].data))+bgsize, it);
  assert(sendBuf.size() == static_cast<std::vector<double>::size_type>(datasize));
  assert(cdSend.size() == chargeSend.size());
  assert(cdSend.size() == bitGridSend.size());
  int sizes[2] = { datasize, cdSend.size() };
  MPI_Send(sizes, 2, MPI_INTEGER, destNode, tag, myworld);
  MPI_Send(&sendBuf[0], sendBuf.size(), MPI_DOUBLE, destNode, tag, myworld);
#else
  int datasize = cdSend.size();
  assert(cdSend.size() == chargeSend.size());
  assert(cdSend.size() == bitGridSend.size());
  assert(sizeof(cdSend[0]) == sizeof(double)*3);
  assert(sizeof(bitGridSend[0]) == sizeof(int)*3);
  MPI_Send(&datasize, 1, MPI_INTEGER, destNode, tag, myworld);
  MPI_Send(&cdSend[0], datasize*3, MPI_DOUBLE, destNode, tag, myworld);
  MPI_Send(&chargeSend[0], datasize, MPI_DOUBLE, destNode, tag, myworld);
  MPI_Send(&bitGridSend[0], datasize*3, MPI_INTEGER, destNode, tag, myworld);
#endif
}

void MPIPMEMethod::recvPositionsAndCharges(int destNode, int tag)
{
  int sizes[2];
  MPI_Status status;
  int count;

  std::cout << "MPIPMEMethod::recvPositionsAndCharges" << std::endl;
#if PACKSEND
  MPI_Recv(sizes, 2, MPI_INTEGER, destNode, tag, myworld, &status);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  assert(count==2);
  recvBuf.resize(sizes[0]);
  MPI_Recv(&recvBuf[0], sizes[0], MPI_DOUBLE, destNode, tag, myworld, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  assert(count==sizes[0]);

  int cdsize = sizes[1]*(Position::Dim);
  int bgsize = (sizeof(BitGrid)*sizes[1]+sizeof(double)-1)
      /sizeof(double);
  int databgsize = (bgsize*sizeof(double)+sizeof(BitGrid)-1)/sizeof(BitGrid);

  assert(cdPME.size() == chargePME.size());
  assert(cdPME.size() == bitGridPME.size());
  int currentSize = cdPME.size();
  int afterSize = currentSize + sizes[1];
  int afterBgSize = currentSize + databgsize;
  cdPME.resize(afterSize);
  double *p = &recvBuf[0];
  std::copy(p, p+cdsize, &(cdPME[currentSize][0]));
  p += cdsize;

  chargePME.resize(afterSize);
  std::copy(p, p+sizes[1], &chargePME[currentSize]);
  p += sizes[1];

  bitGridPME.resize(afterBgSize);
  std::copy(p, p+bgsize, &(bitGridPME[currentSize].data));
  bitGridPME.resize(afterSize);
#else
  int datasize;
  MPI_Recv(&datasize, 1, MPI_INTEGER, destNode, tag, myworld, &status);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  assert(count==1);
  assert(cdPME.size() == chargePME.size());
  assert(cdPME.size() == bitGridPME.size());
  assert(sizeof(cdPME[0]) == sizeof(double)*3);
  assert(sizeof(bitGridPME[0]) == sizeof(int)*3);
  int currentSize = cdPME.size();
  cdPME.resize(currentSize+datasize);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  MPI_Recv(&(cdPME[currentSize][0]), datasize*3, MPI_DOUBLE, destNode, tag,
           myworld, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  assert(count==datasize*3);
  chargePME.resize(currentSize+datasize);
  MPI_Recv(&chargePME[currentSize], datasize, MPI_DOUBLE, destNode, tag,
           myworld, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  assert(count==datasize);
  bitGridPME.resize(currentSize+datasize);
  MPI_Recv(&bitGridPME[currentSize], datasize*3, MPI_INTEGER, destNode, tag,
           myworld, &status);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  assert(count==datasize*3);
#endif
}

#ifdef USE_MODULE_LS_COMM
void MPIPMEMethod::communicateForces()
{
  if (!pmeInterface->isPMENode()) return;
  for (int planIndex=0;planIndex < sendBackPlans.size();++planIndex) {
    SendPlan& plan = sendBackPlans[planIndex];
    if (withDipoleMoment) {
      recvDipoleMoment = 0.0;
    }
    setForces(planIndex);
    fcPME.swap(fcKeep);
    bitGridPME.swap(bitGridKeep);
    if (!plan.keepAll) {
      if (plan.sendOnly) {
        sendForces(plan.destNode, plan.tag);
      }
      else if (plan.sendFirst) {
        sendForces(plan.destNode, plan.tag);
        recvForces(plan.destNode, plan.tag);
      }
      else {
        recvForces(plan.destNode, plan.tag);
        sendForces(plan.destNode, plan.tag);
      }
    }
    else if (plan.recvOnly) {
      recvForces(plan.destNode, plan.tag);
    }
    if (withDipoleMoment) {
      dipoleMoment += recvDipoleMoment;
    }
  }
  if (nodeNum > 1) {
    fcKeep.swap(fcPME);
    fcPME.assign(particleNum, Force());
    for (int i=0;i<fcKeep.size();++i) {
      fcPME.at(bitGridPME[i].particleIndex) += fcKeep[i];
    }
  }
  itfcPME = fcPME.begin();
}

void MPIPMEMethod::communicateForces(SpaceVector<double>& _dipoleMoment)
{
  withDipoleMoment = true;
  dipoleMoment = _dipoleMoment;
  communicateForces();
  _dipoleMoment = dipoleMoment;
}
#endif

void MPIPMEMethod::sendForces(int destNode, int tag)
{
#if PACKSEND
  int fcsize = fcSend.size()*(Force::Dim);
  int bgsize = (sizeof(BitGrid)*bitGridSend.size()+sizeof(double)-1)
      /sizeof(double);
  int datasize = fcsize + bgsize;

  if (withDipoleMoment) {
    datasize += 3;
  }
  sendBuf.resize(datasize);
  std::vector<double>::iterator it = sendBuf.begin();
  it = std::copy(&(fcSend[0][0]), &(fcSend[0][0])+fcsize, it);
  it = std::copy(&(bitGridSend[0].data), (&(bitGridSend[0].data))+bgsize, it);
  if (withDipoleMoment) {
    it = std::copy(&dipoleMoment[0], &dipoleMoment[0]+3, it);
  }
  assert(sendBuf.size() == static_cast<std::vector<double>::size_type>(datasize));
  assert(fcSend.size() == bitGridSend.size());
  int sizes[2] = { datasize, fcSend.size() };
  MPI_Send(sizes, 2, MPI_INTEGER, destNode, tag, myworld);
  MPI_Send(&sendBuf[0], sendBuf.size(), MPI_DOUBLE, destNode, tag, myworld);
#else
  int datasize = fcSend.size();
  assert(fcSend.size() == bitGridSend.size());
  assert(sizeof(fcSend[0]) == sizeof(double)*3);
  assert(sizeof(bitGridSend[0]) == sizeof(int)*3);
  MPI_Send(&datasize, 1, MPI_INTEGER, destNode, tag, myworld);
  MPI_Send(&fcSend[0], datasize*3, MPI_DOUBLE, destNode, tag, myworld);
  MPI_Send(&bitGridSend[0], datasize*3, MPI_INTEGER, destNode, tag, myworld);
#endif
}

void MPIPMEMethod::recvForces(int destNode, int tag)
{
  MPI_Status status;
  int count;
#if PACKSEND
  int sizes[2];
  MPI_Recv(sizes, 2, MPI_INTEGER, destNode, tag, myworld, &status);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  assert(count==2);
  recvBuf.resize(sizes[0]);
  MPI_Recv(&recvBuf[0], sizes[0], MPI_DOUBLE, destNode, tag, myworld, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  assert(count==sizes[0]);

  int fcsize = sizes[1]*(Force::Dim);
  int bgsize = (sizeof(BitGrid)*sizes[1]+sizeof(double)-1)
      /sizeof(double);
  int databgsize = (bgsize*sizeof(double)+sizeof(BitGrid)-1)/sizeof(BitGrid);

  assert(fcPME.size() == bitGridPME.size());
  int currentSize = fcPME.size();
  int afterSize = currentSize + sizes[1];
  int afterBgSize = currentSize + databgsize;
  fcPME.resize(afterSize);
  double *p = &recvBuf[0];
  std::copy(p, p+fcsize, &(fcPME[currentSize][0]));
  p += fcsize;

  bitGridPME.resize(afterBgSize);
  std::copy(p, p+bgsize, &(bitGridPME[currentSize].data));
  bitGridPME.resize(afterSize);
  if (withDipoleMoment) {
    p += bgsize;
    std::copy(p, p+3, &recvDipoleMoment[0]);
  }
#else
  int datasize;
  MPI_Recv(&datasize, 1, MPI_INTEGER, destNode, tag, myworld, &status);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  assert(count==1);
  assert(fcPME.size() == bitGridPME.size());
  assert(sizeof(fcPME[0]) == sizeof(double)*3);
  assert(sizeof(bitGridPME[0]) == sizeof(int)*3);
  int currentSize = fcPME.size();
  fcPME.resize(currentSize+datasize);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  MPI_Recv(&(fcPME[currentSize][0]), datasize*3, MPI_DOUBLE, destNode, tag,
           myworld, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  assert(count==datasize*3);
  bitGridPME.resize(currentSize+datasize);
  MPI_Recv(&bitGridPME[currentSize], datasize*3, MPI_INTEGER, destNode, tag,
           myworld, &status);
  MPI_Get_count(&status, MPI_INTEGER, &count);
  assert(count==datasize*3);
#endif
}

#endif
