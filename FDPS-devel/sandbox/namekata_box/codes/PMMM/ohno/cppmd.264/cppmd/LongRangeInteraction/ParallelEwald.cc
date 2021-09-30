#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H

#include <mpi.h>
#include "Common.h"
#include "EwaldInterface.h"


using namespace EwaldModule;

#ifdef USE_MPI
class MPIEwald: public ParallelEwald {
    public:
	MPIEwald(EwaldInterface* _ewaldInterface):
            ewaldInterface(_ewaldInterface),
	    dftsendbuf(),
	    numbuf(0), cntbuf(0), cdbuf(), chargebuf(), flgbuf(),
            fcbuf(1), recvbuf(1) {}
	~MPIEwald() {}

	void finalize();
	void aggregateDFTResult(std::vector<DFTResult>& dftr);
	void deliverData(const std::vector<Position>& cd,
			 const std::vector<double>& charge);
	void deliverDataComplete();
	void aggregateForce(std::vector<Force>& fc);
	void aggregate(double& potentialEnergy);
	void calcDFT(EwaldInterface::Context pContext,
                     EwaldMethod* p);
	void calcIDFT(EwaldInterface::Context pContext,
                     EwaldMethod* p);
    private:
        EwaldInterface* ewaldInterface;
	std::vector<DFTResult> dftsendbuf;
	int numbuf;
	int cntbuf;
	std::vector<Position> cdbuf;
	std::vector<double> chargebuf;
  std::vector<unsigned int> flgbuf;
	std::vector<Force> fcbuf;
	std::vector<Force> recvbuf;
};
#endif

ParallelEwald::Ptr EwaldModule::createParallelEwald
                              (EwaldInterface *ewaldInterface)
{
#ifdef USE_MPI
    return ParallelEwald::Ptr(new MPIEwald(ewaldInterface));
#else
    return ParallelEwald::Ptr(new ParallelEwald());
#endif
}

#ifdef USE_MPI

void MPIEwald::aggregateDFTResult(std::vector<DFTResult>& dftr)
{
#define DFTR_N 2
    if (ewaldInterface->isEwaldNode()) {
	assert(sizeof(DFTResult)==sizeof(double)*DFTR_N);
	dftsendbuf.assign(dftr.begin(), dftr.end());
        /*
        int n;
        MPI_Comm_size(ewaldInterface->dftrCommunicator(),&n);
        printf("MPIEwald::aggregateDFTResult MPI_Allreduce %d\n",n);
        */
	MPI_Allreduce(&dftsendbuf[0], &dftr[0], dftsendbuf.size()*DFTR_N,
	    MPI_DOUBLE, MPI_SUM, ewaldInterface->dftrCommunicator());
    }
}

void MPIEwald::deliverData(const std::vector<Position>& cd,
				  const std::vector<double>& charge)
{
    if (ewaldInterface->isWaveDevided()) {
	assert(cd.size() == charge.size());
	assert(sizeof(Position) == sizeof(double)*Position::Dim);
	if (ewaldInterface->isWaveParent()) {
	    cdbuf.resize(numbuf+cd.size());
	    std::copy(cd.begin(), cd.end(), &cdbuf[numbuf]);
	    chargebuf.resize(numbuf+charge.size());
	    std::copy(charge.begin(), charge.end(), &chargebuf[numbuf]);
	    numbuf += cd.size();
	}
    }
}

void MPIEwald::deliverDataComplete()
{
    if (ewaldInterface->isWaveDevided()) {
	//std::cout << "DDC " << numbuf << std::endl;
	MPI_Bcast(&numbuf, 1, MPI_INTEGER,
                  ewaldInterface->waveParentRank(),
                  ewaldInterface->waveCommunicator());
	//std::cout << "DDCB " << numbuf << std::endl;
	cdbuf.resize(numbuf);
	MPI_Bcast(&cdbuf[0], numbuf*Position::Dim, MPI_DOUBLE,
                  ewaldInterface->waveParentRank(),
                  ewaldInterface->waveCommunicator());
	chargebuf.resize(numbuf);
	MPI_Bcast(&chargebuf[0], numbuf, MPI_DOUBLE,
                  ewaldInterface->waveParentRank(),
                  ewaldInterface->waveCommunicator());
    }
}

void MPIEwald::calcDFT(EwaldInterface::Context pContext,
                              EwaldMethod* p)
{
    if (ewaldInterface->isWaveDevided()) {
	assert(sizeof(Force)==sizeof(double)*Force::Dim);
	if (!ewaldInterface->isWaveParent()) {
	    p->calculateDFT(pContext, cdbuf, chargebuf);
	}
	else {
	    recvbuf.resize(numbuf);
	}
    }
}

void MPIEwald::calcIDFT(EwaldInterface::Context pContext,
                               EwaldMethod* p)
{
    if (ewaldInterface->isWaveDevided()) {
	assert(sizeof(Force)==sizeof(double)*Force::Dim);
	fcbuf.assign(numbuf, Force());
	if (!ewaldInterface->isWaveParent()) {
            p->calculateEwaldForce(pContext, cdbuf, chargebuf, fcbuf);
	}
	else {
	    recvbuf.resize(numbuf);
	}
	MPI_Reduce(&fcbuf[0], &recvbuf[0], numbuf*Force::Dim, MPI_DOUBLE,
		   MPI_SUM,
                   ewaldInterface->waveParentRank(),
                   ewaldInterface->waveCommunicator());
	numbuf = 0;
	cntbuf = 0;
	cdbuf.clear();
	chargebuf.clear();
    }
}

void MPIEwald::aggregateForce(std::vector<Position>& fc)
{
    if (ewaldInterface->isWaveDevided()) {
	assert(cntbuf+fc.size() <= recvbuf.size());
	for (std::vector<Position>::size_type i=0;i<fc.size();++i,++cntbuf) {
	    fc[i] += recvbuf[cntbuf];
	}
    }
}

void MPIEwald::aggregate(double& potentialEnergy)
{
    if (ewaldInterface->isWaveDevided()) {
	double sendval = potentialEnergy;
	MPI_Reduce(&sendval, &potentialEnergy, 1, MPI_DOUBLE, MPI_SUM,
                   ewaldInterface->waveParentRank(),
                   ewaldInterface->waveCommunicator());
    }
}

#endif
