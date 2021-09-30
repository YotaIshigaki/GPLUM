#include "LongRangeInteraction.h"


class PMEModuleInterface {
public:
  PMEModuleInterface(int unitid, const LongRangeParameter& _param, MPI_Comm lcomm=MPI_COMM_WORLD) 
  : param(_param), myworld(lcomm) 
  {}

  LongRangeParameter param;
  MPI_Comm myworld;
};

class PMEChargeAssign : public ChargeAssign {
public:
  PMEChargeAssign(int unitid, const LongRangeParameter& _param,
		  PMEModuleInterface *pmemi)
    : pmemoduleinterface(*pmemi)
  {
  }

  template<class PA, class GPA>
  void assign(PA& particlearray,
              const std::vector<ParticleRange>& selfrange,
              GPA& ghost,
              const std::vector<ParticleRange>& ghostrange,
              GridData& gridcharge,
              const std::vector<ParticleRange>& self_selfenergy_range,
              const std::vector<ParticleRange>& ghost_selfenergy_range);
  template<class PA, class GPA>
  void backinterpolate(PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
  template<class PA, class GPA>
  void addenergy(PA& particlearray,
                 const std::vector<ParticleRange>& selfrange,
                 GPA& ghost,
                 const std::vector<ParticleRange>& ghostrange,
                 double& energy);

  PMEModuleInterface pmemoduleinterface;
};

class PMEPoissonSolver : public PoissonSolver {
public:
  PMEPoissonSolver(int unitid,
		   const LongRangeParameter& _param,
		   PMEModuleInterface* pmemi);
  void solvePoisson(GridData& gridcharge,
                    GridData& gridpotential, double& energy);
  void solvePoisson(GridData& gridcharge,
                    GridData& gridpotential, double& energy, double &virial);
};


typedef LongRangeInteraction<PMEChargeAssign, PMEPoissonSolver, PMEModuleInterface> PMELongRangeInteraction;
