// -*- mode: C++; -*-
#ifndef LONGRANGEINTERACTION_H
#define LONGRANGEINTERACTION_H

#include <ctime>
#include "Common.h"
#include "ParticleInfo.h"
#include "Geometry.h"
#include "LongRangeParameter.h"

//#define LONGRANGE_DT
#ifdef LONGRANGE_DT
#include <typeinfo>
#include <cxxabi.h>
#endif

#ifdef K_PA
#include "fjcoll.h"
#endif

//! grid data
struct GridData{
  double *gridvalue;
  size_t size;
  //  int inittype;

  GridData()
    : gridvalue()
  {
    //    inittype = 0;
    //    std::cout << "construct GridData " << gridvalue  << std::endl;
  }

  GridData(SpaceVector<int> num)
      : size(num.x*num.y*num.z)
  {
    //    inittype = 1;
    //    std::cout << "no init gridvalue " << gridvalue << std::endl;
    gridvalue = new double[num.x*num.y*num.z];
    //    std::cout << "construct GridData with num " << gridvalue << std::endl;
  }

  GridData(const GridData &gd)
      : gridvalue(),
        size(gd.size)
  {
    //    inittype = 2;
    //    std::cout << "no init gridvalue " << gridvalue << std::endl;
    gridvalue = new double[size];
    for(size_t i=0; i<size;i++){
      gridvalue[i] = gd.gridvalue[i];
    }
    //    std::cout << "construct GridData with gd " << gridvalue << std::endl;
  }

  GridData& operator=(const GridData& gd)
  {
    if(this != &gd){
      size = gd.size;
      //      std::cout << "inittype " << inittype << std::endl;
      //      std::cout << "destruct buffer in = operator " << gridvalue << std::endl;
      delete [] gridvalue;
      //      std::cout << "inittype " << inittype << std::endl;
      gridvalue = new double[size];
      for(size_t i=0; i<size;i++){
        gridvalue[i] = gd.gridvalue[i];
      }
    }
    return *this;
  }

  /* not work ?
   */
  ~GridData(){
    //    std::cout << "destruct GridData" << gridvalue << std::endl;
    delete [] gridvalue;
  }
};

//! charge assign and backinterpolate for Particle-Mesh
/*! particle charge assign to grid points
  and backinterpolate grid point force to particle
*/

class ChargeAssign {
 public:
  int unit_identifier;

  ChargeAssign(){
    // std::cout << "construct ChargeAssign" << std::endl;
  }

  ChargeAssign(int unitid) : unit_identifier(unitid) {}

  ChargeAssign(int unitid, const LongRangeParameter& _param) 
      : unit_identifier(unitid) {}

  virtual ~ChargeAssign() {}

#if 0
  template<class PA, class GPA>
  void assign(PA& particlearray,
              const std::vector<ParticleRange>& selfrange,
              GPA& ghost,
              const std::vector<ParticleRange>& ghostrange,
              GridData& gridcharge) {
  }
#endif
  /*!
    ghost and ghostrange may not be used in backinterpolate
    because this node should not calculate force for ghost
  */
  template<class PA, class GPA>
  void backinterpolate(PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential) {
  }
#if 0
  template<class PA, class GPA>
  void backinterpolate(const LongRangeParameter& param,
                       PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential) {
  }
#endif
  template<class PA, class GPA>
  void addenergy(PA& particlearray,
                 const std::vector<ParticleRange>& selfrange,
                 GPA& ghost,
                 const std::vector<ParticleRange>& ghostrange,
                 double& energy) {
  }

  enum ChargeAssignType{
    Dummy, RealSpace
  };
};


//! Possson Solver
class PoissonSolver {
 public:
  int unit_identifier;

  PoissonSolver(){}

  PoissonSolver(int unitid) : unit_identifier(unitid) {}

  PoissonSolver(int unitid, const LongRangeParameter& _param) : unit_identifier(unitid) {}

  virtual ~PoissonSolver() {}

  virtual void solvePoisson(GridData& gridcharge,
                            GridData& gridpotential, double& energy, double &virial){
  }
  /*
  virtual void solvePoisson(const LongRangeParameter& param,
                            GridData& gridcharge,
                            GridData& gridpotential, double& energy){
  }
  */

  enum PoissonSolverType{
    Dummy, RealSpace
  };
};



//! LongRangeInteraction, exs. PME, LGM, ...
template<class CHARGEASSIGN, class POISSONSOLVER, class MODULEINTERFACE=void>
class LongRangeInteraction {
 public:
  SpaceVector<int> number_of_grid_point;
  GridData gridcharge;
  GridData gridpotential;
  CHARGEASSIGN chargeassign;
  POISSONSOLVER poissonsolver;
  int unit_identifier;

  std::vector<ParticleRange> selfrange;
  std::vector<ParticleRange> ghostrange;

  std::vector<ParticleRange> self_selfenergy_range;
  std::vector<ParticleRange> ghost_selfenergy_range;
  /*
  LongRangeInteraction()
  {
    //    std::cout << "construct LongRangeInteraction no arg " << std::endl;
  }
  */

  LongRangeInteraction(int unitid, SpaceVector<int> grid_num )
      : number_of_grid_point(grid_num),
        gridcharge(number_of_grid_point),
        gridpotential(number_of_grid_point),
        chargeassign(unitid),
        poissonsolver(unitid),
        unit_identifier(unitid)
  {
    lr_time_clear();
    //    std::cout << "construct LongRangeInteraction with grid_num " << std::endl;
  }
  LongRangeInteraction(int unitid, const LongRangeParameter& param)
      : number_of_grid_point(param.grid_num),
        gridcharge(number_of_grid_point),
        gridpotential(number_of_grid_point),
        chargeassign(unitid, param),
        poissonsolver(unitid, param),
        unit_identifier(unitid)
  {
    lr_time_clear();
    //    std::cout << "construct LongRangeInteraction with param " << std::endl;
  }

  LongRangeInteraction(int unitid, const LongRangeParameter& param, MODULEINTERFACE* interface)
      : number_of_grid_point(param.grid_num),
        gridcharge(number_of_grid_point),
        gridpotential(number_of_grid_point),
        chargeassign(unitid, param, interface),
        poissonsolver(unitid, param, interface),
        unit_identifier(unitid)
  {
    lr_time_clear();
    //    std::cout << "construct LongRangeInteraction with param " << std::endl;
  }

  ~LongRangeInteraction() {
    //    delete &gridpotential;
    //    delete &gridcharge;
  }

  //! select charged particle range from typerage indexed by index and store to range
  /*
    poisson solver requires only charged particle (that in range ljcoulmb and coulomb), exclude lj only particles (that in range lj)
  */

  void initialize();

  inline void generate_particlerange(const std::vector<TypeRange>& typerange,
                                     const std::vector<int>& index,
                                     std::vector<ParticleRange>& range)
  {
    range.resize(index.size());
    for(size_t si=0;si<index.size(); ++si){
      //      std::cout << " typerange[index[" << si << "] " << index[si] << "] " << typerange[index[si]].end - typerange[index[si]].begin << std::endl;
      range[si].begin = typerange[index[si]].ljcoulomb.begin;
      range[si].end = typerange[index[si]].coulomb.end;
    }
  }
  template<class PA, class GPA>
  void calcForce(PA& particlearray,
                 const std::vector<TypeRange>& typerange,
                 const std::vector<int>& self_longset_index,
                 const std::vector<int>& self_selfenergycell_index,
                 GPA& ghost,
                 const std::vector<TypeRange>& ghosttyperange,
                 const std::vector<int>& ghost_longset_index,
                 const std::vector<int>& ghost_selfenergycell_index,
                 double& energy, double &virial)
  {
    lr_time_start();
#ifdef K_PA
    start_collection("longparticlerange");
#endif
    generate_particlerange(typerange,self_longset_index,selfrange);
    generate_particlerange(ghosttyperange,ghost_longset_index,ghostrange);
    generate_particlerange(typerange,self_selfenergycell_index,self_selfenergy_range);
    generate_particlerange(ghosttyperange,ghost_selfenergycell_index,ghost_selfenergy_range);
#ifdef K_PA
    stop_collection("longparticlerange");
#endif
    //    std::cout << "number of longrange target set (self, ghost) " << selfrange.size() << " " << ghostrange.size() << std::endl;
    lr_section_time(time_range);
#ifdef K_PA
    start_collection("longassign");
#endif
    chargeassign.assign(particlearray,selfrange,ghost,ghostrange,gridcharge,self_selfenergy_range,ghost_selfenergy_range);
#ifdef K_PA
    stop_collection("longassign");
#endif
    lr_section_time(time_assign);
#ifdef K_PA
    start_collection("longpoissonsolver");
#endif
    poissonsolver.solvePoisson(gridcharge, gridpotential, energy, virial);
#ifdef K_PA
    stop_collection("longpoissonsolver");
#endif
    lr_section_time(time_solve);
#ifdef K_PA
    start_collection("longbackinterpolate");
#endif
    chargeassign.backinterpolate(particlearray,selfrange,
                                 ghost,ghostrange,
                                 gridpotential);
#ifdef K_PA
    stop_collection("longbackinterpolate");
#endif
    lr_section_time(time_back);
#ifdef DUMP_LONG_ENERGY
    std::cout << "LongRangeInteraction without self energy " << energy << std::endl;
#endif
    chargeassign.addenergy(particlearray,selfrange,ghost,ghostrange, energy);
#ifdef DUMP_LONG_ENERGY
    std::cout << "LongRangeInteraction energy " << energy << std::endl;
#endif
    lr_section_time(time_energy);
    lr_time_end();
  }

  timespec time_all;
  timespec time_range;
  timespec time_assign;
  timespec time_solve;
  timespec time_back;
  timespec time_energy;
  timespec lr_fast;
  timespec lr_last;
  timespec lr_current;

  void lr_time_clear()
  {
    timespec zero;
    zero.tv_sec = 0;
    zero.tv_nsec = 0;
    time_all = zero;
    time_range = zero;
    time_assign = zero;
    time_solve = zero;
    time_back = zero;
    time_energy = zero;
    lr_fast = zero;
    lr_last = zero;
    lr_current = zero;
  }

#ifdef LONGRANGE_DT
  inline 
  void lr_time_print(const std::string& name, const timespec& time)
  {
    std::cout << name << " " << long(time.tv_sec)*1000000+time.tv_nsec/1000 << " micro sec" << std::endl;
    //  std::cout << name << " " << time.tv_sec << "." << time.tv_nsec << " sec.nsec" << std::endl;
  }
  inline
  void lr_time_dump()
  {
    std::cout << "Grid size " << gridcharge.size <<  " " << gridpotential.size << std::endl;
    lr_time_print("all",time_all);
    lr_time_print("range",time_range);
    lr_time_print("assign",time_assign);
    lr_time_print("solve",time_solve);
    lr_time_print("back",time_back);
    lr_time_print("energy",time_energy);
  }
  inline
  void lr_time_start()
  {
    clock_gettime(CLOCK_REALTIME,&lr_last);
    lr_fast = lr_last;
  }
  inline
  void lr_time_end()
  {
    clock_gettime(CLOCK_REALTIME,&lr_current);
    time_all.tv_sec += (lr_current.tv_sec - lr_fast.tv_sec);
    time_all.tv_nsec += (lr_current.tv_nsec - lr_fast.tv_nsec);
  }
  inline
  void lr_section_time(timespec &section_time)
  {
    clock_gettime(CLOCK_REALTIME,&lr_current);
    //    lr_time_print("current",lr_current);
    section_time.tv_sec += (lr_current.tv_sec - lr_last.tv_sec);
    section_time.tv_nsec += (lr_current.tv_nsec - lr_last.tv_nsec);
    lr_last = lr_current;
  }
#else
  inline void lr_time_print(const std::string& name, const timespec& time){}
  inline void lr_time_dump(){}
  inline void lr_time_start(){}
  inline void lr_time_end(){}
  inline void lr_section_time(timespec &section_time){}
#endif

};

#endif
