#include <cstring>
#include "CalcForce.h"

#include "Log.h"
#include "Timer.h"

#ifdef CPPMD_ENABLE_LONGRANGE
# ifdef SAMPLE_LONG
#  ifdef FPCOLL
#   include <fj_tool/fjsamp.h>
#   define SAMPLE_LONG_START() fpcoll_start()
#   define SAMPLE_LONG_STOP() fpcoll_stop()
#  else
#   include <fj_tool/fipp.h>
#   define SAMPLE_LONG_START() fipp_start()
#   define SAMPLE_LONG_STOP() fipp_stop()
#  endif
# else
#  define SAMPLE_LONG_START()
#  define SAMPLE_LONG_STOP()
# endif
#endif

#ifdef FAPP_CALCFORCE
# include <fj_tool/fapp.h>
#endif
#ifdef K_PA
#include "fjcoll.h"
#endif

//! Calculate force and potential energy
// Input : Position, Charge, Atomtype, ID
// Output : Force, total Energy

//template<class CBPARTICLE>
CalcForce::CalcForce(std::vector<int>& sset,
                     const MPI_Comm long_comm)
  : unit_identifier(0),
    short_id(0),
    shortset(sset), 
    self_shorttarget_set(),
    longrangeparam(),
    //    longrange(unit_identifier,SpaceVector<int>(1,1,1)),
    longrange(unit_identifier,longrangeparam),
#ifdef CPPMD_ENABLE_LONGRANGE
 #ifdef CPPMD_ENABLE_FMM
  #ifdef CPPMD_ENABLE_MR3EXAFMM
    mr3fmm(unit_identifier,longrangeparam,long_comm),
  #else
    fmmlongrange(unit_identifier,longrangeparam,long_comm),
  #endif
 #elif defined(CPPMD_ENABLE_PMMM)
    pmmmlocal(),
    pmmmlongrange(unit_identifier,longrangeparam,long_comm),
 #else
  #ifdef CPPMD_ENABLE_EWALD
    ewaldinterface(new EwaldModule::EwaldInterface(unit_identifier,longrangeparam,long_comm)),
    ewaldlongrange(unit_identifier,longrangeparam,ewaldinterface),
  #else
   #ifdef STMEWALD
    stmelongrange(unit_identifier,longrangeparam,long_comm),
   #else
    #ifdef CPPMD_ENABLE_OLDPME
    pmeinterface(new PMEModule::PMEInterface(unit_identifier,longrangeparam,long_comm)),
    pmelongrange(unit_identifier,longrangeparam,pmeinterface),
    #else
    pmemoduleinterface(new PMEModuleInterface(unit_identifier,longrangeparam,long_comm)),
    pmelongrange(unit_identifier,longrangeparam,pmemoduleinterface),
    #endif
   #endif
  #endif
 #endif
#endif  // CPPMD_ENABLE_LONGRANGE
    cltype(ShortRange::NoCoulomb),
    cutoff(0.0)
{
  //    bondlist.clear();
  clear_shorttarget_set();
  cutoff2 = cutoff*cutoff;
  ShortRange::set_Shift_Function(cutoff2);
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
#ifdef PRE_CALC_SF
  set_shift_coefficient_(&cutoff);
#endif
#endif

  register_timers();
#ifdef CPPMD_ENABLE_PMMM
  pmmmlocal.register_timers();
  pmmmlongrange.register_timers();
#endif
}

CalcForce::CalcForce(std::vector<int> sset, const ShortRange::CoulombType _cltype,
                     const LongRangeParameter lrp,
                     const int unitid, const int shortid, const double cf,
#ifdef USE_PAIRLIST
                     const double pm,
                     const int num_list,
#endif
                     const MPI_Comm long_comm)
  : unit_identifier(unitid), short_id(shortid),
    shortset(sset), self_shorttarget_set(),  
    bondlist(), 
    longrangeparam(lrp),
    //    longrange(unitid,SpaceVector<int>(1,1,1)), 
    longrange(unitid,longrangeparam), 
#ifdef CPPMD_ENABLE_LONGRANGE
# ifdef CPPMD_ENABLE_FMM
#  ifdef CPPMD_ENABLE_MR3EXAFMM
    mr3fmm(unit_identifier,longrangeparam,long_comm),
#  else
    fmmlongrange(unit_identifier,longrangeparam,long_comm),
#  endif
# elif defined(CPPMD_ENABLE_PMMM)
  pmmmlocal(),
    pmmmlongrange(unit_identifier,longrangeparam,long_comm),
# else
#  ifdef CPPMD_ENABLE_EWALD
    ewaldinterface(new EwaldModule::EwaldInterface(unit_identifier,longrangeparam,long_comm)),
    ewaldlongrange(unit_identifier,longrangeparam,ewaldinterface),
#  else
   #ifdef STMEWALD
    stmelongrange(unit_identifier,longrangeparam,long_comm),
   #else
    #ifdef CPPMD_ENABLE_OLDPME
    pmeinterface(new PMEModule::PMEInterface(unit_identifier,longrangeparam,long_comm)),
    pmelongrange(unit_identifier,longrangeparam,pmeinterface),
    #else
    pmemoduleinterface(new PMEModuleInterface(unit_identifier,longrangeparam,long_comm)),
    pmelongrange(unit_identifier,longrangeparam,pmemoduleinterface),
    #endif
#   endif
#  endif
# endif
#endif  // CPPMD_ENABLE_LONGRANGE
    cltype(_cltype), cutoff(cf)
{
  //    bondlist.clear();
  if(DebugLog::verbose>1){
    std::cout << "num set " << sset.size() << std::endl; 
  }
  clear_shorttarget_set();
  cutoff2 = cutoff*cutoff;
  ShortRange::set_Shift_Function(cutoff2);
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
#ifdef PRE_CALC_SF
  set_shift_coefficient_(&cutoff);
#endif
#endif
#ifdef USE_PAIRLIST
  if(DebugLog::verbose>1){
    std::cout << "new pairlist " << std::endl; 
  }
  pmargin2 = (cutoff+pm)*(cutoff+pm)-cutoff2;
  allocate_list(num_list);
  if(DebugLog::verbose>1){
    std::cout << "newed pairlist " << std::endl; 
  }
#endif

  register_timers();
#ifdef CPPMD_ENABLE_PMMM
  pmmmlocal.register_timers();
  pmmmlongrange.register_timers();
#endif
}

CalcForce::CalcForce(std::vector<int> sset,
                     std::vector< std::vector<int> > targetset, 
                     std::vector< std::pair<int,int> > ghostpairs, 
                     const ShortRange::CoulombType _cltype, 
                     const LongRangeParameter lrp,
                     const int unitid, const int shortid, const double cf,
#ifdef USE_PAIRLIST
                     const double pm,
                     const int num_list,
#endif
                     const MPI_Comm long_comm,
                     bool expire_reverse)
  : unit_identifier(unitid), short_id(shortid),
    shortset(sset),   
    bondlist(), 
    longrangeparam(lrp),
    //    longrange(unitid,SpaceVector<int>(1,1,1)),
    longrange(unitid,longrangeparam),
#ifdef CPPMD_ENABLE_LONGRANGE
# ifdef CPPMD_ENABLE_FMM
#  ifdef CPPMD_ENABLE_MR3EXAFMM
    mr3fmm(unit_identifier,longrangeparam,long_comm),
#  else
    fmmlongrange(unit_identifier,longrangeparam,long_comm),
#  endif
# elif defined(CPPMD_ENABLE_PMMM)
		 pmmmlocal(),
    pmmmlongrange(unit_identifier,longrangeparam,long_comm),
# else
#  ifdef CPPMD_ENABLE_EWALD
    ewaldinterface(new EwaldModule::EwaldInterface(unit_identifier,longrangeparam,long_comm)),
    ewaldlongrange(unit_identifier,longrangeparam,ewaldinterface),
#  else
   #ifdef STMEWALD
    stmelongrange(unit_identifier,longrangeparam,long_comm),
   #else
    #ifdef CPPMD_ENABLE_OLDPME
    pmeinterface(new PMEModule::PMEInterface(unit_identifier,longrangeparam,long_comm)),
    pmelongrange(unit_identifier,longrangeparam,pmeinterface),
    #else
    pmemoduleinterface(new PMEModuleInterface(unit_identifier,longrangeparam,long_comm)),
    pmelongrange(unit_identifier,longrangeparam,pmemoduleinterface),
    #endif
   #endif
#  endif
# endif
#endif  // CPPMD_ENABLE_LONGRANGE
    cltype(_cltype), ghost_pair_set(ghostpairs), cutoff(cf)
{
  //  printf("CalcForce(...) %d\n",unit_identifier);
  //    bondlist.clear();
  // TODO make self_shorttarget_set and ghost_shorttarget_set from targetset
  clear_shorttarget_set();
  generate_shorttarget_set(targetset,expire_reverse);
  cutoff2 = cutoff*cutoff;
  ShortRange::set_Shift_Function(cutoff2);
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
#ifdef PRE_CALC_SF
  set_shift_coefficient_(&cutoff);
#endif
#endif
#ifdef USE_PAIRLIST
  if(DebugLog::verbose>1){
    std::cout << "new pairlist " << std::endl; 
  }
  pmargin2 = (cutoff+pm)*(cutoff+pm)-cutoff2;
  allocate_list(num_list);
  if(DebugLog::verbose>1){
    std::cout << "newed pairlist " << std::endl; 
  }
#endif

  register_timers();
#ifdef CPPMD_ENABLE_PMMM
  pmmmlocal.register_timers();
  pmmmlongrange.register_timers();
#endif
}

void 
CalcForce::register_timers()
{
#ifdef TIMER_DETAIL
  timer_short=PerfCounter::add_target(std::string("CalcShort"));
  timer_long=PerfCounter::add_target(std::string("CalcLong"));
  timer_cb=PerfCounter::add_target(std::string("CalcBond"));
  timer_cellpair=PerfCounter::add_target(std::string("CellPair"));
# ifdef CPPMD_ENABLE_PMMM
  timer_pmmm_pm_top_half=PerfCounter::add_target(std::string("PMMM_PM_TOP"));
  timer_pmmm_pm_comm=PerfCounter::add_target(std::string("PMMM_PM_COMM"));
  timer_pmmm_pm_bottom_half=PerfCounter::add_target(std::string("PMMM_PM_BOTTOM"));
  timer_pmmm_mm_prep=PerfCounter::add_target(std::string("PMMM_MM_PREP"));
  timer_pmmm_mm_calc=PerfCounter::add_target(std::string("PMMM_MM_CALC"));
  timer_pmmm_mm_post=PerfCounter::add_target(std::string("PMMM_MM_POST"));
# endif
#endif
}

void 
CalcForce::set_cutoff(double cf)
{
  cutoff = cf;
  cutoff2 = cutoff*cutoff;
  ShortRange::set_Shift_Function(cutoff2);
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
#ifdef PRE_CALC_SF
  set_shift_coefficient_(&cutoff);
#endif
#endif
}
#ifdef USE_PAIRLIST
void 
CalcForce::set_pmargin(double pm)
{
  pmargin2 = (cutoff+pm)*(cutoff+pm)-cutoff2;
}

void
CalcForce::allocate_list(size_t num_particle)
{
#ifdef INDEX_PAIRLIST
#ifndef SHARE_LIST
  pid = new int[num_particle][MAX_PAIR];
  plj = new int[num_particle][MAX_PAIR];
  iid = new int[num_particle];
#endif
  npair = new int[num_particle];
  fi = new double[num_particle][3];

  selfpid = new int[num_particle][MAX_PAIR];
  selfplj = new int[num_particle][MAX_PAIR];
  selfnpair = new int[num_particle];
  selfiid = new int[num_particle];
  selffi = new double[num_particle][3];
  if((DebugLog::verbose>1)||(short_id==0)){
    printf("Use Index PairList\n");
  }
#else
#ifndef SHARE_LIST
  pid = new int[num_particle][MAX_PAIR];
  pos = new double[num_particle][MAX_PAIR][3];
  charge = new double[num_particle][MAX_PAIR];
  plj = new double[num_particle][MAX_PAIR][4];
  iid = new int[num_particle];
  posi = new double[num_particle][3];
  chargei = new double[num_particle];
#endif
  npair = new int[num_particle];
  fi = new double[num_particle][3];

  selfpid = new int[num_particle][MAX_PAIR];
  selfpos = new double[num_particle][MAX_PAIR][3];
  selfcharge = new double[num_particle][MAX_PAIR];
  selfplj = new double[num_particle][MAX_PAIR][4];
  selfnpair = new int[num_particle];

  selfiid = new int[num_particle];
  selfposi = new double[num_particle][3];
  selfchargei = new double[num_particle];
  selffi = new double[num_particle][3];
  if((DebugLog::verbose>1)||(short_id==0)){
    printf("Use Value PairList\n");
  }
#endif
}

#endif

void 
CalcForce::clear_shorttarget_set()
{
  if(DebugLog::verbose>1){
    std::cout << "num set " << shortset.size() << std::endl;
  }
  self_shorttarget_set.resize(shortset.size());
  self_shorttarget_index.resize(shortset.size());
  ghost_shorttarget_set.resize(shortset.size());
  ghost_shorttarget_index.resize(shortset.size());
  for(size_t s=0;s<shortset.size();s++){
    self_shorttarget_set[s].clear();
    self_shorttarget_index[s].clear();
    ghost_shorttarget_set[s].clear();
    ghost_shorttarget_index[s].clear();
  }
    
}

void 
CalcForce::generate_shorttarget_set(std::vector< std::vector<int> >& targetset,
                                    bool expire_reverse)
{
  for(size_t s=0;s<shortset.size();s++){
    for(size_t ti=0;ti<targetset[s].size();ti++){
      int tid = targetset[s][ti];
      size_t i;
      for(i=0;i<shortset.size();i++){if(shortset[i]==tid)break;}
      if(i<shortset.size()){
        if((expire_reverse==false)||(s<i)){
          self_shorttarget_set[s].push_back(tid);
          self_shorttarget_index[s].push_back(i);
        }
      }else{
        ghost_shorttarget_set[s].push_back(tid);
      }
    }
#ifdef USE_PAIRLIST
    //! ShortRangeInteractionSet always calculate forces of the particle pairs in one own cell. The self_shorttarget_index[s] build upper block dose not include the cell indexed "s". But PairList  calculate just only self_shorttarget_index[s], so append "s" to self_shorttarget_index[s] if absent.
    {
      size_t ti;
      for(ti=0;ti<self_shorttarget_index[s].size();ti++){
        if(size_t(self_shorttarget_index[s][ti])==s)break;
      }
      if(ti>=self_shorttarget_index[s].size()){
        self_shorttarget_index[s].push_back(s);
      }else{
      }
    }
#endif
  }
}

#ifdef CPPMD_ENABLE_FMM
void 
CalcForce::clear_fmmtarget_set()
{
  if(DebugLog::verbose>1){
    std::cout << "num set " << shortset.size() << std::endl;
  }
  self_fmmtarget_set.resize(shortset.size());
  self_fmmtarget_index.resize(shortset.size());
  ghost_fmmtarget_set.resize(shortset.size());
  ghost_fmmtarget_index.resize(shortset.size());
  for(size_t s=0;s<shortset.size();s++){
    self_fmmtarget_set[s].clear();
    self_fmmtarget_index[s].clear();
    ghost_fmmtarget_set[s].clear();
    ghost_fmmtarget_index[s].clear();
  }
    
}

void 
CalcForce::generate_fmmtarget_set(std::vector< std::vector<int> >& fmmtargetset,
				  bool expire_reverse)
{
  for(size_t s=0;s<shortset.size();s++){
    for(size_t ti=0;ti<fmmtargetset[s].size();ti++){
      int tid = fmmtargetset[s][ti];
      size_t i;
      for(i=0;i<shortset.size();i++){if(shortset[i]==tid)break;}
      if(i<shortset.size()){
        if((expire_reverse==false)||(s<i)){
          self_fmmtarget_set[s].push_back(tid);
          self_fmmtarget_index[s].push_back(i);
        }
      }else{
        ghost_fmmtarget_set[s].push_back(tid);
      }
    }
    {
      size_t ti;
      for(ti=0;ti<self_fmmtarget_index[s].size();ti++){
        if(size_t(self_fmmtarget_index[s][ti])==s)break;
      }
      if(ti>=self_fmmtarget_index[s].size()){
        self_fmmtarget_index[s].push_back(s);
      }else{
      }
    }
  }
}
#endif

#ifdef CPPMD_ENABLE_PMMM
void 
CalcForce::clear_pmmmtarget_set()
{
  if(DebugLog::verbose>1){
    std::cout << "num set " << shortset.size() << std::endl;
  }
  self_pmmmtarget_set.resize(shortset.size());
  self_pmmmtarget_index.resize(shortset.size());
  ghost_pmmmtarget_set.resize(shortset.size());
  ghost_pmmmtarget_index.resize(shortset.size());
  for(size_t s=0;s<shortset.size();s++){
    self_pmmmtarget_set[s].clear();
    self_pmmmtarget_index[s].clear();
    ghost_pmmmtarget_set[s].clear();
    ghost_pmmmtarget_index[s].clear();
  }
    
}

void 
CalcForce::generate_pmmmtarget_set(std::vector< std::vector<int> >& pmmmtargetset,
				  bool expire_reverse)
{
  for(size_t s=0;s<shortset.size();s++){
    for(size_t ti=0;ti<pmmmtargetset[s].size();ti++){
      int tid = pmmmtargetset[s][ti];
      size_t i;
      for(i=0;i<shortset.size();i++){if(shortset[i]==tid)break;}
      if(i<shortset.size()){
        if((expire_reverse==false)||(s<i)){
          self_pmmmtarget_set[s].push_back(tid);
          self_pmmmtarget_index[s].push_back(i);
        }
      }else{
        ghost_pmmmtarget_set[s].push_back(tid);
      }
    }
    {
      size_t ti;
      for(ti=0;ti<self_pmmmtarget_index[s].size();ti++){
        if(size_t(self_pmmmtarget_index[s][ti])==s)break;
      }
      if(ti>=self_pmmmtarget_index[s].size()){
        self_pmmmtarget_index[s].push_back(s);
      }else{
      }
    }
  }
}
#endif

void 
CalcForce::convert_setid_to_index(std::map<int,int>& id_to_index,
                                  const std::vector<int>& id,
                                  std::vector<int>& index)
{
  index.resize(id.size());
  for(size_t s=0;s<id.size();s++){
    index[s] = id_to_index[id[s]];
  }
}

void 
CalcForce::convert_setid_to_index(
                  std::map<int,int>& id_to_index,
                  const std::vector< std::pair<int,int> >& pairid,
                  std::vector< std::pair<int,int> >& pairindex)
{
  pairindex.resize(pairid.size());
  for(size_t s=0;s<pairid.size();s++){
    pairindex[s] = std::make_pair(id_to_index[pairid[s].first],
                                  id_to_index[pairid[s].second]);
  }
}

void
CalcForce::set_ghost_longset_index(const std::map<int,int>& id_to_index)
{
  int n=0;
  for(std::vector<int>::size_type i=0;i<ghostlongset.size();i++){
    int id = ghostlongset[i];
    std::map<int,int>::const_iterator it = id_to_index.find(id);
    if(it!=id_to_index.end()){
      ghost_longset_index.push_back(it->second);
      n++;
    }else{
      std::cout << "Not found long target ID " << id << std::endl;
    }
  }
  for(std::vector<int>::size_type i=0;i<ghost_selfenergycell_list.size();i++){
    int id = ghost_selfenergycell_list[i];
    std::map<int,int>::const_iterator it = id_to_index.find(id);
    if(it!=id_to_index.end()){
      ghost_selfenergycell_index.push_back(it->second);
      n++;
    }else{
      std::cout << "Not found selfenergy target ID " << id << std::endl;
    }
  }
}

void
CalcForce::initialize(size_t pa_size, SpaceVector<int> celldiv3d, SpaceVector<int> nodediv3d, OperationSelector operations)
{
#ifdef CPPMD_ENABLE_LONGRANGE
  if(operations.doLongrangecalculation){
# ifdef CPPMD_ENABLE_FMM
#  ifdef CPPMD_ENABLE_MR3EXAFMM
    mr3fmm.initialize(pa_size, celldiv3d, nodediv3d);
#  else
#  endif
# elif defined(CPPMD_ENABLE_PMMM)
    //    pmmmlongrange.initialize();   // initilize in CalcPreparator
# else
#  ifdef CPPMD_ENABLE_EWALD
#  else
#   ifdef STMEWALD
#   else
    //    std::cout << "pmelongrange.initialize()" << std::endl;
    pmelongrange.initialize();
#   endif
#  endif
# endif
  }
#endif  // CPPMD_ENABLE_LONGRANGE
  covalentbond.initialize();
}

void
CalcForce::selfshort_zero()
{
  for(size_t i=0;i<shortforce.size();i++){
    shortforce[i].x = 0.0;
    shortforce[i].y = 0.0;
    shortforce[i].z = 0.0;
  }
  shortenergyself=0.0;
  shortenergy=0.0;
  shortenergyexcluded=0.0;
  shortvirial = 0.0;
}

void
CalcForce::ghostshort_zero()
{
  for(size_t i=0;i<ghostshortforce.size();i++){
    ghostshortforce[i].x = 0.0;
    ghostshortforce[i].y = 0.0;
    ghostshortforce[i].z = 0.0;
  }
}

void
CalcForce::short_zero()
{
  for(size_t i=0;i<shortforce.size();i++){
    shortforce[i].x = 0.0;
    shortforce[i].y = 0.0;
    shortforce[i].z = 0.0;
  }
  for(size_t i=0;i<ghostshortforce.size();i++){
    ghostshortforce[i].x = 0.0;
    ghostshortforce[i].y = 0.0;
    ghostshortforce[i].z = 0.0;
  }
  shortenergyself=0.0;
  shortenergy=0.0;
  shortenergyexcluded=0.0;
  shortvirial = 0.0;
}
void
CalcForce::long_zero()
{
  longenergy=0.0;
  longvirial=0.0;
}

template<typename GPA>
void dump_force_sum(const GPA& ghost)
{
  std::cout << " ghost.force[0] " << ghost.force[0] << std::endl;
  std::cout << " ghost.force[1] " << ghost.force[1] << std::endl;
  Force sum(0.0,0.0,0.0);
  for(int i=0;i<ghost.size();i++){
    sum += ghost.force[i];
  }
  std::cout << " sum of ghost pre long force " << sum << std::endl;
}

template<>
void dump_force_sum(const ParticleArray& ghost)
{
  std::cout << " ghost[0].force " << ghost[0].force << std::endl;
  std::cout << " ghost[1].force " << ghost[1].force << std::endl;
  Force sum(0.0,0.0,0.0);
  for(unsigned int i=0;i<ghost.size();i++){
    sum += ghost[i].force;
  }
  std::cout << " sum of ghost pre long force " << sum << std::endl;
}

template<typename PA>
void dump_force_sum(const PA& particlearray, 
                    const std::vector<TypeRange>& typerangearray)
{
  std::cout << " particlearray.force[0] " << particlearray.force[0] << std::endl;
  std::cout << " particlearray.force[1] " << particlearray.force[1] << std::endl;
  Force sum(0.0,0.0,0.0);
  for(int ci=0;ci<typerangearray.size();ci++){
    for(int i=typerangearray[ci].begin;i<typerangearray[ci].end;i++){
      sum +=  particlearray.force[i];
    }
  }
  std::cout << " sum of shortforce " << sum << std::endl;
}
template<>
void dump_force_sum(const ParticleArray& particlearray, 
                    const std::vector<TypeRange>& typerangearray)
{
  std::cout << " particlearray[0].force " << particlearray[0].force << std::endl;
  std::cout << " particlearray[1].force " << particlearray[1].force << std::endl;
  Force sum(0.0,0.0,0.0);
  for(unsigned int ci=0;ci<typerangearray.size();ci++){
    for(int i=typerangearray[ci].begin;i<typerangearray[ci].end;i++){
      sum +=  particlearray[i].force;
    }
  }
  std::cout << " sum of shortforce " << sum << std::endl;
}
template<>
void dump_force_sum(const ForceArray& forcearray, 
                    const std::vector<TypeRange>& typerangearray)
{
  std::cout << "force[0] " << forcearray[0] << std::endl;
  std::cout << "force[1] " << forcearray[1] << std::endl;
  Force sum(0.0,0.0,0.0);
  for(unsigned int ci=0;ci<typerangearray.size();ci++){
    for(int i=typerangearray[ci].begin;i<typerangearray[ci].end;i++){
      sum +=  forcearray[i];
    }
  }
  std::cout << " sum of shortforce " << sum << std::endl;
}

template<typename PA>
void add_force(PA& particlearray, 
               const std::vector<TypeRange>& typerangearray,
               const ForceArray& force)
{
  //  printf("%f %f\n",particlearray.force[0].x,force[0].x);
  for(unsigned int ci=0;ci<typerangearray.size();ci++){
    for(int i=typerangearray[ci].begin;i<typerangearray[ci].end;i++){
      particlearray.force[i] += force[i];
    }
  }
}
template<>
void add_force(ParticleArray& particlearray, 
               const std::vector<TypeRange>& typerangearray,
               const ForceArray& force)
{
  for(unsigned int ci=0;ci<typerangearray.size();ci++){
    for(int i=typerangearray[ci].begin;i<typerangearray[ci].end;i++){
      particlearray[i].force += force[i];
    }
  }
}

template<typename GPA>
void add_copy_force(ForceArray& ghostforce,
                    GPA& ghost, 
                    const ForceArray& ghostshortforce)
{
  for(size_t fi=0;fi<ghost.size();fi++){
    ghost.force[fi] += ghostshortforce[fi];
    ghostforce[fi] = ghost.force[fi];
  }
}
template<>
void add_copy_force(ForceArray& ghostforce,
                    ParticleArray& ghost, 
                    const ForceArray& ghostshortforce)

{
  for(size_t fi=0;fi<ghost.size();fi++){
    ghost[fi].force += ghostshortforce[fi];
    ghostforce[fi] = ghost[fi].force;
  }
}

template<class PA>
double
CalcForce::calcCorrectionExcludedWater_full_fixed(const PA& particlearray,
				       const WaterList& waterlist,
				       const double cutoff2,
				       const ShortRange::CoulombType clt)
{
#if 1
  WaterList::const_iterator it=waterlist.begin();
  int o_index = it->first;
  int h1_index = it->second.h1;
  int h2_index = it->second.h2;
  double co = getcharge(particlearray,o_index);
  double ch = getcharge(particlearray,h1_index);

  double ce;
  if(clt==ShortRange::ForEwaldReal){
    ce = energy_correction_fixed_water_full<ShortRange::ForEwaldReal>(co,ch,cutoff2);
  }else if(clt==ShortRange::ZeroDipole){
    ce = energy_correction_fixed_water_full<ShortRange::ZeroDipole>(co,ch,cutoff2);
  }
  return ce*waterlist.size();
#else
  double ce = 0.0;
  for(WaterList::const_iterator it=waterlist.begin(); it!=waterlist.end(); it++){
    int o_index = it->first;
    int h1_index = it->second.h1;
    int h2_index = it->second.h2;
    double co = getcharge(particlearray,o_index);
    double ch = getcharge(particlearray,h1_index);

    if(clt==ShortRange::ForEwaldReal){
      ce += energy_correction_fixed_water_full<ShortRange::ForEwaldReal>(co,ch,cutoff2);
    }else if(clt==ShortRange::ZeroDipole){
      ce += energy_correction_fixed_water_full<ShortRange::ZeroDipole>(co,ch,cutoff2);
    }
  }
  return ce;
#endif
}

template
double
CalcForce::calcCorrectionExcludedWater_full_fixed<CombinedParticleArray>(const CombinedParticleArray& particlearray,
				       const WaterList& waterlist,
				       const double cutoff2,
				       const ShortRange::CoulombType clt);

template<class PA>
double
calcCorrectionExcludedWater_full(PA& particlearray,
			    const WaterList& waterlist,
			    double cutoff2,
			    ShortRange::CoulombType clt)
{
  double ce=0.0;
  if(clt==ShortRange::ForEwaldReal){
#ifdef CALCCRECTIONFORCANCELEDFORCE
    ce = forceenergy_correction_excluded_water_full<CombinedParticleArray,ShortRange::ForEwaldReal>(particlearray,waterlist,cutoff2);
#else
    ce = energy_correction_excluded_water_full<CombinedParticleArray,ShortRange::ForEwaldReal>(particlearray,waterlist,cutoff2);
#endif
  }else if(clt==ShortRange::ZeroDipole){
#ifdef CALCCRECTIONFORCANCELEDFORCE
    ce = forceenergy_correction_excluded_water_full<CombinedParticleArray,ShortRange::ZeroDipole>(particlearray,waterlist,cutoff2);
#else
    ce = energy_correction_excluded_water_full<CombinedParticleArray,ShortRange::ZeroDipole>(particlearray,waterlist,cutoff2);
#endif
  }
#ifdef DUMP_ENERGY_CONTENT
  printf("Call energy_correction_excluded_water_full %e\n",ce);
#endif
  return ce;
}

template<class PA>
double
calcCorrectionExcludedWater(PA& particlearray,
			    const WaterList& waterlist,
			    double cutoff2,
			    ShortRange::CoulombType clt)
{
  double ce=0.0;
  if(clt==ShortRange::ForEwaldReal){
#ifdef CALCCRECTIONFORCANCELEDFORCE
    ce = forceenergy_correction_excluded_water<CombinedParticleArray,ShortRange::ForEwaldReal>(particlearray,waterlist,cutoff2);
#else
    ce = energy_correction_excluded_water<CombinedParticleArray,ShortRange::ForEwaldReal>(particlearray,waterlist,cutoff2);
#endif
  }else if(clt==ShortRange::ZeroDipole){
#ifdef CALCCRECTIONFORCANCELEDFORCE
    ce = forceenergy_correction_excluded_water<CombinedParticleArray,ShortRange::ZeroDipole>(particlearray,waterlist,cutoff2);
#else
    ce = energy_correction_excluded_water<CombinedParticleArray,ShortRange::ZeroDipole>(particlearray,waterlist,cutoff2);
#endif
  }
#ifdef DUMP_ENERGY_CONTENT
  printf("Call energy_correction_excluded_water %e\n",ce);
#endif
  return ce;
}

#ifdef CPPMD_ENABLE_PMMM
template<class PA>
void
CalcForce::calcPMMM_top_half(const PA& particlearray, 
			     const std::vector<TypeRange>& typerangearray,
			     OperationSelector operation)
{
#ifdef TIMER_DETAIL
  PerfCounter::start(timer_pmmm_pm_top_half);
#endif
  if(operation.doShortrangecalculation){
    pmmmlocal.calc_local_top_half(particlearray, typerangearray);
  }
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
#ifdef TIMER_DETAIL
  PerfCounter::start(timer_pmmm_mm_prep);
#endif
  if(operation.doLongrangecalculation){
    pmmmlongrange.prep();
  }
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
}

void
CalcForce::make_waterexcludelist(const size_t num,
				 const WaterList& waterlist,
				 std::vector< waterexclude > &wex)
{
  //  int n=0;
  wex.clear();
  wex.resize(num);
  for(WaterList::const_iterator wli = waterlist.begin(); wli != waterlist.end(); wli++ ){
    int wo = wli->first;
    int wh1 = wli->second.h1;
    int wh2 = wli->second.h2;
    wex[wo].exid[0] = wh1;
    wex[wo].exid[1] = wh2;
    wex[wh1].exid[0] = wo;
    wex[wh1].exid[1] = wh2;
    wex[wh2].exid[0] = wo;
    wex[wh2].exid[1] = wh1;
    //    n++;
    //    printf("Found %d water, exclude %d %d %d\n",n,wo,wh1,wh2);
  }

  //  printf("Found %d water for exclude\n",n);
}
#endif

/*!
  @note covalentbond.map_particle must be called at outside of this method.
 */
template<class PA, class GPA>
void 
CalcForce::calcForce_local(PA& particlearray, 
                     const std::vector<TypeRange>& typerangearray, 
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
		     const WaterList& waterlist,
                     GPA& ghost, 
                     const std::vector<int>& ghostsetid, 
                     std::map<int,int>& recvsetid_to_index, 
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<CovalentBondInfo::BondList>& ghostbond,
                     ForceArray& force, 
                     ForceArray& ghostforce, 
                     double& energy,
		     double& virial,
                     OperationSelector operations)
{

#ifdef TIMER_DETAIL
  PerfCounter::start(timer_short);
#endif
  if(operations.doShortrangecalculation) {
#ifdef CPPMD_ENABLE_FMM
# ifdef CPPMD_ENABLE_MR3EXAFMM
#ifdef TIMER_DETAIL
    PerfCounter::start(timer_cellpair);
#endif
#ifdef K_PA
    start_collection("localcellpair");
#endif
    /// DEBUG CODE
    /*
    if(unit_identifier==0){
      for(int ci=0;ci<self_shorttarget_index.size();ci++){
	for(int i=0;i<self_shorttarget_index[ci].size();i++){
	  printf("%d ",self_shorttarget_index[ci][i]);
	}
	printf("\n");
      }
    }
    if(unit_identifier==0){
      for(int ci=0;ci<self_fmmtarget_index.size();ci++){
	for(int i=0;i<self_fmmtarget_index[ci].size();i++){
	  printf("%d ",self_fmmtarget_index[ci][i]);
	}
	printf("\n");
      }
    }
    */
    {
      cellpairsinteraction(particlearray,typerangearray,
			   particlearray,typerangearray, self_fmmtarget_index,
			   shortforce, shortenergy,true, operations);
    }
#ifdef K_PA
    stop_collection("localcellpair");
#endif
#ifdef TIMER_DETAIL
    PerfCounter::stop();
#endif
# endif   // CPPMD_ENABLE_MR3EXAFMM
#endif   // CPPMD_ENABLE_FMM

#ifdef CPPMD_ENABLE_PMMM
#ifdef TIMER_DETAIL
    PerfCounter::start(timer_cellpair);
#endif
#ifdef K_PA
    start_collection("localcellpair");
#endif
    {
# ifdef EXCLUDE_WATERFORCE
      if(operations.doMakePairList){
	make_waterexcludelist(particlearray.size(), waterlist, waterexcludelist);
	//	printf("watrexcludelist.size() %d\n",waterexcludelist.size());
      }
      cellpairsinteraction(particlearray,typerangearray,waterexcludelist,
			   particlearray,typerangearray, self_pmmmtarget_index,
			   shortforce, shortenergy,true, operations);
# else
      cellpairsinteraction(particlearray,typerangearray,
			   particlearray,typerangearray, self_pmmmtarget_index,
			   shortforce, shortenergy,true, operations);
# endif
    }
#ifdef K_PA
    stop_collection("localcellpair");
#endif
#ifdef TIMER_DETAIL
    PerfCounter::stop();
#endif
#endif   // CPPMD_ENABLE_PMMM

#ifdef USE_PAIRLIST
    if(operations.doMakePairList){
      bool mkpl;
#ifdef K_PA
      start_collection("makepairlist");
#endif
# ifdef INDEX_PAIRLIST
      {
#  if !defined(CPPMD_ENABLE_FMM) && !defined(CPPMD_ENABLE_PMMM)
      // make INDEX PairList
       mkpl = makeljcpairlist(selfpid, selfplj, selfnpair, 
                              selfiid, &selfnpl,
                              particlearray,
                              typerangearray,
#ifndef INCLUDE_WATER
                              waterlist,
#endif
                              particlearray,
                              typerangearray, 
                              self_shorttarget_index,
                              cutoff2, pmargin2, true);
#  else   // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
       mkpl = makeljpairlist(selfpid, selfplj, selfnpair, 
                             selfiid, &selfnpl,
                             particlearray,
                             typerangearray,
                             particlearray,
                             typerangearray, 
                             self_shorttarget_index,
                             cutoff2, pmargin2, true);
#  endif  // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
      }
# else // not INDEX PAIRLIST
      {
       mkpl = makeljcpairlist(selfpid, selfpos, selfcharge, selfplj, selfnpair, 
                              selfiid, selfposi, selfchargei, &selfnpl,
                              particlearray, 
                              typerangearray,
                              particlearray, 
                              typerangearray, 
                              self_shorttarget_index,
                              cutoff2, pmargin2, true);
      }
# endif // not INDEX PAIRLIST
      if(!mkpl){
        printf("makeljcpairlist fail\n");
      }
#ifdef K_PA
      stop_collection("makepairlist");
#endif
    }else{
# ifdef INDEX_PAIRLIST
  // do nothing
# else // !INDEX_PAIRLIST
#ifdef K_PA
      start_collection("updatepairlist");
 #endif
      // update Position in PairList
      {
       updatepairlistpos1(selfiid,selfposi,selfnpl,particlearray);
#  ifndef SHARE_LIST
       updatepairlistpos1(iid,posi,npl,particlearray);
 #  endif
       updatepairlistpos(selfpid,selfpos,selfnpair,selfnpl,particlearray);
      }
#ifdef K_PA
      stop_collection("updatepairlist");
#endif
# endif // !INDEX_PAIRLIST
    }
#ifdef FAPP_CALCFORCE
    fapp_start("calclocalplforce",0,1);
#endif
#ifdef K_PA
    start_collection("calclocalplforce");
#endif
    {
      memset((void*)(fi),0,selfnpl*3*sizeof(double));
    }
#ifdef OVERLAP
    MPICHECKER::mpi_checker();
#endif
# ifdef INDEX_PAIRLIST
    {
#  if !defined( CPPMD_ENABLE_FMM) && !defined( CPPMD_ENABLE_PMMM)
      if(cltype==ShortRange::ForEwaldReal){
#ifdef SIMPLE_LJ
       pairlistloopljcfe_ewald(selfpid, selfplj, selfnpair, 
                               selfiid, 
                               particlearray,
                               particlearray,
                               selfnpl,
                               cutoff2, 
                               fi, shortenergy, shortvirial, operations);
#else
       pairlistloopljcfe_ewald_ljshift(selfpid, selfplj, selfnpair, 
				       selfiid, 
				       particlearray,
				       particlearray,
				       selfnpl,
				       cutoff2, 
				       fi, shortenergy, shortvirial, operations);
#endif
      }else if(cltype==ShortRange::MWolf) {
       pairlistloopljcfe_wolf_ljshift(selfpid, selfplj, selfnpair, 
                                      selfiid, 
                                      particlearray,
                                      particlearray,
                                      selfnpl,
                                      cutoff2, 
                                      fi, shortenergy, shortvirial, operations);
      }else if(cltype==ShortRange::ZeroDipole) {
#ifdef SIMPLE_LJ
       pairlistloopljcfe_zerodipole(selfpid, selfplj, selfnpair, 
				    selfiid, 
				    particlearray,
				    particlearray,
				    selfnpl,
				    cutoff2, 
				    fi, shortenergy, shortvirial, operations);
#else
       pairlistloopljcfe_zerodipole_ljshift(selfpid, selfplj, selfnpair, 
					    selfiid, 
					    particlearray,
					    particlearray,
					    selfnpl,
					    cutoff2, 
					    fi, shortenergy, shortvirial, operations);
#endif
       shortenergyself += calc_self_energy_zerodipole(selfiid,particlearray,selfnpl,cutoff2);
      }else {
       pairlistloopljcfe(selfpid, selfplj, selfnpair, 
                         selfiid, 
                         particlearray,
                         particlearray,
                         selfnpl,
                         cutoff2, 
                         fi, shortenergy, shortvirial, operations);
      }
#  else // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
      {
       pairlistloopljfe(selfpid, selfplj, selfnpair, 
                        selfiid, 
                        particlearray,
                        particlearray,
                        selfnpl,
                        cutoff2, 
                        fi, shortenergy, shortvirial, operations);
      }
#  endif // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
    }
# else // ! INDEX_PAIRLIST
    {
      pairlistloopljcfe(selfpos, selfcharge, selfplj, selfnpair, 
                       selfposi, selfchargei, selfnpl,
                       cutoff2, 
                       fi, shortenergy, shortvirial, operations);
    }
# endif
#ifdef K_PA
    stop_collection("calclocalplforce");
#endif
#ifdef FAPP_CALCFORCE
    fapp_stop("calclocalplforce",0,1);
#endif
  }
#endif // USE_PAIRLIST

#ifndef INCLUDE_WARTER
  if(operations.doShortrangecalculation) {
# ifdef FULL_CORRECTTION_EXWATER
#  ifndef CORRECTTION_EXWATER_DYNAMIC
    shortenergyexcluded += calcCorrectionExcludedWater_full_fixed(particlearray,waterlist,cutoff2,cltype);
#  else
    shortenergyexcluded += calcCorrectionExcludedWater_full(particlearray,waterlist,cutoff2,cltype);
#  endif
# else
    shortenergyexcluded += calcCorrectionExcludedWater(particlearray,waterlist,cutoff2,cltype);
# endif
  }
#endif
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
}

template<class PA, class GPA>
void 
CalcForce::calcLongForce(PA& particlearray, 
			 const std::vector<TypeRange>& typerangearray, 
			 GPA& ghost, 
			 const std::vector<TypeRange>& ghosttyperange)
{
#ifdef CPPMD_ENABLE_LONGRANGE
  SAMPLE_LONG_START();
# ifdef CPPMD_ENABLE_FMM
#  ifdef CPPMD_ENABLE_MR3EXAFMM
    mr3fmm.calccoulomb(particlearray, typerangearray, self_longset_index,
		       self_selfenergycell_index,
		       ghost, ghosttyperange, ghost_longset_index,
		       ghost_selfenergycell_index,
		       longenergy);
#  else
    fmmlongrange.calcForce(particlearray, typerangearray, self_longset_index,
			   self_selfenergycell_index,
			   ghost, ghosttyperange, ghost_longset_index,
			   ghost_selfenergycell_index,
			   longenergy);
#  endif
# elif defined(CPPMD_ENABLE_PMMM)
#ifdef TIMER_DETAIL
    PerfCounter::start(timer_pmmm_mm_calc);
#endif
    pmmmlongrange.calcForce();
#ifdef TIMER_DETAIL
    PerfCounter::stop();
#endif
# else
#  ifdef CPPMD_ENABLE_EWALD
    ewaldlongrange.calcForce(particlearray, typerangearray, self_longset_index,
                             self_selfenergycell_index,
                             ghost, ghosttyperange, ghost_longset_index,
                             ghost_selfenergycell_index,
                             longenergy, longvirial);
#  else
#   ifdef STMEWALD
    stmelongrange.calcForce(particlearray, typerangearray, self_longset_index,
			    self_selfenergycell_index,
			    ghost, ghosttyperange, ghost_longset_index,
			    ghost_selfenergycell_index,
			    longenergy, longvirial);
#   else
    pmelongrange.calcForce(particlearray, typerangearray, self_longset_index,
                           self_selfenergycell_index,
                           ghost, ghosttyperange, ghost_longset_index,
                           ghost_selfenergycell_index,
                           longenergy, longvirial);
#   endif
#  endif
# endif
# ifdef MT_LONG
    {
      double mtfactor = (double)MT_LONG;
      for(int ti=0;ti<self_longset_index.size();ti++){
	int t = self_longset_index[ti];
	for(int i=typerangearray[t].begin;i<typerangearray[t].end;i++){
	  getforce(particlearray,i)*=mtfactor;
	}
      }
      for(int ti=0;ti<ghost_longset_index.size();ti++){
	int t = ghost_longset_index[ti];
	for(int i=ghosttyperange[t].begin;i<ghosttyperange[t].end;i++){
	  getforce(ghost,i)*=mtfactor;
	}
      }
    }
# endif
    SAMPLE_LONG_STOP();
#endif
}


/*!
  @note covalentbond.map_particle must be called at outside of this method.
 */
template<class PA, class GPA>
void 
CalcForce::calcForce(PA& particlearray, 
                     const std::vector<TypeRange>& typerangearray, 
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
		     const WaterList& waterlist,
                     GPA& ghost, 
                     const std::vector<int>& ghostsetid, 
                     std::map<int,int>& recvsetid_to_index, 
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<CovalentBondInfo::BondList>& ghostbond,
                     ForceArray& force, 
                     ForceArray& ghostforce, 
                     double& energy,
		     double& virial,
                     OperationSelector operations,
		     bool with_local)
{
  double bondenergy=0.0;
  //    bondlist.clear();
  /*
    TODO
    force and ghostforce are not used.
    every interaction add force to particlearray[*].force and ghost[*].force
    It cause shor/long/bond can't operate in thread prallel, like OpemMP
    make separate force buffer and each interaction return force not in ParticleArray
    For multiple timestep, each force must be hold separately at uncalculated step.
  */

#ifdef TIMER_DETAIL
  PerfCounter::start(timer_short);
#endif
  if(operations.doShortrangecalculation) {
    for(size_t s=0;s<ghost_shorttarget_set.size();s++){
      convert_setid_to_index(recvsetid_to_index,ghost_shorttarget_set[s],
                             ghost_shorttarget_index[s]);
    }
    convert_setid_to_index(recvsetid_to_index, ghost_pair_set,
                           ghost_pair_index);   
# ifdef CPPMD_ENABLE_FMM
#  ifdef CPPMD_ENABLE_MR3EXAFMM
#ifdef TIMER_DETAIL
  PerfCounter::start(timer_cellpair);
#endif
#ifdef K_PA
    start_collection("cellpair");
#endif
    for(size_t s=0;s<ghost_fmmtarget_set.size();s++){
      convert_setid_to_index(recvsetid_to_index,ghost_fmmtarget_set[s],
                             ghost_fmmtarget_index[s]);
    }
    
    {
      if(with_local){
	cellpairsinteraction(particlearray,typerangearray,
			     particlearray,typerangearray, self_fmmtarget_index,
			     shortforce, shortenergy,true, operations);
      }
      cellpairsinteraction(particlearray,typerangearray,
			   ghost, ghosttyperange, ghost_fmmtarget_index,
			   shortforce, shortenergy,false, operations);
    }
#ifdef K_PA
    stop_collection("cellpair");
#endif
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
#  endif
# endif

# ifdef CPPMD_ENABLE_PMMM
#ifdef TIMER_DETAIL
  PerfCounter::start(timer_cellpair);
#endif
#ifdef K_PA
    start_collection("cellpair");
#endif
    for(size_t s=0;s<ghost_pmmmtarget_set.size();s++){
      convert_setid_to_index(recvsetid_to_index,ghost_pmmmtarget_set[s],
                             ghost_pmmmtarget_index[s]);
    }
    
    {
      if(with_local){
# ifdef EXCLUDE_WATERFORCE
	if(operations.doMakePairList){
	  make_waterexcludelist(particlearray.size(), waterlist, waterexcludelist);
	}
	cellpairsinteraction(particlearray,typerangearray,waterexcludelist,
			     particlearray,typerangearray, self_pmmmtarget_index,
			     shortforce, shortenergy,true, operations);
# else
	cellpairsinteraction(particlearray,typerangearray,
			     particlearray,typerangearray, self_pmmmtarget_index,
			     shortforce, shortenergy,true, operations);
# endif
      }
      cellpairsinteraction(particlearray,typerangearray,
			   ghost, ghosttyperange, ghost_pmmmtarget_index,
			   shortforce, shortenergy,false, operations);
    }
#ifdef K_PA
    stop_collection("cellpair");
#endif
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
# endif
#ifdef CPPMD_ENABLE_PMMM
    if(DebugLog::verbose>1)printf("cellpairsinteraction Done\n");
#endif

#ifdef USE_MR3_SHORT
    
# ifndef CPPMD_ENABLE_FMM
    mr3.calcforce(particlearray,typerangearray,
		  self_shorttarget_index,
		  ghost,ghosttyperange, ghost_shorttarget_index,
		  ghost_pair_index,
		  shortforce, ghostshortforce,
		  shortenergyself,shortenergy,cutoff2);
# else
    mr3.calcvdw(particlearray,typerangearray,
		self_shorttarget_index,
		ghost,ghosttyperange, ghost_shorttarget_index,
		ghost_pair_index,
		shortforce, ghostshortforce,
		shortenergyself,shortenergy,cutoff2);
# endif
#else

#ifdef USE_PAIRLIST
    if(operations.doMakePairList){
      bool mkpl=true;
#ifdef K_PA
      start_collection("makepairlist");
#endif
# ifdef INDEX_PAIRLIST
      if(with_local){
#  if !defined( CPPMD_ENABLE_FMM) && !defined (CPPMD_ENABLE_PMMM)
	// make INDEX PairList
	mkpl = makeljcpairlist(selfpid, selfplj, selfnpair, 
			       selfiid, &selfnpl,
			       particlearray,
			       typerangearray,
#ifndef INCLUDE_WATER
			       waterlist,
#endif
			       particlearray,
			       typerangearray, 
			       self_shorttarget_index,
			       cutoff2, pmargin2, true);
#  else
	mkpl = makeljpairlist(selfpid, selfplj, selfnpair, 
			      selfiid, &selfnpl,
			      particlearray,
			      typerangearray,
			      particlearray,
			      typerangearray, 
                             self_shorttarget_index,
                             cutoff2, pmargin2, true);
#  endif
      }
#  ifdef SHARE_LIST
#   if !defined( CPPMD_ENABLE_FMM) && !defined( CPPMD_ENABLE_PMMM)
      mkpl = mkpl&&makeljcpairlist(selfpid, selfplj, npair, 
                                   selfiid, &selfnpl,
                                   particlearray, 
                                   typerangearray,
                                   ghost,
                                   ghosttyperange, ghost_shorttarget_index,
                                   selfnpair,
                                   cutoff2, pmargin2, false);
#   else
      mkpl = mkpl&&makeljpairlist(selfpid, selfplj, npair, 
                                   selfiid, &selfnpl,
                                   particlearray, 
                                   typerangearray,
                                   ghost,
                                   ghosttyperange, ghost_shorttarget_index,
                                   selfnpair,
                                   cutoff2, pmargin2, false);
#   endif
#  else  // Not SHARE_LIST
#   if !defined( CPPMD_ENABLE_FMM ) !defined( CPPMD_ENABLE_PMMM) 
      mkpl = mkpl&&makeljcpairlist(pid, plj, npair, 
                                   iid, &npl,
                                   particlearray, 
                                   typerangearray,
                                   ghost, 
                                   ghosttyperange, ghost_shorttarget_index,
                                   cutoff2, pmargin2, false);
#   else
      mkpl = mkpl&&makeljpairlist(pid, plj, npair, 
                                   iid, &npl,
                                   particlearray, 
                                   typerangearray,
                                   ghost, 
                                   ghosttyperange, ghost_shorttarget_index,
                                   cutoff2, pmargin2, false);
#   endif
#  endif
# else // not INDEX PAIRLIST
      if(with_local){
	mkpl = makeljcpairlist(selfpid, selfpos, selfcharge, selfplj, selfnpair, 
			       selfiid, selfposi, selfchargei, &selfnpl,
			       particlearray, 
			       typerangearray,
			       particlearray, 
			       typerangearray, 
			       self_shorttarget_index,
			       cutoff2, pmargin2, true);
      }
#ifdef SHARE_LIST
      mkpl = mkpl&&makeljcpairlist(selfpid, selfpos, selfcharge, selfplj, npair, 
                                   selfiid, selfposi, selfchargei, &selfnpl,
                                   particlearray, 
                                   typerangearray,
                                   ghost, 
                                   ghosttyperange, ghost_shorttarget_index,
                                   selfnpair,
                                   cutoff2, pmargin2, false);
#else
      mkpl = mkpl&&makeljcpairlist(pid, pos, charge, plj, npair, 
                                   iid, posi, chargei, &npl,
                                   particlearray,
                                   typerangearray,
                                   ghost,
                                   ghosttyperange, ghost_shorttarget_index,
                                   cutoff2, pmargin2, false);
#endif
# endif
      if(!mkpl){
        printf("makeljcpairlist fail\n");
      }
#ifdef K_PA
      stop_collection("makepairlist");
#endif
    }else{
# ifdef INDEX_PAIRLIST
  // do nothing
# else
#ifdef K_PA
      start_collection("updatepairlist");
#endif
      // update Position in PairList
      if(with_local){
	updatepairlistpos1(selfiid,selfposi,selfnpl,particlearray);
#  ifndef SHARE_LIST
	updatepairlistpos1(iid,posi,npl,particlearray);
#  endif
	updatepairlistpos(selfpid,selfpos,selfnpair,selfnpl,particlearray);
      }
#  ifdef SHARE_LIST
      updatepairlistpos(selfpid,selfpos,npair,npl,ghost,selfnpair);
#  else
      updatepairlistpos(pid,pos,npair,npl,ghost);
#  endif
#ifdef K_PA
      stop_collection("updatepairlist");
#endif
# endif
    }
#ifdef FAPP_CALCFORCE
    fapp_start("calcpairlistforce",1,1);
#endif
#ifdef K_PA
    start_collection("calcpairlistforce");
#endif
    if(with_local){
      memset((void*)(fi),0,selfnpl*3*sizeof(double));
    }
# ifdef INDEX_PAIRLIST
    if(with_local){
#  if !defined (CPPMD_ENABLE_FMM ) && !defined (CPPMD_ENABLE_PMMM)
      if(cltype==ShortRange::ForEwaldReal){
#ifdef SIMPLE_LJ
	pairlistloopljcfe_ewald(selfpid, selfplj, selfnpair, 
				selfiid, 
				particlearray,
				particlearray,
				selfnpl,
				cutoff2, 
				fi, shortenergy, shortvirial, operations);
#else
	pairlistloopljcfe_ewald_ljshift(selfpid, selfplj, selfnpair, 
					selfiid, 
					particlearray,
					particlearray,
					selfnpl,
					cutoff2, 
					fi, shortenergy, shortvirial, operations);
#endif
      }else if(cltype==ShortRange::MWolf) {
	pairlistloopljcfe_wolf_ljshift(selfpid, selfplj, selfnpair, 
				       selfiid, 
				       particlearray,
				       particlearray,
				       selfnpl,
				       cutoff2, 
				       fi, shortenergy, shortvirial, operations);
      }else if(cltype==ShortRange::ZeroDipole) {
#ifdef SIMPLE_LJ
	pairlistloopljcfe_zerodipole(selfpid, selfplj, selfnpair, 
				     selfiid, 
				     particlearray,
				     particlearray,
				     selfnpl,
				     cutoff2, 
				     fi, shortenergy, shortvirial, operations);
#else
	pairlistloopljcfe_zerodipole_ljshift(selfpid, selfplj, selfnpair, 
					     selfiid, 
					     particlearray,
					     particlearray,
					     selfnpl,
					     cutoff2, 
					     fi, shortenergy, shortvirial, operations);
#endif
	shortenergyself += calc_self_energy_zerodipole(selfiid,particlearray,selfnpl,cutoff2);
      }else {
	pairlistloopljcfe(selfpid, selfplj, selfnpair, 
			  selfiid, 
			  particlearray,
			  particlearray,
			  selfnpl,
			  cutoff2, 
			  fi, shortenergy, shortvirial, operations);
      }
#  else // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
      {
	pairlistloopljfe(selfpid, selfplj, selfnpair, 
			 selfiid, 
			 particlearray,
			 particlearray,
			 selfnpl,
			 cutoff2, 
			 fi, shortenergy, shortvirial, operations);
      }
#  endif // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
    }
#  ifdef SHARE_LIST
#   if !defined( CPPMD_ENABLE_FMM ) && !defined( CPPMD_ENABLE_PMMM )
    if(cltype==ShortRange::ForEwaldReal){
#ifdef SIMPLE_LJ
      pairlistloopljcfe_ewald(selfpid, selfplj, npair, 
                              selfiid, 
                              particlearray,
                              ghost,
                              selfnpl,
                              selfnpair,
                              cutoff2, 
                              fi, shortenergy, shortvirial, operations);
#else
      pairlistloopljcfe_ewald_ljshift(selfpid, selfplj, npair, 
				      selfiid, 
				      particlearray,
				      ghost,
				      selfnpl,
				      selfnpair,
				      cutoff2, 
				      fi, shortenergy, shortvirial, operations);
#endif
    }else if(cltype==ShortRange::MWolf) {
      pairlistloopljcfe_wolf_ljshift(selfpid, selfplj, npair, 
				     selfiid, 
				     particlearray,
				     ghost,
				     selfnpl,
				     selfnpair,
				     cutoff2, 
				     fi, shortenergy, shortvirial, operations);
    }else if(cltype==ShortRange::ZeroDipole) {
#ifdef SIMPLE_LJ
      pairlistloopljcfe_zerodipole(selfpid, selfplj, npair, 
				   selfiid, 
				   particlearray,
				   ghost,
				   selfnpl,
				   selfnpair,
				   cutoff2, 
				   fi, shortenergy, shortvirial, operations);
#else
      pairlistloopljcfe_zerodipole_ljshift(selfpid, selfplj, npair, 
					   selfiid, 
					   particlearray,
					   ghost,
					   selfnpl,
					   selfnpair,
					   cutoff2, 
					   fi, shortenergy, shortvirial, operations);
#endif
    }else {
      pairlistloopljcfe(selfpid, selfplj, npair, 
                      selfiid, 
                      particlearray,
                      ghost,
                      selfnpl,
                      selfnpair,
                      cutoff2, 
                      fi, shortenergy, shortvirial, operations);
    }
#   else // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
    {
      pairlistloopljfe(selfpid, selfplj, npair, 
                       selfiid, 
                       particlearray,
                       ghost,
                       selfnpl,
                       selfnpair,
                       cutoff2, 
                       fi, shortenergy, shortvirial, operations);
    }
#   endif // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
#  else // !SHARE_LIST
#   if !defined( CPPMD_ENABLE_FMM) !defined( CPPMD_ENABLE_PMMM)
    if(cltype==ShortRange::ForEwaldReal){
#ifdef SIMPLE_LJ
      pairlistloopljcfe_ewald(pid, plj, npair, 
                              iid, 
                              particlearray,
                              ghost,
                              npl,
                              cutoff2, 
                              fi, shortenergy, shortvirial, operations);
#else
      pairlistloopljcfe_ewald_ljshift(pid, plj, npair, 
				      iid, 
				      particlearray,
				      ghost,
				      npl,
				      cutoff2, 
				      fi, shortenergy, shortvirial, operations);
#endif
    }else if(cltype==ShortRange::MWolf) {
      pairlistloopljcfe_wolf_ljshift(pid, plj, npair, 
				     iid, 
				     particlearray,
				     ghost,
				     npl,
				     cutoff2, 
				     fi, shortenergy, shortvirial, operations);
    }else if(cltype==ShortRange::ZeroDipole) {
#ifdef SIMPLE_LJ
      pairlistloopljcfe_zerodipole(pid, plj, npair, 
				   iid, 
				   particlearray,
				   ghost,
				   npl,
				   cutoff2, 
				   fi, shortenergy, shortvirial, operations);
#else
      pairlistloopljcfe_zerodipole_ljshift(pid, plj, npair, 
					   iid, 
					   particlearray,
					   ghost,
					   npl,
					   cutoff2, 
					   fi, shortenergy, shortvirial, operations);
#endif
    }else{
      pairlistloopljcfe(pid, plj, npair, 
                        iid, 
                        particlearray,
                        ghost,
                        npl,
                        cutoff2, 
                        fi, shortenergy, shortvirial, operations);
    }
#   else // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
    pairlistloopljfe(pid, plj, npair, 
                      iid, 
                      particlearray,
                      ghost,
                      npl,
                      cutoff2, 
                      fi, shortenergy, shortvirial, operations);
#   endif // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
#  endif // !SHARE_LIST
# else // ! INDEX_PAIRLIST
    if(with_local){
      pairlistloopljcfe(selfpos, selfcharge, selfplj, selfnpair, 
			selfposi, selfchargei, selfnpl,
			cutoff2, 
			fi, shortenergy, shortvirial, operations);
    }
#  ifdef SHARE_LIST
    pairlistloopljcfe(selfpos, selfcharge, selfplj, npair, 
                      selfposi, selfchargei, selfnpl,
                      selfnpair,
                      cutoff2, 
                      fi, shortenergy, shortvirial, operations);
#  else
    pairlistloopljcfe(pos, charge, plj, npair, 
                      posi, chargei, npl,
                      cutoff2, 
                      fi, shortenergy, shortvirial, operations);
#  endif
# endif
# ifdef SHARE_LIST
    importpairforce(shortforce,selfiid,fi,selfnpl);
# else
    importpairforce(shortforce,iid,fi,npl);
# endif

#ifndef INCLUDE_WARTER
    if(with_local){
      if(operations.doShortrangecalculation) {
# ifdef FULL_CORRECTTION_EXWATER
#  ifndef CORRECTTION_EXWATER_DYNAMIC
	shortenergyexcluded += calcCorrectionExcludedWater_full_fixed(particlearray,waterlist,cutoff2,cltype);
#  else
	shortenergyexcluded += calcCorrectionExcludedWater_full(particlearray,waterlist,cutoff2,cltype);
#  endif
# else
	shortenergyexcluded += calcCorrectionExcludedWater(particlearray,waterlist,cutoff2,cltype);
# endif
      }
    }
#endif

#ifdef K_PA
    stop_collection("calcpairlistforce");
#endif
#ifdef FAPP_CALCFORCE
    fapp_stop("calcpairlistforce",1,1);
#endif
#else  // !USE_PAIRLIST
#ifdef OLDPARTICLE  /// TODO not implemet shortrange for new particle
# if !defined( CPPMD_ENABLE_FMM) !defined( CPPMD_ENABLE_PMMM)
    if(cltype==ShortRange::ForEwaldReal){
      if(cutoff2>0.0){
        ewald_and_lj.allloop(particlearray,typerangearray,
                             self_shorttarget_index,
                             ghost,ghosttyperange, ghost_shorttarget_index,
                             ghost_pair_index,
                             shortforce, ghostshortforce,
                             shortenergyself,shortenergy,cutoff2);
      }else{
        ewald_and_lj.allloop(particlearray,typerangearray,
                             self_shorttarget_index,
                             ghost,ghosttyperange, ghost_shorttarget_index,
                             ghost_pair_index,
                             shortforce, ghostshortforce,
                             shortenergyself,shortenergy);
      }
    }else if(cltype==ShortRange::OriginalCoulomb){
      if(cutoff2>0.0){
        cl_and_lj.allloop(particlearray,typerangearray,
                          self_shorttarget_index,
                          ghost,ghosttyperange,ghost_shorttarget_index,
                          ghost_pair_index,
                          shortforce, ghostshortforce,
                          shortenergyself,shortenergy,cutoff2);
      }else{
        cl_and_lj.allloop(particlearray,typerangearray,
                          self_shorttarget_index,
                          ghost,ghosttyperange,ghost_shorttarget_index,
                          ghost_pair_index,
                          shortforce, ghostshortforce,
                          shortenergyself,shortenergy);
      }
    }else{
      if(cutoff2>0.0){
        lj.allloop(particlearray,typerangearray,self_shorttarget_index,
                   ghost,ghosttyperange,ghost_shorttarget_index,
                   ghost_pair_index,
                   shortforce, ghostshortforce,
                   shortenergyself,shortenergy,cutoff2);
      }else{
        lj.allloop(particlearray,typerangearray,self_shorttarget_index,
                   ghost,ghosttyperange,ghost_shorttarget_index,
                   ghost_pair_index,
                   shortforce, ghostshortforce,
                   shortenergyself,shortenergy);
      }
    }
# else  // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
    {
      if(cutoff2>0.0){
        lj.allloop(particlearray,typerangearray,self_shorttarget_index,
                   ghost,ghosttyperange,ghost_shorttarget_index,
                   ghost_pair_index,
                   shortforce, ghostshortforce,
                   shortenergyself,shortenergy,cutoff2);
      }else{
        lj.allloop(particlearray,typerangearray,self_shorttarget_index,
                   ghost,ghosttyperange,ghost_shorttarget_index,
                   ghost_pair_index,
                   shortforce, ghostshortforce,
                   shortenergyself,shortenergy);
      }
    }
# endif // CPPMD_ENABLE_FMM or CPPMD_ENABLE_PMMM
    // Full Cube, ghostforce of short must be 0
    //! TODO: make allloop without ghostforce  
    if(operations.doReturnShortForce==false){
      ghostshort_zero();
    }
#endif // OLDPARTICLE
#endif // USE_PAIRLIST

#endif // ! USE_MR3
#ifdef CPPMD_ENABLE_PMMM
    if(DebugLog::verbose>1)printf("doShortrangecalculation Done\n");
#endif

  }   /// end operations.doShortrangecalculation
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif

#ifdef CPPMD_ENABLE_LONGRANGE
  SAMPLE_LONG_START();
#ifdef TIMER_DETAIL
  PerfCounter::start(timer_long);
#endif
#ifdef K_PA
  start_collection("calclongrangeforce");
#endif
  
#ifdef CPPMD_ENABLE_PMMM
  if(operations.doShortrangecalculation) {
    if(DebugLog::verbose>1)printf("calc_long_progress\n");
#ifdef TIMER_DETAIL
    PerfCounter::start(timer_pmmm_pm_comm);
#endif
    pmmmlocal.calc_long_progress();
#ifdef TIMER_DETAIL
    PerfCounter::stop();
#endif
    if(DebugLog::verbose>1)printf("calc_long_progress Done, calc_local_bottom_half\n");
#ifdef TIMER_DETAIL
    PerfCounter::start(timer_pmmm_pm_bottom_half);
#endif
    pmmmlocal.calc_local_bottom_half(particlearray, typerangearray, particlearray.force, longenergy);
#ifdef TIMER_DETAIL
    PerfCounter::stop();
#endif
    if(DebugLog::verbose>1)printf("calc_local_bottom_half Done\n");
    longenergy *= 0.5;
  }
#endif
  if(operations.doLongrangecalculation) {
# ifdef CPPMD_ENABLE_FMM
#  ifdef CPPMD_ENABLE_MR3EXAFMM
    mr3fmm.calccoulomb(particlearray, typerangearray, self_longset_index,
		       self_selfenergycell_index,
		       ghost, ghosttyperange, ghost_longset_index,
		       ghost_selfenergycell_index,
		       longenergy);
    /// DEBUG CODE
    /*
    if(unit_identifier==0){
      for(int i=0;i<self_longset_index.size();i++){
	printf("%d ",self_longset_index[i]);
      }
      printf("\n");
      for(int i=0;i<ghost_longset_index.size();i++){
	printf("%d ",ghost_longset_index[i]);
      }
      printf("\n");
    }
    */
#  else
    fmmlongrange.calcForce(particlearray, typerangearray, self_longset_index,
                             self_selfenergycell_index,
                             ghost, ghosttyperange, ghost_longset_index,
                             ghost_selfenergycell_index,
                             longenergy);
#  endif
# elif defined(CPPMD_ENABLE_PMMM)
    if(DebugLog::verbose>1)printf("pmmmlongrange.calcForce\n");
#ifdef TIMER_DETAIL
    PerfCounter::start(timer_pmmm_mm_calc);
#endif
    pmmmlongrange.calcForce();
#ifdef TIMER_DETAIL
    PerfCounter::stop();
#endif
    if(DebugLog::verbose>1)printf("pmmmlongrange.calcForce Done, pmmmlongrange.post\n");
#ifdef TIMER_DETAIL
    PerfCounter::start(timer_pmmm_mm_post);
#endif
    pmmmlongrange.post();
#ifdef TIMER_DETAIL
    PerfCounter::stop();
#endif
    if(DebugLog::verbose>1)printf("pmmmlongrange.post Done\n");
# else
#  ifdef CPPMD_ENABLE_EWALD
    ewaldlongrange.calcForce(particlearray, typerangearray, self_longset_index,
                             self_selfenergycell_index,
                             ghost, ghosttyperange, ghost_longset_index,
                             ghost_selfenergycell_index,
                             longenergy, longvirial);
#  else
#   ifdef STMEWALD
    stmelongrange.calcForce(particlearray, typerangearray, self_longset_index,
                             self_selfenergycell_index,
                             ghost, ghosttyperange, ghost_longset_index,
                             ghost_selfenergycell_index,
			    longenergy, longvirial);
#   else
    pmelongrange.calcForce(particlearray, typerangearray, self_longset_index,
                           self_selfenergycell_index,
                           ghost, ghosttyperange, ghost_longset_index,
                           ghost_selfenergycell_index,
                           longenergy, longvirial);
#   endif
#  endif
# endif
#ifdef MT_LONG
    {
      double mtfactor = (double)MT_LONG;
      for(int ti=0;ti<self_longset_index.size();ti++){
	int t = self_longset_index[ti];
	for(int i=typerangearray[t].begin;i<typerangearray[t].end;i++){
	  getforce(particlearray,i)*=mtfactor;
	}
      }
      for(int ti=0;ti<ghost_longset_index.size();ti++){
	int t = ghost_longset_index[ti];
	for(int i=ghosttyperange[t].begin;i<ghosttyperange[t].end;i++){
	  getforce(ghost,i)*=mtfactor;
	}
      }
    }
#endif
  }
#ifdef K_PA
  stop_collection("calclongrangeforce");
#endif
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
  SAMPLE_LONG_STOP();
#endif  // CPPMD_ENABLE_LONGRANGE

  if (operations.doShortEnergyToHalf) {
    shortenergy *= 0.5;
    shortvirial *= 0.5;
  }
#ifndef INCLUDE_WARTER
  shortenergy += shortenergyexcluded;
#endif
  double cbcse=0.0;
  double cbcsv=0.0;
  if(operations.doCovalentBondcaculation) {
#ifdef TIMER_DETAIL
  PerfCounter::start(timer_cb);
#endif
#ifdef K_PA
    start_collection("covalentbond");
#endif
#if 1
    covalentbond.calcForce(bondlist,bondenergy,cbcse,cbcsv);
    shortenergy += cbcse;
    shortvirial += cbcsv;
    //    if(short_id==0)printf("energy CBcancel %e\n",cbcse);
#else
    covalentbond.calcForce(bondlist,bondenergy,shortenergy,shortvirial);
#endif
#ifdef K_PA
    stop_collection("covalentbond");
#endif
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
  }

  shortenergy += shortenergyself;
  if(operations.doEnergycalculation){
    energy += shortenergy + longenergy + bondenergy;
#ifdef DUMP_ENERGY_CONTENT
    printf("energy (%e + %e + %e) + %e + %e = %e : cb14 %e\n", shortenergy-shortenergyself-cbcse, cbcse, shortenergyself, longenergy, bondenergy,energy, covalentbond.get_cb14energy());
#endif
  }
  if(operations.doVirialcalculation){
    virial += shortvirial + longvirial;  
    //    printf("virial in CalcForce (%e + %e) + %e = %e\n",shortvirial-cbcsv,cbcsv,longenergy,virial);
  }

  add_force(particlearray,typerangearray,shortforce);
  /*
    TODO
    interactions must be set force directory ghostforce
  */
  // for Full cell not require
  ghostforce.resize(ghost.size());
  add_copy_force(ghostforce, ghost, ghostshortforce);

}


#ifdef OLDPARTICLE
template 
void CalcForce::calcForce_local<ParticleArray,ParticleArray>(ParticleArray& particlearray, 
                     const std::vector<TypeRange>& typerangearray, 
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
                     const WaterList& waterlist,
                     ParticleArray& ghost, 
                     const std::vector<int>& ghostsetid, 
                     std::map<int,int>& recvsetid_to_index, 
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<CovalentBondInfo::BondList>& ghostbond,
                     ForceArray& force, 
                     ForceArray& ghostforce, 
		     double& virial,
                     double& energy,
							     OperationSelector operations);
template 
void CalcForce::calcForce<ParticleArray,ParticleArray>(ParticleArray& particlearray, 
                     const std::vector<TypeRange>& typerangearray, 
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
                     const WaterList& waterlist,
                     ParticleArray& ghost, 
                     const std::vector<int>& ghostsetid, 
                     std::map<int,int>& recvsetid_to_index, 
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<CovalentBondInfo::BondList>& ghostbond,
                     ForceArray& force, 
                     ForceArray& ghostforce, 
		     double& virial,
                     double& energy,
							     OperationSelector operations,
							     bool with_local);

#else
#ifdef CPPMD_ENABLE_PMMM
template
void
CalcForce::calcPMMM_top_half(const CombinedParticleArray& particlearray, 
			     const std::vector<TypeRange>& typerangearray,
			     OperationSelector operation);
#endif
template 
void CalcForce::calcForce_local<CombinedParticleArray,GhostParticleArray>(CombinedParticleArray& particlearray, 
                     const std::vector<TypeRange>& typerangearray, 
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
                     const WaterList& waterlist,
                     GhostParticleArray& ghost, 
                     const std::vector<int>& ghostsetid, 
                     std::map<int,int>& recvsetid_to_index, 
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<CovalentBondInfo::BondList>& ghostbond,
                     ForceArray& force, 
                     ForceArray& ghostforce, 
                     double& energy,
		     double& virial,
                     OperationSelector operations);
template 
void CalcForce::calcForce<CombinedParticleArray,GhostParticleArray>(CombinedParticleArray& particlearray, 
                     const std::vector<TypeRange>& typerangearray, 
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                     const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
                     const WaterList& waterlist,
                     GhostParticleArray& ghost, 
                     const std::vector<int>& ghostsetid, 
                     std::map<int,int>& recvsetid_to_index, 
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<CovalentBondInfo::BondList>& ghostbond,
                     ForceArray& force, 
                     ForceArray& ghostforce, 
                     double& energy,
		     double& virial,
								    OperationSelector operations,
								    bool with_local);
#endif
