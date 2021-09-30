#include "CalculationUnit.h"

template<class PA, class GPA>
CalculationUnit<PA,GPA>::CalculationUnit(MPICommunicator _communicator,
                                 MPIReceiverForLong _lrecv,
                                 MPISenderForLong _lsend,
                                 PostProcess _postprocess,
                                 ShortRange::CoulombType cltype,
                                 const LongRangeParameter longrangeparam,
                                 const int unitid, const int shortid,
                                 const SpaceVector<double> bs,
				 const double bms,
                                 const MPI_Comm long_comm,
                                 const Config::TempCtrl& temp_control) 
  : unit_identifier(unitid),
    short_id(shortid),
    boxsize(bs),
    integrator(unitid,shortid,typerangearray,bondlistarray,bondlistarray_idx,bs,bms,temp_control), 
    calcforce(setid, cltype, longrangeparam, unitid, shortid, double(0.0), 
#ifdef USE_PAIRLIST
              double(0.0),
              int(particlearray.size()),
#endif
              long_comm),
    postprocess(_postprocess),
    communicator(_communicator),
    longreceiver(_lrecv),
    longsender(_lsend)
{
  //    std::cout << "construct CalculationUnit()"  << std::endl;
}

template<class PA, class GPA>
CalculationUnit<PA,GPA>::CalculationUnit(PA& _particlearray, 
                                 std::vector<int> sid,
                                 std::vector<TypeRange> typerange,
                                 WaterList _waterlist,
                                 ShakeList _shakelist,
                                 std::vector<CovalentBondInfo::BondList> bondlist,
                                 std::vector<CovalentBondInfo::BondList> bondlist_idx,
                                 const CovalentBondParameterList& cbplist,
                                 std::vector< std::vector<int> > targetset,
                                 std::vector< std::pair<int,int> > ghostpairs,
                                 OperationSelector ope,
                                 MPICommunicator& _communicator,
                                 MPIReceiverForLong& _lrecv,
                                 MPISenderForLong& _lsend,
                                 PostProcess _postprocess,
                                 const ShortRange::CoulombType cltype,
                                 const LongRangeParameter longrangeparam,
                                 const int unitid, const int shortid,
                                 const SpaceVector<double> bs, 
				 const double bms,
                                 const double cutoff,
#ifdef USE_PAIRLIST
                                 const double pmargin,
#endif
                                 int num_t, int num_f,
                                 const MPI_Comm long_comm,
                                 const Config::TempCtrl& temp_control,
                                 bool expire_reverse
                                 ) 
  : unit_identifier(unitid),
    short_id(shortid),
    num_total(num_t),
    num_freedom(num_f),
    particlearray(_particlearray),
    setid(sid),
    typerangearray(typerange),
    waterlist(_waterlist),
    shakelist(_shakelist),
    bondlistarray(bondlist),
    bondlistarray_idx(bondlist_idx),
    boxsize(bs),
    integrator(unitid,shortid,typerangearray,bondlistarray,bondlistarray_idx,bs,bms,temp_control),
    calcforce(setid, targetset, ghostpairs, cltype, longrangeparam, unitid, shortid, cutoff, 
#ifdef USE_PAIRLIST
              pmargin,
              int(particlearray.size()),
#endif
	      long_comm,expire_reverse),
    operations(ope),
    postprocess(_postprocess),
    communicator(_communicator),
    longreceiver(_lrecv),
    longsender(_lsend)
{
  //  printf("CalculationUnit %d \n", unit_identifier);
  //  fflush(stdout);
  //  MPI_Barrier(MPI_COMM_WORLD);
  calcforce.covalentbond.setParameterList(cbplist);
  make_id_to_index(setid,setid_to_index);
  //    std::cout << "construct CalculationUnit("  << unit_identifier << ")" << std::endl;
  //  if(longsender.sprf.size()>0){
  //    std::cout << " longsender.sprf[0].send_requestp  " << _lsend.sprf[0].send_requestp << " " << longsender.sprf[0].send_requestp  << std::endl;
  //  }
}

#ifdef OLDPARTICLE
template
CalculationUnit<ParticleArray,ParticleArray>::CalculationUnit(ParticleArray& _particlearray, 
                                 std::vector<int> sid,
                                 std::vector<TypeRange> typerange,
                                 WaterList _waterlist,
                                 ShakeList _shakelist,
                                 std::vector<CovalentBondInfo::BondList> bondlist,
                                 std::vector<CovalentBondInfo::BondList> bondlist_idx,
                                 const CovalentBondParameterList& cbplist,
                                 std::vector< std::vector<int> > targetset,
                                 std::vector< std::pair<int,int> > ghostpairs,
                                 OperationSelector ope,
                                 MPICommunicator& _communicator,
                                 MPIReceiverForLong& _lrecv,
                                 MPISenderForLong& _lsend,
                                 PostProcess _postprocess,
                                 const ShortRange::CoulombType cltype,
                                 const LongRangeParameter longrangeparam,
                                 const int unitid, const int shortid,
                                 const SpaceVector<double> bs, 
				 const double bms,
                                 const double cutoff,
#ifdef USE_PAIRLIST
                                 const double pmargin,
#endif
                                 int num_t, int num_f,
                                 const MPI_Comm long_comm,
                                 const Config::TempCtrl& temp_control,
                                 bool expire_reverse
                                         );
#else
template
CalculationUnit<CombinedParticleArray,GhostParticleArray>::CalculationUnit(CombinedParticleArray& _particlearray, 
                                 std::vector<int> sid,
                                 std::vector<TypeRange> typerange,
                                 WaterList _waterlist,
                                 ShakeList _shakelist,
                                 std::vector<CovalentBondInfo::BondList> bondlist,
                                 std::vector<CovalentBondInfo::BondList> bondlist_idx,
                                 const CovalentBondParameterList& cbplist,
                                 std::vector< std::vector<int> > targetset,
                                 std::vector< std::pair<int,int> > ghostpairs,
                                 OperationSelector ope,
                                 MPICommunicator& _communicator,
                                 MPIReceiverForLong& _lrecv,
                                 MPISenderForLong& _lsend,
                                 PostProcess _postprocess,
                                 const ShortRange::CoulombType cltype,
                                 const LongRangeParameter longrangeparam,
                                 const int unitid, const int shortid,
                                 const SpaceVector<double> bs, 
				 const double bms,
                                 const double cutoff,
#ifdef USE_PAIRLIST
                                 const double pmargin,
#endif
                                 int num_t, int num_f,
                                 const MPI_Comm long_comm,
                                 const Config::TempCtrl& temp_control,
                                 bool expire_reverse
                                         );
#endif
/*
CalculationUnit& CalculationUnit::operator=(const CalculationUnit& cu)
{
  if(this != &cu){
    unit_identifier = cu.unit_identifier;
    num_total = cu.num_total;
    particlearray = cu.particlearray;
    setid = cu.setid;
    setid_to_index = cu.setid_to_index;
    typerangearray = cu.typerangearray;
    bondlistarray = cu.bondlistarray;
    bondlistarray_idx = cu.bondlistarray_idx;
    targetset = cu.targetset;
    boxsize = cu.boxsize;
    integrator = cu.integrator;
    calcforce = cu.calcforce;
    operations = cu.operations;
    postprocess = cu.postprocess;
    communicator = cu.communicator;
  }
  return *this;
}
*/
//! integrate
template<class PA, class GPA>
bool CalculationUnit<PA,GPA>::startTimeIntegration(double dt, long tmax, 
                                           const ShakeList& sl,
                                           int shake_type, int shake_max_iterate, double shake_tolerance,
                                           int reduce_interval, int printinter,
                                           int dump_crd_inter, int dump_rst_inter, Dump& dump, Dump& dump_rst,
                                           int moveinter,
                                           int& write_restart,
					   long& last_t)
{
  bool result;
  if(unit_identifier==0){
    std::cout << "max step " << tmax << std::endl;
  }
  result = integrator.time_integrations
    (
					particlearray,
				      bondlistarray,
				      bondlistarray_idx,
				      waterlist,
				      shakelist,
				      setid,setid_to_index,
				      calcforce,dt,tmax,operations,communicator,
				      longreceiver, longsender,
				      postprocess,
				      sl,
				      shake_type, shake_max_iterate, shake_tolerance,
				      reduce_interval, printinter, dump_crd_inter, dump_rst_inter, dump, dump_rst,
				      moveinter,
				      num_total,
				      num_freedom,
				      write_restart, last_t);
  return result;
}
#ifdef OLDPARTICLE
template
bool CalculationUnit<ParticleArray,ParticleArray>::startTimeIntegration(double dt, long tmax, 
                                           const ShakeList& sl,
                                           int shake_type, int shake_max_iterate, double shake_tolerance,
                                           int reduce_interval, int printinter,
                                           int dump_crd_inter, int dump_rst_inter, Dump& dump, Dump& dump_rst,
                                           int moveinter,
                                           int& write_restart,
                                                   long& last_t);
#else
template
bool CalculationUnit<CombinedParticleArray,GhostParticleArray>::startTimeIntegration(double dt, long tmax, 
                                           const ShakeList& sl,
                                           int shake_type, int shake_max_iterate, double shake_tolerance,
                                           int reduce_interval, int printinter,
                                           int dump_crd_inter, int dump_rst_inter, Dump& dump, Dump& dump_rst,
                                           int moveinter,
                                           int& write_restart,
                                                   long& last_t);
#endif

template<class PA, class GPA>
size_t CalculationUnit<PA,GPA>::getParticleArray(ParticleArray &pa)
{
  size_t n=0;
  for(size_t tr=0;tr<typerangearray.size();tr++){
    for(int i=typerangearray[tr].begin;i<typerangearray[tr].end;i++){
      pa.push_back(particlearray[i]);
      n++;
    }
  }
  return n;
}


