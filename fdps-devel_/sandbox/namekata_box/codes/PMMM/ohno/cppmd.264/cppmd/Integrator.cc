#ifdef DUMP_MEMSIZE
#include <unistd.h>
#include <fstream>
#include <sstream>
#endif
#include <cstring>
#include "Integrator.h"
#include "UnitParameter.h"

#ifdef SAMPLE_INTEGRATOR
# ifdef FPCOLL
#  include <fj_tool/fjsamp.h>
#  define SAMPLE_START() fpcoll_start()
#  define SAMPLE_STOP() fpcoll_stop()
# else
#  include <fj_tool/fipp.h>
#  define SAMPLE_START() fipp_start()
#  define SAMPLE_STOP() fipp_stop()
# endif
#endif
#ifndef SAMPLE_START
# define SAMPLE_START() 
# define SAMPLE_STOP() 
#endif
#ifdef FAPP
# include <fj_tool/fapp.h>
#endif
#ifdef K_PA
#include "fjcoll.h"
#endif

#define PRINT_DETAIL 0


#ifdef USE_CENTER_OF_MASS
#else
#define MODIFIED_ANDERSEN
#endif
//! integrate
/*! Input  : Mass, Charge, Atomtype, ID
  I/O    : Position,  Velocuty, Force, total Energy
*/

char timernames[12][256];
int timer_id_offset;
int counter_id[12];


template<class IPA, class IGPA>
Integrator<IPA, IGPA>::Integrator()
  :
#ifdef USE_HDF
  hdfdump(HDFDUMP),
  hdfrestore(HDFRESTORE),
#endif
#ifdef BINARY_DUMP
  binarydump(BDUMP),
  binaryrestore(BRESTORE),
#endif
  restore(false)
{
}

template<class IPA, class IGPA>
Integrator<IPA, IGPA>::Integrator(const int unitid,
                       const int shortid,
                       std::vector<TypeRange>& tra,
                       std::vector<CovalentBondInfo::BondList>& ba,
                       std::vector<CovalentBondInfo::BondList>& ba_idx,
                       const SpaceVector<double> bs,
		       const double bms,
                       const Config::TempCtrl& tc) 
    : deltat(1.0),
      unit_identifier(unitid),
      short_id(shortid),
      typerangearray(tra),
      bondlistarray(ba),
      bondlistarray_idx(ba_idx),
      boxsize(bs),
      boxsize_minimum_scale(bms),
      settle(),
      temp_control(tc),
      mti_bond(0),
      mti_short(0),
      mti_long(0),
      pre_calcforce(true),
#ifdef USE_HDF
      hdfdump(HDFDUMP),
      hdfrestore(HDFRESTORE),
#endif
#ifdef BINARY_DUMP
      binarydump(BDUMP),
      binaryrestore(BRESTORE),
#endif
      restore(false)
{
  initial_boxsize = boxsize;
  inv_initial_boxsize.x = 1.0/initial_boxsize.x;
  inv_initial_boxsize.y = 1.0/initial_boxsize.y;
  inv_initial_boxsize.z = 1.0/initial_boxsize.z;


  int ti=0;
  timer_id_offset = timer1;
  timer1 = ti++;
  strncpy(timernames[timer1],"Integrator",255);
  timer_move = ti++;
  strncpy(timernames[timer_move],"MoveParticle",255);
  timer_exp = ti++;
  strncpy(timernames[timer_exp],"ExchangeParticle",255);
#ifdef OVERLAP
  timer_expcalc = ti++;
  strncpy(timernames[timer_expcalc],"ExPCalc",255);
#endif
  timer_calc = ti++;
  strncpy(timernames[timer_calc],"calcForce",255);
  timer_exf = ti++;
  strncpy(timernames[timer_exf],"ExchangeForce",255);
  timer_rede = ti++;
  strncpy(timernames[timer_rede],"ReduceEnergy",255);
  timer_settle = ti++;
  strncpy(timernames[timer_settle],"SETTLE",255);
#ifdef USE_SHAKE
  timer_shake = ti++;
  strncpy(timernames[timer_shake],"SHAKE",255);
#endif
  timer_crd = ti++;
  strncpy(timernames[timer_crd],"DumpCRD",255);
  timer_rst = ti++;
  strncpy(timernames[timer_rst],"DumpRST",255);

#ifdef TIMER_DETAIL
  for(int t=0;t<ti;t++){
    counter_id[t]=PerfCounter::add_target(std::string(timernames[t]));
  }
#endif
}

inline
void timer_start(const int timer_id)
{
#ifdef TIMER_DETAIL
  PerfCounter::start(counter_id[timer_id]);
#endif
#ifdef FAPP
  fapp_start("Integrator",timer_id,1);
#endif
#ifdef K_PA
  start_collection(timernames[timer_id]);
#endif
}

inline
void timer_stop(const int timer_id)
{
#ifdef K_PA
  stop_collection(timernames[timer_id]);
#endif
#ifdef FAPP
  fapp_stop("Integrator",timer_id,1);
#endif
#ifdef TIMER_DETAIL
  PerfCounter::stop();
#endif
}

 
template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::mergerForces(CombinedParticleArray& particlearray) 
{
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      particlearray.force[i] += force[i];
    }
  }
}

template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::mergerForces(ParticleArray& particlearray) 
{
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      particlearray[i].force += force[i];
    }
  }
}

template<typename PA>
inline
Force sum_force(const PA& particlearray, const std::vector<TypeRange>& typerangearray)
{
  Force fc(0.0,0.0,0.0);
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      fc += getforce(particlearray,i);
    }
  }
  return fc;
}

/*
template<typename PA>
inline
void Integrator::reduce(double& ke, double& energy, Force& sf,
                        const double lke, const double le,
                        const PA& particlearray)
{
  double sepa[5];
  double redu[5];
  sepa[0] = lke;
  sepa[1] = le;
  Force lf = sum_force(particlearray);
  sepa[2] = lf.x;
  sepa[3] = lf.y;
  sepa[4] = lf.z;
  MPI_Allreduce(sepa, redu, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ke = redu[0];
  energy = redu[1];
  sf.x = redu[2];
  sf.y = redu[3];
  sf.z = redu[4];
}
*/

template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::reduce(double& ke, double& energy,
                        const double lke, const double le)
{
  double sepa[2];
  double redu[2];
  sepa[0] = lke;
  sepa[1] = le;
  MPI_Allreduce(sepa, redu, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ke = redu[0];
  energy = redu[1];
}

template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::reduce(double& ke, double& energy, double& virial,
                        const double lke, const double le, const double lv)
{
  double sepa[3];
  double redu[3];
  sepa[0] = lke;
  sepa[1] = le;
  sepa[2] = lv;
  MPI_Allreduce(sepa, redu, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ke = redu[0];
  energy = redu[1];
  virial = redu[2];
}

template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::reduce(double& ke, double& energy, double& virial, double &total_mass, Position& center_of_mass,
				   const double lke, const double le, const double lv, const double tm, const Position cm)
{
  double sepa[7];
  double redu[7];
  sepa[0] = lke;
  sepa[1] = le;
  sepa[2] = lv;
  sepa[3] = tm;
  sepa[4] = cm.x;
  sepa[5] = cm.y;
  sepa[6] = cm.z;
  MPI_Allreduce(sepa, redu, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ke = redu[0];
  energy = redu[1];
  virial = redu[2];
  total_mass = redu[3];
  center_of_mass.x = redu[4]/total_mass;
  center_of_mass.y = redu[5]/total_mass;
  center_of_mass.z = redu[6]/total_mass;
}

inline
void reduce_velocity(Velocity &v_sum, Velocity v_sum_local)
{
  double *vs = &(v_sum.x);
  double *vsl = &(v_sum_local.x);
  MPI_Allreduce(vsl, vs, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

inline
void reduce_momentum(SpaceVector<double>& m_sum, 
                                 SpaceVector<double>& am_sum,
                                 const SpaceVector<double> momentum, 
                                 const SpaceVector<double> angular_momentum)
{
  double m[6];
  double sum_m[6];
  m[0] = momentum.x;
  m[1] = momentum.y;
  m[2] = momentum.z;
  m[3] = angular_momentum.x;
  m[4] = angular_momentum.y;
  m[5] = angular_momentum.z;
  MPI_Allreduce(m, sum_m, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  m_sum.x = sum_m[0];
  m_sum.y = sum_m[1];
  m_sum.z = sum_m[2];
  am_sum.x = sum_m[3];
  am_sum.y = sum_m[4];
  am_sum.z = sum_m[5];
}

template<class PA>
void dump_momentum(const int unit_identifier,
                   const PA& particlearray,
                   const std::vector<TypeRange>& typerangearray)
{
  SpaceVector<double> momentum(0.0,0.0,0.0);
  SpaceVector<double> angular_momentum(0.0,0.0,0.0);
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      momentum += getmass(particlearray,i)*getvelocity(particlearray,i);
      angular_momentum += getmass(particlearray,i)*(getpos(particlearray,i) % getvelocity(particlearray,i));
    }
  }
  SpaceVector<double> m_sum;
  SpaceVector<double> am_sum;
  reduce_momentum(m_sum,am_sum,momentum,angular_momentum);
  if(unit_identifier==0){
    std::cout << "Momentum " << m_sum;
    std::cout << " Angular Momentum " << am_sum;
    std::cout << std::endl;
  }
}


void marge_ghostreturnindex(const std::vector<int>& ghost_longset_index,
                            const std::vector<TypeRange>& targettyperange,
                            std::vector<int>& ghostcbindexarray
                            )
{
          // TODO return ghost identified unit is cell for long but
          //      particle index is listed because it is specified by particle for CB.
          // it is uneffcient
          std::set<int> ghostreturnindex;
          {
            for(size_t ci=0;ci<ghost_longset_index.size();ci++){
              int c = ghost_longset_index[ci];
              for(int gi=targettyperange[c].begin;gi<targettyperange[c].end;gi++){
                ghostreturnindex.insert(gi);
              }
            }
            for(size_t gcbi=0;gcbi<ghostcbindexarray.size();gcbi++){
              ghostreturnindex.insert(ghostcbindexarray[gcbi]);
            }
            ghostcbindexarray.resize(ghostreturnindex.size());
            int i=0;
            for(std::set<int>::iterator gri = ghostreturnindex.begin();
                gri != ghostreturnindex.end();
                gri++)
            {
              ghostcbindexarray[i] = *gri;
              i++;
            }
          }
}

template<class PA> inline
void clear_force(PA& particlearray)
{
  for(size_t i=0;i<particlearray.size();i++){
    getforce(particlearray,i).x = 0.0;
    getforce(particlearray,i).y = 0.0;
    getforce(particlearray,i).z = 0.0;
  }
}

template<class PA> inline
void clear_force(PA& particlearray, std::vector<TypeRange>& typerangearray)
{
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getforce(particlearray,i).x = 0.0;
      getforce(particlearray,i).y = 0.0;
      getforce(particlearray,i).z = 0.0;
    }
  }
}

template<> inline
void clear_force(ForceArray& force, std::vector<TypeRange>& typerangearray)
{
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      force[i].x = 0.0;
      force[i].y = 0.0;
      force[i].z = 0.0;
    }
  }
}

/*!
  @brief top half of move particle, select move particle and move inside node
  @param [in,out] particlearray vector of Particle
  @param [in,out] typerangearray vector of typerange
  @param [in,out] bondlistarray vector of BondList
  @param [in,out] waterlist Waterlist
  @param [in,out] shakelist Shakelist
  @param [in,out] communicator MPI communicator, that has buffer for particles move out this node
  @param [in,out] postprocess This instance has methods moving particle
  @param [out] over_move  >0 means some particle move over margin
  @return not overflow cell
  @retval true none cell over flow
  @retval false some cell overflow
 */
template<class PA> inline
bool move_top_half(PA& particlearray,
                   std::vector<TypeRange>& typerangearray,
                   std::vector<CovalentBondInfo::BondList>& bondlistarray,
                   WaterList& waterlist,
		   ShakeList& shakelist,
                   MPICommunicator& communicator,
                   PostProcess& postprocess, 
                   int& over_move)
{
  communicator.all_move_out_particle.clear();
  communicator.all_move_out_cellid.clear();
  communicator.all_move_out_type.clear();
  communicator.all_move_out_bondlistarray.clear();
  bool not_over_flow;
  not_over_flow = postprocess.select_and_move_inside_node(particlearray,
							  typerangearray,
							  bondlistarray,
							  waterlist,
							  shakelist,
							  communicator.setid,
							  communicator.all_move_out_particle,
							  communicator.all_move_out_cellid,
							  communicator.all_move_out_type,
							  communicator.all_move_out_bondlistarray,
							  communicator.setid_to_index,
							  over_move
							  );
  return not_over_flow;
}

/*!
  @brief top half of move particle, select move particle and move inside node
  @param [in,out] particlearray vector of Particle
  @param [in,out] typerangearray vector of typerange
  @param [in,out] bondlistarray vector of BondList
  @param [in,out] waterlist Waterlist
  @param [in,out] shakelist Shakelist
  @param [in,out] communicator MPI communicator, that operate move particles inter node. This has buffer for particles move into this node
  @param [in,out] postprocess This instance has methods merge particles.
  @param [in] sl obsolete
  @return not overflow cell
  @retval true none cell over flow
  @retval false some cell overflow
 */
template<class PA> inline
bool move_bottom_half(PA& particlearray,
		      std::vector<TypeRange>& typerangearray,
		      std::vector<CovalentBondInfo::BondList>& bondlistarray,
		      WaterList& waterlist,
		      ShakeList& shakelist,
		      MPICommunicator& communicator,
		      PostProcess&postprocess,
		      const ShakeList& sl)
{
  //std::cout << " move_particle "  << std::endl;
  communicator.move_particle(particlearray,
                             typerangearray);
  //std::cout << "postprocess.merge "  << std::endl;
  return postprocess.merge(particlearray, 
			   communicator.setid, 
			   communicator.setid_to_index, 
			   typerangearray,
			   bondlistarray,
			   waterlist,
			   shakelist,
			   communicator.move_in_particle, 
			   communicator.move_in_cellid, 
			   communicator.move_in_type,
			   communicator.move_in_bondlistarray,
			   sl);
  
}

template<class IPA, class IGPA>
void Integrator<IPA, IGPA>::rollback_from_backup(IPA& particle,
				      std::vector<TypeRange>& typerangearray,
				      std::vector<CovalentBondInfo::BondList>& bondlistarray,
				      std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
				      WaterList& waterlist,
				      ShakeList& shakelist,
				      long& time)
{
  particle = backup_particle;
  typerangearray = backup_typerangearray;
  bondlistarray = backup_bondlistarray;
  bondlistarray_idx = backup_bondlistarray_idx;
  waterlist = backup_waterlist;
  shakelist = backup_shakelist;
  time = backup_time;
  ri = backup_ri;
  pi = backup_pi;
  dci = backup_dci;
  dri = backup_dri;
  eta_pos = backup_eta_pos;
  eta_vel = backup_eta_vel;
  eta_force = backup_eta_force;
  logv_pos = backup_logv_pos;
  logv_vel = backup_logv_vel;
  logv_force = backup_logv_force;
  boxsize = backup_boxsize;
  volume = backup_volume;
  tljcec = backup_tljcec;
}

inline
bool check_move_over(const int unit_identifier,
		     const bool not_over_flow,
		     const int over_move )
{
  bool move_over=false;

  int over_move_all;
  int over_flag;
  over_flag = over_move + ((not_over_flow)? 0:1);
  MPI_Allreduce((void *)(&over_flag),&over_move_all,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  // Roll back
  if(over_move_all>0){
    if(!not_over_flow){
      printf("cell over flow at node %d\n",unit_identifier);
    }
    if(over_move>0){
      printf("move over margin at node %d\n",unit_identifier);
    }
    if(unit_identifier==0){
      printf("cell over flow or move over of margin %d nodes\n", over_move_all);
    }
    move_over = true;
  }else{
    move_over = false;
  }
  return move_over;
}

/*!
  @brief Check particles move over margin and cells overflow. If move-over or overflow, roll-back last checkpoint, see backup_for_rollback()
  @param [in] not_over_flow >0 indicates over flow at this node, see  move_top_half(), move_bottom_half()
  @param [in] over_move >0 indicates some particles in this node move over margin, see move_top_half() 
  @param [in] max_moveinter maximum interval for checking move particle
  @param [in,out] particlearray vector of Particle
  @param [in,out] typerangearray vector of typerange
  @param [in,out] bondlistarray vector of BondList
  @param [in,out] bondlistarray_idx vector of BondList
  @param [in,out] waterlist 
  @param [in,out] shakelist
  @param [in,out] t current time step
  @param [in,out] moveinter interval for checking move particle
  @note some othre members of Integrator are also rollback
 */
template<class IPA, class IGPA>
inline
bool Integrator<IPA,IGPA>::check_move_over_and_rollback(const bool not_over_flow,
					      const int over_move,
					      const int max_moveinter,
					      IPA& particlearray,
					      std::vector<TypeRange>& typerangearray,
					      std::vector<CovalentBondInfo::BondList>& bondlistarray,
					      std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
					      WaterList& waterlist,
					      ShakeList& shakelist,
					      long& t,
					      int& moveinter)
{
  bool move_over=false;

  int over_move_all;
  int over_flag;
  over_flag = over_move + ((not_over_flow)? 0:1);
  MPI_Allreduce((void *)(&over_flag),&over_move_all,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  // Roll back
  if(over_move_all>0){
    if(!not_over_flow){
      printf("cell over flow at node %d\n",unit_identifier);
    }
    if(over_move>0){
      printf("move over margin at node %d\n",unit_identifier);
    }
    if(unit_identifier==0){
      printf("cell over flow or move over of margin %d nodes\n", over_move_all);
    }
    moveinter = moveinter/2;
    if(t==0){
      if(unit_identifier==0){
        printf("Incorrect move at 1st step\n");
        printf("Give up roll back.\n");
      }
      moveinter = 0;
      move_over = true;
    }else if(moveinter==0){
      if(unit_identifier==0){
        printf("Move interval already shrink 1.\n");
        printf("Give up roll back. Continue with original interval %d, calculation may be incorrect.\n",max_moveinter);
      }
      moveinter = 0;
      rollback_from_backup(particlearray,
			   typerangearray,
			   bondlistarray,
			   bondlistarray_idx,
			   waterlist, 
			   shakelist,
			   t);
      move_over = true;
    }else{
      if(unit_identifier==0){
        printf("Roll Back to %ld step\n",backup_time);
        printf("Move interval change %d\n",moveinter);
      }
      rollback_from_backup(particlearray,
			   typerangearray,
			   bondlistarray,
			   bondlistarray_idx,
			   waterlist, 
			   shakelist,
			   t);
      move_over = true;
    }
  }else{
    move_over = false;
  }
  return move_over;
}

#ifdef USE_HDF
template<class IPA, class IGPA>
int Integrator<IPA, IGPA>::open_hdfdump(const std::string& hdfname)
{
  std::string namebase;
  std::string sufix(".hdf");
  std::string::size_type posdot = hdfname.find_last_of('.');
  if(posdot!=std::string::npos){
    if(posdot<hdfname.length()){
      char c = hdfname.at(posdot+1);
      if((c=='H')||(c=='h')){
	namebase = hdfname.substr(0,posdot);
	sufix = hdfname.substr(posdot);
      }else{
	namebase = hdfname;
      }
    }else{
      namebase = hdfname.substr(0,posdot);
    }
  }else{
    namebase = hdfname;
  }
  char name[256];
  snprintf(name,256,"%s.%05d",namebase.c_str(),short_id);
  std::string filename(name);
  filename.append(sufix);
  hdfdump.init(filename);
  return 0;
}

template<class IPA, class IGPA>
void Integrator<IPA, IGPA>::close_hdfdump()
{
  hdfdump.close();
}

template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::dump_hdf_snapshot(const IPA& particlearray,
					      const std::vector<TypeRange>& typerangearray,
					      const std::vector<CovalentBondInfo::BondList>& bondlistarray,
					      const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
					      const WaterList& waterlist,
					      const ShakeList& shakelist,
					      const double kenergy,
					      const double penergy,
					      const double virial,
					      const double total_penergy,
					      const long t,
					      const long num_total,
					      const int num_rank)
{
  hdfdump.setParameter(deltat,t,ri,pi,dci,dri,num_total);
  std::string ver("alpha");
  HDFSnapshot::HDFSnapshotFile::SnapshotType sst=HDFSnapshot::HDFSnapshotFile::ParallelDump;
  hdfdump.setInfo(short_id,num_rank,&ver,num_total,sst);
  hdfdump.setIntegration(tljcec, kenergy, penergy, total_penergy, virial,
			 eta_pos, eta_vel, eta_force, 
			 logv_pos, logv_vel, logv_force, volume);
  int boxtype = 0;
  SpaceVector<double> angle(90.0,90.0,90.0); 
  hdfdump.setBoxDatatype(boxsize, boxtype, angle);
  hdfdump.setParticle(particlearray,typerangearray);
  hdfdump.setBondlists(bondlistarray);
  hdfdump.setWaterlist(waterlist);
  hdfdump.setShakelist(shakelist);
  hdfdump.dump();
}


template<class IPA, class IGPA>
int Integrator<IPA, IGPA>::open_hdfrestore(const std::string& hdfname)
{
  std::string namebase;
  std::string sufix(".hdf");
  std::string::size_type posdot = hdfname.find_last_of('.');
  if(posdot!=std::string::npos){
    if(posdot<hdfname.length()){
      char c = hdfname.at(posdot+1);
      if((c=='H')||(c=='h')){
	namebase = hdfname.substr(0,posdot);
	sufix = hdfname.substr(posdot);
      }else{
	namebase = hdfname;
      }
    }else{
      namebase = hdfname.substr(0,posdot);
    }
  }else{
    namebase = hdfname;
  }
  char name[256];
  snprintf(name,256,"%s.%05d",namebase.c_str(),short_id);
  std::string filename(name);
  filename.append(sufix);
  return hdfrestore.init(filename);
}

template<class IPA, class IGPA>
void Integrator<IPA, IGPA>::close_hdfrestore()
{
  hdfrestore.close();
}

template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::restore_hdf_snapshot(IPA& particlearray,
						 std::vector<TypeRange>& typerangearray,
						 std::vector<CovalentBondInfo::BondList>& bondlistarray,
						 std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
						 WaterList& waterlist,
						 ShakeList& shakelist,
						 double& kenergy,
						 double& penergy,
						 double& virial,
						 double& total_penergy,
						 long& t,
						 const long num_total,
						 const int num_rank)
{
  hdfrestore.restore();

  int boxtype;
  SpaceVector<double> angle; 
  hdfrestore.getBoxDatatype(boxsize, boxtype, angle);

  hdfrestore.getParticle(particlearray,typerangearray);
  hdfrestore.getBondlists(bondlistarray);
  waterlist.clear();
  waterlist.reverse_list.clear();
  hdfrestore.getWaterlist(waterlist);
  shakelist.clear();
  shakelist.reverse_list.clear();
  hdfrestore.getShakelist(shakelist);

  //  std::string ver("alpha");
  //  HDFSnapshot::HDFSnapshotFile::SnapshotType sst=HDFSnapshot::HDFSnapshotFile::ParallelDump;
  //  hdfrestore.getInfo(short_id,num_rank,&ver,num_total,sst);
  hdfrestore.getIntegration(tljcec, kenergy, penergy, total_penergy, virial,
			    eta_pos, eta_vel, eta_force, 
			    logv_pos, logv_vel, logv_force, volume);
  double dt;
  AtomID nt;
  int ti;
  hdfrestore.getParameter(dt,ti,ri,pi,dci,dri,nt);
  t = ti;
}
#endif

#ifdef BINARY_DUMP

template<class IPA, class IGPA>
int Integrator<IPA, IGPA>::open_binarydump(const std::string& filename)
{
  return  binarydump.init(filename);
}
template<class IPA, class IGPA>
int Integrator<IPA, IGPA>::open_binaryrestore(const std::string& filename)
{
  return binaryrestore.init(filename);
}
template<class IPA, class IGPA>
void Integrator<IPA, IGPA>::close_binarydump()
{
  binarydump.close();
}
template<class IPA, class IGPA>
void Integrator<IPA, IGPA>::close_binaryrestore()
{
  binaryrestore.close();
}


template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::dump_binary_snapshot(const IPA& particlearray,
						 const std::vector<TypeRange>& typerangearray,
						 const std::vector<CovalentBondInfo::BondList>& bondlistarray,
						 const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
						 const WaterList& waterlist,
						 const ShakeList& shakelist,
						 const double kenergy,
						 const double penergy,
						 const double virial,
						 const double total_penergy,
						 const long t)
{
  binarydump.dump_binary_basic(particlearray, typerangearray, bondlistarray, bondlistarray_idx, waterlist, shakelist, t);
  double od[] = {kenergy,
		 penergy,
		 virial,
		 total_penergy,
		 eta_pos, 
		 eta_vel, 
		 eta_force,
		 logv_pos,
		 logv_vel,
		 logv_force,
		 boxsize.x,
		 boxsize.y,
		 boxsize.z,
		 volume,
		 tljcec
  };
  binarydump.dump_binary_optional_double(od, 15);
  int oi[] = {ri,
	      pi,
	      dci,
	      dri
  };
  binarydump.dump_binary_optional_int(oi, 4);
}

template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::restore_binary_snapshot(IPA& particlearray,
						    std::vector<TypeRange>& typerangearray,
						    std::vector<CovalentBondInfo::BondList>& bondlistarray,
						    std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
						    WaterList& waterlist,
						    ShakeList& shakelist,
						    double &kenergy,
						    double &penergy,
						    double &virial,
						    double& total_penergy,
						    long& t)
{
  binaryrestore.restore_binary_basic(particlearray, typerangearray, bondlistarray, bondlistarray_idx, waterlist, shakelist, t);
  double od[15];
  binaryrestore.restore_binary_optional_double(od, 15);
  kenergy = od[0];
  penergy = od[1];
  virial = od[2];
  total_penergy = od[3];
  eta_pos = od[4];
  eta_vel = od[5];
  eta_force = od[6];
  logv_pos = od[7];
  logv_vel = od[8];
  logv_force = od[9];
  boxsize.x = od[10];
  boxsize.y = od[11];
  boxsize.z = od[12];
  volume = od[13];
  tljcec = od[14];
  int oi[4];
  binaryrestore.restore_binary_optional_int(oi, 4);
  ri = oi[0];
  pi = oi[1];
  dci = oi[2];
  dri = oi[3];
  //  printf("%d %d %d %d\n",ri,pi,dci,dri);
}
#endif

/*!
  @brief Make check point. Back up current data.
  @param [in] particlearray vector of Particle
  @param [in] typerangearray vector of typerange
  @param [in] bondlistarray vector of BondList
  @param [in] bondlistarray_idx vector of BondList
  @param [in] waterlist 
  @param [in] shakelist
  @param [in] t current time step
  @note some othre members of Integrator are also backup
 */
template<class IPA, class IGPA>
inline
void Integrator<IPA, IGPA>::backup_for_rollback(const IPA& particlearray,
				     const std::vector<TypeRange>& typerangearray,
				     const std::vector<CovalentBondInfo::BondList>& bondlistarray,
				     const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
				     const WaterList& waterlist,
				     const ShakeList& shakelist,
				     const long t)
{
  backup_particle = particlearray;
  backup_typerangearray = typerangearray;
  backup_bondlistarray = bondlistarray;
  backup_bondlistarray_idx = bondlistarray_idx;
  backup_waterlist = waterlist;
  backup_shakelist = shakelist;
  backup_time = t;
  backup_ri = ri;
  backup_pi = pi;
  backup_dci = dci;
  backup_dri = dri;
  backup_eta_pos = eta_pos;
  backup_eta_vel = eta_vel;
  backup_eta_force = eta_force;
  backup_logv_pos = logv_pos;
  backup_logv_vel = logv_vel;
  backup_logv_force = logv_force;
  backup_boxsize = boxsize;
  backup_volume = volume;
  backup_tljcec = tljcec;
}

template<class PA> inline
void dump_bonds(const int unit_identifier,
                const MPICommunicator& communicator,
                const std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  if((DebugLog::verbose>1)){
    int nr;
    MPI_Comm_size(MPI_COMM_WORLD,&nr);
    for(int i=0;i<nr;i++){
      if(unit_identifier == i){
        std::cout << "rank " << i << " ";
        dump_bondlistarray(communicator.all_move_out_bondlistarray);
        dump_bondlistarray(bondlistarray);
        dump_bondlistarray(communicator.move_in_bondlistarray);
        std::cout << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

template<class PA> inline
void dump_particle_and_bond(const int unit_identifier,
                            const MPICommunicator& communicator,
                            const PA& particlearray,
                            const std::vector<TypeRange>& typerangearray,
                            const std::vector<int>& setid,
                            const double energy,
                            const std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  if((DebugLog::verbose>1)){
    int nr;
    MPI_Comm_size(MPI_COMM_WORLD,&nr);
    for(int i=0;i<nr;i++){
      if(unit_identifier == i){
        std::cout << "rank " << i << " ";
        dump_particle(particlearray,typerangearray,setid,energy);
        dump_bondlistarray(communicator.all_move_out_bondlistarray);
        dump_bondlistarray(bondlistarray);
        dump_bondlistarray(communicator.move_in_bondlistarray);
            
        std::cout << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

template<class PA, class GPA> inline
void dump_atomids(const PA& particlearray, const std::vector<TypeRange>& typerangearray, const std::vector<int>& setid,
                  const GPA& ghost, const std::vector<TypeRange>& targettyperange, const std::vector<int>& recvsetid)
{
  if(DebugLog::verbose>1){
    if(number_of_particle(typerangearray)){
      std::cout << "self ";
      dump_atomid(particlearray,typerangearray,setid);
    }
    if(number_of_particle(targettyperange)){
      std::cout << "ghost ";
      dump_atomid(ghost,targettyperange,recvsetid);
    }
  }
}

template<class PA>
void dump_a_particle(const PA& spa,
		   const std::vector<TypeRange>& typerangearray,
		   const int atomid,
		   const int step,
		   const int myrank)
{
  for(size_t si=0;si<typerangearray.size();si++){
    for(int i=typerangearray[si].begin;i<typerangearray[si].end;i++){
      if(getatomid(spa,i)==atomid){
	Position pos = getpos(spa,i);
	std::cout << "ATOMID " << atomid << " " << pos << " step " << step << " index " << i << " cell " << si << "[" << typerangearray[si].begin << ":" << typerangearray[si].end << "] in rank " << myrank << std::endl;
      }
    }
  }
}


template<class IPA, class IGPA>
inline
void Integrator<IPA,IGPA>::dump_detail(const IPA& particlearray,
                             const WaterList& waterlist,
                             const long t)
{
  for(int uid=0; uid<8; uid++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(unit_identifier == uid){
      printf("\nt = %ld(%d)\n",t,unit_identifier);
      printf("particle (in Integrator)\n");
      printf("water_list(%d)\n",unit_identifier);
      for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++)
        printf("  %d :%d %d\n",it->first, it->second.h1, it->second.h2);
      for(size_t s=0;s<typerangearray.size();s++){
        if(typerangearray[s].begin < typerangearray[s].end){
          printf("s=%ld(%d)\n",s,unit_identifier);
          TypeRange tr = typerangearray[s];
          printf(" (%d, %d) ",tr.begin, tr.end);
          printf(" LJ(%d, %d) ",tr.lj.begin, tr.lj.end);
          printf(" LJC(%d, %d) ",tr.ljcoulomb.begin, tr.ljcoulomb.end);
          printf(" C(%d, %d)\n",tr.coulomb.begin, tr.coulomb.end);
          for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
            printf(" %d (%d) [%d] ",i,getatomtype(particlearray,i), getatomid(particlearray,i));
            printf("position:");
            std::cout << getpos(particlearray,i) << std::endl;
            printf("            velocity:");
            std::cout << getvelocity(particlearray,i) << std::endl;
          }
        }
      }
      fflush(stdout);
      std::cout << std::flush;
      fflush(stdout);
      std::cout << std::flush;
    }
  }
}

template<class PA, class GPA> inline
void exchangeParticleArray_top_half(const PA& particlearray, 
				    const std::vector<TypeRange>& typerangearray,
				    const std::vector<CovalentBondInfo::BondList>& bondlistarray,
				    const std::vector<int>& setid,
				    GPA& ghost,
				    std::vector<ParticleRange>& targetparticlerange,
				    std::vector<int>& recvsetid,
				    std::map<int,int>& recvsetid_to_index,
				    std::vector<TypeRange>& targettyperange,
				    std::vector<CovalentBondInfo::BondList>& targetbond,
				    std::vector<int>& ghostcbindexarray,
				    unsigned int& num_ghost,
				    CalcForce& calcforce,
				    MPICommunicator& communicator,
				    MPIReceiverForLong& longreceiver,
				    MPISenderForLong& longsender,
				    bool& reset_ghost_longset_index,
				    const OperationSelector& operations
				    )
{
  
  communicator.exchangeParticleArraysubset_top_half(particlearray, 
                                           typerangearray, 
                                           bondlistarray, 
                                           ghost,
                                           targetparticlerange, 
                                           recvsetid, 
                                           recvsetid_to_index,
                                           targettyperange, 
                                           targetbond);
#ifdef CPPMD_ENABLE_LONGRANGE
  longreceiver.transferParticle();
  longsender.exchangeParticleArraysubset(particlearray, typerangearray);
#endif  // CPPMD_ENABLE_LONGRANGE
}
template<class PA, class GPA> inline
void exchangeParticleArray_bottom_half(const PA& particlearray, 
				       const std::vector<TypeRange>& typerangearray,
				       const std::vector<CovalentBondInfo::BondList>& bondlistarray,
				       const std::vector<int>& setid,
				       GPA& ghost,
				       std::vector<ParticleRange>& targetparticlerange,
				       std::vector<int>& recvsetid,
				       std::map<int,int>& recvsetid_to_index,
				       std::vector<TypeRange>& targettyperange,
				       std::vector<CovalentBondInfo::BondList>& targetbond,
				       std::vector<int>& ghostcbindexarray,
				       unsigned int& num_ghost,
				       CalcForce& calcforce,
				       MPICommunicator& communicator,
				       MPIReceiverForLong& longreceiver,
				       MPISenderForLong& longsender,
				       bool& reset_ghost_longset_index,
				       const OperationSelector& operations
				       )
{
  
  communicator.exchangeParticleArraysubset_bottom_half(particlearray, 
                                           typerangearray, 
                                           bondlistarray, 
                                           ghost,
                                           targetparticlerange, 
                                           recvsetid, 
                                           recvsetid_to_index,
                                           targettyperange, 
                                           targetbond);
#ifdef CPPMD_ENABLE_LONGRANGE
  longreceiver.exchangeParticleArraysubset(ghost,
                                           targetparticlerange,
                                           recvsetid,
                                           recvsetid_to_index,
                                           targettyperange);
  if(reset_ghost_longset_index){
    calcforce.ghost_longset_index.clear();
    calcforce.set_ghost_longset_index(recvsetid_to_index);
    reset_ghost_longset_index = false;
  }
#endif  // CPPMD_ENABLE_LONGRANGE
  ghostcbindexarray.clear();
  calcforce.ghostshortforce.resize(ghost.size());
  calcforce.covalentbond.map_particle(const_cast<PA&>(particlearray),calcforce.shortforce,
                                      typerangearray,
                                      ghost,calcforce.ghostshortforce,
                                      targettyperange,
                                      calcforce.bondlist,
                                      ghostcbindexarray);
  dump_atomids(particlearray,typerangearray,setid,
               ghost,targettyperange,recvsetid);
  num_ghost = ghost.size();
}
//! TODO separatable long
/*!
  @brief exchange particles for particle move step, include charge, atomtype, typerange
  @param [in] particlearray  Particles in this node
  @param [in] typerangearray typeranges of this node
  @param [in] bondlistarray bondlists of this node
  @param [in] setid identifiers of cell in this node
  @param [out] ghost ghost particles, recievd particles from othen node
  @param [out] targetparticlerange particlerange for ghost
  @param [out] recvsetid identifiers of cell of ghost
  @param [out] recvsetid_to_index map of recvsetid to index of cell of ghost
  @param [out] targettyperange typeranges of ghost
  @param [out] targetbond bond for ghost, not used
  @param [out] ghostcbindexarray indexes of ghost connected covalent bond
  @param [out] num_ghost number of ghost particle
  @param [in,out] calcforce instance of force calculation
  @param [in,out] communicator MPI communicator, that send/receive paticles, etc
  @param [in,out] longreceiver MPI communicator for recieve particles for long-range operation.
  @param [in,out] longsender MPI communicator for send particles for long-range operation.
  @param [in,out] reset_ghost_longset_index necessity of clearing longset index
 */
template<class PA, class GPA> inline
void exchangeParticleArray(const PA& particlearray, 
                           const std::vector<TypeRange>& typerangearray,
                           const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                           const std::vector<int>& setid,
                           GPA& ghost,
                           std::vector<ParticleRange>& targetparticlerange,
                           std::vector<int>& recvsetid,
                           std::map<int,int>& recvsetid_to_index,
                           std::vector<TypeRange>& targettyperange,
                           std::vector<CovalentBondInfo::BondList>& targetbond,
                           std::vector<int>& ghostcbindexarray,
                           unsigned int& num_ghost,
                           CalcForce& calcforce,
                           MPICommunicator& communicator,
                           MPIReceiverForLong& longreceiver,
                           MPISenderForLong& longsender,
                           bool& reset_ghost_longset_index,
			   const OperationSelector& operations
                           )
{
  //  std::cout << " exchangeParticleArraysubset  "  << std::endl;
  communicator.exchangeParticleArraysubset(particlearray, 
                                           typerangearray, 
                                           bondlistarray, 
                                           ghost,
                                           targetparticlerange, 
                                           recvsetid, 
                                           recvsetid_to_index,
                                           targettyperange, 
                                           targetbond);
#ifdef CPPMD_ENABLE_LONGRANGE
  longreceiver.transferParticle();
  longsender.exchangeParticleArraysubset(particlearray, typerangearray);
  longreceiver.exchangeParticleArraysubset(ghost,
                                           targetparticlerange,
                                           recvsetid,
                                           recvsetid_to_index,
                                           targettyperange);
  if(reset_ghost_longset_index){
    calcforce.ghost_longset_index.clear();
    calcforce.set_ghost_longset_index(recvsetid_to_index);
    reset_ghost_longset_index = false;
  }
#endif  // CPPMD_ENABLE_LONGRANGE
  ghostcbindexarray.clear();
  calcforce.ghostshortforce.resize(ghost.size());
  calcforce.covalentbond.map_particle(const_cast<PA&>(particlearray),calcforce.shortforce,
                                      typerangearray,
                                      ghost,calcforce.ghostshortforce,
                                      targettyperange,
                                      calcforce.bondlist,
                                      ghostcbindexarray);
  dump_atomids(particlearray,typerangearray,setid,
               ghost,targettyperange,recvsetid);
  num_ghost = ghost.size();
}

template<class PA, class GPA> inline
void exchangeParticleArray_onlyPosition_top_half(const PA& particlearray, 
						 const std::vector<TypeRange>& typerangearray,
						 const std::vector<CovalentBondInfo::BondList>& bondlistarray,
						 GPA& ghost,
						 const std::vector<ParticleRange>& targetparticlerange,
						 const std::vector<int>& recvsetid,
						 const std::map<int,int>& recvsetid_to_index,
						 const std::vector<TypeRange>& targettyperange,
						 const std::vector<CovalentBondInfo::BondList>& targetbond,
						 const unsigned int& num_ghost,
						 MPICommunicator& communicator,
						 MPIReceiverForLong& longreceiver,
						 MPISenderForLong& longsender,
						 const OperationSelector& operations
						 )
{
  
  communicator.exchangeParticleArraysubset_onlyPosition_top_half(particlearray, 
                                                        typerangearray, 
                                                        bondlistarray, 
                                                        ghost,
                                                        targetparticlerange, 
                                                        recvsetid, 
                                                        recvsetid_to_index,
                                                        targettyperange, 
                                                        targetbond);
#ifdef CPPMD_ENABLE_LONGRANGE
  if(operations.doShortLongCommunication){
    longreceiver.transferParticle();
    longsender.exchangeParticleArraysubset_onlyPosition_top_half(particlearray, typerangearray);
  }
#endif  // CPPMD_ENABLE_LONGRANGE
}

template<class PA, class GPA> inline
void exchangeParticleArray_onlyPosition_bottom_half(const PA& particlearray, 
						    const std::vector<TypeRange>& typerangearray,
						    const std::vector<CovalentBondInfo::BondList>& bondlistarray,
						    GPA& ghost,
						    const std::vector<ParticleRange>& targetparticlerange,
						    const std::vector<int>& recvsetid,
						    const std::map<int,int>& recvsetid_to_index,
						    const std::vector<TypeRange>& targettyperange,
						    const std::vector<CovalentBondInfo::BondList>& targetbond,
						    const unsigned int& num_ghost,
						    MPICommunicator& communicator,
						    MPIReceiverForLong& longreceiver,
						    MPISenderForLong& longsender,
						    const OperationSelector& operations
						    )
{
  
  communicator.exchangeParticleArraysubset_onlyPosition_bottom_half(particlearray, 
                                                        typerangearray, 
                                                        bondlistarray, 
                                                        ghost,
                                                        targetparticlerange, 
                                                        recvsetid, 
                                                        recvsetid_to_index,
                                                        targettyperange, 
                                                        targetbond);
  if(ghost.size()!=num_ghost){
    printf("ghost size change\n");
  }
#ifdef CPPMD_ENABLE_LONGRANGE
  if(operations.doShortLongCommunication){
    longreceiver.exchangeParticleArraysubset_onlyPosition(ghost,
							  targetparticlerange,
							  recvsetid,
							  recvsetid_to_index,
							  targettyperange);
  }
#endif  // CPPMD_ENABLE_LONGRANGE
}

/*!
  @brief exchange positons of particles for particle not move step
  @param [in] particlearray  Particles in this node
  @param [in] typerangearray typeranges of this node
  @param [in] bondlistarray bondlists of this node
  @param [out] ghost ghost particles, recievd particles from othen node
  @param [in] targetparticlerange particlerange for ghost
  @param [in] recvsetid identifiers of cell of ghost
  @param [in] recvsetid_to_index map of recvsetid to index of cell of ghost
  @param [in] targettyperange typeranges of ghost
  @param [in] targetbond bond for ghost, not used
  @param [in] num_ghost number of ghost particle
  @param [in,out] communicator MPI communicator, that send/receive paticles, etc
  @param [in,out] longreceiver MPI communicator for recieve particles for long-range operation.
  @param [in,out] longsender MPI communicator for send particles for long-range operation.
 */
template<class PA, class GPA> inline
void exchangeParticleArray_onlyPosition(const PA& particlearray, 
                                        const std::vector<TypeRange>& typerangearray,
                                        const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                        GPA& ghost,
                                        const std::vector<ParticleRange>& targetparticlerange,
                                        const std::vector<int>& recvsetid,
                                        const std::map<int,int>& recvsetid_to_index,
                                        const std::vector<TypeRange>& targettyperange,
                                        const std::vector<CovalentBondInfo::BondList>& targetbond,
                                        const unsigned int& num_ghost,
                                        MPICommunicator& communicator,
                                        MPIReceiverForLong& longreceiver,
                                        MPISenderForLong& longsender,
					const OperationSelector& operations
                                        )
{
  
  communicator.exchangeParticleArraysubset_onlyPosition(particlearray, 
                                                        typerangearray, 
                                                        bondlistarray, 
                                                        ghost,
                                                        targetparticlerange, 
                                                        recvsetid, 
                                                        recvsetid_to_index,
                                                        targettyperange, 
                                                        targetbond);
  if(ghost.size()!=num_ghost){
    printf("ghost size change\n");
  }
#ifdef CPPMD_ENABLE_LONGRANGE
  if(operations.doShortLongCommunication){
    longreceiver.transferParticle();
    longsender.exchangeParticleArraysubset_onlyPosition(particlearray, typerangearray);
    longreceiver.exchangeParticleArraysubset_onlyPosition(ghost,
							  targetparticlerange,
							  recvsetid,
							  recvsetid_to_index,
							  targettyperange);
  }
#endif  // CPPMD_ENABLE_LONGRANGE
}

template<class GPA> inline
void clear_before_move(GPA& ghost,
                       std::vector<TypeRange>& targettyperange,
                       std::vector<int>& recvsetid,
                       std::vector<CovalentBondInfo::BondList>& bondlistarray,
                       CalcForce& calcforce)
{
  ghost.clear();
  targettyperange.clear();
  recvsetid.clear();
  calcforce.bondlist.clear();
  calcforce.bondlist.merge_bondlistarray(bondlistarray);
  calcforce.bondlist.make_atomidset();
}

template<class GPA> inline
void postprocess_after_exchangeparticle(PostProcess& postprocess,
                                        GPA& ghost,
                                        const std::vector<int>& recvsetid,
                                        const std::map<int,int>& recvsetid_to_index,
                                        const std::vector<TypeRange>& targettyperange,
                                        ForceArray& ghostforce,
                                        double& energy,
					bool& shiftidcheck)
{
  if(ghost.size()>0){
    if(shiftidcheck){
      postprocess.postprocess_receive_with_shiftidcheck(ghost,recvsetid,recvsetid_to_index,targettyperange);
      shiftidcheck = false;
    }else{
      postprocess.postprocess_receive(ghost,recvsetid,recvsetid_to_index,targettyperange);
    }
  }
  clear_force(ghost);
  ghostforce.resize(ghost.size());
  energy = 0.0;
}

inline
void exchangeForceArray(ForceArray& ghostforce,
                        std::vector<ParticleRange>& targetparticlerange,
                        ForceArray& force,
                        const std::vector<int>& ghost_longset_index,
                        const std::vector<TypeRange>& targettyperange,
                        std::vector<int>& ghostcbindexarray,
                        MPICommunicator& communicator,
                        MPIReceiverForLong& longreceiver,
                        MPISenderForLong& longsender,
                        const OperationSelector& operations,
                        const int moved)
{
  if(operations.doExchangeForceIndexed){
    if(moved){
      if(operations.doLongrangecalculation){
        marge_ghostreturnindex(ghost_longset_index,
                               targettyperange,
                               ghostcbindexarray);
      }
      communicator.exchangeForceArraysubset_with_index(ghostforce, 
                                                       targetparticlerange, 
                                                       ghostcbindexarray,
                                                       force);
    }else{
      communicator.exchangeForceArraysubset_indexed(ghostforce, 
                                                    targetparticlerange, 
                                                    ghostcbindexarray,
                                                    force);
    }
  }else{
    communicator.exchangeForceArraysubset(ghostforce, targetparticlerange, 
                                          force);
  }
#ifdef CPPMD_ENABLE_LONGRANGE
  if(operations.doShortLongCommunication){
    longsender.transferForce();
    longreceiver.exchangeForceArraysubset(ghostforce, targetparticlerange);
    longsender.exchangeForceArraysubset(force);
  }
#endif  // CPPMD_ENABLE_LONGRANGE
}

template<class PA> inline
void calc_local_kenergy(const PA& particlearray,
                        const std::vector<TypeRange>& typerangearray,
                        double& ke)
{
  ke = 0.0;
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      ke += getmass(particlearray,i)*getvelocity(particlearray,i).norm2();
    }
  }
  ke *= 0.5;
}

template<class PA> inline
void calc_local_kenergy_virial(PA& particlearray,
			       const std::vector<TypeRange>& typerangearray,
			       const double pairwise_virial,
			       const double water_settle_virial,
			       const double shake_virial,
			       double& ke, double& lv)
{
  ke=0.0;
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      ke += getmass(particlearray,i)*getvelocity(particlearray,i).norm2();
    }
  }
  ke*=0.5;
  lv = pairwise_virial + water_settle_virial + shake_virial;
  //  printf("virial %e + %e + %e = %e\n",pairwise_virial,water_settle_virial,shake_virial,lv);
}

template<class PA> inline
void calc_total_energy(const PA& particlearray,
		       const std::vector<TypeRange>& typerangearray,
		       const double tljcec,
		       const double local_pe,
		       double& local_ke,
		       double& total_pe,
		       double& total_ke)
{
  double ke = 0.0;
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      ke += getmass(particlearray,i)*getvelocity(particlearray,i).norm2();
    }
  }
  double local_kpe[2];
  double total_kpe[2];
  local_kpe[0] = ke*0.5;
  local_kpe[1] = local_pe;
  MPI_Allreduce(local_kpe, total_kpe, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  total_ke = total_kpe[0];
  total_pe = total_kpe[1] + tljcec;
}
template<class PA> inline
void calc_total_energy(const PA& particlearray,
		       const std::vector<TypeRange>& typerangearray,
		       const double tljcec,
		       const double local_pe,
		       double& local_ke,
		       double& total_pe,
		       double& total_ke,
		       Force& force)
{
  double ke = 0.0;
  Force fc(0.0,0.0,0.0);
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      ke += getmass(particlearray,i)*getvelocity(particlearray,i).norm2();
      fc += getforce(particlearray,i);
    }
  }
  double local_kpef[5];
  double total_kpef[5];
  local_kpef[0] = ke*0.5;
  local_kpef[1] = local_pe;
  local_kpef[2] = fc.x;
  local_kpef[3] = fc.y;
  local_kpef[4] = fc.z;
  MPI_Allreduce(local_kpef, total_kpef, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  total_ke = total_kpef[0];
  total_pe = total_kpef[1] + tljcec;
  force.x = total_kpef[2];
  force.y = total_kpef[3];
  force.z = total_kpef[4];
}

inline
void calculate_pressure(const double volume, const double ke_sum,
			const double virial_sum,
			double &pressure)
{
  pressure = (2.0*ke_sum+virial_sum)/(3.0*volume);
#ifdef SIMPLE_CUTOFF
  pressure += ShortRange::LJCutoffEnergyCorrectionV/(volume*volume);
#endif
  /*
    potential energy of cutoff is defined at ShortRangeInteraction.cc
    ShortRange::LJCutoffEnergyCorrectionV/volume;
   */
}

template<class PA, class PVA> inline
void backup_water_posvel(const PA& particlearray,
                         const WaterList& waterlist,
                         PVA& posvelarray)
{
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++){
    getpos(posvelarray,it->first)          =  getpos(particlearray,it->first);
    getvelocity(posvelarray,it->first)     =  getvelocity(particlearray,it->first);
    getpos(posvelarray,it->second.h1)      =  getpos(particlearray,it->second.h1);
    getvelocity(posvelarray,it->second.h1) =  getvelocity(particlearray,it->second.h1);
    getpos(posvelarray,it->second.h2)      =  getpos(particlearray,it->second.h2);
    getvelocity(posvelarray,it->second.h2) =  getvelocity(particlearray,it->second.h2);
  }
}

template<class PA, class PVA> inline
void calc_water_settle_force(const PA& particlearray,
			     const WaterList& waterlist,
			     const PVA& posvelarray,
			     const double deltat,
			     ForceArray& water_settle_force)
{
  double inv_dt = 1.0/deltat;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++){
    water_settle_force[it->first] += (getvelocity(particlearray,it->first) - getvelocity(posvelarray,it->first))*getmass(particlearray,it->first)*inv_dt;
    water_settle_force[it->second.h1] += (getvelocity(particlearray,it->second.h1) - getvelocity(posvelarray,it->second.h1))*getmass(particlearray,it->second.h1)*inv_dt;
    water_settle_force[it->second.h2] += (getvelocity(particlearray,it->second.h2) - getvelocity(posvelarray,it->second.h2))*getmass(particlearray,it->second.h2)*inv_dt;
  }
}

template<class PA, class PVA> inline
void backup_shake_posvel(const PA& particlearray,
                         const ShakeList& shakelist,
                         PVA& posvelarray)
{
  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    getpos(posvelarray,it->first) = getpos(particlearray,it->first);
    getvelocity(posvelarray,it->first) = getvelocity(particlearray,it->first);
    int nh1 = it->second.nh1;
    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      getpos(posvelarray,it->second.h1[n1]) = getpos(particlearray,it->second.h1[n1]);
      getvelocity(posvelarray,it->second.h1[n1]) = getvelocity(particlearray,it->second.h1[n1]);
    }
  }
}

inline
void clear_shake_force(const ShakeList& shakelist,
		       ForceArray& shake_force)
{
  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    shake_force[it->first] = 0.0;
    int nh1 = it->second.nh1;
    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      shake_force[it->second.h1[n1]] = 0.0;
    }
  }
}

template<class PA, class PVA> inline
void calc_shake_force(const PA& particlearray,
		      const ShakeList& shakelist,
		      const PVA& posvelarray,
		      const double deltat,
		      ForceArray& shake_force)
{
  double inv_dt = 1.0/deltat;
  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    shake_force[it->first] += (getvelocity(particlearray,it->first) - getvelocity(posvelarray,it->first))*getmass(particlearray,it->first)*inv_dt;
    int nh1 = it->second.nh1;
    for(int n1=0; n1<nh1; n1++){
      // H atoms of shake bond
      shake_force[it->second.h1[n1]] += (getvelocity(particlearray,it->second.h1[n1]) - getvelocity(posvelarray,it->second.h1[n1]))*getmass(posvelarray,it->second.h1[n1])*inv_dt;
    }
  }
}

template<class PA>
void correct_translation(PA& particlearray, const std::vector<TypeRange>& typerangearray, const double inv_nt, const int node_id)
{
  Velocity v_sum_local(0.0,0.0,0.0);
  for(size_t s=0;s<typerangearray.size();s++){
    Velocity v_sum_c(0.0,0.0,0.0);
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      v_sum_c += getvelocity(particlearray,i);
    }
    v_sum_local += v_sum_c;
  }
  Velocity v_sum(0.0,0.0,0.0);
  reduce_velocity(v_sum,v_sum_local);
  v_sum *= inv_nt;

  if((DebugLog::verbose>0)&&(node_id==0)){
    printf("velocity translation %f %f %f\n",v_sum.x,v_sum.y,v_sum.z);
  }

  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getvelocity(particlearray,i) -= v_sum;
    }
  }
}

/*! Velocity Verlet top half for NVE/NVT
  @brief update velocity(t+deltat/2) by force,(t) position(t+deltat) by velocity(t+deltat/2)
  @param [in,out] particlearray vector of Particle
  @param [in] typerangearray vector of typerange
  @param [in] deltat delta t
 */
template<class IPA, class IGPA> inline
void Integrator<IPA, IGPA>::velocity_verlet_pos_vel(IPA& particlearray,
					 const std::vector<TypeRange>& typerangearray,
					 const double deltat)
{
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      double vf = 0.5*deltat*getinvmass(particlearray,i);
      getvelocity(particlearray,i) += vf*getforce(particlearray,i);
      getpos(particlearray,i) += deltat*getvelocity(particlearray,i);
    }
  }
}
/*! Velocity Verlet top half for NPT
  @brief update velocity(t+deltat/2) by force,(t) position(t+deltat) by velocity(t+deltat/2)
  @param [in,out] particlearray vector of Particle
  @param [in] typerangearray vector of typerange
  @param [in] deltat delta t
  @param [in,out] logv_pos log(volume)
  @param [in,out] boxsize size of box
  @param [in,out] volume volume
 */
template<class IPA, class IGPA> inline
void Integrator<IPA, IGPA>::velocity_verlet_npt_pos_vel(IPA& particlearray,
					     const std::vector<TypeRange>& typerangearray,
					     const double deltat)
{
#define E1 1.0
#define E2 (E1)/6.0
#define E4 (E2)/20.0
#define E6 (E4)/42.0
#define E8 (E6)/72.0 

  double scaling_Vol = exp(0.5*deltat*logv_vel);
  double arg2 = (0.5*deltat*logv_vel)*(0.5*deltat*logv_vel);
  double poly = (((E8*arg2+E6)*arg2+E4)*arg2+E2)*arg2+E1;
  double scaling_VV = deltat*scaling_Vol*poly;
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      double vf = 0.5*deltat*getinvmass(particlearray,i);
      getvelocity(particlearray,i) += vf*getforce(particlearray,i);
      getpos(particlearray,i) *= (scaling_Vol*scaling_Vol);
      getpos(particlearray,i) += getvelocity(particlearray,i)*scaling_VV;
    }
  }
  logv_pos += deltat*logv_vel;
  boxsize *= (scaling_Vol*scaling_Vol);
  if(boxsize.x*inv_initial_boxsize.x<boxsize_minimum_scale){
    printf("Warning boxsize.x %f < initial_boxsize.x %f * %f\n",boxsize.x,inv_initial_boxsize.x,boxsize_minimum_scale);
  }
  volume = boxsize.x*boxsize.y*boxsize.z;
}

template<class IPA, class IGPA> inline
void Integrator<IPA, IGPA>::velocity_verlet_npt_pos_vel_vtcorr(IPA& particlearray,
							       const std::vector<TypeRange>& typerangearray,
							       const double deltat, const double inv_nt)
{

  double scaling_Vol = exp(0.5*deltat*logv_vel);
  double arg2 = (0.5*deltat*logv_vel)*(0.5*deltat*logv_vel);
  double poly = (((E8*arg2+E6)*arg2+E4)*arg2+E2)*arg2+E1;
  double scaling_VV = deltat*scaling_Vol*poly;
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      double vf = 0.5*deltat*getinvmass(particlearray,i);
      getvelocity(particlearray,i) += vf*getforce(particlearray,i);
    }
  }
  correct_translation(particlearray,typerangearray,inv_nt,unit_identifier);
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getpos(particlearray,i) *= (scaling_Vol*scaling_Vol);
      getpos(particlearray,i) += getvelocity(particlearray,i)*scaling_VV;
    }
  }
  logv_pos += deltat*logv_vel;
  boxsize *= (scaling_Vol*scaling_Vol);
  if(boxsize.x*inv_initial_boxsize.x<boxsize_minimum_scale){
    printf("Warning boxsize.x %f < initial_boxsize.x %f * %f\n",boxsize.x,inv_initial_boxsize.x,boxsize_minimum_scale);
  }
  volume = boxsize.x*boxsize.y*boxsize.z;
}


template<class IPA, class IGPA> inline
void Integrator<IPA, IGPA>::velocity_verlet_velocity(IPA& particlearray,
                              const std::vector<TypeRange>& typerangearray,
                              const double deltat)
{
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      double vf = 0.5*deltat*getinvmass(particlearray,i);
      getvelocity(particlearray,i) += vf*getforce(particlearray,i);
    }
  }
}

template<class PA> inline
void velocity_scaling(PA& particlearray,
                      const std::vector<TypeRange>& typerangearray,
                      double& ke_sum, double& temperature,
                      const double ref_temperature)
{
  double scaling_factor2 = ref_temperature/temperature;
  double scaling_factor = sqrt(scaling_factor2);
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getvelocity(particlearray,i) *= scaling_factor;
    }
  }
  ke_sum *= scaling_factor2;
  temperature *= scaling_factor2;
}

/*! Nose Hoover top half
  @param [in,out] particlearray vector of Particle
  @param [in] typerangearray vector of typerange
  @param [in,out] eta_pos position of temperature control
  @param [in,out] eta_vel velocity of temperature control
  @param [in,out] eta_force forece of temperature control
  @param [in,out] ke_sum total kinetic energy
  @param [in] eta_inv_mass mass of temperature control
  @param [in] ref_temperature reference temperatur
  @param [in] inv_unit_energy inverse of unit energy
  @param [in] num_freedom number of freedom
  @param [in] deltat delta t
*/
template<class PA> inline
void nose_hoover_pre(PA& particlearray,
                     const std::vector<TypeRange>& typerangearray,
                     double& eta_pos, double& eta_vel, double& eta_force,
                     double& ke_sum,
                     const double eta_inv_mass,
                     const double ref_temperature,
                     const double inv_unit_energy,
                     const int num_freedom,
                     const double deltat)
{
  eta_vel += 0.25 * deltat * eta_inv_mass * eta_force;
  eta_pos += 0.5 * deltat * eta_vel;

  double scaling_NH = exp(-0.5 * deltat * eta_vel);
  ke_sum *= (scaling_NH * scaling_NH);

  eta_force = (2.0 * ke_sum) - (num_freedom * inv_unit_energy * UnitParameter::Boltzmann * ref_temperature); 
  eta_vel += 0.25 * deltat * eta_inv_mass * eta_force;
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getvelocity(particlearray,i) *= scaling_NH;
    }
  }
}

template<class PA> inline
void nose_hoover_post(PA& particlearray,
                      const std::vector<TypeRange>& typerangearray,
                      double& eta_pos, double& eta_vel, double& eta_force,
                      double& ke_sum, double& temperature,
                      const double eta_inv_mass,
                      const double ref_temperature,
                      const double inv_unit_energy,
                      const int num_freedom,
                      const double deltat)
{
  eta_force = (2.0 * ke_sum) - (num_freedom * inv_unit_energy * UnitParameter::Boltzmann * ref_temperature); 

  eta_vel += 0.25 * deltat * eta_inv_mass * eta_force;
  eta_pos += 0.5 * deltat * eta_vel;

  double scaling_NH = exp(-0.5 * deltat * eta_vel);
  ke_sum *= (scaling_NH * scaling_NH);

  eta_force = (2.0 * ke_sum) - (num_freedom * inv_unit_energy * UnitParameter::Boltzmann * ref_temperature); 
  eta_vel += 0.25 * deltat * eta_inv_mass * eta_force;
  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getvelocity(particlearray,i) *= scaling_NH;
    }
  }
  temperature *= (scaling_NH)*(scaling_NH);
}

/*! Andersen Hoover top half
  @param [in,out] particlearray vector of Particle
  @param [in] typerangearray vector of typerange
  @param [in,out] eta_pos position of temperature control
  @param [in,out] eta_vel velocity of temperature control
  @param [in,out] eta_force forece of temperature control
  @param [in,out] logv_vel velocity of log(volume)
  @param [in,out] logv_force force of log(volume)
  @param [in,out] ke_sum total kinetic energy
  @param [in] eta_inv_mass mass of temperature control
  @param [in] logv_mass mass of log(volume)
  @param [in] logv_mass inverse of logv_mass
  @param [in] ref_temperature reference temperatur
  @param [in] inv_unit_energy inverse of unit energy
  @param [in] num_freedom number of freedom
  @param [in] diff_pressure_volume (difference of pressure)*volume
  @param [in] deltat delta t
*/
template<class PA> inline
void andersen_hoover_pre(PA& particlearray,
                          const std::vector<TypeRange>& typerangearray,
                          double& eta_pos, double& eta_vel, double& eta_force,
                          double& logv_vel, double& logv_force,
                          double& ke_sum,
                          const double eta_inv_mass,
                          const double logv_mass, const double logv_inv_mass,
                          const double ref_temperature,
                          const double inv_unit_energy,
                          const int num_freedom,
                          const double diff_pressure_volume,
                          const double deltat)
{
  double dinvnf = 3.0/num_freedom;
  double onedinvnf = 1.0 + dinvnf;

  eta_vel += eta_force*eta_inv_mass*0.25*deltat;
  eta_pos += eta_vel*0.5*deltat;

  double scaling_VV = exp(-0.125*deltat*eta_vel);
  logv_vel = logv_vel*scaling_VV*scaling_VV + 0.25*deltat*logv_force*logv_inv_mass*scaling_VV;

#ifdef MODIFIED_ANDERSEN
  double scaling_AH = exp(-0.5*deltat*(eta_vel + onedinvnf*logv_vel));
#else
  double scaling_AH = exp(-0.5*deltat*(eta_vel + logv_vel));
#endif
  ke_sum *= scaling_AH*scaling_AH;
#ifdef MODIFIED_ANDERSEN
  //  logv_force = (onedinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// ????
  logv_force = (dinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// p4178 eq(2.9)
#else
  logv_force = 3.0*diff_pressure_volume;  /// p4178 eq(2.7)
#endif

  scaling_VV = exp(-0.125*deltat*eta_vel);
  logv_vel = logv_vel*scaling_VV*scaling_VV + 0.25*deltat*logv_force*logv_inv_mass*scaling_VV;

  eta_force = 2.0*ke_sum + logv_mass*logv_vel*logv_vel
              - (1+num_freedom)*inv_unit_energy*UnitParameter::Boltzmann*ref_temperature;

  eta_vel += eta_force*eta_inv_mass*0.25*deltat;

  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getvelocity(particlearray,i) *= scaling_AH;
    }
  }
}

template<class PA> inline
void andersen_hoover_post(PA& particlearray,
                          const std::vector<TypeRange>& typerangearray,
                          double& eta_pos, double& eta_vel, double& eta_force,
                          double& logv_vel, double& logv_force,
                          double& ke_sum, double& temperature,
                          const double eta_inv_mass,
                          const double logv_mass, const double logv_inv_mass,
                          const double ref_temperature,
                          const double inv_unit_energy,
                          const int num_freedom,
                          const double diff_pressure_volume,
                          const double deltat)
{
  double dinvnf = 3.0/num_freedom;
  double onedinvnf = 1.0 + dinvnf;

  eta_force = 2.0*ke_sum + logv_mass*logv_vel*logv_vel
               - (1+num_freedom)*inv_unit_energy*UnitParameter::Boltzmann*ref_temperature;
#ifdef MODIFIED_ANDERSEN
  //  logv_force = (onedinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// ???
  logv_force = (dinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// p4178 eq(2.9)
#else
  logv_force = (3.0*diff_pressure_volume);  /// p4178 eq(2.7)
#endif
  
  eta_vel += eta_force*eta_inv_mass*0.25*deltat;
  eta_pos += eta_vel*0.5*deltat;

  double scaling_VV = exp(-0.125*deltat*eta_vel);
  logv_vel = logv_vel*scaling_VV*scaling_VV + 0.25*deltat*logv_force*logv_inv_mass*scaling_VV;

#ifdef MODIFIED_ANDERSEN
  double scaling_AH = exp(-0.5*deltat*(eta_vel + onedinvnf*logv_vel));
#else
  double scaling_AH = exp(-0.5*deltat*(eta_vel + logv_vel));
#endif
  ke_sum *= scaling_AH*scaling_AH;
#ifdef MODIFIED_ANDERSEN
  //  logv_force = (onedinvnf*2.0*ke_sum + 3.0*diff_pressure_volume); /// ???
  logv_force = (dinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// p4178 eq(2.9)
#else
  logv_force = (3.0*diff_pressure_volume);  /// p4178 eq(2.7)
#endif

  scaling_VV = exp(-0.125*deltat*eta_vel);
  logv_vel = logv_vel*scaling_VV*scaling_VV + 0.25*deltat*logv_force*logv_inv_mass*scaling_VV;

  eta_force = 2.0*ke_sum + logv_mass*logv_vel*logv_vel
              - (1+num_freedom)*inv_unit_energy*UnitParameter::Boltzmann*ref_temperature;

  eta_vel += eta_force*eta_inv_mass*0.25*deltat;

  for(size_t s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      getvelocity(particlearray,i) *= scaling_AH;
    }
  }

  temperature *= (scaling_AH)*(scaling_AH);
}

double andersen_hoover(double& eta_pos, double& eta_vel, double& eta_force,
		       double& logv_vel, double& logv_force,
		       double& ke_sum, double& temperature,
		       const double eta_inv_mass,
		       const double logv_mass, const double logv_inv_mass,
		       const double ref_temperature,
		       const double inv_unit_energy,
		       const int num_freedom,
		       const double diff_pressure_volume,
		       const double deltat)
{
  double dinvnf = 3.0/num_freedom;
  double onedinvnf = 1.0 + dinvnf;

  eta_force = 2.0*ke_sum + logv_mass*logv_vel*logv_vel
               - (1+num_freedom)*inv_unit_energy*UnitParameter::Boltzmann*ref_temperature;
#ifdef MODIFIED_ANDERSEN
  //  logv_force = (onedinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// ??
  logv_force = (dinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// p4178 eq(2.9)
#else
  logv_force = (3.0*diff_pressure_volume);  /// p4178 eq(2.7)
#endif
  
  eta_vel += eta_force*eta_inv_mass*0.25*deltat;
  eta_pos += eta_vel*0.5*deltat;

  double scaling_VV = exp(-0.125*deltat*eta_vel);
  logv_vel = logv_vel*scaling_VV*scaling_VV + 0.25*deltat*logv_force*logv_inv_mass*scaling_VV;

#ifdef MODIFIED_ANDERSEN
  double scaling_AH = exp(-0.5*deltat*(eta_vel + onedinvnf*logv_vel));
#else
  double scaling_AH = exp(-0.5*deltat*(eta_vel + logv_vel));
#endif
  ke_sum *= scaling_AH*scaling_AH;
#ifdef MODIFIED_ANDERSEN
  //  logv_force = (onedinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// ?
  logv_force = (dinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// p4178 eq(2.9)
#else
  logv_force = (3.0*diff_pressure_volume);  /// p4178 eq(2.7)
#endif

  scaling_VV = exp(-0.125*deltat*eta_vel);
  logv_vel = logv_vel*scaling_VV*scaling_VV + 0.25*deltat*logv_force*logv_inv_mass*scaling_VV;

  eta_force = 2.0*ke_sum + logv_mass*logv_vel*logv_vel
              - (1+num_freedom)*inv_unit_energy*UnitParameter::Boltzmann*ref_temperature;

  eta_vel += eta_force*eta_inv_mass*0.25*deltat;
  
  return scaling_AH;
}

template <class PA> inline
void estimate_1st_water_settle_force(PA& particlearray,
				     const std::vector<TypeRange>& typerangearray,
				     const WaterList& waterlist,
				     const double deltat,
				     Settle& settle,
				     PosVelArray& work_posvel,
				     ForceArray& water_settle_force)
{
  // backup and update velocity by half step;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++){
    double vf = 0.5*deltat*getinvmass(particlearray,it->first);
    double vf1 = 0.5*deltat*getinvmass(particlearray,it->second.h1);
    double vf2 = 0.5*deltat*getinvmass(particlearray,it->second.h2);
    getvelocity(work_posvel,it->first) = getvelocity(particlearray,it->first);
    getvelocity(work_posvel,it->second.h1) = getvelocity(particlearray,it->second.h1);
    getvelocity(work_posvel,it->second.h2) = getvelocity(particlearray,it->second.h2);
    getvelocity(particlearray,it->first) += vf*getforce(particlearray,it->first);
    getvelocity(particlearray,it->second.h1) += vf1*getforce(particlearray,it->second.h1);
    getvelocity(particlearray,it->second.h2) += vf2*getforce(particlearray,it->second.h2);
  }
  // calculate settle and calculate settle force
  settle.constrain_velocity(particlearray, waterlist, typerangearray[0],deltat,water_settle_force);
  // restore velocity
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++){
    getvelocity(particlearray,it->first) = getvelocity(work_posvel,it->first);
    getvelocity(particlearray,it->second.h1) = getvelocity(work_posvel,it->second.h1);
    getvelocity(particlearray,it->second.h2) = getvelocity(work_posvel,it->second.h2);
  }
}

template <class PA> inline
void estimate_1st_water_settle_virial(PA& particlearray,
				     const std::vector<TypeRange>& typerangearray,
				     const WaterList& waterlist,
				     const double deltat,
				     Settle& settle,
				     PosVelArray& work_posvel,
				      double& water_settle_virial)
{
  // backup and update velocity by half step;
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++){
    double vf = 0.5*deltat*getinvmass(particlearray,it->first);
    double vf1 = 0.5*deltat*getinvmass(particlearray,it->second.h1);
    double vf2 = 0.5*deltat*getinvmass(particlearray,it->second.h2);
    getvelocity(work_posvel,it->first) = getvelocity(particlearray,it->first);
    getvelocity(work_posvel,it->second.h1) = getvelocity(particlearray,it->second.h1);
    getvelocity(work_posvel,it->second.h2) = getvelocity(particlearray,it->second.h2);
    getvelocity(particlearray,it->first) += vf*getforce(particlearray,it->first);
    getvelocity(particlearray,it->second.h1) += vf1*getforce(particlearray,it->second.h1);
    getvelocity(particlearray,it->second.h2) += vf2*getforce(particlearray,it->second.h2);
  }
  // calculate settle and calculate settle force
  settle.constrain_velocity(particlearray, waterlist, typerangearray[0],deltat,water_settle_virial);
  // restore velocity
  for(WaterList::const_iterator it=waterlist.begin(); it != waterlist.end(); it++){
    getvelocity(particlearray,it->first) = getvelocity(work_posvel,it->first);
    getvelocity(particlearray,it->second.h1) = getvelocity(work_posvel,it->second.h1);
    getvelocity(particlearray,it->second.h2) = getvelocity(work_posvel,it->second.h2);
  }
}

template <class PA> inline
void estimate_1st_shake_force(PA& particlearray,
			      const ShakeList& shakelist,
			      const CovalentBondParameterList* param_list, 
			      std::vector<CovalentBondInfo::BondList>& bondlistarray,
			      const double deltat,
			      int shake_max_iterate, double shake_tolerance,
			      Settle& settle,
			      PosVelArray& work_posvel,
			      ForceArray& shake_force)
{
  // backup and update velocity by half step;
  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    double vf = 0.5*deltat*getinvmass(particlearray,it->first);
    double vf1[MAX_HBOND];
    int nh1 = it->second.nh1;
    for(int n1=0; n1<nh1; n1++){
      vf1[n1] = 0.5*deltat*getinvmass(particlearray,it->second.h1[n1]);
    }
    getvelocity(work_posvel,it->first) = getvelocity(particlearray,it->first);
    for(int n1=0; n1<nh1; n1++){
      getvelocity(work_posvel,it->second.h1[n1]) = getvelocity(particlearray,it->second.h1[n1]);
    }
    getvelocity(particlearray,it->first) += vf*getforce(particlearray,it->first);
    for(int n1=0; n1<nh1; n1++){
      getvelocity(particlearray,it->second.h1[n1]) += vf1[n1]*getforce(particlearray,it->second.h1[n1]);
    }
  }
  // calculate settle and calculate settle force
  settle.rattle_velocity(particlearray, shakelist, param_list, bondlistarray,
			 deltat, shake_max_iterate, shake_tolerance,
			 shake_force);
  // restore velocity
  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int nh1 = it->second.nh1;
    getvelocity(particlearray,it->first) = getvelocity(work_posvel,it->first);
    for(int n1=0; n1<nh1; n1++){
      getvelocity(particlearray,it->second.h1[n1]) = getvelocity(work_posvel,it->second.h1[n1]);
    }
  }
}
template <class PA> inline
void estimate_1st_shake_virial(PA& particlearray,
			      const ShakeList& shakelist,
			      const CovalentBondParameterList* param_list, 
			      std::vector<CovalentBondInfo::BondList>& bondlistarray,
			      const double deltat,
			      int shake_max_iterate, double shake_tolerance,
			      Settle& settle,
			      PosVelArray& work_posvel,
			      double& shake_virial)
{
  // backup and update velocity by half step;
  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    double vf = 0.5*deltat*getinvmass(particlearray,it->first);
    double vf1[MAX_HBOND];
    int nh1 = it->second.nh1;
    for(int n1=0; n1<nh1; n1++){
      vf1[n1] = 0.5*deltat*getinvmass(particlearray,it->second.h1[n1]);
    }
    getvelocity(work_posvel,it->first) = getvelocity(particlearray,it->first);
    for(int n1=0; n1<nh1; n1++){
      getvelocity(work_posvel,it->second.h1[n1]) = getvelocity(particlearray,it->second.h1[n1]);
    }
    getvelocity(particlearray,it->first) += vf*getforce(particlearray,it->first);
    for(int n1=0; n1<nh1; n1++){
      getvelocity(particlearray,it->second.h1[n1]) += vf1[n1]*getforce(particlearray,it->second.h1[n1]);
    }
  }
  // calculate settle and calculate settle force
  settle.rattle_velocity(particlearray, shakelist, param_list, bondlistarray,
			 deltat, shake_max_iterate, shake_tolerance,
			 shake_virial);
  // restore velocity
  for(ShakeList::const_iterator it=shakelist.begin(); it != shakelist.end(); it++){
    int nh1 = it->second.nh1;
    getvelocity(particlearray,it->first) = getvelocity(work_posvel,it->first);
    for(int n1=0; n1<nh1; n1++){
      getvelocity(particlearray,it->second.h1[n1]) = getvelocity(work_posvel,it->second.h1[n1]);
    }
  }
}

inline
void dump_andersen_hoover(const double ke_sum, const double virial_sum, 
			  const double num_freedom, const double inv_unitEnergy_eV,
			  const double temperature, const double volume,
			  const double eta_pos, const double eta_vel, const double eta_force,
			  const double logv_pos, const double logv_vel, const double logv_force,
			  const double pressure, const SpaceVector<double> boxsize)
{
  printf("ke_sum %f, virial_sum %f,  nkT/volume %f\n", ke_sum, virial_sum,num_freedom*inv_unitEnergy_eV*UnitParameter::Boltzmann*temperature/volume);
  printf("eta_pos %f, eta_vel %f, eta_force %f\n",eta_pos, eta_vel, eta_force);
  printf("logv_pos %f, logv_vel %f, logv_force %f, pressure %f %e Pa %e atm\n",logv_pos, logv_vel, logv_force,pressure,pressure*UnitParameter::unitEnergy_J/(UnitParameter::unitLength*UnitParameter::unitLength*UnitParameter::unitLength),pressure*UnitParameter::unitEnergy_J/(UnitParameter::unitLength*UnitParameter::unitLength*UnitParameter::unitLength)/101325.0);
  std::cout << "boxsize " << boxsize << " volume " << volume << " log(volume)/3.0 "<< log(volume)/3.0 << std::endl;
}

/*!
  @todo correct move timing when restart unaligned time step
  @todo calculate force before main loop except restart
 */
template<class IPA, class IGPA>
bool Integrator<IPA,IGPA>::time_integrations(//ParticleArray& init_particlearray, 
				   IPA& particlearray,
                                   std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
                                   WaterList& waterlist,
                                   ShakeList& shakelist,
                                   std::vector<int>& setid,
                                   std::map<int,int>& setid_to_index,
                                   CalcForce& calcforce, 
                                   double dt, const long tmax, 
                                   OperationSelector operations,
                                   MPICommunicator& communicator,
                                   MPIReceiverForLong& longreceiver,
                                   MPISenderForLong& longsender,
                                   PostProcess& postprocess,
                                   const ShakeList& sl,
                                   int shake_type, int shake_max_iterate, double shake_tolerance,
                                   int reduce_interval, int printinterval,
                                   int dump_crd_inter, int dump_rst_inter, Dump& dump, Dump& dump_rst,
                                   int moveinter,
                                   int num_total,
                                   int num_freedom,
                                   int& write_restart,
				   long& last_t){

#ifdef FPCOLL
#ifdef FPCOLL_INTEGRATOR
  fpcoll_start();
#endif
#endif
  SAMPLE_START();

  timer_start(timer1);

  bool overlap=false;

  IGPA ghost;

  PosVelArray pv_array_prev(particlearray.size());
  ForceArray water_settle_force;

  // Roll Back move particle
  backup_time=0;
  //  int backup_moveinter;
  int max_moveinter;
  max_moveinter = moveinter;
  //  int min_moveinter=5;
  bool rollback=false;
  int over_move;
  bool not_over_flow;
  //  int over_move_all;
  backup_pi=0;
  backup_ri=0;
  backup_dci=0;
  backup_dri=0;

  //  ParticleArray particle_array_prev(particlearray.size());
  //  ParticleArray particle_array_prev_shake(particlearray.size());
  PosVelArray particle_array_prev_shake(particlearray.size());
  ForceArray shake_force;
  double inv_nt = 1.0/num_total;
  double inv_num_freedom = 1.0/num_freedom;
  double inv_boltzmann = 1.0/UnitParameter::Boltzmann;
  double inv_unitEnergy_eV = 1.0/UnitParameter::unitEnergy_eV;
  double ke_sum=0.0;
  double energy_sum=0.0;
  double virial_sum=0.0;
  double temperature=0.0;
  Force force_sum(0.0,0.0,0.0);
#ifdef CHECK_ENERGY
  double lje_sum, se_sum, zde_sum, sse_sum;
  double le_sum, exe_sum;
  double clj_sum, ccl_sum;
#endif
  communicator.resize_receive_data(targetparticlerange);
  volume = boxsize.x*boxsize.y*boxsize.z;
  tljcec = calcTotalLJCutoffEnergyCorrection(volume);
  // TODO fix size of targetbond
  // targetbond.resize(1);
  force.resize(particlearray.size()*2);
  ghostforce.resize(particlearray.size()*2,Force(0.0,0.0,0.0));
  calcforce.bondlist.clear();
  calcforce.bondlist.merge_bondlistarray(bondlistarray);
  pi=0;
  ri=0;
  dci=0;
  dri=0;
  int mi=moveinter-1;
  int moved=1;
  long t=0;
  double ke=0.0;
  double lv=0.0;
  double pairwise_virial=0.0;
  double water_settle_virial=0.0;
  double shake_virial=0.0;
  int mt_bond = mti_bond;
  int mt_short = mti_short;
  int mt_long = mti_long;
  OperationSelector mt_operations;

  double time_stamp, last_time_stamp, lap_time;
  int time_stamp_interval=1000;

#ifdef MT_LONG
  if(MT_LONG>0){
    mti_long = MT_LONG-1;
  }
  if(unit_identifier==0){
    std::cout << "Long Range calculation per " << mti_long+1 << " step " << std::endl;
  }
#endif

  double nxt_dt = dt*reduce_interval;
  if(reduce_interval<0)nxt_dt = dt;

  // constant temerature, pressure
  ref_temperature = temp_control.temperature;
  // nose-hoover
  tau = temp_control.tau_t * 1000.0 *  1e-15 / UnitParameter::unitTime;
  eta_inv_mass = 1.0/(num_freedom * inv_unitEnergy_eV * UnitParameter::Boltzmann * ref_temperature * tau * tau);
  eta_pos = 0.0;
  eta_vel = 0.0;
  eta_force = 0.0;
  // andersen-hoover NPT  MARK E. TUCKERMAN, et.al Molecular Physics, 1996, Vol.87, No.5, 1117-1157
  ref_pressure=101325.0/(UnitParameter::unitEnergy_J)*UnitParameter::unitLength*UnitParameter::unitLength*UnitParameter::unitLength;  //! TODO given by option 
  tauv = temp_control.tau_p * 1000.0 *  1e-15 / UnitParameter::unitTime;
  logv_mass = (num_freedom + 3) * inv_unitEnergy_eV * UnitParameter::Boltzmann * ref_temperature * tauv * tauv;   //! 3 is dimension  p1121 2.5    W=(Nf+d)kT/omega^2
  logv_inv_mass = 1.0/logv_mass;
  logv_pos = log(volume)/3.0;  //! 3.0 is dimension
  logv_vel=0.0;
  logv_force=0.0;

  reset_ghost_longset_index = true;

  std::vector<int> ghostcbindexarray;
  //  CBModule::CBUnordered<AtomID, CBModule::CBInterface_::ForceLocation>::Map ghostcbforcemap;
  unsigned int num_ghost=0;

  bool shiftidcheck = true;

  bool result = true;

  bool restored = false;

  double stored_total_penergy=0.0;

  postprocess.change_fine_margin(moveinter);

#ifdef USE_HDF
  if(restore){
    if(operations.doShortrangecalculation){
      time_stamp = getrealtime();
      if(unit_identifier==0)printf("restore HDF start at %f sec\n",time_stamp);
      last_time_stamp = time_stamp;
      int comm_size;
      MPI_Comm_size(communicator.mpi_comm_short, &comm_size);
      restore_hdf_snapshot(particlearray,
			  typerangearray,bondlistarray,bondlistarray_idx,waterlist,shakelist,ke,energy,lv,stored_total_penergy,t, num_total, comm_size);
      //    printf("energy %d: %f\n",unit_identifier,energy);
      close_hdfrestore();
      time_stamp = getrealtime();
      if(unit_identifier==0)printf("restore HDF end at %f sec, spent %f sec\n",time_stamp-last_time_stamp);
      last_time_stamp = time_stamp;
    }
    restored = true;
    if(unit_identifier==0){
      printf("Restore from HDF dump\n");
      printf("ri pi dci dri : %d %d %d %d\n",ri,pi,dci,dri);
      //      printf("Restored tljcec %e\n",tljcec);
    }
    pre_calcforce = false;
  }
#else
#ifdef BINARY_DUMP
  if(restore){
    if(operations.doShortrangecalculation){
      time_stamp = getrealtime();
      if(unit_identifier==0)printf("restore binary start at %f sec\n",time_stamp);
      last_time_stamp = time_stamp;
      restore_binary_snapshot(particlearray,
			      typerangearray,bondlistarray,bondlistarray_idx,waterlist,shakelist,ke,energy,lv,stored_total_penergy,t);
      //    printf("energy %d: %f\n",unit_identifier,energy);
      close_binaryrestore();
      time_stamp = getrealtime();
      if(unit_identifier==0)printf("restore binary end at %f sec, spent %f sec\n",time_stamp-last_time_stamp);
      last_time_stamp = time_stamp;
    }
    restored = true;
    if(unit_identifier==0){
      printf("Restore from binary dump\n");
      printf("ri pi dci dri : %d %d %d %d\n",ri,pi,dci,dri);
      //      printf("Restored tljcec %e\n",tljcec);
    }
# ifndef DUMP_WITHOUT_FORCE
    pre_calcforce = false;
# endif

  }
#endif
#endif

  calcforce.shortforce.resize(particlearray.size());

  //  operations.doVirialcalculation = true;
  if(temp_control.method == Config::kANDERSEN_HOOVER){
    operations.doVirialcalculation = true;
    water_settle_force.resize(particlearray.size());
    if(shake_type > 0){
      shake_force.resize(particlearray.size());
    }
  }
  
  /*
      {
	int np=0;
	for(int t=0;t<typerangearray.size();t++){
	  np += (typerangearray[t].end-typerangearray[t].begin);
	}
	printf("short_id %d num particle %d\n",unit_identifier,np);
      }
  */
  if(unit_identifier==0){
    printf("Unit of Energy = %e eV = %e J = %e kcal/mol\n",
	   UnitParameter::unitEnergy_eV,UnitParameter::unitEnergy_J,UnitParameter::unitEnergy_kcal_mol);
  }

  //! pre calculation force
  //  if(pre_calcforce&&(tmax>0)){
#ifdef PROF_CHECK
# ifdef PROF_CHECK_VERB
  if(unit_identifier==0){
    printf("MPI_Barrier before pre_calcforce\n");
  }
  fflush(stdout);
# endif
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(pre_calcforce){
    if(unit_identifier==0){
      std::cout << " Calc Force pre integration " << std::endl;
    }
    if(restored){  /// move not required binary restored case
      rollback=false;
    }else
    {
      timer_start(timer_move);
      //      std::cout << "move_top_half  " << unit_identifier  << std::endl;
      if(operations.doShortrangecalculation==true){
	not_over_flow = move_top_half(particlearray, typerangearray, bondlistarray, waterlist, shakelist, communicator, postprocess, over_move);
      }else{
	not_over_flow = true;
	over_move = 0;
      }
      //      std::cout << "check_move_over  " << unit_identifier  << std::endl;
      rollback = check_move_over(unit_identifier, not_over_flow, over_move);
      if(rollback==false){
	//	std::cout << "move_bottom_half  " << unit_identifier  << std::endl;
	if(operations.doShortrangecalculation==true){
	  not_over_flow = move_bottom_half(particlearray, typerangearray, bondlistarray, waterlist, shakelist, communicator, postprocess,sl);
	}else{
	  not_over_flow = true;
	}
	//	std::cout << "check_move_over  " << unit_identifier  << std::endl;
	rollback = check_move_over(unit_identifier, not_over_flow, over_move);
      }
#ifdef PROF_CHECK
# ifdef PROF_CHECK_VERB
    if(unit_identifier==0){
      printf("MPI_Barrier after move_bottom_half\n");
    }
    fflush(stdout);
# endif
    MPI_Barrier(MPI_COMM_WORLD);
# ifdef PROF_CHECK_VERB
    if(unit_identifier==0){
      printf("MPI_Barrier after move_bottom_half done\n");
    }
    fflush(stdout);
# endif
#endif
      timer_stop(timer_move);
      if(DebugLog::verbose>1){
        dump_bondlistarray(bondlistarray);
      }
    }
    if(rollback==false){
      {
	clear_before_move(ghost, targettyperange, recvsetid, bondlistarray, calcforce);
	if(DebugLog::verbose>1){
	  dump_atomid(particlearray,typerangearray,setid);
	}
      }
      clear_force(particlearray,typerangearray);
      clear_force(force,typerangearray);

      timer_start(timer_exp);
      calcforce.covalentbond.set_cancel_short();
      if(operations.doVirialcalculation){
	calcforce.covalentbond.set_cancel_short_virial();
      }
      calcforce.selfshort_zero();
      calcforce.long_zero();
#ifdef CHECK_ENERGY
      {
	ljpair_energy = 0.0;
	shortpair_energy = 0.0;
	zerodipole_energy = 0.0;
	ShortRange::short_coulomb = 0.0;
	ShortRange::short_lj = 0.0;
      }
#endif
      pairwise_virial = 0.0;
#ifdef USE_PAIRLIST
      operations.doMakePairList = true;
#endif
      mt_operations = operations;
#ifdef TEST_CB_ONLY
      mt_operations.doShortrangecalculation = false;
      calcforce.covalentbond.unset_cancel_short();
      if(mt_operations.doVirialcalculation){
	calcforce.covalentbond.unset_cancel_short_virial();
      }
#endif
#ifdef CPPMD_ENABLE_PMMM
      calcforce.calcPMMM_top_half(particlearray, typerangearray,operations);
#endif
#ifdef OVERLAP
      exchangeParticleArray_top_half(particlearray, typerangearray, bondlistarray, setid,
				     ghost, targetparticlerange, 
				     recvsetid, recvsetid_to_index, 
				     targettyperange, targetbond, ghostcbindexarray, num_ghost,
				     calcforce, communicator, longreceiver, longsender,
				     reset_ghost_longset_index,
				     operations);
      timer_stop(timer_exp);
      timer_start(timer_expcalc);
      calcforce.calcForce_local(particlearray, typerangearray, 
			  bondlistarray,bondlistarray_idx, waterlist,
			  ghost, recvsetid, 
			  recvsetid_to_index, targettyperange, targetbond,
			  force, ghostforce, energy, pairwise_virial, mt_operations);
      timer_stop(timer_expcalc);
      timer_start(timer_exp);
      exchangeParticleArray_bottom_half(particlearray, typerangearray, bondlistarray, setid,
					ghost, targetparticlerange, 
					recvsetid, recvsetid_to_index, 
					targettyperange, targetbond, ghostcbindexarray, num_ghost,
					calcforce, communicator, longreceiver, longsender,
					reset_ghost_longset_index,
					operations);
#else
      exchangeParticleArray(particlearray, typerangearray, bondlistarray, setid,
			    ghost, targetparticlerange, 
			    recvsetid, recvsetid_to_index, 
			    targettyperange, targetbond, ghostcbindexarray, num_ghost,
			    calcforce, communicator, longreceiver, longsender,
			    reset_ghost_longset_index,
			    operations);
#endif
#ifdef PROF_CHECK
# ifdef PROF_CHECK_VERB
      if(unit_identifier==0){
	printf("MPI_Barrier after exchangeParticleArray\n");
      }
  fflush(stdout);
# endif
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      timer_stop(timer_exp);

      postprocess_after_exchangeparticle(postprocess,
					 ghost,recvsetid,recvsetid_to_index,targettyperange,
					 ghostforce,energy, shiftidcheck);
      calcforce.ghostshort_zero();

      timer_start(timer_calc);
#ifdef OVERLAP
      calcforce.calcForce(particlearray, typerangearray, 
			  bondlistarray,bondlistarray_idx, waterlist,
			  ghost, recvsetid, 
			  recvsetid_to_index, targettyperange, targetbond,
			  force, ghostforce, energy, pairwise_virial, mt_operations, false);
#else
      calcforce.calcForce(particlearray, typerangearray, 
			  bondlistarray,bondlistarray_idx, waterlist,
			  ghost, recvsetid, 
			  recvsetid_to_index, targettyperange, targetbond,
			  force, ghostforce, energy, pairwise_virial, mt_operations, !overlap);
#endif
#ifdef PROF_CHECK
# ifdef PROF_CHECK_VERB
      if(unit_identifier==0){
	printf("MPI_Barrier after calcForce %d\n",unit_identifier);
      }
  fflush(stdout);
# endif
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      timer_stop(timer_calc);

      if(operations.doExchangeForce){

	timer_start(timer_exf);
	exchangeForceArray(ghostforce, targetparticlerange, force,
			   calcforce.ghost_longset_index,
			   targettyperange, ghostcbindexarray,
			   communicator, longreceiver, longsender, operations, moved);
#ifdef PROF_CHECK
# ifdef PROF_CHECK_VERB
      if(unit_identifier==0){
	printf("MPI_Barrier after exchangeForceArray %d\n",unit_identifier);
      }
  fflush(stdout);
# endif
      MPI_Barrier(MPI_COMM_WORLD);
#endif
        timer_stop(timer_exf);

	mergerForces(particlearray);
      }
    }else{  /// move over, cell over flow
      last_t = 0;
      if(unit_identifier==0){
	printf("move over or cell over flow at initial data\n");
      }
      return false;
    }
    {
      double total_ss,tmp;
      reduce(total_ss,tmp,calcforce.shortenergyself,0.0);
      if(unit_identifier==0){
	printf("Short-range self energy %e\n",total_ss);
      }
    }
  }

#ifdef CORRECT_LEAPFROG_RESTART
  if(unit_identifier==0){
    printf("Correct velocity from restart file by leap-frog\n");
  }
  velocity_verlet_velocity(particlearray,typerangearray,-dt);
#endif

  if(restored){
  }else{
    timer_start(timer_settle);
    // settle for water
    if(operations.doVirialcalculation){
      estimate_1st_water_settle_virial(particlearray,typerangearray,waterlist,dt,settle,
				       pv_array_prev,water_settle_virial);
#ifdef USE_SHAKE
      if(shake_type>0){
	estimate_1st_shake_virial(particlearray, shakelist,
				  calcforce.covalentbond.getParameterList(), bondlistarray,
				  dt, shake_max_iterate, shake_tolerance,
				  settle,
				  particle_array_prev_shake,shake_virial);
      }
#endif
    }
    timer_stop(timer_settle);

    if(operations.doVirialcalculation){
      calc_local_kenergy_virial(particlearray,typerangearray,
				pairwise_virial,water_settle_virial, shake_virial,
				ke,lv);
    }else{
      calc_local_kenergy(particlearray,typerangearray,ke);
    }
  }

#ifdef PROF_CHECK
  /*
      if(unit_identifier==0){
	printf("MPI_Barrier before reduce %d\n",unit_identifier);
      }
  */
      MPI_Barrier(MPI_COMM_WORLD);
#endif
  timer_start(timer_rede);
  if(operations.doVirialcalculation){
    reduce(ke_sum,energy_sum,virial_sum,ke,energy,lv);
  }else{
    reduce(ke_sum,energy_sum,ke,energy);
  }
#ifdef CHECK_ENERGY
  {
    lje_sum = 0.0;
    se_sum = 0.0;
    zde_sum = 0.0;
    sse_sum = 0.0;
    le_sum = 0.0;
    exe_sum = 0.0;
    clj_sum = 0.0;
    ccl_sum = 0.0;
    reduce(se_sum, zde_sum, shortpair_energy,zerodipole_energy);
    reduce(lje_sum, sse_sum, ljpair_energy,calcforce.shortenergyself);
    reduce(le_sum, exe_sum, calcforce.longenergy,calcforce.shortenergyexcluded);
    reduce(ccl_sum,clj_sum,ShortRange::short_coulomb,ShortRange::short_lj);
    if(unit_identifier==0){
      printf("Coulomb energy short %e self %e long %e excluded %e\n",
	     se_sum+zde_sum,sse_sum,le_sum,exe_sum);
      printf("LJ energy %e\n",lje_sum);
      printf("CB cancel Coulomb %e LJ %e\n",ccl_sum,clj_sum);
    }
  }
#endif
  timer_stop(timer_rede);

  if(restored){   /// correct total potential energy with long-range nodes
    energy_sum = stored_total_penergy;
  }else{
    energy_sum+=tljcec;
  }
    
  temperature = 2.0 * ke_sum * UnitParameter::unitEnergy_eV * inv_boltzmann * inv_num_freedom;

  if(operations.doVirialcalculation){
    calculate_pressure(volume,ke_sum,virial_sum,pressure);
  }
  if(restored){
    if(temp_control.method == Config::kANDERSEN_HOOVER){// andersen_hoover
      if(ref_pressure==0.0)ref_pressure=pressure; //<! not given reference pressure 
      diff_pressure_volume = (pressure-ref_pressure)*volume;
    }
  }else{
    // nose_hoover
    if(temp_control.method == Config::kNOSE_HOOVER){
      eta_force = (2.0 * ke_sum) - (num_freedom * inv_unitEnergy_eV * UnitParameter::Boltzmann * ref_temperature); 
    }else if(temp_control.method == Config::kANDERSEN_HOOVER){// andersen_hoover
      eta_force = (2.0 * ke_sum) - ((1+num_freedom) * inv_unitEnergy_eV * UnitParameter::Boltzmann * ref_temperature)
	+ logv_mass*logv_vel*logv_vel;
      if(ref_pressure==0.0)ref_pressure=pressure; //<! not given reference pressure 
      diff_pressure_volume = (pressure-ref_pressure)*volume;
#ifdef MODIFIED_ANDERSEN
      //    double onedinvnf = 1.0+ 3.0/num_freedom; /// ???
      //    logv_force = (onedinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);
      double dinvnf = 3.0/num_freedom;
      logv_force = (dinvnf*2.0*ke_sum + 3.0*diff_pressure_volume);  /// p4178 eq(2.9)
#else
      logv_force = (3.0*diff_pressure_volume);  /// p4178 eq(2.7)
#endif
      if(unit_identifier==0){
	printf("ke_sum %f, virial_sum %f,  nkT/volume %f\n", ke_sum, virial_sum,num_freedom*inv_unitEnergy_eV*UnitParameter::Boltzmann*temperature/volume);
	printf("eta_pos %f, eta_vel %f, eta_force %f\n",eta_pos, eta_vel, eta_force);
	printf("logv_pos %f, logv_vel %f, logv_force %f, (pressure-ref_pressure)*volume %e, ref_pressure %e, pressure %f %e Pa %e atm\n",logv_pos, logv_vel, logv_force, (pressure-ref_pressure)*volume, ref_pressure, pressure,pressure*UnitParameter::unitEnergy_J/(UnitParameter::unitLength*UnitParameter::unitLength*UnitParameter::unitLength),pressure*UnitParameter::unitEnergy_J/(UnitParameter::unitLength*UnitParameter::unitLength*UnitParameter::unitLength)/101325.0);
	std::cout << "boxsize " << boxsize << " volume " << volume << " log(volume)/3.0 "<< log(volume)/3.0 << std::endl;
      }
    }
  }
  if(tmax>0){
    //    t = 0;
    backup_for_rollback(particlearray,
			typerangearray,bondlistarray,bondlistarray_idx,waterlist,shakelist,t);
  }

  last_time_stamp = getrealtime();
  if(unit_identifier==0)printf("Main loop start at %f sec step %d\n",last_time_stamp,t);
  //! main loop
  for(;t<tmax;t++){
    // DBUG
    //    dump_a_particle(particlearray,typerangearray,int(13271),int(t),unit_identifier);

    if(unit_identifier==0){
      if((DebugLog::verbose>1)
         ||((DebugLog::verbose==1)&&(pi==0))){
        std::cout << "step " << t << " E " << ke_sum+energy_sum<< " Ek " << ke_sum << " Ep " << energy_sum 
                  << " T " << temperature;
	if(operations.doVirialcalculation){
	  std::cout << " virial " << virial_sum << " P " << pressure;
	  if(temp_control.method == Config::kANDERSEN_HOOVER){
	    std::cout << " Lx " << boxsize.x;
	  }
	}
	std::cout << std::endl;
      }
    }

#if PRINT_DETAIL
    dump_detail(particlearray, waterlist, t);
#endif

    // copy previous pos and vel to particle_array_prev for water molecules
    backup_water_posvel(particlearray,waterlist,pv_array_prev);

#ifdef USE_SHAKE
    // copy previous pos to particle_array_prev_shake for shake/rattle
    if(shake_type > 0){
      backup_shake_posvel(particlearray,shakelist,particle_array_prev_shake);
    }
#endif

    // nose-hoover
    if(ri==reduce_interval-1){    // TODO fix ri overflow 
      if(temp_control.method == Config::kNOSE_HOOVER){
	nose_hoover_pre(particlearray, typerangearray,
			eta_pos, eta_vel, eta_force, ke_sum,
			eta_inv_mass, ref_temperature, inv_unitEnergy_eV, num_freedom, nxt_dt);
      }else if(temp_control.method == Config::kANDERSEN_HOOVER){
	if((unit_identifier==0)&&(DebugLog::verbose>1)){
	  dump_andersen_hoover(ke_sum, virial_sum, num_freedom, inv_unitEnergy_eV, temperature, volume,
			       eta_pos, eta_vel, eta_force, logv_pos, logv_vel, logv_force,
			       pressure, boxsize);
	}
	andersen_hoover_pre(particlearray, typerangearray,
			    eta_pos, eta_vel, eta_force, logv_vel, logv_force, ke_sum,
			    eta_inv_mass, logv_mass, logv_inv_mass,
			    ref_temperature, inv_unitEnergy_eV, num_freedom, diff_pressure_volume,
			    nxt_dt);
      }
    }

    // integration code
    if(temp_control.method == Config::kANDERSEN_HOOVER){
      if(operations.doCorrectTranslation){
	velocity_verlet_npt_pos_vel_vtcorr(particlearray,typerangearray,dt,inv_nt);
      }else{
	velocity_verlet_npt_pos_vel(particlearray,typerangearray,dt);
      }
      postprocess.change_boxsize(boxsize);
      tljcec = calcTotalLJCutoffEnergyCorrection(volume);
      if((unit_identifier==0)&&(DebugLog::verbose>1)){
	dump_andersen_hoover(ke_sum, virial_sum, num_freedom, inv_unitEnergy_eV, temperature, volume,
			     eta_pos, eta_vel, eta_force, logv_pos, logv_vel, logv_force,
			     pressure, boxsize);
      }
    }else{
      velocity_verlet_pos_vel(particlearray,typerangearray,dt);
    }
      // settle of rattle for water
    timer_start(timer_settle);
    if(temp_control.method == Config::kANDERSEN_HOOVER){
      settle.constrain_position(particlearray, pv_array_prev, waterlist, typerangearray[0],dt,water_settle_force);
    }else{
      settle.constrain_position(particlearray, pv_array_prev, waterlist, typerangearray[0],dt);
    }
    timer_stop(timer_settle);

#ifdef USE_SHAKE
    timer_start(timer_shake);
    if(shake_type > 0){

      //      covalent_bond_parameter_list = calcforce.covalentbond.getParameterList();
      if(temp_control.method == Config::kANDERSEN_HOOVER){
	if(shake_type==1){
	  clear_shake_force(shakelist,shake_force);
	}
	settle.rattle_position(particlearray,
			       particle_array_prev_shake,
			       shakelist,
			       calcforce.covalentbond.getParameterList(),
			       bondlistarray,
			       dt,
			       shake_max_iterate,
			       shake_tolerance,
			       shake_force);
      }else{
	settle.rattle_position(particlearray,
			       particle_array_prev_shake,
			       shakelist,
			       calcforce.covalentbond.getParameterList(),
			       bondlistarray,
			       dt,
			       shake_max_iterate,
			       shake_tolerance);
      }
    }
    timer_stop(timer_shake);
#endif

#ifdef PROF_CHECK
    /*
    if(unit_identifier==0){
      printf("MPI_Barrier after rattle\n");
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    //      select_move_particle(particlearray,);
    mi++;
    if(mi==moveinter){

      timer_start(timer_move);
      not_over_flow = move_top_half(particlearray, typerangearray, bondlistarray, waterlist, shakelist, communicator, postprocess, over_move);
      rollback = check_move_over_and_rollback(not_over_flow, over_move,
					      max_moveinter,
                                              particlearray,
                                              typerangearray,
                                              bondlistarray, bondlistarray_idx, waterlist, shakelist,
                                              t,
					      moveinter);
      if(rollback==false){
        if((unit_identifier==0)&&(DebugLog::verbose>0)){
          printf("Move particle\n");
        }
        // dump_bonds(unit_identifier,communicator,bondlistarray);
        not_over_flow = move_bottom_half(particlearray, typerangearray, bondlistarray, waterlist, shakelist, communicator, postprocess,sl);
	rollback = check_move_over_and_rollback(not_over_flow, over_move,
						max_moveinter,
						particlearray,
						typerangearray,
						bondlistarray, bondlistarray_idx, waterlist, shakelist,
						t,
						moveinter);
      }
      if(rollback==false){
        moveinter = max_moveinter;
        // Back up for roll back
        backup_for_rollback(particlearray,
			    typerangearray,bondlistarray,bondlistarray_idx,waterlist,shakelist,t);
      }else if(moveinter==0){
	if(unit_identifier==0){
	  printf("Can't retry rollback. Break.\n");
	}
	timer_stop(timer_move);
	result = false;
	break;
      }
#ifdef PROF_CHECK
    /*
    if(unit_identifier==0){
      printf("MPI_Barrier after move\n");
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
#endif
      timer_stop(timer_move);

      mi = 0;
      if(rollback){
	postprocess.change_boxsize(boxsize);
        moved = 0;
      }else{
        moved = 1;
      }
    }

    if(moved){
      clear_before_move(ghost, targettyperange, recvsetid, bondlistarray, calcforce);
      if(DebugLog::verbose>1){
        dump_atomid(particlearray,typerangearray,setid);
      }
    }
    clear_force(particlearray,typerangearray);
    clear_force(force,typerangearray);

#ifdef USE_PAIRLIST
    if(moved){
      operations.doMakePairList = true;
    }else{
      operations.doMakePairList = false;
    }
#endif

    mt_operations = operations;
    if(ri!=reduce_interval-1){
      mt_operations.doEnergycalculation = false;
      mt_operations.doVirialcalculation = false;
    }
    if(mt_bond==mti_bond){
      mt_bond = 0;
    }else{
      mt_operations.doCovalentBondcaculation = false;
      mt_bond++;
    }
    if(mt_short==mti_short){
      mt_short = 0;
      calcforce.selfshort_zero();
      calcforce.covalentbond.set_cancel_short();
      if(mt_operations.doVirialcalculation){
	calcforce.covalentbond.set_cancel_short_virial();
      }
    }else{
      mt_operations.doShortrangecalculation = false;
      mt_short++;
      calcforce.covalentbond.unset_cancel_short();
      if(mt_operations.doVirialcalculation){
	calcforce.covalentbond.unset_cancel_short_virial();
      }
    }
    if(mt_long==mti_long){
      mt_long = 0;
      calcforce.long_zero();
    }else{
      mt_operations.doLongrangecalculation = false;
      mt_operations.doShortLongCommunication = false;
      mt_long++;
    }
#ifdef TEST_CB_ONLY
    mt_operations.doShortrangecalculation = false;
    calcforce.covalentbond.unset_cancel_short();
    if(mt_operations.doVirialcalculation){
      calcforce.covalentbond.unset_cancel_short_virial();
    }
#endif

    pairwise_virial = 0.0;

#ifdef PROF_CHECK
    /*
    if(unit_identifier==0){
      printf("MPI_Barrier before exchangeParticleArray\n");
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef CPPMD_ENABLE_PMMM
    if(DebugLog::verbose>1)printf("calcPMMM_top_half\n");
    calcforce.calcPMMM_top_half(particlearray, typerangearray,operations);
    if(DebugLog::verbose>1)printf("calcPMMM_top_half done\n");
#endif
    timer_start(timer_exp);
    if(moved){
#ifdef OVERLAP
      exchangeParticleArray_top_half(particlearray, typerangearray, bondlistarray, setid,
				     ghost, targetparticlerange, 
				     recvsetid, recvsetid_to_index, 
				     targettyperange, targetbond, ghostcbindexarray, num_ghost,
				     calcforce, communicator, longreceiver, longsender,
				     reset_ghost_longset_index,
				     mt_operations);
      timer_stop(timer_exp);
      timer_start(timer_expcalc);
#ifdef CPPMD_ENABLE_PMMM
    if(DebugLog::verbose>1)printf("calcForce_local\n");
#endif
      calcforce.calcForce_local(particlearray, typerangearray, 
			bondlistarray,bondlistarray_idx, waterlist,
                        ghost, recvsetid, 
			recvsetid_to_index, targettyperange, targetbond,
			force, ghostforce, energy, pairwise_virial, mt_operations);
#ifdef CPPMD_ENABLE_PMMM
    if(DebugLog::verbose>1)printf("calcForce_local Done\n");
#endif
      timer_stop(timer_expcalc);
      timer_start(timer_exp);
      exchangeParticleArray_bottom_half(particlearray, typerangearray, bondlistarray, setid,
					ghost, targetparticlerange, 
					recvsetid, recvsetid_to_index, 
					targettyperange, targetbond, ghostcbindexarray, num_ghost,
					calcforce, communicator, longreceiver, longsender,
					reset_ghost_longset_index,
					mt_operations);
#else
      exchangeParticleArray(particlearray, typerangearray, bondlistarray, setid,
                            ghost, targetparticlerange, 
                            recvsetid, recvsetid_to_index, 
                            targettyperange, targetbond, ghostcbindexarray, num_ghost,
                            calcforce, communicator, longreceiver, longsender,
                            reset_ghost_longset_index,
			    mt_operations);
#endif
    }else{
#ifdef OVERLAP
      exchangeParticleArray_onlyPosition_top_half(particlearray, typerangearray, bondlistarray,
						  ghost, targetparticlerange, 
						  recvsetid, recvsetid_to_index, 
						  targettyperange, targetbond, num_ghost,
						  communicator, longreceiver, longsender,
						  mt_operations);
      timer_stop(timer_exp);
      timer_start(timer_expcalc);
      calcforce.calcForce_local(particlearray, typerangearray, 
			bondlistarray,bondlistarray_idx, waterlist,
                        ghost, recvsetid, 
			recvsetid_to_index, targettyperange, targetbond,
			force, ghostforce, energy, pairwise_virial, mt_operations);
      timer_stop(timer_expcalc);
      timer_start(timer_exp);
      exchangeParticleArray_onlyPosition_bottom_half(particlearray, typerangearray, bondlistarray,
						     ghost, targetparticlerange, 
						     recvsetid, recvsetid_to_index, 
						     targettyperange, targetbond, num_ghost,
						     communicator, longreceiver, longsender,
						     mt_operations);
#else
      exchangeParticleArray_onlyPosition(particlearray, typerangearray, bondlistarray,
                                         ghost, targetparticlerange, 
                                         recvsetid, recvsetid_to_index, 
                                         targettyperange, targetbond, num_ghost,
                                         communicator, longreceiver, longsender,
					 mt_operations);
#endif
    }
#ifdef PROF_CHECK
    /*
    if(unit_identifier==0){
      printf("MPI_Barrier after exchangeParticleArray\n");
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    timer_stop(timer_exp);

    postprocess_after_exchangeparticle(postprocess,
                                       ghost,recvsetid,recvsetid_to_index,targettyperange,
                                       ghostforce,energy,shiftidcheck);
    if(mt_short==mti_short){
      mt_short = 0;
      calcforce.ghostshort_zero();
    }


    timer_start(timer_calc);
#ifdef CPPMD_ENABLE_PMMM
    if(DebugLog::verbose>1)printf("calcForce\n");
#endif

#ifdef OVERLAP
    calcforce.calcForce(particlearray, typerangearray, 
			bondlistarray,bondlistarray_idx, waterlist,
                        ghost, recvsetid, 
			recvsetid_to_index, targettyperange, targetbond,
			force, ghostforce, energy, pairwise_virial, mt_operations,false);
#else
    calcforce.calcForce(particlearray, typerangearray, 
			bondlistarray,bondlistarray_idx, waterlist,
                        ghost, recvsetid, 
			recvsetid_to_index, targettyperange, targetbond,
			force, ghostforce, energy, pairwise_virial, mt_operations,!overlap);
#endif
#ifdef CPPMD_ENABLE_PMMM
    if(DebugLog::verbose>1)printf("calcForce Done\n");
#endif
#ifdef PROF_CHECK
    /*
      if(unit_identifier==0){
	printf("MPI_Barrier after calcForce %d\n",unit_identifier);
      }
    */
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    timer_stop(timer_calc);

    if(operations.doExchangeForce){

      timer_start(timer_exf);
      exchangeForceArray(ghostforce, targetparticlerange, force,
                         calcforce.ghost_longset_index,
                         targettyperange, ghostcbindexarray,
                         communicator, longreceiver, longsender, operations, moved);
      timer_stop(timer_exf);

      mergerForces(particlearray);
    }

    // integration code
    velocity_verlet_velocity(particlearray,typerangearray,dt);

    timer_start(timer_settle);
    // settle for water
    if(operations.doVirialcalculation){
      water_settle_virial = 0.0;
      settle.constrain_velocity(particlearray, waterlist, typerangearray[0],dt,water_settle_virial);
    }else{
      settle.constrain_velocity(particlearray, waterlist, typerangearray[0]);
    }
    timer_stop(timer_settle);

#ifdef USE_SHAKE
    timer_start(timer_shake);
    if(shake_type > 1){
      if(operations.doVirialcalculation){
	shake_virial = 0.0;
	settle.rattle_velocity(particlearray,
			       shakelist,
			       calcforce.covalentbond.getParameterList(),
			       bondlistarray,
			       dt,
			       shake_max_iterate,
			       shake_tolerance,
			       shake_virial);
      }else{
	settle.rattle_velocity(particlearray,
			       shakelist,
			       calcforce.covalentbond.getParameterList(),
			       bondlistarray,
			       dt,
			       shake_max_iterate,
			       shake_tolerance);
      }
    }
    timer_stop(timer_shake);
#endif

#ifdef DEBUG_SHORTENERGYEXCLUDE
  {
    double le_sum = 0.0;
    double exe_sum = 0.0;
    reduce(le_sum, exe_sum, calcforce.longenergy,calcforce.shortenergyexcluded);
    if(unit_identifier==0){
      printf("Energy of excluded pair %24.16e\n",exe_sum);
    }
  }
#endif
    //reduce
    ri++;
    if(ri==reduce_interval){    // TODO fix ri overflow 
#ifdef DEBUG_MOMENTUM
      dump_momentum(unit_identifier,particlearray,typerangearray);
#endif
      
      if(operations.doVirialcalculation){
	calc_local_kenergy_virial(particlearray,typerangearray,
				  pairwise_virial,water_settle_virial, shake_virial,
				  ke,lv);
      }else{
	calc_local_kenergy(particlearray,typerangearray,ke);
      }

#ifdef PROF_CHECK
      /*
      if(unit_identifier==0){
	printf("MPI_Barrier before reduce %d\n",unit_identifier);
      }
      */
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      timer_start(timer_rede);
      if(operations.doVirialcalculation){
	//	reduce(ke_sum,energy_sum,virial_sum,ke,energy,lv);
	reduce(ke_sum,energy_sum,virial_sum,ke,energy,lv);
      }else{
	reduce(ke_sum,energy_sum,ke,energy);
      }
      timer_stop(timer_rede);

      energy_sum+=tljcec;
      // calc temperature
      temperature = 2.0 * ke_sum * UnitParameter::unitEnergy_eV * inv_boltzmann * inv_num_freedom;
//      std::cout << "temperature (before) = " << temperature << std::endl;

      if(operations.doVirialcalculation){
	calculate_pressure(volume,ke_sum,virial_sum,pressure);
      }
      // temperature control
      switch (temp_control.method){
        case Config::kNONE:{
          break;
        }
        case Config::kVEL_SCALING:{
          // velocity scaling
          velocity_scaling(particlearray, typerangearray,
                           ke_sum, temperature, ref_temperature);
          break;
        }
        case Config::kNOSE_HOOVER:{
          // nose-hoover
          nose_hoover_post(particlearray, typerangearray,
                           eta_pos, eta_vel, eta_force, ke_sum, temperature,
                           eta_inv_mass, ref_temperature, inv_unitEnergy_eV, num_freedom, nxt_dt);
          break;
        }
        case Config::kLANGEVIN:{
          break;
        }
        case Config::kANDERSEN_HOOVER:{
	  diff_pressure_volume = (pressure-ref_pressure)*volume;
	  andersen_hoover_post(particlearray, typerangearray,
			       eta_pos, eta_vel, eta_force, 
			       logv_vel, logv_force,
			       ke_sum, temperature,
			       eta_inv_mass, logv_mass, logv_inv_mass,
			       ref_temperature, inv_unitEnergy_eV, num_freedom, diff_pressure_volume,
			       nxt_dt);
	  break;
        }
        default: break;
      }
      ri=0;
    }//ri
    
    pi++;
    if(pi==printinterval) {
      pi=0;
    }
    dci++;
    if(dci==dump_crd_inter){
      timer_start(timer_crd);
      if(temp_control.method==Config::kANDERSEN_HOOVER){
	dump.setBoxSize(boxsize);
      }
      dump.GatherDumpAmberCRD(particlearray,typerangearray);
      
      dci=0;
      timer_stop(timer_crd);
    }
    dri++;
    if(dri==dump_rst_inter){
      timer_start(timer_rst);
      if(temp_control.method==Config::kANDERSEN_HOOVER){
	dump_rst.setBoxSize(boxsize);
      }
      dump_rst.GatherDumpAmberCRD(particlearray,typerangearray);

      dri=0;
      timer_stop(timer_rst);
    }
    if(moved)moved=0;

    if((t+1)%time_stamp_interval==0){
      time_stamp = getrealtime();
      lap_time = time_stamp - last_time_stamp;
      if(unit_identifier==0)printf("Time step %d end at %f sec, lap time %f sec / %d step\n",t,time_stamp,lap_time,time_stamp_interval);
      last_time_stamp = time_stamp;
    }
  }//main loop
  time_stamp = getrealtime();
  if(unit_identifier==0)printf("main loop step %d end at %f sec\n",t,time_stamp);
  last_time_stamp = time_stamp;


  calc_total_energy(particlearray,typerangearray,tljcec,energy,ke,energy_sum,ke_sum,force_sum);
  if((restored)&&(tmax<=0)){  /// correct total potential energy with long-range nodes
    energy_sum = stored_total_penergy;
  }
  if(unit_identifier==0){
    temperature = 2.0 * ke_sum * UnitParameter::unitEnergy_eV * inv_boltzmann * inv_num_freedom;
    std::cout << "step " << t << " E " << ke_sum+energy_sum << " Ek " << ke_sum << " Ep " << energy_sum << " T " << temperature;
    if(operations.doVirialcalculation){
      std::cout << " virial " << virial_sum << " P " << pressure;
      if(temp_control.method == Config::kANDERSEN_HOOVER){
	std::cout << " Lx " << boxsize.x;
      }
    }
    std::cout << " F " << force_sum << std::endl;
    std::cout << "Energy(kcal/mol) Etot " << UnitParameter::realizeEnergy_kcal_mol((ke_sum+energy_sum));
    std::cout << " Ek " << UnitParameter::realizeEnergy_kcal_mol(ke_sum);
    std::cout << " Ep " << UnitParameter::realizeEnergy_kcal_mol(energy_sum);
    std::cout << " T " << temperature;
    if(operations.doVirialcalculation){
      std::cout << " P " << pressure*UnitParameter::unitPressure_Pa << " Pa " << pressure*UnitParameter::unitPressure_Pa/101325.0 << " atm";
    }
    std::cout << std::endl;
    std::cout << "Energy(kcal/mol) per atom : Etot " << UnitParameter::realizeEnergy_kcal_mol((ke_sum+energy_sum)*inv_nt);
    std::cout << " Ek " << UnitParameter::realizeEnergy_kcal_mol(ke_sum*inv_nt);
    std::cout << " Ep " << UnitParameter::realizeEnergy_kcal_mol(energy_sum*inv_nt);
    std::cout << std::endl;
    if(temp_control.method == Config::kANDERSEN_HOOVER){
      dump_andersen_hoover(ke_sum, virial_sum, num_freedom, inv_unitEnergy_eV, temperature, volume,
			   eta_pos, eta_vel, eta_force, logv_pos, logv_vel, logv_force,
			   pressure, boxsize);
    }
    std::cout << std::endl;
  }
  if((DebugLog::verbose>2)||(DebugLog::particle_dump>0)){
    std::cout << "step " << t << std::endl;
    dump_particle(particlearray,typerangearray,setid,energy);
  }

  timer_stop(timer1);
#ifdef FPCOLL
#ifdef FPCOLL_INTEGRATOR
  fpcoll_stop();
#endif
#endif
  SAMPLE_STOP();

  last_t = t;
  if((result==false)&&(t>0)){
    if(unit_identifier==0){
      printf("Force write_restart last rollback point %ld\n",t);
    }
    write_restart=1;
  }
#ifdef USE_HDF
  if(operations.doShortrangecalculation){
    time_stamp = getrealtime();
    if(unit_identifier==0)printf("dump HDF start at %f sec\n",time_stamp);
    last_time_stamp = time_stamp;
    int comm_size;
    MPI_Comm_size(communicator.mpi_comm_short, &comm_size);
    dump_hdf_snapshot(particlearray, typerangearray, bondlistarray, bondlistarray_idx,
		      waterlist, shakelist,
		      ke, energy, lv, energy_sum, t, num_total, comm_size);
    close_hdfdump();
    time_stamp = getrealtime();
    if(unit_identifier==0)printf("dump HDF end at %f sec, spent %f sec\n",time_stamp-last_time_stamp);
    last_time_stamp = time_stamp;
  }
#else
#ifdef BINARY_DUMP
  if(operations.doShortrangecalculation){
    time_stamp = getrealtime();
    if(unit_identifier==0)printf("dump binary start at %f sec\n",time_stamp);
    last_time_stamp = time_stamp;
    dump_binary_snapshot(particlearray,
			 typerangearray,bondlistarray,bondlistarray_idx,waterlist,shakelist,ke,energy,lv,energy_sum,t);
    close_binarydump();
    time_stamp = getrealtime();
    if(unit_identifier==0)printf("dump binary end at %f sec, spent %f sec\n",time_stamp-last_time_stamp);
    last_time_stamp = time_stamp;
  }
#endif
#endif
  return result;
}

#ifdef OLDPARTICLE
template class Integrator<ParticleArray,ParticleArray>;
#else
template class Integrator<CombinedParticleArray,GhostParticleArray>;

#endif
