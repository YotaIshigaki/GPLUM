#ifndef PMMMINTERACTION_H
#define PMMMINTERACTION_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <algorithm>

#include "Common.h"
#include "ParticleInfo.h"

#include "pmmm_fmm.h"
#include "M2L_convolution_PBC.h"

#include "MPIParallelLongRangePMMM.h"

#include "Timer.h"

#ifdef PMMM_TIMER_DETAIL
#define start_timer(t) PerfCounter::start(t)
#define stop_timer()  PerfCounter::stop()
#else
#define start_timer(t) 
#define stop_timer()  
#endif

static inline void mm_decompose(std::vector<int> &m_boundary, int num_mm, int order)
{
  int size = (order+1)*(order+1);
  int nm = size/num_mm;
  int rm = size-nm*num_mm;
  int m;
  int t = 0;
  m_boundary.resize(num_mm+1);
  m_boundary[0] = 0;
  for(m=0;m<num_mm;m++){
    m_boundary[m+1] = m_boundary[m] + nm;
    if(m<rm)m_boundary[m+1]++;
  }

}

template<typename real_t=double>
class PMMM_PM {
 public:
  std::vector<MultipoleMoment<> > mm;
  std::vector<LocalExpansion<> > le;
  std::vector<Position> center;
  SpaceVector<double> dipole;
  double quad0;
  int num_cell;
  int order;
 double scale;
 double i_scale;

 MPI_Comm pp_comm;

 MPISenderForLongPMMM<> msender_lreceiver;

 PMMM_PM(int p=5) : order(p)
  {
    Rlm<double> rlm(order);
    rlm.init();
  }

 void set_order(int p){
   order = p;
   Rlm<double> rlm(order);
   rlm.init();
 }

 void set_cell_info(const std::vector<Position> cell_center, const int local_cell_size[3], const double _scale){
   scale = _scale;
   i_scale = 1.0/_scale;
   num_cell = cell_center.size();
   mm.resize(num_cell,order);
   for(int i=0;i<num_cell;i++)mm[i].set_order(order);
   le.resize(num_cell,order);
   for(int i=0;i<num_cell;i++)le[i].set_order(order);
   center.resize(num_cell);
   for(int i=0;i<num_cell;i++)center[i]=cell_center[i]*scale;
   msender_lreceiver.set_size(local_cell_size,(order+1)*(order+1));
 }

 void set_comm_info(const MPI_Comm &short_comm, const MPI_Comm &pmmm_comm, const std::vector<int> &mm_targetid, const std::vector<int> &m_boundary){
   pp_comm = short_comm;
   msender_lreceiver.set_comm(pmmm_comm);
   msender_lreceiver.set_long_target(mm_targetid,m_boundary);
 }

  template<class PA>
  void do_P2M(const PA& particlearray, const std::vector<TypeRange>& typerange){
    int cimax = typerange.size();
#pragma omp parallel for
    for(int ci=0;ci<cimax;ci++){
      mm[ci].clear();
      int i;
      int imin, imax;
      imin = typerange[ci].begin;
      imax = typerange[ci].end;
      for(i=imin;i<imax;i++){
	mm[ci].assign_particle(center[ci],particlearray.poscharge[i].position*scale,particlearray.poscharge[i].charge);
      }
    }
  }
 template<class PA>
  void calc_dipole_for_diplple_correction(const PA& particlearray, const std::vector<TypeRange>& typerange){
   double x = 0.0;
   double y = 0.0;
   double z = 0.0;
   double q = 0.0;
   int cimax = typerange.size();
#pragma omp parallel for reduction(+:x,y,z,q)
   for(int ci=0;ci<cimax;ci++){
      x += mm[ci].buf[3];
      y += mm[ci].buf[1];
      z += mm[ci].buf[2];
      int i;
      int imin, imax;
      imin = typerange[ci].begin;
      imax = typerange[ci].end;
      for(i=imin;i<imax;i++){
	Position dr = particlearray.poscharge[i].position*scale - center[ci];
	q += particlearray.poscharge[i].charge * (dr*dr);
      }
    }
   SpaceVector<double> dp(x,y,z);
    /// Allreduce dipole and quad0
   MPI_Allreduce(&(dp[0]), &(dipole[0]), 3, MPI_DOUBLE, MPI_SUM, pp_comm);
    MPI_Allreduce(&q, &quad0, 1, MPI_DOUBLE, MPI_SUM, pp_comm);
   dipole *= (4.0/3.0)*M_PI;
   {
     if(DebugLog::verbose>1){
       int pr;
       MPI_Comm_rank(pp_comm,&pr);
       if(pr==0){
	 printf("dipole %e %e %e , quad0 %e\n",dipole.x,dipole.y,dipole.z,quad0);
       }
     }
   }
  }

 void pm_to_mm()
 {
   msender_lreceiver.sendM(mm);
 }

 void mm_progress()
 {
   msender_lreceiver.sendM_wait();
   msender_lreceiver.recvL();
 }

 void mm_to_pm()
 {
   msender_lreceiver.recvL_wait(le);
 }

  void dipole_correction(const double alpha){
    int i;
#pragma omp parallel for
    for(i=0;i<num_cell;i++){
      le[i].buf[3] += 2.0*dipole.x;
      le[i].buf[1] -= 2.0*dipole.y;
      le[i].buf[2] += 1.0*dipole.z;
      le[i].buf[0] += ((2.0/3.0)*M_PI)*quad0;
      le[i].buf[0] -= alpha*i_scale*M_2_SQRTPI*mm[i].buf[0];
    }
  }

 template<class PA>
 void do_L2P(const PA& particlearray, const std::vector<TypeRange>& typerange,
	      ForceArray& force, double& energy){
    int cimax = typerange.size();
    double e = 0.0;
    const double scale2 = scale*scale;
#pragma omp parallel for reduction(+:e)
    for(int ci=0;ci<cimax;ci++){
      int i;
      int imin, imax;
      imin = typerange[ci].begin;
      imax = typerange[ci].end;
      for(i=imin;i<imax;i++){
	LocalExpansion<> le1(1);
	le1.assign_from_LE(le[ci],particlearray.poscharge[i].position*scale, center[ci]);
	double phi, ax, ay, az;
	le1.read_phi_and_grad(phi, ax, ay, az);
	force[i].x -= ax*scale2*particlearray.poscharge[i].charge;
	force[i].y -= ay*scale2*particlearray.poscharge[i].charge;
	force[i].z -= az*scale2*particlearray.poscharge[i].charge;
	e+=phi*particlearray.poscharge[i].charge;
      }
    }

    if(DebugLog::verbose>1) {
      double esend[1], esum[1];
      esend[0] = e*scale;
      MPI_Allreduce(esend,esum,1,MPI_DOUBLE,MPI_SUM,pp_comm);
      int r;
      MPI_Comm_rank(pp_comm,&r);
      if(r==0)printf("L2P %24.16e\n",esum[0]);
    }
    energy += e*scale;
  }
 template<class PA>
  void do_L2P_corr(const PA& particlearray, const std::vector<TypeRange>& typerange,
		   ForceArray& force, double& energy,
		   const double csum, const double alpha){
    const double ainv2 = M_PI / (alpha * alpha*i_scale*i_scale);
    const double mc = ((2./3.) * M_PI) * csum;
    int cimax = typerange.size();
    //Debug
    double ed=0.0, ec=0.0;
    double e = 0.0;
    const double scale2 = scale*scale;
#pragma omp parallel for reduction(+:ed,ec)
    for(int ci=0;ci<cimax;ci++){
      int i;
      int imin, imax;
      imin = typerange[ci].begin;
      imax = typerange[ci].end;
      for(i=imin;i<imax;i++){
	const Position dr = particlearray.poscharge[i].position*scale - center[ci];
	force[i] -= (2.0 * mc) * dr*scale2*particlearray.poscharge[i].charge;
	ed += mc * (dr*dr)*particlearray.poscharge[i].charge;
	ec -= (csum*ainv2)*particlearray.poscharge[i].charge;

      }
    }
    e = ed+ec;
    // DEBUG
    //    printf("L2P_corr %e %e : %e\n",ed*scale,ec*scale,e*scale);
    if(DebugLog::verbose>1) {
      double esend[2], esum[2];
      esend[0] = ed*scale;
      esend[1] = ec*scale;
      MPI_Allreduce(esend,esum,2,MPI_DOUBLE,MPI_SUM,pp_comm);
      int r;
      MPI_Comm_rank(pp_comm,&r);
      if(r==0)printf("L2P_corr %24.16e %24.16e : %24.16e\n",esum[0],esum[1],esum[0]+esum[1]);
    }
    energy += e*scale;
  }
 
};

template <typename real_t=double>
class PMMM_MM {
 public:
  std::vector< std::vector<real_t> > mm;
  std::vector< std::vector<real_t> > le;
  M2L_convolution_PBC m2l_conv_pbc;
  GreenFunction_PBC gf;
  int order;
  int mm_length;
  int cell_size[3];
  int num_cell;
 int num_cell_cplx;
 double scale;
 double i_scale;

  int m_begin;
  int m_end;
  int num_local_m;
 int mm_id;

  int timer_mm_conv_forward;
  int timer_mm_conv_exchange;
  int timer_mm_conv_transform;
  int timer_mm_conv_backward;

 MPIReceiverForLongPMMM<> mreceiver_lsender;
 MPIMexchagePMMM<> mmexchanger;

 void set_size(const int global_cell_size[3], const int p, const double _scale){
    cell_size[0] = global_cell_size[0];
    cell_size[1] = global_cell_size[1];
    cell_size[2] = global_cell_size[2];
    mm_length = (p+1)*(p+1);
    num_cell = cell_size[2]*cell_size[1]*cell_size[0];
    num_cell_cplx = cell_size[2]*cell_size[1]*(1+cell_size[0]/2);
    m2l_conv_pbc.set_size(cell_size, p);
    gf.set_size(cell_size, p);
    mmexchanger.set_size(num_cell_cplx);
    scale = _scale;
    i_scale = 1.0/_scale;
  }

 void set_mm_info(const MPI_Comm &mm_comm, const std::vector<int> &m_boundary, 
		  const double alpha, const double cell_length, 
		  const int nmax, const int mmax, const int icut)
  {
    mmexchanger.set_comm(mm_comm);
    mm_id = mmexchanger.mm_id;
    m_begin = m_boundary[mm_id];
    m_end = m_boundary[mm_id+1];
    num_local_m = m_end-m_begin;
    mm.resize(num_local_m);
    for(int m=0;m<num_local_m;m++)mm[m].resize(num_cell);
    le.resize(num_local_m);
    for(int m=0;m<num_local_m;m++)le[m].resize(num_cell);
    m2l_conv_pbc.set_mm_range(m_begin,m_end);
    if((mm_id==0)||(DebugLog::verbose>1)){
      printf("scaled alpha cell_length %e %e\n",alpha*i_scale, cell_length*scale);
    }
#if defined(FFTW3_OPENMP) && defined(_OPENMP)
    fftw_init_threads();
    if((mm_id==0)||(DebugLog::verbose>1)){
      printf("fftw_plan_with_nthreads(%d)\n",omp_get_max_threads());
    }
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    double st = getrealtime();
    gf.gen_gf_r(alpha*i_scale, cell_length*scale, nmax, mmax, icut);
    double et = getrealtime();
    if((mm_id==0)||(DebugLog::verbose>1)){
      printf("gen_gf_r %f sec\n",et-st);
    }
    st = getrealtime();
    gf.gen_gf_k();
    et = getrealtime();
    if((mm_id==0)||(DebugLog::verbose>1)){
      printf("gen_gf_k %f sec\n",et-st);
    }
    mmexchanger.set_m_boundary(m_boundary);
    mmexchanger.set_buf_param();
  }

 void set_pp_info(const MPI_Comm &pmmm_comm, const std::vector<int> &pp_targetid, const std::vector<CellMethodModule::CellRange> &cell_range){
   mreceiver_lsender.set_size(cell_size,num_local_m);
   mreceiver_lsender.set_comm(pmmm_comm,mmexchanger.mm_id);
   mreceiver_lsender.set_target(pp_targetid,cell_range);
 }

  void register_timers() {
#ifdef PMMM_TIMER_DETAIL
    timer_mm_conv_forward = PerfCounter::add_target(std::string("MM_CONV_FORWARD"));
    timer_mm_conv_exchange = PerfCounter::add_target(std::string("MM_CONV_EXCHANGE"));
    timer_mm_conv_transform = PerfCounter::add_target(std::string("MM_CONV_TRANSFORM"));
    timer_mm_conv_backward = PerfCounter::add_target(std::string("MM_CONV_BACKWARD"));
#endif
  }

 void pm_to_mm()
 {
   mreceiver_lsender.recvM();
   mreceiver_lsender.recvM_wait(mm);
 }

 void mm_to_pm()
 {
   mreceiver_lsender.sendL(le);
   mreceiver_lsender.sendL_wait();
 }

  void convolution_PBC(){
    if(DebugLog::verbose>1)printf("forward\n");
    start_timer(timer_mm_conv_forward);
    m2l_conv_pbc.forward(mm);
    stop_timer();

    start_timer(timer_mm_conv_exchange);
    if(mmexchanger.num_long>1){
      if(DebugLog::verbose>2)
	{
	printf("pre mm_k[0]");
	for(int m=0;m<m2l_conv_pbc.mm_k[0].length;m++){
	  printf(" %e",real(m2l_conv_pbc.mm_k[0].buf[m]));
	}
	printf("\n");
	printf("exchangeM %d : %d %d \n",mmexchanger.mm_id,mmexchanger.m_count[mm_id],mmexchanger.m_displs[mm_id]);
	}
      mmexchanger.exchangeM(m2l_conv_pbc.mm_k);
    }
    stop_timer();

    if(DebugLog::verbose>2)
      {
	printf("send_buf[m*]");
	for(int m=0;m<mmexchanger.m_boundary[mmexchanger.mm_id+1]-mmexchanger.m_boundary[mmexchanger.mm_id];m++){
	  printf(" %e",real(mmexchanger.send_buf[m*mmexchanger.total_cell]));
	}
	printf("\n");
	printf("mm_k[0]");
	for(int m=0;m<m2l_conv_pbc.mm_k[0].length;m++){
	  printf(" %e",real(m2l_conv_pbc.mm_k[0].buf[m]));
	}
	printf("\n");
      }
    if(DebugLog::verbose>1)printf("transform\n");
    start_timer(timer_mm_conv_transform);
    m2l_conv_pbc.transform(gf);
    stop_timer();
    if(DebugLog::verbose>1)printf("backward\n");
    start_timer(timer_mm_conv_backward);
    m2l_conv_pbc.backward(le);
    stop_timer();
  }
};

class PMMMLocalInteraction {
 public:
  PMMM_PM<double> pmmm_pm;
  int unit_identifier;
  double alpha;
  MPI_Comm comm;
  int order;
  double scale;
  double i_scale;

  double csum;
  bool charge_neutral;

  int timer_pm_p2m;
  int timer_pm_dipole;
  int timer_pm_pm_to_mm;
  int timer_pm_mm_progress;
  int timer_pm_mm_to_pm;
  int timer_pm_dcorr;
  int timer_pm_l2p;
  int timer_pm_l2p_corr;

 PMMMLocalInteraction(int ui, const LongRangeParameter& param,  MPI_Comm short_comm ) 
   :  pmmm_pm(param.order), unit_identifier(ui), alpha(param.alpha), order(param.order), comm(short_comm),
    csum(0.0), charge_neutral(true)
  {
    calc_scale(param.boxSize);
  }

  PMMMLocalInteraction()
    {
    }

  void set_order(int p){
    pmmm_pm.set_order(p);
  }

  double calc_scale(const Position &box){
    i_scale = std::max(std::max(box.x,box.y),box.z);
    scale = 1.0/i_scale;
    if((unit_identifier==0)||(DebugLog::verbose>1)){
      printf("PMMM scale %f\n",scale);
    }
  }

  void settings(int ui, const LongRangeParameter& param,  MPI_Comm short_comm)
  {
    if((ui==0)||(DebugLog::verbose>1)){
      printf("PMMMLocalInteraction::settings\n");
    }
    pmmm_pm.set_order(param.order);
    unit_identifier = ui;
    alpha = param.alpha;
    order = param.order;
    comm = short_comm;
    csum = 0.0;
    charge_neutral = true;
    calc_scale(param.boxSize);
  }
  
  void initialize(const std::vector<Position> &cell_center, const int local_cell_size[3],
		  const MPI_Comm &pmmm_comm, 
		  const std::vector<int> &mm_targetid, const std::vector<int> &m_boundary){
    pmmm_pm.set_cell_info(cell_center,local_cell_size,scale);
    pmmm_pm.set_comm_info(comm,pmmm_comm,mm_targetid, m_boundary);
  }

  void not_charge_neutral(double charge){
    csum = charge;
    charge_neutral = false;
  }

  void register_timers() {
#ifdef PMMM_TIMER_DETAIL
    timer_pm_p2m = PerfCounter::add_target(std::string("PM_P2M"));
    timer_pm_dipole = PerfCounter::add_target(std::string("PM_DIPOLE"));
    timer_pm_pm_to_mm = PerfCounter::add_target(std::string("PM_PM_TO_MM"));
    timer_pm_mm_progress = PerfCounter::add_target(std::string("PM_MM_PROGRESS"));
    timer_pm_mm_to_pm = PerfCounter::add_target(std::string("PM_MM_TO_PM"));
    timer_pm_dcorr = PerfCounter::add_target(std::string("PM_DCORR"));
    timer_pm_l2p = PerfCounter::add_target(std::string("PM_L2P"));
    timer_pm_l2p_corr = PerfCounter::add_target(std::string("PM_L2P_CORR"));
#endif
  }

  template<class PA>
    void calc_local_top_half(const PA& particlearray, const std::vector<TypeRange>& typerange)
  {
    start_timer(timer_pm_p2m);
    pmmm_pm.do_P2M(particlearray,typerange);
    stop_timer();
    start_timer(timer_pm_pm_to_mm);
    pmmm_pm.pm_to_mm();
    stop_timer();
    start_timer(timer_pm_dipole);
    pmmm_pm.calc_dipole_for_diplple_correction(particlearray,typerange);
    stop_timer();
  }
  
  void calc_long_progress(){
    start_timer(timer_pm_mm_progress);
    pmmm_pm.mm_progress();
    stop_timer();
  }
  
  template<class PA>
    void calc_local_bottom_half(const PA& particlearray, const std::vector<TypeRange>& typerange,
				ForceArray& force, double& energy)
    {
    start_timer(timer_pm_mm_to_pm);
      pmmm_pm.mm_to_pm();
    stop_timer();
    start_timer(timer_pm_dcorr);
      pmmm_pm.dipole_correction(alpha);
    stop_timer();
    start_timer(timer_pm_l2p);
      pmmm_pm.do_L2P(particlearray,typerange, force, energy);
    stop_timer();
    start_timer(timer_pm_l2p_corr);
      if(!charge_neutral){
	pmmm_pm.do_L2P_corr(particlearray,typerange, force, energy, csum, alpha);
      }
    stop_timer();
    }
};


class PMMMLongRangeInteraction {
 public:
  PMMM_MM<double> pmmm_mm;

  int unit_identifier;
  double alpha;
  MPI_Comm comm;
  int order;
  double scale;
  double i_scale;

  int timer_mm_pm_to_mm;
  int timer_mm_mm_to_pm;
  int timer_mm_conv;
  
 PMMMLongRangeInteraction(int ui, const LongRangeParameter& param, MPI_Comm long_comm ) 
   :  unit_identifier(ui), alpha(param.alpha), order(param.order), comm(long_comm)
  {
    calc_scale(param.boxSize);
  }

  PMMMLongRangeInteraction(){
  }

  double calc_scale(const Position &box){
    i_scale = std::max(std::max(box.x,box.y),box.z);
    scale = 1.0/i_scale;
  }
  
  void settings(int ui, const LongRangeParameter& param, const MPI_Comm &long_comm){
    unit_identifier = ui;
    alpha = param.alpha;
    order = param.order;
    comm = long_comm;
    calc_scale(param.boxSize);
  }

  void initialize(const int global_cell_size[3], 
		  const std::vector<int> &m_boundary, 
		  const double cell_length, 
		  const int nmax, const int mmax, const int icut,
		  const MPI_Comm &pmmm_comm, const std::vector<int> &pp_targetid, 
		  const std::vector<CellMethodModule::CellRange> &cell_range){
    pmmm_mm.set_size(global_cell_size, order, scale);
    pmmm_mm.set_mm_info(comm, m_boundary, alpha, cell_length, nmax, mmax, icut);
    pmmm_mm.set_pp_info(pmmm_comm, pp_targetid, cell_range);
  }

  void register_timers() {
#ifdef PMMM_TIMER_DETAIL
    timer_mm_pm_to_mm = PerfCounter::add_target(std::string("MM_PM_TO_MM"));
    timer_mm_mm_to_pm = PerfCounter::add_target(std::string("MM_MM_TO_PM"));
    timer_mm_conv = PerfCounter::add_target(std::string("MM_CONV"));
#endif
    pmmm_mm.register_timers();
  }

  void prep()
  {
    start_timer(timer_mm_pm_to_mm);
    pmmm_mm.pm_to_mm();
    stop_timer();
  }

  void calcForce()
  {
    start_timer(timer_mm_conv);
    pmmm_mm.convolution_PBC();
    stop_timer();
  }

  void post()
  {
    start_timer(timer_mm_mm_to_pm);
    pmmm_mm.mm_to_pm();
    stop_timer();
  }
};

#endif
