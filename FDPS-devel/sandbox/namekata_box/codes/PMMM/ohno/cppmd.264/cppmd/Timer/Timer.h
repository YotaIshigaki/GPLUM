//! Timer and counter
/*
  Measure realtime, usertime and systemtime.
  If USE_PAPI is defined, also measure 5 PAPI_COUNTERs, 
  PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_FP_INS, PAPI_BR_INS, PAPI_L2_TCM.
  10 block/region are counted.
 */
#ifndef TIMER_H
#define TIMER_H

#ifdef USE_PAPI
#include "papi.h"
#endif

#include <sys/time.h>
#include <sys/resource.h>

#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <iostream>

double getrealtime();

//! set result of getrusage to usage and return realtime
/*
  If USE_PAPI is defined, realtime is PAPI_get_real_usec().
  Otherwise, calculate from gettimeofday
 */
long long get_realtime_usec(struct rusage& usage);

namespace PerfCounter {
  extern std::vector< std::set<int> > low_target;  //! nested target
  extern std::vector< std::string > target_name;   //! name of target
  
  extern int num_target;                //! number of target block/region

  //! initialize, nt is number of target block/region
  /* must be call from constructor in this version */
  void init(int nt);

  //! add target name, and return target identifier
  int
  add_target(std::string name);

  //! set EventSet. 
  /* must be call from init in this version */
  void set_event();

  //! stop count, for internal use
  inline
  long long simple_stop();

  //! start count, for internal use
  inline
  long long simple_start();

  //! start count, t is identifier of target block/region. It can be nested.
  /*!
    If nested, previous counter paused untile this count is stopped by stop().
   */
  long long start(int t);

  //! stop count latest(deepest) block.
  /*!
    If target is nested, resume counting upper block.
   */
  long long stop();

  //! print result to stdout.
void print(std::ostream& ofs=std::cout);
void print_time(std::ostream& osf=std::cout);
void print(int& rank, std::string namebase=std::string("timers"));

double * get_time_array(int& num);
void put_time_array(double *times, int &num);
inline const char * get_timer_name(const int timer_id){
  if((timer_id>=0)&&(timer_id<num_target)){
    return target_name[timer_id].c_str();
  }else{
    return (const char *)NULL;
  }
}
}

#endif
