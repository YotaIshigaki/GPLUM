#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include "Timer.h"


inline
long long timeval_to_ll(timeval tv)
{
  return (long long)tv.tv_usec + (long long)(tv.tv_sec)*1000000L;
}

double getrealtime()
{
  double time;

#ifdef NO_CLOCK_GETTIME
  struct timeval tv;
  gettimeofday(&tv,(struct timezone *)NULL);
  time = double(tv.tv_sec) + double(tv.tv_usec)*1.0e-6;
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  time = double(ts.tv_sec) + double(ts.tv_nsec)*1.0e-9;
#endif

  return time;
}

long long get_realtime_usec(struct rusage& usage)
{
  struct timeval tv;

  getrusage(RUSAGE_SELF,&usage);
#ifdef USE_PAPI 
  return PAPI_get_real_usec();
#else
  gettimeofday(&tv,(struct timezone *)NULL);
  return timeval_to_ll(tv);
#endif
}

namespace PerfCounter {
  int retval;
#ifdef USE_PAPI 
#define NUM_EVENTS 5
  int NumEvents;                    //! number of PAPI counter, 5 for x86
  int EventSet;                     //! event set of PAPI. 
  unsigned int *Events_FP;          //! papi variable, internal use
  std::vector< std::vector<long_long> > values; //! papi variable, internal use
  double mflops;
  unsigned int papi_event[NUM_EVENTS];   //! papi events for default
  char papi_event_name[NUM_EVENTS][16];  //! name of papi events for print

  long long vsum[NUM_EVENTS];
#endif
  long long start_time;           //! real start time, internal use
  long long stop_time;            //! real stop  time, internal use
  long long start_utime;          //! user start time, internal use
  long long stop_utime;           //! system stop time, internal use
  long long start_stime;          //! user start time, internal use
  long long stop_stime;           //! system stop time, internal use
  std::vector<long long> time;              //! real time for each target
  std::vector<long long> utime;             //! user time for each target
  std::vector<long long> stime;             //! system time for each target

  long long tsum;
  long long usum;
  long long ssum;

  std::vector< std::set<int> > low_target;  //! nested target
  std::vector< std::string > target_name;   //! name of target
  
  int num_target;                //! number of target block/region

  std::vector<int> target;       //! IDs of targets
  int current_target;            //! target now counted

  struct rusage usage;
}

void PerfCounter::init(int nt)
{
#ifdef USE_PAPI
  NumEvents = PAPI_num_counters();
  std::cout << "Number of PAPI Counter " << NumEvents << std::endl;
  Events_FP = (unsigned int *)malloc(NumEvents*sizeof(Events_FP[0]));
  papi_event[0] = PAPI_TOT_CYC;
  papi_event[1] = PAPI_TOT_INS;
  papi_event[2] = PAPI_FP_INS;
  papi_event[3] = PAPI_BR_INS;
  papi_event[4] = PAPI_L2_TCM;
  strcpy(papi_event_name[0],"PAPI_TOT_CYC");
  strcpy(papi_event_name[1],"PAPI_TOT_INS");
  strcpy(papi_event_name[2],"PAPI_FP_INS");
  strcpy(papi_event_name[3],"PAPI_BR_INS");
  strcpy(papi_event_name[4],"PAPI_L2_TCM");

  EventSet = PAPI_NULL;
  for(int i=0;i<std::min(NumEvents,NUM_EVENTS);i++){
    Events_FP[i] = papi_event[i];
  }
  retval = PAPI_library_init(PAPI_VER_CURRENT);
  if(retval != PAPI_VER_CURRENT){
    printf("PAPI library init error!\n");
  }
  values.reserve(nt);
  set_event();
#endif
  time.reserve(nt);
  utime.reserve(nt);
  stime.reserve(nt);
  low_target.reserve(nt);
  target.reserve(nt);
  num_target=0;
  current_target = -1;
}

int
PerfCounter::add_target(std::string name)
{
  int t;
  for(t=0;t<num_target;t++){
    if(name==target_name[t]){
      return t;
    }
  }
  num_target++;
  target_name.push_back(name);
  time.push_back(0);
  utime.push_back(0);
  stime.push_back(0);
  low_target.resize(num_target);
#ifdef USE_PAPI
  std::vector<long_long> zeroval(NumEvents,0);
  values.push_back(zeroval);
#endif
  return t;
}

void PerfCounter::set_event()
{
#ifdef USE_PAPI
  retval = PAPI_create_eventset(&EventSet);
  if(retval!=PAPI_OK)PAPI_perror(retval,(char *)NULL,0);
  retval = PAPI_add_events(EventSet,(int *)Events_FP,NumEvents);
  if(retval!=PAPI_OK)PAPI_perror(retval,(char *)NULL,0);
#endif
}


long long PerfCounter::simple_stop()
{
#ifdef USE_PAPI
  long_long val[NumEvents];
  PAPI_stop(EventSet,val);
  for(int e=0;e<NumEvents;e++){
    values[current_target][e]+=val[e];
  }
#endif
  stop_time = get_realtime_usec(usage);
  stop_utime = timeval_to_ll(usage.ru_utime);
  stop_stime = timeval_to_ll(usage.ru_stime);
  time[current_target] += (stop_time-start_time);
  utime[current_target] += (stop_utime-start_utime);
  stime[current_target] += (stop_stime-start_stime);
  return stop_time;
}

long long PerfCounter::simple_start()
{
  start_time = get_realtime_usec(usage);
  start_utime = timeval_to_ll(usage.ru_utime);
  start_stime = timeval_to_ll(usage.ru_stime);
#ifdef USE_PAPI
  PAPI_start(EventSet);
#endif
  return start_time;
}

long long PerfCounter::start(int t)
{
  if((t<0)||(t>=num_target)){
    return -1;
  }
  if(target.size()>0){
    simple_stop();
    low_target[current_target].insert(t);
  }
  target.push_back(t);
  current_target = t;
  return simple_start();
}

long long PerfCounter::stop()
{
  if(target.size()>0){
    simple_stop();
    target.pop_back();
    if(target.size()>0){
      current_target = target.back();
      simple_start();
    }else{
      current_target = -1;
    }
  }
  return stop_time;
}

void PerfCounter::print(std::ostream& osf)
{
  for(int t=0;t<num_target;t++){
    osf << target_name[t] << " :";
    for(std::set<int>::iterator lt = low_target[t].begin();
        lt!=low_target[t].end(); ++lt){
      osf << " " << target_name[*lt];
    }
    osf << std::endl;
  }

  osf << "name";
  osf << " real";
  osf << " user";
  osf << " sys";
  osf << " user+sys";
#ifdef USE_PAPI
  for(int e=0;e<NumEvents;e++){
    osf << " " << papi_event_name[e];
  }
#endif
  osf << std::endl;


  tsum = 0;
  usum = 0;
  ssum = 0;
#ifdef USE_PAPI
  for(int e=0;e<NumEvents;e++){
    vsum[e] = 0;
  }
#endif
  for(int t=0;t<num_target;t++){
    tsum+=time[t];
    usum+=utime[t];
    ssum+=stime[t];
#ifdef USE_PAPI
    for(int e=0;e<NumEvents;e++){
      vsum[e]+=values[t][e];
    }
#endif
  }
  osf << "total";
  osf << " " << tsum << " us";
  osf << " " << usum << " us";
  osf << " " << ssum << " us";
  osf << " " << usum+ssum << " us";
#ifdef USE_PAPI
  for(int e=0;e<NumEvents;e++){
    osf << " " << vsum[e] ;
  }
#endif
  osf << std::endl;
  for(int t=0;t<num_target;t++){
    osf << target_name[t];
    osf << " " << time[t] << " us";
    osf << " " << utime[t] << " us";
    osf << " " << stime[t] << " us";
    osf << " " << utime[t]+stime[t] << " us";
#ifdef USE_PAPI
    for(int e=0;e<NumEvents;e++){
      osf << " " << values[t][e] ;
    }
#endif
    osf << std::endl;
  }
}

void PerfCounter::print_time(std::ostream& osf)
{
  osf << "name";
  osf << " real";
  osf << std::endl;


  tsum = 0;
  for(int t=0;t<num_target;t++){
    tsum+=time[t];
  }
  osf << "total";
  osf << " " << tsum << " us";
  osf << std::endl;
  for(int t=0;t<num_target;t++){
    osf << target_name[t];
    osf << " " << time[t] << " us";
    osf << std::endl;
  }
}

void PerfCounter::print(int& rank, std::string namebase)
{
  std::stringstream filename;
  filename << namebase << "." << rank;
  std::string ofsname;
  filename >> ofsname;
  std::ofstream ofs(ofsname.c_str());

  PerfCounter::print(ofs);
}

double * PerfCounter::get_time_array(int& num)
{
  num = time.size();
  double *ta = new double[num];
  for(int i=0;i<num;i++){
    ta[i] = time[i];
  }
  return ta;
}

void PerfCounter::put_time_array(double *times, int &num)
{
  int n = num;
  if(num>static_cast<int>(time.size())){
    n = static_cast<int>(time.size());
  }
  for(int i=0;i<n;i++){
    time[i] = static_cast<long long>(times[i]);
  }
}
