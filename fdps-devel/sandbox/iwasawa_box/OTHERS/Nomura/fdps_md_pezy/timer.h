#ifndef TIMER_H
#define TIMER_H
//copied from HERMITE by Dr. Nitadori

#include <sys/time.h>
struct Profile{
  enum{
    PEZY = 0,
    CPU,
    INTEG,
    DECOMP,
    EXCHANGE,
    OTHER,
    KERNEL,
    TREE,
    TRANSFER,
    TRANSLATE,
    // don't touch the below
    TOTAL,
    MISC,
    NUM_ELEM,
  };
  static const char *name(const int i){
    static const char *strs[NUM_ELEM] = {
      "PEZY      ",
      "CPU       ",
      "Integ     ",
      "Decomp    ",
      "Exchange  ",
      "Other     ",
      "KERNEL    ",
      "Tree      ",
      "Transfer  ",
      "Translate ",
      "Total     ",
      "Misc      ",
    };
    return strs[i];
  }
  static double wtime(){
#if 0
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
#else
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
    Comm::barrier();
#endif //PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
    return MPI::Wtime();
#elif PARTICLE_SIMULATOR_THREAD_PARALLEL
    return omp_get_wtime();
#else
    return clock() / CLOCKS_PER_SEC;
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
#endif
  }
  double time[NUM_ELEM];
  double tprev;
  
  long num_step_prev;
  long num_bstep_prev;
  
  void flush(){
    for(int i=0; i<NUM_ELEM; i++) time[i] = 0.0;
  }
#ifdef TUNING
  void beg(const int elem, const bool reuse = false){
    if(reuse) time[elem] -= tprev;
    else      time[elem] -= wtime();
  }
  void end(const int elem){
    // time[elem] += wtime();
    tprev = wtime();
    time[elem] += tprev;
  }
  void show(){
    time[TREE] = time[PEZY] - (time[KERNEL]+time[TRANSFER]+time[TRANSLATE]);
    for(int i=0;i<NUM_ELEM;i++)
      printf("%s%e\n",name(i),time[i]);
  }
  #else
  void beg(const int elem, const bool reuse = false){};
  void end(const int elem){};
  void show(){}
  #endif
};

#endif
