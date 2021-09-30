#ifndef TIMER_H
#define TIMER_H
//partially copied from HERMITE by Dr. Nitadori
#include <sys/time.h>

inline unsigned long long rdtsc() {
    unsigned long long ret;
    __asm volatile ("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax" : "=a" (ret) :: "%rdx");
    //__asm__ volatile ("rdtsc" : "=A" (ret));
    //fprintf(stderr,"%llu\n",ret);
    return ret;
}

#define TUNING
struct Profile{
  enum{
    FORCE = 0,
    INTEG,
    DECOMP,
    EXCHANGE,
    OTHER,
    // don't touch the below
    TOTAL,
    MISC,
    NUM_ELEM,
  };
  static const char *name(const int i){
    static const char *strs[NUM_ELEM] = {
      "FORCE     ",
      "INTEG     ",
      "DECOMP    ",
      "EXCHANGE  ",
      "Other     ",
      "Total     ",
      "Misc      ",
    };
    return strs[i];
  }
  static double wtime(){
#if 1
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
#else
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
    PS::Comm::barrier();
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
  double tprev[NUM_ELEM];

  unsigned long long clk[NUM_ELEM];
  unsigned long long clkprev[NUM_ELEM];

  long num_step_prev;
  long num_bstep_prev;

  void flush(){
    for(int i=0; i<NUM_ELEM; i++){
      time[i] = 0.0;
      clk[i] = 0;
    }
  }
#ifdef TUNING
  void beg(const int elem, const bool reuse = false){
    __asm("//Timer beg");
#ifdef TIMER_COUNT_BY_CLOCK
    clk[elem]  -= rdtsc();
#else
    if(reuse){
      time[elem] -= tprev[elem];
    }else{
      time[elem] -= wtime();
    }
#endif
  }
  void end(const int elem){
    __asm("//Timer end");
#ifdef TIMER_COUNT_BY_CLOCK
    clkprev[elem] = rdtsc();
    clk[elem]  += clkprev[elem];
#else
    tprev[elem] = wtime();
    time[elem] += tprev[elem];
#endif
  }
  void show(){
    clk[TOTAL] = 0;
    for(int i=0;i<TOTAL;i++){
      clk[TOTAL]  += clk[i];
    }
    for(int i=0;i<NUM_ELEM;i++)
#ifdef TIMER_COUNT_BY_CLOCK
      fprintf(stderr,"%s\t%llu\t%lf\n",
	      name(i),clk[i],(double)clk[i]/(double)clk[TOTAL]*100);
#else
  #if 0
      fprintf(stderr,"%s\t%e\t%5.2f\n",
  	      name(i),time[i],time[i]/time[TOTAL]*100);
  #else
      fprintf(stderr," \t%s",name(i));
    fprintf(stderr,"\n");
    for(int i=0;i<NUM_ELEM;i++)
      fprintf(stderr," \t%e",time[i]);
    fprintf(stderr,"\n");
  #endif
#endif
  }
  #else
  void beg(const int elem, const bool reuse = false){};
  void end(const int elem){};
  void show(){}
  #endif

  Profile(){flush();}
  ~Profile(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    if(PS::Comm::getRank()==0)
#endif
      show();

  }
};

#endif
