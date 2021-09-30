#include "Timer.h"


int tllll;
int tll;
int tlp;

inline
double ll(double a[], double b[], double c[], int n, int r)
{
  PerfCounter::start(tll);
  for(int j=0;j<r;j++){
    for(int i=0;i<n;i++){
      a[i] += b[i]*c[i];
    }
  }
  PerfCounter::stop();
}

inline
double lp(double a[], double b[], double c[], int n, int r)
{
  PerfCounter::start(tlp);
  for(int j=0;j<r;j++){
    for(int i=0;i<n;i++){
      a[i] += b[i];
    }
  }
  PerfCounter::stop();
}


inline
double llll(int l, int r)
{
  int n = l;
  double a[n];
  double b[n];
  double c[n];

  PerfCounter::start(tllll);
  for(int j=0;j<r;j++){
    for(int i=1;i<n;i++){
      a[0] += a[i];
    }
  }
  ll(a,b,c,n,r);
  lp(a,b,c,n,r);
  for(int j=0;j<r;j++){
    for(int i=1;i<n;i++){
      a[0] += a[i];
    }
  }
  PerfCounter::stop();
  return a[0];
}

int
main(int argc, char **argv)
{
  int l=100;
  int r=100;
  int r2=1;

  if(argc>1){
    r = atol(argv[1]);
    if(argc>2){
      r2 = atol(argv[2]);
    } 
  }

  int t=PerfCounter::add_target(std::string("main"));
  tllll=PerfCounter::add_target(std::string("llll"));
  tll=PerfCounter::add_target(std::string("ll"));
  tlp=PerfCounter::add_target(std::string("lp"));

  PerfCounter::start(t);
  
  for(int k=0;k<r2;k++){
    llll(l,r);
  }

  PerfCounter::stop();

  PerfCounter::print();

}
