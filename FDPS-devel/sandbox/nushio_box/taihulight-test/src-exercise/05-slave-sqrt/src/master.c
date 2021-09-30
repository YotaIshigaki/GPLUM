#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "param.h"

extern SLAVE_FUN(func)();

double a[N];
double b[N];
double c[N];

int main() {
  int i;
  printf("hello Sunway TaihuLight\n");

  for (i=0; i<N;++i){
    a[i] = i;
    b[i] = i;
  }

  // for (i=0; i<N;++i){
  //   c[i] = a[i] * b[i];
  // }
  athread_init();
  athread_spawn(func,0);//fflush(NULL);
  athread_join();

  for (i=1; i<N; i=2*i+1){
    printf("%d^2 == %lf\n", i, c[i]);
  }
  athread_halt();
  return 0;
}
