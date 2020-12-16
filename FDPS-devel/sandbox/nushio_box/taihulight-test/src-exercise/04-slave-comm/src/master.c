#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "param.h"

extern SLAVE_FUN(move_r)();
extern SLAVE_FUN(move_c)();

double a[N];
double b[N];
double c[N];

int main() {
  int i;
  athread_init();

  for (i=0; i<N;++i){
    a[i] = i/4;
    b[i] = i/4;
  }

  printf("before:");
  for (i=0; i<N; i++){
    if(i%4==0) printf(" ");
    if(i%32==0) printf("\n");
    printf("%2.0lf", a[i]);
  }
  printf("\n");

  athread_spawn(move_r,0);
  athread_join();

  printf("after move_r:");
  for (i=0; i<N; i++){
    if(i%4==0) printf(" ");
    if(i%32==0) printf("\n");
    printf("%2.0lf", c[i]);
  }
  printf("\n");

  athread_spawn(move_c,0);
  athread_join();

  printf("after move_c:");
  for (i=0; i<N; i++){
    if(i%4==0) printf(" ");
    if(i%32==0) printf("\n");
    printf("%2.0lf", c[i]);
  }
  printf("\n");

  athread_halt();
  return 0;
}
