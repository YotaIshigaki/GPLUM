#include <stdio.h>
#include <athread.h>
#include <fcntl.h>
#define N 4096


//extern SLAVE_FUN(func)();

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

  for (i=0; i<N;++i){
    c[i] = a[i] * b[i];
  }

  for (i=1; i<N; i=2*i+1){
    printf("%d^2 == %lf\n", i, c[i]);
  }

  return 0;
}
