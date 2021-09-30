#include <slave.h>
#include <math.h>
#include <dma.h>

#define N 4096
#define I 64

__thread_local volatile unsigned long get_reply, put_reply;
__thread_local int my_id;
__thread_local double a_slave[I], b_slave[I], c_slave[I];
extern double a[N], b[N], c[N];

void func() {
  int i;
  my_id = athread_get_id(-1);
  get_reply = 0;

  athread_get(PE_MODE, &a[my_id*I], &a_slave[0],I*8,&get_reply,0,0,0);
  athread_get(PE_MODE, &b[my_id*I], &b_slave[0],I*8,&get_reply,0,0,0);
  while(get_reply!=2) {}

  for(i=0;i<I;i++){
    c_slave[i]=a_slave[i]*b_slave[i];
  }

  put_reply=0;
  athread_put(PE_MODE,&c_slave[0],&c[my_id * I],I*8,&put_reply,0,0);
  while(put_reply!=1) {}

}
