#include <slave.h>
#include <math.h>
#include <dma.h>

#include "param.h"

__thread_local volatile unsigned long get_reply, put_reply;
__thread_local doublev4 a_slave[Iv4], b_slave[Iv4], c_slave[Iv4];
extern double a[N], b[N], c[N];

void move_r() {
  int i;
  int my_id = athread_get_id(-1);
  int cid = my_id%8, rid = my_id/8;

  get_reply = 0;

  athread_get(PE_MODE, &a[my_id*I], &a_slave[0],I*8,&get_reply,0,0,0);
  athread_get(PE_MODE, &b[my_id*I], &b_slave[0],I*8,&get_reply,0,0,0);
  while(get_reply!=2) {}

  simd_putr(a_slave[0], (cid+1)%8);
  c_slave[0] = simd_getr(c_slave[0]);

  put_reply=0;
  athread_put(PE_MODE,&c_slave[0],&c[my_id * I],I*8,&put_reply,0,0);
  while(put_reply!=1) {}

}

void move_c() {
  int i;
  int my_id = athread_get_id(-1);
  int cid = my_id%8, rid = my_id/8;

  get_reply = 0;

  athread_get(PE_MODE, &a[my_id*I], &a_slave[0],I*8,&get_reply,0,0,0);
  athread_get(PE_MODE, &b[my_id*I], &b_slave[0],I*8,&get_reply,0,0,0);
  while(get_reply!=2) {}

  simd_putc(b_slave[0], (rid*3)%8);
  c_slave[0] = simd_getc(c_slave[0]);

  put_reply=0;
  athread_put(PE_MODE,&c_slave[0],&c[my_id * I],I*8,&put_reply,0,0);
  while(put_reply!=1) {}

}
