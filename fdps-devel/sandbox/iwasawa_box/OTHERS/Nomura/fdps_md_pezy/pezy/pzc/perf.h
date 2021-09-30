#pragma once
#include"../../common/pzcperf.h"
#define MAX_PERF  4

PZCPerformance gPerf[MAX_PERF];

void pzc_GetPerformance(PZCPerformance*	out,int idx)
{
  int tid = get_tid(); // thread ID (0 - 7)
  int pid = get_pid(); // PE id (depends on work_size)
  if(tid ==0 && pid==0 && idx < MAX_PERF){
    *out = gPerf[ idx ];
    gPerf[ idx ].Clear();
  }
  flush();
}
