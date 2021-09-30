#include <mpi.h>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
#include "Common.h"

static inline Position pbc_minimum_image(const Position &inp){
  return Position(
                  inp.x - round(inp.x),
                  inp.y - round(inp.y),
                  inp.z - round(inp.z));
}

template<typename PA>
static void PP_interact_inner(PA &pa, const std::vector<TypeRange> &ta, const int ti, const int tj, ForceArray &fa, std::vector<double> &phi_direct){
  //  int n = 0;
  for(int i=ta[ti].begin; i<ta[ti].end; i++){
    const Position pi = getpos(pa,i);
    const double chagei = getcharge(pa,i);
    for(int j=ta[tj].begin; j<ta[tj].end; j++){
      const Position pj = getpos(pa,j);
      const double chargej = getcharge(pa,j);
      if(i == j) continue;
      //      printf("PP pair %d %d\n",i,j);
      //      n++;
      const Position  dr  = pbc_minimum_image(pj - pi);
      const double r2  = dr*dr;
      const double ri2 = 1.0 / r2;
      const double ri  = sqrt(ri2);
      const double ri3 = ri * ri2;
      phi_direct[i] += chargej * ri;
      fa[i] -= (chargej * ri3) * dr;
      // puts("evaluated PP");
    }
  }
  //  printf("cell pair %d %d : PP pair %d\n",ti,tj,n);
}

template<typename PA>
void PP_interact_PBC(PA &pa, const std::vector<TypeRange> &ta,
		     ForceArray &fa, std::vector<double> &phi_direct,
		     const int ICUT, const int NX, const int NY, const int NZ)
{
  if(DebugLog::verbose>1)printf("ta.size() %d\n",ta.size());
  if(DebugLog::verbose>1)printf("PP %d %d %d %d\n",ICUT,NX,NY,NZ);
  int i, j, k, ti;
  int nt = NZ*NY*NX;

#pragma omp parallel for
  for(ti=0;ti<nt;ti++){
    k = ti/(NX*NY);
    j = (ti-k*NX*NY)/NX;
    i = (ti-k*NX*NY-j*NY);
    /*
  for(k=0; k<NZ; k++)
    for(j=0; j<NY; j++)
      for(i=0; i<NX; i++)
	{
	  int ti = i+NX*(j+NY*k);
    */
	  int ii, jj, kk;
	  for(kk=k-ICUT; kk<=k+ICUT; kk++)
	    for(jj=j-ICUT; jj<=j+ICUT; jj++)
	      for(ii=i-ICUT; ii<=i+ICUT; ii++)
		{
		  const int iii = (ii+NX)%NX;
		  const int jjj = (jj+NY)%NY;
		  const int kkk = (kk+NX)%NZ;
		  int tj = iii+NX*(jjj+NY*kkk);
		  PP_interact_inner(pa,ta,ti,tj,fa,phi_direct);
		}
	}
}


