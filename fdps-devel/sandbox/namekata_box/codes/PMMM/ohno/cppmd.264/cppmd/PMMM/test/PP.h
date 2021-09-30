#include <mpi.h>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
#include "vector3.h"

#include "cell.h"

static void PP_interact_inner(Cell &ci, const Cell &cj){
	const int ni = ci.plist.size();
	const int nj = cj.plist.size();
	for(int i=0; i<ni; i++){
		Particle &pi = *ci.plist[i];
		for(int j=0; j<nj; j++){
			const Particle &pj = *cj.plist[j];
			if(&pi == &pj) continue;

			const dvec3  dr  = minimum_image(pj.pos - pi.pos);
			const double r2  = dr*dr;
			const double ri2 = 1.0 / r2;
			const double ri  = sqrt(ri2);
			const double ri3 = ri * ri2;
			pi.phi_direct += pj.mass * ri;
			pi.acc_direct += (pj.mass * ri3) * dr;
			// puts("evaluated PP");
		}
	}
}

template <int PFMM, int ICUT, int NX, int NY, int NZ>
void PP_interact_PBC(Cell_FMM<PFMM> cell[NZ][NY][NX])
{
	int i, j, k;
	for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
	{
		int ii, jj, kk;
		for(kk=k-ICUT; kk<=k+ICUT; kk++) for(jj=j-ICUT; jj<=j+ICUT; jj++) for(ii=i-ICUT; ii<=i+ICUT; ii++)
		{
			const int iii = (ii+NX)%NX;
			const int jjj = (jj+NY)%NY;
			const int kkk = (kk+NX)%NZ;
			PP_interact_inner(cell[k][j][i], cell[kkk][jjj][iii]);
		}
	}
}


