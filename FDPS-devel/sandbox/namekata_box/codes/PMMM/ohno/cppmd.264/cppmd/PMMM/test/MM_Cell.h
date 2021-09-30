#pragma once
#include <mpi.h>

#include <cstdio>

#include "M2L_convolution_PBC.h"
#include "MPIParallelPMMM.h"

template<int p, int NX, int NY, int NZ, int NMAX, int MMAX, int ICUT>
struct MM_Cell{
  GreenFunction_PBC<p, NX, NY, NZ, NMAX, MMAX, ICUT> gf;
  MultipoleMoment<p> mm[NZ][NY][NX];
  LocalExpansion <p> le[NZ][NY][NX];
  M2L_convolution_PBC<p,NX,NY,NZ,NMAX,MMAX,ICUT> m2l_conv_pbc;
  MM_MM_Comm<p,NX,NY,NZ> mm_mm_com;

  void initialize(const double alpha, const double cell_length){
    gf.gen_gf_r(alpha, cell_length);
    gf.gen_gf_k();
  }

  void set_mm_mm_com(const MPI::Intracomm &_comm, const std::vector<int> &_mm_list, const int _mm_targets[(p+1)*(p+1)])
  {
    mm_mm_com.set_comm(_comm);
    mm_mm_com.set_mm_list(_mm_list);
    mm_mm_com.set_mm_count(_mm_targets);
    m2l_conv_pbc.set_mm_list(_mm_list);
  }

  void set_mm(const MultipoleMoment<p> m, const int x, const int y, const int z){
    mm[z][y][x] = m;
  }

#if 0
  void set_mm(const MultipoleMoment<p> m[NZ][NY][NX]){
    int i,j,k;
    for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++) mm[k][j][i] = m[k][j][i];
  }
#endif

  LocalExpansion<p> get_le(const int x, const int y, const int z){
    return le[z][y][x];
  }

#if 0
  void get_le(LocalExpansion <p> l[NZ][NY][NX]){
    int i,j,k;
    for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++) l[NZ][NY][NX] = le[NZ][NY][NX];
  }
#endif

  void convolution_PBC(){
    puts("Forward");
    m2l_conv_pbc.forward(mm);
    if(mm_mm_com.mm_size>1){
      puts("Exchange");
      mm_mm_com.exchange_mm(m2l_conv_pbc.mm_k);
    }
    puts("Transform");
    m2l_conv_pbc.transform(gf);
    puts("Backward");
    m2l_conv_pbc.backward(le);
  }
};
