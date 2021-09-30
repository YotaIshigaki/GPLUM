#include <stdio.h>
#include "param.h"

extern Real sU0[SX][SY][SZ], sV0[SX][SY][SZ];
extern Real Uwx[T_MAX][2][SY][SZ], Uwy[T_MAX][SX][2][SZ], Uwz[T_MAX][SX][SY][2];
extern Real Vwx[T_MAX][2][SY][SZ], Vwy[T_MAX][SX][2][SZ], Vwz[T_MAX][SX][SY][2];


// double-buffered simulation state
Real sU[SX][SY][SZ], sV[SX][SY][SZ];
Real sU_1[SX][SY][SZ], sV_1[SX][SY][SZ];


void run_benchmark(){

  printf("Carrying out simulation...\n");

  // set initial condition
  for(int x=0;x<SX;++x) {
    for(int y=0;y<SY;++y) {
      for(int z=0;z<SZ;++z) {
        sU[x][y][z]=sU0[x][y][z];
        sV[x][y][z]=sV0[x][y][z];
      }
    }
  }

  for(int t = 0; t < T_MAX; ++t){
    // load communication values
    for(int x=SX-2;x<SX;++x) {
      for(int y=0;y<SY;++y) {
        for(int z=0;z<SZ;++z) {
          sU[x][y][z] = Uwx[t][x-(SX-2)][y][z];
          sV[x][y][z] = Vwx[t][x-(SX-2)][y][z];
        }
      }
    }

    for(int x=0;x<SX-2;++x) {
      for(int y=SY-2;y<SY;++y) {
        for(int z=0;z<SZ;++z) {
          sU[x][y][z] = Uwy[t][x][y-(SY-2)][z];
          sV[x][y][z] = Vwy[t][x][y-(SY-2)][z];
        }
      }
    }

    for(int x=0;x<SX-2;++x) {
      for(int y=0;y<SY-2;++y) {
        for(int z=SZ-2;z<SZ;++z) {
          sU[x][y][z] = Uwz[t][x][y][z-(SZ-2)];
          sV[x][y][z] = Vwz[t][x][y][z-(SZ-2)];
        }
      }
    }



    // destructively update the state
#define lap(ar, x, y, z)                        \
    ( ar[x][y+1][z+1] + ar[x+2][y+1][z+1]       \
      + ar[x+1][y][z+1] + ar[x+1][y+2][z+1]     \
      + ar[x+1][y+1][z] + ar[x+1][y+1][z+2]     \
        - 6*ar[x+1][y+1][z+1]) / dx / dx

    for(int x=0;x<SX-2;++x) {
      for(int y=0;y<SY-2;++y) {
        for(int z=0;z<SZ-2;++z) {
          Real u=sU[x+1][y+1][z+1] ;
          Real v=sV[x+1][y+1][z+1] ;

          Real du_dt = -Fe * u*v*v + Fu*(1-u) + Du * lap(sU,x,y,z);
          Real dv_dt =  Fe * u*v*v - Fv*v     + Dv * lap(sV,x,y,z);
          sU_1[x][y][z] = u+dt*du_dt;
          sV_1[x][y][z] = v+dt*dv_dt;
        }
      }
    }


    for(int x=0;x<SX-2;++x) {
      for(int y=0;y<SY-2;++y) {
        for(int z=0;z<SZ-2;++z) {
          sU[x][y][z] = sU_1[x][y][z];
          sV[x][y][z] = sV_1[x][y][z];
        }
      }
    }
  }

  // return the final condition
  for(int x=0;x<SX;++x) {
    for(int y=0;y<SY;++y) {
      for(int z=0;z<SZ;++z) {
        sU0[x][y][z]=sU[x][y][z];
        sV0[x][y][z]=sV[x][y][z];
      }
    }
  }

}
