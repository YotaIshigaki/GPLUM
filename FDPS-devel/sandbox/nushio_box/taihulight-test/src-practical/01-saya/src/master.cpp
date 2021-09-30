#include <cmath>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h>

#include "param.h"



Real U[NX][NY][NZ], V[NX][NY][NZ];
Real U_other[NX][NY][NZ], V_other[NX][NY][NZ];
int global_clock;


Real Uwx[T_MAX][2][SY][SZ], Uwy[T_MAX][SX][2][SZ], Uwz[T_MAX][SX][SY][2];
Real Vwx[T_MAX][2][SY][SZ], Vwy[T_MAX][SX][2][SZ], Vwz[T_MAX][SX][SY][2];

Real sU0[SX][SY][SZ], sV0[SX][SY][SZ];


extern "C" {
  void run_benchmark();
}

double wctime() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec*1e-6;
}



void fill_initial_condition() {
  global_clock=0;
  for (int x=0;x<NX;++x) {
    for (int y=0;y<NY;++y) {
      for (int z=0;z<NZ;++z) {
        U[x][y][z] = 1;
        V[x][y][z] = 0;
      }
    }
  }
  int bx = std::max(NX/4,NX/2-8),  ex = std::min(3*NX/4+1,NX/2+8);
  int by = std::max(NY/4,NY/2-8),  ey = std::min(3*NY/4+1,NY/2+8);
  int bz = std::max(NZ/4,NZ/2-8),  ez = std::min(3*NZ/4+1,NZ/2+8);
  for (int x=bx;x<ex;++x){
    for (int y=by;y<ey;++y){
      for (int z=bz;z<ez;++z){
        U[x][y][z] = 0.5;
        V[x][y][z] = 0.25+0.1*sin(x+sqrt(y)+cos(z));
      }
    }
  }
}


inline Real periodic(Real ar[NX][NY][NZ],int x, int y, int z) {
  x = ((x+100*NX)%NX+NX)%NX;
  y = ((y+100*NY)%NY+NY)%NY;
  z = ((z+100*NZ)%NZ+NZ)%NZ;
  //x = (x+NX)%NX;
  //y = (y+NY)%NY;
  //z = (z+NZ)%NZ;
  return ar[x][y][z];
}


void naive_proceed() {
  ++global_clock;

  auto lap = [](Real ar[NX][NY][NZ],int x, int y, int z) {
    auto ret = periodic(ar, x-1, y, z) + periodic(ar, x+1, y, z)
    + periodic(ar, x, y-1, z) + periodic(ar, x, y+1, z)
    + periodic(ar, x, y, z-1) + periodic(ar, x, y, z+1)
    - 6*ar[x][y][z];
    return ret / dx / dx;
  };

  for (int x=0;x<NX;++x) {
    for (int y=0;y<NY;++y) {
      for (int z=0;z<NZ;++z) {
        auto u = U[x][y][z],  v = V[x][y][z];
        auto du_dt = -Fe * u*v*v + Fu*(1-u) + Du * lap(U,x,y,z);
        auto dv_dt =  Fe * u*v*v - Fv*v     + Dv * lap(V,x,y,z);
        U_other[x][y][z] = U[x][y][z] + dt*du_dt;
        V_other[x][y][z] = V[x][y][z] + dt*dv_dt;
      }
    }
  }

  for (int x=0;x<NX;++x) {
    for (int y=0;y<NY;++y) {
      for (int z=0;z<NZ;++z) {
        U[x][y][z]=U_other[x][y][z];
      }
    }
  }
  for (int x=0;x<NX;++x) {
    for (int y=0;y<NY;++y) {
      for (int z=0;z<NZ;++z) {
        V[x][y][z]=V_other[x][y][z];
      }
    }
  }
}

void get_solution_at(int t, int x, int y, int z, Real &u, Real &v) {
  if(global_clock > t) fill_initial_condition();
  while(global_clock < t) naive_proceed();
  u = periodic(U,x,y,z);
  v = periodic(V,x,y,z);
}

int main () {

  fill_initial_condition();
  for(int x=0;x<SX;++x) {
    for(int y=0;y<SY;++y) {
      for(int z=0;z<SZ;++z) {
        double u,v; get_solution_at(0,x,y,z, u,v);
        sU0[x][y][z]=u;
        sV0[x][y][z]=v;
      }
    }
  }

  std::cerr << "Setting up wall values..." << std::endl;
  for(int t = 0;t<T_MAX;++t){
    for(int x=SX-2;x<SX;++x) {
      for(int y=0;y<SY;++y) {
        for(int z=0;z<SZ;++z) {
          double u,v; get_solution_at(t,x+t,y+t,z+t, u,v);
          Uwx[t][x-(SX-2)][y][z] = u;
          Vwx[t][x-(SX-2)][y][z] = v;
        }
      }
    }

    for(int x=0;x<SX;++x) {
      for(int y=SY-2;y<SY;++y) {
        for(int z=0;z<SZ;++z) {
          double u,v; get_solution_at(t,x+t,y+t,z+t, u,v);
          Uwy[t][x][y-(SY-2)][z] = u;
          Vwy[t][x][y-(SY-2)][z] = v;
        }
      }
    }

    for(int x=0;x<SX;++x) {
      for(int y=0;y<SY;++y) {
        for(int z=SZ-2;z<SZ;++z) {
          double u,v; get_solution_at(t,x+t,y+t,z+t, u,v);
          Uwz[t][x][y][z-(SZ-2)] = u;
          Vwz[t][x][y][z-(SZ-2)] = v;
        }
      }
    }
  }


  for(int trial=0;trial<10;++trial) {

    double time_begin = wctime();

    run_benchmark();

    double time_end = wctime();

    double flop = 29.0 * (SX-2)*(SY-2)*(SZ-2) *T_MAX;
    double time_elapse = time_end-time_begin;

    {
      const int t = T_MAX;
      double num=0,den=0;
      for(int x=0;x<SX-2;++x) {
        for(int y=0;y<SY-2;++y) {
          for(int z=0;z<SZ-2;++z) {
            double u,v; get_solution_at(t,x+t,y+t,z+t, u,v);
            num += std::abs(u-sU0[x][y][z]);
            den += 1;
          }
        }
      }
      std::ostringstream msg;
      msg << SX << " " << SY << " " << SZ << " " << T_MAX << " "
          << " t: " << time_elapse << " GFlops: " << flop/time_elapse/1e9<< " error: " << (num/den);
      std::ofstream log_file("benchmark.txt", std::ios::app);
      std::cout << msg.str() << std::endl;
      log_file << msg.str() << std::endl;
    }
  }
}
