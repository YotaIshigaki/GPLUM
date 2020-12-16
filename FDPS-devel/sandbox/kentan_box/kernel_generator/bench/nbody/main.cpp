#include <cstdio>
#include <cstdlib>
#include <assert.h>

#include <particle_simulator.hpp>
#include "user_defined_class.h"

using SPJ = SPJQuad_grav;

#include <unistd.h>
#include <time.h>
double get_dtime(void){
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC,&tp);
  return ((double)(tp.tv_sec)+(double)(tp.tv_nsec)*1e-9);
}

#include "kernel.cpp"
#include "kernel_serial.c"

double random_number(){
  return (double)rand()/((double)RAND_MAX+1);
}

void makeColdUniformSphere(const double mass_glb,
                           const int n_loc,
                           EPJ* epj,
                           const double eng,
                           const int seed) {
  assert(eng < 0.0);
  {
    for(int i = 0; i < n_loc; i++){
      epj[i].mass = mass_glb / n_loc;
      epj[i].eps2 = 0.01f;
      const double radius = 3.0;
      do {
	epj[i].pos.x = (2. * random_number() - 1.) * radius;
	epj[i].pos.y = (2. * random_number() - 1.) * radius;
	epj[i].pos.z = (2. * random_number() - 1.) * radius;
      }while(epj[i].pos*epj[i].pos >= radius * radius);
      //printf("%f %f %f\n",epj[i].pos.x,epj[i].pos.y,epj[i].pos.z);
    }
  }

  double cm_x = 0.0;
  double cm_y = 0.0;
  double cm_z = 0.0;
  double cm_m = 0.0;
  for(int i = 0; i < n_loc; i++){
    cm_x  += epj[i].pos.x * epj[i].mass;
    cm_y  += epj[i].pos.y * epj[i].mass;
    cm_z  += epj[i].pos.z * epj[i].mass;
    cm_m  += epj[i].mass;
  }
  cm_x /= cm_m;
  cm_y /= cm_m;
  cm_z /= cm_m;
  for(int i = 0; i < n_loc; i++){
    epj[i].pos.x -= cm_x;
    epj[i].pos.y -= cm_y;
    epj[i].pos.z -= cm_z;
  }
}

//const int Niter = 10240;
const int Niter = 1024;
const int N = 512;
int main(){
  assert(NJ_EPEP<=N && NJ_EPSP <= N);
  EPI *epi;
  EPJ *epj;
  SPJ *spj;
  Force *force;
  posix_memalign((void**)(&epi), 512, sizeof(EPI)*N);
  posix_memalign((void**)(&epj), 512, sizeof(EPJ)*N);
  posix_memalign((void**)(&spj), 512, sizeof(SPJ)*N);
  posix_memalign((void**)(&force), 512, sizeof(Force)*N);

#if defined(CHECK_FORCE) || defined(CHECK_FORCE_DETAIL)
  makeColdUniformSphere(1.0,N,epj,-0.25,0);
  for(int i=0;i<N;i++){
    //printf("%f %f %f %f %f\n",epj[i].pos.x,epj[i].pos.y,epj[i].pos.z,epj[i].mass,epj[i].eps2);
    epi[i].pos  = epj[i].pos;
    epi[i].eps2 = 0.01f;

    spj[i].pos  = epj[i].pos;
    spj[i].mass = epj[i].mass;
    spj[i].eps2 = 0.01f;
    spj[i].xx = random_number();
    spj[i].yy = random_number();
    spj[i].zz = random_number();
    spj[i].xy = random_number();
    spj[i].yz = random_number();
    spj[i].zx = random_number();
  }
  for(int i=0;i<N;i++){
    force[i].acc.x = force[i].acc.y = force[i].acc.z = force[i].pot = 0.f;
  }
#endif

  for(int i=0;i<100;i++) CalcForceLongEpEp()(epi,N,epj,N,force);
  //for(int i=0;i<100;i++) CalcForceLongEpSp()(epi,N,spj,N,force);

  printf("N=%d, NJ_EPEP=%d, NJ_EPSP=%d\n",N,NJ_EPEP,NJ_EPSP);

  double t = get_dtime();
  for(int i=0;i<Niter;i++) CalcForceLongEpEp()(epi,N,epj,N,force);
  double elapsed_time_epep = get_dtime() - t;
  printf("epep_elapsed_time: %lf sec, %lf Gflops\n",elapsed_time_epep,((unsigned long long)Niter*27*(unsigned long long)N*(unsigned long long)N)/elapsed_time_epep/1e9);

  t = get_dtime();
  for(int i=0;i<Niter;i++) CalcForceLongEpSp()(epi,N,spj,N,force);
  double elapsed_time_epsp = get_dtime() - t;
  printf("epsp_elapsed_time: %f sec, %lf Gflops\n",elapsed_time_epsp,((unsigned long long)Niter*66*(unsigned long long)(N*N))/elapsed_time_epsp/1e9);

#if defined(CHECK_FORCE)||defined(CHECK_FORCE_DETAIL)
  for(int i=0;i<N;i++){
    force[i].acc.x = force[i].acc.y = force[i].acc.z = force[i].pot = 0.f;
  }
  CalcForceLongEpEp()(epi,N,epj,N,force);
  CalcForceLongEpSp()(epi,N,spj,N,force);
  Force* force_orig = (Force*)malloc(sizeof(Force)*N);
  for(int i=0;i<N;i++){
    force_orig[i].acc.x = force_orig[i].acc.y = force_orig[i].acc.z = force_orig[i].pot = 0.f;
  }
  CalcForceLongEpEpOrig(epi,N,epj,N,force_orig);
  CalcForceLongEpSpOrig(epi,N,spj,N,force_orig);

  Force ave,max;
  ave.acc.x  =  ave.acc.y = ave.acc.z = 0.f;
  max.acc.x  =  max.acc.y = max.acc.z = 0.f;
  for(int i=0;i<N;i++){
    float dx = fabs(force[i].acc.x - force_orig[i].acc.x);
    float dy = fabs(force[i].acc.y - force_orig[i].acc.y);
    float dz = fabs(force[i].acc.z - force_orig[i].acc.z);
#ifdef CHECK_FORCE_DETAIL
    printf("%4d %e %e %e, %e %e %e, %e %e %e\n",
	   i,force[i].acc.x,force[i].acc.y,force[i].acc.z,
	   force_orig[i].acc.x,force_orig[i].acc.y,force_orig[i].acc.z,
	   dx,dy,dz);
#endif
    ave.acc.x += dx;
    ave.acc.y += dy;
    ave.acc.z += dz;
    max.acc.x = std::max(dx,ave.acc.x);
    max.acc.y = std::max(dy,ave.acc.y);
    max.acc.z = std::max(dz,ave.acc.z);
  }
  printf("average absolute error of force: %e %e %e\n",ave.acc.x/N,ave.acc.y/N,ave.acc.z/N);
  printf("maximum absolute error of force: %e %e %e\n",max.acc.x,max.acc.y,max.acc.z);
  free(force_orig);
#endif

  free(force);
  free(spj);
  free(epj);
  free(epi);

  return 0;
}

