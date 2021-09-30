#include <particle_simulator.hpp>
#include <water_params.h>

#include "user_defined_class.h"
#include "cuda_pointer.h"
#include "kernel.h"

enum{
  N_WALK_LIMIT = 1200,
  NI_LIMIT     = N_WALK_LIMIT * 1000,
  NJ_LIMIT     = N_WALK_LIMIT * 10000,
  WARP_SIZE = 32,
};

class EpiDev{
 public:
  float4 rm;
  int    w;
  int   id; // type
};
class EpjDev{
 public:
  float4 rm;
  int    id; // type
};
class ForceDev{
 public:
  float4 fcl;
  float4 flj;
  __device__
  void init(){
    fcl.x = fcl.y = fcl.z = fcl.w = 0.0;
    flj.x = flj.y = flj.z = flj.w = 0.0;
  }
};

inline __device__
ForceDev accumulate_lips_sw
(const float4 rmi, const int idi, const int ti,
 const float4 rmj, const int idj, const int tj,
 const float  rc,
 const float  rc2i,
 const float  rclj2,
 ForceDev f
){
  // LIPS constant
  const float alpha2 = 0.19578f*0.19578f;
  const float au[9] = {
       0.0125143224110408f,
      -0.603493863454666f,
      11.7355819865242f,
     -96.296895305654f,
     216.649868508398f,
    -197.409191110696f,
      59.9544311773618f,
      13.9564907382725f,
      -8.66620089071555f
  };
  const float af[9] = {
     2.f*au[0],
     4.f*au[1],
     6.f*au[2],
     8.f*au[3],
    10.f*au[4],
    12.f*au[5],
    14.f*au[6],
    16.f*au[7],
    18.f*au[8]
  };
  const float rc3i = rc2i*rc2i*rc;
  const float bound_ele_pot = rc*rc2i*1.296557;

  // LJ constant
  const float ce12 = 4.f*EPSILON_OXY*powf(SIGMA_OXY,12);
  const float ce06 = 4.f*EPSILON_OXY*powf(SIGMA_OXY, 6);
  const float cf12 = 12.f*ce12;
  const float cf06 =  6.f*ce06;

  const float dx = rmi.x - rmj.x;
  const float dy = rmi.y - rmj.y;
  const float dz = rmi.z - rmj.z;
  const float r2 = ((dx*dx) + dy*dy) + dz*dz;

  // remove molecule itself
  if((idi/3) == (idj/3) && r2 < rclj2) return f;

  const float r2c2 = r2 * rc2i;
  // remove outside cutoff radii
  if(r2c2 > 1.f) return f;
  const float rinv  = rsqrtf(r2);
  const float coef  = r2c2 - alpha2;
  const float coef2 = coef*coef;
  const float utmp = r2c2*(au[0] + r2c2*
			  (au[1] + r2c2*
			  (au[2] + r2c2*
			  (au[3] + r2c2*
			  (au[4] + r2c2*
			  (au[5] + r2c2*
			  (au[6] + r2c2*
			  (au[7] + r2c2*
			  (au[8])))))))));
  const float ftmp =      (af[0] + r2c2*
		          (af[1] + r2c2*
		          (af[2] + r2c2*
		          (af[3] + r2c2*
		          (af[4] + r2c2*
		          (af[5] + r2c2*
		          (af[6] + r2c2*
		          (af[7] + r2c2*
		          (af[8])))))))));
  float fcl = rmj.w * (rinv*rinv*rinv + 0.5f*rc3i*coef2*(6.f*utmp + coef*ftmp));
  float ucl = rmj.w * (rinv - 0.5f*rc2i*rc * coef * coef2 * utmp - bound_ele_pot);
  f.fcl.x += rmi.w*fcl*dx;
  f.fcl.y += rmi.w*fcl*dy;
  f.fcl.z += rmi.w*fcl*dz;
  f.fcl.w += rmi.w*ucl;
  // lj interaction
  if(r2 > rclj2 || ti*tj == 0) return f;
  const float r2i  = rinv*rinv;
  const float r6i  = r2i*r2i*r2i;
  const float r12i = r6i*r6i;
  const float flj  = (cf12 * r12i - cf06 * r6i)*r2i;
  const float ulj  =  ce12 * r12i - ce06 * r6i;
  f.flj.x += flj*dx;
  f.flj.y += flj*dy;
  f.flj.z += flj*dz;
  f.flj.w += ulj;
  return f;
}

inline __device__
float4 warp_reduce(const float4 val){
  float4 ret = val;
#if (__CUDA_ARCH__ >= 300)
  for(int mask=16;mask>0;mask/=2){
    ret.x += __shfl_xor(ret.x,mask);
    ret.y += __shfl_xor(ret.y,mask);
    ret.z += __shfl_xor(ret.z,mask);
    ret.w += __shfl_xor(ret.w,mask);
  }
#else
  __shared__ float4 sh[32];
  sh[threadIdx.x] = val;
  for(int mask=16;mask>0;mask/=2){
    if(threadIdx.x < mask){
      sh[threadIdx.x].x += sh[threadIdx.x+mask].x;
      sh[threadIdx.x].y += sh[threadIdx.x+mask].y;
      sh[threadIdx.x].z += sh[threadIdx.x+mask].z;
      sh[threadIdx.x].w += sh[threadIdx.x+mask].w;
    }
  }
  ret = sh[threadIdx.x];
#endif
  return ret;
}

#if 1
__device__ ForceDev ForceKernel_1walk(
		EpjDev       *jpsh,
		const float4  rmi,
		const int     idi,
		const int     iw,
		const int    *ij_disp,
		const EpjDev *epj,
		ForceDev      f,
		const float   rc,
		const float   rc2i,
		const float   rclj2)

{
    const int tid = threadIdx.x;
    const int j_head = ij_disp[iw  ];
    const int j_tail = ij_disp[iw+1];
    const int ti = (idi%3==0) ? 1:0;
    for(int j=j_head; j<j_tail; j+=WARP_SIZE){
      // __syncthreads();
      jpsh[tid] = ((EpjDev *)(epj + j)) [tid];
      // __syncthreads();
      if(j_tail-j < WARP_SIZE){
	for(int jj=0; jj<j_tail-j; jj++){
	  const float4 rmj = jpsh[jj].rm;
	  const int idj = jpsh[jj].id;
	  const int tj  = (idj%3==0) ? 1:0;
	  f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
	}
      }else{
#pragma unroll
	for(int jj=0; jj<WARP_SIZE; jj++){
	  const float4 rmj = jpsh[jj].rm;
	  const int idj = jpsh[jj].id;
	  const int tj  = (idj%3==0) ? 1:0;
	  f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
	}
      }
    }
    return f;
}

__device__ ForceDev ForceKernel_2walk(
		EpjDev        jpsh[2][WARP_SIZE],
		const float4  rmi,
		const int     idi,
		const int     iw,
		const int     iwalk0,
		const int     iwalk1,
		const int    *ij_disp,
		const EpjDev *epj,
		ForceDev      f,
		const float   rc,
		const float   rc2i,
		const float   rclj2)
{
  const int jbeg0 = ij_disp[iwalk0];
  const int jbeg1 = ij_disp[iwalk1];
  const int jend0 = ij_disp[iwalk0 + 1];
  const int jend1 = ij_disp[iwalk1 + 1];
  const int nj0   = jend0 - jbeg0;
  const int nj1   = jend1 - jbeg1;

  const int nj_longer  = nj0 > nj1 ? nj0 : nj1;
  const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
  const int walk_longer= nj0 > nj1 ? 0 : 1;
  const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;

  const int mywalk = iw==iwalk0 ? 0 : 1;

  const int tid = threadIdx.x;
  const int ti = (idi%3==0) ? 1:0;
  for(int j=0; j<nj_shorter; j+=WARP_SIZE){
    jpsh[0][tid] = ((EpjDev *)(epj + jbeg0 + j)) [tid];
    jpsh[1][tid] = ((EpjDev *)(epj + jbeg1 + j)) [tid];
    if(nj_shorter-j < WARP_SIZE){
      for(int jj=0; jj<nj_shorter-j; jj++){
	const float4 rmj = jpsh[mywalk][jj].rm;
	const int    idj = jpsh[mywalk][jj].id;
	const int    tj  = (idj%3==0)?1:0;
	f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
      }
    }else{
#pragma unroll
      for(int jj=0; jj<WARP_SIZE; jj++){
	const float4 rmj = jpsh[mywalk][jj].rm;
	const int    idj = jpsh[mywalk][jj].id;
	const int    tj  = (idj%3==0)?1:0;
	f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
      }
    }
  }
  for(int j=nj_shorter; j<nj_longer; j+=WARP_SIZE){
    jpsh[0][tid] = ((EpjDev *)(epj + jbeg_longer +  j)) [tid];
    int jrem = nj_longer - j;
    if(jrem < WARP_SIZE){
      for(int jj=0; jj<jrem; jj++){
	if(mywalk == walk_longer){
	  const float4 rmj = jpsh[0][jj].rm;
	  const int    idj = jpsh[0][jj].id;
	  const int    tj  = (idj%3==0)?1:0;
	  f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
	}
      }
    }else{
#pragma unroll
      for(int jj=0; jj<WARP_SIZE; jj++){
	if(mywalk == walk_longer){
	  const float4 rmj = jpsh[0][jj].rm;
	  const int    idj = jpsh[0][jj].id;
	  const int    tj  = (idj%3==0)?1:0;
	  f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
	}
      }
    }
  }
  return f;
}

__device__ ForceDev ForceKernel_multiwalk(const float4  rmi,
					  const int     idi,
					  const int     iw,
					  const int    *ij_disp,
					  const EpjDev *epj,
					  ForceDev      f,
					  const float   rc,
					  const float   rc2i,
					  const float   rclj2)
{
  const int j_head = ij_disp[iw  ];
  const int j_tail = ij_disp[iw+1];

  const int ti = idi%3==0 ? 1:0;
#if 1
  for(int j=j_head; j<j_tail; j++){
    const float4 rmj = epj[j].rm;
    const int    idj = epj[j].id;
    const int    tj  = (idj%3==0)?1:0;
    f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
  }
#else
  int njmin = j_tail - j_head;
  njmin = min(njmin, __shfl_xor(njmin, 1));
  njmin = min(njmin, __shfl_xor(njmin, 2));
  njmin = min(njmin, __shfl_xor(njmin, 4));
  njmin = min(njmin, __shfl_xor(njmin, 8));
  njmin = min(njmin, __shfl_xor(njmin, 16));

  njmin &= 3;;
  for(int j=0; j<njmin; j+=4){
#pragma unroll 4
    for(int jj=0; jj<4; jj++){
      const float4 rmj = epj[j_head + j + jj].rm;
      const int    idj = epj[j_head + j + jj].id;
      const int    tj  = (idj%3==0)?1:0;
      f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
    }
  }
  for(int j=j_head+njmin; j<j_tail; j++){
    const float4 rmj = epj[j].rm;
    const int    idj = epj[j].id;
    const int    tj  = (idj%3==0)?1:0;
    f = accumulate_lips_sw(rmi,idi,ti,rmj,idj,tj,rc,rc2i,rclj2,f);
  }
#endif
  return f;
}

__global__ void kernel_lips_sw(const int      *j_disp,
			       const EpiDev   *epi,
			       const EpjDev   *epj,
			       ForceDev *force,
			       const float    rccl,
			       const float    rclj)
{
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const float rc    = rccl;
    const float rci   = 1.f / rc;
    const float rc2i  = rci*rci;
    const float rclj2 = rclj*rclj;

    const float4 rmi = epi[tid].rm;
    const int    iw  = epi[tid].w;
    const int   idi  = epi[tid].id;

    ForceDev f;
    f.init();

    int t_head = blockDim.x * blockIdx.x;
    int t_tail = t_head + WARP_SIZE - 1;
    int nwalk_in_block = 1 + (epi[t_tail].w - epi[t_head].w);

    __shared__ EpjDev jpsh[2][WARP_SIZE];

    if(1 == nwalk_in_block){
      f = ForceKernel_1walk(jpsh[0], rmi, idi, iw, j_disp, epj, f, rc,rc2i, rclj2);
#if 1
    } else if(2 == nwalk_in_block){
      int iwalk0 = epi[t_head].w;
      int iwalk1 = epi[t_tail].w;
      f = ForceKernel_2walk(jpsh, rmi, idi, iw, iwalk0, iwalk1, j_disp, epj, f, rc,rc2i, rclj2);
#endif
    } else{
      f = ForceKernel_multiwalk(rmi, idi, iw, j_disp, epj, f, rc, rc2i, rclj2);
    }
    force[tid] = f;
}

#else

__global__
void kernel_lips_sw
(const int      *j_disp,
 const EpiDev   *epi,
 const EpjDev   *epj,
       ForceDev *force,
 const float    rccl,
 const float    rclj)
{
  const int tid = threadIdx.x;
  const int i   = blockIdx.x;

  const float rc    = rccl;
  const float rci   = 1.f / rc;
  const float rc2i  = rci*rci;
  const float rclj2 = rclj*rclj;

  const float4 rmi = epi[i].rm;
  const int    iw  = epi[i].w;
  const int    idi = epi[i].id;
  const int    ti  = (idi%3==0) ? 1:0;
  //ForceDev f = force[i];
  ForceDev f;
  f.init();
  for(int j=j_disp[iw]+tid;j<j_disp[iw+1];j+=WARP_SIZE){
    const float4 rmj = epj[j].rm;
    const int idj = epj[j].id;
    const int tj = (idj%3==0)?1:0;
    f = accumulate_lips_sw
      (rmi,idi,ti,
       rmj,idj,tj,
       rc,rc2i,rclj2,f);
  }
  f.fcl = warp_reduce(f.fcl);
  f.flj = warp_reduce(f.flj);
  if(tid==0) force[i] = f;
}
#endif

static cudaPointer<EpiDev>   dev_epi;
static cudaPointer<EpjDev>   dev_epj;
static cudaPointer<ForceDev> dev_force;
static cudaPointer<int>      ij_disp;
static bool init_call = true;

PS::S32 DispatchKernel
(const PS::S32          tag,
 const PS::S32          n_walk,
 const EP              *epi[],
 const PS::S32          n_epi[],
 const EP              *epj[],
 const PS::S32          n_epj[]
 ){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
      int ndevice = 1;
      cudaGetDeviceCount(&ndevice);
      cudaSetDevice(PS::Comm::getRank()%ndevice);
      dev_epi  .allocate(NI_LIMIT);
      dev_epj  .allocate(NJ_LIMIT);
      dev_force.allocate(NI_LIMIT);
      ij_disp  .allocate(N_WALK_LIMIT+2);
      init_call = false;
    }
    ij_disp[0] = 0;
    for(int k=0; k<n_walk; k++){
      ij_disp[k+1] = ij_disp[k] + n_epj[k];
    }
    ij_disp[n_walk+1] = ij_disp[n_walk];
    assert(ij_disp[n_walk] < NJ_LIMIT);
    ij_disp.htod(n_walk + 2);
    int ni_tot = 0;
    int nj_tot = 0;
    int count = 0;
    int nj_ave = 0;
    for(int iw=0; iw<n_walk; iw++){
      PS::F64vec gc = 0.0;
      for(int i=0; i<n_epi[iw]; i++){
	dev_epi[ni_tot].rm.x = epi[iw][i].pos.x;
	dev_epi[ni_tot].rm.y = epi[iw][i].pos.y;
	dev_epi[ni_tot].rm.z = epi[iw][i].pos.z;
	dev_epi[ni_tot].rm.w = epi[iw][i].charge;
	dev_epi[ni_tot].w = iw;
	dev_epi[ni_tot].id= epi[iw][i].id;
	ni_tot++;
      }
      for(int j=0; j<n_epj[iw]; j++){
	dev_epj[nj_tot].rm.x  = epj[iw][j].pos.x;
	dev_epj[nj_tot].rm.y  = epj[iw][j].pos.y;
	dev_epj[nj_tot].rm.z  = epj[iw][j].pos.z;
	dev_epj[nj_tot].rm.w  = epj[iw][j].charge;
	dev_epj[nj_tot].id = epj[iw][j].id;
	nj_tot++;
      }
      //nj_ave += n_epj[iw];
      //count++;
    }
    //printf("nj_ave = %lf\n",(double)nj_ave/count);
    assert(ni_tot < NI_LIMIT);
    assert(nj_tot < NJ_LIMIT);
    int ni_tot_reg = ni_tot;
    if(ni_tot_reg % WARP_SIZE != 0){
      ni_tot_reg /= WARP_SIZE;
      ni_tot_reg++;
      ni_tot_reg *= WARP_SIZE;
    }
    assert(ni_tot_reg < NI_LIMIT);
    for(int i=ni_tot; i<ni_tot_reg; i++){
      dev_epi[i].w = n_walk;
    }
    dev_epi.htod(ni_tot_reg);
    dev_epj.htod(nj_tot);
    //int nblocks  = ni_tot;
    int nblocks  = ni_tot_reg / WARP_SIZE ;
    int nthreads = WARP_SIZE;
    const float rccl = 28.f;
    const float rclj = 4.f * SIGMA_OXY;
    kernel_lips_sw <<<nblocks, nthreads>>>
      (ij_disp, dev_epi,dev_epj, dev_force,rccl,rclj);
    //cudaDeviceSynchronize();
    return 0;
}

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       Force *force[])
{
    int ni_tot = 0;
    for(int k=0; k<n_walk; k++){
      ni_tot += ni[k];
    }
    dev_force.dtoh(ni_tot);

    for(int iw=0; iw<n_walk; iw++){
      for(int i=0; i<ni[iw]; i++){
	force[iw][i].acc = 0.0;
	force[iw][i].pot = 0.0;
      }
    }
    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
      for(int i=0; i<ni[iw]; i++){
	force[iw][i].acc.x += dev_force[n_cnt].fcl.x;
	force[iw][i].acc.y += dev_force[n_cnt].fcl.y;
	force[iw][i].acc.z += dev_force[n_cnt].fcl.z;
	force[iw][i].pot   += 0.5*dev_force[n_cnt].fcl.w;

	force[iw][i].acc.x += dev_force[n_cnt].flj.x;
	force[iw][i].acc.y += dev_force[n_cnt].flj.y;
	force[iw][i].acc.z += dev_force[n_cnt].flj.z;
	force[iw][i].pot   += 0.5*dev_force[n_cnt].flj.w;
	n_cnt++;
      }
    }
    return 0;
}
