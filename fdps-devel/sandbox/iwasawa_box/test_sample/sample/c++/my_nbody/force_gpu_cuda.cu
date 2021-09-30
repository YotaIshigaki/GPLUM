//#include "class.hpp"
//#include "force.hpp"
#include<particle_simulator.hpp>
#include "cuda_pointer.h"
#include "force_gpu_cuda.hpp"

extern "C"{
    double MPI_Wtime(void);
}

//extern double DT_TREE_HOST;

struct CudaTimer{
    cudaEvent_t beg_;
    cudaEvent_t end_;
    float time_;
    CudaTimer() : time_(0.0){}
    void start(){
#ifdef USE_CUDA_TIMER
        cudaEventCreate(&beg_);
        cudaEventCreate(&end_);
        cudaEventRecord(beg_, 0);
#endif
    }
    void stop(){
#ifdef USE_CUDA_TIMER
        cudaEventRecord(end_, 0);
        cudaEventSynchronize(end_);
        cudaEventElapsedTime(&time_, beg_, end_);
#endif
    }
    float getTime(){
#ifdef USE_CUDA_TIMER
        stop();
#endif
        return time_*1e-3;
    }
};

enum{
        //N_STREAM = 40,
        N_STREAM = 8,
        //N_STREAM = 1,
	N_THREAD_GPU = 32,
        //N_WALK_LIMIT = 200,
	N_WALK_LIMIT = 1000,
        //N_WALK_LIMIT = 2000,
        NI_LIMIT_PER_WALK = 17000,
        NJ_LIMIT_PER_WALK = 200000,
        N_LOOP_LIMIT      = 100,
        N_INTERACTION     = 100000000,
	NI_LIMIT            = N_WALK_LIMIT*NI_LIMIT_PER_WALK/2,
        NJ_LIMIT            = 160000000,
        N_ID_TOTAL = 1000000000,
        N_CUM_WALK_LIMIT = 1000000,
};

#ifdef OPTIMIZED_VERSION
struct EpiGPU{
    float3 pos;
    int    id_walk;
    float3 vel;
};

struct EpjGPU{
    float4 posm;
};

struct ForceGPU{
    float4 accp;
    float3 pos;
    float3 vel;
};
#else
struct EpiGPU{
    float3 pos;
    int    id_walk;
};
struct EpjGPU{
    float4 posm;
};
struct ForceGPU{
    float4 accp;
};
#endif

inline __device__ float4 dev_gravity(float  eps2,
                                     float3 ipos,
                                     float4 jposm,
                                     float4 accp){
    float dx = jposm.x - ipos.x;
    float dy = jposm.y - ipos.y;
    float dz = jposm.z - ipos.z;
    
    float r2   = eps2 + dx*dx + dy*dy + dz*dz;
    float rinv = rsqrtf(r2);
    float pij  = jposm.w * rinv;
    float mri3 = rinv*rinv * pij;

    accp.x += mri3 * dx;
    accp.y += mri3 * dy;
    accp.z += mri3 * dz;
    accp.w -= pij;
    
    return accp;
}


////////////////////////
// NO INDEX VERSTION //

#if 0
__global__ void ForceKernel(
		const int2   * ij_disp,
		const EpiGPU * epi,
		const EpjGPU * epj, 
		ForceGPU     * force,
		const float    eps2)
{
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const float3 ipos = epi[tid].pos;
    const int j_head = ij_disp[epi[tid].id_walk  ].y;
    const int j_tail = ij_disp[epi[tid].id_walk+1].y;

    float4 accp = make_float4(0.f, 0.f, 0.f, 0.f);
    for(int j=j_head; j<j_tail; j++){
        float4 jposm = epj[j].posm;
        accp = dev_gravity(eps2, ipos, jposm, accp);
    }
    force[tid].accp = accp;
}
#else
__device__ float4 ForceKernel_1walk(
		float4       *jpsh,
		const float3  ipos,
		const int     id_walk,
		const int2   *ij_disp,
		const EpjGPU *epj,
		float4        accp,
		const float   eps2)
{
    const int tid = threadIdx.x;
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;

	for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
		// __syncthreads();
		jpsh[tid] = ((float4 *)(epj + j)) [tid];
		// __syncthreads();

		if(j_tail-j < N_THREAD_GPU){
			for(int jj=0; jj<j_tail-j; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[jj], accp);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[jj], accp);
			}
		}
	}
	
	return accp;
}

__device__ float4 ForceKernel_2walk(
		float4        jpsh[2][N_THREAD_GPU],
		const float3  ipos,
		const int     id_walk,
		const int     iwalk0,
		const int     iwalk1,
		const int2   *ij_disp,
		const EpjGPU *epj, 
		float4        accp,
		const float   eps2)
{
	const int jbeg0 = ij_disp[iwalk0].y;
	const int jbeg1 = ij_disp[iwalk1].y;
	const int jend0 = ij_disp[iwalk0 + 1].y;
	const int jend1 = ij_disp[iwalk1 + 1].y;
	const int nj0   = jend0 - jbeg0;
	const int nj1   = jend1 - jbeg1;

	const int nj_longer  = nj0 > nj1 ? nj0 : nj1;
	const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
	const int walk_longer= nj0 > nj1 ? 0 : 1;
	const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;

	const int mywalk = id_walk==iwalk0 ? 0 : 1;

        const int tid = threadIdx.x;
	for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
		jpsh[0][tid] = ((float4 *)(epj + jbeg0 + j)) [tid];
		jpsh[1][tid] = ((float4 *)(epj + jbeg1 + j)) [tid];
		if(nj_shorter-j < N_THREAD_GPU){
			for(int jj=0; jj<nj_shorter-j; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[mywalk][jj], accp);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[mywalk][jj], accp);
			}
		}
	}
	for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
		jpsh[0][tid] = ((float4 *)(epj + jbeg_longer +  j)) [tid];
		int jrem = nj_longer - j;
		if(jrem < N_THREAD_GPU){
			for(int jj=0; jj<jrem; jj++){
				if(mywalk == walk_longer)
				accp = dev_gravity(eps2, ipos, jpsh[0][jj], accp);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				if(mywalk == walk_longer)
				accp = dev_gravity(eps2, ipos, jpsh[0][jj], accp);
			}
		}
	}

	return accp;
}

__device__ float4 ForceKernel_multiwalk(
		const float3  ipos,
		const int     id_walk,
		const int2   *ij_disp,
		const EpjGPU *epj, 
		float4        accp,
		const float   eps2){
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;
    for(int j=j_head; j<j_tail; j++){
        float4 jposm = epj[j].posm;
        accp = dev_gravity(eps2, ipos, jposm, accp);
    }
    return accp;
}

__global__ void ForceKernel(
		const int2   * ij_disp,
		const EpiGPU * epi,
		const EpjGPU * epj, 
		ForceGPU     * force,
		const float    eps2){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    float3 ipos    = epi[tid].pos;
    int    id_walk = epi[tid].id_walk;
    float4 accp    = make_float4(0.f, 0.f, 0.f, 0.f);

    int t_head = blockDim.x * blockIdx.x;
    int t_tail = t_head + N_THREAD_GPU - 1;
    int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);

    __shared__ float4 jpsh[2][N_THREAD_GPU];

    if(1 == nwalk_in_block){
        accp = ForceKernel_1walk(jpsh[0], ipos, id_walk, ij_disp, epj, accp, eps2);
    }
    else if(2 == nwalk_in_block){
        // accp = ForceKernel_multiwalk(ipos, id_walk, ij_disp, epj, accp, eps2);
        int iwalk0 = epi[t_head].id_walk;
        int iwalk1 = epi[t_tail].id_walk;
        accp = ForceKernel_2walk(jpsh, ipos, id_walk, iwalk0, iwalk1, ij_disp, epj, accp, eps2);
    }
    else{
        accp = ForceKernel_multiwalk(ipos, id_walk, ij_disp, epj, accp, eps2);
    }
    force[tid].accp = accp;
}
#endif

static cudaPointer<EpiGPU>   dev_epi;
static cudaPointer<EpjGPU>   dev_epj;
static cudaPointer<ForceGPU> dev_force;
static cudaPointer<int2>     ij_disp;
static cudaPointer<int>      dev_id_epj;
static bool init_call = true;

PS::S32 DispatchKernelWithSP(
                             const PS::S32          tag,
                             const PS::S32          n_walk,
                             const EPI          *epi[],
                             const PS::S32          n_epi[],
                             const EPJ          *epj[],
                             const PS::S32          n_epj[],
                             const PS::SPJMonopole *spj[],
                             const PS::S32          n_spj[]){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
		dev_epi.allocate(NI_LIMIT);
		dev_epj.allocate(NJ_LIMIT);
		dev_force.allocate(NI_LIMIT);
		ij_disp.allocate(N_WALK_LIMIT+2);
		init_call = false;
    }
    const float eps2 = FPGrav::eps * FPGrav::eps;
    ij_disp[0].x = 0;
    ij_disp[0].y = 0;
    for(int k=0; k<n_walk; k++){
        ij_disp[k+1].x = ij_disp[k].x + n_epi[k];
        ij_disp[k+1].y = ij_disp[k].y + (n_epj[k] + n_spj[k]);
    }
    ij_disp[n_walk+1] = ij_disp[n_walk];
    assert(ij_disp[n_walk].x < NI_LIMIT);
    assert(ij_disp[n_walk].y < NJ_LIMIT);
    ij_disp.htod(n_walk + 2);
    const int ni_tot = ij_disp[n_walk].x;
    const int nj_tot = ij_disp[n_walk].y;
    int ni_tot_reg = ni_tot;
    ni_tot_reg = ((ni_tot_reg-1)/N_THREAD_GPU + 1)*N_THREAD_GPU;

    double wtime0 = MPI_Wtime();
#pragma omp parallel for
    for(int iw=0; iw<n_walk; iw++){
        int i_dst = ij_disp[iw].x;
        const int n_epi_tmp = n_epi[iw];
        for(int i=0; i<n_epi_tmp; i++, i_dst++){
            const float x = epi[iw][i].pos.x;
            const float y = epi[iw][i].pos.y;
            const float z = epi[iw][i].pos.z;
            dev_epi[i_dst].pos.x = x;
            dev_epi[i_dst].pos.y = y;
            dev_epi[i_dst].pos.z = z;
            dev_epi[i_dst].id_walk = iw;
        }
    }
    WTIME_COPY_IP += MPI_Wtime() - wtime0;
    wtime0 = MPI_Wtime();
#pragma omp parallel for
    for(int iw=0; iw<n_walk; iw++){
        int j_dst = ij_disp[iw].y;
        const int n_epj_tmp = n_epj[iw];
        for(int j=0; j<n_epj_tmp; j++, j_dst++){
            const float x = epj[iw][j].pos.x;
            const float y = epj[iw][j].pos.y;
            const float z = epj[iw][j].pos.z;
            const float m = epj[iw][j].mass;
            dev_epj[j_dst].posm.x  = x;
            dev_epj[j_dst].posm.y  = y;
            dev_epj[j_dst].posm.z  = z;
            dev_epj[j_dst].posm.w  = m;
        }
        const int n_spj_tmp = n_spj[iw];
        for(int j=0; j<n_spj_tmp; j++, j_dst++){
            const float x = spj[iw][j].pos.x;
            const float y = spj[iw][j].pos.y;
            const float z = spj[iw][j].pos.z;
            const float m = spj[iw][j].getCharge();
            dev_epj[j_dst].posm.x  = x;
            dev_epj[j_dst].posm.y  = y;
            dev_epj[j_dst].posm.z  = z;
            dev_epj[j_dst].posm.w  = m;
        }
    }
    WTIME_COPY_JP += MPI_Wtime() - wtime0;
#pragma omp parallel for
    for(int i=ni_tot; i<ni_tot_reg; i++){
        dev_epi[i].id_walk = n_walk;
    }
    //WTIME_COPY_FOR_H2D += MPI_Wtime() - wtime0;
    CudaTimer timer_h2d_ip;
    timer_h2d_ip.start();
    dev_epi.htod(ni_tot_reg);
    WTIME_H2D_IP += timer_h2d_ip.getTime();
    CudaTimer timer_h2d_list;
    timer_h2d_list.start();
    dev_epj.htod(nj_tot);
    WTIME_H2D_LIST += timer_h2d_list.getTime();
    int nblocks  = ni_tot_reg / N_THREAD_GPU;
    int nthreads = N_THREAD_GPU;
    CudaTimer timer_kernel;
    timer_kernel.start();
    ForceKernel <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_force, eps2);
    WTIME_KERNEL += timer_kernel.getTime();
    return 0;
}


////////////////////
// INDEX VERSTION //
#if 1
__device__ float4 ForceKernelIndex_1walk(float4       *jpsh,
                                         const float3  ipos,
                                         const int     id_walk,
                                         const int2   *ij_disp,
                                         const EpjGPU *epj,
                                         const int    *id_epj,
                                         float4        accp,
                                         const float   eps2){
    const int tid = threadIdx.x;
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;
    for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
        int adr  = id_epj[j+tid];
        if(j+tid >= j_tail){
            adr = id_epj[j_tail-1];
        }
        jpsh[tid] = *((float4*)(epj+adr));
        if(j_tail-j < N_THREAD_GPU){
            for(int jj=0; jj<j_tail-j; jj++){
                accp = dev_gravity(eps2, ipos, jpsh[jj], accp);
            }
        }
        else{
#pragma unroll
            for(int jj=0; jj<N_THREAD_GPU; jj++){
                accp = dev_gravity(eps2, ipos, jpsh[jj], accp);
            }
        }
    }
    return accp;
}

__device__ float4 ForceKernelIndex_2walk(float4        jpsh[2][N_THREAD_GPU],
                                         const float3  ipos,
                                         const int     id_walk,
                                         const int     iwalk0,
                                         const int     iwalk1,
                                         const int2   *ij_disp,
                                         const EpjGPU *epj, 
                                         const int    *id_epj,
                                         float4        accp,
                                         const float   eps2){
    const int jbeg0 = ij_disp[iwalk0].y;
    const int jbeg1 = ij_disp[iwalk1].y;
    const int jend0 = ij_disp[iwalk0 + 1].y;
    const int jend1 = ij_disp[iwalk1 + 1].y;
    const int nj0   = jend0 - jbeg0;
    const int nj1   = jend1 - jbeg1;

    const int nj_longer  = nj0 > nj1 ? nj0 : nj1;
    const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
    const int walk_longer= nj0 > nj1 ? 0 : 1;
    const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;

    const int mywalk = id_walk==iwalk0 ? 0 : 1;

    const int tid = threadIdx.x;
    for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
        const int j_tmp = (j + tid < nj_shorter) ? j+tid : nj_shorter-1;
        const int adr0 = id_epj[jbeg0 + j_tmp];
        const int adr1 = id_epj[jbeg1 + j_tmp];
        jpsh[0][tid] = *((float4 *)(epj + adr0));
        jpsh[1][tid] = *((float4 *)(epj + adr1));
        if(nj_shorter-j < N_THREAD_GPU){
            for(int jj=0; jj<nj_shorter-j; jj++){
                accp = dev_gravity(eps2, ipos, jpsh[mywalk][jj], accp);
            }
        }
        else{
#pragma unroll
            for(int jj=0; jj<N_THREAD_GPU; jj++){
                accp = dev_gravity(eps2, ipos, jpsh[mywalk][jj], accp);
            }
        }
    }
    for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
        const int j_tmp = (j + tid < nj_longer) ? j+tid : nj_longer-1;
        const int adr = id_epj[jbeg_longer + j_tmp];
        jpsh[0][tid] = *((float4 *)(epj + adr));
        int jrem = nj_longer - j;
        if(jrem < N_THREAD_GPU){
            for(int jj=0; jj<jrem; jj++){
                if(mywalk == walk_longer)
                    accp = dev_gravity(eps2, ipos, jpsh[0][jj], accp);
            }
        }
        else{
#pragma unroll
            for(int jj=0; jj<N_THREAD_GPU; jj++){
                if(mywalk == walk_longer)
                    accp = dev_gravity(eps2, ipos, jpsh[0][jj], accp);
            }
        }
    }
    return accp;
}

__device__ float4 ForceKernelIndex_multiwalk(const float3  ipos,
                                             const int     id_walk,
                                             const int2   *ij_disp,
                                             const EpjGPU *epj, 
                                             const int    *id_epj,
                                             float4        accp,
                                             const float   eps2){
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;
    for(int j=j_head; j<j_tail; j++){
        const int adr = id_epj[j];
        float4 jposm = epj[adr].posm;
        accp = dev_gravity(eps2, ipos, jposm, accp);
    }
    return accp;
}

__device__ float3 Kick(const float3 vel,
                       const float4 accp,
                       const double dt){
    float3 vel_new;
    vel_new.x = vel.x + accp.x*dt;
    vel_new.y = vel.y + accp.y*dt;
    vel_new.z = vel.z + accp.z*dt;
    return vel_new;
}

__device__ float3 Drift(const float3 pos,
                        const float3 vel,
                        const double dt){
    float3 pos_new;
    pos_new.x = pos.x + vel.x*dt;
    pos_new.y = pos.y + vel.y*dt;
    pos_new.z = pos.z + vel.z*dt;
    return pos_new;
}


// optimized version
__global__ void ForceKernelIndex(const int2   * ij_disp,
                                 const EpiGPU * epi,
                                 const EpjGPU * epj, 
                                 const int * id_epj,
                                 ForceGPU     * force,
                                 const float    eps2){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    float3 ipos    = epi[tid].pos;
    int    id_walk = epi[tid].id_walk;
    float4 accp    = make_float4(0.f, 0.f, 0.f, 0.f);

    int t_head = blockDim.x * blockIdx.x;
    int t_tail = t_head + N_THREAD_GPU - 1;
    int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);

    __shared__ float4 jpsh[2][N_THREAD_GPU];
    if(1 == nwalk_in_block){
        accp = ForceKernelIndex_1walk(jpsh[0], ipos, id_walk, ij_disp, epj, id_epj, accp, eps2);
    } 
    else if(2 == nwalk_in_block){
        int iwalk0 = epi[t_head].id_walk;
        int iwalk1 = epi[t_tail].id_walk;
        accp = ForceKernelIndex_2walk(jpsh, ipos, id_walk, iwalk0, iwalk1, ij_disp, epj, id_epj, 
                                      accp, eps2);
    } 
    else{
        accp = ForceKernelIndex_multiwalk(ipos, id_walk, ij_disp, epj, id_epj, accp, eps2);
    }
    force[tid].accp = accp;
}
#else
// non-optimized version
__global__ void ForceKernelIndex(const int2   * ij_disp,
                                 const EpiGPU * epi,
                                 const EpjGPU * epj,
                                 const int * id_epj, 
                                 ForceGPU     * force,
                                 const float    eps2){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const float3 ipos = epi[tid].pos;
    const int j_head = ij_disp[epi[tid].id_walk  ].y;
    const int j_tail = ij_disp[epi[tid].id_walk+1].y;
    float4 accp = make_float4(0.f, 0.f, 0.f, 0.f);
    for(int j=j_head; j<j_tail; j++){
        int adr      = id_epj[j];
        float4 jposm = epj[adr].posm;
        accp = dev_gravity(eps2, ipos, jposm, accp);
    }
    force[tid].accp = accp;
}
#endif


//////////////////////////////
// INDEX: YES, STREAM: YES  //
static cudaPointer<EpiGPU>   dev_epi_2[N_STREAM];
static cudaPointer<ForceGPU> dev_force_2[N_STREAM];
static cudaPointer<int2>     ij_disp_2[N_STREAM];
static cudaPointer<int>      dev_id_epj_2[N_STREAM];
static cudaStream_t          stream[N_STREAM];
static int n_walk_ar[N_STREAM];
static int n_disp_walk_ar[N_STREAM+1];
//static int ni_tot_reg[N_STREAM];
static int id_epj_ar[N_STREAM+1];
static int id_spj_ar[N_STREAM+1];

//// UTIL ////// 
template<class Tptcl>
void H2dJptcl(const Tptcl pj, 
              const int n_tot,
              cudaPointer<EpjGPU> dev_pj,
              const int n_offset = 0){
    const int n_ave = n_tot/N_STREAM;
    for(int id_stream=0; id_stream<N_STREAM; id_stream++){
        const int n_in_stream = n_ave + ((n_tot%N_STREAM > id_stream) ? 1 : 0);
        const int id_head = n_ave*id_stream + std::min(id_stream, n_tot%N_STREAM);
        const int id_end  = id_head + n_in_stream;
        double wtime1 = MPI_Wtime();
#pragma omp parallel for
        for(int i=id_head; i<id_end; i++){
            const int i_dst = i + n_offset;
            const float x = pj[i].pos.x;
            const float y = pj[i].pos.y;
            const float z = pj[i].pos.z;
            const float m = pj[i].mass;
            dev_pj[i_dst].posm.x  = x;
            dev_pj[i_dst].posm.y  = y;
            dev_pj[i_dst].posm.z  = z;
            dev_pj[i_dst].posm.w  = m;
        }
        //WTIME_COPY_JP += MPI_Wtime() - wtime1;
        CudaTimer timer_h2d_jp;
        timer_h2d_jp.start();
        dev_pj.htod(id_head+n_offset, n_in_stream, stream[id_stream]);
        WTIME_H2D_ALL_PTCL += timer_h2d_jp.getTime();
    }
}


PS::S32 DispatchKernelIndexStream(const PS::S32 tag,
                                  const PS::S32 n_walk,
                                  const EPI ** epi,
                                  const PS::S32 *  n_epi,
                                  const PS::S32 ** id_epj,
                                  const PS::S32 *  n_epj,
                                  const PS::S32 ** id_spj,
                                  const PS::S32 *  n_spj,
                                  const EPJ  * epj,
                                  const PS::S32 n_epj_tot,
                                  const PS::SPJMonopole * spj,
                                  const PS::S32 n_spj_tot,
                                  const bool send_flag){
    //cudaThreadSynchronize();
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
        dev_epj.allocate(NJ_LIMIT);
        dev_epi    .allocate(NI_LIMIT);
        dev_force  .allocate(NI_LIMIT);
        ij_disp    .allocate(N_WALK_LIMIT+2);
        dev_id_epj .allocate(NJ_LIMIT);
        for(int i=0; i<N_STREAM; i++){
            dev_epi_2[i]    .allocate(NI_LIMIT);
            dev_force_2[i]  .allocate(NI_LIMIT);
            ij_disp_2[i]    .allocate(N_WALK_LIMIT+2);
            dev_id_epj_2[i] .allocate(NJ_LIMIT/N_STREAM);
            cudaStreamCreate(&stream[i]);
        }
        init_call = false;
    }
    if(send_flag==true){
        double wtime0 = MPI_Wtime();
        assert(NJ_LIMIT > n_epj_tot+n_spj_tot);
        H2dJptcl(epj, n_epj_tot, dev_epj);
        H2dJptcl(spj, n_spj_tot, dev_epj, n_epj_tot);
        return 0;
    }
    const int n_walk_ave = n_walk/N_STREAM;
    for(int id_stream=0; id_stream<N_STREAM; id_stream++){
        const int n_walk_in_stream = n_walk_ave + ((id_stream < n_walk%N_STREAM) ? 1 : 0);
        const int id_walk_head = n_walk_ave*id_stream + std::min(id_stream, n_walk%N_STREAM);
        const int id_walk_end  = id_walk_head + n_walk_in_stream;
        ij_disp_2[id_stream][0].x = ij_disp_2[id_stream][0].y = 0;
        for(int iw=0; iw<n_walk_in_stream; iw++){
            const int iw_src = iw+id_walk_head;
            ij_disp_2[id_stream][iw+1].x = ij_disp_2[id_stream][iw].x + n_epi[iw_src];
            ij_disp_2[id_stream][iw+1].y = ij_disp_2[id_stream][iw].y + (n_epj[iw_src]+n_spj[iw_src]);
        }
        ij_disp_2[id_stream][n_walk_in_stream+1] = ij_disp_2[id_stream][n_walk_in_stream];
        assert(ij_disp_2[id_stream][n_walk_in_stream].x < NI_LIMIT);
        assert(ij_disp_2[id_stream][n_walk_in_stream].y < NJ_LIMIT);
        ij_disp_2[id_stream].htod(n_walk_in_stream+2, stream[id_stream]);
        
        int ni_tot_reg = ij_disp_2[id_stream][n_walk_in_stream].x;
        ni_tot_reg = ((ni_tot_reg-1)/N_THREAD_GPU + 1)*N_THREAD_GPU;
        double wtime0 = MPI_Wtime();
#pragma omp parallel for
        for(int iw=0; iw<n_walk_in_stream; iw++){
            const int iw_src = id_walk_head + iw;
            const int n_epi_tmp = n_epi[iw_src];
            int i_dst = ij_disp_2[id_stream][iw].x;
            //double wtime1 = MPI_Wtime();
            for(int i=0; i<n_epi_tmp; i++, i_dst++){
                const float x = epi[iw_src][i].pos.x;
                const float y = epi[iw_src][i].pos.y;
                const float z = epi[iw_src][i].pos.z;
                dev_epi_2[id_stream][i_dst].pos.x   = x;
                dev_epi_2[id_stream][i_dst].pos.y   = y;
                dev_epi_2[id_stream][i_dst].pos.z   = z;
                dev_epi_2[id_stream][i_dst].id_walk = iw;
            }
            //double wtime2 = MPI_Wtime();
            //WTIME_COPY_IP += wtime2 - wtime1;
            int j_dst = ij_disp_2[id_stream][iw].y;
            const int n_epj_tmp = n_epj[iw_src];
            for(int j=0; j<n_epj_tmp; j++, j_dst++){
                dev_id_epj_2[id_stream][j_dst] = id_epj[iw_src][j];
            }
            const int n_spj_tmp = n_spj[iw_src];
            for(int j=0; j<n_spj_tmp; j++, j_dst++){
                dev_id_epj_2[id_stream][j_dst] = id_spj[iw_src][j]+n_epj_tot;
            }
            //WTIME_COPY_ID += MPI_Wtime() - wtime2;
        }
        WTIME_COPY_IP_ID += MPI_Wtime() - wtime0;

        const int ni_tot = ij_disp_2[id_stream][n_walk_in_stream].x;
        const int nj_tot = ij_disp_2[id_stream][n_walk_in_stream].y;
        for(int i=ni_tot; i<ni_tot_reg; i++){
            dev_epi_2[id_stream][i].id_walk = n_walk_in_stream;
        }
        //WTIME_COPY_FOR_H2D += MPI_Wtime() - wtime0;

        dev_epi_2[id_stream].htod(ni_tot_reg, stream[id_stream]);
        dev_id_epj_2[id_stream].htod(nj_tot, stream[id_stream]);

        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        const float eps2 = FPGrav::eps * FPGrav::eps;
        CudaTimer timer_kernel;
        timer_kernel.start();
        ForceKernelIndex <<<nblocks, nthreads, 0, stream[id_stream]>>> (ij_disp_2[id_stream], dev_epi_2[id_stream], dev_epj, dev_id_epj_2[id_stream], dev_force_2[id_stream], eps2);
        WTIME_KERNEL += timer_kernel.getTime();
    }
    return 0;
}

PS::S32 RetrieveKernelStream(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const PS::S32 ni[],
                             Force    *force[]){
    //cudaThreadSynchronize();
    static int n_offset[N_WALK_LIMIT+1];
    const int n_walk_ave = n_walk/N_STREAM;
    //for(int i=0; i<N_STREAM+1; i++){ 
    //    n_disp_walk_ar[i] = n_walk_ave*i + std::min(i, n_walk%N_STREAM);
    //}
    for(int id_stream=0; id_stream<N_STREAM; id_stream++){
        //const int id_walk_head = n_disp_walk_ar[id_stream];
        //const int id_walk_end  = n_disp_walk_ar[id_stream+1];
        //const int n_walk_in_stream = n_walk_ar[id_stream];
        const int n_walk_in_stream = n_walk_ave + ((id_stream < n_walk%N_STREAM) ? 1 : 0);
        const int id_walk_head = n_walk_ave*id_stream + std::min(id_stream, n_walk%N_STREAM);
        const int id_walk_end  = id_walk_head + n_walk_in_stream;
        const int ni_tot = ij_disp_2[id_stream][n_walk_in_stream].x;
        dev_force_2[id_stream].dtoh(ni_tot, stream[id_stream]);
        cudaStreamSynchronize(stream[id_stream]);
        n_offset[0] = 0;
        for(int i=0; i<n_walk_in_stream; i++){
            const int iw = id_walk_head + i;
            n_offset[i+1] = n_offset[i] + ni[iw];
        }
        double wtime0 = MPI_Wtime();
#pragma omp parallel for
        for(int iw=id_walk_head; iw<id_walk_end; iw++){
            const int ni_tmp = ni[iw];
            const int n_offset_tmp = n_offset[iw-id_walk_head];
            for(int i=0; i<ni_tmp; i++){
                const int i_src = i + n_offset_tmp;
#ifdef OPTIMIZED_VERSION
                const double px = dev_force_2[id_stream][i_src].pos.x;
                const double py = dev_force_2[id_stream][i_src].pos.y;
                const double pz = dev_force_2[id_stream][i_src].pos.z;
                const double vx = dev_force_2[id_stream][i_src].vel.x;
                const double vy = dev_force_2[id_stream][i_src].vel.y;
                const double vz = dev_force_2[id_stream][i_src].vel.z;
                force[iw][i].pos.x = px;
                force[iw][i].pos.y = py;
                force[iw][i].pos.z = pz;
                force[iw][i].vel.x = vx;
                force[iw][i].vel.y = vy;
                force[iw][i].vel.z = vz;
#else
                const double ax = dev_force_2[id_stream][i_src].accp.x;
                const double ay = dev_force_2[id_stream][i_src].accp.y;
                const double az = dev_force_2[id_stream][i_src].accp.z;
                const double p  = dev_force_2[id_stream][i_src].accp.w;
                force[iw][i].acc.x = ax;
                force[iw][i].acc.y = ay;
                force[iw][i].acc.z = az;
                force[iw][i].pot   = p;
#endif
            }
        }
        WTIME_COPY_FORCE += MPI_Wtime() - wtime0;
    }
    return 0;
}


/////////////////////////////
// INDEX: YES, STREAM: NO  //
PS::S32 DispatchKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const EPI ** epi,
                            const PS::S32 *  n_epi,
                            const PS::S32 ** id_epj,
                            const PS::S32 *  n_epj,
                            const PS::S32 ** id_spj,
                            const PS::S32 *  n_spj,
                            const EPJ  * epj,
                            const PS::S32 n_epj_tot,
                            const PS::SPJMonopole * spj,
                            const PS::S32 n_spj_tot,
                            const bool send_flag){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
        dev_epi    .allocate(NI_LIMIT);
        dev_force  .allocate(NI_LIMIT);
        ij_disp    .allocate(N_WALK_LIMIT+2);
        dev_epj    .allocate(NJ_LIMIT);
        dev_id_epj .allocate(NJ_LIMIT);
        init_call = false;
    }
    if(send_flag==true){
        if(NJ_LIMIT < n_epj_tot+n_spj_tot){
            std::cerr<<"n_epj_tot+n_spj_tot= "<<n_epj_tot+n_spj_tot
                     <<" NJ_LIMIT= "<<NJ_LIMIT
                     <<std::endl;
        }
        assert(NJ_LIMIT > n_epj_tot+n_spj_tot);
        double wtime0 = MPI_Wtime();
#pragma omp parallel for
        for(int i=0; i<n_epj_tot; i++){
            const float x = epj[i].pos.x;
            const float y = epj[i].pos.y;
            const float z = epj[i].pos.z;
            const float m = epj[i].mass;
            dev_epj[i].posm.x  = x;
            dev_epj[i].posm.y  = y;
            dev_epj[i].posm.z  = z;
            dev_epj[i].posm.w  = m;
        }
#pragma omp parallel for
        for(int i=0; i<n_spj_tot; i++){
            const int i_dst = i + n_epj_tot;
            const float x = spj[i].pos.x;
            const float y = spj[i].pos.y;
            const float z = spj[i].pos.z;
            const float m = spj[i].mass;
            dev_epj[i_dst].posm.x  = x;
            dev_epj[i_dst].posm.y  = y;
            dev_epj[i_dst].posm.z  = z;
            dev_epj[i_dst].posm.w  = m;
        }
        WTIME_COPY_JP += MPI_Wtime() - wtime0;

        //CudaTimer timer_send_jp;
        //timer_send_jp.start();
        //dev_epj.htod(n_epj_tot+n_spj_tot);
        //WTIME_SEND_JP += timer_send_jp.getTime();

        CudaTimer timer_h2d_all_ptcl;
        timer_h2d_all_ptcl.start();
        dev_epj.htod(n_epj_tot+n_spj_tot);
        WTIME_H2D_ALL_PTCL += timer_h2d_all_ptcl.getTime();
        return 0;
    }

    ij_disp[0].x = 0;
    ij_disp[0].y = 0;
    for(int k=0; k<n_walk; k++){
        ij_disp[k+1].x = ij_disp[k].x + n_epi[k];
        ij_disp[k+1].y = ij_disp[k].y + (n_epj[k] + n_spj[k]);
    }
    ij_disp[n_walk+1] = ij_disp[n_walk];
    assert(ij_disp[n_walk].x < NI_LIMIT);
    assert(ij_disp[n_walk].y < NJ_LIMIT);
    const int ni_tot = ij_disp[n_walk].x;
    const int nj_tot = ij_disp[n_walk].y;

    ij_disp.htod(n_walk + 2);
    int ni_tot_reg = ij_disp[n_walk].x;
    if(ni_tot_reg % N_THREAD_GPU){
        ni_tot_reg /= N_THREAD_GPU;
        ni_tot_reg++;
        ni_tot_reg *= N_THREAD_GPU;
    }

    double wtime0 = MPI_Wtime();
#pragma omp parallel for
    for(int iw=0; iw<n_walk; iw++){
        int i_dst = ij_disp[iw].x;
        const int n_epi_tmp = n_epi[iw];
        for(int i=0; i<n_epi_tmp; i++, i_dst++){
            const float x = epi[iw][i].pos.x;
            const float y = epi[iw][i].pos.y;
            const float z = epi[iw][i].pos.z;
            dev_epi[i_dst].pos.x = x;
            dev_epi[i_dst].pos.y = y;
            dev_epi[i_dst].pos.z = z;
            dev_epi[i_dst].id_walk = iw;
        }
    }
    double wtime1 = MPI_Wtime();
    WTIME_COPY_IP += wtime1 - wtime0;
#pragma omp parallel for
    for(int iw=0; iw<n_walk; iw++){
        int j_dst = ij_disp[iw].y;
        const int n_epj_tmp = n_epj[iw];
        for(int j=0; j<n_epj_tmp; j++, j_dst++){
            const int id = id_epj[iw][j];
            dev_id_epj[j_dst] = id;
        }
        const int n_spj_tmp = n_spj[iw];
        for(int j=0; j<n_spj_tmp; j++, j_dst++){
            const int id = id_spj[iw][j]+n_epj_tot;
            dev_id_epj[j_dst] = id;
        }
    }
    WTIME_COPY_ID += MPI_Wtime() - wtime1;

#pragma omp parallel for
    for(int i=ni_tot; i<ni_tot_reg; i++){
        dev_epi[i].id_walk = n_walk;
    }

    //WTIME_COPY_FOR_H2D += MPI_Wtime() - wtime0;

    CudaTimer timer_h2d_ip;
    timer_h2d_ip.start();
    dev_epi.htod(ni_tot_reg);
    WTIME_H2D_IP += timer_h2d_ip.getTime();

    CudaTimer timer_h2d_list;
    timer_h2d_list.start();
    dev_id_epj.htod(nj_tot);
    WTIME_H2D_LIST += timer_h2d_list.getTime();

    int nblocks  = ni_tot_reg / N_THREAD_GPU;
    int nthreads = N_THREAD_GPU;
    const float eps2 = FPGrav::eps * FPGrav::eps;

    CudaTimer timer_kernel;
    timer_kernel.start();
    ForceKernelIndex <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_id_epj, dev_force, eps2);
    WTIME_KERNEL += timer_kernel.getTime();

    return 0;
}


PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       Force    *force[]){
    const int ni_tot = ij_disp[n_walk].x;
    //CudaTimer timer_d2h_force;
    //timer_d2h_force.start();
    //dev_force.dtoh(ni_tot);
    //WTIME_RECV_FORCE += timer_d2h_force.getTime();

    CudaTimer timer_d2h_force;
    timer_d2h_force.start();
    dev_force.dtoh(ni_tot);
    WTIME_D2H_FORCE += timer_d2h_force.getTime();

    double wtime0 = MPI_Wtime();
#pragma omp parallel for
    for(int iw=0; iw<n_walk; iw++){
        const int i_head = ij_disp[iw].x;
        const int ni_tmp = ni[iw];
        for(int i=0; i<ni_tmp; i++){
            const int i_src = i+i_head;
            const double ax = dev_force[i_src].accp.x;
            const double ay = dev_force[i_src].accp.y;
            const double az = dev_force[i_src].accp.z;
            const double p = dev_force[i_src].accp.w;
            force[iw][i].acc.x = ax;
            force[iw][i].acc.y = ay;
            force[iw][i].acc.z = az;
            force[iw][i].pot   = p;
        }
    }
    WTIME_COPY_FORCE += MPI_Wtime() - wtime0;
    return 0;
}

/*
PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       ForceGrav    *force[]){
    const int ni_tot = ij_disp[n_walk].x;
    CudaTimer timer_d2h;
    timer_d2h.start();
    dev_force.dtoh(ni_tot);
    WTIME_D2H += timer_d2h.getTime();
    double wtime0 = MPI_Wtime();
#pragma omp parallel for
    for(int iw=0; iw<n_walk; iw++){
        const int i_head = ij_disp[iw].x;
        const int ni_tmp = ni[iw];
        for(int i=0; i<ni_tmp; i++){
            const int i_src = i+i_head;
            const double ax = dev_force[i_src].accp.x;
            const double ay = dev_force[i_src].accp.y;
            const double az = dev_force[i_src].accp.z;
            const double p = dev_force[i_src].accp.w;
            force[iw][i].acc.x = ax;
            force[iw][i].acc.y = ay;
            force[iw][i].acc.z = az;
            force[iw][i].pot   = p;
        }
    }
    WTIME_COPY_FOR_D2H += MPI_Wtime() - wtime0;
    return 0;
}
*/



static cudaPointer<EpjGPU>   dev_epj_2[N_STREAM];

//////////////////////////////
// INDEX: NO, STREAM: YES  //
PS::S32 DispatchKernelStream(const PS::S32          tag,
                             const PS::S32          n_walk,
                             const EPI          *epi[],
                             const PS::S32          n_epi[],
                             const EPJ          *epj[],
                             const PS::S32          n_epj[],
                             const PS::SPJMonopole *spj[],
                             const PS::S32          n_spj[]){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
        dev_epj.allocate(NJ_LIMIT);
        dev_epi    .allocate(NI_LIMIT);
        dev_force  .allocate(NI_LIMIT);
        ij_disp    .allocate(N_WALK_LIMIT+2);
        dev_id_epj .allocate(NJ_LIMIT);
        for(int i=0; i<N_STREAM; i++){
            dev_epi_2[i]    .allocate(NI_LIMIT);
            dev_epj_2[i]    .allocate(NI_LIMIT);
            dev_force_2[i]  .allocate(NI_LIMIT);
            ij_disp_2[i]    .allocate(N_WALK_LIMIT+2);
            dev_id_epj_2[i] .allocate(NJ_LIMIT);
            cudaStreamCreate(&stream[i]);
        }
            init_call = false;
    }
    //cudaThreadSynchronize();
    const int n_walk_ave = n_walk/N_STREAM;
    for(int i=0; i<N_STREAM;   i++){ n_walk_ar[i]  = n_walk_ave   + ((i < n_walk%N_STREAM) ? 1 : 0);}
    for(int i=0; i<N_STREAM+1; i++){ n_disp_walk_ar[i] = n_walk_ave*i + std::min(i, n_walk%N_STREAM);}

    for(int id_stream=0; id_stream<N_STREAM; id_stream++){
        const int n_walk_in_stream   = n_walk_ar[id_stream];
        const int id_walk_head = n_disp_walk_ar[id_stream];
        const int id_walk_end  = n_disp_walk_ar[id_stream+1];
        ij_disp_2[id_stream][0].x = ij_disp_2[id_stream][0].y = 0;
        for(int iw=0; iw<n_walk_in_stream; iw++){
            const int iw_src = iw+id_walk_head;
            ij_disp_2[id_stream][iw+1].x = ij_disp_2[id_stream][iw].x + n_epi[iw_src];
            ij_disp_2[id_stream][iw+1].y = ij_disp_2[id_stream][iw].y + (n_epj[iw_src]+n_spj[iw_src]);
        }
        ij_disp_2[id_stream][n_walk_in_stream+1] = ij_disp_2[id_stream][n_walk_in_stream];
        assert(ij_disp_2[id_stream][n_walk_in_stream].x < NI_LIMIT);
        assert(ij_disp_2[id_stream][n_walk_in_stream].y < NJ_LIMIT);
        ij_disp_2[id_stream].htod(n_walk_in_stream+2, stream[id_stream]);
        
        int ni_tot_reg = ij_disp_2[id_stream][n_walk_in_stream].x;
        if(ni_tot_reg % N_THREAD_GPU){
            ni_tot_reg /= N_THREAD_GPU;
            ni_tot_reg++;
            ni_tot_reg *= N_THREAD_GPU;
        }
        double wtime0 = MPI_Wtime();
#pragma omp parallel for
        for(int iw=id_walk_head; iw<id_walk_end; iw++){
            const int iw_tmp = iw - id_walk_head;
            const int n_epi_tmp = n_epi[iw];
            int i_dst = ij_disp_2[id_stream][iw_tmp].x;
            //double wtime1 = MPI_Wtime();
            for(int i=0; i<n_epi_tmp; i++, i_dst++){
                const float x = epi[iw][i].pos.x;
                const float y = epi[iw][i].pos.y;
                const float z = epi[iw][i].pos.z;
                dev_epi_2[id_stream][i_dst].pos.x   = x;
                dev_epi_2[id_stream][i_dst].pos.y   = y;
                dev_epi_2[id_stream][i_dst].pos.z   = z;
                dev_epi_2[id_stream][i_dst].id_walk = iw_tmp;
            }
            //double wtime2 = MPI_Wtime();
            //WTIME_COPY_IP += wtime2 - wtime1;
            int j_dst = ij_disp_2[id_stream][iw_tmp].y;
            const int n_epj_tmp = n_epj[iw];
            for(int j=0; j<n_epj_tmp; j++, j_dst++){
                const float m = epj[iw][j].mass;
                const float x = epj[iw][j].pos.x;
                const float y = epj[iw][j].pos.y;
                const float z = epj[iw][j].pos.z;
                dev_epj_2[id_stream][j_dst].posm.x = x;
                dev_epj_2[id_stream][j_dst].posm.y = y;
                dev_epj_2[id_stream][j_dst].posm.z = z;
                dev_epj_2[id_stream][j_dst].posm.w = m;
            }
            const int n_spj_tmp = n_spj[iw];
            for(int j=0; j<n_spj_tmp; j++, j_dst++){
                const float m = spj[iw][j].mass;
                const float x = spj[iw][j].pos.x;
                const float y = spj[iw][j].pos.y;
                const float z = spj[iw][j].pos.z;
                dev_epj_2[id_stream][j_dst].posm.x = x;
                dev_epj_2[id_stream][j_dst].posm.y = y;
                dev_epj_2[id_stream][j_dst].posm.z = z;
                dev_epj_2[id_stream][j_dst].posm.w = m;
            }
            //WTIME_COPY_ID += MPI_Wtime() - wtime2;
        }
        WTIME_COPY_IP_JP += MPI_Wtime() - wtime0;
        const int ni_tot = ij_disp_2[id_stream][n_walk_in_stream].x;
        const int nj_tot = ij_disp_2[id_stream][n_walk_in_stream].y;
        for(int i=ni_tot; i<ni_tot_reg; i++){
            dev_epi_2[id_stream][i].id_walk = n_walk_in_stream;
        }
        //WTIME_COPY_FOR_H2D += MPI_Wtime() - wtime0;

        dev_epi_2[id_stream].htod(ni_tot_reg, stream[id_stream]);
        dev_epj_2[id_stream].htod(nj_tot, stream[id_stream]);
        //ni_tot_reg[id_stream] = ni_tot_reg_tmp;

        //int nblocks  = ni_tot_reg[id_stream] / N_THREAD_GPU;
        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        const float eps2 = FPGrav::eps * FPGrav::eps;
        ForceKernel <<<nblocks, nthreads, 0, stream[id_stream]>>> (ij_disp_2[id_stream], dev_epi_2[id_stream], dev_epj_2[id_stream], dev_force_2[id_stream], eps2);
        //CudaTimer timer_kernel;
        //timer_kernel.start();
        //ForceKernelIndex <<<nblocks, nthreads, 0, stream[id_stream]>>> (ij_disp_2[id_stream], dev_epi_2[id_stream], dev_epj, dev_id_epj_2[id_stream], dev_force_2[id_stream], eps2, DT_TREE_HOST);
        //WTIME_KERNEL += timer_kernel.getTime();
    }
    return 0;
}


////////////////////////////////////////
// SEND INDEX ONLY AT THE FIRST STEP  //
extern int CONSTRUCTION_STEP;
static int * id_epj_h[N_STREAM];
static int * id_epj_d[N_STREAM];
static int ID_EPJ_GLB_OFFSET[N_STREAM];

PS::S32 DispatchKernelIndexStream2(const PS::S32 tag,
                                   const PS::S32 n_walk,
                                   const EPI ** epi,
                                   const PS::S32 *  n_epi,
                                   const PS::S32 ** id_epj,
                                   const PS::S32 *  n_epj,
                                   const PS::S32 ** id_spj,
                                   const PS::S32 *  n_spj,
                                   const EPJ  * epj,
                                   const PS::S32 n_epj_tot,
                                   const PS::SPJMonopole * spj,
                                   const PS::S32 n_spj_tot,
                                   const bool send_flag){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
        dev_epj.allocate(NJ_LIMIT);
        //cudaMallocHost((void**)&id_epj_h,  NJ_LIMIT * sizeof(int));
        //cudaMalloc((void**)&id_epj_d,  N_ID_TOTAL * sizeof(int));
        for(int i=0; i<N_STREAM; i++){
            ij_disp_2[i].allocate(N_WALK_LIMIT+2);
            dev_epi_2[i].allocate(NI_LIMIT);
            dev_force_2[i].allocate(NI_LIMIT);
            cudaMallocHost((void**)&id_epj_h[i],  NJ_LIMIT/N_STREAM * sizeof(int));
            cudaMalloc((void**)&id_epj_d[i],  N_ID_TOTAL/N_STREAM * sizeof(int));
            cudaStreamCreate(&stream[i]);
        }
        init_call = false;
    }
    if(send_flag==true){
        double wtime0 = MPI_Wtime();
        assert(NJ_LIMIT > n_epj_tot+n_spj_tot);
        H2dJptcl(epj, n_epj_tot, dev_epj);
        H2dJptcl(spj, n_spj_tot, dev_epj, n_epj_tot);
        for(int i=0; i<N_STREAM; i++){
            ID_EPJ_GLB_OFFSET[i] = 0;
        }
        return 0;
    }
    const int n_walk_ave = n_walk/N_STREAM;
    for(int id_stream=0; id_stream<N_STREAM; id_stream++){
        //std::cerr<<"id_stream= "<<id_stream<<std::endl;
        const int n_walk_in_stream = n_walk_ave + ((id_stream < n_walk%N_STREAM) ? 1 : 0);
        const int id_walk_head = n_walk_ave*id_stream + std::min(id_stream, n_walk%N_STREAM);
        ij_disp_2[id_stream][0].x = 0;
        ij_disp_2[id_stream][0].y = ID_EPJ_GLB_OFFSET[id_stream];
        for(int iw=0; iw<n_walk_in_stream; iw++){
            const int iw_src = iw + id_walk_head;
            ij_disp_2[id_stream][iw+1].x = ij_disp_2[id_stream][iw].x +  n_epi[iw_src];
            ij_disp_2[id_stream][iw+1].y = ij_disp_2[id_stream][iw].y + (n_epj[iw_src]+n_spj[iw_src]);
        }
        ij_disp_2[id_stream][n_walk_in_stream+1] = ij_disp_2[id_stream][n_walk_in_stream];
        ij_disp_2[id_stream].htod(n_walk_in_stream+2, stream[id_stream]);
        int ni_tot_reg = ij_disp_2[id_stream][n_walk_in_stream].x;
        ni_tot_reg = ((ni_tot_reg-1)/N_THREAD_GPU + 1)*N_THREAD_GPU;
        //std::cerr<<"ni_tot_reg= "<<ni_tot_reg<<std::endl;
        double wtime0 = MPI_Wtime();
#pragma omp parallel for
        for(int iw=0; iw<n_walk_in_stream; iw++){
            const int iw_src = id_walk_head + iw;
            const int n_epi_tmp = n_epi[iw_src];
            int i_dst = ij_disp_2[id_stream][iw].x;
            //std::cerr<<"iw_src= "<<iw_src
            //         <<" n_epi_tmp= "<<n_epi_tmp
            //         <<" i_dst= "<<i_dst
            //         <<std::endl;
            for(int ip=0; ip<n_epi_tmp; ip++, i_dst++){
                const float x = epi[iw_src][ip].pos.x;
                const float y = epi[iw_src][ip].pos.y;
                const float z = epi[iw_src][ip].pos.z;
                dev_epi_2[id_stream][i_dst].pos.x   = x;
                dev_epi_2[id_stream][i_dst].pos.y   = y;
                dev_epi_2[id_stream][i_dst].pos.z   = z;
                dev_epi_2[id_stream][i_dst].id_walk = iw;
            }
            if(CONSTRUCTION_STEP==1){
                //double wtime2 = MPI_Wtime();
                //WTIME_COPY_IP += wtime2 - wtime1;
                int j_dst = ij_disp_2[id_stream][iw].y - ID_EPJ_GLB_OFFSET[id_stream];
                //std::cerr<<"j_dst= "<<j_dst<<std::endl;
                const int n_epj_tmp = n_epj[iw_src];
                //std::cerr<<"n_epj_tmp= "<<n_epj_tmp<<std::endl;
                for(int j=0; j<n_epj_tmp; j++, j_dst++){
                    //dev_id_epj_2[id_stream][j_dst] = id_epj[iw_src][j];
                    id_epj_h[id_stream][j_dst] = id_epj[iw_src][j];
                }
                const int n_spj_tmp = n_spj[iw_src];
                //std::cerr<<"n_spj_tmp= "<<n_spj_tmp<<std::endl;
                for(int j=0; j<n_spj_tmp; j++, j_dst++){
                    //dev_id_epj_2[id_stream][j_dst] = id_spj[iw_src][j]+n_epj_tot;
                    id_epj_h[id_stream][j_dst] = id_spj[iw_src][j]+n_epj_tot;
                }
                //WTIME_COPY_ID += MPI_Wtime() - wtime2;
            }
        }
        WTIME_COPY_IP_ID += MPI_Wtime() - wtime0;

        const int ni_tot = ij_disp_2[id_stream][n_walk_in_stream].x;
        const int nj_tot = ij_disp_2[id_stream][n_walk_in_stream].y - ID_EPJ_GLB_OFFSET[id_stream];
        for(int i=ni_tot; i<ni_tot_reg; i++){
            dev_epi_2[id_stream][i].id_walk = n_walk_in_stream;
        }
        //std::cerr<<"ID_EPJ_GLB_OFFSET[id_stream]= "<<ID_EPJ_GLB_OFFSET[id_stream]
        //         <<" ni_tot= "<<ni_tot
        //         <<" nj_tot= "<<nj_tot
        //         <<std::endl;
        //WTIME_COPY_FOR_H2D += MPI_Wtime() - wtime0;
        dev_epi_2[id_stream].htod(ni_tot_reg, stream[id_stream]);
        if(CONSTRUCTION_STEP == 1){
            CUDA_SAFE_CALL(cudaMemcpyAsync(id_epj_d[id_stream]+ID_EPJ_GLB_OFFSET[id_stream], id_epj_h[id_stream], nj_tot*sizeof(int), cudaMemcpyHostToDevice, stream[id_stream]));
        }
        ID_EPJ_GLB_OFFSET[id_stream] += nj_tot;
        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        const float eps2 = FPGrav::eps * FPGrav::eps;
        //std::cerr<<"nblocks= "<<nblocks<<" nthreads= "<<nthreads<<std::endl;
        CudaTimer timer_kernel;
        timer_kernel.start();
        ForceKernelIndex <<<nblocks, nthreads, 0, stream[id_stream]>>> (ij_disp_2[id_stream], dev_epi_2[id_stream], dev_epj, id_epj_d[id_stream], dev_force_2[id_stream], eps2);
        WTIME_KERNEL += timer_kernel.getTime();
    }
    return 0;
}


