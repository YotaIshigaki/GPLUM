//#include "class.hpp"
//#include "force.hpp"
#include<particle_simulator.hpp>
#include "cuda_pointer.h"
#include "force_gpu_cuda.hpp"

enum{
        N_STREAM = 4,
	N_THREAD_GPU = 32,
	N_WALK_LIMIT = 1000,
	NI_LIMIT            = N_WALK_LIMIT*1000,
	NI_LIMIT_PER_STREAM = NI_LIMIT / N_STREAM,
	NJ_LIMIT            = N_WALK_LIMIT*50000,
	NJ_LIMIT_PER_STREAM = NJ_LIMIT / N_STREAM,
};

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
		const float   eps2)
{
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;

#if 1
    for(int j=j_head; j<j_tail; j++){
		float4 jposm = epj[j].posm;
		accp = dev_gravity(eps2, ipos, jposm, accp);
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
			float4 jposm = epj[j_head + j + jj].posm;
			float4 jposm = jpf[jj];
			accp = dev_gravity(eps2, ipos, jposm, accp);
		}
	}
    for(int j=j_head+njmin; j<j_tail; j++){
		float4 jposm = epj[j].posm;
		accp = dev_gravity(eps2, ipos, jposm, accp);
	}
#endif
	return accp;
}

__global__ void ForceKernel(
		const int2   * ij_disp,
		const EpiGPU * epi,
		const EpjGPU * epj, 
		ForceGPU     * force,
		const float    eps2)
{
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
	} else if(2 == nwalk_in_block){
            // accp = ForceKernel_multiwalk(ipos, id_walk, ij_disp, epj, accp, eps2);
            int iwalk0 = epi[t_head].id_walk;
            int iwalk1 = epi[t_tail].id_walk;
            accp = ForceKernel_2walk(jpsh, ipos, id_walk, iwalk0, iwalk1, ij_disp, epj, accp, eps2);
	} else{
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
                             const FPGrav          *epi[],
                             const PS::S32          n_epi[],
                             const FPGrav          *epj[],
                             const PS::S32          n_epj[],
                             const PS::SPJMonopole *spj[],
                             const PS::S32          n_spj[]){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
		dev_epi  .allocate(NI_LIMIT);
		dev_epj  .allocate(NJ_LIMIT);
		dev_force.allocate(NI_LIMIT);
		ij_disp  .allocate(N_WALK_LIMIT+2);
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

    int ni_tot_reg = ij_disp[n_walk].x;
    if(ni_tot_reg % N_THREAD_GPU){
        ni_tot_reg /= N_THREAD_GPU;
        ni_tot_reg++;
        ni_tot_reg *= N_THREAD_GPU;
    }

    int ni_tot = 0;
    int nj_tot = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<n_epi[iw]; i++){
            dev_epi[ni_tot].pos.x = epi[iw][i].pos.x;
            dev_epi[ni_tot].pos.y = epi[iw][i].pos.y;
            dev_epi[ni_tot].pos.z = epi[iw][i].pos.z;
            dev_epi[ni_tot].id_walk = iw;
            ni_tot++;
        }
        for(int j=0; j<n_epj[iw]; j++){
            dev_epj[nj_tot].posm.x  = epj[iw][j].pos.x;
            dev_epj[nj_tot].posm.y  = epj[iw][j].pos.y;
            dev_epj[nj_tot].posm.z  = epj[iw][j].pos.z;
            dev_epj[nj_tot].posm.w  = epj[iw][j].mass;
            nj_tot++;
        }
        for(int j=0; j<n_spj[iw]; j++){
            dev_epj[nj_tot].posm.x  = spj[iw][j].pos.x;
            dev_epj[nj_tot].posm.y  = spj[iw][j].pos.y;
            dev_epj[nj_tot].posm.z  = spj[iw][j].pos.z;
            dev_epj[nj_tot].posm.w  = spj[iw][j].getCharge();
            nj_tot++;
        }
    }
    for(int i=ni_tot; i<ni_tot_reg; i++){
        dev_epi[i].id_walk = n_walk;
    }

    dev_epi.htod(ni_tot_reg);
    dev_epj.htod(nj_tot);

    int nblocks  = ni_tot_reg / N_THREAD_GPU;
    int nthreads = N_THREAD_GPU;
    ForceKernel <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_force, eps2);

    return 0;
}




//////////////////
// for index

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
#if 1
    for(int j=j_head; j<j_tail; j++){
        const int adr = id_epj[j];
        float4 jposm = epj[adr].posm;
        accp = dev_gravity(eps2, ipos, jposm, accp);
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
            float4 jposm = epj[j_head + j + jj].posm;
            float4 jposm = jpf[jj];
            accp = dev_gravity(eps2, ipos, jposm, accp);
        }
    }
    for(int j=j_head+njmin; j<j_tail; j++){
        float4 jposm = epj[j].posm;
        accp = dev_gravity(eps2, ipos, jposm, accp);
    }
#endif
    return accp;
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
    //float mass_tmp = 0.0;
    for(int j=j_head; j<j_tail; j++){
        int adr      = id_epj[j];
        float4 jposm = epj[adr].posm;
        accp = dev_gravity(eps2, ipos, jposm, accp);
        //mass_tmp += jposm.w;
    }
    //if(mass_tmp != 1.0) printf("mass_tmp=%f \n", mass_tmp);
    //assert(mass_tmp == 1.0);
    force[tid].accp = accp;
}
#endif

#if 1


static cudaPointer<EpiGPU>   dev_epi_2[N_STREAM];
static cudaPointer<ForceGPU> dev_force_2[N_STREAM];
static cudaPointer<int2>     ij_disp_2[N_STREAM];
static cudaPointer<int>      dev_id_epj_2[N_STREAM];
static cudaStream_t          stream[N_STREAM];

    #if 0

PS::S32 DispatchKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const FPGrav ** epi,
                            const PS::S32 *  n_epi,
                            const PS::S32 ** id_epj,
                            const PS::S32 *  n_epj,
                            const PS::S32 ** id_spj,
                            const PS::S32 *  n_spj,
                            const FPGrav  * epj,
                            const PS::S32 n_epj_tot,
                            const PS::SPJMonopole * spj,
                            const PS::S32 n_spj_tot,
                            const bool send_flag){
    assert(n_walk <= N_WALK_LIMIT);

    static int id_iptcl_head[N_WALK_LIMIT+1];
    static int id_jptcl_head[N_WALK_LIMIT+1];

    if(init_call){
        dev_epj.allocate(NJ_LIMIT);
        for(int i=0; i<N_STREAM; i++){
            dev_epi_2[i]    .allocate(NI_LIMIT);
            dev_force_2[i]  .allocate(NI_LIMIT);
            ij_disp_2[i]    .allocate(N_WALK_LIMIT+2);
            dev_id_epj_2[i] .allocate(NJ_LIMIT);
            cudaStreamCreate(&stream[i]);
        }
        init_call = false;
    }
    if(send_flag==true){
        assert(NJ_LIMIT > n_epj_tot+n_spj_tot);
#if 1
        /*
        cudaEvent_t beg_send_epj, end_send_epj;
        float time_send_epj_tmp = 0.0;
        cudaEventCreate(&beg_send_epj);
        cudaEventCreate(&end_send_epj);
        cudaEventRecord(beg_send_epj, 0);
        */
        int id_epj_ar[N_STREAM+1];
        const int n_epj_ave = n_epj_tot/N_STREAM;
        int id_spj_ar[N_STREAM+1];
        const int n_spj_ave = n_spj_tot/N_STREAM;
        for(int i=0; i<N_STREAM+1; i++){ 
            id_epj_ar[i] = n_epj_ave*i + std::min(i, n_epj_tot%N_STREAM);
            id_spj_ar[i] = n_spj_ave*i + std::min(i, n_spj_tot%N_STREAM);
        }
        //#pragma omp parallel for num_threads(N_STREAM)
        for(int id_stream=0; id_stream<N_STREAM; id_stream++){
            int id_epj_head = id_epj_ar[id_stream];
            int id_epj_end  = id_epj_ar[id_stream+1];
            int n_epj_tmp = id_epj_end - id_epj_head;
            for(int i=id_epj_head; i<id_epj_end; i++){
                dev_epj[i].posm.x  = epj[i].pos.x;
                dev_epj[i].posm.y  = epj[i].pos.y;
                dev_epj[i].posm.z  = epj[i].pos.z;
                dev_epj[i].posm.w  = epj[i].mass;
            }
            dev_epj.htod(id_epj_head, n_epj_tmp, stream[id_stream]);
            //cudaStreamSynchronize(stream[id_stream]);
        }
        //#pragma omp parallel for num_threads(N_STREAM)
        for(int id_stream=0; id_stream<N_STREAM; id_stream++){
            int id_spj_head = id_spj_ar[id_stream];
            int id_spj_end  = id_spj_ar[id_stream+1];
            int n_spj_tmp = id_spj_end - id_spj_head;
            for(int i=id_spj_head; i<id_spj_end; i++){
                int i_tmp = i + n_epj_tot;
                dev_epj[i_tmp].posm.x  = spj[i].pos.x;
                dev_epj[i_tmp].posm.y  = spj[i].pos.y;
                dev_epj[i_tmp].posm.z  = spj[i].pos.z;
                dev_epj[i_tmp].posm.w  = spj[i].mass;
            }
            dev_epj.htod(id_spj_head+n_epj_tot, n_spj_tmp, stream[id_stream]);
            //cudaStreamSynchronize(stream[id_stream]);
        }
        /*
        cudaEventRecord(end_send_epj, 0);
        cudaEventSynchronize(end_send_epj);
        cudaEventElapsedTime(&time_send_epj_tmp, beg_send_epj, end_send_epj);
        WTIME_SEND_EPJ += time_send_epj_tmp*1e-3;
        */
#else

        cudaEvent_t beg_send_epj, end_send_epj;
        float time_send_epj_tmp = 0.0;
        cudaEventCreate(&beg_send_epj);
        cudaEventCreate(&end_send_epj);
        cudaEventRecord(beg_send_epj, 0);

        //#pragma omp parallel for
        for(PS::S32 i=0; i<n_epj_tot; i++){
            dev_epj[i].posm.x  = epj[i].pos.x;
            dev_epj[i].posm.y  = epj[i].pos.y;
            dev_epj[i].posm.z  = epj[i].pos.z;
            dev_epj[i].posm.w  = epj[i].mass;
        }

        //#pragma omp parallel for
        for(PS::S32 i=0; i<n_spj_tot; i++){
            const int i_tmp = n_epj_tot + i;
            dev_epj[i_tmp].posm.x  = spj[i].pos.x;
            dev_epj[i_tmp].posm.y  = spj[i].pos.y;
            dev_epj[i_tmp].posm.z  = spj[i].pos.z;
            dev_epj[i_tmp].posm.w  = spj[i].mass;
        }
        dev_epj.htod(n_epj_tot+n_spj_tot);

        cudaEventRecord(end_send_epj, 0);
        cudaEventSynchronize(end_send_epj);
        cudaEventElapsedTime(&time_send_epj_tmp, beg_send_epj, end_send_epj);
        WTIME_SEND_EPJ += time_send_epj_tmp*1e-3;

#endif
        return 0;
    }
    else{
        int n_walk_ar[N_STREAM];
        int id_walk_ar[N_STREAM+1];
        int ni_tot_reg[N_STREAM];
        const int n_walk_ave = n_walk/N_STREAM;
        for(int i=0; i<N_STREAM;   i++){ n_walk_ar[i]  = n_walk_ave   + ((i < n_walk%N_STREAM) ? 1 : 0);}
        for(int i=0; i<N_STREAM+1; i++){ id_walk_ar[i] = n_walk_ave*i + std::min(i, n_walk%N_STREAM);}

        id_iptcl_head[0] = id_jptcl_head[0] = 0;
        for(int i=0; i<n_walk; i++){
            id_iptcl_head[i+1] = id_iptcl_head[i] + n_epi[i];
            id_jptcl_head[i+1] = id_jptcl_head[i] + (n_epj[i]+n_spj[i]);
        }

        /*
        cudaEvent_t beg_h2d, end_h2d;
        float time_h2d_tmp = 0.0;
        cudaEventCreate(&beg_h2d);
        cudaEventCreate(&end_h2d);
        cudaEventRecord(beg_h2d, 0);
        */
        for(int id_stream=0; id_stream<N_STREAM; id_stream++){
            int n_walk_tmp   = n_walk_ar[id_stream];
            int id_walk_head = id_walk_ar[id_stream];
            int id_walk_end  = id_walk_ar[id_stream+1];

            ij_disp_2[id_stream][0].x = ij_disp_2[id_stream][0].y = 0;
            int iw_tmp = 0;
            for(int iw=id_walk_head; iw<id_walk_end; iw++, iw_tmp++){
                ij_disp_2[id_stream][iw_tmp+1].x = ij_disp_2[id_stream][iw_tmp].x + n_epi[iw];
                ij_disp_2[id_stream][iw_tmp+1].y = ij_disp_2[id_stream][iw_tmp].y + (n_epj[iw]+n_spj[iw]);
            }
            ij_disp_2[id_stream][n_walk_tmp+1] = ij_disp_2[id_stream][n_walk_tmp];
            assert(ij_disp_2[id_stream][n_walk_tmp].x < NI_LIMIT);
            assert(ij_disp_2[id_stream][n_walk_tmp].y < NJ_LIMIT);
            ij_disp_2[id_stream].htod(n_walk_tmp + 2, stream[id_stream]);

            int ni_tot_reg_tmp = ij_disp_2[id_stream][n_walk_tmp].x;
            if(ni_tot_reg_tmp % N_THREAD_GPU){
                ni_tot_reg_tmp /= N_THREAD_GPU;
                ni_tot_reg_tmp++;
                ni_tot_reg_tmp *= N_THREAD_GPU;
            }

#if 0
            int ni_tot = 0;
            int nj_tot = 0;
            iw_tmp = 0;
            for(int iw=id_walk_head; iw<id_walk_end; iw++, iw_tmp++){
                for(int i=0; i<n_epi[iw]; i++, ni_tot++){
                    dev_epi_2[id_stream][ni_tot].pos.x = epi[iw][i].pos.x;
                    dev_epi_2[id_stream][ni_tot].pos.y = epi[iw][i].pos.y;
                    dev_epi_2[id_stream][ni_tot].pos.z = epi[iw][i].pos.z;
                    dev_epi_2[id_stream][ni_tot].id_walk = iw_tmp;
                }
                for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                    dev_id_epj_2[id_stream][nj_tot] = id_epj[iw][j];
                }
                for(int j=0; j<n_spj[iw]; j++, nj_tot++){
                    dev_id_epj_2[id_stream][nj_tot] = id_spj[iw][j]+n_epj_tot;
                }
            }
            dev_epi_2[id_stream].htod(ni_tot_reg_tmp, stream[id_stream]);
            dev_id_epj_2[id_stream].htod(nj_tot, stream[id_stream]);

#else
    #if 1
#pragma omp parallel for
            for(int iw=id_walk_head; iw<id_walk_end; iw++){
                const int iw_tmp = iw - id_walk_head;
                int i_tmp = id_iptcl_head[iw] - id_iptcl_head[id_walk_head];
                for(int i=0; i<n_epi[iw]; i++, i_tmp++){
                    dev_epi_2[id_stream][i_tmp].pos.x = epi[iw][i].pos.x;
                    dev_epi_2[id_stream][i_tmp].pos.y = epi[iw][i].pos.y;
                    dev_epi_2[id_stream][i_tmp].pos.z = epi[iw][i].pos.z;
                    dev_epi_2[id_stream][i_tmp].id_walk = iw_tmp;
                }
                int j_tmp = id_jptcl_head[iw] - id_jptcl_head[id_walk_head];
                for(int j=0; j<n_epj[iw]; j++, j_tmp++){
                    dev_id_epj_2[id_stream][j_tmp] = id_epj[iw][j];
                }
                for(int j=0; j<n_spj[iw]; j++, j_tmp++){
                    dev_id_epj_2[id_stream][j_tmp] = id_spj[iw][j]+n_epj_tot;
                }
            }
            int ni_tot = id_iptcl_head[id_walk_end] - id_iptcl_head[id_walk_head];
            int nj_tot = id_jptcl_head[id_walk_end] - id_jptcl_head[id_walk_head];
            dev_epi_2[id_stream].htod(ni_tot_reg_tmp, stream[id_stream]);
            dev_id_epj_2[id_stream].htod(nj_tot, stream[id_stream]);

    #else

            int ni_tot = 0;
            int nj_tot = 0;
            iw_tmp = 0;
            for(int iw=id_walk_head; iw<id_walk_end; iw++, iw_tmp++){
                for(int i=0; i<n_epi[iw]; i++, ni_tot++){
                    dev_epi_2[id_stream][ni_tot].pos.x = epi[iw][i].pos.x;
                    dev_epi_2[id_stream][ni_tot].pos.y = epi[iw][i].pos.y;
                    dev_epi_2[id_stream][ni_tot].pos.z = epi[iw][i].pos.z;
                    dev_epi_2[id_stream][ni_tot].id_walk = iw_tmp;
                }
            }
            dev_epi_2[id_stream].htod(ni_tot_reg_tmp, stream[id_stream]);
            for(int iw=id_walk_head; iw<id_walk_end; iw++, iw_tmp++){
                for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                    dev_id_epj_2[id_stream][nj_tot] = id_epj[iw][j];
                }
                for(int j=0; j<n_spj[iw]; j++, nj_tot++){
                    dev_id_epj_2[id_stream][nj_tot] = id_spj[iw][j]+n_epj_tot;
                }
            }
            dev_id_epj_2[id_stream].htod(nj_tot, stream[id_stream]);
    #endif
#endif

            for(int i=ni_tot; i<ni_tot_reg_tmp; i++){
                dev_epi_2[id_stream][i].id_walk = n_walk_tmp;
            }
            ni_tot_reg[id_stream] = ni_tot_reg_tmp;

            //cudaStreamSynchronize(stream[id_stream]);
        }
            /*
        cudaEventRecord(end_h2d, 0);
        cudaEventSynchronize(end_h2d);
        cudaEventElapsedTime(&time_h2d_tmp, beg_h2d, end_h2d);
        WTIME_H2D += time_h2d_tmp*1e-3;
            */

        /*
        cudaEvent_t beg, end;
        float time_kernel_tmp = 0.0;
        cudaEventCreate(&beg);
        cudaEventCreate(&end);
        cudaEventRecord(beg, 0);
        */
        for(int id_stream=0; id_stream<N_STREAM; id_stream++){
            int nblocks  = ni_tot_reg[id_stream] / N_THREAD_GPU;
            int nthreads = N_THREAD_GPU;
            const float eps2 = FPGrav::eps * FPGrav::eps;
            ForceKernelIndex <<<nblocks, nthreads, 0, stream[id_stream]>>> (ij_disp_2[id_stream], dev_epi_2[id_stream], dev_epj, dev_id_epj_2[id_stream], dev_force_2[id_stream], eps2);
            //cudaStreamSynchronize(stream[id_stream]);
        }
        /*
        cudaEventRecord(end, 0);
        cudaEventSynchronize(end);
        cudaEventElapsedTime( &time_kernel_tmp, beg, end);
        WTIME_KERNEL += time_kernel_tmp*1e-3;
        */
        return 0;
    }
}


    #else
PS::S32 DispatchKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const FPGrav ** epi,
                            const PS::S32 *  n_epi,
                            const PS::S32 ** id_epj,
                            const PS::S32 *  n_epj,
                            const PS::S32 ** id_spj,
                            const PS::S32 *  n_spj,
                            const FPGrav  * epj,
                            const PS::S32 n_epj_tot,
                            const PS::SPJMonopole * spj,
                            const PS::S32 n_spj_tot,
                            const bool send_flag){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
        dev_epj.allocate(NJ_LIMIT);
        for(int i=0; i<N_STREAM; i++){
            dev_epi_2[i]    .allocate(NI_LIMIT_PER_STREAM);
            dev_force_2[i]  .allocate(NI_LIMIT_PER_STREAM);
            ij_disp_2[i]    .allocate(N_WALK_LIMIT+2);
            dev_id_epj_2[i] .allocate(NJ_LIMIT_PER_STREAM);
            cudaStreamCreate(&stream[i]);
        }
        init_call = false;
    }
    if(send_flag==true){
        assert(NJ_LIMIT > n_epj_tot+n_spj_tot);
#if 1
        int id_epj_ar[N_STREAM+1];
        const int n_epj_ave = n_epj_tot/N_STREAM;
        int id_spj_ar[N_STREAM+1];
        const int n_spj_ave = n_spj_tot/N_STREAM;
        for(int i=0; i<N_STREAM+1; i++){ 
            id_epj_ar[i] = n_epj_ave*i + std::min(i, n_epj_tot%N_STREAM);
            id_spj_ar[i] = n_spj_ave*i + std::min(i, n_spj_tot%N_STREAM);
        }
        for(int id_stream=0; id_stream<N_STREAM; id_stream++){
            int id_epj_head = id_epj_ar[id_stream];
            int id_epj_end  = id_epj_ar[id_stream+1];
            int n_epj_tmp = id_epj_end - id_epj_head;
            for(int i=id_epj_head; i<id_epj_end; i++){
                dev_epj[i].posm.x  = epj[i].pos.x;
                dev_epj[i].posm.y  = epj[i].pos.y;
                dev_epj[i].posm.z  = epj[i].pos.z;
                dev_epj[i].posm.w  = epj[i].mass;
            }
            dev_epj.htod(id_epj_head, n_epj_tmp, stream[id_stream]);
        }
        for(int id_stream=0; id_stream<N_STREAM; id_stream++){
            int id_spj_head = id_spj_ar[id_stream];
            int id_spj_end  = id_spj_ar[id_stream+1];
            int n_spj_tmp = id_spj_end - id_spj_head;
            for(int i=id_spj_head; i<id_spj_end; i++){
                int i_tmp = i + n_epj_tot;
                dev_epj[i_tmp].posm.x  = spj[i].pos.x;
                dev_epj[i_tmp].posm.y  = spj[i].pos.y;
                dev_epj[i_tmp].posm.z  = spj[i].pos.z;
                dev_epj[i_tmp].posm.w  = spj[i].mass;
            }
            dev_epj.htod(id_spj_head+n_epj_tot, n_spj_tmp, stream[id_stream]);
        }
#else
        for(PS::S32 i=0; i<n_epj_tot; i++){
            dev_epj[i].posm.x  = epj[i].pos.x;
            dev_epj[i].posm.y  = epj[i].pos.y;
            dev_epj[i].posm.z  = epj[i].pos.z;
            dev_epj[i].posm.w  = epj[i].mass;
        }
        int nj_tmp = n_epj_tot;
        for(PS::S32 i=0; i<n_spj_tot; i++, nj_tmp++){
            dev_epj[nj_tmp].posm.x  = spj[i].pos.x;
            dev_epj[nj_tmp].posm.y  = spj[i].pos.y;
            dev_epj[nj_tmp].posm.z  = spj[i].pos.z;
            dev_epj[nj_tmp].posm.w  = spj[i].mass;
        }
        dev_epj.htod(n_epj_tot+n_spj_tot);
#endif

        return 0;
    }
    else{
        int n_walk_ar[N_STREAM];
        int id_walk_ar[N_STREAM+1];
        const int n_walk_ave = n_walk/N_STREAM;
        for(int i=0; i<N_STREAM;   i++){ n_walk_ar[i]  = n_walk_ave   + ((i < n_walk%N_STREAM) ? 1 : 0);}
        for(int i=0; i<N_STREAM+1; i++){ id_walk_ar[i] = n_walk_ave*i + std::min(i, n_walk%N_STREAM);}
        for(int id_stream=0; id_stream<N_STREAM; id_stream++){
            int n_walk_tmp   = n_walk_ar[id_stream];
            int id_walk_head = id_walk_ar[id_stream];
            int id_walk_end  = id_walk_ar[id_stream+1];

            ij_disp_2[id_stream][0].x = ij_disp_2[id_stream][0].y = 0;
            int iw_tmp = 0;
            for(int iw=id_walk_head; iw<id_walk_end; iw++, iw_tmp++){
                ij_disp_2[id_stream][iw_tmp+1].x = ij_disp_2[id_stream][iw_tmp].x + n_epi[iw];
                ij_disp_2[id_stream][iw_tmp+1].y = ij_disp_2[id_stream][iw_tmp].y + (n_epj[iw]+n_spj[iw]);
            }
            ij_disp_2[id_stream][n_walk_tmp+1] = ij_disp_2[id_stream][n_walk_tmp];
            assert(ij_disp_2[id_stream][n_walk_tmp].x < NI_LIMIT_PER_STREAM);
            assert(ij_disp_2[id_stream][n_walk_tmp].y < NJ_LIMIT_PER_STREAM);

            ij_disp_2[id_stream].htod(n_walk_tmp + 2, stream[id_stream]);
            //cudaDeviceSynchronize();

            int ni_tot_reg = ij_disp_2[id_stream][n_walk_tmp].x;
            if(ni_tot_reg % N_THREAD_GPU){
                ni_tot_reg /= N_THREAD_GPU;
                ni_tot_reg++;
                ni_tot_reg *= N_THREAD_GPU;
            }

            int ni_tot = 0;
            int nj_tot = 0;
            iw_tmp = 0;
            for(int iw=id_walk_head; iw<id_walk_end; iw++, iw_tmp++){
                for(int i=0; i<n_epi[iw]; i++, ni_tot++){
                    dev_epi_2[id_stream][ni_tot].pos.x = epi[iw][i].pos.x;
                    dev_epi_2[id_stream][ni_tot].pos.y = epi[iw][i].pos.y;
                    dev_epi_2[id_stream][ni_tot].pos.z = epi[iw][i].pos.z;
                    dev_epi_2[id_stream][ni_tot].id_walk = iw_tmp;
                }
            }
            for(int i=ni_tot; i<ni_tot_reg; i++){
                dev_epi_2[id_stream][i].id_walk = n_walk_tmp;
            }
            dev_epi_2[id_stream].htod(ni_tot_reg, stream[id_stream]);
            for(int iw=id_walk_head; iw<id_walk_end; iw++, iw_tmp++){
                for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                    dev_id_epj_2[id_stream][nj_tot] = id_epj[iw][j];
                }
                for(int j=0; j<n_spj[iw]; j++, nj_tot++){
                    dev_id_epj_2[id_stream][nj_tot] = id_spj[iw][j]+n_epj_tot;
                }
            }
            dev_id_epj_2[id_stream].htod(nj_tot, stream[id_stream]);

            //cudaDeviceSynchronize();

            /*
            cudaDeviceSynchronize();
            dev_epi_2[id_stream].htod(ni_tot_reg);
            cudaDeviceSynchronize();
            cudaDeviceSynchronize();
            dev_id_epj_2[id_stream].htod(nj_tot);
            cudaDeviceSynchronize();
            */

            int nblocks  = ni_tot_reg / N_THREAD_GPU;
            int nthreads = N_THREAD_GPU;

            const float eps2 = FPGrav::eps * FPGrav::eps;

            /*
            float kernel_time_tmp = 0.0;
            cudaEvent_t beg, end;
            cudaEventCreate(&beg);
            cudaEventCreate(&end);
            cudaEventRecord(beg, 0);
            */
            ForceKernelIndex <<<nblocks, nthreads, 0, stream[id_stream]>>> (ij_disp_2[id_stream], dev_epi_2[id_stream], dev_epj, dev_id_epj_2[id_stream], dev_force_2[id_stream], eps2);
            //ForceKernelIndex <<<nblocks, nthreads>>> (ij_disp_2[id_stream], dev_epi_2[id_stream], dev_epj, dev_id_epj_2[id_stream], dev_force_2[id_stream], eps2, id_stream);
            /*
            cudaEventRecord(end, 0);
            cudaEventSynchronize(end);
            cudaEventElapsedTime( &kernel_time_tmp, beg, end);
            WTIME_KERNEL += kernel_time_tmp*1e-3;
            */

            //cudaDeviceSynchronize();
            //cudaStreamSynchronize(stream[id_stream]);
        }
        //cudaDeviceSynchronize();
        return 0;
    }
}

    #endif

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       FPGrav    *force[]){
    int id_walk_ar[N_STREAM+1];
    const int n_walk_ave = n_walk/N_STREAM;
    for(int i=0; i<N_STREAM+1; i++){ id_walk_ar[i] = n_walk_ave*i + std::min(i, n_walk%N_STREAM);}
    
/*
    cudaEvent_t beg_d2h, end_d2h;
    float time_d2h_tmp = 0.0;
    cudaEventCreate(&beg_d2h);
    cudaEventCreate(&end_d2h);
    cudaEventRecord(beg_d2h, 0);
*/
    //#pragma omp parallel for num_threads(N_STREAM)
    for(int id_stream=0; id_stream<N_STREAM; id_stream++){
        int id_walk_head = id_walk_ar[id_stream];
        int id_walk_end  = id_walk_ar[id_stream+1];

        int ni_tot = 0;
        for(int iw=id_walk_head; iw<id_walk_end; iw++){
            ni_tot += ni[iw];
        }
        dev_force_2[id_stream].dtoh(ni_tot, stream[id_stream]);
        cudaStreamSynchronize(stream[id_stream]);

        int n_cnt = 0;
        for(int iw=id_walk_head; iw<id_walk_end; iw++){
            for(int i=0; i<ni[iw]; i++, n_cnt++){
                force[iw][i].acc.x = dev_force_2[id_stream][n_cnt].accp.x;
                force[iw][i].acc.y = dev_force_2[id_stream][n_cnt].accp.y;
                force[iw][i].acc.z = dev_force_2[id_stream][n_cnt].accp.z;
                force[iw][i].pot   = dev_force_2[id_stream][n_cnt].accp.w;
            }
        }
    }
/*
    cudaEventRecord(end_d2h, 0);
    cudaEventSynchronize(end_d2h);
    cudaEventElapsedTime(&time_d2h_tmp, beg_d2h, end_d2h);
    WTIME_D2H += time_d2h_tmp*1e-3;
*/
    return 0;
}

#else
// original
PS::S32 DispatchKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const FPGrav ** epi,
                            const PS::S32 *  n_epi,
                            const PS::S32 ** id_epj,
                            const PS::S32 *  n_epj,
                            const PS::S32 ** id_spj,
                            const PS::S32 *  n_spj,
                            const FPGrav  * epj,
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
        /*
        if(dev_epj.size < n_epj_tot+n_spj_tot){
            dev_epj.free();
            dev_epj.allocate(n_epj_tot+n_spj_tot);
        }
        */
        for(PS::S32 i=0; i<n_epj_tot; i++){
            dev_epj[i].posm.x  = epj[i].pos.x;
            dev_epj[i].posm.y  = epj[i].pos.y;
            dev_epj[i].posm.z  = epj[i].pos.z;
            dev_epj[i].posm.w  = epj[i].mass;
        }
        for(PS::S32 i=0; i<n_spj_tot; i++){
            dev_epj[i+n_epj_tot].posm.x  = spj[i].pos.x;
            dev_epj[i+n_epj_tot].posm.y  = spj[i].pos.y;
            dev_epj[i+n_epj_tot].posm.z  = spj[i].pos.z;
            dev_epj[i+n_epj_tot].posm.w  = spj[i].mass;
        }
        dev_epj.htod(n_epj_tot+n_spj_tot);
        return 0;
    }
    else{
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
        int ni_tot_reg = ij_disp[n_walk].x;
        if(ni_tot_reg % N_THREAD_GPU){
            ni_tot_reg /= N_THREAD_GPU;
            ni_tot_reg++;
            ni_tot_reg *= N_THREAD_GPU;
        }
        int ni_tot = 0;
        int nj_tot = 0;

        for(int iw=0; iw<n_walk; iw++){
            for(int i=0; i<n_epi[iw]; i++){
                dev_epi[ni_tot].pos.x = epi[iw][i].pos.x;
                dev_epi[ni_tot].pos.y = epi[iw][i].pos.y;
                dev_epi[ni_tot].pos.z = epi[iw][i].pos.z;
                dev_epi[ni_tot].id_walk = iw;
                ni_tot++;
            }
            for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                dev_id_epj[nj_tot] = id_epj[iw][j];
            }
            for(int j=0; j<n_spj[iw]; j++, nj_tot++){
                dev_id_epj[nj_tot] = id_spj[iw][j]+n_epj_tot;
            }
        }
        for(int i=ni_tot; i<ni_tot_reg; i++){
            dev_epi[i].id_walk = n_walk;
        }
        dev_epi.htod(ni_tot_reg);
        dev_id_epj.htod(nj_tot);
        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        const float eps2 = FPGrav::eps * FPGrav::eps;

        /*
        cudaEvent_t beg, end;
        float kernel_time_tmp = 0.0;
        cudaEventCreate(&beg);
        cudaEventCreate(&end);
        cudaEventRecord(beg, 0);
        */

        ForceKernelIndex <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_id_epj, dev_force, eps2);
        /*
        cudaEventRecord(end, 0);
        cudaEventSynchronize(end);
        cudaEventElapsedTime( &kernel_time_tmp, beg, end);
        WTIME_KERNEL += kernel_time_tmp*1e-3;
        */
        return 0;
    }
}
PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       FPGrav    *force[])
{
    int ni_tot = 0;
    for(int k=0; k<n_walk; k++){
        ni_tot += ni[k];
    }
    dev_force.dtoh(ni_tot);
    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<ni[iw]; i++, n_cnt++){
            force[iw][i].acc.x = dev_force[n_cnt].accp.x;
            force[iw][i].acc.y = dev_force[n_cnt].accp.y;
            force[iw][i].acc.z = dev_force[n_cnt].accp.z;
            force[iw][i].pot   = dev_force[n_cnt].accp.w;
        }
    }
    return 0;
}
#endif

