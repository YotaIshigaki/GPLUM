
#include "class.hpp"
#include "force.hpp"

#define N_THREAD_GPU 64
const int N_WALK_LIMIT = 1000;
const int NI_LIMIT = N_WALK_LIMIT*1000;
const int NJ_LIMIT = N_WALK_LIMIT*10000;

#if 0
#  include <cutil.h>
#else
#  include <helper_cuda.h>
#  define CUDA_SAFE_CALL checkCudaErrors
#endif

class EpiGPU{
public:
    float2 pos[3];
    int id_walk;
};

class EpjGPU{
public:
    float mass;
    float2 pos[3];
};

class ForceGPU{
public:
    float2 acc[3];
    float2 pot;
};

static float2 float2_split(double x){
    const int shift = 20;
    float2 ret;
    x *= (1<<shift);
    double xi = (int)x;
    double xf = x - xi;
    ret.x = xi * (1./(1<<shift));
    ret.y = xf * (1./(1<<shift));
    return ret;
}

__device__ float2 float2_accum(float2 acc, float x){
    float tmp = acc.x + x;
    acc.y -= (tmp - acc.x) - x;
    acc.x = tmp;
    return acc;
}

__device__ float2 float2_regularize(float2 acc){
    float tmp = acc.x + acc.y;
    acc.y = acc.y -(tmp - acc.x);
    acc.x = tmp;
    return acc;
}

__global__ void ForceKernel(const EpiGPU * epi,
                            const int    * ni_disp,
                            const EpjGPU * epj, 
                            const int    * nj_disp,
                            ForceGPU     * force,
                            const float eps2){
    int id_i = blockDim.x * blockIdx.x + threadIdx.x;
    const EpiGPU & ip = epi[id_i];
    float2 poti;
    float2 acci[3];
    poti = acci[0] = acci[1] = acci[2] = make_float2(0.0, 0.0);
    const int j_head = nj_disp[ip.id_walk];
    const int j_tail = nj_disp[ip.id_walk+1];
    const int nj = j_tail - j_head;
    for(int j=j_head; j<j_tail; j++){
        //int j = j_head + threadIdx.x%nj;
        //for(int jtmp=0; jtmp<nj; jtmp++){
        EpjGPU jp = epj[j];
        const float dx = (jp.pos[0].x - ip.pos[0].x) + (jp.pos[0].y - ip.pos[0].y);
        const float dy = (jp.pos[1].x - ip.pos[1].x) + (jp.pos[1].y - ip.pos[1].y);
        const float dz = (jp.pos[2].x - ip.pos[2].x) + (jp.pos[2].y - ip.pos[2].y); // 9op
        const float r2 = ((eps2 + dx*dx) + dy*dy) + dz*dz; // 15op
        const float r_inv = rsqrtf(r2); // 15op + 1rsqrt
        const float pij = jp.mass * r_inv * (r2 > eps2); 
        const float r2_inv = r_inv * r_inv; 
        const float pij_r3_inv = pij * r2_inv; 
        const float ax = pij_r3_inv * dx; 
        const float ay = pij_r3_inv * dy;
        const float az = pij_r3_inv * dz; //21op + 1rsqrt
        poti = float2_accum(poti, pij);
        acci[0] = float2_accum(acci[0], ax);
        acci[1] = float2_accum(acci[1], ay);
        acci[2] = float2_accum(acci[2], az); // 33op + 1rsqrt
        //j = (j+1)>=j_tail ? j_head : j+1;
    }
    poti = float2_regularize(poti);
    acci[0] = float2_regularize(acci[0]);
    acci[1] = float2_regularize(acci[1]);
    acci[2] = float2_regularize(acci[2]);
    force[id_i].pot = poti;
    force[id_i].acc[0] = acci[0];
    force[id_i].acc[1] = acci[1];
    force[id_i].acc[2] = acci[2];
}

static ForceGPU * force_d;
static ForceGPU * force_h;
static EpiGPU * epi_d;
static EpiGPU * epi_h;
static EpjGPU * epj_d;
static EpjGPU * epj_h;
static int * ni_disp_d;
static int * ni_disp_h;
static int * nj_disp_d;
static int * nj_disp_h;

int DispatchKernelWithSP(const PS::S32 tag,
                         const int    n_walk,
                         const EPIGrav ** epi,
                         const int  *  n_epi,
                         const EPJGrav ** epj,
                         const int  *  n_epj,
                         const PS::SPJMonopole ** spj,
                         const int  *  n_spj){

    static bool first = true;
    assert(n_walk <= N_WALK_LIMIT);
    if(first){
        CUDA_SAFE_CALL( cudaMalloc(     (void**)&ni_disp_d,  (N_WALK_LIMIT+1)*sizeof(int) ) );
        CUDA_SAFE_CALL( cudaMalloc(     (void**)&nj_disp_d,  (N_WALK_LIMIT+1)*sizeof(int) ) );
        CUDA_SAFE_CALL( cudaMallocHost( (void**)&ni_disp_h,  (N_WALK_LIMIT+1)*sizeof(int) ) );
        CUDA_SAFE_CALL( cudaMallocHost( (void**)&nj_disp_h,  (N_WALK_LIMIT+1)*sizeof(int) ) );
        CUDA_SAFE_CALL( cudaMalloc( (void**)&epi_d,       NI_LIMIT*sizeof(EpiGPU) ) );
        CUDA_SAFE_CALL( cudaMalloc( (void**)&epj_d,       NJ_LIMIT*sizeof(EpjGPU) ) );
        CUDA_SAFE_CALL( cudaMalloc( (void**)&force_d,     NI_LIMIT*sizeof(ForceGPU) ) );
        CUDA_SAFE_CALL( cudaMallocHost( (void**)&epi_h,   NI_LIMIT*sizeof(EpiGPU) ) );
        CUDA_SAFE_CALL( cudaMallocHost( (void**)&epj_h,   NJ_LIMIT*sizeof(EpjGPU) ) );
        CUDA_SAFE_CALL( cudaMallocHost( (void**)&force_h, NI_LIMIT*sizeof(ForceGPU) ) );
        first = false;
    }
    const float eps2 = EPIGrav::eps * EPIGrav::eps;
    //CUDA_SAFE_CALL( cudaMalloc(     (void**)&ni_disp_d,     (n_walk+1)*sizeof(int) ) );
    //CUDA_SAFE_CALL( cudaMalloc(     (void**)&nj_disp_d,     (n_walk+1)*sizeof(int) ) );
    //CUDA_SAFE_CALL( cudaMallocHost( (void**)&ni_disp_h,     (n_walk+1)*sizeof(int) ) );
    //CUDA_SAFE_CALL( cudaMallocHost( (void**)&nj_disp_h,     (n_walk+1)*sizeof(int) ) );

    ni_disp_h[0] = nj_disp_h[0] = 0;
    for(int i=0; i<n_walk; i++){
        ni_disp_h[i+1] = ni_disp_h[i] + n_epi[i];
        nj_disp_h[i+1] = nj_disp_h[i] + n_epj[i] + n_spj[i];
    }
    int ni_tot = ni_disp_h[n_walk];
    const int ni_tot_reg = ni_disp_h[n_walk] + ( (ni_tot%N_THREAD_GPU != 0) ? (N_THREAD_GPU - (ni_tot%N_THREAD_GPU)) : 0);
    //CUDA_SAFE_CALL( cudaMalloc( (void**)&epi_d, ni_tot_reg*sizeof(EpiGPU) ) );
    //CUDA_SAFE_CALL( cudaMalloc( (void**)&epj_d, nj_disp_h[n_walk]*sizeof(EpjGPU) ) );
    //CUDA_SAFE_CALL( cudaMalloc( (void**)&force_d,     ni_tot_reg*sizeof(ForceGPU) ) );
    //CUDA_SAFE_CALL( cudaMallocHost( (void**)&epi_h, ni_tot_reg*sizeof(EpiGPU) ) );
    //CUDA_SAFE_CALL( cudaMallocHost( (void**)&epj_h, nj_disp_h[n_walk]*sizeof(EpjGPU) ) );
    //CUDA_SAFE_CALL( cudaMallocHost( (void**)&force_h, ni_tot_reg*sizeof(ForceGPU) ) );
    assert(ni_tot_reg <= NI_LIMIT);
    assert(nj_disp_h[n_walk] <= NJ_LIMIT);
    ni_tot = 0;
    int nj_tot = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int ip=0; ip<n_epi[iw]; ip++){
            epi_h[ni_tot].pos[0]  = float2_split(epi[iw][ip].pos.x);
            epi_h[ni_tot].pos[1]  = float2_split(epi[iw][ip].pos.y);
            epi_h[ni_tot].pos[2]  = float2_split(epi[iw][ip].pos.z);
            epi_h[ni_tot].id_walk = iw;
            force_h[ni_tot].acc[0] = force_h[ni_tot].acc[1] 
                = force_h[ni_tot].acc[2] = force_h[ni_tot].pot = make_float2(0.0, 0.0);
            ni_tot++;
        }
        for(int jp=0; jp<n_epj[iw]; jp++){
            epj_h[nj_tot].mass    = epj[iw][jp].mass;
            epj_h[nj_tot].pos[0]  = float2_split(epj[iw][jp].pos.x);
            epj_h[nj_tot].pos[1]  = float2_split(epj[iw][jp].pos.y);
            epj_h[nj_tot].pos[2]  = float2_split(epj[iw][jp].pos.z);
            nj_tot++;
        }
        for(int jp=0; jp<n_spj[iw]; jp++){
            epj_h[nj_tot].mass    = spj[iw][jp].getCharge();
            epj_h[nj_tot].pos[0]  = float2_split(spj[iw][jp].getPos().x);
            epj_h[nj_tot].pos[1]  = float2_split(spj[iw][jp].getPos().y);
            epj_h[nj_tot].pos[2]  = float2_split(spj[iw][jp].getPos().z);
            nj_tot++;
        }
    }
    for(int ip=ni_tot; ip<ni_tot_reg; ip++){
        epi_h[ni_tot].pos[0]  = epi_h[ni_tot].pos[1]  = epi_h[ni_tot].pos[2]  = make_float2(0.0, 0.0);
        epi_h[ni_tot].id_walk = 0;
        force_h[ni_tot].acc[0] = force_h[ni_tot].acc[1] 
            = force_h[ni_tot].acc[2] = force_h[ni_tot].pot = make_float2(0.0, 0.0);
    }
    CUDA_SAFE_CALL( cudaMemcpy(epi_d, epi_h, ni_tot_reg*sizeof(EpiGPU), cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(epj_d, epj_h, nj_tot*sizeof(EpjGPU), cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(ni_disp_d, ni_disp_h, (n_walk+1)*sizeof(int), cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(nj_disp_d, nj_disp_h, (n_walk+1)*sizeof(int), cudaMemcpyHostToDevice) );
    const int n_grid = ni_tot_reg/N_THREAD_GPU + ((ni_tot_reg%N_THREAD_GPU == 0) ? 0 : 1);
    dim3 size_grid(n_grid, 1, 1);
    dim3 size_thread(N_THREAD_GPU, 1, 1);
    ForceKernel<<<size_grid, size_thread>>> (epi_d, ni_disp_d, epj_d, nj_disp_d, force_d, float(eps2));

    return 0;
}

int RetrieveKernel(const PS::S32 tag,
                   const PS::S32    n_walk,
                   const PS::S32 *  ni,
                   ForceGrav     ** force){

    int ni_tot = 0;
    for(int i=0; i<n_walk; i++){
        ni_tot += ni[i];
    }
    //const int ni_tot_reg = ni_disp_h[n_walk] + ( (ni_tot%N_THREAD_GPU != 0) ? (N_THREAD_GPU - (ni_tot%N_THREAD_GPU)) : 0);
    //CUDA_SAFE_CALL( cudaMemcpy(force_h, force_d,      ni_tot_reg*sizeof(ForceGPU), cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy(force_h, force_d,      ni_tot*sizeof(ForceGPU), cudaMemcpyDeviceToHost) );
    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int ip=0; ip<ni[iw]; ip++){
            force[iw][ip].acc.x = (double)force_h[n_cnt].acc[0].x + (double)force_h[n_cnt].acc[0].y;
            force[iw][ip].acc.y = (double)force_h[n_cnt].acc[1].x + (double)force_h[n_cnt].acc[1].y;
            force[iw][ip].acc.z = (double)force_h[n_cnt].acc[2].x + (double)force_h[n_cnt].acc[2].y;
            force[iw][ip].pot   = (double)force_h[n_cnt].pot.x    + (double)force_h[n_cnt].pot.y;
            force[iw][ip].pot *= -1.0;
            n_cnt++;
        }
    }

    /*
    CUDA_SAFE_CALL( cudaFreeHost(force_h) );
    CUDA_SAFE_CALL( cudaFree(force_d) );
    CUDA_SAFE_CALL( cudaFree(epi_d) );
    CUDA_SAFE_CALL( cudaFreeHost(epi_h) );
    CUDA_SAFE_CALL( cudaFree(epj_d) );
    CUDA_SAFE_CALL( cudaFreeHost(epj_h) );
    CUDA_SAFE_CALL( cudaFree(ni_disp_d) );
    CUDA_SAFE_CALL( cudaFreeHost(ni_disp_h) );
    CUDA_SAFE_CALL( cudaFree(nj_disp_d) );
    CUDA_SAFE_CALL( cudaFreeHost(nj_disp_h) );
    */
    return 0;
}
