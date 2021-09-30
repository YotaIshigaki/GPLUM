#include<particle_simulator.hpp>
#include "cuda_pointer.h"
#include "force_gpu.h"

enum{
    N_THREAD_GPU = 32,
    N_WALK_LIMIT = 1000,
    NI_LIMIT     = N_WALK_LIMIT*1000,
    NJ_LIMIT     = N_WALK_LIMIT*10000,
};

static cudaPointer<int2>    ij_disp_cu;
static cudaPointer<int>     id_epj_cu;

static const float SMTH = 1.2;
static const float SUP_RAD = 2.5;
static const float MY_PI = 3.14159265358979323846;
static const int DIM = 3.0;
static const float MY_C_CFL = 0.3;

inline __device__ float calcW(float3 dr, float h){
    float H = SUP_RAD*h;
    float r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
    float s = sqrt(r2) / H;
    float a = (1.0-s) > 0.0 ? (1.0-s) : 0.0;
    float b = (0.5-s) > 0.0 ? (0.5-s) : 0.0;
    a = a*a*a;
    b = b*b*b;
    float r_value = a - 4.0*b;
    r_value *= 16.0 / MY_PI / (H * H * H);
    return r_value;
}

inline __device__ float3 calcGradW(float3 dr, float h){
    float H = SUP_RAD*h;
    float r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
    float r_abs = sqrt(r2);
    float s = r_abs / H;
    float a = (1.0-s) > 0.0 ? (1.0-s) : 0.0;
    float b = (0.5-s) > 0.0 ? (0.5-s) : 0.0;
    a = a*a;
    b = b*b;
    float r_value = -3.0*a + 12.0*b;
    r_value *= 16.0 / MY_PI / (H * H * H);
    float3 ret;
    ret.x = dr.x*r_value / (r_abs * H + 1.0e-6 * h);
    ret.y = dr.y*r_value / (r_abs * H + 1.0e-6 * h);
    ret.z = dr.z*r_value / (r_abs * H + 1.0e-6 * h);
    return ret;
}

///////////
// HYDRO //
struct HydroEpi{
    float3 pos;
    float3 vel;
    float  smth;
    float  dens;
    float  pres;
    float  snds;
    float  Bal;
    int    id_walk;
};

struct HydroEpj{
    float4 posm;
    float3 vel;
    float  dens;
    float  pres;
    float  smth;
    float  snds;
    float  Bal;
};
struct HydroRes{
    float3 acc;
    float  eng_dot;
    float  dt;
};

static cudaPointer<HydroEpi> hydro_epi_cu;
static cudaPointer<HydroEpj> hydro_epj_cu;
static cudaPointer<HydroRes> hydro_res_cu;

inline __device__ void dev_calc_hydro(const float3 ipos,
                                      const float3 ivel,
                                      const float  ipres,
                                      const float  idens,
                                      const float  isnds,
                                      const float  iBal,
                                      const float  ismth,
                                      const float4 jposm,
                                      const float3 jvel,
                                      const float  jpres,
                                      const float  jdens,
                                      const float  jsnds,
                                      const float  jBal,
                                      const float  jsmth,
                                      float3 & iacc,
                                      float  & ieng_dot,
                                      float & v_sig_max){
    float mj = jposm.w;
    float3 dr;
    dr.x = -(jposm.x - ipos.x);
    dr.y = -(jposm.y - ipos.y);
    dr.z = -(jposm.z - ipos.z);
    float r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
    float r_abs = sqrt(r2);

    float3 dv;
    dv.x = -(jvel.x - ivel.x);
    dv.y = -(jvel.y - ivel.y);
    dv.z = -(jvel.z - ivel.z);

    float dvdr = dv.x*dr.x + dv.y*dr.y + dv.z*dr.z;
    float w_ij = (dvdr < 0.0) ? dvdr/r_abs : 0.0;
    float v_sig = isnds + jsnds - 3.0 * w_ij;
    v_sig_max = (v_sig > v_sig_max) ? v_sig : v_sig_max;
    float AV = -0.5*v_sig*w_ij / (0.5 * (idens + jdens)) * 0.5 * (iBal + jBal);

    float3 igradw = calcGradW(dr, ismth);
    float3 jgradw = calcGradW(dr, jsmth);
    float3 gradw;
    gradw.x = (igradw.x+jgradw.x)*0.5;
    gradw.y = (igradw.y+jgradw.y)*0.5;
    gradw.z = (igradw.z+jgradw.z)*0.5;

    float coef_acc = mj * (ipres / (idens * idens) + jpres / (jdens * jdens) + AV);
    iacc.x -= coef_acc * gradw.x;
    iacc.y -= coef_acc * gradw.y;
    iacc.z -= coef_acc * gradw.z;
    float dvgradw = dv.x*gradw.x + dv.y*gradw.y + dv.z*gradw.z;
    ieng_dot += mj * (ipres / (idens * idens) + 0.5 * AV) * dvgradw;

}

__global__ void HydroKernel(const int2   * ij_disp,
                            const HydroEpi * epi,
                            const HydroEpj * epj,
                            const int * id_epj, 
                            HydroRes     * res){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const float3 ipos = epi[tid].pos;
    const float3 ivel = epi[tid].vel;
    const float ipres = epi[tid].pres;
    const float idens  = epi[tid].dens;
    const float isnds  = epi[tid].snds;
    const float iBal   = epi[tid].Bal;
    const float ismth  = epi[tid].smth;
    const int j_head = ij_disp[epi[tid].id_walk  ].y;
    const int j_tail = ij_disp[epi[tid].id_walk+1].y;
    float  ieng_dot = 0.0;
    float3 iacc;
    iacc.x = iacc.y = iacc.z = 0.0;
    float v_sig_max = -1.0;
    for(int j=j_head; j<j_tail; j++){
        int adr      = id_epj[j];
        float4 jposm = epj[adr].posm;
        float3 jvel  = epj[adr].vel;
        float jpres  = epj[adr].pres;
        float jdens  = epj[adr].dens;
        float jsnds  = epj[adr].snds;
        float jBal   = epj[adr].Bal;
        float jsmth  = epj[adr].smth;
        dev_calc_hydro(ipos,  ivel, ipres, idens, isnds, iBal, ismth,
                       jposm, jvel, jpres, jdens, jsnds, jBal, jsmth,
                       iacc,  ieng_dot, v_sig_max);
    }
    HydroRes res_tmp;
    res_tmp.acc.x = iacc.x;
    res_tmp.acc.y = iacc.y;
    res_tmp.acc.z = iacc.z;
    res_tmp.eng_dot = ieng_dot;
    res_tmp.dt = MY_C_CFL * 2.0 * ismth / v_sig_max;
    res[tid] = res_tmp;
}

PS::S32 CalcHydroDispatch(const PS::S32 tag,
                         const PS::S32 n_walk,
                         const EPI::Hydro ** epi,
                         const PS::S32 *  n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32 *  n_epj,
                         const EPJ::Hydro * epj,
                         const PS::S32 n_epj_tot,
                         const bool send_flag){
    static bool init_call = true;
    if(init_call){
        hydro_epi_cu.allocate(NI_LIMIT);
        hydro_epj_cu.allocate(NJ_LIMIT);
        id_epj_cu.allocate(NJ_LIMIT);
        hydro_res_cu.allocate(NI_LIMIT);
        ij_disp_cu.allocate(N_WALK_LIMIT+2);
        init_call = false;
    }
    if(send_flag==true){
        for(PS::S32 i=0; i<n_epj_tot; i++){
            hydro_epj_cu[i].posm.x  = epj[i].pos.x;
            hydro_epj_cu[i].posm.y  = epj[i].pos.y;
            hydro_epj_cu[i].posm.z  = epj[i].pos.z;
            hydro_epj_cu[i].posm.w  = epj[i].mass;
            hydro_epj_cu[i].vel.x  = epj[i].vel.x;
            hydro_epj_cu[i].vel.y  = epj[i].vel.y;
            hydro_epj_cu[i].vel.z  = epj[i].vel.z;
            hydro_epj_cu[i].dens   = epj[i].dens;
            hydro_epj_cu[i].pres   = epj[i].pres;
            hydro_epj_cu[i].smth   = epj[i].smth;
            hydro_epj_cu[i].snds   = epj[i].snds;
            hydro_epj_cu[i].Bal    = epj[i].Bal;
        }
        hydro_epj_cu.htod(n_epj_tot);
        return 0;
    }
    else{
        ij_disp_cu[0].x = 0;
        ij_disp_cu[0].y = 0;
        for(int k=0; k<n_walk; k++){
            ij_disp_cu[k+1].x = ij_disp_cu[k].x + n_epi[k];
            ij_disp_cu[k+1].y = ij_disp_cu[k].y + n_epj[k];
        }
        ij_disp_cu[n_walk+1] = ij_disp_cu[n_walk];
        assert(ij_disp_cu[n_walk].x < NI_LIMIT);
        assert(ij_disp_cu[n_walk].y < NJ_LIMIT);
        ij_disp_cu.htod(n_walk + 2);
        int ni_tot_reg = ij_disp_cu[n_walk].x;
        if(ni_tot_reg % N_THREAD_GPU){
            ni_tot_reg /= N_THREAD_GPU;
            ni_tot_reg++;
            ni_tot_reg *= N_THREAD_GPU;
        }
        int ni_tot = 0;
        int nj_tot = 0;
        for(int iw=0; iw<n_walk; iw++){
            for(int i=0; i<n_epi[iw]; i++, ni_tot++){
                hydro_epi_cu[ni_tot].pos.x = epi[iw][i].pos.x;
                hydro_epi_cu[ni_tot].pos.y = epi[iw][i].pos.y;
                hydro_epi_cu[ni_tot].pos.z = epi[iw][i].pos.z;
                hydro_epi_cu[ni_tot].vel.x = epi[iw][i].vel.x;
                hydro_epi_cu[ni_tot].vel.y = epi[iw][i].vel.y;
                hydro_epi_cu[ni_tot].vel.z = epi[iw][i].vel.z;
                hydro_epi_cu[ni_tot].smth  = epi[iw][i].smth;
                hydro_epi_cu[ni_tot].dens  = epi[iw][i].dens;
                hydro_epi_cu[ni_tot].pres  = epi[iw][i].pres;
                hydro_epi_cu[ni_tot].snds  = epi[iw][i].snds;
                hydro_epi_cu[ni_tot].Bal   = epi[iw][i].Bal;
                hydro_epi_cu[ni_tot].id_walk = iw;
            }
            for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                id_epj_cu[nj_tot] = id_epj[iw][j];
            }
        }
        for(int i=ni_tot; i<ni_tot_reg; i++){
            hydro_epi_cu[i].id_walk = n_walk;
        }
        hydro_epi_cu.htod(ni_tot_reg);
        id_epj_cu.htod(nj_tot);
        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        HydroKernel <<<nblocks, nthreads>>> (ij_disp_cu,  hydro_epi_cu, 
                                             hydro_epj_cu, id_epj_cu, 
                                             hydro_res_cu);
        return 0;
    }
}


//////////
// DRVT //
struct DrvtEpi{
    float3 pos;
    float3 vel;
    float  smth;
    float    dens;
    int    id_walk;
};

struct DrvtEpj{
    float4 posm;
    float3 vel;
};

struct DrvtRes{
    float div_v;
    float3 rot_v;
};

static cudaPointer<DrvtEpi> drvt_epi_cu;
static cudaPointer<DrvtEpj> drvt_epj_cu;
static cudaPointer<DrvtRes> drvt_res_cu;



inline __device__ void dev_calc_drvt(const float3 ipos,
                                     const float3 ivel,
                                     const float  ismth,
                                     const float4 jposm,
                                     const float3 jvel,
                                     float  & div_v,
                                     float3 & rot_v){
    float mj = jposm.w;
    float3 dr;
    dr.x = -(jposm.x - ipos.x);
    dr.y = -(jposm.y - ipos.y);
    dr.z = -(jposm.z - ipos.z);

    float3 dv;
    dv.x = -(jvel.x - ivel.x);
    dv.y = -(jvel.y - ivel.y);
    dv.z = -(jvel.z - ivel.z);

    float3 gradw = calcGradW(dr, ismth);
    div_v   += mj * (dv.x*gradw.x + dv.y*gradw.y + dv.z*gradw.z);
    rot_v.x += mj*(dv.y*gradw.z - dv.z*gradw.y);
    rot_v.y += mj*(dv.z*gradw.x - dv.x*gradw.z);
    rot_v.z += mj*(dv.x*gradw.y - dv.y*gradw.x);
}

__global__ void DrvtKernel(const int2   * ij_disp,
                           const DrvtEpi * epi,
                           const DrvtEpj * epj,
                           const int * id_epj, 
                           DrvtRes     * res){

    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const float3 ipos = epi[tid].pos;
    const float3 ivel = epi[tid].vel;
    const float ismth  = epi[tid].smth;
    const float idens  = epi[tid].dens;
    const int j_head = ij_disp[epi[tid].id_walk  ].y;
    const int j_tail = ij_disp[epi[tid].id_walk+1].y;
    float  div_v = 0.0;
    float3 rot_v;
    rot_v.x = rot_v.y = rot_v.z = 0.0;
    for(int j=j_head; j<j_tail; j++){
        int adr      = id_epj[j];
        float4 jposm = epj[adr].posm;
        float3 jvel  = epj[adr].vel;
        dev_calc_drvt(ipos, ivel, ismth, jposm, jvel, div_v, rot_v);
    }
    DrvtRes res_tmp;
    res_tmp.div_v   = div_v   / idens;
    res_tmp.rot_v.x = rot_v.x / idens;
    res_tmp.rot_v.y = rot_v.y / idens;
    res_tmp.rot_v.z = rot_v.z / idens;
    res[tid] = res_tmp;
}


PS::S32 CalcDrvtDispatch(const PS::S32 tag,
                         const PS::S32 n_walk,
                         const EPI::Drvt ** epi,
                         const PS::S32 *  n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32 *  n_epj,
                         const EPJ::Drvt * epj,
                         const PS::S32 n_epj_tot,
                         const bool send_flag){
    static bool init_call = true;
    if(init_call){
        drvt_epi_cu.allocate(NI_LIMIT);
        drvt_epj_cu.allocate(NJ_LIMIT);
        id_epj_cu.allocate(NJ_LIMIT);
        drvt_res_cu.allocate(NI_LIMIT);
        ij_disp_cu.allocate(N_WALK_LIMIT+2);
        init_call = false;
    }
    if(send_flag==true){
        for(PS::S32 i=0; i<n_epj_tot; i++){
            drvt_epj_cu[i].posm.x  = epj[i].pos.x;
            drvt_epj_cu[i].posm.y  = epj[i].pos.y;
            drvt_epj_cu[i].posm.z  = epj[i].pos.z;
            drvt_epj_cu[i].posm.w  = epj[i].mass;
        }
        drvt_epj_cu.htod(n_epj_tot);
        return 0;
    }
    else{
        ij_disp_cu[0].x = 0;
        ij_disp_cu[0].y = 0;
        for(int k=0; k<n_walk; k++){
            ij_disp_cu[k+1].x = ij_disp_cu[k].x + n_epi[k];
            ij_disp_cu[k+1].y = ij_disp_cu[k].y + n_epj[k];
        }
        ij_disp_cu[n_walk+1] = ij_disp_cu[n_walk];
        assert(ij_disp_cu[n_walk].x < NI_LIMIT);
        assert(ij_disp_cu[n_walk].y < NJ_LIMIT);
        ij_disp_cu.htod(n_walk + 2);
        int ni_tot_reg = ij_disp_cu[n_walk].x;
        if(ni_tot_reg % N_THREAD_GPU){
            ni_tot_reg /= N_THREAD_GPU;
            ni_tot_reg++;
            ni_tot_reg *= N_THREAD_GPU;
        }
        int ni_tot = 0;
        int nj_tot = 0;
        for(int iw=0; iw<n_walk; iw++){
            for(int i=0; i<n_epi[iw]; i++, ni_tot++){
                drvt_epi_cu[ni_tot].pos.x = epi[iw][i].pos.x;
                drvt_epi_cu[ni_tot].pos.y = epi[iw][i].pos.y;
                drvt_epi_cu[ni_tot].pos.z = epi[iw][i].pos.z;
                drvt_epi_cu[ni_tot].vel.x = epi[iw][i].vel.x;
                drvt_epi_cu[ni_tot].vel.y = epi[iw][i].vel.y;
                drvt_epi_cu[ni_tot].vel.z = epi[iw][i].vel.z;
                drvt_epi_cu[ni_tot].smth  = epi[iw][i].smth;
                drvt_epi_cu[ni_tot].dens  = epi[iw][i].dens;
                drvt_epi_cu[ni_tot].id_walk = iw;
            }
            for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                id_epj_cu[nj_tot] = id_epj[iw][j];
            }
        }
        for(int i=ni_tot; i<ni_tot_reg; i++){
            drvt_epi_cu[i].id_walk = n_walk;
        }
        drvt_epi_cu.htod(ni_tot_reg);
        id_epj_cu.htod(nj_tot);
        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        DrvtKernel <<<nblocks, nthreads>>> (ij_disp_cu,  drvt_epi_cu, 
                                            drvt_epj_cu, id_epj_cu, 
                                            drvt_res_cu);
        return 0;
    }
}

//////////
// DENS //

struct DensEpi{
    float4 posm;
    float  smth;
    int    id_walk;
};

struct DensEpj{
    float4 posm;
};

struct DensRes{
    float smth;
    float dens;
};

static cudaPointer<DensEpi> dens_epi_cu;
static cudaPointer<DensEpj> dens_epj_cu;
static cudaPointer<DensRes> dens_res_cu;

inline __device__ float dev_calc_dens(const float4 iposm,
                                      const float   ismth,
                                      const float4 & jposm,
                                      float dens){
    float mj = jposm.w;
    float3 dr;
    dr.x = -(jposm.x - iposm.x);
    dr.y = -(jposm.y - iposm.y);
    dr.z = -(jposm.z - iposm.z);
    dens += mj * calcW(dr, ismth);
    return dens;
}

__global__ void DensKernel(const int2   * ij_disp,
                           const DensEpi * epi,
                           const DensEpj * epj,
                           const int * id_epj, 
                           DensRes     * res){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const float4 iposm = epi[tid].posm;
    const float ismth = epi[tid].smth;
    const int j_head = ij_disp[epi[tid].id_walk  ].y;
    const int j_tail = ij_disp[epi[tid].id_walk+1].y;
    float dens = 0.0;
    for(int j=j_head; j<j_tail; j++){
        int adr      = id_epj[j];
        float4 jposm = epj[adr].posm;
        dens = dev_calc_dens(iposm, ismth, jposm, dens);
    }
    DensRes res_tmp;
    res_tmp.dens = dens;
    float mi = iposm.w;
    res_tmp.smth = SMTH * pow(mi / res_tmp.dens, 1.0f/(float)(DIM));
    res[tid] = res_tmp;
}

PS::S32 CalcDensDispatch(const PS::S32 tag,
                         const PS::S32 n_walk,
                         const EPI::Dens ** epi,
                         const PS::S32 *  n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32 *  n_epj,
                         const EPJ::Dens * epj,
                         const PS::S32 n_epj_tot,
                         const bool send_flag){
    static bool init_call = true;
    if(init_call){
        dens_epi_cu.allocate(NI_LIMIT);
        dens_epj_cu.allocate(NJ_LIMIT);
        id_epj_cu.allocate(NJ_LIMIT);
        dens_res_cu.allocate(NI_LIMIT);
        ij_disp_cu.allocate(N_WALK_LIMIT+2);
        init_call = false;
    }
    if(send_flag==true){
        for(PS::S32 i=0; i<n_epj_tot; i++){
            dens_epj_cu[i].posm.x  = epj[i].pos.x;
            dens_epj_cu[i].posm.y  = epj[i].pos.y;
            dens_epj_cu[i].posm.z  = epj[i].pos.z;
            dens_epj_cu[i].posm.w  = epj[i].mass;
        }
        dens_epj_cu.htod(n_epj_tot);
        return 0;
    }
    else{
        ij_disp_cu[0].x = 0;
        ij_disp_cu[0].y = 0;
        for(int k=0; k<n_walk; k++){
            ij_disp_cu[k+1].x = ij_disp_cu[k].x + n_epi[k];
            ij_disp_cu[k+1].y = ij_disp_cu[k].y + n_epj[k];
        }
        ij_disp_cu[n_walk+1] = ij_disp_cu[n_walk];
        assert(ij_disp_cu[n_walk].x < NI_LIMIT);
        assert(ij_disp_cu[n_walk].y < NJ_LIMIT);
        ij_disp_cu.htod(n_walk + 2);
        int ni_tot_reg = ij_disp_cu[n_walk].x;
        if(ni_tot_reg % N_THREAD_GPU){
            ni_tot_reg /= N_THREAD_GPU;
            ni_tot_reg++;
            ni_tot_reg *= N_THREAD_GPU;
        }
        int ni_tot = 0;
        int nj_tot = 0;
        for(int iw=0; iw<n_walk; iw++){
            for(int i=0; i<n_epi[iw]; i++, ni_tot++){
                dens_epi_cu[ni_tot].posm.x = epi[iw][i].pos.x;
                dens_epi_cu[ni_tot].posm.y = epi[iw][i].pos.y;
                dens_epi_cu[ni_tot].posm.z = epi[iw][i].pos.z;
                dens_epi_cu[ni_tot].posm.w = epi[iw][i].mass;
                dens_epi_cu[ni_tot].smth  = epi[iw][i].smth;
                dens_epi_cu[ni_tot].id_walk = iw;
            }
            for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                id_epj_cu[nj_tot] = id_epj[iw][j];
            }
        }
        for(int i=ni_tot; i<ni_tot_reg; i++){
            dens_epi_cu[i].id_walk = n_walk;
        }
        dens_epi_cu.htod(ni_tot_reg);
        id_epj_cu.htod(nj_tot);
        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        //std::cerr<<"before dispatch"<<std::endl;
        DensKernel <<<nblocks, nthreads>>> (ij_disp_cu,  dens_epi_cu, 
                                            dens_epj_cu, id_epj_cu, 
                                            dens_res_cu);
        //std::cerr<<"after dispatch"<<std::endl;
        return 0;
    }
}

/////////////
// GRAVITY //
struct GravEpi{
    float3 pos;
    int    id_walk;
};

struct GravEpj{
    float4 posm;
};

struct GravRes{
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

static cudaPointer<GravEpi> grav_epi_cu;
static cudaPointer<GravEpj> grav_epj_cu;
static cudaPointer<GravRes> grav_res_cu;
__global__ void GravKernel(const int2   * ij_disp,
                           const GravEpi * epi,
                           const GravEpj * epj,
                           const int * id_epj, 
                           GravRes     * force,
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

PS::S32 CalcGravDispatch(const PS::S32 tag,
                         const PS::S32 n_walk,
                         const EPI::Grav ** epi,
                         const PS::S32 *  n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32 *  n_epj,
                         const PS::S32 ** id_spj,
                         const PS::S32 *  n_spj,
                         const EPJ::Grav  * epj,
                         const PS::S32 n_epj_tot,
                         const PS::SPJMonopole * spj,
                         const PS::S32 n_spj_tot,
                         const bool send_flag){
    assert(n_walk <= N_WALK_LIMIT);
    static bool init_call = true;
    if(init_call){
        grav_epi_cu.allocate(NI_LIMIT);
        grav_res_cu.allocate(NI_LIMIT);
        grav_epj_cu.allocate(NJ_LIMIT);
        ij_disp_cu .allocate(N_WALK_LIMIT+2);
        id_epj_cu  .allocate(NJ_LIMIT);
        init_call = false;
    }
    if(send_flag==true){
        for(PS::S32 i=0; i<n_epj_tot; i++){
            grav_epj_cu[i].posm.x  = epj[i].pos.x;
            grav_epj_cu[i].posm.y  = epj[i].pos.y;
            grav_epj_cu[i].posm.z  = epj[i].pos.z;
            grav_epj_cu[i].posm.w  = epj[i].mass;
        }
        for(PS::S32 i=0; i<n_spj_tot; i++){
            grav_epj_cu[i+n_epj_tot].posm.x  = spj[i].pos.x;
            grav_epj_cu[i+n_epj_tot].posm.y  = spj[i].pos.y;
            grav_epj_cu[i+n_epj_tot].posm.z  = spj[i].pos.z;
            grav_epj_cu[i+n_epj_tot].posm.w  = spj[i].mass;
        }
        grav_epj_cu.htod(n_epj_tot+n_spj_tot);
        return 0;
    }
    else{
        ij_disp_cu[0].x = 0;
        ij_disp_cu[0].y = 0;
        for(int k=0; k<n_walk; k++){
            ij_disp_cu[k+1].x = ij_disp_cu[k].x + n_epi[k];
            ij_disp_cu[k+1].y = ij_disp_cu[k].y + (n_epj[k] + n_spj[k]);
        }
        ij_disp_cu[n_walk+1] = ij_disp_cu[n_walk];
        assert(ij_disp_cu[n_walk].x < NI_LIMIT);
        assert(ij_disp_cu[n_walk].y < NJ_LIMIT);
        ij_disp_cu.htod(n_walk + 2);
        int ni_tot_reg = ij_disp_cu[n_walk].x;
        if(ni_tot_reg % N_THREAD_GPU){
            ni_tot_reg /= N_THREAD_GPU;
            ni_tot_reg++;
            ni_tot_reg *= N_THREAD_GPU;
        }
        int ni_tot = 0;
        int nj_tot = 0;
        for(int iw=0; iw<n_walk; iw++){
            for(int i=0; i<n_epi[iw]; i++){
                grav_epi_cu[ni_tot].pos.x = epi[iw][i].pos.x;
                grav_epi_cu[ni_tot].pos.y = epi[iw][i].pos.y;
                grav_epi_cu[ni_tot].pos.z = epi[iw][i].pos.z;
                grav_epi_cu[ni_tot].id_walk = iw;
                ni_tot++;
            }
            for(int j=0; j<n_epj[iw]; j++, nj_tot++){
                id_epj_cu[nj_tot] = id_epj[iw][j];
            }
            for(int j=0; j<n_spj[iw]; j++, nj_tot++){
                id_epj_cu[nj_tot] = id_spj[iw][j]+n_epj_tot;
            }
        }
        for(int i=ni_tot; i<ni_tot_reg; i++){
            grav_epi_cu[i].id_walk = n_walk;
        }
        grav_epi_cu.htod(ni_tot_reg);
        id_epj_cu.htod(nj_tot);
        int nblocks  = ni_tot_reg / N_THREAD_GPU;
        int nthreads = N_THREAD_GPU;
        //const float eps2 = FPGrav::eps * FPGrav::eps;
        const float eps2 = epi[0][0].getEps2();
        GravKernel <<<nblocks, nthreads>>> (ij_disp_cu, grav_epi_cu, grav_epj_cu, id_epj_cu, grav_res_cu, eps2);
        return 0;
    }
}


//////////////
// retrieve //

template PS::S32 RetrieveKernel(const PS::S32 tag,
                                const PS::S32 n_walk,
                                const PS::S32 * ni,
                                RESULT::Dens ** force);

template PS::S32 RetrieveKernel(const PS::S32 tag,
                                const PS::S32 n_walk,
                                const PS::S32 * ni,
                                RESULT::Grav ** force);

template PS::S32 RetrieveKernel(const PS::S32 tag,
                                const PS::S32 n_walk,
                                const PS::S32 * ni,
                                RESULT::Drvt ** force);

template PS::S32 RetrieveKernel(const PS::S32 tag,
                                const PS::S32 n_walk,
                                const PS::S32 * ni,
                                RESULT::Hydro ** force);

template<class Tforce>
inline void CopyForceFromDevice(const PS::S32 ni_tot){assert(0);}

template<class Tforce>
inline void CopyForce(Tforce *dst[],
                      const int i_src, 
                      const int iw, 
                      const int i_dst){assert(0);}


// Dens //
template<>
inline void CopyForceFromDevice<RESULT::Dens>(const PS::S32 ni_tot){
    dens_res_cu.dtoh(ni_tot);
}

template<>
inline void CopyForce<RESULT::Dens>(RESULT::Dens *dst[],
                                    const int i_src, 
                                    const int iw, 
                                    const int i_dst){
    dst[iw][i_dst].smth = dens_res_cu[i_src].smth;
    dst[iw][i_dst].dens = dens_res_cu[i_src].dens;
}

// Drvt //
template<>
inline void CopyForceFromDevice<RESULT::Drvt>(const PS::S32 ni_tot){
    drvt_res_cu.dtoh(ni_tot);
}

template<>
inline void CopyForce<RESULT::Drvt>(RESULT::Drvt *dst[],
                                    const int i_src, 
                                    const int iw, 
                                    const int i_dst){
    dst[iw][i_dst].div_v = drvt_res_cu[i_src].div_v;
    dst[iw][i_dst].rot_v.x = drvt_res_cu[i_src].rot_v.x;
    dst[iw][i_dst].rot_v.y = drvt_res_cu[i_src].rot_v.y;
    dst[iw][i_dst].rot_v.z = drvt_res_cu[i_src].rot_v.z;
}

// Hydr //
template<>
inline void CopyForceFromDevice<RESULT::Hydro>(const PS::S32 ni_tot){
    hydro_res_cu.dtoh(ni_tot);
}

template<>
inline void CopyForce<RESULT::Hydro>(RESULT::Hydro *dst[],
                                     const int i_src, 
                                     const int iw, 
                                     const int i_dst){
    dst[iw][i_dst].acc.x = hydro_res_cu[i_src].acc.x;
    dst[iw][i_dst].acc.y = hydro_res_cu[i_src].acc.y;
    dst[iw][i_dst].acc.z = hydro_res_cu[i_src].acc.z;
    dst[iw][i_dst].eng_dot = hydro_res_cu[i_src].eng_dot;
    dst[iw][i_dst].dt = hydro_res_cu[i_src].dt;
}

// Grav //
template<>
inline void CopyForceFromDevice<RESULT::Grav>(const PS::S32 ni_tot){
    grav_res_cu.dtoh(ni_tot);
}

template<>
inline void CopyForce<RESULT::Grav>(RESULT::Grav *dst[],
                                    const int i_src, 
                                    const int iw, 
                                    const int i_dst){
    dst[iw][i_dst].acc.x = grav_res_cu[i_src].accp.x;
    dst[iw][i_dst].acc.y = grav_res_cu[i_src].accp.y;
    dst[iw][i_dst].acc.z = grav_res_cu[i_src].accp.z;
    dst[iw][i_dst].pot   = grav_res_cu[i_src].accp.w;
}

template<class Tforce>
PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       Tforce ** force){
    //std::cerr<<"Before Retrieve"<<std::endl;
    int ni_tot = 0;
    for(int k=0; k<n_walk; k++){
        ni_tot += ni[k];
    }
    CopyForceFromDevice<Tforce>(ni_tot);
    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<ni[iw]; i++, n_cnt++){
            CopyForce<Tforce>(force, n_cnt, iw, i);
        }
    }
    //std::cerr<<"After Retrieve"<<std::endl;
    return 0;
}




