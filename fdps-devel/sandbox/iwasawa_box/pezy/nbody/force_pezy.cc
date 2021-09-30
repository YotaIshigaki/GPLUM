#include "force_multiwalk.hpp"

void InitializeDEVICE(){
    cl_int ret:
    cl_platform_id platform_id = NULL;
    cl_uint ret_num_platforms;
    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    std::cerr<<"platform_id="<<platform_id<<" ret_num_platforms="<<ret_num_platforms<<" ret="<<ret<<std::endl;
    cl_device_id device_id = NULL;
    cl_uint ret_num_devices;
    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
    std::cerr<<"device_id="<<device_id<<" ret_num_devices="<<ret_num_devices<<" ret="<<ret<<std::endl;
    cl_context context = NULL;
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

    cl_command_queue command_queue = NULL;
    command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
}

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
        for(int i=0; i<ni[iw]; i++){
            force[iw][i].acc.x = dev_force[n_cnt].accp.x;
            force[iw][i].acc.y = dev_force[n_cnt].accp.y;
            force[iw][i].acc.z = dev_force[n_cnt].accp.z;
            force[iw][i].pot   = dev_force[n_cnt].accp.w;
            n_cnt++;
        }
    }
    return 0;
}
v
