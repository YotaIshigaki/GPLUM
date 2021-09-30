
//#define DEBUGMODE_DIRECT

#define NO_FORCE_CALCULATION

#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/time.h>
#include<particle_simulator.hpp>

#include "user-defined.hpp"

inline double GetWtime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec) + 1.0e-6*(tv.tv_usec);
}


void GetTimeProfileMax(const PS::TimeProfile & tp, const PS::S32 rank, PS::TimeProfile & tp_max, PS::S32 rank_max[]){
    PS::S32 id = 0;
    PS::Comm::getMaxValue(tp.collect_sample_particle, rank, tp_max.collect_sample_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain, rank, tp_max.decompose_domain, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_particle, rank, tp_max.exchange_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_particle_local_tree, rank, tp_max.set_particle_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_root_cell, rank, tp_max.set_root_cell, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_local_tree, rank, tp_max.make_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_moment_local_tree, rank, tp_max.calc_moment_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_LET_1st, rank, tp_max.make_LET_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_1st, rank, tp_max.exchange_LET_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_LET_2nd, rank, tp_max.make_LET_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_2nd, rank, tp_max.exchange_LET_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_particle_global_tree, rank, tp_max.set_particle_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_global_tree, rank, tp_max.make_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_moment_global_tree, rank, tp_max.calc_moment_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_force, rank, tp_max.calc_force, rank_max[id++]);

    PS::Comm::getMaxValue(tp.make_local_tree_tot, rank, tp_max.make_local_tree_tot, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_global_tree_tot, rank, tp_max.make_global_tree_tot, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_tot, rank, tp_max.exchange_LET_tot, rank_max[id++]);

    PS::Comm::getMaxValue(tp.morton_sort_local_tree, rank, tp_max.morton_sort_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.link_cell_local_tree, rank, tp_max.link_cell_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.morton_sort_global_tree, rank, tp_max.morton_sort_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.link_cell_global_tree, rank, tp_max.link_cell_global_tree, rank_max[id++]);

    PS::Comm::getMaxValue(tp.calc_force__core, rank, tp_max.calc_force__core, rank_max[id++]);
}

#define MULTI_WALK
#ifdef MULTI_WALK
const int N_WALK_LIMIT = 200;
#endif

double WTIME_TOTAL;
double WTIME_TOTAL_OFFSET;
double WTIME_FORCE;
double WTIME_FORCE_OFFSET;
double WTIME_SEND_IP;
double WTIME_SEND_IP_OFFSET;
double WTIME_SEND_JP;
double WTIME_SEND_JP_OFFSET;
double WTIME_RECV;
double WTIME_RECV_OFFSET;
PS::S64 N_IP_SEND;
PS::S64 N_JP_SEND;

const int NI_LIMIT = 8192 * N_WALK_LIMIT;
const int NJ_LIMIT = 131072 * N_WALK_LIMIT;
const int N_THREAD_MAX = 8192;
#include"class_device.hpp"
int j_disp_h[N_WALK_LIMIT+2];
EpiDev epi_h[NI_LIMIT];
EpjDev epj_h[NJ_LIMIT];
ForceDev force_h[NI_LIMIT];


void ClearProfile(){
    WTIME_TOTAL = WTIME_TOTAL_OFFSET = 0.0;
    WTIME_FORCE = WTIME_FORCE_OFFSET = 0.0;
    WTIME_SEND_IP = WTIME_SEND_IP_OFFSET = 0.0;
    WTIME_SEND_JP = WTIME_SEND_JP_OFFSET = 0.0;
    WTIME_RECV = WTIME_RECV_OFFSET = 0.0;
    N_IP_SEND = N_JP_SEND = 0;
}

#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif

#ifdef ENABLE_GPU_CUDA
#define MULTI_WALK
#include"force_multiwalk.hpp"
#endif

#ifdef ENABLE_PEZY
#define MULTI_WALK

#include<pzcl/pzcl_ocl_wrapper.h>
#include<PZSDKHelper.h>


///////////////////////////////////////////////////////////
static bool LoadFile( const char* name, size_t size,  char* pData )
{
	FILE* fp = fopen( name, "rb");
	if( fp == NULL )
	{
		printf("can not open %s\n", name);
		return false;
	}

	if( size == 0 || pData == NULL)
	{
		printf("invalid params %s\n", __FUNCTION__);
		return false;
	}

	size_t size_ret = fread(pData, sizeof(char), size, fp);
	fclose(fp);

	if( size_ret != size )
	{
		printf("can not read requested size\n");
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////
static size_t GetFileSize(const char* name)
{
	FILE* fp = fopen( name, "rb");
	if( fp == NULL )
	{
		printf("can not open %s", name);
		return 0;
	}
	fseek(fp, 0, SEEK_END);

	size_t size = ftell(fp);
	fclose(fp);

	return size;
}

/**
 * Crate program object
 */
cl_program CreateProgram(
	cl_context context,
	std::vector<cl_device_id> &device_id_lists,
	const char* bin_name
	)
{
	cl_program program = NULL;
	char* pBin = NULL;
	cl_int result;

	size_t sizeFile = GetFileSize( bin_name );
	if( sizeFile == 0 )
	{
		goto leaving;
	}

	PZSDK_ALIGNED_ALLOC(pBin, sizeFile, 8 /*8 byte alignment*/);
	if( pBin == NULL)
	{
		printf("out of host memory\n");
		goto leaving;
	}

	if( !LoadFile(bin_name, sizeFile, pBin))
	{
		goto leaving;
	}

	{
		const unsigned char* listBin[1];
		listBin[0] = (unsigned char*)pBin;
		cl_int binary_status = CL_SUCCESS;
		size_t length = sizeFile;
		program = clCreateProgramWithBinary( context, (cl_uint)device_id_lists.size(), &device_id_lists[0], &length, listBin, &binary_status, &result);
        /*
        std::cerr<<"PS::Comm::getRank()= "<<PS::Comm::getRank()<<std::endl;
        std::cerr<<"(cl_uint)device_id_lists.size()="<<(cl_uint)device_id_lists.size()<<std::endl;
        std::cerr<<"device_id_lists[0]="<<device_id_lists[0]<<std::endl;
        std::cerr<<"length= "<<length<<std::endl;
        */
	}
	if(program == NULL)
	{
		printf("clCreateProgramWithBinary failed, %d\n", result);
		goto leaving;
	}
leaving:
	if(pBin)
	{
		PZSDK_ALIGNED_FREE(pBin);
	}
	return program;
}

const int N_DEV_MAX_PER_NODE = 4;



cl_mem j_disp_d = NULL;
cl_mem epi_d = NULL;
cl_mem epj_d = NULL;
cl_mem force_d = NULL;

cl_platform_id platform_id = NULL;

//cl_device_id device_id[N_DEV_MAX*2];
static cl_device_id * device_id = NULL;
//cl_uint ret_num_devices;
cl_context context = NULL;
cl_command_queue command_queue = NULL;
cl_program program = NULL;
cl_kernel kernel = NULL;
cl_event event_force, event_send, event_recv;

std::vector<cl_device_id> device_id_lists;

void InitializeDEVICE(){
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    //device_id = new cl_device_id[N_DEV_MAX_PER_NODE];
    device_id = new cl_device_id[n_proc];
    cl_int ret;
    cl_uint num_platforms;
    ret = clGetPlatformIDs(1, &platform_id, &num_platforms);
    //std::cerr<<"platform_id= "<<platform_id<<" num_platforms= "<<num_platforms<<" ret= "<<ret<<std::endl;
    cl_uint num_devices;
    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, n_proc, device_id, &num_devices);
    //ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, N_DEV_MAX_PER_NODE, device_id, &num_devices);
    //ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 0, device_id, &num_devices);
    /*
    for(int r=0; r<n_proc; r++){
        if(my_rank == r){
            std::cerr<<"my_rank= "<<my_rank<<" platform_id="<<platform_id<<" num_devices= "<<num_devices<<std::endl;
            //for(int i=0; i<N_DEV_MAX_PER_NODE; i++){
            for(int i=0; i<n_proc; i++){
                std::cerr<<"device_id[i]= "<<device_id[i]<<std::endl;
            }
        }
        PS::Comm::barrier();
        usleep(100000);
    }
    */
    int rank_tmp = my_rank % num_devices;
    //context = clCreateContext(NULL, 1, &device_id[my_rank], NULL, NULL, &ret);
    //command_queue = clCreateCommandQueue(context, device_id[my_rank], 0, &ret);
    context = clCreateContext(NULL, 1, &device_id[rank_tmp], NULL, NULL, &ret);
    command_queue = clCreateCommandQueue(context, device_id[rank_tmp], 0, &ret);

    j_disp_d = clCreateBuffer(context, CL_MEM_READ_WRITE, (N_WALK_LIMIT+2)*sizeof(int),   NULL, &ret);
    epi_d   = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT*sizeof(EpiDev),   NULL, &ret);
    epj_d   = clCreateBuffer(context, CL_MEM_READ_WRITE, NJ_LIMIT*sizeof(EpjDev),   NULL, &ret);
    force_d = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT*sizeof(ForceDev), NULL, &ret);


    //device_id_lists.push_back(device_id[my_rank]);
    device_id_lists.push_back(device_id[rank_tmp]);
    program = CreateProgram(context, device_id_lists, "./kernel.sc32/kernel.pz");

    if(program == NULL){
        std::cerr<<"can't create program"<<std::endl;
        exit(1);
    }
    kernel = clCreateKernel(program, "ForceKernel", &ret);
    if(kernel == NULL){
        std::cerr<<"can't create kernel"<<std::endl;
        exit(1);
    }
    ret = clSetKernelArg(kernel, 6, sizeof(int),    (void*)&N_THREAD_MAX);

}

PS::F64 WTIME_SET_OFFSET = 0.0;
PS::F64 WTIME_SET = 0.0;

PS::S32 GetQuantizedValue(const PS::S32 val, const PS::S32 div){
    
#ifdef DEBUGMODE_DIRECT
    return 128;
#else
    return div*((val-1)/div+1);
#endif
}
#endif // ENABLE_PEZY

const int NI_PIPE = 4;

PS::S32 DispatchKernelWithSP(
                             const PS::S32          tag,
                             const PS::S32          n_walk,
                             const EPI          *epi[],
                             const PS::S32          n_epi[],
                             const EPJ          *epj[],
                             const PS::S32          n_epj[],
                             const SPJ         *spj[],
                             const PS::S32          n_spj[]){
/*
PS::S32 DispatchKernelWithSP(
                             const PS::S32          tag,
                             const PS::S32          n_walk,
                             const EPI          *epi[],
                             const PS::S32          n_epi[],
                             const EPJ          *epj[],
                             const PS::S32          n_epj[],
                             const PS::SPJMonopole *spj[],
                             const PS::S32          n_spj[]){
*/
/*
PS::S32 DispatchKernelWithSP(
                             const PS::S32          tag,
                             const PS::S32          n_walk,
                             const FPGrav          *epi[],
                             const PS::S32          n_epi[],
                             const FPGrav          *epj[],
                             const PS::S32          n_epj[],
                             const PS::SPJMonopole *spj[],
                             const PS::S32          n_spj[]){
*/
/*
PS::S32 DispatchKernelWithSP(
                             const PS::S32          tag,
                             const PS::S32          n_walk,
                             const FPGrav          *epi[],
                             const PS::S32          n_epi[],
                             const EPJ          *epj[],
                             const PS::S32          n_epj[],
                             const PS::SPJMonopole *spj[],
                             const PS::S32          n_spj[]){
*/
#ifndef NO_FORCE_CALCULATION
    assert(n_walk <= N_WALK_LIMIT);
    ::WTIME_TOTAL_OFFSET = GetWtime();
    cl_int ret;
    const float eps2 = FPGrav::eps * FPGrav::eps;
    PS::S32 ni_tot = 0;
    j_disp_h[0] = 0;
    for(int k=0; k<n_walk; k++){
        //ni_tot += n_epi[k];
        ni_tot += GetQuantizedValue(n_epi[k], NI_PIPE);
        assert(GetQuantizedValue(n_epi[k], NI_PIPE) % NI_PIPE == 0);
        j_disp_h[k+1] = j_disp_h[k] + (n_epj[k] + n_spj[k]);
    }
    j_disp_h[n_walk+1] = j_disp_h[n_walk];
    assert(ni_tot < NI_LIMIT);
    assert(j_disp_h[n_walk] < NJ_LIMIT);
    ret = clEnqueueWriteBuffer(command_queue, j_disp_d, CL_TRUE, 0, (n_walk+2)*sizeof(int), j_disp_h, 0, NULL, NULL);

    ni_tot = 0;
    int nj_tot = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<n_epi[iw]; i++){
            epi_h[ni_tot].px = epi[iw][i].pos.x;
            epi_h[ni_tot].py = epi[iw][i].pos.y;
            epi_h[ni_tot].pz = epi[iw][i].pos.z;
            epi_h[ni_tot].id_walk = iw;
            ni_tot++;
        }
        for(int i=n_epi[iw]; i<GetQuantizedValue(n_epi[iw], NI_PIPE); i++){
            epi_h[ni_tot].px = 0.0;
            epi_h[ni_tot].py = 0.0;
            epi_h[ni_tot].pz = 0.0;
            epi_h[ni_tot].id_walk = iw;
            ni_tot++;
        }
        for(int j=0; j<n_epj[iw]; j++){
            epj_h[nj_tot].px  = epj[iw][j].pos.x;
            epj_h[nj_tot].py  = epj[iw][j].pos.y;
            epj_h[nj_tot].pz  = epj[iw][j].pos.z;
            epj_h[nj_tot].mass  = epj[iw][j].mass;
            nj_tot++;
        }
        for(int j=0; j<n_spj[iw]; j++){
            epj_h[nj_tot].px  = spj[iw][j].pos.x;
            epj_h[nj_tot].py  = spj[iw][j].pos.y;
            epj_h[nj_tot].pz  = spj[iw][j].pos.z;
            epj_h[nj_tot].mass  = spj[iw][j].getCharge();
            nj_tot++;
        }
    }
    N_IP_SEND += ni_tot;
    N_JP_SEND += nj_tot;
    ::WTIME_SEND_IP_OFFSET = GetWtime();
    ret = clEnqueueWriteBuffer(command_queue, epi_d,   CL_TRUE, 0, (ni_tot)*sizeof(EpiDev),  epi_h, 0, NULL, NULL);
    ::WTIME_SEND_IP += GetWtime() - WTIME_SEND_IP_OFFSET;
    ::WTIME_SEND_JP_OFFSET = GetWtime();
    ret = clEnqueueWriteBuffer(command_queue, epj_d,   CL_TRUE, 0, (nj_tot)*sizeof(EpjDev),  epj_h, 0, NULL, NULL);
    ::WTIME_SEND_JP += GetWtime() - WTIME_SEND_JP_OFFSET;

    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&j_disp_d);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&epi_d);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&epj_d);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*)&force_d);
    ret = clSetKernelArg(kernel, 4, sizeof(float),  (void*)&eps2);
    ret = clSetKernelArg(kernel, 5, sizeof(int),    (void*)&ni_tot);
    assert(ni_tot % NI_PIPE == 0 );
    size_t work_size = N_THREAD_MAX;
    ::WTIME_FORCE_OFFSET = GetWtime();
#if 1
    ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &work_size, NULL, 0, NULL, &event_force);
    //ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
#else
    for(int ith=0; ith<8; ith++){
        for(int ipr=0; ipr<N_THREAD_MAX/8; ipr++){
            cpu_ForceKernel(j_disp_h, epi_h, epj_h, force_h, eps2, ni_tot, N_THREAD_MAX, ith, ipr);
        }
    }
#endif
#endif //NO_FORCE_CALCULATION
    return 0;
}


PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       Force    *force[]){
    /*
PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       FPGrav    *force[]){
    */
#ifndef NO_FORCE_CALCULATION
    clWaitForEvents(1, &event_force);
    ::WTIME_FORCE += GetWtime() - WTIME_FORCE_OFFSET;
    cl_int ret;
    int ni_tot = 0;
    for(int k=0; k<n_walk; k++){
        //ni_tot += ni[k];
        ni_tot += GetQuantizedValue(ni[k], NI_PIPE);
    }
    assert(ni_tot % NI_PIPE == 0 );
    ::WTIME_RECV_OFFSET = GetWtime();
    ret = clEnqueueReadBuffer(command_queue, force_d, CL_TRUE, 0, ni_tot*sizeof(ForceDev), force_h, 0, NULL, &event_recv);
    clWaitForEvents(1, &event_recv);
    ::WTIME_RECV += GetWtime() - WTIME_RECV_OFFSET;
    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<ni[iw]; i++){
            force[iw][i].acc.x = force_h[n_cnt].ax;
            force[iw][i].acc.y = force_h[n_cnt].ay;
            force[iw][i].acc.z = force_h[n_cnt].az;
            force[iw][i].pot   = force_h[n_cnt].pot;
            n_cnt++;
        }
        n_cnt += GetQuantizedValue(ni[iw], NI_PIPE) - ni[iw];
    }
    ::WTIME_TOTAL += GetWtime() - WTIME_TOTAL_OFFSET;
    assert(n_cnt == ni_tot );
#endif //NO_FORCE_CALCULATION
    return 0;
}






void MakePlummerModel(const double mass_glb,
                      const long long int n_glb,
                      const long long int n_loc,
                      double *& mass,
                      PS::F64vec *& pos,
                      PS::F64vec *& vel,
                      const double eng = -0.25,
                      const int seed = 0){

    assert(eng < 0.0);
    static const double PI = atan(1.0) * 4.0;
    const double r_cutoff = 22.8 / (-3.0 * PI * mass_glb * mass_glb / (64.0 * -0.25)); // 22.8 is cutoff in Nbody units
    //const double r_cutoff = 22.8 * 0.25;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];

    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        double r_tmp = 9999.9;
        while(r_tmp > r_cutoff){ 
            double m_tmp = mt.genrand_res53();
            r_tmp = 1.0 / sqrt( pow(m_tmp, (-2.0/3.0)) - 1.0);
        }
        double phi = 2.0 * PI * mt.genrand_res53();
        double cth = 2.0 * (mt.genrand_real2() - 0.5);
        double sth = sqrt(1.0 - cth*cth);
        pos[i][0] = r_tmp * sth * cos(phi);
        pos[i][1] = r_tmp * sth * sin(phi);
        pos[i][2] = r_tmp * cth;
        while(1){
            const double v_max = 0.1;
            const double v_try = mt.genrand_res53();
            const double v_crit = v_max * mt.genrand_res53();
            if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                const double ve = sqrt(2.0) * pow( (r_tmp*r_tmp + 1.0), -0.25 );
                phi = 2.0 * PI * mt.genrand_res53();
                cth = 2.0 * (mt.genrand_res53() - 0.5);
                sth = sqrt(1.0 - cth*cth);
                vel[i][0] = ve * v_try * sth * cos(phi);
                vel[i][1] = ve * v_try * sth * sin(phi);
                vel[i][2] = ve * v_try * cth;
                break;
            }
        }
    }

    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    double  cm_mass = 0.0;
    for(long long int i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(long long int i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(long long int i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }

    double r_max_sq = -1.0;
    for(int i=0; i<n_loc; i++){
        if(r_max_sq < pos[i] * pos[i]){
            r_max_sq = pos[i] * pos[i];
        }
    }
    //std::cout<<"r_max= "<<sqrt(r_max_sq)<<std::endl;
}


template<class Tpsys>
void setParticlesPlummer(Tpsys & psys,
                         const PS::S64 n_glb,
                         PS::S32 & n_loc){
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
#if 1
    if(my_rank == 0){
        n_loc = n_glb;
    }
    else{
        n_loc = 0;
    }
#else
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
#endif
    psys.setNumberOfParticleLocal(n_loc);
    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(long long int i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

void makeColdUniformSphere(const PS::F64 mass_glb,
                           const PS::S64 n_glb,
                           const PS::S64 n_loc,
                           PS::F64 *& mass,
                           PS::F64vec *& pos,
                           PS::F64vec *& vel,
                           const PS::F64 eng = -0.25,
                           const PS::S32 seed = 0) {
    
    assert(eng < 0.0);
    {
        PS::MTTS mt;
        mt.init_genrand(0);
        for(PS::S32 i = 0; i < n_loc; i++){
            mass[i] = mass_glb / n_glb;
            const PS::F64 radius = 3.0;
            do {
                pos[i][0] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][1] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][2] = (2. * mt.genrand_res53() - 1.) * radius;
            }while(pos[i] * pos[i] >= radius * radius);
            vel[i][0] = 0.0;
            vel[i][1] = 0.0;
            vel[i][2] = 0.0;
        }
    }

    PS::F64vec cm_pos  = 0.0;
    PS::F64vec cm_vel  = 0.0;
    PS::F64    cm_mass = 0.0;
    for(PS::S32 i = 0; i < n_loc; i++){
        cm_pos  += mass[i] * pos[i];
        cm_vel  += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(PS::S32 i = 0; i < n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }
}

template<class Tpsys>
void setParticlesColdUniformSphere(Tpsys & psys,
                                   const PS::S32 n_glb,
                                   PS::S32 & n_loc) {

    n_loc = n_glb; 
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64    * mass = new PS::F64[n_loc];
    PS::F64vec * pos  = new PS::F64vec[n_loc];
    PS::F64vec * vel  = new PS::F64vec[n_loc];
    const PS::F64 m_tot = 1.0;
    const PS::F64 eng   = -0.25;
    makeColdUniformSphere(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    for(PS::S32 i = 0; i < n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos  = pos[i];
        psys[i].vel  = vel[i];
        psys[i].id   = i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

template<class Tpsys>
void kick(Tpsys & system,
          const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i = 0; i < n; i++) {
        system[i].pos  += system[i].vel * dt;
    }
}

template<class Tpsys>
void calcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
    if(clear){
        etot = ekin = epot = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * (system[i].pot + system[i].mass / FPGrav::eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    //std::cerr<<"ekin_loc="<<ekin_loc<<std::endl;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
#else
    etot = etot_loc;
    epot = epot_loc;
    ekin = ekin_loc;
#endif
}

void printHelp() {
    std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
    std::cerr<<"s: time_step (default: 1.0 / 128.0)"<<std::endl;
    std::cerr<<"d: dt_diag (default: 1.0 / 8.0)"<<std::endl;
    std::cerr<<"D: dt_snap (default: 1.0)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

void makeOutputDirectory(char * dir_name) {
    struct stat st;
    if(stat(dir_name, &st) != 0) {
        PS::S32 ret_loc = 0;
        PS::S32 ret     = 0;
        if(PS::Comm::getRank() == 0)
            ret_loc = mkdir(dir_name, 0777);
        PS::Comm::broadcast(&ret_loc, ret);
        if(ret == 0) {
            if(PS::Comm::getRank() == 0)
                fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
        } else {
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
            PS::Abort();
        }
    }
}

PS::F64 FPGrav::eps = 1.0/32.0;
//PS::F32 FPGrav::eps = 1.0/32.0;

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);

    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F32 dt = 1.0 / 128.0;
    PS::F32 dt_diag = dt;
    //PS::F32 time_end = dt;
    //PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 time_end = 10.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_tot = 1024;
    PS::S32 c;
    sprintf(dir_name,"./result");
    opterr = 0;
    while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 't':
            theta = atof(optarg);
            std::cerr << "theta= " << theta << std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr << "time_end= " << time_end << std::endl;
            break;
        case 's':
            dt = atof(optarg);
            std::cerr << "time_step= " << dt << std::endl;
            break;
        case 'd':
            dt_diag = atof(optarg);
            std::cerr << "dt_diag= " << dt_diag << std::endl;
            break;
        case 'D':
            dt_snap = atof(optarg);
            std::cerr << "dt_snap= " << dt_snap << std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr << "n_leaf_limit= " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr << "n_group_limit= " << n_group_limit << std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            std::cerr << "n_tot= " << n_tot << std::endl;
            break;
        case 'h':
            if(PS::Comm::getRank() == 0) {
                printHelp();
            }
            PS::Finalize();
            return 0;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
                printHelp();
            }
            PS::Abort();
        }
    }
#ifdef DEBUGMODE_DIRECT
    n_tot = 128;
    n_group_limit = n_tot*2;
    theta = 0.0;
#endif


    makeOutputDirectory(dir_name);

    std::ofstream fout_eng;
    char sout_de[1024];
    sprintf(sout_de, "%s/t-de.dat", dir_name);
    std::cerr << sout_de << std::endl;
    fout_eng.open(sout_de);

    if(PS::Comm::getRank() == 0) {
        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(1000);
    PS::S32 n_loc    = 0;
    PS::F32 time_sys = 0.0;
#if 1
    setParticlesPlummer(system_grav, n_tot, n_loc);
#else
    if(PS::Comm::getRank() == 0) {
        setParticlesColdUniformSphere(system_grav, n_tot, n_loc);
        
    } else {
        system_grav.setNumberOfParticleLocal(n_loc);
    }
#endif
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);

    dinfo.collectSampleParticle(system_grav);

    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);

    n_loc = system_grav.getNumberOfParticleLocal();

    PS::F64 ekin_tmp = 0.0;
    for(int i=0; i<n_loc; i++){
        ekin_tmp += 0.5*system_grav[i].mass*system_grav[i].vel*system_grav[i].vel;
    }
    //std::cerr<<"n_loc="<<n_loc<<" ekin_tmp="<<ekin_tmp<<std::endl;

#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FPGrav::eps);
#endif


#ifdef ENABLE_PEZY
    InitializeDEVICE();
#endif

    PS::TreeForForce<PS::SEARCH_MODE_LONG, Force, EPI, EPJ, Moment, Moment, SPJ> tree_grav;
    //PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;
    //PS::TreeForForceLong<Force, FPGrav, FPGrav>::Monopole tree_grav;
    //PS::TreeForForceLong<Force, FPGrav, EPJ>::Monopole tree_grav;
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    tree_grav.clearCounterAll();
    ClearProfile();
#ifdef MULTI_WALK
    const PS::S32 n_walk_limit = N_WALK_LIMIT;
    const PS::S32 tag_max = 1;
    //std::cerr<<"n_walk_limit="<<n_walk_limit<<std::endl;
    tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                RetrieveKernel,
                                                tag_max,
                                                system_grav,
                                                dinfo,
                                                n_walk_limit);
#else
    tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       system_grav,
                                       dinfo);
#endif



    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    //std::cerr<<"system_grav[0].pot="<<system_grav[0].pot<<"system_grav[0].acc.x= "<<system_grav[0].acc.x<<std::endl;
    calcEnergy(system_grav, Etot0, Ekin0, Epot0);
    PS::S64 n_op = 29;
    PS::S64 n_interaction = tree_grav.getNumberOfInteractionEPEPGlobal() + tree_grav.getNumberOfInteractionEPSPGlobal();
    if(PS::Comm::getRank() == 0){
        std::cerr<<"Ekin0= "<<Ekin0<<" Epot0= "<<Epot0<<" Etot0= "<<Etot0<<std::endl;
        /*
        std::cerr<<"WTIME_FORCE= "<<WTIME_FORCE<<std::endl;
        std::cerr<<"n_interaction="<<n_interaction<<std::endl;
        std::cerr<<"speed= "<< (double)( n_interaction * n_op) / WTIME_FORCE / 1e12<<" [Tflops]"<<std::endl;
        */
        /*
        for(PS::S32 i=0; i<n_loc; i++){
            std::cout<<"system_grav[i].mass= "<<system_grav[i].mass
                     <<" system_grav[i].pos= "<<system_grav[i].pos
                     <<" system_grav[i].acc= "<<system_grav[i].acc
                     <<" system_grav[i].pot= "<<system_grav[i].pot
                     <<std::endl;
        }
        */
    }

#if 1
    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S32 id_snap = 0;
    PS::F64 wtime_offset = 0.0;
    PS::F64 wtime_tot = 0.0;
    PS::F64 weight = 0.0;
    PS::S32 n_loop_measure = 0;
    PS::F64 wtime_kick = 0.0;
    PS::F64 wtime_kick_offset = 0.0;
    PS::F64 wtime_drift = 0.0;
    PS::F64 wtime_drift_offset = 0.0;
    PS::F64 wtime_ex_ptcl = 0.0;
    PS::F64 wtime_ex_ptcl_offset = 0.0;
    PS::F64 wtime_tree = 0.0;
    PS::F64 wtime_tree_offset = 0.0;
    PS::N_CALL_TRAVERSE = 0;
    //while(time_sys < time_end){
    while(n_loop <= 16){
        if( (time_sys >= time_snap) || ( (time_sys + dt) - time_snap ) > (time_snap - time_sys) ){
            char filename[256];
            sprintf(filename, "%s/%04d.dat", dir_name, id_snap++);
            FileHeader header;
            header.time   = time_sys;
            header.n_body = system_grav.getNumberOfParticleGlobal();
            system_grav.writeParticleAscii(filename, header);
            time_snap += dt_snap;
        }

        calcEnergy(system_grav, Etot1, Ekin1, Epot1);

        PS::S64 n_interaction = tree_grav.getNumberOfInteractionEPEPGlobal() + tree_grav.getNumberOfInteractionEPSPGlobal();
        PS::S64 n_walk = tree_grav.getNumberOfWalkGlobal();
        PS::S64 n_glb = system_grav.getNumberOfParticleGlobal();
        PS::TimeProfile tp_grav = tree_grav.getTimeProfile();
        weight = tp_grav.calc_force;
        //if( (time_sys >= time_diag) || ( (time_sys + dt) - time_diag ) > (time_diag - time_sys) ){
        if(n_loop % 8 == 0 && n_loop > 0){
            const PS::S32 n_profile = 30;
            PS::S32 rank_dinfo_max[n_profile], rank_system_max[n_profile], rank_dens_max[n_profile], 
                rank_hydr_max[n_profile], rank_grav_max[n_profile];
            PS::TimeProfile tp_dinfo_max, tp_system_max, tp_dens_max, tp_hydr_max, tp_grav_max;
            GetTimeProfileMax(tp_grav, PS::Comm::getRank(), tp_grav_max, rank_grav_max);
            if(PS::Comm::getRank() == 0){
                fout_eng << time_sys << "   " << (Etot1 - Etot0) / Etot0 << std::endl;
                fprintf(stderr, "time: %10.7f energy error: %+e\n",
                        time_sys, (Etot1 - Etot0) / Etot0);
                time_diag += dt_diag;
                //if(time_sys+dt >= time_end){
                if(n_loop == 16){
                    std::cout<<"n_proc= "<<PS::Comm::getNumberOfProc()
                             <<" n_glb= "<<n_glb
                             <<" wtime_tot= "<<wtime_tot/n_loop_measure
                             <<" (n_interaction*n_op)= "<<(double)(n_interaction*n_op)/1e9/n_loop_measure
                             <<" [Gops] "
                             <<" <ni>= "<<(double)n_glb/(n_walk/n_loop_measure) // 10, 11
                             <<" <nj>= "<<(double)n_interaction/n_glb/n_loop_measure // 12, 13
                             <<" speed= "<<(PS::F64)(n_interaction*n_op)/wtime_tot/1e12
                             <<" [Tflops] "
                             <<" n_walk/n_proc= "<<n_walk/PS::Comm::getNumberOfProc()/n_loop_measure // 17, 18
                             <<" n_group_limit= "<<n_group_limit   // 19, 20
                             <<" WTIME_TOTAL= "<<WTIME_TOTAL/n_loop_measure
                             <<" WTIME_FORCE= "<<WTIME_FORCE/n_loop_measure
                             <<" WTIME_SEDN_IP= "<<WTIME_SEND_IP/n_loop_measure
                             <<" WTIME_SEDN_JP= "<<WTIME_SEND_JP/n_loop_measure
                             <<" WTIME_RECV= "<<WTIME_RECV/n_loop_measure // 29, 30
                             <<" N_IP_SEND= "<<N_IP_SEND/n_loop_measure
                             <<" N_JP_SEND= "<<N_JP_SEND/n_loop_measure
                             <<" speed(pezy_total)= "<<(PS::F64)(n_interaction*n_op)/WTIME_TOTAL/1e12
                             <<" [Tflops] "
                             <<" send_ip_to_pezy= "<<(double)(N_IP_SEND*sizeof(EpiDev)) / WTIME_SEND_IP / 1e9 // 38, 39
                             <<" [GB/sec] "
                             <<" send_jp_to_pezy= "<<(double)(N_JP_SEND*sizeof(EpjDev)) / WTIME_SEND_JP / 1e9 // 41, 42
                             <<" [GB/sec] "
                             <<" recv_from_pezy= "<<(double)(N_IP_SEND*sizeof(ForceDev)) / WTIME_RECV / 1e9   //44, 45
                             <<" [GB/sec] "
                             <<" calc_force= "<<tp_grav.calc_force/n_loop_measure<<" max= "<<tp_grav_max.calc_force/n_loop_measure // 47-50
                             <<" make_local_tree_tot= "<<tp_grav.make_local_tree_tot/n_loop_measure<<" max= "<<tp_grav_max.make_local_tree_tot/n_loop_measure // 51-54
                             <<" make_global_tree_tot= "<<tp_grav.make_global_tree_tot/n_loop_measure<<" max= "<<tp_grav_max.make_global_tree_tot/n_loop_measure // 55-58
                             <<" exchange_LET_tot= "<<tp_grav.exchange_LET_tot/n_loop_measure<<" max= "<<tp_grav_max.exchange_LET_tot/n_loop_measure             // 59-62
                             <<" morton_sort_local_tree= "<<tp_grav.morton_sort_local_tree/n_loop_measure<<" max= "<<tp_grav_max.morton_sort_local_tree/n_loop_measure // 63-66
                             <<" link_cell_local_tree= "<<tp_grav.link_cell_local_tree/n_loop_measure<<" max= "<<tp_grav_max.link_cell_local_tree/n_loop_measure    // 67-70
                             <<" morton_sort_global_tree= "<<tp_grav.morton_sort_global_tree/n_loop_measure<<" max= "<<tp_grav_max.morton_sort_global_tree/n_loop_measure // 71-74
                             <<" link_cell_global_tree= "<<tp_grav.link_cell_global_tree/n_loop_measure<<" max= "<<tp_grav_max.link_cell_global_tree/n_loop_measure  // 75-78
                             <<" wtime_kick= "<<wtime_kick   // 79,80
                             <<" wtime_drift= "<<wtime_drift   // 81,82
                             <<" wtime_ex_ptcl= "<<wtime_ex_ptcl // 83,84
                             <<" wtime_tree= "<<wtime_tree // 85,86
                             <<" calc_force__core= "<<tp_grav.calc_force__core/n_loop_measure<<" max= "<<tp_grav.calc_force__core/n_loop_measure // 87-90
                             <<" calc_force__core__walk_tree= "<<tp_grav.calc_force__core__walk_tree/n_loop_measure<<" max= "<<tp_grav.calc_force__core__walk_tree/n_loop_measure // 91-94
                             <<" n_leaf_limit= "<<n_leaf_limit // 95-96
                             <<" N_CALL_TRAVERSE= "<<(PS::F64)(PS::N_CALL_TRAVERSE)/n_loop_measure
                             <<std::endl;
                }
                std::cerr<<"n_proc= "<<PS::Comm::getNumberOfProc()
                             <<" n_glb= "<<n_glb
                             <<" wtime_tot= "<<wtime_tot/n_loop_measure
                             <<" (n_interaction*n_op)= "<<(double)(n_interaction*n_op)/1e9/n_loop_measure
                             <<" [Gops] "
                             <<" <ni>= "<<(double)n_glb/(n_walk/n_loop_measure)
                             <<" <nj>= "<<(double)n_interaction/n_glb/n_loop_measure
                             <<" speed= "<<(PS::F64)(n_interaction*n_op)/wtime_tot/1e12
                             <<" [Tflops] "
                             <<" n_walk/n_proc= "<<n_walk/PS::Comm::getNumberOfProc()/n_loop_measure
                             <<" n_group_limit= "<<n_group_limit
                             <<" WTIME_TOTAL= "<<WTIME_TOTAL/n_loop_measure
                             <<" WTIME_FORCE= "<<WTIME_FORCE/n_loop_measure
                             <<" WTIME_SEDN_IP= "<<WTIME_SEND_IP/n_loop_measure
                             <<" WTIME_SEDN_JP= "<<WTIME_SEND_JP/n_loop_measure
                             <<" WTIME_RECV= "<<WTIME_RECV/n_loop_measure
                             <<" N_IP_SEND= "<<N_IP_SEND/n_loop_measure
                             <<" N_JP_SEND= "<<N_JP_SEND/n_loop_measure
                             <<" speed(pezy_total)= "<<(PS::F64)(n_interaction*n_op)/WTIME_TOTAL/1e12
                             <<" [Tflops] "
                             <<" send_ip_to_pezy= "<<(double)(N_IP_SEND*sizeof(EpiDev)) / WTIME_SEND_IP / 1e9
                             <<" [GB/sec] "
                             <<" send_jp_to_pezy= "<<(double)(N_JP_SEND*sizeof(EpjDev)) / WTIME_SEND_JP / 1e9
                             <<" [GB/sec] "
                             <<" recv_from_pezy= "<<(double)(N_IP_SEND*sizeof(ForceDev)) / WTIME_RECV / 1e9
                             <<" [GB/sec] "
                             <<" calc_force= "<<tp_grav.calc_force/n_loop_measure<<" max= "<<tp_grav_max.calc_force/n_loop_measure
                             <<" make_local_tree_tot= "<<tp_grav.make_local_tree_tot/n_loop_measure<<" max= "<<tp_grav_max.make_local_tree_tot/n_loop_measure
                             <<" make_global_tree_tot= "<<tp_grav.make_global_tree_tot/n_loop_measure<<" max= "<<tp_grav_max.make_global_tree_tot/n_loop_measure
                             <<" exchange_LET_tot= "<<tp_grav.exchange_LET_tot/n_loop_measure<<" max= "<<tp_grav_max.exchange_LET_tot/n_loop_measure
                             <<" morton_sort_local_tree= "<<tp_grav.morton_sort_local_tree/n_loop_measure<<" max= "<<tp_grav_max.morton_sort_local_tree/n_loop_measure
                             <<" link_cell_local_tree= "<<tp_grav.link_cell_local_tree/n_loop_measure<<" max= "<<tp_grav_max.link_cell_local_tree/n_loop_measure
                             <<" morton_sort_global_tree= "<<tp_grav.morton_sort_global_tree/n_loop_measure<<" max= "<<tp_grav_max.morton_sort_global_tree/n_loop_measure
                             <<" link_cell_global_tree= "<<tp_grav.link_cell_global_tree/n_loop_measure<<" max= "<<tp_grav_max.link_cell_global_tree/n_loop_measure
                             <<" calc_force__core= "<<tp_grav.calc_force__core/n_loop_measure<<" max= "<<tp_grav.calc_force__core/n_loop_measure
                             <<" calc_force__core__walk_tree= "<<tp_grav.calc_force__core__walk_tree/n_loop_measure<<" max= "<<tp_grav.calc_force__core__walk_tree/n_loop_measure
                             <<" n_leaf_limit= "<<n_leaf_limit
                             <<" N_CALL_TRAVERSE= "<<(PS::F64)(PS::N_CALL_TRAVERSE)/n_loop_measure
                             <<std::endl;
                std::cerr<<std::endl;
            } 
            tree_grav.clearCounterAll();
            ClearProfile();
            wtime_tot = 0.0;
            wtime_kick = 0.0;
            wtime_drift = 0.0;
            wtime_ex_ptcl = 0.0;
            wtime_tree = 0.0;
            n_loop_measure = 0;
            PS::N_CALL_TRAVERSE = 0;
        }

        PS::Comm::barrier();
        wtime_offset = PS::GetWtime();
        kick(system_grav, dt * 0.5);
        PS::Comm::barrier();
        wtime_kick += PS::GetWtime() - wtime_offset;

        time_sys += dt;
        wtime_drift_offset = PS::GetWtime();
        drift(system_grav, dt);
        PS::Comm::barrier();
        wtime_drift += PS::GetWtime() - wtime_drift_offset;
        
        if( n_loop < 2){
            dinfo.decomposeDomainAll(system_grav);
        }
        else if(n_loop < 8){
            dinfo.collectSampleParticle(system_grav, true, weight);
            dinfo.decomposeDomainMultiStep();
        }

        PS::Comm::barrier();
        wtime_ex_ptcl_offset = PS::GetWtime();
        system_grav.exchangeParticle(dinfo);
        PS::Comm::barrier();
        wtime_ex_ptcl += PS::GetWtime() - wtime_ex_ptcl_offset;

        wtime_tree_offset = PS::GetWtime();
#ifdef MULTI_WALK
        tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                    RetrieveKernel,
                                                    tag_max,
                                                    system_grav,
                                                    dinfo,
                                                    n_walk_limit,
                                                    true);
#else
        tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                           CalcGravity<PS::SPJMonopole>,
                                           system_grav,
                                           dinfo);
#endif
        PS::Comm::barrier();
        wtime_tree += PS::GetWtime() - wtime_tree_offset;
        wtime_tot += PS::GetWtime() - wtime_offset;
        kick(system_grav, dt * 0.5);
        n_loop++;
        n_loop_measure++;
    }
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

#endif
    PS::Finalize();
    return 0;
}
