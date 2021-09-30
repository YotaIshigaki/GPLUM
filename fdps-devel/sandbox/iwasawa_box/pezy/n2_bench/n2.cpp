#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#include<chrono>
#include<sys/time.h>
#include "user-defined.hpp"

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
#include"class_device.hpp"

inline double GetWtime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec) + 1.0e-6*(tv.tv_usec);
}

double wtime_send  = 0.0;
double wtime_recv  = 0.0;
double wtime_force = 0.0;

struct Particle{
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
};

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

const int N_DEV_MAX = 4;
const int NI_LIMIT = 1000000;
const int NJ_LIMIT = 1000000;
const int N_THREAD_MAX = 8192;

EpiDev epi_h[NI_LIMIT];
EpjDev epj_h[NJ_LIMIT];
ForceDev force_h[NI_LIMIT];

//cl_mem j_disp_d = NULL;
cl_mem epi_d = NULL;
cl_mem epj_d = NULL;
cl_mem force_d = NULL;

cl_platform_id platform_id = NULL;
cl_uint ret_num_platforms;
cl_device_id device_id[N_DEV_MAX];
cl_uint ret_num_devices;
cl_context context = NULL;
cl_command_queue command_queue = NULL;
cl_program program = NULL;
cl_kernel kernel = NULL;
cl_event event_force;

std::vector<cl_device_id> device_id_lists;

void InitializeDEVICE(){
    const PS::S32 my_rank = PS::Comm::getRank();
    std::cerr<<"my_rank="<<my_rank<<std::endl;
    cl_int ret;
    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    std::cerr<<"platform_id="<<platform_id<<" ret_num_platforms="<<ret_num_platforms<<" ret="<<ret<<std::endl;

    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, N_DEV_MAX, device_id, &ret_num_devices);
    std::cerr<<"device_id[my_rank]="<<device_id[my_rank]<<" ret_num_devices="<<ret_num_devices<<" ret="<<ret<<std::endl;

    context = clCreateContext(NULL, 1, &device_id[my_rank], NULL, NULL, &ret);

    command_queue = clCreateCommandQueue(context, device_id[my_rank], 0, &ret);

    epi_d   = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT*sizeof(EpiDev),   NULL, &ret);
    epj_d   = clCreateBuffer(context, CL_MEM_READ_WRITE, NJ_LIMIT*sizeof(EpjDev),   NULL, &ret);
    force_d = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT*sizeof(ForceDev), NULL, &ret);

    device_id_lists.push_back(device_id[my_rank]);
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
    ret = clSetKernelArg(kernel, 5, sizeof(int), (void*)&N_THREAD_MAX);
}

PS::F64 WTIME_SET_OFFSET = 0.0;
PS::F64 WTIME_SET = 0.0;


void calcForce(const PS::F64 mass[], 
               const PS::F64vec pos[], 
               const PS::S64 n_tot, 
               const PS::F32 eps2,
               PS::F64vec acc[],
               PS::F64 pot[]){
    assert(NI_LIMIT > n_tot && NJ_LIMIT > n_tot);
    cl_int ret;
    for(int i=0; i<n_tot; i++){
        epj_h[i].mass  = mass[i];
        epi_h[i].px = epj_h[i].px = pos[i].x;
        epi_h[i].py = epj_h[i].py = pos[i].y;
        epi_h[i].pz = epj_h[i].pz = pos[i].z;
    }
    PS::F64 t0 = GetWtime();
    ret = clEnqueueWriteBuffer(command_queue, epi_d,   CL_TRUE, 0, (n_tot)*sizeof(EpiDev),  epi_h, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, epj_d,   CL_TRUE, 0, (n_tot)*sizeof(EpjDev),  epj_h, 0, NULL, NULL);
    wtime_send += GetWtime() - t0;

    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&epi_d);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&epj_d);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&force_d);
    ret = clSetKernelArg(kernel, 3, sizeof(float),  (void*)&eps2);
    ret = clSetKernelArg(kernel, 4, sizeof(int),    (void*)&n_tot);
    size_t work_size = N_THREAD_MAX;
    t0 = GetWtime();
    //ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &work_size, NULL, 0, NULL, &event_force);
    //clWaitForEvents(1, &event_force);
    ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
    wtime_force += GetWtime() - t0;
    t0 = GetWtime();
    ret = clEnqueueReadBuffer(command_queue, force_d, CL_TRUE, 0, n_tot*sizeof(ForceDev), force_h, 0, NULL, NULL);
    wtime_recv += GetWtime() - t0;
    for(int i=0; i<n_tot; i++){
        acc[i].x = force_h[i].ax;
        acc[i].y = force_h[i].ay;
        acc[i].z = force_h[i].az;
        pot[i]   = force_h[i].pot;
    }
}


#endif

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
    std::cerr<<"r_max= "<<sqrt(r_max_sq)<<std::endl;
}


template<class Tpsys>
void setParticlesPlummer(Tpsys & psys,
                         const PS::S64 n_glb,
                         PS::S64 & n_loc){  

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
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
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
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
    std::cerr<<"ekin_loc="<<ekin_loc<<std::endl;
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

struct Energy{
    PS::F64 tot_;
    PS::F64 kin_;
    PS::F64 pot_;
    Energy(){
        tot_ = kin_ = pot_ = 0.0;
    }
    void clear(){
        tot_ = kin_ = pot_ = 0.0;
    }
    void calc(const PS::F64 mass[], const PS::F64vec vel[], const PS::F64 pot[], const PS::S32 n){
        clear();
        for(PS::S32 i=0; i<n; i++){
            kin_ += 0.5 * mass[i] * vel[i] * vel[i];
            pot_ += 0.5 * mass[i] * pot[i];
        }
        tot_ = kin_ + pot_;
    }
    void dump(){
        std::cerr<<"tot="<<tot_<<" kin="<<kin_<<" pot="<<pot_<<std::endl;
    }
};

void calcEnergy(){

}

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);


    PS::F32 time_end = 10.0;
    PS::F32 dt = 1.0 / 128.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_tot = 1024;
    PS::S32 c;
    sprintf(dir_name,"./result");
    opterr = 0;
    while((c=getopt(argc,argv,"i:o:d:D:T:N:hs:")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr << "time_end = " << time_end << std::endl;
            break;
        case 's':
            dt = atof(optarg);
            std::cerr << "time_step = " << dt << std::endl;
            break;
        case 'd':
            dt_diag = atof(optarg);
            std::cerr << "dt_diag = " << dt_diag << std::endl;
            break;
        case 'D':
            dt_snap = atof(optarg);
            std::cerr << "dt_snap = " << dt_snap << std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            std::cerr << "n_tot = " << n_tot << std::endl;
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
    InitializeDEVICE();
    const PS::F64 m_tot = 1.0;

    PS::F64 * mass = NULL;
    PS::F64vec * pos = NULL;
    PS::F64vec * vel = NULL;
    PS::F64vec * acc = new PS::F64vec[n_tot];
    PS::F64 * pot = new PS::F64[n_tot];
    MakePlummerModel(m_tot, n_tot, n_tot, mass, pos, vel);
    const PS::F32 eps = 4.0 / n_tot;
    const PS::F32 eps2 = eps * eps;

    calcForce(mass, pos, n_tot, eps2, acc, pot);
    PS::S32 n_loop_max = 10;
    //std::chrono::system_clock::time_point wtime_start = std::chrono::system_clock::now();
    double wtime_start = GetWtime();
    for(PS::S32 i=0; i<n_loop_max; i++){
        calcForce(mass, pos, n_tot, eps2, acc, pot);
    }
    double wtime_end = GetWtime();
    PS::F64 wtime_exe = (wtime_end - wtime_start) / n_loop_max;
    wtime_send /= n_loop_max;
    wtime_recv /= n_loop_max;
    wtime_force /= n_loop_max;
    //std::chrono::system_clock::time_point wtime_end = std::chrono::system_clock::now();
    //PS::F64 time_exe = (std::chrono::duration_cast<std::chrono::microseconds>(wtime_end - wtime_start).count()) / n_loop_max; // micro seconds

    for(PS::S32 i=0; i<n_tot; i++) pot[i] += mass[i] / eps;
#if 0
    if(n_tot == 1024){
        for(PS::S32 i=0; i<n_tot; i++){
            std::cout<<"i= "<<i<<" acc[i]= "<<acc[i]<<" pot[i]= "<<pot[i]<<std::endl;
        }
    }
#endif

    Energy eng;
    eng.calc(mass, vel, pot, n_tot);
    eng.dump();
    PS::S64 op_interaction = 38;

    PS::F64 speed = n_tot * n_tot * op_interaction / wtime_exe;
    std::cout<<"n_tot= "<<n_tot<<" wtime_exe= "<<wtime_exe
             <<" wtime_send= "<<wtime_send
             <<" wtime_recv= "<<wtime_recv
             <<" wtime_force= "<<wtime_force
             <<" speed= "<<speed*1e-9<<" [Glops]"<<std::endl;
    PS::Finalize();
    return 0;
}
