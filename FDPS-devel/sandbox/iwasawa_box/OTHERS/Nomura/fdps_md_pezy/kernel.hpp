#pragma once


#define MULTI_WALK

#include<pzcl/pzcl_ocl_wrapper.h>
#include<PZSDKHelper.h>
#include"./pezy/class_device.hpp"

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
cl_program CreateProgram(cl_context context,
			 std::vector<cl_device_id> &device_id_lists,
			 const char* bin_name
			 ){
  
  cl_program program = NULL;
  char* pBin = NULL;
  cl_int result;
  
  size_t sizeFile = GetFileSize( bin_name );
  if( sizeFile == 0 )
    {
      goto leaving;
    }
  
  PZSDK_ALIGNED_ALLOC(pBin, sizeFile, 8 /*8 byte alignment*/);
  if( pBin == NULL){
    printf("out of host memory\n");
    goto leaving;
  }

  if( !LoadFile(bin_name, sizeFile, pBin)){
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
const int N_WALK_LIMIT = 200;
const int NI_LIMIT = 8192 * N_WALK_LIMIT;
const int NJ_LIMIT = 131072 * N_WALK_LIMIT;
const int N_THREAD_MAX = 8192;

int j_disp_h[N_WALK_LIMIT+2];
EpiDev epi_h[NI_LIMIT];
EpjDev epj_h[NJ_LIMIT];

ForceDev force_h[NI_LIMIT];
//ForceDev buf_h[N_THREAD_MAX];
ForceDev buf_h[2*NI_LIMIT];

cl_mem j_disp_d = NULL;
cl_mem epi_d = NULL;
cl_mem epj_d = NULL;
cl_mem force_d = NULL;
cl_mem buf_d = NULL;

cl_platform_id platform_id = NULL;
cl_uint ret_num_platforms;
cl_device_id device_id[N_DEV_MAX];
cl_uint ret_num_devices;
cl_context context = NULL;
cl_command_queue command_queue = NULL;
cl_program program = NULL;
cl_kernel kernel = NULL;
cl_kernel kernelPerf = NULL;
std::vector<cl_device_id> device_id_lists;

#define NPEZY 4

// ------- Profiler --------
#include<common/pzclutil.h>

//#define CACHE_PROFILE
#ifdef CACHE_PROFILE
void printProfileCache(pzcl_profile_cache & profile)
{
  fprintf(stderr,"ReadRequest  = %d\n", profile.read_request);
  fprintf(stderr,"ReadHit      = %d\n", profile.read_hit);
  fprintf(stderr,"WriteRequest = %d\n", profile.write_request);
  fprintf(stderr,"WriteHit     = %d\n", profile.write_hit);
}

void printProfileCacheStatistics(pzcl_profile_cache_stats & stats)
{
  fprintf(stderr,"ReadHitRate  = %10.10lf\n", stats.read_hit_rate);
  fprintf(stderr,"WriteHitRate = %10.10lf\n", stats.write_hit_rate);
}
#endif

//#define KERNEL_PROFILE
#ifdef KERNEL_PROFILE
PZCPerformance GetPerformance(int index){
  PZCPerformance perf;
  cl_int result = 0;
  cl_mem mem = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(PZCPerformance), NULL, &result );
  if( mem == NULL){
    printf("clCreateBuffer failed, %d", result);
  }

  // Setup Kernel argument
  result = clSetKernelArg(kernelPerf, 0, sizeof(cl_mem), (void*)&mem);
  result = clSetKernelArg(kernelPerf, 1, sizeof(int), &index);

  size_t global_work_size = 128;
  if( (result = clEnqueueNDRangeKernel(command_queue, kernelPerf, 1 /* must be 1D*/, NULL, &global_work_size, NULL, 0, NULL, NULL ) ) != CL_SUCCESS ){
      printf("clEnqueueNDRangeKernel failed, %d\n", result);
  }

  result = clEnqueueReadBuffer( command_queue, mem, CL_TRUE, 0, sizeof(PZCPerformance), &perf, 0, NULL, NULL);
  if( result != CL_SUCCESS){
    printf("clEnqueueReadBuffer failed %d\n", result);
  }
  return perf;
}
#endif
// ------- Profiler --------

void InitializeDEVICE(){
    const PS::S32 my_rank = PS::Comm::getRank();
    //std::cout<<"my_rank="<<my_rank<<std::endl;
    cl_int ret;
    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    std::cerr<<"platform_id="<<platform_id<<" ret_num_platforms="<<ret_num_platforms<<" ret="<<ret<<std::endl;

    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, N_DEV_MAX, device_id, &ret_num_devices);
    std::cerr<<"device_id[" << my_rank << "]=" << device_id[my_rank%NPEZY] <<" ret_num_devices="<<ret_num_devices<<" ret="<<ret<<std::endl;

    //context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
    context = clCreateContext(NULL, 1, &device_id[my_rank%NPEZY], NULL, NULL, &ret);
    printf("Context returns %d at rank %d\n",ret, my_rank);

    //command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
    command_queue = clCreateCommandQueue(context, device_id[my_rank%NPEZY], 0, &ret);

    j_disp_d = clCreateBuffer(context, CL_MEM_READ_WRITE, (N_WALK_LIMIT+2)*sizeof(int),   NULL, &ret);
    epi_d    = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT*sizeof(EpiDev),   NULL, &ret);
    epj_d    = clCreateBuffer(context, CL_MEM_READ_WRITE, NJ_LIMIT*sizeof(EpjDev),   NULL, &ret);
    force_d  = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT*sizeof(ForceDev), NULL, &ret);
    //buf_d  = clCreateBuffer(context, CL_MEM_READ_WRITE, N_THREAD_MAX*sizeof(ForceDev), NULL, &ret);
    buf_d  = clCreateBuffer(context, CL_MEM_READ_WRITE, 2*NI_LIMIT*sizeof(ForceDev), NULL, &ret);

    //device_id_lists.push_back(device_id);
    device_id_lists.push_back(device_id[my_rank%NPEZY]);
    program = CreateProgram(context, device_id_lists, "./pezy/kernel.sc32/kernel.pz");
    if(program == NULL){
      std::cerr<<"can't create program"<<std::endl;
      exit(1);
    }
    kernel = clCreateKernel(program, "ForceKernel", &ret);
    if(kernel == NULL){
      std::cerr<<"can't create kernel"<<std::endl;
      exit(1);
    }
    kernelPerf = clCreateKernel(program, "GetPerformance", &ret);
    if(kernel == NULL){
      std::cerr<<"can't create kernelPerf"<<std::endl;
      exit(1);
    }

#if 1  // use local memory or not
    typedef CL_API_ENTRY cl_int (CL_API_CALL *pfnPezyExtSetPerThreadStackSize) (pzcl_kernel kernel, size_t size);
    pfnPezyExtSetPerThreadStackSize clExtSetPerThreadStackSize = (pfnPezyExtSetPerThreadStackSize)clGetExtensionFunctionAddress("pezy_set_per_thread_stack_size");
    assert(clExtSetPerThreadStackSize);
    //size_t per_thread_stack = 0x500;
    size_t per_thread_stack = 0x400;
    ret = clExtSetPerThreadStackSize(kernel, per_thread_stack);
    assert(ret == CL_SUCCESS);
#endif
#if 0 // write memory without cache read
    typedef CL_API_ENTRY pzcl_int (CL_API_CALL *pfnPezyExtSetCacheWriteBuffer)(pzcl_context context, size_t index, bool enable);
    pfnPezyExtSetCacheWriteBuffer clExtSetCacheWriteBuffer = (pfnPezyExtSetCacheWriteBuffer)clGetExtensionFunctionAddress("pezy_set_cache_writebuffer");
    if (!clExtSetCacheWriteBuffer) fprintf(stderr,"error: write memory setting\n");
    ret = clExtSetCacheWriteBuffer(context, 0, CL_TRUE);
    assert(ret == CL_SUCCESS);
#endif
}

PS::F64 WTIME_SET_OFFSET = 0.0;
PS::F64 WTIME_SET = 0.0;

PS::S32 DispatchKernelWithSP(const PS::S32          tag,
                             const PS::S32          n_walk,
                             const EPILJ             *epi[],
                             const PS::S32          n_epi[],
                             const EPJLJ             *epj[],
                             const PS::S32          n_epj[]
                             ){
  cl_int ret;
  assert(n_walk <= N_WALK_LIMIT);
  //const float eps2 = FPGrav::eps * FPGrav::eps;
  PS::S32 ni_tot = 0;
  timer.beg(Profile::TRANSLATE);
  j_disp_h[0] = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += n_epi[k];
    j_disp_h[k+1] = j_disp_h[k] + n_epj[k];
  }
  j_disp_h[n_walk+1] = j_disp_h[n_walk];

  assert(ni_tot < NI_LIMIT);
  assert(j_disp_h[n_walk] < NJ_LIMIT);
  timer.end(Profile::TRANSLATE);

  timer.beg(Profile::TRANSFER);
  ret = clEnqueueWriteBuffer(command_queue, j_disp_d, CL_TRUE, 0, (n_walk+2)*sizeof(int), j_disp_h, 0, NULL, NULL);
  timer.end(Profile::TRANSFER);

  timer.beg(Profile::TRANSLATE);
  ni_tot = 0;
  int nj_tot = 0;
#if 1 //SoA(0) or AoS(1)
  for(int iw=0; iw<n_walk; iw++){
    //printf("%d %d\n",n_epi[iw],n_epj[iw]);
    for(int i=0; i<n_epi[iw]; i++){
      epi_h[ni_tot].px      = epi[iw][i].pos.x;
      epi_h[ni_tot].py      = epi[iw][i].pos.y;
      epi_h[ni_tot].pz      = epi[iw][i].pos.z;
      epi_h[ni_tot].id_walk = iw;
      ni_tot++;
    }
    for(int j=0; j<n_epj[iw]; j++){
      epj_h[nj_tot].px     = epj[iw][j].pos.x;
      epj_h[nj_tot].py     = epj[iw][j].pos.y;
      epj_h[nj_tot].pz     = epj[iw][j].pos.z;
      //epj_h[nj_tot].mass   = epj[iw][j].mass;
      nj_tot++;
    }
  }
#else
  float *epi_tmp = (float*)epi_h;
  float *epj_tmp = (float*)epj_h;
  for(int iw=0; iw<n_walk; iw++){
    ni_tot+=n_epi[iw];
    nj_tot+=n_epj[iw];
  }
  int icount=0,jcount=0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<n_epi[iw]; i++){
      epi_tmp[icount+0*ni_tot] = epi[iw][i].pos.x;
      epi_tmp[icount+1*ni_tot] = epi[iw][i].pos.y;
      epi_tmp[icount+2*ni_tot] = epi[iw][i].pos.z;
      *((int*)&epi_tmp[icount+3*ni_tot]) = iw;
      //printf("%d %d %d %d %d\n",icount,iw,j_disp_h[iw],j_disp_h[iw+1],j_disp_h[iw+1]-j_disp_h[iw]);
      icount++;
    }
    for(int j=0; j<n_epj[iw]; j++){
      epj_tmp[jcount+0*nj_tot] = epj[iw][j].pos.x;
      epj_tmp[jcount+1*nj_tot] = epj[iw][j].pos.y;
      epj_tmp[jcount+2*nj_tot] = epj[iw][j].pos.z;
      //epj_h[nj_tot].mass = epj[iw][j].mass;
      jcount++;
    }
  }
#endif // SoA or AoS
  timer.end(Profile::TRANSLATE);

  timer.beg(Profile::TRANSFER);
  ret = clEnqueueWriteBuffer(command_queue, epi_d,   CL_TRUE, 0, (ni_tot)*sizeof(EpiDev),  epi_h, 0, NULL, NULL);
  ret = clEnqueueWriteBuffer(command_queue, epj_d,   CL_TRUE, 0, (nj_tot)*sizeof(EpjDev),  epj_h, 0, NULL, NULL);
  timer.end(Profile::TRANSFER);

  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&j_disp_d);
  ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&epi_d);
  ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&epj_d);
  ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*)&force_d);
  ret = clSetKernelArg(kernel, 4, sizeof(int),    (void*)&ni_tot);
  //ret = clSetKernelArg(kernel, 5, sizeof(int),    (void*)&nj_tot);
  ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void*)&buf_d);

#ifdef CACHE_PROFILE
      pfnPezyExtSetProfile clSetProfile = (pfnPezyExtSetProfile)clGetExtensionFunctionAddress("pezy_set_profile");
    assert(clSetProfile);

    pfnPezyExtGetProfileCache clGetProfileCache = (pfnPezyExtGetProfileCache)clGetExtensionFunctionAddress("pezy_get_profile_cache");
    assert(clGetProfileCache);

    pfnPezyExtGetProfileCacheStatistics clGetProfileCacheStatistics = (pfnPezyExtGetProfileCacheStatistics)clGetExtensionFunctionAddress("pezy_get_profile_cache_statistics");
    assert(clGetProfileCacheStatistics);

    pzcl_profile_cache profile;
    pzcl_profile_cache_stats stats = {sizeof(stats)};
    pzcl_int level = 0;
    size_t id = 0;
    ret = clSetProfile(context, 0, CL_TRUE);
    if (ret != CL_SUCCESS) {
        printf("clSetProfile failed %d\n", ret);
    }
#endif //CACHE_PROFILE
  size_t work_size = N_THREAD_MAX;
  timer.beg(Profile::KERNEL);
  ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
#ifdef TUNING
  clFinish(command_queue);

#ifdef CACHE_PROFILE
  level = PZCL_EXT_PROFILE_CACHE_L1; //L2,L3
  ret = clGetProfileCache(context, 0, level, id, &profile);
  if (ret != CL_SUCCESS) {
    printf("clGetProfileCache failed %d\n", ret);
  }
  ret = clGetProfileCacheStatistics(context, 0, level, &stats);
  if (ret != CL_SUCCESS) {
    printf("clGetProfileCacheStatistics failed %d\n", ret);
  }
  fprintf(stderr,"L1:\n");
  printProfileCache(profile);
  printProfileCacheStatistics(stats);

  level = PZCL_EXT_PROFILE_CACHE_L2; //L2,L3
  ret = clGetProfileCache(context, 0, level, id, &profile);
  if (ret != CL_SUCCESS) {
    printf("clGetProfileCache failed %d\n", ret);
  }
  ret = clGetProfileCacheStatistics(context, 0, level, &stats);
  if (ret != CL_SUCCESS) {
    printf("clGetProfileCacheStatistics failed %d\n", ret);
  }
  fprintf(stderr,"L2:\n");
  printProfileCache(profile);
  printProfileCacheStatistics(stats);

  level = PZCL_EXT_PROFILE_CACHE_L3; //L2,L3
  ret = clGetProfileCache(context, 0, level, id, &profile);
  if (ret != CL_SUCCESS) {
    printf("clGetProfileCache failed %d\n", ret);
  }
  ret = clGetProfileCacheStatistics(context, 0, level, &stats);
  if (ret != CL_SUCCESS) {
    printf("clGetProfileCacheStatistics failed %d\n", ret);
  }
  fprintf(stderr,"L3:\n");
  printProfileCache(profile);
  printProfileCacheStatistics(stats);
#endif // CACHE_PROFILE
#ifdef KERNEL_PROFILE
  static bool isPerf = true;
  if(isPerf){
    PZCPerformance perf = GetPerformance(0);

    unsigned int perf_count = perf.Perf();
    unsigned int stall      = perf.Stall();
    unsigned int wait       = perf.Wait();
    const double clock = 690e6;
    double sec  = perf_count / clock;

    printf("Perf : %d (%e sec)\n", perf_count, (perf_count != 0) ? sec : 0.0);
    printf("Stall: %d (%5.3f %%)\n", stall, (perf_count != 0) ? stall * 100.0 / perf_count: 0.0 );
    printf("Wait : %d (%5.3f %%)\n", wait, (perf_count != 0) ? wait * 100.0 / perf_count: 0.0 );
    isPerf = false;
  }
#endif //KERNEL_PROFILE
#endif //TUNING
  timer.end(Profile::KERNEL);
  return 0;
}

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       ForceLJ   *force[]
		       ){
  cl_int ret;
  int ni_tot = 0;
  timer.beg(Profile::TRANSFER);
  for(int k=0; k<n_walk; k++){
    ni_tot += ni[k];
  }
  ret = clEnqueueReadBuffer(command_queue, force_d, CL_TRUE, 0, ni_tot*sizeof(ForceDev), force_h, 0, NULL, NULL);
  ret = clEnqueueReadBuffer(command_queue, buf_d, CL_TRUE, 0, 2*ni_tot*sizeof(ForceDev), buf_h, 0, NULL, NULL);
  //ret = clEnqueueReadBuffer(command_queue, buf_d, CL_TRUE, 0, N_THREAD_MAX*sizeof(ForceDev), buf_h, 0, NULL, NULL);
  timer.end(Profile::TRANSFER);
  timer.beg(Profile::TRANSLATE);
  int n_cnt = 0;
  for(int iw=0; iw<n_walk; iw++){
    //if(iw<512) printf("%d: %e %e %e %e\n",n_cnt,force_h[n_cnt].ax,force_h[n_cnt].ay,force_h[n_cnt].az,force_h[n_cnt].pot);
    for(int i=0; i<ni[iw]; i++){
      force[iw][i].acc.x = force_h[n_cnt].ax;
      force[iw][i].acc.y = force_h[n_cnt].ay;
      force[iw][i].acc.z = force_h[n_cnt].az;
      force[iw][i].pot   = force_h[n_cnt].pot;
      n_cnt++;
    }
  }
  //printf("%e: %e + %e = %e\n",force_h[0].ax,buf_h[0].ax,buf_h[ni_tot].ax,buf_h[0].ax+buf_h[ni_tot].ax);
  timer.end(Profile::TRANSLATE);
  return 0;
}

PS::S32 DispatchKernelWithSPcpu(const PS::S32          tag,
				const PS::S32          n_walk,
				const EPILJ             *epi[],
				const PS::S32          n_epi[],
				const EPJLJ             *epj[],
				const PS::S32          n_epj[]
				){
  assert(n_walk <= N_WALK_LIMIT);
  cl_int ret;
  //const float eps2 = FPGrav::eps * FPGrav::eps;
  PS::S32 ni_tot = 0;
  j_disp_h[0] = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += n_epi[k];
    j_disp_h[k+1] = j_disp_h[k] + n_epj[k];
  }
  j_disp_h[n_walk+1] = j_disp_h[n_walk];
  assert(ni_tot < NI_LIMIT);
  assert(j_disp_h[n_walk] < NJ_LIMIT);
  int index = 0;
  const float cut_off2 = 4.5f*4.5f;
  //const float cut_off2 = 9.f*9.f;
  //const float cut_off2 = 13.5f*13.5f;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<n_epi[iw]; i++){
      double ax, ay, az, pot;
      ax = ay = az = pot = 0.0f;
      const EPILJ ip = epi[iw][i];
      const int j_head = j_disp_h[iw];
      const int j_tail = j_disp_h[iw+1];
      for(int j=0; j<n_epj[iw]; j++){
	const EPJLJ jp = epj[iw][j];
	const float dx = ip.pos.x - jp.pos.x;
	const float dy = ip.pos.y - jp.pos.y;
	const float dz = ip.pos.z - jp.pos.z;
	const float r2 = (dx*dx + dy*dy) + dz*dz;
#ifdef TUNING
	flop += 8;
#endif
	if(r2 <= cut_off2 && r2!=0.0){
	  const float r2_inv = 1.f / r2;// 1 + 2 + 1 + 3 + 4 = 11 flop
	  const float r6_inv  = r2_inv*r2_inv*r2_inv;
	  const float r12_inv = r6_inv*r6_inv;
	  pot += 4.0f*(r12_inv - r6_inv);
	  const float dphi = (48.f*r12_inv - 24.f*r6_inv)*r2_inv; 
	  ax  += dphi * dx; // 6 flop
	  ay  += dphi * dy;
	  az  += dphi * dz;
#ifdef TUNING
	  flop += 17;
#endif
	}
      }
      force_h[index].pot = pot;
      force_h[index].ax = ax;
      force_h[index].ay = ay;
      force_h[index].az = az;
      index++;
    }
  }
}

PS::S32 RetrieveKernelcpu(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       ForceLJ   *force[]
		       ){
  int ni_tot = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += ni[k];
  }
  int n_cnt = 0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<ni[iw]; i++){
      force[iw][i].acc.x = force_h[n_cnt].ax;
      force[iw][i].acc.y = force_h[n_cnt].ay;
      force[iw][i].acc.z = force_h[n_cnt].az;
      force[iw][i].pot   = force_h[n_cnt].pot;
      n_cnt++;
    }
  }
  return 0;
}

struct CalcForceEpEpPEZY{
  void operator () (const EPILJ * ep_i,
                    const PS::S32 n_ip,
                    const EPJLJ * ep_j,
                    const PS::S32 n_jp,
                    ForceLJ * force){
    cl_int ret;
    timer.beg(Profile::TRANSLATE);
    for(int i=0;i<n_ip;i++){
      epi_h[i].px = ep_i[i].pos.x;
      epi_h[i].py = ep_i[i].pos.y;
      epi_h[i].pz = ep_i[i].pos.z;
      epi_h[i].id_walk = ep_i[i].id;
    }
    for(int i=0;i<n_jp;i++){
      epj_h[i].px      = ep_j[i].pos.x;
      epj_h[i].py      = ep_j[i].pos.y;
      epj_h[i].pz      = ep_j[i].pos.z;
      //epj_h[i].id_walk = ep_j[i].id;
    }
    timer.end(Profile::TRANSLATE);

    timer.beg(Profile::TRANSFER);
    ret = clEnqueueWriteBuffer(command_queue, epi_d,   CL_TRUE, 0, (n_ip)*sizeof(EpiDev),  epi_h, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, epj_d,   CL_TRUE, 0, (n_jp)*sizeof(EpjDev),  epj_h, 0, NULL, NULL);
    timer.end(Profile::TRANSFER);

    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&epi_d);
    ret = clSetKernelArg(kernel, 1, sizeof(int),    (void*)&n_ip);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&epj_d);
    ret = clSetKernelArg(kernel, 3, sizeof(int),    (void*)&n_jp);
    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem),    (void*)&force_d);

#ifdef CACHE_PROFILE
    pfnPezyExtSetProfile clSetProfile = (pfnPezyExtSetProfile)clGetExtensionFunctionAddress("pezy_set_profile");
    assert(clSetProfile);

    pfnPezyExtGetProfileCache clGetProfileCache = (pfnPezyExtGetProfileCache)clGetExtensionFunctionAddress("pezy_get_profile_cache");
    assert(clGetProfileCache);

    pfnPezyExtGetProfileCacheStatistics clGetProfileCacheStatistics = (pfnPezyExtGetProfileCacheStatistics)clGetExtensionFunctionAddress("pezy_get_profile_cache_statistics");
    assert(clGetProfileCacheStatistics);

    pzcl_profile_cache profile;
    pzcl_profile_cache_stats stats = {sizeof(stats)};
    pzcl_int level = 0;
    size_t id = 0;
    ret = clSetProfile(context, 0, CL_TRUE);
    if (ret != CL_SUCCESS) {
      printf("clSetProfile failed %d\n", ret);
    }
#endif

    size_t work_size = N_THREAD_MAX;
    timer.beg(Profile::KERNEL);
    ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
#ifdef TUNING
    clFinish(command_queue);
#ifdef CACHE_PROFILE
    ret = clGetProfileCache(context, 0, level, id, &profile);
    if (ret != CL_SUCCESS) {
      printf("clGetProfileCache failed %d\n", ret);
    }
    ret = clGetProfileCacheStatistics(context, 0, level, &stats);
    if (ret != CL_SUCCESS) {
      printf("clGetProfileCacheStatistics failed %d\n", ret);
    }
    printProfileCache(profile);
    printProfileCacheStatistics(stats);
#endif //CACHE_PROFILE
#endif //TUNING
    timer.end(Profile::KERNEL);

    timer.beg(Profile::TRANSFER);
    ret = clEnqueueReadBuffer(command_queue, force_d, CL_TRUE, 0, n_ip*sizeof(ForceDev), force_h, 0, NULL, NULL);
    timer.end(Profile::TRANSFER);

    timer.beg(Profile::TRANSLATE);
    for(int i=0; i<n_ip;i++){
      force[i].acc.x = force_h[i].ax;
      force[i].acc.y = force_h[i].ay;
      force[i].acc.z = force_h[i].az;
      force[i].pot   = force_h[i].pot;
    }
    timer.end(Profile::TRANSLATE);
  }
  CalcForceEpEpPEZY(const PS::F32 _cut_off2) : cut_off2(_cut_off2) {}
private:
  PS::F32 cut_off2;
};
