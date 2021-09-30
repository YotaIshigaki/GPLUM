#pragma once

#ifdef ENABLE_PEZY
#define MULTI_WALK

#include<pzcl/pzcl_ocl_wrapper.h>
#include<PZSDKHelper.h>
#include"./pezy/class_device.hpp"

#define NPEZY 4

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
const int N_WALK_LIMIT = 200;
const int NI_LIMIT = 8192 * N_WALK_LIMIT;
const int NJ_LIMIT = 131072 * N_WALK_LIMIT;
const int N_THREAD_MAX = 8192;

int j_disp_h[N_WALK_LIMIT+2];
EpiDev epi_h[NI_LIMIT];
EpjDev epj_h[NJ_LIMIT];
ForceDev force_h[NI_LIMIT];

cl_mem j_disp_d = NULL;
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
std::vector<cl_device_id> device_id_lists;

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
    ret = clSetKernelArg(kernel, 5, sizeof(int),    (void*)&N_THREAD_MAX);
}

PS::F64 WTIME_SET_OFFSET = 0.0;
PS::F64 WTIME_SET = 0.0;

PS::S32 DispatchKernel(const PS::S32          tag,
		       const PS::S32          n_walk,
		       const EP                *epi[],
		       const PS::S32          n_epi[],
		       const EP                *epj[],
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
  
  ret = clEnqueueWriteBuffer(command_queue, j_disp_d, CL_TRUE, 0, (n_walk+2)*sizeof(int), j_disp_h, 0, NULL, NULL);
  
  ni_tot = 0;
  int nj_tot = 0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<n_epi[iw]; i++){
      epi_h[ni_tot].pos.x   = (float)epi[iw][i].pos.x;
      epi_h[ni_tot].pos.y   = (float)epi[iw][i].pos.y;
      epi_h[ni_tot].pos.z   = (float)epi[iw][i].pos.z;
      epi_h[ni_tot].id      = epi[iw][i].id;
      epi_h[ni_tot].walk    = iw;
      ni_tot++;
    }
    for(int j=0; j<n_epj[iw]; j++){
      epj_h[nj_tot].pos.x  = epj[iw][j].pos.x;
      epj_h[nj_tot].pos.y  = epj[iw][j].pos.y;
      epj_h[nj_tot].pos.z  = epj[iw][j].pos.z;
      epj_h[nj_tot].id     = epj[iw][j].id;
      epj_h[nj_tot].charge = epj[iw][j].charge;
      nj_tot++;
    }
  }
  ret = clEnqueueWriteBuffer(command_queue, epi_d,   CL_TRUE, 0, (ni_tot)*sizeof(EpiDev),  epi_h, 0, NULL, NULL);
  ret = clEnqueueWriteBuffer(command_queue, epj_d,   CL_TRUE, 0, (nj_tot)*sizeof(EpjDev),  epj_h, 0, NULL, NULL);

  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&j_disp_d);
  ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&epi_d);
  ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&epj_d);
  ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*)&force_d);
  ret = clSetKernelArg(kernel, 4, sizeof(int),    (void*)&ni_tot);

  size_t work_size = N_THREAD_MAX;
#if 1
  ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
#else
  for(int ith=0; ith<8; ith++){
    for(int ipr=0; ipr<N_THREAD_MAX/8; ipr++){
      cpu_ForceKernel(j_disp_h, epi_h, epj_h, force_h, ni_tot, N_THREAD_MAX, ith, ipr);
    }
  }
#endif
  return 0;
}

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       Force        *force[]
		       ){
  cl_int ret;
  int ni_tot = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += ni[k];
  }
  ret = clEnqueueReadBuffer(command_queue, force_d, CL_TRUE, 0, ni_tot*sizeof(ForceDev), force_h, 0, NULL, NULL);
  int n_cnt = 0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<ni[iw]; i++){
      force[iw][i].acc.x = force_h[n_cnt].f.x;
      force[iw][i].acc.y = force_h[n_cnt].f.y;
      force[iw][i].acc.z = force_h[n_cnt].f.z;
      force[iw][i].pot   = force_h[n_cnt].u;
      n_cnt++;
    }
  }
  return 0;
}

#endif
