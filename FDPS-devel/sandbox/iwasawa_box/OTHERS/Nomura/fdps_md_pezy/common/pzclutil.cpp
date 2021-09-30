/**
 *		PZCL Utilities
 *
 *	Copyright (c) 2012 PEZY Computing, K.K.
 */
#include "pzclutil.h"


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
cl_program PZCLUtil::CreateProgram(
	cl_context context,
	std::vector<cl_device_id> &device_id_lists,
	const char* bin_name
	)
{
	cl_program program = NULL;
	char* pBin = NULL;

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
		cl_int result;

		program = clCreateProgramWithBinary( context, (cl_uint)device_id_lists.size(), &device_id_lists[0], &length, listBin, &binary_status, &result);
		if(program == NULL)
		{
			printf("clCreateProgramWithBinary failed, %d\n", result);
			goto leaving;
		}
	}
leaving:
	if(pBin)
	{
		PZSDK_ALIGNED_FREE(pBin);
	}

	return program;
}

/** 
 * SetKernelArgs
 */
cl_int PZCLUtil::SetKernelArgs( cl_kernel kernel, const std::vector<KernelArg> &vArg)
{
	cl_int result = CL_SUCCESS;
	
	size_t count = vArg.size();
	for(size_t i = 0; i < count; i++)
	{
		if( (result = clSetKernelArg(kernel, (cl_uint)i, vArg[i].size, vArg[i].value ) ) != CL_SUCCESS)
		{
			printf("clSetKernelArg ARG%d failed, %d", i, result);
			return result;
		}
	}
	return result;
}

/**
 * getTime
 */
cl_ulong PZCLUtil::getTime( cl_event event )
{
	cl_ulong start = 0;
	cl_ulong end   = 0;
	cl_int result;
	
	if((result = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL)) != CL_SUCCESS )
	{
		printf("clGetEventProfilingInfo failed, %d\n", result);
		return 0;
	}
	if((result = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL)) != CL_SUCCESS )
	{
		printf("clGetEventProfilingInfo failed, %d\n", result);
		return 0;
	}

	return (end - start);
}



