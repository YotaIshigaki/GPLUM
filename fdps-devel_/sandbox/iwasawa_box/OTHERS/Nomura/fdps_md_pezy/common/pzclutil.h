/*
 *  PZCL Utilities
 *
 *	Copyright (c) 2012 PEZY Computing, K.K.
 */
#pragma once

#ifdef _MSC_VER
	#pragma warning( disable:4996)
#endif

#include <PZSDKHelper.h>
#include <pzcl/pzcl_ocl_wrapper.h> // OpenCL wrapper

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>

#include "pzcperf.h"

namespace PZCLUtil
{
	/**
	 * @brief Create program object
	 */
	cl_program CreateProgram(cl_context context, std::vector<cl_device_id> &device_id_lists, const char* bin_name);

	/**
	 * Container for the kernel argument
	 */
	struct KernelArg
	{
		KernelArg( size_t inSize, void* inValue)
		 : size(inSize), value(inValue)
		{
		}
		size_t size;
		void*  value;
	};

	/**
	 * @brief Set kernel arguments
	 */
	 cl_int SetKernelArgs( cl_kernel kernel, const std::vector<KernelArg> &vArg);

	 /**
	  * @brief getTime
	  * @note  CL_QUEUE_PROFILING_ENABLE
	  */
	 cl_ulong getTime( cl_event event );

}

