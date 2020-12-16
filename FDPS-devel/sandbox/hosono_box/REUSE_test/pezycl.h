#pragma once
#include <pzcl/pzcl_ocl_wrapper.h>
#include <PZSDKHelper.h>

const int N_WALK_LIMIT = 200;
const int NI_LIMIT = 10000 * N_WALK_LIMIT;
const int NJ_LIMIT = 10000 * N_WALK_LIMIT;
const int N_THREAD_MAX = 8192;

struct PezyDevice{
	static const int NumDeviceMax = 4;
	//on device vars.
	struct{
		cl_mem    j_disp;
		cl_mem    epi;
		cl_mem    epj;
		cl_mem    force;
		cl_kernel kernel;
	}dens, drvt, hydr, grav;
	//device info.
	cl_platform_id   platform_id;
	cl_uint          num_of_platforms, num_of_devices;
	cl_device_id     device_id[NumDeviceMax];
	cl_context       context;
	cl_command_queue cmd_queue;
	cl_program       program;

	std::vector<cl_device_id> device_id_list;

	void initialize();
};


