#include <sys/time.h>
#include <particle_simulator.hpp>
#include "EoS.h"
#include "class_platform.hpp"
#include "pzc/class_device.hpp"
#include "pezycl.h"

PezyDevice device;

double wtime_WB = 0;
double wtime_SKA = 0;
double wtime_NDRK = 0;
double wtime_RB = 0;
double wtime_CP = 0;

struct dens_host_t{
	int j_disp[N_WALK_LIMIT+2];
	Dens::EpiDev epi[NI_LIMIT];
	Dens::EpjDev epj[NJ_LIMIT];
	Dens::ForceDev force[NI_LIMIT];
} dens_host;

struct drvt_host_t{
	int j_disp[N_WALK_LIMIT+2];
	Drvt::EpiDev epi[NI_LIMIT];
	Drvt::EpjDev epj[NJ_LIMIT];
	Drvt::ForceDev force[NI_LIMIT];
} drvt_host;

struct hydr_host_t{
	int j_disp[N_WALK_LIMIT+2];
	Hydr::EpiDev epi[NI_LIMIT];
	Hydr::EpjDev epj[NJ_LIMIT];
	Hydr::ForceDev force[NI_LIMIT];
} hydr_host;

struct grav_host_t{
	int j_disp[N_WALK_LIMIT+2];
	Grav::EpiDev epi[NI_LIMIT];
	Grav::EpjDev epj[NJ_LIMIT];
	Grav::ForceDev force[NI_LIMIT];
} grav_host;

static bool LoadFile(const char* name, size_t size, char* pData){
	FILE* fp = fopen(name, "rb");
	if(fp == NULL){
		printf("can not open %s\n", name);
		return false;
	}

	if(size == 0 || pData == NULL){
		printf("invalid params %s\n", __FUNCTION__);
		return false;
	}

	size_t size_ret = fread(pData, sizeof(char), size, fp);
	fclose(fp);

	if(size_ret != size){
		printf("can not read requested size\n");
		return false;
	}
	return true;
}

static size_t GetFileSize(const char* name){
	FILE* fp = fopen(name, "rb");
	if(fp == NULL){
		printf("can not open %s", name);
		return 0;
	}
	fseek(fp, 0, SEEK_END);
	size_t size = ftell(fp);
	fclose(fp);
	return size;
}

cl_program CreateProgram(cl_context context, std::vector<cl_device_id> &device_id_lists, const char* bin_name){
	cl_program program = NULL;
	char* pBin = NULL;
	cl_int result;

	size_t sizeFile = GetFileSize(bin_name);
	if(sizeFile == 0){
		goto leaving;
	}

	PZSDK_ALIGNED_ALLOC(pBin, sizeFile, 8 /*8 byte alignment*/);
	if(pBin == NULL){
		printf("out of host memory\n");
		goto leaving;
	}

	if(!LoadFile(bin_name, sizeFile, pBin)){
		goto leaving;
	}

	{
		const unsigned char* listBin[1];
		listBin[0] = (unsigned char*)pBin;
		cl_int binary_status = CL_SUCCESS;
		size_t length = sizeFile;

		program = clCreateProgramWithBinary(context, (cl_uint)device_id_lists.size(), &device_id_lists[0], &length, listBin, &binary_status, &result);
	}
	if(program == NULL){
		printf("clCreateProgramWithBinary failed, %d\n", result);
		goto leaving;
	}
leaving:
	if(pBin){
		PZSDK_ALIGNED_FREE(pBin);
	}
	return program;
}


PS::S32 DensDispatchKernel(const PS::S32 tag, const PS::S32 n_walk, const EPI::Dens* epi[], const PS::S32 Nepi[], const EPJ::Dens* epj[], const PS::S32 Nepj[]){
	assert(n_walk <= N_WALK_LIMIT);
	cl_int return_code;
	PS::S32 ni_tot = 0;
	dens_host.j_disp[0] = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += Nepi[k];
		dens_host.j_disp[k + 1] = dens_host.j_disp[k] + Nepj[k];
	}
	dens_host.j_disp[n_walk + 1] = dens_host.j_disp[n_walk];
	assert(ni_tot < NI_LIMIT);
	assert(dens_host.j_disp[n_walk] < NJ_LIMIT);

	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.dens.j_disp, CL_TRUE, 0, (n_walk + 2) * sizeof(int), dens_host.j_disp, 0, NULL, NULL);
	ni_tot = 0;
	int nj_tot = 0;
	for(int iw = 0 ; iw < n_walk ; ++ iw){
		for(int i = 0 ; i < Nepi[iw] ; ++ i){
			dens_host.epi[ni_tot].rx      = epi[iw][i].pos.x;
			dens_host.epi[ni_tot].ry      = epi[iw][i].pos.y;
			dens_host.epi[ni_tot].rz      = epi[iw][i].pos.z;
			dens_host.epi[ni_tot].mass    = epi[iw][i].mass;
			dens_host.epi[ni_tot].smth    = epi[iw][i].smth;
			dens_host.epi[ni_tot].id_walk = iw;
			++ ni_tot;
		}
		for(int j = 0 ; j < Nepj[iw] ; ++ j){
			dens_host.epj[nj_tot].rx   = epj[iw][j].pos.x;
			dens_host.epj[nj_tot].ry   = epj[iw][j].pos.y;
			dens_host.epj[nj_tot].rz   = epj[iw][j].pos.z;
			dens_host.epj[nj_tot].mass = epj[iw][j].mass;
			dens_host.epj[nj_tot].smth = epj[iw][j].smth;
			++ nj_tot;
		}
	}
	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.dens.epi, CL_TRUE, 0, (ni_tot) * sizeof(Dens::EpiDev), dens_host.epi, 0, NULL, NULL);
	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.dens.epj, CL_TRUE, 0, (nj_tot) * sizeof(Dens::EpjDev), dens_host.epj, 0, NULL, NULL);

	return_code = clSetKernelArg(device.dens.kernel, 0, sizeof(cl_mem), (void*)&device.dens.j_disp);
	return_code = clSetKernelArg(device.dens.kernel, 1, sizeof(cl_mem), (void*)&device.dens.epi);
	return_code = clSetKernelArg(device.dens.kernel, 2, sizeof(cl_mem), (void*)&device.dens.epj);
	return_code = clSetKernelArg(device.dens.kernel, 3, sizeof(cl_mem), (void*)&device.dens.force);
	return_code = clSetKernelArg(device.dens.kernel, 4, sizeof(int)   , (void*)&ni_tot);

	size_t work_size = N_THREAD_MAX;
	return_code = clEnqueueNDRangeKernel(device.cmd_queue, device.dens.kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
	return 0;
}

PS::S32 DrvtDispatchKernel(const PS::S32 tag, const PS::S32 n_walk, const EPI::Drvt* epi[], const PS::S32 Nepi[], const EPJ::Drvt* epj[], const PS::S32 Nepj[]){
	assert(n_walk <= N_WALK_LIMIT);
	cl_int return_code;
	PS::S32 ni_tot = 0;
	drvt_host.j_disp[0] = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += Nepi[k];
		drvt_host.j_disp[k + 1] = drvt_host.j_disp[k] + Nepj[k];
	}
	drvt_host.j_disp[n_walk + 1] = drvt_host.j_disp[n_walk];
	assert(ni_tot < NI_LIMIT);
	assert(drvt_host.j_disp[n_walk] < NJ_LIMIT);

	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.drvt.j_disp, CL_TRUE, 0, (n_walk + 2) * sizeof(int), drvt_host.j_disp, 0, NULL, NULL);
	ni_tot = 0;
	int nj_tot = 0;
	for(int iw = 0 ; iw < n_walk ; ++ iw){
		for(int i = 0 ; i < Nepi[iw] ; ++ i){
			drvt_host.epi[ni_tot].rx      = epi[iw][i].pos.x;
			drvt_host.epi[ni_tot].ry      = epi[iw][i].pos.y;
			drvt_host.epi[ni_tot].rz      = epi[iw][i].pos.z;
			drvt_host.epi[ni_tot].vx      = epi[iw][i].vel.x;
			drvt_host.epi[ni_tot].vy      = epi[iw][i].vel.y;
			drvt_host.epi[ni_tot].vz      = epi[iw][i].vel.z;
			drvt_host.epi[ni_tot].dens    = epi[iw][i].dens;
			drvt_host.epi[ni_tot].smth    = epi[iw][i].smth;
			drvt_host.epi[ni_tot].id_walk = iw;
			++ ni_tot;
		}
		for(int j = 0 ; j < Nepj[iw] ; ++ j){
			drvt_host.epj[nj_tot].rx   = epj[iw][j].pos.x;
			drvt_host.epj[nj_tot].ry   = epj[iw][j].pos.y;
			drvt_host.epj[nj_tot].rz   = epj[iw][j].pos.z;
			drvt_host.epj[nj_tot].vx   = epj[iw][j].vel.x;
			drvt_host.epj[nj_tot].vy   = epj[iw][j].vel.y;
			drvt_host.epj[nj_tot].vz   = epj[iw][j].vel.z;
			drvt_host.epj[nj_tot].mass = epj[iw][j].mass;
			drvt_host.epj[nj_tot].smth = epj[iw][j].smth;
			++ nj_tot;
		}
	}
	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.drvt.epi, CL_TRUE, 0, (ni_tot) * sizeof(Drvt::EpiDev), drvt_host.epi, 0, NULL, NULL);
	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.drvt.epj, CL_TRUE, 0, (nj_tot) * sizeof(Drvt::EpjDev), drvt_host.epj, 0, NULL, NULL);

	return_code = clSetKernelArg(device.drvt.kernel, 0, sizeof(cl_mem), (void*)&device.drvt.j_disp);
	return_code = clSetKernelArg(device.drvt.kernel, 1, sizeof(cl_mem), (void*)&device.drvt.epi);
	return_code = clSetKernelArg(device.drvt.kernel, 2, sizeof(cl_mem), (void*)&device.drvt.epj);
	return_code = clSetKernelArg(device.drvt.kernel, 3, sizeof(cl_mem), (void*)&device.drvt.force);
	return_code = clSetKernelArg(device.drvt.kernel, 4, sizeof(int)   , (void*)&ni_tot);

	size_t work_size = N_THREAD_MAX;
	return_code = clEnqueueNDRangeKernel(device.cmd_queue, device.drvt.kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
	return 0;
}


PS::S32 HydrDispatchKernel(const PS::S32 tag, const PS::S32 n_walk, const EPI::Hydr* epi[], const PS::S32 Nepi[], const EPJ::Hydr* epj[], const PS::S32 Nepj[]){
	assert(n_walk <= N_WALK_LIMIT);
	cl_int return_code;
	PS::S32 ni_tot = 0;
	hydr_host.j_disp[0] = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += Nepi[k];
		hydr_host.j_disp[k + 1] = hydr_host.j_disp[k] + Nepj[k];
	}
	hydr_host.j_disp[n_walk + 1] = hydr_host.j_disp[n_walk];
	assert(ni_tot < NI_LIMIT);
	assert(hydr_host.j_disp[n_walk] < NJ_LIMIT);

	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.hydr.j_disp, CL_TRUE, 0, (n_walk + 2) * sizeof(int), hydr_host.j_disp, 0, NULL, NULL);
	ni_tot = 0;
	int nj_tot = 0;
	for(int iw = 0 ; iw < n_walk ; ++ iw){
		for(int i = 0 ; i < Nepi[iw] ; ++ i){
			hydr_host.epi[ni_tot].rx      = epi[iw][i].pos.x;
			hydr_host.epi[ni_tot].ry      = epi[iw][i].pos.y;
			hydr_host.epi[ni_tot].rz      = epi[iw][i].pos.z;
			hydr_host.epi[ni_tot].vx      = epi[iw][i].vel.x;
			hydr_host.epi[ni_tot].vy      = epi[iw][i].vel.y;
			hydr_host.epi[ni_tot].vz      = epi[iw][i].vel.z;
			hydr_host.epi[ni_tot].dens    = epi[iw][i].dens;
			hydr_host.epi[ni_tot].pres    = epi[iw][i].pres;
			hydr_host.epi[ni_tot].snds    = epi[iw][i].snds;
			hydr_host.epi[ni_tot].smth    = epi[iw][i].smth;
			hydr_host.epi[ni_tot].Bal     = epi[iw][i].Bal;
			hydr_host.epi[ni_tot].id_walk = iw;
			hydr_host.epi[ni_tot].grad_smth = epi[iw][i].grad_smth;
			++ ni_tot;
		}
		for(int j = 0 ; j < Nepj[iw] ; ++ j){
			hydr_host.epj[nj_tot].rx   = epj[iw][j].pos.x;
			hydr_host.epj[nj_tot].ry   = epj[iw][j].pos.y;
			hydr_host.epj[nj_tot].rz   = epj[iw][j].pos.z;
			hydr_host.epj[nj_tot].vx   = epj[iw][j].vel.x;
			hydr_host.epj[nj_tot].vy   = epj[iw][j].vel.y;
			hydr_host.epj[nj_tot].vz   = epj[iw][j].vel.z;
			hydr_host.epj[nj_tot].dens = epj[iw][j].dens;
			hydr_host.epj[nj_tot].pres = epj[iw][j].pres;
			hydr_host.epj[nj_tot].snds = epj[iw][j].snds;
			hydr_host.epj[nj_tot].mass = epj[iw][j].mass;
			hydr_host.epj[nj_tot].smth = epj[iw][j].smth;
			hydr_host.epj[nj_tot].Bal  = epj[iw][j].Bal;
			hydr_host.epj[nj_tot].grad_smth = epj[iw][j].grad_smth;
			++ nj_tot;
		}
	}
	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.hydr.epi, CL_TRUE, 0, (ni_tot) * sizeof(Hydr::EpiDev), hydr_host.epi, 0, NULL, NULL);
	return_code = clEnqueueWriteBuffer(device.cmd_queue, device.hydr.epj, CL_TRUE, 0, (nj_tot) * sizeof(Hydr::EpjDev), hydr_host.epj, 0, NULL, NULL);

	return_code = clSetKernelArg(device.hydr.kernel, 0, sizeof(cl_mem), (void*)&device.hydr.j_disp);
	return_code = clSetKernelArg(device.hydr.kernel, 1, sizeof(cl_mem), (void*)&device.hydr.epi);
	return_code = clSetKernelArg(device.hydr.kernel, 2, sizeof(cl_mem), (void*)&device.hydr.epj);
	return_code = clSetKernelArg(device.hydr.kernel, 3, sizeof(cl_mem), (void*)&device.hydr.force);
	return_code = clSetKernelArg(device.hydr.kernel, 4, sizeof(int)   , (void*)&ni_tot);

	size_t work_size = N_THREAD_MAX;
	return_code = clEnqueueNDRangeKernel(device.cmd_queue, device.hydr.kernel, 1, NULL, &work_size, NULL, 0, NULL, NULL);
	return 0;
}

PS::S32 GravDispatchKernel(const PS::S32 tag, const PS::S32 n_walk, const EPI::Grav* epi[], const PS::S32 Nepi[], const EPJ::Grav* epj[], const PS::S32 Nepj[], const PS::SPJMonopole* spj[], const PS::S32 Nspj[]){
	assert(n_walk <= N_WALK_LIMIT);
	cl_int return_code;
	PS::S32 ni_tot = 0;
	grav_host.j_disp[0] = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += Nepi[k];
		grav_host.j_disp[k + 1] = grav_host.j_disp[k] + (Nepj[k] + Nspj[k]);
	}
	grav_host.j_disp[n_walk + 1] = grav_host.j_disp[n_walk];
	assert(ni_tot < NI_LIMIT);
	assert(grav_host.j_disp[n_walk] < NJ_LIMIT);

	ni_tot = 0;
	int nj_tot = 0;
	{
		double time = get_dtime();
		for(int iw = 0 ; iw < n_walk ; ++ iw){
			for(int i = 0 ; i < Nepi[iw] ; ++ i){
				grav_host.epi[ni_tot].rx   = epi[iw][i].pos.x;
				grav_host.epi[ni_tot].ry   = epi[iw][i].pos.y;
				grav_host.epi[ni_tot].rz   = epi[iw][i].pos.z;
				grav_host.epi[ni_tot].eps2 = epi[iw][i].eps2;
				grav_host.epi[ni_tot].id_walk = iw;
				++ ni_tot;
			}
			for(int j = 0 ; j < Nepj[iw] ; ++ j){
				grav_host.epj[nj_tot].rx   = epj[iw][j].pos.x;
				grav_host.epj[nj_tot].ry   = epj[iw][j].pos.y;
				grav_host.epj[nj_tot].rz   = epj[iw][j].pos.z;
				grav_host.epj[nj_tot].mass = epj[iw][j].mass;
				++ nj_tot;
			}
			for(int j = 0 ; j < Nspj[iw] ; ++ j){
				grav_host.epj[nj_tot].rx   = spj[iw][j].pos.x;
				grav_host.epj[nj_tot].ry   = spj[iw][j].pos.y;
				grav_host.epj[nj_tot].rz   = spj[iw][j].pos.z;
				grav_host.epj[nj_tot].mass = spj[iw][j].mass;
				++ nj_tot;
			}
		}
		wtime_CP += ((double)(get_dtime() - time));
	}
	{
		double time = get_dtime();
		cl_event event[3];
		return_code = clEnqueueWriteBuffer(device.cmd_queue, device.grav.j_disp, CL_TRUE, 0, (n_walk + 2) * sizeof(int), grav_host.j_disp, 0, NULL, &event[0]);
		return_code = clEnqueueWriteBuffer(device.cmd_queue, device.grav.epi, CL_TRUE, 0, (ni_tot) * sizeof(Grav::EpiDev), grav_host.epi, 0, NULL, &event[1]);
		return_code = clEnqueueWriteBuffer(device.cmd_queue, device.grav.epj, CL_TRUE, 0, (nj_tot) * sizeof(Grav::EpjDev), grav_host.epj, 0, NULL, &event[2]);
		clWaitForEvents(1, &event[0]);
		clWaitForEvents(1, &event[1]);
		clWaitForEvents(1, &event[2]);
		wtime_WB += ((double)(get_dtime() - time));
	}
	{
		double time = get_dtime();
		return_code = clSetKernelArg(device.grav.kernel, 0, sizeof(cl_mem), (void*)&device.grav.j_disp);
		return_code = clSetKernelArg(device.grav.kernel, 1, sizeof(cl_mem), (void*)&device.grav.epi);
		return_code = clSetKernelArg(device.grav.kernel, 2, sizeof(cl_mem), (void*)&device.grav.epj);
		return_code = clSetKernelArg(device.grav.kernel, 3, sizeof(cl_mem), (void*)&device.grav.force);
		return_code = clSetKernelArg(device.grav.kernel, 4, sizeof(int)   , (void*)&ni_tot);
		wtime_SKA += (double)(get_dtime() - time);
	}
	size_t work_size = N_THREAD_MAX;
	{
		double time = get_dtime();
		cl_event event;
		return_code = clEnqueueNDRangeKernel(device.cmd_queue, device.grav.kernel, 1, NULL, &work_size, NULL, 0, NULL, &event);
		clWaitForEvents(1, &event);
		wtime_NDRK += (double)(get_dtime() - time);
	}
	return 0;
}


PS::S32 DensRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32 ni[], RESULT::Dens* force[]){
	cl_int return_code;
	int ni_tot = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += ni[k];
	}
	return_code = clEnqueueReadBuffer(device.cmd_queue, device.dens.force, CL_TRUE, 0, ni_tot * sizeof(Dens::ForceDev), dens_host.force, 0, NULL, NULL);
	int cnt = 0;
	for(int w = 0 ; w < n_walk ; ++ w){
		for(int i = 0;  i < ni[w] ; ++ i){
			force[w][i].dens = dens_host.force[cnt].dens;
			++ cnt;
		}
	}
	return 0;
}

PS::S32 DrvtRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32 ni[], RESULT::Drvt* force[]){
	cl_int return_code;
	int ni_tot = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += ni[k];
	}
	return_code = clEnqueueReadBuffer(device.cmd_queue, device.drvt.force, CL_TRUE, 0, ni_tot * sizeof(Drvt::ForceDev), drvt_host.force, 0, NULL, NULL);
	int cnt = 0;
	for(int w = 0 ; w < n_walk ; ++ w){
		for(int i = 0;  i < ni[w] ; ++ i){
			force[w][i].div_v   = drvt_host.force[cnt].div_v;
			force[w][i].rot_v.x = drvt_host.force[cnt].rot_vx;
			force[w][i].rot_v.y = drvt_host.force[cnt].rot_vy;
			force[w][i].rot_v.z = drvt_host.force[cnt].rot_vz;
			++ cnt;
		}
	}
	return 0;
}

PS::S32 HydrRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32 ni[], RESULT::Hydr* force[]){
	cl_int return_code;
	int ni_tot = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += ni[k];
	}
	return_code = clEnqueueReadBuffer(device.cmd_queue, device.hydr.force, CL_TRUE, 0, ni_tot * sizeof(Hydr::ForceDev), hydr_host.force, 0, NULL, NULL);
	int cnt = 0;
	for(int w = 0 ; w < n_walk ; ++ w){
		for(int i = 0;  i < ni[w] ; ++ i){
			force[w][i].acc.x   = hydr_host.force[cnt].ax;
			force[w][i].acc.y   = hydr_host.force[cnt].ay;
			force[w][i].acc.z   = hydr_host.force[cnt].az;
			force[w][i].eng_dot = hydr_host.force[cnt].eng_dot;
			force[w][i].dt      = hydr_host.force[cnt].dt;
			++ cnt;
		}
	}
	return 0;
}

PS::S32 GravRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32 ni[], RESULT::Grav* force[]){
	cl_int return_code;
	int ni_tot = 0;
	for(int k = 0 ; k < n_walk ; ++ k){
		ni_tot += ni[k];
	}
	{
		double time = get_dtime();
		cl_event event;
		return_code = clEnqueueReadBuffer(device.cmd_queue, device.grav.force, CL_TRUE, 0, ni_tot * sizeof(Hydr::ForceDev), grav_host.force, 0, NULL, &event);
		clWaitForEvents(1, &event);
		wtime_RB += (double)(get_dtime() - time);
	}
	int cnt = 0;
	for(int w = 0 ; w < n_walk ; ++ w){
		for(int i = 0;  i < ni[w] ; ++ i){
			force[w][i].acc.x   = grav_host.force[cnt].ax;
			force[w][i].acc.y   = grav_host.force[cnt].ay;
			force[w][i].acc.z   = grav_host.force[cnt].az;
			force[w][i].pot     = grav_host.force[cnt].pot;
			++ cnt;
		}
	}
	return 0;
}

void PezyDevice::initialize(){
	cl_int return_code;
	int rank = PS::Comm::getRank();
	return_code = clGetPlatformIDs(1, &platform_id, &num_of_platforms);
	return_code = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, NumDeviceMax, device_id, &num_of_devices);
	context     = clCreateContext(NULL, 1, &device_id[rank % NumDeviceMax], NULL, NULL, &return_code);
	cmd_queue   = clCreateCommandQueue(context, device_id[rank % NumDeviceMax], 0, &return_code);

	dens.j_disp = clCreateBuffer(context, CL_MEM_READ_WRITE, (N_WALK_LIMIT + 2) * sizeof(int),  NULL, &return_code);
	dens.epi    = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Dens::EpiDev),   NULL, &return_code);
	dens.epj    = clCreateBuffer(context, CL_MEM_READ_WRITE, NJ_LIMIT * sizeof(Dens::EpjDev),   NULL, &return_code);
	dens.force  = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Dens::ForceDev), NULL, &return_code);

	drvt.j_disp = clCreateBuffer(context, CL_MEM_READ_WRITE, (N_WALK_LIMIT + 2) * sizeof(int),  NULL, &return_code);
	drvt.epi    = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Drvt::EpiDev),   NULL, &return_code);
	drvt.epj    = clCreateBuffer(context, CL_MEM_READ_WRITE, NJ_LIMIT * sizeof(Drvt::EpjDev),   NULL, &return_code);
	drvt.force  = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Drvt::ForceDev), NULL, &return_code);

	hydr.j_disp = clCreateBuffer(context, CL_MEM_READ_WRITE, (N_WALK_LIMIT + 2) * sizeof(int),  NULL, &return_code);
	hydr.epi    = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Hydr::EpiDev),   NULL, &return_code);
	hydr.epj    = clCreateBuffer(context, CL_MEM_READ_WRITE, NJ_LIMIT * sizeof(Hydr::EpjDev),   NULL, &return_code);
	hydr.force  = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Hydr::ForceDev), NULL, &return_code);

	grav.j_disp = clCreateBuffer(context, CL_MEM_READ_WRITE, (N_WALK_LIMIT + 2) * sizeof(int),  NULL, &return_code);
	grav.epi    = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Grav::EpiDev),   NULL, &return_code);
	grav.epj    = clCreateBuffer(context, CL_MEM_READ_WRITE, NJ_LIMIT * sizeof(Grav::EpjDev),   NULL, &return_code);
	grav.force  = clCreateBuffer(context, CL_MEM_READ_WRITE, NI_LIMIT * sizeof(Grav::ForceDev), NULL, &return_code);

	device_id_list.push_back(device_id[rank % NumDeviceMax]);
	program = CreateProgram(context, device_id_list, "./kernel.sc32/kernel.pz");
	if(program == NULL){
		std::cerr << "can't create program" << std::endl;
		exit(1);
	}

	dens.kernel = clCreateKernel(program, "DensityKernel", &return_code);
	if(dens.kernel == NULL){
		std::cerr << "can't create kernel" << std::endl;
		exit(1);
	}
	return_code = clSetKernelArg(dens.kernel, 4, sizeof(int), (void*)&N_THREAD_MAX);

	drvt.kernel = clCreateKernel(program, "DerivativeKernel", &return_code);
	if(drvt.kernel == NULL){
		std::cerr << "can't create kernel" << std::endl;
		exit(1);
	}
	return_code = clSetKernelArg(drvt.kernel, 4, sizeof(int), (void*)&N_THREAD_MAX);
	
	hydr.kernel = clCreateKernel(program, "HydroKernel", &return_code);
	if(hydr.kernel == NULL){
		std::cerr << "can't create kernel" << std::endl;
		exit(1);
	}
	return_code = clSetKernelArg(hydr.kernel, 4, sizeof(int), (void*)&N_THREAD_MAX);

	grav.kernel = clCreateKernel(program, "GravityKernel", &return_code);
	if(hydr.kernel == NULL){
		std::cerr << "can't create kernel" << std::endl;
		exit(1);
	}
	return_code = clSetKernelArg(grav.kernel, 4, sizeof(int), (void*)&N_THREAD_MAX);
}


