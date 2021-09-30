#include <particle_simulator.hpp>
#include "./pzc/class_device.hpp"
#include "class_platform.hpp"
#include "kernel.h"

struct kernel_t_device{
	__device__ static real pow8(const real x){
		real x2 = x  * x;
		real x4 = x2 * x2;
		return x4 * x4;
	}
	__device__ static real plus(const real x){
		return (x > 0.0f) ? x : 0.0f;
	}
	__device__ static real pow7(const real x){
		real x2 = x * x;
		real x4 = x2 * x2;
		return x4 * x2 * x;
	}
	//W
	__device__ real W(const real r, const real h) const{
		const real H = supportRadius() * h;
		const real s = r / H;
		real r_value;
		r_value = (1.0f + s * (8.0f + s * (25.0f + s * (32.0f)))) * pow8(plus(1.0f - s));
		r_value *= (1365.f/64.f) / (H * H * H * M_PI);
		return r_value;
	}
	//gradW
	__device__ real gradW(const real r, const real h) const{
		const real H = supportRadius() * h;
		const real s = r / H;
		real r_value;
		r_value = pow7(plus(1.0f - s)) * (plus(1.0f - s) * (8.0f + s * (50.0f + s * (96.0f))) - 8.0f * (1.0f + s * (8.0f + s * (25.0f + s * (32.0f)))));
		r_value *= (1365.f/64.f) / (H * H * H * M_PI);
		return r_value / (H + 0.01f * h);
	}
	__device__ static real supportRadius(){
		return 3.5f;
	}
};

struct dens_host_t{
	int *ni_displc_d, *nj_displc_d, *ni_displc_h, *nj_displc_h;
	Dens::EpiDev *epi_d, *epi_h;
	Dens::EpjDev *epj_d, *epj_h;
	Dens::ForceDev *res_d, *res_h;
}dens_host;

struct drvt_host_t{
	int *ni_displc_d, *nj_displc_d, *ni_displc_h, *nj_displc_h;
	Drvt::EpiDev *epi_d, *epi_h;
	Drvt::EpjDev *epj_d, *epj_h;
	Drvt::ForceDev *res_d, *res_h;
}drvt_host;

struct hydr_host_t{
	int *ni_displc_d, *nj_displc_d, *ni_displc_h, *nj_displc_h;
	Hydr::EpiDev *epi_d, *epi_h;
	Hydr::EpjDev *epj_d, *epj_h;
	Hydr::ForceDev *res_d, *res_h;
}hydr_host;


__global__ void deviceCalcDensity(const Dens::EpiDev *epi, const int *ni_displc, const Dens::EpjDev *epj, const int *nj_displc, Dens::ForceDev *dens){
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	kernel_t_device kernel;
	const Dens::EpiDev& ith = epi[id];
	const int j_head = nj_displc[ith.id_walk];
	const int j_tail = nj_displc[ith.id_walk + 1];
	real dens_buf = 0.0f;
	for(int j = j_head ; j < j_tail ; ++ j){
		const Dens::EpjDev& jth = epj[j];
		float3 dr = make_float3(jth.rx - ith.rx, jth.ry - ith.ry, jth.rz - ith.rz);
		const real r = sqrtf(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
		dens_buf += jth.mass * kernel.W(r, ith.smth);
	}
	dens[id].dens = dens_buf;
}

__global__ void deviceCalcDerivative(const Drvt::EpiDev *epi, const int *ni_displc, const Drvt::EpjDev *epj, const int *nj_displc, Drvt::ForceDev *force){
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	kernel_t_device kernel;
	const Drvt::EpiDev& ith = epi[id];
	const int j_head = nj_displc[ith.id_walk];
	const int j_tail = nj_displc[ith.id_walk + 1];
	float4 force_buf = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	for(int j = j_head ; j < j_tail ; ++ j){
		const Drvt::EpjDev& jth = epj[j];
		float3 dr = make_float3(jth.rx - ith.rx, jth.ry - ith.ry, jth.rz - ith.rz);
		float3 dv = make_float3(jth.vx - ith.vx, jth.vy - ith.vy, jth.vz - ith.vz);
		const real r = sqrtf(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) + 1.0e-4f;
		const real rinv = 1.0f / r;
		const real drdv = dr.x * dv.x + dr.y * dv.y + dr.z * dv.z;

		const real jth_mass_ith_abs_gradW = jth.mass * kernel.gradW(r, ith.smth) * rinv;
		force_buf.w += - drdv * jth_mass_ith_abs_gradW;
		force_buf.x += - (dr.y * dv.z - dr.z * dv.y) * jth_mass_ith_abs_gradW;
		force_buf.y += - (dr.z * dv.x - dr.x * dv.z) * jth_mass_ith_abs_gradW;
		force_buf.z += - (dr.x * dv.y - dr.y * dv.x) * jth_mass_ith_abs_gradW;
	}
	force[id].rot_vx = force_buf.x;
	force[id].rot_vy = force_buf.y;
	force[id].rot_vz = force_buf.z;
	force[id].div_v  = force_buf.w;
}

__global__ void deviceCalcHydroForce(const Hydr::EpiDev *epi, const int *ni_displc, const Hydr::EpjDev *epj, const int *nj_displc, Hydr::ForceDev *force){
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	kernel_t_device kernel;
	const Hydr::EpiDev& ith = epi[id];
	const int j_head = nj_displc[ith.id_walk];
	const int j_tail = nj_displc[ith.id_walk + 1];


	float v_sig_max = 0.0f;
	float4 force_buf = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	const real ith_pres_over_dens2 = ith.pres / (ith.dens * ith.dens);
	for(int j = j_head ; j < j_tail ; ++ j){
		const Hydr::EpjDev& jth = epj[j];
		float3 dr = make_float3(jth.rx - ith.rx, jth.ry - ith.ry, jth.rz - ith.rz);
		float3 dv = make_float3(jth.vx - ith.vx, jth.vy - ith.vy, jth.vz - ith.vz);
		const real r  = sqrtf(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) + 1.0e-4f;
		const real rinv = 1.0f / r;
		const real drdv = dr.x * dv.x + dr.y * dv.y + dr.z * dv.z;
		//AV
		const real w_ij = (drdv < 0.0f) ? drdv * rinv : 0.0f;
		const real v_sig = ith.snds + jth.snds - 3.0f * w_ij;
		v_sig_max = (v_sig_max < v_sig) ? v_sig : v_sig_max;
		const real AV = - 0.5f * v_sig * w_ij / (0.5f * (ith.dens + jth.dens));
		const real ith_abs_gradW = kernel.gradW(r, ith.smth);
		const real jth_abs_gradW = kernel.gradW(r, jth.smth);
		const real abs_gradW_rinv = 0.5f * (ith_abs_gradW + jth_abs_gradW) * rinv;
		const float3 gradW = make_float3(abs_gradW_rinv * dr.x, abs_gradW_rinv * dr.y, abs_gradW_rinv * dr.z);
		const real acc = jth.mass * (ith_pres_over_dens2 + jth.pres / (jth.dens * jth.dens) + AV);
		force_buf.x -= acc * gradW.x;
		force_buf.y -= acc * gradW.y;
		force_buf.z -= acc * gradW.z;
		force_buf.w += jth.mass * (ith_pres_over_dens2 + 0.5f * AV) * (dv.x * gradW.x + dv.y * gradW.y + dv.z * gradW.z);
	}
	force[id].ax = force_buf.x;
	force[id].ay = force_buf.y;
	force[id].az = force_buf.z;
	force[id].eng_dot = force_buf.w;
	force[id].dt = PARAM::C_CFL * 2.0f * ith.smth / v_sig_max;
}

int DensDispatchKernel(const PS::S32 tag, const int n_walk, const EPI::Dens* epi[], const int* n_epi, const EPJ::Dens* epj[], const int* n_epj){
	static bool isFirst = true;
	if(isFirst == true){
		std::cout << "Alloc Cuda Vars.." << std::endl;
		(cudaMalloc    ((void**)&dens_host.ni_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&dens_host.nj_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&dens_host.ni_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&dens_host.nj_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&dens_host.epi_d, NI_LIMIT * sizeof(Dens::EpiDev)));
		(cudaMalloc    ((void**)&dens_host.epj_d, NJ_LIMIT * sizeof(Dens::EpjDev)));
		(cudaMalloc    ((void**)&dens_host.res_d, NI_LIMIT * sizeof(Dens::ForceDev)));
		(cudaMallocHost((void**)&dens_host.epi_h, NI_LIMIT * sizeof(Dens::EpiDev)));
		(cudaMallocHost((void**)&dens_host.epj_h, NJ_LIMIT * sizeof(Dens::EpjDev)));
		(cudaMallocHost((void**)&dens_host.res_h, NI_LIMIT * sizeof(Dens::ForceDev)));
		isFirst = false;
	}
	dens_host.ni_displc_h[0] = dens_host.nj_displc_h[0] = 0;
	for(std::size_t i = 0; i < n_walk ; ++ i){
		dens_host.ni_displc_h[i+1] = dens_host.ni_displc_h[i] + n_epi[i];
		dens_host.nj_displc_h[i+1] = dens_host.nj_displc_h[i] + n_epj[i];
	}
	const PS::S32 ni_total = dens_host.ni_displc_h[n_walk];
	if(ni_total >= NI_LIMIT){
		std::cout << ni_total << " >= " << NI_LIMIT << std::endl;
		assert(ni_total < NI_LIMIT);
	}
	const int ni_total_reg = dens_host.ni_displc_h[n_walk] + ((ni_total % N_THREAD_GPU != 1) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
	//make data for device on host
	int cnt = 0;
	int cnt_j = 0;
	for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
		for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
			dens_host.epi_h[cnt].rx = epi[walk][i].pos.x;
			dens_host.epi_h[cnt].ry = epi[walk][i].pos.y;
			dens_host.epi_h[cnt].rz = epi[walk][i].pos.z;
			dens_host.epi_h[cnt].mass = epi[walk][i].mass;
			dens_host.epi_h[cnt].smth = epi[walk][i].smth;
			dens_host.epi_h[cnt].id_walk = walk;
			++ cnt;
		}
		for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
			dens_host.epj_h[cnt_j].rx = epj[walk][j].pos.x;
			dens_host.epj_h[cnt_j].ry = epj[walk][j].pos.y;
			dens_host.epj_h[cnt_j].rz = epj[walk][j].pos.z;
			dens_host.epj_h[cnt_j].mass = epj[walk][j].mass;
			dens_host.epj_h[cnt_j].smth = epj[walk][j].smth;
			++ cnt_j;
		}
	}

	(cudaMemcpy(dens_host.epi_d, dens_host.epi_h, ni_total_reg * sizeof(Dens::EpiDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(dens_host.epj_d, dens_host.epj_h, cnt_j * sizeof(Dens::EpjDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(dens_host.ni_displc_d, dens_host.ni_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
	(cudaMemcpy(dens_host.nj_displc_d, dens_host.nj_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

	const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
	dim3 size_grid(n_grid, 1, 1);
	dim3 size_thread(N_THREAD_GPU, 1, 1);
	deviceCalcDensity<<<size_grid, size_thread>>> (dens_host.epi_d, dens_host.ni_displc_d, dens_host.epj_d, dens_host.nj_displc_d, dens_host.res_d);
	return 0;
}

int DensRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, RESULT::Dens* force[]){
	int ni_tot = 0;
	for(int i = 0 ; i < n_walk ; ++ i){
		ni_tot += ni[i];
	}
	(cudaMemcpy(dens_host.res_h, dens_host.res_d, ni_tot * sizeof(Dens::ForceDev), cudaMemcpyDeviceToHost));
	int cnt = 0;
	for(int walk = 0 ; walk < n_walk ; ++ walk){
		for(int i = 0 ; i < ni[walk] ; ++ i){
			force[walk][i].dens = dens_host.res_h[cnt].dens;
			++ cnt;
		}
	}
	return 0;
}

int DrvtDispatchKernel(const PS::S32 tag, const int n_walk, const EPI::Drvt* epi[], const int* n_epi, const EPJ::Drvt* epj[], const int* n_epj){
	static bool isFirst = true;
	if(isFirst == true){
		std::cout << "Alloc Cuda Vars.." << std::endl;
		(cudaMalloc    ((void**)&drvt_host.ni_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&drvt_host.nj_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&drvt_host.ni_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&drvt_host.nj_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&drvt_host.epi_d, NI_LIMIT * sizeof(Drvt::EpiDev)));
		(cudaMalloc    ((void**)&drvt_host.epj_d, NJ_LIMIT * sizeof(Drvt::EpjDev)));
		(cudaMalloc    ((void**)&drvt_host.res_d, NI_LIMIT * sizeof(Drvt::ForceDev)));
		(cudaMallocHost((void**)&drvt_host.epi_h, NI_LIMIT * sizeof(Drvt::EpiDev)));
		(cudaMallocHost((void**)&drvt_host.epj_h, NJ_LIMIT * sizeof(Drvt::EpjDev)));
		(cudaMallocHost((void**)&drvt_host.res_h, NI_LIMIT * sizeof(Drvt::ForceDev)));
		isFirst = false;
	}
	drvt_host.ni_displc_h[0] = drvt_host.nj_displc_h[0] = 0;
	for(std::size_t i = 0; i < n_walk ; ++ i){
		drvt_host.ni_displc_h[i+1] = drvt_host.ni_displc_h[i] + n_epi[i];
		drvt_host.nj_displc_h[i+1] = drvt_host.nj_displc_h[i] + n_epj[i];
	}
	const PS::S32 ni_total = drvt_host.ni_displc_h[n_walk];
	if(ni_total >= NI_LIMIT){
		std::cout << ni_total << " >= " << NI_LIMIT << std::endl;
		assert(ni_total < NI_LIMIT);
	}
	if(drvt_host.nj_displc_h[n_walk] >= NJ_LIMIT){
		std::cout << drvt_host.nj_displc_h[n_walk] << " >= " << NJ_LIMIT << std::endl;
		assert(drvt_host.nj_displc_h[n_walk] < NJ_LIMIT);
	}

	const int ni_total_reg = drvt_host.ni_displc_h[n_walk] + ((ni_total % N_THREAD_GPU != 1) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
	//make data for device on host
	int cnt = 0;
	int cnt_j = 0;
	for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
		for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
			drvt_host.epi_h[cnt].rx = epi[walk][i].pos.x;
			drvt_host.epi_h[cnt].ry = epi[walk][i].pos.y;
			drvt_host.epi_h[cnt].rz = epi[walk][i].pos.z;
			drvt_host.epi_h[cnt].vx = epi[walk][i].vel.x;
			drvt_host.epi_h[cnt].vy = epi[walk][i].vel.y;
			drvt_host.epi_h[cnt].vz = epi[walk][i].vel.z;
			drvt_host.epi_h[cnt].dens = epi[walk][i].dens;
			drvt_host.epi_h[cnt].smth = epi[walk][i].smth;
			drvt_host.epi_h[cnt].id_walk = walk;
			++ cnt;
		}
		for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
			drvt_host.epj_h[cnt_j].rx = epj[walk][j].pos.x;
			drvt_host.epj_h[cnt_j].ry = epj[walk][j].pos.y;
			drvt_host.epj_h[cnt_j].rz = epj[walk][j].pos.z;
			drvt_host.epj_h[cnt_j].vx = epj[walk][j].vel.x;
			drvt_host.epj_h[cnt_j].vy = epj[walk][j].vel.y;
			drvt_host.epj_h[cnt_j].vz = epj[walk][j].vel.z;
			drvt_host.epj_h[cnt_j].mass = epj[walk][j].mass;
			drvt_host.epj_h[cnt_j].smth = epj[walk][j].smth;
			++ cnt_j;
		}
	}

	(cudaMemcpy(drvt_host.epi_d, drvt_host.epi_h, ni_total_reg * sizeof(Drvt::EpiDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(drvt_host.epj_d, drvt_host.epj_h, cnt_j * sizeof(Drvt::EpjDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(drvt_host.ni_displc_d, drvt_host.ni_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
	(cudaMemcpy(drvt_host.nj_displc_d, drvt_host.nj_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

	const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
	dim3 size_grid(n_grid, 1, 1);
	dim3 size_thread(N_THREAD_GPU, 1, 1);
	deviceCalcDerivative<<<size_grid, size_thread>>> (drvt_host.epi_d, drvt_host.ni_displc_d, drvt_host.epj_d, drvt_host.nj_displc_d, drvt_host.res_d);
	return 0;
}

int DrvtRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, RESULT::Drvt* force[]){
	int ni_tot = 0;
	for(int i = 0 ; i < n_walk ; ++ i){
		ni_tot += ni[i];
	}
	(cudaMemcpy(drvt_host.res_h, drvt_host.res_d, ni_tot * sizeof(Drvt::ForceDev), cudaMemcpyDeviceToHost));
	int cnt = 0;
	for(int walk = 0 ; walk < n_walk ; ++ walk){
		for(int i = 0 ; i < ni[walk] ; ++ i){
			force[walk][i].div_v   = drvt_host.res_h[cnt].div_v;
			force[walk][i].rot_v.x = drvt_host.res_h[cnt].rot_vx;
			force[walk][i].rot_v.y = drvt_host.res_h[cnt].rot_vy;
			force[walk][i].rot_v.z = drvt_host.res_h[cnt].rot_vz;
			++ cnt;
		}
	}
	return 0;
}

int HydrDispatchKernel(const PS::S32 tag, const int n_walk, const EPI::Hydr** epi, const int* n_epi, const EPJ::Hydr** epj, const int* n_epj){
	static bool isFirst = true;
	if(isFirst == true){
		std::cout << "Alloc Cuda Vars.." << std::endl;
		(cudaMalloc    ((void**)&hydr_host.ni_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&hydr_host.nj_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&hydr_host.ni_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&hydr_host.nj_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&hydr_host.epi_d, NI_LIMIT * sizeof(Hydr::EpiDev)));
		(cudaMalloc    ((void**)&hydr_host.epj_d, NJ_LIMIT * sizeof(Hydr::EpjDev)));
		(cudaMalloc    ((void**)&hydr_host.res_d, NI_LIMIT * sizeof(Hydr::ForceDev)));
		(cudaMallocHost((void**)&hydr_host.epi_h, NI_LIMIT * sizeof(Hydr::EpiDev)));
		(cudaMallocHost((void**)&hydr_host.epj_h, NJ_LIMIT * sizeof(Hydr::EpjDev)));
		(cudaMallocHost((void**)&hydr_host.res_h, NI_LIMIT * sizeof(Hydr::ForceDev)));
		isFirst = false;
	}
	hydr_host.ni_displc_h[0] = hydr_host.nj_displc_h[0] = 0;
	for(std::size_t i = 0; i < n_walk ; ++ i){
		hydr_host.ni_displc_h[i+1] = hydr_host.ni_displc_h[i] + n_epi[i];
		hydr_host.nj_displc_h[i+1] = hydr_host.nj_displc_h[i] + n_epj[i];
	}
	const PS::S32 ni_total = hydr_host.ni_displc_h[n_walk];
	const int ni_total_reg = hydr_host.ni_displc_h[n_walk] + ((ni_total % N_THREAD_GPU != 0) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
	//make data for device on host
	int cnt = 0;
	int cnt_j = 0;
	for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
		for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
			hydr_host.epi_h[cnt].rx      = epi[walk][i].pos.x;
			hydr_host.epi_h[cnt].ry      = epi[walk][i].pos.y;
			hydr_host.epi_h[cnt].rz      = epi[walk][i].pos.z;
			hydr_host.epi_h[cnt].vx      = epi[walk][i].vel.x;
			hydr_host.epi_h[cnt].vy      = epi[walk][i].vel.y;
			hydr_host.epi_h[cnt].vz      = epi[walk][i].vel.z;
			hydr_host.epi_h[cnt].dens    = epi[walk][i].dens;
			hydr_host.epi_h[cnt].pres    = epi[walk][i].pres;
			hydr_host.epi_h[cnt].snds    = epi[walk][i].snds;
			hydr_host.epi_h[cnt].smth    = epi[walk][i].smth;
			hydr_host.epi_h[cnt].Bal     = epi[walk][i].Bal;
			hydr_host.epi_h[cnt].id_walk = walk;
			hydr_host.epi_h[cnt].grad_smth = epi[walk][i].grad_smth;
			++ cnt;
		}
		for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
			hydr_host.epj_h[cnt_j].rx   = epj[walk][j].pos.x;
			hydr_host.epj_h[cnt_j].ry   = epj[walk][j].pos.y;
			hydr_host.epj_h[cnt_j].rz   = epj[walk][j].pos.z;
			hydr_host.epj_h[cnt_j].vx   = epj[walk][j].vel.x;
			hydr_host.epj_h[cnt_j].vy   = epj[walk][j].vel.y;
			hydr_host.epj_h[cnt_j].vz   = epj[walk][j].vel.z;
			hydr_host.epj_h[cnt_j].dens = epj[walk][j].dens;
			hydr_host.epj_h[cnt_j].pres = epj[walk][j].pres;
			hydr_host.epj_h[cnt_j].snds = epj[walk][j].snds;
			hydr_host.epj_h[cnt_j].mass = epj[walk][j].mass;
			hydr_host.epj_h[cnt_j].smth = epj[walk][j].smth;
			hydr_host.epj_h[cnt_j].Bal  = epj[walk][j].Bal;
			hydr_host.epj_h[cnt_j].grad_smth = epj[walk][j].grad_smth;
			++ cnt_j;
		}
	}

	(cudaMemcpy(hydr_host.epi_d, hydr_host.epi_h, ni_total_reg * sizeof(Hydr::EpiDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(hydr_host.epj_d, hydr_host.epj_h, cnt_j * sizeof(Hydr::EpjDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(hydr_host.ni_displc_d, hydr_host.ni_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
	(cudaMemcpy(hydr_host.nj_displc_d, hydr_host.nj_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

	const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
	dim3 size_grid(n_grid, 1, 1);
	dim3 size_thread(N_THREAD_GPU, 1, 1);
	deviceCalcHydroForce<<<size_grid, size_thread>>> (hydr_host.epi_d, hydr_host.ni_displc_d, hydr_host.epj_d, hydr_host.nj_displc_d, hydr_host.res_d);
	return 0;
}

int HydrRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, RESULT::Hydr** force){
	int ni_tot = 0;
	for(int i = 0 ; i < n_walk ; ++ i){
		ni_tot += ni[i];
	}
	(cudaMemcpy(hydr_host.res_h, hydr_host.res_d, ni_tot * sizeof(Hydr::ForceDev), cudaMemcpyDeviceToHost));
	int cnt = 0;
	for(int walk = 0 ; walk < n_walk ; ++ walk){
		for(int i = 0 ; i < ni[walk] ; ++ i){
			force[walk][i].acc.x = hydr_host.res_h[cnt].ax;
			force[walk][i].acc.y = hydr_host.res_h[cnt].ay;
			force[walk][i].acc.z = hydr_host.res_h[cnt].az;
			force[walk][i].eng_dot = hydr_host.res_h[cnt].eng_dot;
			force[walk][i].dt = hydr_host.res_h[cnt].dt;
			++ cnt;
		}
	}
	return 0;
}

