#include <helper_cuda.h>
#include <particle_simulator.hpp>
#include "math.h"
#include "class.h"
#include "param.h"
#define CUDA_SAFE_CALL checkCudaErrors

const int N_THREAD_GPU = 64;
const int N_WALK_LIMIT = 1000;
const int NI_LIMIT     = 1000 * N_WALK_LIMIT;
const int NJ_LIMIT     = 10000 * N_WALK_LIMIT;

namespace DEVICE{
	namespace EPI{
		class Dens{
			public:
			PS::F32vec pos;
			PS::F32 mass;
			PS::F32 smth;
			PS::S32 walk_id;
			void copyFromHost(const ::EPI::Dens& epi, const PS::S32 _walk_id){
				pos = epi.pos;
				mass = epi.mass;
				smth = epi.smth;
				walk_id = _walk_id;
			}
			void inputDammyData(){
				pos.x = +1.0/0.0;
				pos.y = -1.0/0.0;
				pos.z =  0.0/0.0;
			}
		};
		class Hydro{
			public:
			PS::F32vec pos;
			PS::F32vec vel;
			PS::F32 dens;
			PS::F32 pres;
			PS::F32 snds;
			PS::F32 smth;
			PS::S32 walk_id;
			void copyFromHost(const ::EPI::Hydro& epi, const PS::S32 _walk_id){
				pos = epi.pos;
				vel = epi.vel;
				dens = epi.dens;
				pres = epi.pres;
				snds = epi.snds;
				smth = epi.smth;
				walk_id = _walk_id;
			}
			void inputDammyData(){
				pos.x = +1.0/0.0;
				pos.y = -1.0/0.0;
				pos.z =  0.0/0.0;
			}
		};
	};
	namespace EPJ{
		class Dens{
			public:
			PS::F32vec pos;
			PS::F32 mass;
			PS::F32 smth;
			void copyFromHost(const ::EPJ::Dens& epj){
				pos = epj.pos;
				mass = epj.mass;
				smth = epj.smth;
			}
		};
		class Hydro{
			public:
			PS::F32vec pos;
			PS::F32vec vel;
			PS::F32 dens;
			PS::F32 mass;
			PS::F32 smth;
			PS::F32 pres;
			PS::F32 snds;
			void copyFromHost(const ::EPJ::Hydro& epj){
				pos = epj.pos;
				vel = epj.vel;
				dens = epj.dens;
				mass = epj.mass;
				smth = epj.smth;
				pres = epj.pres;
				snds = epj.snds;
			}
		};
	};
	namespace RESULT{
		class Dens{
			public:
			PS::F32 dens;
			PS::F32 smth;
		};
		class Hydro{
			public:
			PS::F32vec acc;
			PS::F32 eng_dot;
			PS::F32 dt;
		};
	};
}

//Wendland C6
struct device_kernel_t{
	__device__ static float pi(){
		return atan(1.0) * 4.0;
	}
	__device__ static float pow8(const float x){
		float x2 = x  * x;
		float x4 = x2 * x2;
		return x4 * x4;
	}
	__device__ static float plus(const float x){
		return (x > 0) ? x : 0;
	}
	__device__ static float pow7(const float x){
		float x2 = x * x;
		float x4 = x2 * x2;
		return x4 * x2 * x;
	}
	//W
	__device__ float W(const PS::F32 r, const PS::F32 h) const{
		const float H = supportRadius() * h;
		const float s = r / H;
		float r_value;
		r_value = (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))) * pow8(plus(1.0 - s));
		r_value *= (1365./64.) / (H * H * H * pi());
		return r_value;
	}
	//gradW
	__device__ PS::F32 abs_gradW(const PS::F32 r, const PS::F32 h) const{
		const float H = supportRadius() * h;
		const float s = r / H;
		float r_value;
		r_value = pow7(plus(1.0 - s)) * (plus(1.0 - s) * (8.0 + s * (50.0 + s * (96.0))) - 8.0 * (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))));
		r_value *= (1365./64.) / (H * H * H * pi());
		return r_value / (H  + 1.0e-6 * h);
	}
	__device__ static float supportRadius(){
		return 3.5;
	}
};

namespace CalcDensity{
	//displacement
	static int *ni_displc_d, *nj_displc_d, *ni_displc_h, *nj_displc_h;
	static DEVICE::EPI::Dens *epi_d, *epi_h;
	static DEVICE::EPJ::Dens *epj_d, *epj_h;
	static DEVICE::RESULT::Dens *res_d, *res_h;
	__global__ void deviceCalcDensity(const DEVICE::EPI::Dens *epi, const int *ni_displc, const DEVICE::EPJ::Dens *epj, const int *nj_displc, DEVICE::RESULT::Dens *dens){
		const int id = blockDim.x * blockIdx.x + threadIdx.x;
		device_kernel_t kernel;
		const DEVICE::EPI::Dens& ith = epi[id];
		const int j_head = nj_displc[ith.walk_id];
		const int j_tail = nj_displc[ith.walk_id + 1];
		//const int nj = j_tail - j_head;
		float dens_buf = 0;
		for(int j = j_head ; j < j_tail ; ++ j){
			const DEVICE::EPJ::Dens& jth = epj[j];
			const PS::F32 dx = jth.pos.x - ith.pos.x;
			const PS::F32 dy = jth.pos.y - ith.pos.y;
			const PS::F32 dz = jth.pos.z - ith.pos.z;
			const PS::F32 r  = sqrt(dx * dx + dy * dy + dz * dz);
			dens_buf += jth.mass * kernel.W(r, ith.smth);
		}
		assert(PARAM::Dim == 3);
		dens[id].smth = PARAM::SMTH * cbrt(ith.mass / dens_buf);
		dens[id].dens = dens_buf;
	}

	int DispatchKernel(const PS::S32 tag, const int n_walk, const EPI::Dens** epi, const int* n_epi, const EPJ::Dens** epj, const int* n_epj){
		static bool isFirst = true;
		if(isFirst == true){
			std::cout << "Alloc Cuda Vars.." << std::endl;
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&ni_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&nj_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&ni_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&nj_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&epi_d, NI_LIMIT * sizeof(DEVICE::EPI::Dens)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&epj_d, NJ_LIMIT * sizeof(DEVICE::EPJ::Dens)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&res_d, NI_LIMIT * sizeof(DEVICE::RESULT::Dens)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&epi_h, NI_LIMIT * sizeof(DEVICE::EPI::Dens)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&epj_h, NJ_LIMIT * sizeof(DEVICE::EPJ::Dens)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&res_h, NI_LIMIT * sizeof(DEVICE::RESULT::Dens)));
			isFirst = false;
		}
		ni_displc_h[0] = nj_displc_h[0] = 0;
		for(std::size_t i = 0; i < n_walk ; ++ i){
			ni_displc_h[i+1] = ni_displc_h[i] + n_epi[i];
			nj_displc_h[i+1] = nj_displc_h[i] + n_epj[i];
		}
		const PS::S32 ni_total = ni_displc_h[n_walk];
		const int ni_total_reg = ni_displc_h[n_walk] + ((ni_total % N_THREAD_GPU != 0) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
		//make data for device on host
		int cnt = 0;
		int cnt_j = 0;
		for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
			for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
				epi_h[cnt].copyFromHost(epi[walk][i], walk);
				++ cnt;
			}
			for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
				epj_h[cnt_j].copyFromHost(epj[walk][j]);
				++ cnt_j;
			}
		}
		for(std::size_t i = cnt ; i < ni_total_reg ; ++ i){
			epi_h[i].inputDammyData();
		}

		CUDA_SAFE_CALL(cudaMemcpy(epi_d, epi_h, ni_total_reg * sizeof(DEVICE::EPI::Dens), cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(epj_d, epj_h, cnt_j * sizeof(DEVICE::EPJ::Dens), cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(ni_displc_d, ni_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(nj_displc_d, nj_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

		const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
		dim3 size_grid(n_grid, 1, 1);
		dim3 size_thread(N_THREAD_GPU, 1, 1);
		deviceCalcDensity<<<size_grid, size_thread>>> (epi_d, ni_displc_d, epj_d, nj_displc_d, res_d);
		return 0;
	}

	int RetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, RESULT::Dens** dens){
		int ni_tot = 0;
		for(int i = 0 ; i < n_walk ; ++ i){
			ni_tot += ni[i];
		}
		CUDA_SAFE_CALL(cudaMemcpy(res_h, res_d, ni_tot * sizeof(DEVICE::RESULT::Dens), cudaMemcpyDeviceToHost));
		int cnt = 0;
		for(int walk = 0 ; walk < n_walk ; ++ walk){
			for(int i = 0 ; i < ni[walk] ; ++ i){
				dens[walk][i].dens = res_h[cnt].dens;
				dens[walk][i].smth = res_h[cnt].smth;
				++ cnt;
			}
		}
		return 0;
	}
};

namespace CalcHydroForce{
	//displacement
	static int *ni_displc_d, *nj_displc_d, *ni_displc_h, *nj_displc_h;
	static DEVICE::EPI::Hydro *epi_d, *epi_h;
	static DEVICE::EPJ::Hydro *epj_d, *epj_h;
	static DEVICE::RESULT::Hydro *res_d, *res_h;

	__global__ void deviceCalcHydroForce(const DEVICE::EPI::Hydro *epi, const int *ni_displc, const DEVICE::EPJ::Hydro *epj, const int *nj_displc, DEVICE::RESULT::Hydro *force){
		const int id = blockDim.x * blockIdx.x + threadIdx.x;
		device_kernel_t kernel;
		const DEVICE::EPI::Hydro& ith = epi[id];
		const int j_head = nj_displc[ith.walk_id];
		const int j_tail = nj_displc[ith.walk_id + 1];
		float v_sig_max = 0.0;
		PS::F32 acc_x, acc_y, acc_z;
		PS::F32 eng_dot = 0;
		acc_x = acc_y = acc_z = 0;
		for(int j = j_head ; j < j_tail ; ++ j){
			const DEVICE::EPJ::Hydro& jth = epj[j];
			const PS::F32 dx = ith.pos.x - jth.pos.x;
			const PS::F32 dy = ith.pos.y - jth.pos.y;
			const PS::F32 dz = ith.pos.z - jth.pos.z;
			const PS::F32 r  = sqrt(dx * dx + dy * dy + dz * dz);
			const PS::F32 du = ith.vel.x - jth.vel.x;
			const PS::F32 dv = ith.vel.y - jth.vel.y;
			const PS::F32 dw = ith.vel.z - jth.vel.z;
			const PS::F32 xv_inner = dx * du + dy * dv + dz * dw;
			//AV
			const PS::F32 w_ij = (xv_inner < 0) ? xv_inner / r : 0;
			const PS::F32 v_sig = ith.snds + jth.snds - 3.0 * w_ij;
			v_sig_max = (v_sig_max < v_sig) ? v_sig : v_sig_max;
			const PS::F32 AV = - 0.5 * v_sig * w_ij / (0.5 * (ith.dens + jth.dens));
			//
			const PS::F32 ith_abs_gradW = kernel.abs_gradW(r, ith.smth);
			const PS::F32 jth_abs_gradW = kernel.abs_gradW(r, jth.smth);
			const PS::F32 abs_gradW = 0.5 * (ith_abs_gradW + jth_abs_gradW);
			const float4 gradW = (r > 0) ? (float4){abs_gradW * dx / r, abs_gradW * dy / r, abs_gradW * dz / r, 0.0} : (float4){0.0, 0.0, 0.0, 0.0};
			acc_x   -= jth.mass * (ith.pres / (ith.dens * ith.dens) + jth.pres / (jth.dens * jth.dens) + AV) * gradW.x;
			acc_y   -= jth.mass * (ith.pres / (ith.dens * ith.dens) + jth.pres / (jth.dens * jth.dens) + AV) * gradW.y;
			acc_z   -= jth.mass * (ith.pres / (ith.dens * ith.dens) + jth.pres / (jth.dens * jth.dens) + AV) * gradW.z;
			eng_dot += jth.mass * (ith.pres / (ith.dens * ith.dens) + 0.5 * AV) * (du * gradW.x + dv * gradW.y + dw * gradW.z);
		}
		assert(PARAM::Dim == 3);
		force[id].acc.x = acc_x;
		force[id].acc.y = acc_y;
		force[id].acc.z = acc_z;
		force[id].eng_dot = eng_dot;
		force[id].dt = PARAM::C_CFL * 2.0 * ith.smth / v_sig_max;
	}

	int DispatchKernel(const PS::S32 tag, const int n_walk, const EPI::Hydro** epi, const int* n_epi, const EPJ::Hydro** epj, const int* n_epj){
		static bool isFirst = true;
		if(isFirst == true){
			std::cout << "Alloc Cuda Vars.." << std::endl;
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&ni_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&nj_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&ni_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&nj_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&epi_d, NI_LIMIT * sizeof(DEVICE::EPI::Hydro)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&epj_d, NJ_LIMIT * sizeof(DEVICE::EPJ::Hydro)));
			CUDA_SAFE_CALL(cudaMalloc    ((void**)&res_d, NI_LIMIT * sizeof(DEVICE::RESULT::Hydro)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&epi_h, NI_LIMIT * sizeof(DEVICE::EPI::Hydro)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&epj_h, NJ_LIMIT * sizeof(DEVICE::EPJ::Hydro)));
			CUDA_SAFE_CALL(cudaMallocHost((void**)&res_h, NI_LIMIT * sizeof(DEVICE::RESULT::Hydro)));
			isFirst = false;
		}
		ni_displc_h[0] = nj_displc_h[0] = 0;
		for(std::size_t i = 0; i < n_walk ; ++ i){
			ni_displc_h[i+1] = ni_displc_h[i] + n_epi[i];
			nj_displc_h[i+1] = nj_displc_h[i] + n_epj[i];
		}
		const PS::S32 ni_total = ni_displc_h[n_walk];
		const int ni_total_reg = ni_displc_h[n_walk] + ((ni_total % N_THREAD_GPU != 0) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
		//make data for device on host
		int cnt = 0;
		int cnt_j = 0;
		for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
			for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
				epi_h[cnt].copyFromHost(epi[walk][i], walk);
				++ cnt;
			}
			for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
				epj_h[cnt_j].copyFromHost(epj[walk][j]);
				++ cnt_j;
			}
		}
		for(std::size_t i = cnt ; i < ni_total_reg ; ++ i){
			epi_h[i].inputDammyData();
		}

		CUDA_SAFE_CALL(cudaMemcpy(epi_d, epi_h, ni_total_reg * sizeof(DEVICE::EPI::Hydro), cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(epj_d, epj_h, cnt_j * sizeof(DEVICE::EPJ::Hydro), cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(ni_displc_d, ni_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(nj_displc_d, nj_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

		const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
		dim3 size_grid(n_grid, 1, 1);
		dim3 size_thread(N_THREAD_GPU, 1, 1);
		deviceCalcHydroForce<<<size_grid, size_thread>>> (epi_d, ni_displc_d, epj_d, nj_displc_d, res_d);
		return 0;
	}
	int RetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, RESULT::Hydro** force){
		int ni_tot = 0;
		for(int i = 0 ; i < n_walk ; ++ i){
			ni_tot += ni[i];
		}
		CUDA_SAFE_CALL(cudaMemcpy(res_h, res_d, ni_tot * sizeof(DEVICE::RESULT::Hydro), cudaMemcpyDeviceToHost));
		int cnt = 0;
		for(int walk = 0 ; walk < n_walk ; ++ walk){
			for(int i = 0 ; i < ni[walk] ; ++ i){
				force[walk][i].acc = res_h[cnt].acc;
				force[walk][i].eng_dot = res_h[cnt].eng_dot;
				force[walk][i].dt = res_h[cnt].dt;
				++ cnt;
			}
		}
		return 0;
	}

};

