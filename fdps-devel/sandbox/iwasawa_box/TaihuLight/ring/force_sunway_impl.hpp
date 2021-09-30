#pragma once

#include<particle_simulator.hpp>
#include "user-defined.hpp"

extern PS::F64 MASS_DUST_TOTAL;
extern PS::F64vec CM_POS_DUST_TOTAL;

extern MPI_Comm MY_MPI_COMM_SUB;
extern MPI_Comm MY_MPI_COMM_1D;
extern PS::S32 MY_RANK_SUB;
extern PS::S32 N_LOOP_GLB;

#ifdef SUNWAY
extern "C"{
   #include <athread.h>
   #include "cpe_func.h"
   extern void SLAVE_FUN(ForceKernelSunWay1st)(void*);
   extern void SLAVE_FUN(ForceKernelSunWay2nd)(void*);
}
Force_Kernel_Pars cpe_pars; // Force_Kernel_Pars is defined in cpe_func.h
#endif
#ifdef MPE_FORCE_KERNEL
#include "mpe_func.h"
MPE_pars mpe_pars;
#endif

PS::F64 wtime_copy_all_epj   = 0.0;
PS::F64 wtime_copy_all_spj   = 0.0;
PS::F64 wtime_calc_pj    = 0.0;
PS::F64 wtime_calc_epj   = 0.0;
PS::F64 wtime_calc_spj   = 0.0;
PS::F64 wtime_copy_ptcl_to_mm = 0.0;
PS::F64 wtime_copy_force_from_mm = 0.0;
PS::F64 wtime_kick = 0.0;
PS::F64 wtime_drift = 0.0;
PS::F64 wtime_calc_energy = 0.0;

PS::F64 ENERGY_TOT, ENERGY_KIN, ENERGY_POT;
bool DISPATCH_MODE;


extern PS::F64 DT_TREE;
extern Planet PLANET;
extern PS::S32 N_SATELLITE;
extern PS::ReallocatableArray<Satellite> SATELLITE;
extern PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE;
extern PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE_LOC;
extern std::ofstream fout_wtime;
extern std::ofstream fout_debug;

TreeType * TREE_POINTER;

#if !defined(SUNWAY) || !defined(SUNWAY_FORCE_KERNEL)
const PS::S64 N_EPJ_LIMIT = 20000000;
const PS::S64 N_SPJ_LIMIT = 10000000;

//Epi EPI_ARRAY[N_EPI_LIMIT];
Epj EPJ_ARRAY[N_EPJ_LIMIT];
PS::SPJMonopoleScatter SPJ_ARRAY[N_SPJ_LIMIT];
#endif
const PS::S64 N_WALK_LIMIT = 500;
void * EPI_POINTER;
#if 0
const PS::S64 N_EPI_LIMIT = 2000000;
Force FORCE_ARRAY[N_EPI_LIMIT];
#endif

struct LapTimer{ 
	enum{
		STEP_PER_LOOP = 64,
	};
	enum{ // sokutei kukan
		BEGIN_OF_STEP = 0,
		BEFORE_REDUCE,
		AFTER_REDUCE,
		END_OF_STEP,

		NUM_MEASURE, // = 4
	};

	// This must be called after MPI_Init()
	static void initialize(){
		init_time() = MPI_Wtime();
	}

	static int &step_counter(){
		static int cnt;
		return cnt;
	}

	static void measure(const int step, const int place){
		const double t = MPI_Wtime() - init_time();
		at(step, place) = t;
	}

	static void summary(FILE *fp=stdout){
		struct double_int{
			double d;
			int    i;
		};
		double_int max_time[STEP_PER_LOOP][NUM_MEASURE]; // 1 KiB from stack
		double_int min_time[STEP_PER_LOOP][NUM_MEASURE];
		double     avr_time[STEP_PER_LOOP][NUM_MEASURE];
		
		double &src = at(0,0);
		const int nwords = STEP_PER_LOOP * NUM_MEASURE;
		const int rank = MPI::COMM_WORLD.Get_rank();
		const int size = MPI::COMM_WORLD.Get_size();

		const int root = 0;

		double_int src2[STEP_PER_LOOP][NUM_MEASURE];
		for(int i=0; i<STEP_PER_LOOP; i++){
			for(int j=0; j<NUM_MEASURE; j++){
				src2[i][j].d = at(i, j);
				src2[i][j].i = rank;
			}
		}
		
		// int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
		//MPI_Reduce(src2, max_time, nwords, MPI_DOUBLE_INT, MPI_MAXLOC, root, MPI_COMM_WORLD);
		//MPI_Reduce(src2, min_time, nwords, MPI_DOUBLE_INT, MPI_MINLOC, root, MPI_COMM_WORLD);
		//MPI_Reduce(&src, avr_time, nwords, MPI_DOUBLE,     MPI_SUM,    root, MPI_COMM_WORLD);
		MPI_Allreduce(src2, max_time, nwords, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
		MPI_Allreduce(src2, min_time, nwords, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
		MPI_Allreduce(&src, avr_time, nwords, MPI_DOUBLE,     MPI_SUM,    MPI_COMM_WORLD);

		//for(int i=0; i<nwords; i++){
		//	avr_time[0][i] *= 1.0/size;
		//}
		for(int i=0; i<STEP_PER_LOOP; i++){
			for(int j=0; j<NUM_MEASURE; j++){
				avr_time[i][j] *= 1.0/size;
			}
		}

		if(rank == root){
			static const char *names[] = {
				"begin_of_step :",
				"before_reduce :",
				"after_reduce  :",
				"end_of_step   :",
			};
			for(int step=0; step<STEP_PER_LOOP; step++){
				printf("step %d:\n", step);
				for(int m=0; m<NUM_MEASURE; m++){
                                     fprintf(fp, "%s max=%24.16e (rank %d), min=%24.16e (rank %d), avr=%24.16e, max-min=%24.16e\n",
                                             names[m],
                                             max_time[step][m].d, max_time[step][m].i,
                                             min_time[step][m].d, min_time[step][m].i,
                                             avr_time[step][m],
                                             max_time[step][m].d - min_time[step][m].d);
                                             //0.0);
				}
				fprintf(fp, "\n");
			}
			fflush(fp);
		}
	}

private:
	static double &init_time(){
		static double t;
		return t;
	}

	static double &at(const int step, const int place){
		static double time[STEP_PER_LOOP][NUM_MEASURE];
		return time[step][place];
	}
};

#if 0 // SAMPLE USAGE

static void randam_wait(){
	int n = lrand48() >> 16;
	for(int i=0; i<n; i++){
		asm volatile ("nop");
	}
}


int main(int ac, char **av){
	srand48(20170415);

	randam_wait();

	MPI_Init(&ac, &av);

	LapTimer::initialize();

	for(int s=0; s<LapTimer::STEP_PER_LOOP; s++){
		LapTimer::measure(s, LapTimer::BEGIN_OF_STEP);
		randam_wait();
		LapTimer::measure(s, LapTimer::BEFORE_REDUCE);
		MPI_Barrier(MPI_COMM_WORLD);
		LapTimer::measure(s, LapTimer::AFTER_REDUCE);
		randam_wait();
		LapTimer::measure(s, LapTimer::END_OF_STEP);
	}
	LapTimer::summary();

	MPI_Finalize();

	return 0;
}
#endif

template<class Tptcl, class Tptclforce, class Tsat, class Tsatforce>
void CalcEnergy(const Tptcl ptcl[],
                const Tptclforce ptcl_force[],
                const PS::S32 n_ptcl,
                const Tsat satellite[],
                const Tsatforce satellite_force[],
                const PS::S32 n_sat,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot){
    etot = ekin = epot = 0.0;
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    for(PS::S32 i=0; i<n_ptcl; i++){
        ekin_loc += ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
        epot_loc += ptcl[i].mass * ptcl_force[i].pot;
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
    //std::cerr<<"ekin= "<<ekin<<" epot= "<<epot<<" etot= "<<etot<<std::endl;
    for(PS::S32 i=0; i<n_sat; i++){
        ekin += 0.5 * satellite[i].mass * satellite[i].vel * satellite[i].vel;
        epot += 0.5 * satellite[i].mass * satellite_force[i].pot;
    }
    etot = ekin + epot;
}

template<class Tpi, class Tpj, class Tforce>
void CalcGravPair(const Tpi & pi,
                  const Tpj & pj,
                  Tforce & force,
                  const PS::F64 eps_sq){
    const PS::F64vec rij = pi.pos - pj.pos;
    const PS::F64 mj = pj.mass;
    const PS::F64 r_sq = rij * rij + eps_sq;
    const PS::F64 r_inv = 1.0 / sqrt(r_sq);
    const PS::F64 pij = mj * r_inv;
    const double mri3 = r_inv*r_inv*pij;
    force.acc -= mri3 * rij;
    force.pot -= pij;
}

template<class Tpi, class Tpj, class Tforce>
void CalcGravEx(const Tpi & pi,
                const Tpj & pj,
                Tforce & force,
                const PS::F64 eps_sq){
    const PS::F64vec rij = pi.pos - pj.pos;
    const PS::F64 mj = pj.mass;
    const PS::F64 r_sq = rij * rij + eps_sq;
    const PS::F64 r_inv = 1.0 / sqrt(r_sq);
    const PS::F64 pij = mj * r_inv;
    const double mri3 = r_inv*r_inv*pij;
    force.acc -= mri3 * rij;
    force.pot -= 2.0*pij; // factor 2 is just a trick
}

template<class Tpi, class Tpj, class Tforce>
void CalcGravFormArray(const Tpi & pi,
                       const Tpj & pj,
                       Tforce & force,
                       const PS::F64 eps_sq){
    const PS::S32 nj = pj.size();
    for(PS::S32 j=0; j<nj; j++){
        const PS::F64 mj = pj[j].mass;
        const PS::F64vec rij = pi.pos - pj[j].pos;
        const PS::F64 r_sq = rij * rij + eps_sq;
        const PS::F64 r_inv = 1.0 / sqrt(r_sq);
        const PS::F64 pij = mj * r_inv;
        const double mri3 = r_inv*r_inv*pij;
        force.acc -= mri3 * rij;
        force.pot -= pij;
    }
}

template<class Tptcl>
void CalcCmPosVelLoc(const Tptcl ptcl[],
                     const PS::S64 n,
                     PS::F64 & cm_mass,
                     PS::F64vec & cm_pos,
                     PS::F64vec & cm_vel){
    cm_mass = 0.0;
    cm_pos = cm_vel = 0.0;
    for(PS::S64 i=0; i<n; i++){
        cm_mass += ptcl[i].mass;
        cm_pos  += ptcl[i].mass*ptcl[i].pos;
        cm_vel  += ptcl[i].mass*ptcl[i].vel;
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
}

template<class Tptcl>
void CalcCmPosVelGlb(const Tptcl ptcl[],
                     const PS::S64 n,
                     PS::F64 & cm_mass,
                     PS::F64vec & cm_pos,
                     PS::F64vec & cm_vel){
    PS::F64 cm_mass_loc = 0.0;
    PS::F64vec cm_pos_loc = 0.0;
    PS::F64vec cm_vel_loc = 0.0;
    for(PS::S64 i=0; i<n; i++){
        cm_mass_loc += ptcl[i].mass;
        cm_pos_loc += ptcl[i].mass*ptcl[i].pos;
        cm_vel_loc += ptcl[i].mass*ptcl[i].vel;
    }
    cm_mass = PS::Comm::getSum(cm_mass_loc);
    cm_pos = PS::Comm::getSum(cm_pos_loc);
    cm_vel = PS::Comm::getSum(cm_vel_loc);
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
}

template<class Tptcl>
void CalcCmPosFromList(const Tptcl ptcl[],
                       const PS::S64 id_head,
                       const PS::S64 id_tail,
                       const PS::S32 adr[],
                       PS::F64 & cm_mass,
                       PS::F64vec & cm_pos){
    for(PS::S64 i=id_head; i<id_tail; i++){
        cm_mass += ptcl[ adr[i] ].mass;
        cm_pos  += ptcl[ adr[i] ].mass*ptcl[ adr[i] ].pos;
    }
}

template<class Tepi, class Tepj, class Tspj>
PS::S32 DispatchKernelWithSP(
                             const PS::S32   n_walk,
                             Tepi      epi[],
                             const PS::S32   n_epi[],
                             const PS::S32   n_disp_epi[],
                             const PS::S32   adr_epj[],
                             const PS::S32   n_epj[],
                             const PS::S32   n_disp_epj[],
                             const PS::S32   adr_spj[],
                             const PS::S32   n_spj[],
                             const PS::S32   n_disp_spj[],
                             const Tepj      epj[],
                             const Tspj      spj[]){
#ifdef DEBUG_CHECK_SUM
    struct CheckSumList{
        PS::U64 n_epj_tot;
        PS::U64 n_spj_tot;
        PS::U64 id_sum;
        PS::U64 adr_sum;
    } check_sum_list_loc = {n_disp_epj[n_walk], n_disp_spj[n_walk], 0, 0};
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=n_disp_epj[ig]; jp<n_disp_epj[ig+1]; jp++) {
            check_sum_list_loc.adr_sum+= (PS::U64)adr_epj[jp];
            check_sum_list_loc.id_sum += (PS::U64)epj[adr_epj[jp]].id;
        }
    }
    CheckSumList check_sum_list_glb;
    MPI_Allreduce(&check_sum_list_loc, &check_sum_list_glb, 4, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0){
        std::cerr<<"CHECK SUM (INTERACTION LIST)"<<std::endl;
        std::cerr<<"n_epj_tot_glb= "<<check_sum_list_glb.n_epj_tot<<std::endl;
        std::cerr<<"n_spj_tot_glb= "<<check_sum_list_glb.n_spj_tot<<std::endl;
        std::cerr<<"adr_sum_glb= "<<check_sum_list_glb.adr_sum<<std::endl;
        std::cerr<<"id_sum_glb= "<<check_sum_list_glb.id_sum<<std::endl;
    }
    PS::Comm::barrier();
#endif
    
#ifdef DUMP_SATELLITE
    if(PS::Comm::getRank()==PS::CHECK_RANK_MPI){
        for(PS::S32 i=0; i<N_SATELLITE; i++){
            std::cerr<<"i= "<<i
                     <<" SATELLITE[i].pos= "<<SATELLITE[i].pos
                     <<" vel= "<<SATELLITE[i].vel
                     <<std::endl;
        }
    }
#endif    
    
#ifdef DEBUG_DISPATCH_KERNEL
   // for debug 
    PS::F64 cm_mass_ref = 0.0;;
    PS::F64vec cm_pos_ref = 0.0;
    CalcCmPosGlb(epi, n_disp_epi[n_walk], cm_mass_ref, cm_pos_ref);
    PS::Comm::barrier();
    if(PS::Comm::getRank()==PS::CHECK_RANK_MPI){
        std::cerr<<"cm_mass_ref= "<<cm_mass_ref<<" cm_pos_ref= "<<cm_pos_ref
                 <<" @DispatchKernelWithSP"
                 <<std::endl;
    }
    PS::Comm::barrier();
    bool flag_conserve_cm = true;
    for(PS::S32 i=0; i<n_walk; i++){
        PS::F64 cm_mass = 0.0;
        PS::F64vec cm_pos = 0.0;
        CalcCmPosFromList(epj, n_disp_epj[i], n_disp_epj[i]+n_epj[i], adr_epj, cm_mass, cm_pos);
        CalcCmPosFromList(spj, n_disp_spj[i], n_disp_spj[i]+n_spj[i], adr_spj, cm_mass, cm_pos);
        cm_pos /= cm_mass;
        if(PS::Comm::getRank()==PS::CHECK_RANK_MPI){
            std::cerr<<"i= "<<i
                     <<" cm_pos= "<<cm_pos
                     <<" cm_mass= "<<cm_mass
                     <<" rank= "<<PS::Comm::getRank()
                     <<std::endl;
        }
        /*
        if( fabs(cm_mass_ref-cm_mass)/cm_mass_ref > 1e-10
            || fabs(cm_pos.x-cm_pos_ref.x) > 1e-10
            || fabs(cm_pos.y-cm_pos_ref.y) > 1e-10){
            std::cerr<<"rank= "<<PS::Comm::getRank()
                     <<" i= "<<i
                     <<" cm_mass= "<<cm_mass
                     <<" cm_mass_ref= "<<cm_mass_ref
                     <<" cm_pos= "<<cm_pos
                     <<" cm_pos_ref= "<<cm_pos_ref
                     <<std::endl;
            flag_conserve_cm = false;
        }
        */
    }
    /*
    if(flag_conserve_cm){
        //std::cerr<<"conserve CM rank= "<<PS::Comm::getRank()<<std::endl;
    }
    else{
        std::cerr<<"NOT conserve CM rank= "<<PS::Comm::getRank()<<std::endl;
        assert(0);
    }
    */
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"OK1 @DispatchKernelWithSP"<<std::endl;
    }
    //exit(1);
#endif

#ifdef CHECK_NAN
    bool flag_exit = false;
    PS::S32 num_nan_loc,num_nan;
    PS::S32 n_epi_tot,n_epj_tot,n_spj_tot;
    static int num_called=0;
    num_called++;
    if (PS::Comm::getRank() == 0)
        std::cout << "num_called = " << num_called << std::endl;

    n_epi_tot = 0;
    n_epj_tot = 0;
    n_spj_tot = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
       n_epi_tot += n_epi[ig];
       n_epj_tot += n_epj[ig];
       n_spj_tot += n_spj[ig];
    }
    PS::S32 n_epi_tot_min,n_epi_tot_max;
    PS::S32 n_epj_tot_min,n_epj_tot_max;
    PS::S32 n_spj_tot_min,n_spj_tot_max;
    PS::S32 n_jp_tot_min,n_jp_tot_max;
    PS::S32 ntmp;
    MPI::COMM_WORLD.Allreduce(&n_epi_tot,&n_epi_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&n_epi_tot,&n_epi_tot_max,1,MPI::INT,MPI::MAX);
    MPI::COMM_WORLD.Allreduce(&n_epj_tot,&n_epj_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&n_epj_tot,&n_epj_tot_max,1,MPI::INT,MPI::MAX);
    MPI::COMM_WORLD.Allreduce(&n_spj_tot,&n_spj_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&n_spj_tot,&n_spj_tot_max,1,MPI::INT,MPI::MAX);
    ntmp = n_epj_tot + n_spj_tot;
    MPI::COMM_WORLD.Allreduce(&ntmp,&n_jp_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&ntmp,&n_jp_tot_max,1,MPI::INT,MPI::MAX);
    if (PS::Comm::getRank() == 0) {
        std::cout << "n_epi_tot_min = " << n_epi_tot_min << std::endl;
        std::cout << "n_epi_tot_max = " << n_epi_tot_max << std::endl;
        std::cout << "n_epj_tot_min = " << n_epj_tot_min << std::endl;
        std::cout << "n_epj_tot_max = " << n_epj_tot_max << std::endl;
        std::cout << "n_spj_tot_min = " << n_spj_tot_min << std::endl;
        std::cout << "n_spj_tot_max = " << n_spj_tot_max << std::endl;
        std::cout << "n_jp_tot_min  = " << n_jp_tot_min  << std::endl;
        std::cout << "n_jp_tot_max  = " << n_jp_tot_max  << std::endl;
    }

    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epi1] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            if ((std::isfinite(epj[adr_src].pos.x) != true) ||
                (std::isfinite(epj[adr_src].pos.y) != true) ||
                (std::isfinite(epj[adr_src].pos.z) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epj] " << num_nan << std::endl;
    
    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_spj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            if ((std::isfinite(spj[adr_src].pos.x) != true) ||
                (std::isfinite(spj[adr_src].pos.y) != true) ||
                (std::isfinite(spj[adr_src].pos.z) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-spj] " << num_nan << std::endl;


    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true) ||
            (std::isfinite(epi[i].vel.x) != true) ||
            (std::isfinite(epi[i].vel.y) != true) ||
            (std::isfinite(epi[i].vel.z) != true) ||
            (std::isfinite(epi[i].mass) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "A) [nan-epi1] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            if ((std::isfinite(epj[adr_src].pos.x) != true) ||
                (std::isfinite(epj[adr_src].pos.y) != true) ||
                (std::isfinite(epj[adr_src].pos.z) != true) ||
                (std::isfinite(epj[adr_src].vel.x) != true) ||
                (std::isfinite(epj[adr_src].vel.y) != true) ||
                (std::isfinite(epj[adr_src].vel.z) != true) ||
                (std::isfinite(epj[adr_src].mass) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "A) [nan-epj] " << num_nan << std::endl;
    
    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_spj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            if ((std::isfinite(spj[adr_src].pos.x) != true) ||
                (std::isfinite(spj[adr_src].pos.y) != true) ||
                (std::isfinite(spj[adr_src].pos.z) != true) ||
                (std::isfinite(spj[adr_src].mass) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "A) [nan-spj] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        if ((std::isfinite(SATELLITE[i].pos.x) != true) ||
            (std::isfinite(SATELLITE[i].pos.y) != true) ||
            (std::isfinite(SATELLITE[i].pos.z) != true) ||
            (std::isfinite(SATELLITE[i].vel.x) != true) ||
            (std::isfinite(SATELLITE[i].vel.y) != true) ||
            (std::isfinite(SATELLITE[i].vel.z) != true) ||
            (std::isfinite(SATELLITE[i].mass)  != true)) {
            num_nan_loc++;
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "A) [nan-sat] " << num_nan << std::endl;
    
    //if (num_nan_loc > 0) {
    //    std::cout << "myrank = " << PS::Comm::getRank() << " "
    //              << "ntot = " << n_epi_tot + ntmp << std::endl;
    //}
#endif // CHECKNAN
    //******** debug [end] *********

#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
    // [*** Note ****]
    //   Long-Wang's force kernel is implemented here.
    //   cpe_pars below is a global variable defined at nbody.cpp.
    //* Set force kernel parameters
    cpe_pars.epi        = (EpiMM*) epi;
    cpe_pars.epj        = (EpjMM*) epj;
    cpe_pars.spj        = (SpjMM*) spj;
    cpe_pars.sat        = (EpiMM*) SATELLITE.getPointer();
    cpe_pars.force_sat  = (ForceMM*) SATELLITE_FORCE_LOC.getPointer();
    cpe_pars.n_sat      = N_SATELLITE;
    cpe_pars.n_disp_epi = const_cast<int*>(n_disp_epi);
    cpe_pars.n_disp_epj = const_cast<int*>(n_disp_epj);
    cpe_pars.n_disp_spj = const_cast<int*>(n_disp_spj);
    cpe_pars.n_epi      = const_cast<int*>(n_epi);
    cpe_pars.n_epj      = const_cast<int*>(n_epj);
    cpe_pars.n_spj      = const_cast<int*>(n_spj);
    cpe_pars.adr_epj    = const_cast<int*>(adr_epj);
    cpe_pars.adr_spj    = const_cast<int*>(adr_spj);
    cpe_pars.n_walk     = n_walk;
    #ifdef SUNWAY_FORCE_KERNEL_CHECK
    cpe_pars.forcei     = (ForceMM*)force_mm;
    #endif

#ifdef DEBUG_FORCE_OUTPUT
    //******** output [start] **********
    //if (PS::Comm::getRank() == 27) {
    //if (PS::Comm::getRank() == 0) {
    if (PS::Comm::getRank() == 26) {
        for(PS::S32 jp=0; jp<10; jp++){
            std::cerr<<"epj.pos= "<<epj[jp].pos
                     <<" epj.mass= "<<epj[jp].mass
                     <<" epj.vel= "<<epj[jp].vel
                     <<" epj.id= "<<epj[jp].id
                     <<std::endl;
        }
        //for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 ig=0; ig<1; ig++) {
            for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
                const PS::S32 adr_des = n_disp_epj[ig]+jp;
                const PS::S32 adr_src = adr_epj[adr_des];
                std::cerr<<"epj[adr_src].pos= "<<epj[adr_src].pos
                         <<" epj[adr_src].mass= "<<epj[adr_src].mass
                         <<" epj[adr_src].vel= "<<epj[adr_src].vel
                         <<" epj[adr_src].id= "<<epj[adr_src].id
                         <<std::endl;
            }
        }
        
        std::ofstream ofs;
        std::stringstream fnum;
        std::string fname;
        fnum << std::setfill('0') << std::setw(5) << PS::Comm::getRank();
        fname = "./output/intlst" + fnum.str() + ".dat";
        ofs.open(fname.c_str(), std::ios::binary|std::ios::trunc);
        //* Compute the sizes of arrays
        PS::S32 size_of_adr_epj=0, size_of_adr_spj=0;
        PS::S32 num_epi=0, num_epj=0, num_spj=0;
        for (PS::S32 ig=0; ig<n_walk; ig++) {
            //** epi
            num_epi += n_epi[ig];
            //** epj
            for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
                const PS::S32 adr_des = n_disp_epj[ig]+jp;
                const PS::S32 adr_src = adr_epj[adr_des];
                if (adr_des > size_of_adr_epj) size_of_adr_epj=adr_des;
                if (adr_src > num_epj) num_epj=adr_src;
            }
            //** spj
            for (PS::S32 jp=0; jp<n_spj[ig]; jp++) {
                const PS::S32 adr_des = n_disp_spj[ig]+jp;
                const PS::S32 adr_src = adr_spj[adr_des];
                if (adr_des > size_of_adr_spj) size_of_adr_spj=adr_des;
                if (adr_src > num_spj) num_spj=adr_src;
            }
        }
        size_of_adr_epj++;
        size_of_adr_spj++;
        num_epj++;
        num_spj++;
        //* Information required for the indirect access
        ofs.write((char *)&n_walk, sizeof(PS::S32));
        ofs.write((char *)&n_disp_epi[0], sizeof(PS::S32) * n_walk);
        ofs.write((char *)&n_disp_epj[0], sizeof(PS::S32) * n_walk);
        ofs.write((char *)&n_disp_spj[0], sizeof(PS::S32) * n_walk);
        ofs.write((char *)&n_epi[0], sizeof(PS::S32) * n_walk);
        ofs.write((char *)&n_epj[0], sizeof(PS::S32) * n_walk);
        ofs.write((char *)&n_spj[0], sizeof(PS::S32) * n_walk);
        ofs.write((char *)&size_of_adr_epj, sizeof(PS::S32));
        ofs.write((char *)&size_of_adr_spj, sizeof(PS::S32));
        ofs.write((char *)&adr_epj[0], sizeof(PS::S32) * size_of_adr_epj);
        ofs.write((char *)&adr_spj[0], sizeof(PS::S32) * size_of_adr_spj);
        //* Particle data
        size_t bsize_epi = sizeof(Epi);
        size_t bsize_epj = sizeof(Epj);
        size_t bsize_spj = sizeof(SPJMonopoleScatterSW);
        ofs.write((char *)&num_epi, sizeof(PS::S32));
        ofs.write((char *)&num_epj, sizeof(PS::S32));
        ofs.write((char *)&num_spj, sizeof(PS::S32));
        ofs.write((char *)&epi[0], bsize_epi * num_epi);
        ofs.write((char *)&epj[0], bsize_epj * num_epj);
        ofs.write((char *)&spj[0], bsize_spj * num_spj);
        //* Satellite data
        size_t bsize_sat = sizeof(Planet);
        ofs.write((char *)&N_SATELLITE, sizeof(PS::S32));
        ofs.write((char *)SATELLITE.getPointer(), bsize_sat * N_SATELLITE);
        ofs.close();
    }
    //******** output [end] **********
#endif

#ifdef DEBUG_DISPATCH_KERNEL
    /*
    if(PS::Comm::getRank()==59 && N_LOOP_GLB==2){
        for(PS::S32 i=0; i<n_walk; i++){
            std::cerr<<"i= "<<i<<" n_walk= "<<n_walk
                     <<" adr_n_walk[i]= "<<PS::adr_n_walk[i]
                     <<" n_walk_cpe[i]= "<<PS::n_walk_cpe[i]
                     <<" n_disp_walk[i]= "<<PS::n_disp_walk[i]
                     <<std::endl;
            for(PS::S32 j=n_disp_epi[i]; j<n_disp_epi[i]+n_epi[i]; j++){
                std::cerr<<"j= "<<j
                         <<" epi[j].id= "<<epi[j].id
                         <<" epi[j].mass= "<<epi[j].mass
                         <<" epi[j].pos= "<<epi[j].pos
                         <<std::endl;
            }
            for(PS::S32 j=n_disp_epj[i]; j<n_disp_epj[i]+n_epj[i]; j++){
                PS::S32 adr = adr_epj[j];
                std::cerr<<"j= "<<j
                         <<" adr= "<<adr
                         <<" epj[adr].mass= "<<epj[adr].mass
                         <<" epj[adr].pos= "<<epj[adr].pos
                         <<std::endl;
            }
            for(PS::S32 j=n_disp_spj[i]; j<n_disp_spj[i]+n_spj[i]; j++){
                PS::S32 adr = adr_spj[j];
                std::cerr<<"j= "<<j
                         <<" adr= "<<adr
                         <<" spj[adr].mass= "<<spj[adr].mass
                         <<" spj[adr].pos= "<<spj[adr].pos
                         <<std::endl;
            }        
        }
    }
    */
#endif

#if 0
    {
    PS::Comm::barrier();
    PS::F64 cm_mass_loc = 0.0;
    PS::F64vec cm_pos_loc = 0.0;
    PS::F64vec cm_vel_loc = 0.0;
    CalcCmPosVelLoc(epi, n_disp_epi[n_walk], cm_mass_loc, cm_pos_loc, cm_vel_loc);
    PS::F64 cm_mass_glb = 0.0;
    PS::F64vec cm_pos_glb = 0.0;
    PS::F64vec cm_vel_glb = 0.0;
    CalcCmPosVelGlb(epi, n_disp_epi[n_walk], cm_mass_glb, cm_pos_glb, cm_vel_glb);
    PS::Comm::barrier();
    PS::F64    * cm_mass_array = new PS::F64[PS::Comm::getNumberOfProc()];
    PS::F64vec * cm_pos_array = new PS::F64vec[PS::Comm::getNumberOfProc()];
    PS::F64vec * cm_vel_array = new PS::F64vec[PS::Comm::getNumberOfProc()];
    PS::Comm::gather(&cm_mass_loc, 1, cm_mass_array);
    PS::Comm::gather(&cm_pos_loc, 1, cm_pos_array);
    PS::Comm::gather(&cm_vel_loc, 1, cm_vel_array);
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"before 1st foce"<<std::endl;
        std::cerr<<" cm_mass_glb= "<<cm_mass_glb
                 <<" cm_pos_glb= "<<cm_pos_glb
                 <<" cm_vel_glb= "<<cm_vel_glb
                 <<std::endl;
        for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++){
            printf("i= %d, cm_mass= %a, pos= %a %a %a, vel= %a %a %a \n",
                   i, cm_mass_array[i],
                   cm_pos_array[i].x, cm_pos_array[i].y, cm_pos_array[i].z,
                   cm_vel_array[i].x, cm_vel_array[i].y, cm_vel_array[i].z);
        }
    }
    delete [] cm_mass_array;
    delete [] cm_pos_array;
    delete [] cm_vel_array;
    PS::Comm::barrier();
    }
#endif

#ifdef CHECK_ANGULAR_MOMENTUM_2
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0) std::cerr<<"before 1st force"<<std::endl;
    angmom_test0(epi, n_disp_epi[n_walk]);
    PS::Comm::barrier();
#endif

#ifdef DEBUG_PRINT_DISPATCH_KERNEL
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"before 1st force @DispatchKernelWithSP"
                 <<std::endl;
    }
#endif
    
    PS::F64 starttime = MPI::Wtime();
    __real_athread_spawn((void*)slave_ForceKernelSunWay1st,(void*)&cpe_pars);
    athread_join();
    PS::wtime_calc_force += MPI::Wtime() - starttime;

#ifdef DEBUG_PRINT_DISPATCH_KERNEL
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"after 1st force @DispatchKernelWithSP"
                 <<std::endl;
    }
#endif
    
#ifdef CHECK_ANGULAR_MOMENTUM_2
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0) std::cerr<<"after 1st force"<<std::endl;
    angmom_test0(epi, n_disp_epi[n_walk]);
    PS::Comm::barrier();
#endif

#if 0
    {
    PS::Comm::barrier();
    PS::F64 cm_mass_loc = 0.0;
    PS::F64vec cm_pos_loc = 0.0;
    PS::F64vec cm_vel_loc = 0.0;
    CalcCmPosVelLoc(epi, n_disp_epi[n_walk], cm_mass_loc, cm_pos_loc, cm_vel_loc);
    PS::F64 cm_mass_glb = 0.0;
    PS::F64vec cm_pos_glb = 0.0;
    PS::F64vec cm_vel_glb = 0.0;
    CalcCmPosVelGlb(epi, n_disp_epi[n_walk], cm_mass_glb, cm_pos_glb, cm_vel_glb);
    PS::Comm::barrier();
    PS::F64    * cm_mass_array = new PS::F64[PS::Comm::getNumberOfProc()];
    PS::F64vec * cm_pos_array = new PS::F64vec[PS::Comm::getNumberOfProc()];
    PS::F64vec * cm_vel_array = new PS::F64vec[PS::Comm::getNumberOfProc()];
    PS::Comm::gather(&cm_mass_loc, 1, cm_mass_array);
    PS::Comm::gather(&cm_pos_loc, 1, cm_pos_array);
    PS::Comm::gather(&cm_vel_loc, 1, cm_vel_array);
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"after 1st foce"<<std::endl;
        std::cerr<<" cm_mass_glb= "<<cm_mass_glb
                 <<" cm_pos_glb= "<<cm_pos_glb
                 <<" cm_vel_glb= "<<cm_vel_glb
                 <<std::endl;
        for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++){
            printf("i= %d, cm_mass= %a, pos= %a %a %a, vel= %a %a %a \n",
                   i, cm_mass_array[i],
                   cm_pos_array[i].x, cm_pos_array[i].y, cm_pos_array[i].z,
                   cm_vel_array[i].x, cm_vel_array[i].y, cm_vel_array[i].z);
        }
    }
    delete [] cm_mass_array;
    delete [] cm_pos_array;
    delete [] cm_vel_array;
    PS::Comm::barrier();
    }
#endif

    
#ifdef CHECK_NAN
    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true) ||
            (std::isfinite(epi[i].vel.x) != true) ||
            (std::isfinite(epi[i].vel.y) != true) ||
            (std::isfinite(epi[i].vel.z) != true) ||
            (std::isfinite(epi[i].mass) != true)) {
            num_nan_loc++;
            flag_exit = true;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "B) [nan-epi1] " << num_nan << std::endl;


    
    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            if ((std::isfinite(epj[adr_src].pos.x) != true) ||
                (std::isfinite(epj[adr_src].pos.y) != true) ||
                (std::isfinite(epj[adr_src].pos.z) != true) ||
                (std::isfinite(epj[adr_src].vel.x) != true) ||
                (std::isfinite(epj[adr_src].vel.y) != true) ||
                (std::isfinite(epj[adr_src].vel.z) != true) ||
                (std::isfinite(epj[adr_src].mass) != true)) {
                num_nan_loc++;
                flag_exit = true;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "B) [nan-epj] " << num_nan << std::endl;
    
    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_spj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            if ((std::isfinite(spj[adr_src].pos.x) != true) ||
                (std::isfinite(spj[adr_src].pos.y) != true) ||
                (std::isfinite(spj[adr_src].pos.z) != true) ||
                (std::isfinite(spj[adr_src].mass) != true)) {
                num_nan_loc++;
                flag_exit = true;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "B) [nan-spj] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        if ((std::isfinite(SATELLITE[i].pos.x) != true) ||
            (std::isfinite(SATELLITE[i].pos.y) != true) ||
            (std::isfinite(SATELLITE[i].pos.z) != true) ||
            (std::isfinite(SATELLITE[i].vel.x) != true) ||
            (std::isfinite(SATELLITE[i].vel.y) != true) ||
            (std::isfinite(SATELLITE[i].vel.z) != true) ||
            (std::isfinite(SATELLITE[i].mass)  != true)) {
            num_nan_loc++;
            flag_exit = true;
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "B) [nan-sat] " << num_nan << std::endl;

    num_nan_loc = 0;
    for(PS::S32 i=0; i<N_SATELLITE*64; i++){
        if( (std::isfinite(SATELLITE_FORCE_LOC[i].acc.x) != true) ||
            (std::isfinite(SATELLITE_FORCE_LOC[i].acc.y) != true) ||
            (std::isfinite(SATELLITE_FORCE_LOC[i].acc.z) != true) ||
            (std::isfinite(SATELLITE_FORCE_LOC[i].pot) != true) ){
            std::cerr<<"rank="<<PS::Comm::getRank()<<" i= "<<i<<std::endl;
            /*
            if(PS::Comm::getRank() == 0) std::cerr<<"rank0: i= "<<i<<std::endl;
            if(PS::Comm::getRank() == 1) std::cerr<<"rank1: i= "<<i<<std::endl;
            */
            num_nan_loc++;
            flag_exit = true;
        }
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "B) [nan-satforce] " << num_nan << std::endl;

    PS::Comm::barrier();
    if(flag_exit) exit(1);
#endif
    
    //* Reduce satellite force 
    PS::S32 s = LapTimer::step_counter();
    LapTimer::measure(s, LapTimer::BEFORE_REDUCE);
#ifdef MULTI_ALLREDUCE
    /*
    double t0 = PS::GetWtime();
    MPI_Allreduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                   (double*)(SATELLITE_FORCE.getPointer()),
                   4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double t1 = PS::GetWtime();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"A) wtime allreuce: "<<t1-t0<<std::endl;
        for(PS::S32 i=0; i<SATELLITE.size(); i++){
            std::cerr<<"A) i= "<<i<<" SATELLITE_FORCE[i].acc= "<<SATELLITE_FORCE[i].acc<<std::endl;
        }
    }
    */

    //t0 = PS::GetWtime();
    MPI_Reduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                (double*)(SATELLITE_FORCE.getPointer()),
                4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, 0, MY_MPI_COMM_SUB);
    if(MY_RANK_SUB == 0){
        MPI_Allreduce( (double*)(SATELLITE_FORCE.getPointer()),
                       (double*)(SATELLITE_FORCE_LOC.getPointer()),
                       4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MY_MPI_COMM_1D);
    }
    MPI_Bcast((double*)(SATELLITE_FORCE_LOC.getPointer()), 4*SATELLITE.size(),
              MPI_DOUBLE, 0, MY_MPI_COMM_SUB);
    for(PS::S32 i=0; i<SATELLITE.size(); i++){
        SATELLITE_FORCE[i] = SATELLITE_FORCE_LOC[i];
    }
    //t1 = PS::GetWtime();
    /*
    if (PS::Comm::getRank() == 0){
        std::cerr<<"B) wtime allreuce: "<<t1-t0<<std::endl;
        for(PS::S32 i=0; i<SATELLITE.size(); i++){
            std::cerr<<"B) i= "<<i<<" SATELLITE_FORCE[i].acc= "<<SATELLITE_FORCE[i].acc<<std::endl;
        }
    }
    */
    //exit(1);
    
#else
    
    MPI_Allreduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                   (double*)(SATELLITE_FORCE.getPointer()),
                   4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    #if 0
    {
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0){
            std::cerr<<"after allreduce, my_rank= "<<PS::Comm::getRank()<<std::endl;
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_pos= %a %a %a \n",
                       i, SATELLITE[i].pos.x, SATELLITE[i].pos.y, SATELLITE[i].pos.z);
            }
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_vel= %a %a %a \n",
                       i, SATELLITE[i].vel.x, SATELLITE[i].vel.y, SATELLITE[i].vel.z);
            }
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_force= %a %a %a \n",
                       i, SATELLITE_FORCE[i].acc.x, SATELLITE_FORCE[i].acc.y, SATELLITE_FORCE[i].acc.z);
            }
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_force_loc= %a %a %a \n",
                       i, SATELLITE_FORCE_LOC[i].acc.x, SATELLITE_FORCE_LOC[i].acc.y, SATELLITE_FORCE_LOC[i].acc.z);
            }
        }
        PS::Comm::barrier();
    }
    #endif
    
#endif
    
#ifdef DEBUG_PRINT_DISPATCH_KERNEL
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"after allreduce force on sat @DispatchKernelWithSP"
                 <<std::endl;
    }
#endif
    
    LapTimer::measure(s, LapTimer::AFTER_REDUCE);

#ifdef CHECK_NAN
    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true) ||
            (std::isfinite(epi[i].vel.x) != true) ||
            (std::isfinite(epi[i].vel.y) != true) ||
            (std::isfinite(epi[i].vel.z) != true) ||
            (std::isfinite(epi[i].mass) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "C) [nan-epi1] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            if ((std::isfinite(epj[adr_src].pos.x) != true) ||
                (std::isfinite(epj[adr_src].pos.y) != true) ||
                (std::isfinite(epj[adr_src].pos.z) != true) ||
                (std::isfinite(epj[adr_src].vel.x) != true) ||
                (std::isfinite(epj[adr_src].vel.y) != true) ||
                (std::isfinite(epj[adr_src].vel.z) != true) ||
                (std::isfinite(epj[adr_src].mass) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "C) [nan-epj] " << num_nan << std::endl;
    
    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_spj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            if ((std::isfinite(spj[adr_src].pos.x) != true) ||
                (std::isfinite(spj[adr_src].pos.y) != true) ||
                (std::isfinite(spj[adr_src].pos.z) != true) ||
                (std::isfinite(spj[adr_src].mass) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "C) [nan-spj] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        if ((std::isfinite(SATELLITE[i].pos.x) != true) ||
            (std::isfinite(SATELLITE[i].pos.y) != true) ||
            (std::isfinite(SATELLITE[i].pos.z) != true) ||
            (std::isfinite(SATELLITE[i].vel.x) != true) ||
            (std::isfinite(SATELLITE[i].vel.y) != true) ||
            (std::isfinite(SATELLITE[i].vel.z) != true) ||
            (std::isfinite(SATELLITE[i].mass)  != true)) {

            num_nan_loc++;
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "C) [nan-sat] " << num_nan << std::endl;
    //exit(1);
#endif
    
    /*
#ifdef CHECK_NAN
    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epi C] " << num_nan << std::endl;
#endif
    */
    
    //* Compute satellite-satellite force and perform kick-drift for satellite
    cpe_pars.sat        = (EpiMM*) SATELLITE.getPointer();
    cpe_pars.force_sat  = (ForceMM*) SATELLITE_FORCE.getPointer();
    cpe_pars.n_sat      = N_SATELLITE; 
    __real_athread_spawn((void*)slave_ForceKernelSunWay2nd,(void*)&cpe_pars);
    athread_join();

#ifdef CHECK_ANGULAR_MOMENTUM_2
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0) std::cerr<<"after 2nd force"<<std::endl;
    PS::Comm::barrier();
    angmom_test0(epi, n_disp_epi[n_walk]);
#endif


#if 0
    {
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0){
            std::cerr<<"after 2nd force, my_rank= "<<PS::Comm::getRank()<<std::endl;
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_pos= %a %a %a \n",
                       i, SATELLITE[i].pos.x, SATELLITE[i].pos.y, SATELLITE[i].pos.z);
            }
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_vel= %a %a %a \n",
                       i, SATELLITE[i].vel.x, SATELLITE[i].vel.y, SATELLITE[i].vel.z);
            }
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_force= %a %a %a \n",
                       i, SATELLITE_FORCE[i].acc.x, SATELLITE_FORCE[i].acc.y, SATELLITE_FORCE[i].acc.z);
            }
            for(PS::S32 i=0; i<SATELLITE.size(); i++){
                printf("i= %d, sat_force_loc= %a %a %a \n",
                       i, SATELLITE_FORCE_LOC[i].acc.x, SATELLITE_FORCE_LOC[i].acc.y, SATELLITE_FORCE_LOC[i].acc.z);
            }
        }
        PS::Comm::barrier();
    }
    //exit(1);
#endif
    
#ifdef DEBUG_PRINT_DISPATCH_KERNEL
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) {
        std::cerr<<"after ForceKernelSunWay2nd @DispatchKernelWithSP"
                 <<std::endl;
    }
#endif
    
#ifdef CHECK_NAN

    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true) ||
            (std::isfinite(epi[i].vel.x) != true) ||
            (std::isfinite(epi[i].vel.y) != true) ||
            (std::isfinite(epi[i].vel.z) != true) ||
            (std::isfinite(epi[i].mass) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "D) [nan-epi1] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            if ((std::isfinite(epj[adr_src].pos.x) != true) ||
                (std::isfinite(epj[adr_src].pos.y) != true) ||
                (std::isfinite(epj[adr_src].pos.z) != true) ||
                (std::isfinite(epj[adr_src].vel.x) != true) ||
                (std::isfinite(epj[adr_src].vel.y) != true) ||
                (std::isfinite(epj[adr_src].vel.z) != true) ||
                (std::isfinite(epj[adr_src].mass) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "D) [nan-epj] " << num_nan << std::endl;
    
    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_spj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            if ((std::isfinite(spj[adr_src].pos.x) != true) ||
                (std::isfinite(spj[adr_src].pos.y) != true) ||
                (std::isfinite(spj[adr_src].pos.z) != true) ||
                (std::isfinite(spj[adr_src].mass) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "D) [nan-spj] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        if ((std::isfinite(SATELLITE[i].pos.x) != true) ||
            (std::isfinite(SATELLITE[i].pos.y) != true) ||
            (std::isfinite(SATELLITE[i].pos.z) != true) ||
            (std::isfinite(SATELLITE[i].vel.x) != true) ||
            (std::isfinite(SATELLITE[i].vel.y) != true) ||
            (std::isfinite(SATELLITE[i].vel.z) != true) ||
            (std::isfinite(SATELLITE[i].mass)  != true)) {
            num_nan_loc++;
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "D) [nan-sat] " << num_nan << std::endl;
    
    /*
    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epi d] " << num_nan << std::endl;
    */
#endif
    
    //* Copy
    EPI_POINTER = epi;

    //* For test-run
    //fout_wtime.close();
    //athread_halt();
    //PS::Finalize();
    //std::exit(0);
#elif defined(MPE_FORCE_KERNEL)
    EpiMM* sat_mpe = (EpiMM*)SATELLITE.getPointer();
    ForceMM* force_sat_mpe = (ForceMM*)SATELLITE_FORCE_LOC.getPointer();
    memset(force_sat_mpe,0,N_SATELLITE*sizeof(Force));
    std::cerr<<"Force calculation on MPE\n";
    //      wtime_copy0 = PS::GetWtime();
    for(int iw=0; iw<n_walk; iw++){
      Force* force_mpe = &FORCE_ARRAY[n_disp_epi[iw]];
      memset(force_mpe,0,n_epi[iw]*sizeof(Force));
      //      std::cerr<<"iw="<<iw<<" I-J Ni="<<n_epi[iw]<<" Nj="<<n_epj[iw]<<"\n";
      calcDirectGrav((EpiMM*)&(epi[n_disp_epi[iw]]), n_epi[iw], 
                     (EpjMM*)epj, n_epj[iw],(int*)&adr_epj[n_disp_epj[iw]],
                     (ForceMM*)force_mpe, mpe_pars.r_ep, mpe_pars.r_ep, mpe_pars.kappa, mpe_pars.eta);
      //      std::cerr<<"iw="<<iw<<" I-SP\n";
      calcSPGrav((EpiMM*)&(epi[n_disp_epi[iw]]),n_epi[iw],
                 (SpjMM*)spj, n_spj[iw],(int*)&adr_spj[n_disp_spj[iw]],
                 (ForceMM*)force_mpe);
      //      std::cerr<<"iw="<<iw<<" I-SAT\n";
      calcDirectGrav((EpiMM*)&(epi[n_disp_epi[iw]]), n_epi[iw],
                     (EpjMM*)sat_mpe, N_SATELLITE, NULL,
                     (ForceMM*)force_mpe, mpe_pars.r_ep, mpe_pars.r_sat, mpe_pars.kappa, mpe_pars.eta, false);
      //      std::cerr<<"iw="<<iw<<" SAT-I\n";
      calcDirectGrav((EpiMM*)sat_mpe, N_SATELLITE,
                     (EpjMM*)&(epi[n_disp_epi[iw]]), n_epi[iw], NULL,
                     (ForceMM*)force_sat_mpe, mpe_pars.r_ep, mpe_pars.r_sat, mpe_pars.kappa, mpe_pars.eta, false);
      //      std::cerr<<"iw="<<iw<<" I-P\n";
      calcForcePlanet((EpiMM*)&(epi[n_disp_epi[iw]]), n_epi[iw], PLANET.mass, (ForceMM*)force_mpe);
    }
    //    std::cerr<<"Sat MPI Allreduce\n";
    MPI_Allreduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                   (double*)(SATELLITE_FORCE.getPointer()),
                   4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    //    std::cerr<<"Sat-Sat\n";
    calcDirectGrav((EpiMM*)sat_mpe, N_SATELLITE,
                   (EpjMM*)sat_mpe, N_SATELLITE, NULL,
                   (ForceMM*)force_sat_mpe, mpe_pars.r_ep, mpe_pars.r_sat, mpe_pars.kappa, mpe_pars.eta,false);
    //    std::cerr<<"Sat-Planet\n";
    calcForcePlanet(sat_mpe, N_SATELLITE, PLANET.mass, force_sat_mpe);

    //      wtime_copy0 = PS::GetWtime() - wtime_copy0;
    //      fprintf(stderr,"host force time per loop=%e \n",wtime_copy0);

    //* Copy
    EPI_POINTER = epi;
    
#else
    PS::F64 wtime_offset_out = PS::GetWtime();
    const PS::F64 eps_sq = Tepi::eps * Tepi::eps;
    assert(N_EPI_LIMIT > n_disp_epi[n_walk]);
    assert(N_EPJ_LIMIT > n_disp_epj[n_walk]);
    assert(N_SPJ_LIMIT > n_disp_spj[n_walk]);
    const PS::S32 n_epi_tot = n_disp_epi[n_walk];
    
    #ifndef CANCEL_FORCE    
    /////////////////////
    // force on satellite
    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        //SATELLITE_FORCE[ip].clear();
        SATELLITE_FORCE_LOC[ip].clear();
        for(PS::S32 jp=0; jp<n_epi_tot; jp++){
            CalcGravPair(SATELLITE[ip], epi[jp], SATELLITE_FORCE_LOC[ip], eps_sq);
        }
    }
    MPI_Allreduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                   (double*)(SATELLITE_FORCE.getPointer()),
                   4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        for(PS::S32 jp=0; jp<SATELLITE.size(); jp++){
            if(ip == jp) continue;
            CalcGravPair(SATELLITE[ip], SATELLITE[jp], SATELLITE_FORCE[ip], eps_sq);
        }
        CalcGravEx(SATELLITE[ip], PLANET, SATELLITE_FORCE[ip], 0.0);
    }
    // force on satellite
    /////////////////////
    #endif //CANCEL_FORCE

    /////////////////////
    // force on particle
    for(PS::S32 ig=0; ig<n_walk; ig++){
        PS::F64 cm_mass = 0.0;
        PS::F64vec cm_pos = 0.0;
        //PS::F64 wtime_offset = PS::GetWtime();
        for(PS::S32 jp=0; jp<n_epj[ig]; jp++){
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            EPJ_ARRAY[adr_des].id = epj[adr_src].id;
            EPJ_ARRAY[adr_des].mass = epj[adr_src].mass;
            EPJ_ARRAY[adr_des].pos = epj[adr_src].pos;
            cm_mass += EPJ_ARRAY[adr_des].mass;
            cm_pos += EPJ_ARRAY[adr_des].mass*EPJ_ARRAY[adr_des].pos;
        }
        //wtime_copy_all_epj += PS::GetWtime() - wtime_offset;
        //wtime_offset = PS::GetWtime();
        for(PS::S32 jp=0; jp<n_spj[ig]; jp++){
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            SPJ_ARRAY[adr_des].mass = spj[adr_src].mass;
            SPJ_ARRAY[adr_des].pos  = spj[adr_src].pos;
            cm_mass += SPJ_ARRAY[adr_des].mass;
            cm_pos += SPJ_ARRAY[adr_des].mass*SPJ_ARRAY[adr_des].pos;
        }
        //wtime_copy_all_spj += PS::GetWtime() - wtime_offset;        
        //if(PS::Comm::getRank() == 0) std::cerr<<"mass= "<<mass<<std::endl;
        //wtime_offset = PS::GetWtime();
        for(PS::S32 ip=0; ip<n_epi[ig]; ip++){
            cm_mass = 0.0;
            cm_pos = 0.0;
            const PS::S32 adr_i = n_disp_epi[ig] + ip;
            const Tepi & iptcl = epi[adr_i];
            Force & iforce = FORCE_ARRAY[adr_i];
            iforce.acc = 0.0;
            iforce.pot = 0.0;
            
    #ifndef CANCEL_FORCE
            for(PS::S32 jp=0; jp<n_epj[ig]; jp++){
                const PS::S32 adr_j = n_disp_epj[ig] + jp;
                cm_mass += EPJ_ARRAY[adr_j].mass;
                cm_pos += EPJ_ARRAY[adr_j].mass*EPJ_ARRAY[adr_j].pos;
                if(iptcl.id == EPJ_ARRAY[adr_j].id) continue;
                /*
                if(iptcl.id == 3881){
                    PS::F64vec rij = (iptcl.pos - EPJ_ARRAY[adr_j].pos);
                    std::cerr<<"|rij(epj)|= "<<sqrt( rij*rij )<<std::endl;
                }
                */
                CalcGravPair(iptcl, EPJ_ARRAY[adr_j], iforce, eps_sq);
            }
            for(PS::S32 jp=0; jp<n_spj[ig]; jp++){
                const PS::S32 adr_j = n_disp_spj[ig] + jp;
                cm_mass += SPJ_ARRAY[adr_j].mass;
                cm_pos += SPJ_ARRAY[adr_j].mass*SPJ_ARRAY[adr_j].pos;
                if(iptcl.id == 3881){
                    PS::F64vec rij = (iptcl.pos - SPJ_ARRAY[adr_j].pos);
                    //std::cerr<<"SPJ_ARRAY[adr_j].mass= "<<SPJ_ARRAY[adr_j].mass<<" |rij(spj)|= "<<sqrt( rij*rij )<<std::endl;
                }                
                CalcGravPair(iptcl, SPJ_ARRAY[adr_j], iforce, eps_sq);
            }
            for(PS::S32 jp=0; jp<SATELLITE.size(); jp++){
                CalcGravPair(iptcl, SATELLITE[jp], iforce, eps_sq);
            }
    #endif //CANCEL_FORCE
    #ifndef CANCEL_PLANET_FORCE            
            CalcGravEx(iptcl, PLANET, iforce, 0.0);
    #endif
            //if(iptcl.id == 3881) std::cerr<<"cm_mass= "<<cm_mass<<" cm_pos= "<<cm_pos<<std::endl;
        }
        //wtime_calc_pj += PS::GetWtime() - wtime_offset;
    }
    //* Copy
    EPI_POINTER = epi;
    //wtime_dispatch += PS::GetWtime() - wtime_offset_out;
#endif // End of SUNWAY_FORCE_KERNEL

    //****** debug[start] *******
#ifdef CHECK_NAN
    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epi2] " << num_nan << std::endl;
#endif
    //******* debug[end] *******


}

template<class Tforce>
PS::S32 RetrieveKernel(const PS::S32 ni_tot,
                       Tforce force[]){
    /*
    //PS::F64 wtime_offset = PS::GetWtime();
    for(PS::S32 i=0; i<ni_tot; i++){
        force[i].acc = FORCE_ARRAY[i].acc;
        force[i].pot = FORCE_ARRAY[i].pot;
    }
    //wtime_retrieve += PS::GetWtime() - wtime_offset;
    */
    return 0;
}


// origingal hozon
template<class Tepi, class Tepj, class Tspj>
PS::S32 DispatchKernelWithSPKickDrift(
                                      const PS::S32   n_walk,
                                      Tepi      epi[],
                                      const PS::S32   n_epi[],
                                      const PS::S32   n_disp_epi[],
                                      const PS::S32   adr_epj[],
                                      const PS::S32   n_epj[],
                                      const PS::S32   n_disp_epj[],
                                      const PS::S32   adr_spj[],
                                      const PS::S32   n_spj[],
                                      const PS::S32   n_disp_spj[],
                                      const Tepj      epj[],
                                      const Tspj      spj[]){
    #if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
    DispatchKernelWithSP(n_walk, epi, n_epi, n_disp_epi,
                         adr_epj, n_epj, n_disp_epj,
                         adr_spj, n_spj, n_disp_spj,
                         epj, spj);
    //if(N_LOOP_GLB == 2) exit(1);
    #else
    DispatchKernelWithSP(n_walk, epi, n_epi, n_disp_epi,
                         adr_epj, n_epj, n_disp_epj,
                         adr_spj, n_spj, n_disp_spj,
                         epj, spj);
    const PS::S32 n_epi_tot = n_disp_epi[n_walk];
    const PS::F64 dt_kick = (DISPATCH_MODE == 0) ? DT_TREE : DT_TREE*0.5;
    PS::F64 wtime_offset = PS::GetWtime();
    for(PS::S32 ip=0; ip<n_epi_tot; ip++){
        epi[ip].vel += dt_kick * FORCE_ARRAY[ip].acc;
    }
    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        SATELLITE[ip].vel += dt_kick * SATELLITE_FORCE[ip].acc;
    }
    wtime_kick += PS::GetWtime() - wtime_offset;
    if(DISPATCH_MODE == 1){
        wtime_offset = PS::GetWtime();
        CalcEnergy(epi, FORCE_ARRAY, n_epi_tot,
                   SATELLITE.getPointer(), SATELLITE_FORCE.getPointer(), N_SATELLITE,
                   ENERGY_TOT, ENERGY_KIN, ENERGY_POT);
        wtime_calc_energy += PS::GetWtime() - wtime_offset;
        wtime_offset = PS::GetWtime();
        for(PS::S32 ip=0; ip<n_epi_tot; ip++){
            epi[ip].vel += dt_kick * FORCE_ARRAY[ip].acc;
        }
        for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
            SATELLITE[ip].vel += dt_kick * SATELLITE_FORCE[ip].acc;
        }
        wtime_kick += PS::GetWtime() - wtime_offset;
    }
    // TO DO; collision procedure
    wtime_offset = PS::GetWtime();    
    for(PS::S32 ip=0; ip<n_epi_tot; ip++){
        epi[ip].pos += DT_TREE * epi[ip].vel;
        //EPI_ARRAY[ip].pos = epi[ip].pos;
    }
    EPI_POINTER = epi;
    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        SATELLITE[ip].pos += DT_TREE * SATELLITE[ip].vel;
    }
    wtime_drift += PS::GetWtime() - wtime_offset;
    #endif
}



template<class Tsys, class Tsat>
void CalcEnergy(const Tsys & ptcl,
                const PS::ReallocatableArray<Tsat> & sat,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot){
    etot = ekin = epot = 0.0;
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    PS::F64    m_pla = 1.0;
    PS::F64vec r_pla(0.0);
    const PS::S32 n_ptcl = ptcl.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_ptcl; i++){
        ekin_loc += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
        PS::F64vec dr = ptcl[i].pos - r_pla;
        PS::F64 r_sq = dr*dr;
        PS::F64 r_inv = 1.0 / sqrt(r_sq);
        epot_loc -= ptcl[i].mass * m_pla * r_inv;
    }
    if(PS::Comm::getRank() == 0){
        const PS::S32 n_sat = sat.size();
        for(PS::S32 i=0; i<n_sat; i++){
            ekin_loc += 0.5 * sat[i].mass * sat[i].vel * sat[i].vel;
            PS::F64vec dr = sat[i].pos - r_pla;
            PS::F64 r_sq = dr*dr;
            PS::F64 r_inv = 1.0 / sqrt(r_sq);
            epot_loc -= sat[i].mass * m_pla * r_inv;
        }
    }
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
    etot = ekin + epot;
}
