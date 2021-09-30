#pragma once

#include<iostream>
#include<fstream>
#include<functional>
#include<algorithm>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include<mpi.h>
#endif

namespace ParticleSimulator{

    template<S32 DIM>
    void SetNumberOfDomainMultiDimension(S32 np[], S32 rank[]){
        for(S32 i=0; i<DIMENSION_LIMIT; i++){
            np[i] = 1;
            rank[i] = 1;
        }
        std::vector<S32> npv;
        npv.resize(DIM);
        S32 np_tmp = Comm::getNumberOfProc();
        for(S32 d=DIM, cid=0; cid<DIM-1; d--, cid++){
            S32 tmp = (S32)pow(np_tmp+0.000001, (1.0/d)*1.000001 );
            while(np_tmp%tmp){
                tmp--;
            }
            npv[cid] = tmp;
            np_tmp /= npv[cid];
        }
        npv[DIM-1] = np_tmp;
        S32 rank_tmp = Comm::getRank();
        std::sort(npv.begin(), npv.end(), std::greater<S32>());
        for(S32 i=DIM-1; i>=0; i--){
            np[i] = npv[i];
            rank[i] = rank_tmp % np[i];
            rank_tmp /= np[i];
        }
    }

    class DomainInfo{
    private:
        TimeProfile time_profile_;
        ReallocatableArray<F64vec> pos_sample_tot_;
        ReallocatableArray<F64vec> pos_sample_loc_;

        F64ort * pos_domain_;
        F64ort * pos_domain_temp_;
        
        S32 * n_smp_array_;
        S32 * n_smp_disp_array_;

        F32 coef_ema_;
        S32 target_number_of_sample_particle_;
        S32 number_of_sample_particle_tot_;
        S32 number_of_sample_particle_loc_;
        S32 n_domain_[DIMENSION_LIMIT]; // in 2-dim, n_domain_[2] is always 1.

        F64ort pos_root_domain_;
        
        bool first_call_by_initialize;
        bool first_call_by_decomposeDomain;

        S32 boundary_condition_;
        bool periodic_axis_[DIMENSION_LIMIT]; // in 2-dim, periodic_axis_[2] is always false.


#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // NEW
        MPI_Comm comm_1d_[DIMENSION_LIMIT];
        MPI_Comm comm_sub_[DIMENSION_LIMIT];
        int rank_1d_[DIMENSION_LIMIT];
        int n_proc_1d_[DIMENSION_LIMIT];
        int rank_sub_[DIMENSION_LIMIT];
        int n_proc_sub_[DIMENSION_LIMIT];
#endif
        void calculateBoundaryOfDomainX(const S32 np,
					const F64vec pos_sample[],
					const S32 istart,
					const S32 iend,
					F64 & xlow,
					F64 & xhigh) {
            if(istart == 0) xlow  = pos_root_domain_.low_.x;
	    else xlow  = 0.5 * (pos_sample[istart-1].x + pos_sample[istart].x);
            if(iend == np - 1) xhigh = pos_root_domain_.high_.x;
	    else xhigh = 0.5 * (pos_sample[iend].x + pos_sample[iend+1].x);
        }
        void calculateBoundaryOfDomainY(const S32 np,
					const F64vec pos_sample[],
					const S32 istart,
					const S32 iend,
					F64 & xlow,
					F64 & xhigh) {
            if(istart == 0) xlow  = pos_root_domain_.low_.y;
	    else xlow  = 0.5 * (pos_sample[istart-1].y + pos_sample[istart].y);
            if(iend == np - 1) xhigh = pos_root_domain_.high_.y;
	    else xhigh = 0.5 * (pos_sample[iend].y + pos_sample[iend+1].y);
        }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        void calculateBoundaryOfDomainZ(const S32 np,
					const F64vec pos_sample[],
					const S32 istart,
					const S32 iend,
					F64 & xlow,
					F64 & xhigh) {
            if(istart == 0) xlow  = pos_root_domain_.low_.z;
	    else xlow  = 0.5 * (pos_sample[istart-1].z + pos_sample[istart].z);
            if(iend == np - 1) xhigh = pos_root_domain_.high_.z;
	    else xhigh = 0.5 * (pos_sample[iend].z + pos_sample[iend+1].z);
        }	
#endif
    public:
#ifdef TEST_VARIADIC_TEMPLATE
        void DomainInfoDummyFunc(){}
#endif
        TimeProfile getTimeProfile() const {
            return time_profile_;
        }
        void clearTimeProfile(){
            time_profile_.clear();
        }
        DomainInfo() {
            //std::cerr<<"const"<<std::endl;
            first_call_by_initialize = true;
            first_call_by_decomposeDomain = true;
            pos_sample_tot_.setAllocMode(MemoryAllocMode::Pool);
            pos_sample_loc_.setAllocMode(MemoryAllocMode::Pool);
            periodic_axis_[0] = periodic_axis_[1] = false;
            pos_root_domain_.low_.x  = -LARGE_FLOAT;
            pos_root_domain_.high_.x = LARGE_FLOAT;
            pos_root_domain_.low_.y  = -LARGE_FLOAT;
            pos_root_domain_.high_.y = LARGE_FLOAT;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            periodic_axis_[2] = false;
            pos_root_domain_.low_.z  = -LARGE_FLOAT;
            pos_root_domain_.high_.z = LARGE_FLOAT;
#endif
            boundary_condition_ = BOUNDARY_CONDITION_OPEN;
        }
        ~DomainInfo() {
            //std::cerr<<"dest a"<<std::endl;
            //MemoryPool::dump();
            delete [] n_smp_array_;
            delete [] n_smp_disp_array_;
            delete [] pos_domain_;
            delete [] pos_domain_temp_;
            pos_sample_tot_.freeMem();
            pos_sample_loc_.freeMem();
            //delete [] pos_sample_tot_;
            //delete [] pos_sample_loc_;
        }

        void initialize(const F32 coef_ema = FDPS_DFLT_VAL_COEF_EMA){
            //std::cerr<<"initialize a"<<std::endl;
            if( coef_ema < 0.0 || coef_ema > 1.0){
                PARTICLE_SIMULATOR_PRINT_ERROR("The smoothing factor of an exponential moving average is must between 0 and 1.");
                std::cerr<<"The smoothing factor of an exponential moving average is must between 0 and 1."<<std::endl;
                Abort(-1);
            }
            assert(first_call_by_initialize);

            first_call_by_initialize = false;
            //pos_sample_tot_ = NULL;
            //pos_sample_loc_ = NULL;
	    
            n_smp_array_ = new S32[Comm::getNumberOfProc()];
            
            n_smp_disp_array_ = new S32[Comm::getNumberOfProc() + 1];

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_      = new F64ort[Comm::getNumberOfProc()];
            pos_domain_temp_ = new F64ort[Comm::getNumberOfProc()];
#else
            pos_domain_     = new F64ort[1];
            pos_domain_temp_= new F64ort[1];
#endif

            coef_ema_ = coef_ema;
            target_number_of_sample_particle_ = 0;
            number_of_sample_particle_tot_ = 0;
            number_of_sample_particle_loc_ = 0;

            //S32 rank_tmp[DIMENSION];
            S32 rank_tmp[DIMENSION_LIMIT];
            SetNumberOfDomainMultiDimension<DIMENSION>(n_domain_, rank_tmp);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            // NEW
            int rank_glb = Comm::getRank();
            for(S32 d=DIMENSION-1; d>=0; d--){
                /*
                rank_1d_[d] = rank_glb % n_domain_[d];
                rank_glb /= n_domain_[d];
                std::cerr<<"check 1"<<std::endl;
                std::cerr<<"rank_glb= "<<rank_glb
                         <<" rank_1d_[d]= "<<rank_1d_[d]
                         <<" n_domain_[d]= "<<n_domain_[d]
                         <<std::endl;
                MPI_Comm_split(MPI_COMM_WORLD, rank_1d_[d], rank_glb, comm_sub_+d);
                std::cerr<<"check 2"<<std::endl;
                MPI_Comm_size(comm_1d_[d], n_proc_1d_+d);
                std::cerr<<"check 3"<<std::endl;
                MPI_Comm_rank(comm_sub_[d], rank_sub_+d);
                std::cerr<<"check 4"<<std::endl;
                MPI_Comm_split(MPI_COMM_WORLD, rank_sub_[d], rank_glb, comm_1d_+d);
                std::cerr<<"check 5"<<std::endl;
                MPI_Comm_size(comm_sub_[d], n_proc_sub_+d);
                std::cerr<<"check 6"<<std::endl;
                */
                rank_1d_[d] = rank_glb % n_domain_[d];
                rank_glb /= n_domain_[d];
                MPI_Comm_split(MPI_COMM_WORLD, rank_1d_[d], rank_glb, comm_sub_+d);
                MPI_Comm_rank(comm_sub_[d], rank_sub_+d);
                MPI_Comm_split(MPI_COMM_WORLD, rank_sub_[d], rank_glb, comm_1d_+d);
                MPI_Comm_size(comm_sub_[d], n_proc_sub_+d);
                MPI_Comm_size(comm_1d_[d], n_proc_1d_+d);
                MPI_Comm_rank(comm_1d_[d], rank_1d_+d);                
            }
	    for(S32 d=DIMENSION-1; d>=0; d--){
		Comm::setRankMultiDim(d, rank_tmp[d]);
		Comm::setNumberOfProcMultiDim(d, n_domain_[d]);
	    }
            /*
            if(Comm::getRank()==3){
                for(int d=0; d<3; d++){
                    std::cerr<<"n_proc_1d_[d]= "<<n_proc_1d_[d]
                             <<" rank_1d_[d]= "<<rank_1d_[d]
                             <<" n_proc_sub_[d]= "<<n_proc_sub_[d]
                             <<" rank_sub_[d]= "<<rank_sub_[d]
                             <<std::endl;
                }
            }
            */
#endif
        }

        void setNumberOfDomainMultiDimension(const S32 nx, const S32 ny, const S32 nz=1){
            S32 n_proc = Comm::getNumberOfProc();
            if(n_proc != nx*ny*nz){
                PARTICLE_SIMULATOR_PRINT_ERROR("devided number of domains is not consistent with total processe number");
                Abort(-1);
            }
            n_domain_[0] = nx;
            n_domain_[1] = ny;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            n_domain_[2] = 1;
#else
            n_domain_[2] = nz;
#endif
	    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            int rank_glb = Comm::getRank();
	    // make comunicator
            for(S32 d=DIMENSION-1; d>=0; d--){
                rank_1d_[d] = rank_glb % n_domain_[d];
                rank_glb /= n_domain_[d];
                MPI_Comm_split(MPI_COMM_WORLD, rank_1d_[d], rank_glb, comm_sub_+d);
                MPI_Comm_rank(comm_sub_[d], rank_sub_+d);
                MPI_Comm_split(MPI_COMM_WORLD, rank_sub_[d], rank_glb, comm_1d_+d);
                MPI_Comm_size(comm_sub_[d], n_proc_sub_+d);
                MPI_Comm_size(comm_1d_[d], n_proc_1d_+d);
                MPI_Comm_rank(comm_1d_[d], rank_1d_+d);
            }
	    int rank_tmp = Comm::getRank();
	    for(S32 d=DIMENSION-1; d>=0; d--){
		Comm::setRankMultiDim(d, rank_tmp%n_domain_[d]);
		rank_tmp /= n_domain_[d];
	    }
	    Comm::setNumberOfProcMultiDim(0, nx);
	    Comm::setNumberOfProcMultiDim(1, ny);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    Comm::setNumberOfProcMultiDim(2, ny);
#endif //PARTICLE_SIMULATOR_TWO_DIMENSION
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL

            /*
            if(Comm::getRank()==4){
                for(int d=0; d<3; d++){
                    std::cerr<<"n_proc_1d_[d]= "<<n_proc_1d_[d]
                             <<" rank_1d_[d]= "<<rank_1d_[d]
                             <<" n_proc_sub_[d]= "<<n_proc_sub_[d]
                             <<" rank_sub_[d]= "<<rank_sub_[d]
                             <<std::endl;
                }
            }
            */
        }

        void setDomain(const S32 nx, const S32 ny, const S32 nz=1){
            setNumberOfDomainMultiDimension(nx, ny, nz);
        }

        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys,
                                   const bool clear,
                                   const F32 weight) {
            F64 time_offset = GetWtime();
            if(psys.getFirstCallByDomainInfoCollectSampleParticle()) {
                ReallocatableArray<F64vec> temp_loc(target_number_of_sample_particle_, target_number_of_sample_particle_, MemoryAllocMode::Pool);
                for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
                    temp_loc[i] = pos_sample_loc_[i];
                target_number_of_sample_particle_ += psys.getTargetNumberOfSampleParticle();
                pos_sample_tot_.freeMem();
                pos_sample_loc_.freeMem();
                pos_sample_tot_.resizeNoInitialize(target_number_of_sample_particle_);
                pos_sample_loc_.reserve(target_number_of_sample_particle_);
                for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
                    pos_sample_loc_[i] = temp_loc[i];
                temp_loc.freeMem();
            }
            if(clear) {
                number_of_sample_particle_loc_ = 0;
            }
#if 1
            const auto number_of_sample_particle = psys.getNumberOfSampleParticleLocal(weight);
            pos_sample_loc_.increaseSize(number_of_sample_particle);
            psys.getSampleParticle2(number_of_sample_particle, &pos_sample_loc_[number_of_sample_particle_loc_]);
#else
            // original.
            // if used this, remove pos_sample_loc_.freeMem();
            S32 number_of_sample_particle = 0;
            psys.getSampleParticle(number_of_sample_particle, &pos_sample_loc_[number_of_sample_particle_loc_], weight);
#endif
            number_of_sample_particle_loc_ += number_of_sample_particle;
            //pos_sample_loc_.dump();
            pos_sample_loc_.resizeNoInitialize(number_of_sample_particle_loc_);
            time_profile_.collect_sample_particle += GetWtime() - time_offset;
            return;
        }

        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys,
                                   const bool clear){
            const F32 wgh = psys.getNumberOfParticleLocal();
            collectSampleParticle(psys, clear, wgh);
        }

        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys){
            const F32 wgh = psys.getNumberOfParticleLocal();
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
        }

        void decomposeDomain() {
            F64 time_offset = GetWtime();
            // ****** collect sample particles to process 0. ****** 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    S32 nproc  = Comm::getNumberOfProc();
            S32 myrank = Comm::getRank();
#ifdef __HPC_ACE__
            Comm::allGather(&number_of_sample_particle_loc_, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            pos_sample_tot_.resizeNoInitialize(n_smp_disp_array_[nproc]);
            Comm::allGatherV(&pos_sample_loc_[0], &number_of_sample_particle_loc_[0], pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#else //__HPC_ACE__

            //Comm::gather(&number_of_sample_particle_loc_, 1, n_smp_array_);
            Comm::allGather(&number_of_sample_particle_loc_, 1, n_smp_array_);
            
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            pos_sample_tot_.resizeNoInitialize(n_smp_disp_array_[nproc]);
            Comm::allGatherV(&pos_sample_loc_[0], number_of_sample_particle_loc_, &pos_sample_tot_[0], n_smp_array_, n_smp_disp_array_);
            
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#endif //__HPC_ACE__


            // ****************************************************
            // *** decompose domain *******************************
            if(myrank == 0) {
                S32 * istart = new S32[nproc];
                S32 * iend   = new S32[nproc];
                // --- x direction --------------------------
		//std::sort(pos_sample_tot_, pos_sample_tot_+number_of_sample_particle_tot_, Cmpvec(&F64vec::x));
                std::sort(&pos_sample_tot_[0], &pos_sample_tot_[number_of_sample_particle_tot_], Cmpvec(&F64vec::x));
                for(S32 i = 0; i < nproc; i++) {
                    istart[i] = ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
                    if(i > 0)
                        iend[i-1] = istart[i] - 1;
                }
                iend[nproc-1] = number_of_sample_particle_tot_ - 1;
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    F64 x0 = 0.0;
		    F64 x1 = 0.0;
                    //calculateBoundaryOfDomain(number_of_sample_particle_tot_, pos_sample_tot_, 0, istart[ix0], iend[ix1-1], x0, x1);
		    calculateBoundaryOfDomainX(number_of_sample_particle_tot_, &pos_sample_tot_[0], istart[ix0], iend[ix1-1], x0, x1);
                    for(S32 i = ix0; i < ix1; i++) {
                        pos_domain_temp_[i].low_.x  = x0;
                        pos_domain_temp_[i].high_.x = x1;
                    }
                }
                // ------------------------------------------
                // --- y direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    //sortCoordinateOfSampleParticle(pos_sample_tot_, istart[ix0], iend[ix1-1], 1);
		    //std::sort(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), Cmpvec(&F64vec::y));
                    std::sort(&pos_sample_tot_[istart[ix0]], &pos_sample_tot_[(iend[ix1-1]+1)], Cmpvec(&F64vec::y));
/*
                    for(S32 i=istart[ix0]; i<iend[ix1-1]+1; i++){
                        if(pos_sample_tot_[i+1].y < pos_sample_tot_[i].y){
                            std::cout<<"y sort is wrong: i="<<i<<std::endl;
                            std::cout<<"pos_sample_tot_[i+1]="<<pos_sample_tot_[i+1]<<std::endl;
                            std::cout<<"pos_sample_tot_[i]="<<pos_sample_tot_[i]<<std::endl;
                        }
                    }
*/
                    S32 number_of_sample_particle_tot_y = iend[ix1-1] - istart[ix0] + 1;
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        //F64 y0, y1;
			F64 y0 = 0.0;
			F64 y1 = 0.0;
                        //std::cout<<"ix="<<ix<<" ix0="<<ix0<<" ix1="<<ix1<<" iy="<<iy<<" iy0="<<iy0<<" iy1="<<iy1<<" y0="<<y0<<" y1="<<y1<<std::endl;
                        //calculateBoundaryOfDomain(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], 1, istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
			//calculateBoundaryOfDomainY(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
                        calculateBoundaryOfDomainY(number_of_sample_particle_tot_y, &pos_sample_tot_[istart[ix0]], istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
                        for(S32 i = iy0; i < iy1; i++) {
                            pos_domain_temp_[i].low_.y  = y0;
                            pos_domain_temp_[i].high_.y = y1;
                        }
                    }
                }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- z direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        //sortCoordinateOfSampleParticle(pos_sample_tot_, istart[iy0], iend[iy1-1], 2);
			//std::sort(pos_sample_tot_+istart[iy0], pos_sample_tot_+(iend[iy1-1]+1), Cmpvec(&F64vec::z));
                        std::sort(&pos_sample_tot_[istart[iy0]], &pos_sample_tot_[(iend[iy1-1]+1)], Cmpvec(&F64vec::z));
                        S32 number_of_sample_particle_tot_z = iend[iy1-1] - istart[iy0] + 1;
                        for(S32 iz = 0; iz < n_domain_[2]; iz++) {
                            S32 iz0 = iy0 + iz;
                            //F64 z0, z1;
			    F64 z0 = 0.0;
			    F64 z1 = 0.0;						    
                            //calculateBoundaryOfDomain(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], 2, istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
			    //calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                            calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, &pos_sample_tot_[istart[iy0]], istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                            pos_domain_temp_[iz0].low_.z  = z0;
                            pos_domain_temp_[iz0].high_.z = z1;
                        }
                    }
                }
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- process first ------------------------
                if(first_call_by_decomposeDomain) {
                    //first_call_by_decomposeDomain = false;
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                        pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                    }
                } else {
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                        pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                        
                    }
                }
                // ------------------------------------------
                delete [] istart;
                delete [] iend;
            }
            // ****************************************************
            // *** broad cast pos_domain_ *************************
            MPI_Bcast(pos_domain_, nproc, GetDataType<F64ort>(), 0, MPI_COMM_WORLD);
            if(first_call_by_decomposeDomain) {
                first_call_by_decomposeDomain = false;            
                MPI_Bcast(&first_call_by_decomposeDomain, 1, GetDataType<bool>(), 0, MPI_COMM_WORLD);
            }
            //std::cout<<"end of bcast: "<<"time: "<<GetWtime() - Tbegin<<std::endl;
            //Comm::broadcast(pos_domain_, nproc);
            // ****************************************************
#else       // PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_[0] = pos_root_domain_;
#endif     // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"pos_root_domain_="<<pos_root_domain_<<std::endl;
            std::cout<<"pos_domain_[Comm::getRank()]="<<pos_domain_[Comm::getRank()]<<std::endl;
#endif

            pos_sample_tot_.freeMem();
            pos_sample_loc_.freeMem();
            
            time_profile_.decompose_domain += GetWtime() - time_offset;
        }

//#if __cplusplus <= 201112L
#ifdef TEST_VARIADIC_TEMPLATE

    #if 0 // variadic template version
        F32 sample_weight_;
        bool flag_sample_weight_;
        bool flag_first_sample_;
        // variadic template version
        template<class Tail>
        void checkWeightLast(const Tail & tail, std::true_type){
            // sample_weight is defined
            sample_weight_ = tail;
            flag_sample_weight_ = true;
        }
        template<class Tail>
        void checkWeightLast(Tail & tail, std::false_type){
            // sample_weight is not defined
            flag_sample_weight_ = false;
        }
        template<class Tail>
        void checkWeight(const Tail & tail){
            checkWeightLast(tail, typename std::is_convertible<Tail, F32>::type() );
        }
        template<class Head, class... Tail>
        void checkWeight(const Head & head, const Tail&... tail){
            checkWeight(tail...);
        }
        template<class Tail>
        void collectSampleParticleLast(const Tail & tail, std::true_type){
            // last argument is not a particle_system class
            // do nothing
            //std::cerr<<"weight is defined"<<std::endl;
        }
        template<class Tail>
        void collectSampleParticleLast(Tail & tail, std::false_type){
            // last argument is a particle_system class
            //std::cerr<<"no weight"<<std::endl;
            collectSampleParticle(tail, false);
        }
        template<class Tail>
        void collectSampleParticleRecursive(Tail & tail){
            collectSampleParticleLast(tail, typename std::is_convertible<Tail, F32>::type() );
        }
        template<class Head, class ... Tail>
        void collectSampleParticleRecursive(Head & head, Tail& ... tail){
            bool flag_clear = false;
            if(!flag_sample_weight_){sample_weight_ = head.getNumberOfParticleLocal();}
            if(flag_first_sample_){
                flag_clear = true;
                flag_first_sample_ = false;
            }
            collectSampleParticle(head, flag_clear, sample_weight_);
            collectSampleParticleRecursive(tail ...);
        }
        template<class ... Args>
        void decomposeDomainAll(Args & ... args){
            checkWeight(args ...);
            flag_first_sample_ = true;
            collectSampleParticleRecursive(args...);
            decomposeDomain();
        }
    #endif // variadic template version        
        // tuple version
    #if 0
        // use class
        template<int N, class Ttpl, class Head, class ... Tail>
        class CollectSampleParticleTplRecursive{
        public:
            void operator () (Ttpl & sys_tpl, F32 wgh, bool & flag_first_sample, DomainInfo * dinfo){
                bool flag_clear = false;
                if(flag_first_sample){
                    flag_clear = true;
                    flag_first_sample = false;
                }
                dinfo->collectSampleParticle(*(std::get<N>(sys_tpl)), flag_clear, wgh);
                CollectSampleParticleTplRecursive<N+1, Ttpl, Tail...>()(sys_tpl, wgh, flag_first_sample, dinfo);
            }
            void operator () (Ttpl & sys_tpl, bool & flag_first_sample, DomainInfo * dinfo){
                bool flag_clear = false;
                if(flag_first_sample){
                    flag_clear = true;
                    flag_first_sample = false;
                }
                F32 wgh = std::get<N>(sys_tpl)->getNumberOfParticleLocal();
                dinfo->collectSampleParticle(*(std::get<N>(sys_tpl)), flag_clear, wgh);
                CollectSampleParticleTplRecursive<N+1, Ttpl, Tail...>()(sys_tpl, flag_first_sample, dinfo);
            }
        };
        template<int N, class Ttpl, class Tail>
        class CollectSampleParticleTplRecursive<N, Ttpl, Tail>{
        public:
            void operator () (Ttpl & sys_tpl, F32 wgh, bool & flag_first_sample, DomainInfo * dinfo){}
            void operator () (Ttpl & sys_tpl, bool & flag_first_sample, DomainInfo * dinfo){}
        };
        template<class ... Args>
        void decomposeDomainAll(std::tuple<Args...> & sys_tpl,
                                F32 wgh){
            flag_first_sample_ = true;
            CollectSampleParticleTplRecursive<0, std::tuple<Args...>, Args... >()(sys_tpl, wgh, flag_first_sample_, this);
            decomposeDomain();
        }
        template<class ... Args>
        void decomposeDomainAll(std::tuple<Args...> & sys_tpl){
            flag_first_sample_ = true;
            CollectSampleParticleTplRecursive<0, std::tuple<Args...>, Args... >()(sys_tpl, flag_first_sample_, this);
            decomposeDomain();
        }
    #else
        template<int N, class Ttpl>
        void collectSampleParticleTplRecursive(Ttpl & sys_tpl, F32 wgh){}
        template<int N, class Ttpl, class Head, class ... Tail>
        void collectSampleParticleTplRecursive(Ttpl & sys_tpl, F32 wgh){
            bool flag_clear = false;
            if(N==0){ flag_clear = true; }
            collectSampleParticle(std::get<N>(sys_tpl), flag_clear, wgh);
            collectSampleParticleTplRecursive<N+1, Ttpl, Tail...>(sys_tpl, wgh);
        }
        template<int N, class Ttpl>
        void collectSampleParticleTplRecursive(Ttpl & sys_tpl){}
        template<int N, class Ttpl, class Head, class ... Tail>
        void collectSampleParticleTplRecursive(Ttpl & sys_tpl){
            bool flag_clear = false;
            if(N==0){ flag_clear = true; }
            F32 wgh = std::get<N>(sys_tpl).getNumberOfParticleLocal();
            collectSampleParticle(std::get<N>(sys_tpl), flag_clear, wgh);
            collectSampleParticleTplRecursive<N+1, Ttpl, Tail...>(sys_tpl);
        }
        template<class ... Args>
        void decomposeDomainAll(std::tuple<Args...> & sys_tpl,
                                F32 wgh){
            collectSampleParticleTplRecursive<0, std::tuple<Args...>, Args... >(sys_tpl, wgh);
            decomposeDomain();
        }
        template<class ... Args>
        void decomposeDomainAll(std::tuple<Args...> & sys_tpl){
            collectSampleParticleTplRecursive<0, std::tuple<Args...>, Args... >(sys_tpl);
            decomposeDomain();
        }
    #endif
        /*
        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys,
                                const F32 wgh){
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }
        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys){
            const F32 wgh = psys.getNumberOfParticleLocal();
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }
        */
#endif
        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys,
                                const F32 wgh){
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }

        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys){
            const F32 wgh = psys.getNumberOfParticleLocal();
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }
        
        void getRootDomain(FILE *fp) {
            fprintf(fp, "%+e %+e %+e\n",
                pos_root_domain_.low_[0],
                pos_root_domain_.low_[1],
                pos_root_domain_.low_[2]);
            fprintf(fp, "%+e %+e %+e\n",
                pos_root_domain_.high_[0],
                pos_root_domain_.high_[1],
                pos_root_domain_.high_[2]);

            return;
        }


#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        void getSampleParticleLocal(FILE *fp) {
            for(S32 i = 0; i < number_of_sample_particle_loc_; i++) {
                fprintf(fp, "%+e %+e %+e\n", pos_sample_loc_[i].x, pos_sample_loc_[i].y, pos_sample_loc_[i].z);
            }            
            return;
        }

        void getSampleParticleTotal(FILE *fp) {
            for(S32 i = 0; i < number_of_sample_particle_tot_; i++) {
                fprintf(fp, "%+e %+e %+e\n", pos_sample_tot_[i].x, pos_sample_tot_[i].y, pos_sample_tot_[i].z);
            }            
            return;
        }
#endif
	
        void getPosDomainTotal(FILE *fp) {
            S32 nproc = Comm::getNumberOfProc();
            for(S32 i = 0; i < nproc; i++) {
                for(S32 k = 0; k < DIMENSION; k++)
                    fprintf(fp, "%+e ", pos_domain_[i].low_[k]);
                for(S32 k = 0; k < DIMENSION; k++)
                    fprintf(fp, "%+e ", pos_domain_[i].high_[k]);
                fprintf(fp, "\n");
            }
        }

        S32 * getPointerOfNDomain(){return n_domain_;};

        /* AT_DEBUG
        F32ort * getPointerOfPosDomain(){return pos_domain_;};
        */
        F64ort * getPointerOfPosDomain(){return pos_domain_;};

        // A. Tanikawa need this method for Particle Mesh...
        S32 getNDomain(const S32 dim) const {return n_domain_[dim];}
        S32 getRank1D(const S32 rank_glb, const S32 dim) const {
            S32 rank_1d;
            if(dim==0)
                rank_1d = rank_glb / (n_domain_[1]*n_domain_[2]);
            else if(dim==1)
                rank_1d = (rank_glb / n_domain_[2]) % n_domain_[1];
            else if(dim==2)
                rank_1d = rank_glb % n_domain_[2];
            return rank_1d;
        }
        

        /* AT_DEBUG
        // for DEBUG
        F32vec & getPosSample(const S32 id=0){return pos_sample_tot_[id];}
        // for DEBUG // A. Tanikawa would not like to delete this method...
        F32ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
        // for DEBUG
        void setPosDomain(const S32 id, const F32ort & pos){ pos_domain_[id] = pos;}
        */
        // for DEBUG
        F64vec & getPosSample(const S32 id=0){return pos_sample_tot_[id];}
        // for DEBUG // A. Tanikawa would not like to delete this method...
        F64ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
        F64ort getPosDomain(const F64ort & box, const S32 id=0) const {
            F64ort ret = pos_domain_[id];
            if(periodic_axis_[0] == false){
                if(getRank1D(id, 0)==0){
                    ret.low_.x = box.low_.x;
                }
                if(getRank1D(id, 0)==n_domain_[0]-1){
                    ret.high_.x = box.high_.x;
                }
            }
            if(periodic_axis_[1] == false){
                if(getRank1D(id, 1)==0){
                    ret.low_.y = box.low_.y;
                }
                if(getRank1D(id, 1)==n_domain_[1]-1){
                    ret.high_.y = box.high_.y;
                }
            }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            if(periodic_axis_[2] == false){
                if(getRank1D(id, 2)==0){
                    ret.low_.z = box.low_.z;
                }
                if(getRank1D(id, 2)==n_domain_[2]-1){
                    ret.high_.z = box.high_.z;
                }
            }
#endif
            return ret;
        }
        // for DEBUG
        void setPosDomain(const S32 id, const F64ort & pos){ pos_domain_[id] = pos;}

        void setBoundaryCondition(enum BOUNDARY_CONDITION bc){
            /*
            std::cerr<<"bc= "<<bc
                     <<" BOUNDARY_CONDITION_PERIODIC_XY= "<<BOUNDARY_CONDITION_PERIODIC_XY
                     <<std::endl;
            */
            boundary_condition_ = bc;
            if(DIMENSION == 2 && 
               (bc == BOUNDARY_CONDITION_PERIODIC_XYZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_XZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_YZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_Z ) ){
                throw "PS_ERROR: in setBoundaryCondition(enum BOUNDARY_CONDITION) \n boundary condition is incompatible with DIMENSION";
            }
            if(bc == BOUNDARY_CONDITION_PERIODIC_X) periodic_axis_[0] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_Y) periodic_axis_[1] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_Z) periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XY) periodic_axis_[0] = periodic_axis_[1] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XZ) periodic_axis_[0] = periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_YZ) periodic_axis_[1] = periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XYZ) periodic_axis_[0] = periodic_axis_[1] = periodic_axis_[2] = true;
        }

        S32 getBoundaryCondition() const { return boundary_condition_; }

        /* AT_DEBUG
        void setPosRootDomain(const F32vec & low, const F32vec & high){
        */
        void setPosRootDomain(const F64vec & low, const F64vec & high){
            for(S32 k=0; k<DIMENSION; k++){
                if(low[k] > high[k]){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The coodinate of the root domain is inconsistent.");
                    std::cerr<<"The coordinate of the low vertex of the rood domain="<<low<<std::endl;
                    std::cerr<<"The coordinate of the high vertex of the rood domain="<<high<<std::endl;
                    Abort(-1);
                }
            }

            for(S32 i=0; i<DIMENSION; i++){
                if( periodic_axis_[i] == false ) continue;
                if(low[i] < high[i]){
                    pos_root_domain_.low_[i] = low[i];
                    pos_root_domain_.high_[i] = high[i];
                }
            }
        }

        const F64ort getPosRootDomain() const { return pos_root_domain_; }

        void getPeriodicAxis(bool pa[]) const {
            for(S32 i=0; i<DIMENSION; i++) pa[i] = periodic_axis_[i];
        }

        template<class Tpsys>
        bool checkCollectSampleParticleSubset(Tpsys & psys);
        template<class Tpsys>
        bool checkCollectSampleParticleAverage(Tpsys & psys);
        template<class Tpsys>
        bool checkDecomposeDomain(Tpsys & psys);

        S32 whereToGoOpen(const F64 x,
                          const F64 coord_sep[], // size is n_proc+1
                          const S32 my_rank,
                          const S32 shift = 1){
            if(coord_sep[my_rank] <= x && x < coord_sep[my_rank+1]){
                return my_rank;
            }
            auto rank = my_rank;
            if(x < coord_sep[rank]){
                do{
                    rank -= shift;
                }while(x < coord_sep[rank]);
            }
            else if(x >= coord_sep[rank+1]){
                do{
                    rank += shift;
                }while(x >= coord_sep[rank+1]);
            }
            return rank;
        }

        S32 whereToGoPeriodic(const F64 x,
                              const F64 coord_sep[], // size is n_proc+1
                              const S32 my_rank,
                              const F64 len_peri,
                              const S32 n_proc,
                              const S32 shift = 1){
            if(coord_sep[my_rank] <= x && x < coord_sep[my_rank+1]){
                return my_rank;
            }
            auto rank = my_rank;
            auto x_tmp = (x < coord_sep[rank]) ?
                (coord_sep[rank]-x < x+len_peri-coord_sep[rank+1] ? x : x+len_peri) :
                (x-coord_sep[rank+1]< coord_sep[rank]-(x-len_peri) ? x : x-len_peri);
            if(x_tmp < coord_sep[rank]){
                if(x_tmp != x) rank = n_proc-1;
                else rank -= shift;
                while(x < coord_sep[rank]){
                    rank -= shift;
                };
            }
            else if(x_tmp >= coord_sep[rank+1]){
                if(x_tmp != x) rank = 0;
                else rank += shift;
                while(x >= coord_sep[rank+1]){
                    rank += shift;
                };
            }
            return rank;
        }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        
        void determineRoughPartition(ReallocatableArray<F64> & coord, // in. not const, becuase the order is temporally changed in Randomsampling.
                                     ReallocatableArray<F64> & coord_partition, // out
                                     const S32 smp_frac, // in
                                     PS_MPI_COMM comm, // in
                                     bool print = false){
            const auto my_rank = Comm::getRank(comm);
            const auto n_proc  = Comm::getNumberOfProc(comm);
            auto n_coord     = coord.size();
            auto n_coord_smp = ((n_coord/smp_frac) >= 1) ? n_coord/smp_frac : 1;
            if(n_coord <= 0) n_coord_smp = 0;
            ReallocatableArray<F64> coord_smp(n_coord_smp, n_coord_smp, MemoryAllocMode::Pool);
            RandomSampling(&coord[0], &coord_smp[0], n_coord, n_coord_smp,
                           [](F64 & src, F64 & dst){dst = src;});
            ReallocatableArray<S32> n_coord_smp_ar(n_proc, n_proc, MemoryAllocMode::Pool);
            ReallocatableArray<S32> n_disp_coord_smp_ar(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
#if 0
            Comm::allGather(&n_coord_smp, 1, &n_coord_smp_ar[0], comm);
            n_disp_coord_smp_ar[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_coord_smp_ar[i+1] = n_disp_coord_smp_ar[i] + n_coord_smp_ar[i];
            }
            const auto n_coord_smp_glb = n_disp_coord_smp_ar[n_proc];
#else
            Comm::gather(&n_coord_smp, 1, &n_coord_smp_ar[0], comm);
            if(my_rank==0){
                n_disp_coord_smp_ar[0] = 0;
                for(S32 i=0; i<n_proc; i++){
                    n_disp_coord_smp_ar[i+1] = n_disp_coord_smp_ar[i] + n_coord_smp_ar[i];
                }
            }
            const auto n_coord_smp_glb = (my_rank == 0) ? n_disp_coord_smp_ar[n_proc] : 0;
#endif
            ReallocatableArray<F64> coord_smp_glb(n_coord_smp_glb, n_coord_smp_glb, MemoryAllocMode::Pool);
            Comm::gatherV(&coord_smp[0], n_coord_smp, &coord_smp_glb[0],
                          &n_coord_smp_ar[0], &n_disp_coord_smp_ar[0], comm);
            if(my_rank==0){
                std::sort(coord_smp_glb.frontP(), coord_smp_glb.endP());
                if(Comm::getRank()==0){
                    //std::cerr<<"coord_smp_glb.size()= "<<coord_smp_glb.size()<<std::endl;
                    //for(S32 i=0; i<coord_smp_glb.size(); i++){
                    //    std::cerr<<"i= "<<i<<" coord_smp_glb[i]= "<<coord_smp_glb[i]<<std::endl;
                    //}
                }
                for(S32 i=1; i<n_proc; i++){
                    const auto id = (n_coord_smp_glb/n_proc)*i;
                    if(id == 0 || id == n_coord_smp_glb-1)
                        coord_partition[i] = coord_smp_glb[id];
                    else
                        coord_partition[i] = (coord_smp_glb[id]+coord_smp_glb[id+1])*0.5 ;
                }
            }
            Comm::broadcast(&coord_partition[0], n_proc+1, 0, comm);
        }

        template<typename T, typename Tgetter>
        void exchangeSample(ReallocatableArray<T> & smp, // in and out
                            ReallocatableArray<F64> & coord_partition, // in
                            const S32 n_proc, // in
                            const S32 my_rank, // in
                            const F64 len_peri, // in
                            const bool flag_periodic, // in
                            PS_MPI_COMM comm, // in
                            Tgetter getter
                            ){
            const auto n_smp = smp.size();
            ReallocatableArray<S32> n_send(n_proc, n_proc, MemoryAllocMode::Pool);
            ReallocatableArray<S32> rank_send(n_smp, n_smp, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Request> req_recv(n_proc, n_proc, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Status> status(n_proc, n_proc, MemoryAllocMode::Pool);
            for(S32 i=0; i<n_proc; i++)
                n_send[i] = 0;
            S32 dnp_max = 0;
            S32 i_loc = n_smp-1;
            if(!flag_periodic){
                for(S32 i=n_smp-1; i>=0; i--){
                    const auto rank = whereToGoOpen(getter(smp[i]), &coord_partition[0], my_rank);
                    rank_send[i] = rank;
                    n_send[rank]++;
                    if(rank != my_rank){
                        std::swap(smp[i], smp[i_loc]);
                        rank_send[i_loc] = rank;
                        i_loc--;
                    }
                    const auto dnp = std::abs(rank - my_rank);
                    dnp_max = std::max(dnp_max, dnp);
                }
            }
            else{
                for(S32 i=n_smp-1; i>=0; i--){
                    const auto rank = whereToGoPeriodic(getter(smp[i]), &coord_partition[0], my_rank, len_peri, n_proc);
                    rank_send[i] = rank;
                    n_send[rank]++;
                    if(rank != my_rank){
                        std::swap(smp[i], smp[i_loc]);
                        rank_send[i_loc] = rank;
                        i_loc--;
                    }
                    const auto dnp = std::abs(rank - my_rank);
                    dnp_max = std::max(dnp_max, dnp);
                }
            }
            const auto n_remaine = i_loc+1;

            //if(my_rank==0){
            //    std::cerr<<"n_remaine= "<<n_remaine<<std::endl;
            //    for(S32 i=0; i<smp.size(); i++){
            //        std::cerr<<"A) i= "<<i<<"rank_send[i]= "<<rank_send[i]<<" smp[i]= "<<getter(smp[i])<<std::endl;
            //    }
            //}
            
            //std::cerr<<"my_rank= "<<my_rank
            //         <<" n_smp= "<<n_smp
            //         <<" n_remaine= "<<n_remaine
            //         <<" dnp_max= "<<dnp_max
            //         <<std::endl;
            //for(S32 i=0; i<n_proc; i++){
            //    std::cerr<<"n_send[i]= "<<n_send[i]<<std::endl;
            //}
            ReallocatableArray<S32> n_disp_send(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_send[i+1] = n_disp_send[i] + n_send[i];
                n_send[i] = 0; // reuse
            }
            ReallocatableArray<T> smp_send(n_disp_send[n_proc], n_disp_send[n_proc], MemoryAllocMode::Pool);
            for(S32 i=n_remaine; i<n_smp; i++){
                const auto rank = rank_send[i];
                const auto loc  = n_send[rank] + n_disp_send[rank];
                smp_send[loc] = smp[i];
                n_send[rank]++;
            }
            const auto dnp_max_glb = Comm::getMaxValue(dnp_max, comm);
            //std::cerr<<"my_rank= "<<my_rank
            //         <<" dnp_max_glb= "<<dnp_max_glb
            //         <<std::endl;
            S32 dnp_head = -dnp_max_glb;
            S32 dnp_tail = dnp_max_glb;
            if(dnp_max_glb*2+1 > n_proc){
                dnp_head = -n_proc / 2;
                dnp_tail = (n_proc%2==0) ? (n_proc/2 - 1) : n_proc/2;
            }
            //std::cerr<<"Comm::getRank()= "<<Comm::getRank()
            //         <<" dnp_head= "<<dnp_head
            //         <<" dnp_tail= "<<dnp_tail
            //         <<std::endl;
            ReallocatableArray<S32> n_recv(n_proc, n_proc, MemoryAllocMode::Pool);
            for(S32 i=0; i<n_proc; i++){ n_recv[i] = 0; }
            if(1){
                S32 n_cnt = 0;
                for(S32 i=dnp_head; i<=dnp_tail; i++){
                    auto rank = (my_rank+n_proc+i) % n_proc;
                    if(my_rank == rank) continue;
                    MPI_Irecv(&n_recv[rank], 1, GetDataType<S32>(), rank, 0, comm, &req_recv[n_cnt]);
                    n_cnt++;
                }
                for(S32 i=dnp_head; i<=dnp_tail; i++){
                    auto rank = (my_rank+n_proc+i) % n_proc;
                    if(my_rank == rank) continue;
                    MPI_Send(&n_send[rank], 1, GetDataType<S32>(), rank, 0, comm);
                }
                MPI_Waitall(n_cnt, &req_recv[0], &status[0]);
            }
            else{
                // use all to all
                Comm::allToAll(&n_send[0], 1, &n_recv[0], comm);
            }
            
            //std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
            //for(S32 i=0; i<n_proc; i++){
            //    std::cerr<<"i= "<<i<<" n_send[i]= "<<n_send[i]
            //             <<" n_recv[i]= "<<n_recv[i]
            //             <<" n_remaine= "<<n_remaine<<std::endl;
            //}


            
            ReallocatableArray<S32> n_disp_recv(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++)
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            smp.resizeNoInitialize(n_disp_recv[n_proc] + n_remaine);
            //if(2*dnp_max_glb < n_proc){
            if(1){
                S32 n_cnt = 0;
                for(S32 i=dnp_head; i<=dnp_tail; i++){
                    const auto rank = (my_rank+n_proc+i) % n_proc;
                    const auto adr = n_disp_recv[rank] + n_remaine;
                    MPI_Irecv(&smp[adr], n_recv[rank], GetDataType<T>(), rank, 1, comm, &req_recv[n_cnt]);
                    n_cnt++;
                }
                for(S32 i=dnp_head; i<=dnp_tail; i++){
                    const auto rank = (my_rank+n_proc+i) % n_proc;
                    const auto adr  = n_disp_send[rank];
                    MPI_Send(&smp_send[adr], n_send[rank], GetDataType<T>(), rank, 1, comm);
                }                
                MPI_Waitall(n_cnt, &req_recv[0], &status[0]);
            }
            else{
                // alltoallv
                Comm::allToAllV(&smp_send[0], &n_send[0], &n_disp_send[0], &smp[n_remaine], &n_recv[0], &n_disp_recv[0], comm);
            }
            //std::cerr<<"my_rank= "<<my_rank<<" smp.size()= "<<smp.size()<<std::endl;
            //for(S32 i=0; i<smp.size(); i++){
            //    std::cerr<<"i= "<<i<<" smp[i]= "<<getter(smp[i])<<std::endl;
            //}
        }

        
        void determinePartition(ReallocatableArray<F64> & coord, // in
                                ReallocatableArray<F64> & coord_partition_tmp, // in
                                ReallocatableArray<F64> & coord_partition, // out
                                PS_MPI_COMM comm, // in
                                bool print = false){
            const auto my_rank = Comm::getRank(comm);
            const auto n_proc  = Comm::getNumberOfProc(comm);            
            const auto n_smp_loc = coord.size();
            ReallocatableArray<S32> n_smp(n_proc, n_proc, MemoryAllocMode::Pool);
            Comm::allGather(&n_smp_loc, 1, &n_smp[0], comm);
            //std::cerr<<"my_rank= "<<my_rank<<std::endl;
            //for(S32 i=0; i<n_proc; i++){
            //    std::cerr<<"i= "<<i<<" n_smp[i]= "<<n_smp[i]<<std::endl;
            //}
            
            ReallocatableArray<S32> n_disp_smp(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
            n_disp_smp[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_smp[i+1] = n_disp_smp[i] + n_smp[i];
            }
            const auto n_tot = n_disp_smp[n_proc];
            const auto n_ave   = n_tot / n_proc;
            const auto n_amari = n_tot % n_proc;
            //std::cerr<<"n_tot= "<<n_tot
            //         <<" n_ave= "<<n_ave
            //         <<" n_amari= "<<n_amari
            //         <<std::endl;
            ReallocatableArray<F64> coord_partition_send(n_proc-1, 0, MemoryAllocMode::Pool); // whitout eadge
            ReallocatableArray<S32> n_partition_recv(n_proc, n_proc, MemoryAllocMode::Pool);
            for(S32 i=0; i<n_proc; i++)
                n_partition_recv[i] = 0;
            S32 i_domain = 0;
            S32 offset = n_disp_smp[my_rank];
            for(S32 i=1; i<n_proc; i++){
                const auto index_partition = n_ave*i + ((i<n_amari) ? 1 : 0);
                while( !(n_disp_smp[i_domain] <= index_partition && index_partition < n_disp_smp[i_domain+1]) ){
                    i_domain++;
                }
                //if(1){
                //    if(Comm::getRank() == 2){
                //        std::cerr<<"i= "<<i<<" index_partition= "<<index_partition
                //                 <<" i_domain= "<<i_domain
                //                 <<" my_rank= "<<my_rank
                //                 <<" head= "<<n_disp_smp[i_domain]
                //                 <<" tail= "<<n_disp_smp[i_domain+1]-1
                //                 <<std::endl;
                //    }
                //}
                if(my_rank == i_domain){
                    const auto head = n_disp_smp[i_domain];
                    const auto tail = n_disp_smp[i_domain+1]-1;
                    if(index_partition==head){
                        //coord_partition_send[n_partition_recv[i_domain]] = coord_partition_tmp[i-1];
                        coord_partition_send[n_partition_recv[i_domain]] = coord_partition_tmp[i_domain];
                    }
                    else if(index_partition==tail){
                        //coord_partition_send[n_partition_recv[i_domain]] = coord_partition_tmp[i];
                        coord_partition_send[n_partition_recv[i_domain]] = coord_partition_tmp[i_domain+1];
                    }
                    else{
                        //if(Comm::getRank() == 2){
                        //    std::cerr<<"n_partition_recv[i_domain]= "<<n_partition_recv[i_domain]
                        //             <<" index_partition-offset= "<<index_partition-offset
                        //             <<" coord[index_partition-offset]= "<<coord[index_partition-offset]
                        //             <<" coord[index_partition-offset+1]= "<<coord[index_partition-offset+1]
                        //             <<std::endl;
                        //}
                        coord_partition_send[n_partition_recv[i_domain]] = (coord[index_partition-offset]+coord[index_partition-offset+1])*0.5;
                    }
                }
                n_partition_recv[i_domain]++;
            }
            auto n_partition_send = n_partition_recv[my_rank];
            ReallocatableArray<S32> n_disp_partition_recv(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
            n_disp_partition_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_partition_recv[i+1] = n_disp_partition_recv[i] + n_partition_recv[i];
            }
            
            Comm::gatherV(&coord_partition_send[0], n_partition_send, &coord_partition[1], &n_partition_recv[0], &n_disp_partition_recv[0], comm);
            // NOTE: in the 3rd argument, coord_partition[1] means we do not calculate eadge coordinate.
        }

        template<typename T0, typename T1, typename Tgetter>
        void gatherSampleRandomly(ReallocatableArray<T0> & smp, // in
                                  ReallocatableArray<T1> & smp_gat, // out
                                  const S32 n_smp_gat, // in 
                                  PS_MPI_COMM comm, // in
                                  Tgetter getter
                                  ){
            const auto my_rank = Comm::getRank(comm);
            const auto n_proc  = Comm::getNumberOfProc(comm);
            const auto n_smp = smp.size();
            ReallocatableArray<S32> n_smp_ar(n_proc,   n_proc, MemoryAllocMode::Pool);
            ReallocatableArray<S32> n_disp_smp_ar(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
            Comm::allGather(&n_smp, 1, &n_smp_ar[0], comm);
            n_disp_smp_ar[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_smp_ar[i+1] = n_disp_smp_ar[i] + n_smp_ar[i];
            }
            const auto n_smp_glb = n_disp_smp_ar[n_proc];
            std::mt19937 mt(1);
            std::uniform_int_distribution<S32> dist(0, n_smp_glb-1);
            ReallocatableArray<S32> index(n_smp_gat, n_smp_gat, MemoryAllocMode::Pool);
            for(S32 i=0; i<n_smp_gat; i++){
                S32 j = dist(mt);
                index[i] = j;
            }
            S32 n_smp_send = 0;
            for(S32 i=0; i<n_smp_gat; i++){
                if(n_disp_smp_ar[my_rank] <= index[i] && n_disp_smp_ar[my_rank+1] > index[i]){
                    n_smp_send++;
                }
            }
            n_smp_send = std::min(n_smp_send, n_smp);
            const auto n_smp_send_glb = Comm::getSum(n_smp_send, comm);
            assert(n_smp_send_glb == n_smp_gat);
            ReallocatableArray<T1> smp_send(n_smp_send, n_smp_send, MemoryAllocMode::Pool);
            RandomSampling(&smp[0], &smp_send[0], n_smp, n_smp_send, getter);

            // reuse n_smp_ar and n_disp_smp_ar
            Comm::allGather(&n_smp_send, 1, &n_smp_ar[0], comm);
            n_disp_smp_ar[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_smp_ar[i+1] = n_disp_smp_ar[i] + n_smp_ar[i];
            }
            Comm::gatherV(&smp_send[0], n_smp_send, &smp_gat[0], &n_smp_ar[0], &n_disp_smp_ar[0], comm);
        }

        void setEadgeBoundary(F64 & x_low, F64 & x_high, const bool flag_periodic){
            if(!flag_periodic){
                //x_low  = std::numeric_limits<F64>::lowest();
                //x_high = std::numeric_limits<F64>::max();
                x_low  = -LARGE_DOUBLE;
                x_high =  LARGE_DOUBLE;
            }
            else{
                x_low = pos_root_domain_.low_.x;
                x_high = pos_root_domain_.high_.x;
            }
        }
#endif
        
        void decomposeDomainMultiStep2(const bool flag_use_old_domain=false){
            //Comm::barrier();
            //sleep(1);
            F64 time_offset = GetWtime();
            const auto full_len_root_domain = pos_root_domain_.high_ - pos_root_domain_.low_;
            const auto my_rank = Comm::getRank();
            const auto n_proc  = Comm::getNumberOfProc();
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL // ifNdef
            pos_domain_[0] = pos_root_domain_;
#else
            const auto my_rank_x = rank_1d_[0];
            const auto n_proc_x  = n_proc_1d_[0];
            const auto my_rank_y = rank_1d_[1];
            const auto n_proc_y  = n_proc_1d_[1];
            const auto my_rank_z = rank_1d_[2];
            const auto n_proc_z  = n_proc_1d_[2];
            ////////////////
            /// X DIRECTION
            const auto n_smp_sub_x = Comm::getSum(pos_sample_loc_.size(), comm_sub_[0])/(n_proc_y*n_proc_z) ;
            ReallocatableArray<F64> x_smp_sub_in_slab(n_smp_sub_x, n_smp_sub_x, MemoryAllocMode::Pool);
            gatherSampleRandomly(pos_sample_loc_, x_smp_sub_in_slab, n_smp_sub_x, comm_sub_[0],
                                 [](F64vec & src, F64 & dst){dst = src.x;});
            // calculate x_partition
            ReallocatableArray<F64> x_partition(n_proc_x+1, n_proc_x+1, MemoryAllocMode::Pool);
            if(Comm::getRank(comm_sub_[0]) == 0){
                ReallocatableArray<F64> x_partition_tmp(n_proc_x+1, n_proc_x+1, MemoryAllocMode::Pool);
                if(flag_use_old_domain){
                    // reuse old domain
                    x_partition_tmp[0] = pos_domain_[0].low_.x;
                    for(S32 i=1; i<n_proc_x+1; i++){
                        const auto rank_src = (i-1)*n_proc_y*n_proc_z;
                        x_partition_tmp[i] = pos_domain_[rank_src].high_.x;
                    }
                }
                else{
                    // rougly calc postion of partition in x direction
                    determineRoughPartition(x_smp_sub_in_slab, x_partition_tmp, 2*n_proc_x, comm_1d_[0], false);
                    setEadgeBoundary(x_partition_tmp[0], x_partition_tmp[n_proc_x], periodic_axis_[0]);
                }
                exchangeSample(x_smp_sub_in_slab, x_partition_tmp, n_proc_x, my_rank_x, full_len_root_domain[0], periodic_axis_[0], comm_1d_[0], [](const F64 & src)->F64{return src;});
                std::sort(x_smp_sub_in_slab.frontP(), x_smp_sub_in_slab.endP());
#ifdef PS_DEBUG_DD
                for(S32 i=1; i<x_smp_sub_in_slab.size(); i++){
                    assert(x_partition_tmp[my_rank_x] <=  x_smp_sub_in_slab[i]);
                    assert(x_partition_tmp[my_rank_x+1] > x_smp_sub_in_slab[i]);
                    assert(x_smp_sub_in_slab[i-1] < x_smp_sub_in_slab[i]);
                }
#endif
                determinePartition(x_smp_sub_in_slab, x_partition_tmp, x_partition, comm_1d_[0]);
                setEadgeBoundary(x_partition[0], x_partition[n_proc_x], periodic_axis_[0]);
            }
            ////////////////////////////////////////
            // exchange samples in x directino for determination of y partition
            // broadcast x_partition to all processes
            Comm::broadcast(&x_partition[0], x_partition.size());
            // exchange samples so that each particles belong to own apropriate domain.
            exchangeSample(pos_sample_loc_, x_partition, n_proc_x, my_rank_x, full_len_root_domain[0], periodic_axis_[0], comm_1d_[0], [](const F64vec & src)->F64{return src.x;});
#ifdef PS_DEBUG_DD
            for(S32 i=0; i<pos_sample_loc_.size(); i++){
                assert(pos_sample_loc_[i].x < x_partition[my_rank_x+1]);
                assert(pos_sample_loc_[i].x >= x_partition[my_rank_x]);
            }
#endif

            ////////////////
            /// Y DIRECTION
            const auto n_smp_sub_y = Comm::getSum(pos_sample_loc_.size(), comm_1d_[2])/n_proc_z;
            ReallocatableArray<F64> y_smp_sub_in_rod(n_smp_sub_y, n_smp_sub_y, MemoryAllocMode::Pool);
            gatherSampleRandomly(pos_sample_loc_, y_smp_sub_in_rod, n_smp_sub_y, comm_1d_[2],
                                 [](F64vec & src, F64 & dst){dst = src.y;});
            // calc y_partition
            ReallocatableArray<F64> y_partition(n_proc_y+1, n_proc_y+1, MemoryAllocMode::Pool);
            if(my_rank_z == 0){
                ReallocatableArray<F64> y_partition_tmp(n_proc_y+1, n_proc_y+1, MemoryAllocMode::Pool);
                if(flag_use_old_domain){
                    // reuse old domain
                    y_partition_tmp[0] = pos_domain_[my_rank_x*n_proc_y*n_proc_z].low_.y;
                    for(S32 i=1; i<n_proc_y+1; i++){
                        const auto rank_src = my_rank_x*n_proc_y*n_proc_z + (i-1)*n_proc_z;
                        y_partition_tmp[i] = pos_domain_[rank_src].high_.y;
                    }
                }
                else{
                    // rougly calc postion of partition in y direction
                    determineRoughPartition(y_smp_sub_in_rod, y_partition_tmp, 2*n_proc_y, comm_1d_[1], true);
                    setEadgeBoundary(y_partition_tmp[0], y_partition_tmp[n_proc_y], periodic_axis_[1]);
                }
                exchangeSample(y_smp_sub_in_rod, y_partition_tmp, n_proc_y, my_rank_y, full_len_root_domain[1], periodic_axis_[1], comm_1d_[1], [](const F64 & src)->F64{return src;});

#ifdef PS_DEBUG_DD
                for(S32 i=1; i<y_smp_sub_in_rod.size(); i++){
                    assert(y_partition_tmp[my_rank_y] <=  y_smp_sub_in_rod[i]);
                    assert(y_partition_tmp[my_rank_y+1] > y_smp_sub_in_rod[i]);
                }
#endif
                std::sort(y_smp_sub_in_rod.frontP(), y_smp_sub_in_rod.endP());
#ifdef PS_DEBUG_DD
                for(S32 i=1; i<y_smp_sub_in_rod.size(); i++){
                    assert(y_smp_sub_in_rod[i-1] < y_smp_sub_in_rod[i]);
                }
#endif
                determinePartition(y_smp_sub_in_rod, y_partition_tmp, y_partition, comm_1d_[1], true);
                setEadgeBoundary(y_partition[0], y_partition[n_proc_y], periodic_axis_[1]);
            }
            
            ////////////////////////////////////////
            // exchange samples in y directino for determination of z partition
            // broadcast y_partition to processes with the same idx (y-z plane)
            Comm::broadcast(&y_partition[0], y_partition.size(), 0, comm_sub_[0]);
            exchangeSample(pos_sample_loc_, y_partition, n_proc_y, my_rank_y, full_len_root_domain[1], periodic_axis_[1], comm_1d_[1], [](const F64vec & src)->F64{return src.y;});
#ifdef PS_DEBUG_DD
            for(S32 i=0; i<pos_sample_loc_.size(); i++){
                assert(pos_sample_loc_[i].y < y_partition[my_rank_y+1]);
                assert(pos_sample_loc_[i].y >= y_partition[my_rank_y]);
            }
#endif

            ////////////////
            /// Z DIRECTION
            S32 n_smp_sub_z = pos_sample_loc_.size();
            ReallocatableArray<F64> z_smp_sub(n_smp_sub_z, n_smp_sub_z, MemoryAllocMode::Pool);
            for(S32 i=0; i<n_smp_sub_z; i++)
                z_smp_sub[i] = pos_sample_loc_[i].z;
            ReallocatableArray<F64> z_partition(n_proc_z+1, n_proc_z+1, MemoryAllocMode::Pool);
            // rougly calc postion of partition in z direction
            ReallocatableArray<F64> z_partition_tmp(n_proc_z+1, n_proc_z+1, MemoryAllocMode::Pool);
            if(flag_use_old_domain){
                // reuse old domain
                z_partition_tmp[0] = pos_domain_[my_rank_x*n_proc_y*n_proc_z+my_rank_y*n_proc_z].low_.z;
                for(S32 i=1; i<n_proc_z+1; i++){
                    const auto rank_src = my_rank_x*n_proc_y*n_proc_z + my_rank_y*n_proc_z + (i-1);
                    z_partition_tmp[i] = pos_domain_[rank_src].high_.z;
                }
            }
            else{
                // rougly calc postion of partition in y direction
                determineRoughPartition(z_smp_sub, z_partition_tmp, 2*n_proc_z, comm_1d_[2], true);
                setEadgeBoundary(z_partition_tmp[0], z_partition_tmp[n_proc_z], periodic_axis_[2]);
            }
            exchangeSample(z_smp_sub, z_partition_tmp, n_proc_z, my_rank_z, full_len_root_domain[2], periodic_axis_[2], comm_1d_[2], [](const F64 & src)->F64{return src;});
            
#ifdef PS_DEBUG_DD
            for(S32 i=1; i<z_smp_sub.size(); i++){
                assert(z_partition_tmp[my_rank_z] <=  z_smp_sub[i]);
                assert(z_partition_tmp[my_rank_z+1] > z_smp_sub[i]);
            }
#endif
            std::sort(z_smp_sub.frontP(), z_smp_sub.endP());
#ifdef PS_DEBUG_DD
            for(S32 i=1; i<z_smp_sub.size(); i++){
                assert(z_partition_tmp[my_rank_z] <=  z_smp_sub[i]);
                assert(z_partition_tmp[my_rank_z+1] > z_smp_sub[i]);
                assert(z_smp_sub[i-1] < z_smp_sub[i]);
            }
#endif
            determinePartition(z_smp_sub, z_partition_tmp, z_partition, comm_1d_[2], true);
            setEadgeBoundary(z_partition[0], z_partition[n_proc_z], periodic_axis_[2]);
            ////////////////
            // braod cast coordinates of partition
            ReallocatableArray<F64> x_partition_recv(n_proc_x+1, n_proc_x+1, MemoryAllocMode::Pool);
            if(my_rank==0){
                for(S32 i=0; i<n_proc_x+1; i++){                
                    x_partition_recv[i] = x_partition[i];
                }
            }
            Comm::broadcast(&x_partition_recv[0], n_proc_x+1, 0, MPI_COMM_WORLD);

            ReallocatableArray<F64> y_partition_recv((n_proc_y+1)*n_proc_x, (n_proc_y+1)*n_proc_x, MemoryAllocMode::Pool);
            Comm::gather(&y_partition[0], n_proc_y+1, &y_partition_recv[0], comm_1d_[0]);
            Comm::broadcast(&y_partition_recv[0], (n_proc_y+1)*n_proc_x, 0, MPI_COMM_WORLD);
            
            ReallocatableArray<F64> z_partition_recv((n_proc_z+1)*n_proc_y*n_proc_x, (n_proc_z+1)*n_proc_y*n_proc_x, MemoryAllocMode::Pool);
            Comm::gather(&z_partition[0], n_proc_z+1, &z_partition_recv[0], comm_sub_[2]);
            Comm::broadcast(&z_partition_recv[0], (n_proc_z+1)*n_proc_y*n_proc_x, 0, MPI_COMM_WORLD);

            ////////////////
            // set pos_domain_
            for(S32 i=0; i<n_proc_x; i++){
                for(S32 j=0; j<n_proc_y; j++){
                    for(S32 k=0; k<n_proc_z; k++){
                        const auto rank = i*n_proc_y*n_proc_z + j*n_proc_z + k;
                        pos_domain_temp_[rank].low_.x  = x_partition_recv[i];
                        pos_domain_temp_[rank].high_.x = x_partition_recv[i+1];
                        pos_domain_temp_[rank].low_.y  = y_partition_recv[i*(n_proc_y+1) + j];
                        pos_domain_temp_[rank].high_.y = y_partition_recv[i*(n_proc_y+1) + j+1];
                        pos_domain_temp_[rank].low_.z  = z_partition_recv[i*(n_proc_z+1)*n_proc_y + j*(n_proc_z+1) + k];
                        pos_domain_temp_[rank].high_.z = z_partition_recv[i*(n_proc_z+1)*n_proc_y + j*(n_proc_z+1) + k+1];
                        
                    }
                }
            }

            // set pos_domain by using exponential moving average
	    if(first_call_by_decomposeDomain) {
                first_call_by_decomposeDomain = false;
                for(S32 i = 0; i < n_proc; i++) {
                    pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                    pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                }
	    } else {
                for(S32 i = 0; i < n_proc; i++) {
                    pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                    pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                }
	    }

            
            //if(my_rank==1){
            //    for(S32 i=0; i<n_proc; i++){
            //        std::cerr<<"i= "<<i<<" pos_domain_temp_[i]= "<<pos_domain_temp_[i]<<std::endl;
            //    }
            //}
            
            /*
            ReallocatableArray<S32vec, 1> rank_3d(n_proc, n_proc, MemoryAllocMode::Pool);
            for(S32 i=0; i<n_proc; i++){
                rank_3d[i].x = i / (n_proc_y*n_proc_z);
                rank_3d[i].y = (i / n_proc_z) % n_proc_y ;
                rank_3d[i].z = i % n_proc_z;
            }

            for(S32 i=0; i<n_proc; i++){
                S32 rank_x = rank_3d[i].x;
                assert(pos_domain_temp_[i].low_.x  == x_partition_recv[rank_x]);
                assert(pos_domain_temp_[i].high_.x == x_partition_recv[rank_x+1]);
            }

            for(S32 i=0; i<n_proc; i++){
                S32 rank_x = rank_3d[i].x;
                S32 rank_y = rank_3d[i].y;
                assert(pos_domain_temp_[i].low_.y  == y_partition_recv[(rank_x*(n_proc_y+1))+rank_y]);
                assert(pos_domain_temp_[i].high_.y == y_partition_recv[(rank_x*(n_proc_y+1))+rank_y+1]);
            }
            */
            /*
            for(S32 i=0; i<n_proc; i++){
                S32 rank_x = rank_3d[i].x;
                S32 rank_y = rank_3d[i].y;
                S32 rank_z = rank_3d[i].z;
                assert(pos_domain_temp_[i].low_.z  == z_partition_recv[(rank_x*(n_proc_y+1))+rank_z]);
                assert(pos_domain_temp_[i].high_.z == z_partition_recv[(rank_x*(n_proc_y+1))+rank_z+1]);
            }
            */
#endif
            pos_sample_tot_.freeMem();
            pos_sample_loc_.freeMem();
	    time_profile_.decompose_domain = GetWtime() - time_offset;            
        }
    };
}

