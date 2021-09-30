#pragma once
#include<mpi.h>
#include"../../../src/ps_defs.hpp"

namespace ParticleSimulator{

    template<int DIM>
    void DivideProc(int np[],
                    int rank[],
                    const int nproc,
                    const int myrank){
        std::vector<int> npv;
        npv.resize(DIM);
        int np_tmp = nproc;
        for(int d=DIM, cid=0; cid<DIM-1; d--, cid++){
            int tmp = (int)pow(np_tmp+0.000001, (1.0/d)*1.000001 );
            while(np_tmp%tmp){
                tmp--;
            }
            npv[cid] = tmp;
            np_tmp /= npv[cid];
        }
        npv[DIM-1] = np_tmp;
        int rank_tmp = myrank;
        std::sort(npv.begin(), npv.end(), std::greater<int>());
        for(int i=DIM-1; i>=0; i--){
            np[i] = npv[i];
            rank[i] = rank_tmp % np[i];
            rank_tmp /= np[i];
        }
    }
    
    template<class T, int DIM_COMM>
    class CommForAllToAll{
    private:
        int rank_glb_;
        int n_proc_glb_;
        MPI_Comm comm_glb_;
        int rank_1d_[DIM_COMM];
        int n_proc_1d_[DIM_COMM];
        MPI_Comm comm_1d_[DIM_COMM];
        ReallocatableArray<T> val_send_glb_;
        int * n_recv_disp_glb_;
        int * n_send_disp_glb_;
        int * n_send_glb_;
        int * n_recv_1d_;
        int * n_send_1d_;
        int * n_recv_disp_1d_;
        int * n_send_disp_1d_;
    public:
        void dumpRank(){
            for(int i=0; i<DIM_COMM; i++){
                int n_tmp = 0;
                MPI_Comm_size(comm_1d_[i], &n_tmp);
                std::cout<<"n_proc_1d_["<<i<<"]="<<n_proc_1d_[i]<<" n_tmp="<<n_tmp<<std::endl;
	    }
            for(int i=0; i<DIM_COMM; i++){
                int r_tmp = 0;
                MPI_Comm_rank(comm_1d_[i], &r_tmp);
                std::cout<<"rank_1d_["<<i<<"]="<<rank_1d_[i]<<" r_tmp="<<r_tmp<<std::endl;
            }
        }
	
        CommForAllToAll(MPI_Comm comm = MPI_COMM_WORLD){
            comm_glb_ = comm;
            MPI_Comm_rank(comm_glb_, &rank_glb_);
            MPI_Comm_size(comm_glb_, &n_proc_glb_);
            n_recv_disp_glb_ = new int[n_proc_glb_ + 1];
            n_send_disp_glb_ = new int[n_proc_glb_ + 1];
            n_send_glb_ = new int[n_proc_glb_];
            n_recv_1d_ = new int[n_proc_glb_];
            n_send_1d_ = new int[n_proc_glb_];
            n_recv_disp_1d_ = new int[n_proc_glb_ + 1];
            n_send_disp_1d_ = new int[n_proc_glb_ + 1];
            DivideProc<DIM_COMM>(n_proc_1d_, rank_1d_, n_proc_glb_, rank_glb_);
            //int dim_max = n_proc_1d_[0]; // as radix
            int dim_max = -1;
            for(int i=0; i<DIM_COMM; i++){
                if(dim_max < n_proc_1d_[i]){
                    dim_max = n_proc_1d_[i];
                }
            }
            int split_color = 0;
            int factor = 1;
	    
            for(int d=0; d<DIM_COMM; d++){
                split_color += rank_1d_[d] * factor;
                factor *= dim_max;
            }
            for(int d=DIM_COMM-1; d>=0; d--){
                factor = rank_1d_[d];
                for(int d0=0; d0<d; d0++){
                    factor *= dim_max;
                }
                MPI_Comm_split(comm_glb_, split_color-factor, rank_glb_, comm_1d_+d);
            }
        } // Constructor
	
	void execute(const ReallocatableArray<T> & val_send,
		     ReallocatableArray<T> & val_recv,
		     const int cnt){
	    const int n_recv_tot = cnt * n_proc_glb_;
	    val_recv.resizeNoInitialize( n_recv_tot );
	    val_send_glb_.resizeNoInitialize( n_recv_tot );
	    for(int i=0; i<n_recv_tot; i++){
		val_recv[i] = val_send[i];	
	    }
	    for(int d=DIM_COMM-1; d>=0; d--){
		const int radix = n_proc_glb_ / n_proc_1d_[d];
		for(int ib=0; ib<n_proc_glb_; ib++){
		    const int id_send = cnt * ( (ib % radix) * n_proc_1d_[d] + ib / radix );
		    const int offset = ib * cnt;
		    for(int i=0; i<cnt; i++){
			val_send_glb_[offset + i] = val_recv[id_send + i];
		    }
		}
		MPI_Alltoall(val_send_glb_.getPointer(), cnt*radix, GetDataType<T>(),
			     val_recv.getPointer(), cnt*radix, GetDataType<T>(), comm_1d_[d]);
	    }
	}
	
	void executeV(const ReallocatableArray<T> & val_send,
		      ReallocatableArray<T> & val_recv,
		      const int n_send[],
		      int n_recv[]){
	    int cnt = 0;
	    val_recv.reserveAtLeast( val_send.capacity() );
	    val_recv.clearSize();
	    for(int ib=0; ib<n_proc_glb_; ib++){
		n_recv[ib] = n_send[ib];
		for(int ip=0; ip<n_recv[ib]; ip++, cnt++){
		    val_recv.pushBackNoCheck(val_send[cnt]);
		}
	    }

	    for(int d=DIM_COMM-1; d>=0; d--){
		int radix = n_proc_glb_ / n_proc_1d_[d];
		n_recv_disp_glb_[0] = 0;
		for(int i=0; i<n_proc_glb_; i++){
		    n_recv_disp_glb_[i+1] = n_recv_disp_glb_[i] + n_recv[i];
		}
		val_send_glb_.clearSize();
		for(int ib0=0; ib0<n_proc_glb_; ib0++){
		    int id_send = (ib0 % radix) * n_proc_1d_[d] + ib0 / radix;
		    n_send_glb_[ib0] = n_recv[id_send];
		    int offset = n_recv_disp_glb_[id_send];
		    val_send_glb_.reserveEmptyAreaAtLeast(n_send_glb_[ib0]);
		    for(int ib1=0; ib1<n_send_glb_[ib0]; ib1++){
			val_send_glb_.pushBackNoCheck(val_recv[ib1 + offset]);
		    }
		}
		MPI_Alltoall(n_send_glb_, radix, MPI_INT,
			     n_recv, radix, MPI_INT, comm_1d_[d]);
		n_send_disp_1d_[0] = n_recv_disp_1d_[0] = 0;
		for(int ib0=0; ib0<n_proc_1d_[d]; ib0++){
		    n_send_1d_[ib0] = n_recv_1d_[ib0] = 0;
		    int offset = ib0 * radix;
		    for(int ib1=0; ib1<radix; ib1++){
			n_send_1d_[ib0] += n_send_glb_[offset + ib1];
			n_recv_1d_[ib0] += n_recv[offset + ib1];
		    }
		    n_send_disp_1d_[ib0+1] = n_send_disp_1d_[ib0] + n_send_1d_[ib0];
		    n_recv_disp_1d_[ib0+1] = n_recv_disp_1d_[ib0] + n_recv_1d_[ib0];
		}
		val_recv.resizeNoInitialize( n_recv_disp_1d_[n_proc_1d_[d]] );
		MPI_Alltoallv(val_send_glb_.getPointer(), n_send_1d_, n_send_disp_1d_, GetDataType<T>(),
			      val_recv.getPointer(), n_recv_1d_, n_recv_disp_1d_, GetDataType<T>(), comm_1d_[d]);
	    }
	    //n_recv_tot = n_recv_disp_1d_[ n_proc_1d_[0] ];

	}

	
    };
}
