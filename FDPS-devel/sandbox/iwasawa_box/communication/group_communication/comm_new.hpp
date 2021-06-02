#ifndef __COMM_NEW_HPP__
#define __COMM_NEW_HPP__
#include<vector>
#include<algorithm>
#include<cmath>
#include"mpi.h"
#include"omp.h"

template<class T> inline MPI::Datatype GetDataType(){return MPI::INT;};
template<> inline MPI::Datatype GetDataType<int>(){return MPI::INT;}
template<> inline MPI::Datatype GetDataType<unsigned int>(){return MPI::UNSIGNED;}
template<> inline MPI::Datatype GetDataType<long long int>(){return MPI::LONG;}
template<> inline MPI::Datatype GetDataType<unsigned long long int>(){return MPI::UNSIGNED_LONG;}
template<> inline MPI::Datatype GetDataType<float>(){return MPI::FLOAT;}
template<> inline MPI::Datatype GetDataType<double>(){return MPI::DOUBLE;}



template<int DIM>
void MakeDivision(int np[],
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

template<class T, int DIM>
class Comm{
private:
    int my_rank_;
    int n_proc_;
    MPI::Intracomm comm_glb_;
    int rank_[DIM];
    int np_[DIM];
    MPI::Intracomm comm_[DIM];
    std::vector<MPI::Request> req_send_;
    std::vector<MPI::Request> req_recv_;
    std::vector<T> val_send_;
    std::vector<T> val_recv_;
    std::vector<T> val_buff_;
    int factor_ = nproc/np[DIM-1];
public:
    Comm(){}
    Comm(MPI::Intracomm _comm_glb){
        comm_glb_ = _comm_glb;
        my_rank_ = comm_glb_.Get_rank();
        n_proc_ = comm_glb_.Get_size();
        req_send_.reserve[n_proc_];
        req_recv_.reserve[n_proc_];
        val_send_.reserve[n_proc_*10+10000];
        val_recv_.reserve[n_proc_*10+10000];
        val_buff_.reserve[n_proc_*10+10000];


        MakeDivision<DIM>(np, rank, nproc, myrank);
        int dim_max = np[0];
        int dim_tmp = dim_max/DIM;
        if(rank[0]%DIM > dim_max/DIM){

        }

        int split_key = 0;
        int factor = 1;
        req_send = new MPI::Request[nproc];
        req_recv = new MPI::Request[nproc];
        for(int d=0; d<DIM; d++){
            split_key += rank[d]*factor;
            factor *= dim_max;
        }
        for(int d=DIM-1; d>=0; d--){
            factor = rank[d];
            for(int d0=0; d0<d; d0++){
                factor *= dim_max;
            }
            comm[d] = comm_glb.Split(split_key-factor, myrank);
        }
    }

    void Init(MPI::Intracomm _comm_glb){
        comm_glb = _comm_glb;
        myrank = comm_glb.Get_rank();
        nproc = comm_glb.Get_size();
        MakeDivision<int, DIM>(np, rank, nproc, myrank);
        int dim_max = np[0];
        int split_key = 0;
        int factor = 1;
        for(int d=0; d<DIM; d++){
            split_key += rank[d]*factor;
            factor *= dim_max;
        }
        for(int d=DIM-1; d>=0; d--){
            factor = rank[d];
            for(int d0=0; d0<d; d0++){
                factor *= dim_max;
            }
            comm[d] = comm_glb.Split(split_key-factor, myrank);
        }
    }

    void GetRank(int rank_out[DIM]){
        for(int i=0; i<DIM; i++){
            rank_out[i] = rank[i];
        }
    }

    void GetSize(int np_out[DIM]){
        for(int i=0; i<DIM; i++){
            np_out[i] = np[i];
        }
    }

    void GetDim(int & dim_out){
        dim_out = DIM;
    }

    void AllToAll_s(const int & cnt){
        for(int i=0; i<nproc-1; i++){
            int id_send = (i + myrank + 1) % nproc;
            int id_recv = (nproc + myrank - i - 1) % nproc;
            comm.Sendrecv(val_send+(cnt*id_send), cnt, type, id_send, myrank,
                          val_recv+(cnt*id_recv), cnt, type, id_recv, id_recv);
        }
        for(int i=0; i<cnt; i++){
            val_recv[cnt*myrank+i] = val_send[cnt*myrank+i];
        }
    }

    // in this function, normal alltoall is used.
    void AllToAll( const int & cnt ){
#ifdef __OPENMP__
#pragma omp parallel
#endif
        {
#ifdef __OPENMP__
#pragma omp for
#endif
            for(int ib0=0; ib0<n_proc_; ib0++){
                int id_send = cnt*( (ib0%factor_)*np[DIM-1] + ib0/factor_ );
                for(int ib1=0; ib1<cnt; ib1++){
                    val_buff_[ib0*cnt+ib1].push_back(val_send_[id_send+ib1]);
                }
            }
#ifdef __OPENMP__
#pragma omp master
#endif
            comm_[DIM-1].Alltoall(val_buf_, (int)(cnt*factor_), type, val_recv_, (int)(cnt*factor), type);
            for(int d=DIM-2; d>=0; d--){
                factor_ = n_proc/np[d];
#ifdef __OPENMP__
#pragma omp for
#endif
                for(int ib0=0; ib0<nproc; ib0++){
                    int id_send = cnt*( (ib0%factor)*np[d] + ib0/factor );
                    for(int ib1=0; ib1<cnt; ib1++){
                        val_buf[ib0*cnt+ib1] = val_recv[id_send+ib1];
                    }
                }
#ifdef __OPENMP__
#pragma omp master
#endif
                comm[d].Alltoall(val_buf, (int)(cnt*factor), type, val_recv, (int)(cnt*factor), type);
            }
        }
    }


}


#endif
