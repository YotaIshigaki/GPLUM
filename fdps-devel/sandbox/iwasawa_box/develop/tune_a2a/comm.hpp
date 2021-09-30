
extern std::ofstream fout_debug;

void GetRank1D(int rank_1d[],
               const int n_proc_1d[],
               const int my_rank,
               const int dim=3){
    int rank_tmp = my_rank;
    for(int d=dim-1; d>=0; d--){
        rank_1d[d] = rank_tmp % n_proc_1d[d];
        rank_tmp /= n_proc_1d[d];
    }
}

int GetRankGlobal(const int n_proc_1d[],
                  const int rank_1d[],
                  const int dim){
    int rank = (rank_1d[0] + n_proc_1d[0]) % n_proc_1d[0];
    for(int d=1; d<dim; d++){
        rank *= n_proc_1d[d];
        rank += (rank_1d[d] + n_proc_1d[d]) % n_proc_1d[d];
    }
    return rank;
}

void SplitComm(MPI_Comm comm_1d[],
               const int n_proc_1d[],
               const int rank_1d[],
               const int my_rank,
               const int dim=3){
    int color = rank_1d[0];
    int radix = 1;
    for(int d=1; d<dim; d++){
        radix *= n_proc_1d[d-1];
        color += rank_1d[d]*radix;
    }
// processes with the same rank_1d[0] and rank_1d[1] are in the same new communicator.
    radix = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color-rank_1d[0], my_rank, comm_1d+0);
    for(int d=1; d<dim; d++){
        radix *= n_proc_1d[d-1];
        MPI_Comm_split(MPI_COMM_WORLD, color-rank_1d[d]*radix, my_rank, comm_1d+d);
    }
    int nx, ny, nz;
    int rx, ry, rz;
    MPI_Comm_size(comm_1d[0], &nx);
    MPI_Comm_rank(comm_1d[0], &rx);
    if(dim > 1){
        MPI_Comm_size(comm_1d[1], &ny);
        MPI_Comm_rank(comm_1d[1], &ry);
    }
    if(dim > 2){
        MPI_Comm_size(comm_1d[2], &nz);
        MPI_Comm_rank(comm_1d[2], &rz);
    }
}

template<class T>
class AlltoallSystem{
private:
    int dim_;
    int rank_;
    int n_proc_;
    int * n_proc_1d_;
    int * rank_1d_;
    MPI_Comm * comm_1d_;
    MPI_Request * req_send_;
    MPI_Request * req_recv_;
    MPI_Status * stat_send_;
    MPI_Status * stat_recv_;
    PS::ReallocatableArray<T> val_send_buf_;
    int * n_send_buf_;
    int * n_send_1d_;
    int * n_recv_1d_;
    int * n_send_disp_1d_;
    int * n_recv_disp_1d_;
    int size_crit_; // critical value of size
public:
    AlltoallSystem(const int n_proc_1d[], const int dim){
        rank_ = MPI::COMM_WORLD.Get_rank();
        n_proc_ = MPI::COMM_WORLD.Get_size();
        dim_ = dim;
        n_proc_1d_ = new int[dim_];
        for(int d=0; d<dim_; d++) n_proc_1d_[d] = n_proc_1d[d];
        req_send_ = new MPI_Request[n_proc_];
        req_recv_ = new MPI_Request[n_proc_];
        stat_send_ = new MPI_Status[n_proc_];
        stat_recv_ = new MPI_Status[n_proc_];
        rank_1d_ = new int[dim_];
        comm_1d_ = new MPI_Comm[dim_];
        n_send_buf_ = new int[n_proc_];
        n_send_1d_ = new int[n_proc_];
        n_recv_1d_ = new int[n_proc_];
        n_send_disp_1d_ = new int[n_proc_+1];
        n_recv_disp_1d_ = new int[n_proc_+1];
        fout_debug<<"rank_= "<<rank_<<std::endl;
        for(int d=0; d<dim_; d++){
            fout_debug<<"n_proc_1d_["<<d<<"]= "<<n_proc_1d_[d]<<"   ";
        }
        fout_debug<<std::endl;
        GetRank1D(rank_1d_, n_proc_1d_, rank_, dim_);
        for(int d=0; d<dim_; d++){
            fout_debug<<"rank_1d_["<<d<<"]= "<<rank_1d_[d]<<"   ";
        }
        fout_debug<<std::endl;
        SplitComm(comm_1d_, n_proc_1d_, rank_1d_, rank_, dim_);
        for(int i=0; i<dim_; i++){
            int n = 0;
            MPI_Comm_size(comm_1d_[i], &n);
            fout_debug<<"size_1d[i]= "<<n<<"  ";
        }
        fout_debug<<std::endl;
        for(int i=0; i<dim_; i++){
            int n = 0;
            MPI_Comm_rank(comm_1d_[i], &n);
            fout_debug<<"rank_1d[i]= "<<n<<"  ";
        }
        fout_debug<<std::endl;
        size_crit_ = 0;
    }

    const int * getRank1D(){
        return rank_1d_;
    }

    const int * getNProc1D(){
        return n_proc_1d_;
    }

    void allToAllNormal(T * val_send, T * val_recv, const int n){
        MPI_Alltoall(val_send, n, PS::GetDataType<T>(),
                     val_recv, n, PS::GetDataType<T>(),
                     MPI_COMM_WORLD);
    }

    void allToAllVNormal(T * val_send, int * n_send, int * n_send_disp,
                         T * val_recv, int * n_recv, int * n_recv_disp){
        MPI_Alltoallv(val_send, n_send, n_send_disp, PS::GetDataType<T>(),
                      val_recv, n_recv, n_recv_disp, PS::GetDataType<T>(),
                      MPI_COMM_WORLD);
    }

    void allToAllMD(T * val_send, T * val_recv, const int size){
        // multi dimension
        int radix = n_proc_ / n_proc_1d_[0];
        MPI_Alltoall(val_send, size*radix, PS::GetDataType<T>(),
                     val_recv, size*radix, PS::GetDataType<T>(), 
                     comm_1d_[0]);
/*
        if(PS::Comm::getRank() == 0){
            for(int i=0; i<n_proc_*size; i++){
                std::cout<<"i= "<<i<<" val_send[i]= "<<val_send[i]<<std::endl;
            }
            for(int i=0; i<n_proc_*size; i++){
                std::cout<<"i= "<<i<<" val_recv[i]= "<<val_recv[i]<<std::endl;
            }
        }
*/
        if(dim_ > 1){
            val_send_buf_.resizeNoInitialize(size*n_proc_);
            int n_cnt = 0;
            for(int d=1; d<dim_; d++){
#if 0
                n_cnt = 0;
                radix = n_proc_ / n_proc_1d_[d-1];
                for(int i=0; i<radix; i++){
                    int ioffset = size*i;
                    for(int j=0; j<n_proc_1d_[d-1]; j++){
                        int joffset = radix*size*j;
                        for(int k=0; k<size; k++){
                            val_send_buf_[n_cnt++] = val_recv[ioffset+joffset+k];
                        }
                    }
                }
#else
                n_cnt = 0;
                radix = n_proc_ / n_proc_1d_[d-1];
#pragma omp parallel for
                for(int i=0; i<radix; i++){
                    int ioffset = i*size;
                    int ioffset_dst = i*n_proc_1d_[d-1]*size;
                    for(int j=0; j<n_proc_1d_[d-1]; j++){
                        int joffset = radix*size*j;
                        int joffset_dst = j*size;
                        for(int k=0; k<size; k++){
                            //std::cerr<<"n_cnt= "<<n_cnt<<" ioffset_dst+joffset_dst+k= "<<ioffset_dst+joffset_dst+k<<std::endl;
                            val_send_buf_[ioffset_dst+joffset_dst+k] = val_recv[ioffset+joffset+k];
                            //assert(n_cnt == ioffset_dst+joffset_dst+k);
                            //n_cnt++;
                        }
                    }
                }
#endif
                radix = n_proc_ / n_proc_1d_[d];
                PS::Comm::barrier();
                MPI_Alltoall(val_send_buf_.getPointer(), size*radix, PS::GetDataType<T>(),
                             val_recv,                   size*radix, PS::GetDataType<T>(), 
                             comm_1d_[d]);
/*
                if(PS::Comm::getRank() == 0){
                    for(int i=0; i<n_proc_*size; i++){
                        std::cout<<"i= "<<i<<" val_send_buf_[i]= "<<val_send_buf_[i]<<std::endl;
                    }
                    for(int i=0; i<n_proc_*size; i++){
                        std::cout<<"i= "<<i<<" val_recv[i]= "<<val_recv[i]<<std::endl;
                    }
                }
*/
            }
#if 0
            n_cnt = 0;
            radix = n_proc_ / n_proc_1d_[dim_-1];
            for(int i=0; i<radix; i++){
                int ioffset = size*i;
                for(int j=0; j<n_proc_1d_[dim_-1]; j++){
                    int joffset = radix*size*j;
                    for(int k=0; k<size; k++){
                        val_send_buf_[n_cnt++] = val_recv[ioffset+joffset+k];
                    }
                }
            }
#else
            n_cnt = 0;
            radix = n_proc_ / n_proc_1d_[dim_-1];
#pragma omp parallel for
            for(int i=0; i<radix; i++){
                int ioffset = i*size;
                int ioffset_dst = i*n_proc_1d_[dim_-1]*size;
                for(int j=0; j<n_proc_1d_[dim_-1]; j++){
                    int joffset = radix*size*j;
                    int joffset_dst = j*size;
                    for(int k=0; k<size; k++){
                        //std::cerr<<"n_cnt= "<<n_cnt<<"ioffset_dst+joffset_dst+k= "<<ioffset_dst+joffset_dst+k<<std::endl;
                        val_send_buf_[ioffset_dst+joffset_dst+k] = val_recv[ioffset+joffset+k];
                        //assert(n_cnt == ioffset_dst+joffset_dst+k);
                        //n_cnt++;
                    }
                }
            }
#endif
#pragma omp parallel for
            for(int i=0; i<n_proc_*size; i++){
                val_recv[i] = val_send_buf_[i];
            }
        }
    }

    void allToAllI(T * val_send, T * val_recv, const int n){
        int n_cnt_send = 0;
        int n_cnt_recv = 0;
        for(int i=1; i<n_proc_; i++){
            int id = (rank_+n_proc_+i) % n_proc_;
            req_send_[n_cnt_send++] = MPI::COMM_WORLD.Isend(val_send+n*id, n, PS::GetDataType<T>(), id, 0);
            req_recv_[n_cnt_recv++] = MPI::COMM_WORLD.Irecv(val_recv+n*id, n, PS::GetDataType<T>(), id, 0);
        }
        MPI_Waitall(n_cnt_send, req_send_, stat_send_);
        MPI_Waitall(n_cnt_recv, req_recv_, stat_recv_);
    }


    void allToAllSwitch(T * val_send, T * val_recv, const int size){
        if(size < size_crit_){
            this->allToAllMD(val_send, val_recv, size);
        }
        else{
            this->allToAllNormal(val_send, val_recv, size);
        }
    }

    void allToAllVMD(T * val_send, int * n_send, int * n_send_disp,
                     T * val_recv, int * n_recv, int * n_recv_disp){
        // multi dimension
/*
        for(int i=0; i<n_proc_; i++){
            fout_debug<<"n_send[i]= "<<n_send[i]
                      <<" n_send_disp[i]= "<<n_send_disp[i]
                      <<std::endl;
        }
*/
        int radix = n_proc_ / n_proc_1d_[0];
        PS::Comm::barrier();
        MPI_Alltoall(n_send, radix, PS::GetDataType<int>(),
                     n_recv, radix, PS::GetDataType<int>(),
                     comm_1d_[0]);
        n_recv_disp[0] = 0;
        for(int i=0; i<n_proc_; i++) n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        for(int i=0; i<n_proc_1d_[0]; i++){
            n_send_1d_[i] = n_recv_1d_[i] = 0;
            int ioffset = radix * i;
            for(int j=0; j<radix; j++){
                n_send_1d_[i] += n_send[ioffset+j];
                n_recv_1d_[i] += n_recv[ioffset+j];
            }
        }
        n_send_disp_1d_[0] = n_recv_disp_1d_[0] = 0;
        for(int i=0; i<n_proc_1d_[0]; i++){
            n_send_disp_1d_[i+1] = n_send_disp_1d_[i] + n_send_1d_[i];
            n_recv_disp_1d_[i+1] = n_recv_disp_1d_[i] + n_recv_1d_[i];
        }
        PS::Comm::barrier();
        MPI_Alltoallv(val_send, n_send_1d_, n_send_disp_1d_, PS::GetDataType<T>(),
                      val_recv, n_recv_1d_, n_recv_disp_1d_, PS::GetDataType<T>(),
                      comm_1d_[0]);
/*
        fout_debug<<"alltoallvmd"<<std::endl;
        for(int i=0; i<n_send_disp_1d_[n_proc_1d_[0]]; i++){
            fout_debug<<"i= "<<i<<" val_send[i]= "<<val_send[i]<<std::endl;
        }
        for(int i=0; i<n_recv_disp_1d_[n_proc_1d_[0]]; i++){
            fout_debug<<"i= "<<i<<" val_recv[i]= "<<val_recv[i]<<std::endl;
        }
*/
        if(dim_ > 1){
            for(int d=1; d<dim_; d++){
                val_send_buf_.resizeNoInitialize(n_recv_disp[n_proc_]);
                radix = n_proc_ / n_proc_1d_[d-1];
                int n_cnt_0 = 0;
                int n_cnt_1 = 0;
                for(int i=0; i<radix; i++){
                    for(int j=0; j<n_proc_1d_[d-1]; j++){
                        n_send_buf_[n_cnt_0++] = n_recv[j*radix+i];
                        int offset = n_recv_disp[j*radix+i];
                        for(int k=0; k<n_recv[j*radix+i]; k++, n_cnt_1++){
                            val_send_buf_[n_cnt_1] = val_recv[offset+k];
                        }
                    }
                }
                radix = n_proc_ / n_proc_1d_[d];
                PS::Comm::barrier();
                MPI_Alltoall(n_send_buf_, radix, PS::GetDataType<int>(),
                             n_recv,      radix, PS::GetDataType<int>(),
                             comm_1d_[d]);
                n_recv_disp[0] = 0;
                for(int i=0; i<n_proc_; i++) n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
                for(int i=0; i<n_proc_1d_[d]; i++){
                    n_send_1d_[i] = n_recv_1d_[i] = 0;
                    int ioffset = radix * i;
                    for(int j=0; j<radix; j++){
                        n_send_1d_[i] += n_send_buf_[ioffset+j];
                        n_recv_1d_[i] += n_recv[ioffset+j];
                    }
                }
                n_send_disp_1d_[0] = n_recv_disp_1d_[0] = 0;
                for(int i=0; i<n_proc_1d_[d]; i++){
                    n_send_disp_1d_[i+1] = n_send_disp_1d_[i] + n_send_1d_[i];
                    n_recv_disp_1d_[i+1] = n_recv_disp_1d_[i] + n_recv_1d_[i];
                }
                PS::Comm::barrier();
                MPI_Alltoallv(val_send_buf_.getPointer(), n_send_1d_, n_send_disp_1d_, PS::GetDataType<T>(),
                              val_recv,                   n_recv_1d_, n_recv_disp_1d_, PS::GetDataType<T>(),
                              comm_1d_[d]);
/*
                if(PS::Comm::getRank() == 0){
                    fout_debug<<"alltoallvmd, d="<<d<<std::endl;
                    for(int i=0; i<n_send_disp_1d_[n_proc_1d_[d]]; i++){
                        fout_debug<<"i= "<<i<<" val_send_buf_[i]= "<<val_send_buf_[i]<<std::endl;
                    }
                    for(int i=0; i<n_recv_disp_1d_[n_proc_1d_[d]]; i++){
                        fout_debug<<"i= "<<i<<" val_recv[i]= "<<val_recv[i]<<std::endl;
                    }
                }
*/
            }
            val_send_buf_.resizeNoInitialize(n_recv_disp[n_proc_]);
            radix = n_proc_ / n_proc_1d_[dim_-1];
            int n_cnt_0 = 0;
            int n_cnt_1 = 0;
            for(int i=0; i<radix; i++){
                for(int j=0; j<n_proc_1d_[dim_-1]; j++){
                    n_send_buf_[n_cnt_0++] = n_recv[j*radix+i];
                    int offset = n_recv_disp[j*radix+i];
                    for(int k=0; k<n_recv[j*radix+i]; k++){
                        val_send_buf_[n_cnt_1++] = val_recv[offset+k];
                    }
                }
            }
            n_recv_disp[0] = 0;
            for(int i=0; i<n_proc_; i++) n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            for(int i=0; i<n_recv_disp[n_proc_]; i++){
                val_recv[i] = val_send_buf_[i];
            }
        }
    }

    void allToAllVMS(T * val_send, int * n_send, int * n_send_disp,
                     T * val_recv, int * n_recv, int * n_recv_disp,
                     const int size_msg_a2a = 1){
        // combination of alltoall and alltoallv
        static PS::ReallocatableArray<T>   val_send_1st;
        static PS::ReallocatableArray<T>   val_recv_1st;
        val_send_1st.resizeNoInitialize(n_proc_*size_msg_a2a);
        val_recv_1st.resizeNoInitialize(n_proc_*size_msg_a2a);
#pragma omp parallel for
        for(int i=0; i<n_proc_; i++){
            const int n = std::min(size_msg_a2a, n_send[i]);
            const int offset = i*size_msg_a2a;
            const int offset_src = n_send_disp[i];
            for(int j=0; j<n; j++){
                val_send_1st[offset+j] = val_send[offset_src+j];
            }
            for(int j=n; j<size_msg_a2a; j++){
                val_send_1st[offset+j] = 0;
            }
        }

        this->allToAllSwitch(val_send_1st.getPointer(), val_recv_1st.getPointer(), size_msg_a2a);
#pragma omp parallel for
        for(int i=0; i<n_proc_; i++){
            int n = std::min(size_msg_a2a, n_send[i]);
            const int offset_dst = n_recv_disp[i];
            const int offset_src = i*size_msg_a2a;
            for(int j=0; j<n; j++){
                val_recv[offset_dst+j] = val_recv_1st[offset_src+j];
            }
        }
        int n_cnt_send = 0;
        int n_cnt_recv = 0;
        for(int i=1; i<n_proc_; i++){
            int id = (rank_+n_proc_+i) % n_proc_;
            if(n_send[id] > size_msg_a2a){
                req_send_[n_cnt_send++] = MPI::COMM_WORLD.Isend(val_send+n_send_disp[id]+size_msg_a2a, n_send[id]-size_msg_a2a, PS::GetDataType<T>(), id, 0);
            }
            if(n_recv[id] > size_msg_a2a){
                req_recv_[n_cnt_recv++] = MPI::COMM_WORLD.Irecv(val_recv+n_recv_disp[id]+size_msg_a2a, n_recv[id]-size_msg_a2a, PS::GetDataType<T>(), id, 0);
            }
        }
        MPI_Waitall(n_cnt_send, req_send_, stat_send_);
        MPI_Waitall(n_cnt_recv, req_recv_, stat_recv_);
    }

    void allToAllVI(T * val_send, int * n_send, int * n_send_disp,
                    T * val_recv, int * n_recv, int * n_recv_disp){
        int n_cnt_send = 0;
        int n_cnt_recv = 0;
        for(int i=1; i<n_proc_; i++){
            int id = (rank_+n_proc_+i) % n_proc_;
            req_send_[n_cnt_send++] = MPI::COMM_WORLD.Isend(val_send+n_send_disp[id], n_send[id], PS::GetDataType<T>(), id, 0);
            req_recv_[n_cnt_recv++] = MPI::COMM_WORLD.Irecv(val_recv+n_recv_disp[id], n_recv[id], PS::GetDataType<T>(), id, 0);
        }
        MPI_Waitall(n_cnt_send, req_send_, stat_send_);
        MPI_Waitall(n_cnt_recv, req_recv_, stat_recv_);
    }

    void initLoop(int size=16, int n_loop=20){
        T * val_send = new T[size*n_proc_];
        T * val_recv = new T[size*n_proc_];
        int * n_send = new int[n_proc_];
        int * n_recv = new int[n_proc_];
        int * n_send_disp = new int[n_proc_+1];
        int * n_recv_disp = new int[n_proc_+1];
        n_send_disp[0] = n_recv_disp[0] = 0;
        for(int i=0; i<n_proc_; i++){
            n_send[i] = n_recv[i] = size;
            n_send_disp[i+1] = n_send_disp[i] + n_send[i];
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        }
        for(int i=0; i<n_loop; i++){
            allToAllNormal(val_send, val_recv, size);
            allToAllMD(val_send, val_recv, size);
            allToAllVNormal(val_send, n_send, n_send_disp, val_recv, n_recv, n_recv_disp);
            allToAllVMD(val_send, n_send, n_send_disp, val_recv, n_recv, n_recv_disp);
            allToAllVMS(val_send, n_send, n_send_disp, val_recv, n_recv, n_recv_disp, size);
        }
        delete [] val_send;
        delete [] val_recv;
        delete [] n_send;
        delete [] n_recv;
        delete [] n_send_disp;
        delete [] n_recv_disp;
    }

    void setSizeCrit(const int size_max){
        const int size = 1;
        size_crit_ = 0;
        T * val_send = new T[size_max*n_proc_];
        T * val_recv = new T[size_max*n_proc_];
        for(int i=0; i<size_max*n_proc_; i++){
            val_send[i] = val_recv[i] = 0;
        }
        PS::F64 wtime_normal = 0.0;
        PS::F64 wtime_md = 0.0;
        for(int s=size; s<size_max; s*=2){
            for(int l=0; l<5; l++){
                PS::F64 wtime_offset = PS::GetWtime();
                allToAllNormal(val_send, val_recv, s);
                PS::F64 wtime_tmp = PS::GetWtime() - wtime_offset;
                if(l==0) wtime_normal = wtime_tmp;
                else wtime_normal = std::min(wtime_normal, wtime_tmp);
            }
            wtime_normal = PS::Comm::getMinValue(wtime_normal);
            for(int l=0; l<5; l++){
                PS::F64 wtime_offset = PS::GetWtime();
                allToAllMD(val_send, val_recv, s);
                PS::F64 wtime_tmp = PS::GetWtime() - wtime_offset;
                if(l==0) wtime_md = wtime_tmp;
                else wtime_md = std::min(wtime_md, wtime_tmp);
            }
            wtime_md = PS::Comm::getMinValue(wtime_md);
            if(wtime_md > wtime_normal){
                size_crit_ = s;
                break;
            }
        }
        fout_debug<<"size_crit_= "<<size_crit_<<std::endl;
        delete [] val_send;
        delete [] val_recv;
    }
};


