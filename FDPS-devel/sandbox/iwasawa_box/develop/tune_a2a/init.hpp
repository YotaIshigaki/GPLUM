extern std::ofstream fout_debug;

template<class T>
void SetValueUni(T val_send[],
                 const int n_proc,
                 const int size){
    int radix = 100000; // it works # of nodes less than 10000
    for(int i=0; i<n_proc; i++){
        for(int j=0; j<size; j++){
            val_send[i*size+j] = i*radix + PS::Comm::getRank();
        }
    }
}


template<class T>
void SetValueUnUni(T val_send[], 
                   int n_send[], 
                   int n_send_disp[],
                   const int n_proc_1d[],
                   const int n_proc, 
                   const int size_close,
                   const int size_distant,
                   const int dim, 
                   const int n_neighbor_1d){
    const int dim_max = 3;
    int rank_1d[dim_max];
    GetRank1D(rank_1d, n_proc_1d, PS::Comm::getRank(), dim);
    int rank_ngb_1d[dim_max];
    for(int i=0; i<n_proc; i++) n_send[i] = size_distant;
    if(dim == 3){
        for(int i=-n_neighbor_1d; i<=n_neighbor_1d; i++){
            rank_ngb_1d[0] = (rank_1d[0]+n_proc_1d[0]+i) % n_proc_1d[0];
            for(int j=-n_neighbor_1d; j<=n_neighbor_1d; j++){
                rank_ngb_1d[1] = (rank_1d[1]+n_proc_1d[1]+j) % n_proc_1d[1];
                for(int k=-n_neighbor_1d; k<=n_neighbor_1d; k++){
                    rank_ngb_1d[2] = (rank_1d[2]+n_proc_1d[2]+k) % n_proc_1d[2];
                    n_send[GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)] = size_close;
                    //fout_debug<<"GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)= "<<GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)<<std::endl;
                }
            }
        }
    }
    else if(dim == 2){
        for(int i=-n_neighbor_1d; i<=n_neighbor_1d; i++){
            rank_ngb_1d[0] = (rank_1d[0]+n_proc_1d[0]+i) % n_proc_1d[0];
            for(int j=-n_neighbor_1d; j<=n_neighbor_1d; j++){
                rank_ngb_1d[1] = (rank_1d[1]+n_proc_1d[1]+j) % n_proc_1d[1];
                rank_ngb_1d[2] = 0;
                n_send[GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)] = size_close;
                //fout_debug<<"GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)= "<<GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)<<std::endl;
            }
        }
    }
    else if(dim == 1){
        for(int i=-n_neighbor_1d; i<=n_neighbor_1d; i++){
            rank_ngb_1d[0] = (rank_1d[0]+n_proc_1d[0]+i) % n_proc_1d[0];
            rank_ngb_1d[1] = 0;
            rank_ngb_1d[2] = 0;
            n_send[GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)] = size_close;
            //fout_debug<<"GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)= "<<GetRankGlobal(n_proc_1d, rank_ngb_1d, dim)<<std::endl;
        }
    }
    int radix = 100000; // it works # of nodes less than 10000
    n_send_disp[0] = 0;
    for(int i=0; i<n_proc; i++){
        n_send_disp[i+1] = n_send_disp[i] + n_send[i];
        //fout_debug<<"n_send[i]= "<<n_send[i]<<std::endl;
        for(int j=0; j<n_send[i]; j++){
            val_send[n_send_disp[i]+j] = i*radix + PS::Comm::getRank();
        }
    }
    //for(int i=0; i<n_send_disp[n_proc]; i++) fout_debug<<"i= "<<i<<" val_send[i]="<<val_send[i]<<std::endl;
}

