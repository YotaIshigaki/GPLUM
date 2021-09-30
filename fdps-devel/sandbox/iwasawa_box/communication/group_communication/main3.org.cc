#include<iostream>
#include<cmath>
#include<vector>
#include<algorithm>
#include"mpi.h"
extern "C"{
#include<mpi-ext.h>
}

using namespace std;

inline double get_wtime(){
	return MPI::Wtime()*1e3; //[ms]
}

template<int DIM>
void MakeDivision(int np[],   int rank[],
                  const int nproc,   const int myrank){
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


int main(int argc, char *argv[]){
    MPI::Init(argc, argv);

    int myrank = MPI::COMM_WORLD.Get_rank();
    int nproc = MPI::COMM_WORLD.Get_size();

    int fj_dim = 0;
    FJMPI_Topology_get_dimension(&fj_dim);
    std::cout<<"fj_dim="<<fj_dim<<std::endl;
    int fj_dim_x = 0;
    int fj_dim_y = 0;
    int fj_dim_z = 0;
    FJMPI_Topology_get_shape(&fj_dim_x, &fj_dim_y, &fj_dim_z);
    std::cout<<"fj_dim_x="<<fj_dim_x<<", fj_dim_y="<<fj_dim_y<<", fj_dim_z="<<fj_dim_z<<std::endl;
    int fj_rank_x = 0;
    int fj_rank_y = 0;
    int fj_rank_z = 0;
    FJMPI_Topology_rank2x(myrank, &fj_rank_x);
    std::cout<<"fj_rank_x="<<fj_rank_x<<std::endl;
    FJMPI_Topology_rank2xy(myrank, &fj_rank_x, &fj_rank_y);
    std::cout<<"fj_rank_x="<<fj_rank_x<<", fj_rank_y="<<fj_rank_y<<std::endl;
    FJMPI_Topology_rank2xyz(myrank, &fj_rank_x, &fj_rank_y, &fj_rank_z);
    std::cout<<"fj_rank_x="<<fj_rank_x<<", fj_rank_y="<<fj_rank_y<<", fj_rank_z="<<fj_rank_z<<std::endl;
    int fj_rank_a = 0;
    int fj_rank_b = 0;
    int fj_rank_c = 0;
    FJMPI_Topology_sys_rank2xyzabc(myrank, &fj_rank_x, &fj_rank_y, &fj_rank_z, &fj_rank_a, &fj_rank_b, &fj_rank_c);
    std::cout<<"fj_rank_x="<<fj_rank_x<<", fj_rank_y="<<fj_rank_y<<", fj_rank_z="<<fj_rank_z<<", fj_rank_a="<<fj_rank_a<<", fj_rank_b="<<fj_rank_b<<", fj_rank_c="<<fj_rank_c<<std::endl;

    const int DIM = 3;
    int * np = new int[DIM];
    int * rank = new int[DIM];
    MakeDivision<DIM>(np, rank, nproc, myrank);

    std::cout<<"myrank="<<myrank<<", np[0]="<<np[0]<<", np[1]="<<np[1]<<", np[2]="<<np[2]<<std::endl;

    MPI::Intracomm * comm = new MPI::Intracomm[DIM];
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
        comm[d] = MPI::COMM_WORLD.Split(split_key-factor, myrank);
    }

    const int ncnt = 128;
    int * vals_x = new int[ncnt*np[0]];
    int * vals_y = new int[ncnt*np[1]];
    int * vals_z = new int[ncnt*np[2]];
    int * valr_x = new int[ncnt*np[0]];
    int * valr_y = new int[ncnt*np[1]];
    int * valr_z = new int[ncnt*np[2]];
    for(int i=0; i<np[0]*ncnt; i++){
        vals_x[i] = myrank;
        valr_x[i] = 0;
    }
    for(int i=0; i<np[1]*ncnt; i++){
        vals_y[i] = myrank;
        valr_y[i] = 0;
    }
    for(int i=0; i<np[2]*ncnt; i++){
        vals_z[i] = myrank;
        valr_z[i] = 0;
    }

    double Tisend = 0.0;
    double Tisend2 = 0.0;
    double Ta2a = 0.0;
    double Toffset = 0.0;

    MPI::Request * req_send = new MPI::Request[np[0]+np[1]+np[2]];
    MPI::Request * req_recv = new MPI::Request[np[0]+np[1]+np[2]];
    int loop_max = 10000;

    for(int loop=0; loop<100; loop++){
        int nsend = 0;
        for(int i=1; i<np[0]; i++){
            int dst_x = (rank[0] + i)%np[0];
            int src_x = (rank[0] - i + np[0])%np[0];
            int tag_send_x = (rank[0] < dst_x) ? rank[0] : dst_x;
            int tag_recv_x = (rank[0] < src_x) ? rank[0] : src_x;
            req_send[nsend] = comm[0].Isend(vals_x+(ncnt*dst_x), ncnt, MPI::INT, dst_x, tag_send_x);
            req_recv[nsend] = comm[0].Irecv(valr_x+(ncnt*src_x), ncnt, MPI::INT, src_x, tag_recv_x);
            nsend++;
            int dst_y = (rank[1] + i)%np[1];
            int src_y = (rank[1] - i + np[1])%np[1];
            int tag_send_y = (rank[1] < dst_y) ? rank[1] : dst_y;
            int tag_recv_y = (rank[1] < src_y) ? rank[1] : src_y;
            req_send[nsend] = comm[1].Isend(vals_y+(ncnt*dst_y), ncnt, MPI::INT, dst_y, tag_send_y);
            req_recv[nsend] = comm[1].Irecv(valr_y+(ncnt*src_y), ncnt, MPI::INT, src_y, tag_recv_y);
            nsend++;
            int dst_z = (rank[2] + i)%np[2];
            int src_z = (rank[2] - i + np[2])%np[2];
            int tag_send_z = (rank[2] < dst_z) ? rank[2] : dst_z;
            int tag_recv_z = (rank[2] < src_z) ? rank[2] : src_z;
            req_send[nsend] = comm[2].Isend(vals_z+(ncnt*dst_z), ncnt, MPI::INT, dst_z, tag_send_z);
            req_recv[nsend] = comm[2].Irecv(valr_z+(ncnt*src_z), ncnt, MPI::INT, src_z, tag_recv_z);
            nsend++;
        }
        MPI::Request::Waitall(nsend, req_send);
        MPI::Request::Waitall(nsend, req_recv);
    }

    Toffset = get_wtime();
    for(int loop=0; loop<loop_max; loop++){
        int nsend = 0;
        for(int i=1; i<np[0]; i++){
            int dst_x = (rank[0] + i)%np[0];
            int src_x = (rank[0] - i + np[0])%np[0];
            int tag_send_x = (rank[0] < dst_x) ? rank[0] : dst_x;
            int tag_recv_x = (rank[0] < src_x) ? rank[0] : src_x;
            req_send[nsend] = comm[0].Isend(vals_x+(ncnt*dst_x), ncnt, MPI::INT, dst_x, tag_send_x);
            req_recv[nsend] = comm[0].Irecv(valr_x+(ncnt*src_x), ncnt, MPI::INT, src_x, tag_recv_x);
            nsend++;
            int dst_y = (rank[1] + i)%np[1];
            int src_y = (rank[1] - i + np[1])%np[1];
            int tag_send_y = (rank[1] < dst_y) ? rank[1] : dst_y;
            int tag_recv_y = (rank[1] < src_y) ? rank[1] : src_y;
            req_send[nsend] = comm[1].Isend(vals_y+(ncnt*dst_y), ncnt, MPI::INT, dst_y, tag_send_y);
            req_recv[nsend] = comm[1].Irecv(valr_y+(ncnt*src_y), ncnt, MPI::INT, src_y, tag_recv_y);
            nsend++;
            int dst_z = (rank[2] + i)%np[2];
            int src_z = (rank[2] - i + np[2])%np[2];
            int tag_send_z = (rank[2] < dst_z) ? rank[2] : dst_z;
            int tag_recv_z = (rank[2] < src_z) ? rank[2] : src_z;
            req_send[nsend] = comm[2].Isend(vals_z+(ncnt*dst_z), ncnt, MPI::INT, dst_z, tag_send_z);
            req_recv[nsend] = comm[2].Irecv(valr_z+(ncnt*src_z), ncnt, MPI::INT, src_z, tag_recv_z);
            nsend++;
        }
        MPI::Request::Waitall(nsend, req_send);
        MPI::Request::Waitall(nsend, req_recv);
    }
    Tisend2 = get_wtime() - Toffset;


    for(int loop=0; loop<100; loop++){
        int nsend = 0;
        for(int i=1; i<np[0]; i++){
            int dst = (rank[0] + i)%np[0];
            int src = (rank[0] - i + np[0])%np[0];
            int tag_send = (rank[0] < dst) ? rank[0] : dst;
            int tag_recv = (rank[0] < src) ? rank[0] : src;
            req_send[nsend] = comm[0].Isend(vals_x+(ncnt*dst), ncnt, MPI::INT, dst, tag_send);
            req_recv[nsend] = comm[0].Irecv(valr_x+(ncnt*src), ncnt, MPI::INT, src, tag_recv);
            nsend++;
        }
        for(int i=1; i<np[1]; i++){
            int dst = (rank[1] + i)%np[1];
            int src = (rank[1] - i + np[1])%np[1];
            int tag_send = (rank[1] < dst) ? rank[1] : dst;
            int tag_recv = (rank[1] < src) ? rank[1] : src;
            req_send[nsend] = comm[1].Isend(vals_y+(ncnt*dst), ncnt, MPI::INT, dst, tag_send);
            req_recv[nsend] = comm[1].Irecv(valr_y+(ncnt*src), ncnt, MPI::INT, src, tag_recv);
            nsend++;
        }
        for(int i=1; i<np[2]; i++){
            int dst = (rank[2] + i)%np[2];
            int src = (rank[2] - i + np[2])%np[2];
            int tag_send = (rank[2] < dst) ? rank[2] : dst;
            int tag_recv = (rank[2] < src) ? rank[2] : src;
            req_send[nsend] = comm[2].Isend(vals_z+(ncnt*dst), ncnt, MPI::INT, dst, tag_send);
            req_recv[nsend] = comm[2].Irecv(valr_z+(ncnt*src), ncnt, MPI::INT, src, tag_recv);
            nsend++;
        }
        MPI::Request::Waitall(nsend, req_send);
        MPI::Request::Waitall(nsend, req_recv);
    }

    Toffset = get_wtime();
    for(int loop=0; loop<loop_max; loop++){
        int nsend = 0;
        for(int i=1; i<np[0]; i++){
            int dst = (rank[0] + i)%np[0];
            int src = (rank[0] - i + np[0])%np[0];
            int tag_send = (rank[0] < dst) ? rank[0] : dst;
            int tag_recv = (rank[0] < src) ? rank[0] : src;
            req_send[nsend] = comm[0].Isend(vals_x+(ncnt*dst), ncnt, MPI::INT, dst, tag_send);
            req_recv[nsend] = comm[0].Irecv(valr_x+(ncnt*src), ncnt, MPI::INT, src, tag_recv);
            nsend++;
        }
        for(int i=1; i<np[1]; i++){
            int dst = (rank[1] + i)%np[1];
            int src = (rank[1] - i + np[1])%np[1];
            int tag_send = (rank[1] < dst) ? rank[1] : dst;
            int tag_recv = (rank[1] < src) ? rank[1] : src;
            req_send[nsend] = comm[1].Isend(vals_y+(ncnt*dst), ncnt, MPI::INT, dst, tag_send);
            req_recv[nsend] = comm[1].Irecv(valr_y+(ncnt*src), ncnt, MPI::INT, src, tag_recv);
            nsend++;
        }
        for(int i=1; i<np[2]; i++){
            int dst = (rank[2] + i)%np[2];
            int src = (rank[2] - i + np[2])%np[2];
            int tag_send = (rank[2] < dst) ? rank[2] : dst;
            int tag_recv = (rank[2] < src) ? rank[2] : src;
            req_send[nsend] = comm[2].Isend(vals_z+(ncnt*dst), ncnt, MPI::INT, dst, tag_send);
            req_recv[nsend] = comm[2].Irecv(valr_z+(ncnt*src), ncnt, MPI::INT, src, tag_recv);
            nsend++;
        }
        MPI::Request::Waitall(nsend, req_send);
        MPI::Request::Waitall(nsend, req_recv);
    }
    Tisend = get_wtime() - Toffset;

    for(int loop=0; loop<100; loop++){
        comm[0].Alltoall(vals_x, ncnt, MPI::INT, valr_x, ncnt, MPI::INT);
        comm[1].Alltoall(vals_y, ncnt, MPI::INT, valr_y, ncnt, MPI::INT);
        comm[2].Alltoall(vals_z, ncnt, MPI::INT, valr_z, ncnt, MPI::INT);
    }

    Toffset = get_wtime();
    for(int loop=0; loop<loop_max; loop++){
        comm[0].Alltoall(vals_x, ncnt, MPI::INT, valr_x, ncnt, MPI::INT);
        comm[1].Alltoall(vals_y, ncnt, MPI::INT, valr_y, ncnt, MPI::INT);
        comm[2].Alltoall(vals_z, ncnt, MPI::INT, valr_z, ncnt, MPI::INT);
    }
    Ta2a = get_wtime() - Toffset;





/*
    if(myrank == 0){
        std::cout<<"x"<<std::endl;
        for(int i=0; i<np[0]*ncnt; i++){
            std::cout<<"valr_x["<<i<<"]="<<valr_x[i]<<std::endl;
        }
        std::cout<<"y"<<std::endl;
        for(int i=0; i<np[1]*ncnt; i++){
            std::cout<<"valr_y["<<i<<"]="<<valr_y[i]<<std::endl;
        }
        std::cout<<"z"<<std::endl;
        for(int i=0; i<np[2]*ncnt; i++){
            std::cout<<"valr_z["<<i<<"]="<<valr_z[i]<<std::endl;
        }
    }
*/

    std::cout<<"Tisend="<<Tisend/loop_max<<std::endl;
    std::cout<<"Ta2a="<<Ta2a/loop_max<<std::endl;
    std::cout<<"Tisend2="<<Tisend2/loop_max<<std::endl;

    MPI::Finalize();

    return 0;
}
