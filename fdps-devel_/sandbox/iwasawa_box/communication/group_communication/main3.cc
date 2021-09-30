#include<iostream>
#include<cmath>
#include<vector>
#include<algorithm>
#include"mpi.h"
extern "C"{
#include<mpi-ext.h>
}

/*
2D
2   5   8  11
1   4   7  10
0   3   6   9
0 node
send
X: (6)0->6, 0->7, 0->8,   (9)0->9, 0->10, 0->11
Y: (0)0->0, 0->3,   (1)0->1, 0->4,   (2)0->2, 0->5
recv
X: 6->0, 6->1, 6->2,   9->0, 9->1, 9->2 
Y: 0->0, 0->3,  1->0, 1->3,  2->0, 2->3
send2
X: 0->0, 1->0, 2->0,   0->3, 1->3, 2->3
y: 6->0, 9->0, 6->1, 9->1, 6->2, 9->2 
recv2
X: (3)3->0, 4->0, 5->0, 
Y: (1)7->0, 10->0,  (2)8->0, 11->0, 

1 node
send
X: (7)1->6, 1->7, 1->8,   (10)1->9, 1->10, 1->11
Y: (0)1->0, 1->3,   (1)1->1, 1->4,   (2)1->2, 1->5
recv
X: 7->0, 7->1, 7->2,   10->0, 10->1, 10->2 
Y: 0->1, 0->4,  1->1, 1->4,  2->1, 2->4
3 node
send
X: (6)0->6, 0->7, 0->8,   (9)0->9, 0->10, 0->11
Y: (3)3->0, 3->3,   (4)3->1, 3->4,  (5)3->2, 3->5
recv
X: 6->3, 6->4, 6->5,   9->3, 9->4, 9->5
Y: 3->0, 3->3,   4->0, 4->3,  5->0, 5->3

6 node
send
X: (0)6->0, 6->1, 6->2,   (3) 6->3, 6->4, 6->5
Y: (6)6->6, 6->9,  (7)6->7, 6->10,  (8)6->8, 6->11
recv
X: 0->6, 0->7, 0->8,   3->6, 3->7, 3->8
Y: 6->6, 6->9,  7->6, 7->9,  8->6, 8->9

3D
6*4*3 = 72

9  21  33  45  57  69
6  18  30  42  54  66
3  15  27  39  51  63
0  12  24  36  48  60


0 node
24: (0,24), (0.25), (0.26), ..(0.35) nsend = nz*ny = 12
36: (0,36), (0.37), (0.38), ..(0.47) nsend = nz*ny = 12
3:  (0,48), (0,49), (0.50), ..(0.
6:
9:
1:
2:

0 node 
X
to 48 node 0->48, 0->49, 0->50,   0->51, 0->54, 0->57
to 60 node 0->60, 0->63, 0->66, 0->69
Y
to 3 node  0->27, 0->39
to 6 node  0->30, 0->42
Z
to 1 node  0->1, 
to 3 node  0->2

recv
X
from 24->0, 24->1, 24->2,   24->3, 24->4, 
from 36->0, 36->1, 36->2,   36->3, 36->4, 
Y

Z
 */

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
    //std::cerr<<"check 1"<<std::endl;
    int myrank = MPI::COMM_WORLD.Get_rank();
    //std::cerr<<"check 2"<<std::endl;
    int nproc = MPI::COMM_WORLD.Get_size();
    //std::cerr<<"check 3"<<std::endl;
/*
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
*/


    const int DIM = 2;
    int * np = new int[DIM];
    int * rank = new int[DIM];
    MakeDivision<DIM>(np, rank, nproc, myrank);
    if(myrank == 0){
        for(int d=0; d<DIM; d++){
            std::cout<<"np["<<d<<"]="<<np[d]<<"   ";
        }
        std::cout<<std::endl;

    }

    const int ncnt = 128;
    int * vals = new int[ncnt*nproc];
    int * valr = new int[ncnt*nproc];
    for(int i=0; i<nproc; i++){
        for(int j=0; j<ncnt; j++){
            vals[i*ncnt+j] = (myrank+1)*1000+(i+1);
            valr[i*ncnt+j] = 0;
        }
    }
    int dim_x[DIM];
    for(int i=0; i<DIM; i++){
        dim_x[i] = np[0]/DIM;
        if( i < np[0]%DIM ){
            dim_x[i]++;
        }
    }
    int dim_x_disp[DIM+1];
    dim_x_disp[0] = 0;
    for(int i=0; i<DIM; i++){
        dim_x_disp[i+1] = dim_x_disp[i] + dim_x[i];
    }
    if(myrank == 0){
        for(int i=0; i<DIM; i++){
            std::cerr<<"dim_x["<<i<<"]="<<dim_x[i]<<"   ";
        }
        std::cerr<<std::endl;
    }
    int * id_start_x = new int[DIM];
    id_start_x[0] = 0;
    for(int i=0; i<DIM; i++){
        id_start_x[i] = dim_x_disp[i]*nproc/np[0];
    }
    if(myrank == 0){
        for(int i=0; i<DIM; i++){
            std::cerr<<"id_start_x["<<i<<"]="<<id_start_x[i]<<"   ";
        }
        std::cerr<<std::endl;
    }

    int mydomain_x = 0;
    for(int i=0; i<DIM-1; i++){
        //if(rank[0] < dim_x_disp[mydomain_x+1]) break;
        if(myrank < id_start_x[mydomain_x+1]) break;
        mydomain_x++;
    }
    std::cerr<<"rank[0]="<<rank[0]<<", mydomain_x="<<mydomain_x<<std::endl;

    int nsend = 0;
    int dst_domain_x = (mydomain_x != DIM-1) ? mydomain_x+1 : 0;

    int * vals_buf[DIM];
    int * valr_buf[DIM];
    for(int i=0; i<DIM; i++){
        vals_buf[i] = new int[nproc*ncnt];
        valr_buf[i] = new int[nproc*ncnt];
    }
    for(int i=0; i<DIM; i++){
        for(int j=0; j<nproc/np[0]
        vals_buf[i][j] = vals[id_start_x[dst_domain_x]+j];
    }
    for(int ix=0; ix<dim_x[dst_domain_x]; ix++){
        req_send[nsend] = MPI::COMM_WORLD.Isend(vals_buf, ncnt*np[1], MPI::INT, dst, tag_send);
        nsend++;
    }


#if 0
    int nsend = 0;
    for(int ix=0; ix<dim_x[]; ix++){
        int i2x = 0;
        while(dim_x_disp[i2x] < rank[0]){
            i2x++;
        }
        req_send[nsend] = MPI::COMM_WORLD.Isend(vals, ncnt*np[1], MPI::INT, dst, tag_send);
    }

    for(int d=0; d<DIM; d++){
        int dst = start_x[]
        req_send[nsend] = MPI::COMM_WORLD.Isend(vals, ncnt*np[1], MPI::INT, dst, tag_send);
        //req_send[nsend] = MPI::COMM_WORLD.Isend(vals, ncnt*np[1], MPI::INT, dst, tag_send);
    }



/*
    int dim_x = np[0]/DIM;
    if(myrank < np[0]%DIM) dim_x++;
    for(int ix=0; ix<dim_x; ix++){
        MPI::COMM_WORLD.Isend
    }
*/
    for(int ix=0; ix<dim_x; ix++){
        int dst = offset_x;
        int src = dst;
        req_send[nsend] = MPI::COMM_WORLD.Isend(vals, ncnt*np[1], MPI::INT, dst, tag_send);
    }
    for(int iy=0; iy<np[1]; iy++){
        req_send[nsend] = MPI::COMM_WORLD.Isend(vals, ncnt, MPI::INT, dst, tag_send);
            //req_send[nsend] = comm[0].Isend(vals_x+(ncnt*dst_x), ncnt, MPI::INT, dst_x, tag_send_x);
    }
#endif

#if 0
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

#endif



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
/*
    std::cout<<"Tisend="<<Tisend/loop_max<<std::endl;
    std::cout<<"Ta2a="<<Ta2a/loop_max<<std::endl;
    std::cout<<"Tisend2="<<Tisend2/loop_max<<std::endl;
*/

    MPI::Finalize();

    return 0;
}
