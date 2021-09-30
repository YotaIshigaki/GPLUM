#include<iostream>
#include"mpi.h"
#include"comm.hpp"


using namespace std;

inline double get_wtime(){
	return MPI::Wtime()*1e3; // [ms]
}


#ifdef BIT32
typedef int PS_INT;
MPI::Datatype MPI_PS_INT = MPI::INT;
#elif BIT64
typedef long long int PS_INT;
MPI::Datatype  MPI_PS_INT = MPI::LONG;
#endif

int main(int argc, char *argv[]){
    MPI::Init(argc, argv);
    int myrank = MPI::COMM_WORLD.Get_rank();
    int nproc = MPI::COMM_WORLD.Get_size();

    int ncnt = 8;
    int * val_send = new int[nproc*ncnt];
    int * val_recv = new int[nproc*ncnt];
    int * val_buf = new int[nproc*ncnt];
    for(int i=0; i<nproc; i++){
        for(int j=0; j<ncnt; j++){
            val_send[i*ncnt+j] = (myrank+1)*10000+j;
        }
    }

    MPI::Request * req_send = new MPI::Request[nproc];
    MPI::Request * req_recv = new MPI::Request[nproc];

    //AllToAll_1dring(val_send, val_recv, val_buf, ncnt, MPI::INT, MPI::COMM_WORLD);
    //AllToAll_i(val_send, val_recv, req_send, req_recv, ncnt, MPI::INT, MPI::COMM_WORLD);
    AllToAll_s(val_send, val_recv, ncnt, MPI::INT, MPI::COMM_WORLD);

    if(myrank == 1){
        for(int i=0; i<nproc*ncnt; i++){
            std::cout<<"i="<<i<<", val_recv="<<val_recv[i]<<std::endl;
        }
    }

    MPI::Finalize();
    return 0;
}
