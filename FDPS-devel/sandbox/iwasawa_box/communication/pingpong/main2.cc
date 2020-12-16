#include<iostream>
#include<cassert>
#include"mpi.h"

using namespace std;

inline double get_wtime(){
	return MPI::Wtime()*1e3; // [ms]
}

int main(int argc, char *argv[]){

    MPI::Init(argc, argv);
    int myrank = MPI::COMM_WORLD.Get_rank();
    int nproc = MPI::COMM_WORLD.Get_size();
    //assert(nproc == 5);
    //const int powmax=15; // 1-16K(=2^14)
    const int powmax = 21; // 1-1M(=2^20)
    const unsigned int N = 1<<(powmax-1);
    char a[N];
    int cnt = 1;
    int loopmax = 10000;
    double t0;
    double T[powmax];
    MPI::Request * req_send = new MPI::Request[nproc];
    MPI::Request * req_recv = new MPI::Request[nproc];
    for(int p=0; p<powmax; p++){
        if(myrank == 0){
            for(int i=0; i<100; i++){
                req_send[0] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 1, 1);
                //req_send[1] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 2, 2);
                //req_send[2] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 3, 3);
                //req_send[3] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 4, 4);
                //MPI::Request::Waitall(4, req_send);
                MPI::Request::Waitall(1, req_send);
            }
        }
        else{
            for(int i=0; i<100; i++){
                req_recv[0] = MPI::COMM_WORLD.Irecv(a, cnt, MPI::CHAR, 0, myrank);
                MPI::Request::Waitall(1, req_recv);
            }
        }
        MPI::COMM_WORLD.Barrier();
        t0 = get_wtime();
        if(myrank == 0){
            for(int i=0; i<loopmax; i++){
                req_send[0] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 1, 1);
                //req_send[1] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 2, 2);
                //req_send[2] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 3, 3);
                //req_send[3] = MPI::COMM_WORLD.Isend(a, cnt, MPI::CHAR, 4, 4);
                //MPI::Request::Waitall(4, req_send);
                MPI::Request::Waitall(1, req_send);
            }
        }
        else{
            for(int i=0; i<loopmax; i++){
                req_recv[0] = MPI::COMM_WORLD.Irecv(a, cnt, MPI::CHAR, 0, myrank);
                MPI::Request::Waitall(1, req_recv);
            }
        }
        MPI::COMM_WORLD.Barrier();
        T[p] = get_wtime() - t0;
        T[p] /= loopmax;
        cnt *= 2;
    }

    if(myrank == 0){
        cnt = 1;
        for(int p=0; p<powmax; p++){
            std::cout<<cnt<<"   "<<T[p]<<std::endl;
            cnt *= 2;
        }
    }

    MPI::Finalize();

    return 0;
}
