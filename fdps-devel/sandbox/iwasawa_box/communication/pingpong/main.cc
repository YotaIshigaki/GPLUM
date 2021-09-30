#include<iostream>
#include"mpi.h"

using namespace std;

inline double get_wtime(){
    return MPI::Wtime()*1e3; // [ms]
}

int main(int argc, char *argv[]){

    MPI::Init(argc, argv);
    int myrank = MPI::COMM_WORLD.Get_rank();
    int nproc = MPI::COMM_WORLD.Get_size();

    const int powmax = 21; // 1-1M(=2^20)
    const unsigned int N = 1<<(powmax-1);
    char a[N]; 
    int cnt = 1;
    int loopmax = 10000;
    double t0;
    double T[powmax];
    for(int p=0; p<powmax; p++){
    if(myrank == 0){
	for(int i=0; i<100; i++){
	    MPI::COMM_WORLD.Send(a, cnt, MPI::CHAR, 1, cnt);
	    MPI::COMM_WORLD.Recv(a, cnt, MPI::CHAR, 1, cnt);
	}
    }
    else{
	for(int i=0; i<100; i++){
	    MPI::COMM_WORLD.Recv(a, cnt, MPI::CHAR, 0, cnt);
	    MPI::COMM_WORLD.Send(a, cnt, MPI::CHAR, 0, cnt);
	}
    }
    MPI::COMM_WORLD.Barrier();
    t0 = get_wtime();
    if(myrank == 0){
	for(int i=0; i<loopmax; i++){
	    MPI::COMM_WORLD.Send(a, cnt, MPI::CHAR, 1, cnt);
	    MPI::COMM_WORLD.Recv(a, cnt, MPI::CHAR, 1, cnt);
	}
    }
    else{
	for(int i=0; i<loopmax; i++){
	    MPI::COMM_WORLD.Recv(a, cnt, MPI::CHAR, 0, cnt);
	    MPI::COMM_WORLD.Send(a, cnt, MPI::CHAR, 0, cnt);
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
