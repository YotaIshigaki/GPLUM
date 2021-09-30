
#include<iostream>
#include<random>
#include"particle_simulator.hpp"
//#include"comm.hpp"

//namespace PS = ParticleSimulator;

int main(int argc, char *argv[]){
    PS::Initialize(argc, argv);
    int target_rank = atoi(argv[1]);
    if(PS::Comm::getRank() == 0){
	std::cout<<"target_rank="<<target_rank<<std::endl;
    }
    const int DIM_COMM = 2;
    PS::CommForAllToAll<double, DIM_COMM> comm_double(MPI_COMM_WORLD);
    int rank = PS::Comm::getRank();
    int n_proc = PS::Comm::getNumberOfProc();
    int * n_send = new int[n_proc];
    int * n_recv = new int[n_proc];
    std::mt19937 mt(PS::Comm::getRank());
    std::uniform_int_distribution<int> dist(0, 10);
    int n_send_tot = 0;
    for(int i=0; i<n_proc; i++){
	n_send[i] = dist(mt);
        n_send_tot += n_send[i];
    }
    int n_send_offset = dist(mt);
    int n_recv_offset = dist(mt);
    
    int cnt = 0;
    double radix = 100000;
    PS::ReallocatableArray<double> val_send;
    val_send.reserve(n_send_tot+n_send_offset);
    for(int i=0; i<n_send_offset; i++){
        val_send.push_back(-1.0);
    }
    if(PS::Comm::getRank() == target_rank){
        std::cout<<"n_send_offset= "<<n_send_offset<<std::endl;
    }
    val_send.resizeNoInitialize(n_send_offset);
    for(int i=0; i<n_proc; i++){
        if(PS::Comm::getRank() == target_rank){
            std::cout<<"send rank= "<<i<<" n_send[i]= "<<n_send[i]<<std::endl;
        }
	for(int j=0; j<n_send[i]; j++, cnt++){
	    val_send.push_back(radix * (rank+1) + (i+1));
	    if(PS::Comm::getRank() == target_rank){
	        std::cout<<"val_send[cnt+n_send_offset]="<<val_send[cnt+n_send_offset]<<std::endl;
	    }
	}
    }
    std::cout<<std::endl;
    
    PS::ReallocatableArray<double> val_recv;
    val_recv.resizeNoInitialize(n_recv_offset);
    for(int i=0; i<n_recv_offset; i++){
        val_recv.push_back(-10.0);
    }
    comm_double.executeV(val_send, val_recv, n_send, n_recv, n_send_offset, n_recv_offset);
    cnt = 0;
    if(PS::Comm::getRank() == target_rank){
        std::cout<<"n_recv_offset= "<<n_recv_offset<<std::endl;
	for(int i=0; i<n_proc; i++){
            std::cout<<"recv_rank= "<<i<<" n_recv[i]= "<<n_recv[i]<<std::endl;
	    for(int j=0; j<n_recv[i]; j++, cnt++){
                std::cout<<"val_recv[cnt+n_recv_offset]= "<<val_recv[cnt+n_recv_offset]<<std::endl;
	    }
	}
    }
    
    PS::Finalize();
    return 0;
}
