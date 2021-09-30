#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<unistd.h>
#include<vector>
#include<algorithm>
#include<mpi.h>
#include<particle_simulator.hpp>
#include"comm.hpp"
#include"init.hpp"
#include<mpi-ext.h>


inline double get_wtime(){
    MPI::COMM_WORLD.Barrier();
    return MPI::Wtime()*1e3; //[ms]
    //return MPI::Wtime(); //[sec]
}

class Wtime{
private:
    double cum_;
    double offset_;
    int n_split_;
public:
    double split_;
    double ave_;
    double max_;
    double min_;
    Wtime(){
        max_ = -1.0;
        min_ = 1e8;
        ave_ = cum_ = split_ = 0.0;
        n_split_ = 0;
    }
    void clear(){
        max_ = -1.0;
        min_ = 1e8;
        ave_ = cum_ = split_ = 0.0;
        n_split_ = 0;
    }
    void start(){
        offset_ = get_wtime();
    }
    void stop(){
        n_split_++;
        split_ = get_wtime() - offset_;
        cum_ += split_;
        ave_ = cum_ / n_split_;
        max_ = std::max(max_, split_);
        min_ = std::min(min_, split_);
    }
    void dump(const std::string str, std::ostream & fout=std::cout){
        fout<<"***** "<<str<<" ****"<<std::endl;
        fout<<"ave_="<<ave_<<" max_="<<max_<<" min_="<<min_<<std::endl;
    }
};

std::ofstream fout_debug;

int main(int argc, char *argv[]){
    PS::Initialize(argc, argv);
    char sout_debug[1024];
    sprintf(sout_debug, "./debug_%05d.dat", PS::Comm::getRank());
    fout_debug.open(sout_debug);
    const int dim_comm_max = 6;
    int n_proc_1d_comm[dim_comm_max];
    int dim_comm = 0;
    const int dim_sys_max = 3;
    int n_proc_1d_sys[dim_sys_max];
    int dim_sys = 0;
    int c;

    while((c=getopt(argc,argv,"a:b:c:x:y:z:X:Y:Z:h")) != -1){
        switch(c){
        case 'a':
            n_proc_1d_sys[0] = atoi(optarg);
            dim_sys++;
            break;
        case 'b':
            n_proc_1d_sys[1] = atoi(optarg);
            dim_sys++;
            break;
        case 'c':
            n_proc_1d_sys[2] = atoi(optarg);
            dim_sys++;
            break;
        case 'x':
            n_proc_1d_comm[0] = atoi(optarg);
            dim_comm++;
            break;
        case 'y':
            n_proc_1d_comm[1] = atoi(optarg);
            dim_comm++;
            break;
        case 'z':
            n_proc_1d_comm[2] = atoi(optarg);
            dim_comm++;
            break;
        case 'X':
            n_proc_1d_comm[3] = atoi(optarg);
            dim_comm++;
            break;
        case 'Y':
            n_proc_1d_comm[4] = atoi(optarg);
            dim_comm++;
            break;
        case 'Z':
            n_proc_1d_comm[5] = atoi(optarg);
            dim_comm++;
            break;
        case 'h':
            std::cerr<<"a: n_proc_1d_sys (x)"<<std::endl;
            std::cerr<<"b: n_proc_1d_sys (y)"<<std::endl;
            std::cerr<<"c: n_proc_1d_sys (z)"<<std::endl;
            std::cerr<<"x: n_proc_1d_comm (x0)"<<std::endl;
            std::cerr<<"y: n_proc_1d_comm (x1)"<<std::endl;
            std::cerr<<"z: n_proc_1d_comm (x2)"<<std::endl;
            std::cerr<<"X: n_proc_1d_comm (x3)"<<std::endl;
            std::cerr<<"Y: n_proc_1d_comm (x4)"<<std::endl;
            std::cerr<<"Z: n_proc_1d_comm (x5)"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }
    fout_debug<<"PS::Comm::getRank()= "<<PS::Comm::getRank()<<std::endl;
    for(int i=0; i<dim_sys; i++){
        fout_debug<<"n_proc_1d_sys[i]= "<<n_proc_1d_sys[i]<<"  ";
    }
    fout_debug<<std::endl;
    for(int i=0; i<dim_comm; i++){
        fout_debug<<"n_proc_1d_comm[i]= "<<n_proc_1d_comm[i]<<"  ";
    }
    fout_debug<<std::endl;

    int n_proc = PS::Comm::getNumberOfProc();
    double wtime_begin = get_wtime();
    int n_loop = 5;
    const int min_power = 0;
    //const int min_power = 1;
    //const int min_power = 8;
    int max_power_tmp = 10;
    if(n_proc >= 32768 && max_power_tmp > 12) max_power_tmp = 12;

    const int max_power = max_power_tmp;
    const int min_size = 1<<min_power; // array size. NOT Byte. (Bybe : sizeof(int)*array_size)
    const int max_size = 1<<max_power; // array size. NOT Byte
    const int size_factor = 2;
    MPI::COMM_WORLD.Barrier();
    int * val_send_uni = new int[max_size*n_proc+100];
    int * val_recv_uni = new int[max_size*n_proc+100];
    const int n_neighbor_node_1d = 1; // 1 (3d:27), 2(3d:125)
    const int n_neighbor_node = (n_neighbor_node_1d*2+1)*(n_neighbor_node_1d*2+1)*(n_neighbor_node_1d*2+1);
    int * val_send = new int[max_size*n_proc+size_factor*max_size*n_neighbor_node+100];
    int * val_recv = new int[max_size*n_proc+size_factor*max_size*n_neighbor_node+100];
    int * n_send = new int[n_proc];
    int * n_recv = new int[n_proc];
    int * n_send_disp = new int[n_proc+1];
    int * n_recv_disp = new int[n_proc+1];

    AlltoallSystem<int> a2a_system(n_proc_1d_comm, dim_comm);
    a2a_system.initLoop(1, 20);
    a2a_system.initLoop(max_size, 1); // to assign buffer
    a2a_system.setSizeCrit(max_size);
    Wtime wtime_a2a_normal;
    Wtime wtime_a2a_md;
    Wtime wtime_a2a_i;
    Wtime wtime_a2av_normal;
    Wtime wtime_a2av_ms;
    Wtime wtime_a2av_md;
    Wtime wtime_a2av_i;
    fout_debug<<"wtime_cum="<<get_wtime() - wtime_begin<<std::endl;
    for(int size=min_size; size <= max_size; size *= 2){
    //for(int size=min_size; size <= min_size; size *= 2){
        fout_debug<<"******* size= "<<size<<" size_close= "<<size*size_factor<<" *******"<<std::endl;
        for(int loop=0; loop<n_loop; loop++){
            SetValueUni(val_send_uni, n_proc, size);
            SetValueUnUni(val_send, n_send, n_send_disp, n_proc_1d_sys, n_proc, 
                          size*size_factor, size, dim_sys, n_neighbor_node_1d);

            // normal alltoall
            a2a_system.allToAllNormal(val_send_uni, val_recv_uni, size);
            MPI_Barrier(MPI_COMM_WORLD);
            wtime_a2a_normal.start();
            a2a_system.allToAllNormal(val_send_uni, val_recv_uni, size);
            wtime_a2a_normal.stop();
/*
            fout_debug<<"a2a_normal"<<std::endl;
            for(PS::S32 i=0; i<n_proc*size; i++){
                fout_debug<<"i= "<<i
                          <<" val_send_uni[i]= "<<val_send_uni[i]
                          <<" val_recv_uni[i]= "<<val_recv_uni[i]
                          <<std::endl;
            }
            fout_debug<<std::endl;
*/

            // multi dimensional alltoall
            a2a_system.allToAllMD(val_send_uni, val_recv_uni, size);
            MPI_Barrier(MPI_COMM_WORLD);
            wtime_a2a_md.start();
            a2a_system.allToAllMD(val_send_uni, val_recv_uni, size);
            wtime_a2a_md.stop();


            a2a_system.allToAllI(val_send_uni, val_recv_uni, size);
            MPI_Barrier(MPI_COMM_WORLD);
            wtime_a2a_i.start();
            a2a_system.allToAllI(val_send_uni, val_recv_uni, size);
            wtime_a2a_i.stop();
/*
            fout_debug<<"a2a_md"<<std::endl;
            for(PS::S32 i=0; i<n_proc*size; i++){
                fout_debug<<"i= "<<i
                          <<" val_send_uni[i]= "<<val_send_uni[i]
                          <<" val_recv_uni[i]= "<<val_recv_uni[i]
                          <<std::endl;
            }
            fout_debug<<std::endl;
*/

            // normal alltoall
            a2a_system.allToAllNormal(n_send, n_recv, 1);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(int i=0; i<n_proc; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i]; 
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i]; 
            }
            a2a_system.allToAllVNormal(val_send, n_send, n_send_disp,
                                       val_recv, n_recv, n_recv_disp);
            MPI_Barrier(MPI_COMM_WORLD);
            wtime_a2av_normal.start();
            a2a_system.allToAllVNormal(val_send, n_send, n_send_disp,
                                       val_recv, n_recv, n_recv_disp);
            wtime_a2av_normal.stop();
/*
            fout_debug<<"a2av_normal, n_recv_disp[n_proc]="<<n_recv_disp[n_proc]<<std::endl;
            for(PS::S32 i=0; i<n_recv_disp[n_proc]; i++){
                fout_debug<<"i= "<<i<<" val_recv[i]= "<<val_recv[i]<<std::endl;
            }
            fout_debug<<std::endl;
*/


            a2a_system.allToAllVMD(val_send, n_send, n_send_disp,
                                   val_recv, n_recv, n_recv_disp);
            MPI_Barrier(MPI_COMM_WORLD);
            wtime_a2av_md.start();
            a2a_system.allToAllVMD(val_send, n_send, n_send_disp,
                                   val_recv, n_recv, n_recv_disp);
            wtime_a2av_md.stop();
/*
            fout_debug<<"a2av_md, n_recv_disp[n_proc]="<<n_recv_disp[n_proc]<<std::endl;
            for(PS::S32 i=0; i<n_recv_disp[n_proc]; i++){
                fout_debug<<"i= "<<i<<" val_recv[i]= "<<val_recv[i]<<std::endl;
            }
            fout_debug<<std::endl;
*/


            //std::cerr<<"check e"<<std::endl;
            a2a_system.allToAllNormal(n_send, n_recv, 1);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(int i=0; i<n_proc; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i]; 
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i]; 
            }
            a2a_system.allToAllVMS(val_send, n_send, n_send_disp,
                                   val_recv, n_recv, n_recv_disp, size);
            MPI_Barrier(MPI_COMM_WORLD);
            wtime_a2av_ms.start();
            a2a_system.allToAllVMS(val_send, n_send, n_send_disp,
                                   val_recv, n_recv, n_recv_disp, size);
            wtime_a2av_ms.stop();
/*
            fout_debug<<"a2av_ms, n_recv_disp[n_proc]="<<n_recv_disp[n_proc]<<std::endl;
            for(PS::S32 i=0; i<n_recv_disp[n_proc]; i++){
                fout_debug<<"i= "<<i<<" val_recv[i]= "<<val_recv[i]<<std::endl;
            }
            fout_debug<<std::endl;
*/


            a2a_system.allToAllNormal(n_send, n_recv, 1);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(int i=0; i<n_proc; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i]; 
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i]; 
            }
            a2a_system.allToAllVI(val_send, n_send, n_send_disp,
                                  val_recv, n_recv, n_recv_disp);
            MPI_Barrier(MPI_COMM_WORLD);
            wtime_a2av_i.start();
            a2a_system.allToAllVI(val_send, n_send, n_send_disp,
                                  val_recv, n_recv, n_recv_disp);
            wtime_a2av_i.stop();


        }

        wtime_a2a_normal.dump("wtime_a2a_normal", fout_debug);
        wtime_a2a_normal.clear();
        wtime_a2a_md.dump("wtime_a2a_md", fout_debug);
        wtime_a2a_md.clear();
        wtime_a2a_i.dump("wtime_a2a_i", fout_debug);
        wtime_a2a_i.clear();
        fout_debug<<std::endl;
        wtime_a2av_normal.dump("wtime_a2av_normal", fout_debug);
        wtime_a2av_normal.clear();
        wtime_a2av_md.dump("wtime_a2av_md", fout_debug);
        wtime_a2av_md.clear();
        wtime_a2av_ms.dump("wtime_a2av_ms", fout_debug);
        wtime_a2av_ms.clear();
        wtime_a2av_i.dump("wtime_a2av_i", fout_debug);
        wtime_a2av_i.clear();

        fout_debug<<std::endl;
        fout_debug<<std::endl;
    }

    PS::Finalize();
    return 0;
}
