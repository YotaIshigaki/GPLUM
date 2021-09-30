#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<algorithm>
#include<cassert>
#include<mpi.h>
//#include<particle_simulator.hpp>
//#include"comm.hpp"
//#include"mpi.h"

double GetWtime(){
    int err = MPI_Barrier(MPI_COMM_WORLD);
    return MPI_Wtime(); //[sec]
}

class WTime{
    static const int MAX_SIZE = 100;
    int n_cnt_;
    double split_[MAX_SIZE];
    double wtime_beg_;
    double wtime_end_;
public:
    WTime(): n_cnt_(0), wtime_beg_(0.0), wtime_end_(0.0){}
    void clear(){
        for(int i=0; i<n_cnt_; i++){ split_[i] = 0.0; }
        n_cnt_ = 0;
    }
    void start(){
        wtime_beg_ = GetWtime();
    }
    void stop(){
        wtime_end_ = GetWtime();
        split_[n_cnt_++] = wtime_end_ - wtime_beg_;
    }
    double getMax(){
        assert(n_cnt_ > 0);
        double wtime_max = split_[0];
        for(int i=1; i<n_cnt_; i++){ 
            if(wtime_max < split_[i]) wtime_max = split_[i];
        }
        return wtime_max;
    }
    double getMin(){
        assert(n_cnt_ > 0);
        double wtime_min = split_[0];
        for(int i=1; i<n_cnt_; i++){ 
            if(wtime_min > split_[i]) wtime_min = split_[i];
        }
        return wtime_min;
    }
    double getAvg(){
        assert(n_cnt_ > 0);
        double wtime_ave = 0.0;
        for(int i=0; i<n_cnt_; i++){ 
            wtime_ave += split_[i];
        }
        wtime_ave -= (getMax() + getMin());
        return wtime_ave / (n_cnt_-2);
    }
    double getNCnt(){
        return n_cnt_;
    }
};

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int my_rank, n_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    const int min_exp = 1;
    const int max_exp = 15;
    const int min_size = 1<<min_exp;
    const int max_size = 1<<max_exp;
    const int n_try = 20;
    char * send_buf = new char[max_size];
    char * recv_buf = new char[max_size*n_proc];

#if 0
    WTime time_allgather;
    for(int i=0; i<n_try; i++){
        MPI_Allgather(send_buf, max_size, MPI_CHAR,
                      recv_buf, max_size, MPI_CHAR, MPI_COMM_WORLD);
    }
    for(int size=min_size; size < max_size; size *= 2){
        for(int i=0; i<n_try; i++){
            MPI_Allgather(send_buf, size, MPI_CHAR,
                          recv_buf, size, MPI_CHAR, MPI_COMM_WORLD);
        }
        for(int i=0; i<n_try; i++){
            time_allgather.start();
            MPI_Allgather(send_buf, size, MPI_CHAR,
                          recv_buf, size, MPI_CHAR, MPI_COMM_WORLD);
            time_allgather.stop();
        }
        if(my_rank == 0){
            std::cout<<"n_proc= "<<n_proc
                     <<" size(B)= "<<size*sizeof(char)
                     <<" time_allgather.getAvg()= "<<time_allgather.getAvg()
                     <<" time_allgather.getNCnt()= "<<time_allgather.getNCnt()
                     <<" time_allgather.getMax()= "<<time_allgather.getMax()
                     <<" time_allgather.getMin()= "<<time_allgather.getMin()
                     <<std::endl;
        }
        time_allgather.clear();
    }
#endif

#ifdef GATHER
    WTime time_gather;
    for(int i=0; i<n_try; i++){
        MPI_Gather(send_buf, max_size, MPI_CHAR,
                   recv_buf, max_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    for(int size=min_size; size < max_size; size *= 2){
        for(int i=0; i<n_try; i++){
            MPI_Gather(send_buf, size, MPI_CHAR,
                       recv_buf, size, MPI_CHAR, 0, MPI_COMM_WORLD);
        }
        for(int i=0; i<n_try; i++){
            time_gather.start();
            MPI_Gather(send_buf, size, MPI_CHAR,
                       recv_buf, size, MPI_CHAR, 0, MPI_COMM_WORLD);
            time_gather.stop();
        }
        if(my_rank == 0){
            std::cout<<"n_proc= "<<n_proc
                     <<" size(B)= "<<size*sizeof(char)
                     <<" time_gather.getAvg()= "<<time_gather.getAvg()
                     <<" time_gather.getNCnt()= "<<time_gather.getNCnt()
                     <<" time_gather.getMax()= "<<time_gather.getMax()
                     <<" time_gather.getMin()= "<<time_gather.getMin()
                     <<std::endl;
        }
        time_gather.clear();
    }
#endif
    MPI_Finalize();
}

#if 0
inline double get_wtime(){
    MPI::COMM_WORLD.MPI_Barrier();
    return MPI::Wtime()*1e3; //[ms]
}

int main(int argc, char *argv[]){
    //MPI::Init(argc, argv);
    PS::Initialize(argc, argv);
    int myrank = MPI::COMM_WORLD.Get_rank();
    int n_proc = MPI::COMM_WORLD.Get_size();

    //MPI::COMM_WORLD.Barrier();
    //if(myrank == 0) std::cout<<"n_proc="<<n_proc<<std::endl;

    int * n_send = new int[n_proc];
    int * n_recv = new int[n_proc];
    int * n_send_disp = new int[n_proc+1];
    int * n_recv_disp = new int[n_proc+1];
    const int min_power = 2;
    //const int max_power = 20;
    //const int max_power = 18;
    //const int max_power = 16;
    const int max_power = 15;
    const int min_size = 1<<min_power; // Array size. NOT Byte
    const int max_size = 1<<max_power; // Array size. NOT Byte
    const int n_size = max_power - min_power + 1;

    MPI::COMM_WORLD.Barrier();
    if(myrank == 0) std::cout<<"max_size="<<max_size<<std::endl;
    double * wtime_atav_max = new double[n_size];
    double * wtime_agv_max = new double[n_size];
    double * wtime_ag_max = new double[n_size];
    double * wtime_atav_2d_max = new double[n_size];
    double * wtime_atav_min = new double[n_size];
    double * wtime_agv_min = new double[n_size];
    double * wtime_ag_min = new double[n_size];
    double * wtime_atav_2d_min = new double[n_size];
    double * wtime_atav_ave = new double[n_size];
    double * wtime_agv_ave = new double[n_size];
    double * wtime_ag_ave = new double[n_size];
    double * wtime_atav_2d_ave = new double[n_size];
    for(int i=0; i<n_size; i++){
        wtime_atav_max[i] = wtime_atav_2d_max[i] = wtime_agv_max[i] = wtime_ag_max[i] = -1e6;
        wtime_atav_min[i] = wtime_atav_2d_min[i] = wtime_agv_min[i] = wtime_ag_min[i] = 1e6;
        wtime_atav_ave[i] = wtime_atav_2d_ave[i] = wtime_agv_ave[i] = wtime_ag_ave[i] = 0.0;
    }

    int max_size_array = max_size*n_proc;
    ParticleSimulator::ReallocatableArray<char> val_send;
    ParticleSimulator::ReallocatableArray<char> val_recv;
    val_send.reserve(max_size_array);
    val_recv.reserve(max_size_array);
    val_send.resizeNoInitialize(max_size_array);
    val_recv.resizeNoInitialize(max_size_array);
    for(int i=0; i<max_size_array; i++){
        val_send[i]='a';
    }

    //MPI::COMM_WORLD.Barrier();
    //if(myrank == 0) std::cout<<"check a"<<std::endl;

// initialize
    //Comm<2> comm_2d(MPI::COMM_WORLD);
    PS::CommForAllToAll<char, 2> comm_2d(MPI::COMM_WORLD);
    n_send_disp[0] = n_recv_disp[0] = 0;
    for(int i=0; i<n_proc; i++){
        n_send[i] = max_size;
        n_recv[i] = max_size;
        n_send_disp[i+1] = n_send_disp[i] + n_send[i];
        n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
    }
    MPI::COMM_WORLD.Barrier();
    if(myrank == 0) std::cout<<"check a1"<<std::endl;
    for(int i=0; i<20; i++){
        MPI::COMM_WORLD.Alltoallv(val_send.getPointer(), n_send, n_send_disp, MPI::CHAR,
                                  val_recv.getPointer(), n_recv, n_recv_disp, MPI::CHAR);
    }
    //MPI::COMM_WORLD.Barrier();
    //if(myrank == 0) std::cout<<"check b"<<std::endl;
    for(int i=0; i<20; i++){
        MPI::COMM_WORLD.Allgatherv(val_send.getPointer(), max_size, MPI::CHAR,
                                   val_recv.getPointer(), n_recv, n_recv_disp, MPI::CHAR);
    }
    //MPI::COMM_WORLD.Barrier();
    //if(myrank == 0) std::cout<<"check c"<<std::endl;
    int n_recv_tot_tmp = 0;
    for(int i=0; i<20; i++){
        //comm_2d.AllToAllV(val_send, val_recv, n_recv_tot_tmp, n_send, n_recv);
	//comm_2d.AllToAllV(val_send, val_send_buf, val_recv, n_recv_tot_tmp, n_send, n_recv);
	comm_2d.executeV(val_send, val_recv, n_send, n_recv);
    }
    //MPI::COMM_WORLD.Barrier();
    //if(myrank == 0) std::cout<<"check d"<<std::endl;
// end of initialize

////////////
// alltoallv
////////////
    int n_cnt = 0;
    for(int size=max_size; size >= min_size; size /= 2){
        n_send_disp[0] = n_recv_disp[0] = 0;
        for(int i=0; i<n_proc; i++){
            n_send[i] = size;
            n_recv[i] = size;
            n_send_disp[i+1] = n_send_disp[i] + n_send[i];
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        }

        for(int i=0; i<20; i++){
            MPI::COMM_WORLD.Alltoallv(val_send.getPointer(), n_send, n_send_disp, MPI::CHAR,
                                      val_recv.getPointer(), n_recv, n_recv_disp, MPI::CHAR);
        }
        for(int i=0; i<10; i++){
            double wtime_offset = get_wtime();
            MPI::COMM_WORLD.Alltoallv(val_send.getPointer(), n_send, n_send_disp, MPI::CHAR,
                                      val_recv.getPointer(), n_recv, n_recv_disp, MPI::CHAR);
            double wtime_tmp = get_wtime() - wtime_offset;
            if(wtime_atav_min[n_cnt] > wtime_tmp){
                wtime_atav_min[n_cnt] = wtime_tmp;
            }
            if(wtime_atav_max[n_cnt] < wtime_tmp){
                wtime_atav_max[n_cnt] = wtime_tmp;
            }
            wtime_atav_ave[n_cnt] += wtime_tmp;
        }
        wtime_atav_ave[n_cnt] *= 0.1;
        n_cnt++;
    }
    //if(myrank == 0) std::cout<<"check e"<<std::endl;

////////////
// allgatherv
////////////
    n_cnt = 0;
    //for(unsigned int size=1; size<max_size; size=size<<1){
    //for(unsigned int size=min_size; size<max_size; size=size<<1){
    for(int size=max_size; size >= min_size; size /= 2){
        n_recv_disp[0] = 0;
        for(int i=0; i<n_proc; i++){
            n_recv[i] = size;
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        }
        for(int i=0; i<20; i++){
            MPI::COMM_WORLD.Allgatherv(val_send.getPointer(), size, MPI::CHAR,
                                       val_recv.getPointer(), n_recv, n_recv_disp, MPI::CHAR);
        }
        //double wtime_offset = get_wtime();
        for(int i=0; i<10; i++){
            double wtime_offset = get_wtime();
            MPI::COMM_WORLD.Allgatherv(val_send.getPointer(), size, MPI::CHAR,
                                       val_recv.getPointer(), n_recv, n_recv_disp, MPI::CHAR);
            double wtime_tmp = get_wtime() - wtime_offset;
            if(wtime_agv_min[n_cnt] > wtime_tmp){
                wtime_agv_min[n_cnt] = wtime_tmp;
            }
            if(wtime_agv_max[n_cnt] < wtime_tmp){
                wtime_agv_max[n_cnt] = wtime_tmp;
            }
            wtime_agv_ave[n_cnt] += wtime_tmp;
        }
        wtime_agv_ave[n_cnt] *= 0.1;
        n_cnt++;
    }


////////////
// allgather
////////////
    n_cnt = 0;
    //for(unsigned int size=1; size<max_size; size=size<<1){
    //for(unsigned int size=min_size; size<max_size; size=size<<1){
    for(int size=max_size; size >= min_size; size /= 2){
        n_recv_disp[0] = 0;
        for(int i=0; i<n_proc; i++){
            n_recv[i] = size;
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        }
        for(int i=0; i<20; i++){
            MPI::COMM_WORLD.Allgather(val_send.getPointer(), size, MPI::CHAR,
                                      val_recv.getPointer(), size, MPI::CHAR);
        }
        //double wtime_offset = get_wtime();
        for(int i=0; i<10; i++){
            double wtime_offset = get_wtime();
            MPI::COMM_WORLD.Allgather(val_send.getPointer(), size, MPI::CHAR,
                                      val_recv.getPointer(), size, MPI::CHAR);
            double wtime_tmp = get_wtime() - wtime_offset;
            if(wtime_ag_min[n_cnt] > wtime_tmp){
                wtime_ag_min[n_cnt] = wtime_tmp;
            }
            if(wtime_ag_max[n_cnt] < wtime_tmp){
                wtime_ag_max[n_cnt] = wtime_tmp;
            }
            wtime_ag_ave[n_cnt] += wtime_tmp;
        }
        wtime_ag_ave[n_cnt] *= 0.1;
        n_cnt++;
    }

    //if(myrank == 0) std::cout<<"check f"<<std::endl;

////////////
// alltoallv (multi)
////////////
    n_cnt = 0;
    //for(unsigned int size=1; size<max_size; size=size<<1){
    //for(unsigned int size=min_size; size<max_size; size=size<<1){
    for(int size=max_size; size >= min_size; size /= 2){
        n_send_disp[0] = n_recv_disp[0] = 0;
        for(int i=0; i<n_proc; i++){
            n_send[i] = size;
            n_recv[i] = size;
            n_send_disp[i+1] = n_send_disp[i] + n_send[i];
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        }
        val_send.resizeNoInitialize(n_send_disp[n_proc]);
        val_recv.resizeNoInitialize(n_recv_disp[n_proc]);
        int n_recv_tot = 0;
        for(int i=0; i<20; i++){
            comm_2d.executeV(val_send, val_recv, n_send, n_recv);
        }
        for(int i=0; i<10; i++){
            double wtime_offset = get_wtime();
            comm_2d.executeV(val_send, val_recv, n_send, n_recv);
            double wtime_tmp = get_wtime() - wtime_offset;
            if(wtime_atav_2d_min[n_cnt] > wtime_tmp){
                wtime_atav_2d_min[n_cnt] = wtime_tmp;
            }
            if(wtime_atav_2d_max[n_cnt] < wtime_tmp){
                wtime_atav_2d_max[n_cnt] = wtime_tmp;
            }
            wtime_atav_2d_ave[n_cnt] += wtime_tmp;
        }
        wtime_atav_2d_ave[n_cnt] *= 0.1;
        n_cnt++;
    }
    std::cout<<std::setprecision(15);
    if(myrank == 0){
        n_cnt = 0;
        std::cout<<"# size alltoallv(ave)  alltoallv(max)   alltoallv(min)  ";
        std::cout<<"allgatherv(ave)  allgatherv(max)  allgatherv(min)   ";
        std::cout<<"allgather(ave)  allgather(max)  allgather(min)   ";
        std::cout<<"allgatherv(mod ave)  allgatherv(mod max)  allgatherv(mod min)"<<std::endl;
        int size_tmp = min_size;
        for(int i=n_size-1; i>=0; i--){
            //std::cout<<size_tmp<<"   "<<wtime_atav[i]*0.1<<"   "<<wtime_atav_2d[i]*0.1<<"   "<<wtime_agv[i]*0.1<<std::endl;
            std::cout<<size_tmp
                     <<"   "<<wtime_atav_ave[i]<<"   "<<wtime_atav_max[i]<<"   "<<wtime_atav_min[i]
                     <<"   "<<wtime_agv_ave[i]<<"   "<<wtime_agv_max[i]<<"   "<<wtime_agv_min[i]
                     <<"   "<<wtime_ag_ave[i]<<"   "<<wtime_ag_max[i]<<"   "<<wtime_ag_min[i]
                     <<"   "<<wtime_atav_2d_ave[i]<<"   "<<wtime_atav_2d_max[i]<<"   "<<wtime_atav_2d_min[i]<<std::endl;
            size_tmp *= 2;
        }
    }
    PS::Finalize();
    return 0;
}
#endif
