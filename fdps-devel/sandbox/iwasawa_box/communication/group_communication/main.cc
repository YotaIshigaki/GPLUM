
/*
np=512, 4byte
1D
Tcom(Ncom=512, 4byte) = 1.4e-6 * 512 = 7.2e-4 [sec]
2D
Nsend = 4*16*32 + 4*32*16
Ncom = 32+16 = 48
Tcom = 1.8e-6*32 + 4.1e-6*16 = 1.23e-4 [sec]
3D
Nsend = 4*64*8 + 4*64*8 + 4*64*8
Ncom = 8+8+8 = 24
Tcom(Ncom=24, 256byte) = 4.2e-6 * 24 = 1.0e-4 [sec]

 */

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
/*
    Comm<int, 1> comm1d(MPI::COMM_WORLD);
    Comm<int, 2> comm2d(MPI::COMM_WORLD);
    Comm<int, 3> comm3d(MPI::COMM_WORLD);
    Comm<int, 4> comm4d(MPI::COMM_WORLD);
    Comm<int, 5> comm5d(MPI::COMM_WORLD);
    Comm<int, 6> comm6d(MPI::COMM_WORLD);
*/

    Comm<1> comm1d(MPI::COMM_WORLD);
    Comm<2> comm2d(MPI::COMM_WORLD);
    Comm<3> comm3d(MPI::COMM_WORLD);
    Comm<4> comm4d(MPI::COMM_WORLD);
    Comm<5> comm5d(MPI::COMM_WORLD);
    Comm<6> comm6d(MPI::COMM_WORLD);

    //int n_ag = 13; // all gather 13 means 4B*(1-4096)
    //int n_ar = 13; // all reduce 13 means 4B*(1-4096)
    int n_a2a = 8; // alltoall 8 means 4B*(1-128)
    //int n_a2a = 6; // alltoall 8 means 4B*(1-32)

    int N = n_a2a;

    int Ndim = 6+1;
    double ** Txd = new double*[Ndim];
    for(int i=0; i<Ndim; i++){
        Txd[i] = new double[N];
    }

    //int ** val32_send = new int*[N];
    //int ** val32_recv = new int*[N];
    PS_INT ** val32_send = new PS_INT*[N];
    PS_INT ** val32_send_buf = new PS_INT*[N];
    PS_INT ** val32_send_buf2 = new PS_INT*[N];
    PS_INT ** val32_recv = new PS_INT*[N];
    int * cnt = new int[N];
    cnt[0] = 1;
    for(int i=0; i<N-1; i++){
        cnt[i+1] = cnt[i]*2;
    }

    for(int i=0; i<N; i++){
        //val32_send[i] = new PS_INT[cnt[i]];
        val32_send[i] = new PS_INT[cnt[i]*nproc];
        val32_send_buf[i] = new PS_INT[cnt[i]*nproc];
        val32_send_buf2[i] = new PS_INT[cnt[i]*nproc];
        val32_recv[i] = new PS_INT[cnt[i]*nproc];
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<nproc; j++){
            for(int k=0; k<cnt[i]; k++){
                //val32_send[i][k] = myrank+1;
                val32_send[i][j*cnt[i]+k] = myrank+1;
                val32_recv[i][j*cnt[i]+k] = 0;
            }
        }
    }

    //PS_INT * val32_send_buf = new PS_INT[cnt[N-1]*nproc];

    int loopmax = 10000;
    //int loopmax = 1000;
    //int loopmax = 1;
    double t0, t1;

// alltoall old version
    for(int i=0; i<N; i++){
        for(int d=0; d<Ndim; d++){
            if(d==0){
                for(int j=0; j<100; j++){
                    MPI::COMM_WORLD.Alltoall(val32_send[i], cnt[i], MPI_PS_INT, val32_recv[i], cnt[i], MPI_PS_INT);
                }
            }
            if(d==0){
                for(int j=0; j<100; j++){
                    MPI::COMM_WORLD.Alltoall(val32_send[i], cnt[i], MPI_PS_INT, val32_recv[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    MPI::COMM_WORLD.Alltoall(val32_send[i], cnt[i], MPI_PS_INT, val32_recv[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==1){
                for(int j=0; j<100; j++){
                    comm1d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm1d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==2){
                for(int j=0; j<100; j++){
                    comm2d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm2d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==3){
                for(int j=0; j<100; j++){
                    comm3d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm3d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==4){
                for(int j=0; j<100; j++){
                    comm4d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm4d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==5){
                for(int j=0; j<100; j++){
                    comm5d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm5d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==6){
                for(int j=0; j<100; j++){
                    comm6d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm6d.AllToAll(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            MPI::COMM_WORLD.Barrier();
            Txd[d][i] = (t1 - t0)/loopmax;
            //Txd[d][i] = (get_wtime() - t0)/loopmax;
        }
    }

#ifdef BIT32
    std::cout<<"32bit"<<std::endl;
#elif BIT64
    std::cout<<"64bit"<<std::endl;
#endif
#ifdef __OPENMP__
    std::cout<<"using OPENMP"<<std::endl;
#else
    std::cout<<"not using OPENMP"<<std::endl;
#endif
    std::cout<<"nproc="<<nproc<<std::endl;
    std::cout<<"loopmax="<<loopmax<<std::endl;
    std::cout<<"alltoall(old)"<<std::endl;
    for(int i=0; i<N; i++){
        std::cout<<cnt[i];
        for(int d=0; d<Ndim; d++){
            std::cout<<"   "<<Txd[d][i];
        }
        std::cout<<std::endl;
    }


// alltoall sendrecv version
    for(int i=0; i<N; i++){
        for(int j=0; j<nproc; j++){
            for(int k=0; k<cnt[i]; k++){
                val32_send[i][k] = myrank + 1;
                val32_recv[i][j*cnt[i]+k] = 0;
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int d=0; d<Ndim; d++){
            if(d==0){
                for(int j=0; j<100; j++){
                    AllToAll_s(val32_send[i], val32_recv[i], cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
            }
            if(d==0){
                for(int j=0; j<100; j++){
                    AllToAll_s(val32_send[i], val32_recv[i], cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    AllToAll_s(val32_send[i], val32_recv[i], cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
                t1 = get_wtime(); 
            }
            else if(d==1){
                for(int j=0; j<100; j++){
                    comm1d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm1d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==2){
                for(int j=0; j<100; j++){
                    comm2d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm2d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==3){
                for(int j=0; j<100; j++){
                    comm3d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm3d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==4){
                for(int j=0; j<100; j++){
                    comm4d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm4d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==5){
                for(int j=0; j<100; j++){
                    comm5d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm5d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==6){
                for(int j=0; j<100; j++){
                    comm6d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm6d.AllToAll4(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            MPI::COMM_WORLD.Barrier();
            Txd[d][i] = (t1 - t0)/loopmax;
        }
    }

#ifdef BIT32
    std::cout<<"32bit"<<std::endl;
#elif BIT64
    std::cout<<"64bit"<<std::endl;
#endif
#ifdef __OPENMP__
    std::cout<<"using OPENMP"<<std::endl;
#else
    std::cout<<"not using OPENMP"<<std::endl;
#endif
    std::cout<<"nproc="<<nproc<<std::endl;
    std::cout<<"loopmax="<<loopmax<<std::endl;
    std::cout<<"alltoall(new sendrecv)"<<std::endl;
    for(int i=0; i<N; i++){
        std::cout<<cnt[i];
        for(int d=0; d<Ndim; d++){
            std::cout<<"   "<<Txd[d][i];
        }
        std::cout<<std::endl;
    }


// alltoall isend irecv
    MPI::Request * req_send = new MPI::Request[nproc];
    MPI::Request * req_recv = new MPI::Request[nproc];
    for(int i=0; i<N; i++){
        for(int j=0; j<nproc; j++){
            for(int k=0; k<cnt[i]; k++){
                val32_send[i][k] = myrank + 1;
                val32_recv[i][j*cnt[i]+k] = 0;
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int d=0; d<Ndim; d++){
            if(d==0){
                for(int j=0; j<100; j++){
                    AllToAll_i(val32_send[i], val32_recv[i], req_send, req_recv, cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
            }
            if(d==0){
                for(int j=0; j<100; j++){
                    AllToAll_i(val32_send[i], val32_recv[i], req_send, req_recv, cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    AllToAll_i(val32_send[i], val32_recv[i], req_send, req_recv, cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
                t1 = get_wtime(); 
            }
            else if(d==1){
                for(int j=0; j<100; j++){
                    comm1d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm1d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==2){
                for(int j=0; j<100; j++){
                    comm2d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm2d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==3){
                for(int j=0; j<100; j++){
                    comm3d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm3d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==4){
                for(int j=0; j<100; j++){
                    comm4d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm4d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==5){
                for(int j=0; j<100; j++){
                    comm5d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm5d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==6){
                for(int j=0; j<100; j++){
                    comm6d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm6d.AllToAll3(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            MPI::COMM_WORLD.Barrier();
            Txd[d][i] = (t1 - t0)/loopmax;
        }
    }

#ifdef BIT32
    std::cout<<"32bit"<<std::endl;
#elif BIT64
    std::cout<<"64bit"<<std::endl;
#endif
#ifdef __OPENMP__
    std::cout<<"using OPENMP"<<std::endl;
#else
    std::cout<<"not using OPENMP"<<std::endl;
#endif
    std::cout<<"nproc="<<nproc<<std::endl;
    std::cout<<"loopmax="<<loopmax<<std::endl;
    std::cout<<"alltoall(new isend)"<<std::endl;
    for(int i=0; i<N; i++){
        std::cout<<cnt[i];
        for(int d=0; d<Ndim; d++){
            std::cout<<"   "<<Txd[d][i];
        }
        std::cout<<std::endl;
    }



// alltoall 1 dring
    PS_INT * val32_large_buf = new PS_INT[cnt[N-1]*nproc*Ndim];
    //for(int i=0; i<N; i++){
    for(int i=0; i<5; i++){
        for(int j=0; j<nproc; j++){
            for(int k=0; k<cnt[i]; k++){
                val32_send[i][k] = myrank + 1;
                val32_recv[i][j*cnt[i]+k] = 0;
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int d=0; d<Ndim; d++){
            if(d==0){
                for(int j=0; j<100; j++){
                    MPI::COMM_WORLD.Alltoall(val32_send[i], cnt[i], MPI_PS_INT, val32_recv[i], cnt[i], MPI_PS_INT);
                }
            }
            if(d==0){
                for(int j=0; j<100; j++){
                    AllToAll_1dring(val32_send[i], val32_recv[i], val32_large_buf, cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    AllToAll_1dring(val32_send[i], val32_recv[i], val32_large_buf, cnt[i], MPI_PS_INT, MPI::COMM_WORLD);
                }
                t1 = get_wtime(); 
            }
            else if(d==1){
                for(int j=0; j<100; j++){
                    comm1d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm1d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==2){
                for(int j=0; j<100; j++){
                    comm2d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm2d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==3){
                for(int j=0; j<100; j++){
                    comm3d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm3d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==4){
                for(int j=0; j<100; j++){
                    comm4d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm4d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==5){
                for(int j=0; j<100; j++){
                    comm5d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm5d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            else if(d==6){
                for(int j=0; j<100; j++){
                    comm6d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime();
                for(int j=0; j<loopmax; j++){
                    comm6d.AllToAll2(val32_send[i], val32_recv[i], val32_send_buf[i], val32_large_buf, cnt[i], MPI_PS_INT);
                }
                t1 = get_wtime();
            }
            MPI::COMM_WORLD.Barrier();
            Txd[d][i] = (t1 - t0)/loopmax;
        }
    }

#ifdef BIT32
    std::cout<<"32bit"<<std::endl;
#elif BIT64
    std::cout<<"64bit"<<std::endl;
#endif
#ifdef __OPENMP__
    std::cout<<"using OPENMP"<<std::endl;
#else
    std::cout<<"not using OPENMP"<<std::endl;
#endif
    std::cout<<"nproc="<<nproc<<std::endl;
    std::cout<<"alltoall(new 1dring)"<<std::endl;

    for(int i=0; i<N; i++){
        std::cout<<cnt[i];
        for(int d=0; d<Ndim; d++){
            std::cout<<"   "<<Txd[d][i];
        }
        std::cout<<std::endl;
    }


#ifdef ALLTOALLV
// Alltoallv
    int nrecv_tot = 0;
    int ** nsend = new int*[N];
    int ** nsend_buf = new int*[N];
    int ** nrecv = new int*[N];
    int ** nsend_disp = new int*[N];
    int ** nsend_disp_buf = new int*[N];
    int ** nrecv_disp = new int*[N];
    int ** nsend_1d = new int*[N];
    int ** nrecv_1d = new int*[N];
    int ** nsend_disp_1d = new int*[N];
    int ** nrecv_disp_1d = new int*[N];
    for(int i=0; i<N; i++){
        nsend[i] = new int[nproc];
        nsend_buf[i] = new int[nproc];
        nrecv[i] = new int[nproc];
        nsend_disp[i] = new int[nproc+1];
        nsend_disp_buf[i] = new int[nproc+1];
        nrecv_disp[i] = new int[nproc+1];
        nsend_1d[i] = new int[nproc];
        nrecv_1d[i] = new int[nproc];
        nsend_disp_1d[i] = new int[nproc+1];
        nrecv_disp_1d[i] = new int[nproc+1];
    }

    for(int i=0; i<N; i++){
        nsend_disp[i][0] = nrecv_disp[i][0] = 0;
        for(int j=0; j<nproc; j++){
            nsend[i][j] = nrecv[i][j] = cnt[i];
            nsend_disp[i][j+1] = nsend_disp[i][j] + nsend[i][j];
            nrecv_disp[i][j+1] = nrecv_disp[i][j] + nrecv[i][j];
            for(int k=0; k<cnt[i]; k++){
                val32_send[i][k] = myrank+1;
                val32_recv[i][j*cnt[i]+k] = 0;
            }
        }
    }



    for(int i=0; i<N; i++){
        for(int d=0; d<Ndim; d++){
            if(d==0){
                for(int j=0; j<100; j++)
                    MPI::COMM_WORLD.Alltoallv(val32_send[i], nsend[i], nsend_disp[i], MPI_PS_INT, 
                                              val32_recv[i], nrecv[i], nrecv_disp[i], MPI_PS_INT);
            }
            if(d==0){
                for(int j=0; j<100; j++){
                    MPI::COMM_WORLD.Alltoallv(val32_send[i], nsend[i], nsend_disp[i], MPI_PS_INT, 
                                              val32_recv[i], nrecv[i], nrecv_disp[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    MPI::COMM_WORLD.Alltoallv(val32_send[i], nsend[i], nsend_disp[i], MPI_PS_INT, 
                                              val32_recv[i], nrecv[i], nrecv_disp[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==1){
                for(int j=0; j<100; j++){
                    comm1d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i],
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm1d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i],
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==2){
                for(int j=0; j<100; j++){
                    comm2d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i],
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm2d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i],
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==3){
                for(int j=0; j<100; j++){
                    comm3d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm3d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==4){
                for(int j=0; j<100; j++){
                    comm4d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm4d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==5){
                for(int j=0; j<100; j++){
                    comm5d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm5d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            else if(d==6){
                for(int j=0; j<100; j++){
                    comm6d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                MPI::COMM_WORLD.Barrier();
                t0 = get_wtime(); 
                for(int j=0; j<loopmax; j++){
                    comm6d.AllToAllV(val32_send[i], val32_send_buf[i], val32_recv[i], nrecv_tot, nsend[i], nsend_buf[i], nrecv[i], nsend_disp[i], nsend_disp_buf[i], nrecv_disp[i], 
                                     nsend_1d[i], nrecv_1d[i], nsend_disp_1d[i], nrecv_disp_1d[i], MPI_PS_INT);
                }
                t1 = get_wtime(); 
            }
            MPI::COMM_WORLD.Barrier();
            Txd[d][i] = (t1 - t0)/loopmax;
        }

    }

#if 0
    if(myrank == 0){
        for(int i=0; i<N; i++){
           //for(int i=0; i<1; i++){
            for(int j=0; j<nproc*cnt[i]; j++){
                std::cout<<"val32_recv[i][j]="<<val32_recv[i][j]<<std::endl;
            }
            std::cout<<std::endl;
        }
    }
#endif

    std::cout<<std::endl;
    std::cout<<"alltoallv"<<std::endl;
    for(int i=0; i<N; i++){
        std::cout<<cnt[i];
        for(int d=0; d<Ndim; d++){
            std::cout<<"   "<<Txd[d][i];
        }
        std::cout<<std::endl;
    }
#endif //ALLTOALLV

#if 0
// Allreduce
    for(int i=0; i<N; i++){
        for(int j=0; j<nproc; j++){
            for(int k=0; k<cnt[i]; k++){
                val32_send[i][k] = myrank + 1;
                val32_recv[i][j*cnt[i]+k] = 0;
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int d=0; d<Ndim; d++){
            if(d==0)
                for(int j=0; j<100; j++)
                    MPI::COMM_WORLD.Allreduce(val32_send[i], val32_recv[i], cnt[i], MPI_PS_INT, MPI::SUM);
            MPI::COMM_WORLD.Barrier();
            t0 = get_wtime(); 
            if(d==0)
                for(int j=0; j<loopmax; j++)
                    MPI::COMM_WORLD.Allreduce(val32_send[i], val32_recv[i], cnt[i], MPI_PS_INT, MPI::SUM);
            else if(d==1)
                for(int j=0; j<loopmax; j++)
                    comm1d.AllReduce(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI::SUM, MPI_PS_INT);
            else if(d==2)
                for(int j=0; j<loopmax; j++)
                    comm2d.AllReduce(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI::SUM, MPI_PS_INT);
            else if(d==3)
                for(int j=0; j<loopmax; j++)
                    comm3d.AllReduce(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI::SUM, MPI_PS_INT);
            else if(d==4)
                for(int j=0; j<loopmax; j++)
                    comm4d.AllReduce(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI::SUM, MPI_PS_INT);
            else if(d==5)
                for(int j=0; j<loopmax; j++)
                    comm5d.AllReduce(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI::SUM, MPI_PS_INT);
            else if(d==6)
                for(int j=0; j<loopmax; j++)
                    comm6d.AllReduce(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI::SUM, MPI_PS_INT);
            MPI::COMM_WORLD.Barrier();
            Txd[d][i] = (get_wtime() - t0)/loopmax;
/*
            if(myrank==0){
                std::cout<<"i="<<i<<", d="<<d<<std::endl;
                for(int ii=0; ii<cnt[i]; ii++){
                    std::cout<<"val32_recv[i][ii]="<<val32_recv[i][ii]<<std::endl;
                }
            }
*/
        }
    }

    //if(myrank==0){
        std::cout<<std::endl;
        std::cout<<"allreduce"<<std::endl;
        for(int i=0; i<N; i++){
            std::cout<<cnt[i];
            for(int d=0; d<Ndim; d++){
                std::cout<<"   "<<Txd[d][i];
            }
            std::cout<<std::endl;
        }
        //}

#endif // allreduce
#if 0

// Allgather
    for(int i=0; i<N; i++){
        for(int j=0; j<nproc; j++){
            for(int k=0; k<cnt[i]; k++){
                val32_send[i][k] = myrank + 1;
                val32_recv[i][j*cnt[i]+k] = 0;
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int d=0; d<Ndim; d++){
            if(d==0)
                for(int j=0; j<100; j++)
                    MPI::COMM_WORLD.Allgather(val32_send[i], cnt[i], MPI_PS_INT, val32_recv[i], cnt[i], MPI_PS_INT);
            MPI::COMM_WORLD.Barrier();
            t0 = get_wtime(); 
            if(d==0)
                for(int j=0; j<loopmax; j++)
                    MPI::COMM_WORLD.Allgather(val32_send[i], cnt[i], MPI_PS_INT, val32_recv[i], cnt[i], MPI_PS_INT);
            else if(d==1)
                for(int j=0; j<loopmax; j++)
                    comm1d.AllGather(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
            else if(d==2)
                for(int j=0; j<loopmax; j++)
                    comm2d.AllGather(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
            else if(d==3)
                for(int j=0; j<loopmax; j++)
                    comm3d.AllGather(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
            else if(d==4)
                for(int j=0; j<loopmax; j++)
                    comm4d.AllGather(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
            else if(d==5)
                for(int j=0; j<loopmax; j++)
                    comm5d.AllGather(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);
            else if(d==6)
                for(int j=0; j<loopmax; j++)
                    comm6d.AllGather(val32_send[i], val32_recv[i], val32_send_buf[i], cnt[i], MPI_PS_INT);

            MPI::COMM_WORLD.Barrier();
            Txd[d][i] = (get_wtime() - t0)/loopmax;
/*
            if(myrank==0){
                std::cout<<"i="<<i<<", d="<<d<<std::endl;
                for(int ii=0; ii<cnt[i]*nproc; ii++){
                    std::cout<<"val32_recv[i][ii]="<<val32_recv[i][ii]<<std::endl;
                }
            }
*/
        }
    }

    //if(myrank==0){
    std::cout<<std::endl;
    std::cout<<"allgather"<<std::endl;
    for(int i=0; i<N; i++){
        std::cout<<cnt[i];
        for(int d=0; d<Ndim; d++){
            std::cout<<"   "<<Txd[d][i];
        }
        std::cout<<std::endl;
    }
    //}
#endif // allgather




/*
    for(int i=0; i<N; i++){
        delete [] val32_send[i];
        delete [] val32_recv[i];

        delete [] nsend[i];
        delete [] nrecv[i];
        delete [] nsend_disp[i];
        delete [] nrecv_disp[i];
        delete [] nsend_1d[i];
        delete [] nrecv_1d[i];
        delete [] nsend_disp_1d[i];
        delete [] nrecv_disp_1d[i];
    }
    delete [] val32_send;
    delete [] val32_recv;

    delete [] nsend;
    delete [] nrecv;
    delete [] nsend_disp;
    delete [] nrecv_disp;

    delete [] nsend_1d;
    delete [] nrecv_1d;
    delete [] nsend_disp_1d;
    delete [] nrecv_disp_1d;

    delete [] cnt;
    for(int i=0; i<Ndim; i++){
        delete [] Txd[i];
    }
    delete [] Txd;
*/

    MPI::Finalize();
    return 0;
}
