#include<iostream>
#include<cstdlib>
#include<cassert>
#include"mpi.h"

extern "C"{
    #include<athread.h>
    #include<cpe_func.h>
    void SLAVE_FUN(CopyIndirect)(void *);
    void SLAVE_FUN(CopyIndirectInverse)(void *);
    void SLAVE_FUN(CopyDirect)(void *);
    void SLAVE_FUN(CopyStride)(void *);
    void SLAVE_FUN(CopyStrideByHand)(void *);
    //void SLAVE_FUN(CopyIndirectInverse2)(void *);    
}

static inline unsigned long rpcc(){
    unsigned long time;
    asm("rtc %0" : "=r"(time) : );
    return time;
}

class Moment{
public:
    double mass;
    double pos[3];
    double ort0[6];
    double ort1[6];
};

class TreeCell{
public:
    int x0;
    int x1;
    int adr;
    int x3;
    Moment mom;
};

class SuperParticle{
public:
    double mass;
    double pos[3];
};

int main(){
    MPI_Init(0,0);
    athread_init();
    unsigned long arg[7];
    double wtime_st, wtime_ed, wtime;
    const int n_total=1000000;
    TreeCell * tc = new TreeCell[n_total];
    SuperParticle * sp_src = new SuperParticle[n_total];
    SuperParticle * sp_dst = new SuperParticle[n_total];
    int * adr = new int[n_total];
    for(int i=0; i<n_total; i++){
        adr[i] = rand()%n_total;
        //adr[i] = (i+5000)%n_total;
        //adr[i] = i;
        //tc[i].adr = (i+5000)%n_total;
        //tc[i].adr = i;
        tc[i].adr = rand()%n_total;
        //tc[i].mom.mass = 0.125*((double)i)+1.0;
        tc[i].mom.mass = ((double)i);
        tc[i].mom.pos[0] = 1.0;
        tc[i].mom.pos[1] = 2.0;
        tc[i].mom.pos[2] = 3.0;
        sp_src[i].mass   = tc[i].mom.mass;
        sp_src[i].pos[0] = tc[i].mom.pos[0];
        sp_src[i].pos[1] = tc[i].mom.pos[1];
        sp_src[i].pos[2] = tc[i].mom.pos[2];
        sp_dst[i].mass   = -1.0*tc[i].mom.mass;
        sp_dst[i].pos[0] = -1.0*tc[i].mom.pos[0];
        sp_dst[i].pos[1] = -1.0*tc[i].mom.pos[1];
        sp_dst[i].pos[2] = -1.0*tc[i].mom.pos[2];
    }

#if 1
    ///////////////////////
    //INDIRECT COPY 
    std::cerr<<"sp_dst[0].mass(a)= "<<sp_dst[0].mass<<std::endl;
    std::cerr<<"sp_dst[1].mass(a)= "<<sp_dst[1].mass<<std::endl;
    std::cerr<<"sp_dst[2].mass(a)= "<<sp_dst[2].mass<<std::endl;
    //unsigned long arg[5];
    arg[0] = (unsigned long)(n_total);
    arg[1] = (unsigned long)(&sp_src[0]);
    arg[2] = (unsigned long)(&sp_dst[0]);
    arg[3] = (unsigned long)(sizeof(SuperParticle));
    arg[4] = (unsigned long)(adr);
    wtime_st = MPI_Wtime();
    __real_athread_spawn((void*)slave_CopyIndirect, arg);
    athread_join();
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(INDIRECT CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(INDIRECT CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;    
    std::cerr<<"adr[0]= "<<adr[0]<<" sp_dst[adr[0]].mass(B)= "<<sp_dst[adr[0]].mass<<std::endl;
    std::cerr<<"adr[1]= "<<adr[1]<<" sp_dst[adr[1]].mass(B)= "<<sp_dst[adr[1]].mass<<std::endl;
    std::cerr<<"adr[2]= "<<adr[2]<<" sp_dst[adr[2]].mass(B)= "<<sp_dst[adr[2]].mass<<std::endl;
    /*
    for(int i=0; i<n_total; i++){
        if(sp_dst[adr[i]].mass != sp_src[i].mass){
            std::cerr<<"i= "<<i<<" adr[i]= "<<adr[i]<<std::endl;
            std::cerr<<"sp_dst[adr[i]].mass= "<<sp_dst[adr[i]].mass<<" sp_src[i].mass= "<<sp_src[i].mass<<std::endl;
        }
        assert(sp_dst[adr[i]].mass == sp_src[i].mass);
    }
    */
    wtime_st = MPI_Wtime();    
    for(int i=0; i<n_total; i++){
        sp_dst[adr[i]] = sp_src[i];
    }
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(INDIRECT MPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(INDIRECT MPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;    
    //INDIRECT COPY 
    ///////////////////////    
#endif
    
#if 1
    ///////////////////////
    //INDIRECT COPY INVERSE
    //unsigned long arg[5];
    arg[0] = (unsigned long)(n_total);
    arg[1] = (unsigned long)(&sp_src[0]);
    arg[2] = (unsigned long)(&sp_dst[0]);
    arg[3] = (unsigned long)(sizeof(SuperParticle));
    arg[4] = (unsigned long)(adr);
    wtime_st = MPI_Wtime();
    __real_athread_spawn((void*)slave_CopyIndirectInverse, arg);
    athread_join();
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(INDIRECT INVERSE CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(INDIRECT INVERSE CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    
    wtime_st = MPI_Wtime();
    for(int i=0; i<n_total; i++){
        sp_dst[i] = sp_src[adr[i]];
    }
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(INDIRECT INVERSE MPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(INDIRECT INVERSE MPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;    
    /*
    for(int i=0; i<n_total; i++){
        if(sp_dst[i].mass != sp_src[adr[i]].mass){
            std::cerr<<"i= "<<i<<" adr[i]= "<<adr[i]<<std::endl;
            std::cerr<<"sp_dst[i].mass= "<<sp_dst[i].mass<<" sp_src[adr[i]].mass= "<<sp_src[adr[i]].mass<<std::endl;
        }
        assert(sp_src[adr[i]].mass == sp_dst[i].mass);
    }
    */
    //INDIRECT COPY INVERSE
    ///////////////////////
#endif

#if 1
    ///////////////////////
    //DIRECT_COPY
    arg[0] = (unsigned long)(n_total);
    arg[1] = (unsigned long)( &sp_src[0] );
    arg[2] = (unsigned long)( &sp_dst[0] );
    arg[3] = (unsigned long)( sizeof(SuperParticle) );
    wtime_st = MPI_Wtime();
    __real_athread_spawn((void*)slave_CopyDirect, arg);
    athread_join();
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(DIRECT CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(DIRECT CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    wtime_st = MPI_Wtime();
    for(int i=0; i<n_total; i++){
        sp_dst[i] = sp_src[i];
    }
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(DIRECT MPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(DIRECT MPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;    
    //DIRECT_COPY
    ///////////////////////
#endif    

#if 1
    ///////////////////////
    //STRIDE COPY
    //unsigned long arg[5];
    arg[0] = (unsigned long)(n_total);
    arg[1] = (unsigned long)( ((int*)&tc[0])+4 );
    arg[2] = (unsigned long)( &sp_dst[0] );
    arg[3] = (unsigned long)( sizeof(SuperParticle) );
    arg[4] = (unsigned long)( sizeof(int)*4 + sizeof(double)*12 );
    wtime_st = MPI_Wtime();
    __real_athread_spawn((void*)slave_CopyStride, arg);
    athread_join();
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(STRIDE, stride=112B, data=32B CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(STRIDE, stride=112B, data=32B) CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;

    arg[0] = (unsigned long)(n_total);
    arg[1] = (unsigned long)( ((int*)&tc[0])+4 );
    arg[2] = (unsigned long)( &sp_dst[0] );
    arg[3] = (unsigned long)( sizeof(SuperParticle) );
    arg[4] = (unsigned long)( sizeof(int)*4 + sizeof(double)*12 );
    wtime_st = MPI_Wtime();
    __real_athread_spawn((void*)slave_CopyStrideByHand, arg);
    athread_join();
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(STRIDE(by hand), stride=112B, data=32B CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(STRIDE(by hand), stride=112B, data=32B CPE)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    

    wtime_st = MPI_Wtime();
    for(int i=0; i<n_total; i++){
        sp_dst[i].mass = tc[i].mom.mass;
        sp_dst[i].pos[0]  = tc[i].mom.pos[0];
        sp_dst[i].pos[1]  = tc[i].mom.pos[1];
        sp_dst[i].pos[2]  = tc[i].mom.pos[2];
    }
    wtime_ed = MPI_Wtime();
    std::cerr<<"band width(STRIDE MPE, stride=112B, data=32B)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;
    std::cout<<"band width(STRIDE MPE, stride=112B, data=32B)= "<<sizeof(SuperParticle)*n_total / (wtime_ed-wtime_st) / 1e9 <<"[GB/s]"<<" wtime= "<<wtime_ed-wtime_st<<std::endl;    
    /*
    std::cerr<<"sp_dst[0].mass(b)= "<<sp_dst[0].mass<<std::endl;
    std::cerr<<"sp_dst[1].mass(b)= "<<sp_dst[1].mass<<std::endl;
    std::cerr<<"sp_dst[2].mass(b)= "<<sp_dst[2].mass<<std::endl;
    for(int i=0; i<n_total; i++){
        assert(sp_dst[i].mass == tc[i].mom.mass);
        assert(sp_dst[i].pos[0] == tc[i].mom.pos[0]);
        assert(sp_dst[i].pos[1] == tc[i].mom.pos[1]);
        assert(sp_dst[i].pos[2] == tc[i].mom.pos[2]);
    }
    */
    //STRIDE COPY
    ///////////////////////
#endif

#if 0
    ///////////////////////
    //COPY INDIRECT INVERSE
    //unsigned long arg[7];
    arg[0] = (unsigned long)(n_total);
    arg[1] = (unsigned long)(&sp_src[0]);
    arg[2] = (unsigned long)(&sp_dst[0]);
    arg[3] = (unsigned long)(sizeof(sp_src[0]));
    arg[4] = (unsigned long)(((int*)&tc[0])+2); // offset
    arg[5] = (unsigned long)(sizeof(int)); // type of address
    arg[6] = (unsigned long)(sizeof(tc[0])-sizeof(int)); // stride
    double wtime_st = MPI_Wtime();
    __real_athread_spawn((void*)slave_CopyIndirectInverse2, arg);
    athread_join();
    double wtime_ed = MPI_Wtime();
    double wtime = wtime_ed - wtime_st;
    std::cerr<<"band width= "<<sizeof(SuperParticle)*n_total / wtime / 1e9 <<"[GB/s]"<<std::endl;        
    std::cerr<<"tc[0].adr= "<<adr[0]<<" sp_dst[0].mass= "<<sp_dst[0].mass<<" sp_src[tc[0].adr].mass= "<<sp_src[tc[0].adr].mass<<std::endl;
    std::cerr<<"tc[1].adr= "<<adr[1]<<" sp_dst[1].mass= "<<sp_dst[1].mass<<" sp_src[tc[1].adr].mass= "<<sp_src[tc[1].adr].mass<<std::endl;
    std::cerr<<"tc[2].adr= "<<adr[2]<<" sp_dst[2].mass= "<<sp_dst[2].mass<<" sp_src[tc[2].adr].mass= "<<sp_src[tc[2].adr].mass<<std::endl;
    for(int i=0; i<n_total; i++){
        int adr = tc[i].adr;
        if(sp_dst[i].mass   != sp_src[adr].mass){
            std::cerr<<"i= "<<i<<" adr= "<<adr<<std::endl;
            std::cerr<<"sp_dst[i].mass= "<<sp_dst[i].mass<<" sp_src[adr].mass= "<<sp_src[adr].mass<<std::endl;
        }
        assert(sp_dst[i].mass   == sp_src[adr].mass);
        assert(sp_dst[i].pos[0] == sp_src[adr].pos[0]);
        assert(sp_dst[i].pos[1] == sp_src[adr].pos[1]);
        assert(sp_dst[i].pos[2] == sp_src[adr].pos[2]);
    }
    //COPY INDIRECT INVERSE STRIDE    
    ///////////////////////
#endif
    
    athread_halt();

    return 0;
    
}
