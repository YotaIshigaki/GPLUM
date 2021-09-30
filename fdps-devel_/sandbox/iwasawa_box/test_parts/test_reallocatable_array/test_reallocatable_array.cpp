#include<iostream>
#include<unistd.h>
#include<particle_simulator.hpp>

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv, 1000000000);
    PS::ReallocatableArray<PS::F64> f64(10, 10, PS::MemoryAllocMode::Pool);
    SET_CALLER_INFO(f64);
    f64.dump();
#if 0
    //const int n_sml = 5;
    const int n_mid = 10;
    //const int n_big = 100;
    const int num = 1000000;

    PS::ReallocatableArray<PS::F64> * f64;
    try{
	//PS::ReallocatableArray<PS::F64> * f64 = new PS::ReallocatableArray<PS::F64>[num];
	f64 = new PS::ReallocatableArray<PS::F64>[num];
    }
    catch(...){
	std::cerr<<"first"<<std::endl;
    }
    std::cerr<<"second"<<std::endl;
    /*
    f64[0].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[1].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[2].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[3].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[4].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[5].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[6].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[7].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    PS::MemoryPool::dump();

    std::cerr<<"remove 1, 3, 5"<<std::endl;
    f64[1].freeMem();
    f64[3].freeMem();
    f64[5].freeMem();
    PS::MemoryPool::dump();

    std::cerr<<"insert mid (POOL)"<<std::endl;
    f64[8].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    PS::MemoryPool::dump();

    std::cerr<<"insert mid (STACK)"<<std::endl;
    f64[9].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    PS::MemoryPool::dump();

    std::cerr<<"insert mid * 3 (POOL)"<<std::endl;
    f64[10].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[11].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[12].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    PS::MemoryPool::dump();

    std::cerr<<"remove 9 (STACK)"<<std::endl;
    f64[9].freeMem();
    PS::MemoryPool::dump();
    
    std::cerr<<"insert stack->pool->stack"<<std::endl;
    f64[13].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    f64[14].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    f64[15].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    PS::MemoryPool::dump();


    std::cerr<<"insert stack->stack->stack"<<std::endl;
    f64[16].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    f64[17].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    f64[18].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    PS::MemoryPool::dump();

    std::cerr<<"remove 17"<<std::endl;
    f64[17].freeMem();
    PS::MemoryPool::dump();

    std::cerr<<"remove 18"<<std::endl;
    f64[18].freeMem();
    PS::MemoryPool::dump();

    std::cerr<<"clear all pool"<<std::endl;
    PS::MemoryPool::clearAll();
    PS::MemoryPool::dump();
    for(auto i=0; i<num; i++){
	if(f64[i].getPointer() != nullptr){
	    std::cerr<<"i= "<<i<<std::endl;
	}
	assert(f64[i].getPointer() == nullptr);
    }
    */
    std::cerr<<"third"<<std::endl;
    auto t0 = PS::GetWtime();
    for(auto i=0; i<num; i++){
	f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Default);
    }
    auto t1 = PS::GetWtime();
    std::cerr<<"third"<<std::endl;
    /*
    PS::MemoryPool::clearAll();
    auto t2 = PS::GetWtime();
    for(auto i=0; i<num; i++){
	//f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
	f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    }
    auto t3 = PS::GetWtime();

    PS::MemoryPool::clearAll();
    auto t4 = PS::GetWtime();
    for(auto i=0; i<num; i++){
	//f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
	//f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
	f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Default);
    }
    auto t5 = PS::GetWtime();
    PS::MemoryPool::clearAll();
    auto t6 = PS::GetWtime();
    for(auto i=0; i<num; i++){
	//f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Pool);
	f64[i].initialize(n_mid, n_mid, PS::MemoryAllocMode::Stack);
    }
    auto t7 = PS::GetWtime();    
    std::cerr<<"t1-t0= "<<t1-t0
	     <<" t2-t1= "<<t2-t1
	     <<" t3-t2= "<<t3-t2
	     <<" t4-t3= "<<t4-t3
	     <<" t5-t4= "<<t5-t4
	     <<" t6-t5= "<<t6-t5
	     <<" t7-t6= "<<t7-t6
	     <<std::endl;
    f64[0].dump();
    */
    //SET_CALLER_INFO(f64[0].initialize(0, 0, PS::MemoryAllocMode::Stack));
    //SET_CALLER_INFO(PS::ReallocatableArray<PS::F64> a[10]);
    //SET_CALLER_INFO(std::vector<PS::ReallocatableArray<PS::S32>> rank_tmp(12, PS::ReallocatableArray<PS::S32>(PS::MemoryAllocMode::Pool)));
    //SET_CALLER_INFO(printf("hello %c", 10));
    try{
	PS::ReallocatableArray<PS::F64> a;
    }
    catch(int & l){
	std::cerr<<"catch A"<<std::endl;
	l = __LINE__;
    }
    catch(...){
	std::cerr<<"catch B"<<std::endl;
    }
    /*
    std::cerr<<"remove big"<<std::endl;
    f64_a.freeMem();
    PS::MemoryPool::dump();
    
    std::cerr<<"initialize big"<<std::endl;
    f64_e.initialize(n_big, n_big, PS::MemoryAllocMode::Pool);
    PS::MemoryPool::dump();
    
    std::cerr<<"remove big"<<std::endl;
    f64_c.freeMem();
    PS::MemoryPool::dump();
    
    std::cerr<<"insert mid"<<std::endl;
    PS::ReallocatableArray<PS::F64> f64_g(n_mid, n_mid, PS::MemoryAllocMode::Pool);
    PS::MemoryPool::dump();

    std::cerr<<"remove mid"<<std::endl;
    f64_g.freeMem();
    PS::MemoryPool::dump();

    std::cerr<<"insert big"<<std::endl;
    PS::ReallocatableArray<PS::F64> f64_h(n_big, n_big, PS::MemoryAllocMode::Pool);
    PS::MemoryPool::dump();
    */
#endif    
    PS::Finalize();
}
