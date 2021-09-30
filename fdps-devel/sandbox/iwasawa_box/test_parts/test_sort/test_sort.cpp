
#if 0

#include<algorithm>
#include<iostream>
#include<particle_simulator.hpp>
#include"gtest/gtest.h"

int * ARGC;
char *** ARGV;

#if defined(GOOGLE_TEST)

class MPIEnvironment : public ::testing::Environment{
public:
    virtual void SetUp() {
        PS::Initialize(*ARGC, *ARGV);
    }
    virtual void TearDown() {
        PS::Finalize();
    }
    virtual ~MPIEnvironment() {}
};

/*
TEST(test_sample_sort_case, test_sample_sort){
    std::mt19937_64 gen;
    std::uniform_int_distribution<PS::U64> rand_ulong;
    std::uniform_int_distribution<PS::U32> rand_uint;
    std::uniform_real_distribution<PS::F64> rand_double;
    for(auto loop=0; loop<100; loop++){
        const auto n = rand_uint(gen) % 1000000;
        PS::F64 * x = new PS::F64[n];
        for(auto i=0; i<n; i++){
            x[i] = rand_double(gen);
        }        
        PS::SampleSortOmp(x, x+n, n, [](const PS::F64 & l, const PS::F64 & r){return l < r;});
        for(auto i=0; i<n-1; i++){
            EXPECT_LE(x[i], x[i+1]);
        }
        delete [] x;
    }    
}
*/
/*
TEST(test_sample_sort_key_case, test_sample_sort_key){
    std::mt19937_64 gen;
    std::uniform_int_distribution<PS::U64> rand_ulong;
    std::uniform_int_distribution<PS::U32> rand_uint;
    for(auto loop=0; loop<10; loop++){
        const auto n = rand_uint(gen) % 1000000;

#if defined(USE_128BIT_KEY)
        PS::KeyT * x = new PS::KeyT[n];
        for(auto i=0; i<n; i++){
            x[i].hi_ = rand_ulong(gen);
            x[i].lo_ = rand_ulong(gen);
        }
#endif
        PS::SampleSortOmp(x, x+n, n, [](const PS::KeyT & l, const PS::KeyT & r){return l < r;});
        for(auto i=0; i<n-1; i++){
            EXPECT_LE(x[i], x[i+1]);
        }
        delete [] x;
    }
}
*/

/*
TEST(test_mrg_sort_impl_case, test_mrg_sort_impl){
    std::mt19937_64 gen;
    std::uniform_int_distribution<PS::S32> rand_int(-50, 50);
    //std::uniform_int_distribution<PS::S32> rand_int(-100, -10);
    //std::uniform_int_distribution<PS::S32> rand_int(10, 100);

    const auto n0 = 100;
    PS::S32 * x0 = new PS::S32[n0];
    const auto n1 = 50;
    PS::S32 * x1 = new PS::S32[n1];
    PS::S32 * x2 = new PS::S32[n0+n1];
    for(int loop=0; loop<3; loop++){
        std::cerr<<"**************"<<std::endl;
        std::cerr<<"loop= "<<loop<<std::endl;
        if(loop==0){
            for(auto i=0; i<n0; i++){
                x0[i] = rand_int(gen);
            }
            for(auto i=0; i<n1; i++){
                x1[i] = rand_int(gen);
            }
        }
        else if(loop==1){
            for(auto i=0; i<n0; i++){
                x0[i] = 1000+i;
            }
            for(auto i=0; i<n1; i++){
                x1[i] = 100+i;
            }            
        }
        else if(loop==2){
            for(auto i=0; i<n0; i++){
                x0[i] = 100+i;
            }
            for(auto i=0; i<n1; i++){
                x1[i] = 1000+i;
            }            
        }
        std::sort(x0, x0+n0, [](const PS::S32 & l, const PS::S32 & r){return l < r;});
        std::sort(x1, x1+n1, [](const PS::S32 & l, const PS::S32 & r){return l < r;});
        for(int i=0; i<n0; i++)
            std::cerr<<"i= "<<i<<" x0[i]= "<<x0[i]<<std::endl;
        for(int i=0; i<n1; i++)
            std::cerr<<"i= "<<i<<" x1[i]= "<<x1[i]<<std::endl;
        PS::MergeSortOmpImpl2(x0, x0+n0, x1, x1+n1, x2, [](const PS::S32 & l, const PS::S32 & r){return l < r;});
        for(auto i=0; i<n0-1; i++){
            EXPECT_LE(x0[i], x0[i+1]);
        }
        for(auto i=0; i<n1-1; i++){
            EXPECT_LE(x1[i], x1[i+1]);
        }
        for(auto i=0; i<n0+n1-1; i++){
            EXPECT_LE(x2[i], x2[i+1]);
        }
    }
    delete [] x0;
    delete [] x1;
    delete [] x2;
    
}
*/

/*
TEST(test_mrg_sort_impl_case, test_mrg_sort_impl){
    std::mt19937_64 gen;
    std::uniform_int_distribution<PS::S32> rand_int(-500, 500);
    for(int loop=0; loop<10; loop++){
        const auto n = 1000;
        PS::S32 * x = new PS::S32[n];
        for(auto i=0; i<n; i++){
            x[i] = rand_int(gen);
        }        
        PS::MergeSortOmp(x, x+n, n, [](const PS::S32 & l, const PS::S32 & r){return l < r;});
        for(auto i=0; i<n-1; i++){
            std::cerr<<"i= "<<i<<" x[i]= "<<x[i]<<std::endl;
            EXPECT_LE(x[i], x[i+1]);
        }
        delete [] x;
    }
}
*/


TEST(test_sort_performance_case, test_sort_peformance){
    std::mt19937_64 gen;
    std::uniform_int_distribution<PS::U64> rand_ulong;
    std::uniform_int_distribution<PS::U32> rand_uint;
    constexpr int n = 1000000;
    constexpr int n0 = 100;
    double t0, t_mrg, t_smp, t_sir, t_rad;
    //constexpr int n = 100;
    PS::TreeParticle * x_org = new PS::TreeParticle[n];
    for(auto i=0; i<n; i++){
        x_org[i].adr_ptcl_ = i;
        x_org[i].key_.hi_ = rand_ulong(gen);
#if defined(USE_128BIT_KEY)
        x_org[i].key_.lo_ = rand_ulong(gen);
#endif
    }
    PS::RadixSort<PS::KeyT, 8> rs;

    PS::TreeParticle * x_sir = new PS::TreeParticle[n];
    //PS::TreeParticle * x_smp = new PS::TreeParticle[n];
    PS::ReallocatableArray<PS::TreeParticle> x_smp(n, n, PS::MemoryAllocMode::Pool);
    PS::TreeParticle * x_rad = new PS::TreeParticle[n];
    PS::TreeParticle * x_rad_buf = new PS::TreeParticle[n];
    //PS::TreeParticle * x_mrg = new PS::TreeParticle[n];
    PS::ReallocatableArray<PS::TreeParticle> x_mrg(n, n, PS::MemoryAllocMode::Pool);
    for(auto i=0; i<n; i++){
        x_mrg[i] = x_sir[i] = x_smp[i] = x_rad[i] = x_org[i];
    }

    // start measurment

    for(int i=0; i<10; i++){
        std::sort(x_sir+n0, x_sir+n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
        for(int j=0; j<n; j++){
            x_sir[i] = x_org[i];
        }
    }

    t0 = PS::GetWtime();
    PS::MergeSortOmp(x_mrg, n0, n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
    t_mrg = PS::GetWtime() - t0;
    std::cerr<<"t_mrg= "<<t_mrg<<std::endl;
    
    t0 = PS::GetWtime();
    std::sort(x_sir+n0, x_sir+n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
    t_sir = PS::GetWtime() - t0;
    std::cerr<<"t_sir= "<<t_sir<<std::endl;    

    t0 = PS::GetWtime();
    PS::SampleSortOmp(x_smp, n0, n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
    t_smp = PS::GetWtime() - t0;
    std::cerr<<"t_smp= "<<t_smp<<std::endl;
    
    t0 = PS::GetWtime();
    rs.lsdSort(x_rad, x_rad_buf, n0, n-1);
    t_rad = PS::GetWtime() - t0;
    std::cerr<<"t_rad= "<<t_rad<<std::endl;

    for(auto i=0; i<n; i++){
        x_mrg[i] = x_sir[i] = x_smp[i] = x_rad[i] = x_org[i];
    }

    t0 = PS::GetWtime();
    PS::MergeSortOmp(x_mrg, n0, n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
    t_mrg = PS::GetWtime() - t0;
    std::cerr<<"t_mrg= "<<t_mrg<<std::endl;
    
    t0 = PS::GetWtime();
    std::sort(x_sir+n0, x_sir+n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
    t_sir = PS::GetWtime() - t0;
    std::cerr<<"t_sir= "<<t_sir<<std::endl;    

    t0 = PS::GetWtime();
    PS::SampleSortOmp(x_smp, n0, n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
    t_smp = PS::GetWtime() - t0;
    std::cerr<<"t_smp= "<<t_smp<<std::endl;
    
    t0 = PS::GetWtime();
    rs.lsdSort(x_rad, x_rad_buf, n0, n-1);
    t_rad = PS::GetWtime() - t0;
    std::cerr<<"t_rad= "<<t_rad<<std::endl;
    


    
    
    for(auto i=0; i<n0; i++){
        EXPECT_EQ(x_sir[i].getKey(), x_org[i].getKey());
        EXPECT_EQ(x_smp[i].getKey(), x_org[i].getKey());
        EXPECT_EQ(x_mrg[i].getKey(), x_org[i].getKey());
    }
    for(auto i=n0; i<n-1; i++){
        EXPECT_LE(x_sir[i].getKey(), x_sir[i+1].getKey());
        EXPECT_LE(x_smp[i].getKey(), x_smp[i+1].getKey());
        EXPECT_LE(x_mrg[i].getKey(), x_mrg[i+1].getKey());
        //EXPECT_LE(x_rad[i].getKey(), x_rad[i+1].getKey());
    }
    /*    
    std::cerr<<"t_sir= "<<t_sir
             <<" t_smp= "<<t_smp
             <<" t_rad= "<<t_rad
             <<" t_mrg= "<<t_mrg
             <<std::endl;
    */
}

#endif


int main(int argc, char * argv[]){
    ARGC = &argc;
    ARGV = &argv;
    
#if defined(GOOGLE_TEST)
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
    auto ret = RUN_ALL_TESTS();

#else
    PS::Initialize(*ARGC, *ARGV);
#endif

    
    std::mt19937_64 gen;
    std::uniform_int_distribution<PS::U64> rand_ulong;
    std::uniform_int_distribution<PS::U32> rand_uint;
    int n = atoi(argv[1]);
    constexpr int n0 = 0;
    double t0, t_mrg, t_smp, t_sir, t_rad;
    PS::TreeParticle * x_org = new PS::TreeParticle[n];
    for(auto i=0; i<n; i++){
        x_org[i].adr_ptcl_ = i;
        x_org[i].key_.hi_ = rand_ulong(gen);
#if defined(USE_128BIT_KEY)
        x_org[i].key_.lo_ = rand_ulong(gen);
#elif defined(USE_96BIT_KEY)
        x_org[i].key_.lo_ = rand_uint(gen);
#endif
    }
    PS::RadixSort<PS::KeyT, 8> rs;
    PS::TreeParticle * x_sir = new PS::TreeParticle[n];
    PS::ReallocatableArray<PS::TreeParticle> x_smp(n, n, PS::MemoryAllocMode::Pool);
    PS::TreeParticle * x_rad = new PS::TreeParticle[n];
    PS::TreeParticle * x_rad_buf = new PS::TreeParticle[n];
    PS::ReallocatableArray<PS::TreeParticle> x_mrg(n, n, PS::MemoryAllocMode::Pool);
    for(auto i=0; i<n; i++){
        x_mrg[i] = x_sir[i] = x_smp[i] = x_rad[i] = x_org[i];
    }

    std::cerr<<"n= "<<n<<" sizeof(PS::TreeParticle)= "<<sizeof(PS::TreeParticle)<<" n_thread= "<<PS::Comm::getNumberOfThread()<<std::endl;
    for(int loop=0; loop<5; loop++){
        std::cerr<<"******************"<<std::endl;
        t0 = PS::GetWtime();
        std::sort(x_sir+n0, x_sir+n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
        t_sir = PS::GetWtime() - t0;
        std::cerr<<"t_sir= "<<t_sir<<std::endl;

        t0 = PS::GetWtime();
        rs.lsdSort(x_rad, x_rad_buf, n0, n-1);
        t_rad = PS::GetWtime() - t0;
        std::cerr<<"t_rad= "<<t_rad<<std::endl;

        t0 = PS::GetWtime();
        PS::SampleSortOmp(x_smp, n0, n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
        t_smp = PS::GetWtime() - t0;
        std::cerr<<"t_smp= "<<t_smp<<std::endl;
        
        t0 = PS::GetWtime();
        PS::MergeSortOmp(x_mrg, n0, n, [](const PS::TreeParticle & l, const PS::TreeParticle & r){return l.getKey() < r.getKey();});
        t_mrg = PS::GetWtime() - t0;
        std::cerr<<"t_mrg= "<<t_mrg<<std::endl;

        for(auto i=n0; i<n-1; i++){
            assert(x_sir[i].getKey() <= x_sir[i+1].getKey());
            assert(x_rad[i].getKey() <= x_rad[i+1].getKey());
            assert(x_smp[i].getKey() <= x_smp[i+1].getKey());
            assert(x_mrg[i].getKey() <= x_mrg[i+1].getKey());
        }        
        for(auto i=0; i<n; i++){
            x_mrg[i] = x_sir[i] = x_smp[i] = x_rad[i] = x_org[i];
        }
    }
    return 0;
}
#endif

    
#include<iostream>
#include<random>
#include<particle_simulator.hpp>

std::mt19937_64 MT32;
std::mt19937_64 MT64;

struct KIPair{
    PS::KeyT key_;
    PS::U32  id_;
    PS::KeyT getKey() const {
        return key_;
    }
    void dump(std::ostream & fout=std::cout) const {
        fout<<std::oct<<"key_="<<key_<<std::dec<<std::endl;
        fout<<"id_="<<id_<<std::endl;
    }
};


inline PS::KeyT GetRandomKey(){
    PS::KeyT ret;
    ret.hi_ = (PS::U64)MT64();
#if defined (USE_96BIT_KEY)
    ret.lo_ = (PS::U32)MT32();
#elif defined (USE_128BIT_KEY)
    ret.lo_ = (PS::U64)MT64();
#endif
    return ret;
}

inline void CopyOldVal(PS::ReallocatableArray<KIPair> & src, PS::ReallocatableArray<KIPair> & dst){
    for(auto i=0; i<src.size(); i++){
        dst[i] = src[i];
    }
}

int main(int argc, char * argv[]){

    PS::Initialize(argc, argv);

    PS::F64 wtime[3][5]; // 0:radix sort, 1:sample sort, 2:merge sort
    
    const PS::S32 n[] = {1000, 10000, 100000, 1000000, 10000000};
    std::cout<<"sizeof(KIPair)= "<<sizeof(KIPair)<<std::endl;
    PS::ReallocatableArray<KIPair> ki[5];
    PS::ReallocatableArray<KIPair> ki_org[5];
    for(auto i=0; i<5; i++){
        ki[i].resizeNoInitialize(n[i]);
        ki_org[i].resizeNoInitialize(n[i]);
        for(auto j=0; j<n[i]; j++){
            ki_org[i][j].key_ = GetRandomKey();
            ki_org[i][j].id_  = j;
            ki[i][j] = ki_org[i][j];
        }
        CopyOldVal(ki_org[i], ki[i]);
    }

    for(auto i=0; i<5; i++){
        std::cerr<<"i= "<<i<<std::endl;
        for(auto j=0; j<10; j++){
            ki_org[i][j].dump();
        }
    }
    
    PS::ReallocatableArray<KIPair> buf;
    buf.resizeNoInitialize(10000000);
    PS::RadixSort<PS::KeyT, 8> rs;
    for(auto i=0; i<5; i++){
        std::cerr<<"i= "<<i<<std::endl;
        const PS::F64 t0 = PS::GetWtime();
        rs.lsdSort(ki[i].getPointer(), buf.getPointer(), 0, ki[i].size()-1);
        wtime[0][i] = PS::GetWtime() - t0;
        for(auto j=0; j<ki[i].size()-1;  j++){
            assert(ki[i][j].getKey() <= ki[i][j].getKey());
        }
    }

    for(auto i=0; i<5; i++){
        CopyOldVal(ki_org[i], ki[i]);
        std::cerr<<"i= "<<i<<std::endl;
        const PS::F64 t0 = PS::GetWtime();
        SampleSortOmp(ki[i], 0, n[i], [](const auto & l, const auto & r){return l.key_ < r.key_;});
        wtime[1][i] = PS::GetWtime() - t0;
        for(auto j=0; j<ki[i].size()-1;  j++){
            assert(ki[i][j].getKey() <= ki[i][j].getKey());
        }
    }
    
    for(auto i=0; i<5; i++){
        CopyOldVal(ki_org[i], ki[i]);
        std::cerr<<"i= "<<i<<std::endl;
        const PS::F64 t0 = PS::GetWtime();
        MergeSortOmp(ki[i], 0, n[i], [](const auto & l, const auto & r){return l.key_ < r.key_;});
        wtime[2][i] = PS::GetWtime() - t0;
        for(auto j=0; j<ki[i].size()-1;  j++){
            assert(ki[i][j].getKey() <= ki[i][j].getKey());
        }
    }

    for(auto i=0; i<3; i++){
        for(auto j=0; j<5; j++){
            std::cerr<<"i= "<<i<<" j= "<<j<<" wtime[i][j]= "<<wtime[i][j]<<std::endl;
        }
    }

    return 0;
}
