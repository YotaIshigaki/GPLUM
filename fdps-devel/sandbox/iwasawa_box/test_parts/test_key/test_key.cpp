#define USE_96BIT_KEY
//#define USE_128BIT_KEY
//#define USE_STD_SORT

#include<iostream>
#include<unistd.h>
#include<bitset>
#include<particle_simulator.hpp>
#include"gtest/gtest.h"



int * ARGC;
char *** ARGV;


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


TEST(test_dd_prriodic_case, test_periodic_dd){
    PS::F64 full_len = 1.0;
    PS::F64vec center = 0.5;
    PS::F64 half_len = full_len * 0.5;
    PS::MortonKey::initialize(half_len, center);
    PS::U64 n_limit = (((PS::U64)1)<<((PS::U64)PS::TREE_LEVEL_LIMIT));
    std::cerr<<"n_limit= "<<n_limit<<std::endl;
    PS::F64 dx = full_len / n_limit;
    std::cerr<<"dx= "<<dx<<std::endl;
    std::mt19937_64 gen;
    std::uniform_int_distribution<int> rand_int(0, 7);
    PS::F64vec pos;
    PS::KeyT key;
    /*    
    //PS::F64vec pos = center + PS::F64vec(dx*0.9);
    PS::F64vec pos = center + PS::F64vec(half_len - dx*0.5);
    PS::KeyT key = PS::MortonKey::getKey(pos);

    std::cerr<<"pos= "<<pos<<std::endl;
    std::cerr<<"key= "<<std::bitset<sizeof(PS::U64)*8>(key.hi_)
             <<" "<<std::bitset<sizeof(PS::U64)*8>(key.lo_)
             <<std::endl;    
    std::cerr<<"key= ";
    key.dump(std::cerr);
    for(auto i=0; i<=PS::TREE_LEVEL_LIMIT; i++){
        std::cerr<<"i= "<<i
                 <<" PS::MortonKey::getCellID(i, key)= "<<PS::MortonKey::getCellID(i, key)
                 <<" PS::MortonKey::getPosTreeCell(i, key)= "<<PS::MortonKey::getPosTreeCell(i, key)
                 <<std::endl;
    }
    PS::F64vec box_center = center + PS::F64vec(dx*0.9);
    PS::F64vec box_hlen = PS::F64vec(dx*0.5);
    PS::F64ort box(box_center-box_hlen, box_center+box_hlen);
    std::cerr<<"box= "<<box
             <<" PS::MortonKey::getCorrespondingTreeCell(box)= "<<PS::MortonKey::getCorrespondingTreeCell(box)
             <<std::endl;;
    */

    PS::S32 shift_id[PS::TREE_LEVEL_LIMIT];
    const PS::S32 n = 10000;
    for(PS::S32 loop=0; loop<n; loop++){
        //if(loop%(loop/100+1)==0) std::cerr<<"loop= "<<loop<<std::endl;
        if(loop%(n/100+1)==0) std::cerr<<"loop= "<<loop<<std::endl;
        pos = center;
        shift_id[0] = 0;
        PS::F64 len_tmp = half_len;
        for(PS::S32 i=1; i<PS::TREE_LEVEL_LIMIT; i++){
            shift_id[i] = rand_int(gen);
            pos += PS::SHIFT_CENTER[shift_id[i]]*len_tmp;
            len_tmp *= 0.5;
        }
        key = PS::MortonKey::getKey(pos);
        for(PS::S32 i=0; i<PS::TREE_LEVEL_LIMIT; i++){
            EXPECT_EQ(shift_id[i],  PS::MortonKey::getCellID(i, key));
        }
    }

}

TEST(test_shift_case, test_shift){
    std::mt19937_64 gen;
    std::uniform_int_distribution<PS::U64> rand_ulong;
    std::uniform_int_distribution<PS::U32> rand_uint;
    for(auto i=0; i<1; i++){
#if defined(USE_128BIT_KEY)
        PS::KeyT key(rand_ulong(gen), rand_ulong(gen));
        for(int i=0; i<128; i++){
            PS::KeyT key_new = key >> i;
            std::cerr<<"key= "<<std::bitset<sizeof(PS::U64)*8>(key_new.hi_)
                     <<" "<<std::bitset<sizeof(PS::U64)*8>(key_new.lo_)
                     <<std::endl;
        }

        std::cerr<<std::endl;
        for(int i=0; i<128; i++){
            PS::KeyT key_new = key << i;
            std::cerr<<"key= "<<std::bitset<sizeof(PS::U64)*8>(key_new.hi_)
                     <<" "<<std::bitset<sizeof(PS::U64)*8>(key_new.lo_)
                     <<std::endl;
        }

#elif defined(USE_96BIT_KEY)
        PS::KeyT key(rand_ulong(gen), rand_uint(gen));
        std::cerr<<"key= "<<std::bitset<sizeof(PS::U64)*8>(key.hi_)
                     <<" "<<std::bitset<sizeof(PS::U32)*8>(key.lo_)
                     <<std::endl;        
        for(int i=0; i<96; i++){
            PS::KeyT key_new = key >> i;
            std::cerr<<"key_new= "<<std::bitset<sizeof(PS::U64)*8>(key_new.hi_)
                     <<" "<<std::bitset<sizeof(PS::U32)*8>(key_new.lo_)
                     <<std::endl;
        }
        std::cerr<<std::endl;
        for(int i=0; i<96; i++){
            PS::KeyT key_new = key << i;
            std::cerr<<"key= "<<std::bitset<sizeof(PS::U64)*8>(key_new.hi_)
                     <<" "<<std::bitset<sizeof(PS::U32)*8>(key_new.lo_)
                     <<std::endl;
        }
#endif
    }
}

int main(int argc, char * argv[]){
    ARGC = &argc;
    ARGV = &argv;
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
    auto ret = RUN_ALL_TESTS();
    
    return ret;
}

/*

int main(){
    std::cerr<<std::setprecision(15);
    //PS::F64 full_len = 1.0;
    //PS::F64vec center = 0.5;
    PS::F64 full_len = 8.0*atan(1.0);
    PS::F64vec center = 4.0*atan(1.0);
    PS::F64 half_len = full_len * 0.5;

    PS::MortonKey::initialize(half_len, center);
    PS::U64 n_limit = (((PS::U64)1)<<((PS::U64)PS::TREE_LEVEL_LIMIT));
    std::cerr<<"n_limit= "<<n_limit<<std::endl;
    PS::F64 dx = full_len / n_limit;
    std::cerr<<"dx= "<<dx<<std::endl;
    std::mt19937_64 gen;
    std::uniform_int_distribution<int> rand_int(0, 7);
    std::cerr<<"dx= "<<dx<<std::endl;
    PS::F64vec pos = center + PS::F64vec(dx);
    PS::KeyT key = PS::MortonKey::getKey(pos);
    std::cerr<<"pos= "<<pos
             <<std::oct<<" key= "<<key
             <<std::endl;
    bool flag_err = false;
    PS::S32 shift_id[PS::TREE_LEVEL_LIMIT];
    const PS::S32 n = 1000000;
    for(PS::S32 loop=0; loop<n; loop++){
        if(loop%(loop/100+1)==0) std::cerr<<"loop= "<<loop<<std::endl;
        pos = center;
        shift_id[0] = 0;
        PS::F64 len_tmp = half_len;
        for(PS::S32 i=1; i<PS::TREE_LEVEL_LIMIT; i++){
            shift_id[i] = rand_int(gen);
            pos += PS::SHIFT_CENTER[shift_id[i]]*len_tmp;
            len_tmp *= 0.5;
        }
        key = PS::MortonKey::getKey(pos);

        //std::cerr<<"pos= "<<pos<<" cellid= "<<std::dec;
        //for(PS::S32 i=0; i<PS::TREE_LEVEL_LIMIT; i++){
        //    std::cerr<<std::dec<<shift_id[i];
        //}
        //std::cerr<<std::oct<<" key= "<<key<<std::endl;
        //std::cerr<<std::dec;

        for(PS::S32 i=0; i<PS::TREE_LEVEL_LIMIT; i++){
            if(shift_id[i] != PS::MortonKey::getCellID(i, key)){
                std::cerr<<"loop= "<<loop
                         <<" shift_id["<<i<<"]= "<<shift_id[i]
                         <<" cellid= "
                         <<PS::MortonKey::getCellID(i, key)
                         <<std::endl;
                flag_err = true;
            }

        }
    }
    if(flag_err) std::cerr<<"FAIL"<<std::endl;
    else std::cerr<<"PASS"<<std::endl;
    
    return 0;
}


 */
