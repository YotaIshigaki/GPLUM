#define KEY_3D
#define USE_96BIT_KEY

#include<iostream>
#include<iomanip>
#include<random>
#include<ps_defs.hpp>
#include<key.hpp>

namespace PS = ParticleSimulator;

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
        /*
        std::cerr<<"pos= "<<pos<<" cellid= "<<std::dec;
        for(PS::S32 i=0; i<PS::TREE_LEVEL_LIMIT; i++){
            std::cerr<<std::dec<<shift_id[i];
        }
        std::cerr<<std::oct<<" key= "<<key<<std::endl;
        std::cerr<<std::dec;
        */
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
