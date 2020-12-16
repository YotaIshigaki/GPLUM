//***************************************************************************************
//  test program for "atsevb_FDPS_intraCoefTable.hpp"
//
//  compile: $ g++ cpp_test_intraCoefTable.cpp -std=c++11 -I/(PS_DIR)/src/
//
//***************************************************************************************

#include <iostream>

#include <particle_simulator.hpp>
#include "../md_ext_multi_table.hpp"

enum class Type : int {
  a1,
  a2,
  a3,
  a4,
  a5,
  
  ENUM_N_TOT
};

enum class Type2 : int {
  b1,
  b2,
  b3,
  b4,
  b5,
  
  ENUM_N_TOT
};

enum class Type3 : int {
  c1,
  c2,
  c3,
  c4,
  c5,
  
  ENUM_N_TOT
};

int main() {
    
    //--- default: the number of Type <= 4, value range of "Type" is 0~99.
    MD_EXT::MultiTable<double, Type> test;
    MD_EXT::MultiTable<double, Type, Type> test2;
    MD_EXT::MultiTable<double, Type, Type2, Type3> test3;
    MD_EXT::MultiTable<double, Type, Type, Type, Type> test4;
    
    //--- reserve table size
    int size =   static_cast<int>(Type::ENUM_N_TOT)
               * static_cast<int>(Type::ENUM_N_TOT);
    test2.reserve(size);
        size =   static_cast<int>(Type::ENUM_N_TOT)
               * static_cast<int>(Type2::ENUM_N_TOT)
               * static_cast<int>(Type3::ENUM_N_TOT);
    test3.reserve(size);
    
    //--- key generation test
    std::cout << std::endl << "key generation test" << std::endl;
    PS::S64 key;
    key = test.test_key_gen_(Type::a3);
    std::cout << "test(a3): key = " << key << std::endl;
    
    key = test2.test_key_gen_(Type::a2, Type::a4);
    std::cout << "test2(a2, a4): key = " << key << std::endl;
    
    key = test3.test_key_gen_(Type::a5, Type2::b3, Type3::c1);
    std::cout << "test3(a5, b3, c1): key = " << key << std::endl;
    
    key = test4.test_key_gen_(Type::a4, Type::a2, Type::a1, Type::a5);
    std::cout << "test4(a4, a2, a1, a5): key = " << key << std::endl;
    
    
    //--- table function test
    std::cout << std::endl << "table function test" << std::endl;
    double coef  =   1.0;
    double coef2 =  10.0;
    double coef3 = 100.0;
    test2.set(coef,  Type::a2, Type::a2);
    test2.set(coef2, Type::a1, Type::a5);
    test2.set(coef3, Type::a4, Type::a5);
    
    double result;
    result = test2.get(Type::a2, Type::a2);
    std::cout << "Result(a2,a2) = " << result << std::endl;
    result = test2.get(Type::a1, Type::a5);
    std::cout << "Result(a1,a5) = " << result << std::endl;
    result = test2.get(Type::a4, Type::a5);
    std::cout << "Result(a4,a5) = " << result << std::endl;
    
    
    
    std::cout << std::endl << "error case samples are executed in below." << std::endl << std::endl;
    //--- error case sample, undefined value
    // result = test2.get(Type::a1, Type::a1);
    // std::cout << "Result(a1,a1) = " << result << std::endl;
    
    //--- error case sample, too much key types.
    // MD_EXT::MultiTable<double, Type, Type, Type, Type, Type> too_much_types;
    
    //--- error case sample, invalid Key type.
    //  key = test2.test_key_gen_(Type::a1, Type2::b1);  // the keys of test2 is defined as (Type, Type).
    
    //--- error case sample, invalid type number
    // key = test2.test_key_gen_(Type::a1);                      // invalid type number
    // key = test2.test_key_gen_(Type::a1, Type::a1, Type::a1);  // invalid type number
    
    //--- error case sample, key value is out of range. (compile with -DMULTI_TABLE_WINDOW_TEST)
    // key = test2.test_key_gen_(Type::a1, static_cast<Type>(-99));
    // std::cout << "test(a1, 100): key = " << key << std::endl;
    return 0;
}


