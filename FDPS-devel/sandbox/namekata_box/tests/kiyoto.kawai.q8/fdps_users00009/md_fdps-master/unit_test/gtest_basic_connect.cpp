//=======================================================================================
//  This is unit test of MD_EXT::basic_connect.
//     module location: ./generic_ext/md_ext_basic_connect.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include <basic_connect.hpp>


//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST(BasicConnect, edit){
    MD_EXT::basic_connect<int, 4> ct;

    EXPECT_TRUE(ct.empty());
    EXPECT_EQ(ct.max_size(), 4);
    EXPECT_EQ(ct.size()    , 0);

    for(int i=0; i<4; ++i){
        ct.add(i);
    }
    EXPECT_EQ(ct.size(), 4);

    EXPECT_EQ(ct.find(-1), ct.end());
    EXPECT_NE(ct.find(0) , ct.end());
    EXPECT_NE(ct.find(1) , ct.end());
    EXPECT_NE(ct.find(2) , ct.end());
    EXPECT_NE(ct.find(3) , ct.end());
    EXPECT_EQ(ct.find(4) , ct.end());

    ct.remove(2);
    EXPECT_EQ(ct.size()  , 3);
    EXPECT_EQ(ct.find(2) , ct.end());
}

TEST(BasicConnect, combine){
    MD_EXT::basic_connect<int, 4> ct1;
    MD_EXT::basic_connect<int, 4> ct2;
    MD_EXT::basic_connect<int, 8> ct3;

    ct1.add(1);
    ct1.add(2);
    ct1.add(3);

    ct2.add(3);
    ct2.add(4);
    ct2.add(5);

    ct3 += ct1;
    ct3 += ct2;
    EXPECT_EQ(ct3.size(), 5);
    EXPECT_NE(ct3.find(1) , ct3.end());
    EXPECT_NE(ct3.find(2) , ct3.end());
    EXPECT_NE(ct3.find(3) , ct3.end());
    EXPECT_NE(ct3.find(4) , ct3.end());
    EXPECT_NE(ct3.find(5) , ct3.end());

    EXPECT_THROW(ct1 += ct2, std::length_error);

    MD_EXT::basic_connect<int, 6> ct_s1;
    EXPECT_NO_THROW(ct_s1 = ct3);

    MD_EXT::basic_connect<int, 2> ct_s2;
    EXPECT_THROW(   ct_s2 = ct3, std::length_error);

}

#include "gtest_main.hpp"
