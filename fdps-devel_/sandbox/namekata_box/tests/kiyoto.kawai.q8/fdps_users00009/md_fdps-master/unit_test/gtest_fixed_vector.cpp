//=======================================================================================
//  This is unit test of MD_EXT::fixed_vector.
//     module location: ./generic_ext/fixed_vector.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include "fixed_vector.hpp"


//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST(FixedVector, init){
    MD_EXT::fixed_vector<int, 8> vec1;
    MD_EXT::fixed_vector<int, 8> vec2{1,2,3,4};

    EXPECT_TRUE(vec1.empty());
    EXPECT_EQ(vec1.max_size(), 8);
    EXPECT_EQ(vec1.size()    , 0);

    EXPECT_FALSE(vec2.empty());
    EXPECT_EQ(vec2.max_size(), 8);
    EXPECT_EQ(vec2.size()    , 4);

    vec1 = {1,2,3,4,5};
    EXPECT_FALSE(vec1.empty());
    EXPECT_EQ(vec1.max_size(), 8);
    EXPECT_EQ(vec1.size()    , 5);

    vec2.clear();
    EXPECT_TRUE(vec2.empty());
    EXPECT_EQ(vec2.max_size(), 8);
    EXPECT_EQ(vec2.size()    , 0);

    MD_EXT::fixed_vector<int, 6> vec3{vec1};
    EXPECT_NE(vec3.max_size(), vec1.max_size());
    EXPECT_EQ(vec3.size(), vec1.size());
    for(size_t i=0; i<vec3.size(); ++i){
        EXPECT_EQ(vec3[i], vec1[i]);
    }

    auto init_list = {0,1,2,3,4,5,6,7,8,9,10};
    EXPECT_THROW(vec1 = init_list, std::length_error);

    vec1 = {1,2,3,4};
    MD_EXT::fixed_vector<int, 6> vec_s1;
    MD_EXT::fixed_vector<int, 2> vec_s2;
    EXPECT_NO_THROW(vec_s1 = vec1);
    EXPECT_THROW(   vec_s2 = vec1, std::length_error);
}

TEST(FixedVector, access){
    MD_EXT::fixed_vector<int, 8> vec{0,1,2,3,4,5,6};

    int i = 0;
    for(const auto& e : vec){
        EXPECT_EQ(e, i);
        ++i;
    }

    MD_EXT::fixed_vector<int, 8> vec2 = vec;
    EXPECT_EQ(vec2.size(), vec.size());
    for(size_t j=0; j<vec2.size(); ++j){
        EXPECT_EQ(vec2[j], vec.at(j));
    }

    EXPECT_THROW(vec2.at(7), std::out_of_range);
    EXPECT_NO_THROW(vec2[7]);

    for(size_t i=vec2.size(); i<vec2.max_size(); ++i){
        vec2.push_back(i);
    }
    EXPECT_THROW(vec2.push_back(8), std::length_error);
}

TEST(FixedVector, edit){
    MD_EXT::fixed_vector<int, 8> vec{0,1,2};

    vec.push_back(3);
    EXPECT_EQ(vec.size(), 4);
    EXPECT_EQ(vec.at(0) , 0);
    EXPECT_EQ(vec.at(1) , 1);
    EXPECT_EQ(vec.at(2) , 2);
    EXPECT_EQ(vec.at(3) , 3);

    MD_EXT::fixed_vector<int, 8> vec2;
    for(int i=0; i<4; ++i){
        vec2.push_back(i);
    }

    EXPECT_EQ(vec2.size(), vec.size());
    for(size_t i=0; i<vec2.size(); ++i){
        EXPECT_EQ(vec2.at(i), vec.at(i));
    }

    vec.pop_back();
    EXPECT_EQ(vec.size(), 3);
    EXPECT_THROW(vec.at(3), std::out_of_range);

    MD_EXT::fixed_vector<int, 4> vec_s;
    vec_s = vec;
    EXPECT_TRUE(vec_s == vec);

    vec_s.clear();
    EXPECT_TRUE(vec_s.empty());
    EXPECT_NO_THROW(vec_s.pop_back());
}

TEST(FixedVector, resize){
    MD_EXT::fixed_vector<int, 6> vec1;

    vec1.clear();
    vec1.resize(2);
    EXPECT_EQ(vec1.size(), 2);

    vec1.clear();
    vec1.resize(3, -1);
    EXPECT_EQ(vec1.size(), 3);
    for(const auto e : vec1){
        EXPECT_EQ(e, -1);
    }

    vec1.resize(6, 9);
    for(size_t i=0; i<vec1.size(); ++i){
        if(i < 3){
            EXPECT_EQ(vec1.at(i), -1);
        } else {
            EXPECT_EQ(vec1.at(i), 9);
        }
    }

    EXPECT_THROW(vec1.resize(32), std::length_error);
}

TEST(FixedVector, compare){
    MD_EXT::fixed_vector<int, 10> vec1 = {1,1,1,2,1};
    MD_EXT::fixed_vector<int, 8>  vec2 = {1,2,1,1,1};

    EXPECT_TRUE(vec1 < vec2);
    EXPECT_TRUE(vec2 > vec1);

    MD_EXT::fixed_vector<int, 8> vec_cp{vec1};
    EXPECT_TRUE(vec_cp <= vec1);
    EXPECT_TRUE(vec1   <= vec2);
    EXPECT_TRUE(vec_cp >= vec1);
    EXPECT_TRUE(vec2   >= vec1);
    EXPECT_TRUE(vec_cp == vec1);

    vec1.swap(vec2);
    EXPECT_TRUE(vec1 > vec2);
}

#include "gtest_main.hpp"
