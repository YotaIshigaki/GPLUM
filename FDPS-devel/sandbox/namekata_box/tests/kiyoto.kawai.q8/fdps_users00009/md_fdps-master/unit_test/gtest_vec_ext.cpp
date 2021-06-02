//=======================================================================================
//  This is unit test of PS::Vector3<T> and external function.
//     module location: ./generic_ext/vec_ext.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include "vec_ext.hpp"



//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST(Vector3, constructor){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);

    ASSERT_EQ(v0, PS::F64vec(0.0, 0.0, 0.0) );
    ASSERT_EQ(v1, PS::F64vec(2.0, 2.0, 2.0) );
}

TEST(Vector3, assign){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);

    v0 = v1;
    ASSERT_EQ(v0, PS::F64vec(2.0, 2.0, 2.0) );
}

TEST(Vector3, operatorPlus){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);
    PS::F64vec v2 = PS::F64vec(1.0, 0.0, -3.0);

    v0 = v1 + v2;
    ASSERT_EQ(v0, PS::F64vec(3.0, 2.0, -1.0) );
}

TEST(Vector3, operatorPlusEQ){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);
    PS::F64vec v2 = PS::F64vec(1.0, 0.0, -3.0);

    v0 = v1;
    v0 += v2;
    ASSERT_EQ(v0, PS::F64vec(3.0, 2.0, -1.0) );
}

TEST(Vector3, operatorMinus){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);
    PS::F64vec v2 = PS::F64vec(1.0, 0.0, -3.0);

    v0 = v1 - v2;
    ASSERT_EQ(v0, PS::F64vec(1.0, 2.0, 5.0) );
}

TEST(Vector3, operatorMinusEQ){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);
    PS::F64vec v2 = PS::F64vec(1.0, 0.0, -3.0);

    v0 = v1;
    v0 -= v2;
    ASSERT_EQ(v0, PS::F64vec(1.0, 2.0, 5.0) );
}

TEST(Vector3, MultipleScalar){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);
    PS::F64vec v2 = PS::F64vec(1.0, 0.0, -3.0);

    v0 = v1*(-10.0);
    ASSERT_EQ(v0, PS::F64vec(-20.0, -20.0, -20.0) );
    v0 = v2*(2.0);
    ASSERT_EQ(v0, PS::F64vec(2.0, 0.0, -6.0) );
}

TEST(Vector3, InnerProduct){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);
    PS::F64vec v2 = PS::F64vec(1.0, 0.0, -3.0);

    v0 = v1*v2;
    ASSERT_EQ(v0, PS::F64(-4.0) );
    v0 = VEC_EXT::dot(v1, v2);
    ASSERT_EQ(v0, PS::F64(-4.0) );
}

TEST(Vector3, CrossProduct){
    PS::F64vec v0;
    PS::F64vec v1 = PS::F64vec(2.0);
    PS::F64vec v2 = PS::F64vec(1.0, 0.0, -3.0);

    v0 = v1^v2;
    ASSERT_EQ(v0, PS::F64vec(-6.0, 8.0, -2.0) );
    v0 = VEC_EXT::cross(v1, v2);
    ASSERT_EQ(v0, PS::F64vec(-6.0, 8.0, -2.0) );
}

TEST(Vector3, RotateXAxis){
    const PS::F64 pi         = 3.141592653589793;
    const PS::F64 eps        = 1.e-10;
    const PS::F64 _1_sqrt2   = 1.0/std::sqrt(2.0);
    const PS::F64 _45deg_rad = 1.0/4.0*pi;

    PS::F64vec v0, v_ref;
    PS::F64vec v1 = PS::F64vec(0.0, 1.0, 0.0);
    PS::F64vec v2 = PS::F64vec(0.0, _1_sqrt2, _1_sqrt2);

    v0    = VEC_EXT::rot_x(v1, _45deg_rad);
    v_ref = PS::F64vec(0.0, _1_sqrt2, _1_sqrt2);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_x(v2, _45deg_rad);
    v_ref = PS::F64vec(0.0, 0.0, 1.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_x(v1, 3.0*_45deg_rad);
    v_ref = PS::F64vec(0.0, -_1_sqrt2, _1_sqrt2);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_x(v2, 3.0*_45deg_rad);
    v_ref = PS::F64vec(0.0, -1.0, 0.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
}

TEST(Vector3, RotateYAxis){
    const PS::F64 pi         = 3.141592653589793;
    const PS::F64 eps        = 1.e-10;
    const PS::F64 _1_sqrt2   = 1.0/std::sqrt(2.0);
    const PS::F64 _45deg_rad = 1.0/4.0*pi;

    PS::F64vec v0, v_ref;
    PS::F64vec v1 = PS::F64vec(1.0, 0.0, 0.0);
    PS::F64vec v2 = PS::F64vec(_1_sqrt2, 0.0, _1_sqrt2);

    v0    = VEC_EXT::rot_y(v1, _45deg_rad);
    v_ref = PS::F64vec(_1_sqrt2, 0.0, _1_sqrt2);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_y(v2, _45deg_rad);
    v_ref = PS::F64vec(1.0, 0.0, 0.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_y(v1, 3.0*_45deg_rad);
    v_ref = PS::F64vec(-_1_sqrt2, 0.0, _1_sqrt2);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_y(v2, 3.0*_45deg_rad);
    v_ref = PS::F64vec(0.0, 0.0, -1.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
}

TEST(Vector3, RotateZAxis){
    const PS::F64 pi         = 3.141592653589793;
    const PS::F64 eps        = 1.e-10;
    const PS::F64 _1_sqrt2   = 1.0/std::sqrt(2.0);
    const PS::F64 _45deg_rad = 1.0/4.0*pi;

    PS::F64vec v0, v_ref;
    PS::F64vec v1 = PS::F64vec(1.0, 0.0, 0.0);
    PS::F64vec v2 = PS::F64vec(_1_sqrt2, _1_sqrt2, 0.0);

    v0    = VEC_EXT::rot_z(v1, _45deg_rad);
    v_ref = PS::F64vec(_1_sqrt2, _1_sqrt2, 0.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_z(v2, _45deg_rad);
    v_ref = PS::F64vec(0.0, 1.0, 0.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_z(v1, 3.0*_45deg_rad);
    v_ref = PS::F64vec(-_1_sqrt2, _1_sqrt2, 0.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;

    v0    = VEC_EXT::rot_z(v2, 3.0*_45deg_rad);
    v_ref = PS::F64vec(-1.0, 0.0, 0.0);
    v0 -= v_ref;
    ASSERT_TRUE(v0.x < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.y < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
    ASSERT_TRUE(v0.z < eps) << "v0 = " << v0 << " | v_ref = " << v_ref;
}

#include "gtest_main.hpp"
