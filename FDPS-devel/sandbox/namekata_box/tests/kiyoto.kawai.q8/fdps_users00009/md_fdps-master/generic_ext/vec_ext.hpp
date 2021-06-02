/**************************************************************************************************/
/**
* @file  vec_ext.hpp
* @brief extention tools for PS::Vector3<T>.
*/
/**************************************************************************************************/
#pragma once

#include <cmath>

#include <particle_simulator.hpp>


/**
* @brief extention tools for PS::Vector3<T>.
*/
namespace VEC_EXT {

    //! @brief rotate on axis X.
    template <class T>
    PS::Vector3<T> rot_x(const PS::Vector3<T> &v, const PS::F64 &r) {
        PS::Vector3<T> result;
        T sin_tmp = std::sin(r);
        T cos_tmp = std::cos(r);
        result.x = v.x;
        result.y = cos_tmp*v.y - sin_tmp*v.z;
        result.z = sin_tmp*v.y + cos_tmp*v.z;
        return result;
    }

    //! @brief rotate on axis Y.
    template <class T>
    PS::Vector3<T> rot_y(const PS::Vector3<T> &v, const PS::F64 &r) {
        PS::Vector3<T> result;
        T sin_tmp = std::sin(r);
        T cos_tmp = std::cos(r);
        result.x = sin_tmp*v.z + cos_tmp*v.x;
        result.y = v.y;
        result.z = cos_tmp*v.z - sin_tmp*v.x;
        return result;
    }

    //! @brief rotate on axis Z.
    template <class T>
    PS::Vector3<T> rot_z(const PS::Vector3<T> &v, const PS::F64 &r) {
        PS::Vector3<T> result;
        T sin_tmp = std::sin(r);
        T cos_tmp = std::cos(r);
        result.x = cos_tmp*v.x - sin_tmp*v.y;
        result.y = sin_tmp*v.x + cos_tmp*v.y;
        result.z = v.z;
        return result;
    }

    //! @brief dot product (syntactic sugar of vec*vec).
    template <typename T>
    T dot(const PS::Vector3<T>& v1,
          const PS::Vector3<T>& v2){
        return v1*v2;
    }

    //! @brief cross product (syntactic sugar of vec^vec).
    template <typename T>
    PS::Vector3<T> cross(const PS::Vector3<T>& v1,
                         const PS::Vector3<T>& v2){
        return v1^v2;
    }

    //! @brief vector norm^2.
    template <typename T>
    T sq(const PS::Vector3<T>& v){
        return v*v;
    }

    //! @brief vector norm.
    template <typename T>
    T norm(const PS::Vector3<T>& v){
        return std::sqrt(v*v);
    }
}
