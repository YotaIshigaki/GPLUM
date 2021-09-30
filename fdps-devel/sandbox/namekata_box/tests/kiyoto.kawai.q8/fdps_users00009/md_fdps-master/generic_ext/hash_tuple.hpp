/**************************************************************************************************/
/**
* @file  hash_tuple.hpp
* @brief specialization of std::hash<T> for std::tuple<...>.
*/
/**************************************************************************************************/
#pragma once

#include <tuple>
#include <cstdint>
#include <cassert>


/**
* @brief specialization of std::hash<T> for std::tuple<...>.
* @details ref: https://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set
*/
namespace hash_tuple {

    template<class T>
    struct hash_internal{
        size_t operator () (T const &t) const {
            return std::hash<T>()(t);
        }
    };

    //! @brief combine multiple hash values to one hash value.
    //! @details the "hash_combine()" is same to boost::hash_combine()
    template<class T>
    void hash_combine(size_t &seed, T const &v){
        //--- phi = (1 + sqrt(5)) / 2,   2^32 / phi = 0x9e3779b9
        seed ^= hash_internal<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

    template<class Tuple, size_t Index = std::tuple_size<Tuple>::value-1>
    struct HashValueImpl{
        static void apply(size_t &seed, Tuple const &tuple){
            HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
            hash_combine(seed, std::get<Index>(tuple));
        }
    };

    template<class Tuple>
    struct HashValueImpl<Tuple, 0>{
        static void apply(size_t &seed, Tuple const &tuple){
            hash_combine(seed, std::get<0>(tuple));
        }
    };

    //! @brief specialize interface for std::hash<std::tuple<...>>.
    template<typename Tuple>
    struct hash_func {
        size_t operator () (Tuple const &tt) const {
            size_t seed = 0;
            HashValueImpl<Tuple>::apply(seed, tt);
            return seed;
        }
    };
}
