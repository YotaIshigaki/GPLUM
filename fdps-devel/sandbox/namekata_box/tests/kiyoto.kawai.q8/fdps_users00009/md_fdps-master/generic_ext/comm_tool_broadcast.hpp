/**************************************************************************************************/
/**
* @file  comm_tool_broadcast.hpp
* @brief STL container wrapper for PS::Comm::broadcast()
*/
/**************************************************************************************************/
#pragma once

#include <cassert>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>

#include <particle_simulator.hpp>

#include "comm_tool_SerDes.hpp"


//--- implementation in "ps_defs.hpp"
//        ///////////////////////
//        // MPI BCAST WRAPPER //
//        // new functions 10 Feb 2015
//        template<class T>
//        static inline void broadcast(T * val, const int n, const int src=0){
//#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
//            MPI::COMM_WORLD.Bcast(val, n, GetDataType<T>(), src);
//#else
//            // NOP
//#endif
//        }


namespace COMM_TOOL {

    //=====================
    //  wrapper functions
    //=====================

    /**
    * @brief wrapper for static size class.
    * @param[in,out] v target. MUST NOT contain pointer member.
    * @param[in] root ID of source process.
    */
    template <typename T>
    void broadcast(T &v, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** bc_<T>" << std::endl;
        #endif
        PS::Comm::broadcast(&v, 1, root);
    }

    /**
    * @brief specialization for std::vector<T>.
    * @param[in,out] vec target. class T MUST NOT contain pointer member.
    * @param[in] root ID of source process.
    * @details adapting size of the vector between all process and broadcast.
    */
    template <typename T>
    void broadcast(std::vector<T> &vec, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == root )  std::cout << " *** bc_vec<T>  len=" << vec.size() << std::endl;
        #endif
        //--- pass vector length
        size_t size = vec.size();
        broadcast(size ,root);

        if(PS::Comm::getRank() != root){
            vec.clear();
            vec.resize(size);
        }
        if(size == 0) return;

        //--- transport vector
        PS::Comm::broadcast(&vec[0], size, root);
    }

    /**
    * @brief specialization for std::vector<std::string>.
    * @param[in,out] vec_str target.
    * @param[in] root ID of source process.
    * @details transcode std::vector<std::string> to std::vector<char> and call broadcast() for vector<T>.
    */
    void broadcast(std::vector<std::string> &vec_str, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == root )  std::cout << " *** bc_vec<string>  len=" << vec_str.size() << std::endl;
        #endif

        std::vector<char> char_buff;
        serialize_vector_string(vec_str, char_buff);

        broadcast(char_buff, root);

        deserialize_vector_string(char_buff, vec_str);
    }

    /**
    * @brief specialization for std::vector<std::vector<T>>.
    * @param[in,out] vec_vec target.
    * @param[in] root ID of source process.
    * @details devide std::vector<std::vector<T>> into std::vector<T> and std::vector<index>, then call broadcast() for std::vector<T>.
    * @details class T accepts std::vector<>, std::string, or user-defined class WITHOUT pointer member.
    * @details If 3 or more nested std::vector<...> is passed, this function will call itself recurcively.
    */
    template <typename T>
    void broadcast(std::vector<std::vector<T>> &vec_vec, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == root ) std::cout << " *** bc_vec_vec<T>  len=" << vec_vec.size() << std::endl;
        #endif

        std::vector<T>      vec_data;
        std::vector<size_t> vec_index;

        serialize_vector_vector(vec_vec,
                                vec_data, vec_index);

        broadcast(vec_data,  root);
        broadcast(vec_index, root);

        deserialize_vector_vector(vec_data, vec_index,
                                  vec_vec             );
    }

    /**
    * @brief specialization for std::vector<std::pair<Ta, Tb>>.
    * @param[in,out] vec_pair target.
    * @param[in] root ID of source process.
    * @details devide std::vector<std::pair<Ta, Tb>> into std::vector<Ta> and std::vector<Tb>, then call broadcast() for std::vector<T>.
    * @details class Ta and Tb accept std::vector<>, std::string, or user-defined class WITHOUT pointer member.
    */
    template <typename Ta, typename Tb>
    void broadcast(std::vector<std::pair<Ta, Tb>> &vec_pair, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == root )  std::cout << " *** bc_vec<pair<Ta, Tb>  len=" << vec_pair.size() << std::endl;
        #endif
        //--- devide first and second
        std::vector<Ta> vec_1st;
        std::vector<Tb> vec_2nd;

        split_vector_pair(vec_pair,
                          vec_1st, vec_2nd);

        broadcast(vec_1st, root);
        broadcast(vec_2nd, root);

        combine_vector_pair(vec_1st, vec_2nd,
                            vec_pair        );
    }

    /**
    * @brief specialize for STL container with insert().
    * @param[in,out] c target container.
    * @param[in] root ID of source process.
    * @details convert the container into std::vector<Telem>, then call broadcast().
    * @details value of class Tc accepts std::vector<>, std::string, or user-defined class WITHOUT pointer member.
    */
    template<typename Tc, typename Telem>
    void broadcast_insert_container(Tc &c, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == root )  std::cout << " *** bc_insert<T>  len=" << c.size() << std::endl;
        #endif
        //--- buffer in std::vector
        std::vector<Telem> buff;
        buff.reserve(c.size());
        for(auto &e : c){
            buff.push_back(e);
        }

        //--- broadcast buffer
        broadcast(buff, root);

        //--- load from buffer
        c.clear();
        for(auto &elem : buff){
            c.insert(elem);
        }
    }

    //! @brief wrapper for std::unordered_map<>
    template <class Key, class T, class Hash, class Pred, class Allocator>
    void broadcast(std::unordered_map<Key, T, Hash, Pred, Allocator> &map, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == root )  std::cout << " *** bc_unordered_map<T>  len=" << map.size() << std::endl;
        #endif
        broadcast_insert_container<std::unordered_map<Key, T, Hash, Pred, Allocator>,
                                   std::pair<Key, T>
                                   >(map, root);
    }

    //! @brief wrapper for std::unordered_multimap<>
    template <class Key, class T, class Hash, class Pred, class Allocator>
    void broadcast(std::unordered_multimap<Key, T, Hash, Pred, Allocator> &map, const PS::S32 &root){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == root )  std::cout << " *** bc_unordered_multimap<T>  len=" << map.size() << std::endl;
        #endif
        broadcast_insert_container<std::unordered_multimap<Key, T, Hash, Pred, Allocator>,
                                   std::pair<Key, T>
                                   >(map, root);
    }

}
