/**************************************************************************************************/
/**
* @file  comm_tool_SerDes.hpp
* @brief STL container serializer/deserializer for COMM_TOOL
*/
/**************************************************************************************************/
#pragma once

#include <cassert>
#include <vector>
#include <string>
#include <utility>

#include <particle_simulator.hpp>


namespace COMM_TOOL {

    //=====================
    //  support functions
    //=====================

    /**
    * @brief serialize nested std::vector<T> into 2 std::vector<>s.
    */
    template <typename Tdata, typename Tindex>
    void serialize_vector_vector(const std::vector<std::vector<Tdata>> &vec_vec,
                                       std::vector<Tdata>              &vec_data,
                                       std::vector<Tindex>             &vec_index){

        vec_data.clear();
        vec_index.clear();

        size_t count = 0;
        for(auto &vec : vec_vec){
            vec_index.emplace_back(count);  //start point
            for(auto &elem : vec){
                vec_data.emplace_back(elem);
                ++count;
            }
        }
        vec_index.push_back(count);  // terminater
    }

    /**
    * @brief deserialize nested std::vector<T> from 2 std::vector<>s.
    */
    template <typename Tdata, typename Tindex>
    void deserialize_vector_vector(const std::vector<Tdata>              &vec_data,
                                   const std::vector<Tindex>             &vec_index,
                                         std::vector<std::vector<Tdata>> &vec_vec){

        vec_vec.clear();
        vec_vec.resize(vec_index.size()-1);

        for(size_t i=0; i<vec_index.size()-1; ++i){
            vec_vec.at(i).clear();
            size_t index_begin = vec_index.at(i);
            size_t index_end   = vec_index.at(i+1);
            for(size_t j=index_begin; j<index_end; ++j){
                vec_vec.at(i).emplace_back(vec_data.at(j));
            }
        }
    }

    /**
    * @brief split nested std::vector<std::pair<Ta, Tb>> into 2 std::vector<>s.
    */
    template <typename Ta, typename Tb>
    void split_vector_pair(const std::vector<std::pair<Ta, Tb>> &vec_pair,
                                 std::vector<Ta>                &vec_1st,
                                 std::vector<Tb>                &vec_2nd  ){

        vec_1st.clear();
        vec_2nd.clear();

        for(auto &p : vec_pair){
            vec_1st.emplace_back(p.first);
            vec_2nd.emplace_back(p.second);
        }
    }

    /**
    * @brief split nested std::vector<std::pair<Ta, Tb>> into 2 std::vector<>s.
    */
    template <typename Ta, typename Tb>
    void combine_vector_pair(const std::vector<Ta>                &vec_1st,
                             const std::vector<Tb>                &vec_2nd,
                                   std::vector<std::pair<Ta, Tb>> &vec_pair){

        assert(vec_1st.size() == vec_2nd.size());
        vec_pair.clear();

        for(size_t i=0; i<vec_1st.size(); ++i){
            vec_pair.emplace_back( std::make_pair( vec_1st.at(i),
                                                   vec_2nd.at(i) ) );
        }
    }

    /**
    * @brief serialize std::vector<std::string> into std::vector<char>.
    */
    void serialize_vector_string(const std::vector<std::string> &vec_str,
                                       std::vector<char>        &vec_char){

        vec_char.clear();
        for(auto &str : vec_str){
            size_t len = str.size();
            for(size_t index=0; index<len; ++index){
                vec_char.emplace_back(str[index]);
            }
            //--- add char terminator
            vec_char.emplace_back('\0');
        }
    }

    /**
    * @brief deserialize std::vector<char> into std::vector<std::string>.
    */
    void deserialize_vector_string(const std::vector<char>        &vec_char,
                                         std::vector<std::string> &vec_str  ){

        vec_str.clear();
        std::string str_elem;
        str_elem.clear();
        for(auto &c : vec_char){
            if( c != '\0' ){
                str_elem.push_back(c);
            } else {
                vec_str.push_back(str_elem);
                str_elem.clear();
            }
        }
    }

}
