//----------------------------------------------------------------------------------------
//  This file is external enum functions for user-defined enum classes.
//----------------------------------------------------------------------------------------
#pragma once

#include <string>
#include <tuple>
#include <stdexcept>

#include "enum_model.hpp"
#include "enum_coef.hpp"


namespace ENUM {

    //--- convert std::tuple<...> to std::string
    template<typename T>
    std::string whatis_string_Impl(T const &v){
        return std::to_string(v);
    }

    template<typename Tuple, size_t Index = std::tuple_size<Tuple>::value-1>
    struct whatis_Impl{
        static void apply(std::string &str, Tuple const &tuple){
            whatis_Impl<Tuple, Index-1>::apply(str, tuple);
            str += ", " + whatis_string_Impl(std::get<Index>(tuple));
        }
    };

    template<typename Tuple>
    struct whatis_Impl<Tuple, 0>{
        static void apply(std::string &str, Tuple const &tuple){
            str = whatis_string_Impl(std::get<0>(tuple));
        }
    };

    template<typename Tuple>
    inline std::string what(Tuple const &tt){
        std::string str{};
        whatis_Impl<Tuple>::apply(str, tt);
        return "(" + str + ")";
    }

}
