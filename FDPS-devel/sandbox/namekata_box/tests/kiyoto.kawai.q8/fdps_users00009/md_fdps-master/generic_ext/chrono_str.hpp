/**************************************************************************************************/
/**
* @file  chrono_str.hpp
* @brief displaying wrapper for std::chrono.
*/
/**************************************************************************************************/
#pragma once

#include <string>
#include <chrono>

/**
* @brief wrappers for displaying std::chrono results..
*/
namespace chrono_str {

    template<class T>
    std::string to_str_2digit(const T &n){
        std::string tmp = "";
        if(n < 10){
            tmp = "0";
        } else {
            tmp = "";
        }
        tmp += std::to_string(n);

        return tmp;
    }

    /**
    * @brief get string of "milliseconds".
    * @param[in] t duration of std::chrono.
    * @return std::string result. time as "milliseconds".
    */
    template<class Tchrono_value>
    std::string to_str_ms(const Tchrono_value &t){
        return std::to_string( std::chrono::duration_cast<std::chrono::milliseconds>(t).count() );
    }

    /**
    * @brief get string of "sec.msec".
    * @param[in] t duration of std::chrono.
    * @return std::string result. time as "seconds.milliseconds".
    */
    template<class Tchrono_value>
    std::string to_str_s_ms(const Tchrono_value &t){
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(t);
        auto num = sec.count();
        std::string tmp;

        tmp = to_str_2digit(num) + "." + to_str_ms(t - sec);

        return tmp;
    }

    /**
    * @brief get string of "min:sec.msec".
    * @param[in] t duration of std::chrono.
    * @return std::string result. time as "minutes:seconds.milliseconds".
    */
    template<class Tchrono_value>
    std::string to_str_m_s_ms(const Tchrono_value &t){
        auto min = std::chrono::duration_cast<std::chrono::minutes>(t);
        auto num = min.count();
        std::string tmp;

        tmp = to_str_2digit(num) + ":" + to_str_s_ms(t - min);

        return tmp;
    }

    /**
    * @brief get string of "hours:min:sec.msec".
    * @param[in] t duration of std::chrono.
    * @return std::string result. time as "hours:minutes:seconds.milliseconds".
    */
    template<class Tchrono_value>
    std::string to_str_h_m_s_ms(const Tchrono_value &t){
        auto hour = std::chrono::duration_cast<std::chrono::hours>(t);
        auto num  = hour.count();
        std::string tmp;

        tmp = to_str_2digit(num) + ":" + to_str_m_s_ms(t - hour);

        return tmp;
    }

};
