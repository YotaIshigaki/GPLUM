/**************************************************************************************************/
/**
* @file  str_tool.hpp
* @brief tools for reading parameter files.
*/
/**************************************************************************************************/
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>


/**
* @brief tools for reading parameter files.
*/
namespace STR_TOOL {

    /**
    * @brief remove "CR" code from back of std::string.
    * @param[in, out] str std::string value after devided by "LF" code.
    * @details countermeasure for the difference of "\n" in Linux (LF) and Windows (CR+LF).
    */
    inline void removeCR(std::string &str){
        if(!str.empty() && str.back() == static_cast<char>(13)){
            str.pop_back();
        }
    }

    /**
    * @brief split function for std::string.
    * @param[in] str std::string value.
    * @param[in] delim delimiter.
    * @return std::vector<std::string> result.
    */
    inline std::vector<std::string> split(const std::string &str,
                                          const std::string &delim){

        assert(delim.size() == 1);  // use const std::string instead of const char
        std::istringstream ss{str};

        std::string item;
        std::vector<std::string> result;
        while (std::getline(ss, item, delim.at(0))) {
            if(!item.empty()){
                result.push_back(item);
            }
        }
        return result;
    }
    inline std::vector<std::string> split(const std::string &str,
                                          const char        *delim){
        return split(str, std::string{delim});
    }
    inline std::vector<std::string> split(const std::string &str,
                                          const char         delim){
        return split(str, std::string{delim});
    }

    /**
    * @brief checking str is integer of not.
    * @param[in] str std::string value.
    * @return bool result. True: str is integer. False: str is not integer.
    */
    inline bool isInteger(const std::string &str){
        if(str.find_first_not_of("-+0123456789 \t") != std::string::npos) {
            return false;
        } else {
            return true;
        }
    }

    /**
    * @brief checking str is float of not.
    * @param[in] str std::string value.
    * @return bool result. True: str is float. False: str is not float.
    */
    inline bool isNumeric(const std::string &str){
        if(str.find_first_not_of("-+0123456789.Ee \t") != std::string::npos) {
            return false;
        } else {
            return true;
        }
    }
}
