#pragma once
/* C++ headers */
#include "common.h"

namespace debug_utilities {

    extern std::uint64_t nstep;
    extern std::string filename;
    extern std::ofstream fout;
    
    void fopen();
    void fclose();

}
namespace dbg_utils = debug_utilities;
