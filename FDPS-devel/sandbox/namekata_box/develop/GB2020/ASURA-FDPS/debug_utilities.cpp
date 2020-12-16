/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "debug_utilities.h"

namespace debug_utilities {

    std::uint64_t nstep {0};
    std::string filename;
    std::ofstream fout;
    
    void fopen() {
        // Set file name
        std::stringstream ss;
        ss << "debug" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
        filename = ss.str();
        // Check the existence of the file
        std::ifstream ifs(filename);
        const bool is_exist = ifs.is_open();
        // Open the file
        if (is_exist) {
            fout.open(filename.c_str(), std::ios::app);
            auto clock = std::chrono::system_clock::now();
            std::time_t time = std::chrono::system_clock::to_time_t(clock);
            fout << "Alert: resumed at " << std::ctime(&time);
        } else {
            fout.open(filename.c_str(), std::ios::trunc);
        }
    }
    
    void fclose() {
        fout.close();
    }

}
