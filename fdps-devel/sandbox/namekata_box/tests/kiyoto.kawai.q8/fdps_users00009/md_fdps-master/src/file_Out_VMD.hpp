//***************************************************************************************
//  This is record file I/O routine.
//    This code is used by "md_fdps_main.cpp"
//***************************************************************************************
#pragma once

#include <cstdio>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "file_IO_defs.hpp"


namespace FILE_IO {

    class VMDFileManager {
    private:
        _Impl::FileManager_base mngr;

    public:
        void init(const std::string &dir,
                  const PS::S64      start     = 0,
                  const PS::S64      interval  = 1,
                  const std::string &name_pre  = "vmd",
                  const std::string &name_post = ".pdb"){
            this->mngr.init(start, interval, dir, name_pre, name_post);
        }
        VMDFileManager() = default;
        VMDFileManager(const std::string &dir,
                       const PS::S64      start     = 0,
                       const PS::S64      interval  = 1,
                       const std::string &name_pre  = "vmd",
                       const std::string &name_post = ".pdb"){
            this->mngr.init(start, interval, dir, name_pre, name_post);
        }

        //--- VMD data output
        template <class T_FP, class Tstat>
        void record(      PS::ParticleSystem<T_FP> &psys,
                    const Tstat                    &stat){

            //--- output cycle
            const auto i_step = stat.get_istep();
            if(  i_step < this->mngr.get_start()   ||
                (i_step % this->mngr.get_interval() ) != 0 ) return;

            const std::string file_name = this->mngr.get_file_name(i_step);
            if(PS::Comm::getRank() == 0) std::cout << "  output " << file_name << std::endl;

            //--- write data
            psys.writeParticleAscii( file_name.c_str(), &T_FP::write_VMD_ascii );
        }

    };


    namespace VMD {

        template <class Tptcl>
        void write_atom(FILE *fp, const Tptcl &atom){
            std::ostringstream oss;
            const PS::F64vec   pos_real = Normalize::realPos( atom.getPos() );

            //--- VMD output
            oss.str("");   // initialize buffer
            oss << "ATOM  "                                                       // 1~6   "ATOM  "
                << std::setw(5) << atom.getAtomID()                               // 7~11  atom id
                << " "                                                            // 12    space
                << std::setw(4) << ENUM::what( atom.getAtomType() )               // 13~16 atom name
                << " "                                                            // 17    alternate location indicator
                << std::setw(3) << atom.getResidue()                              // 18~20 residue name
                << " "                                                            // 21    chain identifier
                << std::setw(4) << atom.getMolID()                                // 22~25 residue sequence number
                << " "                                                            // 26    code for insertion of residues
                << "   "                                                          // 27~29 space
                << " " << std::setw(7) << std::to_string(pos_real.x).substr(0,7)  // 30~37 X position
                << " " << std::setw(7) << std::to_string(pos_real.y).substr(0,7)  // 38~45 Y position
                << " " << std::setw(7) << std::to_string(pos_real.z).substr(0,7)  // 46~53 Z position
                << "\n";

            std::fprintf(fp, oss.str().c_str());
        }

    }

}
