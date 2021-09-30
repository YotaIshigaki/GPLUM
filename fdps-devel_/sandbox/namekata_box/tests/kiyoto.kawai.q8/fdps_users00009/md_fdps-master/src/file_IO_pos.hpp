//***************************************************************************************
//  This is file I/O for trajecctory data.
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

    //--- header class for pos file
    class PosHeader {
      public:
        //--- flag
        PS::S32     i_proc = -1;

        //--- system property
        PS::S32     n_atom;  // must be PS::S32. if n_atom becoms large, use distributed file format.
        PS::S64     i_step;
        PS::F64     time;
        PS::F64     time_trj;
        PS::F64vec  box;

        //--- definition for tags
        const std::string tag_n_atom   = "n_atom:";
        const std::string tag_i_step   = "i_step:";
        const std::string tag_time     = "time:";
        const std::string tag_time_trj = "time_trj:";
        const std::string tag_box      = "box:";

        const std::string header_mark = "atom_data:";

        void writeAscii(FILE *fp) const {
            std::ostringstream oss;
            oss << std::hexfloat;
            oss << this->tag_n_atom   << "\t" << this->n_atom   << "\n"
                << this->tag_i_step   << "\t" << this->i_step   << "\n"
                << this->tag_time     << "\t" << this->time     << "\t" << "[fs]" << "\n"
                << this->tag_time_trj << "\t" << this->time_trj << "\t" << "[fs]" << "\n"
                << this->tag_box      << "\t" << this->box.x
                <<                       "\t" << this->box.y
                <<                       "\t" << this->box.z << "\n"
                << "\n";
            oss << this->header_mark << "\n";   // mark for end of header data
            oss << "AtomID"   << "\t"
                << "MolID"    << "\t"
                << "AtomType" << "\t"
                << "MolType"  << "\t"
                << "Pos_x"    << "\t"
                << "Pos_y"    << "\t"
                << "Pos_z"    << "\t"
                << "Trj_x"    << "\t"
                << "Trj_y"    << "\t"
                << "Trj_z"    << "\n";
            std::fprintf(fp, oss.str().c_str());
        }
        PS::S32 readAscii(FILE *fp){
            std::map<std::string, std::vector<std::string>> info_map;

            TOOL::load_Ascii_header_info(this->header_mark, "\t", 20, fp, info_map);
            TOOL::line_to_str_list(fp, "\t");  // skip element type discription line

            if(info_map[this->tag_n_atom].size()   < 2 ||
               info_map[this->tag_i_step].size()   < 2 ||
               info_map[this->tag_time].size()     < 2 ||
               info_map[this->tag_time_trj].size() < 2 ||
               info_map[this->tag_box].size()      < 4   ){
                //--- invalid header info
                std::ostringstream oss;
                oss << "invalid header information: format error." << "\n";
                const std::vector<std::string> key_list{ this->tag_n_atom,
                                                         this->tag_i_step,
                                                         this->tag_time,
                                                         this->tag_time_trj,
                                                         this->tag_box,   };
                for(const auto& key : key_list){
                    oss << "  key= " << key << ": string=";
                    for(const auto& s : info_map[key]){
                        oss << " " << s;
                    }
                    oss << "\n";
                }
                throw std::invalid_argument(oss.str());
            }

            this->n_atom   = std::stoi( info_map[this->tag_n_atom].at(1) );
            this->i_step   = std::stoi( info_map[this->tag_i_step].at(1) );
            this->time     = std::stod( info_map[this->tag_time].at(1) );
            this->time_trj = std::stod( info_map[this->tag_time_trj].at(1) );
            this->box      = PS::F64vec{ std::stod( info_map[this->tag_box].at(1) ),
                                         std::stod( info_map[this->tag_box].at(2) ),
                                         std::stod( info_map[this->tag_box].at(3) ) };

            this->i_proc = PS::Comm::getRank();
            return this->n_atom;
        }
    };

    class PosFileManager {
    private:
        _Impl::FileManager_base mngr;

    public:
        void init(const std::string &dir,
                  const PS::S64      start     = 0,
                  const PS::S64      interval  = 1,
                  const std::string &name_pre  = "pos",
                  const std::string &name_post = ".tsv"){
            this->mngr.init(start, interval, dir, name_pre, name_post);
        }
        PosFileManager() = default;
        PosFileManager(const std::string &dir,
                       const PS::S64      start     = 0,
                       const PS::S64      interval  = 1,
                       const std::string &name_pre  = "pos",
                       const std::string &name_post = ".tsv"){
            this->mngr.init(start, interval, dir, name_pre, name_post);
        }

        template <class T_FP, class Tstat>
        void record(PS::ParticleSystem<T_FP> &psys,
                    Tstat                    &stat){

            //--- disable flag
            if(stat.get_pos_interval() <= 0) return;

            //--- output cycle
            const PS::S64 i_step = stat.get_istep();
            if(  i_step < this->mngr.get_start()   ||
                (i_step % this->mngr.get_interval()  ) != 0 ) return;

            const std::string file_name = this->mngr.get_file_name( stat.get_istep() );
            if(PS::Comm::getRank() == 0) std::cout << "  output " << file_name << std::endl;

            //--- create header
            PosHeader header;
            header.n_atom   = psys.getNumberOfParticleGlobal();
            header.i_step   = stat.get_istep();
            header.time     = stat.get_time();
            header.time_trj = stat.get_trj_time();
            header.box      = Normalize::getBoxSize();

            //--- write data
            psys.writeParticleAscii( file_name.c_str(), header, &T_FP::write_pos_ascii );

            //--- clear trajectory data
            const PS::S64 n_local = psys.getNumberOfParticleLocal();
            for(PS::S64 i=0; i<n_local; ++i){
                psys[i].clearTrj();
            }
            stat.clear_trj_time();
        }

        template <class T_FP, class Tstat>
        void load(PS::ParticleSystem<T_FP> &psys,
                  Tstat                    &stat){

            const PS::S64 i_step = stat.get_istep();

            const std::string file_name = this->mngr.get_file_name( stat.get_istep() );
            if(PS::Comm::getRank() == 0) std::cout << "  load " << file_name << std::endl;

            PosHeader header;
            psys.readParticleAscii( file_name.c_str(), header, &T_FP::read_pos_ascii );

            //--- broadcast header information
            PS::S32 data_proc = PS::Comm::getMaxValue(header.i_proc);
            COMM_TOOL::broadcast(header.i_step  , data_proc);
            COMM_TOOL::broadcast(header.time    , data_proc);
            COMM_TOOL::broadcast(header.time_trj, data_proc);
            COMM_TOOL::broadcast(header.box     , data_proc);

            if(header.i_step != i_step){
                std::ostringstream oss;
                oss.str("");
                oss << "invalid header information: value error." << "\n"
                    << "   stat.get_istep() = " << i_step << "\n"
                    << "   header.i_step    = " << header.i_step << "\n";
                throw std::logic_error(oss.str());
            }

            stat.set_time(    header.time);
            stat.set_trj_time(header.time_trj);
            Normalize::setBoxSize(header.box);
            Normalize::broadcast_boxSize();
        }
    };

    namespace Pos {

        //--- atom I/O implementation
        template <class Tptcl>
        void write_atom(FILE *fp, const Tptcl &atom){

            //--- ID & position data output
            std::ostringstream oss;
            oss.str("");   // initialize buffer
            oss << std::hexfloat;
            oss << atom.getAtomID()                 << "\t"
                << atom.getMolID()                  << "\t"
                << ENUM::what( atom.getAtomType() ) << "\t"
                << ENUM::what( atom.getMolType()  ) << "\t"
                << atom.getPos().x << "\t"
                << atom.getPos().y << "\t"
                << atom.getPos().z << "\t"
                << atom.getTrj().x << "\t"
                << atom.getTrj().y << "\t"
                << atom.getTrj().z << "\n";

            std::fprintf(fp, oss.str().c_str());
        }

        template <class Tptcl>
        void read_atom(FILE *fp, Tptcl &atom){
            const auto str_list = TOOL::line_to_str_list(fp, "\t");
            if(str_list.size() < 10){
                throw std::invalid_argument("invalid data field.");
            }

            //--- ID & position data input
            atom.setAtomID( std::stoi(str_list[0]) );
            atom.setMolID(  std::stoi(str_list[1]) );
            atom.setAtomType( str_list[2] );
            atom.setMolType(  str_list[3] );
            atom.setPos( PS::F32vec{ std::stof(str_list[4]),
                                     std::stof(str_list[5]),
                                     std::stof(str_list[6]) } );
            atom.clearTrj();
            atom.addTrj( PS::F32vec{ std::stof(str_list[7]),
                                     std::stof(str_list[8]),
                                     std::stof(str_list[9]) } );

        }

    }

}
