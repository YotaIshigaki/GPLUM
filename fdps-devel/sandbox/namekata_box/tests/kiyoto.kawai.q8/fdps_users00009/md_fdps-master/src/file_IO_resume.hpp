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

    //--- header class for resume file
    template <class Tstate>
    class ResumeHeader {
      public:
        //--- flag
        PS::S32     i_proc = -1;

        //--- system property
        PS::S32     n_atom;  // must be PS::S32. if n_atom becomes large, use distributed file format.
        PS::S64     i_step;
        PS::F64     time;
        PS::F64     time_trj;
        PS::F64vec  box;

        //--- ext_sys_controller property
        Tstate ext_sys_state;

        //--- definition for tags
        const std::string tag_n_atom           = "n_atom:";
        const std::string tag_i_step           = "i_step:";
        const std::string tag_time             = "time:";
        const std::string tag_time_trj         = "time_trj:";
        const std::string tag_box              = "box:";
        const std::string tag_ext_sys_state    = "ext_sys_state:";
        const std::string tag_ext_sys_particle = "ext_sys_particle:";

        const std::string header_mark = "atom_data:";

        void writeAscii(FILE *fp) const {
            std::ostringstream oss;
            oss <<std::hexfloat;
            oss << this->tag_n_atom   << "\t" << this->n_atom   << "\n"
                << this->tag_i_step   << "\t" << this->i_step   << "\n"
                << this->tag_time     << "\t" << this->time     << "\t" << "[normalized]" << "\n"
                << this->tag_time_trj << "\t" << this->time_trj << "\t" << "[normalized]" << "\n"
                << this->tag_box      << "\t" << this->box.x
                <<                       "\t" << this->box.y
                <<                       "\t" << this->box.z << "\n"
                << "\n";
            oss << this->tag_ext_sys_state
                << "\t" << "n_chain="  << "\t" << this->ext_sys_state.n_chain
                << "\t" << "n_rep="    << "\t" << this->ext_sys_state.n_rep
                << "\t" << "n_nys="    << "\t" << this->ext_sys_state.n_nys
                << "\t" << "NVT_freq=" << "\t" << this->ext_sys_state.NVT_freq
                << "\t" << "NPT_freq=" << "\t" << this->ext_sys_state.NPT_freq
                << "\t" << "v_press="  << "\t" << this->ext_sys_state.v_press << "\n";
            oss << this->tag_ext_sys_particle;
            oss << "\t" << "w_coef=";
            for(const auto e : this->ext_sys_state.w_coef){
                oss << "\t" << e;
            }
            oss << "\t" << "x_nhc=";
            for(const auto e : this->ext_sys_state.x_nhc){
                oss << "\t" << e;
            }
            oss << "\t" << "v_nhc=";
            for(const auto e : this->ext_sys_state.v_nhc){
                oss << "\t" << e;
            }
            oss << "\n";
            oss << "\n";
            oss << this->header_mark << "\n"   // mark for end of header data
                << "AtomID"   << "\t"
                << "MolID"    << "\t"
                << "AtomType" << "\t"
                << "MolType"  << "\t"
                << "Pos_x"    << "\t"
                << "Pos_y"    << "\t"
                << "Pos_z"    << "\t"
                << "Mass"     << "\t"
                << "Vel_x"    << "\t"
                << "Vel_y"    << "\t"
                << "Vel_z"    << "\t"
                << "Trj_x"    << "\t"
                << "Trj_y"    << "\t"
                << "Trj_z"    << "\t"
                << "Charge"   << "\t"
                << "VDW_D"    << "\t"
                << "VDW_R"    << "\t"
                << "bond"     << "\n";
            std::fprintf(fp, oss.str().c_str());
        }
        PS::S32 readAscii(FILE *fp){
            std::map<std::string, std::vector<std::string>> info_map;

            TOOL::load_Ascii_header_info(this->header_mark, "\t", 20, fp, info_map);
            TOOL::line_to_str_list(fp, "\t"); // skip element discription line

            if(info_map[this->tag_n_atom].size()   < 2 ||
               info_map[this->tag_i_step].size()   < 2 ||
               info_map[this->tag_time].size()     < 2 ||
               info_map[this->tag_time_trj].size() < 2 ||
               info_map[this->tag_box].size()      < 4 ||
               info_map[this->tag_ext_sys_state].size() < 13){
                //--- invalid header info
                std::ostringstream oss;
                oss << "invalid header information: format error." << "\n";
                const std::vector<std::string> key_list{ this->tag_n_atom,
                                                         this->tag_i_step,
                                                         this->tag_time,
                                                         this->tag_time_trj,
                                                         this->tag_box,
                                                         this->tag_ext_sys_state, };
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

            this->ext_sys_state.n_chain  = std::stoi( info_map[this->tag_ext_sys_state].at(2)  );
            this->ext_sys_state.n_rep    = std::stoi( info_map[this->tag_ext_sys_state].at(4)  );
            this->ext_sys_state.n_nys    = std::stoi( info_map[this->tag_ext_sys_state].at(6)  );
            this->ext_sys_state.NVT_freq = std::stod( info_map[this->tag_ext_sys_state].at(8)  );
            this->ext_sys_state.NPT_freq = std::stod( info_map[this->tag_ext_sys_state].at(10) );
            this->ext_sys_state.v_press  = std::stod( info_map[this->tag_ext_sys_state].at(12) );

            const PS::S32 n_chain   = this->ext_sys_state.n_chain;
            const PS::S32 n_nys     = this->ext_sys_state.n_nys;
            const PS::S32 st_w_coef = 2;
            const PS::S32 st_x_nhc  = st_w_coef + n_nys   + 1;
            const PS::S32 st_v_nhc  = st_x_nhc  + n_chain + 1;
            if(info_map[this->tag_ext_sys_particle].size() < static_cast<size_t>(st_v_nhc + n_chain) ){
                std::ostringstream oss;
                oss << "invalid header information: data field is not enough." << "\n"
                    << "    len for 'w_coef': " << n_nys   << "\n"
                    << "    len for '*_nhc' : " << n_chain << "\n";
                oss << "    data of '" << this->tag_ext_sys_particle << "' =";
                for(const auto& s : info_map[this->tag_ext_sys_particle]){
                    oss << " " << s;
                }
                oss << "\n";
                throw std::invalid_argument(oss.str());
            }
            this->ext_sys_state.w_coef.clear();
            this->ext_sys_state.x_nhc.clear();
            this->ext_sys_state.v_nhc.clear();
            for(PS::S32 i=0; i<n_nys; ++i){
                this->ext_sys_state.w_coef.push_back( std::stod( info_map[this->tag_ext_sys_particle].at(st_w_coef + i) ) );
            }
            for(PS::S32 i=0; i<n_chain; ++i){
                this->ext_sys_state.x_nhc.push_back( std::stod( info_map[this->tag_ext_sys_particle].at(st_x_nhc + i) ) );
            }
            for(PS::S32 i=0; i<n_chain; ++i){
                this->ext_sys_state.v_nhc.push_back( std::stod( info_map[this->tag_ext_sys_particle].at(st_v_nhc + i) ) );
            }

            this->i_proc = PS::Comm::getRank();
            return this->n_atom;
        }
    };

    class ResumeFileManager {
    private:
        _Impl::FileManager_base mngr;

    public:
        void init(const std::string &dir,
                  const PS::S64      start     = 0,
                  const PS::S64      interval  = 1,
                  const std::string &name_pre  = "resume",
                  const std::string &name_post = ".tsv"){
            this->mngr.init(start, interval, dir, name_pre, name_post);
        }
        ResumeFileManager() = default;
        ResumeFileManager(const std::string &dir,
                          const PS::S64      start     = 0,
                          const PS::S64      interval  = 1,
                          const std::string &name_pre  = "resume",
                          const std::string &name_post = ".tsv"){
            this->mngr.init(start, interval, dir, name_pre, name_post);
        }

        //--- resume data output
        template <class T_FP, class Tstat, class Tcontroller>
        void record(      PS::ParticleSystem<T_FP> &psys,
                    const Tstat                    &stat,
                    const Tcontroller              &controller){

            //--- disable flag
            if(stat.get_resume_interval() <= 0) return;

            //--- output cycle
            const PS::S64 i_step = stat.get_istep();
            if(  i_step < this->mngr.get_start()   ||
                (i_step % this->mngr.get_interval() ) != 0 ) return;

            const std::string file_name = this->mngr.get_file_name(i_step);
            if(PS::Comm::getRank() == 0) std::cout << "  output " << file_name << std::endl;

            //--- create file name
            ResumeHeader<typename Tcontroller::State> header;
            header.n_atom   = psys.getNumberOfParticleGlobal();
            header.i_step   = stat.get_istep();
            header.time     = stat.get_time_raw();
            header.time_trj = stat.get_trj_time_raw();
            header.box      = Normalize::getBoxSize();
            header.ext_sys_state = controller.get_resume();

            //--- write data
            psys.writeParticleAscii( file_name.c_str(), header, &T_FP::write_resume_ascii );
        }

        template <class T_FP, class Tstat, class Tcontroller>
        void load(PS::ParticleSystem<T_FP> &psys,
                  Tstat                    &stat,
                  Tcontroller              &controller){

            const PS::S64 i_step = stat.get_istep();

            const std::string file_name = this->mngr.get_file_name(i_step);
            if(PS::Comm::getRank() == 0) std::cout << "  load " << file_name << std::endl;

            ResumeHeader<typename Tcontroller::State> header;
            psys.readParticleAscii( file_name.c_str(), header, &T_FP::read_resume_ascii );

            //--- broadcast header information
            PS::S32 data_proc = PS::Comm::getMaxValue(header.i_proc);
            COMM_TOOL::broadcast(header.i_step       , data_proc);
            COMM_TOOL::broadcast(header.time         , data_proc);
            COMM_TOOL::broadcast(header.time_trj     , data_proc);
            COMM_TOOL::broadcast(header.box          , data_proc);
            COMM_TOOL::broadcast(header.ext_sys_state, data_proc);

            if(header.i_step != i_step){
                std::ostringstream oss;
                oss.str("");
                oss << "invalid header information." << "\n"
                    << "   stat.get_istep() = " << i_step << "\n"
                    << "   header.i_step    = " << header.i_step << "\n";
                throw std::logic_error(oss.str());
            }

            stat.set_time_raw(    header.time);
            stat.set_trj_time_raw(header.time_trj);
            Normalize::setBoxSize(header.box);
            Normalize::broadcast_boxSize();
            controller.load_resume(header.ext_sys_state);
        }

    };

    namespace Resume {

        template <class Tptcl>
        void write_atom(FILE *fp, const Tptcl &atom){

            //--- FP property output
            std::ostringstream oss;
            oss.str("");
            oss << std::hexfloat;
            oss << atom.getAtomID()                 << "\t"
                << atom.getMolID()                  << "\t"
                << ENUM::what( atom.getAtomType() ) << "\t"
                << ENUM::what( atom.getMolType()  ) << "\t"
                << atom.getPos().x  << "\t"
                << atom.getPos().y  << "\t"
                << atom.getPos().z  << "\t"
                << atom.getMass()   << "\t"
                << atom.getVel().x  << "\t"
                << atom.getVel().y  << "\t"
                << atom.getVel().z  << "\t"
                << atom.getTrj().x  << "\t"
                << atom.getTrj().y  << "\t"
                << atom.getTrj().z  << "\t"
                << atom.getCharge() << "\t"
                << atom.getVDW_D()  << "\t"
                << atom.getVDW_R()  << "\t";

            //--- connection information
            oss << atom.bond.size();
            for(const auto b : atom.bond){
                oss << "\t" << b;
            }
            oss << "\n";

            std::fprintf(fp, oss.str().c_str());
        }

        template <class Tptcl>
        void read_atom(FILE *fp, Tptcl &atom){
            const size_t field_len = 18;
            const auto str_list = TOOL::line_to_str_list(fp, "\t");
            if(str_list.size() < field_len){
                std::ostringstream oss;
                oss << "invalid data field." << "\n"
                    << "  str_list.size() = " << str_list.size() << ", must be >= 15" << "\n";
                throw std::invalid_argument("invalid data field.");
            }

            //--- FP property input
            atom.setAtomID( std::stoi(str_list[0]) );
            atom.setMolID(  std::stoi(str_list[1]) );
            atom.setAtomType( str_list[2] );
            atom.setMolType(  str_list[3] );
            atom.setPos( PS::F32vec{ std::stof(str_list[4]),
                                     std::stof(str_list[5]),
                                     std::stof(str_list[6]) } );
            atom.setMass( std::stof(str_list[7]) );
            atom.setVel( PS::F32vec{ std::stof(str_list[8]),
                                     std::stof(str_list[9]),
                                     std::stof(str_list[10]) } );
            atom.clearTrj();
            atom.addTrj( PS::F32vec{ std::stof(str_list[11]),
                                     std::stof(str_list[12]),
                                     std::stof(str_list[13]) } );
            atom.setCharge( std::stof(str_list[14]) );
            atom.setVDW_D(  std::stof(str_list[15]) );
            atom.setVDW_R(  std::stof(str_list[16]) );

            //--- connection information
            const PS::S32 n_bond = std::stoi(str_list[17]);
            if( n_bond          < 0 ||
                str_list.size() < field_len + static_cast<size_t>(n_bond) ){
                std::ostringstream oss;
                oss << "invalid data of bond information." << "\n"
                    << "  n_bond          = " << n_bond
                    << ", must be >= 0." << "\n"
                    << "  str_list.size() = " << str_list.size()
                    << ", must be >= "   << field_len + n_bond << "\n";
                throw std::invalid_argument(oss.str());
            }
            for(PS::S32 i=0; i<n_bond; ++i){
                atom.bond.add( std::stoi(str_list[field_len+i]) );
            }
        }

    }

}
