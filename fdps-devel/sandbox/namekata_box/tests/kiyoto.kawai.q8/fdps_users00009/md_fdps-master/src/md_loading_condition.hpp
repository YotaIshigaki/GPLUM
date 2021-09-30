//***************************************************************************************
//  This program is loading model parameter function.
//***************************************************************************************
#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <cassert>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "atom_class.hpp"
#include "md_coef_table.hpp"
#include "md_setting.hpp"
#include "ext_sys_control.hpp"


//--- file loading mode
enum class CONDITION_LOAD_MODE : int {
    time_step,
    tree,
    cut_off,
    ext_sys,
    ext_sys_sequence,
    record,

    molecule,
    box,
    ex_radius,
};

//--- std::string converter for enum
namespace ENUM {

    static const std::map<std::string, CONDITION_LOAD_MODE> table_str_CONDITION_LOAD_MODE{
        {"TIMESTEP"        , CONDITION_LOAD_MODE::time_step        },
        {"TREE"            , CONDITION_LOAD_MODE::tree             },
        {"CUT_OFF"         , CONDITION_LOAD_MODE::cut_off          },
        {"EXT_SYS"         , CONDITION_LOAD_MODE::ext_sys          },
        {"EXT_SYS_SEQUENCE", CONDITION_LOAD_MODE::ext_sys_sequence },
        {"RECORD"          , CONDITION_LOAD_MODE::record           },

        {"MOLECULE"        , CONDITION_LOAD_MODE::molecule         },
        {"BOX"             , CONDITION_LOAD_MODE::box              },
        {"EX_RADIUS"       , CONDITION_LOAD_MODE::ex_radius        },
    };

    static const std::map<CONDITION_LOAD_MODE, std::string> table_CONDITION_LOAD_MODE_str{
        {CONDITION_LOAD_MODE::time_step       , "TIMESTEP"         },
        {CONDITION_LOAD_MODE::tree            , "TREE"             },
        {CONDITION_LOAD_MODE::cut_off         , "CUT_OFF"          },
        {CONDITION_LOAD_MODE::ext_sys         , "EXT_SYS"          },
        {CONDITION_LOAD_MODE::ext_sys_sequence, "EXT_SYS_SEQUENCE" },
        {CONDITION_LOAD_MODE::record          , "RECORD"           },

        {CONDITION_LOAD_MODE::molecule        , "MOLECULE"         },
        {CONDITION_LOAD_MODE::box             , "BOX"              },
        {CONDITION_LOAD_MODE::ex_radius       , "EX_RADIUS"        },
    };

    CONDITION_LOAD_MODE which_CONDITION_LOAD_MODE(const std::string &str){
        if(table_str_CONDITION_LOAD_MODE.find(str) != table_str_CONDITION_LOAD_MODE.end()){
            return table_str_CONDITION_LOAD_MODE.at(str);
        } else {
            std::cerr << "  CONDITION_LOAD_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in CONDITION_LOAD_MODE.");
        }
    }

    std::string what(const CONDITION_LOAD_MODE &e){
        if(table_CONDITION_LOAD_MODE_str.find(e) != table_CONDITION_LOAD_MODE_str.end()){
            return table_CONDITION_LOAD_MODE_str.at(e);
        } else {
            using type_base = typename std::underlying_type<CONDITION_LOAD_MODE>::type;
            std::cerr << "  CONDITION_LOAD_MODE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in CONDITION_LOAD_MODE.");
        }
    }

}

inline std::ostream& operator << (std::ostream& s, const CONDITION_LOAD_MODE &e){
    s << ENUM::what(e);
    return s;
}


namespace System {

    //--- load settings from file
    void loading_sequence_condition(const std::string         &file_name,
                                          EXT_SYS::Sequence   &sequence,
                                          EXT_SYS::Controller &controller){
        std::ifstream file_sys{file_name};
        std::string line;

        if(file_sys.fail()) throw std::ios_base::failure("file: " + file_name + " was not found.");

        //--- temporary value
        PS::S64 n_chain, n_rep, n_nys;
        PS::F64 NVT_freq, NPT_freq;

        CONDITION_LOAD_MODE mode = CONDITION_LOAD_MODE::time_step;
        while ( getline(file_sys, line) ) {
            STR_TOOL::removeCR(line);
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            std::string mark     = "@<CONDITION>";
            size_t      mark_len = mark.size();
            if( str_list[0].substr(0,mark_len) == mark ){
                mode = ENUM::which_CONDITION_LOAD_MODE( str_list[0].substr(mark_len) );
                continue;
            }

            //--- loading data
            switch (mode) {
                case CONDITION_LOAD_MODE::time_step:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "i_start"){
                        System::profile.nstep_st = std::stoi(str_list[1]);
                        System::profile.istep    = System::profile.nstep_st;
                    }
                    if( str_list[0] == "i_end") System::profile.nstep_ed = std::stoi(str_list[1]);
                    if( str_list[0] == "dt")    System::profile.dt       = Unit::to_norm_time( std::stof(str_list[1]) );
                break;

                case CONDITION_LOAD_MODE::tree:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "n_leaf_limit")  System::profile.n_leaf_limit  = std::stoi(str_list[1]);
                    if( str_list[0] == "coef_ema")      System::profile.coef_ema      = std::stof(str_list[1]);
                    if( str_list[0] == "theta")         System::profile.theta         = std::stof(str_list[1]);
                    if( str_list[0] == "n_group_limit") System::profile.n_group_limit = std::stoi(str_list[1]);
                    if( str_list[0] == "cycle_dinfo")   System::profile.cycle_dinfo   = std::stoi(str_list[1]);
                break;

                case CONDITION_LOAD_MODE::cut_off:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "LJ")    System::profile.cut_off_LJ    = std::stof(str_list[1]);
                    if( str_list[0] == "intra") System::profile.cut_off_intra = std::stof(str_list[1]);
                break;

                case CONDITION_LOAD_MODE::ext_sys:
                    if( str_list.size() < 2) continue;

                    //--- load controller settings
                    if( str_list[0] == "n_chain"  ) n_chain  = std::stoi(str_list[1]);
                    if( str_list[0] == "n_rep"    ) n_rep    = std::stoi(str_list[1]);
                    if( str_list[0] == "n_nys"    ) n_nys    = std::stoi(str_list[1]);
                    if( str_list[0] == "NVT_freq" ) NVT_freq = std::stof(str_list[1]);
                    if( str_list[0] == "NPT_freq" ) NPT_freq = std::stof(str_list[1]);

                    //--- load default control setting
                    if( str_list[0] == "default" ){
                        sequence.setDefault( EXT_SYS::Setting( ENUM::which_EXT_SYS_MODE( str_list[1] ),
                                                               0,   // period. no meaning for default setting.
                                                               std::stof(str_list[2]),
                                                               std::stof(str_list[3]) ) );
                    }
                break;

                case CONDITION_LOAD_MODE::ext_sys_sequence:
                    if( str_list.size() < 4) continue;

                    //--- load control sequence
                    sequence.addSetting( EXT_SYS::Setting( ENUM::which_EXT_SYS_MODE( str_list[0] ),
                                                           std::stoi(str_list[1]),
                                                           std::stof(str_list[2]),
                                                           std::stof(str_list[3]) ) );
                break;

                case CONDITION_LOAD_MODE::record:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "pos_interval")    System::profile.pos_interval    = stoi(str_list[1]);
                    if( str_list[0] == "pos_start")       System::profile.pos_start       = stoi(str_list[1]);
                    if( str_list[0] == "resume_interval") System::profile.resume_interval = stoi(str_list[1]);
                    if( str_list[0] == "resume_start")    System::profile.resume_start    = stoi(str_list[1]);
                    if( str_list[0] == "pdb_interval")    System::profile.pdb_interval    = stoi(str_list[1]);
                    if( str_list[0] == "pdb_start")       System::profile.pdb_start       = stoi(str_list[1]);

                    if( str_list[0] == "eng_interval")  System::profile.eng_interval  = stoi(str_list[1]);
                    if( str_list[0] == "eng_start")     System::profile.eng_start     = stoi(str_list[1]);
                    if( str_list[0] == "prop_interval") System::profile.prop_interval = stoi(str_list[1]);
                    if( str_list[0] == "prop_start")    System::profile.prop_start    = stoi(str_list[1]);
                break;

                default:
                    std::cerr << "  file: " << file_name << std::endl;
                    throw std::invalid_argument("undefined loading mode.");
            }
        }

        //--- initialize ext_sys controller
        controller.init(n_chain,
                        n_rep,
                        n_nys,
                        NVT_freq,
                        NPT_freq);
    }

    void loading_molecular_condition(const std::string &file_name){
        std::ifstream file_sys{file_name};
        std::string line;

        if(file_sys.fail()) throw std::ios_base::failure("file: " + file_name + " was not found.");

        CONDITION_LOAD_MODE mode = CONDITION_LOAD_MODE::molecule;
        while ( getline(file_sys, line) ) {
            STR_TOOL::removeCR(line);
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            std::string mark     = "@<CONDITION>";
            size_t      mark_len = mark.size();
            if( str_list[0].substr(0,mark_len) == mark ){
                mode = ENUM::which_CONDITION_LOAD_MODE( str_list[0].substr(mark_len) );
                continue;
            }

            //--- loading data
            switch (mode) {
                case CONDITION_LOAD_MODE::molecule:
                    if( str_list.size() < 2) continue;

                    System::model_list.push_back( std::make_pair( ENUM::which_MolName(str_list[0]),
                                                                  std::stoi(str_list[1])           ) );
                break;

                case CONDITION_LOAD_MODE::box:
                    if( str_list.size() < 3) continue;

                    Normalize::setBoxSize( PS::F32vec{ std::stof(str_list[0]),
                                                       std::stof(str_list[1]),
                                                       std::stof(str_list[2]) } );
                break;

                case CONDITION_LOAD_MODE::ex_radius:
                    if( str_list.size() < 2) continue;

                    System::profile.ex_radius = std::stof(str_list[0]);
                    System::profile.try_limit = std::stoi(str_list[1]);
                break;

                default:
                    std::cerr << "  file: " << file_name << std::endl;
                    throw std::invalid_argument("undefined loading mode: " + ENUM::what(mode));
            }
        }
        //--- reach EOF

    }

}
