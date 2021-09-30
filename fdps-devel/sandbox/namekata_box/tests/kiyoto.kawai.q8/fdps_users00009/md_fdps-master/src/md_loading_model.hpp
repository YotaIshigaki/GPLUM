//***************************************************************************************
//  This is loading model parameter function.
//***************************************************************************************
#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include <molecular_dynamics_ext.hpp>

#include "md_defs.hpp"
#include "atom_class.hpp"
#include "md_coef_table.hpp"


//--- tag definitions
namespace MODEL {
    namespace DEFS{
        static const std::string ext_mol2_file  = ".mol2";
        static const std::string ext_param_file = ".param";

        static const std::string data_tag_mol2  = "@<TRIPOS>";
        static const std::string data_tag_param = "@<PARAM>";

        static const std::string scaling_LJ_tag      = "scaling_LJ";
        static const std::string scaling_coulomb_tag = "scaling_coulomb";
    }
}

//--- file loading mode
enum class MOL2_LOAD_MODE : int {
    header,
    atom,
    bond,
    angle,
    torsion,
    SUBSTRUCTURE,
    scaling,
};

//--- std::string converter for enum
namespace ENUM {

    static const std::map<std::string, MOL2_LOAD_MODE> table_str_MOL2_LOAD_MODE{
        {"MOLECULE"    , MOL2_LOAD_MODE::header      },
        {"ATOM"        , MOL2_LOAD_MODE::atom        },
        {"BOND"        , MOL2_LOAD_MODE::bond        },
        {"ANGLE"       , MOL2_LOAD_MODE::angle       },
        {"TORSION"     , MOL2_LOAD_MODE::torsion     },
        {"SUBSTRUCTURE", MOL2_LOAD_MODE::SUBSTRUCTURE},
        {"SCALING"     , MOL2_LOAD_MODE::scaling  },
    };

    static const std::map<MOL2_LOAD_MODE, std::string> table_MOL2_LOAD_MODE_str{
        {MOL2_LOAD_MODE::header      , "MOLECULE"    },
        {MOL2_LOAD_MODE::atom        , "ATOM"        },
        {MOL2_LOAD_MODE::bond        , "BOND"        },
        {MOL2_LOAD_MODE::angle       , "ANGLE"       },
        {MOL2_LOAD_MODE::torsion     , "TORSION"     },
        {MOL2_LOAD_MODE::SUBSTRUCTURE, "SUBSTRUCTURE"},
        {MOL2_LOAD_MODE::scaling     , "SCALING"     },
    };

    MOL2_LOAD_MODE which_MOL2_LOAD_MODE(const std::string &str){
        if(table_str_MOL2_LOAD_MODE.find(str) != table_str_MOL2_LOAD_MODE.end()){
            return table_str_MOL2_LOAD_MODE.at(str);
        } else {
            std::cerr << "  MOL2_LOAD_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in MOL2_LOAD_MODE.");
        }
    }

    std::string what(const MOL2_LOAD_MODE &e){
        if(table_MOL2_LOAD_MODE_str.find(e) != table_MOL2_LOAD_MODE_str.end()){
            return table_MOL2_LOAD_MODE_str.at(e);
        } else {
            using type_base = typename std::underlying_type<MOL2_LOAD_MODE>::type;
            std::cerr << "  MOL2_LOAD_MODE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in MOL2_LOAD_MODE.");
        }
    }

}

inline std::ostream& operator << (std::ostream& s, const MOL2_LOAD_MODE &e){
    s << ENUM::what(e);
    return s;
}


namespace MODEL {

    template <class Tptcl>
    void loading_mol2_atom(const std::string              &model_name,
                           const std::vector<std::string> &str_list,
                                 std::vector<Tptcl>       &atom_list,
                                 PS::S32                  &atom_count){

        //--- check format: 9 parameters, str_list[0] must be digit.
        if( str_list.size() < 9) return;
        if( !STR_TOOL::isInteger(str_list[0]) ) return;

        atom_count++;
        if( atom_count != std::stoi(str_list[0]) ){
            std::cerr << "  atom_count = " << atom_count << "\n"
                      << "  ID in file = " << str_list[0] << std::endl;
            throw std::invalid_argument("atom ID in *.mol2 file must be consective noumbers to start with 1.");
        }

        Atom_FP atom_tmp;
        atom_tmp.setAtomID(std::stoi(str_list[0])-1);
        atom_tmp.setMolID(-1);   // this is parameter table. not real atom.
        atom_tmp.setAtomType(str_list[1]);
        atom_tmp.setMolType(model_name);
        atom_tmp.setPos( PS::F64vec{std::stod(str_list[2]),
                                    std::stod(str_list[3]),
                                    std::stod(str_list[4]) } );
        atom_tmp.setCharge( Unit::coef_coulomb*std::stod(str_list[8]) );
        atom_list.push_back(atom_tmp);
    }

    template <class Tptcl>
    void loading_mol2_bond(const std::vector<std::string> &str_list,
                                 std::vector<Tptcl>       &atom_list,
                                 PS::S32                  &bond_count){

        //--- check format: 3 parameters, str_list[0], [1], and [2] must be integer.
        if( str_list.size() < 3) return;
        if( !(STR_TOOL::isInteger(str_list[0]) &&
              STR_TOOL::isInteger(str_list[1]) &&
              STR_TOOL::isInteger(str_list[2]) ) ) return;

        bond_count++;
        if( bond_count != std::stoi(str_list[0]) ){
            std::cerr << "  bond_count = " << bond_count << "\n"
                      << "  ID in file = " << str_list[0] << std::endl;
            throw std::invalid_argument("bond ID in *.mol2 file must be consective noumbers to start with 1.");
        }

        MD_DEFS::ID_type id_i = std::stoi(str_list[1]) - 1;  // ID is started from 0 in this program.
        MD_DEFS::ID_type id_j = std::stoi(str_list[2]) - 1;

        atom_list.at(id_i).bond.add(id_j);
        atom_list.at(id_j).bond.add(id_i);
    }

    template <class Tptcl>
    void loading_mol2_file(const std::string        &model_name,
                           const std::string        &file_name,
                                 std::vector<Tptcl> &atom_list){

        std::ifstream file_mol2{file_name};
        if(file_mol2.fail()) throw std::ios_base::failure("file: " + file_name + " was not found.");

        //--- loading mol2 file
        PS::S32 line_count = 0;
        PS::S32 atom_count = 0;
        PS::S32 bond_count = 0;

        PS::S32 line_index = -1;  // illigal value
        PS::S32 n_atom     = -1;
        PS::S32 n_bond     = -1;
        MOL2_LOAD_MODE mode = MOL2_LOAD_MODE::header;

        std::string line;
        while ( getline(file_mol2, line) ){
            ++line_count;

            try{
                STR_TOOL::removeCR(line);
                std::vector<std::string> str_list = STR_TOOL::split(line, " ");

                //--- skip empty or comment line(start as "!" or "//")
                if(str_list.empty()) continue;
                if(str_list[0].empty()) continue;
                if(str_list[0].substr(0,1) == "!" ||
                   str_list[0].substr(0,2) == "//") continue;

                //--- header information
                const std::string mark     = DEFS::data_tag_mol2;
                size_t            mark_len = mark.size();
                if(str_list[0].substr(0,mark_len) == mark){
                    MOL2_LOAD_MODE mode_past = mode;
                    mode = ENUM::which_MOL2_LOAD_MODE( str_list[0].substr(mark_len) );

                    if(mode == MOL2_LOAD_MODE::header){
                        line_index = line_count + 2;
                        if(mode_past != MOL2_LOAD_MODE::header){
                            throw std::length_error("*.mol2 file must contain single molecule only.");
                        }
                    } else if(mode == MOL2_LOAD_MODE::atom){
                        atom_list.clear();
                    } else if(mode == MOL2_LOAD_MODE::bond){
                        //--- clear all bonds
                        for(auto& atom : atom_list){
                            atom.bond.clear();
                        }
                    } else if(mode == MOL2_LOAD_MODE::SUBSTRUCTURE){
                        //--- ignore SUBSTRUCTURE information.
                        return;
                    }

                    continue;
                }

                //--- loading data
                switch (mode) {
                    case MOL2_LOAD_MODE::header:
                        if(line_count == line_index){
                            n_atom = std::stoi(str_list[0]);
                            n_bond = std::stoi(str_list[1]);
                        }
                    break;

                    case MOL2_LOAD_MODE::atom:
                        loading_mol2_atom(model_name, str_list, atom_list, atom_count);
                    break;

                    case MOL2_LOAD_MODE::bond:
                        loading_mol2_bond(str_list, atom_list, bond_count);
                    break;

                    default:
                        throw std::invalid_argument("undefined loading mode.");
                }

            } catch(...){
                std::cerr << "at line " << line_count << ": " << line << std::endl;
                throw;
            }
        }

        //--- check file consistency
        if(n_atom != atom_count ||
           n_bond != bond_count ){
            std::cerr << "  n_atom    = " << n_atom << " / atom_count    = " << atom_count << "\n"
                      << "  n_bond    = " << n_bond << " / bond_count    = " << bond_count << std::endl;
            throw std::invalid_argument("number of defined parameters are not match with header.");
        }
    }

    template<class Ttable_inter, class Ttable_res>
    void loading_param_atom(const std::string              &model_name,
                            const std::vector<std::string> &str_list,
                                  Ttable_inter             &atom_table,
                                  Ttable_res               &residue_table){
        MODEL::CoefAtom coef_atom;
        MODEL::KeyAtom  key_atom;

        //--- check format: 5 parameters, str_list[0] must be integer.
        if( str_list.size() < 5 ) return;

        if( str_list[1].size() > 3 ) throw std::length_error("residue name must be <= 3 characters.");

        coef_atom.mass   = 1.e-3*std::stod(str_list[2]);  // convert [g/mol] -> [kg/mol]
        coef_atom.charge = 0.0;                           // load from *.mol2 file.
        coef_atom.vdw_d  = std::stod(str_list[3]);
        coef_atom.vdw_r  = std::stod(str_list[4]);

        //--- convert into normalized value
        coef_atom.mass   = coef_atom.mass/Unit::mass_C;
        coef_atom.vdw_d  = std::sqrt(coef_atom.vdw_d);
        coef_atom.vdw_r  = coef_atom.vdw_r*0.5;

        key_atom = std::make_tuple( ENUM::which_MolName(model_name),
                                    ENUM::which_AtomName(str_list[0]) );

        //--- check duplication
        if( atom_table.find(key_atom) != atom_table.cend() ){
            std::cerr << "WARNING: 'coef_atom' of " + ENUM::what(key_atom) + " is overloaded." << std::endl;
        }

        atom_table[key_atom]    = coef_atom;
        residue_table[key_atom] = str_list[1];
    }

    template<class Ttable_bond>
    void loading_param_bond(const std::string              &model_name,
                            const std::vector<std::string> &str_list,
                                   Ttable_bond             &bond_table){
        IntraFuncForm  func_form;
        CoefBond       coef_bond;
        MODEL::KeyBond key_bond;

        //--- check format: total 6 parameters
        //--- set void value in the case of IntraFuncForm == "none".
        if( str_list.size() < 3) return;

        //--- clear coef
        coef_bond.form = IntraFuncForm::none;
        coef_bond.r0   = 1.0;
        coef_bond.k    = 0.0;
        coef_bond.a    = 0.0;

        func_form = ENUM::which_IntraFuncForm(str_list[2]);
        if(func_form == IntraFuncForm::none){
            //--- do nothing
        } else if(func_form == IntraFuncForm::harmonic){
            if( str_list.size() < 5) throw std::invalid_argument("invalid format for harmonic bond.");
            coef_bond.form = func_form;
            coef_bond.r0   = std::stod(str_list[3]);
            coef_bond.k    = std::stod(str_list[4]);
        } else if(func_form == IntraFuncForm::anharmonic){
            if( str_list.size() < 6) throw std::invalid_argument("invalid format for anharmonic bond.");
            coef_bond.form = func_form;
            coef_bond.r0   = std::stod(str_list[3]);
            coef_bond.k    = std::stod(str_list[4]);
            coef_bond.a    = std::stod(str_list[5]);
        } else {
            throw std::invalid_argument("undefined form of bond potential.");
        }

        key_bond = std::make_tuple( ENUM::which_MolName(model_name),
                                    ENUM::which_AtomName(str_list[0]),
                                    ENUM::which_AtomName(str_list[1]) );
        //--- check duplication
        if( bond_table.find(key_bond) != bond_table.cend() ){
            std::cerr << "WARNING: 'coef_bond' of " << ENUM::what(key_bond) << " is overloaded." << std::endl;
        }
        bond_table[key_bond] = coef_bond;

        //--- add reverse shape
        key_bond = reverse_shape(key_bond);
        bond_table[key_bond] = coef_bond;
    }

    template<class Ttable_angle>
    void loading_param_angle(const std::string              &model_name,
                             const std::vector<std::string> &str_list,
                                   Ttable_angle             &angle_table){
        IntraFuncForm   func_form;
        CoefAngle       coef_angle;
        MODEL::KeyAngle key_angle;

        //--- check format: total 6 parameters
        //--- set void value in the case of IntraFuncForm == "none".
        if( str_list.size() < 4) return;

        //--- clear coef
        coef_angle.form = IntraFuncForm::none;
        coef_angle.theta0 = 0.0;
        coef_angle.k      = 0.0;

        func_form = ENUM::which_IntraFuncForm(str_list[3]);
        if( func_form == IntraFuncForm::none ){
            //--- do nothing
        } else if(func_form == IntraFuncForm::harmonic){
            if( str_list.size() < 6 ) throw std::invalid_argument("invalid format for angle potential parameters.");

            coef_angle.form   = ENUM::which_IntraFuncForm(str_list[3]);
            coef_angle.theta0 = std::stod(str_list[4])*Unit::pi/180.0;                           // [degree] -> [rad]
            coef_angle.k      = std::stod(str_list[5])/std::pow(std::sin(coef_angle.theta0),2);  // make [kcal/mol rad^2]/sin(theta0)^2
        } else {
            throw std::invalid_argument("undefined form of angle potential.");
        }

        key_angle = std::make_tuple( ENUM::which_MolName(model_name),
                                     ENUM::which_AtomName(str_list[0]),
                                     ENUM::which_AtomName(str_list[1]),
                                     ENUM::which_AtomName(str_list[2]) );

        //--- check duplication
        if( angle_table.find(key_angle) != angle_table.cend() ){
            std::cerr << "WARNING: 'coef_angle' of " << ENUM::what(key_angle) << " is overloaded." << std::endl;
        }
        angle_table[key_angle] = coef_angle;

        //--- add reverse shape
        key_angle = reverse_shape(key_angle);
        angle_table[key_angle] = coef_angle;
    }

    template<class Ttable_torsion>
    void loading_param_torsion(const std::string              &model_name,
                               const std::vector<std::string> &str_list,
                                     Ttable_torsion           &torsion_table){
        IntraFuncForm     func_form;
        CoefTorsion       coef_torsion;
        MODEL::KeyTorsion key_torsion;

        //--- check format: total 9 parameters,
        //--- set void value in the case of IntraFuncForm == "none".
        if( str_list.size() < 6) return;

        //--- clear coef
        coef_torsion.form   = IntraFuncForm::none;
        coef_torsion.theta0 = 0.0;
        coef_torsion.k      = 0.0;
        coef_torsion.k2     = 0.0;
        coef_torsion.k3     = 0.0;
        coef_torsion.n_min  = 0;

        func_form = ENUM::which_IntraFuncForm(str_list[5]);
        if( func_form == IntraFuncForm::none ){
            // do nothing
        } else if( func_form == IntraFuncForm::cos ){
            if( str_list.size() < 9 ){
                throw std::invalid_argument("invalid format for cos form torsion potential.");
            }
            if( !STR_TOOL::isInteger(str_list[8]) ){
                throw std::invalid_argument("'n_min' value must be a integer.");
            }
            coef_torsion.form   = func_form;
            coef_torsion.theta0 = std::stod(str_list[6])*Unit::pi/180.0;   // [degree] -> [rad]
            coef_torsion.k      = std::stod(str_list[7]);
            coef_torsion.n_min  = std::stoi(str_list[8]);
        } else if( func_form == IntraFuncForm::OPLS_3 ){
            if( str_list.size() < 9 ){
                throw std::invalid_argument("invalid format for OPLS_3 form torsion potential.");
            }
            coef_torsion.form   = func_form;
            coef_torsion.k      = 0.5*std::stod(str_list[6]);
            coef_torsion.k2     = 0.5*std::stod(str_list[7]);
            coef_torsion.k3     = 0.5*std::stod(str_list[8]);
        } else {
            throw std::invalid_argument("undefined form of torsion potential.");
        }

        key_torsion = std::make_tuple( ENUM::which_MolName(model_name),
                                       ENUM::which_TorsionShape(str_list[0]),
                                       ENUM::which_AtomName(str_list[1]),
                                       ENUM::which_AtomName(str_list[2]),
                                       ENUM::which_AtomName(str_list[3]),
                                       ENUM::which_AtomName(str_list[4])     );

        //--- check duplication
        if( torsion_table.find(key_torsion) != torsion_table.cend() ){
            std::cerr << "WARNING: 'coef_torsion' of " << ENUM::what(key_torsion) << " is overloaded." << std::endl;
        }
        torsion_table[key_torsion] = coef_torsion;

        //--- add other shape of same pair
        if(        ENUM::which_TorsionShape(str_list[0]) == TorsionShape::dihedral ){
            key_torsion = reverse_shape(key_torsion);
            torsion_table[key_torsion] = coef_torsion;
        } else if( ENUM::which_TorsionShape(str_list[0]) == TorsionShape::improper ){
            const auto key_list = other_shape(key_torsion);
            for(const auto& key : key_list){
                torsion_table[key] = coef_torsion;
            }
        }
    }

    template <class Ttable_scaling>
    void loading_param_scaling(const std::string              &model_name,
                               const std::vector<std::string> &str_list,
                                     Ttable_scaling           &scaling_table){

        //--- check format: total 2 parameters in line.
        if(str_list.size() < 2) return;

        //--- get reference
        MolName key_scaling = ENUM::which_MolName(model_name);
        auto &coef_scaling  = scaling_table[key_scaling];

        //--- resize mask list
        size_t order = 0;
        for(size_t i=1; i<str_list.size(); ++i){
            if( !STR_TOOL::isNumeric(str_list[i]) )break;
            ++order;
        }
        if(order > coef_scaling.size()) coef_scaling.resize(order);

        //--- set value
        for(size_t i=0; i<order; ++i){
            if(       str_list[0] == DEFS::scaling_LJ_tag){
                coef_scaling[i].scale_LJ      = std::stof(str_list[i+1]);
            } else if(str_list[0] == DEFS::scaling_coulomb_tag){
                coef_scaling[i].scale_coulomb = std::stof(str_list[i+1]);
            }
        }
    }

    template<class Ttable_atom, class Ttable_res,
             class Ttable_bond, class Ttable_angle, class Ttable_torsion,
             class Ttable_scaling>
    void loading_param_file(const std::string    &model_name,
                            const std::string    &file_name,
                                  Ttable_atom    &atom_table,
                                  Ttable_res     &residue_table,
                                  Ttable_bond    &bond_table,
                                  Ttable_angle   &angle_table,
                                  Ttable_torsion &torsion_table,
                                  Ttable_scaling &scaling_table ){

        std::ifstream file_para{file_name};
        if(file_para.fail()) throw std::ios_base::failure("file: " + file_name + ".para was not found.");

        //--- loading para file
        MOL2_LOAD_MODE mode = MOL2_LOAD_MODE::header;

        std::string line;
        PS::S32     line_count = 0;
        while (getline(file_para, line)){
            ++line_count;

            try{
                STR_TOOL::removeCR(line);
                std::vector<std::string> str_list = STR_TOOL::split(line, " ");

                //--- skip empty or comment line(start as "!" or "//")
                if(str_list.empty()) continue;
                if(str_list[0].empty()) continue;
                if(str_list[0].substr(0,1) == "!" ||
                   str_list[0].substr(0,2) == "//"  ) continue;

                //--- header information
                const std::string mark     = DEFS::data_tag_param;
                size_t            mark_len = mark.size();
                if( str_list[0].substr(0,mark_len) == mark ){
                    mode = ENUM::which_MOL2_LOAD_MODE( str_list[0].substr(mark_len) );
                    continue;
                }

                //--- loading data
                switch (mode) {
                    case MOL2_LOAD_MODE::header:
                        continue;
                    break;

                    case MOL2_LOAD_MODE::atom:
                        loading_param_atom(model_name, str_list, atom_table, residue_table);
                    break;

                    case MOL2_LOAD_MODE::bond:
                        loading_param_bond(model_name, str_list, bond_table);
                    break;

                    case MOL2_LOAD_MODE::angle:
                        loading_param_angle(model_name, str_list, angle_table);
                    break;

                    case MOL2_LOAD_MODE::torsion:
                        loading_param_torsion(model_name, str_list, torsion_table);
                    break;

                    case MOL2_LOAD_MODE::scaling:
                        loading_param_scaling(model_name, str_list, scaling_table);
                    break;

                    default:
                        throw std::invalid_argument("undefined loading mode.");
                }

            } catch(...){
                std::cerr << "at line " << line_count << ": " << line << std::endl;
                throw;
            }
        }

        //--- check scaling parameter is set or not
        if( scaling_table.count( ENUM::which_MolName(model_name) ) == 0 ){
            throw std::invalid_argument("scaling parameter was not found.");
        }
    }

    template<class Tptcl, class TCoefTable>
    void loading_model_parameter(const std::string        &model_name,
                                       std::vector<Tptcl> &atom_list,
                                       TCoefTable         &coef_table){

        std::string file_name;

        //--- temporary table for intermolecular parameter
        std::unordered_map< KeyAtom,
                            CoefAtom,
                            hash_tuple::hash_func<KeyAtom>> atom_table;

        //--- setting the load factor of std::unordered_map
        atom_table.max_load_factor(0.7);

        //--- loading ****.mol2 file
        if(model_name.find_first_of("/") != std::string::npos){
            file_name = model_name + DEFS::ext_mol2_file;
        } else {
            file_name = MD_DEFS::model_dir + model_name + DEFS::ext_mol2_file;
        }
        try{
            loading_mol2_file(model_name,
                              file_name,
                              atom_list );
        }
        catch(...){
            std::cerr << "ERROR at loading file: " << file_name << std::endl;
            throw;
        }

        //--- loading ****.param file
        if(model_name.find_first_of("/") != std::string::npos){
            file_name = model_name + DEFS::ext_param_file;
        } else {
            file_name = MD_DEFS::model_dir + model_name + DEFS::ext_param_file;
        }
        try{
            loading_param_file(model_name,
                               file_name,
                               atom_table,
                               coef_table.residue,
                               coef_table.bond,
                               coef_table.angle,
                               coef_table.torsion,
                               coef_table.mask_scaling);
        }
        catch(...){
            std::cerr << "ERROR at loading file: " << file_name << std::endl;
            throw;
        }

        //--- copy VDW coef to atom_list from inter_table
        std::vector<KeyAtom> missed_in_param;
        missed_in_param.clear();
        for(size_t index=0; index<atom_list.size(); index++){
            KeyAtom key_atom = std::make_tuple( atom_list[index].getMolType(),
                                                atom_list[index].getAtomType() );
            if(atom_table.count(key_atom) == 1){
                auto coef_atom = atom_table[key_atom];
                atom_list[index].setMass(  coef_atom.mass );
                atom_list[index].setVDW_D( coef_atom.vdw_d );
                atom_list[index].setVDW_R( coef_atom.vdw_r );
            } else {
                missed_in_param.push_back(key_atom);
            }
        }

        if(missed_in_param.size() > 0){
            std::ostringstream oss;
            oss << "\n";
            for(const auto& key_atom : missed_in_param){
                oss << "  element parameter of " << ENUM::what(std::get<1>(key_atom))
                    << " was not defined in '" << file_name
                    << "' file." << "\n";
            }
            throw std::invalid_argument(oss.str());
        }

        //--- check consistency
        missed_in_param.clear();
        for(auto itr = atom_table.begin(); itr != atom_table.end(); itr++){
            if(std::get<0>(itr->first) != ENUM::which_MolName(model_name)) continue;

            bool find_flag = false;
            for(auto &atom : atom_list){
                if(std::make_tuple(atom.getMolType(), atom.getAtomType()) == itr->first){
                    find_flag = true;
                    break;
                }
            }

            if( !find_flag ){
                missed_in_param.push_back(itr->first);
            }
        }

        if(missed_in_param.size() > 0){
            std::ostringstream oss;
            oss << "\n";
            for(const auto& key_atom : missed_in_param){
                oss << "  atom: " << ENUM::what( std::get<1>(key_atom) )
                    << " is not defined in "
                    << model_name + ".mol2 file."  << "\n";
            }
            throw std::invalid_argument(oss.str());
        }
    }


}
