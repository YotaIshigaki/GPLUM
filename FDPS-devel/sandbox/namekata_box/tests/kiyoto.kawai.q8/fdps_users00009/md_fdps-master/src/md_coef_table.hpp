//***************************************************************************************
//  This is the coefficient data table class for calculate interactions.
//***************************************************************************************
#pragma once

#include <string>
#include <sstream>
#include <unordered_map>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "unit.hpp"
#include "md_enum.hpp"
#include "md_defs.hpp"
#include "atom_class_base.hpp"

//--- grobal object of parameter table
namespace MODEL {

    //--- parameter for intermolecular interactions
    struct CoefAtom {
      public:
        PS::F32 mass;
        PS::F32 charge;
        PS::F32 vdw_d;
        PS::F32 vdw_r;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "mass   : "        << std::setw(12) << this->mass
                                        << "  | real value: " << std::setw(12) << this->mass*Unit::mass_C         << "\n";
            oss << std::setw(shift + 9) << "charge : "        << std::setw(12) << this->charge
                                        << "  | real value: " << std::setw(12) << this->charge/Unit::coef_coulomb << "\n";
            oss << std::setw(shift + 9) << "vdw_d  : "        << std::setw(12) << this->vdw_d
                                        << "  | real value: " << std::setw(12) << this->vdw_d*this->vdw_d         << "\n";
            oss << std::setw(shift + 9) << "vdw_r  : "        << std::setw(12) << this->vdw_r
                                        << "  | real value: " << std::setw(12) << this->vdw_r*2.0                 << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str(shift);
        }
    };

    //--- parameter for intramolecular interactions
    struct CoefBond {
      public:
        IntraFuncForm form;
        PS::F32       r0, k, a;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 7) << "form : "                         << this->form << "\n";
            oss << std::setw(shift + 7) << "r0   : " << std::setprecision(8) << this->r0   << "\n";
            oss << std::setw(shift + 7) << "k    : " << std::setprecision(8) << this->k    << "\n";
            oss << std::setw(shift + 7) << "a    : " << std::setprecision(8) << this->a    << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str(shift);
        }
    };

    struct CoefAngle {
      public:
        IntraFuncForm form;
        PS::F32       theta0, k;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "form   : "                         << this->form   << "\n";
            oss << std::setw(shift + 9) << "k      : " << std::setprecision(8) << this->k      << "\n";
            oss << std::setw(shift + 9) << "theta0 : " << std::setprecision(8) << this->theta0 << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str();
        }
    };

    struct CoefTorsion {
      public:
        IntraFuncForm form;
        PS::S32       n_min;
        PS::F32       k, k2, k3;
        PS::F32       theta0;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "form   : "                         << this->form   << "\n";
            oss << std::setw(shift + 9) << "n_min  : " << std::setprecision(8) << this->n_min  << "\n";
            oss << std::setw(shift + 9) << "k      : " << std::setprecision(8) << this->k      << "\n";
            oss << std::setw(shift + 9) << "k2     : " << std::setprecision(8) << this->k2     << "\n";
            oss << std::setw(shift + 9) << "k3     : " << std::setprecision(8) << this->k3     << "\n";
            oss << std::setw(shift + 9) << "theta0 : " << std::setprecision(8) << this->theta0 << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str();
        }
    };

    //--- for residue information
    //------ key = (model_name, atom_name)
    using KeyAtom = std::tuple<MolName,
                               AtomName>;

    //--- for intramolecular parameter
    //------ key = (model_name, i_atom, j_atom)
    using KeyBond = std::tuple<MolName,
                               AtomName,
                               AtomName>;

    //------ key = (model_name, i_atom, j_atom, k_atom)
    using KeyAngle = std::tuple<MolName,
                                AtomName,
                                AtomName,
                                AtomName>;

    //------ key = (model_name, torsion_shape, i_atom, j_atom, k_atom, l_atom)
    using KeyTorsion = std::tuple<MolName,
                                  TorsionShape,
                                  AtomName,
                                  AtomName,
                                  AtomName,
                                  AtomName    >;


    //--- transform intra pair key
    KeyBond reverse_shape(const KeyBond &key){
        return KeyBond{ std::get<0>(key),
                        std::get<2>(key),
                        std::get<1>(key) };
    }
    KeyAngle reverse_shape(const KeyAngle &key){
        return KeyAngle{ std::get<0>(key),
                         std::get<3>(key),
                         std::get<2>(key),
                         std::get<1>(key) };
    }
    KeyTorsion reverse_shape(const KeyTorsion &key){
        if(std::get<1>(key) != TorsionShape::dihedral) throw std::invalid_argument("inverseShape() must be used for dihedral shape only.");
        return KeyTorsion{ std::get<0>(key),
                           std::get<1>(key),
                           std::get<5>(key),
                           std::get<4>(key),
                           std::get<3>(key),
                           std::get<2>(key) };
    }
    MD_EXT::fixed_vector<KeyTorsion, 3> other_shape(const KeyTorsion &key){
        if(std::get<1>(key) != TorsionShape::improper) throw std::invalid_argument("otherShape() must be used for improper shape only.");

        MD_EXT::fixed_vector<KeyTorsion, 3> result;
        result.clear();

        //--- swap edge
        result.push_back( KeyTorsion{ std::get<0>(key),
                                      std::get<1>(key),
                                      std::get<5>(key),
                                      std::get<3>(key),
                                      std::get<4>(key),
                                      std::get<2>(key) } );

        //--- another axis
        result.push_back( KeyTorsion{ std::get<0>(key),
                                      std::get<1>(key),
                                      std::get<2>(key),
                                      std::get<4>(key),
                                      std::get<3>(key),
                                      std::get<5>(key) } );

        //--- another axis + swap edge
        result.push_back( KeyTorsion{ std::get<0>(key),
                                      std::get<1>(key),
                                      std::get<5>(key),
                                      std::get<4>(key),
                                      std::get<3>(key),
                                      std::get<2>(key) } );
        return result;
    }

    namespace _Impl {
        class CoefTable {
        public:
            std::unordered_map< MolName,
                                MD_DEFS::MaskList,
                                std::hash<MolName> > mask_scaling;

            std::unordered_map< KeyAtom,
                                std::string,
                                hash_tuple::hash_func<KeyAtom>> residue;

            std::unordered_map< KeyBond,
                                CoefBond,
                                hash_tuple::hash_func<KeyBond>> bond;

            std::unordered_map< KeyAngle,
                                CoefAngle,
                                hash_tuple::hash_func<KeyAngle>> angle;

            std::unordered_map< KeyTorsion,
                                CoefTorsion,
                                hash_tuple::hash_func<KeyTorsion>> torsion;

            CoefTable(){
                const PS::F32 factor = 0.7;
                this->mask_scaling.max_load_factor(factor);
                this->residue.max_load_factor(factor);
                this->residue.max_load_factor(factor);
                this->bond.max_load_factor(factor);
                this->angle.max_load_factor(factor);
                this->torsion.max_load_factor(factor);
            }
            void broadcast(const PS::S32 root = 0){
                COMM_TOOL::broadcast(this->mask_scaling, root);
                COMM_TOOL::broadcast(this->residue     , root);
                COMM_TOOL::broadcast(this->bond        , root);
                COMM_TOOL::broadcast(this->angle       , root);
                COMM_TOOL::broadcast(this->torsion     , root);
            }
            void clear(){
                this->mask_scaling.clear();
                this->residue.clear();
                this->residue.clear();
                this->bond.clear();
                this->angle.clear();
                this->torsion.clear();
            }
        };
    }
    static _Impl::CoefTable coef_table;


    //--- display interface
    template<class Tptcl>
    void print_model_template(const std::vector<Tptcl> &atom_list){
        std::ostringstream oss;

        for(size_t index=0; index<atom_list.size(); ++index){
            const auto& atom = atom_list.at(index);

            oss << "  index = " << index << "\n";
            oss << atom.str();
            oss << "\n";
        }
        oss << "\n";

        std::cout << oss.str();
    }

    //--- for Coef****::to_str() in md_coef_table.hpp.
    template<class Ttable>
    void print_coef_table(const Ttable &table, const MolName &model){
        std::string str;
        size_t count = 0;
        for(auto &coef : table){
            if( std::get<0>(coef.first) != model ) continue;

            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += coef.second.to_str(6);
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }

    template<class Ttable>
    void print_coef_table(const Ttable &table){

        std::string str;
        size_t count = 0;
        for(auto &coef : table){
            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += coef.second.to_str(6);
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }

    //--- overload for std::unordered_map<key__, std::string>
    template<class Tkey, class Hash, class Pred, class Allocator>
    void print_coef_table( const std::unordered_map<Tkey, std::string, Hash, Pred, Allocator> &table,
                           const MolName &model ){
        std::string str;
        size_t count = 0;
        for(auto &coef : table){
            if( std::get<0>(coef.first) != model ) continue;

            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += std::string(6, ' ') + coef.second;
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }

    template<class Tkey, class Hash, class Pred, class Allocator>
    void print_coef_table( const std::unordered_map<Tkey, std::string, Hash, Pred, Allocator> &table ){
        std::string str;
        size_t count = 0;
        for(auto &coef : table){
            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += std::string(6, ' ') + coef.second;
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }

}
