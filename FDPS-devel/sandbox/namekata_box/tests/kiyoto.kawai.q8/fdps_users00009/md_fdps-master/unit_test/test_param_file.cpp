//***************************************************************************************
//  This is ***.param file tester.
//***************************************************************************************

#include <tuple>
#include <iostream>
#include <unordered_map>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

//--- user defined headers
//------ definition of data set
#include "unit.hpp"
#include "md_enum.hpp"
#include "atom_class.hpp"
#include "md_coef_table.hpp"
//------ initializer
#include "initialize.hpp"


enum class ParamCheck : int {
    not_defined,
    not_used,
    OK,
};

template <class Tpsys, class Ttree>
void checkIntraForceParam(const Tpsys &test_atom, Ttree &tree){
    std::unordered_map<MODEL::KeyBond,    ParamCheck, hash_tuple::hash_func<MODEL::KeyBond>>    bond_param_check;
    std::unordered_map<MODEL::KeyAngle,   ParamCheck, hash_tuple::hash_func<MODEL::KeyAngle>>   angle_param_check;
    std::unordered_map<MODEL::KeyTorsion, ParamCheck, hash_tuple::hash_func<MODEL::KeyTorsion>> torsion_param_check;

    //--- initialize check flag
    for(auto &bond_param : MODEL::coef_table.bond){
        MODEL::KeyBond key_bond = bond_param.first;
        bond_param_check[key_bond] = ParamCheck::not_used;
    }
    for(auto &angle_param : MODEL::coef_table.angle){
        MODEL::KeyAngle key_angle = angle_param.first;
        angle_param_check[key_angle] = ParamCheck::not_used;
    }
    for(auto &torsion_param : MODEL::coef_table.torsion){
        MODEL::KeyTorsion key_torsion = torsion_param.first;
        torsion_param_check[key_torsion] = ParamCheck::not_used;
    }

    //--- check existance of parameters.
    for(PS::S32 i=0; i<test_atom.getNumberOfParticleLocal(); ++i){
        const auto  ep_i          = test_atom[i];
        const auto& bond_list     = ep_i.bond;
        const auto& angle_list    = ep_i.angle_list();
        const auto& dihedral_list = ep_i.dihedral_list();
        const auto& improper_list = ep_i.improper_list();

        //--- check bond parameters
        for(size_t i=0; i<bond_list.size(); ++i){
            MD_DEFS::ID_type id_j  = bond_list.at(i);
            EP_intra*        ptr_j = tree.getEpjFromId(id_j);
            if(ptr_j == nullptr) throw std::invalid_argument("EPJ is not found.");
            const auto       ep_j = *ptr_j;

            MODEL::KeyBond  key_bond = std::make_tuple(ep_i.getMolType(),
                                                       ep_i.getAtomType(),
                                                       ep_j.getAtomType());

            if(MODEL::coef_table.bond.count(key_bond) == 0){
                //--- param was not found
                bond_param_check[key_bond] = ParamCheck::not_defined;
            } else {
                //--- param was found
                bond_param_check[key_bond] = ParamCheck::OK;
            }
        }

        //--- check angle parameters
        for(size_t i=0; i<angle_list.size(); ++i){
            MD_DEFS::ID_type id_i = std::get<0>( angle_list.at(i) );
            MD_DEFS::ID_type id_j = std::get<1>( angle_list.at(i) );
            MD_DEFS::ID_type id_k = std::get<2>( angle_list.at(i) );

            EP_intra* ptr_i = tree.getEpjFromId(id_i);
            EP_intra* ptr_j = tree.getEpjFromId(id_j);
            EP_intra* ptr_k = tree.getEpjFromId(id_k);

            if(ptr_i == nullptr) throw std::invalid_argument("EPI is not found.");
            if(ptr_j == nullptr) throw std::invalid_argument("EPJ is not found.");
            if(ptr_k == nullptr) throw std::invalid_argument("EPK is not found.");

            const auto ep_i = *ptr_i;
            const auto ep_j = *ptr_j;
            const auto ep_k = *ptr_k;

            MODEL::KeyAngle  key_angle = std::make_tuple(ep_i.getMolType(),
                                                         ep_i.getAtomType(),
                                                         ep_j.getAtomType(),
                                                         ep_k.getAtomType() );

            if(MODEL::coef_table.angle.count(key_angle) == 0){
                //--- param was not found
                angle_param_check[key_angle] = ParamCheck::not_defined;
            } else {
                //--- param was found
                angle_param_check[key_angle] = ParamCheck::OK;
            }
        }

        //--- check torsion parameters
        //------ dihedral torsion potential
        for(size_t i=0; i<dihedral_list.size(); ++i){
            MD_DEFS::ID_type id_i = std::get<0>( dihedral_list.at(i) );
            MD_DEFS::ID_type id_j = std::get<1>( dihedral_list.at(i) );
            MD_DEFS::ID_type id_k = std::get<2>( dihedral_list.at(i) );
            MD_DEFS::ID_type id_l = std::get<3>( dihedral_list.at(i) );

            EP_intra* ptr_i = tree.getEpjFromId(id_i);
            EP_intra* ptr_j = tree.getEpjFromId(id_j);
            EP_intra* ptr_k = tree.getEpjFromId(id_k);
            EP_intra* ptr_l = tree.getEpjFromId(id_l);

            if(ptr_i == nullptr) throw std::invalid_argument("EPI is not found.");
            if(ptr_j == nullptr) throw std::invalid_argument("EPJ is not found.");
            if(ptr_k == nullptr) throw std::invalid_argument("EPK is not found.");
            if(ptr_l == nullptr) throw std::invalid_argument("EPL is not found.");

            const auto ep_i = *ptr_i;
            const auto ep_j = *ptr_j;
            const auto ep_K = *ptr_k;
            const auto ep_l = *ptr_l;

            MODEL::KeyTorsion  key_torsion = std::make_tuple(ep_i.getMolType(),
                                                             TorsionShape::dihedral,
                                                             ep_i.getAtomType(),
                                                             ep_j.getAtomType(),
                                                             ep_K.getAtomType(),
                                                             ep_l.getAtomType() );

            if(MODEL::coef_table.torsion.count(key_torsion) == 0){
                //--- param was not found
                torsion_param_check[key_torsion] = ParamCheck::not_defined;
            } else {
                //--- param was found
                torsion_param_check[key_torsion] = ParamCheck::OK;
            }
        }
        //------ improper torsion potential
        for(size_t i=0; i<improper_list.size(); ++i){
            MD_DEFS::ID_type id_i = std::get<0>( improper_list.at(i) );
            MD_DEFS::ID_type id_j = std::get<1>( improper_list.at(i) );
            MD_DEFS::ID_type id_k = std::get<2>( improper_list.at(i) );
            MD_DEFS::ID_type id_l = std::get<3>( improper_list.at(i) );

            EP_intra* ptr_i = tree.getEpjFromId(id_i);
            EP_intra* ptr_j = tree.getEpjFromId(id_j);
            EP_intra* ptr_k = tree.getEpjFromId(id_k);
            EP_intra* ptr_l = tree.getEpjFromId(id_l);

            if(ptr_i == nullptr) throw std::invalid_argument("EPI is not found.");
            if(ptr_j == nullptr) throw std::invalid_argument("EPJ is not found.");
            if(ptr_k == nullptr) throw std::invalid_argument("EPK is not found.");
            if(ptr_l == nullptr) throw std::invalid_argument("EPL is not found.");

            const auto ep_i = *ptr_i;
            const auto ep_j = *ptr_j;
            const auto ep_K = *ptr_k;
            const auto ep_l = *ptr_l;

            MODEL::KeyTorsion  key_torsion = std::make_tuple(ep_i.getMolType(),
                                                             TorsionShape::improper,
                                                             ep_i.getAtomType(),
                                                             ep_j.getAtomType(),
                                                             ep_K.getAtomType(),
                                                             ep_l.getAtomType() );

            if(MODEL::coef_table.torsion.count(key_torsion) == 0){
                //--- param was not found
                torsion_param_check[key_torsion] = ParamCheck::not_defined;
            } else {
                //--- param was found
                torsion_param_check[key_torsion] = ParamCheck::OK;
            }
        }
    }

    //--- report results
    //------ bond parameters
    std::cout << "\n"
              << "  " << MODEL::coef_table.bond.size() << " parameters for bond." << std::endl;
    for(const auto &bond : bond_param_check){
        if(         bond.second == ParamCheck::not_defined){
            std::cout << "    ERROR: bond parameter of " << ENUM::what(bond.first) << " was not defined." << std::endl;
        } else if ( bond.second == ParamCheck::not_used){
            //--- must be used by both side
            std::cout << "    WARNING: bond parameter of " << ENUM::what(bond.first) << " was not used." << std::endl;
        }
    }
    //------ angle parameters
    std::cout << "  " << MODEL::coef_table.angle.size() << " parameters for angle." << std::endl;
    for(const auto &angle : angle_param_check){
        if(         angle.second == ParamCheck::not_defined){
            std::cout << "    ERROR: angle parameter of " << ENUM::what(angle.first) << " was not defined." << std::endl;
        } else if ( angle.second == ParamCheck::not_used){
            //--- check opposite side is used or not
            auto opposite_key = angle.first;
            std::get<1>(opposite_key) = std::get<3>(angle.first);
            std::get<3>(opposite_key) = std::get<1>(angle.first);
            if( angle_param_check.count(opposite_key) == 0 ){
                std::cout << "    WARNING: angle parameter of " << ENUM::what(angle.first) << " was not defined. (opposite side)" << std::endl;
            } else if( angle_param_check[opposite_key] == ParamCheck::not_used ){
                std::cout << "    WARNING: angle parameter of " << ENUM::what(angle.first) << " was not used." << std::endl;
            }
        }
    }
    //------ torsion parameters
    std::cout << "  " << MODEL::coef_table.torsion.size() << " parameters for torsion." << std::endl;
    for(const auto &torsion : torsion_param_check){
        if(         torsion.second == ParamCheck::not_defined ){
            std::cout << "    ERROR: torsion parameter of " << ENUM::what(torsion.first) << " was not defined." << std::endl;
        } else if ( torsion.second == ParamCheck::not_used ){
            //--- check opposite side is used or not
            auto opposite_key = torsion.first;
            std::get<2>(opposite_key) = std::get<5>(torsion.first);
            std::get<5>(opposite_key) = std::get<2>(torsion.first);
            if( torsion_param_check.count(opposite_key) == 0 ){
                std::cout << "    WARNING: torsion parameter of " << ENUM::what(torsion.first) << " was not defined. (opposite side)" << std::endl;
            } else if( torsion_param_check[opposite_key] == ParamCheck::not_used ){
                std::cout << "    WARNING: torsion parameter of " << ENUM::what(torsion.first) << " was not used." << std::endl;
            }
        }
    }

}


void test_init(){
    if(PS::Comm::getRank() != 0) return;

    //--- system parameters
    System::profile.coef_ema      = 0.3;
    System::profile.theta         = 0.5;
    System::profile.n_leaf_limit  = 8;
    System::profile.n_group_limit = 64;
    System::profile.cycle_dinfo   = 1;

    //--- set loop condition
    System::profile.dt       = 1.0;
    System::profile.istep    = 0;
    System::profile.nstep_st = 0;
    System::profile.nstep_ed = 1;

    //--- for cut_off radius
    System::profile.cut_off_LJ    = 12.0;
    System::profile.cut_off_intra = 12.0;

    //--- for installing molecule at initialize
    System::profile.ex_radius = 3.5;
    System::profile.try_limit = 100;

    //--- set domain size
    Normalize::setBoxSize( PS::F32vec{ 100.0,
                                       100.0,
                                       100.0 } );

}


int main(int argc, char *argv[]) {

    PS::Initialize(argc, argv);

    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> test_atom;

    test_atom.initialize();
    test_atom.setNumberOfParticleLocal(0);

    test_init();

    if(PS::Comm::getRank() == 0) {
        //--- molecular model loading test
        std::cout << "TEST: check consistensy between ***.mol2 file and ***.param file..." << std::endl;

        if(argc <= 1){
            std::cout << "    no input." << std::endl
                      << "    usage:  $ ./test_model.x [model_file_name]";
        }

        System::model_list.clear();
        for(int i=1; i<argc; ++i){
            MolName model = ENUM::which_MolName(argv[i]);
            System::model_list.push_back( std::make_pair(model, 0) );
        }
        System::model_template.resize(System::model_list.size());

        for(size_t i=0; i<System::model_list.size(); ++i){

            //--- clear model data
            System::model_template.at(i).clear();
            MODEL::coef_table.clear();

            //--- load test model
            MolName     model      = System::model_list.at(i).first;
            std::string model_name = ENUM::what(model);
            std::cout << "\n"
                      << "  checking model file: " << model_name << "\n" << std::endl;

            MODEL::loading_model_parameter(model_name,
                                           System::model_template.at(i),
                                           MODEL::coef_table            );

            //--- insert test model
            test_atom.setNumberOfParticleLocal(0);
            for(size_t j=0; j<System::model_list.size(); ++j){
                size_t n = 0;
                if(j == i){
                    n = 1;
                }
                System::model_list.at(j).second = n;
            }
            Initialize::InitParticle(test_atom,
                                     System::model_list,
                                     System::model_template,
                                     System::get_ex_radius(),
                                     System::get_try_limit(),
                                     0.0);

            //--- tree data for use as temporary
            EP_intra::setRcut( Normalize::normCutOff(9.0) );
            PS::TreeForForceShort<ForceIntra<PS::F32>,
                                  EP_intra,
                                  EP_intra>::Scatter test_tree;
            test_tree.initialize(System::model_template.at(i).size(), 0.5, 8, 64);
            test_tree.calcForceAll( IntraPair::dummy_func{},
                                    test_atom,
                                    dinfo );

            //--- make intra pair list
            struct GetBond {
                MD_EXT::basic_connect<MD_DEFS::ID_type,
                                      MD_DEFS::max_bond> operator () (const AtomConnect &atom){
                    return atom.bond;
                }
            };
            IntraPair::IntraMaskMaker<  MD_DEFS::ID_type, GetBond> intra_mask_maker;
            IntraPair::AngleListMaker<  MD_DEFS::ID_type, GetBond> angle_list_maker;
            IntraPair::TorsionListMaker<MD_DEFS::ID_type, GetBond> torsion_list_maker;
            for(PS::S64 i=0; i<test_atom.getNumberOfParticleLocal(); ++i){
                test_atom[i].clear_intra_list();
                intra_mask_maker(  test_atom[i], test_tree, MODEL::coef_table.mask_scaling.at(model), test_atom[i].mask_list());
                angle_list_maker(  test_atom[i], test_tree, test_atom[i].angle_list());
                torsion_list_maker(test_atom[i], test_tree, test_atom[i].dihedral_list(), test_atom[i].improper_list());
            }

            //--- check intraforce coef table.
            checkIntraForceParam(test_atom, test_tree);
        }

        std::cout << "\n  the test finished." << std::endl;
    }

    PS::Finalize();
    return 0;
}
