//***************************************************************************************
//  This routine is initializer for "md_fdps_main.cpp".
//***************************************************************************************
#pragma once

#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <memory>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "unit.hpp"
#include "md_defs.hpp"
#include "md_coef_table.hpp"
#include "atom_class.hpp"
#include "intra_pair.hpp"
#include "md_loading_condition.hpp"
#include "md_loading_model.hpp"

namespace Initialize {

    //--- adjast in [0.0,1.0) space
    template <class Tf>
    PS::Vector3<Tf> periodicAdjustNorm_round(const PS::Vector3<Tf> pos_norm){
        PS::Vector3<Tf> pos_new;
        pos_new.x = pos_norm.x - std::floor(pos_norm.x);
        pos_new.y = pos_norm.y - std::floor(pos_norm.y);
        pos_new.z = pos_norm.z - std::floor(pos_norm.z);
        return pos_new;
    }

    template <class Tf>
    PS::Vector3<Tf> periodicAdjustReal_round(const PS::Vector3<Tf> pos_real){
        PS::Vector3<Tf> pos_new;
        pos_new = Normalize::normPos(pos_new);
        pos_new = periodicAdjustNorm_round(pos_new);
        pos_new = Normalize::realPos(pos_new);
        return pos_new;
    }

    //--- building mask for collision in initialize
    template <class Tptcl>
    void make_collision_mask(      std::vector<MD_DEFS::MaskList> &collision_mask,
                             const std::vector<Tptcl>             &model_template,
                             const PS::F64                         ex_r_real      ){

        collision_mask.clear();
        const PS::F64 ex_r2 = ex_r_real*ex_r_real;
        for(size_t i=0; i<model_template.size(); ++i){
            collision_mask.push_back( MD_DEFS::MaskList{} );
            collision_mask.at(i).clear();
            for(size_t j=0; j<model_template.size(); ++j){
                PS::F64vec r_vec = model_template.at(j).getPos() - model_template.at(i).getPos();
                PS::F64    r2    = r_vec*r_vec;
                if(r2 <= ex_r2){
                    collision_mask.at(i).push_back( (MD_DEFS::IntraMask{}).setId( model_template.at(j).getId() ) );
                }
            }
        }
    }

    template <class Tpsys>
    PS::S64 check_excluded_vol(const PS::F64vec         pos_new,
                               const Tpsys             &psys,
                               const PS::F64            r_ex_norm,
                               const PS::S64            st,
                               const PS::S64            end,
                               const MD_DEFS::MaskList &mask_list ){

        PS::F64 r2_ex = r_ex_norm*r_ex_norm;
        for(PS::S64 i=st; i<end; ++i){
            //--- check mask
            if( MD_DEFS::isFind_mask(mask_list, psys[i].getAtomID()) ) continue;

            PS::F64vec r_vec = psys[i].getPos() - pos_new;
                       r_vec = Normalize::periodicAdjustNorm(r_vec);
            PS::F64    r2    = r_vec*r_vec;

            if(r2 <= r2_ex){
                return i;
            }
        }

        return end;
    }

    template <class Tpsys, class Tindex>
    typename std::vector<Tindex>::const_iterator check_excluded_vol(const PS::F64vec           pos_new,
                                                                    const Tpsys               &psys,
                                                                    const PS::F64              r_ex_norm,
                                                                    const std::vector<Tindex> &index_list,
                                                                    const MD_DEFS::MaskList   &mask_list  ){

        PS::F64 r2_ex = r_ex_norm*r_ex_norm;
        for(auto itr = index_list.begin(); itr != index_list.end(); ++itr){
            //--- check mask
            const Tindex index = *itr;
            if( MD_DEFS::isFind_mask(mask_list, psys[index].getAtomID()) ) continue;

            PS::F64vec r_vec = psys[index].getPos() - pos_new;
                       r_vec = Normalize::relativePosAdjustNorm(r_vec);
            PS::F64    r2    = r_vec*r_vec;

            if(r2 <= r2_ex){
                return itr;
            }
        }

        return index_list.end();
    }

    template <class Trand, class Tdist, class Tblz>
    PS::F64vec calc_vel(const PS::F64 &temperature,
                        const PS::F64 &mass,
                              Trand   &rand,
                              Tdist   &dist,
                              Tblz    &blz         ){

        //--- intensity
        PS::F64    dev = std::sqrt( 2.0*(temperature/Unit::norm_temp)/mass );
        PS::F64vec v   = 0.0;

        v.x = dev*blz.gen( PS::F64(dist(rand)) );

        //--- direction
        v = VEC_EXT::rot_z(v, Unit::pi*PS::F64(dist(rand)) );  // rotate in z-y-z Euler angle
        v = VEC_EXT::rot_y(v, Unit::pi*PS::F64(dist(rand)) );
        v = VEC_EXT::rot_z(v, Unit::pi*PS::F64(dist(rand)) );

        return v;
    }

    template <class Tpsys,
              class Tptcl,
              class Trand, class Tdist, class Tblz>
    void install_molecule(      Tpsys                          &psys,
                          const std::vector<Tptcl>             &model_template,
                          const PS::F64                         ex_r_norm,
                                MD_EXT::CellIndex<size_t>      &cell_index_inter,
                                MD_EXT::CellIndex<size_t>      &cell_index_intra,
                          const size_t                          try_limit,
                          const std::vector<MD_DEFS::MaskList> &collision_mask_template,
                          const PS::F64                         temperature,
                                Trand                          &rand,
                                Tdist                          &dist,
                                Tblz                           &blz,
                          const MD_DEFS::ID_type                atom_inclement,
                          const MD_DEFS::ID_type                mol_inclement){

        size_t try_count = 0;

        //--- set molecule properties
        for(size_t i=0; i<model_template.size(); ++i){
            MD_DEFS::ID_type atom_id = atom_inclement + i;
            psys[atom_id].copyFromModelTemplate(atom_inclement,
                                                mol_inclement,
                                                model_template.at(i) );
        }

        //--- make collision mask from template
        std::vector<MD_DEFS::MaskList> collision_mask;
        collision_mask.resize(model_template.size());
        for(size_t i=0; i<model_template.size(); ++i){
            collision_mask.at(i).clear();
            for(const auto& m : collision_mask_template.at(i)){
                MD_DEFS::IntraMask mask_tmp;
                collision_mask.at(i).push_back( mask_tmp.setId(m.id + atom_inclement) );
            }
        }

        //--- setting position
        for(; try_count<=try_limit; ++try_count){
            PS::F64vec pos_root = { dist(rand), dist(rand), dist(rand) };  // normalized
            PS::F64    rot_1    = Unit::pi*PS::F64(dist(rand));
            PS::F64    rot_2    = Unit::pi*PS::F64(dist(rand));
            PS::F64    rot_3    = Unit::pi*PS::F64(dist(rand));

            //--- set new position of atom (root pos + local_pos)
            cell_index_intra.clear();
            bool pos_flag = true;
            for(size_t i=0; i<model_template.size(); ++i){
                PS::F64vec pos_local = Normalize::normPos( model_template.at(i).getPos() );

                pos_local = VEC_EXT::rot_z(pos_local, rot_1);  // rotate in z-y-z Euler angle
                pos_local = VEC_EXT::rot_y(pos_local, rot_2);
                pos_local = VEC_EXT::rot_z(pos_local, rot_3);

                PS::F64vec pos_norm = pos_root + pos_local;
                pos_norm = periodicAdjustNorm_round(pos_norm);

                //--- check inter molecular collision
                auto index_list = cell_index_inter.get_data_list(pos_norm, ex_r_norm);
                auto itr = check_excluded_vol(pos_norm,
                                              psys,
                                              ex_r_norm,
                                              index_list,
                                              MD_DEFS::MaskList{} );
                if( itr != index_list.end() ){
                    pos_flag = false;
                    break;
                }

                //--- check intra molecular collision
                index_list = cell_index_intra.get_data_list(pos_norm, ex_r_norm);
                itr = check_excluded_vol(pos_norm,
                                         psys,
                                         ex_r_norm,
                                         index_list,
                                         collision_mask.at(i) );
                if( itr != index_list.end() ){
                    pos_flag = false;
                    break;
                }

                //--- add new pos into psys
                psys[atom_inclement + i].setPos(pos_norm);

                //--- record index and pos into cell index for intramolecular collision check
                cell_index_intra.add(pos_norm, atom_inclement + i);
            }

            if(try_count == try_limit){
                std::ostringstream oss;
                oss << "box size is too small or molecules are too much." << "\n"
                    << "   n_try = "       << try_count << "\n"
                    << "     n_atom_mol: " << model_template.size()
                    << " / MolID: "        << mol_inclement << "\n";

                throw std::invalid_argument(oss.str());
            }

            if( !pos_flag ) {
                continue;  // retry
            } else {
                break;     // position was set
            }
        }

        //--- set velocity
    //    PS::F64 mass = 0.0;
    //    for(size_t i=0; i<model_template.size(); ++i){
    //        mass += model_template.at(i).getMass();
    //    }
    //    PS::F64vec v = calc_vel(temperature, mass, rand, dist, blz);
    //    for(size_t i=0; i<model_template.size(); ++i){
    //        psys[atom_inclement + i].setVel(v);
    //    }
        for(size_t i=0; i<model_template.size(); ++i){
            psys[atom_inclement + i].setVel( calc_vel(temperature, psys[i].getMass(), rand, dist, blz) );
        }

        #ifdef TEST_MOL_INSTALL
            std::cout << "   n_try = "       << try_count << "\n"
                      << "     n_atom_mol: " << model_template.size()
                      << " / MolID: "        << mol_inclement << std::endl;
        #endif
    }


    //--- initialize PS::ParticleSystem<Tptcl>
    template <class Tpsys>
    void InitParticle(      Tpsys                                    &psys,
                      const std::vector<std::pair<MolName, PS::S64>> &model_list,
                      const std::vector<std::vector<Atom_FP>>        &model_template,
                      const PS::F64                                   ex_r_real,
                      const PS::S64                                   try_limit,
                      const PS::F64                                   temperature){

        if( PS::Comm::getRank() != 0 ) return;

        std::cout << " Initialize: temperature = " << temperature << " [K]" << std::endl;
        if(temperature < 0.0) throw std::invalid_argument("temperature must be > 0");

        //--- get array size
        PS::S64 n_total = 0;
        for(size_t index=0; index<model_list.size(); ++index){
            PS::S64 n_atom_mol = model_template.at(index).size();
            PS::F64 n_insert   = model_list.at(index).second;
            n_total += n_atom_mol*n_insert;

            #ifndef TEST_MOL_INSTALL
                if(n_insert == 0) continue;
            #endif

            std::cout << " Initialize: model_name = " << model_list.at(index).first << "\n"
                      << "               n_atom_mol = " << std::setw(10) << n_atom_mol
                      <<              " *n_insert = "   << std::setw(10) << n_insert   << std::endl;
        }

        //--- allocate arrays
        psys.setNumberOfParticleLocal(n_total);

        std::cout << " Initialize: total atoms = " << n_total << std::endl;

        //--- initialize random number & distribution generator
        constexpr int                    seed = std::pow(2, 19) + 1;
        std::mt19937_64                  mt(seed);
        std::uniform_real_distribution<> dist(0.0, 1.0);   // [0.0, 1.0)
        MD_EXT::boltzmann_dist           blz_dist;

        //--- initialize cell index neigbor list
        MD_EXT::CellIndex<size_t> cell_index_inter;
        MD_EXT::CellIndex<size_t> cell_index_intra;
        cell_index_inter.init( PS::F32vec{0.0, 0.0, 0.0},
                               PS::F32vec{1.0, 1.0, 1.0},
                               Normalize::normCutOff(ex_r_real) );
        cell_index_intra.init( PS::F32vec{0.0, 0.0, 0.0},
                               PS::F32vec{1.0, 1.0, 1.0},
                               Normalize::normCutOff(ex_r_real) );

        //--- input molecule
        MD_DEFS::ID_type atom_inclement = 0;
        MD_DEFS::ID_type mol_inclement  = 0;
        std::vector<MD_DEFS::MaskList> collision_mask;
        for(size_t index=0; index<model_list.size(); ++index){
            //--- make collision mask for the model
            make_collision_mask(collision_mask,
                                model_template.at(index),
                                ex_r_real);

            //--- insert a molecule in the system
            for(PS::S64 num=0; num<model_list.at(index).second; ++num){
                #ifdef TEST_MOL_INSTALL
                    std::cout << " installing " << ENUM::what(model_list.at(index).first) << std::endl;
                #endif

                install_molecule(psys,
                                 model_template.at(index),
                                 Normalize::normCutOff(ex_r_real), // ex_r
                                 cell_index_inter,
                                 cell_index_intra,
                                 try_limit,                        // try_limit
                                 collision_mask,                   // intramolecular collision mask
                                 temperature,                      // temperature
                                 mt, dist, blz_dist,
                                 atom_inclement,
                                 mol_inclement  );

                for(size_t i=0; i<model_template.at(index).size(); ++i){
                    const PS::S64 id = atom_inclement + i;
                    cell_index_inter.add(psys[id].getPos(), id);
                }

                atom_inclement += model_template.at(index).size();
                ++mol_inclement;
            }
        }
    }
}
