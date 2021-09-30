//***************************************************************************************
//  This is high level function for intramolecular force calculation.
//***************************************************************************************
#pragma once

#include <tuple>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "md_coef_table.hpp"
#include "ff_intra_force_func.hpp"


namespace FORCE {

    template<class Tepi, class Ttree,
             class Tforce>
    void calcForceIntra_IA(const Tepi   &ep_i,
                                 Ttree  &tree,
                                 Tforce &force_IA){

        //--- bond potential
        for(const auto& id_j : ep_i.bond){
            const auto* ptr_j = tree.getEpjFromId(id_j);
            IntraPair::check_nullptr_ptcl(ptr_j, ep_i.getAtomID(), id_j);

            MODEL::KeyBond  key_bond = std::make_tuple(ep_i.getMolType(),
                                                       ep_i.getAtomType(),
                                                       (*ptr_j).getAtomType());

            if(MODEL::coef_table.bond.count(key_bond) != 1){
                std::ostringstream oss;
                oss << "key_bond = " << ENUM::what(key_bond) << " is not defined in MODEL::coef_table.bond." << "\n";
                throw std::invalid_argument(oss.str());
            }

            const auto bond_prm = MODEL::coef_table.bond.at(key_bond);
            switch (bond_prm.form) {
                case IntraFuncForm::none:
                    //--- free potential
                    continue;
                break;

                case IntraFuncForm::anharmonic:
                    calcBondForce_anharmonic_IJ(Normalize::realPos(ep_i.getPos()),
                                                Normalize::realPos((*ptr_j).getPos()),
                                                bond_prm,
                                                force_IA);
                break;

                case IntraFuncForm::harmonic:
                    calcBondForce_harmonic_IJ(Normalize::realPos(ep_i.getPos()),
                                              Normalize::realPos((*ptr_j).getPos()),
                                              bond_prm,
                                              force_IA);
                break;

                default:
                    throw std::invalid_argument("undefined function form: " + ENUM::what(bond_prm.form) + " at bond force.");
            }
        }

        //--- angle potential
        for(const auto& angle_set : ep_i.angle_list()){
            const auto id_i   = std::get<0>(angle_set);
            const auto id_j   = std::get<1>(angle_set);
            const auto id_k   = std::get<2>(angle_set);
            const auto id_tgt = ep_i.getAtomID();

            const auto* ptr_i = tree.getEpjFromId(id_i);
            const auto* ptr_j = tree.getEpjFromId(id_j);
            const auto* ptr_k = tree.getEpjFromId(id_k);

            IntraPair::check_nullptr_ptcl(ptr_i, ep_i.getAtomID(), id_i);
            IntraPair::check_nullptr_ptcl(ptr_j, ep_i.getAtomID(), id_j);
            IntraPair::check_nullptr_ptcl(ptr_k, ep_i.getAtomID(), id_k);

            MODEL::KeyAngle  key_angle = std::make_tuple(ep_i.getMolType(),
                                                         (*ptr_i).getAtomType(),
                                                         (*ptr_j).getAtomType(),
                                                         (*ptr_k).getAtomType() );

            if(MODEL::coef_table.angle.count(key_angle) != 1){
                std::ostringstream oss;
                oss << "key_angle = " << ENUM::what(key_angle) << " is not defined in MODEL::coef_table.angle." << "\n";
                throw std::invalid_argument(oss.str());
            }

            const auto angle_prm = MODEL::coef_table.angle.at(key_angle);
            switch (angle_prm.form){
                case IntraFuncForm::none:
                    //--- free potential
                    continue;
                break;

                case IntraFuncForm::harmonic:
                    calcAngleForce_harmonic_IJK(Normalize::realPos( (*ptr_i).getPos() ),
                                                Normalize::realPos( (*ptr_j).getPos() ),
                                                Normalize::realPos( (*ptr_k).getPos() ),
                                                id_i, id_j, id_k,
                                                id_tgt, angle_prm,
                                                force_IA);
                break;

                default:
                    throw std::invalid_argument("undefined function form: " + ENUM::what(angle_prm.form) + " at angle force");
            }
        }

        //--- dihedral torsion potential
        for(const auto& dihedral_set : ep_i.dihedral_list()){
            const auto id_i   = std::get<0>(dihedral_set);
            const auto id_j   = std::get<1>(dihedral_set);
            const auto id_k   = std::get<2>(dihedral_set);
            const auto id_l   = std::get<3>(dihedral_set);
            const auto id_tgt = ep_i.getAtomID();

            const auto* ptr_i = tree.getEpjFromId(id_i);
            const auto* ptr_j = tree.getEpjFromId(id_j);
            const auto* ptr_k = tree.getEpjFromId(id_k);
            const auto* ptr_l = tree.getEpjFromId(id_l);

            IntraPair::check_nullptr_ptcl(ptr_i, ep_i.getAtomID(), id_i);
            IntraPair::check_nullptr_ptcl(ptr_j, ep_i.getAtomID(), id_j);
            IntraPair::check_nullptr_ptcl(ptr_k, ep_i.getAtomID(), id_k);
            IntraPair::check_nullptr_ptcl(ptr_l, ep_i.getAtomID(), id_l);

            MODEL::KeyTorsion  key_torsion = std::make_tuple(ep_i.getMolType(),
                                                             TorsionShape::dihedral,
                                                             (*ptr_i).getAtomType(),
                                                             (*ptr_j).getAtomType(),
                                                             (*ptr_k).getAtomType(),
                                                             (*ptr_l).getAtomType() );

            if(MODEL::coef_table.torsion.count(key_torsion) != 1){
                std::ostringstream oss;
                oss << "key_torsion = " << ENUM::what(key_torsion) << " is not defined in MODEL::coef_table.torsion." << "\n";
                throw std::invalid_argument(oss.str());
            }

            const auto torsion_prm = MODEL::coef_table.torsion.at(key_torsion);
            switch (torsion_prm.form){
                case IntraFuncForm::none:
                    //--- free potential
                    continue;
                break;

                case IntraFuncForm::cos:
                    calcTorsionForce_harmonic_IJKL(Normalize::realPos( (*ptr_i).getPos() ),
                                                   Normalize::realPos( (*ptr_j).getPos() ),
                                                   Normalize::realPos( (*ptr_k).getPos() ),
                                                   Normalize::realPos( (*ptr_l).getPos() ),
                                                   id_i, id_j, id_k, id_l,
                                                   id_tgt, torsion_prm,
                                                   force_IA);
                break;

                case IntraFuncForm::OPLS_3:
                    calcTorsionForce_OPLS_3rd_IJKL(Normalize::realPos( (*ptr_i).getPos() ),
                                                   Normalize::realPos( (*ptr_j).getPos() ),
                                                   Normalize::realPos( (*ptr_k).getPos() ),
                                                   Normalize::realPos( (*ptr_l).getPos() ),
                                                   id_i, id_j, id_k, id_l,
                                                   id_tgt, torsion_prm,
                                                   force_IA);
                break;

                default:
                    throw std::invalid_argument("undefined potential type: " + ENUM::what(torsion_prm.form) + "at dihedral torsion force");
            }
        }

        //--- improper torsion potential
        for(const auto& improper_set : ep_i.improper_list()){
            const auto id_i   = std::get<0>(improper_set);
            const auto id_j   = std::get<1>(improper_set);
            const auto id_k   = std::get<2>(improper_set);
            const auto id_l   = std::get<3>(improper_set);
            const auto id_tgt = ep_i.getAtomID();

            const auto* ptr_i = tree.getEpjFromId(id_i);
            const auto* ptr_j = tree.getEpjFromId(id_j);
            const auto* ptr_k = tree.getEpjFromId(id_k);
            const auto* ptr_l = tree.getEpjFromId(id_l);

            IntraPair::check_nullptr_ptcl(ptr_i, ep_i.getAtomID(), id_i);
            IntraPair::check_nullptr_ptcl(ptr_j, ep_i.getAtomID(), id_j);
            IntraPair::check_nullptr_ptcl(ptr_k, ep_i.getAtomID(), id_k);
            IntraPair::check_nullptr_ptcl(ptr_l, ep_i.getAtomID(), id_l);

            MODEL::KeyTorsion  key_torsion = std::make_tuple(ep_i.getMolType(),
                                                             TorsionShape::improper,
                                                             (*ptr_i).getAtomType(),
                                                             (*ptr_j).getAtomType(),
                                                             (*ptr_k).getAtomType(),
                                                             (*ptr_l).getAtomType() );

            if(MODEL::coef_table.torsion.count(key_torsion) != 1){
                std::ostringstream oss;
                oss << "key_torsion = " << ENUM::what(key_torsion) << " is not defined in MODEL::coef_table.torsion." << "\n";
                throw std::invalid_argument(oss.str());
            }

            const auto torsion_prm = MODEL::coef_table.torsion.at(key_torsion);
            switch (torsion_prm.form){
                case IntraFuncForm::none:
                    //--- free potential
                    continue;
                break;

                case IntraFuncForm::cos:
                    calcTorsionForce_harmonic_IJKL(Normalize::realPos( (*ptr_i).getPos() ),
                                                   Normalize::realPos( (*ptr_j).getPos() ),
                                                   Normalize::realPos( (*ptr_k).getPos() ),
                                                   Normalize::realPos( (*ptr_l).getPos() ),
                                                   id_i, id_j, id_k, id_l,
                                                   id_tgt, torsion_prm,
                                                   force_IA);
                break;

                case IntraFuncForm::OPLS_3:
                    calcTorsionForce_OPLS_3rd_IJKL(Normalize::realPos( (*ptr_i).getPos() ),
                                                   Normalize::realPos( (*ptr_j).getPos() ),
                                                   Normalize::realPos( (*ptr_k).getPos() ),
                                                   Normalize::realPos( (*ptr_l).getPos() ),
                                                   id_i, id_j, id_k, id_l,
                                                   id_tgt, torsion_prm,
                                                   force_IA);
                break;

                default:
                    throw std::invalid_argument("undefined potential type: " + ENUM::what(torsion_prm.form) + "at improper torsion force");
            }
        }
    }


    //--- template function for FDPS interface
    template <class Tforce,
              class Ttree, class Tpsys, class Tdinfo>
    void calcForceIntra(Ttree  &tree,
                        Tpsys  &atom,
                        Tdinfo &dinfo){

        //--- get neighbor EP_intra information (do not calculate force)
        tree.calcForceAll( IntraPair::dummy_func{}, atom, dinfo);

        const PS::S64 n_local = atom.getNumberOfParticleLocal();

        //--- calculate intramolecular force
        for(PS::S64 i=0; i<n_local; ++i){

            //--- skip no-bonded atoms
            if( atom[i].bond.size() == 0 ) continue;

            Tforce force_IA;
                   force_IA.clear();
            calcForceIntra_IA(atom[i],
                              tree,
                              force_IA);
            atom[i].copyForceIntra(force_IA);
        }
    }

}
