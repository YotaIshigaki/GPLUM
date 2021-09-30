/**************************************************************************************************/
/**
* @file  md_force.hpp
* @brief force calculater interface class.
*/
/**************************************************************************************************/
#pragma once

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>
#include <molecular_dynamics_ext.hpp>

#include "atom_class.hpp"
#include "md_coef_table.hpp"
#include "ff_intra_force.hpp"
#include "ff_inter_force.hpp"
#include "md_setting.hpp"


/**
* @brief force calculater interface class.
*/
class CalcForce {
private:
    //--- FDPS object
    PS::PM::ParticleMesh pm;
    PS::TreeForForceShort<ForceInter<PS::F32>, EP_inter, EP_inter>::Scatter tree_inter;
    PS::TreeForForceShort<ForceIntra<PS::F32>, EP_intra, EP_intra>::Scatter tree_intra;

    //--- intra pair list maker
    struct GetBond {
        MD_EXT::basic_connect<MD_DEFS::ID_type,
                              MD_DEFS::max_bond> operator () (const AtomConnect &atom){
            return atom.bond;
        }
    };
    IntraPair::IntraMaskMaker<  MD_DEFS::ID_type, GetBond> intra_mask_maker;
    IntraPair::AngleListMaker<  MD_DEFS::ID_type, GetBond> angle_list_maker;
    IntraPair::TorsionListMaker<MD_DEFS::ID_type, GetBond> torsion_list_maker;

public:
    void init(const PS::S64 &n_total){
        tree_inter.initialize(n_total,
                              System::profile.theta,
                              System::profile.n_leaf_limit,
                              System::profile.n_group_limit);
        tree_intra.initialize(n_total,
                              System::profile.theta,
                              System::profile.n_leaf_limit,
                              System::profile.n_group_limit);
    }

    /**
    * @brief update cutoff length in normalized space.
    */
    void setRcut(){
        EP_inter::setRcut_LJ(      Normalize::normCutOff( System::get_cut_off_LJ() ) );
        EP_inter::setRcut_coulomb( Normalize::normCutOff_PM() );

        EP_intra::setRcut( Normalize::normCutOff( System::get_cut_off_intra() ) );

        //--- check cut off length
        if(EP_inter::getRcut_LJ() >= 0.5 ||
           EP_inter::getRcut_LJ() <= 0.0 ){
            std::ostringstream oss;
            oss << "RSearch for LJ must be in range of (0.0, 0.5) at normalized space." << "\n"
                << "    EP_inter::getRcut_LJ() = " << EP_inter::getRcut_LJ() << "\n";
            throw std::length_error(oss.str());
        }
        if(EP_intra::getRSearch() >= 0.5 ||
           EP_intra::getRSearch() <= 0.0 ){
            std::ostringstream oss;
            oss << "RSearch for intramolecuar force must be in range of (0.0, 0.5) at normalized space." << "\n"
                << "    EP_intra::getRSearch() = " << EP_intra::getRSearch() << "\n";
            throw std::length_error(oss.str());
        }
    }

    /**
    * @brief update intra pair list at each atom.
    */
    template <class Tpsys, class Tdinfo, class Tmask>
    void update_intra_pair_list(      Tpsys  &atom,
                                      Tdinfo &dinfo,
                                const Tmask  &mask_table){
	PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) 
 	   std::cout << "1st call of update_intra_pair_list()" << std::endl;
        //--- get neighbor EP_intra information (do not calculate force)
        this->setRcut();
        this->tree_intra.calcForceAll( IntraPair::dummy_func{},
                                       atom,
                                       dinfo );

        PS::S32 n_local = atom.getNumberOfParticleLocal();
        for(PS::S32 i=0; i<n_local; ++i){
            atom[i].clear_intra_list();
            const auto& mask_param = mask_table.at(atom[i].getMolType());

            intra_mask_maker(  atom[i], this->tree_intra, mask_param, atom[i].mask_list());
            angle_list_maker(  atom[i], this->tree_intra, atom[i].angle_list());
            torsion_list_maker(atom[i], this->tree_intra, atom[i].dihedral_list(), atom[i].improper_list());
        }

    }

    /**
    * @brief update force on atom.
    */
    template <class Tpsys, class Tdinfo>
    void update_force(Tpsys  &atom,
                      Tdinfo &dinfo){

	PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) 
 	   std::cout << "1st call of update_force()" << std::endl;
        //--- clear force
        PS::S64 n_local = atom.getNumberOfParticleLocal();
        for(PS::S64 i=0; i<n_local; ++i){
            atom[i].clear();
        }

        this->setRcut();

        //=================
        //* PM part
        //=================
        pm.setDomainInfoParticleMesh(dinfo);
        pm.setParticleParticleMesh(atom, true);   // clear previous charge information
        pm.calcMeshForceOnly();

        //--- get potential and field
        for(PS::S64 i=0; i<n_local; ++i){
            PS::F32vec pos = atom[i].getPos();

            atom[i].addFieldCoulomb( Normalize::realPMForce(     -pm.getForce(pos)     ) );
            atom[i].addPotCoulomb(   Normalize::realPMPotential( -pm.getPotential(pos) ) );
        }

        //=================
        //* PP part
        //=================
        //--- calculate force
        this->tree_inter.calcForceAll(FORCE::calcForceShort<ForceInter<PS::F32>, EP_inter, EP_inter>,
                                      atom,
                                      dinfo);
        for(PS::S64 i=0; i<n_local; ++i){
            const auto& result = tree_inter.getForce(i);
            atom[i].addFieldCoulomb( result.getFieldCoulomb() );
            atom[i].addPotCoulomb(   result.getPotCoulomb() );
            atom[i].addForceLJ(  result.getForceLJ()  );
            atom[i].addPotLJ(    result.getPotLJ()    );
            atom[i].addVirialLJ( result.getVirialLJ() );
        }

        //=================
        //* Intra force part
        //=================
        FORCE::calcForceIntra<ForceIntra<PS::F32>>(this->tree_intra, atom, dinfo);
    }
};
