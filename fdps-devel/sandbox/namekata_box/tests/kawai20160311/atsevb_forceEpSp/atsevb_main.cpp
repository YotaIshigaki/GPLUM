//***************************************************************************************
//  This program is the main routine of "aTS-EVB" code with FDPS.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************

//--- general headers
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<algorithm>
#include<unistd.h>
#include<sys/stat.h>

//--- FDPS headers
#include<particle_simulator.hpp>
#include<particle_mesh.hpp>

//--- external library for MD
#include<molecular_dynamics_ext.hpp>

//--- user defined headers
//------ definition of data set
#include "atsevb_unit.hpp"
#include "atsevb_enum_const.hpp"
#include "atsevb_FDPS_atom_class.hpp"
#include "atsevb_FDPS_support_data_class.hpp"
//------ setting data value and intramolecular interaction function
#include "atsevb_FDPS_set_water.hpp"
//#include "atsevb_FDPS_set_EVB.hpp"
//#include "atsevb_FDPS_set_naf.hpp"
//#include "atsevb_FDPS_set_metal_ion.hpp"
//------ calculate interaction
#include "atsevb_FDPS_intra_force.hpp"
#include "atsevb_FDPS_inter_mol_force.hpp"
//------ kick & drift
#include "atsevb_atom_move.hpp"
//------ external system control
#include "atsevb_FDPS_ext_sys_control.hpp"
//------ file I/O
#include "atsevb_FDPS_fileIO.hpp"
//------ initialize
#include "atsevb_initialize_system.hpp"


//--- main routine
int main(int argc, char* argv[]){
    //std::cout << std::fixed << std::setprecision(7);
    std::cout << std::scientific << std::setprecision(7);
    //std::cout << std::setprecision(15);
    
    //--- initialize FDPS
    PS::Initialize(argc, argv);
	
    //--- display total threads for FDPS
    if(PS::Comm::getRank() == 0){
        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }
	
    //--- call system condition object
    SystemSetting *setting = SystemSetting::getInstance();
    
    //--- make the domain information data
    PS::DomainInfo dinfo;
    //------ setting system condition
    Init_Condition(dinfo);
    
    
    //--- make the particle data in FDPS
    PS::ParticleSystem<FP> system;
    system.initialize();
    //------ install molecules
    Init_Molecule(system);
    
    
    //--- exchange particle
    dinfo.collectSampleParticle(system);
    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);
    PS::S32 n_loc = system.getNumberOfParticleLocal();
    PS::S32 n_tot = PS::Comm::getSum(n_loc);
    
    //--- make the tree structure
  //  PS::TreeForForceShort<ForceIntra, EPI_intra, EPJ_intra>::Scatter tree_intra;
  //  PS::TreeForForceShort<ForceLJ, EPI_LJ, EPJ_LJ>::Scatter tree_LJ;
  //  PS::TreeForForceShort<ForceCoulomb, EPI_coulomb, EPJ_coulomb>::Scatter tree_coulomb;
    PS::TreeForForceShort<ForceIntra, EPI_intra, EPJ_intra>::Gather tree_intra;
    PS::TreeForForceShort<ForceLJ, EPI_LJ, EPJ_LJ>::Gather tree_LJ;
  //  PS::TreeForForceShort<ForceCoulomb, EPI_coulomb, EPJ_coulomb>::Gather tree_coulomb;
    PS::TreeForForceLong<ForceCoulomb, EPI_coulomb, EPJ_coulomb>::MonopoleWithCutoff tree_coulomb;
    tree_intra.initialize(n_tot,
                          setting->getTheta(),
                          setting->getNLeafLimit(),
                          setting->getNGroupLimit());
    tree_LJ.initialize(n_tot,
                       setting->getTheta(),
                       setting->getNLeafLimit(),
                       setting->getNGroupLimit());
    tree_coulomb.initialize(n_tot,
                            setting->getTheta(),
                            setting->getNLeafLimit(),
                            setting->getNGroupLimit());
    
    //--- make the Particle Mesh object
    PS::PM::ParticleMesh pm;
    
    //--- calculate interaction
    //------ direct intermolecular force
  //  tree_intra.calcForceAllAndWriteBack(CalcForceIntra<EPI_intra, EPJ_intra>,
  //                                      system, dinfo);
    tree_LJ.calcForceAllAndWriteBack(CalcForceLJ<EPI_LJ, EPJ_LJ>,
                                     system, dinfo);
  //  tree_coulomb.calcForceAllAndWriteBack(CalcForceCoulomb_perm<EPI_coulomb, EPJ_coulomb>,
  //                                        system, dinfo);
    tree_coulomb.calcForceAllAndWriteBack(CalcForceCoulomb_perm<EPI_coulomb, EPJ_coulomb>,
                                          CalcForceCoulomb_EpSp<EPI_coulomb>,
                                          system, dinfo);
    //------ reciprocal coulomb force (PM)
    pm.calcForceAllAndWriteBack(system, dinfo);
    
    //--- record reference value of energy
    calcEnergy(system);
    EXT_SYS_CONTROL::EnergyBuf::getInstance()->setRef();  // set reference energy value
    
    //--- main loop
    if(PS::Comm::getRank() == 0) std::cout << "main loop start!" << std::endl;
    //for(istep = nstep_st; istep <= nstep_ed; istep++){
    while( setting->isLoopContinue() ){
        
        //--- output record
        FILE_IO::recordPos(system, setting->getIstep() );
        FILE_IO::recordVMD(system, setting->getIstep() );
        FILE_IO::recordResume(system, setting->getIstep() );
        
        //--- energy log
        calcEnergy(system);
        FILE_IO::recordEnergy(setting->getIstep());
        
        //--- debug pos trace
      //  FILE_IO::recordTrajecP0(system, setting.getIstep() );
        
        //--- kick
        Kick(system, setting->get_dt()*0.5 );
        
        //--- drift
        Drift(system, setting->get_dt() );
        
        //--- refresh domain info
        if( setting->isDinfoUpdate() ){
            #ifdef DINFO_TEST
                if(PS::Comm::getRank() == 0){
                    std::cerr << "Istep :" << setting->getIstep()
                              << " refresh domain infomation." << std::endl;
                }
            #endif
            dinfo.decomposeDomainAll(system);
        }
        //--- exchange particle
        system.exchangeParticle(dinfo);
        n_loc = system.getNumberOfParticleLocal();
        
        
        //--- calculate intermolecular force in FDPS
        //------ clear force in FP
        for(PS::S32 i=0; i<n_loc; i++){
            system[i].clear();
        }
    //    tree_intra.calcForceAllAndWriteBack(CalcForceIntra<EPI_intra, EPJ_intra>,
    //                                        system, dinfo);
        tree_LJ.calcForceAllAndWriteBack(CalcForceLJ<EPI_LJ, EPJ_LJ>,
                                         system, dinfo);
    //    tree_coulomb.calcForceAllAndWriteBack(CalcForceCoulomb_perm<EPI_coulomb, EPJ_coulomb>,
    //                                          system, dinfo);
        tree_coulomb.calcForceAllAndWriteBack(CalcForceCoulomb_perm<EPI_coulomb, EPJ_coulomb>,
                                              CalcForceCoulomb_EpSp<EPI_coulomb>,
                                              system, dinfo);
        pm.calcForceAllAndWriteBack(system, dinfo);
        
        //--- kick
        Kick(system, setting->get_dt()*0.5 );
        
        //--- nest step
        setting->StepNext();
    }
    if(PS::Comm::getRank() == 0) std::cout << "main loop ends!" << std::endl;
    
    //--- finalize FDPS
    PS::Finalize();
    return 0;
}

