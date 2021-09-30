//***************************************************************************************
//  This program is the main routine of "aTS-EVB" code with FDPS.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************

#include<iostream>
#include<fstream>
#include<cmath>
#include<unistd.h>
#include<sys/stat.h>

//--- FDPS headers
#include<particle_simulator.hpp>
#include<particle_mesh.hpp>

//--- user defined headers
#include "atsevb_enum_const.hpp"
#include "atsevb_normalize_pos.hpp"
// #include "atsevb_FDPS_f90interface.hpp"
#include "atsevb_FDPS_mask_list.hpp"
#include "atsevb_FDPS_force_inter.hpp"
#include "atsevb_FDPS_atom_class.hpp"

//--- main routine
int main(int argc, char* argv[]){
    std::cout<<std::setprecision(15);
    
    //--- initialize FDPS
    PS::Initialize(argc, argv);
    PS::F64 theta = 0.5;
    const PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 cycle_dinfo = 1;      // the time step period for refresh domain information.
    
            boxdh_x=40;
            boxdh_y=40;
            boxdh_z=40;
            setNormalizePosition();
    
    PS::S32 istep;
    PS::S32 nstep_st;
    PS::S32 nstep_ed;
    
    //--- cutoff radius
    EPI_LJ::RSearch_LJ = 12.0;                       // 12 angm. real space
    EPI_coulomb::RSearch_coulomb = 3.0/SIZE_OF_MESH; //  3 mesh. normalized space
    
    //--- display total threads for FDPS
    if(PS::Comm::getRank() == 0){
        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }
    
//    //--- setting the system in Fortran routine
//    if(PS::Comm::getRank() == 0){
//        pre_init_sys();
//    }
//    istep    = getISTEP();
//    nstep_ed = getNSTEP();
    
    //--- make the system data in FDPS
    PS::ParticleSystem<FP> system;
    system.initialize();
    
    //--- make the intramolecular mask list
    intraMaskList_initialize();
    
    //--- for test
    nstep_st = 0;
    istep    = nstep_st;
    nstep_ed = 20;
    PS::S32 n_tot;
    PS::S32 n_loc;
    
    //--- SetParticle
    if(PS::Comm::getRank() == 0){
        //--- for test particle, set in rank 0.
        n_loc = 2;
        //n_loc = 100;
        
        system.setNumberOfParticleLocal(n_loc);
        
        system[0].type    = solvent;
        system[0].atom_id =  0;
        system[0].mol_id  =  0;
		#if 0
        system[0].pos.x   =  10.0;
        system[0].pos.y   =  10.0;
        system[0].pos.z   =  10.0;
		#endif
        system[0].pos.x   =  0.1;
        system[0].pos.y   =  0.25;
        system[0].pos.z   =  0.75;
        system[0].charge  =  1.0;
        system[0].VDW_D   =  0.05;
        system[0].VDW_R   =  2.5;
        //normPos(system[0].pos);      // convert to normalized space
        
        system[1].type    = solvent;
        system[1].atom_id =  1;
        system[1].mol_id  =  1;
        //system[1].pos.x   =  10.5;
        //system[1].pos.y   =  10.0;
        //system[1].pos.z   =  10.0;
        system[1].pos.x   =  0.72;
        system[1].pos.y   =  0.3;
        system[1].pos.z   =  0.6;
        system[1].charge  = -1.0;
        system[1].VDW_D   =  0.05;
        system[1].VDW_R   =  2.5;
        //normPos(system[1].pos);      // convert to normalized space
        
      //  //--- dummy particle, set in Rank == 0.
      //  PS::MTTS mt;
      //  mt.init_genrand(0);
      //  for(PS::S32 i=2; i<n_loc; i++){
      //      system[i].type    = solvent;
      //      system[i].atom_id =  i;
      //      system[i].mol_id  =  i;
      //      
      //      //--- randam distribution for dummy particle
      //      system[i].pos.x   =  (2.0*mt.genrand_res53() - 1.0)*boxdh_x;
      //      system[i].pos.y   =  (2.0*mt.genrand_res53() - 1.0)*boxdh_y;
      //      system[i].pos.z   =  (2.0*mt.genrand_res53() - 1.0)*boxdh_z;
      //      
      //      system[i].charge  =  0.0;  // dummy: no interaction
      //      system[i].VDW_D   =  0.0;  // dummy: no interaction
      //      system[i].VDW_R   =  2.5;
      //  }
        
    } else {
      n_loc = 0;
        
      //  //--- dummy particle, set in rank >= 1.
      //  n_loc = 100;
      //  system.setNumberOfParticleLocal(n_loc);
      //  
      //  PS::MTTS mt;  // Mersenne Twister built in FDPS
      //  mt.init_genrand(PS::Comm::getRank());
      //  
      //  PS::S32 pid = PS::Comm::getRank();
      //  
      //  for(PS::S32 i=0; i<n_loc; i++){
      //      system[i].type    = solvent;
      //      system[i].atom_id =  i + n_loc*pid;
      //      system[i].mol_id  =  i + n_loc*pid;
      //      
      //      //--- randam distribution for dummy particle
      //      system[i].pos.x   =  (2.0*mt.genrand_res53() - 1.0)*boxdh_x;
      //      system[i].pos.y   =  (2.0*mt.genrand_res53() - 1.0)*boxdh_y;
      //      system[i].pos.z   =  (2.0*mt.genrand_res53() - 1.0)*boxdh_z;
      //      
      //      system[i].charge  =  0.0;  // dummy: no interaction
      //      system[i].VDW_D   =  0.0;  // dummy: no interaction
      //      system[i].VDW_R   =  2.5;
      //  }
    }
    
    //--- get total summation of n_loc
    n_tot = PS::Comm::getSum(n_loc);
    
    //--- make the domain information data
    const PS::F64 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    
    //------ set boundary condition
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    
    //------ set system box size
    dinfo.setPosRootDomain(PS::F64vec( 0.0, 0.0, 0.0), PS::F64vec( 1.0, 1.0, 1.0) );  // normalized space
    
    dinfo.collectSampleParticle(system);
    dinfo.decomposeDomain();
    
    system.exchangeParticle(dinfo);
    n_loc = system.getNumberOfParticleLocal();
    
    //--- make the tree structure
    //PS::TreeForForceShort<ForceScatter,
    //                      EPI_scatter,
    //                      EPJ_scatter>::Scatter tree_direct;
    //tree_direct.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    PS::TreeForForceShort<ForceLJ, EPI_LJ, EPJ_LJ>::Scatter tree_LJ;
    PS::TreeForForceShort<ForceCoulomb, EPI_coulomb, EPJ_coulomb>::Scatter tree_coulomb;
    tree_LJ.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    tree_coulomb.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    if(PS::Comm::getRank() == 0) std::cout << "TEST1" << std::endl;
    
    //--- make the Particle Mesh function object
    PS::PM::ParticleMesh pm;
    if(PS::Comm::getRank() == 0) std::cout << "TEST2" << std::endl;
    
    //--- calculate interaction
    //------ direct intermolecular force
    //tree_direct.calcForceAllAndWriteBack(CalcForceScatter<EPI_scatter, EPJ_scatter>, system, dinfo);
    tree_LJ.calcForceAllAndWriteBack(CalcForceLJ<EPI_LJ, EPJ_LJ>,
                                     system, dinfo);
    tree_coulomb.calcForceAllAndWriteBack(CalcForceCoulomb_perm<EPI_coulomb, EPJ_coulomb>,
                                          system, dinfo);
    //------ reciprocal coulomb force (PM)
    pm.calcForceAllAndWriteBack(system, dinfo);
    
    
    //--- thermodynamic treatment in Fortran routine
//    if(PS::Comm::getRank() == 0){
//        post_init_sys();
//    }
    
    std::cout << "main loop start!" << std::endl;
    for(istep = nstep_st; istep < nstep_ed; istep++){
        
        //--- clear force in FP
        for(PS::S32 i=0; i<n_loc; i++){
            system[i].clear();
        }
        
//        //--- pre-force treatment in Fortran routine
//        if(PS::Comm::getRank() == 0){
//          pre_force_treat(&istep);    
//        }
//        
        //--- test move
        for(PS::S32 i=0; i<n_loc; i++){
            if(system[i].atom_id == 1){
                realPos(system[i].pos);
                system[i].pos.x += 1.0;
                normPos(system[i].pos);
            }
        }
        
        //--- refresh domain info
        if( ((istep - nstep_st) % cycle_dinfo) == 0 ){
#ifdef DINFO_TEST
            std::cout << "Istep :" << istep << " refresh domain infomation." << std::endl;
#endif
            dinfo.decomposeDomainAll(system);
        }
        system.exchangeParticle(dinfo);
        n_loc = system.getNumberOfParticleLocal();
        
        //--- calculate intermolecular force in FDPS
        //tree_direct.calcForceAllAndWriteBack(CalcForceScatter<EPI_scatter, EPJ_scatter>, system, dinfo);
        tree_LJ.calcForceAllAndWriteBack(CalcForceLJ<EPI_LJ, EPJ_LJ>,
                                         system, dinfo);
        tree_coulomb.calcForceAllAndWriteBack(CalcForceCoulomb_perm<EPI_coulomb, EPJ_coulomb>,
                                              system, dinfo);
        pm.calcForceAllAndWriteBack(system, dinfo);
        
        //--- output test result
        PS::F64  init_v = 1000000.0;
        PS::F64  g_pos[n_tot];
        PS::F64  buff_pos[n_tot];
        for(PS::S32 i=0; i<n_tot; i++){
            g_pos[i] = 0.0;
            buff_pos[i] = init_v;
        }
        for(PS::S32 i=0; i<n_loc; i++){
            for(PS::S32 j=0; j<n_tot; j++){
                if(system[i].atom_id == j){
                    buff_pos[j] = system[i].pos.x;
                }
            }
        }
        for(PS::S32 i=0; i<n_tot; i++){
            if(buff_pos[i] != init_v){
                g_pos[i] = buff_pos[i];
            }
        }
        for(PS::S32 i=1; i<PS::Comm::getNumberOfProc(); i++){
            PS::Comm::broadcast(buff_pos, n_tot, i);
            for(PS::S32 i=0; i<n_tot; i++){
                if(buff_pos[i] != init_v){
                    g_pos[i] = buff_pos[i];
                }
            }
        }
        PS::Comm::broadcast(g_pos, n_tot, 0);
        
        
    //    if(PS::Comm::getRank() == 0){
    //        std::cout << "rPos : " << double(system[1].pos.x - system[0].pos.x) << std::endl;
    //    }
        for(PS::S32 i=0; i<n_loc; i++){
            if(system[i].atom_id == 0){
                std::cout << " Fx_coulomb : " << system[i].force_coulomb.x << std::endl;
                std::cout << "pot_coulomb : " << system[i].pot_coulomb     << std::endl;
                std::cout << " Fx_LJ      : " << system[i].force_LJ.x << std::endl;
                std::cout << "pot_LJ      : " << system[i].pot_LJ     << std::endl;
                std::cout << " Fx_PM      : " << system[i].force_PM.x << std::endl;
            }
        }
        
#ifdef VIRIAL_TEST
        checkVirial_inter(system, n_loc);
#endif
//
//        //--- post-force treatment in Fortran routine
//        if(PS::Comm::getRank() == 0){
//          post_force_treat(&istep);
//        }
    }
    std::cout << "main loop ends!" << std::endl;
    
    //--- finalize FDPS
    PS::Finalize();
    return 0;
}

