//***************************************************************************************
//  This is initializer module for md_fdps.
//***************************************************************************************

#include <cmath>
#include <iostream>
#include <vector>

#include <particle_simulator.hpp>

#include "ext_sys_control.hpp"
#include "md_coef_table.hpp"
#include "file_IO_pos.hpp"
#include "file_IO_resume.hpp"
#include "file_Out_VMD.hpp"
#include "atom_class.hpp"
#include "initialize.hpp"


int main(int argc, char *argv[]) {
    PS::Initialize(argc, argv);

    //--- display total threads for FDPS
    if(PS::Comm::getRank() == 0){
        std::cout << "Number of processes          : " << PS::Comm::getNumberOfProc()   << "\n"
                  << "Number of threads per process: " << PS::Comm::getNumberOfThread() << std::endl;
    }

    //--- make particle system object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;

    atom.initialize();
    atom.setNumberOfParticleLocal(0);

    //--- make ext_sys controller object
    EXT_SYS::Sequence   ext_sys_sequence;
    EXT_SYS::Controller ext_sys_controller;

    //--- load settings
    if(PS::Comm::getRank() == 0){
        //--- display normalizing units
        Unit::print_unit();

        System::loading_sequence_condition(MD_DEFS::condition_sequence_file,
                                           ext_sys_sequence,
                                           ext_sys_controller );
        System::loading_molecular_condition(MD_DEFS::condition_molecule_file);

        System::model_template.resize( System::model_list.size() );

        for(size_t i=0; i<System::model_list.size(); ++i){
            MODEL::loading_model_parameter(ENUM::what(System::model_list.at(i).first),
                                           System::model_template.at(i),
                                           MODEL::coef_table                          );
        }
        System::print_profile();
        ext_sys_controller.print();
        ext_sys_sequence.print();

        PS::F64 init_temperature = ext_sys_sequence.getSetting(0).temperature;
        Initialize::InitParticle(atom,
                                 System::model_list,
                                 System::model_template,
                                 System::get_ex_radius(),
                                 System::get_try_limit(),
                                 init_temperature);

        //--- warning
        if(System::get_istep() != 0){
            std::cout << "\n"
                      << "WARNING: the 't_start' is not 0." << "\n"
                      << "         this code will generate the resume file for t_start = 0." << "\n"
                      << "         check your setting file: '" << MD_DEFS::condition_sequence_file << "'" << "\n"
                      << std::endl;
        }
    }

    //--- send settings to all MPI processes
    System::broadcast_profile(0);
    MODEL::coef_table.broadcast(0);
    ext_sys_sequence.broadcast(0);
    ext_sys_controller.broadcast(0);

    //--- devide atom particle in MPI processes
    System::InitDinfo(dinfo);
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    std::cout << "proc = " << PS::Comm::getRank() << " / atoms = " << atom.getNumberOfParticleLocal() << std::endl;

    //--- output resume file
    System::profile.istep = 0;
    FILE_IO::VMDFileManager    vmd_file_mngr{    MD_DEFS::VMD_data_dir   , 0, 1};
    FILE_IO::PosFileManager    pos_file_mngr{    MD_DEFS::pos_data_dir   , 0, 1};
    FILE_IO::ResumeFileManager resume_file_mngr{ MD_DEFS::resume_data_dir, 0, 1};

    vmd_file_mngr.record(atom, System::profile);
    pos_file_mngr.record(atom, System::profile);
    resume_file_mngr.record(atom, System::profile, ext_sys_controller);

    PS::Finalize();
    return 0;
}
