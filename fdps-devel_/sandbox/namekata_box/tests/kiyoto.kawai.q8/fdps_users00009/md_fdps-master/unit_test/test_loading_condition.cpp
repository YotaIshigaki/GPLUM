//***************************************************************************************
//  This is unit test of loading system conditions.
//***************************************************************************************

#include <iostream>

#include <particle_simulator.hpp>

#include "md_defs.hpp"
#include "ext_sys_control.hpp"
#include "md_loading_condition.hpp"

int main(int argc, char *argv[]) {
    PS::Initialize(argc, argv);

    PS::DomainInfo dinfo;

    //--- make ext_sys controller object
    EXT_SYS::Sequence   ext_sys_sequence;
    EXT_SYS::Controller ext_sys_controller;

    if(PS::Comm::getRank() == 0) {
        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());

        //--- initialize
        std::cout << std::endl;
        std::cout << "--- settings loaded in process: " << PS::Comm::getRank() << std::endl;
        std::cout << std::endl;

        System::loading_sequence_condition(MD_DEFS::condition_sequence_file,
                                           ext_sys_sequence,
                                           ext_sys_controller );
        System::loading_molecular_condition(MD_DEFS::condition_molecule_file);
    }


    System::broadcast_profile();
    ext_sys_sequence.broadcast(0);
    ext_sys_controller.broadcast(0);
    std::cout << "--- sync in MPI broadcast ---" << std::endl;

    System::InitDinfo(dinfo);

    //--- display settings
    if(PS::Comm::getRank() == 1) {
        std::cout << std::endl;
        System::print_profile();
        ext_sys_controller.print();
        ext_sys_sequence.print();
    }


    PS::Finalize();
    return 0;
}
