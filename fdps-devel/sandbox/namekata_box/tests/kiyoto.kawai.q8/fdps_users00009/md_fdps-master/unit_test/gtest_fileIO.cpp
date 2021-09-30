//***************************************************************************************
//  This program is unit test of LJ interaction.
//***************************************************************************************

#include <gtest/gtest.h>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

//--- external library for MD
#include <molecular_dynamics_ext.hpp>

//--- user defined headers
#include "unit.hpp"
#include "md_enum.hpp"
#include "md_defs.hpp"
#include "md_coef_table.hpp"
#include "file_IO_pos.hpp"
#include "file_IO_resume.hpp"
#include "file_Out_VMD.hpp"
#include "atom_class.hpp"
#include "ext_sys_control.hpp"
#include "md_force.hpp"
#include "observer.hpp"
#include "initialize.hpp"


namespace TEST_DEFS {
    const std::string test_sequence_file = "./unit_test/ref/test_sequence.inp";
    const std::string test_molecule_file = "./unit_test/ref/test_molecule.inp";

    const std::string test_dir = "./test_bin";
};


std::string float_to_str(const PS::F32 f){
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(10) << f;
    return oss.str();
}

std::string float_to_str(const PS::F64 f){
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(18) << f;
    return oss.str();
}

template <class Tf>
::testing::AssertionResult float_complete_eq(const Tf lhs, const Tf rhs){
    if(lhs == rhs){
        return ::testing::AssertionSuccess();
    } else {
        std::ostringstream oss;
        oss << "\n"
            << "  lhs = " << float_to_str(lhs) << "\n"
            << "  rhs = " << float_to_str(rhs) << "\n";
        return ::testing::AssertionFailure() << oss.str();
    }
}

template <class Tpsys, class Tatom>
void allGatherAtom(Tpsys              &psys,
                   std::vector<Tatom> &vec_atom){

    std::vector<Tatom>              vec_local;
    std::vector<std::vector<Tatom>> vec_recv;

    const PS::S64 n_total = psys.getNumberOfParticleGlobal();

    PS::S64 n_local = psys.getNumberOfParticleLocal();
    vec_local.clear();
    for(PS::S32 i=0; i<n_local; ++i){
        vec_local.push_back( psys[i] );
    }

    COMM_TOOL::allGather(vec_local, vec_recv);
    vec_atom.clear();
    for(const auto& vec : vec_recv){
        for(const auto& atom : vec){
            vec_atom.push_back(atom);
        }
    }

    EXPECT_EQ(vec_atom.size(), n_total);
}

void SetUp_atom_data(PS::DomainInfo              &dinfo,
                     PS::ParticleSystem<Atom_FP> &atom,
                     CalcForce                   &force,
                     EXT_SYS::Sequence           &ext_sys_sequence,
                     EXT_SYS::Controller         &ext_sys_controller){

    //--- temporary data
    Observer::Energy eng;

    atom.initialize();
    atom.setNumberOfParticleLocal(0);

    System::model_list.clear();
    System::model_template.clear();

    MODEL::coef_table.clear();

    //--- load settings / model parameters
    if(PS::Comm::getRank() == 0){
        System::loading_sequence_condition(TEST_DEFS::test_sequence_file,
                                           ext_sys_sequence,
                                           ext_sys_controller);
        System::loading_molecular_condition(TEST_DEFS::test_molecule_file);

        System::model_template.resize( System::model_list.size() );

        for(size_t i=0; i<System::model_list.size(); ++i){
            MODEL::loading_model_parameter(ENUM::what(System::model_list.at(i).first),
                                           System::model_template.at(i),
                                           MODEL::coef_table                          );
        }

        PS::F64 init_temperature = ext_sys_sequence.getSetting(0).temperature;
        Initialize::InitParticle(atom,
                                 System::model_list,
                                 System::model_template,
                                 System::get_ex_radius(),
                                 System::get_try_limit(),
                                 init_temperature);
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

    //--- initialize force culculator
    PS::S64 n_total = atom.getNumberOfParticleGlobal();
    force.init(n_total);

    force.update_intra_pair_list(atom, dinfo, MODEL::coef_table.mask_scaling);
    force.update_force(atom, dinfo);

    //--- apply force integrator
    while( System::isLoopContinue() ){
        eng.getEnergy(atom);

        PS::S64 n_rigid = 0;
        ext_sys_controller.apply(n_rigid,
                                 ext_sys_sequence.getSetting( System::get_istep() ),
                                 System::get_dt(),
                                 atom,
                                 eng);

        ext_sys_controller.kick(0.5*System::get_dt(), atom);
        ext_sys_controller.drift(System::get_dt(), atom);

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);
        force.update_intra_pair_list(atom, dinfo, MODEL::coef_table.mask_scaling);

        //--- calculate intermolecular force in FDPS
        force.update_force(atom, dinfo);

        ext_sys_controller.kick(0.5*System::get_dt(), atom);

        //--- nest step
        System::StepNext();
    }
}


TEST(TestFileIO, pos){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    EXT_SYS::Sequence   ext_sys_sequence;
    EXT_SYS::Controller ext_sys_controller;

    PS::F64vec box_size;

    SetUp_atom_data(dinfo, atom, force, ext_sys_sequence, ext_sys_controller);

    //--- file output
    if(PS::Comm::getRank() == 0) std::cout << "  test directory: " << TEST_DEFS::test_dir << std::endl;
    System::profile.istep        = -1;
    System::profile.pos_start    = -1;
    System::profile.pos_interval =  1;

    //--- initialize by init() function
    FILE_IO::PosFileManager pos_file_mngr;
    pos_file_mngr.init( TEST_DEFS::test_dir,
                        System::profile.get_pos_start(),
                        System::profile.get_pos_interval() );
    const PS::F64 time_stamp = System::profile.get_trj_time();  //--- trj time will be cleared in PosFileManager::record()
    pos_file_mngr.record(atom, System::profile);

    box_size = Normalize::getBoxSize();

    //--- data from file
    System::Profile profile_file;
    PS::ParticleSystem<Atom_FP> atom_file;

    atom_file.initialize();
    atom_file.setNumberOfParticleLocal(0);

    //--- input from file
    profile_file.istep = -1;
    pos_file_mngr.load(atom_file, profile_file);

    std::cout << "  Proc = " << PS::Comm::getRank()
              << ", n_atom = " << atom_file.getNumberOfParticleLocal()
              << " / " << atom_file.getNumberOfParticleGlobal() << std::endl;

    //--- check header info
    EXPECT_EQ(System::get_istep(), profile_file.istep);
    EXPECT_TRUE( float_complete_eq(System::profile.get_time(), profile_file.get_time()  ) );
    EXPECT_TRUE( float_complete_eq(time_stamp,                 profile_file.get_trj_time()  ) );
    EXPECT_TRUE( float_complete_eq(box_size.x,                 Normalize::getBoxSize().x) );
    EXPECT_TRUE( float_complete_eq(box_size.y,                 Normalize::getBoxSize().y) );
    EXPECT_TRUE( float_complete_eq(box_size.z,                 Normalize::getBoxSize().z) );

    //--- check atom info
    std::vector<Atom_FP> v_atom;
    std::vector<Atom_FP> v_atom_file;

    allGatherAtom(atom,      v_atom);
    allGatherAtom(atom_file, v_atom_file);

    if(PS::Comm::getRank() == 0){
        std::sort(v_atom.begin(),
                  v_atom.end(),
                  [](const Atom_FP &a, const Atom_FP &b){ return ( a.getAtomID() < b.getAtomID() ); });
        std::sort(v_atom_file.begin(),
                  v_atom_file.end(),
                  [](const Atom_FP &a, const Atom_FP &b){ return ( a.getAtomID() < b.getAtomID() ); });

        ASSERT_EQ(v_atom_file.size(), v_atom.size());
        for(size_t i=0; i<v_atom.size(); ++i){
            EXPECT_EQ(v_atom_file[i].getAtomID()  , v_atom[i].getAtomID()  ) << "i = " << i;
            EXPECT_EQ(v_atom_file[i].getMolID()   , v_atom[i].getMolID()   ) << "i = " << i;
            EXPECT_EQ(v_atom_file[i].getAtomType(), v_atom[i].getAtomType()) << "i = " << i;
            EXPECT_EQ(v_atom_file[i].getMolType() , v_atom[i].getMolType() ) << "i = " << i;
            EXPECT_TRUE( float_complete_eq(v_atom_file[i].getPos().x, v_atom[i].getPos().x) ) << "i = " << i;
            EXPECT_TRUE( float_complete_eq(v_atom_file[i].getPos().y, v_atom[i].getPos().y) ) << "i = " << i;
            EXPECT_TRUE( float_complete_eq(v_atom_file[i].getPos().z, v_atom[i].getPos().z) ) << "i = " << i;
        }
    }
}

TEST(TestFileIO, resume){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    EXT_SYS::Sequence   ext_sys_sequence;
    EXT_SYS::Controller ext_sys_controller;

    PS::F64vec box_size;

    SetUp_atom_data(dinfo, atom, force, ext_sys_sequence, ext_sys_controller);

    //--- file output
    if(PS::Comm::getRank() == 0) std::cout << "  test directory: " << TEST_DEFS::test_dir << std::endl;
    System::profile.istep           = -1;
    System::profile.resume_start    = -1;
    System::profile.resume_interval =  1;

    //--- initialize by constructor
    FILE_IO::ResumeFileManager resume_file_mngr{ TEST_DEFS::test_dir,
                                                 System::profile.get_resume_start(),
                                                 System::profile.get_resume_interval() };
    resume_file_mngr.record(atom, System::profile, ext_sys_controller);

    box_size = Normalize::getBoxSize();

    //--- data from file
    System::Profile   profile_file;
    PS::ParticleSystem<Atom_FP> atom_file;
    EXT_SYS::Controller         controller_file;

    atom_file.initialize();
    atom_file.setNumberOfParticleLocal(0);

    //--- input from file
    profile_file.istep = -1;
    resume_file_mngr.load(atom_file, profile_file, controller_file);

    std::cout << "  Proc = " << PS::Comm::getRank()
              << ", n_atom = " << atom_file.getNumberOfParticleLocal()
              << " / " << atom_file.getNumberOfParticleGlobal() << std::endl;

    //--- check header info
    EXPECT_EQ(System::get_istep(), profile_file.istep);
    EXPECT_TRUE( float_complete_eq(System::profile.get_time()    , profile_file.get_time()  ) );
    EXPECT_TRUE( float_complete_eq(System::profile.get_trj_time(), profile_file.get_trj_time()  ) );
    EXPECT_TRUE( float_complete_eq(box_size.x,                     Normalize::getBoxSize().x) );
    EXPECT_TRUE( float_complete_eq(box_size.y,                     Normalize::getBoxSize().y) );
    EXPECT_TRUE( float_complete_eq(box_size.z,                     Normalize::getBoxSize().z) );

    //------ ext_sys info
    const auto stat_ref  = ext_sys_controller.get_resume();
    const auto stat_file = controller_file.get_resume();
    EXPECT_EQ(stat_file.n_chain , stat_ref.n_chain );
    EXPECT_EQ(stat_file.n_rep   , stat_ref.n_rep   );
    EXPECT_EQ(stat_file.n_nys   , stat_ref.n_nys   );
    EXPECT_TRUE( float_complete_eq(stat_file.NVT_freq, stat_ref.NVT_freq) );
    EXPECT_TRUE( float_complete_eq(stat_file.NPT_freq, stat_ref.NPT_freq) );
    EXPECT_TRUE( float_complete_eq(stat_file.v_press , stat_ref.v_press ) );
    EXPECT_EQ(stat_file.w_coef  , stat_ref.w_coef   );
    EXPECT_EQ(stat_file.x_nhc   , stat_ref.x_nhc   );
    EXPECT_EQ(stat_file.v_nhc   , stat_ref.v_nhc   );

    //--- check atom info
    std::vector<Atom_FP> v_atom;
    std::vector<Atom_FP> v_atom_file;
    allGatherAtom(atom,      v_atom);
    allGatherAtom(atom_file, v_atom_file);
    std::sort(v_atom.begin(),
              v_atom.end(),
              [](const Atom_FP &a, const Atom_FP &b){ return ( a.getAtomID() < b.getAtomID() ); });
    std::sort(v_atom_file.begin(),
              v_atom_file.end(),
              [](const Atom_FP &a, const Atom_FP &b){ return ( a.getAtomID() < b.getAtomID() ); });

    ASSERT_EQ(v_atom_file.size(), v_atom.size());
    for(size_t i=0; i<v_atom.size(); ++i){
        EXPECT_EQ(v_atom_file[i].getAtomID()  , v_atom[i].getAtomID()  );
        EXPECT_EQ(v_atom_file[i].getMolID()   , v_atom[i].getMolID()   );
        EXPECT_EQ(v_atom_file[i].getAtomType(), v_atom[i].getAtomType());
        EXPECT_EQ(v_atom_file[i].getMolType() , v_atom[i].getMolType() );
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getPos().x , v_atom[i].getPos().x  ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getPos().y , v_atom[i].getPos().y  ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getPos().z , v_atom[i].getPos().z  ) ) << " i = " << i;

        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getMass()  , v_atom[i].getMass()   ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getVel().x , v_atom[i].getVel().x  ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getVel().y , v_atom[i].getVel().y  ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getVel().z , v_atom[i].getVel().z  ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getCharge(), v_atom[i].getCharge() ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getVDW_D() , v_atom[i].getVDW_D()  ) ) << " i = " << i;
        EXPECT_TRUE( float_complete_eq(v_atom_file[i].getVDW_R() , v_atom[i].getVDW_R()  ) ) << " i = " << i;
    }
}



#include "gtest_main_mpi.hpp"
