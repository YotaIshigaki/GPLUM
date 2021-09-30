//***************************************************************************************
//  This program is unit test of coulomb interaction.
//***************************************************************************************

#include <gtest/gtest.h>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

//--- external library for MD
#include <molecular_dynamics_ext.hpp>

//--- user defined headers
//------ definition of data set
#include "unit.hpp"
#include "md_enum.hpp"
#include "md_defs.hpp"
#include "atom_class.hpp"
#include "md_coef_table.hpp"
//------ calculate interaction
#include "md_force.hpp"

//--- common tool of unit test for force
#include "gtest_force_common.hpp"


namespace TEST_DEFS {
    const PS::S32 n_atom = 2;
    const PS::S64 n_loop = 1000;

    const MD_DEFS::ID_type id_tgt = 1;

    const PS::S32 data_field_len = 8;

    const PS::F64 eps_abs = 1.e-4;
    const PS::F64 eps_rel = 1.e-3;

    const std::string test_log_file{"test_bin/force_coulomb_log.dat"};
    const std::string test_ref_file{"unit_test/ref/force_coulomb_ref.dat"};
};

void test_model_setting(){
    auto& mask_table = MODEL::coef_table.mask_scaling;

    mask_table[MolName::AA_Ar] = MD_DEFS::MaskList{};
}

template <class Tptcl>
void test_atom_setting(Tptcl &atom){
    atom.setNumberOfParticleLocal(TEST_DEFS::n_atom);

    //--- set initial position & parameters
    //------ for bond test
    for(PS::S32 i=0; i<TEST_DEFS::n_atom; ++i){
    atom[i].setAtomID(i);
    atom[i].setAtomType( AtomName::Ar );
    atom[i].setMolType(  MolName::AA_Ar );
    atom[i].setMolID(i);
    atom[i].setCharge( 0.0 );
    atom[i].setVDW_R( 0.5*3.81637 );
    atom[i].setVDW_D( 0.0 );
    atom[i].clear();
    }

    atom[0].setCharge( 1.0*Unit::coef_coulomb);
    atom[1].setCharge(-1.0*Unit::coef_coulomb);

    atom[0].setPos( PS::F64vec(0.5) );
    atom[1].setPos( PS::F64vec(0.5) + PS::F64vec(0.001, 0.0, 0.0) );
}

template <class Tpsys>
void test_move(Tpsys &atom){

    for(PS::S64 i=0; i<atom.getNumberOfParticleLocal(); ++i){
        PS::F64vec pos_tmp = atom[i].getPos();

        if(atom[i].getAtomID() == 1){
            //--- move for 2-body potential test (LJ, coulomb, bond)
            pos_tmp.x += 0.998/PS::F64(TEST_DEFS::n_loop);         //  setting test distance range (normalized)
        }

        pos_tmp = Normalize::periodicAdjustNorm(pos_tmp);
        atom[i].setPos(pos_tmp);
    }
}

class ForceData {
public:
    PS::S32    count;
    PS::F32vec pos;
    PS::F32    potential;
    PS::F32vec force;
};
std::ostream& operator << (std::ostream &s, const ForceData &d){
    s << std::setw(10) << d.count << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.pos.x << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.pos.y << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.pos.z << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.potential << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.force.x << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.force.y << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.force.z << "\n";
    return s;
}

void read_ref_data(const std::vector<std::string> &str_list,
                         std::vector<ForceData>   &data_list){

    if(str_list.size() < TEST_DEFS::data_field_len) return;
    if( !STR_TOOL::isInteger(str_list[0]) ) return;
    if( !STR_TOOL::isNumeric(str_list[1]) ) return;
    if( !STR_TOOL::isNumeric(str_list[2]) ) return;
    if( !STR_TOOL::isNumeric(str_list[3]) ) return;
    if( !STR_TOOL::isNumeric(str_list[4]) ) return;
    if( !STR_TOOL::isNumeric(str_list[5]) ) return;
    if( !STR_TOOL::isNumeric(str_list[6]) ) return;
    if( !STR_TOOL::isNumeric(str_list[7]) ) return;

    ForceData tmp;
    tmp.count = std::stoi(str_list[0]);
    tmp.pos   = PS::F32vec{ std::stof(str_list[1]),
                            std::stof(str_list[2]),
                            std::stof(str_list[3]) };
    tmp.potential = std::stof(str_list[4]);
    tmp.force     = PS::F32vec{ std::stof(str_list[5]),
                                std::stof(str_list[6]),
                                std::stof(str_list[7]) };

    data_list.push_back(tmp);
}

template <class Tptcl, class Tlog>
void test_record(const Tptcl   &atom,
                 const PS::S32  count,
                       Tlog    &logger){

    Atom_FP buf;
    PS::S32 data_proc = -1;

    for(PS::S32 i=0; i<atom.getNumberOfParticleLocal(); ++i){
        if(atom[i].getAtomID() == TEST_DEFS::id_tgt){
            buf       = atom[i];
            data_proc = PS::Comm::getRank();
        }
    }
    data_proc = PS::Comm::getMaxValue(data_proc);
    COMM_TOOL::broadcast(buf, data_proc);

    logger.push_back( ForceData{ count,
                                 Normalize::realPos( buf.getPos() - PS::F32vec{0.5, 0.5, 0.5} ),
                                 buf.getCharge()*buf.getPotCoulomb(),
                                 buf.getCharge()*buf.getFieldCoulomb() } );
}


class TestForceCoulomb :
public ::testing::Test {
public:
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    virtual void SetUp(){
        //--- initialize
        test_init(TEST_DEFS::n_loop);

        force_log.clear();
        force_ref.clear();

        atom.initialize();
        atom.setNumberOfParticleLocal(0);

        if(PS::Comm::getRank() == 0){
            test_model_setting();
            test_atom_setting(atom);
        }

        execute_force_calc(atom, dinfo, force, force_log, force_ref);

        //--- write test result
        write_log_file(force_log, TEST_DEFS::test_log_file);

        //--- load reference result
        load_log_file(force_ref, TEST_DEFS::test_ref_file);
    }
};

TEST_F(TestForceCoulomb, force){
    EXPECT_EQ(force_log.size(), force_ref.size());
    const PS::S32 n = std::min(force_log.size(), force_ref.size());

    for(PS::S32 i=0; i<n; ++i){
        EXPECT_EQ(force_log[i].count, force_ref[i].count) << " i= " << i;
        EXPECT_FLOAT_EQ(force_log[i].pos.x, force_ref[i].pos.x) << " i= " << i;
        EXPECT_FLOAT_EQ(force_log[i].pos.y, force_ref[i].pos.y) << " i= " << i;
        EXPECT_FLOAT_EQ(force_log[i].pos.z, force_ref[i].pos.z) << " i= " << i;
        EXPECT_TRUE( float_relative_eq(force_log[i].potential, force_ref[i].potential, TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
        EXPECT_TRUE( float_relative_eq(force_log[i].force.x  , force_ref[i].force.x  , TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
        EXPECT_TRUE( float_relative_eq(force_log[i].force.y  , force_ref[i].force.y  , TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
        EXPECT_TRUE( float_relative_eq(force_log[i].force.z  , force_ref[i].force.z  , TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
    }
}


#include "gtest_main_mpi.hpp"
