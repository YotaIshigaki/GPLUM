//***************************************************************************************
//  This program is unit test of bond interaction.
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
//------ loading model parameter
#include "md_loading_model.hpp"

//--- common tool of unit test for force
#include "gtest_force_common.hpp"


namespace TEST_DEFS {
    const PS::S32 n_atom = 2;
    const PS::S64 n_loop = 1000;
    const PS::F64 range  = 3.2;

    const MD_DEFS::ID_type id_tgt = 1;

    const PS::S32 data_field_len = 8;

    const std::string harmonic_log_file{"test_bin/force_bond_harmonic_log.dat"};
    const std::string harmonic_ref_file{"unit_test/ref/force_bond_harmonic_ref.dat"};

    const std::string anharmonic_log_file{"test_bin/force_bond_anharmonic_log.dat"};
    const std::string anharmonic_ref_file{"unit_test/ref/force_bond_anharmonic_ref.dat"};
};

template <class Tpsys>
void test_move(Tpsys &atom){

    for(PS::S64 i=0; i<atom.getNumberOfParticleLocal(); ++i){
        PS::F64vec pos_tmp = Normalize::realPos( atom[i].getPos() );

        if(atom[i].getAtomID() == 1){
            //--- move for 2-body potential test (LJ, coulomb, bond)
            pos_tmp.x += TEST_DEFS::range/PS::F64(TEST_DEFS::n_loop);
        }

        pos_tmp = Normalize::normPos(pos_tmp);
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

void check_result(const std::vector<ForceData> &result,
                  const std::vector<ForceData> &ref    ){

    EXPECT_EQ(result.size(), ref.size());
    const PS::S32 n = std::min(result.size(), ref.size());

    for(PS::S32 i=0; i<n; ++i){
        EXPECT_EQ(result[i].count, ref[i].count) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].pos.x    , ref[i].pos.x    ) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].pos.y    , ref[i].pos.y    ) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].pos.z    , ref[i].pos.z    ) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].potential, ref[i].potential) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].force.x  , ref[i].force.x  ) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].force.y  , ref[i].force.y  ) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].force.z  , ref[i].force.z  ) << " i= " << i;
    }
}

template <class Tptcl, class Tdata>
void test_record(const Tptcl              &atom,
                 const PS::S32             count,
                       std::vector<Tdata> &logger){

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

    logger.push_back( ForceData{count,
                                Normalize::realPos( buf.getPos() - PS::F32vec{0.5, 0.5, 0.5}),
                                buf.getPotBond(),
                                buf.getForceIntra()} );
}


class TestForceBond :
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
    }
};

TEST_F(TestForceBond, harmonic){
    if(PS::Comm::getRank() == 0){
        //--- model parameters
        std::vector<std::string> param_line = {"Ow", "Hw", "harmonic", "1.012", "1059.162", "0.0"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line,
                                  MODEL::coef_table.bond);
        MODEL::coef_table.mask_scaling[MolName::AA_wat_SPC_Fw] = MD_DEFS::MaskList{};

        //--- atom settings
        atom.setNumberOfParticleLocal(TEST_DEFS::n_atom);
        for(PS::S32 i=0; i<TEST_DEFS::n_atom; ++i){
            atom[i].setAtomID(i);
            atom[i].setMolID(i);
            atom[i].setMolType( MolName::AA_wat_SPC_Fw );
            atom[i].setCharge( 0.0 );
            atom[i].setVDW_R( 3.0 );
            atom[i].setVDW_D( 0.0 );
            atom[i].clear();
        }
        atom[0].setAtomType( AtomName::Ow );
        atom[1].setAtomType( AtomName::Hw );

        atom[0].bond.add(1);
        atom[1].bond.add(0);

        atom[0].setPos( PS::F64vec(0.5) );
        atom[1].setPos( PS::F64vec(0.5) + Normalize::normPos( PS::F64vec(0.05, 0.0, 0.0) ) );
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(force_log, TEST_DEFS::harmonic_log_file);
    load_log_file( force_ref, TEST_DEFS::harmonic_ref_file);

    check_result(force_log, force_ref);
}


TEST_F(TestForceBond, anharmonic){
    if(PS::Comm::getRank() == 0){
        //--- model parameters
        std::vector<std::string> param_line = {"Ow", "Hw", "anharmonic", "0.995", "116.09", "2.287"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_aSPC_Fw),
                                  param_line,
                                  MODEL::coef_table.bond);
        MODEL::coef_table.mask_scaling[MolName::AA_wat_aSPC_Fw] = MD_DEFS::MaskList{};

        //--- atom settings
        atom.setNumberOfParticleLocal(TEST_DEFS::n_atom);
        for(PS::S32 i=0; i<TEST_DEFS::n_atom; ++i){
            atom[i].setAtomID(i);
            atom[i].setMolID(i);
            atom[i].setMolType( MolName::AA_wat_aSPC_Fw );
            atom[i].setCharge( 0.0 );
            atom[i].setVDW_R( 3.0 );
            atom[i].setVDW_D( 0.0 );
            atom[i].clear();
        }
        atom[0].setAtomType( AtomName::Ow );
        atom[1].setAtomType( AtomName::Hw );

        atom[0].bond.add(1);
        atom[1].bond.add(0);

        atom[0].setPos( PS::F64vec(0.5) );
        atom[1].setPos( PS::F64vec(0.5) + Normalize::normPos( PS::F64vec(0.05, 0.0, 0.0) ) );
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(force_log, TEST_DEFS::anharmonic_log_file);
    load_log_file( force_ref, TEST_DEFS::anharmonic_ref_file);

    check_result(force_log, force_ref);
}


#include "gtest_main_mpi.hpp"
