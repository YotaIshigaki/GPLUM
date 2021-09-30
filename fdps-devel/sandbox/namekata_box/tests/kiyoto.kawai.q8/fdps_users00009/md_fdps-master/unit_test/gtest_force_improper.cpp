//***************************************************************************************
//  This program is unit test of improper torsion interaction.
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
    const PS::S32 n_atom = 4;
    const PS::S64 n_loop = 1000;
    const PS::F32 range  = 2.0*Unit::pi;

    const MD_DEFS::ID_type id_tgt = 0;

    const PS::S32 data_field_len = 6;

    const std::string harmonic_log_file{"test_bin/force_improper_harmonic_log.dat"};
    const std::string harmonic_ref_file{"unit_test/ref/force_improper_harmonic_ref.dat"};

    const std::string OPLS_log_file{"test_bin/force_improper_OPLS_3_log.dat"};
    const std::string OPLS_ref_file{"unit_test/ref/force_improper_OPLS_3_ref.dat"};
};

template <class Tpsys>
void test_move(Tpsys &atom){

    for(PS::S64 i=0; i<atom.getNumberOfParticleLocal(); ++i){
        PS::F64vec pos_tmp = Normalize::realPos( atom[i].getPos() );

        if(atom[i].getAtomID() == 2 ||
           atom[i].getAtomID() == 3 ){
            //--- move for dihedral torsion test
            auto pos_local = pos_tmp - Normalize::realPos( PS::F64vec{0.5} );
                 pos_local = VEC_EXT::rot_y(pos_local, TEST_DEFS::range/PS::F64(TEST_DEFS::n_loop));
            pos_tmp = pos_local + Normalize::realPos( PS::F64vec{0.5} );
        }

        pos_tmp = Normalize::normPos(pos_tmp);
        pos_tmp = Normalize::periodicAdjustNorm(pos_tmp);
        atom[i].setPos(pos_tmp);
    }
}

class ForceData {
public:
    PS::S32    count;
    PS::F32    degree;
    PS::F32    potential;
    PS::F32vec force;
};
std::ostream& operator << (std::ostream &s, const ForceData &d){
    s << std::setw(10) << d.count << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.degree    << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.potential << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.force.x   << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.force.y   << " "
      << std::setw(15) << std::scientific << std::setprecision(8) << d.force.z   << "\n";
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

    ForceData tmp;
    tmp.count     = std::stoi(str_list[0]);
    tmp.degree    = std::stof(str_list[1]);
    tmp.potential = std::stof(str_list[2]);
    tmp.force     = PS::F32vec{ std::stof(str_list[3]),
                                std::stof(str_list[4]),
                                std::stof(str_list[5]) };

    data_list.push_back(tmp);
}

void check_result(const std::vector<ForceData> &result,
                  const std::vector<ForceData> &ref    ){

    EXPECT_EQ(result.size(), ref.size());
    const PS::S32 n = std::min(result.size(), ref.size());

    for(PS::S32 i=0; i<n; ++i){
        EXPECT_EQ(result[i].count, ref[i].count) << " i= " << i;
        EXPECT_FLOAT_EQ(result[i].degree   , ref[i].degree   ) << " i= " << i;
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
                                static_cast<PS::F32>(180.0*(TEST_DEFS::range/Unit::pi)*PS::F32(count)/PS::F32(TEST_DEFS::n_loop)),
                                buf.getPotTorsion(),
                                buf.getForceIntra()} );
}


class TestForceImproper :
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

TEST_F(TestForceImproper, harmonic){
    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.bond.clear();
        MODEL::coef_table.angle.clear();
        MODEL::coef_table.torsion.clear();

        std::vector<std::string> param_line;
        param_line = {"CT", "CT", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_C6H5_CH3),
                                  param_line,
                                  MODEL::coef_table.bond);

        param_line = {"CT", "CT", "CT", "none"};
        MODEL::loading_param_angle(ENUM::what(MolName::AA_C6H5_CH3),
                                   param_line,
                                   MODEL::coef_table.angle);

        param_line = {"improper", "CT", "CT", "CT", "CT", "cos", "0.0", "3.0", "3"};
        MODEL::loading_param_torsion(ENUM::what(MolName::AA_C6H5_CH3),
                                     param_line,
                                     MODEL::coef_table.torsion);

        MODEL::coef_table.mask_scaling[MolName::AA_C6H5_CH3] = MD_DEFS::MaskList{};

        //--- atom settings
        atom.setNumberOfParticleLocal(TEST_DEFS::n_atom);
        for(PS::S32 i=0; i<TEST_DEFS::n_atom; ++i){
            atom[i].setAtomID(i);
            atom[i].setMolID(i);
            atom[i].setAtomType( AtomName::CT );
            atom[i].setMolType(  MolName::AA_C6H5_CH3 );
            atom[i].setCharge( 0.0 );
            atom[i].setVDW_R( 3.0 );
            atom[i].setVDW_D( 0.0 );
            atom[i].clear();
        }

        atom[0].bond.add(1);
        atom[1].bond.add(0);
        atom[1].bond.add(2);
        atom[1].bond.add(3);
        atom[2].bond.add(1);
        atom[3].bond.add(1);

        PS::F32vec R_ij = Normalize::normPos( PS::F32vec{1.0, 0.0, 0.0} );

        atom[0].setPos( PS::F64vec{0.5} + R_ij );
        atom[1].setPos( PS::F64vec{0.5} );
        atom[2].setPos( PS::F64vec{0.5} + VEC_EXT::rot_z(R_ij,  Unit::pi*(2.0/3.0) ) );
        atom[3].setPos( PS::F64vec{0.5} + VEC_EXT::rot_z(R_ij, -Unit::pi*(2.0/3.0) ) );
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(force_log, TEST_DEFS::harmonic_log_file);
    load_log_file( force_ref, TEST_DEFS::harmonic_ref_file);

    check_result(force_log, force_ref);
}


TEST_F(TestForceImproper, OPLS3){
    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.bond.clear();
        MODEL::coef_table.angle.clear();
        MODEL::coef_table.torsion.clear();

        std::vector<std::string> param_line;
        param_line = {"CT", "CT", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_C6H5_CH3),
                                  param_line,
                                  MODEL::coef_table.bond);

        param_line = {"CT", "CT", "CT", "none"};
        MODEL::loading_param_angle(ENUM::what(MolName::AA_C6H5_CH3),
                                   param_line,
                                   MODEL::coef_table.angle);

        param_line = {"improper", "CT", "CT", "CT", "CT", "OPLS_3", "1.714", "-0.157", "0.279"};
        MODEL::loading_param_torsion(ENUM::what(MolName::AA_C6H5_CH3),
                                     param_line,
                                     MODEL::coef_table.torsion);

        MODEL::coef_table.mask_scaling[MolName::AA_C6H5_CH3] = MD_DEFS::MaskList{};

        //--- atom settings
        atom.setNumberOfParticleLocal(TEST_DEFS::n_atom);
        for(PS::S32 i=0; i<TEST_DEFS::n_atom; ++i){
            atom[i].setAtomID(i);
            atom[i].setMolID(i);
            atom[i].setAtomType( AtomName::CT );
            atom[i].setMolType(  MolName::AA_C6H5_CH3 );
            atom[i].setCharge( 0.0 );
            atom[i].setVDW_R( 3.0 );
            atom[i].setVDW_D( 0.0 );
            atom[i].clear();
        }

        atom[0].bond.add(1);
        atom[1].bond.add(0);
        atom[1].bond.add(2);
        atom[1].bond.add(3);
        atom[2].bond.add(1);
        atom[3].bond.add(1);

        PS::F32vec R_ij = Normalize::normPos( PS::F32vec{1.0, 0.0, 0.0} );

        atom[0].setPos( PS::F64vec{0.5} + R_ij );
        atom[1].setPos( PS::F64vec{0.5} );
        atom[2].setPos( PS::F64vec{0.5} + VEC_EXT::rot_z(R_ij,  Unit::pi*(2.0/3.0) ) );
        atom[3].setPos( PS::F64vec{0.5} + VEC_EXT::rot_z(R_ij, -Unit::pi*(2.0/3.0) ) );
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(force_log, TEST_DEFS::OPLS_log_file);
    load_log_file( force_ref, TEST_DEFS::OPLS_ref_file);

    check_result(force_log, force_ref);
}


#include "gtest_main_mpi.hpp"
