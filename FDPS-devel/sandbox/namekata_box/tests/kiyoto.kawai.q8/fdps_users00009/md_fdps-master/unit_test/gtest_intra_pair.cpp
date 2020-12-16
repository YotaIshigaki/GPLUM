//=======================================================================================
//  This is unit test of IntraPair:: manager tools.
//     module location: ./src/intra_pair.hpp
//=======================================================================================

#include <gtest/gtest.h>

#include <ostream>
#include <tuple>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "md_defs.hpp"

namespace TEST_DEFS {

    constexpr size_t max_sp_pair     = 2;
    constexpr size_t max_all_connect = MD_DEFS::max_bond + max_sp_pair;

};

//--- support function
template <class T>
std::string to_string(const std::vector<T> &vec){
    std::ostringstream oss;
    bool isFirst = true;
    for(const auto& e : vec){
        if( !isFirst ) oss << " ,";
        oss << e;
        isFirst = false;
    }
    return oss.str();
}

template <class T>
::testing::AssertionResult isFoundInVector(const T &elem, const std::vector<T> &vec){
    const auto itr = std::find(vec.begin(), vec.end(), elem);
    if(itr == vec.end()){
        std::ostringstream oss;
        oss << "\n"
            << "    element " << elem << " was not found.\n"
            << "    vector: " << to_string(vec) << "\n";
        return ::testing::AssertionFailure() << oss.str();
    } else {
        return ::testing::AssertionSuccess();
    }
}

template <class T>
::testing::AssertionResult compare_2_vectors(const std::vector<T> &test_tgt,
                                             const std::vector<T> &reference,
                                             const std::string    &comment){

    std::vector<T> same_elem{},
                   not_found_elem{},
                   missed_elem{};

    for(const auto& e : test_tgt){
        const auto itr = std::find(reference.begin(), reference.end(), e);
        if(itr == reference.end()){
            not_found_elem.push_back(e);
        } else {
            same_elem.push_back(e);
        }
    }

    for(const auto& r : reference){
        const auto itr = std::find(test_tgt.begin(), test_tgt.end(), r);
        if(itr == test_tgt.end()){
            missed_elem.push_back(r);
        } else {
            //same_elem.push_back(r);
        }
    }

    //--- judge
    if(test_tgt.size() == reference.size() &&
       not_found_elem.size() == 0 &&
       missed_elem.size()    == 0    ){
        return ::testing::AssertionSuccess();
    } else {
        std::ostringstream oss;
        oss << "\n"
            << "  size:"
            << "    target= "  << test_tgt.size()
            << ", reference= " << reference.size() << "\n";
        oss << "  not found in reference:";
        if(not_found_elem.size() == 0){
            oss << " none." << "\n";
        } else {
            oss << "\n   " << to_string(not_found_elem) << "\n";
        }
        oss << "  missed in target:";
        if(missed_elem.size() == 0){
            oss << " none." << "\n";
        } else {
            oss << "\n   " << to_string(missed_elem) << "\n";
        }
        oss << "  same element:";
        if(same_elem.size() == 0){
            oss << " none." << "\n";
        } else {
            oss << "\n   " << to_string(same_elem) << "\n";
        }
        oss << "  " << comment << "\n";
        return ::testing::AssertionFailure() << oss.str();
    }
}

//--- atom data sample
class Atom {
public:
    MD_DEFS::ID_type id{-1};
    PS::F64vec       pos{0.0, 0.0, 0.0};

    MD_EXT::basic_connect<MD_DEFS::ID_type,
                          MD_DEFS::max_bond> bond;

    MD_DEFS::ID_type getId()  const { return this->id;  }
    PS::F64vec       getPos() const { return this->pos; }

    void setPos(const PS::F64vec &pos_new){ this->pos = pos_new; }

    void initData(const MD_DEFS::ID_type id){
        this->id = id;

        int x_mod = id % 2;
        int y_mod = id % 4;
        int z_mod = id % 8;

        this->pos = {0.49, 0.49, 0.49};
        if(x_mod >  0) this->pos.x += 0.02;
        if(y_mod >= 2) this->pos.y += 0.02;
        if(z_mod >= 4) this->pos.z += 0.02;
    }

    template <class Tptcl>
    void copyFromFP(const Tptcl &ptcl){
        this->id   = ptcl.getId();
        this->pos  = ptcl.getPos();
        this->bond = ptcl.bond;
    }
    void clear(){}

    static PS::F32 R_cut;
    static PS::F32 getRSearch() { return Atom::R_cut; }

    //--- static data for intra pair table
    static std::unordered_map< MD_DEFS::ID_type,
                               std::tuple<MD_DEFS::MaskList,
                                          MD_DEFS::AngleList,
                                          MD_DEFS::TorsionList,
                                          MD_DEFS::TorsionList> > intra_pair_table;

    MD_DEFS::MaskList&    mask_list()    { return std::get<0>(this->intra_pair_table[this->getId()]); }
    MD_DEFS::AngleList&   angle_list()   { return std::get<1>(this->intra_pair_table[this->getId()]); }
    MD_DEFS::TorsionList& dihedral_list(){ return std::get<2>(this->intra_pair_table[this->getId()]); }
    MD_DEFS::TorsionList& improper_list(){ return std::get<3>(this->intra_pair_table[this->getId()]); }

    void clear_intra_list(){
        this->mask_list().clear();
        this->angle_list().clear();
        this->dihedral_list().clear();
        this->improper_list().clear();
    }
};
PS::F32 Atom::R_cut = 0.3;
std::unordered_map< MD_DEFS::ID_type,
                    std::tuple<MD_DEFS::MaskList,
                               MD_DEFS::AngleList,
                               MD_DEFS::TorsionList,
                               MD_DEFS::TorsionList> > Atom::intra_pair_table;

//--- common data
template <class Tptcl>
struct DataBasic {

    //--- test target data
    PS::DomainInfo           dinfo;
    PS::ParticleSystem<Tptcl> atom;

    typename PS::TreeForForceShort<Tptcl, Tptcl, Tptcl>::Scatter tree;

    MD_DEFS::MaskList mask_param{ MD_DEFS::IntraMask(-1, 0.0, 0.0),
                                  MD_DEFS::IntraMask(-1, 0.1, 0.1),
                                  MD_DEFS::IntraMask(-1, 0.5, 0.5) };

    //--- table maker object
    struct GetBond {
        MD_EXT::basic_connect<MD_DEFS::ID_type,
                              MD_DEFS::max_bond> operator () (const Atom &atom){
            return atom.bond;
        }
    };
    IntraPair::IntraMaskMaker<  MD_DEFS::ID_type, GetBond> intra_mask_maker;
    IntraPair::AngleListMaker<  MD_DEFS::ID_type, GetBond> angle_list_maker;
    IntraPair::TorsionListMaker<MD_DEFS::ID_type, GetBond> torsion_list_maker;

    //--- reference result
    std::vector<std::vector<MD_DEFS::IntraMask>>  ref_mask_result;
    std::vector<std::vector<MD_DEFS::AngleSet>>   ref_angle_result;
    std::vector<std::vector<MD_DEFS::TorsionSet>> ref_dihedral_result;
    std::vector<std::vector<MD_DEFS::TorsionSet>> ref_improper_result;

    void init_FDPS_obj(const PS::S64 n){
        dinfo.initialize(0.3);
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
        dinfo.setPosRootDomain(PS::F32vec{0.0, 0.0, 0.0}, PS::F32vec{1.0, 1.0, 1.0});

        atom.initialize();
        atom.setNumberOfParticleLocal(0);

        Atom::intra_pair_table.clear();
        Atom::intra_pair_table.max_load_factor(0.7);

        tree.initialize(n, 0.5, 8, 64);
    }

    void make_intra_list(){
        //--- get neighbor EP_intra information (do not calculate force)
        this->tree.calcForceAll( IntraPair::dummy_func{},
                                 atom,
                                 dinfo );

        PS::S32 n_local = atom.getNumberOfParticleLocal();
        for(PS::S32 i=0; i<n_local; ++i){
            atom[i].clear_intra_list();
            intra_mask_maker(  atom[i], tree, mask_param, atom[i].mask_list());
            angle_list_maker(  atom[i], tree, atom[i].angle_list());
            torsion_list_maker(atom[i], tree, atom[i].dihedral_list(), atom[i].improper_list());
        }
    }

    void check_result(){
        PS::S32 n_local = atom.getNumberOfParticleLocal();
        PS::S32 i_proc  = PS::Comm::getRank();
        for(PS::S32 i=0; i<n_local; ++i){
            const auto id = atom[i].getId();
            std::ostringstream oss;
            oss << "Proc: " << i_proc << ", ID: " << id;

            EXPECT_TRUE( compare_2_vectors( atom[i].mask_list(),
                                            ref_mask_result.at(id),
                                            "mask_list, "     + oss.str() ) );
            EXPECT_TRUE( compare_2_vectors( atom[i].angle_list(),
                                            ref_angle_result.at(id),
                                            "angle_list, "    + oss.str() ) );
            EXPECT_TRUE( compare_2_vectors( atom[i].dihedral_list(),
                                            ref_dihedral_result.at(id),
                                            "dihedral_list, " + oss.str() ) );
            EXPECT_TRUE( compare_2_vectors( atom[i].improper_list(),
                                            ref_improper_result.at(id),
                                            "improper_list, " + oss.str() ) );
        }
    }
};

struct IntraPairBasic :
    public DataBasic<Atom>,
    public ::testing::Test {
    public:

    virtual void SetUp(){
        const PS::S64 n = 9;
        init_FDPS_obj(n);

        if(PS::Comm::getRank() == 0){
            atom.setNumberOfParticleLocal(n);
            for(PS::S64 i=0; i<n; ++i){
                atom[i].initData(i);
            }

            //--- define test connection structure
            //
            //    6
            //    |
            //  0-1-2-3-4-5
            //      |
            //      7-8

            atom[0].bond.add(1);

            atom[1].bond.add(0);
            atom[1].bond.add(2);
            atom[1].bond.add(6);

            atom[2].bond.add(1);
            atom[2].bond.add(3);
            atom[2].bond.add(7);

            atom[3].bond.add(2);
            atom[3].bond.add(4);

            atom[4].bond.add(3);
            atom[4].bond.add(5);

            atom[5].bond.add(4);

            atom[6].bond.add(1);

            atom[7].bond.add(2);
            atom[7].bond.add(8);

            atom[8].bond.add(7);
        }

        //--- reference result
        ref_mask_result.resize(n);
        ref_angle_result.resize(n);
        ref_dihedral_result.resize(n);
        ref_improper_result.resize(n);

        //------ for mask
        auto mask1 = mask_param.at(0);
        auto mask2 = mask_param.at(1);
        auto mask3 = mask_param.at(2);

        ref_mask_result.at(0) = { mask1.setId(1),
                                  mask2.setId(2),
                                  mask2.setId(6),
                                  mask3.setId(3),
                                  mask3.setId(7) };

        ref_mask_result.at(1) = { mask1.setId(0),
                                  mask1.setId(2),
                                  mask1.setId(6),
                                  mask2.setId(3),
                                  mask2.setId(7),
                                  mask3.setId(4),
                                  mask3.setId(8) };

        ref_mask_result.at(2) = { mask1.setId(1),
                                  mask1.setId(3),
                                  mask1.setId(7),
                                  mask2.setId(0),
                                  mask2.setId(6),
                                  mask2.setId(4),
                                  mask2.setId(8),
                                  mask3.setId(5), };

        ref_mask_result.at(3) = { mask1.setId(2),
                                  mask1.setId(4),
                                  mask2.setId(1),
                                  mask2.setId(5),
                                  mask2.setId(7),
                                  mask3.setId(0),
                                  mask3.setId(6),
                                  mask3.setId(8) };

        ref_mask_result.at(4) = { mask1.setId(3),
                                  mask1.setId(5),
                                  mask2.setId(2),
                                  mask3.setId(1),
                                  mask3.setId(7), };

        ref_mask_result.at(5) = { mask1.setId(4),
                                  mask2.setId(3),
                                  mask3.setId(2), };

        ref_mask_result.at(6) = { mask1.setId(1),
                                  mask2.setId(0),
                                  mask2.setId(2),
                                  mask3.setId(3),
                                  mask3.setId(7), };

        ref_mask_result.at(7) = { mask1.setId(2),
                                  mask1.setId(8),
                                  mask2.setId(1),
                                  mask2.setId(3),
                                  mask3.setId(0),
                                  mask3.setId(4),
                                  mask3.setId(6), };

        ref_mask_result.at(8) = { mask1.setId(7),
                                  mask2.setId(2),
                                  mask3.setId(1),
                                  mask3.setId(3) };

        //------ for angle
        ref_angle_result.at(0) = { MD_DEFS::AngleSet{0, 1, 2},
                                   MD_DEFS::AngleSet{0, 1, 6} };

        ref_angle_result.at(1) = { MD_DEFS::AngleSet{0, 1, 2},
                                   MD_DEFS::AngleSet{0, 1, 6},
                                   MD_DEFS::AngleSet{1, 2, 3},
                                   MD_DEFS::AngleSet{1, 2, 7},
                                   MD_DEFS::AngleSet{2, 1, 6} };

        ref_angle_result.at(2) = { MD_DEFS::AngleSet{2, 1, 0},
                                   MD_DEFS::AngleSet{2, 1, 6},
                                   MD_DEFS::AngleSet{1, 2, 3},
                                   MD_DEFS::AngleSet{1, 2, 7},
                                   MD_DEFS::AngleSet{2, 3, 4},
                                   MD_DEFS::AngleSet{2, 7, 8},
                                   MD_DEFS::AngleSet{3, 2, 7} };

        ref_angle_result.at(3) = { MD_DEFS::AngleSet{3, 2, 1},
                                   MD_DEFS::AngleSet{3, 2, 7},
                                   MD_DEFS::AngleSet{2, 3, 4},
                                   MD_DEFS::AngleSet{3, 4, 5} };

        ref_angle_result.at(4) = { MD_DEFS::AngleSet{4, 3, 2},
                                   MD_DEFS::AngleSet{3, 4, 5} };

        ref_angle_result.at(5) = { MD_DEFS::AngleSet{5, 4, 3} };

        ref_angle_result.at(6) = { MD_DEFS::AngleSet{6, 1, 0},
                                   MD_DEFS::AngleSet{6, 1, 2} };

        ref_angle_result.at(7) = { MD_DEFS::AngleSet{7, 2, 1},
                                   MD_DEFS::AngleSet{7, 2, 3},
                                   MD_DEFS::AngleSet{2, 7, 8} };

        ref_angle_result.at(8) = { MD_DEFS::AngleSet{8, 7, 2} };

        //------ for dihedral torsion
        ref_dihedral_result.at(0) = { MD_DEFS::TorsionSet{0, 1, 2, 3},
                                      MD_DEFS::TorsionSet{0, 1, 2, 7} };

        ref_dihedral_result.at(1) = { MD_DEFS::TorsionSet{1, 2, 3, 4},
                                      MD_DEFS::TorsionSet{1, 2, 7, 8},
                                      MD_DEFS::TorsionSet{0, 1, 2, 3},
                                      MD_DEFS::TorsionSet{0, 1, 2, 7},
                                      MD_DEFS::TorsionSet{6, 1, 2, 3},
                                      MD_DEFS::TorsionSet{6, 1, 2, 7} };

        ref_dihedral_result.at(2) = { MD_DEFS::TorsionSet{2, 3, 4, 5},
                                      MD_DEFS::TorsionSet{3, 2, 1, 0},
                                      MD_DEFS::TorsionSet{7, 2, 1, 0},
                                      MD_DEFS::TorsionSet{3, 2, 1, 6},
                                      MD_DEFS::TorsionSet{7, 2, 1, 6},
                                      MD_DEFS::TorsionSet{1, 2, 3, 4},
                                      MD_DEFS::TorsionSet{7, 2, 3, 4},
                                      MD_DEFS::TorsionSet{1, 2, 7, 8},
                                      MD_DEFS::TorsionSet{3, 2, 7, 8} };

        ref_dihedral_result.at(3) = { MD_DEFS::TorsionSet{3, 2, 1, 0},
                                      MD_DEFS::TorsionSet{3, 2, 1, 6},
                                      MD_DEFS::TorsionSet{3, 2, 7, 8},
                                      MD_DEFS::TorsionSet{2, 3, 4, 5},
                                      MD_DEFS::TorsionSet{4, 3, 2, 1},
                                      MD_DEFS::TorsionSet{4, 3, 2, 7} };

        ref_dihedral_result.at(4) = { MD_DEFS::TorsionSet{4, 3, 2, 1},
                                      MD_DEFS::TorsionSet{4, 3, 2, 7},
                                      MD_DEFS::TorsionSet{5, 4, 3, 2} };

        ref_dihedral_result.at(5) = { MD_DEFS::TorsionSet{5, 4, 3, 2} };

        ref_dihedral_result.at(6) = { MD_DEFS::TorsionSet{6, 1, 2, 3},
                                      MD_DEFS::TorsionSet{6, 1, 2, 7} };

        ref_dihedral_result.at(7) = { MD_DEFS::TorsionSet{7, 2, 1, 0},
                                      MD_DEFS::TorsionSet{7, 2, 1, 6},
                                      MD_DEFS::TorsionSet{7, 2, 3, 4},
                                      MD_DEFS::TorsionSet{8, 7, 2, 1},
                                      MD_DEFS::TorsionSet{8, 7, 2, 3} };

        ref_dihedral_result.at(8) = { MD_DEFS::TorsionSet{8, 7, 2, 1},
                                      MD_DEFS::TorsionSet{8, 7, 2, 3} };

        //------ for improper tortion
        ref_improper_result.at(0) = { MD_DEFS::TorsionSet{0, 1, 2, 6},
                                      MD_DEFS::TorsionSet{0, 1, 6, 2},
                                      MD_DEFS::TorsionSet{2, 0, 1, 6} };

        ref_improper_result.at(1) = { MD_DEFS::TorsionSet{2, 1, 0, 6},
                                      MD_DEFS::TorsionSet{0, 1, 2, 6},
                                      MD_DEFS::TorsionSet{0, 1, 6, 2},
                                      MD_DEFS::TorsionSet{3, 1, 2, 7},
                                      MD_DEFS::TorsionSet{1, 2, 3, 7},
                                      MD_DEFS::TorsionSet{1, 2, 7, 3} };

        ref_improper_result.at(2) = { MD_DEFS::TorsionSet{3, 2, 1, 7},
                                      MD_DEFS::TorsionSet{1, 2, 3, 7},
                                      MD_DEFS::TorsionSet{1, 2, 7, 3},
                                      MD_DEFS::TorsionSet{0, 2, 1, 6},
                                      MD_DEFS::TorsionSet{2, 1, 0, 6},
                                      MD_DEFS::TorsionSet{2, 1, 6, 0} };

        ref_improper_result.at(3) = { MD_DEFS::TorsionSet{1, 3, 2, 7},
                                      MD_DEFS::TorsionSet{3, 2, 1, 7},
                                      MD_DEFS::TorsionSet{3, 2, 7, 1} };

        ref_improper_result.at(4) = {};
        ref_improper_result.at(5) = {};

        ref_improper_result.at(6) = { MD_DEFS::TorsionSet{0, 6, 1, 2},
                                      MD_DEFS::TorsionSet{6, 1, 0, 2},
                                      MD_DEFS::TorsionSet{6, 1, 2, 0} };

        ref_improper_result.at(7) = { MD_DEFS::TorsionSet{1, 7, 2, 3},
                                      MD_DEFS::TorsionSet{7, 2, 1, 3},
                                      MD_DEFS::TorsionSet{7, 2, 3, 1} };

        ref_improper_result.at(8) = {};
    }

};

struct IntraPairCircular :
    public DataBasic<Atom>,
    public ::testing::Test {
    public:

    virtual void SetUp(){
        const PS::S64 n = 8;
        init_FDPS_obj(n);

        if(PS::Comm::getRank() == 0){
            atom.setNumberOfParticleLocal(n);
            for(PS::S64 i=0; i<n; ++i){
                atom[i].initData(i);
            }
            //--- circulation structure test
            //
            //  0-3
            //  | |
            //  1-2-4-7
            //      | |
            //      5-6

            atom[0].bond.add(1);
            atom[0].bond.add(3);

            atom[1].bond.add(0);
            atom[1].bond.add(2);

            atom[2].bond.add(1);
            atom[2].bond.add(3);
            atom[2].bond.add(4);

            atom[3].bond.add(0);
            atom[3].bond.add(2);

            atom[4].bond.add(2);
            atom[4].bond.add(5);
            atom[4].bond.add(7);

            atom[5].bond.add(4);
            atom[5].bond.add(6);

            atom[6].bond.add(5);
            atom[6].bond.add(7);

            atom[7].bond.add(4);
            atom[7].bond.add(6);
        }

        //--- reference result
        ref_mask_result.resize(n);
        ref_angle_result.resize(n);
        ref_dihedral_result.resize(n);
        ref_improper_result.resize(n);

        //------ for mask
        auto mask1 = mask_param.at(0);
        auto mask2 = mask_param.at(1);
        auto mask3 = mask_param.at(2);

        ref_mask_result.at(0) = { mask1.setId(1),
                                  mask1.setId(3),
                                  mask2.setId(2),
                                  mask3.setId(4) };

        ref_mask_result.at(1) = { mask1.setId(0),
                                  mask1.setId(2),
                                  mask2.setId(3),
                                  mask2.setId(4),
                                  mask3.setId(5),
                                  mask3.setId(7) };

        ref_mask_result.at(2) = { mask1.setId(1),
                                  mask1.setId(3),
                                  mask1.setId(4),
                                  mask2.setId(0),
                                  mask2.setId(5),
                                  mask2.setId(7),
                                  mask3.setId(6) };

        ref_mask_result.at(3) = { mask1.setId(0),
                                  mask1.setId(2),
                                  mask2.setId(1),
                                  mask2.setId(4),
                                  mask3.setId(5),
                                  mask3.setId(7) };

        ref_mask_result.at(4) = { mask1.setId(2),
                                  mask1.setId(5),
                                  mask1.setId(7),
                                  mask2.setId(1),
                                  mask2.setId(3),
                                  mask2.setId(6),
                                  mask3.setId(0) };

        ref_mask_result.at(5) = { mask1.setId(4),
                                  mask1.setId(6),
                                  mask2.setId(2),
                                  mask2.setId(7),
                                  mask3.setId(1),
                                  mask3.setId(3) };

        ref_mask_result.at(6) = { mask1.setId(5),
                                  mask1.setId(7),
                                  mask2.setId(4),
                                  mask3.setId(2) };

        ref_mask_result.at(7) = { mask1.setId(4),
                                  mask1.setId(6),
                                  mask2.setId(2),
                                  mask2.setId(5),
                                  mask3.setId(1),
                                  mask3.setId(3) };

        //------ for angle
        ref_angle_result.at(0) = { MD_DEFS::AngleSet{0, 1, 2},
                                   MD_DEFS::AngleSet{0, 3, 2},
                                   MD_DEFS::AngleSet{1, 0, 3} };

        ref_angle_result.at(1) = { MD_DEFS::AngleSet{1, 0, 3},
                                   MD_DEFS::AngleSet{1, 2, 3},
                                   MD_DEFS::AngleSet{1, 2, 4},
                                   MD_DEFS::AngleSet{0, 1, 2} };

        ref_angle_result.at(2) = { MD_DEFS::AngleSet{2, 1, 0},
                                   MD_DEFS::AngleSet{2, 3, 0},
                                   MD_DEFS::AngleSet{2, 4, 5},
                                   MD_DEFS::AngleSet{2, 4, 7},
                                   MD_DEFS::AngleSet{1, 2, 3},
                                   MD_DEFS::AngleSet{1, 2, 4},
                                   MD_DEFS::AngleSet{3, 2, 4} };

        ref_angle_result.at(3) = { MD_DEFS::AngleSet{3, 0, 1},
                                   MD_DEFS::AngleSet{3, 2, 1},
                                   MD_DEFS::AngleSet{3, 2, 4},
                                   MD_DEFS::AngleSet{0, 3, 2} };

        ref_angle_result.at(4) = { MD_DEFS::AngleSet{4, 2, 1},
                                   MD_DEFS::AngleSet{4, 2, 3},
                                   MD_DEFS::AngleSet{4, 5, 6},
                                   MD_DEFS::AngleSet{4, 7, 6},
                                   MD_DEFS::AngleSet{2, 4, 5},
                                   MD_DEFS::AngleSet{2, 4, 7},
                                   MD_DEFS::AngleSet{5, 4, 7} };

        ref_angle_result.at(5) = { MD_DEFS::AngleSet{5, 4, 2},
                                   MD_DEFS::AngleSet{5, 4, 7},
                                   MD_DEFS::AngleSet{5, 6, 7},
                                   MD_DEFS::AngleSet{4, 5, 6} };

        ref_angle_result.at(6) = { MD_DEFS::AngleSet{6, 5, 4},
                                   MD_DEFS::AngleSet{6, 7, 4},
                                   MD_DEFS::AngleSet{5, 6, 7} };

        ref_angle_result.at(7) = { MD_DEFS::AngleSet{7, 4, 2},
                                   MD_DEFS::AngleSet{7, 4, 5},
                                   MD_DEFS::AngleSet{7, 6, 5},
                                   MD_DEFS::AngleSet{4, 7, 6} };

        //------ for dihedral torsion
        ref_dihedral_result.at(0) = { MD_DEFS::TorsionSet{0, 1, 2, 3},
                                      MD_DEFS::TorsionSet{0, 1, 2, 4},
                                      MD_DEFS::TorsionSet{0, 3, 2, 1},
                                      MD_DEFS::TorsionSet{0, 3, 2, 4},
                                      MD_DEFS::TorsionSet{1, 0, 3, 2},
                                      MD_DEFS::TorsionSet{3, 0, 1, 2} };

        ref_dihedral_result.at(1) = { MD_DEFS::TorsionSet{1, 0, 3, 2},
                                      MD_DEFS::TorsionSet{1, 2, 3, 0},
                                      MD_DEFS::TorsionSet{1, 2, 4, 5},
                                      MD_DEFS::TorsionSet{1, 2, 4, 7},
                                      MD_DEFS::TorsionSet{0, 1, 2, 3},
                                      MD_DEFS::TorsionSet{0, 1, 2, 4},
                                      MD_DEFS::TorsionSet{2, 1, 0, 3} };

        ref_dihedral_result.at(2) = { MD_DEFS::TorsionSet{2, 1, 0, 3},
                                      MD_DEFS::TorsionSet{2, 3, 0, 1},
                                      MD_DEFS::TorsionSet{2, 4, 5, 6},
                                      MD_DEFS::TorsionSet{2, 4, 7, 6},
                                      MD_DEFS::TorsionSet{3, 2, 1, 0},
                                      MD_DEFS::TorsionSet{4, 2, 1, 0},
                                      MD_DEFS::TorsionSet{1, 2, 3, 0},
                                      MD_DEFS::TorsionSet{4, 2, 3, 0},
                                      MD_DEFS::TorsionSet{1, 2, 4, 5},
                                      MD_DEFS::TorsionSet{3, 2, 4, 5},
                                      MD_DEFS::TorsionSet{1, 2, 4, 7},
                                      MD_DEFS::TorsionSet{3, 2, 4, 7} };

        ref_dihedral_result.at(3) = { MD_DEFS::TorsionSet{3, 0, 1, 2},
                                      MD_DEFS::TorsionSet{3, 2, 1, 0},
                                      MD_DEFS::TorsionSet{3, 2, 4, 5},
                                      MD_DEFS::TorsionSet{3, 2, 4, 7},
                                      MD_DEFS::TorsionSet{2, 3, 0, 1},
                                      MD_DEFS::TorsionSet{0, 3, 2, 1},
                                      MD_DEFS::TorsionSet{0, 3, 2, 4} };

        ref_dihedral_result.at(4) = { MD_DEFS::TorsionSet{4, 2, 1, 0},
                                      MD_DEFS::TorsionSet{4, 2, 3, 0},
                                      MD_DEFS::TorsionSet{4, 5, 6, 7},
                                      MD_DEFS::TorsionSet{4, 7, 6, 5},
                                      MD_DEFS::TorsionSet{5, 4, 2, 1},
                                      MD_DEFS::TorsionSet{7, 4, 2, 1},
                                      MD_DEFS::TorsionSet{5, 4, 2, 3},
                                      MD_DEFS::TorsionSet{7, 4, 2, 3},
                                      MD_DEFS::TorsionSet{2, 4, 5, 6},
                                      MD_DEFS::TorsionSet{7, 4, 5, 6},
                                      MD_DEFS::TorsionSet{2, 4, 7, 6},
                                      MD_DEFS::TorsionSet{5, 4, 7, 6} };

        ref_dihedral_result.at(5) = { MD_DEFS::TorsionSet{5, 4, 2, 1},
                                      MD_DEFS::TorsionSet{5, 4, 2, 3},
                                      MD_DEFS::TorsionSet{5, 4, 7, 6},
                                      MD_DEFS::TorsionSet{5, 6, 7, 4},
                                      MD_DEFS::TorsionSet{6, 5, 4, 2},
                                      MD_DEFS::TorsionSet{6, 5, 4, 7},
                                      MD_DEFS::TorsionSet{4, 5, 6, 7} };

        ref_dihedral_result.at(6) = { MD_DEFS::TorsionSet{6, 5, 4, 2},
                                      MD_DEFS::TorsionSet{6, 5, 4, 7},
                                      MD_DEFS::TorsionSet{6, 7, 4, 2},
                                      MD_DEFS::TorsionSet{6, 7, 4, 5},
                                      MD_DEFS::TorsionSet{7, 6, 5, 4},
                                      MD_DEFS::TorsionSet{5, 6, 7, 4} };

        ref_dihedral_result.at(7) = { MD_DEFS::TorsionSet{7, 4, 2, 1},
                                      MD_DEFS::TorsionSet{7, 4, 2, 3},
                                      MD_DEFS::TorsionSet{7, 4, 5, 6},
                                      MD_DEFS::TorsionSet{7, 6, 5, 4},
                                      MD_DEFS::TorsionSet{6, 7, 4, 2},
                                      MD_DEFS::TorsionSet{6, 7, 4, 5},
                                      MD_DEFS::TorsionSet{4, 7, 6, 5} };

        //------ for improper tortion
        ref_improper_result.at(0) = {};

        ref_improper_result.at(1) = { MD_DEFS::TorsionSet{1, 2, 3, 4},
                                      MD_DEFS::TorsionSet{1, 2, 4, 3},
                                      MD_DEFS::TorsionSet{3, 1, 2, 4} };

        ref_improper_result.at(2) = { MD_DEFS::TorsionSet{3, 2, 1, 4},
                                      MD_DEFS::TorsionSet{1, 2, 3, 4},
                                      MD_DEFS::TorsionSet{1, 2, 4, 3},
                                      MD_DEFS::TorsionSet{5, 2, 4, 7},
                                      MD_DEFS::TorsionSet{2, 4, 5, 7},
                                      MD_DEFS::TorsionSet{2, 4, 7, 5} };

        ref_improper_result.at(3) = { MD_DEFS::TorsionSet{1, 3, 2, 4},
                                      MD_DEFS::TorsionSet{3, 2, 1, 4},
                                      MD_DEFS::TorsionSet{3, 2, 4, 1} };

        ref_improper_result.at(4) = { MD_DEFS::TorsionSet{1, 4, 2, 3},
                                      MD_DEFS::TorsionSet{4, 2, 1, 3},
                                      MD_DEFS::TorsionSet{4, 2, 3, 1},
                                      MD_DEFS::TorsionSet{5, 4, 2, 7},
                                      MD_DEFS::TorsionSet{2, 4, 5, 7},
                                      MD_DEFS::TorsionSet{2, 4, 7, 5} };

        ref_improper_result.at(5) = { MD_DEFS::TorsionSet{5, 4, 2, 7},
                                      MD_DEFS::TorsionSet{5, 4, 7, 2},
                                      MD_DEFS::TorsionSet{2, 5, 4, 7} };

        ref_improper_result.at(6) = {};

        ref_improper_result.at(7) = { MD_DEFS::TorsionSet{7, 4, 2, 5},
                                      MD_DEFS::TorsionSet{7, 4, 5, 2},
                                      MD_DEFS::TorsionSet{2, 7, 4, 5} };
    }

};


struct IntraPairException :
    public DataBasic<Atom>,
    public ::testing::Test {
    public:

    virtual void SetUp(){
        const PS::S64 n = 10;
        init_FDPS_obj(n);

        if(PS::Comm::getRank() == 0){
            atom.setNumberOfParticleLocal(n);
            for(PS::S64 i=0; i<n; ++i){
                atom[i].initData(i);
            }
            //--- throw exception test
            //
            //  1-1    (connected itself)

            atom[1].bond.add(1);
        }
    }

    void make_intra_list(){
        //--- get neighbor EP_intra information (do not calculate force)
        this->tree.calcForceAll( IntraPair::dummy_func{},
                                 atom,
                                 dinfo );

        PS::S32 n_local = atom.getNumberOfParticleLocal();
        for(PS::S32 i=0; i<n_local; ++i){
            atom[i].clear_intra_list();

            if(atom[i].getId() == 1){
                EXPECT_THROW( intra_mask_maker(  atom[i], tree, mask_param, atom[i].mask_list()),
                              std::invalid_argument);
                EXPECT_THROW( angle_list_maker(  atom[i], tree, atom[i].angle_list()),
                              std::invalid_argument);
                EXPECT_THROW( torsion_list_maker(atom[i], tree, atom[i].dihedral_list(), atom[i].improper_list()),
                              std::invalid_argument);
            }

        }
    }

};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(IntraPairBasic, singleProc){
    make_intra_list();
    check_result();
}

TEST_F(IntraPairBasic, withMPI){
    //--- reduce the atom data for all processes.
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    make_intra_list();
    check_result();
}

TEST_F(IntraPairCircular, singleProc){
    make_intra_list();
    check_result();
}

TEST_F(IntraPairCircular, withMPI){
    //--- reduce the atom data for all processes.
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    make_intra_list();
    check_result();
}

TEST_F(IntraPairException, singleProc){
    make_intra_list();
}

TEST_F(IntraPairException, withMPI){
    //--- reduce the atom data for all processes.
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    make_intra_list();
}


//--- test for SP-connection
class Atom_SP:
public Atom {
public:

    MD_EXT::basic_connect<MD_DEFS::ID_type, 2> sp_pair_A;
    MD_EXT::basic_connect<MD_DEFS::ID_type, 2> sp_pair_B;

    void copyFromFP(const Atom_SP &ptcl){
        this->id   = ptcl.getId();
        this->pos  = ptcl.getPos();

        this->bond      = ptcl.bond;
        this->sp_pair_A = ptcl.sp_pair_A;
        this->sp_pair_B = ptcl.sp_pair_B;
    }

    static std::unordered_map< MD_DEFS::ID_type,
                               MD_DEFS::MaskList> state_all_mask;
    static std::unordered_map< MD_DEFS::ID_type,
                               std::tuple<MD_DEFS::AngleList,
                                          MD_DEFS::TorsionList,
                                          MD_DEFS::TorsionList> > state_A_pair,
                                                                  state_B_pair;

    MD_DEFS::MaskList&    state_all_mask_list()  { return this->state_all_mask[this->getId()]; }

    MD_DEFS::AngleList&   state_A_angle_list()   { return std::get<0>(this->state_A_pair[this->getId()]); }
    MD_DEFS::TorsionList& state_A_dihedral_list(){ return std::get<1>(this->state_A_pair[this->getId()]); }
    MD_DEFS::TorsionList& state_A_improper_list(){ return std::get<2>(this->state_A_pair[this->getId()]); }

    MD_DEFS::AngleList&   state_B_angle_list()   { return std::get<0>(this->state_B_pair[this->getId()]); }
    MD_DEFS::TorsionList& state_B_dihedral_list(){ return std::get<1>(this->state_B_pair[this->getId()]); }
    MD_DEFS::TorsionList& state_B_improper_list(){ return std::get<2>(this->state_B_pair[this->getId()]); }

    void clear_intra_list(){
        this->state_all_mask_list().clear();

        this->state_A_angle_list().clear();
        this->state_A_dihedral_list().clear();
        this->state_A_improper_list().clear();

        this->state_B_angle_list().clear();
        this->state_B_dihedral_list().clear();
        this->state_B_improper_list().clear();
    }
};
std::unordered_map< MD_DEFS::ID_type,
                    MD_DEFS::MaskList > Atom_SP::state_all_mask;
std::unordered_map< MD_DEFS::ID_type,
                    std::tuple<MD_DEFS::AngleList,
                               MD_DEFS::TorsionList,
                               MD_DEFS::TorsionList> > Atom_SP::state_A_pair,
                                                       Atom_SP::state_B_pair;

//--- sample data for spacial connection case
struct IntraPairSPconnect :
public DataBasic<Atom_SP>,
public ::testing::Test {
public:

    //--- table maker for SP_connect
    struct GetAllConnect {
        MD_EXT::basic_connect<MD_DEFS::ID_type,
                              TEST_DEFS::max_all_connect> operator () (const Atom_SP &atom){
            return atom.bond + atom.sp_pair_A + atom.sp_pair_B;
        }
    };
    struct GetConnect_state_A {
        MD_EXT::basic_connect<MD_DEFS::ID_type,
                              TEST_DEFS::max_all_connect> operator () (const Atom_SP &atom){
            return atom.bond + atom.sp_pair_A;
        }
    };
    struct GetConnect_state_B {
        MD_EXT::basic_connect<MD_DEFS::ID_type,
                              TEST_DEFS::max_all_connect> operator () (const Atom_SP &atom){
            return atom.bond + atom.sp_pair_B;
        }
    };
    IntraPair::IntraMaskMaker<  MD_DEFS::ID_type, GetAllConnect>      make_intra_mask_all;
    IntraPair::AngleListMaker<  MD_DEFS::ID_type, GetConnect_state_A> make_angle_list_state_A;
    IntraPair::TorsionListMaker<MD_DEFS::ID_type, GetConnect_state_A> make_torsion_list_state_A;
    IntraPair::AngleListMaker<  MD_DEFS::ID_type, GetConnect_state_B> make_angle_list_state_B;
    IntraPair::TorsionListMaker<MD_DEFS::ID_type, GetConnect_state_B> make_torsion_list_state_B;

    //--- reference result for sp_pair
    std::vector<std::vector<MD_DEFS::IntraMask>>  ref_mask_all;
    std::vector<std::vector<MD_DEFS::AngleSet>>   ref_angle_state_A;
    std::vector<std::vector<MD_DEFS::TorsionSet>> ref_dihedral_state_A;
    std::vector<std::vector<MD_DEFS::TorsionSet>> ref_improper_state_A;
    std::vector<std::vector<MD_DEFS::AngleSet>>   ref_angle_state_B;
    std::vector<std::vector<MD_DEFS::TorsionSet>> ref_dihedral_state_B;
    std::vector<std::vector<MD_DEFS::TorsionSet>> ref_improper_state_B;

    void make_intra_list(){
        //--- get neighbor EP_intra information (do not calculate force)
        this->tree.calcForceAll( IntraPair::dummy_func{},
                                 atom,
                                 dinfo );

        PS::S32 n_local = atom.getNumberOfParticleLocal();
        for(PS::S32 i=0; i<n_local; ++i){
            atom[i].clear_intra_list();

            make_intra_mask_all( atom[i], tree, mask_param, atom[i].state_all_mask_list() );

            make_angle_list_state_A(   atom[i], tree, atom[i].state_A_angle_list() );
            make_torsion_list_state_A( atom[i], tree,
                                       atom[i].state_A_dihedral_list(),
                                       atom[i].state_A_improper_list() );

            make_angle_list_state_B(   atom[i], tree, atom[i].state_B_angle_list() );
            make_torsion_list_state_B( atom[i], tree,
                                       atom[i].state_B_dihedral_list(),
                                       atom[i].state_B_improper_list() );
        }
    }

    void check_result(){
        PS::S32 n_local = atom.getNumberOfParticleLocal();
        PS::S32 i_proc  = PS::Comm::getRank();
        for(PS::S32 i=0; i<n_local; ++i){
            const auto id = atom[i].getId();
            std::ostringstream oss;
            oss << "Proc: " << i_proc << ", ID: " << id;

            EXPECT_TRUE( compare_2_vectors( atom[i].state_all_mask_list(),
                                            ref_mask_all.at(id),
                                            "all_mask_list, " + oss.str() ) );

            EXPECT_TRUE( compare_2_vectors( atom[i].state_A_angle_list(),
                                            ref_angle_state_A.at(id),
                                            "state_A_angle_list, " + oss.str() ) );
            EXPECT_TRUE( compare_2_vectors( atom[i].state_A_dihedral_list(),
                                            ref_dihedral_state_A.at(id),
                                            "state_A_dihedral_list, " + oss.str() ) );
            EXPECT_TRUE( compare_2_vectors( atom[i].state_A_improper_list(),
                                            ref_improper_state_A.at(id),
                                            "state_A_improper_list, " + oss.str() ) );

            EXPECT_TRUE( compare_2_vectors( atom[i].state_B_angle_list(),
                                            ref_angle_state_B.at(id),
                                            "state_B_angle_list, " + oss.str() ) );
            EXPECT_TRUE( compare_2_vectors( atom[i].state_B_dihedral_list(),
                                            ref_dihedral_state_B.at(id),
                                            "state_B_dihedral_list, " + oss.str() ) );
            EXPECT_TRUE( compare_2_vectors( atom[i].state_B_improper_list(),
                                            ref_improper_state_B.at(id),
                                            "state_B_improper_list, " + oss.str() ) );
        }
    }

    virtual void SetUp(){
        const PS::S64 n = 9;
        init_FDPS_obj(n);

        if(PS::Comm::getRank() == 0){
            atom.setNumberOfParticleLocal(n);
            for(PS::S64 i=0; i<n; ++i){
                atom[i].initData(i);
            }

            //--- define test connection structure
            //
            //  0-1
            //    |
            //  2-3=(A)=4=(B)=5-6
            //                |
            //                7-8
            //
            //  =(A)= or =(B)= is special connection.

            atom[0].bond.add(1);

            atom[1].bond.add(0);
            atom[1].bond.add(3);

            atom[2].bond.add(3);

            atom[3].bond.add(1);
            atom[3].bond.add(2);
            atom[3].sp_pair_A.add(4);

            atom[4].sp_pair_A.add(3);
            atom[4].sp_pair_B.add(5);

            atom[5].bond.add(6);
            atom[5].bond.add(7);
            atom[5].sp_pair_B.add(4);

            atom[6].bond.add(5);

            atom[7].bond.add(5);
            atom[7].bond.add(8);

            atom[8].bond.add(7);
        }

        //--- reference result
        ref_mask_all.resize(n);
        ref_angle_state_A.resize(n);
        ref_dihedral_state_A.resize(n);
        ref_improper_state_A.resize(n);
        ref_angle_state_B.resize(n);
        ref_dihedral_state_B.resize(n);
        ref_improper_state_B.resize(n);

        //------ for mask
        auto mask1 = mask_param.at(0);
        auto mask2 = mask_param.at(1);
        auto mask3 = mask_param.at(2);

        ref_mask_all.at(0) = { mask1.setId(1),
                               mask2.setId(3),
                               mask3.setId(2),
                               mask3.setId(4) };

        ref_mask_all.at(1) = { mask1.setId(0),
                               mask1.setId(3),
                               mask2.setId(2),
                               mask2.setId(4),
                               mask3.setId(5) };

        ref_mask_all.at(2) = { mask1.setId(3),
                               mask2.setId(1),
                               mask2.setId(4),
                               mask3.setId(0),
                               mask3.setId(5) };

        ref_mask_all.at(3) = { mask1.setId(1),
                               mask1.setId(2),
                               mask1.setId(4),
                               mask2.setId(0),
                               mask2.setId(5),
                               mask3.setId(6),
                               mask3.setId(7) };

        ref_mask_all.at(4) = { mask1.setId(3),
                               mask1.setId(5),
                               mask2.setId(1),
                               mask2.setId(2),
                               mask2.setId(6),
                               mask2.setId(7),
                               mask3.setId(0),
                               mask3.setId(8) };

        ref_mask_all.at(5) = { mask1.setId(4),
                               mask1.setId(6),
                               mask1.setId(7),
                               mask2.setId(3),
                               mask2.setId(8),
                               mask3.setId(1),
                               mask3.setId(2) };

        ref_mask_all.at(6) = { mask1.setId(5),
                               mask2.setId(4),
                               mask2.setId(7),
                               mask3.setId(3),
                               mask3.setId(8) };

        ref_mask_all.at(7) = { mask1.setId(5),
                               mask1.setId(8),
                               mask2.setId(4),
                               mask2.setId(6),
                               mask3.setId(3) };

        ref_mask_all.at(8) = { mask1.setId(7),
                               mask2.setId(5),
                               mask3.setId(4),
                               mask3.setId(6) };

        //------ for angle
        //--------- state A
        ref_angle_state_A.at(0) = { MD_DEFS::AngleSet{0, 1, 3} };

        ref_angle_state_A.at(1) = { MD_DEFS::AngleSet{1, 3, 2},
                                    MD_DEFS::AngleSet{1, 3, 4},
                                    MD_DEFS::AngleSet{0, 1, 3} };

        ref_angle_state_A.at(2) = { MD_DEFS::AngleSet{2, 3, 1},
                                    MD_DEFS::AngleSet{2, 3, 4} };

        ref_angle_state_A.at(3) = { MD_DEFS::AngleSet{3, 1, 0},
                                    MD_DEFS::AngleSet{1, 3, 2},
                                    MD_DEFS::AngleSet{1, 3, 4},
                                    MD_DEFS::AngleSet{2, 3, 4} };

        ref_angle_state_A.at(4) = { MD_DEFS::AngleSet{4, 3, 1},
                                    MD_DEFS::AngleSet{4, 3, 2} };

        ref_angle_state_A.at(5) = { MD_DEFS::AngleSet{5, 7, 8},
                                    MD_DEFS::AngleSet{6, 5, 7} };

        ref_angle_state_A.at(6) = { MD_DEFS::AngleSet{6, 5, 7} };

        ref_angle_state_A.at(7) = { MD_DEFS::AngleSet{7, 5, 6},
                                    MD_DEFS::AngleSet{5, 7, 8} };

        ref_angle_state_A.at(8) = { MD_DEFS::AngleSet{8, 7, 5} };

        //--------- state B
        ref_angle_state_B.at(0) = { MD_DEFS::AngleSet{0, 1, 3} };

        ref_angle_state_B.at(1) = { MD_DEFS::AngleSet{1, 3, 2},
                                    MD_DEFS::AngleSet{0, 1, 3} };

        ref_angle_state_B.at(2) = { MD_DEFS::AngleSet{2, 3, 1} };

        ref_angle_state_B.at(3) = { MD_DEFS::AngleSet{3, 1, 0},
                                    MD_DEFS::AngleSet{1, 3, 2} };

        ref_angle_state_B.at(4) = { MD_DEFS::AngleSet{4, 5, 6},
                                    MD_DEFS::AngleSet{4, 5, 7} };

        ref_angle_state_B.at(5) = { MD_DEFS::AngleSet{5, 7, 8},
                                    MD_DEFS::AngleSet{6, 5, 7},
                                    MD_DEFS::AngleSet{6, 5, 4},
                                    MD_DEFS::AngleSet{7, 5, 4} };

        ref_angle_state_B.at(6) = { MD_DEFS::AngleSet{6, 5, 7},
                                    MD_DEFS::AngleSet{6, 5, 4} };

        ref_angle_state_B.at(7) = { MD_DEFS::AngleSet{7, 5, 6},
                                    MD_DEFS::AngleSet{7, 5, 4},
                                    MD_DEFS::AngleSet{5, 7, 8} };

        ref_angle_state_B.at(8) = { MD_DEFS::AngleSet{8, 7, 5} };

        //------ for dihedral torsion
        //--------- state A
        ref_dihedral_state_A.at(0) = { MD_DEFS::TorsionSet{0, 1, 3, 2},
                                       MD_DEFS::TorsionSet{0, 1, 3, 4} };

        ref_dihedral_state_A.at(1) = { MD_DEFS::TorsionSet{0, 1, 3, 2},
                                       MD_DEFS::TorsionSet{0, 1, 3, 4} };

        ref_dihedral_state_A.at(2) = { MD_DEFS::TorsionSet{2, 3, 1, 0} };

        ref_dihedral_state_A.at(3) = { MD_DEFS::TorsionSet{2, 3, 1, 0},
                                       MD_DEFS::TorsionSet{4, 3, 1, 0} };

        ref_dihedral_state_A.at(4) = { MD_DEFS::TorsionSet{4, 3, 1, 0} };

        ref_dihedral_state_A.at(5) = { MD_DEFS::TorsionSet{6, 5, 7, 8} };

        ref_dihedral_state_A.at(6) = { MD_DEFS::TorsionSet{6, 5, 7, 8} };

        ref_dihedral_state_A.at(7) = { MD_DEFS::TorsionSet{8, 7, 5, 6} };

        ref_dihedral_state_A.at(8) = { MD_DEFS::TorsionSet{8, 7, 5, 6} };

        //--------- state B
        ref_dihedral_state_B.at(0) = { MD_DEFS::TorsionSet{0, 1, 3, 2} };

        ref_dihedral_state_B.at(1) = { MD_DEFS::TorsionSet{0, 1, 3, 2} };

        ref_dihedral_state_B.at(2) = { MD_DEFS::TorsionSet{2, 3, 1, 0} };

        ref_dihedral_state_B.at(3) = { MD_DEFS::TorsionSet{2, 3, 1, 0} };

        ref_dihedral_state_B.at(4) = { MD_DEFS::TorsionSet{4, 5, 7, 8} };

        ref_dihedral_state_B.at(5) = { MD_DEFS::TorsionSet{6, 5, 7, 8},
                                       MD_DEFS::TorsionSet{4, 5, 7, 8} };

        ref_dihedral_state_B.at(6) = { MD_DEFS::TorsionSet{6, 5, 7, 8} };

        ref_dihedral_state_B.at(7) = { MD_DEFS::TorsionSet{8, 7, 5, 6},
                                       MD_DEFS::TorsionSet{8, 7, 5, 4} };

        ref_dihedral_state_B.at(8) = { MD_DEFS::TorsionSet{8, 7, 5, 6},
                                       MD_DEFS::TorsionSet{8, 7, 5, 4} };

        //------ for improper torsion
        //--------- state A
        ref_improper_state_A.at(0) = {};

        ref_improper_state_A.at(1) = { MD_DEFS::TorsionSet{2, 1, 3, 4},
                                       MD_DEFS::TorsionSet{1, 3, 2, 4},
                                       MD_DEFS::TorsionSet{1, 3, 4, 2} };

        ref_improper_state_A.at(2) = { MD_DEFS::TorsionSet{1, 2, 3, 4},
                                       MD_DEFS::TorsionSet{2, 3, 1, 4},
                                       MD_DEFS::TorsionSet{2, 3, 4, 1} };

        ref_improper_state_A.at(3) = { MD_DEFS::TorsionSet{2, 3, 1, 4},
                                       MD_DEFS::TorsionSet{1, 3, 2, 4},
                                       MD_DEFS::TorsionSet{1, 3, 4, 2} };

        ref_improper_state_A.at(4) = { MD_DEFS::TorsionSet{1, 4, 3, 2},
                                       MD_DEFS::TorsionSet{4, 3, 1, 2},
                                       MD_DEFS::TorsionSet{4, 3, 2, 1} };

        ref_improper_state_A.at(5) = {};
        ref_improper_state_A.at(6) = {};
        ref_improper_state_A.at(7) = {};
        ref_improper_state_A.at(8) = {};

        //--------- state B
        ref_improper_state_B.at(0) = {};
        ref_improper_state_B.at(1) = {};
        ref_improper_state_B.at(2) = {};
        ref_improper_state_B.at(3) = {};

        ref_improper_state_B.at(4) = { MD_DEFS::TorsionSet{6, 4, 5, 7},
                                       MD_DEFS::TorsionSet{4, 5, 6, 7},
                                       MD_DEFS::TorsionSet{4, 5, 7, 6} };

        ref_improper_state_B.at(5) = { MD_DEFS::TorsionSet{7, 5, 6, 4},
                                       MD_DEFS::TorsionSet{6, 5, 7, 4},
                                       MD_DEFS::TorsionSet{6, 5, 4, 7} };

        ref_improper_state_B.at(6) = { MD_DEFS::TorsionSet{7, 6, 5, 4},
                                       MD_DEFS::TorsionSet{6, 5, 7, 4},
                                       MD_DEFS::TorsionSet{6, 5, 4, 7} };

        ref_improper_state_B.at(7) = { MD_DEFS::TorsionSet{6, 7, 5, 4},
                                       MD_DEFS::TorsionSet{7, 5, 6, 4},
                                       MD_DEFS::TorsionSet{7, 5, 4, 6} };

        ref_improper_state_B.at(8) = {};
    }
};

TEST_F(IntraPairSPconnect, singleProc){
    make_intra_list();
    check_result();
}

TEST_F(IntraPairSPconnect, withMPI){
    //--- reduce the atom data for all processes.
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    make_intra_list();
    check_result();
}

#include "gtest_main_mpi.hpp"
