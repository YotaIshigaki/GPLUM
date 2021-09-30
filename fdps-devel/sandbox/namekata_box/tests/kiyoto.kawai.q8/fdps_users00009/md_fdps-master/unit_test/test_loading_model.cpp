//***************************************************************************************
//  This is unit test of loading model parameter.
//***************************************************************************************

#include <tuple>
#include <iostream>

#include <particle_simulator.hpp>

#include "md_enum.hpp"
#include "atom_class.hpp"
#include "md_coef_table.hpp"
#include "md_loading_condition.hpp"
#include "md_loading_model.hpp"

int main(int argc, char *argv[]) {

    PS::Initialize(argc, argv);

    if(PS::Comm::getRank() == 0) {
        //--- molecular model loading test
        std::cout << "TEST: molecular model loading test..." << std::endl;

        if(argc <= 1){
            std::cout << "    no input." << std::endl
                      << "    usage:  $ ./test_model.x [model_file_name]";
        }

        System::model_list.clear();
        for(int i=0; i<argc-1; i++){

            std::string model_name = argv[i+1];
            std::cout << std::endl;
            std::cout << "  model file: " << model_name << std::endl;

            MolName model = ENUM::which_MolName(model_name);
            System::model_list.push_back( std::make_pair(model, 0) );
            System::model_template.push_back( std::vector<Atom_FP>{} );

            MODEL::coef_table.clear();
            MODEL::loading_model_parameter(model_name,
                                           System::model_template.at(i),
                                           MODEL::coef_table            );

            //--- show result
            std::cout << endl;
            std::cout << "result: atom_list" << std::endl;
            std::cout << "  \"MolID\" was set as illigal value(-1). this is model template." << std::endl;
            MODEL::print_model_template( System::model_template.at(i));


            std::cout << "result: residue_table" << std::endl;
            MODEL::print_coef_table( MODEL::coef_table.residue,
                                     ENUM::which_MolName(model_name) );


            std::cout << "result: bond_table" << std::endl;
            MODEL::print_coef_table( MODEL::coef_table.bond,
                                     ENUM::which_MolName(model_name) );


            std::cout << "result: angle_table" << std::endl;
            MODEL::print_coef_table( MODEL::coef_table.angle,
                                     ENUM::which_MolName(model_name) );


            std::cout << "result: torsion_table" << std::endl;
            MODEL::print_coef_table( MODEL::coef_table.torsion,
                                     ENUM::which_MolName(model_name) );

            std::cout << "result: intra mask coefficient" << std::endl;
            const auto& mask_param = MODEL::coef_table.mask_scaling.at( ENUM::which_MolName(model_name) );
            int order = 0;
            for(const auto& mask : mask_param){
                ++order;
                std::cout << "   order: " << order << "\n"
                          << "      mask_LJ      = " << mask.scale_LJ      << "\n"
                          << "      mask_coulomb = " << mask.scale_coulomb << std::endl;
            }
            std::cout << "  " << mask_param.size() << " parameters were set." << "\n" << std::endl;
        }

        std::cout << "    the test succeeded." << std::endl;
    }

    PS::Finalize();
    return 0;
}
