//=======================================================================================
//  This is common routine for Force test.
//=======================================================================================

#include <cassert>


::testing::AssertionResult float_relative_eq(const double &lhs,
                                             const double &rhs,
                                             const double &abs_err,
                                             const double &relative_err){

    assert(abs_err      >= 0.0);
    assert(relative_err >= 0.0);

    std::ostringstream oss;
    oss << "\n"
        << "  lhs = " << lhs << "\n"
        << "  rhs = " << rhs << "\n";

    const double diff = lhs - rhs;
    if(std::abs(diff) > abs_err){
        oss << "    diff = " << diff << " > " << "absolute error = " << abs_err << "\n";
    } else {
        return ::testing::AssertionSuccess();
    }

    double r_diff = 0.0;
    if(lhs == 0.0){
        r_diff = diff/rhs;
    } else {
        r_diff = diff/lhs;
    }

    if(std::abs(r_diff) > relative_err ){
        oss << "    relative diff =" << r_diff << " > " << "relative error = " << relative_err << "\n";
    } else {
        return ::testing::AssertionSuccess();
    }

    return ::testing::AssertionFailure() << oss.str();
}

::testing::AssertionResult float_relative_eq(const double &lhs,
                                             const double &rhs,
                                             const double &relative_err){
    return float_relative_eq(lhs, rhs, relative_err, relative_err);
}

template <class Tlog>
void write_log_file(const Tlog &data_list, const std::string &file_name){
    if(PS::Comm::getRank() != 0) return;

    std::ofstream file{file_name};
    for(const auto& data : data_list){
        file << data;
    }
    file.close();
}

template <class Tlog>
void load_log_file(Tlog &data_list, const std::string &file_name){
    data_list.clear();

    if(PS::Comm::getRank() == 0){
        std::ifstream file_ref{file_name};
        if(file_ref.fail()){
            throw std::ios_base::failure("reference data: " + file_name + " was not found.");
        }

        std::string line;
        std::vector<std::string> str_list;
        while( getline(file_ref, line) ){
            STR_TOOL::removeCR(line);
            str_list = STR_TOOL::split(line, " ");
            read_ref_data(str_list, data_list);
        }
    }
    COMM_TOOL::broadcast(data_list, 0);
}

void test_init(const PS::S64 n_step){
    if(PS::Comm::getRank() != 0) return;

    //--- system parameters
    System::profile.coef_ema      = 0.3;
    System::profile.theta         = 0.5;
    System::profile.n_leaf_limit  = 8;
    System::profile.n_group_limit = 64;

    System::profile.cut_off_LJ    = 12.0;
    System::profile.cut_off_intra =  9.0;

    //--- set loop condition
    System::profile.dt       = 1.0;
    System::profile.istep    = 0;
    System::profile.nstep_st = 0;
    System::profile.nstep_ed = n_step;

    //--- set domain size
    Normalize::setBoxSize( PS::F32vec{ 40.0,
                                       40.0,
                                       40.0 } );
}

template <class Tptcl, class Tdinfo, class Tforce,
          class Tdata>
void execute_force_calc(Tptcl              &atom,
                        Tdinfo             &dinfo,
                        Tforce             &force,
                        std::vector<Tdata> &force_log,
                        std::vector<Tdata> &force_ref ){

    //--- sync settings
    System::broadcast_profile(0);
    MODEL::coef_table.broadcast(0);
    System::InitDinfo(dinfo);

    //--- split domain & particle
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    //--- initialize force calculator
    const PS::S64 n_total = atom.getNumberOfParticleGlobal();
    force.init(n_total);

    //--- calculate force
    PS::S32 record_count = 0;
    while( System::isLoopContinue() ){
        //--- exchange particle
        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        //--- calculate intermolecular force in FDPS
        force.update_intra_pair_list(atom, dinfo, MODEL::coef_table.mask_scaling);
        force.update_force(atom, dinfo);

        //--- recording
        test_record(atom, record_count, force_log);
        ++record_count;

        //--- move
        test_move(atom);

        //--- nest step
        System::StepNext();
    }
}
