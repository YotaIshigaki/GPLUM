#pragma once

std::vector<std::string> split_str(std::string str,
                                   std::vector<char> del)
{
    PS::S32 i = 0;
    PS::S32 n_del = del.size();
    std::vector<std::string> list;
    list.clear();
    
    do {
        std::string word;
        while (1){
            bool flag = false;
            for ( PS::S32 ii=0; ii<n_del; ii++ ) if ( str[i]==del.at(ii) ) flag = true;
            if ( str[i]=='\0' || !flag ) break;
            i++;
        }
        PS::S32 j = 0;
        PS::S32 i_start = i;
        while (1){
            bool flag = true;
            for ( PS::S32 ii=0; ii<n_del; ii++ ) {
                if ( str[i]==del.at(ii) ) flag = false;
                break;
            }
            if ( str[i]=='\0' || !flag ) break;
            i++; j++;
        }
        word = str.substr(i_start,j);
        list.push_back(word);
    } while ( str[i]!='\0' );
    
    return list;
}

bool isPowerOf2(PS::F64 p0)
{
    bool check = false;
    //PS::F64 p = std::abs(p0);
    PS::F64 p = p0; 
    
    if ( 0. < p && p < 1. ){
        while ( p < 1. ) p *= 2.;
        if ( p == 1. ) check = true;
    } else if (  1. < p ) {
        while ( 1. < p ) p *= 0.5;
        if ( p == 1. ) check = true;
    } else if ( p == 1. ){
        check = true;
    }

    return check;
}

PS::S32 readParameter(const char * param_file,
                      char * init_file,
                      bool & bHeader,
                      char * output_dir,
                      bool & bRestert,
                      bool & makeInit,
                      PS::F64 & coef_ema,
                      PS::S32 & nx,
                      PS::S32 & ny,
                      PS::F64 & theta,
                      PS::S32 & n_leaf_limit,
                      PS::S32 & n_group_limit,
                      PS::S32 & n_smp_ave,
                      PS::F64 & t_end,
                      PS::F64 & dt_snap,
                      PS::F64 & dt_snap_tmp,
                      PS::F64 & r_max,
                      PS::F64 & r_min,
                      PS::S32 & seed,
                      PS::S32 & reset_step)
{
    std::ifstream ifs(param_file);
    std::string line;
    std::string name;
    std::string value;
    std::vector<char> del{' ','\t','='};
    //std::vector<char>().swap(del);
    //del.push_back(' '); del.push_back('\t'); del.push_back('=');

    const PS::F64 L_MKS = 149597870700;
    const PS::F64 L_CGS = 14959787070000;
    const PS::F64 M_MKS = 1.9884e30;
    const PS::F64 M_CGS = 1.9884e33;
    const PS::F64 T     = 365.25*24.*60.*60./(2.*MY_PI);
    
    if ( ifs.fail() ) {
        errorMessage("The parameter file has FAILED to be successfully opened");
        return 1;
    }
    
    while ( getline(ifs, line) ){
        std::vector<std::string> list;
        if ( line[0]=='#' || line[0]=='\0' ){
            continue;
        } else {
            list = split_str(line, del);
            if ( list.size() < 2 ) continue;
        }

        name = list.at(0);
        value = list.at(1);
        if ( name == "seed" ){
            seed = std::atoi(value.c_str());
            
        } else if ( name == "init_file" ){
            sprintf(init_file,"%s",value.c_str());
            
        } else if ( name == "Header" ){
            bHeader = std::atoi(value.c_str()) > 0;
            
        } else if ( name == "output_dir" ){
            sprintf(output_dir,"%s",value.c_str());
            
        } else if ( name == "Restart" ){
            bRestert = std::atoi(value.c_str()) > 0;

        } else if ( name == "makeInit" ){
            makeInit = std::atoi(value.c_str()) > 0;

        } else if ( name == "n_init" ){
            SolidDisk::n_init = std::atoi(value.c_str());

        } else if ( name == "m_init" ){
            SolidDisk::m_init = getvalue(value, M_MKS, M_CGS);

        } else if ( name == "p" ){
            SolidDisk::p = getvalue(value, 1., 1.);

        } else if ( name == "f_dust" ){
            SolidDisk::f_dust = getvalue(value, 1., 1.);

        } else if ( name == "eta_ice" ){
            SolidDisk::eta_ice = getvalue(value, 1., 1.);

        } else if ( name == "a_in" ){
            SolidDisk::a_in = getvalue(value, L_MKS, L_CGS);

        } else if ( name == "a_out" ){
            SolidDisk::a_out = getvalue(value, L_MKS, L_CGS);

        } else if ( name == "a_ice" ){
            SolidDisk::a_ice = getvalue(value, L_MKS, L_CGS);

        } else if ( name == "ecc_hill" ){
            SolidDisk::ecc_hill = getvalue(value, 1., 1.);
            
        } else if ( name == "inc_hill" ){
            SolidDisk::inc_hill = getvalue(value, 1., 1.);

        } else if ( name == "alpha_gas" ){
            GasDisk::alpha_gas = getvalue(value, 1., 1.);

        } else if ( name == "beta_gas" ){
            GasDisk::beta_gas = getvalue(value, 1., 1.);

        } else if ( name == "f_gas" ){
            GasDisk::f_gas = getvalue(value, 1., 1.);

        } else if ( name == "tau_gas" ){
            GasDisk::tau_gas = getvalue(value, T, T);

        } else if ( name == "C_d" ){
            GasDisk::C_d = getvalue(value, 1., 1.);

        } else if ( name == "mu" ){
            GasDisk::mu = getvalue(value, 1., 1.);

        } else if ( name == "coef_ema" ){
            coef_ema = getvalue(value, 1., 1.);
            
        } else if ( name == "nx" ){
            nx = std::atoi(value.c_str());
            
        } else if ( name == "ny" ){
            ny = std::atoi(value.c_str());
            
        } else if ( name == "theta" ){
            theta = getvalue(value, 1., 1.);
            
        } else if ( name == "n_leaf_limit" ){
            n_leaf_limit = std::atoi(value.c_str());
            
        } else if ( name == "n_group_limit" ){
            n_group_limit = std::atoi(value.c_str());

        } else if ( name == "n_smp_ave" ){
            n_smp_ave = std::atoi(value.c_str());
       
        } else if ( name == "t_end" ){
            t_end = getvalue(value, T, T);
            
        } else if ( name == "dt_snap" ){
            dt_snap = getvalue(value, T, T);

        } else if ( name == "dt_snap_tmp" ){
            dt_snap_tmp = getvalue(value, T, T);
            
        } else if ( name == "dt_tree" ){
            FP_t::dt_tree = getvalue(value, T, T);
            
        } else if ( name == "dt_min" ){
            FP_t::dt_min = getvalue(value, T, T);
            
        } else if ( name == "eta" ){
            FP_t::eta = getvalue(value, 1., 1.);
            
        } else if ( name == "eta_0" ){
            FP_t::eta_0 = getvalue(value, 1., 1.);
            
        } else if ( name == "eta_sun" ){
            FP_t::eta_sun = getvalue(value, 1., 1.);
            
        } else if ( name == "eta_sun0" ){
            FP_t::eta_sun0 = getvalue(value, 1., 1.);
            
        } else if ( name == "alpha" ){
            PS::F64 alpha = getvalue(value, 1., 1.);
            FP_t::alpha2 = alpha*alpha;
            
        } else if ( name == "m_sun" ){
            FP_t::m_sun = getvalue(value, M_MKS, M_CGS);

        } else if ( name == "dens" ){
            FP_t::dens = getvalue(value, M_MKS/(L_MKS*L_MKS*L_MKS), M_CGS/(L_CGS*L_CGS*L_CGS));
            
        } else if ( name == "eps" ){
            PS::F64 eps = getvalue(value, L_MKS, L_CGS);
            FP_t::eps2 = eps*eps;

        } else if ( name == "eps_sun" ){
            PS::F64 eps_sun = getvalue(value, L_MKS, L_CGS);
            FP_t::eps2_sun = eps_sun*eps_sun;
            
        } else if ( name == "R_cut0" ){
            FP_t::R_cut0 = getvalue(value, 1., 1.);

        } else if ( name == "R_cut1" ){
            FP_t::R_cut1 = getvalue(value, 1., 1.);
            
        } else if ( name == "R_search0" ){
            FP_t::R_search0 = getvalue(value, 1., 1.);
            
        } else if ( name == "R_search1" ){
            FP_t::R_search1 = getvalue(value, 1., 1.);
            
#ifdef USE_RE_SEARCH_NEIGHBOR
        } else if ( name == "R_search2" ){
            FP_t::R_search2 = getvalue(value, 1., 1.);
            
        } else if ( name == "R_search3" ){
            FP_t::R_search3 = getvalue(value, 1., 1.);
#endif
#ifdef MERGE_BINARY
        } else if ( name == "R_merge" ){
            FP_t::R_merge = getvalue(value, 1., 1.);
#endif
#ifdef CONSTANT_RANDOM_VELOCITY
        } else if ( name == "v_disp" ){
            FP_t::v_disp = getvalue(value, L_MKS/T, L_CGS/T);
#endif
        } else if ( name == "gamma" ){
            FP_t::setGamma(getvalue(value, 1., 1.));

        } else if ( name == "r_cut_min" ){
            FP_t::r_cut_min = getvalue(value, L_MKS, L_CGS);
            
        } else if ( name == "r_cut_max" ){
            FP_t::r_cut_max = getvalue(value, L_MKS, L_CGS);

        } else if ( name == "p_cut" ){
            FP_t::p_cut = getvalue(value, 1., 1.);
            
        } else if ( name == "r_max" ){
            r_max = getvalue(value, L_MKS, L_CGS);
            
        } else if ( name == "r_min" ){
            r_min = getvalue(value, L_MKS, L_CGS);
            
        } else if ( name == "f" ){
            //HardSystem::f = getvalue(value, 1., 1.);
            //Collision0::f = getvalue(value, 1., 1.);
            FP_t::increase_factor = getvalue(value, 1., 1.);

        } else if ( name == "m_min" ){
           Collision0::m_min = getvalue(value, M_MKS, M_CGS);
            
        } else if ( name == "seed" ){
            seed = std::atoi(value.c_str());
        } else if ( name == "reset_step" ){
            reset_step = std::atoi(value.c_str());
        }
        Collision::readParameter(name, value);
    }
    ifs.close();

#ifdef TEST_PTCL
    FP_t::eps2 = std::max(FP_t::eps2, (PS::F64)std::numeric_limits<PS::F32>::min());
#endif

    return 0;
}


PS::S32 checkParameter(char * init_file,
                       bool bHeader,
                       char * output_dir,
                       bool bRestert,
                       bool makeInit,
                       PS::F64 coef_ema,
                       PS::S32 nx,
                       PS::S32 ny,
                       PS::F64 theta,
                       PS::S32 n_leaf_limit,
                       PS::S32 n_group_limit,
                       PS::S32 n_smp_ave,
                       PS::F64 t_end,
                       PS::F64 dt_snap,
                       PS::F64 dt_snap_tmp,
                       PS::F64 r_max,
                       PS::F64 r_min,
                       PS::S32 seed,
                       PS::S32 reset_step)
{
    if ( !isPowerOf2(FP_t::dt_tree) ){
        errorMessage("dt_tree has NOT been set to be power of 2",
                     "(dt_tree =" + std::to_string(FP_t::dt_tree ) + ")." );
        return 1;
    } else if ( !isPowerOf2(FP_t::dt_min) ){
        errorMessage("dt_min has NOT been set to be power of 2",
                     "(dt_min = " + std::to_string(FP_t::dt_min) + ").");
        return 1;
    } else if ( FP_t::dt_min > 0.5*FP_t::dt_tree ){
        errorMessage("dt_min has NOT been set to be satisfy dt_min <= 0.5 dt_tree",
                     "(dt_min = " + std::to_string(FP_t::dt_min)
                     + ", dt_tree = " + std::to_string(FP_t::dt_tree) + ").");
        return 1;
    } else if ( fmod(dt_snap, FP_t::dt_tree) != 0 ){
        errorMessage("dt_snap has NOT been set to be multiples of dt_tree",
                     "(dt_snap = " + std::to_string(dt_snap)
                     + ", dt_tree = " + std::to_string(FP_t::dt_tree) + ").");
        return 1;

    } else if ( fmod(dt_snap_tmp, FP_t::dt_tree) != 0 ){
        errorMessage("dt_snap_tmp has NOT been set to be multiples of dt_tree",
                     "(dt_snap_tmp = " + std::to_string(dt_snap)
                     + ", dt_tree = " + std::to_string(FP_t::dt_tree) + ").");
        return 1;
    } else if ( FP_t::dt_tree > FP_t::dt_min * pow2(51) ){
        errorMessage("dt_min has NOT been set to satisfy dt_tree * 2^-51 <= dt_min",
                     "(dt_min = " + std::to_string(FP_t::dt_min)
                     + ", dt_tree = " + std::to_string(FP_t::dt_tree) + ").");
        return 1;
    } else if ( makeInit && SolidDisk::n_init == 0 && SolidDisk::m_init == 0. ){
        errorMessage("Both n_init & m_init have NOT been set.");
        return 1;
    } else if ( FP_t::r_cut_min > FP_t::r_cut_max && FP_t::r_cut_max > 0 ){
        errorMessage("r_cut_max & r_cut_min have NOT been set to satisfy r_cut_max >= r_cut_min",
                     "(r_cut_max = " + std::to_string(FP_t::r_cut_max)
                     + ", r_cut_min = " + std::to_string(FP_t::r_cut_min) + ").");
        return 1;
    } else if ( FP_t::R_cut0 < 0. ){
        errorMessage("R_cut0 has NOT been set to satisfy R_cut0 >= 0",
                     "(R_cut0 = " + std::to_string(FP_t::R_cut0) + ").");
        return 1;
    } else if ( FP_t::R_cut1 < 0. ){
        errorMessage("R_cut1 has NOT been set to satisfy R_cut1 >= 0",
                     "(R_cut1 = " + std::to_string(FP_t::R_cut1) + ").");
        return 1;
    } else if ( FP_t::R_search0 < 1. ){
        errorMessage("R_search0 has NOT been set to satisfy R_search0 >= 1",
                     "(R_search0 = " + std::to_string(FP_t::R_search0) + ").");
        return 1;
    } else if ( FP_t::R_search1 < 0. ){
        errorMessage("R_search1 has NOT been set to satisfy R_search1 >= 0",
                     "(R_search1 = " + std::to_string(FP_t::R_search1) + ").");
        return 1;
#ifdef USE_RE_SEARCH_NEIGHBOR
    } else if ( FP_t::R_search2 < 1. ){
        errorMessage("R_search2 has NOT been set to satisfy R_search2 >= 1",
                     "(R_search2 = " + std::to_string(FP_t::R_search2) + ").");
        return 1;
    } else if ( FP_t::R_search3 < 0. ){
        errorMessage("R_search3 has NOT been set to satisfy R_search0 >= 0",
                     "(R_search3 = " + std::to_string(FP_t::R_search3) + ").");
        return 1;
#endif
    } else if ( r_max < r_min ){
        errorMessage("r_max & r_min have NOT been set to satisfy r_max >= r_min",
                     "(r_max = "+ std::to_string(r_max)
                     + ", r_min = " + std::to_string(r_min) + ").");
        return 1;
    } else if ( nx * ny != PS::Comm::getNumberOfProc() ) {
        errorMessage("nx & ny have NOT been set to satisfy nx * ny have been equal to the number of process",
                     "(nx = "+ std::to_string(nx)
                     + ", ny = " + std::to_string(ny) + ").");
        return 1;
    }

    return 0;
}

void showParameter(char * init_file,
                   char * dir_name,
                   bool makeInit,
                   PS::F64 time_sys,
                   PS::F64 coef_ema,
                   PS::S32 nx,
                   PS::S32 ny,
                   PS::F64 theta,
                   PS::S32 n_leaf_limit,
                   PS::S32 n_group_limit,
                   PS::S32 n_smp_ave,
                   PS::F64 t_end,
                   PS::F64 dt_snap,
                   PS::F64 dt_snap_tmp,
                   PS::F64 r_max,
                   PS::F64 r_min,
                   PS::S32 seed,
                   PS::S32 reset_step)
{
    const PS::F64 L = 14959787070000;
    const PS::F64 M = 1.9884e33;
    const PS::F64 T = 365.25*24.*60.*60./(2.*MY_PI);
    
    if ( PS::Comm::getRank() == 0 ){
        std::cout << "Number Of Processes:\t" << PS::Comm::getNumberOfProc() << std::endl;
        std::cout << "Number Of Threads Per Process:\t" << PS::Comm::getNumberOfThread() << std::endl;
        std::cout << std::endl;

        std::cout << "#################### Parameters ####################" << std::endl
            //std::cout << std::setprecision(15);
                  << "Initial File:\t" << init_file << std::endl
                  << "Output Directory:\t" << dir_name << std::endl
                  << "seed          = " << seed << std::endl;
        if ( makeInit ){
            std::cout << std::scientific << std::setprecision(15)
                      << "n_init        = " << SolidDisk::n_init << std::endl
                      << "m_init        = " << SolidDisk::m_init << "\t(" << SolidDisk::m_init*M << " g)" << std::endl
                      << "p             = " << SolidDisk::p << std::endl
                      << "f_dust        = " << SolidDisk::f_dust << std::endl
                      << "eta_ice       = " << SolidDisk::eta_ice << std::endl
                      << "a_in          = " << SolidDisk::a_in << std::endl
                      << "a_out         = " << SolidDisk::a_out << std::endl
                      << "a_ice         = " << SolidDisk::a_ice << std::endl
                      << "ecc_hill      = " << SolidDisk::ecc_hill << std::endl
                      << "inc_hill      = " << SolidDisk::inc_hill << std::endl;
        }
#ifdef GAS_DRAG
        std::cout << std::scientific << std::setprecision(15)
                  << "alpha_gas     = " << GasDisk::alpha_gas << std::endl
                  << "beta_gas      = " << GasDisk::beta_gas << std::endl
                  << "f_gas         = " << GasDisk::f_gas << std::endl
                  << "tau_gas       = " << GasDisk::tau_gas << "\t(" << GasDisk::tau_gas/(2.*MY_PI) << " year)" << std::endl
                  << "C_d           = " << GasDisk::C_d << std::endl
                  << "mu            = " << GasDisk::mu << std::endl;
#endif
        std::cout << std::fixed << std::setprecision(5)
                  << "coef_ema      = " << coef_ema << std::endl
                  << "nx, ny        = " << nx << ", " << ny << std::endl
                  << "reset_step    = " << reset_step << std::endl
                  << "theta         = " << theta << std::endl
                  << "n_leaf_limit  = " << n_leaf_limit << std::endl
                  << "n_group_limit = " << n_group_limit << std::endl
                  << "n_smp_ave     = " << n_smp_ave << std::endl
                  << std::scientific << std::setprecision(15)
                  << "t_begin       = " << time_sys << "\t(" << time_sys/(2.*MY_PI) << " year)" << std::endl
                  << "t_end         = " << t_end << "\t(" << t_end/(2.*MY_PI) << " year)" << std::endl
                  << "dt_snap       = " << dt_snap << "\t(" << dt_snap/(2.*MY_PI) << " year)" << std::endl
                  << "dt_snap_tmp   = " << dt_snap_tmp << "\t(" << dt_snap_tmp/(2.*MY_PI) << " year)" << std::endl
                  << "dt_tree       = " << FP_t::dt_tree << "\t(2^" << (PS::S32)std::log2(FP_t::dt_tree) << ", " << FP_t::dt_tree/(2.*MY_PI) << " year)" << std::endl
                  << "dt_min        = " << FP_t::dt_min << "\t(2^" << (PS::S32)std::log2(FP_t::dt_min) << ", " << FP_t::dt_min/(2.*MY_PI) << " year)" << std::endl
                  << "eta           = " << FP_t::eta  << std::endl
                  << "eta_0         = " << FP_t::eta_0 << std::endl
                  << "eta_sun       = " << FP_t::eta_sun << std::endl
                  << "eta_sun0      = " << FP_t::eta_sun0 << std::endl
                  << "alpha         = " << sqrt(FP_t::alpha2) << std::endl
                  << "m_sun         = " << FP_t::m_sun << "\t(" << FP_t::m_sun*M << " g)" << std::endl
                  << "dens          = " << FP_t::dens << "\t(" << FP_t::dens*M/(L*L*L) << " g/cm^3)"<< std::endl
                  << "eps           = " << sqrt(FP_t::eps2) << "\t(" << sqrt(FP_t::eps2)*L << " cm)"<< std::endl
                  << "eps_sun       = " << sqrt(FP_t::eps2_sun) << "\t(" << sqrt(FP_t::eps2_sun)*L << " cm)"<< std::endl
                  << std::fixed << std::setprecision(5)
                  << "R_cut0        = " << FP_t::R_cut0 << std::endl
                  << "R_cut1        = " << FP_t::R_cut1 << std::endl
                  << "R_search0     = " << FP_t::R_search0 << std::endl
                  << "R_search1     = " << FP_t::R_search1 << std::endl
#ifdef USE_RE_SEARCH_NEIGHBOR
                  << "R_search2     = " << FP_t::R_search2 << std::endl
                  << "R_search3     = " << FP_t::R_search3 << std::endl
#endif
#ifdef MERGE_BINARY
                  << "R_merge       = " << FP_t::R_merge << std::endl
#endif
#ifdef CONSTANT_RANDOM_VELOCITY
                  << "v_disp        = " << FP_t::v_disp << "\t(" << FP_t::v_disp*L/T << " cm/s)"<< std::endl
#endif
                  << "gamma         = " << FP_t::gamma << std::endl
                  << std::scientific << std::setprecision(15)
                  << "r_cut_max     = " << FP_t::r_cut_max << "\t(" << FP_t::r_cut_max*L << " cm)"<< std::endl
                  << "r_cut_min     = " << FP_t::r_cut_min << "\t(" << FP_t::r_cut_min*L << " cm)"<< std::endl
                  << "p_cut         = " << FP_t::p_cut << std::endl
                  << std::fixed << std::setprecision(5)
                  << "r_max         = " << r_max << std::endl
                  << "r_min         = " << r_min << std::endl
                  << "f             = " << FP_t::increase_factor << std::endl
                  << std::scientific << std::setprecision(15)
                  << "m_min         = " << Collision0::m_min << "\t(" << Collision0::m_min*M << " g)" << std::endl;

        Collision::showParameter();
        
        std::cout << std::resetiosflags(std::ios_base::floatfield)<<std::setprecision(8)
                  << "####################################################" << std::endl;

        char sout_param[256];
        std::ofstream fout_param;
        sprintf(sout_param,"%s/param.dat", dir_name);
        if ( time_sys == 0. ) {
            fout_param.open(sout_param, std::ios::out);
        } else {
            fout_param.open(sout_param, std::ios::app);
        }
        fout_param << "####################################################" << std::endl;

#ifdef USE_QUAD
        fout_param << "Use Quadrupole For Tree" << std::endl;
#else
        fout_param << "Use Monopole For Tree" << std::endl;
#endif

#ifdef __AVX512DQ__
        fout_param << "Use AVX512DQ" << std::endl;
#elif defined(__AVX2__)
        fout_param << "Use AVX2" << std::endl;
#endif
        
#ifdef USE_INDIVIDUAL_CUTOFF
        fout_param << "Use Individual CutOff" << std::endl;
#else
        fout_param << "Use Shared CutOff Radii & Search Radii" << std::endl;
#endif
#ifdef GAS_DRAG
        fout_param << "Use Gas Drag" << std::endl;
#endif
#if defined(KOMINAMI)
        fout_param << "Use Kominami Model" << std::endl;
#elif defined(CHAMBERS)
        fout_param << "Use Chambers Model" << std::endl;
#endif
#ifdef ISOTROPIC
        fout_param << "Use Isotropic Method to Set Random Velocity" << std::endl;
#endif 
        fout_param << std::endl;

        fout_param << "Number Of Processes:\t" << PS::Comm::getNumberOfProc() << std::endl;
        fout_param << "Number Of Threads Per Process:\t" << PS::Comm::getNumberOfThread() << std::endl;
        fout_param << std::endl;

        fout_param << "Initial File:\t" << init_file << std::endl
                   << "Output Directory:\t" << dir_name << std::endl
                   << "seed          = " << seed << std::endl;
        if ( makeInit ){
            fout_param << std::scientific << std::setprecision(15)
                       << "n_init        = " << SolidDisk::n_init << std::endl
                       << "m_init        = " << SolidDisk::m_init << "\t(" << SolidDisk::m_init*M << " g)" << std::endl
                       << "p             = " << SolidDisk::p << std::endl
                       << "f_dust        = " << SolidDisk::f_dust << std::endl
                       << "eta_ice       = " << SolidDisk::eta_ice << std::endl
                       << "a_in          = " << SolidDisk::a_in << std::endl
                       << "a_out         = " << SolidDisk::a_out << std::endl
                       << "a_ice         = " << SolidDisk::a_ice << std::endl
                       << "ecc_hill      = " << SolidDisk::ecc_hill << std::endl
                       << "inc_hill      = " << SolidDisk::inc_hill << std::endl;
        }
        #ifdef GAS_DRAG
        fout_param << std::scientific << std::setprecision(15)
                   << "alpha_gas     = " << GasDisk::alpha_gas << std::endl
                   << "beta_gas      = " << GasDisk::beta_gas << std::endl
                   << "f_gas         = " << GasDisk::f_gas << std::endl
                   << "tau_gas       = " << GasDisk::tau_gas << "\t(" << GasDisk::tau_gas/(2.*MY_PI) << " year)" << std::endl
                   << "C_d           = " << GasDisk::C_d << std::endl
                   << "mu            = " << GasDisk::mu << std::endl;
#endif
        fout_param << std::fixed << std::setprecision(5)
                   << "coef_ema      = " << coef_ema << std::endl
                   << "nx, ny        = " << nx << ", " << ny << std::endl
                   << "reset_step    = " << reset_step << std::endl
                   << "theta         = " << theta << std::endl
                   << "n_leaf_limit  = " << n_leaf_limit << std::endl
                   << "n_group_limit = " << n_group_limit << std::endl
                   << "n_smp_ave     = " << n_smp_ave << std::endl
                   << std::scientific << std::setprecision(15)
                   << "t_begin       = " << time_sys << "\t(" << time_sys/(2.*MY_PI) << " year)" << std::endl
                   << "t_end         = " << t_end << "\t(" << t_end/(2.*MY_PI) << " year)" << std::endl
                   << "dt_snap       = " << dt_snap << "\t(" << dt_snap/(2.*MY_PI) << " year)" << std::endl
                   << "dt_snap_tmp   = " << dt_snap_tmp << "\t(" << dt_snap_tmp/(2.*MY_PI) << " year)" << std::endl
                   << "dt_tree       = " << FP_t::dt_tree << "\t(2^" << (PS::S32)std::log2(FP_t::dt_tree) << ", " << FP_t::dt_tree/(2.*MY_PI) << " year)" << std::endl
                   << "dt_min        = " << FP_t::dt_min << "\t(2^" << (PS::S32)std::log2(FP_t::dt_min) << ", " << FP_t::dt_min/(2.*MY_PI) << " year)" << std::endl
                   << "eta           = " << FP_t::eta << std::endl
                   << "eta_0         = " << FP_t::eta_0 << std::endl
                   << "eta_sun       = " << FP_t::eta_sun << std::endl
                   << "eta_sun0      = " << FP_t::eta_sun0 << std::endl
                   << "alpha         = " << sqrt(FP_t::alpha2) << std::endl
                   << "m_sun         = " << FP_t::m_sun << "\t(" << FP_t::m_sun*M << " g)" << std::endl
                   << "dens          = " << FP_t::dens << "\t(" << FP_t::dens*M/(L*L*L) << " g/cm^3)"<< std::endl
                   << "eps           = " << sqrt(FP_t::eps2) << "\t(" << sqrt(FP_t::eps2)*L << " cm)"<< std::endl
                   << "eps_sun       = " << sqrt(FP_t::eps2_sun) << "\t(" << sqrt(FP_t::eps2_sun)*L << " cm)"<< std::endl
                   << std::fixed << std::setprecision(5)
                   << "R_cut0        = " << FP_t::R_cut0 << std::endl
                   << "R_cut1        = " << FP_t::R_cut1 << std::endl
                   << "R_search0     = " << FP_t::R_search0 << std::endl
                   << "R_search1     = " << FP_t::R_search1 << std::endl
#ifdef USE_RE_SEARCH_NEIGHBOR
                   << "R_search2     = " << FP_t::R_search2 << std::endl
                   << "R_search3     = " << FP_t::R_search3 << std::endl
#endif
#ifdef MERGE_BINARY
                   << "R_merge       = " << FP_t::R_merge << std::endl
#endif
#ifdef CONSTANT_RANDOM_VELOCITY
                   << "v_disp        = " << FP_t::v_disp << "\t(" << FP_t::v_disp*L/T << " cm/s)"<< std::endl
#endif
                   << "gamma         = " << FP_t::gamma << std::endl
                   << std::scientific << std::setprecision(15)
                   << "r_cut_max     = " << FP_t::r_cut_max << "\t(" << FP_t::r_cut_max*L << " cm)"<< std::endl
                   << "r_cut_min     = " << FP_t::r_cut_min << "\t(" << FP_t::r_cut_min*L << " cm)"<< std::endl
                   << "p_cut         = " << FP_t::p_cut << std::endl
                   << std::fixed << std::setprecision(5)
                   << "r_max         = " << r_max << std::endl
                   << "r_min         = " << r_min << std::endl
                   << "f             = " << FP_t::increase_factor << std::endl
                   << std::scientific << std::setprecision(15)
                   << "m_min         = " << Collision0::m_min << "\t(" << Collision0::m_min*M << " g)" << std::endl;

        Collision::showParameter(fout_param);
        fout_param.close();
    }
}

PS::S32 makeOutputDirectory(char * dir_name)
{
    struct stat st;
    if ( stat(dir_name, &st) != 0 ) {
        PS::S32 ret_loc = 0;
        PS::S32 ret     = 0;
        if ( PS::Comm::getRank() == 0 ) ret_loc = mkdir(dir_name, 0777);
        PS::Comm::broadcast(&ret_loc, ret);
        if (ret == 0) {
            if ( PS::Comm::getRank() == 0 ) {
                successMessage("The directory \"" + (std::string)dir_name  + "\" is successfully made.");
            }
        } else {
            if ( PS::Comm::getRank() == 0 ) {
                errorMessage("The directory \"" + (std::string)dir_name  + "\" has FAILED to be successfully made.");
                return 1;
            }
        }
    } else {
        if(PS::Comm::getRank() == 0){
            successMessage("The directory \"" + (std::string)dir_name  + "\" have already existed.");
        }
    }

    return 0;
}

PS::S32 getLastSnap(char * dir_name,
                    char * lastsnap_name)
{
    DIR * dir;
    dirent * entry;
    dir = opendir(dir_name);
    if ( dir == NULL ) {
        if(PS::Comm::getRank() == 0){
            errorMessage("The directory \"" + (std::string)dir_name  + "\" does NOT exist.");
        }
        return 1;
    }

    PS::S32 lastnumber = -1;
    do {
        entry = readdir(dir);
        if (entry != NULL){
            char  filename[64];
            char  head[16];
            strcpy(filename, entry->d_name);
            strncpy(head, filename, 4);
            if ( strcmp(head,"snap") == 0 ){
                char  number[16];
                strncpy(number, filename+4, 6);
		lastnumber = std::max(lastnumber, std::atoi(number));
                if ( lastnumber == std::atoi(number) ) {
		  //lastnumber = std::atoi(number);
                    sprintf(lastsnap_name, "%s/%s", dir_name, filename);
                }
            }
        }
    } while (entry != NULL);
    if ( lastnumber < 0 ) {
        if(PS::Comm::getRank() == 0) {
            errorMessage("NO snapshot file exists in the directory \"" + (std::string)dir_name  + "\".");
        }
        return 1;
    }
    
    return 0;
}
