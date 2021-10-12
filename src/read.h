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

class Parameter{
public:
    static char param_file[256];
    
    static char init_file[256];
    static char output_dir[256];
    bool bHeader;
    bool bRestart;
    
    bool makeInit;
    
    PS::F64 coef_ema;
    PS::S32 nx;
    PS::S32 ny;
    
    PS::F64 theta;
    PS::S32 n_leaf_limit;
    PS::S32 n_group_limit;
    PS::S32 n_smp_ave;
    
    PS::F64 t_end;
    PS::F64 dt_snap;
    PS::F64 dt_snap_tmp;
    
    PS::F64 r_max;
    PS::F64 r_min;
    
    PS::S32 seed;
    PS::S32 reset_step;

    PS::F64 wtime_max;

    Parameter(){
        sprintf(param_file, "parameter.dat");
    
        sprintf(init_file, "INIT_3000.dat");
        sprintf(output_dir, "OUTPUT");
        bHeader  = false;
        bRestart = false;

        makeInit = false;

        coef_ema = 0.3;
        nx = 0;
        ny = 0;
    
        theta         = 0.5;
        n_leaf_limit  = 8;
        n_group_limit = 256;
        n_smp_ave     = 100;
    
        t_end       = 1.;
        dt_snap     = pow2(-5);
        dt_snap_tmp = pow2(-5);

        r_max = 40.;
        r_min = 0.1;

        seed = 1;
        reset_step = 1024;

        wtime_max = 0.;
    }

    PS::S32 readParameter();
    PS::S32 checkParameter();
    void showParameter(PS::F64 time_sys,
                       std::ostream & fout);
};

char Parameter::param_file[256];
char Parameter::init_file[256];
char Parameter::output_dir[256];

PS::S32 Parameter::readParameter()
{
    std::ifstream ifs(param_file);
    std::string line;
    std::string name;
    std::string value;
    std::vector<char> del{' ','\t','='};
    //std::vector<char>().swap(del);
    //del.push_back(' '); del.push_back('\t'); del.push_back('=');
    
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
        if ( name == "seed" ) seed = std::atoi(value.c_str());
        if ( name == "init_file" ) sprintf(init_file,"%s",value.c_str());
        if ( name == "Header" ) bHeader = std::atoi(value.c_str()) > 0;
        if ( name == "output_dir" ) sprintf(output_dir,"%s",value.c_str());
        if ( name == "Restart" ) bRestart = std::atoi(value.c_str()) > 0;
        if ( name == "makeInit" ) makeInit = std::atoi(value.c_str()) > 0;
        if ( name == "coef_ema" ) coef_ema = getvalue(value, 1., 1.);
        if ( name == "nx" ) nx = std::atoi(value.c_str());
        if ( name == "ny" ) ny = std::atoi(value.c_str());
        if ( name == "theta" ) theta = getvalue(value, 1., 1.);
        if ( name == "n_leaf_limit" ) n_leaf_limit = std::atoi(value.c_str());
        if ( name == "n_group_limit" ) n_group_limit = std::atoi(value.c_str());
        if ( name == "n_smp_ave" ) n_smp_ave = std::atoi(value.c_str());
        if ( name == "t_end" ) t_end = getvalue(value, T_MKS, T_CGS);
        if ( name == "dt_snap" ) dt_snap = getvalue(value, T_MKS, T_CGS);
        if ( name == "dt_snap_tmp" ) dt_snap_tmp = getvalue(value, T_MKS, T_CGS);
        if ( name == "r_max" ) r_max = getvalue(value, L_MKS, L_CGS);
        if ( name == "r_min" ) r_min = getvalue(value, L_MKS, L_CGS);
        if ( name == "seed" ) seed = std::atoi(value.c_str());
        if ( name == "reset_step" ) reset_step = std::atoi(value.c_str());
        
        FP_t::readParameter(name, value);
        SolidDisk::readParameter(name, value);
        GasDisk::readParameter(name, value);
        Collision::readParameter(name, value);
    }
    ifs.close();

    return 0;
}

PS::S32 Parameter::checkParameter()
{
    if ( PS::Comm::getRank() == 0 ){
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
            errorMessage("nx & ny have NOT been set to satisfy nx * ny = number of process",
                         "(nx = "+ std::to_string(nx)
                         + ", ny = " + std::to_string(ny) + ").");
            return 1;
        } else if ( N_PROC_64 * 64 < PS::Comm::getNumberOfProc() ) {
            errorMessage("N_PROC_64 have NOT been set to satisfy N_PROC_64 * 64 < number of process. Please chack neighbor.h .",
                         "(nx = "+ std::to_string(nx)
                         + ", ny = " + std::to_string(ny) + ").");
            return 1;
        }
    }
        
    return 0;
}

void Parameter::showParameter(PS::F64 time_sys,
                              std::ostream & fout = std::cout)
{    
    fout << "Initial File:\t" << init_file << std::endl
         << "Output Directory:\t" << output_dir << std::endl
         << "seed          = " << seed << std::endl;
    if ( makeInit ) SolidDisk::showParameter(fout);
#ifdef GAS_DRAG
    GasDisk::showParameter(fout);
#endif
    fout << std::fixed << std::setprecision(5)
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
         << std::fixed << std::setprecision(5)
         << "r_max         = " << r_max << std::endl
         << "r_min         = " << r_min << std::endl;
    FP_t::showParameter(fout);         
    Collision::showParameter(fout);
}


PS::S32 readParameter(Parameter & param)
{
    PS::S32 i;
    if ( PS::Comm::getRank() == 0 ){
        i = param.readParameter();
    }
    PS::Comm::broadcast(&param, 1);
    PS::Comm::broadcast(&Parameter::param_file[0], 256);
    PS::Comm::broadcast(&Parameter::init_file[0], 256);
    PS::Comm::broadcast(&Parameter::output_dir[0], 256);

    FP_t::broadcastParameter();
    SolidDisk::broadcastParameter();
    GasDisk::broadcastParameter();
    Collision::broadcastParameter();
    
    PS::Comm::broadcast(&i, 1);
    return i;
}

PS::S32 checkParameter(Parameter & param)
{
    return param.checkParameter();
}

void showParameter(Parameter & param,
                   PS::F64 time_sys)
{
    const PS::F64 L = 14959787070000;
    const PS::F64 M = 1.9884e33;
    const PS::F64 T = 365.25*24.*60.*60./(2.*MY_PI);
    
    //if ( PS::Comm::getRank() == 0 ){
    if ( PS::Comm::getRank() == PS::Comm::getNumberOfProc()-1 ){
        std::cout << "Number Of Processes:\t" << PS::Comm::getNumberOfProc() << std::endl;
        std::cout << "Number Of Threads Per Process:\t" << PS::Comm::getNumberOfThread() << std::endl;
        std::cout << std::endl;

        std::cout << "#################### Parameters ####################" << std::endl;
        param.showParameter(time_sys);
        Collision::showParameter();    
        std::cout << std::resetiosflags(std::ios_base::floatfield)<<std::setprecision(8)
                  << "####################################################" << std::endl;

        char sout_param[256];
        std::ofstream fout_param;
        sprintf(sout_param,"./%s/param.dat", param.output_dir);
        if ( time_sys == 0. ) {
            fout_param.open(sout_param, std::ios::out);
        } else {
            fout_param.open(sout_param, std::ios::app);
        }
        fout_param << "####################################################" << std::endl;

#ifdef USE_QUAD
        fout_param << "Use Quadrupole Moment For Tree" << std::endl;
#else
        fout_param << "Use Monopole Moment For Tree" << std::endl;
#endif

#ifdef __AVX512DQ__
        fout_param << "Use AVX512DQ" << std::endl;
#elif defined(__AVX2__)
        fout_param << "Use AVX2" << std::endl;
#endif
        
#ifdef USE_INDIVIDUAL_CUTOFF
        fout_param << "Use Individual Cut-off" << std::endl;
#else
        fout_param << "Use Shared Cut-off Radii & Search Radii" << std::endl;
#endif
#ifdef GAS_DRAG
        fout_param << "Use Gas Drag" << std::endl;
#endif
#ifdef ISOTROPIC
        fout_param << "Use Isotropic Random Velocity" << std::endl;
#endif 
        fout_param << std::endl;

        fout_param << "Number Of Processes:\t" << PS::Comm::getNumberOfProc() << std::endl;
        fout_param << "Number Of Threads Per Process:\t" << PS::Comm::getNumberOfThread() << std::endl;
        fout_param << std::endl;

        param.showParameter(time_sys, fout_param);
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
