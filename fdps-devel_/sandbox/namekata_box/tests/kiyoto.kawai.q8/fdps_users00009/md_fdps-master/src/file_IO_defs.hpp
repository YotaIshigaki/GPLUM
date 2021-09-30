//***************************************************************************************
//  This is record file I/O routine.
//    This code is used by "md_fdps_main.cpp"
//***************************************************************************************
#pragma once

#include <cstdio>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>


namespace FILE_IO {

    namespace TOOL {

        //--- directory manager
        void makeOutputDirectory(const std::string &dir_name){
            struct stat st;
            if(stat(dir_name.c_str(), &st) != 0) {
                PS::S32 ret_loc = 0;
                PS::S32 ret     = 0;
                if(PS::Comm::getRank() == 0){
                    ret_loc = mkdir(dir_name.c_str(), 0777);
                }
                PS::Comm::broadcast(&ret_loc, ret);
                if(ret == 0) {
                    if(PS::Comm::getRank() == 0){
                        std::cout << "Directory " << dir_name << " is successfully made." << std::endl;
                    }
                } else {
                    std::cerr << "Directory " << dir_name << " fails to be made." << std::endl;
                    PS::Abort();
                }
            }
        }

        //--- line reader
        std::vector<std::string> line_to_str_list(FILE *fp, const std::string &delim){
            std::vector<std::string> result;
            const PS::S32 len = 1024;
            char buf[len] = {0};
            if( std::fgets(buf, len, fp) != NULL ){
                std::string line{buf};

                //--- remove "\n" for std::fgets
                if(line.back() == '\n'){
                    line.pop_back();
                }

                STR_TOOL::removeCR(line);
                result = STR_TOOL::split(line, delim);
            } else {
                result.clear();
            }
            return result;
        }

        //--- Ascii header reader
        void load_Ascii_header_info(const std::string &mark,
                                    const std::string &delim,
                                    const PS::S32      n_limit,
                                          FILE        *fp,
                                          std::map<std::string, std::vector<std::string>> &info_map){
            info_map.clear();
            assert(n_limit > 0);

            PS::S32 n_try = 0;
            while( info_map.count(mark) == 0 ){
                ++n_try;
                const auto str_list = line_to_str_list(fp, delim);
                if(str_list.size() > 0){
                    info_map[str_list[0]] = str_list;
                }

                if(n_try >= n_limit){
                    std::ostringstream oss;
                    oss << "the mark of header end, '" << mark << "' was not found." << "\n"
                        << "   info map:" << "\n";
                    for(const auto& e : info_map){
                        oss << "      key = " << e.first << ", str_list =";
                        for(const auto& s : e.second){
                            oss << " " << s;
                        }
                        oss << "\n";
                    }

                    throw std::logic_error(oss.str());
                }
            }
        }

    }

    namespace _Impl {

        class FileManager_base {
        private:
            PS::S64 start    = 0;
            PS::S64 interval = 0;

            std::string data_dir  = ".";
            std::string name_pre  = "data";
            std::string name_post = ".dat";

        public:
            void init(const PS::S64 start,
                      const PS::S64 interval,
                      const std::string &dir,
                      const std::string &name_pre,
                      const std::string &name_post){
                this->start    = start;
                this->interval = interval;
                this->data_dir = dir;
                this->name_pre  = name_pre;
                this->name_post = name_post;

                TOOL::makeOutputDirectory(this->data_dir);
            }

            PS::S64 get_start()    const { return this->start;    }
            PS::S64 get_interval() const { return this->interval; }

            std::string get_dir() const { return this->data_dir; }

            std::string get_file_name(const PS::S64 i_step,
                                      const size_t  n_digit = 10) const {
                std::ostringstream oss;
                oss << this->data_dir
                    << "/"
                    << this->name_pre
                    << std::setfill('0') << std::setw(n_digit) << i_step
                    << this->name_post;
                return oss.str();
            }
        };

    }
}
