//***************************************************************************************
//  This is record file I/O routine.
//    This code is used by "atsevb_main.cpp"
//***************************************************************************************
#pragma once

#include<climits>
#include<cstdio>
#include<string>
#include<fstream>

namespace FILE_IO {
    
    //--- prototype
    void makeOutputDirectory(char * dir_name);
    
    //--- I/O mode
    static IO_MODE io_mode = IO_MODE::pos;
    
    //--- output cycle setting
    static PS::S32 cycle_eng_log = std::numeric_limits<PS::S32>::max();
    static PS::S32 cycle_pos     = std::numeric_limits<PS::S32>::max();
    static PS::S32 cycle_resume  = std::numeric_limits<PS::S32>::max();
    static PS::S32 cycle_vmd     = std::numeric_limits<PS::S32>::max();
    
    //--- output directry setting
    char dir_name_pos[1024];
    char dir_name_resume[1024];
    char dir_name_vmd[1024];
    
    //--- log file setting
    std::ofstream fout_energy;
    char out_eng_log[1024];
    
    //--- debug log file
    std::ofstream fout_trace;
    
    //--- header class for file I/O
    class FileHeader {
      public:
        PS::S64 n_atom;
        PS::S64 i_step;
        char    buf[512];
        
        PS::S64 readAscii(FILE* fp){
            std::fscanf(fp, "%lld\n", &i_step);
            std::fscanf(fp, "%lld\n", &n_atom);
            std::fgets(buf, sizeof(buf), fp);
            return n_atom;
        }
        void writeAscii(FILE* fp) const {
            std::fprintf(fp, "%lld\n", i_step);
            std::fprintf(fp, "%lld\n", n_atom);
            std::fprintf(fp, buf);
        }
    };
    
    //--- directry manager
    void makeOutputDirectory(char * dir_name){
        struct stat st;
        if(stat(dir_name, &st) != 0) {
            PS::S32 ret_loc = 0;
            PS::S32 ret     = 0;
            if(PS::Comm::getRank() == 0)
                ret_loc = mkdir(dir_name, 0777);
            PS::Comm::broadcast(&ret_loc, ret);
            if(ret == 0) {
                if(PS::Comm::getRank() == 0)
                    std::fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
            } else {
                std::fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
                PS::Abort();
            }
        }
    }
    
    //--- initialize
    void Init(){
        std::sprintf(dir_name_pos,    "./posdata");
        std::sprintf(dir_name_resume, "./resume");
        std::sprintf(dir_name_vmd,    "./vmd");
        makeOutputDirectory(dir_name_pos);
        makeOutputDirectory(dir_name_resume);
        makeOutputDirectory(dir_name_vmd);
        
        std::sprintf(out_eng_log, "./energy_tot.dat");
        if(PS::Comm::getRank() == 0) {
            std::cerr << out_eng_log << std::endl;
            fout_energy.open(out_eng_log);
            fout_energy << std::scientific;
            fout_energy << "istep"
                        << "   eng_tot"
                        << "   eng_kin"
                        << "   eng_pot"
                        << "   pot_LJ"
                        << "   pot_coulomb"
                        << "   pot_intra"  << endl;   // add header
            
            //--- debug log
            fout_trace.open("./trace_p0.dat");
        }
    }
    
    //--- energy log output
    void recordEnergy(PS::S64 i_step){
        
        //--- perform in Rank=0 only
        if( PS::Comm::getRank() != 0 ) return;
        
        //--- output cycle
        if( (i_step % cycle_eng_log) != 0 ) return;
        EXT_SYS_CONTROL::EnergyBuf *energy_buf = EXT_SYS_CONTROL::EnergyBuf::getInstance();
        
        fout_energy << i_step << "  " << energy_buf->getEng()
                              << "  " << energy_buf->getEngKin()
                              << "  " << energy_buf->getEngPot()
                              << "  " << energy_buf->getPotLJ()
                              << "  " << energy_buf->getPotCoulomb()
                              << "  " << energy_buf->getPotIntra() << std::endl;
        
        std::fprintf(stderr, "Istep: %10d energy error: %+e\n",
                     i_step, energy_buf->getR_error() );
    }
    
    //--- pos data output
    template<class Tpsys>
    void recordPos(Tpsys & system,
                   const PS::S64 i_step){
        
        //--- output cycle
        if( (i_step % cycle_pos) != 0 ) return;
        
        char filename[256];
        std::sprintf(filename, "%s/pos%010d.dat", dir_name_pos, i_step);
        FileHeader header;
        header.i_step = i_step;
        header.n_atom = system.getNumberOfParticleGlobal();
        std::sprintf(header.buf, "AtomId\tMolID\tAtomType\tMolType\tpos_x\tpos_y\tpos_z\n");
        
        //--- write data
        io_mode = IO_MODE::pos;
        system.writeParticleAscii(filename, header);
    }
    
    //--- resume data output
    template<class Tpsys>
    void recordResume(Tpsys & system,
                      const PS::S64 i_step){
        
        //--- output cycle
        if( (i_step % cycle_resume) != 0 ) return;
        
        //--- common file name
        char filename[256];
      //  std::sprintf(filename, "%s/resume%010d", dir_name_resume, i_step);
        std::sprintf(filename, "%s/resume%010d.dat", dir_name_resume, i_step);
        FileHeader header;
        header.i_step = i_step;
        header.n_atom = system.getNumberOfParticleLocal();
        std::sprintf(header.buf, "under developping now");
        
      //  //--- distributed file format
      //  char format[256];
      //  std::sprintf(format, "$s_%03d_%03d.dat");  // filename, total MPI Rank, process MPI Rank.
        
        //--- write data
        io_mode = IO_MODE::resume;
      //  system.writeParticleAscii(filename, format);
        system.writeParticleAscii(filename);
    }
    
    //--- VMD data output
    template<class Tpsys>
    void recordVMD(Tpsys & system,
                   const PS::S64 i_step){
        
        //--- output cycle
        if( (i_step % cycle_vmd) != 0 ) return;
        
        char filename[256];
        std::sprintf(filename, "%s/vmd%010d.pdb", dir_name_vmd, i_step);
        
        //--- write data
        io_mode = IO_MODE::VMD;
        system.writeParticleAscii(filename);
    }
    
    
  //  //--- debug output --------------------------------------------------------
  //  template<class Tpsys>
  //  void recordTrajecP0(Tpsys & system,
  //                      const PS::S64 i_step){
  //      
  //      const PS::S32 n_loc = system.getNumberOfParticleLocal();
  //      PS::S32 n_glb;
  //      FP * fp;
  //      PS::AllGatherParticle(fp, n_glb, &system[0], n_loc);  // DANGER. *fp must be delete manually.
  //      
  //      if(PS::Comm::getRank() != 0) return;
  //      
  //      for(PS::S32 i=0; i<n_glb; i++){
  //          if(system[i].getAtomID() == 0){
  //              PS::F64vec pos_tmp = system[i].getPos();
  //                         pos_tmp = Normalize::realPos(pos_tmp);
  //              fout_trace << i_step
  //                         << "  " << system[i].getAtomID()
  //                         << "  " << pos_tmp.x
  //                         << "  " << pos_tmp.y
  //                         << "  " << pos_tmp.z << endl;
  //          }
  //      }
  //      delete [] fp;  // DO NOT forget delete. it will cause memory leak!
  //  }
  //  //-------------------------------------------------------------------------
    
}
    
//--- file I/O interface in FP class
//------ write function at 1 particle in "system"
void FP::writeAscii(FILE* fp) const {
    
    //--- value buffer
    PS::F64vec pos_buff = Normalize::realPos( this->getPos() );
    PS::F64vec vel_buff = Normalize::realPos( this->getVel() );
    
    //--- position & ID data output
    if(FILE_IO::io_mode == IO_MODE::pos){
        std::fprintf(fp, "%6d\t%6d\t%6d\t%6d\t%6e\t%6e\t%6e\n",
                     this->getAtomID(), this->getMolID(),
                     static_cast<int>(this->getAtomType()),
                     static_cast<int>(this->getMolType()),
                     pos_buff.x, pos_buff.y, pos_buff.z);
            
        return; // return or throw error code.
    }
    
    //--- resume data output
    if(FILE_IO::io_mode == IO_MODE::resume){
        
        return; // return or throw error code.
    }
    
    //--- VMD output
    if(FILE_IO::io_mode == IO_MODE::VMD){
        int  type = static_cast<int>(this->getAtomType());
        std::string name;
        std::string res;
        name = ENUM_TAG::name_tag.at(type);
        res  = ENUM_TAG::residue_tag.at(type);
        
        char buf[128];
        char add_buff[32];
        
        std::sprintf(buf, "ATOM  ");                        //  1~6  "ATOM  "
        std::sprintf(add_buff, "%5d", this->getAtomID() );  //  7~11 atom id
        std::strcat(buf, add_buff );
        std::strcat(buf, " ");                              // 12    space
        std::strcat(buf, name.c_str() );                    // 13~16 atom name
        std::strcat(buf, " ");                              // 17    alternate location indicator
        std::strcat(buf, res.c_str() );                     // 18~20 residue name
        std::strcat(buf, " ");                              // 21    chain identifier
        std::sprintf(add_buff, "%4d", this->getMolID() );   // 22~25 residue sequence number
        std::strcat(buf, add_buff );
        std::strcat(buf, " ");                              // 26    code for insertion of residues
        std::strcat(buf, "   ");                            // 27~29 space
        std::sprintf(add_buff, " %7.3f", pos_buff.x );      // 30~37 X position [angstrom]
        std::strcat(buf, add_buff );
        std::sprintf(add_buff, " %7.3f", pos_buff.y );      // 38~45 X position [angstrom]
        std::strcat(buf, add_buff );
        std::sprintf(add_buff, " %7.3f", pos_buff.z );      // 46~53 X position [angstrom]
        std::strcat(buf, add_buff );
        std::strcat(buf, "\n");    // end of line
        
        std::fprintf(fp, buf);
        
        return; // return or throw error code.
    }
    
    throw std::invalid_argument("invalid I/O mode.");
}
void FP::readAscii(FILE* fp){
    
    //--- pos & ID data input  (for analysis program)
    if(FILE_IO::io_mode == IO_MODE::pos){
        
        
        return; // return or throw error code.
    }
    
    //--- resume data input
    if(FILE_IO::io_mode == IO_MODE::resume){
        
        
        return; // return or throw error code.
    }
    
    throw std::invalid_argument("invalid I/O mode.");
}



