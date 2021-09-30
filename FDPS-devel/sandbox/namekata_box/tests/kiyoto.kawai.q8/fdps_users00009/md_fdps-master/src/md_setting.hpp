//***************************************************************************************
//  global data for system settings.
//***************************************************************************************
#pragma once

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "atom_class.hpp"
#include "md_coef_table.hpp"


namespace System {

    //--- setting data class: DO NOT contain pointer or container.
    class Profile {
    private:
        PS::F64 time     = 0.0;
        PS::F64 time_trj = 0.0;

    public:
        //--- for time step
        PS::S64 istep    = -1;
        PS::S64 nstep_st = -1;
        PS::S64 nstep_ed = -1;
        PS::F64 dt       = 0.0;

        //--- for Tree
        PS::F64 coef_ema      = -1.0;
        PS::F64 theta         = -1.0;
        PS::S32 n_leaf_limit  = -1;
        PS::S32 n_group_limit = -1;
        PS::S32 cycle_dinfo   = -1;

        //--- for cut_off radius
        PS::F32 cut_off_LJ    = -1.0;
        PS::F32 cut_off_intra = -1.0;

        //--- for installing molecule at initialize
        PS::F32 ex_radius = -1.0;
        PS::S32 try_limit = -1;

        //--- for recording data
        PS::S64 pos_interval    = std::numeric_limits<PS::S64>::max();
        PS::S64 pos_start       = std::numeric_limits<PS::S64>::max();
        PS::S64 resume_interval = std::numeric_limits<PS::S64>::max();
        PS::S64 resume_start    = std::numeric_limits<PS::S64>::max();
        PS::S64 pdb_interval    = std::numeric_limits<PS::S64>::max();
        PS::S64 pdb_start       = std::numeric_limits<PS::S64>::max();

        PS::S64 eng_interval  = std::numeric_limits<PS::S64>::max();
        PS::S64 eng_start     = std::numeric_limits<PS::S64>::max();
        PS::S64 prop_interval = std::numeric_limits<PS::S64>::max();
        PS::S64 prop_start    = std::numeric_limits<PS::S64>::max();

    public:

        PS::S64 get_istep() const { return this->istep; }
        PS::F64 get_dt()    const { return this->dt;    }

        //--- simulation time
        PS::F64 get_time_raw()                const { return this->time; }
        void    set_time_raw(const PS::F64 t)       { this->time = t;    }

        PS::F64 get_time()                const { return Unit::to_real_time( this->get_time_raw() ); }
        void    set_time(const PS::F64 t)       { this->set_time_raw( Unit::to_norm_time(t) );       }

        //--- trajectory integration time
        PS::F64 get_trj_time_raw()                const { return this->time_trj; }
        void    set_trj_time_raw(const PS::F64 t)       { this->time_trj = t;    }

        PS::F64 get_trj_time()                const { return Unit::to_real_time( this->get_trj_time_raw() ); }
        void    set_trj_time(const PS::F64 t)       { this->set_trj_time_raw( Unit::to_norm_time(t) );       }

        void    clear_trj_time() { this->time_trj = 0.0; }

        //--- cut off
        PS::F64 get_cut_off_intra() const { return this->cut_off_intra; }
        PS::F64 get_cut_off_LJ()    const { return this->cut_off_LJ;    }

        //--- for initializer
        PS::F64 get_ex_radius() const { return this->ex_radius;     }
        PS::F64 get_try_limit() const { return this->try_limit;     }

        //--- for record
        PS::S64 get_eng_start()  const { return this->eng_start;  }
        PS::S64 get_prop_start() const { return this->prop_start; }

        PS::S64 get_eng_interval()  const { return this->eng_interval;  }
        PS::S64 get_prop_interval() const { return this->prop_interval; }

        PS::S64 get_pos_start()    const { return this->pos_start;    }
        PS::S64 get_VMD_start()    const { return this->pdb_start;    }
        PS::S64 get_resume_start() const { return this->resume_start; }

        PS::S64 get_pos_interval()    const { return this->pos_interval;    }
        PS::S64 get_VMD_interval()    const { return this->pdb_interval;    }
        PS::S64 get_resume_interval() const { return this->resume_interval; }

        void step_next(){
            assert(this->istep >= 0  );
            assert(this->dt    >  0.0);
            assert(this->time  >= 0.0);
            ++(this->istep);
            this->time     += this->dt;
            this->time_trj += this->dt;
        }

        //--- timing control
        bool is_dinfo_update() const {
            assert(this->cycle_dinfo > 0);
            return ( ((this->istep - this->nstep_st) % this->cycle_dinfo) == 0 );
        }
        bool is_loop_continue() const {
            assert(this->istep    >= 0);
            assert(this->nstep_ed >= 0);
            return ( this->istep <= this->nstep_ed );
        }
    };
    //--- Global profile object.
    static Profile profile;

    //--- settings for initialize particle
    std::vector<std::pair<MolName, PS::S64>> model_list;
    std::vector<std::vector<Atom_FP>>        model_template;


    //--- sysc settings in MPI processes
    void broadcast_profile(const PS::S32 root = 0){
        COMM_TOOL::broadcast(profile, root);

        Normalize::broadcast_boxSize(root);

        COMM_TOOL::broadcast(model_list,     root);
        COMM_TOOL::broadcast(model_template, root);
    }


    //--- input setting to dinfo
    template<class Tdinfo>
    void InitDinfo(Tdinfo &dinfo){
        dinfo.initialize(profile.coef_ema);
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
        dinfo.setPosRootDomain(PS::F32vec{0.0, 0.0, 0.0}, PS::F32vec{1.0, 1.0, 1.0});  // fixed size for using PS::ParticleMesh
    }

    //--- input setting to tree
    template<class Ttree>
    void InitTree(const size_t &n_ptcl,
                        Ttree   &tree){
        tree.initialize(n_ptcl,
                        profile.theta,
                        profile.n_leaf_limit,
                        profile.n_group_limit);
    }

    //--- global interface
    void StepNext(){ profile.step_next(); }

    bool isDinfoUpdate() { return profile.is_dinfo_update();  }
    bool isLoopContinue(){ return profile.is_loop_continue(); }

    PS::S64 get_istep(){ return profile.get_istep(); }
    PS::F64 get_dt()   { return profile.get_dt();    }
    PS::F64 get_time() { return profile.get_time();  }

    PS::F64 get_cut_off_intra() { return profile.get_cut_off_intra(); }
    PS::F64 get_cut_off_LJ()    { return profile.get_cut_off_LJ();    }

    PS::F64 get_ex_radius()    { return profile.get_ex_radius();     }
    PS::S64 get_try_limit()    { return profile.get_try_limit();     }

    PS::S64 get_eng_start()    { return profile.get_eng_start();     }
    PS::S64 get_prop_start()   { return profile.get_prop_start();    }

    PS::S64 get_eng_interval() { return profile.get_eng_interval();  }
    PS::S64 get_prop_interval(){ return profile.get_prop_interval(); }

    PS::S64 get_VMD_start()    { return profile.get_VMD_start();    }
    PS::S64 get_pos_start()    { return profile.get_pos_start();    }
    PS::S64 get_resume_start() { return profile.get_resume_start(); }

    PS::S64 get_VMD_interval()    { return profile.get_VMD_interval();    }
    PS::S64 get_pos_interval()    { return profile.get_pos_interval();    }
    PS::S64 get_resume_interval() { return profile.get_resume_interval(); }


    //--- display system profiles
    void print_profile(){
        std::ostringstream oss;

        oss << "system conditions:\n";
        oss << "  Time steps:\n";
        oss << "    istep    = " << profile.get_istep() << "\n";
        oss << "    nstep_st = " << profile.nstep_st    << "\n";
        oss << "    nstep_ed = " << profile.nstep_ed    << "\n";
        oss << "    dt       = " << profile.get_dt()    << "\n";
        oss << "             = " << Unit::to_real_time( profile.get_dt() ) << " [fs]\n";
        oss << "\n";

        oss << "  Tree settings:\n";
        oss << "    n_leaf_limit  = " << profile.n_leaf_limit  << "\n";
        oss << "    coef_ema      = " << profile.coef_ema      << "\n";
        oss << "    theta         = " << profile.theta         << "\n";
        oss << "    n_group_limit = " << profile.n_group_limit << "\n";
        oss << "    cycle_dinfo   = " << profile.cycle_dinfo   << "\n";
        oss << "\n";

        oss << "  Cut_off settings:\n";
        oss << "    cut_off_intra   = " << std::setw(9) << std::setprecision(7) << profile.get_cut_off_intra() << " [angstrom]\n";
        oss << "    cut_off_LJ      = " << std::setw(9) << std::setprecision(7) << profile.get_cut_off_LJ()    << " [angstrom]\n";
        oss << "    cut_off_coulomb :  fixed as " << Normalize::normCutOff_PM() << " at normalized space.\n";
        oss << "\n";

        oss << "  Init molecule settings:\n";
        oss << "    box       = (" << Normalize::getBoxSize().x
            << ", "                << Normalize::getBoxSize().y
            << ", "                << Normalize::getBoxSize().z << ") [angstrom]\n";
        oss << "    ex_radius = "  << std::setw(9) << std::setprecision(7) << profile.get_ex_radius() << " [angstrom]\n";
        oss << "    try_limit = "  << std::setw(9) <<                         profile.get_try_limit() << " [times]\n";
        oss << "\n";

        oss << "model conditions:\n";
        if(System::model_list.size() != 0){
            for(auto m : System::model_list){
                oss << "    model_name = " << m.first  << "\n";
                oss << "    n_molecule = " << m.second << "\n";
                oss << "\n";
            }
        } else {
            oss << "    free particle.\n";
            oss << "\n";
        }

        oss << "recording conditions:\n";
        oss << "              " << std::setw(15) << "start    " << " | "
                                << std::setw(15) << "interval  \n";
        oss << "    pos     : " << std::setw(15) << profile.get_pos_start()       << " | "
                                << std::setw(15) << profile.get_pos_interval()    << "\n";
        oss << "    resume  : " << std::setw(15) << profile.get_resume_start()    << " | "
                                << std::setw(15) << profile.get_resume_interval() << "\n";
        oss << "    VMD     : " << std::setw(15) << profile.get_VMD_start()       << " | "
                                << std::setw(15) << profile.get_VMD_interval()    << "\n";
        oss << "\n";
        oss << "    energy  : " << std::setw(15) << profile.get_eng_start()     << " | "
                                << std::setw(15) << profile.get_eng_interval()  << "\n";
        oss << "    property: " << std::setw(15) << profile.get_prop_start()    << " | "
                                << std::setw(15) << profile.get_prop_interval() << "\n";
        oss << "\n";

        //--- output
        std::cout << oss.str() << std::flush;
    }

}
