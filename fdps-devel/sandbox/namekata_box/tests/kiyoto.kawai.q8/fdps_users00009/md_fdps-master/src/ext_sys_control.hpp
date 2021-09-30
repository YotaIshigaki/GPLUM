//***************************************************************************************
//  This is external system control functions.
//***************************************************************************************
#pragma once

#include <cmath>
#include <vector>
#include <sstream>
#include <limits>
#include <iomanip>
#include <stdexcept>
#include <cassert>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "unit.hpp"
#include "atom_move.hpp"


enum class EXT_SYS_MODE : int {
    NVE,
    NVT,
    NPT,
};

extern int debug_flag;

namespace ENUM {

    static const std::map<EXT_SYS_MODE, std::string> table_EXT_SYS_MODE_str{
        {EXT_SYS_MODE::NVE, "NVE"},
        {EXT_SYS_MODE::NVT, "NVT"},
        {EXT_SYS_MODE::NPT, "NPT"},
    };

    static const std::map<std::string, EXT_SYS_MODE> table_str_EXT_SYS_MODE{
        {"NVE", EXT_SYS_MODE::NVE},
        {"NVT", EXT_SYS_MODE::NVT},
        {"NPT", EXT_SYS_MODE::NPT},
    };

    std::string what(const EXT_SYS_MODE &e){
        if(table_EXT_SYS_MODE_str.find(e) != table_EXT_SYS_MODE_str.end()){
            return table_EXT_SYS_MODE_str.at(e);
        } else {
            using type_base = typename std::underlying_type<EXT_SYS_MODE>::type;
            std::cerr << "  EXT_SYS_MODE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in EXT_SYS_MODE.");
        }
    }

    EXT_SYS_MODE which_EXT_SYS_MODE(const std::string &str){
        if(table_str_EXT_SYS_MODE.find(str) != table_str_EXT_SYS_MODE.end()){
            return table_str_EXT_SYS_MODE.at(str);
        } else {
            std::cerr << "  EXT_SYS_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in EXT_SYS_MODE.");
        }
    }
}

inline std::ostream& operator << (std::ostream& s, const EXT_SYS_MODE &e){
    s << ENUM::what(e);
    return s;
}


namespace EXT_SYS {

    //--- sequence manager
    class Setting {
    public:
        EXT_SYS_MODE mode        = EXT_SYS_MODE::NVE;
        PS::S64      period      = 0;
        PS::F64      temperature = 300.0;
        PS::F64      pressure    = 0.0;

        Setting() = default;
        Setting(const Setting &s) = default;
        Setting(const EXT_SYS_MODE &mode,
                const PS::S64      &period,
                const PS::F64      &temperature,
                const PS::F64      &pressure){
            assert(period      >= 0);
            assert(temperature >= 0.0);
            assert(pressure    >= 0.0);
            this->mode        = mode;
            this->period      = period;
            this->temperature = temperature;
            this->pressure    = pressure;
        }
    };

    class Sequence {
    private:
        Setting              set_default;
        std::vector<Setting> set_list;

    public:
        void setDefault(const Setting &set){ this->set_default = set; }
        void addSetting(const Setting &set){ this->set_list.push_back(set); }

        void clear(){
            this->set_default = Setting(EXT_SYS_MODE::NVE, 0, 300.0, 0.0);
            this->set_list.clear();
        }

        Setting getSetting(const PS::S64 &istep) const {
            assert(istep >= 0);

            PS::S64 step_count = 0;
            for(auto set : this->set_list){
                step_count += set.period;
                if(step_count >= istep) return set;
            }
            return set_default;
        }

        void broadcast(const PS::S32 root = 0){
            COMM_TOOL::broadcast(this->set_default, root);
            COMM_TOOL::broadcast(this->set_list,    root);
        }

        void print() const {
            std::ostringstream oss;

            oss << "external system control sequence:\n";
            oss << std::setw(26) << " "
                << std::setw(6)  << "mode " << " "
                << std::setw(17) << "temperature[K] " << " "
                << std::setw(17) << "pressure[MPa]" << "\n";
            oss << std::setw(26) << "  default      : "
                << std::setw(6)  << this->set_default.mode << " "
                << std::setw(17) << std::fixed << this->set_default.temperature << " "
                << std::setw(17) << std::fixed << this->set_default.pressure
                << "\n\n";

            PS::S64 step_count = 0;
            for(auto set : this->set_list){
                oss << std::setw(10) << step_count << " ~ "
                    << std::setw(10) << step_count + set.period << " : "

                    << std::setw(6) << set.mode << " "
                    << std::setw(17) << std::fixed << set.temperature << " "
                    << std::setw(17) << std::fixed << set.pressure
                    << "\n";

                step_count += set.period;
            }
            oss << "\n";

            std::cout << oss.str() << std::flush;
        }
    };

    //--- controller
    constexpr PS::S32 max_n_chain = 5;
    constexpr PS::S32 max_n_nys   = 5;

    class Controller {
    public:
        //--- property data set
        struct State {
            PS::S32 n_chain;
            PS::S32 n_rep;
            PS::S32 n_nys;

            PS::F32 NVT_freq;
            PS::F32 NPT_freq;

            PS::F64 v_press;

            MD_EXT::fixed_vector<PS::F64, max_n_nys  > w_coef;
            MD_EXT::fixed_vector<PS::F64, max_n_chain> x_nhc;
            MD_EXT::fixed_vector<PS::F64, max_n_chain> v_nhc;
        };

    private:
        bool  init_flag = false;
        State state;

        template <class Tpsys, class Teng>
        void nose_hoover_chain(const PS::S64 &n_rigid_local,
                               const PS::F64 &dt,
                               const PS::F64 &norm_tgt_temp,
                                     Tpsys   &psys,
                                     Teng    &eng          );

        template <class Tpsys, class Teng>
        void nose_hoover_chain_andersen(const PS::S64 &n_rigid_local,
                                        const PS::F64 &dt,
                                        const PS::F64 &norm_tgt_temp,
                                        const PS::F64 &norm_tgt_press,
                                              Tpsys   &psys,
                                              Teng    &eng          );

        template <class Tarray>
        void nhc_calc_q_mass(const PS::S64 &n_deg_free,
                             const PS::F64 &norm_tgt_temp,
                                   Tarray  &q_mass,
                                   PS::F64 &g_nkt         ) const;

        template <class Tarray>
        void nhc_calc_w_dti(const PS::F64 &dt,
                                  Tarray  &w_dti);

        template <class Tarray>
        void nhc_update_vel_forward(const Tarray  &g_nhc,
                                    const PS::F64 &w_dt  );

        template <class Tarray>
        void nhc_update_vel_backward(      Tarray  &g_nhc,
                                     const Tarray  &q_mass,
                                     const PS::F64 &w_dt,
                                     const PS::F64 &norm_tgt_temp);

        void andersen_update_v_press(const PS::F64 &w_dt,
                                     const PS::F64 &g_press);

        template <class Tarray>
        PS::F64 nhc_calc_eng_thermostat(const PS::S64 &n_deg_free,
                                        const PS::F64 &norm_tgt_temp,
                                        const Tarray  &q_mass        ) const;

        void update_volume(const PS::F64 &dt);

    public:
        void init(const PS::S32 &n_chain,
                  const PS::S32 &n_rep,
                  const PS::S32 &n_nys,
                  const PS::F64 &NVT_freq,
                  const PS::F64 &NPT_freq);

        template <class Tpsys, class Teng>
        void apply(const PS::S64 &n_rigid_local,
                   const Setting &setting,
                   const PS::F64 &dt,
                         Tpsys   &psys,
                         Teng    &eng     );

        template <class Tpsys>
        void kick(const PS::F64 &dt,
                        Tpsys   &psys);
        template <class Tpsys>
        PS::F64 drift(const PS::F64 &dt,
                            Tpsys   &psys);

        Controller operator = (const State &resume){
            this->state = resume;
            return *this;
        }
        Controller operator = (const Controller &rhs){
            this->state = rhs.get_resume();
            return *this;
        }

        void broadcast(const PS::S32 root = 0){
            COMM_TOOL::broadcast(this->state,     root);
            COMM_TOOL::broadcast(this->init_flag, root);
        }

        State get_resume() const {
            return this->state;
        }
        void load_resume(const State &resume){
            this->state = resume;
            this->broadcast();
        }

        void print() const {
            std::ostringstream oss;

            oss << "external system controller:\n";
            if(this->init_flag){
                oss << "  n_chain  = " << this->state.n_chain  << "\n"
                    << "  n_rep    = " << this->state.n_rep    << "\n"
                    << "  n_nys    = " << this->state.n_nys    << "\n"
                    << "  NVT_freq = " << std::scientific << this->state.NVT_freq << "  (normalized)" << "\n"
                    << "  NPT_freq = " << std::scientific << this->state.NPT_freq << "  (normalized)" << "\n";
            } else {
                oss << "  --- not initialized ---\n";
            }
            oss << "\n";

            std::cout << oss.str() << std::flush;
        }
    };


    //--- implementation for public member
    void Controller::init(const PS::S32 &n_chain,
                          const PS::S32 &n_rep,
                          const PS::S32 &n_nys,
                          const PS::F64 &NVT_freq,
                          const PS::F64 &NPT_freq){

        assert(1 <= n_chain && n_chain <= max_n_chain);
        assert(0 <= n_rep);
        assert(1.e12 <= NVT_freq && NVT_freq <= 1.e13);
        assert(5.e10 <= NPT_freq && NPT_freq <= 5.e11);

        this->state.n_chain  = n_chain;
        this->state.n_rep    = n_rep;
        this->state.n_nys    = n_nys;
        this->state.NVT_freq = NVT_freq*2.0*Unit::pi*Unit::norm_time;
        this->state.NPT_freq = NPT_freq*2.0*Unit::pi*Unit::norm_time;

        switch (n_nys) {
            case 3:
                this->state.w_coef.clear();
                this->state.w_coef.push_back( 1.0/(2.0 - std::pow(2.0, 1.0/3.0)) );
                this->state.w_coef.push_back( 1.0 - 2.0*this->state.w_coef.at(0) );
                this->state.w_coef.push_back( this->state.w_coef.at(0) );
            break;

            case 5:
                this->state.w_coef.clear();
                this->state.w_coef.push_back( 1.0/(4.0 - std::pow(4.0, 1.0/3.0)) );
                this->state.w_coef.push_back( this->state.w_coef.at(0) );
                this->state.w_coef.push_back( 1.0 - 4.0*this->state.w_coef.at(0) );
                this->state.w_coef.push_back( this->state.w_coef.at(0) );
                this->state.w_coef.push_back( this->state.w_coef.at(0) );
            break;

            default:
                std::cerr << "  n_nys = " << n_nys << std::endl;
                throw std::invalid_argument("undefined value of 'n_nys'");
        }

        this->state.v_press = 0.0;

        this->state.v_nhc.resize(n_chain);
        this->state.x_nhc.resize(n_chain);
        std::fill(this->state.v_nhc.begin(), this->state.v_nhc.end(), 0.0);
        std::fill(this->state.x_nhc.begin(), this->state.x_nhc.end(), 0.0);

        this->init_flag = true;
    }

    template <class Tpsys, class Teng>
    void Controller::apply(const PS::S64 &n_rigid_local,
                           const Setting &setting,
                           const PS::F64 &dt,
                                 Tpsys   &psys,
                                 Teng    &eng     ){

        assert(this->init_flag);
        assert(n_rigid_local >= 0);
        assert(dt > 0.0);

        PS::F64 norm_tgt_temp  = setting.temperature/Unit::norm_temp;
        PS::F64 norm_tgt_press = (setting.pressure*1.e6)/Unit::norm_press;  // [MPa] -> [Pa]

        //--- control action
        switch (setting.mode) {
            case EXT_SYS_MODE::NVE:
                this->state.v_press = 0.0;
                eng.ext_sys   = 0.0;
            break;

            case EXT_SYS_MODE::NVT:
                this->nose_hoover_chain(n_rigid_local, dt,
                                        norm_tgt_temp,
                                        psys, eng);
                this->state.v_press = 0.0;
            break;

            case EXT_SYS_MODE::NPT:
                this->nose_hoover_chain_andersen(n_rigid_local, dt,
                                                 norm_tgt_temp, norm_tgt_press,
                                                 psys, eng);
                this->update_volume(dt);
            break;

            default:
                std::cerr << "  mode = " << ENUM::what(setting.mode) << std::endl;
                throw std::invalid_argument("undefined mode of EXT_SYS::Controller");
        }
    }

    template <class Tpsys>
    void Controller::kick(const PS::F64 &dt,
                                Tpsys   &psys){
        ATOM_MOVE::kick(dt, psys);
    }
    template <class Tpsys>
    PS::F64 Controller::drift(const PS::F64 &dt,
                                    Tpsys   &psys){
        constexpr PS::F64 e2 = 1.0/6.0;
        constexpr PS::F64 e4 = e2/20.0;
        constexpr PS::F64 e6 = e4/42.0;
        constexpr PS::F64 e8 = e6/72.0;

        PS::F64 p    = this->state.v_press*dt;
        PS::F64 a    = std::exp(p);
        PS::F64 a2   = a*a;
        PS::F64 p2   = p*p;
        PS::F64 poly = 1.0 + p2*(e2 + p2*(e4 + p2*(e6 + p2*e8)));
        PS::F64 b    = a*poly*dt;

        PS::F64 max_move = 0.0;
        const PS::S64 n_local  = psys.getNumberOfParticleLocal();
        for(PS::F64 i=0; i<n_local; ++i){
            if (debug_flag == 1 && i == 675) {
                std::cout << std::scientific;
                std::cout << "pos (bef. drift) = " 
                          << std::setprecision(std::numeric_limits<float >::max_digits10)
                          << psys[i].getPos() << std::endl;
            }
            PS::F64vec move    = psys[i].getVel()*b;
            PS::F64vec pos_new = psys[i].getPos() + Normalize::normDrift(move*a2);
            if (debug_flag == 1 && i == 675) {
                std::cout << "pos_new (aft. drift) = " 
                          << std::setprecision(std::numeric_limits<double >::max_digits10)
                          << pos_new << std::endl;
                debug_flag = 2;
            }
            if (debug_flag == 2) {
                       //pos_new = Normalize::periodicAdjustNorm(pos_new);
            }
            psys[i].addTrj(move);
            psys[i].setPos(pos_new);

            if (debug_flag == 2 && i == 675) {
                std::cout << "pos_new (aft. periodicAdjustNorm) = " 
                          << std::setprecision(std::numeric_limits<double >::max_digits10)
                          << pos_new << std::endl;
                std::cout << "pos (aft. periodicAdjustNorm) = " 
                          << std::setprecision(std::numeric_limits<float >::max_digits10)
                          << psys[i].getPos() << std::endl;
                debug_flag = 1;
            }
#if 0
            //--- check pos_new
            if ((psys[i].getPos().x < 0.0) || (1.0 <= psys[i].getPos().x) ||
                (psys[i].getPos().y < 0.0) || (1.0 <= psys[i].getPos().y) ||
                (psys[i].getPos().z < 0.0) || (1.0 <= psys[i].getPos().z)) {
                std::cout<<"**** pos_error @drift() ****"<<std::endl;
		std::cout<<"rank = "<<PS::Comm::getRank()<<std::endl;
		std::cout<<"id   = "<<psys[i].getId()<<std::endl;
		std::cout<<"pos  = "<<psys[i].getPos()<<std::endl;
            }
#endif

            //--- check largest move at step
            move     = Normalize::normDrift(move);
            max_move = std::max(max_move, move*move);
        }
        max_move = std::sqrt(max_move);

        if(max_move > 0.5){
            std::cerr << "  max_move = " << max_move << std::endl;
            throw std::logic_error("atoms speed runaway");
        }
        return max_move;
    }

    //--- implementation for private member
    template <class Tpsys, class Teng>
    void Controller::nose_hoover_chain(const PS::S64 &n_rigid_local,
                                       const PS::F64 &dt,
                                       const PS::F64 &norm_tgt_temp,
                                             Tpsys   &psys,
                                             Teng    &eng          ){

        PS::S64 n_local    = psys.getNumberOfParticleLocal();
        PS::S64 n_deg_free = 3*n_local - n_rigid_local;
                n_deg_free = PS::Comm::getSum(n_deg_free);

        PS::F64 eng_kinetic = 0.0;
        for(PS::S64 i=0; i<n_local; ++i){
            PS::F64vec v = psys[i].getVel();
            eng_kinetic += psys[i].getMass()*(v*v);   // virial = 2.0*E_kinetic
        }
        eng_kinetic = PS::Comm::getSum(eng_kinetic);

        //--- set local constant
        MD_EXT::fixed_vector<PS::F64, max_n_chain> g_nhc;
        MD_EXT::fixed_vector<PS::F64, max_n_chain> q_mass;
        MD_EXT::fixed_vector<PS::F64, max_n_nys>   w_dti;

        g_nhc.resize( this->state.n_chain);
        q_mass.resize(this->state.n_chain);
        w_dti.resize( this->state.n_nys);
        std::fill(g_nhc.begin(),  g_nhc.end(),  0.0);
        std::fill(q_mass.begin(), q_mass.end(), 0.0);
        std::fill(w_dti.begin(),  w_dti.end(),  0.0);

        PS::F64 g_nkt = 0.0;
        this->nhc_calc_q_mass(n_deg_free, norm_tgt_temp,
                              q_mass, g_nkt);

        this->nhc_calc_w_dti(dt, w_dti);

        PS::F64 scale = 1.0;
        g_nhc.at(0) = (eng_kinetic - g_nkt)/q_mass.at(0);

        //--- interaction with nose hoover chain
        for(PS::S32 i_rep=0; i_rep<this->state.n_rep; ++i_rep){
            for(PS::S32 i_nys=0; i_nys<this->state.n_nys; ++i_nys){
                //--- update thermostat velocities (heat source -> atom)
                this->nhc_update_vel_forward(g_nhc, w_dti.at(i_nys));

                //--- update particle velocity
                PS::F64 factor = std::exp(-w_dti.at(i_nys)*this->state.v_nhc.at(0));
                scale          = scale*factor;
                eng_kinetic    = eng_kinetic*(factor*factor);

                //--- update thermostat position
                for(PS::S32 i_chain=0; i_chain<this->state.n_chain; ++i_chain){
                    this->state.x_nhc.at(i_chain) += this->state.v_nhc.at(i_chain)*w_dti.at(i_nys);
                }

                //--- update force
                g_nhc.at(0) = (eng_kinetic - g_nkt)/q_mass.at(0);

                //--- update thermostat velocities (atom -> heat source)
                this->nhc_update_vel_backward(g_nhc, q_mass, w_dti.at(i_nys), norm_tgt_temp);
            }
        }

        //--- update velocity of atoms
        COMM_TOOL::broadcast(scale, 0);
        COMM_TOOL::broadcast(eng_kinetic, 0);
        for(PS::S64 i=0; i<n_local; ++i){
            psys[i].setVel( psys[i].getVel()*scale );
        }

        //--- energy of thermostat
        PS::F64 eng_nhc = this->nhc_calc_eng_thermostat(n_deg_free, norm_tgt_temp, q_mass);

        COMM_TOOL::broadcast(eng_nhc, 0);

        eng.kin     = 0.5*eng_kinetic;  // virial -> E_kinetic
        eng.ext_sys = eng_nhc;
    }

    template <class Tpsys, class Teng>
    void Controller::nose_hoover_chain_andersen(const PS::S64 &n_rigid_local,
                                                const PS::F64 &dt,
                                                const PS::F64 &norm_tgt_temp,
                                                const PS::F64 &norm_tgt_press,
                                                      Tpsys   &psys,
                                                      Teng    &eng          ){

        PS::S64 n_local    = psys.getNumberOfParticleLocal();
        PS::S64 n_deg_free = 3*n_local - n_rigid_local;
                n_deg_free = PS::Comm::getSum(n_deg_free);

        PS::F64 eng_kinetic = 0.0;
        for(PS::S64 i=0; i<n_local; ++i){
            PS::F64vec v = psys[i].getVel();
            eng_kinetic += psys[i].getMass()*(v*v);   // virial = 2.0*E_kinetic
        }
        eng_kinetic = PS::Comm::getSum(eng_kinetic);

        //--- set local constant
        MD_EXT::fixed_vector<PS::F64, max_n_chain> g_nhc;
        MD_EXT::fixed_vector<PS::F64, max_n_chain> q_mass;
        MD_EXT::fixed_vector<PS::F64, max_n_nys>   w_dti;

        g_nhc.resize( this->state.n_chain);
        q_mass.resize(this->state.n_chain);
        w_dti.resize( this->state.n_nys);
        std::fill(g_nhc.begin(),  g_nhc.end(),  0.0);
        std::fill(q_mass.begin(), q_mass.end(), 0.0);
        std::fill(w_dti.begin(),  w_dti.end(),  0.0);

        PS::F64 dummy = 0.0;
        this->nhc_calc_q_mass(n_deg_free, norm_tgt_temp,
                              q_mass, dummy);

        this->nhc_calc_w_dti(dt, w_dti);

        PS::F64 v_mass         = PS::F64(n_deg_free + 3)*norm_tgt_temp/std::pow(this->state.NPT_freq, 2);
        PS::F64 n_deg_free_inv = 1.0 + 3.0/PS::F64(n_deg_free);
        PS::F64 press_internal = (eng.virial.x + eng.virial.y + eng.virial.z);
        PS::F64 g_n1kt         = PS::F64(n_deg_free + 1)*norm_tgt_temp;

        PS::F64 scale   = 1.0;
        g_nhc.at(0)     = (eng_kinetic + v_mass*std::pow(this->state.v_press, 2) - g_n1kt)/q_mass.at(0);
        PS::F64 g_press = (   n_deg_free_inv*eng_kinetic
                            + press_internal
                            - 3.0*norm_tgt_press*Normalize::getVol()     )/v_mass;

        //--- interaction with nose hoover chain
        for(PS::S32 i_rep=0; i_rep<this->state.n_rep; ++i_rep){
            for(PS::S32 i_nys=0; i_nys<this->state.n_nys; ++i_nys){
                //--- update thermostat velocities (heat source -> atom)
                this->nhc_update_vel_forward(g_nhc, w_dti.at(i_nys));

                //--- update barostat velocity
                this->andersen_update_v_press(w_dti.at(i_nys), g_press);

                //--- update particle velocity
                PS::F64 factor = std::exp( -w_dti.at(i_nys)*(this->state.v_nhc.at(0) + n_deg_free_inv*this->state.v_press) );
                scale          = scale*factor;
                eng_kinetic    = eng_kinetic*(factor*factor);

                g_press = (  n_deg_free_inv*eng_kinetic
                           + press_internal
                           - 3.0*norm_tgt_press*Normalize::getVol() )/v_mass;

                //--- update thermostat position
                for(PS::S32 i_chain=0; i_chain<this->state.n_chain; ++i_chain){
                    this->state.x_nhc.at(i_chain) += this->state.v_nhc.at(i_chain)*w_dti.at(i_nys);
                }

                //--- update barostat velocity
                this->andersen_update_v_press(w_dti.at(i_nys), g_press);

                //--- update force
                g_nhc.at(0) = (eng_kinetic + v_mass*std::pow(this->state.v_press, 2) - g_n1kt)/q_mass.at(0);

                //--- update thermostat velocities (atom -> heat source)
                this->nhc_update_vel_backward(g_nhc, q_mass, w_dti.at(i_nys), norm_tgt_temp);
            }
        }

        //--- update velocity of atoms
        COMM_TOOL::broadcast(scale, 0);
        COMM_TOOL::broadcast(eng_kinetic, 0);
        for(PS::S64 i=0; i<n_local; ++i){
            psys[i].setVel( psys[i].getVel()*scale );
        }

        //--- energy of thermostat
        PS::F64 eng_nhc  = this->nhc_calc_eng_thermostat(n_deg_free, norm_tgt_temp, q_mass);
                eng_nhc += norm_tgt_temp*this->state.x_nhc.at(0);
                eng_nhc += 0.5*v_mass*std::pow(this->state.v_press, 2) + norm_tgt_press*Normalize::getVol();

        eng.kin     = 0.5*eng_kinetic;  // virial -> E_kinetic
        eng.ext_sys = eng_nhc;
    }

    template <class Tarray>
    void Controller::nhc_calc_q_mass(const PS::S64 &n_deg_free,
                                     const PS::F64 &norm_tgt_temp,
                                           Tarray  &q_mass,
                                           PS::F64 &g_nkt         ) const {

        PS::F64 NVT_freq_sqinv = 1.0/std::pow(this->state.NVT_freq, 2);
        g_nkt                  = PS::F64(n_deg_free)*norm_tgt_temp;

        q_mass.at(0) = g_nkt*NVT_freq_sqinv;
        for(PS::S32 i=1; i<this->state.n_chain; ++i){
            q_mass.at(i) = norm_tgt_temp*NVT_freq_sqinv;
        }
    }

    template <class Tarray>
    void Controller::nhc_calc_w_dti(const PS::F64 &dt,
                                          Tarray  &w_dti){

        for(PS::S32 i=0; i<this->state.n_nys; ++i){
            w_dti.at(i) = 0.5*this->state.w_coef.at(i)*dt/PS::F64(this->state.n_rep);
        }
    }

    template <class Tarray>
    void Controller::nhc_update_vel_forward(const Tarray  &g_nhc,
                                            const PS::F64 &w_dt  ){

        this->state.v_nhc.back() += g_nhc.back()*0.5*w_dt;

        for(PS::S32 next=this->state.n_chain-2; next>=0; --next){
            PS::S32 previous = next + 1;
            PS::F64 factor   = std::exp(-0.25*w_dt*this->state.v_nhc.at(previous));
            this->state.v_nhc.at(next) = ( this->state.v_nhc.at(next)*factor
                                           + g_nhc.at(next)*0.5*w_dt )*factor;
        }
    }

    template <class Tarray>
    void Controller::nhc_update_vel_backward(      Tarray  &g_nhc,
                                             const Tarray  &q_mass,
                                             const PS::F64 &w_dt,
                                             const PS::F64 &norm_tgt_temp){

        for(PS::S32 i_chain=0; i_chain<this->state.n_chain-1; ++i_chain){
            PS::F64 factor                = std::exp(-0.25*w_dt*this->state.v_nhc.at(i_chain+1));
            this->state.v_nhc.at(i_chain) = ( this->state.v_nhc.at(i_chain)*factor
                                              +           g_nhc.at(i_chain)*0.5*w_dt )*factor;
            g_nhc.at(i_chain+1) = ( q_mass.at(i_chain)*std::pow(this->state.v_nhc.at(i_chain), 2)
                                    - norm_tgt_temp )
                                  /q_mass.at(i_chain+1);
        }

        this->state.v_nhc.back() += g_nhc.back()*0.5*w_dt;
    }

    inline void Controller::andersen_update_v_press(const PS::F64 &w_dt,
                                                    const PS::F64 &g_press){
        PS::F64 factor      = std::exp(-0.25*w_dt*this->state.v_nhc.at(0));
        this->state.v_press = (this->state.v_press*factor + 0.5*w_dt*g_press)*factor;
    }

    template <class Tarray>
    PS::F64 Controller::nhc_calc_eng_thermostat(const PS::S64 &n_deg_free,
                                                const PS::F64 &norm_tgt_temp,
                                                const Tarray  &q_mass        ) const {

        PS::F64 eng = PS::F64(n_deg_free)*norm_tgt_temp*this->state.x_nhc.at(0);
        for(PS::S32 i=0; i<this->state.n_chain; ++i){
            eng += 0.5*q_mass.at(i)*std::pow(this->state.v_nhc.at(i), 2);
            eng += norm_tgt_temp*this->state.x_nhc.at(i);
        }
        return eng;
    }

    inline void Controller::update_volume(const PS::F64 &dt){
        assert(this->init_flag);

        PS::F64vec box_size = Normalize::getBoxSize();

        //--- cubic box only
        assert(box_size.x == box_size.y &&
               box_size.x == box_size.z   );

        PS::F64 x_press  = std::log(box_size.x*box_size.y*box_size.z)/3.0;
                x_press += this->state.v_press*dt;
                box_size = std::pow(std::exp(3.0*x_press), 1.0/3.0);

        Normalize::setBoxSize( box_size );
    }
}
