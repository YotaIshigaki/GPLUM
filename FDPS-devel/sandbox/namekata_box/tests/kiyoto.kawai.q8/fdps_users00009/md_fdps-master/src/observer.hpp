//***************************************************************************************
//  This is observer functions.
//***************************************************************************************
#pragma once

#include <vector>
#include <string>
#include <climits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "unit.hpp"


namespace Observer {

    class FilePrinter {
    private:
        PS::S32       rank = 0;
        std::string   name;
        std::ofstream ofs;

    public:
        void file_init(const std::string &name, const PS::S32 &rank){
            if( this->ofs.is_open() ){ this->ofs.close(); }
            this->rank = rank;
            this->name = name;

            if(PS::Comm::getRank() != this->rank) return;
            this->ofs.open(name, std::ios::trunc);
        }

        void print(const std::string &str){
            if(PS::Comm::getRank() != this->rank) return;
            this->ofs << str;
            this->ofs.flush();   // flush write buffer at every call
        }

        std::string file_name() const { return this->name; }
        PS::S32     getRank()   const { return this->rank; }
    };

    class Energy {
    private:
        FilePrinter printer;
        PS::S64     start   = std::numeric_limits<PS::S64>::max();

    public:
        PS::F64 bond;
        PS::F64 angle;
        PS::F64 torsion;

        PS::F64 vdw;
        PS::F64 coulomb;

        PS::F64    kin;
        PS::F64    ext_sys;
        PS::F64vec virial;

        PS::F64 density;
        PS::F64 n_atom;    // buffer for property calculation

        //--- manipulator
        void clear(){
            this->bond    = 0.0;
            this->angle   = 0.0;
            this->torsion = 0.0;

            this->vdw     = 0.0;
            this->coulomb = 0.0;

            this->kin     = 0.0;
            this->ext_sys = 0.0;
            this->virial  = 0.0;

            this->density = 0.0;
            this->n_atom  = 0.0;
        }

    private:
        void _mad(const Energy &rv, const PS::F64 sign = 1.0){
            this->bond    += sign*rv.bond;
            this->angle   += sign*rv.angle;
            this->torsion += sign*rv.torsion;

            this->vdw     += sign*rv.vdw;
            this->coulomb += sign*rv.coulomb;

            this->kin     += sign*rv.kin;
            this->ext_sys += sign*rv.ext_sys;
            this->virial  += sign*rv.virial;

            this->density += sign*rv.density;
            this->n_atom  += sign*rv.n_atom;
        }
        void _assign(const Energy &rv){
            this->clear();
            this->_mad(rv);
        }
        void _multiple(const PS::F64 &r){
            Energy buf;
            buf._mad(*this, r);
            this->_assign(buf);
        }

    public:
        //--- constructor
        Energy(){
            this->clear();
        }
        Energy(const Energy &rv){
            this->_assign(rv);
        }

        //--- operator
        Energy& operator = (const Energy &rv){
            this->_assign(rv);
            return *this;
        }
        Energy  operator +  (const Energy &rv) { this->_mad(rv,  1.0); return *this; }
        Energy  operator -  (const Energy &rv) { this->_mad(rv, -1.0); return *this; }
        Energy& operator += (const Energy &rv) { this->_mad(rv,  1.0); return *this; }
        Energy& operator -= (const Energy &rv) { this->_mad(rv, -1.0); return *this; }

        Energy  operator *  (const PS::F64 &r) { this->_multiple(r);     return *this; }
        Energy  operator /  (const PS::F64 &r) { this->_multiple(1.0/r); return *this; }
        Energy& operator *= (const PS::F64 &r) { this->_multiple(r);     return *this; }
        Energy& operator /= (const PS::F64 &r) { this->_multiple(1.0/r); return *this; }

        PS::F64 inter() const {
            return   this->vdw
                   + this->coulomb;
        }
        PS::F64 intra() const {
            return   this->bond
                   + this->angle
                   + this->torsion;
        }
        PS::F64 kinetic() const {
            return this->kin;
        }
        PS::F64 total() const {
            return   this->inter()
                   + this->intra()
                   + this->kinetic()
                   + this->ext_sys;
        }

        //--- sampling from PS::ParticleSystem<FP>
        template <class Tpsys>
        void getEnergy(const Tpsys &psys){
            Energy  buf;
                    buf.clear();
            PS::S64 n_local  = psys.getNumberOfParticleLocal();

            for(PS::S64 i=0; i<n_local; ++i){
                buf.bond    += psys[i].getPotBond();
                buf.angle   += psys[i].getPotAngle();
                buf.torsion += psys[i].getPotTorsion();

                buf.vdw     += psys[i].getPotLJ();
                buf.coulomb += psys[i].getPotCoulomb()*psys[i].getCharge();

                PS::F64vec v  = psys[i].getVel();
                buf.kin      += 0.5*psys[i].getMass()*(v*v);

            //    buf.ext_sys  = 0.0;

                buf.virial  += psys[i].getVirial();

                buf.density += psys[i].getMass();   // total mass
            //    buf.n_atom   = 0;
            }
            buf.n_atom = PS::F64(n_local);

            //--- sumation
            buf.bond    = PS::Comm::getSum(buf.bond);
            buf.angle   = PS::Comm::getSum(buf.angle);
            buf.torsion = PS::Comm::getSum(buf.torsion);
            buf.vdw     = PS::Comm::getSum(buf.vdw);
            buf.coulomb = PS::Comm::getSum(buf.coulomb);
            buf.kin     = PS::Comm::getSum(buf.kin);
        //    buf.ext_sys = PS::Comm::getSum(buf.ext_sys);
            buf.virial  = PS::Comm::getSum(buf.virial);
            buf.density = PS::Comm::getSum(buf.density);
            buf.n_atom  = PS::Comm::getSum(buf.n_atom);

            buf.ext_sys = this->ext_sys;

            buf.density = buf.density*Normalize::getVolInv();   // mass -> density
            *this = buf;
        }

        std::string header() const {
            std::ostringstream oss;
            const size_t word_length = 16;

            oss << std::left << std::setw(12) << "  Time[fs]";
            oss << std::left << std::setw(12) << "  step";

            oss << std::left << std::setw(word_length) << "  Total";
            oss << std::left << std::setw(word_length) << "  kinetic";
            oss << std::left << std::setw(word_length) << "  bond";
            oss << std::left << std::setw(word_length) << "  angle";
            oss << std::left << std::setw(word_length) << "  torsion";
            oss << std::left << std::setw(word_length) << "  vdw";
            oss << std::left << std::setw(word_length) << "  coulomb";
            oss << std::left << std::setw(word_length) << "  virial_x";
            oss << std::left << std::setw(word_length) << "  virial_y";
            oss << std::left << std::setw(word_length) << "  virial_z";
            oss << std::left << std::setw(word_length) << "  ext_sys";

            oss << "\n";

            return oss.str();
        }
        template <class Tprof>
        std::string str(const Tprof &prof) const {
            std::ostringstream oss;
            const size_t word_length = 16;
            const PS::F64 n_atom_inv = 1.0/this->n_atom;

            oss << std::setw(12) << std::setprecision(8) << prof.get_time();
            oss << std::setw(12) << prof.get_istep();

            oss << std::scientific;

            //--- printing value: energy per 1 atom.
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->total();
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->kinetic();
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->bond;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->angle;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->torsion;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->vdw;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->coulomb;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->virial.x;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->virial.y;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->virial.z;
            oss << std::setw(word_length) << std::setprecision(8) << n_atom_inv*this->ext_sys;

            oss << "\n";

            return oss.str();
        }
        void display(const PS::S64 &i_step) const {
            if(PS::Comm::getRank() != 0) return;

            std::cout << "step " << i_step;

            std::cout << " | energy total= " << this->total();
            std::cout << " kin= "            << this->kinetic();
            std::cout << " intra= "          << this->intra();
            std::cout << " inter= "          << this->inter();

            std::cout << std::endl;
        }

        void file_init(const std::string &name,
                       const PS::S64      start = 0,
                       const PS::S32      proc  = 0){
            this->printer.file_init(name, proc);
            this->printer.print( this->header() );
            this->start = start;
        }
        template <class Tprof>
        void record(const Tprof &prof) {
            if(prof.get_istep() < this->start) return;
            this->printer.print( this->str(prof) );
        }
    };

    class Property {
    private:
        FilePrinter printer;
        PS::S64     start   = std::numeric_limits<PS::S64>::max();

    public:
        PS::F64 temperature;
        PS::F64 density;
        PS::F64 pressure;

        PS::F64 press_temperature;
        PS::F64 press_virial;

        //--- manipulator
        void clear(){
            this->temperature = 0.0;
            this->density     = 0.0;
            this->pressure    = 0.0;

            this->press_temperature = 0.0;
            this->press_virial      = 0.0;
        }

    private:
        void _mad(const Property &rv, const PS::F64 sign = 1.0){
            this->temperature += sign*rv.temperature;
            this->density     += sign*rv.density;
            this->pressure    += sign*rv.pressure;

            this->press_temperature += sign*rv.press_temperature;
            this->press_virial      += sign*rv.press_virial;
        }
        void _assign(const Property &rv){
            this->clear();
            this->_mad(rv);
        }
        void _multiple(const PS::F64 &r){
            Property buf;
            buf._mad(*this, r);
            this->_assign(buf);
        }

    public:
        //--- constructor
        Property() = default;
        Property(const Property& rv){
            this->_assign(rv);
        }

        //--- operator
        Property& operator = (const Property &rv){
            this->_assign(rv);
            return *this;
        }
        Property  operator +  (const Property &rv) { this->_mad(rv,  1.0); return *this; }
        Property  operator -  (const Property &rv) { this->_mad(rv, -1.0); return *this; }
        Property& operator += (const Property &rv) { this->_mad(rv,  1.0); return *this; }
        Property& operator -= (const Property &rv) { this->_mad(rv, -1.0); return *this; }

        Property  operator *  (const PS::F64 &r)   { this->_multiple(r);     return *this; }
        Property  operator /  (const PS::F64 &r)   { this->_multiple(1.0/r); return *this; }
        Property& operator *= (const PS::F64 &r)   { this->_multiple(r);     return *this; }
        Property& operator /= (const PS::F64 &r)   { this->_multiple(1.0/r); return *this; }

        template <class Teng>
        void getProperty(const Teng &eng){
            //--- temperature [K]
            this->temperature = eng.kinetic()*(2.0/(3.0*eng.n_atom - 1.0))*Unit::norm_temp;

            //--- density [g/cm^3] = 10^-3 [kg/m^3]
            this->density     = eng.density*Unit::norm_dens*1.e-3;

            //--- pressure [Pa]
            this->press_temperature = Normalize::getVolInv()*Unit::norm_press*eng.kinetic()*(2.0/3.0);
            this->press_virial      = Normalize::getVolInv()*Unit::norm_press*(  eng.virial.x
                                                                               + eng.virial.y
                                                                               + eng.virial.z)/3.0;

            this->pressure = this->press_temperature
                           + this->press_virial;
        }

        std::string header(){
            std::ostringstream oss;

            oss << std::left << std::setw(12) << "  Time[fs]";
            oss << std::left << std::setw(12) << "  step";

            oss << std::left << std::setw(17) << "  Temperature[K]";
            oss << std::left << std::setw(17) << "  density[g/cm^3]";
            oss << std::left << std::setw(17) << "  pressure[Pa]";
            oss << std::left << std::setw(17) << "  press_temp[Pa]";
            oss << std::left << std::setw(17) << "  press_virial[Pa]";

            oss << "\n";

            return oss.str();
        }
        template <class Tprof>
        std::string str(const Tprof &prof){
            std::ostringstream oss;

            oss << std::setw(12) << std::setprecision(8) << prof.get_time();
            oss << std::setw(12) << prof.get_istep();

            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->temperature;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->density;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->pressure;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->press_temperature;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->press_virial;

            oss << "\n";

            return oss.str();
        }
        void display(const PS::S64 &i_step){
            if(PS::Comm::getRank() != 0) return;

            std::cout << "step " << i_step
                      << " | T= "    << this->temperature << "[K]"
                      <<   " Dens= " << this->density     << "[g/cm^3]"
                      <<   " P= "    << this->pressure    << "[Pa]";
            std::cout << std::endl;
        }

        void file_init(const std::string &name,
                       const PS::S64      start = 0,
                       const PS::S32      proc  = 0){
            this->printer.file_init(name, proc);
            this->printer.print( this->header() );
            this->start = start;
        }
        template <class Tprof>
        void record(const Tprof &prof) {
            if(prof.get_istep() < this->start) return;
            this->printer.print( this->str(prof) );
        }
    };


    template <class Tprop>
    class MovingAve {
    private:
        Tprop ref;
        Tprop sum;

        PS::S64 start = std::numeric_limits<PS::S64>::max();
        PS::S64 cycle = 1;
        PS::S64 count = 0;

    public:
        MovingAve<Tprop>(){
            this->ref.clear();
            this->sum.clear();
        }
        void clear(){
            this->sum.clear();
            this->count = 0;
        }
        void clearReference(){
            this->ref.clear();
        }
        void file_init(const std::string &name,
                       const PS::S64      start = 0,
                       const PS::S64      cycle = 1,
                       const PS::S32      rank  = 0 ){
            assert(start >= 0);
            assert(cycle >  0);
            assert(rank  >= 0);
            this->sum.file_init(name, 0, rank);
            this->start = start;
            this->cycle = cycle;

            this->clear();
            this->clearReference();
        }

        void setReference(const Tprop &prop){
            this->ref = prop;
        }
        template <class Tprofile>
        void record(const Tprofile &profile,
                    const Tprop    &prop   ){

            if(profile.get_istep() < this->start) return;

            //--- add data
            this->sum += prop;
            ++(this->count);

            if(this->count < this->cycle) return;

            //--- output avarage value
            this->sum /= static_cast<PS::F64>(this->count);
            this->sum -= this->ref;
            this->sum.record(profile);

            this->sum.clear();
            this->count = 0;
        }
    };
}
