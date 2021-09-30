//***************************************************************************************
//  This is low level functiona for intramolecular force calculation.
//***************************************************************************************
#pragma once

#include <cmath>
#include <cassert>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "ff_inter_force_func.hpp"


//--- calculate intramolecular potential --------------------------------------------

namespace FORCE {

    namespace _Impl{
        void check_max_bond_length(const PS::F64 r){
            if(r > 4.0){
                std::ostringstream oss;
                oss << "bond length is too long.\n"
                    << "  r = " << r << " must be < 4.0." << "\n";
                throw std::invalid_argument(oss.str());
            }
        }
    }

    //------ harmonic bond potential
    template <class Tcoef, class Tforce>
    void calcBondForce_harmonic_IJ(const PS::F64vec &pos_i,
                                   const PS::F64vec &pos_j,
                                   const Tcoef      &coef,
                                         Tforce     &force_i){

        PS::F64vec R_ij = Normalize::relativePosAdjustReal(pos_i - pos_j);
        PS::F64    r2   = R_ij*R_ij;
        PS::F64    r    = sqrt(r2);

        //--- assert
        #ifndef NDEBUG
        _Impl::check_max_bond_length(r);
        #endif

        //--- potential
        PS::F64 r_diff = r - coef.r0;
        force_i.addPotBond( 0.5*0.5*coef.k*r_diff*r_diff );

        //--- force
        PS::F64vec F_ij = ( -coef.k*r_diff/r )*R_ij;
        force_i.addForceIntra( F_ij );

        //--- virial
        force_i.addVirialIntra( calcVirialEPI(F_ij, R_ij) );
    }

    //------ anharmonic bond potential
    template <class Tcoef, class Tforce>
    void calcBondForce_anharmonic_IJ(const PS::F64vec  &pos_i,
                                     const PS::F64vec  &pos_j,
                                     const Tcoef       &coef,
                                           Tforce      &force_i){

        PS::F64vec R_ij = Normalize::relativePosAdjustReal(pos_i - pos_j);
        PS::F64    r2   = R_ij*R_ij;
        PS::F64    r    = sqrt(r2);

        //--- assert
        #ifndef NDEBUG
        _Impl::check_max_bond_length(r);
        #endif

        //--- potential
        constexpr PS::F64 factor = 7.0/12.0;
        PS::F64 ar  = coef.a*(r - coef.r0);
        PS::F64 ar2 = ar*ar;
        force_i.addPotBond(    0.5*coef.k*(  1.0 - ar + factor*ar2)*ar2 );

        //--- force
        PS::F64    ebp   = -coef.a*coef.k*( (2.0 - 3.0*ar) + 4.0*factor*ar2 )*ar/r;
        PS::F64vec F_ij = ebp*R_ij;
        force_i.addForceIntra( F_ij );

        //--- virial
        force_i.addVirialIntra( calcVirialEPI(F_ij, R_ij) );
    }

    //------ harminic angle potential  (i-j-k form: "j" must be center)
    template <class Tid, class Tcoef, class Tforce>
    void calcAngleForce_harmonic_IJK(const PS::F64vec &pos_i,
                                     const PS::F64vec &pos_j,
                                     const PS::F64vec &pos_k,
                                     const Tid        &id_i,
                                     const Tid        &id_j,
                                     const Tid        &id_k,
                                     const Tid        &id_tgt,
                                     const Tcoef      &coef,
                                           Tforce     &force_tgt){

        PS::F64vec R_ij = Normalize::relativePosAdjustReal(pos_i - pos_j);
        PS::F64vec R_kj = Normalize::relativePosAdjustReal(pos_k - pos_j);

        PS::F64 r_a = VEC_EXT::norm(R_ij);
        PS::F64 r_b = VEC_EXT::norm(R_kj);

        PS::F64 in_prod  = R_ij*R_kj;
        PS::F64 r_ab_inv = 1.0/(r_a*r_b);
        PS::F64 cos_tmp  = in_prod*r_ab_inv;
        PS::F64 diff     = cos_tmp - std::cos(coef.theta0);

        //--- potential
        force_tgt.addPotAngle( (1.0/3.0)*0.5*coef.k*diff*diff );

        //--- force
        PS::F64 acoef     = -coef.k*diff;
        PS::F64 r_ab2_inv = r_ab_inv*r_ab_inv;

        PS::F64vec F_ij = (acoef*r_ab2_inv)*( (r_b*r_b*cos_tmp)*R_ij - (r_a*r_b)*R_kj );
        PS::F64vec F_kj = (acoef*r_ab2_inv)*( (r_a*r_a*cos_tmp)*R_kj - (r_a*r_b)*R_ij );

        //--- select output
        if(       id_tgt == id_i){
            force_tgt.addForceIntra(-F_ij);
            force_tgt.addVirialIntra( -calcVirialEPI(F_ij, R_ij) );
        } else if(id_tgt == id_j) {
            force_tgt.addForceIntra(F_ij + F_kj);
            force_tgt.addVirialIntra(  calcVirialEPI(F_ij, R_ij)
                                     + calcVirialEPI(F_kj, R_kj) );
        } else if(id_tgt == id_k) {
            force_tgt.addForceIntra(-F_kj);
            force_tgt.addVirialIntra( -calcVirialEPI(F_kj, R_kj) );
        } else {
            std::ostringstream oss;
            oss << "id_tgt is not match to id_i, id_j, or id_k in angle pair." << "\n"
                << "   id_tgt = " << id_tgt << "\n"
                << "        i = " << id_i
                <<       "  j = " << id_j
                <<       "  k = " << id_k << "\n";
            throw std::invalid_argument(oss.str());
        }
    }

    //------ support functions for torsion force
    struct TorsionLocalVec{
        PS::F64vec R_ij;
        PS::F64vec R_ik;
        PS::F64vec R_jl;
        PS::F64vec R_kj;
        PS::F64vec R_kl;

        PS::F64vec R_a, R_b;
        PS::F64    r_a, r_b, r_a_inv, r_b_inv;

        PS::F64    phi;
        PS::F64    sin_p, cos_p;
    };
    inline TorsionLocalVec calcTorsionForce_LocalVec(const PS::F64vec &pos_i,
                                                     const PS::F64vec &pos_j,
                                                     const PS::F64vec &pos_k,
                                                     const PS::F64vec &pos_l ){

        TorsionLocalVec result;

        result.R_ij = Normalize::relativePosAdjustReal(pos_i - pos_j);
        result.R_ik = Normalize::relativePosAdjustReal(pos_i - pos_k);
        result.R_jl = Normalize::relativePosAdjustReal(pos_j - pos_l);
        result.R_kj = Normalize::relativePosAdjustReal(pos_k - pos_j);
        result.R_kl = Normalize::relativePosAdjustReal(pos_k - pos_l);

        result.R_a     = VEC_EXT::cross(result.R_ij, result.R_kj);
        result.r_a     = VEC_EXT::norm(result.R_a);
        result.r_a_inv = 1.0/result.r_a;

        result.R_b     = VEC_EXT::cross(result.R_kj, result.R_kl);
        result.r_b     = VEC_EXT::norm(result.R_b);
        result.r_b_inv = 1.0/result.r_b;

        PS::F64 cos_tmp = VEC_EXT::dot(result.R_a, result.R_b)*(result.r_a_inv*result.r_b_inv);
        if(cos_tmp < -1.0) cos_tmp = -1.0;
        if(cos_tmp >  1.0) cos_tmp =  1.0;

        PS::F32 sign = VEC_EXT::dot(result.R_kj, VEC_EXT::cross(result.R_a, result.R_b));
        if(sign < 0.0){
            sign = -1.0;
        } else {
            sign = 1.0;
        }

        result.phi   = -sign*std::abs(std::acos(cos_tmp));
        result.sin_p = std::sin(result.phi);
        result.cos_p = std::cos(result.phi);

        return result;
    }

    struct TorsionForceIntensity {
        PS::F64 eng, f;

        TorsionForceIntensity        operator +  (const TorsionForceIntensity & rv){
            this->eng += rv.eng;
            this->f   += rv.f;
            return *this;
        }
        const TorsionForceIntensity& operator += (const TorsionForceIntensity & rv){
            this->eng += rv.eng;
            this->f   += rv.f;
            return *this;
        }
    };
    template <typename Tlocal_vec>
    inline TorsionForceIntensity calcTorsionForce_intensity(const Tlocal_vec & local_vec,
                                                            const PS::F64    & k,
                                                            const PS::F64    & theta0,
                                                            const PS::S32    & n_min){

        TorsionForceIntensity result;

        //--- skip for none coefficients.
        if(k == 0.0){
            result.eng = 0.0;
            result.f   = 0.0;
            return result;
        }

        PS::F64 cos_eq = std::cos(theta0);
        PS::F64 sin_tmp, sin_1, sin_2;
        switch (n_min) {
            case 1:
                sin_tmp = cos_eq;
            break;

            case 2:
                sin_tmp = 2.0*local_vec.cos_p*cos_eq;
            break;

            case 3:
                sin_tmp = (-4.0*local_vec.sin_p*local_vec.sin_p + 3.0)*cos_eq;
            break;

            case 4:
                sin_tmp = 4.0*local_vec.cos_p*(2.0*local_vec.cos_p*local_vec.cos_p - 1.0)*cos_eq;
            break;

            case 6:
                sin_1   = 4.0*local_vec.cos_p*(2.0*local_vec.cos_p*local_vec.cos_p - 1.0)*std::cos(2.0*local_vec.phi);
                sin_2   = 2.0*local_vec.cos_p*std::cos(4.0*local_vec.phi);
                sin_tmp = (sin_1 + sin_2)*cos_eq;
            break;

            default:
                sin_tmp = std::sin( PS::F64(n_min)*local_vec.phi - theta0 )/local_vec.sin_p;
        }
        PS::F64 force_tmp = -0.5*k*PS::F64(n_min)*sin_tmp;

        result.eng = (1.0/4.0)*0.5*k*(1.0 - std::cos( PS::F64(n_min)*local_vec.phi - theta0 ));
        result.f   = force_tmp;

        return result;
    }

    //------ harminic torsion potential  (i-jk-l form)
    //  shape:      i
    //              |
    //              j---k     ---- axis ----
    //              |   |
    //   (improper) l   l (dihedral)
    template <class Tid, class Tcoef, class Tforce>
    void calcTorsionForce_harmonic_IJKL(const PS::F64vec &pos_i,
                                        const PS::F64vec &pos_j,
                                        const PS::F64vec &pos_k,
                                        const PS::F64vec &pos_l,
                                        const Tid        &id_i,
                                        const Tid        &id_j,
                                        const Tid        &id_k,
                                        const Tid        &id_l,
                                        const Tid        &id_tgt,
                                        const Tcoef      &coef,
                                              Tforce     &force_tgt){

        if( std::abs(coef.theta0)            <= 1.e-5 ||
            std::abs(coef.theta0 - Unit::pi) <= 1.e-5 ){
            //--- pass checking
        } else {
            throw std::invalid_argument("theta0 = " + std::to_string(coef.theta0)
                                        + ", must be 0.0 or pi (0.0 or 180.0 in degree).");
        }

        auto local_vec = calcTorsionForce_LocalVec(pos_i, pos_j, pos_k, pos_l);
        auto Force_tmp = calcTorsionForce_intensity(local_vec,
                                                    coef.k,
                                                    coef.theta0,
                                                    coef.n_min  );

        force_tgt.addPotTorsion(Force_tmp.eng);

        PS::F64vec F_a = (  local_vec.R_b*local_vec.r_b_inv
                          - local_vec.R_a*local_vec.r_a_inv*local_vec.cos_p)*local_vec.r_a_inv;
        PS::F64vec F_b = (  local_vec.R_a*local_vec.r_a_inv
                          - local_vec.R_b*local_vec.r_b_inv*local_vec.cos_p)*local_vec.r_b_inv;

        PS::F64vec F_tmp;
        if(       id_tgt == id_i){
            F_tmp = Force_tmp.f*VEC_EXT::cross(F_a, local_vec.R_kj);
            force_tgt.addForceIntra(F_tmp);
            force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, local_vec.R_ij) );
        } else if(id_tgt == id_j){
            F_tmp = Force_tmp.f*(  VEC_EXT::cross( F_a, local_vec.R_ik)
                                 + VEC_EXT::cross(-F_b, local_vec.R_kl) );
            force_tgt.addForceIntra(F_tmp);
            //--- virial = 0.
        } else if(id_tgt == id_k){
            F_tmp = Force_tmp.f*(  VEC_EXT::cross(-F_a, local_vec.R_ij)
                                 + VEC_EXT::cross( F_b, local_vec.R_jl) );
            force_tgt.addForceIntra(F_tmp);
            force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, local_vec.R_kj) );
        } else if(id_tgt == id_l){
            F_tmp = Force_tmp.f*VEC_EXT::cross(F_b, local_vec.R_kj);
            force_tgt.addForceIntra(F_tmp);
            force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, -local_vec.R_jl) );
        } else {
            std::ostringstream oss;
            oss << "id_tgt is not match to id_i, id_j, id_k, or id_l in torsion pair." << "\n"
                << "   id_tgt = " << id_tgt << "\n"
                << "        i = " << id_i
                <<       ", j = " << id_j
                <<       ", k = " << id_k
                <<       ", l = " << id_l << "\n";
            throw std::invalid_argument(oss.str());
        }
    }


    //------ OPLS_AA 3rd order torsion potential  (i-jk-l form)
    //  shape:      i
    //              |
    //              j---k     ---- axis ----
    //              |   |
    //   (improper) l   l (dihedral)
    template <class Tid, class Tcoef, class Tforce>
    void calcTorsionForce_OPLS_3rd_IJKL(const PS::F64vec &pos_i,
                                        const PS::F64vec &pos_j,
                                        const PS::F64vec &pos_k,
                                        const PS::F64vec &pos_l,
                                        const Tid        &id_i,
                                        const Tid        &id_j,
                                        const Tid        &id_k,
                                        const Tid        &id_l,
                                        const Tid        &id_tgt,
                                        const Tcoef      &coef,
                                              Tforce     &force_tgt){

        auto local_vec = calcTorsionForce_LocalVec(pos_i, pos_j, pos_k, pos_l);

        auto Force_tmp  = calcTorsionForce_intensity(local_vec,
                                                     coef.k,  Unit::pi, 1);
             Force_tmp += calcTorsionForce_intensity(local_vec,
                                                     coef.k2, 0.0,      2);
             Force_tmp += calcTorsionForce_intensity(local_vec,
                                                     coef.k3, Unit::pi, 3);

        force_tgt.addPotTorsion(Force_tmp.eng);

        PS::F64vec F_a = (  local_vec.R_b*local_vec.r_b_inv
                          - local_vec.R_a*local_vec.r_a_inv*local_vec.cos_p)*local_vec.r_a_inv;
        PS::F64vec F_b = (  local_vec.R_a*local_vec.r_a_inv
                          - local_vec.R_b*local_vec.r_b_inv*local_vec.cos_p)*local_vec.r_b_inv;

        PS::F64vec F_tmp;
        if(       id_tgt == id_i){
            F_tmp = Force_tmp.f*VEC_EXT::cross(F_a, local_vec.R_kj);
            force_tgt.addForceIntra(F_tmp);
            force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, local_vec.R_ij) );
        } else if(id_tgt == id_j){
            F_tmp = Force_tmp.f*(  VEC_EXT::cross( F_a, local_vec.R_ik)
                                 + VEC_EXT::cross(-F_b, local_vec.R_kl) );
            force_tgt.addForceIntra(F_tmp);
            //--- virial = 0.
        } else if(id_tgt == id_k){
            F_tmp = Force_tmp.f*(  VEC_EXT::cross(-F_a, local_vec.R_ij)
                                 + VEC_EXT::cross( F_b, local_vec.R_jl) );
            force_tgt.addForceIntra(F_tmp);
            force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, local_vec.R_kj) );
        } else if(id_tgt == id_l){
            F_tmp = Force_tmp.f*VEC_EXT::cross(F_b, local_vec.R_kj);
            force_tgt.addForceIntra(F_tmp);
            force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, -local_vec.R_jl) );
        } else {
            std::ostringstream oss;
            oss << "id_tgt is not match to id_i, id_j, id_k, or id_l in torsion pair." << "\n"
                << "   id_tgt = " << id_tgt << "\n"
                << "        i = " << id_i
                <<       ", j = " << id_j
                <<       ", k = " << id_k
                <<       ", l = " << id_l << "\n";
            throw std::invalid_argument(oss.str());
        }
    }

}
