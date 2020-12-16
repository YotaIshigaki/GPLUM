/**************************************************************************************************/
/**
* @file  md_ext_normalize.hpp
* @brief normalize interface for position, potential, and force.
*/
/**************************************************************************************************/
#pragma once

#include <cmath>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <iomanip>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>


/**
* @brief normalize interface for position, potential, and force.
*/
extern int debug_flag;

namespace Normalize {

    //! @brief system box size.
    PS::F64vec boxdh_{1.0, 1.0, 1.0};

    //! @brief inverse value of system box size.
    PS::F64vec boxdh_inv_{1.0, 1.0, 1.0};

    //! @brief coefficient for force from PS::PM::ParticleMesh.
    PS::F64vec coefPMForce{1.0, 1.0, 1.0};

    /**
    * @brief volume of a mesh.
    * @details "SIZE_OF_MESH" is defined in "$(PS_DIR)/src/particle_mesh/param_fdps.h".
    * @details The recommended value of "SIZE_OF_MESH" is N^(1/3)/2. (N is the total number of particles)
    */
    PS::F64vec mesh_dSize{ 1.0/SIZE_OF_MESH,
                           1.0/SIZE_OF_MESH,
                           1.0/SIZE_OF_MESH };

   /**
   * @brief set boxdh_inv_. internal use only.
   */
    void setNormalizePosition_() {
        boxdh_inv_.x = 1.0/boxdh_.x;
        boxdh_inv_.y = 1.0/boxdh_.y;
        boxdh_inv_.z = 1.0/boxdh_.z;
    }

    /**
    * @brief set coefPMForce. internal use only.
    */
    void setCoefPM_(){
        coefPMForce.x = 1.0*(boxdh_inv_.x*boxdh_inv_.x);
        coefPMForce.y = 1.0*(boxdh_inv_.y*boxdh_inv_.y);
        coefPMForce.z = 1.0*(boxdh_inv_.z*boxdh_inv_.z); // squared
    }


    //--- access functions (must be call only below functions)
    /**
    * @brief setting the real domain size.
    * @param[in] box size of domain.
    */
    void setBoxSize(const PS::F64vec & box){
        boxdh_ = box;
        setNormalizePosition_();
        setCoefPM_();
    }

    PS::F64vec getBoxSize(){
        return boxdh_;
    }
    PS::F64vec getBoxSizeInv(){
        return boxdh_inv_;
    }
    PS::F64 getVol(){
        return (boxdh_.x*boxdh_.y*boxdh_.z);
    }
    PS::F64 getVolInv(){
        return (boxdh_inv_.x*boxdh_inv_.y*boxdh_inv_.z);
    }

    /**
    * @brief broadcast setting of domain size.
    * @param[in] ref ID of source process.
    */
    void broadcast_boxSize(const PS::S32 &ref = 0){
        PS::Comm::broadcast(&boxdh_, 1, ref);
        setNormalizePosition_();
        setCoefPM_();
    }

    //! @brief convert real pos to normalized pos.
    template <class Tf>
    PS::Vector3<Tf> normPos(const PS::Vector3<Tf> &pos){
        PS::Vector3<Tf> result{ static_cast<Tf>(pos.x*boxdh_inv_.x),
                                static_cast<Tf>(pos.y*boxdh_inv_.y),
                                static_cast<Tf>(pos.z*boxdh_inv_.z) };
        return result;
    }
    //! @brief convert normalized pos to real pos.
    template <class Tf>
    PS::Vector3<Tf> realPos(const PS::Vector3<Tf> &pos){
        PS::Vector3<Tf> result{ static_cast<Tf>(pos.x*boxdh_.x),
                                static_cast<Tf>(pos.y*boxdh_.y),
                                static_cast<Tf>(pos.z*boxdh_.z) };
        return result;
    }
    //! @brief convert normalized PM force to real PM force.
    template <class Tf>
    PS::Vector3<Tf> realPMForce(const PS::Vector3<Tf> &Fpm){
        PS::Vector3<Tf> result{ static_cast<Tf>(Fpm.x*coefPMForce.x),
                                static_cast<Tf>(Fpm.y*coefPMForce.y),
                                static_cast<Tf>(Fpm.z*coefPMForce.z) };
        return result;
    }
    //! @brief convert normalized PM potential to real PM potential.
    //! @details the domain must be cubic (X=Y=Z).
    template <class Tf>
    Tf realPMPotential(const Tf &Fpm){
        //--- treat cubic system only
        if(boxdh_.x != boxdh_.y ||
           boxdh_.x != boxdh_.z ){
            std::ostringstream oss;
            oss << "the function 'Normalize::realPMPotential()' can be used cubic real space only." << "\n"
                << "    real domain = " << boxdh_ << "\n";
            throw std::logic_error(oss.str());
        }
        return boxdh_inv_.x*Fpm;
    }
    //! @brief convert real cutoff length to normalized cutoff length.
    template <class Tf>
    Tf normCutOff(const Tf &realCutOff){
        Tf len_inv = std::max( {boxdh_inv_.x,
                                boxdh_inv_.y,
                                boxdh_inv_.z},
                                std::greater<PS::F64>() );
        return realCutOff*len_inv;
    }
    //! @brief convert normalized cutoff length to real cutoff length.
    template <class Tf>
    Tf realCutOff(const Tf &normCutOff){
        Tf len = std::min( {boxdh_.x,
                            boxdh_.y,
                            boxdh_.z},
                            std::less<PS::F64>() );
        return normCutOff*len;
    }
    //! @brief get cutoff length for short-range interaction with PS::PM::ParticleMesh.
    //! @details cutoff length is 3.0/SIZE_OF_MESH.
    inline PS::F64 normCutOff_PM(){
        PS::F64 r = std::max( {mesh_dSize.x,
                               mesh_dSize.y,
                               mesh_dSize.z},
                              std::greater<PS::F64>() );
        return r*3.0;
    }

    //! @brief adjustment position into normalized periodic boundary system [0.0 ~ 1.0).
    template <class Tf>
    PS::Vector3<Tf> periodicAdjustNorm(const PS::Vector3<Tf> &pos_norm){
        PS::Vector3<Tf> pos_new = pos_norm;
        if (debug_flag == 2) {
            std::cout << std::scientific;
            std::cout << "[in periodicAdjustNorm]" << std::endl;
            std::cout << "   pos_norm       = " 
                      << std::setprecision(std::numeric_limits<float >::max_digits10)
                      << pos_norm << std::endl;
            std::cout << "   pos_new (bef.) = " 
                      << std::setprecision(std::numeric_limits<float >::max_digits10)
                      << pos_new  << std::endl;
        }
        if(pos_new.x < 0.0)  pos_new.x += 1.0;
        if(pos_new.y < 0.0)  pos_new.y += 1.0;
        if(pos_new.z < 0.0)  pos_new.z += 1.0;
        if(pos_new.x >= 1.0) pos_new.x -= 1.0;
        if(pos_new.y >= 1.0) pos_new.y -= 1.0;
        if(pos_new.z >= 1.0) pos_new.z -= 1.0;
        if (debug_flag == 2) {
            std::cout << " pos_new (aft.) = " 
                      << std::setprecision(std::numeric_limits<float >::max_digits10)
                      << pos_new << std::endl;
        }
        return pos_new;
    }
    //! @brief adjustment position into real periodic boundary system [0.0 ~ getBoxSize() ).
    template <class Tf>
    PS::Vector3<Tf> periodicAdjustReal(const PS::Vector3<Tf> &pos_real){
        PS::Vector3<Tf> pos_new = pos_real;
        if(pos_new.x <  0.0     ) pos_new.x += boxdh_.x;
        if(pos_new.y <  0.0     ) pos_new.y += boxdh_.y;
        if(pos_new.z <  0.0     ) pos_new.z += boxdh_.z;
        if(pos_new.x >= boxdh_.x) pos_new.x -= boxdh_.x;
        if(pos_new.y >= boxdh_.y) pos_new.y -= boxdh_.y;
        if(pos_new.z >= boxdh_.z) pos_new.z -= boxdh_.z;
        return pos_new;
    }

    //! @brief adjustment relative position into normalized periodic boundary system [-0.5 ~ 0.5).
    template <class Tf>
    PS::Vector3<Tf> relativePosAdjustNorm(const PS::Vector3<Tf> &pos_norm){
        PS::Vector3<Tf> pos_new = pos_norm;
        if(pos_new.x <  -0.5) pos_new.x += 1.0;
        if(pos_new.y <  -0.5) pos_new.y += 1.0;
        if(pos_new.z <  -0.5) pos_new.z += 1.0;
        if(pos_new.x >=  0.5) pos_new.x -= 1.0;
        if(pos_new.y >=  0.5) pos_new.y -= 1.0;
        if(pos_new.z >=  0.5) pos_new.z -= 1.0;
        return pos_new;
    }
    //! @brief adjustment relative position into real periodic boundary system [ -0.5*getBoxSize() ~ 0.5*getBoxSize() ).
    template <class Tf>
    PS::Vector3<Tf> relativePosAdjustReal(const PS::Vector3<Tf> &pos_real){
        PS::Vector3<Tf> pos_new = pos_real;
        if(pos_new.x <  -0.5*boxdh_.x) pos_new.x += boxdh_.x;
        if(pos_new.y <  -0.5*boxdh_.y) pos_new.y += boxdh_.y;
        if(pos_new.z <  -0.5*boxdh_.z) pos_new.z += boxdh_.z;
        if(pos_new.x >=  0.5*boxdh_.x) pos_new.x -= boxdh_.x;
        if(pos_new.y >=  0.5*boxdh_.y) pos_new.y -= boxdh_.y;
        if(pos_new.z >=  0.5*boxdh_.z) pos_new.z -= boxdh_.z;
        return pos_new;
    }

    //! @brief check the position is out of domain or not
    //! @param[in] pos_norm normalized position.
    //! @return bool True: out of domain. False: exists in domain.
    template <class Tf>
    bool checkPosInSpace(const PS::Vector3<Tf> &pos_norm){
        if(pos_norm.x <  0.0 ||
           pos_norm.y <  0.0 ||
           pos_norm.z <  0.0 ||
           pos_norm.x >= 1.0 ||
           pos_norm.y >= 1.0 ||
           pos_norm.z >= 1.0   ){
            return false;
        } else {
            return true;
        }
    }

    //! @brief convert real drift dispracement to normalized drift.
    template <class Tf>
    PS::Vector3<Tf> normDrift(const PS::Vector3<Tf> &move){
        return normPos(move);
    }

    //------ length converter
    //! @brief convert real X length to normalized X length.
    template <class Tf>
    Tf normXLen(const Tf &x_real){
        Tf x_norm = x_real * boxdh_inv_.x;
        return x_norm;
    }
    //! @brief convert normalized X length to real X length.
    template <class Tf>
    Tf realXLen(const Tf &x_norm){
        Tf x_real = x_norm * boxdh_.x;
        return x_real;
    }
    //! @brief convert real Y length to normalized Y length.
    template <class Tf>
    Tf normYLen(const Tf &y_real){
        Tf y_norm = y_real * boxdh_inv_.y;
        return y_norm;
    }
    //! @brief convert normalized Y length to real Y length.
    template <class Tf>
    Tf realYLen(const Tf &y_norm){
        Tf y_real = y_norm * boxdh_.y;
        return y_real;
    }
    //! @brief convert real Z length to normalized Z length.
    template <class Tf>
    Tf normZLen(const Tf &z_real){
        Tf z_norm = z_real * boxdh_inv_.z;
        return z_norm;
    }
    //! @brief convert normalized Z length to real Z length.
    template <class Tf>
    Tf realZLen(const Tf &z_norm){
        Tf z_real = z_norm * boxdh_.z;
        return z_real;
    }
}
