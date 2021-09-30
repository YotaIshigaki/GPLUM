//***************************************************************************************
//  This program is the normalize interface of position data.
//    convert pos{0~boxdh__x, 0~boxdh__y, 0~boxdh__z} and pos{0~1, 0~1, 0~1}.
//    This code is used by "atsevb_main.cpp"
//***************************************************************************************
#pragma once

#include<cmath>
#include<algorithm>

namespace Normalize {
    
    //--- unit rule
    class REAL_F32 {};
    class REAL_F32vec {};
    class REAL_F64 {};
    class REAL_F64vec {};
    class NORM_F32 {};
    class NORM_F32vec {};
    class NORM_F64 {};
    class NORM_F64vec {};
    
    //--- system box size
    PS::F64vec boxdh_ = (1.0, 1.0, 1.0);
    
    //--- inverse value of system box size
    PS::F64vec boxdh_inv_ = (1.0, 1.0, 1.0);
    
    //--- parameter for Particle Mesh
    //------ mesh size in normalized spase for Particle Mesh
    //------ "SIZE_OF_MESH" is defined in "$(PS_DIR)/src/particle_mesh/param_fdps.h".
    //------ recommended value of "SIZE_OF_MESH" is N^(1/3)/2. (N is the total number of particles)
    PS::F64vec coefPMForce = (1.0, 1.0, 1.0);
    PS::F64vec mesh_dSize = ( 1.0/SIZE_OF_MESH,
                              1.0/SIZE_OF_MESH,
                              1.0/SIZE_OF_MESH );
    
    //--- internal functions
    void setNormalizePosition_() {
        boxdh_inv_[0] = 1.0/boxdh_[0];
        boxdh_inv_[1] = 1.0/boxdh_[1];
        boxdh_inv_[2] = 1.0/boxdh_[2];
    }
    void setCoefPM_(){
        coefPMForce[0] = 1.0*(boxdh_inv_[0]*boxdh_inv_[0]);
        coefPMForce[1] = 1.0*(boxdh_inv_[1]*boxdh_inv_[1]);
        coefPMForce[2] = 1.0*(boxdh_inv_[2]*boxdh_inv_[2]); // squared
    }
    
    //--- access functions (must be call only below functions)
    //------ initialize or setting the real domain size
    void setBoxSize(const PS::F64vec & box){
        boxdh_ = box;
        setNormalizePosition_();
        setCoefPM_();
    }
    
    //------ convert real pos to norm pos
    inline PS::F64vec normPos(const PS::F64vec & pos){
        PS::F64vec tmp = 0.0;
        tmp[0] = pos[0] * boxdh_inv_[0];
        tmp[1] = pos[1] * boxdh_inv_[1];
        tmp[2] = pos[2] * boxdh_inv_[2];
        return tmp;
    }
    //------ convert norm pos to real pos
    inline PS::F64vec realPos(const PS::F64vec & pos){
        PS::F64vec tmp = 0.0;
        tmp[0] = pos[0] * boxdh_[0];
        tmp[1] = pos[1] * boxdh_[1];
        tmp[2] = pos[2] * boxdh_[2];
        return tmp;
    }
    //------ convert norm PM force to real PM force
    inline PS::F64vec realPMForce(const PS::F64vec & Fpm){
        PS::F64vec tmp   = 0.0;
                   tmp[0] = Fpm[0] * coefPMForce[0];
                   tmp[1] = Fpm[1] * coefPMForce[1];
                   tmp[2] = Fpm[2] * coefPMForce[2];
        return tmp;
    }
    //------ convert real cutoff length to norm cutoff length
    inline PS::F64 normCutOff(const PS::F64 & realCutOff){
        PS::F64 len_inv = std::max(boxdh_inv_[0], boxdh_inv_[1]);
                len_inv = std::max(len_inv      , boxdh_inv_[2]);
        return realCutOff*len_inv;
    }
    //------ convert norm cutoff length to real cutoff length
    inline PS::F64 realCutOff(const PS::F64 & normCutOff){
        PS::F64 len = std::min(boxdh_[0], boxdh_[1]);
                len = std::min(len      , boxdh_[2]);
        return normCutOff*len;
    }
    
    //------ adjustment position into normalized periodic boundary system
    inline PS::F64vec periodicAdjust(const PS::F64vec & pos_norm){
        PS::F64vec pos_new = pos_norm;
        //* original
        //if(pos_new[0] < 0.0)  pos_new[0] += 1.0;
        //if(pos_new[1] < 0.0)  pos_new[1] += 1.0;
        //if(pos_new[2] < 0.0)  pos_new[2] += 1.0;
        //if(pos_new[0] >= 1.0) pos_new[0] -= 1.0;
        //if(pos_new[1] >= 1.0) pos_new[1] -= 1.0;
        //if(pos_new[2] >= 1.0) pos_new[2] -= 1.0;
        //* FDPS team
        if(pos_new.x < 0.0)  pos_new.x += 1.0;
        if(pos_new.y < 0.0)  pos_new.y += 1.0;
        if(pos_new.z < 0.0)  pos_new.z += 1.0;
        if(pos_new.x >= 1.0) pos_new.x -= 1.0;
        if(pos_new.y >= 1.0) pos_new.y -= 1.0;
        if(pos_new.z >= 1.0) pos_new.z -= 1.0;
      //  pos_new[0] -= round(pos_new[0] - 0.5000000000001);
      //  pos_new[1] -= round(pos_new[1] - 0.5000000000001);
      //  pos_new[2] -= round(pos_new[2] - 0.5000000000001);
        return pos_new;
    }
    
    //------ drift in normalized space
    inline PS::F64vec normDrift(const PS::F64vec & move){
        return normPos(move);
    }
    
    //------ length converter
    inline PS::F64 normXLen(const PS::F64 & x_real){
        PS::F64 x_norm = x_real * boxdh_inv_[0];
        return x_norm;
    }
    inline PS::F64 realXLen(const PS::F64 & x_norm){
        PS::F64 x_real = x_norm * boxdh_[0];
        return x_real;
    }
    inline PS::F64 normYLen(const PS::F64 & y_real){
        PS::F64 y_norm = y_real * boxdh_inv_[1];
        return y_norm;
    }
    inline PS::F64 realYLen(const PS::F64 & y_norm){
        PS::F64 y_real = y_norm * boxdh_[1];
        return y_real;
    }
    inline PS::F64 normZLen(const PS::F64 & z_real){
        PS::F64 z_norm = z_real * boxdh_inv_[2];
        return z_norm;
    }
    inline PS::F64 realZLen(const PS::F64 & z_norm){
        PS::F64 z_real = z_norm * boxdh_[2];
        return z_real;
    }
}
