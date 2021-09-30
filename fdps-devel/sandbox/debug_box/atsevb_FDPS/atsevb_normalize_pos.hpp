//***************************************************************************************
//  This program is the normalize interface of position data.
//    convert pos{0~boxdh_x, 0~boxdh_y, 0~boxdh_z} and pos{0~1, 0~1, 0~1}.
//    This code is used by "atsevb_main.cpp"
//***************************************************************************************
#pragma once

//--- system box size
PS::F64 boxdh_x;
PS::F64 boxdh_y;
PS::F64 boxdh_z;

//--- inverse value of system box size
PS::F64 boxdh_x_inv;
PS::F64 boxdh_y_inv;
PS::F64 boxdh_z_inv;

//--- mesh size in normalized spase for Particle Mesh
PS::F64 mesh_dx = 1.0/SIZE_OF_MESH;
PS::F64 mesh_dy = 1.0/SIZE_OF_MESH;
PS::F64 mesh_dz = 1.0/SIZE_OF_MESH;

void setNormalizePosition() {
    boxdh_x_inv = 1.0/boxdh_x;
    boxdh_y_inv = 1.0/boxdh_y;
    boxdh_z_inv = 1.0/boxdh_z;
}

inline void normPos(PS::F64vec &pos){
    pos.x = pos.x * boxdh_x_inv;
    pos.y = pos.y * boxdh_y_inv;
    pos.z = pos.z * boxdh_z_inv;
}

inline void realPos(PS::F64vec &pos){
    pos.x = pos.x * boxdh_x;
    pos.y = pos.y * boxdh_y;
    pos.z = pos.z * boxdh_z;
}


//--- use only cubic system (the "boxdh_x = boxdh_y = boxdh_z" case only)
inline void normXLen(PS::F64 &x){
    x = x*boxdh_x_inv;
}

inline void realXLen(PS::F64 &x){
    x = x*boxdh_x;
}

