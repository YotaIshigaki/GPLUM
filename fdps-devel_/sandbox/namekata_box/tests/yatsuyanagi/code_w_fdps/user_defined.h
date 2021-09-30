#pragma once
/* Standard headers */
#include <math.h>
/* FDPS headers */
#include "FDPS_c_if.h"

typedef struct full_particle { //$fdps FP,EPI,EPJ,Force
    //$fdps copyFromForce full_particle (vel,vel)
    //$fdps copyFromFP full_particle (id,id) (wz,wz) (pos,pos)
    //$fdps clear id=keep, wz=keep, pos=keep
    long long id; 	// $fdps id
    double  wz; 		// $fdps charge
    fdps_f64vec pos; // $fdps position
    fdps_f64vec vel; 	// $fdps acc
} Full_particle;

void  calc_BS_ep_ep(Full_particle *ep_i,
                         int n_ip,
                         Full_particle *ep_j,
                         int n_jp,
                         Full_particle *f);

void  calc_BS_ep_sp(Full_particle *ep_i,
                         int n_ip,
                         fdps_spj_monopole *ep_j,
                         int n_jp,
                         Full_particle *f);

