#pragma once


#include "kernelfunc.h"



const PS::F64 W_deld       = W_3d( s_deld );

class CalcHydroForce_3d{
public:
  void operator()(const EP* const epi, const PS::S32 Nip, const EP* const epj, const PS::S32 Njp, Hydro* const hydro){
    // Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> ES;
    for( PS::S32 i = 0; i < Nip; ++i ){
      hydro[i].clear();
      if (epi[i].itype > S_W) continue;//i->only water type

      PS::F64vec pro1(zero), pro2(zero), pro3(zero);
      const PS::F64 densi_sq = epi[i].dens * epi[i].dens;
      const PS::F64 densi_sq_inv = 1.0 / densi_sq;

      // const PS::F64vec sigi1(epi[i].sig.xx, epi[i].sig.xy);
      // const PS::F64vec sigi2(epi[i].sig.xy, epi[i].sig.yy);
      
      for( PS::S32 j = 0 ; j < Njp ; ++j ){
        const PS::F64vec dr = epi[i].pos - epj[j].pos;
        const PS::F64 mj_rhoj = epj[j].mass / epj[j].dens;
        const PS::F64 densj_sq = epj[j].dens * epj[j].dens;
        const PS::F64 densj_sq_inv = 1.0 / densj_sq;
        const PS::F64 deno = 1.0 /( densi_sq * densj_sq );
        const PS::F64 dr_sq = dr * dr;
        const PS::F64 s = sqrt( dr_sq )/ smth_len;
        const PS::F64vec gradW_dr = gradW_3d(dr, s);
        PS::F64vec dv = epj[j].vel - epi[i].vel;// dv = J -I !!!

        /********* boundary particle -> non_slip condition: dv **********/
        if ( epj[j].itype == WALL){
          PS::F64 dj = epj[j].disp.getMax();
          PS::F64 di = 0.0;
          if( fabs(epj[j].disp.x - 0.0) > 0.01*dx){
            di = fabs( fabs(epj[j].disp.x - epi[i].disp.x) - dj);
          }
          else if( fabs(epj[j].disp.y - 0.0) > 0.01*dx){
            di = fabs( fabs(epj[j].disp.y - epi[i].disp.y) - dj);
          }
          else{
            di = fabs( fabs(epj[j].disp.z - epi[i].disp.z) - dj);
          }
          PS::F64 betav = std::min(1.5 , 1.0 + dj / di);
          dv = - betav * epi[i].vel;
        }

        // if (epj[j].itype == UNDER){
        //   PS::F64 dj = -epj[j].pos.y;
        //   PS::F64 di = epi[i].pos.y;
        //   PS::F64 betav = std::min( 1.5 , 1.0 + dj / di);
        //   dv = - betav * epi[i].vel;
        // }
        // else if (epj[j].itype == LEFT){
        //   PS::F64 dj = -epj[j].pos.x;
        //   PS::F64 di = epi[i].pos.x;
        //   PS::F64 betav = std::min( 1.5 , 1.0 + dj / di);
        //   dv = - betav * epi[i].vel;
        // }
        // else if (epj[j].itype == RIGHT){
        //   PS::F64 dj = epj[j].pos.x - WallW;
        //   PS::F64 di = WallW - epi[i].pos.x;
        //   PS::F64 betav = std::min( 1.5 , 1.0 + dj / di);
        //   dv = - betav * epi[i].vel;
        // }
        // else if (epj[j].itype == FRONT){
        //   PS::F64 dj = - epj[j].pos.z;
        //   PS::F64 di = epi[i].pos.z;
        //   PS::F64 betav = std::min( 1.5 , 1.0 + dj / di);
        //   dv = - betav * epi[i].vel;
        // }
        // else if (epj[j].itype == BACK){
        //   PS::F64 dj = epj[j].pos.z - WallD;
        //   PS::F64 di = WallD - epi[i].pos.z;
        //   PS::F64 betav = std::min( 1.5 , 1.0 + dj / di);
        //   dv = - betav * epi[i].vel;
        // }

        /****** calc density *****/
        hydro[i].dens_rate -= epj[j].mass * (dv * gradW_dr);

        /***** calc gradv.<dv/dx> (-> eps_rate ) *****/
        pro1 += (mj_rhoj * dv.x) * gradW_dr;
        pro2 += (mj_rhoj * dv.y) * gradW_dr;
        pro3 += (mj_rhoj * dv.z) * gradW_dr;

        /***** motion equation *****/
        PS::F64mat motion(zero);
        PS::F64mat sigj = epj[j].sig;
        PS::F64mat Rj   = epj[j].R_tensor;
        if ( epj[j].itype > S_W ) sigj = epi[i].sig;// non_slip condition
        if ( epj[j].itype > S_W ) Rj   = epi[i].R_tensor;// non_slip condition

        motion += product_02_3d(densi_sq_inv, epi[i].sig) + product_02_3d(densj_sq_inv, sigj);
        // motion += epi[i].sig / densi_sq + epj[j].sig / densj_sq;

        /***** artificial viscousity ->acc *****/
        const PS::F64vec dvel = -dv;
        const PS::F64 vx = dr * dvel;
        const PS::F64 mrho = 0.5 * (epi[i].dens + epj[j].dens);
        PS::F64 phi_av = ( smth_len * vx )/( dr_sq + 0.01 * smth_len * smth_len );
        PS::F64 PIij = zero;
        if (vx < zero) PIij = ( multi_acij + beta * phi_av) * phi_av / mrho;
        PS::F64mat PI_delta(PIij, PIij, PIij, 0.0, 0.0, 0.0);

        motion -= PI_delta;

        hydro[i].acc -=  epj[j].mass * PIij * gradW_dr;

        /***** tensile instability -> artificial stress *****/
        PS::F64 fn = pow( W_3d( s )/ W_deld , n_AS );

        // motion += product_02_3d( fn, Ri + Rj );
        motion += product_02_3d( fn, epi[i].R_tensor + Rj );

        /***** sum motion eq. ******/
        hydro[i].acc += dot_product_21_3d(motion, gradW_dr) * epj[j].mass;
      }// end loop J

      /***** eps_rate *****/
      hydro[i].eps_rate.xx = pro1.x;
      hydro[i].eps_rate.yy = pro2.y;
      hydro[i].eps_rate.zz = pro3.z;
      hydro[i].eps_rate.xy = 0.5 * ( pro1.y + pro2.x );
      hydro[i].eps_rate.xz = 0.5 * ( pro1.z + pro3.x );
      hydro[i].eps_rate.yz = 0.5 * ( pro2.z + pro3.y );
      /***** calc spin tensor (omg_rate) *****/
      // hydro[i].omg_s = 0.5 * ( pro1.y - pro2.x );//yx = -xy
      hydro[i].omg_rate.x = hydro[i].eps_rate.xy - pro2.x;
      hydro[i].omg_rate.y = hydro[i].eps_rate.xz - pro3.x;
      hydro[i].omg_rate.z = hydro[i].eps_rate.yz - pro3.y;
      
      /***** gravity -> acc *****/
      hydro[i].acc.y -= gravity;
    }// end loop I
  }
};


