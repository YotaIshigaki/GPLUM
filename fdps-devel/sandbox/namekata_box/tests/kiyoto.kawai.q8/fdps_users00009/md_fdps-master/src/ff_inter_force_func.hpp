//***************************************************************************************
//  This program is the intermolecular interactuion of "md_fdps_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

namespace FORCE {

    //--- Cutoff functions  (copy from FDPS-master/sample/c++/p3m/main.cpp)
    inline PS::F64 S2_pcut(const PS::F64 xi) {
       // This is the potential cutoff function where we used Eq.(8.75)
       // in Hockney & Eastwood (1987).

       if (xi <= 1.0) {
          return 1.0 - xi*(208.0
                          +(xi*xi)*(-112.0
                                   +(xi*xi)*(56.0
                                            +xi*(-14.0
                                                +xi*(-8.0
                                                    +3.0*xi)))))/140.0;
       } else if ((1.0 < xi) && (xi < 2.0)) {
          return 1.0 - (12.0
                       +xi*(128.0
                           +xi*(224.0
                               +xi*(-448.0
                                   +xi*(280.0
                                       +xi*(-56.0
                                           +xi*(-14.0
                                               +xi*(8.0
                                                   -xi))))))))/140.0;
       } else {
          return 0.0;
       }
    }

    inline PS::F64 S2_fcut(const PS::F64 xi) {
       // This function returns 1 - R(\xi), where \xi is r/(a/2), a is the
       // scale length of the cutoff function, and R(\xi) is almost the same
       // as the function defined as Eq.(8-72) in Hockney & Eastwood (1987).
       // The only difference is that [1/(r/(a/2))]^2 is factored out
       // in this function from Eq.(8-72).

       if (xi <= 1.0) {
          return 1.0 - (xi*xi*xi)*(224.0
                                  +(xi*xi)*(-224.0
                                           +xi*(70.0
                                               +xi*(48.0-21.0*xi))))/140.0;
       } else if ((1.0 < xi) && (xi < 2.0)) {
          return 1.0 - (12.0
                       +(xi*xi)*(-224.0
                                +xi*(896.0
                                    +xi*(-840.0
                                        +xi*(224.0
                                            +xi*(70.0
                                                +xi*(-48.0+7.0*xi)))))))/140.0;
       } else {
          return 0.0;
       }
    }

    //--- simple functions
    //------ culculate virial value of particle i
    inline PS::F64vec calcVirialEPI(const PS::F64vec &pos, const PS::F64vec &force){
        PS::F64vec tmp = 0.0;
        tmp.x = 0.5 * pos.x * force.x;
        tmp.y = 0.5 * pos.y * force.y;
        tmp.z = 0.5 * pos.z * force.z;
        return tmp;
    }

}
