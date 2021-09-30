//***************************************************************************************
//  This program is the intermolecular interactuion of "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <algorithm>

//--- cutoff function for Particle-Particle Particle-Mesh (P3M)
//------ HPCI 戦略プログラム分野５全体シンポジウム (2012/3/7) 課題４報告 牧野淳一郎 P.30
//------ http://jun.artcompsci.org/talks/akihabara20120307a.pdf
inline PS::F64 cutfunc_p3m(const PS::F64 r){
    //--- 0.0 <= r <= 2.0 when scale = 1.0.
    const PS::F64 scale = 1.0;
    
    const PS::F64 c1 =   8.0/5.0;
    const PS::F64 c2 =  -1.0/2.0;
    const PS::F64 c3 = -12.0/35.0;
    const PS::F64 c4 =   3.0/20.0;
    const PS::F64 c5 =   3.0/35.0;
    const PS::F64 c6 =  18.0/35.0;
    const PS::F64 c7 =   1.0/5.0;
    
    PS::F64 tmp;
    
    //--- -1.0 <= tmp <= 1.0
    tmp = std::max(0.0, r - scale);
    tmp = tmp * tmp;         // square
    tmp = tmp * tmp * tmp;   // 6th powered
    
    tmp = 1.0 + r*r*r*(-c1 + r*r*(c1 + r*(c2 + r*(c3 + r*c4)))) - tmp*(c5 + r*(c6 + r*c7));
    return tmp;
}


//--- Force class
//------ This class has the result of force. (used as temporary data)
//------ This class is the subset of Full Particle class.

//------ LJ interaction with cutoff length (simple cutoff)
class ForceLJ {
public:
    PS::F64vec force_LJ;
    PS::F64      pot_LJ;
    
#ifdef VIRIAL_TEST
    PS::F64vec virt_LJ;
#endif
    
    void clear(){
        force_LJ = 0.0;
          pot_LJ = 0.0;
          
#ifdef VIRIAL_TEST
         virt_LJ = 0.0;
        // std::cout << "virial is cleared in ForceLJ" << std::endl << std::endl;
#endif
    }
};

//------ coulomb interaction
class ForceCoulomb {
public:
    PS::F64vec field_coulomb;
    PS::F64      pot_coulomb;
    
#ifdef VIRIAL_TEST
    PS::F64vec virt_coulomb;
#endif
    
    void clear(){
        field_coulomb = 0.0;
          pot_coulomb = 0.0;
          
#ifdef VIRIAL_TEST
         virt_coulomb = 0.0;
        // std::cout << "virial is cleared in ForceCoulomb" << std::endl << std::endl;
#endif
    }
};

//------ LJ and Coulomb-static in cutoff length
//--------- this function do not use SuperParticle.
class ForceScatter {
public:
    PS::F64vec force_LJ;
    PS::F64      pot_LJ;
    PS::F64vec field_coulomb;
    PS::F64      pot_coulomb;
    
#ifdef VIRIAL_TEST
    PS::F64vec virt_LJ;
    PS::F64vec virt_coulomb;
#endif
    
    void clear(){
        force_LJ      = 0.0;
          pot_LJ      = 0.0;
        field_coulomb = 0.0;
          pot_coulomb = 0.0;
          
#ifdef VIRIAL_TEST
         virt_LJ      = 0.0;
         virt_coulomb = 0.0;
        // std::cout << "virial is cleared in ForceScatter" << std::endl << std::endl;
#endif
    }
};

//--- Template function calculate Particle-Particle interaction
//------ for LJ interaction with simple cutoff
template<class Tepi, class Tepj>
void CalcForceLJ(const Tepi *EPI,
                 const PS::S32 n_EPI,
                 const Tepj *EPJ,
                 const PS::S32 n_EPJ,
                       ForceLJ *force){
    //--- setting for cutoff length
    PS::F64 R2_cut = EPI->getRSearch();
            R2_cut = R2_cut*R2_cut;
    
    for(PS::S32 i=0; i<n_EPI; i++){
        MOL_TYPE   i_type = EPI[i].getType();
        if( i_type == cluster ) continue;
        
        PS::S32    i_mol  = EPI[i].getMolID();
        PS::S32    i_atom = EPI[i].getAtomID();
        PS::F64vec pos_i  = EPI[i].getPos();
        PS::F64    vdw_d  = EPI[i].getVDW_D();
        PS::F64    vdw_r  = EPI[i].getVDW_R();
        PS::F64vec F_LJ_i        = 0.0;
        PS::F64    pot_LJ_i      = 0.0;
        
#ifdef VIRIAL_TEST
        PS::F64vec virt_LJ      = 0.0;
#endif
        
        for(PS::S32 j=0; j<n_EPJ; j++){
            //--- intramolecular mask
            if( EPJ[j].getType() == cluster ) continue;
            if( checkIntraMask(i_type, i_mol, i_atom, EPJ, j) ) continue;
            
            //--- intermolecular interaction
            PS::F64vec r_ij = pos_i - EPJ[j].pos;
                       realPos(r_ij);               // convert to real space
            PS::F64    r2   = r_ij*r_ij;
            if(r2 <= R2_cut){
                PS::F64 r      = sqrt(r2);
                
                //--- VDW part
                PS::F64 vddm   = EPJ[j].getVDW_D();
                PS::F64 vdrm   = EPJ[j].getVDW_R();
                        vddm   = sqrt(vddm*vdw_d);
                 //     vddm   = vddm*vdw_d                    // if VDW_D values are pre-affected "sqrt"
                        vdrm   = 0.5*(vdrm + vdw_r);
                        r2     = 1.0/r2;                       // using "r2" as r2_inv
                PS::F64 sbr6   = vdrm*vdrm*r2;                 // (r0/r)^2
                        sbr6   = sbr6*sbr6*sbr6;               // 6 powered, (r0/r)^6
                        vddm   = vddm*sbr6;                    // using "vddm" as vddm*(r0/r)^6
                        vdrm   = 12.0*vddm*(sbr6-1.0)*r2;      // using "vdrm" as (VDW force)/r
                        F_LJ_i   += vdrm*r_ij;
                        pot_LJ_i += vddm*(sbr6-2.0);
                        
#ifdef VIRIAL_TEST
                        virt_LJ.x  += 0.5*vdrm*r_ij.x*r_ij.x;
                        virt_LJ.y  += 0.5*vdrm*r_ij.y*r_ij.y;
                        virt_LJ.z  += 0.5*vdrm*r_ij.z*r_ij.z;
#endif
                        
            }
        }
        //--- writeback result of i_atom
        force[i].force_LJ      += F_LJ_i;
        force[i].pot_LJ        += pot_LJ_i;
        
#ifdef VIRIAL_TEST
        force[i].virt_LJ       += virt_LJ;
#endif
    }
}


//------ for coulomb interaction (Permanent charge part)
template<class Tepi, class Tepj>
void CalcForceCoulomb_perm(const Tepi *EPI,
                           const PS::S32 n_EPI,
                           const Tepj *EPJ,
                           const PS::S32 n_EPJ,
                                 ForceCoulomb *force){
    PS::F64 r_scale = EPI->getRSearch();
            realXLen(r_scale);              // convert to real space (cubic system)
            r_scale = 2.0/r_scale;          // r_scale = 1.0/(0.5*r_cut)
    
    for(PS::S32 i=0; i<n_EPI; i++){
        MOL_TYPE   i_type = EPI[i].getType();
        if( i_type == cluster ) continue;
        
        PS::S32    i_mol  = EPI[i].getMolID();
        PS::S32    i_atom = EPI[i].getAtomID();
        PS::F64vec pos_i  = EPI[i].getPos();
        PS::F64vec F_coulomb_i   = 0.0;
        PS::F64    pot_coulomb_i = 0.0;
        
#ifdef VIRIAL_TEST
        PS::F64    i_charge = EPI[i].getCharge();
        PS::F64vec virt_coulomb = 0.0;
#endif
        
        for(PS::S32 j=0; j<n_EPJ; j++){
            //--- intramolecular mask
            if( EPJ[j].getType() == cluster ) continue;
            if( checkIntraMask(i_type, i_mol, i_atom, EPJ, j) ) continue;
            
            //--- intermolecular interaction
            PS::F64vec r_ij = pos_i - EPJ[j].pos;
                       realPos(r_ij);                         // convert to real space
            PS::F64    r2   = r_ij*r_ij;
            PS::F64    r    = sqrt(r2);
                       r2   = 1.0/r2;                         // r2 means "r2_inv"
                    
            //--- coulomb static part
            PS::F64   g_p3m = cutfunc_p3m(r*r_scale);         // g_p3m means the value of cutoff function
                      r     = g_p3m*r*r2*EPJ[j].getCharge();  // r means "qj/r" with cutoff
                      r2    = r2*r;                           // r2 means "qj/(r^3)" with cutoff
                      F_coulomb_i   += r2*r_ij;
                      pot_coulomb_i += r;
            
#ifdef VIRIAL_TEST
                      virt_coulomb.x += 0.5*r2*i_charge*r_ij.x*r_ij.x;
                      virt_coulomb.y += 0.5*r2*i_charge*r_ij.y*r_ij.y;
                      virt_coulomb.z += 0.5*r2*i_charge*r_ij.z*r_ij.z;
#endif
        }
        //--- writeback result of i_atom
        force[i].field_coulomb += F_coulomb_i;
        force[i].pot_coulomb   += pot_coulomb_i;
        
#ifdef VIRIAL_TEST
        force[i].virt_coulomb  += virt_coulomb;
#endif
    }
}

//--- calculate Particle-Particle interaction in cutoff distance
template<class Tepi, class Tepj>
void CalcForceScatter(const Tepi *EPI,
                      const PS::S32 n_EPI,
                      const Tepj *EPJ,
                      const PS::S32 n_EPJ,
                            ForceScatter *force){
    PS::F64 r_scale = EPI->getRSearch();
    PS::F64 R2_cut = r_scale*r_scale;
            r_scale = 2.0/r_scale;   // r_scale = 1.0/(0.5*r_cut)
    
    for(PS::S32 i=0; i<n_EPI; i++){
        MOL_TYPE   i_type = EPI[i].getType();
        if( i_type == cluster ) continue;
        
        PS::S32    i_mol  = EPI[i].getMolID();
        PS::S32    i_atom = EPI[i].getAtomID();
        PS::F64vec pos_i  = EPI[i].getPos();
        PS::F64    vdw_d  = EPI[i].getVDW_D();
        PS::F64    vdw_r  = EPI[i].getVDW_R();
        PS::F64vec F_LJ_i        = 0.0;
        PS::F64    pot_LJ_i      = 0.0;
        PS::F64vec F_coulomb_i   = 0.0;
        PS::F64    pot_coulomb_i = 0.0;
        
#ifdef VIRIAL_TEST
        PS::F64    i_charge = EPI[i].getCharge();
        PS::F64vec virt_LJ      = 0.0;
        PS::F64vec virt_coulomb = 0.0;
#endif
        
        for(PS::S32 j=0; j<n_EPJ; j++){
            //--- intramolecular mask
            if( EPJ[j].getType() == cluster ) continue;
            if( checkIntraMask(i_type, i_mol, i_atom, EPJ, j) ) continue;
            
            //--- intermolecular interaction
            PS::F64vec r_ij = pos_i - EPJ[j].pos;
            PS::F64    r2   = r_ij*r_ij;
            if(r2 <= R2_cut){
                PS::F64 r      = sqrt(r2);
                
                //--- VDW part
                PS::F64 vddm   = EPJ[j].getVDW_D();
                PS::F64 vdrm   = EPJ[j].getVDW_R();
                        vddm   = sqrt(vddm*vdw_d);
                 //     vddm   = vddm*vdw_d                    // if VDW_D values are pre-affected "sqrt"
                        vdrm   = 0.5*(vdrm + vdw_r);
                        r2     = 1.0/r2;                       // using "r2" as r2_inv
                PS::F64 sbr6   = vdrm*vdrm*r2;                 // (r0/r)^2
                        sbr6   = sbr6*sbr6*sbr6;               // 6 powered, (r0/r)^6
                        vddm   = vddm*sbr6;                    // using "vddm" as vddm*(r0/r)^6
                        vdrm   = 12.0*vddm*(sbr6-1.0)*r2;      // using "vdrm" as (VDW force)/r
                        F_LJ_i   += vdrm*r_ij;
                        pot_LJ_i += vddm*(sbr6-2.0);
                        
#ifdef VIRIAL_TEST
                        virt_LJ.x  += 0.5*vdrm*r_ij.x*r_ij.x;
                        virt_LJ.y  += 0.5*vdrm*r_ij.y*r_ij.y;
                        virt_LJ.z  += 0.5*vdrm*r_ij.z*r_ij.z;
#endif
                        
                //--- coulomb static part
                        sbr6  = cutfunc_p3m(r*r_scale);        // sbr6 means the value of cutoff function
                        r     = sbr6*r*r2*EPJ[j].getCharge();  // r means "qj/r" with cutoff
                        r2    = r2*r;                          // r2 means "qj/(r^3)" with cutoff
                        F_coulomb_i   += r2*r_ij;
                        pot_coulomb_i += r;
                
                    //    r_inv = r_inv*EPJ[j].getCharge(); // r_inv means "qj/r"
                    //    r2    = r2*r_inv;                 // r2 means "qj/(r^3)"
                    //    F_coulomb_i   += r2*r_ij;
                    //    pot_coulomb_i += r_inv;
                
#ifdef VIRIAL_TEST
                        virt_coulomb.x += 0.5*r2*i_charge*r_ij.x*r_ij.x;
                        virt_coulomb.y += 0.5*r2*i_charge*r_ij.y*r_ij.y;
                        virt_coulomb.z += 0.5*r2*i_charge*r_ij.z*r_ij.z;
#endif
            }
        }
        //--- writeback result of i_atom
        force[i].force_LJ      += F_LJ_i;
        force[i].pot_LJ        += pot_LJ_i;
        force[i].field_coulomb += F_coulomb_i;
        force[i].pot_coulomb   += pot_coulomb_i;
        
#ifdef VIRIAL_TEST
        force[i].virt_LJ       += virt_LJ;
        force[i].virt_coulomb  += virt_coulomb;
        
        //std::cout << "in calcforce" << std::endl;
        //std::cout << "  Atom_ID: " << i_atom << std::endl
        //          << "  virt_LJ     :" << virt_LJ << std::endl
        //          << "  virt_coulomb:" << virt_coulomb <<std::endl;
#endif
    }
}

//--- comparison 2 result of virial
template<class Tfp>
void checkVirial_inter(const Tfp &fp,
                       const PS::S32 n_loc){
    PS::F64vec virt_LJ_atom;
    PS::F64vec virt_coulomb_atom;
    PS::F64vec virt_LJ_total;
    PS::F64vec virt_coulomb_total;
    
    PS::F64vec virt_LJ_atom_loc       = 0.0;
    PS::F64vec virt_coulomb_atom_loc  = 0.0;
    PS::F64vec virt_LJ_total_loc      = 0.0;
    PS::F64vec virt_coulomb_total_loc = 0.0;
    
    
    //--- sumation in local process
    for(PS::S32 i=0; i<n_loc; i++){
        virt_LJ_atom_loc       += 0.5 * fp[i].virt_LJ;
        virt_coulomb_atom_loc  += 0.5 * fp[i].virt_coulomb;
        
        PS::F64vec pos_i = fp[i].pos;
                   realPos(pos_i);     // convert to real space
        
        virt_LJ_total_loc.x      += 0.5 * pos_i.x * fp[i].force_LJ.x;
        virt_LJ_total_loc.y      += 0.5 * pos_i.y * fp[i].force_LJ.y;
        virt_LJ_total_loc.z      += 0.5 * pos_i.z * fp[i].force_LJ.z;
        virt_coulomb_total_loc.x += 0.5 * pos_i.x * fp[i].force_coulomb.x;
        virt_coulomb_total_loc.y += 0.5 * pos_i.y * fp[i].force_coulomb.y;
        virt_coulomb_total_loc.z += 0.5 * pos_i.z * fp[i].force_coulomb.z;
    }
    
    //--- sumation in all process
    virt_LJ_atom       = PS::Comm::getSum(virt_LJ_atom_loc);
    virt_coulomb_atom  = PS::Comm::getSum(virt_coulomb_atom_loc);
    virt_LJ_total      = PS::Comm::getSum(virt_LJ_total_loc);
    virt_coulomb_total = PS::Comm::getSum(virt_coulomb_total_loc);
    
    //--- show result
    //------ for LJ
    PS::F64vec  err;
    PS::F64vec  r_err;
    err   = virt_LJ_total - virt_LJ_atom;
    r_err = 0.0;
    if(virt_LJ_total.x != 0.0) r_err.x = err.x/virt_LJ_total.x;
    if(virt_LJ_total.y != 0.0) r_err.y = err.y/virt_LJ_total.y;
    if(virt_LJ_total.z != 0.0) r_err.z = err.z/virt_LJ_total.z;
    if(PS::Comm::getRank() == 0){
        std::cout << "virt_LJ  error: " << err << std::endl;
        std::cout << "     rel error: " << r_err << std::endl;
    }
    //------ for coulomb
    err   = virt_coulomb_total - virt_coulomb_atom;
    r_err = 0.0;
    if(virt_coulomb_total.x != 0.0) r_err.x = err.x/virt_coulomb_total.x;
    if(virt_coulomb_total.y != 0.0) r_err.y = err.y/virt_coulomb_total.y;
    if(virt_coulomb_total.z != 0.0) r_err.z = err.z/virt_coulomb_total.z;
    if(PS::Comm::getRank() == 0){
        std::cout << "virt_coulomb  error: " << err << std::endl;
        std::cout << "          rel error: " << r_err << std::endl;
    }
}

