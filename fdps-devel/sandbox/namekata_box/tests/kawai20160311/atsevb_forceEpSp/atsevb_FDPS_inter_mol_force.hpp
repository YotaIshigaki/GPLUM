//***************************************************************************************
//  This program is the intermolecular interactuion of "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include<algorithm>
#include<cassert>

//--- cutoff function for Particle-Particle Particle-Mesh (P3M)
//------ HPCI 戦略プログラム分野５全体シンポジウム (2012/3/7) 課題４報告 牧野淳一郎 P.30
//------ http://jun.artcompsci.org/talks/akihabara20120307a.pdf
inline PS::F64 cutfunc_p3m(const PS::F64 r,
                           const bool comp_flag = false){
    //--- 0.0 <= r <= 2.0 when scale = 1.0.
    assert( r >= 0.0 );
    if(r > 2.0) return 0.0;
    
    constexpr PS::F64 scale = 1.0;
    
    constexpr PS::F64 c1 =   8.0/5.0;
    constexpr PS::F64 c2 =  -1.0/2.0;
    constexpr PS::F64 c3 = -12.0/35.0;
    constexpr PS::F64 c4 =   3.0/20.0;
    constexpr PS::F64 c5 =   3.0/35.0;
    constexpr PS::F64 c6 =  18.0/35.0;
    constexpr PS::F64 c7 =   1.0/5.0;
    
    PS::F64 tmp;
    
    //--- -1.0 <= tmp <= 1.0
    tmp = std::max(0.0, r - scale);
    tmp = tmp * tmp;         // square
    tmp = tmp * tmp * tmp;   // 6th powered
    
    PS::F64 r2 = r*r;
    
    //--- complementaly value
    tmp = - r2*r*(-c1 + r2*(c1 + r*(c2 + r*(c3 + r*c4)))) + tmp*(c5 + r*(c6 + r*c7));
    
    if(comp_flag){
        //--- complementaly output
        return tmp;
    } else {
        //--- normal output
        return 1.0 - tmp;
    }
}


//--- intraList check function
//------ return true  : the pair is intramolecular.
//------ return false : the pair is intermolecular.
template<class Tepi, class Tepj>
inline bool isIntraPair(const Tepi *epi,
                        const PS::S32 & i,
                        const Tepj *epj,
                        const PS::S32 & j) {
    
    MOL_TYPE i_type = epi[i].getMolType();
    
    //--- check same molecule or not
    if(i_type            != epj[j].getMolType() ) return false;
    if(epi[i].getMolID() != epj[j].getMolID()   ) return false;
    
    //--- check intramolecular force pair or not
    switch(i_type){
        //------ in the case of solvent type
        case MOL_TYPE::solvent:
            //std::cerr << "    ID:" << epi[i].getAtomID() << " and ID:" << epj[j].getAtomID() 
            //          << " is intra pair." << std::endl;
            return true;
        
        //------ in the case of polymer type
        case MOL_TYPE::polymer: {
            if( epi[i].checkIntraList( epj[j].getAtomID() ) ) return true;
            return false;
        }
        
        //------ in the case of wall type
        case MOL_TYPE::wall: {
            if( epi[i].checkIntraList( epj[j].getAtomID() ) ) return true;
            return false;
        }
        
        default:
        #ifdef INTRA_MASK_TEST
            //--- error: definition of intramolecular pair was not found!
            std::cerr << "ERROR: definition of intramolecular pair was not found!" << std::endl
                      << "    check in the `isIntraPair()` function "
                      <<     "and `intraList` data." << std::endl
                      << "    i_type : " << static_cast<int>(i_type)               << std::endl
                      << "    i_mol  : " << epi[i].getMolID()                      << std::endl
                      << "    i_atom : " << epi[i].getAtomID()                     << std::endl
                      << "    j_type : " << static_cast<int>(epj[j].getMolType() ) << std::endl
                      << "    j_mol  : " << epj[j].getMolID()                      << std::endl
                      << "    j_atom : " << epj[j].getAtomID()                     << std::endl;
            
            throw std::invalid_argument("the pair was undefined.");
        #else       
            //--- in release version, non-defined pair error is ignored.
            return false; 
        #endif
    }
}


//--- simple functions
//------ culculate virial value of particle i
inline PS::F64vec pi_virial(const PS::F64vec pos, const PS::F64vec force){
    PS::F64vec tmp = 0.0;
    tmp.x = 0.5 * pos.x * force.x;
    tmp.y = 0.5 * pos.y * force.y;
    tmp.z = 0.5 * pos.z * force.z;
    return tmp;
}


//--- Template function calculate Particle-Particle interaction
//------ for LJ interaction with simple cutoff
template<class Tepi, class Tepj>
void CalcForceLJ(const Tepi *EPI,
                 const PS::S32 n_EPI,
                 const Tepj *EPJ,
                 const PS::S32 n_EPJ,
                       ForceLJ *force){
    //--- setting for cutoff length
    PS::F64 R2_cut = Normalize::realCutOff( EPI->getRSearch() );  // convert to real space
            R2_cut = R2_cut*R2_cut;
    
    for(PS::S32 i=0; i<n_EPI; i++){
        PS::S32    i_atom = EPI[i].getAtomID();
        PS::S32    i_mol  = EPI[i].getMolID();
        MOL_TYPE   i_type = EPI[i].getMolType();
        PS::F64vec pos_i  = EPI[i].getPos();
        PS::F64    vdw_d  = EPI[i].getVDW_D();
        PS::F64    vdw_r  = EPI[i].getVDW_R();
        PS::F64vec F_LJ_i    = 0.0;
        PS::F64vec virt_LJ_i = 0.0;
        PS::F64    pot_LJ_i  = 0.0;
        
        for(PS::S32 j=0; j<n_EPJ; j++){
            //--- intramolecular mask
            if( i_atom == EPJ[j].getAtomID() ) continue;
            if( isIntraPair(EPI, i, EPJ, j)  ) continue;
            
            //--- intermolecular interaction
            PS::F64vec r_ij = pos_i - EPJ[j].getPos();
                       r_ij = Normalize::realPos(r_ij);        // convert to real space
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
                        F_LJ_i    += vdrm*r_ij;
                        pot_LJ_i  += 0.5*vddm*(sbr6-2.0);      // for double count
                        virt_LJ_i += pi_virial(r_ij, vdrm*r_ij);
            }
        }
        //--- writeback result of i_atom
        force[i].addForceLJ(F_LJ_i);
        force[i].addPotLJ(pot_LJ_i);
        force[i].addVirtLJ(virt_LJ_i);
    }
}


//------ for coulomb interaction (Permanent charge part)
template<class Tepi, class Tepj>
void CalcForceCoulomb_perm(const Tepi *EPI,
                           const PS::S32 n_EPI,
                           const Tepj *EPJ,
                           const PS::S32 n_EPJ,
                                 ForceCoulomb *force){
    
  //  //--- real space version
  //  PS::F64 r_scale = Normalize::realCutOff( EPI->getRSearch() );  // convert to real space (cubic system)
  //          r_scale = 2.0/r_scale;                                 // r_scale = 1.0/(0.5*r_cut)
  //  
  //  for(PS::S32 i=0; i<n_EPI; i++){
  //      PS::S32    i_atom = EPI[i].getAtomID();
  //      PS::S32    i_mol  = EPI[i].getMolID();
  //      MOL_TYPE   i_type = EPI[i].getMolType();
  //      PS::F64vec pos_i  = EPI[i].getPos();
  //      PS::F64vec F_coulomb_i   = 0.0;
  //  //    PS::F64    pot_coulomb_i = 0.0;
  //      
  //      for(PS::S32 j=0; j<n_EPJ; j++){
  //          //--- the case of same atom
  //          if( i_atom == EPJ[j].getAtomID() ) continue;
  //          
  //          //--- distance between EPI and EPJ
  //          PS::F64vec r_ij = pos_i - EPJ[j].getPos();
  //                     r_ij = Normalize::realPos(r_ij);       // convert to real space
  //          PS::F64    r2   = r_ij*r_ij;
  //          PS::F64    r    = sqrt(r2);
  //                     r2   = 1.0/r2;                         // r2 means "r2_inv"
  //                     
  //          //--- cutoff function
  //          PS::F64 g_p3m = cutfunc_p3m(r*r_scale, isIntraPair(EPI, i, EPJ, j) );
  //                                      
  //          //--- coulomb static part
  //              r     = g_p3m*r*r2*EPJ[j].getCharge();  // r  means "qj/r" with cutoff
  //              r2    = r2*r;                           // r2 means "qj/(r^3)" with cutoff
  //              F_coulomb_i   += r2*r_ij;
  //            //  pot_coulomb_i += r;
  //
  //      }
  //      //--- writeback result of i_atom
  //      force[i].addFieldCoulomb(  F_coulomb_i);
  //    //  force[i].addPotCoulomb(  pot_coulomb_i);
  //  }
    
    //--- normalized space version
    PS::F64 r_scale = 2.0/( EPI->getRSearch() );   // r_scale = 1.0/(0.5*r_cut)
    
    for(PS::S32 i=0; i<n_EPI; i++){
        PS::S32    i_atom = EPI[i].getAtomID();
        PS::S32    i_mol  = EPI[i].getMolID();
        MOL_TYPE   i_type = EPI[i].getMolType();
        PS::F64vec i_pos  = EPI[i].getPos();
        PS::F64vec i_field = 0.0;
        
        for(PS::S32 j=0; j<n_EPJ; j++){
            //--- skip same atom
            if( EPJ[j].isSameAtom(i_atom) ) continue;
            
            //--- distance
            PS::F64vec r_ij = i_pos - EPJ[j].getPos();
            PS::F64    r2   = r_ij*r_ij;
            PS::F64    r    = sqrt(r2);
            PS::F64    r2_inv = 1.0/r2;
            
            //--- cufoff funciton
        //    std::cout << "r       =" << r << endl
        //              << "r_scale =" << r_scale << endl
        //              << "R       =" << r*r_scale << endl;
        //    bool comp_flag = isIntraPair(EPI, i, EPJ, j);
        //    std::cout << "EPI=" << EPI[i].getAtomID() << " and EPJ=" << EPJ[j].getAtomID()
        //              << " is " << comp_flag << std::endl;
        //    PS::F64 g_p3m = cutfunc_p3m(r*r_scale, comp_flag );
            //PS::F64 g_p3m = cutfunc_p3m(r*r_scale, isIntraPair(EPI, i, EPJ, j) );
            PS::F64 g_p3m = 1.0;
            
            //--- coulomb permanent part
            PS::F64 qj_r3    = g_p3m*r*r2_inv*r2_inv*EPJ[j].getCharge(); // qj/r^3
                    i_field += qj_r3*r_ij;
        }
        //--- writeback result for i_atom
        force[i].addFieldCoulomb( Normalize::realPMForce(i_field) );
    }
}

template<class Tepi>
void CalcForceCoulomb_EpSp(const Tepi *EPI,
                           const PS::S32 n_EPI,
                           const PS::SPJMonopoleCutoff *SPJ,
                           const PS::S32 n_EPJ,
                                 ForceCoulomb *force){
    
    //--- normalized space version
    PS::F64 r_scale = 2.0/( EPI->getRSearch() );   // r_scale = 1.0/(0.5*r_cut)
    
    for(PS::S32 i=0; i<n_EPI; i++){
        PS::S32    i_atom = EPI[i].getAtomID();
        PS::S32    i_mol  = EPI[i].getMolID();
        MOL_TYPE   i_type = EPI[i].getMolType();
        PS::F64vec i_pos  = EPI[i].getPos();
        PS::F64vec i_field = 0.0;
        
        for(PS::S32 j=0; j<n_EPJ; j++){
            //--- distance
            PS::F64vec r_ij = i_pos - SPJ[j].getPos();
            PS::F64    r2   = r_ij*r_ij;
            PS::F64    r    = sqrt(r2);
            PS::F64    r2_inv = 1.0/r2;
            
            //--- cufoff funciton
            //PS::F64 g_p3m = cutfunc_p3m(r*r_scale, isIntraPair(EPI, i, SPJ, j) );
            PS::F64 g_p3m = 1.0;
            
            //--- coulomb permanent part
            PS::F64 qj_r3    = g_p3m*r*r2_inv*r2_inv*SPJ[j].getCharge(); // qj/r^3
                    i_field += qj_r3*r_ij;
        }
        //--- writeback result for i_atom
        force[i].addFieldCoulomb( Normalize::realPMForce(i_field) );
    }
}

