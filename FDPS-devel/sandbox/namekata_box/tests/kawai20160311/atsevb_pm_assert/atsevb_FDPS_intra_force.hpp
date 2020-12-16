//***************************************************************************************
//  This program is the intramolecular force calculation for "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include<unordered_map>

//--- template function for FDPS interface
template<class Tepi, class Tepj>
void CalcForceIntra(const Tepi* EPI,
                    const PS::S32 n_EPI,
                    const Tepj* EPJ,
                    const PS::S32 n_EPJ,
                          ForceIntra* force){
    
    //--- temporary data
    std::unordered_map<PS::S32, PS::S32> IDtable;
    IDtable.max_load_factor(0.7);
    IDtable.reserve( n_EPJ );
    
    //--- create temporary table for atomID and "j"
    for(PS::S32 j=0; j<n_EPJ; j++){
        IDtable[EPJ[j].getAtomID()] = j;
    }
    
    // check contents in IDtable
    std::unordered_map<PS::S32, PS::S32>::iterator it;
    std::cerr << " # check IDtable:" << std::endl;
    for(it = IDtable.begin(); it != IDtable.end(); it++){
        std::pair<PS::S32, PS::S32> element = *it;
        std::cerr << "   AtomID=" << element.first
                  << "  index=" << element.second << std::endl;
    }
    
    //--- calculate intramolecular force
    for(PS::S32 i=0; i<n_EPI; i++){
        //--- switch intraforce function by enum MOL_TYPE
        switch(EPI[i].getMolType()){
            case MOL_TYPE::solvent:
                SetWater::CalcForceIntra(EPI, i,
                                         EPJ, n_EPJ,
                                         IDtable,
                                         force);
                break;
                
            case MOL_TYPE::cluster:
                //if( EPI[i].EVB.getN() == 2 ){  // at excess proton
                //    SetEVB::CalcForceIntra(EPI, i,
                //                           EPJ, n_EPJ,
                //                           IDtable,
                //                           force);
                //}
                break;
            
            case MOL_TYPE::polymer:
                //SetNaf::CalcForceIntra(EPI, i,
                //                       EPJ, n_EPJ,
                //                       IDtable,
                //                       force);
                break;
            
            case MOL_TYPE::wall:
                // (do nothing)
                break;
            
            default:
                // (do nothing)
                break;
        }
    }
}