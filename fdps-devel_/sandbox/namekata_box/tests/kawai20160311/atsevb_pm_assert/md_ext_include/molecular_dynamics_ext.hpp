//***************************************************************************************
//  This file is include interface for #include<md_ext.hpp>
//  This code is made for Molecular Dynamics (MD) simulation on 
//    Framework for Developping Particle Simulator (FDPS)
//***************************************************************************************
#pragma once

#include "md_ext_normalize.hpp"
#include "md_ext_intra_mask_list.hpp"
#include "md_ext_multi_table.hpp"


namespace MD_EXT {
    //--- presets for MultiTable
    //------ 32bit integer is 9(+1) digit in decimal system.
    //--------- 32bit int key, max 4 key types, key window: 0~99.
    template<class Tvalue, typename ...Tkey>
    using MultiTable_4d2 = MultiTable_base_<int_fast32_t, 4, 100, Tvalue, Tkey...>;
    //--------- 32bit int key, max 3 key types, key window: 0~999.
    template<class Tvalue, typename ...Tkey>
    using MultiTable_3d3 = MultiTable_base_<int_fast32_t, 3, 1000, Tvalue, Tkey...>;
    
    //------ 64bit integer is 18(+1) digit in decimal system.
    //--------- 64bit int key, max 9 key types, key window: 0~99.
    template<class Tvalue, typename ...Tkey>
    using MultiTable_9d2 = MultiTable_base_<int64_t, 9, 100, Tvalue, Tkey...>;
    //--------- 64bit int key, max 6 key types, key window: 0~999.
    template<class Tvalue, typename ...Tkey>
    using MultiTable_6d3 = MultiTable_base_<int64_t, 6, 1000, Tvalue, Tkey...>;
    
    //--- default
    template<class Tvalue, typename ...Tkey>
    using MultiTable = MultiTable_4d2<Tvalue, Tkey...>;
}


