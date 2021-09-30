//***************************************************************************************
//  This program is the intermolecular mask list for "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

//--- mask data class for intramolecular interaction pairs
//------ this data must be shared by manually.
class intraMask{
private:
    //--- member variables
    PS::S32 *mask_list;    // pointer for mask_list array
    PS::S32 *mask_index;   // pointer for mask_index array
    PS::S32  list_len;
    PS::S32  index_len;
    
public:
    //--- constructor: default.
    //--- initializer
    void initialize(const PS::S32 len, const PS::S32 n_atom){
        list_len  = len;
        index_len = n_atom + 1;
        mask_list  = new PS::S32[list_len];    // make mask_list array
        mask_index = new PS::S32[index_len];   // make mask_index array
        if(mask_list  == NULL) failure_make_list();
        if(mask_index == NULL) failure_make_list();
    }
    //--- destructor
    ~intraMask(){
        if(mask_list  != NULL) delete [] mask_list;
        if(mask_index != NULL) delete [] mask_index;
    }
    
    //--- member functions
    void setList(const PS::S32 index, const PS::S32 j_atom){
        this->mask_list[index] = j_atom;
    }
    void setIndex(const PS::S32 i_atom, const PS::S32 index){
        this->mask_index[i_atom] = index;
    }
    PS::S32 mask_st(const PS::S32 i_atom) const {
        return mask_index[i_atom];
    }
    PS::S32 mask_ed(const PS::S32 i_atom) const {
        return mask_index[i_atom+1];
    }
    PS::S32 mask(const PS::S32 i) const {
        return mask_list[i];
    }
    
    //--- broadcast the list data
    void BCastListData(const PS::S32 src = 0){
        PS::Comm::broadcast(mask_list,  list_len,  src);
        PS::Comm::broadcast(mask_index, index_len, src);
    }
    
    //--- check mask
    bool checkMask(const PS::S32 i_atom, const PS::S32 j_atom) const {
        //--- same molecule: search in mask list
        for(PS::S32 j=mask_st(i_atom); j<mask_ed(i_atom); j++){
            //------ found in mask list
            if(mask_list[j] == j_atom){
                return true;
                break;
            }
        }
        //------ not found in mask list
        return false;
    }
    
private:
    //--- error check
    void failure_make_list(){
        std::cout << "Rank:" << PS::Comm::getRank()
                  << "  failure to make mask_list" << std::endl;
        PS::Abort();
    }
};

//--- grobal pointer for intraMask class object
intraMask mask_polymer;
intraMask mask_wall;


//--- intraMask flag check function
//------ return true  : the pair is intramolecular.
//------ return false : the pair is intermolecular.
template<class Tepj>
bool checkIntraMask(const MOL_TYPE i_type,
                    const PS::S32 i_mol,
                    const PS::S32 i_atom,
                    const Tepj &epj,
                    const PS::S32 j) {
    
    //--- check same molecule or not
    if(i_type != epj[j].getType() ) return false;
    if(i_mol != epj[j].getMolID() ) return false;
    
    //--- check intramolecular force pair or not
    //------ in the case of solvent type
    if(i_type == solvent) return true;
        
    //------ in the case of polymer type
    if(i_type == polymer) {
        if( mask_polymer.checkMask(i_atom, epj[j].getAtomID()) ) return true;
        return false;
    }
        
    //------ in the case of wall type
    if(i_type == wall) {
        if( mask_wall.checkMask(i_atom, epj[j].getAtomID()) ) return true;
        return false;
    }
    
    
#ifdef MASK_TEST
    //--- error: definition of intramolecular pair was not found!
    std::cout << "ERROR: definition of intramolecular pair was not found!" << std::endl;
    std::cout << "       check in the `checkIntraMask` function "
              <<        "and `intraMask` data class." << std::endl;
    std::cout << "       i_type : " << i_type << std::endl;
    std::cout << "       i_mol  : " << i_mol  << std::endl;
    std::cout << "       i_atom : " << i_atom << std::endl;
    std::cout << "       j_type : " << epj[j].getType()   << std::endl;
    std::cout << "       j_mol  : " << epj[j].getMolID()  << std::endl;
    std::cout << "       j_atom : " << epj[j].getAtomID() << std::endl;
    
    //--- abort program.
    PS::Abort();
#else
    //--- in release version, non-defined pair error is ignored.
    return false; 
#endif

}


//--- initialize the mask list
void intraMaskList_initialize(){
    PS::S32 mask_len_polymer;
    PS::S32 mask_len_wall;
    PS::S32 n_atom_polymer;
    PS::S32 n_atom_wall;
    
    //--- define the size of mask list object
    if(PS::Comm::getRank() == 0){
        mask_len_polymer = 10;
        mask_len_wall    = 10;
        n_atom_polymer = 1;
        n_atom_wall    = 1;
    }
    PS::Comm::broadcast(&mask_len_polymer, 1, 0);
    PS::Comm::broadcast(&mask_len_wall, 1, 0);
    PS::Comm::broadcast(&n_atom_polymer, 1, 0);
    PS::Comm::broadcast(&n_atom_wall, 1, 0);
    
    //--- initialize the mask list object
    mask_polymer.initialize(mask_len_polymer, n_atom_polymer);
    mask_wall.initialize(mask_len_wall, n_atom_wall);
    
    //--- define mask list data
    if(PS::Comm::getRank() == 0){
        // (input mask list data)
    }
    
    //------ broadcast mask list data
    mask_polymer.BCastListData();
    mask_wall.BCastListData();
}

