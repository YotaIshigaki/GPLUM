//***************************************************************************************
//  This program is the temporary buffer data class for intramolecular force calculation.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once


//--- parameter table for intermolecular interactions -----------------------------------
//------ this data is used for initialize
//------ parameter class for table
class InterForceParam {
  public:
    PS::F64 mass;
    PS::F64 charge;
    PS::F64 VDW_R;
    PS::F64 VDW_D;
};

//------ grobal object of parameter table
namespace Grobal {
    MD_EXT::MultiTable<InterForceParam, ATOM_TYPE> coefTable_inter(static_cast<int>(ATOM_TYPE::ENUM_N_TOT));
}

//--- temporary data buffer class for intramolecular interaction ------------------------
template<std::size_t intra_list_len>
class IntraForceTarget{
  private:
    PS::S32 n;
    std::array<PS::S32,    intra_list_len> id_list;
    std::array<PS::S32,    intra_list_len> index_list;
    std::array<ATOM_TYPE,  intra_list_len> type_list;
    std::array<PS::F64vec, intra_list_len> pos_list;
    
    std::array<PS::F64vec, intra_list_len-1> rij_list;
    std::array<PS::F64,    intra_list_len-1> r_list;
    std::array<PS::F64,    intra_list_len-1> r_inv_list;
    
  public:
    template<class Tepj>
    inline void add(const Tepj* epj,
                    const PS::S32 & j){
        //--- input data
        this->id_list.at(n)    = epj[j].getAtomID();
        this->index_list.at(n) = j;
        this->type_list.at(n)  = epj[j].getAtomType();
        this->pos_list.at(n)   = Normalize::realPos( epj[j].getPos() );
        //--- calculate r and r_inv
        if( n > 0 ){
            this->rij_list.at(n-1)   = Normalize::periodicAdjust(this->pos_list.at(n) - this->pos_list.at(n-1));
            this->r_list.at(n-1)     = sqrt(this->rij_list.at(n-1)*this->rij_list.at(n-1));
            this->r_inv_list.at(n-1) = 1.0/this->r_list.at(n-1);
        }
        this->n++;
    }
    PS::S32 getN() const {
        return this->n;
    }
    inline PS::S32 getID(const PS::S32 & i) const {
        return this->id_list.at(i);
    }
    inline PS::S32 getIndex(const PS::S32 & i) const {
        return this->index_list.at(i);
    }
    inline ATOM_TYPE getAtomType(const PS::S32 & i) const {
        return this->type_list.at(i);
    }
    inline PS::F64vec getRij(const PS::S32 & i) const {
        return this->rij_list.at(i);
    }
    inline PS::F64 getR(const PS::S32 & i) const {
        return this->r_list.at(i);
    }
    inline PS::F64 getRinv(const PS::S32 & i) const {
        return this->r_inv_list.at(i);
    }
    inline void pop(){
        #ifdef ARRAY_CB
            if( n <= 0 ){
                throw std::range_error("invalid action to list. list is empty.");
            }
        #endif
        this->n--;
        if(n>0){
            rij_list.at(n-1)   = PS::F64vec(0.0, 0.0, 0.0);
            r_list.at(n-1)     = 0.0;
            r_inv_list.at(n-1) = 0.0;
        }
    }
    inline bool checkParentNode(const PS::S32 & atom_id){
        return (this->id_list.at(n-1) == atom_id);
    }
    inline bool checkRootNode(const PS::S32 & atom_id){
        return (this->id_list.at(0) == atom_id);
    }
    void clear(){
        this->n = 0;
        this->id_list.fill(-1);
        this->index_list.fill(-1);
        this->r_list.fill(0.0);
    }
};
