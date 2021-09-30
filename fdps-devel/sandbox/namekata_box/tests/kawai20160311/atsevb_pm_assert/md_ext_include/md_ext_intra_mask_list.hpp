//***************************************************************************************
//  This program is the intermolecular mask list for "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include<unordered_map>

namespace MD_EXT {
    
    //--- connection class
    template<class Tid ,PS::S32 size>
    class Connection {
      private:
        std::array<Tid, size> pair;
        PS::S32 n_pair;
        
      public:
        //--- clear
        void clear(){
            n_pair = 0;
            pair.fill(-1);  // set illigal value
        }
        
        //--- constructor
        Connection(){
            this->clear();
        }
      
        //--- access to pair information
        void add(const Tid & atomID){
            pair.at(n_pair) = atomID;
            n_pair++;
        }
        void set(const PS::S32 & i, const Tid & atomID){
          #ifdef ARRAY_CB
            if( !(i < n_pair) ){
                throw std::out_of_range("invalid index. index must be < getN().");
            }
          #endif
            pair.at(i) = atomID;
        }
        inline PS::S32 getN() const {
            return this->n_pair;
        }
        inline Tid get(const PS::S32 & i) const {
          #ifdef ARRAY_CB
            return this->pair.at(i);    // with array boundary chack
          #else
            return this->pair[i];       // without array boundary check
          #endif
        }
    };
    
    //--- mask list for PS::TreeForForce<>::Scatter
    template<class Tid>
    class IntraMaskList {
        private:
          std::vector<Tid> intraMaskList;
          
        public:
          void addIntraList(const Tid & j_atom_id){
              this->intraMaskList.push_back(j_atom_id);
          }
          inline bool checkIntraList(const Tid & j_atom_id) const {
              for(auto & intraPair : this->intraMaskList){
                  if(j_atom_id == intraPair) return true;
              }
              return false;
          }
          void clearIntraList(){
              this->intraMaskList.clear();
          }
          
          //--- data copy
          std::vector<Tid> getIntraMaskList() const {
              return this->intraMaskList;
          }
          template<class Tptcl>
          void copyIntraList(const Tptcl & fp){
              this->intraMaskList = fp.getIntraMaskList();
          }
    };
    
    //--- std::array version intra mask list. fixed list length
    template<class Tid, std::size_t len>
    class IntraMaskList_fix {
      private:
        std::array<Tid, len> intraMaskList;
        PS::S32 n_pair;
        
      public:
          void clearIntraList(){
              this->intraMaskList.fill(-1);
          }
          
          //--- constructor
          IntraMaskList_fix(){
              this->clearIntraList();
              this->n_pair = 0;
          }
          
          void addIntraList(const Tid & j_atom_id){
              this->intraMaskList.at(n_pair) = j_atom_id;
              this->n_pair++;
          }
          
          inline bool checkIntraList(const Tid & j_atom_id) const {
              for(PS::S32 i=0; i<n_pair; i++){
                  if(j_atom_id == intraMaskList[i]) return true;
              }
              return false;
          }
          
          //--- data copy
          std::array<Tid, len> getIntraMaskList() const {
              return this->intraMaskList;
          }
          template<class Tptcl>
          void copyIntraList(const Tptcl & fp){
              this->intraMaskList = fp.getIntraMaskList();
          }
    };

    //--- interface function
    template<class Troot, class Tpsys, class Tid, class Tmap>
    void inputIntraList(Troot *const root_node,   // root node (const pointer)
                        const Tpsys & system,     // FullParticle system
                        const Tid parent_id,      // parent node for tree search
                        const Tid node_id,        // current node for tree search
                        const Tmap & IDtable,     // map table for atom_id and i of system[i]
                        const PS::S32 order = 3){ // bond=1, angle=2, tortion=3
        
        //--- convert atom_id of node to "i"
        PS::S32 i = IDtable.at(node_id);
        
        //--- tmp array for search in current node
        std::vector<Tid> v_tmp;
        system[i].getAllConnect(v_tmp);
        
        //--- add pair into intraList
        for(Tid & j_atom_id : v_tmp){
            if(parent_id == j_atom_id) continue;                 // skip parent node
            if(root_node->getAtomID() == j_atom_id) continue;    // check closed loop
            if(root_node->checkIntraList(j_atom_id) ) continue;  // check double count
            
            root_node->addIntraList(j_atom_id); // add pair into intraList
            
            //--- recursive search in tree data
            if(order > 1) inputIntraList(root_node, system, node_id, j_atom_id, IDtable, order-1);
        }
    }
    
    //--- initialize IntraList class in FP class
    template<class Tid, class Tpsys>
    void Init_IntraList(Tpsys & system, PS::S32 order = 3){
        PS::S32 n_loc = system.getNumberOfParticleLocal();
        
        //--- temporary data
        std::unordered_map<Tid, PS::S32> IDtable;
        
        //--- create temporary table for atomID and "i"
        for(PS::S32 i=0; i<n_loc; i++){
            IDtable[system[i].getAtomID()] = i;
        }
        
        //--- define intraList data
        for(PS::S32 i=0; i<n_loc; i++){
            system[i].clearIntraList();  // must be clear
            system[i].addIntraList( system[i].getAtomID() ); // input itself into intraList
            inputIntraList(&system[i],             // root node (pointer)
                           system,                 // FullParticle system
                           system[i].getAtomID(),  // use system[i].getAtomID()
                           system[i].getAtomID(),  // use system[i].getAtomID()
                           IDtable,                // map table for atom_id and i of system[i]
                           order);                 // order: bond=1, angle=2, tortion=3
        }
    }
}
