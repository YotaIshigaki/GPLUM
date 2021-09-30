//***************************************************************************************
//  This program is the intramolecular force calculation for "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include<unordered_map>
#include<cstdint>

namespace MD_EXT {
    
  //  //--- std::Matrix version. only use in the case of small size of "Ttype" ------------
  //  // not defined yet
  //  
  //  template<class Tvalue>
  //  class CoefInMatrix_ : 
  //    public Tvalue {
  //    public:
  //      bool enable;
  //  };
  //  
  //  //------ this class consume much large memory cost. be carefull.
  //  template<class Tvalue, std::size_t ... dim_size>
  //  class CoefTable_Matrix_base_ {
  //    private:
  //      std::Mattrix<CoefInMatrix_<Tvalue>, sizeof...(Ttype)> table(dim_size...);
  //    public:
  //      void setCoef(const Tvalue &value ,const std::size_t & ... dim_size){
  //          this->table[dim_size...] = value;
  //      }
  //      inline void getCoef(const Tvalue &value ,const std::size_t & ... dim_size){
  //          value.copyCoef(this->table[dim_size...]);
  //      }
  //  };
    
    //--- std::unordered_map version ----------------------------------------------------
    template<class Tkey, std::size_t n_key_max, int_fast32_t window,
             class Tvalue, typename...Ttype>
    class MultiTable_base_ {
      private:
        //------ check number of type.
        static_assert(sizeof...(Ttype) <= n_key_max,
                      "the number of key types is too large.");
        static_assert( 0 < window,
                      "window size must be > 0.");
        std::unordered_map<Tkey, Tvalue> table;
        
        //--- convert key
        //------ terminator
        inline void typecast_(const Tkey & key) const { }
        //------ internal function for key generation
        template<class Head_, typename ... Args_>
        inline void typecast_(Tkey & key, const Head_ & first, const Args_ & ... args) const {
            #ifdef MULTI_TABLE_WINDOW_TEST
                if( !(0 <= static_cast<Tkey>(first) && static_cast<Tkey>(first) < window) ){
                    std::cout << "ERORR: key value must be in range of 0~" << (window-1) << std::endl
                              << "  key: " << static_cast<Tkey>(first) << std::endl;
                    throw std::range_error("ERROR: invalid key value.");
                }
            #endif
            key = key*window + static_cast<Tkey>(first);
            typecast_(key, args...);
        }
        
      public:
        //--- reserve table size
        void reserve(std::size_t size){
            this->table.reserve(size);
        }
        
        //--- constructor
        MultiTable_base_(){
            this->table.clear();
            this->table.max_load_factor(0.7);
        }
        MultiTable_base_(std::size_t size){
            this->table.clear();
            this->table.max_load_factor(0.7);
            this->reserve(size);
        }
        
        //--- input to table
        void set(const Tvalue& value, const Ttype & ...args){
            //--- generate the key
            Tkey key = 0;
            this->typecast_(key, args...);
            
            //--- input
            this->table[key] = value;
        }
        //--- read from table
        inline Tvalue get(const Ttype & ...args) const {
            //--- generate the key
            Tkey key = 0;
            this->typecast_(key, args...);
            
            //--- output
            return this->table.at(key);
        }
        
        //--- check the value is setted or not. true: exist / false: void.
        inline bool find(const Ttype & ...args){
            Tkey key = 0;
            this->typecast_(key, args...);
            
            return ( table.count(key) != 0 );
        }
        
        //--- clear all data
        void clearTable(){
            this->table.clear();
        }
        
        //--- key generation test
        int64_t test_key_gen_(const Ttype & ... args){
            //--- generate the key
            Tkey key = 0;
            this->typecast_(key, args...);
            
            //--- output
            return static_cast<int64_t>(key);
        }
    };
}