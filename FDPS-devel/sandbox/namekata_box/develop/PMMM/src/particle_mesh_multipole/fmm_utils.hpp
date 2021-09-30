#pragma once
#include <cassert>
#include <unordered_map>

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        static inline double factorial(const int i){
            assert(i >= 0);
            return i ? double(i) * factorial(i-1) 
                     : 1.0;
        }

        static inline double factinv(const int i){
            return 1.0 / factorial(i);
        }

        template <typename real_t=double>
        class TableForRlm {
        public:
            std::vector<real_t> tbl_inv;
            std::vector<real_t> tbl_factinv;
        
            TableForRlm() {}
        
            void init(int p) {
                for(int i=0; i<p+1; i++){
                    tbl_inv.push_back(real_t(1.0) / real_t(i));
                }
                for(int i=0; i<2*p+1; i++){
                    assert(factorial(i) > 0);
                    tbl_factinv.push_back(1.0 / factorial(i));
                }
            }
        };

        template <typename real_t=double>
        class TableForSlm {
        public:
            std::vector<real_t> tbl_inv;

            TableForSlm() {}

            void init(int p) {
                for(int i=0; i<p+1; i++){
                    tbl_inv.push_back(real_t(1.0) / real_t(i));
                }
            }
        };

        template <typename Ttbl>
        class TableManager {
        public:
            std::unordered_map<int, int> tbl_map;
            std::vector<Ttbl> tbl_vec;
        
            TableManager() {}
        
            void set(int p) {
                const int adr = tbl_vec.size();
                auto itr = tbl_map.find(p);
                if ( itr != tbl_map.end()) {
                    // do nothing
                } else {
                    Ttbl tmp;
                    tmp.init(p);
                    tbl_vec.push_back(tmp);
                    tbl_map[p] = adr;
                }
            }
        
            Ttbl & get(int p){
                const int adr = tbl_map[p];
                return tbl_vec[adr];
            }
        };

    } // END of namespace of ParticleSimulator
} // END of namespace of ParticleSimulator
