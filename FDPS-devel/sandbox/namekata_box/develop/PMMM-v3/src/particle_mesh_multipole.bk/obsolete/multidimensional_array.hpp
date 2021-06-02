#pragma once
#include <typeinfo>
#include <cassert>
#include <iostream>
#include "particle_mesh_multipole_defs.hpp"

#define SANITY_CHECK_MULTIDIMENSIONAL_ARRAY (1)

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        template <typename T, int DIM>
        class MultidimensionalArray {
        private:
            MultidimensionalArray(const MultidimensionalArray &);
            MultidimensionalArray & operator = (const MultidimensionalArray &);
            T * data_;
            int sizes_[DIM];
            int capacity_;
            int alloc_mode_;
    
        public:
            // Constructors
            MultidimensionalArray() : data_(NULL), capacity_(0) {}
            // Destructor
            ~MultidimensionalArray() {
                if (data_ != NULL) {
                    if (alloc_mode_ == 0) {
                        delete [] data_;
                    } else if (alloc_mode_ == 1) {
                        fftw_free(data_);
                    }
                    data_ = NULL;
                }
            }
    
            void initialize (const int sizes[DIM], const int alloc_mode=0) {
                assert(alloc_mode == 0 || alloc_mode == 1);
                alloc_mode_ = alloc_mode;
                capacity_ = 1;
                for (int i = 0; i < DIM; i++) {
                    sizes_[i] = sizes[i];
                    capacity_ *= sizes[i];
                }
                assert(capacity_ > 0);
                if (alloc_mode_ == 0) {
                    data_ = new T[capacity_];
                } else if (alloc_mode_ == 1) {
                    if (typeid(T) == typeid(fftw_real_t) || typeid(T) == typeid(fftw_cplx_t)) {
                        data_ = (T *) fftw_malloc((size_t) sizeof(T) * capacity_ );
                    } else {
                        std::string errmsg = "You cannot use alloc_mode_ = 1 for data types\n";
                        errmsg += " whose size is not consistent with neither double nor fftw_complex.";
                        PARTICLE_SIMULATOR_PRINT_ERROR(errmsg);
                        Abort(-1);
                    }
                }
            }

            void freeMem() {
                for (int i=0; i<DIM; i++) sizes_[i]=0;
                capacity_ = 0;
                if (data_ != NULL) {
                    if (alloc_mode_ == 0) delete [] data_;
                    else if (alloc_mode_ == 1) fftw_free(data_);
                    data_ = NULL;
                }
            }
    
            int capacity() const { return capacity_; }
    
            const T & operator [] (const int i) const {
#if SANITY_CHECK_MULTIDIMENSIONAL_ARRAY > 0
                assert(i >= 0);
                assert(i < capacity_);
#endif
                return data_[i]; 
            }
    
            T & operator [] (const int i) {
#if SANITY_CHECK_MULTIDIMENSIONAL_ARRAY > 0
                assert(i >= 0);
                assert(i < capacity_);
#endif
                return data_[i]; 
            }
    
            T * data() { return data_; }
    
            const T * data() const { return data_; }
    
            int getDimension() const { return DIM; }
    
            int getSize(const int i) const {
#if SANITY_CHECK_MULTIDIMENSIONAL_ARRAY > 0
                assert(i < DIM);
#endif
                return sizes_[i];
            }
    
            T * getPointer(const int i=0) const { return data_ + i; }
    
            inline const int getArrayIndex(const int idx_md_repr[DIM]) const {
                // md := Multi-Dimension, repr := REPResentation
#if SANITY_CHECK_MULTIDIMENSIONAL_ARRAY > 0
                for (int i = 0; i < DIM; i++) {
                    assert(idx_md_repr[i] >= 0);
                    assert(idx_md_repr[i] < sizes_[i]);
                }
#endif
                int idx = 0, stride = 1;
                for (int i = DIM-1; i >= 0; i--) {
                    idx += stride * idx_md_repr[i];
                    stride *= sizes_[i];
                }
#if SANITY_CHECK_MULTIDIMENSIONAL_ARRAY > 0
                assert(idx < capacity_);
#endif
                return idx;
            }
    
            inline const int getArrayIndex(const int i0,
                                           const int i1,
                                           const int i2) const {
                const int idx = sizes_[2]*sizes_[1]*i0
                              + sizes_[2]*i1
                              + i2;
#if SANITY_CHECK_MULTIDIMENSIONAL_ARRAY > 0
                assert(idx < capacity_);
#endif
                return idx;
            }
    
            inline const int getArrayIndex(const int i0,
                                           const int i1,
                                           const int i2,
                                           const int i3) const {
                const int idx = sizes_[3]*sizes_[2]*sizes_[1]*i0
                              + sizes_[3]*sizes_[2]*i1
                              + sizes_[3]*i2
                              + i3;
#if SANITY_CHECK_MULTIDIMENSIONAL_ARRAY > 0
                assert(idx < capacity_);
#endif
                return idx;
            }
    
            const T & operator () (const int idx_md_repr[DIM]) const {
                const int idx = getArrayIndex(idx_md_repr);
                return data_[idx]; 
            }
    
            const T & operator () (const int i0, 
                                   const int i1,
                                   const int i2) const {
                 const int idx = getArrayIndex(i0,i1,i2);
                 return data_[idx];
            }
    
            const T & operator () (const int i0, 
                                   const int i1,
                                   const int i2,
                                   const int i3) const {
                 const int idx = getArrayIndex(i0,i1,i2,i3);
                 return data_[idx];
            }
    
            T & operator () (const int idx_md_repr[DIM]) {
                const int idx = getArrayIndex(idx_md_repr);
                return data_[idx]; 
            }
    
            T & operator () (const int i0, 
                             const int i1,
                             const int i2) {
                 const int idx = getArrayIndex(i0,i1,i2);
                 return data_[idx];
            }
    
            T & operator () (const int i0, 
                             const int i1,
                             const int i2,
                             const int i3) {
                 const int idx = getArrayIndex(i0,i1,i2,i3);
                 return data_[idx];
            }
    
            size_t getMemSize() const { return capacity_ * sizeof(T); }
    
        };

    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
