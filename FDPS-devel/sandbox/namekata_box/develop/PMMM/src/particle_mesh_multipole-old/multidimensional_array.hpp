#pragma once
#include <typeinfo>
#include <cassert>
#include <iostream>

#define SANITY_CHECK_MULTIDIMENSIONAL_ARRAY (1)

namespace ParticleSimulator {

    template <typename T, int DIM>
    class MultidimensionalArray {
    private:
        MultidimensionalArray(const MultidimensionalArray &);
        MultidimensionalArray & operator = (const MultidimensionalArray &);
        T * data_;
        int sizes_[DIM];
        int capacity_;

    public:
        // Constructors
        MultidimensionalArray() : data_(NULL), capacity_(0) {}
        // Destructor
        ~MultidimensionalArray() {
            if (data_ != NULL) {
                delete [] data_;
                data_ = NULL;
            }
        }

        void initialize (const int sizes[DIM]) {
            capacity_ = 1;
            for (int i = 0; i < DIM; i++) {
                sizes_[i] = sizes[i];
                capacity_ *= sizes[i];
            }
            data_ = new T[capacity_];
        }

        void initialize (const int size0,
                         const int size1,
                         const int size2) {
            sizes_[0] = size0;
            sizes_[1] = size1;
            sizes_[2] = size2;
            capacity_ = size0 * size1 * size2;
            data_ = new T[capacity_];
        }

        void initialize (const int size0,
                         const int size1,
                         const int size2,
                         const int size3) {
            sizes_[0] = size0;
            sizes_[1] = size1;
            sizes_[2] = size2;
            sizes_[3] = size3;
            capacity_ = size0 * size1 * size2 * size3;
            data_ = new T[capacity_];
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

} // END of namespace of ParticleSimulator
