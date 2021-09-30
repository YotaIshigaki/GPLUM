/**************************************************************************************************/
/**
* @file  fixed_vector.hpp
* @brief std::vector<T> like wrapper for std::array<T>.
*/
/**************************************************************************************************/
#pragma once

#include <array>
#include <sstream>
#include <initializer_list>


namespace MD_EXT {

    /**
    * @brief std::vector<T> like wrapper for std::array<T>.
    * @tparam <T>    element type.
    * @tparam <Size> max array size.
    */
    template <class T, size_t Size>
    class fixed_vector {
    protected:
        size_t              _n = 0;
        std::array<T, Size> _array;

    private:
        void _M_range_check(const size_t i) const {
            if(i >= this->_n){
                std::ostringstream oss;
                oss << "_M_range_check::out_of_range. ";
                oss << "i = " << i << ", size() = " << this->_n << "\n";
                throw std::out_of_range(oss.str());
            }
        }
        void _M_capacity_check(const size_t sz) const {
            if(sz > Size){
                std::ostringstream oss;
                oss << "_M_capacity_check::length_error. ";
                oss << "sz = " << sz << ", max_size() = " << Size << "\n";
                throw std::length_error(oss.str());
            }
        }

    public:
        //--- internal types
        using reference              = typename decltype(_array)::reference;
        using const_reference        = typename decltype(_array)::const_reference;
        using iterator               = typename decltype(_array)::iterator;
        using const_iterator         = typename decltype(_array)::const_iterator;
        using size_type              = typename decltype(_array)::size_type;
        using difference_type        = typename decltype(_array)::difference_type;
        using value_type             = typename decltype(_array)::value_type;
        using pointer                = typename decltype(_array)::pointer;
        using const_pointer          = typename decltype(_array)::const_pointer;
        using reverse_iterator       = typename decltype(_array)::reverse_iterator;
        using const_reverse_iterator = typename decltype(_array)::const_reverse_iterator;

        fixed_vector()                    = default;
        fixed_vector(const fixed_vector&) = default;
        fixed_vector& operator = (const fixed_vector& rhs) noexcept {
            this->_n     = rhs._n;
            this->_array = rhs._array;
            return *this;
        }
        //fixed_vector(const fixed_vector&&)        = delete;
        //fixed_vector& operator = (fixed_vector&&) = delete;

        ~fixed_vector() = default;

        //--- iterators
        iterator         begin()  noexcept { return this->_array.begin(); }
        iterator         end()    noexcept { return this->_array.begin() + this->_n; }
        reverse_iterator rbegin() noexcept { return std::reverse_iterator<iterator>(this->end()); }
        reverse_iterator rend()   noexcept { return std::reverse_iterator<iterator>(this->begin()); }

        const_iterator   begin()   const noexcept { return this->_array.cbegin(); }
        const_iterator   end()     const noexcept { return this->_array.cbegin() + this->_n; }
        const_iterator   cbegin()  const noexcept { return this->_array.cbegin(); }
        const_iterator   cend()    const noexcept { return this->_array.cbegin() + this->_n; }
        const_reverse_iterator crbegin() const noexcept { return std::reverse_iterator<const_iterator>(this->cend()); }
        const_reverse_iterator crend()   const noexcept { return std::reverse_iterator<const_iterator>(this->cbegin()); }

        //--- area information
        size_type max_size() const { return this->_array.size(); }
        size_type capacity() const { return this->max_size(); }
        size_type size()     const { return this->_n; }
        bool      empty()    const { return (this->_n == 0); }

        //--- value access
        reference operator [] (const size_type i) { return this->_array[i]; }
        reference at(const size_type i) {
            this->_M_range_check(i);
            return this->_array[i];
        }
        reference   front() { return *(this->begin()); }
        reference   back()  { return *(this->end() - 1); }
        value_type* data()  noexcept { return &(this->front()); }

        const_reference operator [] (const size_type i) const { return this->_array[i]; }
        const_reference at(const size_type i) const {
            this->_M_range_check(i);
            return this->_array[i];
        }
        const_reference   front() const { return *(this->cbegin()); }
        const_reference   back()  const { return *(this->cend() - 1); }
        const value_type* data()  const noexcept { return &(this->front()); }

        //--- edit
        void clear(){ this->_n = 0; }

        void swap(fixed_vector& other) noexcept {
            std::swap(this->_n, other._n);
            for(size_t i=0; i<Size; ++i){
                std::swap(this->_array[i], other._array[i]);
            }
        }
        void push_back(const value_type &u){
            this->_M_capacity_check(this->size()+1);
            this->_array[this->_n] = u;
            ++this->_n;
        }
        void pop_back(){
            if( !(this->empty()) ){
                --this->_n;
            }
        }

        iterator insert(iterator pos, const T& x){
            this->_M_capacity_check(this->size()+1);

            size_t index = std::distance(this->begin(), pos);
            for(auto i = this->_n; i>=index; --i){
                this->_array[i+1] = this->_array[i];
            }
            this->_array[index] = x;
            ++this->_n;

            return this->begin() + index;
        }
        iterator erase(iterator pos){
            size_t index = std::distance(this->begin(), pos);
            if(index >= this->_n) return this->end();

            for(auto i = index; i<this->_n-1; ++i){
                this->_array[i] = this->_array[i+1];
            }
            --this->_n;

            return this->begin() + index;
        }

        //--- area operation
        void resize(size_type sz, const T &c){
            if(sz < this->size()){
                const size_type n = this->size() - sz;
                for(size_type i=0; i<n; ++i){
                    this->pop_back();
                }
            }
            if(sz > this->size()){
                this->_M_capacity_check(sz);
                const size_type n = sz - this->size();
                for(size_type i=0; i<n; ++i){
                    this->push_back(c);
                }
            }
        }
        void resize(size_type sz){
            this->resize(sz, T{});
        }

        //--- input for initializer_list
        fixed_vector& operator = (std::initializer_list<T> list){
            auto n = std::distance(list.begin(), list.end());
            if(n >= (decltype(n))Size) throw std::length_error("fixed_vector::_capacity_check");

            this->clear();
            for(auto itr = list.begin(); itr != list.end(); ++itr){
                this->push_back(*itr);
            }

            return *this;
        }
        fixed_vector(std::initializer_list<T> list){
            *this = list;
        }

        //--- interface for different size fixed_vector
        template <size_t Size_rhs>
        fixed_vector& operator = (const fixed_vector<T, Size_rhs> &rhs){
            if(this->max_size() < rhs.size()) throw std::length_error("fixed_vector::_capacity_check");

            this->_n = rhs.size();
            for(size_t i=0; i<rhs.size(); ++i){
                this->_array[i] = rhs[i];
            }
            return *this;
        }
        template <size_t Size_rhs>
        fixed_vector(const fixed_vector<T,Size_rhs> &rhs){
            *this = rhs;
        }
        template <size_t Size_rhs>
        void swap(fixed_vector<T,Size_rhs> &other){
            if(this->_n     > other.max_size()) throw std::length_error("fixed_vector::_capacity_check");
            if(other.size() > this->max_size()) throw std::length_error("fixed_vector::_capacity_check");

            fixed_vector tmp = *this;
            *this = other;
            other = tmp;
        }
    };


    //--- non-member function for fixed_vector<T,Size>
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator == (const fixed_vector<T,Size_lhs> &x, const fixed_vector<T,Size_rhs> &y){
        if(x.size() != y.size()) return false;
        return std::equal(x.begin(), x.end(), y.begin());
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator != (const fixed_vector<T,Size_lhs> &x, const fixed_vector<T,Size_rhs> &y){
        return !(x == y);
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator < (const fixed_vector<T,Size_lhs> &x, const fixed_vector<T,Size_rhs> &y){
        return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator > (const fixed_vector<T,Size_lhs> &x, const fixed_vector<T,Size_rhs> &y){
        return y < x;
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator <= (const fixed_vector<T,Size_lhs> &x, const fixed_vector<T,Size_rhs> &y){
        return !(x > y);
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator >= (const fixed_vector<T,Size_lhs> &x, const fixed_vector<T,Size_rhs> &y){
        return !(x < y);
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    void swap(fixed_vector<T,Size_lhs> &x, fixed_vector<T,Size_rhs> &y){
        x.swap(y);
    }

}
