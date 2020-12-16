/**************************************************************************************************/
/**
* @file  md_ext_basic_connect.hpp
* @brief connection information data class for PS::ParticleSystem<FP>.
*/
/**************************************************************************************************/
#pragma once

#include "fixed_vector.hpp"

namespace MD_EXT {

    /**
    * @brief connection information for PS::ParticleSystem<FP>.
    * @tparam <T>    ID data type.
    * @tparam <Size> max array size.
    */
    template <class T, size_t Size>
    class basic_connect {
    protected:
        fixed_vector<T, Size> _data;

    public:
        //--- internal types
        using reference              = typename fixed_vector<T, Size>::reference;
        using const_reference        = typename fixed_vector<T, Size>::const_reference;
        using iterator               = typename fixed_vector<T, Size>::iterator;
        using const_iterator         = typename fixed_vector<T, Size>::const_iterator;
        using size_type              = typename fixed_vector<T, Size>::size_type;
        using difference_type        = typename fixed_vector<T, Size>::difference_type;
        using value_type             = typename fixed_vector<T, Size>::value_type;
        using pointer                = typename fixed_vector<T, Size>::pointer;
        using const_pointer          = typename fixed_vector<T, Size>::const_pointer;
        using reverse_iterator       = typename fixed_vector<T, Size>::reverse_iterator;
        using const_reverse_iterator = typename fixed_vector<T, Size>::const_reverse_iterator;

        basic_connect(){ this->_data.clear(); }
        basic_connect(const basic_connect& rhs){
            this->_data = rhs._data;
        }
        basic_connect& operator = (const basic_connect& rhs){
            this->_data = rhs._data;
            return *this;
        }
    //    basic_connect(const basic_connect&&)              = delete;
    //    basic_connect& operator = (basic_connect&&)       = delete;

        ~basic_connect() = default;

        //--- iterators
        const_iterator         begin()   const noexcept { return this->_data.begin();   }
        const_iterator         end()     const noexcept { return this->_data.end();     }
        const_iterator         cbegin()  const noexcept { return this->_data.cbegin();  }
        const_iterator         cend()    const noexcept { return this->_data.cend();    }
        const_reverse_iterator crbegin() const noexcept { return this->_data.crbegin(); }
        const_reverse_iterator crend()   const noexcept { return this->_data.crebd();   }

        //--- area information
        size_type max_size() const { return this->_data.max_size(); }
        size_type capacity() const { return this->_data.capacity(); }
        size_type size()     const { return this->_data.size();     }
        bool      empty()    const { return this->_data.empty();    }

        //--- value access
        const_reference operator [] (const size_type n) const { return this->_data[n];    }
        const_reference           at(const size_type n) const { return this->_data.at(n); }

        const_reference   front() const          { return this->_data.front(); }
        const_reference   back()  const          { return this->_data.back();  }
        const value_type* data()  const noexcept { return this->_data.data();  }

    private:
        iterator _find_impl(const value_type &u) const {
            return (iterator)std::find(this->_data.begin(), this->_data.end(), u);
        }
    public:
        const_iterator find(const value_type &u) const {
            return (const_iterator)this->_find_impl(u);
        }
        int count(const value_type &u) const {
            if(this->find(u) == this->end()){
                return 0;
            } else {
                return 1;
            }
        }

        //--- edit
        void clear(){ this->_data.clear(); }

        //! @detail this class has unique value only. Do nothing when the object has same value with the arg "u".
        void add(const value_type &u){
            auto itr = this->_find_impl(u);
            if(itr != this->_data.end()) return;

            this->_data.push_back(u);
        }
        //! @detail Do nothing when the arg value "u" was not found in this object.
        void remove(const value_type &u){
            auto itr = this->_find_impl(u);
            if(itr == this->_data.end()) return;

            this->_data.erase(itr);
        }

        //--- combine 2 basic_connect objects
        template <size_t Size_rhs>
        basic_connect<T, Size>& operator += (const basic_connect<T, Size_rhs> &rhs){
            for(const auto& e : rhs){
                this->add(e);
            }
            return *this;
        }
        template <size_t Size_rhs>
        basic_connect<T, Size> operator + (const basic_connect<T, Size_rhs> &rhs) const {
            basic_connect<T, Size> result = *this;
            result += rhs;
            return result;
        }

        //--- interface for different size connect_base<T,Size>
        template <size_t Size_rhs>
        basic_connect(const basic_connect<T,Size_rhs>& rhs){
            this->clear();
            for(const auto& e : rhs){
                this->_data.push_back(e);
            }
        }
        template <size_t Size_rhs>
        basic_connect<T,Size>& operator = (const basic_connect<T, Size_rhs> rhs){
            this->clear();
            for(const auto& e : rhs){
                this->_data.push_back(e);
            }
            return *this;
        }
    };


    //--- non-member function for basic_connect<T,Size>
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator == (const basic_connect<T,Size_lhs> &x, const basic_connect<T,Size_rhs> &y){
        return x._data == y._data;
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator != (const basic_connect<T,Size_lhs> &x, const basic_connect<T,Size_rhs> &y){
        return !(x == y);
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator < (const basic_connect<T,Size_lhs> &x, const basic_connect<T,Size_rhs> &y){
        return x._data < y._data;
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator > (const basic_connect<T,Size_lhs> &x, const basic_connect<T,Size_rhs> &y){
        return y < x;
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator <= (const basic_connect<T,Size_lhs> &x, const basic_connect<T,Size_rhs> &y){
        return !(x > y);
    }
    template <class T, size_t Size_lhs, size_t Size_rhs>
    bool operator >= (const basic_connect<T,Size_lhs> &x, const basic_connect<T,Size_rhs> &y){
        return !(x < y);
    }

}
