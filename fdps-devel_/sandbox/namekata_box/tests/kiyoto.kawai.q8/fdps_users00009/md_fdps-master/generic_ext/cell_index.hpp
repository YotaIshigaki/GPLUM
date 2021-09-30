/**************************************************************************************************/
/**
* @file  cell_index.hpp
* @brief neighbor search tool by cell index method.
*/
/**************************************************************************************************/
#pragma once

#include <vector>
#include <utility>
#include <tuple>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <sstream>

#include <particle_simulator.hpp>


namespace MD_EXT {

    //! @brief cell index neigbor list manager for single process region.
    //! @tparam <T> ID type.
    template <class T>
    class CellIndex{
    public:
        using value_type = T;
        class index_type {
        public:
            int x = 0;
            int y = 0;
            int z = 0;

            index_type() = default;
            index_type(const index_type&) = default;
            index_type(const int x, const int y, const int z){
                this->x = x;
                this->y = y;
                this->z = z;
            }
            ~index_type() = default;

            index_type& operator = (const index_type &rhs){
                this->x = rhs.x;
                this->y = rhs.y;
                this->z = rhs.z;
                return *this;
            }

            index_type& operator += (const index_type &rhs){
                this->x += rhs.x;
                this->y += rhs.y;
                this->z += rhs.z;
                return *this;
            }
            index_type& operator -= (const index_type &rhs){
                this->x -= rhs.x;
                this->y -= rhs.y;
                this->z -= rhs.z;
                return *this;
            }
            index_type operator + (const index_type &rhs){
                *this += rhs;
                return *this;
            }
            index_type operator - (const index_type &rhs){
                *this -= rhs;
                return *this;
            }
        };

    private:
        std::vector<std::vector<T>> cell;
        int nx = 1;
        int ny = 1;
        int nz = 1;

        PS::F32vec domain_l{0.0, 0.0, 0.0};
        PS::F32vec domain_u{1.0, 1.0, 1.0};

        PS::F32vec cell_size{1.0, 1.0, 1.0};
        PS::F32vec cell_size_inv{1.0, 1.0, 1.0};

        std::vector<index_type> buff_index;
        std::vector<value_type> buff_value;

        //--- internal implementations
        void _checkPosInDomain(const PS::F32vec &pos) const {
            if(pos.x < this->domain_l.x || pos.x >= this->domain_u.x ||
               pos.y < this->domain_l.y || pos.y >= this->domain_u.y ||
               pos.z < this->domain_l.z || pos.z >= this->domain_u.z   ){
                std::ostringstream oss;
                oss << "pos is outside of domain." << "\n"
                    << "  pos.x = " << pos.x << " | domain.x = " << this->domain_l.x << " ~ " << this->domain_u.x << "\n"
                    << "  pos.y = " << pos.y << " | domain.y = " << this->domain_l.y << " ~ " << this->domain_u.y << "\n"
                    << "  pos.z = " << pos.z << " | domain.z = " << this->domain_l.z << " ~ " << this->domain_u.z << "\n";
                throw std::invalid_argument(oss.str());
            }
        }
        void _updateCellSize(){
            this->cell_size.x = (this->domain_u.x - this->domain_l.x)/static_cast<double>(this->nx);
            this->cell_size.y = (this->domain_u.y - this->domain_l.y)/static_cast<double>(this->ny);
            this->cell_size.z = (this->domain_u.z - this->domain_l.z)/static_cast<double>(this->nz);

            this->cell_size_inv.x = 1.0/this->cell_size.x;
            this->cell_size_inv.y = 1.0/this->cell_size.y;
            this->cell_size_inv.z = 1.0/this->cell_size.z;
        }
        index_type _getIndex(const PS::F32vec &pos) const {
            this->_checkPosInDomain(pos);
            const PS::F32vec pos_local = pos - this->domain_l;
            int ix = std::floor(pos_local.x*this->cell_size_inv.x);
            int iy = std::floor(pos_local.y*this->cell_size_inv.y);
            int iz = std::floor(pos_local.z*this->cell_size_inv.z);
            return index_type{ix, iy, iz};
        }
        size_t _access(const index_type &index) const {
            if(index.x < 0 || index.x >= this->nx ||
               index.y < 0 || index.y >= this->ny ||
               index.z < 0 || index.z >= this->nz   ){
                std::ostringstream oss;
                oss << "index is out of range." << "\n"
                    << "  index.x = " << index.x << " | range.x = 0 ~ " << this->nx << "\n"
                    << "  index.y = " << index.y << " | range.y = 0 ~ " << this->ny << "\n"
                    << "  index.z = " << index.z << " | range.z = 0 ~ " << this->nz << "\n";
                throw std::length_error(oss.str());
            }
            return index.x
                 + index.y*(this->nx)
                 + index.z*(this->nx*this->ny);
        }

    public:
        //--- initializer
        void clear(){
            for(auto& c : this->cell){
                c.clear();
            }
        }
        void initDomain(const PS::F32vec &lower, const PS::F32vec &upper){
            if(lower.x >= upper.x ||
               lower.y >= upper.y ||
               lower.z >= upper.z  ){
                std::ostringstream oss;
                oss << "each element of upper must be larger than that of lower." << "\n"
                    << "  lower.x = " << lower.x << " < upper.x = " << upper.x << "\n"
                    << "  lower.y = " << lower.y << " < upper.y = " << upper.y << "\n"
                    << "  lower.z = " << lower.z << " < upper.z = " << upper.z << "\n";
                throw std::invalid_argument(oss.str());
            }

            this->domain_l = lower;
            this->domain_u = upper;

            this->_updateCellSize();
            this->clear();
        }
        void initIndex(const int nx, const int ny, const int nz){
            if(nx <= 0 || ny <= 0 || nz <= 0){
                std::ostringstream oss;
                oss << "index size must be > 0" << "\n"
                    << "  nx = " << nx << "\n"
                    << "  ny = " << ny << "\n"
                    << "  nz = " << nz << "\n";
                throw std::invalid_argument(oss.str());
            }

            this->nx = nx;
            this->ny = ny;
            this->nz = nz;

            int len = nx*ny*nz;
            this->cell.resize(len);

            this->_updateCellSize();
            this->clear();
        }
        void init(const PS::F32vec &lower,
                  const PS::F32vec &upper,
                  const PS::F32     r_search){
            if(r_search <= 0.0){
                std::ostringstream oss;
                oss << "r_search = " << r_search << " must be > 0" << "\n";
                throw std::invalid_argument(oss.str());
            }
            if(r_search*2.0 >= (upper.x - lower.x) ||
               r_search*2.0 >= (upper.y - lower.y) ||
               r_search*2.0 >= (upper.z - lower.z)   ){
                std::ostringstream oss;
                oss << "r_search must be < 0.5*(x, y, or z size of the domain)" << "\n"
                    << "  r_search = " << r_search << "\n"
                    << "  domain_size.x = " << upper.x - lower.x << "\n"
                    << "  domain_size.y = " << upper.y - lower.y << "\n"
                    << "  domain_size.z = " << upper.z - lower.z << "\n";
                throw std::invalid_argument(oss.str());
            }

            int nx = std::max(4, static_cast<int>( (upper.x - lower.x)/r_search ) + 2 );
            int ny = std::max(4, static_cast<int>( (upper.y - lower.y)/r_search ) + 2 );
            int nz = std::max(4, static_cast<int>( (upper.z - lower.z)/r_search ) + 2 );

            this->initDomain(lower, upper);
            this->initIndex(nx, ny, nz);
        }

        //--- internal properties
        std::pair<PS::F32vec,
                  PS::F32vec> get_domain()                     const { return std::make_pair(this->domain_l, this->domain_u); }
        PS::F32vec            get_cell_size()                  const { return this->cell_size;                                }
        index_type            get_grid_size()                  const { return index_type{this->nx, this->ny, this->nz};       }
        index_type            get_index(const PS::F32vec &pos) const { return this->_getIndex(pos);                           }

        //--- cell accessor
        std::vector<T>& operator () (const index_type &i_xyz){
            return this->cell[ this->_access(i_xyz) ];
        }
        std::vector<T>& operator () (const int ix, const int iy, const int iz){
            return this->cell[ this->_access( index_type{ix, iy, iz} ) ];
        }

        const std::vector<T>& operator () (const index_type &i_xyz) const {
            return this->cell[ this->_access(i_xyz) ];
        }
        const std::vector<T>& operator () (const int ix, const int iy, const int iz) const {
            return this->cell[ this->_access( index_type{ix, iy, iz} ) ];
        }

        //--- data interface
        void add(const PS::F32vec pos, const T &data){
            index_type index = this->_getIndex(pos);
            this->cell[ this->_access(index) ].push_back(data);
        }
        void get_index_list(const PS::F32vec              &pos,
                            const PS::F32                  r_search,
                                  std::vector<index_type> &list     ) const {
            list.clear();
            int range_x = static_cast<int>(r_search*this->cell_size_inv.x) + 1;
            int range_y = static_cast<int>(r_search*this->cell_size_inv.y) + 1;
            int range_z = static_cast<int>(r_search*this->cell_size_inv.z) + 1;
            index_type sight = this->_getIndex(pos);

            for(int iz=sight.z-range_z; iz<=sight.z+range_z; ++iz){
                for(int iy=sight.y-range_y; iy<=sight.y+range_y; ++iy){
                    for(int ix=sight.x-range_x; ix<=sight.x+range_x; ++ix){
                        list.push_back( index_type{ (ix + this->nx)%(this->nx),
                                                    (iy + this->ny)%(this->ny),
                                                    (iz + this->nz)%(this->nz) } );
                    }
                }
            }

        }
        const std::vector<index_type>& get_index_list(const PS::F32vec &pos,
                                                      const PS::F32     r_search){
            this->buff_index.clear();
            int range_x = static_cast<int>(r_search*this->cell_size_inv.x) + 1;
            int range_y = static_cast<int>(r_search*this->cell_size_inv.y) + 1;
            int range_z = static_cast<int>(r_search*this->cell_size_inv.z) + 1;
            index_type sight = this->_getIndex(pos);

            for(int iz=sight.z-range_z; iz<=sight.z+range_z; ++iz){
                for(int iy=sight.y-range_y; iy<=sight.y+range_y; ++iy){
                    for(int ix=sight.x-range_x; ix<=sight.x+range_x; ++ix){
                        this->buff_index.push_back( index_type{ (ix + this->nx)%(this->nx),
                                                                (iy + this->ny)%(this->ny),
                                                                (iz + this->nz)%(this->nz) } );
                    }
                }
            }
            return this->buff_index;
        }

        void get_data_list(const PS::F32vec              &pos,
                           const PS::F32                  r_search,
                                 std::vector<index_type> &index_list,
                                 std::vector<value_type> &value_list ) const {
            index_list.clear();
            value_list.clear();
            this->get_index_list(pos, r_search, index_list);
            for(const auto& index : index_list){
                for(const auto& elem : this->cell[ this->_access(index) ] ){
                    value_list.push_back(elem);
                }
            }
        }
        const std::vector<value_type>& get_data_list(const PS::F32vec &pos,
                                                     const PS::F32     r_search){
            this->buff_value.clear();
            const auto& index_list = this->get_index_list(pos, r_search);
            for(const auto& index : index_list){
                for(const auto& elem : this->cell[ this->_access(index) ] ){
                    this->buff_value.push_back(elem);
                }
            }
            return this->buff_value;
        }

    };

}
