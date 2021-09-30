#pragma once

#include<iostream>
#include<iomanip>
#include"vector3.hpp"
#include"ps_defs.hpp"

namespace ParticleSimulator{

    /*
      MSB is always 0.
      next 3bits represent octant index of 8 cells with level 1
     */
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    template <typename Tkey>
    class MortonKey{};

    template <>
    class MortonKey<TagKeyNormal>{
    private:
        enum{
            kLevMax = 30,
        };
        F64 half_len_;
        F64vec center_;
        F64 normalized_factor_;
        MortonKey(){};
        ~MortonKey(){};
        MortonKey(const MortonKey &){};
        MortonKey & operator = (const MortonKey &);
        static MortonKey & getInstance(){
            static MortonKey inst;
            return inst;
        }
        static U64 separateBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | _s<<32) &  0x00000000ffffffff;  //0000 0000 0000 0000 0000 0000 0000 0000 1111 1111 1111 1111 1111 1111 1111 1111 
            _s = (_s | _s<<16) &  0x0000ffff0000ffff;  //0000 0000 0000 0000 1111 1111 1111 1111 0000 0000 0000 0000 1111 1111 1111 1111 
            _s = (_s | _s<<8) &  0x00ff00ff00ff00ff;  //0000 0000 1111 1111 0000 0000 1111 1111
            _s = (_s | _s<<4) &   0x0f0f0f0f0f0f0f0f;  //0000 1111 0000 1111 0000 1111
            _s = (_s | _s<<2) &   0x3333333333333333;  //00 11 00 11 00 11
            return (_s | _s<<1) & 0x5555555555555555;  //0101 0101
        }
    public:
        static void initialize(const F64 half_len,
                               const F64vec & center,
                               const F64vec pos_ref,
                               const S32vec idx_ref,
                               const S32    lev_ref,
                               const F64vec wid_ref) {
            getInstance().half_len_ = half_len;
            getInstance().center_ = center;
            getInstance().normalized_factor_ = (1.0/(half_len*2.0))*(1<<kLevMax);
        }
        static U64 getKey(const F64vec & pos){
            const F64vec cen = getInstance().center_;
            const F64 hlen = getInstance().half_len_;
            const F64 nfactor = getInstance().normalized_factor_;
            const U64 nx = (U64)( (pos.x - cen.x + hlen) * nfactor);
            const U64 ny = (U64)( (pos.y - cen.y + hlen) * nfactor);
            //std::cerr<<"cen="<<cen<<" hlen="<<hlen<<std::endl;
            //std::cerr<<"nx="<<nx<<" ny="<<ny<<std::endl;
            return ( getInstance().separateBit(nx)<<1 | getInstance().separateBit(ny) );
        }
        static S32 getCellID(const S32 lev, const U64 mkey){
            U64 s = mkey >> ( (kLevMax - lev) * 2 );
            return (s & 0x3);
        }
        static S32vec getPMCellID(const F64vec pos) {
            const U64 mkey = 0;
            return getPMCellID(mkey);
        }
        static S32vec getPMCellID(const U64 mkey) {
            // just return a dummy value.
            return S32vec(-1, -1, -1);
        }
    };

    template <>
    class MortonKey<TagKeyMeshBased>{
    private:
        enum{
            kLevMax = 30,
        };
        F64vec pos_ref_;
        F64vec wid_ref_;
        F64vec wid_btm_;
        S32vec idx_ref_;
        S32    lev_ref_;
        MortonKey(){};
        ~MortonKey(){};
        MortonKey(const MortonKey &);
        MortonKey & operator = (const MortonKey &);
        static MortonKey & getInstance(){
            static MortonKey inst;
            return inst;
        }
        static U64 separateBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | _s<<32) &  0x00000000ffffffff;  //0000 0000 0000 0000 0000 0000 0000 0000 1111 1111 1111 1111 1111 1111 1111 1111 
            _s = (_s | _s<<16) &  0x0000ffff0000ffff;  //0000 0000 0000 0000 1111 1111 1111 1111 0000 0000 0000 0000 1111 1111 1111 1111 
            _s = (_s | _s<<8) &  0x00ff00ff00ff00ff;  //0000 0000 1111 1111 0000 0000 1111 1111
            _s = (_s | _s<<4) &   0x0f0f0f0f0f0f0f0f;  //0000 1111 0000 1111 0000 1111
            _s = (_s | _s<<2) &   0x3333333333333333;  //00 11 00 11 00 11
            return (_s | _s<<1) & 0x5555555555555555;  //0101 0101
        }
        static U64 shrinkBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | (_s >>  1)) & 0x3333333333333333ULL;
            _s = (_s | (_s >>  2)) & 0x0f0f0f0f0f0f0f0fULL;
            _s = (_s | (_s >>  4)) & 0x00ff00ff00ff00ffULL;
            _s = (_s | (_s >>  8)) & 0x0000ffff0000ffffULL;
            _s = (_s | (_s >> 16)) & 0x00000000ffffffffULL;
            return _s;
        }
    public:
        static void initialize(const F64 half_len,
                               const F64vec & center,
                               const F64vec pos_ref,
                               const S32vec idx_ref,
                               const S32    lev_ref,
                               const F64vec wid_ref) {
            getInstance().pos_ref_ = pos_ref;
            getInstance().idx_ref_ = idx_ref;
            getInstance().wid_ref_ = wid_ref;
            getInstance().lev_ref_ = lev_ref;
            getInstance().wid_btm_ = wid_ref / (double)(1<<(kLevMax - lev_ref));
        }
        static U64 getKey(F64vec pos) {
            const F64vec pos_ref = getInstance().pos_ref_;
            const S32vec idx_ref = getInstance().idx_ref_;
            const F64vec wid_ref = getInstance().wid_ref_;
            const F64vec wid_btm = getInstance().wid_btm_;
            const S32    lev_ref = getInstance().lev_ref_;
            // Compute cell index 
            F64vec idx_f;
            S32 ix, iy;
            idx_f.x = (pos.x - pos_ref.x) / wid_ref.x;
            idx_f.y = (pos.y - pos_ref.y) / wid_ref.y;
            assert(fabs(idx_f.x) < (F64)std::numeric_limits<S32>::max());
            assert(fabs(idx_f.y) < (F64)std::numeric_limits<S32>::max());
            if (idx_f.x >= 0.0) ix = (S32) idx_f.x;
            else ix = ((S32) idx_f.x) - 1;
            if (idx_f.y >= 0.0) iy = (S32) idx_f.y;
            else iy = ((S32) idx_f.y) - 1;
            S32 ix_btm, iy_btm;
            if (kLevMax > lev_ref) {
                // Compute sub-cell index in cell (ix,iy,iz).
                F64vec pos_cell;
                pos_cell.x = pos_ref.x + wid_ref.x * ix;
                pos_cell.y = pos_ref.y + wid_ref.y * iy;
                // [Notes]
                //    In principal, we can choose pos_cell and wid_btm
                //    different from the those we are currently adopting
                //    because it is guaranteed that the particle of interest
                //    belongs to the cell to which it should belong.
                //    The only condition imposed on sub-cell indices is that
                //    the order of the particles belonging to that cell is
                //    preserved.
                idx_f.x = (pos.x - pos_cell.x)/ wid_btm.x;
                idx_f.y = (pos.y - pos_cell.y)/ wid_btm.y;
                if (idx_f.x >= 0.0) ix_btm = (S32) idx_f.x;
                else ix_btm = ((S32) idx_f.x) - 1;
                if (idx_f.y >= 0.0) iy_btm = (S32) idx_f.y;
                else iy_btm = ((S32) idx_f.y) - 1;
                const S32 maxval = (1<<(kLevMax - lev_ref)) - 1;
                if (ix_btm < 0) ix_btm = 0;
                if (ix_btm > maxval) ix_btm = maxval;
                if (iy_btm < 0) iy_btm = 0;
                if (iy_btm > maxval) iy_btm = maxval;
            } else {
                ix_btm = iy_btm = 0;
            }
            U64 nx = (U64)((ix + idx_ref.x) * (1<<(kLevMax - lev_ref)) + ix_btm);
            U64 ny = (U64)((iy + idx_ref.y) * (1<<(kLevMax - lev_ref)) + iy_btm);
            return ( getInstance().separateBit(nx)<<1 | getInstance().separateBit(ny) );
        }
        static S32 getCellID(const S32 lev, const U64 mkey){
            U64 s = mkey >> ( (kLevMax - lev) * 2 );
            return (s & 0x3);
        }
        static S32vec getPMCellID(const F64vec pos) {
            const U64 mkey = getKey(pos);
            return getPMCellID(mkey);
        }
        static S32vec getPMCellID(const U64 mkey) {
            const S32vec idx_ref = getInstance().idx_ref_;
            const S32    lev_ref = getInstance().lev_ref_;
            const U64 x = (mkey & 0xaaaaaaaaaaaaaaaaULL)>>1; // Repeat of pattern 10
            const U64 y =  mkey & 0x5555555555555555ULL    ; // Repeat of pattern 01
            U64 ix = (shrinkBit(x)>>(kLevMax - lev_ref)) - idx_ref.x;
            U64 iy = (shrinkBit(y)>>(kLevMax - lev_ref)) - idx_ref.y;
            return S32vec((S32)ix, (S32)iy);
        }
    };

#else
    template <typename Tkey>
    class MortonKey{};

    template <>
    class MortonKey<TagKeyNormal>{
    private:
        enum{
            kLevMax = 21,
        };
        F64 half_len_;
        F64vec center_;
        F64 normalized_factor_;
        MortonKey(){};
        ~MortonKey(){};
        MortonKey(const MortonKey &);
        MortonKey & operator = (const MortonKey &);
        static MortonKey & getInstance(){
            static MortonKey inst;
            return inst;
        }
        static U64 separateBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | _s<<32) & 0xffff00000000ffff; //11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
            _s = (_s | _s<<16) & 0x00ff0000ff0000ff; //00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
            _s = (_s | _s<<8) & 0xf00f00f00f00f00f;  //1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
            _s = (_s | _s<<4) & 0x30c30c30c30c30c3;  //11 00 00 11 00 11 00 00 11
            return (_s | _s<<2) & 0x9249249249249249;  //1 0 0 1 0 0 1 0 0 1 0 0 1
        }
    public:
        static void initialize(const F64 half_len,
                               const F64vec & center,
                               const F64vec & pos_ref,
                               const S32vec & idx_ref,
                               const S32    lev_ref,
                               const F64vec & wid_ref){
            getInstance().half_len_ = half_len;
            getInstance().center_ = center;
            getInstance().normalized_factor_ = (1.0/(half_len*2.0))*(1<<kLevMax);
        }
        static U64 getKey(F64vec pos){
            const F64vec cen = getInstance().center_;
            const F64 hlen = getInstance().half_len_;
            const F64 nfactor = getInstance().normalized_factor_;
            U64 nx = (U64)( (pos.x - cen.x + hlen) * nfactor);
            U64 ny = (U64)( (pos.y - cen.y + hlen) * nfactor);
            U64 nz = (U64)( (pos.z - cen.z + hlen) * nfactor);
            return ( getInstance().separateBit(nx)<<2 | getInstance().separateBit(ny)<<1 | getInstance().separateBit(nz) );
        }
        static S32 getCellID(const S32 lev, const U64 mkey){
            U64 s = mkey >> ( (kLevMax - lev) * 3 );
            return (s & 0x7);
        }
        static S32vec getPMCellID(const F64vec pos) {
            const U64 mkey = 0;
            return getPMCellID(mkey);
        }
        static S32vec getPMCellID(const U64 mkey) {
            // just return a dummy value.
            return S32vec(-1, -1, -1);
        }
    };
    
    template <>
    class MortonKey<TagKeyMeshBased>{
    private:
        enum{
            kLevMax = 21,
        };
        F64vec pos_ref_;
        F64vec wid_ref_;
        F64vec wid_btm_;
        S32vec idx_ref_;
        S32    lev_ref_;
        MortonKey(){};
        ~MortonKey(){};
        MortonKey(const MortonKey &);
        MortonKey & operator = (const MortonKey &);
        static MortonKey & getInstance(){
            static MortonKey inst;
            return inst;
        }
        static U64 separateBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | _s<<32) & 0xffff00000000ffff; //11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
            _s = (_s | _s<<16) & 0x00ff0000ff0000ff; //00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
            _s = (_s | _s<<8) & 0xf00f00f00f00f00f;  //1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
            _s = (_s | _s<<4) & 0x30c30c30c30c30c3;  //11 00 00 11 00 11 00 00 11
            return (_s | _s<<2) & 0x9249249249249249;  //1 0 0 1 0 0 1 0 0 1 0 0 1
        }
        static U64 shrinkBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | _s>>2)  & 0x30c30c30c30c30c3ULL; // Repeat of pattern 00 00 11
            _s = (_s | _s>>4)  & 0xf00f00f00f00f00fULL; // Repeat of pattern 0000 0000 1111
            _s = (_s | _s>>8)  & 0x00ff0000ff0000ffULL; // Repeat of pattern 00000000 00000000 11111111
            _s = (_s | _s>>16) & 0xffff00000000ffffULL; // Repeat of pattern 0^{16} 0^{16} 1^{16},
                                                        // where x^{n} reprensets that x continues n times.
            _s = (_s | _s>>32) & 0x00000000ffffffffULL; // Repaat of pattern 0^{32} 0^{32} 1^{32}
            return _s;
        }
    public:
        static void initialize(const F64 half_len,
                               const F64vec & center,
                               const F64vec & pos_ref,
                               const S32vec & idx_ref,
                               const S32    lev_ref,
                               const F64vec & wid_ref) {
            getInstance().pos_ref_ = pos_ref;
            getInstance().idx_ref_ = idx_ref;
            getInstance().wid_ref_ = wid_ref;
            getInstance().lev_ref_ = lev_ref;
            getInstance().wid_btm_ = wid_ref / (double)(1<<(kLevMax - lev_ref));
        }
        static U64 getKey(F64vec pos) {
            const F64vec pos_ref = getInstance().pos_ref_;
            const S32vec idx_ref = getInstance().idx_ref_;
            const F64vec wid_ref = getInstance().wid_ref_;
            const F64vec wid_btm = getInstance().wid_btm_;
            const S32    lev_ref = getInstance().lev_ref_;
            // Compute cell index 
            S32vec idx_pm = GetCellIDMeasuredInUnitCell(pos_ref, wid_ref, pos);
            S32vec idx_btm;
            if (kLevMax > lev_ref) {
                // Compute sub-cell index in cell idx_pm.
                F64vec pos_cell;
                pos_cell.x = pos_ref.x + wid_ref.x * idx_pm.x;
                pos_cell.y = pos_ref.y + wid_ref.y * idx_pm.y;
                pos_cell.z = pos_ref.z + wid_ref.z * idx_pm.z;
                // [Notes]
                //    In principal, we can choose pos_cell and wid_btm
                //    different from the those we are currently adopting
                //    because it is guaranteed that the particle of interest
                //    belongs to the cell to which it should belong.
                //    The only condition imposed on sub-cell indices is that
                //    the order of the particles belonging to that cell is
                //    preserved.
                idx_btm = GetCellIDMeasuredInUnitCell(pos_cell, wid_btm, pos);
                const S32 maxval = (1<<(kLevMax - lev_ref)) - 1;
                if (idx_btm.x < 0) idx_btm.x = 0;
                if (idx_btm.x > maxval) idx_btm.x = maxval;
                if (idx_btm.y < 0) idx_btm.y = 0;
                if (idx_btm.y > maxval) idx_btm.y = maxval;
                if (idx_btm.z < 0) idx_btm.z = 0;
                if (idx_btm.z > maxval) idx_btm.z = maxval;
            } else {
                idx_btm.x = idx_btm.y = idx_btm.z = 0;
            }
            U64 nx = (U64)((idx_pm.x + idx_ref.x) * (1<<(kLevMax - lev_ref)) + idx_btm.x);
            U64 ny = (U64)((idx_pm.y + idx_ref.y) * (1<<(kLevMax - lev_ref)) + idx_btm.y);
            U64 nz = (U64)((idx_pm.z + idx_ref.z) * (1<<(kLevMax - lev_ref)) + idx_btm.z);
#if 0
            //---- a code to test getPMCellID [start]
            U64 mkey = (getInstance().separateBit(nx)<<2 |
                        getInstance().separateBit(ny)<<1 | 
                        getInstance().separateBit(nz));
            S32vec idx = getPMCellID(mkey);
            std::cout << "(ix, iy, iz) [org.]    = " << ix << ", " << iy << ", " << iz << std::endl;
            std::cout << "(ix, iy, iz) [decoded] = " << idx.x << ", " << idx.y << ", " << idx.z << std::endl;
            if ((ix != idx.x) || (iy != idx.y) || (iz != idx.z)) {
                Abort(-1);
            }
            //---- a code to test getPMCellID [end]
#endif
            return ( getInstance().separateBit(nx)<<2 |
                     getInstance().separateBit(ny)<<1 | 
                     getInstance().separateBit(nz) );
        }
        static S32 getCellID(const S32 lev, const U64 mkey){
            U64 s = mkey >> ( (kLevMax - lev) * 3 );
            return (s & 0x7);
        }
        static S32vec getPMCellID(const F64vec pos) {
            const U64 mkey = getKey(pos);
            return getPMCellID(mkey);
        }
        static S32vec getPMCellID(const U64 mkey) {
            const S32vec idx_ref = getInstance().idx_ref_;
            const S32    lev_ref = getInstance().lev_ref_;
            const U64 x = (mkey & 0x4924924924924924ULL)>>2; // Repeat of pattern 100
            const U64 y = (mkey & 0x2492492492492492ULL)>>1; // Repeat of pattern 010
            const U64 z =  mkey & 0x1249249249249249ULL;     // Repeat of pattern 001
            U64 ix = (shrinkBit(x)>>(kLevMax - lev_ref)) - idx_ref.x;
            U64 iy = (shrinkBit(y)>>(kLevMax - lev_ref)) - idx_ref.y;
            U64 iz = (shrinkBit(z)>>(kLevMax - lev_ref)) - idx_ref.z;
            return S32vec((S32)ix, (S32)iy, (S32)iz);
        }
    };
#endif
    
}

