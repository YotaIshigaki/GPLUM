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
    class MortonKey{
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
                               const F64vec & center=0.0){
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
    };


#else
    class MortonKey{
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
                               const F64vec & center=0.0){
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
        
        static F64ort getCorrespondingTreeCell(const F64ort & box){
            const U64 low  = getKey(box.low_);
            const U64 high = getKey(box.high_);
            const U64 tmp = high^low;
            S32 lev = 0;
            for(S32 i=TREE_LEVEL_LIMIT-1; i>=0; i--){
                if( ((tmp >> i*3) & 0x7) != 0){
                    lev = TREE_LEVEL_LIMIT-i-1;
                    break;
                }
            }
            F64ort ret;
            const F64vec cen = getInstance().center_;
            const F64 hlen = getInstance().half_len_;
            if(lev==0){
                ret.low_.x  = cen.x - hlen;
                ret.low_.y  = cen.y - hlen;
                ret.low_.z  = cen.z - hlen;
                ret.high_.x = cen.x + hlen;
                ret.high_.y = cen.y + hlen;
                ret.high_.z = cen.z + hlen;
            }
            else{
                F64vec cen_new = cen;
                F64 hlen_new = hlen;
                for(S32 i=0; i<lev; i++){
                    S32 id = getCellID(i+1, low);
                    cen_new +=  hlen_new*SHIFT_CENTER[id];
                    hlen_new *= 0.5;
                }
                ret.low_.x  = cen_new.x - hlen_new;
                ret.low_.y  = cen_new.y - hlen_new;
                ret.low_.z  = cen_new.z - hlen_new;
                ret.high_.x = cen_new.x + hlen_new;
                ret.high_.y = cen_new.y + hlen_new;
                ret.high_.z = cen_new.z + hlen_new;
            }
            return ret;
        }
        
    };
    
#endif
    
}

