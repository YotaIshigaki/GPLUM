// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>

//#define USE_NONDETERMINISTIC_SEED
//#define CHECK_CALC_POS_ROOT_CELL

inline PS::U64 clearSignBit(const PS::U64 in) {
    return in & 0x7fffffffffffffff; // 01111111 11111111 11111111 11111111 11111111 11111111 11111111 11111111

}
inline void clearSignBit(PS::U64 *in) {
    *in &= 0x7fffffffffffffff; // 01111111 11111111 11111111 11111111 11111111 11111111 11111111 11111111
}

inline PS::U64 clearExponentBits(const PS::U64 in) {
    return in & 0x800fffffffffffff; // 10000000 00001111 11111111 11111111 11111111 11111111 11111111 11111111
}
inline void clearExponentBits(PS::U64 *in) {
    *in &= 0x800fffffffffffff; // 10000000 00001111 11111111 11111111 11111111 11111111 11111111 11111111
}

inline PS::U64 clearFractionBits(const PS::U64 in) {
    return in & 0xfff0000000000000; // 11111111 11110000 00000000 00000000 00000000 00000000 00000000 00000000
}
inline void clearFractionBits(PS::U64 *in) {
    *in &= 0xfff0000000000000; // 11111111 11110000 00000000 00000000 00000000 00000000 00000000 00000000
}

inline PS::U64 getSignBit(const PS::U64 in) {
    return in >> 63;
}
inline PS::U64 getExponentBits(const PS::U64 in) {
    return (in << 1) >> 53;
}
inline PS::U64 getFractionBits(const PS::U64 in) {
    return (in << 12) >> 12;
}

inline PS::F64 getNextSmallValue(PS::F64 in) {
    // See https://ja.wikipedia.org/wiki/IEEE_754#64ビット倍精度
    typedef union {
        PS::F64 f64;
        PS::U64 u64;
    } BitData;
    BitData tmp, ret;
    tmp.f64 = in;
    PS::U64 sign = getSignBit(tmp.u64);
    PS::U64 fraction = getFractionBits(tmp.u64);
    PS::U64 exponent = getExponentBits(tmp.u64);
    PS::U64 one_frac = fraction | 0x0010000000000000;
    std::cout << "original = " << PS::GetBinString(tmp.u64) << std::endl;
    std::cout << "fraction = " << PS::GetBinString(fraction) << std::endl;
    std::cout << "exponent = " << PS::GetBinString(exponent) << std::endl;
    std::cout << "1.frac   = " << PS::GetBinString(one_frac) << std::endl;
    ret.f64 = 0.0;
    if ((exponent > 0) || (fraction > 0)) {
        // In this case, non-zero value is input.
        if (sign) {
            // In this case, the input value is negative.
            if (exponent > 0) {
                // In this case, the input value is a normalized number.
                one_frac++;
                if (one_frac & 0x0020000000000000) {
                    // In this case, 54th bit becomes 1.
                    // In order to make a returned value be a normalized number,
                    // we must adjust exponent.
                    std::cout << "passing through Route (1-a)..." << std::endl;
                    exponent++;
                    one_frac >>= 1;
                    ret.u64 = (sign << 63) | (exponent << 52) | (one_frac & (~0x0010000000000000));
                } else {
                    // In this case, 54 and 53th bits keep 01.
                    std::cout << "passing through Route (1-b)..." << std::endl;
                    ret.u64 = (sign << 63) | (exponent << 52) | (one_frac & (~0x0010000000000000));
                }
            } else {
                // In this case, the input value is a denormalized number.
                std::cout << "passing through Route (1-c)..." << std::endl;
                fraction++;
                if (fraction & 0x0010000000000000) {
                    ret.u64 = (sign << 63) | 0x0010000000000000;
                } else {
                    ret.u64 = (sign << 63) | fraction;
                }
            }
        } else {
            // In this case, the input value is positive.
            if (exponent > 0) {
                // In this case, the input value is a normalized number. 
                one_frac--;
                if (one_frac & 0x0010000000000000) {
                   // In this case, 53th bit is still 1.
                   // We do not have to take care the exponent.
                   std::cout << "passing through Route (2-a)..." << std::endl;
                   ret.u64 = (exponent << 52) | (one_frac & (~0x0010000000000000));
                   std::cout << "ret      = " << PS::GetBinString(ret.u64) << std::endl;
                } else {
                   // In this case, 53th bit becomes 0.
                   // In order to make a returned value be a normalized number, 
                   // we must adjust expoent.
                   std::cout << "passing through Route (2-b)..." << std::endl;
                   PS::S32 nshift = 0;
                   while (!(one_frac & 0x0010000000000000)) {
                      nshift++;
                      one_frac <<= 1;
                   }
                   assert(nshift == 1);
                   exponent -= nshift;
                   one_frac &= 0x000fffffffffffff;
                   ret.u64 = (exponent << 52) | one_frac;
                   std::cout << "ret      = " << PS::GetBinString(ret.u64) << std::endl;
                }
            } else {
                // In this case, the input value is a denormalized number.
                std::cout << "passing through Route (2-c)..." << std::endl; 
                fraction--;
                ret.u64 = fraction;
                std::cout << "ret      = " << PS::GetBinString(ret.u64) << std::endl;
            }
        }
    } else {
        // In this case, the input value is zero. The function must return a negative
        // number whose absolute value is the minimum. The sign bit of a returned value
        // is forced to be 1 (the minus sign) and its fraction bits is also 1.
        // (A returned value is the denormalized number)
        std::cout << "passing through Route (3)..." << std::endl; 
        ret.u64 = 1 | 0x8000000000000000;
        std::cout << "ret      = " << PS::GetBinString(ret.u64) << std::endl;
    }
    return ret.f64;
}

PS::S32 getMinExpOfTwoGE(const PS::S32 in) { // GE := equal to or greater than
    assert(in > 0);
    PS::S32 tmp = 1, lev = 0;
    while (tmp < in) {
        tmp *= 2;
        lev++;
    }
    return lev;
}

PS::S32 getMinExpOfTwoGT(const PS::S32 in) { // GT := greater than
    assert(in > 0);
    PS::S32 tmp = 1, lev = 0;
    while (tmp <= in) {
        tmp *= 2;
        lev++;
    }
    return lev;
}

void calcPosRootCell(const PS::S32vec nc_root_domain,
                     const PS::F64ort pos_root_domain,
                     const PS::S32 icut,
                     PS::S32 level,
                     PS::F64ort pos_root_cell) {
    // Cell width calculated by pos_root_domain & nc_root_domain
    PS::F64vec cwid_root_domain;
    cwid_root_domain.x = pos_root_domain.getFullLength().x / nc_root_domain.x;
    cwid_root_domain.y = pos_root_domain.getFullLength().y / nc_root_domain.y;
    cwid_root_domain.z = pos_root_domain.getFullLength().z / nc_root_domain.z;
#ifdef CHECK_CALC_POS_ROOT_CELL
    std::cout << "cwid_root_domain = " << cwid_root_domain << std::endl;
#endif

    // Compute level & nc_root_cell
    PS::S32vec nc_root_cell_reqmin, lv_root_cell;
    nc_root_cell_reqmin.x = nc_root_domain.x + 2*icut;
    nc_root_cell_reqmin.y = nc_root_domain.y + 2*icut;
    nc_root_cell_reqmin.z = nc_root_domain.z + 2*icut;
#ifdef CHECK_CALC_POS_ROOT_CELL
    std::cout << "nc_root_cell_reqmin = " << nc_root_cell_reqmin << std::endl;
#endif
    lv_root_cell.x = getMinExpOfTwoGT(nc_root_cell_reqmin.x);
    lv_root_cell.y = getMinExpOfTwoGT(nc_root_cell_reqmin.y);
    lv_root_cell.z = getMinExpOfTwoGT(nc_root_cell_reqmin.z);
    level = std::max(lv_root_cell.x, std::max(lv_root_cell.y, lv_root_cell.z));
#ifdef CHECK_CALC_POS_ROOT_CELL
    std::cout << "level = " << level << std::endl;
#endif
    PS::S32vec nc_root_cell;
    nc_root_cell.x = nc_root_cell.y = nc_root_cell.z = std::pow(2,level);
#ifdef CHECK_CALC_POS_ROOT_CELL
    std::cout << "nc_root_cell = " << nc_root_cell << std::endl;
#endif

    // Compute idx_root_domain & pos_root_cell
    PS::S32vec idx_root_domain;
    idx_root_domain = (nc_root_cell - nc_root_cell_reqmin)/2 + icut;
#ifdef CHECK_CALC_POS_ROOT_CELL
    std::cout << "idx_root_domain = " << idx_root_domain << std::endl;
#endif
    // (1) initial guess of pos_root_cell
    pos_root_cell.low_.x = pos_root_domain.low_.x - cwid_root_domain.x * idx_root_domain.x;
    pos_root_cell.low_.y = pos_root_domain.low_.y - cwid_root_domain.y * idx_root_domain.y;
    pos_root_cell.low_.z = pos_root_domain.low_.z - cwid_root_domain.z * idx_root_domain.z;
    pos_root_cell.high_.x = pos_root_cell.low_.x + cwid_root_domain.x * nc_root_cell.x;
    pos_root_cell.high_.y = pos_root_cell.low_.y + cwid_root_domain.y * nc_root_cell.y;
    pos_root_cell.high_.z = pos_root_cell.low_.z + cwid_root_domain.z * nc_root_cell.z;
#ifdef CHECK_CALC_POS_ROOT_CELL
    std::cout << "pos_root_cell (ini. guess) = " << pos_root_cell << std::endl;    
#endif
    // (2) correct pos_root_cell
    const PS::S32 itermax = 2;
    PS::S32 iter;
    PS::F64vec cwid_root_cell;
    PS::F64vec pos, pos_ov_cwid_low, pos_ov_cwid_high;
    PS::F64vec shift;
    for (iter=0; iter < itermax; iter++) {
        std::cout << "iter = " << iter << std::endl;
        // update pos_root_cell
        pos_root_cell.low_.x -= cwid_root_cell.x * shift.x;
        pos_root_cell.low_.y -= cwid_root_cell.y * shift.y;
        pos_root_cell.low_.z -= cwid_root_cell.z * shift.z;
        pos_root_cell.high_.x = pos_root_cell.low_.x + cwid_root_domain.x * nc_root_cell.x;
        pos_root_cell.high_.y = pos_root_cell.low_.y + cwid_root_domain.y * nc_root_cell.y;
        pos_root_cell.high_.z = pos_root_cell.low_.z + cwid_root_domain.z * nc_root_cell.z;
        // cwid_root_cell 
        cwid_root_cell.x = pos_root_cell.getFullLength().x / nc_root_cell.x;
        cwid_root_cell.y = pos_root_cell.getFullLength().y / nc_root_cell.y;
        cwid_root_cell.z = pos_root_cell.getFullLength().z / nc_root_cell.z;
#ifdef CHECK_CALC_POS_ROOT_CELL
        std::cout << "cwid_root_cell = " << cwid_root_cell << std::endl;
#endif
        if ((cwid_root_cell.x != cwid_root_domain.x) ||
            (cwid_root_cell.y != cwid_root_domain.y) ||
            (cwid_root_cell.z != cwid_root_domain.z)) {
            std::cerr << "cwid_root_cell != cwid_root_domain occurs!!" << std::endl;
            std::cerr << "cwid_root_cell   = " << cwid_root_cell << std::endl;
            std::cerr << "cwid_root_domain = " << cwid_root_domain << std::endl;
            std::exit(EXIT_FAILURE);
        }
        // pos_low_ov_cwid & pos_high_ov_cwid
        pos = pos_root_domain.low_ - pos_root_cell.low_;
        pos_ov_cwid_low.x = pos.x/cwid_root_cell.x;
        pos_ov_cwid_low.y = pos.y/cwid_root_cell.y;
        pos_ov_cwid_low.z = pos.z/cwid_root_cell.z;
        pos = pos_root_domain.high_ - pos_root_cell.low_;
        pos_ov_cwid_high.x = pos.x/cwid_root_cell.x;
        pos_ov_cwid_high.y = pos.y/cwid_root_cell.y;
        pos_ov_cwid_high.z = pos.z/cwid_root_cell.z;
        shift.x = pos_ov_cwid_low.x - idx_root_domain.x;
        shift.y = pos_ov_cwid_low.y - idx_root_domain.y;
        shift.z = pos_ov_cwid_low.z - idx_root_domain.z;
        if (shift.x < 0) shift.x = - 2*shift.x;
        else shift.x = 0.0;
        if (shift.y < 0) shift.y = - 2*shift.y;
        else shift.y = 0.0;
        if (shift.z < 0) shift.z = - 2*shift.z;
        else shift.z = 0.0;
        // check 
        if ((shift.x == 0.0) && (shift.y == 0.0) && (shift.z == 0.0))  break;
    }
    PS::S32vec idx ((PS::S32)pos_ov_cwid_low.x,
                    (PS::S32)pos_ov_cwid_low.y,
                    (PS::S32)pos_ov_cwid_low.z);
#ifdef CHECK_CALC_POS_ROOT_CELL
    std::cout << "idx (recalc.) = " << idx << std::endl;
#endif
    if ((idx.x != idx_root_domain.x) ||
        (idx.y != idx_root_domain.y) ||
        (idx.z != idx_root_domain.z)) {
        PS::F64vec ipart, fpart;
        std::cerr << "idx != idx_root_domain occurs!!" << std::endl;
        std::cerr << "idx              = " << idx << std::endl;
        std::cerr << "idx_root_domain  = " << idx_root_domain << std::endl;
        std::cerr << "pos_ov_cwid_low  = " << pos_ov_cwid_low << std::endl;
        std::cerr << "  x comp.        = " << PS::GetBinString((PS::U64 *)&pos_ov_cwid_low.x) << std::endl;
        std::cerr << "  y comp.        = " << PS::GetBinString((PS::U64 *)&pos_ov_cwid_low.y) << std::endl;
        std::cerr << "  z comp.        = " << PS::GetBinString((PS::U64 *)&pos_ov_cwid_low.z) << std::endl;
        fpart.x = std::modf(pos_ov_cwid_low.x, &ipart.x);
        fpart.y = std::modf(pos_ov_cwid_low.y, &ipart.y);
        fpart.z = std::modf(pos_ov_cwid_low.z, &ipart.z);
        std::cerr << "  integral part  = " << ipart << std::endl;
        std::cerr << "  fraction part  = " << fpart << std::endl;
        std::cerr << "pos_ov_cwid_high = " << pos_ov_cwid_high << std::endl;
        std::cerr << "  x comp.        = " << PS::GetBinString((PS::U64 *)&pos_ov_cwid_high.x) << std::endl;
        std::cerr << "  y comp.        = " << PS::GetBinString((PS::U64 *)&pos_ov_cwid_high.y) << std::endl;
        std::cerr << "  z comp.        = " << PS::GetBinString((PS::U64 *)&pos_ov_cwid_high.z) << std::endl;
        fpart.x = std::modf(pos_ov_cwid_high.x, &ipart.x);
        fpart.y = std::modf(pos_ov_cwid_high.y, &ipart.y);
        fpart.z = std::modf(pos_ov_cwid_high.z, &ipart.z);
        std::cerr << "  integral part  = " << ipart << std::endl;
        std::cerr << "  fraction part  = " << fpart << std::endl;
        std::exit(EXIT_FAILURE);
    }

  
    
}

void calcTargetBox(const PS::F64ort pos_root_cell,
                   const PS::S32 level,
                   const PS::F64ort target_box_in,
                   PS::F64ort target_box) {

}

void calcPosPMMCell(const PS::S32vec n_cell_pmm,
                    const PS::F64ort pos_root_domain,
                    const PS::S32vec idx,
                    PS::F64ort pos_cell) {

}


int main(int argc, char *argv[]) {
    std::cout << std::setprecision(15);

#ifdef USE_NONDETERMINISTIC_SEED
    std::random_device rnd;
    std::mt19937_64 mt(rnd());
#else
    const PS::S32 size_seed_vec = 10;
    std::vector<uint32_t> rnd_seed_vec(size_seed_vec);
    for (PS::S32 i=0; i<size_seed_vec; i++) rnd_seed_vec[i]=i;
    std::seed_seq rnd_seed(rnd_seed_vec.begin(), rnd_seed_vec.end());
    std::mt19937_64 mt(rnd_seed);
#endif
    std::uniform_int_distribution<PS::S32> dist_i(4, 512);
    std::uniform_real_distribution<PS::F64> dist_r(0.0, 1.0);
      
#if 0 
    const PS::S32 nmax = 4;
    for (PS::S32 n = 0; n < nmax; n++) {
        std::cout << "n = " << n << std::endl;

        // Set computational box (real image)
        const PS::S32 nc1 (dist_i(mt));
        const PS::S32vec nc_root_domain (PS::S32vec(nc1,nc1,nc1));
        const PS::F64vec low (-dist_r(mt), -dist_r(mt), -dist_r(mt));
        const PS::F64 d = 2.0*dist_r(mt);
        const PS::F64vec high (PS::F64vec(low.x + d, low.y + d, low.z + d));
        const PS::F64ort pos_root_domain (low, high);
        std::cout << "   nc_root_domain = " << nc_root_domain << std::endl;
        std::cout << "   pos_root_domain = " << pos_root_domain << std::endl;

        // Calculate level & pos_root_cell
        const PS::S32 icut = 2;
        PS::S32 level;
        PS::F64ort pos_root_cell;
        calcPosRootCell(nc_root_domain,
                        pos_root_domain,
                        icut,
                        level,
                        pos_root_cell);
        
    }
#endif

    const PS::U64 one = 1;
    std::cout << "one = " << PS::GetBinString(one) << std::endl;
    PS::U64 tmp = std::numeric_limits<unsigned long long int>::max();
    //std::cout << "tmp = " << PS::GetBinString(tmp) << std::endl;
    //clearSignBit(&tmp);
    //clearExponentBits(&tmp);
    //clearFractionBits(&tmp);
    //std::cout << "tmp = " << PS::GetBinString(tmp) << std::endl;

    //const PS::F64 x = 1.2;
    const PS::F64 x = 1.0;
    const PS::F64 y = getNextSmallValue(x);
    std::cout << "x     = " << x << std::endl;
    std::cout << "y     = " << y << std::endl;
    std::cout << "x - y = " << x - y << std::endl;


    return 0;
}
