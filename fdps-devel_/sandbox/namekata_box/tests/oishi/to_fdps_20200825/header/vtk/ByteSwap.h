/*
 * File:   ByteSwap.h
 * Ver1:   by Lalith @ U. of Tokyo
 * Ver2:   Copied from IES_SRA and simplified for FDPS-SPH by J. Chen @AICS Nov 29, 2017
 */

#ifndef NS_BYTESWAP_H_INCLUDED
#define NS_BYTESWAP_H_INCLUDED
#include <iostream>
#include <climits>
#include <cstdlib>
////////////////////////////////////////////////////////////////////////////////
//
// --- function
//
#ifndef Swap2Bytes
#define Swap2Bytes(data)  ( (((data)>>8)& 0x00FF) | (((data)<<8)& 0xFF00) )
#endif
#ifndef Swap4Bytes
#define Swap4Bytes(data)  ( (((data)>>24)& 0x000000FF) | (((data)>> 8)& 0x0000FF00) | \
                            (((data)<< 8)& 0x00FF0000) | (((data)<<24)& 0xFF000000) )
#endif
#ifndef Swap8Bytes
#define Swap8Bytes(data)  ( (((data) >> 56) & 0x00000000000000FF) | (((data) >> 40) & 0x000000000000FF00) | \
                            (((data) >> 24) & 0x0000000000FF0000) | (((data) >> 8 ) & 0x00000000FF000000) | \
                            (((data) << 8 ) & 0x000000FF00000000) | (((data) << 24) & 0x0000FF0000000000) | \
                            (((data) << 40) & 0x00FF000000000000) | (((data) << 56) & 0xFF00000000000000) )
#endif

#undef SIZE_64_UINT
#if      UINT_MAX   == 0xffffffffffffffff
 #define  SIZE_64_UINT  unsigned int
#elif    ULONG_MAX  == 0xffffffffffffffff
 #define  SIZE_64_UINT  unsigned long
#elif    ULLONG_MAX == 0xffffffffffffffff
 #define  SIZE_64_UINT  unsigned long long
#else
 #error  no 64bit primitive integer
#endif

#undef SIZE_32_UINT
#if      UINT_MAX   == 0xffffffff
 #define  SIZE_32_UINT  unsigned int
#elif    ULONG_MAX  == 0xffffffff
 #define  SIZE_32_UINT  unsigned long
#elif    ULLONG_MAX == 0xffffffff
 #define  SIZE_32_UINT  unsigned long long
#else
 #error  no 32bit primitive integer
#endif

#undef SIZE_16_UINT
#if      USHRT_MAX   == 0xffff
 #define  SIZE_16_UINT  unsigned short
#elif    UINT_MAX    == 0xffff
 #define  SIZE_16_UINT  unsigned int
#else
 #error  no 16bit primitive integer
#endif

#ifndef  SWAP_FLOAT
#define  SWAP_FLOAT(x)  (*(SIZE_32_UINT *)&(x)=Swap4Bytes(*(SIZE_32_UINT *)&(x)))
#endif
#ifndef  SWAP_DOUBLE
#define  SWAP_DOUBLE(x) (*(SIZE_64_UINT *)&(x)=Swap8Bytes(*(SIZE_64_UINT *)&(x)))
#endif

#ifndef  SWAP_SHORT
#if      USHRT_MAX == 0xffff
 #define SWAP_SHORT(x)  (*(SIZE_16_UINT *)&(x)=Swap2Bytes(*(SIZE_16_UINT *)&(x)))
#else
 #error  sizeof(short) must be 2
#endif
#endif

#ifndef  SWAP_INT
#if      UINT_MAX == 0xffffffffffffffff
 #define SWAP_INT(x)    (*(SIZE_64_UINT *)&(x)=Swap8Bytes(*(SIZE_64_UINT *)&(x)))
#elif    UINT_MAX == 0xffffffff
 #define SWAP_INT(x)    (*(SIZE_32_UINT *)&(x)=Swap4Bytes(*(SIZE_32_UINT *)&(x)))
#elif    UINT_MAX == 0xffff
 #define SWAP_INT(x)    (*(SIZE_16_UINT *)&(x)=Swap2Bytes(*(SIZE_16_UINT *)&(x)))
#else
 #error  sizeof(int) must be 2, 4 or 8
#endif
#endif

#ifndef  SWAP_LONG
#if      ULONG_MAX == 0xffffffffffffffff
 #define SWAP_LONG(x)   (*(SIZE_64_UINT *)&(x)=Swap8Bytes(*(SIZE_64_UINT *)&(x)))
#elif    ULONG_MAX == 0xffffffff
 #define SWAP_LONG(x)   (*(SIZE_32_UINT *)&(x)=Swap4Bytes(*(SIZE_32_UINT *)&(x)))
#else
 #error  sizeof(long) must be 4 or 8
#endif
#endif

#ifdef ULLONG_MAX
#ifndef  SWAP_LLONG
#if      ULLONG_MAX == 0xffffffffffffffff
 #define SWAP_LLONG(x)  (*(SIZE_64_UINT *)&(x)=Swap8Bytes(*(SIZE_64_UINT *)&(x)))
#else
 #error  sizeof(long long) must be 8
#endif
#endif
#endif

namespace
{
  inline int GetByteOrder() // 1:big, 0:little
  { int i( 1 ); char *s = (char *)&i; return( (s == 0x00)? 1:0 ); }
  
  template<typename T> inline void SwapByte( T &val )
    { std::cerr << "error: SwapByte is not specialized for this data type" << std::endl; exit(1); }
  template <> inline void SwapByte(float  &val){ SWAP_FLOAT (val); }
  template <> inline void SwapByte(double &val){ SWAP_DOUBLE(val); }
  template <> inline void SwapByte(         char  &val){ }
  template <> inline void SwapByte(unsigned char  &val){ }
  template <> inline void SwapByte(         short &val){ SWAP_SHORT(val); }
  template <> inline void SwapByte(unsigned short &val){ SWAP_SHORT(val); }
  template <> inline void SwapByte(         int   &val){ SWAP_INT  (val); }
  template <> inline void SwapByte(unsigned int   &val){ SWAP_INT  (val); }
  template <> inline void SwapByte(         long  &val){ SWAP_LONG (val); }
  template <> inline void SwapByte(unsigned long  &val){ SWAP_LONG (val); }
  #ifdef ULLONG_MAX
  template <> inline void SwapByte(         long long &val){ SWAP_LLONG(val); }
  template <> inline void SwapByte(unsigned long long &val){ SWAP_LLONG(val); }
  #endif
  // N: native
  // B: big
  // L: little
  template<typename T> inline void SwapByteN2B( T &val ){ if( GetByteOrder() == 0 ) SwapByte(val); }
  template<typename T> inline void SwapByteN2L( T &val ){ if( GetByteOrder() == 1 ) SwapByte(val); }
  template<typename T> inline void SwapByteB2N( T &val ){ if( GetByteOrder() == 0 ) SwapByte(val); }
  template<typename T> inline void SwapByteL2N( T &val ){ if( GetByteOrder() == 1 ) SwapByte(val); }
  
  template<typename T> inline void SwapBytes   ( size_t size, T *arr ){ for(size_t i=0; i<size; ++i) SwapByte(arr[i]); }
  template<typename T> inline void SwapBytesN2B( size_t size, T *arr ){ if( GetByteOrder() == 0 ) SwapBytes(size, arr); }
  template<typename T> inline void SwapBytesN2L( size_t size, T *arr ){ if( GetByteOrder() == 1 ) SwapBytes(size, arr); }
  template<typename T> inline void SwapBytesB2N( size_t size, T *arr ){ if( GetByteOrder() == 0 ) SwapBytes(size, arr); }
  template<typename T> inline void SwapBytesL2N( size_t size, T *arr ){ if( GetByteOrder() == 1 ) SwapBytes(size, arr); }
  
  template<typename T> inline void SwapBytes   ( std::vector<T> &vct ){ for(typename std::vector<T>::iterator i=vct.begin(),iend=vct.end(); i!=iend; ++i) SwapByte(*i); }
  template<typename T> inline void SwapBytesN2B( std::vector<T> &vct ){ if( GetByteOrder() == 0 ) SwapBytes(vct); }
  template<typename T> inline void SwapBytesN2L( std::vector<T> &vct ){ if( GetByteOrder() == 1 ) SwapBytes(vct); }
  template<typename T> inline void SwapBytesB2N( std::vector<T> &vct ){ if( GetByteOrder() == 0 ) SwapBytes(vct); }
  template<typename T> inline void SwapBytesL2N( std::vector<T> &vct ){ if( GetByteOrder() == 1 ) SwapBytes(vct); }
  
  #if __cplusplus >= 201103L
  template<typename T, size_t N> inline void SwapBytes   ( std::array<T, N> &arr ){ for(T &t : arr ) SwapByte(t); }
  template<typename T, size_t N> inline void SwapBytesN2B( std::array<T, N> &arr ){ if( GetByteOrder() == 0 ) SwapBytes(arr); }
  template<typename T, size_t N> inline void SwapBytesN2L( std::array<T, N> &arr ){ if( GetByteOrder() == 1 ) SwapBytes(arr); }
  template<typename T, size_t N> inline void SwapBytesB2N( std::array<T, N> &arr ){ if( GetByteOrder() == 0 ) SwapBytes(arr); }
  template<typename T, size_t N> inline void SwapBytesL2N( std::array<T, N> &arr ){ if( GetByteOrder() == 1 ) SwapBytes(arr); }
  #endif
}

////////////////////////////////////////////////////////////////////////////////
#endif /* BYTESWAP_H */
