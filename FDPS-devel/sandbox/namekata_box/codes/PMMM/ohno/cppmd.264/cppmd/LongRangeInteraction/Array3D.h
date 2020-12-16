#ifndef CPPMD_CPPMD_LONGRANGEINTERACTION_ARRAY3D_H_
#define  CPPMD_CPPMD_LONGRANGEINTERACTION_ARRAY3D_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H

#if defined(HAVE_ARRAY)
# include <array>
#elif defined(HAVE_TR1_ARRAY)
# include <tr1/array>
#elif defined(HAVE_BOOST_ARRAY_HPP)
# include <boost/array.hpp>
#else  // HAVE_ARRAY
# error Array3D array
#endif  // HAVE_ARRAY

#ifdef HAVE_BOOST_MULTI_ARRAY_HPP
# include <boost/multi_array.hpp>
#else  // HAVE_BOOST_MULTI_ARRAY_HPP
# error Array3D multi_array
#endif  // HAVE_BOOST_MULTI_ARRAY_HPP

#ifdef HAVE_ARRAY
namespace tr1 {}
namespace std { using namespace tr1; }
#endif  // HAVE_ARRAY

namespace PMEModule {

template <typename T>
struct Array3D {
#ifdef HAVE_BOOST_MULTI_ARRAY_HPP
  typedef boost::multi_array<T, 3> Array;
  typedef boost::multi_array_ref<T, 3> ArrayRef;
#endif  // HAVE_BOOST_MULTI_ARRAY_HPP

#if defined(HAVE_ARRAY)
  typedef std::array<typename ArrayRef::index, 3> Dims;
#elif defined(HAVE_TR1_ARRAY)
  typedef std::tr1::array<typename ArrayRef::index, 3> Dims;
#elif defined(HAVE_BOOST_ARRAY_HPP)
  typedef boost::array<typename ArrayRef::index, 3> Dims;
#endif  // HAVE_ARRAY

  static Array* createArray(int nx, int ny, int nz) {
    return new Array(boost::extents[nx][ny][nz]);
  }
  static ArrayRef* createArray(T* data, int nx, int ny, int nz) {
    return new ArrayRef(data, boost::extents[nx][ny][nz]);
  }
};

}  // namespace PMEModule

#endif  // CPPMD_CPPMD_LONGRANGEINTERACTION_ARRAY3D_H_
