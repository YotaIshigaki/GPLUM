#ifndef CPPMD_CPPMD_COVALENTBOND_CBUNORDERED_H_
#define CPPMD_CPPMD_COVALENTBOND_CBUNORDERED_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H

#if defined(HAVE_UNORDERED_MAP)
# include <unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP)
# include <tr1/unordered_map>
#elif defined(HAVE_BOOST_UNORDERED_MAP_HPP)
# include <boost/unordered_map.hpp>
#else  // HAVE_UNORDERED_MAP
# include <map>
#endif  // HAVE_UNORDERED_MAP

#ifdef HAVE_UNORDERED_MAP
namespace tr1 {}
namespace std { using namespace tr1; }
#endif  // HAVE_UNORDERD_MAP

namespace CBModule {

template <typename T1, typename T2>
struct CBUnordered {
#if defined(HAVE_UNORDERED_MAP)
  typedef std::unordered_map<T1,T2> Map;
#elif defined(HAVE_TR1_UNORDERED_MAP)
  typedef std::tr1::unordered_map<T1,T2> Map;
#elif defined(HAVE_BOOST_UNORDERED_MAP_HPP)
  typedef boost::unordered_map<T1,T2> Map;
#else  // HAVE_UNORDERED_MAP
  typedef std::map<T1,T2> Map;
#endif  // HAVE_UNORDERED_MAP
};

}  // namespace CBModule

#endif  // CPPMD_CPPMD_COVALENTBOND_CBUNORDERED_H_
