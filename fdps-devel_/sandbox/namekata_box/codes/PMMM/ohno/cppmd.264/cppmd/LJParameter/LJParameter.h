#ifndef LJPARAMETER_H
#define LJPARAMETER_H

#ifdef USE_BOOST_MULTI_ARRAY
# include <boost/multi_array.hpp>
#endif  // USE_BOOST_MULTI_ARRAY
#include "Common.h"

#include <cassert>
#include <climits>
#include <string>

//! LJMixParameter
typedef struct LJMixParameter{
  double potFirst;
  double potSecond;
  double forceFirst;
  double forceSecond;
}LJMixParameter;

#ifdef USE_BOOST_MULTI_ARRAY
template <class T>
class array2d : public boost::multi_array<T,2> {
 public:
  array2d()
  {
  }

  explicit array2d(int num)
  {
    array2d::resize(boost::extents[num][num]);
  }

  array2d(int num_i, int num_j)
  {
    array2d::resize(boost::extents[num_i][num_j]);
  }

  boost::multi_array<T,2>& resize2d(int num_i, int num_j)
  {
    return array2d::resize(boost::extents[num_i][num_j]);
  }
};
#else  // USE_BOOST_MULTI_ARRAY
template <typename T>
struct array2d {
  T **array;
  int size_i;
  int size_j;

  array2d()
      : size_i(0), size_j(0)
  {
  }

  explicit array2d(int num)
      : size_i(num), size_j(num)
  {
    assert(0 <= num && num <= INT_MAX);
    resize2d(size_i,size_j);
  }


  array2d(int num_i, int num_j)
      : size_i(num_i), size_j(num_j)
  {
    assert(0 <= num_i && num_i <= INT_MAX);
    assert(0 <= num_j && num_j <= INT_MAX);
    resize2d(size_i,size_j);
  }

  array2d(const array2d& oa)
  {
    resize2d(oa.size_i,oa.size_j);
    for(int i=0;i<size_i;i++){
      for(int j=0;j<size_j;j++){
        array[i][j] = oa.array[i][j];
      }
    }
  }

  array2d& operator=(const array2d& oa)
  {
    resize2d(oa.size_i,oa.size_j);
    for(int i=0;i<size_i;i++){
      for(int j=0;j<size_j;j++){
        array[i][j] = oa.array[i][j];
      }
    }
    return *this;
  }

  ~array2d()
  {
    for(int i=0;i<size_i;i++){
      delete [] array[i];
    }
    delete [] array;
  }

  T * operator [](const int i){
    return array[i];
  }

  inline size_t size()
  {
    return size_i*size_j;
  }

  inline size_t resize2d(int num_i, int num_j)
  {
    for(int i=0;i<size_i;i++){
      delete [] array[i];
    }
    delete [] array;
    size_i = num_i;
    size_j = num_j;
    array = new T*[size_i];
    for(int i=0;i<size_i;i++){
      array[i] = new T[size_j];
    }
    return size();
  }
};
#endif  // USE_BOOST_MULTI_ARRAY
typedef array2d<std::string> LJMixName;
typedef array2d<LJMixParameter> LJMixparameterArray;

//! LJ Parameter
class LJParameter {
 public:
  LJParameter() : name(), sigma(0), epsilon(0), potFirst(0), potSecond(0), 
                  forceFirst(0), forceSecond(0) {}
  LJParameter(const std::string _name, const double _sigma, const double _epsilon);
  //  void setSigmaEpsilon(const double _sigma, const double _epsilon);
  const std::string& getName() const { return name;}
  double getSigma() const { return sigma;}
  double getEpsilon() const { return epsilon;}
  double getPotentialFirst() const { return potFirst;}
  double getPotentialSecond() const { return potSecond;}
  double getForceFirst() const { return forceFirst;}
  double getForceSecond() const { return forceSecond;}
 private:
  void setLocalParm();
  std::string name;
  double sigma,epsilon;
  double potFirst,potSecond;
  double forceFirst,forceSecond;
};

struct LJParameterCoulomb {
  array2d<LJParameter> ljparameter;

  LJParameterCoulomb(int num) : ljparameter(num,num){}
};

struct LJParameterEwaldReal {
  array2d<LJParameter> ljparameter;
  double alpha;
  double alpha2;

  LJParameterEwaldReal(int num) : ljparameter(num,num) {}
};

namespace LJParameterStuff {
  const LJParameter mixParameter(const LJParameter& p1, const LJParameter& p2);
}

#endif
