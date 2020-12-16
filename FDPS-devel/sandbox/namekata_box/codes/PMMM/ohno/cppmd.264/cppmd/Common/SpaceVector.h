#ifndef SPACEVECTOR_H
#define SPACEVECTOR_H

#include <cmath>
#if 1
# include <cstdio>
#endif
#include <iostream>

template<typename T>
class SpaceVector {
 public:
  SpaceVector() : x(T(0)), y(T(0)), z(T(0)) {}
  explicit SpaceVector(const T _v) : x(_v),y(_v),z(_v) {}
  SpaceVector(const T _x, const T _y, const T _z) : x(_x),y(_y),z(_z){}
  ~SpaceVector(){};

  enum {
    Dim = 3
  };

  void set(T _x, T _y, T _z){ x = _x; y = _y; z = _z;}
  void setX(T _x){ x = _x;}
  void setY(T _y){ y = _y;}
  void setZ(T _z){ z = _z;}
  void addX(T _x){ x += _x;}
  void addY(T _y){ y += _y;}
  void addZ(T _z){ z += _z;}

  const T getX() const {return x;}
  const T getY() const {return y;}
  const T getZ() const {return z;}
  const T get(const int i) const { return (&x)[i];}

  const T getComponentMax() const {
    return std::max<T>(x,std::max<T>(y,z));
  }
  const T getComponentMin() const {
    return std::min<T>(x,std::min<T>(y,z));
  }
  
  const T norm2() const {
    return (*this)*(*this);
  }
  const T norm() const {
    return sqrt((*this)*(*this));
  }
  const T angle(const SpaceVector<T> &v) const {
    return acos((*this)*v/(norm()*v.norm()));
  }

  T& operator [](const int i){
    return static_cast<T*>(&x)[i];
  }
  const T operator [](const int i) const {
    return (&x)[i];
  }
  const SpaceVector<T> operator-() const{
    return SpaceVector<T> (-x, -y, -z);
  }
  const SpaceVector<T> operator+(const SpaceVector<T> &v) const{
    return SpaceVector<T> (x+v.x, y+v.y, z+v.z);
  }
  const SpaceVector<T> operator-(const SpaceVector<T> &v) const{
    return SpaceVector<T> (x-v.x, y-v.y, z-v.z);
  }
  const SpaceVector<T> operator*(const T &s) const{
    return SpaceVector<T> (x*s, y*s, z*s);
  }
  friend const SpaceVector<T> operator*(const T &s, const SpaceVector<T> &v){
    return SpaceVector<T>(v*s);
  }

  const T operator*(const SpaceVector<T> &v) const {
    return (x*v.x + y*v.y + z*v.z);
  }
  const SpaceVector<T> operator% (const SpaceVector<T> &v) const{
    return SpaceVector<T> (y*v.z - z*v.y, 
                           z*v.x - x*v.z,
                           x*v.y - y*v.x);
  }
  const SpaceVector<T> operator/ (const T &s) const{
    T r = T(1)/s;
    return (*this)*r;
  }
  const SpaceVector<T>& operator= (const SpaceVector<T> &v){
    x = v.x; y = v.y; z = v.z;
#ifdef FJ_MAGIC
#else  // FJ_MAGIC
    return *this;
#endif  // FJ_MAGIC
  }
  const SpaceVector<T>& operator= (const T v){
    x = v; y = v; z = v;
    return *this;
  }

  const SpaceVector<T> operator- (){
    return SpaceVector<T> (-x, -y, -z);
  }
  SpaceVector<T>& operator+= (const SpaceVector<T> &v){
    *this = *this + v;
    return *this;
  }
  SpaceVector<T>& operator-= (const SpaceVector<T> &v){
    *this = *this - v;
    return *this;
  }
  SpaceVector<T>& operator*= (const T &s){
    *this = *this * s;
    return *this;
  }
  SpaceVector<T>& operator/= (const T &s){
    *this = *this / s;
    return *this;
  }

  bool operator==(const SpaceVector<T> &v) const {
    return x == v.x && y == v.y && z == v.z;
  }
  bool operator!=(const SpaceVector<T> &v) const {
    return x != v.x || y != v.y || z != v.z;
  }
  bool operator<(const SpaceVector<T> &v) const {
    return x < v.x && y < v.y && z < v.z;
  }
  bool operator<=(const SpaceVector<T> &v) const {
    return x <= v.x && y <= v.y && z <= v.z;
  }
  bool operator>(const SpaceVector<T> &v) const {
    return x > v.x && y > v.y && z > v.z;
  }
  bool operator>=(const SpaceVector<T> &v) const {
    return x >= v.x && y >= v.y && z >= v.z;
  }

  SpaceVector<T>& operator()(const T _x, const T _y, const T _z){
    x = _x; y = _y; z = _z;
    return *this;
  }
  SpaceVector<T>& operator()(const T s){
    x = s; y = s; z = s;
    return *this;
  }
  friend std::ostream& operator<< (std::ostream &os, const SpaceVector<T> &v){
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
  }
  friend std::istream& operator>> (std::istream &is, const SpaceVector<T> &v){
    is >> v.x >> v.y >> v.z;
    return is;
  }

#if 1
  friend void fprint(FILE *ofs, const SpaceVector<T> &v) {
    fprintf(ofs, "(%f,%f,%f)", v.x, v.y, v.z);
  }
  friend void print(const SpaceVector<T> &v) {
    printf("(%f,%f,%f)", v.x, v.y, v.z);
  }
#endif

  T x,y,z;
};

#endif
