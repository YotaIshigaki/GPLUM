#ifndef VEC_HH_
#define VEC_HH_ 1

#include <istream>
#include <ostream>

template<int N, typename T>
class Vec {
 private:
  T array[N];

 public:
  Vec() {}
  Vec(const T &u)   { for (int i = 0; i < N; ++i) array[i] = u; }
  Vec(const Vec &v) { for (int i = 0; i < N; ++i) array[i] = v[i]; }
  Vec(const T p[N]) { for (int i = 0; i < N; ++i) array[i] = p[i]; }
  ~Vec() {}

  static const int size = N;
  typedef T value_type;

  template<typename U>
  operator Vec<N, U>() const {
    Vec<N, U> v;
    for (int i = 0; i < N; ++i) v[i] = static_cast<U>(array[i]);
    return v;
  }
  operator T*()             { return array; }
  operator const T*() const { return array; }

  T &operator[](int i)             { return array[i]; }
  const T &operator[](int i) const { return array[i]; }

#define VEC_HH__DEFINE_COMPARISON_OP(op, negop) \
  friend bool operator op(const Vec &x, const Vec &y) { \
    for (int i = 0; i < N; ++i) if (x[i] negop y[i]) return false; \
    return true; \
  }
  VEC_HH__DEFINE_COMPARISON_OP(<, >=)
  VEC_HH__DEFINE_COMPARISON_OP(<=, >)
  VEC_HH__DEFINE_COMPARISON_OP(==, !=)
  VEC_HH__DEFINE_COMPARISON_OP(>=, <)
  VEC_HH__DEFINE_COMPARISON_OP(>, <=)
#undef VEC_HH__DEFINE_COMPARISON_OP
  friend bool operator !=(const Vec &x, const Vec &y) {
    return !(x == y);
  }

#define VEC_HH__DEFINE_ASSIGN_OP(op) \
  Vec &operator op(const T &u) { \
    for (int i = 0; i < N; ++i) array[i] op u; \
    return *this; \
  } \
  Vec &operator op(const Vec &v) { \
    for (int i = 0; i < N; ++i) array[i] op v[i]; \
    return *this; \
  }
  VEC_HH__DEFINE_ASSIGN_OP(=)
  VEC_HH__DEFINE_ASSIGN_OP(+=)
  VEC_HH__DEFINE_ASSIGN_OP(-=)
  VEC_HH__DEFINE_ASSIGN_OP(*=)
  VEC_HH__DEFINE_ASSIGN_OP(/=)
  VEC_HH__DEFINE_ASSIGN_OP(%=)
  VEC_HH__DEFINE_ASSIGN_OP(^=)
  VEC_HH__DEFINE_ASSIGN_OP(&=)
  VEC_HH__DEFINE_ASSIGN_OP(|=)
  VEC_HH__DEFINE_ASSIGN_OP(>>=)
  VEC_HH__DEFINE_ASSIGN_OP(<<=)
#undef VEC_HH__DEFINE_ASSIGN_OP

#define VEC_HH__DEFINE_BINARY_OP(op) \
  friend Vec operator op(const Vec &x, const Vec &y) { \
    Vec v; \
    for (int i = 0; i < N; ++i) v[i] = (x[i] op y[i]); \
    return v; \
  } \
  friend Vec operator op(const T &u, const Vec &y) { \
    Vec v; \
    for (int i = 0; i < N; ++i) v[i] = (u op y[i]); \
    return v; \
  } \
  friend Vec operator op(const Vec &x, const T &u) { \
    Vec v; \
    for (int i = 0; i < N; ++i) v[i] = (x[i] op u); \
    return v; \
  }
  VEC_HH__DEFINE_BINARY_OP(+)
  VEC_HH__DEFINE_BINARY_OP(-)
  VEC_HH__DEFINE_BINARY_OP(*)
  VEC_HH__DEFINE_BINARY_OP(/)
  VEC_HH__DEFINE_BINARY_OP(%)

  VEC_HH__DEFINE_BINARY_OP(^)
  VEC_HH__DEFINE_BINARY_OP(&)
  VEC_HH__DEFINE_BINARY_OP(|)
  VEC_HH__DEFINE_BINARY_OP(>>)
  VEC_HH__DEFINE_BINARY_OP(<<)

  VEC_HH__DEFINE_BINARY_OP(&&)
  VEC_HH__DEFINE_BINARY_OP(||)
#undef VEC_HH__DEFINE_BINARY_OP

#define VEC_HH__DEFINE_UNARY_OP(op) \
  Vec operator op() const { \
    Vec v; \
    for (int i = 0; i < N; ++i) v[i] = op array[i]; \
    return v; \
  }
  VEC_HH__DEFINE_UNARY_OP(~)
  VEC_HH__DEFINE_UNARY_OP(!)
  VEC_HH__DEFINE_UNARY_OP(-)
  VEC_HH__DEFINE_UNARY_OP(+)
#undef VEC_HH__DEFINE_UNARY_OP

#define VEC_HH__DEFINE_PREFIX_OP(op) \
  Vec &operator op() { \
    for (int i = 0; i < N; ++i) op array[i]; \
    return *this; \
  }
  VEC_HH__DEFINE_PREFIX_OP(++)
  VEC_HH__DEFINE_PREFIX_OP(--)
#undef VEC_HH__DEFINE_PREFIX_OP

#define VEC_HH__DEFINE_POSTFIX_OP(op) \
  Vec operator op(int) { \
    Vec v; \
    for (int i = 0; i < N; ++i) v[i] = array[i] op; \
    return v; \
  }
  VEC_HH__DEFINE_POSTFIX_OP(++)
  VEC_HH__DEFINE_POSTFIX_OP(--)
#undef VEC_HH__DEFINE_POSTFIX_OP

  friend std::ostream &operator<<(std::ostream &os, const Vec &v) {
    if (N > 0) {
      os << v[0];
      for (int i = 1; i < N; ++i) { os << " " << v[i]; }
    }
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Vec &v) {
    if (N > 0) {
      is >> v[0];
      for (int i = 1; i < N; ++i) { is.ignore(); is >> v[i]; }
    }
    return is;
  }
};

typedef Vec<3,int> V3i;
typedef Vec<3,double> V3d;

typedef Vec<4,int> V4i;
typedef Vec<4,double> V4d;

#endif  // VEC_HH_
