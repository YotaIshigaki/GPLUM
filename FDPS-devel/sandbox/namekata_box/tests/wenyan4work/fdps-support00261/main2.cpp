#include <iostream>

typedef double F64;

//getRserach があるかないかを判定する構造体（クラス）
template<typename T>
struct HasgetRsearchMethod
{
   template<typename U, F64(U::*)() > struct SFINAE {};
   template<typename U>
   static char Test(SFINAE<U, &U::getRsearch> *);
   template<typename U> static int Test(...);
   static const bool value = sizeof(Test<T>(0)) == sizeof(char);
};

// getRsearch がある時にそれを使う関数
template<class T>
F64 mygetRsearch(T t, std::true_type)
{
   return t.getRsearch();
}

// getRsearch がない時になにか返す関数
template<class T>
F64 mygetRsearch(T t, std::false_type)
{
   return 0.0;
}

// getRsearch があるかないかで呼びわける関数
template<class T>
F64 mygetRsearch(T t)
{
   return mygetRsearch(t, std::integral_constant<bool, HasgetRsearchMethod<T>::value>());
}


//getRsearchメソッドのある/ないHogeクラス
class HogeWithgetRsearch
{
public:
   F64 rsearchval;
   void setRsearch(F64 val) { rsearchval=val; }
   F64 getRsearch() { return rsearchval; }
};
class Hoge{};

int main() 
{
   HogeWithgetRsearch p1, p2;
   Hoge p3;
   p1.setRsearch(1.2);
   p2.setRsearch(2.0);
   std::cout << "p1="<<mygetRsearch(p1)<<std::endl;
   std::cout << "p2="<<mygetRsearch(p2)<<std::endl;
   std::cout << "p3="<<mygetRsearch(p3)<<std::endl;
   return 0;
}
