#include <iostream>
#include <typeinfo>
#include <type_traits>

static std::string getNameByTypeInfo(std::type_info const& iTypeInfo)
{
    return iTypeInfo.name();
}
#define TYPENAME(dType) getNameByTypeInfo(typeid(dType))

struct SEARCH_MODE_LONG {};
struct SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE {};

template
<
    typename T,
    typename std::enable_if< std::is_same<T, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>::value >::type* = nullptr
>
void foo(T)
{
    std::cout << "[0] foo<" << TYPENAME(T) << ">\n";
}
 
template
<
    typename T,
    typename std::enable_if < !std::is_same<T, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>::value >::type* = nullptr
>
void foo(T)
{
    std::cout << "[1] foo<" << TYPENAME(T) << ">\n";
}
 
int main()
{
    std::cout << std::is_same<SEARCH_MODE_LONG, SEARCH_MODE_LONG>::value << std::endl;
    std::cout << std::is_same<SEARCH_MODE_LONG, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>::value << std::endl;
    std::cout << typeid( std::enable_if<true, SEARCH_MODE_LONG>::type* ).name() << std::endl;
    std::cout << typeid( std::enable_if<true>::type* ).name() << std::endl;
    std::cout << typeid( std::enable_if< std::is_same<SEARCH_MODE_LONG, SEARCH_MODE_LONG>::value >::type* ).name() << std::endl;

    foo(123);
    foo(1.0);
    SEARCH_MODE_LONG sm_long;
    SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE sm_long_pmm;
    foo(sm_long);
    foo(sm_long_pmm);
    return 0;
}
