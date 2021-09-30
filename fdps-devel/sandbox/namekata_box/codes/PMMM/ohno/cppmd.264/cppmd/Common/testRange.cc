#include <iostream>
#include "Common.h"

template<int T>
ParticleRange getrange(TypeRange& tr)
{
  if(T==Coulomb){
    return tr.coulomb;
  }else if(T==LJCoulomb){
    return tr.ljcoulomb;
  }else if(T==LJ){
    return tr.lj;
  }else{
    return tr;
  }
}

int
main()
{
  TypeRange tr;

  tr.begin = 0;
  tr.end = 20;
  tr.lj.begin = 0;
  tr.lj.end = 5;
  tr.ljcoulomb.begin = 5;
  tr.ljcoulomb.end = 10;
  tr.coulomb.begin = 10;
  tr.coulomb.end = 20;

  ParticleRange all = getrange<ALL>(tr);
  ParticleRange lj = getrange<LJ>(tr);
  ParticleRange ljc = getrange<LJCoulomb>(tr);
  ParticleRange cl = getrange<Coulomb>(tr);

  std::cout << " " << lj.begin << " " << lj.end;
  std::cout << " " << ljc.begin << " " << ljc.end;
  std::cout << " " << cl.begin << " " << cl.end;
  std::cout << " " << all.begin << " " << all.end;
  std::cout << std::endl;
}
