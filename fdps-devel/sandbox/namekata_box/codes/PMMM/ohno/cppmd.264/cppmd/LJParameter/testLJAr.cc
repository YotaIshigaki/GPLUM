#include <iostream>
#include "LJAr.h"

int
main()
{

  
  LJAr lja;

  std::cout << lja.getSigma() << std::endl;
  std::cout << lja.getEpsilon() << std::endl;

  std::cout << ljmixAr[0][0].potFirst << std::endl;
  std::cout << ljmixAr[0][0].potSecond << std::endl;
  std::cout << ljmixAr[0][0].forceFirst << std::endl;
  std::cout << ljmixAr[0][0].forceSecond << std::endl;
  
}
