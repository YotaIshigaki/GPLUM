/* Standard headers */
#include <iostream>
/* User-defined headers */
#include "sub.h"

int main(int argc, char *argv[]) {

   glbvar = 1.0e0;
   glbobj.var = 2.0e0;
   glbobj.Initialize(argc,argv);

   return 0;
}
