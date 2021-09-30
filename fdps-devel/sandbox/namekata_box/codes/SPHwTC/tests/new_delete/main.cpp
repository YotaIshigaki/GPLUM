#include <cstdio>
#include <iostream>


int main(int argc, char* argv[])
{
   //* Local variables
   double *x,*y,*z;

   //* Memory allocation
   std::cout << "Allocate memory!" << std::endl;
   x = new double[16];
   y = new double[16];
   z = new double[16];

   //* Release emmory
   std::cout << "Release memory!" << std::endl;
   delete[] x,y,z;

   return 0;
}
