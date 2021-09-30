#include <cstdio>
#include <cstdlib>
#include <iostream>


int main(int argc, char* argv[]) {

   bool flag;
   double x=1.0e0;
   double y=2.0e0;

   flag = (x > y);
   std::cout << "flag = " << flag << std::endl;

   return 0;
}
