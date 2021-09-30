#include <cstdio>
#include <cstdlib>
#include <iostream>

void sub();

int main(int argc, char* argv[]) {

   sub();
   sub();
   sub();

   return 0;

}


void sub() {
   //* Static variables
   static int ndump=0;

   //* Output ndump
   std::cout << "ndump = " << ndump << std::endl;

   //* Update ndump
   ndump++;

}
