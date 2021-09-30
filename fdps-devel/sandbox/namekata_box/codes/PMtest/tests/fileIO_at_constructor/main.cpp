#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

class Calc_Conditions {
   public:
      const double x;
      Calc_Conditions();
};

Calc_Conditions::Calc_Conditions() {
   //* Local variables
   std::string buf;
   std::ifstream ifs;

   ifs.open("params.dat",std::ios::in);
   if (ifs.fail()) {
      std::cout << "cannot open params.dat" << std::endl;
      std::exit(1);
   }
   while(std::getline(ifs,buf)) {
      std::cout << buf << std::endl;
   }
   ifs.close();
   this->x = 1.0e0;
}


int main(int argc, char* argv[]) {
   //* Local variables
   Calc_Conditions calc_condi;

   //std::cout << calc_condi.x << std::endl;


   return 0;
}
